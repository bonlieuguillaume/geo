"""SAR acquisition opportunities over an area of interest, from orbit parameters.

Given a satellite orbit — a TLE, a NORAD id whose latest TLE is fetched from
CelesTrak, mean orbital elements for a satellite not flying yet, or the
operator's ephemerides — list the passes during which a side-looking radar can
image an area of interest: the zero-Doppler time over the area, the orbit
direction, the side the area lies on, and the incidence angle on the ground,
with how much of the area falls in the instrument's usable incidence range.

Propagation is SGP4 (``sgp4`` package) for TLEs and elements; ephemerides are
interpolated. Everything else — frames, zero-Doppler search, incidence, swath
edges — is numpy here.

Usable as a library — `satellite` builds the orbit, `find_passes` returns a
DataFrame, `swath_footprint` the imaged strip — or through the ipyleaflet +
ipywidgets interface returned by `build_ui`, which ``tasking.ipynb`` drives.
"""

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import shapely
from sgp4.api import WGS72, Satrec
from shapely import wkt as shapely_wkt
from shapely.geometry import LineString, Polygon, mapping, shape
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

# WGS84 ellipsoid: the ground
WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563
WGS84_B = WGS84_A * (1 - WGS84_F)
WGS84_E2 = WGS84_F * (2 - WGS84_F)
OMEGA_EARTH = 7.292115146706979e-5  # rad/s

# WGS72 constants: the ones SGP4 mean elements are defined with
MU_WGS72 = 398600.8  # km^3/s^2
RE_WGS72 = 6378.135  # km
J2_WGS72 = 0.001082616

# Rate of the mean Sun, which a sun-synchronous orbit plane turns at
SUN_RATE = 2 * np.pi / (365.2421897 * 86400.0)  # rad/s

JD_UNIX_EPOCH = 2440587.5
JD_SGP4_EPOCH = 2433281.5  # 1949-12-31 00:00 UT, origin of sgp4init's epoch
JD_J2000 = 2451545.0

CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php"
# Beyond a few days from its epoch, a TLE drifts by several km along track
TLE_MAX_AGE_DAYS = 3.0

# Sentinel-1 IW swath: the default usable incidence range
DEFAULT_INCIDENCE_RANGE = (29.1, 46.0)
LOOK_SIDES = ("right", "left", "both")
DIRECTIONS = ("ascending", "descending")

PASS_COLUMNS = [
    "time_utc", "direction", "side", "incidence_ref", "incidence_min",
    "incidence_max", "coverage", "off_nadir_ref", "slant_range_km",
    "local_time", "sat_lon", "sat_lat", "sat_alt_km", "heading_deg",
]


# --- Area of interest ---------------------------------------------------------

def _geometry_from_geojson(obj):
    """Build a single shapely geometry from a parsed GeoJSON object."""
    kind = obj.get("type")
    if kind == "FeatureCollection":
        geoms = [shape(f["geometry"]) for f in obj.get("features", []) if f.get("geometry")]
    elif kind == "Feature":
        geoms = [shape(obj["geometry"])]
    elif kind == "GeometryCollection":
        geoms = [shape(g) for g in obj.get("geometries", [])]
    elif kind:
        geoms = [shape(obj)]
    else:
        raise ValueError("Unrecognised GeoJSON: no 'type' member")
    if not geoms:
        raise ValueError("GeoJSON holds no geometry")
    return geoms[0] if len(geoms) == 1 else unary_union(geoms)


def parse_aoi(aoi):
    """Parse an area of interest into a shapely geometry (lon/lat, EPSG:4326).

    Parameters
    ----------
    aoi : shapely geometry | dict | str | Path
        A shapely geometry is returned as is; a dict is read as GeoJSON; a
        string is the path of a WKT or GeoJSON file if such a file exists,
        otherwise inline WKT or inline GeoJSON.
    """
    if isinstance(aoi, BaseGeometry):
        return aoi
    if isinstance(aoi, dict):
        return _geometry_from_geojson(aoi)

    text = str(aoi).strip()
    if not text:
        raise ValueError("empty area of interest")
    if "\n" not in text and len(text) < 4096:
        try:
            if Path(text).is_file():
                text = Path(text).read_text(encoding="utf-8").strip()
        except (OSError, ValueError):
            pass  # not a usable path, treat the string as a geometry

    if text.startswith("{"):
        return _geometry_from_geojson(json.loads(text))
    try:
        return shapely_wkt.loads(text)
    except Exception as exc:  # shapely raises several unrelated types here
        raise ValueError(f"not a WKT geometry, not GeoJSON: {text[:60]!r}") from exc


def unwrap_antimeridian(geom):
    """Move an AOI that crosses the antimeridian to continuous longitudes.

    Shapely knows nothing of the sphere: a polygon from 179° to -179° spans
    358° for it. Wider than 180°, an AOI is taken as crossing the
    antimeridian, and its negative longitudes are moved by +360° — the same
    points on the ground, beyond 180°.
    """
    geom = parse_aoi(geom)
    minx, _, maxx, _ = geom.bounds
    if maxx - minx <= 180.0:
        return geom
    return shapely.transform(
        geom, lambda c: np.column_stack([np.where(c[:, 0] < 0, c[:, 0] + 360.0, c[:, 0]), c[:, 1]])
    )


def aoi_points(aoi, spacing_km=5.0, max_points=400):
    """Ground points an AOI is checked at: a reference point, a grid, the boundary.

    Incidence grows monotonically with ground range, so over an area its
    extremes lie on the boundary: sampling the boundary bounds the incidence
    over the whole area. The grid inside measures which share of the area a
    pass covers.

    Parameters
    ----------
    aoi : anything `parse_aoi` accepts
    spacing_km : float
        Spacing of the points, along the boundary and in the grid.
    max_points : int
        Cap on the boundary points and on the grid points, each; the spacing
        grows for large areas.

    Returns
    -------
    points : numpy.ndarray
        (n, 2) lon/lat in degrees. Row 0 is the reference point of the AOI (a
        point inside it), the one passes are timed on; rows 1 to `n_area` the
        grid; the rest the boundary.
    n_area : int
        Number of grid points (0 for an AOI smaller than the spacing, whose
        area is then represented by the reference point alone).
        Longitudes may exceed 180° for an AOI across the antimeridian.
    """
    geom = unwrap_antimeridian(aoi)
    inner = geom.representative_point()
    step = spacing_km / 111.0

    grid = np.empty((0, 2))
    if geom.area > 0:
        grid_step = max(step, np.sqrt(geom.area / max_points))
        minx, miny, maxx, maxy = geom.bounds
        gx, gy = np.meshgrid(
            np.arange(minx + grid_step / 2, maxx, grid_step),
            np.arange(miny + grid_step / 2, maxy, grid_step),
        )
        inside = shapely.contains_xy(geom, gx, gy)
        grid = np.column_stack([gx[inside], gy[inside]])

    line = geom.boundary if geom.area > 0 else geom
    if line.length / step > max_points:
        step = line.length / max_points
    coords = shapely.get_coordinates(shapely.segmentize(line, step)) if step > 0 else np.empty((0, 2))
    coords = np.unique(np.round(coords, 9), axis=0)
    return np.vstack([[inner.x, inner.y], grid, coords]), len(grid)


# --- Time and frames ----------------------------------------------------------

def _utc(value):
    """Any time (str, datetime, date, Timestamp) as a naive UTC datetime64[ns]."""
    if isinstance(value, np.datetime64):
        return value.astype("datetime64[ns]")
    ts = pd.Timestamp(value)
    if ts.tzinfo is not None:
        ts = ts.tz_convert("UTC").tz_localize(None)
    return ts.to_datetime64().astype("datetime64[ns]")


def _julian(times):
    """Split datetime64 times into the (jd, fraction) pair SGP4 takes."""
    ns = np.atleast_1d(np.asarray(times, dtype="datetime64[ns]")).astype(np.int64)
    day_ns = 86_400 * 10**9
    days = ns // day_ns
    return JD_UNIX_EPOCH + days.astype(float), (ns - days * day_ns) / day_ns


def _gmst(jd, fr):
    """Greenwich mean sidereal time (IAU 1982, the one TEME is tied to), radians.

    UT1 is taken as UTC: at most 0.9 s of Earth rotation, 400 m at the equator.
    """
    t = (jd - JD_J2000 + fr) / 36525.0
    sec = 67310.54841 + (876600.0 * 3600 + 8640184.812866) * t + 0.093104 * t**2 - 6.2e-6 * t**3
    return np.radians((sec % 86400.0) / 240.0)


def _rot(axis, angle):
    """Frame rotations R1/R2/R3 by `angle` (rad, array), (n, 3, 3)."""
    c, s = np.cos(angle), np.sin(angle)
    o, z = np.ones_like(c), np.zeros_like(c)
    rows = {
        1: [[o, z, z], [z, c, s], [z, -s, c]],
        2: [[c, z, -s], [z, o, z], [s, z, c]],
        3: [[c, s, z], [-s, c, z], [z, z, o]],
    }[axis]
    return np.moveaxis(np.array(rows), [0, 1], [-2, -1])


def _precession_nutation(jd, fr):
    """J2000 to true-of-date rotation (n, 3, 3) and equation of the equinoxes (rad).

    IAU 1976 precession and the four main IAU 1980 nutation terms: ~30 m on
    the satellite position, against the full 106-term series. UTC stands
    for TT (69 s, nothing at this level).
    """
    t = (jd - JD_J2000 + fr) / 36525.0
    arcsec = np.pi / (180.0 * 3600.0)
    zeta = (2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3) * arcsec
    z = (2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3) * arcsec
    theta = (2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3) * arcsec
    precession = _rot(3, -z) @ _rot(2, theta) @ _rot(3, -zeta)

    node = np.radians(125.04452 - 1934.136261 * t)    # Moon's ascending node
    sun = np.radians(280.4665 + 36000.7698 * t)       # Sun's mean longitude
    moon = np.radians(218.3165 + 481267.8813 * t)     # Moon's mean longitude
    dpsi = (-17.20 * np.sin(node) - 1.32 * np.sin(2 * sun) - 0.23 * np.sin(2 * moon)
            + 0.21 * np.sin(2 * node)) * arcsec
    deps = (9.20 * np.cos(node) + 0.57 * np.cos(2 * sun) + 0.10 * np.cos(2 * moon)
            - 0.09 * np.cos(2 * node)) * arcsec
    eps = (84381.448 - 46.8150 * t - 0.00059 * t**2 + 0.001813 * t**3) * arcsec
    nutation = _rot(1, -(eps + deps)) @ _rot(3, -dpsi) @ _rot(1, eps)
    return nutation @ precession, dpsi * np.cos(eps + deps)


INERTIAL_FRAMES = ("TEME", "EME2000")


def inertial_to_ecef(times, r, v, frame="TEME"):
    """Inertial position/velocity to Earth-fixed, relative to the ground.

    Parameters
    ----------
    times : array of datetime64
        UTC times.
    r, v : numpy.ndarray
        (n, 3), in any unit (velocity in that unit per second).
    frame : {"TEME", "EME2000"}
        TEME is SGP4's frame, rotated by GMST alone. EME2000 (J2000, and GCRF
        within a metre) goes through precession and nutation first.

    Returns
    -------
    r, v : numpy.ndarray
        (n, 3) ECEF, same units. UT1 is taken as UTC and polar motion is
        ignored: a few hundred metres at most.
    """
    jd, fr = _julian(times)
    theta = _gmst(jd, fr)
    if frame == "EME2000":
        pn, eqeq = _precession_nutation(jd, fr)
        r, v = np.einsum("nij,nj->ni", pn, r), np.einsum("nij,nj->ni", pn, v)
        theta = theta + eqeq  # true sidereal time, from true-of-date axes
    elif frame != "TEME":
        raise ValueError(f"inertial frame must be one of {INERTIAL_FRAMES}, not {frame!r}")
    c, s = np.cos(theta), np.sin(theta)
    xe, ye = c * r[:, 0] + s * r[:, 1], -s * r[:, 0] + c * r[:, 1]
    r_ecef = np.column_stack([xe, ye, r[:, 2]])
    # Velocity seen from the rotating Earth: v - omega x r
    v_ecef = np.column_stack([
        c * v[:, 0] + s * v[:, 1] + OMEGA_EARTH * ye,
        -s * v[:, 0] + c * v[:, 1] - OMEGA_EARTH * xe,
        v[:, 2],
    ])
    return r_ecef, v_ecef


def propagate(sat, times):
    """Satellite position and velocity in the Earth-fixed frame.

    Parameters
    ----------
    sat : sgp4.api.Satrec | Ephemeris
    times : array of datetime64
        UTC times.

    Returns
    -------
    r, v : numpy.ndarray
        (n, 3) position (m) and velocity relative to the ground (m/s), ECEF.
    """
    times = np.atleast_1d(np.asarray(times, dtype="datetime64[ns]"))
    if isinstance(sat, Ephemeris):
        return sat.state(times)
    jd, fr = _julian(times)
    err, r, v = sat.sgp4_array(jd, fr)
    if np.any(err):
        bad = int(np.flatnonzero(err)[0])
        raise RuntimeError(
            f"SGP4 error {int(err[bad])} at jd {jd[bad] + fr[bad]:.5f}: decayed orbit, "
            "or elements used far from their epoch"
        )
    return inertial_to_ecef(times, r * 1e3, v * 1e3, "TEME")


def geodetic_to_ecef(lon, lat, height=0.0):
    """WGS84 lon/lat (degrees) and ellipsoidal height (m) to ECEF (m), (n, 3)."""
    lon, lat = np.radians(lon), np.radians(lat)
    n = WGS84_A / np.sqrt(1 - WGS84_E2 * np.sin(lat) ** 2)
    return np.column_stack([
        (n + height) * np.cos(lat) * np.cos(lon),
        (n + height) * np.cos(lat) * np.sin(lon),
        (n * (1 - WGS84_E2) + height) * np.sin(lat),
    ])


def ecef_to_geodetic(xyz):
    """ECEF (m), (n, 3), to WGS84 lon, lat (degrees) and height (m)."""
    x, y, z = np.atleast_2d(xyz).T
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1 - WGS84_E2))
    for _ in range(5):
        n = WGS84_A / np.sqrt(1 - WGS84_E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1 - WGS84_E2 * n / (n + h)))
    n = WGS84_A / np.sqrt(1 - WGS84_E2 * np.sin(lat) ** 2)
    return np.degrees(np.arctan2(y, x)), np.degrees(lat), p / np.cos(lat) - n


def _up(lon, lat):
    """Ellipsoid normal (unit, ECEF) at lon/lat in degrees, (n, 3)."""
    lon, lat = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def _dot(a, b):
    return np.einsum("ij,ij->i", a, b)


def _unit(a):
    return a / np.linalg.norm(a, axis=1, keepdims=True)


# --- Orbit --------------------------------------------------------------------

def parse_tle(text):
    """Split a TLE (two lines, or three with a name first) into (name, line1, line2)."""
    lines = [ln.strip() for ln in str(text).strip().splitlines() if ln.strip()]
    l1 = next((ln for ln in lines if ln.startswith("1 ")), None)
    l2 = next((ln for ln in lines if ln.startswith("2 ")), None)
    if l1 is None or l2 is None:
        raise ValueError("not a TLE: no line starting with '1 ' and '2 '")
    i = lines.index(l1)
    name = lines[i - 1].lstrip("0 ").strip() if i > 0 else ""
    return name, l1, l2


def fetch_tle(norad_id, timeout=30):
    """Latest TLE of a satellite from CelesTrak (anonymous), as (name, line1, line2)."""
    resp = requests.get(
        CELESTRAK_URL, params={"CATNR": int(norad_id), "FORMAT": "TLE"}, timeout=timeout
    )
    resp.raise_for_status()
    if "1 " not in resp.text:
        raise ValueError(f"CelesTrak has no TLE for NORAD id {norad_id}: {resp.text.strip()[:80]!r}")
    return parse_tle(resp.text)


def _mean_sun_ra(jd):
    """Right ascension of the mean Sun (its mean longitude), radians."""
    return np.radians((280.460 + 0.9856474 * (jd - JD_J2000)) % 360.0)


def _sso_inclination(a_km, ecc=0.0):
    """Inclination (rad) that makes the J2 node drift follow the mean Sun."""
    n = np.sqrt(MU_WGS72 / a_km**3)
    cos_i = -SUN_RATE * (1 - ecc**2) ** 2 / (1.5 * n * J2_WGS72 * (RE_WGS72 / a_km) ** 2)
    if abs(cos_i) > 1:
        raise ValueError(f"no sun-synchronous orbit at a = {a_km:.0f} km")
    return np.arccos(cos_i)


def _node_rate(a_km, inc, ecc=0.0):
    """J2 secular drift of the ascending node, rad/s."""
    n = np.sqrt(MU_WGS72 / a_km**3)
    return -1.5 * n * J2_WGS72 * (RE_WGS72 / a_km) ** 2 * np.cos(inc) / (1 - ecc**2) ** 2


def _nodal_period(sat, days=1.0, step_s=10.0):
    """Mean time between ascending nodes over `days` from the epoch, seconds."""
    t = np.arange(0.0, days * 86400.0 + step_s, step_s)
    jd = np.full(t.shape, sat.jdsatepoch)
    _, r, _ = sat.sgp4_array(jd, sat.jdsatepochF + t / 86400.0)
    z = r[:, 2]
    i = np.flatnonzero((z[:-1] < 0) & (z[1:] >= 0))
    tc = t[i] + step_s * -z[i] / (z[i + 1] - z[i])
    return (tc[-1] - tc[0]) / (len(tc) - 1)


def satrec_from_elements(
    epoch,
    *,
    altitude_km=None,
    repeat_cycle=None,
    eccentricity=0.0,
    inclination_deg=None,
    raan_deg=None,
    ltan_hours=None,
    arg_perigee_deg=0.0,
    mean_anomaly_deg=0.0,
    bstar=0.0,
    satnum=99999,
):
    """Build an SGP4 satellite from mean orbital elements.

    For a satellite with no TLE yet. The size of the orbit comes from exactly
    one of `altitude_km` or `repeat_cycle`, its orientation from exactly one of
    `raan_deg` or `ltan_hours`.

    Parameters
    ----------
    epoch : str | datetime
        UTC time the elements hold at.
    altitude_km : float, optional
        Mean altitude: mean semi-major axis minus the WGS84 equatorial radius.
        Fixes the period to ~1e-3 only — mission altitudes are rounded figures
        — which over days shifts the passes by minutes. Use `repeat_cycle`
        when the satellite is on a repeat ground track.
    repeat_cycle : (days, revolutions), optional
        Repeat ground track, e.g. ``(12, 175)`` for Sentinel-1: the nodal
        period is solved so the track repeats exactly, in SGP4 itself.
    eccentricity : float
    inclination_deg : float, optional
        ``None`` makes the orbit sun-synchronous at its altitude.
    raan_deg : float, optional
        Right ascension of the ascending node at the epoch (TEME).
    ltan_hours : float, optional
        Local mean solar time at the ascending node, e.g. 18.0 for
        Sentinel-1 — the natural way to set a sun-synchronous orbit.
    arg_perigee_deg, mean_anomaly_deg : float
        Position on the orbit at the epoch. For a circular orbit only their
        sum counts — the argument of latitude: 0 puts the satellite at the
        ascending node at the epoch. These set *when* the passes happen.
    bstar : float
        SGP4 drag term (1/earth radii); 0 for no drag.
    satnum : int

    Returns
    -------
    sgp4.api.Satrec
    """
    if (altitude_km is None) == (repeat_cycle is None):
        raise ValueError("give exactly one of altitude_km and repeat_cycle")
    if (raan_deg is None) == (ltan_hours is None):
        raise ValueError("give exactly one of raan_deg and ltan_hours")

    jd = sum(_julian(_utc(epoch)))[0]
    raan = (
        np.radians(raan_deg) if raan_deg is not None
        else _mean_sun_ra(jd) + np.radians((ltan_hours - 12.0) * 15.0)
    ) % (2 * np.pi)
    ecc = float(eccentricity)

    def inclination(a_km):
        return _sso_inclination(a_km, ecc) if inclination_deg is None else np.radians(inclination_deg)

    def build(n_rad_s, inc):
        sat = Satrec()
        sat.sgp4init(
            WGS72, "i", satnum, jd - JD_SGP4_EPOCH, bstar, 0.0, 0.0, ecc,
            np.radians(arg_perigee_deg), inc, np.radians(mean_anomaly_deg),
            n_rad_s * 60.0, raan,  # sgp4init wants rad/min
        )
        return sat

    if altitude_km is not None:
        a_target = WGS84_A / 1e3 + altitude_km
        n = np.sqrt(MU_WGS72 / a_target**3)
        # SGP4 reads its mean motion as Kozai's: rescale until its own mean
        # semi-major axis is the one asked for
        for _ in range(8):
            sat = build(n, inclination(a_target))
            ratio = sat.a * RE_WGS72 / a_target
            n *= ratio**1.5
            if abs(ratio - 1) < 1e-10:
                break
        return build(n, inclination(a_target))

    days, revs = repeat_cycle
    a = WGS84_A / 1e3 + 700.0
    n = np.sqrt(MU_WGS72 / a**3)
    for _ in range(12):
        inc = inclination(a)
        node_rate = SUN_RATE if inclination_deg is None else _node_rate(a, inc, ecc)
        # The track repeats when `revs` nodal periods last `days` turns of the
        # Earth under the orbit plane
        target = days * 2 * np.pi / (OMEGA_EARTH - node_rate) / revs
        sat = build(n, inc)
        measured = _nodal_period(sat)
        n *= measured / target
        a = (MU_WGS72 / n**2) ** (1 / 3)
        if abs(measured - target) < 1e-4:
            break
    return build(n, inclination(a))


# --- Ephemerides --------------------------------------------------------------

# Offset to add to a time in each scale to get UTC, valid since the 2017 leap
# second (TAI - UTC = 37 s)
TIME_SCALES = {"UTC": 0.0, "GPS": -18.0, "TAI": -37.0, "TT": -69.184}


def _frame_kind(name):
    """Map a frame name onto "ITRF" (Earth-fixed), "EME2000" or "TEME"."""
    key = str(name).upper().replace("-", "").replace("_", "")
    if key.startswith("ITRF") or key.startswith("IGS") or key in (
        "ECEF", "ECF", "EFG", "ITRS", "EARTHFIXED", "WGS84"
    ):
        return "ITRF"
    if key in ("EME2000", "J2000", "GCRF", "ICRF"):
        return "EME2000"
    if key == "TEME":
        return "TEME"
    raise ValueError(f"unsupported reference frame {name!r}: Earth-fixed, EME2000/GCRF or TEME")


def _lagrange(tk, yk, t, order, derivative=False):
    """Lagrange interpolation of yk(tk) at t on a sliding window of `order` points.

    tk (n,) and t (m,) in seconds; yk (n, 3). Returns (m, 3), or its time
    derivative with `derivative`.
    """
    i = np.searchsorted(tk, t)
    start = np.clip(i - order // 2, 0, len(tk) - order)
    idx = start[:, None] + np.arange(order)
    x, dt = tk[idx], t[:, None] - tk[idx]
    w = np.zeros_like(dt)
    for j in range(order):
        others = [k for k in range(order) if k != j]
        if not derivative:
            w[:, j] = np.prod([dt[:, k] / (x[:, j] - x[:, k]) for k in others], axis=0)
            continue
        for m in others:
            rest = [k for k in others if k != m]
            w[:, j] += np.prod([dt[:, k] / (x[:, j] - x[:, k]) for k in rest], axis=0) / (
                x[:, j] - x[:, m]
            )
    return np.einsum("mj,mjc->mc", w, yk[idx])


class Ephemeris:
    """An orbit given as state vectors, interpolated: the operator's ephemerides.

    Stored Earth-fixed, in m and m/s, velocity relative to the ground.
    Positions are interpolated with 8-point Lagrange polynomials — the CCSDS
    default, sub-metre at a 60 s sampling —, velocities likewise when given,
    otherwise taken from the derivative of the position polynomial.

    Parameters
    ----------
    times : array of datetime64 or str
        UTC times.
    positions : array (n, 3)
        Metres.
    velocities : array (n, 3), optional
        Metres per second, in `frame`.
    frame : str
        "ITRF" (any Earth-fixed name), "EME2000" (J2000, GCRF) or "TEME".
    name : str
    """

    ORDER = 8
    # A query between two samples further apart than this many median steps
    # falls in a gap of the ephemerides
    MAX_GAP = 5.0

    def __init__(self, times, positions, velocities=None, frame="ITRF", name=""):
        times = np.array([_utc(t) for t in times], dtype="datetime64[ns]")
        times, first = np.unique(times, return_index=True)  # sorted, duplicates dropped
        r = np.asarray(positions, float)[first]
        v = None if velocities is None else np.asarray(velocities, float)[first]
        if len(times) < self.ORDER:
            raise ValueError(f"{len(times)} state vectors: at least {self.ORDER} are needed")

        self.t0 = times[0]
        tk = (times - self.t0) / np.timedelta64(1, "s")
        kind = _frame_kind(frame)
        if v is None:
            # Differentiated in the frame the positions are given in: from
            # Earth-fixed positions, the velocity is relative to the ground
            v = _lagrange(tk, r, tk, self.ORDER, derivative=True)
        if kind != "ITRF":
            r, v = inertial_to_ecef(times, r, v, kind)

        self.times, self.tk, self.r, self.v = times, tk, r, v
        self.frame, self.name = frame, name
        self.step = float(np.median(np.diff(tk)))

    def __repr__(self):
        return (
            f"Ephemeris({self.name or 'unnamed'}, {len(self.times)} states, "
            f"{self.start:%Y-%m-%d %H:%M:%S} to {self.end:%Y-%m-%d %H:%M:%S} UTC, "
            f"step {self.step:g} s, from {self.frame})"
        )

    @property
    def start(self):
        return pd.Timestamp(self.times[0])

    @property
    def end(self):
        return pd.Timestamp(self.times[-1])

    def _seconds(self, times):
        return (np.asarray(times, dtype="datetime64[ns]") - self.t0) / np.timedelta64(1, "s")

    def covers(self, times):
        """Whether each time can be interpolated: inside the span, not in a gap."""
        t = self._seconds(times)
        i = np.clip(np.searchsorted(self.tk, t), 1, len(self.tk) - 1)
        gap = self.tk[i] - self.tk[i - 1]
        return (t >= -1e-3) & (t <= self.tk[-1] + 1e-3) & (gap <= self.MAX_GAP * self.step)

    def state(self, times):
        """Interpolated ECEF position (m) and ground-relative velocity (m/s), (n, 3)."""
        ok = self.covers(times)
        if not ok.all():
            at = pd.Timestamp(np.asarray(times, dtype="datetime64[ns]")[~ok][0])
            raise ValueError(
                f"{at} UTC is outside the ephemerides ({self.start} to {self.end} UTC) "
                "or in a gap of them"
            )
        t = self._seconds(times)
        return _lagrange(self.tk, self.r, t, self.ORDER), _lagrange(self.tk, self.v, t, self.ORDER)


def _parse_time(text, scale="UTC"):
    """An ephemeris time (ISO, or year-day-of-year 2026-278T...) as UTC datetime64."""
    text = text.strip().removeprefix("UTC=").rstrip("Z")
    doy = re.match(r"^(\d{4})-(\d{3})T(.*)$", text)
    if doy:
        t = np.datetime64(f"{doy[1]}-01-01T{doy[3]}", "ns") + np.timedelta64(int(doy[2]) - 1, "D")
    else:
        t = np.datetime64(text, "ns")
    return t + np.timedelta64(int(round(TIME_SCALES[scale.upper()] * 1e9)), "ns")


def _read_oem(text):
    """CCSDS OEM, KVN text form: km and km/s, frame and time system from the metadata."""
    times, states = [], []
    frame = scale = name = None
    in_meta = in_covariance = False
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("COMMENT"):
            continue
        if s in ("META_START", "META_STOP"):
            in_meta = s == "META_START"
            continue
        if s in ("COVARIANCE_START", "COVARIANCE_STOP"):
            in_covariance = s == "COVARIANCE_START"
            continue
        if in_covariance:
            continue
        if in_meta or "=" in s:
            key, _, value = (p.strip() for p in s.partition("="))
            if key == "REF_FRAME":
                if frame is not None and _frame_kind(value) != _frame_kind(frame):
                    raise ValueError("OEM segments in different frames are not supported")
                frame = value
            elif key == "TIME_SYSTEM":
                if value.upper() not in TIME_SCALES:
                    raise ValueError(f"unsupported OEM time system {value!r}")
                scale = value
            elif key == "OBJECT_NAME":
                name = value
            continue
        parts = s.split()
        if len(parts) >= 7:
            times.append(_parse_time(parts[0], scale or "UTC"))
            states.append([float(p) for p in parts[1:7]])
    if not states:
        raise ValueError("no state vector in the OEM file")
    states = np.array(states) * 1e3
    return Ephemeris(times, states[:, :3], states[:, 3:], frame or "EME2000", name or "")


def _read_eof(text):
    """ESA Earth Explorer orbit file (Sentinel-1 AUX_POEORB/RESORB/PREORB): Earth-fixed, m."""
    root = ET.fromstring(text)
    times, states = [], []
    for osv in root.iter("OSV"):
        times.append(_parse_time(osv.findtext("UTC")))
        states.append([float(osv.findtext(k)) for k in ("X", "Y", "Z", "VX", "VY", "VZ")])
    if not states:
        raise ValueError("no OSV element in the orbit file")
    states = np.array(states)
    name = root.findtext(".//Mission") or ""
    return Ephemeris(times, states[:, :3], states[:, 3:], "ITRF", name)


def _read_table(text, frame, time_scale):
    """CSV/whitespace table: a time column, x y z, optionally vx vy vz; m or km guessed."""
    from io import StringIO

    table = pd.read_csv(StringIO(text), sep=None, engine="python", comment="#")
    table.columns = [c.strip().lower() for c in table.columns]
    time_col = next((c for c in table.columns if c in ("time", "utc", "epoch", "date", "datetime")), None)
    if time_col is None or not {"x", "y", "z"} <= set(table.columns):
        raise ValueError("a table needs a time (or utc/epoch) column and x, y, z columns")
    r = table[["x", "y", "z"]].to_numpy(float)
    v = table[["vx", "vy", "vz"]].to_numpy(float) if {"vx", "vy", "vz"} <= set(table.columns) else None
    # An orbit radius is ~7000 km: the magnitude tells the unit
    scale = 1.0 if np.median(np.linalg.norm(r, axis=1)) > 1e5 else 1e3
    times = [_parse_time(str(t), time_scale) for t in table[time_col]]
    return Ephemeris(times, r * scale, None if v is None else v * scale, frame, "")


def read_ephemeris(path, frame="ITRF", time_scale="UTC"):
    """Read ephemerides from a file, its format told from its content.

    Parameters
    ----------
    path : str | Path
        - CCSDS OEM (KVN text): frame, time system and units from the file;
        - ESA orbit file (``.EOF``, Sentinel-1 POEORB/RESORB/PREORB);
        - a table (CSV, or whitespace-separated): a ``time`` column, ``x y z``
          and optionally ``vx vy vz``, in m or km (told from the magnitude).
    frame : str
        Frame of a table: "ITRF" (Earth-fixed), "EME2000" or "TEME". Files
        that state their frame ignore it.
    time_scale : {"UTC", "GPS", "TAI", "TT"}
        Time scale of a table.

    Returns
    -------
    Ephemeris
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    head = text.lstrip()[:2000]
    if "CCSDS_OEM_VERS" in head:
        if head.startswith("<"):
            raise ValueError("OEM in XML form is not supported: ask for the KVN (text) form")
        return _read_oem(text)
    if head.startswith("<"):
        return _read_eof(text)
    return _read_table(text, frame, time_scale)


def satellite(tle=None, norad_id=None, elements=None, ephemeris=None):
    """Build the satellite from exactly one orbit source.

    Parameters
    ----------
    tle : str, optional
        Two lines, or three with a name first.
    norad_id : int, optional
        The latest TLE is fetched from CelesTrak.
    elements : dict, optional
        Keyword arguments of `satrec_from_elements`.
    ephemeris : str | Path | Ephemeris, optional
        A file `read_ephemeris` reads, or an `Ephemeris` already built.

    Returns
    -------
    sgp4.api.Satrec | Ephemeris
    """
    if sum(x is not None for x in (tle, norad_id, elements, ephemeris)) != 1:
        raise ValueError("give exactly one of tle, norad_id, elements and ephemeris")
    if ephemeris is not None:
        return ephemeris if isinstance(ephemeris, Ephemeris) else read_ephemeris(ephemeris)
    if elements is not None:
        return satrec_from_elements(**elements)
    _, l1, l2 = fetch_tle(norad_id) if norad_id is not None else parse_tle(tle)
    return Satrec.twoline2rv(l1, l2)


def epoch_of(sat):
    """Epoch of a satellite's elements — start of ephemerides —, as a UTC Timestamp."""
    if isinstance(sat, Ephemeris):
        return sat.start
    days = sat.jdsatepoch + sat.jdsatepochF - JD_UNIX_EPOCH
    return pd.Timestamp(np.datetime64(round(days * 86400e6), "us"))


def _ephemeris_summary(eph):
    """`orbit_summary` of ephemerides, measured on the state vectors themselves."""
    r, v = eph.r, eph.v
    _, _, h = ecef_to_geodetic(r)
    # Back to an inertial velocity for the orbit plane: v + omega x r
    v_in = v + OMEGA_EARTH * np.column_stack([-r[:, 1], r[:, 0], np.zeros(len(r))])
    mom = np.cross(r, v_in)
    inclination = np.degrees(np.arccos(mom[:, 2] / np.linalg.norm(mom, axis=1))).mean()

    z = r[:, 2]
    i = np.flatnonzero((z[:-1] < 0) & (z[1:] >= 0))
    frac = -z[i] / (z[i + 1] - z[i])
    t_node = eph.tk[i] + frac * (eph.tk[i + 1] - eph.tk[i])
    period = np.nan
    if len(t_node) > 1:
        # Count the revolutions rather than the nodes: a gap may hide some
        revs = np.round((t_node[-1] - t_node[0]) / np.median(np.diff(t_node)))
        period = (t_node[-1] - t_node[0]) / revs
    ltan = np.nan
    if len(t_node):
        node = r[i] + frac[:, None] * (r[i + 1] - r[i])
        lon = np.degrees(np.arctan2(node[:, 1], node[:, 0]))
        hours = ((eph.t0 - eph.t0.astype("datetime64[D]")) / np.timedelta64(1, "h") + t_node / 3600.0)
        angle = np.radians((hours + lon / 15.0) * 15.0)  # mean solar time, as an angle
        ltan = np.degrees(np.arctan2(np.sin(angle).mean(), np.cos(angle).mean())) / 15.0 % 24.0
    return {
        "epoch": eph.start,
        "end": eph.end,
        "mean_altitude_km": h.mean() / 1e3,
        "inclination_deg": inclination,
        "nodal_period_min": period / 60.0,
        "revs_per_day": 86400.0 / period,
        "ltan_hours": ltan,
    }


def orbit_summary(sat):
    """A few readable figures of an orbit, for a sanity check.

    From elements, the altitude is the mean semi-major axis minus the
    equatorial radius; from ephemerides, the mean height above the ellipsoid
    (~10 km more on a polar orbit, the Earth being flattened).
    """
    if isinstance(sat, Ephemeris):
        return _ephemeris_summary(sat)
    jd = sat.jdsatepoch + sat.jdsatepochF
    ltan = (12.0 + np.degrees(sat.nodeo - _mean_sun_ra(jd)) / 15.0) % 24.0
    period = _nodal_period(sat)
    return {
        "epoch": epoch_of(sat),
        "mean_altitude_km": sat.a * RE_WGS72 - WGS84_A / 1e3,
        "inclination_deg": np.degrees(sat.inclo),
        "eccentricity": sat.ecco,
        "nodal_period_min": period / 60.0,
        "revs_per_day": 86400.0 / period,
        "ltan_hours": ltan,
    }


# --- Passes -------------------------------------------------------------------

def _doppler(sat, t_ns, targets):
    """(r - target) . v at float ns times; zero at closest approach (zero-Doppler)."""
    r, v = propagate(sat, np.round(t_ns).astype(np.int64).astype("datetime64[ns]"))
    return _dot(r - targets, v)


def _refine_zero_doppler(sat, targets, ta, tb, iterations=4):
    """Zero-Doppler time of each target in [ta, tb], where f goes from - to +.

    Secant steps on f = (r - T) . v, nearly linear over a few seconds.
    """
    a = ta.astype("datetime64[ns]").astype(np.int64).astype(float)
    b = tb.astype("datetime64[ns]").astype(np.int64).astype(float)
    t0, t1 = a, b
    f0, f1 = _doppler(sat, a, targets), _doppler(sat, b, targets)
    for _ in range(iterations):
        denom = np.where(f1 != f0, f1 - f0, 1.0)
        t2 = np.clip(np.where(f1 != f0, t1 - f1 * (t1 - t0) / denom, t1), a, b)
        t0, f0 = t1, f1
        t1, f1 = t2, _doppler(sat, t2, targets)
    return np.round(t1).astype(np.int64).astype("datetime64[ns]")


def _geometry(r, v, targets, ups):
    """Incidence (deg), side, slant range (m), off-nadir (deg) per satellite/target pair."""
    los = r - targets  # ground to satellite
    rng = np.linalg.norm(los, axis=1)
    incidence = np.degrees(np.arccos(np.clip(_dot(los, ups) / rng, -1, 1)))
    right = np.cross(v, r)  # velocity x up points to the right of the track
    side = np.where(_dot(-los, right) > 0, "right", "left")
    off_nadir = np.degrees(
        np.arccos(np.clip(_dot(-r, -los) / (np.linalg.norm(r, axis=1) * rng), -1, 1))
    )
    return incidence, side, rng, off_nadir


def _heading(lon, lat, v):
    """Ground-track heading (deg from north, clockwise) of an ECEF velocity."""
    lo, la = np.radians(lon), np.radians(lat)
    east = np.column_stack([-np.sin(lo), np.cos(lo), np.zeros_like(lo)])
    north = np.column_stack([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    return np.degrees(np.arctan2(_dot(v, east), _dot(v, north))) % 360.0


def find_passes(
    sat,
    aoi,
    start,
    end,
    *,
    look_side="right",
    incidence_range=DEFAULT_INCIDENCE_RANGE,
    directions=DIRECTIONS,
    height_m=0.0,
    only_feasible=True,
    step_s=10.0,
    spacing_km=5.0,
):
    """List the passes over an AOI and the SAR geometry of each.

    A pass is timed when the AOI's reference point is at zero-Doppler — the
    closest approach, when a SAR images it. Every boundary point of the AOI is
    then solved for its own zero-Doppler time, which gives the incidence and
    side it is seen at.

    Parameters
    ----------
    sat : sgp4.api.Satrec | Ephemeris
        From `satellite`.
    aoi : anything `parse_aoi` accepts
    start, end : str | datetime
        UTC window; with ephemerides, clipped to their span.
    look_side : {"right", "left", "both"}
        Side the radar looks to; "both" for a satellite that can roll.
    incidence_range : (near, far) in degrees, or None
        Usable incidence range of the instrument; None for no limit.
    directions : iterable of {"ascending", "descending"}
    height_m : float
        Ellipsoidal height of the AOI.
    only_feasible : bool
        Keep only the passes that can image part of the AOI; False keeps
        every pass with the satellite above the AOI's horizon.
    step_s : float
        Time step of the coarse search.
    spacing_km : float
        Spacing of the AOI boundary points.

    Returns
    -------
    pandas.DataFrame
        One row per pass, `PASS_COLUMNS`. ``side`` is the side of the track
        the AOI lies on ("both" when the track crosses it); ``coverage`` the
        share of the AOI area imageable — on an allowed side, within
        `incidence_range` —, counted on a grid of `spacing_km`;
        ``incidence_ref``, ``off_nadir_ref`` and ``slant_range_km`` are at the
        reference point; ``local_time`` is the
        mean solar time there. ``attrs`` holds the window searched — clipped
        to the span of ephemerides, whose gaps are skipped — and the gap
        between the window and the epoch of the elements.
    """
    if look_side not in LOOK_SIDES:
        raise ValueError(f"look_side must be one of {LOOK_SIDES}")
    directions = (directions,) if isinstance(directions, str) else tuple(directions)
    if not set(directions) <= set(DIRECTIONS):
        raise ValueError(f"directions must be among {DIRECTIONS}")
    near, far = (0.0, 90.0) if incidence_range is None else incidence_range

    lonlat, n_area = aoi_points(aoi, spacing_km)
    # Points the covered share of the area is counted on
    area = slice(1, n_area + 1) if n_area else slice(0, 1)
    targets = geodetic_to_ecef(lonlat[:, 0], lonlat[:, 1], height_m)
    ups = _up(lonlat[:, 0], lonlat[:, 1])
    ref, ref_up = targets[:1], ups[:1]

    # Every AOI point reaches zero-Doppler within its along-track distance
    # from the reference point, at ~7 km/s
    extent = np.linalg.norm(targets - ref, axis=1).max()
    half = extent / 6000.0 + 20.0
    offsets = (np.arange(-half, half + 5.0, 5.0) * 1e9).astype("timedelta64[ns]")

    t0, t1 = _utc(start), _utc(end)
    if t1 <= t0:
        raise ValueError("end must be after start")
    if isinstance(sat, Ephemeris):
        # Leave room for the per-pass search window at both ends
        margin = np.timedelta64(int((half + 10.0) * 1e9), "ns")
        t0, t1 = max(t0, sat.times[0] + margin), min(t1, sat.times[-1] - margin)
        if t1 <= t0:
            raise ValueError(f"window outside the ephemerides, which run {sat.start} to {sat.end} UTC")
    step = np.timedelta64(int(step_s * 1e9), "ns")
    times = np.arange(t0, t1 + step, step)
    if isinstance(sat, Ephemeris):
        times = times[sat.covers(times)]  # no search across a gap

    # Closest approaches of the reference point, kept above its horizon
    r, v = propagate(sat, times)
    f = _dot(r - ref, v)
    i = np.flatnonzero((f[:-1] < 0) & (f[1:] >= 0) & (np.diff(times) == step))
    i = i[_dot(r[i] - ref, np.repeat(ref_up, len(i), 0)) > 0]
    t_ref = _refine_zero_doppler(sat, np.repeat(ref, len(i), 0), times[i], times[i + 1])

    rows = []
    for tc in t_ref:
        local = tc + offsets
        if isinstance(sat, Ephemeris) and not sat.covers(local).all():
            continue  # pass too close to a gap to be solved over the whole AOI
        r, v = propagate(sat, local)
        f = _dot(r, v)[:, None] - v @ targets.T  # (time, target)
        cross = (f[:-1] < 0) & (f[1:] >= 0)
        found = cross.any(axis=0)
        k = cross.argmax(axis=0)[found]
        tz = _refine_zero_doppler(sat, targets[found], local[k], local[k + 1])
        rz, vz = propagate(sat, tz)
        inc, side, _, _ = _geometry(rz, vz, targets[found], ups[found])

        allowed = np.ones(len(inc), bool) if look_side == "both" else side == look_side
        ok = np.zeros(len(targets), bool)
        ok[found] = allowed & (inc >= near) & (inc <= far) & (inc < 90.0)
        coverage = ok[area].mean()

        # Reference point (row 0) is always found: its time is tc itself
        r0, v0 = propagate(sat, np.array([tc]))
        inc0, _, rng0, off0 = _geometry(r0, v0, ref, ref_up)
        lon, lat, alt = ecef_to_geodetic(r0)
        sides = set(side)
        hours = (tc - tc.astype("datetime64[D]")) / np.timedelta64(1, "h")
        solar = (hours + lonlat[0, 0] / 15.0) % 24.0
        rows.append({
            "time_utc": pd.Timestamp(tc),
            "direction": "ascending" if v0[0, 2] > 0 else "descending",
            "side": sides.pop() if len(sides) == 1 else "both",
            "incidence_ref": inc0[0],
            "incidence_min": inc.min(),
            "incidence_max": inc.max(),
            "coverage": coverage,
            "off_nadir_ref": off0[0],
            "slant_range_km": rng0[0] / 1e3,
            "local_time": f"{int(solar):02d}:{int(solar % 1 * 60):02d}",
            "sat_lon": lon[0],
            "sat_lat": lat[0],
            "sat_alt_km": alt[0] / 1e3,
            "heading_deg": _heading(lon, lat, v0)[0],
        })

    passes = pd.DataFrame(rows, columns=PASS_COLUMNS)
    passes = passes[passes.direction.isin(directions)]
    if only_feasible:
        passes = passes[passes.coverage > 0]
    passes = passes.reset_index(drop=True)

    epoch = epoch_of(sat)
    gap = max(abs(pd.Timestamp(t0) - epoch), abs(pd.Timestamp(t1) - epoch))
    passes.attrs.update({
        "window": (pd.Timestamp(t0), pd.Timestamp(t1)),  # clipped to ephemerides
        "epoch": epoch,
        "epoch_gap_days": gap / pd.Timedelta(days=1),
        "look_side": look_side,
        "incidence_range": incidence_range,
        "aoi_reference": tuple(lonlat[0]),
    })
    return passes


# --- Swath on the ground ------------------------------------------------------

def _ray_to_ellipsoid(r, d, height=0.0):
    """First hit of rays r + t d (ECEF, m) on the ellipsoid raised by `height`; NaN if none."""
    a, b = WGS84_A + height, WGS84_B + height
    scale = np.array([1.0, 1.0, a / b])
    p, q = r * scale, d * scale
    qa, qb, qc = _dot(q, q), 2 * _dot(p, q), _dot(p, p) - a**2
    disc = qb**2 - 4 * qa * qc
    t = (-qb - np.sqrt(np.where(disc >= 0, disc, np.nan))) / (2 * qa)
    hit = r + t[:, None] * d
    normal = _unit(hit / np.array([a**2, a**2, b**2]))
    return hit, normal


def _ground_at_incidence(r, v, side, incidence_deg, height=0.0, iterations=6):
    """Ground point seen at a given incidence, in the zero-Doppler plane, (n, 3) ECEF."""
    up = _unit(r)
    vh = _unit(v)
    # Zero-Doppler plane: perpendicular to the ground-relative velocity
    e1 = _unit(-up - _dot(-up, vh)[:, None] * vh)
    e2 = _unit(np.cross(v, r)) * (1.0 if side == "right" else -1.0)
    target = np.radians(incidence_deg)

    def incidence(alpha):
        d = np.cos(alpha)[:, None] * e1 + np.sin(alpha)[:, None] * e2
        hit, normal = _ray_to_ellipsoid(r, d, height)
        return np.arccos(np.clip(_dot(normal, _unit(r - hit)), -1, 1)), hit

    # Spherical first guess, then secant steps on the ellipsoid
    rn = np.linalg.norm(r, axis=1)
    a0 = np.arcsin(np.clip(WGS84_A / rn * np.sin(target), -1, 1))
    a1 = a0 * 1.01
    g0 = incidence(a0)[0] - target
    g1 = incidence(a1)[0] - target
    for _ in range(iterations):
        denom = np.where(g1 != g0, g1 - g0, 1.0)
        a2 = np.where(g1 != g0, a1 - g1 * (a1 - a0) / denom, a1)
        a0, g0 = a1, g1
        a1 = a2
        g1 = incidence(a1)[0] - target
    return incidence(a1)[1]


def swath_footprint(
    sat,
    time,
    side="right",
    incidence_range=DEFAULT_INCIDENCE_RANGE,
    half_length_s=60.0,
    step_s=5.0,
    height_m=0.0,
    center_lon=None,
):
    """Ground strip imaged between two incidence angles, around a time.

    Parameters
    ----------
    sat : sgp4.api.Satrec
    time : str | datetime
        Centre of the strip, typically a pass's ``time_utc``.
    side : {"right", "left"}
    incidence_range : (near, far) in degrees
    half_length_s : float
        Half the strip length in seconds of flight (~7 km/s).
    step_s : float
    height_m : float
    center_lon : float, optional
        Longitude to draw near, typically the AOI's: the strip is shifted by
        whole turns to lie within 180° of it.

    Returns
    -------
    swath : shapely.geometry.Polygon
        The strip, lon/lat. Longitudes are continuous along it, so they may
        leave [-180, 180] across the antimeridian.
    track : shapely.geometry.LineString
        The nadir track over the same time span, lon/lat, same convention.
    """
    tc = _utc(time)
    offsets = (np.arange(-half_length_s, half_length_s + step_s, step_s) * 1e9).astype("timedelta64[ns]")
    r, v = propagate(sat, tc + offsets)
    near = _ground_at_incidence(r, v, side, incidence_range[0], height_m)
    far = _ground_at_incidence(r, v, side, incidence_range[1], height_m)
    ring = np.vstack([
        np.column_stack(ecef_to_geodetic(near)[:2]),
        np.column_stack(ecef_to_geodetic(far)[:2])[::-1],
    ])
    track = np.column_stack(ecef_to_geodetic(r)[:2])
    # No 360° jump along the strip, and the track on the same turn as the
    # strip: its first point is the near edge's first point's neighbour
    ring[:, 0] = np.degrees(np.unwrap(np.radians(ring[:, 0])))
    track[:, 0] = np.degrees(np.unwrap(np.radians(track[:, 0])))
    track[:, 0] += 360.0 * np.round((ring[0, 0] - track[0, 0]) / 360.0)
    if center_lon is not None:
        shift = 360.0 * np.round((center_lon - ring[:, 0].mean()) / 360.0)
        ring[:, 0] += shift
        track[:, 0] += shift
    return Polygon(ring), LineString(track)


# --- Notebook interface -------------------------------------------------------

def _pass_label(row):
    return (
        f"{row.time_utc:%Y-%m-%d %H:%M:%S}  {row.direction[:4].replace('asce', 'asc '):4}  {row.side:5}  "
        f"inc {row.incidence_min:5.1f}-{row.incidence_max:5.1f}°  "
        f"cover {row.coverage:4.0%}  LT {row.local_time}"
    )


def _optional_float(text):
    text = str(text).strip()
    return float(text) if text else None


class TaskingUI:
    """Webmap interface: AOI, orbit, sensor, window — then passes and their swath.

    Everything is also reachable from Python: `aoi` (shapely geometry), `sat`
    (the last orbit used, a Satrec or an Ephemeris) and `passes` (the last
    result, a DataFrame).
    """

    def __init__(self, center=(46.5, 2.5), zoom=6, start=None, days=12):
        # Imported here so that the library functions above work without the
        # widget packages installed
        import ipywidgets as w
        from ipyleaflet import DrawControl, GeoJSON, LayersControl, Map, basemaps

        self.aoi = None
        self.sat = None
        self.passes = None
        empty = {"type": "FeatureCollection", "features": []}

        # --- map ---
        self.map = Map(
            center=center, zoom=zoom, basemap=basemaps.OpenStreetMap.Mapnik,
            scroll_wheel_zoom=True, layout=w.Layout(height="480px", width="100%"),
        )
        aoi_style = {"color": "#d62728", "weight": 2, "fillOpacity": 0.05}
        self.draw = DrawControl(
            polygon={"shapeOptions": aoi_style}, rectangle={"shapeOptions": aoi_style},
            circle={}, circlemarker={}, marker={}, polyline={}, edit=False, remove=False,
        )
        self.draw.on_draw(self._on_draw)
        self.map.add(self.draw)
        self.swath_layer = GeoJSON(
            data=empty, name="Swath", style={"color": "#1f77b4", "weight": 2, "fillOpacity": 0.15}
        )
        self.track_layer = GeoJSON(
            data=empty, name="Nadir track", style={"color": "#555555", "weight": 2, "dashArray": "6"}
        )
        self.aoi_layer = GeoJSON(data=empty, name="AOI", style=aoi_style)
        for layer in (self.swath_layer, self.track_layer, self.aoi_layer):
            self.map.add(layer)
        self.map.add(LayersControl(position="topright"))

        # --- AOI text ---
        self.aoi_text = w.Textarea(
            placeholder="WKT or GeoJSON, lon/lat — or draw on the map",
            layout=w.Layout(width="100%", height="60px"),
        )
        self.use_text_btn = w.Button(description="Use this AOI", icon="check")
        self.use_text_btn.on_click(lambda _: self._guard(self._use_text))

        # --- orbit ---
        self.source = w.ToggleButtons(
            options=["TLE", "NORAD id", "Elements", "Ephemeris"], value="TLE"
        )
        self.tle_text = w.Textarea(
            placeholder="TLE: two lines, or three with the name first",
            layout=w.Layout(width="100%", height="70px"),
        )
        self.norad = w.IntText(description="NORAD id", value=39634)
        text = lambda desc, value="", ph="": w.Text(  # noqa: E731
            description=desc, value=value, placeholder=ph, style={"description_width": "130px"}
        )
        self.el = {
            "epoch": text("Epoch (UTC)", ph="2026-10-05T00:00:00"),
            "altitude_km": text("Altitude (km)", ph="or repeat cycle"),
            "repeat_days": text("Repeat days", ph="e.g. 12"),
            "repeat_revs": text("Repeat revs", ph="e.g. 175"),
            "inclination_deg": text("Inclination (°)", ph="blank: sun-synchronous"),
            "ltan_hours": text("LTAN (h)", ph="e.g. 18 — or RAAN"),
            "raan_deg": text("RAAN (°)", ph="or LTAN"),
            "eccentricity": text("Eccentricity", "0"),
            "arg_latitude_deg": text("Arg. latitude (°)", "0", "0: at asc. node at epoch"),
        }
        self.eph_path = w.Text(
            description="File", placeholder="OEM, Sentinel-1 .EOF, or CSV (time, x, y, z[, vx, vy, vz])",
            layout=w.Layout(width="70%"),
        )
        self.eph_frame = w.Dropdown(
            description="CSV frame", options=["ITRF", "EME2000", "TEME"], value="ITRF",
            tooltip="Frame of a CSV table; OEM and EOF files state their own",
        )
        self.orbit_stack = w.Stack(
            [self.tle_text, self.norad, w.GridBox(
                list(self.el.values()),
                layout=w.Layout(grid_template_columns="repeat(3, 1fr)"),
            ), w.HBox([self.eph_path, self.eph_frame])],
            selected_index=0,
        )
        w.jslink((self.source, "index"), (self.orbit_stack, "selected_index"))

        # --- sensor and window ---
        self.look_side = w.Dropdown(description="Look side", options=list(LOOK_SIDES), value="right")
        self.near = w.FloatText(description="Near inc. (°)", value=DEFAULT_INCIDENCE_RANGE[0])
        self.far = w.FloatText(description="Far inc. (°)", value=DEFAULT_INCIDENCE_RANGE[1])
        self.start = w.DatePicker(description="From", value=start or pd.Timestamp.now("UTC").date())
        self.days = w.IntText(description="Days", value=days)
        self.direction = w.Dropdown(description="Direction", options=["any", *DIRECTIONS], value="any")
        self.only_feasible = w.Checkbox(description="Imageable passes only", value=True, indent=False)
        self.run_btn = w.Button(description="Find passes", icon="search", button_style="primary")
        self.run_btn.on_click(lambda _: self._guard(self._run))
        self.status = w.HTML("Draw an AOI on the map or paste one below, then set the orbit.")

        # --- results ---
        self.list = w.Select(options=[], rows=12, layout=w.Layout(width="100%"))
        self.list.observe(lambda _: self._guard(self._on_select), names="value")
        self.details = w.HTML("&nbsp;")

        self.widget = w.VBox([
            self.map,
            w.HBox([self.aoi_text, self.use_text_btn]),
            self.source,
            self.orbit_stack,
            w.HBox([self.look_side, self.near, self.far]),
            w.HBox([self.start, self.days, self.direction]),
            w.HBox([self.only_feasible, self.run_btn]),
            self.status,
            self.list,
            self.details,
        ])

    def _ipython_display_(self):
        from IPython.display import display
        display(self.widget)

    # --- helpers ---

    def _guard(self, action):
        """Run a handler and show its error: widget callbacks swallow exceptions."""
        try:
            action()
        except Exception as exc:  # noqa: BLE001 — anything must reach the user
            self.status.value = f'<span style="color:#d62728">{type(exc).__name__}: {exc}</span>'

    def _set_aoi(self, geom, fit=True):
        # Shown as computed: an AOI across the antimeridian stays in one piece
        geom = unwrap_antimeridian(geom)
        self.aoi = geom
        self.aoi_layer.data = {
            "type": "FeatureCollection",
            "features": [{"type": "Feature", "properties": {}, "geometry": mapping(geom)}],
        }
        if fit:
            minx, miny, maxx, maxy = geom.bounds
            self.map.fit_bounds([[miny, minx], [maxy, maxx]])

    def _elements(self):
        el = {k: wdg.value.strip() for k, wdg in self.el.items()}
        if not el["epoch"]:
            raise ValueError("elements need an epoch")
        repeat = None
        if el["repeat_days"] or el["repeat_revs"]:
            repeat = (float(el["repeat_days"]), int(el["repeat_revs"]))
        return {
            "epoch": el["epoch"],
            "altitude_km": _optional_float(el["altitude_km"]),
            "repeat_cycle": repeat,
            "eccentricity": _optional_float(el["eccentricity"]) or 0.0,
            "inclination_deg": _optional_float(el["inclination_deg"]),
            "raan_deg": _optional_float(el["raan_deg"]),
            "ltan_hours": _optional_float(el["ltan_hours"]),
            "mean_anomaly_deg": _optional_float(el["arg_latitude_deg"]) or 0.0,
        }

    def _satellite(self):
        if self.source.value == "TLE":
            return satellite(tle=self.tle_text.value)
        if self.source.value == "NORAD id":
            return satellite(norad_id=self.norad.value)
        if self.source.value == "Elements":
            return satellite(elements=self._elements())
        path = self.eph_path.value.strip().strip('"')
        if not path:
            raise ValueError("no ephemeris file given")
        return read_ephemeris(path, frame=self.eph_frame.value)

    # --- handlers ---

    def _on_draw(self, control, action, geo_json):
        if action != "created":
            return
        geom = shape(geo_json["geometry"])
        control.clear()  # keep a single AOI, shown through aoi_layer
        self._set_aoi(geom, fit=False)
        self.aoi_text.value = geom.wkt
        self.status.value = f"AOI set from the map ({geom.geom_type})."

    def _use_text(self):
        geom = parse_aoi(self.aoi_text.value)
        self._set_aoi(geom)
        self.status.value = f"AOI set from text ({geom.geom_type})."

    def _run(self):
        if self.aoi is None:
            raise ValueError("no AOI: draw one on the map or paste one and click 'Use this AOI'")
        self.status.value = "Computing…"
        self.sat = self._satellite()
        start = pd.Timestamp(self.start.value)
        self.passes = find_passes(
            self.sat, self.aoi, start, start + pd.Timedelta(days=self.days.value),
            look_side=self.look_side.value,
            incidence_range=(self.near.value, self.far.value),
            directions=DIRECTIONS if self.direction.value == "any" else self.direction.value,
            only_feasible=self.only_feasible.value,
        )
        self.list.options = [
            (_pass_label(row), i) for i, row in enumerate(self.passes.itertuples(index=False))
        ]
        s = orbit_summary(self.sat)
        gap = self.passes.attrs["epoch_gap_days"]
        warn = (
            f' <span style="color:#d62728">— window up to {gap:.1f} days from the epoch: '
            "expect km-level along-track error</span>"
            if gap > TLE_MAX_AGE_DAYS and self.source.value in ("TLE", "NORAD id") else ""
        )
        if isinstance(self.sat, Ephemeris):
            w0, w1 = self.passes.attrs["window"]
            span = f"ephemerides {s['epoch']:%Y-%m-%d %H:%M} to {s['end']:%Y-%m-%d %H:%M}, " \
                   f"searched {w0:%Y-%m-%d %H:%M} to {w1:%Y-%m-%d %H:%M}"
        else:
            span = f"epoch {s['epoch']:%Y-%m-%d %H:%M}"
        self.status.value = (
            f"{len(self.passes)} pass(es). Orbit: {span}, "
            f"alt {s['mean_altitude_km']:.1f} km, inc {s['inclination_deg']:.2f}°, "
            f"{s['revs_per_day']:.4f} rev/day, LTAN {s['ltan_hours']:.2f} h{warn}"
        )

    def _on_select(self):
        if self.passes is None or self.list.value is None:
            return
        row = self.passes.iloc[self.list.value]
        sides = ["right", "left"] if row.side == "both" else [row.side]
        if self.look_side.value != "both":
            sides = [self.look_side.value]
        features = []
        for side in sides:
            swath, track = swath_footprint(
                self.sat, row.time_utc, side, (self.near.value, self.far.value),
                center_lon=self.aoi.representative_point().x,
            )
            features.append({"type": "Feature", "properties": {"side": side}, "geometry": mapping(swath)})
        self.swath_layer.data = {"type": "FeatureCollection", "features": features}
        self.track_layer.data = {
            "type": "FeatureCollection",
            "features": [{"type": "Feature", "properties": {}, "geometry": mapping(track)}],
        }
        self.details.value = (
            f"<b>{row.time_utc:%Y-%m-%d %H:%M:%S} UTC</b> — {row.direction}, AOI on the {row.side} "
            f"side. Incidence {row.incidence_ref:.2f}° at the reference point "
            f"({row.incidence_min:.2f}–{row.incidence_max:.2f}° over the AOI), off-nadir "
            f"{row.off_nadir_ref:.2f}°, slant range {row.slant_range_km:.1f} km, heading "
            f"{row.heading_deg:.1f}°, coverage {row.coverage:.0%}."
        )


def build_ui(**kwargs):
    """Build the webmap interface; see `TaskingUI` for the keyword arguments."""
    return TaskingUI(**kwargs)
