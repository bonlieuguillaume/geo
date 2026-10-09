# pie_echo_wall — SAR tasking

When can a side-looking radar image an area of interest, and at which
incidence? From a satellite's orbit, list the passes over the AOI with their
zero-Doppler time, orbit direction, the side of the track the AOI lies on, the
incidence over the AOI and the share of it inside the instrument's swath; on
the map, the imaged strip of any pass.

- `tasking.py` — the logic: orbit (`satellite`, `satrec_from_elements`,
  `fetch_tle`, `read_ephemeris`, `orbit_summary`), passes (`find_passes`), swath on the ground
  (`swath_footprint`), plus a notebook interface (`build_ui`, ipyleaflet +
  ipywidgets). Usable as a library without the widget packages.
- `tasking.ipynb` — the webmap interface, plus the calls from plain Python. A
  thin driver of the module, *not* a mirror of it.
- `test_tasking.py` — the same as a script, every parameter at the top.

## Orbit input

Exactly one of:

| Source | When | Accuracy |
| --- | --- | --- |
| **TLE** (2 lines, or 3 with a name) | satellite flying | km-level, growing a few km/day away from the TLE epoch: keep the window within ~3 days of it |
| **NORAD id** | same, latest TLE fetched from CelesTrak (anonymous) | same |
| **Mean elements** | satellite not flying yet | as good as the elements: see below |
| **Ephemerides** | the operator's predicted orbit | tens to hundreds of metres over days; the window is clipped to the file's span |

### Ephemerides

The best input for a flying satellite: the operator's flight dynamics fit the
orbit to onboard GNSS and propagate it with full force models, **including
the manoeuvres planned when the file was made** (not later ones, such as an
unplanned collision avoidance: refresh the file as the operator does, often
daily). Not to be confused with the *reference* orbit, the ideal track a
mission is kept to, nor with *restituted* orbits (Sentinel-1 POEORB, cm-level
but after the fact — useful to check the tool, not to plan).

`read_ephemeris` tells the format from the content:

| Format | Frame | Units / time |
| --- | --- | --- |
| CCSDS OEM, KVN text form (not XML) | from `REF_FRAME`: Earth-fixed (ITRF…), EME2000/J2000/GCRF, TEME | km, km/s; `TIME_SYSTEM` UTC, GPS, TAI or TT |
| ESA orbit file `.EOF` (Sentinel-1 POEORB, RESORB, PREORB) | Earth-fixed | m, m/s, UTC |
| Table, CSV or whitespace: `time, x, y, z[, vx, vy, vz]` | `frame` argument | m or km (told from the magnitude); `time_scale` argument |

States are interpolated with 8-point Lagrange polynomials (CCSDS default):
millimetre-level at a 60 s sampling. Velocities missing from a table come
from the derivative of the position polynomial. A query in a gap of more than
5 sampling steps is refused, and passes too close to a gap are skipped. A
manoeuvre inside the file is smoothed over a few samples.

EME2000 goes to Earth-fixed through IAU 1976 precession and the main IAU 1980
nutation terms: checked against Vallado's example 3-15, 16 m from ITRF with
UT1−UTC applied — UT1−UTC itself, ignored here, adds up to ~50 m today.

Mean elements (`satrec_from_elements`):

- **epoch** — UTC time they hold at.
- size: **`altitude_km`** (mean semi-major axis minus 6378.137 km) **or
  `repeat_cycle=(days, revs)`**. Prefer the repeat cycle for a repeat-track
  mission (Sentinel-1: `(12, 175)`): the nodal period is solved so the track
  repeats exactly. A rounded altitude puts the period off by ~1e-3, minutes of
  drift in the pass times after a few days.
- orientation: **`ltan_hours`** (local time at the ascending node, 18.0 for
  Sentinel-1) **or `raan_deg`**.
- **`inclination_deg`** — leave `None` for a sun-synchronous orbit, computed
  from the altitude.
- **`mean_anomaly_deg`** (+ `arg_perigee_deg`) — where the satellite is on its
  orbit at the epoch; for a circular orbit their sum is the argument of
  latitude, 0 = at the ascending node. This sets *when* the passes happen; with
  a wrong value, the geometry of the passes stays right but their dates and
  times do not.
- `eccentricity` (0), `bstar` (0, no drag).

## Sensor

- **`look_side`** — `"right"` (Sentinel-1), `"left"`, or `"both"` for a
  satellite that rolls to either side.
- **`incidence_range`** — usable incidence `(near, far)` in degrees, default
  Sentinel-1 IW `(29.1, 46.0)`; `None` for no limit.

## Method

1. **Propagation** — SGP4 (`sgp4` package) in TEME, rotated to Earth-fixed by
   GMST, or ephemerides interpolated (see above); velocity taken relative to
   the ground. UT1−UTC and polar motion ignored: tens of metres today, a few
   hundred at most — far below a TLE's own error (~1 km at its epoch, growing
   1–3 km/day, mostly along track: a time shift of a few seconds, the
   incidence moving by ~0.06° per km across track).
2. **Passes** — the AOI is sampled: a reference point inside, a grid
   (`spacing_km`, 5 km by default), the boundary. A pass is a closest approach
   of the reference point — the **zero-Doppler** time, when a SAR images it —
   with the satellite above its horizon: the zero crossings of
   `(r_sat − r_target) · v_sat` on a 10 s grid, refined by secant steps.
3. **Geometry** — every AOI point is solved for its own zero-Doppler time, then
   - incidence = angle between the ellipsoid normal and the line of sight,
   - side = sign of the target along `v × up`, right of the track or left,
   - off-nadir angle and slant range.

   Incidence grows with ground range, so its extremes over the AOI are on the
   boundary: `incidence_min/max` bound it over the whole area. `coverage` is
   the share of grid points on an allowed side and within the incidence range.
4. **Swath** — for the map, the ground points at near and far incidence in the
   zero-Doppler plane, ray-traced to the ellipsoid, along ±60 s of flight.

## Limits

- Ellipsoid incidence: the local incidence on slopes needs a DEM.
- Opportunities, not a schedule: duty cycle, conflicts between requests and
  the operator's acquisition plan are not modelled. For Sentinel-1 itself, the
  ESA acquisition plans say what *will* be acquired.
- AOI across the antimeridian: handled by taking any AOI wider than 180° in
  longitude as crossing it (`unwrap_antimeridian` moves its negative
  longitudes by +360°, the same points on the ground). An AOI genuinely wider
  than 180° is therefore not supported. The swath is drawn with continuous
  longitudes, beyond 180° if need be.

## Setup

From conda-forge, in the `geo` env: `sgp4`, `numpy`, `pandas`, `shapely`,
`requests`; for the interface `ipyleaflet`, `ipywidgets`.
