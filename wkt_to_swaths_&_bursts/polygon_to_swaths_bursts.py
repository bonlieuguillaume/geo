"""Find the Sentinel-1 SLC swaths and bursts intersecting a polygon.

The area of interest is given as WKT or GeoJSON, either as a literal string or
as a path to a file; the format is detected automatically.

Usable both as a library (import the functions below) and as a command-line
tool (see `python polygon_to_swaths_bursts.py --help`).
"""

import argparse
import fnmatch
import json
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely import wkt as shapely_wkt
from shapely.affinity import translate
from shapely.errors import ShapelyError
from shapely.geometry import MultiPolygon, Polygon, shape
from shapely.ops import unary_union

DEFAULT_COARSE_MARGIN = 0.02  # degrees, ~2 km


def _geometry_from_geojson(obj):
    """Build a single shapely geometry from a parsed GeoJSON object.

    Accepts a bare geometry, a Feature, a FeatureCollection or a
    GeometryCollection; several features are merged into one geometry.
    """
    kind = obj.get("type")
    if kind == "FeatureCollection":
        geoms = [
            shape(f["geometry"]) for f in obj.get("features", []) if f.get("geometry")
        ]
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


def parse_polygon(polygon):
    """Parse an area of interest given as WKT or GeoJSON, literal or file path.

    The format is detected automatically:
      - an existing file path: its content is read, then parsed as below;
      - a dict, or a string starting with '{': GeoJSON;
      - anything else: WKT.

    Coordinates are expected in lon/lat (EPSG:4326), the convention of both
    Sentinel-1 annotations and GeoJSON (RFC 7946).
    """
    if isinstance(polygon, (str, Path)):
        text = str(polygon).strip()
        # A path has no newline and stays short: cheap enough to probe
        if "\n" not in text and len(text) < 4096:
            try:
                if Path(text).is_file():
                    text = Path(text).read_text(encoding="utf-8").strip()
            except (OSError, ValueError):
                pass  # not a usable path, treat the string as a geometry
        polygon = text

    if isinstance(polygon, dict):
        return _geometry_from_geojson(polygon)
    if polygon.startswith("{"):
        return _geometry_from_geojson(json.loads(polygon))
    return shapely_wkt.loads(polygon)


def _read_annotation_files(slc_path):
    """Read the annotation XML files of an SLC product (.SAFE or .zip).

    Returns a list of (file_name, file_bytes). Only product annotations are
    kept (annotation/s1*.xml), not calibration/ nor rfi/.
    """
    slc_path = Path(slc_path)
    if not slc_path.exists():
        raise FileNotFoundError(slc_path)

    pattern = "*/annotation/s1*.xml"
    if slc_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(slc_path) as zf:
            names = [
                n for n in zf.namelist()
                if fnmatch.fnmatch(n, pattern)
                and "/calibration/" not in n
                and "/rfi/" not in n
            ]
            return [(Path(n).name, zf.read(n)) for n in sorted(names)]

    files = sorted((slc_path / "annotation").glob("s1*.xml"))
    if not files:
        raise ValueError(f"No annotation file found in {slc_path}")
    return [(f.name, f.read_bytes()) for f in files]


def _unwrap_antimeridian(coords):
    """Fix a ring straddling the antimeridian (±180°).

    If the longitudes span more than 180°, the ring must cross the ±180 line
    (a genuine ring cannot be that wide): shift negative longitudes by +360 to
    get a continuous polygon in a [0, 360) frame.
    """
    lons = [lon for lon, _ in coords]
    if max(lons) - min(lons) > 180:
        coords = [(lon + 360 if lon < 0 else lon, lat) for lon, lat in coords]
    return coords


def _unwrap_polygon(poly):
    """Apply the antimeridian fix to every ring of a polygon."""
    return Polygon(
        _unwrap_antimeridian(list(poly.exterior.coords)),
        [_unwrap_antimeridian(list(ring.coords)) for ring in poly.interiors],
    )


def _unwrap_geometry(geom):
    """Apply the antimeridian fix to a (multi)polygon; other types pass through."""
    if geom.geom_type == "Polygon":
        return _unwrap_polygon(geom)
    if geom.geom_type == "MultiPolygon":
        return MultiPolygon([_unwrap_polygon(p) for p in geom.geoms])
    return geom


def _burst_footprints_from_annotation(xml_bytes):
    """Extract the burst footprints from one annotation file.

    Returns a list of dicts {swath, polarisation, burst, geometry}.
    """
    root = ET.fromstring(xml_bytes)
    swath = root.findtext(".//adsHeader/swath")
    polarisation = root.findtext(".//adsHeader/polarisation")
    lines_per_burst = int(root.findtext(".//swathTiming/linesPerBurst"))
    n_bursts = len(root.findall(".//swathTiming/burstList/burst"))

    pts = np.array([
        (
            int(p.findtext("line")),
            int(p.findtext("pixel")),
            float(p.findtext("latitude")),
            float(p.findtext("longitude")),
        )
        for p in root.findall(".//geolocationGridPoint")
    ])
    grid_lines = np.unique(pts[:, 0])

    def ring_points(target_line, reverse=False):
        # Grid points on the row closest to target_line
        line = grid_lines[np.argmin(np.abs(grid_lines - target_line))]
        sel = pts[pts[:, 0] == line]
        sel = sel[np.argsort(sel[:, 1])]  # sort by pixel (range)
        if reverse:
            sel = sel[::-1]
        return [(lon, lat) for _, _, lat, lon in sel]

    bursts = []
    for k in range(n_bursts):
        first_line = k * lines_per_burst
        last_line = (k + 1) * lines_per_burst
        coords = ring_points(first_line) + ring_points(last_line, reverse=True)
        coords = _unwrap_antimeridian(coords)
        bursts.append({
            "swath": swath,
            "polarisation": polarisation,
            "burst": k + 1,  # 1-based numbering within the product
            "geometry": Polygon(coords),
        })
    return bursts


def load_burst_footprints(slc_path):
    """Footprints of every burst in an SLC product, as a GeoDataFrame (EPSG:4326).

    Footprints are identical across polarisations (VV/VH), so only one
    polarisation per swath is kept.
    """
    records = []
    seen_swaths = set()
    for _, xml_bytes in _read_annotation_files(slc_path):
        bursts = _burst_footprints_from_annotation(xml_bytes)
        if bursts and bursts[0]["swath"] not in seen_swaths:
            seen_swaths.add(bursts[0]["swath"])
            records.extend(bursts)
    return gpd.GeoDataFrame(records, crs="EPSG:4326")


def get_intersecting_bursts(
    slc_path, polygon, coarse=False, coarse_margin=DEFAULT_COARSE_MARGIN
):
    """Main entry point: swaths and bursts intersecting the polygon.

    Parameters
    ----------
    slc_path : str | Path
        Path to the SLC product (.SAFE directory or .zip archive).
    polygon : str | Path | dict
        Area of interest as WKT or GeoJSON, given as a literal string, a
        parsed GeoJSON dict, or a path to a file holding either. The format
        is detected automatically; coordinates are lon/lat (EPSG:4326).
    coarse : bool
        False (default): strict test against the edge-matched footprints.
        True: footprints are dilated by coarse_margin before the test —
        favours recall: an AOI close to a burst seam or a swath edge also
        returns the neighbouring bursts, whose valid data extend beyond the
        edge-matched footprints.
    coarse_margin : float
        Dilation margin in degrees (0.02 ~ 2 km). Ignored when coarse=False.

    Returns
    -------
    hits : GeoDataFrame
        The intersecting bursts (columns swath, burst, geometry).
    summary : dict
        {swath: [burst numbers]}, e.g. {"IW1": [3, 4], "IW2": [3, 4, 5]}.
    """
    aoi = parse_polygon(polygon)
    # Same antimeridian fix as for the footprints
    aoi = _unwrap_geometry(aoi)
    # Fixed footprints may live in a [0, 360) frame: test the polygon against
    # all three possible frames.
    test_geom = unary_union(
        [aoi, translate(aoi, xoff=360), translate(aoi, xoff=-360)]
    )

    footprints = load_burst_footprints(slc_path)
    test_footprints = (
        footprints.buffer(coarse_margin) if coarse else footprints.geometry
    )
    # The result keeps the original, undilated geometries
    hits = footprints[test_footprints.intersects(test_geom)].copy()
    summary = (
        hits.groupby("swath")["burst"].apply(lambda s: sorted(s.tolist())).to_dict()
    )
    return hits, summary


DESCRIPTION = """\
Find the Sentinel-1 SLC swaths and bursts intersecting a polygon.

The area of interest is read as WKT or GeoJSON, passed either inline or as a
file path; the format is detected automatically. Coordinates must be lon/lat
(EPSG:4326).

How it works: no SAR library is needed, everything comes from the product
annotation XML files (annotation/s1*.xml inside the .SAFE). Each file holds
swathTiming/linesPerBurst, which splits the image into fixed-size azimuth
tiles, and a geolocation grid whose rows sit on the burst boundaries. The
footprint of a burst is rebuilt from the grid points of its first and last
line, then tested for intersection against the polygon with shapely.

Bursts are numbered from 1 within each swath, as in SNAP's TOPSAR-Split.
Polygons straddling the antimeridian (+/-180 deg) are handled.
"""

EPILOG = """\
examples:
  # inline WKT, strict test, human-readable output
  python polygon_to_swaths_bursts.py S1B_IW_SLC__1SDV_....SAFE "POLYGON ((2.2 48.8, 2.5 48.8, 2.5 49.0, 2.2 49.0, 2.2 48.8))"

  # AOI read from a GeoJSON file, recall-oriented test, JSON output
  python polygon_to_swaths_bursts.py product.zip aoi.geojson --coarse --json

  # write the footprints of the selected bursts as GeoJSON
  python polygon_to_swaths_bursts.py product.zip aoi.wkt --geojson hits.geojson
"""


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="polygon_to_swaths_bursts.py",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "slc_path",
        help="path to the Sentinel-1 SLC product (.SAFE directory or .zip archive)",
    )
    parser.add_argument(
        "polygon",
        help="area of interest: a WKT or GeoJSON string, or a path to a file "
        "holding either (format detected automatically), in lon/lat EPSG:4326",
    )
    parser.add_argument(
        "--coarse",
        action="store_true",
        help="favour recall: dilate the footprints before the intersection test, "
        "so an AOI close to a burst seam or a swath edge also returns the "
        "neighbouring bursts (their valid data overlap by ~1-2 km)",
    )
    parser.add_argument(
        "--coarse-margin",
        type=float,
        default=DEFAULT_COARSE_MARGIN,
        metavar="DEG",
        help="dilation margin in degrees used by --coarse "
        f"(default: {DEFAULT_COARSE_MARGIN}, i.e. ~2 km)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="print the result as JSON ({swath: [burst numbers]}) instead of text",
    )
    parser.add_argument(
        "--geojson",
        metavar="PATH",
        help="also write the footprints of the intersecting bursts to this "
        "GeoJSON file",
    )
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        hits, summary = get_intersecting_bursts(
            args.slc_path,
            args.polygon,
            coarse=args.coarse,
            coarse_margin=args.coarse_margin,
        )
    except (FileNotFoundError, ValueError, KeyError, ShapelyError) as exc:
        parser.exit(2, f"error: {exc}\n")

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    elif summary:
        print("Intersecting swaths:", ", ".join(sorted(summary)))
        for swath, bursts in sorted(summary.items()):
            print(f"  {swath}: bursts {', '.join(str(b) for b in bursts)}")
    else:
        print("No intersecting burst found.")

    if args.geojson:
        hits.to_file(args.geojson, driver="GeoJSON")
        print(f"Footprints written to {args.geojson}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
