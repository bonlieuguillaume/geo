"""Serve the AOI → Sentinel-1 webmap in the browser: ``python serve.py --open``.

Standard library only. The server hands out viewer.html and answers three JSON
routes with the module's functions — parse the AOI, search the CDSE STAC
catalogue, write the path file. Everything else — the map, the drawing, the
list, the selection — happens in the page, with Leaflet and Leaflet.draw.

Usage:
    python serve.py                       # http://localhost:8050
    python serve.py --open                # and open the browser
    python serve.py --port 9000 --path-file C:/data/list.txt

The default port is 8050 on purpose — 8000 is taken by the Moonfleet viewer,
and both may run at the same time.

The page's defaults — path file, path style, map view, dates — come from the
command line; the path file stays editable in the page.
"""

import argparse
import json
import math
import sys
import threading
import webbrowser
from datetime import date, timedelta
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pandas as pd
from shapely.geometry import mapping

# The module sits next to this script, whatever the current directory
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from aoi_to_slc import (  # noqa: E402
    DEFAULT_MAX_ITEMS,
    DEFAULT_STYLE,
    S3_PATH_STYLES,
    format_s3_path,
    parse_aoi,
    search_products,
    write_path_file,
)

VIEWER = HERE / "viewer.html"

# Filled by main() from the command line, handed to the page by /api/config
CONFIG = {}


# --- JSON routes --------------------------------------------------------------

def _clean(value):
    """A cell of the results table as something json.dumps accepts."""
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if isinstance(value, float) and math.isnan(value):
        return None
    if hasattr(value, "item"):  # numpy scalar
        return value.item()
    return value


def _feature(row):
    """One row of `search_products` as a GeoJSON feature, path included."""
    props = {k: _clean(v) for k, v in row._asdict().items() if k != "geometry"}
    props["path"] = format_s3_path(props["s3_key"], CONFIG["style"]) if props["s3_key"] else None
    return {"type": "Feature", "properties": props, "geometry": mapping(row.geometry)}


def api_config(_body):
    return CONFIG


def api_parse_aoi(body):
    geom = parse_aoi(body.get("text", ""))
    return {"geometry": mapping(geom), "type": geom.geom_type}


def api_search(body):
    for key in ("aoi", "start", "end"):
        if not body.get(key):
            raise ValueError(f"{key} is missing")
    products = search_products(
        body["aoi"], body["start"], body["end"],
        product_type=body.get("product_type") or "SLC",
        mode=body.get("mode") or None,
        orbit_direction=body.get("orbit_direction") or None,
        platforms=body.get("platforms") or None,
        max_items=CONFIG["max_items"],
    )
    print(f"  search: {len(products)} product(s)"
          + (" (truncated)" if products.attrs["truncated"] else ""), flush=True)
    return {
        "type": "FeatureCollection",
        "features": [_feature(row) for row in products.itertuples(index=False)],
        "truncated": products.attrs["truncated"],
    }


def api_write(body):
    products = body.get("products") or []
    if not products:
        raise ValueError("nothing selected")
    path = Path(body.get("path_file") or CONFIG["path_file"])
    if not path.is_absolute():
        path = HERE / path
    aoi = body.get("aoi")
    written = write_path_file(products, path, style=CONFIG["style"], aoi=aoi)
    aoi_name = f"{written.stem}_aoi.geojson" if aoi else None
    print(f"  write: {len(products)} path(s) -> {written}" + (f" (+ {aoi_name})" if aoi_name else ""),
          flush=True)
    return {"path": str(written), "aoi_path": aoi_name, "count": len(products)}


ROUTES = {
    "GET": {"/api/config": api_config},
    "POST": {"/api/parse_aoi": api_parse_aoi, "/api/search": api_search, "/api/write": api_write},
}


# --- Server -------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    """viewer.html at the root, JSON under /api, nothing else."""

    def log_message(self, *args):
        pass  # the routes print what matters; request logs would drown it

    def do_GET(self):
        if self.path in ("/", "/viewer.html"):
            self._send_file(VIEWER, "text/html; charset=utf-8")
        elif self.path in ROUTES["GET"]:
            self._call(ROUTES["GET"][self.path], {})
        else:
            self.send_error(404)

    def do_POST(self):
        route = ROUTES["POST"].get(self.path)
        if route is None:
            self.send_error(404)
            return
        length = int(self.headers.get("Content-Length") or 0)
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
        except json.JSONDecodeError:
            self._send_json({"error": "invalid JSON body"}, 400)
            return
        self._call(route, body)

    def _call(self, route, body):
        # Any error goes back to the page as text, where the status line shows
        # it; the server never dies on a bad AOI or a catalogue hiccup
        try:
            self._send_json(route(body))
        except Exception as exc:  # noqa: BLE001
            self._send_json({"error": f"{type(exc).__name__}: {exc}"}, 400)

    def _send_json(self, payload, status=200):
        self._send(json.dumps(payload).encode("utf-8"), "application/json", status)

    def _send_file(self, path, content_type):
        self._send(path.read_bytes(), content_type)

    def _send(self, data, content_type, status=200):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)


class Server(ThreadingHTTPServer):
    """One server per port. HTTPServer sets SO_REUSEADDR, which on Windows lets
    a second process bind a port already in use — both then run and requests
    land on either. Without it, the second start fails, as it should."""

    allow_reuse_address = False


def main(argv=None):
    parser = argparse.ArgumentParser(description="AOI → Sentinel-1 webmap, served locally")
    parser.add_argument("--port", type=int, default=8050, help="port to serve on (default: 8050)")
    parser.add_argument("--open", action="store_true", help="open the browser once the server is up")
    parser.add_argument(
        "--path-file", default=str(HERE / "path_files" / "products.txt"),
        help="default path file, editable in the page (default: path_files/products.txt here)",
    )
    parser.add_argument(
        "--style", choices=list(S3_PATH_STYLES), default=DEFAULT_STYLE,
        help="form of each line: mount = /eodata/..., s3 = s3://eodata/..., key = eodata/... "
        f"(default: {DEFAULT_STYLE})",
    )
    parser.add_argument("--center", nargs=2, type=float, default=(46.5, 2.5), metavar=("LAT", "LON"),
                        help="initial map centre (default: 46.5 2.5)")
    parser.add_argument("--zoom", type=int, default=6, help="initial zoom (default: 6)")
    parser.add_argument("--days", type=int, default=30, help="initial date range: the last N days (default: 30)")
    parser.add_argument("--max-items", type=int, default=300,
                        help="stop listing after that many products (default: 300)")
    args = parser.parse_args(argv)

    end = date.today()
    CONFIG.update(
        path_file=args.path_file, style=args.style,
        center=list(args.center), zoom=args.zoom,
        start=(end - timedelta(days=args.days)).isoformat(), end=end.isoformat(),
        max_items=args.max_items,
    )

    url = f"http://localhost:{args.port}/"
    try:
        httpd = Server(("", args.port), Handler)
    except OSError as exc:
        # Another server — this one already running, or another project's —
        # holds the port; the raw error is unreadable on Windows
        print(f"  port {args.port} is already in use ({exc.strerror}). "
              "Pick another one with --port, e.g. --port 8051", file=sys.stderr)
        return 1
    with httpd:
        print(f"  AOI → Sentinel-1 webmap at {url}")
        print(f"  path file: {args.path_file} ({args.style} style)")
        print("  Ctrl+C to stop")
        if args.open:
            # A short delay so the server is listening when the page loads
            threading.Timer(0.5, lambda: webbrowser.open(url)).start()
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n  Server stopped.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
