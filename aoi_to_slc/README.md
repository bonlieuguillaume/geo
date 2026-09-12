# aoi_to_slc

From an area of interest to a list of Sentinel-1 products to download: a webmap
in a notebook to draw or paste the AOI, list the SLC / GRD scenes covering it,
tick some and write their S3 paths to a text file, one per line. A stripped-down
Copernicus Browser that stops where the download starts: the path file is what
a downloader takes, with an output folder, to fetch the products. Nothing is
downloaded here and no credentials are needed.

- `aoi_to_slc.py` — everything: the search (`search_products`), the path file
  (`write_path_file`, `s3_paths`, `format_s3_path`), the AOI helpers
  (`parse_aoi`, `save_aoi`) and the interface (`build_ui`). Usable as a library
  without the widget packages.
- `aoi_to_slc.ipynb` — a thin driver of the module: parameters, the interface,
  and the same calls from plain Python. It does *not* mirror the module.
- `path_files/` — where the path files land by default (created on demand).

## Backend: the CDSE STAC catalogue, hit directly

Search hits `https://stac.dataspace.copernicus.eu/v1/search` with `requests`
(no `pystac-client`): a POST per page with the AOI as `intersects`, the date
range, a CQL2 filter on `sar:instrument_mode`, `sat:orbit_state` and
`platform`, sorted newest first, following the `next` links up to `max_items`.
Every STAC item exposes each file of the product as an asset with an
`s3://eodata/...` href; the `.SAFE` key is read from the manifest asset and is
what goes in the path file, in one of three styles:

| `style` | line written |
| --- | --- |
| `"mount"` (default) | `/eodata/Sentinel-1/SAR/IW_SLC__1S/2026/09/11/<product>.SAFE` |
| `"s3"` | `s3://eodata/Sentinel-1/SAR/IW_SLC__1S/2026/09/11/<product>.SAFE` |
| `"key"` | `eodata/Sentinel-1/SAR/IW_SLC__1S/2026/09/11/<product>.SAFE` |

**GRD means the COG variant.** CDSE distributes every GRD scene twice: the
original SAFE (`IW_GRDH_1S/.../..._614A.SAFE`) and a Cloud-Optimised one
(`IW_GRDH_1S-COG/.../..._DA9F_COG.SAFE`) — same values and same annotation
XML, measurements as Zstandard-compressed COG (~30 % lighter), always online
where an original older than a year may sit in cold storage. The two are
distinct products with different checksum suffixes, so one cannot be derived
from the other. The STAC catalogue lists the COG only, and this tool goes with
it: SNAP reads it from version 10 with the preprocessing graph unchanged. The
name difference matters only when cross-referencing by name with catalogues
that ignore the COG (ASF, HyP3). SLC products have no such variant.

Should the original GRD ever be needed, the CDSE OData catalogue lists both
variants (`productType` `IW_GRDH_1S` / `IW_GRDH_1S-COG`) with their `S3Path`:
a rewrite of `search_products` behind the same interface, nothing else.

**Why not `asf_search`.** Same Sentinel-1 scenes, but ASF only hands out its
own HTTPS zip URLs — no key in the CDSE `eodata` bucket, which is what the
downloader wants. `asf_search` stays the right tool for bursts as products and
InSAR baselines, which is the `asf/` folder's business.

## Setup

All dependencies come from conda-forge:

```
conda install -c conda-forge ipyleaflet ipywidgets requests geopandas shapely
```

No account, no token: the catalogue search is anonymous.

## Usage

**Notebook.** Fill the parameters cell (path file, path style, map view), run
the interface cell, draw or paste the AOI, search, tick, *Write S3 paths*. The
box under the list previews the lines before they are written. The file is
overwritten each time — one file is one selection — and the AOI is written next
to it as `<name>_aoi.geojson` unless unticked. See the notebook's opening cell.

**Python.**

```python
from aoi_to_slc import search_products, write_path_file

products = search_products(
    "POLYGON ((2.2 48.8, 2.5 48.8, 2.5 49.0, 2.2 49.0, 2.2 48.8))",
    "2026-08-01", "2026-09-12",
    product_type="SLC",          # or "GRD" (COG)
    mode="IW",                   # "IW", "EW", "SM", None
    orbit_direction="descending",
    platforms=["S1C", "S1D"],
)
# GeoDataFrame, newest first: name, datetime, platform, product_type, mode,
# orbit, relative_orbit, absolute_orbit, polarisations, size_gb, s3_key,
# geometry. products.attrs["truncated"] tells if max_items was hit.

chosen = products[products.relative_orbit == 110]
write_path_file(chosen, "path_files/paris_orbit110.txt", style="mount", aoi=aoi)
```

`parse_aoi` accepts a shapely geometry, a GeoJSON dict, inline WKT or GeoJSON,
or the path of a WKT / GeoJSON file. Coordinates are lon/lat (EPSG:4326).

## Limits and things to know

- The CDSE front-end answers HTTP 429 to bursts of requests; the search retries
  with an exponential back-off. Listing a few hundred products takes a little
  longer, nothing more.
- Sentinel-1 SLC bursts as products are not covered: whole products only here.
  Which bursts of a product cover the AOI is `polygon_to_swaths_bursts`'s job.
- The interface needs the Jupyter widgets front-end (ipyleaflet, ipywidgets).
  It works in VS Code; if the map stays blank, JupyterLab is the safe bet.
