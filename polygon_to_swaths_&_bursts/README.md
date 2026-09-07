# polygon → swaths & bursts (Sentinel-1 SLC)

Given a Sentinel-1 SLC product and an area of interest, return the sub-swaths and
the bursts that the AOI intersects.

---

## 1. The problem

A Sentinel-1 IW SLC product covers roughly 250 × 170 km and weighs 4–8 GB. An area
of interest — a mining site, a city, a landslide — usually covers a few square
kilometres of it. Processing the whole product to study that patch wastes hours of
computation and tens of gigabytes of disk.

Every SAR toolbox therefore offers a burst selection step (`TOPSAR-Split` in SNAP,
its equivalents in ISCE or GAMMA), but they all expect the answer to be already
known: *which* sub-swath, and *which* burst numbers inside it. Finding them by hand
means opening the product in a GUI, eyeballing the burst outlines against the AOI,
and writing down indices — slow, error-prone, and impossible to script over a stack
of a hundred dates.

This tool answers that question programmatically:

> input: a path to an SLC product + a polygon (WKT or GeoJSON)
> output: `{"IW1": [1, 2], "IW2": [2, 3]}`

It relies only on packages available from conda-forge (numpy, shapely, geopandas) —
no SNAP, no GDAL-heavy SAR stack, and it never opens the multi-gigabyte image data.

---

## 2. How Sentinel-1 IW products are organised

**Sub-swaths.** In its main acquisition mode (IW, Interferometric Wide swath), the
Sentinel-1 radar does not stare at a single strip. It cycles its antenna across
three adjacent sub-swaths — **IW1** (closest to nadir), **IW2**, **IW3** — which
together span the ~250 km swath. The wide-area mode EW works the same way with five
sub-swaths, EW1 to EW5. Adjacent sub-swaths overlap by roughly 1–2 km so that the
coverage has no gap.

**Bursts.** Within a sub-swath, the antenna is also steered *along* the track
(the TOPS technique), which chops the acquisition into a series of short takes
called **bursts**. Each burst images about 20 km along track, and a standard product
slice holds around nine of them per sub-swath. Consecutive bursts overlap by roughly
a kilometre — a deliberate margin that lets the deburst step stitch them without a
seam, and that interferometric processing exploits for azimuth calibration (ESD).

**On disk.** A `.SAFE` product is a directory (often distributed zipped):

```
S1B_IW_SLC__1SDV_20170804T215105_..._B333.SAFE/
├── measurement/            one GeoTIFF per sub-swath × polarisation (the pixels)
├── annotation/             one XML per sub-swath × polarisation (the metadata)
│   ├── calibration/
│   └── rfi/                radio-frequency-interference reports
└── manifest.safe
```

Inside one measurement file, the bursts are stacked vertically as **tiles of
identical size**: `linesPerBurst` azimuth lines each. The valid data does not fill
its tile entirely — the edges, where the synthetic aperture is incomplete, are set
to zero (*black-fill*), and the annotation records the real extent line by line in
`firstValidSample` / `lastValidSample`.

---

## 3. The method

Everything needed is in the annotation XML; the image data is never read.

**Step 1 — Collect the annotation files.** `annotation/s1*.xml`, read directly from
inside the `.zip` when the product is still archived (no extraction). The
`calibration/` and `rfi/` sub-directories are skipped: they share the same header
structure but carry no geometry.

**Step 2 — Deduplicate polarisations.** There is one annotation file per sub-swath
*and* per polarisation (6 files for a dual-pol IW product). VV and VH image exactly
the same ground, so only the first polarisation met for each sub-swath is kept.

**Step 3 — Cut the image into bursts.** The XML gives `swathTiming/linesPerBurst`
and the list of bursts. Burst *k* (0-based internally) occupies azimuth lines
`[k · linesPerBurst, (k+1) · linesPerBurst]`.

**Step 4 — Turn line numbers into ground coordinates.** The `geolocationGrid`
element is a list of tie points, each mapping an image position to the ground:
`(line, pixel) → (latitude, longitude)`. The decisive property is that **the rows of
this grid are placed on the burst boundaries** — so no interpolation is required:

```
line=0     •  •  •  •  •  …  •      ← top edge of burst 1
line=1500  •  •  •  •  •  …  •      ← bottom of burst 1 / top of burst 2
line=3000  •  •  •  •  •  …  •      ← bottom of burst 2 / top of burst 3
line=4499  •  •  •  •  •  …  •      ← last line of the image
```

**Step 5 — Rebuild each burst footprint.** For burst *k*, take the grid row at its
first line from left to right, then the grid row at its last line from right to
left, and close the ring:

```
first line, increasing pixel  →   A ─ B ─ C ─ D ─ E
                                  │               │
last line, decreasing pixel   ←   J ─ I ─ H ─ G ─ F
```

Reversing the second row is what keeps the ring from crossing itself. Keeping every
tie point rather than just the four corners matters too: a burst on the ground is
not a rectangle, and the ~21 points per edge let the polygon follow its real
curvature.

The row is selected as the one *closest* to the theoretical boundary rather than
strictly equal to it — a safeguard, since the last row of the grid sits at
`numberOfLines − 1`, not at `numberOfLines`.

**Step 6 — Parse the area of interest.** WKT or GeoJSON, passed inline or as a file
path; the format is detected automatically. A GeoJSON `FeatureCollection` with
several features is merged into a single geometry. Coordinates are expected in
lon/lat (EPSG:4326), the convention of both the Sentinel-1 annotations and
GeoJSON (RFC 7946).

**Step 7 — Handle the antimeridian.** Longitude jumps from +180 to −180 in the
middle of the Pacific. A footprint straddling that line would, taken literally,
produce a polygon wrapping the wrong way around the globe. When a ring's longitudes
span more than 180° — impossible for a genuine burst — negative longitudes are
shifted by +360 to make it continuous again. The AOI gets the same treatment, and
the intersection is tested against the AOI plus its ±360° copies so that both frames
match whichever way each geometry was unwrapped.

**Step 8 — Intersect and summarise.** Footprints go into a GeoDataFrame
(EPSG:4326), shapely tests them against the AOI, and the hits are grouped into
`{swath: [burst numbers]}`. Burst numbers are 1-based within each sub-swath, the
same convention as SNAP's `TOPSAR-Split`.

### Strict mode vs. coarse mode

The footprints built at step 5 are **edge-matched**: two consecutive bursts share
their boundary row, so they touch without overlapping. The real valid data does
overlap by about a kilometre, and asymmetrically — the upstream burst extends past
the seam, while the downstream burst starts with a few hundred metres of black-fill
*after* it.

Consequently, an AOI whose edge falls near a seam may need a burst that the strict
test does not return. That is what `coarse=True` (`--coarse` on the command line) is
for: footprints are dilated by `coarse_margin` (0.02° ≈ 2 km) **for the test only**,
so neighbouring bursts are also returned. The geometries in the result stay
undilated. Use it whenever the AOI is not comfortably inside a single burst; the
cost of splitting one extra burst is small, and the deburst step handles the
duplicated strip cleanly.

### Known limits

- Footprints derive from the geolocation grid, not from `firstValidSample` /
  `lastValidSample`; they are accurate to roughly a kilometre near burst edges.
- Burst numbers are local to the product — they are not the ESA global burst IDs.
- The AOI must be in lon/lat WGS84; a polygon in a projected CRS (Lambert-93, UTM)
  must be reprojected first, and `POLYGON((lon lat, …))` is longitude-first.
- An AOI legitimately wider than 180° of longitude would be mistaken for an
  antimeridian crossing.

---

## 4. Usage

### Install

```bash
conda env update -n geo -f ../env_light.yml     # numpy, shapely, geopandas, folium
```

### Command line

```bash
# inline WKT, strict test
python polygon_to_swaths_bursts.py product.SAFE "POLYGON ((2.2 48.8, 2.5 48.8, 2.5 49.0, 2.2 49.0, 2.2 48.8))"

# AOI read from a GeoJSON file, recall-oriented, JSON output
python polygon_to_swaths_bursts.py product.zip aoi.geojson --coarse --json

# also export the footprints of the selected bursts
python polygon_to_swaths_bursts.py product.zip aoi.wkt --geojson hits.geojson
```

```
Intersecting swaths: IW1, IW2
  IW1: bursts 1, 2
  IW2: bursts 2, 3
```

| Option | Effect |
| --- | --- |
| `--coarse` | dilate footprints before the test (favours recall) |
| `--coarse-margin DEG` | dilation margin in degrees, default `0.02` (~2 km) |
| `--json` | print `{swath: [bursts]}` as JSON instead of text |
| `--geojson PATH` | write the footprints of the selected bursts to a file |

`python polygon_to_swaths_bursts.py --help` prints the same description as this
section.

### Python

```python
from polygon_to_swaths_bursts import get_intersecting_bursts, load_burst_footprints

hits, summary = get_intersecting_bursts(
    "product.SAFE",
    "aoi.geojson",        # or inline WKT / GeoJSON, or a dict
    coarse=True,          # optional, default False
    coarse_margin=0.02,   # optional
)
# summary -> {"IW1": [1, 2], "IW2": [2, 3]}
# hits    -> GeoDataFrame (swath, polarisation, burst, geometry)

footprints = load_burst_footprints("product.SAFE")   # every burst, for inspection
```

### Notebook

`polygon_to_swaths_&_bursts.ipynb` holds the same functions cell by cell, plus a
folium map that draws every burst of the product — the selected ones in red, the
others in grey, the AOI in blue. It is the quickest way to check a result visually.
If the map does not render in VS Code, save it and open it in a browser:
`m.save("map.html")`.

### Files

| File | Role |
| --- | --- |
| `polygon_to_swaths_bursts.py` | library + command-line interface |
| `polygon_to_swaths_&_bursts.ipynb` | same functions, interactive, with a map |
