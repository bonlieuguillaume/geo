# geo

Repo of small, independent geospatial tools. Each subfolder is one self-contained
feature, unrelated to the others — work inside the relevant subfolder and read its
own README before changing anything.

## Environment

- conda/miniforge env named `geo`, Python 3.11. No single env file: each tool
  folder carries its own `env_<tool>.yml` (`polygon_to_swaths_bursts/`, `asf/` —
  the latter covers both tools), all declaring the same env `geo`.
- Update it with `mamba env update -n geo -f <tool>/env_<tool>.yml`. The user
  maintains these files and the env themselves: hand them the conda-forge package
  names to add, do not edit the yml files.
- **Every dependency must be installable from conda-forge.** This is a hard
  requirement for all tools here: no pip-only packages, no heavyweight SAR stacks
  (SNAP, ISCE, GAMMA) — the tools reimplement what they need from product metadata.

## Working rules

- The user runs their own terminal, conda env and tests on real data. Do not install
  packages, modify the env, or run commands against their SAR products unless they
  explicitly ask. Ad-hoc checks belong in the scratchpad directory.
- Windows. The "Miniforge Prompt" is `cmd.exe` (no `PS` in the prompt); VS Code's
  integrated terminal is PowerShell. Shell quoting differs between the two — matters
  when documenting CLI examples.
- **Never overwrite the user's own values.** Parameter cells hold what they chose:
  input paths, AOI, mode, resolution, output directories. `NotebookEdit` rewrites a
  whole cell, so editing one line of such a cell silently restores every other value
  from whenever the cell was last read — and the user's edits since are lost. Re-read
  the cell immediately before writing it, carry the current values across verbatim,
  and if a value cannot be confirmed, ask instead of guessing.

## Conventions

- Code, comments and docstrings in English (numpy-style docstrings); conversation
  with the user in French.
- When a tool exists as both a notebook and a module, both hold the same functions:
  any change must be applied to the two in the same edit.
- A notebook whose opening markdown cell describes what each cell does keeps that
  description in sync: adding, removing or reordering a cell means updating the
  table in the same edit. A stale walkthrough is worse than none.

## Maintaining this file

Keep it current as the repo evolves: add a line under Tools for a new tool, and
record a decision or constraint here once it is durable and not derivable from the
code. Be selective — this file is loaded into context every session, so prefer
rewriting an existing line over appending a new one, and drop what no longer helps.

## Tools

- `polygon_to_swaths_bursts/` — Sentinel-1 SLC: which sub-swaths and bursts a polygon
  intersects, read from the product annotation XML without touching the image data.
  Module + CLI, a mirrored notebook, and a `README.md` detailing the algorithm,
  its accuracy limits and its usage.
- `asf/` — turns Sentinel-1 acquisitions into gamma0 RTC VV/VH GeoTIFFs on a
  common grid, by submitting jobs to **ASF HyP3** and regridding the results
  locally. Two notebooks, because the two HyP3 job types share no options:
  `asf_gamma0.ipynb` for whole scenes (`RTC_GAMMA`, SLC or GRD, resolution and
  radiometry configurable) and `asf_gamma0_burst.ipynb` for bursts
  (`OPERA_RTC_S1`, no options at all, 30 m, one co-pol burst id per job). OTB, ISCE2 and ISCE3 were all ruled out first: OTB has no
  terrain flattening and no Windows conda-forge build, the ISCE family is
  Linux-only and does not calibrate radiometrically. Imports
  `polygon_to_swaths_bursts` through `sys.path` (the one cross-folder dependency
  in this repo): renaming or moving that folder breaks this notebook.
- `aoi_to_slc/` — webmap in a notebook (ipyleaflet + ipywidgets): draw or paste
  an AOI, list the Sentinel-1 SLC/GRD scenes covering it, tick some, write their
  S3 paths (`/eodata/Sentinel-1/SAR/.../<product>.SAFE`, the downloader's
  convention; `s3://` and bare-key forms optional) to a text file in
  `path_files/`, one per line, plus the AOI as `<name>_aoi.geojson`. **No download
  here**: the user's downloader, in another repo, takes that file and an output
  folder. Search = CDSE STAC hit directly with `requests`, anonymous, no
  credentials. Not asf_search: no `eodata` paths there. **GRD = the COG variant**,
  the only one in the CDSE STAC — a distinct product from the original GRD (other
  checksum suffix, `IW_GRDH_1S-COG` folder), accepted deliberately since the
  user's chain runs SNAP 13 (COG readable from SNAP 10). If originals are ever
  needed, CDSE OData lists both with `S3Path`: rewrite `search_products` only.
  The module holds everything and the notebook only drives it — not a mirror.