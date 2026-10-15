# parallel_download

Download the Sentinel-1 products of a path file from the Copernicus Data Space
(CDSE) with [rclone](https://rclone.org/), in parallel, into a folder. This is
the step after [`aoi_to_slc`](../aoi_to_slc/README.md), which writes the path
file. No notebook: it is a command.

- `parallel_download.py` — `parallel_download(list_path, dest_dir, ...)`, the
  filter-file builder (`build_filter_lines`, `build_filter_file`), the rclone
  call (`run_rclone_copy`) and the flattening of the downloaded tree
  (`flatten_safe_dirs`), plus the command line.

Adapted from the `download_products` feature of the rosarium repository: same
logic, but the output folder is an argument instead of rosarium's fixed
`data/raw/<folder>/`, and there is no progress hook for a webmap.

## How it works

1. The path file — one product per line in any of the three forms
   `aoi_to_slc` writes (`/eodata/...`, `s3://eodata/...`, `eodata/...`) —
   loses the products whose `.SAFE` is already at the root of the output
   folder (nothing is run when none is left), then the rest is turned into an
   rclone **filter file**, `filter.txt`, next to it: one
   `+ /Sentinel-1/SAR/.../<product>.SAFE/**` line per product and a final
   `- *` that excludes everything else. Lines already in filter syntax are
   kept, so a converted list can be fed again.
2. One `rclone copy cdse:eodata . --filter-from filter.txt` runs in the output
   folder, with `--transfers 8 --multi-thread-streams 8`: every product at
   once, several streams per file. rclone walks only the listed folders, so
   the bucket-wide filter costs nothing.
3. rclone reproduces the bucket's tree (`Sentinel-1/SAR/IW_SLC__1S/2026/09/11/…`);
   every `*.SAFE` folder is moved up to the output folder and the emptied tree
   pruned. A `.SAFE` already present is left alone (warning).

The download resumes: rclone skips files already complete, so an interrupted
run is simply launched again. And a product downloaded once is not downloaded
again: without step 1's check it would be, since rclone looks for it at its
bucket path, which the flattening emptied. A `.SAFE` reaches the root of the
folder only after the whole copy succeeded, so its presence there means
complete.

## Configuring rclone

rclone is not a Python package but a single binary; it comes from conda-forge
(package `rclone`) and is on the PATH as soon as the env is activated.

What does need doing once, and only once, is the **remote**: the binary knows
nothing of CDSE until it is given your S3 keys, which come from
<https://eodata-s3keysmanager.dataspace.copernicus.eu/>. Either run
`rclone config` (interactive) or paste this section into the file that
`rclone config file` points to:

```
[cdse]
type = s3
provider = Other
access_key_id = <your access key>
secret_access_key = <your secret key>
endpoint = https://eodata.dataspace.copernicus.eu
region = default
```

That config file is **per user, not per environment** (`%APPDATA%\rclone\` on
Windows): if the remote was already set up for rosarium, it works here as is.

The remote name (`cdse`) and bucket (`eodata`) are the default of `--remote`;
another name is passed as `--remote name:eodata`. Check with
`rclone lsd cdse:eodata` — it should list `Sentinel-1`, `Sentinel-2`, …

## Usage

From this folder:

```
python parallel_download.py ../aoi_to_slc/path_files/products.txt D:/s1/paris
python parallel_download.py list.txt out --transfers 4 --multi-thread-streams 4
python parallel_download.py --help
```

| Argument | Default | Effect |
| --- | --- | --- |
| `LIST` | — | the path file to download |
| `OUTPUT` | — | folder receiving the `.SAFE` (created if missing) |
| `--transfers N` | `8` | files copied in parallel |
| `--multi-thread-streams N` | `8` | streams per large file |
| `--remote REMOTE:BUCKET` | `cdse:eodata` | the rclone remote |

From Python:

```python
from parallel_download import parallel_download

parallel_download("../aoi_to_slc/path_files/products.txt", "D:/s1/paris")  # -> Path("D:/s1/paris")
```

## Dependencies

Standard library only on the Python side. rclone comes from conda-forge; its
CDSE remote is configured once (above).
