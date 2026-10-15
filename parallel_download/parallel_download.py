"""Parallel download of Sentinel-1 .SAFE products from CDSE via rclone.

Given a path file — the one `aoi_to_slc` writes, one `/eodata/...` product
path per line — this builds an rclone filter file next to it, runs a single
parallelised `rclone copy` covering every product at once, then flattens the
resulting date/mission tree so every .SAFE folder ends up directly in the
output folder. Products whose .SAFE is already there are left out of the copy.

    python parallel_download.py ../aoi_to_slc/path_files/products.txt D:/s1/paris

rclone comes from conda-forge, but its remote for the CDSE S3 endpoint is
configured once by the user (`rclone config`, named `cdse` by default — see
the README).

Adapted from the `download_products` feature of the rosarium repository,
without its fixed ``data/raw/<folder>`` layout and its webmap progress hook.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_REMOTE = "cdse:eodata"
DEFAULT_TRANSFERS = 8
DEFAULT_MULTI_THREAD_STREAMS = 8


def _normalize_remote_path(raw):
    """Turn a raw S3-style product path into a path relative to the `eodata` bucket root.

    Accepts the three forms `aoi_to_slc` can write (`/eodata/...`,
    `s3://eodata/...`, `eodata/...`).
    """
    raw = raw.strip()
    if raw.lower().startswith("s3://"):
        raw = raw[5:]
    parts = [p for p in raw.split("/") if p and p != "**"]
    if parts and parts[0].lower() == "eodata":
        parts = parts[1:]
    return "/" + "/".join(parts)


def _raw_path_to_filter_line(raw):
    return f"+ {_normalize_remote_path(raw)}/**"


def product_name(raw):
    """The ``<product>.SAFE`` folder name a path-file line points to, or None.

    Works on the three path forms and on lines already in filter syntax
    (``+ /Sentinel-1/.../X.SAFE/**``).
    """
    for part in raw.strip().lstrip("+- ").split("/"):
        if part.upper().endswith(".SAFE"):
            return part
    return None


def missing_lines(raw_lines, dest_dir):
    """The path-file lines whose product is not yet in ``dest_dir``.

    A .SAFE lands at the root of the folder only once the whole copy has
    succeeded (``flatten_safe_dirs`` runs after rclone), so its presence
    there means complete. Without this, a second run would copy it again:
    rclone looks for it at its bucket path (``Sentinel-1/SAR/...``), which
    the flattening emptied.
    """
    dest_dir = Path(dest_dir)
    return [
        raw for raw in raw_lines
        if raw.strip() and not (
            (name := product_name(raw)) and (dest_dir / name).is_dir()
        )
    ]


def build_filter_lines(raw_lines):
    """Convert path-file lines into rclone filter syntax.

    Lines already in filter syntax (starting with "+" or "-") are kept as-is,
    so the function is safe to re-run on a list that was already converted.
    A trailing "- *" catch-all (excluding everything else) is (re-)appended once.
    """
    lines = []
    seen = set()
    for raw in raw_lines:
        raw = raw.strip()
        if not raw:
            continue
        if raw.startswith("+") or raw.startswith("-"):
            if raw.replace(" ", "") == "-*":
                continue
            line = raw
        else:
            line = _raw_path_to_filter_line(raw)
        if line not in seen:
            seen.add(line)
            lines.append(line)

    if not lines:
        raise ValueError("No product paths found to build the rclone filter from.")

    lines.append("- *")
    return lines


def build_filter_file(list_path, raw_lines=None):
    """Write the filter.txt of list_path next to it.

    ``raw_lines`` replaces the content of list_path when given (the lines
    left once the products already downloaded are dropped).
    """
    list_path = Path(list_path)
    if raw_lines is None:
        raw_lines = list_path.read_text(encoding="utf-8").splitlines()

    filter_lines = build_filter_lines(raw_lines)

    filter_path = list_path.parent / "filter.txt"
    filter_path.write_text("\n".join(filter_lines) + "\n", encoding="utf-8")

    return filter_path


def run_rclone_copy(
    filter_path,
    dest_dir,
    transfers=DEFAULT_TRANSFERS,
    multi_thread_streams=DEFAULT_MULTI_THREAD_STREAMS,
    remote=DEFAULT_REMOTE,
):
    """Run a single parallelised rclone copy of every product matched by filter_path into dest_dir.

    rclone draws its live progress in the console (``--progress``).
    """
    rclone = shutil.which("rclone")
    if rclone is None:
        raise FileNotFoundError(
            "rclone not found on the PATH: it comes from conda-forge (package "
            "`rclone`) — activate the geo env, or add the package to it. The CDSE "
            "remote itself is configured once with `rclone config` (see the README "
            "of parallel_download)."
        )
    cmd = [
        rclone, "copy", remote, ".",
        "--filter-from", str(Path(filter_path).resolve()),
        "--transfers", str(transfers),
        "--multi-thread-streams", str(multi_thread_streams),
        "--progress",
    ]
    subprocess.run(cmd, cwd=str(dest_dir), check=True)


def remove_empty_dirs(root):
    """Remove every directory under root left empty after flattening (bottom-up, non-destructive)."""
    root = Path(root)
    for p in sorted(root.rglob("*"), key=lambda q: len(q.parts), reverse=True):
        if p.is_dir():
            try:
                p.rmdir()
            except OSError:
                pass  # not empty (or in use) — leave it alone


def flatten_safe_dirs(dest_dir):
    """Move every *.SAFE directory found anywhere under dest_dir up to its root, then prune the empty tree."""
    dest_dir = Path(dest_dir)
    safe_dirs = sorted(p for p in dest_dir.rglob("*.SAFE") if p.is_dir())

    moved = []
    for safe_dir in safe_dirs:
        target = dest_dir / safe_dir.name
        if safe_dir == target:
            moved.append(target)
            continue
        if target.exists():
            print(f"[WARN] {target} already exists, skipping {safe_dir}")
            continue
        safe_dir.rename(target)
        moved.append(target)

    remove_empty_dirs(dest_dir)
    return moved


def parallel_download(
    list_path,
    dest_dir,
    transfers=DEFAULT_TRANSFERS,
    multi_thread_streams=DEFAULT_MULTI_THREAD_STREAMS,
    remote=DEFAULT_REMOTE,
):
    """Download every Sentinel-1 product listed in a path file into a folder.

    Products whose .SAFE is already at the root of that folder are skipped
    (see `missing_lines`); when none is left, rclone is not run at all.

    Parameters
    ----------
    list_path : str | Path
        Path file: one product per line, `/eodata/...`, `s3://eodata/...` or
        `eodata/...` (what `aoi_to_slc` writes). The rclone filter file is
        written next to it, as ``filter.txt``.
    dest_dir : str | Path
        Folder receiving the .SAFE products (created if missing, reused if
        present).
    transfers : int
        rclone ``--transfers``: files copied in parallel.
    multi_thread_streams : int
        rclone ``--multi-thread-streams``: streams per large file.
    remote : str
        rclone ``remote:bucket`` to copy from.

    Returns
    -------
    Path
        ``dest_dir``, holding the flattened .SAFE products.
    """
    list_path = Path(list_path)
    if not list_path.is_file():
        raise FileNotFoundError(f"path file not found: {list_path}")
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    raw_lines = [line for line in list_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    todo = missing_lines(raw_lines, dest_dir)
    if len(todo) < len(raw_lines):
        print(f"{len(raw_lines) - len(todo)} product(s) already in {dest_dir}, skipped")
    if not todo:
        print("Nothing left to download")
        return dest_dir

    filter_path = build_filter_file(list_path, todo)
    print(f"Filter file written to {filter_path}")

    run_rclone_copy(filter_path, dest_dir, transfers, multi_thread_streams, remote)

    moved = flatten_safe_dirs(dest_dir)
    print(f"{len(moved)} .SAFE product(s) available in {dest_dir}")

    return dest_dir


# --- Command-line interface ---------------------------------------------------

def _build_parser():
    prog = "python parallel_download.py"
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Parallel download of Sentinel-1 .SAFE products from CDSE via rclone.\n\n"
            "Reads a path file (what aoi_to_slc writes: one /eodata/... product path\n"
            "per line), builds an rclone filter file next to it, runs a single\n"
            "parallelised `rclone copy` for every product at once, then flattens the\n"
            "resulting mission/date tree so every .SAFE ends up directly in OUTPUT.\n"
            "Products already there are skipped.\n\n"
            "rclone comes from conda-forge; what it needs is a remote configured once\n"
            "for the CDSE S3 endpoint, with your own keys from\n"
            "https://eodata-s3keysmanager.dataspace.copernicus.eu/ (`rclone config`;\n"
            "the default remote name is `cdse`, bucket `eodata`)."
        ),
        epilog=(
            "examples:\n"
            f"  {prog} ../aoi_to_slc/path_files/products.txt D:/s1/paris\n"
            f"  {prog} list.txt out --transfers 4 --multi-thread-streams 4   # gentler on the connection"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("list", metavar="LIST", help="path file to download")
    parser.add_argument("output", metavar="OUTPUT", help="folder receiving the .SAFE products")
    parser.add_argument("--transfers", type=int, default=DEFAULT_TRANSFERS, metavar="N",
                        help=f"[optional] rclone --transfers value (default: {DEFAULT_TRANSFERS})")
    parser.add_argument("--multi-thread-streams", type=int, default=DEFAULT_MULTI_THREAD_STREAMS, metavar="N",
                        help=f"[optional] rclone --multi-thread-streams value (default: {DEFAULT_MULTI_THREAD_STREAMS})")
    parser.add_argument("--remote", default=DEFAULT_REMOTE, metavar="REMOTE:BUCKET",
                        help=f"[optional] rclone remote to copy from (default: {DEFAULT_REMOTE!r})")
    return parser


def main(argv=None):
    args = _build_parser().parse_args(argv)
    parallel_download(
        args.list,
        args.output,
        transfers=args.transfers,
        multi_thread_streams=args.multi_thread_streams,
        remote=args.remote,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
