# geo

Repo of small, independent geospatial tools. Each subfolder is one self-contained
feature, unrelated to the others — work inside the relevant subfolder and read its
own README before changing anything.

## Environment

- conda/miniforge env named `geo`, Python 3.11, declared in `env_light.yml`.
- Update it with `mamba env update -n geo -f env_light.yml`.
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

## Conventions

- Code, comments and docstrings in English (numpy-style docstrings); conversation
  with the user in French.
- When a tool exists as both a notebook and a module, both hold the same functions:
  any change must be applied to the two in the same edit.

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