# geo
Repository to develop and test separate and unrelated features.

- `polygon_to_swaths_bursts/` — which Sentinel-1 SLC swaths and bursts a polygon intersects.
- `asf/` — gamma0 RTC from Sentinel-1 scenes or bursts, through ASF HyP3.
- `aoi_to_slc/` — webmap to draw an AOI, list the Sentinel-1 scenes covering it and write their CDSE paths to a file.
- `parallel_download/` — download the products of such a path file from CDSE, in parallel, with rclone.
- `pie_echo_wall/` — SAR tasking: passes of a satellite over an AOI, with their incidence, from its orbit.
