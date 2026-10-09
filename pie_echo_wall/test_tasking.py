"""Passes of a SAR satellite over an AOI, from plain Python — fill in the parameters.

Run from this folder: ``python test_tasking.py``. The interactive version is
``tasking.ipynb``.
"""

import pandas as pd

import tasking

# --- Orbit: fill exactly one of the four -------------------------------------
# ELEMENTS is used when the other three are None.

# 1. A TLE, two lines (or three with the name first)
TLE = None
# TLE = """
# 1 39634U 14016A   26277.50000000  .00000123  00000-0  35000-4 0  9991
# 2 39634  98.1820 284.1234 0001300  85.1234 274.9876 14.59198500660000
# """

# 2. A NORAD id: the latest TLE is fetched from CelesTrak
NORAD_ID = None  # e.g. 39634 for Sentinel-1A

# 3. The operator's ephemerides: CCSDS OEM (text), Sentinel-1 .EOF, or a CSV
# table (time, x, y, z[, vx, vy, vz], m or km; its frame in EPHEMERIS_FRAME).
# The search window is clipped to the span of the file.
EPHEMERIS = None  # e.g. r"C:\path\to\predicted_orbit.oem"
EPHEMERIS_FRAME = "ITRF"  # CSV only: "ITRF", "EME2000" or "TEME"

# 4. Mean elements, for a satellite with no TLE yet. Size: altitude_km OR
# repeat_cycle; orientation: ltan_hours OR raan_deg; inclination_deg None =
# sun-synchronous; mean_anomaly_deg = argument of latitude at the epoch for a
# circular orbit (0: at the ascending node).
ELEMENTS = dict(
    epoch="2026-10-05T00:00:00",
    repeat_cycle=(12, 175),
    ltan_hours=18.0,
    mean_anomaly_deg=0.0,
)

# --- Sensor -------------------------------------------------------------------

LOOK_SIDE = "right"                  # "right", "left" or "both" (satellite able to roll)
INCIDENCE_RANGE = (29.1, 46.0)       # usable incidence, near/far, degrees; None = no limit

# --- Area and window ----------------------------------------------------------

AOI = "POLYGON((2.2 48.7, 2.6 48.7, 2.6 49.0, 2.2 49.0, 2.2 48.7))"  # WKT, GeoJSON or a file path
START = "2026-10-05"
END = "2026-10-17"
DIRECTIONS = ("ascending", "descending")


if __name__ == "__main__":
    if EPHEMERIS is not None:
        sat = tasking.read_ephemeris(EPHEMERIS, frame=EPHEMERIS_FRAME)
        print(sat)
    else:
        sat = tasking.satellite(
            tle=TLE, norad_id=NORAD_ID,
            elements=ELEMENTS if TLE is None and NORAD_ID is None else None,
        )
    for key, value in tasking.orbit_summary(sat).items():
        print(f"{key:>18}: {value:.4f}" if isinstance(value, float) else f"{key:>18}: {value}")

    passes = tasking.find_passes(
        sat, AOI, START, END,
        look_side=LOOK_SIDE, incidence_range=INCIDENCE_RANGE, directions=DIRECTIONS,
    )
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(passes.to_string(float_format="%.2f"))
    if EPHEMERIS is not None:
        w0, w1 = passes.attrs["window"]
        print(f"Searched {w0} to {w1} UTC, the part of the window the ephemerides cover.")
    from_tle = EPHEMERIS is None and (TLE is not None or NORAD_ID is not None)
    if from_tle and passes.attrs["epoch_gap_days"] > tasking.TLE_MAX_AGE_DAYS:
        print(f"Window up to {passes.attrs['epoch_gap_days']:.1f} days from the TLE epoch: "
              "expect km-level along-track error.")
