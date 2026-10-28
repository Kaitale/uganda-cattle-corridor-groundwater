"""ERA5-Land potential evaporation and 2 m temperature for Uganda, 1998-2026.
Uses the same CDS request pattern as KaiEWS\\AM_EWS_HEC_RTS\\scripts\\am_fetch_history.py,
which already works on this computer (needs ~/.cdsapirc, already present for AM_EWS).
Run from this folder:  python get_era5_pev.py
ERA5-Land accumulations at 00:00 hold the total for the previous day."""
import os, time, cdsapi
OUT = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "era5")
os.makedirs(OUT, exist_ok=True)
AREA = [4.3, 29.5, -1.5, 35.1]          # N, W, S, E  (all of Uganda)
c = cdsapi.Client()
for y in range(1998, 2027):
    target = os.path.join(OUT, f"era5land_pev_t2m_{y}.nc")
    if os.path.exists(target) and os.path.getsize(target) > 5000:
        print(y, "cached"); continue
    req = dict(variable=["potential_evaporation", "total_precipitation", "2m_temperature"], area=AREA, year=str(y),
               month=[f"{m:02d}" for m in range(1, 13)], day=[f"{d:02d}" for d in range(1, 32)],
               time=["00:00"], data_format="netcdf", download_format="unarchived")
    for k in range(3):
        try:
            c.retrieve("reanalysis-era5-land", req, target + ".part")
            os.replace(target + ".part", target); print(y, "ok"); break
        except Exception as e:
            print(y, "retry", str(e)[:100]); time.sleep(60)
