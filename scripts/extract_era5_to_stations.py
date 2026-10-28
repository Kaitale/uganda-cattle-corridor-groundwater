"""Extract daily ERA5-Land potential evaporation and precipitation (mm/day) at each located station."""
import glob, os, pandas as pd, xarray as xr
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "data")
st = pd.read_csv(os.path.join(R, "raw", "Stations.csv")).dropna(subset=["latitude", "longitude"])
files = sorted(glob.glob(os.path.join(R, "raw", "era5", "era5land_pev_t2m_*.nc")))
ds = xr.open_mfdataset(files, combine="by_coords")
t = "valid_time" if "valid_time" in ds.dims else "time"
pev, tp = {}, {}
for _, r in st.iterrows():
    pt = ds.sel(latitude=float(r.latitude), longitude=float(r.longitude), method="nearest")
    for var, out in (("pev", pev), ("tp", tp)):
        s = pt[var].to_series(); s.index = pd.to_datetime(s.index) - pd.Timedelta(days=1)  # 00:00 = previous day total
        out[r.station_name] = (s.abs() * 1000) if var == "pev" else (s * 1000).clip(lower=0)
pev = pd.DataFrame(pev)
# ERA5-Land potential evaporation is strongly overestimated over the tropics (about 15-20 mm/day here).
# Keep its daily and seasonal pattern but rescale each station to a long-term mean of 4.5 mm/day
# (about 1,650 mm/yr, typical reference PET for Uganda). Pastas' factor f absorbs any remaining scale error.
pev = pev * (4.5 / pev.mean())
pev.round(2).to_csv(os.path.join(R, "processed", "era5_pev_stations.csv"))
pd.DataFrame(tp).round(2).to_csv(os.path.join(R, "processed", "era5_tp_stations.csv"))
print("ERA5 extracted for", len(pev), "stations")
