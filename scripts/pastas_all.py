"""Pastas recharge models (rain + ERA5 evaporation) for every located station with >= 3 years of levels."""
import os, re, pandas as pd, numpy as np, pastas as ps, matplotlib, warnings
matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); ps.set_log_level("ERROR")
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "data"); FIG = os.path.join(H, "..", "figures", "pastas")
os.makedirs(FIG, exist_ok=True)
st = pd.read_csv(os.path.join(R, "raw", "Stations.csv")).dropna(subset=["latitude", "longitude"])
pev = pd.read_csv(os.path.join(R, "processed", "era5_pev_stations.csv"), index_col=0, parse_dates=True)
tp = pd.read_csv(os.path.join(R, "processed", "era5_tp_stations.csv"), index_col=0, parse_dates=True)
rows = []
for _, s in st.iterrows():
    f = os.path.join(R, "raw", "levels", s.data_file)
    try:
        d = pd.read_csv(f, encoding="latin1", on_bad_lines="skip")
        dc = d.columns[0]; d[dc] = pd.to_datetime(d[dc], errors="coerce"); d = d.dropna(subset=[dc]).set_index(dc).sort_index()
        d = d[~d.index.duplicated()]
        h = -pd.to_numeric(d.iloc[:, 0], errors="coerce").dropna()
        h = h[(h - h.median()).abs() < 5 * h.std()]
        if len(h) < 365 * 3 or (h.index[-1] - h.index[0]).days < 365 * 3:
            continue
        name = s.station_name
        p = tp[name].dropna(); e = pev[name].dropna()
        rain_src = "ERA5-Land"
        if d.shape[1] > 1 and re.search("rf|rain", str(d.columns[1]), re.I):
            rf = pd.to_numeric(d.iloc[:, 1], errors="coerce").clip(lower=0)
            if rf.notna().sum() > 0.5 * len(h):          # station gauge good enough: use it, gap-fill with ERA5
                p = rf.asfreq("D").combine_first(p).dropna(); rain_src = "station + ERA5 fill"
        h = h.loc[max(h.index[0], p.index[0], e.index[0]):min(h.index[-1], p.index[-1], e.index[-1])]
        ml = ps.Model(h, name=name)
        rm = ps.RechargeModel(p / 1000, e / 1000, ps.Gamma(), name="recharge", recharge=ps.rch.Linear())
        ml.add_stressmodel(rm)
        ml.solve(report=False)
        par = ml.get_parameters("recharge")
        rows.append(dict(station=name, wmz=s.wmz, purpose=s.purpose, start=h.index[0].date(), end=h.index[-1].date(),
                         n=len(h), rain=rain_src, EVP=round(ml.stats.evp(), 1), RMSE_m=round(ml.stats.rmse(), 2),
                         t50_days=round(rm.rfunc.get_tmax(par, cutoff=0.5)), t95_days=round(rm.rfunc.get_tmax(par, cutoff=0.95)),
                         evap_factor_f=round(float(ml.parameters.loc["recharge_f", "optimal"]), 2)))
        ml.plots.results(figsize=(10, 7)); plt.savefig(os.path.join(FIG, f"{name}.png".replace("/", "-")), dpi=100, bbox_inches="tight"); plt.close("all")
        print("ok", name, rows[-1]["EVP"])
    except Exception as ex:
        print("fail", s.station_name, str(ex)[:120])
out = pd.DataFrame(rows).sort_values("EVP", ascending=False)
out.to_csv(os.path.join(R, "processed", "pastas_recharge_summary.csv"), index=False)
print(out.to_string())
