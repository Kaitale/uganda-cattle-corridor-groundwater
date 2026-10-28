"""Pastas recharge + linear trend models, SGI drought index and monthly drought thresholds per station."""
import os, re, pandas as pd, numpy as np, pastas as ps, matplotlib, warnings
matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); ps.set_log_level("ERROR")
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "data")
FT = os.path.join(H, "..", "figures", "pastas_trend"); FD = os.path.join(H, "..", "figures", "drought")
os.makedirs(FT, exist_ok=True); os.makedirs(FD, exist_ok=True)
st = pd.read_csv(os.path.join(R, "raw", "Stations.csv")).dropna(subset=["latitude", "longitude"])
pev = pd.read_csv(os.path.join(R, "processed", "era5_pev_stations.csv"), index_col=0, parse_dates=True)
tp = pd.read_csv(os.path.join(R, "processed", "era5_tp_stations.csv"), index_col=0, parse_dates=True)
rows, thr_rows = [], []
for _, s in st.iterrows():
    name = s.station_name
    try:
        d = pd.read_csv(os.path.join(R, "raw", "levels", s.data_file), encoding="latin1", on_bad_lines="skip")
        dc = d.columns[0]; d[dc] = pd.to_datetime(d[dc], errors="coerce"); d = d.dropna(subset=[dc]).set_index(dc).sort_index()
        d = d[~d.index.duplicated()]
        dep = pd.to_numeric(d.iloc[:, 0], errors="coerce").dropna()           # depth below ground, m
        dep = dep[(dep - dep.median()).abs() < 5 * dep.std()]
        if len(dep) < 365 * 3 or (dep.index[-1] - dep.index[0]).days < 365 * 3:
            continue
        h = -dep
        p = tp[name].dropna(); e = pev[name].dropna(); rain_src = "ERA5-Land"
        if d.shape[1] > 1 and re.search("rf|rain", str(d.columns[1]), re.I):
            rf = pd.to_numeric(d.iloc[:, 1], errors="coerce").clip(lower=0)
            if rf.notna().sum() > 0.5 * len(h):
                p = rf.asfreq("D").combine_first(p).dropna(); rain_src = "station + ERA5 fill"
        h = h.loc[max(h.index[0], p.index[0], e.index[0]):min(h.index[-1], p.index[-1], e.index[-1])]
        ml = ps.Model(h, name=name)
        rm = ps.RechargeModel(p / 1000, e / 1000, ps.Gamma(), name="recharge", recharge=ps.rch.Linear())
        ml.add_stressmodel(rm)
        ml.add_stressmodel(ps.LinearTrend(start=h.index[0], end=h.index[-1], name="trend"))
        ml.solve(report=False)
        par = ml.get_parameters("recharge")
        slope = float(ml.parameters.loc["trend_a", "optimal"]) * 365.25       # m/yr (head; + = rising)
        # --- drought indicators from observed monthly depth ---
        mon = dep.resample("MS").mean().dropna()
        sgi = ps.stats.sgi(-mon)                                              # standardised groundwater index (head)
        thr = mon.groupby(mon.index.month).quantile([0.5, 0.8, 0.9]).unstack()
        thr.columns = ["normal_P50_mbgl", "moderate_P80_mbgl", "severe_P90_mbgl"]
        for m, r in thr.iterrows():
            thr_rows.append(dict(station=name, month=m, **r.round(2).to_dict()))
        yrs = (mon.index - mon.index[0]).days / 365.25
        from scipy.stats import theilslopes, kendalltau
        obs_slope = -theilslopes(mon.values, yrs)[0]                       # head change, m/yr (+ = rising)
        mk_p = kendalltau(yrs, mon.values)[1]
        run = (sgi < -1).astype(int); longest = int(run.groupby((run != run.shift()).cumsum()).sum().max())
        last = sgi.index[-1]; sgi_last = float(sgi.iloc[-1])
        drought_months = int((sgi < -1).sum()); severe = int((sgi < -1.5).sum())
        status = "severe drought" if sgi_last < -1.5 else "drought" if sgi_last < -1 else "below normal" if sgi_last < -0.5 else "normal or above"
        rows.append(dict(station=name, wmz=s.wmz, purpose=s.purpose, start=h.index[0].date(), end=h.index[-1].date(), rain=rain_src,
                         EVP=round(ml.stats.evp(), 1), RMSE_m=round(ml.stats.rmse(), 2),
                         t50_days=round(rm.rfunc.get_tmax(par, cutoff=0.5)), t95_days=round(rm.rfunc.get_tmax(par, cutoff=0.95)),
                         obs_trend_m_per_yr=round(obs_slope, 3), trend_p_value=round(mk_p, 3), model_trend_m_per_yr=round(slope, 3), longest_drought_months=longest, months=len(sgi), drought_months_SGI_below_m1=drought_months,
                         severe_months_SGI_below_m1_5=severe, last_month=last.strftime("%Y-%m"), last_SGI=round(sgi_last, 2), status=status))
        ml.plots.results(figsize=(10, 7)); plt.savefig(os.path.join(FT, f"{name}.png"), dpi=100, bbox_inches="tight"); plt.close("all")
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.bar(sgi.index, sgi.clip(upper=0), width=25, color="#b2182b"); ax.bar(sgi.index, sgi.clip(lower=0), width=25, color="#2166ac")
        for y, c in [(-1, "#ef8a62"), (-1.5, "#b2182b")]: ax.axhline(y, ls="--", c=c, lw=1)
        ax.set_title(f"{name}: Standardised Groundwater Index (SGI), monthly"); ax.set_ylabel("SGI"); ax.set_ylim(-3, 3)
        plt.tight_layout(); plt.savefig(os.path.join(FD, f"{name}_SGI.png"), dpi=100); plt.close("all")
        print("ok", name, rows[-1]["EVP"], slope)
    except Exception as ex:
        print("fail", name, str(ex)[:150])
out = pd.DataFrame(rows).sort_values("EVP", ascending=False)
out.to_csv(os.path.join(R, "processed", "pastas_trend_drought_summary.csv"), index=False)
pd.DataFrame(thr_rows).to_csv(os.path.join(R, "processed", "drought_thresholds_monthly.csv"), index=False)
print(out.to_string())
