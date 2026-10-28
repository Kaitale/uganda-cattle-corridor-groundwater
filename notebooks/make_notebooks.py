"""Generate the three workflow notebooks (run from the repo root or notebooks/)."""
import os, nbformat as nbf
H = os.path.dirname(os.path.abspath(__file__))
def nb(cells):
    n = nbf.v4.new_notebook(); n.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    n.cells = [nbf.v4.new_markdown_cell(c[1:]) if c.startswith("#MD") else nbf.v4.new_code_cell(c) for c in cells]; return n
SETUP = """import os, sys, warnings
import numpy as np, pandas as pd, matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
ROOT = os.path.abspath(".." if os.path.basename(os.getcwd()) == "notebooks" else ".")
DATA = os.path.join(ROOT, "data"); sys.path.insert(0, os.path.join(ROOT, "scripts"))"""
n1 = nb([
"#MD# 1. National borehole database and district summaries\nStatistics of the MWE national borehole database (12,529 records) and district medians on the 2026 OCHA/UBOS boundaries. Raw data are not in the repository; see `data/README.md`.",
SETUP,
"""a = pd.read_excel(os.path.join(DATA, "raw", "boreholes", "BOREHOLE DATA.xlsx"))
print(f"{len(a):,} records, {a.District.nunique()} district names")
a[["TotalDrilledDepth", "DepthToBedRock", "StaticWaterLevel", "ConstantDischargeYield"]].apply(pd.to_numeric, errors="coerce").describe().round(1)""",
"""import borehole_stats; borehole_stats.main()""",
"#MD## Distributions",
"""fig, axs = plt.subplots(1, 4, figsize=(12, 2.8))
for ax, (c, hi) in zip(axs, [("TotalDrilledDepth", 150), ("DepthToBedRock", 100), ("StaticWaterLevel", 60), ("ConstantDischargeYield", 10)]):
    x = pd.to_numeric(a[c], errors="coerce"); x = x[x.between(0, hi)]
    ax.hist(x, 40, color="#4a7aa8"); ax.axvline(x.median(), c="k", ls="--"); ax.set_title(f"{c}\\nmedian {x.median():.1f}", fontsize=9)
plt.tight_layout()""",
"#MD## District medians (2026 boundaries)\n`make_maps.py` joins records to districts by name and writes `data/processed/district_borehole_summary.csv`.",
"""d = pd.read_csv(os.path.join(DATA, "processed", "district_borehole_summary.csv"))
print(d.groupby("adm1_name")[["median_yield_m3h", "median_bedrock_m"]].median().round(1))
d.dropna().sort_values("median_bedrock_m", ascending=False).head(10)""",
"""from IPython.display import Image
Image(os.path.join(ROOT, "figures", "paper", "fig06_district_boreholes.png"), width=900)"""])
n2 = nb([
"#MD# 2. Groundwater level models, trends and drought indicators\nPastas recharge model (rain − f·evaporation, Gamma response) with a linear trend, Theil–Sen / Mann–Kendall trends on monthly levels, the Standardised Groundwater Index (SGI) and monthly percentile drought thresholds. The full run for all stations is `scripts/pastas_trend_drought.py`; this notebook walks through one station.",
SETUP + "\nimport pastas as ps\nfrom scipy.stats import theilslopes, kendalltau\nps.set_log_level('ERROR')",
"""STATION = "Rwonyo"
st = pd.read_csv(os.path.join(DATA, "raw", "Stations.csv")).set_index("station_name")
d = pd.read_csv(os.path.join(DATA, "raw", "levels", st.loc[STATION, "data_file"]), encoding="latin1", on_bad_lines="skip")
dc = d.columns[0]; d[dc] = pd.to_datetime(d[dc], errors="coerce"); d = d.dropna(subset=[dc]).set_index(dc).sort_index(); d = d[~d.index.duplicated()]
dep = pd.to_numeric(d.iloc[:, 0], errors="coerce").dropna(); dep = dep[(dep - dep.median()).abs() < 5 * dep.std()]
dep.plot(figsize=(10, 3), title=f"{STATION}: depth to water (m bgl)").invert_yaxis()""",
"#MD## Forcing: ERA5-Land rain and evaporation at the station",
"""tp = pd.read_csv(os.path.join(DATA, "processed", "era5_tp_stations.csv"), index_col=0, parse_dates=True)[STATION].dropna()
pev = pd.read_csv(os.path.join(DATA, "processed", "era5_pev_stations.csv"), index_col=0, parse_dates=True)[STATION].dropna()
h = -dep.loc[max(dep.index[0], tp.index[0]):min(dep.index[-1], tp.index[-1])]
ml = ps.Model(h, name=STATION)
ml.add_stressmodel(ps.RechargeModel(tp / 1000, pev / 1000, ps.Gamma(), name="recharge", recharge=ps.rch.Linear()))
ml.add_stressmodel(ps.LinearTrend(h.index[0], h.index[-1], name="trend"))
ml.solve(report=False)
print(f"EVP {ml.stats.evp():.1f}%, RMSE {ml.stats.rmse():.2f} m")
ml.plots.results(figsize=(10, 7));""",
"#MD## Observed trend (Theil–Sen, Mann–Kendall)",
"""mon = dep.resample("MS").mean().dropna()
yrs = (mon.index - mon.index[0]).days / 365.25
slope = -theilslopes(mon.values, yrs)[0]; p = kendalltau(yrs, mon.values).pvalue
print(f"trend {slope:+.3f} m/yr (positive = rising), Mann-Kendall p = {p:.3g}")""",
"#MD## Drought: SGI and monthly thresholds",
"""sgi = ps.stats.sgi(-mon)
ax = sgi.plot(kind="line", figsize=(10, 2.8)); ax.axhline(-1, ls="--", c="orange"); ax.axhline(-1.5, ls="--", c="r"); ax.set_title("SGI")
thr = mon.groupby(mon.index.month).quantile([0.5, 0.8, 0.9]).unstack().round(2)
thr.columns = ["normal_P50", "moderate_P80", "severe_P90"]; thr""",
"#MD## All stations\nResults written by `pastas_trend_drought.py`:",
"""s = pd.read_csv(os.path.join(DATA, "processed", "pastas_trend_drought_summary.csv"))
s[["station", "EVP", "t50_days", "obs_trend_m_per_yr", "trend_p_value", "last_SGI", "status"]].round(2).sort_values("EVP", ascending=False)"""])
n3 = nb([
"#MD# 3. Aquifer geometry from drilling logs: regolith maps and GemPy\nDrilling logs from the Aswa catchment are classified into cover, weathered/fractured rock and fresh rock (`litho_prep.py`), interpolated to maps (`make_maps.py`) and modelled in 3D with GemPy for Kitgum–Pader (`gempy_kitgum.py`).",
SETUP,
"""import litho_prep
for s in ["Top soil", "Lateritic clay", "Weathered granite", "Fractured gneiss", "Fresh granite"]:
    print(f"{s:20s} -> {litho_prep.classify(s)}")""",
"""p = pd.read_csv(os.path.join(DATA, "processed", "aswa_layer_picks_dem.csv"))
p = p[(p.fresh_top >= p.cover_base) & (p.fresh_top < 120)]
print(len(p), "consistent logs"); p[["cover_base", "fresh_top", "weathered_thk", "ground_elev_m"]].describe().round(1)""",
"""p.groupby("district")[["cover_base", "fresh_top", "weathered_thk"]].median().round(1).sort_values("fresh_top")""",
"""from IPython.display import Image
Image(os.path.join(ROOT, "figures", "paper", "fig07_aswa_regolith_maps.png"), width=900)""",
"#MD## GemPy model\nRun `python scripts/gempy_kitgum.py` (several minutes). The W–E section is masked above the DEM and more than 6 km from any borehole.",
"""Image(os.path.join(ROOT, "figures", "paper", "fig08_gempy_kitgum_section.png"), width=900)"""])
n4 = nb([
"#MD# 4. Cattle corridor: wells in 2D map view and 3D\nLocates national-database boreholes, interpolates static water level, bedrock depth and yield in map view, and plots wells with water-table and bedrock surfaces in 3D. Full script: `scripts/cattle_corridor.py`.",
SETUP,
"""b = pd.read_csv(os.path.join(DATA, "processed", "boreholes_located.csv"))
src = open(os.path.join(ROOT, "scripts", "cattle_corridor.py")).read(); i = src.index("CC = ["); exec(src[i:src.index("]", i) + 1])
b["corridor"] = b.adm2_name.isin(CC)
print(f"located boreholes: {len(b):,}; in the cattle corridor ({len(CC)} districts): {b.corridor.sum():,}")
b.groupby("corridor")[["swl", "bedrock", "depth", "yield"]].median().round(1)""",
"""pd.read_csv(os.path.join(DATA, "processed", "corridor_district_summary.csv")).sort_values("median_swl_m", ascending=False).head(12)""",
"""from IPython.display import Image, IFrame
Image(os.path.join(ROOT, "figures", "paper", "fig03_corridor_mapview.png"), width=900)""",
"""Image(os.path.join(ROOT, "figures", "paper", "fig05_corridor_3d_wells.png"), width=1000)""",
"#MD## Groundwater atlas: Uganda and the corridor sub-regions\n`scripts/regional_maps.py` zooms from Uganda to each sub-region (as for the Aswa catchment).",
"""print(pd.read_csv(os.path.join(DATA, "processed", "region_summary.csv")).round(2).to_string(index=False))""",
"""for f in ["figB1_uganda_mapview.png", "figB2_central_mapview.png", "figB3_southwest_mapview.png", "figB4_teso_mapview.png", "figB5_karamoja_mapview.png"]:
    display(Image(os.path.join(ROOT, "figures", "paper", f), width=900))""",
"#MD## Full-page maps\nOne map per area and variable is stored in `figures/maps/<area>/` (Uganda, cattle corridor, central, south-west, Teso, Karamoja).",
"""import glob
for d in sorted(glob.glob(os.path.join(ROOT, "figures", "maps", "*"))):
    print(os.path.basename(d), [os.path.basename(f) for f in sorted(glob.glob(os.path.join(d, "*.png")))])
display(Image(os.path.join(ROOT, "figures", "maps", "central", "central_static_water_level.png"), width=800))""",
"#MD## Interactive 3D\n`figures/interactive/kai_cattle_corridor_3d.html` (KAI Water & Energy showcase): three areas, hover over a well for its ID, district, layer type and depth interval. `corridor_3d_wells.html` is the plain version."])
for name, n in [("01_borehole_database.ipynb", n1), ("02_pastas_trends_drought.ipynb", n2), ("03_regolith_gempy.ipynb", n3), ("04_cattle_corridor_2d_3d.ipynb", n4)]:
    nbf.write(n, os.path.join(H, name)); print("wrote", name)
