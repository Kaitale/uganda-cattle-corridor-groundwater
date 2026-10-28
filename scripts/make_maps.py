"""National maps with Uganda administrative boundaries (OCHA/UBOS COD-AB, HDX): Fig. 5, 6 and S1.

Needs data/raw/admin/uga_admin{0,1,2}.shp. Produces:
  figS1_station_trends.png    national stations coloured by observed level trend
  fig07_aswa_regolith_maps.png IDW depth to fresh rock and weathered thickness, Aswa
  fig06_district_boreholes.png district medians of yield and depth to bedrock
"""
import os, re, warnings, numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from shapely.geometry import Point
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data")
F = os.path.join(H, "..", "figures", "paper"); os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9})
A = os.path.join(D, "raw", "admin")
adm0 = gpd.read_file(os.path.join(A, "uga_admin0.shp"))
adm1 = gpd.read_file(os.path.join(A, "uga_admin1.shp"))
adm2 = gpd.read_file(os.path.join(A, "uga_admin2.shp"))
lakes = gpd.read_file(os.path.join(A, "Lakes.shp")).to_crs(4326)
lakes = lakes[lakes.geometry.area > 1e-4]  # drop ponds < ~1 km²
LK = dict(color="#cfe3f3", edgecolor="#8fb4d6", lw=0.3)

def base(ax, crs=None, labels=False):
    a0, a1, a2 = (adm0, adm1, adm2) if crs is None else (adm0.to_crs(crs), adm1.to_crs(crs), adm2.to_crs(crs))
    (lakes if crs is None else lakes.to_crs(crs)).plot(ax=ax, zorder=0.5, **LK)
    a2.boundary.plot(ax=ax, color="#c8c8c8", lw=0.3, zorder=1)
    a1.boundary.plot(ax=ax, color="#8a8a8a", lw=0.6, zorder=1)
    a0.boundary.plot(ax=ax, color="k", lw=0.9, zorder=2)
    return a2

# ---------- Fig 1: stations ----------
st = pd.read_csv(os.path.join(D, "raw", "Stations.csv"))
sm = pd.read_csv(os.path.join(D, "processed", "pastas_trend_drought_summary.csv"))
m = sm.merge(st[["station_name", "latitude", "longitude"]], left_on="station", right_on="station_name")
fig, ax = plt.subplots(figsize=(6, 6.2))
base(ax)
for _, r in adm1.iterrows():
    c = r.geometry.representative_point(); ax.text(c.x, c.y, r.adm1_name.upper(), fontsize=8, color="#9a9a9a", ha="center", style="italic", zorder=1)
allst = st.dropna(subset=["latitude"])
ax.scatter(allst.longitude, allst.latitude, s=12, c="#bbbbbb", edgecolor="#777", lw=0.3, label="Record < 3 years", zorder=3)
sig = m.trend_p_value < 0.05
sc = ax.scatter(m.longitude, m.latitude, c=m.obs_trend_m_per_yr.clip(-0.5, 0.5), cmap="RdBu", vmin=-0.5, vmax=0.5,
                s=np.where(sig, 75, 30), edgecolor="k", lw=0.4, zorder=4)
for _, r in m.iterrows():
    if r.station in ["Rwonyo", "Soroti", "Kitgum", "Pader", "Bombo Lukole", "Isingiro", "Masindi", "Kangole", "Galilaya", "Mubende 2"]:
        off = {"Isingiro": (6, -11), "Galilaya": (6, 4), "Rwonyo": (6, 3)}.get(r.station, (4, 3))
        ax.annotate(r.station, (r.longitude, r.latitude), xytext=off, textcoords="offset points", fontsize=7, zorder=5)
cb = plt.colorbar(sc, ax=ax, shrink=0.6); cb.set_label("Level trend, m/yr (blue rising, red falling)")
ax.scatter([], [], s=75, c="w", edgecolor="k", label="Significant, p < 0.05"); ax.scatter([], [], s=30, c="w", edgecolor="k", label="Not significant")
ax.legend(loc="lower right", fontsize=7, frameon=True, framealpha=0.9)
ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)"); ax.set_aspect("equal")
ax.set_title("Groundwater monitoring stations and observed level trends (districts and regions, 2026)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(F, "figS1_station_trends.png"), dpi=200); plt.close()

# ---------- Fig 5: Aswa regolith ----------
U = "EPSG:32636"
p = pd.read_csv(os.path.join(D, "processed", "aswa_layer_picks_dem.csv")).dropna(subset=["fresh_top"])
p = p[(p.fresh_top >= p.cover_base) & (p.fresh_top < 120) & (p.x > 1e5) & (p.y > 1e5)]
pg = gpd.GeoDataFrame(p, geometry=gpd.points_from_xy(p.x, p.y), crs=U)
a2u = adm2.to_crs(U)
pg = gpd.sjoin(pg, a2u[["adm2_name", "geometry"]], predicate="within").drop(columns="index_right")
cnt = pg.adm2_name.value_counts(); dist = a2u[a2u.adm2_name.isin(cnt[cnt >= 10].index)]
mask = dist.union_all()
x0, y0, x1, y1 = mask.bounds; g = 2000
gx, gy = np.meshgrid(np.arange(x0, x1, g), np.arange(y0, y1, g))
tr = cKDTree(pg[["x", "y"]].values); dd, ii = tr.query(np.c_[gx.ravel(), gy.ravel()], k=12)
w = 1 / np.maximum(dd, 50) ** 2
inside = np.array([mask.contains(Point(a, b)) for a, b in zip(gx.ravel(), gy.ravel())]) & (dd[:, 0] < 15000)
fig, axs = plt.subplots(1, 2, figsize=(11, 5.4))
for ax, col, lab, cm, vmax in [(axs[0], "fresh_top", "Depth to fresh rock, m", "YlOrBr", 50), (axs[1], "weathered_thk", "Weathered/fractured zone thickness, m", "PuBu", 30)]:
    v = pg[col].values[ii]; z = (w * v).sum(1) / w.sum(1); z[~inside] = np.nan
    im = ax.pcolormesh(gx / 1000, gy / 1000, z.reshape(gx.shape), cmap=cm, vmin=0, vmax=vmax, shading="auto", zorder=0)
    lakes.to_crs(U).plot(ax=ax, zorder=0.5, **LK)
    a2u.boundary.plot(ax=ax, color="#9a9a9a", lw=0.4, zorder=1)
    adm0.to_crs(U).boundary.plot(ax=ax, color="k", lw=0.9, zorder=2)
    ax.scatter(pg.x / 1000, pg.y / 1000, s=1.2, c="k", alpha=0.35, zorder=3)
    for _, r in dist.iterrows():
        c = r.geometry.representative_point(); ax.text(c.x / 1000, c.y / 1000, r.adm2_name, fontsize=6.5, ha="center", zorder=4,
                                                         bbox=dict(fc="w", ec="none", alpha=0.6, pad=0.5))
    ax.set_xlim(x0 / 1000 - 5, x1 / 1000 + 5); ax.set_ylim(y0 / 1000 - 5, y1 / 1000 + 5); ax.set_aspect("equal")
    ax.set_xlabel("Easting, km (UTM 36N)"); plt.colorbar(im, ax=ax, shrink=0.75, label=lab)
for ax in axs: ax.set_ylabel("Northing, km")
# section line for Fig 6 (GemPy model, Kitgum–Pader, row 30 of 60)
k = pg[pg.district.isin(["KITGUM", "PADER"])]
k = k[(k.cover_base > 0) & (k.fresh_top < 100)]
c2 = (k.x // 2000).astype(int).astype(str) + "_" + (k.y // 2000).astype(int).astype(str)
s = k.groupby(c2)[["x", "y"]].median()
ys = np.linspace(s.y.min() - 1000, s.y.max() + 1000, 60)[30]
for ax in axs:
    ax.plot([k.x.min() / 1000 - 5, k.x.max() / 1000 + 5], [ys / 1000] * 2, "r-", lw=1.4, zorder=5)
for ax in axs: ax.text(k.x.min() / 1000 - 4, ys / 1000 + 4, "Section, Fig. 6", color="r", fontsize=7)
fig.suptitle(f"Aswa catchment regolith from {len(pg):,} borehole logs (inverse-distance interpolation, 2 km grid)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(F, "fig07_aswa_regolith_maps.png"), dpi=200); plt.close()
print("fig5 boreholes", len(pg), "districts", len(dist))

# ---------- Fig 7: district borehole medians ----------
b = pd.read_excel(os.path.join(D, "raw", "boreholes", "BOREHOLE DATA.xlsx"))
norm = lambda s: re.sub(r"[^a-z]", "", str(s).lower().replace("district", "").replace("city", ""))
b["key"] = b.District.map(norm)
for c in ["ConstantDischargeYield", "DepthToBedRock", "TotalDrilledDepth"]:
    b[c] = pd.to_numeric(b[c], errors="coerce")
b.loc[~b.ConstantDischargeYield.between(0, 50), "ConstantDischargeYield"] = np.nan
b.loc[~b.DepthToBedRock.between(0, 200), "DepthToBedRock"] = np.nan
agg = b.groupby("key").agg(n=("District", "size"), yld=("ConstantDischargeYield", "median"), bed=("DepthToBedRock", "median")).reset_index()
a2 = adm2.copy(); a2["key"] = a2.adm2_name.map(norm)
j = a2.merge(agg, on="key", how="left"); j.loc[j.n < 10, ["yld", "bed"]] = np.nan
matched = b.key.isin(a2.key).mean()
fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
for ax, col, lab, cm, vr in [(axs[0], "yld", "Median yield, m³/h", "viridis", (0.5, 4)), (axs[1], "bed", "Median depth to bedrock, m", "YlOrBr", (5, 50))]:
    j.plot(column=col, ax=ax, cmap=cm, vmin=vr[0], vmax=vr[1], edgecolor="#ffffff", lw=0.25,
           missing_kwds=dict(color="#eeeeee", hatch="///", edgecolor="#cccccc", lw=0.2), legend=True,
           legend_kwds=dict(shrink=0.65, label=lab))
    lakes.plot(ax=ax, zorder=2, **LK)
    adm1.boundary.plot(ax=ax, color="#555", lw=0.6, zorder=3); adm0.boundary.plot(ax=ax, color="k", lw=0.9, zorder=3)
    ax.set_aspect("equal"); ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
fig.suptitle(f"National borehole database by district (median of districts with ≥10 records; hatched: no or few data). {matched:.0%} of records matched to 2026 districts", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(F, "fig06_district_boreholes.png"), dpi=200); plt.close()
j[["adm2_name", "adm2_pcode", "adm1_name", "n", "yld", "bed"]].rename(columns=dict(n="boreholes", yld="median_yield_m3h", bed="median_bedrock_m")).to_csv(
    os.path.join(D, "processed", "district_borehole_summary.csv"), index=False)
print("matched", round(matched, 3), "districts with data", int(j.yld.notna().sum()))
print("top yield", j.nlargest(5, "yld")[["adm2_name", "yld"]].values.tolist())
print("low yield", j.nsmallest(5, "yld")[["adm2_name", "yld"]].values.tolist())
