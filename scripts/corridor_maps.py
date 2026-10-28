"""Refined cattle corridor map-view figures.

  figures/paper/figC0_corridor_mapview_classes.png   2 x 2 panel (alternative to fig03)
  figures/paper/figC1..C4_corridor_*.png               one full-page map per variable for the whole corridor
Interpolation: inverse distance (k = 15, power 1.5) on a 1.5 km grid, limited to 12 km from a well,
then a light Gaussian smoothing that respects the mask; classed colour scales.
"""
import os, warnings, numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from shapely import contains_xy as contains
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data"); F = os.path.join(H, "..", "figures", "paper")
src = open(os.path.join(H, "cattle_corridor.py")).read(); i = src.index("CC = ["); exec(src[i:src.index("]", i) + 1])
U = "EPSG:32636"; A = os.path.join(D, "raw", "admin")
adm0 = gpd.read_file(os.path.join(A, "uga_admin0.shp")).to_crs(U); adm1 = gpd.read_file(os.path.join(A, "uga_admin1.shp")).to_crs(U)
adm2 = gpd.read_file(os.path.join(A, "uga_admin2.shp")).to_crs(U)
lakes = gpd.read_file(os.path.join(A, "Lakes.shp")).to_crs(U); lakes = lakes[lakes.area > 2e6]
cc = adm2[adm2.adm2_name.isin(CC)]; ccu = cc.union_all()
b = pd.read_csv(os.path.join(D, "processed", "boreholes_located.csv")); bc = b[b.adm2_name.isin(CC)].copy()
st = pd.read_csv(os.path.join(D, "raw", "Stations.csv")).dropna(subset=["latitude"])
st = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.longitude, st.latitude), crs=4326).to_crs(U)
st = st[contains(ccu, st.geometry.x.values, st.geometry.y.values)]
x0, y0, x1, y1 = ccu.bounds; g = 1500
gx, gy = np.meshgrid(np.arange(x0, x1, g), np.arange(y0, y1, g)); inside = contains(ccu, gx, gy)
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
with rasterio.open(os.path.join(D, "raw", "dem", "uganda_dem_220m_utm36n.tif")) as r:
    src_ = r.read(1).astype("float32"); src_[(src_ < 300) | (src_ > 5000)] = np.nan
    DEM = np.full(gx.shape, np.nan, "float32")
    reproject(src_, DEM, src_transform=r.transform, src_crs=r.crs, dst_transform=from_origin(x0 - g / 2, y0 - g / 2, g, -g),
              dst_crs=r.crs, src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.average)
# group stations sharing a site (e.g. Mubende 1-4) for labels
st["lab"] = st.station_name.str.replace(r"\s*\(.*\)", "", regex=True)
grp = st.groupby([st.geometry.x.round(-2), st.geometry.y.round(-2)])
SL = []
for (xx, yy), q in grp:
    names = sorted(q.lab.unique())
    base = sorted(set(n.rstrip(" 0123456789") for n in names))
    lab = names[0] if len(names) == 1 else (f"{base[0]} {min(n[-1] for n in names)}–{max(n[-1] for n in names)}" if len(base) == 1 else " / ".join(names))
    SL.append((xx, yy, lab))

def surface(col, log=False):
    v = bc[col].values.astype(float); ok = np.isfinite(v); px, py, pv = bc.x.values[ok], bc.y.values[ok], v[ok]
    if log: pv = np.log10(np.clip(pv, 0.1, 50))
    d, i = cKDTree(np.c_[px, py]).query(np.c_[gx.ravel(), gy.ravel()], k=15)
    w = 1 / np.maximum(d, 200) ** 1.5; z = ((w * pv[i]).sum(1) / w.sum(1)).reshape(gx.shape)
    m = (d[:, 0] <= 12000).reshape(gx.shape) & inside
    a = gaussian_filter(np.where(m, z, 0), 1.0); wt = gaussian_filter(m.astype(float), 1.0)
    z = np.where(m, a / np.maximum(wt, 1e-6), np.nan)
    return (10 ** z if log else z), ok.sum()

SPECS = {
 "swl": dict(title="Static water level", unit="m below ground", cmap="YlGnBu", bounds=[0, 5, 10, 15, 20, 30, 40, 60]),
 "bedrock": dict(title="Depth to bedrock", unit="m below ground", cmap="YlOrBr", bounds=[0, 10, 20, 30, 40, 50, 70]),
 "yield": dict(title="Borehole yield", unit="m³/h", cmap="viridis", bounds=[0, 0.5, 1, 1.5, 2, 3, 5, 10], log=True),
 "wt": dict(title="Water-table elevation", unit="m above sea level", cmap="Spectral_r", bounds=list(range(1000, 1601, 50))),
}
Z = {}
for k, s in SPECS.items():
    if k == "wt":
        z, n = surface("swl"); Z[k] = (DEM - z, n)
    else:
        Z[k] = surface(k, s.get("log", False))

def scalebar(ax, length_km=50):
    xa, xb = ax.get_xlim(); ya, yb = ax.get_ylim(); X = xa + 0.70 * (xb - xa); Y = ya + 0.03 * (yb - ya); L = length_km * 1000
    for j in range(2):
        ax.add_patch(Rectangle((X + j * L / 2, Y), L / 2, 0.006 * (yb - ya), fc="k" if j == 0 else "w", ec="k", lw=0.6, zorder=6))
    for j, t in enumerate([0, length_km // 2, length_km]):
        ax.text(X + j * L / 2, Y + 0.012 * (yb - ya), f"{t}", ha="center", fontsize=7, zorder=6)
    ax.text(X + L + 0.01 * (xb - xa), Y, "km", fontsize=7, va="bottom", zorder=6)
    ax.annotate("N", xy=(xb - 0.05 * (xb - xa), yb - 0.06 * (yb - ya)), xytext=(xb - 0.05 * (xb - xa), yb - 0.13 * (yb - ya)),
                ha="center", fontsize=10, fontweight="bold", arrowprops=dict(arrowstyle="-|>", color="k", lw=1.5), zorder=6)

def draw(ax, k, big=False):
    s = SPECS[k]; z, n = Z[k]; cmap = plt.get_cmap(s["cmap"], len(s["bounds"])); norm = BoundaryNorm(s["bounds"], cmap.N, extend="max")
    cc.plot(ax=ax, color="#f2f2f2", zorder=0.1)
    im = ax.pcolormesh(gx, gy, z, cmap=cmap, norm=norm, shading="auto", zorder=0.3, rasterized=True)
    if k == "wt":
        cs = ax.contour(gx, gy, gaussian_filter(np.nan_to_num(z, nan=np.nanmean(z)), 1.5) * np.where(np.isfinite(z), 1, np.nan),
                        levels=range(1000, 1700, 100), colors="k", linewidths=0.5, zorder=3)
        ax.clabel(cs, fontsize=6 if big else 5, fmt="%d m")
    lakes.plot(ax=ax, color="#cfe3f3", edgecolor="#8fb4d6", lw=0.3, zorder=2)
    adm2.boundary.plot(ax=ax, color="#bdbdbd", lw=0.25, zorder=2.5); cc.boundary.plot(ax=ax, color="#6e6e6e", lw=0.45, zorder=2.6)
    adm0.boundary.plot(ax=ax, color="k", lw=0.9, zorder=2.7)
    ax.scatter(bc.x, bc.y, s=1.6 if big else 0.6, c="k", alpha=0.45, lw=0, zorder=3.5)
    ax.scatter(st.geometry.x, st.geometry.y, s=40 if big else 16, marker="^", c="#d7301f", edgecolor="w", lw=0.6, zorder=5)
    ax.set_xlim(x0 - 8000, x1 + 8000); ax.set_ylim(y0 - 8000, y1 + 8000)
    if big:
        from adjustText import adjust_text
        T = []
        for _, r in cc.iterrows():
            c = r.geometry.representative_point()
            T.append(ax.text(c.x, c.y, r.adm2_name, fontsize=6.6, ha="center", va="center", zorder=4.5, color="#1f1f1f", style="italic"))
        for xx, yy, lab in SL:
            T.append(ax.text(xx, yy, lab, fontsize=6.8, color="#a3150a", fontweight="bold", zorder=5.5))
        adjust_text(T, ax=ax, expand=(1.15, 1.3), arrowprops=dict(arrowstyle="-", color="#a3150a", lw=0.4), only_move={"text": "xy"})
    ax.set_xlim(x0 - 8000, x1 + 8000); ax.set_ylim(y0 - 8000, y1 + 8000); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    for sp in ax.spines.values(): sp.set_color("#888")
    cb = plt.colorbar(im, ax=ax, shrink=0.55 if big else 0.62, pad=0.01, ticks=s["bounds"] if k != "wt" else s["bounds"][::2])
    cb.set_label(f"{s['title']}, {s['unit']}", fontsize=9 if big else 7.5); cb.ax.tick_params(labelsize=8 if big else 6.5)
    scalebar(ax, 100 if not big else 50)
    return n

# chapter 2 x 2
fig, axs = plt.subplots(2, 2, figsize=(10.5, 10.2))
for ax, k, lab in zip(axs.ravel(), ["swl", "bedrock", "yield", "wt"], "abcd"):
    n = draw(ax, k); ax.set_title(f"({lab}) {SPECS[k]['title']} (n = {n:,} wells)", fontsize=9.5, loc="left")
fig.suptitle("Cattle corridor in map view: borehole data interpolated between wells (1.5 km grid, ≤ 12 km from a well); red triangles: monitoring stations", fontsize=9.5)
plt.tight_layout(rect=(0, 0, 1, 0.97)); plt.savefig(os.path.join(F, "figC0_corridor_mapview_classes.png"), dpi=220); plt.close()

# annex full pages
names = {"swl": "C1_corridor_static_water_level", "bedrock": "C2_corridor_depth_to_bedrock", "yield": "C3_corridor_yield", "wt": "C4_corridor_water_table"}
for k, nm in names.items():
    fig, ax = plt.subplots(figsize=(8.6, 8.6))
    n = draw(ax, k, big=True)
    ax.set_title(f"Uganda cattle corridor: {SPECS[k]['title'].lower()} ({SPECS[k]['unit']}), interpolated from {n:,} boreholes", fontsize=10.5)
    fig.text(0.02, 0.012, "Inverse-distance interpolation (k = 15, power 1.5), 1.5 km grid, blank > 12 km from a borehole. Dots: boreholes; red triangles: monitoring stations.\n"
            "Boundaries: OCHA/UBOS 2026; projection UTM 36N." + (" Water table = DEM − interpolated static level; contours every 100 m." if k == "wt" else ""),
            fontsize=7.2, color="#444", va="bottom")
    plt.tight_layout(rect=(0, 0.03, 1, 1)); plt.savefig(os.path.join(F, f"fig{nm}.png"), dpi=220, bbox_inches="tight", pad_inches=0.15); plt.close()
print("corridor maps done")
