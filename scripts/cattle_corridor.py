"""Cattle corridor: well maps, 2D map-view interpolation and 3D well plots.

Inputs : data/processed/boreholes_located.csv (borehole_locate.py), pastas_trend_drought_summary.csv,
         data/raw/admin/*.shp, data/raw/dem/uganda_dem_220m_utm36n.tif
Outputs: figures/paper/fig01_cattle_corridor.png      corridor, stations and trends, boreholes
         figures/paper/fig03_corridor_mapview.png     interpolated SWL, bedrock depth, yield, water-table contours
         figures/paper/fig05_corridor_3d_wells.png    3D wells with ground, water-table and bedrock surfaces
         figures/interactive/corridor_3d_wells.html  interactive 3D (plotly)
         data/processed/corridor_district_summary.csv
"""
import os, warnings, numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from shapely import contains_xy as contains
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data")
F = os.path.join(H, "..", "figures", "paper"); FI = os.path.join(H, "..", "figures", "interactive")
os.makedirs(F, exist_ok=True); os.makedirs(FI, exist_ok=True)
plt.rcParams.update({"font.size": 9})
U = "EPSG:32636"
# Cattle corridor districts (2026 boundaries): the semi-arid rangeland belt from Isingiro/Ntungamo through
# the central rangelands and Teso to Karamoja (after MAAIF and Mwiru 2015, updated to new districts)
CC = ["Isingiro", "Ntungamo", "Mbarara", "Kiruhura", "Kazo", "Ibanda", "Lyantonde", "Ssembabule", "Rakai", "Kyotera", "Lwengo",
      "Gomba", "Mubende", "Kassanda", "Kiboga", "Kyankwanzi", "Nakaseke", "Luwero", "Nakasongola", "Kayunga", "Kiryandongo",
      "Masindi", "Apac", "Kwania", "Amolatar", "Kaberamaido", "Kalaki", "Soroti", "Serere", "Katakwi", "Amuria", "Kapelebyong",
      "Kumi", "Ngora", "Bukedea", "Nakapiripirit", "Nabilatuk", "Amudat", "Napak", "Moroto", "Kotido", "Abim", "Kaabong", "Karenga"]
A = os.path.join(D, "raw", "admin")
adm0 = gpd.read_file(os.path.join(A, "uga_admin0.shp")).to_crs(U); adm1 = gpd.read_file(os.path.join(A, "uga_admin1.shp")).to_crs(U)
adm2 = gpd.read_file(os.path.join(A, "uga_admin2.shp")).to_crs(U)
lakes = gpd.read_file(os.path.join(A, "Lakes.shp")).to_crs(U); lakes = lakes[lakes.area > 1e6]
LK = dict(color="#cfe3f3", edgecolor="#8fb4d6", lw=0.3)
cc = adm2[adm2.adm2_name.isin(CC)]; ccu = cc.union_all()
b = pd.read_csv(os.path.join(D, "processed", "boreholes_located.csv"))
b["in_cc"] = b.adm2_name.isin(CC)
b["wt_elev"] = b.ground_elev_m - b.swl; b["rock_elev"] = b.ground_elev_m - b.bedrock
bc = b[b.in_cc].copy()
km = lambda v: np.asarray(v) / 1000

def idw(px, py, pv, gx, gy, k=10, maxd=12000, power=2):
    ok = np.isfinite(pv); tr = cKDTree(np.c_[px[ok], py[ok]]); v = pv[ok]
    d, i = tr.query(np.c_[gx.ravel(), gy.ravel()], k=k); w = 1 / np.maximum(d, 100) ** power
    z = (w * v[i]).sum(1) / w.sum(1); z[d[:, 0] > maxd] = np.nan
    return z.reshape(gx.shape)

def frame(ax, xlim=None, ylim=None, labels=True):
    lakes.plot(ax=ax, zorder=0.5, **LK)
    adm2.boundary.plot(ax=ax, color="#d0d0d0", lw=0.25, zorder=1); adm0.boundary.plot(ax=ax, color="k", lw=0.8, zorder=2)
    ax.set_aspect("equal")
    if xlim: ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")

# ---------- Fig 1: corridor overview ----------
st = pd.read_csv(os.path.join(D, "raw", "Stations.csv")).dropna(subset=["latitude"])
sm = pd.read_csv(os.path.join(D, "processed", "pastas_trend_drought_summary.csv"))
st = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.longitude, st.latitude), crs=4326).to_crs(U)
st = st.merge(sm[["station", "obs_trend_m_per_yr", "trend_p_value", "status"]], left_on="station_name", right_on="station", how="left")
st["in_cc"] = contains(ccu, st.geometry.x, st.geometry.y)
fig, ax = plt.subplots(figsize=(6.4, 6.6))
cc.plot(ax=ax, color="#f3e3c3", zorder=0.2)
frame(ax)
ax.scatter(km(b.x[~b.in_cc]) * 1000, km(b.y[~b.in_cc]) * 1000, s=0.6, c="#9a9a9a", alpha=0.5, zorder=2.5, label=f"Borehole, outside corridor ({(~b.in_cc).sum():,})")
ax.scatter(bc.x, bc.y, s=1.2, c="#7a4b16", alpha=0.7, zorder=2.6, label=f"Borehole, cattle corridor ({len(bc):,})")
m = st.dropna(subset=["obs_trend_m_per_yr"]); sig = m.trend_p_value < 0.05
sc = ax.scatter(m.geometry.x, m.geometry.y, c=m.obs_trend_m_per_yr.clip(-0.5, 0.5), cmap="RdBu", vmin=-0.5, vmax=0.5,
                s=np.where(sig, 85, 35), edgecolor="k", lw=0.5, zorder=4, marker="o")
n = st[st.obs_trend_m_per_yr.isna()]
ax.scatter(n.geometry.x, n.geometry.y, s=22, marker="^", c="w", edgecolor="k", lw=0.5, zorder=4, label="Station, < 3 years of data")
for _, r in m.iterrows():
    if r.station_name in ["Rwonyo", "Isingiro", "Mubende 2", "Kiboga", "Nakasongola", "Bombo Lukole", "Soroti", "Katakwi", "Kangole", "Rakai", "Mbarara"]:
        ax.annotate(r.station_name, (r.geometry.x, r.geometry.y), xytext=(5, 2), textcoords="offset points", fontsize=6.5, zorder=5,
                    bbox=dict(fc="w", ec="none", alpha=0.6, pad=0.3))
cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.02); cb.set_label("Groundwater level trend, m/yr\n(large: p < 0.05; blue rising, red falling)", fontsize=8)
lg = ax.legend(loc="lower right", fontsize=6.5, framealpha=0.9)
for h, sz in zip(lg.legend_handles, [12, 12, 30]): h.set_sizes([sz])
ax.set_title(f"Uganda's cattle corridor ({len(cc)} districts, shaded): monitoring stations and boreholes", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(F, "fig01_cattle_corridor.png"), dpi=200); plt.close()
print("stations in corridor", int(st.in_cc.sum()), "of", len(st), "; with trends", int((st.in_cc & st.obs_trend_m_per_yr.notna()).sum()))

# ---------- Fig 4: 2D map view ----------
x0, y0, x1, y1 = ccu.bounds; g = 2500
gx, gy = np.meshgrid(np.arange(x0, x1, g), np.arange(y0, y1, g)); inside = contains(ccu, gx, gy)
with rasterio.open(os.path.join(D, "raw", "dem", "uganda_dem_220m_utm36n.tif")) as r:
    DEMG = np.array([v[0] for v in r.sample(zip(gx.ravel(), gy.ravel()))], float).reshape(gx.shape)
DEMG[(DEMG < 300) | (DEMG > 5000)] = np.nan
from scipy.ndimage import uniform_filter
panels = [("swl", "Static water level, m bgl", "YlGnBu", (0, 40)), ("bedrock", "Depth to bedrock, m", "YlOrBr", (0, 60)),
          ("yield", "Yield, m³/h", "viridis", (0, 5)), ("wt_elev", "Water-table elevation, m a.s.l.", "terrain", None)]
fig, axs = plt.subplots(2, 2, figsize=(10.5, 10.5))
for ax, (col, lab, cmap, vr) in zip(axs.ravel(), panels):
    v = bc[("swl" if col == "wt_elev" else col)].values.astype(float); ok = np.isfinite(v)
    if col == "yield": v = np.log10(np.clip(v, 0.1, 50))       # interpolate log-yield
    z = idw(bc.x.values[ok], bc.y.values[ok], v[ok], gx, gy); z[~inside] = np.nan
    if col == "wt_elev": z = DEMG - z                            # water table = ground (DEM) - interpolated SWL
    if col == "yield": z = 10 ** z
    cc.plot(ax=ax, color="#f4f4f4", zorder=0.1)
    if vr is None: vr = (np.nanpercentile(z, 2), np.nanpercentile(z, 98))
    im = ax.pcolormesh(gx, gy, z, cmap=cmap, vmin=vr[0], vmax=vr[1], shading="auto", zorder=0.3)
    frame(ax); ax.set_xlim(x0 - 5000, x1 + 5000); ax.set_ylim(y0 - 5000, y1 + 5000)
    if col == "wt_elev":
        cs = ax.contour(gx, gy, z, levels=np.arange(900, 1700, 100), colors="k", linewidths=0.4, zorder=3)
        ax.clabel(cs, fontsize=5, fmt="%d")
    ax.scatter(bc.x[ok], bc.y[ok], s=0.8, c="k", alpha=0.45, zorder=3.5)
    plt.colorbar(im, ax=ax, shrink=0.6, pad=0.01).set_label(lab, fontsize=8)
    ax.set_title(lab + f" (n = {ok.sum():,})", fontsize=9)
fig.suptitle("Cattle corridor, map view: wells interpolated by inverse distance (2.5 km grid, ≤12 km from a well).\n"
             "Water-table elevation = DEM − interpolated static water level (100 m contours).", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(F, "fig03_corridor_mapview.png"), dpi=200); plt.close()

# district summary
ds = bc.groupby("adm2_name").agg(boreholes=("x", "size"), median_swl_m=("swl", "median"), median_bedrock_m=("bedrock", "median"),
                                  median_yield_m3h=("yield", "median"), median_depth_m=("depth", "median")).round(2)
ds.to_csv(os.path.join(D, "processed", "corridor_district_summary.csv"))
print(bc[["swl", "bedrock", "yield", "depth"]].median().round(2).to_dict(), "low-yield share <1 m3/h", round((bc["yield"] < 1).mean(), 3))

# ---------- Fig 5: 3D wells in depth below ground ----------
from mpl_toolkits.mplot3d import Axes3D  # noqa
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
BLOCKS = [("Central corridor: Nakasongola and Luwero", ["Nakasongola", "Luwero"], ["Nakasongola", "Bombo Lukole"]),
          ("South-west corridor: Isingiro, Mbarara and Kiruhura", ["Isingiro", "Mbarara", "Kiruhura", "Kazo"], ["Rwonyo", "Isingiro", "Mbarara"])]
from scipy.ndimage import gaussian_filter
def smooth(Z, sig=1.2):
    m = np.isfinite(Z); a = gaussian_filter(np.where(m, Z, 0), sig); w_ = gaussian_filter(m.astype(float), sig)
    return np.where(m, a / np.maximum(w_, 1e-6), np.nan)
def block(ax, dists, stns, title, g3=2000, maxd=12000, nmax=170):
    w = bc[bc.adm2_name.isin(dists)].dropna(subset=["depth"]).copy()
    reg = cc[cc.adm2_name.isin(dists)].union_all()
    X, Y = np.meshgrid(np.arange(*reg.bounds[0::2], g3), np.arange(*reg.bounds[1::2], g3)); inr = contains(reg, X, Y)
    S = smooth(idw(w.x.values, w.y.values, w.swl.values, X, Y, k=15, maxd=maxd, power=1.5)); S[~inr] = np.nan
    B = smooth(idw(w.x.values, w.y.values, w.bedrock.values, X, Y, k=15, maxd=maxd, power=1.5)); B[~inr] = np.nan
    ax.plot_surface(km(X), km(Y), np.where(inr, 0, np.nan), color="#d9c49a", alpha=0.25, linewidth=0)
    ax.plot_surface(km(X), km(Y), -S, color="#2f7bd8", alpha=0.55, linewidth=0)
    ax.plot_surface(km(X), km(Y), -B, color="#707070", alpha=0.55, linewidth=0)
    for _, r in w.sample(min(nmax, len(w)), random_state=2).iterrows():
        xb, yb = r.x / 1000, r.y / 1000; bd = r.bedrock if np.isfinite(r.bedrock) else np.nan
        if np.isfinite(bd):
            ax.plot([xb, xb], [yb, yb], [0, -bd], c="#b07a2a", lw=1.0); ax.plot([xb, xb], [yb, yb], [-bd, -r.depth], c="#2b2b2b", lw=0.7)
        else:
            ax.plot([xb, xb], [yb, yb], [0, -r.depth], c="#9a9a9a", lw=0.7)
        if np.isfinite(r.swl): ax.scatter([xb], [yb], [-r.swl], s=4, c="#0b4ea2", depthshade=False)
    for n in stns:
        q = st[st.station_name == n]
        if len(q):
            gx_, gy_ = q.geometry.x.iloc[0] / 1000, q.geometry.y.iloc[0] / 1000
            ax.plot([gx_, gx_], [gy_, gy_], [0, 25], c="r", lw=1.5); ax.text(gx_, gy_, 28, n, color="r", fontsize=7)
    ax.set_zlim(-120, 30); ax.set_box_aspect((np.ptp(X), np.ptp(Y), 0.45 * max(np.ptp(X), np.ptp(Y))))
    ax.set_xlabel("Easting, km", fontsize=7, labelpad=-2); ax.set_ylabel("Northing, km", fontsize=7, labelpad=-2)
    ax.set_zlabel("Depth below ground, m", fontsize=7); ax.tick_params(labelsize=6, pad=0)
    ax.set_title(f"{title} ({min(nmax, len(w))} of {len(w)} wells shown)", fontsize=8.5); ax.view_init(elev=24, azim=-58)
    return w, X, Y, S, B
fig = plt.figure(figsize=(12, 6.2))
res = []
for k, (t, d_, s_) in enumerate(BLOCKS):
    res.append(block(fig.add_subplot(1, 2, k + 1, projection="3d"), d_, s_, t))
fig.legend(handles=[Patch(color="#2f7bd8", alpha=0.6, label="Water table (interpolated static level)"), Patch(color="#707070", alpha=0.6, label="Top of bedrock (interpolated)"),
                    Line2D([], [], c="#b07a2a", lw=2, label="Well in regolith"), Line2D([], [], c="#2b2b2b", lw=2, label="Well in bedrock"),
                    Line2D([], [], ls="", marker="o", c="#0b4ea2", ms=4, label="Static water level"), Line2D([], [], c="r", lw=2, label="Monitoring station")],
           loc="lower center", ncol=6, fontsize=7.5, frameon=False)
fig.suptitle("Cattle corridor wells in 3D: regolith aquifer, water table and top of bedrock, interpolated between wells (depth below ground)", fontsize=9.5)
plt.subplots_adjust(left=0, right=1, top=0.92, bottom=0.08, wspace=0.02)
plt.savefig(os.path.join(F, "fig05_corridor_3d_wells.png"), dpi=200); plt.close()
w, X, Y, S, B = res[0]
G = np.where(np.isfinite(S), 0.0, np.nan); WT = -S; RK = -B; ry0, ry1, rx0, rx1 = Y.min(), Y.max(), X.min(), X.max()
w["ground_elev_m"] = 0.0
# interactive 3D
try:
    import plotly.graph_objects as go
    tr = [go.Surface(x=km(X), y=km(Y), z=G, colorscale=[[0, "#c9b48a"], [1, "#c9b48a"]], opacity=0.35, showscale=False, name="Ground"),
          go.Surface(x=km(X), y=km(Y), z=WT, colorscale=[[0, "#3b7dd8"], [1, "#3b7dd8"]], opacity=0.55, showscale=False, name="Water table"),
          go.Surface(x=km(X), y=km(Y), z=RK, colorscale=[[0, "#6d6d6d"], [1, "#6d6d6d"]], opacity=0.7, showscale=False, name="Top of bedrock")]
    for seg, colr, nm in [((0, "bedrock"), "#b07a2a", "Regolith"), (("bedrock", "depth"), "#333333", "Bedrock section")]:
        xs, ys, zs = [], [], []
        for _, r in w.dropna(subset=["bedrock"]).iterrows():
            top = r.ground_elev_m - (0 if seg[0] == 0 else r[seg[0]]); bot = r.ground_elev_m - r[seg[1]]
            xs += [r.x / 1000] * 2 + [None]; ys += [r.y / 1000] * 2 + [None]; zs += [top, bot, None]
        tr.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=colr, width=4), name=nm))
    ww = w.dropna(subset=["swl"])
    tr.append(go.Scatter3d(x=km(ww.x), y=km(ww.y), z=ww.ground_elev_m - ww.swl, mode="markers", marker=dict(size=2.5, color="#1f6fd1"),
                           name="Static water level", text=ww.SourceID + " (" + ww.adm2_name + ")"))
    f3 = go.Figure(tr); f3.update_layout(title="Cattle corridor, Uganda (Nakasongola and Luwero): wells, water table and top of bedrock, depth below ground (m)",
                                         scene=dict(aspectmode="manual", aspectratio=dict(x=1, y=(ry1 - ry0) / (rx1 - rx0), z=0.5),
                                                    xaxis_title="Easting km (UTM 36N)", yaxis_title="Northing km", zaxis_title="depth below ground, m"),
                                         margin=dict(l=0, r=0, t=40, b=0), legend=dict(x=0, y=1))
    f3.write_html(os.path.join(FI, "corridor_3d_wells.html"), include_plotlyjs="cdn")
except ImportError:
    print("plotly not installed: interactive 3D skipped")
print("3D wells", len(w))
