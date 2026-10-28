"""Map-view groundwater atlas: Uganda and the four cattle corridor sub-regions.

For each area: static water level, depth to bedrock, yield and water-table elevation, interpolated
between located boreholes (inverse distance, k = 15, power 1.5, light Gaussian smoothing), with
districts, lakes, monitoring stations, scale bar and north arrow.

Outputs (figures/paper):
  fig03_corridor_mapview.png     whole corridor, 2 x 2 (chapter)
  figB1_uganda_mapview.png       Uganda, 2 x 2 (report annex)
  figB2..B5_<region>_mapview.png  corridor sub-regions, 2 x 2 (report annex)
  figures/maps/<area>/<area>_<variable>.png   one full-page map per area and variable (24 maps)
  data/processed/region_summary.csv
"""
import os, warnings, numpy as np, pandas as pd, geopandas as gpd, rasterio, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from shapely import contains_xy as contains
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from adjustText import adjust_text
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data"); F = os.path.join(H, "..", "figures", "paper")
src = open(os.path.join(H, "cattle_corridor.py")).read(); i = src.index("CC = ["); exec(src[i:src.index("]", i) + 1])
SUB = {"Central rangelands": ["Nakasongola", "Luwero", "Nakaseke", "Kiboga", "Kyankwanzi", "Mubende", "Kassanda", "Kayunga", "Gomba", "Masindi", "Kiryandongo"],
       "South-west rangelands": ["Isingiro", "Mbarara", "Kiruhura", "Kazo", "Ntungamo", "Ibanda", "Lyantonde", "Ssembabule", "Rakai", "Kyotera", "Lwengo"],
       "Teso and Kyoga plains": ["Soroti", "Serere", "Katakwi", "Amuria", "Kapelebyong", "Kumi", "Ngora", "Bukedea", "Kaberamaido", "Kalaki", "Apac", "Kwania", "Amolatar"],
       "Karamoja": ["Nakapiripirit", "Nabilatuk", "Amudat", "Napak", "Moroto", "Kotido", "Abim", "Kaabong", "Karenga"]}
U = "EPSG:32636"; A = os.path.join(D, "raw", "admin")
adm0 = gpd.read_file(os.path.join(A, "uga_admin0.shp")).to_crs(U); adm1 = gpd.read_file(os.path.join(A, "uga_admin1.shp")).to_crs(U)
adm2 = gpd.read_file(os.path.join(A, "uga_admin2.shp")).to_crs(U)
lakes = gpd.read_file(os.path.join(A, "Lakes.shp")).to_crs(U); lakes = lakes[lakes.area > 1e6]
b = pd.read_csv(os.path.join(D, "processed", "boreholes_located.csv"))
st = pd.read_csv(os.path.join(D, "raw", "Stations.csv")).dropna(subset=["latitude"])
st = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.longitude, st.latitude), crs=4326).to_crs(U)
st["lab"] = st.station_name.str.replace(r"\s*\(.*\)", "", regex=True)
SL = []
for (xx, yy), q in st.groupby([st.geometry.x.round(-2), st.geometry.y.round(-2)]):
    names = sorted(q.lab.unique()); base = sorted(set(n.rstrip(" 0123456789") for n in names))
    SL.append((xx, yy, names[0] if len(names) == 1 else (f"{base[0]} {min(n[-1] for n in names)}–{max(n[-1] for n in names)}" if len(base) == 1 else " / ".join(names))))
SL = pd.DataFrame(SL, columns=["x", "y", "lab"])
with rasterio.open(os.path.join(D, "raw", "dem", "uganda_dem_220m_utm36n.tif")) as r:
    DEMSRC = r.read(1).astype("float32"); DEMSRC[(DEMSRC < 300) | (DEMSRC > 5000)] = np.nan; DT, DC = r.transform, r.crs

SPECS = {"swl": dict(title="Static water level", unit="m below ground", cmap="YlGnBu", bounds=[0, 5, 10, 15, 20, 30, 40, 60]),
         "bedrock": dict(title="Depth to bedrock", unit="m below ground", cmap="YlOrBr", bounds=[0, 10, 20, 30, 40, 50, 70]),
         "yield": dict(title="Borehole yield", unit="m³/h", cmap="viridis", bounds=[0, 0.5, 1, 1.5, 2, 3, 5, 10]),
         "wt": dict(title="Water-table elevation", unit="m a.s.l.", cmap="Spectral_r", bounds=list(range(600, 1701, 100)))}

class Area:
    def __init__(self, name, districts, g, maxd, pad):
        self.name = name
        self.poly = adm0.union_all() if districts is None else adm2[adm2.adm2_name.isin(districts)].union_all()
        self.dist = adm2 if districts is None else adm2[adm2.adm2_name.isin(districts)]
        x0, y0, x1, y1 = self.poly.bounds; self.ext = (x0 - pad, x1 + pad, y0 - pad, y1 + pad)
        self.gx, self.gy = np.meshgrid(np.arange(x0, x1, g), np.arange(y0, y1, g)); self.inside = contains(self.poly, self.gx, self.gy)
        self.g, self.maxd = g, maxd
        dem = np.full(self.gx.shape, np.nan, "float32")
        reproject(DEMSRC, dem, src_transform=DT, src_crs=DC, dst_transform=from_origin(x0 - g / 2, y0 - g / 2, g, -g), dst_crs=DC,
                  src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.average)
        self.dem = dem
        x0e, x1e, y0e, y1e = self.ext
        self.b = b[b.x.between(x0e - maxd, x1e + maxd) & b.y.between(y0e - maxd, y1e + maxd)]      # wells incl. neighbours
        self.bin = b[contains(self.poly, b.x.values, b.y.values)]                                     # wells inside
        self.Z = {k: self.surface("swl" if k == "wt" else k, log=(k == "yield")) for k in SPECS}
        self.Z["wt"] = (self.dem - self.Z["wt"][0], self.Z["wt"][1])

    def surface(self, col, log=False):
        v = self.b[col].values.astype(float); ok = np.isfinite(v); px, py, pv = self.b.x.values[ok], self.b.y.values[ok], v[ok]
        if log: pv = np.log10(np.clip(pv, 0.1, 50))
        d, i = cKDTree(np.c_[px, py]).query(np.c_[self.gx.ravel(), self.gy.ravel()], k=15)
        w = 1 / np.maximum(d, 200) ** 1.5; z = ((w * pv[i]).sum(1) / w.sum(1)).reshape(self.gx.shape)
        m = (d[:, 0] <= self.maxd).reshape(self.gx.shape) & self.inside
        sg = 1.0 if self.g <= 1500 else 0.8
        a = gaussian_filter(np.where(m, z, 0), sg); wt = gaussian_filter(m.astype(float), sg)
        z = np.where(m, a / np.maximum(wt, 1e-6), np.nan)
        return (10 ** z if log else z), int(np.isfinite(self.bin[col]).sum())

def scalebar(ax, km):
    xa, xb = ax.get_xlim(); ya, yb = ax.get_ylim(); X = xa + 0.66 * (xb - xa); Y = ya + 0.035 * (yb - ya); L = km * 1000; hgt = 0.008 * (yb - ya)
    ax.add_patch(Rectangle((xa + 0.64 * (xb - xa), ya + 0.01 * (yb - ya)), L + 0.09 * (xb - xa), 0.07 * (yb - ya), fc="w", ec="none", alpha=0.85, zorder=6))
    for j in range(2): ax.add_patch(Rectangle((X + j * L / 2, Y), L / 2, hgt, fc="k" if j == 0 else "w", ec="k", lw=0.6, zorder=7))
    for j, t in enumerate([0, km // 2, km]): ax.text(X + j * L / 2, Y + 2.2 * hgt, f"{t}", ha="center", fontsize=6.5, zorder=7)
    ax.text(X + L + 0.012 * (xb - xa), Y, "km", fontsize=6.5, zorder=7)
    ax.annotate("", xy=(xb - 0.05 * (xb - xa), yb - 0.04 * (yb - ya)), xytext=(xb - 0.05 * (xb - xa), yb - 0.11 * (yb - ya)),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4), zorder=7)
    ax.text(xb - 0.05 * (xb - xa), yb - 0.135 * (yb - ya), "N", ha="center", va="center", fontsize=9, fontweight="bold", zorder=7)

def draw(ax, A_, k, labels="districts", sb=50):
    s = SPECS[k]; z, n = A_.Z[k]; bd = s["bounds"]
    if k == "wt":   # classes adapted to the area
        lo, hi = np.nanpercentile(z, 3), np.nanpercentile(z, 97)
        step = next(st_ for st_ in (10, 20, 25, 50, 100, 200) if (hi - lo) / st_ <= 14)
        bd = list(np.arange(np.floor(lo / step) * step, np.ceil(hi / step) * step + step, step))
    ext = "both" if k == "wt" else "max"; cmap = plt.get_cmap(s["cmap"], len(bd) + (1 if ext == "both" else 0)); norm = BoundaryNorm(bd, cmap.N, extend=ext)
    A_.dist.plot(ax=ax, color="#f2f2f2", zorder=0.1)
    im = ax.pcolormesh(A_.gx, A_.gy, z, cmap=cmap, norm=norm, shading="auto", zorder=0.3, rasterized=True)
    if k == "wt":
        zz = gaussian_filter(np.nan_to_num(z, nan=np.nanmean(z)), 1.5) * np.where(np.isfinite(z), 1, np.nan)
        cs = ax.contour(A_.gx, A_.gy, zz, levels=bd[::2], colors="k", linewidths=0.45, zorder=3); ax.clabel(cs, fontsize=5.5, fmt="%d")
    lakes.plot(ax=ax, color="#cfe3f3", edgecolor="#8fb4d6", lw=0.3, zorder=2)
    adm2.boundary.plot(ax=ax, color="#b5b5b5", lw=0.3, zorder=2.5); A_.dist.boundary.plot(ax=ax, color="#5e5e5e", lw=0.5, zorder=2.6)
    adm0.boundary.plot(ax=ax, color="k", lw=0.9, zorder=2.7)
    ax.scatter(A_.bin.x, A_.bin.y, s=1.2 if labels != "regions" else 0.4, c="k", alpha=0.5, lw=0, zorder=3.5)
    x0, x1, y0, y1 = A_.ext; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    vis = SL[SL.x.between(x0, x1) & SL.y.between(y0, y1)]
    ax.scatter(vis.x, vis.y, s=26, marker="^", c="#d7301f", edgecolor="w", lw=0.6, zorder=5)
    T = []
    if labels == "districts":
        for _, r in A_.dist.iterrows():
            c = r.geometry.representative_point(); T.append(ax.text(c.x, c.y, r.adm2_name, fontsize=6.3, ha="center", va="center", style="italic", color="#1f1f1f", zorder=4.5))
        for _, r in vis.iterrows(): T.append(ax.text(r.x, r.y, r.lab, fontsize=6.3, color="#a3150a", fontweight="bold", zorder=5.5))
    elif labels == "regions":
        for _, r in adm1.iterrows():
            c = r.geometry.representative_point(); T.append(ax.text(c.x, c.y, r.adm1_name.upper(), fontsize=7, color="#555", ha="center", fontweight="bold", zorder=4.5))
    if T: adjust_text(T, ax=ax, expand=(1.12, 1.25), arrowprops=dict(arrowstyle="-", color="#a3150a", lw=0.35), only_move={"text": "xy"})
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    for sp in ax.spines.values(): sp.set_color("#888")
    if labels == "subregions":
        for nm, ds in SUB.items():
            g_ = adm2[adm2.adm2_name.isin(ds)].union_all()
            gpd.GeoSeries([g_], crs=U).boundary.plot(ax=ax, color="#b2182b", lw=1.1, ls="--", zorder=4)
            c = g_.representative_point()
            T.append(ax.text(c.x, c.y, nm, fontsize=7.5, fontweight="bold", color="#b2182b", ha="center", zorder=6,
                             bbox=dict(fc="w", ec="#b2182b", lw=0.5, alpha=0.85, pad=1.5)))
        adjust_text(T, ax=ax, only_move={"text": "xy"})
    cb = plt.colorbar(im, ax=ax, shrink=0.62, pad=0.01, ticks=bd if len(bd) < 10 else bd[::2])
    cb.set_label(f"{s['title']}, {s['unit']}", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    scalebar(ax, sb)
    ax.set_xlabel(""); ax.set_ylabel("")
    return n

def page(A_, fname, title, labels, sb):
    x0, x1, y0, y1 = A_.ext; asp = (x1 - x0) / (y1 - y0)
    if asp >= 1: w_in = 4.7; h_in = w_in / asp
    else: h_in = 5.2; w_in = h_in * asp
    fig, axs = plt.subplots(2, 2, figsize=(2 * (w_in + 1.25), 2 * (h_in + 0.45) + (0.95 if "\n" in title else 0.6)))
    for ax, k, lab in zip(axs.ravel(), ["swl", "bedrock", "yield", "wt"], "abcd"):
        n = draw(ax, A_, k, labels, sb); ax.set_title(f"({lab}) {SPECS[k]['title']}, {SPECS[k]['unit']} (n = {n:,})", fontsize=9, loc="left", pad=4)
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.997, linespacing=1.5)
    fig.text(0.01, 0.004, "Inverse-distance interpolation (k = 15, power 1.5) with light smoothing; blank where no borehole lies within "
             f"{A_.maxd / 1000:.0f} km. Dots: boreholes; red triangles: monitoring stations. Water table = DEM − interpolated static level. "
             "Boundaries: OCHA/UBOS 2026; UTM 36N.", fontsize=6.5, color="#444")
    plt.tight_layout(rect=(0, 0.015, 1, 0.975), h_pad=2.4, w_pad=0.6)
    plt.savefig(os.path.join(F, fname), dpi=200); plt.close()

rows = []
def summarise(name, A_):
    d = A_.bin
    rows.append(dict(area=name, boreholes=len(d), median_swl_m=d.swl.median(), median_bedrock_m=d.bedrock.median(), median_depth_m=d.depth.median(),
                     median_yield_m3h=d["yield"].median(), share_yield_lt1=(d["yield"] < 1).mean(), share_swl_gt30=(d.swl > 30).mean(),
                     stations=int(contains(A_.poly, st.geometry.x.values, st.geometry.y.values).sum())))

# corridor (chapter figure)
Ac = Area("Cattle corridor", CC, 1500, 12000, 10000); summarise("Cattle corridor", Ac)
page(Ac, "fig03_corridor_mapview.png", "Uganda's cattle corridor in map view (44 districts; dashed: sub-regions)", "subregions", 100)
# Uganda
Au = Area("Uganda", None, 2000, 12000, 10000); summarise("Uganda", Au)
page(Au, "figB1_uganda_mapview.png", f"Uganda: groundwater from {len(Au.bin):,} located boreholes", "regions", 100)
ZOOM = []
for j, (nm, ds) in enumerate(SUB.items(), start=2):
    Ar = Area(nm, ds, 1000, 10000, 8000); summarise(nm, Ar); ZOOM.append((nm, Ar))
    sb = 50
    page(Ar, f"figB{j}_{nm.split()[0].lower().replace('-', '')}_mapview.png", f"Cattle corridor – {nm}\n" + "\u2009" + ", ".join(ds), "districts", sb)
    print("done", nm, len(Ar.bin))
ZOOM.append(None)
pd.DataFrame(rows).round(3).to_csv(os.path.join(D, "processed", "region_summary.csv"), index=False)
print(pd.DataFrame(rows).round(2).to_string())

# ---------- chapter Fig. 4: zoom into the four sub-regions (static water level and yield) ----------
from matplotlib.gridspec import GridSpec
Zs = [z for z in ZOOM if z]
fig = plt.figure(figsize=(11, 15.5))
gs = GridSpec(4, 2, figure=fig, height_ratios=[(z[1].ext[3] - z[1].ext[2]) / (z[1].ext[1] - z[1].ext[0]) for z in Zs])
for r_, (nm, Ar) in enumerate(Zs):
    for c_, k in enumerate(["swl", "yield"]):
        ax = fig.add_subplot(gs[r_, c_]); n = draw(ax, Ar, k, "districts", 50)
        ax.set_title(f"({'abcdefgh'[2 * r_ + c_]}) {nm}: {SPECS[k]['title'].lower()}, {SPECS[k]['unit']} (n = {n:,})", fontsize=8.5, loc="left")
fig.suptitle("Zoom into the cattle corridor sub-regions: static water level (left) and borehole yield (right)", fontsize=11, fontweight="bold", y=0.997)
plt.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.2)
plt.savefig(os.path.join(F, "fig04_subregions_zoom.png"), dpi=190); plt.close()
print("zoom figure done")

# ---------- single full-page maps: one per area and variable ----------
FM = os.path.join(H, "..", "figures", "maps")
AREAS = [("uganda", "Uganda", Au, "regions", 100), ("cattle_corridor", "Cattle corridor", Ac, "subregions", 100)] + \
        [(nm.split()[0].lower().replace("-", ""), f"Cattle corridor – {nm}", Ar, "districts", 50) for nm, Ar in Zs]
VARS = {"swl": "static_water_level", "bedrock": "depth_to_bedrock", "yield": "borehole_yield", "wt": "water_table_elevation"}
for slug, title, A_, lab, sb in AREAS:
    os.makedirs(os.path.join(FM, slug), exist_ok=True)
    x0, x1, y0, y1 = A_.ext; asp = (x1 - x0) / (y1 - y0)
    w_in = 8.6; h_in = min(10.2, w_in / asp * 0.9 + 0.8)
    if asp < 0.8: h_in = 10.2; w_in = min(8.6, (h_in - 1.3) * asp / 0.9 + 1.6)
    for k, vn in VARS.items():
        fig, ax = plt.subplots(figsize=(w_in, h_in))
        n = draw(ax, A_, k, lab, sb)
        ax.set_title(f"{title}: {SPECS[k]['title'].lower()} ({SPECS[k]['unit']}), from {n:,} boreholes", fontsize=10.5, loc="left")
        import textwrap
        note = ("Inverse-distance interpolation (k = 15, power 1.5), light smoothing; blank where no borehole lies within "
                f"{A_.maxd / 1000:.0f} km. Dots: boreholes; red triangles: monitoring stations." + (" Water table = DEM − interpolated static level." if k == "wt" else "")
                + " Boundaries: OCHA/UBOS 2026; UTM 36N.")
        plt.tight_layout()
        ax.text(0, -0.015, textwrap.fill(note, int(w_in * 19)), transform=ax.transAxes, fontsize=6.3, color="#444", va="top")
        plt.savefig(os.path.join(FM, slug, f"{slug}_{vn}.png"), dpi=220, bbox_inches="tight", pad_inches=0.12); plt.close()
    print("single maps", slug)
