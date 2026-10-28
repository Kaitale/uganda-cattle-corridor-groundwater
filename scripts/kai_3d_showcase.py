"""KAI Water & Energy showcase: interactive 3D groundwater models (Plotly).

Three areas, switched with buttons:
  1. Central cattle corridor (Nakasongola, Luwero)      national borehole database
  2. South-west cattle corridor (Isingiro ... Kazo)      national borehole database
  3. Aswa catchment, Kitgum and Pader                    drilling logs with lithology layers
Hover over any well section for the well ID, district, layer type and depth interval.
Output: figures/interactive/kai_cattle_corridor_3d.html
"""
import os, sys, numpy as np, pandas as pd, geopandas as gpd, plotly.graph_objects as go
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from shapely import contains_xy as contains
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data"); FI = os.path.join(H, "..", "figures", "interactive")
sys.path.insert(0, H); from litho_prep import classify
os.makedirs(FI, exist_ok=True)
adm2 = gpd.read_file(os.path.join(D, "raw", "admin", "uga_admin2.shp")).to_crs(32636)
b = pd.read_csv(os.path.join(D, "processed", "boreholes_located.csv"))
st = pd.read_csv(os.path.join(D, "raw", "Stations.csv")).dropna(subset=["latitude"])
st = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.longitude, st.latitude), crs=4326).to_crs(32636)
COL = {"Regolith (soil, clay, weathered rock)": "#c58b3c", "Bedrock": "#2b2b2b", "Cover (soil, clay, laterite)": "#d9a65b",
       "Weathered / fractured rock (aquifer)": "#2f8fc4", "Fresh rock": "#3d3d3d", "Unclassified": "#a8a8a8"}
f1 = lambda v, u="m": "n/a" if pd.isna(v) else f"{v:.1f} {u}"

def idw(px, py, pv, gx, gy, maxd=12000):
    d, i = cKDTree(np.c_[px, py]).query(np.c_[gx.ravel(), gy.ravel()], k=15); w = 1 / np.maximum(d, 200) ** 1.5
    z = ((w * pv[i]).sum(1) / w.sum(1)).reshape(gx.shape); m = (d[:, 0] <= maxd).reshape(gx.shape)
    a = gaussian_filter(np.where(m, z, 0), 1.2); wt = gaussian_filter(m.astype(float), 1.2)
    return np.where(m, a / np.maximum(wt, 1e-6), np.nan)

def well_traces(segs, x0, y0, vis):
    """segs: id, district, x, y, top, bot, layer, desc, info -> one line trace per layer type, hover per vertex."""
    out = []
    for layer, colr in COL.items():
        q = segs[segs.layer == layer]
        if q.empty: continue
        xs, ys, zs, tx = [], [], [], []
        for r in q.itertuples():
            t = f"<b>{r.id}</b> · {r.district}<br><b>{layer}</b>: {r.top:.0f}–{r.bot:.0f} m" + (f"<br><i>{r.desc}</i>" if r.desc else "") + f"<br>{r.info}"
            xs += [(r.x - x0) / 1000] * 2 + [None]; ys += [(r.y - y0) / 1000] * 2 + [None]; zs += [-r.top, -r.bot, None]; tx += [t, t, None]
        out.append(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=colr, width=5), name=layer, visible=vis, text=tx,
                                hovertemplate="%{text}<extra></extra>"))
    return out

def surf(X, Y, x0, y0, Z, cs, lo, hi, name, title, ypos, vis):
    return go.Surface(x=(X - x0) / 1000, y=(Y - y0) / 1000, z=-Z, surfacecolor=Z, colorscale=cs, cmin=lo, cmax=hi, opacity=0.75, name=name, showlegend=True,
                      colorbar=dict(title=title, x=1.0, len=0.42, y=ypos, thickness=14), visible=vis, hovertemplate=name + " at %{surfacecolor:.1f} m depth<extra></extra>")

def stations(reg, x0, y0, vis):
    q = st[contains(reg, st.geometry.x.values, st.geometry.y.values)]
    return go.Scatter3d(x=(q.geometry.x - x0) / 1000, y=(q.geometry.y - y0) / 1000, z=[12] * len(q), mode="markers+text", marker=dict(size=6, color="#d7301f", symbol="diamond"),
                        text=q.station_name, textposition="top center", textfont=dict(color="#a3150a", size=11), name="Monitoring station", visible=vis,
                        hovertemplate="Monitoring station %{text}<extra></extra>")

traces, groups = [], []
# --- corridor blocks (national database)
for k, (title, ds) in enumerate([("Central cattle corridor: Nakasongola and Luwero", ["Nakasongola", "Luwero"]),
                                 ("South-west cattle corridor: Isingiro, Mbarara, Kiruhura, Kazo", ["Isingiro", "Mbarara", "Kiruhura", "Kazo"])]):
    vis = k == 0; n0 = len(traces)
    reg = adm2[adm2.adm2_name.isin(ds)].union_all(); x0, y0, x1, y1 = reg.bounds
    X, Y = np.meshgrid(np.arange(x0, x1, 2000), np.arange(y0, y1, 2000)); inr = contains(reg, X, Y)
    w = b[b.adm2_name.isin(ds)].dropna(subset=["depth"]); ws = w.dropna(subset=["swl"]); wb = w.dropna(subset=["bedrock"])
    S = idw(ws.x.values, ws.y.values, ws.swl.values, X, Y); S[~inr] = np.nan
    B = idw(wb.x.values, wb.y.values, wb.bedrock.values, X, Y); B[~inr] = np.nan
    traces += [surf(X, Y, x0, y0, S, "Blues", 0, 45, "Water table", "Water table,<br>m below ground", 0.77, vis),
               surf(X, Y, x0, y0, B, [[0, "#ececec"], [1, "#4d4d4d"]], 0, 70, "Top of bedrock", "Top of bedrock,<br>m below ground", 0.27, vis)]
    rows = []
    for r in w.itertuples():
        info = f"water level {f1(r.swl)} · yield {f1(r._14, 'm³/h')} · drilled {f1(r.depth)}"
        if pd.notna(r.bedrock) and r.bedrock < r.depth:
            rows += [(r.SourceID, r.adm2_name, r.x, r.y, 0, r.bedrock, "Regolith (soil, clay, weathered rock)", "", info),
                     (r.SourceID, r.adm2_name, r.x, r.y, r.bedrock, r.depth, "Bedrock", "", info)]
        else:
            rows.append((r.SourceID, r.adm2_name, r.x, r.y, 0, r.depth, "Unclassified", "", info))
    traces += well_traces(pd.DataFrame(rows, columns=["id", "district", "x", "y", "top", "bot", "layer", "desc", "info"]), x0, y0, vis)
    traces.append(go.Scatter3d(x=(ws.x - x0) / 1000, y=(ws.y - y0) / 1000, z=-ws.swl, mode="markers", marker=dict(size=2.8, color="#0b4ea2"), name="Water level in well",
                               visible=vis, text=[f"<b>{a}</b> · {c}<br>static water level {s_:.1f} m" for a, c, s_ in zip(ws.SourceID, ws.adm2_name, ws.swl)],
                               hovertemplate="%{text}<extra></extra>"))
    traces.append(stations(reg, x0, y0, vis))
    groups.append((title, n0, len(traces), f"{len(w)} boreholes from the national database", (x1 - x0) / 1000, (y1 - y0) / 1000))
# --- Aswa: drilling logs with lithology
LAB = {"cover": "Cover (soil, clay, laterite)", "weathered": "Weathered / fractured rock (aquifer)", "fresh": "Fresh rock"}
t = pd.read_excel(os.path.join(D, "raw", "boreholes", "Litho_strati.xlsx"), "Aswa cat 1933-2012")
t["x"] = pd.to_numeric(t.UTME, errors="coerce"); t["y"] = pd.to_numeric(t.UTMNS, errors="coerce")
ds = ["Kitgum", "Pader"]; reg = adm2[adm2.adm2_name.isin(ds)].union_all(); x0, y0, x1, y1 = reg.bounds
t = t.dropna(subset=["x", "y"]); t = t[contains(reg, t.x.values, t.y.values)]
rows = []
for r in t.itertuples(index=False):
    rr = r._asdict(); prev = 0.0; seq = []
    for i in range(1, 16):
        to = pd.to_numeric(rr.get(f"Lithology{i}To"), errors="coerce"); desc = rr.get(f"Lithology{i}Describe")
        if pd.isna(to) or to <= prev: continue
        c = classify(desc); seq.append((prev, float(to), LAB.get(c, "Unclassified"), str(desc).strip().capitalize() if pd.notna(desc) else ""))
        prev = float(to)
    for a, bb, lay, de in seq:
        rows.append((rr["SourceID"], str(rr["District"]).title(), rr["x"], rr["y"], a, bb, lay, de, f"drilled {prev:.0f} m"))
segs = pd.DataFrame(rows, columns=["id", "district", "x", "y", "top", "bot", "layer", "desc", "info"])
pk = pd.read_csv(os.path.join(D, "processed", "aswa_layer_picks.csv")); pk = pk[contains(reg, pk.x.values, pk.y.values)]
pk = pk[(pk.fresh_top >= pk.cover_base) & (pk.fresh_top < 100)]
X, Y = np.meshgrid(np.arange(x0, x1, 2000), np.arange(y0, y1, 2000)); inr = contains(reg, X, Y)
CB = idw(pk.x.values, pk.y.values, pk.cover_base.values, X, Y, 10000); CB[~inr] = np.nan
FT = idw(pk.dropna(subset=["fresh_top"]).x.values, pk.dropna(subset=["fresh_top"]).y.values, pk.dropna(subset=["fresh_top"]).fresh_top.values, X, Y, 10000); FT[~inr] = np.nan
vis = False; n0 = len(traces)
traces += [surf(X, Y, x0, y0, CB, [[0, "#f6e3bf"], [1, "#b07a2a"]], 0, 30, "Base of cover", "Base of cover,<br>m below ground", 0.77, vis),
           surf(X, Y, x0, y0, FT, [[0, "#ececec"], [1, "#4d4d4d"]], 0, 60, "Top of fresh rock", "Top of fresh rock,<br>m below ground", 0.27, vis)]
traces += well_traces(segs, x0, y0, vis); traces.append(stations(reg, x0, y0, vis))
groups.append(("Aswa catchment: Kitgum and Pader", n0, len(traces), f"{segs.id.nunique()} drilling logs with {len(segs):,} logged layers", (x1 - x0) / 1000, (y1 - y0) / 1000))

buttons = []
for title, a0, a1, sub, wx, wy in groups:
    v = [a0 <= i < a1 for i in range(len(traces))]
    buttons.append(dict(label=title.split(":")[0], method="update",
                        args=[{"visible": v}, {"title.text": f"{title} · {sub}", "scene.aspectratio": dict(x=1, y=wy / wx, z=0.45)}]))
t0, _, _, sub0, wx0, wy0 = groups[0]
fig = go.Figure(traces)
fig.update_layout(title=dict(text=f"{t0} · {sub0}", x=0.02, y=0.97, font=dict(size=15)),
                  updatemenus=[dict(type="buttons", direction="right", x=0.02, y=1.0, xanchor="left", yanchor="bottom", buttons=buttons, pad=dict(t=4))],
                  scene=dict(aspectmode="manual", aspectratio=dict(x=1, y=wy0 / wx0, z=0.45), xaxis_title="km east", yaxis_title="km north",
                             zaxis=dict(title="depth below ground, m", range=[-130, 20]), camera=dict(eye=dict(x=-1.3, y=-1.5, z=0.9))),
                  legend=dict(x=0.01, y=0.88, bgcolor="rgba(255,255,255,0.75)", font=dict(size=11)), hoverlabel=dict(bgcolor="white", font_size=12),
                  margin=dict(l=0, r=0, t=80, b=0), height=740)
html = fig.to_html(include_plotlyjs="cdn", full_html=False, config=dict(displaylogo=False, responsive=True))
page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Groundwater in 3D</title>
<style>:root{{--ink:#12263a;--muted:#5b6b7b;--bg:#f7f9fb;--card:#fff;--head:#0e2a47}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{--ink:#e7eef6;--muted:#a9b8c7;--bg:#0f1720;--card:#ffffff;--head:#0b1f35}}}}
body{{margin:0;font-family:Inter,Segoe UI,Arial,sans-serif;background:var(--bg);color:var(--ink)}}
header{{padding:18px 24px;background:var(--head);color:#fff}}header h1{{margin:0;font-size:22px}}header p{{margin:4px 0 0;color:#bcd3ea;font-size:14px}}
main{{max-width:1240px;margin:0 auto;padding:16px}}.card{{background:var(--card);border-radius:10px;box-shadow:0 1px 4px rgba(0,0,0,.08);padding:8px;overflow:hidden}}
.note{{color:var(--muted);font-size:13px;line-height:1.55;margin:12px 4px}}footer{{text-align:center;color:var(--muted);font-size:12px;padding:16px}}</style></head>
<body><header><h1>KAI Water &amp; Energy · Groundwater in 3D</h1><p>Water, foreseen. Boreholes, layers, water table and bedrock in Uganda's drylands, built from national records with open tools.</p></header>
<main><div class="card">{html}</div>
<p class="note">Drag to rotate, scroll to zoom, and hover over a well for its ID, district, layer type and depth interval (in the Aswa logs, also the driller's description). Click a legend entry to hide or show a layer; use the buttons to switch area. Surfaces are interpolated between boreholes (inverse distance, 2 km grid) and shown as depth below ground; red diamonds are national monitoring stations.
Data: national borehole database, drilling logs and monitoring network, Ministry of Water and Environment, Uganda; boundaries OCHA/UBOS. Method and code: open workflow (Pastas, GemPy, GeoPandas, Plotly).</p></main>
<footer>KAI Water &amp; Energy, Kampala · hydrological forecasting · groundwater and flood digital twins · water-utility software · irrigation · subsurface modelling</footer></body></html>"""
open(os.path.join(FI, "kai_cattle_corridor_3d.html"), "w").write(page)
print("written", round(len(page) / 1e6, 2), "MB;", [(g[0], g[3]) for g in groups])
