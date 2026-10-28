"""Realistic 3D groundwater model of three regions of Uganda in three.js (GroundWater Analyst), one tab per region.

Regions: central cattle corridor (Nakasongola, Luwero), south-west cattle corridor (Isingiro, Mbarara, Kiruhura,
Kazo) and the Aswa catchment (Kitgum, Pader, drilling logs with lithology). No borehole IDs are written.
Fills scripts/templates/corridor_3d_realistic_template.html and writes figures/interactive/corridor_3d_realistic.html
(open it with ?region=central|southwest|aswa, or use the tabs). The KAI website splits the regions into one JSON each.
"""
import json, os, re, sys
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, uniform_filter
from shapely import contains_xy
from shapely.geometry import box

U = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..') + os.sep
adm2 = gpd.read_file(U + 'data/raw/admin/uga_admin2.shp').to_crs(32636)
lakes = gpd.read_file(U + 'data/raw/admin/Lakes.shp').to_crs(32636)
b = pd.read_csv(U + 'data/processed/boreholes_located.csv')
st = pd.read_csv(U + 'data/raw/Stations.csv').dropna(subset=['latitude'])
stg = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st.longitude, st.latitude), crs=4326).to_crs(32636)
with rasterio.open(U + 'data/raw/dem/uganda_dem_220m_utm36n.tif') as r:
    DEMA = r.read(1).astype(float); DT = r.transform
    if r.nodata is not None: DEMA[DEMA == r.nodata] = np.nan
DEMA = uniform_filter(np.where(np.isnan(DEMA), np.nanmedian(DEMA), DEMA), 3)
# named lakes: name and approximate centre (WGS84 -> UTM 36N, checked against the lake polygons)
LAKE_NAMES = [('Lake Kyoga', 472191, 165797), ('Lake Mburo', 271851, -70044), ('Lake Nakivale', 269999, -86635),
              ('Lake Kachera', 288549, -62672), ('Lake Kijanebalola', 307109, -68195)]


def classify(text):
    """Driller's description -> cover / weathered / fresh (same rules as scripts/litho_prep.py)."""
    s = str(text).lower().strip()
    if s in ('', 'nan'): return None
    if re.search(r"fresh|hard rock|grey rock|massive|gneiss$|^granite$|^rock$|basement|quartzite", s): return 'fresh'
    if re.search(r"weather|fractur|decompos|saprol|sand.*rock|jointed", s): return 'weathered'
    if re.search(r"soil|clay|laterit|murram|sand|silt|gravel|loam", s): return 'cover'
    if re.search(r"granite|gneiss|rock", s): return 'fresh'
    return None


def r1(v): return None if v is None or pd.isna(v) else round(float(v), 1)
def med(s): s = pd.Series(s, dtype=float).dropna(); return round(float(s.median()), 1) if len(s) else None


def region(key, ds, step, title, sub, aswa=False):
    reg = adm2[adm2.adm2_name.isin(ds)]; U_ = reg.union_all()
    x0, y0, x1, y1 = reg.total_bounds; pad = 3000; x0 -= pad; y0 -= pad; x1 += pad; y1 += pad
    nx = int((x1 - x0) // step) + 1; ny = int((y1 - y0) // step) + 1
    X, Y = np.meshgrid(x0 + np.arange(nx) * step, y0 + np.arange(ny) * step)
    rows, cols = rasterio.transform.rowcol(DT, X.ravel(), Y.ravel())
    rows = np.clip(rows, 0, DEMA.shape[0] - 1); cols = np.clip(cols, 0, DEMA.shape[1] - 1)
    D = DEMA[rows, cols].reshape(X.shape)
    bx = box(x0, y0, x1, y1)
    # lakes: mask, flatten each to its own level, label the named ones
    lmask = np.zeros(X.shape, bool); lab = []
    for g in lakes[lakes.intersects(bx)].geometry:
        gi = g.intersection(bx); m = contains_xy(gi, X, Y)
        if m.sum() == 0: continue
        D[m] = np.percentile(D[m], 10); lmask |= m
        if m.sum() > 25:
            for n, lx, ly in LAKE_NAMES:
                if gi.distance(gpd.points_from_xy([lx], [ly])[0]) < 4000:
                    p = gi.representative_point() if n != 'Lake Kyoga' else gi.centroid
                    lab.append([round((p.x - x0) / 1000, 2), round((p.y - y0) / 1000, 2), n])
    near = b[(b.x > x0 - 20000) & (b.x < x1 + 20000) & (b.y > y0 - 20000) & (b.y < y1 + 20000)]

    def idw(px, py, pv, k=15, p=1.5):
        d, i = cKDTree(np.c_[px, py]).query(np.c_[X.ravel(), Y.ravel()], k=k); w = 1 / np.maximum(d, 200) ** p
        return gaussian_filter(((w * pv[i]).sum(1) / w.sum(1)).reshape(X.shape), 1.5)

    def samp(G, px, py):
        fx = (px - x0) / step; fy = (py - y0) / step
        i0 = np.clip(np.floor(fx).astype(int), 0, nx - 2); j0 = np.clip(np.floor(fy).astype(int), 0, ny - 2); tx = fx - i0; ty = fy - j0
        return G[j0, i0] * (1 - tx) * (1 - ty) + G[j0, i0 + 1] * tx * (1 - ty) + G[j0 + 1, i0] * (1 - tx) * ty + G[j0 + 1, i0 + 1] * tx * ty

    ws = near.dropna(subset=['swl']); ws = ws[(ws.swl >= 0) & (ws.swl < 120)]
    S = idw(ws.x.values, ws.y.values, ws.swl.values)
    out = dict(key=key, title=title, sub=sub, aswa=aswa)
    if not aswa:
        wb = near.dropna(subset=['bedrock']); wb = wb[(wb.bedrock > 0) & (wb.bedrock < 200)]
        B = idw(wb.x.values, wb.y.values, wb.bedrock.values); CB = None
        inw = contains_xy(U_, b.x.values, b.y.values)
        w = b[inw].dropna(subset=['depth']).copy().rename(columns={'yield': 'yld'}); w = w[(w.depth > 0) & (w.depth < 300)]
        w['gz'] = samp(D, w.x.values, w.y.values)
        wells = []
        for r in w.itertuples():
            br = r1(r.bedrock) if pd.notna(r.bedrock) and 0 < r.bedrock <= r.depth else None
            sw = r1(r.swl) if pd.notna(r.swl) and 0 <= r.swl <= r.depth else None
            segs = [[0, br, 0], [br, r1(r.depth), 1]] if br is not None else [[0, r1(r.depth), 9]]
            wells.append([round((r.x - x0) / 1000, 3), round((r.y - y0) / 1000, 3), round(float(r.gz), 1), r1(r.depth), br, sw, r1(r.yld), r.adm2_name, segs])
        n_swl = sum(q[5] is not None for q in wells); n_br = sum(q[4] is not None for q in wells)
        out['stats'] = [[len(wells), 'boreholes'], [f"{med([q[5] for q in wells])} m", 'median water level'],
                        [f"{med([q[4] for q in wells])} m", 'median depth to bedrock'], [f"{med(w.yld)} m³/h", 'median yield']]
        names = ' and '.join([', '.join(ds[:-1]), ds[-1]]) if len(ds) > 1 else ds[0]
        out['method'] = (f"Water table and top of bedrock are interpolated between boreholes (inverse distance, {int(step)} m grid, smoothed) "
                         f"from {len(ws)} water levels and {len(wb)} bedrock depths in and around {names}, and hung below the national DEM. "
                         "Depths below ground are exaggerated separately from the terrain relief. The weathered-zone and rock textures are schematic. "
                         "Flow particles follow the slope of the interpolated water-table heads: indicative direction only, not velocity. "
                         f"Of {len(wells)} boreholes shown, {n_swl} have a static water level and {n_br} a depth to bedrock.")
        out['legend'] = 'corridor'
    else:
        t = pd.read_excel(U + 'data/raw/boreholes/Litho_strati.xlsx', 'Aswa cat 1933-2012')
        t['x'] = pd.to_numeric(t.UTME, errors='coerce'); t['y'] = pd.to_numeric(t.UTMNS, errors='coerce')
        t = t.dropna(subset=['x', 'y'])
        # layer picks from every log within 20 km (interpolation), logs drawn only inside the two districts
        CLS = {'cover': 2, 'weathered': 3, 'fresh': 4, None: 9}
        picks, wells, nlay = [], [], 0
        bsw = b.dropna(subset=['swl']).copy(); bsw['k'] = bsw.SourceID.astype(str).str.replace(r'\s+', '', regex=True).str.upper()
        swl_of = bsw.groupby('k').swl.first().to_dict(); yld_of = b.assign(k=b.SourceID.astype(str).str.replace(r'\s+', '', regex=True).str.upper()).dropna(subset=['yield']).groupby('k')['yield'].first().to_dict()
        tin = t[(t.x > x0 - 20000) & (t.x < x1 + 20000) & (t.y > y0 - 20000) & (t.y < y1 + 20000)]
        inside = contains_xy(U_, tin.x.values, tin.y.values)
        for (_, r), ins in zip(tin.iterrows(), inside):
            prev, seq = 0.0, []
            for i in range(1, 16):
                to = pd.to_numeric(r.get(f'Lithology{i}To'), errors='coerce'); desc = r.get(f'Lithology{i}Describe')
                if pd.isna(to) or to <= prev: continue
                c = classify(desc); seq.append([round(prev, 1), round(float(to), 1), CLS[c], (str(desc).strip().capitalize()[:48] if pd.notna(desc) else '')]); prev = float(to)
            if not seq: continue
            known = [s for s in seq if s[2] != 9]
            cb = max([s[1] for s in known if s[2] == 2], default=0.0)
            ft = min([s[0] for s in known if s[2] == 4], default=np.nan)
            if known: picks.append((r.x, r.y, cb, ft))
            if ins:
                k = re.sub(r'\s+', '', str(r.SourceID)).upper()
                sw = swl_of.get(k); sw = r1(sw) if sw is not None and 0 <= sw <= prev else None
                wells.append([round((r.x - x0) / 1000, 3), round((r.y - y0) / 1000, 3), 0.0, r1(prev), r1(ft) if pd.notna(ft) else None, sw,
                              r1(yld_of.get(k)), str(r.District).title(), seq]); nlay += len(seq)
        pk = pd.DataFrame(picks, columns=['x', 'y', 'cb', 'ft'])
        pcb = pk[pk.cb < 80]; pft = pk.dropna(subset=['ft']); pft = pft[(pft.ft >= pft.cb) & (pft.ft < 100)]
        CB = idw(pcb.x.values, pcb.y.values, pcb.cb.values); B = idw(pft.x.values, pft.y.values, pft.ft.values)
        B = np.maximum(B, CB + 1.0)
        gz = samp(D, np.array([w[0] * 1000 + x0 for w in wells]), np.array([w[1] * 1000 + y0 for w in wells]))
        for w, g in zip(wells, gz): w[2] = round(float(g), 1)
        ins_swl = ws[contains_xy(U_, ws.x.values, ws.y.values)]
        out['stats'] = [[len(wells), 'drilling logs'], [f"{nlay:,}", 'logged layers'],
                        [f"{med([max([s[1] for s in w[8] if s[2] == 2], default=0) for w in wells])} m", 'median base of cover'],
                        [f"{med([w[4] for w in wells])} m", 'median top of fresh rock']]
        out['method'] = (f"Each drilling log's layers are classed from the driller's description as cover (soil, clay, laterite), weathered or fractured rock "
                         f"(the aquifer) or fresh rock. Base of cover and top of fresh rock are interpolated (inverse distance, {int(step)} m grid, smoothed) "
                         f"from {len(pcb)} and {len(pft)} logs in and around Kitgum and Pader; the water table from {len(ws)} static water levels in the national "
                         f"borehole database (median {med(ins_swl.swl)} m in the two districts). Depths below ground are exaggerated separately from the terrain relief. "
                         "Flow particles follow the slope of the interpolated water-table heads: indicative direction only, not velocity. "
                         f"{len(wells)} logs with {nlay:,} layers are drawn; {sum(w[5] is not None for w in wells)} could be matched to a water level.")
        out['legend'] = 'aswa'
    H = D - S; BR = D - B
    dmask = contains_xy(U_, X, Y)
    stq = stg[(stg.geometry.x > x0) & (stg.geometry.x < x1) & (stg.geometry.y > y0) & (stg.geometry.y < y1)]
    q = lambda G, dec=1: [round(float(v), dec) for v in G.ravel()]
    def lines_of(geom):
        res = []
        for g in ([geom] if geom.geom_type == 'Polygon' else list(geom.geoms)):
            c = np.array(g.exterior.simplify(300).coords); res.append([[round((px - x0) / 1000, 2), round((py - y0) / 1000, 2)] for px, py in c])
        return res
    out.update(nx=nx, ny=ny, step=step / 1000, w=round((nx - 1) * step / 1000, 3), h=round((ny - 1) * step / 1000, 3),
               dem=q(D, 0), wt=q(H), br=q(BR), cb=(q(D - CB) if CB is not None else None),
               lake=''.join('1' if v else '0' for v in lmask.ravel()), dist=''.join('1' if v else '0' for v in dmask.ravel()),
               wells=wells, outlines={n: lines_of(reg[reg.adm2_name == n].union_all()) for n in ds},
               stations=[[round((g.x - x0) / 1000, 3), round((g.y - y0) / 1000, 3), str(n)] for g, n in zip(stq.geometry, stq.station_name)],
               lakes=lab, utm0=[round(x0), round(y0)])
    print(key, 'grid', nx, ny, nx * ny, 'wells', len(wells), 'lakes', [l[2] for l in lab], 'stations', len(out['stations']), out['stats'])
    return out


R = [region('central', ['Nakasongola', 'Luwero'], 600.0, 'Groundwater beneath Nakasongola and Luwero',
            'Central cattle corridor. Borehole records turned into a 3D model: water table, top of bedrock and aquifer layers under real terrain.'),
     region('southwest', ['Isingiro', 'Mbarara', 'Kiruhura', 'Kazo'], 650.0, 'Groundwater beneath Isingiro, Mbarara, Kiruhura and Kazo',
            'South-west cattle corridor. Borehole records turned into a 3D model: water table, top of bedrock and aquifer layers under real terrain.'),
     region('aswa', ['Kitgum', 'Pader'], 800.0, 'Rock layers beneath Kitgum and Pader',
            'Aswa catchment. Drilling logs read layer by layer: cover, the weathered and fractured rock aquifer, and fresh rock.', aswa=True)]
data = dict(order=[r['key'] for r in R], tabs={'central': 'Central corridor', 'southwest': 'South-west corridor', 'aswa': 'Aswa catchment'},
            regions={r['key']: r for r in R})
s = json.dumps(data, separators=(',', ':'), ensure_ascii=False)
tpl = open(os.path.join(U, 'scripts', 'templates', 'corridor_3d_realistic_template.html'), encoding='utf-8').read()
out = os.path.join(U, 'figures', 'interactive', 'corridor_3d_realistic.html')
open(out, 'w', encoding='utf-8').write(tpl.replace('__DATA__', s)); print('written', out, round(len(s) / 1e6, 2), 'MB data')
