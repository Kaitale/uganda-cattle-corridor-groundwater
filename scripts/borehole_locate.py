"""Locate national-database boreholes: test UTM 36N / 36S / 35N readings of UTME/UTMNS and keep the
reading that falls inside the borehole's recorded district (2026 boundaries).

Output: data/processed/boreholes_located.csv (lon, lat, x/y in EPSG:32636, ground elevation, key fields)
"""
import os, re, warnings, numpy as np, pandas as pd, geopandas as gpd, rasterio
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data")
norm = lambda s: re.sub(r"[^a-z]", "", str(s).lower().replace("district", "").replace("city", ""))
ALIAS = {"ssembabule": "sembabule", "sembabule": "sembabule", "kasanda": "kassanda", "luweero": "luwero"}


def main():
    a = pd.read_excel(os.path.join(D, "raw", "boreholes", "BOREHOLE DATA.xlsx"))
    a["E"] = pd.to_numeric(a.UTME, errors="coerce"); a["N"] = pd.to_numeric(a.UTMNS, errors="coerce")
    a["key"] = a.District.map(norm).replace(ALIAS)
    adm2 = gpd.read_file(os.path.join(D, "raw", "admin", "uga_admin2.shp")); adm2["key"] = adm2.adm2_name.map(norm).replace(ALIAS)
    a = a[a.E.between(1e5, 9e5) & a.N.between(1, 1e7)].copy()
    best = pd.Series(index=a.index, dtype=object)
    # readings: 36N, 36S (10,000 km false northing), 36N with the southern northing written without sign, 35N
    for name, epsg, sign in (("36N", 32636, 1), ("36S", 32736, 1), ("36N-south", 32636, -1), ("35N", 32635, 1)):
        g = gpd.GeoDataFrame(a[["key"]], geometry=gpd.points_from_xy(a.E, sign * a.N), crs=epsg).to_crs(4326)
        j = gpd.sjoin(g, adm2[["key", "adm2_name", "adm1_name", "geometry"]], predicate="within", how="left")
        j = j[~j.index.duplicated()]
        ok = (j.key_left == j.key_right) & best.isna()
        best[ok] = name
        a.loc[ok, "lon"] = j.geometry.x[ok]; a.loc[ok, "lat"] = j.geometry.y[ok]; a.loc[ok, "adm2_name"] = j.adm2_name[ok]; a.loc[ok, "adm1_name"] = j.adm1_name[ok]
    a["crs_read"] = best; a = a.dropna(subset=["lon"])
    u = gpd.GeoSeries(gpd.points_from_xy(a.lon, a.lat), crs=4326).to_crs(32636)
    a["x"], a["y"] = u.x.values, u.y.values
    with rasterio.open(os.path.join(D, "raw", "dem", "uganda_dem_220m_utm36n.tif")) as r:
        z = np.array([v[0] for v in r.sample(zip(a.x, a.y))], float)
    z[(z < 300) | (z > 5000)] = np.nan; a["ground_elev_m"] = z
    for c, n, lo, hi in [("TotalDrilledDepth", "depth", 5, 300), ("DepthToBedRock", "bedrock", 0, 200), ("StaticWaterLevel", "swl", 0, 150),
                         ("ConstantDischargeYield", "yield", 0, 50), ("1st Water Strike", "strike1", 0, 300)]:
        v = pd.to_numeric(a[c], errors="coerce"); a[n] = v.where(v.between(lo, hi))
    out = a[["SourceID", "District", "adm2_name", "adm1_name", "lon", "lat", "x", "y", "crs_read", "ground_elev_m", "depth", "bedrock", "swl", "yield", "strike1"]]
    out.to_csv(os.path.join(D, "processed", "boreholes_located.csv"), index=False)
    print(f"located {len(out):,} of {len(pd.read_excel(os.path.join(D, 'raw', 'boreholes', 'BOREHOLE DATA.xlsx'))):,}; by CRS", out.crs_read.value_counts().to_dict())


if __name__ == "__main__":
    main()
