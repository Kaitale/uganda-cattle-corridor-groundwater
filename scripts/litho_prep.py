"""Classify Aswa drilling logs into cover / weathered / fresh rock and pick layer depths.

Input : data/raw/boreholes/Litho_strati.xlsx (sheet 'Aswa cat 1933-2012', MWE)
        data/raw/dem/aswa_dem_110m_utm36n.tif (optional, ground elevation)
Output: data/processed/aswa_layer_picks.csv, data/processed/aswa_layer_picks_dem.csv
"""
import os, re, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data")


def classify(text):
    """Map a driller's lithology description to 'cover', 'weathered' or 'fresh' (None if unknown)."""
    s = str(text).lower().strip()
    if s in ("", "nan"):
        return None
    if re.search(r"fresh|hard rock|grey rock|massive|gneiss$|^granite$|^rock$|basement|quartzite", s):
        return "fresh"
    if re.search(r"weather|fractur|decompos|saprol|sand.*rock|jointed", s):
        return "weathered"
    if re.search(r"soil|clay|laterit|murram|sand|silt|gravel|loam", s):
        return "cover"
    if re.search(r"granite|gneiss|rock", s):
        return "fresh"
    return None


def picks(row, n=15):
    """Return depth to base of cover, top of fresh rock and weathered thickness for one log."""
    prev, seq = 0.0, []
    for i in range(1, n + 1):
        to = pd.to_numeric(row.get(f"Lithology{i}To"), errors="coerce"); c = classify(row.get(f"Lithology{i}Describe"))
        if pd.isna(to) or c is None:
            continue
        seq.append((prev, to, c)); prev = to
    if not seq:
        return None
    fresh = min([a for a, b, c in seq if c == "fresh"], default=np.nan)
    return dict(depth=prev, cover_base=max([b for a, b, c in seq if c == "cover"], default=0),
                fresh_top=fresh, weathered_thk=sum(b - a for a, b, c in seq if c == "weathered"), regolith=fresh)


def main():
    t = pd.read_excel(os.path.join(D, "raw", "boreholes", "Litho_strati.xlsx"), "Aswa cat 1933-2012")
    rows = []
    for _, r in t.iterrows():
        p = picks(r)
        if p:
            rows.append(dict(id=r["SourceID"], district=r["District"], x=pd.to_numeric(r["UTME"], errors="coerce"),
                             y=pd.to_numeric(r["UTMNS"], errors="coerce"), **p))
    d = pd.DataFrame(rows).dropna(subset=["x", "y"])
    d = d[d.x.between(300000, 800000) & d.y.between(200000, 450000)]          # UTM 36N window, Aswa
    os.makedirs(os.path.join(D, "processed"), exist_ok=True)
    d.to_csv(os.path.join(D, "processed", "aswa_layer_picks.csv"), index=False)
    print(f"logs read {len(t)}, with lithology {len(rows)}, located {len(d)}")
    print(d[["depth", "cover_base", "fresh_top", "weathered_thk"]].describe().round(1).to_string())
    dem = os.path.join(D, "raw", "dem", "aswa_dem_110m_utm36n.tif")
    if os.path.exists(dem):
        import rasterio
        with rasterio.open(dem) as src:
            z = np.array([v[0] for v in src.sample(zip(d.x, d.y))], dtype=float)
            z[(z == src.nodata) | (z < 0)] = np.nan
        d["ground_elev_m"] = z
        d.dropna(subset=["ground_elev_m"]).to_csv(os.path.join(D, "processed", "aswa_layer_picks_dem.csv"), index=False)
        print("with DEM elevation", int(d.ground_elev_m.notna().sum()))


if __name__ == "__main__":
    main()
