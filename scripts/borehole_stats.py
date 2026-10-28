"""Summary statistics of the national borehole database (depths, levels, yield, water quality).

Input : data/raw/boreholes/BOREHOLE DATA.xlsx (MWE national borehole database)
Output: data/processed/borehole_national_stats.csv
"""
import os, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, "..", "data")
FIELDS = {"TotalDrilledDepth": (0, 300), "DepthToBedRock": (0, 200), "StaticWaterLevel": (0, 150),
          "ConstantDischargeYield": (0, 50)}
LIMITS = {"flouride": 1.5, "iron": 0.3}          # WHO guideline values, mg/L


def main():
    a = pd.read_excel(os.path.join(D, "raw", "boreholes", "BOREHOLE DATA.xlsx"))
    print(f"records {len(a):,}, districts {a.District.nunique()}")
    rows = []
    for c, (lo, hi) in FIELDS.items():
        x = pd.to_numeric(a[c], errors="coerce"); x = x[x.between(lo, hi)]
        rows.append(dict(variable=c, n=len(x), p25=x.quantile(.25), median=x.median(), p75=x.quantile(.75)))
    y = pd.to_numeric(a.ConstantDischargeYield, errors="coerce"); y = y[y.between(0, 50)]
    rows.append(dict(variable="share_yield_ge_3m3h", n=len(y), median=(y >= 3).mean()))
    rows.append(dict(variable="share_yield_lt_0.5m3h", n=len(y), median=(y < 0.5).mean()))
    for c, lim in LIMITS.items():
        x = pd.to_numeric(a[c], errors="coerce").dropna(); x = x[x >= 0]
        rows.append(dict(variable=f"share_{c}_gt_{lim}", n=len(x), median=(x > lim).mean()))
    out = pd.DataFrame(rows).round(3)
    out.to_csv(os.path.join(D, "processed", "borehole_national_stats.csv"), index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
