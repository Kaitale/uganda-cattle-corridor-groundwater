"""Run the full workflow in order.

    python scripts/run_all.py              # analysis + figures (uses ERA5 already extracted)
    python scripts/run_all.py --era5       # also download ERA5-Land and extract at stations
    python scripts/run_all.py --no-gempy   # skip the GemPy model (slow, needs gempy)
"""
import argparse, os, subprocess, sys, time
H = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser(); p.add_argument("--era5", action="store_true"); p.add_argument("--no-gempy", action="store_true")
a = p.parse_args()
steps = (["get_era5_pev.py", "extract_era5_to_stations.py"] if a.era5 else []) + [
    "borehole_stats.py",        # national borehole statistics
    "borehole_locate.py",       # locate national-database boreholes (UTM zone/hemisphere checks)
    "litho_prep.py",            # Aswa lithology picks (+ DEM elevation)
    "pastas_trend_drought.py",  # Pastas + trend, SGI, drought thresholds
] + ([] if a.no_gempy else ["gempy_kitgum.py"]) + [
    "cattle_corridor.py",       # Fig. 1, 4: corridor map, 3D wells (+ interactive HTML)
    "regional_maps.py",         # Fig. 3 and atlas B1-B5: Uganda and corridor sub-regions in map view
    "corridor_maps.py",         # C0-C4: full-page corridor maps per variable
    "kai_3d_showcase.py",       # interactive 3D (KAI showcase page)
    "make_framework.py",        # Fig. 2: framework diagram
    "make_maps.py",             # Fig. 5, 6 (district and Aswa maps)
    "make_paper_figures.py",    # Fig. 7, 8, 9, A1
]
for s in steps:
    t = time.time(); print(f"\n=== {s}", flush=True)
    r = subprocess.run([sys.executable, os.path.join(H, s)], cwd=H)
    if r.returncode:
        sys.exit(f"{s} failed (exit {r.returncode})")
    print(f"--- {s} done in {time.time() - t:.0f} s")
print("\nAll steps finished. Tables: data/processed, figures: figures/paper")
