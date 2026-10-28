"""Framework diagram: data -> models -> indicators -> decisions (groundwater digital twin)."""
import os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
H = os.path.dirname(os.path.abspath(__file__)); F = os.path.join(H, "..", "figures", "paper")
cols = [("DATA", "#e8eef6", "#2c5d8f", ["Telemetry and manual\nwater levels (48 stations)", "National borehole database\n(12,529 records)", "Drilling logs\n(5,052, Aswa)", "ERA5-Land rain and\nevaporation (1998–2026)", "DEM, admin boundaries,\nlakes"]),
        ("MODELS", "#eaf3ea", "#2f7a3d", ["Quality control and\nwell location (UTM checks)", "Time-series models\n(Pastas: recharge + trend)", "2D map view and\n3D well interpolation", "3D geological model\n(GemPy)", "Trend tests\n(Theil–Sen, Mann–Kendall)"]),
        ("INDICATORS", "#fbf0df", "#a46a12", ["Standardised Groundwater\nIndex (SGI)", "Monthly drought thresholds\n(P50, P80, P90)", "Level trend, m/yr", "Aquifer geometry: regolith,\nwater table, bedrock", "Response time and\nrecharge sensitivity"]),
        ("DECISIONS", "#f6e6ea", "#9b2c48", ["Drought bulletins and\nanticipatory action triggers", "Abstraction permits and\nsafe pumping limits", "Borehole siting and\ndrilling depth", "Network design and\ndatum maintenance", "Water for production:\nlivestock and irrigation"])]
fig, ax = plt.subplots(figsize=(11, 5.6)); ax.set_xlim(0, 4); ax.set_ylim(-0.85, 5.75); ax.axis("off")
for i, (t, fc, ec, items) in enumerate(cols):
    ax.add_patch(FancyBboxPatch((i + 0.05, -0.05), 0.9, 5.45, boxstyle="round,pad=0.02,rounding_size=0.05", fc=fc, ec=ec, lw=1.2))
    ax.text(i + 0.5, 5.18, t, ha="center", va="center", fontsize=10.5, fontweight="bold", color=ec)
    for j, s in enumerate(items):
        y = 4.35 - j * 0.98
        ax.add_patch(FancyBboxPatch((i + 0.11, y - 0.36), 0.78, 0.72, boxstyle="round,pad=0.01,rounding_size=0.04", fc="white", ec=ec, lw=0.7))
        ax.text(i + 0.5, y, s, ha="center", va="center", fontsize=7.6)
    if i < 3:
        ax.add_patch(FancyArrowPatch((i + 0.96, 2.4), (i + 1.04, 2.4), arrowstyle="-|>", mutation_scale=18, color="#444", lw=1.5))
ax.add_patch(FancyArrowPatch((3.5, -0.1), (0.5, -0.1), connectionstyle="arc3,rad=-0.08", arrowstyle="-|>", mutation_scale=16, color="#777", lw=1.2, ls="--"))
ax.text(2, -0.68, "Feedback: decisions on network, abstraction and data standards improve the next data cycle (annual rerun, open code on GitHub and Zenodo)",
        ha="center", fontsize=8, color="#555", style="italic")
ax.text(2, 5.62, "SCIENCE  +  TECHNOLOGY  +  POLICY: an open groundwater digital twin for drylands", ha="center", fontsize=11, fontweight="bold")
plt.tight_layout(); plt.savefig(os.path.join(F, "fig02_framework.png"), dpi=200); plt.close()
