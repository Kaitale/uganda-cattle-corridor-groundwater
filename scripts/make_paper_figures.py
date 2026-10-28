"""Non-map figures: Fig. 7 (GemPy section copy), 8 (Pastas Rwonyo copy), 9 (SGI panel), A1 (borehole statistics).

Run after pastas_trend_drought.py and gempy_kitgum.py."""
import os, pandas as pd, numpy as np, matplotlib, warnings
matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
H=os.path.dirname(os.path.abspath(__file__)); R=os.path.join(H,"..","data"); F=os.path.join(H,"..","figures","paper"); os.makedirs(F,exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
st=pd.read_csv(os.path.join(R,"raw","Stations.csv"))
sm=pd.read_csv(os.path.join(R,"processed","pastas_trend_drought_summary.csv"))
# Fig 2: national borehole statistics
a=pd.read_excel(os.path.join(R,"raw","boreholes","BOREHOLE DATA.xlsx"))
cols=[("TotalDrilledDepth","Drilled depth (m)",(0,150)),("DepthToBedRock","Depth to bedrock (m)",(0,100)),("StaticWaterLevel","Static water level (m bgl)",(0,60)),("ConstantDischargeYield","Yield (m³/h)",(0,10))]
fig,axs=plt.subplots(1,4,figsize=(10,2.6))
for ax,(c,l,xl) in zip(axs,cols):
    x=pd.to_numeric(a[c],errors="coerce").dropna(); x=x[(x>=xl[0])&(x<=xl[1])]
    ax.hist(x,bins=40,color="#4a7aa8"); ax.axvline(x.median(),c="k",ls="--",lw=1)
    ax.set_xlabel(l); ax.set_title(f"n = {len(x):,}, median {x.median():.1f}",fontsize=8)
axs[0].set_ylabel("Boreholes")
plt.tight_layout(); plt.savefig(os.path.join(F,"figA1_borehole_statistics.png"),dpi=200); plt.close()
# Fig 5 and 6: copy model figures
import shutil
for src,dst in [(("figures","gempy_kitgum_section.png"),"fig08_gempy_kitgum_section.png"),(("figures","pastas_trend","Rwonyo.png"),"fig09_pastas_rwonyo.png")]:
    p=os.path.join(H,"..",*src)
    if os.path.exists(p): shutil.copy(p,os.path.join(F,dst))
# Fig 7: SGI panel
import pastas as ps
pick=["Rwonyo","Nakasongola","Soroti","Kangole","Bombo Lukole"]
fig,axs=plt.subplots(len(pick),1,figsize=(9,7.5),sharex=True)
for ax,n in zip(axs,pick):
    s=st[st.station_name==n].iloc[0]
    d=pd.read_csv(os.path.join(R,"raw","levels",s.data_file),encoding="latin1",on_bad_lines="skip")
    dc=d.columns[0]; d[dc]=pd.to_datetime(d[dc],errors="coerce"); d=d.dropna(subset=[dc]).set_index(dc).sort_index(); d=d[~d.index.duplicated()]
    dep=pd.to_numeric(d.iloc[:,0],errors="coerce").dropna(); dep=dep[(dep-dep.median()).abs()<5*dep.std()]
    mon=dep.resample("MS").mean().dropna(); sgi=ps.stats.sgi(-mon)
    ax.bar(sgi.index,sgi.clip(upper=0),width=25,color="#b2182b"); ax.bar(sgi.index,sgi.clip(lower=0),width=25,color="#2166ac")
    ax.axhline(-1,ls="--",c="#ef8a62",lw=0.8); ax.axhline(-1.5,ls="--",c="#b2182b",lw=0.8)
    ax.set_ylim(-2.6,2.6); ax.set_ylabel("SGI"); ax.text(0.01,0.85,n,transform=ax.transAxes,fontweight="bold")
plt.tight_layout(); plt.savefig(os.path.join(F,"fig10_sgi_panel.png"),dpi=200); plt.close()
print("figures done")
