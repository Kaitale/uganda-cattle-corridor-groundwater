"""GemPy implicit model of the regolith (cover / weathered-fractured / fresh rock) for Kitgum and Pader.

Input : data/processed/aswa_layer_picks_dem.csv (from litho_prep.py), data/raw/dem/aswa_dem_110m_utm36n.tif
Output: figures/gempy_kitgum_section.png (W-E section, masked above the DEM and >6 km from boreholes)
"""
import pandas as pd, numpy as np, gempy as gp, rasterio, matplotlib, warnings, os, tempfile
matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
H=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(H,"..","data")
d=pd.read_csv(os.path.join(D,"processed","aswa_layer_picks_dem.csv")).dropna(subset=["fresh_top","ground_elev_m"])
d=d[(d.district.isin(["KITGUM","PADER"]))&(d.fresh_top>=d.cover_base)&(d.cover_base>0)&(d.fresh_top<100)]
d["z_cover"]=d.ground_elev_m-d.cover_base; d["z_fresh"]=d.ground_elev_m-d.fresh_top
# thin to 1 per 2 km cell to keep kriging tractable
d["cell"]=(d.x//2000).astype(int).astype(str)+"_"+(d.y//2000).astype(int).astype(str)
s=d.groupby("cell").median(numeric_only=True).reset_index()
print("boreholes",len(d),"cells",len(s))
xmin,xmax,ymin,ymax=s.x.min()-1000,s.x.max()+1000,s.y.min()-1000,s.y.max()+1000
zmin=float(s.z_fresh.quantile(0.01)-40); zmax=float(s.ground_elev_m.quantile(0.99)+20)
pts=pd.concat([pd.DataFrame(dict(X=s.x,Y=s.y,Z=s.z_cover,surface="base_cover")),pd.DataFrame(dict(X=s.x,Y=s.y,Z=s.z_fresh,surface="top_fresh"))])
cx,cy=s.x.median(),s.y.median()
ori=pd.DataFrame(dict(X=[cx,cx],Y=[cy,cy],Z=[s.z_cover.median(),s.z_fresh.median()],G_x=[0,0],G_y=[0,0],G_z=[1,1],surface=["base_cover","top_fresh"]))
T=tempfile.mkdtemp(); fp,fo=os.path.join(T,"points.csv"),os.path.join(T,"orientations.csv")
pts.to_csv(fp,index=False); ori.to_csv(fo,index=False)
geo=gp.create_geomodel(project_name="Kitgum_Pader",extent=[xmin,xmax,ymin,ymax,zmin,zmax],resolution=[60,60,120],
   importer_helper=gp.data.ImporterHelper(path_to_surface_points=fp,path_to_orientations=fo))
gp.map_stack_to_surfaces(geo,{"Cover":("base_cover",),"Weathered":("top_fresh",)})
gp.compute_model(geo)
lb=geo.solutions.raw_arrays.lith_block.reshape(60,60,120)
zs=np.linspace(zmin,zmax,120); xs=np.linspace(xmin,xmax,60); ys=np.linspace(ymin,ymax,60)
j=30
with rasterio.open(os.path.join(D,"raw","dem","aswa_dem_110m_utm36n.tif")) as r:
    topo=np.array([v[0] for v in r.sample([(x,ys[j]) for x in xs])])
sec=lb[:,j,:].T.astype(float)
from scipy.spatial import cKDTree
tr=cKDTree(d[["x","y"]].values)
for i,x in enumerate(xs):
    sec[zs>topo[i],i]=np.nan
    if tr.query([x,ys[j]])[0]>6000: sec[:,i]=np.nan
fig,ax=plt.subplots(figsize=(11,4))
from matplotlib.colors import ListedColormap
ax.pcolormesh(xs/1000,zs,sec,cmap=ListedColormap(["#d9b38c","#8fbcd4","#7a7a7a"]),shading="auto",vmin=1,vmax=3)
ax.plot(xs/1000,topo,"k",lw=1)
near=d[(d.y-ys[j]).abs()<3000]
for _,b in near.iterrows():
    ax.plot([b.x/1000]*2,[b.ground_elev_m,b.ground_elev_m-b.depth],c="k",lw=0.6)
ax.set_xlim(near.x.min()/1000-5, near.x.max()/1000+5); ax.set_ylim(zs.min(), zmax)
ax.set_xlabel("Easting, km (UTM 36N)"); ax.set_ylabel("Elevation, m a.s.l.")
ax.set_title(f"Kitgum–Pader regolith model (GemPy), W–E section at N {ys[j]/1000:.0f} km: cover (tan), weathered/fractured aquifer (blue), fresh rock (grey)",fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(H,"..","figures","gempy_kitgum_section.png"),dpi=110)
u,c=np.unique(lb,return_counts=True); print(dict(zip(u,c)))
