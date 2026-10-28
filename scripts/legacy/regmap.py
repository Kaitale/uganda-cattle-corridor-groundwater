import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
d=pd.read_csv('gwout/aswa_layer_picks.csv').dropna(subset=['fresh_top'])
d=d[(d.fresh_top>=d.cover_base)&(d.fresh_top<100)]
xy=d[['x','y']].values; tree=cKDTree(xy)
gx,gy=np.meshgrid(np.linspace(d.x.min(),d.x.max(),220),np.linspace(d.y.min(),d.y.max(),180))
q=np.c_[gx.ravel(),gy.ravel()]; dist,idx=tree.query(q,k=12)
w=1/np.maximum(dist,50)**2
def idw(v): return (np.sum(w*v[idx],1)/np.sum(w,1)).reshape(gx.shape)
mask=(dist[:,0]>8000).reshape(gx.shape)
reg=idw(d.fresh_top.values); reg[mask]=np.nan
cov=idw(d.cover_base.values); cov[mask]=np.nan
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
for a,z,t in [(ax[0],reg,'Depth to fresh rock (regolith thickness), m'),(ax[1],reg-cov,'Weathered/fractured zone thickness, m')]:
    im=a.pcolormesh(gx/1000,gy/1000,z,cmap='viridis',shading='auto',vmin=0,vmax=50); a.scatter(d.x/1000,d.y/1000,s=2,c='k',alpha=.3)
    a.set_title(t,fontsize=10); a.set_xlabel('Easting, km (UTM 36N)'); a.set_aspect('equal'); plt.colorbar(im,ax=a,shrink=.8)
ax[0].set_ylabel('Northing, km')
fig.suptitle(f'Aswa catchment: {len(d)} borehole logs, inverse-distance interpolation (preliminary)',fontsize=11)
plt.tight_layout(); plt.savefig('gwout/aswa_regolith_maps.png',dpi=110)
print(len(d), np.nanmedian(reg), np.nanmedian(reg-cov))
