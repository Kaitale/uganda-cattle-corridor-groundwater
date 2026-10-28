import pandas as pd, numpy as np, gempy as gp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, warnings
warnings.filterwarnings('ignore')
d=pd.read_csv('gwout/aswa_layer_picks.csv').dropna(subset=['fresh_top'])
d=d[(d.cover_base>0)&(d.fresh_top>=d.cover_base)]
s=d.sample(min(300,len(d)),random_state=1)
xmin,xmax,ymin,ymax=s.x.min()-2000,s.x.max()+2000,s.y.min()-2000,s.y.max()+2000
pts=pd.concat([pd.DataFrame(dict(X=s.x,Y=s.y,Z=-s.cover_base,surface='base_cover')),pd.DataFrame(dict(X=s.x,Y=s.y,Z=-s.fresh_top,surface='top_fresh'))])
cx,cy=(xmin+xmax)/2,(ymin+ymax)/2
ori=pd.DataFrame(dict(X=[cx,cx],Y=[cy,cy],Z=[-13,-24],G_x=[0,0],G_y=[0,0],G_z=[1,1],surface=['base_cover','top_fresh']))
pts.to_csv('gwout/gp_points.csv',index=False); ori.to_csv('gwout/gp_ori.csv',index=False)
geo=gp.create_geomodel(project_name='Aswa_regolith',extent=[xmin,xmax,ymin,ymax,-140,0],resolution=[50,50,40],
    importer_helper=gp.data.ImporterHelper(path_to_surface_points='gwout/gp_points.csv',path_to_orientations='gwout/gp_ori.csv'))
gp.map_stack_to_surfaces(geo,{'Regolith':('base_cover','top_fresh')})
gp.compute_model(geo)
import gempy_viewer as gpv
p=gpv.plot_2d(geo,direction=['y'],cell_number='mid',show_data=True,show_boundaries=False,show=False)
plt.gcf().set_size_inches(10,4); plt.savefig('gwout/gempy_section.png',dpi=110,bbox_inches='tight')
import numpy as np
lb=geo.solutions.raw_arrays.lith_block.reshape(50,50,40)
print('done',len(s),np.unique(lb,return_counts=True))
