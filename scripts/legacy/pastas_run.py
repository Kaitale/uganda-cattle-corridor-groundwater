import pandas as pd, pastas as ps, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, warnings, sys
warnings.filterwarnings('ignore'); ps.set_log_level("ERROR")
D='cc/Cattle corridor/'
st={'Soroti':('Soroti.csv',),'Kangole':('824007 Kangole.csv',),'Ntungamo':('GW_Ntungamo.csv',),'Moroto':('GW_824011_Moroto.csv',),'Nakasongola':('NAKASONGOLA.csv',),'Serere':('Serere.csv',),'Isingiro':('GW_814005_Isingiro.csv',),'Rwebisengo':('Rwebisengo.csv',)}
res=[]
for name,(f,) in st.items():
    d=pd.read_csv(D+f,encoding='latin1'); dc=d.columns[0]; d[dc]=pd.to_datetime(d[dc],errors='coerce'); d=d.dropna(subset=[dc]).set_index(dc).sort_index()
    d=d[~d.index.duplicated()]
    h=pd.to_numeric(d.iloc[:,0],errors='coerce').dropna()
    h=-h  # m below ground -> negative head
    r=pd.to_numeric(d.iloc[:,1],errors='coerce').clip(lower=0)
    r=r.asfreq('D')
    # use rain period only; fill missing rain days with 0 inside coverage window of rain
    rv=r.dropna(); r=r.loc[rv.index.min():rv.index.max()].fillna(0)
    h=h.loc[r.index.min():r.index.max()]
    h=h[(h-h.median()).abs()<5*h.std()]
    if len(h)<300: print(name,'skip',len(h)); continue
    ml=ps.Model(h,name=name)
    sm=ps.StressModel(r/1000,ps.Gamma(),name='rain',settings='prec')
    ml.add_stressmodel(sm); ml.add_noisemodel(ps.ArNoiseModel()) if hasattr(ps,'ArNoiseModel') else None
    try:
        ml.solve(report=False)
    except Exception as e: print(name,e); continue
    ev=ml.stats.evp(); rmse=ml.stats.rmse()
    t50=sm.rfunc.get_tmax(ml.get_parameters('rain'),cutoff=0.5); t95=sm.rfunc.get_tmax(ml.get_parameters('rain'),cutoff=0.95)
    gain=sm.rfunc.gain(ml.get_parameters('rain'))
    res.append([name,h.index.min().date(),h.index.max().date(),len(h),round(ev,1),round(rmse,2),round(t50),round(t95),round(gain,1)])
    if name=='Soroti':
        fig=ml.plots.results(figsize=(10,7)); plt.savefig('gwout/pastas_soroti.png',dpi=110,bbox_inches='tight'); plt.close('all')
print(pd.DataFrame(res,columns=['station','start','end','n','EVP%','RMSE m','t50 d','t95 d','gain m per m rain']).to_string())
pd.DataFrame(res,columns=['station','start','end','n','EVP%','RMSE m','t50 d','t95 d','gain m per m rain']).to_csv('gwout/pastas_summary.csv',index=False)
