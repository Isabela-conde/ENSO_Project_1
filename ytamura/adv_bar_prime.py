#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
datadir="/Volumes/Raid2/tamura/GODAS/data"
u_da=xr.open_dataset(f"{datadir}/ucur.merged.nc").sel(lat=slice(-5,5)).ucur
v_da=xr.open_dataset(f"{datadir}/vcur.merged.nc").sel(lat=slice(-6,6)).vcur
tmp_da=xr.open_dataset(f"{datadir}/pottmp.merged.nc").sel(lat=slice(-5,5)).pottmp
dtdx=xr.open_dataarray(f"{datadir}/dtdx.merged.nc")
dtdy=xr.open_dataarray(f"{datadir}/dtdy.merged.nc")

ulon,ulat=u_da.lon.values,u_da.lat.values
tlon,tlat=tmp_da.lon.values,tmp_da.lat.values
# %%
p1=[np.datetime64("1980-01-01"),np.datetime64("2000-12-01")]
p2=[np.datetime64("2001-01-01"),np.datetime64("2021-12-01")]
#%%
# calculate monthly climatology and anomalies
# bar: climatology, prm: anomaly
yr1,yr2=1980,2000
t1,t2=np.datetime64(f"{yr1}-01-01"),np.datetime64(f"{yr2}-02-01")
ubar=u_da.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
vbar=v_da.groupby('time.month').mean('time')
dtdxbar=dtdx.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
dtdybar=dtdy.sel(time=slice(t1,t2)).groupby('time.month').mean('time')

uprm=u_da.groupby('time.month') - ubar
vprm=v_da.groupby('time.month') - vbar
dtdxprm=dtdx.groupby('time.month') - dtdxbar
dtdyprm=dtdy.groupby('time.month') - dtdybar
#%%
ubar_dtdxprm=ubar*dtdxprm
ubar_dtdxprm.name="ubar dt'dx"
ubar_dtdxprm.attrs["units"]="K /s"

vbar_dtdyprm=vbar*dtdyprm
vbar_dtdyprm.name="vbar dt'dy"
vbar_dtdyprm.attrs["units"]="K /s"

uprm_dtdxbar=uprm*dtdxbar
uprm_dtdxbar.name="u' dt_bardx"
uprm_dtdxbar.attrs["units"]="K /s"
vprm_dtdybar=vprm*dtdybar
vprm_dtdybar.name="v' dtdybar"
#%%
# Vertical integration
dz=10 # in meters, (10m for GODAS)
zsum_ubar_dtdxprm=ubar_dtdxprm.sum("level")*dz  
zsum_uprm_dtdxbar=uprm_dtdxbar.sum("level")*dz
zsum_vbar_dtdyprm=vbar_dtdyprm.sum("level")*dz
zsum_vprm_dtdybar=vprm_dtdybar.sum("level")*dz
# %%
lon1_n3,lon2_n3,lat1_n3,lat2_n3=210,270,-5,5
lon1_n4,lon2_n4,lat1_n4,lat2_n4=160,210,-5,5
# %%
dlat=np.deg2rad(1/3)
dlon=np.deg2rad(1)
rad_earth=6371e3
dy = (rad_earth * dlat)
dx = (rad_earth * np.cos(np.deg2rad(ulat)) * dlon)
# Sum over the region
# Niño 3
ubar_dtdxprm_n3=dy*zsum_ubar_dtdxprm\
    .sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)).sum(["lat","lon"])
uprm_dtdxbar_n3=dy*zsum_uprm_dtdxbar\
    .sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)).sum(["lat","lon"])
vbar_dtdyprm_n3=dx*zsum_vbar_dtdyprm\
    .sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)).sum(["lat","lon"])
vprm_dtdybar_n3=dx*zsum_vprm_dtdybar\
    .sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)).sum(["lat","lon"])
# Niño 4
ubar_dtdxprm_n4=dy*zsum_ubar_dtdxprm\
    .sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)).sum(["lat","lon"])
uprm_dtdxbar_n4=dy*zsum_uprm_dtdxbar\
    .sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)).sum(["lat","lon"])
vbar_dtdyprm_n4=dx*zsum_vbar_dtdyprm\
    .sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)).sum(["lat","lon"])
vprm_dtdybar_n4=dx*zsum_vprm_dtdybar\
    .sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)).sum(["lat","lon"])