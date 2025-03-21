#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import tamfunc as tf
#%%
datadir="/Volumes/Raid2/tamura/GODAS/data"
u_da=xr.open_dataset(f"{datadir}/ucur.merged.nc").sel(lat=slice(-20,20)).ucur
v_da=xr.open_dataset(f"{datadir}/vcur.merged.nc").sel(lat=slice(-20,20)).vcur
w_da=xr.open_dataset(f"{datadir}/dzdt.merged.nc").sel(lat=slice(-20,20)).dzdt
tmp_da=xr.open_dataset(f"{datadir}/pottmp.merged.nc").sel(lat=slice(-20,20)).pottmp
dtdx=xr.open_dataarray(f"{datadir}/dtdx.merged.nc")
dtdy=xr.open_dataarray(f"{datadir}/dtdy.merged.nc")
dtdz=xr.open_dataarray(f"{datadir}/dtdz.merged.nc")

ulon,ulat=u_da.lon.values,u_da.lat.values
tlon,tlat=tmp_da.lon.values,tmp_da.lat.values
wlev=w_da.level.values
tlev=tmp_da.level.values

datadir="/Volumes/Raid2/tamura/GODAS/mhb"
u_dtdx=xr.open_dataarray(f"{datadir}/u_dtdx.nc")
v_dtdy=xr.open_dataarray(f"{datadir}/v_dtdy.nc")
w_dtdz=xr.open_dataarray(f"{datadir}/w_dtdz.nc")
# %%
yr1 = 1980
yr2 = 2000

ubar_dtdxprm_zm = xr.open_dataarray(f"{datadir}/ubar_dtdxprm_zm_{yr1}-{yr2}.nc")
uprm_dtdxbar_zm = xr.open_dataarray(f"{datadir}/uprm_dtdxbar_zm_{yr1}-{yr2}.nc")
vbar_dtdyprm_zm = xr.open_dataarray(f"{datadir}/vbar_dtdyprm_zm_{yr1}-{yr2}.nc")
vprm_dtdybar_zm = xr.open_dataarray(f"{datadir}/vprm_dtdybar_zm_{yr1}-{yr2}.nc")
wbar_dtdzprm_zm = xr.open_dataarray(f"{datadir}/wbar_dtdzprm_zm_{yr1}-{yr2}.nc")
wprm_dtdzbar_zm = xr.open_dataarray(f"{datadir}/wprm_dtdzbar_zm_{yr1}-{yr2}.nc")
# %%
u_dtdx_prm = tf.rm_monthlyclim(u_dtdx)
v_dtdy_prm = tf.rm_monthlyclim(v_dtdy)
w_dtdz_prm = tf.rm_monthlyclim(w_dtdz)
#%%
# Temperature tendency
tmp_prm_zm=tf.rm_monthlyclim(tmp_da.mean("level"))
dtdt=tmp_prm_zm.diff("time",label="lower")
# %%
lon1_n3,lon2_n3,lat1_n3,lat2_n3=210,270,-5,5
lon1_n4,lon2_n4,lat1_n4,lat2_n4=160,210,-5,5
# %%
# dlat=np.deg2rad(1/3)
# dlon=np.deg2rad(1)
# rad_earth=6371e3
# dy = (rad_earth * dlat)
# dx = (rad_earth * np.cos(xr.DataArray(np.deg2rad(ulat),coords={"lat":ulat})) * dlon)
cos_lat=np.cos(xr.DataArray(np.deg2rad(ulat),coords={"lat":ulat}))
wcos_lat=cos_lat/cos_lat.median()
def wcoslat_amean(var,wcos_lat=wcos_lat,mdims=("lat","lon")):
    return (var*wcos_lat).mean(mdims)
# Average over the region
# Niño 3
ubar_dtdxprm_n3=wcoslat_amean(ubar_dtdxprm_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
uprm_dtdxbar_n3=wcoslat_amean(uprm_dtdxbar_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
vbar_dtdyprm_n3=wcoslat_amean(vbar_dtdyprm_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
vprm_dtdybar_n3=wcoslat_amean(vprm_dtdybar_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
wbar_dtdzprm_n3=wcoslat_amean(wbar_dtdzprm_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
wprm_dtdzbar_n3=wcoslat_amean(wprm_dtdzbar_zm.sel(lat=slice(lat1_n3,lat2_n3),lon=slice(lon1_n3,lon2_n3)))
# Niño 4
ubar_dtdxprm_n4=wcoslat_amean(ubar_dtdxprm_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
uprm_dtdxbar_n4=wcoslat_amean(uprm_dtdxbar_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
vbar_dtdyprm_n4=wcoslat_amean(vbar_dtdyprm_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
vprm_dtdybar_n4=wcoslat_amean(vprm_dtdybar_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
wbar_dtdzprm_n4=wcoslat_amean(wbar_dtdzprm_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
wprm_dtdzbar_n4=wcoslat_amean(wprm_dtdzbar_zm.sel(lat=slice(lat1_n4,lat2_n4),lon=slice(lon1_n4,lon2_n4)))
#%%
ep_yrs=[1982,1986,1987,1991,1997]
