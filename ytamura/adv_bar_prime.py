#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
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
# %%
p1=[np.datetime64("1980-01-01"),np.datetime64("2000-12-01")]
p2=[np.datetime64("2001-01-01"),np.datetime64("2021-12-01")]
#%%
u_dtdx=u_da*dtdx
v_dtdy=v_da*dtdy
w_dtdz=w_da*dtdz
#%%
# calculate monthly climatology and anomalies
# bar: climatology, prm: anomaly
yr1,yr2=1980,2000
t1,t2=np.datetime64(f"{yr1}-01-01"),np.datetime64(f"{yr2}-02-01")
ubar=u_da.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
vbar=v_da.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
wbar=w_da.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
dtdxbar=dtdx.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
dtdybar=dtdy.sel(time=slice(t1,t2)).groupby('time.month').mean('time')
dtdzbar=dtdz.sel(time=slice(t1,t2)).groupby('time.month').mean('time')

uprm=u_da.groupby('time.month') - ubar
vprm=v_da.groupby('time.month') - vbar
wprm=w_da.groupby('time.month') - wbar
dtdxprm=dtdx.groupby('time.month') - dtdxbar
dtdyprm=dtdy.groupby('time.month') - dtdybar
dtdzprm=dtdz.groupby('time.month') - dtdzbar
#%%
ubar_dtdxprm = ubar * dtdxprm.groupby("time.month")
ubar_dtdxprm.name = "ubar_dtdxprm"
ubar_dtdxprm.attrs["units"] = "K /s"

vbar_dtdyprm = vbar * dtdyprm.groupby("time.month")
vbar_dtdyprm.name = "vbar_dtdyprm"
vbar_dtdyprm.attrs["units"] = "K /s"

wbar_dtdzprm = wbar * dtdzprm.groupby("time.month")
wbar_dtdzprm.name = "wbar_dtdzprm"
wbar_dtdzprm.attrs["units"] = "K /s"

uprm_dtdxbar = uprm.groupby("time.month") * dtdxbar
uprm_dtdxbar.name = "uprm_dtdxbar"
uprm_dtdxbar.attrs["units"] = "K /s"

vprm_dtdybar = vprm.groupby("time.month") * dtdybar
vprm_dtdybar.name = "vprm_dtdybar"
vprm_dtdybar.attrs["units"] = "K /s"

wprm_dtdzbar = wprm.groupby("time.month") * dtdzbar
wprm_dtdzbar.name = "wprm_dtdzbar"
wprm_dtdzbar.attrs["units"] = "K /s"
#%%
# Vertical mean, not integration
# all layer has the same thickness
dz=10 # in meters, (10m for GODAS)
dsec=3600*24*30 # in seconds, (30 days)
ubar_dtdxprm_zm=ubar_dtdxprm.mean("level")*dsec#*dz  
uprm_dtdxbar_zm=uprm_dtdxbar.mean("level")*dsec#*dz
vbar_dtdyprm_zm=vbar_dtdyprm.mean("level")*dsec#*dz
vprm_dtdybar_zm=vprm_dtdybar.mean("level")*dsec#*dz
wbar_dtdzprm_zm=wbar_dtdzprm.mean("level")*dsec#*dz
wprm_dtdzbar_zm=wprm_dtdzbar.mean("level")*dsec#*dz
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
# %%
# Save DataArrays to netCDF files
output_dir = "/Volumes/Raid2/tamura/GODAS/mhb"  # Use the data directory as the output directory

variables_n3 = [ubar_dtdxprm_n3, uprm_dtdxbar_n3, vbar_dtdyprm_n3, vprm_dtdybar_n3, wbar_dtdzprm_n3, wprm_dtdzbar_n3]
variables_n4 = [ubar_dtdxprm_n4, uprm_dtdxbar_n4, vbar_dtdyprm_n4, vprm_dtdybar_n4, wbar_dtdzprm_n4, wprm_dtdzbar_n4]
variables_zm = [ubar_dtdxprm_zm, uprm_dtdxbar_zm, vbar_dtdyprm_zm, vprm_dtdybar_zm, wbar_dtdzprm_zm, wprm_dtdzbar_zm]

for var in variables_zm:
    var.to_netcdf(f"{output_dir}/{var.name}_zm_{yr1}-{yr2}.nc")

# for var in variables_n4:
#     var.to_netcdf(f"{output_dir}/{var.name}_n4_{yr1}-{yr2}.nc")
#%%
u_dtdx.to_netcdf(f"{output_dir}/u_dtdx.nc")
v_dtdy.to_netcdf(f"{output_dir}/v_dtdy.nc")
w_dtdz.to_netcdf(f"{output_dir}/w_dtdz.nc")
# %%==================================================================================
# Check data
import tamdraw as dr # pip install git+https://github.com/y-tamura/tamdraw.git
import cartopy.crs as ccrs
#%%===================================================================================
t1,t2=np.datetime64("1998-10-01"),np.datetime64("1998-11-01")
var=5e5*dtdx[:,0]#ubar_dtdxprm
nax=var.sel(time=slice(t1,t2)).time.size
ncol=2; nrow=(nax-1)//ncol+1
fig,axes=plt.subplots(nrow,ncol,figsize=(4*ncol,1.6*nrow),
                      subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)},
                      layout="constrained")
fig.suptitle(var.name)
for it,t in enumerate(var.sel(time=slice(t1,t2)).time.values):
    ax=np.ravel(axes)[it]
    c=dr.axplot_hrz_field_double(ax,var.sel(time=t),tmp_da.sel(time=t)[0]-273,
                               -1,1,.2, 0,40,2, 120,280,-20,20, 30,10,
                               title=str(t)[:7])
fig.colorbar(c, ax=axes, orientation='horizontal',aspect=40, shrink=.6)
#%%
t1,t2=np.datetime64("1998-08-01"),np.datetime64("1998-10-01")
var=wprm_dtdzbar_zm#ubar_dtdxprm_zm
nax=var.sel(time=slice(t1,t2)).time.size
ncol=2; nrow=(nax-1)//ncol+1
fig,axes=plt.subplots(nrow,ncol,figsize=(4*ncol,1.6*nrow),
                      subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)},
                      layout="constrained")
fig.suptitle(var.name)
for it,t in enumerate(var.sel(time=slice(t1,t2)).time.values):
    ax=np.ravel(axes)[it]
    c=dr.axplot_hrz_field_double(ax,var.sel(time=t),tmp_da.sel(time=t)[0]-273,
                               -2,2,.5, 0,40,2, 120,280,-20,20, 30,10,
                               title=str(t)[:7])
fig.colorbar(c, ax=axes, orientation='horizontal',aspect=40, shrink=.6)

# %%
