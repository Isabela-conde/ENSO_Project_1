#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
datadir="/Volumes/Raid2/tamura/GODAS/data"
u_da=xr.open_dataset(f"{datadir}/ucur.merged.nc").sel(lat=slice(-5,5)).ucur
tmp_da=xr.open_dataset(f"{datadir}/pottmp.merged.nc").sel(lat=slice(-5,5)).pottmp

ulon,ulat=u_da.lon.values,u_da.lat.values
tlon,tlat=tmp_da.lon.values,tmp_da.lat.values
# %%
def xr_xdev(var,lat,lon,xdim="lon",rad_earth=6371e3):
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    dlat=lat_rad[1]-lat_rad[0]
    dlon=lon_rad[1]-lon_rad[0]
    
    dx = (rad_earth * np.cos(lat_rad) * dlon * xr.ones_like(var[xdim])).astype(np.float32)
    # return (var.shift({xdim:-1}) - var.shift({xdim:1})) / (2 * dx)
    return (var.shift({xdim:-1}) - var) / (dx)

def xr_ydev(var,lat,ydim="lat",rad_earth=6371e3):
    lat_rad = np.deg2rad(lat)
    dlat=lat_rad[1]-lat_rad[0]
    
    dy = (rad_earth * dlat * xr.ones_like(var[ydim])).astype(np.float32)
    return (var.shift({ydim:-1}) - var.shift({ydim:1})) / (2 * dy)
#%%--------------------------------------------------------
# zonal advection
#----------------------------------------------------------
dtdx=xr_xdev(tmp_da,tmp_da.lat,tmp_da.lon)
dtdx["lon"]=ulon[1:]
dtdx=dtdx.interp(lat=ulat)
dtdx.name="dtdx"
dtdx.attrs["units"]="K /m"

# udtdx=u_da*dtdx
# udtdx.name="udtdx"
# udtdx.attrs["units"]="K /s"
#%%
dtdx.to_netcdf(f"{datadir}/dtdx.merged.nc")
# udtdx.to_netcdf(f"{datadir}/udtdx.merged.nc")
#%%
def rm_bar(var):
    clim = var.groupby('time.month').mean('time')
    return var.groupby('time.month') - clim
#%%
dtdx_bar=dtdx.groupby('time.month').mean('time')
dtdx_prm=dtdx.groupby('time.month') - dtdx_bar
u_bar=u_da.groupby('time.month').mean('time')
u_prm=u_da.groupby('time.month') - u_bar
#%%
ubar_dtdxprm=u_bar*dtdx_prm
ubar_dtdxprm.name="ubar dt'dx"
ubar_dtdxprm.attrs["units"]="K /s"

uprm_dtdxbar=u_prm*dtdx_bar
uprm_dtdxbar.name="u' dt_bardx"
#%%
# Vertical integration
dz=10 # in meters, (10m for GODAS)
zsum_ubar_dtdxprm=ubar_dtdxprm.sum("level")*dz
zsum_uprm_dtdxbar=uprm_dtdxbar.sum("level")*dz
zsum_ubar_dtdxprm.to_netcdf(f"{datadir}/zsum_ubar_dtdxprm.merged.nc")
zsum_uprm_dtdxbar.to_netcdf(f"{datadir}/zsum_uprm_dtdxbar.merged.nc")
# %%----------------------------------------------------------
w_da=xr.open_dataset(f"{datadir}/dzdt.merged.nc").sel(lat=slice(-5,5)).dzdt
wlev=w_da.level.values
tlev=tmp_da.level.values
# %%
def xr_zdev(var,dlev,zdim="level"):
    dvardz= -(var.shift({zdim:-1}) - var) / dlev
    return dvardz
# %%
dlev=tmp_da.level.shift(level=-1)-tmp_da.level
dtdz=xr_zdev(tmp_da,dlev,zdim="level")
dtdz["level"]=wlev[:]
# %%
