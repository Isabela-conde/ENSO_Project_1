#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
datadir="/Volumes/Raid2/tamura/GODAS/data"
v_da=xr.open_dataset(f"{datadir}/vcur.merged.nc").sel(lat=slice(-6,6)).vcur
tmp_da=xr.open_dataset(f"{datadir}/pottmp.merged.nc").sel(lat=slice(-6,6)).pottmp

ulon,ulat=v_da.lon.values,v_da.lat.values
tlon,tlat=tmp_da.lon.values,tmp_da.lat.values
# %%
def xr_xdev(var,lat,lon,xdim="lon",rad_earth=6371e3):
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    dlon=lon_rad[1]-lon_rad[0]
    
    dx = (rad_earth * np.cos(lat_rad) * dlon * xr.ones_like(var[xdim])).astype(np.float32)
    # return (var.shift({xdim:-1}) - var.shift({xdim:1})) / (2 * dx)
    return (var.shift({xdim:-1}) - var) / (dx)

def xr_ydev(var,lat,ydim="lat",rad_earth=6371e3):
    lat_rad = np.deg2rad(lat)
    # dlat=lat_rad[1]-lat_rad[0]
    dlat=lat_rad.shift({ydim:-1})-lat_rad
    
    # dy = (rad_earth * dlat * xr.ones_like(var[ydim])).astype(np.float32)
    dy = (rad_earth * dlat).astype(np.float32)
    return (var.shift({ydim:-1}) - var) / (dy)

def xr_ydev(var,lat,ydim="lat",rad_earth=6371e3):
    lat_rad = np.deg2rad(lat)
    # dlat=lat_rad[1]-lat_rad[0]
    dlat=lat_rad.shift({ydim:-1})-lat_rad
    
    # dy = (rad_earth * dlat * xr.ones_like(var[ydim])).astype(np.float32)
    dy = (rad_earth * dlat).astype(np.float32)
    return (var.shift({ydim:-1}) - var) / (dy)
#%%--------------------------------------------------------
# meridional advection
#----------------------------------------------------------
dtdy=xr_ydev(tmp_da,tmp_da.lat)
#%%
dtdy["lat"]=ulat[:]
dtdy=dtdy.interp(lon=ulon)
dtdy.name="dtdy"
dtdy.attrs["units"]="K /m"
#%%
dtdy.to_netcdf(f"{datadir}/dtdy.merged.nc")
# udtdy.to_netcdf(f"{datadir}/udtdy.merged.nc")
#%%
def rm_bar(var):
    clim = var.groupby('time.month').mean('time')
    return var.groupby('time.month') - clim
#%%
dtdy_bar=dtdy.groupby('time.month').mean('time')
dtdy_prm=dtdy.groupby('time.month') - dtdy_bar
v_bar=v_da.groupby('time.month').mean('time')
v_prm=v_da.groupby('time.month') - v_bar
#%%
vbar_dtdyprm=v_bar*dtdy_prm
vbar_dtdyprm.name="vbar dt'dy"
vbar_dtdyprm.attrs["units"]="K /s"

vprm_dtdybar=v_prm*dtdy_bar
vprm_dtdybar.name="v' dt_bardy"
#%%
# Vertical integration
dz=10 # in meters, (10m for GODAS)
zsum_vbar_dtdyprm=vbar_dtdyprm.sum("level")*dz
zsum_vprm_dtdybar=vprm_dtdybar.sum("level")*dz
zsum_vbar_dtdyprm.to_netcdf(f"{datadir}/zsum_vbar_dtdyprm.merged.nc")
zsum_vprm_dtdybar.to_netcdf(f"{datadir}/zsum_vprm_dtdybar.merged.nc")
# %%
