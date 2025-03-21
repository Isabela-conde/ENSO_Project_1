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

ulon,ulat=u_da.lon.values,u_da.lat.values
tlon,tlat=tmp_da.lon.values,tmp_da.lat.values
wlev=w_da.level.values
tlev=tmp_da.level.values
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
def xr_zdev(var,dlev,zdim="level"):
    dvardz= -(var.shift({zdim:-1}) - var) / dlev
    return dvardz
#%%--------------------------------------------------------
# zonal gradient
#----------------------------------------------------------
dtdx=xr_xdev(tmp_da,tmp_da.lat,tmp_da.lon)
dtdx["lon"]=ulon[1:] # convert tlon to match ulon
dtdx=dtdx.interp(lat=ulat) # interpolate tlat to ulat
dtdx.name="dtdx"
dtdx.attrs["units"]="K /m"
#%%
dtdx.to_netcdf(f"{datadir}/dtdx.merged.nc")
#%%--------------------------------------------------------
# meridional gradient
#----------------------------------------------------------
dtdy=xr_ydev(tmp_da,tmp_da.lat)
#%%
dtdy["lat"]=ulat[:] # convert tlat to match ulat
dtdy=dtdy.interp(lon=ulon) # interpolate tlon to ulon
dtdy.name="dtdy"
dtdy.attrs["units"]="K /m"
#%%
dtdy.to_netcdf(f"{datadir}/dtdy.merged.nc")
# %%--------------------------------------------------------
# vertical gradient
#----------------------------------------------------------
dlev=tmp_da.level.shift(level=-1)-tmp_da.level
dtdz=xr_zdev(tmp_da,dlev,zdim="level")[:,:-1]
dtdz["level"]=wlev[:]
dtdz=dtdz.interp(lat=w_da.lat,method="linear")
dtdz.name="dtdz"
dtdz.attrs["units"]="K/m"
# %%
dtdz.to_netcdf(f"{datadir}/dtdz.merged.nc")
# %%
