# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 15:37:21 2022

@author: sauda
"""
#Investigating Hurricane Harvey Rainfall with ERA5

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
from netCDF4 import Dataset
import glob


from geocat.viz import util as gvutil
from geocat.viz import cmaps as gvcmaps

file_h = 'real_assign_4_data.nc'
folder_h = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/4/'
filename = folder_h+file_h
ds_h = xr.open_dataset(filename, decode_times=True)

mslp = ds_h.msl.values/100
time = ds_h.time.values
lats = ds_h.latitude.values
lons = ds_h.longitude.values
tp = ds_h.tp.values
[a,b,c] = np.where(tp==np.amax(tp))

mslp2 = ds_h.msl.isel(time=a).drop('time')/100 #time='2019-01-01'
mslp_plot = np.squeeze(mslp)
tp_2 = ds_h.tp.isel(time=a).drop('time')*1000 #time='2019-01-01'
tp_2 = np.squeeze(tp_2)
ds = ds_h

#extracting the day, month & year
months = ds_h.time.dt.month.values
years = ds_h.time.dt.year.values
days = ds_h.time.dt.day.values
hours = ds_h.time.dt.hour.values
dates = np.transpose(np.stack((years,months,days,hours)))
dates_plot = dates[a,:]

#Map of tp at the time of peak rainfall

cmap = gvcmaps.precip4_11lev #Using the color map to define colors for contours 

plt.figure(figsize=(10,8))

which_plot = np.where((dates[:,0]==2017) & (dates[:,1]==8)& (dates[:,2]==26) & (dates[:,3]==6))
#To convert the tuple file into int64
plot_idx = which_plot[0][0]
plot_dat = var[plot_idx,:,:]
#TRIM DATA TO 2D IN ORDER TO MAKE A MAP
plot_date = dates[plot_idx,:]
plot_str = str(plot_date[1])+ '-' +str(plot_date[2]) + '-'+str(plot_date[0])+' '+str(plot_date[3])+'Z'
ax = plt.axes(projection=ccrs.PlateCarree())

ax.set_extent([np.min(ds.longitude),np.max(ds.longitude),np.min(ds.latitude),np.max(ds.latitude)], 
              crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE, linewidth=0.5,
               edgecolor='black', facecolor='None')

states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lines',
                                                scale='50m', facecolor='none')

tp_lev = np.arange(0,30,5) 
mslp_lev = np.arange(980,1020,4) 

p_plot = tp_2.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
                      levels=tp_lev, extend='both', add_colorbar=False,
                      add_labels=False)

plt.colorbar(p_plot, ax=ax, ticks=tp_lev, orientation='horizontal', pad=0.075)

msl_plot = mslp_plot.plot.contour(ax=ax, transform=ccrs.PlateCarree(),
                      levels=mslp_lev, colors='red', linewidths=0.5,
                      add_labels=True)

ax.add_feature(states_provinces, edgecolor='gray', zorder=2)

ax.clabel(msl_plot, fmt='%d', levels=mslp_lev)

[
     txt.set_bbox(dict(facecolor='yellow', edgecolor='yellow', pad=1))
     for txt in mslp_plot.labelTexts
]

gvutil.set_titles_and_labels(ax, maintitle='Total Precipitation/MSLP at Surface'+plot_str)
gvutil.add_lat_lon_ticklabels(ax)

ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

gvutil.add_major_minor_ticks(ax, x_minor_per_major=3, y_minor_per_major=5,
                             labelsize=5)