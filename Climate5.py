# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:48:05 2022

@author: sauda
"""
#Creating contour maps showing the mean monthly values for JANUARY and for JULY

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
from netCDF4 import Dataset

from geocat.viz import util as gvutil
from geocat.viz import cmaps as gvcmaps

file_mid = 'mid_term_data.nc'
folder_mid = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/Mid/'
filename_mid = folder_mid+file_mid
ds_mid = xr.open_dataset(filename_mid, decode_times=True)
print(ds_mid)

#extracting the variables
d2m = ds_mid.d2m.values 
ds_mid.variables

time = ds_mid.time.values
lat = ds_mid.latitude.values
lon = ds_mid.longitude.values

# Extracting grid points for Texas state
var_gridlat_t = d2m[:,((lat>25.59)&(lat<36.53)),:] 
var_grid_t = var_gridlat_t[:,:,((lon<-92.99)&(lon>-107.6082))]
var_mean_grid_t = np.mean(var_grid_t, axis = 0)

#Extracting maximum & minimum value of d2m
d2m_max_t= np.amax(var_mean_grid_t)
d2m_min_t= np.amin(var_mean_grid_t)

#extracting the day, month & year
months = ds_mid.time.dt.month.values
years = ds_mid.time.dt.year.values
days = ds_mid.time.dt.day.values
hours = ds_mid.time.dt.hour.values
#Keeping all these 4 variables to the same dataset
dates_all = np.transpose(np.stack((years,months,days,hours)))

#Time for the maximum d2m
time_max_t_idx = np.where(d2m_max_t)
time_max_t = dates_all[time_max_t_idx]

unique_months = np.unique(months)
[a,b,c] = np.shape(d2m)
monthly_daily_mean = np.empty((len(unique_months),b,c))

range(len(unique_months))
#range(0, 12) = 0 to 11
#finding months 
for i in range(len(unique_months)):
    curr_month = unique_months[i]
    find_months = np.where(months==curr_month)[0]
    #taking the mean only for Jan exclusing othe rmonths
    #monthly_daily_mean = np.mean(prcp_all[find_months,:,:],axis=0)- the loop overwrites the values everytime we run so see below
    monthly_daily_mean[i,:,:] = np.mean(d2m[find_months,:,:],axis=0) 
    #axis = 0 means taking the mean for Jan

#jan = np.where((months == 1))
#july = np.where((months == 7))
#jan_data = d2m[jan]
#july_data = d2m[july]


#mean monthly values for JANUARY and JULY

plot_jan = monthly_daily_mean[0,:]
plot_july = monthly_daily_mean[6,:]

#below we plot
#For Jan
fig = plt.figure(figsize=(16,7),constrained_layout=True)
ax = plt.axes(projection=ccrs.PlateCarree())

# Cartopy packages for making maps within python
gvutil.set_axes_limits_and_ticks(ax,
                                 xlim=(min(lon),max(lon)),
                                 ylim=(min(lat),max(lat)))
                                 
gvutil.add_major_minor_ticks(ax,labelsize=10)

gvutil.add_lat_lon_ticklabels(ax)

colormap = "OrRd"

space = 2
min_con = 250
max_con = 280

levels = np.arange(min_con,max_con,space)

kwargs = dict(levels=levels,add_colorbar=False,
              transform=ccrs.PlateCarree(), extend="both")

fillplot = ax.contourf(lon,lat,plot_jan,cmap=colormap,**kwargs)
cb = fig.colorbar(fillplot, orientation='horizontal',
                  ticks=np.arange(min_con,max_con,space),
                  label='Temperature (K)', shrink=0.75,extendrect=True,
                  extendfrac = 'auto', pad=0.04)
gvutil.set_titles_and_labels(ax,maintitle="2 metre dewpoint temperature in January")
ax.coastlines(linewidth=1.5,alpha=0.6)
states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m')
ax.add_feature(states, linewidth=0.75, edgecolor='k', facecolor='none', zorder=15)

out_folder = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/Mid/'
plot_out_str = 'D2M_Jan'
plt.savefig(out_folder+plot_out_str+'.png',dpi=1000)

plt.show()  

#For July
fig = plt.figure(figsize=(16,7),constrained_layout=True)
ax = plt.axes(projection=ccrs.PlateCarree())
np.amin(plot_july)
np.amax(plot_july)
# Cartopy packages for making maps within python
gvutil.set_axes_limits_and_ticks(ax,
                                 xlim=(min(lon),max(lon)),
                                 ylim=(min(lat),max(lat)))
                                 
gvutil.add_major_minor_ticks(ax,labelsize=10)

gvutil.add_lat_lon_ticklabels(ax)

colormap = "OrRd"

space = 2
min_con = 260
max_con = 300

levels = np.arange(min_con,max_con,space)

kwargs = dict(levels=levels,add_colorbar=False,
              transform=ccrs.PlateCarree(), extend="both")

fillplot = ax.contourf(lon,lat,plot_july,cmap=colormap,**kwargs)
cb = fig.colorbar(fillplot, orientation='horizontal',
                  ticks=np.arange(min_con,max_con,space),
                  label='Temperature (K)', shrink=0.75,extendrect=True,
                  extendfrac = 'auto', pad=0.04)
gvutil.set_titles_and_labels(ax,maintitle="2 metre dewpoint temperature in July")
ax.coastlines(linewidth=1.5,alpha=0.6)
states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m')
ax.add_feature(states, linewidth=0.75, edgecolor='k', facecolor='none', zorder=15)


out_folder = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/Mid/'
plot_out_str = 'D2M_July'
plt.savefig(out_folder+plot_out_str+'.png',dpi=1000)
plt.show()  

