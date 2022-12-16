# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:27:17 2022

@author: sauda
"""
#Mapping precipitation during Hurricane Florence

#This module is dedicated to mapping better than ArcGIS.
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
#Cartopy produces maps and other geospatial data analyses.
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import geocat.viz.util as gvutil

#GeoCAT-viz is a collection of utility functions to make plotting with Matplotlib and Cartopy

file = 'ERA5_hourly.nc'
#folder = '/home/craig/Documents/'
folder = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/3/'
filename3 = folder+file

#open_dataset opens the file with read-only access. 
#If True, decode times encoded in the standard NetCDF datetime format into datetime objects.
ds = xr.open_dataset(filename3, decode_times=True)
print(ds)
#Extracting the precipitation variables
var = ds.tp.values
#var2 = ds.sst.values
#IVT = the vertically integrated horizontal water vapor transport 
time = ds.time.values
lats = ds.latitude.values
lons = ds.longitude.values

#extracting the day, month & year
#dt can be used to access the values of the series as datetimelike 
years = ds.time.dt.year.values
months = ds.time.dt.month.values
days = ds.time.dt.day.values
hours = ds.time.dt.hour.values
dates = np.transpose(np.stack((years,months,days,hours)))

#LOCATING EXACT POSITION OF A CERTAIN DAY'S DATA (14/09/2018) 
which_plot = np.where((dates[:,0]==2018) & (dates[:,1]==9)& (dates[:,2]==14) & (dates[:,3]==12))
#To convert the tuple file into int64
plot_idx = which_plot[0][0]
plot_data = var[plot_idx,:,:]
#TRIM DATA TO 2D IN ORDER TO MAKE A MAP
plot_date = dates[plot_idx,:]
plot_str = str(plot_date[1])+ '-' +str(plot_date[2]) + '-'+str(plot_date[0])+' '+str(plot_date[3])+'Z'

#below we plot
fig = plt.figure(figsize=(16,7),constrained_layout=True)
ax = plt.axes(projection=ccrs.PlateCarree())

# Cartopy packages for making maps within python
gvutil.set_axes_limits_and_ticks(ax,
                                 xlim=(min(lons),max(lons)),
                                 ylim=(min(lats),max(lats)))
                                 
gvutil.add_major_minor_ticks(ax,labelsize=10)

gvutil.add_lat_lon_ticklabels(ax)

colormap = "Blues"

space = 0.005
min_con = 0
max_con = 0.03

levels = np.arange(min_con,max_con,space)

kwargs = dict(levels=levels,add_colorbar=False,
              transform=ccrs.PlateCarree(), extend="both")

fillplot = ax.contourf(lons,lats,plot_data,cmap=colormap,**kwargs)

cb = fig.colorbar(fillplot, orientation='horizontal',
                  ticks=np.arange(min_con,max_con,space),
                  label='TP', shrink=0.75,extendrect=True,
                  extendfrac = 'auto', pad=0.04)
gvutil.set_titles_and_labels(ax,maintitle="Hurricane Florence TP"+plot_str)
ax.coastlines(linewidth=1.5,alpha=0.6)
states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m')
ax.add_feature(states, linewidth=0.75, edgecolor='k', facecolor='none', zorder=15)

out_folder = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/3/'
plot_out_str = 'hurr_Florence_'+plot_str
plt.savefig(out_folder+plot_str+'.png',dpi=1000)
plt.show()

#MEAN of all of your timesteps 
var_mean = np.mean(var,axis=0)
fig = plt.figure(figsize=(16,7),constrained_layout=True)
ax = plt.axes(projection=ccrs.PlateCarree())

# Cartopy packages for making maps within python
gvutil.set_axes_limits_and_ticks(ax,
                                 xlim=(min(lons),max(lons)),
                                 ylim=(min(lats),max(lats)))
                                 
gvutil.add_major_minor_ticks(ax,labelsize=10)

gvutil.add_lat_lon_ticklabels(ax)

colormap = "BuPu"

space = 0.0025
min_con = 0
max_con = 0.0275

levels = np.arange(min_con,max_con,space)

kwargs = dict(levels=levels,add_colorbar=False,
              transform=ccrs.PlateCarree(), extend="both")

fillplot = ax.contourf(lons,lats,var_mean,cmap=colormap,**kwargs)

cb = fig.colorbar(fillplot, orientation='horizontal',
                  ticks=np.arange(min_con,max_con,space),
                  label='TP', shrink=0.75,extendrect=True,
                  extendfrac = 'auto', pad=0.04)
gvutil.set_titles_and_labels(ax,maintitle="Hurricane Florence Mean TP")
ax.coastlines(linewidth=1.5,alpha=0.6)
states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m')
ax.add_feature(states, linewidth=0.75, edgecolor='k', facecolor='none', zorder=15)

out_folder = 'G:/My Drive/01_Spring 2022/Climate Data/Assignments/3/'
plot_out_str = 'hurr_Florence_Mean'
plt.savefig(out_folder+plot_out_str+'.png',dpi=1000)
plt.show()      