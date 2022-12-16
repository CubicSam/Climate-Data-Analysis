# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 09:53:20 2022

@author: sauda
"""
#what % of the data exceeded the cape value greater than 1000 in a particular location?
import xarray as xr 
import numpy as np
file = 'era5_cape.nc'
folder = 'G:/My Drive/01_Spring 2022/Climate Data//Lab/Lab 02/'
filename = folder+file


ds = xr. open_dataset(filename,decode_times=True)
print(ds)

#extracting the 3 variables
cape = ds.cape.values 

lat = ds.latitude.values
lon = ds.longitude.values
time = ds.time.values

#extracting the day, month & year
months = ds.time.dt.month.values
years = ds.time.dt.year.values
days = ds.time.dt.day.values
hours = ds.time.dt.hour.values
#Keeping all these 4 variables to the same dataset
dates = np.transpose(np.stack((years,months,days,hours)))

# Extracting grid points for Alabama state
var_gridlat_a = cape[:,((lat>30.4)&(lat<35.00)),:] 
var_grid_a = var_gridlat_a[:,:,((lon<360-85.02)&(lon>360-88.39))]
var_mean_grid_a = np.mean(var_grid_a)
var_mean_grid_a = np.mean(var_grid_a, axis = 0)

#Maximum Value
var_max_a= np.amax(var_mean_grid_a)

time_max_a_idx = np.where(var_mean_grid_a == np.amax(var_mean_grid_a)) 
time_max_a = dates[time_max_a_idx]

# make a new variable below to select the timestep with max CAPE
var_max_cape = var_grid_a[3,:,:]
cape_2860= np.where(var_max_cape>=2860)
cape_exceed_a = len(cape_2860[0])
#what % of the data exceeded the value in a particular location?
[d,e,f] = np.shape(var_grid_a)
#just need e and f here because those are the lat/lon variables (number of grid points)
total_num_grids_a = e*f
percent_exceed_a = (cape_exceed_a/total_num_grids_a)*100
