# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 06:51:01 2022

@author: sauda
"""
#Calculating seasonal statistics with 2-D Data in Python

import numpy as np
#mathematical operations
import pandas as pd
#data analysis
import matplotlib.pyplot as plt
#data visualization
import seaborn as sns
#better data visualization

data=pd.read_csv("G:/My Drive/01_Spring 2022/Climate Data/Assignments/1/Ashville_Tmax&Prec.csv")
var = np.array(data.TMAX)
date = pd.DatetimeIndex(data.DATE)
month = date.month
winter = np.where((month == 12)|(month == 1)|(month == 2))
summer = np.where((month == 6)|(month == 7) | (month == 8))
spring = np.where ((month == 3) | (month == 4) | (month == 5))
fall = np.where ((month == 9) | (month == 10) | (month == 11))

#Extracting the Tmax for the seasons
var = np.array(data.TMAX)
var_max = np.max(var)
winter_data = var[winter]
summer_data = var[summer]
spring_data = var[spring]
fall_data = var[fall]
#Summer stat of TMax
summer_mean = np.mean(summer_data)
summer_median = np.median(summer_data)
summer_std = np.std(summer_data)
summer_max = np.max(summer_data)
#Winter stat of TMax
winter_mean = np.mean(winter_data)
winter_median = np.median(winter_data)
winter_std = np.std(winter_data)
winter_max = np.max(winter_data)
#Fall stat of TMax
fall_mean = np.mean(fall_data)
fall_median = np.median(fall_data)
fall_std = np.std(fall_data)
fall_max = np.max(fall_data)
#Spring stat for TMax
spring_mean = np.mean(spring_data)
spring_median = np.median(spring_data)
spring_std = np.std(spring_data)
spring_max = np.max(spring_data)

#Plotting
#TMax for Winter (Simple Plotting)
out_folder = "G:/My Drive/01_Spring 2022/Climate Data/Assignments/1"

#simple line plot
plt.rcParams["font.family"] = ' Times New Roman'
plt.plot(winter_data,color='m')
plt.ylabel('TMax')
plt.xlabel('Day')
plt.title('Asheville Max Temperature in Winter (2017-2021)')
plt.savefig(out_folder+"ashv_maxt_lineplot.png")

#Extracting the Precipitation for seasons
varp = np.array(data.PRCP)
varp_max = np.max(varp)
winter_pre = varp[winter]
summer_pre = varp[summer]
spring_pre = varp[spring]
fall_pre = varp[fall]
#Summer stat of Precip
summer_mean_p = np.mean(summer_pre)
summer_median_p = np.median(summer_pre)
summer_std_p = np.std(summer_pre)
summer_max_p = np.max(summer_pre)
#Winter stat of Precip
winter_mean_p = np.mean(winter_pre)
winter_median_p = np.median(winter_pre)
winter_std_p = np.std(winter_pre)
winter_max_p = np.max(winter_pre)
#Fall stat of Precip
fall_mean_p = np.mean(fall_pre)
fall_median_p = np.median(fall_pre)
fall_std_p = np.std(fall_pre)
fall_max_p = np.max(fall_pre)
#Spring stat for Precip
spring_mean_p = np.mean(spring_pre)
spring_median_p = np.median(spring_pre)
spring_std_p = np.std(spring_pre)
spring_max_p = np.max(spring_pre)

#For histogram
plt.rcParams["font.family"]= 'Times New Roman'
sns.set_palette('muted')
sns.set_style('dark')
kwargs = dict(alpha=0.9,linewidth=2)
min_con = 0
max_con=  1
num_bins = 10
plt.hist(fall_pre, **kwargs,color='r',label='Fall')
plt.hist(spring_pre, **kwargs,color='b',label='Spring')
plt.title('Fall & Spring Precipitation in Asheville (2017-2021)')
plt.ylabel('Frequency')
plt.xlabel('Precipitation')
plt.legend()
plt.savefig(out_folder+"ashv_prec_spr_fall.jpg")
plt.close()
