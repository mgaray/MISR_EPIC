# AERONET_V3_NetCDF_Daily_AOD15_Site_MERRA_Linear_Regression_Plot_01b.py
#
# This is a Python 3.8.9 code to read the NetCDF data generated by the
# AERONET_V3_DailyAverages_AOD15_Site_2_NetCDF_02.py code and the MERRA-2 data
# provided by Kyo Lee (JPL) and generate a linear regression plot with overall statistics.
#
# Creation Date: 2021-04-24
# Last Modified: 2021-04-24
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

import datetime as dt
import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mpdates
from netCDF4 import Dataset
import numpy as np
import os
import time

def main():  # Main code

# Time the code

    start_time = time.time()

# Set the paths

    aeropath = '/Users/mgaray/Desktop/CODING/PYTHON/PY38/APR21/AERONET/DATA/'
    datapath = '/Users/mgaray/Desktop/DATA/MERRA/'
    figpath = '/Users/mgaray/Desktop/CODING/PYTHON/PY38/APR21/MERRA/FIGS/'

# Set the AERONET product name

#    site_name = 'Banizoumbou'
#    site_name = 'Ilorin'
#    site_name = 'IER_Cinzana'
#    site_name = 'Ouagadougou'
#    site_name = 'LAMTO-STATION'
    site_name = 'Koforidua_ANUC'
    aero_name = 'AERONET_V3_DailyAverages_AOD15_'+site_name+'_2021_04_17_v02.nc'

# Set the MERRA-2 product name and epoch
    
    merra_name = "MERRA2_aerosol_variables_over_Africa_daily*"
    merra_epoch_str = "2000-02-29T00:00:00"
    merra_epoch = np.datetime64(merra_epoch_str)
  
# Set the base output file name  
  
    outname = site_name+'_MERRA_AOD_linear_plot_v01b.png'

### Find and read the correct AERONET NetCDF files

# Change to the data directory

    os.chdir(aeropath)

# Get the file list

    file_list = glob.glob(aero_name)
    
    num_files = len(file_list)

# Choose the correct file

    inputName = file_list[0]    

# Tell user location in process

    print('Reading: ',inputName)
    
# Open the NetCDF file

    rootgrp = Dataset(inputName, 'r', format='NETCDF4')
    
# Read the data

    lat = rootgrp.variables['Latitude'][:]
    lon = rootgrp.variables['Longitude'][:]
    year_raw = rootgrp.variables['Year'][:]
    month_raw = rootgrp.variables['Month'][:]
    day_raw = rootgrp.variables['Day'][:]
    aod_raw = rootgrp.variables['AOD_550nm'][:]
    
# Close the NetCDF file

    rootgrp.close() 
    
# Convert AERONET data to numpy arrays for analysis

    aero_year = np.array(year_raw)
    aero_month = np.array(month_raw)
    aero_day = np.array(day_raw)
        
    aero_aod = np.array(aod_raw)
    
# Extract the site latitude and longitude

    aero_lat = lat[0]
    aero_lon = lon[0]
   
# Get the number of unique years

    unique_years = np.unique(aero_year)
    num_years = len(aero_year)
    print("CHECK")
    print(num_years)
    print(unique_years)
    
# Generate the appropriate time string format

    time_raw = []

    num_found = len(aero_year)
    
    for i in np.arange(num_found):
    
        time_str = "{:}-".format(int(aero_year[i]))
        time_str = time_str+"{:02}-".format(int(aero_month[i]))
        time_str = time_str+"{:02}".format(int(aero_day[i]))
        base_time = time_str  

        time_str = base_time+"T12:00:00"  # Set to 12:00:00 

# Store the time string

        time_raw.append(time_str)

# Convert the date string to a numpy array

    aero_dates = np.array(time_raw,dtype='datetime64[D]') 
    
### Find and read the correct MERRA-2 NetCDF files

# Change to the data directory

    os.chdir(datapath)

# Get the file list

    file_list = glob.glob(merra_name)
    
    num_files = len(file_list)

# Choose the correct file

    inputName = file_list[0]    

# Tell user location in process

    print('Reading: ',inputName)
    
# Open the NetCDF file

    rootgrp = Dataset(inputName, 'r', format='NETCDF4')
    
# Read the data

    lat_raw = rootgrp.variables['lat'][:]
    lon_raw = rootgrp.variables['lon'][:]
    tim = rootgrp.variables['time'][:]
    aod_raw = rootgrp.variables['TOTEXTTAU'][:] 
    
# Close the NetCDF file

    rootgrp.close() 

# Convert MERRRA-2 data to numpy arrays for analysis

    lat = np.array(lat_raw)
    lon = np.array(lon_raw)

    day = np.array(tim)
    merra_aod = np.array(aod_raw)
    
# Convert the MERRA-2 time using the epoch

    merra_dates = merra_epoch + day.astype('timedelta64[D]')

### FIND THE GRID CELL CONTAINING THE AERONET SITE

    lat_ind = np.arange(len(lat))
    lon_ind = np.arange(len(lon))
    dlat = abs(lat-aero_lat)
    dlon = abs(lon-aero_lon)
    
    min_lat = np.amin(dlat)
    lat_match = (min_lat == dlat)
    lat_index = np.squeeze(lat_ind[lat_match])
    
    min_lon = np.amin(dlon)
    lon_match = (min_lon == dlon)
    lon_index = np.squeeze(lon_ind[lon_match])
    
### Extract the data from this grid cell

    merra_aod_grid = merra_aod[:,lat_index,lon_index]
    
### MATCH TO THE AERONET DATA

    num_aero = len(aero_aod)
    aero_good = np.zeros(num_aero)
    merra_good = np.zeros(num_aero)
    dates_good = np.zeros_like(aero_dates)
    
    for i in np.arange(num_aero):
    
        this_aod = aero_aod[i]
        this_date = aero_dates[i]
        
        delta = np.abs(merra_dates - this_date)
        min_delta = np.amin(delta)
        delta_int = int(min_delta/np.timedelta64(1,'D'))
        if(delta_int > 1):
            continue
            
        keep = (delta == min_delta)
        
        dates_good[i] = this_date
        aero_good[i] = this_aod
        merra_good[i] = np.squeeze(merra_aod_grid[keep])
        
### PLOT THE RESULTS

# Set the minimum and maximum AOD for plotting

    aod_plot_min = 0.00
    aod_plot_max = 4.5
    
    aod_plot_ticks = 8 # Usually 1 more than you think
    aod_plot_step = 0.25
    
    max_val = aod_plot_max
    
# Set the plot area (using the concise format)

    fig, ax = plt.subplots(figsize=(6.0,6.0), dpi=120)

# Extract the valid data to plot

    good = ((aero_good > 0) & (merra_good > 0))
    ref_aod = aero_good[good]  # Set AERONET as the reference dataset
    test_aod = merra_good[good]  # Set MERRA as the comparison dataset

# Plot the data

    ax.scatter(ref_aod,test_aod,marker='o',color='black',s=5)

# Plot the one-to-one line

    ax.plot([aod_plot_min,aod_plot_max], [aod_plot_min,aod_plot_max], color="0.25", lw=1)

# Plot the envelopes

    dummy_aod = np.logspace(-4,1,num=100)
    up1_aod = 1.20*dummy_aod
    up2_aod = dummy_aod+0.05
    upper_aod = np.maximum(up1_aod,up2_aod)

    lo1_aod = 0.80*dummy_aod
    lo2_aod = dummy_aod-0.05
    lower_aod = np.minimum(lo1_aod,lo2_aod)

    ax.plot(dummy_aod,lower_aod,color="0.75", lw=1)
    ax.plot(dummy_aod,upper_aod,color="0.75", lw=1)

# Axes and labels

    title_text = site_name
    ax.set_title(title_text,size=20)

    ax.set_xlim(aod_plot_min,aod_plot_max)
    ax.set_xticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5])
    ax.set_xticklabels(['0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5'],size=15)
    ax.set_xlabel('AERONET AOD',size=18)
    
    ax.set_ylim(aod_plot_min,aod_plot_max)
    ax.set_yticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5])
    ax.set_yticklabels(['0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5'],size=15)
    ax.set_ylabel('MERRA-2 AOD',size=18)
    
    ax.grid(True)
    
# Additional Text

    x_pos = 0.25
    y_pos1 = 3.0
    y_pos2 = y_pos1 + (aod_plot_max/25.)
    y_pos3 = y_pos2 + (aod_plot_max/25.)
    y_pos4 = y_pos3 + (aod_plot_max/25.)
    y_pos5 = y_pos4 + (aod_plot_max/25.)
    y_pos6 = y_pos5 + (aod_plot_max/25.)
    
    count = len(test_aod)
    out_text = 'N = '+str(count)
    plt.text(x_pos,y_pos5,out_text,fontsize=15) # Count

    temp = np.corrcoef(ref_aod,test_aod)
    be_r = temp[0,1]
    out_text = 'r = '+"{0:.4f}".format(be_r)
    plt.text(x_pos,y_pos4,out_text,fontsize=15) # Correlation coefficient

    rmse = np.sqrt(((test_aod - ref_aod) ** 2).mean())
    out_text = 'RMSE = '+"{0:.4f}".format(rmse)
    plt.text(x_pos,y_pos3,out_text,fontsize=15) # Root mean squared error

    diff = test_aod - ref_aod
    bias = np.mean(diff)
    out_text = 'Bias = '+"{0:.4f}".format(bias)
    plt.text(x_pos,y_pos2,out_text,fontsize=15) # Bias

    offset = np.ones_like(ref_aod)*0.05
    inner = np.absolute(diff) < np.maximum(offset,ref_aod*0.2)
    in_frac = (np.sum(inner)/(1.0*count))*100.0
    out_text = 'Percent In = '+"{0:.2f}".format(in_frac)
    plt.text(x_pos,y_pos1,out_text,fontsize=15) # Percent in envelope

# Tight layout    

    plt.tight_layout()

# Save the data
    
    os.chdir(figpath)
    plt.savefig(outname,dpi=300)
    print("Writing: "+outname)

    plt.show()  
        
# Print the time

    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION

def Haversine_Distance(lat1,lon1,lat_arr,lon_arr):
### Distance Calculation Based on the Haversine Formula
# Creation Date: 2015-05-13
# Last Modified: 2015-05-13
# By Michael J. Garay
# Michael.J.Garay@jpl.nasa.gov
#
# Note: This follows the formula from http://williams.best.vwh.net/avform.htm#Dist
# But see a discussion on the Earth-radius at 
# http://www.cs.nyu.edu/visual/home/proj/tiger/gisfaq.html
#
# Input: lat1 = First latitude, single element(degrees)
#        lon1 = First longitude, single element
#        lat2 = Second latitude, array of values
#        lon2 = Second longitude, array of values
#
# Output: Returns an array of distances (km)

# Convert lat/lon to radians

    rat1 = lat1*np.pi/180.0
    ron1 = lon1*np.pi/180.0
    
    rat2 = lat_arr*np.pi/180.0
    ron2 = lon_arr*np.pi/180.0

# Calculate the distance using the Haversine Formula

    d = 2.0*np.arcsin(np.sqrt((np.sin((rat2-rat1)/2))**2 +
      np.cos(rat2)*np.cos(rat1)*(np.sin((ron2-ron1)/2))**2))

# Convert to kilometers
    
    dist = 6371.0 * d
    
    return dist

### END Haversine_Distance

if __name__ == '__main__':
    main()    