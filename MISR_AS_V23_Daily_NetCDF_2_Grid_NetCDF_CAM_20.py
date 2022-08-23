# MISR_AS_V23_Daily_NetCDF_2_Grid_NetCDF_CAM_20.py
#
# This is a Python 3.8.11 code to read a set of days of MISR V23 aerosol files, 
# extract the data globally, grid the data onto the CAM model grid, and store the data 
# in a single NetCDF file for each day.
#
# NOTE: The date iteration comes from 
# https://stackoverflow.com/questions/1060279/iterating-through-a-range-of-dates-in-python
#
# The rolling variance comes from B. P. Welford (Technometrics, 1962) and is referenced 
# here:
# https://www.johndcook.com/blog/standard_deviation/
#
# Creation Date: 2021-09-06
# Last Modified: 2021-09-06
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

from datetime import date, timedelta
import glob
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os
import time

def main():  # Main code

# Print a blank line

    print()

# Set the overall timer

    all_start_time = time.time()  

### START USER INPUTS ###
           
# Set the paths

    basepath = '/Volumes/NEW_MISR_AEROSOL/DATA/MISR_AEROSOL/2014/'
    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY38/SEP21/CAM/DATA/'

# Set the start and end dates
# NOTE: Format is start_date = date(YYYY, MM, DD)

    start_date = date(2014, 2, 1)
    end_date = date(2015, 1, 1)

# Set the grid resolution
# NOTE: The CAM grid is specified by the number of cells

    num_lon = 288.0
    num_lat = 192.0
    
    lon_min = 0.0
    lon_max = 360.0  # Degrees east
    
    lat_min = -90.0
    lat_max = 90.0  # Degrees north

### END USER INPUTS ###  

# Set the  grids for the appropriate resolution

    delt_lon = lon_max - lon_min
    delt_lat = lat_max - lat_min
  
    lon_grid = delt_lon/num_lon
    lat_grid = delt_lat/(num_lat-1.)  # Eliminate one to center on poles

# Calculate the center latitude/longitude
# NOTE: These correspond to the CAM output received from Yuan Wang
    
    glon = np.arange(num_lon)*lon_grid+lon_min
    glat = np.arange(num_lat)*lat_grid+lat_min
    
    lon_off = lon_grid/2.0  # Longitude offset
    lat_off = lat_grid/2.0  # Latitude offset
    
    num_lon = len(glon)
    num_lat = len(glat)
    
# Create 2-D lat/lon arrays from the 1-D arrays

    lat2D, lon2D = np.meshgrid(glat,glon)

# Calculate the temporal arrays

    num_time = 24  # Hourly
    gtime = np.arange(num_time)

### LARGE LOOP OVER THE DATES IN THE SELECTION

    for single_date in daterange(start_date,end_date):

        temp = single_date.strftime("%Y-%m-%d")
        start_date_str = temp+"T00:00:00"
        print(start_date_str)
        
# Parse the start date string

        hold = start_date_str.split('T')
        dhold = hold[0].split('-')
        thold = hold[1].split(':')
    
        start_year_str = dhold[0]
        start_month_str = dhold[1]
        start_day_str = dhold[2]
        start_YYYMMDD_str = hold[0]
    
        start_hour_str = thold[0]
        start_minute_str = thold[1]
        start_second_str = thold[2]
        start_HHMMSS_str = hold[1]

# Generate the output file

        outfile = "MISR_AS_V23_CAM_Grid_"+start_YYYMMDD_str+"_288x192_v20.nc"

# Generate the arrays to store the output data
# NOTE: This needs to be in the loop to zero out the arrays

        aod_best_tot = np.zeros((num_lon,num_lat,num_time))
        log_aod_best_tot = np.zeros((num_lon,num_lat,num_time))
        aod_best_tot_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_best_mk_tot = np.zeros((num_lon,num_lat,num_time))
        aod_best_sk_tot = np.zeros((num_lon,num_lat,num_time))
    
        aod_all_tot = np.zeros((num_lon,num_lat,num_time))
        log_aod_all_tot = np.zeros((num_lon,num_lat,num_time))
        aod_all_tot_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_all_mk_tot = np.zeros((num_lon,num_lat,num_time))
        aod_all_sk_tot = np.zeros((num_lon,num_lat,num_time))
    
        aod_best_land = np.zeros((num_lon,num_lat,num_time))
        log_aod_best_land = np.zeros((num_lon,num_lat,num_time))
        aod_best_land_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_best_mk_land = np.zeros((num_lon,num_lat,num_time))
        aod_best_sk_land = np.zeros((num_lon,num_lat,num_time))
    
        aod_all_land = np.zeros((num_lon,num_lat,num_time))
        log_aod_all_land = np.zeros((num_lon,num_lat,num_time))
        aod_all_land_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_all_mk_land = np.zeros((num_lon,num_lat,num_time))
        aod_all_sk_land = np.zeros((num_lon,num_lat,num_time))
    
        aod_best_water = np.zeros((num_lon,num_lat,num_time))
        log_aod_best_water = np.zeros((num_lon,num_lat,num_time))
        aod_best_water_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_best_mk_water = np.zeros((num_lon,num_lat,num_time))
        aod_best_sk_water = np.zeros((num_lon,num_lat,num_time))
    
        aod_all_water = np.zeros((num_lon,num_lat,num_time))
        log_aod_all_water = np.zeros((num_lon,num_lat,num_time))
        aod_all_water_cnt = np.zeros((num_lon,num_lat,num_time))
        aod_all_mk_water = np.zeros((num_lon,num_lat,num_time))
        aod_all_sk_water = np.zeros((num_lon,num_lat,num_time))

## Identify the orbits associated with the requested date
    
        start_date = np.datetime64(start_date_str)  # Convert the time string to datetime64
        first_orbit = misr_datetime64_to_orbit(start_date) 
    
        end_date = start_date+np.timedelta64(1,'D')
        last_orbit = misr_datetime64_to_orbit(end_date) 

# Get the number of orbits and generate a list containing the orbit numbers
# NOTE: Because of the way orbits map (poorly) to dates, we include orbits
#       prior to the first and after the last
    
        orbit_list = np.arange(first_orbit-2,last_orbit+2,1)  # Inclusive

### FIND THE AEROSOL DATA

# Change to the correct directory

        os.chdir(basepath)

# Loop over the orbits

        for orbit in orbit_list:
    
            orbit_str = str(orbit)

# Search for the individual MISR files

            search_str = 'MISR_AM1_AS_AEROSOL_*'+orbit_str+'*.nc'
            dum_list = glob.glob(search_str)
        
            if(len(dum_list) < 1):
                print("***MISSING ORBIT***")
                print(orbit_str)
                continue
        
# Select the file

            misr_file = dum_list[0]
            
# Tell user location in process

            print("Reading: "+misr_file)
        
# Open the NetCDF file

            rootgrp = Dataset(misr_file, 'r', format='NETCDF4')

# Choose the appropriate subgroup

            subgrp = rootgrp.groups['4.4_KM_PRODUCTS']

# Assign the variables to (masked) arrays

            lat_masked = subgrp.variables['Latitude'][:]
            lon_masked = subgrp.variables['Longitude'][:]
        
            year_masked = subgrp.variables['Year'][:]
            month_masked = subgrp.variables['Month'][:]
            day_masked = subgrp.variables['Day'][:]
            
            hour_masked = subgrp.variables['Hour'][:]
            minute_masked = subgrp.variables['Minute'][:]
               
            aod_masked = subgrp.variables['Aerosol_Optical_Depth'][:]
            
            lwr_masked = subgrp.variables['Land_Water_Retrieval_Type'][:]
            
# Choose the AUXILIARY subsubgroup

            subsubgrp = subgrp.groups['AUXILIARY']
        
            aod_raw_masked = subsubgrp.variables['Aerosol_Optical_Depth_Raw'][:]
            lwr_raw_masked = subsubgrp.variables['Land_Water_Retrieval_Type_Raw'][:]
                             
# Close the NetCDF file

            rootgrp.close()
        
# Convert the masked arrays to numpy arrays that can be properly analyzed

            lat = ma.filled(lat_masked,fill_value=-9999.0)
            lon = ma.filled(lon_masked,fill_value=-9999.0)
        
            year = ma.filled(year_masked,fill_value=-9999.0)
            month = ma.filled(month_masked,fill_value=-9999.0)
            day = ma.filled(day_masked,fill_value=-9999.0)
            
            hour = ma.filled(hour_masked,fill_value=-9999.0)
            minute = ma.filled(minute_masked,fill_value=-9999.0)
    
            aod = ma.filled(aod_masked,fill_value=-9999.0)
            lwr = ma.filled(lwr_masked,fill_value=253)
            
            aod_raw = ma.filled(aod_raw_masked,fill_value=-9999.0)
            lwr_raw = ma.filled(lwr_raw_masked,fill_value=253)

# Extract the locations with valid raw retrievals
# NOTE: There are far fewer locations with valid AOD (raw) retrievals than locations seen
#       by the instrument

            raw = (aod_raw > 0.0)
        
            lat_all = lat[raw]
            lon_all = lon[raw]
        
            year_all = year[raw]
            month_all = month[raw]
            day_all = day[raw]
        
            hour_all = hour[raw]
            minute_all = minute[raw]
        
            aod_best = aod[raw]
            lwr_best = lwr[raw]
        
            aod_all = aod_raw[raw]
            lwr_all = lwr_raw[raw]
        
            num_all = len(aod_all)
        
            if(num_all < 1):
                print("NO VALID DATA")
                continue

# LOOP THROUGH THE DATA

            for i in range(num_all):
        
# Generate a time string

                this_year = year_all[i]
                this_month = month_all[i]
                this_day = day_all[i]
            
                this_hour = hour_all[i]
                this_minute = minute_all[i]
            
                time_str = "{:}-".format(int(this_year))
                time_str = time_str+"{:02}-".format(int(this_month))
                time_str = time_str+"{:02}".format(int(this_day))
                base_time = time_str  

                time_str = base_time+"T"+"{:02}:".format(int(this_hour))
                time_str = time_str+"{:02}:".format(int(this_minute))  
                time_str = time_str+"00"  # No seconds
            
# Convert the time string to the correct time format

                this_date = np.datetime64(time_str)
               
# Parse the time to find relevant temporal matches

                delta_hour = (this_date - start_date)/np.timedelta64(1,'h')
            
# Convert the time difference to a grid index
# NOTE: We add 0.5 hours to center on the hour itself

                time_ind = int(np.floor(delta_hour+0.5))  # Convert output to an INT
            
                if(time_ind < 0):
                    continue
                
                if(time_ind > 23):
                    continue

# Calculate the other grid values
     
                this_lat = lat_all[i]
                this_lon = lon_all[i]
                
# Convert lat/lon to a grid index

                lat_ind = int((this_lat - lat_min + lat_off)/lat_grid)
                lon_ind = int((this_lon - lon_min + lon_off)/lon_grid)

## RAW VALUES - ALL
                
# Extract the RAW values

                this_aod = aod_all[i]
                                
# Store the RAW values

                temp = aod_all_tot[lon_ind,lat_ind,time_ind]
                aod_all_tot[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                temp = log_aod_all_tot[lon_ind,lat_ind,time_ind]
                log_aod_all_tot[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                temp = aod_all_tot_cnt[lon_ind,lat_ind,time_ind]
                aod_all_tot_cnt[lon_ind,lat_ind,time_ind] = temp+1
                
# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                k = temp  # Current count in bin
                if(k < 1):
                    aod_all_mk_tot[lon_ind,lat_ind,time_ind] = this_aod
                    aod_all_sk_tot[lon_ind,lat_ind,time_ind] = 0.0
                else:
                    old_m = aod_all_mk_tot[lon_ind,lat_ind,time_ind]
                    old_s = aod_all_sk_tot[lon_ind,lat_ind,time_ind]
                    new_m = old_m + (this_aod-old_m)/k
                    new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                    aod_all_mk_tot[lon_ind,lat_ind,time_ind] = new_m
                    aod_all_sk_tot[lon_ind,lat_ind,time_ind] = new_s

## GET THE TYPE OF RETRIEVAL
## 0 - DARK WATER
## 1 - HET SURF

                this_type = lwr_all[i]

# Store the RAW DARK WATER values

                if(this_type == 0):

                    temp = aod_all_water[lon_ind,lat_ind,time_ind]
                    aod_all_water[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                    temp = log_aod_all_water[lon_ind,lat_ind,time_ind]
                    log_aod_all_water[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                    temp = aod_all_water_cnt[lon_ind,lat_ind,time_ind]
                    aod_all_water_cnt[lon_ind,lat_ind,time_ind] = temp+1
                
# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                    k = temp  # Current count in bin
                    if(k < 1):
                        aod_all_mk_water[lon_ind,lat_ind,time_ind] = this_aod
                        aod_all_sk_water[lon_ind,lat_ind,time_ind] = 0.0
                    else:
                        old_m = aod_all_mk_water[lon_ind,lat_ind,time_ind]
                        old_s = aod_all_sk_water[lon_ind,lat_ind,time_ind]
                        new_m = old_m + (this_aod-old_m)/k
                        new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                        aod_all_mk_water[lon_ind,lat_ind,time_ind] = new_m
                        aod_all_sk_water[lon_ind,lat_ind,time_ind] = new_s    
                        
# Store the RAW HET SURF values

                if(this_type == 1):

                    temp = aod_all_land[lon_ind,lat_ind,time_ind]
                    aod_all_land[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                    temp = log_aod_all_land[lon_ind,lat_ind,time_ind]
                    log_aod_all_land[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                    temp = aod_all_land_cnt[lon_ind,lat_ind,time_ind]
                    aod_all_land_cnt[lon_ind,lat_ind,time_ind] = temp+1
                
# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                    k = temp  # Current count in bin
                    if(k < 1):
                        aod_all_mk_land[lon_ind,lat_ind,time_ind] = this_aod
                        aod_all_sk_land[lon_ind,lat_ind,time_ind] = 0.0
                    else:
                        old_m = aod_all_mk_land[lon_ind,lat_ind,time_ind]
                        old_s = aod_all_sk_land[lon_ind,lat_ind,time_ind]
                        new_m = old_m + (this_aod-old_m)/k
                        new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                        aod_all_mk_land[lon_ind,lat_ind,time_ind] = new_m
                        aod_all_sk_land[lon_ind,lat_ind,time_ind] = new_s                

## USER VALUES - ALL
            
# Extract the USER values

                this_aod = aod_best[i]
                if(this_aod <= 0.):
                    continue
                                
# Store the USER values

                temp = aod_best_tot[lon_ind,lat_ind,time_ind]
                aod_best_tot[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                temp = log_aod_best_tot[lon_ind,lat_ind,time_ind]
                log_aod_best_tot[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                temp = aod_best_tot_cnt[lon_ind,lat_ind,time_ind]
                aod_best_tot_cnt[lon_ind,lat_ind,time_ind] = temp+1

# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                k = temp
                if(k < 1):
                    aod_best_mk_tot[lon_ind,lat_ind,time_ind] = this_aod
                    aod_best_sk_tot[lon_ind,lat_ind,time_ind] = 0.0
                else:
                    old_m = aod_best_mk_tot[lon_ind,lat_ind,time_ind]
                    old_s = aod_best_sk_tot[lon_ind,lat_ind,time_ind]
                    new_m = old_m + (this_aod-old_m)/k
                    new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                    aod_best_mk_tot[lon_ind,lat_ind,time_ind] = new_m
                    aod_best_sk_tot[lon_ind,lat_ind,time_ind] = new_s

## GET THE TYPE OF RETRIEVAL
## 0 - DARK WATER
## 1 - HET SURF

                this_type = lwr_best[i]

# Store the RAW DARK WATER values

                if(this_type == 0):

                    temp = aod_best_water[lon_ind,lat_ind,time_ind]
                    aod_best_water[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                    temp = log_aod_best_water[lon_ind,lat_ind,time_ind]
                    log_aod_best_water[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                    temp = aod_best_water_cnt[lon_ind,lat_ind,time_ind]
                    aod_best_water_cnt[lon_ind,lat_ind,time_ind] = temp+1
                
# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                    k = temp  # Current count in bin
                    if(k < 1):
                        aod_best_mk_water[lon_ind,lat_ind,time_ind] = this_aod
                        aod_best_sk_water[lon_ind,lat_ind,time_ind] = 0.0
                    else:
                        old_m = aod_best_mk_water[lon_ind,lat_ind,time_ind]
                        old_s = aod_best_sk_water[lon_ind,lat_ind,time_ind]
                        new_m = old_m + (this_aod-old_m)/k
                        new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                        aod_best_mk_water[lon_ind,lat_ind,time_ind] = new_m
                        aod_best_sk_water[lon_ind,lat_ind,time_ind] = new_s    
                        
# Store the RAW HET SURF values

                if(this_type == 1):

                    temp = aod_best_land[lon_ind,lat_ind,time_ind]
                    aod_best_land[lon_ind,lat_ind,time_ind] = temp+this_aod
            
                    temp = log_aod_best_land[lon_ind,lat_ind,time_ind]
                    log_aod_best_land[lon_ind,lat_ind,time_ind] = temp+np.log10(this_aod)
                
                    temp = aod_best_land_cnt[lon_ind,lat_ind,time_ind]
                    aod_best_land_cnt[lon_ind,lat_ind,time_ind] = temp+1
                
# Calculate the running variance
# NOTE: This uses the approach from 
# https://www.johndcook.com/blog/standard_deviation/

                    k = temp  # Current count in bin
                    if(k < 1):
                        aod_best_mk_land[lon_ind,lat_ind,time_ind] = this_aod
                        aod_best_sk_land[lon_ind,lat_ind,time_ind] = 0.0
                    else:
                        old_m = aod_best_mk_land[lon_ind,lat_ind,time_ind]
                        old_s = aod_best_sk_land[lon_ind,lat_ind,time_ind]
                        new_m = old_m + (this_aod-old_m)/k
                        new_s = old_s + (this_aod-old_m)*(this_aod-new_m)
                        aod_best_mk_land[lon_ind,lat_ind,time_ind] = new_m
                        aod_best_sk_land[lon_ind,lat_ind,time_ind] = new_s                
                                                              
### SAVE THE DATA TO A NETCDF FILE

# Create 2-D arrays containing NaNs

        mean_aod_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_nan = np.full((num_lon,num_lat,num_time),np.nan)

        mean_aod_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        
        mean_aod_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_land_nan = np.full((num_lon,num_lat,num_time),np.nan)

        mean_aod_raw_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_raw_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_raw_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_raw_nan = np.full((num_lon,num_lat,num_time),np.nan)
        
        mean_aod_raw_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_raw_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_raw_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_raw_water_nan = np.full((num_lon,num_lat,num_time),np.nan)
        
        mean_aod_raw_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        gmean_aod_raw_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        count_aod_raw_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
        stdev_aod_raw_land_nan = np.full((num_lon,num_lat,num_time),np.nan)
    
# Change to the data directory

        os.chdir(datapath)
    
# Open the file
        
        n_id = Dataset(outfile,'w')
    
# Define the dimensions

        lon = n_id.createDimension('lon',num_lon) # Longitude
        lat = n_id.createDimension('lat',num_lat) # Latitude
        hour = n_id.createDimension('hour',num_time) # Time

# Define the output variables (mapping)

        out01 = n_id.createVariable('lon','f4',('lon'),zlib=True)
        out01.units = 'degrees_east'
        out02 = n_id.createVariable('lat','f4',('lat'),zlib=True)
        out02.units = 'degrees_north'
        out03 = n_id.createVariable('hour','f4',('hour'),zlib=True)
        out03.units = 'hour'

# Define the output variables

        out04 = n_id.createVariable('Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out05 = n_id.createVariable('Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out06 = n_id.createVariable('Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out07 = n_id.createVariable('Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)

        out08 = n_id.createVariable('Water_Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out09 = n_id.createVariable('Water_Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out10 = n_id.createVariable('Water_Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out11 = n_id.createVariable('Water_Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)

        out12 = n_id.createVariable('Land_Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out13 = n_id.createVariable('Land_Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out14 = n_id.createVariable('Land_Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out15 = n_id.createVariable('Land_Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)

        out16 = n_id.createVariable('Raw_Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out17 = n_id.createVariable('Raw_Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out18 = n_id.createVariable('Raw_Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out19 = n_id.createVariable('Raw_Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)

        out20 = n_id.createVariable('Raw_Water_Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out21 = n_id.createVariable('Raw_Water_Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out22 = n_id.createVariable('Raw_Water_Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out23 = n_id.createVariable('Raw_Water_Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)
        
        out24 = n_id.createVariable('Raw_Land_Aerosol_Optical_Depth_Mean','f4',('lon','lat','hour'),zlib=True)
        out25 = n_id.createVariable('Raw_Land_Aerosol_Optical_Depth_Geometric_Mean','f4',('lon','lat','hour'),zlib=True)
        out26 = n_id.createVariable('Raw_Land_Aerosol_Optical_Depth_Count','f4',('lon','lat','hour'),zlib=True)
        out27 = n_id.createVariable('Raw_Land_Aerosol_Optical_Depth_StDev','f4',('lon','lat','hour'),zlib=True)

# Define the output variables (2-D Mapping)
  
        out28 = n_id.createVariable('Lon2D','f4',('lon','lat'),zlib=True)
        out29 = n_id.createVariable('Lat2D','f4',('lon','lat'),zlib=True)
    
# Put the data into the output (dimensions)

        out01[:] = glon
        out02[:] = glat
        out03[:] = gtime

# Calculate the mean USER values - ALL

        denom = np.copy(aod_best_tot_cnt)*1.0
        denom[aod_best_tot_cnt == 0] = 1.0
        good = aod_best_tot_cnt > 0
        mean_best_aod = aod_best_tot/denom
        gmean_best_aod = log_aod_best_tot/denom
    
        mean_aod_nan[good] = mean_best_aod[good]
        gmean_aod_nan[good] = 10**gmean_best_aod[good]
        count_aod_nan[good] = aod_best_tot_cnt[good]
        
        keep = aod_best_tot_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_nan[keep] = np.sqrt(aod_best_sk_tot[keep]/(aod_best_tot_cnt[keep]-1.0))

# Calculate the mean USER values - WATER

        denom = np.copy(aod_best_water_cnt)*1.0
        denom[aod_best_water_cnt == 0] = 1.0
        good = aod_best_water_cnt > 0
        mean_aod_water = aod_best_water/denom
        gmean_aod_water = log_aod_best_water/denom
    
        mean_aod_water_nan[good] = mean_aod_water[good]
        gmean_aod_water_nan[good] = 10**gmean_aod_water[good]
        count_aod_water_nan[good] = aod_best_water_cnt[good]
        
        keep = aod_best_water_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_water_nan[keep] = np.sqrt(aod_best_sk_water[keep]/(aod_best_water_cnt[keep]-1.0))

# Calculate the mean USER values - LAND

        denom = np.copy(aod_best_land_cnt)*1.0
        denom[aod_best_land_cnt == 0] = 1.0
        good = aod_best_land_cnt > 0
        mean_aod_land = aod_best_land/denom
        gmean_aod_land = log_aod_best_land/denom
    
        mean_aod_land_nan[good] = mean_aod_land[good]
        gmean_aod_land_nan[good] = 10**gmean_aod_land[good]
        count_aod_land_nan[good] = aod_best_land_cnt[good]
        
        keep = aod_best_land_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_land_nan[keep] = np.sqrt(aod_best_sk_land[keep]/(aod_best_land_cnt[keep]-1.0))

# Calculate the mean RAW values - ALL

        denom = np.copy(aod_all_tot_cnt)*1.0
        denom[aod_all_tot_cnt == 0] = 1.0
        good = aod_all_tot_cnt > 0
        mean_raw_aod = aod_all_tot/denom
        gmean_raw_aod = log_aod_all_tot/denom
    
        mean_aod_raw_nan[good] = mean_raw_aod[good]
        gmean_aod_raw_nan[good] = 10**gmean_raw_aod[good]
        count_aod_raw_nan[good] = aod_all_tot_cnt[good]
        
        keep = aod_all_tot_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_raw_nan[keep] = np.sqrt(aod_all_sk_tot[keep]/(aod_all_tot_cnt[keep]-1.0))

# Calculate the mean RAW values - WATER

        denom = np.copy(aod_all_water_cnt)*1.0
        denom[aod_all_water_cnt == 0] = 1.0
        good = aod_all_water_cnt > 0
        mean_raw_aod_water = aod_all_water/denom
        gmean_raw_aod_water = log_aod_all_water/denom
    
        mean_aod_raw_water_nan[good] = mean_raw_aod_water[good]
        gmean_aod_raw_water_nan[good] = 10**gmean_raw_aod_water[good]
        count_aod_raw_water_nan[good] = aod_all_water_cnt[good]
        
        keep = aod_all_water_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_raw_water_nan[keep] = np.sqrt(aod_all_sk_water[keep]/(aod_all_water_cnt[keep]-1.0))

# Calculate the mean RAW values - LAND

        denom = np.copy(aod_all_land_cnt)*1.0
        denom[aod_all_land_cnt == 0] = 1.0
        good = aod_all_land_cnt > 0
        mean_raw_aod_land = aod_all_land/denom
        gmean_raw_aod_land = log_aod_all_land/denom
    
        mean_aod_raw_land_nan[good] = mean_raw_aod_land[good]
        gmean_aod_raw_land_nan[good] = 10**gmean_raw_aod_land[good]
        count_aod_raw_land_nan[good] = aod_all_land_cnt[good]
        
        keep = aod_all_land_cnt > 1  # Do not calculate variance for a single value
        stdev_aod_raw_land_nan[keep] = np.sqrt(aod_all_sk_land[keep]/(aod_all_land_cnt[keep]-1.0))
        
# Store the data

        out04[:] = mean_aod_nan
        out05[:] = gmean_aod_nan
        out06[:] = count_aod_nan
        out07[:] = stdev_aod_nan
        
        out08[:] = mean_aod_water_nan
        out09[:] = gmean_aod_water_nan
        out10[:] = count_aod_water_nan
        out11[:] = stdev_aod_water_nan
        
        out12[:] = mean_aod_land_nan
        out13[:] = gmean_aod_land_nan
        out14[:] = count_aod_land_nan
        out15[:] = stdev_aod_land_nan
    
        out16[:] = mean_aod_raw_nan
        out17[:] = gmean_aod_raw_nan
        out18[:] = count_aod_raw_nan
        out19[:] = stdev_aod_raw_nan
        
        out20[:] = mean_aod_raw_water_nan
        out21[:] = gmean_aod_raw_water_nan
        out22[:] = count_aod_raw_water_nan
        out23[:] = stdev_aod_raw_water_nan
        
        out24[:] = mean_aod_raw_land_nan
        out25[:] = gmean_aod_raw_land_nan
        out26[:] = count_aod_raw_land_nan
        out27[:] = stdev_aod_raw_land_nan
        
# Put the 2-D Mapping arrays into the file

        out28[:] = lon2D
        out29[:] = lat2D  
    
# Close the NetCDF file
        
        n_id.close()
    
# Print the time

    all_end_time = time.time()
    print("Total elapsed time was %g seconds" % (all_end_time - all_start_time))

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION


def daterange(start_date, end_date):
### Iterate over a range of dates
# Taken from 
# https://stackoverflow.com/questions/1060279/iterating-through-a-range-of-dates-in-python
#
# Input: start_date = date to start, format date(YYYY, MM, DD)
#        end_date = date to end, format date(YYYY, MM, DD)
# Output: List of dates in string YYYY-MM-DD

    for n in range(int((end_date-start_date).days)):
        yield start_date + timedelta(n)


### END daterange

def misr_datetime64_to_orbit(this_datetime64):
### Convert a Numpy datetime64 object to a MISR orbit
# This is the inverse of the misr_orbit_to_datetime function.
# Creation Date: 2021-05-12
# Last Modified: 2021-05-12
# By Michael J. Garay
# Michael.J.Garay@jpl.nasa.gov
# 
# Input: this_datetime64 = datetime64 object
# Output: Orbit (int) = MISR orbit number

# First MISR orbit: Dec. 18.4, 1999
    
    first_datetime_str = "1999-12-18T09:36:00"
    first_datetime64 = np.datetime64(first_datetime_str)

# Compute the elapsed time in seconds

    delta_sec = (this_datetime64 - first_datetime64)/np.timedelta64(1,'s')
#    delta_sec = (this_datetime64 - first_datetime64)/np.timedelta64(1,'ms')

# Calculate the number of orbits this is using the number of seconds per orbit
# NOTE: There was a small change to the number of seconds per orbit to get consistency
#       with the results from the MISR orbit/date tool.
    
    num_orbits = delta_sec/5933.0
#    num_orbits = delta_sec/5933140.  # Original value in ms
#    num_orbits = delta_sec/5933000.  # For calculation in ms
    now = int(1+num_orbits)
    
# Return result

    return now

### END misr_datetime_to_orbit

if __name__ == '__main__':
    main()