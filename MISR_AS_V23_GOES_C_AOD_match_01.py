# MISR_AS_V23_GOES_C_AOD_match_01.py
#
# This is a Python 3.6.10 code to read MISR aerosol data
# then generate AOD regression information for GOES-CONUS sites in the Western US
# during FIREX-AQ.
#
# Creation Date: 2020-02-21
# Last Modified: 2020-02-21
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

from datetime import datetime
import fnmatch
import glob
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os
import time

def main():  # Main code

# Set the overall timer

    all_start_time = time.time()

# Set the file exist (0 to start, 1 after initial run)
# NOTE: The NetCDF filename is hardcoded

    exist = 0

# Set some data limits

    max_lat = 49.0
    min_lat = 31.0

# Set the time window

    time_window = 2.5  #  Time in minutes
    
# Set the folder

    this_folder = '2019_08_15'

# Set the paths

    outpath = '/Users/mgaray/Desktop/CODING/PYTHON/PY36/FEB20/FIREX/DATA/'
    misrpath = '/Users/mgaray/Desktop/DATA/FIREX/MISR/'+this_folder+'/'
    goespath = '/Users/mgaray/Desktop/DATA/FIREX/GOES/'+this_folder+'/'

# Set the output file

    outfile = 'MISR_GOES_AOD_'+this_folder+'_v01.nc'

# Set the product names

    misr_name = 'F13_0023' # New Product
 
# Set the GOES filename

    goes_name = 'OR_ABI-L2-AODC' # CONUS Product

# Set up dummy arrays to store the matched data

    misr_lat_match_raw = []
    misr_lon_match_raw = []
    
    goes_lat_match_raw = []
    goes_lon_match_raw = []
    
    dist_match_raw = []
      
    misr_aod_match_raw = []
    misr_raw_match_raw = []
    misr_lwr_match_raw = []
      
    goes_aod_match_raw = []
    goes_dqf_match_raw = []

### GET THE MISR DATA

    os.chdir(misrpath)
        
# Parse the folder to get some information

    hold = this_folder.split('_')
    this_year = int(hold[0])
    this_month = int(hold[1])
    this_day = int(hold[2])
                        
# Choose the correct file

    file_list = glob.glob('*'+misr_name+'.nc')

    num_files = len(file_list)
        
    if(num_files != 1):
        print()
        print("***ERROR***")
        print("MISSING MISR Files...")
        print(misr_name_1)
        print("***ERROR***")
        print()
        print(error)

    misr_file = file_list[0]
    
# Tell user location in process

    print("Reading: "+misr_file)
        
# Parse the filename to get information
    
    words = misr_file.split('.')
    fn = words[0]
    temp = fn.split('_')
    path_str = temp[4]
    misr_path_number = int(path_str[1:4])
    orbit_str = temp[5]
    misr_orbit_number = float(orbit_str[1:7])
                     
# Open the NetCDF file

    rootgrp = Dataset(misr_file, 'r', format='NETCDF4')

# Choose the appropriate subgroup

    subgrp = rootgrp.groups['4.4_KM_PRODUCTS']

# Assign the variables to an array

    lat_raw = subgrp.variables['Latitude'][:]
    lon_raw = subgrp.variables['Longitude'][:]
        
    year_raw = subgrp.variables['Year'][:]
    month_raw = subgrp.variables['Month'][:]
    day_raw = subgrp.variables['Day'][:]
        
    hour_raw = subgrp.variables['Hour'][:]
    minute_raw = subgrp.variables['Minute'][:]
        
    aod_raw = subgrp.variables['Aerosol_Optical_Depth'][:]
    lwr_raw = subgrp.variables['Land_Water_Retrieval_Type'][:]

# Choose the AUXILIARY subsubgroup
# NOTE: This is a typo in the subgroup name in this version of the software

    subsubgrp = subgrp.groups['AUXILIARY']
        
    aod_raw_raw = subsubgrp.variables['Aerosol_Optical_Depth_Raw'][:]

# Close the NetCDF file

    rootgrp.close()

# Convert arrays to numpy arrays for easier handling

    misr_lat = np.copy(lat_raw)
    misr_lon = np.copy(lon_raw)
        
    year = np.copy(year_raw)
    month = np.copy(month_raw)
    day = np.copy(day_raw)
        
    hour = np.copy(hour_raw)
    minute = np.copy(minute_raw)
        
    misr_aod = np.copy(aod_raw)
    misr_aod_raw = np.copy(aod_raw_raw)
    misr_lwr = np.copy(lwr_raw)
    
### CHECK FOR GOES DATA WITHIN THE TIME/SPACE WINDOW

# Find valid MISR lat/lon values in the subset

    good = ((misr_lat >= min_lat) & (misr_lat <= max_lat) & (misr_aod_raw > 0.))

# Extract the values into a subset

    misr_lat_good = misr_lat[good]
    misr_lon_good = misr_lon[good]
        
    year_test = year[good]
    month_test = month[good]
    day_test = day[good]
        
    hour_test = hour[good]
    minute_test = minute[good]
        
    misr_aod_good = misr_aod[good]
    misr_aod_raw_good = misr_aod_raw[good]
    misr_lwr_good = misr_lwr[good]

# Convert the time information to numpy datetime objects

    num_test = len(misr_lon_good)
        
    time_raw = []
        
    for i in range(num_test):
        
        time_str = "{:}-".format(int(year_test[i]))
        time_str = time_str+"{:02}-".format(int(month_test[i]))
        time_str = time_str+"{:02}".format(int(day_test[i]))
            
        time_str = time_str+"T"
        time_str = time_str+"{:02}:".format(int(hour_test[i]))
        time_str = time_str+"{:02}".format(int(minute_test[i]))
        
# Store the time string

        time_raw.append(time_str)    

# Convert the timestring to a numpy datetime object

    misr_time = np.array(time_raw,dtype='datetime64[m]')

# Choose the middle time from the array (more or less)
        
    mean_test = int(num_test/2.)
    misr_mean_time = misr_time[mean_test]
    misr_mean_time_str = time_raw[mean_test]
    misr_mean_year = int(misr_mean_time_str[0:4])
    misr_mean_month = int(misr_mean_time_str[5:7])
    misr_mean_day = int(misr_mean_time_str[8:10])
    misr_mean_hour = int(misr_mean_time_str[11:13])
    misr_mean_minute = int(misr_mean_time_str[14:16])
    
    print(misr_mean_time)
    
# Change directory to the GOES path

    os.chdir(goespath)
        
# Choose the correct files

    file_list = glob.glob(goes_name+'*.nc')

    num_files = len(file_list)
        
    if(num_files < 1):
        print()
        print("***ERROR***")
        print("MISSING GOES Files...")
        print(goes_path)
        print(goes_name)
        print("***ERROR***")
        print()

# Loop through files

    misr_count = 0
    
    for file_loop in range(num_files):

        goes_file = file_list[file_loop]
    
# Parse the filename to get information
    
        words = goes_file.split('_')
        temp = words[3]
        goes_start_hr = int(temp[8:10])
        goes_start_min = int(temp[10:12])
        temp = words[4]
        goes_end_hr = int(temp[8:10])
        goes_end_min = int(temp[10:12])  
        
# Convert the time information into a numpy datetime object

        time_str = "{:}-".format(this_year)
        time_str = time_str+"{:02}-".format(this_month)
        time_str = time_str+"{:02}".format(this_day)
            
        time_str = time_str+"T"
        time_str = time_str+"{:02}:".format(goes_start_hr)
        time_str = time_str+"{:02}".format(goes_start_min)

# Convert the timestring to a numpy datetime object

        goes_time = np.array(time_str,dtype='datetime64[m]')

## Find the GOES/MISR matches within the time window (in minutes)

        tdiff = np.abs(np.divide((misr_mean_time-goes_time),np.timedelta64(1,'m')))        

# Test if any GOES/MISR matches fall within the time window
    
        found = (tdiff <= time_window)
        
        if(np.amax(found) == True):

# If time matches are found, extract the GOES data

            print("Reading: "+goes_file)

# Open the NetCDF file

            rootgrp = Dataset(goes_file, 'r', format='NETCDF4')

# Assign the variables to an array

            aod_raw = rootgrp.variables['AOD'][:]
            dqf_raw = rootgrp.variables['DQF'][:]
            x_raw = rootgrp.variables['x'][:]
            y_raw = rootgrp.variables['y'][:]

# Get the scale factors and add offset for AOD
    
            aod_sf = getattr(rootgrp.variables['AOD'],'scale_factor')
            aod_ao = getattr(rootgrp.variables['AOD'],'add_offset')
    
# Get the projection information

            r_eq = getattr(rootgrp.variables['goes_imager_projection'],'semi_major_axis')
            f_inv = getattr(rootgrp.variables['goes_imager_projection'],'inverse_flattening')
            r_pol = getattr(rootgrp.variables['goes_imager_projection'],'semi_minor_axis')
            pph = getattr(rootgrp.variables['goes_imager_projection'],'perspective_point_height')
            l_0 = getattr(rootgrp.variables['goes_imager_projection'],'longitude_of_projection_origin')
    
# Close the NetCDF file

            rootgrp.close()

# Convert the read arrays to numpy arrays

            aod_bare = np.array(aod_raw)
            goes_dqf = np.array(dqf_raw)
            x_dim = np.array(x_raw)
            y_dim = np.array(y_raw)

# AOD: Scale and add the offset

            goes_aod = aod_bare*aod_sf+aod_ao

# Make x and y into 2-D arrays with the same shape as AOD

            x, y = np.meshgrid(x_dim,y_dim)

### NAVIGATE THE DATA USING THE INFORMATION IN THE PUG-L2+-vol5.pdf
# Section 4.2.8.1

# Calculate the missing navigation parameters

            e = np.sqrt((r_eq**2-r_pol**2)/r_eq**2)
            H = pph+r_eq
            ratio_r2 =  (r_eq/r_pol)**2  # This term is always used in the squared form 
            c = H**2-r_eq**2

# Calculate the trig functions of x and y (these are in radians)

            sin_x = np.sin(x)
            cos_x = np.cos(x)
            sin_y = np.sin(y)
            cos_y = np.cos(y)
    
# Calculate the navigation equations

            a = sin_x**2+cos_x**2*(cos_y**2+ratio_r2*sin_y**2)
            b = -2.0*H*cos_x*cos_y
            r_s = (-1.0*b-np.sqrt(b**2-4.0*a*c))/(2.0*a)
            s_x = r_s*cos_x*cos_y
            s_y = -1.0*r_s*sin_x
            s_z = r_s*cos_x*sin_y
    
            tan_phi = ratio_r2*(s_z/np.sqrt((H-s_x)**2+s_y**2))
            phi = np.arctan(tan_phi)
    
            tan_paren = s_y/(H-s_x)
            paren = np.arctan(tan_paren)
            lamb = np.deg2rad(l_0)-paren  # Note the need to convert l_0 to radians
            
            goes_lat = np.rad2deg(phi)
            goes_lon = np.rad2deg(lamb)
            
# Loop through the valid MISR data and match to GOES

            num_good = len(misr_aod_good)
            
            for inner in np.arange(num_good):
            
                misr_lat_test = misr_lat_good[inner]
                misr_lon_test = misr_lon_good[inner]
                misr_aod_test = misr_aod_good[inner]
                misr_aod_raw_test = misr_aod_raw_good[inner]
                misr_lwr_test = misr_lwr_good[inner]
                
# Use the GOES INVERSE OPERATION (PUG-L2+-vol5.pdf Section 4.2.8.2) to find the match for
# this location
                    
                phi = np.deg2rad(misr_lat_test)
                lam = np.deg2rad(misr_lon_test)
                tan_phi_c = np.tan(phi)/ratio_r2
                phi_c = np.arctan(tan_phi_c)
                lam_c = lam-np.deg2rad(l_0)
    
                cos_phi = np.cos(phi_c)
                sin_phi = np.sin(phi_c)
                cos_lam = np.cos(lam_c)
                sin_lam = np.sin(lam_c)
    
                r_c = r_pol/np.sqrt(1.0-e**2*cos_phi**2)
                s_x = H - r_c*cos_phi*cos_lam
                s_y = -1.0*r_c*cos_phi*sin_lam
                s_z = r_c*sin_phi
    
# Test for visibility

                left_test = H*(H-s_x)
                right_test = s_y**2+ratio_r2*s_z**2
    
                if(right_test > left_test):
                    print("Point not visible")
                    continue
        
# Calculate y and x (in radians)

                tan_yr = s_z/s_x
                yr = np.arctan(tan_yr)
                sin_xr = -1.0*s_y/(np.sqrt(s_x**2+s_y**2+s_z**2))
                xr = np.arcsin(sin_xr)
    
# Find the index on the grid

                del_x = np.abs(xr-x_dim)
                x_val = np.argmin(del_x)
                del_y = np.abs(yr-y_dim)
                y_val = np.argmin(del_y)
                    
# Extract the latitude and longitude

                goes_lat_found = goes_lat[y_val,x_val]
                goes_lon_found = goes_lon[y_val,x_val]
                goes_aod_found = goes_aod[y_val,x_val]
                goes_dqf_found = goes_dqf[y_val,x_val]
                goes_dist_found = Haversine_Distance(misr_lat_test,
                        misr_lon_test,goes_lat[y_val,x_val],goes_lon[y_val,x_val]) 

# Store the data
                
                misr_lat_match_raw.append(misr_lat_test)
                misr_lon_match_raw.append(misr_lon_test)
    
                goes_lat_match_raw.append(goes_lat_found)
                goes_lon_match_raw.append(goes_lon_found)
                dist_match_raw.append(goes_dist_found)
                    
                goes_aod_match_raw.append(goes_aod_found)
                goes_dqf_match_raw.append(goes_dqf_found)
      
                misr_aod_match_raw.append(misr_aod_test)
                misr_raw_match_raw.append(misr_aod_raw_test)
                misr_lwr_match_raw.append(misr_lwr_test)

                misr_count = misr_count+1
            
        else:
            print("No Matching Data Found in Time Window")
            continue
                
### OUTPUT THE DATA

# Change to the correct directory

    os.chdir(outpath)
            
# Extract the data

    if(misr_count > 0):
        
# Convert the raw data to numpy arrays
            
        goes_lat_match = np.array(goes_lat_match_raw)
        goes_lon_match = np.array(goes_lon_match_raw)
        
        misr_lat_match = np.array(misr_lat_match_raw)
        misr_lon_match = np.array(misr_lon_match_raw)
        
        dist_match = np.array(dist_match_raw)
            
        goes_aod_match = np.array(goes_aod_match_raw)
        goes_dqf_match = np.array(goes_dqf_match_raw)
        
        misr_aod_match = np.array(misr_aod_match_raw)
        misr_raw_match = np.array(misr_raw_match_raw)
        misr_lwr_match = np.array(misr_lwr_match_raw)

# Open the file
        
        n_id = Dataset(outfile,'w')
        
# Define the dimensions

        xdim = n_id.createDimension('xdim') # Unlimited  

# Define the output variables

        out01 = n_id.createVariable('GOES_Longitude','f4',('xdim',),zlib=True)
        out02 = n_id.createVariable('GOES_Latitude','f4',('xdim',),zlib=True)
        
        out03 = n_id.createVariable('MISR_Longitude','f4',('xdim',),zlib=True)
        out04 = n_id.createVariable('MISR_Latitude','f4',('xdim',),zlib=True)
        
        out05 = n_id.createVariable('Distance','f4',('xdim',),zlib=True)
                
        out06 = n_id.createVariable('GOES_AOD_550nm','f4',('xdim',),zlib=True)
        out07 = n_id.createVariable('GOES_DQF','f4',('xdim',),zlib=True)
        
        out08 = n_id.createVariable('MISR_AOD_550nm','f4',('xdim',),zlib=True)
        out09 = n_id.createVariable('MISR_AOD_RAW_550nm','f4',('xdim',),zlib=True)
        out10 = n_id.createVariable('Land_Water_Retrieval_Type','f4',('xdim',),zlib=True)

# Put the data into the output
             
        out01[:] = goes_lat_match
        out02[:] = goes_lon_match
        
        out03[:] = misr_lat_match
        out04[:] = misr_lon_match
        
        out05[:] = dist_match
                
        out06[:] = goes_aod_match
        out07[:] = goes_dqf_match
        
        out08[:] = misr_aod_match
        out09[:] = misr_raw_match
        out10[:] = misr_lwr_match

# Close the NetCDF file
        
        n_id.close()

# Print the time

    all_end_time = time.time()
    print("Total elapsed time was %g seconds" % (all_end_time - all_start_time))

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
