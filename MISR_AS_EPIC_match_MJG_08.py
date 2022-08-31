# MISR_AS_EPIC_match_MJG_08.py
#
# This is a Python 3.10.6 code to read a MISR V23 aerosol file and an associated EPIC
# MAIAC aerosol file, extract the data, and find matches. 
#
# Creation Date: 2022-08-31
# Last Modified: 2022-08-31
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

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

# Set the paths
# NOTE: misrpath = Location of the MISR aerosol data
#       epicpath = Location of the EPIC MAIAC data
#       datapath = Location of the output data file

    misrpath = '/Users/mgaray/Desktop/DATA/MISR_FIRES/'
    epicpath = '/Users/mgaray/Desktop/DATA/EPIC_FIRES/'
    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY310/AUG22/MISR_EPIC/DATA/'

# Set the MISR orbit
  
    misr_orbit = '110564'  # Mullen Fire
    start_block = 58 # 1-based
    end_block = 62 # 1-based

# Set the EPIC data
    
    epic_spec = '20200930181817'
    
### FIND THE MISR AEROSOL DATA

# Change to the correct directory

    os.chdir(misrpath)

# Search for the individual MISR files

    search_str = 'MISR_AM1_AS_AEROSOL_*'+misr_orbit+'*.nc'
    dum_list = glob.glob(search_str)
        
    if(len(dum_list) < 1):
        print("***ERROR***")
        print("MISSING ORBIT")
        print(misr_orbit)
        print("***ERROR***")
        print()
        print(error)
        
# Select the file

    misr_file = dum_list[0]
            
# Tell user location in process

    print("Reading: "+misr_file)
        
# Open the NetCDF file

    rootgrp = Dataset(misr_file, 'r', format='NETCDF4')
    
# Choose the appropriate subgroup

    subgrp = rootgrp.groups['4.4_KM_PRODUCTS']

# Assign the variables to (masked) arrays

    block_num = subgrp.variables['Block_Number'][:]

    lat_masked = subgrp.variables['Latitude'][:]
    lon_masked = subgrp.variables['Longitude'][:]
        
    year_masked = subgrp.variables['Year'][:]
    month_masked = subgrp.variables['Month'][:]
    day_masked = subgrp.variables['Day'][:]
            
    hour_masked = subgrp.variables['Hour'][:]
    minute_masked = subgrp.variables['Minute'][:]
               
    aod_masked = subgrp.variables['Aerosol_Optical_Depth'][:]
    ssc_masked = subgrp.variables['Spectral_AOD_Scaling_Coeff'][:]
            
# Choose the AUXILIARY subsubgroup

    subsubgrp = subgrp.groups['AUXILIARY']
        
    aod_raw_masked = subsubgrp.variables['Aerosol_Optical_Depth_Raw'][:]
    ssa_raw_masked = subsubgrp.variables['Single_Scattering_Albedo_446nm_Raw'][:]
    ssc_raw_masked = subsubgrp.variables['Spectral_AOD_Scaling_Coeff_Raw'][:]
    
    arsf_masked = subsubgrp.variables['Aerosol_Retrieval_Screening_Flags'][:]
    
    lrm_masked = subsubgrp.variables['Lowest_Residual_Mixture'][:]
                             
# Close the NetCDF file

    rootgrp.close()
    
# Convert the masked arrays to numpy arrays that can be properly analyzed

    misr_lat = ma.filled(lat_masked,fill_value=-9999.0)
    misr_lon = ma.filled(lon_masked,fill_value=-9999.0)
        
    year = ma.filled(year_masked,fill_value=-9999.0)
    month = ma.filled(month_masked,fill_value=-9999.0)
    day = ma.filled(day_masked,fill_value=-9999.0)
            
    hour = ma.filled(hour_masked,fill_value=-9999.0)
    minute = ma.filled(minute_masked,fill_value=-9999.0)
    
    misr_aod_550 = ma.filled(aod_masked,fill_value=-9999.0)
    misr_aod_scale = ma.filled(ssc_masked,fill_value=-9999.0)
    
    misr_aod_raw_550 = ma.filled(aod_masked,fill_value=-9999.0)
    misr_aod_raw_scale = ma.filled(ssc_masked,fill_value=-9999.0)
    misr_ssa_raw_446 = ma.filled(aod_masked,fill_value=-9999.0)
    
    misr_flags = ma.filled(arsf_masked,fill_value=253)
    misr_mix = ma.filled(lrm_masked,fill_value=-9999.0)

### FIND THE EPIC MAIAC AEROSOL DATA

# Change to the correct directory

    os.chdir(epicpath)

# Search for the individual MISR files

    search_str = 'EPIC_MAIAC_speciation_*'+epic_spec+'.nc'
    dum_list = glob.glob(search_str)
        
    if(len(dum_list) < 1):
        print("***ERROR***")
        print("MISSING ORBIT")
        print(epic_spec)
        print("***ERROR***")
        print()
        print(error)
        
# Select the file

    epic_file = dum_list[0]
            
# Tell user location in process

    print("Reading: "+epic_file)
    
# Open the NetCDF file

    rootgrp = Dataset(epic_file, 'r', format='NETCDF4')
    
# Choose the appropriate subgroup

    subgrp = rootgrp.groups['Geolocation_fields']

# Assign the variables to (masked) arrays

    lat_masked = subgrp.variables['latitude'][:]
    lon_masked = subgrp.variables['longitude'][:]
    
# Choose the appropriate subgroup

    subgrp = rootgrp.groups['Height_1km']
    
    aod_1km_masked = subgrp.variables['AOD_443_1km'][:]
    aod_1km_sf = getattr(subgrp.variables['AOD_443_1km'],'scale_factor')
    ssa_1km_masked = subgrp.variables['SSA_443_1km'][:]
    ssa_1km_sf = getattr(subgrp.variables['SSA_443_1km'],'scale_factor')
    
# Choose the appropriate subgroup

    subgrp = rootgrp.groups['Height_4km']
    
    aod_4km_masked = subgrp.variables['AOD_443_4km'][:]
    aod_4km_sf = getattr(subgrp.variables['AOD_443_4km'],'scale_factor')
    ssa_4km_masked = subgrp.variables['SSA_443_4km'][:]
    ssa_4km_sf = getattr(subgrp.variables['SSA_443_4km'],'scale_factor')

# Close the NetCDF file

    rootgrp.close()
    
# Convert the masked arrays to numpy arrays that can be properly analyzed
# NOTE: Do NOT apply the scale factors

    epic_lat = ma.filled(lat_masked,fill_value=-9999.0)
    epic_lon = ma.filled(lon_masked,fill_value=-9999.0)
    
    epic_aod_1km = ma.filled(aod_1km_masked,fill_value=-9999.0)
    epic_ssa_1km = ma.filled(ssa_1km_masked,fill_value=-9999.0)
    
    epic_aod_4km = ma.filled(aod_4km_masked,fill_value=-9999.0)
    epic_ssa_4km = ma.filled(ssa_4km_masked,fill_value=-9999.0)
    
# Get the sizes of the MISR and EPIC arrays

    print("MISR Sizes")
    print(np.shape(misr_aod_550))
    
    print("EPIC Sizes")
    print(np.shape(epic_aod_1km))

### SUBSET THE MISR DATA ON THE CORRECT BLOCK RANGE
# NOTE: We do this so that we can restrict the MISR/EPIC comparison to plumes
#       of interest

    block_range = (block_num >= start_block) & (block_num <= end_block)
    block_extract = block_num[block_range]
    block_min_ind = np.argwhere(block_num == np.amin(block_extract)) # Extract index
    block_max_ind = np.argwhere(block_num == np.amax(block_extract)) # Extract index
    
# Convert indices to array range for the aerosol datasets by multiplying by 32

    full_min_ind = int(block_min_ind*32)
    full_max_ind = int((block_max_ind+1)*32-1)  # Get the next block and subtract one
    
# Extract the subset of the MISR data

    misr_lat_sub = misr_lat[full_min_ind:full_max_ind,:]
    misr_lon_sub = misr_lon[full_min_ind:full_max_ind,:]
    
    misr_aod_raw_550_sub = misr_aod_raw_550[full_min_ind:full_max_ind,:]
    misr_aod_raw_scale_sub = misr_aod_raw_scale[full_min_ind:full_max_ind,:]
    misr_ssa_raw_446_sub = misr_ssa_raw_446[full_min_ind:full_max_ind,:]
    
### GET THE LAT/LON INTERSECTION OF THE MISR AND EPIC DATA

    good = ((misr_lat_sub > -999.) & (misr_lon_sub > -999.))
    
    lat_min = np.amin(misr_lat_sub[good])
    lat_max = np.amax(misr_lat_sub[good])
    
    lon_min = np.amin(misr_lon_sub[good])
    lon_max = np.amax(misr_lon_sub[good])

    inside = ((epic_lat >= lat_min) & (epic_lat <= lat_max) & 
        (epic_lon >= lon_min) & (epic_lon <= lon_max))
        
    epic_lat_sub = epic_lat[inside]
    epic_lon_sub = epic_lon[inside]
    
    epic_aod_1km_sub = epic_aod_1km[inside]
    epic_ssa_1km_sub = epic_ssa_1km[inside]
    
    epic_aod_4km_sub = epic_aod_4km[inside]
    epic_ssa_4km_sub = epic_ssa_4km[inside]
    
### CHOOSE ONLY THOSE MISR PIXELS WITH VALID AOD DATA TO MATCH

    good = (misr_aod_raw_550_sub > 0.0)
    misr_aod_raw_500_good = misr_aod_raw_550_sub[good]
    misr_aod_raw_scale_good = misr_aod_raw_scale_sub[good]
    misr_ssa_raw_446_good = misr_ssa_raw_446_sub[good]
    
    misr_lat_good = misr_lat_sub[good]
    misr_lon_good = misr_lon_sub[good]
    
### CHOOSE ONLY THOSE EPIC PIXELS WITH VALID AOD DATA TO MATCH

    good = ((epic_aod_1km_sub > 0.0) & (epic_aod_4km_sub > 0.0))  # This may be redundant
    
    epic_aod_1km_good = epic_aod_1km_sub[good]
    epic_ssa_1km_good = epic_ssa_1km_sub[good]
    
    epic_aod_4km_good = epic_aod_4km_sub[good]
    epic_ssa_4km_good = epic_ssa_4km_sub[good]
    
    epic_lat_good = epic_lat_sub[good]
    epic_lon_good = epic_lon_sub[good]

    print("MISR")
    print(np.shape(misr_lat))
    print(np.shape(misr_lat_sub))
    print(np.shape(misr_lat_good))
    
    print("EPIC")
    print(np.shape(epic_lat))
    print(np.shape(epic_lat_sub))
    print(np.shape(epic_lat_good))
    print()
    
### MATCH EPIC DATA TO MISR
# NOTE: This is just the nearest neighbor match to the center

    num_epic = len(epic_lat_good)
    print("Matching {:d} MAIAC/EPIC pixels".format(num_epic))
    
# Set arrays to store the data

    epic_aod_1km_match = np.zeros(num_epic)
    epic_ssa_1km_match = np.zeros(num_epic)
    epic_aod_4km_match = np.zeros(num_epic)
    epic_ssa_4km_match = np.zeros(num_epic)
    
    misr_aod_match = np.zeros(num_epic)
    misr_ssa_match = np.zeros(num_epic)
    
    distance_match = np.zeros(num_epic)
    
# Loop over the EPIC data    
    
    for loop in range(num_epic):

       epic_lat = epic_lat_good[loop]  # Note: I'm overwriting the original epic_lat array
       epic_lon = epic_lon_good[loop]  # Note: I'm overwriting the original epic_lon array
       
       misr_dist = Haversine_Distance(epic_lat,epic_lon,misr_lat_good,misr_lon_good)
       
       min_dist = np.amin(misr_dist)  # This is the distance in km
       
       min_ind = np.argwhere(misr_dist == min_dist)  # Get the index
       
       misr_aod_550 = misr_aod_raw_500_good[min_ind[0]]
       misr_aod_scale = misr_aod_raw_scale_good[min_ind[0]]
       misr_ssa = misr_ssa_raw_446_good[min_ind[0]]

# Use the scale factors to calculate the MISR AOD at 443 nm to compare to MAIAC retrievals
       
       misr_sf1 = misr_aod_scale[0,0]
       misr_sf2 = misr_aod_scale[0,1]
       misr_sf3 = misr_aod_scale[0,2]
       
#       misr_550 = misr_sf1*(0.550**2)+misr_sf2*(0.550)+misr_sf3  # Check
       misr_443 = misr_sf1*(0.443**2)+misr_sf2*(0.443)+misr_sf3

# Store the matched data

       epic_aod_1km_match[loop] = epic_aod_1km_good[loop]
       epic_ssa_1km_match[loop] = epic_ssa_1km_good[loop]
       
       epic_aod_4km_match[loop] = epic_aod_4km_good[loop]
       epic_ssa_4km_match[loop] = epic_ssa_4km_good[loop]
       
       misr_aod_match[loop] = misr_443
       misr_ssa_match[loop] = misr_ssa[0]
       
       distance_match[loop] = min_dist
       
    print(np.amin(distance_match))
    print(np.amax(distance_match))
  
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