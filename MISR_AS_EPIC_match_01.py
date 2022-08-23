# MISR_AS_EPIC_match_01.py
#
# This is a Python 3.9.13 code to read a MISR V23 aerosol file and an associated EPIC
# MAIAC aerosol file, extract the data, and find matches. 
#
# Creation Date: 2022-08-23
# Last Modified: 2022-08-23
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

    misrpath = '/Volumes/JPL_2022/DATA/MISR_FIRES/'
    epicpath = '/Volumes/JPL_2022/DATA/EPIC_FIRES/'
    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY39/AUG22/MISR_EPIC/DATA/'

# Set the MISR orbit
  
    misr_orbit = '110564'  # Mullen Fire
    
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

# Print the time

    all_end_time = time.time()
    print("Total elapsed time was %g seconds" % (all_end_time - all_start_time))

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION


if __name__ == '__main__':
    main()    