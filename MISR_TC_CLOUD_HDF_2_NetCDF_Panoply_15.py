# MISR_TC_CLOUD_HDF_2_NetCDF_Panoply_15.py
#
# This is a Python 3.9.13 code to read a MISR TC_CLOUD file and associated AGP file for
# the same path and generate a navigated NetCDF output file that can be visualized
# using Panoply. 
#
# Creation Date: 2022-08-04
# Last Modified: 2022-08-12
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

import glob
from netCDF4 import Dataset
import numpy as np
import os
from pyhdf.HDF import *
from pyhdf.SD import SD, SDC
import time

def main():  # Main code

# Print a blank line

    print()

# Set the overall timer

    all_start_time = time.time()  

# Set the paths
# NOTE: basepath = Location of the TC_CLOUD data
#       agppath = Location of the AGP data
#       datapath = Location of the output data file(s)

    basepath = '/Volumes/JPL_2022/DATA/MISR_FIRES/'
    agppath = '/Volumes/JPL_2022/DATA/AGP/'
    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY39/AUG22/TC_CLOUD/DATA/'

# Set the MISR path and orbit of interest and the start and end blocks
    
#    misr_path = 'P035'  # Silver Creek Fire
#    misr_orbit = '099686'
#    start_block = 52 # 1-based
#    end_block = 60 # 1-based
    
#    misr_path = 'P044'  # Carr Fire (1)
#    misr_orbit = '098973'
#    start_block = 55 # 1-based
#    end_block = 61 # 1-based
 
#    misr_path = 'P045'  # Carr Fire (2)
#    misr_orbit = '099075'
#    start_block = 55 # 1-based
#    end_block = 62 # 1-based
    
#    misr_path = 'P044'  # Carr Fire (3)
#    misr_orbit = '099206'
#    start_block = 50 # 1-based
#    end_block = 61 # 1-based
    
#    misr_path = 'P043'  # Donnell Fire
#    misr_orbit = '099104'
#    start_block = 55 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P046'  # Taylor Creek Fire
#    misr_orbit = '099177'
#    start_block = 43 # 1-based
#    end_block = 60 # 1-based
    
#    misr_path = 'P042'  # Roosters Comb Fire (1)
#    misr_orbit = '093410'
#    start_block = 56 # 1-based
#    end_block = 58 # 1-based

#    misr_path = 'P040'  # Roosters Comb Fire (2)
#    misr_orbit = '093439'
#    start_block = 57 # 1-based
#    end_block = 58 # 1-based
   
#    misr_path = 'P034'  # Mullen Fire (1)
#    misr_orbit = '110535'
#    start_block = 57 # 1-based
#    end_block = 58 # 1-based
  
#    misr_path = 'P032'  # Mullen Fire (2)
#    misr_orbit = '110564'
#    start_block = 58 # 1-based
#    end_block = 62 # 1-based
    
#    misr_path = 'P035'  # Mullen Fire (3) + Friend
#    misr_orbit = '110637'
#    start_block = 56 # 1-based
#    end_block = 65 # 1-based
    
#    misr_path = 'P033'  # Mullen Fire (4)
#    misr_orbit = '110666'
#    start_block = 56 # 1-based
#    end_block = 61 # 1-based
    
#    misr_path = 'P044'  # Evans Canyon Fire
#    misr_orbit = '110157'
#    start_block = 53 # 1-based
#    end_block = 55 # 1-based
    
#    misr_path = 'P040'  # Pioneer Fire (1)
#    misr_orbit = '088313'
#    start_block = 55 # 1-based
#    end_block = 56 # 1-based
    
#    misr_path = 'P041'  # Pioneer Fire (2)
#    misr_orbit = '088415'
#    start_block = 53 # 1-based
#    end_block = 57 # 1-based
    
#    misr_path = 'P039'  # Pioneer Fire (3)
#    misr_orbit = '088444'
#    start_block = 53 # 1-based
#    end_block = 57 # 1-based
    
#    misr_path = 'P043'  # Windy Fire (1)
#    misr_orbit = '115647'
#    start_block = 56 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P041'  # Windy Fire (2)
#    misr_orbit = '115676'
#    start_block = 54 # 1-based
#    end_block = 62 # 1-based
    
#    misr_path = 'P044'  # Windy Fire (3)
#    misr_orbit = '115749'
#    start_block = 53 # 1-based
#    end_block = 61 # 1-based
    
#    misr_path = 'P042'  # Windy Fire (4)
#    misr_orbit = '115778'
#    start_block = 59 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P043'  # Windy Fire (5)
#    misr_orbit = '115880'
#    start_block = 60 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P041'  # Windy Fire (6)
#    misr_orbit = '115909'
#    start_block = 54 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P034'  # E Troublesome Fire (1)
#    misr_orbit = '110768'
#    start_block = 58 # 1-based
#    end_block = 65 # 1-based
    
#    misr_path = 'P026'  # E Troublesome Fire (2)
#    misr_orbit = '110884'
#    start_block = 51 # 1-based
#    end_block = 56 # 1-based
    
#    misr_path = 'P045'  # Williams Flats Fire (1)
#    misr_orbit = '104434'
#    start_block = 51 # 1-based
#    end_block = 52 # 1-based
    
#    misr_path = 'P043'  # Williams Flats Fire (2)
#    misr_orbit = '104463'
#    start_block = 47 # 1-based
#    end_block = 53 # 1-based
    
#    misr_path = 'P038'  # Little Bear + Sheridan Fires
#    misr_orbit = '104652'
#    start_block = 60 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P046'  # Milepost 97 Fire (1)
#    misr_orbit = '104303'
#    start_block = 53 # 1-based
#    end_block = 58 # 1-based
    
#    misr_path = 'P044'  # Milepost 97 Fire (2)
#    misr_orbit = '104332'
#    start_block = 53 # 1-based
#    end_block = 58 # 1-based
    
#    misr_path = 'P045'  # Kincade Fire
#    misr_orbit = '105599'
#    start_block = 56 # 1-based
#    end_block = 61 # 1-based
    
#    misr_path = 'P043'  # Caples Fire
#    misr_orbit = '105395'
#    start_block = 56 # 1-based
#    end_block = 63 # 1-based
    
#    misr_path = 'P039'  # Horse River Fire
#    misr_orbit = '087279'
#    start_block = 45 # 1-based
#    end_block = 52 # 1-based
    
#    misr_path = 'P043'  # Camp Fire (1)
#    misr_orbit = '100502'
#    start_block = 58 # 1-based
#    end_block = 64 # 1-based
    
#    misr_path = 'P044'  # Camp Fire (2)
#    misr_orbit = '100604'
#    start_block = 58 # 1-based
#    end_block = 62 # 1-based
    
#    misr_path = 'P042'  # Camp Fire (3)
#    misr_orbit = '100633'
#    start_block = 57 # 1-based
#    end_block = 62 # 1-based
    
#    misr_path = 'P045'  # Pearl Hill Fire
#    misr_orbit = '110259'
#    start_block = 51 # 1-based
#    end_block = 61 # 1-based
    
#    misr_path = 'P046'  # Alberta, Canada
#    misr_orbit = '114555'
#    start_block = 35 # 1-based
#    end_block = 58 # 1-based
    
#    misr_path = 'P037'  # Saskatchewan
#    misr_orbit = '114569'
#    start_block = 41 # 1-based
#    end_block = 45 # 1-based
    
#    misr_path = 'P041'  # BC (1)
#    misr_orbit = '093541'
#    start_block = 46 # 1-based
#    end_block = 60 # 1-based
    
#    misr_path = 'P046'  # BC (2)
#    misr_orbit = '093818'
#    start_block = 48 # 1-based
#    end_block = 58 # 1-based
    
#    misr_path = 'P044'  # BC (3)
#    misr_orbit = '093847'
#    start_block = 49 # 1-based
#    end_block = 56 # 1-based
    
#    misr_path = 'P051'  # BC (4)
#    misr_orbit = '093862'
#    start_block = 41 # 1-based
#    end_block = 50 # 1-based
    
#    misr_path = 'P047'  # BC (5)
#    misr_orbit = '099279'
#    start_block = 43 # 1-based
#    end_block = 58 # 1-based
    
    misr_path = 'P042'  # Oak Fire
    misr_orbit = '120205'
    start_block = 57 # 1-based
    end_block = 61 # 1-based


    
# Get the software version number to help track issues

    hold = os.path.basename(__file__)
    words = hold.split('_')
    temp = words[len(words)-1]  # Choose the last element
    hold = temp.split('.')
    vers = hold[0]

# Generate the output file

    outfile = 'MISR_AM1_TC_CLOUD_'+misr_path+'_O'+misr_orbit+'_Panoply_v'+vers+'.nc'

### GET THE MISR AGP DATA (NAVIGATION AND ALTITUDE)

# Change to the correct directory

    os.chdir(agppath)

# Get the correct file

    search_str = 'MISR_AM1_AGP*'+misr_path+'*.hdf'

    file_list = glob.glob(search_str)

# Set the filename

    inputName = file_list[0]

# Tell user location in process

    print("Reading: "+inputName)

# Open the file

    hdf = SD(inputName, SDC.READ)

# Read the data fields

    var01 = hdf.select('GeoLatitude')
    var02 = hdf.select('GeoLongitude')
    var03 = hdf.select('AveSceneElev')

    agp_lat = var01.get()
    agp_lon = var02.get()
    agp_ele = var03.get()

# Close the file

    hdf.end()
    
### GET THE MISR TC_CLOUD DATA

# Change to the correct directory

    os.chdir(basepath)

    search_str = 'MISR_AM1_TC_CLOUD_'+misr_path+'_O'+misr_orbit+'*.hdf'
    dum_list = glob.glob(search_str)
    raw_list = np.array(dum_list)  # Convert to a numpy array

# Get the number of files
    
    num_files = len(raw_list)
    
    if(num_files < 1):
        print()
        print("***ERROR***")
        print("CANNOT FIND MISR TC_CLOUD FILE")
        print(search_str)
        print("***ERROR***")
        print()
        print(error)

# Set the filename

    inputName = raw_list[0]

# Tell user location in process

    print("Reading: "+inputName)

# Open the file

    hdf = SD(inputName, SDC.READ)

# Read the data fields

    var01 = hdf.select('CloudTopHeight_WithoutWindCorrection')
    var02 = hdf.select('StereoDerivedCloudMask_WithoutWindCorrection')
    
    var03 = hdf.select('CloudTopHeight')
    var04 = hdf.select('StereoDerivedCloudMask')
    
    var05 = hdf.select('CloudTopHeightOfMotion')
    var06 = hdf.select('MotionQualityIndicator')
    var07 = hdf.select('CloudMotionEastward')
    var08 = hdf.select('CloudMotionNorthward')
    var09 = hdf.select('MotionDerivedCloudMask')
    
    cth_wwc = var01.get()
    sdcm_wwc = var02.get()
    
    cth = var03.get()
    sdcm = var04.get()
    
    cthom = var05.get()
    motqa = var06.get()
    cme = var07.get()
    cmn = var08.get()
    mdcm = var09.get()

# Close the file

    hdf.end()

### EXTRACT MULTIPLE BLOCKS

    s_blk = start_block - 1 # 0-based
    num_blk = end_block - start_block + 1 # Count number of blocks +1 for loop control

# Extract the navigation data for the first block
    
    this_lat = agp_lat[s_blk,:,:]
    this_lon = agp_lon[s_blk,:,:]
    this_ele = agp_ele[s_blk,:,:]
    
# Extract the stereo data for the first block

    this_cth_wwc = cth_wwc[s_blk,:,:]
    this_sdcm_wwc = sdcm_wwc[s_blk,:,:]
    
    this_cth = cth[s_blk,:,:]
    this_sdcm = sdcm[s_blk,:,:]
    
    small_cthom = cthom[s_blk,:,:]
    small_motqa = motqa[s_blk,:,:]
    small_cme = cme[s_blk,:,:]
    small_cmn = cmn[s_blk,:,:]
    small_mdcm = mdcm[s_blk,:,:]
    
# Calculate the wind speed

    small_wspd = np.sqrt(small_cme*small_cme+small_cmn*small_cmn)
    
# Reduce the resolution of the surface elevation dataset to match the winds
# NOTE: The factor is 16 in both dimensions (1.1 km -> 17.6 km)

    small_y = small_wspd.shape[0]
    small_x = small_wspd.shape[1]
    
    full_y = this_ele.shape[0]
    full_x = this_ele.shape[1]
    
    factor_y = int(full_y/small_y)  # Conversion factor
    factor_x = int(full_x/small_x)  # Conversion factor

    small_ele = np.squeeze(this_ele.reshape([small_y,factor_y,
        small_x,factor_x]).mean(3).mean(1))
    
# Resize the wind data to match the navigation
# Note: The factor is 16 in both dimensions (17.6 -> 1.1 km)
    
    this_cthom = np.repeat(np.repeat(small_cthom,16,axis=0),16,axis=1)
    this_motqa = np.repeat(np.repeat(small_motqa,16,axis=0),16,axis=1)
    this_wspd = np.repeat(np.repeat(small_wspd,16,axis=0),16,axis=1)
    this_mdcm = np.repeat(np.repeat(small_mdcm,16,axis=0),16,axis=1)
    this_wele = np.repeat(np.repeat(small_ele,16,axis=0),16,axis=1)

# Loop through the blocks
    
    for loop in range(1,num_blk):
    
        this_blk = s_blk + loop
    
# Extract the navigation data

        temp_lat = agp_lat[this_blk,:,:]
        temp_lon = agp_lon[this_blk,:,:]
        temp_ele = agp_ele[this_blk,:,:]
        
# Extract the stereo data

        temp_cth_wwc = cth_wwc[this_blk,:,:]
        temp_sdcm_wwc = sdcm_wwc[this_blk,:,:]
    
        temp_cth = cth[this_blk,:,:]
        temp_sdcm = sdcm[this_blk,:,:]
    
        small_cthom = cthom[this_blk,:,:]
        small_motqa = motqa[this_blk,:,:]
        small_cme = cme[this_blk,:,:]
        small_cmn = cmn[this_blk,:,:]
        small_mdcm = mdcm[this_blk,:,:]
        
# Calculate the wind speed

        small_wspd = np.sqrt(small_cme*small_cme+small_cmn*small_cmn)

# Reduce the resolution of the surface elevation dataset to match the winds
# NOTE: The factor is 16 in both dimensions (1.1 km -> 17.6 km)

        small_y = small_wspd.shape[0]
        small_x = small_wspd.shape[1]
    
        full_y = temp_ele.shape[0]
        full_x = temp_ele.shape[1]
    
        factor_y = int(full_y/small_y)  # Conversion factor
        factor_x = int(full_x/small_x)  # Conversion factor

        small_ele = np.squeeze(temp_ele.reshape([small_y,factor_y,
            small_x,factor_x]).mean(3).mean(1))
    
# Resize the wind data to match the navigation
# Note: The factor is 16 in both dimensions (17.6 -> 1.1 km)
    
        temp_cthom = np.repeat(np.repeat(small_cthom,16,axis=0),16,axis=1)
        temp_motqa = np.repeat(np.repeat(small_motqa,16,axis=0),16,axis=1)
        temp_wspd = np.repeat(np.repeat(small_wspd,16,axis=0),16,axis=1)
        temp_mdcm = np.repeat(np.repeat(small_mdcm,16,axis=0),16,axis=1)
        temp_wele = np.repeat(np.repeat(small_ele,16,axis=0),16,axis=1)
        
# Stack the data into a 2-D arrays

        this_lat = np.vstack((this_lat,temp_lat))
        this_lon = np.vstack((this_lon,temp_lon))
        this_ele = np.vstack((this_ele,temp_ele))
        
        this_cth_wwc = np.vstack((this_cth_wwc,temp_cth_wwc))
        this_sdcm_wwc = np.vstack((this_sdcm_wwc,temp_sdcm_wwc))
        
        this_cth = np.vstack((this_cth,temp_cth))
        this_sdcm = np.vstack((this_sdcm,temp_sdcm))
        
        this_cthom = np.vstack((this_cthom,temp_cthom))
        this_motqa = np.vstack((this_motqa,temp_motqa))
        this_wspd = np.vstack((this_wspd,temp_wspd))
        this_mdcm = np.vstack((this_mdcm,temp_mdcm))
        this_wele = np.vstack((this_wele,temp_wele))
    
# Get the dimensions

    hold = np.shape(this_lat)
    num_x = hold[0]
    num_y = hold[1]
    
# Rescale the SDCM
# Original: 0 = no retrieval, 1 = high confidence cloud, 2 = low confidence cloud,
#           3 = low confidence clear, 4 = high confidence clear
#      New: 0 = no retrieval/clear, 1 = cloud

    temp_sdcm_wwc = np.copy(this_sdcm_wwc)    

    this_sdcm_wwc = np.copy(temp_sdcm_wwc)
    this_sdcm_wwc[temp_sdcm_wwc == 2] = 1
    this_sdcm_wwc[temp_sdcm_wwc > 2] = 0
    
    temp_sdcm = np.copy(this_sdcm) 
    
    this_sdcm = np.copy(temp_sdcm)
    this_sdcm[temp_sdcm==2] = 1
    this_sdcm[temp_sdcm>2] = 0
    
    temp_mdcm = np.copy(this_mdcm) 
    
    this_mdcm = np.copy(temp_sdcm)
    this_mdcm[temp_mdcm==2] = 1
    this_mdcm[temp_mdcm>2] = 0

# Generate a new motion derived cloud mask

    wind_height_above_surf = this_cthom - this_wele
    this_wind_height_above_surf = np.copy(wind_height_above_surf)
    this_wind_height_above_surf[wind_height_above_surf <= 500.0] = 0
    this_wind_height_above_surf[wind_height_above_surf > 500.0] = 1
    
# Eliminate high wind speeds

    temp_wspd = np.copy(this_wspd) 
    this_wspd = np.copy(temp_wspd)
    this_wspd[temp_wspd>100.0] = 0.0
    
### SAVE THE DATA TO A NETCDF FILE

# Change to the figure directory

    os.chdir(datapath)

# Open the file
        
    n_id = Dataset(outfile,'w')
    
# Define the dimensions

    X_Dim = n_id.createDimension('X_Dim',num_x) # SOM-x
    Y_Dim = n_id.createDimension('Y_Dim',num_y) # SOM-y

# Define the output variables (2-D Mapping)
  
    out01 = n_id.createVariable('lon','f4',('X_Dim','Y_Dim'),zlib=True)
    out01.units = 'degrees_east'
    
    out02 = n_id.createVariable('lat','f4',('X_Dim','Y_Dim'),zlib=True)
    out02.units = 'degrees_north'

# Define the output variables

    out03 = n_id.createVariable('Elevation','f4',('X_Dim','Y_Dim'),zlib=True)
    
    out04 = n_id.createVariable('Cloud Top Height No Wind AGL','f4',('X_Dim','Y_Dim'),zlib=True)
    out05 = n_id.createVariable('Stereo Derived Cloud Mask No Wind','f4',('X_Dim','Y_Dim'),zlib=True)
    
    out06 = n_id.createVariable('Cloud Top Height AGL','f4',('X_Dim','Y_Dim'),zlib=True)
    out07 = n_id.createVariable('Stereo Derived Cloud Mask','f4',('X_Dim','Y_Dim'),zlib=True)
    
    out08 = n_id.createVariable('Height of Motion AGL','f4',('X_Dim','Y_Dim'),zlib=True)
    out09 = n_id.createVariable('Motion Quality','f4',('X_Dim','Y_Dim'),zlib=True)
    out10 = n_id.createVariable('Wind Speed','f4',('X_Dim','Y_Dim'),zlib=True)

    out12 = n_id.createVariable('Reduced Resolution Elevation','f4',('X_Dim','Y_Dim'),zlib=True)
    out13 = n_id.createVariable('Derived Motion Derived Cloud Mask','f4',('X_Dim','Y_Dim'),zlib=True)
  
# Put the data into the output (dimensions)

    out01[:] = this_lon
    out02[:] = this_lat
    
    out03[:] = this_ele
    
    out04[:] = this_cth_wwc - this_ele
    out05[:] = this_sdcm_wwc
    
    out06[:] = this_cth - this_ele
    out07[:] = this_sdcm
    
    out08[:] = this_cthom - this_wele
    out09[:] = this_motqa
    out10[:] = this_wspd
    
    out12[:] = this_wele
    out13[:] = this_wind_height_above_surf
    
# Close the NetCDF file
        
    n_id.close() 
    
# Print the time

    all_end_time = time.time()
    print("Total elapsed time was %g seconds" % (all_end_time - all_start_time))

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION


if __name__ == '__main__':
    main()