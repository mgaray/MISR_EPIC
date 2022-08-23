# MISR_TERRAIN_HDF_2_NetCDF_Panoply_09.py
#
# This is a Python 3.9.13 code to read a MISR TERRAIN file and associated AGP file for
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
from pyhdf.V import *
from pyhdf.VS import *
import time

def main():  # Main code

# Print a blank line

    print()

# Set the overall timer

    all_start_time = time.time()  

# Set the paths
# NOTE: basepath = Location of the TERRAIN data
#       agppath = Location of the AGP data
#       datapath = Location of the output data file(s)

    basepath = '/Volumes/JPL_2022/DATA/MISR_FIRES/'
    agppath = '/Volumes/JPL_2022/DATA/AGP/'
    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY39/AUG22/TERRAIN/DATA/'

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

    outfile = 'MISR_AM1_TERRAIN_'+misr_path+'_O'+misr_orbit+'_Panoply_v'+vers+'.nc'

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

    agp_lat = var01.get()
    agp_lon = var02.get()

# Close the file

    hdf.end()
    
### GET THE MISR TERRAIN DATA (Df Camera)

# Change to the correct directory

    os.chdir(basepath)

    search_str = 'MISR_AM1_GRP_TERRAIN_GM_'+misr_path+'_O'+misr_orbit+'*DF*.hdf'
    dum_list = glob.glob(search_str)
    raw_list = np.array(dum_list)  # Convert to a numpy array

# Get the number of files
    
    num_files = len(raw_list)
    
    if(num_files < 1):
        print()
        print("***ERROR***")
        print("CANNOT FIND MISR TERRAIN FILE")
        print(search_str)
        print("***ERROR***")
        print()
        print(error)

# Set the filename

    inputName = raw_list[0]

# Tell user location in process

    print("Reading: "+inputName)
    
# Get the scale factors from the EOS Grid data
# Note: The reference numbers were found by parsing the full VDATA description
#       using the code MISR_L1B2_NAV_vis_tool_08.py

    f = HDF(inputName)
    vs = f.vstart()
    v = f.vgstart()

# BlueBand

    vg = vs.attach(223) # Scale factor
    blue_sf_raw = vg.read()
    blue_sf = np.squeeze(blue_sf_raw)
    vg.detach()

    vg = vs.attach(691) # std_solar_wgted_height
    blue_E0_raw = vg.read()
    vg.detach()

    vg = vs.attach(692) # SunDistanceAU
    MISR_AU_raw = vg.read()
    vg.detach()

# GreenBand

    vg = vs.attach(224) # Scale factor
    green_sf_raw = vg.read()
    green_sf = np.squeeze(green_sf_raw)
    vg.detach()

    vg = vs.attach(693) # std_solar_wgted_height
    green_E0_raw = vg.read()
    vg.detach()

# RedBand

    vg = vs.attach(225) # Scale factor
    red_sf_raw = vg.read()
    red_sf = np.squeeze(red_sf_raw)
    vg.detach()

    vg = vs.attach(695) # std_solar_wgted_height
    red_E0_raw = vg.read()
    vg.detach()

# NIRBand

    vg = vs.attach(226) # Scale factor
    nir_sf_raw = vg.read()
    nir_sf = np.squeeze(nir_sf_raw)
    vg.detach()

    vg = vs.attach(697) # std_solar_wgted_height
    nir_E0_raw = vg.read()
    vg.detach()

    v.end()
    vs.end()
    f.close()

# Open the file

    hdf = SD(inputName, SDC.READ)

# Read the data fields

    var01 = hdf.select('Blue Radiance/RDQI')
    var02 = hdf.select('Green Radiance/RDQI')
    var03 = hdf.select('Red Radiance/RDQI')
    var04 = hdf.select('NIR Radiance/RDQI')
    var05 = hdf.select('BlueConversionFactor')
    var06 = hdf.select('GreenConversionFactor')
    var07 = hdf.select('RedConversionFactor')
    var08 = hdf.select('NIRConversionFactor')
    var09 = hdf.select('SolarZenith')

    blue_raw = var01.get()
    green_raw = var02.get()
    red_raw = var03.get()
    nir_raw = var04.get()
    blue_cf_raw = var05.get()
    green_cf_raw = var06.get()
    red_cf_raw = var07.get()
    nir_cf_raw = var08.get()
    sza_raw = var09.get()

# Close the file

    hdf.end()
 
### EXTRACT MULTIPLE BLOCKS

    s_blk = start_block - 1 # 0-based
    num_blk = end_block - start_block + 1 # Count number of blocks +1 for loop control

# Extract the navigation data for the first block
    
    this_lat = agp_lat[s_blk,:,:]
    this_lon = agp_lon[s_blk,:,:]
    
# Extract the camera data for the first block

    blue = blue_raw[s_blk,:,:]
    green = green_raw[s_blk,:,:]
    red = red_raw[s_blk,:,:]
    nir = nir_raw[s_blk,:,:]
    
    blue_cf = blue_cf_raw[s_blk,:,:]
    green_cf = green_cf_raw[s_blk,:,:]
    red_cf = red_cf_raw[s_blk,:,:]
    nir_cf = nir_cf_raw[s_blk,:,:] 
    
# Strip off RDQI to get scaled DNs

    blue_dn = np.fix(blue/4.0)
    green_dn = np.fix(green/4.0)
    red_dn = np.fix(red/4.0)
    nir_dn = np.fix(nir/4.0)
    
# Mask the data outside the swath and with high RDQI
# Note: Based on the MISR documentation, it's not necessary to
#       calculate the RDQI independently since a value of 16380
#       indicates RDQI out of range.

    out_rng = 16378
        
    blue_dn[blue_dn >= out_rng] = 0.0
    green_dn[green_dn >= out_rng] = 0.0
    red_dn[red_dn >= out_rng] = 0.0
    nir_dn[nir_dn >= out_rng] = 0.0
    
# Convert to radiance

    blue_rad = blue_dn * blue_sf
    green_rad = green_dn * green_sf
    red_rad = red_dn * red_sf
    nir_rad = nir_dn * nir_sf
    
# Convert to BRF
# Note: The conversion factor has be resized by a factor of 64 in both
#       dimensions for the red band, but only 16 for the other bands

    blue_to_brf = np.repeat(np.repeat(blue_cf,16,axis=0),16,axis=1)
    green_to_brf = np.repeat(np.repeat(green_cf,16,axis=0),16,axis=1)
    red_to_brf = np.repeat(np.repeat(red_cf,64,axis=0),64,axis=1)
    nir_to_brf = np.repeat(np.repeat(nir_cf,16,axis=0),16,axis=1)
    
    blue_brf = blue_rad * blue_to_brf
    green_brf = green_rad * green_to_brf
    red_full = red_rad * red_to_brf
    nir_brf = nir_rad * nir_to_brf

# Store the results

    this_blue = np.copy(blue_brf)
    this_green = np.copy(green_brf)
    full_red = np.copy(red_full)
    this_nir = np.copy(nir_brf)

# Set missing data to zero (this is TERRAIN data)

    this_blue[blue_brf < 0.0] = 0.0
    this_green[green_brf < 0.0] = 0.0
    full_red[red_full < 0.0] = 0.0
    this_nir[nir_brf < 0.0] = 0.0
    
    this_blue[blue_brf > 1.0] = 0.0
    this_green[green_brf > 1.0] = 0.0
    full_red[red_full > 1.0] = 0.0
    this_nir[nir_brf > 1.0] = 0.0
    
# Resize the red band to match the others
# Note: We choose to use blue here as the reference size

    blue_y = blue_brf.shape[0]
    blue_x = blue_brf.shape[1]
    
    red_y = red_full.shape[0]
    red_x = red_full.shape[1]
    
    factor_y = int(red_y/blue_y)  # Conversion factor
    factor_x = int(red_x/blue_x)  # Conversion factor
    
    this_red = np.squeeze(full_red.reshape([blue_y,factor_y,
        blue_x,factor_x]).mean(3).mean(1))
    
# Loop through the blocks
    
    for loop in range(1,num_blk):
    
        this_blk = s_blk + loop
    
# Extract the navigation data

        temp_lat = agp_lat[this_blk,:,:]
        temp_lon = agp_lon[this_blk,:,:]
        
# Extract the camera data for the first block

        blue = blue_raw[this_blk,:,:]
        green = green_raw[this_blk,:,:]
        red = red_raw[this_blk,:,:]
        nir = nir_raw[this_blk,:,:]
    
        blue_cf = blue_cf_raw[this_blk,:,:]
        green_cf = green_cf_raw[this_blk,:,:]
        red_cf = red_cf_raw[this_blk,:,:]
        nir_cf = nir_cf_raw[this_blk,:,:] 
    
# Strip off RDQI to get scaled DNs

        blue_dn = np.fix(blue/4.0)
        green_dn = np.fix(green/4.0)
        red_dn = np.fix(red/4.0)
        nir_dn = np.fix(nir/4.0)
    
# Mask the data outside the swath and with high RDQI
# Note: Based on the MISR documentation, it's not necessary to
#       calculate the RDQI independently since a value of 16380
#       indicates RDQI out of range.

        out_rng = 16378
        
        blue_dn[blue_dn >= out_rng] = 0.0
        green_dn[green_dn >= out_rng] = 0.0
        red_dn[red_dn >= out_rng] = 0.0
        nir_dn[nir_dn >= out_rng] = 0.0
    
# Convert to radiance

        blue_rad = blue_dn * blue_sf
        green_rad = green_dn * green_sf
        red_rad = red_dn * red_sf
        nir_rad = nir_dn * nir_sf
    
# Convert to BRF
# Note: The conversion factor has be resized by a factor of 64 in both
#       dimensions for the red band, but only 16 for the other bands

        blue_to_brf = np.repeat(np.repeat(blue_cf,16,axis=0),16,axis=1)
        green_to_brf = np.repeat(np.repeat(green_cf,16,axis=0),16,axis=1)
        red_to_brf = np.repeat(np.repeat(red_cf,64,axis=0),64,axis=1)
        nir_to_brf = np.repeat(np.repeat(nir_cf,16,axis=0),16,axis=1)
    
        blue_brf = blue_rad * blue_to_brf
        green_brf = green_rad * green_to_brf
        red_full = red_rad * red_to_brf
        nir_brf = nir_rad * nir_to_brf

# Store the results

        temp_blue = np.copy(blue_brf)
        temp_green = np.copy(green_brf)
        full_red = np.copy(red_full)
        temp_nir = np.copy(nir_brf)
        
# Set missing data to zero (this is TERRAIN data)

        temp_blue[blue_brf < 0.0] = 0.0
        temp_green[green_brf < 0.0] = 0.0
        full_red[red_full < 0.0] = 0.0
        temp_nir[nir_brf < 0.0] = 0.0
        
        temp_blue[blue_brf > 1.0] = 0.0
        temp_green[green_brf > 1.0] = 0.0
        full_red[red_full > 1.0] = 0.0
        temp_nir[nir_brf > 1.0] = 0.0

# Resize the red band to match the others
# Note: We choose to use blue here as the reference size

        blue_y = blue_brf.shape[0]
        blue_x = blue_brf.shape[1]
    
        red_y = red_full.shape[0]
        red_x = red_full.shape[1]
    
        factor_y = int(red_y/blue_y)  # Conversion factor
        factor_x = int(red_x/blue_x)  # Conversion factor
    
        temp_red = np.squeeze(full_red.reshape([blue_y,factor_y,
            blue_x,factor_x]).mean(3).mean(1))

# Stack the data into a 2-D arrays

        this_lat = np.vstack((this_lat,temp_lat))
        this_lon = np.vstack((this_lon,temp_lon))
        
        this_blue = np.vstack((this_blue,temp_blue))
        this_green = np.vstack((this_green,temp_green))
        this_red = np.vstack((this_red,temp_red))
        this_nir = np.vstack((this_nir,temp_nir))
    
# Get the dimensions

    hold = np.shape(this_lat)
    num_x = hold[0]
    num_y = hold[1]
    
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

    out03 = n_id.createVariable('Blue_Band','f4',('X_Dim','Y_Dim'),zlib=True)
    out04 = n_id.createVariable('Green_Band','f4',('X_Dim','Y_Dim'),zlib=True)
    out05 = n_id.createVariable('Red_Band','f4',('X_Dim','Y_Dim'),zlib=True)
    out06 = n_id.createVariable('NIR_Band','f4',('X_Dim','Y_Dim'),zlib=True)
  
# Put the data into the output (dimensions)

    out01[:] = this_lon
    out02[:] = this_lat
    
    out03[:] = this_blue
    out04[:] = this_green
    out05[:] = this_red
    out06[:] = this_nir
    
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