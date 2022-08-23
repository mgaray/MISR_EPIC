# MISR_AS_NetCDF_ghost_vis_tool_01.py
#
# This is a Python 2.7.10 code to read a MISR NetCDF AS_AEROSOL file, find 
# the associated L1B2 nadir file, and generate mapped images.
#
# Creation Date: 2016-02-19
# Last Modified: 2016-02-19
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

from __future__ import print_function # Makes 2.7 behave like 3.3
import fnmatch
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import os
from pyhdf.HDF import *
from pyhdf.SD import SD, SDC
from pyhdf.V import *
from pyhdf.VS import *
import time

import img_scale # Extra package (img_scale.py needs to be in execution directory)

def main():  # Main code

# Set the overall timer

    all_start_time = time.time()

# Set the pixel sizes for 1.1 km data

    x_size = 512
    y_size = 128

# Set the paths

    aeropath = '/Volumes/Makalu/DATA/POOPEX/2008_06_24/'
    agppath = '/Volumes/Makalu/DATA/POOPEX/2008_06_24/'
    imgpath = '/Volumes/Makalu/DATA/POOPEX/2008_06_24/'
    figpath = '/Users/mgaray/Desktop/CODING/PYTHON/PY27/FEB16/AEROSOL_VIS/FIGS/'

# Set the MISR product to test

    misr_name1 = '22' # 17.6 km Standard Product
    misr_name2 = '22b24-71' # 4.4 km Product

# Set the MISR orbit to test
    
    misr_path = 'P041'
    misr_orbit = '45310'    

# Set the block range

#    start_block = 60 # Start block 1-based
#    num_block = 3
    
#    start_block = 64 # Start block 1-based
#    num_block = 3
    
    start_block = 62 # Start block 1-based
    num_block = 3    

# Start the timer

    start_time = time.time()    

# Generate the output file name
    
    out_name = 'MISR_AOD_O'+misr_orbit+'_{}'.format(start_block)
    out_name = out_name+'_{}'.format(num_block)+'_07.png'

### Get the MISR AGP file for navigation

    os.chdir(agppath)
    search_str = 'MISR*'+misr_path+'*.hdf'
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

    lat_raw = var01.get()
    lon_raw = var02.get()

# Close the file

    hdf.end()

# Print the time

    end_time = time.time()
    print("Time to Read AGP data was %g seconds" % (end_time - start_time))
    
    
### Get the L1B2 Ellipsoid data for the An camera

# Start the timer

    start_time = time.time()

    os.chdir(imgpath)
    search_str = 'MISR*'+misr_path+'*'+misr_orbit+'*0024.hdf'
    file_list = glob.glob(search_str)

# Set the filename

    inputName = file_list[0]

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

# Print the time

    end_time = time.time()
    print("Time to Read MISR L1B2 data was %g seconds" % (end_time - start_time))


### Read the AEROSOL Data

# Start the timer

    start_time = time.time()

## Get the first MISR Aerosol File

# Change to the correct directory

#    basepath = aeropath+'V'+misr_name1+'/'
    basepath = aeropath

    os.chdir(basepath)

# Get the correct file

    search_str = 'MISR_AM1_AS_AEROSOL*'+misr_orbit+'*'+misr_name1+'.hdf'

    file_list = glob.glob(search_str)

# Set the filename

    inputName = file_list[0]

# Tell user location in process

    print("Reading: "+inputName)

# Open the file

    hdf = SD(inputName, SDC.READ)

# Read the data fields

    var01 = hdf.select('AlgTypeFlag')
    var02 = hdf.select('AerRetrSuccFlag')

    var03 = hdf.select('RegBestEstimateSpectralOptDepth')

    alg_type_01 = var01.get()
    succ_flag_01 = var02.get()

    rbe_aod_01 = var03.get()

# Close the file

    hdf.end()

## Get the second MISR Aerosol File

# Change to the correct directory

#    basepath = aeropath+'V'+misr_name2+'/'
    basepath = aeropath

    os.chdir(basepath)

# Get the correct file

    search_str = 'MISR_AM1_AS_AEROSOL*'+misr_orbit+'*'+misr_name2+'*.nc'

    file_list = glob.glob(search_str)

# Set the filename

    inputName = file_list[0]

# Tell user location in process

    print("Reading: "+inputName)
    
# Open the NetCDF file

    rootgrp = Dataset(inputName, 'r', format='NETCDF4')

# Choose the appropriate subgroup

    subgrp = rootgrp.groups['Regional']

# Assign the variables to an array

    alg_type_02 = subgrp.variables['Land_Water_Retrieval_Type_Flag'][:]
    succ_flag_02 = subgrp.variables['Retrieval_Success_Flag'][:]
    rbe_aod_02 = subgrp.variables['Aerosol_Optical_Depth'][:]
    clear_flag_02 = subgrp.variables['Clear_Flag_Fraction'][:]
    lat_02 = subgrp.variables['Latitude'][:]
    lon_02 = subgrp.variables['Longitude'][:]

# Close the NetCDF file

    rootgrp.close()

### Plot the image for multiple blocks

# Change to the figure directory

    os.chdir(figpath)

# Extract the navigation information to get corners for the entire image

    lat = lat_raw[start_block-1,:,:]
    lon = lon_raw[start_block-1,:,:]
    lat_max = np.amax(lat)
    lon_max = np.amax(lon)

    lat = lat_raw[start_block-1+num_block-1,:,:]
    lon = lon_raw[start_block-1+num_block-1,:,:]
    lat_min = np.amin(lat)
    lon_min = np.amin(lon)

# Set the plot area

    fig = plt.figure(figsize=(12,6), dpi=120)

### Set second subplot

    ax1 = fig.add_subplot(1,2,1)

## Set the plot

#    ax1.set_title('V22 (17.6 km)')
    ax1.set_title('17.6 km Resolution')

# Draw basemap

    m = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max,
        projection='cyl',resolution='i')
    m.drawmapboundary(fill_color='0.3')

## Loop over blocks

    for i in range(num_block):

        block = start_block + i - 1 # Block 0-based
        
# Extract the data

        lat = lat_raw[block,:,:]
        lon = lon_raw[block,:,:]
        alg_v22 = alg_type_01[block,:,:]
        succ_v22 = succ_flag_01[block,:,:]
        aod = rbe_aod_01[block,:,:,1]*1. # Green band
        
# Process the success flag as a mask (if successful set success to 1, otherwise 0)
# NOTE: Change to handle V22 format

        succ = np.copy(succ_v22)
        succ[succ_v22 != 7] = 0.0
        succ[succ_v22 == 7] = 1.0
        
# Process the algorithm type flag as a mask (if Dark Water set to 1, otherwise 0)

        dark_water = np.copy(alg_v22)
        dark_water[alg_v22 != 1] = 0.0
        dark_water[alg_v22 == 1] = 1.0
        
# Reset the success flag to require both success and Dark Water

        succ_water = succ*dark_water    

# Resize the data to match the navigation
# Note: The factor is 16 in both dimensions (17.6 -> 1.1 km)

        succ_full = np.repeat(np.repeat(succ,16,axis=0),16,axis=1)
        aod_full = np.repeat(np.repeat(aod,16,axis=0),16,axis=1)

# Plot data

        img = succ_full*aod_full # Only plot successful retrievals
        img_min = np.amin(img)
        img_max = np.amax(img)
        print(img_min)
        print(img_max)
        plot_min = 0.00
        plot_max = 1.0
        
        im = m.pcolormesh(lon,lat,img,shading='flat',cmap=plt.cm.CMRmap_r,
            latlon=True,vmin=plot_min,vmax=plot_max)

# Add the coastlines

    coast_color = 'green'
    m.drawcoastlines(color=coast_color,linewidth=0.8)
    
#    m.fillcontinents(color=coast_color)

# Add the colorbar

    aod_plot_ticks = 6 # Usually 1 more than you think
    aod_plot_step = 0.2
    cb = m.colorbar(im,"bottom", size="10%", pad="5%",
        ticks=np.arange(aod_plot_ticks)*aod_plot_step)
        
    cb.set_label('Green Band AOD')
    
### Set third subplot

    ax1 = fig.add_subplot(1,2,2)

## Set the plot

#    ax1.set_title('V'+misr_name2+' (4.4 km)')
    ax1.set_title('4.4 km Resolution')

# Draw basemap

    m = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max,
        projection='cyl',resolution='i')
    m.drawmapboundary(fill_color='0.3')

# Process the data

    lat_new = lat_02
    lon_new = lon_02

    succ_new = succ_flag_02
    alg_new = alg_type_02
    aod = rbe_aod_02*1.
    clear_flag = clear_flag_02*1.
        
# Process the success flag as a mask (if successful set success to 1, otherwise 0)

    succ = np.copy(succ_new)
    succ[succ_new != 1] = 0.0
    succ[succ_new == 1] = 1.0
    
# Process the algorithm type flag as a mask (if Dark Water set to 1, otherwise 0)

    dark_water = np.copy(alg_new)
    dark_water[alg_new != 0] = 0.0
    dark_water[alg_new == 0] = 1.0

# Reset the success flag to require both success and Dark Water

    succ_water = succ*dark_water    
        
# Process the clear flag fraction

#    thresh = 0.8
#    cff = np.copy(clear_flag)
#    cff[clear_flag > thresh] = 1.0
#    cff[clear_flag <= thresh] = 0.0      

# Plot data

    img = succ*aod.filled(0.0) # Only plot successful retrievals
    img_min = np.amin(img)
    img_max = np.amax(img)
    print(img_min)
    print(img_max)
    plot_min = 0.00
#    plot_max = 1.00
    
# Calculate the array index corresponding to the minimum latitude on the map

    idx1 = (np.abs(lat_new - lat_min)).argmin()
    idx1t = np.unravel_index(idx1,np.shape(lat_new))

# Calculate the array index corresponding to the maximum latitude on the map
    
    idx2 = (np.abs(lat_new - lat_max)).argmin()
    idx2t = np.unravel_index(idx2,np.shape(lat_new))
 
    max_x = idx1t[0]
    min_x = idx2t[0]
 
    lon_plot = lon_new[min_x:max_x,:]
    lat_plot = lat_new[min_x:max_x,:]
    img_plot = img[min_x:max_x,:]
    aod_plot = aod[min_x:max_x,:]
    
#    im = m.pcolormesh(lon_plot,lat_plot,img_plot,shading='flat',cmap=plt.cm.CMRmap_r,
#        latlon=True,vmin=plot_min,vmax=plot_max) 

# Reproject the map coordinates to handle an issue with pcolormesh   
        
    x, y = m(lon_plot,lat_plot)    
    im = m.pcolormesh(x,y,img_plot,shading='flat',cmap=plt.cm.CMRmap_r,
        vmin=plot_min,vmax=plot_max)     

# Add the coastlines

    m.drawcoastlines(color=coast_color,linewidth=0.8)
#    m.fillcontinents(color=coast_color)

# Add the colorbar

#    aod_plot_ticks = 6 # Usually 1 more than you think
#    aod_plot_step = 0.05
    cb = m.colorbar(im,"bottom", size="10%", pad="5%",
        ticks=np.arange(aod_plot_ticks)*aod_plot_step)
        
    cb.set_label('Green Band AOD')    

# Save the figure

    plt.savefig(out_name,dpi=300)  

# Show the figure

    plt.show()

# Print the time

    all_end_time = time.time()
    print("Total elapsed time was %g seconds" % (all_end_time - all_start_time))

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION


if __name__ == '__main__':
    main()    