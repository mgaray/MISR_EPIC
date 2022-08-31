# GOES_C_MISR_AOD_plot_01.py
#
# This is a Python 2.7.10 code to read the NetCDF file generated by the code
# MISR_AS_V23_GOES_C_AOD_match*.py and generate plots with overall statistics.
#
# Creation Date: 2020-02-21
# Last Modified: 2020-02-21
#
# by Michael J. Garay
# (Michael.J.Garay@jpl.nasa.gov)

# Import packages

import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os

def main():  # Main code

# Set the acceptable AERONET standard deviation

    acc_std = 1000.0
    acc_chi = 1000.0

# Set the minimum and maximum AOD for plotting

    aod_plot_min = 0.00
    
#    aod_plot_max = 0.30
#    aod_plot_max = 0.60
#    aod_plot_max = 1.20
    aod_plot_max = 1.40
    
    log_plot_ref = 1.00
    
#    aod_plot_ticks = 7 # Usually 1 more than you think
    aod_plot_ticks = 8 # Usually 1 more than you think
    
#    aod_plot_step = 0.05
#    aod_plot_step = 0.1
    aod_plot_step = 0.2

# Set the paths

    datapath = '/Users/mgaray/Desktop/CODING/PYTHON/PY36/FEB20/FIREX/DATA/'
    figpath = '/Users/mgaray/Desktop/CODING/PYTHON/PY36/FEB20/FIREX/FIGS/'

# Set the GOES product

    misr_name = 'OR_ABI-L2-AODC' # CONUS Product

# Set the filename

    data_name = "MISR_GOES_C_AOD_2019_08_06_v01.nc"
    
# Set the base file name

    out_base = "_"+misr_name+"_2019_08_06_raw_stats_01.png"

### Find and read the correct NetCDF file

# Change directory to the basepath and get the file list

    os.chdir(datapath)
    file_list = glob.glob(data_name)

# Choose the first file

    inputName = file_list[0]

# Tell user location in process

    print('Reading: ',inputName)

# Open the NetCDF file

    rootgrp = Dataset(inputName, 'r', format='NETCDF4')

# Assign the variables to arrays
               
    goes_lon = rootgrp.variables['GOES_Longitude'][:]
    goes_lat = rootgrp.variables['GOES_Latitude'][:]
    dist = rootgrp.variables['Distance'][:]
                
    goes_aod = rootgrp.variables['GOES_AOD_550nm'][:]
    goes_dqf = rootgrp.variables['GOES_DQF'][:]
    
    misr_lon = rootgrp.variables['MISR_Longitude'][:]
    misr_lat = rootgrp.variables['MISR_Latitude'][:]
                
    misr_aod = rootgrp.variables['MISR_AOD_550nm'][:]
    misr_raw = rootgrp.variables['MISR_AOD_RAW_550nm'][:]
    misr_lwr = rootgrp.variables['Land_Water_Retrieval_Type'][:]
   
# Close the NetCDF file

    rootgrp.close()
    
### Set the plotting values

#    test_best_raw = misr_aod_v22
#    test_low_raw = misr_low_aod_v22
#    ref_aod_raw = aero_aod_558_close
#    ref_std = aero_std_558
    
### FIRST PLOT: LINEAR REGRESSION

    max_val = aod_plot_max

# Set the plot area
# NOTE: The base plot size is 6 x 6, so a 2 row, 3 column set would be 18 x 12

    plt.figure(figsize=(12,6), dpi=120)

## Linear plot (Best Estimate)

#    good = ((goes_dqf == 0) & (misr_aod > 0))
    good = ((goes_dqf == 0) & (misr_raw > 0))
#    ref_aod = misr_aod[good]
    ref_aod = misr_raw[good]
    test_aod = goes_aod[good]

    plt.subplot(1, 2, 1)
    plt.scatter(ref_aod,test_aod,marker='o',color='black',s=10)   
    plt.title(misr_name+" High Quality")

# Plot the one-to-one line

    plt.plot([0.0,max_val], [0.0,max_val], color="k", lw=1)

# Plot the envelopes

    dummy_aod = np.logspace(-4,1,num=100)
    up1_aod = 1.20*dummy_aod
    up2_aod = dummy_aod+0.05
    upper_aod = np.maximum(up1_aod,up2_aod)

    lo1_aod = 0.80*dummy_aod
    lo2_aod = dummy_aod-0.05
    lower_aod = np.minimum(lo1_aod,lo2_aod)

    plt.plot(dummy_aod,lower_aod,color="0.75", lw=1)
    plt.plot(dummy_aod,upper_aod,color="0.75", lw=1)

# Set the limits and axis labels

    plt.xlim(0.0,max_val)
    plt.ylim(0.0,max_val)

    plt.xlabel('MISR AOD')
    plt.ylabel('GOES AOD')

    plt.grid(True)

# Include some text on the Best Estimate Figure

    x_pos = (aod_plot_max/2.) + (aod_plot_max/10.)
    y_pos1 = (aod_plot_max/10.)
    y_pos2 = y_pos1 + (aod_plot_max/30.)
    y_pos3 = y_pos2 + (aod_plot_max/30.)
    y_pos4 = y_pos3 + (aod_plot_max/30.)
    y_pos5 = y_pos4 + (aod_plot_max/30.)
    y_pos6 = y_pos5 + (aod_plot_max/30.)
    
    plt.text(x_pos,y_pos6,'High Quality',fontsize=12) # Version

    count = len(test_aod)
    out_text = 'N = '+str(count)
    plt.text(x_pos,y_pos5,out_text,fontsize=10) # Count

    temp = np.corrcoef(ref_aod,test_aod)
    be_r = temp[0,1]
    out_text = 'r = '+"{0:.4f}".format(be_r)
    plt.text(x_pos,y_pos4,out_text,fontsize=10) # Correlation coefficient

    rmse = np.sqrt(((test_aod - ref_aod) ** 2).mean())
    out_text = 'RMSE = '+"{0:.4f}".format(rmse)
    plt.text(x_pos,y_pos3,out_text,fontsize=10) # Root mean squared error

    diff = test_aod - ref_aod
    bias = np.mean(diff)
    out_text = 'Bias = '+"{0:.4f}".format(bias)
    plt.text(x_pos,y_pos2,out_text,fontsize=10) # Bias

    offset = np.ones_like(ref_aod)*0.05
    inner = np.absolute(diff) < np.maximum(offset,ref_aod*0.2)
    in_frac = (np.sum(inner)/(1.0*count))*100.0
    out_text = 'Percent In = '+"{0:.2f}".format(in_frac)
    plt.text(x_pos,y_pos1,out_text,fontsize=10) # Percent in envelope

## Linear plot (Lowest Residual)

#    good = ((goes_dqf <= 1) & (misr_aod > 0))
#    ref_aod = misr_aod[good]
    good = ((goes_dqf <= 1) & (misr_raw > 0))
    ref_aod = misr_raw[good]
    test_aod = goes_aod[good]

    plt.subplot(1, 2, 2)
    plt.scatter(ref_aod,test_aod,marker='o',color='black',s=10)   
    plt.title(misr_name+" Medium-High Quality")

# Plot the one-to-one line

    plt.plot([0.0,max_val], [0.0,max_val], color="k", lw=1)

# Plot the envelopes

    dummy_aod = np.logspace(-4,1,num=100)
    up1_aod = 1.20*dummy_aod
    up2_aod = dummy_aod+0.05
    upper_aod = np.maximum(up1_aod,up2_aod)

    lo1_aod = 0.80*dummy_aod
    lo2_aod = dummy_aod-0.05
    lower_aod = np.minimum(lo1_aod,lo2_aod)

    plt.plot(dummy_aod,lower_aod,color="0.75", lw=1)
    plt.plot(dummy_aod,upper_aod,color="0.75", lw=1)

# Set the limits and axis labels

    plt.xlim(0.0,max_val)
    plt.ylim(0.0,max_val)

    plt.xlabel('MISR AOD')
    plt.ylabel('GOES AOD')

    plt.grid(True)

# Include some text on the Lowest Residual Figure

    plt.text(x_pos,y_pos6,'Medium-High',fontsize=12) # Version

    count = len(test_aod)
    out_text = 'N = '+str(count)
    plt.text(x_pos,y_pos5,out_text,fontsize=10) # Count

    temp = np.corrcoef(ref_aod,test_aod)
    be_r = temp[0,1]
    out_text = 'r = '+"{0:.4f}".format(be_r)
    plt.text(x_pos,y_pos4,out_text,fontsize=10) # Correlation coefficient

    rmse = np.sqrt(((test_aod - ref_aod) ** 2).mean())
    out_text = 'RMSE = '+"{0:.4f}".format(rmse)
    plt.text(x_pos,y_pos3,out_text,fontsize=10) # Root mean squared error

    diff = test_aod - ref_aod
    bias = np.mean(diff)
    out_text = 'Bias = '+"{0:.4f}".format(bias)
    plt.text(x_pos,y_pos2,out_text,fontsize=10) # Bias

    offset = np.ones_like(ref_aod)*0.05
    inner = np.absolute(diff) < np.maximum(offset,ref_aod*0.2)
    in_frac = (np.sum(inner)/(1.0*count))*100.0
    out_text = 'Percent In = '+"{0:.2f}".format(in_frac)
    plt.text(x_pos,y_pos1,out_text,fontsize=10) # Percent in envelope

# Save the figure

    os.chdir(figpath)
    outname = 'AOD_Mean_Line_Regression'+out_base
    plt.savefig(outname,dpi=120)

### SECOND PLOT: LOG REGRESSION

# Set the plot area
# NOTE: The base plot size is 6 x 6, so a 2 row, 3 column set would be 18 x 12

    plt.figure(figsize=(12,6), dpi=120)

## Log plot (Best Estimate)

#    good = ((goes_dqf == 0) & (misr_aod > 0))
#    ref_aod = misr_aod[good]
    good = ((goes_dqf == 0) & (misr_raw > 0))
    ref_aod = misr_raw[good]
    test_aod = goes_aod[good]

    plt.subplot(1, 2, 1)
    plt.scatter(ref_aod,test_aod,marker='o',color='black',s=10)   
    plt.title(misr_name+" High Quality")

# Plot the one-to-one line

    plt.plot([0.001,10.0], [0.001,10.0], color="k", lw=1)

# Plot the envelopes

    dummy_aod = np.logspace(-4,1,num=100)
    up1_aod = 1.20*dummy_aod
    up2_aod = dummy_aod+0.05
    upper_aod = np.maximum(up1_aod,up2_aod)

    lo1_aod = 0.80*dummy_aod
    lo2_aod = dummy_aod-0.05
    lower_aod = np.minimum(lo1_aod,lo2_aod)

    plt.plot(dummy_aod,lower_aod,color="0.75", lw=1)
    plt.plot(dummy_aod,upper_aod,color="0.75", lw=1)

# Set the limits and axis labels

    plt.xlim(0.001,10.0)
    plt.ylim(0.001,10.0)

    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel('MISR AOD')
    plt.ylabel('GOES AOD')

    plt.grid(True)

# Include some text on the Best Estimate Figure

    x_pos = (log_plot_ref/3.) + (log_plot_ref/30.)
    y_pos1 = (log_plot_ref/300.)
    y_pos2 = y_pos1 + (log_plot_ref/800.)
    y_pos3 = y_pos2 + (log_plot_ref/600.)
    y_pos4 = y_pos3 + (log_plot_ref/450.)
    y_pos5 = y_pos4 + (log_plot_ref/400.)
    y_pos6 = y_pos5 + (log_plot_ref/250.)
    
    plt.text(x_pos,y_pos6,'High Quality',fontsize=12) # Version

    count = len(test_aod)
    out_text = 'N = '+str(count)
    plt.text(x_pos,y_pos5,out_text,fontsize=10) # Count

    temp = np.corrcoef(ref_aod,test_aod)
    be_r = temp[0,1]
    out_text = 'r = '+"{0:.4f}".format(be_r)
    plt.text(x_pos,y_pos4,out_text,fontsize=10) # Correlation coefficient

    rmse = np.sqrt(((test_aod - ref_aod) ** 2).mean())
    out_text = 'RMSE = '+"{0:.4f}".format(rmse)
    plt.text(x_pos,y_pos3,out_text,fontsize=10) # Root mean squared error

    diff = test_aod - ref_aod
    bias = np.mean(diff)
    out_text = 'Bias = '+"{0:.4f}".format(bias)
    plt.text(x_pos,y_pos2,out_text,fontsize=10) # Bias

    offset = np.ones_like(ref_aod)*0.05
    inner = np.absolute(diff) < np.maximum(offset,ref_aod*0.2)
    in_frac = (np.sum(inner)/(1.0*count))*100.0
    out_text = 'Percent In = '+"{0:.2f}".format(in_frac)
    plt.text(x_pos,y_pos1,out_text,fontsize=10) # Percent in envelope

## Log plot (Lowest Residual)

#    good = ((goes_dqf <= 1) & (misr_aod > 0))
#    ref_aod = misr_aod[good]
    good = ((goes_dqf <= 1) & (misr_raw > 0))
    ref_aod = misr_raw[good]
    test_aod = goes_aod[good]

    plt.subplot(1, 2, 2)
    plt.scatter(ref_aod,test_aod,marker='o',color='black',s=10)   
    plt.title(misr_name+" Medium-High Quality")

# Plot the one-to-one line

    plt.plot([0.001,10.0], [0.001,10.0], color="k", lw=1)

# Plot the envelopes

    dummy_aod = np.logspace(-4,1,num=100)
    up1_aod = 1.20*dummy_aod
    up2_aod = dummy_aod+0.05
    upper_aod = np.maximum(up1_aod,up2_aod)

    lo1_aod = 0.80*dummy_aod
    lo2_aod = dummy_aod-0.05
    lower_aod = np.minimum(lo1_aod,lo2_aod)

    plt.plot(dummy_aod,lower_aod,color="0.75", lw=1)
    plt.plot(dummy_aod,upper_aod,color="0.75", lw=1)

# Set the limits and axis labels

    plt.xlim(0.001,10.0)
    plt.ylim(0.001,10.0)

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('MISR AOD')
    plt.ylabel('GOES AOD')

    plt.grid(True)

# Include some text on the Lowest Residual Figure

    plt.text(x_pos,y_pos6,'Medium-High',fontsize=12) # Version

    count = len(test_aod)
    out_text = 'N = '+str(count)
    plt.text(x_pos,y_pos5,out_text,fontsize=10) # Count

    temp = np.corrcoef(ref_aod,test_aod)
    be_r = temp[0,1]
    out_text = 'r = '+"{0:.4f}".format(be_r)
    plt.text(x_pos,y_pos4,out_text,fontsize=10) # Correlation coefficient

    rmse = np.sqrt(((test_aod - ref_aod) ** 2).mean())
    out_text = 'RMSE = '+"{0:.4f}".format(rmse)
    plt.text(x_pos,y_pos3,out_text,fontsize=10) # Root mean squared error

    diff = test_aod - ref_aod
    bias = np.mean(diff)
    out_text = 'Bias = '+"{0:.4f}".format(bias)
    plt.text(x_pos,y_pos2,out_text,fontsize=10) # Bias

    offset = np.ones_like(ref_aod)*0.05
    inner = np.absolute(diff) < np.maximum(offset,ref_aod*0.2)
    in_frac = (np.sum(inner)/(1.0*count))*100.0
    out_text = 'Percent In = '+"{0:.2f}".format(in_frac)
    plt.text(x_pos,y_pos1,out_text,fontsize=10) # Percent in envelope

# Save the figure

    os.chdir(figpath)
    outname = 'AOD_Mean_Log_Regression'+out_base
    plt.savefig(outname,dpi=120)

# Show the plot

    plt.show()

# Tell user completion was successful

    print("\nSuccessful Completion\n")

### END MAIN FUNCTION


if __name__ == '__main__':
    main()