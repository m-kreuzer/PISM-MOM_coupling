#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:47:22 2022

@author: huiskamp
"""
# Create default runoff file for MOM6 that does not include fluxes around Antarctica.
# These are to be provided by PISM as part of the MOM-PISM-PICO coupling.

# The netCDF file we're altering is simply a copy of the defaul runoff for MOM6,
# i.e. runoff_newgrid.nc
# Similar to the default runoff, fluxes around Antarctica are time-invariant.

from netCDF4 import Dataset as CDF
import time

runoff_file = CDF('/p/projects/climber3/huiskamp/MOM6-examples/ice_ocean_SIS2/SIS2_coarse/new_grid/basal_test/melt_test_files/runoff_default.nc','r+')

runoff = runoff_file.variables['runoff'][:,:,:].data

for i in range(11):
    for j in range(runoff.shape[2]):
        if (runoff[:,i,j] != -9999).all() and (runoff[:,i,j] != 0).all():
            runoff[:,i,j] = 0
            
runoff_file.variables['runoff'][:] = runoff[:]

runoff_file.history = "This is the default runoff file for MOM6-PISM-PICO coupling created using the script runoff_default.py. \
    It was created on Created on " + time.ctime(time.time()) + " by Willem Huiskamp"
    
runoff_file.close()