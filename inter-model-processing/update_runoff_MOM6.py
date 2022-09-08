#  Copyright (C) 2019-2021 PISM-MOM_coupling authors, see AUTHORS file
#
#  This file is part of PISM-MOM_coupling
#
#  PISM-MOM_coupling is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PISM-MOM_coupling is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PISM-MOM_coupling.  If not, see <https://www.gnu.org/licenses/>.

""" Update the river runoff file for MOM6 with newly calculated values from PISM-PICO

usage: ./update_runoff_MOM6.py -i fluxes_file -o out_file [-v]

Runoff fields regridded from the PISM to the MOM grid are read in and applied 
to the model's default runoff file. This default runoff file should be specificially
created for use with PISM-PICO coupling and should have *no* runoff around the 
Antarctic margin.

Arguments:
    -i fluxes_file
        a netCDF file with variables to be processed
    -o out_file
       file to store processed fields which serve as PISM/PICO input 
    -v (optional)
        print verbose output        

This script was created as a processing tool for preparing output of the land ice model 
PISM-PICO as a forcing field for the ocean model MOM6-SIS2. This was done in the scope 
of coupling PISM to the climate model POEM at PIK.

"""
import sys
import argparse
from netCDF4 import Dataset as CDF
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                description=
                ("Update MOM6 runoff file with values from PISM-PICO.")
            )
    parser.add_argument('-i', '--PISM', 
                        action="store", 
                        dest="fluxes_file",
                        required=True, 
                        help="file with runoff data from PISM")
    parser.add_argument('-o', '--runoff', 
                        action="store", 
                        dest="out_file",
                        required=True, 
                        help="file with updated runoff for MOM6")
    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()
    
    if args.verbose:
        print("Running", sys.argv[0])
        print()
        
    # Open runoff from PISM
    if args.verbose:
        print("Opening PISM runoff file: " + args.fluxes_file)
    try:
        ds_ice = CDF(args.fluxes_file)
    except:
        print("fluxes_file '", args.fluxes_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    # Open runoff file
    if args.verbose:
        print("Opening MOM6 runoff file: " + args.out_file)
    try:
        ds_runoff = CDF(args.out_file)
    except:
        print("Runoff file can't be found! Exiting.")
        sys.exit(1)
        
    runoff_ice = ds_ice.variables['mass_flux_surf_runoff'][:,:,:].data
    runoff_new = ds_runoff.variables['runoff'][:,:,:].data
        
    # Apply ice runoff to all months in reference runoff file.
    for i in range(runoff_ice.shape[1]):
        for j in range(runoff_ice.shape[2]):
            if runoff_ice[0,i,j] != 0:
                runoff_new[:,i,j] = runoff_ice[0,i,j]

    ds_runoff.variables['runoff'][:] = runoff_new[:]

    ds_ice.close()
    ds_runoff.close()
    
    ##########################################################################
    #Test code#
    
    # from netCDF4 import Dataset as CDF
       
    # ds_ice = CDF('/p/tmp/kreuzer/coupled_PISM_MOM/experiments/MOM5_PISM_16km_spinup_run01/x_PISM-to-MOM/017860.fluxes.nc','r')
    # ds_runoff = CDF('/p/projects/climber3/huiskamp/POEM/work/PISM-MOM_coupling/runoff_newgrid.nc','r+')
        
    # runoff_ice = ds_ice.variables['mass_flux_surf_runoff'][:,:,:].data
    # runoff_new = ds_runoff.variables['runoff'][:,:,:].data

    
    # for i in range(runoff_ice.shape[1]):
    #     for j in range(runoff_ice.shape[2]):
    #         if runoff_ice[0,i,j] != 0:
    #             runoff_new[:,i,j] = runoff_ice[0,i,j]
                
    # ds_runoff.variables['runoff'][:] = runoff_new[:]
    
    # ds_ice.close()
    # ds_runoff.close()
                
    
        







