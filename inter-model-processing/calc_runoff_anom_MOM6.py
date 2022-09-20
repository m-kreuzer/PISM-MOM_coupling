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

""" Calculate the anomalous runoff required for sea level change and 
return the value required for prcme_adj in the MOM6 data_table.

usage: ./update_runoff_MOM6.py -i fluxes_file -g grid_data [-v]

returns: Float of the global mean mass flux value to be applied to the ocean [kg/m2/s]

This script requires an experiment that is outputting grid information for both 
grid cell area and the ocean mask. This is typically the file ocean_static.nc

Arguments:
    -i fluxes_file
        a netCDF file with variables to be processed
    -g grid_data
    	A diagnostic file containing cell area and ocean mask.
    -v (optional)
        print verbose output        

This script was created as a processing tool for preparing output of the land ice model 
PISM-PICO as a forcing field for the ocean model MOM6-SIS2. This was done in the scope 
of coupling PISM to the climate model POEM at PIK.

"""

import sys
import argparse
from netCDF4 import Dataset as CDF
import numpy as np
import numpy.ma as ma
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                description=
                ("Calculate anomalous mass flux from PISM to MOM6")
            )
    parser.add_argument('-i', '--PISM', 
                        action="store", 
                        dest="fluxes_file",
                        required=True, 
                        help="file with runoff data from PISM")
    parser.add_argument('-g', '--grid', 
                        action="store", 
                        dest="grid_file",
                        required=True, 
                        help="file with static grid data for MOM6")
    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()
    
    if args.verbose:
        print("Running", sys.argv[0])
        print()
        
    # Open fluxes from PISM
    if args.verbose:
        print("Opening PISM runoff file: " + args.fluxes_file)
    try:
        ice_file = CDF(args.fluxes_file)
    except:
        print("fluxes_file '", args.fluxes_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    # Open grid file
    if args.verbose:
        print("Opening MOM6 grid file: " + args.out_file)
    try:
        grid_file = CDF(args.grid_file)
    except:
        print("Grid file can't be found! Exiting.")
        sys.exit(1)
        
    # Extract vars.    
    omask = grid_file.variables['wet'][:,:].data
    area  = grid_file.variables['area_t'][:,:].data

    mass_tot = ice_file.variables['mass_flux'][:,:].data


    # Calculate anomalous mass flux & total ocean area (I need extra info for surface mass balance I think?)
    area_m = ma.masked_array(area,mask=omask)
    area_tot = area_m.sum() 

    # Redistribute






    return 

    ds_ice.close()
    ds_grid.close()