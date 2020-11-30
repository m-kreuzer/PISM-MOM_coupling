#!/usr/bin/env python3

#  Copyright (C) 2019, 2020 PISM-MOM_coupling authors, see AUTHORS file
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

""" Aggregate PICO/PISM input fields from regridded MOM output.

usage: ./regriddedMOM-to-PISM_processing.py -i field_file -b basin_file \
            -e edge_file -f var1 [var2 ...] [-d basin_shelf_depth_file] \
            -o out_file [-t] [-v]

Regridded variables from ocean model MOM5 to landice model PISM are processed 
and aggregated to be used as 2d boundary conditions for the PICO component (PISM).
Empty grid cells of variables are filled with mean of all (precomputed) edge 
values [missing/non-missing] located in the same PICO basin. Averaging 
and filling is done for each basin in every layer and each variable. Then 
vertical interpolation is done for every grid point to the mean bedrock depth 
underneath the ice shelf for each basin. These values are read from file or 
use a default depth if none is given.

Arguments:
    -i field_file
        a netCDF file with variables to be processed
    -b basin_file
        a netCDF file with variable 'basins' describing the basins
        of the PICO model as well as the variable 'pico_contshelf_mask'
        with the continental shelf mask
    -e edge_file
        a netCDF file with variable 'field_edge' which masks all
        grid cells with empty neighbouring cells. This file was precomputed
        by the script 'find_edge.py'
    -f var1 [var2 ...]
        list of variable names whose missing cells should be filled
    -d basin_shelf_depth_file (optional)
        file with variable 'mean_shelf_topg' which stores basin shelf depths 
        to determine vertical layer of ocean boundary condition input to 
        PISM/PICO. Gets computed by script 'PISM-to-MOM_processing.py'.
    -o out_file
       file to store processed fields which serve as PISM/PICO input 
    -t (optional)
        print script time statistics
    -v (optional)
        print verbose output        


This script requires the ouput of script find_edge.py and optionally output of
PISM-to-MOM_processing.py

This script was created as a processing tool for preparing output of ocean 
model MOM5 as boundary condition input to the landice model PISM/PICO. This was 
done in the scope of coupling PISM to the climate model POEM at PIK.

"""

import sys
import os
import time as t
import copy as cp
import argparse
import numpy as np
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


if __name__ == "__main__":
   

    parser = argparse.ArgumentParser(
                description=
                ("Fill empty grid cells of netCDF variables with appropriate "
                 "basin edge mean values."),
                epilog=
                ("Empty grid cells of variables in VAR_FILL are filled with "
                 "mean of all values located in the same PICO basin (read from "
                 "BASIN_FILE) and which have empty grid cells in their direct "
                 "neighbourhood (4-point stencil). These cells are the so "
                 "called 'edge' of missing/non-missing values and are read "
                 "from EDGE_FILE which is precomputed by 'find_edge.py'. "
                 "Averaging and filling is done for each basin in every layer "
                 "and each variable.")
            )

    parser.add_argument('-i', '--input', 
                        action="store", 
                        dest="field_file",
                        required=True, 
                        help="file with variables to be processed")
    parser.add_argument('-e', '--edges', 
                        action="store", 
                        dest="edge_file",
                        required=True, 
                        help="file storing precalculated edge cells of        \
                                FIELD_FILE")
    parser.add_argument('-b', '--basins', 
                        action="store", 
                        dest="basin_file",
                        required=True, 
                        help="file with variable 'basins' from PICO model")
    parser.add_argument('-f', '--fill', 
                        action="store", 
                        dest="var_fill",
                        required=True, 
                        nargs='+',
                        help="list of variable names whose missing cells     \
                                should be filled")
    parser.add_argument('-d', '--depth',
                        action="store",
                        dest="basin_shelf_depth_file",
                        required=False,
                        help=("file with variable 'mean_shelf_topg' which "
                              "stores basin shelf depths to determine vertical "
                              "layer of ocean boundary condition input to "
                              "PISM/PICO. Gets computed by script "
                              "'PISM_postprocessing.py'."))
    parser.add_argument('-o', '--output', 
                        action="store", 
                        dest="out_file",
                        required=True, 
                        help="file to store processed fields")
    parser.add_argument('-t', '--time', 
                        action="store_true", 
                        help="print script timings")
    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()



    # -------------------- general setup --------------------  
    t_main_start = t.time()
     
    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()
       
    # preparing data structure to manage fields to be filled
    d = dict({'field_in':None, 'field_out':None , 'field_out_basin_mean': None,
                'dim':None})
    fields = dict({f: cp.deepcopy(d)  for f in args.var_fill})
    

    # -------------------- read basin mask --------------------  
    if args.verbose:
        print("... reading basin mask from file '" + args.basin_file + "'")
    t_read_start = t.time()
    t_read_bmask_start = t.time()
    
    try:
        nc_b = CDF(args.basin_file, 'r')
    except:
        print("BASIN_FILE '", args.basin_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    # read basin array    
    try:
        basins = np.squeeze(nc_b.variables['basins'][:])
        basins_dtype = nc_b.variables['basins'].datatype
        basins_dim   = nc_b.variables['basins'].dimensions
        basins_dict  = nc_b.variables['basins'].__dict__
    except:
        print("Variable 'basins' can't be read from BASIN_FILE '" 
                + args.basin_file + "'!")
    try:
        # mask values:
        #   0 = grounded ice or open ocean
        #   1 = continental shelf but ignored
        #   2 = continental shelf as used in PICO
        contshelf_mask = np.squeeze(nc_b.variables['pico_contshelf_mask'][:])
        contshelf_mask_dtype = nc_b.variables['pico_contshelf_mask'].datatype
        contshelf_mask_dim   = nc_b.variables['pico_contshelf_mask'].dimensions
        contshelf_mask_dict  = nc_b.variables['pico_contshelf_mask'].__dict__
    except:
        print("Variable 'pico_contshelf_mask' can't be read from BASIN_FILE '" 
                + args.basin_file + "'!")

    # check dimensions
    basin_ndim = len(basins.shape)
    if basin_ndim != 2 and basin_ndim != 3:
        raise ValueError( str("basin field is of dimension " + \
                            str( basin_ndim ) + ". Expected: 2 or 3.") )
    if basin_ndim == 3:
        # cut of time dimension and take first time slice
        basins = basins[0,:,:]

    contshelf_ndim = len(contshelf_mask.shape)
    if contshelf_ndim != 2 and contshelf_ndim != 3:
        raise ValueError( str("pico continental shelf mask field is of dimension " + \
                            str( contshelf_ndim ) + ". Expected: 2 or 3.") )
    if contshelf_ndim == 3:
        # cut of time dimension and take first time slice
        contshelf_mask = contshelf_mask[0,:,:]

    assert_str = ("non matching dimensions of variables 'basins' and "
                  "'pico_contshelf_mask' in BASIN_FILE '{}'.")
    assert basins.shape == contshelf_mask.shape, \
                assert_str.format(args.basin_file)
        
    # create list with all occuring basins
    basin_vals = np.unique(basins)
    basin_list = basin_vals[~basin_vals.mask].data
    # remove basin 0
    basin_list = np.delete(basin_list, np.where(basin_list==0) ) 
    # make sure datatype is integer
    basin_list = basin_list.astype(int)


    nc_b.close()
    t_read_bmask_end = t.time()

    # -------------------- read edge mask -------------------- 
    t_read_edge_start = t.time()

    if args.verbose:
        print("... reading edge mask from file '" + args.edge_file + "'")

    try:
        nc_e = CDF(args.edge_file, 'r')
    except:
        print("EDGE_FILE '", args.edge_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    # read edge array    
    #  -> important to read in as bool, otherwise type is int >> bad performance
    try:
        field_edge = nc_e.variables['field_edge'][:].astype(bool)
    except:
        print("Variable 'field_edge' can't be read from EDGE_FILE '" 
                + args.edge_file + "'!")
    
    # check basin dimension
    field_edge_ndim = len(field_edge.shape)
    if field_edge_ndim != 3:
        raise ValueError( str("edge field is of dimension " + \
                            str( field_edge_ndim ) + ". Expected: 3.") )
    
    nc_e.close()
    t_read_edge_end = t.time()


    # -------------------- read variables -------------------- 
    #   -> read file with regridded variables to be processed
    t_read_infile_start = t.time()
    
    if args.verbose:
        print("... reading VAR_FILL fields '" + str(args.var_fill) + 
                "' from FIELD_FILE '" + args.field_file + "'")

    try:
        nc_src = CDF(args.field_file, 'r')
    except:
        print("FIELD_FILE '" + args.field_file + "' can't be found! Exiting.")
        sys.exit(1)
        
    # extract vertical axis
    ocean_z = nc_src.variables['st_ocean'][:]   # units: dbars (interpreting as m)
    ocean_z = -1 * ocean_z                      # positive = upwards
    
    # read field arrays 
    for f in fields.keys():
        try:
            fields[f]['field_in'] = np.squeeze( nc_src.variables[f][:] )
            #fields[f]['field_use'] = np.squeeze( nc_src.variables[f][:] )
        except:
            err_str = "Variable '" + f + "' can't be read from FIELD_FILE '" \
                        + args.field_file + "'!"
            raise ValueError( str(err_str) )
        fields[f]['dim'] = nc_src.variables[f].dimensions
        field_ndim = len(fields[f]['field_in'].shape)

        if ( field_ndim != 3) :
            if (field_ndim == 4 and fields[f]['dim'][0] == 'time'):
                # select last timeslice
                fields[f]['field_in'] = np.squeeze(fields[f]['field_in'][-1,:,:,:])
            else:
                err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
                        args.field_file + "' is of dimension " + \
                        str(field_ndim) + ". Expected: 3 (z,y,x) or 4 (time,z,y,x)."
                raise ValueError( str(err_str) )

        fields[f]['field_out'] = np.zeros(shape=fields[f]['field_in'].shape[-2:])
        fields[f]['field_out'][:] = np.nan
        fields[f]['field_out_basin_mean'] = np.zeros(shape=fields[f]['field_in'].shape[-2:])
        fields[f]['field_out_basin_mean'][:] = np.nan
        print(fields[f]['field_out'].shape)

    nc_src.close()
    t_read_infile_end = t.time()
   
    

    # -------------------- read basin depths -------------------- 
    #   -> read file with regridded variables to be processed
    
    # input file for shelf depth is given -> read it
    if args.basin_shelf_depth_file is not None:
        
        t_read_depthfile_start = t.time()
        if args.verbose:
            print("... reading basin shelf topography from BASIN_SHELF_DEPTH_FILE '" 
                  + args.basin_shelf_depth_file + "'")
    
        try:
            nc_src = CDF(args.basin_shelf_depth_file, 'r')
        except:
            print("BASIN_SHELF_DEPTH_FILE '" + args.basin_shelf_depth_file + 
                  "' can't be found! Exiting.")
            sys.exit(1)
        
        n_basin_name = nc_src.dimensions['n_basin'].name
        n_basin_size = nc_src.dimensions['n_basin'].size
        
        shelf_depth_basin_list  = nc_src.variables['basin'][:]
        shelf_depth_basin_dtype = nc_src.variables['basin'].datatype
        shelf_depth_basin_dim   = nc_src.variables['basin'].dimensions
        shelf_depth_basin_dict  = nc_src.variables['basin'].__dict__
        
        mean_shelf_topg         = np.squeeze(nc_src.variables['mean_shelf_topg'][:])
        mean_shelf_topg_dtype   = nc_src.variables['mean_shelf_topg'].datatype
        mean_shelf_topg_dim     = nc_src.variables['mean_shelf_topg'].dimensions
        mean_shelf_topg_dict    = nc_src.variables['mean_shelf_topg'].__dict__

        nc_src.close()
        t_read_depthfile_end = t.time()
        
    else:        
        # set default values

        if args.verbose:
            print("... no BASIN_SHELF_DEPTH_FILE given! \n"
                  "\t -> using default depth -500m for all basins")
    
        n_basin_name = 'n_basin'
        n_basin_size = len(basin_vals)
        
        shelf_depth_basin_list  = basin_list
        #shelf_depth_basin_list  = np.array(np.arange(1,20))  #temorary fix
        shelf_depth_basin_dtype = int
        shelf_depth_basin_dim   = ('n_basin',)
        shelf_depth_basin_dict  = \
            dict({'long_name':  "list of valid PISM/PICO basins"})
        
        mean_shelf_topg         = np.ones_like(shelf_depth_basin_list) * -500
        mean_shelf_topg_dtype   = float
        mean_shelf_topg_dim     = ('n_basin',)
        mean_shelf_topg_dict    = \
            dict({'long_name':  "mean basin topography of ice shelf areas", \
                  'units':      "m", \
                  'axis':       "Z", \
                  'positive':   "up"})        
    
    t_read_end = t.time()


    ###  -------------------- fill missing field values -------------------- 
    #  -> empty cells are filled with the mean of all values on the edge 
    #       from missing to non-missing values. This is done for each basin
    #       in each layer in each variable/field.

    t_fill_start = t.time()

    if args.verbose:
        print('... fill empty grid cells')

    # iterate through fields to be filled
    for f in fields.keys():
        if args.verbose: 
            print('\t > ', f)
        
        # check field dimension
        field = fields[f]['field_in']
        field_shape = field.shape
        depth_levels = field_shape[0]

        if (field_shape[-2:] != basins.shape[-2:]) :
            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
                        args.field_file + "' has shape " + \
                        str(field_shape[-2:]) + " which does not match the " + \
                        "shape of 'basins' " + str(basins.shape) + \
                        " of BASIN_FILE '" + args.basin_file + "'!" 
            raise ValueError( str(err_str) ) 
        
        if (field_shape[-2:] != field_edge.shape[-2:]) :
            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
                        args.field_file + "' has shape " + \
                        str(field_shape[-2:]) + " which does not match the " + \
                        "shape of 'field_edge' " + str(field_edge.shape[-2:]) + \
                        " of EDGE_FILE '" + args.edge_file + "'!" 
            raise ValueError( str(err_str) ) 
        
        # calculate basin mean of edge cells and apply for empty cells
        for z in range(depth_levels):
            for b in basin_vals:

                ### create masks
                # mask of current basin
                m__basin = (basins==b)
                # mask of field edge AND current basin
                m__fedge_A_basin = field_edge[z,:] & m__basin
                # mask of missing fields AND current basin
                m__fmiss_A_basin = (field[z,:].mask==True) & m__basin

                ### calculate mean and write to array 
                mean_field_basin_edge = field[z,m__fedge_A_basin].mean()
                fields[f]['field_in'][z,m__fmiss_A_basin] = mean_field_basin_edge

    t_fill_end = t.time()

    ### --------------- depth interpolation ---------------------------------
    t_interp_start = t.time()    

    if args.verbose:
        print('... interpolate basin shelf depths')
    
    # check whether basin lists are matching
    assert_str = ("non matching basins between BASIN_FILE '{}' and "
                  "BASIN_SHELF_DEPTH_FILE '{}'.")
    assert set(basin_list) == set(shelf_depth_basin_list.data), \
                assert_str.format(args.basin_file, args.basin_shelf_depth_file)

    # iterate basins 
    for b_idx, b_val in enumerate(shelf_depth_basin_list.data):
        if args.verbose:
            print('\t > basin ', b_val, ' / ', shelf_depth_basin_list.max())
            
        # depth of current basin: mean_shelf_topg[idx]
        # find higher and lower ocean levels 
        z_idx_closest = np.abs(ocean_z - mean_shelf_topg[b_idx]).argmin()
        if (ocean_z[z_idx_closest] - mean_shelf_topg[b_idx]) < 0:
            z_idx_lower = z_idx_closest
            z_idx_higher = z_idx_closest - 1
            # special case for (mean_shelf_topg >= uppermost ocean level)
            if z_idx_closest == 0:
                z_idx_lower = 0
                z_idx_higher = 0
        else:
            z_idx_lower = z_idx_closest + 1
            z_idx_higher = z_idx_closest 
        
        
        # define sampling points
        z_l = ocean_z[z_idx_lower]      # z_lower
        z_h = ocean_z[z_idx_higher]     # z_higher
        z_i = mean_shelf_topg[b_idx]    # z_interpolate
        dz = z_h - z_l                  # delta(z_l,z_h)
        
        # iterate fields
        for f in fields.keys():
            bf_l = fields[f]['field_in'][z_idx_lower, basins==b_val]     # basin field lower
            bf_h = fields[f]['field_in'][z_idx_higher, basins==b_val]    # basin field higher
            bf_i = np.zeros_like(bf_l)      # basin field interpolate
            bf_i[:] = np.nan
            
            # use higher values if lower values not present
            if np.ma.is_masked(bf_l):
                bf_i[:] = bf_h[:]
            # use uppermost ocean level when mean_shelf_topg >= ocean_z[0]
            elif ( (z_idx_lower == 0) and (z_idx_higher == 0) ):
                bf_i[:] = bf_h[:]
            # regular interpolation case
            else:
                # compute linear interpolation for each basin grid point between 
                #   higher and lower depth
                for bf_idx, bf_val in enumerate(bf_l):
                    bf_i[bf_idx] =      (z_h - z_i) / dz * bf_l[bf_idx] \
                                      + (z_i - z_l) / dz * bf_h[bf_idx]
            # store interpolated basin field
            fields[f]['field_out'][basins==b_val] = bf_i
    

    t_interp_end = t.time()    
    

    ### --------------------------- basin averaging -------------------------
    t_basin_ave_start = t.time()
    if args.verbose:
        print('... calculate basin means')

    # iterate through fields to be filled
    for f in fields.keys():
        if args.verbose: 
            print('\t > ', f)
        
        # check field dimension
        field = fields[f]['field_out']
        field_shape = field.shape

        for b in basin_vals:
            # mask of current basin
            m__basin = (basins==b)
            # mask of contshelf condition == True
            m__contshelf = (contshelf_mask==2)
            # mask of current basin & contshelf
            m__basin_contshelf = m__basin & m__contshelf

            # calculate mean and write to array 
            basin_mean = field[m__basin].mean()
            #fields[f]['field_out_basin_mean'][m__basin] = basin_mean
            fields[f]['field_out_basin_mean'][m__basin_contshelf] = basin_mean

    t_basin_ave_end = t.time()
    

    ### -------------------- write result to output file -------------------- 

    ###   read input file 
    #   -> read file with regridded variables to be processed (again)
    #       to create identical output file (incl. correct dimensions
    #       and attributes) but with modified fields

    if args.verbose:
        print("... writing output file '" + args.out_file + "'")

    t_write_outfile_start = t.time()
    
    

    try:
        nc_src = CDF(args.field_file, 'r')
    except:
        print("FIELD_FILE '" + args.field_file + "' can't be found! Exiting.")
        sys.exit(1)

    ### create file for output 
    nc_dst = CDF(args.out_file, "w", format='NETCDF4')

    # copy general dimensions from field input file
    for name, dimension in nc_src.dimensions.items():
        nc_dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # create dimension n_basins
    nc_dst.createDimension(n_basin_name, n_basin_size)

    ### copy and modify global attributes
    glob_dict = nc_src.__dict__
    glob_dict['filename'] = os.path.basename(args.out_file)
    glob_dict['title'] = ("MOM output variables interpolated to PISM grid "
                          "with missing values filled and interpolated to "
                          "correct basin input depth for PISM")
    # remove old attributes from MOM grid
    try:
        del glob_dict['grid_type']
        del glob_dict['grid_tile']
    except:
        pass

    # modify history string
    cmd_line = ' '.join(sys.argv)
    histstr = t.asctime() + ': ' + cmd_line + "\n "

    if 'history' in glob_dict.keys():
        glob_dict['history'] = histstr + glob_dict['history']
    elif 'History' in glob_dict.keys():
        glob_dict['History'] = histstr + glob_dict['History']
    else:
        glob_dict['history'] = histstr
    
    # set global attributes
    nc_dst.setncatts(glob_dict)

    ### copy variables from field_file incl attributes & dimensions
    var_copy = ['lon','lon_bnds','lat','lat_bnds','x','y']
    for name, var in nc_src.variables.items():
        if name in var_copy:
            # create variables with correct datatype and dimensions
            nc_dst.createVariable(name, var.datatype, var.dimensions)
            # copy variable attributes
            nc_dst[name].setncatts(nc_src[name].__dict__)
            # copy variable data
            nc_dst[name][:] = nc_src[name][:]
    for name, var in nc_src.variables.items():
        if name in ['time','time_bnds']:
            # create variables with correct datatype and dimensions
            nc_dst.createVariable(name, var.datatype, var.dimensions)
            # copy variable attributes
            nc_dst[name].setncatts(nc_src[name].__dict__)
            # copy variable data
            if name == 'time':
                nc_dst['time'][:] = nc_src['time'][-1]
            if name == 'time_bnds':
                nc_dst['time_bnds'][:] = nc_src['time_bnds'][-2]

    
    # write depth condensed fields     
    for name in fields:
        if name == 'temp':
            name_out = 'theta_ocean'
            field_unit = 'Celsius'
        elif name == 'salt':
            name_out = 'salinity_ocean'
            field_unit = 'g/kg'
        else:
            name_out = name + '_dcon'   # depth condensed
            field_unit = 'unknown'
        
        # remove depth dimension but keep time
        dim_out = list(nc_src.variables['temp'].dimensions)
        dim_out.remove("st_ocean")

        for out_field in ['field_out', 'field_out_basin_mean']:
            if out_field == 'field_out_basin_mean':
                name_out += '_basin_mean'
            nc_dst.createVariable(name_out, 
                                  nc_src.variables[name].datatype, 
                                  #nc_src.variables['temp'].dimensions[-2:])
                                  dim_out)
            var_dict = nc_src[name].__dict__
            var_dict['long_name'] += (" depth_condensed: linear depth interpolated "
                                      "for basin shelf depth")
            if out_field == 'field_out_basin_mean':
                var_dict['long_name'] += (" -- mean value for continental shelf"
                                          " region of each basin")
            var_dict['units'] = field_unit
            nc_dst[name_out].setncatts(var_dict)
            field_out = fields[name][out_field]
            nc_dst[name_out][:] = np.array(field_out)[np.newaxis,:] # added time dimension

    # write basins
    nc_dst.createVariable('basins', 
                          basins_dtype, 
                          basins_dim[-2:])
    nc_dst['basins'].setncatts(basins_dict)
    nc_dst['basins'][:] = basins[:]
    
    # write basin list
    nc_dst.createVariable('basin_list', 
                          shelf_depth_basin_dtype, 
                          shelf_depth_basin_dim)
    nc_dst['basin_list'].setncatts(shelf_depth_basin_dict)
    nc_dst['basin_list'][:] = shelf_depth_basin_list[:] 
    
    # write basin mean_shelf_topg depth
    nc_dst.createVariable('mean_shelf_topg', 
                          mean_shelf_topg_dtype, 
                          mean_shelf_topg_dim)
    nc_dst['mean_shelf_topg'].setncatts(mean_shelf_topg_dict)
    nc_dst['mean_shelf_topg'][:] = np.array(mean_shelf_topg)[np.newaxis,:] # added time dimension


    nc_dst.close()
    nc_src.close()

    t_write_outfile_end = t.time()
    t_main_end = t.time()





    # -------------------- performance -------------------- 

    if args.verbose | args.time:
        t_main          = t_main_end            - t_main_start
#        t_read_bmask    = t_read_bmask_end      - t_read_bmask_start
#        t_read_edge     = t_read_edge_end       - t_read_edge_start
#        t_read_infile   = t_read_infile_end     - t_read_infile_start
        t_read_files    = t_read_end            - t_read_start
        t_fill          = t_fill_end            - t_fill_start
        t_interp        = t_interp_end          - t_interp_start        
        t_basin_ave     = t_basin_ave_end       - t_basin_ave_start
        t_write_outfile = t_write_outfile_end   - t_write_outfile_start

        format_total = "{:<15} \t\t {:9.2f} s \t {:6.2f} %"
        format_sub   = "\t{:<15} \t {:9.2f} s \t {:6.2f} %"

        print()
        print('{:-^58}'.format(' elapsed time '))
        print(format_total.format('total', t_main, t_main/t_main*100))
        print('{:.^58}'.format(''))
#        print(format_sub.format('read basin mask', t_read_bmask, 
#                                    t_read_bmask/t_main*100))
#        print(format_sub.format('read edge mask', t_read_edge, 
#                                    t_read_edge/t_main*100))
#        print(format_sub.format('read field file', t_read_infile, 
#                                    t_read_infile/t_main*100))
        print(format_sub.format('read files', t_read_files, 
                                    t_read_files/t_main*100))
        print(format_sub.format('fill empty cells', t_fill, 
                                    t_fill/t_main*100))
        print(format_sub.format('interpolate depths', t_interp, 
                                    t_interp/t_main*100))
        print(format_sub.format('basin average', t_basin_ave,
                                    t_basin_ave/t_main*100))
        print(format_sub.format('write output file', t_write_outfile, 
                                    t_write_outfile/t_main*100))
        print('{:.^58}'.format(''))
        print()

