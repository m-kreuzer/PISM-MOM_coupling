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
import time
import copy as cp
import argparse
import numpy as np
import xarray as xr
import cftime
## for debugging
#import code 
#import matplotlib.pyplot as plt
#import warnings
##warnings.filterwarnings('error')
##warnings.simplefilter('always')


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
    t_main_start = time.time()
     
    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()
       

    # -------------------- read basin mask --------------------  
    if args.verbose:
        print("... reading basin mask from file '" + args.basin_file + "'")
    t_read_start = time.time()
    t_read_bmask_start = time.time()

    
    try:
        ds_basins = xr.open_dataset(args.basin_file)
    except:
        print("BASIN_FILE '", args.basin_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    basins = ds_basins['basins']
    # mask values:
    #   0 = grounded ice or open ocean
    #   1 = continental shelf but ignored
    #   2 = continental shelf as used in PICO
    contshelf_mask = ds_basins['pico_contshelf_mask']

    # check dimensions
    basin_ndim = len(basins.dims)
    if basin_ndim != 2 and basin_ndim != 3:
        raise ValueError( str("basin field has dimensions " + \
                            str( basins.dims ) \
                            + ". Expected: 2 or 3 dimensions.") )
    if basin_ndim == 3:
        # cut of time dimension and take first time slice
        basins = basins.isel(time=0)

    contshelf_ndim = len(contshelf_mask.dims)
    if contshelf_ndim != 2 and contshelf_ndim != 3:
        raise ValueError( str("pico continental shelf mask has dimension " + \
                            str( contshelf_mask.dims ) + \
                            ". Expected: 2 or 3 dimensions.") )
    if contshelf_ndim == 3:
        # cut of time dimension and take first time slice
        contshelf_mask = contshelf_mask.isel(time=0)

    assert_str = ("non matching dimensions of variables 'basins' and "
                  "'pico_contshelf_mask' in BASIN_FILE '{}'.")
    assert basins.shape == contshelf_mask.shape, \
                assert_str.format(args.basin_file)
        
    # create list with all occuring basins
    basin_vals = np.unique(basins)
    # remove basin 0
    basin_list = np.delete(basin_vals, np.where(basin_vals==0) ) 
    # make sure datatype is integer
    basin_list = basin_list.astype(int)

    t_read_bmask_end = time.time()

    # -------------------- read edge mask -------------------- 
    t_read_edge_start = time.time()

    if args.verbose:
        print("... reading edge mask from file '" + args.edge_file + "'")

    try:
        ds_field_edge = xr.open_dataset(args.edge_file)
        #nc_e = CDF(args.edge_file, 'r')
    except:
        print("EDGE_FILE '", args.edge_file, "' can't be found! Exiting.")
        sys.exit(1)
        
    # read edge array    
    #  -> important to read in as bool, otherwise type is int >> bad performance
    try:
        field_edge = ds_field_edge['field_edge'].astype(bool)
    except:
        print("Variable 'field_edge' can't be read from EDGE_FILE '" 
                + args.edge_file + "'!")
    
    # check basin dimension
    field_edge_ndim = len(field_edge.dims)
    if field_edge_ndim != 3:
        raise ValueError( str("edge field has dimensions " + \
                            str( field_edge.dims ) + \
                            ". Expected: 3 dimensions.") )
    
    t_read_edge_end = time.time()


    # -------------------- read variables -------------------- 
    #   -> read file with regridded variables to be processed
    t_read_infile_start = time.time()
    
    if args.verbose:
        print("... reading VAR_FILL fields '" + str(args.var_fill) + 
                "' from FIELD_FILE '" + args.field_file + "'")

    try:
        ds_field = xr.open_dataset(args.field_file)
    except:
        print("FIELD_FILE '" + args.field_file + "' can't be found! Exiting.")
        sys.exit(1)
    
    # create condensed data structure (no depth) for output
    ds_field_out = ds_field.isel(st_ocean=0, drop=True)
    for f in args.var_fill:
        ds_field_out[f] = ds_field_out[f]*np.nan
        # copy attributes
        for k in ds_field[f].attrs.keys():
            if 'name' in k:
                ds_field_out[f].attrs[k] = ds_field[f].attrs[k]
        ds_field_out[f].attrs['units'] = ds_field[f].attrs['units']

    # extract vertical axis
    ocean_z = ds_field['st_ocean']              # units: dbars (interpreting as m)
    ocean_z = -1 * ocean_z                      # positive = upwards
    
    # read field arrays 
    for f in args.var_fill:
        field_ndim = len(ds_field[f].shape)

        if ( field_ndim != 4) :
            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
                    args.field_file + "' has dimensions " + \
                    str(ds_field[f].dims) + ". Expected: (time,z,y,x)."
            raise ValueError( str(err_str) )

        with xr.set_options(keep_attrs=True):
            ds_field_out = ds_field_out.assign({f'{f}_basin_mean': ds_field_out[f]*np.nan}) 


    t_read_infile_end = time.time()

    # -------------------- read basin depths -------------------- 
    #   -> read file with regridded variables to be processed
    
    # input file for shelf depth is given -> read it
    if args.basin_shelf_depth_file is not None:
        
        t_read_depthfile_start = time.time()
        if args.verbose:
            print("... reading basin shelf topography from BASIN_SHELF_DEPTH_FILE '" 
                  + args.basin_shelf_depth_file + "'")
    
        try:
            ds_basin_shelf_depth = xr.open_dataset(args.basin_shelf_depth_file)
        except:
            print("BASIN_SHELF_DEPTH_FILE '" + args.basin_shelf_depth_file + 
                  "' can't be found! Exiting.")
            sys.exit(1)
        
        n_basin_size = ds_basin_shelf_depth.dims['n_basin']
        shelf_depth_basin_list  = ds_basin_shelf_depth['basin'].data
        mean_shelf_topg         = ds_basin_shelf_depth['mean_shelf_topg'].isel(time=0).data

        t_read_depthfile_end = time.time()
        
    else:        
        ## set default values

        #if args.verbose:
        #    print("... no BASIN_SHELF_DEPTH_FILE given! \n"
        #          "\t -> using default depth -500m for all basins")
    
        #n_basin_name = 'n_basin'
        #n_basin_size = len(basin_vals)
        #
        #shelf_depth_basin_list  = basin_list
        ##shelf_depth_basin_list  = np.array(np.arange(1,20))  #temorary fix
        #shelf_depth_basin_dtype = int
        #shelf_depth_basin_dim   = ('n_basin',)
        #shelf_depth_basin_dict  = \
        #    dict({'long_name':  "list of valid PISM/PICO basins"})
        #
        #mean_shelf_topg         = np.ones_like(shelf_depth_basin_list) * -500
        #mean_shelf_topg_dtype   = float
        #mean_shelf_topg_dim     = ('n_basin',)
        #mean_shelf_topg_dict    = \
        #    dict({'long_name':  "mean basin topography of ice shelf areas", \
        #          'units':      "m", \
        #          'axis':       "Z", \
        #          'positive':   "up"})        

        raise ValueError("No BASIN_SHELF_DEPTH_FILE given (-b). Exiting! ")
    
    t_read_end = time.time()


    ###  -------------------- fill missing field values -------------------- 
    #  -> empty cells are filled with the mean of all values on the edge 
    #       from missing to non-missing values. This is done for each basin
    #       in each layer in each variable/field.

    t_fill_start = time.time()

    if args.verbose:
        print('... fill empty grid cells')

    ### --- working but very inefficient implementation with xr.where() --- ###

    #basins.load()
    #field_edge.load()
    ## iterate through time
    #for t, t_val in enumerate(ds_field.time):
    #    ds_field.isel(time=t).load()

    #    # iterate through fields to be filled
    #    for f in args.var_fill:
    #        if args.verbose: 
    #            print('\t > ', f)
    #        
    #        # check field dimension
    #        field = ds_field[f].sel(time=t_val)
    #        field_shape = field.shape
    #        depth_levels = field_shape[0]

    #        if (field_shape[-2:] != basins.shape[-2:]) :
    #            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
    #                        args.field_file + "' has shape " + \
    #                        str(field_shape[-2:]) + " which does not match the " + \
    #                        "shape of 'basins' " + str(basins.shape) + \
    #                        " of BASIN_FILE '" + args.basin_file + "'!" 
    #            raise ValueError( str(err_str) ) 
    #        
    #        if (field_shape[-2:] != field_edge.shape[-2:]) :
    #            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
    #                        args.field_file + "' has shape " + \
    #                        str(field_shape[-2:]) + " which does not match the " + \
    #                        "shape of 'field_edge' " + str(field_edge.shape[-2:]) + \
    #                        " of EDGE_FILE '" + args.edge_file + "'!" 
    #            raise ValueError( str(err_str) ) 
    #        
    #        # calculate basin mean of edge cells and apply for empty cells
    #        for z,z_val in enumerate(field.st_ocean):
    #            f_z =  ds_field[f].isel(time=t, st_ocean=z)
    #            fe_z = field_edge.sel(st_ocean=z_val)

    #            for b in basin_vals:
    #                ### create masks
    #                ## mask of current basin
    #                m__basin = (basins==b)
    #                ## mask of field edge AND current basin
    #                m__fedge_A_basin = fe_z & m__basin
    #                ## mask of missing fields AND current basin
    #                m__fmiss_A_basin = np.isnan(f_z) & m__basin

    #                ### calculate mean and write to array 
    #                mean_field_basin_edge = f_z.where(m__fedge_A_basin==True).mean()
    #                with ds_field[f] as fs:
    #                    mask_assign = ( (fs.coords['time']==t_val) & 
    #                                    (fs.coords['st_ocean'] == z_val) & 
    #                                    (m__fmiss_A_basin==True)        )
    #                ds_field[f] = xr.where(mask_assign, mean_field_basin_edge, ds_field[f])

    # iterate through fields to be filled
    for f in args.var_fill:
        if args.verbose: 
            print(f'\t field: {f}')

        # extract data from xarray Dataset, do the filling and put it back in
        da_data = ds_field[f].data      # dim: (time, z, y, x)

        # iterate through time
        for t, t_val in enumerate(ds_field.time):
            
            # check field dimension
            field = da_data[t,:]
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
            for z,z_val in enumerate(ds_field[f].st_ocean):
                for b in basin_vals:
                    ### create masks
                    # mask of current basin
                    m__basin = (basins.data==b)
                    # mask of field edge AND current basin
                    m__fedge_A_basin = field_edge.data[z,:] & m__basin
                    # mask of missing fields AND current basin
                    m__fmiss_A_basin = np.isnan(da_data[t,z,:]) & m__basin

                    ### calculate mean and write to array 
                    mean_field_basin_edge = np.nanmean(da_data[t,z,m__fedge_A_basin])
                    da_data[t,z,m__fmiss_A_basin] = mean_field_basin_edge

        # apparently modifications of da_data are reflected in ds_field[f] and the line below is not needed
        #ds_field[f'{f}_fill'] = xr.DataArray(da_data, coords=ds_field[f].coords)

            
    t_fill_end = time.time()

    ### --------------- depth interpolation ---------------------------------
    t_interp_start = time.time()    
    #FIXME: don't know how to use xarray's internal interpolation method 
    #       to assign edge nan's with last valid value. Sticking to old 
    #       implementation for now.

    if args.verbose:
        print('... interpolate basin shelf depths')
    
    # check whether basin lists are matching
    assert_str = ("non matching basins between BASIN_FILE '{}' and "
                  "BASIN_SHELF_DEPTH_FILE '{}'.")
    assert set(basin_list) == set(shelf_depth_basin_list), \
                assert_str.format(args.basin_file, args.basin_shelf_depth_file)


    # iterate basins 
    for b_idx, b_val in enumerate(shelf_depth_basin_list):
        if args.verbose:
            print('\t > basin ', b_val, ' / ', shelf_depth_basin_list.max())
            
        # depth of current basin: mean_shelf_topg[idx]
        # find higher and lower ocean levels 
        z_idx_closest = np.abs(ocean_z.data - mean_shelf_topg[b_idx]).argmin()
        if (ocean_z.data[z_idx_closest] - mean_shelf_topg[b_idx]) < 0:
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
        z_l = ocean_z.data[z_idx_lower]      # z_lower
        z_h = ocean_z.data[z_idx_higher]     # z_higher
        z_i = mean_shelf_topg[b_idx]    # z_interpolate
        dz = z_h - z_l                  # delta(z_l,z_h)
        
        # iterate fields
        for f in args.var_fill:
            # extract data from xarray Dataset, do the filling and put it back in
            da_data = ds_field_out[f].data      # dim: (time, y, x)

            # iterate time
            for t,t_val in enumerate(ds_field[f].time):
                bf_l = ds_field[f].data[t,z_idx_lower, basins==b_val]     # basin field lower
                bf_h = ds_field[f].data[t,z_idx_higher, basins==b_val]    # basin field higher
                bf_i = np.zeros_like(bf_l)      # basin field interpolate
                bf_i[:] = np.nan
                
                # use higher values if lower values not present
                if any(np.isnan(bf_l)):
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
                da_data[t,basins.data==b_val] = bf_i
                #fields[f]['field_out'][basins==b_val] = bf_i
    

    t_interp_end = time.time()    
    

    ### --------------------------- basin averaging -------------------------
    t_basin_ave_start = time.time()
    if args.verbose:
        print('... calculate basin means')

    # iterate through fields to be filled
    for f in args.var_fill:
        if args.verbose: 
            print('\t > ', f)
        
        for t, t_val in enumerate(ds_field_out.time):
            da_data = ds_field_out[f'{f}_basin_mean'].data

            for b in basin_vals:
                # mask of current basin
                m__basin = (basins.data==b)
                # mask of contshelf condition == True
                m__contshelf = (contshelf_mask.data==2)
                # mask of current basin & contshelf
                m__basin_contshelf = m__basin & m__contshelf

                # calculate mean and write to array 
                basin_mean = ds_field_out[f].data[t,m__basin_contshelf].mean()
                da_data[t,m__basin_contshelf] = basin_mean
                #fields[f]['field_out_basin_mean'][m__basin_contshelf] = basin_mean

    t_basin_ave_end = time.time()
    

    ### -------------------- write result to output file -------------------- 

    ###   read input file 
    #   -> read file with regridded variables to be processed (again)
    #       to create identical output file (incl. correct dimensions
    #       and attributes) but with modified fields

    if args.verbose:
        print("... writing output file '" + args.out_file + "'")

    t_write_outfile_start = time.time()
    
    
    #### global attributes
    ds_field_out.attrs['filename'] = os.path.basename(args.out_file)
    ds_field_out.attrs['title'] = ("MOM output variables interpolated to PISM grid "
                          "with missing values filled and interpolated to "
                          "correct basin input depth for PISM")
    # remove old attributes from MOM grid
    try:
        del ds_field_out.attrs['grid_type']
        del ds_field_out.attrs['grid_tile']
    except:
        pass

    # modify history string
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "

    if 'history' in ds_field_out.attrs:
        ds_field_out.attrs['history'] = histstr + ds_field_out.attrs['history']
    elif 'History' in ds_field_out.attrs:
        ds_field_out.attrs['History'] = histstr + ds_field_out.attrs['History']
    else:
        ds_field_out.attrs['history'] = histstr
    

    
    ### variables    
    # add some metadata 
    for f in args.var_fill:
        for v in ds_field_out.data_vars:
            if f in v:
                ds_field_out[v].attrs['long_name'] += \
                    " depth_condensed: linear depth interpolated for basin shelf depth"
                if '_basin_mean' in v:
                    ds_field_out[v].attrs['long_name'] += \
                        "; mean value for continental shelf region of each basin"

    # change temperature unit as required by PISM
    for v in ds_field_out.data_vars:
        if 'temp' in v:
            if ds_field_out[v].attrs['units'] == 'degrees C':
                ds_field_out[v].attrs['units'] = 'Celsius' 
        if 'salt' in v:
            if ds_field_out[v].attrs['units'] == 'psu':
                ds_field_out[v].attrs['units'] = 'g/kg' 

    # rename variables
    variable_rename = { 'temp':'theta_ocean', 
                        'salt':'salinity_ocean',
                        'temp_basin_mean':'theta_ocean_basin_mean',
                        'salt_basin_mean':'salinity_ocean_basin_mean'} 
    ds_field_out = ds_field_out.rename_vars(variable_rename)


    # add basins and mean_shelf_depth for output
    ds_field_out = ds_field_out.assign({'basins':basins})
    ds_field_out = xr.merge([ds_field_out, ds_basin_shelf_depth.drop_dims('time')])
    
    ### time bounds (required by PISM-PICO for time series input)
    if 'time_bnds' in ds_field_out.variables:
        # rename time bounds
        ds_field_out = ds_field_out.rename_vars({'time_bnds':'time_bounds'}) 
        ds_field_out['time'].attrs['bounds'] = 'time_bounds'
    elif ('time_bounds' not in ds_field_out.variables) and \
        (ds_field_out.time.size > 1) :
        # create bounds when data is timeseries (not required for single timestamp)
        ds_field_out = ds_field_out.assign(time_bounds=cp.deepcopy( 
                                              ds_field_out.time.expand_dims({'nv':2}, axis=1)))
        ds_field_out['time'].attrs['bounds'] = 'time_bounds'
    
        for t,t_val in enumerate(ds_field_out.get_index('time')):
            ds_field_out['time_bounds'].loc[dict(time=t_val, nv=0)] = cftime.DatetimeNoLeap(t_val.year, 1, 1)
            ds_field_out['time_bounds'].loc[dict(time=t_val, nv=1)] = cftime.DatetimeNoLeap(t_val.year+1, 1, 1)

    # time units as required by PISM-PICO
    encoding_dict = dict({'time': {'units':'seconds since 0001-01-01 00:00:00',
                                   'dtype':'double'}})
    if 'time_bounds' in ds_field_out.variables:
        encoding_dict['time_bounds'] = {'units':'seconds since 0001-01-01 00:00:00',
                                        'dtype':'double'}

    ### write output
    ds_field_out.to_netcdf(args.out_file, encoding=encoding_dict)

    t_write_outfile_end = time.time()
    t_main_end = time.time()




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
