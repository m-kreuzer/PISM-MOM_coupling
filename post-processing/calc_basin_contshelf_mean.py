#!/usr/bin/env python3

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

""" Calculate basin mean values of input fields

usage: ./calc_basin_mean.py -i field_file -b basin_file -f var1 [var2 ...] \
                            [-c] -o out_file [-t] [-v]

Given fields are averaged for each PISM/PICO basin in continental shelf mask
region.

Arguments:
    -i, --input field_file
        a netCDF file with variables to be processed (PISM grid)
    -b, --basins basin_file
        a netCDF file with variable 'basins' describing the basins
        of the PICO model as well as the variable 'pico_contshelf_mask'
        with the continental shelf mask
    -f, --fill  var1 [var2 ...]
        list of variable names whose missing cells should be filled
    -o, --output out_file
       file to store processed fields which serve as PISM/PICO input 
    -t, --time (optional)
        print script time statistics
    -v, --verbose (optional)
        print verbose output        


This script was created as a post-processing tool for preparing output of ocean 
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
from tqdm import tqdm
#import matplotlib.pyplot as plt
#import warnings
##warnings.filterwarnings('error')
##warnings.simplefilter('always')


if __name__ == "__main__":
   

    parser = argparse.ArgumentParser(
                description=
                ("Calculate basin mean values of input fields"),
                epilog=
                ("Given fields are averaged for each PISM/PICO basin. If "
                 "continental shelf mask is given, only this region can be used "
                 "for averaging.")
            )

    parser.add_argument('-i', '--input', 
                        action="store", 
                        dest="field_file",
                        required=True, 
                        help="file with variables to be processed")
    parser.add_argument('-b', '--basins', 
                        action="store", 
                        dest="basin_file",
                        required=True, 
                        help="file with variable 'basins' from PICO model, "
                             "optionally also includes variable "
                             "'pico_contshelf_mask'")
    parser.add_argument('-f', '--fill', 
                        action="store", 
                        dest="var_fill",
                        required=True, 
                        nargs='+',
                        help="list of variable names whose missing cells     \
                                should be filled")
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
    
    # read field arrays 
    for f in args.var_fill:
        field_ndim = len(ds_field[f].shape)


        if field_ndim==2:
            data_has_time_dim=False
        elif field_ndim==3:
            data_has_time_dim=True
        else:
            err_str = "VAR_FILL variable '" + f + "' in FIELD_FILE '" + \
                    args.field_file + "' has dimensions " + \
                    str(ds_field[f].dims) + ". Expected: (time,y,x)."
            raise ValueError( str(err_str) )

    t_read_infile_end = time.time()
    t_read_end = time.time()

    ### ------------------- create data structure for output -------------------

    if data_has_time_dim:
        ds_basin_mean = xr.Dataset(
                            coords={"time": ds_field.time,
                                    "basin": basin_list})
                            
        if 'bounds' in ds_field.time.attrs:
            tb = ds_field.time.attrs['bounds']
            try:
                ds_basin_mean[tb] = ds_field[tb] 
            except KeyError:
                print(f'Warning: Variable "{tb}" cannot be found in {args.field_file}. Resume without copying {tb} to output file.')

    else:
        ds_basin_mean = xr.Dataset(coords={"basin": basin_list})

    ds_basin_mean['basin'].attrs['name'] = 'basin' 
    ds_basin_mean['basin'].attrs['reference'] = 'basin number according to ' +\
            f'basin mask definied in {args.basin_file}'

    for f in args.var_fill:
        if data_has_time_dim:
            ds_basin_mean[f] = (('time','basin'), 
                                np.zeros(
                                    (len(ds_basin_mean['time']),
                                    len(ds_basin_mean['basin']))
                                    )*np.nan
                                ) 
            ds_basin_mean[f'{f}_allcontshelf_mean'] = \
                    ('time', np.zeros(len(ds_basin_mean['time']))*np.nan)
        else:
            ds_basin_mean[f] = (('basin'), 
                                np.zeros(len(ds_basin_mean['basin']))*np.nan
                                )
            ds_basin_mean[f'{f}_allcontshelf_mean'] = np.nan

        # copy attributes
        for k in ds_field[f].attrs.keys():
            if ('name' in k) or ('units' in k):
                ds_basin_mean[f].attrs[k] = ds_field[f].attrs[k]
                ds_basin_mean[f'{f}_allcontshelf_mean'].attrs[k] = ds_field[f].attrs[k]

    ### --------------------------- basin averaging -------------------------
    t_basin_ave_start = time.time()
    if args.verbose:
        print('... calculate basin means')

    # iterate through fields to be filled
    for f in args.var_fill:
        if args.verbose: 
            print('\t > ', f)

        if data_has_time_dim:
            for t, t_val in enumerate(tqdm(ds_field.time, desc=f'{f}')):#, file=sys.stdout, position=0, leave=True)):
                da_data = ds_basin_mean[f].data

                # cont shelf mean for each basin
                for b, b_val in enumerate(basin_list):
                    # mask of current basin
                    m__basin = (basins.data==b_val)
                    # mask of contshelf condition == True
                    m__contshelf = (contshelf_mask.data==2)
                    # mask of current basin & contshelf
                    m__basin_contshelf = m__basin & m__contshelf

                    # calculate mean and write to array 
                    basin_select = ds_field[f].data[t,m__basin_contshelf]
                    if basin_select.size > 0:
                        da_data[t,b] = basin_select.mean()
                    else:
                        da_data[t,b] = np.nan

                # cont shelf mean for all basins
                basin_select = ds_field[f].data[t,m__contshelf]
                if basin_select.size > 0:
                    ds_basin_mean[f'{f}_allcontshelf_mean'].data[t] = basin_select.mean()
                else:
                    ds_basin_mean[f'{f}_allcontshelf_mean'].data[t] = np.nan


        else:   # data has no time dimension
            da_data = ds_basin_mean[f].data

            # cont shelf mean for each basin
            for b, b_val in enumerate(basin_list):
                # mask of current basin
                m__basin = (basins.data==b_val)
                # mask of contshelf condition == True
                m__contshelf = (contshelf_mask.data==2)
                # mask of current basin & contshelf
                m__basin_contshelf = m__basin & m__contshelf

                # calculate mean and write to array 
                basin_select = ds_field[f].data[m__basin_contshelf]
                if basin_select.size > 0:
                    da_data[b] = basin_select.mean()
                else:
                    da_data[b] = np.nan

            # cont shelf mean for all basins
            basin_select = ds_field[f].data[m__contshelf]
            if basin_select.size > 0:
                ds_basin_mean[f'{f}_allcontshelf_mean'].data = basin_select.mean()
            else:
                ds_basin_mean[f'{f}_allcontshelf_mean'].data = np.nan

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
    ds_basin_mean.attrs['filename'] = os.path.basename(args.out_file)
    ds_basin_mean.attrs['title'] = ("PISM-PICO continental shelf mask average of each basin")
    # remove old attributes from MOM grid
    try:
        del ds_basin_mean.attrs['grid_type']
        del ds_basin_mean.attrs['grid_tile']
    except:
        pass

    # modify history string
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "

    if 'history' in ds_basin_mean.attrs:
        ds_basin_mean.attrs['history'] = histstr + ds_basin_mean.attrs['history']
    elif 'History' in ds_basin_mean.attrs:
        ds_basin_mean.attrs['History'] = histstr + ds_basin_mean.attrs['History']
    else:
        ds_basin_mean.attrs['history'] = histstr
    

    
    ### variables    
    # add some metadata 
    for f in args.var_fill:
        for v in ds_basin_mean.data_vars:
            if f in v:
                if 'long_name' in ds_basin_mean[v].attrs:
                    if 'allcontshelf_mean' in v:
                        ds_basin_mean[v].attrs['long_name'] += \
                            "; mean value for all continental shelf regions"
                    else:
                        ds_basin_mean[v].attrs['long_name'] += \
                            "; mean value for continental shelf region of each basin"


    ## time units as required by PISM-PICO
    #encoding_dict = dict({'time': {'units':'seconds since 0001-01-01 00:00:00',
    #                               'dtype':'double'}})
    encoding_dict = dict({'basin': {'dtype':'i4'}})
    #if 'time_bounds' in ds_field_out.variables:
    #    encoding_dict['time_bounds'] = {'units':'seconds since 0001-01-01 00:00:00',
    #                                    'dtype':'double'}

    ### write output

    if data_has_time_dim:
        ds_basin_mean.to_netcdf(args.out_file, unlimited_dims='time', encoding=encoding_dict)
    else:
        ds_basin_mean.to_netcdf(args.out_file, encoding=encoding_dict)

    t_write_outfile_end = time.time()
    t_main_end = time.time()




    # -------------------- performance -------------------- 

    if args.verbose | args.time:
        t_main          = t_main_end            - t_main_start
        t_read_files    = t_read_end            - t_read_start
        t_basin_ave     = t_basin_ave_end       - t_basin_ave_start
        t_write_outfile = t_write_outfile_end   - t_write_outfile_start

        format_total = "{:<15} \t\t {:9.2f} s \t {:6.2f} %"
        format_sub   = "\t{:<15} \t {:9.2f} s \t {:6.2f} %"

        print()
        print('{:-^58}'.format(' elapsed time '))
        print(format_total.format('total', t_main, t_main/t_main*100))
        print('{:.^58}'.format(''))
        print(format_sub.format('read files', t_read_files, 
                                    t_read_files/t_main*100))
        print(format_sub.format('basin average', t_basin_ave,
                                    t_basin_ave/t_main*100))
        print(format_sub.format('write output file', t_write_outfile, 
                                    t_write_outfile/t_main*100))
        print('{:.^58}'.format(''))
