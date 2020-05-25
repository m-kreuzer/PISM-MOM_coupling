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

""" Identifies cells that have neighbours with missing values. 

usage: ./find_edge.py -i infile.nc -f field_name -e var_exclude -o outfile.nc \
        [-t] [-v]

Edge grid cells of a netcdf 3d variable are identified by this script and 
written to a netCDF output file. A cell is referred as 'edge' if at least one 
of it's direct neighbours (4 point stencil) has no value.

Arguments:
    -i infile 
        a netCDF file with variable 'field' to be processed
    -f field_name
        a variable in 'infile' which is examined
    -e var_exclude_1 [var_exclude_2] [...]
        variable name(s) from 'infile' not to be copied to 'outfile'
    -o outfile
        name of the output netCDF file
    -t (optional)
        print script time statistics
    -v (optional)
        print verbose output
        
The output file will be used by the script regriddedMOM-to-PISM_processing.py

This script was created as a preprocessing tool for regridding ocean output 
fields to landice grid. It is designed for the context of coupling the ice 
sheet model PISM/PICO to the climate model POEM at PIK.

"""

import sys
import time as t
import argparse
import numpy as np
import collections as col
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


def mark_edge(data, mask):
    """Marks edge cells of 2d masked array 'data' in 2d field 'mask'"""

    # iterate colums
    for j in range(data.shape[0]):
        # iterate rows
        for i in range(data.shape[1]):

            # check if cell is masked
            if data.mask[j,i] == False:
                mask[j,i] = check_neighbor(data,j,i)

    return mask
    

        
def check_neighbor(data, col, row):
    """Returns true if one of the data point's neighbours is masked.
    
    Left, right, front and behind neighbours (if existing) are checked for 
    given column and row indices on 2d field 'data'. 
    """

    col_p1 = col+1
    col_m1 = col-1
    row_p1 = row+1
    row_m1 = row-1

    ### correction for domain edges (prevents out of bounds)
    # rightmost column
    if col >= data.shape[0]-1:
        col_p1 = col
    # leftmost column
    if col <= 0:
        col_m1 = col
    # uppermost row
    if row <= 0:
        row_m1 = row    
    # lowermost row
    if row >= data.shape[1]-1:
        row_p1 = row
        #print("row_p1", row_p1)

    return data.mask[col_p1,row] + data.mask[col_m1,row] + \
            data.mask[col,row_p1] + data.mask[col,row_m1]

    


if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(
                description=
                "Identifies cells that have neighbours with missing values",
                epilog=
                "Edge grid cells of 3d variable FIELD in netCDF INFILE are     \
                identified and written to netCDF file OUTFILE. A cell is       \
                referred as 'edge' if at least one of it's direct neighbours   \
                has no value. VAR_EXCLUDE variables are not copied from        \
                INFILE to OUTFILE."
            )

    parser.add_argument('-i', '--input', action="store", dest="infile",
                        required=True, 
                        help="input file with variable to be processed")
    parser.add_argument('-o', '--output', action="store", dest="outfile",
                        required=True, 
                        help="file to store calculated edge cells")
    parser.add_argument('-f', '--field', action="store", dest="field",
                        required=True, 
                        help="name of variable to be examined")
    parser.add_argument('-e', '--exclude', action="store", dest="var_exclude",
                        required=True, nargs='+',
                        help="list of variable names not to copy to output \
                                file")
    parser.add_argument('-t', '--time', action="store_true", 
                        help="print script timings")
    parser.add_argument('-v', '--verbose', action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()




    # -------------------- general setup --------------------  
    t_main_start = t.time()

    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()
        
    # a list of possible x,y-dimensions names
    xdims = ['x', 'x1']
    ydims = ['y', 'y1']
    
    # -------------------- read input file --------------------  
    #   -> read file with regridded variables to be processed
    if args.verbose:
        print("... reading input file '" + args.infile + "'")

    t_read_infile_start = t.time()    

    try:
        nc_src = CDF(args.infile, 'r')
    except:
        print("INFILE '" + args.infile + "' can't be found! Exiting.")
        sys.exit(1)
        
    # assign x,y dimension
    for dim in xdims:
        if dim in list(nc_src.dimensions.keys()):
            xdim = dim
    for dim in ydims:
        if dim in list(nc_src.dimensions.keys()):
            ydim = dim
    
    # read field array 
    try:
        field = np.squeeze(nc_src.variables[args.field][:])
    except:
        print("Variable FIELD '" + args.field + "' can't be read from file '" +
                args.infile + "'!")

    field_dim = nc_src.variables[args.field].dimensions
    
    # check field dimension
    field_ndim = len(field.shape)
    if (field_ndim != 3) :
        err_str = "regridded field is of dimension " + str(field_ndim) + \
                    ". Expected: 3."
        raise ValueError( str(err_str) )

    t_read_infile_end = t.time()

    # -------------------- identify edge cells --------------------  
    #   -> determine which cells of field are edge cells
   
    t_calc_edge_start = t.time()

    if args.verbose:
        print("... calculating edge cells of given field '" + args.field + "'")

    field_edge = np.zeros_like(field.data, dtype=bool)
    depth_levels = field.shape[0]

    # compute field edges for all vertical levels
    for z in range(depth_levels):
        if args.verbose:
            print("     level: ", z+1, "/", depth_levels)
        mark_edge( field[z,:], field_edge[z, :] )


    t_calc_edge_end = t.time()

    # -------------------- write output --------------------  
    t_write_outfile_start = t.time()

    if args.verbose:
        print("... writing output file '" + args.outfile + "'")

    # create file for output
    nc_dst = CDF(args.outfile, "w", format='NETCDF4')

    # copy general dimensions from infile
    for name, dimension in nc_src.dimensions.items():
        nc_dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # copy coordinate variables from infile incl attributes & dimensions
    for name, var in nc_src.variables.items():
        if name not in args.var_exclude:
            nc_dst.createVariable(name, var.datatype, var.dimensions)
            # copy variable attributes
            nc_dst[name].setncatts(nc_src[name].__dict__)
            nc_dst[name][:] = nc_src[name][:]

    # create and write field_edge variable
    if field_dim[0] == 'time' :
        nc_dst.createVariable('field_edge', 'b', field_dim[1:])
    else:
        nc_dst.createVariable('field_edge', 'b', field_dim)
    var_dict = col.OrderedDict([
         ('long_name', 'mask determing edge between missing/non-missing values'),
         ('valid_range', np.array([0, 1], dtype=np.int32)),
         ('fill_value', -1),
         ('coordinates', 'lat lon')])
    nc_dst['field_edge'].setncatts(var_dict)
    nc_dst.variables['field_edge'][:] = field_edge[:]

    ### write history attribute 
    # concat all commandline arguments to string (separated by ' ')   
    cmd_line = ' '.join(sys.argv)
    nc_dst.history = t.asctime() + ': ' + cmd_line + "\n "

    nc_dst.close()
    nc_src.close()


    t_write_outfile_end = t.time()
    t_main_end = t.time()





    # -------------------- performance -------------------- 

    if args.verbose | args.time:
        t_main          = t_main_end            - t_main_start
        t_read_infile   = t_read_infile_end     - t_read_infile_start
        t_calc_edge     = t_calc_edge_end       - t_calc_edge_start
        t_write_outfile = t_write_outfile_end   - t_write_outfile_start

        format_total = "{:<15} \t\t {:9.2f} s \t {:6.2f} %"
        format_sub   = "\t{:<15} \t {:9.2f} s \t {:6.2f} %"

        print()
        print('{:-^58}'.format(' elapsed time '))
        print(format_total.format('total', t_main, t_main/t_main*100))
        print('{:.^58}'.format(''))
        print(format_sub.format('read input file', t_read_infile, 
                                    t_read_infile/t_main*100))
        print(format_sub.format('calculate edges', t_calc_edge, 
                                    t_calc_edge/t_main*100))
        print(format_sub.format('write output file', t_write_outfile, 
                                    t_write_outfile/t_main*100))
        print('{:.^58}'.format(''))
        print()

