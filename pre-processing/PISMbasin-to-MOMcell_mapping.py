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

""" Creates mapping between PISM/PICO basins and MOM ocean cells.

usage: ./PISMbasin-to-MOMcell_mapping.py -b basin_file -m mom_file \
            -o out_file [-t] [-v]

PISM/PICO basins are read from file as well as the geometry of the MOM grid 
structure. The southern ocean cells around Antarctica are identified. The 
closest PISM grid cell is calculated to each ocean edge cell. A mapping of 
ocean edge cells with corresponding PICO basins is created along with it's 
fraction attribution to the specified basin. The computed fields are stored to 
an output file based on the MOM ocean grid.  

Arguments:
    -b basin_file
        a netCDF file with PISM/PICO variable 'basins'
    -m mom_file
        a MOM output file (netCDF) with variable 'temp' and 'area_t'
    -o out_file
        name of netCDF output file to store mapping of PISM basins to MOM cells
    -t (optional)
        print script time statistics
    -v (optional)
        print verbose output

The output file will be used by the script PISM-to-MOM_processing.py 

This script was created as a preprocessing tool for distributing the flux output
of the landice model PISM/PICO to the grid of ocean model MOM5. This was done
in the scope of coupling PISM to the climate model POEM at PIK.

"""


import sys
import os
import numpy as np
import copy
import collections as col
import time
import argparse
try:
    import netCDF4
    from netCDF4 import Dataset as CDF
except:
    raise ImportError("netCDF4 is not installed!")


def mark_edge(data, mask, j_limit):
    """Marks edge cells of 2d masked array 'data' in 2d field 'mask'
    
    'j_limit' is the maximum index to use for iterating 2d array 'data' in 
    first dimension.
    """
    # iterate colums
    #for j in range(data.shape[0]):
    for j in range(j_limit):
        # iterate rows
        for i in range(data.shape[1]):
            # check if cell exists
            if data.mask[j,i] == False:
                mask[j,i] = check_neighbor(data,j,i)
                    
    return mask
    

        
def check_neighbor(data, col, row):
    """Returns true if one of the data point's neighbours is masked.
    
    All eight neighbours (if existing) are checked for 
    given column and row indices on 2d field 'data'. 
    """
    # save neighbor indices
    col_p1 = col+1
    col_m1 = col-1
    row_p1 = row+1
    row_m1 = row-1
 
    ## correction for domain edges (prevents out of bounds)
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

    return data.mask[col_p1,row_p1]         \
            + data.mask[col_p1,row]         \
            + data.mask[col_p1,row_m1]      \
            + data.mask[col,row_m1]         \
            + data.mask[col_m1,row_m1]      \
            + data.mask[col_m1,row]         \
            + data.mask[col_m1,row_p1]      \
            + data.mask[col,row_p1] 
            
def get_closest_pism_cell(ocean_lat, ocean_lon, pism_lat_array, pism_lon_array):
    """ returns PISM cell indices closest to given MOM cell center coordinates
    
    'ocean_lat', 'ocean_lon'            - coordinates of one MOM cell
    'pism_lat_array', 'pism_lon_array'  - 2d latitude/longitude arrays of PISM grid
    """
    
    # define 2D potential field for latitude & longitude
    lat_pot = np.abs(pism_lat_array - ocean_lat)
    lon_pot = np.abs(pism_lon_array - ocean_lon)
    # minimize cumulative potential
    index_1d = (lat_pot + lon_pot).argmin()
    # transform to 2d index 
    return np.unravel_index(index_1d, pism_lat_array.shape)


if __name__ == "__main__":

    # -------------------- argument parser --------------------------
    parser = argparse.ArgumentParser(
                description=
                "Creates mapping between PISM/PICO basins and MOM ocean cells.",
                epilog= ("PISM/PICO basins are read from file as well as the "
                         "geometry of the MOM grid structure. The southern "
                         "ocean cells around Antarctica are identified. The "
                         "closest PISM grid cell is calculated to each ocean "
                         "edge cell. A mapping of ocean edge cells with "
                         "corresponding PICO basins is stored along with it's "
                         "fraction attributing to the specified basin. The "
                         "computed fields are stored to an output file based "
                         "on the MOM ocean grid.")
            )

    parser.add_argument('-b', '--basins', action="store", dest="basin_file",
                        required=True, 
                        help=("PISM output file with PICO variable 'basins'"))
    parser.add_argument('-m', '--mom', action="store", dest="MOM_file",
                        required=True, 
                        help=("MOM output file with coordinates 'xt_ocean', "
                              "'yt_ocean', 'st_ocean' and variables 'temp'"))
    parser.add_argument('-o', '--output', action="store", dest="out_file", 
                        required=True, 
                        help="file to store basin - ocean cell mapping")
    parser.add_argument('-l', '--limit', action='store', dest="southern_limit",
                        default=-60, type=float, 
                        help=("southern limit of ocean cells to consider for "
                              "basin mapping, units: degN"))
    parser.add_argument('-t', '--time', action="store_true", 
                        help="print script timings")
    parser.add_argument('-v', '--verbose', action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args() 
    
    
    # -------------------- general setup --------------------------
    t_main_start = time.time()
     
    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()
        
    # a list of possible x,y-dimensions names
    xdims = ['x', 'x1']
    ydims = ['y', 'y1']    
    
    ### ---------- read basin mask - BEGIN -----------------------------------
    t_read_files_start = time.time() 
    if args.verbose:
        print(" - reading basin mask from " + args.basin_file )
    try:
        nc_fh = CDF(args.basin_file, 'r')
    except:
        s = ("PISM basin file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.basin_file) )
        
    # assign x,y dimension
    for dim in xdims:
        if dim in list(nc_fh.dimensions.keys()):
            xdim = dim
    for dim in ydims:
        if dim in list(nc_fh.dimensions.keys()):
            ydim = dim
    
    # coordinate variable in x,y-direction
    basin_x = nc_fh.variables[xdim][:]
    basin_y = nc_fh.variables[ydim][:]
    basin_lat = nc_fh.variables['lat'][:]
    basin_lon = nc_fh.variables['lon'][:]
     # shift ocean longitudes to range [-180, 180] degE
    basin_lon_s = copy.deepcopy(basin_lon)
    basin_lon_s[basin_lon_s < -180] +=360
    basin_lon_s[basin_lon_s >  180] -=360

    # read basin array    
    pism_basins = np.squeeze(nc_fh.variables['basins'][:])#.astype(np.int32)
    # check basin dimension
    basin_ndim = len(pism_basins.shape)
    if basin_ndim != 2:
        s = ("Variable 'basins' from file '{}' is of dimension {}. Expected: 2")
        raise ValueError( s.format(args.basin_file, basin_ndim))     

    nc_fh.close()
    ### ---------- read basin mask - END -------------------------------------
    
        
    ### ---------- read ocean grid - BEGIN -----------------------------------
    if args.verbose:
        print(" - reading ocean grid from " + args.MOM_file )
    try:
        nc_fh = CDF(args.MOM_file, 'r')
    except:
        s = ("ocean file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.MOM_file) )
        
    # assign x,y dimension
    xdim = 'xt_ocean'
    ydim = 'yt_ocean'
    zdim = 'st_ocean'

    # deactivate mask for NC file input
    # (wrong valid range for longitude leads to missing values)    
    nc_fh.set_auto_mask(False)
    
    # coordinate variable in x,y-direction
    ocean_x =   nc_fh.variables[xdim][:]
    ocean_y =   nc_fh.variables[ydim][:]
    ocean_z =   nc_fh.variables[zdim][:]
    ocean_lat = nc_fh.variables['geolat_t'][:]
    ocean_lon = nc_fh.variables['geolon_t'][:]
    
    # activate mask for NC file input again
    nc_fh.set_auto_mask(True)
    
    # shift ocean longitudes to range [-180, 180] degE
    ocean_lon_s = copy.deepcopy(ocean_lon)
    ocean_lon_s[ocean_lon_s < -180] +=360
    ocean_lon_s[ocean_lon_s >  180] -=360

    # read tracer area
    oc_area = np.squeeze(nc_fh.variables['area_t'][:])
    # check area dimension
    ocean_a_ndim = len(oc_area.shape)
    if ocean_a_ndim != 2:
        s = ("Ocean variable 'area_t' from file '{}' is of dimension {}. "
             "Expected: 2")
        raise ValueError( s.format(args.MOM_file, ocean_a_ndim))

    # read temperature
    oc_temp = np.squeeze(nc_fh.variables['temp'][:])#.astype(np.int32)
    # check temp dimension
    ocean_t_ndim = len(oc_temp.shape)
    if ocean_t_ndim != 3:
        s = ("Ocean variable 'temp' from file '{}' is of dimension {}. "
             "Expected: 3")
        raise ValueError( s.format(args.MOM_file, ocean_t_ndim))
        
    oc_nlat = oc_temp.shape[1]
    oc_nlon = oc_temp.shape[2]
        
    nc_fh.close()

    t_read_files_end = time.time() 
    ### ---------- read ocean grid - END -------------------------------------
    
    t_process_start = time.time()
    if args.verbose:
        print(" - computing mapping of southern ocean cells and PISM grid cells")
    
    # create data structure for ocean cell information
    oc_south_edge = dict({'mask':None,              \
                          'pism_i':None,            \
                          'pism_j':None,            \
                          'pism_basin':None,        \
                          'pism_basin_ratio':None })
    
    # mask              - whether ocean cell belongs to southern ring around 
    #                       Antarctica or not
    #                       mask = True : cell DOES NOT belong to southern ring
    #                       mask = False: cell DOES belong to southern ring
    # pism_i, pism_j    - corresponding indices for the ocean cell center on
    #                       the PISM grid
    # pism_basin        - corresponding PISM/PICO basin for the ocean cell 
    #                       center with indices pism_i, pism_j on the PISM grid
    # pism_basin_ratio  - fraction of pism_basin flux to be routed in specific
    #                       ocean cell
    
    
    ### identify southern ocean edge cells around Antarctuca
    #southern_limit = -60 # unit: degN
    j_limit = np.abs(ocean_lat[:,0] - args.southern_limit).argmin()
    oc_south_edge['mask']= np.zeros((oc_nlat, oc_nlon), dtype=bool)
    #oc_south_edge['mask']= np.zeros_like(oc_temp[0,:].data, dtype=bool)
    mark_edge(oc_temp[0,:], oc_south_edge['mask'], j_limit)    
    # flip mask --> True: no southern edge cell; False: is southern edge cell
    oc_south_edge['mask'] = ~oc_south_edge['mask']


    # initialize rest of datastructure
    oc_south_edge['pism_i'] = np.ma.masked_array(oc_south_edge['mask']*-1,    \
                                                 mask= oc_south_edge['mask'], \
                                                 dtype= int)
    oc_south_edge['pism_j'] = np.ma.masked_array(oc_south_edge['mask']*-1,    \
                                                 mask= oc_south_edge['mask'], \
                                                 dtype= int)
    oc_south_edge['pism_basin'] = np.ma.masked_array(oc_south_edge['mask']*-1,\
                                                 mask= oc_south_edge['mask'], \
                                                 dtype= int)
    oc_south_edge['pism_basin_ratio'] = np.ma.masked_array(                   \
                                                 oc_south_edge['mask']*-1,    \
                                                 mask= oc_south_edge['mask'], \
                                                 dtype= float)
    
    # find corresponding PISM indices for ocean edge cell centers via coordinates
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_south_edge['mask'][j,i] == False:
                pism_index = get_closest_pism_cell(ocean_lat[j,i],   \
                                                   ocean_lon_s[j,i], \
                                                   basin_lat, \
                                                   basin_lon_s)                              
                oc_south_edge['pism_i'][j,i] = pism_index[0]
                oc_south_edge['pism_j'][j,i] = pism_index[1]

    
    ## create PISM field for verification
    #ocean_edge_on_pism_grid = np.empty_like(pism_basins.data)
    #ocean_edge_on_pism_grid[:] = np.nan
    #
    #fill_val = -40
    #for j in range(oc_nlat):
    #    for i in range(oc_nlon):
    #        if oc_south_edge['mask'][j,i] == False:
    #            ocean_edge_on_pism_grid[oc_south_edge['pism_i'][j,i], \
    #                                    oc_south_edge['pism_j'][j,i]] = fill_val

    
    # identify corresponding basin for each ocean edge cell 
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_south_edge['mask'][j,i] == False:
                 oc_south_edge['pism_basin'][j,i] = \
                                 pism_basins[oc_south_edge['pism_i'][j,i], \
                                             oc_south_edge['pism_j'][j,i]]


    ### check whether all PISM basins have at least one corresponding MOM cell
    mom_basin_list_tmp, mom_basin_count_tmp = \
                    np.unique(oc_south_edge['pism_basin'], return_counts=True)
    mom_basin_list    = mom_basin_list_tmp[~mom_basin_list_tmp.mask].data
    mom_basin_count   = mom_basin_count_tmp[~mom_basin_list_tmp.mask]
    mom_basin_area    = np.zeros_like(mom_basin_list, dtype=float)

    #pism_basin_list = np.unique(pism_basins).data
    pism_basin_list = np.unique(pism_basins)
    # remove basin 0
    pism_basin_list = np.delete(pism_basin_list, np.where(pism_basin_list==0) ) 
    
    # check whether all elements of basin_lists are matching
    assert set(mom_basin_list) == set(pism_basin_list), \
            'not all basins on PISM grid have a corresponding MOM cell!'
            
    # calculate cummulative area of ocean cells associated with each basin
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_south_edge['mask'][j,i] == False:
                basin = oc_south_edge['pism_basin'][j,i]
                mom_basin_area[np.where(mom_basin_list == basin)]  += oc_area[j,i]


    # write basin ratio for each ocean edge cell
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_south_edge['mask'][j,i] == False:
                basin = oc_south_edge['pism_basin'][j,i]
                basin_area = mom_basin_area[np.where(mom_basin_list == basin)]
                oc_south_edge['pism_basin_ratio'][j,i] =  oc_area[j,i] / basin_area
                # old version: ratio only by number of ocean cells per basin
                #count = mom_basin_count[np.where(mom_basin_list == basin)]
                #oc_south_edge['pism_basin_ratio'][j,i] =  1 / count[0]

    # check whether ratios of ocean cells in same basin add up to 1
    mom_basin_ratio_check = np.empty_like(mom_basin_area)*np.nan
    for b in mom_basin_list:
        basin_mask = (oc_south_edge['pism_basin'] == b)
        basin_ratio_sum = np.sum(oc_south_edge['pism_basin_ratio'][basin_mask])
        assert np.isclose(1, basin_ratio_sum), \
                'MOM cell area ratios in basin %s are not summing up to 1 but to %f' \
                % (b, basin_ratio_sum)
                      
    t_process_end = time.time()
    
    ### ---------- writing output  -----------------------------------
    t_write_file_start = time.time()
    if args.verbose:
        print(" - writing output to file " + args.out_file )    


        
    ### write oc_south_edge information to output file
    dim_copy = ['xt_ocean','yt_ocean']
    var_copy = ['geolat_t', 'geolon_t', 'xt_ocean', 'yt_ocean']
    
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "

    
    with CDF(args.MOM_file, 'r') as src, CDF(args.out_file, "w") as dst:
        # copy global attributes all at once via dictionary
        glob_dict = src.__dict__
        glob_dict['filename'] = os.path.basename(args.out_file)
        glob_dict['title'] = ("mapping of PISM/PICO basins to MOM ocean cells "
                              "at southern domain edge")
        
        if 'history' in glob_dict.keys():
            glob_dict['history'] = histstr + glob_dict['history']
        elif 'History' in glob_dict.keys():
            glob_dict['History'] = histstr + glob_dict['History']
        else:
            glob_dict['history'] = histstr
        
        dst.setncatts(glob_dict)
        
        # copy dimensions
        for name, dimension in src.dimensions.items():
            if name in dim_copy:
                dst.createDimension(name, len(dimension) )
                
        # copy variables
        for name, var in src.variables.items():
            if name in var_copy:
                x = dst.createVariable(name, var.datatype, var.dimensions)
                # fix wrong valid range attribute in geolon_t
                if name == 'geolon_t':
                    d = src[name].__dict__
                    d['valid_range'][0] = -360
                    dst[name].setncatts(d)
                    dst[name][:] = ocean_lon[:]
                else:
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                    dst[name][:] = src[name][:]
                
        ### write new variables    
        x = dst.createVariable('basin', 'i4', ('yt_ocean','xt_ocean'),
                fill_value=oc_south_edge['pism_basin'].fill_value)
        var_dict = col.OrderedDict([
             ('long_name', 'corresponding PISM/PICO basin of cell center'),
             ('valid_range', np.array([0, oc_south_edge['pism_basin'].max()], dtype=np.int32)),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['basin'].setncatts(var_dict)
        dst['basin'][:] = oc_south_edge['pism_basin'][:]
         
        x = dst.createVariable('basin_ratio', 'f8', ('yt_ocean','xt_ocean'),
                fill_value=oc_south_edge['pism_basin_ratio'].fill_value)
        var_dict = col.OrderedDict([
             ('long_name', ('ratio of corresponding PISM/PICO basin total '
                            'flux value to be mapped to cell')),
             ('valid_range', np.array([0, 1], dtype=np.int32)),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['basin_ratio'].setncatts(var_dict)
        dst['basin_ratio'][:] = oc_south_edge['pism_basin_ratio'][:]
        
        x = dst.createVariable('pism_i', 'i', ('yt_ocean','xt_ocean'),
                fill_value=oc_south_edge['pism_i'].fill_value)
        var_dict = col.OrderedDict([
             ('long_name', 'index of closest PISM grid cell to ocean cell center'),
             ('valid_range', np.array([0, basin_x.size], dtype=np.int32)),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['pism_i'].setncatts(var_dict)
        dst['pism_i'][:] = oc_south_edge['pism_i'][:]
        
        x = dst.createVariable('pism_j', 'i', ('yt_ocean','xt_ocean'),
                fill_value=oc_south_edge['pism_j'].fill_value)
        var_dict = col.OrderedDict([
             ('long_name', 'index of closest PISM grid cell to ocean cell center'),
             ('valid_range', np.array([0, basin_y.size], dtype=np.int32)),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['pism_j'].setncatts(var_dict)
        dst['pism_j'][:] = oc_south_edge['pism_j'][:]
        
    
    t_write_file_end = time.time()
    t_main_end = time.time()    
    
    
    # -------------------- performance -------------------- 
        
    if args.verbose | args.time:
        t_main          = t_main_end            - t_main_start
        t_read_files    = t_read_files_end      - t_read_files_start
        t_process       = t_process_end         - t_process_start
        t_write_file    = t_write_file_end      - t_write_file_start


        format_total = "{:<15} \t\t {:9.2f} s \t {:6.2f} %"
        format_sub   = "\t{:<15} \t {:9.2f} s \t {:6.2f} %"

        print()
        print('{:-^58}'.format(' elapsed time '))
        print(format_total.format('total', t_main, t_main/t_main*100))
        print('{:.^58}'.format(''))
        print(format_sub.format('read input files', t_read_files, 
                                    t_read_files/t_main*100))
        print(format_sub.format('process', t_process, 
                                    t_process/t_main*100))
        print(format_sub.format('write output file', t_write_file, 
                                    t_write_file/t_main*100))
        print('{:.^58}'.format(''))
        print()
        
