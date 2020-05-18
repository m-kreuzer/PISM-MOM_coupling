#!/usr/bin/env python3

""" Redistributing mass and energy fluxes from PISM/PICO to MOM cells.

useage: ./PISM-to-MOM_processing -o PISM_output_file -e PISM_snap_file 
            -m PISM_MOM_mapping_file -a MOM_file -f fluxes_out_file 
            -d basin_shelf_depth_file [-t] [-v]

Mass and energy fluxes are computed from PISM snapshot file (from multiple 
variables). Conversion to total fluxes per PISM grid cell is done assuming 
uniform grid cell area.
Fluxes are aggregated in PICO basins and distributed to southern ocean edge 
cells via previous computed mapping. Fluxes on the ocean grid are again converted
into unit fluxes per area via division with corresponding MOM grid cell area.
The output file can be used by the FMS data overwrite mechanism to put PISM 
fluxes to ocean/sea-ice surface.

Arguments:
    -o PISM_output_file
        PISM output file with flux variables 'mask', 'ice_area_specific_volume' 
        and 'topg'
    -e PISM_snap_file
        input file from PISM snapshot file with flux variables
        'surface_runoff_flux_accumulator', 
        'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator' 
        and 'tendency_of_ice_amount_due_to_discharge_accumulator'.
        Caution: variables '*_accumulator' contain accumulated fluxes over 
        interval specified in '*_time_since_reset'. 
    -m PISM_MOM_mapping_file
        input file with PICO basin to MOM cell mapping from 
        PISMbasin-to-MOMcell_mapping script
    -a MOM_file
        input file with MOM grid area variable 'area_t'
    -f fluxes_out_file
        file to store processed fluxes which serve as MOM_input
    -d basin_shelf_depth_file
        file to store basin shelf depths which determine vertical layer of 
        ocean boundary condition input to PISM/PICO
    -t (optional)
        print script time statistics
    -v (optional)
        print verbose output
        
This script requires the ouput of script PISMbasin-to-MOMcell-mapping.py 

This script was created as a processing tool for distributing the flux output
of the landice model PISM/PICO to the grid of ocean model MOM5. This was done
in the scope of coupling PISM to the climate model POEM at PIK.

"""

import sys
import os
import numpy as np
import copy as cp
import collections as col
import time
import argparse
try:
    import netCDF4
    from netCDF4 import Dataset as CDF
except:
    raise ImportError("netCDF4 is not installed!")
    
    
    
__author__ = "Moritz Kreuzer"
__copyright__ = "Copyright 2019"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.2"
__maintainer__ = "Moritz Kreuzer"
__email__ = "kreuzer@pik-potsdam.de"
__status__ = "Prototype"


    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                description=
                "Redistributing mass and energy fluxes from PISM/PICO to MOM cells.",
                epilog=
                ("Mass and energy fluxes are computed from PISM output file. "
                 "Conversion to total fluxes per PISM grid cell is done "
                 "assuming uniform grid cell area. Fluxes are aggregated in "
                 "PICO basins and distributed to southern ocean edge cells "
                 "via previous computed mapping. The output file can be used "
                 "by the FMS data overwrite mechanism.")
            )
    parser.add_argument('-o', '--output', 
                        action="store", 
                        dest="PISM_output_file", 
                        required=True, 
                        help=("PISM output file with flux variables "
                              "'mask', 'ice_area_specific_volume' and "
                              "'topg'"))
    parser.add_argument('-e', '--snap_output', 
                        action="store", 
                        dest="PISM_snap_file", 
                        required=True, 
                        help=("PISM snapshot output file with flux variables "
                              "'surface_runoff_flux_accumulator', "
                              "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'" 
                              " and 'tendency_of_ice_amount_due_to_discharge_accumulator'"))
    parser.add_argument('-m', '--mapping', 
                        action="store", 
                        dest="PISM_MOM_mapping_file",
                        required=True, 
                        help=("file with PICO basin to MOM cell mapping from "
                              "PISM_basin_to_MOM_cell_mapping script"))
    parser.add_argument('-a', '--area', 
                        action="store", 
                        dest="MOM_file",
                        required=True, 
                        help="MOM file with area variable 'area_t'")
    parser.add_argument('-f', '--flux-out', 
                        action="store", 
                        dest="PISM_to_MOM_fluxes_file", 
                        required=True, 
                        help=("file to store processed fluxes which serve as "
                              "MOM input"))
    parser.add_argument('-d', '--depth-out', 
                        action="store", 
                        dest="basin_shelf_depth_file", 
                        required=True, 
                        help=("file to store basin shelf depths which determine"
                              "vertical layer of ocean boundary condition input"
                              "to PISM/PICO"))
    parser.add_argument('-t', '--time',
                        action="store_true", 
                        help="print script timings")
    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        help="increase output verbosity")
    args = parser.parse_args()    
    
    
    # -------------------- general setup --------------------------
    
    t_main_start = time.time()
     
    if args.verbose:
        print("Running", sys.argv[0])
        print(" -> verbose output = True")
        print()
    
    
    # constant definitions
    seconds_p_year = 365*24*60*60               # for no-leap year calendar!
    latent_heat_of_fusion = 3.34 * 1e5          # unit: J/kg   
        # -> see: Reese et al 2018: Antarctic sub-shelf melt rates via PICO
        # consistent with HLF in mom5.0.2/src/shared/constants/constants.F90
    
    # a list of possible x,y-dimensions names to read nc-files
    xdims = ['x', 'x1']
    ydims = ['y', 'y1']        
    
    ### ---------- read PISM snapshot - BEGIN ----------------------------------
    t_read_files_start = time.time()
    
    if args.verbose:
        print(" - reading PISM snapshot variables from " + args.PISM_snap_file )
    try:
        nc_fh = CDF(args.PISM_snap_file, 'r')
    except:
        s = ("PISM snapshot file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_snap_file) )
        
    # assign x,y dimension
    for dim in xdims:
        if dim in list(nc_fh.dimensions.keys()):
            xdim = dim
    for dim in ydims:
        if dim in list(nc_fh.dimensions.keys()):
            ydim = dim
    
    # coordinate variable in x,y-direction
    pism_x = nc_fh.variables[xdim][:]
    pism_y = nc_fh.variables[ydim][:]
    #pism_lat = nc_fh.variables['lat'][:]
    #pism_lon = nc_fh.variables['lon'][:]
    ## transform ocean longitudes to range [-180, 180] degE
    #basin_lon2 = copy.deepcopy(basin_lon)
    #basin_lon2[basin_lon2 < -180] +=360
    #basin_lon2[basin_lon2 >  180] -=360

    # read time axis incl dimensions and variables
    pism_snap_time = nc_fh.variables['time'][:]
    pism_snap_time_n = len(pism_snap_time)
    pism_snap_time_dict = nc_fh['time'].__dict__
    if pism_snap_time_n != 2:
        s = ("PISM snap variables have {} timestamps. "
             "Expected: 2")
        raise ValueError( s.format(pism_snap_time_n)) 
 
    #pism_snap_time_bounds = nc_fh.variables['time_bounds'][:]
    #pism_snap_time_bounds_n = len(pism_snap_time_bounds)
    #pism_snap_time_bounds_dict = nc_fh['time_bounds'].__dict__
    
    #pism_snap_nv = nc_fh.dimensions['nv'].size

    ### read PISM variables concerning mass flux from ice to ocean
    #   all in units: [kg/m^2]
    pism_bmf = np.squeeze(nc_fh.variables['basal_mass_flux_floating_accumulator'][1])
    pism_bmf_dtype = nc_fh.variables['basal_mass_flux_floating_accumulator'].dtype
    pism_bmf_ndim = len(pism_bmf.shape)
    if pism_bmf_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_bmf_ndim ) + ". Expected: 2.") )

    pism_surf_runoff = np.squeeze(nc_fh.variables['surface_runoff_flux_accumulator'][1])
    pism_surf_runoff_dtype = nc_fh.variables['surface_runoff_flux_accumulator'].dtype
    pism_surf_runoff_ndim = len(pism_surf_runoff.shape)
    if pism_surf_runoff_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_surf_runoff_ndim ) + ". Expected: 2.") )

    pism_tend_bmf = np.squeeze(nc_fh.variables['tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'][1])
    pism_tend_bmf_dtype = nc_fh.variables['tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'].dtype
    pism_tend_bmf_ndim = len(pism_tend_bmf.shape)
    if pism_tend_bmf_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_tend_bmf_ndim ) + ". Expected: 2.") )

    pism_tend_discharge = np.squeeze(nc_fh.variables['tendency_of_ice_amount_due_to_discharge_accumulator'][1])
    pism_tend_discharge_dtype = nc_fh.variables['tendency_of_ice_amount_due_to_discharge_accumulator'].dtype
    pism_tend_discharge_ndim = len(pism_tend_discharge.shape)
    if pism_tend_discharge_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_tend_discharge_ndim ) + ". Expected: 2.") )

    ### read time span for PISM accumulation variables in snapshot file
    #   time in unit: [s]
    pism_snaptime_raw = np.squeeze(nc_fh.variables['basal_mass_flux_floating_time_since_reset'][:])
    pism_snaptime_dtype = nc_fh.variables['basal_mass_flux_floating_time_since_reset'].dtype
    pism_snaptime_ndim = len(pism_snaptime_raw.shape)
    if pism_snaptime_ndim != 1:
        raise ValueError( str("snapshot time field is of dimension " + \
                            str( pism_snaptime_ndim ) + ". Expected: 1.") )
    # select second timestep of two (beginning and end of run)
    if len(pism_snaptime_raw) != 2:
        raise ValueError( str("snapshot time field has " + \
                            str( len(pism_snaptime_raw) ) + "entries. Expected: 2.") )
    else:
        pism_snaptime = pism_snaptime_raw[1]

    # read PISM basins
    pism_basins = np.squeeze(nc_fh.variables['basins'][:])
    pism_basins_ndim = len(pism_basins.shape)
    if pism_basins_ndim == 3:
        # cut of time dimension and take first time slice
        pism_basins = pism_basins[0,:,:]
    
    # read reporting interval, unit: [years]
    d = nc_fh['pism_config'].__dict__
    pism_snap_times__str = d['output.snapshot.times']           
    
    nc_fh.close()
    ### ---------- read PISM snapshot - END -----------------------------------


    ### ---------- read PISM output - BEGIN -----------------------------------   
    if args.verbose:
        print(" - reading PISM output from " + args.PISM_output_file )
    try:
        nc_fh = CDF(args.PISM_output_file, 'r')
    except:
        s = ("PISM output file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_output_file) )
        
    pism_mask = np.squeeze(nc_fh.variables['mask'][:])
    pism_iasv = np.squeeze(nc_fh.variables['ice_area_specific_volume'][:])
    pism_topg = np.squeeze(nc_fh.variables['topg'][:])
    
    nc_fh.close()
    ### ---------- read PISM output - END ------------------------------------
           
    ### ---------- read PISM-MOM mapping - BEGIN -----------------------------
    if args.verbose:
        print(" - reading PISM to MOM mapping file " + args.PISM_MOM_mapping_file )
    try:
        nc_fh = CDF(args.PISM_MOM_mapping_file, 'r')
    except:
        s = ("PISM to MOM mapping file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_MOM_mapping_file) )
        

    # deactivate mask for NC file input
    # (wrong valid range for longitude leads to missing values)    
    nc_fh.set_auto_mask(False)
    
    # coordinate variable in x,y-direction
    ocean_x =   nc_fh.variables['xt_ocean'][:]
    ocean_y =   nc_fh.variables['yt_ocean'][:]
    ocean_lat = nc_fh.variables['geolat_t'][:]
    ocean_lon = nc_fh.variables['geolon_t'][:]
    #ocean_area = nc_fh.variables['area_t'][:]
    
    # activate mask for NC file input again
    nc_fh.set_auto_mask(True)
    
#    # shift ocean longitudes to range [-180, 180] degE
#    ocean_lon_s = copy.deepcopy(ocean_lon)
#    ocean_lon_s[ocean_lon_s < -180] +=360
#    ocean_lon_s[ocean_lon_s >  180] -=360

    # read basin info arrays    
    oc_edge_basin = np.squeeze(nc_fh.variables['basin'][:])
    ocean_ndim = len(oc_edge_basin.shape)
    if ocean_ndim != 2:
        s = ("Variable 'south_edge_basin' from file '{}' is of dimension {}."
             "Expected: 2")
        raise ValueError( s.format(ocean_ndim)) 
    
    oc_edge_basin_ratio = np.squeeze(nc_fh.variables['basin_ratio'][:])
    ocean_ndim = len(oc_edge_basin_ratio.shape)
    if ocean_ndim != 2:
        s = ("Variable 'oc_edge_basin_ratio' is of dimension {}. Expected: 2")
        raise ValueError( s.format(ocean_ndim)) 
        
    oc_nlat = oc_edge_basin.shape[0]
    oc_nlon = oc_edge_basin.shape[1]
    
    nc_fh.close()
    ### ---------- read PISM-MOM mapping - END -------------------------------
   
    ### ------------- read MOM area - BEGIN ---------------------------------
    if args.verbose:
        print(" - reading MOM file " + args.MOM_file )
    try:
        nc_fh = CDF(args.MOM_file, 'r')
    except:
        s = ("MOM file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.MOM_file) )
        
    # deactivate mask for NC file input
    # (wrong valid range for longitude leads to missing values)    
    nc_fh.set_auto_mask(False)
    
    # coordinate variable in x,y-direction
    ocean_x2 =   nc_fh.variables['xt_ocean'][:]
    ocean_y2 =   nc_fh.variables['yt_ocean'][:]
    ocean_lat2 = nc_fh.variables['geolat_t'][:]
    ocean_lon2 = nc_fh.variables['geolon_t'][:]
    # read area variable
    ocean_area = nc_fh.variables['area_t'][:]   # units: m^2
    
    # activate mask for NC file input again
    nc_fh.set_auto_mask(True)
    
    assert_str = ("non matching {{}}-coordinates between PISM-to-MOM-mapping "
                  "file '{}' and MOM file '{}'."
                  ).format(args.PISM_MOM_mapping_file, args.MOM_file)
    assert (ocean_x == ocean_x2).all(), assert_str.format('x')
    assert (ocean_y == ocean_y2).all(), assert_str.format('y')
    assert (ocean_lat == ocean_lat2).all(), assert_str.format('latitude')
    assert (ocean_lon == ocean_lon2).all(), assert_str.format('longitude')  
    
    nc_fh.close()
    t_read_files_end = time.time()
    ### --------------- read MOM area - END ----------------------------------
    
    ### ---------- start general processing ----------------------------------
    t_process_start = time.time()
    
    ### extract reporting interval of PISM flux 
    #   -> time over which flux was aggregated, unit: years
    pism_snap_times = pism_snap_times__str.split(':')
    if len(pism_snap_times)==3:
        reporting_interval = float(pism_snap_times[1])
    elif len(pism_snap_times)==2:
        reporting_interval = float(pism_snap_times[1]) - float(pism_snap_times[0])
    elif len(pism_snap_times)==1:
        reporting_interval = float(pism_snap_times[0])
    else:
        s = ("Cannot identify PISM reporting interval! "
             "PISM snapshot-output time interval has {} items. "
             "Required are 1, 2 or 3.")
        raise( ValueError( s.format(len(pism_snap_times)) ) )
        
    ### calculate cell area
    pism_dx = np.diff(pism_x)[0]    # unit: m
    pism_dy = np.diff(pism_y)[0]    # unit: m
    # uniform area corresponds to PISM internal area representation
    pism_cell_area_uniform = pism_dx*pism_dy    # unit: m^2 
    
    # TODO?: calculate real cell areas based on projection with PROJ4
    # proj4_str = ("+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 "
    #              "+proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0 ")
    
    ### create subset of topography for grid cells either floating or at edge
    #   of ice shield
    pism_floating_mask = (pism_mask==3) | (pism_iasv!=0)
    pism_shelf_topg = np.ma.array(pism_topg, mask=~pism_floating_mask)
    
    
    # aggregate mass from ice to ocean for mass & energy flux calculations
    #  positive corresponds to ice gain
    #  unit[pism_massflux*] = kg/m^2
    pism_massflux = -pism_surf_runoff + pism_tend_bmf + pism_tend_discharge
    pism_massflux_energy = pism_tend_bmf + pism_tend_discharge


    ### ------------- conversion of variables ----------------
    
    ### mass flux ice to ocean 
    #  -> unit[pism_massflux] :         kg/m^2
    #  -> unit[pism_massflux_total] :   kg/s
    # positive mass flux corresponds to transfer from ice to ocean
    pism_massflux_total = -1 * pism_massflux * pism_cell_area_uniform \
                                    / pism_snaptime
    
    ### heatflux
    #  -> unit[latent_heat_of_fusion] : J/kg
    #  -> unit[pism_heatflux_total] :   J/s = W
    # heat flux PISM to Ocean: should be negative
    pism_heatflux_total = pism_massflux_energy * pism_cell_area_uniform \
                            / pism_snaptime * latent_heat_of_fusion
     

    
    ### ------------- aggregation of fluxes per basin ----------------
    if args.verbose:
        print(" - aggregating PISM output to basins ")
    
    # create list of PISM basins
    pism_basin_list_tmp = np.unique(pism_basins)
    pism_basin_list = pism_basin_list_tmp[~pism_basin_list_tmp.mask].data
    # remove basin 0
    pism_basin_list = np.delete(pism_basin_list, np.where(pism_basin_list==0) ) 
    # make sure datatype is integer
    pism_basin_list = pism_basin_list.astype(int)
    
    # create datastructre for basin cumulated fluxes
    pism_basin_dummy = np.zeros_like(pism_basin_list, dtype=pism_tend_bmf_dtype)
    pism_basin_flux = dict({'mass_total':  cp.deepcopy(pism_basin_dummy),   \
                            'heat_total':  cp.deepcopy(pism_basin_dummy) })
    # create datastructre for basin topography depth
    pism_basin_shelf_depth = np.zeros_like(pism_basin_list, dtype=np.float64)
    pism_basin_shelf_depth[:] = np.nan
    
    
    for idx, val in enumerate(pism_basin_list):
         # cumulate PISM output flux for each basin
        pism_basin_flux['mass_total'][idx] = \
                                    np.sum(pism_massflux_total[pism_basins==val])
        pism_basin_flux['heat_total'][idx] = \
                                    np.sum(pism_heatflux_total[pism_basins==val])
        # calculate basin mean topography for shelf ice
        basin_mean_depth = np.mean(pism_shelf_topg[pism_basins==val])
        if (basin_mean_depth is np.ma.masked):
            # default depth in meters
            pism_basin_shelf_depth[idx] = -500
        else:
            pism_basin_shelf_depth[idx] = np.mean(pism_shelf_topg[pism_basins==val])

                            
    # create output structure on MOM grid
    oc_dummy = np.zeros_like(oc_edge_basin, dtype=pism_tend_bmf_dtype)
    oc_edge_flux = dict({'mass_total': cp.deepcopy(oc_dummy),  \
                         'mass':       cp.deepcopy(oc_dummy),  \
                         'heat_total': cp.deepcopy(oc_dummy),  \
                         'heat':       cp.deepcopy(oc_dummy) })
    
    # distribute basin fluxes to MOM cells
    if args.verbose:
        print(" - distributing basin data to ocean grid cells ")
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_edge_basin.mask[j,i] == False:
                list_index = np.where( pism_basin_list==oc_edge_basin[j,i] )
                oc_edge_flux['mass_total'][j,i] = \
                                    pism_basin_flux['mass_total'][list_index] * \
                                    oc_edge_basin_ratio[j,i] 
                oc_edge_flux['heat_total'][j,i] = \
                                    pism_basin_flux['heat_total'][list_index] * \
                                    oc_edge_basin_ratio[j,i] 


    # convert fluxes in MOM cells from total to area-relative fluxes
    # units of oc_edge_flux['*']:
    #   mass_total:     kg/s
    #   mass:           kg/s/m^2
    #   heat_total:     J/s = W 
    #   heat:           W/m^2
    oc_edge_flux['mass'] = oc_edge_flux['mass_total'] / ocean_area
    oc_edge_flux['heat'] = oc_edge_flux['heat_total'] / ocean_area
    
    # conservation check
    pism_mf_cum = np.sum( np.float128(pism_massflux_total) )
    oc_mf_cum =   np.sum( np.float128(oc_edge_flux['mass_total']) )
    
    pism_ef_cum = np.sum( np.float128(pism_heatflux_total) )
    oc_ef_cum =   np.sum( np.float128(oc_edge_flux['heat_total']) )
    
    error_rate_mass =  (pism_mf_cum - oc_mf_cum) / pism_mf_cum 
    error_rate_energy =  (pism_ef_cum - oc_ef_cum) / pism_ef_cum 

    if args.verbose:
        print(' - relative conservation error')
        print('\tmass: \t\t', error_rate_mass)
        print('\tenergy: \t', error_rate_energy)
    
    
    t_process_end = time.time()
    
    ### ---------------------- save fluxes to file ---------------------------
    #   write redistributed flux variables to file PISM_to_MOM_fluxes_file
    t_write_file_start = time.time()
    
    if args.verbose:
        print(" - write fluxes on ocean grid to file ", args.PISM_to_MOM_fluxes_file)

    dim_copy = ['xt_ocean','yt_ocean']
    var_copy = ['geolat_t', 'geolon_t', 'xt_ocean', 'yt_ocean']
    
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "

    
    with CDF(args.PISM_MOM_mapping_file, 'r') as src,   \
         CDF(args.PISM_to_MOM_fluxes_file, "w") as dst:
        # copy global attributes all at once via dictionary
        glob_dict = src.__dict__
        glob_dict['filename'] = os.path.basename(args.PISM_to_MOM_fluxes_file)
        glob_dict['title'] = 'heat and mass fluxes from PISM on MOM grid'
        
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
                
        # create time dimension
        dst.createDimension('time', None)
        #dst.createDimension('nv', pism_snap_nv)
        
        # write time variable
        dst.createVariable('time', np.double, ("time",) )
        dst['time'].setncatts(pism_snap_time_dict)
        dst['time'][:] = pism_snap_time[1]
        
        #dst.createVariable('time_bounds', np.double, ("time","nv",) )
        #dst['time_bounds'].setncatts(pism_snap_time_bounds_dict)
        #dst['time_bounds'][:] = pism_snap_time_bounds
        
        
        # copy variables
        for name, variable in src.variables.items():
            if name in var_copy:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
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
        if pism_tend_bmf_dtype == 'float32':
            nc_dtype = 'f4'
        elif pism_tend_bmf_dtype == 'float64':
            nc_dtype = 'f8'
        else:
            s = 'pism_tend_bmf_dtype is "{}". Only "float32" and "float64" are allowed.'
            raise ValueError(s.format(pism_tend_bmf_dtype))
             
        x = dst.createVariable('mass_flux', \
                               nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average mass flux from PISM diagnostic output variables"
                            "'surface_runoff_flux_accumulator', "
                            "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'"
                            " and 'tendency_of_ice_amount_due_to_discharge_accumulator' "
                            "in reporting interval")),
             ('units', 'kg/m^2/s'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['mass_flux'].setncatts(var_dict)
        dst['mass_flux'][0,:] = oc_edge_flux['mass'][:].data
         
        x = dst.createVariable('heat_flux', nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average heat flux calculated from PISM diagnostic output variables"
                            "'surface_runoff_flux_accumulator', "
                            "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'"
                            " and 'tendency_of_ice_amount_due_to_discharge_accumulator' "
                            "in reporting interval")),
             ('units', 'W/m^2'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['heat_flux'].setncatts(var_dict)        
        dst['heat_flux'][0,:] = oc_edge_flux['heat'][:].data


    ### ---------------------- save basin depths file ---------------------------
    #   write basin mean topography of ice shelf areas to basin_shelf_depth_file 
    if args.verbose:
        print(" - write basin mean topography of ice shelf areas to file ",
                  args.basin_shelf_depth_file)
    
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "
   
    with CDF(args.basin_shelf_depth_file, "w") as dst:
        # create dictionary for global attributes
        glob_dict = dict({  \
            'filename': os.path.basename( args.basin_shelf_depth_file),
            'title':    'basin mean topography of ice shelf areas',
            'history':  histstr
            })
        dst.setncatts(glob_dict)
        
        # create basin dimension
        dst.createDimension('n_basin', len(pism_basin_list) )

        # create variables
        x = dst.createVariable('basin', int, 'n_basin')
        var_dict = col.OrderedDict([
             ('long_name', "list of valid PISM/PICO basins")])
        dst['basin'].setncatts(var_dict)      
        dst['basin'][:] = pism_basin_list[:]

        x = dst.createVariable('mean_shelf_topg', float, 'n_basin')
        var_dict = col.OrderedDict([
             ('long_name', "mean basin topography of ice shelf areas"),
             ('units', 'm'),
             ('axis', 'Z'),
             ('positive', 'up'),
             ('fill_value', netCDF4._netCDF4.default_fillvals['f4'])])
        dst['mean_shelf_topg'].setncatts(var_dict)
        dst['mean_shelf_topg'][:] = pism_basin_shelf_depth[:]


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

