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

""" Processing PISM output as valid MOM input

usage: ./PISM-to-MOM_processing -o PISM_output_file -e PISM_extra_file
            -m PISM_MOM_mapping_file -a MOM_file -f fluxes_out_file 
            -b basin_shelf_topg_depth_file [-t] [-v]

Mass and energy fluxes to be passed from ice model PISM to ocean model MOM5/6
are computed from PISM extra output file (from multiple variables). Optionally
an ice to ocean reference runoff mass flux is computed from PISM's surface
accumulation to represent the equilibrium runoff under present surface forcing.
Conversion to total fluxes per PISM grid cell is done assuming uniform grid
cell area.  Fluxes are aggregated in PICO basins and distributed to southern
ocean edge cells via previous computed mapping. Fluxes on the ocean grid are
again converted into unit fluxes per area via division with corresponding MOM
grid cell area.  The output file can be used by the FMS data overwrite
mechanism to put PISM fluxes to ocean/sea-ice surface.  Additionally the basin
mean topography depth is computed and stored, which is used to select the depth
of temperature and salinity extraction from 3d ocean output it the
regriddedMOM-to-PISM processing routine. Also optionally the basin mean depth
of ice shelf fronts are computed to determine the input depth of basal melt
fluxes into the ocean model.

Arguments:
    -o, --output PISM_output_file
        PISM output file with flux variables 'mask', 'ice_area_specific_volume' 
        and 'topg'
    -e, --extra-output PISM_extra_file
        input file from PISM extra-output with flux variables
        'surface_runoff_flux', 'tendency_of_ice_amount_due_to_basal_mass_flux' 
        and 'tendency_of_ice_amount_due_to_discharge'.
        Caution: variables contain accumulated fluxes over reporting interval,
        which is the interval of writing out the extra variables. 
    -m, --mapping PISM_MOM_mapping_file
        input file with PICO basin to MOM cell mapping from 
        PISMbasin-to-MOMcell_mapping script
    -a, --area MOM_file
        input file with MOM grid area variable 'area_t'
    -f, --flux-out fluxes_out_file
        file to store processed fluxes which serve as MOM_input
    -b, --topg-depth-out basin_shelf_topg_depth_file
        file to store basin shelf depths which determine vertical layer of 
        ocean boundary condition input to PISM/PICO
    -s, --shelf-depth-out basin_shelf_front_depth_file (optional)
        file to store basin shelf frontal depths which determine the vertical layer of 
        basal melt input into MOM
    --density-ice density_ice (optional)
        set ice density used for frontal shelf depth calculation
    --density-ocean density_ocean (optional)
        set ocean density used for frontal shelf depth calculation
    -r, --runoff-reference-out runoff_reference_file (optional)
        file to store the ice to ocean runoff reference which is computed from
        PISM's surface accumulation and used to calculate the sea level
        changing fraction of the ice to ocean runoff fluxes stored by
        --flux-out
    -t, --time (optional)
        print script time statistics
    -v, --verbose (optional)
        print verbose output
        
This script requires the output of script PISMbasin-to-MOMcell-mapping.py 

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
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
                description=(
"Processing PISM output as valid MOM input. "
"Mass and energy fluxes to be passed from ice model PISM to ocean model MOM5/6 "
"are computed from PISM extra output file (from multiple variables). Optionally "
"an ice to ocean reference runoff mass flux is computed from PISM's surface "
"accumulation to represent the equilibrium runoff under present surface forcing. "
"Conversion to total fluxes per PISM grid cell is done assuming uniform grid "
"cell area.  Fluxes are aggregated in PICO basins and distributed to southern "
"ocean edge cells via previous computed mapping. Fluxes on the ocean grid are "
"again converted into unit fluxes per area via division with corresponding MOM "
"grid cell area.  The output file can be used by the FMS data overwrite "
"mechanism to put PISM fluxes to ocean/sea-ice surface.  Additionally the basin "
"mean topography depth is computed and stored, which is used to select the depth "
"of temperature and salinity extraction from 3d ocean output it the "
"regriddedMOM-to-PISM processing routine. Also optionally the basin mean depth "
"of ice shelf fronts are computed to determine the input depth of basal melt "
"fluxes into the ocean model. "),
                epilog=(
"This script requires the output of script PISMbasin-to-MOMcell-mapping.py "
" "
"This script was created as a processing tool for distributing the flux output "
"of the landice model PISM/PICO to the grid of ocean model MOM5. This was done "
"in the scope of coupling PISM to the climate model POEM at PIK. "))
    parser.add_argument('-o', '--output', 
                        action="store", 
                        dest="PISM_output_file", 
                        required=True, 
                        help=("PISM output file with flux variables "
                              "'mask', 'ice_area_specific_volume' and "
                              "'topg'"))
    parser.add_argument('-e', '--extra-output', 
                        action="store", 
                        dest="PISM_extra_file", 
                        required=True, 
                        help=("PISM output file with flux variables "
                              "'surface_runoff_flux', "
                              "'tendency_of_ice_amount_due_to_basal_mass_flux'" 
                              " and 'tendency_of_ice_amount_due_to_discharge'"))
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
    parser.add_argument('-b', '--topg-depth-out', 
                        action="store", 
                        dest="basin_shelf_topg_depth_file", 
                        required=True, 
                        help=("file to store basin shelf topography depths "
                            "which determine vertical layer of ocean boundary "
                            "condition input to PISM/PICO"))
    parser.add_argument('-s', '--shelf-depth-out', 
                        action="store", 
                        dest="basin_shelf_front_depth_file", 
                        required=False, 
                        help=("file to store basin shelf frontal depths which "
                            "determine the vertical layer basal melt input "
                            "into MOM. Requires PISM extra output variables "
                            "'pism_shelf_mask', 'pism_box_mask' and 'thk'"))
    parser.add_argument('--density-ice', 
                        action="store", 
                        dest="density_ice", 
                        required=False, 
                        help=("set density of ice [kg/m^3] (default: 910). Used "
                              "for computation of frontal shelf depths (when "
                              "--shelf-depth-out is set)."))
    parser.add_argument('--density-ocean', 
                        action="store", 
                        dest="density_ocean", 
                        required=False, 
                        help=("set density of ocean [kg/m^3] (default: 1028). Used "
                              "for computation of frontal shelf depths (when "
                              "--shelf-depth-out is set)."))
    parser.add_argument('-r', '--runoff-reference-out',
                        action="store",
                        dest="runoff_reference_file",
                        required=False,
                        help=("file to store the ice to ocean runoff reference "
                              "which is computed from PISM's surface accumulation "
                              "and used to calculate the sea level changing "
                              "fraction of the ice to ocean runoff fluxes "
                              "stored by --flux-out"))
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
    dens_ice = 910                              # kg/m^3
    dens_ocn = 1028                             # kg/m^3

    if args.density_ice:
        dens_ice = float(args.density_ice)
    if args.density_ocean:
        dens_ocn = float(args.density_ocean)
    
    # a list of possible x,y-dimensions names to read nc-files
    xdims = ['x', 'x1']
    ydims = ['y', 'y1']        
    
    ### ---------- read PISM extra - BEGIN ----------------------------------
    t_read_files_start = time.time()
    
    if args.verbose:
        print(" - reading PISM extra variables from " + args.PISM_extra_file )
    try:
        nc_fh = CDF(args.PISM_extra_file, 'r')
    except:
        s = ("PISM extra file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_extra_file) )
        
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
    pism_extra_time = nc_fh.variables['time'][:]
    pism_extra_time_n = len(pism_extra_time)
    pism_extra_time_dict = nc_fh['time'].__dict__
    if pism_extra_time_n != 1:
        s = ("PISM extra variables have {} timestamps. "
             "Expected: 1")
        raise ValueError( s.format(pism_extra_time_n)) 
 
    #pism_extra_time_bounds = nc_fh.variables['time_bounds'][:]
    #pism_extra_time_bounds_n = len(pism_extra_time_bounds)
    #pism_extra_time_bounds_dict = nc_fh['time_bounds'].__dict__
    
    #pism_extra_nv = nc_fh.dimensions['nv'].size

    #### read PISM variables concerning mass flux from ice to ocean
    ##   all in units: [kg/m^2]
    #pism_bmf = np.squeeze(nc_fh.variables['basal_mass_flux_floating'][:])
    #pism_bmf_dtype = nc_fh.variables['basal_mass_flux_floating'].dtype
    #pism_bmf_ndim = len(pism_bmf.shape)
    #if pism_bmf_ndim != 2:
    #    raise ValueError( str("flux field is of dimension " + \
    #                        str( pism_bmf_ndim ) + ". Expected: 2.") )

    varname = 'tendency_of_ice_amount_due_to_basal_mass_flux'
    pism_tend_bmf = np.squeeze(nc_fh.variables[varname][:])
    pism_tend_bmf_dtype = nc_fh.variables[varname].dtype
    pism_tend_bmf_ndim = len(pism_tend_bmf.shape)
    if pism_tend_bmf_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_tend_bmf_ndim ) + ". Expected: 2.") )

    varname = 'tendency_of_ice_amount_due_to_discharge'
    pism_tend_discharge = np.squeeze(nc_fh.variables[varname][:])
    pism_tend_discharge_dtype = nc_fh.variables[varname].dtype
    pism_tend_discharge_ndim = len(pism_tend_discharge.shape)
    if pism_tend_discharge_ndim != 2:
        raise ValueError( str("flux field is of dimension " + \
                            str( pism_tend_discharge_ndim ) + ". Expected: 2.") )

     # if no surface runoff variable, initialize field with zeros  
    if 'surface_runoff_flux' in nc_fh.variables:
        pism_surf_runoff = np.squeeze(nc_fh.variables['surface_runoff_flux'][:])
        pism_surf_runoff_dtype = nc_fh.variables['surface_runoff_flux'].dtype
        pism_surf_runoff_ndim = len(pism_surf_runoff.shape)
        if pism_surf_runoff_ndim != 2:
            raise ValueError( str("flux field is of dimension " + \
                                str( pism_surf_runoff_ndim ) + ". Expected: 2.") )
    else:
        pism_surf_runoff = np.zeros_like(pism_tend_discharge)

    if args.runoff_reference_file:
        try:
            pism_surf_accum = np.squeeze(nc_fh.variables['surface_accumulation_flux'][:])
            pism_surf_accum_dtype = nc_fh.variables['surface_accumulation_flux'].dtype
            pism_surf_accum_ndim = len(pism_surf_accum.shape)
            if pism_surf_accum_ndim != 2:
                raise ValueError( str("flux field is of dimension " + \
                                str( pism_surf_accum_ndim ) + ". Expected: 2.") )
        except KeyError:
            raise KeyError("PISM's extra output does not have variable "
                "'surface_accumulation_flux', which is needed for calculation "
                "of ice to ocean runoff reference (-r/--runoff-reference-out flag)")


    # read PISM basins
    pism_basins = np.squeeze(nc_fh.variables['basins'][:])
    pism_basins_ndim = len(pism_basins.shape)
    if pism_basins_ndim == 3:
        # cut of time dimension and take first time slice
        pism_basins = pism_basins[0,:,:]
    
    # read PICO contshelf mask
    pism_contshelf_mask = np.squeeze(nc_fh.variables['pico_contshelf_mask'][:])
    pism_contshelf_mask_ndim = len(pism_contshelf_mask.shape)
    if pism_contshelf_mask_ndim == 3:
        # cut of time dimension and take first time slice
        pism_contshelf_mask = pism_contshelf_mask[0,:,:]

    if args.basin_shelf_front_depth_file:
        pico_shelf_mask = np.squeeze(nc_fh.variables['pico_shelf_mask'][:])
        pico_shelf_mask_ndim = len(pico_shelf_mask.shape)
        if pico_shelf_mask_ndim == 3:
            # cut of time dimension and take first time slice
            pico_shelf_mask = pico_shelf_mask[0,:,:]

        pico_box_mask = np.squeeze(nc_fh.variables['pico_box_mask'][:])
        pico_box_mask_ndim = len(pico_box_mask.shape)
        if pico_box_mask_ndim == 3:
            # cut of time dimension and take first time slice
            pico_box_mask = pico_box_mask[0,:,:]

        pism_thk = np.squeeze(nc_fh.variables['thk'][:])
        pism_thk_dtype = nc_fh.variables['thk'].dtype
        pism_thk_ndim = len(pism_thk.shape)
        if pism_thk_ndim == 3:
            # cut of time dimension and take first time slice
            pism_thk = pism_thk[0,:,:]

    # read reporting interval, unit: [years]
    d = nc_fh['pism_config'].__dict__
    pism_extra_times__str = d['output.extra.times']           
    
    nc_fh.close()
    ### ---------- read PISM extra - END -----------------------------------


    ### ---------- read PISM output - BEGIN -----------------------------------   
    if args.verbose:
        print(" - reading PISM output from " + args.PISM_output_file )
    try:
        nc_fh = CDF(args.PISM_output_file, 'r')
    except:
        s = ("PISM output file '{}' can't be found! ")
        raise FileNotFoundError( s.format(args.PISM_output_file) )
        
    pism_mask = np.squeeze(nc_fh.variables['mask'][:])
    #pism_iasv = np.squeeze(nc_fh.variables['ice_area_specific_volume'][:])
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
    pism_extra_times = pism_extra_times__str.split(':')
    if len(pism_extra_times)==3:
        reporting_interval = float(pism_extra_times[1])
    elif len(pism_extra_times)==2:
        reporting_interval = float(pism_extra_times[1]) - float(pism_extra_times[0])
    else:
        s = ("Cannot identify PISM reporting interval! "
             "PISM extra-output time interval has {} items. "
             "Required are 2 or 3.")
        raise( ValueError( s.format(len(pism_extra_times)) ) )
        
    ### calculate cell area
    pism_dx = np.diff(pism_x)[0]    # unit: m
    pism_dy = np.diff(pism_y)[0]    # unit: m
    # uniform area corresponds to PISM internal area representation
    pism_cell_area_uniform = pism_dx*pism_dy    # unit: m^2 
    
    # TODO?: calculate real cell areas based on projection with PROJ4
    # proj4_str = ("+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 "
    #              "+proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0 ")
    
    #### create subset of topography for grid cells either floating or at edge
    ##   of ice shield
    #pism_floating_mask = (pism_mask==3) | (pism_iasv!=0)
    #pism_shelf_topg = np.ma.array(pism_topg, mask=~pism_floating_mask)
    ### create subset of topography for continental shelf grid cells
    #   pism/pico continental shelf mask 
    #       <=> continental shelf (topg above threshold) AND no land ice
    #    0 = False
    #    1 = True, but not relevant
    #    2 = True and relevant
    pism_contshelf_topg = np.ma.array(pism_topg, mask=~(pism_contshelf_mask==2))
    
    
    # aggregate mass from ice to ocean for mass & energy flux calculations
    #  positive corresponds to ice gain
    #  unit[pism_massflux*] = kg/m^2/year
    pism_massflux = {}
    pism_massflux["mass_net"] = -pism_surf_runoff + pism_tend_bmf + \
                                    pism_tend_discharge
    pism_massflux["mass_surf_runoff"] = -pism_surf_runoff
    pism_massflux["mass_basal_melt"] = pism_tend_bmf
    pism_massflux["mass_calving"] = pism_tend_discharge

    pism_massflux_energy = {}
    pism_massflux_energy["energy_net"] = pism_tend_bmf + pism_tend_discharge
    pism_massflux_energy["energy_basal_melt"] = pism_tend_bmf 
    pism_massflux_energy["energy_calving"] = pism_tend_discharge


    ### ------------- conversion of variables ----------------
    
    ### mass flux ice to ocean 
    #  -> unit[pism_massflux] :         kg/m^2/year
    #  -> unit[pism_massflux_total] :   kg/s
    # positive mass flux corresponds to transfer from ice to ocean
    pism_massflux_total = {}.fromkeys(pism_massflux.keys(), None)
    for k in pism_massflux_total.keys():
        pism_massflux_total[k] = -1 * pism_massflux[k] * pism_cell_area_uniform \
                                    / seconds_p_year 
    
    ### heatflux
    #  -> unit[latent_heat_of_fusion] : J/kg
    #  -> unit[pism_heatflux_total] :   J/s = W
    # heat flux PISM to Ocean: should be negative
    pism_heatflux_total = {}.fromkeys(pism_massflux_energy.keys(), None)
    for k in pism_heatflux_total.keys():
        pism_heatflux_total[k] = cp.deepcopy(pism_massflux_energy[k] * pism_cell_area_uniform \
                                    / seconds_p_year * latent_heat_of_fusion)

    if args.runoff_reference_file:
        pism_surf_accum_total = pism_surf_accum * pism_cell_area_uniform \
                                    / seconds_p_year

    # combine dictioniaries
    if args.runoff_reference_file:
        pism_fluxes_total = {**pism_massflux_total,**pism_heatflux_total,
                'surf_accumulation': pism_surf_accum_total}
    else:
        pism_fluxes_total = {**pism_massflux_total,**pism_heatflux_total}
    
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
    pism_basin_flux = {}.fromkeys(pism_fluxes_total.keys())
    for k in pism_basin_flux.keys():
        pism_basin_flux[k] = cp.deepcopy(pism_basin_dummy)

    # create datastructre for basin topography depth
    pism_basin_shelf_topg_depth = np.zeros_like(pism_basin_list, dtype=np.float64)
    pism_basin_shelf_topg_depth[:] = np.nan
    
    
    for idx, val in enumerate(pism_basin_list):
        # cumulate PISM output flux for each basin
        for k in pism_basin_flux.keys():
            pism_basin_flux[k][idx] = np.sum(pism_fluxes_total[k][pism_basins==val])

        # calculate basin mean topography for shelf ice
        basin_mean_depth = np.mean(pism_contshelf_topg[pism_basins==val])
        if (basin_mean_depth is np.ma.masked):
            # default depth in meters
            pism_basin_shelf_topg_depth[idx] = -500
        else:
            pism_basin_shelf_topg_depth[idx] = basin_mean_depth

    # create output structure on MOM grid
    oc_dummy = np.zeros_like(oc_edge_basin, dtype=pism_tend_bmf_dtype)
    oc_edge_flux = {}.fromkeys(pism_fluxes_total.keys())
    for k in oc_edge_flux.keys():
        oc_edge_flux[k] = cp.deepcopy(oc_dummy)
    
    # distribute basin fluxes to MOM cells
    if args.verbose:
        print(" - distributing basin data to ocean grid cells ")
    for j in range(oc_nlat):
        for i in range(oc_nlon):
            if oc_edge_basin.mask[j,i] == False:
                list_index = np.where( pism_basin_list==oc_edge_basin[j,i] )
                for k in oc_edge_flux.keys():
                    oc_edge_flux[k][j,i] = pism_basin_flux[k][list_index] * \
                                            oc_edge_basin_ratio[j,i] 

    # convert fluxes in MOM cells from total to area-relative fluxes
    # units of oc_edge_flux:
    #   mass*:      kg/s
    #   energy*:    J/s = W 
    # units of oc_edge_flux_per_area:
    #   mass*:      kg/s/m^2
    #   energy*:    W/m^2
    oc_edge_flux_per_area = cp.deepcopy(oc_edge_flux)
    for k in oc_edge_flux_per_area.keys():
        oc_edge_flux_per_area[k] /= ocean_area

    # conservation check
    pism_mf_cum = np.sum( np.float128(pism_massflux_total['mass_net']) )
    oc_mf_cum =   np.sum( np.float128(oc_edge_flux['mass_net']) )
    
    pism_ef_cum = np.sum( np.float128(pism_heatflux_total['energy_net']) )
    oc_ef_cum =   np.sum( np.float128(oc_edge_flux['energy_net']) )
    
    error_rate_mass =  (pism_mf_cum - oc_mf_cum) / pism_mf_cum 
    error_rate_energy =  (pism_ef_cum - oc_ef_cum) / pism_ef_cum 

    if args.verbose:
        print(' - relative conservation error')
        print('\tmass: \t\t', error_rate_mass)
        print('\tenergy: \t', error_rate_energy)
    


    ### --------------- calculate frontal ice shelf draft depth ----------------
    if args.basin_shelf_front_depth_file:
        if args.verbose:
            print(" - calculating frontal ice shelf draft depth ")
        shelf_front_box_depth          = np.zeros_like(pism_basins,     dtype=np.float64)
        shelf_front_box_depth_basin    = np.zeros_like(pism_basin_list, dtype=np.float64)
        shelf_front_box_depth_ocean    = np.zeros_like(oc_edge_basin,   dtype=np.float64)
        shelf_front_box_depth[:]       = np.nan
        shelf_front_box_depth_basin[:] = np.nan
        shelf_front_box_depth_ocean[:] = np.nan

        shelf_list = np.unique(pico_shelf_mask)
        shelf_list = shelf_list[shelf_list>0]

        for s in shelf_list:
            # create mask for last PICO box in shelf
            m__shelf = (pico_shelf_mask == s)
            list_boxes = np.unique(pico_box_mask[m__shelf])
            if np.all(list_boxes.mask):
                print(f"    shelf {int(s):4} (cells: {np.sum(m__shelf)}, avg "
                        f"thk: {np.mean(pism_thk[m__shelf]):>6.1f}m) has no "
                        f"PICO box value. Skipping.")
                continue
            else:
                box_max = int(list_boxes.max())
            m__box_max = (pico_box_mask == box_max)           
            m__shelf_box_max = m__shelf & m__box_max
            # compute depth below water surface
            shelf_front_box_depth[m__shelf_box_max] = pism_thk[m__shelf_box_max] * dens_ice/dens_ocn

        # aggregate depths as mean per basin
        for b in pism_basin_list:
            m__basin = (pism_basins == b)
            shelf_front_box_depth_basin[b-1] = np.nanmean(shelf_front_box_depth[m__basin])

        # map basin depth from PISM to MOM grid
        #   multiplied by -1 as positive axis pointing upwards
        for j in range(oc_nlat):
            for i in range(oc_nlon):
                if oc_edge_basin.mask[j,i] == False:
                    list_index = np.where( pism_basin_list==oc_edge_basin[j,i] )
                    shelf_front_box_depth_ocean[j,i] = -1*shelf_front_box_depth_basin[list_index] 

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
        #dst.createDimension('nv', pism_extra_nv)
        
        # write time variable
        dst.createVariable('time', np.double, ("time",) )
        dst['time'].setncatts(pism_extra_time_dict)
        dst['time'][:] = pism_extra_time
        
        #dst.createVariable('time_bounds', np.double, ("time","nv",) )
        #dst['time_bounds'].setncatts(pism_extra_time_bounds_dict)
        #dst['time_bounds'][:] = pism_extra_time_bounds
        
        
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
             
        #### ---- mass flux variables ---- 
        x = dst.createVariable('mass_flux', \
                               nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average mass flux from PISM diagnostic output variables"
                            "'surface_runoff_flux', "
                            "'tendency_of_ice_amount_due_to_basal_mass_flux' and "
                            "'tendency_of_ice_amount_due_to_discharge' "
                            "in reporting interval")),
             ('units', 'kg/m^2/s'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['mass_flux'].setncatts(var_dict)
        dst['mass_flux'][0,:] = oc_edge_flux_per_area['mass_net'][:].data
         
        x = dst.createVariable('mass_flux_surf_runoff', \
                               nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average mass flux from PISM diagnostic output variable "
                            "'surface_runoff_flux_accumulator' "
                            "in reporting interval")),
             ('units', 'kg/m^2/s'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['mass_flux_surf_runoff'].setncatts(var_dict)
        dst['mass_flux_surf_runoff'][0,:] = oc_edge_flux_per_area['mass_surf_runoff'][:].data

        x = dst.createVariable('mass_flux_basal_melt', \
                               nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average mass flux from PISM diagnostic output variable "
                            "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator' "
                            "in reporting interval")),
             ('units', 'kg/m^2/s'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['mass_flux_basal_melt'].setncatts(var_dict)
        dst['mass_flux_basal_melt'][0,:] = oc_edge_flux_per_area['mass_basal_melt'][:].data

        x = dst.createVariable('mass_flux_calving', \
                               nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average mass flux from PISM diagnostic output variable "
                            "'tendency_of_ice_amount_due_to_discharge_accumulator' "
                            "in reporting interval")),
             ('units', 'kg/m^2/s'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['mass_flux_calving'].setncatts(var_dict)
        dst['mass_flux_calving'][0,:] = oc_edge_flux_per_area['mass_calving'][:].data


        #### ---- heat flux variables ---- 
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
        dst['heat_flux'][0,:] = oc_edge_flux_per_area['energy_net'][:].data

        x = dst.createVariable('heat_flux_basal_melt', nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average heat flux calculated from PISM diagnostic output variable "
                            "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator' "
                            "in reporting interval")),
             ('units', 'W/m^2'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['heat_flux_basal_melt'].setncatts(var_dict)        
        dst['heat_flux_basal_melt'][0,:] = oc_edge_flux_per_area['energy_basal_melt'][:].data

        x = dst.createVariable('heat_flux_calving', nc_dtype, ('time','yt_ocean','xt_ocean'))
        var_dict = col.OrderedDict([
             ('long_name', ("average heat flux calculated from PISM diagnostic output variable "
                            " 'tendency_of_ice_amount_due_to_discharge_accumulator' "
                            "in reporting interval")),
             ('units', 'W/m^2'),
             ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
             ('reporting_interval', reporting_interval ),
             ('reporting_interval_units', 'years'),
             ('cell_methods', 'time: point'),
             ('coordinates', 'geolon_t geolat_t')])
        dst['heat_flux_calving'].setncatts(var_dict)
        dst['heat_flux_calving'][0,:] = oc_edge_flux_per_area['energy_calving'][:].data

    # copy area variable to file
    with CDF(args.MOM_file, 'r') as src,   \
         CDF(args.PISM_to_MOM_fluxes_file, "a") as dst:
        # copy area_t variable
        for name, variable in src.variables.items():
            if name in ['area_t']:
                #code.interact(local=locals())
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]


    ### ------------------ save runoff reference to file -----------------------
    #   write redistributed ice to ocean reference runoff (surface accumulation)
    #   to file runoff_reference_file
    if args.runoff_reference_file:
        if args.verbose:
            print(" - write reference runoff to file ",
                      args.runoff_reference_file)
        with CDF(args.PISM_MOM_mapping_file, 'r') as src,   \
             CDF(args.runoff_reference_file, "w") as dst:
            # copy global attributes all at once via dictionary
            glob_dict = src.__dict__
            glob_dict['filename'] = os.path.basename(args.runoff_reference_file)
            glob_dict['title'] = "ice to ocean reference runoff fluxes computed " \
                "from PISM extra output variable 'surface_accumulation_flux' " \
                "and redistributed to MOM grid"

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
            #dst.createDimension('nv', pism_extra_nv)

            # write time variable
            dst.createVariable('time', np.double, ("time",) )
            dst['time'].setncatts(pism_extra_time_dict)
            dst['time'][:] = pism_extra_time

            #dst.createVariable('time_bounds', np.double, ("time","nv",) )
            #dst['time_bounds'].setncatts(pism_extra_time_bounds_dict)
            #dst['time_bounds'][:] = pism_extra_time_bounds


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
            if pism_surf_accum_dtype == 'float32':
                nc_dtype = 'f4'
            elif pism_surf_accum_dtype == 'float64':
                nc_dtype = 'f8'
            else:
                s = 'pism_surf_accum_dtype is "{}". Only "float32" and "float64" are allowed.'
                raise ValueError(s.format(pism_tend_bmf_dtype))

            #### ---- runoff reference variable ----
            x = dst.createVariable('mass_flux', nc_dtype, \
                                   ('time','yt_ocean','xt_ocean'))
            var_dict = col.OrderedDict([
                 ('long_name', ("ice to ocean reference runoff flux computed "
                                "as the mean of PISM extra output variable "
                                "'surface_accumulation_flux' over reporting "
                                "interval and redistributed to MOM grid")),
                 ('units', 'kg/m^2/s'),
                 ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
                 ('reporting_interval', reporting_interval ),
                 ('reporting_interval_units', 'years'),
                 ('cell_methods', 'time: point'),
                 ('coordinates', 'geolon_t geolat_t')])
            dst['mass_flux'].setncatts(var_dict)
            dst['mass_flux'][0,:] = oc_edge_flux_per_area['surf_accumulation'][:].data


        with CDF(args.MOM_file, 'r') as src,   \
             CDF(args.runoff_reference_file, "a") as dst:
            # copy area_t variable
            for name, variable in src.variables.items():
                if name in ['area_t']:
                    #code.interact(local=locals())
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name].setncatts(src[name].__dict__)
                    dst[name][:] = src[name][:]


    ### ---------------------- save basin topography depths file ---------------------------
    #   write basin mean topography of ice shelf areas to basin_shelf_topg_depth_file 
    if args.verbose:
        print(" - write basin mean topography of ice shelf areas to file ",
                  args.basin_shelf_topg_depth_file)
    
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "
   
    with CDF(args.basin_shelf_topg_depth_file, "w") as dst:
        # create dictionary for global attributes
        glob_dict = dict({  \
            'filename': os.path.basename( args.basin_shelf_topg_depth_file),
            'title':    'basin mean topography of ice shelf areas',
            'history':  histstr
            })
        dst.setncatts(glob_dict)
        
        # create dimensions
        dst.createDimension('time', None)
        dst.createDimension('n_basin', len(pism_basin_list) )

        # create variables
        x = dst.createVariable('basin', int, 'n_basin')
        var_dict = col.OrderedDict([
             ('name', "PISM-PICO basin"),
             ('long_name', "list of valid PISM-PICO basins")])
        dst['basin'].setncatts(var_dict)      
        dst['basin'][:] = pism_basin_list[:]

        dst.createVariable('time', np.double, ("time",) )
        dst['time'].setncatts(pism_extra_time_dict)
        dst['time'][:] = pism_extra_time

        x = dst.createVariable('mean_shelf_topg', float, ('time','n_basin'))
        var_dict = col.OrderedDict([
             ('long_name', "mean basin topography of ice shelf areas"),
             ('units', 'm'),
             ('axis', 'Z'),
             ('positive', 'up'),
             ('fill_value', netCDF4._netCDF4.default_fillvals['f4'])])
        dst['mean_shelf_topg'].setncatts(var_dict)
        dst['mean_shelf_topg'][0,:] = pism_basin_shelf_topg_depth[:]


    ### ---------------------- save shelf depth file ---------------------------
    #   write mean ice shelf frontal depth to basin_shelf_front_depth_file
    if args.basin_shelf_front_depth_file:
        if args.verbose:
            print(" - write basin mean shelf front depth to file ",
                      args.basin_shelf_front_depth_file)
        with CDF(args.PISM_MOM_mapping_file, 'r') as src,   \
             CDF(args.basin_shelf_front_depth_file, "w") as dst:
            # copy global attributes all at once via dictionary
            glob_dict = src.__dict__
            glob_dict['filename'] = os.path.basename(args.basin_shelf_front_depth_file)
            glob_dict['title'] = "ice shelf front depth computed as last PICO" \
                "box mean" 

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
            #dst.createDimension('nv', pism_extra_nv)

            # write time variable
            dst.createVariable('time', np.double, ("time",) )
            dst['time'].setncatts(pism_extra_time_dict)
            dst['time'][:] = pism_extra_time

            #dst.createVariable('time_bounds', np.double, ("time","nv",) )
            #dst['time_bounds'].setncatts(pism_extra_time_bounds_dict)
            #dst['time_bounds'][:] = pism_extra_time_bounds


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
            if pism_thk_dtype == 'float32':
                nc_dtype = 'f4'
            elif pism_thk_dtype == 'float64':
                nc_dtype = 'f8'
            else:
                raise ValueError(f'pism_thk_dtype is {pism_thk_dtype}. '
                        'Only "float32" and "float64" are allowed.')

            x = dst.createVariable('shelf_front_depth', nc_dtype, 
                                   ('time','yt_ocean','xt_ocean'))
            var_dict = col.OrderedDict([
                 ('long_name', ("front depth of ice shelf draft, computed as "
                                "the average over the last PICO box in every "
                                "shelf, aggregated as mean for each basin")),
                 ('units', 'm'),
                 ('fill_value', netCDF4._netCDF4.default_fillvals[nc_dtype]),
                 ('reporting_interval', reporting_interval ),
                 ('reporting_interval_units', 'years'),
                 ('cell_methods', 'time: point'),
                 ('coordinates', 'geolon_t geolat_t')])
            dst['shelf_front_depth'].setncatts(var_dict)
            dst['shelf_front_depth'][0,:] = shelf_front_box_depth_ocean[:].data


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

