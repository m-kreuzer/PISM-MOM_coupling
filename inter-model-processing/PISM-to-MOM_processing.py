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
of temperature and salinity extraction from 3d ocean output by the
regriddedMOM-to-PISM processing routine. Also optionally the basin mean depth
of ice shelf fronts are computed to determine the input depth of basal melt
fluxes into the ocean model.

Arguments:
    -e, --extra-output PISM_extra_file
        input file from PISM extra-output with flux variables 'mask', 'topg',
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
import time
import argparse
from tqdm import tqdm
import xarray as xr
import warnings
# debug
import code

def check_ndims(da, valid_ndims=[2,3], ds_source=None):
    '''checks whether xr.DataArray has the expected dimensions 
    and raises error instead'''
    
    n_dims = len(da.dims)
    ndims_match = False
    for d in valid_ndims:
        if d==n_dims:
            ndims_match = True
            break
    
    if ndims_match==False:
        if type(ds_source)==str:
            err_path = f"(path: '{ds_source}') "
        else:
            err_path = ''
        err = f"variable '{da.name}' {err_path}has {n_dims} dimension(s)."
        if len(valid_ndims) == 1:
            err += f" Expected: {valid_ndims[0]}"
        else:
            err += f" Expected one of: {valid_ndims}"
        raise ValueError(err)

def check_ndims_ds(ds, check_dims_dict):
    '''checks whether a number of variables in xr.Dataset have expected dimensions 
    
    ds: xr.Dataset,
        dataset to check for given variables
    check_dims_dict: dict,
        dictionary holding variable names with list of valid number of dimensions to check
    '''
    for v, valid_ndims in check_dims.items():
        try:
            check_ndims(ds[v], 
                        valid_ndims=valid_ndims, 
                        ds_source=ds.encoding['source'])
        except ValueError:
            raise
        except KeyError:
            warnings.warn(f"variable '{v}' not found in '{ds.encoding['source']}'.")
    

def check_var_match(ds1, ds2, var):
    '''checks whether given variable is identical between datasets'''
    assert_str = f"non matching variable '{var}' between " \
                 f"file '{ds1.encoding['source']}' and file '{ds2.encoding['source']}'."
    assert (ds1[var] == ds2[var]).all(), assert_str

### https://gist.github.com/bsolomon1124/44f77ed2f15062c614ef6e102bc683a5
class DeprecateAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        warnings.warn("Argument %s is deprecated and is *ignored*." % self.option_strings)
        delattr(namespace, self.dest)

def mark_deprecated_help_strings(parser, prefix="DEPRECATED"):
    for action in parser._actions:
        if isinstance(action, DeprecateAction):
            h = action.help
            if h is None:
                action.help = prefix
            else:
                action.help = prefix + ": " + h
    

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
"of temperature and salinity extraction from 3d ocean output by the "
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
                        #action="store", 
                        action=DeprecateAction,
                        dest="PISM_output_file", 
                        required=True, 
                        help=("Using PISM output variables 'mask' "
                              "and 'topg' now from PISM extra output file"))
    parser.add_argument('-e', '--extra-output', 
                        action="store", 
                        dest="PISM_extra_file", 
                        required=True, 
                        help=("PISM output file with flux variables "
                              "'mask', 'topg', "
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
    
    mark_deprecated_help_strings(parser)
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
    #xdims = ['x', 'x1']
    #ydims = ['y', 'y1']        
    
    ### ---------- read PISM extra - BEGIN ----------------------------------
    t_read_files_start = time.time()
    
    if args.verbose:
        print(" - reading PISM extra variables from " + args.PISM_extra_file )
    
    pism_extra = xr.open_dataset(args.PISM_extra_file, use_cftime=True)
    
    check_dims = {'tendency_of_ice_amount_due_to_basal_mass_flux':[3],
                  'tendency_of_ice_amount_due_to_discharge':[3],
                  'surface_runoff_flux':[3],
                  'surface_accumulation_flux':[3],
                  'basins':[3],
                  'pico_contshelf_mask':[3],
                  'pico_shelf_mask':[3],
                  'pico_box_mask':[3],
                  'thk':[3],
                 }
    
    check_ndims_ds(pism_extra, check_dims_dict=check_dims)
              
    if args.runoff_reference_file:
        if 'surface_accumulation_flux' not in pism_extra.variables:
            raise KeyError("variable 'surface_accumulation_flux' not found in PISM's "
                           f"extra file (path '{pism_extra.encoding['source']}'). "
                           "Required for calculation of ice to ocean runoff reference "
                           "(-r/--runoff-reference-out flag).")
            
    ### ---------- read PISM extra - END -----------------------------------

    # ### ---------- read PISM output - BEGIN -----------------------------------   
    # DEPRECATED, using 'mask' and 'topg' variables now from PISM extra output
    # if args.verbose:
    #     print(" - reading PISM output from " + args.PISM_output_file )
    # pism_out = xr.open_dataset(args.PISM_output_file, use_cftime=True)
    # ### ---------- read PISM output - END ------------------------------------

    ### ---------- read PISM-MOM mapping - BEGIN -----------------------------
    if args.verbose:
        print(" - reading PISM to MOM mapping file " + args.PISM_MOM_mapping_file )
    
    pism_mom_mapping = xr.open_dataset(args.PISM_MOM_mapping_file, use_cftime=True)
    
    check_dims = {'basin':[2],
                  'basin_ratio':[2],
                 }
    
    check_ndims_ds(pism_mom_mapping, check_dims_dict=check_dims)
    ### ---------- read PISM-MOM mapping - END -------------------------------
    
    
    ### ------------- read MOM area - BEGIN ---------------------------------
    if args.verbose:
        print(" - reading MOM file " + args.MOM_file )
    
    mom_out = xr.open_dataset(args.MOM_file, use_cftime=True)
    
    check_var_match(mom_out, pism_mom_mapping, 'xt_ocean')
    check_var_match(mom_out, pism_mom_mapping, 'yt_ocean')
    check_var_match(mom_out, pism_mom_mapping, 'geolat_t')
    check_var_match(mom_out, pism_mom_mapping, 'geolon_t')
    
    t_read_files_end = time.time()
    ### --------------- read MOM area - END ----------------------------------
    
    
    ### ---------- start general processing ----------------------------------
    t_process_start = time.time()
    
    
    ### calculate cell area
    pism_dx = pism_extra['x'].diff(dim='x').data[0]    # unit: m
    pism_dy = pism_extra['y'].diff(dim='y').data[0]    # unit: m
    # uniform area corresponds to PISM internal area representation
    pism_cell_area_uniform = pism_dx*pism_dy    # unit: m^2 
    
    # TODO?: calculate real cell areas based on projection with PROJ4
    # proj4_str = ("+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 "
    #              "+proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0 ")
    
    ### create subset of topography for continental shelf grid cells
    #   pism/pico continental shelf mask 
    #       <=> continental shelf (topg above threshold) AND no land ice
    #    0 = False
    #    1 = True, but not relevant
    #    2 = True and relevant
    pism_contshelf_topg = pism_extra['topg'].where(pism_extra['pico_contshelf_mask']==2, np.nan)
    
    # aggregate mass from ice to ocean for mass & energy flux calculations
    #  positive corresponds to ice gain
    #  unit[pism_massflux*] = kg/m^2/year
    
    pism_tend_bmf        = pism_extra['tendency_of_ice_amount_due_to_basal_mass_flux']
    pism_tend_discharge  = pism_extra['tendency_of_ice_amount_due_to_discharge']
    if 'surface_runoff_flux' in pism_extra.variables:
        pism_surf_runoff = pism_extra['surface_runoff_flux']
    else:
        # initialise field with 0
        pism_surf_runoff = pism_extra['tendency_of_ice_amount_due_to_discharge'] * 0
        pism_surf_runoff.attrs = {'long_name': 'surface runoff, averaged over the reporting interval',
                                  'standard_name': 'surface_runoff_flux',
                                  'units': pism_extra['tendency_of_ice_amount_due_to_discharge'].attrs['units']}
        pism_surf_runoff = pism_surf_runoff.rename('surface_runoff_flux')
        
    if args.runoff_reference_file:
        pism_reference_runoff = -pism_extra['surface_accumulation_flux']
    
    pism_flux = xr.Dataset()
    pism_flux["mass_flux"] = -pism_surf_runoff + pism_tend_bmf + \
                             pism_tend_discharge
    pism_flux["mass_flux"].attrs['units'] = pism_tend_bmf.attrs['units']
    pism_flux["mass_flux"].attrs['long_name'] = \
        ("average mass flux from PISM diagnostic output variables"
         "'surface_runoff_flux', 'tendency_of_ice_amount_due_to_basal_mass_flux' "
         "and 'tendency_of_ice_amount_due_to_discharge' in reporting interval")
        
    pism_flux["mass_flux_surf_runoff"] = -pism_surf_runoff
    pism_flux["mass_flux_surf_runoff"].attrs['long_name'] = \
        ("average mass flux from PISM diagnostic output variable "
         "'surface_runoff_flux_accumulator' in reporting interval")
    
    pism_flux["mass_flux_basal_melt"] = pism_tend_bmf
    pism_flux["mass_flux_basal_melt"].attrs['long_name'] = \
        ("average mass flux from PISM diagnostic output variable "
         "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator' "
         "in reporting interval")
    
    pism_flux["mass_flux_calving"] = pism_tend_discharge
    pism_flux["mass_flux_calving"].attrs['long_name'] = \
        ("average mass flux from PISM diagnostic output variable "
         "'tendency_of_ice_amount_due_to_discharge_accumulator' "
         "in reporting interval")
    
    if args.runoff_reference_file:
        pism_flux["mass_ref_runoff"] = pism_reference_runoff
        pism_flux["mass_ref_runoff"].attrs['long_name'] = \
            ("ice to ocean reference runoff flux computed as the mean of "
             "PISM extra output variable 'surface_accumulation_flux' over "
             "reporting interval and redistributed to MOM grid")
    
    pism_flux["heat_flux"] = pism_tend_bmf + pism_tend_discharge
    pism_flux["heat_flux"].attrs['units'] = pism_tend_bmf.attrs['units']
    pism_flux["heat_flux"].attrs['long_name'] = \
        ("average heat flux calculated from PISM diagnostic output variables"
         "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator'"
         " and 'tendency_of_ice_amount_due_to_discharge_accumulator' "
         "in reporting interval")
    
    pism_flux["heat_flux_basal_melt"] = pism_tend_bmf 
    pism_flux["heat_flux_basal_melt"].attrs['long_name'] = \
        ("average heat flux calculated from PISM diagnostic output variable "
         "'tendency_of_ice_amount_due_to_basal_mass_flux_accumulator' "
         "in reporting interval")
    
    pism_flux["heat_flux_calving"] = pism_tend_discharge
    pism_flux["heat_flux_calving"].attrs['long_name'] = \
        ("average heat flux calculated from PISM diagnostic output variable "
         " 'tendency_of_ice_amount_due_to_discharge_accumulator' "
         "in reporting interval")
    
    ### ------------- conversion of variables ----------------
    for v in pism_flux.variables:
        if 'mass' in v:
            assert pism_flux[v].attrs['units'] == "kg m-2 year-1", \
                    f"variable {v} is not in units 'kg m-2 year-1'"
            attrs_sav = cp.deepcopy(pism_flux[v].attrs)
            pism_flux[v] = -1 * pism_flux[v] * pism_cell_area_uniform / seconds_p_year
            pism_flux[v].attrs = attrs_sav
            pism_flux[v].attrs['units'] = "kg/s"
            pism_flux[v].attrs['sign'] = \
                "positive mass flux corresponds to transfer from ice to ocean"
            
        if 'heat' in v:
            assert pism_flux[v].attrs['units'] == "kg m-2 year-1", \
                    f"variable {v} is not in units 'kg m-2 year-1'"
            attrs_sav = cp.deepcopy(pism_flux[v].attrs)
            pism_flux[v] = pism_flux[v] * pism_cell_area_uniform / seconds_p_year \
                            * latent_heat_of_fusion
            pism_flux[v].attrs = attrs_sav
            pism_flux[v].attrs['units'] = "W"
            pism_flux[v].attrs['sign'] = \
                "positive heat flux corresponds to adding heat to ocean"
            
        if 'coordinates' in pism_flux[v].attrs:
            pism_flux[v].attrs.pop('coordinates')
    
    ### ------------- aggregation of fluxes per basin ----------------
    if args.verbose:
        print(" - aggregating PISM output to basins ")
    
    # cumulate PISM output flux for each basin
    pism_basin_flux_time_list = []
    for t in tqdm(pism_extra.time, desc='cumulating PISM output basin wise (iterating time)', disable=not(args.verbose)):
        ds_tmp = pism_flux.sel(time=t).groupby(pism_extra['basins'].sel(time=t)).sum(..., keep_attrs=True)
        pism_basin_flux_time_list.append(ds_tmp.expand_dims(dim='time'))
    pism_basin_flux = xr.merge(pism_basin_flux_time_list)
    
    # calculate basin mean topography for shelf ice
    #  and fill NaNs with default depth (in meters)
    basin_mean_depth_time_list = []
    for t in tqdm(pism_extra.time, desc='calculate basin mean topography for shelf ice (iterating time)', disable=not(args.verbose)):
        pism_basins = pism_extra['basins'].sel(time=t)
        pism_basins_nan = pism_basins.where(pism_basins!=0,np.nan)
        da_tmp = pism_contshelf_topg.sel(time=t).groupby(pism_basins_nan).mean(keep_attrs=True)
        da_tmp = da_tmp.fillna(-500)
        basin_mean_depth_time_list.append(da_tmp.expand_dims(dim='time'))
    pism_basin_shelf_topg_depth = xr.merge(basin_mean_depth_time_list)
    
    pism_basin_shelf_topg_depth = pism_basin_shelf_topg_depth.rename_dims({'basins':'n_basin'})
    pism_basin_shelf_topg_depth = pism_basin_shelf_topg_depth.rename({'basins':'basin'})
    pism_basin_shelf_topg_depth['basin'].attrs['name']      = 'PISM-PICO basin'
    pism_basin_shelf_topg_depth['basin'].attrs['long_name'] = 'list of valid PISM-PICO basins'
    
    pism_basin_shelf_topg_depth = pism_basin_shelf_topg_depth.rename({'topg':'mean_shelf_topg'})
    attrs = {}
    attrs['units'] = pism_basin_shelf_topg_depth['mean_shelf_topg'].attrs['units']
    attrs['long_name'] = 'mean basin topography of ice shelf areas'
    attrs['axis'] = 'Z'
    attrs['positive'] = 'up'
    pism_basin_shelf_topg_depth['mean_shelf_topg'].attrs = attrs
    
    
    # distribute basin fluxes to MOM cells
    oc_edge_flux = xr.Dataset()
    
    
    pism_mom_mapping_basin_time = pism_mom_mapping['basin'].broadcast_like(pism_extra.time)
    
    for v in pism_basin_flux.data_vars:
        oc_edge_flux[v] = pism_mom_mapping_basin_time * np.nan
    
    for v in oc_edge_flux.data_vars:
        for b in pism_basin_flux.basins:
            oc_edge_flux[v].data = xr.where(pism_mom_mapping_basin_time==b, 
                                            pism_basin_flux[v].sel(basins=b), 
                                            oc_edge_flux[v]).data
    
        oc_edge_flux[v] *= pism_mom_mapping['basin_ratio']
        oc_edge_flux[v].attrs = pism_basin_flux[v].attrs
    
        
    # remove unwanted basin coordinate
    #oc_edge_flux = oc_edge_flux.drop('basins')
    
    # convert fluxes in MOM cells from total to area-relative fluxes
    assert mom_out['area_t'].attrs['units']=='m^2', \
            f"variable 'area_t' is not in units 'm^2'"
    oc_edge_flux_per_area = oc_edge_flux / mom_out['area_t']
    
    ### ------------- conversion of variables ----------------
    for v in oc_edge_flux_per_area.variables:
        #print(v)
        if 'mass' in v:
            oc_edge_flux_per_area[v].attrs = oc_edge_flux[v].attrs
            assert oc_edge_flux_per_area[v].attrs['units'] == "kg/s", \
                    f"variable {v} is not in units 'kg/s'"
            oc_edge_flux_per_area[v].attrs['units'] = "kg/s/m^2"
            
        if 'heat' in v:
            oc_edge_flux_per_area[v].attrs = oc_edge_flux[v].attrs
            assert oc_edge_flux_per_area[v].attrs['units'] == "W", \
                    f"variable {v} is not in units 'W'"
            oc_edge_flux_per_area[v].attrs['units'] = "W/m^2"
    
    oc_edge_flux_per_area['area_t'] = mom_out['area_t']
    
    # conservation check
    pism_mf_cum = pism_flux['mass_flux'].astype(np.float128).sum(dim=('x','y'))
    oc_mf_cum   = oc_edge_flux['mass_flux'].astype(np.float128).sum(dim=('xt_ocean','yt_ocean'))
    
    pism_ef_cum = pism_flux['heat_flux'].astype(np.float128).sum(dim=('x','y'))
    oc_ef_cum   = oc_edge_flux['heat_flux'].astype(np.float128).sum(dim=('xt_ocean','yt_ocean'))
    
    error_rate_mass   =  (pism_mf_cum - oc_mf_cum) / pism_mf_cum 
    error_rate_energy =  (pism_ef_cum - oc_ef_cum) / pism_ef_cum 
    
    if args.verbose:
        print(' - max. relative conservation error')
        print(f'\tmass:   {np.max(np.abs(error_rate_mass)).data:4.2e}')
        print(f'\tenergy: {np.max(np.abs(error_rate_energy)).data:4.2e}')
    
    ## --------------- calculate frontal ice shelf draft depth ----------------
    if args.basin_shelf_front_depth_file:
        if args.verbose:
            print(" - calculating frontal ice shelf draft depth ")
    
            shelf_front_box_depth = pism_extra['pico_box_mask'] * np.nan
    
        for t in tqdm(pism_extra.time, desc='on ice grid (iterating time)', disable=not(args.verbose)):
            pico_shelf_mask = pism_extra['pico_shelf_mask'].sel(time=t)
            pico_box_mask   = pism_extra['pico_box_mask'].sel(time=t)
            pism_thk        = pism_extra['thk'].sel(time=t)
    
            shelf_list = np.unique(pico_shelf_mask)
            shelf_list = shelf_list[shelf_list>0]
    
            for s in shelf_list:
                m__shelf = (pico_shelf_mask==s)
                shelf_boxes = pico_box_mask.where(m__shelf)
                list_boxes = np.unique(shelf_boxes)
                if np.all(np.isnan(list_boxes)):
                    if args.verbose:
                        n_cells = np.sum(~np.isnan(shelf_boxes)).data
                        thk_mean = np.mean(pism_extra['thk'].where(~np.isnan(shelf_boxes), np.nan)).data
                        if n_cells >0:
                            print(f"    {t.dt.strftime('%Y-%M-%d').data}, shelf {s:6.3f} (cells: {n_cells}, avg " 
                                    f"thk: {thk_mean:>6.1f}m) has no PICO box value. Skipping.")
                    continue
                else:
                    box_max = int(np.nanmax(list_boxes))
    
                m__box_max = (pico_box_mask == box_max)
    
                m__shelf_box_max = (m__shelf & m__box_max)
                m__depth = (pism_thk.where(m__shelf_box_max, np.nan) * dens_ice/dens_ocn).to_masked_array()
    
                shelf_front_box_depth.sel(time=t).data[m__shelf_box_max] = m__depth[~m__depth.mask]
    
                
        # aggregate depths as mean per basin        
        shelf_front_box_depth_basin_time_list = []
        for t in tqdm(pism_extra.time, desc='calculating mean depths per basin', disable=not(args.verbose)):
            da_tmp = shelf_front_box_depth.sel(time=t).groupby(pism_extra['basins'].sel(time=t)).mean(...)
            shelf_front_box_depth_basin_time_list.append(da_tmp.expand_dims(dim='time'))
        shelf_front_box_depth_basin = xr.merge(shelf_front_box_depth_basin_time_list).rename({'pico_box_mask':'shelf_front_depth'})
    
    
        attrs = {}
        attrs['units'] = pism_thk.attrs['units']
        attrs['long_name'] = ("front depth of ice shelf draft, computed as "
                              "the average over the last PICO box in every "
                              "shelf, aggregated as mean for each basin")
        attrs['axis'] = 'Z'
        attrs['positive'] = 'down'
        shelf_front_box_depth_basin.attrs = attrs
    
        # map basin depth from PISM to MOM grid
        shelf_front_box_depth_ocean = xr.Dataset()
    
        pism_mom_mapping_basin_time = pism_mom_mapping['basin'].broadcast_like(pism_extra.time)
        v = 'shelf_front_depth'
        shelf_front_box_depth_ocean[v] = pism_mom_mapping_basin_time * np.nan
    
        for b in tqdm(shelf_front_box_depth_basin.basins, desc='redistribution on ocean grid (iterating basins)', disable=not(args.verbose)):
            shelf_front_box_depth_ocean[v] = xr.where(pism_mom_mapping_basin_time==b, 
                                                           shelf_front_box_depth_basin[v].sel(basins=b), 
                                                           shelf_front_box_depth_ocean[v])
    
        shelf_front_box_depth_ocean['time'].attrs = pism_extra.time.attrs
    
        # remove unwanted basin coordinate
        #shelf_front_box_depth_ocean = shelf_front_box_depth_ocean.drop('basins')
    
    t_process_end = time.time()
    
    ### ---------------------- save fluxes to file ---------------------------
    #   write redistributed flux variables to file PISM_to_MOM_fluxes_file
    t_write_file_start = time.time()
    
    if 'mass_ref_runoff' in oc_edge_flux_per_area:
        # move reference runoff to different Dataset
        mass_ref_runoff = xr.Dataset()
        mass_ref_runoff['mass_flux'] = oc_edge_flux_per_area['mass_ref_runoff']
        mass_ref_runoff['area_t']    = oc_edge_flux_per_area['area_t']
        oc_edge_flux_per_area = oc_edge_flux_per_area.drop('mass_ref_runoff')
    
        
    if args.verbose:
        print(" - write fluxes on ocean grid to file ", args.PISM_to_MOM_fluxes_file)
    
    cmd_line = ' '.join(sys.argv)
    histstr = time.asctime() + ': ' + cmd_line + "\n "
    
    glob_attrs = cp.deepcopy(pism_mom_mapping.attrs)
    glob_attrs['filename'] = os.path.basename(args.PISM_to_MOM_fluxes_file)
    glob_attrs['title'] = 'heat and mass fluxes from PISM on MOM grid'
    
    if 'history' in glob_attrs.keys():
        glob_attrs['history'] = histstr + glob_attrs['history']
    elif 'History' in glob_attrs.keys():
        glob_attrs['History'] = histstr + glob_attrs['History']
    else:
        glob_attrs['history'] = histstr
    
    oc_edge_flux_per_area.attrs = glob_attrs
    oc_edge_flux_per_area.to_netcdf(args.PISM_to_MOM_fluxes_file,
                                   encoding={"time":      {"dtype": "float", "units":f"days since 0001-01-01 00:00:00"}})
    
    ### ------------------ save runoff reference to file -----------------------
    #   write redistributed ice to ocean reference runoff (surface accumulation)
    #   to file runoff_reference_file
    if args.runoff_reference_file:
        if args.verbose:
            print(" - write reference runoff to file ",
                      args.runoff_reference_file)
            
    
        glob_attrs = cp.deepcopy(pism_mom_mapping.attrs)
        glob_attrs['filename'] = os.path.basename(args.runoff_reference_file)
        glob_attrs['title'] =  "ice to ocean reference runoff fluxes computed " \
                    "from PISM extra output variable 'surface_accumulation_flux' " \
                    "and redistributed to MOM grid"
    
        if 'history' in glob_attrs.keys():
            glob_attrs['history'] = histstr + glob_attrs['history']
        elif 'History' in glob_attrs.keys():
            glob_attrs['History'] = histstr + glob_attrs['History']
        else:
            glob_attrs['history'] = histstr
    
        mass_ref_runoff.attrs = glob_attrs
        mass_ref_runoff.to_netcdf(args.runoff_reference_file,
                                  encoding={"time":      {"dtype": "float", "units":f"days since 0001-01-01 00:00:00"}})
    
    ### ---------------------- save basin topography depths file ---------------------------
    #   write basin mean topography of ice shelf areas to basin_shelf_topg_depth_file 
    if args.verbose:
        print(" - write basin mean topography of ice shelf areas to file ",
                  args.basin_shelf_topg_depth_file)
        
    glob_attrs = {}
    glob_attrs['filename'] = os.path.basename(args.basin_shelf_topg_depth_file)
    glob_attrs['title']    =  'basin mean topography of ice shelf areas'
    glob_attrs['history']  = histstr
    
    pism_basin_shelf_topg_depth.attrs = glob_attrs
    pism_basin_shelf_topg_depth.to_netcdf(args.basin_shelf_topg_depth_file,
                                          encoding={"time":      {"dtype": "float", "units":f"days since 0001-01-01 00:00:00"}})
    
    
    ### ---------------------- save shelf depth file ---------------------------
    #   write mean ice shelf frontal depth to basin_shelf_front_depth_file
    if args.basin_shelf_front_depth_file:
        if args.verbose:
            print(" - write basin mean shelf front depth to file ",
                      args.basin_shelf_front_depth_file)
    
        glob_attrs = cp.deepcopy(pism_mom_mapping.attrs)
        glob_attrs['filename'] = os.path.basename(args.basin_shelf_front_depth_file)
        glob_attrs['title'] =  "ice shelf front depth computed as last PICO" \
                    "box mean" 
    
        if 'history' in glob_attrs.keys():
            glob_attrs['history'] = histstr + glob_attrs['history']
        elif 'History' in glob_attrs.keys():
            glob_attrs['History'] = histstr + glob_attrs['History']
        else:
            glob_attrs['history'] = histstr
    
        shelf_front_box_depth_ocean.attrs = glob_attrs
        shelf_front_box_depth_ocean.to_netcdf(args.basin_shelf_front_depth_file,
                                                   encoding={"time":      {"dtype": "float", "units":f"days since 0001-01-01 00:00:00"}})
        
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
    
