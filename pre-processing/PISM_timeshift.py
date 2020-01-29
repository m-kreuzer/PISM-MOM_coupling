#!/usr/bin/env python3

""" Shift timestamp of PISM restart file to align with POEM/MOM timestamp.

usage: ./PISM_timeshift.py -p pism_restart_file -c climate_restart_file \
        -s pism_restart_timeshift_path_file [-t] [-v]

The restart file of ice model PISM is modified such as the timestamp is the 
same as the corresponding climate model (POEM/MOM). To do so the 
pism_restart_file and the climate_restart_file have to be given as arguments.
The name of the modified restart file is written to the 
pism_restart_timeshift_path_file.

Arguments:
    -p pism_restart_file
        a netCDF file with PISM restart conditions 
    -c climate_restart_file
        file named coupler.res from INPUT/ directory in GFDL ocean/climate
        model structure 
    -s pism_restart_timeshift_path_file
        file name to which the path of modified PISM restart file is written
    -t (optional)
        print script time statistics
    -v (optional)
        print verbose output
        
This script was created as a preprocessing tool for coupling the ice 
sheet model PISM/PICO to the climate model POEM at PIK.

"""

import sys
import os
import time as t
import argparse
import numpy as np
import collections as col
try:
    from netCDF4 import Dataset as CDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

__author__ = "Moritz Kreuzer"
__copyright__ = "Copyright 2020"
__credits__ = ["", ""]
__license__ = "GPLv3"
__version__ = "0.0.3"
__maintainer__ = "Moritz Kreuzer"
__email__ = "kreuzer@pik-potsdam.de"
__status__ = "Prototype"



if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(
            description=
                "Shift timestamp of PISM restart file to align with POEM/MOM    \
                timestamp",
                epilog=
                "The restart file of ice model PISM is modified such as the     \
                timestamp is the same as the corresponding climate model        \
                (POEM/MOM). To do so the pism_restart_file and the              \
                climate_restart_file have to be given as arguments.             \
                The name & path of the modified restart file is written to the  \
                pism_restart_timeshift_path_file."
            )


    parser.add_argument('-p', '--pism-restart', action="store", 
                        dest="pism_restart_file", required=True, 
                        help="a netCDF file with PISM restart conditions")
    parser.add_argument('-c', '--climate-restart', action="store", 
                        dest="climate_restart_file", required=True, 
                        help="file named coupler.res from INPUT/ directory in \
                              GFDL ocean/climate model structure")
    parser.add_argument('-s', '--pism-restart-timeshift-path', action="store", 
                        dest="pism_restart_timeshift_path_file", required=True, 
                        help="file name to which the path of modified PISM \
                                restart file is written")
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
        
    # --------------- read PISM restart file ----------------
    #   -> read PISM restart file to extract timestamp
    if args.verbose:
        print("... reading PISM restart file '" + args.pism_restart_file + "'")
    t_read_pismrestartfile_start = t.time()    

    try:
        nc_pism = CDF(args.pism_restart_file, 'r')
    except:
        print("pism_restart_file '" + args.pism_restart_file + \
                "' can't be found! Exiting.")
        sys.exit(1)
        
    # read field array 
    try:
        nc_pism.variables['time'].set_auto_mask(False)
        pism_time_raw = np.squeeze(nc_pism.variables['time'][:])
        #pism_time__sec = pism_time_raw.data
        pism_time__sec = pism_time_raw
        print("pism_time__sec: ", pism_time__sec)
        print("pism_time_raw: ", pism_time_raw)
        print("type(pism_time_raw): ", type(pism_time_raw))
    except:
        print("Variable 'time' can't be read from file '" +
                args.pism_restart_file + "'!")

    # check field dimension
    n_pism_time = pism_time__sec.size
    if (n_pism_time != 1) :
        err_str = "PISM restart file has more than one timestamp. Ambigous!"
        raise ValueError( str(err_str) )

    nc_pism.close()
    t_read_pismrestartfile_end = t.time()

    # --------------- read MOM restart file ----------------
    #   -> read MOM restart file to extract timestamp
    if args.verbose:
        print("... reading MOM restart file '" + args.climate_restart_file + "'")
    t_read_climrestartfile_start = t.time()    

    try:
        fh_clim_time = open(args.climate_restart_file, 'r')
    except:
        print("climate_restart_file '" + args.climate_restart_file + \
                "' can't be found! Exiting.")
        sys.exit(1)
        
    # read field array 
    try:
        mom_time_file = fh_clim_time.readlines()
        tmp_str = mom_time_file[2]     # extract 3rd line 
        mom_time_array_str = np.array(' '.join(tmp_str.split()).split(' ')[0:5])
    except:
        print("POEM/MOM timestamp can't be read from file '" +
                args.climate_restart_file + "'!")

    # check time stamp content 
    for s in mom_time_array_str:
        if( s.isdigit() == False) :
            err_str = "Time stamp from climate model does not only contain numbers!"
            raise ValueError( str(err_str) )

    fh_clim_time.close()
    t_read_climrestartfile_end = t.time()

    # -------------------- calculate MOM timestamp in seconds --------------------  
    if args.verbose:
        print("... caluating MOM timestamp in seconds")
    t_MOM_time_sec_start = t.time()

    monthly_days = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    seconds_per_year = 365*24*60*60
    seconds_per_day = 24*60*60
    seconds_per_hour = 60*60

    mom_time_array = mom_time_array_str.astype(float)
    mom_yrs, mom_months, mom_days, mom_hours, mom_secs = mom_time_array

    if args.verbose:
        print()
        print("--- MOM time stamp ---")
        print(" year:  \t",mom_yrs)
        print(" month: \t",mom_months)
        print(" day:   \t",mom_days)
        print(" hours: \t",mom_hours)
        print(" secs:  \t",mom_secs)
        print()
    
    # calculate passed seconds for running year, month, day and hour
    mom_time_yrs__sec = mom_yrs * seconds_per_year 
    if( mom_months > 1): 
        mom_time_months__sec = (monthly_days[0:int(mom_months)-1].sum() ) \
                                    * seconds_per_day
    else:
        mom_time_months__sec = 0
    mom_time_days__sec = (mom_days - 1) * seconds_per_day
    mom_time_hours__sec = mom_hours * seconds_per_hour
    mom_time_secs__sec = mom_secs

    #print(" years_sec: \t",mom_time_yrs__sec)
    #print(" months_sec:\t",mom_time_months__sec)
    #print(" days_sec:  \t",mom_time_days__sec)
    #print(" hours_sec: \t",mom_time_hours__sec)
    #print(" secs_sec:  \t",mom_time_secs__sec)
    #print()
    mom_time__sec = mom_time_yrs__sec + mom_time_months__sec + \
            mom_time_days__sec + mom_time_hours__sec + \
            mom_time_secs__sec

    if args.verbose:
        print("time stamps [seconds since 0000/1/1 - 00:00:00]:")
        print("-----------------------------------------------")
        print(" {:<10}{:>15.1f} s".format("MOM", mom_time__sec))
        print(" {:<10}{:>15.1f} s".format("PISM", pism_time__sec))
        print()

    t_MOM_time_sec_end = t.time()

    # -------------------- do timeshift of PISM restart file --------------------  
    if args.verbose:
        print("... shifting timestamp of PISM restart file with NCO ")
    t_timeshift_start = t.time()

    time_shift__sec = mom_time__sec - pism_time__sec
    if args.verbose:
        print(" PISM time shift (sec): ", time_shift__sec )
        print()

    # MOM time stamp in format YYYYMMDD
    MOM_timestamp_print ='{:04n}{:02n}{:02n}'.format(mom_time_array[0],\
                                                     mom_time_array[1],\
                                                     mom_time_array[2])

    # create new file name including MOM timestamp
    pism_restart_file_path, pism_restart_file_extension = \
            os.path.splitext(args.pism_restart_file)
    pism_restart_file_shift = pism_restart_file_path + '.' + MOM_timestamp_print \
                                + pism_restart_file_extension
    
    # use NCO to shift time variable
    try:
        os.system('ncap2 -s time+=' + str(time_shift__sec) + ' ' \
                + args.pism_restart_file + ' -O -o ' + pism_restart_file_shift )
    except:
        print("ncap2 command failed! Exiting.")
        sys.exit(1)


    t_timeshift_end = t.time()

    # -------------- write location of shifted PISM restart file -------------- 
    if args.verbose:
        print("... writing location of shifted PISM restart file to '" + \
                args.pism_restart_timeshift_path_file + "' ")
    t_write_outfile_start = t.time()

    try:
        fh_outfile = open(args.pism_restart_timeshift_path_file, 'w')
    except:
        print("Cant't write to pism_restart_timeshift_path_file '" + \
                args.pism_restart_timeshift_path_file + " Exiting.")
        sys.exit(1)


    fh_outfile.write("# relative path to time shifted PISM restart file written by PISM_timeshift.py \n")
    fh_outfile.write(pism_restart_file_shift)
    fh_outfile.close()


    t_write_outfile_end = t.time()
    t_main_end = t.time()

    # -------------------- performance -------------------- 

    if args.verbose | args.time:
        t_main              = t_main_end            - t_main_start
        t_read_pismrestartfile   = t_read_pismrestartfile_end     - \
                                                t_read_pismrestartfile_start
        t_read_climrestartfile   = t_read_climrestartfile_end     - \
                                                t_read_climrestartfile_start
        t_MOM_time_sec      = t_MOM_time_sec_end    - t_MOM_time_sec_start
        t_timeshift         = t_timeshift_end       - t_timeshift_start
        t_write_outfile     = t_write_outfile_end   - t_write_outfile_start

        format_total = "{:<35}{:9.2f} s \t {:6.2f} %"
        format_sub   = "  {:<33}{:9.2f} s \t {:6.2f} %"

        print()
        print('{:-^58}'.format(' elapsed time '))
        print(format_total.format('total', t_main, t_main/t_main*100))
        print('{:.^58}'.format(''))
        print(format_sub.format('read PISM restart file', t_read_pismrestartfile, 
                                    t_read_pismrestartfile/t_main*100))
        print(format_sub.format('read MOM restart file', t_read_climrestartfile, 
                                    t_read_climrestartfile/t_main*100))
        print(format_sub.format('calculate MOM timestamp in sec', t_MOM_time_sec, 
                                    t_MOM_time_sec/t_main*100))
        print(format_sub.format('timeshift PISM restart file', t_timeshift, 
                                    t_timeshift/t_main*100))
        print(format_sub.format('write outfile', t_write_outfile, 
                                    t_write_outfile/t_main*100))
        print('{:.^58}'.format(''))
        print()

