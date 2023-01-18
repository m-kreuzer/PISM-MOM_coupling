
import sys
import os
import stat
import shutil
import jinja2
import collections
import distutils.dir_util as dist
import subprocess
import warnings
import glob

import helpers 
import settings 


from contextlib import contextmanager

# redirect stderr/stdout functions copied from
# https://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python/22434262#22434262
def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

def merged_stderr_stdout():  # $ exec 2>&1
    return stdout_redirected(to=sys.stdout, stdout=sys.stderr)


def create_script_from_template(settings, template_file,
                      experiment=settings.experiment):

    if not os.path.exists(settings.experiment_dir):
        os.makedirs(settings.experiment_dir)

    # make jinja aware of templates
    template_path = os.path.join(settings.project_root,"templates")
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_path),
                                   trim_blocks=True, 
                                   lstrip_blocks=True
                                   )

    template = jinja_env.get_template(template_file)
    out = template.render(settings=settings)

    template_write_path = os.path.realpath(os.path.join(settings.experiment_dir,
                            "create_experiment","templates",template_file))
    fname = os.path.join(os.path.dirname(template_write_path), 
                template_file.replace(".jinja2",""))

    with open(fname, 'w') as f: f.write(out)

    exec_files = ["pism_prerun_script.sh","pism_run_script.sh"]
    if os.path.basename(fname) in exec_files:
        os.chmod(fname, os.stat(fname).st_mode | stat.S_IEXEC)

    print('   - created {} from template'.format(os.path.basename(fname)))


def get_pism_config_as_dict(settings):

    """ find the configs in standard settings.pism_config_file """

    pism_configs = {}

    for l in open(settings.pism_config_file,"r"):

        if "pism_config:" in l:
            key,val =  [s.strip() for s in l.split("=",1)]
            pism_configs[key.replace("pism_config:","")] = val.strip(";")

    return pism_configs


def check_if_override_is_in_config(settings,pism_config_dict):

    for k in settings.override_params:

        if k not in pism_config_dict.keys():
            raise ValueError(
                '{} is not in {}'.format(k, settings.pism_config_file))


def copy_from_template(settings, filename,
                       experiment=settings.experiment):

    """ just a dummy for testing. nothing is modified. """

    #settings.experiment_dir = os.path.join(settings.working_dir,
    #                              experiment)
    shutil.copy(os.path.join("templates",filename), settings.experiment_dir)
    print(os.path.join(settings.experiment_dir, filename), "copied.")


def create_run(settings=settings, experiment=settings.experiment):
    
    # copy template structure to new experiment location
    try:
        shutil.copytree(settings.coupl_template_dir, settings.experiment_dir, symlinks=True)
    except OSError as error:
        print(error)
        print("Choose a different experiment name or remove " + settings.experiment_dir)
        sys.exit(1)

    for d in ['x_MOM-to-PISM','x_PISM-to-MOM']:
        dir_path = os.path.join(settings.experiment_dir,d) 
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
    print(f" > created experiment directory {settings.experiment_dir}")

    # create main coupling script from template
    create_script_from_template(settings, "run_coupled.sh.jinja2")

    # prepare PISM subdirectory
    print(f" > set up PISM")
    PISM_folders = ['initdata', 'prerun', 'results']
    for f in PISM_folders:
        fpath = os.path.join(settings.pism_exp_dir, f)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
            print(f"   - created directory PISM/{f}")

    # copy PISM binary to experiment dir
    if not os.path.exists(settings.pism_exp_bin_dir):
        os.makedirs(settings.pism_exp_bin_dir)
    shutil.copy2(settings.pism_sys_bin, settings.pism_exp_bin_dir)
    print(f"   - copied PISM binary {settings.pism_sys_bin} to PISM/bin")

    # copy PISM input files to PISM/initdata/
    pism_input_files_to_copy = [
            settings.pism_infile_path, 
            settings.pism_atm_data_path, 
            settings.pism_ocn_data_path,
            settings.pism_ocnkill_data_path]
    if settings.pism_use_atm_anomaly_file:
        pism_input_files_to_copy.append(settings.pism_atm_anomaly_data_path)
    if settings.pism_use_atm_lapse_rate_file:
        pism_input_files_to_copy.append(settings.pism_atm_lapse_rate_data_path)
    for f in pism_input_files_to_copy:
        shutil.copy2(f, os.path.join(settings.pism_exp_dir, 'initdata'))
        print(f"   - copied PISM input file {f} to PISM/initdata")

    # adapt PISM input files
    if settings.pism_use_atm_anomaly_file:
        if hasattr(settings, 'pism_atm_anomaly_time_shift_years'):
            cmd = f"ncdump -h {settings.pism_atm_anomaly_data_path} | grep 'time:units'"
            time_units_str = str(subprocess.check_output(cmd, shell=True))
            if 'years' in time_units_str:
                time_shift = settings.pism_atm_anomaly_time_shift_years
            elif 'days' in time_units_str:
                time_shift = settings.pism_atm_anomaly_time_shift_years * 365
            elif 'seconds' in time_units_str:
                time_shift = settings.pism_atm_anomaly_time_shift_years * 365*24*3600
            else:
                raise ValueError(f"Cannot read time units of "
                    f"{settings.pism_atm_anomaly_data_path} to do time shift.")

            # shift time, time_bnds
            path_file = os.path.join(settings.pism_exp_dir, 'initdata',settings.pism_atm_anomaly_file)
            path_tmp = os.path.join(settings.pism_exp_dir, 'initdata',settings.pism_atm_anomaly_file + '.tmp')
            shutil.move(path_file, path_tmp)
            cmd = f"ncap2 -s 'time+={time_shift}ll; time_bnds+={time_shift}ll' {path_tmp} {path_file}"
            subprocess.call(cmd, shell=True)
            os.remove(path_tmp)
            print(f"   - applied time shift to PISM/initdata/{settings.pism_atm_anomaly_file} ({settings.pism_atm_anomaly_time_shift_years} years)")



    # prepare pism config_override file
    pism_config_dict = get_pism_config_as_dict(settings)
    check_if_override_is_in_config(settings, pism_config_dict)
    create_script_from_template(settings, "config_override.cdl.jinja2")
    cmd = "ncgen3 "+os.path.join(settings.pism_exp_dir,"initdata","config_override.cdl") \
            +" -o "+os.path.join(settings.pism_exp_dir,"initdata","config_override.nc")
    os.system(cmd)
    print("   - created PISM/config_override.nc from PISM/config_override.cdl")

    # create pism run scripts from template
    create_script_from_template(settings, "pism_prerun_script.sh.jinja2")
    create_script_from_template(settings, "pism_run_script.sh.jinja2")



    ## prepare POEM subdirectory
    print(f" > set up POEM")
    # delete all content first
    for filename in os.listdir(settings.poem_exp_dir):
        file_path = os.path.join(settings.poem_exp_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')
    # copy from template
    #shutil.copytree(settings.poem_template_dir, settings.poem_exp_dir, 
    #        dirs_exist_ok=True, symlinks=True)
    dist.copy_tree(settings.poem_template_dir, settings.poem_exp_dir, 
            preserve_symlinks=1, update=1, verbose=1)
    print(f"   - copied POEM template {settings.poem_template_dir} to POEM")

    # a new coupled run (not restarting from a previous coupled run)
    if settings.coupled_restart==False:
        # copy MOM restart files
        poem_input_dir = os.path.join(settings.poem_exp_dir,'INPUT')
        poem_restart_files_dir = settings.poem_restart_files_dir
        if os.path.exists(poem_restart_files_dir):
            with helpers.cd( str(poem_restart_files_dir) ):
                cmd = f"for i in *.res*; do rm {poem_input_dir}/$i; "\
                        f"cp -a $i {poem_input_dir} ; done"
                subprocess.call(cmd, shell=True)
            print(f"   - copied MOM restart files from {poem_restart_files_dir} "\
                    "to POEM/INPUT")
        else:
            warnings.warn(f"WARNING: path {poem_restart_files_dir} does not "\
                    f"exist! Need to copy MOM restart files to INPUT dir by "\
                    f"hand...")

    # set up data_table
    with helpers.cd( str(settings.poem_exp_dir) ):
        cmd = f"ln -sf {settings.poem_data_table_dummy} data_table-dummy"
        subprocess.call(cmd, shell=True)
    if settings.poem_data_table_replace is not {}:
        for pattern, replace_str in settings.poem_data_table_replace.items():
            #with helpers.cd( str(settings.poem_exp_dir) ):
            #    shutil.copy2(settings.poem_data_table_dummy, f"{settings.poem_data_table_dummy}.tmp")
            with helpers.cd( str(settings.poem_exp_dir) ):
                cmd = f'sed "s/{pattern}/{replace_str}/g" -i {settings.poem_data_table_dummy}'
                subprocess.call(cmd, shell=True)

    if settings.poem_copy_forcing_data:
        if not os.path.exists(settings.poem_forcing_data_target_path):
            os.mkdir(settings.poem_forcing_data_target_path)
        for f in glob.glob(settings.poem_forcing_data_source_path):
            shutil.copy2(f, settings.poem_forcing_data_target_path)
        #shutil.copy2(settings.poem_forcing_data_source_path,
        #             settings.poem_forcing_data_target_path)
        print(f"   - copied MOM forcing files "
              f"{settings.poem_forcing_data_source_path} to "
              f"{settings.poem_forcing_data_target_path}")

        # shift forcing data time
        if hasattr(settings, 'poem_forcing_time_shift_years'):
            file_paths = os.path.join(settings.poem_forcing_data_target_path,
                    settings.poem_forcing_data_source_pattern)
            cmd = f"ncdump -h {glob.glob(file_paths)[0]} | grep 'time:units'"
            time_units_str = str(subprocess.check_output(cmd, shell=True))
            if 'years' in time_units_str:
                time_shift = settings.poem_forcing_time_shift_years
            elif 'days' in time_units_str:
                time_shift = settings.poem_forcing_time_shift_years * 365
            elif 'seconds' in time_units_str:
                time_shift = settings.poem_forcing_time_shift_years * 365*24*3600
            else:
                raise ValueError(f"Cannot read time units of "
                    f"{glob.glob(file_paths)[0]} to do time shift.")

            # shift time, time_bnds
            for f in glob.glob(file_paths):
                file_ext = f.split('.')[-1]
                file_base = ".".join(f.split('.')[0:-1])
                f_shift = ".".join([file_base,'shift',file_ext])
                cmd = f"ncap2 -s 'time+={time_shift}ll; time_bnds+={time_shift}ll' {f} {f_shift}"
                subprocess.call(cmd, shell=True)
                print(f"   - applied time shift to {f} ({settings.poem_forcing_time_shift_years} years)")

    print(f" > set up additional coupling files")


    if (settings.do_ocean_tracer_anomaly==True and
        settings.use_ocean_tracer_anomaly_from_prev_run==True):
        # copy ocean tracer anomaly reference file from given path
        if os.path.exists(settings.ocean_tracer_anomaly_reference_path):
            shutil.copy2(settings.ocean_tracer_anomaly_reference_path,
                os.path.join(settings.experiment_dir, 'x_MOM-to-PISM'))
            print(f"   - copied ocean tracer anomaly reference file "\
                    f"{settings.ocean_tracer_anomaly_reference_path} "\
                    f"to compute ocean tracer anomalies "\
                    f"to same reference like in other run")
        else:
            warnings.warn(f"path {settings.ocean_tracer_anomaly_reference_path} "\
                f"does not exist!")

    if (settings.do_ocean_sealevel_anomaly==True and
        settings.use_ocean_sealevel_anomaly_from_prev_run==True):
        # copy ocean sealevel anomaly reference file from given path
        if os.path.exists(settings.ocean_sealevel_anomaly_reference_path):
            shutil.copy2(settings.ocean_sealevel_anomaly_reference_path,
                os.path.join(settings.experiment_dir, 'x_MOM-to-PISM'))
            print(f"   - copied ocean sealevel anomaly reference file "\
                    f"{settings.ocean_sealevel_anomaly_reference_path} "\
                    f"to compute ocean sealevel anomalies "\
                    f"to same reference like in other run")
        else:
            warnings.warn(f"path {settings.ocean_sealevel_anomaly_reference_path} "\
                f"does not exist!")

    # a new coupled run (not restarting from a previous coupled run)
    if settings.coupled_restart==False:
        # copy inital PISM-to-MOM fluxes for the first coupling iteration
        shutil.copy2(settings.pism_to_mom_flux_init_path,
                os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
        print(f"   - copied initial PISM-to-MOM flux file "\
                f"{settings.pism_to_mom_flux_init_file} from "\
                f"{settings.pism_to_mom_flux_init_path} for ice to ocean "\
                f"fluxes in first coupling iteration")

        # TODO: copy MOM RESTART files

    # restarting from a previous coupled run
    if settings.coupled_restart==True :
        # copy PISM-to-MOM fluxes file from previous run (last iteration)
        shutil.copy2(settings.pism_to_mom_flux_restart_path, 
                os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
        print(f"   - copied PISM-to-MOM flux file "\
                f"{settings.pism_to_mom_flux_restart_file} from "\
                f"{settings.restart_dir} to restart from previous run")

        if settings.use_prescribed_pico_input_depth==False:
            # copy PICO input depth restart file from previous run (last iteration)
            shutil.copy2(settings.pico_input_depth_restart_path, 
                    os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
            print(f"   - copied PICO input depth restart file "\
                    f"{settings.pico_input_depth_restart_file} from "\
                    f"{settings.restart_dir} to restart from previous run")

        if (settings.use_prescribed_basal_melt_input_depth==False and \
                settings.insert_basal_melt_at_depth==True):
            # copy basal melt input depth restart file from previous run (last iteration)
            shutil.copy2(settings.basal_melt_input_depth_restart_path, 
                    os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
            print(f"   - copied basal melt input depth restart file "\
                    f"{settings.basal_melt_input_depth_restart_file} from "\
                    f"{settings.restart_dir} to restart from previous run")


        # copy MOM restart files from previous run
        poem_input_dir = os.path.join(settings.poem_exp_dir,'INPUT')
        poem_restart_files_dir = os.path.join(settings.restart_dir,'POEM/INPUT')
        if os.path.exists(poem_restart_files_dir):
            with helpers.cd( str(poem_restart_files_dir) ):
                cmd = f"for i in *.res*; do rm {poem_input_dir}/$i; "\
                        f"cp -a $i {poem_input_dir} ; done"
                subprocess.call(cmd, shell=True)
            print(f"   - copied MOM restart files from {poem_restart_files_dir} "\
                    "to POEM/INPUT")
        else:
            warnings.warn(f"WARNING: path {poem_restart_files_dir} does not "\
                    f"exist! Need to copy MOM restart files to INPUT dir by "\
                    f"hand...")
    # [end] if settings.coupled_restart==True

    # copying PISM runoff reference file for calculation of ice to ocean runoff
    # with sea level impact
    if settings.do_runoff_slc == True:
        if settings.runoff_reference_surf_accum == True:
            if settings.coupled_restart == True:
                # copy ice-to-ocean runoff reference file from previous run
                if os.path.exists(settings.runoff_reference_restart_path):
                    shutil.copy2(settings.runoff_reference_restart_path,
                        os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
                    print(f"   - copied PISM to MOM runoff reference file "\
                            f"{settings.runoff_reference_restart_path} "\
                            f"from previous run (computed from PISM surface "\
                            f"accumulation flux) to identify the part of ice "\
                            f"to ocean runoff which changes sea level in the "\
                            f"ocean")
                else:
                    warnings.warn(f"WARNING: tried to copy PISM runoff "\
                            f"reference file, but path "\
                            f"{settings.runoff_reference_restart_path} "\
                            f"does not exist!")

        else:   # runoff_reference_surf_accum == False
            # copy pre-computed ice-to-ocean runoff reference file
            if os.path.exists(settings.runoff_reference_path):
                shutil.copy2(settings.runoff_reference_path,
                    os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
                print(f"   - copied pre-computed PISM runoff reference file "\
                        f"{settings.runoff_reference_path} to identify to "\
                        f"part of ice to ocean runoff which changes sea level "\
                        f"in the ocean")
            else:
                warnings.warn(f"WARNING: tried to copy PISM runoff reference file, "\
                        f"but path {settings.runoff_reference_path} does "\
                        f"not exist!")


    # copying prescribed PICO input depth (if used)
    if settings.use_prescribed_pico_input_depth==True:
        # copy prescribed PICO input depth file 
        shutil.copy2(settings.prescribed_pico_input_depth_path, 
                os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
        print(f"   - copied prescribed PICO input depth restart file "\
                f"{settings.pico_input_depth_restart_path}")

    # copying prescribed  basal melt depth (if used)
    if (settings.use_prescribed_basal_melt_input_depth==True and \
            settings.insert_basal_melt_at_depth==True):
        # copy prescribed basal melt input depth file 
        shutil.copy2(settings.prescribed_basal_melt_input_depth_path, 
                os.path.join(settings.experiment_dir, 'x_PISM-to-MOM'))
        print(f"   - copied prescribed basal melt input depth restart file "\
                f"{settings.basal_melt_input_depth_restart_path}")


if __name__ == "__main__":

    # redirect standard output to logfile
    orig_stdout = sys.stdout
    logfile = "create_run.log"
    print(f"writing script output to '{logfile}'")
    print("  -> check there for possible errors")
    with open(logfile, "w") as f:
        #sys.stdout = f
        #try:
        #    create_run()
        #finally:
        #    sys.stdout = orig_stdout
        with stdout_redirected(to=f), merged_stderr_stdout():
        # copied from https://stackoverflow.com/questions/6796492/temporarily-redirect-stdout-stderr
            create_run()


    # copy logfile to experiment location
    logfile_dst = os.path.join(settings.experiment_dir, logfile)
    shutil.copy2(logfile, logfile_dst)
    print(f"script ended sucessfully; logfile copied to {logfile_dst}")

