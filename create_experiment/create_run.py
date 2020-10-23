
import sys
import os
import stat
import shutil
import jinja2
import collections
import distutils.dir_util as dist

import settings 



def create_script_from_template(settings, template_file,
                      experiment=settings.experiment):

    if not os.path.exists(settings.experiment_dir):
        os.makedirs(settings.experiment_dir)

    # make jinja aware of templates
    template_path = os.path.join(settings.project_root,"templates")
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=template_path))

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
    print(" > created experiment directory " + settings.experiment_dir)

    # create main coupling script from template
    create_script_from_template(settings, "run_coupled.sh.jinja2")

    # prepare PISM subdirectory
    PISM_folders = ['initdata', 'prerun', 'results']
    for f in PISM_folders:
        fpath = os.path.join(settings.pism_exp_dir, f)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
            print("   - created directory PISM/"+f)

    # copy PISM binary to experiment dir
    if not os.path.exists(settings.pism_exp_bin_dir):
        os.makedirs(settings.pism_exp_bin_dir)
    shutil.copy2(settings.pism_sys_bin, settings.pism_exp_bin_dir)
    print("   - copied PISM binary {} to PISM/bin".format(settings.pism_sys_bin))

    # copy PISM input files to PISM/initdata/
    pism_input_files_to_copy = [
            settings.pism_infile_path, 
            settings.pism_atm_data_path, 
            settings.pism_ocn_data_path,
            settings.pism_ocnkill_data_path]
    for f in pism_input_files_to_copy:
        shutil.copy2(f, os.path.join(settings.pism_exp_dir, 'initdata'))
        print("   - copied PISM input file {} to PISM/initdata".format(f))

    #shutil.copy2(settings.pism_infile_path, os.path.join(settings.pism_exp_dir, 'initdata'))
    #print("   - copied PISM input file {} to PISM/initdata".format(settings.pism_infile_path))
    #shutil.copy2(settings.pism_atm_data_path, os.path.join(settings.pism_exp_dir, 'initdata'))
    #print("   - copied PISM input file {} to PISM/initdata".format(settings.pism_atm_data_path))
    #shutil.copy2(settings.pism_ocn_data_path, os.path.join(settings.pism_exp_dir, 'initdata'))
    #print("   - copied PISM input file {} to PISM/initdata".format(settings.pism_ocn_data_path))
    #shutil.copy2(settings.pism_ocnkill_data_path, os.path.join(settings.pism_exp_dir, 'initdata'))
    #print("   - copied PISM input file {} to PISM/initdata".format(settings.pism_ocnkill_data_path))


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
    # delete all content first
    for filename in os.listdir(settings.poem_exp_dir):
        file_path = os.path.join(settings.poem_exp_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
    # copy from template
    #shutil.copytree(settings.poem_template_dir, settings.poem_exp_dir, dirs_exist_ok=True, symlinks=True)
    dist.copy_tree(settings.poem_template_dir, settings.poem_exp_dir, preserve_symlinks=1, update=1, verbose=1)
    print("   - copied POEM template {} to POEM".format(settings.poem_template_dir))



if __name__ == "__main__":

    # redirect standard output to logfile
    orig_stdout = sys.stdout
    logfile = "create_run.log"
    print("writing script output to '{}'".format(logfile))
    print("  -> check there for possible errors")
    with open(logfile, "w") as f:
        sys.stdout = f
        try:
            create_run()
        finally:
            sys.stdout = orig_stdout

    # copy logfile to experiment location
    logfile_dst = os.path.join(settings.experiment_dir, logfile)
    shutil.copy2(logfile, logfile_dst)
    print("script ended sucessfully; logfile copied to {}".format(logfile_dst))

