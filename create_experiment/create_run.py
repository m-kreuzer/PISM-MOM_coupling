
import os
import stat
import shutil
import jinja2
import collections
import settings


def write_pism_script(settings, template_file,
                      experiment=settings.experiment,
                      printe=False):

    experiment_dir = os.path.join(settings.pism_experiments_dir,
                                  experiment)

    if not os.path.exists(experiment_dir):
        os.makedirs(experiment_dir)

    # make jinja aware of templates
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(
        searchpath=os.path.join(settings.project_root,"templates")))

    template = jinja_env.get_template(template_file)
    out = template.render(settings=settings)

    fname = os.path.join(experiment_dir,
                         template_file.replace(".jinja2",""))

    with open(fname, 'w') as f: f.write(out)

    if template_file == "pism_run.sh.jinja2":
        os.chmod(fname, os.stat(fname).st_mode | stat.S_IEXEC)

    if printe:
        print("##", experiment_dir)

    print(fname, "written.")


def get_pism_config_as_dict(settings):

    """ find the configs in standard settings.pism_config_file """

    pism_configs = {}

    for l in open(settings.pism_config_file,"r"):

        if "pism_config:" in l and "_doc" not in l:
            key,val =  [s.strip() for s in l.split("=")]
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

    experiment_dir = os.path.join(settings.pism_experiments_dir,
                                  experiment)
    shutil.copy(os.path.join("templates",filename), experiment_dir)
    print(os.path.join(experiment_dir, filename), "copied.")


if __name__ == "__main__":

    write_pism_script(settings, "pism_run.sh.jinja2",printe=True)
    write_pism_script(settings, "submit.sh.jinja2")
    write_pism_script(settings, "prepare_restart.sh.jinja2")

    pism_config_dict = get_pism_config_as_dict(settings)
    check_if_override_is_in_config(settings, pism_config_dict)
    write_pism_script(settings, "config_override.cdl.jinja2")
