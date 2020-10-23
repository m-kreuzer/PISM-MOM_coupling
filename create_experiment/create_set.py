
import os
import sys
import shutil
import itertools
import collections
import hashlib
import pandas as pd
import create_run as cr
import settings

import code # debug

ensemble_parent_dir = ""

def create_set():

    iterables = settings.iterables
    param_iterables = settings.param_iterables
    
    # our ensemble table needs to contain keys and values of the iterables
    index_from_iterables = [[k]+[k+"_value"] for k in iterables.keys()]
    flat_list = [item for sublist in index_from_iterables for item in sublist]
    
    # collect things to vary, combine iterables and param_iterables
    values_to_vary = collections.OrderedDict()
    values_to_vary.update(param_iterables)
    
    for k in iterables.keys():
        values_to_vary[k] = iterables[k].keys()
    
    
    # find all combinations of iterables and param_iterables
    combinations = list(itertools.product(*values_to_vary.values()))
    labels = values_to_vary.keys()
    
    
    # create unique hash for each combination
    ensemble_table = pd.DataFrame(columns=list(param_iterables.keys())+flat_list)
    for i,pc in enumerate(combinations):
        em_id = hashlib.sha224(str(pc).encode('utf-8')).hexdigest()[0:8]
        for j,l in enumerate(labels):
            ensemble_table.loc[em_id,l] = list(pc)[j]
    
    # split up topg_to_phi to original parameters.
    if "topg_to_phi" in ensemble_table.columns:
        for em in ensemble_table.index:
            tphip = ensemble_table.loc[em,"topg_to_phi"]
            ensemble_table.loc[em,
            "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"] = tphip[0]
            ensemble_table.loc[em,
            "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"] = tphip[1]
            ensemble_table.loc[em,
            "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"] = tphip[2]
            ensemble_table.loc[em,
            "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max"] = tphip[3]
    
        ensemble_table = ensemble_table.drop("topg_to_phi",axis=1)
    
    # also fill values from iterables dict
    for k in iterables.keys():
        for i,v in enumerate(iterables[k].keys()):
            ensemble_table.loc[ensemble_table[k] == v,k+"_value"] = iterables[k][v]
    
    
    # create ensemble parent directory
    global ensemble_parent_dir
    ensemble_parent_dir = os.path.join(settings.working_dir, settings.experiment + "_ensemble")
    
    if not os.path.exists(ensemble_parent_dir):
        os.makedirs(ensemble_parent_dir)
    else:
        print("Directory {} exists. Choose a different experiment name or remove directory".format(ensemble_parent_dir))
        sys.exit(1)

    
    # bring ensemble table entries to settings and write experiments
    for ind in ensemble_table.index:
    
        for col in ensemble_table.columns:
            if col in iterables:
                # case of iterables, e.g., oceanfiles, should be written to pism_run file
                settings.__dict__[col] = ensemble_table.loc[ind,col+"_value"]
            elif any(iterable in col for iterable in iterables):
                # case of oceanfile_value which should be simply ignored, should not be written to config_override
                pass 
            else:
                # case of standard parameters, should be written to config_override
                settings.override_params[col] = ensemble_table.loc[ind,col]
    
        # adapt experiment directory in settings and update all derived path definitions
        experiment = settings.experiment+"_"+ind
        new_experiment_dir = os.path.join(ensemble_parent_dir,experiment)
        settings.__dict__['experiment_dir'] = new_experiment_dir
        settings.__dict__['poem_exp_dir'] = os.path.join(new_experiment_dir, 'POEM')
        settings.__dict__['pism_exp_dir'] = os.path.join(new_experiment_dir, 'PISM')
        settings.__dict__['pism_exp_bin_dir'] = os.path.join(new_experiment_dir, 'PISM', 'bin')
        settings.__dict__['pism_exp_bin'] = os.path.join(new_experiment_dir, 'PISM', 'bin', settings.pism_exec)
        
        # build run directory 
        cr.create_run(settings=settings, experiment=experiment)


    table_path = os.path.join(ensemble_parent_dir,"ensemble_hash_table.txt")
    ensemble_table.to_csv( table_path, sep=" ", index_label="hash")
    print("Wrote ensemble table to "+ table_path)

if __name__ == "__main__":

    # redirect standard output to logfile
    orig_stdout = sys.stdout
    logfile = "create_set.log"
    print("writing script output to '{}'".format(logfile))
    print("  -> check there for possible errors")
    with open(logfile, "w") as f:
        sys.stdout = f
        try:
            create_set()
        finally:
            sys.stdout = orig_stdout

    # copy logfile to experiment location
    shutil.copy2(logfile, ensemble_parent_dir)
    print("script ended sucessfully; logfile copied to {}".format(ensemble_parent_dir))


