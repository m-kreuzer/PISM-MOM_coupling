###
# This script is designed to continue a run from a number of promising runs at given years. 
# There are hand-selected for example using pism-ensemble-analysis and written to a file
# as indicated in settings.py 
# The file should contain the hash and year of the runs to be continued
#

import os
import itertools
import collections
import hashlib
import pandas as pd
import create_run as cr
import settings

# iterables = settings.iterables
# param_iterables = settings.param_iterables

# # our ensemble table needs to contain keys and values of the iterables
# index_from_iterables = [[k]+[k+"_value"] for k in iterables.keys()]
# flat_list = [item for sublist in index_from_iterables for item in sublist]

# # collect things to vary, combine iterables and param_iterables
# values_to_vary = collections.OrderedDict()
# values_to_vary.update(param_iterables)

# for k in iterables.keys():
#     values_to_vary[k] = iterables[k].keys()

# # find all combinations of iterables and param_iterables
# combinations = list(itertools.product(*values_to_vary.values()))
# labels = values_to_vary.keys()

print('Full ensemlbe: '+os.path.join("sets",settings.source_ensemble_table))

ensemble_table = pd.read_csv(os.path.join("sets",settings.source_ensemble_table),
                             sep='\s+',index_col=0)

print('Runs to be continued are defined in:'+settings.runs_to_continue)
runs_to_continue = pd.read_csv(settings.runs_to_continue,
                             sep='\s+')


# ensemble_table = pd.DataFrame(columns=list(param_iterables.keys())+flat_list)

# for i,pc in enumerate(combinations):
#     em_id = hashlib.sha224(str(pc).encode('utf-8')).hexdigest()[0:8]
#     for j,l in enumerate(labels):
#         ensemble_table.loc[em_id,l] = list(pc)[j]

# split up topg_to_phi to original parameters.
# if "topg_to_phi" in ensemble_table.columns:
#     for em in ensemble_table.index:
#         tphip = ensemble_table.loc[em,"topg_to_phi"]
#         ensemble_table.loc[em,
#         "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"] = tphip[0]
#         ensemble_table.loc[em,
#         "basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"] = tphip[1]
#         ensemble_table.loc[em,
#         "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"] = tphip[2]
#         ensemble_table.loc[em,
#         "basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max"] = tphip[3]

#     ensemble_table = ensemble_table.drop("topg_to_phi",axis=1)

# also fill values from iterables dict
# for k in iterables.keys():
#     for i,v in enumerate(iterables[k].keys()):
#         ensemble_table.loc[ensemble_table[k] == v,k+"_value"] = iterables[k][v]

# bring ensemble table entries to settings and write experiments

# if settings.hashes_to_run == "all":
#     hashes = ensemble_table.index
# else:
#     hashes = hashes_to_run

for ehash, year in runs_to_continue.values:
    
    for col in ensemble_table.columns:

        # TODO: ensure that this option works.
        if col in settings.iterables:
            settings.__dict__[col] = ensemble_table.loc[ehash,col+"_value"]
        else:
            settings.override_params[col] = ensemble_table.loc[ehash,col]

    experiment = settings.experiment+"_"+ehash+"_"+str(year)

    # TODO make this more general
    settings.infile = settings.get_infile_to_continue(ehash, year)
    # print(settings.infile)

    cr.write_pism_script(settings, "pism_run.sh.jinja2",
                      experiment=experiment)
    cr.write_pism_script(settings, "submit.sh.jinja2",
                      experiment=experiment)
    cr.write_pism_script(settings, "config_override.cdl.jinja2",
                          experiment=experiment)
    cr.write_pism_script(settings, "prepare_restart.sh.jinja2",
                          experiment=experiment)


ensemble_table.to_csv(os.path.join("sets",settings.experiment+".txt"),
                      sep=" ", index_label="hash")
print("Wrote ensemble table to", os.path.join("sets",settings.experiment+".txt"))

