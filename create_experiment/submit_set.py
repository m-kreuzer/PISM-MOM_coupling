"""
Submit the ensemble of PISM runs, as created with
create_set.py, to the batch system of cluster PIK or SUPERMUC.
"""

import os
import pandas as pd
import subprocess
import settings

# Set this to True if you want to continue previous runs, when submitting this the first time, set to False
continued_submission = False

def submit(ens_member_name, ens_path):

    ens_member_path = os.path.join(ens_path,ens_member_name)
    print(ens_member_path)
    subprocess.check_call("cd "+ens_member_path+" && "+settings.submit_command,
                         shell=True)

if continued_submission:

    runs_to_continue = pd.read_csv(settings.runs_to_continue,
                             sep='\s+')

    for ehash, year in runs_to_continue.values:
        ens_member_name = settings.experiment+"_"+ehash+"_"+str(year)
        print(ens_member_name)
        submit(ens_member_name)


if not continued_submission:

    ensemble_parent_dir = os.path.join(settings.working_dir, settings.experiment + "_ensemble")
    ensemble_hash_table_path = os.path.join(ensemble_parent_dir,"ensemble_hash_table.txt")
    ensemble_table = pd.read_csv(ensemble_hash_table_path, sep='\s+',index_col=0)

    for ehash in ensemble_table.index[:]:
        ens_member_name = settings.experiment+"_"+ehash
        submit(ens_member_name=ens_member_name, ens_path=ensemble_parent_dir)

