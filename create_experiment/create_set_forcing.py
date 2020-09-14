
import os
import itertools
import collections
import hashlib
import pandas as pd
import create_run as cr
import settings


print('Full ensemlbe: '+os.path.join("sets",settings.source_ensemble_table))

ensemble_table = pd.read_csv(os.path.join("sets",settings.source_ensemble_table),
                             sep='\s+',index_col=0)

print('Runs to be continued are defined in:'+settings.runs_for_forcing)
runs_for_forcing = pd.read_csv(settings.runs_for_forcing,
                             sep='\s+')

print(runs_for_forcing.values)

for ehash, year in runs_for_forcing.values: # FIXME this works currently only for 1 entry!
	## get parameters used in the run for forcing and replace in settings
	#for col in ensemble_table.columns:
	#	# TODO: ensure that this option works.
	#	if col in settings.iterables:
	#		settings.__dict__[col] = ensemble_table.loc[ehash,col+"_value"]
	#	elif any(iterable in col for iterable in iterables):
	#		# case of oceanfile_value which should be simply ignored, should not be written to config_override
	#		pass 
	#	else:
	#		settings.override_params[col] = ensemble_table.loc[ehash,col]

	settings.infile_forcing = settings.get_infile_to_continue(str(ehash), year)
	print(settings.infile_forcing)
	# set newly defined iterables
	# !Attention, model parameters might be overwritten
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

	# bring ensemble table entries to settings and write experiments
	for ind in ensemble_table.index:
		for col in ensemble_table.columns:
			if col in iterables:
				settings.__dict__[col] = ensemble_table.loc[ind,col+"_value"]
			elif any(iterable in col for iterable in iterables):
				# case of oceanfile_value which should be simply ignored, should not be written to config_override
				pass 
			else:
				settings.override_params[col] = ensemble_table.loc[ind,col]

		experiment = settings.experiment+"_"+ind
	
		# settings.infile = os.path.join(settings.iinfile_forcing+ind,"snapshots_2300.000.nc")
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

