## pism-run - create one or many pism run scripts

pism-run allows you to create the bash scripts that define your pism simulation.
You can define ensemble parameters and create a set of scripts where only
these parameters vary.
As this is git tracked, you can easily see the changes you did to your simulations.
Separate settings for different clusters allow easy switch between machines
without code edits.

### Usage

Edit the `settings.py`, `supermuc_settings.py` and `pikcluster_settings.py` files.

Take special care of naming your experiment through `experiment` in settings.py.
You may overwrite older run scripts if this is not changed.

`python create_run.py` will create a single run.

`python create_set.py` will create a ensemble of runs.
Currently, `iterables` and `param_iterables` from `settings.py` are put together in all possible ways through python's `itertools.product`. A file in the subfolder `sets/`
with the experiment name will be created, which holds the information of each ensemble member.

`python submit_set.py` will submit the set of runs you created before.

You may read the ensemble file in python with

 `pd.read_csv("your_ens_file.txt", sep='\s+", index_col=0)`

Hashes identify your run and relate the set of parameters to the
script folder for running pism. The hash is part of each folder name of
an ensemble.

Not all changes to the run scripts can be done via `settings.py`. You may
edit the templates in `templates/` for heavier tweaks.

### License

pism-run is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License, or (at your option) any later version.

pism-run is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You find a copy of the GNU General Public License along with pism-run in the file LICENSE.txt;
