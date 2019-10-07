# layeredMaterialToolKit
This is a tiny scripting tool to support annoying preparation 
of DFT calculations especially in layered materials.

## requirement
ASE > 3.18.1, 
numpy,
spglib

## install
This program is still constructing status.
Please use virtualenv & be careful to install.  
```
git clone https://github.com/eminamitani/layeredMaterialToolKit.git
cd cd layeredMaterialToolKit/
pip install -e .
```

## required files
Normaly, three files are required to setup calculation condition.
The name of files are not specified, here I add some explanation on sample files including this repository.


input_data.txt : text file for dictionary of Quantum-espresso calculation

config.txt: text file for dictionary of the other set-ups

PBS_header.txt: text file contains the header section of jobscript.  

## memo
In the constructor of espresso class, the information of pseudopotetials in target directory is parsed automatically and stored in `self.pseudopotentials`.

## example
see MoS2-sample for setting the input_data, config, PBS_header etc.

## usage
`generate_different_lattice_constant_inputs` module to create input files to manual optimization for lattice constant.
(to avoid trap local minimum in vc-relax)
This function creates several directories to run pw.x calculation and store the results.
`run.sh` is also created to submit the job in cluster machines.
The headers of job schedulers and mpi commands are specified through `PBS_header` and `config` files.

`gather_optimize_results` module to gather the calculation results in the directories generated by `generate_different_lattice_constant_inputs`.
Used to determine most stable structures & generate scf calculation input with the optimized structure.


`create_band_input` module to generate band structure input for `pw.x` & `bands.x`.
See `band-test/band_test.py`

`create_phonon_input` module to generate phonon calculation input. To reduce the calculation time, here I assume that phonon calculations are split into several jobs.
The number of jobs are determined by `njobs` in phonon_config.
This module generates `njobs` directories and stores input files for each directory.
This module also generates `runall.sh` script to submit all jobs.
See `phonon_test/phonon_test.py`.

`gather_phonon_info` module to gather the calculation results stored in several directories generated by 
`create_phonon_input`. In order to calculate phonon band structure, this module 
also generate input files for `q2r.x` &
  `matdyn.x`. The number of q-point in phonon band calculation are given by `'nbandkpts'` in phonon_condig.
 If you set `'calcRun':True` in phonon_config file, `q2r.x` & `matdyn.x` call directory from this python script.
 If you set `'wfClean':True` , the wavefuntion files are removed to save your storage.
 See `phonon_test/phonon_post.py`.

## example of phonon_config
```
{
'nscfkpts':[24,24,1],
'nphmesh':[6,6,1],
'njobs':4,
'fildvscf':'dvscf',
'tr2_ph':1.0e-16,
'ldisp':True,
'epsil':True,
'trans':True,
'search_sym':False,
'nbandkpts':100,
'wfClean':True,
'calcRun':False,
}
```
