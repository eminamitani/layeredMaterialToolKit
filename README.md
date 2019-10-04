# layeredMaterialToolKit
This is a tiny scripting tool to support annoying preparation 
of DFT calculations especially in layered materials.

## requirement
ASE > 3.18.1, 
numpy

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
optimize_lattice_constant: create input files to manual optimization for lattice constant.
(to avoid trap local minimum in vc-relax)
This function creates several directories to run pw.x calculation and store the results.
`run.sh` is also created to submit the job in cluster machines.
The headers of job schedulers and mpi commands are specified through `PBS_header` and `config` files.