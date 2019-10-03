from ase.build import mx2
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
import sys
import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../' )
from espresso import *
import os

input_data = "input_data.txt"
pseudo_dir = "./pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.gather_optimize_results('log_lattice_optim')


