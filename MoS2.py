from ase.build import mx2
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
from espresso import *
import os

input_data = "input_data.txt"
pseudo_dir = "./pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.optimize_lattice_constant(initial_lattice_constant=3.3,range=5,samples=20,kpts=(12,12,1),vacuum=25.0)

