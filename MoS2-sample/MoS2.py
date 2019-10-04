import sys
import pathlib
#current_dir = pathlib.Path(__file__).resolve().parent
#sys.path.append( str(current_dir) + '/../' )
from layeredMaterialToolKit.espresso import *

input_data = "input_data.txt"
pseudo_dir = "./pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.optimize_lattice_constant(initial_lattice_constant=3.1,range=5,samples=20,kpts=(16,16,1),vacuum=25.0)

