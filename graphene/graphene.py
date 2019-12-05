from layeredMaterialToolKit.espresso import *

input_data = "input_data.txt"
pseudo_dir = "./pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.generate_different_lattice_constant_inputs(initial_lattice_constant=3.1, range=5, samples=2, kpts=(16, 16, 1), vacuum=15.0)