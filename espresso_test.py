
from espresso import espresso

input_data = "input_data.txt"
pseudo_dir = "./pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
