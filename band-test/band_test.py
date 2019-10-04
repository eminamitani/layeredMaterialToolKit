import sys
import pathlib
current_dir = pathlib.Path(__file__).resolve().parent
sys.path.append( str(current_dir) + '/../' )
from espresso import *
import os

input_data = "input_data.txt"
pseudo_dir = "../pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.create_band_input(scf_file='MoS2-optimized.pwi', band_config='band_config.txt')