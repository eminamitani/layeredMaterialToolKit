from layeredMaterialToolKit.espresso import *

input_data = "input_data.txt"
pseudo_dir = "../pseudo"
config = "config.txt"
PBS_header="PBS_header.txt"

esp = espresso(input_data=input_data, pseudo_dir=pseudo_dir, config=config, PBS_header=PBS_header)
esp.gather_phonon_info(scf_file='MoS2-scf.pwi', phonon_config_file='phonon_config.txt')