# -*- coding: utf-8 -*-
import os
import sys
import json
import yaml
from adcp import ADCP
from functions import latest_files, all_files, log, select_parameters


with open("scripts/input_python.yaml", "r") as f:
    directories = yaml.load(f, Loader=yaml.FullLoader) # Path to the data directories L0, L1, L2
    
for directory in directories.values():
    if not os.path.exists(directory):
        os.makedirs(directory) # Create the directories if they do not exist

with open('scripts/parameters.json', 'r') as f:
    parameter_dict = json.load(f) # Parameters related to the settings of the ADCP for each deployment period 

if len(sys.argv) == 1: # Reprocess all the L0 files
    files_300 = all_files(directories["Level0_dir"], "RDI300")
    files_600 = all_files(directories["Level0_dir"], "RDI600")
    log("Reprocessing complete dataset.")

elif len(sys.argv) == 2 and str(sys.argv[1]) == "live": # Live data: process recent L0 files only
    files_300 = latest_files(directories["Level0_dir"], "RDI300")
    files_600 = latest_files(directories["Level0_dir"], "RDI600")
    log("Live processing recent data.")
else:
    raise ValueError()

for file in files_300:
    a = ADCP()
    p = select_parameters(file, parameter_dict) # Get the parameters corresponding to the specific deployment period
    if a.read_data(file, transducer_depth=p["transducer_depth"], bottom_depth=p["bottom_depth"], up=p["up"]): # Read the raw data  
        a.quality_flags("./quality_assurance.json") # Flag the data based on quality checks
        a.export(os.path.join(directories["Level1_dir"], "RDI300"), "L1", output_period="file", overwrite_file=True) # Create Level 1 file
        a.mask_data() # Replace flagged data by nan 
        a.derive_variables(p["rotate_velocity"]) # Compute additional variables to add to Level 2
        a.export(os.path.join(directories["Level2_dir"], "RDI300"), "L2", output_period="file", overwrite_file=True) # Create Level 2 file

for file in files_600:
    a = ADCP()
    p = select_parameters(file, parameter_dict)
    if a.read_data(file, transducer_depth=p["transducer_depth"], bottom_depth=p["bottom_depth"], up=p["up"]):
        a.quality_flags("./quality_assurance.json")
        a.export(os.path.join(directories["Level1_dir"], "RDI600"), "L1", output_period="file", overwrite_file=True)
        a.mask_data()
        a.derive_variables(p["rotate_velocity"])
        a.export(os.path.join(directories["Level2_dir"], "RDI600"), "L2", output_period="file", overwrite_file=True)
