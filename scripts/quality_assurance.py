import yaml
import netCDF4
import numpy as np
import xarray as xr
from functions import log, advanced_quality_flags
import glob
from adcp import ADCP_data

def run_quality_assurance(files, p):
    first = True
    for file in files: 
        dataset = xr.open_dataset(file, decode_times =False)
        if first:
            combined_dataset = dataset
            first = False
        else:
            if max(np.abs(combined_dataset["depth"].values-dataset["depth"].values))>1:
                log("RDI moved from more than 1m, depth grid can no longer be trusted")
            dataset["depth"] = combined_dataset["depth"]
            combined_dataset = xr.combine_by_coords([combined_dataset, dataset], combine_attrs='override')

    log("Apply advance quality checks to data")
    advanced_dataset = advanced_quality_flags(combined_dataset, json_path="quality_assurance.json")
    dataset.close()

    log("Update NetCDF files with advanced QA")
    for file in files:
        dset = netCDF4.Dataset(file, 'r+')
        idx = np.where((advanced_dataset["time"] >= dset["time"][0]) & (advanced_dataset["time"] <= dset["time"][-1]))[0]
        for var in dset.variables:
            if "_qual" in var:
                if advanced_dataset[var][:].ndim==1:
                    dset[var][:][np.array(advanced_dataset[var][:][idx], dtype=bool)] = 1
                else:
                    if advanced_dataset[var][:].ndim==2:
                        dset[var][:][np.array(advanced_dataset[var][:][:,idx], dtype=bool)]
                    else: 
                        log("didn't recognized number of dimensions")
        dset.close()
    log("Updating L2 data with advanced QA")
    for file in files:
        ADCP = ADCP_data()
        ADCP.read_nc(file)
        ADCP.mask_data()
        ADCP.derive_variables(p["rotate_velocity"])
        new_folder = file[:file.rfind("L1")]
        ADCP.to_netcdf(new_folder.replace("Level1", "Level2"), "L2")

log("Performing advanced quality check")
with open("scripts/input_python.yaml", "r") as f:
    directories = yaml.load(f, Loader=yaml.FullLoader)

files600 = glob.glob(directories["Level1_dir"]+"/RDI600/"+"*.nc")
files300 = glob.glob(directories["Level1_dir"]+"/RDI300/"+"*.nc")
files600.sort()
files300.sort()

p300 = {"bandwidth": "300", "transducer_depth": 8, "bottom_depth": 110, "up": "False", "irt": 6*24, "rotate_velocity": 0, "cabled": True}
p600 = {"bandwidth": "600", "transducer_depth": 8, "bottom_depth": 110, "up": "True", "irt": 6*24, "rotate_velocity": 0, "cabled": True}


run_quality_assurance(files600, p600)
run_quality_assurance(files300, p300)