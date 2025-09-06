# -*- coding: utf-8 -*-
import os
import sys
import yaml
import json
import time
import argparse
import requests
from instruments import ADCP
from general.functions import logger, files_in_directory
from functions import retrieve_new_files, select_parameters

def main(server=False, logs=False):
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if logs:
        log = logger(os.path.join(repo, "logs/adcp"))
    else:
        log = logger()
    log.initialise("Processing LéXPLORE ADCP data")
    directories = {f: os.path.join(repo, "data", f) for f in ["Level0", "Level1", "Level2"]}
    for directory in directories:
        os.makedirs(directories[directory], exist_ok=True)
    edited_files = []

    log.begin_stage("Collecting mooring parameters")
    with open(os.path.join(repo, 'notes/parameters.json'), 'r') as f:
        parameter_dict = json.load(f)
    log.end_stage()

    log.begin_stage("Collecting inputs")
    if server:
        log.info("Processing files from sftp server")
        if not os.path.exists(os.path.join(repo, "creds.json")):
            raise ValueError("Credential file required to retrieve live data from the fstp server.")
        with open(os.path.join(repo, "creds.json"), 'r') as f:
            creds = json.load(f)
        files = retrieve_new_files(directories["Level0"], creds, server_location=["data/ADCP_300", "data/ADCP_600"], filetype=".LTA")
        edited_files = edited_files + files
    else:
        files = files_in_directory(directories["Level0"])
        files.sort()
        log.info("Reprocessing complete dataset from {}".format(directories["Level0"]))
    log.end_stage()

    files = [files[10]]

    log.begin_stage("Processing data")
    for file in files:
        sensor = ADCP(log=log)
        p = select_parameters(file, parameter_dict)
        if sensor.read_data(file, transducer_depth=p["transducer_depth"], bottom_depth=p["bottom_depth"], cabled=p["cabled"], up=p["up"]):
            sensor.quality_flags(envass_file=os.path.join(repo, "notes/quality_assurance.json"), adcp_file=os.path.join(repo, 'notes/quality_specific_adcp.json'))
            edited_files.extend(sensor.export(os.path.join(directories["Level1"], "RDI" + p["bandwidth"]), "L1_ADCP_", output_period="weekly", overwrite=True))
            sensor.mask_data()
            sensor.derive_variables(p["rotate_velocity"])
            edited_files.extend(sensor.export(os.path.join(directories["Level2"], "RDI" + p["bandwidth"]), "L2_ADCP_", output_period="weekly", overwrite=True))
    log.end_stage()

    return edited_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--server', '-s', help="Collect and process new files from FTP server", action='store_true')
    parser.add_argument('--logs', '-l', help="Write logs to file", action='store_true')
    args = vars(parser.parse_args())
    main(server=args["server"], logs=args["logs"])