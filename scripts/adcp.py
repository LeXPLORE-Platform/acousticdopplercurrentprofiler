# -*- coding: utf-8 -*-
import os
import netCDF4
import numpy as np
import xarray as xr
import dolfyn as dlfn
import distutils.util
from functions import *
from datetime import datetime
from envass import qualityassurance
from dateutil.relativedelta import relativedelta
from general_functions import GenericInstrument
from quality_checks_adcp import *


class ADCP(GenericInstrument):
    def __init__(self, *args, **kwargs):
        super(ADCP, self).__init__(*args, **kwargs)
        self.general_attributes = {
            "institution": "EPFL",
            "source": "ADCP",
            "references": "LéXPLORE commun instruments camille.minaudo@epfl.ch>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Data from ADCP on Lexplore Platform in Lake Geneva",
            "title": "LéXPLORE ADCP Velocities",
            "Instrument": "ADCP",
        }
        
        self.dimensions = {
            'time': {'dim_name': 'time', 'dim_size': None},
            'depth': {'dim_name': 'depth', 'dim_size': None}
        }
        
        self.variables = {
            'time': {'var_name': 'time', 'dim': ('time',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'time'},
            'depth': {'var_name': 'depth', 'dim': ('depth',), 'unit': 'm', 'long_name': 'nominal depth'},
            'u': {'var_name': 'u', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'long_name': 'eastern velocity'},
            'v': {'var_name': 'v', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'long_name': 'northern velocty'},
            'w': {'var_name': 'w', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'long_name': 'upward velocty'},
            'temp': {'var_name': 'temp', 'dim': ('time',), 'unit': 'degC', 'long_name': 'temperature', },
            'eu': {'var_name': 'vel_err', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'long_name': 'velocity error', },
            'echo1': {'var_name': 'echo1', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 1 echo'},
            'echo2': {'var_name': 'echo2', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 2 echo'},
            'echo3': {'var_name': 'echo3', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 3 echo'},
            'echo4': {'var_name': 'echo4', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 4 echo'},
            'corr1': {'var_name': 'corr1', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 1 correlation', },
            'corr2': {'var_name': 'corr2', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 2 correlation', },
            'corr3': {'var_name': 'corr3', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 3 correlation', },
            'corr4': {'var_name': 'corr4', 'dim': ('depth', 'time'), 'unit': '-', 'long_name': 'Beam 4 correlation', },
            'prcnt_gd1': {'var_name': 'prcnt_gd1', 'dim': ('depth', 'time'), 'unit': '%', 'long_name': 'Percentage Good 1', },
            'prcnt_gd2': {'var_name': 'prcnt_gd2', 'dim': ('depth', 'time'), 'unit': '%', 'long_name': 'Percentage Good 2', },
            'prcnt_gd3': {'var_name': 'prcnt_gd3', 'dim': ('depth', 'time'), 'unit': '%', 'long_name': 'Percentage Good 3', },
            'prcnt_gd4': {'var_name': 'prcnt_gd4', 'dim': ('depth', 'time'), 'unit': '%', 'long_name': 'Percentage Good 4', },
            'battery': {'var_name': 'battery', 'dim': ('time',), 'unit': '-', 'long_name': 'Battery level'},
            'heading': {'var_name': 'heading', 'dim': ('time',), 'unit': 'deg', 'long_name': 'Heading'},
            'roll': {'var_name': 'roll', 'dim': ('time',), 'unit': 'deg', 'long_name': 'Roll'},
            'pitch': {'var_name': 'pitch', 'dim': ('time',), 'unit': 'deg', 'long_name': 'Pitch'},
        }

        self.derived_variables = {
            'mU': {'var_name': 'mU', 'dim': ('time',), 'unit': 'm s-1', 'long_name': 'modulus of depth-averaged velocity'},
            'mdir': {'var_name': 'mdir', 'dim': ('time',), 'unit': 'deg', 'long_name': 'direction (anticlockwise from east) of depth-averaged velocity'},
            'Sv': {'var_name': 'Sv', 'dim': ('depth', 'time'), 'unit': 'dB', 'long_name': 'absolute backscatter'}
        }

        self.data = {}

    def read_data(self, file, transducer_depth, bottom_depth=110., cabled=False, up=False, **kwargs):
        """
        Read the ADCP data and store it in an ADCP object.

        Parameters:
            file (str): path and ADCP filename (e.g., .LTA file)
            transducer_depth (float): depth of the ADCP [m]
            bottom_depth (float): total lake depth [m] (default: lake depth at the LéXPLORE platform)
            cabled (bool): =True is the ADCP is cable-linked, =False otherwise
            up (bool): =True is the ADCP is upward-looking, =False if the ADCP is downward-looking
    
        Additional arguments (**kwargs):
            start_date (str, %Y%m%d %H:%M format): starting time of the period to extract
            end_date (str, %Y%m%d %H:%M format): end time of the period to extract
            
        Returns:
            True if the data was correctly read, False otherwise
        """
        
        log("Parsing data from {}.".format(file))

        try:
            dlfn_data = dlfn.read(file)
            #date = np.array([mplt_datetime(t) for t in dlfn_data.mpltime], dtype='datetime64[s]') # dolfyn version <1.0
            date= dlfn_data.time.data.astype('datetime64[s]')
            time = date.astype('int')

            if not isinstance(up, bool):
                up = bool(distutils.util.strtobool(up))
            if not isinstance(cabled, bool):
                cabled = bool(distutils.util.strtobool(cabled))

            # Define the measurement period either from correlation>20% or from specified start and end dates:
                
            #idx = np.where(np.nanmean(np.nanmean(dlfn_data.signal.corr/255.*100, axis=0), axis=0) > 20)[0] # dolfyn version <1.0
            idx = np.where(np.nanmean(np.nanmean(dlfn_data.corr/255.*100, axis=0), axis=0) > 20)[0]
            start = date[idx[[0, -1]]]

            date_subset = []
            if "start_date" in kwargs:
                start_date = datetime.strptime(kwargs["start_date"], "%Y%m%d %H:%M")
                date_subset.append(max([start_date, start[0]]))
            else:
                date_subset.append(start[0])
            if "end_date" in kwargs:
                end_date = datetime.strptime(kwargs["end_date"], "%Y%m%d %H:%M")
                date_subset.append(min([end_date, start[1]]))
            else:
                date_subset.append(start[1])

            time_subset = np.array(date_subset).astype('int')
            idx_subset = (time_subset[0] < time) & (time < time_subset[1])
            #dlfn_data = dlfn_data.subset[idx_subset] # dolfyn version <1.0
            # dlfn_data.date = date[idx_subset]
            # dlfn_data.timestamp = time[idx_subset]
            # dlfn_data.transducer_depth = transducer_depth
            
            # Create a new dataset with the subset data:
            dlfn_subset=xr.Dataset(coords={key: value.values for key, value in dlfn_data.coords.items() if key != "time"})
            dlfn_subset=dlfn_subset.expand_dims({"time":dlfn_data.time.values[idx_subset]})
            dlfn_subset=dlfn_subset.assign(timestamp=(["time"],time[idx_subset]))
            for var in dlfn_data.variables:
                dimvar=list(dlfn_data[var].dims)
                if "time" in dimvar:
                    print(var)
                    dlfn_subset=dlfn_subset.assign({var:(dimvar,dlfn_data[var].isel(time=idx_subset).values)})
            #self.general_attributes["Er"] = np.nanmin(dlfn_data.signal.echo.astype(float)) # dolfyn version <1.0
            self.general_attributes["Er"] = np.nanmin(dlfn_subset.amp.values.astype(float))
            self.general_attributes["cabled"] = str(cabled)
            self.general_attributes["up"] = str(up)
            self.general_attributes["bottom_depth"] = bottom_depth
            self.general_attributes["transducer_depth"] = transducer_depth
            # self.general_attributes["beam_angle"] = float(dlfn_data.config.beam_angle) # dolfyn version <1.0
            # self.general_attributes["xmit_length"] = float(dlfn_data.config.xmit_pulse) # transmit pulse length used for backscattering calculation, dolfyn version <1.0
            # self.general_attributes["beam_freq"] = float(dlfn_data.config.beam_freq_khz) # dolfyn version <1.0
            self.general_attributes["xmit_length"] = float(dlfn_data.attrs["transmit_pulse_m"]) # transmit pulse length [m] used for backscattering calculation 
            self.general_attributes["beam_angle"] = float(dlfn_data.attrs["beam_angle"]) 
            self.general_attributes["blank_dist"] = float(dlfn_data.attrs["blank_dist"])
            self.general_attributes["beam_freq"] = float(dlfn_data.attrs["freq"])
            
            if up:
                z0 = transducer_depth - dlfn_subset.range.values
            else:
                z0 = transducer_depth + dlfn_subset.range.values
            #echo = dlfn_data.signal.echo.astype(float) # dolfyn version <1.0
            echo = dlfn_subset.amp.values.astype(float)
            self.data["r"] = z0/np.cos(self.general_attributes["beam_angle"]*np.pi/180.)
            self.data["eu"] = dlfn_subset.vel[3, :, :].values # Velocity error
            #self.data["corr"] = dlfn_data.signal.corr.astype(float)/255. # dolfyn version <1.0
            #self.data["prcnt_gd"] = dlfn_data.signal.prcnt_gd.astype(float) # dolfyn version <1.0
            # self.data["heading"] = dlfn_data.orient.raw.heading # dolfyn version <1.0
            # self.data["roll"] = dlfn_data.orient.raw.roll # dolfyn version <1.0
            # self.data["pitch"] = dlfn_data.orient.raw.pitch # dolfyn version <1.0
            # self.data["time"] = np.array(dlfn_data.timestamp, dtype='float') # dolfyn version <1.0
            self.data["corr"] = dlfn_subset.corr.values.astype(float)/255.
            self.data["corr1"] = self.data["corr"][0, :, :]
            self.data["corr2"] = self.data["corr"][1, :, :]
            self.data["corr3"] = self.data["corr"][2, :, :]
            self.data["corr4"] = self.data["corr"][3, :, :]
            self.data["prcnt_gd"] = dlfn_subset.prcnt_gd.values.astype(float)
            self.data["prcnt_gd1"] = self.data["prcnt_gd"][0, :, :]
            self.data["prcnt_gd2"] = self.data["prcnt_gd"][1, :, :]
            self.data["prcnt_gd3"] = self.data["prcnt_gd"][2, :, :]
            self.data["prcnt_gd4"] = self.data["prcnt_gd"][3, :, :]
            self.data["heading"] = dlfn_subset.heading.values
            self.data["roll"] = dlfn_subset["roll"].values
            self.data["pitch"] = dlfn_subset.pitch.values
            self.data["time"] = dlfn_subset.timestamp.values.astype(float)
            self.data["depth"] = z0
            # self.data["u"] = dlfn_data.u # dolfyn version <1.0
            # self.data["v"] = dlfn_data.v # dolfyn version <1.0
            # self.data["w"] = dlfn_data.w # dolfyn version <1.0
            self.data["u"] = dlfn_subset.vel[0, :, :].values # Eastward velocity
            self.data["v"] = dlfn_subset.vel[1, :, :].values # Northward velocity
            self.data["w"] = dlfn_subset.vel[2, :, :].values # Upward velocity
            self.data["echo"] = echo
            self.data["echo1"] = echo[0, :, :]
            self.data["echo2"] = echo[1, :, :]
            self.data["echo3"] = echo[2, :, :]
            self.data["echo4"] = echo[3, :, :]
            # self.data["battery"] = dlfn_data.sys.adc[1, :] # dolfyn version <1.0 # dolfyn version <1.0
            # self.data["temp"] = dlfn_data.env.temperature_C # dolfyn version <1.0 # dolfyn version <1.0
            # Battery and temperature_C not available in the new dolfyn version
            self.data["battery"]=np.full(self.data["roll"].shape,np.nan)
            self.data["temp"]=np.full(self.data["roll"].shape,np.nan)

            return True
        except:
            log("Failed to process {}.".format(file))
            return False

    def quality_flags(self, envass_file = './quality_assurance.json', adcp_file='./quality_specific_adcp.json', simple=True):

        log("Performing quality assurance")
        log("1. ADCP-specific quality checks",indent=1) # Additional ADCP tests on the velocity matrix, increases qa with a base-2 approach (check#2 returns 0 or 2, chech#3 returns 0 or 4, etc.)
        quality_adcp_dict = json_converter(json.load(open(adcp_file))) # Load parameters related to simple and advanced quality checks
        
        varname0=quality_adcp_dict["variables"][0]
        qa_adcp=init_flag_adcp(np.array(self.data[varname0])) # Initial qa array (zero values)
        
        if "interface" in quality_adcp_dict["tests"].keys():
            if self.general_attributes['up']=='True': # Upward looking: surface detection
                qa_adcp=qa_adcp_interface_top(qa_adcp,self.data["depth"],self.general_attributes['transducer_depth'],beam_angle=self.general_attributes["beam_angle"])
            else: # Downward looking: sediment detection
                qa_adcp=qa_adcp_interface_bottom(qa_adcp,self.data["depth"],self.general_attributes['transducer_depth'],self.general_attributes['bottom_depth'],beam_angle=self.general_attributes["beam_angle"])
        
        if "corr" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_corr(qa_adcp,self.data["corr1"],self.data["corr2"],self.data["corr3"],self.data["corr4"],corr_threshold=quality_adcp_dict["tests"]["corr"]["corr_threshold"])
        
        if "PG14" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_PG14(qa_adcp,self.data["prcnt_gd1"],self.data["prcnt_gd4"],percentage_threshold=quality_adcp_dict["tests"]["PG14"]["percentage_threshold"])
        
        if "PG3" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_PG3(qa_adcp,self.data["prcnt_gd3"],percentage_threshold=quality_adcp_dict["tests"]["PG3"]["percentage_threshold"])
        
        if "velerror" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_velerror(qa_adcp,self.data["eu"],vel_threshold=quality_adcp_dict["tests"]["velerror"]["vel_threshold"])
        
        if "tilt" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_tilt(qa_adcp,self.data["roll"],self.data["pitch"],tilt_threshold=quality_adcp_dict["tests"]["tilt"]["tilt_threshold"])
        
        if "corrstd" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_corrstd(qa_adcp,self.data["corr1"],self.data["corr2"],self.data["corr3"],self.data["corr4"],std_threshold=quality_adcp_dict["tests"]["corrstd"]["std_threshold"])
        
        if "echodiff" in quality_adcp_dict["tests"].keys():
            qa_adcp=qa_adcp_echodiff(qa_adcp,self.data["echo1"],self.data["echo2"],self.data["echo3"],self.data["echo4"],diff_threshold=quality_adcp_dict["tests"]["echodiff"]["diff_threshold"])
        
        log("2. envass quality checks",indent=1) # Corresponds to quality check #1: qa is 0 (all good) or 1 (flagged)
        quality_assurance_dict = json_converter(json.load(open(envass_file))) # Load parameters related to simple and advanced quality checks
        
        for key, values in self.variables.copy().items():
            if (key in quality_assurance_dict):
                name = key + "_qual" # Create a new variable _qual to flag the data
                self.variables[name] = {'var_name': name, 'dim': values["dim"],
                                        'unit': '0 = nothing to report, 1 = more investigation',
                                        'long_name': name, }
                if key in quality_adcp_dict["variables"] and self.data[key].shape==qa_adcp.shape:
                    prior_qa=qa_adcp
                else:
                    prior_qa=np.zeros(self.data[key].shape)
                if simple: # Simple quality check only
                    self.data[name] = prior_qa+qualityassurance(np.array(self.data[key]), np.array(self.data["time"]), **quality_assurance_dict[key]["simple"])
                else:
                    quality_assurance_all = dict(quality_assurance_dict[key]["simple"], **quality_assurance_dict[key]["advanced"])
                    self.data[name] = prior_qa+qualityassurance(np.array(self.data[key]), np.array(self.data["time"]), **quality_assurance_all)
               


    def derive_variables(self, rotate_velocity):
        log("Computing derived variables.", indent=1)
        self.variables.update(self.derived_variables)

        log("Computes mean velocities", indent=2)
        self.data["mu"] = np.nanmean(self.data["u"], axis=0)
        self.data["mv"] = np.nanmean(self.data["v"], axis=0)
        self.data["mU"] = (self.data["mu"] ** 2 + self.data["mv"] ** 2) ** 0.5
        self.data["mdir"] = np.arctan2(self.data["mv"], self.data["mu"]) * 180 / np.pi
        self.data["mdir"][self.data["mdir"] < 0] = 360 + self.data["mdir"][self.data["mdir"] < 0]

        log("Compute rotate velocity", indent=2)
        self.data["u"], self.data["v"] = perform_rotate_velocity(self.data["u"], self.data["v"], rotate_velocity)

        log("Smooth data with moving average filter", indent=2)
        self.data["u"] = moving_average_filter(self.data["u"])
        self.data["v"] = moving_average_filter(self.data["v"])
        self.data["w"] = moving_average_filter(self.data["w"])

        log("Absolute backscatter", indent=2)
        # self.data["Sv"] = absolute_backscatter(self.data["echo"], self.data["temp"],
        #                                        self.general_attributes["beam_freq"],
        #                                        self.general_attributes["beam_angle"],
        #                                        bool(distutils.util.strtobool(self.general_attributes["cabled"])),
        #                                        self.data["depth"], self.data["r"],
        #                                        self.general_attributes["xmit_length"], self.data["battery"],
        #                                        self.general_attributes["Er"])
        self.data["Sv"] = np.mean(self.data["echo"], axis=0)