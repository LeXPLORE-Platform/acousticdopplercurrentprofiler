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


class ADCP():
    def __init__(self):
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
            'time': {'var_name': 'time', 'dim': ('time',), 'unit': 'seconds since 1970-01-01 00:00:00', 'longname': 'time'},
            'depth': {'var_name': 'depth', 'dim': ('depth',), 'unit': 'm', 'longname': 'nominal depth'},
            'u': {'var_name': 'u', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'longname': 'eastern velocity'},
            'v': {'var_name': 'v', 'dim': ('depth', 'time'), 'unit': 'm s-1', 'longname': 'northern velocty'},
            'temp': {'var_name': 'temp', 'dim': ('time',), 'unit': 'degC', 'longname': 'temperature', },
            'echo1': {'var_name': 'echo1', 'dim': ('depth', 'time'), 'unit': '-', 'longname': 'Beam 1 echo'},
            'echo2': {'var_name': 'echo2', 'dim': ('depth', 'time'), 'unit': '-', 'longname': 'Beam 2 echo'},
            'echo3': {'var_name': 'echo3', 'dim': ('depth', 'time'), 'unit': '-', 'longname': 'Beam 3 echo'},
            'echo4': {'var_name': 'echo4', 'dim': ('depth', 'time'), 'unit': '-', 'longname': 'Beam 4 echo'},
            'battery': {'var_name': 'battery', 'dim': ('time',), 'unit': '-', 'longname': 'Battery level'},
            'heading': {'var_name': 'heading', 'dim': ('time',), 'unit': 'deg', 'longname': 'Heading'},
            'roll': {'var_name': 'roll', 'dim': ('time',), 'unit': 'deg', 'longname': 'Roll'},
            'pitch': {'var_name': 'pitch', 'dim': ('time',), 'unit': 'deg', 'longname': 'Pitch'},
        }

        self.derived_variables = {
            'mU': {'var_name': 'mU', 'dim': ('time',), 'unit': 'm s-1', 'longname': 'modulus of depth-averaged velocity'},
            'mdir': {'var_name': 'mdir', 'dim': ('time',), 'unit': 'deg', 'longname': 'direction (anticlockwise from east) of depth-averaged velocity'},
            'Sv': {'var_name': 'Sv', 'dim': ('depth', 'time'), 'unit': '%', 'longname': 'absolute backscatter'}
        }

        self.data = {}

    def read_data(self, file, transducer_depth=0.65, bottom_depth=110., cabled=False, up=False, **kwargs):
        log("Parsing data from {}.".format(file))
        dlfn_data = dlfn.read(file)
        date = np.array([mplt_datetime(t) for t in dlfn_data.mpltime], dtype='datetime64[s]')
        time = date.astype('int')

        if not isinstance(up, bool):
            up = bool(distutils.util.strtobool(up))
        if not isinstance(cabled, bool):
            cabled = bool(distutils.util.strtobool(cabled))

        idx = np.where(np.nanmean(np.nanmean(dlfn_data.signal.corr/255.*100, axis=0), axis=0) > 20)[0]
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
        dlfn_data = dlfn_data.subset[idx_subset]
        dlfn_data.date = date[idx_subset]
        dlfn_data.timestamp = time[idx_subset]
        dlfn_data.transducer_depth = transducer_depth

        self.general_attributes["Er"] = np.nanmin(dlfn_data.signal.echo.astype(float))
        self.general_attributes["cabled"] = str(cabled)
        self.general_attributes["up"] = str(up)
        self.general_attributes["bottom_depth"] = bottom_depth
        self.general_attributes["transducer_depth"] = transducer_depth
        self.general_attributes["beam_angle"] = float(dlfn_data.config.beam_angle)
        self.general_attributes["xmit_length"] = float(dlfn_data.config.xmit_pulse)
        self.general_attributes["beam_freq"] = float(dlfn_data.config.beam_freq_khz)

        if up:
            z0 = transducer_depth - dlfn_data.range
        else:
            z0 = transducer_depth + dlfn_data.range
        echo = dlfn_data.signal.echo.astype(float)
        self.data["r"] = z0/np.cos(self.general_attributes["beam_angle"]*np.pi/180.)
        self.data["eu"] = dlfn_data.vel[3, :, :]
        self.data["corr"] = dlfn_data.signal.corr.astype(float)/255.
        self.data["prcnt_gd"] = dlfn_data.signal.prcnt_gd.astype(float)
        self.data["heading"] = dlfn_data.orient.raw.heading
        self.data["roll"] = dlfn_data.orient.raw.roll
        self.data["pitch"] = dlfn_data.orient.raw.pitch
        self.data["time"] = np.array(dlfn_data.timestamp, dtype='float')
        self.data["depth"] = z0
        self.data["u"] = dlfn_data.u
        self.data["v"] = dlfn_data.v
        self.data["w"] = dlfn_data.w
        self.data["echo"] = echo
        self.data["echo1"] = echo[0, :, :]
        self.data["echo2"] = echo[1, :, :]
        self.data["echo3"] = echo[2, :, :]
        self.data["echo4"] = echo[3, :, :]
        self.data["battery"] = dlfn_data.sys.adc[1, :]
        self.data["temp"] = dlfn_data.env.temperature_C

    def to_netcdf(self, folder, title, time_label="time", output_period="profile", mode='a', grid=False):
        
        if not os.path.exists(folder):
            os.makedirs(folder)

        if grid:
            variables = self.grid_variables
            dimensions = self.grid_dimensions
            data = self.grid
        else: 
            variables = self.variables
            dimensions = self.dimensions
            data = self.data

        time_arr = data[time_label]
        dt_min = datetime.utcfromtimestamp(np.nanmin(time_arr))
        dt_max = datetime.utcfromtimestamp(np.nanmax(time_arr))

        if type(output_period) == int: 
            start = (dt_min - dt.timedelta(days=dt_min.weekday())).replace(hour=0, minute=0, second=0, microsecond=0)
            td = dt.timedelta(days=output_period)
        elif output_period == "daily": 
            start = (dt_min - dt.timedelta(days=dt_min.weekday())).replace(hour=0, minute=0, second=0, microsecond=0)
            td = dt.timedelta(days=1)
        elif output_period == "weekly":
            start = (dt_min - dt.timedelta(days=dt_min.weekday())).replace(hour=0, minute=0, second=0, microsecond=0)
            td = dt.timedelta(weeks=1)
        elif output_period == "monthly":
            start = dt_min.replace(day=1, hour=0, minute=0, second=0)
            td = relativedelta(months=+1)
        elif output_period == "yearly":
            start = dt_min.replace(month=1, day=1, hour=0, minute=0, second=0)
            td = relativedelta(year=+1)
        elif output_period == "profile":
            start = dt_min
            td = dt_max-dt_min
        else:
            log("Output periods {} not defined.".format(output_period))



        while start < dt_max:
            end = start + td

            s = datetime.timestamp(start)
            e = datetime.timestamp(end)

            filename = "{}_{}.nc".format(title, start.strftime('%Y%m%d_%H%M%S'))
            out_file = os.path.join(folder, filename)
            log("Writing data from {} until {} to NetCDF file {}".format(start, end, filename), 2)

            if os.path.isfile(out_file):
                nc = netCDF4.Dataset(out_file, mode=mode, format='NETCDF4')

                nc_time_arr = np.array(nc.variables[time_label][:])

                if (np.all(np.isin(time_arr, nc_time_arr))) and (mode!="r+"):
                    log("Duplicated run, no data added", 2)
                    nc.close()
                    start = start + td
                    continue
                else:

                    valid_time = (time_arr >= s) & (time_arr < e)
                    valid_duplicates = ~np.isin(time_arr, nc_time_arr)
                    valid = np.logical_and(valid_time, valid_duplicates)
                    combined_time = np.append(nc_time_arr, time_arr[valid])
                    order = np.argsort(combined_time)
                    
                    nc_copy = copy_variables(nc.variables)

                    for key, values in variables.items():
                        nc_data = np.array(nc.variables[key][:])
                        var = data[key]
                        if "time" in values["dim"]:
                            if var.ndim==1:
                                combined = np.append(nc_data, var[valid])
                                out = combined[order]
                                nc.variables[key][:] = out
                            if var.ndim==2:
                                combined = np.hstack((nc_copy[key][:], var[:,valid]))
                                out = combined[:,order]
                                nc.variables[key][:] = out
                    nc.close()
            else:

                valid_time = (time_arr >= s) & (time_arr < e)
                nc = netCDF4.Dataset(out_file, mode='w', format='NETCDF4')

                for key in self.general_attributes:
                    setattr(nc, key, self.general_attributes[key])

                for key, values in dimensions.items():
                    nc.createDimension(values['dim_name'], values['dim_size'])

                for key, values in variables.items():
                    var = nc.createVariable(values["var_name"], np.float64, values["dim"], fill_value=np.nan)
                    var.units = values["unit"]
                    var.long_name = values["longname"]
                    if  time_label in values["dim"]: 
                        if var.ndim==1:
                            var[:] = data[key][valid_time]
                        if var.ndim==2:
                            var[:] = data[key][:,valid_time]
                    else: var[:] = data[key][:]
                nc.close()

            start = start + td

    def quality_flags(self, file_path = './quality_assurance.json', simple=True):

        log("Performing quality assurance")
        qa_adcp = self.quality_flags_adcp(mincor=0.3, maxcor=0.7, minpcg=50, errorfactor=1)

        quality_assurance_dict = json_converter(json.load(open(file_path)))

        for key, values in self.variables.copy().items():
            if (key in quality_assurance_dict):
                name = key + "_qual"
                self.variables[name] = {'var_name': name, 'dim': values["dim"],
                                        'unit': '0 = nothing to report, 1 = more investigation',
                                        'longname': name, }
                if simple:
                    self.data[name] = qualityassurance(np.array(self.data[key]), np.array(self.data["time"]), **quality_assurance_dict[key]["simple"])
                else:
                    quality_assurance_all = dict(quality_assurance_dict[key]["simple"], **quality_assurance_dict[key]["advanced"])
                    self.data[name] = qualityassurance(np.array(self.data[key]), np.array(self.data["time"]), **quality_assurance_all)
                if key == "u" or key == "v":
                    self.data[name][qa_adcp[key] > 0] = 1
                    log("Specific ADCP check done for {}".format(key))

    def quality_flags_adcp(self, mincor, maxcor, minpcg, errorfactor):
        d1, d2, d3 = self.data["corr"].shape
        quality_flag = np.full((d2, d3), False)
        for j in range(d3):
            for i in range(d2):
                flg = 0
                for k in range(d1):
                    cr = self.data["corr"][k, i, j]
                    if (cr < mincor) or (cr > maxcor):
                        flg += 1
                pcg = self.data["prcnt_gd"][3, i, j]
                if pcg < minpcg:
                    flg += 1
                modW = np.sqrt(self.data["u"][i, j]**2 + self.data["v"][i, j]**2 + self.data["w"][i, j]**2)**0.5
                if modW < errorfactor*np.abs(self.data["eu"][i, j]):
                    flg += 1
                if flg > 0:
                    quality_flag[i, j] = True
        qa_spe = {"u": quality_flag, "v": quality_flag}
        return qa_spe

    def mask_data(self):
        for var in self.variables:
            if var + "_qual" in self.data:
                idx = self.data[var + "_qual"][:] > 0
                self.data[var][idx] = np.nan

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

        log("Absolute backscatter", indent=2)
        self.data["Sv"] = absolute_backscatter(self.data["echo"], self.data["temp"],
                                               self.general_attributes["beam_freq"],
                                               self.general_attributes["beam_angle"],
                                               bool(distutils.util.strtobool(self.general_attributes["cabled"])),
                                               self.data["depth"], self.data["r"],
                                               self.general_attributes["xmit_length"], self.data["battery"],
                                               self.general_attributes["Er"])
