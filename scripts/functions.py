import os
import copy
import json
import numpy as np
from datetime import datetime
from envass import qualityassurance


def fixed_grid_resample_guide(data, grid):
    resample = []
    for g in grid:
        for j in range(len(data)):
            if data[j] > g or j >= len(data) - 1:
                resample.append({"index": False})
                break
            elif data[j] <= g < data[j + 1]:
                itp = (g - data[j]) / (data[j + 1] - data[j])
                resample.append({"index": j, "interpolation": itp})
                break
    return resample

def resample(guide, data):
    out = []
    for i in range(len(guide) - 1):
        if not guide[i]["index"]:
            out.append(np.nan)
        else:
            value = ((data[guide[i]["index"] + 1] - data[guide[i]["index"]]) * guide[i]["interpolation"]) + data[guide[i]["index"]]
            out.append(value)
    out.append(np.nan)
    return out
    
def rotation_matrix_2d(alpha):
    M = np.zeros((2,2))
    M[0,0] = np.cos(alpha)
    M[0,1] = -np.sin(alpha)
    M[1,0] = np.sin(alpha)
    M[1,1] = np.cos(alpha)
    return M

def perform_rotate_velocity(u,v, alpha = 40):
    #rotates the velocities by an angle alpha (degrees)
    for index in np.ndindex(u.shape):
        M = rotation_matrix_2d(alpha/180*np.pi)
        V = np.array([[u[index],v[index]]]).T
        V2 = np.dot(M,V)
        u[index] = V2[0]
        v[index] = V2[1]
    return u,v 

def moving_average_filter(array, m=3, n=7, valid_entries=3):
    ## Filtering velocities by a moving average filter, 
    ## keeping only elements that do not relay on zero-padding

    # Define filter size: 4m vertical (m=3 element), 1h horizontal (n=7, elements)
    y, x = array.shape
    y = y - m + 1
    x = x - n + 1

    u_smoothed = np.full(array.shape, np.nan)
    
    # Define minimum number of valid entries needed to perform the smoothing
    valid_entries = 3
    
    # Filter velocities
    for i in range(y):
        for j in range(x):
            
            # Mask data to only sum the valid entries later
            masked_u = np.ma.masked_invalid(array[i:i+m, j:j+n])

            # Count number of valid entries
            count_u = len(masked_u.compressed())

            # Perform smoothing only if minimum number of valid entries is satisfied, otherwise set to NaN
            if count_u >= valid_entries:
                u_smoothed[i][j] = np.sum(masked_u)/count_u
            else:
                u_smoothed[i][j] = np.nan

            
    #self.data["depth"] = self.z0[1:-1]
    #self.data["time"] = self.time[3:-3]
    # Adapt time and depth array
    #self.z0 = self.z0[1:-1]
    #self.time = self.time[3:-3]
    #self.date = self.date[3:-3]
    
    # Adapt variables for calculating backscattering later
    #self.r = self.r[1:-1]
    #self.temp = self.temp[3:-3]
    #self.batery = self.battery[3:-3]
    #self.echo = self.echo[:, 1:-1, 3:-3]
    return u_smoothed

def absolute_backscatter(echo, temp, beam_freq, beam_angle, cabled, z0, r, xmit_length, battery, Er,
                                kc = 0.45, PLOT = False, msv = -80, Msv = -55, dsv = 0.5):
    """
    Computes absolute backscatter Sv [dB] from echo [counts] for RDI Workhorse Signature.

    Parameters:
        echo (4*n_bins*n_time np.array of floats): echo data for the 4 beams [counts]
        temp (n_time np.array of floats): temperature [°C]
        beam_freq (float): acoustic frequency provided by dolfyn [kHz]
        beam_angle (float): angle between the transmitted signal and the ADCP z-axis, provided by dolfyn [°]
        cabled (boolean): = True if cabled ADCP
        z0 (n_bins np.array of floats): depth values, positive downward [m]
        r (n_bins np.array of floats): radial distances r=depth/cos(beam_angle) [m]
        xmit_length (float): transmit pulse length provided by dolfyn [m]
        battery (n_time np.array of floats): battery level [counts]
        Er (float): noise echo in absence of any signal [counts]
        kc (float): conversion factor from counts to decibels [dB/count]
        PLOT (boolean): = True to make plot, NOT USED HERE
        msv (int): NOT USED HERE
        Msv (int): NOT USED HERE
        dsv (float): NOT USED HERE
        
        
    Returns:
        mean_Sv (n_bins*n_time np.array of floats): mean backscatter between the four beams [dB]
        
    """
    mbat = 32.
    if beam_freq == 300:
        C = -140.87
        if not cabled:
            Pdbw = 14.0
        else:
            Pdbw = 17.5
        alpha = 0.025
    elif beam_freq == 600:
        C = -139.09
        if not cabled:
            Pdbw = 9.0
        else:
            Pdbw = 12.5
        alpha = 0.098
                
    R = r + 0.25*(z0[1]-z0[0])/np.cos(beam_angle*np.pi/180.)
    Sv = np.full(echo.shape, np.nan)
    shp = echo.shape
    for i in range(shp[0]):
        for j in range(shp[1]):
            for k in range(shp[2]):
                Sv[i,j,k] = C +10*np.log10((temp[k]+273.16)*R[j]**2) - 10*np.log10(xmit_length) - Pdbw  + 2.*alpha*R[j]+ 10.*np.log10( 10**(kc*(echo[i,j,k]-Er)/10.)-1 )
                if not cabled:
                    Sv[i,j,k] -= 20*np.log10(battery[k]/np.nanmean(battery))
    mean_Sv = np.mean(Sv, axis=0)
    return mean_Sv
def finds_surface_1prof(echo, itime, irt = 100, factor = 1.5, PLOT = False):
    #This function just finds the surface for a given time step and given accoustic beam
    #It probably needs to be improved
    
    #uses the maximum echo
    d1,d2,d3 = echo.shape
    istart = itime-irt//2
    iend = itime+irt//2


    if istart<0:
        iend = iend - istart
        istart = 0
    if iend>d3:
        ic = iend-d3
        iend = d3
        istart = istart-ic

        
    itt = np.arange(istart,iend)
    ECHO = np.nanmedian(np.nanmedian(echo[:,:,itt], axis = 0), axis = 1)
    #ECHO = np.nanmedian(self.echo[:,:,itime], axis = 0)
    imin = np.argmin(ECHO)
    if imin >= ECHO.size-1:
        imin = 0
    minECHO = np.min(ECHO)
    ECHO = ECHO - minECHO
    maxECHO = np.max(ECHO[imin:])
    imax = np.argmax(ECHO[imin:])+imin
    isurfA = np.where( (ECHO[:imax+1]<maxECHO/factor) )[0]
    if isurfA.size>0:
        isurfA = isurfA[-1]
    else:
        isurfA = imax
        
    #uses the maximum change in echo
    diffECHO = np.diff(ECHO)
    maxdiffECHO = np.max(diffECHO[imin:])
    imaxdiff = np.where(diffECHO==maxdiffECHO)[0][0]
    isurfB = np.where( diffECHO<=maxdiffECHO)[0]
    if imaxdiff ==0:
        imaxdiff = diffECHO.size
    isurfB = isurfB[:imaxdiff][-1]
    
    
    #gets the minimum of the three
    isurf = int(max([isurfA, isurfB+1]))
    
    return isurf

def finds_surface_timeseries(echo, range, bottom_depth, up, irt = 100, factor = 1.5, PLOT = False):
    #looks where is the surface in each step, to be sure of the depth where the instrument is
    #in some of the deployements the instrument mooved to a deeper part of the lake and the surface was lost
    #i dont know how to deal with this yet
    print("Finding surfaces")
    d1,d2,d3 = echo.shape
    irt = min([irt,d3])
    isurf = np.full( d3, 0 ) #index of the surface for each beam
    rsurf = np.full( d3, np.nan) #distance from transducer to the surface
    z = np.full( (d2,d3), np.nan ) #this variable is the actual depth of each bin at each time step
    watercol = np.full ( (d2,d3), True ) #flags the data that is underwater
    for j in range(d3):
        print("Step %d of %d"%(j+1,d3))
        isurf[j] = finds_surface_1prof( itime = j, irt = irt, factor = factor,PLOT = PLOT)
        watercol[isurf[j]+1:,j] = False
        rsurf[j] = range[isurf[j]]
        if not up:
                r = (bottom_depth - rsurf[j])+range
        else:
                r = rsurf[j]-range
        z[:,j] = r #*np.cos(20*np.pi/180.) #this was incorrect
    return isurf, rsurf, z, watercol

def copy_variables(variables_dict):
    var_dict = dict()
    for var in variables_dict:
        var_dict[var] = variables_dict[var][:]
    nc_copy = copy.deepcopy(var_dict)
    return nc_copy

def mplt_datetime(t):
    return datetime.utcfromtimestamp((t - 719163) * 24 * 60 * 60)
    
def advanced_quality_flags(df, json_path="quality_assurance.json"):
    """
        input :
            - df is a dataframe of level 1B where basic check have been performed
            - json path: path for the advanced quality check json file, produced by the jupyter notebook
        output:
            - dictionnary where the dataframe is stored with updated advanced quality checks
        """
    quality_assurance_dict = json.load(open(json_path))
    var_name = quality_assurance_dict.keys()
    advanced_df = df.copy()
    for var in var_name:
        if quality_assurance_dict[var]:
            qa = qualityassurance(np.array(df[var]), np.array(df["time"]), **quality_assurance_dict[var]["advanced"])
            advanced_df[var + "_qual"].values[np.array(qa, dtype=bool)] = 1
    return advanced_df
    
def json_converter(qa):
    for keys in qa.keys():
        try:
            if qa[keys]["simple"]["bounds"][0] == "-inf":
                qa[keys]["simple"]["bounds"][0] = -np.inf
            if qa[keys]["simple"]["bounds"][1] == "inf":
                qa[keys]["simple"]["bounds"][1] = np.inf
        except:pass
    try:
        if qa["time"]["simple"]["bounds"][1] == "now":
            qa["time"]["simple"]["bounds"][1] = datetime.now().timestamp()
        return qa
    except:
        return qa
    
def log(str, indent=0, start=False):
    if start:
        out = "\n" + str + "\n"
        with open("log.txt", "w") as file:
            file.write(out + "\n")
    else:
        out = datetime.now().strftime("%H:%M:%S.%f") + (" " * 3 * (indent + 1)) + str
        with open("log.txt", "a") as file:
            file.write(out + "\n")
    print(out)


def error(str):
    out = datetime.now().strftime("%H:%M:%S.%f") + "   ERROR: " + str
    with open("log.txt", "a") as file:
        file.write(out + "\n")
    raise ValueError(str)


def find_closest_index(arr, value):
    return min(range(len(arr)), key=lambda i: abs(arr[i] - value))


def is_number(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return True


def isnt_number(n):
    try:
        float(n)
    except ValueError:
        return True
    else:
        return False


def select_parameters(file, parameters):
    if "RDI300" in file:
        instrument = "300"
    elif "RDI600" in file:
        instrument = "600"
    else:
        raise ValueError('Unrecognised file input string')
    try:
        dt = datetime.strptime(os.path.normpath(file).split(os.path.sep)[-2], '%Y%m%d')
        print(parameters)
        for i in range(len(parameters)):
            start = datetime.strptime(parameters[i]["start"], '%Y%m%d')
            if parameters[i]["end"] == "now":
                end = datetime.now()
            else:
                end = datetime.strptime(parameters[i]["end"], '%Y%m%d')
            if start <= dt <= end:
                return parameters[i]["data"][instrument]
        raise ValueError("Couldn't find parameters for files time period.")
    except:
        raise ValueError('Unable to parse data from filename')


def latest_files(indir, folder):
    dirs = next(os.walk(os.path.join(indir, folder)))[1]
    dirs = list(filter(filter_dir, dirs))
    if len(dirs) > 0:
        dirs.sort()
        latest_dir = dirs[-1]
        files = os.listdir(os.path.join(indir, folder, latest_dir))
        files = list(filter(filter_file, files))
        files = [os.path.join(indir, folder, latest_dir, i) for i in files]
        files.sort()
        return files
    else:
        return []


def all_files(indir, folder):
    filelist = []
    for root, dirs, files in os.walk(os.path.join(indir, folder)):
        for file in files:
            if file.endswith('.LTA'):
                filelist.append(os.path.join(root, file))
    filelist.sort()
    return filelist

def filter_dir(dir):
    if len(dir) == 8 and str(int(dir)) == dir:
        return True
    else:
        return False


def filter_file(file):
    if file.endswith(".LTA"):
        return True
    else:
        return False