import os
import copy
import numpy as np
from datetime import datetime

def adcp_flag_indices(data_dict,gen_att,quality_adcp_dict,instrument_parameters=dict()):
    """
    Get the flag indices for all the ADCP-specific quality checks.

    Parameters:
        data_dict (dict): data variables from the ADCP object
        gen_att (dict): general attributes from the ADCP object
        instrument_parameters (dict): instrument and period-specific parameters (from "parameters.json" file)

    Returns:
        qa_adcp (numpy array of floats): flag indices with a base 2 format
    """
    
    varname0=quality_adcp_dict["variables"][0]
    qa_adcp=init_flag_adcp(np.array(data_dict[varname0])) # Initial qa array (zero values)
    
    if "interface" in quality_adcp_dict["tests"].keys():
        if gen_att['up']=='True': # Upward looking: surface detection
            qa_adcp=qa_adcp_interface_top(qa_adcp,data_dict["depth"],gen_att['transducer_depth'],beam_angle=gen_att["beam_angle"])
        else: # Downward looking: sediment detection
            qa_adcp=qa_adcp_interface_bottom(qa_adcp,data_dict["depth"],gen_att['transducer_depth'],gen_att['bottom_depth'],beam_angle=gen_att["beam_angle"])
    
    if "corr" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"corr")
        qa_adcp=qa_adcp_corr(qa_adcp,data_dict["corr1"],data_dict["corr2"],data_dict["corr3"],data_dict["corr4"],corr_threshold=test_param["corr_threshold"])
    
    if "PG14" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"PG14")
        qa_adcp=qa_adcp_PG14(qa_adcp,data_dict["prcnt_gd1"],data_dict["prcnt_gd4"],percentage_threshold=test_param["percentage_threshold"])
    
    if "PG3" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"PG3")
        qa_adcp=qa_adcp_PG3(qa_adcp,data_dict["prcnt_gd3"],percentage_threshold=test_param["percentage_threshold"])
    
    if "velerror" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"velerror")
        qa_adcp=qa_adcp_velerror(qa_adcp,data_dict["eu"],vel_threshold=test_param["vel_threshold"])
    
    if "tilt" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"tilt")
        qa_adcp=qa_adcp_tilt(qa_adcp,data_dict["roll"],data_dict["pitch"],tilt_threshold=test_param["tilt_threshold"])
    
    if "corrstd" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"corrstd")
        qa_adcp=qa_adcp_corrstd(qa_adcp,data_dict["corr1"],data_dict["corr2"],data_dict["corr3"],data_dict["corr4"],std_threshold=test_param["std_threshold"])
    
    if "echodiff" in quality_adcp_dict["tests"].keys():
        test_param=get_test_param(instrument_parameters,quality_adcp_dict,"echodiff")
        qa_adcp=qa_adcp_echodiff(qa_adcp,data_dict["echo1"],data_dict["echo2"],data_dict["echo3"],data_dict["echo4"],diff_threshold=test_param["diff_threshold"])
        
    return qa_adcp
    
def get_test_param(p,quality_adcp_dict,test_name):
    """
    Get the threshold value for the given quality test, either from the instrument parameters (if present) or from the default values specified in quality_adcp_dict.

    Parameters:
        p (dict): instrument and period-specific parameters (from "parameters.json" file)
        quality_adcp_dict (dict): parameters of the ADCP-specific quality checks (from "quality_specific_adcp.json" file)
        test_name (string): name of the quality test
    Returns:
        test_param (dict): parameters of the specified quality test
    """
    
    if test_name in p.keys():
        test_param=p[test_name]
    else:
        test_param=quality_adcp_dict["tests"][test_name]
    
    return test_param

def init_flag_adcp(data_array):
    """
    Initialize the flag array.

    Parameters:
        data_array (np.array of floats): data array to which the quality assurance is applied.
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(data_array.shape)
    return flags

def qa_adcp_interface_top(prior_flags,depthval,depthADCP,beam_angle,flag_nb=2):
    """
    Flag data affected by the surface (sidelobe interference).

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        depthval (np.array of floats): positive depth values, already corrected for the ADCP location [m]
        depthADCP (float): depth of the ADCP [m]
        beam_angle (float): angle of the beams of the ADCP to compute sidelobe interference [°]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    dist_sidelobe=depthADCP*(1-np.cos(beam_angle*np.pi/180))
    flags=np.zeros(prior_flags.shape)
    indtop=np.where(depthval<dist_sidelobe)[0][0]-1 # Add one cell to be conservative, depth values are decreasing for upward ADCP
    flags[indtop:,:]=flag_nb
    flags=flags+prior_flags

    return flags

def qa_adcp_interface_bottom(prior_flags,depthval,depthADCP,depth_bottom,beam_angle,flag_nb=2):
    """
    Flag data affected by the sediment interface (sidelobe interference).

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        depthval (np.array of floats): positive depth values, already corrected for the ADCP location [m]
        depthADCP (float): depth of the ADCP [m]
        depth_bottom (float): depth of the sediment interface [m]
        beam_angle (float): angle of the beams of the ADCP to compute sidelobe interference [°]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    dist_sidelobe=(depth_bottom-depthADCP)*(1-np.cos(beam_angle*np.pi/180))
    flags=np.zeros(prior_flags.shape)
    indbot=np.where(depthval>(depth_bottom-dist_sidelobe))[0][0]-1 # Add one cell to be conservative
    flags[indbot:,:]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_corr(prior_flags,corr1,corr2,corr3,corr4,corr_threshold=64,flag_nb=2**2):
    """
    Flag data if correlation<threshold for at least one beam.

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        corr1 (np.array of floats): correlation values of beam 1 [0-1]
        corr2 (np.array of floats): correlation values of beam 2 [0-1]
        corr3 (np.array of floats): correlation values of beam 3 [0-1]
        corr4 (np.array of floats): correlation values of beam 4 [0-1]
        corr_threshold (float): minimum correlation for each beam [counts]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    corr_threshold=corr_threshold/255 # Between 0 and 1
    flags=np.zeros(prior_flags.shape)
    flags[(corr1<corr_threshold)|(corr2<corr_threshold)|(corr3<corr_threshold)|(corr4<corr_threshold)]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_PG14(prior_flags,PG1,PG4,percentage_threshold=25,flag_nb=2**3):
    """
    Flag data if PG1+PG4<threshold (i.e., low percentage of pings are used).

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        PG1 (np.array of floats): Percentage Good 1, percentage of samples/pings of the ensemble that are calculated with 3 beams because one beam was flagged “bad” [%]
        PG4 (np.array of floats): Percentage Good 4, percentage of samples/pings of the ensemble that are calculated with 4 beams [%]
        percentage_threshold (float): minimum percentage for PG1+PG4 [%]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    flags[PG1+PG4<percentage_threshold]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_PG3(prior_flags,PG3,percentage_threshold=25,flag_nb=2**4):
    """
    Flag data if PG3>threshold (i.e., high percentage of pings have been removed).

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        PG3 (np.array of floats): Percentage Good 3, percentage of samples/pings of the ensemble that have been removed by the software due to low correlation, high velocity error or fish echo [%]
        percentage_threshold (float): maximum percentage for PG3 [%]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    flags[PG3>percentage_threshold]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_velerror(prior_flags,vel_error,vel_threshold=0.05,flag_nb=2**5):
    """
    Flag data if velocity error>threshold.

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        vel_error (np.array of floats): velocity error [m/s]
        vel_threshold (float): maximum velocity error allowed [m/s]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    flags[np.abs(vel_error)>vel_threshold]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_tilt(prior_flags,roll,pitch,tilt_threshold=15,flag_nb=2**6):
    """
    Flag data if tilt angle>threshold.

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        roll (np.array of floats): roll angle time series [°]
        pitch (np.array of floats): pitch angle time series [°]
        tilt_threshold (float): maximum tilt allowed [°]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    flags[:,(np.abs(roll)>tilt_threshold)|(np.abs(pitch)>tilt_threshold)]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_corrstd(prior_flags,corr1,corr2,corr3,corr4,std_threshold=0.01,flag_nb=2**7):
    """
    Flag data if the std of the 4 correlations is > threshold (indicates that the beams measured something different, which is not normal)

    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        corr1 (np.array of floats): correlation values of beam 1 [0-1]
        corr2 (np.array of floats): correlation values of beam 2 [0-1]
        corr3 (np.array of floats): correlation values of beam 3 [0-1]
        corr4 (np.array of floats): correlation values of beam 4 [0-1]
        std_threshold (float): maximum value allowed for the standard deviation of the 4 correlations [0-1]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    corr_std=np.std(np.stack([corr1,corr2,corr3,corr4],axis=0), axis=0)
    flags[corr_std>std_threshold]=flag_nb
    flags=flags+prior_flags 
    return flags

def qa_adcp_echodiff(prior_flags,echo1,echo2,echo3,echo4,diff_threshold=30,flag_nb=2**8):
    """
    Flag data if there is a vertical difference in echo > threshold (obstacle)
    Parameters:
        prior_flags (np.array of ints): flags array with existing flags 
        echo1 (np.array of floats): echo values of beam 1 [counts]
        echo2 (np.array of floats): echo values of beam 2 [counts]
        echo3 (np.array of floats): echo values of beam 3 [counts]
        echo4 (np.array of floats): echo values of beam 4 [counts]
        diff_threshold (float): maximum value allowed for the echo different between two bins [counts]
        flag_nb (int): index of the flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    flags=np.zeros(prior_flags.shape)
    decho1=np.concatenate((np.zeros((1,echo1.shape[1])),np.diff(echo1,axis=0)),axis=0)
    decho2=np.concatenate((np.zeros((1,echo2.shape[1])),np.diff(echo2,axis=0)),axis=0)
    decho3=np.concatenate((np.zeros((1,echo3.shape[1])),np.diff(echo3,axis=0)),axis=0)
    decho4=np.concatenate((np.zeros((1,echo4.shape[1])),np.diff(echo4,axis=0)),axis=0)
    flags[(decho1>diff_threshold)|(decho2>diff_threshold)|(decho3>diff_threshold)|(decho4>diff_threshold)]=flag_nb
    flags=flags+prior_flags 
    return flags
