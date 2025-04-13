import os
import copy
import numpy as np
from datetime import datetime

def init_flag_adcp(data_array,prior_flags=False):
    """
    Initialize the flag array.

    Parameters:
        data_array (np.array of floats): data array to which the quality assurance is applied.
        prior_flags (np.array of ints): flags array with existing flags or False if no existing flag
    Returns:
        flags (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """
    
    if prior_flags:
        flags=prior_flags.copy()
    else:
        flags=np.zeros(data_array.shape,dtype=bool)
    return flags

def qa_adcp_interface(a, prior_flags=False,flag_nb=1):
    """
    Flag data affected by the surface or the ground (sidelobe interference).

    Parameters:
        a (ADCP object): ADCP object containing ADCP data.
        prior_flags (np.array of ints): flags array with existing flags or False if no existing flag
        flag_nb (int): index of the flag
    Returns:
        flag (np.array of ints): data array where flagged data is shown with a number>0 and 0 = no flag.
    """

    flags = init_flag_adcp(a.echo1, prior_flags)
    
    # a.press # To be created
    # a.general_attributes['up']
    # a.general_attributes['bottom_depth']
    # a.general_attributes['transducer_depth']

    
    return flags