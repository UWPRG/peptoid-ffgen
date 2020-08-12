import numpy as np
import os
import pandas as pd

##### General Utilities #####
# Helper functions for core.py
#
#
# coterminal/coterminal_rev are used to deal with weird gaussian outputs
def coterminal(cv_deg):
    """
    Takes an input CV range from -180 to 180 and converts to 0 to 360

    Inputs
    cv_deg : array
        Array of degree values to change
        * Warning will only deal with range from -180 to 180
    """
    cv_new = []
    for cv in cv_deg:
        if cv < 0:
            cv+=360
            cv_new.append(cv)
        else:
            cv_new.append(cv)
    cv_new = np.array(cv_new)
    return cv_new

def coterminal_rev(cv_deg):
    """
    Inputs
    Takes an input CV range including numbers between 360->540  and
    converts them to -180 to 180. Does not handle neg numbers

    cv_deg : array
        Array of degree values to change
    """
    cv_new = []
    for cv in cv_deg:
        if cv > 180:
            cv-=360
            cv_new.append(cv)
        else:
            cv_new.append(cv)
    cv_new = np.array(cv_new)
    return cv_new


def scan_to_dataframe(scan_coord, scan_energy, min_to_zero=True):
    """
    Saves output of parse_logfile to dataframe structure with new units

    Inputs
    --------
    log2xyz_outfile : txt file
        Output of log2xyz function
    coord : array
        Array describing scan coordinates (not printed out from log2xyz)
    min_to_zero : bool
        Will normalize energy units by minimum energy in scan

    Returns
    --------
    scan_df : pandas dataframe
        dataframe containing energy and scan coordinates in different forms
    """
    ht_to_kcal = 627.509
    ht_to_kj = 2625.50

    if min_to_zero:
        scan_energy = scan_energy-np.min(scan_energy)
    else:
        raise Exception("Proceed with caution, energies are not normalized")

    scan_df = pd.DataFrame(scan_energy, columns=['hartree'])
    convert_to_360 = u.coterminal(scan_coord)
    scan_df['coord'] = convert_to_360
    scan_df['coord_rev'] = u.coterminal_rev(convert_to_360)
    scan_df['kj'] = scan_df['hartree']*ht_to_kj
    scan_df['kcal'] = scan_df['hartree']*ht_to_kcal

    return scan_df

def check_path_exists(file):
    directory = file.split(sep='/')[0] # only works on first dir in path
    if not os.path.exists(directory):
        os.makedirs(directory)
