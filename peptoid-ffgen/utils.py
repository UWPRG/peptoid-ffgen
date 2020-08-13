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

    Parameters
    ----------
    cv_deg : array
        Array of degree values to change
        * Warning will only deal with range from -180 to 180

    Returns
    ----------
    cv_new : array
        degree + 360
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
    Parameters
    ----------
    Takes an input CV range including numbers between 360->540  and
    converts them to -180 to 180. Does not handle neg numbers

    cv_deg : array
        Array of degree values to change

    Returns
    ----------
    cv_new : array
        degree + 360
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

    Parameters
    ----------
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
        print("Proceed with caution, energies are not normalized to any value")

    scan_df = pd.DataFrame(scan_energy, columns=['hartree'])
    convert_to_360 = coterminal(scan_coord)
    scan_df['coord'] = convert_to_360
    scan_df['coord_rev'] = coterminal_rev(convert_to_360)
    scan_df['kj'] = scan_df['hartree']*ht_to_kj
    scan_df['kcal'] = scan_df['hartree']*ht_to_kcal

    return scan_df


def mdscan_to_dataframe(scan_coord, md_energy, min_to_zero=True):
    """
    Saves output of run_md w/ Plumed = False (default kJ) to dataframe
    structure with new units

    Parameters
    ----------
    log2xyz_outfile : txt file
        Output of log2xyz function
    coord : array
        Array describing scan coordinates (not printed out from log2xyz)
    min_to_zero : bool
        Will normalize energy units by minimum energy in scan

    Returns
    --------
    md_df : pandas dataframe
        dataframe containing energy and scan coordinates in different forms
    """
    ht_to_kcal = 627.509
    ht_to_kj = 2625.50
    kj_to_kcal = 1/4.184

    if min_to_zero:
        md_energy = md_energy-np.min(md_energy)
    else:
        print("Proceed with caution, energies are not normalized to any value")

    md_df = pd.DataFrame(md_energy, columns=['kj'])
    convert_to_360 = coterminal(scan_coord)
    md_df['coord'] = convert_to_360
    md_df['coord_rev'] = coterminal_rev(convert_to_360)
    md_df['hartree'] = md_df['kj']/ht_to_kj
    md_df['kcal'] = md_df['kj']*kj_to_kcal

    return md_df

def check_path_exists(file):
    """
    Creates directory of path to dir doesnt exist (Uses a "/" delimiter so only
    does this for first entry - needs to be updated)

    Parameters
    ----------
    file : str
        path to file (ex: xyz/xyz2 will create the directory xyz/ from working
        path)
    """
    if not os.path.exists(file):
        os.makedirs(file)

class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
