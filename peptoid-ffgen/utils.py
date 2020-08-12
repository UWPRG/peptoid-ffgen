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


def xyz_to_scan(log2xyz_outfile, coord):
    """
    Pulls energies from output of log2xyz (v1) and saves to data frame
    with different energy units
    *** Probably will be rewritten for updated code ***

    Inputs
    log2xyz_outfile : txt file
        Output of log2xyz function
    coord : array
        Array describing scan coordinates (not printed out from log2xyz)
    """
    ht_to_kcal = 627.509
    ht_to_kj = 2625.50

    with open(log2xyz_outfile) as f:
        f_lines = []
        scan_E = []
        for ln in f:
            if ln.startswith("i"):
                scan_E.append(float(ln.split()[5]))

    print(type(scan_E))
    scan = pd.DataFrame(scan_E, columns=['hartree'])
    convert_to_360 = coterminal(coord)
    scan['coord'] = convert_to_360
    scan['coord_rev'] = coterminal_rev(convert_to_360)
    scan['kj'] = scan['hartree']*ht_to_kj
    scan['kcal'] = scan['hartree']*ht_to_kcal

    return scan

def check_path_exists(file):
    directory = file.split(sep='/')[0] # only works on first dir in path
    if not os.path.exists(directory):
        os.makedirs(directory)
