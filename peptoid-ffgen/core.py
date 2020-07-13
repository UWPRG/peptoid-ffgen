import pandas as pd
import numpy as np

def coterminal(cv_deg):
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
    cv_new = []
    for cv in cv_deg:
        if cv > 180:
            cv-=360
            cv_new.append(cv)
        else:
            cv_new.append(cv)
    cv_new = np.array(cv_new)
    return cv_new

def xyz_to_scan(file_name, coord):
    ht_to_kcal = 627.509
    ht_to_kj = 2625.50
    
    with open(file_name) as f:
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
