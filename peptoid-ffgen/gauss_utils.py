import os.path
import re
import sys
import argparse
import numpy as np

def fortrandouble(x):
    """Converts string of Fortran double scientific notation to python float.

    Example
    -------
    >>> val = '-0.12345D+03'
    >>> print(fortrandouble(val))
    -123.45
    """
    return float(x.replace('D', 'E'))


def check_post_hartree_fock(post_hartree_fock):
    """Check if post-HF method is being used and which one.

    **SA - SHOULD THIS BE UPDATED? in case B3LYP is used? **

    Parameters
    ----------
    post_hartree_fock : str
        String containing post-HF level (MP2 or MP3) else leave blank

    Returns
    -------
    using_post_HF : bool
        Boolean that denotes whether PHF method is used
    """
    post_HF_options = {'MP2': 'EUMP2', 'MP3': 'EUMP3'}
    if post_hartree_fock:
        assert post_hartree_fock in post_HF_options.keys(), \
            "{} not available. ".format(post_hartree_fock) +\
            "Must choose from: {} or leave blank".format(post_HF_options.keys())
        return True
    else:
        return False

def read_log_into_lists(logfile, post_hartree_fock, no_energy=False):
    """Read the log file into lists.

    Paramters
    ---------
    post_hartree_fock : str
	       String of post-HF level
    logfile : str
	   String with  full name of log file including extension
    no_energy : bool
        Enabling turns off recording of energy (Default = False)

    Returns
    -------
    cartcoords : array
        cartesian coordinates (x,y,z)
    step_num : list of lists
        **ADD DESCRIPTION**
    energy : list of lists
        **ADD DESCRIPTION**
    """
    post_HF_options = {'MP2': 'EUMP2', 'MP3': 'EUMP3'}
    cartcoords = []
    frame_num = 0
    step_num = []
    step_pattern = r'\d+'
    energy = []
    scf_pattern =r'([-]?\d+\.\d+)'
    mp_pattern = r'([-]?\d+\.\d+[DE][+-]\d+)'
    with open(logfile, 'r') as logf:
        for line in logf:
            if 'Input orientation' in line.strip():
                # Throw away the 4 header lines.
                for _ in range(4):
                    _ = next(logf)
                cartcoords.append([])
                line = logf.readline()
                while not '-----' in line:
                    cartcoords[frame_num].append(line.split())
                    line = logf.readline()
                frame_num += 1
            elif 'Step number' in line.strip():
                nums = re.findall(step_pattern, line)
                step_num.append([int(n) for n in nums])
            elif not no_energy:
                # Set parameters for grepping the energy from log file
                if post_hartree_fock:
                    grep_string = post_HF_options[post_hartree_fock]
                    ener_pattern = mp_pattern
                else:
                    grep_string = 'SCF Done'
                    ener_pattern = scf_pattern
                if grep_string in line.strip():
                    ener_string = re.findall(ener_pattern, line)[-1]
                    ener_float = fortrandouble(ener_string)
                    energy.append(ener_float)
            else:
                pass


    # Prepare data for writing to file ----------------------------------------
    # Delete first and third columns from coordinates.
    cartcoords = np.array(cartcoords, dtype=np.float)
    cartcoords = np.delete(cartcoords, [0, 2] ,axis=2)
    return cartcoords, step_num, energy


def determine_and_save_frames(step_num, full=False, scan=True):
    """
    Determine which frames to write to file

    ** SA - this function isn't called right now **

    Parameters
    ----------
    scan : bool
	   log file corresponds to scan (Default: True)
    full : bool
	   option to save all frames (including non-optimized steps) (True) or to
       only save optimized coordinates (False) (Default: False)
    step_num : list of lists
        ** ADD DESCRIPTION **
    Returns
    -------
    save_frames : list

    """
    #Determine which frames to write to file.
    if scan and not full:
        # `save_frames` is a list of frame indices to be saved, not actual
        # frame numbers.
        save_frames = []
        for i, step in enumerate(step_num[:-1]):
            # Only pull the last frame of each scan point.
            if step[2] != step_num[i + 1][2]:
                save_frames.append(i)
            else:
                pass
        # Save last frame as well.
        save_frames.append(len(step_num) - 1)
    else:
        if full:
            # Save all frames, even if it is a scan calculation.
            pass
        else:
            # Check if step numbers repeat (i.e., if calculation is a scan).
            unique_step_num = np.unique([step[-2] for step in step_num])
            # Allow for 1 or 2 extra frames associated with freq calculations.
            if len(step_num) - len(unique_step_num) < 2:
                pass
            else:
                raise Exception("Calculation is likely a scan. Recheck "
                                "arguments.\nIf so, use scan== True")
        save_frames = range(len(step_num))
    return save_frames


def write_to_file(cartcoords, step_num, energy, save_frames, outfile="", overwrite=True, last_frame=False, no_energy=False):
    """
    Write to file

    Paramters
    ---------
    cartcoords : array
        cartesian coordinates (x,y,z)
    step_num : list of lists
        ** ADD DESCRIPTION **
    energy : array
        ** ADD DESCRIPTION **
    save_frames : list
        ** ADD DESCRIPTION **
    outfile : str
	   str that contains name of outputfile without extension (Default: "" no extension)
    overwrite : bool
        denotes whether or not output file overwrites previous file of same name (Default = True)
    last_frame: bool
        denotes whether to save only the last frame (Default = False)
    no_energy : bool
        Enabling turns off recording of energy (Default = False)

    """
    element_dict = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
                    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 16: 'S'}
    if outfile:
        if outfile.endswith('.xyz'):
            pass
        else:
            outfile = outfile + '.xyz'
        outputfile = outfile
    else:
        outputfile = os.path.splitext(logfile)[0] + '.xyz'

    if overwrite:
        # Do not change output file name.
        pass
    else:
        # Add version number without overwriting.
        counter = 0
        while os.path.exists(outputfile):
            if counter > 99:
                raise Exception("Too many versions.")
            else:
                pass
            if outfile:
                outputfile = os.path.splitext(outfile)[0] +\
                    '{}.xyz'.format(counter)
            else:
                outputfile = os.path.splitext(logfile)[0] +\
                    '{}.xyz'.format(counter)
            counter += 1

    with open(outputfile, 'w') as outf:
        if last_frame:
            save_frames = [save_frames[-1]]
        else:
            pass

        for f in save_frames:
            outf.write('{}\n'.format(len(cartcoords[f])))
            if no_energy:
                outf.write('i = {0:3d}\n'.format(step_num[f][-2]))
            else:
                outf.write('i = {0:3d}, E = {1: >17.12f}\n'
                           .format(step_num[f][-2], energy[f]))
            for line in cartcoords[f]:
                element = element_dict[int(line[0])]
                outf.write(' {0:s}\t\t\t{1: 2.5f}  {2: 2.5f}  {3: 2.5f}\n'\
                           .format(element, *line[1:]))

    print('Output written to {}'.format(outputfile))
    return


def write_to_multiple_files(outfile, save_frames, cartcoords, step_num, energy, no_energy=False):
    """
    Write to multiple files

    Paramters
    ---------
    outfile : str
	   str that contains name of file without extension
    no_energy : bool
        Enabling turns off recording of energy (Default = False)
        ** SA - DO WE NEED THIS? **
    save_frames : list
        list of frame indices
    cartcoords : array
        xyz coordinates
    step_num : list of lists
       list of lists containing current and max iterations
    energy : array
        energy of each frame

    """
    element_dict = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
                    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 16: 'S'}
    scan_energy = []

    for i, f in enumerate(save_frames):
        #Append a number to outfile to note which frame
        single_frame_file = outfile + str(i) + '.xyz'
        with open(single_frame_file, 'w') as outf:
            scan_energy.append(energy[f])
            #outf.write('{}\n'.format(len(cartcoords[f])))
            #if no_energy:
            #    outf.write('i = {0:3d}\n'.format(step_num[f][-2]))
            #else:
            #    outf.write('i = {0:3d}, E = {1: >17.12f}\n'
            #               .format(step_num[f][-2], energy[f]))
            for line in cartcoords[f]:
                element = element_dict[int(line[0])]
                outf.write(' {0:s}\t\t\t{1: 2.5f}  {2: 2.5f}  {3: 2.5f}\n'\
                           .format(element, *line[1:]))

    return scan_energy

def get_scan_coord(logfile, scan_atoms, scan_dihedral=True):
    """
    Gets values for scan

    logfile: str
        string with  full name of log file including extension
    scan_atoms : list
        list of atoms used for scan
    scan_dihedral : bool
        indicates if a dihedral scan was performed (Default: True)

    """
    if not scan_dihedral:
        raise Exception("Currently only handles scans for dihedral (scan_dihedral must = True)")
    if len(scan_atoms) is not 4:
        raise Exception("Need 4 atoms to describe a dihedral")

    pattern=f'D({scan_atoms[0]},{scan_atoms[1]},{scan_atoms[2]},{scan_atoms[3]})'

    scan_coord = []

    with open(logfile, 'r') as logf:
        for line in logf:
            if pattern in line:
                if 'Scan' not in line:
                    scan_coord.append(float(line.split()[3]))

    return scan_coord
