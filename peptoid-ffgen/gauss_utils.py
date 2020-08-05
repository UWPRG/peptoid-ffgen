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

    Parameters
    ----------
    post_hartree_fock : str
        String of post-HF level

    Returns
    -------
    using_post_HF : bool
    """
    post_HF_options = {'MP2': 'EUMP2', 'MP3': 'EUMP3'}
    if post_hartree_fock:
        assert post_hartree_fock in post_HF_options.keys(), \
            "{} not available. ".format(post_hartree_fock) +\
            "Must choose from: {}".format(post_HF_options.keys())
        return True
    else:
        return False

def read_log_into_lists(logfile, post_hartree_fock):
    """Read the log file into lists.

    Paramters
    post_hartree_fock : str
	String of post-HF level
    logfile : str
	String with  full name of log file including extension

    ---------

    Returns
     cartcoords : array
        array with modifictions (1st and 3rd columns deleted from coordinates).
        step_num : list of lists
        energy : list of lists
    ------

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
    return cartcoords, step_num

def determine_and_save_frames(scan, full, step_num):
    """Determine which frames to write to file.Save all frames, even if it is a scan calculation.
    Parameters
    scan : bool
	boolean on scan calculation (Default = True)
    full : bool
	boolean on saving all frames (Defuault = True)
    step_num : list of lists
	 
    ----------
    
    Returns
    save_frames : list
    -------

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
                                "arguments.\nIf so, add -s (or --scan) flag "
                                "to function call.")
        save_frames = range(len(step_num))
    return save_frames

def write_to_file(outfile, overwrite, last_frame, save_frames, no_energy, cartcoords, step_num, energy):
    """Write to file
    Paramters
    outfile : str
	str that contains name of file without extension 
    overwrite : bool
        denotes whether or not output file overwrites previous file of same name (Default = True)
    last_frame: bool
        denotes whether to save only the last frame (Default = False)
    no_energy : bool
        Enabling turns off recording of energy (Default = False)
    save_frames : list
    cartcoords : array
    step_num : list of lists
    energy : array
    ---------

    Returns
    -------
     
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

def write_to_multiple_files(outfile, save_frames, no_energy, cartcoords, step_num, energy):
    """Write to file
    Paramters
    ---------
    outfile : str
	str that contains name of file without extension
    no_energy : bool
        Enabling turns off recording of energy (Default = False)
    save_frames : list
    cartcoords : array
    step_num : list of lists
    energy : array


    """
    element_dict = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
                    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 16: 'S'}
    """
    if outfile:
        if outfile.endswith('.xyz'):
            pass
        else:
            outfile = outfile + '.xyz'
        outputfile = outfile
    else:
        outputfile = os.path.splitext(logfile)[0] + '.xyz'

    """


    for i, f in enumerate(save_frames):
        #Append a number to outfile to note which frame
        single_frame_file = outfile + str(i) + '.xyz'
        with open(single_frame_file, 'w') as outf:
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
    
    return

#Wrapper functiom
def parse_logfile(scan, logfile, overwrite, full, no_energy, last_frame, post_hartree_fock, outfile):

    check_post_hartree_fock(post_hartree_fock)

    cartcoords, step_num, energy = read_log_into_lists(logfile, post_hartree_fock)

    save_frames = determine_and_save_frames(scan, full, step_num)

    if multiple_files:
        write_to_multiple_files(outfile, save_frames, no_energy, cartcoords, step_num, energy)
    else: 
        write_to_file(outfile, overwrite, last_frame, save_frames, no_energy, cartcoords, step_num, energy)
