import os.path
import re
import sys
import argparse
import numpy as np

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
