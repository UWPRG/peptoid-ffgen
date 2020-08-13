import fileinput
import os
import pandas as pd
import subprocess
import shutil
from gauss_utils import *

from utils import *

##### Main Functions #####

# Parse through Gaussian log file
def parse_logfile(logfile, output_dir, scan_atoms, post_hartree_fock="",
                  scan_dihedral=True, multiple_files=True, full=False, scan=True,
                  outfile="", overwrite=True, last_frame=False, no_energy=False):
    """Parsing the logfile, wrapper function.
    Parameters
    ----------
    logfile : string
        string with  full name of log file including extension
    output_dir : str
        str that contains name of output directory
    post_hartree_fock : str
        post-HF level used (example='MP2' or 'MP3' else leave blank)
    scan_atoms : list
        list of atoms used for scan
    scan_dihedral : bool
        indicates if a dihedral scan was performed (Default: True)
    multiple_files : bool
        option to write coordinates to seperate files or single file
        (Default = True; writes each scan coordinate to separate file)
    full : bool
        denotes saving all frames (True) or optimized frames (False)
        (Default = False)
    no_energy : bool
        Enabling turns off recording of energy (Default = False)
    last_frame: bool
        denotes whether to save only the last frame (Default = False)
    overwrite : bool
        denotes whether or not output file overwrites previous file of same
        name (Default = True)
    scan : bool
        Is log file for a scan calculation (Default = True)

    Returns
    ---------
    scan_coord : list
        list of scan coordinates
    energy : list
        qm energy from run (in Hartree)
    """
    check_path_exists(output_dir)

    check_post_hartree_fock(post_hartree_fock)

    cartcoords, step_num, energy = read_log_into_lists(logfile, post_hartree_fock, no_energy)

    save_frames = determine_and_save_frames(step_num, full, scan)

    if multiple_files:
        scan_energy = write_to_multiple_files(output_dir+outfile, save_frames,
                                   cartcoords, step_num, energy)
    else:
        write_to_file(output_dir+outfile, overwrite, last_frame, save_frames,
                         cartcoords, step_num, energy)

    scan_coord = get_scan_coord(logfile, scan_atoms)

    return scan_coord, scan_energy # need this


### Convert Gaussian xyz to Gromacs pdb ###

def parse_xyz(xyz_dir, hash_table, skel_file, output_dir, file_index=1,
            file_type='pdb', multiple_files=True):
    """
    A wrapper for gen_pdb to be run on multiple xyz files that are all located
    in the same dir

    Parameters
    ----------
    xyz_dir : str
        Directory where all processed xyz files are located
    hash_table : dictionary
        2D dictionary that maps the row name in xyz file to the correct
        atom name and residue names
    skel_file : PDB file
        Contains sample PDB file with residues in correct order (coordinates
        must be present but do not have to be correct)
    output_dir : str
        Where new PDB files will be output to
    file_index : int
        Used to index file names in the case of iterating over many files
        (i.e. if editing all the xyz files in a scan). Default to 1, can pass
        variable in a for loop format
    file_type: str
        Determines type of output, only handles PDB files for now
    """
    all_files = len([name for name in os.listdir(xyz_dir) if os.path.isfile(os.path.join(xyz_dir, name))])

    for i in range(all_files):
        xyz_file=xyz_dir+str(i)+".xyz"
        gen_pdb(xyz_file, hash_table=hash_table, skel_file=skel_file,
                output_dir=output_dir, file_index=i, file_type='pdb')


def gen_pdb(xyz_file, hash_table, skel_file, output_dir, file_index,
            file_type='pdb'):
    """
    Reads in an xyz file and returns a formatted pdb file

    Parameters
    ----------
    xyz_file : xyz coordinate file
        Output of log2xyz
    hash_table : dictionary
        2D dictionary that maps the row name in xyz file to the correct
        atom name and residue names
    skel_file : PDB file
        Contains sample PDB file with residues in correct order (coordinates
        must be present but do not have to be correct)
    output_dir : str
        Where new PDB files will be output to
    file_index : int
        Used to index file names in the case of iterating over many files
        (i.e. if editing all the xyz files in a scan). Default to 1, can pass
        variable in a for loop format
    file_type: str
        Determines type of output, only handles PDB files for now
    """

    check_path_exists(output_dir)

    if file_type == 'pdb':
        f_type = '.pdb'
    else:
        raise Exception("Sorry this code only vibes w/ PDB files right now")

    # load xyz file
    xyz_names = ['element', 'x','y','z']
    df_xyz = pd.read_csv(xyz_file, delim_whitespace=True, names=xyz_names)
    #print(df_xyz)

    # load skel file
    pdb_col = ['atom','atom_number','atom_name','residue_name','residue_number',
               'x','y','z','occupancy','mass']
    df_skel = pd.read_csv(skel_file, delim_whitespace=True, skiprows=4,
                          skipfooter=2, names=pdb_col, engine='python')

    # initalize strings to write too
    header = ""
    write_string = ""
    footer = ""

    # iterate over each row in skel file
    for i, row in df_skel.iterrows():
        value = [row['atom_name'], row['residue_name']] # iterates through skel.pdb to get correct order of atoms/residues
        #print(value)

        # find atom_name and residue_name that map to index using a reverse
        # hash table lookup
        for k,v in hash_table.items():
            #print (v)
            if v == value:
                #print("v == value", v, value, k)
                parser_index = k - 1

                # write to PDB (requires perfect formatting) --------------------
                # PDBs have specific formatting and alignment rules
                # These PDBs are the bareminimum format to be read by gromacs, be aware
                # of potential issues for PDB formatting refer to :
                # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

                spacer = " "

                # for ease in writing out print statement...
                f_a = df_skel['atom'][i]
                f_a_num = df_skel['atom_number'][i]
                f_a_name = df_skel['atom_name'][i]
                f_res_name = df_skel['residue_name'][i]
                f_chainid = 1*spacer # *Use spacer for chainID column for now*
                f_res_num = df_skel['residue_number'][i]
                f_x = df_xyz['x'][parser_index] # new coord
                f_y = df_xyz['y'][parser_index] # new coord
                f_z = df_xyz['z'][parser_index] # new coord
                f_occ = df_skel['occupancy'][i]
                f_mass = df_skel['mass'][i]

                # hard coded spacers -------------
                # spacers are hard-coded if they represent a space in the PDB file
                s1=2*spacer # b/w ATOM -> atom serial num
                # b/w atom serial num and atom name
                if len(f_a_name) == 4:
                    s2=1*spacer
                else:
                    s2=2*spacer
                s3=1*spacer # alternate location indicator
                s4=1*spacer # b/w residue name and chain identifier
                s5=1*spacer # insertion of residues
                s6=3*spacer # between insertion and x coordinate

                write_string+=(f"{f_a}{s1}{f_a_num:>5}{s2}{f_a_name:>3}{s3}{f_res_name:>3}{s4}{f_chainid}{f_res_num:>4}{s5}{s6}{f_x:8.3f}{f_y:8.3f}{f_z:8.3f}{f_occ:6.2f}{f_mass:6.2f}\n")

    with open(skel_file) as r:
        h_lines = r.readlines()[0:4]
        for h in h_lines:
            header+=str(h)

    with open(skel_file) as r:
        f_lines = r.readlines()[-2:]
        for f in f_lines:
            footer+=str(f)

    #print(header+write_string+footer)

    # write to file
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir+str(file_index)+f_type, "w") as outfile:
        outfile.write(header+write_string+footer)

    return


### Functions for running MD on GROMACS ###

def run_md(scan_coord, scan_atoms, mpirun, plumed, md_dir='md/'):
    """
    Wrapper for running energy minimization for each scan coordinate. Only runs
    with Gromacs.

    Parameters
    ----------
    scan_coord : list
        list of scan coordinates
    md_dir : str
        directory name where md will run and be stored
    """
    md_energy = []

    # setup MD
    setup_runs(scan_coord)

    # run MD

    for i, s in enumerate(scan_coord):
        top_dir = os.getcwd()+"/"
        i_pdb = 'pdb/'+str(i)+'.pdb'

        # copy forcefield
        src = top_dir+"skel/"+"charmm36-nov2016.ff/"
        dst = top_dir+"md/"+str(i)+"/charmm36-nov2016.ff"
        if not os.path.exists(dst):
            shutil.copytree(src,dst)

        # cd into directory and do some md
        with cd(top_dir+"md/"+str(i)):
            pdb2gmx(mpirun, input_pdb=str(i)+'.pdb', output_gro=str(i)+'.gro')
            editconf(mpirun, input_gro=str(i)+'.gro', output_gro='box.gro')
            make_ndx(mpirun, scan_atoms, index_input='box.gro') # how do i make this w/o running MD?
            grompp(mpirun, mdp='em.mdp', coord='box.gro', tpr='em.tpr')
            #print(i, "grompp")
            mdrun(mpirun, deffnm='em', plumed=False)
            print(f'finished scan {i}')

    # grep from colvar files or run gmx_energy for md_energy
            if plumed:
                pass
                ## still need to script this lmao ** ##
                ## basically run plumed driver and extract energy corresponding to scan CV
            else:
                energy(mpirun, input_f='em.edr', output_xvg='potential.xvg')
                with open('potential.xvg') as f:
                    for line in f:
                        pass
                    last_line = line
                    #print(last_line.split())
                    md_energy.append(float(last_line.split()[1]))

    return md_energy


def setup_runs(scan_coord, md_dir='md/'):
    """
    Generates folders for md runs based on number of scan coordinates

    Parameters
    ----------
    scan_coord : list
        list of scan coordinates
    md_dir : str
        directory name where md will run and be stored
    """
    for i, s in enumerate(scan_coord):
        check_path_exists(md_dir)
        i_scan_dir = md_dir+str(i)+"/"
        check_path_exists(i_scan_dir)
        i_pdb = 'pdb/'+str(i)+'.pdb'

        shutil.copyfile(i_pdb, i_scan_dir+str(i)+'.pdb')
        shutil.copyfile('skel/em.mdp', i_scan_dir+'em.mdp')
        shutil.copyfile('skel/cat_pdb2gmx.txt', i_scan_dir+'cat_pdb2gmx.txt')
        shutil.copyfile('skel/cat_make_ndx.txt', i_scan_dir+'cat_make_ndx.txt')
        shutil.copyfile('skel/cat_energy.txt', i_scan_dir+'cat_energy.txt')

        ## ** JOSHUA REWRITE (or supplement) w/ python funcs to be more modular (vs relying on user setup) ** ##
        shutil.copyfile('skel/plumed.dat', i_scan_dir+'plumed.dat')
        plumed_find = 'XXX'
    with fileinput.FileInput(i_scan_dir+'plumed.dat', inplace=True) as f_plumed:
            for line in f_plumed:
                print(line.replace(plumed_find, str(scan_coord[i])))
        ## ** end rewrite section ** ##
    return

## Functions for calculating error ##
def calc_error(d1, d2,rss_method):
    """
    Calculates root mean square error between two data arrays

    Parameters
    ----------
    d1 : pandas dataframe column
        list of energies
    d2 : pandas dataframe column
        list of energies
    rss_method: string
        Can evaluate error in 3 ways:
        'min_to_zero' : sets min energy from d1 and d2 to zero
        'average_energy' : normalizes energies to average energy
        'KLDIV' : takes into account a weighted average
    """
    if rss_method == 'min_to_zero':
        d1 = d1 - d1.min()
        d2 = d2 - d2.min()

         # error metric
        rmse = np.sqrt(((d1 - d2)**2).mean())

    elif rss_method == 'average_energy':
        # get on same scale
        d1 = d1 - d1.min()
        d2 = d2 - d2.min()

        # normalize by average energy - I think I'm doing this wrong check w Jim
        d1 = np.mean(d1)
        d2 = np.mean(d2)

        # error metric
        rmse = np.sqrt(((d1 - d2)**2).mean())

    elif rss_method == 'KLDIV': # Kullback Liebler divergence
        # get on same scale
        d1 = d1 - d1.min()
        d2 = d2 - d2.min()

        # calculate weight
        kb = 0.008314463 # kJ/mol
        T = 300 # Kelvin
        beta = 1/(kb*T) # 2.49 kj/mol K at 300K
        weight = np.exp(-beta*d1)

        # error metric
        rmse = np.sqrt((((d1 - d2)*weight)**2).mean())

    else:
        print (f"choose from ['min_to_zero', 'average_energy', or 'KLDIV']")

    return rmse
