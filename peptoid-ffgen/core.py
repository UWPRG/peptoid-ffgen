import os
import pandas as pd
import gauss_utils as gu
import utils as u

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
    u.check_path_exists(output_dir)

    gu.check_post_hartree_fock(post_hartree_fock)

    cartcoords, step_num, energy = gu.read_log_into_lists(logfile, post_hartree_fock, no_energy)

    save_frames = gu.determine_and_save_frames(step_num, full, scan)

    if multiple_files:
        scan_energy = gu.write_to_multiple_files(output_dir+outfile, save_frames,
                                   cartcoords, step_num, energy)
    else:
        gu.write_to_file(output_dir+outfile, overwrite, last_frame, save_frames,
                         cartcoords, step_num, energy)

    scan_coord = gu.get_scan_coord(logfile, scan_atoms)

    return scan_coord, scan_energy # need this

### Convert Gaussian xyz to Gromacs pdb ###

def parse_xyz(xyz_dir, hash_table, skel_file, output_dir, file_index=1,
            file_type='pdb', multiple_files=True):
    """
    A wrapper for gen_pdb to be run on multiple xyz files that are all located
    in the same dir

    Inputs
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

    Inputs
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

    u.check_path_exists(output_dir)

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
        value = [row['atom_name'], row['residue_name']]

        # find atom_name and residue_name that map to index using a reverse
        # hash table lookup
        for k,v in hash_table.items():
            #print (k,v)
            if v == value:
                #print(k)
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
    with open(output_dir+str(file_index)+"_names"+f_type, "w") as outfile:
        outfile.write(header+write_string+footer)

    return

### Functions for running MD on GROMACS ###

def run_md(engine="GROMACS"):
    # recode wrapper script here

    # symbolic link charmm directory

    #

    return

def grompp(mpirun, mdp, coord, tpr, index="index.ndx", topol="topol.top",
           maxwarn=1, np=1, logfile="stdout.log"):
    """
    Python wrapper for gmx grompp

    Inputs
    mpirun : Bool
        Is this a multi-node run or not (gmx vs gmx_mpi), if True must specify
        number of processes (np)
    mdp : str
        Filename of .mdp file
    coord : str
        Filename of .gro or .pdb file
    tpr : str
        Filename of output of grompp
    index : str
        Filename of index file (Default index.ndx)
    topol : str
        Filename of topology file (Default topol.top)
    maxwarn : int
        Maximum number of acceptable warnings when grompping
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + np + "gmx_mpi "
    elif mpirun == False:
        mpi = "gmx "
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "grompp", "-f", mdp, "-p",
                topol, "-c", coord, "-o", tpr, "-n", index,
                "-maxwarn", maxwarn]

    subprocess.run(commands, stdout=logfile)
    return

def mdrun(mpirun, deffnm, plumed, np=1):
    """
    Python wrapper for gmx mdrun -deffnm

    Inputs
    mpirun : Bool
        Is this a multi-node run or not (gmx vs gmx_mpi), if True must specify
        number of processes (np)
    deffnm : str
         File names for md run
    plumed : str
        name of plumed file
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + np + "gmx_mpi "
    elif mpirun == False:
        mpi = "gmx "
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "mdrun", "-deffnm", deffnm, "-plumed", plumed]

    subprocess.run(commands, stdout=logfile)
    return
