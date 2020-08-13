def grompp(mpirun, mdp, coord, tpr, index="index.ndx", topol="topol.top",
           maxwarn=1, np=1):
    """
    Python wrapper for gmx grompp

    Parameters
    ----------
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
    mpirun : bool
        Is this a multi-node run or not gmx (False) vs gmx_mpi (Default: True)
        number of processes (np)
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + "gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "grompp", "-f", mdp, "-p",
                topol, "-c", coord, "-o", tpr, "-n", index,
                "-maxwarn", str(maxwarn)]

    subprocess.run(commands)
    return


def mdrun(mpirun, deffnm, plumed=False, plumed_file='plumed.dat', np=1):
    """
    Python wrapper for gmx mdrun -deffnm

    Parameters
    ----------
    deffnm : str
         File names for md run
    mpirun : bool
        Is this a multi-node run or not gmx (False) vs gmx_mpi (Default: True)
        number of processes (np)
    plumed : bool
        Turns plumed on/off
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + " gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    if plumed:
        commands = [mpi, "mdrun", "-deffnm", deffnm, "-plumed", plumed_file]
    else:
        commands = [mpi, "mdrun", "-deffnm", deffnm]

    subprocess.run(commands)
    return

def make_ndx(mpirun, scan_atoms, index_input, np=1):
    """
    Python wrapper for gmx make_ndx

    Parameters
    ----------
    index_input : str
        file (with extension) used for gmx make_ndx
    mpirun : bool
        Is this a multi-node run or not gmx (False) vs gmx_mpi (Default: True)
        number of processes (np)
    scan_coordinates : array of ints
        Indicates atoms involved in scan. Need this for running MD
    np : int
        Number of processes to run mpirun on (Default 1 for non-mpi run)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + " gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "make_ndx", "-f", index_input]
    pipe_command=["cat","cat_make_ndx.txt"] # pick ff in working directory and TIP3P ** this should not be hardcoded lmao FIX **
    ps = subprocess.Popen(pipe_command, stdout=subprocess.PIPE)
    output = subprocess.check_output(commands, stdin=ps.stdout)
    subprocess.run(commands)

    # append scan atoms to end of file
    if len(scan_atoms) is not 4:
        raise Exception("Need 4 atoms to describe a dihedral")

    w1 = "[ SCAN ]\n"
    w2 = f'\t{scan_atoms[0]}\t{scan_atoms[1]}\t{scan_atoms[2]}\t{scan_atoms[3]}\n'
    f1 = open("index.ndx", "a")  # append mode
    f1.write(w1+w2)
    f1.close()

    return


def pdb2gmx(mpirun, input_pdb, output_gro, np=1):
    """
    Python wrapper for gmx pdb2gmx

    Parameters
    ----------
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + " gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "pdb2gmx", "-f", input_pdb, "-o", output_gro]
    pipe_command=["cat","cat_pdb2gmx.txt"] # pick ff in working directory and TIP3P ** this should not be hardcoded lmao FIX **
    ps = subprocess.Popen(pipe_command, stdout=subprocess.PIPE)
    output = subprocess.check_output(commands, stdin=ps.stdout)

    return

def editconf(mpirun, input_gro, output_gro, distance=1.3, np=1):
    """
    Python wrapper for gmx editconf

    Parameters
    ----------
    distance : int
        distance between solute and box (default 1.3)
    """

    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + " gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "editconf", "-f", input_gro, "-o", output_gro, '-c', '-d', str(distance)]
    subprocess.run(commands) # selects FF in current directory

    return


def energy(mpirun, input_f, output_xvg, np=1):
    """
    Python wrapper for gmx energy

    Parameters
    ----------

    """
    if mpirun == True:
        mpi = "mpirun " + "np " + str(np) + " gmx_mpi"
    elif mpirun == False:
        mpi = "gmx"
    else:
        print ("mpirun only takes bool as input")

    commands = [mpi, "energy", "-f", input_f, "-o", output_xvg]
    pipe_command=["cat","cat_energy.txt"] # pick ff in working directory and TIP3P ** this should not be hardcoded lmao FIX **
    ps = subprocess.Popen(pipe_command, stdout=subprocess.PIPE)
    output = subprocess.check_output(commands, stdin=ps.stdout)

    return
