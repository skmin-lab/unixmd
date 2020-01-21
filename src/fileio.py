from __future__ import division
from misc import au_to_A
import os, textwrap, datetime

def touch_file(molecule, calc_coupling, propagation, unixmd_dir, SH_chk):
    """ Open the UNIXMD output files
    """
    write_init_base(molecule, unixmd_dir)
    if (calc_coupling):
        write_init_coupling(molecule, propagation, unixmd_dir, SH_chk)

    # Print UNIXMD version
    cur_time = datetime.datetime.now()
    cur_time = cur_time.strftime("%Y-%m-%d %H:%M:%S")
    prog_info = textwrap.dedent(f"""\
    {"-" * 68}

    {"UNI-XMD version 20.1":>43s}

    {"< Developers >":>40s}
    {" " * 4}Seung Kyu Min,  In Seong Lee,  Jong-Kwon Ha,  Daeho Han,
    {" " * 4}Kicheol Kim,  Tae-In Kim,  Sung Wook Moon

    {"-" * 68}

    {" " * 4}Please cite UNI-XMD as follows:
    {" " * 4}This is article

    {" " * 4}UNI-XMD begins on {cur_time}
    """)
    print (prog_info, flush=True)

def write_init_base(molecule, unixmd_dir):
    """ Header for non-coupling output files
    """
    # energy information file header
    tmp = f'{"#":5s}{"Step":9s}{"Kinetic(H)":15s}{"Potential(H)":15s}{"Total(H)":15s}' + \
        "".join([f'E({ist})(H){"":8s}' for ist in range(molecule.nst)])
    typewriter(tmp, unixmd_dir, "MDENERGY")

def write_init_coupling(molecule, propagation, unixmd_dir, SH_chk):
    """ Header for coupling output files
    """
    # surface hopping running state and hopping probabilities file header
    if (SH_chk):
        tmp = f'{"#":5s}{"Step":8s}{"Running State":10s}'
        typewriter(tmp, unixmd_dir, "SHSTATE")

        tmp = f'{"#":5s}{"Step":12s}' + "".join([f'Prob({ist}){"":8s}' for ist in range(molecule.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB")

    # BO coefficents, densities file header
    if (propagation == "density"):
        tmp = f'{"#":5s} Density Matrix: population Re; see the manual for detail orders'
        typewriter(tmp, unixmd_dir, "BOPOP")
        tmp = f'{"#":5s} Density Matrix: coherence Re-Im; see the manual for detail orders'
        typewriter(tmp, unixmd_dir, "BOCOH")
    elif (propagation == "coefficient"):
        tmp = f'{"#":5s} BO State Coefficients: state Re-Im; see the manual for detail orders'
        typewriter(tmp, unixmd_dir, "BOCOEF")
    else:
        raise ValueError (f"Unknown propagation type! (current input: {propagation})")

    # NACME file header
    tmp = f'{"#":5s}Non-Adiabatic Coupling Matrix Eliments: off-diagonal'
    typewriter(tmp, unixmd_dir, "NACME")

def write_md_output(molecule, calc_coupling, istep, propagation, unixmd_dir):
    """ Writting output files
    """
    # write MOVIE.xyz file including positions and velocities
    tmp = f'{molecule.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)'
    typewriter(tmp, unixmd_dir, "MOVIE.xyz")
    for iat in range(molecule.nat):
        tmp = f'{molecule.symbols[iat]:5s}' + \
            "".join([f'{molecule.pos[iat, isp] * au_to_A:15.8f}' for isp in range(molecule.nsp)]) \
            + "".join([f"{molecule.vel[iat, isp]:15.8f}" for isp in range(molecule.nsp)])
        typewriter(tmp, unixmd_dir, "MOVIE.xyz")

    # write MDENERGY file including several energy information
    tmp = f'{istep + 1:9d}{molecule.ekin:15.8f}{molecule.epot:15.8f}{molecule.etot:15.8f}' \
        + "".join([f'{states.energy:15.8f}' for states in molecule.states])
    typewriter(tmp, unixmd_dir, "MDENERGY")

    if (calc_coupling):
        # write BOCOEF, BOPOP, BOCOH files
        if (propagation == "density"):
            tmp = f'{istep + 1:9d}' + "".join([f'{molecule.rho.real[ist, ist]:15.8f}' for ist in range(molecule.nst)])
            typewriter(tmp, unixmd_dir, "BOPOP")
            tmp = f'{istep + 1:9d}' + "".join([f"{molecule.rho.real[ist, jst]:15.8f}{molecule.rho.imag[ist, jst]:15.8f}" \
                for ist in range(molecule.nst) for jst in range(ist + 1, molecule.nst)])
            typewriter(tmp, unixmd_dir, "BOCOH")
        elif (propagation == "coefficient"):
            tmp = f'{istep + 1:9d}' + "".join([f'{states.coef.real:15.8f}{states.coef.imag:15.8f}' \
                for states in molecule.states])
            typewriter(tmp, unixmd_dir, "BOCOEF")

        # write NACME file
        tmp = f'{istep + 1:10d}' + "".join([f'{molecule.nacme[ist, jst]:15.8f}' \
            for ist in range(molecule.nst) for jst in range(ist + 1, molecule.nst)])
        typewriter(tmp, unixmd_dir, "NACME")

def write_final_xyz(molecule, istep, unixmd_dir):
    """ Write final positions and velocities
    """
    # write FINAL.xyz file including positions and velocities
    tmp = f'{molecule.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)'
    typewriter(tmp, unixmd_dir, "FINAL.xyz")
    for iat in range(molecule.nat):
        tmp = f'{molecule.symbols[iat]:5s}' + \
            "".join([f'{molecule.pos[iat, isp] * au_to_A:15.8f}' for isp in range(molecule.nsp)]) \
            + "".join([f"{molecule.vel[iat, isp]:15.8f}" for isp in range(molecule.nsp)])
        typewriter(tmp, unixmd_dir, "FINAL.xyz")

def typewriter(string, dir_name, filename):
    """ funtion to write any string in dir_name/filename
    """
    tmp_name = os.path.join(dir_name, filename)
    with open(tmp_name, "a") as f:
        f.write(string + "\n")



