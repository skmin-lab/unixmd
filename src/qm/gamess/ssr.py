from __future__ import division
from qm.gamess.gamess import GAMESS
from misc import data, call_name, au_to_A
import textwrap

class SSR(GAMESS):
    """ Class for SSR method of GAMESS

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string functional: Exchange-correlation functional information
        :param integer active_space: Active space for SSR calculation
        :param double shift: Level shifting value in REKS SCF iterations
        :param boolean l_state_interactions: Include state-interaction terms to SA-REKS
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of GAMESS, check $VERNO
    """
    def __init__(self, molecule, basis_set="6-31g*", memory="50", functional="bhhlyp", \
        active_space=2, shift=0.3, l_state_interactions=False, qm_path="./", nthreads=1, version="00"):
        # Initialize GAMESS common variables
        super(SSR, self).__init__(basis_set, memory, qm_path, nthreads, version)

        # Initialize GAMESS SSR variables
        self.functional = functional

        self.active_space = active_space
        if not (self.active_space in [2]):
            error_message = "Invalid active space for SSR!"
            error_vars = f"active_space = {self.active_space}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.active_space == 2):
            if (molecule.nst > 2):
                error_message = "SSR(3,2) calculation is not implemented in current version of pyUNIxMD!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            elif (molecule.nst == 1):
                error_message = "GAMESS does not support single-state REKS calculation!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.shift = shift
        self.l_state_interactions = l_state_interactions

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            error_message = "Only closed shell configuration supported, check charge!"
            error_vars = f"Molecule.nelec = {int(molecule.nelec)}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # SSR can compute the gradient of several states simultaneously.
        molecule.l_nacme = False
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from SSR method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        self.get_input(molecule, bo_list)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list):
        """ Generate GAMESS input files: gamess.inp.xx

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Geometry Block
        input_geom = ""
        for iat in range(molecule.nat_qm):
            if (iat > 0):
                input_geom += " " * 8
            input_geom += f"{molecule.symbols[iat]:4s}"
            sym_index = list(data.keys()).index(molecule.symbols[iat])
            input_geom += f"{sym_index:4.1f}"
            input_geom += "".join([f"{i:14.8f}" for i in molecule.pos[iat] * au_to_A])
            if (iat < molecule.nat_qm - 1):
                input_geom += "\n"

        input_data = textwrap.indent(textwrap.dedent(f"""\
        $data
        gamess
        C1
        {input_geom}
        $end
        """), " ")

        # Basis Block
        # TODO : Currently, only 6-31g* is possible
        if not (self.basis_set in ["6-31g*"]):
            error_message = "Invalid basis sets for GAMESS!"
            error_vars = f"basis_set = {self.basis_set}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        input_basis = textwrap.indent(textwrap.dedent(f"""\
        $basis
        gbasis=n31 ngauss=6 ndfunc=1
        $end
        """), " ")

        # System Block
        input_system = textwrap.indent(textwrap.dedent(f"""\
        $system
        mwords={self.memory}
        $end
        """), " ")

        # Control Block: calculate gradient option
        runtyp = "gradient"

        input_control = textwrap.indent(textwrap.dedent(f"""\
        $contrl
        scftyp=REKS runtyp={runtyp} dfttyp={self.functional}
        units=angs icharg={molecule.charge} mult=1
        $end
        """), " ")

        # Energy functional and NAC options
        if (self.l_state_interactions):
            # SSR state
            reks_type = 1
            if (self.calc_coupling):
                # SHXF, SH, Eh
                self.nac = "Yes"
            else:
                # BOMD
                self.nac = "No"
        else:
            # SA-REKS state
            reks_type = 0
            if (self.calc_coupling):
                # SHXF, SH, Eh
                error_message = "SA-REKS calculation does not provide NACVs!"
                error_vars = f"l_state_interactions = {self.l_state_interactions}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                # BOMD
                self.nac = "No"

        # First calculation
        # REKS Block: calculate gradient for target SA-REKS or SSR state
        input_reks = textwrap.indent(textwrap.dedent(f"""\
        $reks
        rexType={reks_type} rexTarget={bo_list[0] + 1} rexShift={self.shift}
        $end
        """), " ")

        # Make 'gamess.inp.1' file
        input_gamess = ""

        input_gamess += input_control
        input_gamess += input_system
        input_gamess += input_basis
        input_gamess += input_reks
        input_gamess += input_data

        # Write 'gamess.inp.1' file
        file_name = "gamess.inp.1"
        with open(file_name, "w") as f:
            f.write(input_gamess)

        # For evaluation of NACVs, two additional calculations are required
        if (self.nac == "Yes"):
            # Second calculation
            # TODO : To reduce computational cost, reading eigenvectors is essential
            # REKS Block: calculate gradient for target SA-REKS state
            input_reks = textwrap.indent(textwrap.dedent(f"""\
            $reks
            rexType=0 rexTarget={bo_list[0] + 1} rexShift={self.shift}
            $end
            """), " ")

            # Make 'gamess.inp.2' file
            input_gamess = ""

            input_gamess += input_control
            input_gamess += input_system
            input_gamess += input_basis
            input_gamess += input_reks
            input_gamess += input_data

            # Write 'gamess.inp.2' file
            file_name = "gamess.inp.2"
            with open(file_name, "w") as f:
                f.write(input_gamess)

            # Third calculation
            # TODO : To reduce computational cost, reading eigenvectors is essential
            # REKS Block: calculate gradient for averaged functional
            input_reks = textwrap.indent(textwrap.dedent(f"""\
            $reks
            rexType=0 rexTarget=0 rexShift={self.shift}
            $end
            """), " ")

            # Make 'gamess.inp.3' file
            input_gamess = ""

            input_gamess += input_control
            input_gamess += input_system
            input_gamess += input_basis
            input_gamess += input_reks
            input_gamess += input_data

            # Write 'gamess.inp.3' file
            file_name = "gamess.inp.3"
            with open(file_name, "w") as f:
                f.write(input_gamess)


