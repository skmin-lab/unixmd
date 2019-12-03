from __future__ import division
import numpy as np
import random, os, shutil
from mqc.mqc import MQC
from fileio import unixmd_init, write_md_output, write_final_xyz, typewriter
from misc import eps
from mqc.el_prop.el_propagator import *

""" THERMOSTAT NOT IMPLEMENTED YET FOR SHXF!!!!!!!!!!!!!
    DENSITY PROPAGATRO NOT IMPLEMENTED YET FOR SHXF!!!!!!!!!!!!!
"""

class Auxiliary_Molecule(object):
    """ auxiliary trajectory class
    """
    def __init__(self, molecule):
        # Initialize auxiliary molecule
        self.pos = []
        self.vel = []
        self.vel_old = []

        for ist in range(molecule.nst):
            self.pos.append(molecule.pos)
            self.vel.append(molecule.vel)
            self.vel_old.append(molecule.vel)

        self.pos = np.array(self.pos)
        self.vel = np.array(self.vel)
        self.vel_old = np.array(self.vel_old)

        self.force_hop = False


class SHXF(MQC):
    """ decoherence-indeced surface hopping based on exact factorization dynamics
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, threshold=0.01, wsigma=0.1):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_adjnac)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        self.l_hop = False

        # Initialize XF related variables
        self.l_coh = []
        self.l_first = []
        for ist in range(molecule.nst):
            self.l_coh.append(False)
            self.l_first.append(False)
        self.phase = np.array(np.zeros((molecule.nst, molecule.nat, molecule.nsp)))
        self.tot_E = np.array(np.zeros((molecule.nst)))
        self.pos_old = np.array(molecule.pos)
        self.threshold = threshold
        self.wsigma = wsigma

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(molecule)

    def run(self, molecule, theory, thermostat, input_dir="./", \
        save_QMlog=False, save_scr=True):
        """ run MQC dynamics according to decoherence-induced surface hopping dynamics
        """
        # TODO: current SHXF with thermostat Not Implemented -> should be removed
        #thermotype = thermostat.__class__.__name__
        #if (thermotype != "none"):
        #    raise ValueError ("Thermostat for SHXF Dynamics Not Implemented")
        # set directory information
        input_dir = os.path.expanduser(input_dir)
        base_dir = os.path.join(os.getcwd(), input_dir)

        unixmd_dir = os.path.join(base_dir, "md")
        if (os.path.exists(unixmd_dir)):
            shutil.rmtree(unixmd_dir)
        os.makedirs(unixmd_dir)

        QMlog_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(QMlog_dir)):
            shutil.rmtree(QMlog_dir)
        if (save_QMlog):
            os.makedirs(QMlog_dir)

        # initialize unixmd
        os.chdir(base_dir)
        bo_list = [self.rstate]
        theory.calc_coupling = True
        unixmd_init(molecule, theory.calc_coupling, self.propagation, \
            unixmd_dir, SH_chk=True)

        # calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.hop_prob(molecule, -1, unixmd_dir)
        self.hop_check(molecule, bo_list, -1, unixmd_dir)
        if (self.l_hop):
            self.evaluate_hop(molecule)

        self.update_energy(molecule)

        self.check_decoherence(molecule)
        self.check_coherence(molecule)
        self.aux_propagator(molecule)
        self.get_phase(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, \
            self.propagation, unixmd_dir)

        # main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            molecule.backup_bo()
            theory.get_bo(molecule, base_dir, istep, bo_list, calc_force_only=False)

            if (not molecule.l_nacme):
                molecule.adjust_nac()

            self.cl_update_velocity(molecule)

            if (not molecule.l_nacme):
                molecule.get_nacme()

            self.el_propagator(molecule)

            self.hop_prob(molecule, istep, unixmd_dir)
            self.hop_check(molecule, bo_list, istep, unixmd_dir)
            if (self.l_hop):
                self.evaluate_hop(molecule)
                if (theory.re_calc):
                    theory.get_bo(molecule, base_dir, istep, bo_list, calc_force_only=True)

            thermostat.run(molecule, self)

            self.update_energy(molecule)

            self.check_decoherence(molecule)
            self.check_coherence(molecule)
            self.aux_propagator(molecule)
            self.get_phase(molecule)

            write_md_output(molecule, theory.calc_coupling, istep, \
                self.propagation, unixmd_dir)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, istep, unixmd_dir)

        # delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(base_dir, "md/scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def hop_prob(self, molecule, istep, unixmd_dir):
        """ Routine to calculate hopping probabilities
        """
        # reset surface hopping variables
        self.l_hop = False
        self.force_hop = False
        self.rstate_old = self.rstate
        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        accum = 0.

        if (molecule.rho.real[self.rstate, self.rstate] < eps):
            self.force_hop = True

        for ist in range(molecule.nst):
            if (ist != self.rstate):
                if (self.force_hop):
                    self.prob[ist] = molecule.rho.real[ist, ist] / (1. - eps) 
                else: 
                    self.prob[ist] = - 2. * molecule.rho.real[ist, self.rstate] * \
                        molecule.nacme[ist, self.rstate] * self.dt / molecule.rho.real[self.rstate, self.rstate]

                if (self.prob[ist] < 0.):
                    self.prob[ist] = 0.
                accum += self.prob[ist]
            self.acc_prob[ist + 1] = accum
        psum = self.acc_prob[molecule.nst]
 
        if (psum > 1.):
            self.prob /= psum
            self.acc_prob /= psum

        # write SHPROB file
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(molecule.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB")

    def hop_check(self, molecule, bo_list, istep, unixmd_dir):
        """ Routine to check hopping occurs with random number
        """
        rand = random.random()
        for ist in range(molecule.nst):
            if (ist == self.rstate):
                continue
            if (rand > self.acc_prob[ist] and rand <= self.acc_prob[ist + 1]):
                self.l_hop = True
                self.rstate = ist
                bo_list[0] = self.rstate

        # write SHSTATE file
        tmp = f'{istep + 1:9d}{"":14s}{self.rstate}'
        typewriter(tmp, unixmd_dir, "SHSTATE")

    def evaluate_hop(self, molecule):
        """ Routine to evaluate hopping and velocity rescaling
        """
        pot_diff = molecule.states[self.rstate].energy - molecule.states[self.rstate_old].energy
        if (molecule.ekin < pot_diff):
            if (not self.force_hop):
                self.l_hop = False
                self.rstate = self.rstate_old
        else:
            if (molecule.ekin < eps):
                raise ValueError ("Too small kinetic energy!")
            fac = 1. - pot_diff / molecule.ekin
            molecule.vel *= np.sqrt(fac)
            # update kinetic energy
            molecule.update_kinetic()

    def calculate_force(self, molecule):
        """ Routine to calculate the forces
        """
        self.rforce = np.copy(molecule.states[self.rstate].force)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in surface hopping
        """
        # update kinetic energy
        molecule.update_kinetic()
        molecule.epot = molecule.states[self.rstate].energy
        molecule.etot = molecule.epot + molecule.ekin

    def check_coherence(self, molecule):
        """ Check coherence and reset density
        """
        count = 0
        for ist in range(molecule.nst):
            rho = molecule.rho.real[ist, ist]
            if (rho > self.upper_th or rho < self.lower_th):
                self.l_coh[ist] = False
            else:
                if (self.l_coh[ist]):
                    self.l_first[ist] = False
                else:
                    self.l_first[ist] = True
                self.l_coh[ist] = True
                count += 1

        if (count < 2):
            for ist in range(molecule.nst):
                self.l_coh[ist] = False
                self.l_first[ist] = False

#        for ist in range(molecule.nst):
#            if (self.l_coh[ist]):
#                print(f" TSHXF DECO T           {ist + 1}", flush=True)
#            else:
#                print(f" TSHXF DECO F           {ist + 1}", flush=True)

    def check_decoherence(self, molecule):
        """ Check coherence and reset density
        """
        if (self.l_hop):
            for ist in range(molecule.nst):
                self.l_coh[ist] = False
                self.l_first[ist] = False
        else:
            for ist in range(molecule.nst):
                if (self.l_coh[ist]):
                    rho = molecule.rho.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(molecule, ist)
                        return

    def aux_propagator(self, molecule):
        """ Update auxiliary positions/velocities
        """
        self.pos_old = np.copy(molecule.pos)
        # Get auxiliary position
        for ist in range(molecule.nst):
            if (self.l_coh[ist] and ist != self.rstate and not self.l_first[ist]):
                self.aux.pos[ist] += self.aux.vel[ist] * self.dt
            else:
                self.aux.pos[ist] = molecule.pos
#            print(f"AUX_POS {self.aux.pos[ist, 0, 0]:15.8f}            {ist + 1}", flush=True)

        # Get auxiliary velocity
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(molecule.nst):
            if (self.l_coh[ist] and ist != self.rstate):
                if (self.l_first[ist]):
                    self.tot_E[ist] = molecule.ekin + molecule.states[ist].energy
                    alpha = 1.
                    self.aux.vel[ist] = molecule.vel
                else:
                    alpha = self.tot_E[ist] - molecule.states[ist].energy
                    if (alpha < eps):
                        alpha = 0.
                    else:
                        alpha /= molecule.ekin
                        alpha = np.sqrt(alpha)
                    self.aux.vel[ist] = molecule.vel * alpha
            else:
                self.aux.vel[ist] = molecule.vel
#            print(f"AUX_VEL {self.aux.vel[ist, 0, 0]:15.8f}            {ist + 1}", flush=True)

    def set_decoherence(self, molecule, one_st):
        """ Reset densities in case of molecule got decohered
        """
        self.phase = np.zeros((molecule.nst, molecule.nat, molecule.nsp))
        for ist in range(molecule.nst):
            self.l_coh[ist] = False
            self.l_first[ist] = False
            if (ist == one_st):
                molecule.states[ist].coef /= np.sqrt(molecule.rho[ist, ist].real)
                molecule.rho[ist, ist] = 1. + 0.j
            else:
                molecule.states[ist].coef = 0. + 0.j
 
    def get_phase(self, molecule):
        """ Phase calculation routine 
        """
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.phase[ist] = 0.
                else:
                    for iat in range(molecule.nat):
                        for isp in range(molecule.nsp):
                            self.phase[ist, iat, isp] += molecule.mass[iat] * \
                            (self.aux.vel[ist, iat, isp] - self.aux.vel_old[ist, iat, isp])
#            print(f"NABPH {self.phase[ist, 0, 0]:15.8f}            {ist + 1}", flush=True)

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix
        """
        if (self.propagation == "coefficient"):
            el_coef_xf(self, molecule)
        elif (self.propagation == "density"):
            raise ValueError ("Density Propagation Not Implemented")
        else:
            raise ValueError ("Other Propagators Not Implemented")


