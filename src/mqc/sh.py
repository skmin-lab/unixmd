from __future__ import division
import numpy as np
import random, os, shutil
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz, typewriter
from misc import eps, au_to_K
from mqc.el_prop.el_propagator import *
 
class SH(MQC):    
    """ surface hopping dynamics
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_adjnac)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate  
        
        self.rand = 0.0
        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        self.l_hop = False
        self.force_hop = False

    def run(self, molecule, theory, thermostat, input_dir="./", \
        save_QMlog=False, save_scr=True):
        """ run MQC dynamics according to surface hopping dynamics
        """
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

        touch_file(molecule, theory.calc_coupling, self.propagation, \
            unixmd_dir, SH_chk=True)
        self.print_init(molecule, theory, thermostat)

        # calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.hop_prob(molecule, -1, unixmd_dir)
        self.hop_check(molecule, bo_list, -1, unixmd_dir)
        if (self.l_hop):
            self.evaluate_hop(molecule)

        self.update_energy(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, \
            self.propagation, unixmd_dir)

        self.print_step(molecule, -1)

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

            self.print_step(molecule, istep)
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
        self.rand = random.random()
        for ist in range(molecule.nst):
            if (ist == self.rstate):
                continue
            if (self.rand > self.acc_prob[ist] and self.rand <= self.acc_prob[ist + 1]):
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

#        a = 0.
#        b = 0.
#        for iat in range(molecule.nat):
#            for isp in range(molecule.nsp):
#                a += 0.5 * molecule.mass[iat] * molecule.nac[rstate_old, self.rstate, iat, isp] * molecule.nac[rstate_old, self.rstate, iat, isp]
#                b += molecule.mass[iat] * molecule.nac[rstate_old, self.rstate, iat, isp] * molecule.vel[iat, isp]
#
#        pot_diff = molecule.states[self.rstate].energy - molecule.states[rstate_old].energy
#        det = b ** 2. - 4. * a * pot_diff
#
#        if (det < 0.):
#            self.l_hop = False
#            self.rstate = self.rstate_old
#        else:
#            if(b < 0.):
#                x = (- b - np.sqrt(det)) / a * 0.5
#            else:
#                x = (- b + np.sqrt(det)) / a * 0.5
#            for iat in range(molecule.nat):
#                molecule.vel[iat, :] += x * molecule.nac[rstate_old, self.rstate, iat, :]

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

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix
        """
        if (self.propagation == "coefficient"):
            el_coef(self.nesteps, self.dt, molecule)
        elif (self.propagation == "density"):
            el_rho(self.nesteps, self.dt, molecule)
        else:
            raise ValueError ("Other propagators Not Implemented")

    #def print_step(self, molecule, istep, debug=0):
    def print_step(self, molecule, istep):
        """ Routine to print each steps infomation about dynamics
        """
        if(istep == -1):
            INFO = f"INFO {istep+1:7d}{self.rstate:9d}{0.0:14.8f}({self.rstate}->{self.rstate}){self.rand:15.8f}"
        else:
            max_prob = max(self.prob)
            hstate = np.where(self.prob == max_prob)[0][0]
            INFO = f"INFO {istep+1:7d}{self.rstate:9d}{max_prob:14.8f}({self.rstate}->{hstate}){self.rand:15.8f}"
        
        INFO += f"{molecule.ekin:16.8f}{molecule.epot:16.8f}{molecule.etot:16.8f}"

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        INFO += f"{ctemp:18.6f}"
        
        norm = 0.0
        for ist in range(molecule.nst):
            norm += molecule.rho.real[ist,ist]
        INFO += f"{norm:18.8f}"
        print(INFO, flush=True)

        if(self.rstate != self.rstate_old):
            print(f"Hopping {self.rstate_old} -> {self.rstate}", flush=True)
        
        if(self.force_hop):
            print(f"Force hop {self.rstate_old} -> {self.rstate}", flush=True)
        
        #if(debug=1):
        debug1 = f"DEBUG1 {istep+1:5d}"
        for ist in range(molecule.nst):
            debug1 += f"{molecule.states[ist].energy:15.8f}"
        print(debug1, flush=True)

        debug2 = f"DEBUG2 {istep+1:5d}"
        for ist in range(molecule.nst): 
          debug2 += f"{self.prob[ist]:14.8f}({self.rstate}->{ist})"
        print(debug2, flush=True)
