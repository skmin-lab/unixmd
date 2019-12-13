from __future__ import division
import os, shutil, re
from bo.columbus.columbus import columbus
import numpy as np
import textwrap
from bo.columbus.colbasis import *
from misc import data, au_to_A, A_to_au, amu_to_au

class CASSCF(columbus):
    """ Class for Columbus CASSCF method
                 bomd | sh | eh | nac | re_calc
        CASSCF :  o     o    o     T      T
    """
    def __init__(self, molecule, basis_set="6-31g*", memory="500", \
        active_elec=2, active_orb=2, qm_path="/opt/Columbus7.0/Columbus", nthreads=1, version=7.0):
        # Initialize Columbus common variables
        super().__init__(basis_set, memory, qm_path, nthreads, version)


        # TODO : periodic option?
#        self.max_iter = max_iter
#        self.scf_en_tol = scf_en_tol
#        self.scf_grad_tol = scf_grad_tol
#        self.scf_step_tol = scf_step_tol
        self.active_elec = active_elec
        self.active_orb = active_orb
#        self.cpscf_grad_tol = cpscf_grad_tol

        # calculate number of frozen, closed and occ orbitals in CASSCF method
        # no positive frozen core orbitals in CASSCF
        self.frozen_orb = 0
        self.closed_orb = int((int(molecule.nelec) - self.active_elec) / 2)
        self.occ_orb = int(self.closed_orb + self.active_orb)

        os.environ["COLUMBUS"] = qm_path
        # set 'l_nacme' with respect to the computational method
        # CASSCF can produce NACs, so we do not need to get NACME from CIoverlap
        # CASSCF can compute the gradient of several states simultaneously,
        #        but self.re_calc is set to be true to reduce cost.
        molecule.l_nacme = False
        self.re_calc = True

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from Columbus
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, restart):
        """ Generate Columbus input files: base reference file need for qm in base_directory/bo/columbus/ref_col
        """
        #generate geom
        geom = ""
        for iat in range(molecule.nat):
            anum =   list(data.keys()).index(f"{molecule.symbols[iat]}")
            tmp = f' {molecule.symbols[iat]:5s}{anum:7.2f}' + "".join([f'{molecule.pos[iat, isp]:15.8f}' for isp in range(molecule.nsp)])  \
             + f'{molecule.mass[iat]/amu_to_au:15.8f}' + "\n" 
            geom += tmp
   
        # write geom file
        file_name = "geom"
        with open(file_name, "w") as f:
            f.write(geom)

        if (molecule.nelec % 2 == 0):
            closeshell = "yes"
            DOCC1 = int(int(molecule.nelec) / 2)
        else:
            closeshell = "no"
            DOCC1 = int((int(molecule.nelec) - 1 ) / 2  )

        nelec = int(molecule.nelec) 
   
        stdin = f"\ny\n1\nn\nno\n2\n{closeshell}\n{DOCC1}\nyes\nno\n\n"
        stdin += f"3\nn\n3\n1\n{nelec}\n1\n1\n0\n0\n{self.closed_orb}\n{self.active_orb}\nn\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        stdin += f"5\n1\n1\n2\n3\n11\n1\ny\n\n3\ny\n8\n4\n7\n\n"
        
        f = open("geom", "r")
        geomsym = f.readlines()
        f.close
        
        asymbol = []
        for i in range(len(geomsym)):
          tmp = geomsym[i].split()
          if tmp[0] in asymbol:
            continue
          else:
            asymbol.append(tmp[0])
        asymbol.sort()
        
        shutil.copy(os.path.join(self.qm_path, "prepinp"), "prepinp_copy")
        
        f = open("prepinp_copy", "r")
        prepinp = f.readlines()
        f.close
        
        indexnu = prepinp.index("  foreach $el ( keys %sumformula ){\n")
        prepinp[indexnu]  = "  foreach $el (sort keys %sumformula ){\n"
        
        with open("prepinp_fix", "w") as f:
            for i in range(len(prepinp)):
                f.write(prepinp[i])
        os.system("chmod 755 prepinp_fix")
        
        basisnu = [1] * len(asymbol)
        
        if (self.basis_set == "cc-pvdz"):
            for i in range(len(asymbol)):
                basisnu[i] = cc_pvdz[asymbol[i]]
        elif (self.basis_set == "cc-pvtz"):
            for i in range(len(asymbol)):
                basisnu[i] = cc_pvtz[asymbol[i]]
        elif (self.basis_set == "cc-pvqz"):
            for i in range(len(asymbol)):
                basisnu[i] = cc_pvqz[asymbol[i]]
        elif (self.basis_set == "3-21g*"):
            for i in range(len(asymbol)):
                basisnu[i] = t_21gs[asymbol[i]]
        elif (self.basis_set == "3-21+g*"):
            for i in range(len(asymbol)):
                basisnu[i] = t_21pgs[asymbol[i]]
        elif (self.basis_set == "6-31g"):
            for i in range(len(asymbol)):
                basisnu[i] = s_31g[asymbol[i]]
        elif (self.basis_set == "6-31g*"):
            for i in range(len(asymbol)):
                basisnu[i] = s_31gs[asymbol[i]]
        elif (self.basis_set == "6-31+g*"):
            for i in range(len(asymbol)):
                basisnu[i] = s_31pgs[asymbol[i]]
        elif (self.basis_set == "6-311g*"):
            for i in range(len(asymbol)):
                basisnu[i] = s_311gs[asymbol[i]]
        elif (self.basis_set == "6-311g+*"):
            for i in range(len(asymbol)):
                basisnu[i] = s_311pgs[asymbol[i]]
        else:
            raise ValueError ("Basis set not yet implemented in Columbus input: add manually (colbasis.py)")
        
        prepin = f"1\nc1\ngeom\n\n"
        for i in range(len(asymbol)):
            if (basisnu[i] == 0):
                raise ValueError (f"Data not found in Columbus : Atom {asymbol[i]} with  basis {basis}")
            prepin += f"{basisnu[i]}\n"
        prepin += "y\n\n" 
        
        with open("prepin", "w") as f:
            f.write(prepin)
        
        os.system("./prepinp_fix < prepin")
        
        with open("stdin", "w") as f:
            f.write(stdin)
        os.system(f"{self.qm_path}/colinp < stdin")
       
        f = open("mcscfin", "r")
        mcscfin = f.readlines()
        f.close
        
        mcs_length = len(mcscfin)
        target_line = mcs_length - 3
        new_mcs = ""

        for i in range (target_line):
            new_mcs += mcscfin[i]
        new_mcs += f"  NAVST(1) = {molecule.nst},\n"
        for i in range (molecule.nst):
            new_mcs += f"  WAVST(1,{i+1})=1 ,\n"
        new_mcs += " &end\n"
        
        os.rename("mcscfin", "mcscfin.old")
        
        with open("mcscfin", "w") as f:
            f.write(new_mcs)
        
        transmomin = "MCSCF\n"
        for i in range (molecule.nst):
            for j in range (i):
                transmomin += f"1  {i+1}  1  {j+1}\n"
        with open("transmomin", "w") as f:
            f.write(transmomin)
        
        shutil.copy("daltcomm", "daltcomm.new")


    def run_QM(self, base_dir, istep, bo_list):
        """ run columbus calculation and save the output files
        """
        # run dynamics
        command = f"{self.qm_path}/runc -m {self.memory} > runls"
        os.system(command)
          
        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("runls", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, restart):
        """ read the output files to get BO data
        """
        # read 'log' file
        file_name = "LISTINGS/mcscfsm.sp"
        with open(file_name, "r") as f:
            log_out = f.read()

        # energy
        if (not restart):
            for states in molecule.states:
                states.energy = 0.

            tmp_e = 'total\senergy[=]\s*[-]\d+[.]\d+'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
#            energy = energy.astype(float)

            for ist in range(molecule.nst):
                tmp = energy[ist].split()
                molecule.states[ist].energy = float(tmp[2])

        # force
        if (not restart):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:

            os.system(f"sed s/D/E/g GRADIENTS/cartgrd.drt1.state{ist+1}.sp > tmp")
            with open("tmp", "r") as f:
                log_out = f.read() 

            tmp_f ='\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (not restart and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):

                        file_name = f"GRADIENTS/cartgrd.nad.drt1.state{jst+1}.drt1.state{ist+1}.sp"
                        with open(file_name, "r") as f:
                            log_out = f.read()

                        tmp_c =  '\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[0])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, 3, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]


