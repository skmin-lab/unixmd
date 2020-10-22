

     ******************************************
     **    PROGRAM:              MCPC        **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


 original author: Daniel Robertson, FSU
 later revisions: Ron Shepard, ANL;
                  Michal Dallos, University Vienna



 This Version of Program mcpc is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



   ******  File header section  ******

 Headers form the restart file:
    Hermit Integral Program : SIFS version  odin1             11:04:44.739 15-Oct-20
     title                                                                          


   ******  DRT info section  ******

 Informations for the DRT no.  1
 Header form the DRT file: 
     title                                                                          
 Molecular symmetry group:   sym1
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    4
 Number of active electrons:   6
 Total number of CSFs:        10

   ***  Informations from the DRT number:   1

 
 Symmetry orbital summary:
 Symm.blocks:         1
 Symm.labels:         a  

 List of doubly occupied orbitals:
  1 a    2 a    3 a    4 a    5 a  

 List of active orbitals:
  6 a    7 a    8 a    9 a  


   ******  MCSCF convergence information:  ******

 MCSCF convergence criteria were satisfied.

 mcscf energy=   -94.1577090136    nuclear repulsion=    38.8015478431
 demc=             0.0000000003    wnorm=                 0.0000007657
 knorm=            0.0000000370    apxde=                 0.0000000000


 MCSCF calculation performmed for   1 symmetry.

 State averaging:
 No,  ssym, navst, wavst
  1    a      3   0.3333 0.3333 0.3333

 Input the DRT No of interest: [  1]:
In the DRT No.: 1 there are  3 states.

 Which one to take? [  1]:
 The CSFs for the state No  2 of the symmetry  a   will be printed
 according to the following print options :

 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: csfs will be printed based on coefficient magnitudes.

 input the coefficient threshold (end with 0.) [ 0.0000]:
 List of active orbitals:
  6 a    7 a    8 a    9 a  

   csf       coeff       coeff**2    step(*)
  -----  ------------  ------------  ------------
      2 -0.6796301197  0.4618970995  3312
      4  0.6627638381  0.4392559051  3132
      5 -0.2066105185  0.0426879064  3123
      1 -0.2000325413  0.0400130176  3330
      8 -0.0961927216  0.0092530397  1323
      7 -0.0491168555  0.0024124655  1332
      6 -0.0485829642  0.0023603044  3033
      3  0.0345460339  0.0011934285  3303
     10  0.0247759328  0.0006138468  0333
      9 -0.0176914233  0.0003129865  1233

 input the coefficient threshold (end with 0.) [ 0.0000]:
 1) print csf info by sorted index number.
 2) print csf info by contribution threshold.
 3) print csf info by csf number.
 4) set additional print options.
 5) print the entire sorted csf vector.
 6) print the entire csf vector.
 7) print the mcscf molecular orbitals.
 8) print the mcscf natural orbitals and occupation numbers.
 9) export wave function files for cioverlap (all states).
 0) end.

 input menu number [  0]: