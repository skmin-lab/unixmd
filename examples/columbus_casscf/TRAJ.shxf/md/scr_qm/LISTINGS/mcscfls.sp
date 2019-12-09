

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
        65536000 of real*8 words (  500.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-4,
   tol(2)=1.e-4,
   tol(1)=1.e-8,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2,30
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,8,20,1,9,20
   NAVST(1) = 2,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /home/swmoon/work/unixmd-fixcol/src/TRAJ.shxf/md/scr_qm/W
 ORK

 Integral file header information:
 Hermit Integral Program : SIFS version  odin              17:46:14.039 09-Dec-19

 Core type energy values:
 energy( 1)=  3.278508429890E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   32.785084299


   ******  Basis set information:  ******

 Number of irreps:                  1
 Total number of basis functions:  36

 irrep no.              1
 irrep label           A  
 no. of bas.fcions.    36


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-04. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-04. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  1 DRT.

 DRT  first state   no.of aver.states   weights
  1   ground state          2             0.500 0.500

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        3

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       8(  8)    20
       1       9(  9)    20

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).
 30:  Compute mcscf (transition) density matrices


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a  
 Total number of electrons:   16
 Spin multiplicity:            1
 Number of active orbitals:    2
 Number of active electrons:   2
 Total number of CSFs:         3
 

 faar:   0 active-active rotations allowed out of:   1 possible.


 Number of active-double rotations:        14
 Number of active-active rotations:         0
 Number of double-virtual rotations:      189
 Number of active-virtual rotations:       54
 lenbfsdef=                131071  lenbfs=                   729
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #           6
 class  2 (pq|ri):         #          42
 class  3 (pq|ia):         #         567
 class  4 (pi|qa):         #         756
 class  5 (pq|ra):         #         162
 class  6 (pq|ij)/(pi|qj): #         231
 class  7 (pq|ab):         #        1134
 class  8 (pa|qb):         #        2187
 class  9 p(bp,ai)         #       10206
 class 10p(ai,jp):        #        2646
 class 11p(ai,bj):        #       20412

 Size of orbital-Hessian matrix B:                    36354
 Size of the orbital-state Hessian matrix C:           1542
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:          37896


 Source of the initial MO coeficients:

 Input MO coefficient file: /home/swmoon/work/unixmd-fixcol/src/TRAJ.shxf/md/scr_qm/WORK
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     2
 ciiter=   3 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8329489719     -110.6180332708        0.0000000000        0.0000010000
    2       -77.7057488483     -110.4908331472        0.0000000000        0.0000010000
    3       -77.6152113695     -110.4002956684        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  5.409698176517922E-002
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99297611 pnorm= 0.0000E+00 rznorm= 1.8940E-06 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.542173  -22.457267   -2.065445   -1.608897   -1.176766   -1.153595   -1.090895

 qvv(*) eigenvalues. symmetry block  1
     0.502739    0.590323    0.625226    0.706336    0.890479    1.366406    1.508853    1.608801    1.651956    1.749155
     1.826224    2.182658    2.270695    2.300338    2.453568    2.507621    2.931790    3.496664    3.853062    4.004054
     4.210634    4.646536    4.915389    5.138526    5.262747    5.628552    5.977105

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -77.7693489101 demc= 7.7769E+01 wnorm= 4.3278E-01 knorm= 1.1832E-01 apxde= 1.6953E-02    *not conv.*     

               starting mcscf iteration...   2

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8841966920     -110.6692809909        0.0000000000        0.0000100000
    2       -77.6938263447     -110.4789106436        0.0000000000        0.0000100000
    3       -77.6865604989     -110.4716447978        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.205053661879966E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99061050 pnorm= 0.0000E+00 rznorm= 5.0263E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.737604  -22.198458   -2.070572   -1.621803   -1.277775   -1.095715   -1.015544

 qvv(*) eigenvalues. symmetry block  1
     0.489603    0.569592    0.661055    0.692516    0.899007    1.364260    1.524066    1.602963    1.653955    1.772506
     1.829885    2.149709    2.261190    2.337775    2.380567    2.589372    2.985287    3.509723    3.763894    4.094041
     4.280704    4.631191    4.839065    5.201586    5.306065    5.644168    6.015564

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -77.7890115183 demc= 1.9663E-02 wnorm= 5.7640E-01 knorm= 1.3671E-01 apxde= 2.7596E-02    *not conv.*     

               starting mcscf iteration...   3

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8881577607     -110.6732420596        0.0000000000        0.0000100000
    2       -77.7996343822     -110.5847186811        0.0000000000        0.0000100000
    3       -77.6434864794     -110.4285707783        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.836850088361920E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99457483 pnorm= 0.0000E+00 rznorm= 5.1950E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.709898  -22.240713   -2.080452   -1.635506   -1.296127   -1.107871   -1.018899

 qvv(*) eigenvalues. symmetry block  1
     0.486832    0.562410    0.657113    0.698269    0.906173    1.361001    1.527409    1.604765    1.659132    1.773498
     1.834063    2.150271    2.262232    2.330771    2.377349    2.580681    2.981805    3.508077    3.759054    4.095539
     4.288625    4.626672    4.839273    5.201643    5.303756    5.647944    6.017875

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -77.8438960715 demc= 5.4885E-02 wnorm= 3.0695E-01 knorm= 1.0402E-01 apxde= 9.0153E-03    *not conv.*     

               starting mcscf iteration...   4

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8867695464     -110.6718538453        0.0000000000        0.0000100000
    2       -77.8197323172     -110.6048166161        0.0000000000        0.0000100000
    3       -77.5599775201     -110.3450618190        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.349154381462732E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99985912 pnorm= 0.0000E+00 rznorm= 9.3747E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.589070  -22.366437   -2.073180   -1.630067   -1.271007   -1.120613   -1.039268

 qvv(*) eigenvalues. symmetry block  1
     0.501936    0.567665    0.663094    0.689088    0.907721    1.364384    1.543365    1.613519    1.642769    1.761220
     1.824014    2.170644    2.276531    2.338345    2.380765    2.555271    2.971975    3.510078    3.800180    4.058387
     4.269691    4.637828    4.867303    5.184982    5.283989    5.644906    6.006641

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -77.8532509318 demc= 9.3549E-03 wnorm= 1.0793E-02 knorm= 1.6785E-02 apxde= 6.0629E-05    *not conv.*     

               starting mcscf iteration...   5

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8868174854     -110.6719017843        0.0000000000        0.0000010000
    2       -77.8198065762     -110.6048908751        0.0000000000        0.0000010000
    3       -77.5528840254     -110.3379683243        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.946969730770716E-005
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99999994 pnorm= 0.0000E+00 rznorm= 8.0599E-07 rpnorm= 0.0000E+00 noldr=  5 nnewr=  5 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.585182  -22.369829   -2.072595   -1.629469   -1.268992   -1.120784   -1.039970

 qvv(*) eigenvalues. symmetry block  1
     0.502837    0.567985    0.663410    0.688827    0.907713    1.364651    1.543374    1.614097    1.642112    1.760883
     1.824173    2.171539    2.277402    2.337833    2.382347    2.554418    2.971625    3.510205    3.801941    4.057167
     4.269072    4.638369    4.868862    5.184443    5.283406    5.645106    6.006492

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=    -77.8533120308 demc= 6.1099E-05 wnorm= 3.1576E-04 knorm= 3.3293E-04 apxde= 3.4432E-08    *not conv.*     

               starting mcscf iteration...   6

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8868082370     -110.6718925359        0.0000000000        0.0000010000
    2       -77.8198159060     -110.6049002049        0.0000000000        0.0000010000
    3       -77.5528162754     -110.3379005743        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.086607953396653E-006
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    4

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 5.0864E-05 rznorm= 6.0807E-07 rpnorm= 1.2362E-07 noldr=  4 nnewr=  4 nolds=  1 nnews=  1
 

 fdd(*) eigenvalues. symmetry block  1
   -22.585108  -22.369903   -2.072586   -1.629463   -1.268960   -1.120789   -1.039985

 qvv(*) eigenvalues. symmetry block  1
     0.502849    0.567987    0.663414    0.688821    0.907711    1.364654    1.543373    1.614112    1.642095    1.760878
     1.824175    2.171552    2.277415    2.337819    2.382371    2.554399    2.971614    3.510206    3.801972    4.057144
     4.269060    4.638377    4.868888    5.184432    5.283395    5.645107    6.006486

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=    -77.8533120715 demc= 4.0684E-08 wnorm= 5.6693E-05 knorm= 8.0114E-05 apxde= 1.5993E-09    *not conv.*     

               starting mcscf iteration...   7

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     1, naopsy(1) =    36, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore=  65530061

 inoutp: segmentation information:
 in-core transformation space,   avcinc =  65360681
 address segment size,           sizesg =  65209195
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=    65391914 available sort2 space, avcisx=    65392166

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   1 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.8868079543     -110.6718922532        0.0000000000        0.0000010000
    2       -77.8198161919     -110.6049004908        0.0000000000        0.0000010000
    3       -77.5528148515     -110.3378991504        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.525922686627879E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 6.7128E-07 rpnorm= 1.0743E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.585107  -22.369903   -2.072585   -1.629463   -1.268959   -1.120789   -1.039985

 qvv(*) eigenvalues. symmetry block  1
     0.502850    0.567987    0.663415    0.688821    0.907710    1.364654    1.543373    1.614114    1.642093    1.760878
     1.824175    2.171552    2.277415    2.337819    2.382371    2.554399    2.971613    3.510206    3.801973    4.057144
     4.269060    4.638377    4.868889    5.184431    5.283395    5.645107    6.006486

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    7 emc=    -77.8533120731 demc= 1.5993E-09 wnorm= 6.0207E-07 knorm= 4.7795E-08 apxde= 1.3981E-14    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.500 total energy=      -77.886807954, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.500 total energy=      -77.819816192, rel. (eV)=   1.822939
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.49984423
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.50015577     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
         97 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      
 Computing the requested mcscf (transition) density matrices (flag 30)
 Reading mcdenin ...
 Number of density matrices (ndens):                     3
 Number of unique bra states (ndbra):                     2
 qind: F
 (Transition) density matrices:
 d1(*) written to the 1-particle density matrix file.
         97 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st02                                           
 d1(*) written to the 1-particle density matrix file.
         97 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st01                                           
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
         97 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st02-st01                                      

          state spec. NOs: DRT 1, State  2
 *** warning *** small active-orbital occupation. i=  1 nocc= 3.1494E-04
 *** warning *** large active-orbital occupation. i=  2 nocc= 1.9997E+00

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99968506
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.00031494     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000

          state spec. NOs: DRT 1, State  1

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.01645851
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.98354149     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.000291   1.997274   0.648385   0.472184   0.034083   0.028995
     1_ p       0.000417   0.000002   0.095104   0.098371   0.056222   0.750774
     1_ d      -0.000030   0.000000   0.001997   0.002467   0.002687   0.014506
     2_ s       1.998198   0.000584   0.919247   0.407145   0.000092   0.009673
     2_ p      -0.000000   0.000635   0.037265   0.327870   1.096880   0.611917
     2_ d       0.000001  -0.000043   0.004784   0.005405   0.011354   0.005770
     3_ s       0.000050   0.000786   0.081570   0.148001   0.000644   0.010734
     3_ s       0.000015   0.000748   0.049605   0.156439   0.001007   0.414711
     3_ s       0.000529   0.000020   0.073624   0.252922   0.277049   0.139769
     3_ s       0.000530  -0.000007   0.088419   0.129196   0.519983   0.013152
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002646   0.191894   0.000177   0.000000   0.000000   0.000000
     1_ p       0.915884   1.175178   0.009577   0.000000   0.000000   0.000000
     1_ d       0.017423  -0.000134   0.005347   0.000000   0.000000   0.000000
     2_ s       0.006523   0.005462   0.000208   0.000000   0.000000   0.000000
     2_ p       0.300002   0.023369   0.437942   0.000000   0.000000   0.000000
     2_ d       0.005441   0.012806   0.000098   0.000000   0.000000   0.000000
     3_ s       0.562960   0.007604   0.026495   0.000000   0.000000   0.000000
     3_ s       0.169746   0.013046   0.020037   0.000000   0.000000   0.000000
     3_ s       0.019110   0.018605   0.000224   0.000000   0.000000   0.000000
     3_ s       0.000264   0.052015   0.000052   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.375929   3.347132   0.838844   0.825353   0.781852   0.803603
      p         3.101529   2.835880   0.000000   0.000000   0.000000   0.000000
      d         0.044263   0.045616   0.000000   0.000000   0.000000   0.000000
    total       6.521721   6.228627   0.838844   0.825353   0.781852   0.803603
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.000291   1.997274   0.648384   0.472185   0.034083   0.028995
     1_ p       0.000417   0.000002   0.095104   0.098370   0.056222   0.750776
     1_ d      -0.000030   0.000000   0.001997   0.002467   0.002687   0.014506
     2_ s       1.998198   0.000584   0.919249   0.407143   0.000092   0.009672
     2_ p      -0.000000   0.000635   0.037265   0.327872   1.096880   0.611912
     2_ d       0.000001  -0.000043   0.004784   0.005405   0.011354   0.005770
     3_ s       0.000050   0.000786   0.081570   0.148001   0.000644   0.010733
     3_ s       0.000015   0.000748   0.049604   0.156438   0.001007   0.414717
     3_ s       0.000529   0.000020   0.073624   0.252921   0.277049   0.139768
     3_ s       0.000530  -0.000007   0.088419   0.129196   0.519982   0.013152
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002646   0.256101   0.000000   0.000000   0.000000   0.000000
     1_ p       0.915883   1.566320   0.000006   0.000000   0.000000   0.000000
     1_ d       0.017423  -0.000162   0.000003   0.000000   0.000000   0.000000
     2_ s       0.006524   0.007287   0.000000   0.000000   0.000000   0.000000
     2_ p       0.300006   0.031266   0.000276   0.000000   0.000000   0.000000
     2_ d       0.005441   0.017054   0.000000   0.000000   0.000000   0.000000
     3_ s       0.562961   0.009767   0.000017   0.000000   0.000000   0.000000
     3_ s       0.169741   0.017940   0.000013   0.000000   0.000000   0.000000
     3_ s       0.019111   0.024720   0.000000   0.000000   0.000000   0.000000
     3_ s       0.000264   0.069392   0.000000   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.439960   3.348749   0.814529   0.810223   0.787742   0.820929
      p         3.483100   2.406111   0.000000   0.000000   0.000000   0.000000
      d         0.038891   0.049766   0.000000   0.000000   0.000000   0.000000
    total       6.961952   5.804625   0.814529   0.810223   0.787742   0.820929
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.000291   1.997274   0.648384   0.472185   0.034083   0.028995
     1_ p       0.000417   0.000002   0.095104   0.098370   0.056222   0.750776
     1_ d      -0.000030   0.000000   0.001997   0.002467   0.002687   0.014506
     2_ s       1.998198   0.000584   0.919249   0.407143   0.000092   0.009672
     2_ p      -0.000000   0.000635   0.037265   0.327872   1.096880   0.611912
     2_ d       0.000001  -0.000043   0.004784   0.005405   0.011354   0.005770
     3_ s       0.000050   0.000786   0.081570   0.148001   0.000644   0.010733
     3_ s       0.000015   0.000748   0.049604   0.156438   0.001007   0.414717
     3_ s       0.000529   0.000020   0.073624   0.252921   0.277049   0.139768
     3_ s       0.000530  -0.000007   0.088419   0.129196   0.519982   0.013152
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002646   0.057282   0.070756   0.000000   0.000000   0.000000
     1_ p       0.915883   0.423451   0.379731   0.000000   0.000000   0.000000
     1_ d       0.017423   0.004879   0.005706   0.000000   0.000000   0.000000
     2_ s       0.006524   0.001929   0.002125   0.000000   0.000000   0.000000
     2_ p       0.300006   0.449463   0.441617   0.000000   0.000000   0.000000
     2_ d       0.005441   0.005040   0.003713   0.000000   0.000000   0.000000
     3_ s       0.562961   0.040946   0.017468   0.000000   0.000000   0.000000
     3_ s       0.169741   0.007911   0.040302   0.000000   0.000000   0.000000
     3_ s       0.019111   0.009198   0.003742   0.000000   0.000000   0.000000
     3_ s       0.000264   0.016359   0.018382   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.311898   3.345515   0.863159   0.840484   0.775962   0.786278
      p         2.719957   3.265648   0.000000   0.000000   0.000000   0.000000
      d         0.049634   0.041466   0.000000   0.000000   0.000000   0.000000
    total       6.081489   6.652629   0.863159   0.840484   0.775962   0.786278
 

 Total number of electrons:   16.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 02 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s      -0.000000  -0.000001  -0.000000   0.000002  -0.000001  -0.000000
     1_ p      -0.000015   0.000027  -0.000001   0.000007  -0.000001   0.000000
     1_ d      -0.000009   0.000018  -0.000000   0.000023  -0.000005  -0.000000
     2_ s      -0.000000   0.000015   0.000000   0.000000  -0.000005  -0.000000
     2_ p      -0.000000   0.000053  -0.000001   0.000000  -0.000003  -0.000000
     2_ d       0.000000   0.000081  -0.000028   0.000000  -0.000109  -0.000032
     3_ s       0.000000   0.000009  -0.000002  -0.000000  -0.000002  -0.000004
     3_ s      -0.000000   0.000010  -0.000001  -0.000000  -0.000000  -0.000000
     3_ s       0.000000   0.000003  -0.000003   0.000000  -0.000002  -0.000005
     3_ s      -0.000000  -0.000000  -0.000000   0.000000   0.000000  -0.000002
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s      -0.000020  -0.000002   0.000000  -0.000007   0.000000  -0.000000
     1_ p      -0.000100  -0.000031   0.000001  -0.000004   0.000002  -0.000000
     1_ d      -0.000312   0.000061   0.000003  -0.000040  -0.000059  -0.000002
     2_ s      -0.000017   0.000011  -0.000000  -0.000001  -0.000138  -0.000001
     2_ p      -0.000002   0.000145   0.000002   0.000000  -0.000000  -0.000000
     2_ d       0.000000   0.000186   0.000009   0.000000  -0.000255  -0.000003
     3_ s      -0.000003   0.000112   0.000004   0.000000  -0.000000  -0.000001
     3_ s      -0.000002  -0.000008   0.000001   0.000000  -0.000001  -0.000000
     3_ s       0.000000  -0.000007   0.000004  -0.000000   0.000004  -0.000004
     3_ s      -0.000000   0.000005   0.000000  -0.000000   0.000001  -0.000001
 
   ao class      13A        14A        15A        16A        17A        18A  
     1_ s       0.000001   0.000000  -0.000000   0.000000   0.000000  -0.000000
     1_ p       0.000010   0.000000   0.000000   0.000000  -0.000000   0.000000
     1_ d       0.000053  -0.000000   0.000000   0.000000  -0.000006   0.000000
     2_ s       0.000000  -0.000000   0.000000   0.000000  -0.000005   0.000000
     2_ p       0.000001  -0.000000   0.000000   0.000000  -0.000004   0.000000
     2_ d      -0.000000  -0.000000   0.000000   0.000000  -0.000017   0.000001
     3_ s      -0.000000  -0.000000   0.000000  -0.000000  -0.000000   0.000001
     3_ s       0.000000   0.000000   0.000000   0.000000  -0.000000  -0.000000
     3_ s       0.000000   0.000000   0.000000  -0.000000  -0.000000   0.000004
     3_ s      -0.000000  -0.000000   0.000000   0.000000  -0.000000   0.000001
 
   ao class      19A        20A        21A        22A        23A        24A  
     1_ s      -0.000009   0.000000   0.000000   0.000000   0.000000  -0.000000
     1_ p      -0.000023  -0.000006  -0.000000   0.000001  -0.000000  -0.000000
     1_ d      -0.000026  -0.000067   0.000000   0.000002  -0.000004  -0.000001
     2_ s      -0.000000  -0.000033   0.000000   0.000000  -0.000002   0.000000
     2_ p      -0.000001  -0.000003   0.000000  -0.000000  -0.000006  -0.000000
     2_ d      -0.000001  -0.000113   0.000001   0.000000  -0.000001  -0.000003
     3_ s       0.000000   0.000001   0.000000  -0.000000  -0.000000  -0.000000
     3_ s       0.000000  -0.000002   0.000000   0.000000  -0.000000  -0.000000
     3_ s       0.000000  -0.000013   0.000000   0.000000   0.000000  -0.000001
     3_ s      -0.000000   0.000000   0.000000  -0.000000  -0.000000  -0.000001
 
   ao class      25A        26A        27A        28A        29A        30A  
     1_ s       0.000002  -0.000000  -0.000000  -0.000087   0.000000  -0.000001
     1_ p       0.000001   0.000000   0.000000  -0.000076  -0.000049  -0.000049
     1_ d       0.000005  -0.000006  -0.000000  -0.000159  -0.000040  -0.000002
     2_ s       0.000000  -0.000005  -0.000000  -0.000001  -0.000102  -0.000016
     2_ p       0.000000  -0.000008   0.000000  -0.000002  -0.000002  -0.000000
     2_ d       0.000000  -0.000041   0.000008  -0.000002  -0.000055  -0.000019
     3_ s      -0.000000  -0.000001   0.000001  -0.000002  -0.000002  -0.000007
     3_ s      -0.000000  -0.000004  -0.000000   0.000000  -0.000010  -0.000005
     3_ s       0.000000  -0.000003   0.000002  -0.000000  -0.000009  -0.000097
     3_ s       0.000000  -0.000000   0.000000  -0.000001   0.000003  -0.000074
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s       0.000032   0.000016   0.000003  -0.000159   0.000003   0.000000
     1_ p       0.000062   0.000071   0.000001  -0.000269  -0.000011  -0.000027
     1_ d       0.000086   0.000106   0.000000  -0.003310  -0.000273  -0.000001
     2_ s       0.000001   0.000312   0.000170  -0.000008  -0.000331  -0.000002
     2_ p       0.000049   0.000072   0.000007  -0.000042  -0.000003  -0.000002
     2_ d       0.000018   0.000215   0.000694  -0.000058  -0.000095  -0.000037
     3_ s       0.000005  -0.000017  -0.000002   0.000018  -0.000001  -0.000005
     3_ s      -0.000009   0.000090   0.000017   0.000011  -0.000046  -0.000003
     3_ s       0.000016  -0.000001   0.000101  -0.000039  -0.000035  -0.000009
     3_ s       0.000007   0.000002   0.000000  -0.000000   0.000001  -0.000013


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s        -0.000227  -0.000157   0.000101   0.000038  -0.000092  -0.000072
      p        -0.000482   0.000248   0.000000   0.000000   0.000000   0.000000
      d        -0.003967   0.000348   0.000000   0.000000   0.000000   0.000000
    total      -0.004676   0.000440   0.000101   0.000038  -0.000092  -0.000072
 

 Total number of electrons:   -0.00426005

 !timer: mcscf                           cpu_time=     0.232 walltime=     0.232
