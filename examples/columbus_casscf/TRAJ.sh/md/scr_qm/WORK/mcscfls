

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


 input integral file : /home/swmoon/work/unixmd-fixcol/src/TRAJ.sh/md/scr_qm/WOR
 K/a

 Integral file header information:
 Hermit Integral Program : SIFS version  odin              16:39:39.058 09-Dec-19

 Core type energy values:
 energy( 1)=  3.250102702831E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   32.501027028


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

 Input MO coefficient file: /home/swmoon/work/unixmd-fixcol/src/TRAJ.sh/md/scr_qm/WORK/m
 

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

 trial vector  1 is unit matrix column     1
 ciiter=   3 noldhv=  3 noldv=  3

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -77.9198331106     -110.4208601389        0.0000000000        0.0000010000
    2       -77.8006845065     -110.3017115348        0.0000000000        0.0000010000
    3       -77.6387696573     -110.1397966856        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.205150261319313E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.98206757 pnorm= 0.0000E+00 rznorm= 2.9054E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.576825  -22.486751   -2.102944   -1.571240   -1.248569   -1.118577   -1.046156

 qvv(*) eigenvalues. symmetry block  1
     0.483288    0.556343    0.590590    0.670549    0.872394    1.375821    1.465153    1.546169    1.664401    1.715988
     1.848042    2.206282    2.240877    2.369470    2.413618    2.442099    2.764159    3.490577    3.644602    3.975544
     4.370141    4.496624    4.783972    5.151786    5.273594    5.495178    5.928996

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -77.8602588086 demc= 7.7860E+01 wnorm= 1.7641E-01 knorm= 1.8853E-01 apxde= 1.1110E-02    *not conv.*     

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
    1*      -77.9189084533     -110.4199354816        0.0000000000        0.0000100000
    2       -77.8246950513     -110.3257220796        0.0000000000        0.0000100000
    3       -77.6748859004     -110.1759129287        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.176780658640676E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99995304 pnorm= 0.0000E+00 rznorm= 6.5749E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.493330  -22.474303   -2.086084   -1.555286   -1.237698   -1.109898   -1.030685

 qvv(*) eigenvalues. symmetry block  1
     0.486908    0.561321    0.597649    0.675316    0.884798    1.388281    1.471698    1.560336    1.667965    1.734950
     1.865953    2.217060    2.249010    2.385123    2.422691    2.455728    2.773048    3.510386    3.663086    4.002907
     4.390620    4.518989    4.804401    5.169776    5.289228    5.513627    5.953347

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -77.8718017523 demc= 1.1543E-02 wnorm= 9.4142E-03 knorm= 9.6909E-03 apxde= 2.9260E-05    *not conv.*     

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
    1*      -77.9190006928     -110.4200277211        0.0000000000        0.0000010000
    2       -77.8246619998     -110.3256890281        0.0000000000        0.0000010000
    3       -77.6747816221     -110.1758086504        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.765473299140537E-005
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999998 pnorm= 0.0000E+00 rznorm= 8.8407E-07 rpnorm= 0.0000E+00 noldr=  5 nnewr=  5 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.493759  -22.470816   -2.085066   -1.554707   -1.236279   -1.108586   -1.029964

 qvv(*) eigenvalues. symmetry block  1
     0.487449    0.561897    0.598176    0.675902    0.885211    1.388607    1.472696    1.560719    1.665232    1.734986
     1.866611    2.218196    2.249944    2.385710    2.423625    2.456427    2.773997    3.510936    3.663961    4.002657
     4.391169    4.519552    4.804855    5.170611    5.290223    5.514824    5.953983

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -77.8718313463 demc= 2.9594E-05 wnorm= 3.0124E-04 knorm= 2.0685E-04 apxde= 2.0355E-08    *not conv.*     

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
    1*      -77.9190048646     -110.4200318929        0.0000000000        0.0000010000
    2       -77.8246578720     -110.3256849003        0.0000000000        0.0000010000
    3       -77.6748025105     -110.1758295388        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.857897572504278E-006
 Total number of micro iterations:    4

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-1.00000000 pnorm= 0.0000E+00 rznorm= 3.8860E-07 rpnorm= 0.0000E+00 noldr=  4 nnewr=  4 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.493749  -22.470811   -2.085063   -1.554706   -1.236276   -1.108576   -1.029963

 qvv(*) eigenvalues. symmetry block  1
     0.487453    0.561899    0.598180    0.675905    0.885212    1.388608    1.472699    1.560719    1.665218    1.734987
     1.866612    2.218200    2.249952    2.385713    2.423631    2.456433    2.774001    3.510937    3.663963    4.002664
     4.391169    4.519554    4.804858    5.170613    5.290228    5.514827    5.953985

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -77.8718313683 demc= 2.2042E-08 wnorm= 3.0863E-05 knorm= 1.8558E-05 apxde= 1.7123E-10    *not conv.*     

               starting mcscf iteration...   5

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
    1*      -77.9190050342     -110.4200320626        0.0000000000        0.0000010000
    2       -77.8246577028     -110.3256847311        0.0000000000        0.0000010000
    3       -77.6748040268     -110.1758310551        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.861915502948886E-007
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    2

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 6.2220E-07 rpnorm= 3.0371E-07 noldr=  2 nnewr=  2 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.493748  -22.470811   -2.085063   -1.554706   -1.236276   -1.108576   -1.029963

 qvv(*) eigenvalues. symmetry block  1
     0.487454    0.561898    0.598180    0.675905    0.885212    1.388608    1.472699    1.560719    1.665218    1.734987
     1.866612    2.218200    2.249952    2.385713    2.423631    2.456433    2.774002    3.510937    3.663963    4.002664
     4.391169    4.519554    4.804858    5.170613    5.290229    5.514827    5.953985

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    5 emc=    -77.8718313685 demc= 1.9106E-10 wnorm= 3.8895E-06 knorm= 2.1159E-06 apxde= 2.4482E-12    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.500 total energy=      -77.919005034, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.500 total energy=      -77.824657703, rel. (eV)=   2.567323
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.47858014
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.52141986     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
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

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.06502293
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.93497707     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000

          state spec. NOs: DRT 1, State  1

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.95651754
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.04348246     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
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
     1_ s       1.996096   0.002953   0.518955   0.689733   0.001279   0.061041
     1_ p      -0.000015   0.000463   0.084307   0.105326   0.184993   0.093372
     1_ d       0.000001  -0.000033   0.003967   0.003903   0.004220   0.007377
     2_ s       0.002801   1.994846   1.039096   0.274123   0.004436   0.000189
     2_ p       0.000346   0.000009   0.042694   0.308256   1.014962   1.003787
     2_ d      -0.000038  -0.000000   0.001380   0.005495   0.006939   0.015951
     3_ s       0.000134   0.000634   0.124616   0.027057   0.069647   0.401090
     3_ s       0.000653   0.000009   0.037865   0.306843   0.092776   0.040077
     3_ s       0.000022   0.000619   0.070724   0.168544   0.062699   0.373099
     3_ s       0.000000   0.000500   0.076395   0.110721   0.558049   0.004017
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002072   0.245139  -0.000446   0.000000   0.000000   0.000000
     1_ p       0.646279   1.056164   0.483938   0.000000   0.000000   0.000000
     1_ d       0.014048   0.000655   0.000249   0.000000   0.000000   0.000000
     2_ s       0.007289   0.008176   0.000760   0.000000   0.000000   0.000000
     2_ p       0.744078   0.043504   0.010897   0.000000   0.000000   0.000000
     2_ d       0.011367   0.007276   0.003900   0.000000   0.000000   0.000000
     3_ s       0.158769   0.007572   0.017324   0.000000   0.000000   0.000000
     3_ s       0.288775   0.065113   0.000055   0.000000   0.000000   0.000000
     3_ s       0.113746   0.009214   0.001817   0.000000   0.000000   0.000000
     3_ s       0.013578   0.035767   0.002926   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.516822   3.331715   0.806841   0.832166   0.800484   0.801954
      p         2.654826   3.168534   0.000000   0.000000   0.000000   0.000000
      d         0.034387   0.052270   0.000000   0.000000   0.000000   0.000000
    total       6.206035   6.552520   0.806841   0.832166   0.800484   0.801954
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       1.996096   0.002953   0.518955   0.689733   0.001279   0.061041
     1_ p      -0.000015   0.000463   0.084307   0.105326   0.184993   0.093371
     1_ d       0.000001  -0.000033   0.003967   0.003903   0.004220   0.007377
     2_ s       0.002800   1.994847   1.039096   0.274122   0.004436   0.000189
     2_ p       0.000346   0.000009   0.042694   0.308257   1.014962   1.003787
     2_ d      -0.000038  -0.000000   0.001380   0.005495   0.006939   0.015951
     3_ s       0.000134   0.000634   0.124616   0.027056   0.069648   0.401089
     3_ s       0.000653   0.000009   0.037865   0.306843   0.092776   0.040077
     3_ s       0.000022   0.000619   0.070724   0.168544   0.062699   0.373099
     3_ s       0.000000   0.000500   0.076395   0.110721   0.558049   0.004017
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002072   0.084980   0.079610   0.000000   0.000000   0.000000
     1_ p       0.646279   0.901267   0.744413   0.000000   0.000000   0.000000
     1_ d       0.014048   0.000102   0.000771   0.000000   0.000000   0.000000
     2_ s       0.007289   0.003725   0.003263   0.000000   0.000000   0.000000
     2_ p       0.744078   0.010113   0.038171   0.000000   0.000000   0.000000
     2_ d       0.011367   0.007235   0.005243   0.000000   0.000000   0.000000
     3_ s       0.158769   0.006940   0.029759   0.000000   0.000000   0.000000
     3_ s       0.288775   0.021355   0.022525   0.000000   0.000000   0.000000
     3_ s       0.113746   0.000438   0.008699   0.000000   0.000000   0.000000
     3_ s       0.013578   0.028866   0.002523   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.436719   3.329767   0.818645   0.810878   0.798590   0.794650
      p         2.760404   3.162418   0.000000   0.000000   0.000000   0.000000
      d         0.034357   0.053572   0.000000   0.000000   0.000000   0.000000
    total       6.231480   6.545757   0.818645   0.810878   0.798590   0.794650
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       1.996096   0.002953   0.518955   0.689733   0.001279   0.061041
     1_ p      -0.000015   0.000463   0.084307   0.105326   0.184993   0.093371
     1_ d       0.000001  -0.000033   0.003967   0.003903   0.004220   0.007377
     2_ s       0.002800   1.994847   1.039096   0.274122   0.004436   0.000189
     2_ p       0.000346   0.000009   0.042694   0.308257   1.014962   1.003787
     2_ d      -0.000038  -0.000000   0.001380   0.005495   0.006939   0.015951
     3_ s       0.000134   0.000634   0.124616   0.027056   0.069648   0.401089
     3_ s       0.000653   0.000009   0.037865   0.306843   0.092776   0.040077
     3_ s       0.000022   0.000619   0.070724   0.168544   0.062699   0.373099
     3_ s       0.000000   0.000500   0.076395   0.110721   0.558049   0.004017
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.002072   0.324844  -0.000048   0.000000   0.000000   0.000000
     1_ p       0.646279   1.394090   0.040434   0.000000   0.000000   0.000000
     1_ d       0.014048   0.000915   0.000020   0.000000   0.000000   0.000000
     2_ s       0.007289   0.010821   0.000063   0.000000   0.000000   0.000000
     2_ p       0.744078   0.059656   0.000862   0.000000   0.000000   0.000000
     2_ d       0.011367   0.009548   0.000327   0.000000   0.000000   0.000000
     3_ s       0.158769   0.011684   0.001408   0.000000   0.000000   0.000000
     3_ s       0.288775   0.086457  -0.000002   0.000000   0.000000   0.000000
     3_ s       0.113746   0.012786   0.000138   0.000000   0.000000   0.000000
     3_ s       0.013578   0.045718   0.000280   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.596925   3.333663   0.795038   0.853454   0.802377   0.809258
      p         2.549247   3.174651   0.000000   0.000000   0.000000   0.000000
      d         0.034418   0.050969   0.000000   0.000000   0.000000   0.000000
    total       6.180590   6.559283   0.795038   0.853454   0.802377   0.809258
 

 Total number of electrons:   16.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 02 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s      -0.000000   0.000000  -0.000000   0.000002  -0.000000   0.000000
     1_ p      -0.000000  -0.000000   0.000000   0.000006   0.000000  -0.000000
     1_ d      -0.000000  -0.000003   0.000001   0.000012  -0.000000  -0.000000
     2_ s      -0.000000  -0.000000   0.000000   0.000001  -0.000000  -0.000002
     2_ p      -0.000000  -0.000001   0.000000  -0.000000  -0.000001  -0.000001
     2_ d      -0.000000  -0.000001   0.000002   0.000000  -0.000001  -0.000032
     3_ s      -0.000000  -0.000000   0.000002   0.000001   0.000000   0.000000
     3_ s      -0.000000   0.000000   0.000001  -0.000000  -0.000000  -0.000002
     3_ s       0.000000  -0.000001   0.000000   0.000000   0.000000  -0.000000
     3_ s      -0.000000   0.000000   0.000000   0.000000   0.000000  -0.000002
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000002  -0.000000   0.000000   0.000001   0.000000  -0.000000
     1_ p       0.000002  -0.000000  -0.000000   0.000004   0.000000   0.000000
     1_ d       0.000020  -0.000006   0.000000   0.000004  -0.000002   0.000000
     2_ s       0.000001  -0.000000  -0.000000   0.000004  -0.000001   0.000000
     2_ p       0.000001  -0.000001   0.000000   0.000000  -0.000001   0.000000
     2_ d       0.000003  -0.000001   0.000001   0.000004  -0.000001   0.000001
     3_ s       0.000000  -0.000000   0.000007   0.000000  -0.000000   0.000000
     3_ s      -0.000000  -0.000000   0.000001  -0.000000  -0.000001   0.000000
     3_ s       0.000000  -0.000001  -0.000000   0.000000  -0.000000   0.000000
     3_ s       0.000000  -0.000000   0.000002   0.000000  -0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
     1_ s       0.000003  -0.000000   0.000000  -0.000000   0.000000  -0.000000
     1_ p       0.000000  -0.000000  -0.000000  -0.000000   0.000000  -0.000000
     1_ d       0.000001  -0.000001  -0.000000  -0.000000   0.000000  -0.000000
     2_ s       0.000001  -0.000000  -0.000000  -0.000000   0.000000  -0.000000
     2_ p       0.000000  -0.000000  -0.000000  -0.000000   0.000000   0.000000
     2_ d       0.000001  -0.000000  -0.000002  -0.000000   0.000000  -0.000000
     3_ s       0.000000  -0.000000   0.000000  -0.000000  -0.000000  -0.000000
     3_ s       0.000000  -0.000000  -0.000000  -0.000000   0.000000  -0.000000
     3_ s       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
     3_ s      -0.000000  -0.000000  -0.000001  -0.000000   0.000000  -0.000000
 
   ao class      19A        20A        21A        22A        23A        24A  
     1_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000   0.000000
     1_ p      -0.000000   0.000000   0.000000   0.000000  -0.000000   0.000000
     1_ d      -0.000000   0.000004   0.000001   0.000000  -0.000000   0.000000
     2_ s      -0.000000   0.000000  -0.000000   0.000000  -0.000000  -0.000000
     2_ p      -0.000000   0.000001   0.000000  -0.000000  -0.000000   0.000000
     2_ d      -0.000000   0.000001   0.000001   0.000000  -0.000000   0.000000
     3_ s       0.000000  -0.000000   0.000000  -0.000000  -0.000000   0.000000
     3_ s      -0.000000   0.000000   0.000000  -0.000000  -0.000000   0.000000
     3_ s      -0.000000   0.000001   0.000000   0.000000  -0.000000   0.000000
     3_ s      -0.000000   0.000000   0.000000   0.000000  -0.000000   0.000000
 
   ao class      25A        26A        27A        28A        29A        30A  
     1_ s       0.000000   0.000000  -0.000000  -0.000010   0.000026   0.000003
     1_ p       0.000001  -0.000000   0.000000  -0.000007  -0.000184  -0.000022
     1_ d       0.000001  -0.000001   0.000000  -0.000020  -0.001796  -0.000020
     2_ s       0.000000  -0.000001  -0.000000  -0.000003  -0.001186  -0.000011
     2_ p       0.000000  -0.000002   0.000000   0.000000  -0.000100  -0.000001
     2_ d       0.000001  -0.000002   0.000003  -0.000001  -0.001106  -0.000054
     3_ s       0.000000  -0.000000   0.000003  -0.000000   0.000015  -0.000025
     3_ s       0.000000  -0.000001   0.000003   0.000000  -0.000057  -0.000010
     3_ s       0.000000  -0.000000   0.000002   0.000000  -0.000035  -0.000008
     3_ s       0.000000  -0.000000   0.000001   0.000000  -0.000033  -0.000004
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s       0.000002  -0.000001  -0.000002  -0.000000   0.000002  -0.000001
     1_ p       0.000016  -0.000113  -0.000005  -0.000008  -0.000003  -0.000020
     1_ d       0.000187  -0.000990  -0.000001  -0.000074  -0.002527  -0.000000
     2_ s       0.000007  -0.000140  -0.000013  -0.000001  -0.000428   0.000005
     2_ p       0.000003  -0.000097   0.000000   0.000000  -0.000062  -0.000000
     2_ d       0.000013  -0.000039  -0.000048  -0.000011  -0.000480  -0.000026
     3_ s       0.000000  -0.000129  -0.000000  -0.000001   0.000001  -0.000057
     3_ s      -0.000000  -0.000047  -0.000033  -0.000001  -0.000016  -0.000006
     3_ s       0.000001   0.000002  -0.000033  -0.000006  -0.000119  -0.000045
     3_ s       0.000004  -0.000001   0.000003  -0.000000   0.000001  -0.000013


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         0.000027  -0.001769  -0.000184  -0.000169  -0.000243  -0.000042
      p        -0.000331  -0.000260   0.000000   0.000000   0.000000   0.000000
      d        -0.005212  -0.001774   0.000000   0.000000   0.000000   0.000000
    total      -0.005516  -0.003803  -0.000184  -0.000169  -0.000243  -0.000042
 

 Total number of electrons:   -0.00995572

 !timer: mcscf                           cpu_time=     0.176 walltime=     0.176
