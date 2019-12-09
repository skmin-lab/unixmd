

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


 input integral file : /home/swmoon/work/unixmd-fixcol/src/TRAJ.bomd/md/scr_qm/W
 ORK

 Integral file header information:
 Hermit Integral Program : SIFS version  odin              18:49:55.859 09-Dec-19

 Core type energy values:
 energy( 1)=  3.280297367831E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   32.802973678


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

 Input MO coefficient file: /home/swmoon/work/unixmd-fixcol/src/TRAJ.bomd/md/scr_qm/WORK
 

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
    1*      -78.0464579815     -110.8494316598        0.0000000000        0.0000010000
    2       -77.6669561941     -110.4699298725        0.0000000000        0.0000010000
    3       -77.4719971784     -110.2749708567        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.828928178178175E-002
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99839792 pnorm= 0.0000E+00 rznorm= 9.5375E-06 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.535186  -22.530198   -2.082473   -1.577304   -1.285430   -1.184988   -1.012714

 qvv(*) eigenvalues. symmetry block  1
     0.506411    0.550634    0.590649    0.762836    0.959098    1.321763    1.526352    1.527724    1.714520    1.736106
     1.859907    2.160793    2.243863    2.342171    2.369324    2.626873    2.899281    3.419744    3.635635    4.194911
     4.333154    4.525501    4.750230    5.133656    5.345286    5.936941    6.054489

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -77.8567070878 demc= 7.7857E+01 wnorm= 1.4631E-01 knorm= 5.6583E-02 apxde= 1.6873E-03    *not conv.*     

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
    1*      -78.0447371046     -110.8477107829        0.0000000000        0.0000100000
    2       -77.6721178554     -110.4750915337        0.0000000000        0.0000100000
    3       -77.4811873737     -110.2841610520        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.482387416701491E-004
 Total number of micro iterations:   11

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99956195 pnorm= 0.0000E+00 rznorm= 6.3186E-07 rpnorm= 0.0000E+00 noldr= 11 nnewr= 11 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.495025  -22.490279   -2.057258   -1.568791   -1.274237   -1.166793   -1.003992

 qvv(*) eigenvalues. symmetry block  1
     0.510596    0.558347    0.594170    0.765944    0.964497    1.334050    1.537867    1.547507    1.723221    1.754340
     1.876286    2.176103    2.260041    2.349548    2.376678    2.644523    2.906198    3.442090    3.658716    4.216473
     4.354966    4.544637    4.773412    5.153983    5.366992    5.956093    6.079973

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -77.8584274800 demc= 1.7204E-03 wnorm= 1.9859E-03 knorm= 2.9596E-02 apxde= 1.4596E-05    *not conv.*     

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
    1*      -78.0444842814     -110.8474579597        0.0000000000        0.0000010000
    2       -77.6724008359     -110.4753745143        0.0000000000        0.0000010000
    3       -77.4806503121     -110.2836239904        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.116393838034815E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99999950 pnorm= 0.0000E+00 rznorm= 5.5406E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.494291  -22.490000   -2.056901   -1.568574   -1.274026   -1.166103   -1.003541

 qvv(*) eigenvalues. symmetry block  1
     0.510549    0.558439    0.594131    0.765995    0.964445    1.334197    1.538005    1.547683    1.723286    1.754461
     1.876440    2.176290    2.260228    2.349742    2.376612    2.644740    2.906317    3.442401    3.658944    4.216742
     4.355215    4.544801    4.773655    5.154300    5.367222    5.956290    6.080214

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -77.8584425587 demc= 1.5079E-05 wnorm= 8.9312E-04 knorm= 9.9808E-04 apxde= 9.3754E-08    *not conv.*     

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
    1*      -78.0444660493     -110.8474397276        0.0000000000        0.0000010000
    2       -77.6724193450     -110.4753930233        0.0000000000        0.0000010000
    3       -77.4805928589     -110.2835665372        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  6.378517615465427E-005
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99999997 pnorm= 0.0000E+00 rznorm= 8.4556E-07 rpnorm= 0.0000E+00 noldr=  7 nnewr=  7 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.494296  -22.490085   -2.056922   -1.568585   -1.274041   -1.166106   -1.003540

 qvv(*) eigenvalues. symmetry block  1
     0.510540    0.558430    0.594119    0.765987    0.964433    1.334186    1.537990    1.547666    1.723274    1.754450
     1.876421    2.176283    2.260199    2.349780    2.376542    2.644721    2.906309    3.442378    3.658922    4.216723
     4.355195    4.544779    4.773632    5.154282    5.367199    5.956268    6.080186

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -77.8584426971 demc= 1.3843E-07 wnorm= 5.1028E-04 knorm= 2.4688E-04 apxde= 2.4741E-08    *not conv.*     

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
    1*      -78.0444628260     -110.8474365043        0.0000000000        0.0000010000
    2       -77.6724226458     -110.4753963241        0.0000000000        0.0000010000
    3       -77.4805817201     -110.2835553984        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.661214475213587E-005
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99999996 pnorm= 5.9876E-04 rznorm= 5.4533E-07 rpnorm= 3.0576E-20 noldr=  7 nnewr=  7 nolds=  2 nnews=  2
 

 fdd(*) eigenvalues. symmetry block  1
   -22.494283  -22.490111   -2.056927   -1.568587   -1.274045   -1.166108   -1.003541

 qvv(*) eigenvalues. symmetry block  1
     0.510540    0.558429    0.594114    0.765985    0.964432    1.334183    1.537987    1.547663    1.723272    1.754448
     1.876417    2.176286    2.260191    2.349808    2.376509    2.644718    2.906307    3.442373    3.658919    4.216719
     4.355192    4.544774    4.773628    5.154278    5.367195    5.956264    6.080181

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=    -77.8584427359 demc= 3.8801E-08 wnorm= 2.9290E-04 knorm= 2.8462E-04 apxde= 1.8735E-08    *not conv.*     

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
    1*      -78.0444605946     -110.8474342729        0.0000000000        0.0000010000
    2       -77.6724249147     -110.4753985930        0.0000000000        0.0000010000
    3       -77.4805722960     -110.2835459743        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.204501907734334E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 4.8038E-07 rpnorm= 3.0865E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.494261  -22.490140   -2.056929   -1.568588   -1.274048   -1.166108   -1.003541

 qvv(*) eigenvalues. symmetry block  1
     0.510540    0.558428    0.594110    0.765984    0.964431    1.334182    1.537985    1.547661    1.723270    1.754446
     1.876415    2.176292    2.260182    2.349847    2.376466    2.644716    2.906306    3.442371    3.658918    4.216717
     4.355191    4.544772    4.773626    5.154276    5.367192    5.956261    6.080178

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=    -77.8584427547 demc= 1.8736E-08 wnorm= 5.7636E-07 knorm= 2.2295E-08 apxde= 6.9866E-15    *not conv.*     

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
    1*      -78.0444605960     -110.8474342743        0.0000000000        0.0000010000
    2       -77.6724249133     -110.4753985916        0.0000000000        0.0000010000
    3       -77.4805722948     -110.2835459731        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  6.009668769112967E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 3.8948E-07 rpnorm= 4.0995E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -22.494261  -22.490140   -2.056929   -1.568588   -1.274048   -1.166108   -1.003541

 qvv(*) eigenvalues. symmetry block  1
     0.510540    0.558428    0.594110    0.765984    0.964431    1.334182    1.537985    1.547661    1.723270    1.754446
     1.876415    2.176292    2.260182    2.349847    2.376466    2.644716    2.906306    3.442371    3.658918    4.216717
     4.355191    4.544772    4.773626    5.154276    5.367192    5.956261    6.080178

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    7 emc=    -77.8584427547 demc= 4.2633E-14 wnorm= 4.8077E-07 knorm= 1.7080E-08 apxde= 3.9644E-15    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.500 total energy=      -78.044460596, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.500 total energy=      -77.672424913, rel. (eV)=  10.123610
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.46333857
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.53666143     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
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
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.00877845
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.99122155     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000

          state spec. NOs: DRT 1, State  1

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.92666151
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.07333849     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
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
     1_ s       0.249774   1.748872   0.821791   0.464297  -0.000010   0.006783
     1_ p       0.000335   0.000120   0.068020   0.191351   0.626272   0.704234
     1_ d      -0.000018  -0.000026   0.005314   0.004164   0.000600   0.007283
     2_ s       1.748878   0.249792   0.805328   0.480874   0.000271   0.010177
     2_ p      -0.000011   0.000459   0.071958   0.167539   0.600798   0.719911
     2_ d       0.000016  -0.000060   0.005306   0.003929   0.000237   0.007785
     3_ s       0.000093   0.000398   0.060933   0.164901   0.208120   0.088154
     3_ s       0.000083   0.000379   0.053663   0.184119   0.189224   0.175964
     3_ s       0.000425   0.000034   0.053599   0.180072   0.174651   0.185777
     3_ s       0.000425   0.000032   0.054088   0.158755   0.199837   0.093932
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000108  -0.000006   0.001811   0.000000   0.000000   0.000000
     1_ p       0.475083   0.715078   0.265610   0.000000   0.000000   0.000000
     1_ d       0.016593   0.007589   0.002172   0.000000   0.000000   0.000000
     2_ s      -0.000051   0.000097   0.000062   0.000000   0.000000   0.000000
     2_ p       0.499810   0.729708   0.262471   0.000000   0.000000   0.000000
     2_ d       0.016475   0.007450   0.002130   0.000000   0.000000   0.000000
     3_ s       0.288088   0.001103   0.001109   0.000000   0.000000   0.000000
     3_ s       0.195154   0.000171   0.000806   0.000000   0.000000   0.000000
     3_ s       0.206242  -0.000010   0.000336   0.000000   0.000000   0.000000
     3_ s       0.302497   0.002157   0.000154   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.293419   3.295428   0.812900   0.799564   0.801128   0.811876
      p         3.046102   3.052643   0.000000   0.000000   0.000000   0.000000
      d         0.043670   0.043269   0.000000   0.000000   0.000000   0.000000
    total       6.383191   6.391340   0.812900   0.799564   0.801128   0.811876
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.249773   1.748873   0.821791   0.464297  -0.000010   0.006783
     1_ p       0.000335   0.000120   0.068020   0.191351   0.626272   0.704234
     1_ d      -0.000018  -0.000026   0.005314   0.004164   0.000600   0.007283
     2_ s       1.748879   0.249791   0.805328   0.480874   0.000271   0.010177
     2_ p      -0.000011   0.000459   0.071958   0.167539   0.600798   0.719911
     2_ d       0.000016  -0.000060   0.005306   0.003929   0.000237   0.007785
     3_ s       0.000093   0.000398   0.060933   0.164901   0.208120   0.088154
     3_ s       0.000083   0.000379   0.053663   0.184119   0.189224   0.175964
     3_ s       0.000425   0.000034   0.053599   0.180072   0.174651   0.185777
     3_ s       0.000425   0.000032   0.054088   0.158755   0.199837   0.093932
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000108   0.000840   0.002516   0.000000   0.000000   0.000000
     1_ p       0.475083  -0.072532   1.046227   0.000000   0.000000   0.000000
     1_ d       0.016593   0.009510  -0.000192   0.000000   0.000000   0.000000
     2_ s      -0.000051   0.000061   0.000120   0.000000   0.000000   0.000000
     2_ p       0.499810   1.066081  -0.068456   0.000000   0.000000   0.000000
     2_ d       0.016475  -0.000150   0.009129   0.000000   0.000000   0.000000
     3_ s       0.288088   0.002734   0.000110   0.000000   0.000000   0.000000
     3_ s       0.195154   0.000228   0.001381   0.000000   0.000000   0.000000
     3_ s       0.206242   0.000363   0.000257   0.000000   0.000000   0.000000
     3_ s       0.302497   0.001644   0.000130   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.294970   3.295450   0.813531   0.800196   0.801422   0.811339
      p         3.039110   3.058088   0.000000   0.000000   0.000000   0.000000
      d         0.043227   0.042667   0.000000   0.000000   0.000000   0.000000
    total       6.377306   6.396206   0.813531   0.800196   0.801422   0.811339
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.249773   1.748873   0.821791   0.464297  -0.000010   0.006783
     1_ p       0.000335   0.000120   0.068020   0.191351   0.626272   0.704234
     1_ d      -0.000018  -0.000026   0.005314   0.004164   0.000600   0.007283
     2_ s       1.748879   0.249791   0.805328   0.480874   0.000271   0.010177
     2_ p      -0.000011   0.000459   0.071958   0.167539   0.600798   0.719911
     2_ d       0.000016  -0.000060   0.005306   0.003929   0.000237   0.007785
     3_ s       0.000093   0.000398   0.060933   0.164901   0.208120   0.088154
     3_ s       0.000083   0.000379   0.053663   0.184119   0.189224   0.175964
     3_ s       0.000425   0.000034   0.053599   0.180072   0.174651   0.185777
     3_ s       0.000425   0.000032   0.054088   0.158755   0.199837   0.093932
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000108   0.000007   0.000247   0.000000   0.000000   0.000000
     1_ p       0.475083   0.951775   0.035906   0.000000   0.000000   0.000000
     1_ d       0.016593   0.009904   0.000300   0.000000   0.000000   0.000000
     2_ s      -0.000051   0.000129   0.000008   0.000000   0.000000   0.000000
     2_ p       0.499810   0.950475   0.036260   0.000000   0.000000   0.000000
     2_ d       0.016475   0.009895   0.000288   0.000000   0.000000   0.000000
     3_ s       0.288088   0.001428   0.000152   0.000000   0.000000   0.000000
     3_ s       0.195154   0.000236   0.000110   0.000000   0.000000   0.000000
     3_ s       0.206242  -0.000014   0.000046   0.000000   0.000000   0.000000
     3_ s       0.302497   0.002827   0.000022   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.291868   3.295406   0.812268   0.798932   0.800834   0.812414
      p         3.053095   3.047198   0.000000   0.000000   0.000000   0.000000
      d         0.044114   0.043871   0.000000   0.000000   0.000000   0.000000
    total       6.389076   6.386475   0.812268   0.798932   0.800834   0.812414
 

 Total number of electrons:   16.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 02 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s      -0.000000   0.000000  -0.000000   0.000000   0.000000   0.000000
     1_ p      -0.000000   0.000000  -0.000000   0.000031  -0.000003  -0.000000
     1_ d      -0.000000   0.000001  -0.000002   0.000002  -0.000000  -0.000000
     2_ s      -0.000000   0.000001  -0.000000   0.000000  -0.000000  -0.000000
     2_ p      -0.000000   0.000001  -0.000001  -0.000000  -0.000000  -0.000000
     2_ d      -0.000000   0.000003  -0.000000   0.000000  -0.000166  -0.000000
     3_ s       0.000000   0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     3_ s      -0.000000   0.000000  -0.000000  -0.000000   0.000000  -0.000000
     3_ s      -0.000000   0.000000  -0.000000   0.000001   0.000000  -0.000000
     3_ s      -0.000000   0.000000  -0.000000   0.000000  -0.000000   0.000000
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s      -0.000000  -0.000000   0.000000  -0.000000   0.000000   0.000000
     1_ p      -0.000004  -0.000000  -0.000000  -0.000004   0.000000   0.000000
     1_ d      -0.000000  -0.000000  -0.000000  -0.000001   0.000003   0.000000
     2_ s       0.000000  -0.000001  -0.000000  -0.000000   0.000000   0.000000
     2_ p      -0.000000   0.000000  -0.000000  -0.000000  -0.000000   0.000000
     2_ d      -0.000000  -0.000000  -0.000000  -0.000000   0.000001   0.000000
     3_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000   0.000000
     3_ s       0.000000  -0.000000  -0.000001   0.000000   0.000000   0.000000
     3_ s       0.000000   0.000000  -0.000000  -0.000000   0.000000  -0.000000
     3_ s      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
     1_ s      -0.000000   0.000000  -0.000000  -0.000001   0.000000  -0.000000
     1_ p      -0.000000  -0.000000  -0.000000  -0.000000   0.000000   0.000000
     1_ d      -0.000000  -0.000001  -0.000000  -0.000000   0.000002  -0.000000
     2_ s       0.000000  -0.000001   0.000000  -0.000000   0.000000  -0.000000
     2_ p      -0.000000  -0.000000  -0.000000  -0.000000   0.000003  -0.000000
     2_ d      -0.000000  -0.000000  -0.000002  -0.000000   0.000004  -0.000000
     3_ s      -0.000000  -0.000000   0.000000  -0.000000   0.000000   0.000000
     3_ s       0.000000   0.000000  -0.000000  -0.000000   0.000000  -0.000000
     3_ s      -0.000000  -0.000000   0.000000  -0.000000   0.000000  -0.000000
     3_ s      -0.000000  -0.000000  -0.000000   0.000000  -0.000000   0.000000
 
   ao class      19A        20A        21A        22A        23A        24A  
     1_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000000
     1_ p       0.000000   0.000000   0.000000   0.000000  -0.000000  -0.000000
     1_ d       0.000000  -0.000000  -0.000000   0.000000  -0.000001  -0.000000
     2_ s       0.000000  -0.000000  -0.000000   0.000000  -0.000000   0.000000
     2_ p       0.000000  -0.000000   0.000000   0.000000  -0.000000   0.000000
     2_ d       0.000000  -0.000004  -0.000000   0.000000  -0.000000  -0.000000
     3_ s      -0.000000  -0.000000   0.000000  -0.000000  -0.000000  -0.000000
     3_ s      -0.000000   0.000000  -0.000000  -0.000000   0.000000   0.000000
     3_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000  -0.000000
     3_ s       0.000000  -0.000000   0.000000   0.000000   0.000000  -0.000000
 
   ao class      25A        26A        27A        28A        29A        30A  
     1_ s      -0.000002  -0.000000  -0.000000   0.000000  -0.000000  -0.000000
     1_ p      -0.000001  -0.000000  -0.000000   0.000025   0.000000  -0.000259
     1_ d      -0.000001  -0.000000  -0.000000   0.000002   0.000017  -0.000002
     2_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000   0.000009
     2_ p      -0.000000  -0.000000  -0.000000  -0.000000   0.000000  -0.000004
     2_ d      -0.000001  -0.000001  -0.000000   0.000000   0.000000  -0.000688
     3_ s       0.000000   0.000000   0.000000  -0.000000   0.000000   0.000003
     3_ s       0.000000  -0.000000  -0.000000  -0.000000   0.000000   0.000002
     3_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000  -0.000104
     3_ s       0.000000  -0.000000   0.000000   0.000000   0.000000   0.000026
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s      -0.000000  -0.000001   0.000005  -0.000000  -0.000000   0.000000
     1_ p      -0.000009  -0.000041  -0.000482  -0.000015  -0.000004   0.000002
     1_ d      -0.000000   0.000000  -0.000061  -0.000001  -0.000002   0.000000
     2_ s      -0.000000  -0.000479  -0.000719   0.000000  -0.000060   0.000001
     2_ p      -0.000000  -0.000000  -0.000004  -0.000000  -0.000000  -0.000000
     2_ d      -0.000000  -0.000014  -0.000329  -0.000000  -0.000000  -0.000000
     3_ s      -0.000000   0.000001   0.000007  -0.000000  -0.000000  -0.000000
     3_ s       0.000000  -0.000012  -0.000000   0.000001  -0.000008   0.000012
     3_ s      -0.000000   0.000005   0.000008  -0.000002   0.000001   0.000000
     3_ s      -0.000001   0.000001  -0.000000   0.000000  -0.000001  -0.000000


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         0.000002  -0.001250   0.000010  -0.000008  -0.000092   0.000024
      p        -0.000764  -0.000006   0.000000   0.000000   0.000000   0.000000
      d        -0.000044  -0.001196   0.000000   0.000000   0.000000   0.000000
    total      -0.000806  -0.002453   0.000010  -0.000008  -0.000092   0.000024
 

 Total number of electrons:   -0.00332524

 !timer: mcscf                           cpu_time=     0.233 walltime=     0.233
