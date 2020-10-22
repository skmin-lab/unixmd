

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
       131072000 of real*8 words ( 1000.00 MB) of work space has been allocated.

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
   FCIORB=  1,6,20,1,7,20,1,8,20,1,9,20
   NAVST(1) = 3,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
   WAVST(1,3)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /work/swmoon/Col.25235/TRAJ_shxf/md/scr_qm/WORK/aoints   
    

 Integral file header information:
 Hermit Integral Program : SIFS version  odin1             11:42:40.666 15-Oct-20

 Core type energy values:
 energy( 1)=  3.880154784307E+01, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =   38.801547843


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
  1   ground state          3             0.333 0.333 0.333

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        4

 orbital coefficients are optimized for the ground state (nstate=0).

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       6(  6)    20
       1       7(  7)    20
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
 Number of active orbitals:    4
 Number of active electrons:   6
 Total number of CSFs:        10
 

 faar:   0 active-active rotations allowed out of:   6 possible.


 Number of active-double rotations:        20
 Number of active-active rotations:         0
 Number of double-virtual rotations:      135
 Number of active-virtual rotations:      108
 lenbfsdef=                131071  lenbfs=                   729
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #          55
 class  2 (pq|ri):         #         200
 class  3 (pq|ia):         #        1350
 class  4 (pi|qa):         #        2160
 class  5 (pq|ra):         #        1080
 class  6 (pq|ij)/(pi|qj): #         400
 class  7 (pq|ab):         #        3780
 class  8 (pa|qb):         #        7290
 class  9 p(bp,ai)         #       14580
 class 10p(ai,jp):        #        2700
 class 11p(ai,bj):        #       10935

 Size of orbital-Hessian matrix B:                    37915
 Size of the orbital-state Hessian matrix C:           7890
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:          45805


 Source of the initial MO coeficients:

 Input MO coefficient file: /work/swmoon/Col.25235/TRAJ_shxf/md/scr_qm/WORK/mocoef      
 

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     1
 ciiter=   5 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3907987650     -133.1923466081        0.0000000000        0.0000010000
    2       -94.0547002060     -132.8562480491        0.0000000000        0.0000010000
    3       -93.9934954779     -132.7950433210        0.0000000000        0.0000010000
    4       -93.9371125647     -132.7386604078        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.859309694450713E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.98567853 pnorm= 0.0000E+00 rznorm= 3.3248E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.842830  -23.210044   -3.129606   -2.369772   -2.085437

 qvv(*) eigenvalues. symmetry block  1
    -0.045637    0.041942    0.087340    0.257371    0.395431    0.977853    0.998069    1.046577    1.176586    1.380066
     1.388732    1.596787    1.804799    1.824596    1.844385    2.456322    2.510513    2.916126    3.035722    3.439214
     3.930718    4.104932    4.415812    4.763849    4.871664    5.477139    5.490990
 *** warning *** large active-orbital occupation. i=  4 nocc= 1.9996E+00

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=    -94.1463314830 demc= 9.4146E+01 wnorm= 1.4874E-01 knorm= 1.6864E-01 apxde= 6.4548E-03    *not conv.*     

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3918751306     -133.1934229737        0.0000000000        0.0000100000
    2       -94.0652112783     -132.8667591214        0.0000000000        0.0000100000
    3       -94.0035970573     -132.8051449004        0.0000000000        0.0000100000
    4       -93.9468410121     -132.7483888551        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.790230851505196E-003
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.94887087 pnorm= 0.0000E+00 rznorm= 5.9796E-06 rpnorm= 0.0000E+00 noldr=  7 nnewr=  7 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.801249  -23.163774   -3.104450   -2.352615   -2.065015

 qvv(*) eigenvalues. symmetry block  1
    -0.036837    0.051292    0.090984    0.263712    0.402512    0.993741    1.000051    1.057828    1.184380    1.397492
     1.401661    1.612089    1.817322    1.844521    1.858283    2.473069    2.524490    2.936370    3.053689    3.458460
     3.949387    4.123727    4.437185    4.782854    4.890521    5.500091    5.510947
 *** warning *** large active-orbital occupation. i=  4 nocc= 1.9994E+00

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=    -94.1535611554 demc= 7.2297E-03 wnorm= 2.2322E-02 knorm= 3.1566E-01 apxde= 4.8348E-04    *not conv.*     

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3924296716     -133.1939775147        0.0000000000        0.0000010000
    2       -94.0660665974     -132.8676144405        0.0000000000        0.0000010000
    3       -94.0047808628     -132.8063287059        0.0000000000        0.0000010000
    4       -93.8930565493     -132.6946043924        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.050496746495724E-003
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.90443651 pnorm= 0.0000E+00 rznorm= 2.5051E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.797430  -23.151628   -3.088566   -2.315555   -2.057615

 qvv(*) eigenvalues. symmetry block  1
    -0.035721    0.052058    0.090551    0.263679    0.402462    0.994709    0.999824    1.058001    1.183816    1.399550
     1.402933    1.612378    1.817753    1.846137    1.859405    2.474363    2.525378    2.937592    3.054838    3.459396
     3.950315    4.124501    4.438334    4.783953    4.891788    5.501742    5.512344

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=    -94.1544257106 demc= 8.6456E-04 wnorm= 1.6404E-02 knorm= 4.2661E-01 apxde= 4.4152E-04    *not conv.*     

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3932468562     -133.1947946993        0.0000000000        0.0000010000
    2       -94.0664441920     -132.8679920351        0.0000000000        0.0000010000
    3       -94.0058408755     -132.8073887186        0.0000000000        0.0000010000
    4       -93.8054834697     -132.6070313128        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.052084124601761E-003
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.94586728 pnorm= 0.0000E+00 rznorm= 4.6865E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.793350  -23.152085   -3.036262   -2.188396   -2.044812

 qvv(*) eigenvalues. symmetry block  1
    -0.035232    0.052319    0.090134    0.263514    0.402226    0.995001    0.999917    1.057836    1.183287    1.400494
     1.404092    1.611972    1.817571    1.847215    1.859985    2.475208    2.525701    2.938000    3.055297    3.459683
     3.950495    4.124471    4.438619    4.784423    4.892270    5.502441    5.512881

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=    -94.1551773079 demc= 7.5160E-04 wnorm= 8.4167E-03 knorm= 3.2455E-01 apxde= 2.4621E-04    *not conv.*     

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3934118679     -133.1949597110        0.0000000000        0.0000010000
    2       -94.0668113508     -132.8683591939        0.0000000000        0.0000010000
    3       -94.0062769526     -132.8078247957        0.0000000000        0.0000010000
    4       -93.7929878764     -132.5945357195        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.610680590891884E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99510977 pnorm= 0.0000E+00 rznorm= 5.7252E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.791181  -23.154970   -2.999810   -2.112616   -2.000290

 qvv(*) eigenvalues. symmetry block  1
    -0.035058    0.052412    0.090018    0.263491    0.402175    0.995073    0.999975    1.057789    1.183177    1.400846
     1.404617    1.611863    1.817541    1.847739    1.860305    2.475606    2.525915    2.938220    3.055582    3.459875
     3.950621    4.124505    4.438787    4.784674    4.892536    5.502797    5.513127

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=    -94.1555000571 demc= 3.2275E-04 wnorm= 3.6885E-03 knorm= 9.8775E-02 apxde= 2.8501E-05    *not conv.*     

               starting mcscf iteration...   6

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3933493864     -133.1948972295        0.0000000000        0.0000010000
    2       -94.0669830593     -132.8685309024        0.0000000000        0.0000010000
    3       -94.0062819484     -132.8078297915        0.0000000000        0.0000010000
    4       -93.7936659502     -132.5952137933        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.312469838545728E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99860522 pnorm= 0.0000E+00 rznorm= 7.0399E-07 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790444  -23.156211   -2.992636   -2.111027   -1.978878

 qvv(*) eigenvalues. symmetry block  1
    -0.035030    0.052417    0.090037    0.263496    0.402191    0.995102    1.000026    1.057774    1.183196    1.400912
     1.404748    1.611873    1.817563    1.847879    1.860423    2.475723    2.525995    2.938289    3.055718    3.459943
     3.950671    4.124545    4.438865    4.784758    4.892618    5.502922    5.513193

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    6 emc=    -94.1555381314 demc= 3.8074E-05 wnorm= 1.8500E-03 knorm= 5.2798E-02 apxde= 8.5827E-06    *not conv.*     

               starting mcscf iteration...   7

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3932847255     -133.1948325686        0.0000000000        0.0000010000
    2       -94.0671759359     -132.8687237789        0.0000000000        0.0000010000
    3       -94.0062122426     -132.8077600857        0.0000000000        0.0000010000
    4       -93.7955436273     -132.5970914703        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.730700924756957E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99672578 pnorm= 0.0000E+00 rznorm= 5.0261E-07 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790299  -23.157045   -2.996278   -2.124527   -1.974008

 qvv(*) eigenvalues. symmetry block  1
    -0.035067    0.052384    0.090098    0.263508    0.402228    0.995119    1.000095    1.057780    1.183268    1.400855
     1.404715    1.611910    1.817610    1.847856    1.860435    2.475714    2.525996    2.938287    3.055759    3.459957
     3.950691    4.124569    4.438897    4.784761    4.892632    5.502933    5.513177

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    7 emc=    -94.1555576347 demc= 1.9503E-05 wnorm= 1.3846E-03 knorm= 8.0856E-02 apxde= 1.6420E-05    *not conv.*     

               starting mcscf iteration...   8

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3931865355     -133.1947343785        0.0000000000        0.0000010000
    2       -94.0675140380     -132.8690618810        0.0000000000        0.0000010000
    3       -94.0060984388     -132.8076462819        0.0000000000        0.0000010000
    4       -93.7991887687     -132.6007366118        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.276387380112861E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99138888 pnorm= 0.0000E+00 rznorm= 9.6478E-07 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790407  -23.158229   -3.006920   -2.148568   -1.970497

 qvv(*) eigenvalues. symmetry block  1
    -0.035170    0.052317    0.090223    0.263543    0.402306    0.995135    1.000207    1.057821    1.183432    1.400688
     1.404563    1.612007    1.817711    1.847729    1.860393    2.475629    2.525968    2.938269    3.055773    3.459963
     3.950725    4.124607    4.438932    4.784728    4.892626    5.502887    5.513120

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    8 emc=    -94.1555996708 demc= 4.2036E-05 wnorm= 1.8211E-03 knorm= 1.3095E-01 apxde= 4.2289E-05    *not conv.*     

               starting mcscf iteration...   9

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3930097084     -133.1945575514        0.0000000000        0.0000010000
    2       -94.0681804797     -132.8697283228        0.0000000000        0.0000010000
    3       -94.0059542908     -132.8075021339        0.0000000000        0.0000010000
    4       -93.8073876127     -132.6089354558        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  3.852401610710454E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.97458111 pnorm= 0.0000E+00 rznorm= 6.9140E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.790920  -23.159728   -3.026329   -2.188961   -1.965513

 qvv(*) eigenvalues. symmetry block  1
    -0.035408    0.052182    0.090490    0.263629    0.402479    0.995154    1.000447    1.057941    1.183799    1.400305
     1.404169    1.612246    1.817933    1.847408    1.860263    2.475397    2.525904    2.938242    3.055764    3.459973
     3.950808    4.124692    4.438998    4.784637    4.892602    5.502755    5.512983

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    9 emc=    -94.1557148263 demc= 1.1516E-04 wnorm= 3.0819E-03 knorm= 2.2403E-01 apxde= 1.3085E-04    *not conv.*     

               starting mcscf iteration...  10

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3927164830     -133.1942643260        0.0000000000        0.0000010000
    2       -94.0696091792     -132.8711570222        0.0000000000        0.0000010000
    3       -94.0059480074     -132.8074958505        0.0000000000        0.0000010000
    4       -93.8300496461     -132.6315974892        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.462469324490928E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.93756637 pnorm= 0.0000E+00 rznorm= 3.7082E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.792767  -23.160435   -3.058981   -2.253882   -1.956693

 qvv(*) eigenvalues. symmetry block  1
    -0.036051    0.051861    0.091181    0.263872    0.402945    0.995172    1.001139    1.058316    1.184779    1.399298
     1.403046    1.612949    1.818478    1.846499    1.859900    2.474717    2.525750    2.938236    3.055692    3.460034
     3.951085    4.124947    4.439177    4.784390    4.892568    5.502391    5.512622

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   10 emc=    -94.1560912232 demc= 3.7640E-04 wnorm= 5.9700E-03 knorm= 3.4781E-01 apxde= 4.0478E-04    *not conv.*     

               starting mcscf iteration...  11

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3925387254     -133.1940865685        0.0000000000        0.0000010000
    2       -94.0719381926     -132.8734860357        0.0000000000        0.0000010000
    3       -94.0066435756     -132.8081914186        0.0000000000        0.0000010000
    4       -93.8789384147     -132.6804862578        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.265475635229796E-003
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.95525487 pnorm= 0.0000E+00 rznorm= 4.9306E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.798134  -23.156570   -3.094315   -2.324758   -1.941603

 qvv(*) eigenvalues. symmetry block  1
    -0.037666    0.051056    0.092914    0.264512    0.404119    0.995048    1.003065    1.059237    1.187285    1.396831
     1.400167    1.614776    1.819409    1.844134    1.859464    2.472982    2.525366    2.938301    3.055444    3.460303
     3.951849    4.125630    4.439652    4.783825    4.892558    5.501488    5.511769

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   11 emc=    -94.1570401645 demc= 9.4894E-04 wnorm= 1.0124E-02 knorm= 2.9578E-01 apxde= 4.0625E-04    *not conv.*     

               starting mcscf iteration...  12

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3926509738     -133.1941988168        0.0000000000        0.0000010000
    2       -94.0729439754     -132.8744918185        0.0000000000        0.0000010000
    3       -94.0072730817     -132.8088209248        0.0000000000        0.0000010000
    4       -93.9032288479     -132.7047766909        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  9.391205621189203E-004
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99181049 pnorm= 0.0000E+00 rznorm= 4.5118E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.799763  -23.154079   -3.104004   -2.349526   -1.931141

 qvv(*) eigenvalues. symmetry block  1
    -0.038461    0.050692    0.093839    0.264858    0.404700    0.995049    1.004043    1.059757    1.188616    1.395675
     1.398793    1.615703    1.819752    1.843152    1.859661    2.472246    2.525236    2.938486    3.055402    3.460555
     3.952316    4.126008    4.439983    4.783645    4.892674    5.501172    5.511502

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   12 emc=    -94.1576226770 demc= 5.8251E-04 wnorm= 7.5130E-03 knorm= 1.2772E-01 apxde= 7.1806E-05    *not conv.*     

               starting mcscf iteration...  13

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3927279322     -133.1942757752        0.0000000000        0.0000010000
    2       -94.0730060683     -132.8745539114        0.0000000000        0.0000010000
    3       -94.0073852190     -132.8089330621        0.0000000000        0.0000010000
    4       -93.9049205178     -132.7064683609        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.583442118414107E-004
 Total number of micro iterations:    9

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99968590 pnorm= 0.0000E+00 rznorm= 5.0469E-07 rpnorm= 0.0000E+00 noldr=  9 nnewr=  9 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.797871  -23.155267   -3.103833   -2.350709   -1.928607

 qvv(*) eigenvalues. symmetry block  1
    -0.038375    0.050799    0.093735    0.264829    0.404602    0.995148    1.003842    1.059798    1.188435    1.395848
     1.399095    1.615519    1.819796    1.843625    1.859643    2.472468    2.525352    2.938623    3.055432    3.460552
     3.952307    4.125914    4.439965    4.783716    4.892714    5.501299    5.511649

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   13 emc=    -94.1577064065 demc= 8.3730E-05 wnorm= 2.0668E-03 knorm= 2.5062E-02 apxde= 2.4406E-06    *not conv.*     

               starting mcscf iteration...  14

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3927483558     -133.1942961988        0.0000000000        0.0000010000
    2       -94.0729943029     -132.8745421460        0.0000000000        0.0000010000
    3       -94.0073843334     -132.8089321764        0.0000000000        0.0000010000
    4       -93.9047896507     -132.7063374938        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.585964332266391E-005
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99999923 pnorm= 0.0000E+00 rznorm= 9.0565E-07 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.797207  -23.155791   -3.103678   -2.350421   -1.928324

 qvv(*) eigenvalues. symmetry block  1
    -0.038342    0.050825    0.093687    0.264814    0.404564    0.995160    1.003760    1.059780    1.188355    1.395908
     1.399229    1.615427    1.819793    1.843784    1.859635    2.472551    2.525379    2.938641    3.055436    3.460532
     3.952277    4.125863    4.439942    4.783733    4.892713    5.501332    5.511685

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   14 emc=    -94.1577089974 demc= 2.5909E-06 wnorm= 3.6688E-04 knorm= 1.2373E-03 apxde= 1.4171E-08    *not conv.*     

               starting mcscf iteration...  15

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   4 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3927519129     -133.1942997560        0.0000000000        0.0000010000
    2       -94.0729929451     -132.8745407882        0.0000000000        0.0000010000
    3       -94.0073821816     -132.8089300247        0.0000000000        0.0000010000
    4       -93.9047974753     -132.7063453184        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.905620097030913E-006
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    4

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 2.2165E-04 rznorm= 6.6292E-07 rpnorm= 5.7368E-07 noldr=  4 nnewr=  4 nolds=  3 nnews=  3
 

 fdd(*) eigenvalues. symmetry block  1
   -31.797086  -23.155887   -3.103657   -2.350405   -1.928315

 qvv(*) eigenvalues. symmetry block  1
    -0.038336    0.050828    0.093680    0.264811    0.404558    0.995161    1.003747    1.059775    1.188342    1.395918
     1.399255    1.615411    1.819792    1.843811    1.859635    2.472567    2.525383    2.938642    3.055439    3.460530
     3.952271    4.125855    4.439938    4.783737    4.892713    5.501339    5.511690

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=   15 emc=    -94.1577090132 demc= 1.5853E-08 wnorm= 6.3245E-05 knorm= 1.7766E-05 apxde= 3.4383E-10    *not conv.*     

               starting mcscf iteration...  16

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

 module tranlib: workspace lcore= 131066078

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 130896698
 address segment size,           sizesg = 130745212
 number of in-core blocks,       nincbk =         1
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      58626 transformed 1/r12    array elements were written in      11 records.


 mosort: allocated sort2 space, avc2is=   130927931 available sort2 space, avcisx=   130928183

   5 trial vectors read from nvfile (unit= 29).
 ciiter=   3 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*      -94.3927526565     -133.1943004996        0.0000000000        0.0000010000
    2       -94.0729927700     -132.8745406131        0.0000000000        0.0000010000
    3       -94.0073816141     -132.8089294572        0.0000000000        0.0000010000
    4       -93.9048022431     -132.7063500862        0.0000000000        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  9.570920924216618E-008
 performing all-state projection
 Total number of micro iterations:    1

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 1.00000000 pnorm= 0.0000E+00 rznorm= 9.3189E-07 rpnorm= 5.3475E-09 noldr=  1 nnewr=  1 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -31.797061  -23.155907   -3.103653   -2.350405   -1.928316

 qvv(*) eigenvalues. symmetry block  1
    -0.038334    0.050828    0.093678    0.264810    0.404556    0.995161    1.003745    1.059774    1.188339    1.395921
     1.399260    1.615408    1.819792    1.843816    1.859636    2.472571    2.525384    2.938643    3.055440    3.460529
     3.952270    4.125853    4.439938    4.783738    4.892713    5.501340    5.511692

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=   16 emc=    -94.1577090136 demc= 3.4461E-10 wnorm= 7.6567E-07 knorm= 3.6972E-08 apxde= 1.3425E-14    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.333 total energy=      -94.392752657, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.333 total energy=      -94.072992770, rel. (eV)=   8.701113
   DRT #1 state # 3 wt 0.333 total energy=      -94.007381614, rel. (eV)=  10.486484
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A  
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99195941     1.67602946     1.62103743
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.71097371     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22        MO   23        MO   24
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   25        MO   26        MO   27        MO   28        MO   29        MO   30        MO   31        MO   32
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   33        MO   34        MO   35        MO   36
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      
 Computing the requested mcscf (transition) density matrices (flag 30)
 Reading mcdenin ...
 Number of density matrices (ndens):                     4
 Number of unique bra states (ndbra):                     3
 qind: F
 (Transition) density matrices:
 d1(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st01                                           
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st02-st01                                      
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st03-st01                                      
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        180 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st03-st02                                      

          state spec. NOs: DRT 1, State  1
 *** warning *** large active-orbital occupation. i=  3 nocc= 1.9993E+00
 *** warning *** large active-orbital occupation. i=  4 nocc= 2.0000E+00

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99996128     1.99929619     1.94884220
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.05190033     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
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
     1_ s       0.000547   1.999414   0.315470   0.809338   0.163258   0.000868
     1_ p       0.000483   0.000003   0.098467   0.031602   0.700305   0.234679
     1_ d      -0.000081   0.000001   0.010161   0.003579   0.014472   0.002475
     2_ s       1.998441  -0.000039   1.404151   0.206273  -0.001678   0.000070
     2_ p       0.000000  -0.000270   0.040259   0.499985   0.646765   1.227564
     2_ d       0.000000  -0.000104   0.004499   0.006448   0.004271   0.002321
     3_ s      -0.000002   0.000486   0.008365   0.130577   0.183605   0.041658
     3_ s      -0.000006   0.000498   0.009398   0.122224   0.160103   0.044447
     3_ s       0.000305   0.000006   0.057582   0.081252   0.052039   0.237226
     3_ s       0.000312   0.000005   0.051647   0.108723   0.076860   0.200650
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000015   0.000031   0.000322   0.000000   0.000000   0.000000
     1_ p       0.834299   0.193956   0.626512   0.000000   0.000000   0.000000
     1_ d       0.012727   0.015493   0.001872   0.000000   0.000000   0.000000
     2_ s      -0.000024   0.000144   0.000011   0.000000   0.000000   0.000000
     2_ p       0.224850   1.335757   0.072007   0.000000   0.000000   0.000000
     2_ d       0.010430   0.003533   0.004835   0.000000   0.000000   0.000000
     3_ s       0.236026   0.024028   0.000457   0.000000   0.000000   0.000000
     3_ s       0.266304   0.025745   0.001072   0.000000   0.000000   0.000000
     3_ s       0.048992   0.012700   0.001823   0.000000   0.000000   0.000000
     3_ s       0.042410   0.009651   0.002063   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.289263   3.607349   0.625201   0.629785   0.491925   0.492321
      p         2.720306   4.046918   0.000000   0.000000   0.000000   0.000000
      d         0.060700   0.036233   0.000000   0.000000   0.000000   0.000000
    total       6.070269   7.690500   0.625201   0.629785   0.491925   0.492321
 

 Total number of electrons:   16.00000000

 Mulliken population for:
DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.000547   1.999414   0.315464   0.809337   0.163265   0.000716
     1_ p       0.000483   0.000003   0.098467   0.031599   0.700308   0.044385
     1_ d      -0.000081   0.000001   0.010161   0.003579   0.014472   0.005651
     2_ s       1.998441  -0.000039   1.404157   0.206267  -0.001678   0.000037
     2_ p       0.000000  -0.000270   0.040257   0.500000   0.646752   1.394435
     2_ d       0.000000  -0.000104   0.004499   0.006448   0.004271   0.005714
     3_ s      -0.000002   0.000486   0.008365   0.130574   0.183609   0.006275
     3_ s      -0.000006   0.000498   0.009398   0.122221   0.160106   0.005651
     3_ s       0.000305   0.000006   0.057583   0.081252   0.052038   0.291104
     3_ s       0.000312   0.000005   0.051648   0.108724   0.076857   0.245994
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000081   0.000247   0.000020   0.000000   0.000000   0.000000
     1_ p       1.191837   0.572048   0.036569   0.000000   0.000000   0.000000
     1_ d       0.014559   0.011947   0.000249   0.000000   0.000000   0.000000
     2_ s       0.000007   0.000147   0.000001   0.000000   0.000000   0.000000
     2_ p       0.009534   1.357809   0.014375   0.000000   0.000000   0.000000
     2_ d       0.011073   0.005524   0.000267   0.000000   0.000000   0.000000
     3_ s       0.343402   0.001156   0.000089   0.000000   0.000000   0.000000
     3_ s       0.384521   0.000706   0.000160   0.000000   0.000000   0.000000
     3_ s       0.023869  -0.000419   0.000075   0.000000   0.000000   0.000000
     3_ s       0.020413  -0.000323   0.000096   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         3.289091   3.607339   0.673953   0.683255   0.505813   0.503726
      p         2.675699   3.962894   0.000000   0.000000   0.000000   0.000000
      d         0.060539   0.037693   0.000000   0.000000   0.000000   0.000000
    total       6.025328   7.607925   0.673953   0.683255   0.505813   0.503726
 

 Total number of electrons:   16.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 02 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
 
   ao class       7A         8A         9A        10A        11A        12A  
 
   ao class      13A        14A        15A        16A        17A        18A  
 
   ao class      19A        20A        21A        22A        23A        24A  
 
   ao class      25A        26A        27A        28A        29A        30A  
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s       0.000000  -0.000000  -0.000006  -0.000001   0.000001   0.000000
     1_ p       0.000000  -0.000003  -0.000000  -0.000040   0.000000  -0.000000
     1_ d      -0.000000  -0.000000  -0.000004  -0.000097   0.000000   0.000225
     2_ s       0.000000  -0.000000  -0.000002  -0.000028   0.000001   0.000000
     2_ p      -0.000000  -0.000001  -0.000001  -0.000016   0.000002   0.000002
     2_ d      -0.000000  -0.000002  -0.000000  -0.000003   0.000000   0.000043
     3_ s      -0.000000  -0.000000  -0.000005   0.000000   0.000000  -0.000000
     3_ s      -0.000000  -0.000000  -0.000004  -0.000000   0.000001   0.000012
     3_ s      -0.000000  -0.000000   0.000000   0.000001   0.000001   0.000000
     3_ s      -0.000000  -0.000001  -0.000001  -0.000001  -0.000000   0.000002


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s        -0.000005  -0.000029  -0.000005   0.000008   0.000002  -0.000001
      p        -0.000043  -0.000014   0.000000   0.000000   0.000000   0.000000
      d         0.000123   0.000037   0.000000   0.000000   0.000000   0.000000
    total       0.000075  -0.000007  -0.000005   0.000008   0.000002  -0.000001
 

 Total number of electrons:    0.00007243

 Off-diagonal Mulliken population for:
DRT 1,state 03 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s      -0.000022   0.000000  -0.000000  -0.000000  -0.000000   0.000000
     1_ p      -0.000008   0.000000  -0.000005  -0.000002   0.000000   0.000022
     1_ d      -0.000001   0.000000  -0.000018  -0.000000  -0.000001   0.000002
     2_ s      -0.000005  -0.000000  -0.000000  -0.000000  -0.000000   0.000008
     2_ p      -0.000034   0.000000  -0.000002  -0.000002  -0.000000   0.000010
     2_ d      -0.000036   0.000000  -0.000001  -0.000000  -0.000001   0.000033
     3_ s      -0.000001   0.000001  -0.000000  -0.000000  -0.000000   0.000001
     3_ s       0.000000   0.000002   0.000000  -0.000000  -0.000000  -0.000002
     3_ s       0.000001   0.000000  -0.000000  -0.000001  -0.000001  -0.000002
     3_ s       0.000000   0.000001  -0.000000  -0.000000  -0.000000   0.000012
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000010   0.000000  -0.000000   0.000015  -0.000000  -0.000000
     1_ p       0.000002  -0.000000  -0.000000   0.000009   0.000000  -0.000001
     1_ d       0.000030  -0.000000  -0.000000   0.000001   0.000000  -0.000029
     2_ s       0.000006  -0.000000  -0.000000  -0.000000   0.000000  -0.000001
     2_ p       0.000012  -0.000000  -0.000001   0.000000   0.000000  -0.000006
     2_ d       0.000002  -0.000000  -0.000000   0.000000   0.000000  -0.000001
     3_ s      -0.000000  -0.000000  -0.000000  -0.000000  -0.000000  -0.000000
     3_ s      -0.000000   0.000000  -0.000000   0.000000   0.000000   0.000000
     3_ s      -0.000000  -0.000000  -0.000001   0.000000  -0.000000  -0.000000
     3_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000000
 
   ao class      13A        14A        15A        16A        17A        18A  
     1_ s       0.000000  -0.000000   0.000000  -0.000001  -0.000005   0.000000
     1_ p       0.000000  -0.000000   0.000000  -0.000000  -0.000009   0.000000
     1_ d       0.000000  -0.000001  -0.000000  -0.000000  -0.000022   0.000000
     2_ s       0.000000  -0.000000  -0.000000  -0.000000  -0.000008   0.000000
     2_ p       0.000000  -0.000000  -0.000000  -0.000000  -0.000009   0.000000
     2_ d       0.000001  -0.000000  -0.000000  -0.000000  -0.000025   0.000000
     3_ s       0.000000   0.000000  -0.000000  -0.000000   0.000000   0.000000
     3_ s       0.000000  -0.000000  -0.000000  -0.000000   0.000000  -0.000000
     3_ s       0.000000  -0.000000  -0.000000  -0.000000   0.000000   0.000000
     3_ s       0.000000  -0.000000  -0.000000  -0.000000  -0.000001   0.000000
 
   ao class      19A        20A        21A        22A        23A        24A  
     1_ s       0.000000   0.000001  -0.000000  -0.000001  -0.000000   0.000000
     1_ p       0.000000   0.000028  -0.000000  -0.000000  -0.000002   0.000000
     1_ d       0.000000   0.000009  -0.000000  -0.000001  -0.000017   0.000001
     2_ s       0.000000   0.000013  -0.000000  -0.000000  -0.000002   0.000000
     2_ p       0.000000   0.000031  -0.000000  -0.000001  -0.000001   0.000000
     2_ d       0.000000   0.000111  -0.000000  -0.000000  -0.000001   0.000000
     3_ s       0.000000   0.000009  -0.000000  -0.000001  -0.000000   0.000000
     3_ s       0.000001   0.000001   0.000000  -0.000021   0.000000   0.000000
     3_ s       0.000000   0.000033  -0.000000  -0.000007   0.000000  -0.000000
     3_ s       0.000000   0.000004  -0.000000   0.000000  -0.000001  -0.000000
 
   ao class      25A        26A        27A        28A        29A        30A  
     1_ s       0.000000  -0.000002   0.000000   0.000003  -0.000002   0.000001
     1_ p      -0.000000  -0.000000   0.000000   0.000053  -0.000000   0.000014
     1_ d       0.000000  -0.000001  -0.000000   0.000047  -0.000003   0.000000
     2_ s       0.000000   0.000000  -0.000000   0.000057   0.000000   0.000000
     2_ p       0.000004  -0.000000  -0.000000   0.000003  -0.000001   0.000000
     2_ d       0.000025  -0.000000  -0.000000   0.000006  -0.000001  -0.000000
     3_ s       0.000002  -0.000000  -0.000000   0.000002  -0.000000   0.000000
     3_ s       0.000000  -0.000000  -0.000005  -0.000000  -0.000002  -0.000000
     3_ s       0.000000  -0.000000  -0.000005   0.000008  -0.000002   0.000000
     3_ s       0.000001   0.000000   0.000000   0.000003  -0.000000   0.000000
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s      -0.000000  -0.000000  -0.000000  -0.000000   0.000032  -0.000000
     1_ p      -0.000000  -0.000000  -0.000000   0.000000   0.000005  -0.000000
     1_ d      -0.000000  -0.000001  -0.000000  -0.000000   0.000000  -0.000002
     2_ s      -0.000000  -0.000000  -0.000000   0.000000   0.000000  -0.000000
     2_ p      -0.000000  -0.000001  -0.000000  -0.000000   0.000002  -0.000000
     2_ d      -0.000001  -0.000000  -0.000000  -0.000000   0.000002  -0.000000
     3_ s      -0.000000   0.000000  -0.000000  -0.000001  -0.000001  -0.000000
     3_ s      -0.000000   0.000000  -0.000000  -0.000000   0.000001  -0.000000
     3_ s      -0.000000  -0.000000  -0.000000   0.000000   0.000001  -0.000000
     3_ s       0.000000   0.000000  -0.000000  -0.000000   0.000002   0.000000


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         0.000028   0.000066   0.000009  -0.000027   0.000025   0.000023
      p         0.000105   0.000005   0.000000   0.000000   0.000000   0.000000
      d        -0.000006   0.000110   0.000000   0.000000   0.000000   0.000000
    total       0.000128   0.000181   0.000009  -0.000027   0.000025   0.000023
 

 Total number of electrons:    0.00033891

 Off-diagonal Mulliken population for:
DRT 1,state 03 - DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A   partial gross atomic populations
   ao class       1A         2A         3A         4A         5A         6A  
     1_ s       0.000000   0.000000  -0.000000  -0.000054  -0.100451  -0.000124
     1_ p       0.000000   0.000000  -0.000000   0.000027  -0.017111  -0.013043
     1_ d       0.000000   0.000000  -0.000000  -0.000029  -0.015019  -0.079271
     2_ s       0.000000  -0.000000  -0.000000  -0.000004  -0.000629  -0.000000
     2_ p       0.000000  -0.000000  -0.000000   0.261641   0.000210   0.000000
     2_ d       0.000000  -0.000000  -0.000000   1.246223  -0.001835  -0.073159
     3_ s      -0.000000  -0.000000  -0.000000   0.045072  -0.045390  -0.005443
     3_ s       0.000000  -0.000000  -0.000000  -0.000020  -0.029432  -0.003053
     3_ s       0.000000   0.000000   0.000000  -0.000582  -0.000025  -0.000317
     3_ s       0.000000  -0.000000   0.000000   0.360945  -0.010963  -0.000000
 
   ao class       7A         8A         9A        10A        11A        12A  
     1_ s       0.000000  -0.000001   0.109529   0.000000   0.000000  -0.000008
     1_ p       0.000000   0.000000   0.593572  -0.000000   0.000000  -0.000000
     1_ d      -0.000000   0.000003   2.686737   0.008547  -0.068334   0.000000
     2_ s      -0.000000   0.000002   0.178221   0.000000  -0.000517   0.000000
     2_ p       0.596942  -0.055546   0.148099  -0.000000   0.000000  -0.033480
     2_ d       3.252927  -2.744769   0.200471  -0.000000   0.000000  -0.000004
     3_ s       0.000000  -0.259837   0.000201   0.008717   0.000000  -0.000000
     3_ s      -0.000000  -0.503230   0.000127  -0.000001  -0.000000  -0.000000
     3_ s       0.000000  -1.130252  -0.000000  -0.000000  -0.071788   0.000021
     3_ s      -0.000000  -2.464286  -0.000000   0.000000  -0.000031  -0.003295
 
   ao class      13A        14A        15A        16A        17A        18A  
     1_ s      -0.012043   0.008255   0.000528  -0.002692   0.000276   0.000000
     1_ p      -0.001884   0.023163   0.010737   0.000012   0.000802   0.000000
     1_ d      -0.205984   0.021273   0.000296  -0.002067   0.000000  -0.000000
     2_ s       0.005186   0.000231   0.000085  -0.000552   0.000000  -0.263782
     2_ p      -0.000798  -0.000088   0.000409  -0.000571   0.000000  -0.000000
     2_ d      -0.000474  -0.000231   0.000332  -0.219689   0.000000  -0.000000
     3_ s      -0.002089  -0.000002   0.002748  -0.008924   0.000000  -0.000000
     3_ s      -0.003086  -0.000002   0.004359  -0.009423   0.000000   0.000000
     3_ s      -0.011525  -0.000003   0.003488  -0.003085   0.000000   0.000000
     3_ s      -0.005233   0.043767   0.004699  -0.000233  -0.000000  -0.000000
 
   ao class      19A        20A        21A        22A        23A        24A  
     1_ s      -0.000000   0.000000   0.000001  -0.000329  -0.000001   0.000001
     1_ p      -0.000000   0.000000   0.000001  -0.001130  -0.000001   0.000004
     1_ d      -2.211219   0.000000   0.000002  -0.000870  -0.000003   0.000001
     2_ s       0.000038   0.000000   0.000002  -0.000047  -0.000000   0.000000
     2_ p       0.000445   0.000000   0.000001  -0.000377  -0.000000  -0.000000
     2_ d      -0.001095   0.000002   0.000012  -0.000672  -0.000001   0.000026
     3_ s      -0.000004   0.000000   0.000000  -0.000094   0.000000   0.000000
     3_ s       0.000011   0.000000   0.000000  -0.000042   0.000000   0.000000
     3_ s       0.000003   0.000000   0.000000  -0.000643  -0.000000   0.000000
     3_ s       0.000016  -0.000000  -0.000000  -0.000003  -0.000000   0.000005
 
   ao class      25A        26A        27A        28A        29A        30A  
     1_ s       0.000001   0.000000   0.000001  -0.000000   0.000000   0.000000
     1_ p       0.000048   0.000049   0.000000  -0.000000  -0.000000   0.000000
     1_ d       0.000016   0.000008   0.000000  -0.000000  -0.000000   0.000010
     2_ s       0.000001  -0.000001   0.000000  -0.000000  -0.000000   0.000001
     2_ p       0.000041   0.000001   0.000000  -0.000000  -0.000000   0.000001
     2_ d       0.000008   0.000004   0.000000  -0.000001  -0.000000   0.000011
     3_ s      -0.000004  -0.000001   0.000000  -0.000000  -0.000000   0.000000
     3_ s       0.000001   0.000003   0.000000   0.000000  -0.000000   0.000000
     3_ s       0.000002   0.000006   0.000000  -0.000000  -0.000001   0.000000
     3_ s      -0.000002   0.000001   0.000000  -0.000000  -0.000000   0.000000
 
   ao class      31A        32A        33A        34A        35A        36A  
     1_ s       0.000000   0.000000  -0.000000   0.000010   0.000000  -0.000000
     1_ p      -0.000000  -0.000000  -0.000001   0.000004  -0.000001   0.000000
     1_ d      -0.000000  -0.000000  -0.000002   0.000019  -0.000002   0.000000
     2_ s      -0.000000  -0.000000  -0.000000   0.000001  -0.000000   0.000000
     2_ p      -0.000001   0.000000  -0.000001   0.000008   0.000000   0.000000
     2_ d      -0.000000  -0.000001  -0.000004   0.000018  -0.000001   0.000000
     3_ s       0.000000  -0.000000  -0.000000   0.000000  -0.000000   0.000000
     3_ s      -0.000000   0.000000  -0.000000   0.000000   0.000000   0.000000
     3_ s       0.000000  -0.000003  -0.000000   0.000000  -0.000000   0.000000
     3_ s      -0.000000  -0.000001   0.000000   0.000003  -0.000002   0.000000


                        gross atomic populations
     ao            1_         2_         3_         3_         3_         3_
      s         0.002899  -0.081765  -0.265052  -0.543787  -1.214703  -2.074613
      p         0.595246   0.916936   0.000000   0.000000   0.000000   0.000000
      d         0.134111   1.658098   0.000000   0.000000   0.000000   0.000000
    total       0.732256   2.493269  -0.265052  -0.543787  -1.214703  -2.074613
 

 Total number of electrons:   -0.87262960

 !timer: mcscf                           cpu_time=     0.466 walltime=     0.465
