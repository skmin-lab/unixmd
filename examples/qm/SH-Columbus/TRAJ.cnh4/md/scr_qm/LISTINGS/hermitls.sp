
 Work memory size (LMWORK) :   131072000 =1000   megabytes.

 Default basis set library used :
        /sphome/kedziora/dalton/basis/                              


    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    $$$$$$$$$$$  DALTON - An electronic structure program  $$$$$$$$$$$
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

               This is output from DALTON (beta-version 0.9) 

                          Principal authors:

            Trygve Helgaker,     University of Oslo,        Norway 
            Hans Joergen Jensen, University of Odense,      Denmark
            Poul Joergensen,     University of Aarhus,      Denmark
            Henrik Koch,         University of Aarhus,      Denmark
            Jeppe Olsen,         University of Lund,        Sweden 
            Hans Aagren,         University of Linkoeping,  Sweden 

                          Contributors:

            Torgeir Andersen,    University of Oslo,        Norway 
            Keld L. Bak,         University of Copenhagen,  Denmark
            Vebjoern Bakken,     University of Oslo,        Norway 
            Ove Christiansen,    University of Aarhus,      Denmark
            Paal Dahle,          University of Oslo,        Norway 
            Erik K. Dalskov,     University of Odense,      Denmark
            Thomas Enevoldsen,   University of Odense,      Denmark
            Asger Halkier,       University of Aarhus,      Denmark
            Hanne Heiberg,       University of Oslo,        Norway 
            Dan Jonsson,         University of Linkoeping,  Sweden 
            Sheela Kirpekar,     University of Odense,      Denmark
            Rika Kobayashi,      University of Aarhus,      Denmark
            Alfredo S. de Meras, Valencia University,       Spain  
            Kurt Mikkelsen,      University of Aarhus,      Denmark
            Patrick Norman,      University of Linkoeping,  Sweden 
            Martin J. Packer,    University of Sheffield,   UK     
            Kenneth Ruud,        University of Oslo,        Norway 
            Trond Saue,          University of Oslo,        Norway 
            Peter Taylor,        San Diego Superc. Center,  USA    
            Olav Vahtras,        University of Linkoeping,  Sweden

                                             Release Date:  August 1996
------------------------------------------------------------------------


      
     NOTE:
      
     This is an experimental code for the evaluation of molecular
     properties using (MC)SCF/CC wave functions. The authors accept
      no responsibility for the performance of the code or for the
     correctness of the results.
      
     The code (in whole or part) is not to be reproduced for further
     distribution without the written permission of T. Helgaker,
     H. J. Aa. Jensen or P. Taylor.
      
     If results obtained with this code are published, an
     appropriate citation would be:
      
     T. Helgaker, H. J. Aa. Jensen, P.Joergensen, H. Koch,
     J. Olsen, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
     O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
     A. Halkier, H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi,
     A. S. de Meras, K. V. Mikkelsen, P. Norman, M. J. Packer,
     K. Ruud, T.Saue, P. R. Taylor, and O. Vahtras:
     DALTON, an electronic structure program"



     ******************************************
     **    PROGRAM:              DALTON      **
     **    PROGRAM VERSION:      5.4.0.0     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************



 <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>




 Default print level:        0

    Integral sections will be executed
    Starting in Integral Section -



 *************************************************************************
 ****************** Output from HERMIT input processing ******************
 *************************************************************************



 Default print level:        2


 Calculation of one- and two-electron Hamiltonian integrals.


 The following one-electron property integrals are calculated:

          - overlap integrals
          - Cartesian multipole moment integrals of orders 4 and lower
          - electronic angular momentum around the origin


 Changes of defaults for READIN:
 -------------------------------


 Maximum number of primitives per integral block :   10



 *************************************************************************
 ****************** Output from READIN input processing ******************
 *************************************************************************



  Title Cards
  -----------

                                                                          
                                                                          


                      SYMGRP:Point group information
                      ------------------------------

Point group: C1 

   * Character table

        |  E 
   -----+-----
    A   |   1

   * Direct product table

        | A  
   -----+-----
    A   | A  


  Atoms and basis sets
  --------------------

  Number of atom types:     3
  Total number of atoms:    6

  label    atoms   charge   prim    cont     basis   
  ----------------------------------------------------------------------
  C  1        1       6      27      14      [10s4p1d|3s2p1d]                             
  N  1        1       7      27      14      [10s4p1d|3s2p1d]                             
  H  1        1       1       4       2      [4s|2s]                                      
  H  2        1       1       4       2      [4s|2s]                                      
  H  3        1       1       4       2      [4s|2s]                                      
  H  4        1       1       4       2      [4s|2s]                                      
  ----------------------------------------------------------------------
  ----------------------------------------------------------------------
  total:      6      17      70      36

  Spherical harmonic basis used.
  Threshold for integrals:  1.00D-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates: 18


   1   C  1     x      2.8609140600
   2            y      4.2451482300
   3            z      4.2535066900

   4   N  1     x      5.2753030500
   5            y      4.0903880100
   6            z      4.0984674100

   7   H  1     x      1.9803792300
   8            y      5.4371881700
   9            z      5.6541816700

  10   H  2     x      1.7726110400
  11            y      3.2497601400
  12            z      2.8674415200

  13   H  3     x      6.3648569700
  14            y      5.3657881400
  15            z      5.0174784100

  16   H  4     x      6.0884235600
  17            y      2.6765319200
  18            z      3.1170715100



   Interatomic separations (in Angstroms):
   ---------------------------------------

            C  1        N  1        H  1        H  2        H  3        H  4

   C  1    0.000000
   N  1    1.282887    0.000000
   H  1    1.079079    2.055677    0.000000
   H  2    1.071028    2.014416    1.877940    0.000000
   H  3    1.988263    1.012148    2.344806    2.907536    0.000000
   H  4    1.991904    1.007286    2.943204    2.307670    1.748693    0.000000




  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2                           distance
                  ------     ------                           --------
  bond distance:    N  1       C  1                           1.282887
  bond distance:    H  1       C  1                           1.079079
  bond distance:    H  2       C  1                           1.071028
  bond distance:    H  3       N  1                           1.012148
  bond distance:    H  4       N  1                           1.007286


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3                   angle
                  ------     ------     ------                   -----
  bond angle:       H  1       C  1       N  1                 120.750
  bond angle:       H  2       C  1       N  1                 117.408
  bond angle:       H  2       C  1       H  1                 121.716
  bond angle:       H  3       N  1       C  1                 119.605
  bond angle:       H  4       N  1       C  1                 120.385
  bond angle:       H  4       N  1       H  3                 119.979


  Nuclear repulsion energy :   38.801547843067


  Orbital exponents and contraction coefficients
  ----------------------------------------------


  C  1   1s    1     3047.524900    0.0018  0.0000  0.0000
   seg. cont.  2      457.369510    0.0140  0.0000  0.0000
               3      103.948690    0.0688  0.0000  0.0000
               4       29.210155    0.2322  0.0000  0.0000
               5        9.286663    0.4679  0.0000  0.0000
               6        3.163927    0.3623  0.0000  0.0000
               7        7.868272    0.0000 -0.1193  0.0000
               8        1.881288    0.0000 -0.1609  0.0000
               9        0.544249    0.0000  1.1435  0.0000
              10        0.168714    0.0000  0.0000  1.0000

  C  1   2px  11        7.868272    0.0690  0.0000
   seg. cont. 12        1.881288    0.3164  0.0000
              13        0.544249    0.7443  0.0000
              14        0.168714    0.0000  1.0000

  C  1   2py  15        7.868272    0.0690  0.0000
   seg. cont. 16        1.881288    0.3164  0.0000
              17        0.544249    0.7443  0.0000
              18        0.168714    0.0000  1.0000

  C  1   2pz  19        7.868272    0.0690  0.0000
   seg. cont. 20        1.881288    0.3164  0.0000
              21        0.544249    0.7443  0.0000
              22        0.168714    0.0000  1.0000

  C  1   3d2- 23        0.800000    1.0000

  C  1   3d1- 24        0.800000    1.0000

  C  1   3d0  25        0.800000    1.0000

  C  1   3d1+ 26        0.800000    1.0000

  C  1   3d2+ 27        0.800000    1.0000

  N  1   1s   28     4173.511000    0.0018  0.0000  0.0000
   seg. cont. 29      627.457900    0.0140  0.0000  0.0000
              30      142.902100    0.0686  0.0000  0.0000
              31       40.234330    0.2322  0.0000  0.0000
              32       12.820210    0.4691  0.0000  0.0000
              33        4.390437    0.3605  0.0000  0.0000
              34       11.626358    0.0000 -0.1150  0.0000
              35        2.716280    0.0000 -0.1691  0.0000
              36        0.772218    0.0000  1.1459  0.0000
              37        0.212031    0.0000  0.0000  1.0000

  N  1   2px  38       11.626358    0.0676  0.0000
   seg. cont. 39        2.716280    0.3239  0.0000
              40        0.772218    0.7409  0.0000
              41        0.212031    0.0000  1.0000

  N  1   2py  42       11.626358    0.0676  0.0000
   seg. cont. 43        2.716280    0.3239  0.0000
              44        0.772218    0.7409  0.0000
              45        0.212031    0.0000  1.0000

  N  1   2pz  46       11.626358    0.0676  0.0000
   seg. cont. 47        2.716280    0.3239  0.0000
              48        0.772218    0.7409  0.0000
              49        0.212031    0.0000  1.0000

  N  1   3d2- 50        0.800000    1.0000

  N  1   3d1- 51        0.800000    1.0000

  N  1   3d0  52        0.800000    1.0000

  N  1   3d1+ 53        0.800000    1.0000

  N  1   3d2+ 54        0.800000    1.0000

  H  1   1s   55       18.731137    0.0335  0.0000
   seg. cont. 56        2.825394    0.2347  0.0000
              57        0.640122    0.8138  0.0000
              58        0.161278    0.0000  1.0000

  H  2   1s   59       18.731137    0.0335  0.0000
   seg. cont. 60        2.825394    0.2347  0.0000
              61        0.640122    0.8138  0.0000
              62        0.161278    0.0000  1.0000

  H  3   1s   63       18.731137    0.0335  0.0000
   seg. cont. 64        2.825394    0.2347  0.0000
              65        0.640122    0.8138  0.0000
              66        0.161278    0.0000  1.0000

  H  4   1s   67       18.731137    0.0335  0.0000
   seg. cont. 68        2.825394    0.2347  0.0000
              69        0.640122    0.8138  0.0000
              70        0.161278    0.0000  1.0000


  Contracted Orbitals
  -------------------

   1  C  1    1s       1     2     3     4     5     6
   2  C  1    1s       7     8     9
   3  C  1    1s      10
   4  C  1    2px     11    12    13
   5  C  1    2py     15    16    17
   6  C  1    2pz     19    20    21
   7  C  1    2px     14
   8  C  1    2py     18
   9  C  1    2pz     22
  10  C  1    3d2-    23
  11  C  1    3d1-    24
  12  C  1    3d0     25
  13  C  1    3d1+    26
  14  C  1    3d2+    27
  15  N  1    1s      28    29    30    31    32    33
  16  N  1    1s      34    35    36
  17  N  1    1s      37
  18  N  1    2px     38    39    40
  19  N  1    2py     42    43    44
  20  N  1    2pz     46    47    48
  21  N  1    2px     41
  22  N  1    2py     45
  23  N  1    2pz     49
  24  N  1    3d2-    50
  25  N  1    3d1-    51
  26  N  1    3d0     52
  27  N  1    3d1+    53
  28  N  1    3d2+    54
  29  H  1    1s      55    56    57
  30  H  1    1s      58
  31  H  2    1s      59    60    61
  32  H  2    1s      62
  33  H  3    1s      63    64    65
  34  H  3    1s      66
  35  H  4    1s      67    68    69
  36  H  4    1s      70




  Symmetry Orbitals
  -----------------

  Number of orbitals in each symmetry:        36


  Symmetry  A  ( 1)

    1     C  1     1s         1
    2     C  1     1s         2
    3     C  1     1s         3
    4     C  1     2px        4
    5     C  1     2py        5
    6     C  1     2pz        6
    7     C  1     2px        7
    8     C  1     2py        8
    9     C  1     2pz        9
   10     C  1     3d2-      10
   11     C  1     3d1-      11
   12     C  1     3d0       12
   13     C  1     3d1+      13
   14     C  1     3d2+      14
   15     N  1     1s        15
   16     N  1     1s        16
   17     N  1     1s        17
   18     N  1     2px       18
   19     N  1     2py       19
   20     N  1     2pz       20
   21     N  1     2px       21
   22     N  1     2py       22
   23     N  1     2pz       23
   24     N  1     3d2-      24
   25     N  1     3d1-      25
   26     N  1     3d0       26
   27     N  1     3d1+      27
   28     N  1     3d2+      28
   29     H  1     1s        29
   30     H  1     1s        30
   31     H  2     1s        31
   32     H  2     1s        32
   33     H  3     1s        33
   34     H  3     1s        34
   35     H  4     1s        35
   36     H  4     1s        36

  Symmetries of electric field:  A  (1)  A  (1)  A  (1)

  Symmetries of magnetic field:  A  (1)  A  (1)  A  (1)


 Copy of input to READIN
 -----------------------

INTGRL                                                                          
                                                                                
                                                                                
s   3    0           0.10D-14                                                   
       6.0    1    3    3    2    1                                             
C  1   2.860914060000000   4.245148230000000   4.253506690000000       *        
H   6   1                                                                       
       3047.52490000         0.00183470                                         
        457.36951000         0.01403730                                         
        103.94869000         0.06884260                                         
         29.21015500         0.23218440                                         
          9.28666300         0.46794130                                         
          3.16392700         0.36231200                                         
H   3   1                                                                       
          7.86827240        -0.11933240                                         
          1.88128850        -0.16085420                                         
          0.54424930         1.14345640                                         
H   1   1                                                                       
          0.16871440         1.00000000                                         
H   3   1                                                                       
          7.86827240         0.06899910                                         
          1.88128850         0.31642400                                         
          0.54424930         0.74430830                                         
H   1   1                                                                       
          0.16871440         1.00000000                                         
H   1   1                                                                       
          0.80000000         1.00000000                                         
       7.0    1    3    3    2    1                                             
N  1   5.275303050000000   4.090388010000000   4.098467410000000       *        
H   6   1                                                                       
       4173.51100000         0.00183480                                         
        627.45790000         0.01399500                                         
        142.90210000         0.06858700                                         
         40.23433000         0.23224100                                         
         12.82021000         0.46907000                                         
          4.39043700         0.36045500                                         
H   3   1                                                                       
         11.62635800        -0.11496100                                         
          2.71628000        -0.16911800                                         
          0.77221800         1.14585200                                         
H   1   1                                                                       
          0.21203130         1.00000000                                         
H   3   1                                                                       
         11.62635800         0.06758000                                         
          2.71628000         0.32390700                                         
          0.77221800         0.74089500                                         
H   1   1                                                                       
          0.21203130         1.00000000                                         
H   1   1                                                                       
          0.80000000         1.00000000                                         
       1.0    4    1    2                                                       
H  1   1.980379230000000   5.437188170000000   5.654181670000000       *        
H  2   1.772611040000000   3.249760140000000   2.867441520000000       *        
H  3   6.364856970000000   5.365788140000000   5.017478410000000       *        
H  4   6.088423560000000   2.676531920000000   3.117071510000000       *        
H   3   1                                                                       
         18.73113700         0.03349460                                         
          2.82539370         0.23472695                                         
          0.64012170         0.81375733                                         
H   1   1                                                                       
          0.16127780         1.00000000                                         




 ************************************************************************
 ************************** Output from HERONE **************************
 ************************************************************************

found     496 non-vanashing overlap integrals
found     666 non-vanashing nuclear attraction integrals
found     504 non-vanashing kinetic energy integrals






 found     530 non-vanashing integrals ( typea=  1 typeb=  0)
 found     533 non-vanashing integrals ( typea=  1 typeb=  1)
 found     525 non-vanashing integrals ( typea=  1 typeb=  2)


 found     544 non-vanashing integrals ( typea=  1 typeb=  3)
 found     581 non-vanashing integrals ( typea=  1 typeb=  4)
 found     574 non-vanashing integrals ( typea=  1 typeb=  5)
 found     550 non-vanashing integrals ( typea=  1 typeb=  6)
 found     582 non-vanashing integrals ( typea=  1 typeb=  7)
 found     536 non-vanashing integrals ( typea=  1 typeb=  8)


 found     551 non-vanashing integrals ( typea=  1 typeb=  9)
 found     593 non-vanashing integrals ( typea=  1 typeb= 10)
 found     592 non-vanashing integrals ( typea=  1 typeb= 11)
 found     594 non-vanashing integrals ( typea=  1 typeb= 12)
 found     651 non-vanashing integrals ( typea=  1 typeb= 13)
 found     589 non-vanashing integrals ( typea=  1 typeb= 14)
 found     558 non-vanashing integrals ( typea=  1 typeb= 15)
 found     603 non-vanashing integrals ( typea=  1 typeb= 16)
 found     598 non-vanashing integrals ( typea=  1 typeb= 17)
 found     549 non-vanashing integrals ( typea=  1 typeb= 18)


 found     556 non-vanashing integrals ( typea=  1 typeb= 19)
 found     599 non-vanashing integrals ( typea=  1 typeb= 20)
 found     596 non-vanashing integrals ( typea=  1 typeb= 21)
 found     598 non-vanashing integrals ( typea=  1 typeb= 22)
 found     664 non-vanashing integrals ( typea=  1 typeb= 23)
 found     596 non-vanashing integrals ( typea=  1 typeb= 24)
 found     599 non-vanashing integrals ( typea=  1 typeb= 25)
 found     665 non-vanashing integrals ( typea=  1 typeb= 26)
 found     656 non-vanashing integrals ( typea=  1 typeb= 27)
 found     589 non-vanashing integrals ( typea=  1 typeb= 28)
 found     565 non-vanashing integrals ( typea=  1 typeb= 29)
 found     610 non-vanashing integrals ( typea=  1 typeb= 30)
 found     610 non-vanashing integrals ( typea=  1 typeb= 31)
 found     607 non-vanashing integrals ( typea=  1 typeb= 32)
 found     556 non-vanashing integrals ( typea=  1 typeb= 33)


 found     542 non-vanashing integrals ( typea=  2 typeb=  6)
 found     533 non-vanashing integrals ( typea=  2 typeb=  7)
 found     546 non-vanashing integrals ( typea=  2 typeb=  8)




 ************************************************************************
 ************************** Output from TWOINT **************************
 ************************************************************************


 Number of two-electron integrals written:    212533 (95.7%)
 Kilobytes written:                             3408




 >>>> Total CPU  time used in HERMIT:   0.07 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds

- End of Integral Section
