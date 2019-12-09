
 Work memory size (LMWORK) :    65536000 =500   megabytes.

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

  Number of atom types:     2
  Total number of atoms:    6

  label    atoms   charge   prim    cont     basis   
  ----------------------------------------------------------------------
  C  1        1       6      27      14      [10s4p1d|3s2p1d]                             
  C  2        1       6      27      14      [10s4p1d|3s2p1d]                             
  H  1        1       1       4       2      [4s|2s]                                      
  H  2        1       1       4       2      [4s|2s]                                      
  H  3        1       1       4       2      [4s|2s]                                      
  H  4        1       1       4       2      [4s|2s]                                      
  ----------------------------------------------------------------------
  ----------------------------------------------------------------------
  total:      6      16      70      36

  Spherical harmonic basis used.
  Threshold for integrals:  1.00D-15


  Cartesian Coordinates
  ---------------------

  Total number of coordinates: 18


   1   C  1     x     -1.4165094100
   2            y      0.0261996700
   3            z      0.0532323000

   4   C  2     x      1.2752627400
   5            y      0.0914142300
   6            z     -0.0648086900

   7   H  1     x     -1.7192232000
   8            y      0.4351672400
   9            z     -1.9063976400

  10   H  2     x     -2.2276067600
  11            y     -1.8500448200
  12            z     -0.0135828800

  13   H  3     x      2.4993930200
  14            y      1.2465773600
  15            z      1.0236847600

  16   H  4     x      2.0770538800
  17            y     -1.1764990900
  18            z     -1.4247174600



   Interatomic separations (in Angstroms):
   ---------------------------------------

            C  1        C  2        H  1        H  2        H  3        H  4

   C  1    0.000000
   C  2    1.426211    0.000000
   H  1    1.071376    1.869169    0.000000
   H  2    1.082246    2.119483    1.593113    0.000000
   H  3    2.230429    1.060695    2.751746    3.040329    0.000000
   H  4    2.105819    1.071479    2.197277    2.423554    1.836507    0.000000




  Bond distances (angstroms):
  ---------------------------

                  atom 1     atom 2                           distance
                  ------     ------                           --------
  bond distance:    C  2       C  1                           1.426211
  bond distance:    H  1       C  1                           1.071376
  bond distance:    H  2       C  1                           1.082246
  bond distance:    H  3       C  2                           1.060695
  bond distance:    H  4       C  2                           1.071479


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3                   angle
                  ------     ------     ------                   -----
  bond angle:       H  1       C  1       C  2                  95.857
  bond angle:       H  2       C  1       C  2                 114.637
  bond angle:       H  2       C  1       H  1                  95.417
  bond angle:       H  3       C  2       C  1                 126.878
  bond angle:       H  4       C  2       C  1                 114.189
  bond angle:       H  4       C  2       H  3                 118.932


  Nuclear repulsion energy :   32.785084298900


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

  C  2   1s   28     3047.524900    0.0018  0.0000  0.0000
   seg. cont. 29      457.369510    0.0140  0.0000  0.0000
              30      103.948690    0.0688  0.0000  0.0000
              31       29.210155    0.2322  0.0000  0.0000
              32        9.286663    0.4679  0.0000  0.0000
              33        3.163927    0.3623  0.0000  0.0000
              34        7.868272    0.0000 -0.1193  0.0000
              35        1.881288    0.0000 -0.1609  0.0000
              36        0.544249    0.0000  1.1435  0.0000
              37        0.168714    0.0000  0.0000  1.0000

  C  2   2px  38        7.868272    0.0690  0.0000
   seg. cont. 39        1.881288    0.3164  0.0000
              40        0.544249    0.7443  0.0000
              41        0.168714    0.0000  1.0000

  C  2   2py  42        7.868272    0.0690  0.0000
   seg. cont. 43        1.881288    0.3164  0.0000
              44        0.544249    0.7443  0.0000
              45        0.168714    0.0000  1.0000

  C  2   2pz  46        7.868272    0.0690  0.0000
   seg. cont. 47        1.881288    0.3164  0.0000
              48        0.544249    0.7443  0.0000
              49        0.168714    0.0000  1.0000

  C  2   3d2- 50        0.800000    1.0000

  C  2   3d1- 51        0.800000    1.0000

  C  2   3d0  52        0.800000    1.0000

  C  2   3d1+ 53        0.800000    1.0000

  C  2   3d2+ 54        0.800000    1.0000

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
  15  C  2    1s      28    29    30    31    32    33
  16  C  2    1s      34    35    36
  17  C  2    1s      37
  18  C  2    2px     38    39    40
  19  C  2    2py     42    43    44
  20  C  2    2pz     46    47    48
  21  C  2    2px     41
  22  C  2    2py     45
  23  C  2    2pz     49
  24  C  2    3d2-    50
  25  C  2    3d1-    51
  26  C  2    3d0     52
  27  C  2    3d1+    53
  28  C  2    3d2+    54
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
   15     C  2     1s        15
   16     C  2     1s        16
   17     C  2     1s        17
   18     C  2     2px       18
   19     C  2     2py       19
   20     C  2     2pz       20
   21     C  2     2px       21
   22     C  2     2py       22
   23     C  2     2pz       23
   24     C  2     3d2-      24
   25     C  2     3d1-      25
   26     C  2     3d0       26
   27     C  2     3d1+      27
   28     C  2     3d2+      28
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
                                                                                
                                                                                
s   2    0           0.10D-14                                                   
       6.0    2    3    3    2    1                                             
C  1  -1.416509410000000   0.026199670000000   0.053232300000000       *        
C  2   1.275262740000000   0.091414230000000  -0.064808690000000       *        
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
       1.0    4    1    2                                                       
H  1  -1.719223200000000   0.435167240000000  -1.906397640000000       *        
H  2  -2.227606760000000  -1.850044820000000  -0.013582880000000       *        
H  3   2.499393020000000   1.246577360000000   1.023684760000000       *        
H  4   2.077053880000000  -1.176499090000000  -1.424717460000000       *        
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
found     498 non-vanashing kinetic energy integrals






 found     525 non-vanashing integrals ( typea=  1 typeb=  0)
 found     524 non-vanashing integrals ( typea=  1 typeb=  1)
 found     520 non-vanashing integrals ( typea=  1 typeb=  2)


 found     538 non-vanashing integrals ( typea=  1 typeb=  3)
 found     570 non-vanashing integrals ( typea=  1 typeb=  4)
 found     568 non-vanashing integrals ( typea=  1 typeb=  5)
 found     538 non-vanashing integrals ( typea=  1 typeb=  6)
 found     568 non-vanashing integrals ( typea=  1 typeb=  7)
 found     526 non-vanashing integrals ( typea=  1 typeb=  8)


 found     538 non-vanashing integrals ( typea=  1 typeb=  9)
 found     584 non-vanashing integrals ( typea=  1 typeb= 10)
 found     586 non-vanashing integrals ( typea=  1 typeb= 11)
 found     584 non-vanashing integrals ( typea=  1 typeb= 12)
 found     646 non-vanashing integrals ( typea=  1 typeb= 13)
 found     574 non-vanashing integrals ( typea=  1 typeb= 14)
 found     538 non-vanashing integrals ( typea=  1 typeb= 15)
 found     586 non-vanashing integrals ( typea=  1 typeb= 16)
 found     574 non-vanashing integrals ( typea=  1 typeb= 17)
 found     526 non-vanashing integrals ( typea=  1 typeb= 18)


 found     538 non-vanashing integrals ( typea=  1 typeb= 19)
 found     586 non-vanashing integrals ( typea=  1 typeb= 20)
 found     586 non-vanashing integrals ( typea=  1 typeb= 21)
 found     584 non-vanashing integrals ( typea=  1 typeb= 22)
 found     664 non-vanashing integrals ( typea=  1 typeb= 23)
 found     586 non-vanashing integrals ( typea=  1 typeb= 24)
 found     586 non-vanashing integrals ( typea=  1 typeb= 25)
 found     664 non-vanashing integrals ( typea=  1 typeb= 26)
 found     652 non-vanashing integrals ( typea=  1 typeb= 27)
 found     574 non-vanashing integrals ( typea=  1 typeb= 28)
 found     538 non-vanashing integrals ( typea=  1 typeb= 29)
 found     586 non-vanashing integrals ( typea=  1 typeb= 30)
 found     586 non-vanashing integrals ( typea=  1 typeb= 31)
 found     574 non-vanashing integrals ( typea=  1 typeb= 32)
 found     526 non-vanashing integrals ( typea=  1 typeb= 33)


 found     510 non-vanashing integrals ( typea=  2 typeb=  6)
 found     510 non-vanashing integrals ( typea=  2 typeb=  7)
 found     512 non-vanashing integrals ( typea=  2 typeb=  8)




 ************************************************************************
 ************************** Output from TWOINT **************************
 ************************************************************************


 Number of two-electron integrals written:    212533 (95.7%)
 Kilobytes written:                             3408




 >>>> Total CPU  time used in HERMIT:   0.08 seconds
 >>>> Total wall time used in HERMIT:   0.00 seconds

- End of Integral Section
