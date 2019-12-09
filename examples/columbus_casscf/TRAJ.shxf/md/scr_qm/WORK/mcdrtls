
 program "mcdrt 4.1 a3"

 distinct row table specification and csf
 selection for mcscf wavefunction optimization.

 programmed by: ron shepard

 version date: 17-oct-91


 This Version of Program mcdrt is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCDRT       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 expanded keystroke file:
 /home/swmoon/work/unixmd-fixcol/src/TRAJ.shxf/md/scr_qm/WORK/mcdrtky            
 
 input the spin multiplicity [  0]: spin multiplicity:    1    singlet 
 input the total number of electrons [  0]: nelt:     16
 input the number of irreps (1-8) [  0]: nsym:      1
 enter symmetry labels:(y,[n]) enter 1 labels (a4):
 enter symmetry label, default=   1
 input the molecular spatial symmetry (irrep 1:nsym) [  0]: spatial symmetry is irrep number:      1
 
 input the list of doubly-occupied orbitals (sym(i),rmo(i),i=1,ndot):
 number of doubly-occupied orbitals:      7
 number of inactive electrons:     14
 number of active electrons:      2
 level(*)        1   2   3   4   5   6   7
 symd(*)         1   1   1   1   1   1   1
 slabel(*)     a   a   a   a   a   a   a  
 doub(*)         1   2   3   4   5   6   7
 
 input the active orbitals (sym(i),rmo(i),i=1,nact):
 nact:      2
 level(*)        1   2
 syml(*)         1   1
 slabel(*)     a   a  
 modrt(*)        8   9
 input the minimum cumulative occupation for each active level:
  a   a  
    8   9
 input the maximum cumulative occupation for each active level:
  a   a  
    8   9
 slabel(*)     a   a  
 modrt(*)        8   9
 occmin(*)       0   2
 occmax(*)       2   2
 input the minimum b value for each active level:
  a   a  
    8   9
 input the maximum b value for each active level:
  a   a  
    8   9
 slabel(*)     a   a  
 modrt(*)        8   9
 bmin(*)         0   0
 bmax(*)         2   2
 input the step masks for each active level:
 modrt:smask=
   8:1111   9:1111
 input the number of vertices to be deleted [  0]: number of vertices to be removed (a priori):      0
 number of rows in the drt:      5
 are any arcs to be manually removed?(y,[n])
 nwalk=       3
 input the range of drt levels to print (l1,l2):
 levprt(*)       0   2

 level  0 through level  2 of the drt:

 row lev a b syml lab rmo  l0  l1  l2  l3 isym xbar   y0    y1    y2    xp     z

   5   2 1 0   1 a     9    0   0   0   0   1     1     0     0     0     3     0
 ........................................

   2   1 1 0   1 a     8    5   0   0   0   1     1     0     0     0     1     0

   3   1 0 1   1 a     8    0   0   5   0   1     1     1     1     0     1     1

   4   1 0 0   1 a     8    0   0   0   5   1     1     1     1     1     1     2
 ........................................

   1   0 0 0   0       0    4   3   0   2   1     3     2     1     1     1     0
 ........................................

 initial csf selection step:
 total number of walks in the drt, nwalk=       3
 keep all of these walks?(y,[n]) individual walks will be generated from the drt.
 apply orbital-group occupation restrictions?(y,[n]) apply reference occupation restrictions?(y,[n]) manually select individual walks?(y,[n])
 step-vector based csf selection complete.
        3 csfs selected from       3 total walks.

 beginning step-vector based csf selection.
 enter [step_vector/disposition] pairs:

 enter the active orbital step vector, (-1/ to end):

 step-vector based csf selection complete.
        3 csfs selected from       3 total walks.

 beginning numerical walk selection:
 enter positive walk numbers to add walks,
 negative walk numbers to delete walks, and zero to end.

 input walk number (0 to end) [  0]:
 final csf selection complete.
        3 csfs selected from       3 total walks.
  drt construction and csf selection complete.
 
 input a title card, default=mdrt2_title
  title                                                                         
  
 input a drt file name, default=mcdrtfl
 drt and indexing arrays written to file:
 /home/swmoon/work/unixmd-fixcol/src/TRAJ.shxf/md/scr_qm/WORK/mcdrtfl            
 
 write the drt file?([y],n) include step(*) vectors?([y],n) drt file is being written...


   List of selected configurations (step vectors)


   CSF#     1    3 0
   CSF#     2    1 2
   CSF#     3    0 3
