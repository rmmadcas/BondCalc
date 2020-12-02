IMPLICIT NONE
REAL (KIND=r):: a0, a1, a2, a3, a4, a5, a6, a7
REAL (KIND=r):: a8, a9, a10, a11, a12, a13, a14
INTEGER:: t1(8) 
CHARACTER(LEN=8) date1
CHARACTER(LEN=10) time1
CHARACTER(LEN=5) zona

CALL DATE_AND_TIME (VALUES=t1, DATE=date1, ZONE=zona, TIME=time1)

OPEN(30,file="Output_BondCalc.dat",STATUS='REPLACE', ACTION='WRITE')

WRITE(30,'(A28,I2.2,A1,I2.2,A1,I4,A4,I2.2,A1,I2.2,A1,I2.2)') "Program compiled and run on ",&
        t1(3),"/",t1(2),"/",t1(1)," at ", t1(5),":",t1(6),":",t1(7)

WRITE(30,*)" "
WRITE(30,*)" INPUT"

WRITE(30,*)"  ===================================================== "
WRITE(30,*)" "

WRITE(30,'(A40,I15)')" Number of components", number_of_pdb
SELECT CASE (type_of_molecules(1))
    CASE ("wate")
    WRITE(30,'(A40,A30)')" Component 1", "water"

    CASE ("alco")
    WRITE(30,'(A40,A30)')" Component 1", "alcohol"

    CASE ("diol")
    WRITE(30,'(A40,A30)')" Component 1", "diol"

    CASE ("ring")
    WRITE(30,'(A40,A30)')" Component 1", "aromatic ring"

    CASE ("ammo")
    WRITE(30,'(A40,A30)')" Component 1", "ammonia"

    CASE ("gene")
    WRITE(30,'(A40,A30)')" Component 1", "generic molecule"
END SELECT
   
IF (number_of_pdb.eq.2) THEN
  SELECT CASE (type_of_molecules(2))
    CASE ("wate")
    WRITE(30,'(A40,A30)')" Component 2", "water"

    CASE ("alco")
    WRITE(30,'(A40,A30)')" Component 2", "alcohol"

    CASE ("diol")
    WRITE(30,'(A40,A30)')" Component 2", "diol"

    CASE ("ammo")
    WRITE(30,'(A40,A30)')" Component 2", "ammonia"

    CASE ("ring")
    WRITE(30,'(A40,A30)')" Component 2", "aromatic ring"

    CASE ("gene")
    WRITE(30,'(A40,A30)')" Component 2", "generic molecule"
  END SELECT
END IF


WRITE(30,*)" "
WRITE(30,'(A40,I15)')" Number of movies", nmovies


WRITE(30,*)" "
WRITE(30,*)" RESULTS"

WRITE(30,*)"  ===================================================== "
WRITE(30,*)" "

WRITE(30,'(A40)')" Component 1"
WRITE(30,'(A80,A20,F8.3,A5,F8.3)')" Average number of bonds per molecule","", avgnHB1(nmovies),&
        " +/-", sqrt((disp1)/(nmovies-1)-avgnHB1(nmovies)*avgnHB1(nmovies))

WRITE(30,*)" "
OPEN(51,FILE="Ave_bonds1.dat",STATUS='OLD', ACTION='READ')
READ(51,*)          !comentarios
READ(51,*)

READ(51,*) a0, a1, a2, a3, a4, a5, a6

WRITE(30,'(A80,A20,7F8.3)')" Percentage of molecules with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 1", &
       "", a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
IF ((type_of_molecules(1).eq."wate").or.(type_of_molecules(1).eq."diol").or.(type_of_molecules(1).eq."alco")) THEN
    READ(51,*) a0, a1, a2, a3, a4, a5, a6
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of oxygen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 1", &
       "", a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
    READ(51,*) a0, a1, a2, a3, a4, a5, a6
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of hydrogen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 1", &
       "", a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
END IF
CLOSE (UNIT=51)

WRITE(30,*)" "
WRITE(30,*)" "

IF (number_of_pdb.eq.2) THEN
 WRITE(30,'(A40)')" Component 2"
 WRITE(30,'(A80,A20,F8.3,A5,F8.3)')" Average number of bonds per molecule","", avgnHB2(nmovies),&
    " +/-", sqrt((disp2)/(nmovies-1)-avgnHB2(nmovies)*avgnHB2(nmovies))

 OPEN(52,FILE="Ave_bonds2.dat",STATUS='OLD', ACTION='READ')
 READ(52,*)          !comentarios
 READ(52,*)

 READ(52,*) a0, a1, a2, a3, a4, a5, a6

 WRITE(30,'(A80,A20,7F8.3)')" Percentage of molecules with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 2", &
  "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
 IF ((type_of_molecules(2).eq."wate").or.(type_of_molecules(2).eq."diol").or.(type_of_molecules(2).eq."alco")) THEN
    READ(52,*) a0, a1, a2, a3, a4, a5, a6
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of oxygen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 2", &
     "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
    READ(52,*) a0, a1, a2, a3, a4, a5, a6
                              
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of hydrogen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in Component 2", &
     "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
 END IF
 CLOSE (UNIT=52)

 WRITE(30,*)" "
 WRITE(30,*)" "

 WRITE(30,'(A40)')" Mixture "
 WRITE(30,'(A80,A20,F8.3,A5,F8.3)')" Average number of bonds per molecule","", avgnHBm(nmovies),&
   " +/-", sqrt((dispm)/(nmovies-1)-avgnHBm(nmovies)*avgnHBm(nmovies))
                   
 OPEN(53,FILE="Ave_bondsMix.dat",STATUS='OLD', ACTION='READ')
 READ(53,*)          !comentarios
 READ(53,*)

 READ(53,*) a0, a1, a2, a3, a4, a5, a6
 WRITE(30,'(A80,A20,7F8.3)')" Percentage of molecules with 0, 1, 2, 3, 4, 5 and 6 bonds in the mixture", &
  "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
 IF ((type_of_molecules(1).eq."wate").or.(type_of_molecules(1).eq."diol").or.(type_of_molecules(1).eq."alco")) THEN
    READ(53,*) a0, a1, a2, a3, a4, a5, a6
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of oxygen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in the mixture", &
     "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
    READ(53,*) a0, a1, a2, a3, a4, a5, a6
    WRITE(30,'(A80,A20,7F8.3)')" Percentage of hydrogen atoms with 0, 1, 2, 3, 4, 5 and 6 bonds in the mixture", &
     "",a0*100, a1*100, a2*100, a3*100, a4*100, a5*100, a6*100
 END IF
 CLOSE(UNIT=53)
 WRITE(30,*)" "
END IF


IF (clusters.eq."y") THEN
 WRITE(30,*)" ===================================================== "
 WRITE(30,*)" "
 WRITE(30,*)" Clusters have been calculated"
 WRITE(30,*)" "
 WRITE(30,*)" ===================================================== "
END IF

IF (SDF.eq."y") THEN
 WRITE(30,*)" ===================================================== "
 WRITE(30,*)" "
 WRITE(30,*)" Spatial Distribution Function has been calculated"
 WRITE(30,*)" "
 WRITE(30,*)" ===================================================== "
END IF

IF (life_time.eq."y") THEN
 WRITE(30,*)" "
 WRITE(30,*)" Life Time Results"
 WRITE(30,*)" ===================================================== "

 WRITE(30,*)" These data have to be taken with care."
 WRITE(30,*)" Check that the showed parameters fit the experimental data found in the files ",&
  "'lifetime.dat' using the curve f(x)=a*exp(-x/b)"

 WRITE(30,*)" The parameter 'b' must be multiplied by the step between movies to obtain the lifetime."
 WRITE(30,*)" Component 1 "

END IF

CLOSE(UNIT=30)
