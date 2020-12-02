IMPLICIT NONE
CHARACTER(LEN=4):: mole1, mole2
INTEGER:: flag1, flag2

!*********************************pures 1************************************
!*****************************pura de agua*********************************
IF((mole1.eq."wate").and.(mole2.eq."wate").and.(flag1.eq.1).and.(flag2.eq.1)) THEN
   dOOmin=dOOw
   dOHmin=dOHw
   theta_min=pi/6

   CALL distance(molecule_water1(k,i)%O,molecule_water1(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_water1(k,i)%O,molecule_water1(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_water1(k,j)%O,molecule_water1(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water1(k,j)%O,molecule_water1(k,j)%H1,molecule_water1(k,j)%H2,molecule_water1(k,i)%O,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          END IF
        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_water1(k,i)%O,molecule_water1(k,j)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_water1(k,j)%O,molecule_water1(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water1(k,j)%O,molecule_water1(k,j)%H1,molecule_water1(k,j)%H2,molecule_water1(k,i)%O,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      CALL distance(molecule_water1(k,j)%O,molecule_water1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.3
        CALL distance(molecule_water1(k,i)%O,molecule_water1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.3
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water1(k,i)%O,molecule_water1(k,i)%H1,molecule_water1(k,i)%H2,molecule_water1(k,j)%O,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.3
      END IF !(dOH.lt.dOHmin)           !condition 2.3

      CALL distance(molecule_water1(k,j)%O,molecule_water1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.4
        CALL distance(molecule_water1(k,i)%O,molecule_water1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.4
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1
 
         
          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water1(k,i)%O,molecule_water1(k,i)%H1,molecule_water1(k,i)%H2,molecule_water1(k,j)%O,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.4
      END IF !(dOH.lt.dOHmin)           !condition 2.4

   END IF !dOO.LT.dOOmin							!condition 1
END IF

!*****************************pura de alcohol*******************************
IF((mole1.eq."alco").and.(mole2.eq."alco").and.(flag1.eq.1).and.(flag2.eq.1)) THEN
dOOmin=dOOal
dOHmin=dOHal
theta_min=pi/6

   CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1
      CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,j)%H,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_alcohol1(k,j)%O,molecule_alcohol1(k,j)%H,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_alcohol1(k,j)%O,molecule_alcohol1(k,j)%H,molecule_alcohol1(k,i)%O,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          END IF

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol1(k,j)%O,molecule_alcohol1(k,i)%H,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_alcohol1(k,j)%O,1)
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2
      
   END IF !dOO.LT.dOOmin							!condition 1
END IF

!*****************************pura de diol*******************************
IF((mole1.eq."diol").and.(mole2.eq."diol").and.(flag1.eq.1).and.(flag2.eq.1)) THEN
dOOmin=dOOd
dOHmin=dOHd
theta_min=pi/6

   CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.1

      CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol1(k,j)%O1,molecule_diol1(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,j)%O1,molecule_diol1(k,j)%H1,molecule_diol1(k,i)%O1,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          END IF

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,j)%O1,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_diol1(k,j)%O1,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

   END IF !dOO.LT.dOOmin							!condition 1.1
   
   CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.2

      CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol1(k,j)%O2,molecule_diol1(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,j)%O2,molecule_diol1(k,j)%H2,molecule_diol1(k,i)%O1,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          END IF

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,j)%O2,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_diol1(k,j)%O2,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

   END IF !dOO.LT.dOOmin							!condition 1.2
   
   CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.3

      CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol1(k,j)%O1,molecule_diol1(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,j)%O1,molecule_diol1(k,j)%H1,molecule_diol1(k,i)%O2,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "5"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "5"
          END IF

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,j)%O1,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,molecule_diol1(k,j)%O1,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "6"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "6"
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2
      
   END IF !dOO.LT.dOOmin							!condition 1.3
   
   CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.4

      CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol1(k,j)%O2,molecule_diol1(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,j)%O2,molecule_diol1(k,j)%H2,molecule_diol1(k,i)%O2,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "7"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "7"
          END IF

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,j)%O2,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolec1(k)=nHBpmolec1(k)+2

          nn1=nn1+1; d1=d1+1          
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1
          
          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "8"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "8"
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,molecule_diol1(k,j)%O2,1)
          END IF

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

   END IF !dOO.LT.dOOmin							!condition 1.4
   
END IF

!*****************************pura de los aromaticos*************************
IF((mole1.eq."ring").and.(mole2.eq."ring").and.(flag1.eq.1).and.(flag2.eq.1)) THEN 
   dCD=dFFr
   dangleCD=0.5
   dFE=dFEr
   dangleFE=0.86602540378
 
   CALL distance(molecule_ring1(k,i)%centro, molecule_ring1(k,j)%centro,dist)

   IF(dist.lt.dCD) THEN
      CALL vectorial_product(molecule_ring1(k,i)%v_nor(1),molecule_ring1(k,i)%v_nor(2),molecule_ring1(k,i)%v_nor(3),&
           molecule_ring1(k,j)%v_nor(1),molecule_ring1(k,j)%v_nor(2),molecule_ring1(k,j)%v_nor(3),modu)
           modu=ABS(modu)
      IF(modu.lt.dangleCD)THEN  
              nHBpmolec1(k)=nHBpmolec1(k)+2

              nn1=nn1+1; d1=d1+1              
              
              enl1(k,i)=enl1(k,i)+1
              enl1(k,j)=enl1(k,j)+1
                
              IF (life_time.eq."y") THEN
               matrix1(k,i,j)=matrix1(k,i,j)+1
              END IF  

              IF (SDF.eq."y") THEN
                CALL SDF_ring(molecule_ring1(k,j)%centro,molecule_ring1(k,j)%C1,molecule_ring1(k,j)%C2,molecule_ring1(k,j)%C3,&
                     molecule_ring1(k,j)%C4,molecule_ring1(k,j)%C5,molecule_ring1(k,j)%C6,molecule_ring1(k,i)%centro,1)
              END IF

              WRITE(17,*)  k, i, j, nHBpmolec1(k), "1"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          END IF
      END IF
   END IF

   IF(dist.lt.dFE) THEN
           CALL vectorial_product(molecule_ring1(k,i)%v_nor(1),molecule_ring1(k,i)%v_nor(2),molecule_ring1(k,i)%v_nor(3),&
                molecule_ring1(k,j)%v_nor(1),molecule_ring1(k,j)%v_nor(2),molecule_ring1(k,j)%v_nor(3),modu)
           modu=ABS(modu)
           IF(modu.gt.dangleFE)THEN
              nHBpmolec1(k)=nHBpmolec1(k)+2

              nn1=nn1+1; d1=d1+1
              
              enl1(k,i)=enl1(k,i)+1
              enl1(k,j)=enl1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

              IF (SDF.eq."y") THEN
                CALL SDF_ring(molecule_ring1(k,j)%centro,molecule_ring1(k,j)%C1,molecule_ring1(k,j)%C2,molecule_ring1(k,j)%C3,&
                     molecule_ring1(k,j)%C4,molecule_ring1(k,j)%C5,molecule_ring1(k,j)%C6,molecule_ring1(k,i)%centro,1)
              END IF

              WRITE(17,*)  k, i, j, nHBpmolec1(k), "2"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          END IF
           END IF
   END IF
END IF

!*****************************pura de los amoniaco*************************
IF((mole1.eq."ammo").and.(mole2.eq."ammo").and.(flag1.eq.1).and.(flag2.eq.1)) THEN
   dOOmin=dNN
   dOHmin=dNH

   CALL distance(molecule_ammo1(k,i)%N,molecule_ammo1(k,j)%N,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo1(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,j)%N,molecule_ammo1(k,j)%H1,molecule_ammo1(k,j)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_ammo1(k,i)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "1"
          END IF
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo1(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,j)%N,molecule_ammo1(k,j)%H1,molecule_ammo1(k,j)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_ammo1(k,i)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "2"
          END IF
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo1(k,j)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,i)=enl1(k,i)+1
          enl1(k,j)=enl1(k,j)+1

          enlO1(k,i)=enlO1(k,i)+1
          enlH1(k,j)=enlH1(k,j)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,j)%N,molecule_ammo1(k,j)%H1,molecule_ammo1(k,j)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_ammo1(k,i)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "3"
          END IF
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo1(k,i)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,j)=enl1(k,j)+1
          enl1(k,i)=enl1(k,i)+1

          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_ammo1(k,j)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "4"
          END IF
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo1(k,i)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,j)=enl1(k,j)+1
          enl1(k,i)=enl1(k,i)+1

          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_ammo1(k,j)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "5"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "5"
          END IF
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo1(k,i)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolec1(k)=nHBpmolec1(k)+2
          
          nn1=nn1+1; d1=d1+1
          
          enl1(k,j)=enl1(k,j)+1
          enl1(k,i)=enl1(k,i)+1

          enlO1(k,j)=enlO1(k,j)+1
          enlH1(k,i)=enlH1(k,i)+1

          IF (life_time.eq."y") THEN
           matrix1(k,i,j)=matrix1(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_ammo1(k,j)%N,1)
          END IF

          WRITE(17,*)  k, i, j, nHBpmolec1(k), dOHin, "6"
          IF (number_of_pdb.eq.2) THEN
                  WRITE(19,*)  k, i, j, nHBpmolec1(k), dOHin, "6"
          END IF
      END IF !dOH.lt.dOHmin
   END IF !dOO.LT.dOOmin							!condition 1
END IF
