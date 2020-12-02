!*****************************mezcla de diol-agua******************************************
IF(((mole1.eq."diol").and.(mole2.eq."wate".and.(flag1.eq.1).and.(flag2.eq.2)))) THEN
dOOmin=dOOwd
dOHmin=dOHwd
theta_min=pi/6

CALL distance(molecule_diol1(k,i)%O1,molecule_water2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_diol1(k,i)%O1,molecule_water2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H2,molecule_water2(k,j)%H2,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      
      CALL distance(molecule_diol1(k,i)%O1,molecule_water2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_water2(k,j)%O,molecule_water2(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H2,molecule_water2(k,j)%H2,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min) 	!condition 3.2
      END IF !(dOH.lt.dOHmin)		!condition 2.2
      
      CALL distance(molecule_water2(k,j)%O,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.3
        CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.3
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"

        END IF !(theta.lt.theta_min)    !condition 3.3
      END IF !(dOH.lt.dOHmin)           !condition 2.3

   END IF !dOO.LT.dOOmin							!condition 1


CALL distance(molecule_diol1(k,i)%O2,molecule_water2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_diol1(k,i)%O2,molecule_water2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      
      CALL distance(molecule_diol1(k,i)%O2,molecule_water2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_water2(k,j)%O,molecule_water2(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "5"

        END IF !(theta.lt.theta_min) 	!condition 3.2
      END IF !(dOH.lt.dOHmin)		!condition 2.2
      
      CALL distance(molecule_water2(k,j)%O,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.3
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.3
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "6"

        END IF !(theta.lt.theta_min)    !condition 3.3
      END IF !(dOH.lt.dOHmin)           !condition 2.3

   END IF !dOO.LT.dOOmin							!condition 1
END IF


!*****************************mezcla de diol-alcohol*************************************
IF((mole1.eq."diol").and.(mole2.eq."alcohol").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
dOOmin=dOOald
dOHmin=dOHald
theta_min=pi/6

   CALL distance(molecule_diol1(k,i)%O1,molecule_alcohol2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.1

      CALL distance(molecule_diol1(k,i)%O1,molecule_alcohol2(k,j)%H,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol2(k,j)%O,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2
      
   END IF !dOO.LT.dOOmin							!condition 1.1
   
   CALL distance(molecule_diol1(k,i)%O2,molecule_alcohol2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.3

      CALL distance(molecule_diol1(k,i)%O2,molecule_alcohol2(k,j)%H,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol2(k,j)%O,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      
   END IF !dOO.LT.dOOmin							!condition 1.3ule_diol1(k,i)%H2,dOHin)         
END IF

!*****************************mezcla de dioles*******************************
IF((mole1.eq."diol").and.(mole2.eq."diol").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
dOOmin=dOOdd
dOHmin=dOHdd
theta_min=pi/6

   CALL distance(molecule_diol1(k,i)%O1,molecule_diol2(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.1

      CALL distance(molecule_diol1(k,i)%O1,molecule_diol2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol2(k,j)%O1,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_diol2(k,j)%O1,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      
   END IF !dOO.LT.dOOmin							!condition 1.1

   
   CALL distance(molecule_diol1(k,i)%O1,molecule_diol2(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.2

      CALL distance(molecule_diol1(k,i)%O1,molecule_diol2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol2(k,j)%O2,molecule_diol1(k,i)%H1,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol1(k,i)%H1,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      
   END IF !dOO.LT.dOOmin							!condition 1.2
   
      CALL distance(molecule_diol1(k,i)%O2,molecule_diol2(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.3

      CALL distance(molecule_diol1(k,i)%O2,molecule_diol2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "5"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol2(k,j)%O1,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,molecule_diol2(k,j)%O1,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "6"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      
   END IF !dOO.LT.dOOmin							!condition 1.3ule_diol1(k,i)%H2,dOHin)         
   
      CALL distance(molecule_diol1(k,i)%O2,molecule_diol2(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.4

      CALL distance(molecule_diol1(k,i)%O2,molecule_diol2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
        CALL distance(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "7"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1ule_diol1(k,i)%H2,dOHin)         
        CALL bend

      CALL distance(molecule_diol2(k,j)%O2,molecule_diol1(k,i)%H2,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_diol1(k,i)%O2,molecule_diol1(k,i)%H2,dOHin)         
        CALL bend
        IF (theta.lt.theta_min) THEN                                            !condition 3.2
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1          
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1
          
          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol2(k,i)%O2,molecule_diol2(k,i)%H2,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "8"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2
   END IF !dOO.LT.dOOmin							!condition 1.4
   
END IF

!*********************** diol-ammoniaco*****************************
IF((mole1.eq."diol").and.(mole2.eq."ammo").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dONdam
   dOHmin=dOHdam

   CALL distance(molecule_diol1(k,i)%O1,molecule_ammo2(k,j)%N,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_diol1(k,i)%O1,molecule_ammo2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%O1,molecule_ammo2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

         IF (SDF.eq."y") THEN
          CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                  molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O1,4)
         END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%O1,molecule_ammo2(k,j)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF
  
          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O1,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%H1,molecule_ammo2(k,j)%N,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1

          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O1,molecule_diol2(k,i)%H1,molecule_ammo1(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1


   CALL distance(molecule_diol1(k,i)%O2,molecule_ammo2(k,j)%N,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_diol1(k,i)%O2,molecule_ammo2(k,j)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%O2,molecule_ammo2(k,j)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%O2,molecule_ammo2(k,j)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF
  
          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_diol1(k,i)%O2,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol1(k,i)%H2,molecule_ammo2(k,j)%N,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2

          nn3=nn3+1; d3=d3+1

          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_alcohol(molecule_diol1(k,i)%O2,molecule_diol2(k,i)%H2,molecule_ammo2(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1
END IF
