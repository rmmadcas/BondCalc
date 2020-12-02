!*********************************alcohol-agua*****************************
IF(((mole1.eq."alco").and.(mole2.eq."wate".and.(flag1.eq.1).and.(flag2.eq.2)))) THEN
dOOmin=dOOwal
dOHmin=dOHwal
theta_min=pi/6

CALL distance(molecule_alcohol1(k,i)%O,molecule_water2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_water2(k,j)%H1,dOH)         
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
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      
      CALL distance(molecule_alcohol1(k,i)%O,molecule_water2(k,j)%H2,dOH)         
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
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min) 	!condition 3.2
      END IF !(dOH.lt.dOHmin)		!condition 2.2
      
      CALL distance(molecule_water2(k,j)%O,molecule_alcohol1(k,i)%H,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.3
        CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,dOHin)         
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
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"

        END IF !(theta.lt.theta_min)    !condition 3.3
      END IF !(dOH.lt.dOHmin)           !condition 2.3

   END IF !dOO.LT.dOOmin							!condition 1
END IF

!*****************************mezcla de alcoholes*******************************
IF((mole1.eq."alco").and.(mole2.eq."alco").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
dOOmin=dOOalal
dOHmin=dOHalal
theta_min=pi/6

   CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol2(k,j)%H,dOH)         
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
           CALL SDF_alcohol(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol2(k,j)%O,molecule_alcohol1(k,i)%H,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,dOHin)         
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
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2
      
   END IF !dOO.LT.dOOmin							!condition 1

END IF

!*****************************mezcla de alcohol-diol************************************
IF((mole1.eq."alcohol").and.(mole2.eq."diol").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
dOOmin=dOOald
dOHmin=dOHald
theta_min=pi/6

   CALL distance(molecule_alcohol1(k,i)%O,molecule_diol2(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_diol2(k,j)%H1,dOH)         
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
           CALL SDF_alcohol(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol2(k,j)%O1,molecule_alcohol1(k,i)%H,dOH)
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.2
        CALL distance(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,dOHin)         
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
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_diol2(k,j)%O1,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

   END IF !dOO.LT.dOOmin							!condition 1.1
   
   CALL distance(molecule_alcohol1(k,i)%O,molecule_diol2(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1.2

      CALL distance(molecule_alcohol1(k,i)%O,molecule_diol2(k,j)%H2,dOH)         
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
           CALL SDF_alcohol(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"

        END IF !(theta.lt.theta_min) 	!condition 3.1
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_diol2(k,j)%O2,molecule_alcohol1(k,i)%H,dOH)
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
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"

        END IF !(theta.lt.theta_min)    !condition 3.2
      END IF !(dOH.lt.dOHmin)           !condition 2.2

      
   END IF !dOO.LT.dOOmin							!condition 1.2
   
END IF ! alcohol-diol   

!*********************** alcohol-ammoniaco*****************************
IF((mole1.eq."alcohol").and.(mole2.eq."ammo").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dONalam
   dOHmin=dOHalam

   CALL distance(molecule_alcohol1(k,i)%O,molecule_ammo2(k,j)%N,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_ammo2(k,j)%H1,dOH)         
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
                  molecule_ammo2(k,j)%H3,molecule_alcohol1(k,i)%O,4)
         END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_ammo2(k,j)%H2,dOH)         
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
                  molecule_ammo2(k,j)%H3,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_alcohol1(k,i)%O,molecule_ammo2(k,j)%H3,dOH)         
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
                   molecule_ammo2(k,j)%H3,molecule_alcohol1(k,i)%O,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_alcohol1(k,i)%H,molecule_ammo2(k,j)%N,dOH)
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
           CALL SDF_alcohol(molecule_alcohol1(k,i)%O,molecule_alcohol1(k,i)%H,molecule_ammo2(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1
END IF

