!**********************************amoniaco-agua**********************************
IF((mole1.eq."ammo").and.(mole2.eq."water").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dONwam
   dOHmin=dOHwam

   CALL distance(molecule_ammo1(k,i)%N,molecule_water2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_water2(k,j)%H1,dOH)         
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
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,i)%N,molecule_water2(k,j)%H2,dOH)         
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
           CALL SDF_water(molecule_water2(k,j)%O,molecule_water2(k,j)%H1,molecule_water2(k,j)%H2,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_ammo1(k,i)%H1,molecule_water2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H2,molecule_water2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H3,molecule_water2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,j)%H3,molecule_water2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1
END IF


!**********************************amoniaco-alcohol**********************************
IF((mole1.eq."ammo").and.(mole2.eq."alcohol").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dONalam
   dOHmin=dOHalam

   CALL distance(molecule_ammo1(k,i)%N,molecule_alcohol2(k,j)%O,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_alcohol2(k,j)%H,dOH)         
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
           CALL SDF_alcohol(molecule_alcohol2(k,j)%O,molecule_alcohol2(k,j)%H,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_ammo1(k,i)%H1,molecule_alcohol2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H2,molecule_alcohol2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H3,molecule_alcohol2(k,j)%O,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_alcohol2(k,j)%O,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1
END IF


!*********************** ammoniaco-diol*****************************
IF((mole1.eq."ammo").and.(mole2.eq."diol").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dONdam
   dOHmin=dOHdam

   CALL distance(molecule_ammo1(k,i)%N,molecule_diol2(k,j)%O1,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_diol2(k,j)%H1,dOH)         
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
          CALL SDF_alcohol(molecule_diol2(k,j)%O1,molecule_diol2(k,j)%H1,molecule_ammo1(k,i)%N,4)
         END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_ammo1(k,i)%H1,molecule_diol2(k,j)%O1,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O1,3)
          END IF
 
          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H2,molecule_diol2(k,j)%O1,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O1,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O1,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O1,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1

   CALL distance(molecule_ammo1(k,i)%N,molecule_diol2(k,j)%O2,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_diol2(k,j)%H2,dOH)         
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
          CALL SDF_alcohol(molecule_diol2(k,j)%O2,molecule_diol2(k,j)%H2,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_ammo1(k,i)%H1,molecule_diol2(k,j)%O2,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H2,molecule_diol2(k,j)%O2,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

      CALL distance(molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O2,dOH)
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
           CALL SDF_ammo(molecule_ammo1(k,i)%N,molecule_ammo1(k,i)%H1,molecule_ammo1(k,i)%H2,&
                   molecule_ammo1(k,i)%H3,molecule_diol2(k,j)%O2,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)           !condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1
END IF

!*****************************mixture de ammo*********************************
IF((mole1.eq."ammo").and.(mole2.eq."ammo").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
   dOOmin=dNN
   dOHmin=dNH

   CALL distance(molecule_ammo1(k,i)%N,molecule_ammo2(k,j)%N,dOO)
   IF (dOO.lt.dOOmin) THEN                                                      !condition 1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo2(k,j)%H1,dOH)         
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
                   molecule_ammo2(k,j)%H3,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "1"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo2(k,j)%H2,dOH)         
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
                   molecule_ammo2(k,j)%H3,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "2"
    
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,i)%N,molecule_ammo2(k,j)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,i)=enlm(k,i)+1
          enlm(k,j)=enlm(k,j)+1

          enlOm(k,i)=enlOm(k,i)+1
          enlHm(k,j)=enlHm(k,j)+1

          IF (life_time.eq."y") THEN
           matrix2(k,i,j)=matrix2(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,j)%N,molecule_ammo2(k,j)%H1,molecule_ammo2(k,j)%H2,&
                   molecule_ammo2(k,j)%H3,molecule_ammo1(k,i)%N,4)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "3"
      END IF !(dOH.lt.dOHmin)		!condition 2.1


      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo2(k,i)%H1,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,j)=enlm(k,j)+1
          enlm(k,i)=enlm(k,i)+1

          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,i)%N,molecule_ammo2(k,i)%H1,molecule_ammo2(k,i)%H2,&
                   molecule_ammo2(k,i)%H3,molecule_ammo1(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "4"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo2(k,i)%H2,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,j)=enlm(k,j)+1
          enlm(k,i)=enlm(k,i)+1

          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,i)%N,molecule_ammo2(k,i)%H1,molecule_ammo2(k,i)%H2,&
                   molecule_ammo2(k,i)%H3,molecule_ammo1(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "5"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

      CALL distance(molecule_ammo1(k,j)%N,molecule_ammo2(k,i)%H3,dOH)         
      IF (dOH.lt.dOHmin) THEN                                                   !condition 2.1
          nHBpmolecm(k)=nHBpmolecm(k)+2
          
          nn3=nn3+1; d3=d3+1
          
          enlm(k,j)=enlm(k,j)+1
          enlm(k,i)=enlm(k,i)+1

          enlOm(k,j)=enlOm(k,j)+1
          enlHm(k,i)=enlHm(k,i)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

          IF (SDF.eq."y") THEN
           CALL SDF_ammo(molecule_ammo2(k,i)%N,molecule_ammo2(k,i)%H1,molecule_ammo2(k,i)%H2,&
                   molecule_ammo2(k,i)%H3,molecule_ammo1(k,j)%N,3)
          END IF

          WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), dOHin, "6"
      END IF !(dOH.lt.dOHmin)		!condition 2.1

   END IF !dOO.LT.dOOmin							!condition 1

END IF

