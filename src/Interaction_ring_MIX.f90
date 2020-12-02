!*****************************mezcla de aromáticos******************************************
IF((mole1.eq."ring").and.(mole2.eq."ring").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
  dCD=dFFrr
  dangleCD=0.5
  dFE=dFErr
  dangleFE=0.86602540378

  CALL distance(molecule_ring1(k,i)%centro,molecule_ring2(k,j)%centro,dist)

  IF(dist.lt.dCD)THEN
     CALL vectorial_product(molecule_ring1(k,i)%v_nor(1),molecule_ring1(k,i)%v_nor(2),molecule_ring1(k,i)%v_nor(3),&
          molecule_ring2(k,j)%v_nor(1),molecule_ring2(k,j)%v_nor(2),molecule_ring2(k,j)%v_nor(3),modu)
          modu=ABS(modu)
     IF(modu.lt.dangleCD)THEN  
             nHBpmolecm(k)=nHBpmolecm(k)+2

             nn3=nn3+1; d3=d3+1             
             
             enlm(k,i)=enlm(k,i)+1
             enlm(k,j)=enlm(k,j)+1
             
          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF             
             
             WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), "1"
     END IF
  END IF

  IF(dist.lt.dFE)THEN
          CALL vectorial_product(molecule_ring1(k,i)%v_nor(1),molecule_ring1(k,i)%v_nor(2),molecule_ring1(k,i)%v_nor(3),&
               molecule_ring2(k,j)%v_nor(1),molecule_ring2(k,j)%v_nor(2),molecule_ring2(k,j)%v_nor(3),modu)
          modu=ABS(modu)
          IF(modu.gt.dangleFE)THEN
             nHBpmolecm(k)=nHBpmolecm(k)+2

             nn3=nn3+1; d3=d3+1             
             
             enlm(k,i)=enlm(k,i)+1
             enlm(k,j)=enlm(k,j)+1

          IF (life_time.eq."y") THEN
           matrixm(k,i,j)=matrixm(k,i,j)+1
          END IF

             WRITE(19,*)  k, i, maxnum1+j, nHBpmolecm(k), "2"
          END IF
  END IF

END IF
