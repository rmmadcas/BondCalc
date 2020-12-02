IMPLICIT NONE
CHARACTER(LEN=4):: mole1, mole2
INTEGER:: flag1, flag2, ii, jj, kk, a, b, aa, bb
REAL(KIND=r):: va(3), vb(3), vaa(3), vbb(3)
REAL(KIND=r):: dab, daabb

!*********************************1 MOVIE************************************
!*********************************pures 1************************************
IF((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.1).and.(flag2.eq.1)) THEN

   DO kk=1,NumInteraction1-1, 2
        a=list_int(kk,1)
        b=list_int(kk,2)
        aa=list_int(kk+1,1)
        bb=list_int(kk+1,2)

        va(1)=molecule_gene1(k,i)%A(a,1)
        va(2)=molecule_gene1(k,i)%A(a,2)
        va(3)=molecule_gene1(k,i)%A(a,3)

        vb(1)=molecule_gene1(k,j)%A(b,1)
        vb(2)=molecule_gene1(k,j)%A(b,2)
        vb(3)=molecule_gene1(k,j)%A(b,3)

        vaa(1)=molecule_gene1(k,i)%A(aa,1)
        vaa(2)=molecule_gene1(k,i)%A(aa,2)
        vaa(3)=molecule_gene1(k,i)%A(aa,3)

        vbb(1)=molecule_gene1(k,j)%A(bb,1)
        vbb(2)=molecule_gene1(k,j)%A(bb,2)
        vbb(3)=molecule_gene1(k,j)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_Interaction1(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_Interaction1(aa,bb)) THEN    !condition 2
            nHBpmolec1(k)=nHBpmolec1(k)+2
 
            nn1=nn1+1; d1=d1+1
          
            enl1(k,i)=enl1(k,i)+1
            enl1(k,j)=enl1(k,j)+1

            IF (life_time.eq."y") THEN
             matrix1(k,i,j)=matrix1(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene1(k,j)%A,molecule_gene1(k,i)%A(a,1),molecule_gene1(k,i)%A(a,2),molecule_gene1(k,i)%A(a,3),1,dab)
            END IF

            WRITE(17,*)  k, i, j, nHBpmolec1(k), "1"
          END IF !aabb   condtion 2
        END IF !dab      conditio 1

        va(1)=molecule_gene1(k,j)%A(a,1)
        va(2)=molecule_gene1(k,j)%A(a,2)
        va(3)=molecule_gene1(k,j)%A(a,3)

        vb(1)=molecule_gene1(k,i)%A(b,1)
        vb(2)=molecule_gene1(k,i)%A(b,2)
        vb(3)=molecule_gene1(k,i)%A(b,3)

        vaa(1)=molecule_gene1(k,j)%A(aa,1)
        vaa(2)=molecule_gene1(k,j)%A(aa,2)
        vaa(3)=molecule_gene1(k,j)%A(aa,3)

        vbb(1)=molecule_gene1(k,i)%A(bb,1)
        vbb(2)=molecule_gene1(k,i)%A(bb,2)
        vbb(3)=molecule_gene1(k,i)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_Interaction1(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_Interaction1(aa,bb)) THEN    !condition 2
            nHBpmolec1(k)=nHBpmolec1(k)+2
          
            nn1=nn1+1; d1=d1+1
          
            enl1(k,i)=enl1(k,i)+1
            enl1(k,j)=enl1(k,j)+1

            IF (life_time.eq."y") THEN
             matrix1(k,i,j)=matrix1(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene1(k,i)%A,molecule_gene1(k,j)%A(b,1),molecule_gene1(k,j)%A(b,2),molecule_gene1(k,j)%A(b,3),1,dab)
            END IF

            WRITE(17,*)  k, i, j, nHBpmolec1(k), "2"
          END IF !aabb   condtion 2
        END IF !dab      conditio 1

   END DO                               !kk
END IF   !IF ((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.1).and.(flag2.eq.1)) THEN

!*********************************pures 2************************************
IF((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.2).and.(flag2.eq.2)) THEN

   DO kk=NumInteraction1+1,NumInteraction1+NumInteraction2-1,2
        a=list_int(kk,1)
        b=list_int(kk,2)
        aa=list_int(kk+1,1)
        bb=list_int(kk+1,2)
        
        va(1)=molecule_gene2(k,i)%A(a,1)
        va(2)=molecule_gene2(k,i)%A(a,2)
        va(3)=molecule_gene2(k,i)%A(a,3)

        vb(1)=molecule_gene2(k,j)%A(b,1)
        vb(2)=molecule_gene2(k,j)%A(b,2)
        vb(3)=molecule_gene2(k,j)%A(b,3)

        vaa(1)=molecule_gene2(k,i)%A(aa,1)
        vaa(2)=molecule_gene2(k,i)%A(aa,2)
        vaa(3)=molecule_gene2(k,i)%A(aa,3)

        vbb(1)=molecule_gene2(k,j)%A(bb,1)
        vbb(2)=molecule_gene2(k,j)%A(bb,2)
        vbb(3)=molecule_gene2(k,j)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_Interaction2(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_Interaction2(aa,bb)) THEN    !condition 2
            nHBpmolec2(k)=nHBpmolec2(k)+2
          
            nn2=nn2+1; d2=d2+1
          
            enl2(k,i)=enl2(k,i)+1
            enl2(k,j)=enl2(k,j)+1

            IF (life_time.eq."y") THEN
             matrix2(k,i,j)=matrix2(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene2(k,j)%A,molecule_gene2(k,i)%A(a,1),molecule_gene2(k,i)%A(a,2),molecule_gene2(k,i)%A(a,3),2,dab)
            END IF

            WRITE(18,*)  k, i, j, nHBpmolec2(k), "1"
          END IF !aabb   condtiion 2
        END IF !dab      condition 1

        va(1)=molecule_gene2(k,j)%A(a,1)
        va(2)=molecule_gene2(k,j)%A(a,2)
        va(3)=molecule_gene2(k,j)%A(a,3)

        vb(1)=molecule_gene2(k,i)%A(b,1)
        vb(2)=molecule_gene2(k,i)%A(b,2)
        vb(3)=molecule_gene2(k,i)%A(b,3)

        vaa(1)=molecule_gene2(k,j)%A(aa,1)
        vaa(2)=molecule_gene2(k,j)%A(aa,2)
        vaa(3)=molecule_gene2(k,j)%A(aa,3)

        vbb(1)=molecule_gene2(k,i)%A(bb,1)
        vbb(2)=molecule_gene2(k,i)%A(bb,2)
        vbb(3)=molecule_gene2(k,i)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_Interaction2(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_Interaction2(aa,bb)) THEN    !condition 2
            nHBpmolec2(k)=nHBpmolec2(k)+2
          
            nn2=nn2+1; d2=d2+1
          
            enl2(k,i)=enl2(k,i)+1
            enl2(k,j)=enl2(k,j)+1

            IF (life_time.eq."y") THEN
             matrix2(k,i,j)=matrix2(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene2(k,i)%A,molecule_gene2(k,j)%A(b,1),molecule_gene2(k,j)%A(b,2),molecule_gene2(k,j)%A(b,3),2,dab)
            END IF

            WRITE(18,*)  k, i, j, nHBpmolec2(k), "2"
          END IF !aabb   condtiion 2
        END IF !dab      condition 1

   END DO                               !kk
END IF   !IF ((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.2).and.(flag2.eq.2)) THEN
!*********************************MIXTURE**********************************
!*********************************gene-gene*****************************
IF((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.1).and.(flag2.eq.2)) THEN

   DO kk=NumInteraction1+NumInteraction2+1,NumInteraction1+NumInteraction2+NumInteractionMix-1,2
        a=list_int(kk,1)
        b=list_int(kk,2)
        aa=list_int(kk+1,1)
        bb=list_int(kk+1,2)

        va(1)=molecule_gene1(k,i)%A(a,1)
        va(2)=molecule_gene1(k,i)%A(a,2)
        va(3)=molecule_gene1(k,i)%A(a,3)

        vb(1)=molecule_gene2(k,j)%A(b,1)
        vb(2)=molecule_gene2(k,j)%A(b,2)
        vb(3)=molecule_gene2(k,j)%A(b,3)

        vaa(1)=molecule_gene1(k,i)%A(aa,1)
        vaa(2)=molecule_gene1(k,i)%A(aa,2)
        vaa(3)=molecule_gene1(k,i)%A(aa,3)

        vbb(1)=molecule_gene2(k,j)%A(bb,1)
        vbb(2)=molecule_gene2(k,j)%A(bb,2)
        vbb(3)=molecule_gene2(k,j)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_InteractionMix(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_InteractionMix(aa,bb)) THEN    !condition 2
            nHBpmolecm(k)=nHBpmolecm(k)+2
          
            nn3=nn3+1; d3=d3+1
          
            enlm(k,i)=enlm(k,i)+1
            enlm(k,j)=enlm(k,maxnum1+j)+1

            IF (life_time.eq."y") THEN
             matrixm(k,i,j)=matrixm(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene2(k,j)%A,molecule_gene1(k,i)%A(b,1),molecule_gene1(k,i)%A(b,2),molecule_gene1(k,i)%A(b,3),4,dab)
            END IF

           WRITE(19,*)  k, i, j, nHBpmolecm(k), "1"
          END IF !aabb   condtiion 2
        END IF !dab      condition 1

        va(1)=molecule_gene1(k,j)%A(a,1)
        va(2)=molecule_gene1(k,j)%A(a,2)
        va(3)=molecule_gene1(k,j)%A(a,3)

        vb(1)=molecule_gene2(k,i)%A(b,1)
        vb(2)=molecule_gene2(k,i)%A(b,2)
        vb(3)=molecule_gene2(k,i)%A(b,3)

        vaa(1)=molecule_gene1(k,j)%A(aa,1)
        vaa(2)=molecule_gene1(k,j)%A(aa,2)
        vaa(3)=molecule_gene1(k,j)%A(aa,3)

        vbb(1)=molecule_gene2(k,i)%A(bb,1)
        vbb(2)=molecule_gene2(k,i)%A(bb,2)
        vbb(3)=molecule_gene2(k,i)%A(bb,3)

        CALL distance(va,vb,dab)
        IF (dab.lt.Generic_InteractionMix(a,b)) THEN          !condition 1 

          CALL distance(vaa,vbb,daabb)         
          IF (daabb.lt.Generic_InteractionMix(aa,bb)) THEN    !condition 2
            nHBpmolecm(k)=nHBpmolecm(k)+2
          
            nn3=nn3+1; d3=d3+1
          
            enlm(k,i)=enlm(k,i)+1
            enlm(k,j)=enlm(k,j)+1

            IF (life_time.eq."y") THEN
             matrixm(k,i,j)=matrixm(k,i,j)+1
            END IF

            IF (SDF.eq."y") THEN
CALL SDF_gene(molecule_gene1(k,i)%A,molecule_gene1(k,j)%A(b,1),molecule_gene1(k,j)%A(b,2),molecule_gene1(k,j)%A(b,3),3,dab)
            END IF

            WRITE(19,*)  k, i, j, nHBpmolecm(k), "2"
          END IF !aabb   condtiion 2
        END IF !dab      condition 1

   END DO                               !kk

END IF! IF((mole1.eq."gene").and.(mole2.eq."gene").and.(flag1.eq.1).and.(flag2.eq.2)) THEN
