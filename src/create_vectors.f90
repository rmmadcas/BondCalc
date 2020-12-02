IMPLICIT NONE
INTEGER::flag, lineas, i, j, k, num, maxi(nmovies)
CHARACTER(LEN=30):: rubbish1
CHARACTER(LEN=24):: rubbish2
SELECT CASE (flag)
       CASE (1)
        DO k=1, nmovies
          DO j=1, maxi(k)
             molecule_ring1(k,j)%Centro(1)=(molecule_ring1(k,j)%C1(1)+molecule_ring1(k,j)%C2(1)+molecule_ring1(k,j)%C3(1)&
                     +molecule_ring1(k,j)%C4(1)+molecule_ring1(k,j)%C5(1)+molecule_ring1(k,j)%C6(1))/6
             molecule_ring1(k,j)%Centro(2)=(molecule_ring1(k,j)%C1(2)+molecule_ring1(k,j)%C2(2)+molecule_ring1(k,j)%C3(2)&
                     +molecule_ring1(k,j)%C4(2)+molecule_ring1(k,j)%C5(2)+molecule_ring1(k,j)%C6(2))/6
             molecule_ring1(k,j)%Centro(3)=(molecule_ring1(k,j)%C1(3)+molecule_ring1(k,j)%C2(3)+molecule_ring1(k,j)%C3(3)&
                     +molecule_ring1(k,j)%C4(3)+molecule_ring1(k,j)%C5(3)+molecule_ring1(k,j)%C6(3))/6
          

            !Buscamos el vector normal, para ello buscamos dos vectores (1,3) y (1,5) y los multiplicamos vectorialmente
             vec_aux_a(1)=molecule_ring1(k,j)%C3(1)-molecule_ring1(k,j)%C1(1)
             vec_aux_a(2)=molecule_ring1(k,j)%C3(2)-molecule_ring1(k,j)%C1(2)
             vec_aux_a(3)=molecule_ring1(k,j)%C3(3)-molecule_ring1(k,j)%C1(3)

             vec_aux_b(1)=molecule_ring1(k,j)%C1(1)-molecule_ring1(k,j)%C5(1)
             vec_aux_b(2)=molecule_ring1(k,j)%C1(2)-molecule_ring1(k,j)%C5(2)
             vec_aux_b(3)=molecule_ring1(k,j)%C1(3)-molecule_ring1(k,j)%C5(3)

             molecule_ring1(k,j)%v_nor(1)=vec_aux_b(3)*vec_aux_a(2)-vec_aux_b(2)*vec_aux_a(3)
             molecule_ring1(k,j)%v_nor(2)=vec_aux_b(1)*vec_aux_a(3)-vec_aux_b(3)*vec_aux_a(1)
             molecule_ring1(k,j)%v_nor(3)=vec_aux_b(2)*vec_aux_a(1)-vec_aux_b(1)*vec_aux_a(2)

             norm=(molecule_ring1(k,j)%v_nor(1)**2)+(molecule_ring1(k,j)%v_nor(2)**2)+(molecule_ring1(k,j)%v_nor(3)**2)
             norm=sqrt(norm)

             molecule_ring1(k,j)%v_nor(1)=molecule_ring1(k,j)%v_nor(1)/norm
             molecule_ring1(k,j)%v_nor(2)=molecule_ring1(k,j)%v_nor(2)/norm
             molecule_ring1(k,j)%v_nor(3)=molecule_ring1(k,j)%v_nor(3)/norm
          END DO
        END DO

       CASE (2)
        DO k=1, nmovies
          DO j=1, maxi(k)
             molecule_ring2(k,j)%Centro(1)=(molecule_ring2(k,j)%C1(1)+molecule_ring2(k,j)%C2(1)+molecule_ring2(k,j)%C3(1)&
                     +molecule_ring2(k,j)%C4(1)+molecule_ring2(k,j)%C5(1)+molecule_ring2(k,j)%C6(1))/6
             molecule_ring2(k,j)%Centro(2)=(molecule_ring2(k,j)%C1(2)+molecule_ring2(k,j)%C2(2)+molecule_ring2(k,j)%C3(2)&
                     +molecule_ring2(k,j)%C4(2)+molecule_ring2(k,j)%C5(2)+molecule_ring2(k,j)%C6(2))/6
             molecule_ring2(k,j)%Centro(3)=(molecule_ring2(k,j)%C1(3)+molecule_ring2(k,j)%C2(3)+molecule_ring2(k,j)%C3(3)&
                     +molecule_ring2(k,j)%C4(3)+molecule_ring2(k,j)%C5(3)+molecule_ring2(k,j)%C6(3))/6


            !Buscamos el vector normal, para ello buscamos dos vectores (1,3) y (1,5) y los multiplicamos vectorialmente
             vec_aux_a(1)=molecule_ring2(k,j)%C3(1)-molecule_ring2(k,j)%C1(1)
             vec_aux_a(2)=molecule_ring2(k,j)%C3(2)-molecule_ring2(k,j)%C1(2)
             vec_aux_a(3)=molecule_ring2(k,j)%C3(3)-molecule_ring2(k,j)%C1(3)

             vec_aux_b(1)=molecule_ring2(k,j)%C1(1)-molecule_ring2(k,j)%C5(1)
             vec_aux_b(2)=molecule_ring2(k,j)%C1(2)-molecule_ring2(k,j)%C5(2)
             vec_aux_b(3)=molecule_ring2(k,j)%C1(3)-molecule_ring2(k,j)%C5(3)

             molecule_ring2(k,j)%v_nor(1)=vec_aux_b(3)*vec_aux_a(2)-vec_aux_b(2)*vec_aux_a(3)
             molecule_ring2(k,j)%v_nor(2)=vec_aux_b(1)*vec_aux_a(3)-vec_aux_b(3)*vec_aux_a(1)
             molecule_ring2(k,j)%v_nor(3)=vec_aux_b(2)*vec_aux_a(1)-vec_aux_b(1)*vec_aux_a(2)

             norm=(molecule_ring2(k,j)%v_nor(1)**2)+(molecule_ring2(k,j)%v_nor(2)**2)+(molecule_ring2(k,j)%v_nor(3)**2)
             norm=sqrt(norm)

             molecule_ring2(k,j)%v_nor(1)=molecule_ring2(k,j)%v_nor(1)/norm
             molecule_ring2(k,j)%v_nor(2)=molecule_ring2(k,j)%v_nor(2)/norm
             molecule_ring2(k,j)%v_nor(3)=molecule_ring2(k,j)%v_nor(3)/norm
          END DO
        END DO
END SELECT
