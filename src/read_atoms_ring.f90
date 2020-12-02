 IMPLICIT NONE
   INTEGER::flag, lineas, i, j, k, num, maxi(nmovies)
   CHARACTER(LEN=30):: rubbish1
   CHARACTER(LEN=24):: rubbish2
   SELECT CASE (flag)
     CASE (1)
        filename="movie1.pdb"
        lineas=n1
        num=11
     CASE (2)
        filename="movie2.pdb"
        lineas=n2
        num=12
   END SELECT

   OPEN(num,FILE=filename,STATUS='OLD', ACTION='READ')
   DO k=1, nmovies
      READ(num,*)       !MODEL
      READ(num,*)       !CRYST1
      DO j=1, maxi(k)
        SELECT CASE (flag)
          CASE (1)
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C1(1), molecule_ring1(k,j)%C1(2),&
           molecule_ring1(k,j)%C1(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C2(1), molecule_ring1(k,j)%C2(2),&
           molecule_ring1(k,j)%C2(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C3(1), molecule_ring1(k,j)%C3(2),&
           molecule_ring1(k,j)%C3(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C4(1), molecule_ring1(k,j)%C4(2),&
           molecule_ring1(k,j)%C4(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C5(1), molecule_ring1(k,j)%C5(2),&
           molecule_ring1(k,j)%C5(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring1(k,j)%C6(1), molecule_ring1(k,j)%C6(2),&
           molecule_ring1(k,j)%C6(3), rubbish2
          CASE (2)
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C1(1), molecule_ring2(k,j)%C1(2),&
           molecule_ring2(k,j)%C1(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C2(1), molecule_ring2(k,j)%C2(2),&
           molecule_ring2(k,j)%C2(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C3(1), molecule_ring2(k,j)%C3(2),&
           molecule_ring2(k,j)%C3(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C4(1), molecule_ring2(k,j)%C4(2),&
           molecule_ring2(k,j)%C4(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C5(1), molecule_ring2(k,j)%C5(2),&
           molecule_ring2(k,j)%C5(3), rubbish2
           READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_ring2(k,j)%C6(1), molecule_ring2(k,j)%C6(2),&
           molecule_ring2(k,j)%C6(3), rubbish2
        END SELECT
      END DO
      READ(num,*)  !ENDMDL
   END DO
!   IF (nmovies.ne.k) STOP 'PROGRAM FAILED'
   CLOSE (UNIT=num)
