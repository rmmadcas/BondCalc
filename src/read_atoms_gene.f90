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
           DO l=1,NumMol(1)
                   READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_gene1(k,j)%A(l,1), molecule_gene1(k,j)%A(l,2),&
                   molecule_gene1(k,j)%A(l,3), rubbish2
           END DO
          CASE (2)
           DO l=1,NumMol(2)
                   READ(num,'(A30,3F8.3,A24)')rubbish1, molecule_gene2(k,j)%A(l,1), molecule_gene2(k,j)%A(l,2),&
                   molecule_gene2(k,j)%A(l,3), rubbish2
           END DO
        END SELECT
      END DO
      READ(num,*)  !ENDMDL
   END DO

!   IF (nmovies.ne.k) STOP 'PROGRAM FAILED'

   CLOSE (UNIT=num)
