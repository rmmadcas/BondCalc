IMPLICIT NONE
INTEGER:: a, aux
INTEGER:: i,j,k, ii, jj

OPEN(31,file="timelife1.dat",STATUS='NEW', ACTION='WRITE')
WRITE(31,*) "#n_of_Bonds1, frecuency"
WRITE(31,*)"#========================================="

old=matrix1
new=matrix1
a=0
DO k=1, nmovies-1
    DO i=1, maxi1(k)
        DO j=i+1, maxi1(k)
            IF (old(k,i,j).eq.0) THEN
                new(k+1,i,j) = 0
            END IF
        END DO !j 
    END DO !i
    
    aux=0
    DO ii=1, maxi1(k)
        DO jj=ii+1, maxi1(k)
            aux=aux+new(k,ii,jj)
        END DO !j
    END DO !i
    WRITE(31,*)k-a, aux
    old=new
    IF (aux.eq.0) THEN
        old=matrix1
        new=matrix1
        WRITE(31,*)k-a+1, 0       !hacemos esto para ajustar un poco mejor la cola de la distribucion
        WRITE(31,*)k-a+2, 0
        WRITE(31,*)k-a+3, 0
        WRITE(31,*)k-a+4, 0
        WRITE(31,*)k-a+5, 0
        a=k
    END IF
END DO !k

CLOSE (UNIT=31)

IF (number_of_pdb.eq.2) THEN

 DEALLOCATE (old, new)

 OPEN(32,file="timelife2.dat",STATUS='NEW', ACTION='WRITE')
 WRITE(32,*) "#n_of_Bonds2, frecuency"
 WRITE(32,*)"#========================================="

 ALLOCATE (old(1:nmovies,1:maxnum2,1:maxnum2))
 ALLOCATE (new(1:nmovies,1:maxnum2,1:maxnum2))
  
 old=matrix2
 new=matrix2
 a=0
 DO k=1, nmovies-1
     DO i=1, maxi2(k)
         DO j=i+1, maxi2(k)
             IF (old(k,i,j).eq.0) THEN
                 new(k+1,i,j) = 0
             END IF
         END DO !j 
     END DO !i
     
     aux=0
     DO ii=1, maxi2(k)
         DO jj=ii+1, maxi2(k)
             aux=aux+new(k,ii,jj)
         END DO !j
     END DO !i
     WRITE(32,*)k-a, aux
     old=new
     IF (aux.eq.0) THEN
         old=matrix2
         new=matrix2
         WRITE(32,*)k-a+1, 0       !hacemos esto para ajustar un poco mejor la cola de la distribucion
         WRITE(32,*)k-a+2, 0
         WRITE(32,*)k-a+3, 0
         WRITE(32,*)k-a+4, 0
         WRITE(32,*)k-a+5, 0
         a=k
     END IF
 END DO !k

 CLOSE (UNIT=32)


 DEALLOCATE (old, new)

 OPEN(33,file="timelifeMix.dat",STATUS='NEW', ACTION='WRITE')
 WRITE(33,*) "#n_of_BondsMix, frecuency"
 WRITE(33,*)"#========================================="

 ALLOCATE (old(1:nmovies,1:maxnum1,1:maxnum2))
 ALLOCATE (new(1:nmovies,1:maxnum1,1:maxnum2))
  
 old=matrixm
 new=matrixm
 a=0
 DO k=1, nmovies-1
     DO i=1, maxi1(k)
         DO j=i+1, maxi2(k)
             IF (old(k,i,j).eq.0) THEN
                 new(k+1,i,j) = 0
             END IF
         END DO !j 
     END DO !i
     
     aux=0
     DO ii=1, maxi1(k)
         DO jj=ii+1, maxi2(k)
             aux=aux+new(k,ii,jj)
         END DO !j
     END DO !i
     WRITE(33,*)k-a, aux
     old=new
     IF (aux.eq.0) THEN
         old=matrixm
         new=matrixm
         WRITE(33,*)k-a+1, 0       !hacemos esto para ajustar un poco mejor la cola de la distribucion
         WRITE(33,*)k-a+2, 0
         WRITE(33,*)k-a+3, 0
         WRITE(33,*)k-a+4, 0
         WRITE(33,*)k-a+5, 0
         a=k
     END IF
 END DO !k

 CLOSE (UNIT=33)

END IF
