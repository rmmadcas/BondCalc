IMPLICIT NONE
INTEGER::n(0:6), nO(0:6), nH(0:6)
INTEGER::i,j,k, flag, maxi(nmovies)
REAL (KIND=r):: f0(nmovies),f1(nmovies),f2(nmovies),f3(nmovies),f4(nmovies),f5(nmovies),f6(nmovies)
REAL (KIND=r):: f0O(nmovies),f1O(nmovies),f2O(nmovies),f3O(nmovies),f4O(nmovies),f5O(nmovies),f6O(nmovies)
REAL (KIND=r):: f0H(nmovies),f1H(nmovies),f2H(nmovies),f3H(nmovies),f4H(nmovies),f5H(nmovies),f6H(nmovies)


IF (flag.eq.1) THEN
 maxi=maxi1 
   OPEN(23,file="data_bonds1.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(23,'(A42)') "#Movie molecule  TotalBonds O-Bonds H-Bonds Total-O-H "
   WRITE(23,'(A42)') "#========================================="

   OPEN(26,file="Ave_bonds1.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(26,'(A42)') "#f0, f1, f2, f3, f4, f5, f6, SUM (fi)"
   WRITE(26,'(A42)') "#========================================="

ELSE IF (flag.eq.2) THEN
 maxi=maxi2
   OPEN(24,file="data_bonds2.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(24,'(A42)') "#Movie molecule  TotalBonds O-Bonds H-Bonds Total-O-H "
   WRITE(24,'(A42)') "#========================================="

   OPEN(27,file="Ave_bonds2.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(27,'(A42)') "#f0, f1, f2, f3, f4, f5, f6, SUM (fi)"
   WRITE(27,'(A42)') "#========================================="

ELSE IF (flag.eq.3) THEN
 maxi=maxi1+maxi2
   OPEN(25,file="data_bondsMix.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(25,'(A42)') "#Movie molecule  TotalBonds O-Bonds H-Bonds Total-O-H "
   WRITE(25,'(A42)') "#========================================="

   OPEN(28,file="Ave_bondsMix.dat",STATUS='NEW', ACTION='WRITE')
   WRITE(28,'(A42)') "#f0, f1, f2, f3, f4, f5, f6, SUM (fi)"
   WRITE(28,'(A42)') "#========================================="

END IF

DO k=1, nmovies
n=0 ; nO=0; nH=0
 DO i=1, maxi(k)
  IF (flag.eq.1) THEN 
       WRITE(23,*) k, i, enl1(k,i), enlO1(k,i),enlH1(k,i),enl1(k,i)-enlO1(k,i)-enlH1(k,i)
       SELECT CASE (enl1(k,i))
       CASE (0)
               n(0)=n(0)+1
       CASE (1)
               n(1)=n(1)+1
       CASE (2)
               n(2)=n(2)+1
       CASE (3)
               n(3)=n(3)+1
       CASE (4)
               n(4)=n(4)+1
       CASE (5)
               n(5)=n(5)+1
       CASE (6)
               n(6)=n(6)+1
       END SELECT
      
       SELECT CASE (enlO1(k,i))
       CASE (0)
               nO(0)=nO(0)+1
       CASE (1)
               nO(1)=nO(1)+1
       CASE (2)
               nO(2)=nO(2)+1
       CASE (3)
               nO(3)=nO(3)+1
       CASE (4)
               nO(4)=nO(4)+1
       CASE (5)
               nO(5)=nO(5)+1
       CASE (6)
               nO(6)=nO(6)+1
       END SELECT

       SELECT CASE (enlH1(k,i))
       CASE (0)
               nH(0)=nH(0)+1
       CASE (1)
               nH(1)=nH(1)+1
       CASE (2)
               nH(2)=nH(2)+1
       CASE (3)
               nH(3)=nH(3)+1
       CASE (4)
               nH(4)=nH(4)+1
       CASE (5)
               nH(5)=nH(5)+1
       CASE (6)
               nH(6)=nH(6)+1
       END SELECT


  ELSE IF (flag.eq.2) THEN 
       WRITE(24,*) k, i, enl2(k,i), enlO2(k,i),enlH2(k,i),enl2(k,i)-enlO2(k,i)-enlH2(k,i)
       SELECT CASE (enl2(k,i))
       CASE (0)
               n(0)=n(0)+1
       CASE (1)
               n(1)=n(1)+1
       CASE (2)
               n(2)=n(2)+1
       CASE (3)
               n(3)=n(3)+1
       CASE (4)
               n(4)=n(4)+1
       CASE (5)
               n(5)=n(5)+1
       CASE (6)
               n(6)=n(6)+1
       END SELECT
      
       SELECT CASE (enlO2(k,i))
       CASE (0)
               nO(0)=nO(0)+1
       CASE (1)
               nO(1)=nO(1)+1
       CASE (2)
               nO(2)=nO(2)+1
       CASE (3)
               nO(3)=nO(3)+1
       CASE (4)
               nO(4)=nO(4)+1
       CASE (5)
               nO(5)=nO(5)+1
       CASE (6)
               nO(6)=nO(6)+1
       END SELECT

       SELECT CASE (enlH2(k,i))
       CASE (0)
               nH(0)=nH(0)+1
       CASE (1)
               nH(1)=nH(1)+1
       CASE (2)
               nH(2)=nH(2)+1
       CASE (3)
               nH(3)=nH(3)+1
       CASE (4)
               nH(4)=nH(4)+1
       CASE (5)
               nH(5)=nH(5)+1
       CASE (6)
               nH(6)=nH(6)+1
       END SELECT

  ELSE IF (flag.eq.3) THEN 
       WRITE(25,*) k, i, enlm(k,i), enlOm(k,i),enlHm(k,i),enlm(k,i)-enlOm(k,i)-enlHm(k,i)
       SELECT CASE (enlm(k,i))
       CASE (0)
               n(0)=n(0)+1
       CASE (1)
               n(1)=n(1)+1
       CASE (2)
               n(2)=n(2)+1
       CASE (3)
               n(3)=n(3)+1
       CASE (4)
               n(4)=n(4)+1
       CASE (5)
               n(5)=n(5)+1
       CASE (6)
               n(6)=n(6)+1
       END SELECT
      
       SELECT CASE (enlOm(k,i))
       CASE (0)
               nO(0)=nO(0)+1
       CASE (1)
               nO(1)=nO(1)+1
       CASE (2)
               nO(2)=nO(2)+1
       CASE (3)
               nO(3)=nO(3)+1
       CASE (4)
               nO(4)=nO(4)+1
       CASE (5)
               nO(5)=nO(5)+1
       CASE (6)
               nO(6)=nO(6)+1
       END SELECT

       SELECT CASE (enlHm(k,i))
       CASE (0)
               nH(0)=nH(0)+1
       CASE (1)
               nH(1)=nH(1)+1
       CASE (2)
               nH(2)=nH(2)+1
       CASE (3)
               nH(3)=nH(3)+1
       CASE (4)
               nH(4)=nH(4)+1
       CASE (5)
               nH(5)=nH(5)+1
       CASE (6)
               nH(6)=nH(6)+1
       END SELECT
  END IF
 END DO !i

 f0(k)=1.*n(0)/maxi(k)
 f1(k)=1.*n(1)/maxi(k)
 f2(k)=1.*n(2)/maxi(k)
 f3(k)=1.*n(3)/maxi(k)
 f4(k)=1.*n(4)/maxi(k)
 f5(k)=1.*n(5)/maxi(k)
 f6(k)=1.*n(6)/maxi(k)

 f0O(k)=1.*nO(0)/maxi(k)
 f1O(k)=1.*nO(1)/maxi(k)
 f2O(k)=1.*nO(2)/maxi(k)
 f3O(k)=1.*nO(3)/maxi(k)
 f4O(k)=1.*nO(4)/maxi(k)
 f5O(k)=1.*nO(5)/maxi(k)
 f6O(k)=1.*nO(6)/maxi(k)

 f0H(k)=1.*nH(0)/maxi(k)
 f1H(k)=1.*nH(1)/maxi(k)
 f2H(k)=1.*nH(2)/maxi(k)
 f3H(k)=1.*nH(3)/maxi(k)
 f4H(k)=1.*nH(4)/maxi(k)
 f5H(k)=1.*nH(5)/maxi(k)
 f6H(k)=1.*nH(6)/maxi(k)

END DO !k

IF (flag.eq.1) THEN
 j=26
ELSE IF (flag.eq.2) THEN
 j=27
ELSE IF (flag.eq.3) THEN
 j=28
END IF

WRITE(j,*) SUM(f0)/nmovies, SUM(f1)/nmovies, SUM(f2)/nmovies, SUM(f3)/nmovies, SUM(f4)/nmovies, SUM(f5)/nmovies,& 
        SUM(f6)/nmovies,(SUM(f0)+SUM(f1)+SUM(f2)+SUM(f3)+SUM(f4)+SUM(f5)+SUM(f6))/nmovies

WRITE(j,*) SUM(f0O)/nmovies, SUM(f1O)/nmovies, SUM(f2O)/nmovies, SUM(f3O)/nmovies, SUM(f4O)/nmovies, SUM(f5O)/nmovies,&
        SUM(f6O)/nmovies, (SUM(f0O)+SUM(f1O)+SUM(f2O)+SUM(f3O)+SUM(f4O)+SUM(f5O)+SUM(f6O))/nmovies

WRITE(j,*) SUM(f0H)/nmovies, SUM(f1H)/nmovies, SUM(f2H)/nmovies, SUM(f3H)/nmovies, SUM(f4H)/nmovies, SUM(f5H)/nmovies,& 
        SUM(f6H)/nmovies, (SUM(f0H)+SUM(f1H)+SUM(f2H)+SUM(f3H)+SUM(f4H)+SUM(f5H)+SUM(f6H))/nmovies

CLOSE (UNIT=j)
