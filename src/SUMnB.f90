IMPLICIT NONE
REAL (KIND=r)::na 

OPEN(20,file="AVG-nBpermolec1.dat",STATUS='NEW', ACTION='WRITE')
WRITE(20,'(A42)') "#Movie  nB_per_molec1.  AVG_nHB_per_molec1. Dispersion1."
WRITE(20,'(A42)') "#========================================="
avgnHB1=0
sumnHB1=0
disp1=0
IF (number_of_pdb.eq.2) THEN
 OPEN(21,file="AVG-nBpermolec2.dat",STATUS='NEW', ACTION='WRITE')
 WRITE(21,'(A42)') "#Movie  nB_per_molec2.  AVG_nHB_per_molec2. Dispersion2."
 WRITE(21,'(A42)') "#========================================="
 avgnHB2=0
 sumnHB2=0
 disp2=0

 OPEN(22,file="AVG-nBpermolecMix.dat",STATUS='NEW', ACTION='WRITE')
 WRITE(22,'(A42)') "#Movie  nB_per_molecm.  AVG_nHB_per_molecm. DispersionMix."
 WRITE(22,'(A42)') "#========================================="
 avgnHBm=0
 sumnHBm=0
 dispm=0
END IF

DO k=2,nmovies

   na=INT(maxi1(k))
   sumnHB1(k)=sumnHB1(k-1)+nHBpmolec1(k)/na
   avgnHB1(k)=sumnHB1(k)/(k-1.0)
   disp1=disp1+nHBpmolec1(k)/na*nHBpmolec1(k)/na
   WRITE(20,'(I5,F14.4,2F18.4)') k, nHBpmolec1(k)/na, avgnHB1(k), sqrt((disp1)/(k-1.0)-avgnHB1(k)*avgnHB1(k))

 IF (number_of_pdb.eq.2) THEN

   na=INT(maxi2(k))
   sumnHB2(k)=sumnHB2(k-1)+nHBpmolec2(k)/na
   avgnHB2(k)=sumnHB2(k)/(k-1.0)
   disp2=disp2+nHBpmolec2(k)/na*nHBpmolec2(k)/na
   WRITE(21,'(I5,F14.4,2F18.4)') k, nHBpmolec2(k)/na, avgnHB2(k), sqrt((disp2)/(k-1.0)-avgnHB2(k)*avgnHB2(k))

   na=INT(maxi1(k)+maxi2(k))
   sumnHBm(k)=sumnHBm(k-1)+nHBpmolecm(k)/na
   avgnHBm(k)=sumnHBm(k)/(k-1.0)
   dispm=dispm+nHBpmolecm(k)/na*nHBpmolecm(k)/na
   WRITE(22,'(I5,F14.4,2F18.4)') k, nHBpmolecm(k)/na, avgnHBm(k), sqrt((dispm)/(k-1.0)-avgnHBm(k)*avgnHBm(k))
 END IF

END DO


CLOSE (UNIT=20)
IF (number_of_pdb.eq.2) THEN
   CLOSE (UNIT=21)
   CLOSE (UNIT=22)
END IF
