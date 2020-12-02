IMPLICIT NONE

REAL(KIND=r)::Density(141,141,141)         !irá desde -7 hasta 7 pasando por 0
REAL(KIND=r)::xxx, yyy, zzz
CHARACTER(LEN=30):: rubbish1
CHARACTER(LEN=24):: rubbish2
REAL(KIND=r)::step
INTEGER::i, j, k, l 
INTEGER::intX, intY, intZ, intvalue

step=0.1

OPEN(121,file="SDF1.pdb",STATUS='OLD', ACTION='READ')
OPEN(131,file="DensityProfile1.pdb",STATUS='REPLACE', ACTION='WRITE')

Density=0
DO l=1, sdf1
        READ(121,'(A30,3F8.3,A24)') rubbish1, xxx, yyy, zzz, rubbish2
        intX=INT((xxx*10+70)+1)      
        intY=INT((yyy*10+70)+1)      
        intZ=INT((zzz*10+70)+10)
        Density(intX,intY,intZ)=Density(intX,intY,intZ)+1
END DO
        Density=Density/MAXVAL(Density)

DO i=1,141
  DO j=1, 141
    DO k=1, 141
      intvalue=int(Density(i,j,k)*1000)
      IF ((intvalue.ge.100).and.(intvalue.lt.200)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  B   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           B"
      ELSE IF ((intvalue.ge.200).and.(intvalue.lt.300)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  D   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           D"
      ELSE IF ((intvalue.ge.300).and.(intvalue.lt.400)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  E   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           E"
      ELSE IF ((intvalue.ge.400).and.(intvalue.lt.500)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  F   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           F"
      ELSE IF ((intvalue.ge.500).and.(intvalue.lt.600)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  G   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           G"
      ELSE IF ((intvalue.ge.600).and.(intvalue.lt.700)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  J   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           J"
      ELSE IF ((intvalue.ge.700).and.(intvalue.lt.800)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  K   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           K"
      ELSE IF ((intvalue.ge.800).and.(intvalue.lt.900)) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  L   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           L"
      ELSE IF (intvalue.ge.900) THEN
WRITE(131,'(A31,3F7.3,A24)')"ATOM      2  M   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           M"
      END IF
    END DO 
  END DO
END DO
CLOSE(UNIT=121)

IF (number_of_pdb.eq.2) THEN

       OPEN(122,file="SDF2.pdb",STATUS='OLD', ACTION='READ')
       OPEN(132,file="DensityProfile2.dat",STATUS='REPLACE', ACTION='WRITE')

       Density=0
       DO l=1, sdf2
        READ(122,'(A30,3F8.3,A24)') rubbish1, xxx, yyy, zzz, rubbish2
        intX=INT((xxx*10+70)+1)      
        intY=INT((yyy*10+70)+1)      
        intZ=INT((zzz*10+70)+1)
        Density(intX,intY,intZ)=Density(intX,intY,intZ)+1
       END DO
       Density=Density/MAXVAL(Density)

       DO i=1, 141
         DO j=1, 141
           DO k=1, 141
      intvalue=int(Density(i,j,k)*1000)
      IF ((intvalue.ge.100).and.(intvalue.lt.200)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  B   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           B"
      ELSE IF ((intvalue.ge.200).and.(intvalue.lt.300)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  D   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           D"
      ELSE IF ((intvalue.ge.300).and.(intvalue.lt.400)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  E   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           E"
      ELSE IF ((intvalue.ge.400).and.(intvalue.lt.500)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  F   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           F"
      ELSE IF ((intvalue.ge.500).and.(intvalue.lt.600)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  G   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           G"
      ELSE IF ((intvalue.ge.600).and.(intvalue.lt.700)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  J   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           J"
      ELSE IF ((intvalue.ge.700).and.(intvalue.lt.800)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  K   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           K"
      ELSE IF ((intvalue.ge.800).and.(intvalue.lt.900)) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  L   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           L"
      ELSE IF (intvalue.ge.900) THEN
WRITE(132,'(A31,3F7.3,A24)')"ATOM      2  M   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           M"
      END IF
           END DO
         END DO
       END DO
       CLOSE(UNIT=122)


       OPEN(123,file="SDF1-2.pdb",STATUS='OLD', ACTION='READ')
       OPEN(133,file="DensityProfile1-2.dat",STATUS='REPLACE', ACTION='WRITE')

       Density=0
       DO l=1, sdf12
        READ(123,'(A30,3F8.3,A24)') rubbish1, xxx, yyy, zzz, rubbish2
        intX=INT((xxx*10+70)+1)      
        intY=INT((yyy*10+70)+1)      
        intZ=INT((zzz*10+70)+1)
        Density(intX,intY,intZ)=Density(intX,intY,intZ)+1
       END DO
        Density=Density/MAXVAL(Density)

       DO i=1, 141
         DO j=1, 141
           DO k=1, 141
      intvalue=int(Density(i,j,k)*1000)
      IF ((intvalue.ge.100).and.(intvalue.lt.200)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  B   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           B"
      ELSE IF ((intvalue.ge.200).and.(intvalue.lt.300)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  D   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           D"
      ELSE IF ((intvalue.ge.300).and.(intvalue.lt.400)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  E   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           E"
      ELSE IF ((intvalue.ge.400).and.(intvalue.lt.500)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  F   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           F"
      ELSE IF ((intvalue.ge.500).and.(intvalue.lt.600)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  G   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           G"
      ELSE IF ((intvalue.ge.600).and.(intvalue.lt.700)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  J   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           J"
      ELSE IF ((intvalue.ge.700).and.(intvalue.lt.800)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  K   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           K"
      ELSE IF ((intvalue.ge.800).and.(intvalue.lt.900)) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  L   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           L"
      ELSE IF (intvalue.ge.900) THEN
WRITE(133,'(A31,3F7.3,A24)')"ATOM      2  M   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           M"
      END IF
           END DO
         END DO
       END DO
       CLOSE(UNIT=123)

       OPEN(124,file="SDF2-1.pdb",STATUS='OLD', ACTION='READ')
       OPEN(134,file="DensityProfile2-1.dat",STATUS='REPLACE', ACTION='WRITE')

       Density=0
       DO l=1, sdf21
        READ(124,'(A30,3F8.3,A24)') rubbish1, xxx, yyy, zzz, rubbish2
        intX=INT((xxx*10+70)+1)      
        intY=INT((yyy*10+70)+1)      
        intZ=INT((zzz*10+70)+1)
        Density(intX,intY,intZ)=Density(intX,intY,intZ)+1
       END DO
       Density=Density/MAXVAL(Density)

       DO i=1, 141
         DO j=1, 141
           DO k=1, 141
      intvalue=int(Density(i,j,k)*1000)
      IF ((intvalue.ge.100).and.(intvalue.lt.200)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  B   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           B"
      ELSE IF ((intvalue.ge.200).and.(intvalue.lt.300)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  D   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           D"
      ELSE IF ((intvalue.ge.300).and.(intvalue.lt.400)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  E   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           E"
      ELSE IF ((intvalue.ge.400).and.(intvalue.lt.500)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  F   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           F"
      ELSE IF ((intvalue.ge.500).and.(intvalue.lt.600)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  G   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           G"
      ELSE IF ((intvalue.ge.600).and.(intvalue.lt.700)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  J   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           J"
      ELSE IF ((intvalue.ge.700).and.(intvalue.lt.800)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  K   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           K"
      ELSE IF ((intvalue.ge.800).and.(intvalue.lt.900)) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  L   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           L"
      ELSE IF (intvalue.ge.900) THEN
WRITE(134,'(A31,3F7.3,A24)')"ATOM      2  M   MOL           ", step*i-7, step*j-7, step*k-7, "  1.00  0.00           M"
      END IF
           END DO
         END DO
       END DO
       CLOSE(UNIT=124)
END IF
