IMPLICIT NONE
REAL(KIND=r):: a1(3), a2(3),a3(3),a4(3),a5(3),a6(3),aC(3),ao(3)
REAL(KIND=r)::Cref(3),C1ref(3), C2ref(3), C3ref(3), C4ref(3), C5ref(3), C6ref(3)
REAL(KIND=r)::o(3)
REAL(KIND=r)::xxx, yyy, zzz
REAL(KIND=r)::phi, theta
INTEGER::flag

C1ref=a1
C2ref=a2
C3ref=a3
C4ref=a4
C5ref=a5
C6ref=a6
Cref=aC
o=ao

!1º trasladamos el sistema para colocar O_ref en el eje de coordenadas
C1ref=C1ref-Cref
C2ref=C2ref-Cref
C3ref=C3ref-Cref
C4ref=C4ref-Cref
C5ref=C5ref-Cref
C6ref=C6ref-Cref
o=o-Cref
Cref=0

!2º giramos el sistema para colocar C3 en el plano xy; es decir C1ref(3)=0
phi=DATAN(-C3ref(3)/C3ref(2))
!C1
xxx=C1ref(1) ; yyy=C1ref(2) ; zzz=C1ref(3)
C1ref(1)=xxx ; C1ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C1ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C2
xxx=C2ref(1) ; yyy=C2ref(2) ; zzz=C2ref(3)
C2ref(1)=xxx ; C2ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C2ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C3
xxx=C3ref(1) ; yyy=C3ref(2) ; zzz=C3ref(3)
C3ref(1)=xxx ; C3ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C3ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C4
xxx=C4ref(1) ; yyy=C4ref(2) ; zzz=C4ref(3)
C4ref(1)=xxx ; C4ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C4ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C5
xxx=C5ref(1) ; yyy=C5ref(2) ; zzz=C5ref(3)
C5ref(1)=xxx ; C5ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C5ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C6
xxx=C6ref(1) ; yyy=C6ref(2) ; zzz=C6ref(3)
C6ref(1)=xxx ; C6ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C6ref(3)=yyy*sin(phi)+zzz*cos(phi)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(phi)-zzz*sin(phi) ; o(3)=yyy*sin(phi)+zzz*cos(phi)


!3º giramos el sistema para colocar C3 en el eje  x; es decir C3ref(2)=0
theta=DATAN(-C3ref(2)/C3ref(1))

!C1
xxx=C1ref(1) ; yyy=C1ref(2) ; zzz=C1ref(3)
C1ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C1ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C1ref(3)=zzz
!C2
xxx=C2ref(1) ; yyy=C2ref(2) ; zzz=C2ref(3)
C2ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C2ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C2ref(3)=zzz
!C3
xxx=C3ref(1) ; yyy=C3ref(2) ; zzz=C3ref(3)
C3ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C3ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C3ref(3)=zzz
!C4
xxx=C4ref(1) ; yyy=C4ref(2) ; zzz=C4ref(3)
C4ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C4ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C4ref(3)=zzz
!C5
xxx=C5ref(1) ; yyy=C5ref(2) ; zzz=C5ref(3)
C5ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C5ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C5ref(3)=zzz
!C6
xxx=C6ref(1) ; yyy=C6ref(2) ; zzz=C6ref(3)
C6ref(1)=xxx*cos(theta)-yyy*sin(theta) ; C6ref(2)=xxx*sin(theta)+yyy*cos(theta) ; C6ref(3)=zzz
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx*cos(theta)-yyy*sin(theta) ; o(2)=xxx*sin(theta)+yyy*cos(theta) ; o(3)=zzz



!4º giramos el sistema para colocar C1 en el plano  xy; es decir C1ref(3)=0
phi=DATAN(-C1ref(3)/C1ref(2))
!C1
xxx=C1ref(1) ; yyy=C1ref(2) ; zzz=C1ref(3)
C1ref(1)=xxx ; C1ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C1ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C2
xxx=C2ref(1) ; yyy=C2ref(2) ; zzz=C2ref(3)
C2ref(1)=xxx ; C2ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C2ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C3
xxx=C3ref(1) ; yyy=C3ref(2) ; zzz=C3ref(3)
C3ref(1)=xxx ; C3ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C3ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C4
xxx=C4ref(1) ; yyy=C4ref(2) ; zzz=C4ref(3)
C4ref(1)=xxx ; C4ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C4ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C5
xxx=C5ref(1) ; yyy=C5ref(2) ; zzz=C5ref(3)
C5ref(1)=xxx ; C5ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C5ref(3)=yyy*sin(phi)+zzz*cos(phi)
!C6
xxx=C6ref(1) ; yyy=C6ref(2) ; zzz=C6ref(3)
C6ref(1)=xxx ; C6ref(2)=yyy*cos(phi)-zzz*sin(phi) ; C6ref(3)=yyy*sin(phi)+zzz*cos(phi)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(phi)-zzz*sin(phi) ; o(3)=yyy*sin(phi)+zzz*cos(phi)

SELECT CASE (flag)
CASE (1)
        IF(first1.eq."y") THEN
         OPEN(110+flag,file="Reference_molecule1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  C   MOL           ", C1ref(1), C1ref(2), C1ref(3), "  1.00  0.00           C"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  C   MOL           ", C2ref(1), C2ref(2), C2ref(3), "  1.00  0.00           C"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  C   MOL           ", C3ref(1), C3ref(2), C3ref(3), "  1.00  0.00           C"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      4  C   MOL           ", C4ref(1), C4ref(2), C4ref(3), "  1.00  0.00           C"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      5  C   MOL           ", C5ref(1), C5ref(2), C5ref(3), "  1.00  0.00           C"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      6  C   MOL           ", C6ref(1), C6ref(2), C6ref(3), "  1.00  0.00           C"
         CLOSE(UNIT=110+flag)
        first1="n"
        END IF
CASE (2)
        IF(first2.eq."y") THEN
         OPEN(110+flag,file="Reference_molecule2.pdb",STATUS='REPLACE', ACTION='WRITE') 
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      1  C   MOL           ", C1ref(1), C1ref(2), C1ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      2  C   MOL           ", C2ref(1), C2ref(2), C2ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      3  C   MOL           ", C3ref(1), C3ref(2), C3ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      4  C   MOL           ", C4ref(1), C4ref(2), C4ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      5  C   MOL           ", C5ref(1), C5ref(2), C5ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      6  C   MOL           ", C6ref(1), C6ref(2), C6ref(3), "  1.00  0.00           C"
         CLOSE(UNIT=110+flag)
        first2="n"
        END IF
CASE (3)
        IF(first3.eq."y") THEN
         OPEN(110+flag,file="Reference_molecule1-2.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      1  C   MOL           ", C1ref(1), C1ref(2), C1ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      2  C   MOL           ", C2ref(1), C2ref(2), C2ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      3  C   MOL           ", C3ref(1), C3ref(2), C3ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      4  C   MOL           ", C4ref(1), C4ref(2), C4ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      5  C   MOL           ", C5ref(1), C5ref(2), C5ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      6  C   MOL           ", C6ref(1), C6ref(2), C6ref(3), "  1.00  0.00           C"
         CLOSE(UNIT=110+flag)
        first3="n"
        END IF
CASE (4)
        IF(first4.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule2-1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      1  C   MOL           ", C1ref(1), C1ref(2), C1ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      2  C   MOL           ", C2ref(1), C2ref(2), C2ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      3  C   MOL           ", C3ref(1), C3ref(2), C3ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      4  C   MOL           ", C4ref(1), C4ref(2), C4ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      5  C   MOL           ", C5ref(1), C5ref(2), C5ref(3), "  1.00  0.00           C"
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      6  C   MOL           ", C6ref(1), C6ref(2), C6ref(3), "  1.00  0.00           C"
         CLOSE(UNIT=110+flag)
        first4="n"
        END IF
END SELECT

IF (o(1)*o(1)+o(2)*o(2)+o(3)*o(3).lt.dCD*dCD) THEN
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      7  F   MOL           ", o(1), o(2), o(3), "  1.00  0.00           F"
 SELECT CASE (flag)
 CASE (1)
  sdf1=sdf1+1
 CASE (2)
  sdf2=sdf2+1
 CASE (3)
  sdf12=sdf12+1
 CASE (4)
  sdf21=sdf21+1
 END SELECT
END IF
