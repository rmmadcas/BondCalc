IMPLICIT NONE

REAL(KIND=r)::N1(3), H1(3), H2(3), H3(3), Odon(3)
REAL(KIND=r)::Nref(3), H1ref(3), H2ref(3), H3ref(3)
REAL(KIND=r)::o(3), centro(3)
REAL(KIND=r)::xxx, yyy, zzz
REAL(KIND=r)::phi, theta, alpha
INTEGER::flag

Nref=N1
H1ref=H1
H2ref=H2
H3ref=H3
o=Odon
centro=(H1+H2+H3)/3
!1º trasladamos el sistema para colocar O_ref en el eje de coordenadas
H1ref=H1ref-Nref
H2ref=H2ref-Nref
H3ref=H3ref-Nref
o=o-Nref
centro=centro-Nref
Nref=0

!2º giramos el sistema para colocar Centro en el plano xy; es decir centro(3)=0
phi=DATAN(-centro(3)/centro(2))
!H1
xxx=H1ref(1) ; yyy=H1ref(2) ; zzz=H1ref(3)
H1ref(1)=xxx ; H1ref(2)=yyy*cos(phi)-zzz*sin(phi) ; H1ref(3)=yyy*sin(phi)+zzz*cos(phi)
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx ; H2ref(2)=yyy*cos(phi)-zzz*sin(phi) ; H2ref(3)=yyy*sin(phi)+zzz*cos(phi)
!H3
xxx=H3ref(1) ; yyy=H3ref(2) ; zzz=H3ref(3)
H3ref(1)=xxx ; H3ref(2)=yyy*cos(phi)-zzz*sin(phi) ; H3ref(3)=yyy*sin(phi)+zzz*cos(phi)
!centro
xxx=centro(1) ; yyy=centro(2) ; zzz=centro(3)
centro(1)=xxx ; centro(2)=yyy*cos(phi)-zzz*sin(phi) ; centro(3)=yyy*sin(phi)+zzz*cos(phi)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(phi)-zzz*sin(phi) ; o(3)=yyy*sin(phi)+zzz*cos(phi)

!3º giramos el sistema para colocar centro en el eje x; es decir, centro(2)=0
theta=DATAN(-centro(2)/centro(1))
!H1
xxx=H1ref(1) ; yyy=H1ref(2) ; zzz=H1ref(3)
H1ref(1)=xxx*cos(theta)-yyy*sin(theta) ; H1ref(2)=xxx*sin(theta)+yyy*cos(theta) ; H1ref(3)=zzz
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx*cos(theta)-yyy*sin(theta) ; H2ref(2)=xxx*sin(theta)+yyy*cos(theta) ; H2ref(3)=zzz
!H3
xxx=H3ref(1) ; yyy=H3ref(2) ; zzz=H3ref(3)
H3ref(1)=xxx*cos(theta)-yyy*sin(theta) ; H3ref(2)=xxx*sin(theta)+yyy*cos(theta) ; H3ref(3)=zzz
!centro
xxx=centro(1) ; yyy=centro(2) ; zzz=centro(3)
centro(1)=xxx*cos(theta)-yyy*sin(theta) ; centro(2)=xxx*sin(theta)+yyy*cos(theta) ; centro(3)=zzz
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx*cos(theta)-yyy*sin(theta) ; o(2)=xxx*sin(theta)+yyy*cos(theta) ; o(3)=zzz

!4º giramos el sistema para colocar H1 de forma que H1(2)=0
alpha=DATAN(-H2ref(2)/H2ref(3))
!H1
xxx=H1ref(1) ; yyy=H1ref(2) ; zzz=H1ref(3)
H1ref(1)=xxx ; H1ref(2)=yyy*cos(alpha)-zzz*sin(alpha) ; H1ref(3)=yyy*sin(alpha)+zzz*cos(alpha)
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx ; H2ref(2)=yyy*cos(alpha)-zzz*sin(alpha) ; H2ref(3)=yyy*sin(alpha)+zzz*cos(alpha)
!H3
xxx=H3ref(1) ; yyy=H3ref(2) ; zzz=H3ref(3)
H3ref(1)=xxx ; H3ref(2)=yyy*cos(alpha)-zzz*sin(alpha) ; H3ref(3)=yyy*sin(alpha)+zzz*cos(alpha)
!centro
xxx=centro(1) ; yyy=centro(2) ; zzz=centro(3)
centro(1)=xxx ; centro(2)=yyy*cos(alpha)-zzz*sin(alpha) ; centro(3)=yyy*sin(alpha)+zzz*cos(alpha)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(alpha)-zzz*sin(alpha) ; o(3)=yyy*sin(alpha)+zzz*cos(alpha)

IF(H1ref(1).lt.0) THEN
        H1ref(1)=-H1ref(1)
        H2ref(1)=-H2ref(1)
        H3ref(1)=-H3ref(1)
        o(1)=-o(1)
END IF
IF(H2ref(2).lt.0) THEN
        H2ref(2)=-H2ref(2)
        H2ref(3)=-H2ref(3)
        o(2)=-o(2)
END IF

SELECT CASE (flag)
CASE (1)
        IF(first1.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  N   MOL           ", Nref(1), Nref(2), Nref(3), "  1.00  0.00           N"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      4  H   MOL           ", H3ref(1), H3ref(2), H3ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first1="n"
        END IF
CASE (2)
        IF(first2.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule2.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  N   MOL           ", Nref(1), Nref(2), Nref(3), "  1.00  0.00           N"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      4  H   MOL           ", H3ref(1), H3ref(2), H3ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first2="n"
        END IF
CASE (3)
        IF(first3.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule1-2.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  N   MOL           ", Nref(1), Nref(2), Nref(3), "  1.00  0.00           N"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      4  H   MOL           ", H3ref(1), H3ref(2), H3ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first3="n"
        END IF
CASE (4)
        IF(first4.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule2-1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  N   MOL           ", Nref(1), Nref(2), Nref(3), "  1.00  0.00           N"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      4  H   MOL           ", H3ref(1), H3ref(2), H3ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first4="n"
        END IF
END SELECT

IF (o(1)*o(1)+o(2)*o(2)+o(3)*o(3).lt.dOOmin*dOOmin) THEN
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
