IMPLICIT NONE

REAL(KIND=r)::O1(3), H1(3), H2(3), Odon(3)
REAL(KIND=r)::Oref(3), H1ref(3), H2ref(3)
REAL(KIND=r)::o(3)
REAL(KIND=r)::xxx, yyy, zzz
REAL(KIND=r)::phi, theta, alpha
INTEGER::flag
Oref=O1
H1ref=H1
H2ref=H2
o=Odon
!1º trasladamos el sistema para colocar O_ref en el eje de coordenadas
H1ref=H1ref-Oref
H2ref=H2ref-Oref
o=o-Oref
Oref=0

!2º giramos el sistema para colocar H1 en el plano xy; es decir H1ref(3)=0
phi=DATAN(-H1ref(3)/H1ref(2))
!H1
xxx=H1ref(1) ; yyy=H1ref(2) ; zzz=H1ref(3)
H1ref(1)=xxx ; H1ref(2)=yyy*cos(phi)-zzz*sin(phi) ; H1ref(3)=yyy*sin(phi)+zzz*cos(phi)
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx ; H2ref(2)=yyy*cos(phi)-zzz*sin(phi) ; H2ref(3)=yyy*sin(phi)+zzz*cos(phi)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(phi)-zzz*sin(phi) ; o(3)=yyy*sin(phi)+zzz*cos(phi)

!3º giramos el sistema para colocar H1 en el eje x; es decir, H1ref(2)=0
theta=DATAN(-H1ref(2)/H1ref(1))
!H1
xxx=H1ref(1) ; yyy=H1ref(2) ; zzz=H1ref(3)
H1ref(1)=xxx*cos(theta)-yyy*sin(theta) ; H1ref(2)=xxx*sin(theta)+yyy*cos(theta) ; H1ref(3)=zzz
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx*cos(theta)-yyy*sin(theta) ; H2ref(2)=xxx*sin(theta)+yyy*cos(theta) ; H2ref(3)=zzz
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx*cos(theta)-yyy*sin(theta) ; o(2)=xxx*sin(theta)+yyy*cos(theta) ; o(3)=zzz

!4º giramos el sistema para colocar H2 en el plano xy
alpha=DATAN(-H2ref(3)/H2ref(2))
!H2
xxx=H2ref(1) ; yyy=H2ref(2) ; zzz=H2ref(3)
H2ref(1)=xxx ; H2ref(2)=yyy*cos(alpha)-zzz*sin(alpha) ; H2ref(3)=yyy*sin(alpha)+zzz*cos(alpha)
!o
xxx=o(1) ; yyy=o(2) ; zzz=o(3)
o(1)=xxx ; o(2)=yyy*cos(alpha)-zzz*sin(alpha) ; o(3)=yyy*sin(alpha)+zzz*cos(alpha)

IF(H1ref(1).lt.0) THEN
        H1ref(1)=-H1ref(1)
        H2ref(1)=-H2ref(1)
        o(1)=-o(1)
END IF
IF(H2ref(2).lt.0) THEN
        H2ref(2)=-H2ref(2)
        o(2)=-o(2)
END IF


SELECT CASE (flag)
CASE (1)
        IF(first1.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  O   MOL           ", Oref(1), Oref(2), Oref(3), "  1.00  0.00           O"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first1="n"
        END IF
CASE (2)
        IF(first2.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule2.pdb",STATUS='REPLACE', ACTION='WRITE') 
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  O   MOL           ", Oref(1), Oref(2), Oref(3), "  1.00  0.00           O"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first2="n"
        END IF
CASE (3)
        IF(first3.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule1-2.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  O   MOL           ", Oref(1), Oref(2), Oref(3), "  1.00  0.00           O"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
 CLOSE(UNIT=110+flag)
        first3="n"
        END IF
CASE (4)
        IF(first4.eq."y") THEN
 OPEN(110+flag,file="Reference_molecule2-1.pdb",STATUS='REPLACE', ACTION='WRITE')
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      1  O   MOL           ", Oref(1), Oref(2), Oref(3), "  1.00  0.00           O"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      2  H   MOL           ", H1ref(1), H1ref(2), H1ref(3), "  1.00  0.00           H"
 WRITE(110+flag,'(A31,3F7.3,A24)')"ATOM      3  H   MOL           ", H2ref(1), H2ref(2), H2ref(3), "  1.00  0.00           H"
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
