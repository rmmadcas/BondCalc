IMPLICIT NONE

INTEGER::flag, lineas, i, j, k, l, Num

REAL(KIND=r):: Acopy(50,3), A(50,3)
REAL(KIND=r)::Aref(50,3), Ocopy(3), firstA(3)
REAL(KIND=r)::o(3), vn(3), Vcop(3), Vref(3), Vncop(3), Vnref(3)
REAL(KIND=r)::xxx, yyy, zzz, modu
REAL(KIND=r)::alp,dab

CHARACTER(LEN=30):: rubbish1
CHARACTER(LEN=24):: rubbish2
Acopy=A
Ocopy(1)=xxx
Ocopy(2)=yyy
Ocopy(3)=zzz
Aref=0


!leemos la geometria de la molécula
IF ((flag.eq.1).or.(flag.eq.3)) THEN
   OPEN(40,FILE="Geometry_molecule1.pdb",STATUS='OLD', ACTION='READ')
   Num=NUmMol(1)
 DO l=1,Num
        READ(40,'(A30,3F8.3,A24)')rubbish1, Aref(l,1), Aref(l,2), Aref(l,3), rubbish2
 END DO
 CLOSE (UNIT=40)
ELSE
   OPEN(41,FILE="Geometry_molecule2.pdb",STATUS='OLD', ACTION='READ')
   Num=NUmMol(2)
 DO l=1,Num
        READ(41,'(A30,3F8.3,A24)')rubbish1, Aref(l,1), Aref(l,2), Aref(l,3), rubbish2
 END DO
 CLOSE (UNIT=41)
END IF

firstA(1)=Acopy(1,1) ; firstA(2)=Acopy(1,2) ; firstA(3)=Acopy(1,3) 

!colocamos la molecula de foma que el átomo A1 coicida con el origen
DO l=1, Num
        Acopy(l,1)=Acopy(l,1)-firstA(1)
        Acopy(l,2)=Acopy(l,2)-firstA(2)
        Acopy(l,3)=Acopy(l,3)-firstA(3)
END DO
        Ocopy(1)=Ocopy(1)-firstA(1)
        Ocopy(2)=Ocopy(2)-firstA(2)
        Ocopy(3)=Ocopy(3)-firstA(3)

!ponemos el átomo Acopy2 alineado con el átomo Aref(2). Para ello hemos de rotar respecto a un eje perpendicular
! a ambos vectores, el formado por los átomos 1 y 2 del Aref de referencia y 1y2 del Acopy
!buscamos el ángulo que forma los vectores Aref(1,i)Aref(2,i) con el vector Acopy(1,i)Acopy(2,i). Como hemos puesto
!el primer átomo en el orgien, los vectores tendrán el valor del segundo átomo. Usamos el producto escalar

alp=-DACOS(((Aref(2,1)-Aref(1,1))*Acopy(2,1)+(Aref(2,2)-Aref(1,2))*Acopy(2,2)+(Aref(2,3)-Aref(1,3))*Acopy(2,3))&
 /(sqrt((Aref(2,1)-Aref(1,1))**2+(Aref(2,2)-Aref(1,2))**2+(Aref(2,3)-Aref(1,3))**2)*&
 sqrt(Acopy(2,1)**2+Acopy(2,2)**2+Acopy(2,3)**2)))

!buscamos ahora el vector normal a estos dos vectores. Para ello usamos el producto vectorial

vn(1)=(Aref(2,2)-Aref(1,2))*Acopy(2,3)-(Aref(2,3)-Aref(1,3))*Acopy(2,2)
vn(2)=(Aref(2,3)-Aref(1,3))*Acopy(2,1)-(Aref(2,1)-Aref(1,1))*Acopy(2,3)
vn(3)=(Aref(2,1)-Aref(1,1))*Acopy(2,2)-(Aref(2,2)-Aref(1,2))*Acopy(2,1)

modu=sqrt(vn(1)**2+vn(2)**2+vn(3)**2)
vn(1)=vn(1)/modu
vn(2)=vn(2)/modu
vn(3)=vn(3)/modu

!ahora hacemos una rotación de un angulo alpha respecto al vector vn
DO l=2, Num
xxx=Acopy(l,1) ; yyy=Acopy(l,2) ; zzz=Acopy(l,3)

Acopy(l,1)=(cos(alp)+vn(1)**2*(1-cos(alp)))*xxx+&
           (vn(1)*vn(2)*(1-cos(alp))-vn(3)*sin(alp))*yyy+&
           (vn(1)*vn(3)*(1-cos(alp))+vn(2)*sin(alp))*zzz

Acopy(l,2)=(vn(2)*vn(1)*(1-cos(alp))+vn(3)*sin(alp))*xxx+&
           (cos(alp)+vn(2)**2*(1-cos(alp)))*yyy+&
           (vn(2)*vn(3)*(1-cos(alp))-vn(1)*sin(alp))*zzz

Acopy(l,3)=(vn(3)*vn(1)*(1-cos(alp))-vn(2)*sin(alp))*xxx+&
           (vn(3)*vn(2)*(1-cos(alp))+vn(1)*sin(alp))*yyy+&
           (cos(alp)+vn(3)**2*(1-cos(alp)))*zzz
END DO

xxx=Ocopy(1) ; yyy=Ocopy(2) ; zzz=Ocopy(3)

Ocopy(1)=(cos(alp)+vn(1)**2*(1-cos(alp)))*xxx+&
         (vn(1)*vn(2)*(1-cos(alp))-vn(3)*sin(alp))*yyy+&
         (vn(1)*vn(3)*(1-cos(alp))+vn(2)*sin(alp))*zzz

Ocopy(2)=(vn(2)*vn(1)*(1-cos(alp))+vn(3)*sin(alp))*xxx+&
         (cos(alp)+vn(2)**2*(1-cos(alp)))*yyy+&
         (vn(2)*vn(3)*(1-cos(alp))-vn(1)*sin(alp))*zzz

Ocopy(3)=(vn(3)*vn(1)*(1-cos(alp))-vn(2)*sin(alp))*xxx+&
         (vn(3)*vn(2)*(1-cos(alp))+vn(1)*sin(alp))*yyy+&
         (cos(alp)+vn(3)**2*(1-cos(alp)))*zzz

!ya hemos alineado el vector formado por el primer y segundo átomo con la referencia, ahora giremos el tercero para colocarlo.
!para ello buscamos el angulo que forman las proyecciones del vector 1 y 3 con las proyecciones del 1y3 de la referencia
!sobre el plano que define el vetor formado los atomos 1 y 2. Buscamos las proyecciones

!el vector normal que define el plano de rotación será el formado por los átomos 1 y 2
vn(1)=Acopy(2,1)
vn(2)=Acopy(2,2)
vn(3)=Acopy(2,3)

modu=sqrt(vn(1)**2+vn(2)**2+vn(3)**2)
vn(1)=vn(1)/modu
vn(2)=vn(2)/modu
vn(3)=vn(3)/modu

!vn es el vector que une los átomos uno y dos. Es el vector que usamos como eje para realizaremos las rotaciones
!necesitamos ahora las proyecciones de los vectores formados por los átomos 1 y 3 de ambos sitemas (ref y copy)
!la proyección sobre el plano será el vector menos la proyección sobre el vector vn


Vcop(1)=Acopy(3,1)-Acopy(1,1) ; Vcop(2)=Acopy(3,2)-Acopy(1,2) ; Vcop(3)=Acopy(3,3)-Acopy(1,3)
Vref(1)=Aref(3,1)-Aref(1,1) ; Vref(2)=Aref(3,2)-Aref(1,2) ; Vref(3)=Aref(3,3)-Aref(1,3) 

!si multiplicamos vectorialmente dichos vectores por el vector normal obtendremos dos vectores contenidos en el plano pero que guardan
!los ángulos entre ellos se conservarán

Vnref(1)=(Vref(1)*vn(1)+Vref(2)*vn(2)+Vref(3)*vn(3))*vn(1)
Vnref(2)=(Vref(1)*vn(1)+Vref(2)*vn(2)+Vref(3)*vn(3))*vn(2)
Vnref(3)=(Vref(1)*vn(1)+Vref(2)*vn(2)+Vref(3)*vn(3))*vn(3)

!restamos Vnref a Vref para obtener la proyeccion sobre el plano
Vref(1)=Vref(1)-Vnref(1)
Vref(2)=Vref(2)-Vnref(2)
Vref(3)=Vref(3)-Vnref(3)
!hacemos lo mismo para Vcopy

Vncop(1)=(Vcop(1)*vn(1)+Vcop(2)*vn(2)+Vcop(3)*vn(3))*vn(1)
Vncop(2)=(Vcop(1)*vn(1)+Vcop(2)*vn(2)+Vcop(3)*vn(3))*vn(2)
Vncop(3)=(Vcop(1)*vn(1)+Vcop(2)*vn(2)+Vcop(3)*vn(3))*vn(3)

!restamos Vnref a Vcop para obtener la proyeccion sobre el plano
Vcop(1)=Vcop(1)-Vncop(1)
Vcop(2)=Vcop(2)-Vncop(2)
Vcop(3)=Vcop(3)-Vncop(3)

!buscamos el angulo alpha entre Vref y Vcop

alp=-DACOS((Vref(1)*Vcop(1)+Vref(2)*Vcop(2)+Vref(3)*Vcop(3))&
 /(sqrt(Vref(1)**2+Vref(2)**2+Vref(3)**2)*sqrt(Vcop(1)**2+Vcop(2)**2+Vcop(3)**2)))

!ahora hacemos una rotación de un angulo alpha respecto al vector vn
DO l=3, Num
xxx=Acopy(l,1) ; yyy=Acopy(l,2) ; zzz=Acopy(l,3)

Acopy(l,1)=(cos(alp)+vn(1)**2*(1-cos(alp)))*xxx+&
           (vn(1)*vn(2)*(1-cos(alp))-vn(3)*sin(alp))*yyy+&
           (vn(1)*vn(3)*(1-cos(alp))+vn(2)*sin(alp))*zzz

Acopy(l,2)=(vn(2)*vn(1)*(1-cos(alp))+vn(3)*sin(alp))*xxx+&
           (cos(alp)+vn(2)**2*(1-cos(alp)))*yyy+&
           (vn(2)*vn(3)*(1-cos(alp))-vn(1)*sin(alp))*zzz

Acopy(l,3)=(vn(3)*vn(1)*(1-cos(alp))-vn(2)*sin(alp))*xxx+&
           (vn(3)*vn(2)*(1-cos(alp))+vn(1)*sin(alp))*yyy+&
           (cos(alp)+vn(3)**2*(1-cos(alp)))*zzz
END DO

xxx=Ocopy(1) ; yyy=Ocopy(2) ; zzz=Ocopy(3)

Ocopy(1)=(cos(alp)+vn(1)**2*(1-cos(alp)))*xxx+&
         (vn(1)*vn(2)*(1-cos(alp))-vn(3)*sin(alp))*yyy+&
         (vn(1)*vn(3)*(1-cos(alp))+vn(2)*sin(alp))*zzz

Ocopy(2)=(vn(2)*vn(1)*(1-cos(alp))+vn(3)*sin(alp))*xxx+&
         (cos(alp)+vn(2)**2*(1-cos(alp)))*yyy+&
         (vn(2)*vn(3)*(1-cos(alp))-vn(1)*sin(alp))*zzz

Ocopy(3)=(vn(3)*vn(1)*(1-cos(alp))-vn(2)*sin(alp))*xxx+&
         (vn(3)*vn(2)*(1-cos(alp))+vn(1)*sin(alp))*yyy+&
         (cos(alp)+vn(3)**2*(1-cos(alp)))*zzz


IF ((Ocopy(1)**2+Ocopy(2)**2+Ocopy(3)**2).lt.(dab**2)) THEN
 DO l=1, Num
        Acopy(l,1)=Acopy(l,1)+Aref(1,1)
        Acopy(l,2)=Acopy(l,2)+Aref(1,2)
        Acopy(l,3)=Acopy(l,3)+Aref(1,3)
 END DO
        Ocopy(1)=Ocopy(1)+Aref(1,1)
        Ocopy(2)=Ocopy(2)+Aref(1,2)
        Ocopy(3)=Ocopy(3)+Aref(1,3)
 WRITE(100+flag,'(A31,3F7.3,A24)')"ATOM      7  F   MOL           ", o(1), o(2), o(3), "  1.00  0.00           F"
END IF
