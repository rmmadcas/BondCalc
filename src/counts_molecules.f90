IMPLICIT NONE
INTEGER::flag, model, a, i, maxi(nmovies), lineas, num
CHARACTER(LEN=10)::filename
CHARACTER(LEN=4):: type_mol
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
maxi=0
model=0
a=0
DO i=1, lineas
        READ(num,'(A4)')atom
        IF(ATOM.eq."MODE") model=model+1
        IF(ATOM.eq."ATOM") a=a+1
        IF(ATOM.eq."ENDM") THEN
                maxi(model)=a
                a=0
        END IF
END DO
SELECT CASE (type_mol)
  CASE ("wate")                !numero de moleculas de agua en cada movie
        maxi=maxi/3
  CASE ("alco")              !numero de moleculas de alchol en cada movie
        maxi=maxi/2
  CASE ("diol")                 !numero de moleculas de diol en cada movie
        maxi=maxi/4 
  CASE ("ammo")                 !numero de moleculas de amoniaco en cada movie
        maxi=maxi/4 
  CASE ("ring")                 !numero de anillos aromaticos en cada movie
        maxi=maxi/6
  CASE ("gene")                 !numero de moleculas genéricas en cada movie
        maxi=maxi/NumMol(flag)
END SELECT

!IF (nmovies.ne.model) STOP 'PROGRAM FAILED'

CLOSE (UNIT=num)
