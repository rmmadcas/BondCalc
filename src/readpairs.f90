IMPLICIT NONE
ini=0; inj=0
w=1
DO i=a,b
       ini(w)=initot(i)
       inj(w)=injtot(i)
       w=w+1
END DO
