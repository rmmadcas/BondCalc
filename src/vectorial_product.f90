INTEGER, PARAMETER:: r=SELECTED_REAL_KIND(P=15)
REAL (KIND=r)::ax,ay,az,bx,by,bz
REAL (KIND=r)::a, b, c, modu

a=bz*ay-by*az
b=bx*az-bz*ax
c=by*ax-bx*ay

modu=ABS(SQRT(a*a+b*b+c*c))
