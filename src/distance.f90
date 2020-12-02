INTEGER, PARAMETER:: r=SELECTED_REAL_KIND(P=15)
REAL (KIND=r)::xf,yf,zf,dist
REAL (KIND=r)::r1(3),r2(3)
REAL (KIND=r)::lx, ly, lz

xf=ABS(r1(1)-r2(1))
yf=ABS(r1(2)-r2(2))
zf=ABS(r1(3)-r2(3))

lx=lxxx(k)
ly=lyyy(k)
lz=lzzz(k)

IF (xf.GT.(lx/2.0)) xf=ABS(xf-lx)
IF (yf.GT.(ly/2.0)) yf=ABS(yf-ly)
IF (zf.GT.(lz/2.0)) zf=ABS(zf-lz)

dist=SQRT(xf*xf+yf*yf+zf*zf)
