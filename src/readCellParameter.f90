OPEN(16,FILE="cell_parameters.dat",STATUS='OLD', ACTION='READ')
DO k=1,nmovies    !recorro las movies
        READ(16,*) lxxx(k), lyyy(k), lzzz(k)
END DO
CLOSE (UNIT=16)
