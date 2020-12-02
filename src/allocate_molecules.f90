IMPLICIT NONE
INTEGER::flag
IF(flag.eq.1) THEN
    IF (type_of_molecules(1).eq."wate") ALLOCATE (molecule_water1(1:nmovies,1:MAXVAL(maxi1)))
    IF (type_of_molecules(1).eq."alco") ALLOCATE (molecule_alcohol1(1:nmovies,1:MAXVAL(maxi1)))
    IF (type_of_molecules(1).eq."diol") ALLOCATE (molecule_diol1(1:nmovies,1:MAXVAL(maxi1)))
    IF (type_of_molecules(1).eq."ammo") ALLOCATE (molecule_ammo1(1:nmovies,1:MAXVAL(maxi1)))
    IF (type_of_molecules(1).eq."ring") ALLOCATE (molecule_ring1(1:nmovies,1:MAXVAL(maxi1)))
    IF (type_of_molecules(1).eq."gene") ALLOCATE (molecule_gene1(1:nmovies,1:MAXVAL(maxi1)))
ELSE IF (flag.eq.2) THEN
    IF (type_of_molecules(2).eq."wate") ALLOCATE (molecule_water2(1:nmovies,1:MAXVAL(maxi2)))
    IF (type_of_molecules(2).eq."alco") ALLOCATE (molecule_alcohol2(1:nmovies,1:MAXVAL(maxi2)))
    IF (type_of_molecules(2).eq."diol") ALLOCATE (molecule_diol2(1:nmovies,1:MAXVAL(maxi2)))
    IF (type_of_molecules(2).eq."ammo") ALLOCATE (molecule_ammo2(1:nmovies,1:MAXVAL(maxi2)))
    IF (type_of_molecules(2).eq."ring") ALLOCATE (molecule_ring2(1:nmovies,1:MAXVAL(maxi2)))
    IF (type_of_molecules(2).eq."gene") ALLOCATE (molecule_gene2(1:nmovies,1:MAXVAL(maxi2)))
END IF
