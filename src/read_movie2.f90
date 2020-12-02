INTEGER:: flag
flag=2
IF (type_of_molecules(2).eq."wate") call read_water(flag, maxi2)
IF (type_of_molecules(2).eq."alco") call read_alcohol(flag, maxi2)
IF (type_of_molecules(2).eq."diol") call read_diol(flag, maxi2)
IF (type_of_molecules(2).eq."ring") call read_ring(flag, maxi2)
IF (type_of_molecules(2).eq."gene") call read_gene(flag, maxi2)
IF (type_of_molecules(2).eq."ammo") call read_ammo(flag, maxi2)
