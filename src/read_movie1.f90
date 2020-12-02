INTEGER:: flag
flag=1
IF (type_of_molecules(1).eq."wate") call read_water(flag,maxi1)
IF (type_of_molecules(1).eq."alco") call read_alcohol(flag,maxi1)
IF (type_of_molecules(1).eq."diol") call read_diol(flag,maxi1)
IF (type_of_molecules(1).eq."ring") call read_ring(flag,maxi1)
IF (type_of_molecules(1).eq."gene") call read_gene(flag,maxi1)
IF (type_of_molecules(1).eq."ammo") call read_ammo(flag,maxi1)
