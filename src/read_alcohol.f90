INTEGER::flag, maxi(nmovies)
CALL counts_molecules(flag, maxi, "alco")
CALL allocate_molecules(flag)
CALL read_atoms_alcohol(flag,maxi)
