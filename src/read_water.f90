INTEGER::flag, maxi(nmovies)
CALL counts_molecules(flag, maxi, "wate")
CALL allocate_molecules(flag)
CALL read_atoms_water(flag, maxi)
