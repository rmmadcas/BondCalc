INTEGER::flag, maxi(nmovies)
CALL counts_molecules(flag, maxi, "diol")
CALL allocate_molecules(flag)
CALL read_atoms_diol(flag,maxi)
