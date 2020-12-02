INTEGER::flag, maxi(nmovies)
CALL counts_molecules(flag, maxi, "gene")
CALL allocate_molecules(flag)
CALL read_atoms_gene(flag, maxi)
