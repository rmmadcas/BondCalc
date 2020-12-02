INTEGER::flag, maxi(nmovies)
CALL counts_molecules(flag, maxi, "ammo")
CALL allocate_molecules(flag)
CALL read_atoms_ammo(flag, maxi)
