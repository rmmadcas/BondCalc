    READ*, number_of_pdb                                !numero de pdb que voy a meter, 1 ó 2
    ALLOCATE (type_of_molecules(1:number_of_pdb))
    DO  i=1, number_of_pdb                              !typos de pdb: agua, alcohol, diol, anillos
        READ*, type_of_molecules(i)
    END DO
    DO  i=1, number_of_pdb
        IF (type_of_molecules(i).eq."gene") READ*, NumMol(i) !número de átomos de la molécula genérica
    END DO
    READ*, life_time                                    !calculamos el tiempo de vida medio del enlace?
    READ*, clusters                                     !calculamos los clusters?
    READ*, SDF                                          !calculamos la SDF?
    READ*, n1                                           !lineas del primer archivo
    IF (number_of_pdb.eq.2) READ*, n2                   !lineas del segundo, solo la lee si hay dos movies
    READ*, nmovies                                      !numero de movies del primer archivo (debe coincidir con el segundo)
    IF (type_of_molecules(1).eq."gene") THEN
        READ*, NumInteraction1                          !numero de interacciones propias de la molécula genérica 1
        NumInteraction1=2*NumInteraction1
        IF (number_of_pdb.eq.2) THEN
                READ*, NumInteraction2                  !numero de interacciones propias de la molécula genérica 2
                NumInteraction2=2*NumInteraction2
                READ*, NumInteractionMix                !numero de interacciones propias de la mezcla de moléculas genéricas
                NumInteractionMix=2*NumInteractionMix
        END IF
    END IF
