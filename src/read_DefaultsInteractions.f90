IMPLICIT NONE

CHARACTER(LEN=4):: atom

OPEN(45,FILE="Default_Interactions.def",STATUS='OLD', ACTION='READ')
        READ(45,*)        !coment #Default Interacction (DO NOT MODIFY THE ORDER OF THE PARAMETERS)
        READ(45,*)        !coment #PURES

        READ(45,*)        !coment #water
        READ(45,*) ATOM, ATOM, dOOw
        READ(45,*) ATOM, ATOM, dOHw
                
        READ(45,*)        !coment #alcohol
        READ(45,*) ATOM, ATOM, dOOal
        READ(45,*) ATOM, ATOM, dOHal
                
        READ(45,*)        !coment #diol
        READ(45,*) ATOM, ATOM, dOOd
        READ(45,*) ATOM, ATOM, dOHd
                
        READ(45,*)        !coment #ammonia
        READ(45,*) ATOM, ATOM, dNN
        READ(45,*) ATOM, ATOM, dNH
                
        READ(45,*)        !coment #aromatic ring
        READ(45,*) ATOM, ATOM, dFFr
        READ(45,*) ATOM, ATOM, dFEr
                
        READ(45,*)        !coment #MIXTURE
        READ(45,*)        !coment #water-alcohol
        READ(45,*) ATOM, ATOM, dOOwal
        READ(45,*) ATOM, ATOM, dOHwal
                
        READ(45,*)        !coment #water-ammonia
        READ(45,*) ATOM, ATOM, dONwam
        READ(45,*) ATOM, ATOM, dOHwam
                
        READ(45,*)        !coment #water-diol
        READ(45,*) ATOM, ATOM, dOOwd
        READ(45,*) ATOM, ATOM, dOHwd
                
        READ(45,*)        !coment #alcohol-alcohol
        READ(45,*) ATOM, ATOM, dOOalal
        READ(45,*) ATOM, ATOM, dOHalal
                
        READ(45,*)        !coment #alcohol-diol
        READ(45,*) ATOM, ATOM, dOOald
        READ(45,*) ATOM, ATOM, dOHald
                
        READ(45,*)        !coment #alcohol-ammonia
        READ(45,*) ATOM, ATOM, dONalam
        READ(45,*) ATOM, ATOM, dOHalam
                
        READ(45,*)        !coment #diol-diol
        READ(45,*) ATOM, ATOM, dOOdd
        READ(45,*) ATOM, ATOM, dOHdd
                
        READ(45,*)        !coment #diol-ammonia
        READ(45,*) ATOM, ATOM, dONdam
        READ(45,*) ATOM, ATOM, dOHdam
                
        READ(45,*)        !coment #aromatic-aromatic
        READ(45,*) ATOM, ATOM, dFFrr
        READ(45,*) ATOM, ATOM, dFErr

        CLOSE(UNIT=45)
