IMPLICIT NONE
INTEGER:: a, b, i, k, sumita
CHARACTER(LEN=40):: cosa
        OPEN(29,FILE="Generic_Interactions.def",STATUS='OLD', ACTION='READ')
        i=1
        sumita=NumInteraction1+NumInteraction2+NumInteractionMix
        ALLOCATE(list_int(1:sumita,1:3))
        list_int=0
        Generic_Interaction1=0
        Generic_Interaction2=0
        Generic_InteractionMix=0

        READ(29,*)        !coment #Atom_i         Atom_j          CutOff
        READ(29,*)        !coment #Interactions between molecules of the first pdb

        DO k=1, NumInteraction1/2
                READ(29,*)      !NewInteraction1
                READ(29,*) a, b, Generic_Interaction1(a,b)
                Generic_Interaction1(b,a)=Generic_Interaction1(a,b)
                list_int(i,1)=a
                list_int(i,2)=b
                list_int(i,3)=1
                i=i+1
                READ(29,*) a, b, Generic_Interaction1(a,b)
                Generic_Interaction1(b,a)=Generic_Interaction1(a,b)
                list_int(i,1)=a
                list_int(i,2)=b
                list_int(i,3)=1
                i=i+1
        END DO  

        IF (number_of_pdb.eq.2) THEN
                READ(29,*)              !coment #Interactions between molecules of the second pdb
                DO k=1, NumInteraction2/2
                        READ(29,*)      !NewInteraction2
                        READ(29,*) a, b, Generic_Interaction2(a,b)
                        Generic_Interaction2(b,a)=Generic_Interaction2(a,b)
                        list_int(i,1)=a
                        list_int(i,2)=b
                        list_int(i,3)=2
                        i=i+1
                        READ(29,*) a, b, Generic_Interaction2(a,b)
                        Generic_Interaction2(b,a)=Generic_Interaction2(a,b)
                        list_int(i,1)=a
                        list_int(i,2)=b
                        list_int(i,3)=2
                        i=i+1
                END DO  

                READ(29,*)              !coment #Interactions between mixture of molecules
                DO k=1, NumInteractionMix/2
                        READ(29,*)      !NewInteractionMix
                        READ(29,*) a, b, Generic_InteractionMix(a,b)
                        Generic_InteractionMix(b,a)=Generic_InteractionMix(a,b)
                        list_int(i,1)=a
                        list_int(i,2)=b
                        list_int(i,3)=3
                        i=i+1
                        READ(29,*) a, b, Generic_InteractionMix(a,b)
                        Generic_InteractionMix(b,a)=Generic_InteractionMix(a,b)
                        list_int(i,1)=a
                        list_int(i,2)=b
                        list_int(i,3)=3
                        i=i+1
                END DO  
        END IF
        CLOSE(UNIT=29)
