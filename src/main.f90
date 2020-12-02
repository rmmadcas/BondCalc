PROGRAM BondCalc

    INCLUDE 'declareVAR.f90'
    INCLUDE 'ReadInitialVAR.f90'

    ALLOCATE(maxi1(1:nmovies))                          !numero maximo de moleculas en cada movie
    IF (number_of_pdb.eq.2) ALLOCATE(maxi2(1:nmovies))  !moleculas del segundo, solo la lee si hay dos movies

    ALLOCATE (nHBpmolec1(1:nmovies), avgnHB1(1:nmovies),&
        sumnHB1(1:nmovies))

    IF (number_of_pdb.eq.2) THEN
     ALLOCATE (nHBpmolec2(1:nmovies), avgnHB2(1:nmovies),&
        sumnHB2(1:nmovies))
     ALLOCATE (nHBpmolecm(1:nmovies), avgnHBm(1:nmovies),&
        sumnHBm(1:nmovies))
    END IF


!**************************** iniciamos variables ***********************************************
    OPEN(17,file="All-bonds1.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(17,'(A42)') "#Movie  i  j   nHB/nPi   dOHin   condition"
    
    IF (SDF.EQ."y") THEN
     OPEN(101, file="SDF1.pdb",STATUS='REPLACE', ACTION='WRITE')
     IF (number_of_pdb.eq.2) THEN   
      OPEN(102, file="SDF2.pdb",STATUS='REPLACE', ACTION='WRITE')
      OPEN(103, file="SDF1-2.pdb",STATUS='REPLACE', ACTION='WRITE')
      OPEN(104, file="SDF2-1.pdb",STATUS='REPLACE', ACTION='WRITE')
     END IF
    END IF

    IF (number_of_pdb.eq.2) THEN                                !solo si hay dos movies
        OPEN(18,file="All-bonds2.dat",STATUS='NEW', ACTION='WRITE')
        WRITE(18,'(A42)') "#Movie  i  j   nHB/nPi   dOHin   condition"
        OPEN(19,file="All-bondsMix.dat",STATUS='NEW', ACTION='WRITE')
        WRITE(19,'(A42)') "#Movie  i  j   nHB/nPi   dOHin   condition"
    END IF

    ALLOCATE(lxxx(1:nmovies),lyyy(1:nmovies),lzzz(1:nmovies))
    CALL readCellParameter                             !lee los parametros de celda
    CALL read_movie1                                    !lee la movie 1
    maxnum1=MAXVAL(maxi1)
    ALLOCATE(enl1(1:nmovies,1:maxnum1),enlO1(1:nmovies,1:maxnum1),enlH1(1:nmovies,1:maxnum1))

    IF (life_time.eq."y") THEN
     ALLOCATE (matrix1(1:nmovies,1:maxnum1,1:maxnum1))
     ALLOCATE (old(1:nmovies,1:maxnum1,1:maxnum1))
     ALLOCATE (new(1:nmovies,1:maxnum1,1:maxnum1))
    END IF

    IF (number_of_pdb.eq.2) THEN
     CALL read_movie2            !lee la movie 2, si existe
     maxnum2=MAXVAL(maxi2)
     maxnumMix=maxnum1+maxnum2
     ALLOCATE(enl2(1:nmovies,1:maxnum2),enlO2(1:nmovies,1:maxnum2),enlH2(1:nmovies,1:maxnum2))
     ALLOCATE(enlm(1:nmovies,1:maxnumMix),enlOm(1:nmovies,1:maxnumMix),enlHm(1:nmovies,1:maxnumMix))

     IF (life_time.eq."y") THEN
      ALLOCATE (matrix2(1:nmovies,1:maxnum2,1:maxnum2))
      ALLOCATE (matrixm(1:nmovies,1:maxnum1,1:maxnum2))
     END IF 
    END IF

     nmoviestotal = nmovies
     ALLOCATE (naggregates(1:nmoviestotal), maxilinepmov(1:nmoviestotal))
     ALLOCATE (maxilinepmov1(1:nmoviestotal))
     IF(number_of_pdb.eq.2) ALLOCATE (maxilinepmov2(1:nmoviestotal),maxilinepmov3(1:nmoviestotal))

!**************************** iniciamos variables ***********************************************
    sdf1=0
    sdf2=0
    sdf12=0
    sdf21=0

    nHBpmolec1=0
    enl1=0 
    enlO1=0 
    enlH1=0 
    IF (life_time.eq."y") THEN
     matrix1=0
    END IF

    IF (number_of_pdb.eq.2) THEN
     nHBpmolec2=0
     enl2=0
     enlO2=0
     enlH2=0
     IF (life_time.eq."y") THEN
     matrix2=0
     END IF
     
     nHBpmolecm=0
     enlm=0
     enlOm=0
     enlHm=0
     IF (life_time.eq."y") THEN
      matrixm=0
     END IF
    END IF
    first1="y"
    first2="y"
    first3="y"
    first4="y"


    CALL read_DefaultsInteractions
    IF(type_of_molecules(1).eq."gene") CALL read_Interactions

    DO k=1, nmovies
        d=0; d1=0;d2=0;d3=0
        DO i=1, maxi1(k)-1
                DO j=i+1, maxi1(k)
                   IF(type_of_molecules(1).eq."gene") THEN
                        Call Generic_interaction(type_of_molecules(1), type_of_molecules(1),1,1)
                   ELSE
                        Call Interaction(type_of_molecules(1), type_of_molecules(1),1,1)
                   END IF
                END DO !j=i+1, maxi1(k)
        END DO !i=1, maxi1(k)-1

      IF(number_of_pdb.eq.2) THEN                               !only for 2 pdb files
        DO i=1, maxi2(k)-1
                DO j=i+1, maxi2(k)
                   IF(type_of_molecules(1).eq."gene") THEN
                        Call Generic_interaction(type_of_molecules(2), type_of_molecules(2),2,2)
                   ELSE
                        Call Interaction(type_of_molecules(2), type_of_molecules(2),2,2)
                   END IF     
                END DO !j=i+1, maxi2(k)
        END DO !i=1, maxi2(k)-1

        DO i=1, maxi1(k)
                DO j=1, maxi2(k)
                   IF(type_of_molecules(1).eq."gene") THEN
                        Call Generic_interaction(type_of_molecules(1), type_of_molecules(2),1,2)
                   ELSE
                        Call Interaction(type_of_molecules(1), type_of_molecules(2),1,2)
                   END IF     
                END DO !j=i+1, maxi2(k)
        END DO !i=1, maxi2(k)-1
      END IF                                                    !only for 2 pdb files

    maxilinepmov1(k)=d1-1
    IF(number_of_pdb.eq.2) maxilinepmov2(k)=d2-1
    IF(number_of_pdb.eq.2) maxilinepmov3(k)=d3-1
    END DO ! k=1, nmovies

    CLOSE(UNIT=17)
    IF (number_of_pdb.eq.2) THEN
     CLOSE(UNIT=18)
     CLOSE(UNIT=19)
    END IF

    CALL SUMnB

    CALL nbonds(1)
    IF(number_of_pdb.eq.2) CALL nbonds(2)
    IF(number_of_pdb.eq.2) CALL nbonds(3)


    IF (clusters.eq."y") THEN
     CALL CalcClusters(1)
     IF(number_of_pdb.eq.2) CALL CalcClusters(2)
     IF(number_of_pdb.eq.2) CALL CalcClusters(3)
    END IF

    IF (life_time.eq."y") THEN
     CALL CalcLifeTime
    END IF

    IF (SDF.eq."y") THEN
     CLOSE(UNIT=101) ; CLOSE(UNIT=102) ;CLOSE(UNIT=103) ;CLOSE(UNIT=104)
     CALL colorin3D
    END IF
   
    CALL WriteOutput

!********************************************* rutina que evalua el número de enlaces*****************************************
 CONTAINS

 SUBROUTINE WriteOutput
    INCLUDE 'WriteOutput.f90'
 END SUBROUTINE WriteOutput
 
 SUBROUTINE CalcClusters(flag)
    INCLUDE 'CalcCluster.f90'
 END SUBROUTINE CalcClusters
 
 SUBROUTINE removeclusterfound
    INCLUDE 'removeclusterfound.f90'
 END SUBROUTINE removeclusterfound

 SUBROUTINE readpairs
    INCLUDE 'readpairs.f90'
 END SUBROUTINE readpairs
 
 SUBROUTINE CalcLifeTime
    INCLUDE 'CalcLifeTime.f90'
 END SUBROUTINE CalcLifeTime

 SUBROUTINE nbonds(flag)
    INCLUDE 'nbonds.f90'
 END SUBROUTINE nbonds

  SUBROUTINE SUMnB 
     INCLUDE 'SUMnB.f90'
  END SUBROUTINE SUMnB  

  SUBROUTINE allocate_molecules(flag) 
     INCLUDE 'allocate_molecules.f90'
  END SUBROUTINE allocate_molecules

  SUBROUTINE read_movie1      
     INCLUDE 'read_movie1.f90'
  END SUBROUTINE read_movie1

  SUBROUTINE read_movie2      
     INCLUDE 'read_movie2.f90'
  END SUBROUTINE read_movie2

  SUBROUTINE read_water(flag,maxi)  
     INCLUDE 'read_water.f90'
  END SUBROUTINE read_water
  
  SUBROUTINE read_ammo(flag,maxi)  
     INCLUDE 'read_ammo.f90'
  END SUBROUTINE read_ammo

  SUBROUTINE read_gene(flag,maxi)  
     INCLUDE 'read_gene.f90'
  END SUBROUTINE read_gene

  SUBROUTINE read_alcohol(flag,maxi)
     INCLUDE 'read_alcohol.f90'
  END SUBROUTINE read_alcohol

  SUBROUTINE read_diol(flag,maxi)   
     INCLUDE 'read_diol.f90'
  END SUBROUTINE read_diol

  SUBROUTINE read_ring(flag,maxi)   
     INCLUDE 'read_ring.f90'
  END SUBROUTINE read_ring

  SUBROUTINE counts_molecules(flag, maxi,type_mol)               
     INCLUDE 'counts_molecules.f90'
  END SUBROUTINE counts_molecules

  SUBROUTINE read_atoms_water(flag, maxi)
     INCLUDE 'read_atoms_water.f90'
  END SUBROUTINE read_atoms_water 

  SUBROUTINE read_atoms_gene(flag, maxi)
     INCLUDE 'read_atoms_gene.f90'
  END SUBROUTINE read_atoms_gene 

  SUBROUTINE read_atoms_ammo(flag, maxi)
     INCLUDE 'read_atoms_ammo.f90'
  END SUBROUTINE read_atoms_ammo 

  SUBROUTINE read_atoms_alcohol(flag, maxi)
     INCLUDE 'read_atoms_alcohol.f90'
  END SUBROUTINE read_atoms_alcohol

  SUBROUTINE read_atoms_diol(flag, maxi)
     INCLUDE 'read_atoms_diol.f90'
  END SUBROUTINE read_atoms_diol

  SUBROUTINE read_atoms_ring(flag, maxi)
     INCLUDE 'read_atoms_ring.f90'
  END SUBROUTINE read_atoms_ring

  SUBROUTINE create_vectors(flag, maxi) 
     INCLUDE 'create_vectors.f90'
  END SUBROUTINE create_vectors
 
  SUBROUTINE readCellParameter
     INCLUDE 'readCellParameter.f90'
  END SUBROUTINE readCellParameter

  SUBROUTINE vectorial_product(ax,ay,az,bx,by,bz,modu)
     INCLUDE 'vectorial_product.f90'
  END SUBROUTINE vectorial_product

  SUBROUTINE distance(r1,r2,dist)
     INCLUDE 'distance.f90'
  END SUBROUTINE distance

  SUBROUTINE bend
     INCLUDE 'bend.f90'
  END SUBROUTINE bend

  SUBROUTINE read_Interactions
     INCLUDE 'read_Interactions.f90'
  END SUBROUTINE

  SUBROUTINE read_DefaultsInteractions
     INCLUDE 'read_DefaultsInteractions.f90'
  END SUBROUTINE

  SUBROUTINE Interaction(mole1,mole2,flag1,flag2)
     INCLUDE 'Interaction_Puras1.f90'
     INCLUDE 'Interaction_Puras2.f90'
     INCLUDE 'Interaction_water_MIX.f90'
     INCLUDE 'Interaction_alcohol_MIX.f90'
     INCLUDE 'Interaction_diol_MIX.f90'
     INCLUDE 'Interaction_ammo_MIX.f90'
     INCLUDE 'Interaction_ring_MIX.f90'
  END SUBROUTINE Interaction
     
  SUBROUTINE Generic_Interaction(mole1,mole2,flag1,flag2)
     INCLUDE 'Generic_Interaction.f90'
  END SUBROUTINE Generic_Interaction
     
  SUBROUTINE SDF_water(O1,H1,H2,Odon,flag)
     INCLUDE 'SDF_water.f90'
  END SUBROUTINE SDF_water

  SUBROUTINE SDF_alcohol(aO1,aH1,ao,flag)
     INCLUDE 'SDF_alcohol.f90'
  END SUBROUTINE SDF_alcohol

  SUBROUTINE SDF_ring(aC,a1,a2,a3,a4,a5,a6,ao,flag)
     INCLUDE 'SDF_ring.f90'
  END SUBROUTINE SDF_ring

  SUBROUTINE SDF_ammo(N1,H1,H2,H3,Odon,flag)
     INCLUDE 'SDF_ammo.f90'
  END SUBROUTINE SDF_ammo

  SUBROUTINE SDF_gene(A,xxx,yyy,zzz,flag,dab)
     INCLUDE 'SDF_gene.f90'
  END SUBROUTINE SDF_gene

  SUBROUTINE colorin3D
     INCLUDE 'colorin3D.f90'
  END SUBROUTINE colorin3D


END PROGRAM
