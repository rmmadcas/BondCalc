IMPLICIT NONE
INTEGER, PARAMETER:: r=SELECTED_REAL_KIND(P=15)
REAL, PARAMETER:: pi=4*atan(1.)
INTEGER::i, j, k, a, n, n1, n2, number_of_pdb, m, mm, f
INTEGER::maxnum1, maxnum2, c, maxnumMix
INTEGER::n_atoms, n_molecules, nmovies, n0
CHARACTER (LEN=75):: firstline
CHARACTER (LEN=4):: atom
CHARACTER(LEN=10):: filename
CHARACTER(LEN=4), ALLOCATABLE:: type_of_molecules(:)
CHARACTER:: life_time, clusters, SDF, first1, first2, first3, first4
INTEGER:: sdf1, sdf2, sdf12, sdf21

REAL (KIND=r):: dOHmin, dOOmin, theta_min
REAL (KIND=r):: vec_aux_a(3), vec_aux_b(3), norm
REAL (KIND=r):: dCD, dangleCD, dFE, dangleFE, dist, modu
REAL (KIND=r):: dOO, dOH, dOHin, theta, costheta
REAL (KIND=r), ALLOCATABLE:: nPipmolec(:)
INTEGER, ALLOCATABLE:: maxi1(:), maxi2(:)
REAL (KIND=r), ALLOCATABLE::avgnPi(:), sumnPi(:)
REAL (KIND=r), ALLOCATABLE::lxxx(:), lyyy(:), lzzz(:)
INTEGER, ALLOCATABLE:: enl1(:,:), enlO1(:,:), enlH1(:,:)   !Matriz de enlaces
INTEGER, ALLOCATABLE:: enl2(:,:), enlO2(:,:), enlH2(:,:)   !Matriz de enlaces
INTEGER, ALLOCATABLE:: enlm(:,:), enlOm(:,:), enlHm(:,:)   !Matriz de enlaces
INTEGER, ALLOCATABLE:: matrix1(:,:,:), new(:,:,:), old(:,:,:)
INTEGER, ALLOCATABLE:: matrix2(:,:,:), matrixm(:,:,:)
REAL (KIND=r):: mu, dispersion, acumulate, real_aux!, aux
REAL (KIND=r):: lx, ly, lz
REAL (KIND=r), ALLOCATABLE:: nHBpmolec1(:), avgnHB1(:), sumnHB1(:)
REAL (KIND=r), ALLOCATABLE:: nHBpmolec2(:), avgnHB2(:), sumnHB2(:)
REAL (KIND=r), ALLOCATABLE:: nHBpmolecm(:), avgnHBm(:), sumnHBm(:)
REAL (KIND=r)::disp1, disp2, dispm 

!Variables de los Clusters
INTEGER:: nn, nmoviestotal, nn1, nn2, nn3
INTEGER:: l, h, s, cntaggreg, cntmov, freepairs, b, w
INTEGER:: nlinestotal, d, d1, d2, d3
INTEGER:: maxfreepairs, lineread
INTEGER, ALLOCATABLE:: maxilinepmov(:), maxilinepmov1(:),maxilinepmov2(:),maxilinepmov3(:)
INTEGER, ALLOCATABLE::ini(:), inj(:)
INTEGER, ALLOCATABLE::initot(:), injtot(:)
INTEGER, ALLOCATABLE:: clus(:),naggregates(:)
REAL (KIND=r):: avgaggregates, percaggregates     

!Variables para las interacciones por defecto
REAL (KIND=r)::dOOw, dOHw, dOOal, dOHal, dNN, dNH, dOOd, dOHd, dFFr, dFEr 
REAL (KIND=r)::dOOwal, dOHwal, dONwam, dOHwam, dOOwd, dOHwd 
REAL (KIND=r)::dOOalal, dOHalal, dOOald, dOHald, dONalam, dOHalam 
REAL (KIND=r)::dOOdd, dOHdd, dONdam, dOHdam, dFFrr, dFErr 

!Variables para interacciones genericas
INTEGER::NumMol(2), NumInteraction1, NumInteraction2, NumInteractionMix
INTEGER, ALLOCATABLE::list_int(:,:)
REAL (KIND=r):: Generic_Interaction1(50,50), Generic_Interaction2(50,50), Generic_InteractionMix(50,50)



!Definimos las caracteristicas de los átomos en una estructura
type:: type_water
        REAL (KIND=r):: O(3), H1(3), H2(3)
end type

type:: type_alcohol
        REAL (KIND=r):: O(3), H(3) 
end type

type:: type_diol
        REAL (KIND=r):: O1(3), H1(3)
        REAL (KIND=r):: O2(3), H2(3)
end type

type:: type_ring
        REAL (KIND=r):: C1(3), C2(3), C3(3)
        REAL (KIND=r):: C4(3), C5(3), C6(3)
        REAL (KIND=r):: Centro(3), v_nor(3)
end type

type:: type_generic
        REAL (KIND=r):: A(50,3)
end type

type:: type_ammo
        REAL (KIND=r):: N(3), H1(3), H2(3), H3(3)
end type

type(type_water), allocatable :: molecule_water1(:,:), molecule_water2(:,:) 
type(type_alcohol), allocatable :: molecule_alcohol1(:,:), molecule_alcohol2(:,:)
type(type_diol), allocatable :: molecule_diol1(:,:), molecule_diol2(:,:)
type(type_ammo), allocatable :: molecule_ammo1(:,:), molecule_ammo2(:,:)
type(type_ring), allocatable::molecule_ring1(:,:), molecule_ring2(:,:)
type(type_generic), allocatable:: molecule_gene1(:,:), molecule_gene2(:,:)
