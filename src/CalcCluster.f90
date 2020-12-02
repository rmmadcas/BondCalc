IMPLICIT NONE
INTEGER::flag
REAL (KIND=r):: natomspm

SELECT CASE (flag)
   CASE (1)
    OPEN(61,FILE="out-cluster-details1.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(61,'(A67)') "#        Movie     cluster   freepairs cluster-element cluster-atom"
   
    OPEN(71,FILE="out-cluster-summary1.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(71,'(A53)') "#           Movie       Cluster-nº nº-aggregates"

    OPEN(81,FILE="out-final-average1.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(81,'(A58)') "#   Nº molecules    % Agreggates    Avg nº of aggregated molecules"
   
    OPEN(91,FILE="All-bonds1.dat",STATUS='OLD', ACTION='READ') 
    natomspm=SUM(maxi1)/nmovies
    nn=nn1
    d=d1
    maxilinepmov=maxilinepmov1     

    nlinestotal=nn-1
    maxfreepairs=MAXVAL(maxilinepmov)
    ALLOCATE (ini(1:maxfreepairs), inj(1:maxfreepairs))
    ALLOCATE (initot(1:nlinestotal), injtot(1:nlinestotal))
    ALLOCATE (clus(1:5000))

    READ(91,*)firstline  !comentario
    DO j=1,nlinestotal
         READ(91,*) lineread, initot(j), injtot(j)    !movie, atom i, atom j
    END DO
    CLOSE (UNIT=91)

   CASE (2)
    OPEN(62,FILE="out-cluster-details2.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(62,'(A67)') "#        Movie     cluster   freepairs cluster-element cluster-atom"

    OPEN(72,FILE="out-cluster-summary2.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(72,'(A53)') "#           Movie       Cluster-nº nº-aggregates"

    OPEN(82,FILE="out-final-average2.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(82,'(A58)') "#   Nº molecules    % Agreggates    Avg nº of aggregated molecules"
    
    OPEN(92,FILE="All-bonds2.dat",STATUS='OLD', ACTION='READ')
    natomspm=SUM(maxi2)/nmovies
    nn=nn2
    d=d2
    maxilinepmov=maxilinepmov2
    nlinestotal=nn-1
    maxfreepairs=MAXVAL(maxilinepmov)
    ALLOCATE (ini(1:maxfreepairs), inj(1:maxfreepairs))
    ALLOCATE (initot(1:nlinestotal), injtot(1:nlinestotal))
    ALLOCATE (clus(1:5000))

    READ(92,*)  !comentario
    DO j=1,nlinestotal
         READ(92,*) lineread, initot(j), injtot(j)    !movie, atom i, atom j
    END DO
    CLOSE (UNIT=92)
       
   CASE (3)
    OPEN(63,FILE="out-cluster-detailsMix.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(63,'(A67)') "#        Movie     cluster   freepairs cluster-element cluster-atom"

    OPEN(73,FILE="out-cluster-summaryMix.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(73,'(A53)') "#           Movie       Cluster-nº nº-aggregates"

    OPEN(83,FILE="out-final-averageMix.dat",STATUS='NEW', ACTION='WRITE')
    WRITE(83,'(A58)') "#   Nº molecules    % Agreggates    Avg nº of aggregated molecules"

    OPEN(93,FILE="All-bondsMix.dat",STATUS='OLD', ACTION='READ')
    natomspm=(SUM(maxi1)+SUM(maxi2))/nmovies
    nn=nn3
    d=d3
    maxilinepmov=maxilinepmov3     
    nlinestotal=nn-1
    maxfreepairs=MAXVAL(maxilinepmov)
    ALLOCATE (ini(1:maxfreepairs), inj(1:maxfreepairs))
    ALLOCATE (initot(1:nlinestotal), injtot(1:nlinestotal))
    ALLOCATE (clus(1:5000))

    READ(93,*)  !comentario
    DO j=1,nlinestotal
         READ(93,*) lineread, initot(j), injtot(j)    !movie, atom i, atom j
    END DO
    CLOSE (UNIT=93)

END SELECT

   a=1      !start of movie pairs
   b=maxilinepmov(1)     !end of movie pairs
   
   naggregates=0
   DO cntmov=1,nmoviestotal  !main loop movie  
       CALL readpairs
       IF (b.LT.nlinestotal) THEN
           a=b+1
           b=b+maxilinepmov(cntmov+1)
       END IF
       freepairs=maxilinepmov(cntmov)  !parejas libres = parejas de atomos en la movie en la que estoy
       mm=1
       clus=0
       cntaggreg=0
       DO s=1,125   !main loop frame     
           CALL removeclusterfound     !elimina un cluster lleno
           IF (freepairs.EQ.0) EXIT    !si ya no hay mas parejas libres, sal del bucle
           clus(1)=ini(1); clus(2)=inj(1)      !primer y segundo elemento del cluster son la primera pareja
           nn=1; f=1; mm=2
           c=0; d=1
           DO l=1,250   !loop 1
               IF (d.EQ.mm) EXIT
               d=mm
               DO k=2,freepairs     !loop 2    
                   DO nn=1,c+mm     !loop 3
                       IF (ini(k).EQ.clus(nn)) THEN
                           DO f=1,c+mm     !loop 4.1
                               IF (inj(k).EQ.clus(f)) THEN
                                   GO TO 10
                               END IF
                           END DO       !loop 4.1
                           mm=mm+1; clus(mm)=inj(k)
                           GO TO 10
                       END IF

                       IF (inj(k).EQ.clus(nn)) THEN
                           DO f=1,c+mm     !loop 4.2
                               IF (ini(k).EQ.clus(f)) THEN
                                   GO TO 10
                               END IF
                           END DO       !loop 4.2
                           mm=mm+1; clus(mm)=ini(k)
                           GO TO 10
                       END IF
                   END DO     !loop 3
10                 CONTINUE
               END DO      !loop 2
               
           END DO       !loop 1
          
           DO i=1,c+mm 
               WRITE(60+flag,*) cntmov, s, freepairs, i, clus(i)
           END DO
               naggregates(cntmov)=naggregates(cntmov)+mm
               cntaggreg=cntaggreg+1
               WRITE(70+flag,*) "MOV", cntmov, cntaggreg, mm
       END DO   !main loop frame
       WRITE(70+flag,*) "TOTAL", cntmov, cntaggreg, naggregates(cntmov)
   END DO   !main loop movie
   avgaggregates=REAL(SUM(naggregates))/REAL(nmoviestotal)
   percaggregates=avgaggregates/(natomspm)*100

   WRITE(80+flag,'(F16.5,F16.5,F18.5)') natomspm, percaggregates, avgaggregates
   
   CLOSE (UNIT=60+flag); CLOSE (UNIT=70+flag); CLOSE (UNIT=80+flag)

DEALLOCATE (ini, inj)
DEALLOCATE (initot, injtot)
DEALLOCATE (clus)
