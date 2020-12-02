        h=1
        DO j=1,freepairs
                DO i=1,mm
                        IF ((clus(i).EQ.ini(j)).OR.(clus(i).EQ.inj(j))) THEN
                                GO TO 15
                        END IF
                END DO
                ini(h)=ini(j)
                inj(h)=inj(j)           
                h=h+1
15              CONTINUE
        END DO
        freepairs=h-1
