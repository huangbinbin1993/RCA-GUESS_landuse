!Subroutine for routing the runoff of RCA/GCM 
!Written by Charlotta Pers 2008-03-18
!Calls subroutine: routcatchment
!The subroutine uses the module routmodvar
!Limitations: Only calculates on a daily timestep.
!The program is limited to 1000 outlet points and a maximum subbasin size of 100000 km2. 

SUBROUTINE RIVER_ROUTING(firstcall,savedstate,iyy,imm,idd,acctime,nlon,nlat,nlon2,nlat2, &
     dimstate,noutflow,soilstate,riverstate,lakestate,                                &
     runoff,GCMrunoff,dischargearray)

  USE IFPORT
  USE ROUTMODVAR
  IMPLICIT NONE

  !Declaration of arguments
  LOGICAL, INTENT(IN) :: firstcall                    !First call of subroutine for this run
  LOGICAL, INTENT(IN) :: savedstate                   !Use saved state
  !  INTEGER, INTENT(IN) :: idayno                       !Current timestep for routing routine (day)
  INTEGER, INTENT(IN) :: iyy                          !Current year
  INTEGER, INTENT(IN) :: imm                          !Current month
  INTEGER, INTENT(IN) :: idd                          !Current day
  !  INTEGER, INTENT(IN) :: ntime                        !Number of timesteps (days)
  INTEGER, INTENT(IN) :: acctime                      !Time (s) over which runoff is accumulated (i.e. 86400)
  INTEGER, INTENT(IN) :: nlon                         !RCA grid size
  INTEGER, INTENT(IN) :: nlat                         !RCA grid size
!  INTEGER, INTENT(IN) :: ntile                        !Number of tiles, RCA grid size
  INTEGER, INTENT(IN) :: nlon2                        !GCM grid size
  INTEGER, INTENT(IN) :: nlat2                        !GCM grid size
  INTEGER, INTENT(IN) :: dimstate                     !Size of state variables
  INTEGER,INTENT(OUT) :: noutflow                     !Number of discharge points
  REAL, INTENT(INOUT) :: soilstate(dimstate)          !soil response box (uz), one per subbasin (mm)
  REAL, INTENT(INOUT) :: riverstate(0:12,dimstate)    !river water in translation delay (m3/timestep)
  REAL, INTENT(INOUT) :: lakestate(dimstate)          !water stage of lake (m)
  REAL,INTENT(IN)     :: runoff(nlon,nlat,ntile)      !Runoff from RCA-grid (mm/acctime): idx,idy,tile
  REAL,INTENT(IN)     :: gcmrunoff(nlon2,nlat2)       !Runoff from GCM-grid (mm/acctime): idx,idy
  REAL,INTENT(OUT)    :: dischargearray(1000)         !Discharge at discharge points (m3/s)

  !Declaration of local variables
  INTEGER ID,path1,xlat,ylon
  INTEGER xrid1,xrid2,xrid3,xrid4,xrid5,xrid6
  INTEGER yrid1,yrid2,yrid3,yrid4,yrid5,yrid6
  INTEGER I,J,formatcode,mcols
  INTEGER nsubriver,maxriverrow
  INTEGER ilen,istat
  LOGICAL istat2
  REAL    wgt1,wgt2,wgt3,wgt4,wgt5,wgt6
  REAL    x,flowacc
  REAL    subbasinarea,FRAC1,FRAC2,FRAC3,FRAC4,SEAFRAC
  CHARACTER(LEN=20) text
  CHARACTER(LEN=200) filedir,filename
  CHARACTER(LEN=1000) line
  INTEGER, ALLOCATABLE :: resultid(:)
  INTEGER, ALLOCATABLE :: resultrows(:,:)
  INTEGER, ALLOCATABLE :: riverrows(:,:)
  REAL, ALLOCATABLE :: grid(:,:)
  REAL, ALLOCATABLE :: outflow(:)
  REAL, ALLOCATABLE :: resoutflow(:)

  !------------------------------------------------
  ! FIRST TIMESTEP
  !------------------------------------------------
  !Preparation, should only be done first call of river_routing

  ! write(*,*) ' river_routing start'

  IF(firstcall)THEN	
     !Get result point data
     CALL GETENVAR('DEF_ROUT1',filedir,ilen,istat)  !LINUX-RCA
     IF(istat/=0) RETURN                            !LINUX-RCA
     !    istat2 = SETENVQQ('DEF_ROUT1=..\routing\')    !WINDOWS-TEST
     !    ilen = GETENVQQ('DEF_ROUT1',filedir)          !WINDOWS-TEST
     !     filedir=''
     !     ilen=0
     !    WRITE(*,*) filedir
     OPEN(300,FILE=filedir(1:ilen)//'Point.txt',STATUS='old')
     READ(300,*)  !Comments
     READ(300,*)  nresult, formatcode
     ALLOCATE(resultrows(nresult,2))
     ALLOCATE(riverrows(nresult,2))
     ALLOCATE(resultid(nresult))
     resultid = 0
     nriver = 0
     maxriverrow = 0
     DO i = 1, nresult
        IF(formatcode==1)THEN       !lat lon in decimal degrees
           READ(300,*,END=1) id,text, x,x,x,x,x, resultrows(i,1),resultrows(i,2)
        ELSEIF(formatcode==2)THEN   !lat lon in degrees, min, sec
           READ(300,*,END=1) id,text, x,x,x,x,x,x,x,x,x, resultrows(i,1),resultrows(i,2)
        ELSE
           !error
        ENDIF
        IF(maxriverrow<resultrows(i,2))THEN
           nriver = nriver + 1
           maxriverrow = resultrows(i,2)
           riverrows(nriver,:) = resultrows(i,:)
        ENDIF
        resultid(i)=id
     ENDDO
1    CLOSE(300)
     nsub = maxriverrow

     !     write(*,*) ' river_routing Point file ok nsub= ',nsub

     !Ger river row data
     ALLOCATE(iriverrow(nriver,2))
     iriverrow = riverrows(1:nriver,:)
     DEALLOCATE(riverrows)

     !Get fysiographic data
     OPEN(300, FILE = filedir(1:ilen)//'RoutingData.txt',STATUS='old')  !"GeoData"
     READ(300,602) line                                !Read headings to calculate number of columns
     j=0
     DO i = 1,1000
        IF(line(i:i)==CHAR(32).OR.line(i:i)==CHAR(9).OR.line(i:i)==CHAR(13))THEN   !space, tab or carriage return
        ELSE
           IF(line(i+1:i+1)==CHAR(32).OR.line(i+1:i+1)==CHAR(9).OR.line(i+1:i+1)==CHAR(13)) j=j+1
        ENDIF
     ENDDO
     mcols=j     !number of columns in file

     !    write(*,*) ' river_routing j= ',j

     ALLOCATE(grid(nsub,28))
     xrid2=0;xrid3=0;xrid4=0;xrid5=0;xrid6=0;yrid2=0;yrid3=0;yrid4=0;yrid5=0;yrid6=0;wgt2=0.;wgt3=0.;wgt4=0.;wgt5=0.;wgt6=0.
     I = 0
     DO
        IF(mcols==34)THEN
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,xrid5,yrid5,xrid6,yrid6,wgt1,wgt2,wgt3,wgt4,wgt5,wgt6
           dimcell = 6
        ELSEIF(mcols==31)THEN
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,xrid5,yrid5,wgt1,wgt2,wgt3,wgt4,wgt5
           dimcell = 5
        ELSEIF(mcols==28)THEN
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,wgt1,wgt2,wgt3,wgt4
           dimcell = 4
        ELSEIF(mcols==25)THEN
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,wgt1,wgt2,wgt3
           dimcell = 3
        ELSEIF(mcols==22)THEN
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,xrid2,yrid2,wgt1,wgt2
           dimcell = 2
        ELSE
           READ(300,*, END = 2)ID,path1,x,flowacc,subbasinarea,x,x,FRAC2,FRAC1,SEAFRAC,FRAC3,FRAC4,x,x,xlat,ylon,xrid1,yrid1,wgt1
           dimcell = 1
        ENDIF
        I = I + 1
        grid(I,1) = ID              !ID NUMBER OF GRID CELL
        grid(I,2) = SEAFRAC         !FRACTION OF SEA ON LAND (SHOULD BE ZERO)
        grid(I,4) = path1           !ID NUMBER OF RECIEVING GRID CELL
        grid(I,5) = FRAC1           !FRACTION OF FOREST LAND
        grid(I,6) = FRAC2           !FRACTION OF OPEN LAND
        grid(I,7) = FRAC3          	!FRACTION OF LAKE	LAND
        grid(I,8) = FRAC4          	!FRACTION OF OUTLETLAKE	
        grid(I,9) = subbasinarea   	!SUBBASIN AREA (KM2)	
        grid(I,10) = flowacc      	!CATCHMENT AREA (KM2)	
        grid(I,11) = xrid1          !XID FOR RUNOFF GCM/RCM CELL 1
        grid(I,12) = xrid2          !XID FOR RUNOFF GCM/RCM CELL 2
        grid(I,13) = xrid3          !XID FOR RUNOFF GCM/RCM CELL 3
        grid(I,14) = xrid4          !XID FOR RUNOFF GCM/RCM CELL 4
        grid(I,15) = xrid5          !XID FOR RUNOFF GCM/RCM CELL 5
        grid(I,16) = xrid6          !XID FOR RUNOFF GCM/RCM CELL 6
        grid(I,17) = yrid1          !YID FOR RUNOFF GCM/RCM CELL 1
        grid(I,18) = yrid2          !YID FOR RUNOFF GCM/RCM CELL 2
        grid(I,19) = yrid3          !YID FOR RUNOFF GCM/RCM CELL 3
        grid(I,20) = yrid4          !YID FOR RUNOFF GCM/RCM CELL 4
        grid(I,21) = yrid5          !YID FOR RUNOFF GCM/RCM CELL 5
        grid(I,22) = yrid6          !YID FOR RUNOFF GCM/RCM CELL 6
        grid(I,23) = wgt1           !WEIGHT FOR RUNOFF GCM/RCM CELL 1
        grid(I,24) = wgt2           !WEIGHT FOR RUNOFF GCM/RCM CELL 2
        grid(I,25) = wgt3           !WEIGHT FOR RUNOFF GCM/RCM CELL 3
        grid(I,26) = wgt4           !WEIGHT FOR RUNOFF GCM/RCM CELL 4
        grid(I,27) = wgt5           !WEIGHT FOR RUNOFF GCM/RCM CELL 5
        grid(I,28) = wgt6           !WEIGHT FOR RUNOFF GCM/RCM CELL 6

        if ( i.eq.nsub) goto 2

     ENDDO
2    CLOSE(300)
     IF(nsub /= I)THEN
        WRITE(*,*) 'ERROR: Number of subbasins not the same as number of rows in file'
     ENDIF


     !   write(*,*) ' river_routing input ok'

     !Prepare path variable, save to module variable patharray
     IF(.NOT.ALLOCATED(patharray)) ALLOCATE(patharray(nsub))
     patharray = 0.
     DO I = 1,nsub
        DO J = I+1,nsub
           IF(grid(I,4)==grid(J,1))THEN
              patharray(I)=J
              EXIT
           ENDIF
        ENDDO
     ENDDO

     !Prepare area and tilefraction variables
     IF(.NOT.ALLOCATED(tilefractions)) ALLOCATE(tilefractions(nsub,maxtiles))
     IF(.NOT.ALLOCATED(basinarea)) ALLOCATE(basinarea(nsub))

     ! write(*,*) ' river_routing alloc ok'

     WHERE(grid(:,2)>0)
        tilefractions(:,1)=grid(:,5)+grid(:,2)*grid(:,5)/(grid(:,5)+grid(:,6))    !safe for positive sea fraction
        tilefractions(:,2)=grid(:,6)+grid(:,2)*grid(:,6)/(grid(:,5)+grid(:,6))
     ELSEWHERE
        tilefractions(:,1)=grid(:,5)
        tilefractions(:,2)=grid(:,6)
     ENDWHERE
     tilefractions(:,3:4)=grid(:,7:8)
     !    tilefractions(:,:)=grid(:,5:8)
     basinarea = grid(:,9)

     !Prepare state variables, allocate and initiate to zero
     dimriverq = INT(SQRT(MAXVAL(grid(:,9)))*2150./acctime)+3     !approximation of river delay based on gridsize (1.5/0.7*1000)
     IF(.NOT.ALLOCATED(soilvec)) ALLOCATE(soilvec(nsub))
     IF(.NOT.ALLOCATED(lakewvec)) ALLOCATE(lakewvec(nsub))
     IF(.NOT.ALLOCATED(riverqvec)) ALLOCATE(riverqvec(0:dimriverq,nsub))
     IF(savedstate)THEN
        soilvec = soilstate(1:nsub)
        lakewvec = lakestate(1:nsub)
        riverqvec = riverstate(0:dimriverq,1:nsub)
     ELSE
        soilvec = 0.
        lakewvec = SQRT(meanrunoff * grid(:,10) / (86.400 * 365. * 2.))
        riverqvec = 0.
     ENDIF

     !Prepare output variables, save river result point id to module variable
     IF(.NOT.ALLOCATED(iresoutflow)) ALLOCATE(iresoutflow(nresult))
     iresoutflow(:) = resultrows(1:nresult,2)
     noutflow = nresult
     dischargearray = 0.

     !Prepare for reading runoff. Save id for RGCM grid cells.
     IF(.NOT.ALLOCATED(cmidx)) ALLOCATE(cmidx(nsub,dimcell))
     IF(.NOT.ALLOCATED(cmidy)) ALLOCATE(cmidy(nsub,dimcell))
     IF(.NOT.ALLOCATED(cmwgt)) ALLOCATE(cmwgt(nsub,dimcell))
     cmidx = NINT(grid(:,11:10+dimcell))
     cmidy = NINT(grid(:,17:16+dimcell))
     cmwgt = grid(:,23:22+dimcell)

     !Prepare writing result to file
     filename=filedir(1:ilen)//'allresult_'
     WRITE(filename(ilen+11:ilen+14),'(I4.4)') iyy
     WRITE(filename(ilen+15:ilen+16),'(I2.2)') imm
     WRITE(filename(ilen+17:ilen+18),'(I2.2)') idd
     OPEN(300,FILE=filename(1:ilen+18))
     filename=filedir(1:ilen)//'pointresult_'
     WRITE(filename(ilen+13:ilen+16),'(I4.4)') iyy
     WRITE(filename(ilen+17:ilen+18),'(I2.2)') imm
     WRITE(filename(ilen+19:ilen+20),'(I2.2)') idd
     OPEN(301,FILE=filename(1:ilen+20))
     WRITE(300,601)'yyyy mm dd',(i,i=1,maxriverrow)
     WRITE(301,601)'yyyy mm dd',(resultid(i),i=1,nresult)

     DEALLOCATE(grid)
     DEALLOCATE(resultid)
     DEALLOCATE(resultrows)

  ENDIF


  !Prepare for calculating routing current timestep
  ALLOCATE(outflow(nsub))        !Allocate memory for discharge results
  ALLOCATE(resoutflow(nresult))

  outflow=0

  ! write(*,*) ' start routecatchment   nriver ',nriver

  !------------------------------------------------
  ! LOOP OVER CATCHMENTS
  !------------------------------------------------
  !Calculate the routing; push the runoff of the current timestep through the river
  DO i = 1,nriver
     nsubriver = iriverrow(i,2) - (iriverrow(i,1) - 1)

     ! write(*,*) ' routecatchment  i ',i

     CALL ROUTCATCHMENT(nsubriver,acctime,patharray(iriverrow(i,1):iriverrow(i,2))-(iriverrow(i,1)-1), &
          basinarea(iriverrow(i,1):iriverrow(i,2)),                    &
          tilefractions(iriverrow(i,1):iriverrow(i,2),:),              &
          outflow(iriverrow(i,1):iriverrow(i,2)),soilvec(iriverrow(i,1):iriverrow(i,2)),                  &
          riverqvec(:,iriverrow(i,1):iriverrow(i,2)),lakewvec(iriverrow(i,1):iriverrow(i,2)),             &
          cmidx(iriverrow(i,1):iriverrow(i,2),:),cmidy(iriverrow(i,1):iriverrow(i,2),:),cmwgt(iriverrow(i,1):iriverrow(i,2),:),   &
          nlat,nlon,ntile,nlat2,nlon2,runoff,GCMrunoff)
  ENDDO


  ! write(*,*) ' start routecatchment end '

  !------------------------------------------------
  !Finishing off the routing calculation for the current timestep
  DO i = 1, nresult           !Find the outflow at result points
     DO j = 1, nsub
        IF(j==iresoutflow(i)) resoutflow(i) = outflow(j)
     ENDDO
  ENDDO
  noutflow = nresult          !Set result variables of subroutine
  dischargearray(1:nresult) = resoutflow

  WRITE(300,600) iyy,imm,idd,outflow  !Write the results to file
  WRITE(301,600) iyy,imm,idd,resoutflow

  DEALLOCATE(outflow)         !Deallocate local variables
  DEALLOCATE(resoutflow)

  !Save states
  soilstate(1:nsub) = soilvec
  lakestate(1:nsub) = lakewvec
  riverstate(0:dimriverq,1:nsub) = riverqvec 

  !Format statements
600 FORMAT(I4,x,I2,x,I2,x,1000000(1X,ES9.2))
601 FORMAT(A10,x,1000000(1X,I9))
602 FORMAT(A1000)

  ! write(*,*) ' river_routing end '

END	SUBROUTINE RIVER_ROUTING
