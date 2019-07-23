PROGRAM ROUTINGMODEL

  IMPLICIT NONE

  !Dimensions of delivered RCA runoff
  INTEGER, PARAMETER :: nlon  = 70
  INTEGER, PARAMETER :: nlat  = 78
  INTEGER, PARAMETER :: ntile  = 1		!vilken ordning?
  !Dimensions of delivered GCM runoff
  INTEGER, PARAMETER :: nlon2  = 3
  INTEGER, PARAMETER :: nlat2  = 3

  INTEGER idayno,ntime
  INTEGER iyy,imm,idd
  INTEGER acctime
  REAL RCArunoff(nlon,nlat,ntile)    !Runoff from RCA-grid (mm/acctime): id,tile	
  REAL GCMrunoff(nlon2,nlat2)        !Runoff from GCM-grid (mm/acctime): id
  INTEGER col, xid(5460), yid(5460), date   !for RCArunoff.txt
  INTEGER ndischarge                 !Number of discharge values
  REAL discharge(1000)               !Resulting discharge
  REAL soil(1000)                    !Soil moisture
  REAL river(0:12,1000)              !River water
  REAL lake(1000)                    !Lake water stage

  !Preparations
  ntime = 1948
  ntime = 14
  acctime = 86400     !s

  !Temporary ASCII version of RCArunoff input
  OPEN(10,FILE='rcmgrid_mod.txt', STATUS='old')
  READ(10,*)  !Comments
  DO
     READ(10,*,END=100)  col, xid(col), yid(col)
  ENDDO

100 CLOSE(10)
  OPEN(12,FILE='RCArunoff.txt',STATUS='old')
  READ(12,*)  !Comments
  GCMrunoff = 0.
  soil = 0.
  river = 0.
  lake = 0.

  !First call, preparations
  CALL RIVER_ROUTING(.TRUE.,.FALSE.,1984,9,1,acctime,nlon,nlat,ntile,nlon2,nlat2,1000,ndischarge,soil,river,lake,RCArunoff,GCMrunoff,discharge)

  !Calculate the discharge for every day
  time_loop: &
       DO idayno = 1,ntime

  !Temporary ASCII version of RCArunoff input
  READ(12,*) date, (RCArunoff(xid(col),yid(col),1),col=1,5460)
  iyy=date/10000
  imm=date/100 - iyy*100
  idd=date - iyy*10000 - imm*100

  !Routa RCA/GCM-runoff	
  CALL RIVER_ROUTING(.FALSE.,.FALSE.,iyy,imm,idd,acctime,nlon,nlat,ntile,nlon2,nlat2,1000,ndischarge,soil,river,lake,RCArunoff,GCMrunoff,discharge)

ENDDO time_loop

CLOSE(12)

END	PROGRAM ROUTINGMODEL
