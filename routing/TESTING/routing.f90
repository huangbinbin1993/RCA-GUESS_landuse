!Subroutines for routing the runoff through catchments one time step
!Written by Charlotta Pers 2008-02-14
!Contains subroutines: routcatchment and routgrid
!The subroutines uses the module routmodvar
!

SUBROUTINE ROUTCATCHMENT(nsubs,acctime,paths,area,fractions,outflow1,soil,riverq,lakew,     &
     cmidrx,cmidry,cmidrw,nlat,nlon,ntile,nlat2,nlon2,runoff,GCMrunoff)
  !------------------------------------------------------------------------------------------------------------------	

  USE ROUTMODVAR, ONLY : maxtiles, i_landriver, i_lakeoutflow, dimriverq, dimcell
  IMPLICIT NONE

  !Declaration of arguments
  INTEGER, INTENT(IN) :: nsubs                      !Number of subbasins, routing grid cells
  INTEGER, INTENT(IN) :: acctime                    !Time (s) over which runoff is accumulated (e.g. 86400)
  INTEGER, INTENT(IN) :: paths(nsubs)               !Receiving grid cells
  REAL, INTENT(IN)    :: area(nsubs)                !Area of routing grid cell (km2)
  REAL, INTENT(IN)    :: fractions(nsubs,maxtiles)  !Fractions of forest,open,lakeland and lake tiles
  REAL, INTENT(OUT)   :: outflow1(nsubs)            !Discharge in rivers
  REAL, INTENT(INOUT) :: soil(nsubs)                !Water in soil shallow ground water (m3)
  REAL, INTENT(INOUT) :: lakew(nsubs)               !Water level of lake (m)
  REAL, INTENT(INOUT) :: riverq(0:dimriverq,nsubs)  !Water in river (m3)
  INTEGER,INTENT(IN)  :: cmidrx(nsubs,dimcell)      !subbasin ID for RCA/GCM-grid idx
  INTEGER,INTENT(IN)  :: cmidry(nsubs,dimcell)      !subbasin ID for RCA/GCM-grid idy
  REAL,INTENT(IN)     :: cmidrw(nsubs,dimcell)      !subbasin ID for RCA/GCM-grid weight
  INTEGER, INTENT(IN) :: nlat                       !RCA grid size
  INTEGER, INTENT(IN) :: nlon                       !RCA grid size
  INTEGER, INTENT(IN) :: ntile                      !Number of tiles for RCA runoff
  INTEGER, INTENT(IN) :: nlat2                      !GCM grid size
  INTEGER, INTENT(IN) :: nlon2                      !GCM grid size
  REAL,INTENT(IN)     :: runoff(nlon,nlat,ntile)    !Runoff from RCA-grid (mm/acctime): idx,idy,tile
  REAL,INTENT(IN)     :: gcmrunoff(nlon2,nlat2)     !Runoff from GCM-grid (mm/acctime): idx,idy

  !Declaration of local variables
  INTEGER i,isub
  INTEGER ttwhole
  REAL ttpart,transtime
  REAL lakefrac,landfrac,par(6),riverlength
  REAL irrigation
  REAL gridoutflow
  REAL dataarray(8)
  REAL accinflow(nsubs)               !Accumulated inflow to routing grid cell
  LOGICAL conduct(2)

  !Declaration of model parameters
  REAL, PARAMETER :: velocity   = 0.7          !river velocity (m/s)
  REAL, PARAMETER :: recession  = 0.1          !groundwater recession coefficient (1/d)
  REAL, PARAMETER :: ratingrate = 2.0          !rate of lake rating curve (m/s if ratingexp=2)
  REAL, PARAMETER :: ratingexp  = 2.0          !exponent in lake rating curve (-)

  ! write(*,*) ' routcatchment start '

  !Preparations
  accinflow = 0
  par(1) = recession
  par(5) = ratingrate
  par(6) = ratingexp

  ! write(*,*) ' nsubs ',nsubs

  basin_loop : &
       DO isub = 1, nsubs        !for every subcatchment

  ! write(*,*) ' isub ',isub

  !Set parameter values, subbasin specific
  riverlength = 1.5 * SQRT(area(isub)) * 1000.   !m
  transtime = riverlength / velocity / REAL(acctime)     !delay in timesteps
  ttwhole = INT(transtime)
  ttpart = transtime - REAL(ttwhole)
  par(2) = transtime
  par(3) = REAL(ttwhole)
  par(4) = ttpart

  !Get tiles area fractions
  landfrac = fractions(isub,1) + fractions(isub,2) + fractions(isub,3)
  lakefrac = fractions(isub,4)
  DO i = 1, maxtiles
     dataarray(i) = fractions(isub,i)
  ENDDO

  !Get runoff for this grid cell
  dataarray(5:7) = 0.

  ! write(*,*) ' dimcell ',dimcell

  DO i = 1, dimcell

     ! write(*,*) ' i ',i,cmidrx(isub,i)
     IF ( abs(cmidrw(isub,i)) .gt. 1.e-5) then 

        IF(cmidrx(isub,i)>0)THEN      !RCAgrid
           IF(ntile==1)THEN
              dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * runoff(cmidrx(isub,i),cmidry(isub,i),1)
           ELSE
              dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * runoff(cmidrx(isub,i),cmidry(isub,i),:)
           ENDIF
        ELSEIF(cmidrx(isub,i)<0)THEN   !GCMgrid

           ! write(*,*) ' GCM ',cmidrx(isub,i),cmidry(isub,i)

           dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * GCMrunoff(-cmidrx(isub,i),-cmidry(isub,i))  !same runoff all tiles
        ENDIF
     ENDIF      ! cmidrw > 0
  ENDDO

  !Irrigation not yet included
  irrigation = 0.
  dataarray(8) = irrigation

  !Set type of routing for the grid cell
  conduct = .FALSE.
  IF(landfrac>0) conduct(i_landriver)   = .TRUE.
  IF(lakefrac>0) conduct(i_lakeoutflow) = .TRUE.

  !     write(*,*) ' isub ',isub

  !Transport water through the grid cell
  CALL ROUTGRID(isub,paths(isub),nsubs,acctime,area(isub),dataarray,par,irrigation,accinflow,gridoutflow,soil,riverq,lakew,conduct)

  !     write(*,*) ' acctime ',acctime 

  !Save discharge
  outflow1(isub) = gridoutflow / REAL(acctime)   !outflow = m3/s, gridoutflow=m3/timestep

ENDDO basin_loop

!        write(*,*) ' routcatchment end '

END	SUBROUTINE ROUTCATCHMENT

!------------------------------------------------------------------------------------------------------------------	
SUBROUTINE ROUTGRID(i,ito,nsubs,acctime,gridarea,indata,parvalues,abstract,accinflow,outflow3,soil,riverq,lakew,conduct)
  !------------------------------------------------------------------------------------------------------------------	

USE ROUTMODVAR, ONLY : i_landriver, i_lakeoutflow, dimriverq
IMPLICIT NONE 

!Declaration of arguments
INTEGER, INTENT(IN) :: i                  				!subbasin, routing grid cell
INTEGER, INTENT(IN) :: ito                				!subbasin, recieving routing grid cell
INTEGER, INTENT(IN) :: nsubs              				!Number of subbasins, routing grid cells
INTEGER, INTENT(IN) :: acctime            				!timestep (s)
REAL, INTENT(IN)    :: gridarea           				!area of grid (km2)
REAL, INTENT(IN)    :: indata(8)          				!input data needed by subroutine
REAL, INTENT(IN)    :: parvalues(6)       				!parameter values for subbasin
REAL, INTENT(INOUT) :: abstract           				!allowed abstraction (out) (m3/timestep)
REAL, INTENT(INOUT) :: accinflow(nsubs)   				!accumulated routed flow (m3/timestep)
REAL, INTENT(OUT)   :: outflow3           				!Outflow from routing grid cell (m3/timestep)
REAL, INTENT(INOUT) :: soil(nsubs)        				!Water in soil shallow ground water (m3)
REAL, INTENT(INOUT) :: lakew(nsubs)       				!Water level of lake (m)
REAL, INTENT(INOUT) :: riverq(0:dimriverq,nsubs)  !Water in river (m3)
LOGICAL, INTENT(IN) :: conduct(2)                 !calculations to be made

!Declaration of local variables
INTEGER y,j
REAL    x,qupstream,qin,transq
REAL    ldepthm,outflowm,outflowm3s,rating1,rating2,soilrunoff
REAL    lakearea   !m2
REAL    runoff(3),tilefrac(3),lakefrac

!Preparations
tilefrac = indata(1:3)    !forest,open,lakeland
lakefrac = indata(4)      !outlet lake
runoff   = indata(5:7)

!Routing
!Here begins the routing part. The primary units are m3 and m3/timestep. The waterlevel of the lakes are measured in m.
!Shallow ground water
IF(conduct(i_landriver))THEN
  qin = 0.
  DO j = 1,3
     qin = qin + runoff(j) * gridarea * tilefrac(j) * 1000.      !runoff, m3/tidsteg, assuming runoff=mm/timestep
  ENDDO
  soil(i) = soil(i) + qin                    !m3
  soilrunoff = parvalues(1) * soil(i)        !m3/timesteg   !if timestep>1d we should have some loop here
  soil(i) = soil(i) - soilrunoff
ELSE
  soilrunoff = 0.
ENDIF

!Flow from upstream routing grid cells
qupstream = accinflow(i)   !m3/timestep

!River
IF(conduct(i_landriver))THEN
   !Inflow to river
  qin = soilrunoff + qupstream                       !!m3/timesteg

  !Translation (delay) in river
  riverq(0,i) = qin       !Add new inflow to translation variable
  IF(parvalues(2)>0)THEN
     y = NINT(parvalues(3))
     x = parvalues(4)
     transq = (1.-x)*riverq(y,i) + x*riverq(y+1,i)   !Calculate flow from river after translation
     riverq(1:y+1,i) = riverq(0:y,i)
  ELSE
     transq = qin
  ENDIF
ELSE
  transq = qupstream            !River outflow = transq       !m3/timestep
ENDIF

!Lake
!"Runoff" on lake
qin = transq                                    !m3/timestep
IF(lakefrac<0) qin = qin + runoff(3) * gridarea * (1. - SUM(tilefrac)) * 1000.      !lake runoff, m3/tidsteg

IF(conduct(i_lakeoutflow))THEN
  lakearea = lakefrac * gridarea * 1000000.   !m2

  !Inflow to lake
  lakew(i) = lakew(i) + qin / lakearea        !m

  !Outflow from lake
  ldepthm = 0.                  !Threshold olake depth (m)
  rating1 = parvalues(5)
  rating2 = parvalues(6)
  IF(rating1>0)THEN
     IF(lakew(i)>ldepthm)THEN
        outflowm3s = rating1 * (lakew(i)-ldepthm)**rating2  !General rating curve (m3/s)   !if timestep>1d we should have some loop her
     ELSE
        outflowm3s = 0.
     ENDIF
     outflowm = outflowm3s * REAL(acctime) / lakearea     !m3/s -> m/tidsteg
     IF(outflowm>lakew(i)-ldepthm)THEN     !Check
        IF(lakew(i)>ldepthm)THEN
           outflowm = lakew(i)-ldepthm
        ELSE
           outflowm = 0.
        ENDIF
     ENDIF
  ELSE
     IF(lakew(i)>ldepthm)THEN
        outflowm = lakew(i)-ldepthm
     ELSE
        outflowm = 0.
     ENDIF
  ENDIF
  lakew(i) = lakew(i) - outflowm      !m/timestep
  outflow3 = outflowm * lakearea      !m3/timestep

ELSE    
   !No lake
  outflow3 = qin
ENDIF

!Transport outflow to next grid cell
IF(ito>0)  accinflow(ito) = accinflow(ito) + outflow3

END SUBROUTINE ROUTGRID
