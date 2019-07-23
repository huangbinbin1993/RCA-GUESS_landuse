module routmodvar
  !module for declaration of routing model variables
  use decomp,only:realkind
  implicit none
  private
  !parameter declarations
  integer, parameter :: maxtiles       = 4 !dimension of routing model tile fraction
  integer, parameter :: i_landriver    = 1 !conduct index parameter for land and river translation
  integer, parameter :: i_lakeoutflow  = 2 !conduct index parameter for lake outflow
  real(kind=realkind),    parameter :: meanrunoff = 300._realkind
  !mean runoff global (mm/year)

  !variable declarations
  integer,save:: nsub                !number of subbasins for all rivers in current run
  integer,save:: nriver              !number of rivers in current run
  integer,save:: nresult             !number of subbasins for result in current run
  integer,save:: dimriverq           !dimension of riverqvec, days delay in river
  integer,save:: dimcell             !dimension of cm-arrays, number of rgcm grid cells a subbasin get runoff from
  integer,save, allocatable :: patharray(:)        !river paths; id of recieving routing grid cells
  integer,save, allocatable :: iresoutflow(:)    !river result points; id-subbasinnumber-row of river result points
  integer,save, allocatable :: iriverrow(:,:)       !river; begin and end row in list
  real(kind=realkind),save, allocatable :: tilefractions(:,:)     !tile fractions of routing grid cells; forest, open, lake
  real(kind=realkind),save, allocatable :: basinarea(:)           !subbasin area (km2)
  integer,save, allocatable :: cmidx(:,:)           !xid for gcm/rcm grid cells corresponding to each routing grid cell
  integer,save, allocatable :: cmidy(:,:)           !yid for gcm/rcm grid cells corresponding to each routing grid cell
  real(kind=realkind),save, allocatable :: cmwgt(:,:)              !weight for gcm/rcm grid cells corresponding to each routing grid cell


  !state variable declarations
  real(kind=realkind),save, allocatable ::  soilvec(:)          !soil response box (uz), one per subbasin (mm)
  real(kind=realkind),save, allocatable ::  riverqvec(:,:)      !river water in translation delay (m3/timestep)
  real(kind=realkind),save, allocatable ::  lakewvec(:)         !water stage of lake (m)
  integer, public, save:: ntile=1

  integer, public, save::klon_bound
  integer, public, save::klat_bound
  real(kind=realkind),save::acctime=86400.0_realkind

  integer,public,save::dimstate=1000
  real(kind=realkind),allocatable,public,save:: soilstate(:) !(dimstate) soil response box (uz), one per subbasin (mm)
  real(kind=realkind),allocatable,public,save:: riverstate(:,:)!(0:12,dimstate)!river water in translation delay (m3/timestep)
  real(kind=realkind),allocatable,public,save:: lakestate(:)!(dimstate)   !water stage of lake (m)
  integer,public,save :: noutflow                     !number of discharge points
  real(kind=realkind),allocatable,public,save    :: discharge(:)!1000)  !discharge at discharge points (m3/s)

  character(len=200),save :: filedir
  real(kind=realkind), public,allocatable,save::runoff_bound(:,:)!(klon_bound,klat_bound)

  public river_routing,initiate_routing
contains
  
  subroutine initiate_routing(iyy,imm,idd,savedstate)
    implicit none
    integer, intent(in) :: iyy                          !current year
    integer, intent(in) :: imm                          !current month
    integer, intent(in) :: idd                          !current day
    logical, intent(in) :: savedstate                   !use saved state
    integer::ilen
    integer::formatcode
    integer, allocatable :: resultrows(:,:)
    integer, allocatable :: riverrows(:,:)
    integer, allocatable :: resultid(:)
    real(kind=realkind), allocatable :: grid(:,:)
    integer::maxriverrow
    integer::i,j,mcols
    integer::id
    character(len=20):: text
    character(len=1000):: line
    real(kind=realkind)::x
    integer:: xrid1,xrid2,xrid3,xrid4,xrid5,xrid6
    integer:: yrid1,yrid2,yrid3,yrid4,yrid5,yrid6
    real(kind=realkind)::    wgt1,wgt2,wgt3,wgt4,wgt5,wgt6
    integer:: path1,xlat,ylon
    real(kind=realkind)::    flowacc
    real(kind=realkind)::    subbasinarea,frac1,frac2,frac3,frac4,seafrac
    character(len=200)::filename




    !get result point data
    call read_namrouting(filedir,ilen)
    allocate(runoff_bound(klon_bound,klat_bound))
    allocate(soilstate(dimstate))
    allocate(riverstate(0:12,dimstate))
    allocate(lakestate(dimstate))
    allocate(discharge(dimstate))


    open(300,file=filedir(1:ilen)//'point.txt',status='old')
    read(300,*)  !comments
    read(300,*)  nresult, formatcode
    allocate(resultrows(nresult,2))
    allocate(riverrows(nresult,2))
    allocate(resultid(nresult))
    resultid = 0
    nriver = 0
    maxriverrow = 0
    do i = 1, nresult
       if(formatcode==1)then       !lat lon in decimal degrees
          read(300,*,end=1) id,text, x,x,x,x,x, resultrows(i,1),resultrows(i,2)
       elseif(formatcode==2)then   !lat lon in degrees, min, sec
          read(300,*,end=1) id,text, x,x,x,x,x,x,x,x,x, resultrows(i,1),resultrows(i,2)
       else
          !error
       endif
       if(maxriverrow<resultrows(i,2))then
          nriver = nriver + 1
          maxriverrow = resultrows(i,2)
          riverrows(nriver,:) = resultrows(i,:)
       endif
       resultid(i)=id
    enddo
1   close(300)
    nsub = maxriverrow

    !ger river row data
    allocate(iriverrow(nriver,2))
    iriverrow = riverrows(1:nriver,:)
    deallocate(riverrows)

    !get fysiographic data
    open(300, file = filedir(1:ilen)//'routingdata.txt',status='old')  !"geodata"
    read(300,602) line                                !read headings to calculate number of columns
    j=0
    do i = 1,1000
       if(line(i:i)==char(32).or.line(i:i)==char(9).or.line(i:i)==char(13))then   !space, tab or carriage return
       else
          if(line(i+1:i+1)==char(32).or.line(i+1:i+1)==char(9).or.line(i+1:i+1)==char(13)) j=j+1
       endif
    enddo
    mcols=j     !number of columns in file

    allocate(grid(nsub,28))
    xrid2=0;xrid3=0;xrid4=0;xrid5=0;xrid6=0;yrid2=0;yrid3=0;yrid4=0;yrid5=0;yrid6=0
    wgt2=0._realkind;wgt3=0._realkind;wgt4=0._realkind;wgt5=0._realkind;wgt6=0._realkind
    i = 0
    do
       if(mcols==34)then
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,xrid5,yrid5,xrid6,yrid6,wgt1,wgt2,wgt3,wgt4,wgt5,wgt6
          dimcell = 6
       elseif(mcols==31)then
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,xrid5,yrid5,wgt1,wgt2,wgt3,wgt4,wgt5
          dimcell = 5
       elseif(mcols==28)then
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,xrid4,yrid4,wgt1,wgt2,wgt3,wgt4
          dimcell = 4
       elseif(mcols==25)then
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,xrid2,yrid2,xrid3,yrid3,wgt1,wgt2,wgt3
          dimcell = 3
       elseif(mcols==22)then
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,xrid2,yrid2,wgt1,wgt2
          dimcell = 2
       else
          read(300,*, end = 2)id,path1,x,flowacc,subbasinarea,x,x,frac2,frac1,seafrac,frac3,frac4,x,x,xlat,&
               ylon,xrid1,yrid1,wgt1
          dimcell = 1
       endif
       i = i + 1
       grid(i,1) = real(id,realkind)        !id number of grid cell
       grid(i,2) = seafrac         !fraction of sea on land (should be zero)
       grid(i,4) = real(path1,realkind)           !id number of recieving grid cell
       grid(i,5) = frac1           !fraction of forest land
       grid(i,6) = frac2           !fraction of open land
       grid(i,7) = frac3          	!fraction of lake	land
       grid(i,8) = frac4          	!fraction of outletlake	
       grid(i,9) = subbasinarea   	!subbasin area (km2)	
       grid(i,10) = flowacc      	!catchment area (km2)	
       grid(i,11) = real(xrid1,realkind)          !xid for runoff gcm/rcm cell 1
       grid(i,12) = real(xrid2,realkind)          !xid for runoff gcm/rcm cell 2
       grid(i,13) = real(xrid3,realkind)          !xid for runoff gcm/rcm cell 3
       grid(i,14) = real(xrid4,realkind)          !xid for runoff gcm/rcm cell 4
       grid(i,15) = real(xrid5,realkind)          !xid for runoff gcm/rcm cell 5
       grid(i,16) = real(xrid6,realkind)          !xid for runoff gcm/rcm cell 6
       grid(i,17) = real(yrid1,realkind)          !yid for runoff gcm/rcm cell 1
       grid(i,18) = real(yrid2,realkind)          !yid for runoff gcm/rcm cell 2
       grid(i,19) = real(yrid3,realkind)          !yid for runoff gcm/rcm cell 3
       grid(i,20) = real(yrid4,realkind)          !yid for runoff gcm/rcm cell 4
       grid(i,21) = real(yrid5,realkind)          !yid for runoff gcm/rcm cell 5
       grid(i,22) = real(yrid6,realkind)          !yid for runoff gcm/rcm cell 6
       grid(i,23) = wgt1           !weight for runoff gcm/rcm cell 1
       grid(i,24) = wgt2           !weight for runoff gcm/rcm cell 2
       grid(i,25) = wgt3           !weight for runoff gcm/rcm cell 3
       grid(i,26) = wgt4           !weight for runoff gcm/rcm cell 4
       grid(i,27) = wgt5           !weight for runoff gcm/rcm cell 5
       grid(i,28) = wgt6           !weight for runoff gcm/rcm cell 6

       if ( i==nsub) goto 2

    enddo
2   close(300)
    if(nsub /= i)then
       write(*,*) 'error: number of subbasins not the same as number of rows in file'
    endif


    !prepare path variable, save to module variable patharray
    if(.not.allocated(patharray)) allocate(patharray(nsub))
    patharray = 0
    do i = 1,nsub
       do j = i+1,nsub
          if(grid(i,4)==grid(j,1))then
             patharray(i)=j
             exit
          endif
       enddo
    enddo

    !prepare area and tilefraction variables
    if(.not.allocated(tilefractions)) allocate(tilefractions(nsub,maxtiles))
    if(.not.allocated(basinarea)) allocate(basinarea(nsub))

    where(grid(:,2)>0.0_realkind)
       tilefractions(:,1)=grid(:,5)+grid(:,2)*grid(:,5)/(grid(:,5)+grid(:,6))    !safe for positive sea fraction
       tilefractions(:,2)=grid(:,6)+grid(:,2)*grid(:,6)/(grid(:,5)+grid(:,6))
    elsewhere
       tilefractions(:,1)=grid(:,5)
       tilefractions(:,2)=grid(:,6)
    endwhere
    tilefractions(:,3:4)=grid(:,7:8)
    !    tilefractions(:,:)=grid(:,5:8)
    basinarea = grid(:,9)

    !prepare state variables, allocate and initiate to zero
    dimriverq = nint(sqrt(maxval(grid(:,9)))*2150._realkind/acctime)+3     !approximation of river delay based on gridsize (1.5/0.7*1000)
    if(.not.allocated(soilvec)) allocate(soilvec(nsub))
    if(.not.allocated(lakewvec)) allocate(lakewvec(nsub))
    if(.not.allocated(riverqvec)) allocate(riverqvec(0:dimriverq,nsub))
    if(savedstate)then
       soilvec = soilstate(1:nsub)
       lakewvec = lakestate(1:nsub)
       riverqvec = riverstate(0:dimriverq,1:nsub)
    else
       soilvec = 0._realkind
       lakewvec = sqrt(meanrunoff * grid(:,10) / (86.400_realkind * 365._realkind * 2._realkind))
       riverqvec = 0._realkind
    endif

    !prepare output variables, save river result point id to module variable
    if(.not.allocated(iresoutflow)) allocate(iresoutflow(nresult))
    iresoutflow(:) = resultrows(1:nresult,2)
    noutflow = nresult
    discharge = 0._realkind

    !prepare for reading runoff. save id for rgcm grid cells.
    if(.not.allocated(cmidx)) allocate(cmidx(nsub,dimcell))
    if(.not.allocated(cmidy)) allocate(cmidy(nsub,dimcell))
    if(.not.allocated(cmwgt)) allocate(cmwgt(nsub,dimcell))
    cmidx = nint(grid(:,11:10+dimcell))
    cmidy = nint(grid(:,17:16+dimcell))
    cmwgt = grid(:,23:22+dimcell)

    !prepare writing result to file
    filename=filedir(1:ilen)//'allresult_'
    write(filename(ilen+11:ilen+14),'(i4.4)') iyy
    write(filename(ilen+15:ilen+16),'(i2.2)') imm
    write(filename(ilen+17:ilen+18),'(i2.2)') idd
    open(300,file=filename(1:ilen+18))
    filename=filedir(1:ilen)//'pointresult_'
    write(filename(ilen+13:ilen+16),'(i4.4)') iyy
    write(filename(ilen+17:ilen+18),'(i2.2)') imm
    write(filename(ilen+19:ilen+20),'(i2.2)') idd
    open(301,file=filename(1:ilen+20))
    write(300,601)'yyyy mm dd',(i,i=1,maxriverrow)
    write(301,601)'yyyy mm dd',(resultid(i),i=1,nresult)

    deallocate(grid)
    deallocate(resultid)
    deallocate(resultrows)
601 format(a10,x,1000000(1x,i9))
602 format(a1000)

  end subroutine initiate_routing



  !subroutine for routing the runoff of rca/gcm 
  !written by charlotta pers 2008-03-18
  !calls subroutine: routcatchment
  !the subroutine uses the module routmodvar
  !limitations: only calculates on a daily timestep.
  !the program is limited to 1000 outlet points and a maximum subbasin size of 100000 km2. 

  subroutine river_routing(iyy,imm,idd,runoff,gcmrunoff)
    use decomp
    implicit none

    !declaration of arguments
    integer, intent(in) :: iyy                          !current year
    integer, intent(in) :: imm                          !current month
    integer, intent(in) :: idd                          !current day

    real(kind=realkind),intent(in)::runoff(klon_global,klat_global,ntile)!runoff from rca-grid (mm/acctime): idx,idy,tile
    real(kind=realkind),intent(in)::gcmrunoff(klon_bound,klat_bound)!runoff from gcm-grid (mm/acctime): idx,idy

    !declaration of local variables
    integer i,j
    integer nsubriver
    integer ilen,istat
    logical istat2


    real(kind=realkind), allocatable :: outflow(:)
    real(kind=realkind), allocatable :: resoutflow(:)

    allocate(outflow(nsub))        !allocate memory for discharge results
    allocate(resoutflow(nresult))

    outflow=0.0_realkind

    ! loop over catchments
    !calculate the routing; push the runoff of the current timestep through the river
    do i = 1,nriver
       nsubriver = iriverrow(i,2) - (iriverrow(i,1) - 1)
       call routcatchment(nsubriver,patharray(iriverrow(i,1):iriverrow(i,2))-(iriverrow(i,1)-1), &
            basinarea(iriverrow(i,1):iriverrow(i,2)),                    &
            tilefractions(iriverrow(i,1):iriverrow(i,2),:),              &
            outflow(iriverrow(i,1):iriverrow(i,2)),soilvec(iriverrow(i,1):iriverrow(i,2)),                  &
            riverqvec(:,iriverrow(i,1):iriverrow(i,2)),lakewvec(iriverrow(i,1):iriverrow(i,2)),             &
            cmidx(iriverrow(i,1):iriverrow(i,2),:),cmidy(iriverrow(i,1):iriverrow(i,2),:),&
            cmwgt(iriverrow(i,1):iriverrow(i,2),:),   &
            runoff,gcmrunoff)
    enddo

    !finishing off the routing calculation for the current timestep
    do i = 1, nresult           !find the outflow at result points
       do j = 1, nsub
          if(j==iresoutflow(i)) resoutflow(i) = outflow(j)
       enddo
    enddo
    noutflow = nresult          !set result variables of subroutine
    discharge(1:nresult) = resoutflow

    write(300,600) iyy,imm,idd,outflow  !write the results to file
    write(301,600) iyy,imm,idd,resoutflow

    deallocate(outflow)         !deallocate local variables
    deallocate(resoutflow)

    !save states
    soilstate(1:nsub) = soilvec
    lakestate(1:nsub) = lakewvec
    riverstate(0:dimriverq,1:nsub) = riverqvec 

    !format statements
600 format(i4,x,i2,x,i2,x,1000000(1x,es9.2))


  end subroutine river_routing



  !subroutines for routing the runoff through catchments one time step
  !written by charlotta pers 2008-02-14
  !contains subroutines: routcatchment and routgrid
  !the subroutines uses the module routmodvar
  !

  subroutine routcatchment(nsubs, &
       paths,area,fractions,outflow1,soil,riverq,lakew,     &
       cmidrx,cmidry,cmidrw,runoff,gcmrunoff)

    use decomp
    implicit none

    !declaration of arguments
    integer, intent(in) :: nsubs                      !number of subbasins, routing grid cells
    integer, intent(in) :: paths(nsubs)               !receiving grid cells
    real(kind=realkind), intent(in)    :: area(nsubs)                !area of routing grid cell (km2)
    real(kind=realkind), intent(in)    :: fractions(nsubs,maxtiles)  !fractions of forest,open,lakeland and lake tiles
    real(kind=realkind), intent(out)   :: outflow1(nsubs)            !discharge in rivers
    real(kind=realkind), intent(inout) :: soil(nsubs)                !water in soil shallow ground water (m3)
    real(kind=realkind), intent(inout) :: lakew(nsubs)               !water level of lake (m)
    real(kind=realkind), intent(inout) :: riverq(0:dimriverq,nsubs)  !water in river (m3)
    integer,intent(in)  :: cmidrx(nsubs,dimcell)      !subbasin id for rca/gcm-grid idx
    integer,intent(in)  :: cmidry(nsubs,dimcell)      !subbasin id for rca/gcm-grid idy
    real(kind=realkind),intent(in)     :: cmidrw(nsubs,dimcell)      !subbasin id for rca/gcm-grid weight
    real(kind=realkind),intent(in)     :: runoff(klon_global,klat_global,ntile)    !runoff from rca-grid (mm/acctime): idx,idy,tile
    real(kind=realkind),intent(in)     :: gcmrunoff(klon_bound,klat_bound)     !runoff from gcm-grid (mm/acctime): idx,idy

    !declaration of local variables
    integer i,isub
    integer ttwhole
    real(kind=realkind) ttpart,transtime
    real(kind=realkind) lakefrac,landfrac,par(6),riverlength
    real(kind=realkind) irrigation
    real(kind=realkind) gridoutflow
    real(kind=realkind) dataarray(8)
    real(kind=realkind) accinflow(nsubs)               !accumulated inflow to routing grid cell
    logical conduct(2)

    !declaration of model parameters
    real(kind=realkind), parameter :: velocity   = 0.7_realkind          !river velocity (m/s)
    real(kind=realkind), parameter :: recession  = 0.1_realkind          !groundwater recession coefficient (1/d)
    real(kind=realkind), parameter :: ratingrate = 2.0_realkind          !rate of lake rating curve (m/s if ratingexp=2)
    real(kind=realkind), parameter :: ratingexp  = 2.0_realkind          !exponent in lake rating curve (-)

    ! write(*,*) ' routcatchment start '

    !preparations
    accinflow = 0.0_realkind
    par(1) = recession
    par(5) = ratingrate
    par(6) = ratingexp

    ! write(*,*) ' nsubs ',nsubs

    basin_loop : &
         do isub = 1, nsubs        !for every subcatchment

    ! write(*,*) ' isub ',isub

    !set parameter values, subbasin specific
    riverlength = 1.5_realkind * sqrt(area(isub)) * 1000._realkind   !m
    transtime = riverlength / velocity / real(acctime,realkind)     !delay in timesteps
    ttwhole = nint(transtime)
    ttpart = transtime - real(ttwhole,realkind)
    par(2) = transtime
    par(3) = real(ttwhole,realkind)
    par(4) = ttpart

    !get tiles area fractions
    landfrac = fractions(isub,1) + fractions(isub,2) + fractions(isub,3)
    lakefrac = fractions(isub,4)
    do i = 1, maxtiles
       dataarray(i) = fractions(isub,i)
    enddo

    !get runoff for this grid cell
    dataarray(5:7) = 0._realkind

    ! write(*,*) ' dimcell ',dimcell

    do i = 1, dimcell
       if ( abs(cmidrw(isub,i)) > 1.e-5_realkind) then 
          if(cmidrx(isub,i)>0)then      !rcagrid
             if(ntile==1)then
                dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * runoff(cmidrx(isub,i),cmidry(isub,i),1)
             else
                dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * runoff(cmidrx(isub,i),cmidry(isub,i),:)
             endif
          elseif(cmidrx(isub,i)<0)then   !gcmgrid
             dataarray(5:7) = dataarray(5:7) + cmidrw(isub,i) * gcmrunoff(-cmidrx(isub,i),-cmidry(isub,i))  !same runoff all tiles
          endif
       endif      ! cmidrw > 0
    enddo

    !irrigation not yet included
    irrigation = 0._realkind
    dataarray(8) = irrigation

    !set type of routing for the grid cell
    conduct = .false.
    if(landfrac>0.0_realkind) conduct(i_landriver)   = .true.
    if(lakefrac>0.0_realkind) conduct(i_lakeoutflow) = .true.

    !     write(*,*) ' isub ',isub

    !transport water through the grid cell
    call routgrid(isub,paths(isub),nsubs,area(isub),dataarray,par,irrigation,accinflow,gridoutflow,soil,riverq,lakew,conduct)



    !save discharge
    outflow1(isub) = gridoutflow / real(acctime,realkind)   !outflow = m3/s, gridoutflow=m3/timestep

 enddo basin_loop

 !        write(*,*) ' routcatchment end '

end	subroutine routcatchment

!------------------------------------------------------------------------------------------------------------------	
subroutine routgrid(i,ito,nsubs,gridarea,indata,parvalues,abstract,accinflow,outflow3,soil,riverq,lakew,conduct)
  !------------------------------------------------------------------------------------------------------------------	

  ! use routmodvar, only : i_landriver, i_lakeoutflow, dimriverq
 implicit none 

 !declaration of arguments
 integer, intent(in) :: i                  				!subbasin, routing grid cell
 integer, intent(in) :: ito                				!subbasin, recieving routing grid cell
 integer, intent(in) :: nsubs              				!number of subbasins, routing grid cells

 real(kind=realkind), intent(in)    :: gridarea           				!area of grid (km2)
 real(kind=realkind), intent(in)    :: indata(8)          				!input data needed by subroutine
 real(kind=realkind), intent(in)    :: parvalues(6)       				!parameter values for subbasin
 real(kind=realkind), intent(inout) :: abstract           				!allowed abstraction (out) (m3/timestep)
 real(kind=realkind), intent(inout) :: accinflow(nsubs)   				!accumulated routed flow (m3/timestep)
 real(kind=realkind), intent(out)   :: outflow3           				!outflow from routing grid cell (m3/timestep)
 real(kind=realkind), intent(inout) :: soil(nsubs)        				!water in soil shallow ground water (m3)
 real(kind=realkind), intent(inout) :: lakew(nsubs)       				!water level of lake (m)
 real(kind=realkind), intent(inout) :: riverq(0:dimriverq,nsubs)  !water in river (m3)
 logical, intent(in) :: conduct(2)                 !calculations to be made

 !declaration of local variables
 integer y,j
 real(kind=realkind)    x,qupstream,qin,transq
 real(kind=realkind)    ldepthm,outflowm,outflowm3s,rating1,rating2,soilrunoff
 real(kind=realkind)    lakearea   !m2
 real(kind=realkind)    runoff(3),tilefrac(3),lakefrac

 !preparations
 tilefrac = indata(1:3)    !forest,open,lakeland
 lakefrac = indata(4)      !outlet lake
 runoff   = indata(5:7)

 !routing
 !here begins the routing part. the primary units are m3 and m3/timestep. the waterlevel of the lakes are measured in m.
 !shallow ground water
 if(conduct(i_landriver))then
    qin = 0._realkind
    do j = 1,3
       qin = qin + runoff(j) * gridarea * tilefrac(j) * 1000._realkind      !runoff, m3/tidsteg, assuming runoff=mm/timestep
    enddo
    soil(i) = soil(i) + qin                    !m3
    soilrunoff = parvalues(1) * soil(i)        !m3/timesteg   !if timestep>1d we should have some loop here
    soil(i) = soil(i) - soilrunoff
 else
    soilrunoff = 0._realkind
 endif

 !flow from upstream routing grid cells
 qupstream = accinflow(i)   !m3/timestep

 !river
 if(conduct(i_landriver))then
    !inflow to river
    qin = soilrunoff + qupstream                       !!m3/timesteg

    !translation (delay) in river
    riverq(0,i) = qin       !add new inflow to translation variable
    if(parvalues(2)>0.0_realkind)then
       y = nint(parvalues(3))
       x = parvalues(4)
       transq = (1._realkind-x)*riverq(y,i) + x*riverq(y+1,i)   !calculate flow from river after translation
       riverq(1:y+1,i) = riverq(0:y,i)
    else
       transq = qin
    endif
 else
    transq = qupstream            !river outflow = transq       !m3/timestep
 endif

 !lake
 !"runoff" on lake
 qin = transq                                    !m3/timestep
 if(lakefrac<0.0_realkind) qin = qin + runoff(3) * gridarea * (1._realkind - sum(tilefrac)) * 1000._realkind      !lake runoff, m3/tidsteg

 if(conduct(i_lakeoutflow))then
    lakearea = lakefrac * gridarea * 1000000._realkind   !m2

    !inflow to lake
    lakew(i) = lakew(i) + qin / lakearea        !m

    !outflow from lake
    ldepthm = 0._realkind                 !threshold olake depth (m)
    rating1 = parvalues(5)
    rating2 = parvalues(6)
    if(rating1>0.0_realkind)then
       if(lakew(i)>ldepthm)then
          outflowm3s = rating1 * (lakew(i)-ldepthm)**rating2  !general rating curve (m3/s)   !if timestep>1d we should have some loop her
       else
          outflowm3s = 0._realkind
       endif
       outflowm = outflowm3s * real(acctime,realkind) / lakearea     !m3/s -> m/tidsteg
       if(outflowm>lakew(i)-ldepthm)then     !check
          if(lakew(i)>ldepthm)then
             outflowm = lakew(i)-ldepthm
          else
             outflowm = 0._realkind
          endif
       endif
    else
       if(lakew(i)>ldepthm)then
          outflowm = lakew(i)-ldepthm
       else
          outflowm = 0._realkind
       endif
    endif
    lakew(i) = lakew(i) - outflowm      !m/timestep
    outflow3 = outflowm * lakearea      !m3/timestep

 else    
    !no lake
    outflow3 = qin
 endif

 !transport outflow to next grid cell
 if(ito>0)  accinflow(ito) = accinflow(ito) + outflow3

end subroutine routgrid



subroutine read_namrouting(routing_path,ilen)
  use decomp
  implicit none
  character(len=132)::routing_path
  integer ilen
  character(len=20)namelistfile,inst
  namelist/institute/inst
  namelist/namrouting/klon_bound,klat_bound,ntile,dimstate,acctime
  namelist/rout_path/routing_path
  ilen=1
  if(mype==0)then
     open(57,file='namelists.dat',status='old',form='formatted')
     read(57,nml=namrouting)
     read(57,nml=institute)
     close(57)
     write(6,nml=namrouting)
     namelistfile='gcmpaths.'//trim(inst)
     open(61,file=namelistfile,status='old',form='formatted')
     read(61,nml=rout_path)
     close(61)
  endif

  ilen=index(routing_path,' ')-1
end subroutine read_namrouting



end module routmodvar


