module confys
  use decomp
  implicit none
  private
  real(kind=realkind),public,parameter::pi = 3.141592653589793_realkind 
  real(kind=realkind),public,parameter::omega  = 2.0_realkind*pi/(24.0_realkind*3600.0_realkind)
  real(kind=realkind),public,save::latvap = 2.5003e6_realkind
  real(kind=realkind),public,save::rair   = 2.8704e2_realkind
  real(kind=realkind),public,save::cpair  = 1.0046e3_realkind
  real(kind=realkind),public,save::ccpq   = 0.8593_realkind
  real(kind=realkind),public,save::epsilo = 0.622_realkind
  real(kind=realkind),public,save::gravit = 9.80665_realkind
  real(kind=realkind),public,save::tmelt  = 273.16_realkind
  real(kind=realkind),public,save::latice = 3.336e5_realkind
  real(kind=realkind),public,save::rhos   = 1.3e3_realkind
  real(kind=realkind),public,save::rhoh2o = 1.0e3_realkind
  real(kind=realkind),public,save::solar  = 1.367E3_realkind
  real(kind=realkind),public,save::stebol = 5.67e-8_realkind
  real(kind=realkind),public,save::carman = 0.40_realkind
  real(kind=realkind),public,save::rearth = 6.37e6_realkind
  logical,public,save::lbkf=.true.
  logical,public,save::maximum_random=.false.
  logical,public,save::lcarb=.true.
  logical,public,save::lmulch=.true.
  public conphys,ficefun,t2es,q2e,es2t !,adiagcnd
  interface ficefun
     module procedure ficefun1,ficefun0,ficefun2
  end interface
contains 
  !       saturation water vapour look-up-table
  !     function to compute dew point temperature (bolton's formula)
  !     !! warning: everything in mks units !!
  real(kind=realkind) function t2es(targ)
    implicit none
    real(kind=realkind),intent(in)::targ
    t2es=100.0_realkind*(6.112_realkind*exp(17.67_realkind*(targ-tmelt)/((targ-tmelt)+243.5_realkind)))
  end function t2es

  real(kind=realkind) function q2e(ptarg,qtarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg,qtarg
    q2e=ptarg*qtarg/(epsilo+qtarg*(1.0_realkind-epsilo))
  end function q2e

  real(kind=realkind) function es2t(earg)
    implicit none
    real(kind=realkind),intent(in)::earg
    es2t=tmelt+243.5_realkind*log((max(earg,1.e-14_realkind)/100.0_realkind)/6.112_realkind) &
         /(17.67_realkind-log((max(earg,1.e-14_realkind)/100.0_realkind)/6.112_realkind))
  end function es2t


  subroutine ficefun1(nhor,nlev,kstart,kstop,t,fice)
    ! compute fraction of ice in clouds
    ! call this subroutine from RADIA, VCBR, AKFRAK, and GEMINI

    implicit none

    integer,intent(in):: nhor,nlev,kstart,kstop
    real(kind=realkind),intent(in):: t(nhor,nlev)
    real(kind=realkind),intent(out):: fice(nhor,nlev)
    ! LOCAL
    integer i,k
    real(kind=8)::denom
    real(kind=8)::zbeta
    real(kind=realkind),parameter::tlow=250.16_realkind

    denom=1.0d0/(dble(tmelt-tlow))
    do k=1,nlev
       do i=kstart,kstop
          zbeta= dble((t(i,k)-tlow))*denom
          zbeta = max(0d0,min(1d0,zbeta))
          fice(i,k)=real(1d0-(zbeta*zbeta),realkind)
       enddo
    enddo

    return
  end subroutine ficefun1

  subroutine ficefun2(klon,klat,klev,t,fice)
    ! compute fraction of ice in clouds
    ! call this subroutine from RADIA, VCBR, AKFRAK, and GEMINI

    implicit none

    integer,intent(in):: klon,klat,klev
    real(kind=realkind),intent(in):: t(klon,klat,klev)
    real(kind=realkind),intent(out):: fice(klon,klat,klev)
    ! LOCAL
    integer i,j,k
    real(kind=8)::denom
    real(kind=8)::zbeta
    real(kind=realkind),parameter::tlow=250.16_realkind

    denom=1.0d0/(dble(tmelt-tlow))
    do k=1,klev
       do j=1,klat
          do i=1,klon
             zbeta= dble((t(i,j,k)-tlow))*denom
             zbeta = max(0d0,min(1d0,zbeta))
             fice(i,j,k)=real(1d0-(zbeta*zbeta),realkind)
          enddo
       enddo
    enddo

    return
  end subroutine ficefun2

  subroutine ficefun0(t,fice)
    ! compute fraction of ice in clouds
    ! call this subroutine from RADIA, VCBR, AKFRAK, and GEMINI
    implicit none

    real(kind=realkind),intent(in):: t
    real(kind=realkind),intent(out):: fice
    real(kind=8)::denom
    real(kind=8)::zbeta
    real(kind=realkind),parameter::tlow=250.16_realkind

    denom=1.0d0/(dble(tmelt-tlow))
    zbeta= dble((t-tlow))*denom
    zbeta = max(0d0,min(1d0,zbeta))
    fice=real(1d0-(zbeta*zbeta),realkind)

    return
  end subroutine ficefun0

  subroutine conphys(lecmwf)
    !subroutine conphys:   stefan gollvik feb 1991
    !give values to physical constants

    implicit none
    logical,intent(in)::lecmwf
#ifdef MPI_SRC
#include"mpif.h"
    integer ierr
    real(kind=realkind):: buf(18)
#endif
    namelist/namphys/latvap,rair,cpair,ccpq,epsilo,gravit,tmelt,&
         latice ,rhos,rhoh2o,solar,stebol,carman,rearth,lbkf,maximum_random,&
         lcarb,lmulch 

    if(lecmwf)then
       solar = 1.370E3_realkind
    endif

    if(mype==0) then
       open(57,file='namelists.dat',status='old')
       read(57,nml=namphys)
       close(57)
       write(6,nml=namphys)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namphys)
       close(1)
    endif

#ifdef MPI_SRC
    buf = (/latvap,rair,cpair,ccpq,epsilo,gravit,tmelt,latice,&
         rhos,rhoh2o,solar,stebol,carman,rearth,real(1.0,realkind),real(1.0,realkind),&
         real(1.0,realkind),real(1.0,realkind)/)
    if(.not.lbkf)then
       buf(15) = real(-1.0,realkind)
    endif
    if(.not.maximum_random)then
       buf(16) = real(-1.0,realkind)
    endif
    if(.not.lcarb)then
       buf(17) = real(-1.0,realkind)
    endif
    if(.not.lmulch)then
       buf(18) = real(-1.0,realkind)
    endif
    call MPI_Bcast(buf,18,REALTYPE,0,LOCALCOMM,ierr)
    latvap = buf(1)
    rair   = buf(2)
    cpair  = buf(3)
    ccpq   = buf(4)
    epsilo = buf(5)
    gravit = buf(6)
    tmelt  = buf(7) 
    latice = buf(8)
    rhos   = buf(9)
    rhoh2o = buf(10)
    solar  = buf(11)
    stebol = buf(12)
    carman = buf(13)
    rearth = buf(14)
    lbkf = buf(15) > real(0.0,realkind)
    maximum_random = buf(16) > real(0.0,realkind)
    lcarb = buf(17) > real(0.0,realkind)
    lmulch = buf(18) > real(0.0,realkind)
#endif
    

    return
  end subroutine conphys



!!$  subroutine adiagcnd (nhor,nlev,kstart,kstop,                   &
!!$       dtime,                                      &
!!$       acprcu,acprst,acsncu,acsnst,sevapr,sevapc,  &
!!$       kht,                                        &
!!$       cov2d,prcpcu,prcpst,                 &
!!$       dcwdin,dcwdt,cw,dpf,dtdt,hcucov,totcov,     &
!!$       stcov,stdcov,stcon,stdcon,stccov,           &
!!$       stscov,stcw,stcpnt,sttpnt,stscal)
!!$    !
!!$    ! diagcnd - diagnostics of sundqvist scheme
!!$    !
!!$    !
!!$
!!$    implicit none
!!$    !
!!$    !  ****** input 0-d *********
!!$    !
!!$    !  nhor   - number of horizontal gridpoints (nlon*nlat)
!!$    !  nlev   - number of vertical levels
!!$    !  kstart - starting gridpoint of the call
!!$    !  kstop  - final gridpoint of the call
!!$    !
!!$    integer nhor,nlev,kstart,kstop
!!$    !
!!$    !  dtime  - 2*dt
!!$    !  acprcu - accumulated precipitation from convective clouds
!!$    !  acprst - accumulated precipitation from stratiform clouds
!!$    !  acsncu - accumulated snow from convective clouds
!!$    !  acsnst - accumulated snow from stratiform clouds
!!$    !  sevapr - accumulated  evaporation from precipitation
!!$    !  sevapc - accumulated evaporation of cloud water
!!$    !
!!$    real(kind=realkind) dtime,acprcu,acprst,acsncu,acsnst,sevapr,sevapc
!!$    !
!!$    !  ****** input 1-d (nhor) *********
!!$    !
!!$    !  kht    - top of the convective clouds
!!$    !  prcpcu - convective precipitation (rain and snow)
!!$    !  prcpst - stratiform precipitation
!!$
!!$    !  cov2d  - 2d cloud cover from maximum/random overlapping
!!$    !
!!$    integer kht(nhor)
!!$    !
!!$    real(kind=realkind) prcpcu(nhor),prcpst(nhor),cov2d(nhor)
!!$    !
!!$    !  ****** input 2-d (nhor,nlev) *********
!!$    !
!!$    !  cw     - cloud water content
!!$    !  dcwdin - cloud water tendendy due to other processes than convection
!!$    !           (dynamics, radiation, turbulence fluxes, ..)
!!$    !  dpf    - delta p of the model layers
!!$    !  dtdt   - temperature tendency
!!$    !  dcwdt  - cloud water tendency
!!$    !  totcov - total cloud cover
!!$    !  hcucov - convective cloud cover
!!$    !
!!$    real(kind=realkind) cw(nhor,nlev),dcwdin(nhor,nlev),dpf(nhor,nlev), &
!!$         dtdt(nhor,nlev),dcwdt(nhor,nlev),totcov(nhor,nlev), &
!!$         hcucov(nhor,nlev)
!!$    !
!!$    !  ****** output 1-d ********
!!$    !
!!$    !  stcon  - tendency of the temperature
!!$    !  stdcon - (dtdt)**2
!!$    !  stcov  -
!!$    !  stdcov -
!!$    !  stccov - convective cloud cover
!!$    !  stscov - total cloud cover
!!$    !  stcw   - cloud water
!!$    !  stcpnt - number of columns where convection happen
!!$    !  sttpnt - number of columns with total cloud cover larger than zero
!!$    !  stscal - statistical array
!!$    !  stscal(1) -  accumulate the maximum convective cloud cover of each
!!$    !               column
!!$    !  stscal(2) -  accumulate the total 2d cloud cover
!!$    !  stscal(3) -  accumulate the cloud water vertically integrated
!!$    !  stscal(4) -  total number of points with convection in tha integratio
!!$    !               area
!!$    !  stscal(5) -  total number of points with cloud cover larger than zero
!!$    !  stscal(6) -  accumulate the convective precipitation (in tendency for
!!$    !  stscal(7) -  accumualte the stratiform precipitation (in tendency for
!!$    !  stscal(8) -  total convective precipitation
!!$    !  stscal(9) -  total stratiform precipitation
!!$    !  stscal(10) - total snow
!!$    !  stscal(11) - total evaporation of rain
!!$    !  stscal(12) - total evaporation of cloud water
!!$    !
!!$    real(kind=realkind) stcon(nlev),stdcon(nlev),stcov(nlev),stdcov(nlev), &
!!$         stccov(nlev),stscov(nlev),stcw(nlev),stcpnt(nlev), &
!!$         sttpnt(nlev),stscal(nlev)
!!$    !
!!$    !****************************************************************
!!$    !
!!$    !                       work space:
!!$    !
!!$    !****************************************************************
!!$    !
!!$    !
!!$    ! ********** here follows 1-d integer work arrays (nhor) *********
!!$
!!$    real(kind=realkind) wrk1(nhor),wrk2(nhor)
!!$
!!$    call diagcnd (nhor,nlev,kstart,kstop,                  &
!!$         dtime,                                     &
!!$         acprcu,acprst,acsncu,acsnst,sevapr,sevapc, &
!!$         kht,                                       &
!!$         cov2d,prcpcu,prcpst,                &
!!$         dcwdin,dcwdt,cw,dpf,dtdt,hcucov,totcov,    &
!!$         stcov,stdcov,stcon,stdcon,stccov,          &
!!$         stscov,stcw,stcpnt,sttpnt,stscal,          &
!!$         wrk1,wrk2)
!!$    return
!!$  end subroutine adiagcnd
!!$
!!$
!!$
!!$  subroutine diagcnd (nhor,nlev,kstart,kstop,      &
!!$       dtime,                                       &
!!$       acprcu,acprst,acsncu,acsnst,sevapr,sevapc,   &
!!$       kht,                                         &
!!$       cov2d,prcpcu,prcpst,                  &
!!$       dcwdin,dcwdt,cw,dpf,dtdt,hcucov,totcov,      &
!!$       stcov,stdcov,stcon,stdcon,stccov,            &
!!$       stscov,stcw,stcpnt,sttpnt,stscal,            &
!!$       wrk1,wrk2)
!!$    !     diagcnd - diagnostics of sundqvist scheme
!!$    implicit none
!!$    !     
!!$    !     ****** input 0-d *********
!!$    !     
!!$    !     nhor   - number of horizontal gridpoints (nlon*nlat)
!!$    !     nlev   - number of vertical levels
!!$    !     kstart - starting gridpoint of the call
!!$    !     kstop  - final gridpoint of the call
!!$    !     
!!$    integer nhor,nlev,kstart,kstop
!!$    !     
!!$    !     dtime  - 2*dt
!!$    !     acprcu - accumulated precipitation from convective clouds
!!$    !     acprst - accumulated precipitation from stratiform clouds
!!$    !     acsncu - accumulated snow from convective clouds
!!$    !     acsnst - accumulated snow from stratiform clouds
!!$    !     sevapr - accumulated  evaporation from precipitation
!!$    !     sevapc - accumulated evaporation of cloud water
!!$    !     
!!$    real(kind=realkind) dtime,acprcu,acprst,acsncu,acsnst,sevapr,sevapc
!!$    !     
!!$    !     ****** input 1-d (nhor) *********
!!$    !     
!!$    !     kht    - top of the convective clouds
!!$    !     prcpcu - convective precipitation (rain and snow)
!!$    !     prcpst - stratiform precipitation
!!$
!!$    !     cov2d  - 2d cloud cover from maximum/random overlapping
!!$    !     
!!$    integer kht(nhor)
!!$    !     
!!$    real(kind=realkind) prcpcu(nhor),prcpst(nhor),cov2d(nhor)
!!$    !     
!!$    !     ****** input 2-d (nhor,nlev) *********
!!$    !     
!!$    !     cw     - cloud water content
!!$    !     dcwdin - cloud water tendendy due to other processes than convection
!!$    !     (dynamics, radiation, turbulence fluxes, ..)
!!$    !     dpf    - delta p of the model layers
!!$    !     dtdt   - temperature tendency
!!$    !     dcwdt  - cloud water tendency
!!$    !     totcov - total cloud cover
!!$    !     hcucov - convective cloud cover
!!$    !     
!!$    real(kind=realkind) cw(nhor,nlev),dcwdin(nhor,nlev),dpf(nhor,nlev), &
!!$         dtdt(nhor,nlev),dcwdt(nhor,nlev),totcov(nhor,nlev), &
!!$         hcucov(nhor,nlev)
!!$
!!$    !     ****** output 1-d ********
!!$    !     
!!$    !     stcon  - tendency of the temperature
!!$    !     stdcon - (dtdt)**2
!!$    !     stcov  -
!!$    !     stdcov -
!!$    !     stccov - convective cloud cover
!!$    !     stscov - total cloud cover
!!$    !     stcw   - cloud water
!!$    !     stcpnt - number of columns where convection happen
!!$    !     sttpnt - number of columns with total cloud cover larger than zero
!!$    !     stscal - statistical array
!!$    !     stscal(1) -  accumulate the maximum convective cloud cover of each
!!$    !     column
!!$    !     stscal(2) -  accumulate the total 2d cloud cover
!!$    !     stscal(3) -  accumulate the cloud water vertically integrated
!!$    !     stscal(4) -  total number of points with convection in tha integratio
!!$    !     area
!!$    !     stscal(5) -  total number of points with cloud cover larger than zero
!!$    !     stscal(6) -  accumulate the convective precipitation (in tendency for
!!$    !     stscal(7) -  accumualte the stratiform precipitation (in tendency for
!!$    !     stscal(8) -  total convective precipitation
!!$    !     stscal(9) -  total stratiform precipitation
!!$    !     stscal(10) - total snow
!!$    !     stscal(11) - total evaporation of rain
!!$    !     stscal(12) - total evaporation of cloud water
!!$    real(kind=realkind) stcon(nlev),stdcon(nlev),stcov(nlev),stdcov(nlev), &
!!$         stccov(nlev),stscov(nlev),stcw(nlev),stcpnt(nlev), &
!!$         sttpnt(nlev),stscal(nlev)
!!$    real(kind=realkind) wrk1(nhor),wrk2(nhor)
!!$    integer jl,jk
!!$    real(kind=realkind) zwrk1,zwrk2,zwrk10,zwrk20,zwrk11,zwrk12,zwrk13, &
!!$         zwrk14,zwrk15,zwrk16,zwrk32,zwrk57
!!$
!!$    !     perform diagnostic calculations
!!$    !     the calculated quantities are explained under
!!$    !     "names of parameters and other quantities;
!!$
!!$    do 1930 jl = kstart,kstop
!!$       wrk1(jl) = 0.0_realkind
!!$       wrk2(jl) = 0.0_realkind
!!$1930 enddo
!!$
!!$    zwrk16 =  0.0_realkind
!!$    zwrk57 =  0.0_realkind
!!$
!!$    do 1950 jk = 1,nlev
!!$       zwrk1 =  0.0_realkind
!!$       zwrk2 =  0.0_realkind
!!$       zwrk10 =  0.0_realkind
!!$       zwrk20 =  0.0_realkind
!!$       zwrk11 =  0.0_realkind
!!$       zwrk12 =  0.0_realkind
!!$       zwrk13 =  0.0_realkind
!!$       zwrk14 =  0.0_realkind
!!$       zwrk15 =  0.0_realkind
!!$
!!$       !     zwrk32 = mean cloud water amount for each model layer
!!$       do 1940 jl = kstart,kstop
!!$          zwrk32 = ((dcwdt(jl,jk) + dcwdin(jl,jk)) *dtime + &
!!$               cw(jl,jk))*dpf(jl,jk)/gravit
!!$          if(zwrk32/(dpf(jl,jk)/gravit)<=1.e-6_realkind*cw(jl,jk) &
!!$               .or.zwrk32/(dpf(jl,jk)/gravit)<1.e-13_realkind) then
!!$             zwrk32 = 0.0_realkind
!!$          endif
!!$          if(jk>=kht(jl)) then
!!$             zwrk10 = zwrk10 + dtdt(jl,jk)
!!$             zwrk20 = zwrk20 + dtdt(jl,jk)*dtdt(jl,jk)
!!$          endif
!!$          zwrk1 = zwrk1 + dtdt(jl,jk)
!!$          zwrk2 = zwrk2 + dtdt(jl,jk)*dtdt(jl,jk)
!!$          zwrk11 = zwrk11 + hcucov (jl,jk)
!!$          zwrk12 = zwrk12 + totcov (jl,jk)
!!$          zwrk13 = zwrk13 +  zwrk32
!!$          zwrk14 = zwrk14 + real(int( hcucov(jl,jk) + 0.4999_realkind),realkind)
!!$          zwrk15 = zwrk15 + real(int( totcov(jl,jk) + 0.4999_realkind),realkind)
!!$          wrk2(jl) = max(wrk2(jl), hcucov(jl,jk))
!!$
!!$          if(hcucov(jl,jk)>0.0_realkind .and. abs(wrk1(jl))<1.e-14_realkind) then
!!$             wrk1(jl) = 1.0_realkind
!!$             zwrk16 = zwrk16 + 1.0_realkind
!!$          endif
!!$1940   enddo
!!$
!!$       !     stcon  - tendency of the temperature
!!$       !     stdcon - (dtdt)**2
!!$       !     stccov - convective cloud cover
!!$       !     stscov - total cloud cover
!!$       !     stcw   - cloud water
!!$       !     stcpnt - number of columns where convection happen
!!$       !     sttpnt - number of columns with total cloud cover larger tha
!!$       !     
!!$       stcov(jk) = stcov (jk)+ zwrk10
!!$       stdcov(jk)= stdcov(jk)+ zwrk20
!!$       stcon(jk)  = stcon (jk) + zwrk1
!!$       stdcon(jk) = stdcon(jk) + zwrk2
!!$       stccov(jk) = stccov (jk) + zwrk11
!!$       stscov(jk) = stscov (jk) + zwrk12
!!$       stcw(jk)   = stcw (jk) + zwrk13
!!$       stcpnt(jk) = stcpnt(jk) + zwrk14
!!$       sttpnt(jk) = sttpnt(jk) + zwrk15
!!$       stscal(5) = stscal(5) + zwrk15
!!$       stscal(3) = stscal(3) + stcw(jk)
!!$1950 enddo
!!$
!!$    !     stscal(1) -  accumulate the maximum convective cloud cover of each
!!$    !     column
!!$    !     stscal(2) -  accumulate the total 2d cloud cover
!!$    !     stscal(3) -  accumulate the cloud water vertically integrated
!!$    !     stscal(4) -  total number of points with convection in tha integratio
!!$    !     area
!!$    !     stscal(5) -  total number of points with cloud cover larger than zero
!!$    !     stscal(6) -  accumulate the convective precipitation (in tendency for
!!$    !     stscal(7) -  accumualte the stratiform precipitation (in tendency for
!!$    !     stscal(8) -  total convective precipitation
!!$    !     stscal(9) -  total stratiform precipitation
!!$    !     stscal(10) - total snow
!!$    !     stscal(11) - total evaporation of rain
!!$    !     stscal(12) - total evaporation of cloud water
!!$    !     
!!$    do 1960 jl = kstart, kstop
!!$       stscal(1) = stscal(1) + wrk2(jl)
!!$       stscal(2) = stscal(2) + cov2d(jl)
!!$       stscal(6) = stscal(6) + prcpcu(jl)
!!$       stscal(7) = stscal(7) + prcpst(jl)
!!$1960 enddo
!!$
!!$    stscal(4) = stscal(4) + zwrk16
!!$    stscal(8) = stscal(8) + acprcu
!!$    stscal(9) = stscal(9) + acprst
!!$    stscal(10) = stscal(10) + acsncu+acsnst
!!$    stscal(11) = stscal(11) + sevapr
!!$    stscal(12) = stscal(12) + sevapc
!!$    return
!!$  end subroutine diagcnd



end module confys
