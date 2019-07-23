module condensmod
  use decomp,only:realkind
  implicit none
  private
  real(kind=realkind),public,save::cevapcu(100)
  ! A1.+2 Geoide
  real(kind=realkind),public,save::ra
  real(kind=realkind),public,save::rg

  ! a1.4 thermodynamic gas phase
  real(kind=realkind),public,save:: r
  real(kind=realkind),public,save::rmd
  real(kind=realkind),public,save::rmv
  real(kind=realkind),public,save::rd
  real(kind=realkind),public,save::rv
  real(kind=realkind),public,save::rcpd
  real(kind=realkind),public,save:: retv
  !  thermodynamic transition of phase
  real(kind=realkind),public,save::rlvtt
  real(kind=realkind),public,save::rlstt
  real(kind=realkind),public,save::rlmlt
  real(kind=realkind),public,save::rtt
  real(kind=realkind),public,save::restt !  curve of saturation

  !DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
  real(kind=realkind),public,save::r2es
  real(kind=realkind),public,save::r3les
  real(kind=realkind),public,save::r3ies
  real(kind=realkind),public,save::r4les
  real(kind=realkind),public,save::r4ies
  real(kind=realkind),public,save::r5les
  real(kind=realkind),public,save::r5ies
  real(kind=realkind),public,save::r5alvcp
  real(kind=realkind),public,save::r5alscp
  real(kind=realkind),public,save::ralvdcp
  real(kind=realkind),public,save::ralsdcp
  real(kind=realkind),public,save::ralfdcp
  real(kind=realkind),public,save::rtwat
  real(kind=realkind),public,save::rtice

  real(kind=realkind),public,save::ENTRPEN=1.0E-4_realkind!AVERAGE ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
  real(kind=realkind),public,save::ENTRSCV=3.0E-4_realkind!AVERAGE ENTRAINMENT RATE FOR SHALLOW CONVECTION
  real(kind=realkind),public,save::ENTRMID=1.0E-4_realkind!AVERAGE ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
  real(kind=realkind),public,save::ENTRDD=2.0E-4_realkind!AVERAGE ENTRAINMENT RATE FOR DOWNDRAFTS
  real(kind=realkind),public,save::CMFCTOP=0.33_realkind!RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVE
  real(kind=realkind),public,save::CMFCMAX=1.0_realkind!MAXIMUM MASSFLUX VALUE ALLOWED FOR UPDRAFTS ETC
  real(kind=realkind),public,save::CMFCMIN=1.e-14_realkind!MINIMUM MASSFLUX VALUE (FOR SAFETY)
  real(kind=realkind),public,save::CMFDEPS=0.3_realkind!FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  real(kind=realkind),public,save::RHCDD=1.0_realkind!IS NO LONGER USED ( FORMULATION IMPLIES SATURATION)
  real(kind=realkind),public,save::CPRCON=2.0E-4_realkind!COEFFICIENTS FOR DETERMINING CONVERSION FROM CLOUD WATE
  logical,public,save::lmfpen
  logical,public,save::lmfscv
  logical,public,save::lmfmid
  logical,public,save::lmfdd
  logical,public,save::lmfdudv
  logical,public,save::leevap
  public sucst,suphec,sucond, sucumf
  public foelhm,foealfa,foeewm,foedem,foeldcpm


  public acondens,prcond
contains

  subroutine sucond ( klev , peta )

    !sucond - initialize common yoecnd controlling large-scale cond.
    !     
    !     purpose.
    !     --------
    !     initialize yoecnd, the common that controls the
    !     large-scale condensation routine of the model.


    implicit none
    integer:: jk,klev


    real(kind=realkind):: peta(klev)

    do 101 jk=1,klev
       if(leevap) then
          cevapcu(jk)=1.93e-6_realkind*261.0_realkind*sqrt(1.e3_realkind/(38.3_realkind*0.293_realkind) &
               *sqrt(peta(jk))) * 0.5_realkind/rg
       else
          cevapcu(jk)=0.0_realkind
       endif
101 enddo

    return
  end subroutine sucond

  subroutine sucumf
    implicit none

    entrpen=1.0e-4_realkind

    entrscv=3.0e-4_realkind

    entrmid=1.0e-4_realkind

    entrdd =2.0e-4_realkind

    cmfctop=0.33_realkind

    cmfcmax=1.0_realkind

    cmfcmin=1.e-14_realkind

    cmfdeps=0.3_realkind

    cprcon =2.0e-4_realkind

    rhcdd=1.0_realkind

    lmfpen  =.true.
    lmfscv  =.true.
    lmfmid  =.true.
    lmfdd   =.true.
    lmfdudv =.true.
    leevap  =.true.


    return
  end subroutine sucumf

  subroutine suphec(klev,ahyb,bhyb)
    !suphec - initialises physical constants of uncertain value.
    !               within the e.c.m.w.f. physics package
    !          this routine sets the values for the physical constants used
    !     in the parameterization routines whenever these values are not
    !     known well enough to forbid any tuning or whenever they are
    !     subject to an arbitrary choice of the modeller. these constants
    !     are distributed in common decks *yoexxxx* where xxxx corresponds
    !     to the individual physical parametrization
    use referenceParameters,only:sip0

    implicit none
    real(kind=realkind):: vp00,afull,bfull
    integer:: jk,klev

    real(kind=realkind):: ahyb(klev+1) , bhyb(klev+1) , zeta(klev)
    !          0.1    defining universal constants for ecmwf model
    call sucst
    !          0.2    defining derived constants from universal constants
    r2es=restt*rd/rv
    r3les=17.269_realkind
    r3ies=21.875_realkind
    r4les=35.86_realkind
    r4ies=7.66_realkind
    r5les=r3les*(rtt-r4les)
    r5ies=r3ies*(rtt-r4ies)
    r5alvcp=r5les*rlvtt/rcpd
    r5alscp=r5ies*rlstt/rcpd
    ralvdcp=rlvtt/rcpd
    ralsdcp=rlstt/rcpd
    ralfdcp=rlmlt/rcpd
    rtwat=rtt
    rtice=rtt-23.0_realkind
    !         0.5    define standard atmosphere vertical configuration
    vp00=sip0

    do 52 jk=1,klev
       afull = 0.5_realkind*( ahyb(jk) + ahyb(jk+1) )
       bfull = 0.5_realkind*( bhyb(jk) + bhyb(jk+1) )
       zeta(jk)= afull/vp00 + bfull
52  enddo

    !          2.     setting constants for convection scheme
    call sucumf
    !          3.     setting constants for large-scale condensation scheme
    call sucond(klev,zeta)
    return
  end subroutine suphec

  subroutine sucst()
    use confys
    implicit none
    !
    !sucst  - routine to initialize the constants of the model.
    !
    !     purpose.
    !     --------
    !           initialize and print the common yomcst + initialize
    !         date and time of yomrip.
    !
    !     reference.
    !     ----------
    !        ecmwf research department documentation of the ifs
    !
    !     author.
    !     -------
    !        mats hamrud and philippe courtier  *ecmwf*
    !
    real(kind=realkind):: rkbol,rnavo

    !
    rkbol=1.380658e-23_realkind
    rnavo=6.0221367e+23_realkind
    !     rg=9.80665
    rg=gravit
    !     ra=6371229.
    ra=rearth
    !
    !     ------------------------------------------------------------------
    !
    !*       5.    define thermodynamic constants, gas phase.
    !              ------------------------------------------
    !

    !
    r=rnavo*rkbol
    rmd=28.9644_realkind
    rmv=18.0153_realkind
    !     rd=1000.*r/rmd
    rd=rair
    rv=1000.0_realkind*r/rmv
    !     rcpd=3.5*rd
    rcpd=cpair
    retv=rv/rd-1.0_realkind
    !
    !     ------------------------------------------------------------------
    !
    !*       8.    define thermodynamic constants, transition of phase.
    !              ----------------------------------------------------
    !

    !
    rtt=tmelt
    !     rlvtt=2.5008e+6
    rlvtt=latvap
    !     rlmlt=rlstt-rlvtt
    rlmlt=latice
    !     rlstt=2.8345e+6
    rlstt=latvap+latice
    !
    !     ------------------------------------------------------------------
    !
    !*       9.    saturated vapour pressure.
    !              --------------------------
    !

    !
    restt=611.14_realkind
    !
    !     ------------------------------------------------------------------
    !
    return
  end subroutine sucst


  real(kind=realkind) function foedelta(ptarg) 
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foedelta = max (0.0_realkind,sign(1.0_realkind,ptarg-rtt))
  end function foedelta

  real(kind=realkind) function foeew(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foeew = r2es*exp (    &
         (r3les*foedelta(ptarg)+r3ies*(1.0_realkind-foedelta(ptarg)))*(ptarg-rtt) &
         / (ptarg-(r4les*foedelta(ptarg)+r4ies*(1.0_realkind-foedelta(ptarg)))))
  end function foeew

  real(kind=realkind) function foede(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foede  = &
         (foedelta(ptarg)*r5alvcp+(1.0_realkind-foedelta(ptarg))*r5alscp) &
         / (ptarg-(r4les*foedelta(ptarg)+r4ies*(1.0_realkind-foedelta(ptarg))))**2.0_realkind
  end function foede

  real(kind=realkind) function foedesu(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foedesu  = &
         (foedelta(ptarg)*r5les+(1.0_realkind-foedelta(ptarg))*r5ies) &
         / (ptarg-(r4les*foedelta(ptarg)+r4ies*(1.0_realkind-foedelta(ptarg))))**2.0_realkind
  end function foedesu

  real(kind=realkind) function foelh(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foelh = foedelta(ptarg)*rlvtt + (1.0_realkind-foedelta(ptarg))*rlstt
  end function foelh
  
  real(kind=realkind) function foeldcp(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foeldcp = foedelta(ptarg)*ralvdcp + (1.0_realkind-foedelta(ptarg))*ralsdcp
  end function foeldcp
  
!     foealfa is calculated to distinguish the three cases:
!
!                       foealfa=1            water phase
!                       foealfa=0            ice phase
!                       0 < foealfa < 1      mixed phase
!
!               input : ptarg = temperature
!
  real(kind=realkind) function foealfa(ptarg)
     implicit none
     real(kind=realkind),intent(in)::ptarg
    foealfa = min(1.0_realkind,((max(rtice,min(rtwat,ptarg))-rtice) &
         /(rtwat-rtice))**2.0_realkind)
  end function foealfa


!     pressure of water vapour at saturation
!        input : ptarg = temperature
!
  real(kind=realkind) function foeewm(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foeewm = r2es * &
         (foealfa(ptarg)*exp(r3les*(ptarg-rtt)/(ptarg-r4les))+ &
         (1.0_realkind-foealfa(ptarg))*exp(r3ies*(ptarg-rtt)/(ptarg-r4ies)))
  end function foeewm

  real(kind=realkind) function foedem(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foedem = foealfa(ptarg)*r5alvcp*(1._realkind/(ptarg-r4les)**2_realkind)+ &
                   (1.0_realkind-foealfa(ptarg))*r5alscp*(1._realkind/(ptarg-r4ies)**2_realkind)
  end function foedem

  real(kind=realkind) function foeldcpm(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foeldcpm = foealfa(ptarg)*ralvdcp+(1.0_realkind-foealfa(ptarg))*ralsdcp
  end function foeldcpm

  real(kind=realkind) function foelhm(ptarg)
    implicit none
    real(kind=realkind),intent(in)::ptarg
    foelhm =foealfa(ptarg)*rlvtt+(1.0_realkind-foealfa(ptarg))*rlstt
  end function foelhm



  subroutine acondens(nhor,nlev,kstart,kstop,            &
       dtime,                                  &
       frland,ts,ps,dpsdin,                    &
       bhyb,                                   &
       bfull,                                  &
       ph,gpot,senf,latf,u,v,omf,              &
       t,q,cw,pf,dpf,                          &
       dcwdin,dtdtin,dqdtin,dudtin,dvdtin,     &
       cov2d,cwpath,kht,                       &
       sevapr,sevapc,                          &
       prcpst,stsnow,cusnow,prcpcu,            &
       draindt,dsnowdt,                        &
       dtdt,dqdt,dcwdt,dudt,dvdt,totcov,hcucov,&
       preta,prsl,prsi,prcl,prci)              
    !
    ! Sundqvist scheme
    !
    implicit none
    !
    !
    !

    !
    !  nhor   - number of horizontal gridpoints (nlon*nlat)
    !  nlev   - number of vertical levels
    !  kstart - starting gridpoint of the call
    !  kstop  - final gridpoint of the call
    !
    integer:: nhor,nlev,kstart,kstop
    !
    !  dtime  - 2*dt
    !
    real(kind=realkind):: dtime
    !

    !
    !  ts     - surface temperature
    !  ps     - surface pressure
    !  dpsdin - surface pressure tendendy due to dynamics
    !  senf   - surface sensible heat flux
    !  latf   - surface latent heat flux
    !
    real(kind=realkind):: ts(nhor),ps(nhor),dpsdin(nhor)
    real(kind=realkind):: senf(nhor),latf(nhor)
    !

    !
    !  bhyb   - parameter b of the hybrid coordinate of the half model level
    !
    real(kind=realkind):: bhyb(nlev+1)
    !au
    real(kind=realkind):: bfull(nlev)
    !

    !
    !  t      - temperature
    !  q      - specific humidity
    !  cw     - cloud water content
    !  pf     - pressure at the full model levels
    !  dpf    - delta p of the model layers
    !  gpot   - geopotential
    !  u      - wind x-component
    !  v      - wind y-component
    !  omf    - vertical velocity at full levels
    !  dcwdin - cloud water tendendy due to other processes than convection
    !           (dynamics, radiation, turbulence fluxes, ..)
    !  dtdtin - temperature tendendy equal than dcwdin
    !  dqdtin - q tendendy equal than dcwdin
    !  dudtin - wind x-component tendendy coming from dynamics
    !  dvdtin - wind y-component tendendy coming from dynamics
    !
    real(kind=realkind):: t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),pf(nhor,nlev), dpf(nhor,nlev)
    real(kind=realkind):: gpot(nhor,nlev),u(nhor,nlev),v(nhor,nlev),omf(nhor,nlev)
    real(kind=realkind):: dcwdin(nhor,nlev),dtdtin(nhor,nlev),dqdtin(nhor,nlev),&
         dudtin(nhor,nlev),dvdtin(nhor,nlev)
    !

    !
    !  ph     - pressure at the half model layers
    real(kind=realkind):: ph(nhor,nlev+1)
    !

    !
    !  frland  - fraction land
    !
    real(kind=realkind):: frland(nhor)
    !

    !
    !  cov2d  - 2d cloud cover from maximum/random overlapping
    !  cwpath - cloud water vertically integrated
    !  dsnowdt- tendency of snow
    !  draindt- tendency of liquid precipitation
    !
    real(kind=realkind):: cov2d(nhor),cwpath(nhor),dsnowdt(nhor),draindt(nhor)
    !

    !
    !  dtdt   - temperature tendency
    !  dqdt   - q tendency
    !  dcwdt  - cloud water tendency
    !  dudt   - wind x-component tendency
    !  dvdt   - wind y-component tendency
    !  totcov - total cloud cover
    !  hcucov - convective cloud cover
    !  preta  - precip. rate at model levels
    !
    real(kind=realkind):: dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev),&
         totcov(nhor,nlev),hcucov(nhor,nlev),preta(nhor,nlev)
    real(kind=realkind):: dudt(nhor,nlev),dvdt(nhor,nlev)
    !
    !                    fluxes for arpege interface (nhor,0:nlev):
    !  prsl   - flux of liquid precipitation due to stratiform part
    !  prsi   - flux of solid precipitation due to stratiform part
    !  prcl   - flux of liquid precipitation due to convection
    !  prci   - flux of solid precipitation due to convection

    real(kind=realkind):: prsl(nhor,0:nlev),prsi(nhor,0:nlev),&
         prcl(nhor,0:nlev),prci(nhor,0:nlev)
    !
    !

    !
    !  acprcu - accumulated precipitation from convective clouds
    !  acprst - accumulated precipitation from stratiform clouds
    !  acsncu - accumulated snow from convective clouds
    !  acsnst - accumulated snow from stratiform clouds
    !  sevapr - accumulated evaporation from precipitation
    !  sevapc - accumulated evaporation of cloud water
    !
    real(kind=realkind):: acprcu,acprst,acsncu,acsnst,sevapr,sevapc
    !
    !

    !
    !  prcpst - stratiform precipitation
    !  stsnow - stratiform snow
    !  cusnow - convective snow
    !  prcpcu - convective precipitation (rain and snow)
    !
    real(kind=realkind):: prcpst(nhor),stsnow(nhor),cusnow(nhor),prcpcu(nhor)
    !
    !

    !
    !  kht    - top of the convective clouds
    !
    integer:: kht(nhor)
    !

    !
    !                       work space:
    !

    !
    !

    !
    !  khb    - base of the convective clouds
    !  jwanv  - gridpoints where anvil is allowed to becomes stratiform
    !

    !
    !  jwmask - set = 0 if convection,   otherwise = 1
    !

    !
    !  frland - fraction of land of the grid square
    !
    !

    !
    !  hldcp2 - 1/cpair * latent heat as a function of
    !           temperature by weighting in latent heat of
    !           freezing times probability of ice crystals
    !  dlnpdt - = 1/p*dp/dt at level jk
    !  hdcwcu - tendency of cloud water due to convection
    !  hdcwst - tendency of cloud water due to stratiform
    !  hqsat  - saturation mixing ratio with respect to gridpoint temperatur
    !  hsq    - eps*l**2*qsat/(r*cp*t**2)
    !  hu     - relative humidity
    !  pcond  - condensate from mass flux to stamic
    !  pcondo - condensate from mass flux to stamic not update by stamic
    !



    integer:: khb(nhor),jwanv(nhor)
    integer:: jwmask(nhor,nlev)

    real(kind=realkind):: hldcp2(nhor,nlev),dlnpdt(nhor,nlev),&
         hdcwcu(nhor,nlev),hdcwst(nhor,nlev),&
         hqsat(nhor,nlev),hsq(nhor,nlev),hu(nhor,nlev),&
         pcond(nhor,nlev),pcondo(nhor,nlev)



    integer:: jl,jk,jtyp

    !  setting to zero local arrays

    do jl=kstart,kstop
       khb(jl)=0
       jwanv(jl)=0
    enddo
    !
    do jk=1,nlev
       do jl=kstart,kstop
          jwmask(jl,jk)=0
          hldcp2(jl,jk)=0.0_realkind
          dlnpdt(jl,jk)=0.0_realkind
          hdcwcu(jl,jk)=0.0_realkind
          hdcwst(jl,jk)=0.0_realkind
          hqsat(jl,jk)=0.0_realkind
          hsq(jl,jk)=0.0_realkind
          hu(jl,jk)=0.0_realkind
       enddo
    enddo
    !
    do jk=1,nlev
       do jl=kstart,kstop
          dtdt(jl,jk)=0.0_realkind
          dqdt(jl,jk)=0.0_realkind
          dcwdt(jl,jk)=0.0_realkind
          dudt(jl,jk)=0.0_realkind
          dvdt(jl,jk)=0.0_realkind
          pcond(jl,jk)=0.0_realkind
          pcondo(jl,jk)=0.0_realkind
       enddo
    enddo
    !
    !  do some computations before call the routines
    !
    call prcond (nhor,nlev,kstart,kstop,t,hldcp2)

    !  call to the convective part of the scheme
    !
    call acumastr   (nhor,nlev,kstart,kstop,                &
         dtime,                                   &
         frland,pf,ph,dpf,gpot,ts,senf,latf,hqsat,&
         t,q,cw,u,v,omf,                          &
         dcwdin,dtdtin,dqdtin,dudtin,dvdtin,      &
         jwanv,kht,khb,jwmask,hcucov,hdcwcu,      &
         dcwdt,dtdt,dqdt,dudt,dvdt,pcond)
    !
    !  store pcond to compute dtdt after stamic
    !
    do jk=1,nlev
       do jl=kstart,kstop
          pcondo(jl,jk)=pcond(jl,jk)
       enddo
    enddo
    !

    !  !!  WARNING  !!
    !
    !  dtdt, dqdt, dcwdt are both input and output to stamic,
    !  because they store tendencies due to convection which are
    !  updated due to evaporation in the microphysics part
    !
    !
    !
    !  call to the stratiform and microphysics
    !
    call astamic (nhor,nlev,kstart,kstop,                   &
         dtime,                                     &
         jwanv,kht,khb,                             &
         jwmask,                                    &
         frland,ts,ps,                              &
         bhyb,                                      &
         t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,        &
         hcucov,hdcwcu,dlnpdt,hldcp2,hu,hqsat,hsq,  &
         pcond,dqdt,dcwdt,                          &
         sevapr,sevapc,                             &
         totcov,prcpcu,prcpst,cusnow,stsnow,cov2d,  &
         cwpath,                                    &
         draindt,dsnowdt,                           &
         hdcwst,                                    &
         prcl,prci,prsl,prsi)

    !  precip. rate at model levels
    !
    do jk=1,nlev
       do jl=kstart,kstop
          preta(jl,jk)=prsl(jl,jk)+prsi(jl,jk) +prcl(jl,jk)+prci(jl,jk)
       enddo
    enddo
    !
    !  update dtdt after stamic
    !
    do jk=1,nlev
       do jl=kstart,kstop
          dtdt(jl,jk)=dtdt(jl,jk)-pcondo(jl,jk)+pcond(jl,jk)
       enddo
    enddo
    !
    return
  end subroutine acondens


  subroutine acumastr (nhor,nlev,kstart,kstop,                 &
       dtime,                                   &
       frland,pf,ph,dpf,gpot,ts,senf,latf,hqsat,&
       t,q,cw,u,v,omf,                          &
       dcwdin,dtdtin,dqdtin,dudtin,dvdtin,      &
       jwanv,kctop,kcbot,jwmask,hcucov,hdcwcu,  &
       dcwdt,dtdt,dqdt,dudt,dvdt,pcond)

    use confys
    use escom
    use condsmod
    implicit none
    !
    !  input 0-D
    !
    integer:: nhor,nlev,kstart,kstop
    real(kind=realkind):: dtime
    !
    !  input 1-D
    !
    real(kind=realkind):: frland(nhor),senf(nhor),latf(nhor),ts(nhor)
    !
    !  input 2-D
    !
    real(kind=realkind):: t(nhor,nlev),q(nhor,nlev),u(nhor,nlev),v(nhor,nlev),&
         omf(nhor,nlev),gpot(nhor,nlev),cw(nhor,nlev),        &
         pf(nhor,nlev),ph(nhor,nlev+1),dpf(nhor,nlev),        &
         hqsat(nhor,nlev)                                      
    real(kind=realkind):: dtdtin(nhor,nlev),dqdtin(nhor,nlev),dudtin(nhor,nlev),&
         dvdtin(nhor,nlev),dcwdin(nhor,nlev)
    !
    !  output 1-D
    !
    integer:: jwanv(nhor),kcbot(nhor),kctop(nhor)
    !
    !  output 2-D
    !
    integer:: jwmask(nhor,nlev)
    real(kind=realkind):: dtdt(nhor,nlev),dqdt(nhor,nlev),dudt(nhor,nlev),&
         dvdt(nhor,nlev),dcwdt(nhor,nlev),pcond(nhor,nlev)
    real(kind=realkind):: hcucov(nhor,nlev),hdcwcu(nhor,nlev)
    !
    ! work space:
    !

    !  here follows 0-D integer work numbers
    integer:: jl,jk
    !  here follows 0-D real work numbers
    real(kind=realkind):: zlat,c0,c1,c2,zmfu
    !  here follows 1-D logical work arrays
    logical:: ldland(nhor),ldcum(nhor)
    !  here follows 1-D integer work arrays
    integer:: ktype(nhor)
    !  here follows 1-D real work arrays
    real(kind=realkind):: qhfl(nhor),rain(nhor)
    !  here follows 2-D real work arrays
    real(kind=realkind):: tu(nhor,nlev),qu(nhor,nlev),plu(nhor,nlev),     &
         pld(nhor,nlev),pleen(nhor,nlev),                 &
         plude(nhor,nlev),penth(nhor,nlev),mfu(nhor,nlev),&
         mfd(nhor,nlev)
    real(kind=realkind):: mflxr(nhor,nlev+1),mflxs(nhor,nlev+1)



    ! setting tendencies because cumastr updates them
    !
    do jl=kstart,kstop
       do jk=1,nlev
          dtdt(jl,jk)=dtdtin(jl,jk)
          dqdt(jl,jk)=dqdtin(jl,jk)
          dcwdt(jl,jk)=dcwdin(jl,jk)
          dudt(jl,jk)=dudtin(jl,jk)
          dvdt(jl,jk)=dvdtin(jl,jk)
       enddo
    enddo
    !
    ! logical array for land or sea
    !
    do jl=kstart,kstop
       if(frland(jl)>=0.5_realkind) then
          ldland(jl)=.true.
       else
          ldland(jl)=.false.
       endif
    enddo
    !
    ! moisture surface flux ( latf / zlat )
    !
    do jl=kstart,kstop
       zlat=latvap
       if(ts(jl)<=tmelt) zlat=zlat+latice
       qhfl(jl)=latf(jl)/zlat
    enddo

    call cumastr(kstart,kstop,nhor,1,nlev,         &
         ldland,dtime,                   &
         t,q,cw,u,v,omf,                 &
         hqsat,qhfl,senf,                &
         pf,ph,gpot,                     &
         dtdt,dqdt,dcwdt,dudt,dvdt,      &
         ldcum,ktype,kcbot,kctop,        &
         tu,qu,plu,pld,plude,pleen,penth,&
         mflxr,mflxs,rain,mfu,mfd,pcond)

    ! setting back tendencies because cumastr updates them
    !
    do jl=kstart,kstop
       do jk=1,nlev
          dtdt(jl,jk)=dtdt(jl,jk)-dtdtin(jl,jk)
          dqdt(jl,jk)=dqdt(jl,jk)-dqdtin(jl,jk)
          dcwdt(jl,jk)=dcwdt(jl,jk)-dcwdin(jl,jk)
          dudt(jl,jk)=dudt(jl,jk)-dudtin(jl,jk)
          dvdt(jl,jk)=dvdt(jl,jk)-dvdtin(jl,jk)
       enddo
    enddo

    ! interface for stamic
    !
    do jl=kstart,kstop
       !
       ! setting default value for kctop in agreement with stamic
       !
       if(kctop(jl)==nlev-1) kctop(jl)=nlev+1
       !
       do jk=1,nlev
          !
          ! jwmask=0 if convection was allowed and in sublcoud layers
          !          otherwise is 1
          !
          if(ktype(jl)/=0.and.jk>=kctop(jl))then
             jwmask(jl,jk)=0
          else
             jwmask(jl,jk)=1
          endif
          !
          ! hcucov
          !
          if(ktype(jl)/=0.and.jk>=kctop(jl).and. jk<=kcbot(jl))then
             !
             !  cloudiness as a function of mass flux (mb/h)
             !
             !  Xu K. and Krueger S.K. (1991).-'Evaluation of Cloudiness
             !  Paameterizations Using a Cumulus ensemble Model'.
             !  M.W.R., 119, 342-362.
             !
             zmfu=mfu(jl,jk)*gravit*3600.0_realkind/100.0_realkind
             if(pf(jl,jk)<=55000.0_realkind) then
                c0=0.0708_realkind
                c1=0.0941_realkind
                c2=0.0233_realkind
             elseif(pf(jl,jk)<=80000.0_realkind) then
                c0=0.0722_realkind
                c1=0.0760_realkind
                c2=0.0254_realkind
             else
                c0=0.0337_realkind
                c1=0.0453_realkind
                c2=0.0237_realkind
             endif
             if(zmfu>0.01_realkind) then
                hcucov(jl,jk)=c0+c1*log10(zmfu)+c2*(log10(zmfu))**2.0_realkind
             else
                hcucov(jl,jk)=0._realkind
             endif
             hcucov(jl,jk)=min(hcucov(jl,jk),0.25_realkind+0.5_realkind*q(jl,jk)/hqsat(jl,jk))
             hcucov(jl,jk)=max(hcucov(jl,jk),0.01_realkind)
          else
             hcucov(jl,jk)=0.0_realkind
          endif
          !
          ! hdcwcu
          !
          hdcwcu(jl,jk)=dcwdt(jl,jk)
          !
       enddo
       !
       ! jwanv=1 if stratiform condensation in anvil is allowed
       !         (only for deep convection), otherwise is 0
       !
       if(ktype(jl)==1.or.ktype(jl)==3) then
          if(t(jl,min(kctop(jl),nlev))<=tanvil) then
             jwanv(jl)=1
          else
             jwanv(jl)=0
          endif
       endif
       !
    enddo
    !
    return
  end subroutine acumastr


  subroutine cumastr(kidia,kfdia,klon,ktdia,    klev   &
       ,  ldland,   ptsphy                              &
       ,  pten,     pqen,     plen,     puen,     pven  &
       ,  pvervel,  pqsen,    pqhfl,    pahfs           &
       ,  pap,      paph,     pgeo                      &
       ,  ptent,    ptenq,    ptenl,    ptenu,    ptenv &
       ,  ldcum,    ktype,    kcbot,    kctop           &
       ,  ptu,      pqu,      plu,      pld,      plude &
       ,  pleen,    penth,    pmflxr,   pmflxs,   prain &
       ,  pmfu,     pmfd,     pcond )
    !     
    !cumastr  master routine for cumulus massflux-scheme
    !     
    !     m.tiedtke      e.c.m.w.f.     1986/1987/1989
    !     
    !     
    !     purpose
    !     -------
    !     
    !     this routine computes the physical tendencies of the
    !     prognostic variables t,q,u and v due to convective processes.
    !     processes considered are: convective fluxes, formation of
    !     precipitation, evaporation of falling rain below cloud base,
    !     saturated cumulus downdrafts.
    !     
    !   interface.
    !     ----------
    !     
    !     *cumastr* is called from *cucall*
    !     the routine takes its input from the long-term storage
    !     t,q,u,v,phi and p and moisture tendencies.
    !     it returns its output to the same space
    !     1.modified tendencies of model variables
    !     2.rates of convective precipitation
    !     (used in subroutine surf)
    !     3.cloud base, cloud top and precip for radiation
    !     (used in subroutine cloud)
    !     
    !     method
    !     -------
    !     
    !     parameterization is done using a massflux-scheme.
    !     (1) define constants and parameters
    !     (2) specify values (t,q,qs...) at half levels and
    !     initialize updraft- and downdraft-values in 'cuini'
    !     (3) calculate cloud base in 'cubase'
    !     and specify cloud base massflux from pbl moisture budget
    !     (4) do cloud ascent in 'cuasc' in absence of downdrafts
    !     (5) do downdraft calculations:
    !     (a) determine values at lfs in 'cudlfs'
    !     (b) determine moist descent in 'cuddraf'
    !     (c) recalculate cloud base massflux considering the
    !     effect of cu-downdrafts
    !     (6) do final cloud ascent in 'cuasc'
    !     (7) do final adjusments to convective fluxes in 'cuflx',
    !     do evaporation in subcloud layer
    !     (8) calculate increments of t and q in 'cudtdq'
    !     (9) calculate increments of u and v in 'cududv'
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        end point
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     
    !     input parameters (logical)
    !     
    !     *ldland*       land sea mask (.true. for land)
    !     
    !     input parameters (real)
    !     
    !     *ptsphy*       time step for the physics                       s
    !     *pten*         provisional environment temperature (t+1)       k
    !     *pqen*         provisional environment spec. humidity (t+1)  kg/kg
    !     *plen*         provisional env. cloud water cont. (t+1)      kg/kg
    !     *puen*         provisional environment u-velocity (t+1)       m/s
    !     *pven*         provisional environment v-velocity (t+1)       m/s
    !     *pvervel*      vertical velocity                             pa/s
    !     *pqsen*        environment spec. saturation humidity (t+1)   kg/kg
    !     *pqhfl*        surface moisture flux                       kg/(sm2)
    !     *pahfs*        surface sensible heat flux                   w/m2
    !     *pap*          provisional pressure on full levels             pa
    !     *paph*         provisional pressure on half levels             pa
    !     *pgeo*         geopotential                                  m2/s2
    !     
    !     updated parameters (real):
    !     
    !     *ptent*        temperature tendency                           k/s
    !     *ptenq*        moisture tendency                             kg/(kg s)
    !     *ptenl*        cloud water cont. tendency                    kg/(kg s)
    !     *ptenu*        tendency of u-comp. of wind                    m/s2
    !     *ptenv*        tendency of v-comp. of wind                    m/s2
    !     
    !     output parameters (logical):
    !     
    !     *ldcum*        flag: .true. for convective points
    !     
    !     output parameters (integer):
    !     
    !     *ktype*        type of convection
    !     1 = penetrative convection
    !     2 = shallow convection
    !     3 = midlevel convection
    !     *kcbot*        cloud base level
    !     *kctop*        cloud top level
    !     
    !     output parameters (real):
    !     
    !     *ptu*          temperature in updrafts                         k
    !     *pqu*          spec. humidity in updrafts                    kg/kg
    !     *plu*          liquid water content in updrafts              kg/kg
    !     *pld*          liquid water content in downdrafts            kg/kg
    !     *plude*        detrained liquid water                        kg/(m3*s)
    !     *pleen*        entrained liquid water                        kg/(m3*s)
    !     *penth*        increment of dry static energy                 j/(kg*s)
    !     *pmflxr*       convective rain flux                          kg/(m2*s)
    !     *pmflxs*       convective snow flux                          kg/(m2*s)
    !     *prain*        total precip. produced in conv. updrafts      kg/(m2*s)
    !     (no evaporation in downdrafts)
    !     *pmfu*         massflux updrafts                             kg/(m2*s)
    !     *pmfd*         massflux downdrafts                           kg/(m2*s)
    !     *pcond*        condensate                                     k/s
    !     
    !     externals.
    !     ----------
    !     
    !     cuini:  initializes values at vertical grid used in cu-parametr.
    !     cubase: cloud base calculation for penetr.and shallow convection
    !     cuasc:  cloud ascent for entraining plume
    !     cudlfs: determines values at lfs for downdrafts
    !     cuddraf:does moist descent for cumulus downdrafts
    !     cuflx:  final adjustments to convective fluxes (also in pbl)
    !     cudqdt: updates tendencies for t and q
    !     cududv: updates tendencies for u and v
    !     
    !     switches.
    !     --------
    !     
    !     lmfpen=.true.   penetrative convection is switched on
    !     lmfscv=.true.   shallow convection is switched on
    !     lmfmid=.true.   midlevel convection is switched on
    !     lmfdd=.true.    cumulus downdrafts switched on
    !     lmfdudv=.true.  cumulus friction switched on
    !     
    !     
    !     model parameters (defined in subroutine cuparam)
    !     ------------------------------------------------
    !     entrpen    entrainment rate for penetrative convection
    !     entrscv    entrainment rate for shallow convection
    !     entrmid    entrainment rate for midlevel convection
    !     entrdd     entrainment rate for cumulus downdrafts
    !     cmfctop    relative cloud massflux at level above nonbuoyancy leve
    !     cmfcmax    maximum massflux value allowed for
    !     cmfcmin    minimum massflux value (for safety)
    !     cmfdeps    fractional massflux for downdrafts at lfs
    !     cprcon     coefficient for conversion from cloud water to rain
    !     
    !     reference.
    !     ----------
    !     
    !     paper on massflux scheme (tiedtke,1989)
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)
    !     
    !----------------------------------------------------------------------
    !     

    implicit none


    integer:: kidia,kfdia,jk,jl,ktdia,ikb,itopm2,klon,klev
    real(kind=realkind):: zcons2,ptsphy,zqumqe,zdqmin,zhfl,zdh,zmf2,zmfmax,zal,zqal
    real(kind=realkind):: zalfaw,zalfai,z5ldcp,z4es,zhsat,zgam,zzz,zhhat,zpbmpt,zeps
    real(kind=realkind):: zfac
    real(kind=realkind)::     pten(klon,klev),        pqen(klon,klev), &
         plen(klon,klev),                              &
         puen(klon,klev),        pven(klon,klev),      &
         ptent(klon,klev),       ptenq(klon,klev),     &
         ptenl(klon,klev),                             &
         ptenu(klon,klev),       ptenv(klon,klev),     &
         pqsen(klon,klev),       pgeo(klon,klev),      &
         pap(klon,klev),         paph(klon,klev+1),    &
         pvervel(klon,klev),     pqhfl(klon),          &
         pahfs(klon)
    real(kind=realkind)::     ptu(klon,klev),         pqu(klon,klev),&
         plu(klon,klev),         plude(klon,klev),   &
         pld(klon,klev),         pleen(klon,klev),   &
         pmfu(klon,klev),        pmfd(klon,klev),    &
         penth(klon,klev),       prain(klon),        &
         pmflxr(klon,klev+1),    pmflxs(klon,klev+1),&
         pcond(klon,klev)
    integer::  kcbot(klon),            kctop(klon),     ktype(klon)
    logical::  ldland(klon),           ldcum(klon)
    !     
    real(kind=realkind)::     ztenh(klon,klev),       zqenh(klon,klev),&
         zlenh(klon,klev),                             &
         zgeoh(klon,klev),       zqsenh(klon,klev),    &
         ztd(klon,klev),         zqd(klon,klev),       &
         zmfus(klon,klev),       zmfds(klon,klev),     &
         zmfuq(klon,klev),       zmfdq(klon,klev),     &
         zdmfup(klon,klev),      zdmfdp(klon,klev),    &
         zmful(klon,klev),       zrfl(klon),           &
         zmfdl(klon,klev),                             &
         zuu(klon,klev),         zvu(klon,klev),       &
         zud(klon,klev),         zvd(klon,klev)
    real(kind=realkind)::     zentr(klon), zhcbase(klon), zmfub(klon),zmfub1(klon), &
         zdqpbl(klon),           zdqcv(klon)
    real(kind=realkind)::     zsfl(klon),             zdpmel(klon,klev)
    real(kind=realkind)::     zlglac(klon,klev)
    integer::  ilab(klon,klev),idtop(klon),ictop0(klon),ilwmin(klon)
    logical::  llddraf(klon)
    logical::  llo1

    !     1.           specify constants and parameters


    zcons2=1.0_realkind/(rg*ptsphy)

    !     *    2.           initialize values at vertical grid points in 'cuini'
    call cuini( kidia,    kfdia,    klon,     ktdia,    klev &
         , pten,     pqen,     pqsen,    puen,     pven       &
         , plen,     pvervel,  pgeo,     paph                 &
         , ilwmin,   ilab                                     &
         , ztenh,    zqenh,    zqsenh,   zgeoh,    zlenh      &
         , ptu,      pqu,      ztd,      zqd                  &
         , zuu,      zvu,      zud,      zvd                  &
         , plu,      pld       )

    !     *    3.0          cloud base calculations

    !     *             (a) determine cloud base values in 'cubase'
    call cubase( kidia,    kfdia,    klon,     ktdia,    klev &
         , ztenh,    zqenh,    zgeoh,    paph                  &
         , puen,     pven                                      &
         , ptu,      pqu,      plu,      zuu,      zvu         &
         , ilab,     ldcum,    kcbot  )
    !     
    !     *             (b) determine total moisture convergence and
    !     *                 then decide on type of cumulus convection
    !     -----------------------------------------
    !     
    jk=1
    do 310 jl=kidia,kfdia
       zdqcv(jl) =ptenq(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
       zdqpbl(jl)=0.0_realkind
       idtop(jl)=0
310 enddo
    do 320 jk=2,klev
       do 315 jl=kidia,kfdia
          zdqcv(jl)=zdqcv(jl)+ptenq(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
          if(jk>=kcbot(jl))then
             zdqpbl(jl)=zdqpbl(jl)+ptenq(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
          endif
315    enddo
320 enddo
    do 330 jl=kidia,kfdia
       if(zdqcv(jl)>max(0.0_realkind,-1.1_realkind*pqhfl(jl)*rg)) then
          ktype(jl)=1
       else
          ktype(jl)=2
       endif
330 enddo
    !     
    !     *             (c) determine moisture supply for boundary layer
    !     *                 and determine cloud base massflux ignoring
    !     *                 the effects of downdrafts at this stage
    !     ------------------------------------------
    !     
    do 340 jl=kidia,kfdia
       ikb=kcbot(jl)
       zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
       zdqmin=max(0.01_realkind*zqenh(jl,ikb),1.e-14_realkind)
       if(zdqpbl(jl)>0.0_realkind.and.zqumqe>zdqmin.and.ldcum(jl)) then
          zmfub(jl)=zdqpbl(jl)/(rg*max(zqumqe,zdqmin))
       else
          zmfub(jl)=0.01_realkind
          ldcum(jl)=.false.
       endif
       if(ktype(jl)==2) then
          zhfl=-pahfs(jl)-rlvtt*pqhfl(jl)
          zdh=rcpd*(ptu(jl,ikb)-ztenh(jl,ikb))+rlvtt*zqumqe
          zmf2=zhfl/max(zdh,1.e5_realkind*zdqmin)
          if (zmf2>0.0_realkind) then
             zmfub(jl)=zmf2
          endif
       end if
       zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2
       zmfub(jl)=min(zmfub(jl),zmfmax)
       if(ktype(jl)==1) then
          zentr(jl)=entrpen
       else
          zentr(jl)=entrscv
       endif
340 enddo
    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     *    4.0          determine cloud ascent for entraining plume
    !     -------------------------------------------
    !     
    !     

    !     
    !     *             (a) estimate cloud height for entrainment/detrainment
    !     *                 calculations in cuasc (max.possible cloud height
    !     *                 for non-entraining plume, following a.-s.,1974)
    !     -------------------------------------------------
    !     
    do 410 jl=kidia,kfdia
       ikb=kcbot(jl)
       zhcbase(jl)=rcpd*ptu(jl,ikb)+zgeoh(jl,ikb)+rlvtt*pqu(jl,ikb)
       ictop0(jl)=kcbot(jl)-1
410 enddo
    do 430 jk=klev-1,3,-1
       do 420 jl=kidia,kfdia
          zal=foelhm(ztenh(jl,jk))
          zqal=1.0_realkind/zal
          zalfaw=foealfa(ztenh(jl,jk))
          zalfai=1.0_realkind-zalfaw
          z5ldcp=zalfaw*r5alvcp+zalfai*r5alscp
          z4es=zalfaw*r4les+zalfai*r4ies
          zhsat=rcpd*ztenh(jl,jk)+zgeoh(jl,jk)+zal*zqsenh(jl,jk)
          zgam=z5ldcp*zqsenh(jl,jk)/((1.0_realkind-retv*zqsenh(jl,jk))* &
               (ztenh(jl,jk)-z4es)**2.0_realkind)
          zzz=rcpd*ztenh(jl,jk)*0.608_realkind
          zhhat=zhsat-(zzz+zgam*zzz)/(1.0_realkind+zgam*zzz*zqal)* &
               max(zqsenh(jl,jk)-zqenh(jl,jk),0.0_realkind)
          if(jk<ictop0(jl).and.zhcbase(jl)>zhhat) ictop0(jl)=jk
420    enddo
430 enddo
    !     
    !     *             (b) do ascent in 'cuasc'in absence of downdrafts
    !     --------------------------------------------
    !     
    call cuasc(kidia,    kfdia,    klon,     ktdia,    klev &
         , ptsphy                                            &
         , ztenh,    zqenh,    zlenh,    puen,     pven      &
         , pten,     pqen,     pqsen,    plen                &
         , pgeo,     zgeoh,    pap,      paph                &
         , ptenq,    pvervel,  ilwmin                        &
         , ldland,   ldcum,    ktype,    ilab                &
         , ptu,      pqu,      plu,      zuu,      zvu       &
         , pmfu,     zmfub,    zentr                         &
         , zmfus,    zmfuq,    zmful,    plude,    pleen     &
         , zdmfup,   zlglac                                  &
         , kcbot,    kctop,    ictop0 )
    !     
    !     
    !     *         (c) check cloud depth and change entrainment rate accordingly
    !     calculate precipitation rate (for downdraft calculation)
    !     -----------------------------------------------------
    !     
    !     
    do 440 jl=kidia,kfdia
       zpbmpt=paph(jl,kcbot(jl))-paph(jl,kctop(jl))
       if(ldcum(jl).and.ktype(jl)==1.and.zpbmpt<2.e4_realkind) ktype(jl)=2
       if(ldcum(jl)) ictop0(jl)=kctop(jl)
       if(ktype(jl)==2) zentr(jl)=entrscv
       zrfl(jl)=zdmfup(jl,1)
440 enddo
    do 460 jk=2,klev
       do 450 jl=kidia,kfdia
          zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
450    enddo
460 enddo
    do 480 jk=1,klev
       do 470 jl=kidia,kfdia
          pmfd(jl,jk)=0.0_realkind
          zmfds(jl,jk)=0.0_realkind
          zmfdq(jl,jk)=0.0_realkind
          zmfdl(jl,jk)=0.0_realkind
          zdmfdp(jl,jk)=0.0_realkind
          zdpmel(jl,jk)=0.0_realkind
470    enddo
480 enddo
    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     *    5.0          cumulus downdraft calculations
    !     ------------------------------
    !     
    !     

    !     
    if(lmfdd) then
       !     
       !     *             (a) determine lfs in 'cudlfs'
       !     -------------------------
       !     
       call cudlfs( kidia,    kfdia,    klon,     ktdia,    klev  &
            , kcbot,    kctop,    ldland,   ldcum                  &
            , ztenh,    zqenh,    zlenh,    puen,     pven         &
            , zgeoh,    paph,     ptu,      pqu,      plu          &
            , zuu,      zvu,      zmfub,    zrfl                   &
            , ztd,      zqd,      pld,      zud,      zvd          &
            , pmfd,     zmfds,    zmfdq,    zmfdl,    zdmfdp       &
            , idtop,    llddraf )
       !     
       !     *            (b)  determine downdraft t,q,l and fluxes in 'cuddraf'
       !     -----------------------------------------------
       !     
       call cuddraf( kidia,    kfdia,    klon,     ktdia,    klev &
            , llddraf                                              &
            , ztenh,    zqenh,    zlenh,    puen,     pven         &
            , zgeoh,    paph,     zrfl                             &
            , ztd,      zqd,      pld,      zud,      zvd          &
            , pmfd,     zmfds,    zmfdq,    zmfdl,    zdmfdp )
       !     
       !     *            (c)  recalculate convective fluxes due to effect of
       !     downdrafts on boundary layer moisture budget
       !     --------------------------------------------
       !     
       do 520 jl=kidia,kfdia
          if(llddraf(jl)) then
             ikb=kcbot(jl)
             llo1=pmfd(jl,ikb)<0.0_realkind
             if(pmfd(jl,ikb)<0.0_realkind) then
                zeps=cmfdeps
             else
                zeps=0.0_realkind
             endif
             !     
             !     no interaction between cloud water of the
             !     environment and new cloud
             !     -------------------------   
             zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zeps*(zqd(jl,ikb)+pld(jl,ikb))- &
                  (1.0_realkind-zeps)*(zqenh(jl,ikb))
             !     
             zdqmin=max(0.01_realkind*zqenh(jl,ikb),1.e-14_realkind)
             zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2
             if(zdqpbl(jl)>0.0_realkind.and.zqumqe>zdqmin.and.ldcum(jl) &
                  .and.zmfub(jl)<zmfmax) then
                zmfub1(jl)=zdqpbl(jl)/(rg*max(zqumqe,zdqmin))
             else
                zmfub1(jl)=zmfub(jl)
             endif
             if(ktype(jl)==2) then
                zhfl=-pahfs(jl)-rlvtt*pqhfl(jl)
                zdh=rcpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb)- &
                     (1.0_realkind-zeps)*ztenh(jl,ikb))+rlvtt*zqumqe
                zmf2=zhfl/max(zdh,1.e5_realkind*zdqmin)
                if (zmf2>0.0_realkind) then
                   zmfub1(jl)=zmf2
                endif
             end if
             if(.not.((ktype(jl)==1.or.ktype(jl)==2).and. &
                  abs(zmfub1(jl)-zmfub(jl))<0.2_realkind*zmfub(jl))) then
                zmfub1(jl)=zmfub(jl)
             endif
          end if
520    enddo
       do 540 jk=1,klev
          do 530 jl=kidia,kfdia
             if(llddraf(jl)) then
                zfac=zmfub1(jl)/max(zmfub(jl),1.e-14_realkind)
                pmfd(jl,jk)=pmfd(jl,jk)*zfac
                zmfds(jl,jk)=zmfds(jl,jk)*zfac
                zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
                zmfdl(jl,jk)=zmfdl(jl,jk)*zfac
                zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
             end if
530       enddo
540    enddo
       do 550 jl=kidia,kfdia
          if(llddraf(jl)) zmfub(jl)=zmfub1(jl)
550    enddo
       !     
    end if
    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     *    6.0          determine final cloud ascent for entraining plume
    !     *                 for penetrative convection (type=1),
    !     *                 for shallow to medium convection (type=2)
    !     *                 and for mid-level convection (type=3).
    !     -------------------------------------------------
    !     

    call cuasc( kidia,    kfdia,    klon,     ktdia,    klev &
         , ptsphy                                             &
         , ztenh,    zqenh,    zlenh,    puen,     pven       &
         , pten,     pqen,     pqsen,    plen                 &
         , pgeo,     zgeoh,    pap,      paph                 &
         , ptenq,    pvervel,  ilwmin                         &
         , ldland,   ldcum,    ktype,    ilab                 &
         , ptu,      pqu,      plu,      zuu,      zvu        &
         , pmfu,     zmfub,    zentr                          &
         , zmfus,    zmfuq,    zmful,    plude,    pleen      &
         , zdmfup,   zlglac                                   &
         , kcbot,    kctop,    ictop0 )
    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     *    7.0          determine final convective fluxes in 'cuflx'
    !     ------------------------------------------
    !     

    call cuflx( kidia,    kfdia,    klon,     ktdia,    klev &
         , ptsphy                                             &
         , pten,     pqen,     plen                           &
         , ztenh,    zqenh,    zlenh                          &
         , paph,     zgeoh,    pap,      ldland,   ldcum      &
         , kcbot,    kctop,    idtop,    itopm2               &
         , ktype,    llddraf                                  &
         , pmfu,     pmfd,     zmfus,    zmfds                &
         , zmfuq,    zmfdq,    zmful,    zmfdl,    plude      &
         , pleen,    zdmfup,   zdmfdp,   pqsen,    zlglac     &
         , zdpmel,   pmflxr,   pmflxs,   prain )
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     *    8.0          update tendencies for t,q and l in subroutine cudtdq
    !     --------------------------------------------------
    !     

    call cudtdq( kidia,    kfdia,    klon,     ktdia,    klev      &
         , itopm2,   ldland,   ldcum,    ptsphy                     &
         , paph,     pten,     pmfu,     zlenh                      &
         , zmfus,    zmfds,    zmfuq,    zmfdq                      &
         , zmful,    zmfdl,    zdmfup,   zdmfdp,   zdpmel,   zlglac &
         , ptent,    ptenq,    ptenl,    penth,    pcond)
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     *    9.0          update tendencies for u and v in subroutine cududv
    !     --------------------------------------------------
    !     

    if(lmfdudv) then
       call cududv( kidia,    kfdia,    klon,     ktdia,    klev &
            , itopm2,   ktype,    kcbot,    ldcum,    paph        &
            , puen,     pven,     pmfu,     pmfd                  &
            , zuu,      zud,      zvu,      zvd                   &
            , ptenu,    ptenv     )
    end if
    !     

    !     
    !     
    return
  end subroutine cumastr

  subroutine cuasc( kidia,    kfdia,    klon,     ktdia,    klev &
       , ptsphy                                                   &
       , ptenh,    pqenh,    plenh,    puen,     pven             &
       , pten,     pqen,     pqsen,    plen                       &
       , pgeo,     pgeoh,    pap,      paph                       &
       , ptenq,    pvervel,  klwmin                               &
       , ldland,   ldcum,    ktype,    klab                       &
       , ptu,      pqu,      plu,      puu,      pvu              &
       , pmfu,     pmfub,    pentr                                &
       , pmfus,    pmfuq,    pmful,    plude,    pleen            &
       , pdmfup,   plglac                                         &
       , kcbot,    kctop,    kctop0 )
    !     
    !     this routine does the calculations for cloud ascents
    !     for cumulus parameterization
    !     
    !     m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89
    !     
    !     purpose.
    !     --------
    !     to produce cloud ascents for cu-parametrization
    !     (vertical profiles of t,q,l,u and v and corresponding
    !     fluxes as well as precipitation rates)
    !     
    !     interface
    !     ---------
    !     
    !     this routine is called from *cumastr*.
    !     
    !     method.
    !     --------
    !     lift surface air dry-adiabatically to cloud base
    !     and then calculate moist ascent for
    !     entraining/detraining plume.
    !     entrainment and detrainment rates differ for
    !     shallow and deep cumulus convection.
    !     in case there is no penetrative or shallow convection
    !     check for possibility of mid level convection
    !     (cloud base values calculated in *cubasmc*)
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        end point
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *klwmin*       level of maximum vertical velocity
    !     *ktype*        type of convection
    !     1 = penetrative convection
    !     2 = shallow convection
    !     3 = midlevel convection
    !     *kcbot*        cloud base level
    !     
    !     input parameters (real):
    !     
    !     *ptsphy*       time step for the physics                      s
    !     *ptenh*        env. temperature (t+1) on half levels          k
    !     *pqenh*        env. spec. humidity (t+1) on half levels     kg/kg
    !     *plenh*        env. cloud water c. (t+1) on half levels     kg/kg
    !     *puen*         provisional environment u-velocity (t+1)      m/s
    !     *pven*         provisional environment v-velocity (t+1)      m/s
    !     *pten*         provisional environment temperature (t+1)      k
    !     *pqen*         provisional environment spec. humidity (t+1) kg/kg
    !     *plen*         provisional environment cloud wat. c.  (t+1) kg/kg
    !     *pqsen*        environment spec. saturation humidity (t+1)  kg/kg
    !     *pgeo*         geopotential                                 m2/s2
    !     *pgeoh*        geopotential on half levels                  m2/s2
    !     *pap*          provisional pressure on full levels           pa
    !     *paph*         provisional pressure on half levels           pa
    !     *ptenq*        moisture tendency                            kg/(kg s)
    !     *pvervel*      vertical velocity                            pa/s
    !     
    !     input parameters (logical):
    !     
    !     *ldland*       land sea mask (.true. for land)
    !     *ldcum*        flag: .true. for convective points
    !     
    !     updated parameters (integer):
    !     
    !     *klab*         flag klab=1 for subcloud levels
    !     klab=2 for cloud levels
    !     
    !     updated parameters (real):
    !     
    !     *ptu*          temperature in updrafts                        k
    !     *pqu*          spec. humidity in updrafts                   kg/kg
    !     *plu*          liquid water content in updrafts             kg/kg
    !     *puu*          u-velocity in updrafts                        m/s
    !     *pvu*          v-velocity in updrafts                        m/s
    !     
    !     output parameters (integer):
    !     
    !     *kctop*        cloud top level
    !     *kctop0*       first guess of cloud top level
    !     
    !     output parameters (real):
    !     
    !     *pmfu*         massflux in updrafts                         kg/(m2*s)
    !     *pmfub*        massflux in updrafts at cloud base           kg/(m2*s)
    !     *pentr*        fractional mass entrainment rate              1/m
    !     *pmfus*        flux of dry static energy in updrafts         j/(m2*s)
    !     *pmfuq*        flux of spec. humidity in updrafts           kg/(m2*s)
    !     *pmful*        flux of liquid water in updrafts             kg/(m2*s)
    !     *plude*        detrained liquid water                       kg/(m3*s)
    !     *pleen*        entrained liquid water                       kg/(m3*s)
    !     *pdmfup*       flux difference of precip. in updrafts       kg/(m2*s)
    !     *plglac*       frozen cloud water content                    kg/kg
    !     
    !     
    !     externals
    !     ---------
    !     *cuadjtq* adjust t and q due to condensation in ascent
    !     *cuentr*  calculate entrainment/detrainment rates
    !     *cubasmc* calculate cloud base values for midlevel convection
    !     
    !     reference
    !     ---------
    !     (tiedtke,1989)
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)
    !     

    implicit none

    real(kind=realkind):: zcons2,ptsphy,ztglace,zdphi,zmfmax,zfac,zmftest,zqeen
    real(kind=realkind):: zseen,zscde,zqude,zmfusk,zmfuqk,zmfulk,zalfawu,zalfaiu
    real(kind=realkind):: zalfawd,zalfaid,zbuo,zdnoprc,zprcon,zlnew,zz,zdmfeu
    real(kind=realkind):: zdmfdu,zzdmf
    integer:: jl,ktdia,kfdia,jk,ik,is,kidia,icall,kcum,klon,klev

    real(kind=realkind)::     ptenh(klon,klev),       pqenh(klon,klev), &
         plenh(klon,klev),                              &
         puen(klon,klev),        pven(klon,klev),       &
         pten(klon,klev),        pqen(klon,klev),       &
         plen(klon,klev),                               &
         pgeo(klon,klev),        pgeoh(klon,klev),      &
         pap(klon,klev),         paph(klon,klev+1),     &
         pqsen(klon,klev),       ptenq(klon,klev),      &
         pvervel(klon,klev)
    !     
    real(kind=realkind)::     ptu(klon,klev),         pqu(klon,klev), &
         puu(klon,klev),         pvu(klon,klev),      &
         pmfu(klon,klev),                             &
         pmfub(klon),            pentr(klon),         &
         pmfus(klon,klev),       pmfuq(klon,klev),    &
         plu(klon,klev),         plude(klon,klev),    &
         pleen(klon,klev),                            &
         pmful(klon,klev),       pdmfup(klon,klev),   &
         plglac(klon,klev)
    integer::  klwmin(klon),           ktype(klon), &
         klab(klon,klev),        kcbot(klon),      &
         kctop(klon),            kctop0(klon)
    logical::  ldland(klon),           ldcum(klon)
    !     
    real(kind=realkind)::     zdmfen(klon),           zdmfde(klon), &
         zmfuu(klon),            zmfuv(klon),       &
         zpbase(klon),           zqold(klon)
    real(kind=realkind)::     zdland(klon)
    real(kind=realkind)::     zph(klon)
    logical::  llflag(klon), llsum,  llo1

    !     *    1.           specify parameters

    zcons2=1.0_realkind/(rg*ptsphy)
    ztglace=rtt  -13.0_realkind
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     2.           set default values
    !     ------------------
    !     
    llsum=.false.
    llo1=.false.
    do 210 jl=kidia,kfdia
       zmfuu(jl)=0.0_realkind
       zmfuv(jl)=0.0_realkind
       if(.not.ldcum(jl)) ktype(jl)=0
210 enddo
    do 230 jk=1,klev
       do 220 jl=kidia,kfdia
          plu(jl,jk)=0.0_realkind
          pmfu(jl,jk)=0.0_realkind
          pmfus(jl,jk)=0.0_realkind
          pmfuq(jl,jk)=0.0_realkind
          pmful(jl,jk)=0.0_realkind
          plude(jl,jk)=0.0_realkind
          pleen(jl,jk)=0.0_realkind
          pdmfup(jl,jk)=0.0_realkind
          plglac(jl,jk)=0.0_realkind
          if(.not.ldcum(jl).or.ktype(jl)==3) klab(jl,jk)=0
          if(.not.ldcum(jl).and.paph(jl,jk)<4.e4_realkind) kctop0(jl)=jk
220    enddo
230 enddo
    do 240 jl=kidia,kfdia
       if(ldland(jl)) then
          zdland(jl)=3.0e4_realkind
          zdphi=pgeoh(jl,kctop0(jl))-pgeoh(jl,kcbot(jl))
          if(ptu(jl,kctop0(jl))>=ztglace) zdland(jl)=zdphi
          zdland(jl)=max(3.0e4_realkind,zdland(jl))
          zdland(jl)=min(5.0e4_realkind,zdland(jl))
       endif
240 enddo
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     3.0          initialize values at lifting level
    !     ----------------------------------
    !     

    do 310 jl=kidia,kfdia
       kctop(jl)=klev-1
       if(.not.ldcum(jl)) then
          kcbot(jl)=klev-1
          pmfub(jl)=0.0_realkind
          pqu(jl,klev)=0.0_realkind
       endif
       pmfu(jl,klev)=pmfub(jl)
       pmfus(jl,klev)=pmfub(jl)*(rcpd*ptu(jl,klev)+pgeoh(jl,klev))
       pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
       pmful(jl,klev)=pmfub(jl)*plu(jl,klev)
       if(lmfdudv) then
          zmfuu(jl)=pmfub(jl)*puu(jl,klev)
          zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
       endif
310 enddo
    do 320 jl=kidia,kfdia
       ldcum(jl)=.false.
320 enddo
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     4.           do ascent: subcloud layer (klab=1) ,clouds (klab=2)
    !     by doing first dry-adiabatic ascent and then
    !     by adjusting t,q and l accordingly in *cuadjtq*,
    !     then check for buoyancy and set flags accordingly
    !     -------------------------------------------------
    !     

    do 480 jk=klev-1,3,-1
       !     
       !     specify cloud base values for midlevel convection
       !     in *cubasmc* in case there is not already convection
       !     ----------------------------------------------------
       !     
       ik=jk
       call cubasmc( kidia,    kfdia,    klon,     ktdia,    klev    &
            , ik                                                      &
            , pten,     pqen,     pqsen,    plen,     puen,     pven  &
            , pvervel,  pgeo,     pgeoh,    ldcum,    ktype,    klab  &
            , kcbot,    pmfu,     pmfub,    pentr                     &
            , ptu,      pqu,      plu,      puu,      pvu             &
            , pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    zmfuv)
       !     
       is=0
       do 410 jl=kidia,kfdia
          is=is+klab(jl,jk+1)
          if(klab(jl,jk+1)==0) klab(jl,jk)=0
          if(klab(jl,jk+1)>0) then
             llflag(jl)=.true.
          else
             llflag(jl)=.false.
          endif
          zph(jl)=paph(jl,jk)
          if(ktype(jl)==3.and.jk==kcbot(jl)) then
             zmfmax=(paph(jl,jk)-paph(jl,jk-1))*zcons2
             if(pmfub(jl)>zmfmax) then
                zfac=zmfmax/pmfub(jl)
                pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
                pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
                pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
                pmful(jl,jk+1)=pmful(jl,jk+1)*zfac
                zmfuu(jl)=zmfuu(jl)*zfac
                zmfuv(jl)=zmfuv(jl)*zfac
                pmfub(jl)=zmfmax
             endif
          endif
410    enddo

       if(is>0) llo1=.true.


       !     *                  specify entrainment rates in *cuentr*
       ik=jk
       call cuentr( kidia,    kfdia,    klon,     ktdia,    klev   &
            , ik,       klwmin,   ktype,    kcbot,    kctop0        &
            , ldcum,    llo1                                        &
            , ptenh,    pqenh,    plenh,    ptenq,    paph,     pap &
            , pmfu,     pentr                                       &
            , zpbase,   zdmfen,   zdmfde )
       !     
       !     
       !     
       !     do adiabatic ascent for entraining/detraining plume
       !     ---------------------------------------------------
       !     
       if(llo1) then
          !     
          do 420 jl=kidia,kfdia
             if(llflag(jl)) then
                if(jk<kcbot(jl)) then
                   zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
                   zmfmax=min(zmftest,(paph(jl,jk)-paph(jl,jk-1))*zcons2)
                   zdmfen(jl)=max(zdmfen(jl)-max(zmftest-zmfmax,0.0_realkind),0.0_realkind)
                endif
                zdmfde(jl)=min(zdmfde(jl),0.75_realkind*pmfu(jl,jk+1))
                pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
                zqeen=pqenh(jl,jk+1)*zdmfen(jl)
                zseen=(rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
                zscde=(rcpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
                zqude=pqu(jl,jk+1)*zdmfde(jl)
                pleen(jl,jk)=plenh(jl,jk+1)*zdmfen(jl)
                plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
                zmfusk=pmfus(jl,jk+1)+zseen-zscde
                zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
                !     
                !     entrainment of cloud water is not allowed
                !     -----------------------------------------
                !     
                !     zmfulk=pmful(jl,jk+1)+pleen(jl,jk)-plude(jl,jk)
                zmfulk=pmful(jl,jk+1)             -plude(jl,jk)
                !     
                plu(jl,jk)=zmfulk*(1.0_realkind/max(cmfcmin,pmfu(jl,jk)))
                pqu(jl,jk)=zmfuqk*(1.0_realkind/max(cmfcmin,pmfu(jl,jk)))
                ptu(jl,jk)=(zmfusk*(1.0_realkind/max(cmfcmin,pmfu(jl,jk)))-pgeoh(jl,jk))/rcpd
                ptu(jl,jk)=max(100.0_realkind,ptu(jl,jk))
                ptu(jl,jk)=min(400.0_realkind,ptu(jl,jk))
                zqold(jl)=pqu(jl,jk)
             else
                zqold(jl)=0.0_realkind
             endif
420       enddo
          !     
          !     
          !     do corrections for moist ascent
          !     by adjusting t,q and l in *cuadjtq*
          !     -----------------------------------
          !     
          ik=jk
          icall=1
          call cuadjtq( kidia,    kfdia,    klon,     ktdia,    klev &
               , ik, zph,      ptu,      pqu,      llflag,  icall )
          !     
          do 440 jl=kidia,kfdia
             if(llflag(jl).and.abs(pqu(jl,jk)-zqold(jl))>1.e-14_realkind) then
                zalfawu=foealfa(ptu(jl,jk))
                zalfaiu=1.0_realkind-zalfawu
                zalfawd=foealfa(ptu(jl,jk+1))
                zalfaid=1.0_realkind-zalfawd
                plglac(jl,jk)=plu(jl,jk)*(zalfaiu-zalfaid)
                ptu(jl,jk)=ptu(jl,jk)+ralfdcp*plglac(jl,jk)
             endif
440       enddo
          do 442 jl=kidia,kfdia
             if(llflag(jl).and.abs(pqu(jl,jk)-zqold(jl))>1.e-14_realkind) then
                klab(jl,jk)=2
                plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
                zbuo=ptu(jl,jk)*(1.0_realkind+retv  *pqu(jl,jk))- &
                     ptenh(jl,jk)*(1.0_realkind+retv  *pqenh(jl,jk))
                if(klab(jl,jk+1)==1) zbuo=zbuo+0.5_realkind
                if(zbuo>0.0_realkind.and.pmfu(jl,jk)>=0.1_realkind*pmfub(jl)) then
                   kctop(jl)=jk
                   ldcum(jl)=.true.
                   if(ldland(jl)) then
                      zdnoprc=zdland(jl)
                   else
                      zdnoprc=1.5e4_realkind
                   endif
                   if(zpbase(jl)-paph(jl,jk)<zdnoprc) then
                      zprcon=0.0_realkind
                   else
                      zprcon=cprcon
                   endif
                   zlnew=plu(jl,jk)/(1.0_realkind+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
                   pdmfup(jl,jk)=max(0.0_realkind,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
                   plu(jl,jk)=zlnew
                else
                   klab(jl,jk)=0
                   pmfu(jl,jk)=0.0_realkind
                endif
             endif
442       enddo
          do 455 jl=kidia,kfdia
             if(llflag(jl)) then
                pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
                pmfus(jl,jk)=(rcpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
                pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
             endif
455       enddo
          if(lmfdudv) then
             do 460 jl=kidia,kfdia
                if(llflag(jl)) then
                   if(ktype(jl)==1.or.ktype(jl)==3) then
                      if(abs(zdmfen(jl))<1.e-14_realkind) then
                         zz=3.0_realkind
                      else
                         zz=2.0_realkind
                      endif
                   else
                      if(abs(zdmfen(jl))<1.e-14_realkind) then
                         zz=1.0_realkind
                      else
                         zz=0.0_realkind
                      endif
                   endif
                   zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
                   zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
                   zdmfdu=min(zdmfdu,0.75_realkind*pmfu(jl,jk+1))
                   zmfuu(jl)=zmfuu(jl)+zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
                   zmfuv(jl)=zmfuv(jl)+zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
                   if(pmfu(jl,jk)>0.0_realkind) then
                      puu(jl,jk)=zmfuu(jl)*(1.0_realkind/pmfu(jl,jk))
                      pvu(jl,jk)=zmfuv(jl)*(1.0_realkind/pmfu(jl,jk))
                   endif
                endif
460          enddo
          endif
          !     
       endif
480 enddo
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     5.           determine convective fluxes above non-buoyancy level
    !     ----------------------------------------------------
    !     (note: cloud variables like t,q and l are not
    !     affected by detrainment and are already known
    !     from previous calculations above)
    !     


    do 510 jl=kidia,kfdia
       if(kctop(jl)==klev-1) ldcum(jl)=.false.
       kcbot(jl)=max(kcbot(jl),kctop(jl))
510 enddo
    is=0
    do 520 jl=kidia,kfdia
       if(ldcum(jl)) then
          is=is+1
       endif
520 enddo
    if(is>0) then
       llsum=.true.
       kcum=is
    endif

    if(llsum) then
       do 530 jl=kidia,kfdia
          if(ldcum(jl)) then
             jk=kctop(jl)-1
             zzdmf=cmfctop
             zdmfde(jl)=(1.0_realkind-zzdmf)*pmfu(jl,jk+1)
             plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
             pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
             zlnew=plu(jl,jk)
             pdmfup(jl,jk)=max(0.0_realkind,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
             plu(jl,jk)=zlnew
             pmfus(jl,jk)=(rcpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
             pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
             pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
             plude(jl,jk-1)=pmful(jl,jk)
          endif
530    enddo
       if(lmfdudv) then

          do 540 jl=kidia,kfdia
             if(ldcum(jl)) then
                jk=kctop(jl)-1
                puu(jl,jk)=puu(jl,jk+1)
                pvu(jl,jk)=pvu(jl,jk+1)
             endif
540       enddo
       endif
    endif

    return
  end subroutine cuasc
  subroutine cubase( kidia,    kfdia,    klon,     ktdia,    klev &
       , ptenh,    pqenh,    pgeoh,    paph                          &
       , puen,     pven                                              &
       , ptu,      pqu,      plu,      puu,      pvu                 &
       , klab,     ldcum,    kcbot   )
    !
    !          this routine calculates cloud base values (t and q)
    !          for cumulus parameterization
    !
    !          m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89
    !
    !          purpose.
    !          --------
    !          to produce cloud base values for cu-parametrization
    !
    !          interface
    !          ---------
    !          this routine is called from *cumastr*.
    !          input are environm. values of t,q,p,phi at half levels.
    !          it returns cloud base values and flags as follows;
    !                 klab=1 for subcloud levels
    !                 klab=2 for condensation level
    !
    !          method.
    !          --------
    !          lift surface air dry-adiabatically to cloud base
    !          (non entraining plume,i.e.constant massflux)
    !
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !
    !    *kidia*        start point
    !    *kfdia*        end point
    !    *klon*         number of grid points per packet
    !    *ktdia*        start of the vertical loop
    !    *klev*         number of levels
    !
    !    input parameters (real):
    !
    !    *ptenh*        env. temperature (t+1) on half levels           k
    !    *pqenh*        env. spec. humidity (t+1) on half levels      kg/kg
    !    *pgeoh*        geopotential on half levels                   m2/s2
    !    *paph*         provisional pressure on half levels             pa
    !    *puen*         provisional environment u-velocity (t+1)       m/s
    !    *pven*         provisional environment v-velocity (t+1)       m/s
    !
    !    updated parameters (real):
    !
    !    *ptu*          temperature in updrafts                         k
    !    *pqu*          spec. humidity in updrafts                    kg/kg
    !    *plu*          liquid water content in updrafts              kg/kg
    !    *puu*          u-velocity in updrafts                         m/s
    !    *pvu*          v-velocity in updrafts                         m/s
    !
    !    updated parameters (integer):
    !
    !    *klab*         flag klab=1 for subcloud levels
    !                        klab=2 for cloud levels
    !
    !    output parameters (logical):
    !
    !    *ldcum*        flag: .true. for convective points
    !
    !    output parameters (integer):
    !
    !    *kcbot*        cloud base level
    !
    !
    !          externals
    !          ---------
    !          *cuadjtq* for adjusting t and q due to condensation in ascent
    !
    !          modifications
    !          -------------
    !             92-09-21 : update to cy44      j.-j. morcrette
    !        adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !        mod. for cloud water of envt: 97-04-23 (j. a. garcia-moya, inm)



    implicit none

    integer:: jl,kidia,kfdia,jk,is,ik,icall,ktdia,klon,klev,ikb
    real(kind=realkind):: zbuo,zz



    !
    real(kind=realkind):: ptenh(klon,klev), pqenh(klon,klev),pgeoh(klon,klev),paph(klon,klev+1)
    !
    real(kind=realkind):: ptu(klon,klev),pqu(klon,klev),plu(klon,klev)
    real(kind=realkind):: puen(klon,klev), pven(klon,klev),puu(klon,klev),pvu(klon,klev)
    integer::  klab(klon,klev),        kcbot(klon)
    logical::  ldcum(klon)
    !
    real(kind=realkind)::     zqold(klon)
    real(kind=realkind)::     zph(klon)
    logical::  llflag(klon)
    !
    !
    !
    !----------------------------------------------------------------------
    !
    !     1.           initialize values at lifting level
    !                  ----------------------------------
    !

    !
    do 110 jl=kidia,kfdia
       klab(jl,klev)=1
       kcbot(jl)=klev-1
       ldcum(jl)=.false.
       puu(jl,klev)=puen(jl,klev)*(paph(jl,klev+1)-paph(jl,klev))
       pvu(jl,klev)=pven(jl,klev)*(paph(jl,klev+1)-paph(jl,klev))
110 enddo
    !
    !
    !----------------------------------------------------------------------
    !
    !     2.0          do ascent in subcloud layer,
    !                  check for existence of condensation level,
    !                  adjust t,q and l accordingly in *cuadjtq*,
    !                  check for buoyancy and set flags
    !                  -------------------------------------
    !

    do jk=klev-1,2,-1
       is=0
       do 210 jl=kidia,kfdia
          if(klab(jl,jk+1)==1) then
             is=is+1
             llflag(jl)=.true.
          else
             llflag(jl)=.false.
          endif
          zph(jl)=paph(jl,jk)
210    enddo
       if(is==0) goto 290
       do 220 jl=kidia,kfdia
          if(llflag(jl)) then
             pqu(jl,jk)=pqu(jl,jk+1)
             ptu(jl,jk)=(rcpd*ptu(jl,jk+1)+pgeoh(jl,jk+1)-pgeoh(jl,jk))/rcpd
             zbuo=ptu(jl,jk)*(1.0_realkind+retv  *pqu(jl,jk))- &
                  ptenh(jl,jk)*(1.0_realkind+retv  *pqenh(jl,jk))+0.5_realkind
             if(zbuo>0.0_realkind) klab(jl,jk)=1
             zqold(jl)=pqu(jl,jk)
          end if
220    enddo
       !
       ik=jk
       icall=1
       call cuadjtq( kidia,    kfdia,    klon,     ktdia,    klev &
            , ik , zph,      ptu,      pqu,      llflag,   icall)

       do 240 jl=kidia,kfdia
          if(llflag(jl).and.abs(pqu(jl,jk)-zqold(jl))>1.e-14_realkind) then
             klab(jl,jk)=2
             plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
             zbuo=ptu(jl,jk)*(1.0_realkind+retv  *pqu(jl,jk))- &
                  ptenh(jl,jk)*(1.0_realkind+retv  *pqenh(jl,jk))+0.5_realkind
             if(zbuo>0.0_realkind) then
                kcbot(jl)=jk
                ldcum(jl)=.true.
             end if
          end if
240    enddo
       !
       !             calculate averages of u and v for subcloud ara,.
       !             the values will be used to define cloud base values.
       !
       if(lmfdudv) then
          do 250 jl=kidia,kfdia
             if(jk>=kcbot(jl)) then
                puu(jl,klev)=puu(jl,klev)+puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
                pvu(jl,klev)=pvu(jl,klev)+pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
             end if
250       enddo
       end if
       !
290    continue
    enddo
    !
    if(lmfdudv) then
       do 310 jl=kidia,kfdia
          if(ldcum(jl)) then
             ikb=kcbot(jl)
             zz=1.0_realkind/(paph(jl,klev+1)-paph(jl,ikb))
             puu(jl,klev)=puu(jl,klev)*zz
             pvu(jl,klev)=pvu(jl,klev)*zz
          else
             puu(jl,klev)=puen(jl,klev-1)
             pvu(jl,klev)=pven(jl,klev-1)
          end if
310    enddo
    end if
    !
    return
  end subroutine cubase

  subroutine cuddraf( kidia,    kfdia,    klon,     ktdia,    klev &
       , lddraf                                                       &
       , ptenh,    pqenh,    plenh,    puen,     pven                 &
       , pgeoh,    paph,     prfl                                     &
       , ptd,      pqd,      pld,      pud,      pvd                  &
       , pmfd,     pmfds,    pmfdq,    pmfdl,    pdmfdp )
    !     
    !     this routine calculates cumulus downdraft descent
    !     
    !     m.tiedtke         e.c.m.w.f.    12/86 modif. 12/89
    !     
    !     purpose.
    !     --------
    !     to produce the vertical profiles for cumulus downdrafts
    !     (i.e. t,q,l,u and v and fluxes)
    !     
    !     interface
    !     ---------
    !     
    !     this routine is called from *cumastr*.
    !     input is t,q,p,phi,u,v at half levels.
    !     it returns fluxes of s,q and evaporation rate
    !     and u,v at levels where downdraft occurs
    !     
    !     method.
    !     --------
    !     calculate moist descent for entraining/detraining plume by
    !     a) moving air dry-adiabatically to next level below and
    !     b) correcting for evaporation to obtain saturated state.
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     
    !     input parameters (logical):
    !     
    !     *lddraf*       .true. if downdrafts exist
    !     
    !     input parameters (real):
    !     
    !     *ptenh*        env. temperature (t+1) on half levels          k
    !     *pqenh*        env. spec. humidity (t+1) on half levels     kg/kg
    !     *plenh*        env. cloud wat. co. (t+1) on half levels     kg/kg
    !     *puen*         provisional environment u-velocity (t+1)      m/s
    !     *pven*         provisional environment v-velocity (t+1)      m/s
    !     *pgeoh*        geopotential on half levels                  m2/s2
    !     *paph*         provisional pressure on half levels           pa
    !     
    !     updated parameters (real):
    !     
    !     *prfl*         precipitation rate                      kg/(m2*s)
    !     
    !     output parameters (real):
    !     
    !     *ptd*          temperature in downdrafts                      k
    !     *pqd*          spec. humidity in downdrafts                 kg/kg
    !     *pld*          cloud water content in downdrafts            kg/kg
    !     *pud*          u-velocity in downdrafts                      m/s
    !     *pvd*          v-velocity in downdrafts                      m/s
    !     *pmfd*         massflux in downdrafts                   kg/(m2*s)
    !     *pmfds*        flux of dry static energy in downdrafts   j/(m2*s)
    !     *pmfdq*        flux of spec. humidity in downdrafts     kg/(m2*s)
    !     *pmfdl*        flux of cloud wat. co. in downdrafts     kg/(m2*s)
    !     *pdmfdp*       flux difference of precip. in downdrafts kg/(m2*s)
    !     
    !     externals
    !     ---------
    !     *cuadjtq* for adjusting t and q due to evaporation in
    !     saturated descent
    !     
    !     reference
    !     ---------
    !     (tiedtke,1989)
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)



    implicit none
    integer:: jk,is,jl,kidia,kfdia,itopde,klon,klev,ktdia,icall,ik
    real(kind=realkind):: zentr,zmfdvk,zmfduk,zdmfdp,zbuo,zmfdlk,zmfdqk,zmfdsk,zldde
    real(kind=realkind):: zqdde,zsdde,zleen,zqeen,zseen



    !     
    !     -----------------------------------------------------------------
    !     
    real(kind=realkind)::     ptenh(klon,klev),       pqenh(klon,klev), &
         plenh(klon,klev),                          &
         puen(klon,klev),        pven(klon,klev),   &
         pgeoh(klon,klev),       paph(klon,klev+1)
    !     
    real(kind=realkind)::     ptd(klon,klev),         pqd(klon,klev),   &
         pld(klon,klev),                            &
         pud(klon,klev),         pvd(klon,klev),    &
         pmfd(klon,klev),        pmfds(klon,klev),  &
         pmfdq(klon,klev),       pdmfdp(klon,klev), &
         pmfdl(klon,klev),                          &
         prfl(klon)                                  
    integer::  ktype(klon),            kcbot(klon),       &
         kdtop(klon)                                 
    logical::  lddraf(klon)                                    
    !     
    real(kind=realkind)::     zdmfen(klon),           zdmfde(klon),      &
         zcond(klon)
    real(kind=realkind)::     zph(klon)
    logical::  llo2(klon)
    logical::  llo1
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     1.           calculate moist descent for cumulus downdraft by
    !     (a) calculating entrainment rates, assuming
    !     linear decrease of massflux in pbl
    !     (b) doing moist descent - evaporative cooling
    !     and moistening is calculated in *cuadjtq*
    !     (c) checking for negative buoyancy and
    !     specifying final t,q,u,v and downward fluxes
    !     -------------------------------------------------
    !     

    do  jk=3,klev
       is=0
       do 110 jl=kidia,kfdia
          zph(jl)=paph(jl,jk)
          llo2(jl)=lddraf(jl).and.pmfd(jl,jk-1)<0.0_realkind
          if(llo2(jl)) then
             is=is+1
          endif
110    enddo
       if(is==0) goto 180
       do 122 jl=kidia,kfdia
          if(llo2(jl)) then
             zentr=entrdd*pmfd(jl,jk-1)*rd*ptenh(jl,jk-1)/ &
                  (rg*paph(jl,jk-1))*(paph(jl,jk)-paph(jl,jk-1))
             zdmfen(jl)=zentr
             zdmfde(jl)=zentr
          endif
122    enddo
       itopde=klev-2
       if(jk>itopde) then
          do 124 jl=kidia,kfdia
             if(llo2(jl)) then
                zdmfen(jl)=0.0_realkind
                zdmfde(jl)=pmfd(jl,itopde)*(paph(jl,jk)-paph(jl,jk-1))/ &
                     (paph(jl,klev+1)-paph(jl,itopde))
             endif
124       enddo
       endif
       !     
       do 126 jl=kidia,kfdia
          if(llo2(jl)) then
             pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
             zseen=(rcpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
             zqeen=pqenh(jl,jk-1)*zdmfen(jl)
             zleen=plenh(jl,jk-1)*zdmfen(jl)
             zsdde=(rcpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
             zqdde=pqd(jl,jk-1)*zdmfde(jl)
             zldde=pld(jl,jk-1)*zdmfde(jl)
             zmfdsk=pmfds(jl,jk-1)+zseen-zsdde
             zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde
             zmfdlk=pmfdl(jl,jk-1)+zleen-zldde
             pqd(jl,jk)=zmfdqk*(1.0_realkind/min(-cmfcmin,pmfd(jl,jk)))
             pld(jl,jk)=zmfdlk*(1.0_realkind/min(-cmfcmin,pmfd(jl,jk)))
             ptd(jl,jk)=(zmfdsk*(1.0_realkind/min(-cmfcmin,pmfd(jl,jk)))-&
                  pgeoh(jl,jk))/rcpd
             ptd(jl,jk)=min(400.0_realkind,ptd(jl,jk))
             ptd(jl,jk)=max(100.0_realkind,ptd(jl,jk))
             zcond(jl)=pqd(jl,jk)
          endif
126    enddo
       !     
       ik=jk
       icall=2
       call cuadjtq( kidia,    kfdia,    klon,     ktdia,    klev &
            , ik , zph,      ptd,      pqd,      llo2,     icall )
       !     
       do 150 jl=kidia,kfdia
          if(llo2(jl)) then
             pld(jl,jk)=pld(jl,jk)+(zcond(jl)-pqd(jl,jk))
             zcond(jl)=zcond(jl)-pqd(jl,jk)
             zbuo=ptd(jl,jk)*(1.0_realkind+retv  *pqd(jl,jk))- &
                  ptenh(jl,jk)*(1.0_realkind+retv  *pqenh(jl,jk))
             if(zbuo>=0.0_realkind.or.prfl(jl)<=(pmfd(jl,jk)*zcond(jl))) then
                pmfd(jl,jk)=0.0_realkind
             endif
             pmfds(jl,jk)=(rcpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
             pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
             pmfdl(jl,jk)=pld(jl,jk)*pmfd(jl,jk)
             zdmfdp=-pmfd(jl,jk)*zcond(jl)
             pdmfdp(jl,jk-1)=zdmfdp
             prfl(jl)=prfl(jl)+zdmfdp
          endif
150    enddo
       if(lmfdudv) then
          do 160 jl=kidia,kfdia
             if(llo2(jl).and.pmfd(jl,jk)<0.0_realkind) then
                zmfduk=pmfd(jl,jk-1)*pud(jl,jk-1)+ &
                     zdmfen(jl)*puen(jl,jk-1)-zdmfde(jl)*pud(jl,jk-1)
                zmfdvk=pmfd(jl,jk-1)*pvd(jl,jk-1)+ &
                     zdmfen(jl)*pven(jl,jk-1)-zdmfde(jl)*pvd(jl,jk-1)
                pud(jl,jk)=zmfduk*(1.0_realkind/min(-cmfcmin,pmfd(jl,jk)))
                pvd(jl,jk)=zmfdvk*(1.0_realkind/min(-cmfcmin,pmfd(jl,jk)))
             endif
160       enddo
       endif
       !     
180    continue
    enddo
    return
  end subroutine cuddraf


  subroutine cudlfs(kidia,    kfdia,    klon,     ktdia,    klev, &
       kcbot,    kctop,    ldland,   ldcum,                        &
       ptenh,    pqenh,    plenh,    puen,     pven,               &
       pgeoh,    paph,     ptu,      pqu,      plu,                &
       puu,      pvu,      pmfub,    prfl,                         &
       ptd,      pqd,      pld,      pud,      pvd,                &
       pmfd,     pmfds,    pmfdq,    pmfdl,    pdmfdp,             &
       kdtop,    lddraf)
    !     
    !     this routine calculates level of free sinking for
    !     cumulus downdrafts and specifies t,q,u and v values
    !     
    !     m.tiedtke         e.c.m.w.f.    12/86 modif. 12/89
    !     
    !     purpose.
    !     --------
    !     to produce lfs-values for cumulus downdrafts
    !     for massflux cumulus parameterization
    !     
    !     interface
    !     ---------
    !     this routine is called from *cumastr*.
    !     input are environmental values of t,q,u,v,p,phi
    !     and updraft values t,q,u and v and also
    !     cloud base massflux and cu-precipitation rate.
    !     it returns t,q,u and v values and massflux at lfs.
    !     
    !     method.
    !     
    !     check for negative buoyancy of air of equal parts of
    !     moist environmental air and cloud air.
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        end point
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *kcbot*        cloud base level
    !     *kctop*        cloud top level
    !     
    !     input parameters (logical):
    !     
    !     *ldland*       land sea mask (.true. for land)
    !     *ldcum*        flag: .true. for convective points
    !     
    !     input parameters (real):
    !     
    !     *ptenh*        env. temperature (t+1) on half levels       k
    !     *pqenh*        env. spec. humidity (t+1) on half levels  kg/kg
    !     *plenh*        env. cloud wat. co. (t+1) on half levels  kg/kg
    !     *puen*         provisional environment u-velocity (t+1)   m/s
    !     *pven*         provisional environment v-velocity (t+1)   m/s
    !     *pgeoh*        geopotential on half levels               m2/s2
    !     *paph*         provisional pressure on half levels        pa
    !     *ptu*          temperature in updrafts                     k
    !     *pqu*          spec. humidity in updrafts                kg/kg
    !     *plu*          liquid water content in updrafts          kg/kg
    !     *puu*          u-velocity in updrafts                     m/s
    !     *pvu*          v-velocity in updrafts                     m/s
    !     *pmfub*        massflux in updrafts at cloud base      kg/(m2*s)
    !     
    !     updated parameters (real):
    !     
    !     *prfl*         precipitation rate                        kg/(m2*s)
    !     
    !     output parameters (real):
    !     
    !     *ptd*          temperature in downdrafts                      k
    !     *pqd*          spec. humidity in downdrafts              kg/kg
    !     *pld*          cloud water content in downdrafts         kg/kg
    !     *pud*          u-velocity in downdrafts                   m/s
    !     *pvd*          v-velocity in downdrafts                   m/s
    !     *pmfd*         massflux in downdrafts                    kg/(m2*s)
    !     *pmfds*        flux of dry static energy in downdrafts    j/(m2*s)
    !     *pmfdq*        flux of spec. humidity in downdrafts      kg/(m2*s)
    !     *pmfdl*        flux of cloud wat. co. in downdrafts      kg/(m2*s)
    !     *pdmfdp*       flux difference of precip. in downdrafts  kg/(m2*s)
    !     
    !     output parameters (integer):
    !     
    !     *kdtop*        top level of downdrafts
    !     
    !     output parameters (logical):
    !     
    !     *lddraf*       .true. if downdrafts exist
    !     
    !     externals
    !     ---------
    !     *cuadjtq* for calculating wet bulb t and q at lfs
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)



    implicit none
    integer:: klon,klev,jl,kidia,kfdia,ike,jk,is,ik,icall,ktdia
    real(kind=realkind):: zttest,zqtest,zltest,zbuo,zmftop



    real(kind=realkind)::     ptenh(klon,klev),       pqenh(klon,klev), &
         plenh(klon,klev),                              &
         puen(klon,klev),        pven(klon,klev),       &
         pgeoh(klon,klev),       paph(klon,klev+1),     &
         ptu(klon,klev),         pqu(klon,klev),        &
         puu(klon,klev),         pvu(klon,klev),        &
         plu(klon,klev),                                &
         pmfub(klon),            prfl(klon)
    !     
    real(kind=realkind)::     ptd(klon,klev),         pqd(klon,klev), &
         pld(klon,klev),                              &
         pud(klon,klev),         pvd(klon,klev),      &
         pmfd(klon,klev),        pmfds(klon,klev),    &
         pmfdl(klon,klev),                            &
         pmfdq(klon,klev),       pdmfdp(klon,klev)     
    integer::  kcbot(klon),            kctop(klon),     &
         kdtop(klon)                                   
    logical::  ldland(klon),           ldcum(klon),     &
         lddraf(klon)                                  
    !                                                       
    real(kind=realkind)::     ztenwb(klon,klev),      zqenwb(klon,klev),&
         zlenwb(klon,klev),      zqold(klon,klev),    &
         zcond(klon),            zph(klon)
    logical::  llo2(klon)
    !     
    !----------------------------------------------------------------------
    !     
    !     1.           set default values for downdrafts
    !     ---------------------------------
    !     

    do 110 jl=kidia,kfdia
       lddraf(jl)=.false.
       kdtop(jl)=klev+1
110 enddo
    !     
    if(lmfdd)then
       !     
       !     
       !----------------------------------------------------------------------
       !     
       !     2.           determine level of free sinking by
       !     doing a scan from top to base of cumulus clouds
       !     
       !     for every point and proceed as follows:
       !     
       !     (1) determine wet bulb environmental t and q
       !     (2) do mixing with cumulus cloud air
       !     (3) check for negative buoyancy
       !     
       !     the assumption is that air of downdrafts is mixture
       !     of 50% cloud air + 50% environmental air at wet bulb
       !     temperature (i.e. which became saturated due to
       !     evaporation of rain and cloud water)
       !     ----------------------------------------------------
       !     

       !     
       ike=klev-3
       do 290 jk=3,ike
          !     
          !     
          !     2.1          calculate wet-bulb temperature and moisture
          !     for environmental air in *cuadjtq*
          !     -------------------------------------------
          !     
          is=0
          do 212 jl=kidia,kfdia
             ztenwb(jl,jk)=ptenh(jl,jk)
             zqenwb(jl,jk)=pqenh(jl,jk)
             zlenwb(jl,jk)=plenh(jl,jk)
             zph(jl)=paph(jl,jk)
             llo2(jl)=ldcum(jl).and.prfl(jl)>0.0_realkind.and..not.lddraf(jl).and. &
                  (jk<kcbot(jl).and.jk>kctop(jl))
             if(llo2(jl))then
                is=is+1
             endif
212       enddo
          if(is/=0)then
             !     
             do 214 jl=kidia,kfdia
                zqold(jl,jk)=zqenwb(jl,jk)
214          enddo
             !     
             ik=jk
             icall=2
             call cuadjtq( kidia,    kfdia,    klon,     ktdia,    klev &
                  , ik, zph,      ztenwb,   zqenwb,   llo2,     icall)
             !     
             !     
             !     2.2          do mixing of cumulus and environmental air
             !     and check for negative buoyancy.
             !     then set values for downdraft at lfs.
             !     ----------------------------------------
             !     

             do 222 jl=kidia,kfdia
                if(llo2(jl)) then
                   zlenwb(jl,jk)=zlenwb(jl,jk)+(zqold(jl,jk)-zqenwb(jl,jk))
                   zttest=0.5_realkind*(ptu(jl,jk)+ztenwb(jl,jk))
                   zqtest=0.5_realkind*(pqu(jl,jk)+zqenwb(jl,jk))
                   !     
                   !     environmental cloud water doesn't 
                   !     influence new convection 
                   !     ------------------------ 
                   !     
                   !     zltest=0.5*(plu(jl,jk)+zlenwb(jl,jk))
                   zltest=0.0_realkind
                   !     
                   zbuo=zttest*(1.0_realkind+retv  *zqtest)-ptenh(jl,jk)*&
                        (1.0_realkind+retv  *pqenh(jl,jk))
                   zcond(jl)=zltest+pqenh(jl,jk)-zqenwb(jl,jk)
                   zmftop=-cmfdeps*pmfub(jl)
                   if(zbuo<0.0_realkind.and.prfl(jl)>10.0_realkind*zmftop*zcond(jl)) then
                      kdtop(jl)=jk
                      lddraf(jl)=.true.
                      ptd(jl,jk)=zttest
                      pqd(jl,jk)=zqtest
                      pld(jl,jk)=zltest
                      pmfd(jl,jk)=zmftop
                      pmfds(jl,jk)=pmfd(jl,jk)*(rcpd*ptd(jl,jk)+pgeoh(jl,jk))
                      pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)
                      pmfdl(jl,jk)=pmfd(jl,jk)*pld(jl,jk)
                      pdmfdp(jl,jk-1)=-0.5_realkind*pmfd(jl,jk)*zcond(jl)
                      prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
                   end if
                end if
222          enddo
             if(lmfdudv) then
                do 224 jl=kidia,kfdia
                   if(pmfd(jl,jk)<0.0_realkind) then
                      pud(jl,jk)=0.5_realkind*(puu(jl,jk)+puen(jl,jk-1))
                      pvd(jl,jk)=0.5_realkind*(pvu(jl,jk)+pven(jl,jk-1))
                   end if
224             enddo
             end if
             !     
          endif
290    enddo
       !     
    endif
    !     
    return
  end subroutine cudlfs


  subroutine cudtdq(  kidia,    kfdia,    klon,     ktdia,    klev     &
       ,  ktopm2,   ldland,   ldcum,    ptsphy                    &
       ,  paph,     pten,     pmfu,     plenh                     &
       ,  pmfus,    pmfds,    pmfuq,    pmfdq                     &
       ,  pmful,    pmfdl,    pdmfup,   pdmfdp,   pdpmel,   plglac &
       ,  ptent,    ptenq,    ptenl,    penth ,   pcond)
    !     
    !     
    !cudtdq - updates t and q tendencies, precipitation rates
    !     does global diagnostics
    !     
    !     m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89
    !     
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)
    !     
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     
    !     input parameters (logical):
    !     
    !     *ldland*       land sea mask (.true. for land)
    !     *ldcum*        flag: .true. for convective points
    !     
    !     input parameters (real):
    !     
    !     *ptsphy*       time step for the physics                       s
    !     *paph*         provisional pressure on half levels            pa
    !     *pten*         provisional environment temperature (t+1)       k
    !     *pmfu*         massflux updrafts                         kg/(m2*s)
    !     *plenh*        env. cloud wat. co. (t+1) on half levels  kg/kg
    !     *pmfus*        flux of dry static energy in updrafts      j/(m2*s)
    !     *pmfds*        flux of dry static energy in downdrafts    j/(m2*s)
    !     *pmfuq*        flux of spec. humidity in updrafts        kg/(m2*s)
    !     *pmfdq*        flux of spec. humidity in downdrafts      kg/(m2*s)
    !     *pmful*        flux of liquid water in updrafts          kg/(m2*s)
    !     *pmfdl*        flux of liquid water in downdrafts        kg/(m2*s)
    !     *pdmfup*       flux difference of precip. in updrafts    kg/(m2*s)
    !     *pdmfdp*       flux difference of precip. in downdrafts  kg/(m2*s)
    !     *pdpmel*       change in precip.-fluxes due to melting   kg/(m2*s)
    !     *plglac*       flux of frozen cloud water                kg/(m2*s)
    !     
    !     updated parameters (real):
    !     
    !     *ptent*        temperature tendency                           k/s
    !     *ptenq*        moisture tendency                         kg/(kg s)
    !     *ptenl*        cloud water content tendency              kg/(kg s)
    !     
    !     output parameters (real):
    !     
    !     *penth*        increment of dry static energy             j/(kg*s)
    !     *pcond*        condensate                                     k/s
    !     


    implicit none

    integer:: jk,jl,kidia,kfdia,klon,ktdia,klev,ktopm2
    real(kind=realkind):: zalv,zdtdt,zdqdt,zdldt,ptsphy

    logical::  llo1,llo2,llo3

    real(kind=realkind)::     ptent(klon,klev),       ptenq(klon,klev), &
         ptenl(klon,klev),       pcond(klon,klev),      &
         pten(klon,klev),                               &
         pgeo(klon,klev),        paph(klon,klev+1)       
    real(kind=realkind)::     pmfus(klon,klev),       pmfds(klon,klev),  &
         pmfuq(klon,klev),       pmfdq(klon,klev),      &
         pmful(klon,klev),       pmfdl(klon,klev),      &
         penth(klon,klev),                              &
         pdmfup(klon,klev),      pdmfdp(klon,klev),     &
         pdpmel(klon,klev),      plglac(klon,klev),     &
         pmfu(klon,klev),        plenh(klon,klev)
    logical::  ldland(klon),           ldcum(klon)
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     
    !     *    1.0          incrementation of t and q tendencies
    !     ------------------------------------
    !     


    do 120 jk=1,klev
       do 110 jl=kidia,kfdia
          penth(jl,jk)=0.0_realkind
110    enddo
120 enddo
    !     
    do 150 jk=1,klev
       !     
       if(jk<klev) then
          do 130 jl=kidia,kfdia
             if(ldcum(jl)) then
                zalv=foelhm(pten(jl,jk))
                zdtdt=(rg/(paph(jl,jk+1)-paph(jl,jk)))/rcpd* &
                     (pmfus(jl,jk+1)-pmfus(jl,jk)+            &
                     pmfds(jl,jk+1)-pmfds(jl,jk)              &
                     +rlmlt*plglac(jl,jk)                     &
                     -rlmlt*pdpmel(jl,jk)                     &
                     -zalv*(pmful(jl,jk+1)-pmful(jl,jk)-      &
                     (pdmfup(jl,jk)+pdmfdp(jl,jk))))
                !     
                !     increment of dry static energy
                !     
                penth(jl,jk)= zdtdt*rcpd
                !     
                pcond(jl,jk)=(rg/(paph(jl,jk+1)-paph(jl,jk)))/rcpd* &
                     (rlmlt*plglac(jl,jk)                            &
                     -rlmlt*pdpmel(jl,jk)                            &
                     -zalv*(pmful(jl,jk+1)-pmful(jl,jk)-             &
                     (pdmfup(jl,jk)+pdmfdp(jl,jk))))
                !     
                ptent(jl,jk)=ptent(jl,jk)+zdtdt
                !     
                zdqdt=(rg/(paph(jl,jk+1)-paph(jl,jk)))* &
                     (pmfuq(jl,jk+1)-pmfuq(jl,jk)+       &
                     pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                     pmful(jl,jk+1)-pmful(jl,jk)-        &
                     (pdmfup(jl,jk)+pdmfdp(jl,jk)))
                ptenq(jl,jk)=ptenq(jl,jk)+zdqdt
                !     
                !     cloud water tendency is managed by sundqvist
                !     --------------------------------------------
                !     here is the right place to substract 
                !     environmental compensating subsidence from
                !     pmful (not done in cuflx)
                !     
                !     zdldt=(rg/(paph(jl,jk+1)-paph(jl,jk)))*
                !     &            ((pmful(jl,jk+1)-pmfu(jl,jk+1)*plenh(jl,jk+1))-
                !     &            (pmful(jl,jk)-pmfu(jl,jk)*plenh(jl,jk))+
                !     &            pmfdl(jl,jk+1)-pmfdl(jl,jk)+
                !     &            (pdmfup(jl,jk)+pdmfdp(jl,jk)))
                zdldt=0.0_realkind
                !     
                ptenl(jl,jk)=ptenl(jl,jk)+zdldt
             endif
125          format(i3,2x,7(e12.5,2x))
130       enddo
          !     
       else
          do 140 jl=kidia,kfdia
             if(ldcum(jl)) then
                zalv=foelhm(pten(jl,jk))
                zdtdt=-(rg/(paph(jl,jk+1)-paph(jl,jk)))/rcpd*           &
                     (pmfus(jl,jk)+pmfds(jl,jk)+rlmlt*pdpmel(jl,jk)-zalv* &
                     (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
                ptent(jl,jk)=ptent(jl,jk)+zdtdt
                !     
                !     increment of dry static energy
                !     
                penth(jl,jk)= zdtdt*rcpd
                !     
                pcond(jl,jk)=(rg/(paph(jl,jk+1)-paph(jl,jk)))/rcpd* &
                     (rlmlt*pdpmel(jl,jk)-zalv*                      &
                     (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))      
                !                                                                      
                zdqdt=-(rg/(paph(jl,jk+1)-paph(jl,jk)))*             &
                     (pmfuq(jl,jk)+pmfdq(jl,jk)+                     &
                     (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
                ptenq(jl,jk)=ptenq(jl,jk)+zdqdt
                !     
                !     cloud water tendency is managed by sundqvist
                !     --------------------------------------------
                !     here is the right place to substract 
                !     environmental compensating subsidence from
                !     pmful (not done in cuflx)
                !     
                !     zdldt=-(rg/(paph(jl,jk+1)-paph(jl,jk)))*
                !     1            ((pmful(jl,jk)-pmfu(jl,jk)*plenh(jl,jk))+
                !     1            pmfdl(jl,jk)+(pdmfup(jl,jk)+pdmfdp(jl,jk)))
                zdldt=0.0_realkind
                !     
                ptenl(jl,jk)=ptenl(jl,jk)+zdldt
             endif
140       enddo
       endif
150 enddo
    return
  end subroutine cudtdq



  subroutine cududv( kidia,    kfdia,    klon,     ktdia,    klev &
       , ktopm2,   ktype,    kcbot,    ldcum,    paph              &
       , puen,     pven,     pmfu,     pmfd                        &
       , puu,      pud,      pvu,      pvd                         &
       , ptenu,    ptenv  )
    !     
    !     
    !cududv - updates u and v tendencies,
    !     does global diagnostic of dissipation
    !     
    !     m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89
    !     

    !     ----------
    !     
    !     *cududv* is called from *cumastr*
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *ktype*        type of convection
    !     1 = penetrative convection
    !     2 = shallow convection
    !     3 = midlevel convection
    !     *kcbot*        cloud base level
    !     
    !     input parameters (logical):
    !     
    !     *ldcum*        flag: .true. for convective points
    !     
    !     input parameters (real):
    !     
    !     *paph*         provisional pressure on half levels            pa
    !     *puen*         provisional environment u-velocity (t+1)       m/s
    !     *pven*         provisional environment v-velocity (t+1)       m/s
    !     *pmfu*         massflux updrafts                          kg/(m2*s)
    !     *pmfd*         massflux downdrafts                        kg/(m2*s)
    !     *puu*          u-velocity in updrafts                         m/s
    !     *pud*          u-velocity in downdrafts                       m/s
    !     *pvu*          v-velocity in updrafts                         m/s
    !     *pvd*          v-velocity in downdrafts                       m/s
    !     
    !     updated parameters (real):
    !     
    !     *ptenu*        tendency of u-comp. of wind                    m/s2
    !     *ptenv*        tendency of v-comp. of wind                    m/s2
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of envt: 97-04-23 (j. a. garcia-moya, inm)


    implicit none
    integer:: jk,ktopm2,ik,jl,kidia,kfdia,ikb,klon,ktdia,klev
    real(kind=realkind):: zzp,zdudt,zdvdt



    !     
    real(kind=realkind)::     puen(klon,klev),        pven(klon,klev), &
         ptenv(klon,klev),       ptenu(klon,klev),    &
         paph(klon,klev+1)
    real(kind=realkind)::     puu(klon,klev),         pud(klon,klev), &
         pvu(klon,klev),         pvd(klon,klev),     &
         pmfu(klon,klev),        pmfd(klon,klev)
    integer::  ktype(klon),            kcbot(klon)
    logical::  ldcum(klon)
    !     
    real(kind=realkind)::     zmfuu(klon,klev),       zmfdu(klon,klev),&
         zmfuv(klon,klev),       zmfdv(klon,klev)
    !     
    !----------------------------------------------------------------------
    !     
    !     *    1.0          calculate fluxes and update u and v tendencies
    !     ----------------------------------------------
    !     

    do 120 jk=ktopm2,klev
       ik=jk-1
       do 110 jl=kidia,kfdia
          if(ldcum(jl)) then
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          endif
110    enddo
120 enddo
    do 140 jk=ktopm2,klev
       do 130 jl=kidia,kfdia
          if(ldcum(jl).and.jk>kcbot(jl)) then
             ikb=kcbot(jl)
             zzp=((paph(jl,klev+1)-paph(jl,jk))/ &
                           (paph(jl,klev+1)-paph(jl,ikb)))
             if(ktype(jl)==3) then
                zzp=zzp**2.0_realkind
             endif
             zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
             zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
             zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
             zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
          endif
130    enddo
140 enddo
    !     
    !     
    do 190 jk=ktopm2,klev
       !     
       if(jk<klev) then
          do 160 jl=kidia,kfdia
             if(ldcum(jl)) then
                zdudt=(rg/(paph(jl,jk+1)-paph(jl,jk)))*        &
                     (zmfuu(jl,jk+1)-zmfuu(jl,jk)+ &
                     zmfdu(jl,jk+1)-zmfdu(jl,jk))
                zdvdt=(rg/(paph(jl,jk+1)-paph(jl,jk)))* &
                     (zmfuv(jl,jk+1)-zmfuv(jl,jk)+ &
                     zmfdv(jl,jk+1)-zmfdv(jl,jk))
                ptenu(jl,jk)=ptenu(jl,jk)+zdudt
                ptenv(jl,jk)=ptenv(jl,jk)+zdvdt
             endif
160       enddo
          !     
       else
          do 170 jl=kidia,kfdia
             if(ldcum(jl)) then
                zdudt=-(rg/(paph(jl,jk+1)-paph(jl,jk)))* &
                     (zmfuu(jl,jk)+zmfdu(jl,jk))
                zdvdt=-(rg/(paph(jl,jk+1)-paph(jl,jk)))* &
                     (zmfuv(jl,jk)+zmfdv(jl,jk))
                ptenu(jl,jk)=ptenu(jl,jk)+zdudt
                ptenv(jl,jk)=ptenv(jl,jk)+zdvdt
             endif
170       enddo
       endif
       !     
190 enddo
    !     
    !     
    return
  end subroutine cududv


  subroutine cuflx(kidia,kfdia,klon,ktdia,klev             &
       ,  ptsphy                                            &
       ,  pten,     pqen,     plen                          &
       ,  ptenh,    pqenh,    plenh                         &
       ,  paph,     pgeoh,    pap,      ldland,   ldcum     &
       ,  kcbot,    kctop,    kdtop,    ktopm2              &
       ,  ktype,    lddraf                                  &
       ,  pmfu,     pmfd,     pmfus,    pmfds               &
       ,  pmfuq,    pmfdq,    pmful,    pmfdl,    plude     &
       ,  pleen,    pdmfup,   pdmfdp,   pqsen,    plglac    &
       ,  pdpmel,   pmflxr,   pmflxs,   prain )
    !     
    !     m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89
    !     
    !     purpose
    !     -------
    !     
    !     this routine does the final calculation of convective
    !     fluxes in the cloud layer and in the subcloud layer
    !     
    !     interface
    !     ---------
    !     this routine is called from *cumastr*.
    !     
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        end point
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *kcbot*        cloud base level
    !     *kctop*        cloud top level
    !     *kdtop*        top level of downdrafts
    !     
    !     input parameters (logical):
    !     
    !     *ldland*       land sea mask (.true. for land)
    !     *ldcum*        flag: .true. for convective points
    !     
    !     input parameters (real):
    !     
    !     *ptsphy*       time step for the physics                       s
    !     *pten*         provisional environment temperature (t+1)       k
    !     *pqen*         provisional environment spec. humidity (t+1)  kg/kg
    !     *plen*         provisional environment cloud wat. co. (t+1)  kg/kg
    !     *ptenh*        env. temperature (t+1) on half levels           k
    !     *pqenh*        env. spec. humidity (t+1) on half levels      kg/kg
    !     *plenh*        env. cloud water co. (t+1) on half levels     kg/kg
    !     *paph*         provisional pressure on half levels            pa
    !     *pgeoh*        geopotential on half levels                   m2/s2
    !     *pap*          provisional pressure full levels               pa
    !     
    !     updated parameters (integer):
    !     
    !     *ktype*        set to zero if ldcum=.false.
    !     
    !     updated parameters (logical):
    !     
    !     *lddraf*       set to .false. if ldcum=.false. or kdtop<kctop
    !     
    !     updated parameters (real):
    !     
    !     *pmfu*         massflux in updrafts                      kg/(m2*s)
    !     *pmfd*         massflux in downdrafts                    kg/(m2*s)
    !     *pmfus*        flux of dry static energy in updrafts      j/(m2*s)
    !     *pmfds*        flux of dry static energy in downdrafts    j/(m2*s)
    !     *pmfuq*        flux of spec. humidity in updrafts        kg/(m2*s)
    !     *pmfdq*        flux of spec. humidity in downdrafts      kg/(m2*s)
    !     *pmful*        flux of liquid water in updrafts          kg/(m2*s)
    !     *pmfdl*        flux of liquid water in downdrafts        kg/(m2*s)
    !     *plude*        detrained liquid water                    kg/(m3*s)
    !     *pleen*        entrained liquid water                    kg/(m3*s)
    !     *pdmfup*       flux difference of precip. in updrafts    kg/(m2*s)
    !     *pdmfdp*       flux difference of precip. in downdrafts  kg/(m2*s)
    !     *pqsen*        environment spec. saturation humidity (t+1)   kg/kg
    !     *plglac*       flux of frozen cloud water                kg/(m2*s)
    !     
    !     output parameters (real):
    !     
    !     *pdpmel*       change in precip.-fluxes due to melting   kg/(m2*s)
    !     *pmflxr*       convective rain flux                      kg/(m2*s)
    !     *pmflxs*       convective snow flux                      kg/(m2*s)
    !     *prain*        total precip. produced in conv. updrafts  kg/(m2*s)
    !     (no evaporation in downdrafts)
    !     
    !     externals
    !     ---------
    !     none
    !     
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)
    !     



    implicit none


    integer:: jk,ik,jl,kidia,kfdia,klon,ktdia,klev,itop,ktopm2,ikb
    real(kind=realkind):: ztmst,ptsphy,zcons1,zcons2,zcucov,ztmelp2,zdp,zzp,zfac
    real(kind=realkind):: zsnmlt,ztmsmlt,zalfaw,zrfl,zrnew,zrmin,zrfln,zdrfl,zdenom
    real(kind=realkind):: zpdr,zpds




    !     
    real(kind=realkind)::     pten(klon,klev),        ptenh(klon,klev), &
         pqen(klon,klev),        pqenh(klon,klev),      &
         plen(klon,klev),        plenh(klon,klev),      &
         paph(klon,klev+1),      pgeoh(klon,klev),      &
         pap(klon,klev)                                  
    !                                                         
    real(kind=realkind)::     pmfu(klon,klev),        pmfd(klon,klev),   &
         pmfus(klon,klev),       pmfds(klon,klev),      &
         pmfuq(klon,klev),       pmfdq(klon,klev),      &
         pdmfup(klon,klev),      pdmfdp(klon,klev),     &
         pqsen(klon,klev),       plglac(klon,klev),     &
         pdpmel(klon,klev),      prain(klon),           &
         pmful(klon,klev),       pmfdl(klon,klev),      &
         plude(klon,klev),       pleen(klon,klev),      &
         pmflxr(klon,klev+1),    pmflxs(klon,klev+1)     
    integer::  kcbot(klon),            kctop(klon),       &
         kdtop(klon),            ktype(klon)             
    logical::  lddraf(klon),           ldland(klon),      &
         ldcum(klon)



    !     *             specify constants

    ztmst=ptsphy
    zcons1=rcpd/(rlmlt*rg*ztmst)
    zcons2=1.0_realkind/(rg*ztmst)
    zcucov=0.05_realkind
    ztmelp2=rtt  +2.0_realkind
    !     
    !     
    !     
    !     *    1.0          determine final convective fluxes
    !     ---------------------------------
    !     

    itop=klev
    do 110 jl=kidia,kfdia
       prain(jl)=0.0_realkind
       itop=min(itop,kctop(jl))
       if(.not.ldcum(jl).or.kdtop(jl)<kctop(jl)) lddraf(jl)=.false.
       if(.not.ldcum(jl)) ktype(jl)=0
110 enddo
    ktopm2=itop-2
    do 120 jk=ktopm2,klev
       do 115 jl=kidia,kfdia
          if(ldcum(jl).and.jk>=kctop(jl)-1) then
             pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)* &
                  (rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
             pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
             !     
             !     pmful is also condensate, so environmental 
             !     compensating subsidence can not be added here 
             !     it will be in cudtdq 
             !     -------------------- 
             !     
             !     pmful(jl,jk)=pmful(jl,jk)-pmfu(jl,jk)*plenh(jl,jk)
             !     
             plglac(jl,jk)=pmfu(jl,jk)*plglac(jl,jk)
             if(ldland(jl)) then
                zdp=3.0e4_realkind
             else
                zdp=1.5e4_realkind
             endif
             if(paph(jl,kcbot(jl))-paph(jl,kctop(jl))>=zdp.and. &
                  pqen(jl,jk-1)>0.8_realkind*pqsen(jl,jk-1)) &
                                !     
                                !     entraiment from environment is not allowed
                                !     ------------------------------------------
                                !     
                                !     &      pdmfup(jl,jk-1)=pdmfup(jl,jk-1)+plude(jl,jk-1)-
                                !     &      pleen(jl,jk-1)
                  pdmfup(jl,jk-1)=pdmfup(jl,jk-1)+plude(jl,jk-1)
             !     
             if(lddraf(jl).and.jk>=kdtop(jl)) then
                pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)* &
                     (rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
                pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
                pmfdl(jl,jk)=pmfdl(jl,jk)-pmfd(jl,jk)*plenh(jl,jk)
             else
                pmfd(jl,jk)=0.0_realkind
                pmfds(jl,jk)=0.0_realkind
                pmfdq(jl,jk)=0.0_realkind
                pmfdl(jl,jk)=0.0_realkind
                pdmfdp(jl,jk-1)=0.0_realkind
             end if
          else
             pmfu(jl,jk)=0.0_realkind
             pmfd(jl,jk)=0.0_realkind
             pmfus(jl,jk)=0.0_realkind
             pmfds(jl,jk)=0.0_realkind
             pmfuq(jl,jk)=0.0_realkind
             pmfdq(jl,jk)=0.0_realkind
             pmful(jl,jk)=0.0_realkind
             pmfdl(jl,jk)=0.0_realkind
             plglac(jl,jk)=0.0_realkind
             pdmfup(jl,jk-1)=0.0_realkind
             pdmfdp(jl,jk-1)=0.0_realkind
             plude(jl,jk-1)=0.0_realkind
             pleen(jl,jk-1)=0.0_realkind
          end if
115    enddo
120 enddo
    do 130 jk=ktopm2,klev
       do 125 jl=kidia,kfdia
          if(ldcum(jl).and.jk>kcbot(jl)) then
             ikb=kcbot(jl)
             zzp=((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
             if(ktype(jl)==3) then
                zzp=zzp**2.0_realkind
             endif
             pmfu(jl,jk)=pmfu(jl,ikb)*zzp
             pmfus(jl,jk)=pmfus(jl,ikb)*zzp
             pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
             pmful(jl,jk)=pmful(jl,ikb)*zzp
          end if
125    enddo
130 enddo
    !     
    !     *    2.            calculate rain/snow fall rates
    !     *                  calculate melting of snow
    !     *                  calculate evaporation of precip
    !     -------------------------------
    !     

    !     
    !     initialazing raind and snow fluxes
    !     
    do 210 jk=1,klev+1
       do 205 jl=kidia,kfdia
          pmflxr(jl,jk)=0.0_realkind
          pmflxs(jl,jk)=0.0_realkind
205    enddo
210 enddo
    !     
    do 220 jk=ktopm2,klev
       do 215 jl=kidia,kfdia
          if(ldcum(jl)) then
             prain(jl)=prain(jl)+pdmfup(jl,jk)
             if(pmflxs(jl,jk)>0.0_realkind.and.pten(jl,jk)>ztmelp2) then
                zfac=zcons1*(paph(jl,jk+1)-paph(jl,jk))
                zsnmlt=min(pmflxs(jl,jk),zfac*(pten(jl,jk)-ztmelp2))
                pdpmel(jl,jk)=zsnmlt
                ztmsmlt=pten(jl,jk)-zsnmlt/zfac
                pqsen(jl,jk)=foeewm(ztmsmlt)/pap(jl,jk)
             end if
             zalfaw=foealfa(pten(jl,jk))
             pmflxr(jl,jk+1)=pmflxr(jl,jk)+zalfaw* &
                  (pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)
             pmflxs(jl,jk+1)=pmflxs(jl,jk)+(1.0_realkind-zalfaw)* &
                  (pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk)
             !     if(pten(jl,jk)>rtt  ) then
             !     pmflxr(jl,jk+1)=pmflxr(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)
             !     &       +pdpmel(jl,jk)
             !     pmflxs(jl,jk+1)=pmflxs(jl,jk)-pdpmel(jl,jk)
             !     else
             !     pmflxs(jl,jk+1)=pmflxs(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)
             !     pmflxr(jl,jk+1)=pmflxr(jl,jk)
             !     end if
             if(pmflxr(jl,jk+1)+pmflxs(jl,jk+1)<0.0_realkind) then
                pdmfdp(jl,jk)=-(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
                pmflxr(jl,jk+1)=0.0_realkind
                pmflxs(jl,jk+1)=0.0_realkind
                pdpmel(jl,jk)  =0.0_realkind
             endif
          end if
215    enddo
220 enddo
    do 240 jk=ktopm2,klev
       do 235 jl=kidia,kfdia
          if(ldcum(jl).and.jk>=kcbot(jl)) then
             zrfl=pmflxr(jl,jk)+pmflxs(jl,jk)
             if(zrfl>1.e-20_realkind) then
                zrnew=(max(0.0_realkind,sqrt(zrfl/zcucov)-                     &
                     cevapcu(jk)*(paph(jl,jk+1)-paph(jl,jk))*         &
                     max(0.0_realkind,pqsen(jl,jk)-pqen(jl,jk))))**2.0_realkind*zcucov      
                zrmin=zrfl-zcucov*max(0.0_realkind,0.8_realkind*pqsen(jl,jk)-pqen(jl,jk))&
                     *zcons2*(paph(jl,jk+1)-paph(jl,jk))
                zrnew=max(zrnew,zrmin)
                zrfln=max(zrnew,0.0_realkind)
                zdrfl=min(0.0_realkind,zrfln-zrfl)
                zdenom=1.0_realkind/max(1.e-20_realkind,pmflxr(jl,jk)+pmflxs(jl,jk))
                if(pten(jl,jk)>rtt  ) then
                   zpdr=pdmfdp(jl,jk)
                   zpds=0.0_realkind
                else
                   zpds=pdmfdp(jl,jk)
                   zpdr=0.0_realkind
                endif
                pmflxr(jl,jk+1)=pmflxr(jl,jk)+zpdr &
                     +pdpmel(jl,jk)+zdrfl*pmflxr(jl,jk)*zdenom
                pmflxs(jl,jk+1)=pmflxs(jl,jk)+zpds &
                     -pdpmel(jl,jk)+zdrfl*pmflxs(jl,jk)*zdenom
                pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
             else
                pmflxr(jl,jk+1)=0.0_realkind
                pmflxs(jl,jk+1)=0.0_realkind
                pdmfdp(jl,jk)=0.0_realkind
                pdpmel(jl,jk)=0.0_realkind
             endif
          end if
235    enddo
240 enddo
    !     
    return
  end subroutine cuflx


  subroutine cuini( kidia,    kfdia,    klon,     ktdia,    klev &
       , pten,     pqen,     pqsen,    puen,     pven             &
       , plen,     pvervel,  pgeo,     paph                       &
       , klwmin,   klab                                           &
       , ptenh,    pqenh,    pqsenh,   pgeoh,    plenh            &
       , ptu,      pqu,      ptd,      pqd                        &
       , puu,      pvu,      pud,      pvd                        &
       , plu,      pld    )
    !     
    !     m.tiedtke         e.c.m.w.f.     12/89
    !     
    !     purpose
    !     -------
    !     
    !     this routine interpolates large-scale fields of t,q etc.
    !     to half levels (i.e. grid for massflux scheme),
    !     determines level of maximum vertical velocity
    !     and initializes values for updrafts and downdrafts
    !     
    !     interface
    !     ---------
    !     this routine is called from *cumastr*.
    !     
    !     method.
    !     --------
    !     for extrapolation to half levels see tiedtke(1989)
    !     
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        end point
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     
    !     input parameters (real):
    !     
    !     *pten*         provisional environment temperature (t+1)       k
    !     *pqen*         provisional environment spec. humidity (t+1)  kg/kg
    !     *pqsen*        environment spec. saturation humidity (t+1)   kg/kg
    !     *puen*         provisional environment u-velocity (t+1)       m/s
    !     *pven*         provisional environment v-velocity (t+1)       m/s
    !     *plen*         provisional envir. cloud water cont. (t+1)    kg/kg
    !     *pvervel*      vertical velocity                             pa/s
    !     *pgeo*         geopotential                                  m2/s2
    !     *paph*         provisional pressure on half levels             pa
    !     
    !     output parameters (integer):
    !     
    !     *klwmin*       level of maximum vertical velocity
    !     *klab*         flag klab=1 for subcloud levels
    !     klab=2 for condensation level
    !     
    !     output parameters (real):
    !     
    !     *ptenh*        env. temperature (t+1) on half levels         k
    !     *pqenh*        env. spec. humidity (t+1) on half levels    kg/kg
    !     *plenh*        env. cloud water cont. (t+1) on half levels kg/kg
    !     *pqsenh*       env. spec. saturation humidity (t+1)
    !     on half levels                              kg/kg
    !     *pgeoh*        geopotential on half levels                 m2/s2
    !     *ptu*          temperature in updrafts                       k
    !     *pqu*          spec. humidity in updrafts                  kg/kg
    !     *ptd*          temperature in downdrafts                     k
    !     *pqu*          spec. humidity in downdrafts                kg/kg
    !     *puu*          u-velocity in updrafts                       m/s
    !     *pvu*          v-velocity in updrafts                       m/s
    !     *pud*          u-velocity in downdrafts                     m/s
    !     *pvd*          v-velocity in downdrafts                     m/s
    !     *plu*          liquid water content in updrafts            kg/kg
    !     *pld*          liquid water content in downdrafts          kg/kg
    !     
    !     
    !     externals
    !     ---------
    !     *cuadjtq* to specify qs at half levels
    !     
    !     modifications
    !     -------------
    !     92-09-21 : update to cy44      j.-j. morcrette
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)
    !     

    implicit none
    real(kind=realkind):: zdp,zzs
    integer:: jk,jl,kidia,kfdia,ik,icall,ktdia,klon,klev



    !     
    real(kind=realkind)::     pten(klon,klev),        pqen(klon,klev), &
         puen(klon,klev),        pven(klon,klev),      &
         plen(klon,klev),                              &
         pqsen(klon,klev),       pvervel(klon,klev),   &
         pgeo(klon,klev),        pgeoh(klon,klev),     &
         paph(klon,klev+1),      ptenh(klon,klev),     &
         pqenh(klon,klev),       pqsenh(klon,klev),    &
         plenh(klon,klev)                               
    !                                                        
    real(kind=realkind)::     ptu(klon,klev),         pqu(klon,klev),   &
         ptd(klon,klev),         pqd(klon,klev),       &
         puu(klon,klev),         pvu(klon,klev),       &
         pud(klon,klev),         pvd(klon,klev),       &
         plu(klon,klev),         pld(klon,klev)
    integer::  klab(klon,klev),        klwmin(klon)
    !     
    real(kind=realkind)::     zwmax(klon)
    real(kind=realkind)::     zph(klon)
    logical::  llflag(klon)
    !     
    !----------------------------------------------------------------------
    !     
    !     *    1.           specify large scale parameters at half levels
    !     *                 adjust temperature fields if staticly unstable
    !     *                 find level of maximum vertical velocity
    !     ----------------------------------------------
    !     

    zdp=0.5_realkind
    do 130 jk=2,klev
       do 110 jl=kidia,kfdia
          pgeoh(jl,jk)=pgeo(jl,jk)+(pgeo(jl,jk-1)-pgeo(jl,jk))*zdp
          ptenh(jl,jk)=(max(rcpd*pten(jl,jk-1)+pgeo(jl,jk-1), &
               rcpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))/rcpd
          pqsenh(jl,jk)=pqsen(jl,jk-1)
          zph(jl)=paph(jl,jk)
          llflag(jl)=.true.
110    enddo
       !     
       ik=jk
       icall=0
       call cuadjtq( kidia,    kfdia,    klon,     ktdia,    klev &
            , ik, zph,      ptenh,    pqsenh,   llflag,   icall)
       !     
       do 120 jl=kidia,kfdia
          pqenh(jl,jk)=min(pqen(jl,jk-1),pqsen(jl,jk-1)) &
               +(pqsenh(jl,jk)-pqsen(jl,jk-1))
          pqenh(jl,jk)=max(pqenh(jl,jk),0.0_realkind)
          plenh(jl,jk)=plen(jl,jk-1)+(pqsen(jl,jk-1)-pqsenh(jl,jk))
          plenh(jl,jk)=max(plenh(jl,jk),0.0_realkind)
120    enddo
130 enddo
    !     
    do 140 jl=kidia,kfdia
       ptenh(jl,klev)=(rcpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))/rcpd
       pqenh(jl,klev)=pqen(jl,klev)
       plenh(jl,klev)=plen(jl,klev)
       ptenh(jl,1)=pten(jl,1)
       pqenh(jl,1)=pqen(jl,1)
       plenh(jl,1)=plen(jl,1)
       pgeoh(jl,1)=pgeo(jl,1)
       klwmin(jl)=klev
       zwmax(jl)=0.0_realkind
140 enddo
    !     
    do 160 jk=klev-1,2,-1
       do 150 jl=kidia,kfdia
          zzs=max(rcpd*ptenh(jl,jk)+pgeoh(jl,jk),&
               rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
          ptenh(jl,jk)=(zzs-pgeoh(jl,jk))/rcpd
150    enddo
160 enddo
    !     
    do 190 jk=klev,3,-1
       do 180 jl=kidia,kfdia
          if(pvervel(jl,jk)<zwmax(jl)) then
             zwmax(jl)=pvervel(jl,jk)
             klwmin(jl)=jk
          end if
180    enddo
190 enddo
    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     *    2.0          initialize values for updrafts and downdrafts
    !     *                 ---------------------------------------------
    !     

    do 230 jk=1,klev
       ik=jk-1
       if(jk==1) ik=1
       do 220 jl=kidia,kfdia
          ptu(jl,jk)=ptenh(jl,jk)
          ptd(jl,jk)=ptenh(jl,jk)
          pqu(jl,jk)=pqenh(jl,jk)
          pqd(jl,jk)=pqenh(jl,jk)
          !     
          !     interaction between old a new clouds is not allowed
          !     ---------------------------------------------------
          !     
          !     plu(jl,jk)=plenh(jl,jk)
          !     pld(jl,jk)=plenh(jl,jk)
          plu(jl,jk)=0.0_realkind
          pld(jl,jk)=0.0_realkind
          !     
          puu(jl,jk)=puen(jl,ik)
          pud(jl,jk)=puen(jl,ik)
          pvu(jl,jk)=pven(jl,ik)
          pvd(jl,jk)=pven(jl,ik)
          klab(jl,jk)=0
220    enddo
230 enddo
    !     
    return
  end subroutine cuini



  subroutine cubasmc(kidia,    kfdia,    klon,     ktdia,    klev,  &
       kk,                                                     &
       pten,     pqen,     pqsen,    plen,     puen,     pven, &
       pvervel,  pgeo,     pgeoh,    ldcum,    ktype,    klab, &
       kcbot,    pmfu,     pmfub,    pentr,                    &
       ptu,      pqu,      plu,      puu,      pvu,            &
       pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,    pmfuv)
    !     
    !     m.tiedtke         e.c.m.w.f.     12/89
    !     
    !     purpose.
    !     --------
    !     this routine calculates cloud base values
    !     for midlevel convection
    !     
    !     interface
    !     ---------
    !     
    !     this routine is called from *cuasc*.
    !     input are environmental values t,q etc
    !     it returns cloudbase values for midlevel convection
    !     
    !     method.
    !     --------
    !     s. tiedtke (1989)
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *kk*           actual level
    !     
    !     input parameters (real):
    !     
    !     *pten*         provisional environment temperature (t+1)       k
    !     *pqen*         provisional environment spec. humidity (t+1)  kg/kg
    !     *pqsen*        environment spec. saturation humidity (t+1)   kg/kg
    !     *plen*         provisional environment cloud wat. co. (t+1)  kg/kg
    !     *puen*         provisional environment u-velocity (t+1)       m/s
    !     *pven*         provisional environment v-velocity (t+1)       m/s
    !     *pvervel*      vertical velocity                             pa/s
    !     *pgeo*         geopotential                                  m2/s2
    !     *pgeoh*        geopotential on half levels                   m2/s2
    !     
    !     input parameters (logical):
    !     
    !     *ldcum*        flag: .true. for convective points
    !     
    !     updated parameters (integer):
    !     
    !     *ktype*        type of convection
    !     1 = penetrative convection
    !     2 = shallow convection
    !     3 = midlevel convection
    !     *klab*         flag klab=1 for subcloud levels
    !     klab=2 for cloud levels
    !     *kcbot*        cloud base level
    !     
    !     output parameters (real):
    !     
    !     *pmfu*         massflux in updrafts                     kg/(m2*s)
    !     *pmfub*        massflux in updrafts at cloud base       kg/(m2*s)
    !     *pentr*        fractional mass entrainment rate          1/m
    !     *ptu*          temperature in updrafts                    k
    !     *pqu*          spec. humidity in updrafts               kg/kg
    !     *plu*          liquid water content in updrafts         kg/kg
    !     *puu*          u-velocity in updrafts                    m/s
    !     *pvu*          v-velocity in updrafts                    m/s
    !     *pmfus*        flux of dry static energy in updrafts     j/(m2*s)
    !     *pmfuq*        flux of spec. humidity in updrafts       kg/(m2*s)
    !     *pmful*        flux of liquid water in updrafts         kg/(m2*s)
    !     *pdmfup*       flux difference of precip. in updrafts   kg/(m2*s)
    !     *pmfuu*        flux of u-momentum in updrafts     (m/s)*kg/(m2*s)
    !     *pmfuv*        flux of v-momentum in updrafts     (m/s)*kg/(m2*s)
    !     
    !     externals
    !     ---------
    !     none
    !     
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of environment: 97-04-23 (j. a. garcia-moya, inm)



    implicit none
    integer:: kk,jl,kidia,kfdia,klon,klev,ktdia
    real(kind=realkind):: zzzmb




    real(kind=realkind)::     pten(klon,klev),        pqen(klon,klev), &
         plen(klon,klev),                              &
         puen(klon,klev),        pven(klon,klev),      &
         pqsen(klon,klev),       pvervel(klon,klev),   &
         pgeo(klon,klev),        pgeoh(klon,klev)       
    !                                                        
    real(kind=realkind)::     ptu(klon,klev),         pqu(klon,klev),   &
         puu(klon,klev),         pvu(klon,klev),       &
         plu(klon,klev),         pmfu(klon,klev),      &
         pmfub(klon),            pentr(klon),          &
         pmfus(klon,klev),       pmfuq(klon,klev),     &
         pmful(klon,klev),       pdmfup(klon,klev),    &
         pmfuu(klon),            pmfuv(klon)            
    integer::  ktype(klon),            kcbot(klon),      &
         klab(klon,klev)
    logical::  ldcum(klon)
    !     
    logical::  llo2
    !     
    !     
    !----------------------------------------------------------------------
    !     
    !     *    1.           calculate entrainment and detrainment rates
    !     -------------------------------------------
    !     
    if(lmfmid.and.kk<klev-1.and.kk>klev/2) then
       do 150 jl=kidia,kfdia
          if(.not.ldcum(jl).and.klab(jl,kk+1)==0 &
               .and.pqen(jl,kk)>0.90_realkind*pqsen(jl,kk)) then
             ptu(jl,kk+1)=(rcpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))/rcpd
             pqu(jl,kk+1)=pqen(jl,kk)
             !     
             !     interaction between old clouds 
             !     and new ones is not allowed 
             !     ---------------------------
             !     
             !     plu(jl,kk+1)=plen(jl,kk)
             plu(jl,kk+1)=0.0_realkind
             !     
             zzzmb=max(cmfcmin,-pvervel(jl,kk)/rg)
             zzzmb=min(zzzmb,cmfcmax)
             pmfub(jl)=zzzmb
             pmfu(jl,kk+1)=pmfub(jl)
             pmfus(jl,kk+1)=pmfub(jl)*(rcpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
             pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
             pmful(jl,kk+1)=pmfub(jl)*plu(jl,kk+1)
             pdmfup(jl,kk+1)=0.0_realkind
             kcbot(jl)=kk
             klab(jl,kk+1)=1
             ktype(jl)=3
             pentr(jl)=entrmid
             if(lmfdudv) then
                puu(jl,kk+1)=puen(jl,kk)
                pvu(jl,kk+1)=pven(jl,kk)
                pmfuu(jl)=pmfub(jl)*puu(jl,kk+1)
                pmfuv(jl)=pmfub(jl)*pvu(jl,kk+1)
             endif
          endif
150    enddo
    endif

    return
  end subroutine cubasmc


  subroutine cuentr( kidia,    kfdia,    klon,     ktdia,    klev    &
       , kk,       klwmin,   ktype,    kcbot,    kctop0         &
       , ldcum,    ldwork                                       &
       , ptenh,    pqenh,    plenh,    ptenq,    paph,     pap  &
       , pmfu,     pentr                                        &
       , pcbase,   pdmfen,   pdmfde )
    !     
    !     m.tiedtke         e.c.m.w.f.     12/89
    !     
    !     purpose.
    !     --------
    !     this routine calculates entrainment/detrainment rates
    !     for updrafts in cumulus parameterization
    !     
    !     interface
    !     ---------
    !     
    !     this routine is called from *cuasc*.
    !     input are environmental values t,q etc
    !     and updraft values t,q etc
    !     it returns entrainment/detrainment rates
    !     
    !     method.
    !     --------
    !     s. tiedtke (1989)
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer::):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *kk*           current level
    !     *klwmin*       level of maximum vertical velocity
    !     *ktype*        type of convection
    !     1 = penetrative convection
    !     2 = shallow convection
    !     3 = midlevel convection
    !     *kcbot*        cloud base level
    !     *kctop0*       first guess of cloud top level
    !     
    !     input parameters (logical):
    !     
    !     *ldcum*        flag: .true. for convective points
    !     
    !     input parameters (real):
    !     
    !     *ptenh*        env. temperature (t+1) on half levels         k
    !     *pqenh*        env. spec. humidity (t+1) on half levels    kg/kg
    !     *plenh*        env. cloud water co. (t+1) on half levels   kg/kg
    !     *ptenq*        moisture tendency                          kg/(kg s)
    !     *paph*         provisional pressure on half levels          pa
    !     *pap*          provisional pressure on full levels          pa
    !     *pmfu*         massflux in updrafts                       kg/(m2*s)
    !     *pentr*        fractional mass entrainment rate             1/m
    !     
    !     output parameters (real):
    !     
    !     *pcbase*       pressure at cloud base                       pa
    !     *pdmfen*       entrainment rate                           kg/(m2*s)
    !     *pdmfde*       detrainment rate                           kg/(m2*s)
    !     
    !     externals
    !     ---------
    !     none
    !     
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of envt: 97-04-23 (j. a. garcia-moya, inm)
    !     
    !----------------------------------------------------------------------


    implicit none
    integer:: jk,jl,kidia,kfdia,kk,iklwmin,klon,klev,ktdia
    real(kind=realkind):: zdprho,zentr,zpmid




    real(kind=realkind)::     ptenh(klon,klev),       pqenh(klon,klev), &
         plenh(klon,klev),                              &
         pap(klon,klev),         paph(klon,klev+1),     &
         ptenq(klon,klev),       pmfu(klon,klev),       &
         pentr(klon),            pcbase(klon)            
    integer::  klwmin(klon),           ktype(klon),       &
         kcbot(klon),            kctop0(klon)
    logical::  ldcum(klon)

    real(kind=realkind)::     pdmfen(klon),           pdmfde(klon)

    real(kind=realkind)::     zrho(klon),             zptop(klon)

    logical::  llo1,llo2,ldwork

    !     *    1.           calculate entrainment and detrainment rates

    if(ldwork) then


       do 105 jl=kidia,kfdia
          pdmfen(jl)=0.0_realkind
          pdmfde(jl)=0.0_realkind
          zrho(jl)=paph(jl,kk+1)/(rd*ptenh(jl,kk+1))
          pcbase(jl)=paph(jl,kcbot(jl))
          zptop(jl)=paph(jl,kctop0(jl))
105    enddo



       !     *    1.2          specify entrainment rates for deep clouds

       do 125 jl=kidia,kfdia
          if(ldcum(jl)) then
             zdprho=(paph(jl,kk+1)-paph(jl,kk))/(rg*zrho(jl))
             zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho
             llo1=kk<kcbot(jl)
             if(llo1) pdmfde(jl)=zentr
             zpmid=0.5_realkind*(pcbase(jl)+zptop(jl))
             llo2=llo1.and.ktype(jl)==2.and. &
                  (pcbase(jl)-paph(jl,kk)<0.2e5_realkind.or. &
                  paph(jl,kk)>zpmid)
             if(llo2) pdmfen(jl)=zentr
             iklwmin=max(klwmin(jl),kctop0(jl)+2)
             llo2=llo1.and.(ktype(jl)==1.or.ktype(jl)==3).and. &
                  (kk>=iklwmin.or.pap(jl,kk)>zpmid)
             if(llo2) pdmfen(jl)=zentr
             llo1=pdmfen(jl)>0.0_realkind.and.(ktype(jl)==1.or.ktype(jl)==2)
             if(llo1) then
                zentr=zentr* &
                     (1.0_realkind+3.0_realkind*(1.0_realkind-min(1.0_realkind,(pcbase(jl)-pap(jl,kk))/1.5e4_realkind)))
                pdmfen(jl)=pdmfen(jl)* &
                     (1.0_realkind+3.0_realkind*(1.0_realkind-min(1.0_realkind,(pcbase(jl)-pap(jl,kk))/1.5e4_realkind)))
                pdmfde(jl)=pdmfde(jl)* &
                     (1.0_realkind+3.0_realkind*(1.0_realkind-min(1.0_realkind,(pcbase(jl)-pap(jl,kk))/1.5e4_realkind)))
             endif
             if(llo2.and.pqenh(jl,kk+1)>1.e-5_realkind) &
                  pdmfen(jl)=zentr+max(ptenq(jl,kk),0.0_realkind)/pqenh(jl,kk+1)* &
                  zrho(jl)*zdprho
          endif
125    enddo

    endif

    return
  end subroutine cuentr


  subroutine cuadjtq(kidia,    kfdia,    klon,     ktdia,    klev, &
       kk, psp,      pt,       pq,       ldflag,   kcall)
    !     
    !     
    !     m.tiedtke         e.c.m.w.f.     12/89
    !     
    !     modifications
    !     -------------
    !     d.salmond         cray(uk))      12/8/91
    !     j.j. morcrette    ecmwf          92-09-18   update to cy44
    !     adapted to hirlam: 96-01-23 (j. a. garcia-moya, inm)
    !     mod. for cloud water of envt: 97-04-23 (j. a. garcia-moya, inm)
    !     
    !     purpose.
    !     --------
    !     to produce t,q and l values for cloud ascent
    !     
    !     interface
    !     ---------
    !     this routine is called from subroutines:
    !     *cond*     (t and q at condensation level)
    !     *cubase*   (t and q at condensation level)
    !     *cuasc*    (t and q at cloud levels)
    !     *cuini*    (environmental t and qs values at half levels)
    !     *custrat*  (t and q at condensation level)
    !     input are unadjusted t and q values,
    !     it returns adjusted values of t and q
    !     
    !     parameter     description                                   units
    !     ---------     -----------                                   -----
    !     input parameters (integer):
    !     
    !     *kidia*        start point
    !     *kfdia*        endpoint
    !     *klon*         number of grid points per packet
    !     *ktdia*        start of the vertical loop
    !     *klev*         number of levels
    !     *kk*           level
    !     *kcall*        defines calculation as
    !     kcall=0  env. t and qs in*cuini*
    !     kcall=1  condensation in updrafts  (e.g. cubase, cuasc)
    !     kcall=2  evaporation in downdrafts (e.g. cudlfs,cuddraf)
    !     
    !     input parameters (logical):
    !     
    !     *ldland*       land-sea mask (.true. for land points)
    !     
    !     input parameters (real):
    !     
    !     *psp*          pressure                                        pa
    !     
    !     updated parameters (real):
    !     
    !     *pt*           temperature                                     k
    !     *pq*           specific humidity                             kg/kg
    !     
    !     
    !     externals
    !     ---------
    !     3 lookup tables ( tlucua, tlucub, tlucuc )
    !     for condensation calculations.
    !     the tables are initialised in *suphec*.
    !     
    !----------------------------------------------------------------------
    !     

    implicit none



    integer:: kcall,isum,jl,kidia,kfdia,kk
    integer:: klon,klat,klev,ktdia
    real(kind=realkind):: zqsat,zcor,zcond1
    real(kind=realkind)::  pt(klon,klev), pq(klon,klev), psp(klon)
    logical:: ldflag(klon)
    real(kind=realkind):: zcond(klon), zqp(klon)


    !     2.           calculate condensation and adjust t and q accordingly


    if (kcall==1 ) then

       isum=0
       do 210 jl=kidia,kfdia
          if(ldflag(jl)) then
             zqp(jl)=1.0_realkind/psp(jl)
             !     lucup      it=pt(jl,kk)*1000.
             !     lucup      zqsat=tlucua(it)*zqp(jl)
             !     
             zqsat=foeewm(pt(jl,kk))*zqp(jl)
             !     
             zqsat=min(0.5_realkind,zqsat)
             zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
             zqsat=zqsat*zcor
             !     lucup      zcond(jl)=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
             zcond(jl)=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
             zcond(jl)=max(zcond(jl),0.0_realkind)
             !     lucup      pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
             pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond(jl)
             pq(jl,kk)=pq(jl,kk)-zcond(jl)
             if(abs(zcond(jl))>1.e-14_realkind) isum=isum+1
          else
             zcond(jl)=0.0_realkind
          endif
210    enddo
       !     

       if(isum/=0)then
          !     
          do 220 jl=kidia,kfdia
             if(ldflag(jl).and.abs(zcond(jl))>1.e-14_realkind) then
                !     lucup      it=pt(jl,kk)*1000.
                !     lucup      zqsat=tlucua(it)*zqp(jl)
                !     
                zqsat=foeewm(pt(jl,kk))*zqp(jl)
                !     
                zqsat=min(0.5_realkind,zqsat)
                zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
                zqsat=zqsat*zcor
                !     lucup      zcond1=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
                zcond1=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
                !     lucup      pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
                pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
                pq(jl,kk)=pq(jl,kk)-zcond1
             endif
220       enddo
          !     
       endif

       !     
    endif

    if(kcall==2) then
       !     
       isum=0
       do 310 jl=kidia,kfdia
          if(ldflag(jl)) then
             zqp(jl)=1.0_realkind/psp(jl)
             !     lucup      it=pt(jl,kk)*1000.
             !     lucup      zqsat=tlucua(it)*zqp(jl)
             !     
             zqsat=foeewm(pt(jl,kk))*zqp(jl)
             !     
             zqsat=min(0.5_realkind,zqsat)
             zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
             zqsat=zqsat*zcor
             !     lucup      zcond(jl)=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
             zcond(jl)=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
             zcond(jl)=min(zcond(jl),0.0_realkind)
             !     lucup      pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
             pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond(jl)
             pq(jl,kk)=pq(jl,kk)-zcond(jl)
             if(abs(zcond(jl))>1.e-14_realkind) isum=isum+1
          else
             zcond(jl)=0.0_realkind
          endif
310    enddo
       !     

       if(isum/=0)then
          !     
          do 320 jl=kidia,kfdia
             if(ldflag(jl).and.abs(zcond(jl))>1.e-14_realkind) then
                !     lucup      it=pt(jl,kk)*1000.
                !     lucup      zqsat=tlucua(it)*zqp(jl)
                !     
                zqsat=foeewm(pt(jl,kk))*zqp(jl)
                !     
                zqsat=min(0.5_realkind,zqsat)
                zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
                zqsat=zqsat*zcor
                !     lucup      zcond1=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
                zcond1=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
                !     lucup      pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
                pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
                pq(jl,kk)=pq(jl,kk)-zcond1
             endif
320       enddo
          !     
       endif

       !     
    endif
    !     
    if(kcall==0) then
       !     
       do 410 jl=kidia,kfdia
          zqp(jl)=1.0_realkind/psp(jl)
          !     lucup   it=pt(jl,kk)*1000.
          !     lucup   zqsat=tlucua(it)*zqp(jl)
          !     
          zqsat=foeewm(pt(jl,kk))*zqp(jl)
          !     
          zqsat=min(0.5_realkind,zqsat)
          zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
          zqsat=zqsat*zcor
          !     lucup   zcond1=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          zcond1=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
          !     lucup   pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
          pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk)=pq(jl,kk)-zcond1
410    enddo
       !     
       do 420 jl=kidia,kfdia
          !     lucup   it=pt(jl,kk)*1000.
          !     lucup   zqsat=tlucua(it)*zqp(jl)
          !     
          zqsat=foeewm(pt(jl,kk))*zqp(jl)
          !     
          zqsat=min(0.5_realkind,zqsat)
          zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
          zqsat=zqsat*zcor
          !     lucup   zcond1=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          zcond1=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
          !     lucup   pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
          pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk)=pq(jl,kk)-zcond1
420    enddo



    endif
    !     
    if(kcall==4) then
       !     
       do 510 jl=kidia,kfdia
          zqp(jl)=1.0_realkind/psp(jl)
          !     lucup   it=pt(jl,kk)*1000.
          !     lucup   zqsat=tlucua(it)*zqp(jl)
          !     
          zqsat=foeewm(pt(jl,kk))*zqp(jl)
          !     
          zqsat=min(0.5_realkind,zqsat)
          zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
          zqsat=zqsat*zcor
          !     lucup   zcond(jl)=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          zcond(jl)=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
          !     lucup   pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
          pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond(jl)
          pq(jl,kk)=pq(jl,kk)-zcond(jl)
510    enddo
       !     
       do 520 jl=kidia,kfdia
          !     lucup   it=pt(jl,kk)*1000.
          !     lucup   zqsat=tlucua(it)*zqp(jl)
          !     
          zqsat=foeewm(pt(jl,kk))*zqp(jl)
          !     
          zqsat=min(0.5_realkind,zqsat)
          zcor=1.0_realkind/(1.0_realkind-retv  *zqsat)
          zqsat=zqsat*zcor
          !     lucup   zcond1=(pq(jl,kk)-zqsat)/(1.+zqsat*zcor*tlucub(it))
          zcond1=(pq(jl,kk)-zqsat)/(1.0_realkind+zqsat*zcor*foedem(pt(jl,kk)))
          !     lucup   pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
          pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk)=pq(jl,kk)-zcond1
520    enddo
       !     
    endif
    !     
    !     
    return
  end subroutine cuadjtq


  subroutine prcond (nhor,nlev,kstart,kstop,t,hldcp2)

    !     prcond:  prepare several parameters for the sundqvist scheme

    use confys
    use escom
    use ccons
    implicit none

    integer:: nhor,nlev,kstart,kstop
    real(kind=realkind):: t(nhor,nlev),q(nhor,nlev),pf(nhor,nlev)
    real(kind=realkind):: hldcp2(nhor,nlev)

    integer:: jk,jl,jintt

    logical:: lqsice
    real(kind=realkind):: cfreez,aiw,x,y,xt,riw(nhor),hedldr

    !     lqsice      if = true, then qs is weighted with
    !     ice probability for t < 273 k
    cfreez=0.12_realkind
    lqsice=.true.
    hedldr = epsilo*latice/rair
    aiw=exp(hedldr*(1.0_realkind/tmelt-1.0_realkind/tcir))

    do jk=1,nlev
       do jl=kstart,kstop
          if(lqsice) then 

             !     ice crystal probability
             jintt = max(1, 1 + nint((tmelt - t(jl,jk))/dttabl))
             jintt = min(jintt, 1750)

             !     calculate ice - water saturation ratio  riw
             xt = min(tmelt, t(jl,jk))  
             xt = max(xt, 173.0_realkind)
             x = xt  - tcir
             y = elotci/xt
             !     riw is expanded around tcir (but instead of using 1, 1.0015607 is
             !     put in, to compensate that y above uses xt and not a mean of
             !     t and tcir...?
             riw(jl) = aiw*(1.0015607_realkind+x*y*(1.0_realkind+0.5_realkind*x*y*(1.0_realkind+x*y/3.0_realkind)))
          endif
          !     
          !     as before:
          hldcp2(jl,jk) =(latvap + latice* prbice(jintt))/cpair
       enddo
    enddo
    !     
    return
  end subroutine prcond



  subroutine astamic(nhor,nlev,kstart,kstop,                     &
       dtime,                                       &
       jwanv,kht,khb,                               &
       jwmask,                                      &
       frland,ts,ps,                                &
       bhyb,                                        &
       t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,          &
       hcucov,hdcwcu,dlnpdt,hldcp2,hu,hqsat,hsq,    &
       dtdt,dqdt,dcwdt,                             &
       sevapr,sevapc,                               &
       totcov,prcpcu,prcpst,cusnow,stsnow,cov2d,    &
       cwpath,                                      &
       draindt,dsnowdt,                             &
       hdcwst,                                      &
       prcl,prci,prsl,prsi)
    !
    ! stamic  - stratiform condensation and microphysic parameterization of
    !           the scheme (including cloud water content forecast)

    implicit none
    !
    !
    !

    !
    !  nhor   - number of horizontal gridpoints (nlon*nlat)
    !  nlev   - number of vertical levels
    !  kstart - starting gridpoint of the call
    !  kstop  - final gridpoint of the call
    !
    integer:: nhor,nlev,kstart,kstop
    !
    !  dtime  - 2*dt
    !
    real(kind=realkind):: dtime,conacc
    !

    !
    !  jwanv  - gridpoints where anvil is allowed to becomes stratiform
    !  kht    - top of the convective clouds
    !  khb    - base of the convective clouds
    !
    integer:: jwanv(nhor),kht(nhor),khb(nhor)
    !
    !  frland - fraction of land in the square grid
    !  ts     - surface temperature
    !  ps     - surface pressure
    !
    real(kind=realkind):: frland(nhor),ts(nhor),ps(nhor)
    !

    !
    !  bhyb   - parameter b of the hybrid coordinate of the half model level
    !
    real(kind=realkind):: bhyb(nlev+1)
    !

    !
    !  jwmask - set = 0 if convection,   otherwise = 1
    !
    integer:: jwmask(nhor,nlev)
    !
    !  t      - temperature
    !  q      - specific humidity
    !  cw     - cloud water content
    !  dcwdin - cloud water tendendy due to other processes than convection
    !           (dynamics, radiation, turbulence fluxes, ..)
    !  dtdtin - temperature tendendy equal than dcwdin
    !  dqdtin - q tendendy equal than dcwdin
    !  pf     - pressure at the full model levels
    !  dpf    - delta p of the model layers
    !  dlnpdt - = 1/p*dp/dt at level
    !  hldcp2 - 1/cpair * latent heat as a function of
    !           temperature by weighting in latent heat of
    !           freezing times probability of ice crystals
    !  hu     - relative humidity
    !  hqsat  - saturation mixing ratio with respect to
    !                    gridpoint temperature
    !  hsq    - eps*l**2*qsat/(r*cp*t**2)
    !  hcucov - convective cloud cover
    !  hdcwcu - tendency of cloud water due to convection
    !
    real(kind=realkind):: t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),dcwdin(nhor,nlev),  &
         dtdtin(nhor,nlev),dqdtin(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev),&
         dlnpdt(nhor,nlev),hldcp2(nhor,nlev),hu(nhor,nlev),hsq(nhor,nlev),&
         hqsat(nhor,nlev),hcucov(nhor,nlev),hdcwcu(nhor,nlev)
    !     

    !
    !  dtdt   - temperature tendency
    !  dqdt   - q tendency
    !  dcwdt  - cloud water tendency
    !
    real(kind=realkind):: dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev)
    !

    !
    !  sevapr - accumulated  evaporation from precipitation
    !  sevapc - accumulated evaporation of cloud water
    !
    real(kind=realkind):: sevapr,sevapc
    !

    !
    !  cov2d  - 2d cloud cover from maximum/random overlapping
    !  cwpath - cloud water vertically integrated
    !  draindt- tendency of liquid precipitation
    !  dsnowdt- tendency of snow
    !  prcpcu - convective precipitation (rain and snow)
    !  prcpst - stratiform precipitation
    !  stsnow - stratiform snow
    !  cusnow - convective snow
    !
    real(kind=realkind):: cov2d(nhor),cwpath(nhor),draindt(nhor),dsnowdt(nhor),prcpcu(nhor), &
         prcpst(nhor),cusnow(nhor),stsnow(nhor)
    !

    !
    !  hdcwst - tendency of cloud water due to stratiform
    !  totcov - total cloud cover
    !
    real(kind=realkind):: hdcwst(nhor,nlev),totcov(nhor,nlev)
    !
    !  prsl   - flux of liquid precipitation due to stratiform part
    !  prsi   - flux of solid precipitation due to stratiform part
    !  prcl   - flux of liquid precipitation due to convection
    !  prci   - flux of solid precipitation due to convection
    !
    real(kind=realkind):: prcl(nhor,0:nlev),prci(nhor,0:nlev),prsl(nhor,0:nlev),prsi(nhor,0:nlev)
    !

    !
    !                       work space:
    !

    !
    !

    !
    !  jwqcon - indicator = 1 if stratiform condensation is allowed
    !                    to take place,  else jwqcon = 0
    !  jwssat - index for stratiform condensation
    !
    !

    !
    !  hfoorl - function giving the vertical influence of ocean -
    !                    - land contrast for hu00
    !

    !
    !  huz00  - u00 for every grid point in every vertical level
    !  wrk1   - working array number 1
    !  wrk2   - working array number 2
    !  wrk3   - working array number 3
    !  wrk4   - working array number 4
    !  tprel  - preliminary value of t(t+dt)
    !  qsprel - preliminary value of qs(t+dt)
    !  xfprim - m(tau-1)/b.mr, exponent of the precipitation with coalescenc
    !           effect equation
    !  xfrcoa - fco as in osz94
    !  hfmrx  - actual mr in the microphysics scale
    !  xevap  - square root second term in eq. of precipitation evaporation
    !  evapri - precipitation evaporation out of the cloud
    !  duoorl - u00 difference between land and sea
    !  dustab - increment of u00 due to stability (postulated equation)
    !  vterm  - terminal velocity of precipitation, weighted
    !           from hvsnow and hvterm with resp to ice probab
    !

    !
    !  hdqad  - tendency of vapour due to effects other than condensation
    !  hdqst  - tendency of vapour due to stratiform condensation
    !  hdtst  - tendency of temperature due to stratiform condensation
    !  hdcwad - tendency of cloud water due to other effects
    !           but condensation after evaporation of cloud water
    !  hstcov - contains stratiform cloud cover
    !  wqp    - gridpoint q at (t+dt)
    !  wro    - density of air
    !  hevac  - cloud water evaporation up to the saturation of air
    !  hevapr - rate of precipitation evaporation in the layer
    !  cndsat - condensate

    integer:: jwqcon(nhor),jwssat(nhor)
    real(kind=realkind):: hfoorl(nlev)
    real(kind=realkind):: huz00(nhor),wrk1(nhor),wrk2(nhor),wrk3(nhor),         &
         wrk4(nhor),tprel(nhor),qsprel(nhor),xfprim(nhor),      &
         xfrcoa(nhor),hfmrx(nhor),xevap(nhor),evapri(nhor),     &
         duoorl(nhor),dustab(nhor),vterm(nhor),storcu(nhor)      
    real(kind=realkind):: hdqad(nhor,nlev),hdqst(nhor,nlev),                     &
         hdtst(nhor,nlev),hdcwad(nhor,nlev),hstcov(nhor,nlev),  &
         wqp(nhor,nlev),wro(nhor,nlev),                         &
         hevac(nhor,nlev),hevapr(nhor,nlev),cndsat(nhor,nlev)         

    call stamic (nhor,nlev,kstart,kstop,                        &
         dtime,                                         &
         jwanv,kht,khb,                                 &
         jwmask,                                        &
         frland,ts,ps,bhyb,                             &
         t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,            &
         dlnpdt,hldcp2,hu,hqsat,hsq,hcucov,hdcwcu,      &
         dtdt,dqdt,dcwdt,                               &
         sevapr,sevapc,                                 &
         jwqcon,jwssat,                                 &
         prcpcu,prcpst,cusnow,stsnow,cov2d,cwpath,      &
         draindt,dsnowdt,                               &
         huz00,wrk1,wrk2,wrk3,wrk4,tprel,qsprel,        &
         xfprim,xfrcoa,hfmrx,xevap,evapri,              &
         duoorl,dustab,vterm,hfoorl,storcu,             &
         totcov,hdcwst,hdqad,hdqst,hdtst,               &
         hdcwad,hstcov,wqp,                             &
         wro,hevac,hevapr,cndsat,                       &
         prcl,prci,prsl,prsi)
    return
  end subroutine astamic

  subroutine stamic (nhor,nlev,kstart,kstop,   &
       dtime,                                   &
       jwanv,kht,khb,                           &
       jwmask,                                  &
       frland,ts,ps,bhyb,                       &
       t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,      &
       dlnpdt,hldcp2,hu,hqsat,hsq,hcucov,hdcwcu,&
       dtdt,dqdt,dcwdt,                         &
       sevapr,sevapc,                           &
       jwqcon,jwssat,                           &
       prcpcu,prcpst,cusnow,stsnow,cov2d,cwpath,&
       draindt,dsnowdt,                         &
       huz00,wrk1,wrk2,wrk3,wrk4,tprel,qsprel,  &
       xfprim,xfrcoa,hfmrx,xevap,evapri,        &
       duoorl,dustab,vterm,hfoorl,storcu,       &
       totcov,hdcwst,hdqad,hdqst,hdtst,         &
       hdcwad,hstcov,wqp,                       &
       wro,hevac,hevapr,cndsat,                 &
       prcl,prci,prsl,prsi)

    !     stamic  - stratiform condensation and microphysic parameterization of
    !     the scheme (including cloud water content forecast)


    use confys
    use condsmod
    implicit none
    !     
    !     

    !     
    !     nhor   - number of horizontal gridpoints (nlon*nlat)
    !     nlev   - number of vertical levels
    !     kstart - starting gridpoint of the call
    !     kstop  - final gridpoint of the call
    !     
    integer:: nhor,nlev,kstart,kstop
    !     
    !     dtime  - 2*dt
    !     
    real(kind=realkind):: dtime,conacc
    !     

    !     
    !     jwanv  - gridpoints where anvil is allowed to becomes stratiform
    !     kht    - top of the convective clouds
    !     khb    - base of the convective clouds
    !     
    integer:: jwanv(nhor),kht(nhor),khb(nhor)
    !     
    !     frland - fraction of land in the square grid
    !     ts     - surface temperature
    !     ps     - surface pressure
    !     
    real(kind=realkind):: frland(nhor),ts(nhor),ps(nhor)
    !     

    !     
    !     bhyb   - parameter b of the hybrid coordinate of the half model level
    !     
    real(kind=realkind):: bhyb(nlev+1)


    !     
    !     jwmask - set = 0 if convection,   otherwise = 1
    !     
    integer:: jwmask(nhor,nlev)
    !     
    !     t      - temperature
    !     q      - specific humidity
    !     cw     - cloud water content
    !     dcwdin - cloud water tendendy due to other processes than convection
    !     (dynamics, radiation, turbulence fluxes, ..)
    !     dtdtin - temperature tendendy equal than dcwdin
    !     dqdtin - q tendendy equal than dcwdin
    !     pf     - pressure at the full model levels
    !     dpf    - delta p of the model layers
    !     dlnpdt - = 1/p*dp/dt at level
    !     hldcp2 - 1/cpair * latent heat as a function of
    !     temperature by weighting in latent heat of
    !     freezing times probability of ice crystals
    !     hu     - relative humidity
    !     hqsat  - saturation mixing ratio with respect to
    !     gridpoint temperature
    !     hsq    - eps*l**2*qsat/(r*cp*t**2)
    !     hcucov - convective cloud cover
    !     hdcwcu - tendency of cloud water due to convection
    !     
    real(kind=realkind):: t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),dcwdin(nhor,nlev),      &
         dtdtin(nhor,nlev),dqdtin(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev),&
         dlnpdt(nhor,nlev),hldcp2(nhor,nlev),hu(nhor,nlev),hsq(nhor,nlev),&
         hqsat(nhor,nlev),hcucov(nhor,nlev),hdcwcu(nhor,nlev)
    !     

    !     
    !     dtdt   - temperature tendency
    !     dqdt   - q tendency
    !     dcwdt  - cloud water tendency
    !     
    real(kind=realkind):: dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev)


    !     
    !     sevapr - accumulated  evaporation from precipitation
    !     sevapc - accumulated evaporation of cloud water
    !     
    real(kind=realkind):: sevapr,sevapc
    !     

    !     
    !     jwqcon - indicator = 1 if stratiform condensation is allowed
    !     to take place,  else jwqcon = 0
    !     jwssat - index for stratiform condensation
    !     
    integer:: jwqcon(nhor),jwssat(nhor)
    !     
    !     cov2d  - 2d cloud cover from maximum/random overlapping
    !     cwpath - cloud water vertically integrated
    !     draindt- tendency of liquid precipitation
    !     dsnowdt- tendency of snow
    !     prcpcu - convective precipitation (rain and snow)
    !     prcpst - stratiform precipitation
    !     stsnow - stratiform snow
    !     cusnow - convective snow
    !     
    real(kind=realkind):: cov2d(nhor),cwpath(nhor), draindt(nhor),dsnowdt(nhor), &
         prcpcu(nhor),prcpst(nhor),stsnow(nhor),cusnow(nhor)
    !     
    !     huz00  - u00 for every grid point in every vertical level
    !     wrk1   - working array number 1
    !     wrk2   - working array number 2
    !     wrk3   - working array number 3
    !     wrk4   - working array number 4
    !     tprel  - preliminary value of t(t+dt)
    !     qsprel - preliminary value of qs(t+dt)
    !     xfprim - m(tau-1)/b.mr, exponent of the precipitation with coalescenc
    !     effect equation
    !     xfrcoa - fco as in osz94
    !     hfmrx  - actual mr in the microphysics scale
    !     xevap  - square root second term in eq. of precipitation evaporation
    !     evapri - precipitation evaporation out of the cloud
    !     duoorl - u00 difference between land and sea
    !     dustab - increment of u00 due to stability (postulated equation)
    !     vterm  - terminal velocity of precipitation, weighted
    !     from hvsnow and hvterm with resp to ice probab
    !     storcu - used to distinguish between stratiform, storcu = 0,
    !     and convective, storcu = 1, precipitation
    !     
    real(kind=realkind):: huz00(nhor),wrk1(nhor),wrk2(nhor),wrk3(nhor),tprel(nhor), &
         qsprel(nhor),xfprim(nhor),xfrcoa(nhor),hfmrx(nhor),        &
         wrk4(nhor),vterm(nhor),xevap(nhor),evapri(nhor),           &
         duoorl(nhor),dustab(nhor),storcu(nhor)

    !     hfoorl - function giving the vertical influence of ocean -
    !     - land contrast for hu00
    real(kind=realkind):: hfoorl(nlev)
    !     

    !     
    !     hdcwst - tendency of cloud water due to stratiform
    !     totcov - total cloud cover
    !     
    real(kind=realkind):: hdcwst(nhor,nlev),totcov(nhor,nlev)
    !     
    !     prsl   - flux of liquid precipitation due to stratiform part
    !     prsi   - flux of solid precipitation due to stratiform part
    !     prcl   - flux of liquid precipitation due to convection
    !     prci   - flux of solid precipitation due to convection
    !     
    real(kind=realkind):: prcl(nhor,0:nlev),prci(nhor,0:nlev),prsl(nhor,0:nlev),prsi(nhor,0:nlev)
    !     
    !     hdqad  - tendency of vapour due to effects other than condensation
    !     hdqst  - tendency of vapour due to stratiform condensation
    !     hdtst  - tendency of temperature due to stratiform condensation
    !     hdcwad - tendency of cloud water due to other effects
    !     but condensation after evaporation of cloud water
    !     hstcov - contains stratiform cloud cover
    !     wqp    - gridpoint q at (t+dt)
    !     wro    - density of air
    !     hevac  - cloud water evaporation up to the saturation of air
    !     hevapr - rate of precipitation evaporation in the layer
    !     cndsat - condensate
    !     
    real(kind=realkind):: hdqad(nhor,nlev),hdqst(nhor,nlev),hdtst(nhor,nlev), &
         hdcwad(nhor,nlev),hstcov(nhor,nlev),wqp(nhor,nlev),  &
         wro(nhor,nlev),                                      &
         hevac(nhor,nlev),hevapr(nhor,nlev),cndsat(nhor,nlev)
    !     


    !     work space:
    !     

    !     
    integer:: jl,jk,jintt,nlevm1,nlevp1,nlevp2
    real(kind=realkind):: cstocu(nhor)
    real(kind=realkind):: cbfsno
    real(kind=realkind):: zxpsq,zxp,zxf,zxhj,zxn,zxnpm,dt,dt2,zwrk1,zwrk3,zwrk11, &
         zwrk25,zwrk51,zwrk26,zwrk17,zwrk35,zwrk52,zwrk57,zwrk2,  &
         zwrk6,xsp
    !     references:
    !     
    !     sundqvist, h.; berge, e.; kristjansson, j.e. (1989).-"condensation an
    !     cloud parameterization studies with a mesoscale numerical weather
    !     prediction model". mwr, 117, 1641-1657. later will be sbk89.
    !     
    !     sundqvist, h. (1991).-personal notes. later will be s91.
    !     
    !     sundqvist, h. (1993).-"inclusion of ice phase of hydrometeors in clou
    !     parameterization for mesoscale and largescale models". beit. phys.
    !     atmosph., 66,137-147. later will be s93.
    !     
    !     olofsson, p-o.; sundqvist, h.; zurovac-jevtic, d. (1994).-"documentat
    !     of the routine condens". misu. later will be osz94.
    !     
    !     
    !_______________________________________________________________________
    !     
    !     the tendencies of specific humidity have to be made 2-dimensional
    !     in routine dyn.
    !_______________________________________________________________________
    !     
    !     this routine deals with parameterization of stratiform condensati
    !     and microphysics.
    !     
    !     points 11) - 15) handle parameterization of stratiform condensa-
    !     tion including prediction of cloud water content and
    !     
    !     point 16) handles parameterization of the microphysical processes
    !     for both convective and stratiform cases, comprising
    !     prediction of cloud water content and
    !     evaporation of precipitation
    !     
    !     points 17) - 18) arrange various calculated fields in the routine
    !     
    !     ii)  declarations
    !     iv)  values of constants including derived ones
    !     v)   preparations
    !     
    !-----------------------------------------------------------------------
    !     
    !     11)   modification of basic threshold relative humidity, hu00
    !     12)   general requirement for condensation, jwqcon = 1 or 0
    !     13)   setting of specific hu00 and jwqcon
    !     14)   cloud cover and evaporation of cloud water and precipitation
    !     and a special treatment of cirrus cloud
    !     15)   accession of vapour, and final setting of heating/cooling
    !     
    !     16)   prediction of cloud water content, cw, for both convective and
    !     stratiform cases, employing newton-raphson iterative solution
    !     
    !     17)   accumulate precipitation,final setting of tendencies
    !     and calculate cov2d for output
    !     18)   quantities for ground parameterization in the hirlam-model
    !     
    !-----------------------------------------------------------------------
    !     ii)        declarations
    !--------------------------------------------------------------------

    !     
    !     iv)          values of constants including derived ones
    !     -------------------------------------------------------
    !     
    nlevm1 = nlev-1
    nlevp1 = nlev+1
    nlevp2 = nlev+2
    !     
    !     v)            preparations
    !     ---------------------------
    sevapr = 0.0_realkind
    sevapc = 0.0_realkind
    !     
    do 6 jl = kstart,kstop
       !     
       huz00(jl) = hu00
       prcpcu(jl) = 0.0_realkind
       prcpst(jl) = 0.0_realkind

       !     arpege stuff within subroutine condens
       !     
       prcl(jl,0)=0.0_realkind
       prci(jl,0)=0.0_realkind
       prsl(jl,0)=0.0_realkind
       prsi(jl,0)=0.0_realkind
       !     
       cusnow(jl) = 0.0_realkind
       stsnow(jl) = 0.0_realkind

       !     cw integrated and 2-d cover, max-random overlap
       !     
       cwpath(jl) = 0.0_realkind
       cov2d(jl)  = 1.0_realkind
       !     
       !     terminal velocity, cloudiness storage and u00 difference
       !     
       duoorl(jl) = 0.0_realkind
       vterm(jl)=hvsnow
       cstocu(jl)=hcst
       storcu(jl)=0.0_realkind
       !     
       evapri(jl)=0.0_realkind
       tprel(jl)=0.0_realkind
       qsprel(jl)=0.0_realkind
       wrk1(jl)=0.0_realkind
       wrk2(jl)=0.0_realkind
       wrk3(jl)=0.0_realkind
       wrk4(jl)=0.0_realkind
       xfprim(jl)=0.0_realkind
       xfrcoa(jl)=0.0_realkind
       hfmrx(jl)=0.0_realkind
       xevap(jl)=0.0_realkind
       !     
6   enddo

    !     update the vertical arrays
    !     
    do 20 jk = 1,nlev
       !     
       hfoorl(jk) = 1._realkind
       do 10 jl = kstart,kstop
          !     
          wqp(jl,jk)=q(jl,jk)+dtime*dqdtin(jl,jk)
          !     
          !     arpege stuff within subroutine condens
          !     
          prcl(jl,jk)=0.0_realkind
          prci(jl,jk)=0.0_realkind
          prsl(jl,jk)=0.0_realkind
          prsi(jl,jk)=0.0_realkind
          !     
          hdqad(jl,jk) = dqdtin(jl,jk)
          hdcwad(jl,jk) = dcwdin(jl,jk)
          !     
          hdcwst(jl,jk)= 0.0_realkind
          hevac(jl,jk) = 0.0_realkind
          hevapr(jl,jk)= 0.0_realkind
          hdtst(jl,jk) = 0.0_realkind
          hdqst(jl,jk) = 0.0_realkind
          !     
          hstcov(jl,jk) = 0.0_realkind
          !     
          totcov(jl,jk) = hcucov(jl,jk)
          !     
          if(jk>=kht(jl) .and. jk<=khb(jl) .and.jwmask(jl,jk)==0) then
             cndsat(jl,jk) = dtdt(jl,jk) / hldcp2(jl,jk)
             cndsat(jl,jk) = max(cndsat(jl,jk),0.0_realkind)
          else
             cndsat(jl,jk) = 0.0_realkind
          endif

          !     air density stored in wro
          !     
          wro(jl,jk) = pf(jl,jk)/(rair *t(jl,jk))
          !     
10     enddo
20  enddo
    !     
    !-----------------------------------------------------------------------

    !     11)      modification of basic threshold relative humidity, hu00
    !     
    !-----------------------------------------------------------------------
    !     
    !     ocean or land
    !     
    !     duoorl - u00 difference between land and sea
    !     frland - fraction of land in the square grid
    !     
    do 1110 jl = kstart,kstop
       if (frland(jl) > 0.5_realkind) duoorl(jl) = 0.1_realkind
1110 enddo
    !     

    !     
    !     beginning main vertical loop
    !     
    do 1670 jk = 1, nlev
       !     

       !     
       !-----------------------------------------------------------------------

       !     12)      general requirement for condensation, jwqcon = 1 or 0
       !     
       !-----------------------------------------------------------------------
       !     
       !     xevap  - square root second term in eq. 29 in osz94
       !     jwqcon - indicator = 1 if stratiform condensation is allowed
       !     to take place,  else jwqcon = 0
       !     jwssat - index for stratiform condensation
       !     jwmask - set = 0 if convection,   otherwise = 1
       !     hu00   - u00 default value over sea
       !     hfoorl - function giving the vertical influence of ocean -
       !     - land contrast for hu00
       !     huz00  - u00 for every grid point in every vertical level
       !     
       do 1210 jl = kstart,kstop
          !     
          xevap(jl) = 0.0_realkind
          jwqcon(jl)=jwmask(jl,jk)
          dustab(jl) = 0.0_realkind
          jwssat(jl) = 0
          huz00(jl) = hu00 - duoorl(jl) * hfoorl(jk)
          !     
1210   enddo

       !     wqp    - preliminary value of q(t+dt)
       !     hqsat  - qsat(t)
       !     hdqad  - advection of q
       !     
       do 1215 jl = kstart,kstop
          !     
          if (wqp(jl,jk)*real(jwmask(jl,jk),realkind)/ hqsat(jl,jk)> 1.02_realkind ) then
             jwqcon(jl) = 1
             jwssat(jl) = 1
             hdqad(jl,jk) = (wqp(jl,jk)-min(q(jl,jk),hqsat(jl,jk)))/ dtime
          endif
          !     
1215   enddo
       !     
       !-----------------------------------------------------------------------

       !     13)      setting of specific hu00
       !     
       !-----------------------------------------------------------------------
       !     
       !     
       !     xsp    - bhyb at full levels
       !     
       xsp=1._realkind
       if (jk/=nlev) xsp = 0.5_realkind*(bhyb(jk+1)+bhyb(jk+2))
       !     
       !     let hu00 increase linearly with p for p/ps > 0.9
       !     
       !     u00max - maximum value of u00 allowed
       !     
       do 1304 jl = kstart,kstop
          !     
          if(jwqcon(jl) == 1.and. pf(jl,jk)/ps(jl) > 0.9_realkind)  then
             huz00(jl) = huz00(jl) + (u00max-huz00(jl))*(1.0_realkind-10.0_realkind*(1.0_realkind-pf(jl,jk)/ps(jl)))
          endif
          !     
1304   enddo

       !     let hu00 increase with stability
       !     
       !     wrk1   - adequate temperature
       !     dustab - increment of u00 due to stability (postulated equation)
       !     wrk2   - maximum and actual u00 difference
       !     
       if (jk>=nlevm1) then
          !     
          do 1305 jl = kstart,kstop
             !     
             if(jk==nlevm1) then
                wrk1(jl) = t(jl,jk + 1)
             else
                wrk1(jl) = ts(jl)
             endif
             !     
             dustab(jl) = 0.0_realkind
             !     
1305      enddo
          !     
          if( abs(xsp-0.5_realkind*(bhyb(jk-1)+bhyb(jk))) > 1.e-14_realkind)  then
             !     
             do 1306 jl = kstart,kstop
                !     
                wrk2(jl) = u00max - huz00(jl)
                dustab(jl)=max(wrk2(jl) * (1.85_realkind-0.019_realkind*(t(jl,jk)-277.0_realkind)    &
                     -1.67e-2_realkind*(wrk1(jl)-t(jl,jk-1)) /(xsp-0.5_realkind*(bhyb(jk-1) &
                     +bhyb(jk)))) , 0.0_realkind)* real(jwqcon(jl),realkind)
                !     
1306         enddo
             !     
          endif
          !     
       endif
       !     
       do 1313 jl = kstart,kstop
          huz00(jl) = min((huz00(jl) + dustab(jl)),u00max)
1313   enddo

       !     let u00 increase for t<238k
       !     
       !     wrk1   - postulated equation for the variation of u00 with temperatur
       !     below 238k
       !     wrk2   - increment up to u00max
       !     
       do 1320 jl = kstart,kstop
          !     
          if (t(jl,jk) < 238.0_realkind ) then
             wrk1(jl) = 1.0_realkind+ 0.15_realkind * (238.0_realkind - t(jl,jk))
             wrk2(jl) = (u00max - huz00(jl)) * (1.0_realkind - 1.0_realkind / wrk1(jl))
             huz00(jl) = min((huz00(jl) + wrk2(jl)) , u00max)
          endif
          !     
1320   enddo

       !     t > tcir: cloud cover; but if u<u0, set jwqcon = 0;
       !     
       do 1330 jl = kstart,kstop
          !     
          if (hu(jl,jk)<=huz00(jl).and.jwssat(jl)==0 ) jwqcon(jl) = 0
          !     
1330   enddo
       !     
       !     totcov - total cloud cover
       !     hstcov - contains stratiform cloud cover
       !     
       do 1334 jl = kstart,kstop
          !     
          hstcov(jl,jk) = (1.0_realkind - sqrt((1.0_realkind - hu(jl,jk)) * 0.98_realkind/ &
               (1.0_realkind - huz00(jl)))) * real(jwqcon(jl),realkind)
          hstcov(jl,jk) = max(hstcov(jl,jk),0.0_realkind)
          totcov(jl,jk) = totcov(jl,jk) + hstcov(jl,jk)
          !     
1334   enddo
       !     

       !     condensation calculations
       !     
       !     
       !-----------------------------------------------------------------------
       !     
       !     14)      evaporation of cloud water and precipitation
       !     
       !-----------------------------------------------------------------------
       !     
       !     
       !     possible evaporation of cloud water, proportional to the cloud free a
       !     
       !     
       !     wrk1   -
       !     wrk2   -
       !     hevac  -
       !     hdcwad -
       !     
       do 1405 jl = kstart,kstop
          !     
          wrk1(jl)=max(hdcwad(jl,jk),0.0_realkind)*(1.0_realkind- hstcov(jl,jk)*real(jwmask(jl,jk),realkind) &
               - hcucov(jl,jk)*real(1-jwmask(jl,jk),realkind))                           
          wrk2(jl) = min (wrk1(jl),(((hqsat(jl,jk) - q(jl,jk))/dtime +     &
               hu(jl,jk)* hsq(jl,jk)*dtdtin(jl,jk)/hldcp2(jl,jk)           &
               -hdqad(jl,jk)) / (1.0_realkind + hsq(jl,jk))))
          hevac(jl,jk) = max(wrk2(jl),0.0_realkind)
          !     
1405   enddo

       !     evaporation of precipitation
       !     
       !     wrk2   - air density multiplied by precipitation terminal velocity
       !     wrk3   - minimum value in eq. 29 in osz94
       !     wrk4   - cloud cover square root
       !     xevap  - square root second term in eq. 29 in osz94
       !     cov2d  - cloud cover with maximum/random overlapping (the real cloud
       !     cover is 1-cov2d) (eq. 31 in osz94)
       !     
       if(jk<=2) then
          !     
          do 1414 jl = kstart,kstop
             !     
             wrk2(jl) = wro(jl,jk) * vterm(jl)
             wrk3(jl) = min((dpf(jl,jk)/(gravit*wrk2(jl))),dtime)
             cov2d(jl) = 1.0_realkind - totcov(jl,jk)
             wrk4(jl) = sqrt(1.0_realkind - cov2d(jl))
             !     
             xevap(jl) = 0.5_realkind*hkevap*(1.0_realkind - huz00(jl)*real(jwqcon(jl),realkind) - &
                  hu(jl,jk)* real(1  - jwqcon(jl),realkind)) * sqrt(vterm(jl)) &
                  * wrk3(jl) * wrk4(jl)
             !     
1414      enddo
          !     
       else
          !     
          do 1415 jl = kstart,kstop
             !     
             wrk2(jl) = wro(jl,jk) * vterm(jl)
             wrk3(jl) = min((dpf(jl,jk)/(gravit*wrk2(jl))),dtime)
             cov2d(jl) = cov2d(jl) * (1.0_realkind - max(totcov(jl,jk-1), &
                  totcov(jl,jk))) /                              &
                  (1.0_realkind - min(totcov(jl,jk-1),0.99_realkind))
             cov2d(jl)=min(cov2d(jl),1.0_realkind)
             wrk4(jl) = sqrt(1.0_realkind - cov2d(jl))
             !     
             xevap(jl) = 0.5_realkind*hkevap*(1.0_realkind - huz00(jl)*real(jwqcon(jl),realkind) - &
                  hu(jl,jk)* real(1  - jwqcon(jl),realkind)) * sqrt(vterm(jl)) &
                  * wrk3(jl) * wrk4(jl)
             !     
1415      enddo
          !     
       endif

       !     zxpsq  - square root of pout
       !     zxp    - precipitation out of the layer
       !     wrk1   - gues of the evaporation value
       !     evapri - precipitation evaporation out of the cloud
       !     hevapr - rate of precipitation evaporation in the layer
       !     stsnow - stratiform snow afther the evaporation
       !     prcpst - stratiform precipitation afther the evaporation
       !     
       do 1420 jl = kstart,kstop
          !     
          if(prcpst(jl) > 0.0_realkind .and. jwmask(jl,jk) == 1) then
             zxpsq = sqrt(prcpst(jl)) - xevap(jl)
             zxp = (max (zxpsq,0.0_realkind))**2.0_realkind
             wrk1(jl) = (prcpst(jl) -zxp) * real(1-jwssat(jl),realkind)
             evapri(jl) = (1.0_realkind- hstcov(jl,jk)) * (1.0_realkind- hstcov(jl,jk))*wrk1(jl)
             evapri(jl) = max(evapri(jl),0.0_realkind)
             hevapr(jl,jk) = evapri(jl) * gravit / dpf(jl,jk)
             stsnow(jl)=max(0.0_realkind,(stsnow(jl)-evapri(jl)*stsnow(jl)/prcpst(jl)))
             prcpst(jl) = prcpst(jl) - evapri(jl)
             stsnow(jl) = min(stsnow(jl), prcpst(jl))
          elseif(prcpcu(jl) > 0.0_realkind .and. jwmask(jl,jk) == 0) then
             zxpsq = sqrt(prcpcu(jl)) - xevap(jl)
             zxp = (max (zxpsq,0.0_realkind))**2.0_realkind
             wrk1(jl) = prcpcu(jl) -zxp
             evapri(jl) = (1.0_realkind- hcucov(jl,jk))*(1.0_realkind- hcucov(jl,jk))*wrk1(jl)
             evapri(jl) = max(evapri(jl),0.0_realkind)
             hevapr(jl,jk) = evapri(jl) * gravit / dpf(jl,jk)
             cusnow(jl)=max(0.0_realkind,(cusnow(jl)-evapri(jl)*cusnow(jl)/prcpcu(jl)))
             prcpcu(jl) = prcpcu(jl) - evapri(jl)
             cusnow(jl) = min(cusnow(jl), prcpcu(jl))
             prci(jl,jk)= cusnow(jl)
             prcl(jl,jk)=prcpcu(jl)-cusnow(jl)
             !     
             dtdt(jl,jk) = dtdt(jl,jk) - hevapr(jl,jk) * hldcp2(jl,jk)
             dqdt(jl,jk) = dqdt(jl,jk) + hevapr(jl,jk)
          endif
          !     
1420   enddo
       !     
       !-----------------------------------------------------------------------

       !     t > tcir
       !     
       !     15)         accession of vapour, and final setting of heating/cooling
       !     
       !-----------------------------------------------------------------------
       !     
       !     wrk1   - convergence of available latent heat (m, eq. 3.18 in sbk89)
       !     zxn    - 2.(1.b).(1.u00), working paramter
       !     zxnpm  - denominator of eq. 3.20 in sbk89
       !     wrk2   - working parameter
       !     cndsat - condensate
       !     
       do 1508 jl = kstart,kstop
          !     
          if(jwmask(jl,jk) == 1) then
             wrk1(jl) = (hdqad(jl,jk) - hu(jl,jk) * hsq(jl,jk) * &
                  dtdtin(jl,jk) /  hldcp2(jl,jk) + q(jl,jk) *     &
                  dlnpdt(jl,jk) ) * real(jwqcon(jl),realkind)
             !     
             !-----------------------------------------------------------------------
             !     
             !     this is the proper place for possible modification of
             !     acces due to entrainment in sc-layer
             !     
             !-----------------------------------------------------------------------
             !     
             !     calculation of heating rate composed of release of latent
             !     heat and cooling due to evaporation of precipitation
             !     
             zxn = 2.0_realkind * (1.0_realkind - huz00(jl))*(1.0_realkind-hstcov(jl,jk) + 1.e-3_realkind)
             zxnpm = zxn + cw(jl,jk)/(hstcov(jl,jk) * hqsat(jl,jk)+1.e-6_realkind)
             wrk2(jl) = zxnpm - zxn + hstcov(jl,jk) * zxn
             cndsat(jl,jk) = (wrk2(jl) * wrk1(jl) - zxn *(hevapr(jl,jk)+ &
                  hevac(jl,jk)))/((1.0_realkind+ hu(jl,jk)* hsq(jl,jk)) * zxnpm)    &
                  *real(jwqcon(jl),realkind)
             !     
          endif
          !     
1508   enddo
       !     
       do 1510 jl = kstart,kstop
          !     
          if(jwqcon(jl) == 1) then
             if ((cndsat(jl,jk) + cw(jl,jk)/dtime + &
                  hdcwad(jl,jk))<0.0_realkind.and.jwssat(jl)==0) jwqcon(jl) = 0
          endif
          !     
1510   enddo

       !     final setting of heating/cooling and cloud cover
       !     
       !     hstcov - stratiform cloud cover
       !     wrk1   - preliminary value of q(t+dt)
       !     tprel  - the same for t(t+dt)
       !     qsprel - the same for qs(t+dt)
       !     
       do 1520 jl = kstart,kstop
          if(jwmask(jl,jk) == 1.or.(jwmask(jl,jk)==0.and.&
               jk>khb(jl))) then
             cndsat(jl,jk) = cndsat(jl,jk) * real(jwqcon(jl),realkind) - &
                  real(1-jwqcon(jl),realkind) * (hevac(jl,jk) + hevapr(jl,jk))
             if( jwqcon(jl)==0.and.abs(cw(jl,jk))<1.e-14_realkind) then
                hdtst(jl,jk) = - hevac(jl,jk)*hldcp2(jl,jk)
                hdqst(jl,jk) = hevac(jl,jk)
             endif
             hstcov(jl,jk) = hstcov(jl,jk) * real(jwqcon(jl),realkind)
             totcov(jl,jk) = hstcov(jl,jk)

             !     check if supersaturation may result and if so, adjust
             !     
             if(jwmask(jl,jk)==1) then
                wrk1(jl) = q(jl,jk) + dtime*(hdqad(jl,jk)- cndsat(jl,jk))
                tprel(jl) = t(jl,jk) + dtime*(dtdtin(jl,jk)+ &
                     cndsat(jl,jk) * hldcp2(jl,jk))
                qsprel(jl) = hqsat(jl,jk)+hsq(jl,jk)*(tprel(jl) &
                     -t(jl,jk))/hldcp2(jl,jk) -hqsat(jl,jk)*dlnpdt(jl,jk)*dtime
             endif
          endif
          !     
1520   enddo
       !     
       do 1525 jl = kstart,kstop
          !     
          if(jwmask(jl,jk)==1) then
             if(wrk1(jl)>qsprel(jl))  cndsat(jl,jk) = cndsat(jl,jk)+ &
                  (wrk1(jl)-qsprel(jl))/(dtime*(1.0_realkind+hsq(jl,jk)))
          endif
          !     
1525   enddo

       !     let old cw evaporate if it exists without having condensation
       !     conditions fulfilled
       !     
       do 1530 jl = kstart,kstop
          if((jwqcon(jl) == 0 .and. cw(jl,jk)*real(jwmask(jl,jk),realkind) > 0.0_realkind) &
               .or. (cw(jl,jk)*(1._realkind-real(jwmask(jl,jk),realkind)) > 0.0_realkind &
               .and. jk > khb(jl))) then
             hevac(jl,jk) = max(0.0_realkind, cw(jl,jk) + dtime*hdcwad(jl,jk)) / dtime
             hdcwst(jl,jk) = - hevac(jl,jk)
             hdtst(jl,jk) = - hevac(jl,jk)*hldcp2(jl,jk)
             hdqst(jl,jk) = hevac(jl,jk)
          endif
1530   enddo

       !     let evaporatated cw influence the cloud water tendency in
       !     even in the case jwqcon(jl)=0 and cw(jl,jk)=0
       !     
       do 1540 jl = kstart,kstop
          if(jwqcon(jl) == 0 .and. abs(cw(jl,jk))< 1.e-14_realkind) then
             if(jwmask(jl,jk) == 1 .or. (jwmask(jl,jk) == 0 &
                  .and. jk > khb(jl))) then
                hevac(jl,jk) = max(0.0_realkind,hdcwad(jl,jk))
                hdcwst(jl,jk) = - hevac(jl,jk)
                cndsat(jl,jk) = - hevac(jl,jk)- hevapr(jl,jk)
             endif
          endif
1540   enddo
       !     
       !-----------------------------------------------------------------------

       !     t > tcir
       !     
       !     16)      prediction of cloud water content, cw
       !     for convective and stratiform clouds
       !     
       !-----------------------------------------------------------------------
       !     
       !     calculate factors for coalescence hfcox, freezing hfrezx,
       !     reduction of hmrcu or hmrst at low temps, hfmrx(jl)
       !     
       !     the various quantities in this loop are calculated with
       !     appropriate attention to whether the condensation is
       !     convective or stratiform
       !     
       !     modified probability of ice resulting from ice in precip
       !     coming from above
       !     
       !     jintt  - table index
       !     zwrk1  - ice crystal probability
       !     zwrk2  - vapour press diff water-ice
       !     stsnow - stratiform snow
       !     cusnow - convective snow
       !     prcpcu - convective percipitation
       !     storcu - used to distinguish between stratiform, storcu = 0,
       !     and convective, storcu = 1, precipitation
       !     zwrk3  - modified ice crystal probability (eq. 11 in s93)
       !     zwrk6  - bergeron-findeisen effect function (eq. 13 in s93)
       !     wrk1   - ice crystal probability
       !     wrk2   - vapour press diff water-ice
       !     wrk3   - modified ice crystal probability

       !     
       do 1605 jl=kstart,kstop
          !     
          jintt = max(1,1 + nint((tmelt - t(jl,jk))/dttabl))
          jintt = min(jintt, 1750)
          zwrk1  = prbice(jintt)
          zwrk2  = hdewi(jintt)
          !     
          if(kht(jl) < nlevp1 .and. jk == kht(jl) &
               .and. jwmask(jl,jk)==0) then
             cusnow(jl) = stsnow(jl)
             prcpcu(jl) = prcpst(jl)
             stsnow(jl) = 0.0_realkind
             prcpst(jl) = 0.0_realkind
             storcu(jl) = 1.0_realkind
          endif
          !     
          cusnow(jl) = max(0.0_realkind, cusnow(jl))
          stsnow(jl) = max(0.0_realkind, stsnow(jl))

          !     probicemod = probice+(1-probice)*(snow/totprecip)
          !     
          zwrk3 = zwrk1 + (1.0_realkind-zwrk1) *(storcu(jl) * cusnow(jl) &
               / (prcpcu(jl) + 1.e-12_realkind) + (1.0_realkind-storcu(jl)) * stsnow(jl) &
               / (prcpst(jl) + 1.e-12_realkind))
          vterm(jl) = (1._realkind- zwrk3)*hvterm + zwrk3*hvsnow
          !     
          !     take ice precip rate from above into account in
          !     the b-f effect
          !     
          zwrk17=(storcu(jl)*cusnow(jl)+(1._realkind-storcu(jl))*stsnow(jl))/snoref
          cbfsno = asnow * zwrk17 ** bsnow
          cbfsno = cbfsno / (1.0_realkind+cbfsno)

          !     valud of cbfsno to reproduce old results
          !     
          !     cbfsno = 1.
          !     
          !     f-bf = prbicemod * (1-prbice)*dwi
          !     
          zwrk6=max(0.0_realkind,bfeff(jintt)+(zwrk3-zwrk1)* (1._realkind-zwrk1) * zwrk2)
          !     
          !     zwrk1  - fmr, so probability of ice multiplied by the differenc
          !     between saturated pressure vapor between water and ice
          !     xfprim - m(tau-1)/b-mr, exponent of the precipitation with
          !     coalescence effect equation (eq. 16 in osz94)
          !     zwrk55 - precipitation with the coalescence effect (eq. 16 in o
          !     zwrk2  - fco, function of the coalescence effect (eq. 15 in s93
          !     xfrcoa - fco as in osz94
          !     hfmrx  - actual mr in the microphysics scale
          !     
          zwrk1 = hmroft(jintt)
          if(jwmask(jl,jk) == 1) then
             !     -------- stratiform case -------
             zwrk2 = hmrst
             zwrk51 = coalst
             cstocu(jl) = hcst
          else

             !     -------- convective case -------
             zwrk2 = hmrcu
             zwrk51 = coalcu
             cstocu(jl) = hccu
             if(jk == kht(jl) .and. jwanv(jl) == 1) then
                !     ----- with anvil ------
                zwrk2 = hmrst
                zwrk51 = coalst
                cstocu(jl) = hcst
             endif
          endif
          !     
          hfmrx(jl) = zwrk2
          xfprim(jl) = cw(jl,jk)/((totcov(jl,jk)+1.e-2_realkind) * zwrk1 * zwrk2)
          zwrk51 = zwrk51 * dpf(jl,jk) * totcov(jl,jk) * &
               xfprim(jl)*(1.0_realkind- exp(-xfprim(jl)**2.0_realkind))
          zwrk51 = max(zwrk51,0.0_realkind)
          !     
          zwrk2 = 1.0_realkind + coales*sqrt((prcpcu(jl) + prcpst(jl)) &
               /(totcov(jl,jk)+1.e-2_realkind) + zwrk51)
          !     
          zwrk2  = zwrk2 + cbfeff*zwrk6*cbfsno
          xfrcoa(jl) = zwrk2
          !     
          hfmrx(jl) = hfmrx(jl) * zwrk1 / zwrk2
          !     
          wrk1(jl)=zwrk1
          wrk2(jl)=zwrk2
          wrk3(jl)=zwrk3
          !     
1605   enddo

       !     special treatment for t<238 k
       !     
       do 1610 jl = kstart,kstop
          !     
          if (t(jl,jk) <= 238.0_realkind) then
             !     
             xfrcoa(jl) = xfrcoa(jl) * (1._realkind + (238.0_realkind - t(jl,jk))/2.0_realkind)
             if (t(jl,jk) <= 230.0_realkind) xfrcoa(jl) = wrk2(jl) * 5.0_realkind
             !     
          endif
          !     
1610   enddo
       !     

       !     calculate the fixed part of equ and normalize by b * mr
       !     
       !     wrk1   - independent term in eq.35 in osz94
       !     wrk2   - x final value
       !     zwrk52 - first guess of x/2
       !     zwrk25 - factor x in the newton-raphson solution (eq. 34 in osz94)
       !     zwrk1  - maximum value of x allowed to avoid too small exponentials
       !     zxp    - exp(-x**2)
       !     zxhj   - 1+dt.c0.(1-exp(-x**2))
       !     zxf    - residual value of eq. 35 in osz94
       !     zwrk11 - intermediate computation
       !     
       !     $dir scalar
       do 1628 jl = kstart,kstop
          !     
          wrk1(jl) = ((2.0_realkind*cw(jl,jk)+dtime*(hdcwad(jl,jk) &
               + cndsat(jl,jk)))/(2.0_realkind*(totcov(jl,jk)+1.e-2_realkind)*hfmrx(jl)))
          !     
          !     the conversion rate times 2*dt (dt.rhs in eq 33 in osz94)
          !     
          xfrcoa(jl) = 0.5_realkind * cstocu(jl) * xfrcoa(jl) * dtime
          !     
          !     first guess is the normalized m(t-dt)
          !     
          zwrk25 = cw(jl,jk)/((totcov(jl,jk)+1.e-2_realkind)*hfmrx(jl))

          !     to make m(t+dt)ge.0, the final solution
          !     of ym has to be ym>=m(t-dt)/(2.*b*mr)
          !     
          zwrk52 = 0.5_realkind*zwrk25
          !     
          if(jwqcon(jl) == 1 .or.(jk>=kht(jl) .and. jk<=khb(jl) &
               .and. jwmask(jl,jk)==0)) then
             !     
             !     newton-raphson solution
             !     -----------------------
             !     
             !     ---------------- iteration 1 -----------------------
             zwrk1 = min(zwrk25, 5.0_realkind)
             zxp = exp(-zwrk1*zwrk1)
             zxhj = 1.0_realkind + xfrcoa(jl)*(1.0_realkind-zxp)
             zxf = zwrk25*zxhj - wrk1(jl)
             zwrk11 = zxhj + 2.0_realkind*xfrcoa(jl)*zwrk25*zwrk25*zxp
             zwrk25 = max((zwrk25-zxf/zwrk11), zwrk52)

             !     ---------------- iteration 2 -----------------------
             zwrk1 = min(zwrk25, 5.0_realkind)
             zxp = exp(-zwrk1*zwrk1)
             zxhj = 1.0_realkind + xfrcoa(jl)*(1.0_realkind-zxp)
             zxf = zwrk25*zxhj - wrk1(jl)
             zwrk11 = zxhj + 2.0_realkind*xfrcoa(jl)*zwrk25*zwrk25*zxp
             zwrk25 = max((zwrk25-zxf/zwrk11), zwrk52)

             !     ---------------- iteration 3 -----------------------
             zwrk1 = min(zwrk25, 5.0_realkind)
             zxp = exp(-zwrk1*zwrk1)
             zxhj = 1.0_realkind + xfrcoa(jl)*(1.0_realkind-zxp)
             zxf = zwrk25*zxhj - wrk1(jl)
             zwrk11 = zxhj + 2.0_realkind*xfrcoa(jl)*zwrk25*zwrk25*zxp
             zwrk25 = max((zwrk25-zxf/zwrk11), zwrk52)

             !     ---------------- iteration 4 -----------------------
             zwrk1 = min(zwrk25, 5.0_realkind)
             zxp = exp(-zwrk1*zwrk1)
             zxhj = 1.0_realkind + xfrcoa(jl)*(1.0_realkind-zxp)
             zxf = zwrk25*zxhj - wrk1(jl)
             zwrk11 = zxhj + 2.0_realkind*xfrcoa(jl)*zwrk25*zwrk25*zxp
             zwrk25 = max((zwrk25-zxf/zwrk11), zwrk52)

             !     ---------------- iteration 5 -----------------------
             zwrk1 = min(zwrk25, 5.0_realkind)
             zxp = exp(-zwrk1*zwrk1)
             zxhj = 1.0_realkind + xfrcoa(jl)*(1.0_realkind-zxp)
             zxf = zwrk25*zxhj - wrk1(jl)
             zwrk11 = zxhj + 2.0_realkind*xfrcoa(jl)*zwrk25*zwrk25*zxp
             zwrk25 = max((zwrk25-zxf/zwrk11), zwrk52)

             !     ---------------- end of iteration -----------------
             !     
             wrk2(jl) = zwrk25
          endif
          !     
1628   enddo
       !     
       !     wrk2   - cw(tau+1)
       !     wrk1   - cw(tau+1)
       !     hdcwst - tendency of cloud water due to stratiform condensation
       !     hstcov - stratiform cloud cover
       !     cwpath - cloud water vertically integrated
       !     prcpst - stratiform precipitation
       !     stsnow - stratiform snow
       !     cov2d  - max-random overlap in calculation of cloudcover
       !     which is used in calculation of evaporation of
       !     precipitation
       !     *****   note that the resulting cover then is (1-cov2d);
       !     hence cov2d=1 means cloudfree
       !     *****   note further that on output cov2d = mean cloud cover
       !     cstocu - used to be set equal to either the stratiform hcst, or
       !     the convective hccu
       !     
       do 1635 jl = kstart,kstop
          !     
          wrk1(jl) = 2.0_realkind*wrk2(jl)*(totcov(jl,jk)+1.e-2_realkind)*hfmrx(jl) - cw(jl,jk)
          wrk1(jl) = max(wrk1(jl),0.0_realkind)
          if(wrk1(jl)<=1.e-6_realkind*cw(jl,jk).or.wrk1(jl)<1.e-13_realkind) then
             wrk1(jl) = 0._realkind
             hstcov(jl,jk) = 0.0_realkind
             hcucov(jl,jk) = 0.0_realkind
             totcov(jl,jk) = 0.0_realkind
          endif
          cwpath(jl) = cwpath(jl) + dpf(jl,jk)/gravit * wrk1(jl)
          if(jk>=kht(jl) .and. jk<=khb(jl) .and. jwmask(jl,jk)==0)then
             !     
             hdcwcu(jl,jk)= (wrk1(jl)-cw(jl,jk))/dtime-hdcwad(jl,jk)

             !     rate of convective precipitation
             !     
             zwrk25 = dpf(jl,jk)/gravit * &
                  max((cndsat(jl,jk)-hdcwcu(jl,jk)),0.0_realkind)
             prcpcu(jl) = prcpcu(jl)+ zwrk25
             prcl(jl,jk)=prcpcu(jl)
             cusnow(jl) = min((cusnow(jl) + wrk3(jl)*zwrk25),prcpcu(jl))
             !     
             if(jk > 2) then
                cov2d(jl) = cov2d(jl) * (1.0_realkind - max(hcucov(jl,jk-1), &
                     hcucov(jl,jk))) / &
                     (1.0_realkind - min(hcucov(jl,jk-1),0.99_realkind))
                cov2d(jl)=min(cov2d(jl),1.0_realkind)
             endif
             !     
             wrk4(jl) = sqrt(1.0_realkind - cov2d(jl))
             !     
          else
             !     
             hdcwst(jl,jk)=((wrk1(jl)-cw(jl,jk))/dtime- &
                  hdcwad(jl,jk)) * real(jwqcon(jl),realkind)

             !     rate of stratiform precipitation
             !     
             wrk1(jl) = dpf(jl,jk)/gravit * &
                  max((cndsat(jl,jk)-hdcwst(jl,jk)), 0.0_realkind) * real(jwmask(jl,jk),realkind)
             prcpst(jl) = prcpst(jl)+ wrk1(jl)
             stsnow(jl) = min((stsnow(jl) + wrk3(jl)*wrk1(jl)),prcpst(jl))
             !     
             hdtst(jl,jk) = cndsat(jl,jk)*hldcp2(jl,jk)
             hdqst(jl,jk) = -cndsat(jl,jk)
             !     
             if(abs(cwpath(jl)) < 1.e-14_realkind) then

                !     cov2d max-random overlap, i.e., that cov2d is set=1
                !     
                cov2d(jl) = 1.0_realkind
                wrk4(jl) = 0.0_realkind
             endif
          endif
1635   enddo
       !     
       !     calculate possible melting
       !     
       !     wrk1   - minimum value in eq. 32 in osz94
       !     wrk2   - sqrt of piceout in eq. 32 in osz94
       !     zwrk34 - piceout in eq. 32 in osz94
       !     zwrk53 - amount of ice precipitation melted
       !     stsnow - new value of cusnow
       !     hdtst  - dtdt with the decreasing for melting
       !     
       do 1655 jl = kstart,kstop
          !     
          zwrk25 = storcu(jl) * cusnow(jl) + (1.0_realkind-storcu(jl))*stsnow(jl)
          !     
          if(t(jl,jk) > tmelt .and. zwrk25 > 0.0_realkind) then

             wrk1(jl) = min((dpf(jl,jk)/&
                  (gravit*wro(jl,jk) * hvsnow)),dtime)
             wrk2(jl) = sqrt(zwrk25) - &
                  0.5_realkind*hkmelt*(t(jl,jk)-tmelt) * sqrt(hvsnow) * &
                  wrk1(jl) * wrk4(jl)
             wrk3(jl) = (max(wrk2(jl),0.0_realkind))**2.0_realkind
             zwrk57 = zwrk25 - wrk3(jl)
             !     
             if(jwmask(jl,jk) == 1) then

                !     -------- stratiform case -------
                stsnow(jl) = wrk3(jl)
                hdtst(jl,jk) = hdtst(jl,jk) - hdldcp * zwrk57 * &
                     gravit / dpf(jl,jk)
                !     
             else
                !     -------- convective case -------
                if(jk >= kht(jl) .and. jwmask(jl,jk) == 0 ) then
                   cusnow(jl) = wrk3(jl)
                   dtdt(jl,jk) = dtdt(jl,jk) - hdldcp * zwrk57 * &
                        gravit / dpf(jl,jk)
                endif
                !     
             endif
             !     
          endif
          !     
1655   enddo

       !     accumulate evaporation from cloud water and precipitation
       !     
       do 1665 jl = kstart,kstop
          !     
          sevapr = sevapr + hevapr(jl,jk)*dpf(jl,jk)/gravit
          sevapc = sevapc + hevac(jl,jk)*dpf(jl,jk)/gravit
          prsi(jl,jk) = stsnow(jl)
          prsl(jl,jk)=prcpst(jl)-stsnow(jl)
          !     
1665   enddo
       !     
       !     


       !     the end of the main vertical loop
       !     
1670 enddo
    !     

    !     
    !     
    !-----------------------------------------------------------------------
    !     
    !     17)  final setting of tendencies
    !     and calculation of cov2d for output
    !     
    !------------------------------------------------------------
    !     
    !     
    !     put cov2d to 1. here, see below
    !     
    do 1705 jl = kstart,kstop
       !     
       cov2d(jl) = 1.0_realkind
       !     
1705 enddo
    !     
    !     return calculated accumulated tendencies to the main arrays
    !     
    !     dcwdt  - main array of cw tendency
    !     totcov - cloudiness
    !     dtdt   - total tendency of temperature
    !     dqdt   - total tendency of specific humidity
    !     dcwdt  - total tendency of cloud water
    !     
    do 1720 jk = 1,nlev
       do 1714 jl = kstart,kstop
          !     
          dcwdt(jl,jk) = hdcwcu(jl,jk) + hdcwst(jl,jk)
          dtdt(jl,jk) = dtdt(jl,jk)+hdtst(jl,jk)
          dqdt(jl,jk) = dqdt(jl,jk)+hdqst(jl,jk)
          !     
1714   enddo
       !     
       do 1715 jl = kstart,kstop
          if(abs(dcwdt(jl,jk))<1.e-16_realkind) dcwdt(jl,jk) = 0.0_realkind
1715   enddo
       !     
       if(jk > 1) then
          do 1718 jl = kstart,kstop
             cov2d (jl) = cov2d(jl) * (1.0_realkind - max(totcov(jl,jk-1), &
                  totcov(jl,jk))) /(1.0_realkind- min(totcov(jl,jk-1),0.99_realkind))
1718      enddo
       endif
1720 enddo

    !     arrange cov2d for output:
    !     
    do 1725 jl = kstart,kstop

       cov2d(jl) = 1.0_realkind - cov2d(jl)
       if(abs(cwpath(jl)) < 1.e-14_realkind) cov2d(jl) = 0.0_realkind

1725 enddo
    !     
    !-----------------------------------------------------------------------

    !     18)   quantities for ground parameterization in the hirlam-model
    !     
    !     loop for making draindt and dsnowdt, which are needed
    !     in ground parameterization in hirlam-model
    !     
    !-----------------------------------------------------------------------
    !     
    do 1800 jl = kstart,kstop
       draindt(jl)=draindt(jl)+prcpst(jl)+prcpcu(jl)-stsnow(jl)-cusnow(jl)
       dsnowdt(jl) = dsnowdt(jl) + stsnow(jl) + cusnow(jl)
       !     
1800 enddo
    !     
    return
  end subroutine stamic



end module condensmod
