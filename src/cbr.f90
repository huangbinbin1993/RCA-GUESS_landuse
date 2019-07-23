module cbr
  use confys
  use escom
  use decomp,only:realkind
  implicit none
  private
  real(kind=realkind),save::ct=1.0_realkind
  real(kind=realkind),save::cq=1.0_realkind
  real(kind=realkind),save::cte=3.0_realkind
  real(kind=realkind),save::cu=1.0_realkind       ! everything is put into the length scale formulation
  real(kind=realkind),save::ctvar=0.139_realkind
  real(kind=realkind),save::cqvar=0.139_realkind
  real(kind=realkind),save::cht1=0.069_realkind
  real(kind=realkind),save::cht2=0.069_realkind
  real(kind=realkind),save::cdis=1.0_realkind/(3.75_realkind*3.75_realkind) ! consistent with surface bndary condition

  public vcbr
contains


  subroutine vcbr(nhor,nlev,kstart,kstop,dtime, &
       t,q,u,v,cw,tke,gpot, &
       pf,ph,dpf,dph, &
       ts, ps,taux,tauy,senf,latf,ustar, &
       dtdt ,dqdt, dudt ,dvdt,dcwdt, dtkedt, &
       lesat,pblh,om_cpbl, &
       zfrpbl,zfrtop,zlevtop, &
       z0wind,gustest,gustlow,gustup,gustestuc,&
       zvarcu,wdzf)


    ! ---------------------------------------------------------------------
    !    Routine that performs the computation of the turbulent fluxes
    !    and gives the corresponding tendencies to the prognostic variables.
    !    Routine is based on VCBR (v1.2) by J. Cuxart and J. Calvo (see
    !    references) but modified with respect to the following points.
    !    CHANGE 1) New formulation length scale (see section 3)
    !    CHANGE 2) TKE defined on half levels ! This increases
    !    numerical precision considerably.
    !    To make the change in code relatively simple, we used
    !    tke(nlev) for the surface. Therefore, the grid index of tke
    !    (and tke related variables) is one lower than the usual half
    !    level numbering. Note this carefully !!!
    !    This makes the scheme a bit incompatible with the computations
    !    outside the tke scheme, although probably the errors introduced
    !    are (very) small.
    !    CHANGE 3) TKE equation more implicit. Buoyancy, shear and
    !    dissipation are now all semi-implicit.
    !    CHANGE 4) boundary condition momentum flux (explicitly prescribed
    !    momentum flux) is replaced by implicitly computed one by
    !    prescribing u = 0 (surface) and Km(surf)
    !    Km(surf) is reconstructed from ustar. This change solves
    !    instability problems as encountered (over mountains)
    !    in various NWP runs when the surface momentum flux is large.
    !    CHANGE 5) Solving the diffusion equation has been changed.
    !    The original version did not incorporate density changes,
    !    and leaks about 5-10 % of surface fluxes


    !-----------------------------------------------------------------------
    ! Description:
    !
    !  Main routine of the scheme. It performs the computation of
    !  the turbulent fluxes and gives the corresponding tendencies to
    !  the prognostic variables.
    !
    !  The present version of the scheme computes dry turbulence, though
    !  the scheme is prepared for an easy jump to a moist turbulence scheme.
    !
    !  References:
    !
    !  J.Cuxart, P.Bougeault and Jean-Luc Redelsperger, Jan 2000
    !  A turbulence scheme allowing for mesoscale and large-eddy simulations
    !  Quarterly Journal of the Royal Meteorological Society, pp 1
    !
    !  J.Calvo, J.Cuxart and P.Kallberg:
    !  Implementation of an E-L turbulence scheme in the HIRLAM model.
    !  Submitted to Monthly Weather Review.
    !
    !  De Rooy, 2000: Experiences with the HIRLAM including the CBR scheme.
    !  Hirlam Newsletter No. 35
    !
    !  Lenderink, De Rooy, 2000: A robust mixing length formulation for a
    !  TKE-l turbulence scheme. Hirlam newsletter No. 36
    !
    !-----------------------------------------------------------------------
    !

    !    Input variables: surface fluxes, prognostic variables
    !    Output variables: tendencies of the prognostic variables
    !    Local variables: turbulent fluxes, Prandtl and Schmidt numbers,
    !                     mixing length, working variables



    !
    !     !   basic grid point resolution parameters
    implicit none

    integer nhor,nlev,kstart,kstop



    !-----------------------------------------------------------------------
    !
    !     ! input variables
    !
    real(kind=realkind),intent(in)::t(nhor,nlev), &     ! temperature input
         q(nhor,nlev),   &   ! moisture input
         cw(nhor,nlev),   &  ! cloud water input
         u(nhor,nlev),     & ! u wind input
         v(nhor,nlev),     & ! v wind input
         gpot(nhor,nlev),  & ! virtual potential above ground
         ts(nhor),         & ! surface temperature
         ps(nhor),         & ! surface temperature
         ustar(nhor),      & ! surface friction velocity
         taux(nhor),       & ! x surface stress ((kg/(m s2))
         tauy(nhor),       & ! y surface stress ((kg/(m s2))
         senf(nhor),       & ! surface heat flux (w/m2)
         latf(nhor),       & ! surface latent heat flux (w/m2)
         pblh(nhor),       & ! boundary layer height
         pf(nhor,nlev),    & ! pressure at model levels
         ph(nhor,nlev+1),  & ! pressure at model interface
         dpf(nhor,nlev),   & ! delta-p between model levels
         dph(nhor,nlev+1), & ! delta-p at half levels between full leve
!cau110525
!         along(nhor),      &
!         coslat(nhor),      &
!         sinlat(nhor),      &
!         frland(nhor),      &
         dtime              ! 2 delta-t

    real(kind=realkind),intent(inout)::tke(nhor,nlev) ! TKE input

    !
    !     ! input and output 1d variables
    !
!    real(kind=realkind),intent(inout)::stddif(nlev),stdif(nlev)

    !
    !     ! output variables
    !
    !         dtdtacc(nhor,nlev),& ! potential temperature after vertical
    !                              diffusion. on output tendency for
    !                              actual temperature!
    
    real(kind=realkind),intent(inout):: dudt(nhor,nlev), &   ! u-wind after vertical diffusion
         dvdt(nhor,nlev),  & ! v-wind after vertical diffusion
                                !                               on output tendencies for wind
         dtdt(nhor,nlev),  & ! potential temperature after vertical
                                !                              diffusion. on output tendency for
                                !                              actual temperature!
         dqdt(nhor,nlev), &  ! moisture after vertical diffusion
                                !                               on output tendency for specific hum.
         dtkedt(nhor,nlev), & ! TKE tendency
         dcwdt(nhor,nlev)   ! cw  tendency
    !gjkf
    !
    !     ! work space
    !
    integer i,k,kk
    !
    real(kind=realkind) wthm(nhor, nlev+1), &   ! potential temperature
         wsenf(nhor),      &  ! -senf
         wcflx(nhor),      &  ! -surface moisture flux (kg/m2/s)
         wths(nhor),       &  ! potential surface temperature
         wzh(nhor,nlev+1), &  ! z half levels
         wzf(nhor,nlev+1), &  ! z full levels
         wdzf(nhor,nlev),  &  ! delta-z between half levels
         wdzh(nhor,nlev+1),&  ! delta-z between full levels
         ahtf(nhor,nlev),  &  ! vertical interpolation vectors
         bhtf(nhor,nlev),  &  !              "
         afth(nhor,nlev+1),&  !              "
         bfth(nhor,nlev+1)   !              "
    !
    !     Local variables explicitly related to the CBR scheme
    !
    real(kind=realkind) phi3(nhor,nlev+1), &  ! stability function
         keff(nhor,nlev+1),  & ! basic eddy coefficient sqrt(e)*L
         thl(nhor,nlev),     & ! liquid potential temperature
         totw(nhor,nlev),    & ! total water
         qt(nhor,nlev+1),    &   ! total water (geert, new)
         thliq(nhor,nlev+1), &   ! liquid wat. pot. temp (geert, new)
         qsat(nhor,nlev), &
         jcoef(nhor,nlev), &
         tkefx(nhor,nlev),&
         tlfx(nhor,nlev+1),&
         twfx(nhor,nlev+1),&
         buofx(nhor,nlev+1),&
         trtke(nhor,nlev),&
         source(nhor,nlev), &
         tpvir(nhor,nlev+1), &
         zdeltpv(nhor,nlev), & ! delta-virtual pot. T between full lev.
         zasc(nhor,nlev),zcsc(nhor,nlev), &
         zamom(nhor,nlev),zcmom(nhor,nlev) &
         ,zatke(nhor,nlev),zctke(nhor,nlev)
    real(kind=realkind) zgbuoy(nhor,nlev+1), zgri(nhor,nlev+1), zcloud(nhor,nlev) &
         ,zmixdwh(nhor,nlev+1), zmixuph(nhor,nlev+1) &
         ,zmixdwm(nhor,nlev+1), zmixupm(nhor,nlev+1) &
         ,zmixquadh(nhor,nlev+1), zmixquadm(nhor,nlev+1) &
         ,zmixquadd(nhor,nlev+1), zri(nhor,nlev+1) &
         ,wdzfp(nhor,nlev)    & ! delta-pstar between half levels
         ,wdzhp(nhor,nlev+1)  & ! delta-pstar between full levels
         ,keffp(nhor,nlev+1)   ! basic eddy coefficient sqrt(e)*L
                                !abc  
                                !gjsplit
    real(kind=realkind) zgbuoy_clr(nhor,nlev+1),zgbuoy_cld(nhor,nlev+1)
    real(kind=realkind) zgri_clr(nhor,nlev+1),zgri_cld(nhor,nlev+1)
    real(kind=realkind) zmixquadh_clr(nhor,nlev+1),zmixquadh_cld(nhor,nlev+1)
    real(kind=realkind) zmixquadm_clr(nhor,nlev+1),zmixquadm_cld(nhor,nlev+1)
    real(kind=realkind) zmixquadd_clr(nhor,nlev+1),zmixquadd_cld(nhor,nlev+1)
    !gjsplit


    real(kind=realkind) ztkesq, zzb, zdisl, zktest, zimplmom, dtemp,&
         zdu, zldrt,zldct, zadry , zamoist, zbdry , zbmoist,&
         x, znex, delta, ccr, zathl, zbqt
    !GL
    !gj
    real(kind=realkind) zvarcu(nhor,nlev),st2d(nhor)
    real(kind=realkind) fice(nhor,nlev),latsub,rvair,cpv
    real(kind=realkind) afunc(nhor,nlev),bfunc(nhor,nlev)
    real(kind=realkind) tliq(nhor,nlev),qsliq(nhor,nlev)
    real(kind=realkind) zwrk1,zwrk2,zwrk3,zwrk4,zwrk5,tc
    real(kind=realkind) zwrk22,zwrk33
    real(kind=realkind) hliq(nhor,nlev)
    real(kind=realkind) qlst(nhor,nlev)
    real(kind=realkind) zwrk6,zwrk7,ess,essl
    logical lesat
    !gjgb2001
    integer corr,kinv(nhor),p300(nhor)
    real(kind=realkind) rho(nhor,nlev)
    real(kind=realkind) thlv(nhor,nlev),wthv(nhor,nlev)
    real(kind=realkind) buoy_inv(nhor)
    real(kind=realkind) entfac(nhor)
    real(kind=realkind) wentm(nhor)
    real(kind=realkind) enhf(nhor,nlev),test1(nhor,nlev)
    real(kind=realkind) ffunc(nhor,nlev),enhf2(nhor,nlev),varf(nhor,nlev)
    real(kind=realkind) zpote(nhor), zinte(nhor), zlwork(nhor)

    real(kind=realkind) thlp(nhor,nlev),totwp(nhor,nlev),up(nhor,nlev), &
         vp(nhor,nlev),tkep(nhor,nlev)

    real(kind=realkind) zwork1(nhor,nlev),zwork2(nhor,nlev+1),zwork(nhor), &
         zbeta(nhor),zrho(nhor),zi(nhor),buofxs(nhor)
    real(kind=realkind) cwpart, qpart, qsatpart, tpvpart,etheta, emoist, &
         wstar, lobukhov(nhor), psat, &
         zginv,zkappa,zp0,zlati,zlatw,zlatin
    !gj
    real(kind=realkind) om_cpbl(nhor),zfrtop(nhor),tpv_pbl(nhor), &
         dppbl(nhor),z0wind(nhor),gustest(nhor), &
         gustlow(nhor),gustup(nhor),gustestuc(nhor)
    integer zfrpbl(nhor),zlevtop(nhor),zfree(nhor)
    real(kind=realkind) qtnew(nhor,nlev+1)
    real(kind=realkind) tpp
    !gj060805 Countergradient term
    real(kind=realkind) zdt,zdq,zdfft(nhor),zdffq(nhor)
    real(kind=realkind) estemp(nhor),qstemp(nhor),gamma(nhor),  &
         zcpair,zidqs,dqsdt,zqmin


    !gj
!cau110525
!    real alat

    !       --------------------
    !    1. Preliminary settings
    !       --------------------
    !
    !     define local constants
    !

    zginv=1._realkind/gravit
    zkappa=rair/cpair
    zp0=1.0e5_realkind
    zlati=1.0_realkind/(latvap+latice)
    zlatw=1.0_realkind/latvap
    !gjgb2001
    latsub=latvap+latice
    rvair=461.0_realkind
    cpv=4.0_realkind*rvair
    !gjgb2001
    !gj060805....Deardorff type counter-gradient flux terms
    zdt=-5.0_realkind
    zdq=-5.0_realkind
!cgj280611
    zqmin=1.0e-9_realkind
    !gj060805
    !
    !     Level closest to the ground
    !     (sign of fluxes changed to have positive upwards)
    !
    do i=kstart,kstop
       !
       !gj   Variables needed for free convective velocity in KFCUMULUS
       !gj   These are output from VCBR and not needed in this routine.
       om_cpbl(i)=0.0_realkind
       zfrtop(i)=0.0_realkind
       tpv_pbl(i)=0.0_realkind
       dppbl(i)=0.0_realkind
       zfrpbl(i)=0
       zlevtop(i)=nlev+1
       zfree(i)=1
       !gj
       if(ts(i)<tmelt)then
          zlatin=zlati
       else
          zlatin=zlatw
       endif
       zrho(i)=pf(i,nlev)/rair/t(i,nlev)
       wsenf(i)=-senf(i)/(zrho(i)*cpair)
       wcflx(i)=-latf(i)*zlatin/zrho(i)
       wths(i)=ts(i)*(zp0/ps(i))**zkappa
    enddo
    !
    !     geometric height(wzf), theta and Tv(zwork1) at full levels
    !
    do k=1,nlev
       do i=kstart,kstop
          zcloud(i,k)=0.0_realkind
          wzf(i,k)=gpot(i,k)*zginv
          zwork1(i,k)=t(i,k)* (1.0_realkind + 0.61_realkind*q(i,k))
          wthm(i,k)= t(i,k)*(zp0/pf(i,k))**zkappa
          !gjgb2001
          wthv(i,k)=zwork1(i,k)*(zp0/pf(i,k))**zkappa   !virt pot. temp.
          !gjgb2001
          !gj define fraction of cw that is assumed to be frozen=fice
          tc=t(i,k) - tmelt
          fice(i,k)=max(0._realkind,min(-tc*0.05_realkind,1.0_realkind))
       enddo
    enddo

    do i=kstart,kstop
       zbeta(i)=gravit/wthm(i,nlev)
       wzf(i,nlev+1)=0.0_realkind
    enddo
    !
    !    geometric height (wzh) of the half levels
    !
    !      zwork2 contains the geopotential of the intermediate levels
    do i=kstart,kstop
       zwork2(i,nlev+1)=0.0_realkind
       do k=nlev,2,-1
          zwork2(i,k)=zwork2(i,k+1)+rair*zwork1(i,k) &
               *log(ph(i,k+1)/ph(i,k))
       enddo
       zwork2(i,1)=zwork2(i,2)+(gpot(i,1)-gpot(i,2))/2.0_realkind
    enddo

    do k=1,nlev+1
       do i=kstart,kstop
          wzh(i,k)=zwork2(i,k)*zginv
       enddo
    enddo
    !
    !     thickness between full and half levels (at half and full resp.)
    !
    do k=1,nlev-1
       do i=kstart,kstop
          wdzf(i,k)=wzh(i,k)-wzh(i,k+1)
          wdzh(i,k+1)=wzf(i,k)-wzf(i,k+1)

          wdzfp(i,k)   = ( ph(i,k+1) - ph (i,k) ) *zginv
          wdzhp(i,k+1) = ( pf(i,k+1) - pf (i,k) ) *zginv

       enddo
    enddo


    do i=kstart,kstop
       wdzh(i,1) = 0.0_realkind
       wdzf(i,nlev)=wzh(i,nlev)-wzh(i,nlev+1)
       wdzh(i,nlev+1)=wzf(i,nlev)

       wdzfp(i,nlev) = (ph(i,nlev+1)-ph(i,nlev))*zginv
       wdzhp(i,nlev+1)=(ph(i,nlev+1)-pf(i,nlev))*zginv
       wdzhp(i,1) = 0.0_realkind

    enddo


    do k=2,nlev
       do i=kstart,kstop
          bhtf(i,k)=      (wzh(i,k)-wzf(i,k)) / wdzf(i,k)
          ahtf(i,k)= 1.0_realkind -  bhtf(i,k)
          afth(i,k)=      (wzh(i,k)-wzf(i,k)) / wdzh(i,k)
          bfth(i,k)= 1.0_realkind -  afth(i,k)
          !gjgb2001
          if(ps(i)-pf(i,k)>30000.0_realkind)p300(i)=k  !k 300hpa above ground
          !search depth for KINV
          !gjgb2001
       enddo
    enddo

    do i=kstart,kstop
       afth(i,nlev+1)= 0.0_realkind
       bfth(i,nlev+1)= 1.0_realkind
       bhtf(i,1)=     (wzh(i,1)-wzf(i,1)) / wdzf(i,1)
       ahtf(i,1)=1.0_realkind -  bhtf(i,1)
    enddo
    !
    !

    !     TKE at the ground level (similarity formula)
    !     --------------------------------------------

    !
    !       tke(i,nlev+1)=F(surf fluxes)
    !       e*=Au*+Bw*, A=3.75,B=0 cas estable;
    !                   A=3.75+(-Z/L)**2/3, B=0.2 inestable
    !       L=-(theta0*ustar^3)/(k*g*buofxs)
    !       algoritme:
    !                  Zi= height where e < 0.01
    !                  w*=w*(zi,wth_surf)
    !                  L=L(zi,u*,w*)
    !      Note that cdis is dependent on A, by 
    !      approx. cdis = 1/A^2

    !     the following is a very simple way of computing boundary layer top
    !
    do i=kstart,kstop
       zi(i)=1.0_realkind
       do k=nlev,1,-1
          if (tke(i,k)<0.01_realkind) goto 130
          zi(i)=wzh(i,k)
       enddo
130    continue
    enddo

    do i=kstart,kstop
       buofxs(i)=wsenf(i)*(1.0_realkind+epsilo*q(i,nlev)) + &
            epsilo*wthm(i,nlev)*wcflx(i)
    enddo

    do i=kstart,kstop
       if ( buofxs(i)>=0.0_realkind) then
          wstar=(  zbeta(i)*zi(i)*buofxs(i) ) **(1.0_realkind/3.0_realkind)
          lobukhov(i)=-wthm(i,nlev)*(ustar(i)**3.0_realkind)/(0.4_realkind*gravit*buofxs(i))
          tke(i,nlev)=(3.75_realkind +((-1._realkind)*wzf(i,nlev)/ &
               lobukhov(i))**(2.0_realkind/3.0_realkind)) &
               *(ustar(i)**2_realkind)+0.2_realkind*(wstar**2.0_realkind)
       else
          tke(i,nlev)=3.75_realkind*(ustar(i)**2.0_realkind)
       endif
    enddo

    !
    !     TKE positivity control
    !     ----------------------
    !
    !     To prevent any occurrence of negative values produce by previous
    !     manipulation of the field (for instance, through advection or
    !     horizontal diffusion)
    !
    do k=1,nlev
       do i=kstart,kstop
          tke(i,k)=max( 1.0E-4_realkind,tke(i,k) )
       enddo
    enddo



    !         -------------------------------------------------------------
    !     2.1 Construction of the conservative variables (Theta_l and qtot)
    !         -------------------------------------------------------------
    !
    !     Here the conservative variables for the moist non-precipitating
    !     processes are defined.
    !     At present, the lines making Theta_l and q_tot differnt form Theta
    !     and qv are in comments, waiting for the planned activation of the
    !     moist scheme.

    !     continuation lines active to work with dry conservative var.
!cgj260611 compare thl to that in 7.3
    do k=1,nlev
       do i=kstart,kstop
          totw(i,k)=q(i,k) + cw(i,k)
          qt(i,k) = totw(i,k)
!cgjthliq
          thl(i,k)=wthm(i,k) &
               -(1.0_realkind-fice(i,k))*cw(i,k)*(latvap/ (cpair+qt(i,k)*cpv)) &
               *((zp0/pf(i,k))**zkappa) &
               -(fice(i,k))*cw(i,k)*(latsub/ (cpair+qt(i,k)*cpv)) &
               *((zp0/pf(i,k))**zkappa)
          thliq(i,k)=thl(i,k)
          !gjstatcld
          !gj     Define Liquid water Temperature, qsat(tliq), Latent heat
          !gj     term as function of ice/water ratio and constants a & b
          !gj     in Chaboureau & Bechtold JAS 2002...can use thliq also
          !gj
          tliq(i,k)=t(i,k) &
               -((1.0_realkind-fice(i,k))*cw(i,k)*(latvap/(cpair+qt(i,k)*cpv))) &
               -(fice(i,k)*cw(i,k)*(latsub/(cpair+qt(i,k)*cpv)))
!cau110518
!          if ( tliq(i,k) .lt. 100. .or. tliq(i,k) .gt. 400. )then
!             alat=acos(coslat(i))*180./3.141592654
!             if(sinlat(i).lt.0.)alat=-alat
!             write(6,6123) i,k,along(i),alat,frland(i),tliq(i,k)
!6123         format(' cbr i,k,tliq ',2I4,4F10.2)
!             stop ' cbr '
!          endif
          !gj
          if(lesat)then
             if(cw(i,k)<=0.0_realkind)then
                essl=esatw(tliq(i,k))
             else
                essl=(zcloud(i,k)*esat(tliq(i,k)))+((1.0_realkind-zcloud(i,k))*esatw(tliq(i,k))) 
             endif
          else
             essl=esat(tliq(i,k))
          endif
          !gj

          hliq(i,k)=((cpair+qt(i,k)*cpv)*tliq(i,k))+ &
               ((1._realkind+qt(i,k))*gpot(i,k)) &
               -((1.0_realkind-fice(i,k))*cw(i,k)*latvap) &
               -(fice(i,k)*cw(i,k)*latsub)
          !
          !gjgb2001
          thlv(i,k)=wthv(i,k) &
               -(1.0_realkind-fice(i,k))*cw(i,k)*(latvap/ (cpair+qt(i,k)*cpv)) &
               *((zp0/pf(i,k))**zkappa) &
               -(fice(i,k))*cw(i,k)*(latsub/ (cpair+qt(i,k)*cpv)) &
               *((zp0/pf(i,k))**zkappa)
          !gjgb2001
          qsliq(i,k)=epsilo*essl/(pf(i,k)-0.378_realkind*essl)
          !
          zwrk1=(1.0_realkind-fice(i,k))*latvap + fice(i,k)*latsub
          zwrk2=zwrk1*qsliq(i,k)/(rvair*tliq(i,k)*tliq(i,k))
          zwrk3=1.0_realkind+(zwrk1*zwrk2/(cpair+qt(i,k)*cpv))
          afunc(i,k)=1.0_realkind/zwrk3
          bfunc(i,k)=afunc(i,k)*zwrk2
          !gjstatcld
       enddo
    enddo


    !         -----------------------------
    !     2.2 virtual potential temperature
    !         -----------------------------
    !     For its use in the mixing length and in the buoyancy flux in
    !     the TKE equation.


    !     virtual potential temperature  at full levels
    !     simple formulation of saturation value. This does not have to be
    !     very precise as all the function are not strongly dependent on qsat
    !     Preferably qsat should be "imported" from outside to be consistent
    !     with other parts of HIRLAM
    !     used now is Tetens formula
    do k=1,nlev
       do i=kstart,kstop
          !         rvap/rair=1.61
          tpvir(i,k)=wthm(i,k)*(1.0_realkind+1.61_realkind*q(i,k))/(1.0_realkind+totw(i,k))
          !gj
          if(lesat)then
             if(cw(i,k)<=0.0_realkind)then
                ess=esatw(t(i,k))
             else
                ess=(zcloud(i,k)*esat(t(i,k)))+((1.0_realkind-zcloud(i,k))*esatw(t(i,k))) 
             endif
          else
             ess=esat(t(i,k))
          endif
          !gj
          qsat(i,k)=epsilo*ess/(pf(i,k)-0.378_realkind*ess)
       enddo
    enddo

    do i=kstart,kstop
       !        extrapolate some fields to surface at nlev+1
       wthm(i,nlev+1)= wthm(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - wthm(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
       tpvir(i,nlev+1)= tpvir(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - tpvir(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
       qt(i,nlev+1)= qt(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - qt(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
       thliq(i,nlev+1)= thliq(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - thliq(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
    enddo
    !gj
    !gj   Extract top of PBL from stability and save height and level
    !gj   calculate mean virtual temperature of PBL..needed for
    !gj   later calculation of free convective scale velocity used
    !gj   in KFCUMULUS.f...code not needed in VCBR.f
    !gj
    do k=nlev,1,-1
       do i=kstart,kstop
          if(zfree(i)==1)then
             if((tpvir(i,k+1)-tpvir(i,k))<0.0_realkind)then
                zfree(i)=0
             else
                zfrpbl(i)=1
                if(k==nlev)then
                   zfrtop(i)=0.0_realkind
                   zlevtop(i)=k
                else
                   zfrtop(i)=wzh(i,k)
                   zlevtop(i)=k
                endif
             endif
          endif
       enddo
    enddo
    !gj
    do k=nlev,1,-1
       do i=kstart,kstop
          if(k>=zlevtop(i))then
             if(zlevtop(i)==nlev)then
                tpv_pbl(i)=tpvir(i,k)*dpf(i,k)
                dppbl(i)=dpf(i,k)
             else
                tpv_pbl(i)=tpv_pbl(i)+(tpvir(i,k)*dpf(i,k))
                dppbl(i)=dppbl(i)+dpf(i,k)
             endif
          endif
       enddo
    enddo
    !gj
    do i=kstart,kstop
       if(zfrtop(i)>0.0_realkind)then
          tpv_pbl(i)=tpv_pbl(i)/dppbl(i)
       else
          tpv_pbl(i)=0.0_realkind
       endif
    enddo
    !gj   finished with PBL mean quantity calculations
    !gj   finished with deriving variables needed by KFCUMULUS.f
    !gj
    !gjstatcld
    !     compute brunt vaisala frequency (stability) zgbuoy from
    !     virtual potential temperature gradient
    !     This is a preliminary value for use in calculating the Cuijpers 
    !     constants and thus a moist mixing length. We will use
    !     a simple estimate of the cloud fraction and derive
    !     preliminary "moist" mixing lengths
    !     We will use this and moist Ri No to get preliminary
    !     "moist" mixing lengths. These mixing lengths will then be
    !     used to caluclate the variance of saturation deficit as
    !     described in Bechtold etal MWR 1992 and Bechtold & Chabereau
    !     JAS 2002. Once we have an estimate of the variance, we can
    !     calculate a new, turbulence based cloud fraction. This cloud
    !     fraction is then used to recalculate the Cuijpers constants
    !     and a new estimate of moist stability and moist mixing 
    !     lengths are derived. Once we have this final moist mixing lengths
    !     A final cloud fraction is calculated, using the saturation
    !     deficit variance derived from the finalised moist mixing lengths.
    !     At this point a consistent estimate of the liquid water portion
    !     of the total water content is also made. Presently this variable
    !     does not yet communicate with the rest of the code. The final
    !     cloud cover does and is the active cloud cover for the model 
    !     timestep. Cloud is therefore calculated after advection and
    !     directly before radiation. For greatest consistency use of the
    !     diagnosed cloud water would also be desirable at this point.
    !gjstatcld
    do k = 1,nlev
       do i=kstart,kstop
          if(abs(zcloud(i,k))<1.e-14_realkind)then
             znex = 0.7_realkind
             delta    = 0.02_realkind*qsat(i,k)
             !gj          x =(qt(i,k) - qsat(i,k)) / delta         
             x =afunc(i,k)*(qt(i,k) - qsat(i,k)) / delta         
             ccr = 0.5_realkind + 0.36_realkind * atan( 1.55_realkind * x )          ! bechtold formula
             zcloud(i,k) = min(max(ccr,0.0_realkind),1.0_realkind)
          endif
       enddo
    enddo
    do k = 1, nlev 
       do i=kstart,kstop
          zldrt = latvap / (rair * t(i,k))
          zldct = latvap / (cpair * t(i,k)) 
          zadry = 1.0_realkind + 0.61_realkind*qt(i,k)
          zamoist = (1.0_realkind - qt(i,k) + 1.61_realkind*qsat(i,k)*(1.0_realkind + epsilo*zldrt)) &
               / (1.0_realkind + epsilo * zldrt * zldct * qsat (i,k) ) 
          zbdry =   0.61_realkind * wthm(i,k)
          zbmoist = (zldct * zamoist - 1.0_realkind) *  wthm(i,k)

          !        moist stability

          zathl = (1.0_realkind-zcloud(i,k))*zadry + zcloud(i,k)*zamoist
          zbqt = (1.0_realkind-zcloud(i,k))*zbdry + zcloud(i,k)*zbmoist

          zgbuoy(i,k+1) = gravit/tpvir(i,k)*  &
               (zathl * (thliq(i,k)-thliq(i,k+1)) & 
               + zbqt * (qt(i,k)-qt(i,k+1))        )/wdzh(i,k+1)
          !gj
       enddo
    enddo
    do i=kstart,kstop
       zgbuoy(i,1)=zgbuoy(i,2)
    enddo
    !
    !     add computation Richardson number

    do k=1,nlev-1
       do i=kstart,kstop
          zdu = max(0.1_realkind,((u(i,k)-u(i,k+1))**2.0_realkind + (v(i,k)-v(i,k+1))**2.0_realkind)) &
               /(wdzh(i,k+1)**2.0_realkind)
          !         max as a security for noshear limit; no phys. relevance
          zgri(i,k+1) = zgbuoy(i,k+1) / zdu
       enddo
    enddo
    do i=kstart,kstop
       zgri(i,nlev+1) = zgbuoy(i,nlev+1)*(wzf(i,nlev)**2.0_realkind) &
            / max( u(i,nlev)**2.0_realkind + v(i,nlev)**2.0_realkind ,0.1_realkind) 
       !       max as a security for noshear limit; no phys. relevance
       zgri(i,1) = zgri(i,2)
    enddo
    !gjgb2001
    do i=kstart,kstop
       zgbuoy(i,1) = zgbuoy(i,2)   
       kinv(i)=-999
    enddo
    !ps050920
    !
    !gj  Diagnostic calculation of wind speed gustiness 
    !gj  Not needed for VCBR.f routine.
    !    Calculate gust wind speed
    call  wind_gust(nhor,nlev,kstart,kstop, &
         u,v,tke,pblh, &
         wzf,wzh,wdzh,zgbuoy,z0wind, &
         gustest,gustlow,gustup,gustestuc)

!    call entrain(nhor,nlev,kstart,kstop, &
!         t,cw,wthv,thlv,qt,pf,ph,wthm, &
!         zwork1,qsat,p300,wdzh, &
!         buoy_inv,entfac,kinv)
    !gj
    !gj   Get the new mixing lengths
    !gj
    call mixlen(nhor,nlev,kstart,kstop,&
         wdzf,wzh,tke,zgbuoy,zgri,&
         zmixquadm,zmixquadh,zmixquadd, &
         pf,phi3)
    !
    do i=kstart,kstop
       zmixquadm(i,nlev+1)=0.01_realkind
       zmixquadh(i,nlev+1)=0.01_realkind
       zmixquadd(i,nlev+1)=0.01_realkind
       zmixquadm(i,1)=0.01_realkind
       zmixquadh(i,1)=0.01_realkind
       zmixquadd(i,1)=0.01_realkind
    enddo

    !gj
    !gj   Now use dry mixing lengths to get variance quantities
    !gj   and calculate statistical cloud fraction & water content
    !gj
    call statcld(nhor,nlev,kstart,kstop,wzf,wzh, &
         qt,hliq,qsliq,wdzh,qsat,cw,zi, &
         zmixquadh,zmixquadm,afunc,bfunc, &
         zvarcu,ffunc,varf, &
         zcloud,qlst)
    !gj
    !gj    Perhaps one should calculate qlst at this point also
    !gj    and use it in a recalculation of qt for use in the
    !gj    Cuijpers constants below.??
    !gj    Calculate Cuijpers consatnts and moist zgbuoy & zgri
    !gj
    !gj  Check use of qlst instead of cw or qtnew=totw
    !cg  also check deardorff terms and whether dqdt=totwp-totw works
    !cj  perhaps also think about statcld formulation 
    do k = 1, nlev 
       do i=kstart,kstop
          qtnew(i,k)=q(i,k)+qlst(i,k)
          thliq(i,k)=wthm(i,k) &
!cgjnoqlst               -(1.0_realkind-fice(i,k))*qlst(i,k)*(latvap/ (cpair+q(i,k)*ccpq)) &
               -(1.0_realkind-fice(i,k))*cw(i,k)*(latvap/ (cpair+q(i,k)*ccpq)) &
               *((zp0/pf(i,k))**zkappa) &
!cgjnoqlst               -(fice(i,k))*qlst(i,k)*(latsub/ (cpair+q(i,k)*ccpq)) &
               -(fice(i,k))*cw(i,k)*(latsub/ (cpair+q(i,k)*ccpq)) &
               *((zp0/pf(i,k))**zkappa)
       enddo
    enddo
    do i=kstart,kstop
       !        extrapolate some fields to surface at nlev+1
       qtnew(i,nlev+1)= qtnew(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - qtnew(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
       thliq(i,nlev+1)= thliq(i,nlev)  * wzf(i,nlev-1) / wdzh(i,nlev) &
            - thliq(i,nlev-1)*   wzf(i,nlev) / wdzh(i,nlev)
    enddo
    !gj
    do k = 1, nlev 
       do i=kstart,kstop
          zldrt = latvap / (rair * t(i,k))
          zldct = latvap / (cpair * t(i,k)) 
          zadry = 1.0_realkind + 0.61_realkind*qt(i,k)
!cgjnoqlst          zamoist = (1.0_realkind-qtnew(i,k)+1.61_realkind*qsat(i,k)*(1.0_realkind + epsilo*zldrt)) &
          zamoist = (1.0_realkind-qt(i,k)+1.61_realkind*qsat(i,k)*(1.0_realkind + epsilo*zldrt)) &
               / (1.0_realkind + epsilo * zldrt * zldct * qsat (i,k) ) 
          zbdry =   0.61_realkind * wthm(i,k)
          zbmoist = (zldct * zamoist - 1.0_realkind) *  wthm(i,k)

          !        moist stability

          zathl = (1.0_realkind-zcloud(i,k))*zadry + zcloud(i,k)*zamoist
          zbqt = (1.0_realkind-zcloud(i,k))*zbdry + zcloud(i,k)*zbmoist

          zgbuoy(i,k+1) = gravit/tpvir(i,k)*  &
               (zathl * (thliq(i,k)-thliq(i,k+1)) & 
               + zbqt * (qt(i,k)-qt(i,k+1))        )/wdzh(i,k+1)
       enddo
    enddo
    do i=kstart,kstop
       zgbuoy(i,1)=zgbuoy(i,2)
    enddo
    !gj
    do k=1,nlev-1
       do i=kstart,kstop
          zdu = max(0.1_realkind,((u(i,k)-u(i,k+1))**2.0_realkind + (v(i,k)-v(i,k+1))**2.0_realkind)) &
               /(wdzh(i,k+1)**2.0_realkind)
          !         max as a security for noshear limit; no phys. relevance
          zgri(i,k+1) = zgbuoy(i,k+1) / zdu
       enddo
    enddo
    !
    do i=kstart,kstop
       zgri(i,nlev+1) = zgbuoy(i,nlev+1)*(wzf(i,nlev)**2.0_realkind) &
            / max( u(i,nlev)**2.0_realkind + v(i,nlev)**2.0_realkind ,0.1_realkind) 
       !       max as a security for noshear limit; no phys. relevance
       zgri(i,1) = zgri(i,2)
    enddo
    !gj
    !gj   Now have moist zgbuoy & zgri and dry mixing length defined cloud
    !gj   variance. Use these to define a moist mixing length
    !gj
    call mixlen(nhor,nlev,kstart,kstop, &
         wdzf,wzh,tke,zgbuoy,zgri, &
         zmixquadm,zmixquadh,zmixquadd, &
         pf,phi3)
    !
    do i=kstart,kstop
       zmixquadm(i,nlev+1)=0.01_realkind
       zmixquadh(i,nlev+1)=0.01_realkind
       zmixquadd(i,nlev+1)=0.01_realkind
       zmixquadm(i,1)=0.01_realkind
       zmixquadh(i,1)=0.01_realkind
       zmixquadd(i,1)=0.01_realkind
    enddo
    !gjsplit
    !gj

    !        ------------------------------------------------
    !     5. Turbulent fluxes and values of variables at t+1
    !        ------------------------------------------------
    !
    !    The diffusion equation is solved implicitly in time.
    !    The solver is not yet optimal; small errors may result
    !    from the balance between explicit processes and the
    !    implicit diffusion solver
    !    (see e.g., Lenderink and Holtslag, MWR, 2000, 244-258)
    !    In fact, stability should be computed from full time
    !    level fields, whereas the diffusion solver should use
    !    the updated (intermediate time level) fields.
    !    In the present code this cannot be done as one needs
    !    at least the tendencies due to the explicit processes
    !    or the full and intermediate time level fields.
    !    Errors are probably small for dry turbulence, but may
    !    significant for cloudy BL and high vertical resolution.


    !    The classical method of inverting a tridiagonal matrix is used.
    !    The same routine is called and what differs between variables is
    !    the coefficients of the matrix, passed as arguments to TRIDIAG.
    !    The sources (basically the surface fluxes) are introduced in
    !    the independent term of the System (the Y array). The output of
    !    TRIDIAG is the value of the variable as if it has been only forced
    !    by the turbulence (variables thlp, qtotp, cwp, up, vp). With these
    !    values the tendency is constructed in the paragraph 8.
    !
    zimplmom = 1.0_realkind      ! 1 surface momentum flux implicit
    ! 0 explicit

    !    implicit is used to solve instability problem resulting from
    !    high surface roughness areas (see Hirlam Newsletter 38, 70-73)
    !
    do k=2,nlev
       do i=kstart,kstop
          keff(i,k) = zmixquadm(i,k)*sqrt(tke(i,k-1))
          keff(i,k) = max(keff(i,k),1.e-2_realkind)
       enddo
    enddo
    !
    do i=kstart,kstop
       keff(i,1)=keff(i,2)
       keff(i,nlev+1) = zimplmom*ustar(i)**2.0_realkind*wzf(i,nlev)/ &
            (cu*max(0.0001_realkind,(u(i,nlev)**2.0_realkind + v(i,nlev)**2.0_realkind)**0.5_realkind))
       !       max as a security for noshear limit; no phys. relevance
    enddo

    !gl
    !     make keffp = rho **2 * keff

    do k=2,nlev
       do i=kstart,kstop
          rho(i,k) = 0.5_realkind*( pf(i,k)/rair/t(i,k) +  pf(i,k-1)/rair/t(i,k-1) )
          keffp(i,k) = keff(i,k) * rho(i,k) ** 2.0_realkind
       enddo
    enddo

    do i=kstart,kstop
       keffp(i,nlev+1) = keff(i,nlev+1) * zrho(i) ** 2.0_realkind
       keffp(i,1) = keffp(i,2)
    enddo
!entrain
!    do i=kstart,kstop
!       if(kinv(i)>0)then
!          kk=kinv(i)
!          wentm(i)=entfac(i)*tke(i,kk-1)*sqrt(tke(i,kk-1))/ &
!               (zmixquadm(i,kk)*buoy_inv(i))
!       endif
!    enddo
!entrain

    !
    !     Coefficients for the tridiagonals
    !     ---------------------------------

    do k=nlev,2,-1                           !
       do i=kstart,kstop
          kk=kinv(i)
          if(k==kk)then
             zamom(i,k)= -dtime*(wentm(i)*wdzhp(i,kk)*rho(i,kk))/ &
                  (wdzhp(i,kk)*wdzfp(i,kk))
             zasc(i,k)= -dtime*(wentm(i)*wdzhp(i,kk)*rho(i,kk))/ &
                  (wdzhp(i,kk)*wdzfp(i,kk))
             zatke(i,k)= -dtime* &
                  0.5_realkind*((wentm(i)*wdzhp(i,kk)*rho(i,kk)) &
                  +keffp(i,kk+1))/(wdzhp(i,min(nlev,kk+1))*wdzfp(i,kk))
          else
             zamom(i,k)= -dtime*keffp(i,k)/(wdzhp(i,k)*wdzfp(i,k))
             zasc(i,k)= -dtime*keffp(i,k)*phi3(i,k)/(wdzhp(i,k)*wdzfp(i,k))
             !gl
             zatke(i,k) = -dtime* &      ! defined on tke-levels !
                  0.5_realkind*( keffp(i,k) + keffp(i,k+1) ) / &
                  (wdzhp(i,min(nlev,k+1))*wdzfp(i,k))
          endif
       enddo
    enddo

    do k=nlev-1,1,-1
       do i=kstart,kstop
          kk=kinv(i)
          if(k==kk-1)then
             zcmom(i,k)= -dtime*(wentm(i)*wdzhp(i,kk)*rho(i,kk))/ &
                  (wdzhp(i,kk)*wdzfp(i,kk-1))
             zcsc(i,k)= -dtime*(wentm(i)*wdzhp(i,kk)*rho(i,kk))/ &
                  (wdzhp(i,kk)*wdzfp(i,kk-1))
             zctke(i,k)=-dtime* &
                  0.5_realkind*((wentm(i)*wdzhp(i,kk)*rho(i,kk)) &
                  +keffp(i,kk+1))/ &
                  (wdzhp(i,kk)*wdzfp(i,kk))
             !
          else
             !
             zcmom(i,k)= -dtime*keffp(i,k+1)/(wdzhp(i,k+1)*wdzfp(i,k))
             zcsc(i,k)= -dtime*keffp(i,k+1)*phi3(i,k+1)/(wdzhp(i,k+1)*wdzfp(i,k))
             !gl
             zctke(i,k) = -dtime* &   ! defined on tke-levels !
                  0.5_realkind*( keffp(i,k+1)+keffp(i,k+2) ) / &
                  (wdzhp(i,k+1)*wdzfp(i,k+1))
          endif
          !
          source(i,k)=0.0_realkind
       enddo
    enddo

    do i=kstart,kstop
       zamom(i,1)=0.0_realkind
       zasc(i,1)=0.0_realkind
       zcmom(i,nlev)= -dtime*keffp(i,nlev+1)/  &      ! important
            (wdzhp(i,nlev+1)*wdzfp(i,nlev))
       zcsc(i,nlev)=0.0_realkind
       zatke(i,1)=0.0_realkind
       zctke(i,nlev)=0.0_realkind
    enddo


    !     Source due to surface fluxes and splitted fields at t+1
    !     -------------------------------------------------------

    !     Theta_l

    do i=kstart,kstop
       source(i,nlev)= wsenf(i)/wdzf(i,nlev)
    enddo
    call tridiag(nhor,nlev,kstart,kstop, &
         thl,zasc,zcsc,dtime,ct,source,thlp )!, &
!         tkefx,trtke,zwork1,zwork)

    !     qtot
    
    do i=kstart,kstop
       source(i,nlev)= wcflx(i) / wdzf(i,nlev)
    enddo
    call tridiag(nhor,nlev,kstart,kstop, &
         totw,zasc,zcsc,dtime,cq,source,totwp) 

    !     u-component

    do i=kstart,kstop
       source(i,nlev)= - (1.0_realkind-zimplmom)* taux(i) / (zrho(i) * wdzf(i,nlev))
    enddo
    up = 99999.0_realkind

    call tridiag(nhor,nlev,kstart,kstop, &
         u,zamom,zcmom,dtime,cu,source,up) 


    !     v-component

    do i=kstart,kstop
       source(i,nlev)= - (1.0_realkind-zimplmom)*tauy(i) / (zrho(i) * wdzf(i,nlev))
    enddo
    call tridiag(nhor,nlev,kstart,kstop, &
         v,zamom,zcmom,dtime,cu,source,vp) !, &
!         tkefx,trtke,zwork1,zwork)


    !        -------------------------------
    !     6. Turbulent fluxes  (half levels)
    !        -------------------------------
    !
    ! Here the turbulent fluxes of heat and water are reconstructed at
    ! present for their use within the TKE equation (buoyancy term).
    ! They are consistently computed according to the values at t+1.
    !
    ! When the moist scheme will be installed in the reference system,
    ! in this section, the fluxes, variances and correlations of the
    ! conservatives variables will be computed for the option of the
    ! subgrid-condensation scheme.
    !
    do k=nlev,2,-1
       do i=kstart,kstop

          !         Conservative potential temperature flux
          tlfx(i,k)=-ct*phi3(i,k)*keff(i,k)*(thlp(i,k-1)-thlp(i,k))/wdzh(i,k)

          !         total water flux
          twfx(i,k)=-cq*phi3(i,k)*keff(i,k)*(totwp(i,k-1)-totwp(i,k))/wdzh(i,k)
       enddo
    enddo

    !
    do i=kstart,kstop
       tlfx(i,1)=0.0_realkind
       twfx(i,1)=0.0_realkind
       tlfx(i,nlev+1)=wsenf(i)
       twfx(i,nlev+1)=wcflx(i)
       !gj     Free convective scale velocity for use in KFCUMULUS.f
       !gj     Not needed in VCBR.f
       if(zfrtop(i)>0.0_realkind .and. buofxs(i)>0.0_realkind)then
          zwrk1=gravit*zfrtop(i)/tpv_pbl(i)
          zwrk2=zwrk1*buofxs(i)
          om_cpbl(i)=zwrk2**(1.0_realkind/3.0_realkind) !0.333_realkind
          zdfft(i)=zdt*buofxs(i)/(zfrtop(i)*om_cpbl(i))
          zdffq(i)=zdq*wcflx(i)/(zfrtop(i)*om_cpbl(i))
       else
          om_cpbl(i)=0.0_realkind
          zdfft(i)=0.0_realkind
          zdffq(i)=0.0_realkind
       endif
    enddo
    !
    !

    !        ----------------
    !     7. Evolution of TKE
    !        ----------------
    !
    !     computation shear production

    !gj   Make shear computation more implicit
    do k=nlev,2,-1
       do i=kstart,kstop
          zwork1(i,k-1)=cu*keff(i,k)*((up(i,k-1)-up(i,k))*(u(i,k-1)-u(i,k)) &
               +  (vp(i,k-1)-vp(i,k))*(v(i,k-1)-v(i,k)) ) /(wdzh(i,k)**2.0_realkind)
       enddo
    enddo
!!cgj260611 Geert sets zwork1(i,1) & zwork1(i,nlev) = 0. in hl 7.3

    !     buoyancy production

    do k=nlev-1,1,-1
       do i=kstart,kstop
          buofx(i,k) = -zgbuoy(i,k+1)*ct*keff(i,k+1)*phi3(i,k+1)
       enddo
    enddo

    !     dissipation

    do k=1,nlev-1
       do i=kstart,kstop
          zwork2(i,k)=- cdis*(tke(i,k)**(3.0_realkind/2.0_realkind))/zmixquadd(i,k+1)
       enddo
    enddo

    do i=kstart,kstop
       zwork2(i,nlev)=0.0_realkind
    enddo


    !gl implicit dissipation shear, and buoyancy.

    do k=1,nlev-1
       do i=kstart,kstop

          !         method as in echam4.
          !         basically an implicit (in time) equation for the 
          !         sqrt(tke) is solved. This can be done exactly 
          !         for buoy, shear, and dissipation (given the length scales)
          !         the provisional value tke is put into tkep

          ztkesq = tke(i,k)
          ztkesq  = sqrt(max(ztkesq,1.e-3_realkind))
          zzb =  (zwork1(i,k) + buofx(i,k))/ztkesq 
          zdisl=  zmixquadd(i,k+1)/cdis/dtime
          zktest= 1.0_realkind+(zzb*dtime + sqrt(tke(i,k))*2.0_realkind)/zdisl
          if (zktest<=1.0_realkind) then
             tkep(i,k) = 1e-4_realkind
          else
             tkep(i,k) = max(1e-4_realkind,(zdisl*(sqrt(zktest)-1.0_realkind))**2.0_realkind)
          endif
       enddo
    enddo
    !
    !     source term (buoy,shear and dissipation)
    do k=1,nlev-1
       do i=kstart,kstop
          !gj  Use implicit solver for TKE
          source(i,k) = (tkep(i,k)-tke(i,k))/dtime          ! implicit
          !gj         source(i,k)= buofx(i,k) + zwork1(i,k)+zwork2(i,k)  ! explicit
          !         explicit is better satisfying balance between 
          !         (shear/buoy/diss)  and diffusion, but might 
          !         be numerically somewhat less stable
       enddo
    enddo

    !     Activate if e(u*,w*) is desired
    do i=kstart,kstop
       source(i,nlev)= 0.0_realkind
    enddo

    !     coefficients now zatke,zctke
    call tridiag(nhor,nlev,kstart,kstop,&
         tke,zatke,zctke,dtime,cte,source,tkep)
!, &
!         tkefx,trtke,jcoef,zwork)

    !     positivity control
    do k=1,nlev
       do i=kstart,kstop
          tkep(i,k)=max(1.e-4_realkind,tkep(i,k))
       enddo
    enddo
    !
    !
    !
    !

    !        ---------------------
    !     8. Turbulence tendencies
    !        ---------------------
    ! The turbulence tendencies for the prognostic variables are computed
    ! from the output values of the tridiagonal and the t-1 values.
    ! 

    do k= nlev,1,-1
!cgj
!cgj280611 Get gamma to derive dqsdt etc from vqsatd for use in diagnosing a 
!cgj280611 dcwdt term from the total water tendency
!cgj
       call vqsatd (nhor,kstart,kstop,lesat, &
            t(1,k),pf(1,k),cw(1,k),zcloud(1,k),estemp,qstemp,gamma)
!cgj280611

!cgj010911 Testing old cw treatment in cbr...to get correct back, comment cgj010911_comment
! deleted lines between cgjrca35 and uncomment lines cgj010911
!
       do i=kstart,kstop
!          dcwdt(i,k)=0.0_realkind	!cgj010911_comment
          dudt(i,k) = (up(i,k) - u(i,k))/dtime
          dvdt(i,k) = (vp(i,k) - v(i,k))/dtime
          dtkedt(i,k)=(tkep(i,k)-tke(i,k))/dtime 
!cgjca35
!          dtdt(i,k) =(thlp(i,k) - thl(i,k))/dtime			&
!           +((latvap/(cpair+totwp(i,k)*cpv))				&
!           *((zp0/pf(i,k))**zkappa)*(1.-fice(i,k))*dcwdt(i,k))		&
!           +((latsub/(cpair+totwp(i,k)*cpv))				&
!           *((zp0/pf(i,k))**zkappa)*fice(i,k)*dcwdt(i,k))	
!          dtdt(i,k)=dtdt(i,k)*(pf(i,k)/zp0)**zkappa
!cgjca35
          dtdt(i,k) =(thlp(i,k) - thl(i,k))/dtime 
          dqdt(i,k) =(totwp(i,k) - totw(i,k))/dtime  

          zcpair = cpair*(1._realkind + ccpq*q(i,k))
          dqsdt = cpair*gamma(i)/latvap
          zidqs = 1._realkind/(1._realkind + (latvap + fice(i,k)*latice)  &
                  * dqsdt / zcpair )
          dcwdt(i,k) = zidqs*zcloud(i,k)*(dqdt(i,k)-dqsdt*dtdt(i,k))
!cgj280611 check for possible -dcwdt overshoots
          if((cw(i,k)/dtime + dcwdt(i,k)) < 0._realkind)          &
                              dcwdt(i,k)=-cw(i,k)/dtime
!cgj280611
          dqdt(i,k) = dqdt(i,k) - dcwdt(i,k)
!cgj280611 check for possible q undershoots
          if((q(i,k)/dtime + dqdt(i,k)) < 0._realkind)                   &
                             dqdt(i,k)=(-q(i,k)+zqmin)/dtime
!cgj280611
          dtdt(i,k)=(pf(i,k)/zp0)**zkappa*( dtdt(i,k) +           &
                    ((latvap+fice(i,k)*latice)/zcpair)*dcwdt(i,k) )

          !gj060805
          !gj060805 Add in countergradient term...after Deardorff
          !gj060805
          if(wzh(i,k)<zfrtop(i))then
             dtdt(i,k)=dtdt(i,k)-((keffp(i,k+1)*zdfft(i))/(dtime*wdzfp(i,k)))
             dqdt(i,k)=dqdt(i,k)-((keffp(i,k+1)*zdffq(i))/(dtime*wdzfp(i,k)))
          endif
          !gj060805
       enddo
    enddo
    !



    return        
  end subroutine vcbr


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! ======================================================================

  subroutine tridiag(nhor,nlev,kstart,kstop, &
       pvarm,pa,pc,ptime,ct,psource,pvarp) 

    implicit none
    !
    !     Solves a tridiagonal system of equations. Conceived for a variable
    !     'pvarm' at a full level
    !
    !gl
    !     the standard version of this scheme is designed for a flux
    !     boundary condition at the surface level, but
    !     this routine also is suitable for solving a system with
    !     a zero (value) boundary condition (like for the horizontal wind)
    !     In that case pc(nlev) has to be nonzero, and contains K
    !     in surface layer.

    !     Input variables

    integer,intent(in)::nhor,nlev,kstart,kstop
    real(kind=realkind),intent(in):: pvarm(nhor,nlev)
    real(kind=realkind),intent(in):: pa(nhor,nlev)
    real(kind=realkind),intent(in):: pc(nhor,nlev)
    real(kind=realkind),intent(in):: ptime,ct
    real(kind=realkind),intent(in):: psource(nhor,nlev)

    !     Output variables
    real(kind=realkind),intent(inout)::pvarp(nhor,nlev)

    !     Local variables
    integer:: k,i
    real(kind=realkind):: zb(nhor,nlev),zy(nhor,nlev),zom(nhor,nlev)
    real(kind=realkind):: bet(nhor)


    do k=nlev,1,-1
       do i=kstart,kstop
          zy(i,k)=pvarm(i,k) + ptime*psource(i,k)
          zb(i,k)=1.0_realkind-ct*(pa(i,k)+pc(i,k))
          zom(i,k)=0.0_realkind
       enddo
    enddo

    do i=kstart,kstop
       bet(i)=zb(i,1)
       pvarp(i,1)=zy(i,1)/bet(i)
    enddo

    do k=2,nlev
       do i=kstart,kstop
          zom(i,k)=ct*pc(i,k-1)/bet(i)
          bet(i)=zb(i,k)-ct*pa(i,k)*zom(i,k)
          pvarp(i,k)=(zy(i,k)-ct*pa(i,k)*pvarp(i,k-1))/bet(i)
       enddo
    enddo

    do k=nlev-1,1,-1
       do i=kstart,kstop
          pvarp(i,k)=pvarp(i,k)-zom(i,k+1)*pvarp(i,k+1)
       enddo
    enddo


    return
  end subroutine tridiag

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine mixlen(nhor,nlev,kstart,kstop,&
       wdzf,wzh,tke,zgbuoy,zgri,&
       zmixquadm,zmixquadh,zmixquadd,&
       pf,phi3)

    implicit none

    integer,intent(in)::nhor,nlev,kstart,kstop
    integer::i,k
    !   
    !     input
    !
    real(kind=realkind),intent(in):: wdzf(nhor,nlev),wzh(nhor,nlev+1),tke(nhor,nlev),&
         pf(nhor,nlev), zgbuoy(nhor,nlev+1),zgri(nhor,nlev+1)
    !   
    !     output
    !
    real(kind=realkind),intent(inout):: zmixquadm(nhor,nlev+1), zmixquadh(nhor,nlev+1),&
         zmixquadd(nhor,nlev+1), phi3(nhor,nlev+1)
    !
    !     internal variables
    !
    real(kind=realkind) zbuoythress,zprandt,zalpha,zkap,pi2,zneutr, &
         zfineutr,zfias,zbri,zbuoy_ful,zri_ful, &
         xarm,xarh,xm,xh,dz,x,dzm,dzh,pfunc
    real(kind=realkind) zmix_help,zave,zlengthrh,zlengthrm
    real(kind=realkind) zwz,zfunc,zentr,zentrh,zentrm,zlam,zmix2, &
         zmixh,zmixm,ztkesq, zmixh_stable,zmixm_stable,zmixch,zmixcm
    real(kind=realkind) zmixuph(nhor,nlev+1), zmixupm(nhor,nlev+1), zmixdwh(nhor,nlev+1),&
         zmixdwm(nhor,nlev+1)
    !
    !
    !        --------------------------------
    !     3. Computation of the mixing length
    !        --------------------------------
    !
    !     compute a simple quadratic length scale
    !     this is a preliminary code, which gives a "quadratic" lengthscale
    !     according to 1/zmixquad = 1/zmixup + 1/zmixdw
    !     zmixup: lengthscale from bottom (say ground)
    !        lower boundary: ground or lowest position unstable layer
    !     zmixdw: lengthscale from top (say inversion)
    !        upper boundary: "inversion"
    !     both zmixup and zmixdw increase for unstable and neutral
    !     conditions, and decrease for stable conditions.
    !     this is governed by zbuoythress; small zbuoythress ==>
    !     strong damping stable conditions (inversion), and
    !     increase for unstable conditions.
    !     in neutral conditions the lengthscale zmixup ~ zmixquad ~ zkap*z
    !     near the surface
    !     zprandt: governs difference between lengthscale for momentum
    !              and heat.
    !gj
    !gj   At this point zgbuoy & zgri are based on dry variables
    !gj   we will initially compute dry mixing lengths to get the
    !gj   variance in S from which we calculate zcloud. Once
    !gj   zcloud is calculated we can then calculate the Cuijpers
    !gj   constants and moist values of zgbuoy & zgri from which
    !gj   we will calculate moist mixing lengths and from these
    !gj   a new "moist variance" cloud fraction and liquid water
    !gj

    !     some initial constant and tuning parameters !

    zbuoythress = 0.1e-4_realkind  ! between 0.1e4 - 0.2e-4
    zprandt =  0.6_realkind        ! between 0.2 - 1
    zalpha = 1.0_realkind
    zkap = 0.4_realkind            
    pi2 = 2.0_realkind * atan(1.0_realkind)
    zneutr = 1.0_realkind/(3.75_realkind**0.5_realkind)

    !     tuning constants.
    !     zprandt determines lm/lh convective limit
    !     zfineutr scaling length scale near surface
    !     zfineutr = zneutr *zkap standard value 
    !     (derived in Lenderink 2002)
    !     zbri length scale dependency on Ri determined
    !     from linearization flux prof. relations
    !     zbri ~ 4-8.
    !     see details see Lenderink 2002)
    !     zfias convective limit Kh 
    !     zfias =  zneutr *zkap * 3. / zprandt  standard  

    zprandt =  0.6_realkind  

    zfineutr = zneutr *zkap * 1.3_realkind     ! test: to compensate 
    ! for transport of TKE in
    ! neutral conditions
    zfias =  zneutr *zkap * 3.0_realkind / zprandt        
    zbri = 8.0_realkind

    !     compute bottom up length scale

    do i=kstart,kstop               ! initialization at surface !
       zmixuph(i, nlev+1) = 0.0_realkind
       zmixupm(i, nlev+1) = 0.0_realkind
    enddo

    !     Loop over model levels          ! bottom-up loop !

    do k=nlev,2,-1
       do i=kstart,kstop
          zbuoy_ful = 0.5_realkind*(zgbuoy(i,k+1) + zgbuoy(i,k))
          zri_ful = 0.5_realkind*(zgri(i,k+1) + zgri(i,k)) 

          xarm = zbri * pi2 * zneutr * zkap /(zfias*zprandt-zfineutr)
          xarh = zbri * pi2 * zneutr * zkap /(zfias - zfineutr)

          xm = xarm * zri_ful
          xh = xarh * zri_ful
          dz = wdzf(i,k)

          x = zri_ful
          if (x>0._realkind) then
             dzm = zfineutr*dz - (zfias*zprandt - zfineutr) / pi2*dz*xm
             dzh = zfineutr*dz - (zfias -zfineutr) / pi2*dz*xh
          else
             dzm = zfineutr*dz - (zfias*zprandt - zfineutr) / pi2*dz*atan(xm)
             dzh = zfineutr*dz - (zfias - zfineutr) / pi2*dz*atan(xh)
          endif

          !         bottom up lengthscale

          zmixuph(i,k) =  zmixuph (i,k + 1) + dzh
          zmixuph(i,k) = max(zmixuph(i,k) , 0.01_realkind)
          zmixupm(i,k) =  zmixupm (i,k + 1) + dzm
          zmixupm(i,k) = max(zmixupm(i,k) , 0.01_realkind)

       enddo
    enddo
    !
    !     compute top down length scale
    !
    do i=kstart,kstop           ! initialize at the top of domain
       zmixdwh(i, 1) = 0.0_realkind
       zmixdwm(i, 1) = 0.0_realkind
    enddo

    do k = 2, nlev              ! vertical loop top, down
       do i=kstart,kstop
          !gj Slight correction to zri_full
          zbuoy_ful = 0.5_realkind*(zgbuoy(i,k) + zgbuoy(i,k+1))
          zri_ful = 0.5_realkind*(zgri(i,k) + zgri(i,k+1))
          !gj
          xarm = zbri * pi2 * zneutr * zkap /(zfias*zprandt-zfineutr)
          xarh = zbri * pi2 * zneutr * zkap /(zfias - zfineutr)

          xm = xarm * zri_ful
          xh = xarh * zri_ful
          dz = wdzf(i,k-1)

          x = zri_ful
          if (x>0._realkind) then
             dzm = zfineutr*dz - (zfias*zprandt - zfineutr) / pi2 *dz*xm
             dzh = zfineutr*dz - (zfias -zfineutr) / pi2*dz*xh
          else
             dzm = zfineutr*dz - (zfias*zprandt - zfineutr) / pi2*dz*atan(xm)
             dzh = zfineutr*dz - (zfias - zfineutr) / pi2*dz*atan(xh)
          endif


          !         limit downward length scale close to the surface
          !         to reduce influence of it on zmixquad (final length
          !         scale). This means that the zmixup will determine
          !         scaling close to the surface. 
          
          zmix_help = 75._realkind * exp ( - wzh(i,k) / 500.0_realkind) 
          zmixdwh(i,k) = zmixdwh (i,k-1)  +  dzh             
          zmixdwh(i,k) = max (zmixdwh(i,k),zmix_help) 
          zmixdwm(i,k) = zmixdwm (i,k-1)  +  dzm
          zmixdwm(i,k) = max (zmixdwm(i,k),zmix_help)    

          !         composed length scale of top-down ls and bottom-up ls

          zave = 1.0_realkind
          zlengthrh = 1.0_realkind / zmixdwh(i,k)**zave + 1.0_realkind / zmixuph(i,k)**zave
          zmixquadh(i,k) = 1.0_realkind / zlengthrh**(1.0_realkind/zave)
          zlengthrm = 1.0_realkind / zmixdwm(i,k)**zave + 1.0_realkind / zmixupm(i,k)**zave
          zmixquadm(i,k) = 1.0_realkind / zlengthrm**(1.0_realkind/zave)

       enddo
    enddo                                  ! end vertical loop


    !GL
    !    4. Computation lengthscale for stable conditions.
    !    At this moment basic convective lengthscale is in zmixquadm/h
    !    final lengthscale is put in zmixquadm/h
    !
    zmix_help = 0.0_realkind
    do k = 2, nlev            ! vertical loop top, down
       do i=kstart,kstop
!cgj040711          zentr = 0.1_realkind
!cgj040711          !         included to represent momentum mixing by waves !
!cgj040711          zentrh  = zentr
!cgj040711          !gj          zentrm  = max(zentr, zentr*zmix_help)
!cgj040711          zentrm  = zentr
!cgj040711          !gjteststablemods
!cgj040711
          zwz = wzh (i, k) / 1000._realkind
          zfunc = 2._realkind * zwz * zwz * zwz * exp ( - zwz)
          zfunc = min (1._realkind, max (zfunc, 0._realkind) )
          zmix_help = max (1._realkind, min ( (1._realkind +   &
                          (zfunc * zgri (i, k) * 3.8_realkind) ), 10._realkind) )
!cgj010812          zentr=0.18
          zentr=0.1
          zentrh=zentr
          zentrm=zentr*zmix_help
!cgj040711
          zlam = 75.0_realkind
          zmix2 = zlam / (1._realkind + zlam / (zneutr*zkap  * wzh(i,k))) 

          !         zmixh = max(zmixquadh(i,k),zmix2)                   
          !         zmixm = max(zmixquadm(i,k),zmix2)
          
          zave = 2._realkind    !  interpolation
          zmixh =  (zmix2**zave + zmixquadh(i,k)**zave)**(1.0_realkind/zave) 
          zmixm =  (zmix2**zave + zmixquadm(i,k)**zave)**(1.0_realkind/zave)


          ztkesq = tke(i,k-1)
          ztkesq  = sqrt(max(ztkesq,1e-4_realkind))

          zave = -2._realkind   ! interpolation
          if (zgbuoy(i,k)>0._realkind ) then                         ! stable
             zmixh_stable = zentrh * ztkesq / sqrt(zgbuoy(i,k))
             zmixch = (zmixh_stable**(zave) + zmixh**(zave))**(1.0_realkind/zave) 
             zmixm_stable = zentrm * ztkesq / sqrt(zgbuoy(i,k))
             zmixcm = (zmixm_stable**(zave) + zmixm**(zave))**(1.0_realkind/zave) 
          else
             zmixch = zmixh
             zmixcm = zmixm
          endif

          zmixquadh(i,k) = zmixch
          zmixquadm(i,k) = zmixcm
          zmixquadd(i,k) = zmixcm

          ! make it fit in the CBR scheme
          ! phi3 is defined on half levels consistent with the CBR scheme

          phi3(i,k) = zmixquadh(i,k) / zmixquadm(i,k)

       enddo
    enddo

    do i=kstart,kstop
       phi3(i,nlev+1)=phi3(i,nlev)
       phi3(i,1)=phi3(i,2)
       zmixquadm(i,nlev+1) = 0.01_realkind
       zmixquadm(i,1) = 0.01_realkind
       zmixquadd(i,nlev+1) = 0.01_realkind
       zmixquadd(i,1) = 0.01_realkind
       zmixquadh(i,nlev+1) = 0.01_realkind
       zmixquadh(i,1) = 0.01_realkind

    enddo
    !
    return
  end subroutine mixlen

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine wind_gust(nhor,nlev,kstart,kstop,&
       u,v,tke,pblh, &
       wzf,wzh,wdzh,zgbuoy,z0wind,&
       gustest,gustlow,gustup,gustestuc)
    !
    ! By Maria Nordstrom, 2005
    ! "Estimation of gusty winds in RCA"
    ! Master thesis, Department of Earth Sciences, Meteorology
    ! Uppsala University, ISSN 1650-6553 Nr 101
    ! http://www.met.uu.se/exjobb/exjobb2005/Maria_Nordstrom.pdf
    !
    implicit none
    integer,intent(in):: nhor,nlev,kstart,kstop
    !
    !     ! input variables
    !
    real(kind=realkind),intent(in)::u(nhor,nlev), &     ! u wind input
         v(nhor,nlev),  &    ! v wind input
         tke(nhor,nlev),&    ! TKE input
         pblh(nhor),    &    ! boundary layer height
         z0wind(nhor)       ! roughness length for wind
    !
    !     ! output variables
    !
    real(kind=realkind),intent(inout):: gustest(nhor), &    ! gust wind speed estimate
         gustlow(nhor), &    ! lower bound of gust wind speed estimate
         gustup(nhor), &     ! upper bound of gust wind speed estimate
         gustestuc(nhor)  ! uncorrected gust wind speed estimate



    !     ! work space
    real(kind=realkind)::gustlowuc(nhor)     ! uncorrected lower bound of gust wind speed estimate
    integer i,k
    real(kind=realkind) wzh(nhor,nlev+1), &  ! z half levels
         wzf(nhor,nlev+1), &  ! z full levels
         wdzh(nhor,nlev+1),&  ! delta-z between full levels
         zgbuoy(nhor,nlev+1),&
         fulltke(nhor,nlev), &
         ztke(nhor,nlev+1), &
         zfulltke(nhor,nlev+1)
    !
    real(kind=realkind) zanlev,zm,zk,vvar, inttke,intbuoy,zalpha,zpblh,zz0wind
    !
    !----------------------------------------------------------------------
    !
    !     define local constants
    !
    zanlev=10.0_realkind
    zalpha=1.4_realkind
    !
    do k=1,nlev
       do i=kstart,kstop
          ztke(i,k)=tke(i,k)
          zfulltke(i,k)=0.0_realkind
       enddo
    enddo
    do i=kstart,kstop
       ztke(i,nlev+1)=tke(i,nlev)
       zfulltke(i,nlev+1)=0.0_realkind
    enddo
    !     --------------------------------
    !     wind gust estimate
    !     --------------------------------
    !
    do i=kstart,kstop
       !
       zpblh=max(pblh(i),20.0_realkind)
       zz0wind=min(z0wind(i),1.0_realkind)
       !
       k=nlev
       inttke=0.0_realkind
       intbuoy=0.0_realkind
       do while(pblh(i)>=wzf(i,k))
          zk=(ztke(i,k+1)-ztke(i,k))/(wzh(i,k+1)-wzh(i,k))
          zm=ztke(i,k)-wzh(i,k)*zk
          fulltke(i,k)=wzf(i,k)*zk+zm
          zfulltke(i,k)=fulltke(i,k)
          zfulltke(i,nlev+1)=fulltke(i,nlev)
          inttke=1.0_realkind/(wzf(i,k)-zanlev)*0.5_realkind*(zfulltke(i,k)+zfulltke(i,k+1)) &
               *wdzh(i,k+1)+inttke
          intbuoy=0.5_realkind*(zgbuoy(i,k)*wdzh(i,k)+zgbuoy(i,k+1) &
               *wdzh(i,k+1))*wdzh(i,k+1)+intbuoy
          if(inttke>=intbuoy)then
             gustestuc(i)=max(sqrt(u(i,k)**2.0_realkind+v(i,k)**2.0_realkind),gustestuc(i))
             !ps050920 Include correction
             gustest(i)=gustestuc(i)*( log(10.0_realkind/zz0wind) / &
                  log(0.1_realkind*zpblh/zz0wind))*zalpha
          endif
          k=k-1
       enddo
    enddo       ! for i
    !
    !
    !     --------------------------------
    !     upper bound
    !     --------------------------------
    !
    do i=kstart,kstop
       !
       zpblh=max(pblh(i),20.0_realkind)
       zz0wind=min(z0wind(i),1.0_realkind)
       !
       k=nlev
       do while(pblh(i)>=wzf(i,k))
          gustup(i)=max(sqrt(u(i,k)**2.0_realkind+v(i,k)**2.0_realkind),gustup(i))
          !ps050920 Include correction
          gustup(i)=gustup(i)*( log(10.0_realkind/zz0wind) / &
               log(0.1_realkind*zpblh/zz0wind))*zalpha
          k=k-1
       enddo
    enddo
    !
    !     --------------------------------
    !     lower bound
    !     --------------------------------
    !
    do i=kstart,kstop
       !
       gustlowuc(i)=gustlow(i)
       !
       zpblh=max(pblh(i),20.0_realkind)
       zz0wind=min(z0wind(i),1.0_realkind)
       !
       k=nlev
       intbuoy=0.0_realkind
       do while(pblh(i)>=wzf(i,k))
          zk=(ztke(i,k+1)-ztke(i,k))/(wzh(i,k+1)-wzh(i,k))
          zm=ztke(i,k)-wzh(i,k)*zk
          fulltke(i,k)=wzf(i,k)*zk+zm
          vvar=2.5_realkind/11.0_realkind*fulltke(i,k)
          intbuoy=0.5_realkind*(zgbuoy(i,k)*wdzh(i,k)+zgbuoy(i,k+1) &
               *wdzh(i,k+1))*wdzh(i,k+1)+intbuoy
          if(vvar>=intbuoy)then
             gustlowuc(i)=max(sqrt(u(i,k)**2.0_realkind+v(i,k)**2.0_realkind),gustlowuc(i))
             !ps050920 Include correction
             gustlow(i)=gustlowuc(i)*( log(10.0_realkind/zz0wind) / &
                  log(0.1_realkind*zpblh/zz0wind))*zalpha
          endif
          k=k-1
       enddo
       !
       gustest(i)=max(gustest(i),gustlow(i))
       gustestuc(i)=max(gustestuc(i),gustlowuc(i))
       gustup(i)=max(gustup(i),gustest(i))
       !
    enddo
    !
    return
  end subroutine wind_gust

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine entrain(nhor,nlev,kstart,kstop, &
       t,cw,wthv,thlv,qt,pf,ph,wthm,&
       tv,qsat,p300,wdzh,&
       buoy_inv,entfac,kinv)

    implicit none
    !
    !     input
    !
    integer,intent(in)::nlev,nhor,kstart,kstop,p300(nhor)
    real(kind=realkind),intent(in):: t(nhor,nlev),&
         tv(nhor,nlev),&
         cw(nhor,nlev),&
         wthv(nhor,nlev),&
         thlv(nhor,nlev),&
         qt(nhor,nlev+1),&
         pf(nhor,nlev),&
         ph(nhor,nlev+1),&
         wthm(nhor,nlev),&
         wdzh(nhor,nlev+1),&
         qsat(nhor,nlev)
    !
    !     output
    !
    integer,intent(inout):: kinv(nhor)
    real(kind=realkind),intent(inout):: entfac(nhor)
    !
    integer:: i,k,kk
    !     internal arrays
    !
    real(kind=realkind):: tpvgr(nhor,nlev),&
         tpvgr_test,&
         tpvmax,&
         buoy_inv(nhor),&
         exn_inv,&
         thlv_inv(nhor), &
         qt_inv(nhor), &
         lam_inv,&
         eps_inv,zp0,zkappa, &
         tinv,esat_inv,qsat_inv,gamma_inv,beta_inv, &
         rr_inv,kd_inv,cone_inv, &
         zwrk1,zwrk2,zwrk3,zwrk4,zwrk5, &
         f_inv,evap_inv, &
         cwmin,cpl,a1,a2,omeps,rvair


    !     Find the half layer the inversion is located at and deterime
    !     the evaporative enhancement factor and entrainment factor needed
    !     for calculation of the entrainmant velocity in VCBR code

    a1=0.16_realkind
    a2=30.0_realkind
    cwmin=5.e-6_realkind
    tpvmax=0.0001_realkind
    zkappa=rair/cpair
    zp0=1.e5_realkind
    cpl=cpair/latvap
    omeps=1.0_realkind-epsilo
    rvair=461.0_realkind

    do k=nlev,2,-1
       do i=kstart,kstop
          if(k>=p300(i))then
             tpvgr_test=wthv(i,k-1)-wthv(i,k)
             tpvgr(i,k)=tpvgr_test/wdzh(i,k)
             if(tpvgr(i,k)>tpvmax .and. cw(i,k)>cwmin )then
                if(qt(i,k-1)<qsat(i,k-1) )then
                   tpvmax=tpvgr(i,k)
                   kinv(i)=k
                endif
             endif
          endif
       enddo
    enddo
    !cgjggb2001
    do i=kstart,kstop
       if(kinv(i)>0)then
          kk=kinv(i)
          buoy_inv(i)=gravit*(wthv(i,kk-1)-wthv(i,kk))/wthv(i,nlev)
          exn_inv=(zp0/ph(i,kk))**zkappa
          thlv_inv(i)=thlv(i,kk-1)-thlv(i,kk)
          qt_inv(i)=qt(i,kk-1)-qt(i,kk)
          lam_inv=(0.5_realkind*(cw(i,kk)+cw(i,kk-1)))/(cpl*thlv_inv(i)*exn_inv)
          eps_inv=cpl*(0.5_realkind*(wthm(i,kk-1)+wthm(i,kk)))
          tinv=0.5_realkind*(t(i,kk-1)+t(i,kk))
          esat_inv=esat(tinv)
          qsat_inv=epsilo*esat_inv/(ph(i,kk)-omeps*esat_inv)
          gamma_inv=(qsat_inv*latvap)/(tv(i,kk)*tv(i,kk)*cpl*rvair)
          beta_inv=(1.0_realkind+(1.608_realkind*(gamma_inv*eps_inv)))/(1.0_realkind+gamma_inv)
          rr_inv=-qt_inv(i)/(exn_inv*thlv_inv(i)*cpl)
          kd_inv=1._realkind-(exn_inv*eps_inv*0.608_realkind)
          cone_inv=gamma_inv/(1._realkind+gamma_inv)
          !cgjgb2001
          zwrk1=exn_inv*eps_inv
          zwrk2=(kd_inv-zwrk1)*lam_inv
          zwrk3=zwrk1-(beta_inv*kd_inv)
          zwrk4=rr_inv*zwrk3
          zwrk5=rr_inv*(1._realkind-(kd_inv*cone_inv))
          f_inv=(1._realkind-beta_inv-zwrk2-zwrk4)/((1._realkind-zwrk2)*(cone_inv+zwrk5))
          if(cw(i,kk)>0.0_realkind)then
             evap_inv=lam_inv*f_inv
          else
             evap_inv=0.0_realkind
          endif
          evap_inv=min(25.0_realkind,max(evap_inv,0.0_realkind))
          entfac(i)=a1*(1.0_realkind+a2*evap_inv)
       else
          entfac(i)=0.0_realkind
       endif
       !gjgb2001
    enddo
    !
    return
  end subroutine entrain

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  subroutine statcld(nhor,nlev,kstart,kstop,wzf,wzh, &
       qt,hliq,qsliq,wdzh,qsat,cw,pblh, &
       zmixquadh,zmixquadm,afunc,bfunc, &
       zvarcu,ffunc,varf, &
       zcloud,qlst)

    implicit none
    !
    !    Calculate a cloud fraction & liquid water content from
    !    the total water & moist static energy & turbulent mixing
    !    length, following the general ideas developed for statistical
    !    cloud schemes.
    integer::i,k,ind
    integer,intent(in)::nhor,nlev,kstart,kstop
    !
    !     input
    !      
    real(kind=realkind),intent(in)::qt(nhor,nlev+1),&
         hliq(nhor,nlev),&
         qsliq(nhor,nlev),&
         wdzh(nhor,nlev+1),&
         zmixquadh(nhor,nlev+1),&
         zmixquadm(nhor,nlev+1),&
         afunc(nhor,nlev),&
         bfunc(nhor,nlev),&
         qsat(nhor,nlev),&
         wzf(nhor,nlev),&
         wzh(nhor,nlev+1),&
         zvarcu(nhor,nlev),&
         cw(nhor,nlev),&
         pblh(nhor)
    !
    !     output
    !      
    real(kind=realkind),intent(inout)::zcloud(nhor,nlev),qlst(nhor,nlev)
    !
    !     internal
    !      
    real(kind=realkind) cpm,cpv,rvair,hldz,qtdz,zwrk1,zwrk2,zwrk3,zwrk4, &
         varf(nhor,nlev),x,x2,ccr
    real(kind=realkind) varh(nhor,nlev)
    real(kind=realkind) ffunc(nhor,nlev)
    real(kind=realkind) zeps
    real(kind=realkind) csp,zpbl1,zpbl2
    real(kind=realkind) zfix,zmlen

    zeps=1.e-12_realkind
    rvair=461.0_realkind
    cpv=4.0_realkind*rvair
    csp=0.12_realkind
    zfix=200.0_realkind
    !
    do k=2,nlev
       do i=kstart,kstop
          !gjgb
          zpbl1=max(2000.0_realkind,2.0_realkind*pblh(i))
          zpbl2=zpbl1+1000.0_realkind
          if(wzf(i,k)<zpbl1)then
             zmlen=zmixquadh(i,k)
          elseif(wzf(i,k)>=zpbl1 .and. wzf(i,k)<=zpbl2)then
             zmlen=(((wzf(i,k)-zpbl1)/1000.0_realkind)*zfix) &
                  +(((zpbl2-wzf(i,k))/1000.0_realkind)*zmixquadh(i,k))
          else
             zmlen=zfix
          endif
          !gj
          zcloud(i,k)=0.0_realkind
          ffunc(i,k)=0.0_realkind
          cpm=cpair+(0.5_realkind*(qt(i,k)+qt(i,k-1))*cpv)
          hldz=(hliq(i,k-1)-hliq(i,k))/wdzh(i,k)
          qtdz=(qt(i,k-1)-qt(i,k))/wdzh(i,k)
          zwrk1=afunc(i,k)*afunc(i,k)*qtdz*qtdz
          zwrk2=2._realkind*afunc(i,k)*bfunc(i,k)*hldz*qtdz/cpm
          zwrk3=(bfunc(i,k)*bfunc(i,k)*hldz*hldz)/(cpm*cpm)
          zwrk4=zwrk1-zwrk2+zwrk3
          zwrk4=max(0.0_realkind,zwrk4)
          !gj        varh(i,k)=sqrt(zmixquadh(i,k)*zmixquadh(i,k))*sqrt(zwrk4)
          varh(i,k)=sqrt(zmlen*zmlen)*sqrt(zwrk4)
          varh(i,k)=max(varh(i,k),zeps)
       enddo
    enddo
    !
    do i=kstart,kstop
       varh(i,1)=varh(i,2)
       ffunc(i,1)=ffunc(i,2)
       zcloud(i,1)=zcloud(i,2)
    enddo
    !gj
    do k=2,nlev
       do i=kstart,kstop
          if(k<nlev)then
             zwrk1=0.5_realkind*(zvarcu(i,k)+zvarcu(i,k+1))	!conv var
             zwrk2=0.5_realkind*(varh(i,k)+varh(i,k+1))	!turb var
             zwrk3=0.02_realkind*(qsat(i,k))		!meso var
             if(zwrk1>0.0_realkind)then
                varf(i,k)=zwrk1+zwrk2+zwrk3
             else
                varf(i,k)=zwrk2+zwrk3
             endif
          else
             varf(i,k)=varh(i,k)+(0.02_realkind*qsat(i,nlev))
             !gj          varf(i,k)=varh(i,k)	!use nlev value
          endif
          !gjtest
          !gj
          x=afunc(i,k)*(qt(i,k)-qsliq(i,k))/varf(i,k)
          ccr=0.5_realkind + 0.36_realkind* atan(1.55_realkind*x)
          zcloud(i,k)=min(max(ccr,0.0_realkind),1.0_realkind)
          if(x<0.0_realkind)then
             qlst(i,k)=varf(i,k)*(exp((2.4_realkind*x)-1.0_realkind))
             !gj          ffunc(i,k)=exp(-1.4*x)-1.       
          elseif(x>=0.0_realkind .and. x<=2.0_realkind)then
             qlst(i,k)=varf(i,k)*(exp(-1.0_realkind)+ 0.66_realkind*x + 0.086_realkind*x*x)
             !gj          ffunc(i,k)=0.
          elseif(x>2.0_realkind)then
             qlst(i,k)=x*varf(i,k)
             !gj          ffunc(i,k)=0.
          endif
          if(zcloud(i,k)<=0._realkind)qlst(i,k)=0.0_realkind
          !
          !      Now get func for skewness contribution to buoyancy flux
          !
          !gj          ffunc(i,k)=min(ffunc(i,k),200.)
          !
       enddo
    enddo

    do i=kstart,kstop
       zcloud(i,1)=0.0_realkind
       qlst(i,1)=0.0_realkind
    enddo
    !
    return
  end subroutine statcld

  subroutine vqsatd(plond,kstart,kstop,lesat, &
       t    ,p  ,cwat   ,cloud  ,    &
       es      ,qs      ,gam)
    !-----------------------------------------------------------------------
    !
    ! utility procedure to look up and return saturation vapor pressure from 
    ! precomputed table, calculate and return saturation specific humidity 
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dt).  the same
    ! function as qsatd, but operates on vectors of temperature and pressure
    !
    !----------------------------code history-------------------------------
    !
    ! original version:  j. hack
    ! standardized:      j. rosinski, june 1992
    !                    t. acker, march 1996
    ! reviewed:          j. hack, august 1992
    !
    use confys
    use escom
    implicit none
    !------------------------------arguments--------------------------------
    !
    ! input arguments
    !
    integer plond       ! vector length
    integer kstart      ! start point of vector
    integer kstop       ! end point of vector
    real(kind=realkind) t(plond)       ! temperature
    real(kind=realkind) p(plond)       ! pressure
    real(kind=realkind) cwat(plond)
    real(kind=realkind) cloud(plond)
    ! 
    ! output arguments
    !
    real(kind=realkind) es(plond)   ! saturation vapor pressure
    real(kind=realkind) qs(plond)   ! saturation specific humidity
    real(kind=realkind) gam(plond)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------local variables------------------------------
    !
    logical lflg   ! true if in temperature transition region
    !
    integer i      ! index for vector calculations
    !
    real(kind=realkind) omeps     ! 1. - epsilo
    real(kind=realkind) trinv     ! reciprocal of ttrice (transition range)
    real(kind=realkind) tc        ! temperature (in degrees c)
    real(kind=realkind) weight    ! weight for es transition from water to ice
    real(kind=realkind) hltalt    ! appropriately modified hlat for t derivatives  
    !
    real(kind=realkind) hlatsb    ! hlat weighted in transition region
    real(kind=realkind) hlatvp    ! hlat modified for t changes above 273.16
    real(kind=realkind) tterm     ! account for d(es)/dt in transition region
    real(kind=realkind) desdtm     ! d(es)/dt
    real(kind=realkind) cp        ! heat capacity for dry air
    real(kind=realkind) hlatv     ! latent heat of vaporization
    real(kind=realkind) hlatf     ! latent heat of freezing
    real(kind=realkind) ttrice    ! temperature treshold for ice
    real(kind=realkind) rgasv     ! gas constant for water vapour
    !vo replacing the routine gestbl.f
    real(kind=realkind) pcf(1:5)  ! diff. sat. vapour press. over water and ice
    logical lesat


    hlatv = latvap
    hlatf = latice
    rgasv = 461.50_realkind
    omeps = 1.0_realkind - epsilo
    !gj
    !gj  see comments in findsp
    !gj
    !gjorig      ttrice = 20.
    ttrice = 40._realkind
    !gj      ttrice = 0.
    cp = cpair
    ! pressure over ice for -ttrice < t < 0 (degrees c). note: polynomial
    ! is valid in the range -40 < t < 0 (degrees c).
    !
    !                  --- degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_realkind
    pcf(2) = -5.47288442819e+00_realkind
    pcf(3) = -3.67471858735e-01_realkind
    pcf(4) = -8.95963532403e-03_realkind
    pcf(5) = -7.78053686625e-05_realkind
    !vo ---

    do i=kstart,kstop
       !gj080304
       if(lesat)then
          if(cwat(i)<=0._realkind)then
             es(i)=esatw(t(i))
          else
             es(i)=(cloud(i)*esat(t(i)))+((1._realkind-cloud(i))*esatw(t(i)))
          endif
          !gj080304
       else
          es(i) = esat(t(i))
       endif
       !
       ! saturation specific humidity
       !
       qs(i) = epsilo*es(i)/(p(i) - omeps*es(i))
       !
       ! the following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = min(1.0_realkind,qs(i))
       !
       if (qs(i) < 0.0_realkind) then
          qs(i) = 1.0_realkind
          es(i) = p(i)
       end if
    end do
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_realkind
    if(abs(ttrice)<1.e-14_realkind) goto 10
    trinv = 1.0_realkind/ttrice
    do i=kstart,kstop
       !
       ! weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (c) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition 
       ! range from ice to water also accounting for change of hlatv with t 
       ! above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc>=-ttrice .and. tc<0.0_realkind)
       weight = min(-tc*trinv,1.0_realkind)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_realkind*tc
       if (t(i)<tmelt) then
          hltalt = hlatsb
       else
          hltalt = hlatvp
       end if
       if (lflg) then
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       else
          tterm = 0.0_realkind
       end if
       desdtm  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdtm/(cp*es(i)*(p(i) - omeps*es(i)))
       if(abs(qs(i)-1.0_realkind)<1.e-14_realkind) gam(i) = 0.0_realkind
    end do
    return
    !
    ! no icephs or water to ice transition
    !
10  continue
    do i=kstart,kstop
       !
       ! account for change of hlatv with t above 273.16 where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_realkind*(t(i)-tmelt)
       hlatsb = hlatv + hlatf
       if (t(i)<tmelt) then
          hltalt = hlatsb
       else
          hltalt = hlatvp
       end if
       desdtm  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdtm/(cp*es(i)*(p(i) - omeps*es(i)))
       if(abs(qs(i)-1.0_realkind)<1.e-14_realkind) gam(i) = 0.0_realkind
    end do
    !
    return
    !
  end subroutine vqsatd

end module cbr
