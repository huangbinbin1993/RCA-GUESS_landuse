module comdfb
  use confys
  use decomp

  implicit none
  private
  real(kind=realkind),public,save::ckmin= 0.01_realkind              ! reduced from 0.1 by bab
  real(kind=realkind),public,save::cms1     !computed
  real(kind=realkind),public,save::clog     !computed
  real(kind=realkind),public,save::cchrn
  real(kind=realkind),public,save::chrmax
  integer,public,save::npbl !computed
  integer,public,save::ntopfl = 1

  real(kind=realkind),save::betam=15.0_realkind
  real(kind=realkind),save::betas=5.0_realkind
  real(kind=realkind),save::betah=15.0_realkind
  real(kind=realkind),save::bmax =25.0_realkind               ! max values for b's in interpolation
  real(kind=realkind),save::bmin  =-25.0_realkind              ! min values for b's in interpolation
  real(kind=realkind),save::fak=8.5_realkind
  real(kind=realkind),save::g                   ! gravity
  real(kind=realkind),save::onet = 1.0_realkind/3.0_realkind
  real(kind=realkind),save::fakn=7.2_realkind
  real(kind=realkind),save::ricr                ! critical Richardson number
  real(kind=realkind),save::sffrac
  real(kind=realkind),save::vk                  ! von Karman's constant
  real(kind=realkind),save::ccon
  real(kind=realkind),save::binm
  real(kind=realkind),save::binh

  public inidif2, avdifnl
contains
  subroutine inidif2(nlev,hybf,hybh,cml2)

    !     define constants used in vertical diffusion


    implicit none
    integer:: nlev,k
    real(kind=realkind):: cml2(nlev+1),hybh(nlev+1),hybf(nlev)

    if(mype==0)print *,'The Holtslag scheme is used'

    ckmin = 0.01_realkind   !dfb           ! reduced from 0.1 by bab
    npbl = int(real(nlev,realkind)*0.5_realkind)+1 !dfb
    ntopfl = 1                   !dfb

    clog=hybf(nlev)              !dfb

    cms1=gravit/((hybh(nlev+1)-hybh(nlev)) *rair )*clog*carman**2 !dfb
    cchrn=0.032_realkind*(hybh(nlev+1)-hybh(nlev))*rair*288.15_realkind/(gravit**2_realkind*hybf(nlev)) !dfb
    clog=log(clog)   !dfb
    chrmax=0.50_realkind !dfb


    cml2(1) = 0.0_realkind !I/O
    do 10 k=2,nlev
       cml2(k) = 30.0_realkind**2_realkind
10  enddo
    cml2(nlev+1) = 0._realkind


    fak=8.5_realkind  !comdif2
    fakn=7.2_realkind !comdif2
    betam=15.0_realkind !comdif2
    betas=5.0_realkind !comdif2
    betah=15.0_realkind !comdif2
    sffrac=0.1_realkind !comdif2
    ricr=0.5_realkind !comdif2
    bmax=25.0_realkind !comdif2
    bmin=-25.0_realkind !comdif2

    vk= carman !comdif2
    g = gravit !comdif2
    onet = 1.0_realkind/3.0_realkind !comdif2
    ccon = fak*sffrac*vk !comdif2
    binm = betam*sffrac !comdif2
    binh = betah*sffrac !comdif2

    return
  end subroutine inidif2




  subroutine avdifnl(nhor, nlev, kstart, kstop, dtime,  t, &
       q, cw, totcov, u, v, gpot, heatv, obkl, pf, ph, dpf, dph, ts, ps, &
       t2m, q2m, taux, tauy, senf, latf, ustar, frland, &
       ustarg, &
       pblh, z0oro, dtdt, dqdt, dcwdt, dudt, dvdt )!, &
!       stddif, stdif)
    !     driver routine to compute vertical diffusion of momentum,
    !     moisture and potential temperature.
    !     coefficients are computed locally and passed to individual
    !     parameterizations mvdiff, qvdiff
    !     sg, smhi and bh, knmi november 91.
    !
    !     revised by bent hansen sass, december 1997
    !
    !     revised by niels woetmann nielsen, december-february 1999
    !     to improve treatment of momentum drag over land.

    implicit none  

    integer :: nhor, nlev, kstart, kstop  
    !     ! input variables
    ! temperature input

    real(kind=realkind) :: t(nhor, nlev), q(nhor, nlev), cw(nhor, nlev), totcov( &
         nhor, nlev), u(nhor, nlev), v(nhor, nlev), gpot(nhor, nlev), &
         ts(nhor), ps(nhor), t2m(nhor), q2m(nhor), ustar(nhor), &
         ustarg(nhor), z0oro(nhor), frland(nhor)
!    real(kind=realkind)::frice(nhor)
    real(kind=realkind)::taux(nhor), tauy(nhor), senf(nhor), latf(nhor), pf(nhor, nlev), &
         ph(nhor, nlev + 1), dpf(nhor, nlev), dph(nhor, nlev + 1), &
         pblh(nhor), dtime
    ! moisture input
    ! cloud water input
    ! cloud cover
    ! u wind input
    ! v wind input
    ! virtual potential above ground
    ! surface temperature
    ! surface temperature
    ! temperature at 2 meters height
    ! specific humidity at 2 meters height
    ! surface friction velocity
    ! friction velocity based on heat roughnes
    ! orography roughness
    ! fraction of land area
    ! fraction of ice area(not presently used)
    ! x surface stress(n)
    ! y surface stress(n)
    ! surface heat flux(w/m2)
    ! surface latent heat flux(w/m2)
    ! pressure at model levels
    ! pressure at model interface
    ! delta-p between model levels
    ! delta-p at half levels between full lvls
    ! planetary boundary layer height
    ! 2 delta-t
    !     output variables
    ! u-wind after vertical diffusion
    real(kind=realkind) :: dudt(nhor, nlev), dvdt(nhor, nlev), dtdt(nhor, nlev), &
         dqdt(nhor, nlev), dcwdt(nhor, nlev)
    ! v-wind after vertical diffusion
    !                                on output tendencies for wind
    ! potential temperature after vertical
    !                               diffusion. on output tendency for
    !                               actual temperature!
    ! moisture after vertical diffusion
    !                                on output tendency for specific hum.
    ! cloud water after vertical diffusion,
    !                               output tendency for cloud water
    !     ------------------------------------------------------------------
    !     input, output 1d vertical array
    !     ------------------------------------------------------------------

!    real(kind=realkind):: stddif(nlev), stdif(nlev)
    ! monin obukhovs length

    real(kind=realkind):: obkl(nhor), heatv(nhor)  
    ! virtual surface heat flux
    !     work space

    integer :: jwork(nhor)  
    ! potential temperature

    real(kind=realkind) :: wthm(nhor, nlev), wthx(nhor, nlev), wqmx(nhor, nlev), &
         wcwx(nhor, nlev), wzm(nhor, nlev), wcah(nhor, nlev), wcch( &
         nhor, nlev), wcam(nhor, nlev), wccm(nhor, nlev), wcgh(nhor, &
         nlev + 1), wcgq(nhor, nlev + 1), wkvh(nhor, nlev + 1), wkvm( &
         nhor, nlev + 1), wkvq(nhor, nlev + 1), wkvf(nhor, nlev + 1), &
         wsenf(nhor), wcflx(nhor), wdubot(nhor), wdvbot(nhor), wdtbot( &
         nhor), wdqbot(nhor), wdcwb(nhor), wtseff(nhor)
    ! temperature input + counter gradient
    ! humidity input + counter gradient
    ! cloud water( input to qvdiff )
    ! height above ground
    ! -upper diagonal for heat and spec.hum.
    ! -lower diagonal for heat and spec.hum.
    ! -upper diagonal for momentum
    ! -lower diagonal for momentum
    ! countergradient term for heat
    ! countergradient term for constitue
    ! coefficient for heat
    ! coefficient for momentum
    ! coefficient for moisture
    ! extra work space for momentum diffusn
    ! -senf
    ! -surface moisture flux(kg/m2/s)
    ! lowest layer u change from stress
    ! lowest layer v change from stress
    ! lowest layer t change from heat flux
    ! lowest layer q change from constituent
    ! lowest layer cw change(=0) due to sfc
    ! effective near surface virt. pot. temp
    call vdifnl(nhor, nlev, kstart, kstop, dtime,  t, q, cw, &
         totcov, u, v, gpot, heatv, obkl, pf, ph, dpf, dph, ts, ps, t2m, &
         q2m, taux, tauy, senf, latf, ustar, frland, &
         ustarg, pblh, &
         z0oro, dtdt, dqdt, dcwdt, dudt, dvdt, &
         jwork, wthm, wthx, wqmx, wcwx, wzm, wcah, wcch, wcam, wccm, &
         wcgh, wcgq, wkvh, wkvm, wkvq, wkvf, wsenf, wcflx, wdubot, wdvbot, &
         wdtbot, wdqbot, wdcwb, wtseff)
    return  


  end subroutine avdifnl




  subroutine vdifnl(nhor, nlev, kstart, kstop, dtime,  t, q, &
       cw, totcov, u, v, gpot, heatv, obkl, pf, ph, dpf, dph, ts, ps, &
       t2m, q2m, taux, tauy, senf, latf, ustar, frland, &!frice, 
       ustarg, &
       pblh, z0oro, dtdt, dqdt, dcwdt, dudt, dvdt,   &
       !stddif, stdif, &
       jwork, wthm, wthx, wqmx, wcwx, wzm, wcah, wcch, &
       wcam, wccm, wcgh, wcgq, wkvh, wkvm, wkvq, wkvf, wsenf, wcflx, &
       wdubot, wdvbot, wdtbot, wdqbot, wdcwb, wtseff)
    !     authors: niels woetmann nielsen and bent hansen sass,
    !     dmi december 1997
    !     modified by niels woetmann nielsen(dmi), december 1998
    !    (to be consistent with slfluxo with roughness lengths
    !     for both momentum and heat)

    implicit none  


    integer :: nhor, nlev, kstart, kstop, km1  
    !     ! input variables
    ! temperature input

    real(kind=realkind) :: t(nhor, nlev), q(nhor, nlev), cw(nhor, nlev), totcov( &
         nhor, nlev), u(nhor, nlev), v(nhor, nlev), gpot(nhor, nlev), &
         obkl(nhor), heatv(nhor), ts(nhor), ps(nhor), t2m(nhor), &
         q2m(nhor), ustar(nhor), ustarg(nhor), z0oro(nhor), frland( &
         nhor)
!    real(kind=realkind)::frice(nhor)
    real(kind=realkind)::taux(nhor), tauy(nhor), senf(nhor), &
         latf(nhor), pf(nhor, nlev), ph(nhor, nlev + 1), dpf(nhor, &
         nlev), dph(nhor, nlev + 1), pblh(nhor), dtime
    ! moisture input
    ! cloud water input
    ! cloud cover
    ! u wind input
    ! v wind input
    ! virtual potential above ground
    ! monin obukhovs length
    ! virtual surface heat flux
    ! surface temperature
    ! surface temperature
    ! temperature at 2 meters height
    ! specific humidity at 2 meters height
    ! surface friction velocity
    ! friction velocity based on heat roughn
    ! orography roughness
    ! fraction of land area
    ! fraction of ice area relative to sea f
    ! x surface stress(n)
    ! y surface stress(n)
    ! surface heat flux(w/m2)
    ! surface latent heat flux(w/m2)
    ! pressure at model levels
    ! pressure at model interface
    ! delta-p between model levels
    ! delta-p at half levels between full lv
    ! planetary boundary layer height
    ! 2 delta-t
    !     input,output 1d vertical array

!    real(kind=realkind) ::  stddif(nlev), stdif(nlev)
    !     work space
    integer :: jwork(nhor), i, k, kp  
    integer :: nlevp  

    real(kind=realkind) :: zcmlen, zm, zp, zsstab, zpblh, z, zh, zl, zzh  
    ! potential temperature
    real(kind=realkind) :: wthm(nhor, nlev), wthx(nhor, nlev), wqmx(nhor, nlev), &
         wcwx(nhor, nlev), wzm(nhor, nlev), wcah(nhor, nlev), wcch( &
         nhor, nlev), wcam(nhor, nlev), wccm(nhor, nlev), wcgh(nhor, &
         nlev + 1), wcgq(nhor, nlev + 1), wkvh(nhor, nlev + 1), wkvm( &
         nhor, nlev + 1), wkvq(nhor, nlev + 1), wkvf(nhor, nlev + 1), &
         wsenf(nhor), wcflx(nhor), wdubot(nhor), wdvbot(nhor), wdtbot( &
         nhor), wdqbot(nhor), wdcwb(nhor), wtseff(nhor)
    ! temperature input + counter gradient
    ! humidity input + counter gradient
    ! cloud water( input to qvdiff )
    ! height above ground
    ! -upper diagonal for heat and spec.hum.
    ! -lower diagonal for heat and spec.hum.
    ! -upper diagonal for momentum
    ! -lower diagonal for momentum
    ! countergradient term for heat
    ! countergradient term for constitue
    ! coefficient for heat
    ! coefficient for momentum
    ! coefficient for moisture
    ! extra work space for momentum diffusio
    ! -senf
    ! -surface moisture flux(kg/m2/s)
    ! lowest layer u change from stress
    ! lowest layer v change from stress
    ! lowest layer t change from heat flux
    ! lowest layer q change from constituent
    ! lowest layer cw change(=0) due to sfc
    ! effective near surface virt. pot. temp
    !     ! output variables
    ! u-wind after vertical diffusion

    real(kind=realkind) :: dudt(nhor, nlev), dvdt(nhor, nlev), dtdt(nhor, nlev), &
         dqdt(nhor, nlev), dcwdt(nhor, nlev)
    ! v-wind after vertical diffusion
    !     on output tendencies for wind
    ! potential temperature after vertical
    !     diffusion. on output tendency for
    !     actual temperature!
    ! moisture after vertical diffusion
    !     on output tendency for specific hum.
    ! cloud water after vertical diffusion,
    !     output tendency for cloud water
    !       declare local variables
    real(kind=realkind) :: zginv, zkappa, zp0, zlati, zlatw, zlatin, ztmp1, zdup2, &
         zdz, zrinub, zfunst, zkvn, ztmp3, ztbarp, ztbarm, ztmp2, zthvk, &
         zthvp, ztp1, zrho, zfstam, zfstah, zb, zd, zhalf, zchpbl, zarg, &
         zlamb, zlmix, zmixb, zmixh, zmixd, zfakn, zwstr, zth2, zsl1, &
         zpblk, zslask, zfmt, zfht, zwsc, zpr, z2m, z10m, zsecu
    real(kind=realkind) :: zshala, zshalb, zcl1, zcl2, zcl3, zteta, ztetap, ztetb, &
         ztetbp, ztetc, ztetcp, zdtetr, zonelim, zfrlim
    real(kind=realkind) :: zapbl, zrapbl, zexpz  
    real(kind=realkind) :: zlocal, zdh  

    real(kind=realkind) :: zqon, zqoff, zq1on, zq1off, zq2, zq3, zqqon, zqqoff  
    !     1.1 define local constants
    zginv = 1.0_realkind / gravit  
    zkappa = rair / cpair  
    zp0 = 1.e5_realkind  
    zlati = 1.0_realkind /(latvap + latice)  
    zlatw = 1.0_realkind / latvap  
    zlocal = 2.0_realkind  
    zdh = 1.0_realkind  
    zb = 5.0_realkind  
    zd = 5.0_realkind  
    zmixb = 0.3_realkind  
    zmixh = 0.35_realkind  
    zmixd = 0.35_realkind  
    zfrlim = 0.01_realkind  

    zonelim = 1.0_realkind - zfrlim  
    !     1.2 preparation for mixing length calculation
    zapbl =(1.0_realkind - zmixh) / zmixd  
    zapbl = exp( - zapbl * zapbl)  

    zrapbl = 1.0_realkind /(1.0_realkind - zapbl)  
    nlevp = nlev + 1  
    zcmlen = 30.0_realkind  
    zsecu = 0.00001_realkind  
    z2m = 2.0_realkind  

    z10m = 10.0_realkind  
    !     1.3 constants related to shallow convection
    zshala = 0.0_realkind  


    zshalb = 15.0_realkind  
    !     compute first kinematic surface fluxes( wsenf and wcflx) and
    !     'effective' grid average virtual potential temperature(wtseff)
    !     at a height of 2 meters.
    !     at first 'wtseff' is area average surface temperature,
    !     secondly it is redefined as 'effective' near surface
    !     virtual potential temperature.
    do 100 i = kstart, kstop  
       jwork(i) = 0  

       wtseff(i) = ts(i)  
       if(wtseff(i) <tmelt) then  
          zlatin = zlati  
       else  
          zlatin = zlatw  

       endif
       zrho = gpot(i, nlev) / dph(i, nlev + 1)  
       wsenf(i) = - senf(i) * zrho / cpair  

       wcflx(i) = - latf(i) * zrho * zlatin  
       wzm(i, nlev) = gpot(i, nlev) * zginv  


       wthm(i, nlev) = t(i, nlev) *(zp0 / pf(i, nlev) ) **zkappa  
       !     redefine 'wtseff'
       ztmp1 = wzm(i, nlev) - z2m  
       if(ztmp1>=0.0_realkind) ztmp1 = max(ztmp1, zsecu)  

       if(ztmp1<0.0_realkind) ztmp1 = min(ztmp1, - zsecu)  
       ztmp3 = t2m(i) *(1._realkind - cw(i, nlev) /(zlatin * cpair * t2m(  i) ) )
       ztmp2 = ztmp3 *((zp0 / ps(i) ) **zkappa) *(1.0_realkind + 0.607717_realkind * &
            q2m(i) - cw(i, nlev) )

       wtseff(i) = ztmp2 +(wthm(i, nlev) - ztmp2) *(z10m - z2m) / ztmp1


100 enddo
    !     compute relative height and the
    !     potential temperature.
    do k = nlev - 1, ntopfl, - 1  
       do i = kstart, kstop  
          wzm(i, k) = gpot(i, k) * zginv  
          wthm(i, k) = t(i, k) *(zp0 / pf(i, k) ) **zkappa  
       enddo

    enddo
    !     convert the surface fluxes to lowest level tendencies
    do 90 k = 1, nlev  
       do 85 i = kstart, kstop  
          dudt(i, k) = 0.0_realkind  
          dvdt(i, k) = 0.0_realkind  
          dtdt(i, k) = 0.0_realkind  
          dqdt(i, k) = 0.0_realkind  
          dcwdt(i, k) = 0.0_realkind  
85     enddo

90  enddo
    !     4. compute local exchange coefficients
    !     4.1 initialize exchange coefficients and counter gradient
    !     terms
    do 120 k = 1, nlevp  
       do 110 i = kstart, kstop  
          wkvm(i, k) = 0.0_realkind  
          wkvh(i, k) = 0.0_realkind  
          wkvq(i, k) = 0.0_realkind  
          wkvf(i, k) = 0.0_realkind  
          wcgh(i, k) = 0.0_realkind  
          wcgq(i, k) = 0.0_realkind  
110    enddo
120 enddo
    do 140 k = nlev - 1, ntopfl, - 1  

       do 130 i = kstart, kstop  
          !     4.2 vertical shear squared, min value of
          !    (delta v)**2 prevents zero shear.
          zdup2 =(u(i, k) - u(i, k + 1) ) **2 +(v(i, k) - v(i, &
               k + 1) ) **2
          zdup2 = max(zdup2, 1.e-14_realkind)  
          zdz = wzm(i, k) - wzm(i, k + 1)  

          zdup2 = zdup2 / zdz**2  

          zhalf = 0.5_realkind *(wzm(i, k) + wzm(i, k + 1) )  
          !     4.3  manipulation with pblh(i), and compute
          !     mixing length paramters.
          zpblh = max(pblh(i), zlocal * zcmlen / zmixb)  
          pblh(i) = zpblh  
          zchpbl = zmixh * zpblh  
          zarg =(max(zhalf, zchpbl) - zchpbl) /(zmixd * zpblh)  
          zexpz = exp( - zarg * zarg)  
          zexpz = max(zexpz, zapbl)  
          zchpbl = zmixb * zpblh  
          zlamb = zrapbl *(1.0_realkind - zexpz) *(zcmlen - zapbl * zchpbl) &
               + zchpbl * zexpz
          zlmix = carman * zhalf /(1.0_realkind + carman * zhalf / zlamb)  

          zkvn = zlmix * zlmix * sqrt(zdup2)  
          !     4.4 static stability(use virtual potential
          !     temperature). special treatment of exchange
          !     coefficients in cloudy layers.
          zcl1 = min(totcov(i, k), totcov(i, k + 1) )  
          zcl1 = min(zcl1, zonelim)  
          zcl2 = max(totcov(i, k + 1) - totcov(i, k), 0.0_realkind)  

          zcl3 = max(totcov(i, k) - totcov(i, k + 1), 0.0_realkind)  
          if(t(i, k) <tmelt) then  
             zlatin = zlati  
          else  
             zlatin = zlatw  

          endif
          zteta = wthm(i, k)  

          zthvk = zteta *(1.0_realkind + 0.61_realkind * q(i, k) - cw(i, k) )  
          if(totcov(i, k) >zfrlim) then  
             ztetb = wthm(i, k)  
             ztetc = wthm(i, k) *(1.0_realkind + zshala * cw(i, k) /(zlatin &
                  * cpair * t(i, k) * totcov(i, k) ) )
          else  
             ztetb = wthm(i, k)  
             ztetc = wthm(i, k)  

          endif
          if(t(i, k + 1) <tmelt) then  
             zlatin = zlati  
          else  
             zlatin = zlatw  

          endif
          ztetap = wthm(i, k + 1)  

          zthvp = ztetap *(1.0_realkind + 0.61_realkind * q(i, k + 1) - cw(i, k + 1) &
               )
          if(totcov(i, k + 1) >zfrlim) then  
             ztetbp = wthm(i, k + 1)  
             ztetcp = wthm(i, k + 1) *(1.0_realkind + zshala * cw(i, k + 1) &
                  /(zlatin * cpair * t(i, k + 1) * totcov(i, k + 1) ) )
          else  
             ztetbp = wthm(i, k + 1)  
             ztetcp = wthm(i, k + 1)  

          endif

          zdtetr =((ztetc - ztetcp) *(zcl2 - zcl3) +(zteta - &
               ztetap) *(1.0_realkind - zcl1 - zcl2 - zcl3) ) /(1.0_realkind - zcl1)

          zsstab = 2.0_realkind * gravit * zdtetr /(zdz *(zthvk + zthvp) )  
          !     4.5 richardson number, stable and unstable
          !     modifying functions

          zrinub = zsstab / zdup2  
          if(zrinub>=0.0_realkind) then  
             zfstam = 2.0_realkind * zb * zrinub / sqrt(1.0_realkind + zd * zrinub)  
             zfstam = 1.0_realkind /(1.0_realkind + zfstam)  
             zfstah = 2.0_realkind * zb * zrinub * sqrt(1.0_realkind + zdh * zrinub)  
             zfstah = 1.0_realkind /(1.0_realkind + zfstah)  

             zfunst = 1.0_realkind  
             wkvm(i, k + 1) = zkvn * zfstam *(1.0_realkind - zcl1) + zshalb * &
                  zcl1
             wkvh(i, k + 1) = zkvn * zfstah**(1.0_realkind - zcl1) + zshalb * &
                  zcl1
             wkvq(i, k + 1) = wkvh(i, k + 1)  

          else  
             zfunst = max(1.0_realkind - 18.0_realkind * zrinub, 0.0_realkind)  
             zfstam = sqrt(zfunst)  
             zfstah = zfstam  
             wkvm(i, k + 1) = zkvn * zfstam *(1.0_realkind - zcl1) + zshalb * &
                  zcl1
             wkvh(i, k + 1) = wkvm(i, k + 1)  
             wkvq(i, k + 1) = wkvm(i, k + 1)  

          endif

130    enddo



140 enddo
    !     end of local exchange coefficient calculations
    !     5. determine non local effects including counter gradient
    !     terms.
    !     add the counter gradient terms to potential temperature
    !     and specific in the boundary layer humidity.
    do 770 k = 1, npbl - 1  
       kp = k + 1  
       do 771 i = kstart, kstop  
          zm = wzm(i, nlevp - k)  
          zp = wzm(i, nlevp - kp)  
          if(pblh(i) >zm) then  
             zqon = 1.0_realkind  
             zqoff = 0.0_realkind  
          else  
             zqon = 0.0_realkind  
             zqoff = 1.0_realkind  

          endif
          z = 0.5_realkind *(zm + zp)  

          zh = zqon * z / pblh(i)  
          !     zl must not be zero
          zl = zqon * z / obkl(i) + zqoff  
          if(zh>1.0_realkind) then  
             zq3 = 0.0_realkind  
          else  
             zq3 = 1.0_realkind  

          endif

          zzh = zqon * zq3 *(1.0_realkind - zh) **2  
          if(heatv(i) >0.0_realkind) then  
             zq1on = 1.0_realkind 
             zq1off = 0.0_realkind  
          else  
             zq1on = 0.0_realkind  
             zq1off = 1.0_realkind  

          endif

          zq2 = zqon *(1.0_realkind - zq1on)  
          !     5.1 stable and neutral:
          !     prevent zpblk to become too small in very stable conditions
          if(zl>1.0_realkind) then  
             zqqon = 1.0_realkind  
             zqqoff = 0.0_realkind  
          else  
             zqqon = 0.0_realkind  
             zqqoff = 1.0_realkind  

          endif
          ztmp1 = max(z0oro(i), chrmax)  
          zslask = log((z + ztmp1) / ztmp1)  
          zpblk = log((z + chrmax) / chrmax)  
          zslask = zslask / zpblk  
          zl = zq1on * zl + zq1off *(frland(i) * zslask * zl + &
              (1.0_realkind - frland(i) ) * zl)
          zslask = ustar(i) * pblh(i) * vk * zh * zzh  
          zpblk = zslask /(zq1on + zq1off *(zqqon *(betas + zl) &
               + zqqoff *(1.0_realkind + betas * zl) ) )

          zpblk = zq2 * zpblk  
          !     5.2  correction of km/h/q due to orography roughness.
          !     the correction does not affect the local values
          !     above the pbl.

          ztmp1 = ustarg(i) / ustar(i)  

          wkvm(i, nlevp - k) = zq2 *(zq3 * max(zpblk, ckmin) &
               +(1.0_realkind - zq3) * wkvm(i, nlevp - k) ) +(1.0_realkind - zq2) * wkvm( &
               i, nlevp - k)
          !     5.3   zq3=0 if zm< pblh < 0.5(zm+zp); else zq3=1;
          !     and if zq3=0  km/h/q must be assigned their
          !     local values at z=0.5(zm+zp).
          wkvh(i, nlevp - k) = zq2 * zq3 * ztmp1 * wkvm(i, nlevp - &
               k) +(1.0_realkind - zq2 * zq3) * wkvh(i, nlevp - k)

          wkvq(i, nlevp - k) = wkvh(i, nlevp - k)  

          zq3 = zq3 * zqon * zq1on  
          !     5.4  unstable case(surface layer values)
          zslask = 1.0_realkind - betam * zl  

          zslask = zq3 * zslask +(1.0_realkind - zq3)  
          zfmt = zslask**onet  
          zwsc = ustar(i) * zfmt  
          zpblk =(1.0_realkind - zq3) * zpblk + zq3 *(vk * zwsc * pblh(i) &
               * zh * zzh)
          zsl1 = 1.0_realkind - betah * zl  

          zsl1 = zq3 * zsl1 +(1.0_realkind - zq3)  
          zfht = sqrt(zsl1)  
          zwstr =(abs(heatv(i) ) * g * pblh(i) / wtseff(i) ) ** &
               onet
          zslask = 2.0_realkind * onet  
          zslask = zh**zslask  
          zslask = 33.1_realkind * zslask  
          zfakn = zslask * zwstr / zwsc  
          zth2 = zfakn /(pblh(i) * zwsc)  


          wcgh(i, nlevp - k) = zq3 * wsenf(i) * zth2  
          !     5.5   modification of the prandtl number(zpr) to take into
          !     account the different surface friction velocities for
          !     momentum(ustar) and heat/moisture(ustarg).

          zpr = zq3 *(zfmt /(zfht * ztmp1) + zslask * vk * zh * &
               zwstr / zwsc)


          wkvm(i, nlevp - k) = zq3 * max(zpblk, ckmin) +(1.0_realkind - zq3) &
               * wkvm(i, nlevp - k)
          !     5.6   final adjustment of prandtl number 'zpr'
          !     and exchange coefficients for heat and moisture.
          zpr = zq3 * zpr +(1.0_realkind - zq3)  
          zpr = max(zpr, zfrlim)  
          wkvh(i, nlevp - k) = zq3 * max(zpblk / zpr, ckmin) &
               +(1.0_realkind - zq3) * wkvh(i, nlevp - k)

          wkvq(i, nlevp - k) = wkvh(i, nlevp - k)  
771    enddo


770 enddo
    !     end of non local calculations.
    !     6. preparations for solving the diffusion equations.
    do i = kstart, kstop  
       if(ts(i) <tmelt) then  
          zlatin = zlati  
       else  
          zlatin = zlatw  
       endif
       wsenf(i) = - senf(i)  
       wcflx(i) = - latf(i) * zlatin  
       ztmp1 = dtime * gravit / dpf(i, nlev)  
       wdubot(i) = - taux(i) * ztmp1  
       wdvbot(i) = - tauy(i) * ztmp1  
       wdqbot(i) = wcflx(i) * ztmp1  
       wdcwb(i) = 0.0_realkind  
       wdtbot(i) = wsenf(i) * ztmp1 / cpair  


    enddo
    !    (above the pbl):
    do 160 k = 1, nlev - npbl  
       do 150 i = kstart, kstop  
          wthx(i, k) = wthm(i, k)  
          wqmx(i, k) = q(i, k)  
          wcwx(i, k) = cw(i, k)  
150    enddo


160 enddo
    !    (in the pbl):
    do 180 k = nlev - npbl + 1, nlev  
       do 170 i = kstart, kstop  
          if(k<nlev) then  
             ztbarp = 0.5_realkind *(t(i, k + 1) + t(i, k) )  
          else  
             ztbarp = t(i, nlev)  
          endif
          ztbarm = 0.5_realkind *(t(i, k) + t(i, k - 1) )  
          ztmp1 = dtime * gravit /(rair * dpf(i, k) )  
          wthx(i, k) = wthm(i, k) + ztmp1 *(ph(i, k + 1) * wkvh( &
               i, k + 1) * wcgh(i, k + 1) / ztbarp - ph(i, k) * wkvh(i, &
               k) * wcgh(i, k) / ztbarm)
          wqmx(i, k) = q(i, k) + ztmp1 *(ph(i, k + 1) * wkvh(i, &
               k + 1) * wcgq(i, k + 1) / ztbarp - ph(i, k) * wkvh(i, k) &
               * wcgq(i, k) / ztbarm)
          wcwx(i, k) = cw(i, k)  
170    enddo
180 enddo
    !     &
    !     check for negatives q's and put the original profiles back
    !     -
    do 183 k = nlev - npbl + 1, nlev  
       do 182 i = kstart, kstop  
          if(wqmx(i, k) <0.0_realkind) jwork(i) = 1  
182    enddo

183 enddo
    do 185 k = nlev - npbl + 1, nlev  
       do 184 i = kstart, kstop  
          if(abs(real(jwork(i),realkind)-0.1_realkind)<1.e-14_realkind)then
             wqmx(i,k) = q(i,k)  
          endif
184    enddo

185 enddo
    !     6.1 determine superdiagonal(ca(k)) and subdiagonal(cc(k))
    !     coefficients of the tridiagonal diffusion matrix. the dia-
    !     gonal elements are a combination of ca and cc; they are not
    !     explicitly provided to the solv.
    do 210 k = ntopfl, nlev - 1  
       do 200 i = kstart, kstop  
          ztbarp = 0.5_realkind *(t(i, k + 1) + t(i, k) )  
          ztmp2 = dtime / dph(i, k + 1) *(gravit * ph(i, k + 1) &
               /(ztbarp * rair) ) **2
          wcch(i, k + 1) = wkvh(i, k + 1) * ztmp2 / dpf(i, k + 1)  
          wcah(i, k) = wkvh(i, k + 1) * ztmp2 / dpf(i, k)  
          wccm(i, k + 1) = wkvm(i, k + 1) * ztmp2 / dpf(i, k + 1)  
          wcam(i, k) = wkvm(i, k + 1) * ztmp2 / dpf(i, k)  
200    enddo
210 enddo
    !     &
    !     the last element of the upper diagonal is zero.
    !     -
    do 220 i = kstart, kstop  
       wcah(i, nlev) = 0.0_realkind  
       wcam(i, nlev) = 0.0_realkind  



220 enddo
    !     7. solution fo diffusion equations
    !     7.1 solution for momentum

    call mvdiff(nhor, nlev, kstart, kstop, u, v, wdubot, wdvbot, &
         wcam, wccm, dudt, dvdt, wkvh, wkvm, wkvf)
    !     7.2 solution for specific cloud water

    call qvdiff(nhor, nlev, kstart, kstop, wcwx, wdcwb, wcah, wcch, &
         dcwdt, wkvh, wkvm)
    !     7.3 negative cloud water is adjusted to zero,
    !     and specific humidity is adjusted correspondingly
    !     to conserve total specific humidity.
    do 224 k = 1, nlev  
       do 222 i = kstart, kstop  
          if(dcwdt(i, k) <0.0_realkind) then  
             wqmx(i, k) = q(i, k) - dcwdt(i, k)  
             dcwdt(i, k) = 0.00000_realkind  
          endif
222    enddo

224 enddo
    do 228 k = 2, nlev  
       km1 = k - 1  
       do 226 i = kstart, kstop  
          if(wqmx(i, km1) <0.0_realkind) then  
             wqmx(i, k) = wqmx(i, k) + wqmx(i, km1) * dpf(i, km1) &
                  / dpf(i, k)
             wqmx(i, km1) = 0.0_realkind  
          endif
226    enddo

228 enddo
    !     7.4 solution for specific humidity

    call qvdiff(nhor, nlev, kstart, kstop, wqmx, wdqbot, wcah, wcch, &
         dqdt, wkvh, wkvm)
    !     7.5 solution for potential temperature.

    call qvdiff(nhor, nlev, kstart, kstop, wthx, wdtbot, wcah, wcch, &
         dtdt, wkvh, wkvm)
    !     8. convert the diffused fields back to diffusion
    !     tendencies.
    do 400 k = ntopfl, nlev  
       do 410 i = kstart, kstop  
          dudt(i, k) =(dudt(i, k) - u(i, k) ) / dtime  
          dvdt(i, k) =(dvdt(i, k) - v(i, k) ) / dtime  
          dqdt(i, k) =(dqdt(i, k) - q(i, k) ) / dtime  


          dcwdt(i, k) =(dcwdt(i, k) - cw(i, k) ) / dtime  
          !     9. compute actual temperature from potential
          !     temperature.
          ztp1 = dtdt(i, k) *(pf(i, k) / zp0) **zkappa  
          dtdt(i, k) =(ztp1 - t(i, k) ) / dtime  
410    enddo

400 enddo
    return  
  end subroutine vdifnl

  subroutine qvdiff(nhor, nlev, kstart, kstop, qm1, qflx, ca, cc, &
       qp1, we, wfq)
    !     ! solve vertical diffusion equation for moisture with explicit
    !     ! surface flux.
    !     ! procedure for solution of the implicit equation follows
    !     ! richtmyer and morton(1967,pp198-199).

    implicit none  

    integer :: nhor, nlev, kstart, kstop, jk, jl  
    !     ! input variables
    ! initial moisture

    real(kind=realkind) :: qm1(nhor, nlev), qflx(nhor), ca(nhor, nlev), cc(nhor, &
         nlev)
    ! sfc q flux into lowest model level
    ! -upper diagonal coeff. of tri-diagonal
    ! -lower diagonal coeff. of tri-diagonal
    !     ! output variables
    ! final moisture

    real(kind=realkind) :: qp1(nhor, nlev)  
    !     ! local workspace
    ! terms appearing in soln of tridiagonal



    real(kind=realkind) :: we(nhor, nlev), wfq(nhor, nlev), ztmp  
    ! terms appearing in soln of tridiagonal
    ! temporary workspace
    !     ! calculate e(k) and fq(k).
    !     ! these are terms required in solution of tridiagonal matrix
    !     ! defined by implicit diffusion eqn.
    do 10 jl = kstart, kstop  
       ztmp = 1.0_realkind + ca(jl, ntopfl)  
       we(jl, ntopfl) = ca(jl, ntopfl) / ztmp  
       wfq(jl, ntopfl) = qm1(jl, ntopfl) / ztmp  

10  enddo

    do 30 jk = ntopfl + 1, nlev - 1  
       do 20 jl = kstart, kstop  
          ztmp = 1.0_realkind + ca(jl, jk) + cc(jl, jk) - cc(jl, jk) * we( &
               jl, jk - 1)
          we(jl, jk) = ca(jl, jk) / ztmp  
          wfq(jl, jk) =(qm1(jl, jk) + cc(jl, jk) * wfq(jl, jk - &
               1) ) / ztmp
20     enddo


30  enddo
    !     !  bottom level:(includes  surface fluxes
    do 40 jl = kstart, kstop  
       ztmp = 1.0_realkind + cc(jl, nlev) - cc(jl, nlev) * we(jl, nlev - 1)  
       we(jl, nlev) = 0.0_realkind  
       wfq(jl, nlev) =(qm1(jl, nlev) + qflx(jl) + cc(jl, nlev) &
            * wfq(jl, nlev - 1) ) / ztmp

40  enddo
    !     ! perform back substitution
    do 50 jl = kstart, kstop  
       qp1(jl, nlev) = wfq(jl, nlev)  

50  enddo
    do 70 jk = nlev - 1, ntopfl, - 1  
       do 60 jl = kstart, kstop  
          qp1(jl, jk) = wfq(jl, jk) + we(jl, jk) * qp1(jl, jk + 1)  
60     enddo
70  enddo
    return  



  end subroutine qvdiff


  subroutine mvdiff(nhor, nlev, kstart, kstop, um1, vm1, uflx, &
       vflx, ca, cc, up1, vp1, we, wfu, wfv)
    !     ! vertical momentum diffusion with explicit surface flux.
    !     ! solve the vertical diffusion equation for momentum.
    !     ! procedure for solution of the implicit equation follows
    !     ! richtmyer and morton(1967,pp.198-199)
    implicit none  


    integer :: nhor, nlev, kstart, kstop, jk, jl  
    !     ! input variables
    ! u wind input


    real(kind=realkind) :: um1(nhor, nlev), vm1(nhor, nlev), uflx(nhor), vflx( &
         nhor), ca(nhor, nlev), cc(nhor, nlev)
    ! v wind input
    ! sfc u flux into lowest model level
    ! sfc v flux into lowest model level
    ! -upper diagonal coeff. of tri-diagonal
    ! -lower diagonal coeff. of tri-diagonal
    !     ! output variables
    ! u wind after mvdiff


    real(kind=realkind) :: up1(nhor, nlev), vp1(nhor, nlev)  
    ! v wind after mvdiff
    !     ! local workspace
    ! terms appearing in soln of tridiagonal



    real(kind=realkind) :: we(nhor, nlev), wfu(nhor, nlev), wfv(nhor, nlev), &
         ztmp
    ! terms appearing in soln of tridiagonal
    ! terms appearing in soln of tridiagonal
    ! temporary workspace
    !     ! calculate e(k),fu(k) and fv(k).
    !     ! these are terms required in solution of tridiagonal matrix defin
    !     ! by implicit diffusion eqn.
    do 10 jl = kstart, kstop  
       ztmp = 1.0_realkind + ca(jl, ntopfl)  
       we(jl, ntopfl) = ca(jl, ntopfl) / ztmp  
       wfu(jl, ntopfl) = um1(jl, ntopfl) / ztmp  
       wfv(jl, ntopfl) = vm1(jl, ntopfl) / ztmp  

10  enddo
    do 30 jk = ntopfl + 1, nlev - 1  
       do 20 jl = kstart, kstop  
          ztmp = 1.0_realkind + ca(jl, jk) + cc(jl, jk) - cc(jl, jk) * we( &
               jl, jk - 1)
          we(jl, jk) = ca(jl, jk) / ztmp  
          wfu(jl, jk) =(um1(jl, jk) + cc(jl, jk) * wfu(jl, jk - &
               1) ) / ztmp
          wfv(jl, jk) =(vm1(jl, jk) + cc(jl, jk) * wfv(jl, jk - &
               1) ) / ztmp
20     enddo


30  enddo
    !     !  bottom level:(includes  surface fluxes
    do 40 jl = kstart, kstop  
       ztmp = 1.0_realkind + cc(jl, nlev) - cc(jl, nlev) * we(jl, nlev - 1)  
       we(jl, nlev) = 0.0_realkind  
       wfu(jl, nlev) =(um1(jl, nlev) + uflx(jl) + cc(jl, nlev) &
            * wfu(jl, nlev - 1) ) / ztmp
       wfv(jl, nlev) =(vm1(jl, nlev) + vflx(jl) + cc(jl, nlev) &
            * wfv(jl, nlev - 1) ) / ztmp


40  enddo
    !     ! perform back substitution
    do 50 jl = kstart, kstop  
       up1(jl, nlev) = wfu(jl, nlev)  
       vp1(jl, nlev) = wfv(jl, nlev)  

50  enddo
    do 70 jk = nlev - 1, ntopfl, - 1  
       do 60 jl = kstart, kstop  
          up1(jl, jk) = wfu(jl, jk) + we(jl, jk) * up1(jl, jk + 1)  
          vp1(jl, jk) = wfv(jl, jk) + we(jl, jk) * vp1(jl, jk + 1)  
60     enddo
70  enddo
    return  
  end subroutine mvdiff

end module comdfb
