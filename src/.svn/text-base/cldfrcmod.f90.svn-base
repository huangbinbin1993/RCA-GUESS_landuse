module cldfrcmod
  use decomp,only:stop_program,realkind
  use confys
  use escom

  implicit none
  private

  public acldfrc,axucld
contains

  subroutine acldfrc (nhor, nlev, kstart, kstop, pf, rdph, t, q, &
       cwat, fice, omf, kf_top, kf_base, totcov, oldcov, dpf, mft, oro, &
       snowh, frland, frice, concld, ts, ps, detr, cov2d, shal_cgj, &
       kstep, rad_cloud, zdx, zdy, orogsigm, dzh, rhcrit, lesat, qvap, &
       qliq, qice, gpot,pblh, shalcld, relhum)

    implicit none  
    logical :: lesat  
    integer :: nhor, nlev, kstart, kstop  
    real(kind=realkind) :: snowh (nhor), ts (nhor), ps (nhor)  
    real(kind=realkind) :: frland (nhor), frice (nhor)  
    integer :: oro (nhor)  
    real(kind=realkind) :: t (nhor, nlev), q (nhor, nlev), omf (nhor, nlev)
    real(kind=realkind)::pf(nhor, nlev), dpf (nhor, nlev), rdph (nhor, nlev)
    real(kind=realkind)::mft (nhor, nlev)  , detr (nhor, nlev)
    real(kind=realkind) :: cwat (nhor, nlev)  
    real(kind=realkind) :: dzh (nhor, nlev)  
    real(kind=realkind) :: rhcrit (nhor, nlev)  
    real(kind=realkind) :: cov2d (nhor)  
    real(kind=realkind) :: concld (nhor, nlev), totcov (nhor, nlev)  
    integer :: jl, jk, k700 (nhor), kf_top (nhor), kf_base (nhor)  
    integer :: shal_cgj (nhor)  
    integer :: kstep  
    real(kind=realkind) :: rad_cloud (nhor, nlev)  
    real(kind=realkind) :: shalcld (nhor, nlev)  
    real(kind=realkind) :: zdx (nhor), zdy (nhor)  
    real(kind=realkind) :: orogsigm (nhor)  
    real(kind=realkind) :: qvap (nhor, nlev), qliq (nhor, nlev)  
    real(kind=realkind) :: qice (nhor, nlev)  
    real(kind=realkind) :: relhum (nhor, nlev)  
    real(kind=realkind) :: fice (nhor, nlev)  
    real(kind=realkind) :: oldcov (nhor, nlev)  
!cgjpblh
    real(kind=realkind) :: gpot (nhor, nlev)  
    real(kind=realkind) :: pblh (nhor)  

    !  do some necessary computations
    !  calc rpdeli for cldfrc
    ! determining the level closest to 700mb
    do jl = kstart, kstop  
       k700(jl) = nlev  
       do jk = 1, nlev - 1  
          if(pf(jl, jk)<7.e4_realkind)k700(jl)=jk
       enddo
    enddo
    !
    call cldfrc (kstart, kstop, nhor, nlev, pf, rdph, t, q, cwat, &
         fice, omf, kf_base, kf_top, totcov, oldcov, k700, dpf, mft, oro, &
         snowh, frland, frice, detr, concld, ts, ps, cov2d, shal_cgj, &
         kstep, rad_cloud, zdx, zdy, orogsigm, dzh, rhcrit, lesat, qvap, &
         qliq, qice, gpot, pblh, shalcld, relhum)
    return  
  end subroutine acldfrc

  subroutine cldfrc(kstart  ,kstop   ,plond   ,plev    , &
       pmid  ,rpdeli ,temp  ,q     ,cwat  ,fice,      &
       omga    ,                                      &
       kf_base  ,kf_top  ,cloud,oldcov,k700,pdel,     &
       cmfmc   ,oro     ,                             &
       snowh   ,frland  ,frice  ,                     &
       zdu     ,concld  ,ts      ,ps      ,cov2d,     &
!cgj       zdu     ,concld  ,ts      ,ps      ,ptropo,     &
       shal_cgj,kstep,rad_cloud,                      &
       zdx,zdy,orogsigm,dzh,rhlim,lesat,              &
       qvap,qliq,qice,gpot,pblh,                      &
       shalcld,rh)
    !     
    !     compute cloud fraction using scheme of j.m.slingo, 
    !     as modified by j.j.hack and j.t.kiehl
    !     
    !     this scheme is based on the operational scheme used in the ecmwf model
    !     a full description of its development can be found in slingo (1987),
    !     which appears in the qjrms july issue.  a number of modifications have
    !     been introduced to the original scheme in the following implementation 
    !     
    !---------------------------code history--------------------------------
    !     
    !     original version: based on code from j. slingo 
    !     modified:         j. j. hack, may 1990
    !     modified:         j. t. kiehl, june, 1990
    !     modified:         j. j. hack, january 1991
    !     rewritten:        b. p. briegleb, november, 1991
    !     standardized:     j. rosinski, june 1992
    !     last updated:     j. j. hack, july 1992 
    !     reviewed:         j. kiehl, april 1996
    !     modified          p. rasch, mar 1997



    implicit none
    !------------------------------arguments--------------------------------
    !     
    !     input 0-d
    !     
    integer ::kstart,kstop,plond,plev
    !     
    !     input 1-d
    !     
    integer:: k700(plond)       ! level nearest to 700 mb over ocean
    integer:: oro(plond)        ! land/ocean/seaice flag
    real(kind=realkind) ::snowh(plond)         ! snow depth (liquid water equivalent)
    real(kind=realkind):: ts(plond)            ! surface temperature
    real(kind=realkind):: ps(plond)            ! surface pressure
    !     
    !     input 2-d
    !     
    real(kind=realkind):: pmid(plond,plev)     ! midpoint pressures
    real(kind=realkind):: rpdeli(plond,plev)   ! 1./(pmid(k+1)-pmid(k))
    real(kind=realkind):: temp(plond,plev)     ! temperature
    real(kind=realkind):: q(plond,plev)        ! specific humidity
    real(kind=realkind):: omga(plond,plev)     ! vertical pressure velocity
    real(kind=realkind):: cmfmc(plond,plev)    ! convective mass flux--m sub c
    real(kind=realkind):: pdel(plond,plev)     ! pressure depth of layer
    real(kind=realkind):: zdu(plond,plev)      ! detrainment rate from deep convection
    !     vo +++
    !     
    !     output 1-d
    !     
    real(kind=realkind):: cov2d(plond)         ! 2-d cloud fraction
    real(kind=realkind):: shal2d(plond)        ! 2-d cloud fraction
    real(kind=realkind):: cu2d(plond)          ! 2-d cloud fraction
    !     gj240302
    real(kind=realkind):: cov2d_odp(plond)     ! 2-d cloud fraction
    !     gj240302
    !     vo ---
    !     
    !     output 2-d
    !     
    real(kind=realkind):: cloud(plond,plev)    ! cloud fraction
    real(kind=realkind):: concld(plond,plev)   ! column convective cloud amount
    !     
    !------------------------------work space--------------------------------
    !     
    !     here follows 0-d integer work numbers
    integer:: i,k               ! longitude, level indices
    integer:: kp1
    !     here follows 0-d real work numbers
    real(kind=realkind):: cappa                ! r/cp
    real(kind=realkind):: cld                  ! intermediate scratch variable (low cld)
    real(kind=realkind):: dthdp                ! lapse rate (intermediate variable)
    real(kind=realkind):: premib               ! bottom pressure bound of middle cloud
    real(kind=realkind) ::pretop               ! pressure bounding high cloud
    !     vo      real rhb                   ! intermediate scratch variable
    real(kind=realkind):: rhdif                ! intermediate scratch variable
    real(kind=realkind):: strat                ! intermediate scratch variable
    real(kind=realkind):: bvf(plond,plev)      ! brunt-vaisalla frequency
    real(kind=realkind):: rbvflim              ! bound on inverse of bvf
    real(kind=realkind):: rho                  ! local density (used to calculate bvf)
    real(kind=realkind):: rhlim(plond,plev)    ! local rel. humidity threshold estimate 
    real(kind=realkind):: rhden                ! intermediate scratch variable
    real(kind=realkind):: rhdif2               ! intermediate scratch variable
    !     vo      real pdepth                ! intermediate scratch variable
    !     vo      real stratfac              ! intermediate scratch variable
    real(kind=realkind):: rhminl               ! minimum rh for low stable clouds
    real(kind=realkind):: rhminh               ! minimum rh for high stable clouds
    real(kind=realkind):: coef1                ! coefficient to convert mass flux to mb/d
    real(kind=realkind):: pnot                 ! reference pressure
    !     here follows 1-d logical work arrays
    logical:: lol(plond)        ! region of low level cloud
    logical:: cldbnd(plond)     ! region below high cloud boundary
    !     here follows 1-d integer work arrays
    integer:: kdthdp(plond)
    !     here follows 1-d real work arrays
    real(kind=realkind):: cld8(plond)          ! low cloud fraction estimate
    real(kind=realkind):: cld9(plond)          ! mid and high cloud fraction estimate
    real(kind=realkind):: cck(plond)           ! convective cloud per level (assuming
    !     random overlap in convective layer)
    real(kind=realkind):: dtdpmn(plond)        ! most stable lapse rate below 750 mb
    real(kind=realkind):: mcbar(plond)         ! mean convective scale motion in column
    real(kind=realkind):: dpsum(plond)         ! vertical sum of delta-p (k-1 levels)
    real(kind=realkind):: clrsky(plond)        ! temporary used in random overlap calc
    real(kind=realkind):: thetas(plond)        ! potential temperature surface
    real(kind=realkind):: clc(plond)           ! coloumn convective cloud amount
    !     here follows 2-d real work arrays
    real(kind=realkind):: cldst(plond,plev)    ! cloud fraction
    real(kind=realkind):: dthtdp(plond,plev)   ! lapse rate (d theta/dp) below 750 mb
    real(kind=realkind):: qs(plond,plev)       ! saturation specific humidity
    real(kind=realkind):: rh(plond,plev)       ! relative humidity
    real(kind=realkind):: theta(plond,plev)    ! potential temperature
    real(kind=realkind):: thetav(plond,plev)   ! virtual potential temperature
    real(kind=realkind):: ccldt(plond)
    real(kind=realkind):: zrth,ztv
    integer:: kf_top(plond)
    integer:: kf_base(plond)
    real(kind=realkind):: ess,pss,mcmbh,wrk1,wrk2,rhmax
    real(kind=realkind):: cu0(plond,plev),cu1(plond,plev),cu2(plond,plev)
    real(kind=realkind):: xucldarray(plond,plev)
    !     gjshallow conv
    integer:: shal_cgj(plond)
    integer:: ktropo_bvf(plond),ktropo_wmo(plond),ktt,kk
    integer:: kstep
    real(kind=realkind):: rad_cloud(plond,plev),pdep(plond)
    real(kind=realkind):: cwat(plond,plev)
    real(kind=realkind):: rhcrit(plond,plev)
    real(kind=realkind):: zdx(plond),zdy(plond)
    real(kind=realkind):: orogsigm(plond)
    !     gjkfarm
    real(kind=realkind):: shalcld(plond,plev)
    real(kind=realkind):: shalcld2(plond,plev)
    real(kind=realkind):: qvap(plond,plev)
    real(kind=realkind):: qliq(plond,plev),sr
    real(kind=realkind):: qice(plond,plev)
    real(kind=realkind):: ptropo(plond),stropo(plond)
    real(kind=realkind):: ptropo_out(plond),ttropo(plond),ztropo(plond)
    real(kind=realkind):: zztr(plond),zzcnt,zztint,zdz,zzt,zztav,zwmo
    real(kind=realkind):: lscld(plond,plev),wzf(plond,plev)
    !     gjkfarm
    real(kind=realkind):: cldopd(plond,plev)
    real(kind=realkind):: dzh(plond,plev)
    real(kind=realkind):: cwliq,cwice,lwp,iwp
    !     gj
    real(kind=realkind):: pdepst(plond)
    real(kind=realkind):: fice(plond,plev)
    real(kind=realkind):: oldcov(plond,plev)
    real(kind=realkind):: zfrice,frice(plond),frland(plond)
!cgjnewrh
    real(kind=realkind):: zzm,zginv,zsigma,zwrk1,zwrk2,zwrk3
    real(kind=realkind):: zbvf,zdetr,zcdet 
    real(kind=realkind):: zbvfg(plond),zbvflim(plond)
    real(kind=realkind):: gpot(plond,plev)     ! midpoint geopotential height
    real(kind=realkind):: rhmin(plond,plev)    ! rhminh RH threshold can now vary 
    real(kind=realkind):: pblh(plond)
    logical:: lesat
    !--------------------------------------------------------------------

    !     
    !     set bound for inverse of brunt-vaisalla frequency and minimum relative
    !     humidity thresholds for stable clouds.  these are the principal 
    !     "disposable" parameters for the cloud fraction scheme
    !     
    rbvflim = 1.0_realkind/0.00035_realkind
!cgj020711
    rhminl = 0.7_realkind
    rhminh = 0.8_realkind
!cgj020711
    rhmax  = 0.98_realkind
    pnot = 1.e5_realkind
    mcmbh=0.0_realkind
    zginv=1.0_realkind/gravit
    zwmo=-2._realkind
    zcdet=4.0_realkind/3.0_realkind
    !     
    !     initialize other temporary variables
    !     
    coef1 = gravit*864.0_realkind      ! conversion to millibars/day
    cappa=rair/cpair
    do i=kstart,kstop
       thetas(i)  = ts(i)*(pnot/ps(i))**cappa
       cck(i) = 0.0_realkind    
       clc(i) = 0.0_realkind
       mcbar(i) = 0.0_realkind
       dpsum(i) = 0.0_realkind
!cgjnewrh
       ptropo(i)=5000._realkind
       stropo(i)=5000._realkind/ps(i)
       ktropo_wmo(i)=1
       ktropo_bvf(i)=1
       zbvflim(i)=0._realkind
    enddo

    do k=1,plev
       do i=kstart,kstop
          qs(i,k) = (esat(temp(i,k))*epsilo)/ (pmid(i,k)-  &
                    ((1._realkind-epsilo)*esat(temp(i,k))))
! make theta_v
          ztv = temp(i,k)*(1.0_realkind + 0.61_realkind*q(i,k))
          theta(i,k)  = temp(i,k)*(pnot/pmid(i,k))**cappa
          thetav(i,k) = ztv*(pnot/pmid(i,k))**cappa
!
          rh(i,k)     = q(i,k)/qs(i,k)
          rh(i,k)     = max(rh(i,k),0.0_realkind)
          rh(i,k)     = min(rh(i,k),1.0_realkind)
          cloud(i,k)  = 0.0_realkind
          cldst(i,k)  = 0.0_realkind
          concld(i,k) = 0.0_realkind
          shalcld(i,k) = 0.0_realkind
          shalcld2(i,k) = 0.0_realkind
          xucldarray(i,k)  = 0.0_realkind
          rad_cloud(i,k)=0.0_realkind
          rhcrit(i,k)=0.0_realkind
          wzf(i,k)=gpot(i,k)*zginv
!parameter to modify rhminh to account for thicker layers in a dz metres sense
          zwrk1=min(500._realkind,max(0._realkind,         &
                                     (dzh(i,k)-500._realkind)))
          zwrk2=min(1._realkind,max(0.6_realkind,(1._realkind-       &
                               (0.4_realkind/500._realkind)*zwrk1)))
          rhmin(i,k)=max(0.5_realkind,(rhminh*zwrk2))
!
!v588          zwrk1=min(500._realkind,max(0._realkind,         &
!v588                                     (dzh(i,k)-650._realkind)))
!v588          zwrk2=min(1._realkind,max(0.65_realkind,(1._realkind-       &
!v588                               (0.35_realkind/500._realkind)*zwrk1)))
!v588          rhmin(i,k)=rhminh*zwrk2
!
       enddo
    enddo
!bvf
    do k=2,plev
       do i=kstart,kstop
          rho = pmid(i,k)/(rair*temp(i,k))
          bvf(i,k) = -rho*gravit*gravit*((theta(i,k)-theta(i,k-1))* &
                     rpdeli(i,k))/thetav(i,k)
       enddo
    enddo

!define tropopause as layer just above maximum bvf
     do k=plev-1,2,-1
        do i=kstart,kstop
         zbvfg(i)=bvf(i,k)-bvf(i,k+1)
           if(zbvfg(i)>zbvflim(i))then
             ktropo_bvf(i)=k
             zbvflim(i)=zbvfg(i)
          endif
!
          if(ktropo_wmo(i)==1)then
            zztr(i)=1000._realkind*(temp(i,k-1)-temp(i,k))/(wzf(i,k-1)-wzf(i,k))
            if(zztr(i)>zwmo)then	!lapse rate smaller than 2K/km
             zzcnt=0._realkind
             zztint=0._realkind
              do kk=k-1,2,-1
               if(zzcnt<3000._realkind)then
                zdz=wzf(i,kk-1)-wzf(i,kk)
                zzt=temp(i,kk-1)-temp(i,kk)
                zzcnt=zzcnt+zdz
                zztint=zztint+zzt
                 if(zzcnt>=3000._realkind)then
                   zztav=zztint*1000._realkind/zzcnt
                   if(zztav>zwmo)ktropo_wmo(i)=k-1
                 endif
               endif
              enddo
            endif
          endif
!
        enddo
     enddo
!cgjnewrh
    do i=kstart,kstop
!
      if(ktropo_bvf(i)>1 .and. ktropo_wmo(i)>1)then 
       ktt=max(ktropo_bvf(i),ktropo_wmo(i))
      elseif(ktropo_bvf(i)>1 .and. ktropo_wmo(i)==1)then
       ktt=ktropo_bvf(i)
      elseif(ktropo_bvf(i)==1 .and. ktropo_wmo(i)>1)then
       ktt=ktropo_wmo(i)
      else
       ktt=1
      endif
!
      ptropo(i)=min(35000._realkind,max(7000._realkind,pmid(i,ktt)))
      stropo(i)=min(0.35_realkind,(max(0.07_realkind,(ptropo(i)/ps(i)))))
!
      ztropo(i)=min(25000._realkind,max(10000._realkind,wzf(i,ktt)))
      ttropo(i)=min(250._realkind,max(150._realkind,temp(i,ktt)))
!      ptropo_out(i)=min(20000.,ptropo(i))
    enddo
    !     
    !     calculate mean convective motion throughout column (in units of mb/day)
    !     
    do k=1,plev-1
       do i=kstart,kstop
          mcbar(i) = mcbar(i) + max(cmfmc(i,k+1)*coef1,0.0_realkind)*pdel(i,k)
          dpsum(i) = dpsum(i) + pdel(i,k)
       enddo
    enddo
    !     
    !     gj    add in conv cld calc from xu & krueger mwr 1991,p342
    !     gj
    call cucoeff(plond,plev,kstart,kstop,pmid,pdel,cu0,cu1,cu2)
    !     gj
    !     gj   do not allow conv clouds where shallow conv occuring,
    !     gj   shallow conv clds only come thru large sale calc for now.
    !     gj
    do k = 1,plev-1
       do i = kstart,kstop
          if(shal_cgj(i)==1)then
             if(k<=kf_base(i) .and. k>=kf_top(i))then
!
                sr=(qvap(i,k)+qliq(i,k)+qice(i,k))/qs(i,k)
                if(sr>1.0_realkind .and. rh(i,k)<1.0_realkind)then
                   shalcld(i,k)=(sr-1.0_realkind)/(sr-rh(i,k))
!                   shalcld(i,k)=((sr-1.0_realkind)/(sr-rh(i,k)))**2
                elseif(sr<=1.0_realkind)then
                   shalcld(i,k)=0.0_realkind
                elseif(rh(i,k)>=1.0_realkind)then
                   shalcld(i,k)=1.0_realkind
                endif
!
                mcmbh=cmfmc(i,k)*gravit*36.0_realkind !kg/m2/s to mb/hour
                if(mcmbh>0.0001_realkind .and.&
                     (qvap(i,k)+qliq(i,k)+qice(i,k))>0.0_realkind)then
                   shalcld2(i,k)=cu0(i,k)+cu1(i,k)*log(mcmbh)+ &
                        cu2(i,k)*((log(mcmbh))**2)
                else
                   shalcld2(i,k)=0.0_realkind
                endif
!
             endif
             shalcld(i,k)=shalcld(i,k)+shalcld2(i,k)
          else
!
             mcmbh=cmfmc(i,k)*gravit*36.0_realkind !kg/m2/s to mb/hour
!260911             if(mcmbh>0.001_realkind.and.(qvap(i,k)+qliq(i,k)+qice(i,k))>0.0_realkind)then
             if(mcmbh>0.001_realkind)then
                if(k<=kf_base(i) .and. k>=kf_top(i))then
                   concld(i,k)=cu0(i,k)+cu1(i,k)*log(mcmbh)+ &
                        cu2(i,k)*((log(mcmbh))**2)
                else
                   concld(i,k)=0.0_realkind
                endif
             endif
!
             sr=(qvap(i,k)+qliq(i,k)+qice(i,k))/qs(i,k)
             if(sr>1.0_realkind .and. rh(i,k)<1.0_realkind)then
!cgj180811                xucldarray(i,k)=(sr-1.0_realkind)/(sr-rh(i,k))
                xucldarray(i,k)=((sr-1.0_realkind)/(sr-rh(i,k)))
             elseif(sr<=1.0_realkind)then
                xucldarray(i,k)=0.0_realkind
             elseif(rh(i,k)>=1.0_realkind)then
                xucldarray(i,k)=1.0_realkind
             endif
!
             zsigma=pmid(i,k)/ps(i)
             if (pmid(i,k)<5.e4_realkind .and. zsigma>stropo(i))then
             !cgjdo only detrainment bases convective clouds above 500hpa
!v588                zdetr=min(1._realkind,max(0._realkind,    &
!v588                         (40000._realkind-pmid(i,k))*5.e-5_realkind))
!v588                xucldarray(i,k) = max(xucldarray(i,k), &
!v588                     (min(rh(i,k),min(1.0_realkind,max(0.0_realkind,  &
!v588                     zdu(i,k)*zdetr*5.0e4_realkind)))))
                zdetr=min(8._realkind,max(0._realkind,    &
                         ((50000._realkind-pmid(i,k))*3.e-4_realkind)**zcdet))
                xucldarray(i,k) = max(xucldarray(i,k), &
                     (min(1.0_realkind,max(0.0_realkind,  &
                     zdu(i,k)*zdetr*5.0e4_realkind))))
             endif
!
             concld(i,k)=xucldarray(i,k)+concld(i,k)
             concld(i,k)=max(0.0_realkind,min(1.0_realkind,concld(i,k)))
          endif
!
          if(shalcld(i,k)>0._realkind)then
             if( (concld(i,k)+shalcld(i,k)) > 1.0_realkind)then
                concld(i,k)=concld(i,k)/(concld(i,k)+shalcld(i,k))
                shalcld(i,k)=shalcld(i,k)/(concld(i,k)+shalcld(i,k))
             endif
          endif
!
       enddo
    enddo
    !     
    !     
    !     evaluate effective column-integrated convective cloud cover using
    !     random overlap assumption (for diagnostic purposes only)
    !     
    do i=kstart,kstop
       clrsky(i) = 1.0_realkind
    enddo
    do k=plev,1,-1
       do i=kstart,kstop
          clrsky(i) = clrsky(i)*(1.0_realkind - concld(i,k))
!
          zsigma=pmid(i,k)/ps(i)
          if(zsigma>=0.875_realkind)then
             zwrk2=((zsigma-0.875_realkind)/0.125_realkind)**2
             rhcrit(i,k)=rhmin(i,k)+( (1.-rhmin(i,k))*zwrk2 )
          elseif(zsigma>=stropo(i) .and. zsigma<0.875_realkind)then	!below tropopause above ~875hpa
             rhcrit(i,k)=rhmin(i,k)
          else
             rhcrit(i,k)=rhmax
          endif
!
          if(zsigma>0.9)then
           if(temp(i,k)<273.15_realkind .and. temp(i,k)>233.15_realkind)then
              zwrk1=abs((temp(i,k)-253.15_realkind))-20._realkind
              zwrk2=zwrk1*zwrk1
              zwrk3=max(0._realkind,min(0.08_realkind,(zwrk2*0.08_realkind/400._realkind)))
              rhcrit(i,k)=rhcrit(i,k)-zwrk3
           endif
          endif
       enddo
    enddo
    do i=kstart,kstop
       clc(i) = 1.0_realkind - clrsky(i)
    enddo
    !     
    !     ****** compute layer cloudiness ******
    !     
    premib = 750.e2_realkind
    pretop = 1.0e2_realkind            ! top of cloud layer is at 1 mb
    !     
    !     find most stable level below 750 mb for evaluating stratus regimes
    !     
    do i=kstart,kstop
       dtdpmn(i) = 0.0_realkind     
!cgj211011       kdthdp(i) = 1
       kdthdp(i) = 0
       dthtdp(i,1) = 0.0_realkind
    enddo
    do k=2,plev-2
       do i=kstart,kstop
          if (pmid(i,k)>=premib) then
             dthdp = 100.0_realkind*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
          else
             dthdp = 0.0_realkind
          endif
          if (dthdp<dtdpmn(i)) then
             dtdpmn(i) = dthdp
             kdthdp(i) = k    ! index of interface of max inversion
          endif
          dthtdp(i,k) = dthdp !dthtdp(i,k) & dtdpmn(i) are both set equal to max negative theta jump
			      !negative meaning theta(k)-theta(k-1)
       enddo
    enddo
    do k=plev-1,plev
       do i=kstart,kstop
          if (0.0_realkind<dtdpmn(i)) then
             dtdpmn(i) = 0.0_realkind	!if dtdpmn(i) is positive (theta(k) always > theta(k-1)
          endif
          dthtdp(i,k) = 0.0_realkind
       enddo
    enddo
    !     
    !     bvf => brunt-vaisalla frequency (approx. 1-sided diff.)
    !     this stability measure is used to set a local relative humidity 
    !     threshold when evaluating the fractional area of layered cloud
    !     
    do 10 k=2,plev
       kp1 = min(k + 1,plev)
       do i=kstart,kstop
          if (dthtdp(i,k)>dtdpmn(i)) then
             dthtdp(i,k) = 0.0_realkind
          endif
          cldbnd(i) = pmid(i,k)>=pretop
          lol(i) = pmid(i,k)>=premib
!
          pdep(i)=ps(i)-pmid(i,k)
          if(pdep(i)<=30000.0_realkind)then
             lol(i)=.true.
             cldbnd(i)=.true.
          elseif(pdep(i)>30000.0_realkind .and. pmid(i,k)>=pretop)then
             lol(i)=.false.
             cldbnd(i)=.true.
          else
             lol(i)=.false.
             cldbnd(i)=.false.
          endif
!
! Control use of bvf*rbvflim
!
          zzm=gpot(i,k)*zginv
          if(zzm<=pblh(i))then
            if(temp(i,k)<273.16_realkind .and. zzm<1750._realkind)then
              zbvf=0._realkind
             else
              zbvf=1._realkind
            endif
          else				!above PBL
            if(ptropo(i)<15000._realkind)then	!deep tropics
             zbvf=0._realkind
            else
             zbvf=1._realkind
            endif
!            if(pmid(i,k)>30000._realkind)zbvf=1._realkind
          endif
          bvf(i,k)=min(bvf(i,k),0.0006_realkind)
!
          rhlim(i,k)=rhcrit(i,k)
!     
          if (cldbnd(i)) then
             rhlim(i,k) = rhmax - (1.0_realkind-rhcrit(i,k)) &
                          *(1.0_realkind-min(1.0_realkind,max(0.0_realkind,bvf(i,k)*rbvflim*zbvf)))
             rhden = 1.0_realkind - rhlim(i,k)
!             if(temp(i,k)<228.0_realkind)then
!                wrk1=1.0_realkind+0.15_realkind*(228.0_realkind-temp(i,k))
!                wrk2=(rhmax-rhlim(i,k))*(1.0_realkind-1.0_realkind/wrk1)
!                rhlim(i,k)=min(rhlim(i,k)+wrk2,rhmax)
!             endif
          else
             rhlim(i,k) = rhmax
             rhden = 0.001_realkind
          endif
          rhdif = (rh(i,k) - rhlim(i,k))/rhden
          cld9(i) = min(0.999_realkind,(max(rhdif,0.0_realkind))**2)
          !     
          if (lol(i)) then
             rhlim(i,k) = rhmax - (1.0_realkind-rhcrit(i,k)) &
                         *(1.0_realkind-min(1.0_realkind,max(0.0_realkind,bvf(i,k)*rbvflim*zbvf)))
             rhdif2 = (rh(i,k) - rhlim(i,k))/(1.0_realkind-rhlim(i,k))
             cld8(i) = min(0.999_realkind,(max(rhdif2,0.0_realkind))**2)
          else
             cld8(i) = cld9(i)
          endif
       enddo
       !     
       !     final evaluation of layered cloud fraction
       !     
       do i=kstart,kstop
          if (lol(i)) then
             cloud(i,k) = cld8(i)
          else                ! middle and high level cloud 
             if( cldbnd(i) ) then
                cloud(i,k) = cld9(i)
             else
                cloud(i,k) = 0.0_realkind
             endif
          endif
       enddo
10  enddo                     ! k=2,plev-1
    !     
    !     add in the marine strat
    !     
    do i=kstart,kstop
       if(kdthdp(i)/=0)then
!          pdepst(i)=ps(i)-pmid(i,kdthdp(i))
!          if(rh(i,kdthdp(i))<0.7_realkind .or.pdepst(i)>20000.0_realkind)kdthdp(i)=0
!cgj2408          zfrice=frice(i)*(1.0_realkind-frland(i))
!cgj2408          if(abs(zfrice-1.0_realkind)<1.e-14_realkind)kdthdp(i)=0

       endif
       !     gj080502

       if (kdthdp(i)/=0) then
          k = kdthdp(i)
          kp1 = min(k+1,plev)
          strat = min(1.0_realkind,max(0.0_realkind,(theta(i,k700(i))-thetas(i))*.057_realkind-.5573_realkind))
          !     
          if (oro(i)==0 .and. dthtdp(i,k)<=-0.125_realkind ) then
             cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
             cloud(i,k) = max(cloud(i,k),cldst(i,k))
          else
             cldst(i,k) = 0.0_realkind
          endif
       endif
    enddo
    !     
    !     merge convective and layered cloud fraction for total cloud
    !     
    do k=1,plev
       do i=kstart,kstop
          lscld(i,k)=max(0.0_realkind,min(0.999_realkind,cloud(i,k))) 
          if(shal_cgj(i)==1 .and.shalcld(i,k)>0.0_realkind)then
             cloud(i,k)=shalcld(i,k)
          else
             cloud(i,k) = max(0.0_realkind,min(0.999_realkind,cloud(i,k)))
          endif
          if (rh(i,k)>0.99_realkind) then
             cloud(i,k) = max(0.01_realkind,cloud(i,k))
             lscld(i,k) = max(0.01_realkind,lscld(i,k))
          endif
!
          rad_cloud(i,k)=max(0.0_realkind,min(1.0_realkind,(concld(i,k)+cloud(i,k)) ))
          if((concld(i,k)+cloud(i,k))>1.0_realkind)then
             concld(i,k)=concld(i,k)/(concld(i,k)+cloud(i,k))
             cloud(i,k)=cloud(i,k)/(concld(i,k)+cloud(i,k))
          endif
          !     gj
       enddo
    enddo
!
    do i = kstart,kstop
       shal2d(i)=0.0_realkind
       cu2d(i)=0.0_realkind
       cov2d(i)=1.0_realkind
       cov2d_odp(i)=1.0_realkind
    enddo
    !     
    !     maximum/random overlaping
    !     
    do k=2,plev
       do i=kstart,kstop
          cov2d(i)=cov2d(i) *(1.0_realkind - max(rad_cloud(i,k-1), &
               rad_cloud(i,k))) /(1.0_realkind- min(rad_cloud(i,k-1),0.99_realkind))
          cu2d(i)=cu2d(i)*(1.0_realkind-max(concld(i,k-1),concld(i,k))) &
               /(1.0_realkind- min(concld(i,k-1),0.99_realkind))
          shal2d(i)=shal2d(i)*(1.0_realkind - max(shalcld(i,k-1), &
               shalcld(i,k))) /(1.0_realkind- min(shalcld(i,k-1),0.99_realkind))
       enddo
    enddo

    do i=kstart,kstop
       cov2d(i)=1.0_realkind-cov2d(i)
       cu2d(i)=1.0_realkind-cu2d(i)
       shal2d(i)=1.0_realkind-shal2d(i)
    enddo
    return
  end subroutine cldfrc
  !gj------------------------------------------------------------------------------
  !gj
  subroutine axucld(plond,plev,kstart,kstop,             &
       temp,q,cw,pmid,ps,rhcrit,lesat,lscov, &
       lscloud,concld,                       &
       shalcld,shal_cgj,cov2d)
    implicit none
    !
    !gj   xu-randall jas 1996 approach to calculate cloud fraction as
    !gj   a function of rh,cw & qsat. through a timestream we now 
    !gj   in t=t-1 calculate a slingo cloud fn(rh) this controls 
    !gj   the t=t-1 prediction of cloud water & precipitation.
    !gj   then in t=t we recalculate a cloud fraction here.
    !gj   this is the cloud fraction used in radiation calculations
    !gj   in essence this means that condensation/evaporation of ls
    !gj   cloud water occurs/is controlled by a threshold rh which the
    !gj   slingo cloud is a surrogate for. and the subsequent "radiativel
    !gj   active cloud" is formed after condensation and is controlled by
    !gj   rh,cloud water & qsat.
    !gj
    !gj
    integer:: plond,plev,kstart,kstop,i,k
    integer:: shal_cgj(plond)
    !gj
    real(kind=realkind):: temp(plond,plev)
    real(kind=realkind):: q(plond,plev)
    real(kind=realkind):: cw(plond,plev)
    real(kind=realkind):: pmid(plond,plev)
    real(kind=realkind):: concld(plond,plev)
    real(kind=realkind):: rhcrit(plond,plev)
    real(kind=realkind):: lscov(plond,plev)
    real(kind=realkind):: shalcld(plond,plev)
    real(kind=realkind):: ps(plond)
    !gj
    !gj   output
    !gj
    real(kind=realkind):: rad_cloud(plond,plev)
    real(kind=realkind):: lscloud(plond,plev)
    real(kind=realkind):: cov2d(plond)
    !gj
    real(kind=realkind):: rh(plond,plev)
    real(kind=realkind):: qs(plond,plev)
    real(kind=realkind):: ess,pss
    logical:: lesat
    !


    !gj   calculate saturation specific humidity and relative humidity
    !gj
    do k=1,plev
       do i=kstart,kstop
          !gj080304
          if(lesat)then
             if(cw(i,k)<=0.0_realkind)then
                ess=esatw(temp(i,k))
             else
                ess=(lscov(i,k)*esat(temp(i,k)))+((1.0_realkind-lscov(i,k))*esatw(temp(i,k)))
             endif
             !gj080304
          else
             ess=esat(temp(i,k))
          endif
          pss=pmid(i,k)
          qs(i,k)=(ess*epsilo)/(pss-((1.0_realkind-epsilo)*ess))
          rh(i,k)     = q(i,k)/qs(i,k)
          rh(i,k)     = max(rh(i,k),0.0_realkind)
          rh(i,k)     = min(rh(i,k),1.0_realkind)
          !kw add cloud initialisation
          lscloud(i,k)=0.0_realkind
       enddo
    enddo
    !gj
    call xucld(plond,plev,kstart,kstop,rh,qs,cw, &
         pmid,ps,rhcrit,shalcld,shal_cgj,lscloud)
    !gj
    !
    ! merge convective and layered cloud fraction for total cloud
    !
    do k=1,plev
       do i=kstart,kstop
          !gj
          lscloud(i,k) = max(0.0_realkind,min(0.999_realkind,lscloud(i,k)))
          !gj
          if (rh(i,k)>0.99_realkind) then
             lscloud(i,k) = max(0.01_realkind,lscloud(i,k))
          endif
          !gj
          rad_cloud(i,k)=max(0.0_realkind,min(1.0_realkind,(concld(i,k)+lscloud(i,k)) ))
          !gj
          !gj    it is possible that cloud+concld can become greater than 1.
          !gj    if so scale the clouds so that the total value is 1, whilst
          !gj    retaining the relative proportions of each cloud type.
          !gj
          if((concld(i,k)+lscloud(i,k))>1.0_realkind)then
             concld(i,k)=concld(i,k)/(concld(i,k)+lscloud(i,k))
             lscloud(i,k)=lscloud(i,k)/(concld(i,k)+lscloud(i,k))
          endif
          !gj
       enddo
    enddo
    !vo +++
    !vo adding hirlam algorithm for calculation of 2-dimensional cloud cover
    !vo
    do i = kstart,kstop
       !
       !  put cov2d to 1. here, see below
       !
       cov2d(i)=1.0_realkind
    enddo
    !
    !  maximum/random overlaping
    !
    do k=2,plev
       do i=kstart,kstop
          cov2d(i)=cov2d(i)*(1.0_realkind - max(rad_cloud(i,k-1), &
               rad_cloud(i,k))) /(1.0_realkind- min(rad_cloud(i,k-1),0.99_realkind))
       enddo
    enddo

    do i=kstart,kstop
       cov2d(i)=1.0_realkind-cov2d(i)
    enddo
    !vo ---
    !
    return
  end subroutine axucld
  !gj
  !gj-----------------------------------------------------------------------------------
  !gj
  subroutine xucld(plond,plev,kstart,kstop,rh,qs,cw, &
       pmid,ps,rhcrit,shalcld,shal_cgj,lscloud)
    implicit none
    !gj
    integer:: plond,plev,kstart,kstop,i,k
    integer:: shal_cgj(plond)
    real(kind=realkind):: cw(plond,plev)
    real(kind=realkind):: qs(plond,plev)
    real(kind=realkind):: rh(plond,plev)
    real(kind=realkind):: rhcrit(plond,plev)
    real(kind=realkind):: lscloud(plond,plev)
    real(kind=realkind):: pmid(plond,plev)
    real(kind=realkind):: shalcld(plond,plev)
    real(kind=realkind):: ps(plond)
    real(kind=realkind):: cwxu,qsxu,wrk1,wrk2,wrk3,wrk4,gamma,alpha,pcoeff
    real(kind=realkind):: zwrk1,za,zalfa,zfac
    !gj
    !cc      write(*,*)'entered xucld!'
    za=0.028_realkind
    !gj      gamma=0.49
    !gj      alpha=100.
    pcoeff=0.25_realkind
    zfac=1.e-6_realkind
    !gj

    !gj
    do k=1,plev
       do i=kstart,kstop
          !gj300408
          if(pmid(i,k)<=20000.0_realkind)then
             gamma=0.54_realkind
          else
             gamma=0.54_realkind+((20000.0_realkind-pmid(i,k))*zfac)
          endif
          !gj300408
          !gj040604
          zwrk1=ps(i)-pmid(i,k)
          if(zwrk1<=20000.0_realkind)then
             zalfa=100.0_realkind
          else
             zalfa=100.0_realkind+((zwrk1-20000.0_realkind)*za)
          endif
          alpha=max(100.0_realkind,min(zalfa,650.0_realkind))

          cwxu=cw(i,k)
          qsxu=qs(i,k)
          !ghshallowcloud
          ! europe
          if(rh(i,k)>=rhcrit(i,k) .and. rh(i,k)<1.0_realkind)then
             !gj          if(rh(i,k)<1.0)then
             ! africa tuned
             !          if(rh(i,k)<1.0)then
             wrk1=-alpha*cwxu
             wrk2=((1.0_realkind-rh(i,k))*qsxu)**gamma
             wrk3=wrk1/wrk2
             wrk4=rh(i,k)**pcoeff
             lscloud(i,k)=wrk4*(1.0_realkind-exp(wrk3))
          elseif(rh(i,k)>=1.0_realkind)then
             lscloud(i,k)=1.0_realkind
          endif
          lscloud(i,k)=min(1.0_realkind,max(lscloud(i,k),0.0_realkind))
          !gj
          if(shalcld(i,k)>0.0_realkind)lscloud(i,k)=max(lscloud(i,k),shalcld(i,k))
       enddo
    enddo
    !gj
    return
  end subroutine xucld



  subroutine cucoeff(plond,plev,kstart,kstop,pmid,pdel,c0,c1,c2)
    implicit none
    integer:: i,k,plond,plev,kstart,kstop
    real(kind=realkind):: pmid(plond,plev),pdel(plond,plev),c0(plond,plev), &
         c1(plond,plev),c2(plond,plev)
    do k=1,plev
       do i=kstart,kstop
          if(pdel(i,k)<=7500.0_realkind)then	!"Thin" layers
             if(pmid(i,k)<=50000.0_realkind)then
                c0(i,k)=0.0131_realkind
                c1(i,k)=0.0123_realkind
                c2(i,k)=0.0030_realkind
             elseif(pmid(i,k)<=80000.0_realkind .and. pmid(i,k)>50000.0_realkind)then
                c0(i,k)=0.0363_realkind
                c1(i,k)=0.0226_realkind
                c2(i,k)=0.0054_realkind
             else
                c0(i,k)=0.0116_realkind
                c1(i,k)=0.0103_realkind
                c2(i,k)=0.0025_realkind
             endif
             !gj
          else			!"Thick" layers
             if(pmid(i,k)<=50000.0_realkind)then
                c0(i,k)=0.0708_realkind
                c1(i,k)=0.0941_realkind
                c2(i,k)=0.0233_realkind
             elseif(pmid(i,k)<=80000.0_realkind .and. pmid(i,k)>50000.0_realkind)then
                c0(i,k)=0.0722_realkind
                c1(i,k)=0.0760_realkind
                c2(i,k)=0.0254_realkind
             else
                c0(i,k)=0.0337_realkind
                c1(i,k)=0.0453_realkind
                c2(i,k)=0.0237_realkind
             endif
             !gj
          endif
       enddo
    enddo

    return
  end subroutine cucoeff


  subroutine rhcrak(plond,plev,kstart,kstop,rhminh,p,ps,rhcrit)

    implicit none
    integer:: i,k,plond,plev,kstart,kstop
    real(kind=realkind):: rh1,rh2,rh3,rh4,rh5,rh6,rh7,rh8,rh9,rhminh
    real(kind=realkind):: p(plond,plev),rhcrit(plond,plev),ps(plond),rhmax
    !
    rhmax=0.975_realkind
    !
    !gj
    rh1=rhmax
    rh2=0.50_realkind	!0.50
    rh3=rhminh-0.3_realkind	!0.40
    rh4=rhminh-0.3_realkind	!0.40
    rh5=rhminh-0.3_realkind	!0.40
    rh6=rhminh-0.15_realkind	!0.55
    rh7=rhminh-0.08_realkind	!0.62
    rh8=rhminh-0.05_realkind	!0.65
    rh9=rhminh	!0.7
    !
    !gj
    !gjnew       rh1=rhmax
    !gjnew       rh2=0.70	!0.7
    !gjnew       rh3=rhminh-0.1	!0.7
    !gjnew       rh4=rhminh-0.1	!0.7
    !gjnew       rh5=rhminh-0.1	!0.7
    !gjnew       rh6=rhminh-0.1	!0.7
    !gjnew       rh7=rhminh-0.05	!0.75
    !gjnew       rh8=rhminh-0.05	!0.75
    !gjnew        rh9=rhminh	!0.8
    !
    !gj
    do k=1,plev
       do i=kstart,kstop
          if(p(i,k)<15000.0_realkind)rhcrit(i,k)=rh1
          if(p(i,k)>=15000._realkind .and. p(i,k)<20000.0_realkind)           &
               rhcrit(i,k)=(rh1*(20000._realkind-p(i,k))/5000._realkind)+        &
               (rh2*(p(i,k)-15000._realkind)/5000._realkind)                    
          if(p(i,k)>=20000._realkind .and. p(i,k)<30000.0_realkind)            &
               rhcrit(i,k)=(rh2*(30000._realkind-p(i,k))/10000._realkind)+       &
               (rh3*(p(i,k)-20000._realkind)/10000._realkind)                   
          if(p(i,k)>=30000._realkind .and. p(i,k)<40000.0_realkind)            &
               rhcrit(i,k)=(rh3*(40000._realkind-p(i,k))/10000._realkind)+       &
               (rh4*(p(i,k)-30000._realkind)/10000._realkind)                   
          if(p(i,k)>=40000._realkind .and. p(i,k)<50000.0_realkind)            &
               rhcrit(i,k)=(rh4*(50000._realkind-p(i,k))/10000._realkind)+       &
               (rh5*(p(i,k)-40000._realkind)/10000._realkind)                   
          if(p(i,k)>=50000._realkind .and. p(i,k)<60000.0_realkind)            &
               rhcrit(i,k)=(rh5*(60000._realkind-p(i,k))/10000._realkind)+       &
               (rh6*(p(i,k)-50000._realkind)/10000._realkind)                   
          if(p(i,k)>=60000._realkind .and. p(i,k)<70000.0_realkind)            &
               rhcrit(i,k)=(rh6*(70000._realkind-p(i,k))/10000._realkind)+       &
               (rh7*(p(i,k)-60000._realkind)/10000._realkind)                   
          if(p(i,k)>=70000._realkind .and. p(i,k)<80000.0_realkind)            &
               rhcrit(i,k)=(rh7*(80000._realkind-p(i,k))/10000._realkind)+       &
               (rh8*(p(i,k)-70000._realkind)/10000._realkind)                   
          if(p(i,k)>=80000._realkind .and. p(i,k)<90000.0_realkind)            &
               rhcrit(i,k)=(rh8*(90000._realkind-p(i,k))/10000._realkind)+       &
               (rh9*(p(i,k)-80000._realkind)/10000._realkind)                   
          if(p(i,k)>=90000.0_realkind)rhcrit(i,k)=rhminh                 
          if(p(i,k)/ps(i)>0.9_realkind)then                            
             rhcrit(i,k)=rhcrit(i,k)+(rhmax-rhcrit(i,k))*         &
                  (1.0_realkind-10._realkind*(1.0_realkind-p(i,k)/ps(i)))
             rhcrit(i,k)=min(rhcrit(i,k),rhmax)
          endif
          if(rhcrit(i,k)>1.0_realkind .or. rhcrit(i,k)<=0.0_realkind)then
             call stop_program('stop! rhcrit is wrong!!!!!!!!!!!!!!!!!')
          endif
       enddo
    enddo
    !
    return
  end subroutine rhcrak



end module cldfrcmod
