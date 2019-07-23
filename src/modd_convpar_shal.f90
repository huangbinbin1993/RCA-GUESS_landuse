module modd_convpar_shal
  !!  modd_convpar_shal - declaration of convection constants 
  !!
  !!    purpose
  !!    -------
  !!      the purpose of this declarative module is to declare  the 
  !!      constants in the deep convection parameterization.    
  !!
  !!
  !!  implicit arguments
  !!    ------------------
  !!      none 
  !!
  !!    reference
  !!    ---------
  !!      book2 of documentation of meso-nh (modd_convpar_shal)
  !!          
  !!    author
  !!    ------
  !!      p. bechtold   laboratoire d'aerologie
  !!
  !!    modifications
  !!    -------------
  !!      original    26/03/96                      
  !!   last modified  04/10/98
  !-------------------------------------------------------------------------------

  !       0.   declarations

  use modd_cst
  use modd_convparext
  use confys,only:epsilo
  use decomp,only:realkind
  implicit none 
  private
  real(kind=realkind), save :: xa25= 625.e6_realkind         ! 25 km x 25 km reference grid area

  real(kind=realkind), save :: xcrad =  50._realkind       ! cloud radius 
  real(kind=realkind), save :: xctime_shal= 10800._realkind ! convective adjustment time
  real(kind=realkind), save :: xcdepth = 0.3e3_realkind   ! minimum necessary cloud depth
  real(kind=realkind), save :: xcdepth_d= 3.0e3_realkind   ! maximum allowed cloud thickness
  real(kind=realkind), save :: xdtpert = .2_realkind   ! add small temp perturb. at lcl
  real(kind=realkind), save :: xentr  = 0.03_realkind      ! entrainment constant (m/pa) = 0.2 (m)  

  real(kind=realkind), save :: xzlcl  = 1.5e3_realkind      ! maximum allowed allowed height 
  ! difference between departure level and surface
  real(kind=realkind), save :: xzpbl  = 50.e2_realkind ! minimum mixed layer depth to sustain convection
!  real(kind=realkind), save :: xwtrig      ! constant in vertical velocity trigger


  real(kind=realkind), save :: xnhgam= 1.3333_realkind       ! accounts for non-hydrost. pressure 
  ! in buoyancy term of w equation
  ! = 2 / (1+gamma)
!  real(kind=realkind), save :: xtfrz1= 268.16      ! begin of freezing interval
!  real(kind=realkind), save :: xtfrz2= 248.16      ! end of freezing interval


  real(kind=realkind), save :: xstabt= 0.75_realkind      ! factor to assure stability in  fractional time
  ! integration, routine convect_closure
  real(kind=realkind), save :: xstabc= 0.95_realkind      ! factor to assure stability in cape adjustment,
  !  routine convect_closure
  public  shallow_convection
contains
  subroutine shallow_convection( klon, klev, kidia, kfdia, kbdia, ktdia,        &
       pdtconv, kice, osettadj, ptadjs,               &
       indexcv,gcrow,wstar,                           &  !jyj
       ppabst, pzz,                                   &
       ptt, prvt, prct, prit, pwt,                    &
       ptten, prvten, prcten, priten,                 &
       kcltop, kclbas, pumf,                          &
       qvap, qliq, qice,                              &
       och1conv, kch1, pch1, pch1ten, pscdep            )
    !   ############################################################################
    !
    !! monitor routine to compute all convective tendencies by calls
    !!     of several subroutines.
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the convective
    !!      tendencies. the routine first prepares all necessary grid-scale
    !!      variables. the final convective tendencies are then computed by
    !!      calls of different subroutines.
    !!
    !!
    !!  method
    !!    ------
    !!      we start by selecting convective columns in the model domain through
    !!      the call of routine trigger_funct. then, we allocate memory for the
    !!      convection updraft and downdraft variables and gather the grid scale
    !!      variables in convective arrays. 
    !!      the updraft and downdraft computations are done level by level starting
    !!      at the  bottom and top of the domain, respectively.
    !!      all computations are done on mnh thermodynamic levels. the depth
    !!      of the current model layer k is defined by dp(k)=p(k-1)-p(k)
    !!      
    !!     
    !!
    !!    external
    !!    --------
    !!    convect_trigger_shal
    !!    convect_satmixratio
    !!    convect_updraft_shal
    !!        convect_condens
    !!        convect_mixing_funct
    !!    convect_closure_shal
    !!        convect_closure_thrvlcl
    !!        convect_closure_adjust_shal
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_cst
    !!          xg                   ! gravity constant
    !!          xpi                  ! number pi
    !!          xp00                 ! reference pressure
    !!          xrd, xrv             ! gaz  constants for dry air and water vapor
    !!          xcpd, xcpv           ! specific heat for dry air and water vapor
    !!          xrholw               ! density of liquid water
    !!          xalpw, xbetaw, xgamw ! constants for water saturation pressure
    !!          xtt                  ! triple point temperature
    !!          xlvtt, xlstt         ! vaporization, sublimation heat constant
    !!          xcl, xci             ! specific heat for liquid water and ice
    !!
    !!      module modd_convparext
    !!          jcvexb, jcvext       ! extra levels on the vertical boundaries
    !!
    !!      module modd_convpar
    !!          xa25                 ! reference grid area
    !!          xcrad                ! cloud radius
    !!
    !!         
    !!    reference
    !!    ---------
    !!
    !!      bechtold, 1997 : meso-nh scientific  documentation (31 pp)
    !!      fritsch and chappell, 1980, j. atmos. sci., vol. 37, 1722-1761.
    !!      kain and fritsch, 1990, j. atmos. sci., vol. 47, 2784-2801.
    !!      kain and fritsch, 1993, meteor. monographs, vol. 24, 165-170.
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    26/03/96 
    !!   peter bechtold 15/11/96 replace theta_il by enthalpy
    !!         "        10/12/98 changes for arpege
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !

    !
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    !
    integer,                    intent(in) :: klon     ! horizontal dimension
    integer,                    intent(in) :: klev     ! vertical dimension
    integer,                    intent(in) :: kidia    ! value of the first point in x
    integer,                    intent(in) :: kfdia    ! value of the last point in x
    integer,                    intent(in) :: kbdia    ! vertical  computations start at
    !                                                  ! kbdia that is at least 1
    integer,                    intent(in) :: ktdia    ! vertical computations can be
    ! limited to klev + 1 - ktdia
    ! default=1
    real(kind=realkind),                       intent(in) :: pdtconv  ! interval of time between two
    ! calls of the deep convection
    ! scheme
    integer,                    intent(in) :: kice     ! flag for ice ( 1 = yes, 
    !                0 = no ice )
    logical,                    intent(in) :: osettadj ! logical to set convective
    ! adjustment time by user
    real(kind=realkind),                       intent(in) :: ptadjs   ! user defined adjustment time
    real(kind=realkind), dimension(klon,klev), intent(in) :: ptt      ! grid scale temperature at t
    real(kind=realkind), dimension(klon,klev), intent(in) :: prvt     ! grid scale water vapor "
    real(kind=realkind), dimension(klon,klev), intent(in) :: prct     ! grid scale r_c  "
    real(kind=realkind), dimension(klon,klev), intent(in) :: prit     ! grid scale r_i "
    real(kind=realkind), dimension(klon,klev), intent(in) :: pwt      ! grid scale vertical 
    ! velocity (m/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppabst   ! grid scale pressure at t
    real(kind=realkind), dimension(klon,klev), intent(in) :: pzz      ! height of model layer (m) 
    !   
    real(kind=realkind), dimension(klon,klev), intent(out):: qvap     ! cloud plumes total water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: qliq     ! cloud plumes cloud water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: qice     ! cloud plumes cloud ice   (kg/kg)

    real(kind=realkind), dimension(klon,klev), intent(inout):: ptten  ! convective temperature 
    ! tendency (k/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: prvten ! convective r_v tendency (1/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: prcten ! convective r_c tendency (1/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: priten ! convective r_i tendency (1/s)
    real(kind=realkind), dimension(klon),      intent(in)   :: gcrow  ! jyj land sea mask
    real(kind=realkind), dimension(klon),      intent(in)   :: wstar  ! jyj wstar free convection velocity
    integer, dimension(klon),   intent(inout):: indexcv! convection index (deep=2, shallow=1,none=0)  !jyj
    integer, dimension(klon),   intent(inout):: kcltop ! cloud top level
    integer, dimension(klon),   intent(inout):: kclbas ! cloud base level
    ! they are given a value of
    ! 0 if no convection
    real(kind=realkind), dimension(klon,klev), intent(inout):: pumf   ! updraft mass flux (kg/s m2)
    !
    logical,                    intent(in) :: och1conv ! include tracer transport
    integer,                    intent(in) :: kch1     ! number of species
    real(kind=realkind), dimension(klon,klev,kch1), intent(in) :: pch1! grid scale chemical species
    real(kind=realkind), dimension(klon,klev,kch1), intent(inout):: pch1ten! species conv. tendency (1/s)
    !cgj
    real(kind=realkind), dimension(klon),      intent(inout)   :: pscdep  ! shallow conv cloud depth (m)

    !
    !
    !       0.2   declarations of local fixed memory variables :
    !
    integer  :: itest, iconv, iconv1    ! number of convective columns
    integer  :: iib, iie                ! horizontal loop bounds
    integer  :: ikb, ike                ! vertical loop bounds
    integer  :: iks                     ! vertical dimension
    integer  :: ji, jl                  ! horizontal loop index
    integer  :: jn                      ! number of tracers
    integer  :: jk, jkp, jkm            ! vertical loop index
    integer  :: iftsteps                ! only used for chemical tracers
    real(kind=realkind)     :: zeps, zepsa, zepsb      ! r_d / r_v, r_v / r_d, xcpv / xcpd - zepsa
    real(kind=realkind)     :: zcpord, zrdocp          ! c_p/r_d,  r_d/c_p
    !
    logical, dimension(klon, klev)     :: gtrig3 ! 3d logical mask for convection 
    logical, dimension(klon)           :: gtrig  ! 2d logical mask for trigger test
    real(kind=realkind), dimension(klon,klev)         :: ztht, zsthv, zsthes  ! grid scale theta, theta_v
    real(kind=realkind), dimension(klon)              :: ztime  ! convective time period
    real(kind=realkind), dimension(klon)              :: zwork2, zwork2b ! work array 
    !
    !
    !       0.2   declarations of local allocatable  variables :
    !
    integer, dimension(:),allocatable  :: idpl    ! index for parcel departure level
    integer, dimension(:),allocatable  :: ipbl    ! index for source layer top
    integer, dimension(:),allocatable  :: ilcl    ! index for lifting condensation level 
    integer, dimension(:),allocatable  :: ietl    ! index for zero buoyancy level
    integer, dimension(:),allocatable  :: ictl    ! index for cloud top level
    integer, dimension(:),allocatable  :: ilfs    ! index for level of free sink
    !
    integer, dimension(:), allocatable :: isdpl   ! index for parcel departure level
    integer, dimension(:),allocatable  :: ispbl   ! index for source layer top
    integer, dimension(:), allocatable :: islcl   ! index for lifting condensation level 
    real(kind=realkind), dimension(:), allocatable    :: zsthlcl ! updraft theta at lcl
    real(kind=realkind), dimension(:), allocatable    :: zstlcl  ! updraft temp. at lcl
    real(kind=realkind), dimension(:), allocatable    :: zsrvlcl ! updraft rv at lcl
    real(kind=realkind), dimension(:), allocatable    :: zswlcl  ! updraft w at lcl
    real(kind=realkind), dimension(:), allocatable    :: zszlcl  ! lcl height
    real(kind=realkind), dimension(:), allocatable    :: zsthvelcl! envir. theta_v at lcl
    real(kind=realkind), dimension(:), allocatable    :: zsdxdy  ! grid area (m^2)
    real(kind=realkind), dimension(:), allocatable    :: zswstar ! jyj wstar free convection velocity
    real(kind=realkind), dimension(:), allocatable    :: zsland  ! jyj land sea mask
    !cgj
    real(kind=realkind), dimension(:), allocatable    :: zsscdep  ! shallow conv cloud depth (m)
    !
    ! grid scale variables
    real(kind=realkind), dimension(:,:), allocatable  :: zz      ! height of model layer (m) 
    real(kind=realkind), dimension(:,:), allocatable  :: zpres   ! grid scale pressure
    real(kind=realkind), dimension(:,:), allocatable  :: zdpres  ! pressure difference between 
    ! bottom and top of layer (pa)
    real(kind=realkind), dimension(:,:), allocatable  :: zw      ! grid scale vertical velocity on theta grid
    real(kind=realkind), dimension(:,:), allocatable  :: ztt     ! temperature
    real(kind=realkind), dimension(:,:), allocatable  :: zth     ! grid scale theta     
    real(kind=realkind), dimension(:,:), allocatable  :: zthv    ! grid scale theta_v     
    real(kind=realkind), dimension(:,:), allocatable  :: zthl    ! grid scale enthalpy (j/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zthes, zthest ! grid scale saturated theta_e
    real(kind=realkind), dimension(:,:), allocatable  :: zrw     ! grid scale total water (kg/kg) 
    real(kind=realkind), dimension(:,:), allocatable  :: zrv     ! grid scale water vapor (kg/kg) 
    real(kind=realkind), dimension(:,:), allocatable  :: zrc     ! grid scale cloud water (kg/kg) 
    real(kind=realkind), dimension(:,:), allocatable  :: zri     ! grid scale cloud ice (kg/kg) 
    real(kind=realkind), dimension(:),   allocatable  :: zdxdy   ! grid area (m^2)
    !cgj
    real(kind=realkind), dimension(:),   allocatable  :: zscdep   ! shallow conv cloud depth (m) 
    !
    ! updraft variables
    real(kind=realkind), dimension(:,:), allocatable  :: zumf    ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zuer    ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zudr    ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zuthl   ! updraft enthalpy (j/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zuthv   ! updraft theta_v (k)
    real(kind=realkind), dimension(:,:), allocatable  :: zurw    ! updraft total water (kg/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zurc    ! updraft cloud water (kg/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zuri    ! updraft cloud ice   (kg/kg)
    real(kind=realkind), dimension(:),   allocatable  :: zmflcl  ! cloud base unit mass flux(kg/s) 
    real(kind=realkind), dimension(:),   allocatable  :: zcape   ! available potent. energy     
    real(kind=realkind), dimension(:),   allocatable  :: zthlcl  ! updraft theta at lcl
    real(kind=realkind), dimension(:),   allocatable  :: ztlcl   ! updraft temp. at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zrvlcl  ! updraft rv at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zwlcl   ! updraft w at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zzlcl   ! lcl height
    real(kind=realkind), dimension(:),   allocatable  :: zthvelcl! envir. theta_v at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zwstar  ! jyj wstar free convection velocity
    !
    ! downdraft variables
    real(kind=realkind), dimension(:,:), allocatable  :: zdmf    ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zder    ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zddr    ! downdraft detrainment (kg/s)
    !
    ! closure variables
    real(kind=realkind), dimension(:,:), allocatable  :: zlmass  ! mass of model layer (kg)
    real(kind=realkind), dimension(:),   allocatable  :: ztimec  ! advective time period
    !
    real(kind=realkind), dimension(:,:), allocatable  :: zthc    ! conv. adj. grid scale theta
    real(kind=realkind), dimension(:,:), allocatable  :: zrvc    ! conv. adj. grid scale r_w 
    real(kind=realkind), dimension(:,:), allocatable  :: zrcc    ! conv. adj. grid scale r_c 
    real(kind=realkind), dimension(:,:), allocatable  :: zric    ! conv. adj. grid scale r_i 
    real(kind=realkind), dimension(:,:), allocatable  :: zwsub   ! envir. compensating subsidence (pa/s)
    !
    logical, dimension(:),allocatable  :: gtrig1  ! logical mask for convection    
    logical, dimension(:),allocatable  :: gwork   ! logical work array
    integer, dimension(:),allocatable  :: iindex, ijindex, ijsindex, ijpindex!hor.index
    real(kind=realkind), dimension(:),   allocatable  :: zcph    ! specific heat c_ph 
    real(kind=realkind), dimension(:),   allocatable  :: zlv, zls! latent heat of vaporis., sublim.
    real(kind=realkind)                               :: zes     ! saturation vapor mixng ratio
    real(kind=realkind)                               :: zw1     ! work variable
    !
    ! chemical tracers:
    real(kind=realkind), dimension(:,:,:), allocatable:: zch1    ! grid scale chemical specy (kg/kg)
    real(kind=realkind), dimension(:,:,:), allocatable:: zch1c   ! conv. adjust. chemical specy 1
    real(kind=realkind), dimension(:,:),   allocatable:: zwork3  ! conv. adjust. chemical specy 1
    logical, dimension(:,:,:),allocatable::gtrig4 ! logical mask
    !
    !-------------------------------------------------------------------------------
    !
    !
    !       0.3    compute loop bounds
    !               -------------------
    !
    iib = kidia
    iie = kfdia
    jcvexb = max( 0, kbdia - 1 )
    ikb = 1 + jcvexb 
    iks = klev
    jcvext = max( 0, ktdia - 1)
    ike = iks - jcvext 
    !
    !
    !       0.5    update convective counter ( where kcount > 0 
    !               convection is still active ).
    !               ---------------------------------------------
    !
    gtrig(:) = .false.
    gtrig(iib:iie) = .true.
    !jyj turn off the shallow if deep has been triggered---------------
    do ji = iib, iie
       if (indexcv(ji) == 2 ) gtrig(ji) = .false.
    enddo
    !jyj turn off the shallow if deep has been triggered---------------
    itest = count( gtrig(:) )
    if ( itest == 0 ) return 

    !
    !
    !       0.7    reset convective tendencies to zero if convective
    !               counter becomes negative
    !               -------------------------------------------------
    !
    gtrig3(:,:) = spread( gtrig(:), dim=2, ncopies=iks )
    where ( gtrig3(:,:) ) 
       ptten(:,:)  = 0._realkind
       prvten(:,:) = 0._realkind
       prcten(:,:) = 0._realkind
       priten(:,:) = 0._realkind
       !   puten(:,:)  = 0._realkind
       !   pvten(:,:)  = 0._realkind
       pumf(:,:)   = 0._realkind
    end where
    where ( gtrig(:) ) 
       kcltop(:)  = 0
       kclbas(:)  = 0
       indexcv(:) = 0  !jyj
       !cgj
       pscdep(:)  = 0._realkind
    end where
    if ( och1conv ) then
       allocate( gtrig4(klon,klev,kch1) )
       gtrig4(:,:,:) = spread( gtrig3(:,:), dim=3, ncopies=kch1 )
       where( gtrig4(:,:,:) ) pch1ten(:,:,:) = 0._realkind
       deallocate( gtrig4 )
    end if
    !
    !
    !       1.     initialize  local variables
    !               ----------------------------
    !
    zeps   = xrd / xrv
    zepsa  = xrv / xrd 
    zepsb  = xcpv / xcpd - zepsa
    zcpord = xcpd / xrd
    zrdocp = xrd / xcpd
    !
    !
    !       1.1    set up grid scale theta, theta_v, theta_es 
    !               ------------------------------------------
    !
    ztht(:,:) = 300._realkind
    zsthv(:,:)= 300._realkind
    zsthes(:,:) = 400._realkind
    do jk = ikb, ike
       do ji = iib, iie
          if ( ppabst(ji,jk) > 40.e2_realkind ) then
             ztht(ji,jk)  = ptt(ji,jk) * ( xp00 / ppabst(ji,jk) ) ** zrdocp
             zsthv(ji,jk) = ztht(ji,jk) * ( 1._realkind + zepsa * prvt(ji,jk) ) /              &
                  ( 1._realkind + prvt(ji,jk) + prct(ji,jk) + prit(ji,jk) )
             !
             ! use conservative bolton (1980) formula for theta_e
             ! it is used to compute cape for undilute parcel ascent
             ! for economical reasons we do not use routine convect_satmixratio here
             !
             zes = exp( xalpw - xbetaw / ptt(ji,jk) - xgamw * log( ptt(ji,jk) ) )
             zes = min( 1._realkind, zeps * zes / ( ppabst(ji,jk) - zes ) )
             zsthes(ji,jk) = ptt(ji,jk) * ( ztht(ji,jk) / ptt(ji,jk) ) **             &
                  ( 1._realkind - 0.28_realkind * zes ) * exp( ( 3374.6525_realkind / ptt(ji,jk) - 2.5403_realkind ) &
                  * zes * ( 1._realkind + 0.81_realkind * zes ) )
          end if
       end do
    end do
    !
    !
    !
    !       2.     test for convective columns and determine properties at the lcl
    !               --------------------------------------------------------------
    !
    !       2.1    allocate arrays depending on number of model columns that need
    !               to be tested for convection (i.e. where no convection is present
    !               at the moment.
    !               --------------------------------------------------------------
    !
    allocate( zpres(itest,iks) )
    allocate( zz(itest,iks) )
    allocate( zw(itest,iks) )
    allocate( zth(itest,iks) )
    allocate( zthv(itest,iks) )
    allocate( zthest(itest,iks) )
    allocate( zrv(itest,iks) )
    allocate( zsthlcl(itest) )
    allocate( zstlcl(itest) )
    allocate( zsrvlcl(itest) )
    allocate( zswlcl(itest) )
    allocate( zszlcl(itest) )
    allocate( zsthvelcl(itest) )
    allocate( isdpl(itest) )
    allocate( ispbl(itest) )
    allocate( islcl(itest) )
    allocate( zsdxdy(itest) )
    allocate( zswstar(itest) )   !jyj wstar
    allocate( zsland(itest) )    !jyj
    allocate( gtrig1(itest) )
    allocate( iindex(klon) )
    allocate( ijsindex(itest) )
    !cgj
    allocate( zsscdep(itest) )
    !cgj
    do ji = 1, klon
       iindex(ji) = ji
    end do
    ijsindex(:) = pack( iindex(:), mask=gtrig(:) )
    !
    do jk = ikb, ike
       do ji = 1, itest
          jl = ijsindex(ji)
          zpres(ji,jk)  = ppabst(jl,jk)
          zz(ji,jk)     = pzz(jl,jk)
          zth(ji,jk)    = ztht(jl,jk)
          zthv(ji,jk)   = zsthv(jl,jk)
          zthest(ji,jk) = zsthes(jl,jk)
          zrv(ji,jk)    = max( 0._realkind, prvt(jl,jk) )
          zw(ji,jk)     = pwt(jl,jk)
       end do
    end do
    do ji = 1, itest
       jl = ijsindex(ji)
       zsdxdy(ji)    = xa25
       zsland (ji)   = gcrow(jl)   !jyj
       zswstar(ji)   = wstar(jl)   !jyj wstar
    end do
    !
    !       2.2    compute environm. enthalpy and total water = r_v + r_i + r_c 
    !               and envir. saturation theta_e
    !               ------------------------------------------------------------
    !
    !
    !       2.3    test for convective columns and determine properties at the lcl
    !               --------------------------------------------------------------
    !
    islcl(:) = max( ikb, 2 )   ! initialize dpl pbl and lcl 
    isdpl(:) = ikb
    ispbl(:) = ikb
    !
    call convect_trigger_shal(  itest, klev,                              &
         zpres, zth, zthv, zthest,                 &
         zrv, zw, zz, zsdxdy,zsland,               &   !jyj
         zsthlcl, zstlcl, zsrvlcl, zswlcl, zszlcl, &
         zsthvelcl, zsscdep, islcl, isdpl, ispbl, gtrig1    )
    !
    deallocate( zpres )
    deallocate( zz )
    deallocate( zth )
    deallocate( zthv )
    deallocate( zthest )
    deallocate( zrv )
    deallocate( zw )
    !
    !
    !       3.     after the call of trigger_funct we allocate all the dynamic
    !               arrays used in the convection scheme using the mask gtrig, i.e.
    !               we do calculus only in convective columns. this corresponds to
    !               a gather operation.
    !               --------------------------------------------------------------
    !
    iconv = count( gtrig1(:) )
    if ( iconv == 0 )  then 
       deallocate( zsthlcl )
       deallocate( zstlcl )
       deallocate( zsrvlcl )
       deallocate( zswlcl )
       deallocate( zszlcl )
       deallocate( zsthvelcl )
       deallocate( zsdxdy )
       deallocate( islcl )
       deallocate( isdpl )
       deallocate( ispbl )
       deallocate( gtrig1 )
       deallocate( iindex )
       deallocate( ijsindex )
       deallocate( zswstar )     !jyj wstar
       !cgj
       deallocate( zsscdep )
       return   ! no convective column has been found, exit deep_convection
    endif
    !
    ! vertical index variables
    !
    allocate( idpl(iconv) )
    allocate( ipbl(iconv) )
    allocate( ilcl(iconv) )
    allocate( ictl(iconv) )
    allocate( ietl(iconv) )
    !
    ! grid scale variables
    !
    allocate( zz(iconv,iks) )
    allocate( zpres(iconv,iks) )
    allocate( zdpres(iconv,iks) )
    allocate( ztt(iconv, iks) )
    allocate( zth(iconv,iks) )
    allocate( zthv(iconv,iks) )
    allocate( zthl(iconv,iks) )
    allocate( zthes(iconv,iks) )
    allocate( zrv(iconv,iks) )
    allocate( zrc(iconv,iks) )
    allocate( zri(iconv,iks) )
    allocate( zrw(iconv,iks) )
    allocate( zdxdy(iconv) )
    allocate( zwstar(iconv) )   !jyj wstar
    !cgj
    allocate( zscdep(iconv) )
    !
    ! updraft variables
    !
    allocate( zumf(iconv,iks) )
    allocate( zuer(iconv,iks) )
    allocate( zudr(iconv,iks) )
    allocate( zuthl(iconv,iks) )
    allocate( zuthv(iconv,iks) )
    allocate( zurw(iconv,iks) )
    allocate( zurc(iconv,iks) )
    allocate( zuri(iconv,iks) )
    allocate( zthlcl(iconv) )
    allocate( ztlcl(iconv) )
    allocate( zrvlcl(iconv) )
    allocate( zwlcl(iconv) )
    allocate( zmflcl(iconv) )
    allocate( zzlcl(iconv) )
    allocate( zthvelcl(iconv) )
    allocate( zcape(iconv) )
    !
    ! work variables
    !
    allocate( ijindex(iconv) )
    allocate( ijpindex(iconv) )
    allocate( zcph(iconv) )
    allocate( zlv(iconv) )
    allocate( zls(iconv) )
    !
    !
    !           3.1    gather grid scale and updraft base variables in
    !                   arrays using mask gtrig
    !                   ---------------------------------------------------
    !
    gtrig(:)      = unpack( gtrig1(:), mask=gtrig(:), field=.false. )  
    ijindex(:)    = pack( iindex(:), mask=gtrig(:) )
    !
    do jk = ikb, ike
       do ji = 1, iconv
          jl = ijindex(ji)
          zz(ji,jk)     = pzz(jl,jk)
          zpres(ji,jk)  = ppabst(jl,jk)
          ztt(ji,jk)    = ptt(jl,jk)
          zth(ji,jk)    = ztht(jl,jk)
          zthes(ji,jk)  = zsthes(jl,jk)
          zrv(ji,jk)    = max( 0._realkind, prvt(jl,jk) )
          zrc(ji,jk)    = max( 0._realkind, prct(jl,jk) )
          zri(ji,jk)    = max( 0._realkind, prit(jl,jk) )
          zthv(ji,jk)   = zsthv(jl,jk)
       end do
    end do
    !
    do ji = 1, itest
       ijsindex(ji) = ji	
    end do
    ijpindex(:) = pack( ijsindex(:), mask=gtrig1(:) )
    do ji = 1, iconv
       jl = ijpindex(ji)
       idpl(ji)      = isdpl(jl)
       ipbl(ji)      = ispbl(jl)
       ilcl(ji)      = islcl(jl)
       zthlcl(ji)    = zsthlcl(jl)
       ztlcl(ji)     = zstlcl(jl)
       zrvlcl(ji)    = zsrvlcl(jl)
       zwlcl(ji)     = zswlcl(jl)
       zzlcl(ji)     = zszlcl(jl)
       zthvelcl(ji)  = zsthvelcl(jl)
       zdxdy(ji)     = zsdxdy(jl)
       zwstar(ji)    = zswstar(jl)      !jyj wstar
       !cgj
       zscdep(ji)    = zsscdep(jl)
    end do
    allocate( gwork(iconv) )
    gwork(:)      = pack( gtrig1(:),  mask=gtrig1(:) ) 
    deallocate( gtrig1 )
    allocate( gtrig1(iconv) )
    gtrig1(:)     = gwork(:)
    !                 
    deallocate( gwork )
    deallocate( ijpindex )
    deallocate( isdpl )
    deallocate( ispbl )
    deallocate( islcl )
    deallocate( zsthlcl )
    deallocate( zstlcl )
    deallocate( zsrvlcl )
    deallocate( zswlcl )
    deallocate( zszlcl )
    deallocate( zsthvelcl )
    deallocate( zsdxdy )
    deallocate( zswstar )     !jyj wstar
    !cgj
    deallocate( zsscdep )
    !
    !
    !           3.2    compute pressure difference 
    !                   ---------------------------------------------------
    !
    zdpres(:,ikb) = 0._realkind
    do jk = ikb + 1, ike
       zdpres(:,jk)  = zpres(:,jk-1) - zpres(:,jk)
    end do
    !
    !           3.3   compute environm. enthalpy and total water = r_v + r_i + r_c
    !                  ----------------------------------------------------------
    !
    do jk = ikb, ike, 1
       zrw(:,jk)  = zrv(:,jk) + zrc(:,jk) + zri(:,jk)
       zcph(:)    = xcpd + xcpv * zrw(:,jk)
       zlv(:)     = xlvtt + ( xcpv - xcl ) * ( ztt(:,jk) - xtt ) ! compute l_v
       zls(:)     = xlstt + ( xcpv - xci ) * ( ztt(:,jk) - xtt ) ! compute l_i
       zthl(:,jk) = zcph(:) * ztt(:,jk) + ( 1._realkind + zrw(:,jk) ) * xg * zz(:,jk) &
            - zlv(:) * zrc(:,jk) - zls(:) * zri(:,jk)
    end do
    !
    deallocate( zcph )
    deallocate( zlv )
    deallocate( zls )
    !
    !
    !           4.     compute updraft properties 
    !                   ----------------------------
    !
    !           4.1    set mass flux at lcl ( here a unit mass flux with w = 1 m/s ) 
    !                   -------------------------------------------------------------
    !
    zdxdy(:)  = xa25
    zmflcl(:) = xa25 * 1.e-3_realkind
    zmflcl(:) = xa25 * 0.03_realkind * zwstar(:)  !jyj wstar, 0.03 is the value in grant (2001)
    !
    !
    call convect_updraft_shal( iconv, klev,                                     &
         kice, zpres, zdpres, zz, zthl, zthv, zthes, zrw, &
         zthlcl, ztlcl, zrvlcl, zwlcl, zzlcl, zthvelcl,   & 
         zmflcl, gtrig1, ilcl, idpl, ipbl,                &
         zumf, zuer, zudr, zuthl, zuthv, zurw,            &
         zurc, zuri, zcape, ictl, ietl                    )
    !
    !
    !
    !           4.2    in routine updraft gtrig1 has been set to false when cloud 
    !                   thickness is smaller than 3 km
    !                   -----------------------------------------------------------
    !
    !
    iconv1 = count(gtrig1) 
    !
    if ( iconv1 > 0 )  then
       !
       !       4.3    allocate memory for downdraft variables
       !               ---------------------------------------
       !
       ! downdraft variables
       !
       allocate( zdmf(iconv,iks) )
       allocate( zder(iconv,iks) )
       allocate( zddr(iconv,iks) )
       allocate( ilfs(iconv) )
       allocate( zlmass(iconv,iks) )
       zdmf(:,:) = 0._realkind
       zder(:,:) = 0._realkind
       zddr(:,:) = 0._realkind
       ilfs(:)   = ikb
       do jk = ikb, ike
          zlmass(:,jk)  = zdxdy(:) * zdpres(:,jk) / xg  ! mass of model layer
       end do
       zlmass(:,ikb) = zlmass(:,ikb+1)
       !
       ! closure variables
       !
       allocate( ztimec(iconv) )
       allocate( zthc(iconv,iks) )
       allocate( zrvc(iconv,iks) )
       allocate( zrcc(iconv,iks) )
       allocate( zric(iconv,iks) )
       allocate( zwsub(iconv,iks) )
       !
       !
       !           5.     compute downdraft properties 
       !                   ----------------------------
       !
       ztimec(:) = xctime_shal
       if ( osettadj ) ztimec(:) = ptadjs
       !
       !           7.     determine adjusted environmental values assuming
       !                   that all available buoyant energy must be removed
       !                   within an advective time step ztimec.
       !                   ---------------------------------------------------
       !
       call convect_closure_shal( iconv, klev,                         &
            zpres, zdpres, zz, zdxdy, zlmass,    &
            zthl, zth, zrw, zrc, zri, gtrig1,    &
            zthc, zrvc, zrcc, zric, zwsub,       &
            ilcl, idpl, ipbl, ictl,              &
            zumf, zuer, zudr, zuthl, zurw,       &
            zurc, zuri, zcape, ztimec, iftsteps  )

       !
       !
       !
       !           8.     determine the final grid-scale (environmental) convective 
       !                   tendencies and set convective counter
       !                   --------------------------------------------------------
       !
       !
       !           8.1    grid scale tendencies
       !                   ---------------------
       !
       ! in order to save memory, the tendencies are temporarily stored
       ! in the tables for the adjusted grid-scale values
       !
       do jk = ikb, ike
          zthc(:,jk) = ( zthc(:,jk) - zth(:,jk) ) / ztimec(:)             &
               * ( zpres(:,jk) / xp00 ) ** zrdocp ! change theta in temperature
          zrvc(:,jk) = ( zrvc(:,jk) - zrw(:,jk) + zrc(:,jk) + zri(:,jk) ) &
               / ztimec(:)

          zrcc(:,jk) = ( zrcc(:,jk) - zrc(:,jk) ) / ztimec(:)
          zric(:,jk) = ( zric(:,jk) - zri(:,jk) ) / ztimec(:)
          !
       end do
       !
       !
       !           8.2    apply conservation correction
       !                   -----------------------------
       !
       ! compute vertical integrals
       !
       jkm = maxval( ictl(:) )
       zwork2(:) = 0._realkind
       zwork2b(:) = 0._realkind
       do jk = jkm, ikb+1, -1
          jkp = jk + 1
          do ji = 1, iconv
             zw1 =  zrvc(ji,jk) + zrcc(ji,jk) + zric(ji,jk)
             zwork2(ji) = zwork2(ji) +  zw1 *          & ! moisture
                  (zpres(ji,jk) - zpres(ji,jkp)) / xg
             zw1 = ( xcpd + xcpv * zrw(ji,jk) )* zthc(ji,jk)   - &
                  ( xlvtt + ( xcpv - xcl ) * ( ztt(ji,jk) - xtt ) ) * zrcc(ji,jk) - &
                  ( xlstt + ( xcpv - xcl ) * ( ztt(ji,jk) - xtt ) ) * zric(ji,jk) 
             zwork2b(ji) = zwork2b(ji) + zw1 *         & ! energy
                  (zpres(ji,jk) - zpres(ji,jkp)) / xg
          end do
       end do
       !
       ! budget error (integral must be zero)
       !
       do ji = 1, iconv
          if ( ictl(ji) > 2 ) then
             jkp = ictl(ji) + 1
             zwork2(ji) =  zwork2(ji) * xg / ( zpres(ji,ikb+1) - zpres(ji,jkp) )
             zwork2b(ji) = zwork2b(ji) * xg / ( zpres(ji,ikb+1) - zpres(ji,jkp) )
          end if
       end do
       !
       ! apply uniform correction
       !
       do jk = jkm, ikb+1, -1
          do ji = 1, iconv
             if ( ictl(ji) > 2 .and. jk <= ictl(ji) ) then
                zrvc(ji,jk) = zrvc(ji,jk) - zwork2(ji)                                ! moisture
                zthc(ji,jk) = zthc(ji,jk) - zwork2b(ji) / ( xcpd + xcpv * zrw(ji,jk) )! energy
             end if
          end do
       end do
       !
       !
       !                   cloud base and top levels
       !                   -------------------------
       !
       ilcl(:) = min( ilcl(:), ictl(:) )
       do ji = 1, iconv
          jl = ijindex(ji)
          kcltop(jl) = ictl(ji)
          kclbas(jl) = ilcl(ji)
          indexcv(jl)= 1  !jyj
          !cgj
          pscdep(jl) = zscdep(ji)
       end do
       !
       ! execute a "scatter"= pack command to store the tendencies in
       ! the final 2d tables
       !
       do jk = ikb, ike
          do ji = 1, iconv
             jl = ijindex(ji)
             ptten(jl,jk)   = zthc(ji,jk)
             prvten(jl,jk)  = zrvc(ji,jk)
             prcten(jl,jk)  = zrcc(ji,jk)
             priten(jl,jk)  = zric(ji,jk)
          end do
       end do
       !jyj qvap for output--------------------------------------
       do jk = ikb, ike
          do ji = 1, iconv
             jl = ijindex(ji)

             if( jk<idpl(ji) ) then        !below conv source layer
                qvap(jl,jk)=prvt(jl,jk)
                qliq(jl,jk)=0._realkind
                qice(jl,jk)=0._realkind
             elseif( jk>=idpl(ji) .and. jk<kclbas(jl) ) then
                qvap(jl,jk)=zrvlcl(ji)
                qliq(jl,jk)=0._realkind
                qice(jl,jk)=0._realkind
             else
                qvap(jl,jk)=zurw(ji,jk)-zurc(ji,jk)-zuri(ji,jk)        !in the cloud
                qliq(jl,jk)=zurc(ji,jk)
                qice(jl,jk)=zuri(ji,jk)
             endif

          end do
       end do
       !jyj qvap for output--------------------------------------
       !
       !           8.7    compute convective tendencies for tracers
       !                   ------------------------------------------
       !
       if ( och1conv ) then
          !
          allocate( zch1(iconv,iks,kch1) )
          allocate( zch1c(iconv,iks,kch1) )
          allocate( zwork3(iconv,kch1) )
          !
          do jk = ikb, ike
             do ji = 1, iconv
                jl = ijindex(ji)
                zch1(ji,jk,:) = pch1(jl,jk,:)
             end do
          end do
          !
          call convect_chem_transport( iconv, klev, kch1, zch1, zch1c,          &
               idpl, ipbl, ilcl, ictl, ilfs, ilfs,      &
               zumf, zuer, zudr, zdmf, zder, zddr,      &
               ztimec, zdxdy, zdmf(:,1), zlmass, zwsub, &
               iftsteps )
          !
          do jk = ikb, ike
             do jn = 1, kch1
                zch1c(:,jk,jn) = ( zch1c(:,jk,jn)- zch1(:,jk,jn) ) / ztimec(:)
             end do
          end do
          !
          !
          !           8.8    apply conservation correction
          !                   -----------------------------
          !
          ! compute vertical integrals
          !
          jkm = maxval( ictl(:) )
          zwork3(:,:) = 0._realkind
          do jk = jkm, ikb+1, -1
             jkp = jk + 1
             do ji = 1, iconv
                zwork3(ji,:) = zwork3(ji,:) + zch1c(ji,jk,:) *                    &
                     (zpres(ji,jk) - zpres(ji,jkp)) / xg
             end do
          end do
          !
          ! mass error (integral must be zero)
          !
          do ji = 1, iconv
             jkp = ictl(ji) + 1
             zwork3(ji,:) = zwork3(ji,:) *                                     &
                  xg / ( zpres(ji,ikb+1) - zpres(ji,jkp) )
          end do
          !
          ! apply uniform correction but assure positive mass at each level
          !
          do jk = jkm, ikb+1, -1
             do ji = 1, iconv
                if ( jk <= ictl(ji) ) then
                   zch1c(ji,jk,:) = zch1c(ji,jk,:) - zwork3(ji,:)
                   !        zch1c(ji,jk,:) = max( zch1c(ji,jk,:), -zch1(ji,jk,:)/ztimec(ji) )
                end if
             end do
          end do
          !
          do jk = ikb, ike
             do ji = 1, iconv
                jl = ijindex(ji)
                pch1ten(jl,jk,:) = zch1c(ji,jk,:)
             end do
          end do
       end if
       !
       !
       !           9.     write up- and downdraft mass fluxes
       !                   ------------------------------------
       !
       do jk = ikb, ike
          zumf(:,jk)  = zumf(:,jk) / zdxdy(:) ! mass flux per unit area
       end do
       zwork2(:) = 1._realkind
       do jk = ikb, ike
          do ji = 1, iconv
             jl = ijindex(ji)
             pumf(jl,jk) = zumf(ji,jk) * zwork2(jl)
          end do
       end do
       !
       !
       !           10.    deallocate all local arrays
       !                   ---------------------------
       !
       ! downdraft variables
       !
       deallocate( zdmf )
       deallocate( zder )
       deallocate( zddr )
       deallocate( ilfs )
       deallocate( zlmass )
       !
       !   closure variables
       !
       deallocate( ztimec )
       deallocate( zthc )
       deallocate( zrvc )
       deallocate( zrcc )
       deallocate( zric )
       deallocate( zwsub )
       !
       if ( och1conv ) then
          deallocate( zch1 )
          deallocate( zch1c )
          deallocate( zwork3 )
       end if
       !
    endif
    !
    !    vertical index
    !
    deallocate( idpl )
    deallocate( ipbl )
    deallocate( ilcl )
    deallocate( ictl )
    deallocate( ietl )
    !
    ! grid scale variables
    !
    deallocate( zz )
    deallocate( zpres )
    deallocate( zdpres )
    deallocate( ztt )
    deallocate( zth )
    deallocate( zthv )
    deallocate( zthl )
    deallocate( zthes )
    deallocate( zrw )
    deallocate( zrv )
    deallocate( zrc )
    deallocate( zri )
    deallocate( zdxdy )
    !
    ! updraft variables
    !
    deallocate( zumf )
    deallocate( zuer )
    deallocate( zudr )
    deallocate( zuthl )
    deallocate( zuthv )
    deallocate( zurw )
    deallocate( zurc )
    deallocate( zuri )
    deallocate( zthlcl )
    deallocate( ztlcl )
    deallocate( zrvlcl )
    deallocate( zwlcl )
    deallocate( zzlcl )
    deallocate( zthvelcl )
    deallocate( zmflcl )
    deallocate( zcape )
    deallocate( zwstar)     !jyj wstar
    !cgj
    deallocate( zscdep)
    !
    ! work arrays
    !
    deallocate( iindex )
    deallocate( ijindex )
    deallocate( ijsindex )
    deallocate( gtrig1 )
    !
    !
  end subroutine shallow_convection



  subroutine convect_trigger_shal(  klon, klev,                           &
       ppres, pth, pthv, pthes,              &
       prv, pw, pz, pdxdy,pzland,            &  !jyj
       pthlcl, ptlcl, prvlcl, pwlcl, pzlcl,  &
       pthvelcl, pscdep, klcl, kdpl, kpbl, otrig     ) 
    !     ######################################################################
    !
    !! determine convective columns as well as the cloudy values of theta,
    !!     and qv at the lifting condensation level (lcl) 
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine convective columns
    !!   
    !!
    !!
    !!  method
    !!    ------
    !!      computations are done at every model level starting from bottom.
    !!      the use of masks allows to optimise the inner loops (horizontal loops).
    !!      what we look for is the undermost unstable level at each grid point.
    !!      
    !!     
    !!
    !!    external
    !!    --------
    !!     routine convect_satmixratio
    !!     
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_cst
    !!          xg                 ! gravity constant
    !!          xp00               ! reference pressure
    !!          xrd, xrv           ! gaz  constants for dry air and water vapor
    !!          xcpd               ! cpd (dry air)
    !!          xtt                ! triple point temperature
    !!          xbetaw, xgamw      ! constants for vapor saturation pressure
    !!
    !!      module modd_convpar
    !!          xa25               ! reference grid area
    !!          xzlcl              ! maximum height difference between
    !!                             ! the surface and the dpl
    !!          xzpbl              ! minimum mixed layer depth to sustain convection
    !!          xcdepth            ! minimum necessary cloud depth
    !!          xcdepth_d          ! maximum allowed cloud depth
    !!          xdtpert            ! add small temp peturbation
    !!          xnhgam             ! coefficient for buoyancy term in w eq.
    !!                             ! accounting for nh-pressure
    !!
    !!      module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!    reference
    !!    ---------
    !!
    !!      book2 of documentation ( routine trigger_funct)
    !!      fritsch and chappell (1980), j. atm. sci., vol. 37, 1722-1761.
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  20/03/97  select first departure level
    !!                            that produces a cloud thicker than xcdepth
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    integer, intent(in)                   :: klon      ! horizontal loop index
    integer, intent(in)                   :: klev      ! vertical loop index
    real(kind=realkind), dimension(klon),     intent(in) :: pdxdy     ! grid area
    real(kind=realkind), dimension(klon),     intent(in) :: pzland    ! zsland
    real(kind=realkind), dimension(klon,klev),intent(in) :: pth, pthv ! theta, theta_v
    real(kind=realkind), dimension(klon,klev),intent(in) :: pthes     ! envir. satur. theta_e
    real(kind=realkind), dimension(klon,klev),intent(in) :: prv       ! vapor mixing ratio 
    real(kind=realkind), dimension(klon,klev),intent(in) :: ppres     ! pressure
    real(kind=realkind), dimension(klon,klev),intent(in) :: pz        ! height of grid point (m)
    real(kind=realkind), dimension(klon,klev),intent(in) :: pw        ! vertical velocity
    !
    real(kind=realkind), dimension(klon),     intent(out):: pthlcl    ! theta at lcl
    real(kind=realkind), dimension(klon),     intent(out):: ptlcl     ! temp. at lcl
    real(kind=realkind), dimension(klon),     intent(out):: prvlcl    ! vapor mixing ratio at  lcl
    real(kind=realkind), dimension(klon),     intent(out):: pwlcl     ! parcel velocity at  lcl
    real(kind=realkind), dimension(klon),     intent(out):: pzlcl     ! height at lcl (m)
    real(kind=realkind), dimension(klon),     intent(out):: pthvelcl  ! environm. theta_v at lcl (k)
    !cgj
    real(kind=realkind), dimension(klon),     intent(out):: pscdep     ! shallow conv cloud depth (m)
    logical, dimension(klon),  intent(out):: otrig     ! logical mask for convection 
    integer, dimension(klon),  intent(inout):: klcl    ! contains vert. index of lcl
    integer, dimension(klon),  intent(inout):: kdpl    ! contains vert. index of dpl
    integer, dimension(klon),  intent(inout):: kpbl    ! contains index of source layer top
    !
    !       0.2   declarations of local variables :
    !
    integer :: jkk, jk, jkp, jkm, jl, jkt, jt      ! vertical loop index
    integer :: ji                                  ! horizontal loop index 
    integer :: iie, ikb, ike                       ! horizontal + vertical loop bounds
    real(kind=realkind)    :: zeps, zepsa                         ! r_d / r_v, r_v / r_d 
    real(kind=realkind)    :: zcpord, zrdocp                      ! c_pd / r_d, r_d / c_pd
    !
    real(kind=realkind), dimension(klon) :: zthlcl, ztlcl, zrvlcl, & ! locals for pthlcl,ptlcl
         zwlcl,  zzlcl, zthvelcl  ! prvlcl, ....
    integer, dimension(klon) :: idpl, ipbl, ilcl      ! locals for kdpl, ...
    real(kind=realkind), dimension(klon) :: zplcl    ! pressure at lcl
    real(kind=realkind), dimension(klon) :: zzdpl    ! height of dpl 
    real(kind=realkind), dimension(klon) :: zthvlcl  ! theta_v at lcl = mixed layer value
    real(kind=realkind), dimension(klon) :: ztmix    ! mixed layer temperature
    real(kind=realkind), dimension(klon) :: zevmix   ! mixed layer water vapor pressure 
    real(kind=realkind), dimension(klon) :: zdpthmix, zpresmix ! mixed layer depth and pressure
    real(kind=realkind), dimension(klon) :: zcape    ! convective available energy (m^2/s^2/g)
    real(kind=realkind), dimension(klon) :: ztheul   ! updraft equiv. pot. temperature (k)
    real(kind=realkind), dimension(klon) :: zlv, zcph! specific heats of vaporisation, dry air
    real(kind=realkind), dimension(klon) :: zdp      ! pressure between lcl and model layer
    real(kind=realkind), dimension(klon) :: ztop     ! estimated cloud top (m)
    !integer, dimension(klon) :: itop  ! work array to store highest test layer
    real(kind=realkind), dimension(klon) :: zwork1, zwork2, zwork3    ! work arrays
    logical, dimension(klon) :: gtrig, gtrig2          ! local arrays for otrig
    logical, dimension(klon) :: gwork1                 ! work array
    !
    !jyj dtrh--------------------------------------------------
    integer :: k0
    real(kind=realkind), dimension(klev) :: p00,t00,z00,q00,qes
    real(kind=realkind), dimension(klon) :: dtrh
    real(kind=realkind) :: dlp,qenv,qslcl,rhlcl,dqsdt,ees
    real(kind=realkind) :: aliq,bliq,cliq,dliq
    data aliq,bliq,cliq,dliq/613.3_realkind, 17.502_realkind, 4780.8_realkind, 32.19_realkind/
    !jyj dtrh--------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    !
    !       0.3    compute array bounds
    !               --------------------
    !
    iie = klon
    ikb = 1 + jcvexb 
    ike = klev - jcvext 
    !
    !
    !       1.     initialize local variables
    !               --------------------------
    !
    zeps       = xrd / xrv
    zepsa      = xrv / xrd 
    zcpord     = xcpd / xrd
    zrdocp     = xrd / xcpd
    otrig(:)   = .false.
    idpl(:)    = kdpl(:)
    ipbl(:)    = kpbl(:)
    ilcl(:)    = klcl(:)
    !itop(:)    = ikb
    pwlcl(:)   = 0._realkind
    zwlcl(:)   = 0._realkind
    pthlcl(:)  = 1._realkind
    pthvelcl(:)= 1._realkind
    ptlcl(:)   = 1._realkind
    prvlcl(:)  = 0._realkind
    pwlcl(:)   = 0._realkind
    pzlcl(:)   = pz(:,ikb)
    zzdpl(:)   = pz(:,ikb)
    gtrig2(:)  = .true.
    !cgj
    pscdep(:)  = 0._realkind
    !
    !
    !
    !       1.     determine highest necessary loop test layer
    !              -------------------------------------------
    !
    jt = ike - 2
    do jk = ikb + 1, ike - 2
       ! do ji = 1, iie
       !    if ( pz(ji,jk) - pz(ji,ikb) <= xzlcl ) itop(ji) = jk
       ! end do
       if ( pz(1,jk) - pz(1,ikb) < 5.e3_realkind ) jt = jk 
    end do
    !
    !
    !       2.     enter loop for convection test
    !               ------------------------------
    !
    jkp = minval( idpl(:) ) + 1
    !jkt = maxval( itop(:) )
    jkt = jt
    do jkk = jkp, jkt
       !
       gwork1(:) = zzdpl(:) - pz(:,ikb) < xzlcl
       ! we exit the trigger test when the center of the mixed layer is more
       ! than 1500 m  above soil level.
       where ( gwork1(:) )
          zdpthmix(:) = 0._realkind
          zpresmix(:) = 0._realkind
          zthlcl(:)   = 0._realkind
          zrvlcl(:)   = 0._realkind
          zzdpl(:)    = pz(:,jkk)
          idpl(:)     = jkk
       end where
       !
       !
       !       3.     construct a mixed layer of at least 50 hpa (xzpbl)
       !               ------------------------------------------
       !
       do jk = jkk, ike - 1
          jkm = jk + 1
          do ji = 1, iie     
             if ( gwork1(ji) .and. zdpthmix(ji) < xzpbl ) then
                ipbl(ji)     = jk
                zwork1(ji)   = ppres(ji,jk) - ppres(ji,jkm)
                zdpthmix(ji) = zdpthmix(ji) + zwork1(ji)
                zpresmix(ji) = zpresmix(ji) + ppres(ji,jk) * zwork1(ji)
                zthlcl(ji)   = zthlcl(ji)   + pth(ji,jk)   * zwork1(ji)
                zrvlcl(ji)   = zrvlcl(ji)   + prv(ji,jk)   * zwork1(ji)
             end if
          end do
          if ( minval ( zdpthmix(:) ) >= xzpbl ) exit
       end do
       !
       !
       where ( gwork1(:) )
          !
          zpresmix(:) = zpresmix(:) / zdpthmix(:)
          zthlcl(:)   = zthlcl(:)   / zdpthmix(:) + xdtpert ! add small temp perturb.
          zrvlcl(:)   = zrvlcl(:)   / zdpthmix(:)
          zthvlcl(:)  = zthlcl(:) * ( 1._realkind + zepsa * zrvlcl(:) )                 &
               / ( 1._realkind + zrvlcl(:) )
          !
          !       4.1    use an empirical direct solution ( bolton formula )
          !               to determine temperature and pressure at lcl. 
          !               nota: the adiabatic saturation temperature is not
          !                     equal to the dewpoint temperature
          !               ----------------------------------------------------
          !
          ! 
          ztmix(:)  = zthlcl(:) * ( zpresmix(:) / xp00 ) ** zrdocp 
          zevmix(:) = zrvlcl(:) * zpresmix(:) / ( zrvlcl(:) + zeps )
          zevmix(:) = max( 1.e-8_realkind, zevmix(:) )
          zwork1(:) = log( zevmix(:) / 613.3_realkind )
          ! dewpoint temperature
          zwork1(:) = ( 4780.8_realkind - 32.19_realkind * zwork1(:) ) / ( 17.502_realkind - zwork1(:) ) 
          ! adiabatic saturation temperature
          ztlcl(:)  = zwork1(:) - ( .212_realkind + 1.571e-3_realkind * ( zwork1(:) - xtt )      &
               - 4.36e-4_realkind * ( ztmix(:) - xtt ) ) * ( ztmix(:) - zwork1(:) )
          ztlcl(:)  = min( ztlcl(:), ztmix(:) )
          zplcl(:)  = xp00 * ( ztlcl(:) / zthlcl(:) ) ** zcpord
          !
       end where
       !
       !
       !       4.2    correct ztlcl in order to be completely consistent
       !               with mnh saturation formula
       !               ---------------------------------------------
       !
       call convect_satmixratio( klon, zplcl, ztlcl, zwork1, zlv, zwork2, zcph )
       where( gwork1(:) )
          zwork2(:) = zwork1(:) / ztlcl(:) * ( xbetaw / ztlcl(:) - xgamw ) ! dr_sat/dt
          zwork2(:) = ( zwork1(:) - zrvlcl(:) ) /                              &
               ( 1._realkind + zlv(:) / zcph(:) * zwork2(:) ) 
          ztlcl(:)  = ztlcl(:) - zlv(:) / zcph(:) * zwork2(:)
          !
       end where
       !
       !
       !       4.3    if zrvlcl = prvmix is oversaturated set humidity 
       !               and temperature to saturation values. 
       !               ---------------------------------------------
       !
       call convect_satmixratio( klon, zpresmix, ztmix, zwork1, zlv, zwork2, zcph )
       where( gwork1(:) .and. zrvlcl(:) > zwork1(:) )
          zwork2(:) = zwork1(:) / ztmix(:) * ( xbetaw / ztmix(:) - xgamw ) ! dr_sat/dt
          zwork2(:) = ( zwork1(:) - zrvlcl(:) ) /                              &
               ( 1._realkind + zlv(:) / zcph(:) * zwork2(:) ) 
          ztlcl(:)  = ztmix(:) - zlv(:) / zcph(:) * zwork2(:)
          zrvlcl(:) = zrvlcl(:) - zwork2(:)
          zplcl(:)  = zpresmix(:)
          zthlcl(:) = ztlcl(:) * ( xp00 / zplcl(:) ) ** zrdocp
          zthvlcl(:)= zthlcl(:) * ( 1._realkind + zepsa * zrvlcl(:) )                   &
               / ( 1._realkind + zrvlcl(:) )
       end where
       !
       !
       !        5.1   determine  vertical loop index at the lcl and dpl
       !               --------------------------------------------------
       !
       do jk = jkk, ike - 1
          do ji = 1, iie
             if ( zplcl(ji) <= ppres(ji,jk) .and. gwork1(ji) ) ilcl(ji) = jk + 1
          end do
       end do
       !
       !
       !        5.2   estimate height and environm. theta_v at lcl
       !               --------------------------------------------------
       !
       do ji = 1, iie
          jk   = ilcl(ji)
          jkm  = jk - 1
          zdp(ji)    = log( zplcl(ji) / ppres(ji,jkm) ) /                     &
               log( ppres(ji,jk) / ppres(ji,jkm) )
          zwork1(ji) = pthv(ji,jkm) + ( pthv(ji,jk) - pthv(ji,jkm) ) * zdp(ji) 
          ! we compute the precise value of the lcl
          ! the precise height is between the levels ilcl and ilcl-1.
          zwork2(ji) = pz(ji,jkm) + ( pz(ji,jk) - pz(ji,jkm) ) * zdp(ji)
       end do
       where( gwork1(:) )
          zthvelcl(:) = zwork1(:) 
          zzlcl(:)    = zwork2(:)
       end where
       !        
       !jyj dtrh beg-------------------------------------------------------
       do ji  = 1, iie
          jk  = ilcl(ji)
          jkm = jk - 1

          do k0=jkm,jk
             p00(k0) = ppres(ji,k0)
             t00(k0) = pth(ji,k0) * ( p00(k0) / xp00 ) ** zrdocp
             q00(k0) = prv(ji,k0) / (1._realkind + prv(ji,k0))
             z00(k0) = pz (ji,k0)
             ees = aliq*exp((bliq*t00(k0)-cliq)/(t00(k0)-dliq))
             qes(k0) = epsilo*ees/(p00(k0)-ees)
             q00(k0) = min(qes(k0),q00(k0))
             q00(k0) = max(q00(k0),0.1e-8_realkind)
          end do

          !calculate dlp using z instead of log(p)...
          dlp = (zzlcl(ji)-z00(jkm))/(z00(jk)-z00(jkm))

          !estimate specific humidity at lcl...
          qenv = q00(jkm)+(q00(jk)-q00(jkm))*dlp

          qslcl = qes(jkm)+(qes(jk)-qes(jkm))*dlp
          rhlcl = qenv/qslcl
          dqsdt = zrvlcl(ji)*(cliq-bliq*dliq)/((ztlcl(ji)-dliq)*(ztlcl(ji)-dliq))

          dtrh(ji) = 0._realkind
          !pzland in crcm: land=-1 ocean=0 ice=+1
          !       if(pzland(ji)<0.0 .or. (pzland(ji)==0.0.and.zzlcl(ji)>800.) ) then
          !pzland in scm : land=1  ocean=0 
          if(pzland(ji)>0.0_realkind .or. &
               (abs(pzland(ji))<1.e-14_realkind.and.zzlcl(ji)>800._realkind) ) then
             if(rhlcl>=0.70_realkind .and. rhlcl<=0.90_realkind) then
                dtrh(ji) = 0.20_realkind*(rhlcl-0.7_realkind)*zrvlcl(ji)/dqsdt
             elseif(rhlcl>0.9_realkind) then
                dtrh(ji) = (1._realkind/rhlcl-1._realkind)*zrvlcl(ji)/dqsdt
             endif
          endif

       end do
       !jyj dtrh end-------------------------------------------------------

       !
       !       6.     check to see if cloud is bouyant 
       !               --------------------------------
       !
       !      6.1    compute grid scale vertical velocity perturbation term zwork1
       !               -------------------------------------------------------------
       ! 
       !            !  normalize w grid scale to a 25 km refer. grid
       !    do ji = 1, iie
       !       jk  = ilcl(ji)
       !       jkm = jk - 1 
       !       zwork1(ji) =  ( pw(ji,jkm)  + ( pw(ji,jk) - pw(ji,jkm) ) * zdp(ji) )  &
       !                          * sqrt( pdxdy(ji) / xa25 )
       !                         - 0.02 * zzlcl(ji) / xzlcl ! avoid spurious convection
       !    end do
       !            ! compute sign of normalized grid scale w
       !       zwork2(:) = sign( 1., zwork1(:) ) 
       !       zwork1(:) = xwtrig * zwork2(:) * abs( zwork1(:) ) ** 0.333       &
       !                          * ( xp00 / zplcl(:) ) ** zrdocp
       !
       !       6.2    compute parcel vertical velocity at lcl
       !               ---------------------------------------
       !                   
       !    do ji = 1, iie
       !       jkdl = idpl(ji)
       !       zwork3(ji) = xg * zwork1(ji) * ( zzlcl(ji) - pz(ji,jkdl) )       &
       !                      / ( pthv(ji,jkdl) + zthvelcl(ji) )
       !    end do
       !    where( gwork1(:) )
       !      zwlcl(:)  = 1. + .5 * zwork2(:) * sqrt( abs( zwork3(:) ) ) 
       !      gtrig(:)  = zthvlcl(:) - zthvelcl(:) + zwork1(:) > 0. .and.       &
       !                  zwlcl(:) > 0. 
       !    end where
       zwlcl(:) = 1._realkind
       !
       !jyj dtrh-------------------------------------------
       !calculate new theta and temperature at lcl

       zthlcl(:) = zthlcl(:) + ( dtrh(:) + xdtpert ) &
            * ( xp00 / zplcl(:) ) ** zrdocp
       ztlcl (:) = zthlcl(:) * ( zplcl(:) / xp00 ) ** zrdocp
       !jyj dtrh-------------------------------------------

       !
       !       6.3    look for parcel that produces sufficient cloud depth.
       !               the cloud top is estimated as the level where the cape 
       !               is smaller  than a given value (based on vertical velocity eq.)
       !               --------------------------------------------------------------
       !
       ztheul(:) = ztlcl(:) * ( zthlcl(:) / ztlcl(:) )                       &
            ** ( 1._realkind - 0.28_realkind * zrvlcl(:) )  &
            * exp( ( 3374.6525_realkind / ztlcl(:) - 2.5403_realkind ) *       &
            zrvlcl(:) * ( 1._realkind + 0.81_realkind * zrvlcl(:) ) )
       !
       zcape(:) = 0._realkind
       ztop(:)  = 0._realkind
       zwork3(:)= 0._realkind
       jkm = minval( ilcl(:) )
       do jl = jkm, jt
          jk = jl + 1
          do ji = 1, iie
             zwork1(ji) = ( 2._realkind * ztheul(ji) /                                &
                  ( pthes(ji,jk) + pthes(ji,jl) ) - 1._realkind) * ( pz(ji,jk) - pz(ji,jl) )
             if ( jl < ilcl(ji) ) zwork1(ji) = 0._realkind
             zcape(ji)  = zcape(ji) + zwork1(ji)
             zwork2(ji) = xnhgam * xg * zcape(ji) + 1.05_realkind * zwlcl(ji) * zwlcl(ji)
             ! the factor 1.05 takes entrainment into account
             zwork2(ji) = sign( 1._realkind, zwork2(ji) )
             zwork3(ji) = zwork3(ji) + min(0._realkind, zwork2(ji) )
             zwork3(ji) = max( -1._realkind, zwork3(ji) )
             ! nota, the factors zwork2 and zwork3 are only used to avoid
             ! if and goto statements, the difficulty is to extract only
             ! the level where the criterium is first fullfilled
             ztop(ji)   = pz(ji,jl) * .5_realkind * ( 1._realkind+ zwork2(ji) ) * ( 1._realkind + zwork3(ji) ) + &
                  ztop(ji) * .5_realkind * ( 1._realkind - zwork2(ji) )
          end do
       end do
       !
       !
       zwork2(:) = ztop(:) - zzlcl(:)
       where( zwork2(:)  >= xcdepth  .and. zwork2(:) < xcdepth_d .and. gtrig2(:) )
          gtrig2(:)   = .false.
          otrig(:)    = .true.
          ! otrig(:)    = gtrig(:)     ! we  select the first departure level
          pthlcl(:)   = zthlcl(:)    ! that gives sufficient cloud depth
          prvlcl(:)   = zrvlcl(:)
          ptlcl(:)    = ztlcl(:)
          pwlcl(:)    = zwlcl(:)
          pzlcl(:)    = zzlcl(:)
          pthvelcl(:) = zthvelcl(:)
          kdpl(:)     = idpl(:)
          kpbl(:)     = ipbl(:)
          klcl(:)     = ilcl(:)
          !cgj
          pscdep(:)   = zwork2(:)
       end where

       if ( count(.not.otrig(:) ) == 0 ) exit   !jyj   
       !
    end do
    !
    !
  end subroutine convect_trigger_shal


  subroutine convect_updraft_shal( klon, klev,       &
       kice, ppres, pdpres, pz, pthl, pthv, pthes, prw,&
       pthlcl, ptlcl, prvlcl, pwlcl, pzlcl, pthvelcl,  &
       pmflcl, otrig, klcl, kdpl, kpbl,                &
       pumf, puer, pudr, puthl, puthv, purw,           &
       purc, puri, pcape, kctl, ketl )
    !    ###############################################################################
    !
    !! compute updraft properties from dpl to ctl. 
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine updraft properties
    !!      ( mass flux, thermodynamics, precipitation ) 
    !!
    !!
    !!  method
    !!    ------
    !!      computations are done at every model level starting from bottom.
    !!      the use of masks allows to optimise the inner loops (horizontal loops).
    !!      
    !!     
    !!
    !!    external
    !!    --------
    !!     routine convect_mixing_funct
    !!     routine convect_condens
    !!     
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_cst
    !!          xg                 ! gravity constant
    !!          xp00               ! reference pressure
    !!          xrd, xrv           ! gaz  constants for dry air and water vapor
    !!          xcpd, xcpv, xcl    ! cp of dry air, water vapor and liquid water
    !!          xtt                ! triple point temperature
    !!          xlvtt              ! vaporisation heat at xtt
    !!        
    !!
    !!      module modd_convpar_shal
    !!          xa25               ! reference grid area
    !!          xcrad              ! cloud radius
    !!          xcdepth            ! minimum necessary cloud depth
    !!          xentr              ! entrainment constant
    !!          xnhgam             ! coefficient for buoyancy term in w eq.
    !!                             ! accounting for nh-pressure
    !!          xtfrz1             ! begin of freezing interval
    !!          xtfrz2             ! begin of freezing interval
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_updraft)
    !!      kain and fritsch, 1990, j. atmos. sci., vol.
    !!      kain and fritsch, 1993, meteor. monographs, vol.
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  10/12/97
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    integer, intent(in)                    :: klon  ! horizontal dimension
    integer, intent(in)                    :: klev  ! vertical dimension
    integer, intent(in)                    :: kice  ! flag for ice ( 1 = yes,
    !                0 = no ice )
    real(kind=realkind), dimension(klon,klev), intent(in) :: pthl  ! grid scale enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pthv  ! grid scale theta_v     
    real(kind=realkind), dimension(klon,klev), intent(in) :: pthes ! grid scale saturated theta_e 
    real(kind=realkind), dimension(klon,klev), intent(in) :: prw   ! grid scale total water  
    ! mixing ratio 
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres ! pressure (p)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pdpres! pressure difference between 
    ! bottom and top of layer (pa) 
    real(kind=realkind), dimension(klon,klev), intent(in) :: pz    ! height of model layer (m) 
    real(kind=realkind), dimension(klon),     intent(in) :: pthlcl ! theta at lcl
    real(kind=realkind), dimension(klon),     intent(in) :: ptlcl  ! temp. at lcl
    real(kind=realkind), dimension(klon),     intent(in) :: prvlcl ! vapor mixing ratio at  lcl
    real(kind=realkind), dimension(klon),     intent(in) :: pwlcl  ! parcel velocity at lcl (m/s)
    real(kind=realkind), dimension(klon),     intent(in) :: pmflcl ! cloud  base unit mass flux
    ! (kg/s)
    real(kind=realkind), dimension(klon),     intent(in) :: pzlcl  ! height at lcl (m)
    real(kind=realkind), dimension(klon),     intent(in) :: pthvelcl  ! environm. theta_v at lcl (k)
    logical, dimension(klon),  intent(inout):: otrig! logical mask for convection 
    integer, dimension(klon),  intent(in) :: klcl   ! contains vert. index of lcl
    integer, dimension(klon),  intent(in) :: kdpl   ! contains vert. index of dpl 
    integer, dimension(klon),  intent(in) :: kpbl   !  " vert. index of source layertop
    !
    !
    integer, dimension(klon),  intent(out):: kctl   ! contains vert. index of ctl 
    integer, dimension(klon),  intent(out):: ketl   ! contains vert. index of        &
    !equilibrium (zero buoyancy) level 
    real(kind=realkind), dimension(klon,klev), intent(out):: pumf  ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: puer  ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: pudr  ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: puthl ! updraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: puthv ! updraft theta_v (k)
    real(kind=realkind), dimension(klon,klev), intent(out):: purw  ! updraft total water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: purc  ! updraft cloud water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: puri  ! updraft cloud ice   (kg/kg)
    real(kind=realkind), dimension(klon),     intent(out):: pcape  ! available potent. energy
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike  ! horizontal and vertical loop bounds
    integer :: ji             ! horizontal loop index
    integer :: jk, jkp, jkm, jk1, jk2, jkmin   ! vertical loop index
    real(kind=realkind)    :: zepsa, zcvocd  ! r_v / r_d, c_pv / c_pd 
    real(kind=realkind)    :: zcpord, zrdocp ! c_pd / r_d, r_d / c_pd
    !
    real(kind=realkind), dimension(klon)    :: zut             ! updraft temperature (k)
    real(kind=realkind), dimension(klon)    :: zuw1, zuw2      ! square of updraft vert.
    ! velocity at levels k and k+1
    real(kind=realkind), dimension(klon)    :: ze1,ze2,zd1,zd2 ! fractional entrainm./detrain
    ! rates at levels k and k+1
    real(kind=realkind), dimension(klon)    :: zmixf           ! critical mixed fraction  
    real(kind=realkind), dimension(klon)    :: zcph            ! specific heat c_ph 
    real(kind=realkind), dimension(klon)    :: zlv, zls        ! latent heat of vaporis., sublim.       
    real(kind=realkind), dimension(klon)    :: zurv            ! updraft water vapor at level k+1
    real(kind=realkind), dimension(klon)    :: zpi             ! pi=(p0/p)**(rd/cpd)  
    real(kind=realkind), dimension(klon)    :: ztheul          ! theta_e for undilute ascent
    real(kind=realkind), dimension(klon)    :: zwork1, zwork2, zwork3, zwork4, zwork5,   &
         zwork6          ! work arrays
    integer, dimension(klon) :: iwork           ! wok array
    logical, dimension(klon) :: gwork1, gwork2, gwork4, gwork5 
    ! work arrays
    logical, dimension(klon,klev) :: gwork6     ! work array
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        0.3   set loop bounds
    !              ---------------
    !
    ikb = 1 + jcvexb 
    ike = klev - jcvext 
    iie = klon
    !
    !
    !       1.     initialize updraft properties and local variables
    !               -------------------------------------------------
    !
    zepsa      = xrv / xrd 
    zcvocd     = xcpv / xcpd
    zcpord     = xcpd / xrd
    zrdocp     = xrd / xcpd
    !
    pumf(:,:)  = 0._realkind
    puer(:,:)  = 0._realkind
    pudr(:,:)  = 0._realkind
    puthl(:,:) = 0._realkind
    puthv(:,:) = 0._realkind
    purw(:,:)  = 0._realkind
    purc(:,:)  = 0._realkind
    puri(:,:)  = 0._realkind
    zuw1(:)    = pwlcl(:) * pwlcl(:)
    zuw2(:)    = 0._realkind
    ze1(:)     = 0._realkind
    ze1(:)     = 1._realkind   !jyj
    zd1(:)     = 0._realkind
    pcape(:)   = 0._realkind
    kctl(:)    = ikb
    ketl(:)    = klcl(:)
    gwork2(:)  = .true.
    gwork5(:)  = .true.
    zpi(:)     = 1._realkind
    zwork3(:)  = 0._realkind
    zwork4(:)  = 0._realkind
    zwork5(:)  = 0._realkind
    zwork6(:)  = 0._realkind
    gwork1(:)  = .false.
    gwork4(:)  = .false.
    !
    !
    !       1.1    compute undilute updraft theta_e for cape computations
    !               bolton (1980) formula.
    !               define accurate enthalpy for updraft
    !               -----------------------------------------------------
    !
    ztheul(:) = ptlcl(:) * ( pthlcl(:) / ptlcl(:) ) ** ( 1._realkind - 0.28_realkind * prvlcl(:) )  &
         * exp( ( 3374.6525_realkind / ptlcl(:) - 2.5403_realkind ) *                        &
         prvlcl(:) * ( 1._realkind + 0.81_realkind * prvlcl(:) ) )
    !
    !
    zwork1(:) = ( xcpd + prvlcl(:) * xcpv ) * ptlcl(:)                            &
         + ( 1._realkind + prvlcl(:) ) * xg * pzlcl(:)
    !
    !
    !       2.     set updraft properties between dpl and lcl
    !               ------------------------------------------
    !
    jkp = maxval( klcl(:) )
    jkm = minval( kdpl(:) )
    do jk = jkm, jkp
       do ji = 1, iie
          if ( jk >= kdpl(ji) .and. jk < klcl(ji) ) then
             pumf(ji,jk)  = pmflcl(ji)
             puthl(ji,jk) = zwork1(ji) 
             puthv(ji,jk) = pthlcl(ji) * ( 1._realkind + zepsa * prvlcl(ji) ) /             &
                  ( 1._realkind + prvlcl(ji) )
             purw(ji,jk)  = prvlcl(ji) 
          end if
       end do
    end do
    !
    !
    !       3.     enter loop for updraft computations
    !               ------------------------------------
    !
    jkmin = minval( klcl(:) - 1 ) 
    do jk = max( ikb + 1, jkmin ), ike - 1
       zwork6(:) = 1._realkind
       jkp = jk + 1  
       !
       gwork4(:) = jk >= klcl(:) - 1 
       gwork1(:) = gwork4(:) .and. gwork2(:) ! this mask is used to confine
       ! updraft computations between the lcl and the ctl
       !                                                         
       where( jk == klcl(:) - 1 ) zwork6(:) = 0._realkind ! factor that is used in buoyancy
       ! computation at first level above lcl
       !
       !
       !       4.     estimate condensate, l_v l_i, cph and theta_v at level k+1   
       !               ----------------------------------------------------------
       !
       zwork1(:) = purc(:,jk) 
       zwork2(:) = puri(:,jk) 
       call convect_condens( klon, kice, ppres(:,jkp), puthl(:,jk), purw(:,jk),&
            zwork1, zwork2, pz(:,jkp), gwork1, zut, zurv,     &
            purc(:,jkp), puri(:,jkp), zlv, zls, zcph )
       !
       !
       zpi(:) = ( xp00 / ppres(:,jkp) ) ** zrdocp   
       where ( gwork1(:) )
          !
          puthv(:,jkp) = zpi(:) * zut(:) * ( 1._realkind + zepsa * zurv(:) )           &  
               / ( 1._realkind + purw(:,jk) )     
          !
          !
          !       5.     compute square of vertical velocity using entrainment   
          !               at level k
          !               -----------------------------------------------------
          !    
          zwork3(:) = pz(:,jkp) - pz(:,jk) * zwork6(:) -         &
               ( 1._realkind - zwork6(:) ) * pzlcl(:)          ! level thickness  
          zwork4(:) = pthv(:,jk) * zwork6(:) +                   &
               ( 1._realkind - zwork6(:) ) * pthvelcl(:)
          zwork5(:) = 2._realkind * zuw1(:) * puer(:,jk) / max( .1_realkind, pumf(:,jk) )
          zuw2(:)   = zuw1(:) + zwork3(:) * xnhgam * xg *        & 
               ( ( puthv(:,jk) + puthv(:,jkp) ) /       &
               ( zwork4(:) + pthv(:,jkp) ) - 1._realkind )       & ! buoyancy term
               - zwork5(:)                                  ! entrainment term
          !
          !
          !       6.     update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)  
          !               --------------------------------------------------------
          !
          !                    compute level mean vertical velocity  
          zwork2(:)   = 0.5_realkind *                                                    &
               ( sqrt( max( 1.e-2_realkind, zuw2(:) ) ) +                   &
               sqrt( max( 1.e-2_realkind, zuw1(:) ) ) )          
          !
          !
          !       7.     update r_c, r_i, enthalpy, r_w  for precipitation 
          !               -------------------------------------------------------
          !
          purw(:,jkp)  = purw(:,jk) 
          purc(:,jkp)  = purc(:,jkp)
          puri(:,jkp)  = puri(:,jkp)
          puthl(:,jkp) = puthl(:,jk)
          !    
          zuw1(:)      = zuw2(:)       
          !
       end where
       !
       !
       !       8.     compute entrainment and detrainment using conservative
       !               variables adjusted for precipitation ( not for entrainment)
       !               -----------------------------------------------------------
       !
       !       8.1    compute critical mixed fraction by estimating unknown  
       !               t^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
       !               we determine the zero crossing of the linear curve
       !               evaluating the derivative using zmixf=0.1.
       !               -----------------------------------------------------
       !    
       zmixf(:)  = 0.1_realkind   ! starting value for critical mixed fraction
       zwork1(:) = zmixf(:) * pthl(:,jkp)                                     &
            + ( 1._realkind - zmixf(:) ) * puthl(:,jkp) ! mixed enthalpy
       zwork2(:) = zmixf(:) * prw(:,jkp)                                      &
            + ( 1._realkind - zmixf(:) ) * purw(:,jkp)  ! mixed r_w
       !
       call convect_condens( klon, kice, ppres(:,jkp), zwork1, zwork2,        &
            purc(:,jkp), puri(:,jkp), pz(:,jkp), gwork1, zut,&
            zwork3, zwork4, zwork5, zlv, zls, zcph )
       !        put in enthalpy and r_w and get t r_c, r_i (zut, zwork4-5)
       !        
       ! compute theta_v of mixture
       zwork3(:) = zut(:) * zpi(:) * ( 1._realkind + zepsa * (                         &
            zwork2(:) - zwork4(:) - zwork5(:) ) ) / ( 1._realkind + zwork2(:) )
       ! compute final value of critical mixed fraction using theta_v
       ! of mixture, grid-scale and updraft
       zmixf(:) = max( 0._realkind, puthv(:,jkp) - pthv(:,jkp) ) * zmixf(:) /          &
            ( puthv(:,jkp) - zwork3(:) + 1.e-14_realkind )
       zmixf(:) = max( 0._realkind, min( 1._realkind, zmixf(:) ) )

       !    
       !
       !       8.2     compute final midlevel values for entr. and detrainment    
       !                after call of distribution function
       !                -------------------------------------------------------
       !    
       !
       call convect_mixing_funct ( klon, zmixf, 1, ze2, zd2 )
       !       note: routine mixing_funct returns fractional entrainm/detrainm. rates
       !
       ! zwork1(:) = xentr * pmflcl(:) * pdpres(:,jkp) / xcrad ! rate of env. inflow
       !*mod
       zwork1(:) = xentr * xg / xcrad * pumf(:,jk) * ( pz(:,jkp) - pz(:,jk) )
       ! zwork1(:) = xentr * pumf(:,jk) * pdpres(:,jkp) / xcrad ! rate of env. inflow
       !*mod
       zwork2(:) = 0._realkind
       where ( gwork1(:) ) zwork2(:) = 1._realkind
       where ( puthv(:,jkp) > pthv(:,jkp) )
!cgj250811 as in RCA35
          ze2(:)=max(ze2(:),0.5_realkind);zd2(:)=zd2(:)*1.5_realkind        !jyj
          puer(:,jkp) = 0.5_realkind * zwork1(:) * ( ze1(:) + ze2(:) ) * zwork2(:)
          pudr(:,jkp) = 0.5_realkind * zwork1(:) * ( zd1(:) + zd2(:) ) * zwork2(:)
       elsewhere
          puer(:,jkp) = 0._realkind
          pudr(:,jkp) = zwork1(:) * zwork2(:)
          !   pudr(:,jkp) = zwork1(:) * zwork2(:) * 1.5      !jyj
       end where
       !
       !       8.3     determine equilibrium temperature level
       !                --------------------------------------
       !
       where ( puthv(:,jkp) > pthv(:,jkp) .and. jk > klcl(:) + 1 &   
            .and. gwork1(:) )
          ketl(:) = jkp            ! equilibrium temperature level 
       end where
       !
       !       8.4     if the calculated detrained mass flux is greater than    
       !                the total updraft mass flux, or vertical velocity is
       !                negative, all cloud mass detrains at previous model level,
       !                exit updraft calculations - ctl is attained
       !                -------------------------------------------------------
       !
       where( gwork1(:) )                                                   &
            gwork2(:) = pumf(:,jk) - pudr(:,jkp) > 10._realkind .and. zuw2(:) > 0._realkind        
       where ( gwork2(:) ) kctl(:) = jkp   ! cloud top level
       gwork1(:) = gwork2(:) .and. gwork4(:)
       !
       if ( count( gwork2(:) ) == 0 ) exit           
       !
       !
       !       9.   compute cape for undilute ascent using theta_e and 
       !             theta_es instead of theta_v. this estimation produces 
       !             a significantly larger value for cape than the actual one.
       !             ----------------------------------------------------------
       !
       where ( gwork1(:) )
          !
          zwork3(:)   = pz(:,jkp) - pz(:,jk) * zwork6(:) -                      &
               ( 1._realkind - zwork6(:) ) *  pzlcl(:)              ! level thickness
          zwork2(:)   = pthes(:,jk) + ( 1._realkind - zwork6(:) ) *                      &
               ( pthes(:,jkp) - pthes(:,jk) ) / ( pz(:,jkp) - pz(:,jk) ) *          &
               ( pzlcl(:) - pz(:,jk) ) ! linear interpolation for theta_es at lcl
          ! ( this is only done for model level just above lcl
          !
          zwork1(:) = ( 2._realkind * ztheul(:) ) / ( zwork2(:) + pthes(:,jkp) ) - 1._realkind   
          pcape(:)  = pcape(:) + xg * zwork3(:) * max( 0._realkind, zwork1(:) )
          !
          !
          !       10.   compute final values of updraft mass flux, enthalpy, r_w 
          !              at level k+1    
          !              --------------------------------------------------------
          !    
          pumf(:,jkp)  = pumf(:,jk) - pudr(:,jkp) + puer(:,jkp) 
          pumf(:,jkp)  = max( pumf(:,jkp), 0.1_realkind )
          puthl(:,jkp) = ( pumf(:,jk) * puthl(:,jk) +                              &
               puer(:,jkp) * pthl(:,jk) - pudr(:,jkp) * puthl(:,jk) )  &
               / pumf(:,jkp) 
          purw(:,jkp)  = ( pumf(:,jk) * purw(:,jk) +                               &
               puer(:,jkp) * prw(:,jk) - pudr(:,jkp) * purw(:,jk) )    &
               / pumf(:,jkp) 
          !    
          !
          ze1(:) = ze2(:) ! update fractional entrainment/detrainment
          zd1(:) = zd2(:)
          !
       end where
       !
    end do
    !
    !       12.1    set otrig to false if cloud thickness < 0.5km
    !                or > 3km (deep convection) or cape < 1
    !                ------------------------------------------------
    !
    do ji = 1, iie
       jk  = kctl(ji)
       zwork1(ji) = pz(ji,jk) - pzlcl(ji)
       otrig(ji) = zwork1(ji) >= xcdepth  .and. zwork1(ji) < 3.e3_realkind        &
            .and. pcape(ji) > 1._realkind 
    end do
    where( .not. otrig(:) )
       kctl(:) = ikb 
    end where
    ketl(:) = max( ketl(:), klcl(:) + 2 )
    ketl(:) = min( ketl(:), kctl(:) )
    !
    !
    !       12.2    if the etl and ctl are the same detrain updraft mass   
    !                flux at this level
    !                ------------------------------------------------------- 
    !
    zwork1(:) = 0._realkind
    where ( ketl(:) == kctl(:) ) zwork1(:) = 1._realkind
    !
    do ji = 1, iie
       jk = ketl(ji) 
       pudr(ji,jk)   = pudr(ji,jk) +                                    &
            ( pumf(ji,jk) - puer(ji,jk) )  * zwork1(ji)  
       puer(ji,jk)   = puer(ji,jk) * ( 1._realkind - zwork1(ji) )
       pumf(ji,jk)   = pumf(ji,jk) * ( 1._realkind - zwork1(ji) )
       jkp = kctl(ji) + 1
       puer(ji,jkp)  = 0._realkind ! entrainm/detr rates have been already computed
       pudr(ji,jkp)  = 0._realkind ! at level kctl+1, set them to zero
    end do
    !    
    !       12.3    adjust mass flux profiles, detrainment rates, and   
    !                precipitation fallout rates to reflect linear decrease
    !                in mass flux between the etl and ctl
    !                -------------------------------------------------------        
    ! 
    zwork1(:) = 0._realkind
    jk1 = minval( ketl(:) )
    jk2 = maxval( kctl(:) )

    do jk = jk1, jk2
       do ji = 1, iie
          if( jk > ketl(ji) .and. jk <= kctl(ji) ) then
             zwork1(ji) = zwork1(ji) + pdpres(ji,jk)
          end if
       end do
    end do
    !
    do ji = 1, iie
       jk = ketl(ji) 
       zwork1(ji) = pumf(ji,jk) / max( 1._realkind, zwork1(ji) )
    end do
    !
    do jk = jk1 + 1, jk2
       jkp = jk - 1
       do ji = 1, iie
          if ( jk > ketl(ji) .and. jk <= kctl(ji) ) then
             pudr(ji,jk)  = pdpres(ji,jk) * zwork1(ji)
             pumf(ji,jk)  = pumf(ji,jkp) - pudr(ji,jk)
          end if
       end do
    end do
    !
    !         12.4   set mass flux and entrainment in the source layer.
    !                linear increase throughout the source layer.
    !                -------------------------------------------------------
    !
    !iwork(:) = min( kpbl(:), klcl(:) - 1 )
    iwork(:) = kpbl(:)
    do ji = 1, iie
       jk  = kdpl(ji)
       jkp = iwork(ji)
       !          mixed layer depth
       zwork2(ji) = ppres(ji,jk) - ppres(ji,jkp) + pdpres(ji,jk)
    end do
    !
    jkp = maxval( iwork(:) )
    do jk = jkm, jkp
       do ji = 1, iie
          if ( jk >= kdpl(ji)  .and. jk <= iwork(ji) ) then
             puer(ji,jk) = puer(ji,jk) + pmflcl(ji) * pdpres(ji,jk) / ( zwork2(ji) + 0.1_realkind )
             pumf(ji,jk) = pumf(ji,jk-1) + puer(ji,jk)
          end if
       end do
    end do
    !
    !
    !       13.   if cloud thickness is smaller than  .5 km or > 3 km
    !              no shallow convection is allowed
    !              nota: for technical reasons, we stop the convection
    !                    computations in this case and do not go back to
    !                    trigger_funct to look for the next unstable lcl
    !                    which could produce a thicker cloud.
    !              ---------------------------------------------------
    !
    gwork6(:,:) = spread( otrig(:), dim=2, ncopies=klev )
    where ( .not. gwork6(:,:) )
       pumf(:,:)  = 0._realkind
       pudr(:,:)  = 0._realkind
       puer(:,:)  = 0._realkind
       puthl(:,:) = pthl(:,:)
       purw(:,:)  = prw(:,:)
       purc(:,:)  = 0._realkind
       puri(:,:)  = 0._realkind
    end where
    !
  end subroutine convect_updraft_shal


  subroutine convect_closure_adjust_shal( klon, klev, padj, &
       pumf, pzumf, puer, pzuer, pudr, pzudr  )
    !! uses closure adjustment factor to adjust mass flux and to modify
    !!     precipitation efficiency  when necessary. the computations are
    !!     similar to routine convect_precip_adjust.
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to adjust the mass flux using the
    !!      factor padj computed in convect_closure
    !!
    !!
    !!  method
    !!    ------
    !!      computations are done at every model level starting from bottom.
    !!      the use of masks allows to optimise the inner loops (horizontal loops).
    !!      
    !!
    !!    external
    !!    --------
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!     
    !!    none
    !!
    !!    implicit arguments
    !!    ------------------
    !!
    !!    none
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_closure_adjust)
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    26/03/96 
    !!   last modified  15/11/96
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    !
    integer,                    intent(in) :: klon     ! horizontal dimension
    integer,                    intent(in) :: klev     ! vertical dimension
    real(kind=realkind), dimension(klon),      intent(in) :: padj     ! mass adjustment factor
    !
    !
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pumf  ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzumf ! initial value of  "
    real(kind=realkind), dimension(klon,klev), intent(inout) :: puer  ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzuer ! initial value of  "
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pudr  ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzudr ! initial value of  "
    !
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike                 ! horiz. + vert. loop bounds
    integer :: jk                            ! vertical loop index
    !
    !
    !-------------------------------------------------------------------------------
    !
    !       0.3   compute loop bounds
    !              -------------------
    !
    iie  = klon
    ikb  = 1 + jcvexb 
    ike  = klev - jcvext
    !
    !
    !       1.     adjust mass flux by the factor padj to converge to
    !               specified degree of stabilization
    !               ----------------------------------------------------
    !
    do jk = ikb + 1, ike
       pumf(:,jk)  = pzumf(:,jk)   * padj(:)
       puer(:,jk)  = pzuer(:,jk)   * padj(:)
       pudr(:,jk)  = pzudr(:,jk)   * padj(:)
    end do
    !
  end subroutine convect_closure_adjust_shal


  subroutine convect_closure_shal( klon, klev,                                 &
       ppres, pdpres, pz, pdxdy, plmass,           &
       pthl, pth, prw, prc, pri, otrig1,           &
       pthc, prwc, prcc, pric, pwsub,              &
       klcl, kdpl, kpbl, kctl,                     &
       pumf, puer, pudr, puthl, purw,              &
       purc, puri, pcape, ptimec, kftsteps         )
    !    #######################################################################
    !
    !! uses modified fritsch-chappell closure
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the final adjusted 
    !!     (over a time step ptimec) environmental values of theta_l, r_w, r_c, r_i
    !!      the final convective tendencies can then be evaluated in the main
    !!      routine deep_convect by (pthc-pth)/ptimec
    !!
    !!
    !!  method
    !!    ------
    !!      computations are done at every model level starting from bottom.
    !!      the use of masks allows to optimise the inner loops (horizontal loops).
    !!      
    !!     
    !!
    !!    external
    !!    --------
    !!     
    !!    convect_closure_thrvlcl
    !!    convect_closure_adjust_shal
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_cst
    !!          xg                 ! gravity constant
    !!          xp00               ! reference pressure
    !!          xrd, xrv           ! gaz  constants for dry air and water vapor
    !!          xcpd, xcpv         ! specific heat for dry air and water vapor
    !!          xcl, xci           ! specific heat for liquid water and ice
    !!          xtt                ! triple point temperature
    !!          xlvtt, xlstt       ! vaporization, sublimation heat constant
    !!
    !!      module modd_convpar_shal
    !!          xa25               ! reference grid area
    !!          xstabt             ! stability factor in time integration 
    !!          xstabc             ! stability factor in cape adjustment
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_closure)
    !!      fritsch and chappell, 1980, j. atmos. sci.
    !!      kain and fritsch, 1993, meteor. monographs, vol.
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    26/03/96 
    !!   peter bechtold 15/11/96 change for enthalpie, r_c + r_i tendencies
    !!      tony dore   14/10/96 initialise local variables
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    integer,                   intent(in) :: klon   ! horizontal dimension
    integer,                   intent(in) :: klev   ! vertical dimension
    integer, dimension(klon),  intent(in) :: klcl   ! index lifting condens. level
    integer, dimension(klon),  intent(in) :: kctl   ! index for cloud top level
    integer, dimension(klon),  intent(in) :: kdpl   ! index for departure level 
    integer, dimension(klon),  intent(in) :: kpbl   ! index for top of source layer
    real(kind=realkind), dimension(klon),  intent(inout) :: ptimec ! convection time step 
    real(kind=realkind), dimension(klon),     intent(in) :: pdxdy  ! grid area (m^2)
    real(kind=realkind), dimension(klon,klev),intent(in) :: pthl   ! grid scale enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev),intent(in) :: pth    ! grid scale theta        
    real(kind=realkind), dimension(klon,klev),intent(in) :: prw    ! grid scale total water  
    ! mixing ratio 
    real(kind=realkind), dimension(klon,klev),intent(in) :: prc    ! grid scale r_c 
    real(kind=realkind), dimension(klon,klev),intent(in) :: pri    ! grid scale r_i 
    logical, dimension(klon),  intent(in) :: otrig1 ! logical to keep trace of 
    ! convective arrays modified in updraft
    !   
    !
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres  ! pressure (p)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pdpres ! pressure difference between 
    ! bottom and top of layer (pa)
    real(kind=realkind), dimension(klon,klev), intent(in) :: plmass ! mass of model layer (kg)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pz     ! height of model layer (m) 
    real(kind=realkind), dimension(klon),     intent(in)  :: pcape  ! available potent. energy
    integer,                intent(out)   :: kftsteps! maximum of fract time steps
    ! only used for chemical tracers
    !
    !
    real(kind=realkind), dimension(klon,klev), intent(inout):: pumf  ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: puer  ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pudr  ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: puthl  ! updraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: purw   ! updraft total water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: purc   ! updraft cloud water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: puri   ! updraft cloud ice   (kg/kg)
    !
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pthc  ! conv. adj. grid scale theta
    real(kind=realkind), dimension(klon,klev), intent(out)  :: prwc  ! conv. adj. grid scale r_w 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: prcc  ! conv. adj. grid scale r_c 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pric  ! conv. adj. grid scale r_i 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pwsub ! envir. compensating subsidence(pa/s)
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike  ! horizontal + vertical loop bounds
    integer :: iks            ! vertical dimension
    integer :: jk, jkp, jkmax ! vertical loop index
    integer :: ji             ! horizontal loop index
    integer :: jiter          ! iteration loop index
    integer :: jstep          ! fractional time loop index
    real(kind=realkind)    :: zcpord, zrdocp ! c_pd / r_d, r_d / c_pd
    real(kind=realkind)    :: zcvocd, zepsa  ! c_pv / c_pd, r_v / r_d
    !
    real(kind=realkind), dimension(klon,klev) :: zthlc       ! convectively adjusted 
    ! grid scale enthalpy
    real(kind=realkind), dimension(klon,klev) :: zomg        ! conv. environm. subsidence (pa/s)
    real(kind=realkind), dimension(klon,klev) :: zumf        ! non-adjusted updraft mass flux
    real(kind=realkind), dimension(klon,klev) :: zuer        !   "     updraft  entrainm. rate
    real(kind=realkind), dimension(klon,klev) :: zudr        !   "     updraft  detrainm. rate
    real(kind=realkind), dimension(klon)     :: zadj         ! mass adjustment factor
    real(kind=realkind), dimension(klon)     :: zadjmax      ! limit value for zadj
    real(kind=realkind), dimension(klon)     :: zcape        ! new cape after adjustment
    real(kind=realkind), dimension(klon)     :: ztimec       ! fractional convective time step
    real(kind=realkind), dimension(klon,klev):: ztimc        ! 2d work array for ztimec
    !
    real(kind=realkind), dimension(klon)     :: zthlcl       ! new  theta at lcl
    real(kind=realkind), dimension(klon)     :: zrvlcl       ! new  r_v at lcl
    real(kind=realkind), dimension(klon)     :: zzlcl        ! height of lcl
    real(kind=realkind), dimension(klon)     :: ztlcl        ! temperature at lcl
    real(kind=realkind), dimension(klon)     :: ztelcl       ! envir. temper. at lcl
    real(kind=realkind), dimension(klon)     :: ztheul       ! theta_e for undilute ascent
    real(kind=realkind), dimension(klon)     :: zthes1, zthes2! saturation environm. theta_e
    real(kind=realkind), dimension(klon,klev) :: zthmfin, zthmfout, zrwmfin, zrwmfout
    real(kind=realkind), dimension(klon,klev) :: zrcmfin, zrcmfout, zrimfin, zrimfout
    ! work arrays for environm. compensat. mass flux
    real(kind=realkind), dimension(klon)     :: zpi          ! (p/p00)**r_d/c_pd 
    real(kind=realkind), dimension(klon)     :: zlv          ! latent heat of vaporisation
    real(kind=realkind), dimension(klon)     :: zls          ! latent heat of sublimation 
    real(kind=realkind), dimension(klon)     :: zcph         ! specific heat c_ph
    integer, dimension(klon)  :: itstep       ! fractional convective time step
    integer, dimension(klon)  :: icount       ! timestep counter 
    integer, dimension(klon)  :: ilcl         ! index lifting condens. level
    integer, dimension(klon)  :: iwork1       ! work array
    real(kind=realkind), dimension(klon)     :: zwork1, zwork2, zwork3, zwork4, zwork5
    logical, dimension(klon)  :: gwork1, gwork3! work arrays
    logical, dimension(klon,klev) :: gwork4    ! work array
    !
    !
    !-------------------------------------------------------------------------------
    !
    !       0.2    initialize  local variables
    !               ----------------------------
    !
    !
    ztimc(:,:) = 0._realkind
    zthes2(:) = 0._realkind
    zwork1(:) = 0._realkind 
    zwork2(:) = 0._realkind 
    zwork3(:) = 0._realkind 
    zwork4(:) = 0._realkind 
    zwork5(:) = 0._realkind 
    gwork1(:) = .false.
    gwork3(:) = .false.  
    gwork4(:,:) = .false.  
    ilcl(:)   = klcl(:)
    !
    zcpord    = xcpd / xrd
    zrdocp    = xrd / xcpd
    zcvocd    = xcpv / xcpd 
    zepsa     = xrv / xrd
    !
    zadj(:)   = 1._realkind
    zwork5(:) = 1._realkind
    where( .not. otrig1(:) ) zwork5(:) = 0._realkind 
    !
    !
    !       0.3   compute loop bounds
    !              ------------------- 
    !
    iie    = klon
    ikb    = 1 + jcvexb 
    iks    = klev
    ike    = klev - jcvext 
    jkmax  = maxval( kctl(:) )
    !
    !
    !       2.     save initial mass flux values to be used in adjustment procedure
    !               ---------------------------------------------------------------
    !
    zumf(:,:)  = pumf(:,:)
    zuer(:,:)  = puer(:,:)
    zudr(:,:)  = pudr(:,:)
    zomg(:,:)  = 0._realkind
    pwsub(:,:) = 0._realkind 
    !
    !
    !       3.     compute limits on the closure adjustment factor so that the
    !               inflow in convective drafts from a given layer can't be larger 
    !               than the mass contained in this layer initially.
    !               ---------------------------------------------------------------
    !
    zadjmax(:) = 1000._realkind
    iwork1(:) = ilcl(:)
    jkp = minval( kdpl(:) )
    do jk = jkp, ike
       do ji = 1, iie
          if( jk > kdpl(ji) .and. jk <= iwork1(ji) ) then
             zwork1(ji)  = plmass(ji,jk) / ( ( puer(ji,jk) + 1.e-5_realkind ) * ptimec(ji) )
             zadjmax(ji) = min( zadjmax(ji), zwork1(ji) )
          end if
       end do
    end do
    !
    !
    gwork1(:) = otrig1(:)  ! logical array to limit adjustment to not definitively
    ! adjusted columns
    !
    do jk = ikb, ike
       zthlc(:,jk) = pthl(:,jk) ! initialize adjusted envir. values 
       prwc(:,jk)  = prw(:,jk)
       prcc(:,jk)  = prc(:,jk)
       pric(:,jk)  = pri(:,jk)
       pthc(:,jk)  = pth(:,jk)
    end do
    !
    !
    !
    do jiter = 1, 7  ! enter adjustment loop to assure that all cape is
       ! removed within the advective time interval timec
       !
       ztimec(:) = ptimec(:)
       gwork4(:,:)   = spread( gwork1(:), dim=2, ncopies=iks )
       where( gwork4(:,:) ) pwsub(:,:) = 0._realkind
       zomg(:,:)=0._realkind
       !
       do jk = ikb + 1, jkmax 
          jkp = max( ikb + 1, jk - 1 )
          where ( gwork1(:) .and. jk <= kctl(:) )
             !
             !
             !       4.     determine vertical velocity at top and bottom of each layer
             !               to satisfy mass continuity.
             !               ---------------------------------------------------------------
             ! we compute here domega/dp = - g rho dw/dz = 1/dt
             !
             zwork1(:)   = - ( puer(:,jkp) - pudr(:,jkp) ) / plmass(:,jkp)
             !    
             pwsub(:,jk) = pwsub(:,jkp) - pdpres(:,jk-1) * zwork1(:)
             ! we use pdpres(jk-1) and not jkp in order to have zero subsidence
             ! at the first layer
             !
             !   
             !       5.     compute fractional time step. for stability or 
             !               mass conservation reasons one must split full time step ptimec)
             !               ---------------------------------------------------------------
             !
             zwork1(:) = xstabt * pdpres(:,jkp) / ( abs( pwsub(:,jk) ) + 1.e-14_realkind )
             ! the factor xstabt is used for stability reasons
             ztimec(:) = min( ztimec(:), zwork1(:) ) 
             !
             ! transform vertical velocity in mass flux units
             zomg(:,jk) = pwsub(:,jk) * pdxdy(:) / xg 
          end where
       end do
       !
       !
       where( gwork4(:,:) )
          zthlc(:,:) = pthl(:,:) ! reinitialize adjusted envir. values 
          prwc(:,:)  = prw(:,:)  ! when iteration criterium not attained
          prcc(:,:)  = prc(:,:)
          pric(:,:)  = pri(:,:)
          pthc(:,:)  = pth(:,:)
       end where
       !
       ! 
       !        6. check for mass conservation, i.e. zwork1 > 1.e-2
       !           if mass is not conserved, the convective tendencies
       !           automatically become zero.
       !           ----------------------------------------------------
       !
       do ji = 1, iie
          jk=kctl(ji)
          zwork1(ji) = pudr(ji,jk) * pdpres(ji,jk) / ( plmass(ji,jk) + .1_realkind )    &
               - pwsub(ji,jk)
       end do
       where( gwork1(:) .and. abs( zwork1(:) ) - .01_realkind > 0._realkind)
          gwork1(:) = .false.
          ptimec(:) = 1.e-1_realkind
          zwork5(:) = 0._realkind
       end where
       do jk = ikb, ike
          pwsub(:,jk) = pwsub(:,jk) * zwork5(:)
       end do
       gwork4(:,1:ikb) = .false.
       gwork4(:,iks)   = .false.
       !
       itstep(:) = int( ptimec(:) / ztimec(:) ) + 1 
       ztimec(:) = ptimec(:) / real( itstep(:),realkind ) ! adjust  fractional time step
       ! to be an integer multiple of ptimec
       ztimc(:,:)= spread( ztimec(:), dim=2, ncopies=iks )
       icount(:) = 0
       !
       !
       !
       kftsteps = maxval( itstep(:) )
       do jstep = 1, kftsteps ! enter the fractional time step loop here
          !
          icount(:) = icount(:) + 1
          !
          gwork3(:) =  itstep(:) >= icount(:) .and. gwork1(:) 
          !
          !
          !       7.     assign enthalpy and r_w values at the top and bottom of each
          !               layer based on the sign of w
          !               ------------------------------------------------------------
          !
          zthmfin(:,:)   = 0._realkind
          zrwmfin(:,:)   = 0._realkind
          zrcmfin(:,:)   = 0._realkind
          zrimfin(:,:)   = 0._realkind
          zthmfout(:,:)  = 0._realkind
          zrwmfout(:,:)  = 0._realkind
          zrcmfout(:,:)  = 0._realkind
          zrimfout(:,:)  = 0._realkind
          !
          do jk = ikb + 1, jkmax
             gwork4(:,jk) = gwork3(:) .and. jk <= kctl(:)
             jkp = max( ikb + 1, jk - 1 )
             do ji = 1, iie
                if ( gwork3(ji) ) then
                   !
                   zwork1(ji)       = sign( 1._realkind, zomg(ji,jk) )
                   zwork2(ji)       = 0.5_realkind * ( 1._realkind + zwork1(ji) )
                   zwork1(ji)       = 0.5_realkind * ( 1._realkind - zwork1(ji) )
                   zthmfin(ji,jk)   = - zomg(ji,jk) * zthlc(ji,jkp) * zwork1(ji)
                   zthmfout(ji,jk)  =   zomg(ji,jk) * zthlc(ji,jk)  * zwork2(ji)
                   zthmfin(ji,jkp)  = zthmfin(ji,jkp)  + zthmfout(ji,jk) * zwork2(ji)
                   zthmfout(ji,jkp) = zthmfout(ji,jkp) + zthmfin(ji,jk)  * zwork1(ji)
                   zrwmfin(ji,jk)   = - zomg(ji,jk) * prwc(ji,jkp) * zwork1(ji)
                   zrwmfout(ji,jk)  =   zomg(ji,jk) * prwc(ji,jk)  * zwork2(ji)
                   zrwmfin(ji,jkp)  = zrwmfin(ji,jkp)  + zrwmfout(ji,jk) * zwork2(ji)
                   zrwmfout(ji,jkp) = zrwmfout(ji,jkp) + zrwmfin(ji,jk)  * zwork1(ji)
                   zrcmfin(ji,jk)   = - zomg(ji,jk) * prcc(ji,jkp) * zwork1(ji)
                   zrcmfout(ji,jk)  =   zomg(ji,jk) * prcc(ji,jk)  * zwork2(ji)
                   zrcmfin(ji,jkp)  = zrcmfin(ji,jkp)  + zrcmfout(ji,jk) * zwork2(ji)
                   zrcmfout(ji,jkp) = zrcmfout(ji,jkp) + zrcmfin(ji,jk)  * zwork1(ji)
                   zrimfin(ji,jk)   = - zomg(ji,jk) * pric(ji,jkp) * zwork1(ji)
                   zrimfout(ji,jk)  =   zomg(ji,jk) * pric(ji,jk)  * zwork2(ji)
                   zrimfin(ji,jkp)  = zrimfin(ji,jkp)  + zrimfout(ji,jk) * zwork2(ji)
                   zrimfout(ji,jkp) = zrimfout(ji,jkp) + zrimfin(ji,jk)  * zwork1(ji)
                   !
                end if
             end do
          end do
          !
          where ( gwork4(:,:) )
             !
             !
             !
             !       8.     update the environmental values of enthalpy and r_w at each level
             !               nota: these are the main equations of the scheme
             !
             zthlc(:,:) = zthlc(:,:) + ztimc(:,:) / plmass(:,:) * (      &
                  zthmfin(:,:) + pudr(:,:) * puthl(:,:)        &
                  - zthmfout(:,:) - puer(:,:) * pthl(:,:)   )
             prwc(:,:)  = prwc(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrwmfin(:,:) + pudr(:,:) * purw(:,:)          &
                  - zrwmfout(:,:) - puer(:,:) * prw(:,:)    )    
             prcc(:,:)  = prcc(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrcmfin(:,:) + pudr(:,:) * purc(:,:) - zrcmfout(:,:) -  &
                  puer(:,:) * prc(:,:)    )    
             pric(:,:)  = pric(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrimfin(:,:) + pudr(:,:) * puri(:,:) - zrimfout(:,:) -  & 
                  puer(:,:) * pri(:,:)    )    
          end where
          !
       end do ! exit the fractional time step loop
       !
       ! 
       !          10.    compute final linearized value of theta envir.
       !                  ----------------------------------------------
       !
       do jk = ikb + 1, jkmax
          do ji = 1, iie
             if( gwork1(ji) .and. jk <= kctl(ji) ) then
                zpi(ji)    = ( xp00 / ppres(ji,jk) ) ** zrdocp
                zcph(ji)   = xcpd + prwc(ji,jk) * xcpv
                zwork2(ji) = pth(ji,jk) / zpi(ji)  ! first temperature estimate
                zlv(ji)    = xlvtt + ( xcpv - xcl ) * ( zwork2(ji) - xtt )
                zls(ji)    = xlvtt + ( xcpv - xci ) * ( zwork2(ji) - xtt )
                ! final linearized temperature
                zwork2(ji) = ( zthlc(ji,jk) + zlv(ji) * prcc(ji,jk) + zls(ji) * pric(ji,jk) &
                     - (1._realkind + prwc(ji,jk) ) * xg * pz(ji,jk) ) / zcph(ji)
                zwork2(ji) = max( 180._realkind, min( 340._realkind, zwork2(ji) ) )
                pthc(ji,jk)= zwork2(ji) * zpi(ji) ! final adjusted envir. theta
             end if
          end do
       end do
       !
       !
       !         11.     compute new cloud ( properties at new lcl )
       !                     nota: the computations are very close to
       !                           that in routine trigger_funct
       !                  ---------------------------------------------
       !
       call convect_closure_thrvlcl(  klon, klev,                           &
            ppres, pthc, prwc, pz, gwork1,        &
            zthlcl, zrvlcl, zzlcl, ztlcl, ztelcl, &
            ilcl, kdpl, kpbl )
       !
       !
       ztlcl(:)  = max( 230._realkind, min( 335._realkind, ztlcl(:) ) )  ! set some overflow bounds
       ztelcl(:) = max( 230._realkind, min( 335._realkind, ztelcl(:) ) )
       zthlcl(:) = max( 230._realkind, min( 345._realkind, zthlcl(:) ) )
       zrvlcl(:) = max(   0._realkind, min(   1._realkind, zrvlcl(:) ) )
       !
       !
       !         12.    compute adjusted cape
       !                 ---------------------
       !
       zcape(:)  = 0._realkind
       zpi(:)    = zthlcl(:) / ztlcl(:)
       zpi(:)    = max( 0.95_realkind, min( 1.5_realkind, zpi(:) ) )
       zwork1(:) = xp00 / zpi(:) ** zcpord ! pressure at lcl
       !
       call convect_satmixratio( klon, zwork1, ztelcl, zwork3, zlv, zls, zcph )
       zwork3(:) = min(   .1_realkind, max(   0._realkind, zwork3(:) ) )
       !
       ! compute theta_e updraft undilute
       ztheul(:) = ztlcl(:) * zpi(:) ** ( 1._realkind - 0.28_realkind * zrvlcl(:) )            &
            * exp( ( 3374.6525_realkind / ztlcl(:) - 2.5403_realkind )   &
            * zrvlcl(:) * ( 1._realkind + 0.81_realkind * zrvlcl(:) ) )
       !
       ! compute theta_e saturated environment at lcl
       zthes1(:) = ztelcl(:) * zpi(:) ** ( 1._realkind - 0.28_realkind * zwork3(:) )           &
            * exp( ( 3374.6525_realkind / ztelcl(:) - 2.5403_realkind )  &
            * zwork3(:) * ( 1._realkind + 0.81_realkind * zwork3(:) ) )
       !
       do jk = minval( ilcl(:) ), jkmax
          jkp = jk - 1
          do ji = 1, iie
             zwork4(ji) = 1._realkind
             if ( jk == ilcl(ji) ) zwork4(ji) = 0._realkind
             !
             ! compute theta_e saturated environment and adjusted values
             ! of theta
             !
             gwork3(ji)  = jk >= ilcl(ji) .and. jk <= kctl(ji) .and. gwork1(ji) 
             !
             zpi(ji)     = ( xp00 / ppres(ji,jk) ) ** zrdocp
             zwork2(ji)  = pthc(ji,jk) / zpi(ji)
          end do
          !
          call convect_satmixratio( klon, ppres(:,jk), zwork2, zwork3, zlv, zls, zcph )
          !
          !
          do ji = 1, iie
             if ( gwork3(ji) ) then
                zthes2(ji)  = zwork2(ji) * zpi(ji) ** ( 1._realkind - 0.28_realkind * zwork3(ji) )   &
                     * exp( ( 3374.6525_realkind / zwork2(ji) - 2.5403_realkind ) &
                     * zwork3(ji) * ( 1._realkind + 0.81_realkind * zwork3(ji) ) )
                !
                zwork3(ji)  = pz(ji,jk) - pz(ji,jkp) * zwork4(ji) -                &
                     ( 1._realkind - zwork4(ji) ) * zzlcl(ji)    ! level thickness
                zwork1(ji)  = ( 2._realkind * ztheul(ji) ) / ( zthes1(ji) + zthes2(ji) ) - 1._realkind
                zcape(ji)   = zcape(ji) + xg * zwork3(ji) * max( 0._realkind, zwork1(ji) )
                zthes1(ji)  = zthes2(ji)
             end if
          end do
       end do
       !
       !                                                          
       !         13.     determine mass adjustment factor knowing how much
       !                  cape has been removed.
       !                  -------------------------------------------------
       !
       where ( gwork1(:) )
          zwork1(:) = max( pcape(:) - zcape(:), 0.1_realkind * pcape(:) )
          zwork2(:) = zcape(:) / ( pcape(:) + 1.e-8_realkind )
          !       
          gwork1(:) = zwork2(:) > 0.1_realkind .or. abs(zcape(:))< 1.e-14_realkind ! mask for adjustment
       end where
       !
       where ( abs(zcape(:))< 1.e-14_realkind .and. gwork1(:) )  zadj(:) = zadj(:) * 0.5_realkind
       where ( abs(zcape(:)) > 1.e-14_realkind .and. gwork1(:) )                              &
            zadj(:) = zadj(:) * xstabc * pcape(:) / ( zwork1(:) + 1.e-8_realkind )
       zadj(:) = min( zadj(:), zadjmax(:) )  
       !
       !
       !         13.     adjust mass flux by the factor zadj to converge to
       !                  specified degree of stabilization
       !                 ----------------------------------------------------
       !
       call convect_closure_adjust_shal( klon, klev, zadj,                     &
            pumf, zumf, puer, zuer, pudr, zudr    )
       !
       !
       where ( zadj(:) < 0.05_realkind ) gwork1(:) = .false.   ! jyj mask for adjustment
       if ( count( gwork1(:) ) == 0 ) exit ! exit big adjustment iteration loop
       ! when all columns have reached 
       ! desired degree of stabilization.
       !
    end do  ! end of big adjustment iteration loop
    !
    !
    ! skip adj. total water array  to water vapor
    do jk = ikb, ike
       prwc(:,jk) = max( 0._realkind, prwc(:,jk) - prcc(:,jk) - pric(:,jk) )
    end do
    !
    !
  end subroutine convect_closure_shal

end module modd_convpar_shal
