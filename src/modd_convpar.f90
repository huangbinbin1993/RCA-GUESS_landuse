module modd_convpar
  !!- declaration of convection constants 
  !!
  !!    purpose
  !!    -------
  !      the purpose of this declarative module is to declare  the 
  !      constants in the deep convection parameterization.    
  !
  !!
  !!  implicit arguments
  !!    ------------------
  !!      none 
  !!
  !!    reference
  !!    ---------
  !!      book2 of documentation of meso-nh 
  !!          
  !!    author
  !!    ------
  !!      p. bechtold   laboratoire d'aerologie
  !!
  !!    modifications
  !!    -------------
  !!      original    26/03/96                      
  !!   last modified  15/11/96
  !!         updated  18/10/07 by y. jiao (search !jyj for the details)
  !!                  for deep 1) variation cloud radius xcrad
  !!                           2) variation minimum cloud depth xcdepth              
  !!                           3) dilute updraft
  !!               for shallow 1) additional triggering term dtrh
  !!                           2) cloud bass mass flux closure with wstar
  !!                           3) turn off the shallow once deep being triggered, indexcv
  !
  !       0.   declarations
  use modd_cst
  use modd_convparext
  use decomp,only:realkind
  implicit none 
  private
!cgj050711  real(kind=realkind), save :: xa25=525.e6_realkind        ! 25 km x 25 km reference grid area
  real(kind=realkind), save :: xa25=625.e6_realkind        ! 25 km x 25 km reference grid area
  real(kind=realkind), save :: xcrad=1500._realkind       ! cloud radius 
  real(kind=realkind), save :: xcdepth=3.e3_realkind     ! minimum necessary cloud depth
  real(kind=realkind), save :: xentr=0.03_realkind       ! entrainment constant (m/pa) = 0.2 (m)  
  real(kind=realkind), save :: xzlcl=3.5e3_realkind       ! maximum allowed allowed height 
  ! difference between departure level and surface
  !cgj170812real(kind=realkind), save :: xzpbl=60.e2_realkind       ! minimum mixed layer depth to sustain convection
  real(kind=realkind), save :: xzpbl=40.e2_realkind       ! minimum mixed layer depth to sustain convection
  real(kind=realkind), save :: xwtrig=4.64_realkind      ! constant in vertical velocity trigger

  real(kind=realkind), save :: xnhgam=1.3333_realkind      ! accounts for non-hydrost. pressure 
  ! in buoyancy term of w equation
  ! = 2 / (1+gamma)

  real(kind=realkind), save :: xrhdbc=0.9_realkind      ! relative humidity below cloud in downdraft

!cgj250612  real(kind=realkind), save :: xrconv=0.005_realkind      ! constant in precipitation conversion 
  real(kind=realkind), save :: xrconv=0.015_realkind      ! constant in precipitation conversion 
  real(kind=realkind), save :: xstabt=0.75_realkind      ! factor to assure stability in  fractional time
  ! integration, routine convect_closure
  real(kind=realkind), save :: xstabc=0.95_realkind      ! factor to assure stability in cape adjustment,
  !  routine convect_closure
  real(kind=realkind), save :: xusrdpth=165.e2_realkind    ! pressure thickness used to compute updraft
  ! moisture supply rate for downdraft
  real(kind=realkind), save :: xmeldpth=200.e2_realkind    ! layer (pa) through which precipitation melt is
  ! allowed below  melting level
  real(kind=realkind), save :: xuvdp=0.7_realkind       ! constant for pressure perturb in momentum transport

  public deep_convection
contains
  subroutine deep_convection( klon, klev, kidia, kfdia, kbdia, ktdia,        &
       pdtconv, kice, orefresh, odown, osettadj,      &
       ppabst, pzz, pdxdy, ptimec,                    &
       indexcv,                                       &   !jyj
       ptt, prvt, prct, prit, put, pvt, pwt, pworo,poro,     &
       kcount, ctrig ,ptten, prvten, prcten, priten,         &
       pprlten, pprsten,                              &
       kcltop, kclbas, pprlflx, pprsflx,              &
       pumf, pdmf, pudr, pudrf, pcape,                &
       qvap, qliq, qice,                              &   !jyj output cloud plumes variables
       och1conv, kch1, pch1, pch1ten, pdcdep           )
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
    !!    convect_trigger_funct 
    !!    convect_satmixratio
    !!    convect_updraft
    !!        convect_condens
    !!        convect_mixing_funct
    !!    convect_tstep_pref
    !!    convect_downdraft
    !!    convect_precip_adjust
    !!    convect_closure
    !!        convect_closure_thrvlcl
    !!        convect_closure_adjust
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
    !!
    !!         
    !!    reference
    !!    ---------
    !!
    !!      bechtold et al., 2001, quart. j. roy. meteor. soc. : 
    !!           a mass flux convection scheme for regional and global models.
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
    !!   peter bechtold 04/10/97 replace theta_il by enthalpy
    !!         "        10/12/98 changes for arpege
    !!   dominique paquin uqam suivi des corrections de debordements
    !
    !          "          ouranos avril 2003 
    !                     convect_closure correction pcp < 0 
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !

    implicit none
    !       0.1   declarations of dummy arguments :
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
    logical,                    intent(in) :: orefresh ! refresh or not tendencies
    ! at every call
    logical,                    intent(in) :: odown    ! take or not convective
    ! downdrafts into account
    logical,                    intent(in) :: osettadj ! logical to set convective
    ! adjustment time by user 
    real(kind=realkind), dimension(klon,klev), intent(in) :: ptt      ! grid scale temperature at t
    real(kind=realkind), dimension(klon,klev), intent(in) :: prvt     ! grid scale water vapor "
    real(kind=realkind), dimension(klon,klev), intent(in) :: prct     ! grid scale r_c  "
    real(kind=realkind), dimension(klon,klev), intent(in) :: prit     ! grid scale r_i "
    real(kind=realkind), dimension(klon,klev), intent(in) :: put      ! grid scale horiz. wind u "
    real(kind=realkind), dimension(klon,klev), intent(in) :: pvt      ! grid scale horiz. wind v "
    real(kind=realkind), dimension(klon,klev), intent(in) :: pwt      ! grid scale vertical 
    ! velocity (m/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppabst   ! grid scale pressure at t
    real(kind=realkind), dimension(klon,klev), intent(in) :: pzz      ! height of model layer (m) 
    real(kind=realkind), dimension(klon),      intent(in) :: pdxdy    ! horizontal grid area (m-a2)
    !cgj
    real(kind=realkind), dimension(klon),      intent(in) :: pworo    ! orographic standard deviation
    real(kind=realkind), dimension(klon),      intent(in) :: poro     ! fraction of land
    !cgj
    real(kind=realkind), dimension(klon),      intent(in) :: ptimec   ! value of convective adjustment
    ! time if osettadj=.true.
    !   
    integer, dimension(klon),   intent(inout):: kcount ! convective counter (recompute
    !cgj
    integer, dimension(klon),   intent(inout):: ctrig  ! external control to look for convection
    ! at a grid column (1) or not (0)
    real(kind=realkind), dimension(klon,klev), intent(inout):: ptten  ! convective temperature
    ! tendency (k/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: prvten ! convective r_v tendency (1/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: prcten ! convective r_c tendency (1/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: priten ! convective r_i tendency (1/s)

    real(kind=realkind), dimension(klon,klev), intent(out):: qvap     ! cloud plumes total water (kg/kg) !jyj
    real(kind=realkind), dimension(klon,klev), intent(out):: qliq     ! cloud plumes cloud water (kg/kg) !jyj
    real(kind=realkind), dimension(klon,klev), intent(out):: qice     ! cloud plumes cloud ice   (kg/kg) !jyj

    real(kind=realkind), dimension(klon),      intent(inout):: pprlten! liquid surf. precipitation
    ! tendency (m/s)
    real(kind=realkind), dimension(klon),      intent(inout):: pprsten! solid surf. precipitation
    ! tendency (m/s)
    integer, dimension(klon),   intent(inout):: indexcv! convection index (deep=2, shallow=1,none=0)  !jyj
    integer, dimension(klon),   intent(inout):: kcltop ! cloud top level
    integer, dimension(klon),   intent(inout):: kclbas ! cloud base level
    ! they are given a value of
    ! 0 if no convection
    real(kind=realkind), dimension(klon,klev), intent(inout):: pprlflx! liquid precip flux (m/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pprsflx! solid  precip flux (m/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pumf   ! updraft mass flux (kg/s m2)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pudr   ! detrainment rate (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pudrf  ! detrainment rate (/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pdmf   ! downdraft mass flux (kg/s m2)
    real(kind=realkind), dimension(klon),      intent(inout):: pcape  ! maximum cape (j/kg)
    !
    logical,                    intent(in) :: och1conv ! include tracer transport
    integer,                    intent(in) :: kch1     ! number of species
    real(kind=realkind), dimension(klon,klev,kch1), intent(in) :: pch1! grid scale chemical species
    real(kind=realkind), dimension(klon,klev,kch1), intent(inout):: pch1ten! species conv. tendency (1/s)
    !cgj
    real(kind=realkind), dimension(klon),      intent(inout):: pdcdep ! deep conv cloud depth (m)
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
    real(kind=realkind), dimension(klon,klev)         :: ztht, zsthv, zsthes  ! grid scale theta, 
    ! theta_v, theta_es
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
    integer, dimension(:),allocatable  :: idbl    ! index for downdraft base level  
    integer, dimension(:),allocatable  :: iml     ! melting level  
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
    !cgj
    real(kind=realkind), dimension(:), allocatable    :: zworo   ! orographic standard deviation
    real(kind=realkind), dimension(:), allocatable    :: zoro    ! fraction of land 
    !cgj
    real(kind=realkind), dimension(:), allocatable    :: zscrad  ! cloud radius (m) !jyj xcrad
    real(kind=realkind), dimension(:), allocatable    :: zscdep  ! cloud depth  (m) !jyj xcdepth
    !cgj
    real(kind=realkind), dimension(:), allocatable    :: zsdcdep ! deep conv cloud depth (m)
    !
    ! grid scale variables
    real(kind=realkind), dimension(:,:), allocatable  :: zz      ! height of model layer (m) 
    real(kind=realkind), dimension(:,:), allocatable  :: zpres   ! grid scale pressure
    real(kind=realkind), dimension(:,:), allocatable  :: zdpres  ! pressure difference between 
    ! bottom and top of layer (pa)
    real(kind=realkind), dimension(:,:), allocatable  :: zu      ! grid scale horiz. u component on theta grid
    real(kind=realkind), dimension(:,:), allocatable  :: zv      ! grid scale horiz. v component on theta grid
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
    real(kind=realkind), dimension(:),   allocatable  :: zcrad   ! cloud radius (m) !jyj xcrad
    real(kind=realkind), dimension(:),   allocatable  :: zcdep   ! cloud depth  (m) !jyj xcdepth
    !cgj
    real(kind=realkind), dimension(:),   allocatable  :: zdcdep   ! cloud depth  (m) !jyj xcdepth
    !
    ! updraft variables
    real(kind=realkind), dimension(:,:), allocatable  :: zumf    ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zuer    ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zudr    ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zudrf   ! updraft detrainment (/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zupr    ! updraft precipitation in
    ! flux units (kg water / s)
    real(kind=realkind), dimension(:,:), allocatable  :: zuthl   ! updraft enthalpy (j/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zuthv   ! updraft theta_v (k)
    real(kind=realkind), dimension(:,:), allocatable  :: zurw    ! updraft total water (kg/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zurc    ! updraft cloud water (kg/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zuri    ! updraft cloud ice   (kg/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zurr    ! liquid precipit. (kg/kg)
    ! produced in  model layer
    real(kind=realkind), dimension(:,:), allocatable  :: zurs    ! solid precipit. (kg/kg)
    ! produced in  model layer
    real(kind=realkind), dimension(:),   allocatable  :: zutpr   ! total updraft precipitation (kg/s)
    real(kind=realkind), dimension(:),   allocatable  :: zmflcl  ! cloud base unit mass flux(kg/s) 
    real(kind=realkind), dimension(:),   allocatable  :: zcape   ! available potent. energy     
    real(kind=realkind), dimension(:),   allocatable  :: zthlcl  ! updraft theta at lcl
    real(kind=realkind), dimension(:),   allocatable  :: ztlcl   ! updraft temp. at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zrvlcl  ! updraft rv at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zwlcl   ! updraft w at lcl
    real(kind=realkind), dimension(:),   allocatable  :: zzlcl   ! lcl height
    real(kind=realkind), dimension(:),   allocatable  :: zthvelcl! envir. theta_v at lcl
    !
    ! downdraft variables
    real(kind=realkind), dimension(:,:), allocatable  :: zdmf    ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zder    ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zddr    ! downdraft detrainment (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zdthl   ! downdraft enthalpy (j/kg)
    real(kind=realkind), dimension(:,:), allocatable  :: zdrw    ! downdraft total water (kg/kg)
    real(kind=realkind), dimension(:),   allocatable  :: zmixf   ! mixed fraction at lfs        
    real(kind=realkind), dimension(:),   allocatable  :: ztpr    ! total surf precipitation (kg/s)
    real(kind=realkind), dimension(:),   allocatable  :: zspr    ! solid surf precipitation (kg/s)
    real(kind=realkind), dimension(:),   allocatable  :: zdtevr  ! donwndraft evapor. (kg/s)
    real(kind=realkind), dimension(:),   allocatable  :: zpref   ! precipitation efficiency
    real(kind=realkind), dimension(:,:), allocatable  :: zdtevrf ! donwndraft evapor. (kg/s)
    real(kind=realkind), dimension(:,:), allocatable  :: zprlflx ! liquid precip flux
    real(kind=realkind), dimension(:,:), allocatable  :: zprsflx ! solid precip flux
    !
    ! closure variables
    real(kind=realkind), dimension(:,:), allocatable  :: zlmass  ! mass of model layer (kg)
    real(kind=realkind), dimension(:),   allocatable  :: ztimea  ! advective time period
    real(kind=realkind), dimension(:),   allocatable  :: ztimec, ztimed! time during which convection is
    ! active at grid point (as ztime)
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
    !
    ! chemical tracers:
    real(kind=realkind), dimension(:,:,:), allocatable:: zch1    ! grid scale chemical specy (kg/kg)
    real(kind=realkind), dimension(:,:,:), allocatable:: zch1c   ! conv. adjust. chemical specy 1
    real(kind=realkind), dimension(:,:),   allocatable:: zwork3  ! work arrays 
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
    jcvext = max( 0, ktdia - 1 )
    ike = iks - jcvext 
    !
    !
    !       0.5    update convective counter ( where kcount > 0 
    !               convection is still active ).
    !               ---------------------------------------------
    !
    kcount(iib:iie) = kcount(iib:iie) - 1 
    !
    if ( orefresh ) then
       kcount(:) = 1
       kcount(iib:iie) = 0 ! refresh or not at every call
    end if
    !
    !original bkf..gtrig(:)  = kcount(:) <= 0
    gtrig(:)  = kcount(:) <= 0 .and. ctrig(:) > 0
    !cgj calculate the akfcum filter for convective points in bkfcall.f90
    !cgj and pass in a 1/0 (klon) dimensioned integer field, tis should be 
    !cgj then combined with the kcount statement above in order to define
    !cgj gtrig(:)
    itest = count( gtrig(:) )
    if ( itest == 0 ) return  ! if convection is already active at every grid point
    ! exit deep_convection
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
       pprlflx(:,:)= 0._realkind
       pprsflx(:,:)= 0._realkind
       !   puten(:,:)  = 0._realkind
       !   pvten(:,:)  = 0._realkind
       pumf(:,:)   = 0._realkind
       pdmf(:,:)   = 0._realkind
!cgjudr
       pudr(:,:)   = 0._realkind
       pudrf(:,:)   = 0._realkind
    end where
    where ( gtrig(:) ) 
       pprlten(:) = 0._realkind
       pprsten(:) = 0._realkind
       indexcv(:) = 0  !jyj
       kcltop(:)  = 0
       kclbas(:)  = 0
       pcape(:)   = 0._realkind
       !cgj_realkind
       pdcdep(:)  = 0._realkind
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
    allocate( zscrad(itest) )      !jyj xcrad
    allocate( zscdep(itest) )      !jyj xcdepth
    allocate( gtrig1(itest) )
    allocate( zcape(itest) )
    allocate( iindex(klon) )
    allocate( ijsindex(itest) )
    !cgj
    allocate( zsdcdep(itest) )
    allocate( zworo(itest) )
    allocate( zoro(itest) )
    !
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
       zsdxdy(ji)    = pdxdy(jl)
       !cgj
       zworo(ji)     = pworo(jl)
       zoro(ji)     =  poro(jl)
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
    !
    call convect_trigger_funct( itest, klev,                              &
         zpres, zth, zthv, zthest,                 &
         zrv, zw, zz, zsdxdy,zworo,zoro,           &
         zsthlcl, zstlcl, zsrvlcl, zswlcl, zszlcl, &
         zsthvelcl, islcl, isdpl, ispbl, gtrig1,   &
         zcape, zscrad, zscdep, zsdcdep )   !jyj xcrad and xcdepth
    !
    do ji = 1, itest
       jl = ijsindex(ji)
       pcape(jl) = zcape(ji)
    end do
    !
    deallocate( zpres )
    deallocate( zz )
    deallocate( zth )
    deallocate( zthv )
    deallocate( zthest )
    deallocate( zrv )
    deallocate( zw )
    deallocate( zcape )
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
       deallocate( zscrad )     !jyj xcrad
       deallocate( zscdep )     !jyj xcdepth
       deallocate( islcl )
       deallocate( isdpl )
       deallocate( ispbl )
       deallocate( gtrig1 )
       deallocate( iindex )
       deallocate( ijsindex )
       !cgj
       deallocate( zsdcdep )
       deallocate( zworo )
       deallocate( zoro )
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
    allocate( zu(iconv,iks) )
    allocate( zv(iconv,iks) )
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
    allocate( zcrad(iconv) )      !jyj xcrad
    allocate( zcdep(iconv) )      !jyj xcdepth
    !cgj
    allocate( zdcdep(iconv) )
    !
    ! updraft variables
    !
    allocate( zumf(iconv,iks) )
    allocate( zuer(iconv,iks) )
    allocate( zudr(iconv,iks) )
    allocate( zudrf(iconv,iks) )
    allocate( zupr(iconv,iks) )
    allocate( zuthl(iconv,iks) )
    allocate( zuthv(iconv,iks) )
    allocate( zurw(iconv,iks) )
    allocate( zurc(iconv,iks) )
    allocate( zuri(iconv,iks) )
    allocate( zurr(iconv,iks) )
    allocate( zurs(iconv,iks) )
    allocate( zutpr(iconv) )
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
          zu(ji,jk)     = put(jl,jk)
          zv(ji,jk)     = pvt(jl,jk)
       end do
    end do
    if ( osettadj ) then
       allocate( ztimed(iconv) )
       do ji = 1, iconv
          jl = ijindex(ji)
          ztimed(ji) = ptimec(jl)
       end do
    end if
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
       zcrad(ji)     = zscrad(jl)    !jyj xcrad
       zcdep(ji)     = zscdep(jl)    !jyj xcdepth
       !cgj
       zdcdep(ji)    = zsdcdep(jl)
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
    deallocate( zscrad )      !jyj xcrad
    deallocate( zscdep )      !jyj xcdepth
    !cgj
    deallocate( zsdcdep )
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
    !
    !           4.     compute updraft properties 
    !                   ----------------------------
    !
    !           4.1    set mass flux at lcl ( here a unit mass flux with w = 1 m/s ) 
    !                   -------------------------------------------------------------
    !
    do ji = 1, iconv
       jk = ilcl(ji) - 1
       zmflcl(ji) = zpres(ji,jk) / ( xrd * ztt(ji,jk) *                &
            ( 1._realkind + zeps * zrvlcl(ji) ) ) * xpi * zcrad(ji) * zcrad(ji) 
    end do
    !
    deallocate( zcph )
    deallocate( zlv )
    deallocate( zls )
    !
    !
    call convect_updraft( iconv, klev,                                     &
         kice, zpres, zdpres, zz, zthl, zthv, zthes, zrw, &
         zthlcl, ztlcl, zrvlcl, zwlcl, zzlcl, zthvelcl,   & 
         zmflcl, gtrig1, ilcl, idpl, ipbl,                &
         zumf, zuer, zudr, zuthl, zuthv, zurw,            &
         zurc, zuri, zurr, zurs, zupr,                    &
         zutpr, zcape, ictl, ietl, zcrad, zcdep)    
    !jyj xcrad and xcdepth
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
       allocate( ilfs(iconv) )
       allocate( idbl(iconv) )
       allocate( iml(iconv) )
       allocate( zdmf(iconv,iks) )
       allocate( zder(iconv,iks) )
       allocate( zddr(iconv,iks) )
       allocate( zdthl(iconv,iks) )
       allocate( zdrw(iconv,iks) )
       allocate( zlmass(iconv,iks) )
       do jk = ikb, ike
          zlmass(:,jk)  = zdxdy(:) * zdpres(:,jk) / xg  ! mass of model layer
       end do
       zlmass(:,ikb) = zlmass(:,ikb+1)
       allocate( zmixf(iconv) )
       allocate( ztpr(iconv) )
       allocate( zspr(iconv) )
       allocate( zdtevr(iconv) )
       allocate( zpref(iconv) )
       allocate( zdtevrf(iconv,iks) )
       allocate( zprlflx(iconv,iks) )
       allocate( zprsflx(iconv,iks) )
       !
       ! closure variables
       !
       allocate( ztimea(iconv) )
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
       !           5.1    compute advective time period and precipitation 
       !                   efficiency as a function of mean ambient wind (shear) 
       !                   --------------------------------------------------------
       !
       call convect_tstep_pref( iconv, klev,                          &
            zu, zv, zpres, zz, zdxdy, ilcl, ictl, &
            ztimea, zpref )
       !
       ! exclude convective downdrafts if desired
       if ( .not. odown ) zpref(:) = 1._realkind
       !
       ! compute the period during which convection is active
       ztimec(:) = max( 1800._realkind, min( 3600._realkind, ztimea(:) ) )
       !cgj050711
       !cgj cape adjustment time made inversely proportional to resolution centred on 25km grid box
       !ztimec(:) = ztimec(:)*sqrt(zdxdy(:)/xa25)
       !cgj050711
       ztimec(:) = real( int( ztimec(:) / pdtconv ),realkind ) * pdtconv
       ztimec(:) = max( pdtconv, ztimec(:) ) ! necessary if pdtconv > 1800
       if ( osettadj ) then
          ztimec(:) = max( pdtconv, ztimed(:) )
       end if
       !
       !
       !           5.2    compute melting level
       !                   ----------------------
       !
       iml(:) = ikb
       do jk = ike, ikb, -1
          where( ztt(:,jk) <= xtt )  iml(:) = jk
       end do
       !
       call convect_downdraft( iconv, klev,                               &
            kice, zpres, zdpres, zz, zth, zthes,       & 
            zrw, zrc, zri,                             &
            zpref, ilcl, ictl, ietl,                   &
            zuthl, zurw, zurc, zuri,                   &
            zdmf, zder, zddr, zdthl, zdrw,             &
            zmixf, zdtevr, ilfs, idbl, iml,            &
            zdtevrf,zcrad                              )   !jyj xcrad
       !
       !
       !           6.     adjust up and downdraft mass flux to be consistent
       !                   with precipitation efficiency relation.
       !                   --------------------------------------------------- 
       !
       call convect_precip_adjust( iconv, klev,                              &
            zpres,zumf, zuer, zudr, zupr, zutpr, zurw,&
            zdmf, zder, zddr, zdthl, zdrw,            &
            zpref, ztpr, zmixf, zdtevr,               &
            ilfs, idbl, ilcl, ictl, ietl,             &
            zdtevrf                                   )
       !
       !
       !           7.     determine adjusted environmental values assuming
       !                   that all available buoyant energy must be removed
       !                   within an advective time step ztimec.
       !                   ---------------------------------------------------
       !
       call convect_closure( iconv, klev,                                &
            zpres, zdpres, zz, zdxdy, zlmass,           &
            zthl, zth, zrw, zrc, zri, gtrig1,           &
            zthc, zrvc, zrcc, zric, zwsub,              &
            ilcl, idpl, ipbl, ilfs, ictl, iml,          &
            zumf, zuer, zudr, zuthl, zurw,              &
            zurc, zuri, zupr,                           &
            zdmf, zder, zddr, zdthl, zdrw,              &
            ztpr, zspr, zdtevr,                         &
            zcape, ztimec,                              &
            iftsteps,                                   &
            zdtevrf, zprlflx, zprsflx )
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
          zprlflx(:,jk) = zprlflx(:,jk) / ( xrholw * zdxdy(:) )
          zprsflx(:,jk) = zprsflx(:,jk) / ( xrholw * zdxdy(:) )
          !
       end do
       !
       zprlflx(:,ikb) = zprlflx(:,ikb+1)
       zprsflx(:,ikb) = zprsflx(:,ikb+1)
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
             zwork2(ji) = zwork2(ji) + ( zrvc(ji,jk) + zrcc(ji,jk) + zric(ji,jk) ) *   & ! moisture
                  (zpres(ji,jk) - zpres(ji,jkp)) / xg
             !    zwork2b(ji) = zwork2b(ji) + (                                             & ! energy
             !                                ( xcpd + xcpv * zrw(ji,jk) )* zthc(ji,jk)   - &
             !          ( xlvtt + ( xcpv - xcl ) * ( ztt(ji,jk) - xtt ) ) * zrcc(ji,jk)   - & 
             !          ( xlstt + ( xcpv - xcl ) * ( ztt(ji,jk) - xtt ) ) * zric(ji,jk) ) * & 
             !                                      (zpres(ji,jk) - zpres(ji,jkp)) / xg
          end do
       end do
       !
       ! budget error (compare integral to surface precip.)
       !
       do ji = 1, iconv
          if ( ztpr(ji) > 0._realkind) then
             jkp = ictl(ji) + 1
             zwork2(ji) = ( ztpr(ji) / zdxdy(ji) + zwork2(ji) ) * xg /                 &
                  ( zpres(ji,ikb+1) - zpres(ji,jkp) )
             !    zwork2b(ji) = ( ztpr(ji) / zdxdy(ji) *                                    &
             !       ( xlvtt + ( xcpv - xcl ) * ( ztt(ji,ikb) - xtt ) ) - zwork2b(ji) )     &
             !                                * xg / ( zpres(ji,ikb+1) - zpres(ji,jkp) )
          end if
       end do
       !
       ! apply uniform correction
       !
       do jk = jkm, ikb+1, -1
          do ji = 1, iconv
             if ( ztpr(ji) > 0._realkind .and. jk <= ictl(ji) ) then
                zrvc(ji,jk) = zrvc(ji,jk) - zwork2(ji)                                ! moisture
                !        zthc(ji,jk) = zthc(ji,jk) + zwork2b(ji) / ( xcpd + xcpv * zrw(ji,jk) )! energy
             end if
          end do
       end do
       !
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
             !
             pprlflx(jl,jk) = zprlflx(ji,jk)
             pprsflx(jl,jk) = zprsflx(ji,jk)
          end do
       end do
       !
       !
       !
       !
       !           8.3    convective rainfall tendency
       !                   ----------------------------
       !
       ! liquid and solid surface rainfall tendency in m/s
       ztpr(:)   = ztpr(:) / ( xrholw * zdxdy(:) ) ! total surf precip
       zspr(:)   = zspr(:) / ( xrholw * zdxdy(:) ) ! solid surf precip
       ztpr(:)   = ztpr(:) - zspr(:) ! compute liquid part
       !
       do ji = 1, iconv
          jl = ijindex(ji)
          pprlten(jl) = ztpr(ji)
          pprsten(jl) = zspr(ji)
          !cgj
          pdcdep(jl)  = zdcdep(ji)
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
          indexcv(jl)= 2            !jyj deep convection index
       end do
       !jyj qvap for output--------------------------------------
       do jk = ikb, ike
          do ji = 1, iconv

             jl = ijindex(ji)
             if( jk<idpl(ji) ) then        !below conv source layer
                qvap(jl,jk)=prvt(jl,jk)
                qliq(jl,jk)=0._realkind
                qice(jl,jk)=0._realkind
             elseif( jk>=idpl(ji) .and. jk<kclbas(jl) ) then   !in source mixing layer
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
       !
       !           8.4    set convective counter
       !                   ----------------------
       !
       ! compute convective counter for just activated convective
       ! grid points
       ! if the advective time period is less than specified
       ! minimum for convective period, allow feedback to occur only
       ! during advective time
       !
       ztime(:) = 1._realkind
       zwork2(:) = 0._realkind
       do ji = 1, iconv
          jl = ijindex(ji)
          ztime(jl)  =  ztimec(ji)
          zwork2(jl) =  ztimea(ji)
          zwork2(jl) =  min( zwork2(jl), ztime(jl) )
          zwork2(jl) =  max( zwork2(jl), pdtconv )
          if ( gtrig(jl) )  kcount(jl) = int( zwork2(jl) / pdtconv )
          if ( gtrig(jl) .and. pprlten(jl)<1.e-14_realkind ) kcount(jl) = 0
       end do
       !
       !
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
          call convect_chem_transport( iconv, klev, kch1, zch1, zch1c,      &
               idpl, ipbl, ilcl, ictl, ilfs, idbl,  &
               zumf, zuer, zudr, zdmf, zder, zddr,  &
               ztimec, zdxdy, zmixf, zlmass, zwsub, &
               iftsteps )
          !
          do jk = ikb, ike
             do jn = 1, kch1
                zch1c(:,jk,jn) = ( zch1c(:,jk,jn)- zch1(:,jk,jn) ) / ztimec(:)
             end do
          end do
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
             if ( ztpr(ji) > 0._realkind) then
                jkp = ictl(ji) + 1
                zwork3(ji,:) = zwork3(ji,:) *                                     &
                     xg / ( zpres(ji,ikb+1) - zpres(ji,jkp) )
             end if
          end do
          !
          ! apply uniform correction but assure positive mass at each level
          !
          do jk = jkm, ikb+1, -1
             do ji = 1, iconv
                if ( ztpr(ji) > 0._realkind .and. jk <= ictl(ji) ) then
                   zch1c(ji,jk,:) = zch1c(ji,jk,:) - zwork3(ji,:)
                   zch1c(ji,jk,:) = max( zch1c(ji,jk,:), -zch1(ji,jk,:)/ztimec(ji) )
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
          zdmf(:,jk)  = zdmf(:,jk) / zdxdy(:)
!cgjudr
          zudrf(:,jk) = zudr(:,jk) / zlmass(:,jk)
!cgjudr
       end do
       zwork2(:) = 1._realkind
       where ( pprlten(:)<1.e-14_realkind ) zwork2(:) = 0._realkind
       do jk = ikb, ike
          do ji = 1, iconv
             jl = ijindex(ji)
             pumf(jl,jk) = zumf(ji,jk) * zwork2(jl)
             pdmf(jl,jk) = zdmf(ji,jk) * zwork2(jl)
!cgjudr
             pudr(jl,jk) = zudr(ji,jk)
             pudrf(jl,jk) =zudrf(ji,jk)
!cgjudr
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
       deallocate( zdthl )
       deallocate( zdrw )
       deallocate( zlmass )
       deallocate( zmixf )
       deallocate( ztpr )
       deallocate( zspr )
       deallocate( zdtevr )
       deallocate( zpref )
       deallocate( iml )
       deallocate( ilfs )
       deallocate( idbl )
       deallocate( zdtevrf )
       deallocate( zprlflx )
       deallocate( zprsflx )
       !
       !   closure variables
       !
       deallocate( ztimea )
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
    deallocate( zu )
    deallocate( zv )
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
    deallocate( zcrad )     !jyj xcrad
    deallocate( zcdep )     !jyj xcdepth
    !cgj
    deallocate( zdcdep )
    !
    ! updraft variables
    !
    deallocate( zumf )
    deallocate( zuer )
    deallocate( zudr )
    deallocate( zudrf )
    deallocate( zuthl )
    deallocate( zuthv )
    deallocate( zurw )
    deallocate( zurc )
    deallocate( zuri )
    deallocate( zurr )
    deallocate( zurs )
    deallocate( zupr )
    deallocate( zutpr )
    deallocate( zthlcl )
    deallocate( ztlcl )
    deallocate( zrvlcl )
    deallocate( zwlcl )
    deallocate( zzlcl )
    deallocate( zthvelcl )
    deallocate( zmflcl )
    deallocate( zcape )
    if ( osettadj ) deallocate( ztimed )
    !
    ! work arrays
    !
    deallocate( iindex )
    deallocate( ijindex )
    deallocate( ijsindex )
    deallocate( gtrig1 )
    !
    !
  end subroutine deep_convection


  subroutine convect_updraft( klon, klev,                                      &
       kice, ppres, pdpres, pz, pthl, pthv, pthes, prw, &
       pthlcl, ptlcl, prvlcl, pwlcl, pzlcl, pthvelcl,   &
       pmflcl, otrig, klcl, kdpl, kpbl,                 &
       pumf, puer, pudr, puthl, puthv, purw,            &
       purc, puri, purr, purs, pupr,                    &
       putpr, pcape, kctl, ketl, pcrad, pcdep ) 
    !jyj xcrad and xcdepth

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
    !
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
    real(kind=realkind), dimension(klon),     intent(in) :: pcrad  ! cloud radius   !jyj xcrad
    real(kind=realkind), dimension(klon),     intent(in) :: pcdep  ! cloud depth    !jyj xcdepth
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
    real(kind=realkind), dimension(klon,klev), intent(out):: purr  ! liquid precipit. (kg/kg)
    ! produced in  model layer
    real(kind=realkind), dimension(klon,klev),   intent(out)::purs ! solid precipit. (kg/kg)
    ! produced in  model layer
    real(kind=realkind), dimension(klon,klev),   intent(out)::pupr ! updraft precipitation in
    ! flux units (kg water / s)
    real(kind=realkind), dimension(klon),     intent(out):: putpr  ! total updraft precipitation
    ! in flux units (kg water / s)
    real(kind=realkind), dimension(klon),     intent(out):: pcape  ! available potent. energy
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike  ! horizontal and vertical loop bounds
    integer :: ji             ! horizontal loop index
    integer :: jk, jkp, jkm, jk1, jk2, jkmin  ! vertical loop index
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
    pupr(:,:)  = 0._realkind
    purr(:,:)  = 0._realkind
    purs(:,:)  = 0._realkind
    putpr(:)   = 0._realkind
    zuw1(:)    = pwlcl(:) * pwlcl(:)
    zuw2(:)    = 0._realkind
    ze1(:)     = 1._realkind
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
       !zwork6 is zero where jk=klcl-1
       !
       !
       !       4.     estimate condensate, l_v l_i, cph and theta_v at level k+1   
       !               ----------------------------------------------------------
       !
       zwork1(:) = purc(:,jk) + purr(:,jk)
       zwork2(:) = puri(:,jk) + purs(:,jk)
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
          purr(:,jkp) = 0.5_realkind * ( purc(:,jk) + purc(:,jkp) + puri(:,jk) + puri(:,jkp) )&
               * ( 1._realkind - exp( - xrconv  * zwork3(:) / zwork2(:) ) )
          pupr(:,jkp) = purr(:,jkp) * pumf(:,jk) ! precipitation rate ( kg water / s)
          putpr(:)    = putpr(:) + pupr(:,jkp)   ! total precipitation rate
          zwork2(:)   = purr(:,jkp) / max( 1.e-8_realkind, purc(:,jkp) + puri(:,jkp) )
          purr(:,jkp) = zwork2(:) * purc(:,jkp)          ! liquid precipitation
          purs(:,jkp) = zwork2(:) * puri(:,jkp)          ! solid precipitation
          !
          !
          !       7.     update r_c, r_i, enthalpy, r_w  for precipitation 
          !               -------------------------------------------------------
          !
          purw(:,jkp)  = purw(:,jk) - purr(:,jkp) - purs(:,jkp) 
          purc(:,jkp)  = purc(:,jkp) - purr(:,jkp)
          puri(:,jkp)  = puri(:,jkp) - purs(:,jkp)       
          puthl(:,jkp) = ( xcpd + purw(:,jkp) * xcpv ) * zut(:)                     &
               + ( 1._realkind + purw(:,jkp) ) * xg * pz(:,jkp)                    &
               - zlv(:) * purc(:,jkp) - zls(:) * puri(:,jkp)             
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
       ! zwork1(:) = xentr * pmflcl(:) * pdpres(:,jkp) / pcrad(:) ! rate of env. inflow   !jyj xcrad
       !mod
       zwork1(:) = xentr * xg / pcrad(:) * pumf(:,jk) * ( pz(:,jkp) - pz(:,jk) )        !jyj xcrad
       ! zwork1(:) = xentr * pumf(:,jk) * pdpres(:,jkp) / pcrad(:) ! rate of env. inflow  !jyj xcrad
       !mod
       zwork2(:) = 0._realkind
       where ( gwork1(:) ) zwork2(:) = 1._realkind
       where ( puthv(:,jkp) > pthv(:,jkp) )
       !cgj120811
       !cgj120811 Use buoyancy sorted method for ent/det instead of fixed 0.5 rate, comment next line
       !cgj180811   ze2=.5_realkind; zd2=.5_realkind  ! modif entrainment=detrainment, this avoids
       ! too large mass flux values at upper levels
       !cgj120811 allow Kain 2004 suggested limit for entrainment rates....uncomment next line
       !cgj   ze2(:)=max(ze2(:),0.5)                            !jyj need more test kain2004,eq(4)
       !cgj120811
       !   ze2(:)=max(ze2(:),0.5);zd2(:)=max(zd2(:),0.5)     !jyj need more test kain2004,eq(4)
          puer(:,jkp) = 0.5_realkind * zwork1(:) * ( ze1(:) + ze2(:) ) * zwork2(:)
          pudr(:,jkp) = 0.5_realkind * zwork1(:) * ( zd1(:) + zd2(:) ) * zwork2(:)
       elsewhere
          puer(:,jkp) = 0._realkind
          pudr(:,jkp) = zwork1(:) * zwork2(:)
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
          !
          !jyj dilute beg---------------------------------------------------------- 
          ! cape now calculated by an entrainment influenced dilute ascent
          ! profile rather than the original undilute ascent

          zwork2(:) = pthvelcl(:)* (1._realkind - zwork6(:)) + puthv(:,jk) * zwork6(:)
          zwork1(:) = (zwork2(:)+puthv(:,jkp)) / (pthv(:,jk) + pthv(:,jkp)) - 1._realkind
          !jyj dilute end---------------------------------------------------------- 

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
               / pumf(:,jkp) + puthl(:,jkp) - puthl(:,jk)
          purw(:,jkp)  = ( pumf(:,jk) * purw(:,jk) +                               &
               puer(:,jkp) * prw(:,jk) - pudr(:,jkp) * purw(:,jk) )    &
               / pumf(:,jkp) - purr(:,jkp) - purs(:,jkp)
          !    
          ze1(:) = ze2(:) ! update fractional entrainment/detrainment
          zd1(:) = zd2(:)
          !
       end where
       !
    end do
    !
    !       12.1    set otrig to false if cloud thickness < xcdepth
    !                or cape < 1
    !                ------------------------------------------------
    !
    do ji = 1, iie
       jk  = kctl(ji)
       otrig(ji) = pz(ji,jk) - pzlcl(ji) >= pcdep(ji)               &
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
             ! putpr(ji)    = putpr(ji) - ( purr(ji,jk) + purs(ji,jk) ) * pumf(ji,jkp)      
             putpr(ji)    = putpr(ji) - pupr(ji,jk)
             pudr(ji,jk)  = pdpres(ji,jk) * zwork1(ji)
             pumf(ji,jk)  = pumf(ji,jkp) - pudr(ji,jk)
             pupr(ji,jk)  = pumf(ji,jkp) * ( purr(ji,jk) + purs(ji,jk) )
             putpr(ji)    = putpr(ji) + pupr(ji,jk)
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
    !       13.   if cloud thickness is smaller than  3 km, no
    !              convection is allowed
    !              nota: for technical reasons, we stop the convection
    !                    computations in this case and do not go back to
    !                    trigger_funct to look for the next unstable lcl
    !                    which could produce a thicker cloud.
    !              ---------------------------------------------------
    !
    gwork6(:,:) = spread( otrig(:), dim=2, ncopies=klev )
    where ( .not. otrig(:) ) putpr(:) = 0._realkind
    where ( .not. gwork6(:,:) )
       pumf(:,:)  = 0._realkind
       pudr(:,:)  = 0._realkind
       puer(:,:)  = 0._realkind
       puthl(:,:) = pthl(:,:)
       purw(:,:)  = prw(:,:)
       pupr(:,:)  = 0._realkind
       purc(:,:)  = 0._realkind
       puri(:,:)  = 0._realkind
       purr(:,:)  = 0._realkind
       purs(:,:)  = 0._realkind
    end where
    !
  end subroutine convect_updraft


  subroutine convect_trigger_funct( klon, klev,                           &
       ppres, pth, pthv, pthes,              &
       prv, pw, pz, pdxdy,pworo,poro,        &
       pthlcl, ptlcl, prvlcl, pwlcl, pzlcl,  &
       pthvelcl, klcl, kdpl, kpbl, otrig,    &
       pcape, pcrad, pcdep, pdcdep )                  !jyj xcrad xcdepth

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
    !
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    integer, intent(in)                   :: klon      ! horizontal loop index
    integer, intent(in)                   :: klev      ! vertical loop index
    real(kind=realkind), dimension(klon),     intent(in) :: pdxdy     ! grid area
    !cgj
    real(kind=realkind), dimension(klon),     intent(in) :: pworo     ! orographic standard deviation
    real(kind=realkind), dimension(klon),     intent(in) :: poro      ! fraction of land 
    !cgj
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
    logical, dimension(klon),  intent(out):: otrig     ! logical mask for convection 
    integer, dimension(klon),  intent(inout):: klcl    ! contains vert. index of lcl
    integer, dimension(klon),  intent(inout):: kdpl    ! contains vert. index of dpl
    integer, dimension(klon),  intent(inout):: kpbl    ! contains index of source layer top
    real(kind=realkind), dimension(klon),     intent(out):: pcape     ! cape (j/kg) for diagnostics
    real(kind=realkind), dimension(klon),     intent(out):: pcrad     ! cloud radius (m)    !jyj xcrad
    real(kind=realkind), dimension(klon),     intent(out):: pcdep     ! cloud depth  (m)    !jyj xcdepth
    !cgj
    real(kind=realkind), dimension(klon),     intent(out):: pdcdep     ! deep conv cloud depth  (m) 
    !
    !       0.2   declarations of local variables :
    !
    integer :: jkk, jk, jkp, jkm, jkdl, jl, jkt, jt! vertical loop index
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
    real(kind=realkind), dimension(klon,klev):: zcap ! cape at every level for diagnostics
    !integer, dimension(klon) :: itop  ! work array to store highest test layer
    real(kind=realkind), dimension(klon) :: zwork1, zwork2, zwork3, zwork4    ! work arrays
    logical, dimension(klon) :: gtrig, gtrig2          ! local arrays for otrig
    logical, dimension(klon) :: gwork1                 ! work array

    !
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
    zcap(:,:)  = 0._realkind
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
       if ( pz(1,jk) - pz(1,ikb) < 12.e3_realkind ) jt = jk 
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
       ! than 3500 m  above soil level.
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
       !       3.     construct a mixed layer of at least 60 hpa (xzpbl)
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
          ! zthlcl(:)   = zthlcl(:)   / zdpthmix(:) 
          ! zrvlcl(:)   = zrvlcl(:)   / zdpthmix(:) 
          zthlcl(:)   = zthlcl(:)   / zdpthmix(:) +.3_realkind
          zrvlcl(:)   = zrvlcl(:)   / zdpthmix(:) +1.e-4_realkind
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

       if ( jkk>minval(ilcl(:)) ) exit     !jyj
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
       !
       !       6.     check to see if cloud is bouyant 
       !               --------------------------------
       !
       !      6.1    compute grid scale vertical velocity perturbation term zwork1
       !               -------------------------------------------------------------
       ! 
       !  normalize w grid scale to a 25 km refer. grid
       do ji = 1, iie
          jk  = ilcl(ji)
          jkm = jk - 1 
          !jyj xcrad beg---------------------------------------------
          !kain (2004), journal fo applied meteorology, 43, 170-181. eq(2)

          if(zzlcl(ji) <= 2000._realkind) then
             zwork1(ji)=0.02_realkind*zzlcl(ji)/2000._realkind
          else
             zwork1(ji)=0.02_realkind
          endif

          !calculate dlp using z instead of log(p)...
          zwork2(ji) = (zzlcl(ji)-pz(ji,jkm))/(pz(ji,jk)-pz(ji,jkm))

          zwork1(ji)=(pw(ji,jkm) + (pw(ji,jk) - pw(ji,jkm)) * zwork2(ji)) &
               *sqrt( pdxdy(ji) / xa25 ) - zwork1(ji)                   !jyj
          !cgj
          !cgj  modify zwork1, based on orographic standard deviation woro
          !cgj  this acts to reduce the vertical velocity at zlcl in mountainous
          !cgj  regions
          !cgj
          zwork4(ji)=sqrt(pworo(ji)/70.0_realkind)
          zwork4(ji)=min(2.5_realkind, max(1.0_realkind,zwork4(ji)))
          if(zwork1(ji)>0._realkind)zwork1(ji)=zwork1(ji)*zwork4(ji)
          !cgj
          !kain (2004), journal fo applied meteorology, 43, 170-181. eq(6)
          !cgj tuned for rca3.5
          !cgj010812if(zwork1(ji)<0.0_realkind) then
          !cgj010812   pcrad(ji) = 1000._realkind
          !cgj010812elseif(zwork1(ji)>0.1_realkind) then
          !cgj010812   pcrad(ji) = 2000._realkind
          !cgj010812else
          !cgj010812   pcrad(ji) = 1000._realkind + (1000._realkind*(zwork1(ji)/0.1_realkind))
          !cgj010812endif
          pcrad(ji) = 1500.0_realkind
          !cgj010812if(zwork1(ji)<-0.05_realkind) then
          !cgj010812   pcrad(ji) = 100._realkind
          !cgj010812elseif(zwork1(ji)>=-0.05_realkind .and. zwork1(ji)<=0._realkind) then
          !cgj010812   pcrad(ji) = 300._realkind + (zwork1(ji)*4000._realkind)
          !cgj010812elseif(zwork1(ji)>0.2_realkind) then
          !cgj010812   pcrad(ji) = 2000._realkind
          !cgj010812else
          !cgj010812   !cgjoriginal           pcrad(ji) = 500. * (1.0 + (3.0*(zwork1(ji)/0.1)) )
          !cgj010812   pcrad(ji) = 300._realkind + (1700._realkind*(zwork1(ji)/0.2_realkind))
          !cgj010812endif
          !jyj xcrad end---------------------------------------------

          !cgj        zwork1(ji) =  ( pw(ji,jkm)  + ( pw(ji,jk) - pw(ji,jkm) ) * zdp(ji) )  &
          !cgj                           * sqrt( pdxdy(ji) / xa25 )
          !                         - 0.02 * zzlcl(ji) / xzlcl ! avoid spurious convection
       end do
       ! compute sign of normalized grid scale w
       zwork2(:) = sign( 1._realkind, zwork1(:) ) 
       zwork1(:) = xwtrig * zwork2(:) * abs( zwork1(:) ) **(1.0_realkind/3.0_realkind) &!0.333_realkind       &
            * ( xp00 / zplcl(:) ) ** zrdocp
       !
       !       6.2    compute parcel vertical velocity at lcl
       !               ---------------------------------------
       !                   
       do ji = 1, iie
          jkdl = idpl(ji)
          zwork3(ji) = xg * zwork1(ji) * ( zzlcl(ji) - pz(ji,jkdl) )       &
               / ( pthv(ji,jkdl) + zthvelcl(ji) )
       end do
       where( gwork1(:) )
       !cgjstdrca4   zwlcl(:)  = 1._realkind + .5_realkind * zwork2(:) * sqrt( abs( zwork3(:) ) ) 
          zwlcl(:)  = 1._realkind + .5_realkind * zwork2(:) * zwork4(:) * sqrt( abs( zwork3(:) ) ) 
          !cgj
       !cgjstdrca4   zwlcl(:)  = min(zwlcl(:),3._realkind)
          zwlcl(:)  = min(zwlcl(:),5._realkind)
          !cgj
          gtrig(:)  = zthvlcl(:) - zthvelcl(:) + zwork1(:) > 0._realkind .and.       &
               zwlcl(:) > 0._realkind 
       end where
       !
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
                  ( pthes(ji,jk) + pthes(ji,jl) ) - 1._realkind ) * ( pz(ji,jk) - pz(ji,jl) )
             if ( jl < ilcl(ji) ) zwork1(ji) = 0._realkind
             zcape(ji)  = zcape(ji) + zwork1(ji)
             zcap(ji,jkk) = zcap(ji,jkk) + xg * max( 0._realkind, zwork1(ji) ) ! actual cape
             zwork2(ji) = xnhgam * xg * zcape(ji) + 1.05_realkind * zwlcl(ji) * zwlcl(ji)
             ! the factor 1.05 takes entrainment into account
             zwork2(ji) = sign( 1._realkind, zwork2(ji) )
             zwork3(ji) = zwork3(ji) + min(0._realkind, zwork2(ji) )
             zwork3(ji) = max( -1._realkind, zwork3(ji) )
             ! nota, the factors zwork2 and zwork3 are only used to avoid
             ! if and goto statements, the difficulty is to extract only
             ! the level where the criterium is first fullfilled
             ztop(ji)   = pz(ji,jl) * .5_realkind * ( 1._realkind + zwork2(ji) ) * ( 1._realkind + zwork3(ji) ) + &
                  ztop(ji) * .5_realkind * ( 1._realkind - zwork2(ji) )
          end do
       end do
       !jyj xcdepth beg---------------------
       ! specifying minimum cloud depth as a function of tlcl
       ! kain, 2004, journal fo applied meteorology, 43, 170-181. eq(7)
       do ji = 1, iie
!cgj220811          if (ztlcl(ji)>293.0_realkind)  then
!cgj220811             pcdep(ji) = 4000.0_realkind
!cgj220811             pcdep(ji) = 3000.0_realkind
!cgj220811          else if( ztlcl(ji) <= 293.0_realkind .and. ztlcl(ji) >= 273.0_realkind) then
!cgj220811             pcdep(ji) = 2000.0_realkind + 100._realkind*(ztlcl(ji)-273.0_realkind)
!cgj220811          else if( ztlcl(ji) < 273.0_realkind) then
!cgj220811             pcdep(ji) = 2000.0_realkind
!cgj220811          end if
       !cgj010812 pcdep(ji)=3000.0_realkind
       if(poro(ji)<=0.001_realkind)then
        pcdep(ji)=3500.0_realkind
       else
        pcdep(ji)=500._realkind + (2500._realkind*( (1._realkind - poro(ji))**3))
       endif
!        pcdep(ji)=500._realkind
!cgj060912
       end do
       !jyj xcdepth end----------------------
       !
       !
       where( ztop(:) - zzlcl(:) >= pcdep(:)  .and. gtrig(:) .and. gtrig2(:) )
          gtrig2(:)   = .false.
          otrig(:)    = gtrig(:)     ! we  select the first departure level
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
          pdcdep(:)   = ztop(:) - zzlcl(:)
       end where
       !
       if ( count(.not.otrig(:) ) == 0 ) exit   !jyj
    end do
    !
    do ji = 1, iie
       pcape(ji) = maxval( zcap(ji,:) ) ! maximum cape for diagnostics
    end do
    !
    !
  end subroutine convect_trigger_funct


  subroutine convect_precip_adjust( klon, klev,                        &
       ppres, pumf, puer, pudr,           &
       pupr, putpr, purw,                 &
       pdmf, pder, pddr, pdthl, pdrw,     &
       ppref, ptpr, pmixf, pdtevr,        &
       klfs, kdbl, klcl, kctl, ketl,      &
       pdtevrf )

    !
    !! adjust up- and downdraft mass fluxes to be consistent with the
    !!     mass transport at the lfs given by the precipitation efficiency
    !!     relation. 
    !!
    !!
    !!    purpose                                                       
    !!    -------
    !!      the purpose of this routine is to adjust up- and downdraft mass
    !!      fluxes below the lfs to be consistent with the precipitation
    !!      efficiency relation
    !!
    !!
    !!
    !!  method
    !!    ------
    !!      
    !!
    !!    external
    !!    --------
    !!     none
    !!     
    !!
    !!    implicit arguments
    !!    ------------------
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_precip_adjust)
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  04/10/97
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
    integer,                    intent(in) :: klon  ! horizontal dimension
    integer,                    intent(in) :: klev  ! vertical dimension
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres ! pressure (pa)
    real(kind=realkind), dimension(klon,klev), intent(in) :: purw  ! updraft total water (kg/kg) 
    real(kind=realkind), dimension(klon),      intent(in) :: putpr ! updraft  total precipit. (kg/s
    real(kind=realkind), dimension(klon),      intent(in) :: ppref ! precipitation efficiency
    real(kind=realkind), dimension(klon),      intent(in) :: pmixf ! critical mixed fraction at lcl
    integer, dimension(klon),   intent(in) :: klcl  ! contains vert. index of lcl
    integer, dimension(klon),   intent(in) :: kctl  ! contains vert. index of ctl
    integer, dimension(klon),   intent(in) :: ketl  ! contains vert. index of equilibrium 
    ! (zero buoyancy) level 
    integer, dimension(klon),  intent(inout) :: klfs ! contains vert. index of lfs
    integer, dimension(klon),  intent(inout) :: kdbl ! contains vert. index of dbl
    !
    real(kind=realkind), dimension(klon),      intent(inout) :: pdtevr ! total downdraft evaporation
    ! rate at lfs   
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pdtevrf! downdraft evaporation rate
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pumf   ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: puer   ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pudr   ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pupr   ! updraft  precipit. (kg/s)     
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pdmf   ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pder   ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pddr   ! downdraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pdthl  ! downdraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pdrw   ! downdraft total water (kg/kg)
    !
    real(kind=realkind), dimension(klon),     intent(out)   :: ptpr    ! total precipitation (kg/s) 
    ! = downdraft precipitation
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike        ! horizontal + vertical loop bounds
    integer :: jk, jkt1, jkt2, jkt3 ! vertical loop index
    integer :: ji                   ! horizontal loop index
    !
    integer, dimension(klon) :: iprl
    real(kind=realkind), dimension(klon)    :: zwork1, zwork2, zwork3,     &
         zwork4, zwork5, zwork6 ! work arrays
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        0.3   set loop bounds
    !              ---------------
    !
    ikb  = 1 + jcvexb 
    ike  = klev - jcvext 
    iie  = klon
    jkt1 = maxval( klfs(:) )
    jkt2 = maxval( kctl(:) )
    jkt3 = minval( klcl(:) )
    !
    !
    !        1.    set some output variables for columns where no downdraft 
    !              exists. exit if there is no downdraft at all.
    !              --------------------------------------------------------
    !
    iprl(:) = ikb
    ptpr(:) = 0._realkind
    !
    where ( abs(pdtevr(:)) < 1.e-14_realkind )
       ptpr(:)    = putpr(:)  ! no downdraft evaporation => no downdraft, all
       ! precipitation occurs in updraft
    end where
    if ( count( pdtevr(:) > 1.e-14_realkind ) == 0 ) return ! exit routine if no downdraft exists
    !
    !       2.     the total mass transported from the updraft to the down-  
    !               draft at the lfs must be consistent with the three water
    !               budget terms :
    !               ---------------------------------------------------------
    !
    !       2.1    downdraft evaporation rate at the dbl. the evaporation
    !               rate in downdraft must be consistent with precipitation
    !               efficiency relation.
    !               --------------------------------------------------------
    !
    !
    do ji = 1, iie
       jk = klfs(ji)
       zwork1(ji) = pdtevr(ji) / min( -1.e-1_realkind, pdmf(ji,jk) )
       zwork6(ji) = pdmf(ji,jk)
    end do
    !
    !       2.2    some preliminar computations for downdraft = total 
    !               precipitation rate. the precipitation is evaluated in 
    !               a layer thickness dp=xusrdpth=165 hpa above the lcl.
    !               the difference between updraft precipitation and downdraft
    !               precipitation (updraft supply rate) is used to drive the
    !               downdraft through evaporational cooling.
    !               --------------------------------------------------------
    !
    do ji = 1, iie
       jk = klcl(ji)
       zwork5(ji) = ppres(ji,jk)
    end do
    !
    ptpr(:) = 0._realkind
    do jk = jkt3, jkt2
       where ( jk >= klcl(:) .and. ppres(:,jk) >= zwork5(:) - xusrdpth )
          ptpr(:) = ptpr(:) + pupr(:,jk)
          iprl(:) = jk
       end where
    end do
    iprl(:) = min( ketl(:), iprl(:) )
    !
    do ji = 1, iie
       jk = iprl(ji)
       ptpr(ji) = pumf(ji,jk+1) * purw(ji,jk+1) + ptpr(ji) 
    end do
    !
    ptpr(:) = ppref(:) * min( putpr(:), ptpr(:) )
    zwork4(:) = putpr(:) - ptpr(:) 
    !
    !
    !       2.3    total amount of precipitation that falls out of the up-
    !               draft between the lcl and the lfs.
    !               condensate transfer from up to downdraft at lfs
    !               ---------------------------------------------------------
    !
    zwork5(:) = 0._realkind
    do jk = jkt3, jkt1
       where ( jk >= klcl(:) .and. jk <= klfs(:) )
          zwork5(:) = zwork5(:) +  pupr(:,jk)
       end where
    end do
    !
    do ji = 1, iie
       jk = klfs(ji)
       zwork2(ji) = ( 1._realkind - ppref(ji) ) * zwork5(ji) *                     &
            ( 1._realkind - pmixf(ji) ) / max( 1.e-1_realkind, pumf(ji,jk) )
    end do
    !
    !
    !       2.4    increase the first guess downdraft mass flux to satisfy
    !               precipitation efficiency relation.
    !               if downdraft does not evaporate any water at the dbl for  
    !               the specified relative humidity, or if the corrected mass 
    !               flux at the lfs is positive no downdraft is allowed
    !               ---------------------------------------------------------
    !    
    !
    zwork1(:) = zwork4(:) / ( zwork1(:) + zwork2(:) + 1.e-8_realkind ) 
    zwork2(:) = zwork1(:) / min( -1.e-1_realkind, zwork6(:) ) ! ratio of budget consistent to actual dmf
    !
    zwork3(:) = 1._realkind
    zwork6(:) = 1._realkind
    where ( zwork1(:) > 0._realkind .or. pdtevr(:) < 1._realkind ) 
       kdbl(:)   = ikb
       klfs(:)   = ikb
       pdtevr(:) = 0._realkind 
       zwork2(:) = 0._realkind
       zwork3(:) = 0._realkind
       zwork6(:) = 0._realkind
    end where
    !
    do jk = ikb, jkt1   
       pdmf(:,jk)  = pdmf(:,jk)  * zwork2(:)
       pder(:,jk)  = pder(:,jk)  * zwork2(:)  
       pddr(:,jk)  = pddr(:,jk)  * zwork2(:)  
       pdtevrf(:,jk) = pdtevrf(:,jk)* zwork2(:)  
       pdrw(:,jk)  = pdrw(:,jk)  * zwork3(:)  
       pdthl(:,jk) = pdthl(:,jk) * zwork3(:)  
    end do
    zwork4(:) = zwork2(:)
    !
    !
    !       3.     increase updraft mass flux, mass detrainment rate, and water  
    !               substance detrainment rates to be consistent with the transfer
    !               of the estimated mass from the up- to the downdraft at the lfs
    !               --------------------------------------------------------------
    !
    do ji = 1, iie
       jk = klfs(ji)
       zwork2(ji) = ( 1._realkind - zwork6(ji) ) + zwork6(ji) *                   &
            ( pumf(ji,jk) - ( 1._realkind - pmixf(ji) ) * zwork1(ji) ) / &
            max( 1.e-1_realkind, pumf(ji,jk) )
    end do
    !
    !
    jkt1  = maxval( klfs(:) )  ! value of klfs might have been reset to ikb above
    do jk = ikb, jkt1
       do ji = 1, iie
          if ( jk <= klfs(ji) ) then
             pumf(ji,jk)  = pumf(ji,jk)  * zwork2(ji) 
             puer(ji,jk)  = puer(ji,jk)  * zwork2(ji)
             pudr(ji,jk)  = pudr(ji,jk)  * zwork2(ji)
             pupr(ji,jk)  = pupr(ji,jk)  * zwork2(ji)
          end if
       end do
    end do
    !
    !
    !       4.     increase total = downdraft precipitation and evaporation rate
    !               -------------------------------------------------------------
    !
    where ( pdtevr(:) > 0._realkind )
       pdtevr(:)  = pdtevr(:) * zwork4(:)
       ptpr(:)    = ptpr(:) + ppref(:) * zwork5(:) * ( zwork2(:) - 1._realkind )
    elsewhere
       ptpr(:)    = putpr(:)
    end where
    !
    !
  end subroutine convect_precip_adjust



  subroutine convect_downdraft( klon, klev,                                &
       kice, ppres, pdpres, pz, pth, pthes,       &
       prw, prc, pri,                             &
       ppref, klcl, kctl, ketl,                   &
       puthl, purw, purc, puri,                   &
       pdmf, pder, pddr, pdthl, pdrw,             &
       pmixf, pdtevr, klfs, kdbl, kml,            &
       pdtevrf, pcrad )    !jyj xcrad

    !
    !! compute downdraft properties from lfs to dbl. 
    !!
    !!
    !!    pdrpose                                                       
    !!    -------
    !!      the purpose of this routine is to determine downdraft properties
    !!      ( mass flux, thermodynamics ) 
    !!
    !!
    !!  method
    !!    ------
    !!      computations are done at every model level starting from top.
    !!      the use of masks allows to optimise the inner loops (horizontal loops).
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
    !!
    !!      module modd_cst
    !!          xg                 ! gravity constant
    !!          xpi                ! pi
    !!          xp00               ! reference pressure
    !!          xrd, xrv           ! gaz  constants for dry air and water vapor
    !!          xcpd               ! cpd (dry air)
    !!          xcpv, xcl, xci     ! cp of water vapor, liquid water and ice
    !!          xtt                ! triple point temperature
    !!          xlvtt, xlstt       ! vaporisation/sublimation heat at xtt
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_downdraft)
    !!      kain and fritsch, 1993, meteor. monographs, vol.
    !!
    !!    author
    !!    ------
    !!      p. bechtold        laboratoire d'aerologie 
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  04/10/97
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
    integer,                    intent(in) :: klon  ! horizontal dimension
    integer,                    intent(in) :: klev  ! vertical dimension
    integer,                    intent(in) :: kice  ! flag for ice ( 1 = yes,
    !                0 = no ice )
    real(kind=realkind), dimension(klon,klev), intent(in) :: pth   ! grid scale theta        
    real(kind=realkind), dimension(klon,klev), intent(in) :: pthes ! grid scale saturated theta_e 
    real(kind=realkind), dimension(klon,klev), intent(in) :: prw   ! grid scale total water  
    ! mixing ratio 
    real(kind=realkind), dimension(klon,klev), intent(in) :: prc   ! grid scale r_c (cloud water) 
    real(kind=realkind), dimension(klon,klev), intent(in) :: pri   ! grid scale r_i (cloud ice) 
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres ! pressure (pa)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pdpres! pressure difference between 
    ! bottom and top of layer (pa) 
    real(kind=realkind), dimension(klon,klev), intent(in) :: pz    ! level height (m)
    integer, dimension(klon),   intent(in) :: klcl  ! contains vert. index of lcl
    integer, dimension(klon),   intent(in) :: kctl  ! contains vert. index of ctl 
    integer, dimension(klon),   intent(in) :: ketl  ! contains vert. index of 
    ! equilibrium (zero buoyancy) level 
    integer, dimension(klon),   intent(in) :: kml   ! " vert. index of melting level
    real(kind=realkind), dimension(klon,klev), intent(in) :: puthl ! updraft enthalpy (j/kg)      
    real(kind=realkind), dimension(klon,klev), intent(in) :: purw  ! updraft total water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in) :: purc  ! updraft r_c (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in) :: puri  ! updraft r_i (kg/kg)
    real(kind=realkind), dimension(klon),      intent(in) :: ppref ! precipitation efficiency
    real(kind=realkind), dimension(klon),      intent(in) :: pcrad ! jyj xcrad
    !
    !
    real(kind=realkind), dimension(klon,klev), intent(out):: pdmf   ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: pder   ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: pddr   ! downdraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: pdthl  ! downdraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(out):: pdrw   ! downdraft total water (kg/kg)
    real(kind=realkind), dimension(klon),      intent(out):: pmixf  ! mixed fraction at lfs
    real(kind=realkind), dimension(klon),      intent(out):: pdtevr ! total downdraft evaporation
    ! rate at lfs (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(out):: pdtevrf! downdraft evaporation rate
    integer, dimension(klon),  intent(out):: klfs    ! contains vert. index of lfs 
    integer, dimension(klon),  intent(out):: kdbl    ! contains vert. index of dbl   
    !
    !       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike     ! horizontal + vertical loop bounds
    integer :: jk, jkp, jkm, jkt ! vertical loop index
    integer :: ji, jl            ! horizontal loop index
    integer :: jiter          ! iteration loop index
    real(kind=realkind)    :: zcpord, zrdocp ! c_pd / r_d, r_d / c_pd
    real(kind=realkind)    :: zeps           ! r_d / r_v
    real(kind=realkind)    :: zepsa, zcvocd  ! r_v / r_d, c_pv / c_pd
    !
    integer, dimension(klon) :: iddt      ! top level of detrainm. layer
    real(kind=realkind), dimension(klon)    :: zthe      ! environm. theta_e (k)
    real(kind=realkind), dimension(klon)    :: zdt, zdtp ! downdraft temperature (k)
    real(kind=realkind), dimension(klon)    :: zcph      ! specific heat c_ph 
    real(kind=realkind), dimension(klon)    :: zlv, zls  ! latent heat of vaporis., sublim.       
    real(kind=realkind), dimension(klon)    :: zddt      ! thickness (hpa) of detrainm. layer
    real(kind=realkind), dimension(klon)    :: zpi       ! pi=(p0/p)**(rd/cpd)  
    real(kind=realkind), dimension(klon)    :: zwork1, zwork2, zwork3, zwork4,  &
         zwork5                  ! work arrays 
    logical, dimension(klon) :: gwork1                         ! work array
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        0.3    set loop bounds
    !               ---------------
    !
    iie = klon
    ikb = 1 + jcvexb 
    ike = klev - jcvext 
    !
    !
    !       1.     initialize downdraft properties
    !               -------------------------------
    !
    zcpord     = xcpd / xrd
    zrdocp     = xrd / xcpd
    zeps       = xrd / xrv
    zepsa      = xrv / xrd
    zcvocd     = xcpv / xcpd
    pdmf(:,:)  = 0._realkind
    pder(:,:)  = 0._realkind
    pddr(:,:)  = 0._realkind
    pdrw(:,:)  = 0._realkind
    pdthl(:,:) = 0._realkind
    pdtevr(:)  = 0._realkind
    pmixf(:)   = 0._realkind
    zthe(:)    = 0._realkind
    zddt(:)    = pdpres(:,ikb+2)
    kdbl(:)    = ikb + 1
    klfs(:)    = ikb + 1
    iddt(:)    = kdbl(:) + 1
    !  
    !
    !       2.     determine the lfs by looking for minimum of environmental 
    !               saturated theta_e 
    !               ----------------------------------------------------------
    !
    zwork1(:) = 900._realkind   ! starting value for search of minimum envir. theta_e
    do jk = minval( klcl(:) ) + 2, maxval( ketl(:) )
       do ji = 1, iie
          gwork1(ji) = jk >= klcl(ji) + 2 .and. jk < ketl(ji)  
          if ( gwork1(ji) .and. zwork1(ji) > pthes(ji,jk) ) then
             klfs(ji)   = jk
             zwork1(ji) = min( zwork1(ji), pthes(ji,jk) )
          end if
       end do
    end do
    !
    !
    !       3.     determine the mixed fraction using environmental and updraft
    !               values of theta_e at lfs
    !               ---------------------------------------------------------   
    !
    do ji = 1, iie
       jk = klfs(ji)
       zpi(ji)    = ( xp00 / ppres(ji,jk) ) ** zrdocp
       ! compute updraft theta_e
       zwork3(ji) = purw(ji,jk) - purc(ji,jk) - puri(ji,jk)
       zdt(ji)    = pth(ji,jk) / zpi(ji) 
       zlv(ji)    = xlvtt + ( xcpv - xcl ) * ( zdt(ji) - xtt )                   
       zls(ji)    = xlstt + ( xcpv - xci ) * ( zdt(ji) - xtt )                   
       zcph(ji)   = xcpd + xcpv * purw(ji,jk)
       zdt(ji)    = ( puthl(ji,jk) - ( 1._realkind + purw(ji,jk) ) * xg * pz(ji,jk)       &
            + zlv(ji) * purc(ji,jk) + zls(ji) * puri(ji,jk) ) / zcph(ji)           
       zwork1(ji) = zdt(ji) * zpi(ji) ** ( 1._realkind - 0.28_realkind * zwork3(ji) )              &
            * exp( ( 3374.6525_realkind / zdt(ji) - 2.5403_realkind )     &
            * zwork3(ji) * ( 1._realkind + 0.81_realkind * zwork3(ji) ) )
       ! compute environmental theta_e
       zdt(ji)    = pth(ji,jk) / zpi(ji)
       zlv(ji)    = xlvtt + ( xcpv - xcl ) * ( zdt(ji) - xtt )                   
       zls(ji)    = xlstt + ( xcpv - xci ) * ( zdt(ji) - xtt )                   
       zwork3(ji) = prw(ji,jk) - prc(ji,jk) - pri(ji,jk)
       zcph(ji)   = xcpd + xcpv * prw(ji,jk)
       zwork2(ji) = zdt(ji) * zpi(ji) ** ( 1._realkind - 0.28_realkind * zwork3(ji) )              &
            * exp( ( 3374.6525_realkind / zdt(ji) - 2.5403_realkind )     &
            * zwork3(ji) * ( 1._realkind + 0.81_realkind * zwork3(ji) ) )
       ! compute mixed fraction
       pmixf(ji)  = max( 0._realkind, ( zwork1(ji) - pthes(ji,jk) ) )                   &
            / ( zwork1(ji) - zwork2(ji) + 1.e-14_realkind )
       pmixf(ji)  = max(0._realkind, min( 1._realkind, pmixf(ji) ) )
       zwork4(ji) = ppres(ji,jk)
    end do
    !
    !
    !       4.     estimate the effect of melting on the downdraft  
    !               ---------------------------------------------
    !
    zwork1(:) = 0._realkind
    ! use total solid precipitation
    !do jk = ikb + 1, ike
    !    zwork1(:) = zwork1(:) + purs(:,jk) ! total snow/hail content
    !end do
    !
    do ji = 1, iie
       jk  = klcl(ji)
       jkp = kctl(ji)
       zwork1(ji) = 0.5_realkind * ( purw(ji,jk) - purw(ji,jkp) )
    end do
    !
    ! temperature perturbation due to melting at lfs
    zwork3(:) = 0._realkind
    where( kml(:) > ikb + 2 )
       zwork3(:) = zwork1(:) * ( zls(:) - zlv(:) ) / zcph(:)
       zdt(:)    = zdt(:) - zwork3(:) * real(kice,realkind)
    end where
    !
    !
    !       5.     initialize humidity at lfs as a saturated mixture of
    !               updraft and environmental air
    !               -----------------------------------------------------    
    !
    do ji = 1, iie
       jk = klfs(ji)
       pdrw(ji,jk)  = pmixf(ji) * prw(ji,jk) + ( 1._realkind - pmixf(ji) ) * purw(ji,jk)
       zwork2(ji)   = pdrw(ji,jk) - ( 1._realkind - pmixf(ji) )                          &
            * ( purc(ji,jk) + puri(ji,jk) )
    end do
    !
    !
    !       6.1    determine the dbl by looking for level where the envir.
    !               theta_es at the lfs corrected by melting effects  becomes
    !               larger than envir. value
    !               ---------------------------------------------------------
    !
    ! compute satur. mixing ratio for melting corrected temperature
    call convect_satmixratio( klon, zwork4, zdt, zwork3, zlv, zls, zcph )  
    !
    ! compute envir. saturated theta_e for melting corrected temperature
    zwork1(:) = min( zwork2(:), zwork3(:) )
    zwork3(:) = zwork3(:) * zwork4(:) / ( zwork3(:) + zeps ) ! sat. pressure
    zwork3(:) = log( zwork3(:) / 613.3_realkind )
    ! dewp point temperature
    zwork3(:) = ( 4780.8_realkind - 32.19_realkind * zwork3(:) ) / ( 17.502_realkind - zwork3(:) )
    ! adiabatic saturation temperature
    zwork3(:) = zwork3(:) - ( .212_realkind + 1.571e-3_realkind * ( zwork3(:) - xtt )          &
         - 4.36e-4_realkind * ( zdt(:) - xtt ) ) * ( zdt(:) - zwork3(:) )
    zwork4(:) = sign(0.5_realkind, zwork2(:) - zwork3(:) )
    zdt(:)    = zdt(:) * ( .5_realkind + zwork4(:) ) + ( .5_realkind - zwork4(:) ) * zwork3(:) 
    zwork2(:) = zdt(:) * zpi(:) ** ( 1._realkind - 0.28_realkind * zwork2(:) )                 &
         * exp( ( 3374.6525_realkind / zdt(:) - 2.5403_realkind )     &
         * zwork1(:) * ( 1._realkind + 0.81_realkind * zwork1(:) ) )
    !
    gwork1(:) = .true.
    jkm = maxval( klfs(:) )
    do jk = jkm - 1, ikb + 1, -1
       do ji = 1, iie
          if ( jk < klfs(ji) .and. zwork2(ji) > pthes(ji,jk) .and. gwork1(ji) ) then
             kdbl(ji) = jk
             gwork1(ji) = .false.
          end if
       end do
    end do
    !
    !
    !       7.     define mass flux and entr/detr. rates at lfs
    !               -------------------------------------------
    !
    do ji = 1, iie
       jk = klfs(ji)
       zwork1(ji)  = ppres(ji,jk) /                                            &
            ( xrd * zdt(ji) * ( 1._realkind + zeps * zwork1(ji) ) ) ! density
       pdmf(ji,jk) = - ( 1._realkind - ppref(ji) ) * zwork1(ji) * xpi * pcrad(ji) * pcrad(ji)  !jyj xcrad
       pdthl(ji,jk)= zwork2(ji)   ! theta_l is here actually theta_e
       zwork2(ji)  = pdmf(ji,jk)
       pddr(ji,jk) = 0._realkind
       pder(ji,jk) = - pmixf(ji) * pdmf(ji,jk)
    end do
    !
    !
    !         7.1   downdraft detrainment is assumed to occur in a layer
    !               of 60 hpa, determine top level iddt of this layer
    !               ---------------------------------------------------------
    !
    zwork1(:) = 0._realkind
    do jk = ikb + 2, jkm
       zwork1(:) = zwork1(:) + pdpres(:,jk)
       where ( jk > kdbl(:) .and. zwork1(:) <= xzpbl )
          zddt(:) = zwork1(:) 
          iddt(:) = jk
       end where
    end do
    !
    !
    !       8.     enter loop for downdraft computations. make a first guess
    !               of initial downdraft mass flux. 
    !               in the downdraft computations we use theta_es instead of 
    !               enthalpy as it allows to better take into account evaporation
    !               effects. as the downdraft detrainment rate is zero apart 
    !               from the detrainment layer, we just compute enthalpy 
    !               downdraft from theta_es in this layer.
    !               ----------------------------------------------------------
    !
    !
    zwork5(:) = 0._realkind
    !
    do jk =  jkm - 1, ikb + 1, -1
       jkp = jk + 1
       do ji = 1, iie
          if ( jk < klfs(ji) .and. jk >= iddt(ji) )  then
             pder(ji,jk)  = - zwork2(ji) * xentr * pdpres(ji,jkp) / pcrad(ji)   !jyj xcrad
             ! der and dpres are positive
             pdmf(ji,jk)  = pdmf(ji,jkp) - pder(ji,jk) 
             zpi(ji)      = ( xp00 / ppres(ji,jk) ) ** zrdocp
             zdt(ji)      = pth(ji,jk) / zpi(ji)
             zwork1(ji)   = prw(ji,jk) - prc(ji,jk) - pri(ji,jk)
             zthe(ji)     = zdt(ji) * zpi(ji) ** ( 1._realkind - 0.28_realkind * zwork1(ji) )           &
                  * exp( ( 3374.6525_realkind / zdt(ji) - 2.5403_realkind )         &
                  * zwork1(ji) * ( 1._realkind+ 0.81_realkind* zwork1(ji) ) )
             ! pdthl is here theta_es, later on in this routine this table is
             ! reskipped to enthalpy 
             pdthl(ji,jk) = ( pdthl(ji,jkp) * pdmf(ji,jkp) - zthe(ji) * pder(ji,jk)    &
                  ) / ( pdmf(ji,jk) - 1.e-7_realkind )      
             pdrw(ji,jk)  = ( pdrw(ji,jkp) * pdmf(ji,jkp) - prw(ji,jk) * pder(ji,jk)   &
                  ) / ( pdmf(ji,jk) - 1.e-7_realkind )       
          end if
          if ( jk < iddt(ji) .and. jk >= kdbl(ji) )   then
             jl = iddt(ji)
             pddr(ji,jk)  = - pdmf(ji,jl) * pdpres(ji,jkp) / zddt(ji) 
             pdmf(ji,jk)  = pdmf(ji,jkp) + pddr(ji,jk) 
             pdthl(ji,jk) = pdthl(ji,jkp)
             pdrw(ji,jk)  = pdrw(ji,jkp)
          end if
       end do
    end do
    !
    !
    !       9.     calculate total downdraft evaporation 
    !               rate for given mass flux (between dbl and iddt)
    !               -----------------------------------------------
    !
    pdtevrf(:,:) = 0._realkind
    !
    jkt = maxval( iddt(:) )
    do jk = ikb + 1, jkt
       !
       zpi(:) = ( xp00 / ppres(:,jk) ) ** zrdocp
       zdt(:) = pth(:,jk) / zpi(:)
       !
       !       9.1    determine wet bulb temperature at dbl from theta_e.
       !               the iteration algoritm is similar to that used in
       !               routine convect_condens
       !               --------------------------------------------------
       !
       do jiter = 1, 4
          call convect_satmixratio( klon, ppres(:,jk), zdt, zwork1, zlv, zls, zcph )  
          zdtp(:) = pdthl(:,jk) / ( zpi(:) ** ( 1._realkind - 0.28_realkind * zwork1(:) )         &
               * exp( ( 3374.6525_realkind / zdt(:) - 2.5403_realkind )                 &
               * zwork1(:) * ( 1._realkind + 0.81_realkind * zwork1(:) ) ) )
          zdt(:)  = 0.4_realkind * zdtp(:) + 0.6_realkind * zdt(:) ! force convergence
       end do
       !
       !
       !       9.2    sum total downdraft evaporation rate. no evaporation
       !               if actual humidity is larger than specified one.
       !               -----------------------------------------------------
       !
       !cgj
       zwork4(:) = 1._realkind-0.1_realkind/1000._realkind*(pz(:,jkt)-pz(:,jk))
       !cgj
       zwork2(:) = zwork1(:) / zdt(:) * ( xbetaw / zdt(:) - xgamw ) ! dr_sat/dt
       !cgj   zwork2(:) = zlv(:) / zcph(:) * zwork1(:) * ( 1. - xrhdbc ) /              &
       zwork2(:) = zlv(:) / zcph(:) * zwork1(:) * ( 1._realkind - zwork4(:) ) /              &
            ( 1._realkind + zlv(:) / zcph(:) * zwork2(:) ) ! temperature perturb                                                           ! due to evaporation
       zdt(:)    = zdt(:) + zwork2(:)
       !
       call convect_satmixratio( klon, ppres(:,jk), zdt, zwork3, zlv, zls, zcph )
       !
       !cgj   zwork3(:)    = zwork3(:) * xrhdbc
       zwork3(:)    = zwork3(:) * zwork4(:)
       !cgj
       zwork1(:)    = max( 0._realkind, zwork3(:) - pdrw(:,jk) ) 
       pdtevr(:)    = pdtevr(:) + zwork1(:) * pddr(:,jk) 
       pdtevrf(:,jk)= pdtevrf(:,jk) + zwork1(:) * pddr(:,jk) 
       ! compute enthalpie and humidity in the detrainment layer
       pdrw(:,jk)   = max( pdrw(:,jk), zwork3(:) ) 
       pdthl(:,jk)  = ( ( xcpd + pdrw(:,jk) * xcpv ) * zdt(:)                    &
            + ( 1._realkind + pdrw(:,jk) ) * xg * pz(:,jk) ) 
       !
    end do
    !
    !
    !      12.     if downdraft does not evaporate any water for specified 
    !               relative humidity, no downdraft is allowed
    !               ---------------------------------------------------------
    !
    zwork2(:) = 1._realkind
    where ( pdtevr(:) < 1._realkind .or. klfs(:) == ikb + 1 ) zwork2(:) = 0._realkind
    do jk = ikb, jkm
       kdbl(:)     = kdbl(:) * int( zwork2(:) ) + ( 1 - int( zwork2(:) ) ) * ikb
       klfs(:)     = klfs(:) * int( zwork2(:) ) + ( 1 - int( zwork2(:) ) ) * ikb
       pdmf(:,jk)  = pdmf(:,jk)  * zwork2(:)
       pder(:,jk)  = pder(:,jk)  * zwork2(:) 
       pddr(:,jk)  = pddr(:,jk)  * zwork2(:) 
       zwork1(:)   = real( klfs(:) - jk,realkind )         ! use this to reset thl_d
       zwork1(:)   = max( 0._realkind,min(1._realkind,zwork1(:) ) ) ! and rv_d to zero above lfs
       pdthl(:,jk) = pdthl(:,jk) * zwork2(:) * zwork1(:)
       pdrw(:,jk)  = pdrw(:,jk)  * zwork2(:) * zwork1(:)
       pdtevr(:)   = pdtevr(:)   * zwork2(:)
       pdtevrf(:,jk)= pdtevrf(:,jk) * zwork2(:)
    end do
    !
  end subroutine convect_downdraft




  subroutine convect_closure( klon, klev,                                 &
       ppres, pdpres, pz, pdxdy, plmass,           &
       pthl, pth, prw, prc, pri, otrig1,           &
       pthc, prwc, prcc, pric, pwsub,              &
       klcl, kdpl, kpbl, klfs, kctl, kml,          &
       pumf, puer, pudr, puthl, purw,              &
       purc, puri, pupr,                           &
       pdmf, pder, pddr, pdthl, pdrw,              &
       ptpr, pspr, pdtevr,                         &
       pcape, ptimec,                              &
       kftsteps,                                   &
       pdtevrf, pprlflx, pprsflx                   )

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
    !!    convect_closure_adjust
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
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!      original    26/03/96 
    !!   peter bechtold 04/10/97 change for enthalpie, r_c + r_i tendencies
    !!   dominique paquin uqam suivi des corrections de debordements
    !!                    ouranos avril 2003 correction pcp < 0 

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
    integer, dimension(klon),  intent(in) :: klfs   ! index for level of free sink
    integer, dimension(klon),  intent(in) :: klcl   ! index lifting condens. level
    integer, dimension(klon),  intent(in) :: kctl   ! index for cloud top level
    integer, dimension(klon),  intent(in) :: kdpl   ! index for departure level 
    integer, dimension(klon),  intent(in) :: kpbl   ! index for top of source layer
    integer, dimension(klon),  intent(in) :: kml    ! index for melting level
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
    real(kind=realkind), dimension(klon,klev), intent(inout):: pupr  ! updraft precipitation in
    ! flux units (kg water / s)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: puthl  ! updraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: purw   ! updraft total water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: purc   ! updraft cloud water (kg/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)  :: puri   ! updraft cloud ice   (kg/kg)
    !
    real(kind=realkind), dimension(klon,klev), intent(inout):: pdmf  ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pder  ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout):: pddr  ! downdraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in)   :: pdthl ! downdraft enthalpy (j/kg)
    real(kind=realkind), dimension(klon,klev), intent(in)   :: pdrw  ! downdraft total water (kg/kg)
    real(kind=realkind), dimension(klon),      intent(inout):: ptpr  ! total surf precipitation (kg/s)
    real(kind=realkind), dimension(klon),      intent(out)  :: pspr  ! solid surf precipitation (kg/s)
    real(kind=realkind), dimension(klon),      intent(inout):: pdtevr! donwndraft evapor. (kg/s)
    !
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pthc  ! conv. adj. grid scale theta
    real(kind=realkind), dimension(klon,klev), intent(out)  :: prwc  ! conv. adj. grid scale r_w 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: prcc  ! conv. adj. grid scale r_c 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pric  ! conv. adj. grid scale r_i 
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pwsub ! envir. compensating subsidence(pa/s)
    !
    real(kind=realkind), dimension(klon,klev), intent(inout):: pdtevrf! downdraft evaporation rate
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pprlflx! liquid precip flux
    real(kind=realkind), dimension(klon,klev), intent(out)  :: pprsflx! solid  precip flux
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
    real(kind=realkind), dimension(klon,klev) :: zdmf        !   "   downdraft mass flux
    real(kind=realkind), dimension(klon,klev) :: zder        !   "   downdraft  entrainm. rate
    real(kind=realkind), dimension(klon,klev) :: zddr        !   "   downdraft  detrainm. rate
    real(kind=realkind), dimension(klon)     :: ztpr         !   "   total precipitation
    real(kind=realkind), dimension(klon)     :: zdtevr       !   "   total downdraft evapor. 
    real(kind=realkind), dimension(klon,klev):: zprlflx      !   "   liquid precip flux
    real(kind=realkind), dimension(klon,klev):: zprsflx      !   "   solid  precip flux
    real(kind=realkind), dimension(klon)     :: zprmelt      ! melting of precipitation
    real(kind=realkind), dimension(klon)     :: zprmelto     ! non-adjusted  "
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
    real(kind=realkind), dimension(klon)     :: zlm          ! latent heat of melting
    real(kind=realkind), dimension(klon)     :: zcph         ! specific heat c_ph
    real(kind=realkind), dimension(klon)     :: zmeldpth     ! actual depth of melting layer 
    integer, dimension(klon)  :: itstep       ! fractional convective time step
    integer, dimension(klon)  :: icount       ! timestep counter 
    integer, dimension(klon)  :: ilcl         ! index lifting condens. level
    integer, dimension(klon)  :: iwork1       ! work array
    real(kind=realkind), dimension(klon)     :: zwork1, zwork2, zwork3, zwork4, zwork5
    real(kind=realkind), dimension(klon,klev):: zwork6
    logical, dimension(klon)  :: gwork1, gwork3! work arrays
    logical, dimension(klon,klev) :: gwork4    ! work array

    real(kind=realkind), dimension(klon) :: zuthed1,zuthed2   ! jyj dilute updraft theta_e 
    real(kind=realkind)                  :: zufr              ! jyj dilute updraft fraction
    !
    !
    !-------------------------------------------------------------------------------
    !
    !       0.2    initialize  local variables
    !               ----------------------------
    !
    !
    pspr(:)   = 0._realkind
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
    zdmf(:,:)  = pdmf(:,:)
    zder(:,:)  = pder(:,:)
    zddr(:,:)  = pddr(:,:)
    ztpr(:)    = ptpr(:)
    zdtevr(:)  = pdtevr(:)
    zomg(:,:)  = 0._realkind
    pwsub(:,:) = 0._realkind 
    zprmelt(:) = 0._realkind
    pprlflx(:,:) = 0._realkind
    zprlflx(:,:) = 0._realkind
    pprsflx(:,:) = 0._realkind
    zprsflx(:,:) = 0._realkind
    !
    !
    !       2.1    some preliminar computations for melting of precipitation
    !               used later in section 9 and computation of precip fluxes
    !               precipitation fluxes are updated for melting and evaporation
    !               ---------------------------------------------------------
    !
    !
    zwork1(:) = 0._realkind
    zmeldpth(:) = 0._realkind
    zwork6(:,:) = 0._realkind
    do jk = jkmax + 1, ikb + 1, -1
       ! nota: pupr is total precipitation flux, but the solid, liquid
       !       precipitation is stored in units kg/kg; therefore we compute here
       !       the solid fraction of the total precipitation flux.
       do ji = 1, iie
          zwork2(ji)    = pupr(ji,jk) / ( purc(ji,jk) + puri(ji,jk) + 1.e-8_realkind )
          zprmelt(ji)   = zprmelt(ji) + puri(ji,jk) * zwork2(ji)
          zwork1(ji)    = zwork1(ji) + purc(ji,jk) * zwork2(ji) - pdtevrf(ji,jk)
          zprlflx(ji,jk)= max( 0._realkind, zwork1(ji) )
          zprmelt(ji)   = zprmelt(ji) + min( 0._realkind, zwork1(ji) )
          zprsflx(ji,jk)= zprmelt(ji) 
          if ( kml(ji) >= jk .and. zmeldpth(ji) <= xmeldpth ) then                 
             zpi(ji)    = ( ppres(ji,jk) / xp00 ) ** zrdocp 
             zwork3(ji) = pth(ji,jk) * zpi(ji)            ! temperature estimate
             zlm(ji)    = xlstt + ( xcpv - xci ) * ( zwork3(ji) - xtt ) -       &
                  ( xlvtt + ( xcpv - xcl ) * ( zwork3(ji) - xtt ) ) ! l_s - l_v
             zcph(ji)   = xcpd + xcpv * prw(ji,jk)
             zmeldpth(ji) = zmeldpth(ji) + pdpres(ji,jk)
             zwork6(ji,jk)= zlm(ji) * ptimec(ji) / plmass(ji,jk) * pdpres(ji,jk)
             zomg(ji,jk)= 1._realkind ! at this place only used as work variable
          end if
       end do
       !
    end do
    !
    zwork2(:) = 0._realkind
    do jk = jkmax, ikb + 1, -1
       zwork1(:) = zprmelt(:) * pdpres(:,jk) / max( xmeldpth, zmeldpth(:) )
       zwork2(:) = zwork2(:) + zwork1(:) * zomg(:,jk)
       zprlflx(:,jk) = zprlflx(:,jk) + zwork2(:) 
       zprsflx(:,jk) = zprsflx(:,jk) - zwork2(:)
    end do
    where( zprsflx(:,:) < 1._realkind ) zprsflx(:,:)=0._realkind
    zprmelto(:) = zprmelt(:)
    !
    !
    !       3.     compute limits on the closure adjustment factor so that the
    !               inflow in convective drafts from a given layer can't be larger 
    !               than the mass contained in this layer initially.
    !               ---------------------------------------------------------------
    !
    zadjmax(:) = 1000._realkind
    iwork1(:) = max( ilcl(:), klfs(:) )
    jkp = minval( kdpl(:) )
    do jk = jkp, ike
       do ji = 1, iie
          if( jk > kdpl(ji) .and. jk <= iwork1(ji) ) then
             zwork1(ji)  = plmass(ji,jk) /                                      &
                  ( ( puer(ji,jk) + pder(ji,jk) + 1.e-5_realkind ) * ptimec(ji) )
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
       zthlc(:,:) = pthl(:,:) ! initialize adjusted envir. values 
       prwc(:,:)  = prw(:,:)
       prcc(:,:)  = prc(:,:)
       pric(:,:)  = pri(:,:)
       pthc(:,:)  = pth(:,:)
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
             zwork1(:)   = - ( puer(:,jkp) + pder(:,jkp) -                   &
                  pudr(:,jkp) - pddr(:,jkp) ) / plmass(:,jkp)
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
       where( gwork1(:) .and. abs( zwork1(:) ) - .01_realkind > 0._realkind )
          gwork1(:) = .false.
          ptimec(:) = 1.e-1_realkind
          ztpr(:)   = 0._realkind
          zwork5(:) = 0._realkind
       end where
       do jk = ikb, ike
          pwsub(:,jk) = pwsub(:,jk) * zwork5(:)
          zprlflx(:,jk) = zprlflx(:,jk) * zwork5(:)
          zprsflx(:,jk) = zprsflx(:,jk) * zwork5(:)
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
             !               -----------------------------------------------------------------
             !
             !
             zthlc(:,:) = zthlc(:,:) + ztimc(:,:) / plmass(:,:) * (      &
                  zthmfin(:,:) + pudr(:,:) * puthl(:,:)  +     &
                  pddr(:,:) * pdthl(:,:) - zthmfout(:,:) -     &
                  ( puer(:,:) + pder(:,:) ) * pthl(:,:)   )
             prwc(:,:)  = prwc(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrwmfin(:,:) + pudr(:,:) * purw(:,:)  +       &
                  pddr(:,:) * pdrw(:,:) - zrwmfout(:,:) -       &
                  ( puer(:,:) + pder(:,:) ) * prw(:,:)    )    
             prcc(:,:)  = prcc(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrcmfin(:,:) + pudr(:,:) * purc(:,:) - zrcmfout(:,:) -  &
                  ( puer(:,:) + pder(:,:) ) * prc(:,:)    )    
             pric(:,:)  = pric(:,:) + ztimc(:,:) / plmass(:,:) *  (      &
                  zrimfin(:,:) + pudr(:,:) * puri(:,:) - zrimfout(:,:) -  & 
                  ( puer(:,:) + pder(:,:) ) * pri(:,:)    )    
             !
             !

             !
          end where
          !
       end do ! exit the fractional time step loop
       !
       ! 
       !           9.    allow frozen precipitation to melt over a 200 mb deep layer
       !                  -----------------------------------------------------------
       !
       do jk = jkmax, ikb + 1, -1
          zthlc(:,jk) = zthlc(:,jk) -                                &
               zprmelt(:) * zwork6(:,jk) / max( xmeldpth, zmeldpth(:) )
       end do
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
       !jyj dilute updraft theta_e at lcl
       zuthed1(:) = ztheul(:)   
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
                !jyj dilute beg--------------------------------------------
                !calculate updraft dilute fraction using mass fluxes
                zufr     = (pumf(ji,jkp)-pudr(ji,jk))/(pumf(ji,jkp)-pudr(ji,jk)+puer(ji,jk)+0.000001_realkind)
                zufr     = max(0.0_realkind, min(zufr,1.0_realkind) )

                !update theta_e of updraft dilute
                zuthed2(ji) = zuthed1(ji)*zufr + zthes2(ji)*(1.0_realkind-zufr)

                !the bracket part in eq. (16.29) for calculating cape
                zwork1(ji) = ( zuthed1(ji) + zuthed2(ji) ) / ( zthes1(ji) + zthes2(ji) ) - 1._realkind

                !update theta_e for the next upper layer
                zuthed1(ji)  = zuthed2(ji)

                !jyj dilute end--------------------------------------------
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
       where ( abs(zcape(:)) < 1.e-14_realkind .and. gwork1(:) )  zadj(:) = zadj(:) * 0.5_realkind
       where ( abs(zcape(:)) > 1.e-14_realkind .and. gwork1(:) )                              &
            zadj(:) = zadj(:) * xstabc * pcape(:) / ( zwork1(:) + 1.e-8_realkind )
       zadj(:) = min( zadj(:), zadjmax(:) )  
       !
       !
       !         13.     adjust mass flux by the factor zadj to converge to
       !                  specified degree of stabilization
       !                 ----------------------------------------------------
       !
       call convect_closure_adjust( klon, klev, zadj,                     &
            pumf, zumf, puer, zuer, pudr, zudr,   &
            pdmf, zdmf, pder, zder, pddr, zddr,   &
            zprmelt, zprmelto, pdtevr, zdtevr,    &
            ptpr, ztpr,                           &
            pprlflx, zprlflx, pprsflx, zprsflx    )
       !
       !
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
    ! compute surface solid (ice) precipitation 
    pspr(:) = zprmelt(:) * ( 1._realkind - zmeldpth(:) / xmeldpth )
    pspr(:) = max( 0._realkind, pspr(:) )
    !
    !domi correction pcp < 0 
    ptpr(:) = max( ptpr(:) , pspr(:) )
    !fin correction domi
    !
    !
  end subroutine convect_closure


  subroutine convect_tstep_pref( klon, klev,                           &
       pu, pv, ppres, pz, pdxdy, klcl, kctl, &
       ptimea, ppref )

    !! routine to compute convective advection time step and precipitation 
    !!     efficiency 
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the convective
    !!      advection time step ptimec as a function of the mean ambient 
    !!      wind as well as the precipitation efficiency as a function
    !!      of wind shear and cloud base height.
    !!
    !!
    !!  method
    !!    ------
    !!     
    !!
    !!    external
    !!    --------
    !!     none
    !!
    !!
    !!    implicit arguments
    !!    ------------------
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation 
    !!      fritsch and chappell, 1980, j. atmos. sci.
    !!      kain and fritsch, 1993, meteor. monographs, vol.
    !!
    !!    author
    !!    ------
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  04/10/97
    !-------------------------------------------------------------------------------
    !
    !       0.    declarations
    !              ------------
    !
    implicit none
    !
    !       0.1   declarations of dummy arguments :
    !
    integer, intent(in)                    :: klon   ! horizontal dimension
    integer, intent(in)                    :: klev   ! vertical dimension
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres  ! pressure (pa) 
    real(kind=realkind), dimension(klon,klev), intent(in) :: pu     ! grid scale horiz. wind u 
    real(kind=realkind), dimension(klon,klev), intent(in) :: pv     ! grid scale horiz. wind v
    real(kind=realkind), dimension(klon,klev), intent(in) :: pz     ! height of model layer (m) 
    real(kind=realkind), dimension(klon),      intent(in) :: pdxdy  ! grid area (m^2)
    integer, dimension(klon),   intent(in) :: klcl   ! lifting condensation level index
    integer, dimension(klon),   intent(in) :: kctl   ! cloud top level index
    !
    real(kind=realkind), dimension(klon),      intent(out):: ptimea ! advective time period
    real(kind=realkind), dimension(klon),      intent(out):: ppref  ! precipitation efficiency 
    !
    !
    !       0.2   declarations of local variables klon
    !
    integer :: iie, ikb, ike                      ! horizontal + vertical loop bounds
    integer :: ji                                 ! horizontal loop index
    integer :: jk, jklc, jkp5, jkct               ! vertical loop index
    !
    integer, dimension(klon)  :: ip500       ! index of 500 hpa levels
    real(kind=realkind), dimension(klon)     :: zcbh        ! cloud base height 
    real(kind=realkind), dimension(klon)     :: zwork1, zwork2, zwork3  ! work arrays
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        0.3   set loop bounds
    !              ---------------
    !
    iie = klon
    ikb = 1 + jcvexb 
    ike = klev - jcvext 
    !
    !
    !       1.     determine vertical index for 500 hpa levels 
    !               ------------------------------------------
    !
    !
    ip500(:) = ikb
    do jk = ikb, ike
       where ( ppres(:,jk) >= 500.e2_realkind ) ip500(:) = jk
    end do
    !
    !
    !       2.     compute convective time step 
    !               ----------------------------
    !
    ! compute wind speed at lcl, 500 hpa, ctl

    do ji = 1, iie
       jklc = klcl(ji)
       jkp5 = ip500(ji)
       jkct = kctl(ji)
       zwork1(ji) = sqrt( pu(ji,jklc) * pu(ji,jklc) +           &
            pv(ji,jklc) * pv(ji,jklc)  ) 
       zwork2(ji) = sqrt( pu(ji,jkp5) * pu(ji,jkp5) +           &
            pv(ji,jkp5) * pv(ji,jkp5)  ) 
       zwork3(ji) = sqrt( pu(ji,jkct) * pu(ji,jkct) +           &
            pv(ji,jkct) * pv(ji,jkct)  ) 
    end do
    !
    zwork2(:) = max( 0.1_realkind, 0.5_realkind * ( zwork1(:) + zwork2(:) ) )
    !
    !correction debordement domi
    do ji = 1, iie
       ptimea(:) = sqrt( pdxdy(:) ) / zwork2(:) 
    end do
    !
    !
    !       3.     compute precipitation efficiency 
    !               -----------------------------------
    !
    !       3.1    precipitation efficiency as a function of wind shear
    !               ----------------------------------------------------
    !
    zwork2(:) = sign( 1._realkind, zwork3(:) - zwork1(:) )
    do ji = 1, iie
       jklc = klcl(ji)
       jkct = kctl(ji)
       zwork1(ji) = ( pu(ji,jkct) - pu(ji,jklc) )  *          &
            ( pu(ji,jkct) - pu(ji,jklc) )  +          &
            ( pv(ji,jkct) - pv(ji,jklc) )  *          &
            ( pv(ji,jkct) - pv(ji,jklc) )  
       zwork1(ji) = 1.e3_realkind * zwork2(ji) * sqrt( zwork1(ji) ) /  &
            max( 1.e-2_realkind, pz(ji,jkct) - pz(ji,jklc) )
    end do
    !
    ppref(:)  = 1.591_realkind + zwork1(:) * ( -.639_realkind + zwork1(:) * (        &
         9.53e-2_realkind - zwork1(:) * 4.96e-3_realkind ) ) 
    ppref(:)  = max( .4_realkind, min( ppref(:), .92_realkind ) )
    !
    !       3.2    precipitation efficiency as a function of cloud base height 
    !               ----------------------------------------------------------
    !
    do ji = 1, iie
       jklc = klcl(ji)
       zcbh(ji)   = max( 3._realkind, ( pz(ji,jklc) - pz(ji,ikb) ) * 3.281e-3_realkind ) 
    end do
    zwork1(:) = .9673_realkind + zcbh(:) * ( -.7003_realkind + zcbh(:) * ( .1622_realkind + &
         zcbh(:) *  ( -1.2570e-2_realkind + zcbh(:) * ( 4.2772e-4_realkind -  &
         zcbh(:) * 5.44e-6_realkind ) ) ) )
    zwork1(:) = max( .4_realkind, min( .92_realkind, 1._realkind/ ( 1._realkind + zwork1(:) ) ) )
    !
    !       3.3    mean precipitation efficiency is used to compute rainfall 
    !               ----------------------------------------------------------
    !
    ppref(:) = 0.5_realkind * ( ppref(:) + zwork1(:) )
    !
    !
  end subroutine convect_tstep_pref


end module modd_convpar
