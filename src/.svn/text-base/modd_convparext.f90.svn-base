module modd_convparext
  use modd_cst
  use decomp,only:realkind
  implicit none
  private
  integer,public, save :: jcvexb ! start vertical computations at
  ! 1 + jcvexb = 1 + ( kbdia - 1 )
  integer,public, save :: jcvext ! limit vertical computations to
  ! klev - jcvext = klev - ( ktdia - 1 )
  real(kind=realkind),public, save :: xtfrz1=268.16_realkind      ! begin of freezing interval
  real(kind=realkind),public, save :: xtfrz2=248.16_realkind      ! end of freezing interval
  public convect_mixing_funct,convect_chem_transport,convect_closure_adjust,convect_closure_thrvlcl,&
       convect_condens,convect_satmixratio

contains
  subroutine convect_mixing_funct( klon,                &
       pmixc, kmf, per, pdr ) 
    !     #######################################################
    !
    !!**** determine the area under the distribution function
    !!     kmf = 1 : gaussian  kmf = 2 : triangular distribution function
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the entrainment and
    !!      detrainment rate by evaluating the are under the distribution 
    !!      function. the integration interval is limited by the critical
    !!      mixed fraction pmixc
    !!   
    !!
    !!
    !!**  method
    !!    ------
    !!      use handbook of mathemat. functions by abramowitz and stegun, 1968
    !!      
    !!     
    !!
    !!    external
    !!    --------
    !!      none
    !!     
    !!
    !!    implicit arguments
    !!    ------------------
    !!      none
    !!
    !!
    !!    reference
    !!    ---------
    !!
    !!      book2 of documentation ( routine mixing_funct)
    !!      abramovitz and stegun (1968), handbook of math. functions 
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
    !*       0.    declarations
    !              ------------
    !
    !
    implicit none
    !
    !*       0.1   declarations of dummy arguments :
    !
    integer,               intent(in) :: klon   ! horizontal dimension
    integer,               intent(in) :: kmf    ! switch for dist. function
    real(kind=realkind), dimension(klon), intent(in) :: pmixc  ! critical mixed fraction
    !
    real(kind=realkind), dimension(klon), intent(out):: per    ! normalized entrainment rate
    real(kind=realkind), dimension(klon), intent(out):: pdr    ! normalized detrainment rate
    !
    !*       0.2   declarations of local variables :
    !
    real(kind=realkind)    :: zsigma = 0.166666667_realkind                   ! standard deviation 
    real(kind=realkind)    :: zfe    = 4.931813949_realkind                   ! integral normalization 
    real(kind=realkind)    :: zsqrtp = 2.506628_realkind,  zp  = 0.33267_realkind      ! constants
    real(kind=realkind)    :: za1    = 0.4361836_realkind, za2 =-0.1201676_realkind    ! constants
    real(kind=realkind)    :: za3    = 0.9372980_realkind, zt1 = 0.500498_realkind     ! constants
    real(kind=realkind)    :: ze45   = 0.01111_realkind                       ! constant
    !
    real(kind=realkind), dimension(klon) :: zx, zy, zw1, zw2         ! work variables
    real(kind=realkind)    :: zw11
    !
    !
    !-------------------------------------------------------------------------------
    !
    !       1.     use gaussian function for kmf=1
    !              -------------------------------
    !
    if( kmf == 1 ) then 
       ! zx(:)  = ( pmixc(:) - 0.5 ) / zsigma
       zx(:)  = 6._realkind * pmixc(:) - 3._realkind
       zw1(:) = 1._realkind / ( 1._realkind+ zp * abs ( zx(:) ) )
       zy(:)  = exp( -0.5_realkind * zx(:) * zx(:) )
       zw2(:) = za1 * zw1(:) + za2 * zw1(:) * zw1(:) +                   &
            za3 * zw1(:) * zw1(:) * zw1(:)
       zw11   = za1 * zt1 + za2 * zt1 * zt1 + za3 * zt1 * zt1 * zt1
    endif
    !
    where ( kmf == 1 .and. zx(:) >= 0._realkind )
       per(:) = zsigma * ( 0.5_realkind * ( zsqrtp - ze45 * zw11                 &
            - zy(:) * zw2(:) ) + zsigma * ( ze45 - zy(:) ) )        &
            - 0.5_realkind * ze45 * pmixc(:) * pmixc(:)
       pdr(:) = zsigma*( 0.5_realkind * ( zy(:) * zw2(:) - ze45 * zw11   )       &
            + zsigma * ( ze45 - zy(:) ) )                           &
            - ze45 * ( 0.5_realkind + 0.5_realkind * pmixc(:) * pmixc(:) - pmixc(:) )
    end where
    where ( kmf == 1 .and. zx(:) < 0._realkind ) 
       per(:) = zsigma*( 0.5_realkind * ( zy(:) * zw2(:) - ze45 * zw11   )       &
            + zsigma * ( ze45 - zy(:) ) )                           &
            - 0.5_realkind * ze45 * pmixc(:) * pmixc(:)
       pdr(:) = zsigma * ( 0.5_realkind * ( zsqrtp - ze45 * zw11 - zy(:)         &
            * zw2(:) ) + zsigma * ( ze45 - zy(:) ) )                &
            - ze45 * ( 0.5_realkind + 0.5_realkind * pmixc(:) * pmixc(:) - pmixc(:) )
    end where
    !
    per(:) = per(:) * zfe
    pdr(:) = pdr(:) * zfe
    !
    !
    !       2.     use triangular function kmf=2
    !              -------------------------------
    !
    !     not yet released
    !
    !
  end subroutine convect_mixing_funct


  subroutine convect_chem_transport( klon, klev, kch, pch1, pch1c,       &
       kdpl, kpbl, klcl, kctl, klfs, kdbl, &
       pumf, puer, pudr, pdmf, pder, pddr, &
       ptimec, pdxdy, pmixf, plmass, pwsub,&
       kftsteps )
    !    #######################################################################
    !
    !!**** compute  modified chemical tracer values due to convective event
    !!
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the final adjusted
    !!      environmental values of the chemical tracers
    !!      the final convective tendencies can then be evaluated in the main
    !!      routine deep_convect by (pch1c-pch1)/ptimec
    !!
    !!
    !!**  method
    !!    ------
    !!      identical to the computation of the conservative variables in the
    !!      main deep convection code
    !!
    !!    external
    !!    --------
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_cst
    !!          xg                 ! gravity constant
    !!
    !!     module modd_convparext
    !!          jcvexb, jcvext     ! extra levels on the vertical boundaries
    !!
    !!    author
    !!    ------
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!
    !!      original    11/12/97
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    declarations
    !              ------------
    !

    !
    implicit none
    !
    !*       0.1   declarations of dummy arguments :
    !
    integer,                intent(in) :: klon     ! horizontal dimension
    integer,                intent(in) :: klev     ! vertical dimension
    integer,                intent(in) :: kch      ! number of passive tracers
    !
    real(kind=realkind),dimension(klon,klev,kch),intent(in) :: pch1 ! grid scale tracer concentr.
    real(kind=realkind),dimension(klon,klev,kch),intent(out):: pch1c! conv adjusted tracer concntr.
    !
    integer, dimension(klon), intent(in) :: kdpl   ! index for departure level
    integer, dimension(klon), intent(in) :: kpbl   ! index for top of source layer
    integer, dimension(klon), intent(in) :: klcl   ! index lifting condens. level
    integer, dimension(klon), intent(in) :: kctl   ! index for cloud top level
    integer, dimension(klon), intent(in) :: klfs   ! index for level of free sink
    integer, dimension(klon), intent(in) :: kdbl   ! index for downdraft base level
    !
    real(kind=realkind), dimension(klon,klev), intent(in) :: pumf ! updraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: puer ! updraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pudr ! updraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pdmf ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pder ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(in) :: pddr ! downdraft detrainment (kg/s)
    !
    real(kind=realkind), dimension(klon),     intent(in) :: ptimec! convection time step
    real(kind=realkind), dimension(klon),     intent(in) :: pdxdy ! grid area (m^2)
    real(kind=realkind), dimension(klon),     intent(in) :: pmixf ! mixed fraction at lfs
    real(kind=realkind), dimension(klon,klev),intent(in) :: plmass! mass of model layer (kg)
    real(kind=realkind), dimension(klon,klev),intent(in) :: pwsub ! envir. compensating subsidence(pa/s)
    integer,                intent(in) :: kftsteps  ! maximum fractional time steps
    !
    !
    !*       0.2   declarations of local variables :
    !
    integer :: inch1          ! number of chemical tracers
    integer :: iie, ikb, ike  ! horizontal + vertical loop bounds
    integer :: iks            ! vertical dimension
    integer :: ji             ! horizontal loop index
    integer :: jk, jkp        ! vertical loop index
    integer :: jn             ! chemical tracer loop index
    integer :: jstep          ! fractional time loop index
    integer :: jklc, jkld, jklp, jkmax ! loop index for levels
    !
    real(kind=realkind), dimension(klon,klev)     :: zomg ! compensat. subsidence (pa/s)
    real(kind=realkind), dimension(klon,klev,kch) :: zuch1, zdch1 ! updraft/downdraft values
    real(kind=realkind), dimension(klon)          :: ztimec  ! fractional convective time step
    real(kind=realkind), dimension(klon,klev)     :: ztimc! 2d work array for ztimec
    real(kind=realkind), dimension(klon,klev,kch) :: zch1mfin, zch1mfout
    ! work arrays for environm. compensat. mass
    real(kind=realkind), dimension(klon,kch)      :: zwork1, zwork2, zwork3
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.3   compute loop bounds
    !              -------------------
    !
    inch1  = kch
    iie    = klon
    ikb    = 1 + jcvexb 
    iks    = klev
    ike    = klev - jcvext 
    jkmax  = maxval( kctl(:) )
    !
    !
    !*      2.      updraft computations
    !               --------------------
    !
    zuch1(:,:,:) = 0._realkind
    !
    !*      2.1     initialization  at lcl
    !               ----------------------------------
    !
    do ji = 1, iie
       jklc = klcl(ji)
       jkld = kdpl(ji)
       jklp = kpbl(ji)
       zwork1(ji,:) = .5_realkind * ( pch1(ji,jkld,:) + pch1(ji,jklp,:) )
    end do
    !
    !*      2.2     final updraft loop
    !               ------------------
    !
    do jk = minval( kdpl(:) ), jkmax
       jkp = jk + 1
       !
       do jn = 1, inch1
          do ji = 1, iie
             if ( kdpl(ji) <= jk .and. klcl(ji) > jk )                             &
                  zuch1(ji,jk,jn) = zwork1(ji,jn)
             !
             if ( klcl(ji) - 1 <= jk .and. kctl(ji) > jk ) then
                zuch1(ji,jkp,jn) = zuch1(ji,jk,jn) 
                !if you have reactive i.e. non-passive tracers
                ! update their values here and add the corresponding
                ! sink term in the following equation
                zuch1(ji,jkp,jn) = ( pumf(ji,jk) * zuch1(ji,jk,jn) +              &
                     puer(ji,jkp) * pch1(ji,jk,jn) )  /           & 
                     ( pumf(ji,jkp) + pudr(ji,jkp) + 1.e-7_realkind )
             end if
          end do
       end do
       !
    end do
    !
    !*      3.      downdraft computations
    !               ----------------------
    !
    zdch1(:,:,:) = 0._realkind
    !
    !*      3.1     initialization at the lfs
    !               -------------------------
    !
    zwork1(:,:) = spread( pmixf(:), dim=2, ncopies=inch1 )
    do ji = 1, iie
       jk = klfs(ji)
       zdch1(ji,jk,:) = zwork1(ji,:) * pch1(ji,jk,:) +                          &
            ( 1._realkind - zwork1(ji,:) ) * zuch1(ji,jk,:)
    end do
    !
    !*      3.2     final downdraft loop
    !               --------------------
    !
    do jk = maxval( klfs(:) ), ikb + 1, -1
       jkp = jk - 1
       do jn = 1, inch1
          do ji = 1, iie
             if ( jk <= klfs(ji) .and. jkp >= kdbl(ji) ) then
                zdch1(ji,jkp,jn) = ( zdch1(ji,jk,jn) * pdmf(ji,jk) -              &
                     pch1(ji,jk,jn) *  pder(ji,jkp) ) /           &
                     ( pdmf(ji,jkp) - pddr(ji,jkp) - 1.e-7_realkind ) 
             end if
          end do
       end do
    end do
    !
    !							   
    !*      4.      final closure (environmental) computations
    !               ------------------------------------------
    !
    pch1c(:,ikb:ike,:) = pch1(:,ikb:ike,:) ! initialize adjusted envir. values
    !
    do jk = ikb, ike
       zomg(:,jk) = pwsub(:,jk) * pdxdy(:) / xg ! environmental subsidence
    end do
    !
    ztimec(:) = ptimec(:) / real( kftsteps,realkind ) ! adjust  fractional time step
    ! to be an integer multiple of ptimec
    where ( ptimec(:) < 1._realkind ) ztimec(:) = 0._realkind
    ztimc(:,:)= spread( ztimec(:), dim=2, ncopies=iks )
    !
    zch1mfin(:,:,:)   = 0._realkind
    zch1mfout(:,:,:)  = 0._realkind
    !
    do jstep = 1, kftsteps ! enter the fractional time step loop
       !
       do jk = ikb + 1, jkmax
          jkp = max( ikb + 1, jk - 1 )
	  zwork3(:,:) = spread( zomg(:,jk), dim=2, ncopies=inch1 )
          zwork1(:,:) = sign( 1._realkind, zwork3(:,:) )
          zwork2(:,:) = 0.5_realkind * ( 1._realkind + zwork1(:,:) )
          zwork1(:,:) = 0.5_realkind * ( 1._realkind - zwork1(:,:) )
          zch1mfin(:,jk,:)  = - zwork3(:,:) * pch1c(:,jkp,:) * zwork1(:,:)
          zch1mfout(:,jk,:) =   zwork3(:,:) * pch1c(:,jk,:)  * zwork2(:,:)
          zch1mfin(:,jkp,:) = zch1mfin(:,jkp,:) + zch1mfout(:,jk,:) * zwork2(:,:)
          zch1mfout(:,jkp,:)= zch1mfout(:,jkp,:) + zch1mfin(:,jk,:) * zwork1(:,:)
       end do
       !
       do jn = 1, inch1
          do jk = ikb + 1, jkmax
             pch1c(:,jk,jn) = pch1c(:,jk,jn) + ztimc(:,jk) / plmass(:,jk) *  (    &
                  zch1mfin(:,jk,jn) + pudr(:,jk) * zuch1(:,jk,jn) +       &
                  pddr(:,jk) * zdch1(:,jk,jn) - zch1mfout(:,jk,jn) -      &
                  ( puer(:,jk) + pder(:,jk) ) * pch1(:,jk,jn)    )
             pch1c(:,jk,jn) = max( 0._realkind, pch1c(:,jk,jn) )
          end do
       end do
       !
    end do ! exit the fractional time step loop
    !
    !
  end subroutine convect_chem_transport


  subroutine convect_closure_adjust( klon, klev, padj,                      &
       pumf, pzumf, puer, pzuer, pudr, pzudr, &
       pdmf, pzdmf, pder, pzder, pddr, pzddr, &
       pprmelt, pzprmelt, pdtevr, pzdtevr,    &
       ptpr, pztpr,                           &
       pprlflx, pzprlfl, pprsflx, pzprsfl     )

    !    #########################################################################
    !
    !!**** uses closure adjustment factor to adjust mass flux and to modify
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
    !!**  method
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
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!      original    26/03/96 
    !!   last modified  04/10/97
    !-------------------------------------------------------------------------------
    !
    !*       0.    declarations
    !              ------------
    !

    !
    implicit none
    !
    !*       0.1   declarations of dummy arguments :
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
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pdmf  ! downdraft mass flux (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzdmf ! initial value of  "
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pder  ! downdraft entrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzder ! initial value of  "
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pddr  ! downdraft detrainment (kg/s)
    real(kind=realkind), dimension(klon,klev), intent(inout) :: pzddr ! initial value of  "
    real(kind=realkind), dimension(klon),   intent(inout):: ptpr     ! total precipitation (kg/s)
    real(kind=realkind), dimension(klon),   intent(inout):: pztpr    ! initial value of "
    real(kind=realkind), dimension(klon),   intent(inout):: pdtevr   ! donwndraft evapor. (kg/s)
    real(kind=realkind), dimension(klon),   intent(inout):: pzdtevr  ! initial value of " 
    real(kind=realkind), dimension(klon),   intent(inout):: pprmelt  ! melting of precipitation
    real(kind=realkind), dimension(klon),   intent(inout):: pzprmelt ! initial value of " 
    real(kind=realkind), dimension(klon,klev),intent(inout)  :: pprlflx! liquid precip flux
    real(kind=realkind), dimension(klon,klev),intent(inout)  :: pzprlfl! initial value "
    real(kind=realkind), dimension(klon,klev),intent(inout)  :: pprsflx! solid  precip flux
    real(kind=realkind), dimension(klon,klev),intent(inout)  :: pzprsfl! initial value "
    !
    !
    !*       0.2   declarations of local variables :
    !
    integer :: iie, ikb, ike                 ! horiz. + vert. loop bounds
    integer :: jk                            ! vertical loop index
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.3   compute loop bounds
    !              -------------------
    !
    iie  = klon
    ikb  = 1 + jcvexb 
    ike  = klev - jcvext 
    !
    !
    !*       1.     adjust mass flux by the factor padj to converge to
    !               specified degree of stabilization
    !               ----------------------------------------------------
    !
    pprmelt(:)  = pzprmelt(:)   * padj(:)
    pdtevr(:)   = pzdtevr(:)    * padj(:)
    ptpr(:)     = pztpr(:)      * padj(:)
    !
    do jk = ikb + 1, ike
       pumf(:,jk)  = pzumf(:,jk)   * padj(:)
       puer(:,jk)  = pzuer(:,jk)   * padj(:)
       pudr(:,jk)  = pzudr(:,jk)   * padj(:)
       pdmf(:,jk)  = pzdmf(:,jk)   * padj(:)
       pder(:,jk)  = pzder(:,jk)   * padj(:)
       pddr(:,jk)  = pzddr(:,jk)   * padj(:)
       pprlflx(:,jk) = pzprlfl(:,jk) * padj(:)
       pprsflx(:,jk) = pzprsfl(:,jk) * padj(:)
    end do
    !
  end subroutine convect_closure_adjust


  subroutine convect_closure_thrvlcl( klon, klev,                         &
       ppres, pth, prv, pz, owork1,        &
       pthlcl, prvlcl, pzlcl, ptlcl, ptelcl,&
       klcl, kdpl, kpbl )
    !     ######################################################################
    !
    !!**** determine thermodynamic properties at new lcl
    !!
    !!    purpose
    !!    -------
    !!      the purpose of this routine is to determine the thermodynamic
    !!      properties at the new lifting condensation level lcl
    !!   
    !!
    !!
    !!**  method
    !!    ------
    !!    see convect_trigger_funct
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
    !!          xzlcl              ! lowest allowed pressure difference between
    !!                             ! surface and lcl
    !!          xzpbl              ! minimum mixed layer depth to sustain convection
    !!          xwtrig             ! constant in vertical velocity trigger
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
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  04/10/97
    !-------------------------------------------------------------------------------
    !
    !*       0.    declarations
    !              ------------
    !

    implicit none
    !
    !*       0.1   declarations of dummy arguments :
    !
    integer,                    intent(in) :: klon  ! horizontal dimension
    integer,                    intent(in) :: klev  ! vertical dimension
    real(kind=realkind), dimension(klon,klev), intent(in) :: pth   ! theta
    real(kind=realkind), dimension(klon,klev), intent(in) :: prv   ! vapor mixing ratio 
    real(kind=realkind), dimension(klon,klev), intent(in) :: ppres ! pressure
    real(kind=realkind), dimension(klon,klev), intent(in) :: pz    ! height of grid point (m)
    integer, dimension(klon),   intent(in) :: kdpl  ! contains vert. index of dpl
    integer, dimension(klon),   intent(in) :: kpbl  ! " vert. index of source layer top
    logical, dimension(klon),   intent(in) :: owork1! logical mask 
    !
    real(kind=realkind), dimension(klon),     intent(out):: pthlcl ! theta at lcl
    real(kind=realkind), dimension(klon),     intent(out):: prvlcl ! vapor mixing ratio at  lcl
    real(kind=realkind), dimension(klon),     intent(out):: pzlcl  ! height at lcl (m)
    real(kind=realkind), dimension(klon),     intent(out):: ptlcl  ! temperature at lcl (m)
    real(kind=realkind), dimension(klon),     intent(out):: ptelcl ! environm. temp. at lcl (k)
    integer, dimension(klon),  intent(out):: klcl   ! contains vert. index of lcl
    !
    !*       0.2   declarations of local variables :
    !
    integer :: jk, jkm, jkmin, jkmax      ! vertical loop index
    integer :: ji                         ! horizontal loop index 
    integer :: iie, ikb, ike              ! horizontal + vertical loop bounds
    real(kind=realkind)    :: zeps, zepsa    ! r_d / r_v, r_v / r_d 
    real(kind=realkind)    :: zcpord, zrdocp ! c_pd / r_d, r_d / c_pd
    !
    real(kind=realkind), dimension(klon) :: zplcl    ! pressure at lcl
    real(kind=realkind), dimension(klon) :: ztmix    ! mixed layer temperature
    real(kind=realkind), dimension(klon) :: zevmix   ! mixed layer water vapor pressure 
    real(kind=realkind), dimension(klon) :: zdpthmix, zpresmix ! mixed layer depth and pressure
    real(kind=realkind), dimension(klon) :: zlv, zcph! specific heats of vaporisation, dry air
    real(kind=realkind), dimension(klon) :: zdp      ! pressure between lcl and model layer
    real(kind=realkind), dimension(klon) :: zwork1, zwork2     ! work arrays
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.3    compute array bounds
    !               --------------------
    !
    iie = klon
    ikb = 1 + jcvexb 
    ike = klev - jcvext 
    !
    !
    !*       1.     initialize local variables
    !               --------------------------
    !
    zeps      = xrd / xrv
    zepsa     = xrv / xrd 
    zcpord    = xcpd / xrd
    zrdocp    = xrd / xcpd
    !
    zdpthmix(:) = 0._realkind
    zpresmix(:) = 0._realkind
    pthlcl(:)   = 300._realkind
    ptlcl(:)    = 300._realkind
    ptelcl(:)   = 300._realkind
    prvlcl(:)   = 0._realkind
    pzlcl(:)    = pz(:,ikb)
    ztmix(:)    = 230._realkind
    zplcl(:)    = 1.e4_realkind 
    klcl(:)     = ikb + 1
    !
    !
    !*       2.     construct a mixed layer as in trigger_funct
    !               -------------------------------------------
    !
    jkmax = maxval( kpbl(:) )
    jkmin = minval( kdpl(:) )
    do jk = ikb + 1, jkmax
       jkm = jk + 1
       do ji = 1, iie
          if ( jk >= kdpl(ji) .and. jk <= kpbl(ji) ) then
             !           
             zwork1(ji)   = ppres(ji,jk) - ppres(ji,jkm)
             zdpthmix(ji) = zdpthmix(ji) + zwork1(ji)
             zpresmix(ji) = zpresmix(ji) + ppres(ji,jk) * zwork1(ji)
             pthlcl(ji)   = pthlcl(ji)   + pth(ji,jk)   * zwork1(ji)
             prvlcl(ji)   = prvlcl(ji)   + prv(ji,jk)   * zwork1(ji)
             !
          end if
       end do
    end do
    !
    !
    where ( owork1(:) )
       !
       zpresmix(:) = zpresmix(:) / zdpthmix(:)
       pthlcl(:)   = pthlcl(:)   / zdpthmix(:)
       prvlcl(:)   = prvlcl(:)   / zdpthmix(:)
       !
       !*       3.1    use an empirical direct solution ( bolton formula )
       !               to determine temperature and pressure at lcl.
       !               nota: the adiabatic saturation temperature is not
       !                     equal to the dewpoint temperature
       !               --------------------------------------------------
       !
       !
       ztmix(:)  = pthlcl(:) * ( zpresmix(:) / xp00 ) ** zrdocp
       zevmix(:) = prvlcl(:) * zpresmix(:) / ( prvlcl(:) + zeps )
       zevmix(:) = max( 1.e-8_realkind, zevmix(:) )
       zwork1(:) = log( zevmix(:) / 613.3_realkind )
       ! dewpoint temperature
       zwork1(:) = ( 4780.8_realkind - 32.19_realkind * zwork1(:) ) / ( 17.502_realkind - zwork1(:) ) 
       ! adiabatic saturation temperature
       ptlcl(:)  = zwork1(:) - ( .212_realkind + 1.571e-3_realkind * ( zwork1(:) - xtt )      &
            - 4.36e-4_realkind * ( ztmix(:) - xtt ) ) * ( ztmix(:) - zwork1(:) )
       ptlcl(:)  = min( ptlcl(:), ztmix(:) )
       zplcl(:)  = xp00 * ( ptlcl(:) / pthlcl(:) ) ** zcpord
       !
    end where
    !
    zplcl(:) = min( 2.e5_realkind, max( 10._realkind, zplcl(:) ) ) ! bound to avoid overflow
    !
    !
    !*       3.2    correct ptlcl in order to be completely consistent
    !               with mnh saturation formula
    !               --------------------------------------------------
    !
    call convect_satmixratio( klon, zplcl, ptlcl, zwork1, zlv, zwork2, zcph )
    where( owork1(:) )
       zwork2(:) = zwork1(:) / ptlcl(:) * ( xbetaw / ptlcl(:) - xgamw ) ! dr_sat/dt
       zwork2(:) = ( zwork1(:) - prvlcl(:) ) /                              &
            ( 1._realkind + zlv(:) / zcph(:) * zwork2(:) ) 
       ptlcl(:)  = ptlcl(:) - zlv(:) / zcph(:) * zwork2(:)
       !
    end where
    !
    !
    !*       3.3    if prvlcl is oversaturated set humidity and temperature
    !               to saturation values.
    !               -------------------------------------------------------
    !
    call convect_satmixratio( klon, zpresmix, ztmix, zwork1, zlv, zwork2, zcph )
    where( owork1(:) .and. prvlcl(:) > zwork1(:) )
       zwork2(:) = zwork1(:) / ztmix(:) * ( xbetaw / ztmix(:) - xgamw ) ! dr_sat/dt
       zwork2(:) = ( zwork1(:) - prvlcl(:) ) /                              &
            ( 1._realkind + zlv(:) / zcph(:) * zwork2(:) )
       ptlcl(:)  = ztmix(:) + zlv(:) / zcph(:) * zwork2(:)
       prvlcl(:) = prvlcl(:) - zwork2(:)
       zplcl(:)  = zpresmix(:)
       pthlcl(:) = ptlcl(:) * ( xp00 / zplcl(:) ) ** zrdocp
    end where
    !
    !
    !*        4.1   determine  vertical loop index at the lcl 
    !               -----------------------------------------
    !
    do jk = jkmin, ike - 1
       do ji = 1, iie
          if ( zplcl(ji) <= ppres(ji,jk) .and. owork1(ji) ) then
             klcl(ji)  = jk + 1
             pzlcl(ji) = pz(ji,jk+1)
          end if
       end do
    end do
    !
    !
    !*        4.2   estimate height and environmental temperature at lcl
    !               ----------------------------------------------------
    !
    do ji = 1, iie
       jk   = klcl(ji)
       jkm  = jk - 1
       zdp(ji)     = log( zplcl(ji) / ppres(ji,jkm) ) /                     &
            log( ppres(ji,jk) / ppres(ji,jkm) )
       zwork1(ji)  = pth(ji,jk)  * ( ppres(ji,jk) / xp00 ) ** zrdocp
       zwork2(ji)  = pth(ji,jkm) * ( ppres(ji,jkm) / xp00 ) ** zrdocp
       zwork1(ji)  = zwork2(ji) + ( zwork1(ji) - zwork2(ji) ) * zdp(ji) 
       ! we compute the precise value of the lcl
       ! the precise height is between the levels klcl and klcl-1.
       zwork2(ji) = pz(ji,jkm) + ( pz(ji,jk) - pz(ji,jkm) ) * zdp(ji)
    end do
    where( owork1(:) )
       ptelcl(:) = zwork1(:)
       pzlcl(:)  = zwork2(:)
    end where
    !        
    !
    !
  end subroutine convect_closure_thrvlcl


  subroutine convect_condens( klon,                                           &
       kice, ppres, pthl, prw, prco, prio, pz, owork1, &
       pt, pew, prc, pri, plv, pls, pcph   )
    !     ###########################################################################
    !
    !!**** compute temperature cloud and ice water content from enthalpy and r_w 
    !!
    !!
    !!    purpose
    !!    -------
    !!     the purpose of this routine is to determine cloud condensate
    !!     and to return values for l_v, l_s and c_ph
    !!
    !!
    !!**  method
    !!    ------
    !!     condensate is extracted iteratively 
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
    !!      module modd_cst
    !!          xg                   ! gravity constant
    !!          xalpw, xbetaw, xgamw ! constants for water saturation pressure
    !!          xalpi, xbetai, xgami ! constants for ice saturation pressure
    !!          xp00                 ! reference pressure
    !!          xrd, xrv             ! gaz  constants for dry air and water vapor
    !!          xcpd, xcpv           ! specific heat for dry air and water vapor
    !!          xcl, xci             ! specific heat for liquid water and ice
    !!          xtt                  ! triple point temperature
    !!          xlvtt, xlstt         ! vaporization, sublimation heat constant
    !!
    !!    implicit arguments
    !!    ------------------
    !!      module modd_convpar
    !!          xtfrz1               ! begin of freezing interval
    !!          xtfrz2               ! end of freezing interval
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_condens)
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
    !*       0.    declarations
    !              ------------
    !
    !
    !
    implicit none
    !
    !*       0.1   declarations of dummy arguments :
    !
    integer, intent(in)                :: klon    ! horizontal loop index
    integer, intent(in)                :: kice    ! flag for ice ( 1 = yes,
    !                0 = no ice )
    real(kind=realkind), dimension(klon),   intent(in) :: ppres  ! pressure
    real(kind=realkind), dimension(klon),   intent(in) :: pthl   ! enthalpy (j/kg)
    real(kind=realkind), dimension(klon),   intent(in) :: prw    ! total water mixing ratio  
    real(kind=realkind), dimension(klon),   intent(in) :: prco   ! cloud water estimate (kg/kg)
    real(kind=realkind), dimension(klon),   intent(in) :: prio   ! cloud ice   estimate (kg/kg)
    real(kind=realkind), dimension(klon),   intent(in) :: pz     ! level height (m)
    logical, dimension(klon),intent(in) :: owork1 ! logical mask         
    !
    !
    real(kind=realkind), dimension(klon),   intent(out):: pt     ! temperature   
    real(kind=realkind), dimension(klon),   intent(out):: prc    ! cloud water mixing ratio(kg/kg)
    real(kind=realkind), dimension(klon),   intent(out):: pri    ! cloud ice mixing ratio  (kg/kg)
    real(kind=realkind), dimension(klon),   intent(out):: plv    ! latent heat l_v    
    real(kind=realkind), dimension(klon),   intent(out):: pls    ! latent heat l_s  
    real(kind=realkind), dimension(klon),   intent(out):: pcph   ! specific heat c_ph   
    real(kind=realkind), dimension(klon),   intent(out):: pew    ! water saturation mixing ratio  
    !
    !*       0.2   declarations of local variables klon
    !
    integer :: jiter          ! iteration index
    real(kind=realkind)    :: zeps, zepsa    ! r_d / r_v, 1 / zeps
    real(kind=realkind)    :: zcvocd         ! xcpv / xcpd
    real(kind=realkind)    :: zrdocp         ! r_d / c_pd
    !
    real(kind=realkind), dimension(klon)    :: zei           ! ice saturation mixing ratio
    real(kind=realkind), dimension(klon)    :: zwork1, zwork2, zwork3, zt ! work arrays
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     initialize temperature and exner function
    !               -----------------------------------------
    !
    zrdocp      = xrd / xcpd  
    zeps        = xrd / xrv
    zepsa       = 1._realkind / zeps
    zcvocd      = xcpv / xcpd
    !
    !
    ! make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level
    !
    !! note that the definition of zcph is not the same as used in
    !! routine convect_satmixratio
    pcph(:)   = xcpd + xcpv * prw(:)
    zwork1(:) = ( 1._realkind + prw(:) ) * xg * pz(:)
    pt(:)     = ( pthl(:) + prco(:) * xlvtt + prio(:) * xlstt - zwork1(:) )   &
         / pcph(:)
    pt(:)     = max(180._realkind, min( 330._realkind, pt(:) ) ) ! set overflow bounds in
    ! case that pthl=0     
    !
    !
    !*       2.     enter the iteration loop
    !               ------------------------
    !    
    do jiter = 1,6
       pew(:) = exp( xalpw - xbetaw / pt(:) - xgamw * log( pt(:) ) )
       zei(:) = exp( xalpi - xbetai / pt(:) - xgami * log( pt(:) ) )
       pew(:) = zeps * pew(:) / ( ppres(:) - pew(:) )
       zei(:) = zeps * zei(:) / ( ppres(:) - zei(:) )    
       !
       plv(:)    = xlvtt + ( xcpv - xcl ) * ( pt(:) - xtt ) ! compute l_v
       pls(:)    = xlstt + ( xcpv - xci ) * ( pt(:) - xtt ) ! compute l_i
       !    
       zwork2(:) = ( xtfrz1 - pt(:) ) / ( xtfrz1 - xtfrz2 ) ! freezing interval
       zwork2(:) = max( 0._realkind, min(1._realkind, zwork2(:) ) ) * real( kice,realkind )
       zwork3(:) = ( 1._realkind - zwork2(:) ) * pew(:) + zwork2(:) * zei(:)
       prc(:)    = max( 0._realkind, ( 1._realkind - zwork2(:) ) * ( prw(:) - zwork3(:) ) )
       pri(:)    = max( 0._realkind,  zwork2(:) * ( prw(:) - zwork3(:) ) )
       zt(:)     = ( pthl(:) + prc(:) * plv(:) + pri(:) * pls(:) - zwork1(:) )   &
            / pcph(:)
       pt(:) = pt(:) + ( zt(:) - pt(:) ) * 0.4_realkind  ! force convergence
       pt(:) = max( 175._realkind, min( 330._realkind, pt(:) ) )
    end do
    !
    !
  end subroutine convect_condens

  subroutine convect_satmixratio( klon,                          &
       ppres, pt, pew, plv, pls, pcph )      
    !     ################################################################
    !
    !!**** compute vapor saturation mixing ratio over liquid water
    !!
    !!
    !!    pdrpose
    !!    -------
    !!     the purpose of this routine is to determine saturation mixing ratio
    !!     and to return values for l_v l_s and c_ph
    !!
    !!
    !!**  method
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
    !!      module modd_cst
    !!          xalpw, xbetaw, xgamw ! constants for water saturation pressure
    !!          xrd, xrv             ! gaz  constants for dry air and water vapor
    !!          xcpd, xcpv           ! specific heat for dry air and water vapor
    !!          xcl, xci             ! specific heat for liquid water and ice
    !!          xtt                  ! triple point temperature
    !!          xlvtt, xlstt         ! vaporization, sublimation heat constant
    !!
    !!
    !!    reference
    !!    ---------
    !!
    !!      book1,2 of documentation ( routine convect_satmixratio)
    !!
    !!    author
    !!    ------
    !!      p. bechtold       * laboratoire d'aerologie *
    !!
    !!    modifications
    !!    -------------
    !!      original    07/11/95 
    !!   last modified  04/10/97
    !------------------------- ------------------------------------------------------
    !
    !*       0.    declarations
    !              ------------
    !
    !
    implicit none
    !
    !*       0.1   declarations of dummy arguments :
    !
    !
    integer,                intent(in) :: klon    ! horizontal loop index
    real(kind=realkind), dimension(klon),  intent(in) :: ppres   ! pressure
    real(kind=realkind), dimension(klon),  intent(in) :: pt      ! temperature   
    !
    real(kind=realkind), dimension(klon),  intent(out):: pew     ! vapor saturation mixing ratio
    real(kind=realkind), dimension(klon),  intent(out):: plv     ! latent heat l_v    
    real(kind=realkind), dimension(klon),  intent(out):: pls     ! latent heat l_s  
    real(kind=realkind), dimension(klon),  intent(out):: pcph    ! specific heat c_ph   
    !
    !*       0.2   declarations of local variables :
    !
    real(kind=realkind), dimension(klon)              :: zt      ! temperature   
    real(kind=realkind)    :: zeps           ! r_d / r_v
    !
    !
    !-------------------------------------------------------------------------------
    !
    zeps      = xrd / xrv
    !
    zt(:)     = min( 400._realkind, max( pt(:), 10._realkind ) ) ! overflow bound
    pew(:)    = exp( xalpw - xbetaw / zt(:) - xgamw * log( zt(:) ) )
    pew(:)    = zeps * pew(:) / ( ppres(:) - pew(:) )
    !
    plv(:)    = xlvtt + ( xcpv - xcl ) * ( zt(:) - xtt ) ! compute l_v
    pls(:)    = xlstt + ( xcpv - xci ) * ( zt(:) - xtt ) ! compute l_i
    !    
    pcph(:)   = xcpd + xcpv * pew(:)                     ! compute c_ph 
    !
  end subroutine convect_satmixratio


end module modd_convparext
