module flake
  use decomp,only:realkind
  implicit none
  private

  type,public::lakeType
     integer::lake_types=3
     integer::lake_no_prog=12
     integer::lake_no_tend=12
     integer::lake_no_diag=15
     real(kind=realkind),allocatable,dimension(:,:,:,:)::prog_lakes !(klon,klat,lake_types,lake_no_prog)
     real(kind=realkind),allocatable,dimension(:,:,:,:)::tend_lakes !(klon,klat,lake_types,lake_no_tend)
     real(kind=realkind),allocatable,dimension(:,:,:,:)::diag_lakes !(klon,klat,lake_types,lake_no_diag)
     real(kind=realkind),pointer::frac_lakes(:,:,:)
     real(kind=realkind),pointer::depth_lakes(:,:,:)
  end type lakeType

!  integer, parameter :: ireals    = selected_real_kind (12,200)
   integer, parameter :: ireals=realkind
  ! number of desired significant digits for
  ! real variables
  ! corresponds to 8 byte real variables

  integer, parameter ::iintegers = kind(1)
  ! kind-type parameter of the integer values
  ! corresponds to the default integers


  real(kind=realkind) :: z0ice = 8.e-04_realkind        !  roughness for ice
  real(kind=realkind) :: z0snowice = 0.005_realkind         !  roughness for snow on ice

  integer,public, parameter ::  nlaketype=3_iintegers ! number of lake types - may be declared elsewhere

  real (kind = ireals),public, parameter :: fraclakemin = 0.01_ireals ! minimum lake fraction to be taken into account
  !h_snow_min_flk = 1.0e-5_ireals        , & ! minimum snow thickness [m] - from flake
  ! for rca coupling we change the minimum snow thickness to 1.0e-3
  real (kind = ireals),public, parameter :: h_snow_min_flk = 1.0e-3_ireals         ! minimum snow thickness [m] - from flake
  !h_ice_min_flk  = 1.0e-9_ireals            ! minimum ice thickness [m] - from flake
  ! for rca coupling we change the minimum ice thickness to 1.0e-3
  real (kind = ireals), public,parameter :: h_ice_min_flk  = 1.0e-3_ireals            ! minimum ice thickness [m] - for rca


  !  albedo for water, ice and snow.
  real(kind=ireals),parameter::albedo_water_ref=0.07_ireals     !water
  real(kind=ireals),public,parameter::albedo_whiteice_ref=0.75_ireals   !white ice
  real(kind=ireals),public,parameter::albedo_blueice_ref=0.10_ireals    ! blue ice
  real(kind=ireals),parameter::albedo_drysnow_ref=0.75_ireals    ! dry snow 
  real(kind=ireals),parameter::albedo_meltingsnow_ref=0.10_ireals!melting snow 

  !  empirical parameters.
  real(kind = ireals),public,parameter ::c_albice_mr = 95.6_ireals
  ! constant in the interpolation formula for 
  ! the ice albedo (mironov and ritter 2004)


  real (kind = ireals), parameter ::         &
       c_cbl_1       = 0.17_ireals                   , & ! constant in the cbl entrainment equation
       c_cbl_2       = 1.0_ireals                     , & ! constant in the cbl entrainment equation
       c_sbl_zm_n    = 0.5_ireals                    , & ! constant in the zm1996 equation for the equilibrium sbl depth
       c_sbl_zm_s    = 10.0_ireals                    , & ! constant in the zm1996 equation for the equilibrium sbl depth
       c_sbl_zm_i    = 20.0_ireals                    , & ! constant in the zm1996 equation for the equilibrium sbl depth
                                !  c_relax_h     = 0.1_ireals                    , & ! constant in the relaxation equation for the sbl depth
       c_relax_h     = 0.01_ireals                    , & ! constant in the relaxation equation for the sbl depth
                                !  c_relax_c     = 0.01_ireals                       ! constant in the relaxation equation for the shape factor
       c_relax_c     = 0.003_ireals                       ! constant in the relaxation equation for the shape factor
  ! with respect to the temperature profile in the thermocline

  !  parameters of the shape functions 
  !  indices refer to t - thermocline, s - snow, i - ice,
  !  b1 - upper layer of the bottom sediments, b2 - lower layer of the bottom sediments.
  !  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
  !  at "zeta=0" ad "zeta=1", respectively.
  real (kind = ireals), parameter ::         &
       c_t_min       = 0.5_ireals                    , & ! minimum value of the shape factor c_t (thermocline)
       c_t_max       = 0.8_ireals                    , & ! maximum value of the shape factor c_t (thermocline)
       phi_t_pr0_1   = 40.0_ireals/3.0_ireals   , & ! constant in the expression for the t shape-function derivative 
       phi_t_pr0_2   = 20.0_ireals/3.0_ireals   , & ! constant in the expression for the t shape-function derivative 
       c_tt_1        = 11.0_ireals/18.0_ireals  , & ! constant in the expression for c_tt (thermocline)
       c_tt_2        = 7.0_ireals/45.0_ireals   , & ! constant in the expression for c_tt (thermocline)
       c_b1          = 2.0_ireals/3.0_ireals    , & ! shape factor (upper layer of bottom sediments)
       c_b2          = 3.0_ireals/5.0_ireals    , & ! shape factor (lower layer of bottom sediments)
       phi_b1_pr0    = 2.0_ireals                     , & ! b1 shape-function derivative 
       c_s_lin       = 0.5_ireals                    , & ! shape factor (linear temperature profile in the snow layer)
       phi_s_pr0_lin = 1.0_ireals                    , & ! s shape-function derivative (linear profile) 
       c_i_lin       = 0.5_ireals                    , & ! shape factor (linear temperature profile in the ice layer)
       phi_i_pr0_lin = 1.0_ireals                     , & ! i shape-function derivative (linear profile) 
       phi_i_pr1_lin = 1.0_ireals                     , & ! i shape-function derivative (linear profile) 
       phi_i_ast_mr  = 2.0_ireals                     , & ! constant in the mr2003 expression for i shape function
       c_i_mr        = 1.0_ireals/12.0_ireals   , & ! constant in the mr2003 expression for i shape factor
       h_ice_max     = 3.0_ireals                         ! maximum ice tickness in the mr2003 ice model [m] 

  !  thermodynamic parameters
  real (kind = ireals), parameter ::    &
       tpl_grav          = 9.81_ireals          , & ! acceleration due to gravity [m s^{-2}]
                                !  tpl_t_r           = 277.13_ireals        , & ! temperature of maximum density of fresh water [k  ]
       tpl_t_r           = 277.13_ireals        , & ! temperature of maximum density of fresh water [k], dbg, kk, sd, 21.04.05 
       tpl_t_f           = 273.15_ireals        , & ! fresh water freezing point [k]
       tpl_a_t           = 1.6509e-05_ireals    , & ! constant in the fresh-water equation of state [k^{-2}]
       tpl_rho_w_r       = 1.0e+03_ireals       , & ! maximum density of fresh water [kg m^{-3}]
       tpl_rho_i         = 9.1e+02_ireals       , & ! density of ice [kg m^{-3}]
       tpl_rho_s_min     = 1.0e+02_ireals       , & ! minimum snow density [kg m^{-3}]
       tpl_rho_s_max     = 4.0e+02_ireals       , & ! maximum snow density [kg m^{-3}]
       tpl_gamma_rho_s   = 2.0e+02_ireals       , & ! empirical parameter [kg m^{-4}]  
                                ! in the expression for the snow density 
       tpl_l_f           = 3.3e+05_ireals       , & ! latent heat of fusion [j kg^{-1}]
       tpl_c_w           = 4.2e+03_ireals       , & ! specific heat of water [j kg^{-1} k^{-1}]
       tpl_c_i           = 2.1e+03_ireals       , & ! specific heat of ice [j kg^{-1} k^{-1}]
       tpl_c_s           = 2.1e+03_ireals       , & ! specific heat of snow [j kg^{-1} k^{-1}]
       tpl_kappa_w       = 5.46e-01_ireals      , & ! molecular heat conductivity of water [j m^{-1} s^{-1} k^{-1}]
                                !cps090202  tpl_kappa_i       = 2.29_ireals          , & ! molecular heat conductivity of ice [j m^{-1} s^{-1} k^{-1}]
       tpl_kappa_i       = 1.5_ireals          , & ! molecular heat conductivity of ice [j m^{-1} s^{-1} k^{-1}]
       tpl_kappa_s_min   = 0.2_ireals           , & ! minimum molecular heat conductivity of snow [j m^{-1} s^{-1} k^{-1}]
       tpl_kappa_s_max   = 1.5_ireals           , & ! maximum molecular heat conductivity of snow [j m^{-1} s^{-1} k^{-1}]
       tpl_gamma_kappa_s = 1.3_ireals               ! empirical parameter [j m^{-2} s^{-1} k^{-1}] 
  ! in the expression for the snow heat conductivity 

  !  maximum value of the wave-length bands 
  !  in the exponential decay law for the radiation flux.
  !  a storage for a ten-band approximation is allocated,
  !  although a smaller number of bands can actually be used.
  integer (kind = iintegers), parameter ::   nband_optic_max = 10

  !  define type "opticpar_medium"
  type opticpar_medium
     integer (kind = iintegers)                        ::   & 
          nband_optic                                            ! number of wave-length bands
     real (kind = ireals), dimension (nband_optic_max) ::   & 
          frac_optic                                         , & ! fractions of total radiation flux 
          extincoef_optic                                        ! extinction coefficients                 
  end type opticpar_medium



  ! optical characteristics of water, ice and snow 

  integer (kind = iintegers), private :: i ! help variable(s) ! for pgi - compiler
  ! do loop index

  type (opticpar_medium), parameter :: & ! opticpar_water_ref, optical characteristics of water
       opticpar_water =  opticpar_medium(1_iintegers,                    &
       (/1.0_ireals, (0.0_ireals,i=2_iintegers,nband_optic_max)/),            &
       (/5.0_ireals, (1.e+10_ireals,i=2_iintegers,nband_optic_max)/))

  type (opticpar_medium), parameter ::opticpar_ice = opticpar_medium(1_iintegers,   & ! opticpar_ice_opaque, opaque ice
       (/1.0_ireals, (0.0_ireals,i=2_iintegers,nband_optic_max)/),            &
       (/1.0e+07_ireals, (1.e+10_ireals,i=2_iintegers,nband_optic_max)/)) 

  type (opticpar_medium), parameter ::opticpar_snow = opticpar_medium(1_iintegers, & ! opticpar_snow_opaque, opaque snow
       (/1.0_ireals, (0.0_ireals,i=2_iintegers,nband_optic_max)/),            &
       (/1.0e+07_ireals, (1.e+10_ireals,i=2_iintegers,nband_optic_max)/)) 


  real (kind = ireals), parameter :: depth_bs_lk = 5.0_ireals           ! thickness of the thermally active layer of the bottom sediments [m]
  real (kind = ireals), parameter :: t_bs_lk = 278.13_ireals           ! temperature at the outer edge of 
  ! the thermally active layer of bottom sediments [k],
  ! maximum on the bottom (now is dummy - reference value)
  ! tpl_t_r+1.0_ireals

  ! physical (not of the lake) parameters - better to use other units of rca, but not now :-)
  real (kind = ireals), parameter :: omega_earth     = 7.29e-05_ireals            ! the angular velocity of the earth's rotation [s^{-1}]


  logical, parameter :: lflk_botsed_use   = .true.       ! .true. indicates that the bottom-sediment scheme should be used
  ! to compute the depth penetrated by the thermal wave, 
  ! the temperature at this depth and the bottom heat flux.
  ! otherwise, the heat flux at the water-bottom sediment interface
  ! is set to zero, the depth penetrated by the thermal wave 
  ! is set to a reference value defined below,
  ! and the temperature at this depth is set to 
  ! the temperature of maximum density of the fresh water.

  real (kind = ireals), parameter :: rflk_depth_bs_ref = 5.0_ireals
  ! reference value of the depth of the thermally active
  ! reference value of the depth of the thermally active
  ! layer of bottom sediments [m].
  ! this value is used to (formally) define
  ! the depth penetrated by the thermal wave
  ! in case the bottom-sediment scheme is not used.



  !  "flake" (fresh-water lake) is a lake model capable of predicting the surface temperature 
  !  in lakes of various depth on the time scales from a few hours to a year.
  !  the model is based on a two-layer parametric representation of
  !  the temperature profile, where the structure of the stratified layer between the
  !  upper mixed layer and the basin bottom, the lake thermocline,
  !  is described using the concept of self-similarity of the evolving temperature profile.
  !  the concept was put forward by kitaigorodskii and miropolsky (1970) 
  !  to describe the vertical temperature structure of the oceanic seasonal thermocline.
  !  it has since been successfully used in geophysical applications.
  !  the concept of self-similarity of the evolving temperature profile
  !  is also used to describe the vertical structure of the thermally active upper layer 
  !  of bottom sediments and of the ice and snow cover.
  !
  !  the lake model incorporates the heat budget equations
  !  for the four layers in question, viz., snow, ice, water and bottom sediments,
  !  developed with due regard for the vertically distributed character
  !  of the short-wave radiation heating.
  !  the entrainment equation that incorporates the zilitinkevich (1975) spin-up term
  !  is used to compute the depth of a convectively-mixed layer. 
  !  a relaxation-type equation is used
  !  to compute the wind-mixed layer depth in stable and neutral stratification,
  !  where a multi-limit formulation for the equilibrium mixed-layer depth
  !  proposed by zilitinkevich and mironov (1996)
  !  accounts for the effects of the earth's rotation, of the surface buoyancy flux
  !  and of the static stability in the thermocline.
  !  the equations for the mixed-layer depth are developed with due regard for  
  !  the volumetric character of the radiation heating.
  !  simple thermodynamic arguments are invoked to develop
  !  the evolution equations for the ice thickness and for the snow thickness.
  !  the heat flux through the water-bottom sediment interface is computed,
  !  using a parameterization proposed by golosov et al. (1998).
  !  the heat flux trough the air-water interface 
  !  (or through the air-ice or air-snow interface)
  !  is provided by the driving atmospheric model.
  !
  !  empirical constants and parameters of the lake model
  !  are estimated, using independent empirical and numerical data.
  !  they should not be re-evaluated when the model is applied to a particular lake.
  !  the only lake-specific parameters are the optical characteristics of lake water,
  !  the temperature at the bottom of the thermally active layer
  !  of bottom sediments and the depth of this layer.
  !
  !  a detailed description of the lake model is given in
  !  mironov, d. v., et al., 2003:
  !  parameterization of lakes in numerical weather prediction.
  !  part 1: description of a lake model.
  !  details of the ice-snow model are given in
  !  mironov, d., and b. ritter, 2003:
  !  a thermodynamic model of ice and snow over water bodies.
  !  manuscripts are available from the authors.
  !  contact dmitrii mironov 
  !  german weather service, referat fe14,
  !  frankfurter str. 135, d-63067 offenbach am main, germany.
  !  e-mail: dmitrii.mironov@dwd.de 
  !
  !  apart from serving as a lake module (parameterisation scheme) 
  !  within a three-dimensional nwp or climate model,
  !  "flake" can serve as an independent single-column lake model.
  !  the code is easily re-configured by setting a number of logical switches accordingly.
  !  this is explained in comment lines where appropriate.
  !
  !
  ! current code owner: dwd, dmitrii mironov
  !  phone:  +49-69-8062 2705
  !  fax:    +49-69-8062 3721
  !  e-mail: dmitrii.mironov@dwd.de
  !
  ! history:
  ! version    date       name
  ! ---------- ---------- ----
  ! n.nn       yyyy/mm/dd dmitrii mironov
  !  initial release
  ! !version 1!  !2004.05.11!  katherina kourzeneva
  ! move   h_snow_min_flk and  h_ice_min_flk to 

  ! code description:
  ! language: fortran 90.
  ! software standards: "european standards for writing and
  ! documenting exchangeable fortran 90 code".

  !
  ! declarations:
  !
  ! modules used:



  !  the variables declared below
  !  are accessible to all program units of the module "flake"
  !  and to the driving routines that use "flake".
  !  these are basically the quantities computed by flake.
  !  apart from these quantities, there are a few local scalars 
  !  used by flake routines mainly for security reasons.
  !  these are declared as private.
  !  all variables declared below have a suffix "flk".

  !  flake variables of type real

  !  temperatures at the previous time step ("p") and the updated temperatures ("n") 
  real (kind = ireals) ::           &
       t_mnw_p_flk, t_mnw_n_flk      , & ! mean temperature of the water column [k] 
       t_snow_p_flk, t_snow_n_flk    , & ! temperature at the air-snow interface [k] 
       t_ice_p_flk, t_ice_n_flk      , & ! temperature at the snow-ice or air-ice interface [k] 
       t_wml_p_flk, t_wml_n_flk      , & ! mixed-layer temperature [k] 
       t_bot_p_flk, t_bot_n_flk      , & ! temperature at the water-bottom sediment interface [k] 
       t_b1_p_flk, t_b1_n_flk            ! temperature at the bottom of the upper layer of the sediments [k] 

  !  thickness of various layers at the previous time step ("p") and the updated values ("n") 
  real (kind = ireals) ::           &
       h_snow_p_flk, h_snow_n_flk    , & ! snow thickness [m]
       h_ice_p_flk, h_ice_n_flk      , & ! ice thickness [m]
       h_ml_p_flk, h_ml_n_flk        , & ! thickness of the mixed-layer [m] 
       h_b1_p_flk, h_b1_n_flk            ! thickness of the upper layer of bottom sediments [m] 

  !  the shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
  real (kind = ireals) ::           &
       c_t_p_flk, c_t_n_flk          , & ! shape factor (thermocline)
       c_tt_flk                      , & ! dimensionless parameter (thermocline)
       c_q_flk                       , & ! shape factor with respect to the heat flux (thermocline)
       c_i_flk                       , & ! shape factor (ice)
       c_s_flk                           ! shape factor (snow)

  !  derivatives of the shape functions
  real (kind = ireals) ::           &
       phi_t_pr0_flk                 , & ! d\phi_t(0)/d\zeta   (thermocline)
       phi_i_pr0_flk                 , & ! d\phi_i(0)/d\zeta_i (ice)
       phi_i_pr1_flk                 , & ! d\phi_i(1)/d\zeta_i (ice)
       phi_s_pr0_flk                     ! d\phi_s(0)/d\zeta_s (snow)

  !  heat and radiation fluxes
  real (kind = ireals) ::           &
       q_snow_flk                    , & ! heat flux through the air-snow interface [w m^{-2}]
       q_ice_flk                     , & ! heat flux through the snow-ice or air-ice interface [w m^{-2}]
       q_w_flk                       , & ! heat flux through the ice-water or air-water interface [w m^{-2}]
       q_bot_flk                     , & ! heat flux through the water-bottom sediment interface [w m^{-2}]
       i_atm_flk                     , & ! radiation flux at the lower boundary of the atmosphere [w m^{-2}],
                                ! i.e. the incident radiation flux with no regard for the surface albedo.
       i_snow_flk                    , & ! radiation flux through the air-snow interface [w m^{-2}]
       i_ice_flk                     , & ! radiation flux through the snow-ice or air-ice interface [w m^{-2}]
       i_w_flk                       , & ! radiation flux through the ice-water or air-water interface [w m^{-2}]
       i_h_flk                       , & ! radiation flux through the mixed-layer-thermocline interface [w m^{-2}]
       i_bot_flk                     , & ! radiation flux through the water-bottom sediment interface [w m^{-2}]
       i_intm_0_h_flk                , & ! integral-mean radiation flux over the mixed layer [w m^{-1}]
       i_intm_h_d_flk                , & ! integral-mean radiation flux over the thermocline [w m^{-1}]
       q_star_flk                        ! a generalized heat flux scale [w m^{-2}]

  !  velocity scales
  real (kind = ireals) ::           &
       u_star_w_flk                  , & ! friction velocity in the surface layer of lake water [m s^{-1}]
       w_star_sfc_flk                    ! convective velocity scale, 
  ! using a generalized heat flux scale [m s^{-1}]

  !  the rate of snow accumulation
  real(kind = ireals)::dmsnowdt_flk !the rate of snow accumulation [kg m^{-2} s^{-1}]

  real (kind = ireals), parameter :: &
       h_ml_min_flk   = 1.0e-2_ireals        , & ! minimum mixed-layer depth [m] 
       h_ml_max_flk   = 1.0e+3_ireals        , & ! maximum mixed-layer depth [m] 
       h_b1_min_flk   = 1.0e-3_ireals        , & ! minimum thickness of the upper layer of bottom sediments [m] 
       u_star_min_flk = 1.0e-6_ireals            ! minimum value of the surface friction velocity [m s^{-1}]

  real (kind = ireals), parameter:: &
       c_small_flk    = 1.0e-10_ireals, &          ! a small number
       c_small_4_flk  = 1.0e-4_ireals              ! kk, 070205, for rca:  rca has much smaller accuracy!  


  public lakemasks,slfluxo_surf_lake_ice,lakeinit
contains
  subroutine flake_driver ( depth_w, depth_bs, t_bs, par_coriolis,       &
       extincoef_water_typ,                         &
       del_time, &
       t_sfc_n,    checkout)

    ! Description:
    !
    !  The main driving routine of the FLake model,
    !  where computations are performed.
    !  Advances the surface temperature
    !  and other FLake variables one time step.
    !  At the moment, the Euler explicit scheme is used.
    !
    !  Lines embraced with "!_tmp" contain temporary parts of the code.
    !  These should be removed prior to using FLake in applications.
    !  Lines embraced/marked with "!_dev" may be replaced
    !  as improved parameterizations are developed and tested.
    !  Lines embraced/marked with "!_dm" are DM's comments
    !  that may be helpful to a user.
    !  Lines embraced/marked with "!_dbg" are used 
    !  for debugging purposes only.



    implicit none

    ! Declarations

    !  Input (procedure arguments)
    real (kind = ireals), intent(in) ::   &
         depth_w                           , & ! the lake depth [m]
         depth_bs                          , & ! depth of the thermally active layer of the bottom sediments [m]
         t_bs                              , & ! temperature at the outer edge of 
                                ! the thermally active layer of the bottom sediments [k]
         par_coriolis                      , & ! the coriolis parameter [s^{-1}]
         extincoef_water_typ               , & ! "typical" extinction coefficient of the lake water [m^{-1}],
                                ! used to compute the equilibrium cbl depth
         del_time                           ! the model time step [s]

    ! (equal to either t_ice, t_snow or to t_wml)

    !  output (procedure arguments)

    real(kind = ireals), intent(out) ::t_sfc_n! updated surface temperature [k] 
                                ! (equal to the updated value of either t_ice, t_snow or t_wml)
    real(kind = ireals), intent(out) ::checkout!checking variable - to check something (optional)



    !  local variables of type logical
    logical :: l_ice_create    , & ! switch, true = ice does not exist but should be created
         l_snow_exists   , & ! switch, true = there is snow above the ice
         l_ice_meltabove     ! switch, true = snow/ice melting from above takes place


    !  local variables of type real
    real (kind = ireals) ::    &
         d_t_mnw_dt             , & ! time derivative of t_mnw [k s^{-1}] 
         d_t_ice_dt             , & ! time derivative of t_ice [k s^{-1}] 
         d_t_bot_dt             , & ! time derivative of t_bot [k s^{-1}] 
         d_t_b1_dt              , & ! time derivative of t_b1 [k s^{-1}] 
         d_h_snow_dt            , & ! time derivative of h_snow [m s^{-1}]
         d_h_ice_dt             , & ! time derivative of h_ice [m s^{-1}]
         d_h_ml_dt              , & ! time derivative of h_ml [m s^{-1}]
         d_h_b1_dt              , & ! time derivative of h_b1 [m s^{-1}]
         d_c_t_dt                   ! time derivative of c_t [s^{-1}]

    !  local variables of type real
    real (kind = ireals) ::    &
         n_t_mean               , & ! the mean buoyancy frequency in the thermocline [s^{-1}] 
         zm_h_scale             , & ! the zm96 equilibrium sbl depth scale [m] 
         conv_equil_h_scale         ! the equilibrium cbl depth scale [m]

    !  local variables of type real
    real (kind = ireals) :: &
         h_ice_threshold     , & ! if h_ice<h_ice_threshold, use quasi-equilibrium ice model 
         flk_str_1           , & ! help storage variable
         flk_str_2           , & ! help storage variable
         r_h_icesnow         , & ! dimensionless ratio, used to store intermediate results
         r_rho_c_icesnow     , & ! dimensionless ratio, used to store intermediate results
         r_ti_icesnow        , & ! dimensionless ratio, used to store intermediate results
         r_tstar_icesnow         ! dimensionless ratio, used to store intermediate results


    !------------------------------------------------------------------------------
    !  compute fluxes, using variables from the previous time step.
    !------------------------------------------------------------------------------

    !_dm
    ! at this point, the heat and radiation fluxes, namely,
    ! q_snow_flk, q_ice_flk, q_w_flk, 
    ! i_atm_flk, i_snow_flk, i_ice_flk, i_w_flk, i_h_flk, i_bot_flk,     
    ! the integral-mean radiation flux over the mixed layer, i_intm_0_h_flk, 
    ! and the integral-mean radiation flux over the thermocline, i_intm_h_d_flk, 
    ! should be known.
    ! they are computed within "flake_interface" (or within the driving model)
    ! and are available to "flake_driver"
    ! through the above variables declared in the module "flake".
    ! in case a lake is ice-covered, q_w_flk is re-computed below.
    !_dm

    ! heat flux through the ice-water interface

    if(h_ice_p_flk>=h_ice_min_flk) then    ! ice exists 
       if(h_ml_p_flk<=h_ml_min_flk) then    ! mixed-layer depth is zero, compute flux 
          q_w_flk = -tpl_kappa_w*(t_bot_p_flk-t_wml_p_flk)/depth_w  ! flux with linear t(z) 
          phi_t_pr0_flk = phi_t_pr0_1*c_t_p_flk-phi_t_pr0_2         ! d\phi(0)/d\zeta (thermocline)
          q_w_flk = q_w_flk*max(phi_t_pr0_flk, 1.0_ireals)           ! account for an increased d\phi(0)/d\zeta 
       else                    
          q_w_flk = 0.0_ireals                  ! mixed-layer depth is greater than zero, set flux to zero
       end if
    end if

    ! a generalized heat flux scale 
    q_star_flk = q_w_flk + i_w_flk + i_h_flk - 2.0_ireals*i_intm_0_h_flk

    ! heat flux through the water-bottom sediment interface
    if(lflk_botsed_use) then
       q_bot_flk = -tpl_kappa_w*(t_b1_p_flk-t_bot_p_flk)/max(h_b1_p_flk, h_b1_min_flk)*phi_b1_pr0
    else  
       q_bot_flk = 0.0_ireals   ! the bottom-sediment scheme is not used
    end if


    !------------------------------------------------------------------------------
    !  check if ice exists or should be created.
    !  if so, compute the thickness and the temperature of ice and snow.
    !------------------------------------------------------------------------------

    !_dm
    ! notice that a quasi-equilibrium ice-snow model is used 
    ! to avoid numerical instability when the ice is thin.
    ! this is always the case when new ice is created.
    !_dm

    !_dev
    ! the dependence of snow density and of snow heat conductivity 
    ! on the snow thickness is accounted for parametrically.
    ! that is, the time derivatives of \rho_s and \kappa_s are neglected.
    ! the exception is the equation for the snow thickness 
    ! in case of snow accumulation and no melting, 
    ! where d\rho_s/dt is incorporated.
    ! furthermore, some (presumably small) correction terms incorporating 
    ! the snow density and the snow heat conductivity are dropped out.
    ! those terms may be included as better formulations 
    ! for \rho_s and \kappa_s are available.
    !_dev

    ! default values
    l_ice_create    = .false.  
    l_ice_meltabove = .false.  


    ice_exist: if(h_ice_p_flk<h_ice_min_flk) then   ! ice does not exist 

       l_ice_create = t_wml_p_flk<=(tpl_t_f+c_small_4_flk).and.q_w_flk<0.0_ireals
       if(l_ice_create) then                            ! ice does not exist but should be created
          d_h_ice_dt = -q_w_flk/tpl_rho_i/tpl_l_f                                  
          h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          ! advance h_ice 
          t_ice_n_flk = tpl_t_f + h_ice_n_flk*q_w_flk/tpl_kappa_i/phi_i_pr0_lin    ! ice temperature
          d_h_snow_dt = dmsnowdt_flk/tpl_rho_s_min 
          h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       ! advance h_snow
          phi_i_pr1_flk = phi_i_pr1_lin                                    & 
               + phi_i_ast_mr*min(1.0_ireals, h_ice_n_flk/h_ice_max)       ! d\phi_i(1)/d\zeta_i (ice)
          r_h_icesnow = phi_i_pr1_flk/phi_s_pr0_lin*tpl_kappa_i/flake_snowheatconduct(h_snow_n_flk) &
               * h_snow_n_flk/max(h_ice_n_flk, h_ice_min_flk)
          t_snow_n_flk = t_ice_n_flk + r_h_icesnow*(t_ice_n_flk-tpl_t_f)           ! snow temperature
       else !kk, 30092004 - to avoid unsertainty here
          h_ice_n_flk = h_ice_p_flk
          h_snow_n_flk = h_snow_p_flk
          t_ice_n_flk = t_ice_p_flk
          t_snow_n_flk = t_snow_p_flk
       end if

    else ice_exist                                     ! ice exists

       l_snow_exists = h_snow_p_flk>=h_snow_min_flk   ! check if there is snow above the ice


       melting: if(t_snow_p_flk>=(tpl_t_f-c_small_4_flk)) then  ! t_sfc = t_f, check for melting from above
          ! t_snow = t_ice if snow is absent 
          if(l_snow_exists) then   ! there is snow above the ice
             flk_str_1 = q_snow_flk + i_snow_flk - i_ice_flk        ! atmospheric forcing
             if(flk_str_1>=0.0_ireals) then  ! melting of snow and ice from above
                l_ice_meltabove = .true.
                d_h_snow_dt = (-flk_str_1/tpl_l_f+dmsnowdt_flk)/flake_snowdensity(h_snow_p_flk)
                d_h_ice_dt  = -(i_ice_flk - i_w_flk - q_w_flk)/tpl_l_f/tpl_rho_i 
             end if
          else                     ! no snow above the ice
             flk_str_1 = q_ice_flk + i_ice_flk - i_w_flk - q_w_flk  ! atmospheric forcing + heating from the water
             if(flk_str_1>=0.0_ireals) then  ! melting of ice from above, snow accumulation may occur
                l_ice_meltabove = .true.
                d_h_ice_dt  = -flk_str_1/tpl_l_f/tpl_rho_i 
                d_h_snow_dt = dmsnowdt_flk/tpl_rho_s_min
             end if
          end if
          if(l_ice_meltabove) then  ! melting from above takes place
             h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! advance h_ice
             h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! advance h_snow
             t_ice_n_flk  = tpl_t_f                              ! set t_ice to the freezing point
             t_snow_n_flk = tpl_t_f                              ! set t_snow to the freezing point
          end if

       end if melting

       no_melting: if(.not.l_ice_meltabove) then                 ! no melting from above

          d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
          if(d_h_snow_dt<tpl_rho_s_max) then    ! account for d\rho_s/dt
             flk_str_1 = h_snow_p_flk*tpl_gamma_rho_s/tpl_rho_w_r
             flk_str_1 = flk_str_1/(1.0_ireals-flk_str_1)
          else                                     ! snow density is equal to its maximum value, d\rho_s/dt=0
             flk_str_1 = 0.0_ireals
          end if
          d_h_snow_dt = dmsnowdt_flk/d_h_snow_dt/(1.0_ireals+flk_str_1)       ! snow accumulation
          h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         ! advance h_snow

          phi_i_pr0_flk = h_ice_p_flk/h_ice_max                              ! h_ice relative to its maximum value
          c_i_flk = c_i_lin - c_i_mr*(1.0_ireals+phi_i_ast_mr)*phi_i_pr0_flk  ! shape factor (ice)
          phi_i_pr1_flk = phi_i_pr1_lin + phi_i_ast_mr*phi_i_pr0_flk         ! d\phi_i(1)/d\zeta_i (ice)
          phi_i_pr0_flk = phi_i_pr0_lin - phi_i_pr0_flk                      ! d\phi_i(0)/d\zeta_i (ice)

          h_ice_threshold = max(1.0_ireals, 2.0_ireals*c_i_flk*tpl_c_i*(tpl_t_f-t_ice_p_flk)/tpl_l_f)
          h_ice_threshold = phi_i_pr0_flk/c_i_flk*tpl_kappa_i/tpl_rho_i/tpl_c_i*h_ice_threshold
          h_ice_threshold = sqrt(h_ice_threshold*del_time)                   ! threshold value of h_ice
          h_ice_threshold = min(0.9_ireals*h_ice_max, max(h_ice_threshold, h_ice_min_flk))
          ! h_ice(threshold) < 0.9*h_ice_max

          if(h_ice_p_flk<h_ice_threshold) then  ! use a quasi-equilibrium ice model

             if(l_snow_exists) then   ! use fluxes at the air-snow interface
                flk_str_1 = q_snow_flk + i_snow_flk - i_w_flk
             else                     ! use fluxes at the air-ice interface
                flk_str_1 = q_ice_flk + i_ice_flk - i_w_flk
             end if
             d_h_ice_dt = -(flk_str_1-q_w_flk)/tpl_l_f/tpl_rho_i
             h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         ! advance h_ice
             t_ice_n_flk = tpl_t_f + h_ice_n_flk*flk_str_1/tpl_kappa_i/phi_i_pr0_flk  ! ice temperature

          else                                     ! use a complete ice model

             d_h_ice_dt = tpl_kappa_i*(tpl_t_f-t_ice_p_flk)/h_ice_p_flk*phi_i_pr0_flk
             d_h_ice_dt = (q_w_flk+d_h_ice_dt)/tpl_l_f/tpl_rho_i
             h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time                         ! advance h_ice

             r_ti_icesnow = tpl_c_i*(tpl_t_f-t_ice_p_flk)/tpl_l_f         ! dimensionless parameter
             r_tstar_icesnow = 1.0_ireals - c_i_flk                        ! dimensionless parameter
             if(l_snow_exists) then  ! there is snow above the ice
                r_h_icesnow = phi_i_pr1_flk/phi_s_pr0_lin*tpl_kappa_i/flake_snowheatconduct(h_snow_p_flk) &
                     * h_snow_p_flk/h_ice_p_flk
                r_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*tpl_c_s/tpl_rho_i/tpl_c_i 
                r_tstar_icesnow = r_tstar_icesnow*r_ti_icesnow             ! dimensionless parameter

                flk_str_2 = q_snow_flk+i_snow_flk-i_w_flk                  ! atmospheric fluxes
                flk_str_1  = c_i_flk*h_ice_p_flk + (1.0_ireals+c_s_lin*r_h_icesnow)*r_rho_c_icesnow*h_snow_p_flk
                d_t_ice_dt = -(1.0_ireals-2.0_ireals*c_s_lin)*r_h_icesnow*(tpl_t_f-t_ice_p_flk)             & 
                     * tpl_c_s*dmsnowdt_flk                          ! effect of snow accumulation
             else                    ! no snow above the ice
                r_tstar_icesnow = r_tstar_icesnow*r_ti_icesnow             ! dimensionless parameter
                flk_str_2 = q_ice_flk+i_ice_flk-i_w_flk                    ! atmospheric fluxes
                flk_str_1  = c_i_flk*h_ice_p_flk
                d_t_ice_dt = 0.0_ireals
             end if
             d_t_ice_dt = d_t_ice_dt + tpl_kappa_i*(tpl_t_f-t_ice_p_flk)/h_ice_p_flk*phi_i_pr0_flk       &
                  * (1.0_ireals-r_tstar_icesnow)                     ! add flux due to heat conduction
             d_t_ice_dt = d_t_ice_dt - r_tstar_icesnow*q_w_flk            ! add flux from water to ice
             d_t_ice_dt = d_t_ice_dt + flk_str_2                          ! add atmospheric fluxes
             d_t_ice_dt = d_t_ice_dt/tpl_rho_i/tpl_c_i                    ! total forcing
             d_t_ice_dt = d_t_ice_dt/flk_str_1                            ! dt_ice/dt 
             t_ice_n_flk = t_ice_p_flk + d_t_ice_dt*del_time                          ! advance t_ice
          end if

          phi_i_pr1_flk = min(1.0_ireals, h_ice_n_flk/h_ice_max)          ! h_ice relative to its maximum value
          phi_i_pr1_flk = phi_i_pr1_lin + phi_i_ast_mr*phi_i_pr1_flk     ! d\phi_i(1)/d\zeta_i (ice)
          r_h_icesnow = phi_i_pr1_flk/phi_s_pr0_lin*tpl_kappa_i/flake_snowheatconduct(h_snow_n_flk) &
               *h_snow_n_flk/max(h_ice_n_flk, h_ice_min_flk)
          t_snow_n_flk = t_ice_n_flk + r_h_icesnow*(t_ice_n_flk-tpl_t_f)             ! snow temperature

       end if no_melting


    end if ice_exist

    ! security, limit h_ice by its maximum value
    h_ice_n_flk = min(h_ice_n_flk, h_ice_max) 

    ! security, limit the ice and snow temperatures by the freezing point 
    t_snow_n_flk = min(t_snow_n_flk, tpl_t_f)  
    t_ice_n_flk =  min(t_ice_n_flk,  tpl_t_f) 

    ! security, avoid too low values (these constraints are used for debugging purposes)
    t_snow_n_flk = max(t_snow_n_flk, 73.15_ireals)  
    t_ice_n_flk =  max(t_ice_n_flk,  73.15_ireals)    
    !_tmp

    ! remove too thin ice and/or snow
    if(h_ice_n_flk<h_ice_min_flk)  then        ! check ice
       h_ice_n_flk = 0.0_ireals       ! ice is too thin, remove it, and
       t_ice_n_flk = tpl_t_f         ! set t_ice to the freezing point.
       h_snow_n_flk = 0.0_ireals      ! remove snow when there is no ice, and
       t_snow_n_flk = tpl_t_f        ! set t_snow to the freezing point.
       l_ice_create = .false.        ! "exotic" case, ice has been created but proved to be too thin
    else if(h_snow_n_flk<h_snow_min_flk) then  ! ice exists, check snow
       h_snow_n_flk = 0.0_ireals      ! snow is too thin, remove it, 
       t_snow_n_flk = t_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
    end if


    !------------------------------------------------------------------------------
    !  compute the mean temperature of the water column.
    !------------------------------------------------------------------------------

    if(l_ice_create) q_w_flk = 0.0_ireals     ! ice has just been created, set q_w to zero
    d_t_mnw_dt = (q_w_flk - q_bot_flk + i_w_flk - i_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w
    t_mnw_n_flk = t_mnw_p_flk + d_t_mnw_dt*del_time   ! advance t_mnw
    t_mnw_n_flk = max(t_mnw_n_flk, tpl_t_f)           ! limit t_mnw by the freezing point 


    !------------------------------------------------------------------------------
    !  compute the mixed-layer depth, the mixed-layer temperature, 
    !  the bottom temperature and the shape factor
    !  with respect to the temperature profile in the thermocline. 
    !  different formulations are used, depending on the regime of mixing. 
    !------------------------------------------------------------------------------


    htc_water: if(h_ice_n_flk>=h_ice_min_flk) then    ! ice exists

       t_mnw_n_flk = min(t_mnw_n_flk, tpl_t_r) ! limit the mean temperature under the ice by t_r 
       t_wml_n_flk = tpl_t_f                   ! the mixed-layer temperature is equal to the freezing point 

       if(l_ice_create) then                  ! ice has just been created 
          if(h_ml_p_flk>=depth_w-h_ml_min_flk) then    ! h_ml=d when ice is created 
             h_ml_n_flk = 0.0_ireals                 ! set h_ml to zero 
             c_t_n_flk = c_t_min                    ! set c_t to its minimum value 
          else                                          ! h_ml<d when ice is created 
             h_ml_n_flk = h_ml_p_flk                ! h_ml remains unchanged 
             c_t_n_flk = c_t_p_flk                  ! c_t (thermocline) remains unchanged 
          end if
          t_bot_n_flk = t_wml_n_flk - (t_wml_n_flk-t_mnw_n_flk)/c_t_n_flk/(1.0_ireals-h_ml_n_flk/depth_w)
          ! update the bottom temperature 

       else if(t_bot_p_flk<tpl_t_r) then   ! ice exists and t_bot < t_r, molecular heat transfer 
          h_ml_n_flk = h_ml_p_flk                  ! h_ml remains unchanged 
          c_t_n_flk = c_t_p_flk                    ! c_t (thermocline) remains unchanged 
          t_bot_n_flk = t_wml_n_flk - (t_wml_n_flk-t_mnw_n_flk)/c_t_n_flk/(1.0_ireals-h_ml_n_flk/depth_w)
          ! update the bottom temperature 

       else                                   ! ice exists and t_bot = t_r, convection due to bottom heating 
          t_bot_n_flk = tpl_t_r                      ! t_bot is equal to the temperature of maximum density 
          if(h_ml_p_flk>=c_small_flk) then   ! h_ml > 0 
             c_t_n_flk = c_t_p_flk                     ! c_t (thermocline) remains unchanged 
             h_ml_n_flk = depth_w*(1.0_ireals-(t_wml_n_flk-t_mnw_n_flk)/(t_wml_n_flk-t_bot_n_flk)/c_t_n_flk)
             h_ml_n_flk = max(h_ml_n_flk, 0.0_ireals)   ! update the mixed-layer depth  
          else                                 ! h_ml = 0 
             h_ml_n_flk = h_ml_p_flk                   ! h_ml remains unchanged 
             c_t_n_flk = (t_wml_n_flk-t_mnw_n_flk)/(t_wml_n_flk-t_bot_n_flk) 
             c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min)) ! update the shape factor (thermocline)  
          end if
       end if

       t_bot_n_flk = min(t_bot_n_flk, tpl_t_r)    ! security, limit the bottom temperature by t_r 

    else htc_water                                      ! open water

       ! generalised buoyancy flux scale and convective velocity scale
       flk_str_1 = flake_buoypar(t_wml_p_flk)*q_star_flk/tpl_rho_w_r/tpl_c_w                    
       if(flk_str_1<0.0_ireals) then       
          w_star_sfc_flk = (-flk_str_1*h_ml_p_flk)**(1.0_ireals/3.0_ireals)  ! convection     
       else 
          w_star_sfc_flk = 0.0_ireals                                       ! neutral or stable stratification
       end if

       !_dm
       ! the equilibrium depth of the cbl due to surface cooling with the volumetric heating
       ! is not computed as a solution to the transcendental equation.
       ! instead, an algebraic formula is used
       ! that interpolates between the two asymptotic limits.
       !_dm
       conv_equil_h_scale = -q_w_flk/max(i_w_flk, c_small_flk)
       if(conv_equil_h_scale>0.0_ireals.and.conv_equil_h_scale<1.0_ireals  &
            .and.t_wml_p_flk>tpl_t_r) then    ! the equilibrium cbl depth scale is only used above t_r
          conv_equil_h_scale = sqrt(6.0_ireals*conv_equil_h_scale)               &
               + 2.0_ireals*conv_equil_h_scale/(1.0_ireals-conv_equil_h_scale)
          conv_equil_h_scale = min(depth_w, conv_equil_h_scale/extincoef_water_typ)
       else
          conv_equil_h_scale = 0.0_ireals       ! set the equilibrium cbl depth to zero
       end if

       ! mean buoyancy frequency in the thermocline
       n_t_mean = flake_buoypar(0.5_ireals*(t_wml_p_flk+t_bot_p_flk))*(t_wml_p_flk-t_bot_p_flk)
       if(h_ml_p_flk<=depth_w-h_ml_min_flk) then
          n_t_mean = sqrt(n_t_mean/(depth_w-h_ml_p_flk))  ! compute n                   
       else 
          n_t_mean = 0.0_ireals                            ! h_ml=d, set n to zero
       end if

       ! the rate of change of c_t
       d_c_t_dt = max(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2.0_ireals
       d_c_t_dt = n_t_mean*(depth_w-h_ml_p_flk)**2.0_ireals       &
            / c_relax_c/d_c_t_dt                               ! relaxation time scale for c_t
       d_c_t_dt = (c_t_max-c_t_min)/max(d_c_t_dt, c_small_flk)     ! rate-of-change of c_t 

       ! compute the shape factor and the mixed-layer depth, 
       ! using different formulations for convection and wind mixing

       c_tt_flk = c_tt_1*c_t_p_flk-c_tt_2         ! c_tt, using c_t at the previous time step
       c_q_flk = 2.0_ireals*c_tt_flk/c_t_p_flk     ! c_q using c_t at the previous time step

       mixing_regime: if(flk_str_1<0.0_ireals) then  ! convective mixing 

          checkout=(i_w_flk + i_h_flk - 2.0_ireals*i_intm_0_h_flk)
          c_t_n_flk = c_t_p_flk + d_c_t_dt*del_time                        ! update c_t, assuming dh_ml/dt>0
          c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min))                ! limit c_t 
          d_c_t_dt = (c_t_n_flk-c_t_p_flk)/del_time                        ! re-compute dc_t/dt

          if(h_ml_p_flk<=depth_w-h_ml_min_flk) then       ! compute dh_ml/dt
             if(h_ml_p_flk<=h_ml_min_flk) then    ! use a reduced entrainment equation (spin-up)
                d_h_ml_dt = c_cbl_1/c_cbl_2*max(w_star_sfc_flk, c_small_flk)

                !_dbg
                ! write(*,*) ' flake: reduced entrainment eq. d_time*d_h_ml_dt  = ', d_h_ml_dt*del_time
                ! write(*,*) '         w_*       = ', w_star_sfc_flk
                ! write(*,*) '         \beta*q_* = ', flk_str_1
                !_dbg

             else                                   ! use a complete entrainment equation 
                r_h_icesnow     = depth_w/h_ml_p_flk
                r_rho_c_icesnow = r_h_icesnow-1.0_ireals
                r_ti_icesnow    = c_t_p_flk/c_tt_flk
                r_tstar_icesnow = (r_ti_icesnow/2.0_ireals-1.0_ireals)*r_rho_c_icesnow + 1.0_ireals
                d_h_ml_dt = -q_star_flk*(r_tstar_icesnow*(1.0_ireals+c_cbl_1)-1.0_ireals) - q_bot_flk
                d_h_ml_dt = d_h_ml_dt/tpl_rho_w_r/tpl_c_w                        ! q_* and q_b flux terms
                flk_str_2 = (depth_w-h_ml_p_flk)*(t_wml_p_flk-t_bot_p_flk)*c_tt_2/c_tt_flk*d_c_t_dt 
                d_h_ml_dt = d_h_ml_dt + flk_str_2                                 ! add dc_t/dt term
                flk_str_2 = i_bot_flk + (r_ti_icesnow-1.0_ireals)*i_h_flk - r_ti_icesnow*i_intm_h_d_flk
                flk_str_2 = flk_str_2 + (r_ti_icesnow-2.0_ireals)*r_rho_c_icesnow*(i_h_flk-i_intm_0_h_flk)
                flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
                d_h_ml_dt = d_h_ml_dt + flk_str_2                                 ! add radiation terms
                flk_str_2 = -c_cbl_2*r_tstar_icesnow*q_star_flk/tpl_rho_w_r/tpl_c_w/max(w_star_sfc_flk, c_small_flk)
                flk_str_2 = flk_str_2 + c_t_p_flk*(t_wml_p_flk-t_bot_p_flk)
                d_h_ml_dt = d_h_ml_dt/flk_str_2                                   ! dh_ml/dt = r.h.s.
             end if
             !_dm
             ! notice that dh_ml/dt may appear to be negative  
             ! (e.g. due to buoyancy loss to bottom sediments and/or
             ! the effect of volumetric radiation heating),
             ! although a negative generalized buoyancy flux scale indicates 
             ! that the equilibrium cbl depth has not yet been reached
             ! and convective deepening of the mixed layer should take place.
             ! physically, this situation reflects an approximate character of the lake model.
             ! using the self-similar temperature profile in the thermocline, 
             ! there is always communication between the mixed layer, the thermocline 
             ! and the lake bottom. as a result, the rate of change of the cbl depth
             ! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
             ! in reality, convective mixed-layer deepening may be completely decoupled
             ! from the processes underneath. in order to account for this fact,
             ! the rate of cbl deepening is set to a small value
             ! if dh_ml/dt proves to be negative.
             ! this is "double insurance" however, 
             ! as a negative dh_ml/dt is encountered very rarely.
             !_dm

             !_dbg
             if(d_h_ml_dt<0.0_ireals) then 
             end if

             d_h_ml_dt = max(d_h_ml_dt, c_small_flk)    
             h_ml_n_flk = h_ml_p_flk + d_h_ml_dt*del_time                       ! update h_ml 
             h_ml_n_flk = max(h_ml_min_flk, min(h_ml_n_flk, depth_w))           ! security, limit h_ml
          else                                              ! mixing down to the lake bottom
             h_ml_n_flk = depth_w
          end if

       else mixing_regime                              ! wind mixing


          checkout=(i_w_flk + i_h_flk - 2.0_ireals*i_intm_0_h_flk)
          d_h_ml_dt = max(u_star_w_flk, u_star_min_flk)                        ! the surface friction velocity
          zm_h_scale = (abs(par_coriolis)/c_sbl_zm_n + n_t_mean/c_sbl_zm_i)*d_h_ml_dt**2.0_ireals
          zm_h_scale = zm_h_scale + flk_str_1/c_sbl_zm_s
          zm_h_scale = max(zm_h_scale, c_small_flk)
          zm_h_scale = d_h_ml_dt**3.0_ireals/zm_h_scale 
          zm_h_scale = max(h_ml_min_flk, min(zm_h_scale, h_ml_max_flk))        ! the zm96 sbl depth scale 
          zm_h_scale = max(zm_h_scale, conv_equil_h_scale)                     ! equilibrium mixed-layer depth 

          !_dm 
          ! in order to avoid numerical discretization problems,
          ! an analytical solution to the evolution equation 
          ! for the wind-mixed layer depth is used.
          ! that is, an exponential relaxation formula is applied
          ! over the time interval equal to the model time step.
          !_dm 

          d_h_ml_dt = c_relax_h*d_h_ml_dt/zm_h_scale*del_time
          h_ml_n_flk = zm_h_scale - (zm_h_scale-h_ml_p_flk)*exp(-d_h_ml_dt)    ! update h_ml 
          h_ml_n_flk = max(h_ml_min_flk, min(h_ml_n_flk, depth_w))             ! limit h_ml 
          d_h_ml_dt = (h_ml_n_flk-h_ml_p_flk)/del_time                         ! re-compute dh_ml/dt

          if(h_ml_n_flk<=h_ml_p_flk)           &
               d_c_t_dt = -d_c_t_dt                 ! mixed-layer retreat or stationary state, dc_t/dt<0
          c_t_n_flk = c_t_p_flk + d_c_t_dt*del_time                            ! update c_t
          c_t_n_flk = min(c_t_max, max(c_t_n_flk, c_t_min))                    ! limit c_t 
          d_c_t_dt = (c_t_n_flk-c_t_p_flk)/del_time                            ! re-compute dc_t/dt

       end if mixing_regime

       ! compute the time-rate-of-change of the the bottom temperature, 
       ! depending on the sign of dh_ml/dt 
       ! update the bottom temperature and the mixed-layer temperature

       if(h_ml_n_flk<=depth_w-h_ml_min_flk) then       ! mixing did not reach the bottom 

          if(h_ml_n_flk>h_ml_p_flk) then   ! mixed-layer deepening 
             r_h_icesnow     = h_ml_p_flk/depth_w
             r_rho_c_icesnow = 1.0_ireals-r_h_icesnow 
             r_ti_icesnow    = 0.5_ireals*c_t_p_flk*r_rho_c_icesnow+c_tt_flk*(2.0_ireals*r_h_icesnow-1.0_ireals)
             r_tstar_icesnow = (0.5_ireals+c_tt_flk-c_q_flk)/r_ti_icesnow
             r_ti_icesnow    = (1.0_ireals-c_t_p_flk*r_rho_c_icesnow)/r_ti_icesnow

             d_t_bot_dt = (q_w_flk-q_bot_flk+i_w_flk-i_bot_flk)/tpl_rho_w_r/tpl_c_w
             d_t_bot_dt = d_t_bot_dt - c_t_p_flk*(t_wml_p_flk-t_bot_p_flk)*d_h_ml_dt
             d_t_bot_dt = d_t_bot_dt*r_tstar_icesnow/depth_w                   ! q+i fluxes and dh_ml/dt term

             flk_str_2 = i_intm_h_d_flk - (1.0_ireals-c_q_flk)*i_h_flk - c_q_flk*i_bot_flk
             flk_str_2 = flk_str_2*r_ti_icesnow/(depth_w-h_ml_p_flk)/tpl_rho_w_r/tpl_c_w
             d_t_bot_dt = d_t_bot_dt + flk_str_2                               ! add radiation-flux term

             flk_str_2 = (1.0_ireals-c_tt_2*r_ti_icesnow)/c_t_p_flk
             flk_str_2 = flk_str_2*(t_wml_p_flk-t_bot_p_flk)*d_c_t_dt
             d_t_bot_dt = d_t_bot_dt + flk_str_2                               ! add dc_t/dt term

          else                                ! mixed-layer retreat or stationary state
             d_t_bot_dt = 0.0_ireals                                            ! dt_bot/dt=0
          end if

          t_bot_n_flk = t_bot_p_flk + d_t_bot_dt*del_time                      ! update t_bot  
          t_bot_n_flk = max(t_bot_n_flk, tpl_t_f)           ! security, limit t_bot by the freezing point
          flk_str_2 = (t_bot_n_flk-tpl_t_r)*flake_buoypar(t_mnw_n_flk)
          if(flk_str_2<0.0_ireals) t_bot_n_flk = tpl_t_r  ! security, avoid t_r crossover 
          t_wml_n_flk = c_t_n_flk*(1.0_ireals-h_ml_n_flk/depth_w)
          t_wml_n_flk = (t_mnw_n_flk-t_bot_n_flk*t_wml_n_flk)/(1.0_ireals-t_wml_n_flk)
          t_wml_n_flk = max(t_wml_n_flk, tpl_t_f)           ! security, limit t_wml by the freezing point
       else                                              ! mixing down to the lake bottom 
          h_ml_n_flk = depth_w
          t_wml_n_flk = t_mnw_n_flk
          t_bot_n_flk = t_mnw_n_flk
          c_t_n_flk = c_t_min
       end if

    end if htc_water

    !  compute the depth of the thermally active layer of bottom sediments
    !  and the temperature at this depth.

    use_sediment: if(lflk_botsed_use) then   ! the bottom-sediment scheme is used

       if(h_b1_p_flk>=depth_bs-h_b1_min_flk) then   ! no t(z) maximum (no thermal wave) 
          h_b1_p_flk = 0.0_ireals                       ! set h_b1_p to zero
          t_b1_p_flk = t_bot_p_flk                     ! set t_b1_p to the bottom temperature
       end if

       flk_str_1 = 2.0_ireals*phi_b1_pr0/(1.0_ireals-c_b1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
       h_ice_threshold = sqrt(flk_str_1)                              ! threshold value of h_b1
       h_ice_threshold = min(0.9_ireals*depth_bs, h_ice_threshold)    ! limit h_b1
       flk_str_2 = c_b2/(1.0_ireals-c_b2)*(t_bs-t_b1_p_flk)/(depth_bs-h_b1_p_flk)

       if(h_b1_p_flk<h_ice_threshold) then  ! use a truncated equation for h_b1(t)
          h_b1_n_flk = sqrt(h_b1_p_flk**2.0_ireals+flk_str_1)  ! advance h_b1
          d_h_b1_dt = (h_b1_n_flk-h_b1_p_flk)/del_time          ! re-compute dh_b1/dt
       else                                    ! use a full equation for h_b1(t)
          flk_str_1 = (q_bot_flk+i_bot_flk)/h_b1_p_flk/tpl_rho_w_r/tpl_c_w
          flk_str_1 = flk_str_1 - (1.0_ireals-c_b1)*(t_bot_n_flk-t_bot_p_flk)/del_time
          d_h_b1_dt = (1.0_ireals-c_b1)*(t_bot_p_flk-t_b1_p_flk)/h_b1_p_flk + c_b1*flk_str_2

          if(abs(d_h_b1_dt)<1.e-15_ireals)then
             write(*,*) 'flake_driver, bottom sediments: deviding by zero' 
          endif
          d_h_b1_dt = flk_str_1/d_h_b1_dt
          h_b1_n_flk = h_b1_p_flk + d_h_b1_dt*del_time          ! advance h_b1
       end if
       d_t_b1_dt = flk_str_2*d_h_b1_dt
       t_b1_n_flk = t_b1_p_flk + d_t_b1_dt*del_time            ! advance t_b1

       l_snow_exists = h_b1_n_flk>=depth_bs-h_b1_min_flk                    & ! h_b1 reached depth_bs, or
            .or. h_b1_n_flk<h_b1_min_flk                             & ! h_b1 decreased to zero, or
            .or.(t_bot_n_flk-t_b1_n_flk)*(t_bs-t_b1_n_flk)<=0.0_ireals   ! there is no t(z) maximum
       if(l_snow_exists) then      
          h_b1_n_flk = depth_bs                     ! set h_b1 to the depth of the thermally active layer
          t_b1_n_flk = t_bs                         ! set t_b1 to the climatological temperature 
       end if

    else use_sediment                        ! the bottom-sediment scheme is not used

       h_b1_n_flk = rflk_depth_bs_ref              ! h_b1 is set to a reference value 
       t_b1_n_flk = tpl_t_r                        ! t_b1 is set to the temperature of maximum density

    end if use_sediment

    !  impose additional constraints.

    ! in case of unstable stratification, force mixing down to the bottom
    flk_str_2 = (t_wml_n_flk-t_bot_n_flk)*flake_buoypar(t_mnw_n_flk)
    if(flk_str_2<0.0_ireals) then 
       h_ml_n_flk = depth_w
       t_wml_n_flk = t_mnw_n_flk
       t_bot_n_flk = t_mnw_n_flk
       c_t_n_flk = c_t_min

    end if


    !  update the surface temperature.

    if(h_snow_n_flk>=h_snow_min_flk) then   
       t_sfc_n = t_snow_n_flk                   ! snow exists, use the snow temperature
    else if(h_ice_n_flk>=h_ice_min_flk) then
       t_sfc_n = t_ice_n_flk                    ! ice exists but there is no snow, use the ice temperature
    else 
       t_sfc_n = t_wml_n_flk                    ! no ice-snow cover, use the mixed-layer temperature
    end if

    checkout = q_bot_flk

  end subroutine flake_driver

  real (kind = ireals) function flake_buoypar (t_water)

    ! Description:
    !
    !  Computes the buoyancy parameter,
    !  using a quadratic equation of state for the fresh-water.



    implicit none

    real (kind = ireals), intent(in) :: t_water! water temperature [k]

    flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)


  end function flake_buoypar

  real(kind = ireals) function flake_snowdensity (h_snow)

    ! Description:
    !
    !  Computes the snow density, using an empirical approximation from Heise et al. (2003).

    implicit none

    real (kind = ireals), intent(in) :: h_snow                              ! snow thickness [m]

    !  security. ensure that the expression in () does not become negative at a very large h_snow.
    flake_snowdensity = max( c_small_flk, (1.0_ireals - h_snow*tpl_gamma_rho_s/tpl_rho_w_r) )
    flake_snowdensity = min( tpl_rho_s_max, tpl_rho_s_min/flake_snowdensity )

  end function flake_snowdensity
  real (kind = ireals) function flake_snowheatconduct (h_snow)

    ! Description:
    !  Computes the snow heat conductivity,
    !  using an empirical approximation from Heise et al. (2003).



    implicit none

    real (kind = ireals), intent(in) :: h_snow                              ! snow thickness [m]

    flake_snowheatconduct = flake_snowdensity( h_snow )   ! Compute snow density
    flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
         + h_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )


  end function flake_snowheatconduct

  subroutine slfluxo_surf_lake_ice(nhor,nlev,kstart,kstop, &
       dtime, coslat,lcounttice, &
       t,q,u,v,gpot,dph, &
       ps, &
       tlakesfc, tlakesnow, tlakeice, tlakemean, &
       tlakeml, tlakebot, ctlake, sndepthlake, &
       icedepthlake, mldepth, &
       t_b1lake, h_b1lake, &
       frlake, &
       depthlake, &
       tskin, & 
       draindt, dsnowdt, &
       radf, &
       albedo, &
       scos, &
       sswdn, &
       emskin, &
       z0lake, &     
       t2ms,t2mi, &
       q2ms,q2mi, &
       rh2ms,rh2mi, &
       u10ms,u10mi, &
       v10ms,v10mi, &
       senfs,senfi, &
       latfs,latfi, &
       evaps,evapi, &
       accrunofflake, &
       momfus,momfui, &
       momfvs,momfvi, &
       ustars,ustari, &
       rousea, &
       tsnice,frsnice, &
       t2mlake, q2mlake, rh2mlake, &
       u10mlake, v10mlake, &
       momfulake, momfvlake , &
       senflake, latflake, &     
       sswrlake, slwrlake, &
       ustarlake, &
       tlakesfctend, tlakesnowtend, tlakeicetend, tlakemeantend, &
       tlakemltend, tlakebottend, ctlaketend, sndepthlaketend, &
       icedepthlaketend, mldepthtend, &
       t_b1laketend, h_b1laketend)


    !      "common" atmospheric variables -----------------------------
    !cl    nhor:    number of points in  horizontal loop.
    !cl    nlev:    number of vertical levels.
    !cl    kstart:  start index of horizontal loop.
    !cl    kstop:   stop index of horizontal loop.
    !cl    dtime:   time step (s),(two time the time step in the dynamics)
    !      coslat:  cosine if latitude
    !cl    t    :   temperature (k),(3d field)
    !cl    q    :   specific humidity ,(kg/kg),(3d field)
    !cl    u    :   wind component towards the east ,(m/s), (3d field)
    !cl    v    :   wind component towards the north,(m/s), (3d field)

    !cl    gpot :   geopotential height (m2/s2)
    !cl    dph  :   pressure difference (pa) between half levels.
    !cl    ps   :   surface pressure (pa)
    !c   draindt:   rain intencity
    !c   dsnowdt:   snow intencity
    !cl    scos:    cosine for solar zenith angle
    !cl    sswdn:   global radiation
    !      surface (lake) specific variables ---------------------------
    !      tlakesfc  :   lake surface temperature, 
    !                    including ice and snow (pseudo-prognostic) (k)
    !      tlakesnow :  temperature of snow over the lake ice (k)
    !      tlakeice  :  temperature of lake ice (k)
    !      tlakemean :  mean water lake temperature (k)
    !      tlakeml   : mixed layer lake temperature (k)
    !      tlakebot : lake bottom temperature (k)
    !      ctlake   : chape-factor 
    !      sndepthlake : snow depth over lake ice (m)
    !      icedepthlake : lake ice depth (m)
    !      mldepth : mixed layer depth (m)
    !      t_b1lake : temperature on the extremum on bottom sedoments (k)
    !               : when bottom sediments = off, use dummy value
    !      h_b1lake: depth of the bottom sediments temperature extremum (m)
    !              : when bottom sediments = off, use dummy value
    !      frlake : farction of lakes of each type
    !      depthlake : depth of lake of each type, (m)
    !      surface (not lake) specific variables -----------------------
    !      tskin:   skin surface temperature
    !cl    radf:    grid square average of surface radiation budget
    !1cl   albedo:  grid square average of surface albedo
    !cl    emskin:  grid average emissivity
    !cl    z0lake:  lake surface roughness (m) 
    !      t2mlake: 2 meter temperature (k) over lake (with/without ice)
    !      q2mlake: 2 mer specific humidity (kg/kg) over lake (with/without ice)
    !      rh2mlake: 2 meter relative humidity over lake (with/without ice)
    !cl    u10lake: 10 meter u-wind over lake (with/without ice)
    !cl    v10lake: 10 meter v-wind over lake (with/without ice)
    !      momfulake: momentum flux towards the east (kg/ (ms2)) over lake (with/without ice)
    !      momfvlake: momentum flux towards the north (kg/ (ms2)) over lake (with/without ice)
    !      senflake: sensible heat flux over lake (with/without ice)
    !      latflake : latent heat flux over lake (with/without ice)
    !      sswrlake : surface short wave radiation over lake (with/without ice)
    !      slwrlake : surface long wave radiation over lake (with/without ice)
    !      ustarlake : surface friction velocity (m/s) for lake (with/without ice)
    !      tlakesfctend  : lake surface temperature, including ice and snow 
    !                      (pseudo-prognostic) tendency
    !      tlakesnowtend :  temperature of snow over the lake ice tendency
    !      tlakeicetend  :  temperature of lake ice tendency
    !      tlakemeantend :  mean water lake temperature tendency
    !      tlakemltend   : mixed layer lake temperature tendency
    !      tlakebottend : lake bottom temperature tendency
    !      ctlaketend   : chape-factor tendency
    !      sndepthlaketend : snow depth over lake ice tendency
    !      icedepthlaketend : lake ice depth tendency
    !      mldepthtend : mixed layer depth tendency
    !      t_b1laketend : temperature on the extremum of bottom sediments tendency
    !      h_b1laketend : depth of the tempersature extremum of bottom sediments tendency

    !c     note: latfs,senfs > 0 for fluxes towards surface.


    
    use confys !, only:epsilo,tmelt
    implicit none

    integer, intent (in) :: nhor,nlev,kstart,kstop

    real(kind=realkind), intent (in) :: dtime
    real(kind=realkind), dimension(nhor),intent (in) :: coslat,lcounttice
    real(kind=realkind), dimension(nhor,nlev), intent (in) ::  gpot,q,t,u,v
    real(kind=realkind), dimension(nhor,nlev+1), intent (in) :: dph
    real(kind=realkind), dimension(nhor), intent (in) :: ps, tskin, draindt, dsnowdt, &
         radf, albedo, scos, sswdn, emskin
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: tlakesfc, tlakesnow, &
         tlakeice, tlakemean, & 
         tlakeml, tlakebot, &
         ctlake, sndepthlake, &
         icedepthlake, mldepth, &
         t_b1lake, h_b1lake, &
         frlake, depthlake

    real(kind=realkind), dimension(nhor, nlaketype), intent (inout) :: z0lake 

    real(kind=realkind), dimension(nhor), intent (inout) :: t2ms,t2mi, &
         q2ms,q2mi,rh2ms,rh2mi,u10ms,u10mi,v10ms,v10mi, &
         senfs,senfi,latfs,latfi,evaps,evapi,accrunofflake, &
         momfus,momfui,momfvs,momfvi,ustars,ustari,rousea,tsnice,frsnice

    real(kind=realkind), dimension(nhor, nlaketype), intent (out) :: t2mlake, q2mlake, rh2mlake, &
         u10mlake, v10mlake, momfulake, momfvlake, senflake, latflake, &
         sswrlake, slwrlake, ustarlake

    real(kind=realkind), dimension(nhor, nlaketype), intent (out) :: tlakesfctend, tlakesnowtend, &
         tlakeicetend, tlakemeantend, & 
         tlakemltend, tlakebottend, &
         ctlaketend, sndepthlaketend, &
         icedepthlaketend, mldepthtend, &
         t_b1laketend, h_b1laketend

    real(kind=realkind) zfrsum

    ! locals
    integer ihor, ilaketype ! loop indicies

    real(kind=realkind), dimension(nhor, nlaketype) :: evaplake


    do ilaketype=1,nlaketype ! put zero to all output values 
       do ihor=kstart,kstop 
          u10mlake(ihor,ilaketype)=0.0_realkind 
          v10mlake(ihor,ilaketype)=0.0_realkind  
          momfulake(ihor,ilaketype)=0.0_realkind  
          momfvlake(ihor,ilaketype)=0.0_realkind  
          senflake(ihor,ilaketype)=0.0_realkind  
          latflake(ihor,ilaketype)=0.0_realkind 
          evaplake(ihor,ilaketype)=0.0_realkind 
          sswrlake(ihor,ilaketype)=0.0_realkind 
          slwrlake(ihor,ilaketype)=0.0_realkind 
          ustarlake(ihor,ilaketype)=0.0_realkind 
          tlakeicetend(ihor,ilaketype)=0.0_realkind 
          tlakemeantend(ihor,ilaketype)=0.0_realkind 
          tlakemltend(ihor,ilaketype)=0.0_realkind 
          tlakebottend(ihor,ilaketype)=0.0_realkind 
          ctlaketend(ihor,ilaketype)=0.0_realkind 
          sndepthlaketend(ihor,ilaketype)=0.0_realkind 
          icedepthlaketend(ihor,ilaketype)=0.0_realkind 
          mldepthtend(ihor,ilaketype)=0.0_realkind     
          t_b1laketend(ihor,ilaketype) = 0.0_realkind 
          h_b1laketend(ihor,ilaketype) = 0.0_realkind 
          tlakesfctend(ihor,ilaketype) = 0.0_realkind 
       end do
    end do ! of put zero to all output values loop 

    ! roughness, turbulent fluxes, diagnostic variables 
    ! over lakes free of ice
    call turbfluxlakewater(nhor,nlev,kstart,kstop, &
         dtime, &
         t,q,u,v,gpot,dph, &
         ps, &
         tlakesfc, &
         sndepthlake, &
         icedepthlake, &
         frlake, &
         z0lake, &     
         t2mlake, q2mlake, &
         u10mlake, v10mlake, &
         momfulake, momfvlake , &
         senflake, latflake, evaplake, &     
         ustarlake)


    ! roughness, turbulent fluxes, diagnostic variables 
    ! over lakes with ice and/or snow
    call turbfluxlakeicesnow(nhor,nlev,kstart,kstop, &
         dtime, &
         t,q,u,v,gpot,dph, &
         ps, &
         tlakesfc, &
         sndepthlake, &
         icedepthlake, &
         frlake, &    
         t2mlake, q2mlake, &
         u10mlake, v10mlake, &
         momfulake, momfvlake , &
         senflake, latflake, evaplake, &     
         ustarlake, &
         coslat)

    do ilaketype=1,nlaketype 
       do ihor=kstart,kstop
          ! relative humidity
          if(q2mlake(ihor,ilaketype)>0.0_realkind ) then
             rh2mlake(ihor,ilaketype)=q2e(ps(ihor),q2mlake(ihor,ilaketype))/ &
                  t2es(t2mlake(ihor,ilaketype))
             rh2mlake(ihor,ilaketype)=min(max(rh2mlake(ihor,ilaketype),0.0_realkind ),1.0_realkind )
          else
             rh2mlake(ihor,ilaketype)=0.0_realkind 
          endif
       enddo
    enddo

    do ihor=kstart,kstop
       if(nint(lcounttice(ihor))==3)then
          t2ms(ihor)=0.0_realkind 
          q2ms(ihor)=0.0_realkind 
          rh2ms(ihor)=0.0_realkind 
          u10ms(ihor)=0.0_realkind 
          v10ms(ihor)=0.0_realkind 
          senfs(ihor)=0.0_realkind 
          latfs(ihor)=0.0_realkind 
          evaps(ihor)=0.0_realkind 
          momfus(ihor)=0.0_realkind 
          momfvs(ihor)=0.0_realkind 
          ustars(ihor)=0.0_realkind 
          rousea(ihor)=0.0_realkind 
          tsnice(ihor)=0.0_realkind 
          frsnice(ihor)=0.0_realkind 
          !zfrsum=0.0_realkind 
          zfrsum = sum(frlake(ihor,:))
          do ilaketype=1,nlaketype 
             !             if(frlake(ihor,ilaketype)>=fraclakemin)then !Detta är sant by construction!
             t2ms(ihor)=t2ms(ihor)+frlake(ihor,ilaketype)*t2mlake(ihor,ilaketype)
             q2ms(ihor)=q2ms(ihor)+frlake(ihor,ilaketype)*q2mlake(ihor,ilaketype)
             rh2ms(ihor)=rh2ms(ihor)+frlake(ihor,ilaketype)*rh2mlake(ihor,ilaketype)
             u10ms(ihor)=u10ms(ihor)+frlake(ihor,ilaketype)*u10mlake(ihor,ilaketype)
             v10ms(ihor)=v10ms(ihor)+frlake(ihor,ilaketype)*v10mlake(ihor,ilaketype)
             senfs(ihor)=senfs(ihor)+frlake(ihor,ilaketype)*senflake(ihor,ilaketype)
             latfs(ihor)=latfs(ihor)+frlake(ihor,ilaketype)*latflake(ihor,ilaketype)
             evaps(ihor)=evaps(ihor)+frlake(ihor,ilaketype)*evaplake(ihor,ilaketype)
             momfus(ihor)=momfus(ihor)+frlake(ihor,ilaketype)*momfulake(ihor,ilaketype)
             momfvs(ihor)=momfvs(ihor)+frlake(ihor,ilaketype)*momfvlake(ihor,ilaketype)
             ustars(ihor)=ustars(ihor)+frlake(ihor,ilaketype)*ustarlake(ihor,ilaketype)
             rousea(ihor)=rousea(ihor)+frlake(ihor,ilaketype)*z0lake(ihor,ilaketype)
             tsnice(ihor)=tsnice(ihor)+frlake(ihor,ilaketype)*tlakesfc(ihor,ilaketype)
             frsnice(ihor)=frsnice(ihor)+frlake(ihor,ilaketype)*sndepthlake(ihor,ilaketype)
             !                zfrsum=zfrsum+frlake(ihor,ilaketype)
             !             endif
          enddo
          t2ms(ihor)=t2ms(ihor)/zfrsum
          q2ms(ihor)=q2ms(ihor)/zfrsum
          rh2ms(ihor)=rh2ms(ihor)/zfrsum
          u10ms(ihor)=u10ms(ihor)/zfrsum
          v10ms(ihor)=v10ms(ihor)/zfrsum
          senfs(ihor)=senfs(ihor)/zfrsum
          latfs(ihor)=latfs(ihor)/zfrsum
          evaps(ihor)=evaps(ihor)/zfrsum
          momfus(ihor)=momfus(ihor)/zfrsum
          momfvs(ihor)=momfvs(ihor)/zfrsum
          ustars(ihor)=ustars(ihor)/zfrsum
          rousea(ihor)=rousea(ihor)/zfrsum
          tsnice(ihor)=tsnice(ihor)/zfrsum
          if(frsnice(ihor)>0.0_realkind )then
             frsnice(ihor)=1.0_realkind 
          endif
          !
          t2mi(ihor)=t2ms(ihor)
          q2mi(ihor)=q2ms(ihor)
          rh2mi(ihor)=rh2ms(ihor)
          u10mi(ihor)=u10ms(ihor)
          v10mi(ihor)=v10ms(ihor)
          senfi(ihor)=senfs(ihor)
          latfi(ihor)=latfs(ihor)
          evapi(ihor)=evaps(ihor)
          momfui(ihor)=momfus(ihor)
          momfvi(ihor)=momfvs(ihor)
          ustari(ihor)=ustars(ihor)
          !
          accrunofflake(ihor)=accrunofflake(ihor)+ &
               ((draindt(ihor)+dsnowdt(ihor))*dtime-evaps(ihor))*zfrsum
       endif ! lcounttice=3
    end do

    ! radiation fluxes over lakes
    call radsurflake(nhor,kstart,kstop, & 
         tskin, radf, albedo, scos, & 
         sswdn, emskin, & 
         tlakesfc,frlake,icedepthlake, sndepthlake, &  
         sswrlake, slwrlake )

    ! interface to flake - lake model
    call flakeinterface3d(nhor,kstart,kstop, &
         dtime, coslat, &
         momfulake, momfvlake , senflake, latflake, &
         sswrlake, slwrlake,   &
         tlakesfc, tlakesnow, tlakeice, tlakemean, &
         tlakeml, tlakebot, ctlake, sndepthlake, &
         icedepthlake, mldepth, &
         t_b1lake, h_b1lake, &
         frlake, depthlake, &
         tlakesfctend, tlakesnowtend, tlakeicetend, tlakemeantend, &
         tlakemltend, tlakebottend, ctlaketend, sndepthlaketend, &
         icedepthlaketend, mldepthtend, &
         t_b1laketend, h_b1laketend)

    return 
  end subroutine slfluxo_surf_lake_ice





  subroutine radsurflake(nhor,kstart,kstop, & 
       tskin, radf, albedo, scos, & 
       sswdn, emskin, & 
       tlakesfc,frlake,icedepthlake, sndepthlake, & 
       sswrlake, slwrlake ) 


    !cl    input:

    !
    !      "common" atmospheric variables -----------------------------
    !cl    nhor:    number of points in  horizontal loop.
    !cl    kstart:  start index of horizontal loop.
    !cl    kstop:   stop index of horizontal loop.
    !cl    scos:    cosine for solar zenith angle
    !      surface (not lake) specific variables -----------------------
    !      tskin:   skin surface temperature
    !cl    radf:    grid square average of surface radiation budget
    !1cl   albedo:  grid square average of surface albedo
    !cl    emskin:  grid average emissivity
    !cl    sswdn:   global radiation
    !      surface (lake) specific variables ----------------------------
    !      tlakesfc  :   lake surface temperature, 
    !                    including ice and snow (pseudo-prognostic) (k)
    !      frlake : farction of lakes of each type
    !      icedepthlake : lake ice depth (m)
    !      sndepthlake : snow depth over lake ice (m)

    !cl    output

    !      sswrlake : surface short wave radiation over lake (with/without ice)
    !      slwrlake : surface long wave radiation over lake (with/without ice)
    ! uses --------------------------------------------------------------------------
    use confys
    use comrpar
    implicit none

    ! exchange variables
    integer, intent (in) :: nhor,kstart,kstop
    real(kind=realkind), dimension(nhor), intent (in) ::  tskin, &
         radf, albedo, scos, sswdn, emskin        
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: tlakesfc, &
         icedepthlake, sndepthlake, frlake
    real(kind=realkind), dimension(nhor, nlaketype), intent (out) ::sswrlake, slwrlake

    ! local variables 

    integer ihor, ilaketype ! loop indicies

    real(kind=realkind) zscos, zalb, zrads, zradf, zradl, zalblake, emlake
    real(kind=realkind) zlakewateron, zlakewateroff, &
         zlakeiceon, zlakeiceoff, &
         zlakesnowon, zlakesnowoff
    real(kind=realkind) zalbicenew,zalbice0


    do ilaketype=1,nlaketype 
       do ihor=kstart,kstop 

          ! switchers 
          if(frlake(ihor,ilaketype)<fraclakemin) then ! fraction of lake is too small 
             zlakewateron = 0._realkind
             zlakewateroff = 1._realkind
             zlakeiceon = 0._realkind
             zlakeiceoff = 1._realkind
             zlakesnowon = 0._realkind
             zlakesnowoff = 1._realkind 
          else ! fraction of lake is big enough
             if(icedepthlake(ihor,ilaketype)<h_ice_min_flk) then ! there is no ice
                zlakewateron = 1._realkind
                zlakewateroff = 0._realkind
                zlakeiceon = 0._realkind
                zlakeiceoff = 1._realkind
                zlakesnowon = 0._realkind
                zlakesnowoff = 1._realkind  
             else ! there is ice
                if(sndepthlake(ihor,ilaketype)<h_snow_min_flk) then ! there is no snow on ice
                   zlakewateron = 1._realkind
                   zlakewateroff = 0._realkind
                   zlakeiceon = 1._realkind
                   zlakeiceoff = 0._realkind
                   zlakesnowon = 0._realkind
                   zlakesnowoff = 1._realkind
                else ! there is snow on ice
                   zlakewateron = 1._realkind
                   zlakewateroff = 0._realkind
                   zlakeiceon = 1._realkind
                   zlakeiceoff = 0._realkind
                   zlakesnowon = 1._realkind
                   zlakesnowoff = 0._realkind
                end if ! snow condition
             end if ! ice condition
          end if ! fraction of lake condition

          zscos=0.2_realkind /(1.0_realkind +scos(ihor)) - 0.1_realkind 
          zalb = albedo(ihor) + zscos
          zrads=sswdn(ihor)
          zradf = radf(ihor)-zrads*(1.0_realkind -zalb)
          zalbice0 = exp(-c_albice_mr*(tmelt - tlakesfc(ihor,ilaketype))/tmelt)! sd, 16.03.05      
          zalbicenew = albedo_whiteice_ref*(1.0_realkind  - zalbice0) + albedo_blueice_ref*zalbice0 ! sd, 16.03.05  
          zalblake = zlakewateron*zlakeiceoff*zlakesnowoff*albwater + &
               zlakewateron*zlakeiceon*zlakesnowoff*zalbicenew + &  ! 16.03.05, kk, sd
               zlakewateron*zlakeiceon*zlakesnowon*albsnow +  zscos 
          sswrlake(ihor,ilaketype) = (1.0_realkind -zalblake)*zrads*zlakewateron

          !c  zradl is downward long wave radiation (common for all surfaces)
          zradl = zradf/emskin(ihor) + stebol*tskin(ihor)**4.0_realkind
          emlake = zlakewateron*zlakeiceoff*zlakesnowoff*emwater + &
               zlakewateron*zlakeiceon*zlakesnowoff*emice + &
               zlakewateron*zlakeiceon*zlakesnowon*emsnow
          slwrlake(ihor,ilaketype)=emlake* &
               (zradl-stebol*tlakesfc(ihor,ilaketype)**4.0_realkind)*zlakewateron      

       enddo ! loop over hor. grid points
    enddo ! loop over lake types


    return
  end subroutine radsurflake


  subroutine lakemasks(klon,klat,meco,eco,lcounttice)


    implicit none
    ! forms mask lake/sea  and corrects fractions
    ! not the best algorythm at the moment

    ! exchange variables 
    integer, intent(in) :: klon , klat,meco
    real(kind=realkind),target,dimension(klon,klat,meco), intent(inout) :: eco
    real(kind=realkind),dimension(klon,klat), intent(inout) :: lcounttice
    real(kind=realkind),pointer,dimension(:,:,:)::frlake,depthlake
    real(kind=realkind),pointer,dimension(:,:)::frland,fr_lake

    ! nb!!! lcounttice used here twohold: as a lake-water/sea-water mask and a 
    !       switcher between lake models: 1 - lake is treated as sea
    !                                     2 - probe
    !                                     3 - flake
    !       as the main destination of lcounttice - to be a mask, we use 
    !       environmental variables flakeflag and multilakeflag additionaly 
    !       if multilakeflag is not defined  => we can't use probe and flake in parallel
    !       if flakeflag not defined, we use (sea) or probe, if is defined by lcounttice
    !       if flakeflag is defined, we use flake, and change lcounttice = 3 in appr. points

    ! local variables 
    real(kind=realkind) zfrdiff,zfrsea,zfrtemp,zfrlaketot
    integer i, j, k ! loop indexes


    frlake => eco(:,:,44:46)
    depthlake => eco(:,:,47:49)
    frland => eco(:,:,10)
    fr_lake => eco(:,:,42)

    lcounttice=1.0_realkind 

    do i=1,klon
       do j=1,klat
          ! put frland=0 if frland<0.01
          if(frland(i,j)<fraclakemin.and.frland(i,j)>0.0_realkind )then
             if(fr_lake(i,j)>=fraclakemin)then
                fr_lake(i,j)=fr_lake(i,j)+frland(i,j)
             else
                fr_lake(i,j)=0.0_realkind 
                do k = 1, nlaketype
                   frlake(i,j,k)=0.0_realkind 
                enddo
             endif
             frland(i,j)=0.0_realkind 
          endif

          ! put all lake-water to sea if zfrsea>=0.01
          zfrsea=1.0_realkind -frland(i,j)-fr_lake(i,j)
          if(zfrsea>=fraclakemin)then
             if(fr_lake(i,j)>0.0_realkind )then
                zfrsea=zfrsea+fr_lake(i,j)
                fr_lake(i,j)=0.0_realkind 
                frland(i,j)=1.0_realkind -zfrsea
             endif
          endif

          ! put lake-water to land if fr_lake<0.01
          if(fr_lake(i,j)<fraclakemin.and.fr_lake(i,j)>0.0_realkind )then
             frland(i,j)=frland(i,j)+fr_lake(i,j)
             fr_lake(i,j)=0.0_realkind 
             do k = 1, nlaketype
                frlake(i,j,k)=0.0_realkind 
             enddo
          endif

          ! check if svars-lake are big enough (>=0.01) and sum them together
          zfrlaketot=0.0_realkind 
          do k = 1, nlaketype
             if (frlake(i,j,k)>=fraclakemin) then
                zfrlaketot=zfrlaketot+frlake(i,j,k)
             else
                frlake(i,j,k)=0.0_realkind 
             endif
          enddo

          ! the total fraction of svars-lakes are >=0.01
          ! correct for difference between area of svars-lakes and ecoclimap lakes
          if (zfrlaketot>=fraclakemin) then
             zfrdiff=fr_lake(i,j)-zfrlaketot
             ! there are less svars-lakes than ecoclimap lakes.
             ! thus, increase the deep svars-lakes
             if(zfrdiff>0.0_realkind )then
                frlake(i,j,1)=frlake(i,j,1)+zfrdiff
                depthlake(i,j,1)=max(10.0_realkind ,depthlake(i,j,1))
                ! there are more svars-lakes than ecoclimap lakes.
                ! thus, decrease the svars-lakes
             elseif(zfrdiff<0.0_realkind )then
                do k = 1,nlaketype
                   if(frlake(i,j,k)>0.0_realkind )then
                      zfrtemp=frlake(i,j,k)+zfrdiff
                      frlake(i,j,k)=max(zfrtemp,0.0_realkind )
                      zfrdiff=min(zfrtemp,0.0_realkind )
                   endif
                enddo
             endif
             ! the total fraction of svars-lakes are <0.01
             ! check if there are enough ecoclimap lakes
          else
             if(fr_lake(i,j)>=fraclakemin)then
                frlake(i,j,1)=fr_lake(i,j)
                depthlake(i,j,1)=max(10.0_realkind ,depthlake(i,j,1))
             else
                frlake(i,j,1)=0.0_realkind 
             endif
             do k = 2, nlaketype
                frlake(i,j,k)=0.0_realkind 
             enddo
          endif

          ! check again if svars-lake are big enough (>=0.01) and sum them together
          zfrlaketot=0.0_realkind 

          do k = nlaketype,2,-1
             zfrtemp=0.0_realkind 
             if(frlake(i,j,k)<fraclakemin)then
                zfrtemp=frlake(i,j,k)
                frlake(i,j,k)=0.0_realkind 
             else
                zfrlaketot=zfrlaketot+frlake(i,j,k)
             endif
             frlake(i,j,k-1)=frlake(i,j,k-1)+zfrtemp
          enddo

          if(frlake(i,j,1)<fraclakemin)then
             frlake(i,j,1)=0.0_realkind 
          else
             zfrlaketot=zfrlaketot+frlake(i,j,1)
          endif

          ! set lcounttice=3. if the total fraction of lakes are >=0.01 
          zfrtemp=fr_lake(i,j)-zfrlaketot
          if(zfrlaketot>=fraclakemin)then
             lcounttice(i,j)=3.0_realkind 
             fr_lake(i,j)=zfrlaketot
             frland(i,j)=frland(i,j)+zfrtemp
             !  put the lake depth for deep lakes to maximum 40 m
             depthlake(i,j,1)=min(40.0_realkind ,depthlake(i,j,1))
          else
             fr_lake(i,j)=0.0_realkind 
             frland(i,j)=frland(i,j)+zfrlaketot
          endif

       enddo
    enddo

    return
  end subroutine lakemasks



  subroutine lakeinit(tlakesfc,     & ! lake surface temperature, including ice and snow (pseudo-prognostic) (k), 1, in
       depthlake,    & ! depth of lake, (m), in
       tlakesnow,    & ! temperature of snow over the lake ice (k), 2, out
       tlakeice,     & ! temperature of lake ice (k), 3, out
       tlakemean,    & ! mean water lake temperature (k), 4, out
       tlakeml,      & ! mixed layer lake temperature (k), 5, out
       tlakebot,     & ! lake bottom temperature (k), 6, out
       ctlake,       & ! chape-factor, 7, out
       sndepthlake,  & !snow depth over lake ice (m), 8, out
       icedepthlake, & ! lake ice depth (m), 9, out
       mldepth,      & ! mixed layer depth (m), 10, out
       t_b1lake,     & ! temperature on the extremum on bottom sedoments (k),
                                ! when bottom sediments = off, use dummy value, 11, out
       h_b1lake)       ! depth of the bottom sediments temperature extremum (m) 
    ! when bottom sediments = off, use dummy value), 12, OUT
    ! Initialise lake prognostic variables
    ! NOW: The lake surface temperature  (input) is equal to
    ! deep soil temperature; all lake is mixed!!! 
    ! there is no ice and snow cower (runs from September, Europe)
    !--------------------------------------------------------------------------------------------



    implicit none

    real(kind=realkind), intent(in)  ::tlakesfc 
    real(kind=realkind), intent(in)  ::depthlake 
    real(kind=realkind), intent(out) ::tlakesnow, tlakeice, tlakemean
    real(kind=realkind), intent(out) ::tlakeml,tlakebot, ctlake, sndepthlake
    real(kind=realkind), intent(out) ::icedepthlake, mldepth, t_b1lake, h_b1lake

    ! locals
    real(kind=ireals) :: depthinit=7.0_ireals


    tlakesnow = tpl_t_f 
    tlakeice = tpl_t_f 
    if(depthlake<=depthinit)then 
       tlakeml = tlakesfc
       ctlake = c_t_min
       mldepth = depthlake 
       tlakebot = tlakeml
       tlakemean = tlakeml
    else 
       tlakeml = tlakesfc
       ctlake = c_t_min
       mldepth = depthinit
       tlakebot = tpl_t_r
       tlakemean = tlakeml - ctlake*(1.0_realkind - mldepth/depthlake)*(tlakeml-tlakebot)
    endif
    t_b1lake = t_bs_lk+1._realkind
    h_b1lake = depth_bs_lk
    sndepthlake = 0.0_realkind 
    icedepthlake = 0.0_realkind 

    return
  end subroutine lakeinit



  subroutine flakeinterface3d(nhor,kstart,kstop, dtime, coslat,  momfulake, momfvlake, senflake, latflake, & 
       sswrlake, slwrlake, &
       tlakesfc, tlakesnow, tlakeice, tlakemean,  tlakeml, tlakebot, ctlake, sndepthlake, & 
       icedepthlake, mldepth,  t_b1lake, h_b1lake,  frlake, depthlake, tlakesfctend, tlakesnowtend, tlakeicetend, tlakemeantend, & 
       tlakemltend, tlakebottend, ctlaketend, sndepthlaketend,  icedepthlaketend, mldepthtend, t_b1laketend, h_b1laketend)
!       nstep,along ) 

    !
    !      "Common" atmospheric variables -----------------------------
    !    nhor:    Number of points in  horizontal loop.
    !    kstart:  Start index of horizontal loop.
    !    kstop:   Stop index of horizontal loop.
    !    dtime:   Time step (s),(two time the time step in the dynamics)
    !      coslat:  cosine if latitude
    !      Fluxes from the atmosphere ---------------------------------
    !      momfuLake: momentum flux towards the east (kg/ (ms2)) over lake (with/without ice)
    !      momfvLake: momentum flux towards the north (kg/ (ms2)) over lake (with/without ice)
    !      senfLake: sensible heat flux over lake (with/without ice)
    !      latfLake : latent heat flux over lake (with/without ice)
    !      sswrLake : surface short wave radiation over lake (with/without ice)
    !      slwrLake : surface long wave radiation over lake (with/without ice)
    !     dsnowdt:   Snow intencity
    !      Surface (lake) specific variables ---------------------------
    !      tLakeSfc  :   lake surface temperature, 
    !                    including ice and snow (pseudo-prognostic) (K)
    !      tLakeSnow :  temperature of snow over the lake ice (K)
    !      tLakeIce  :  temperature of lake ice (K)
    !      tLakeMean :  mean water lake temperature (K)
    !      tLakeML   : mixed layer lake temperature (K)
    !      tLakeBot : lake bottom temperature (K)
    !      CTLake   : chape-factor 
    !      SnDepthLake : snow depth over lake ice (m)
    !      IceDepthLake : lake ice depth (m)
    !      MLDepth : mixed layer depth (m)
    !      T_B1Lake : temperature on the extremum on bottom sedoments (K)
    !               : when bottom sediments = off, use dummy value
    !      H_B1Lake: depth of the bottom sediments temperature extremum (m)
    !              : when bottom sediments = off, use dummy value
    !      frLake : farction of lakes of each type
    !      depthLake : depth of lake of each type, (m)

    !    Output

    !   
    !      tLakeSfcTend  : lake surface temperature, 
    !                     : including ice and snow (pseudo-prognostic) tendency
    !      tLakeSnowTend :  temperature of snow over the lake ice tendency
    !      tLakeIceTend  :  temperature of lake ice tendency
    !      tLakeMeanTend :  mean water lake temperature tendency
    !      tLakeMLTend   : mixed layer lake temperature tendency
    !      tLakeBotTend : lake bottom temperature tendency
    !      CTLakeTend   : chape-factor tendency
    !      SnDepthLakeTend : snow depth over lake ice tendency
    !      IceDepthLakeTend : lake ice depth tendency
    !      MLDepthTend : mixed layer depth tendency
    !      T_B1LakeTend : temperature on the extremum of bottom sediments tendency
    !      H_B1LakeTend : depth of the tempersature extremum of bottom sediments tendency




    implicit none 

    ! exchange variables: -----------------------------------------------------------
    integer, intent (in) :: nhor,kstart,kstop
    real(kind=realkind), intent (in) :: dtime
    real(kind=realkind), dimension(nhor),intent (in) :: coslat
!    real(kind=realkind), dimension(nhor), intent (in) :: dsnowdt
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: senflake, latflake,    sswrlake, slwrlake, momfulake, momfvlake 
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: tlakesfc, tlakesnow,  tlakeice, tlakemean, & 
         tlakeml, tlakebot, ctlake, sndepthlake, icedepthlake, mldepth, t_b1lake, h_b1lake,   frlake, depthlake
    real(kind=realkind), dimension(nhor, nlaketype), intent (out) ::  tlakesfctend, tlakesnowtend,  tlakeicetend, tlakemeantend, & 
         tlakemltend, tlakebottend,  ctlaketend, sndepthlaketend, icedepthlaketend, mldepthtend, t_b1laketend, h_b1laketend

    ! local variables 

    integer ihor, ilaketype ! loop indicies
    real (kind = ireals) :: zparcoriolis ! coriolis parameter
    real (kind = ireals) :: zqsum ! flux into the surface (water/ice/snow)
    real (kind = ireals) ::  zqw  ! heat flux into the water [w m^{-2}]
    real (kind = ireals) ::  zqice ! heat flux into the ice [w m^{-2}]
    real (kind = ireals) ::  zqsnow ! heat flux into the snow [w m^{-2}]
    real (kind = ireals) :: zt_sfc_p, zt_sfc_n ! surface temperatues at current and next step
    real (kind = ireals) :: zcheckout        !  checking variable - to check something (optional)
    ! additional - to make the right types
    real (kind = ireals) ::  zdepthlake, zdtime

    real(kind=realkind) zalat

    do ilaketype=1,nlaketype 
       do ihor=kstart,kstop

          ! not yet use switches here (not to debug this part of code just now) 
          if(frlake(ihor,ilaketype)<fraclakemin) then ! fraction of lake is too small - do nothing

             tlakesfctend(ihor, ilaketype) = 0.0_realkind 
             tlakesnowtend(ihor, ilaketype) = 0.0_realkind 
             tlakeicetend(ihor, ilaketype) = 0.0_realkind 
             tlakemeantend(ihor, ilaketype) = 0.0_realkind  
             tlakemltend(ihor, ilaketype) = 0.0_realkind 
             tlakebottend(ihor, ilaketype) = 0.0_realkind 
             ctlaketend(ihor, ilaketype) = 0.0_realkind 
             sndepthlaketend(ihor, ilaketype) = 0.0_realkind 
             icedepthlaketend(ihor, ilaketype) = 0.0_realkind 
             mldepthtend(ihor, ilaketype) = 0.0_realkind 
             t_b1laketend(ihor, ilaketype) = 0.0_realkind 
             h_b1laketend(ihor, ilaketype) = 0.0_realkind   

          else ! fraction of lake is big enough - do everything

             !  calculate coriolis parameter
             zparcoriolis =  2.0_ireals*omega_earth*sin(acos(coslat(ihor)))

             ! set "background _flk" values

             !  take care of fluxes
             i_atm_flk = sswrlake(ihor, ilaketype) 
             u_star_w_flk = sqrt(sqrt(momfulake(ihor, ilaketype) **2.0_realkind + &
                  momfvlake(ihor, ilaketype) **2.0_realkind)/tpl_rho_w_r)   
             !   dmsnowdt_flk = dsnowdt(ihor) ! kk, 03.02.05 - switch off snow!
             dmsnowdt_flk = 0.0_ireals
             !  arrangements necessary for heat flux: 
             zqsum = slwrlake(ihor, ilaketype)  & ! (notice the signs !)
                  + senflake(ihor, ilaketype) + latflake(ihor, ilaketype)
             ! set fluxes depending on kind of surface
             if(icedepthlake(ihor, ilaketype)>=h_ice_min_flk) then  ! ice exists
                if(sndepthlake(ihor, ilaketype)>=h_snow_min_flk) then   ! there is snow above the ice
                   zqsnow = zqsum
                   zqice  = 0.0_ireals
                   zqw    = 0.0_ireals
                else ! no snow above the ice
                   zqsnow = 0.0_ireals
                   zqice  = zqsum
                   zqw    = 0.0_ireals
                end if ! snow condition
             else ! no ice-snow cover
                zqsnow = 0.0_ireals
                zqice  = 0.0_ireals
                zqw    = zqsum
             end if ! ice condition
             q_w_flk = zqw
             q_ice_flk = zqice
             q_snow_flk = zqsnow
             ! prognostic variables 
             t_snow_p_flk = tlakesnow(ihor, ilaketype)     
             t_ice_p_flk  = tlakeice(ihor, ilaketype)         
             t_mnw_p_flk  = tlakemean(ihor, ilaketype)       
             t_wml_p_flk  = tlakeml(ihor, ilaketype)       
             t_bot_p_flk  = tlakebot(ihor, ilaketype)
             t_b1_p_flk   = t_b1lake(ihor, ilaketype)        
             c_t_p_flk    = ctlake(ihor, ilaketype)        
             h_snow_p_flk = sndepthlake(ihor, ilaketype)      
             h_ice_p_flk  = icedepthlake(ihor, ilaketype)       
             h_ml_p_flk   = mldepth(ihor, ilaketype)        
             h_b1_p_flk   = h_b1lake (ihor, ilaketype)

             ! default values for next time step 
             t_snow_n_flk = t_snow_p_flk
             t_ice_n_flk  = t_ice_p_flk         
             t_mnw_n_flk  = t_mnw_p_flk      
             t_wml_n_flk  = t_wml_p_flk      
             t_bot_n_flk  = t_bot_p_flk
             t_b1_n_flk   = t_b1_p_flk        
             c_t_n_flk    = c_t_p_flk       
             h_snow_n_flk = h_snow_p_flk      
             h_ice_n_flk  = h_ice_p_flk       
             h_ml_n_flk   = h_ml_p_flk        
             h_b1_n_flk   = h_b1_p_flk   

             ! set value to be sent 
             zt_sfc_p = tlakesfc(ihor, ilaketype)

             ! make right types (kind=ireals)
             zdepthlake = depthlake(ihor, ilaketype)
             zdtime = dtime

             !  compute short-wave radiation fluxes (positive downward)

             call flake_radflux_modif ( zdepthlake,  &
                  opticpar_water, opticpar_ice, opticpar_snow )
             !  advance flake variables

             zalat=acos(coslat(ihor))*180.0_realkind /3.141592654_realkind 

             call flake_driver(zdepthlake,depth_bs_lk,t_bs_lk,zparcoriolis,&
                  opticpar_water%extincoef_optic(1),  &
                  zdtime, zt_sfc_n,  zcheckout)

             ! make tendencies of advanced variables 
             tlakesfctend(ihor, ilaketype) = (zt_sfc_n - zt_sfc_p) /dtime
             tlakesnowtend(ihor, ilaketype) = (t_snow_n_flk - t_snow_p_flk) /dtime
             tlakeicetend(ihor, ilaketype) = (t_ice_n_flk - t_ice_p_flk) / dtime
             tlakemeantend(ihor, ilaketype) = (t_mnw_n_flk - t_mnw_p_flk ) / dtime 
             tlakemltend(ihor, ilaketype) = (t_wml_n_flk - t_wml_p_flk) / dtime
             tlakebottend(ihor, ilaketype) = (t_bot_n_flk - t_bot_p_flk) / dtime
             ctlaketend(ihor, ilaketype) = (c_t_n_flk - c_t_p_flk) / dtime
             sndepthlaketend(ihor, ilaketype) = (h_snow_n_flk - h_snow_p_flk) / dtime
             icedepthlaketend(ihor, ilaketype) = (h_ice_n_flk - h_ice_p_flk) / dtime
             mldepthtend(ihor, ilaketype) = ( h_ml_n_flk -  h_ml_p_flk) / dtime
             t_b1laketend(ihor, ilaketype) = (t_b1_n_flk - t_b1_p_flk) / dtime
             h_b1laketend(ihor, ilaketype) = (h_b1_n_flk - h_b1_p_flk ) / dtime

          end if ! fraction of lake condition     

       end do ! loop over hor. grid points
    end do ! loop over lake types

    return 
  end subroutine flakeinterface3d


  subroutine flake_radflux_modif ( depth_w,  & 
       opticpar_water, opticpar_ice, opticpar_snow )       

    !------------------------------------------------------------------------------
    !
    ! description:
    !
    !  computes the radiation fluxes 
    !  at the snow-ice, ice-water, air-water, 
    !  mixed layer-thermocline and water column-bottom sediment interfaces,
    !  the integral-mean radiation flux over the mixed layer,
    !  and the integral-mean radiation flux over the thermocline.
    !
    !
    ! current code owner: dwd, dmitrii mironov
    !  phone:  +49-69-8062 2705
    !  fax:    +49-69-8062 3721
    !  e-mail: dmitrii.mironov@dwd.de
    !
    ! history:
    ! version    date       name
    ! ---------- ---------- ----
    ! n.nn       yyyy/mm/dd dmitrii mironov
    !  initial release
    ! !version 1!  !2004.05.11     katherina kourzeneva
    ! make input radiation flux to be with regard of albedo already -
    ! needed to be included into rca model
    ! code description:
    ! language: fortran 90.
    ! software standards: "european standards for writing and
    ! documenting exchangeable fortran 90 code".

    !
    ! declarations:
    !
    ! modules used:




    implicit none


    !
    ! declarations

    !  input (procedure arguments)

    real (kind = ireals), intent(in) ::   &
         depth_w                            ! the lake depth [m]

    type (opticpar_medium), intent(in) :: & 
         opticpar_water                    , & ! optical characteristics of water
         opticpar_ice                      , & ! optical characteristics of ice
         opticpar_snow                         ! optical characteristics of snow 

    integer (kind = iintegers) :: & ! help variable(s)
         i                             ! do loop index


    !  start calculations
    !------------------------------------------------------------------------------

    if(h_ice_p_flk>=h_ice_min_flk) then            ! ice exists
       if(h_snow_p_flk>=h_snow_min_flk) then        ! there is snow above the ice
          i_snow_flk = i_atm_flk
          i_bot_flk = 0.0_ireals
          do i=1_iintegers, opticpar_snow%nband_optic
             i_bot_flk = i_bot_flk +                    & 
                  opticpar_snow%frac_optic(i)*exp(max(-89._realkind,-opticpar_snow%extincoef_optic(i)*h_snow_p_flk))
          end do
          i_ice_flk  = i_snow_flk*i_bot_flk
       else                                           ! no snow above the ice 
          i_snow_flk = i_atm_flk  
          i_ice_flk  = i_atm_flk
       end if
       i_bot_flk = 0.0_ireals
       do i=1_iintegers, opticpar_ice%nband_optic
          i_bot_flk = i_bot_flk +                      & 
               opticpar_ice%frac_optic(i)*exp(max(-89._realkind,-opticpar_ice%extincoef_optic(i)*h_ice_p_flk )) 
       end do
       i_w_flk      = i_ice_flk*i_bot_flk
    else                                             ! no ice-snow cover
       i_snow_flk   = i_atm_flk  
       i_ice_flk    = i_atm_flk
       i_w_flk      = i_atm_flk
    end if

    if(h_ml_p_flk>=h_ml_min_flk) then           ! radiation flux at the bottom of the mixed layer
       i_bot_flk = 0.0_ireals
       do i=1_iintegers, opticpar_water%nband_optic
          i_bot_flk = i_bot_flk +            & 
               opticpar_water%frac_optic(i)*exp(max(-89._realkind,-opticpar_water%extincoef_optic(i)*h_ml_p_flk)) 
       end do
       i_h_flk = i_w_flk*i_bot_flk
    else                                          ! mixed-layer depth is less then a minimum value
       i_h_flk = i_w_flk
    end if

    i_bot_flk = 0.0_ireals                         ! radiation flux at the lake bottom
    do i=1_iintegers, opticpar_water%nband_optic
       i_bot_flk = i_bot_flk +              & 
            opticpar_water%frac_optic(i)*exp(max(-89._realkind,-opticpar_water%extincoef_optic(i)*depth_w)) 
    end do
    i_bot_flk = i_w_flk*i_bot_flk

    if(h_ml_p_flk>=h_ml_min_flk) then           ! integral-mean radiation flux over the mixed layer
       i_intm_0_h_flk = 0.0_ireals
       do i=1_iintegers, opticpar_water%nband_optic
          i_intm_0_h_flk = i_intm_0_h_flk +                                &
               opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
               (1.0_ireals - exp(max(-89._realkind,-opticpar_water%extincoef_optic(i)*h_ml_p_flk)))
       end do
       i_intm_0_h_flk = i_w_flk*i_intm_0_h_flk/h_ml_p_flk
    else
       i_intm_0_h_flk = i_h_flk
    end if

    if(h_ml_p_flk<=depth_w-h_ml_min_flk) then   ! integral-mean radiation flux over the thermocline
       i_intm_h_d_flk = 0.0_ireals 
       do i=1_iintegers, opticpar_water%nband_optic
          i_intm_h_d_flk = i_intm_h_d_flk +                                &
               opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
               ( exp(max(-89._realkind,-opticpar_water%extincoef_optic(i)*h_ml_p_flk))             &
               - exp(max(-89._realkind,-opticpar_water%extincoef_optic(i)*depth_w) ))
       end do
       i_intm_h_d_flk = i_w_flk*i_intm_h_d_flk/(depth_w-h_ml_p_flk)
    else
       i_intm_h_d_flk = i_h_flk
    end if

    ! kk - security
    i_bot_flk = max(i_bot_flk, c_small_flk )

    !------------------------------------------------------------------------------
    !  end calculations


  end subroutine flake_radflux_modif

  subroutine turbfluxlakeicesnow(nhor,nlev,kstart,kstop, &
       dtime, &
       t,q,u,v,gpot,dph, &
       ps, &
       tlakesfc, &
       sndepthlake, &
       icedepthlake, &
       frlake, &
       t2mlake,  q2mlake, &                                 ! input-output
       u10mlake, v10mlake, &
       momfulake, momfvlake , &
       senflake, latflake, evaplake, &     
       ustarlake, &
       coslat)


    !cl    input:

    !
    !      "common" atmospheric variables -----------------------------
    !cl    nhor:    number of points in  horizontal loop.
    !cl    nlev:    number of vertical levels.
    !cl    kstart:  start index of horizontal loop.
    !cl    kstop:   stop index of horizontal loop.
    !cl    dtime:   time step (s),(two time the time step in the dynamics)
    !cl    t    :   temperature (k),(3d field)
    !cl    q    :   specific humidity ,(kg/kg),(3d field)
    !cl    u    :   wind component towards the east ,(m/s), (3d field)
    !cl    v    :   wind component towards the north,(m/s), (3d field)
    !cl    gpot :   geopotential height (m2/s2)
    !cl    dph  :   pressure difference (pa) between half levels.
    !cl    ps   :   surface pressure (pa)
    !      surface (lake) specific variables ------------------------------------------------------
    !      tlakesfc  :   lake surface temperature, including ice and snow 
    !                    (pseudo-prognostic) (k)
    !      sndepthlake : lake snow depth (m)
    !      icedepthlake : lake ice depth (m)
    !      frlake : farction of lakes of each type

    !      input-output

    !cl    z0lake:  lake surface roughness (m) 
    !      t2mlake: 2 meter temperature (k) over lake (with/without ice)
    !      q2mlake: 2 mer specific humidity (kg/kg) over lake (with/without ice)
    !cl    u10lake: 10 meter u-wind over lake (with/without ice)
    !cl    v10lake: 10 meter v-wind over lake (with/without ice)
    !      momfulake: momentum flux towards the east (kg/ (ms2)) over lake (with/without ice)
    !      momfvlake: momentum flux towards the north (kg/ (ms2)) over lake (with/without ice)
    !      senflake: sensible heat flux over lake (with/without ice)
    !      latflake : latent heat flux over lake (with/without ice)
    !      evaplake : evaporation over lake (with/without ice)
    !      ustarlake : surface friction velocity (m/s) for lake (with/without ice)

    ! uses : ---------------------------------------------------------------------------------


    use confys
    use ctun
    use escom
    use comdfb
    implicit none

    ! exchange variables: --------------------------------------------------------------------
    integer, intent (in) :: nhor,nlev,kstart,kstop
    !c
    real(kind=realkind), intent (in) :: dtime
    real(kind=realkind), dimension(nhor,nlev), intent (in) ::  gpot,q,t,u,v
    real(kind=realkind), dimension(nhor,nlev+1), intent (in) :: dph
    real(kind=realkind), dimension(nhor), intent (in) :: ps
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: tlakesfc, &
         sndepthlake, icedepthlake, frlake

    real(kind=realkind), dimension(nhor, nlaketype), intent (inout) :: t2mlake, q2mlake, &
         u10mlake, v10mlake, momfulake, momfvlake, senflake, latflake, &
         evaplake, ustarlake


    real(kind=realkind), dimension(nhor),intent (in) :: coslat
    real(kind=realkind) zalat

    ! local variables ------------------------------------------------------------------------

    integer ihor, ilaketype ! loop indicies

    !  - this mark means code copied from slf_surf_sea_ice.f90

    real(kind=realkind) zlakeiceon, zlakeiceoff, zlakesnowon, zlakesnowoff, &
         zdens,zvirnl,zdup2,zvel,zcrdq, &
         z03,zriq,zrous,zroumin,zslask,zhnlev,z01,z02,zcneut,zstaon, &
         ztsi,zepcr,zstaoff,zria,zdr,zcdrag,zcdrgh,zm,zsecu, &
         zfmx,zqb,zcams3,zcharg,zqam,ztam,zlat,zh,zuneg,zupos, &
         zunlev,zvnlev,zrougl,zustar,zmoin,zsl1,ztotf,zl2,zln2, &
         zl2lim,zust,zvst,zthst,zqstar,zunson,zunsoff,zt2m,zq2m, &
         zln2k,zy,zy2,zrkar,zfrlim, &
         ztdum, &
         zqseff,ztseff,z10m, &
         zqs,zra,zpar,zqlim,zueps,zqc,zqd,zqdh,zcons1, &
         z2m,zsl11,zhm,zq,zrepac,zcrit,zri,zvneg,zvpos,ze, &
         zx10,zu10,zv10,zl10,zln10,zln10k,zrpi4, &
         zicelayer, &
         zzustar,zevap
    
    
    real(kind=realkind) zrdt,zrocg,zd1,zslask1,ztlamda,zct, &
         ztmelt,zsnlayer, &
         rhoice, &
         zrsfl, &
         zwsat,zeps,zsncrit, &
         zk1,zk2,zk,zfrsnasymp,zsfdist, &
         zsnlim,&
         zsnswcrit, &
         ztdumout,zqdumout,zudumout
    !c       additional local variables for calculation of
    !c       heat and momentum flux over rough and smooth sea.
    real(kind=realkind) zcvis,zdcharg, zcm,zch,z1r3 
    real(kind=realkind) znrz2,zhrmax


    !c     1.2: define local constants
    !c     molecular kinematic viscity of air = 'zcvis'
    zcvis=1.5e-5_realkind 
    z1r3=1.0_realkind /3.0_realkind 

    !c     zcm=1./(conm*delta_h*pr**(2/3), conm=2.*zqb/(3.*zqb*zqc).
    !c     zch=1./(conh*delta_h*pr**(-2/3)), conh=3.*zqb/(3*zqb*zqc).
    !c     delta_h=0.17 (deardorff et al., 1969, townsend, 1964).
    !c     pr=0.71 is the prandtl number.

    zfrlim=0.01_realkind 
    zcm=0.9855_realkind 
    zch=0.9363_realkind 
    zhrmax=0.50_realkind 


    zqlim=0.0001_realkind 
    zueps=0.1_realkind 
    zroumin=5.e-5_realkind 
    zqb=5.0_realkind 
    zqc=5.0_realkind 
    zqd=5.0_realkind 
    zqdh=1.0_realkind 
    zrpi4=pi/4.0_realkind 
    zcrdq=1.0_realkind /epsilo-1.0_realkind 
    zcams3=3.0_realkind *zqb*zqc*carman**2.0_realkind
    zcons1=ccpq*cpair
    zhnlev=-rair*288.15_realkind *clog/gravit
    zsl1=2.0_realkind *rair*gravit*carman/cpair
    zsl11=0.61_realkind *cpair
    z2m=2.0_realkind 
    zk1=0.2_realkind 
    zk2=1.e-06_realkind 
    zk=exp(-zk2*dtime)     
    !car010502
    znrz2=zhnlev/z2m
    z10m=10.0_realkind 
    zrkar=1.0_realkind /carman
    zm=4.0_realkind *zrkar
    zhm=zm
    zq=zm
    zrepac=1.0_realkind /(epsilo*acrit)
    ztdum=280.0_realkind 
    ztdumout=99.0_realkind 
    zqdumout=-1.e-2_realkind 
    zudumout=99.0_realkind 
    !car010502nwn
    zcharg= 0.014_realkind /gravit
    zdcharg=0.018_realkind /gravit
    !carnwn
    !c     zcharg=0.032/gravit
    zcrit=acrit
    zepcr=acrit*epsilo
    zsecu=1.e-7_realkind 
    zl2lim=-4.0_realkind 
    zicelayer=0.1_realkind      !  minimum ice thickness where heat cond is neglected
    zwsat=0.1_realkind 
    zeps=1.e-06_realkind 
    zsncrit=0.03_realkind 
    zsnswcrit=0.0001_realkind 
    zsfdist=0.6_realkind 
    zfrsnasymp=0.95_realkind 
    zsnlim=0.0015_realkind 
    zsnlayer=0.15_realkind 
    !cps
    !cps
    !csg031203
    zrdt=1.0_realkind /dtime
    zd1=0.07_realkind 
    rhoice=920.0_realkind       ! ice density
    zrocg=2.05e06_realkind     ! ice volumetric heat capacity (j m-3 k-1)
    ztlamda=2.22_realkind      ! ice thermal conductivity
    ztmelt=273.15_realkind 
    !c
    zrsfl=dtime/rhoh2o
    zct=1.0_realkind /(zd1*zrocg)
    zslask1=0.5_realkind *dtime*zct


    do ilaketype=1,nlaketype ! ------------------------------------
       do ihor=kstart,kstop 

          ! kkps, debug, 070605:
          zalat=acos(coslat(ihor))*180.0_realkind /3.141592654_realkind 
          zqam=max(q(ihor,nlev),zqlim)

          ! switchers ---------------------------------------------------
          if(frlake(ihor,ilaketype)<fraclakemin) then ! fraction of lake is too small 
             zlakeiceon = 0.0_realkind 
             zlakeiceoff = 1.0_realkind 
             zlakesnowon = 0.0_realkind 
             zlakesnowoff = 1.0_realkind 
          else ! fraction of lake is big enough
             if(icedepthlake(ihor,ilaketype)<h_ice_min_flk) then ! there is no ice
                zlakeiceon = 0.0_realkind 
                zlakeiceoff = 1.0_realkind 
                zlakesnowon = 0.0_realkind 
                zlakesnowoff = 1.0_realkind   
             else ! there is ice
                if(sndepthlake(ihor,ilaketype)<h_snow_min_flk) then ! there is no snow on ice
                   zlakeiceon = 1.0_realkind 
                   zlakeiceoff = 0.0_realkind 
                   zlakesnowon = 0.0_realkind 
                   zlakesnowoff = 1.0_realkind 
                else ! there is snow on ice
                   zlakeiceon = 1.0_realkind 
                   zlakeiceoff = 0.0_realkind 
                   zlakesnowon = 1.0_realkind 
                   zlakesnowoff = 0.0_realkind 
                end if ! snow condition
             end if ! ice condition
          end if ! fraction of lake condition
          ! ---------------------------------------------------------

          ztsi = zlakeiceon*min(tlakesfc(ihor,ilaketype),ztmelt) +ztdum*zlakeiceoff

          zdens=dph(ihor,nlev+1)/gpot(ihor,nlev)
          zqseff=zepcr*esati(ztsi)/ &
               (ps(ihor) -(1.0_realkind -epsilo)*zcrit*esati(ztsi))
          zvirnl=1.0_realkind +zcrdq*zqam
          zdup2=max( u(ihor,nlev)**2.0_realkind + v(ihor,nlev)**2.0_realkind,zueps )
          zvel=sqrt(zdup2)
          !c
          !c     ln(hybf(nlev))= ln(pf(nlev))-ln(psm) => z03=-fi/(t*v**2)
          !c     -------------------------------------------------------------
          z03=rair*clog/zdup2
          zriq =-z03*( gpot(ihor,nlev)/(cpair+zcons1*zqam) + &
               zvirnl*t(ihor,nlev) )      

          zrous =zlakeiceon*zlakesnowon*z0snowice + zlakeiceon*zlakesnowoff*z0ice &
               + zlakeiceoff*z0ice ! - dummy - here is needed, if we follow this philosophy :-)
          !c     ----------------------------------------------------------------
          zslask=(zrous+zhnlev)/zrous
          z01=1.0_realkind /(log(zslask))**2.0_realkind
          z02=sqrt(zslask)*z01
          zcneut=z01*(carman**2.0_realkind)
          !c
          !c     -------------------------------------------------------
          !cl    calculation of richardson number ( zri ) for surface
          !cl    temperature ztsi.
          !c     -------------------------------------------------------
          !cps
          ztseff=ztsi
          !c
          zri=zriq + z03 *( 1.0_realkind +zcrdq*zqseff)*ztseff
          !c
          !cps
          !c     zstaon=0 if zri<=0.
          !c
          if( zri>0.0_realkind  ) then
             zstaon=1.0_realkind 
             zstaoff=0.0_realkind 
          else
             zstaon=0.0_realkind 
             zstaoff=1.0_realkind 
          endif
          !c
          zria=abs(zri)
          !c
          !c     stable richardson number > 0.
          !c     (formulae according to louis et al,1982)
          !c
          zdr=2.0_realkind *zqb*zria/sqrt(1.0_realkind  +zqd*zria)
          zcdrag=zcneut/( 1.0_realkind  +zdr )
          zdr=2.0_realkind *zqb*zria*sqrt(1.0_realkind  +zqdh*zria)
          zcdrgh=zcneut/( 1.0_realkind +zdr )
          !c
          !c     unstable richardson number
          !c     (louis et al,1982: ecmwf)
          !c
          zfmx=zqb*zria/(1.0_realkind +z02*zcams3*sqrt(zria))
          zslask=zcneut*(1.0_realkind +3.0_realkind *zfmx)
          zcdrag=zstaoff*zcneut*(1.0_realkind +2.0_realkind *zfmx) +zstaon*zcdrag
          zcdrgh=zstaoff*zslask +zstaon*zcdrgh
          !c
          !c     aerodynamic resistance for heat='zra'
          !c     ( computation not used if zicensoff=1.)
          !c     -------------------------------------
          !c

          zra=zlakeiceon/(zcdrgh*zvel) +zlakeiceoff

          !c
          !c     specific humidity at the lowest model layer.
          !c     --------------------------------------------
          ztam=t(ihor,nlev)+gpot(ihor,nlev)/cpair
          !c

          if(tlakesfc(ihor,ilaketype)>=ztmelt) then
             zlat = latvap
          else
             zlat = latvap+latice
          end if

          !c
          !c     sensible and latent heat fluxes over ice.
          !c     --------------------------------------------------------
          zh=-zdens/zra*cpair*(ztsi-ztam)
          ze=-zdens/zra*(zqseff-zqam)*zlat
          !c evaporation from ice [kg/m2 or mm per time step]
          zevap=-ze/zlat*dtime
          zqs=zqseff  

          !c     calculation of t2m,q2m over ice (and snow on ice - kk) 
          !c     wind components must not be zero.
          !c     --------------------------------------------------------
          if( u(ihor,nlev)<0.0_realkind  ) then
             zuneg=1.0_realkind 
             zupos=0.0_realkind 
          else
             zuneg=0.0_realkind 
             zupos=1.0_realkind 
          endif
          !c
          zunlev=zuneg*min(u(ihor,nlev),-zueps)+ &
               zupos*max(u(ihor,nlev),zueps)
          !c
          if( v(ihor,nlev)<0.0_realkind  ) then
             zvneg=1.0_realkind 
             zvpos=0.0_realkind 
          else
             zvneg=0.0_realkind 
             zvpos=1.0_realkind 
          endif
          !c
          zvnlev=zvneg*min(v(ihor,nlev),-zueps)+ &
               zvpos*max(v(ihor,nlev),zueps)
          !c
          !c     avoid values of rougness > 2m
          !c
          zrougl=min( zrous,z2m )
          !c
          zustar=sqrt(zcdrag*zdup2)
          zustar=max(zustar,0.01_realkind )

          ! update friction velocity.
          ! -----------------------------------------------------------
          zzustar=zustar*zustar
          ustarlake(ihor,ilaketype)=ustarlake(ihor,ilaketype)*ustarlake(ihor,ilaketype)
          ustarlake(ihor,ilaketype)=zlakeiceon*zzustar + &
               zlakeiceoff*ustarlake(ihor,ilaketype)
          ustarlake(ihor,ilaketype)=sqrt(ustarlake(ihor,ilaketype))

          ztotf=zh +zsl11*t(ihor,nlev)*ze/zlat
          zmoin=zsl1*ztotf/( 2.0_realkind *(ps(ihor)-dph(ihor,nlev+1))*zustar**3.0_realkind )
          !c
          zl2=z2m*zmoin
          zl2=max( zl2, zl2lim )
          !c
          zln2=log(z2m/zrougl)
          zln2k=zln2*zrkar
          !c
          !c     compute components 'u*' and 'v*' (ustar and vstar)
          !c
          zust=zustar*zunlev/zvel
          zvst=zustar*zvnlev/zvel
          !c
          !c
          !c     compute 't*','q*' (zthst,zqstar)
          !c     -----------------------------------------------------------
          zthst=zh/(cpair*zdens*zustar)
          zqstar=ze/(zlat*zdens*zustar)
          !c
          !c     compute 't*','q*' (zthst,zqstar)
          !c     -----------------------------------------------------------
          zthst=zh/(cpair*zdens*zustar)
          zqstar=ze/(zlat*zdens*zustar)
          !c
          !c     stable case zunson=0, unstable zunson=1
          !c
          if( ztotf<0.0_realkind  ) then
             zunson=1.0_realkind 
             zunsoff=0.0_realkind 
          else
             zunson=0.0_realkind 
             zunsoff=1.0_realkind 
          endif

          zslask = max(-89._realkind,zlakeiceon* &
               ( zunson + zunsoff*(-zhm*zthst*zl2/ztam) ) &
               + zlakeiceoff)

          !cps
          zt2m=ztseff +zthst*zln2k + &
               (ztam-ztseff)*(1.0_realkind - exp(zslask) )

          zslask = max(-89._realkind,zlakeiceon* &
               ( zunson + zunsoff*(-zq*zqstar*zl2/zqam) ) &
               + zlakeiceoff)

          zq2m=zqs +zqstar*zln2k + &
               (zqam-zqs)*(1.0_realkind - exp(zslask) )
          !c
          !c     zy2 dummy if zunson=0
          !c
          zy=sqrt( 1.0_realkind  - 9.0_realkind *zl2*zunson )
          zy=min( zy, sqrt( 8.0_realkind /zrougl ) -1.0_realkind  )
          zy2=zunson*zy +zunsoff

          zy2 = zlakeiceon*zy2 + zlakeiceoff

          zpar= ( zln2-2.0_realkind *log( 0.5_realkind *(1.0_realkind +zy2)) )*zrkar

          t2mlake(ihor,ilaketype)=zlakeiceon* &
               ( zunson*(ztseff+zthst*zpar) + zunsoff*zt2m ) + &
               zlakeiceoff*t2mlake(ihor,ilaketype)

          ! debug, kk, 070205
          !       if(along(ihor)>62.1 .and. along(ihor)<62.2 )then
          !          if(zalat>56.3 .and. zalat<56.4)then
          !            write(555,*) 'nstep=',nstep, along(ihor),zalat 
          !            write(555,*) t2mlake(ihor,ilaketype)
          !            write(555,*) zlakeiceon, zunson*(ztseff+zthst*zpar), zunsoff*zt2m
          !            write(555,*) zlakeiceoff,t2mlake(ihor,ilaketype) 
          !         endif
          !       endif


          q2mlake(ihor,ilaketype)=zlakeiceon* &
               ( zunson*(zqs+zqstar*zpar) + zunsoff*zq2m ) + &
               zlakeiceoff*q2mlake(ihor,ilaketype)

          !c
          !c---- new by anna to have values of wind over ice
          !c
          zl10=z10m*zmoin
          zln10=log(z10m/zrougl)
          zln10k=zln10*zrkar
          !c
          !c     zu10 etc dummy if zunson=1.0
          !c
          zslask=max(-89._realkind,zunson + zunsoff*(-zm*zust*zl10/zunlev))
          zu10=zust*zln10k + zunlev*(1.0_realkind -exp(zslask))
          zslask=max(-89._realkind,zunson + zunsoff*(-zm*zvst*zl10/zvnlev))
          zv10=zvst*zln10k + zvnlev*(1.0_realkind -exp(zslask))
          !c
          !c     zx10 dummy if zunson=0
          !c
          zx10=zunson*(1.0_realkind -15.0_realkind *zl10) + zunsoff
          zx10=sqrt(sqrt(zx10))
          !c
          zpar=(zln10 - (log( 0.5_realkind *(1.0_realkind +zx10*zx10) )+ 2.0_realkind * &
               ( log( 0.5_realkind *(1.0_realkind +zx10))-atan(zx10)+zrpi4 )) )*zrkar
          !c

          u10mlake(ihor,ilaketype) = zlakeiceon * & 
               ( zunson*zust*zpar + zunsoff*zu10 ) &
               + zlakeiceoff * u10mlake(ihor, ilaketype)
          v10mlake(ihor,ilaketype) = zlakeiceon * &
               ( zunson*zvst*zpar + zunsoff*zv10 ) &
               + zlakeiceoff * v10mlake(ihor, ilaketype)
          momfulake(ihor,ilaketype) = zlakeiceon * &
               zdens*zustar*zustar*u(ihor,nlev)/(zvel+zsecu) &
               + zlakeiceoff * momfulake(ihor, ilaketype)
          momfvlake(ihor,ilaketype) = zlakeiceon * &
               zdens*zustar*zustar*v(ihor,nlev)/(zvel+zsecu) &
               + zlakeiceoff * momfvlake(ihor, ilaketype)
          senflake(ihor,ilaketype) = zlakeiceon * zh &
               + zlakeiceoff * senflake(ihor, ilaketype)
          latflake(ihor,ilaketype) = zlakeiceon * ze &
               + zlakeiceoff * latflake(ihor, ilaketype)  
          evaplake(ihor,ilaketype) = zlakeiceon * zevap &
               + zlakeiceoff * evaplake(ihor, ilaketype)  

       end do ! loop over hor. grid points
    end do ! loop over lake types

    return
  end subroutine turbfluxlakeicesnow


  subroutine turbfluxlakewater(nhor,nlev,kstart,kstop, &
       dtime, &
       t,q,u,v,gpot,dph, &
       ps, &
       tlakesfc, &
       sndepthlake, &
       icedepthlake, &
       frlake, &
       z0lake, &     
       t2mlake, q2mlake, &
       u10mlake, v10mlake, &
       momfulake, momfvlake , &
       senflake, latflake, evaplake, &     
       ustarlake)


    !cl    input:

    !
    !      "common" atmospheric variables -----------------------------
    !cl    nhor:    number of points in  horizontal loop.
    !cl    nlev:    number of vertical levels.
    !cl    kstart:  start index of horizontal loop.
    !cl    kstop:   stop index of horizontal loop.
    !cl    dtime:   time step (s),(two time the time step in the dynamics)
    !cl    t    :   temperature (k),(3d field)
    !cl    q    :   specific humidity ,(kg/kg),(3d field)
    !cl    u    :   wind component towards the east ,(m/s), (3d field)
    !cl    v    :   wind component towards the north,(m/s), (3d field)
    !cl    gpot :   geopotential height (m2/s2)
    !cl    dph  :   pressure difference (pa) between half levels.
    !cl    ps   :   surface pressure (pa)
    !      surface (lake) specific variables ---------
    !      tlakesfc  :   lake surface temperature, including ice and snow
    !                     (pseudo-prognostic) (k)
    !      sndepthlake : lake snow depth (m)
    !      icedepthlake : lake ice depth (m)
    !      frlake : farction of lakes of each type

    !      input-output

    !cl    z0lake:  lake surface roughness (m) 
    !      t2mlake: 2 meter temperature (k) over lake (with/without ice)
    !      q2mlake: 2 mer specific humidity (kg/kg) over lake (with/without ice)
    !cl    u10lake: 10 meter u-wind over lake (with/without ice)
    !cl    v10lake: 10 meter v-wind over lake (with/without ice)
    !      momfulake: momentum flux towards the east (kg/ (ms2)) over lake (with/without ice)
    !      momfvlake: momentum flux towards the north (kg/ (ms2)) over lake (with/without ice)
    !      senflake: sensible heat flux over lake (with/without ice)
    !      latflake : latent heat flux over lake (with/without ice)
    !      evaplake : evaporation over lake (with/without ice)
    !      ustarlake : surface friction velocity (m/s) for lake (with/without ice)

    ! uses : -----------


    use confys
    use ctun
    use escom
    use comdfb
    implicit none 

    ! exchange variables
    integer, intent (in) :: nhor,nlev,kstart,kstop
    real(kind=realkind), intent (in) :: dtime
    real(kind=realkind), dimension(nhor,nlev), intent (in) ::  gpot,q,t,u,v
    real(kind=realkind), dimension(nhor,nlev+1), intent (in) :: dph
    real(kind=realkind), dimension(nhor), intent (in) :: ps
    real(kind=realkind), dimension(nhor, nlaketype), intent (in) :: tlakesfc, &
         sndepthlake, icedepthlake, frlake

    real(kind=realkind), dimension(nhor, nlaketype), intent (inout) :: z0lake 

    real(kind=realkind), dimension(nhor, nlaketype), intent (out) :: t2mlake, q2mlake, &
         u10mlake, v10mlake, momfulake, momfvlake, senflake, latflake, &
         evaplake, ustarlake

    ! local variables --

    integer ihor, ilaketype ! loop indicies




    ! some variables and constants may be not used, some useless ones are deleted - kk

    real(kind=realkind) zlakewateron,zlakewateroff,zlakeiceon,zlakeiceoff, &
         zlakesnowon, zlakesnowoff
    real(kind=realkind) &
         zdens,zvirnl,zdup2,zvel,zcrdq, &
         z03,zriq,zrous,zroumin,zslask,zhnlev,z01,z02,zcneut,zstaon, &
         ztsi,zqsi,zepcr,zstaoff,zria,zdr,zcdrag,zcdrgh,zm,zsecu, &
         zfmx,zqb,zcams3,zcharg,zqam,ztam,zlat,zh,zeeff,zuneg,zupos, &
         zunlev,zvnlev,zrougl,zustar,zmoin,zsl1,ztotf,zl2,zln2, &
         zl2lim,zust,zvst,zthst,zqstar,zunson,zunsoff,zt2m,zq2m, &
         zln2k,zy,zy2,zrkar,zfrlim, &
         ztdum, &
         z10m, &
         zra,zpar,zqlim,zueps,zqc,zqd,zqdh,zcons1, &
         z2m,zsl11,zhm,zq,zrepac,zcrit,zri,zvneg,zvpos, &
         zx10,zu10,zv10,zl10,zln10,zln10k,zrpi4, &
         zicelayer, &
         zevap

    real(kind=realkind) zrocg,zd1,zslask1,ztlamda,zct, &
         ztmelt,zsnlayer, &
         rhoice, &
         zrsfl, &
         zwsat,zeps,zsncrit, &
         zk1,zk2,zfrsnasymp,zsfdist, &
         zsnlim,&
         zsnswcrit, &
         ztdumout,zqdumout,zudumout
    !car010502nwn beg
    !c
    !c
    !c       additional local variables for calculation of
    !c       heat and momentum flux over rough and smooth sea.
    !c      --------------------------------------------------
    real(kind=realkind) zrey,zfrey,zslash,zslasq,zslasm,z01h,z01q, &
         zcneuh,zcneuq,zcdrgq,zfint,zcvis,zustars, &
         zlnzrzh,zlnzrzq,zust2,zchargf,zdcharg, &
         zcm,zch,z1r3,zcmol,zqar
    real(kind=realkind) zln,zyn,zyn2,znrz2,zxn,zhrmax

    !c     1.2: define local constants
    !c     molecular kinematic viscity of air = 'zcvis'
    zcvis=1.5e-5_realkind 
    z1r3=1.0_realkind /3.0_realkind 

    !c     zcm=1./(conm*delta_h*pr**(2/3), conm=2.*zqb/(3.*zqb*zqc).
    !c     zch=1./(conh*delta_h*pr**(-2/3)), conh=3.*zqb/(3*zqb*zqc).
    !c     delta_h=0.17 (deardorff et al., 1969, townsend, 1964).
    !c     pr=0.71 is the prandtl number.

    zfrlim=0.01_realkind 
    zcm=0.9855_realkind 
    zch=0.9363_realkind 
    zhrmax=0.50_realkind 

    zqlim=0.0001_realkind 
    zueps=0.1_realkind 
    zroumin=1.e-6_realkind 
    zqb=5.0_realkind 
    zqc=5.0_realkind 
    zqd=5.0_realkind 
    zqdh=1.0_realkind 
    zrpi4=pi/4.0_realkind 
    zcrdq=1.0_realkind /epsilo-1.0_realkind 
    zcams3=3.0_realkind *zqb*zqc*carman**2.0_realkind
    zcons1=ccpq*cpair
    zhnlev=-rair*288.15_realkind *clog/gravit
    zsl1=2.0_realkind *rair*gravit*carman/cpair
    zsl11=0.61_realkind *cpair
    z2m=2.0_realkind 
    zk1=0.2_realkind 
    zk2=1.e-06_realkind  
    !car010502
    znrz2=zhnlev/z2m
    z10m=10.0_realkind 
    zrkar=1.0_realkind /carman
    zm=4.0_realkind *zrkar
    zhm=zm
    zq=zm
    zrepac=1.0_realkind /(epsilo*acrit)
    ztdum=280.0_realkind 
    ztdumout=99.0_realkind 
    zqdumout=-1.e-2_realkind 
    zudumout=99.0_realkind 
    !car010502nwn
    zcharg= 0.014_realkind /gravit
    zdcharg=0.018_realkind /gravit
    !carnwn
    !c     zcharg=0.032/gravit
    zcrit=acrit
    zepcr=acrit*epsilo
    zsecu=1.e-7_realkind 
    zl2lim=-4.0_realkind 
    zicelayer=0.1_realkind      !  minimum ice thickness where heat cond is neglected
    zwsat=0.1_realkind 
    zeps=1.e-06_realkind 
    zsncrit=0.03_realkind 
    zsnswcrit=0.0001_realkind 
    zsfdist=0.6_realkind 
    zfrsnasymp=0.95_realkind 
    zsnlim=0.0015_realkind 
    zsnlayer=0.15_realkind 
    zd1=0.07_realkind 
    rhoice=920.0_realkind       ! ice density
    zrocg=2.05e06_realkind     ! ice volumetric heat capacity (j m-3 k-1)
    ztlamda=2.22_realkind      ! ice thermal conductivity
    ztmelt=273.1_realkind 
    zrsfl=dtime/rhoh2o
    zct=1.0_realkind /(zd1*zrocg)
    zslask1=0.5_realkind *dtime*zct


    do ilaketype=1,nlaketype ! ------------------------------------
       do ihor=kstart,kstop 

          zqam=max(q(ihor,nlev),zqlim)

          ! switchers --------------------------------------------------- 

          ! ice condition here is necessary to prezerve the value of z0lake 
          ! when lake is covered by ice   
          if(frlake(ihor,ilaketype)<fraclakemin) then ! fraction of lake is too small  
             zlakewateron = 0.0_realkind 
             zlakewateroff = 1.0_realkind 
             zlakeiceon = 0.0_realkind 
             zlakeiceoff = 1.0_realkind 
             zlakesnowon = 0.0_realkind 
             zlakesnowoff = 1.0_realkind   
          else ! fraction of lake is big enough
             if(icedepthlake(ihor,ilaketype)<h_ice_min_flk) then ! no ice on lake
                zlakewateron = 1.0_realkind 
                zlakewateroff = 0.0_realkind 
                zlakeiceon = 0.0_realkind 
                zlakeiceoff =1.0_realkind 
                zlakesnowon = 0.0_realkind 
                zlakesnowoff = 1.0_realkind  
             else ! ice on lake
                if (sndepthlake(ihor,ilaketype)<h_snow_min_flk) then ! there is no snow on ice
                   zlakewateron = 1.0_realkind 
                   zlakewateroff = 0.0_realkind 
                   zlakeiceon = 1.0_realkind 
                   zlakeiceoff =0.0_realkind 
                   zlakesnowon = 0.0_realkind 
                   zlakesnowoff = 1.0_realkind      
                else ! there is snow on ice
                   zlakewateron = 1.0_realkind 
                   zlakewateroff = 0.0_realkind 
                   zlakeiceon = 1.0_realkind 
                   zlakeiceoff = 0.0_realkind 
                   zlakesnowon = 1.0_realkind 
                   zlakesnowoff = 0.0_realkind     
                end if ! snow condition
             end if ! ice condition
          end if ! lake condition


          zdens=dph(ihor,nlev+1)/gpot(ihor,nlev)
          zvirnl=1.0_realkind +zcrdq*zqam
          !c     --------------------------------------------------------
          zdup2=max( u(ihor,nlev)**2.0_realkind + v(ihor,nlev)**2.0_realkind,zueps )
          zvel=sqrt(zdup2)
          z03=rair*clog/(zdup2)
          zriq =-z03*( gpot(ihor,nlev)/(cpair+zcons1*zqam) + &
               zvirnl*t(ihor,nlev) )
          !c     2.05. calculate interpolation function zfint for use in
          !c     the calculation of z0-sea,
          !c     assuming rough sea for v &gt; 5.0 m/s
          !c     and smooth sea for v &lt; 3.0 m/s, and applying
          !c     a square root interpolation between the regimes.
          !c     rough sea: z0-sea=beta*ustars*ustars/gravit.
          !c     smooth sea: z0-sea=0.11*1.5e-5/ustars.
          !c     ------------------------------------------------
          zfint=min((zvel-3.0_realkind )/2.0_realkind ,1.0_realkind )
          zfint=max(zfint,0.0_realkind )
          zfint=sqrt(zfint)
          !c


          !cl    computing new surface roughness for water:
          !c     ---------------------------------------------------------
          zrous=max( z0lake(ihor,ilaketype), zroumin )

          zslask=zhnlev/zrous
          z01=1.0_realkind /(log(zslask))**2.0_realkind
          zcneut=z01*(carman**2.0_realkind)
          !c
          !c     -------------------------------------------------------
          !c     2.06. redefine constant beta in charnock relation if
          !c     1-fraction of sea is greater than zfrlim.
          !c     -------------------------------------------------------
          !c
          zchargf=zcharg

          if(1.0_realkind-frlake(ihor,ilaketype)>zfrlim) zchargf=zcharg+zdcharg

          !c     --------------------------------------------------------
          !cl    2.1 calculation of richardson number ( zri ) for surface
          !cl        temperature  tsea(jl).
          !c

          ztsi=tlakesfc(ihor,ilaketype)
          if(tlakesfc(ihor,ilaketype)<270.0_realkind )ztsi=tmelt
          if(tlakesfc(ihor,ilaketype)>350.0_realkind )ztsi=350.0_realkind 
          ztsi = ztsi*zlakewateron*zlakeiceoff*zlakesnowoff + &
               ztdum*(zlakewateroff+zlakeiceon)
          !c
          zqsi=zepcr*esatw(ztsi)/ &
               (ps(ihor)-(1.0_realkind -epsilo)*zcrit*esatw(ztsi))
          !c

          zri=zriq + z03 *( 1.0_realkind +zcrdq*zqsi)*ztsi
          !c
          !c     zstaon=1 if zri>0.
          !c
          if( zri>0.0_realkind  ) then
             zstaon=1.0_realkind 
             zstaoff=0.0_realkind 
          else
             zstaon=0.0_realkind 
             zstaoff=1.0_realkind 
          endif
          !c
          zria=abs(zri)
          !car010502nwn beg
          !cnwn  calculation of surface drag coefficients
          !c     -----------------------------------------------
          !c     2.15. calculate free convection velocity 'zslask'

          !c     free convection roughness length.
          !c     -----------------------------------------------
          zslask=(zria*zdup2*zcvis/zhnlev)**z1r3
          z02=sqrt(zslask*zhnlev/zcvis)
          !c
          !cnwn  step 1: momentum drag
          !c
          !c     stable richardson number &gt; 0.
          !c     (formulae according to louis et al,1982)
          !c     -----------------------------------------
          zdr=2.0_realkind *zqb*zria/sqrt(1.0_realkind  +zqd*zria)
          zcmol=zcvis/(zhnlev*zvel)
          zcdrag=zcneut/( 1.0_realkind  +zdr )+zcmol
          !c
          !c     unstable richardson number
          !c     (formulae according to louis et al,1982)
          !c     -----------------------------------------
          zfmx=zqb*zria/(1.0_realkind +zcm*z01*z02*zcams3*sqrt(zria))
          zcdrag=zstaoff*zcneut*(1.0_realkind +2.0_realkind *zfmx) +zstaon*zcdrag
          !c
          !c
          !c     step 2: kinematic heat and moisture drag
          !c     ---------------------------------------
          !c
          zust2=zcdrag*zdup2
          zustars=sqrt(zust2)
          !c
          !c
          !c     ---------------------------------------
          !c     2.18. calculate surface reynolds number
          !c     ---------------------------------------
          !c
          zrey=zrous*zustars/zcvis
          !c
          !c     calculate ln(z0m/z0h) and ln(z0m/z0q)
          !c     according to ln(z0m/zoh)=2.48*sqrt(sqrt(re*))-2.
          !c     and          ln(z0m/z0q)=2.28*surt(surt(re*))-2.
          !c     taken from j.r. garratt, 1992: the atmospheric
          !c     boundary layer, p 102.
          !c     ------------------------------------------------
          !c
          zfrey = sqrt(zrey)
          zfrey = sqrt(zfrey)
          zslasm = 1.0_realkind /sqrt(z01)
          !car   zslash = (0.05*zfint+2.43)*zfrey - 2.
          !car   zslasq = zslash + (0.5*zfint-0.70)*zfrey
          !car010427 new coefficient to reduce chn
          zslash = (0.92_realkind *zfint+2.43_realkind )*zfrey - 2.0_realkind 
          zslasq = zslash + (0.08_realkind *zfint-0.70_realkind )*zfrey
          zlnzrzh = zslash+zslasm
          zlnzrzq = zslasq+zslasm
          z01h = 1.0_realkind /(zlnzrzh*zslasm)
          z01q = 1.0_realkind /(zlnzrzq*zslasm)
          !c
          !c     ------------------------------------------
          !c     chn=k*k/(ln(z/z0m)*ln(z/z0h))=
          !c         cmn/(1.+ln(z0m/z0h)/ln(z/z0m))
          !c     ------------------------------------------
          zcneuh = zcneut/(1.0_realkind  + zslash/zslasm)
          zcneuq = zcneut/(1.0_realkind  + zslasq/zslasm)
          !c
          !c     stable richardson number &gt; 0,
          !c     (formulae according to louis et al., 1982)
          !c     ------------------------------------------
          !c
          !car010427new
          zqar=2._realkind
          !c     zdr=2.0*zqb*zria*sqrt(1. +zqdh*zria)
          zdr=2.0_realkind *zqar*zqb*zria*sqrt(1.0_realkind  +zqdh*zria)
          zdr=1.0_realkind /(1.0_realkind  + zdr)
          zcdrgh=zcneuh*zdr+zcmol/0.71_realkind 
          zcdrgq=zcneuq*zdr+zcmol/0.60_realkind 
          !c
          !c     ------------------------------------------
          !c     unstable richardson number
          !c     (formulae according to louis et al., 1982)
          !c     but with a free convection modification.
          !c     ------------------------------------------
          !c
          zdr=zcams3*sqrt(zria)
          !car010427new
          zqar=2.0_realkind 
          !c     zslask=zcneuh*(1.+3.*(zqb*zria/(1.+zch*z01h*z02*zdr)))
          zslask=zcneuh*(1.0_realkind +3.0_realkind *(zqb*zria/ &
               (1.0_realkind +zqar*zch*z01h*z02*zdr)))
          zcdrgh=zstaoff*zslask+zstaon*zcdrgh
          !car010427new
          !c     zslask=zcneuq*(1.+3.*(zqb*zria/(1.+zch*z01q*z02*zdr)))
          zslask=zcneuq*(1.0_realkind +3.0_realkind * &
               (zqb*zria/(1.0_realkind +zqar*zch*z01q*z02*zdr)))
          zcdrgq=zstaoff*zslask+zstaon*zcdrgq
          !c
          !c
          !c     roughness calc for sea:
          !c     -----------------------------------------------------
          !c     calculate z0-sea, using smooth interpolation in
          !c     wind speed between smooth and rough surface.
          !c     calculation of zfint is performed in 2.05 above.
          !c     zchargf is calculated in 2.06.
          !c     -----------------------------------------------------
          !c
          zrous=(1.0_realkind -zfint)*0.11_realkind *zcvis/zustars &
               + zfint*zust2*zchargf 
          ! if ice (or snow) on the lake appears, z0lake is preserved as
          !  z0ice (z0snowice) ,but not previous z0lake  
          z0lake(ihor,ilaketype)= &
               zlakewateron*zlakeiceoff*zlakesnowoff * max(zrous,zroumin) &
               + zlakewateron*zlakeiceon*zlakesnowoff*z0ice &
               + zlakewateron*zlakeiceon*zlakesnowon*z0snowice &
               + zlakewateroff*1.0e-3_realkind 
          !c     ztam is pot. temp. at nlev. it is assumed that
          !c     theta/temp = 1.
          ztam=t(ihor,nlev) +gpot(ihor,nlev)/cpair
          !c     2.2  flux computations
          !cps    zlat = latvap+latice*(0.5-sign(0.5,ztsi-tmelt))
          !c always water!!!
          zlat = latvap
          !c     sensible heat flux 'zh' and latent heat flux 'zeeff'
          !c     are negative upwards
          !c      aerodynamic resistance for heat and moisture='zra'
          !c     ( computation not used if zseaoff=1.)
          zra=zlakewateron*zlakeiceoff/(zcdrgh*zvel) + zlakeiceon + zlakewateroff
          zh=-zdens/zra*cpair*(ztsi-ztam)
          zra=zlakewateron*zlakeiceoff/(zcdrgq*zvel) + zlakeiceon + zlakewateroff
          zeeff = -zdens/zra * (zqsi - zqam) * zlat
          zevap=-zeeff/zlat*dtime
          !c     calculation of t2m,q2m
          !c     wind components must not be zero.
          if( u(ihor,nlev)<0.0_realkind  ) then
             zuneg=1.0_realkind 
             zupos=0.0_realkind 
          else
             zuneg=0.0_realkind 
             zupos=1.0_realkind 
          endif
          !c
          zunlev=zuneg*min(u(ihor,nlev),-zueps)+ &
               zupos*max(u(ihor,nlev),zueps)
          !c
          if( v(ihor,nlev)<0.0_realkind  ) then
             zvneg=1.0_realkind 
             zvpos=0.0_realkind 
          else
             zvneg=0.0_realkind 
             zvpos=1.0_realkind 
          endif
          !c
          zvnlev=zvneg*min(v(ihor,nlev),-zueps)+ &
               zvpos*max(v(ihor,nlev),zueps)
          !c
          !c     avoid values of rougness > 2m
          !c

          zrougl=min( z0lake(ihor,ilaketype),z2m )

          !c     ------------------------------------------------------
          !cc      zustar=sqrt(zcdrag*zdup2)
          !cc      zustar=max(zustar,0.01)
          !car010502 nwn
          zustar=sqrt(zcdrag*zdup2)
          zustar=max(zustar,0.001_realkind )
          !cps051017
          !cps051017      zustar=zustar*1.25
          !cps051017
          !car

          ustarlake(ihor,ilaketype)=zustar*zustar*zlakewateron*zlakeiceoff
          ustarlake(ihor,ilaketype)=sqrt(ustarlake(ihor,ilaketype))

          !c     ---------------------------------------------------------
          ztotf=zh+zsl11*t(ihor,nlev)*zeeff/zlat
          !c     ---------------------------------------------------------
          !c     'zmoin' is the uinverse monin-obukhov length
          !c     ---------------------------------------------------------
          zmoin=zsl1*ztotf/( 2.0_realkind *(ps(ihor)-dph(ihor,nlev+1))*zustar**3.0_realkind )
          !c
          zl2=z2m*zmoin
          !car010502
          zln=zhnlev*zmoin
          zl2=max( zl2, zl2lim )
          zln2=log(z2m/zrougl)
          zln2k=zln2*zrkar
          !c
          !c     compute components 'u*' and 'v*' (ustar and vstar)
          !c
          zust=zustar*zunlev/zvel
          zvst=zustar*zvnlev/zvel
          !c
          !c     temperature-and humidity fluctuations teta*,q*,
          !c     ( 'zthst', 'zqstar' )
          !c
          !car010502 nwn begin
          zthst=zh/(cpair*zdens*zustar*sqrt( zcneuh/zcneut ))
          zqstar=zeeff/(zlat*zdens*zustar*sqrt( zcneuq/zcneut ))
          !c
          !c     unstable case zunson=1.0
          !c
          !c
          if( ztotf<0.0_realkind  ) then
             zunson=1.0_realkind 
             zunsoff=0.0_realkind 
          else
             zunson=0.0_realkind 
             zunsoff=1.0_realkind 
          endif
          !c
          !c     zt2m etc dummy if zunson=1.0
          !c
          zslask=max(-89._realkind,zlakewateron*zlakeiceoff*( zunson + zunsoff* &
               (-zhm*zthst*zl2/ztam) ) + zlakeiceon + zlakewateroff)
          zt2m=ztsi+zthst*zln2k + &
               (ztam-ztsi)*(1.0_realkind - exp(zslask) )
          zslask=max(-89._realkind,zlakewateron*zlakeiceoff*( zunson + zunsoff* &
               (-zq*zqstar*zl2/zqam) ) + zlakeiceon + zlakewateroff)
          zq2m=zqsi+zqstar*zln2k + &
               (zqam-zqsi)*(1.0_realkind - exp(zslask) )
          !c
          !c     zy2 dummy if zunson=0
          !c
          zy=sqrt( 1.0_realkind  -9.0_realkind *zl2*zunson )
          !ps060928      zy=min( zy, sqrt( 8.0/z0lake(ihor,ilaketype) ) -1.0 )
          zy=min( zy, sqrt( 8.0_realkind /zrougl ) -1.0_realkind  )
          zy2=zunson*zy +zunsoff
          zy2=zy2*zlakewateron*zlakeiceoff + zlakeiceon + zlakewateroff
          !car010502 nwn beg
          !c  should it really be ztam/zqam in the calc. of wt2ms etc?
          !c
          zyn=sqrt( 1.0_realkind  -9.0_realkind *zln*zunson )
          zyn2=zyn*zunson + zunsoff
          zyn2=zyn2*zlakewateron*zlakeiceoff + zlakeiceon + zlakewateroff
          !c
          zpar=( (1.0_realkind  +zy2)/(1.0_realkind  +zyn2) )**2.0_realkind
          zpar=zpar*znrz2
          zpar=log(zpar)*zrkar
          !c
          t2mlake(ihor,ilaketype)=( zunson*( ztam -zthst*zpar) + &
               zunsoff*zt2m )*zlakewateron*zlakeiceoff + ztdumout*(zlakeiceon+zlakewateroff)
          q2mlake(ihor,ilaketype)=( zunson*( zqam -zqstar*zpar ) + &
               zunsoff*zq2m )*zlakewateron*zlakeiceoff + zqdumout*(zlakeiceon+zlakewateroff)
          !c
          !car010502 end

          !c
          !c---- new by anna to have values of wind over sea
          !c
          zl10=z10m*zmoin
          zln10=log(z10m/zrougl)
          zln10k=zln10*zrkar
          !c
          !c      zu10 etc dummy if zunson=1.0
          !c
          zslask=max(-89._realkind,zunson + zunsoff*(-zm*zust*zl10/zunlev))
          zu10=zust*zln10k + zunlev*(1.0_realkind -exp(zslask))
          zslask=max(-89._realkind,zunson + zunsoff*(-zm*zvst*zl10/zvnlev))
          zv10=zvst*zln10k + zvnlev*(1.0_realkind -exp(zslask))
          !c
          !c     zx10 dummy if zunson=0
          !c
          zx10=zunson*(1.0_realkind -15.0_realkind *zl10) + zunsoff
          zx10=sqrt(sqrt(zx10))
          !car010502 nwn beg
          !c
          zxn=zunson*( 1.0_realkind  -15.0_realkind *zln) + zunsoff
          zxn=zunson*sqrt( sqrt(zxn) )
          !c
          zpar= -log( z10m*(1.0_realkind  +zxn*zxn)/( zhnlev*(1.0_realkind  +zx10*zx10)) ) - &
               2.0_realkind *log( (1.0_realkind  +zxn)/(1.0_realkind  +zx10) ) + &
               2.0_realkind *( atan(zxn) -atan(zx10) )

          u10mlake(ihor,ilaketype)=zunson*(u(ihor,nlev) -zrkar*zust*zpar) &
               +zunsoff*zu10
          v10mlake(ihor,ilaketype)=zunson*( v(ihor,nlev) -zrkar*zvst*zpar) &
               +zunsoff*zv10
          momfulake(ihor,ilaketype)=zdens*u(ihor,nlev)/( zvel +zsecu ) * &
               zustar*zustar*zlakewateron*zlakeiceoff
          momfvlake(ihor,ilaketype)=zdens*v(ihor,nlev)/( zvel +zsecu ) * &
               zustar*zustar*zlakewateron*zlakeiceoff
          u10mlake(ihor,ilaketype)=u10mlake(ihor,ilaketype)*zlakewateron*zlakeiceoff & 
               + zudumout*zlakewateroff
          v10mlake(ihor,ilaketype)=v10mlake(ihor,ilaketype)*zlakewateron*zlakeiceoff &
               + zudumout*zlakewateroff
          !c----
          senflake(ihor,ilaketype)=zh*zlakewateron*zlakeiceoff
          latflake(ihor,ilaketype)=zeeff*zlakewateron*zlakeiceoff
          evaplake(ihor,ilaketype)=zevap*zlakewateron*zlakeiceoff
          !--------------------------------------------------------------

       end do ! loop over hor. grid points
    end do ! loop over lake types


    return
  end subroutine turbfluxlakewater

end module flake
