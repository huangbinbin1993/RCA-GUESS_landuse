subroutine gemini(klon,klat,klev)
  use boundaryRelaxation, only:npbpts
  use mod_implicit_solver
  use coupling
  use config
  use co2mod
  use comhkp
#ifdef USING_NETCDF
  use netcdfrestart
#else
  use restart
#endif
  use decomp
  use lateralBc
  use surface_bc
  use postprocess
  use escom
  use confys
  use ctun
  use climate
  use util
  use interpol
  use IO
  use physpar
  use mod_implicit_solver
  use sl2tim4
  use runtimeDiagnostics
  use accumulate
  use routmodvar
  use mod_diffh
  use domainmod
  use RCAdomainMod
  use timers
  use calendar
  use gcm
  use derived_types
  use filter_orog
  use ecoclimap
  use flake,only:lakeType
  use surface,only:surfVariables
  implicit none


  integer::halo

  integer,intent(in)::klon,klat,klev 
  integer::iacdg,iacdg2,ksvar

  integer::nstep,i,j,kk

#ifdef MPI_SRC
#include "mpif.h"
#endif


  integer::jx,js,jy,jk

  !     DYNAMIC VARIABLES
  type(atm)::phim,phiz,phip


  real(kind=realkind),pointer,dimension(:,:)::t2ms, t2mi
  real(kind=realkind),pointer,dimension(:,:)::senfs,senfi,latfs,latfi

  real(kind=realkind),allocatable,dimension(:,:,:)::dtdtph,dqdtph 
  real(kind=realkind),allocatable,dimension(:,:,:)::dsdtph
  real(kind=realkind),allocatable,dimension(:,:)::drolddt,dsolddt,accrunoff 
  real(kind=realkind),allocatable,dimension(:,:)::accrunoffopl,accrunofffor
  real(kind=realkind),allocatable,dimension(:,:)::accrunofflake,q2d
  real(kind=realkind),allocatable,dimension(:,:)::accprl,accprc,sswr
  
  real(kind=realkind),allocatable,dimension(:,:)::tsmax,tsmin,t2max,t2min,uv10max
  real(kind=realkind),allocatable,dimension(:,:)::cw2d,tswr,slwr,tlwr
  real(kind=realkind),allocatable,dimension(:,:)::accq2d  
  real(kind=realkind),allocatable,dimension(:,:)::acccw2d,accslwr,accsswr,acctlwr
  real(kind=realkind),allocatable,dimension(:,:)::acctswr,accsnow
  integer::stat

  real(kind=realkind),allocatable,dimension(:,:)::acctcov,acclcov,accmcov,acchcov
  real(kind=realkind),allocatable,dimension(:,:)::tcov,lcov,mcov,hcov
  real(kind=realkind),allocatable,dimension(:,:)::pr1h,prc1h,prl1h,maxdayprecip


  integer,parameter::esvars=27
  real(kind=realkind),allocatable,dimension(:,:,:)::extrem_surf
  !     1,2  :  t2mopsn   max min
  !     3,4  :  t2ms      max min
  !     5,6  :  t2mi      max min
  !     7,8  :  tsns      max min
  !     9,10 :  tc        max min
  !     11,12 :  tsc       max min
  !     13   :  uv10opsn  max
  !     14   :  uv10ms    max
  !     15   :  uv10mi    max
  !     16,17 :  tsnow     max min
  !     18   :  utot10ms  mean 3h
  !     19,20 :  rh2m      max min
  !     21-23 :  gust      max(estimated,lower,upper)
  !     24,25 :  t2mopsnsi max min
  !     26   :  uv10opsnsimax  max
  !     27   :  gustuc      max(uncorrected estimated)

  real(kind=realkind),allocatable,dimension(:,:)::uv10opsn_max,uv10ms_max
  real(kind=realkind),allocatable,dimension(:,:)::uv10mi_max,t2mopsn_max 
  real(kind=realkind),allocatable,dimension(:,:)::t2mopsn_min,t2ms_max 
  real(kind=realkind),allocatable,dimension(:,:)::t2ms_min,t2mi_max,t2mi_min 
  real(kind=realkind),allocatable,dimension(:,:)::tsns_max,tsns_min,tc_max,tc_min 
  real(kind=realkind),allocatable,dimension(:,:)::tsc_max,tsc_min 

  real(kind=realkind),allocatable,dimension(:,:)::slwrsea,sswrsea,slwrice,sswrice 

  !     prognostic variables for the land surface
  integer,parameter::ksvars=49
  real(kind=realkind),allocatable,dimension(:,:,:)::svarsm,svarsz,svarsp

  real(kind=realkind),allocatable,dimension(:,:)::sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn

  !     diagnostic variables for the land surface
  integer,parameter::msvars=121
  integer,allocatable::account_surf(:)
  real(kind=realkind),target,allocatable,dimension(:,:,:)::svar_surf
  real(kind=realkind),allocatable,dimension(:,:,:)::acum_surf
  real(kind=realkind),allocatable,target,dimension(:,:,:)::mean_surf

  !     ecoclimap physiography
  integer,parameter::meco=49
  real(kind=realkind),allocatable:: eco(:,:,:)

  real(kind=realkind),allocatable:: storage_w_ini(:,:)

  type(lakeType)::lakes

  real(kind=realkind),allocatable,dimension(:,:)::acct2ms,accq2m ,acct2mi
  real(kind=realkind),allocatable,dimension(:,:)::accsenfs,accsenfi,acclatfs
  real(kind=realkind),allocatable,dimension(:,:)::acclatfi
  real(kind=realkind),allocatable,dimension(:,:)::oaspr,oassnow
  real(kind=realkind),allocatable,dimension(:,:)::lastsnow
  real(kind=realkind),allocatable,dimension(:,:)::accslwrdown,accsswrdown
  real(kind=realkind),allocatable,dimension(:,:)::accslwrs,accsswrs
  real(kind=realkind),allocatable,dimension(:,:)::accslwri,accsswri
  real(kind=realkind),allocatable,dimension(:,:)::psl,accpsl
  real(kind=realkind),allocatable,dimension(:)  ::bfullOasis

  integer:: oasacccount
  real(kind=realkind),allocatable,dimension(:,:):: slwrs,sswrs  


  logical::  nlinit

  type(surfVariables)::surf_vars

  real(kind=realkind),pointer,dimension(:,:)::u10       ! 10-meter velocity in x-direction
  real(kind=realkind),pointer,dimension(:,:)::v10       ! 10-meter velocity in y-direction
  real(kind=realkind),pointer,dimension(:,:)::q2m       ! 2-meter specific humidity

  real(kind=realkind),pointer,dimension(:,:)::senf      ! sensible heat flux
  real(kind=realkind),pointer,dimension(:,:)::latf      ! latent heat flux
  real(kind=realkind),pointer,dimension(:,:)::momf      ! momentum flux


  real(kind=realkind),allocatable,dimension(:,:):: cov2d              ! 2-dimensional cloud cover
  real(kind=realkind),allocatable,dimension(:,:):: cwpath             ! vertically integrated cloud water cont
  real(kind=realkind),allocatable,dimension(:,:):: pblh               ! boundary layer height

  real(kind=realkind),allocatable,dimension(:,:):: albice             ! albedo over ice from pb
  real(kind=realkind),allocatable,dimension(:,:):: fwat               !heat water to ice in pb

  real(kind=realkind),allocatable,dimension(:,:):: lcounttice         !baltic sea counter (pb)


  real(kind=realkind),allocatable,dimension(:,:,:):: totcov ! total cloud cover
  real(kind=realkind),allocatable,dimension(:,:,:)::  cucov !cumulus cloud cover
  !     = 1 RCAO , = 2 lake PROBE, = 3 lake FLake, =0 else
  integer,allocatable,dimension(:,:)::guessid

  real(kind=realkind),allocatable,dimension(:,:,:)::dtdt_kf,dqdt_kf,dqcdt_kf
  real(kind=realkind),allocatable,dimension(:,:,:)::om_t1,om_t2,om_t3
  real(kind=realkind),allocatable,dimension(:,:,:)::oldcov,zvarcu,cwcu,rad_cloud
  real(kind=realkind),allocatable,dimension(:,:,:)::cwice,cwliq,fice_cj
  real(kind=realkind),allocatable,dimension(:,:)::  raincv_kf,snowcv_kf, umfb
  integer,allocatable,dimension(:,:)::nca,kf_ind,kf_top,kf_base,shal_cgj

  real(kind=realkind),allocatable,dimension(:,:):: ci2d,accci2d,qu2d,accqu2d,qv2d,accqv2d, vimfc2d,accvimfc2d

  real(kind=realkind),allocatable,dimension(:,:)::u10_reg,v10_reg,u10s_reg 
  real(kind=realkind),allocatable,dimension(:,:)::v10s_reg,u10i_reg,v10i_reg

  real(kind=realkind),allocatable,dimension(:,:)::tsea_clim, frice_clim     



  real(kind=realkind),allocatable,dimension(:,:):: ice_thic, accsunny 



  real(kind=realkind),allocatable,dimension(:,:)::tsclim_years

  ! For routing
  logical :: savedstate

  real(kind=realkind), allocatable, dimension(:,:,:):: big_runoff

  real(kind=realkind),allocatable,dimension(:,:):: runoff24

  real(kind=realkind),allocatable,dimension(:,:):: utot10ms, speed10
  

  real(kind=realkind),allocatable,dimension(:,:)::tsea_lake,iceconc_lake, &
       icethic_lake,lcount_lake,origlcounttice,origfrland
  real(kind=realkind) zfrsum,zfrsumi
  integer ilaketype

! OASIS3

  real(kind=realkind),allocatable,dimension(:,:)::momuo_fieldm1,momui_fieldm1,momvo_fieldm1,momvi_fieldm1
!  real(kind=realkind),allocatable,dimension(:,:)::u10s_fieldm1,u10i_fieldm1,v10s_fieldm1,v10i_fieldm1


  !     new rca3, transient+ensembles
  real(kind=realkind),allocatable,dimension(:,:)::accpr,lastpr,accsnoc,accsnol


  !     the following are the default values for the horizontal diffusion
  !     coefficients on the 0.5x0.5 degree grid.
  real(kind=realkind),parameter::ppdif52=5.0e+5_realkind  
  real(kind=realkind),parameter::ppdif54=3.5e+14_realkind 
  real(kind=realkind),parameter::ppdif56=1.0e+24_realkind 
  real(kind=realkind),parameter::ppdth5=0.5_realkind     

  !     definitions for 2d and 3d accumulated diagnostics

  integer,parameter:: macdg=8, macdg2=16
  integer::nwmodg(macdg) 
  integer:: nwmodg2(macdg2) 
  real(kind=realkind),allocatable:: sacdg(:,:,:,:)
  real(kind=realkind),allocatable:: sacdg2(:,:,:) 

  !     data nwmodg / 061, 210, 212, 211, 220, 221, 230, 231/
  !     F     1: *APRETA*    ACC. PRECIP. AT MODEL LEVELS
  !     F     2: *ADTP*      ACC. TOTAL TEND OF TEMPERATURE DUE TO PHSYICS
  !     F     3: *ADTPR*     ACC. TEND OF TEMPERATURE DUE TO RADIATION
  !     F     4: *ADTPT*     ACC. TEND OF TEMPERATURE DUE TO VERTICAL DIFFUSION
  !     F     5: *ADQP*      ACC. TOTAL TEND OF SPEC. HUM. DUE TO PHSYICS
  !     F     6: *ADQPT*     ACC. TEND OF SPEC. HUM DUE TO VERTICAL DIFFUSION
  !     F     7: *ADCP*      ACC. TOTAL TEND OF CLOUD CONDENSATE DUE TO PHSYICS
  !     F     8: *ADCPT*     ACC. TEND OF CLOUD CONDENSATE DUE TO VERTICAL DIFF
  !     
  !     LR    1 : *ACCPRC*    ACC. CONVECTIVE PRECIPITATION
  !     LR    2 : *ACCPRL*    ACC. STRATIFORM PRECIPITATION
  !     LR    3 : *ACSNOC*    ACC. CONVECTIVE SNOW PRECIPITATION
  !     LR    4 : *ACSNOL*    ACC. STRATIFORM SNOW PRECIPITATION
  !     LR    5 : *ASENF*     ACC. SENSIBLE HEAT FLUX
  !     LR    6 : *ALATF*     ACC. LATENT HEAT FLUX
  !     LR    7 : *AMOMFU*    ACC. MOMENTUM FLUX U-COMP
  !     LR    8 : *AMOMFV*    ACC. MOMENTUM FLUX V-COMP
  !     LR    9 : *AEVAP*     ACC. EVAPORATION
  !     LR    10: *ASLNET*    ACC. LONG-WAVE NET RADIATION AT SFC
  !     LR    11: *ASSNET*    ACC. SHORT-WAVE NET RADIATION AT SFC
  !     LR    12: *ASLDN*     ACC. LONG-WAVE DOWNWARD RADIATION AT SFC
  !     LR    13: *ASSDN*     ACC. SHORT-WAVE DOWNWARD RADIATION AT SFC
  !     LR    14: *ATLNET*    ACC. LONG-WAVE NET RADIATION AT TOA
  !     LR    15: *ATSNET*    ACC. SHORT-WAVE NET RADIATION AT TOA
  !     LR    16: *ATSDN*     ACC. SHORT-WAVE DOWNWARD RADIATION AT TOA
  !     F NOTE: The same hard coding also in PHTASK + POSTPP!
#ifdef TIMING
  call timer(0,'all')
#endif

! OASIS3
  if(use_oasis) then
    allocate(momuo_fieldm1(klon,klat),momui_fieldm1(klon,klat),  &
             momvo_fieldm1(klon,klat),momvi_fieldm1(klon,klat))	   
!  if(use_oasis_arctic) then
!    allocate(u10s_fieldm1(klon,klat),u10i_fieldm1(klon,klat),  &
!             v10s_fieldm1(klon,klat),v10i_fieldm1(klon,klat))
!  endif	   
  endif	   

  allocate(tsclim_years(klon,klat))
  allocate(guessid(klon,klat))
  allocate(acct2ms(klon,klat),accq2m(klon,klat) ,acct2mi(klon,klat))
  allocate(accsenfs(klon,klat),accsenfi(klon,klat),acclatfs(klon,klat))
  allocate(acclatfi(klon,klat))
  allocate(oaspr(klon,klat),oassnow(klon,klat))
  allocate(lastsnow(klon,klat))
  allocate(accslwrdown(klon,klat),accsswrdown(klon,klat))
  allocate(accslwrs(klon,klat),accsswrs(klon,klat))
  allocate(accslwri(klon,klat),accsswri(klon,klat))
  allocate(psl(klon,klat),accpsl(klon,klat))
  allocate(bfullOasis(klev+1))
  allocate(runoff24(klon,klat))
  allocate(account_surf(msvars),stat=stat)
  if(stat/=0)stop'account_surf allocation error'
  allocate(dtdtph(klon,klat,klev))
  allocate(dqdtph(klon,klat,klev))
  allocate(dsdtph(klon,klat,klev))
  allocate(drolddt(klon,klat))
  allocate(dsolddt(klon,klat)) 
  allocate(accrunoff(klon,klat))
  allocate(accrunoffopl(klon,klat))
  allocate(accrunofffor(klon,klat))
  allocate(accrunofflake(klon,klat))
  allocate(q2d(klon,klat))
  allocate(accprl(klon,klat))
  allocate(accprc(klon,klat))
  allocate(sswr(klon,klat))

  allocate(tsmax(klon,klat),tsmin(klon,klat),t2max(klon,klat),t2min(klon,klat),&
       uv10max(klon,klat) )                  
  allocate(cw2d(klon,klat),                    &
       tswr(klon,klat),                    &
       slwr(klon,klat),                    &
       tlwr(klon,klat))                    
  allocate(accq2d(klon,klat) ,                &
       acccw2d(klon,klat) ,                &
       accslwr(klon,klat) ,                &
       accsswr(klon,klat) ,                &
       acctlwr(klon,klat) ,                &
       acctswr(klon,klat))                 
!CUW110301beg
  allocate(acctcov(klon,klat) ,                 &
          acclcov(klon,klat) ,                 &
          accmcov(klon,klat) ,                 &
          acchcov(klon,klat) ,                 &
          tcov(klon,klat) ,                    &                 
          lcov(klon,klat) ,                    &                 
          mcov(klon,klat) ,                    &
          hcov(klon,klat)) 
!CUW110301end
  allocate(accsnow(klon,klat))                 
  allocate(extrem_surf(klon,klat,esvars))      
  allocate(uv10opsn_max(klon,klat),            &
       uv10ms_max(klon,klat) ,           &
       uv10mi_max(klon,klat) ,           &
       t2mopsn_max(klon,klat) ,           &
       t2mopsn_min(klon,klat) ,           &
       t2ms_max(klon,klat) ,           &
       t2ms_min(klon,klat) ,           &
       t2mi_max(klon,klat) ,           &
       t2mi_min(klon,klat) ,           &
       tsns_max(klon,klat) ,           &
       tsns_min(klon,klat) ,           &
       tc_max(klon,klat) ,           &
       tc_min(klon,klat) ,           &
       tsc_max(klon,klat) ,           &
       tsc_min(klon,klat))                 

  allocate(slwrsea(klon,klat),                 &
       sswrsea(klon,klat) ,                &
       slwrice(klon,klat) ,                &
       sswrice(klon,klat))                      

  allocate(svarsm(klon,klat,ksvars),           &
       svarsz(klon,klat,ksvars),           &
       svarsp(klon,klat,ksvars))
  allocate(sswdn(klon,klat),                   &
       slwdn(klon,klat) ,                  &
       accsswdn(klon,klat) ,                  &
       accslwdn(klon,klat) ,                  &
       tswdn(klon,klat) ,                  &
       acctswdn(klon,klat))
  allocate(svar_surf(klon,klat,msvars))
  allocate(acum_surf(klon,klat,msvars))
  allocate(mean_surf(klon,klat,msvars))
  allocate(eco(klon,klat,meco))
  allocate(storage_w_ini(klon,klat))

  allocate(lakes%prog_lakes(klon,klat,lakes%lake_types,lakes%lake_no_prog),           &
       lakes%tend_lakes(klon,klat,lakes%lake_types,lakes%lake_no_tend),           &
       lakes%diag_lakes(klon,klat,lakes%lake_types,lakes%lake_no_diag))           
  nullify(lakes%frac_lakes,lakes%depth_lakes)

  allocate( slwrs( klon,klat),sswrs( klon,klat)) 
  allocate(surf_vars%tsm(klon,klat))     
  allocate(surf_vars%tdm(klon,klat))     
  allocate(surf_vars%tsc(klon,klat))     
  allocate(surf_vars%tsea(klon,klat))     
  allocate(surf_vars%swm(klon,klat))     
  allocate(surf_vars%sdm(klon,klat))     
  allocate(surf_vars%swc(klon,klat))     
  allocate(surf_vars%snm(klon,klat))     
  allocate(surf_vars%snc(klon,klat))     
  allocate(surf_vars%rou(klon,klat))     
  allocate(surf_vars%roc(klon,klat))     
  allocate(surf_vars%alb(klon,klat))     

  allocate(surf_vars%fri(klon,klat))     
  allocate(surf_vars%frf(klon,klat))
  allocate(cov2d(klon,klat))
  allocate(cwpath(klon,klat))  
  allocate(pblh(klon,klat))    
  allocate(albice(klon,klat))  
  allocate(fwat(klon,klat))    
  allocate(lcounttice(klon,klat))
  allocate(totcov(klon,klat,klev),cucov(klon,klat,klev))
  allocate(dtdt_kf(klon,klat,klev))
  allocate(dqdt_kf(klon,klat,klev))
  allocate(dqcdt_kf(klon,klat,klev))
  allocate(om_t1(klon,klat,klev))
  allocate(om_t2(klon,klat,klev))
  allocate(om_t3(klon,klat,klev))
  allocate(oldcov(klon,klat,klev))
  allocate(zvarcu(klon,klat,klev))
  allocate(cwcu(klon,klat,klev))
  allocate(rad_cloud(klon,klat,klev))
  allocate(cwice(klon,klat,klev))
  allocate(cwliq(klon,klat,klev))
  allocate(fice_cj(klon,klat,klev))
  allocate(raincv_kf(klon,klat))
  allocate(snowcv_kf(klon,klat))
  allocate(umfb(klon,klat))
  allocate(nca(klon,klat))
  allocate(kf_ind(klon,klat))
  allocate(kf_top(klon,klat))
  allocate(kf_base(klon,klat))
  allocate(shal_cgj(klon,klat))
  allocate(ci2d(klon,klat),accci2d(klon,klat),qu2d(klon,klat),       &
       accqu2d(klon,klat),qv2d(klon,klat),accqv2d(klon,klat),         &
       vimfc2d(klon,klat),accvimfc2d(klon,klat))                       

  allocate(u10_reg(klon,klat),v10_reg(klon,klat), u10s_reg(klon,klat))
  allocate(v10s_reg(klon,klat),u10i_reg(klon,klat), v10i_reg(klon,klat))

  allocate(tsea_clim(klon,klat), frice_clim(klon,klat))


  allocate(ice_thic(klon,klat),accsunny(klon,klat))
  allocate(utot10ms(klon,klat),speed10(klon,klat))
  allocate(tsea_lake(klon,klat),iceconc_lake(klon,klat))
  allocate(icethic_lake(klon,klat),lcount_lake(klon,klat))
  allocate(origlcounttice(klon,klat),origfrland(klon,klat))
  allocate(accpr(klon,klat),lastpr(klon,klat),accsnoc(klon,klat),accsnol(klon,klat))
  allocate(maxdayprecip(klon,klat))
  allocate(pr1h(klon,klat))
  allocate(prc1h(klon,klat))
  allocate(prl1h(klon,klat))



  !     Check inputs

  nwmodg=(/061,210,212,211,220,221,230,231/)
  nwmodg2=(/063,062,078,079,122,121,124,125,057,112,111,115,116,114,113,117/)


  !         read namelist namtun and initiate physical constants
  call contun()
  call conphys(lecmwf)
  call tabdef()


  halo=1
  call initRca(klon,klat,klev)

!cau110622     moved from hlprog
  call getgcm(startTime,klon,klat,RCAdomain)

  !Have to decide if we are to read from a restart file here or else...

  nlinit = .true.


  current_time = startTime !startTime is read from namelist

  call read_nampos()

  call read_namvar(ksvar,iacdg,iacdg2)
  allocate(sacdg(klon,klat,klev,iacdg))
  allocate(sacdg2(klon,klat,iacdg2))

  phim = alloc_atm(klon,klat,klev,ksvar,current_time)
  phiz = alloc_atm(klon,klat,klev,ksvar,current_time)
  phip = alloc_atm(klon,klat,klev,ksvar,current_time)


  call initEcoclimap(current_time,klon,klat,meco,RCAdomain,lcounttice)
  call setInitialEcoclimap(current_time,klon,klat,meco,eco)


  if(use_oasis)then
    call read_namcoupling(current_time,ndtime)
    call couple_prism_define(Nseconds(ndtime),klon,klat,meco,rcaDomain,lcounttice,eco)
  endif

  ! For routing
  savedstate=.false.
  if(use_routing)then
     call initiate_routing(current_time%year,current_time%month,current_time%day,savedstate)
  endif

  call initLateral_bc(startTime,klon,klat,klev,ksvar, nwmosv,nlhumc,RCAdomain)
  
  if(.not.doRestart)then
     call setInitialCondition(phim,nlsimp,nlslan)
  endif
  call initSurface_bc(klon,klat,startTime,eco,meco,RCAdomain)

  call setSurface_bc(klon,klat,startTime,surf_vars%tsea,surf_vars%fri,svarsm(:,:,14),eco,&
       meco,.true.,RCAdomain,lcounttice)

  if(.not.doRestart)then
     call setInitialCondPhys(klon,klat,tsclim_years, &
          svar_surf,svarsm,lakes%prog_lakes,eco,&
          startTime,&
          RCAdomain%west,RCAdomain%south,RCAdomain%dlon, &
          RCAdomain%dlat,RCAdomain%polon,RCAdomain%polat)  
  endif


  if(doRestart)then
     call readDump(phim%ps,phim%lnps,phim%t,phim%u,phim%v,phim%q,phim%cw,phim%edot,phim%svar, &
          phiz%ps,phiz%lnps,phiz%t,phiz%u,phiz%v,phiz%q,phiz%cw,phiz%edot,phiz%svar, &
          phip%ps,phip%lnps,phip%t,phip%u,phip%v,phip%q,phip%cw,phip%edot,phip%svar, &
          dtdtph,dqdtph,dsdtph,drolddt,dsolddt,accsunny,&
          accrunoff,q2d,accrunoffopl,accrunofffor,accrunofflake,&
          accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr,slwr,tlwr,&
          accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri,&
          accq2d,acccw2d,accslwr,accsswr,acctlwr,acctswr,accsnow,&
          sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn,slwrsea,sswrsea,slwrice,sswrice,&
          svarsm,svarsz,svarsp,svar_surf,acum_surf,account_surf,mean_surf,&
          slwrs,sswrs,eco,surf_vars%tdm,surf_vars%tsc,surf_vars%tsea,surf_vars%swm,surf_vars%sdm,surf_vars%swc,surf_vars%snm,&
          surf_vars%snc,surf_vars%rou,surf_vars%roc,surf_vars%alb,surf_vars%fri,surf_vars%frf, &
          totcov,cucov,cov2d,cwpath,pblh,dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,&
          tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
          lcounttice,oldcov,zvarcu,cwcu,raincv_kf,snowcv_kf,umfb,rad_cloud,ci2d,accci2d,&
          cwice,cwliq,qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d,&
          u10_reg,v10_reg,u10s_reg,v10s_reg,u10i_reg,v10i_reg,&
          ice_thic,extrem_surf,sacdg,sacdg2,accpr,lastpr,&
          lakes%prog_lakes,lakes%tend_lakes,&
          accsnol,accsnoc,utot10ms,tsclim_years)

     if(use_oasis)then

        momuo_fieldm1(:,:)=mean_surf(:,:,41)/real(account_surf(41),realkind)
        momui_fieldm1(:,:)=mean_surf(:,:,42)/real(account_surf(42),realkind)
        momvo_fieldm1(:,:)=mean_surf(:,:,47)/real(account_surf(47),realkind)
        momvi_fieldm1(:,:)=mean_surf(:,:,48)/real(account_surf(48),realkind)

       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            momuo_fieldm1,momvo_fieldm1,momuo_fieldm,momvo_fieldm)

       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            momui_fieldm1,momvi_fieldm1,momui_fieldm,momvi_fieldm)

        netheato_fieldm(:,:)=(accsenfs(:,:)+acclatfs(:,:)+accslwrs(:,:))/real(account_surf(41),realkind)
        netheati_fieldm(:,:)=(accsenfi(:,:)+acclatfi(:,:)+accslwri(:,:))/real(account_surf(41),realkind)
        solaro_fieldm(:,:)=accsswrs(:,:)/real(account_surf(41),realkind)
        solari_fieldm(:,:)=accsswri(:,:)/real(account_surf(41),realkind)

        oemp_fieldm(:,:)=-1.0*(accprc(:,:)+accprl(:,:))/10800.

        call set_dqndt(svar_surf,t2ms)

        if(use_oasis_arctic) then
   
          t2ms_fieldm(:,:)=svar_surf(:,:,1)  ! use 2 instead !grid_averaged_T2m
          t2mi_fieldm = svar_surf(:,:,6)   ! 250 105 6	t2mi_i		lake_and/or_sea_ice_and_snow_averaged_T2m
          evas_fieldm = 0

          shfs_fieldm = accsenfs(:,:)/real(account_surf(41),realkind)    ! 
          shfi_fieldm = accsenfi(:,:)/real(account_surf(41),realkind)    ! 
          lhfs_fieldm = acclatfs(:,:)/real(account_surf(41),realkind)    ! 
          lhfi_fieldm = acclatfi(:,:)/real(account_surf(41),realkind)    ! 
          swd_fieldm = accsswdn(:,:)/real(account_surf(41),realkind)   ! accumulated
!          swns_fieldm = accsswr   !slfluxo_surf_sea_ice.incf90:  !     sswrs:  surface short wave radiation over sea
          swns_fieldm = accsswrs(:,:)/real(account_surf(41),realkind)   !slfluxo_surf_sea_ice.incf90:  !     sswrs:  surface short wave radiation over sea
!          lwd_fieldm = accslwdn 
          lwd_fieldm = accslwrdown 
          lwns_fieldm = accslwrs(:,:)/real(account_surf(41),realkind)
!          u10s_fieldm = svar_surf(:,:,17)     ! 250 105 17	u10ms_i		lake_and/or_sea_water_u10                       m_s^{-1}
!          u10i_fieldm = svar_surf(:,:,18)     ! 250 105 18	u10mi_i		lake_and/or_sea_ice_and_snow_averaged_u10
!          v10s_fieldm = svar_surf(:,:,23)     ! 250 105 23	v10ms_i		lake_and/or_sea_water_v10
!          v10i_fieldm = svar_surf(:,:,24)     ! 250 105 24	v10mi_i		lake_and/or_sea_ice_and_snow_averaged_v10    
          u10s_fieldm1 = svar_surf(:,:,13)     ! 250 105 17	u10ms_i		lake_and/or_sea_water_u10                       m_s^{-1}
          u10i_fieldm1 = svar_surf(:,:,13)     ! 250 105 18	u10mi_i		lake_and/or_sea_ice_and_snow_averaged_u10
          v10s_fieldm1 = svar_surf(:,:,19)     ! 250 105 23	v10ms_i		lake_and/or_sea_water_v10
          v10i_fieldm1 = svar_surf(:,:,19)     ! 250 105 24	v10mi_i		lake_and/or_sea_ice_and_snow_averaged_v10    
       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            u10s_fieldm1,v10s_fieldm1,u10s_fieldm,v10s_fieldm)
       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            u10i_fieldm1,v10i_fieldm1,u10i_fieldm,v10i_fieldm)

          rai_fieldm = accprl+accprc   ! mm/3 hours, includes rain and snow
          sno_fieldm = accsnow
          rqa_fieldm =  svar_surf(:,:,11)     !  250 105 11	q2ms_i		lake_and/or_sea_water_q2m
          clo_fieldm = acctcov
          slp_fieldm = 1000.0 
          amus_fieldm = svar_surf(:,:,41)     ! 250 105 41	momfus_i	lake_and/or_sea_water_uw N_m^{-2}
          amvs_fieldm = svar_surf(:,:,47)     ! 250 105 47	momfvs_i	lake_and/or_sea_water_vw N_m^{-2}
          amui_fieldm = svar_surf(:,:,42)     ! 250 105 42	momfui_i	lake_and/or_sea_ice_and_snow_averaged_uw  N_m^{-2}
          amvi_fieldm = svar_surf(:,:,48)     ! 250 105 48	momfvi_i	lake_and/or_sea_ice_and_snow_averaged_vw  N_m^{-2}
       endif    !  use_oasis_arctic



        call oasaccscratch(klon, klat,oasacccount, acct2ms, acct2mi, &
            accsenfs, accsenfi, acclatfs, acclatfi, accslwrs,accslwri,accsswrs,accsswri, &
            accslwrdown, accsswrdown, accq2m, accpsl)

      endif

  endif

  if(.not.dorestart)then
     accsunny = 0.0_realkind 
     totcov = 0.0_realkind 
     cucov  = 0.0_realkind 
     rad_cloud = 0.0_realkind 
     cov2d = 0.0_realkind 
     cwpath = 0.0_realkind 
     pblh = 0.0_realkind 
     accrunoff = 0.0_realkind 
     accrunoffopl = 0.0_realkind 
     accrunofffor = 0.0_realkind 
     accrunofflake = 0.0_realkind 
     accsunny = 0.0_realkind 
     accq2d = 0.0_realkind 
     acccw2d = 0.0_realkind 
     accci2d = 0.0_realkind 
     ci2d = 0.0_realkind 
     qu2d = 0.0_realkind 
     accqu2d = 0.0_realkind 
     qv2d = 0.0_realkind 
     accqv2d = 0.0_realkind 
     vimfc2d = 0.0_realkind 
     accvimfc2d = 0.0_realkind 
     slwr = 0.0_realkind 
     sswr = 0.0_realkind 
     tlwr = 0.0_realkind 
     tswr = 0.0_realkind 
     slwdn = 0.0_realkind 
     sswdn = 0.0_realkind 
     tswdn = 0.0_realkind 
     accslwr = 0.0_realkind 
     accsswr = 0.0_realkind 
     acctlwr = 0.0_realkind 
     acctswr = 0.0_realkind 
     accslwdn = 0.0_realkind 
     accsswdn = 0.0_realkind 
     acctswdn = 0.0_realkind 
     q2d = 0.0_realkind 
     cw2d = 0.0_realkind 
     ice_thic = 0.0_realkind 
     lastpr = 0.0_realkind 
     storage_w_ini = 0.0_realkind 
     svar_surf = 0._realkind 
     acum_surf = 0._realkind 
     extrem_surf = 0._realkind 
     account_surf = 0
     guessid = 0
     lakes%diag_lakes = 0._realkind 
     svarsz = svarsm
     svarsp = svarsm
!CUW110301beg
     acctcov = 0.0_realkind
     acclcov = 0.0_realkind
     accmcov = 0.0_realkind
     acchcov = 0.0_realkind
     tcov = 0.0_realkind
     lcov = 0.0_realkind
     mcov = 0.0_realkind
     hcov = 0.0_realkind
!CUW110301end

     if(use_oasis)then
        accslwrdown = 0.0_realkind 
        accsswrdown = 0.0_realkind 
        accslwrs = 0.0_realkind 
        accsswrs = 0.0_realkind 
        accslwri = 0.0_realkind 
        accsswri = 0.0_realkind 
        accsenfs = 0.0_realkind
        accsenfi = 0.0_realkind
        acclatfs = 0.0_realkind
        acclatfi = 0.0_realkind
        accprc   = 0.0_realkind
        accprl   = 0.0_realkind
        oaspr = 0.0_realkind 
        oassnow = 0.0_realkind 
        lastsnow = 0.0_realkind 
     endif
     sacdg2 = 0.0_realkind 
     sacdg = 0.0_realkind 

     phiz%ps   = phim%ps
     phiz%lnps = phim%lnps
     phiz%u = phim%u
     phiz%v = phim%v
     phiz%t = phim%t
     phiz%q = phim%q
     phiz%cw = phim%cw
     if (nlslan) then
        phiz%edot = phim%edot
     endif
     phiz%svar = phim%svar
     nca = 0
     raincv_kf = 0._realkind 
     snowcv_kf = 0._realkind 
     umfb = 0._realkind 
     shal_cgj = 0
     kf_ind = 0
     kf_base = 0
     kf_top = 0
     dtdt_kf = 0._realkind 
     dqdt_kf = 0._realkind 
     dqcdt_kf = 0._realkind 
     zvarcu = 0._realkind 
     cwcu = 0._realkind 
     oldcov = 0._realkind 
     om_t1 = 0._realkind 
     om_t2 = 0._realkind 
     om_t3 = 0._realkind 
     tsmax=-999.0_realkind 
     tsmin=999.0_realkind 
     t2max=-999.0_realkind 
     t2min=999.0_realkind 
     uv10max=-999.0_realkind 
     maxdayprecip=0._realkind 
     pr1h=0._realkind 
     prc1h=0._realkind 
     prl1h=0._realkind 
     do js=1,esvars/2
        extrem_surf(:,:,1+(js-1)*2) = -999.0_realkind 
        extrem_surf(:,:,2+(js-1)*2) = 999.0_realkind 
     enddo
     storage_w_ini=svar_surf(:,:,81)
  endif



  if (nlsimp)then
     call impini(klon,klat,klev,npbpts,RCAdomain%ahyb,RCAdomain%bhyb)
  endif
  call difhini(1,klon,klat,klev,real(2*Nseconds(ndtime),realkind),npbpts,&
       ksvar,RCAdomain%ahyb,RCAdomain%bhyb)


  nstep = 0
!  call timer(0,'All')
!  call timer(1,'All')
1000 do while( current_time <= stopTime ) !main time stepping loooop

     if(mype==0)then
        if(nstep<100.or.&
             (modulo(current_time%hour,6)==0.and.current_time%min==0.and.current_time%sec==0))then
           print *,' -------------------------------------- '
           write(6,'((a),i7,(a))')'nstep',nstep
           write(6,'(1x,(a),6i5)')'year,month,day,hour,min,sec=', &
                current_time%year,current_time%month,current_time%day,current_time%hour,&
                current_time%min,current_time%sec
        endif
     endif

!!     call mean_and_accumulation(svar_surf,account_surf,mean_surf,acum_surf)

     !     a few local arrayer for later use in gemini
     u10 => svar_surf(:,:,13)
     v10 => svar_surf(:,:,19)
     senf => svar_surf(:,:,25)
     latf => svar_surf(:,:,31)
     momf => svar_surf(:,:,49)

     if(nstep>0)then

       call hilot(klon,klat,svar_surf(:,:,63),tsmax,tsmin) !  tseff
       call hilot(klon,klat,svar_surf(:,:,1),t2max,t2min) !  t2m
       call hilot(klon,klat,svar_surf(:,:,3),extrem_surf(:,:,1),extrem_surf(:,:,2))
       call hilot(klon,klat,svar_surf(:,:,5), extrem_surf(:,:,3),extrem_surf(:,:,4))!  t2ms
       call hilot(klon,klat,svar_surf(:,:,6), extrem_surf(:,:,5),extrem_surf(:,:,6)) !  t2mi   
       call hilot(klon,klat,svarsz(:,:,1),extrem_surf(:,:,7),extrem_surf(:,:,8))  !  tsns
       call hilot(klon,klat,svarsz(:,:,2),extrem_surf(:,:,9),extrem_surf(:,:,10))!  tc
       call hilot(klon,klat,svarsz(:,:,3), extrem_surf(:,:,11),extrem_surf(:,:,12))!  tsc
       call hilot(klon,klat,svarsz(:,:,4), extrem_surf(:,:,16),extrem_surf(:,:,17))!  tsnow
       call hilorh(klon,klat,svar_surf(:,:,83), extrem_surf(:,:,19),extrem_surf(:,:,20))!  rh2m
       call higust(klon,klat,svar_surf(:,:,92),svar_surf(:,:,93),svar_surf(:,:,94),svar_surf(:,:,120),&
            extrem_surf(:,:,21),extrem_surf(:,:,22),extrem_surf(:,:,23),extrem_surf(:,:,27)) ! wind gusts
       call hilot(klon,klat,svar_surf(:,:,107),extrem_surf(:,:,24),extrem_surf(:,:,25))!  t2opsnsi
       call hiuv10(klon,klat,svar_surf(:,:,15),svar_surf(:,:,21),extrem_surf(:,:,13),current_time)!  u10opsn  v10opsn
       call hiuv10(klon,klat,svar_surf(:,:,17),svar_surf(:,:,23),extrem_surf(:,:,14),current_time)!  u10ms  v10ms
       call hiuv10(klon,klat,svar_surf(:,:,18),svar_surf(:,:,24),extrem_surf(:,:,15),current_time)!  u10mi  v10mi
       call hiuv10(klon,klat,svar_surf(:,:,110),svar_surf(:,:,111),extrem_surf(:,:,26),current_time)!  u10opsnsi  v10opsnsi
       call hiuv10(klon,klat,u10,v10,uv10max,current_time)
       do jy=1,klat
          do jx=1,klon
             if ( svar_surf(jx,jy,17) < 90._realkind  ) then
                utot10ms(jx,jy) = utot10ms(jx,jy) + sqrt(svar_surf(jx,jy,17)**2.0_realkind+&
                     svar_surf(jx,jy,23)**2.0_realkind)
             elseif ( svar_surf(jx,jy,18) < 90._realkind  ) then
                utot10ms(jx,jy) = utot10ms(jx,jy) + sqrt(svar_surf(jx,jy,18)**2.0_realkind+ &
                     svar_surf(jx,jy,24)**2.0_realkind)
             else
                utot10ms(jx,jy) = -1._realkind 
             endif
          enddo
       enddo

       call cloud_cover(klon,klat,klev,phim%ps,cucov,totcov, &
                        rad_cloud,tcov,lcov,mcov,hcov)

       call accumulate1(klon,klat,q2d,accq2d)
       call accumulate1(klon,klat,cw2d,acccw2d)
       call accumulate1(klon,klat,slwr,accslwr)
       call accumulate1(klon,klat,sswr,accsswr)
       call accumulate1(klon,klat,tlwr,acctlwr)
       call accumulate1(klon,klat,tswr,acctswr)
       call accumulate1(klon,klat,slwdn,accslwdn)
       call accumulate1(klon,klat,sswdn,accsswdn)
       call accumulate1(klon,klat,tswdn,acctswdn)
       call accumulate1(klon,klat,ci2d,accci2d)
       call accumulate1(klon,klat,qu2d,accqu2d)
       call accumulate1(klon,klat,qv2d,accqv2d)
       call accumulate1(klon,klat,vimfc2d,accvimfc2d)
!CUW110301beg
     call accumulate1(klon,klat,tcov,acctcov)
     call accumulate1(klon,klat,lcov,acclcov)
     call accumulate1(klon,klat,mcov,accmcov)
     call accumulate1(klon,klat,hcov,acchcov)
!CUW110301end

!uh
! calculate maximum hourly precipitation rate

     call accumulate1(klon,klat,sacdg2(:,:,1),prc1h)
     call accumulate1(klon,klat,sacdg2(:,:,2),prl1h)
     pr1h=prl1h+prc1h
     if(mod(nstep, 3600/Nseconds(ndtime))==0)then ! full hour
        call accmean1(klon,klat,pr1h,3600/Nseconds(ndtime))
        maxdayprecip=max(maxdayprecip,pr1h)
        pr1h=0._realkind
        prc1h=0._realkind
        prl1h=0._realkind
     endif

       call acc_surf(klon,klat,msvars,svar_surf,acum_surf)
       call acc4mean_surf(klon,klat,account_surf,msvars, svar_surf,mean_surf)

        if(use_oasis)then
           do jk=1,klev
              bfullOasis(jk) = 0.5_realkind *( RCAdomain%bhyb(jk) + RCAdomain%bhyb(jk+1) )
           enddo
           !     computation of mean sea-level pressure
           call pslcom(nlphys,phim%ps,RCAdomain%fis,phim%t,svar_surf(:,:,63), &
                psl,klon,klat,klev,bfullOasis)

           t2ms=>svar_surf(:,:,5)
           t2mi=>svar_surf(:,:,6)
           senfs=>svar_surf(:,:,29)
           senfi=>svar_surf(:,:,30)
           latfs=>svar_surf(:,:,35)
           latfi=>svar_surf(:,:,36)
           q2m=>svar_surf(:,:,7)
           call oasacc(klon,klat,oasacccount,t2ms,acct2ms,t2mi,acct2mi,  &
                senfs,accsenfs,senfi,accsenfi,latfs,acclatfs,latfi,     &
                acclatfi,slwrs,accslwrs,slwrice,accslwri,sswrs,accsswrs,sswrice,accsswri, &
                slwdn,accslwrdown, &
                sswdn,accsswrdown,q2m,accq2m,psl,accpsl)
        endif

       call accumulate_post_process_surface(nstep,klon,klat,msvars,account_surf,mean_surf,extrem_surf,utot10ms)
       call accumulate_post_process(nstep,klon,klat,accsswr,accslwr,acctswr,acctlwr,accslwdn,accsswdn,acctswdn,accq2d,accci2d,&
            acccw2d,accqu2d,accqv2d,accvimfc2d,acctcov,acclcov,accmcov,acchcov)
       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            u10,v10,u10_reg,v10_reg)


       if(mod(nstep, 21600/Nseconds(ndtime))==0)then !number of timesteps to 6h
          if(use_routing)then
             runoff24 = runoff24 + accrunoff
             if( current_time%hour==18 ) then
                allocate(big_runoff(klon_global,klat_global,1),stat=stat)
                if(stat/=0)then
                   stop 'memory allocation error gemini big_runoff'
                endif
                call colfld(0,big_runoff,runoff24,klon,klat)
                if ( mype==0 ) then
                   call river_routing(current_time%year,current_time%month,current_time%day,&
                        big_runoff,runoff_bound)
                   savedstate = .false.
                endif
                deallocate(big_runoff)
                runoff24 = 0._realkind 
             endif
          endif
       endif

     endif  ! nstep>0

     call writeDump(current_time, & !check is made inside writeDump wheter to write or not
          phim%ps,phim%lnps,phim%t,phim%u,phim%v,phim%q,phim%cw,phim%edot,phim%svar, &
          phiz%ps,phiz%lnps,phiz%t,phiz%u,phiz%v,phiz%q,phiz%cw,phiz%edot,phiz%svar, &
          phip%ps,phip%lnps,phip%t,phip%u,phip%v,phip%q,phip%cw,phip%edot,phip%svar, &
          dtdtph,dqdtph,dsdtph,drolddt,dsolddt,accsunny,&
          accrunoff,q2d,accrunoffopl,accrunofffor,accrunofflake,&
          accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr,slwr,tlwr,&
          accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri,&
          accq2d,acccw2d,accslwr,accsswr,acctlwr,acctswr,accsnow,&
          sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn,slwrsea,sswrsea,slwrice,sswrice,&
          svarsm,svarsz,svarsp,svar_surf,acum_surf,account_surf,mean_surf,&
          slwrs,sswrs,eco, &
          surf_vars%tdm,surf_vars%tsc,surf_vars%tsea,surf_vars%swm,surf_vars%sdm,surf_vars%swc,&
          surf_vars%snm,surf_vars%snc,surf_vars%rou,surf_vars%roc,surf_vars%alb,surf_vars%fri,surf_vars%frf, &
          totcov,cucov,cov2d,cwpath,pblh,dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,&
          tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
          lcounttice,oldcov,zvarcu,cwcu,raincv_kf,snowcv_kf,umfb,rad_cloud,ci2d,accci2d,&
          cwice,cwliq,qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d,&
          u10_reg,v10_reg,u10s_reg,v10s_reg,u10i_reg,v10i_reg,&
          ice_thic,extrem_surf,sacdg,sacdg2,accpr,lastpr,&
          lakes%prog_lakes,lakes%tend_lakes,&
          accsnol,accsnoc,utot10ms,tsclim_years)

     call postproc(klon,klat,klev,ksvar,current_time,phim,  &
          surf_vars,&
          totcov ,cucov  ,cw2d  ,acccw2d  , pblh,          &
          nwmosv,sacdg,iacdg,nwmodg,sacdg2,iacdg2, &
          ci2d,accci2d,rad_cloud ,accqu2d,accqv2d,accvimfc2d, &
          accsunny, ice_thic,  &
          accrunoff, accprl, accprc ,tsmax,tsmin,t2max,t2min,maxdayprecip,     &
          uv10max,accq2d,accslwr,accsswr,acctlwr,acctswr,     &
          accslwdn,accsswdn,acctswdn  ,accsnow  ,&
          accsnoc,accsnol,   &
          tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov, &
          lcounttice ,msvars, svar_surf ,ksvars, svarsz ,esvars, extrem_surf, &
          acum_surf, mean_surf,meco, eco ,lakes)

     if(use_oasis)then
       if(nstep > 0) then
       if(mod(nstep,10800/Nseconds(ndtime))==0)then  !coupling every 3 hours
           do kk=1,msvars
           call accmean_surf(klon, klat, account_surf, msvars, kk, mean_surf)
           enddo

        momuo_fieldm1(:,:)=mean_surf(:,:,41)
        momui_fieldm1(:,:)=mean_surf(:,:,42)
        momvo_fieldm1(:,:)=mean_surf(:,:,47)
        momvi_fieldm1(:,:)=mean_surf(:,:,48)

       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            momuo_fieldm1,momvo_fieldm1,momuo_fieldm,momvo_fieldm)

       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            momui_fieldm1,momvi_fieldm1,momui_fieldm,momvi_fieldm)

        if(use_oasis_arctic) then
        write(6,*)'nstep,momuo_fieldm(10,10) ',nstep,momuo_fieldm(10,10)
        endif  ! use_oasis_arctic

           do kk=1,msvars
           call accscratch4mean_surf(klon, klat,account_surf, msvars, kk, mean_surf)
           enddo

           call  oasaccmean(klon,klat,oasacccount,acct2ms,&
             acct2mi,accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri, &
             accq2m,accpsl,accslwrdown,accsswrdown )

        netheato_fieldm(:,:)=accsenfs(:,:)+acclatfs(:,:)+accslwrs(:,:)
        netheati_fieldm(:,:)=accsenfi(:,:)+acclatfi(:,:)+accslwri(:,:)
        solaro_fieldm(:,:)=accsswrs(:,:)
        solari_fieldm(:,:)=accsswri(:,:)
        call set_dqndt(svar_surf,t2ms)
!wangsy        call set_dqndt(svar_surf,t2ms,klon,klat)
        oemp_fieldm(:,:)=-1.0*(accprc(:,:)+accprl(:,:))/10800.

      if(use_oasis_arctic) then
!rd20120109:
        shfs_fieldm = accsenfs(:,:)    ! 
        shfi_fieldm = accsenfi(:,:)    ! 
        lhfs_fieldm = acclatfs(:,:)    ! 
        lhfi_fieldm = acclatfi(:,:)    ! 
        !t2ms_fieldm(:,:)=svar_surf(:,:,5)
        t2ms_fieldm(:,:)=svar_surf(:,:,1)   ! use 2 instead !grid_averaged_T2m
        t2mi_fieldm = svar_surf(:,:,6)      ! 250 105 6	t2mi_i		lake_and/or_sea_ice_and_snow_averaged_T2m
        evas_fieldm =0
!        swd_fieldm = accsswdn  
        swd_fieldm = accsswrdown 
!        swns_fieldm = accsswr  
        swns_fieldm = accsswrs  
!        lwd_fieldm = accslwdn 
        lwd_fieldm = accslwrdown 
        lwns_fieldm = accslwrs
        !u10s_fieldm = svar_surf(:,:,17)     ! 250 105 17	u10ms_i		lake_and/or_sea_water_u10
        !u10i_fieldm = svar_surf(:,:,18)     ! 250 105 18	u10mi_i		lake_and/or_sea_ice_and_snow_averaged_u10
        !v10s_fieldm = svar_surf(:,:,23)     ! 250 105 23	v10ms_i		lake_and/or_sea_water_v10
        !v10i_fieldm = svar_surf(:,:,24)     ! 250 105 24	v10mi_i		lake_and/or_sea_ice_and_snow_averaged_v10    

        ! svar_surf  instanteneous
        u10s_fieldm1 = svar_surf(:,:,13)     ! 250 105 17	u10ms_i		lake_and/or_sea_water_u10                       m_s^{-1}
        u10i_fieldm1 = svar_surf(:,:,13)     ! 250 105 18	u10mi_i		lake_and/or_sea_ice_and_snow_averaged_u10
        v10s_fieldm1 = svar_surf(:,:,19)     ! 250 105 23	v10ms_i		lake_and/or_sea_water_v10
        v10i_fieldm1 = svar_surf(:,:,19)     ! 250 105 24	v10mi_i		lake_and/or_sea_ice_and_snow_averaged_v10    
       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            u10s_fieldm1,v10s_fieldm1,u10s_fieldm,v10s_fieldm)
       call turn_wind(klon,klat,&
            RCAdomain%west,RCAdomain%south,RCAdomain%polon,RCAdomain%polat,RCAdomain%dlon,RCAdomain%dlat,&
            u10i_fieldm1,v10i_fieldm1,u10i_fieldm,v10i_fieldm)

        rai_fieldm = accprl+accprc
        sno_fieldm = accsnow
        rqa_fieldm =  svar_surf(:,:,11)     !  250 105 11	q2ms_i		lake_and/or_sea_water_q2m
        clo_fieldm = tcov
        slp_fieldm = psl    ! accpsl 
        amus_fieldm = svar_surf(:,:,41)     ! 250 105 41	momfus_i	lake_and/or_sea_water_uw
        amvs_fieldm = svar_surf(:,:,47)     ! 250 105 47	momfvs_i	lake_and/or_sea_water_vw
        amui_fieldm = svar_surf(:,:,42)     ! 250 105 42	momfui_i	lake_and/or_sea_ice_and_snow_averaged_uw 
        amvi_fieldm = svar_surf(:,:,48)     ! 250 105 48	momfvi_i	lake_and/or_sea_ice_and_snow_averaged_vw
        write(6,*)'nstep,t2ms_fieldm(10,10) ',nstep,t2ms_fieldm(10,10)
      endif   ! use_oasis_arctic

           call oasaccscratch(klon, klat,oasacccount, acct2ms, acct2mi, &
            accsenfs, accsenfi, acclatfs, acclatfi, accslwrs,accslwri,accsswrs,accsswri, &
            accslwrdown, accsswrdown, accq2m, accpsl)
       endif
       endif   ! step > 0

!       if(nstep < (nnstep-2)) then
!       call couple_prism_send(klon,klat,Nseconds(ndtime))
!       endif

    endif      ! use_oasis


!!     call scratch_mean_and_accumulation(account_surf,mean_surf,acum_surf)

     call scratch_post_process_surface(nstep,klon,klat,msvars,&
          account_surf,mean_surf,acum_surf,extrem_surf,utot10ms)
     call scratch_post_process(nstep,klon,klat,accsswr,accslwr,acctswr, &
          acctlwr,accslwdn,accsswdn,acctswdn,&
          accq2d,accci2d,acccw2d,accqu2d,accqv2d,accvimfc2d,accsunny, &
          accrunoff,accrunoffopl,accrunofffor,accrunofflake,&
          accsnow,accsnol,accsnoc,accprl,lastpr,accprc,tsmax,t2max,tsmin,t2min,uv10max,maxdayprecip,&
          acctcov,acclcov,accmcov,acchcov)

     call gco2(current_time,nlinit)
     call setEcoclimap(current_time,klon,klat,meco,eco,RCAdomain,lcounttice)


     call sl2tim(nstep,klon,klat,klev,ksvar,RCAdomain,          &
          phim,phiz,phip, &
          surf_vars,&
          totcov,cucov,cov2d,cwpath, &
          pblh,sacdg,iacdg,sacdg2,iacdg2,  &
          nlinit,                                                   &
          drolddt, dsolddt, accsunny,                               &
          accrunoff,                                                &
          accrunoffopl,accrunofffor,accrunofflake,                  &
          q2d, &
          accprl, accprc,                                           &
          !
          slwr, sswr, tlwr, tswr,                                   &
          slwrs, sswrs,                                             &
          slwrice,sswrice,                                          &
          slwdn,sswdn,tswdn,                                        &
          !
          svarsz, svarsp, ksvars,                           &
          ice_thic,                                                 &
          guessid,                                                  &
          lcounttice,                                               &
          accsnow,accsnoc,accsnol,                                  &
          msvars,svar_surf,                                         &
          eco,meco,                                                 &
          tsclim_years,                                             &
          dtdtph, dqdtph, dsdtph,&
          nca,dtdt_kf,dqdt_kf,                                &
          dqcdt_kf,raincv_kf,snowcv_kf,                              &
          umfb,shal_cgj,kf_ind,                                      &
          kf_base,                                                   &
          kf_top,                                                    &
          oldcov,                                                    &
          zvarcu,cwcu,lakes)

     call setLateral_bc(klon,klat,klev,ksvar,phip,nlslan,RCAdomain)

     !     prevent negative values of cloud water and extra scalars
     phip%cw = max(phip%cw,0.0_realkind )
     phip%q = max(phip%q,1.e-7_realkind )
     phip%svar = max(phip%svar,0.0_realkind )

     if ( (nstep==0.and..not.doRestart) .or. epsn<=0.0_realkind ) then
        phim%ps = phiz%ps
        phim%lnps = phiz%lnps
        phim%u = phiz%u
        phim%v = phiz%v
        phim%t = phiz%t
        phim%q = phiz%q
        phim%cw = phiz%cw
        phim%svar = phiz%svar
        svarsm = svarsz
        phim%edot = phiz%edot
     else
        phim%ps   = (1.0_realkind -2.0_realkind *epsn)*phiz%ps+epsn*(phim%ps+phip%ps)
        phim%lnps = log( phim%ps )
        phim%u = (1.0_realkind -2.0_realkind *epsn)*phiz%u+epsn*(phim%u+phip%u)
        phim%v = (1.0_realkind -2.0_realkind *epsn)*phiz%v+epsn*(phim%v+phip%v)
        phim%t = (1.0_realkind -2.0_realkind *epsn)*phiz%t+epsn*(phim%t+phip%t)
        phim%q = (1.0_realkind -2.0_realkind *epsn)*phiz%q+epsn*(phim%q+phip%q)
        phim%cw = phiz%cw
        phim%svar = phiz%svar
        svarsm = svarsz
        phim%edot = phiz%edot
     endif
     accpr = accprl + accprc - lastpr
     lastpr = accprl + accprc

     if(use_oasis)then
        oaspr = oaspr + accpr
        oassnow = oassnow + accsnow - lastsnow
        lastsnow = accsnow
     endif

     svar_surf(:,:,81)=svar_surf(:,:,81)-storage_w_ini

     !         compute and print diagnostics at new timestep
     if(nlstat) then
        if((modulo(current_time%hour,6)==0.and.current_time%min==0.and.current_time%sec==0).or.nstep<=10)then
           call diagnos(klon,klat,klev,real(2*Nseconds(ndtime),realkind),&
                phip%ps,phip%u,phip%v,phip%t,phip%q,phip%cw,phiz%ps,&
                RCAdomain%fis,RCAdomain%ahyb,RCAdomain%bhyb)
        endif
     endif

     phiz%ps = phip%ps
     phiz%lnps = phip%lnps
     phiz%u = phip%u
     phiz%v = phip%v
     phiz%t = phip%t
     phiz%q = phip%q
     phiz%cw = phip%cw
     phiz%svar = phip%svar
     phiz%timestamp = phip%timestamp
     !     gjkf   redistribute running mean omega and move cloud fields
     !     gjkf   into oldcov
     if(nlslan)then
        oldcov = totcov
     endif


     if (nlphys) then
        svarsm = svarsz
        svarsz = svarsp
     endif
!wangsy     if(use_oasis)then
!        call couple_prism_recv(klon,klat,Nseconds(ndtime))    
!     endif

     call setSurface_bc(klon,klat,current_time,surf_vars%tsea,surf_vars%fri,svarsm(:,:,14),&
          eco,meco,.false.,RCAdomain,lcounttice)

     if(use_oasis)then
       if(nstep < (nnstep-2)) then
#if defined(OASIS3)
        seconds_since_begin=nstep*Nseconds(ndtime)
#endif
       call couple_prism_send(klon,klat,Nseconds(ndtime))
!       call couple_prism_recv(klon,klat,Nseconds(ndtime))
       call couple_prism_recv(klon,klat,surf_vars%tsea,surf_vars%fri,svarsm(:,:,14),svar_surf(:,:,62),Nseconds(ndtime))
      endif
     endif


     !This loop should be in flake?
     do jy=1,klat
        do jx=1,klon
           if(lcounttice(jx,jy)>2.5_realkind ) then
              tsea_lake(jx,jy)=0._realkind 
              icethic_lake(jx,jy)=0._realkind 
              iceconc_lake(jx,jy)=0._realkind 
              zfrsum = 0.0_realkind !sum(eco(jx,jy,(43+1):(43+laketypes)))
              do ilaketype=1,lakes%lake_types
                 !if(eco(jx,jy,43+ilaketype)>=0.01_realkind )then this is true by construct
                 tsea_lake(jx,jy)=tsea_lake(jx,jy)+eco(jx,jy,43+ilaketype)*lakes%prog_lakes(jx,jy,ilaketype,1)            
                 icethic_lake(jx,jy)=icethic_lake(jx,jy)+eco(jx,jy,43+ilaketype)*lakes%prog_lakes(jx,jy,ilaketype,9)            
                 if(lakes%prog_lakes(jx,jy,ilaketype,9)>0._realkind )then   
                    iceconc_lake(jx,jy)=iceconc_lake(jx,jy)+eco(jx,jy,43+ilaketype)
                 endif
                 zfrsum=zfrsum+eco(jx,jy,43+ilaketype)
                 !endif
              enddo
              if(abs(zfrsum)<0.0000000001_realkind)then
                 do ilaketype=1,lakes%lake_types
                    print *,eco(jx,jy,43+ilaketype)
                 enddo
                 stop 'zfrsum=0 in gemini'
              endif
              zfrsumi = 1.0_realkind /zfrsum
              surf_vars%tsea(jx,jy) = tsea_lake(jx,jy)*zfrsumi
              ice_thic(jx,jy) = icethic_lake(jx,jy)*zfrsumi
              surf_vars%fri(jx,jy) = iceconc_lake(jx,jy)*zfrsumi
              !     put prognostic ice temperature to lake surface temp
              svarsm(jx,jy,14) = tsea_lake(jx,jy)*zfrsumi
              svarsp(jx,jy,14) = tsea_lake(jx,jy)*zfrsumi
              svarsz(jx,jy,14) = tsea_lake(jx,jy)*zfrsumi
              !     put prognostic snow-ice temperature to lake surface temp
              svarsm(jx,jy,39) = tsea_lake(jx,jy)*zfrsumi
              svarsp(jx,jy,39) = tsea_lake(jx,jy)*zfrsumi
              svarsz(jx,jy,39) = tsea_lake(jx,jy)*zfrsumi
           endif
        enddo
     enddo
     !end flake

     call ficefun(klon,klat,klev,phip%t,fice_cj)
     cwice=phip%cw*fice_cj
     cwliq=phip%cw*(1._realkind -fice_cj)

     call intvert2(klon,klat,klev,phip%q,phip%ps,RCAdomain%ahyb,RCAdomain%bhyb,q2d)
     call intvert2(klon,klat,klev,cwliq,phip%ps,RCAdomain%ahyb,RCAdomain%bhyb,cw2d)
     call intvert2(klon,klat,klev,cwice,phip%ps,RCAdomain%ahyb,RCAdomain%bhyb,ci2d)
     call intvert_moist(klon,klat,klev,phip%q,phip%u,phip%v,phip%ps,qu2d,qv2d,vimfc2d,RCAdomain%ahyb,RCAdomain%bhyb)

     if(use_oasis)then   


     endif

     current_time = current_time + ndtime 

     nlinit = .false. 
     nstep = nstep + 1

  enddo !end main time loop

!  call timer(2,'All')

!this has been called in hlprog.F90, so it is not required in here anymore.
!  if(use_oasis)then
!     
!     call couple_prism_finalize
!  endif

  deallocate(account_surf)
  deallocate(dtdtph)
  deallocate(dqdtph)
  deallocate(dsdtph)
  deallocate(drolddt)
  deallocate(dsolddt) 
  deallocate(accrunoff)
  deallocate(accrunoffopl)
  deallocate(accrunofffor)
  deallocate(accrunofflake)
  deallocate(q2d)
  deallocate(accprl)
  deallocate(accprc)
  deallocate(sswr)

  deallocate(tsmax,tsmin,t2max,t2min, uv10max )                  
  deallocate(cw2d,tswr,slwr,tlwr)                    
  deallocate(accq2d,acccw2d,accslwr,accsswr,acctlwr,acctswr)                 
  deallocate(accsnow)                 
  deallocate(extrem_surf)
  deallocate(uv10opsn_max,uv10ms_max,uv10mi_max,t2mopsn_max, &
       t2mopsn_min,t2ms_max,t2ms_min,t2mi_max,t2mi_min,&
       tsns_max,tsns_min,tc_max,tc_min,tsc_max,tsc_min)                 
  deallocate(slwrsea,sswrsea,slwrice,sswrice)                      
  deallocate(svarsm,svarsz,svarsp)
  deallocate(sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn)
  deallocate(svar_surf)
  deallocate(acum_surf)
  deallocate(mean_surf)
  deallocate(eco)
  deallocate(storage_w_ini)
  deallocate(lakes%prog_lakes, lakes%tend_lakes,lakes%diag_lakes)
  deallocate(slwrs,sswrs)
  deallocate(surf_vars%tsm)     
  deallocate(surf_vars%tdm)     
  deallocate(surf_vars%tsc)     
  deallocate(surf_vars%tsea)     
  deallocate(surf_vars%swm)     
  deallocate(surf_vars%sdm)     
  deallocate(surf_vars%swc)     
  deallocate(surf_vars%snm)     
  deallocate(surf_vars%snc)     
  deallocate(surf_vars%rou)     
  deallocate(surf_vars%roc)     
  deallocate(surf_vars%alb)     
  deallocate(surf_vars%fri)     
  deallocate(surf_vars%frf)
  deallocate(cov2d)  
  deallocate(cwpath)  
  deallocate(pblh)    
  deallocate(albice)  
  deallocate(fwat)    
  deallocate(lcounttice)
  deallocate(totcov,cucov)
  deallocate(dtdt_kf)
  deallocate(dqdt_kf)
  deallocate(dqcdt_kf)
  deallocate(om_t1)
  deallocate(om_t2)
  deallocate(om_t3)
  deallocate(oldcov)
  deallocate(zvarcu)
  deallocate(cwcu)
  deallocate(rad_cloud)
  deallocate(cwice)
  deallocate(cwliq)
  deallocate(fice_cj)
  deallocate(raincv_kf)
  deallocate(snowcv_kf)
  deallocate(umfb)
  deallocate(nca)
  deallocate(kf_ind)
  deallocate(kf_top)
  deallocate(kf_base)
  deallocate(shal_cgj)
  deallocate(ci2d,accci2d,qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d)
  deallocate(u10_reg,v10_reg, u10s_reg)
  deallocate(v10s_reg,u10i_reg, v10i_reg)
  deallocate(tsea_clim, frice_clim)
  deallocate(ice_thic,accsunny)
  deallocate(utot10ms,speed10)
  deallocate(tsea_lake,iceconc_lake)
  deallocate(icethic_lake,lcount_lake)
  deallocate(tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov)
  deallocate(origlcounttice,origfrland)
  deallocate(accpr,lastpr,accsnoc,accsnol)
  deallocate(sacdg)
  deallocate(sacdg2)
  deallocate(maxdayprecip)                  
  deallocate(pr1h)                  
  deallocate(prc1h)                  
  deallocate(prl1h)                  



  return
end subroutine gemini


