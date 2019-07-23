module physpar
  use decomp
  use rcadomainMod
  use comhkp,only:current_time,nlphys
  implicit none
  private
  
  logical,save::namelist_is_read=.false.
  logical,save::lserial=.true.! if .f. run the physics in 'batches' of nhorph
  integer,save::nhorph=500
  integer,save::ntask=1
  integer,public,save::noption(20)=0          ! = option array to select schemes in physics

  type:: phCoefficinets
     real(kind=realkind),allocatable,dimension(:)::emc    !initiated in iniphys/updated in inirad
     real(kind=realkind),allocatable,dimension(:)::csusa  !initiated in iniphys/updated in inirad
     real(kind=realkind),allocatable,dimension(:)::c1er   !initiated in iniphys
     real(kind=realkind),allocatable,dimension(:)::cml2   !initiated in iniphys
     real(kind=realkind),allocatable,dimension(:)::fabso3 !initiated in iniphys/updated in inirad
     real(kind=realkind),allocatable,dimension(:)::cadd   !initiated in iniphys/updated in inirad
     logical::isallocated=.false.
  end type phCoefficinets

  type(phCoefficinets),save::phcoef

  public phcall 
  public read_namprc
contains
  

  subroutine options(koption)  
    !     initialize option array 'noption(20)'
    !     to define which physical parameterizations
    !     will be use d in the physics of hirlam
    !     input parameters:
    !     noption:  number of processes which may be initialized
    !     output parameters:
    !     noption:  integer values identifyiing schemes to
    !               be use d in the physics.
    !     the purpose of initializing the integer array 'noption'
    !     is to define a framework to allow for options
    !     in the physical parameterizations. in reality,
    !     the reference system does not have to support
    !     (many) options. the philosophy is the following:
    !
    !     it is possible to define up to 20 processes in physcis
    !     by means of the option array. noption(j) identifies
    !     process number j. the integer value put into noption(j)
    !     identifies the parameterization scheme to be use d.
    !
    !     in order to implement a given parameterization the
    !     following three steps should be made:
    !        - define integer value noption(j) in subroutine
    !          'options'
    !        - call associated initializing routine in subroutine
    !          'iniphys'
    !        - call the appropriate parameterization from subroutine
    !          'phys' provided that noption(j) has the correct value.
    !     index 'j' is always associated with the same process.
    !     the following definitions have been made:
    !     noption(1) ==> radiation processes
    !     noption(2) ==> condensation processes
    !     noption(3) ==> turbulence processes
    !     noption(4) ==> surface processes
    !     noption(j) , j > 4 have not yet been defined.
    !     regarding the integer values the following principle
    !     applies:
    !     noption(j) = 0 means that process number 'j' will
    !     not be initialized. a warning will be written out in
    !     subroutine 'iniphys' called from 'phcall' in the physics.
    !     it is assumed that a physical parameterization scheme
    !     is always associated with the same prescribed integer
    !     value. at a given instant the default reference system
    !     choice will be defined by one specific array
    !     noption(j)=( k1,k2,k3,k4,...,k20 )
    !     currently, k5=k6=k7=...=k20 =0
    !     predefined options identifying schemes
    !     radiation
    !     noption(1) =1  for savijaervi hirlam radiation scheme
    !     condensation
    !     noption(2) =1  for straco condensation and convection
    !                    scheme
    !
    !     noption(2) =2  for hirlam version of the tiedtke
    !                    mass flux scheme plus the sundqvist
    !                    stratiform condensation scheme
    !
    !     noption(2) =3  for hirlam version of the tiedtke
    !                    mass flux scheme plus the rasch-
    !                    kristjansson condensation scheme
    !
    !     noption(2) =4  for the kain-fritsch convection
    !                    scheme plus the rasch-kristjansson
    !                    condensation scheme
    !     turbulence
    !     noption(3) =1 for the holtslag scheme
    !
    !     noption(3) =2 for the cbr turbulent kinetic
    !                   energy based scheme.
    !     surface processes
    !
    !     noption(4) =1 for the hirlam4.4.3 surface flux
    !                   parameterization ( 'slfluxo' and 'surf')
    !
    !     noption(4) =2 for the isba scheme

    implicit none  
    integer :: koption (20)  
    integer :: j  
    !        radiation
    !        (savijaervi radiation)

    koption (1) = 1  
    !        condensation
    !        (straco scheme)
    koption (2) = 4  
    !
    !        turbulence
    !        (cbr turbulence scheme)
    koption (3) = 2  
    !        surface processes
    !        (hirlam4.3.3 'slfluxo' + 'surf')

    koption (4) = 1  
    !     remaining processes initialized with zero values
    do j = 5, 20  
       koption (j) = 0  
    enddo
  end subroutine options

  subroutine read_namprc()
#ifdef MPI_SRC
    integer::ibuf(22),ierr
#include"mpif.h"
#endif    
    namelist/namprc/nhorph,lserial,noption

    call options(noption)
    if(mype==0)then
       open(57,file='namelists.dat',status='old')
       read(57,nml=namprc)
       close(57)
    endif
#ifdef MPI_SRC
    ibuf(1) = nhorph
    if(lserial)then        !regular input namprc
       ibuf(2)=1
    else 
       ibuf(2)=0
    endif
    ibuf(3:22)=noption
    call mpi_bcast(ibuf,22,mpi_integer,0,localComm,ierr)
    nhorph       = ibuf(1) 
    lserial     = ibuf(2)==1
    noption = ibuf(3:22)
#endif
    if(mype==0)then
       write(6,nml=namprc)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namprc)
       close(1)
    endif

    !     Print out information on the schemes used    
    if(mype==0)then
       write(6,*)' ----------------------------------------'
       write(6,*)' model components:'

       if(nlphys)then
          write(6,*)' physical parametrisations:'
          if( noption(1)==1 )write(6,*)' - savijaervi radiation '
          if( noption(2)==1 )write(6,*)' - straco condensation '
          if( noption(2)==2 )write(6,*)' - tiedtke convection + sundqvist cond.' 
          if( noption(2)==3 )write(6,*)' - tiedtke convection + rasch-kr. cond.' 
          if( noption(2)==4 )write(6,*)' - kain-fr. convection + rasch-kr. cond.'
          if( noption(3)==1 )write(6,*)' - holtslag turbulence scheme'
          if( noption(3)==2 )write(6,*)' - cbr turbulence scheme'
          if( noption(4)==1 )write(6,*)' - hirlam4.4 ''slfluxo+surf'' scheme'
          if( noption(4)==2 )write(6,*)' - isba surface scheme' 
       endif
    endif

  end subroutine read_namprc
  
  subroutine phys(nhor,nlev,kstart,kstop,        &
       kstep, &
       yearc,monthc,dayc,hourc,minc,secc,         &
       dtime,nsvdif,conacc,dtheta,      &
       lsamedt,lnewpht,ldynvd,    &
       t,q,cw,uvel,vvel,omf,        &
       ps,tsea,frice,along,coslat,sinlat,         &
       dpsdin,      &
       cov2d,cwpath,&
       rough,rousea,&
       prcpst,prcpcu,cusnow,stsnow,     &
       drolddt,dsolddt,       &
       accsunny,dtphysh,      &
       orosigm, &
       accrunoff,accrunoffopl,accrunofffor,accrunofflake, & 
                                !    following are prognostic variables stored in arrays svarsm/z/p
                                !    variable description is found below.
                                !    ksvars in gemini should correspond to current number of variables
                                !    in this list.
       tsns2,tsns3,tsns4,tsns5,    &
       tssn2,tssn3,tssn4,tssn5,    &
       tsc2,tsc3,tsc4,tsc5,        &
       tscsn2,tscsn3,tscsn4,tscsn5,&
       snowcan,&
       ticed,ticesn,ticesnd,tsnice,&
       snice,swsnice,rhosnice,snmaxice,      &
       sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,      &
       tsns,tc,tsc,tsnow,tssn,svegopl,svegfor,snopl,snfor,swsn,  &
       rhosn,  &
       snmax,snmaxf,     &
       tice,   &
       tsnc,tscsn,swsnc,rhosnc,    & 
                                !    following are time tendency variables stored in arrays dsvarsdt
                                !    variable description is found below. 
                                !    ksvars in gemini should correspond to current number of variables
                                !    in this list.
       dtsnsdt,dtcdt,dtscdt,dtsndt,dtssndt,dsvegopldt,dsvegfordt, &
       dsnopldt,dsnfordt,dswsndt,drhosndt,dsnmaxdt,dsnmaxfdt,     &
       dticedt, &
       dtsncdt,dtscsndt,dswsncdt,drhosncdt,   &
       dtsns2dt,dtsns3dt,dtsns4dt,dtsns5dt,   &
       dtssn2dt,dtssn3dt,dtssn4dt,dtssn5dt,   &
       dtsc2dt,dtsc3dt,dtsc4dt,dtsc5dt,       &
       dtscsn2dt,dtscsn3dt,dtscsn4dt,dtscsn5dt,         &
       dsnowcandt,        &
       dticeddt,dticesndt,dticesnddt,dtsnicedt,         &
       dsnicedt,dswsnicedt,drhosnicedt,dsnmaxicedt,     &
       dsw1opldt,dsw1fordt,dsw2opldt,dsw2fordt,dsw3opldt,dsw3fordt, & 
                                !    following are diagnostic variables stored in array svar_surf
                                !    variable description is found below.
                                !    msvars in gemini should correspond to current number of variables
                                !    in this list. 
       t2,t2ml,t2mopsn,t2mfor,t2ms,t2mi,  &
       q2,q2ml,q2mopsn,q2mfor,q2ms,q2mi,  &
       u10,&!u10l,&
       u10opsn,&!u10for,&
       u10ms,u10mi,         &
       v10,&!v10l,
       v10opsn,&!v10for,
       v10ms,v10mi,         &
       senf,senfl,senfopsn,senffor,senfs,senfi,     &
       latf,latfl,latfopsn,latffor,latfs,latfi,     &
       momfu,momful,momfuopsn,momfufor,momfus,momfui,         &
       momfv,momfvl,momfvopsn,momfvfor,momfvs,momfvi,         &
       momf,ustar,    &
       frcw,vegopl,frsn,frsnfor,frsnice,  &
       laiopn_int,lai_conif,lai_decid,    &
       albedo,albsnowl,albsnice,albicenl, &
       tseff,tskin,tca,evap,faopet,sn,    &
       tsopsn1,tsopsn2,tsopsn3,tsopsn4,tsopsn5,     &
       tsfor1,tsfor2,tsfor3,tsfor4,tsfor5,&
       swa,lwlai,storage_w,flux_w,        &
       rh2,rh2ml,rh2mopsn,rh2mfor,rh2ms,rh2mi,frsngrid,wevopl,   &
       cloudbot,gustest,gustlow,gustup,tsland1,tsland2, &
       tsland3,tsland4,tsland5,tice1,tice2,mv10,        &
       frdecid,frfor,fr_rain_hrs,emsnowl,     &
       t2mopsnsi,q2mopsnsi,rh2mopsnsi,        &
       u10opsnsi,v10opsnsi,mv10opsnsi,        &
       soilwmm,soilfrwmm, &
       dzsnow,dzsnowopl,dzsnowfor,dzsnowice,  & ! -118
       snowmeltland,  &
       gustestuc,                                & !  120
       lake_types,lake_no_prog,lake_no_tend,lake_no_diag,  &
       frlake, depthlake,   &
       prog_lake,tend_lake,diag_lake, &
       tsclim_years,        & 
                                !    l    tice:        ice surface temperature                         
                                !    l                           
                                !    l tice is updated either outside phys i.e. :                      
                                !    l           1: by values from boundary/climate-data               
                                !    l           2: over the baltic by doesher's model or as 1         
                                !    l                           
                                !    l or it is updated within phys in case of lake-ice, within        
                                !    l the routine surf_sea_ice. 
                                !    l the choice is given by the real array lcounttice(nhor) which    
                                !    l in cases of 0. => dticedt=0.                                    
                                !    l                           
                                !    l note the following shortcoming:                                 
                                !    l in cases where frland < 1 and most of the no-land area is covered                                     
                                !    l by lakes and only little is covered by the baltic, it means that lcounttice                           
                                !    l is 0. and the dticedt=0., i.e. the lake-ice gets the same temperature                                 
                                !    l as that of the baltic. in the worst case, where we have no ice in the baltic, but                     
                                !    l ice in the lakes, we put tice=min(tice,273.15), implying a somewhat erronious tskin,                  
                                !    l and also errors in senf and latf. the area with these errors are in any case, so small                
                                !    l so we neglect it at present.                                    
                                !    l                           
                                !    l    icethick:    thickness of ice from oceanographic model       
                                !    l inout/output:             
                                !    l    albicenl:    albedo of sea/lake-ice (including influence from snow)                                
                                !    l output:                   
                                !    l    dticedt:     time tendency of tice                           
                                !                               
       lcounttice,icethick,      & 
                                !kk  lcounttice used here twohold: as a lake-water/sea-water mask and a
                                !kk  switcher between lake models: 1 - lake is treated a sea           
                                !kk  2 - probe                   
                                !kk  3 - flake                   
                                !kk  as the main destination of lcounttice - to be a mask, we use      
                                !kk  environmental variable flake additionaly                          
                                !    variables for the dynamic vegetation model guess                 
       guessid,    & 
                                !    following is physiographic information from ecoclimap stored in array eco                               
                                !    variable description is found below.                              
                                !    meco in gemini should correspond to current number of variables   
                                !    in this list.              
       lai_t1,lai_t2,lai_t3,z0_t1,z0_t2,     &        ! 1-5
       z0_t3,emis_t1,alb_t1,veg_t1,frland,   &        ! 6-10
       alb_soil,clay,sand,frac_t1,frac_t2,   &        ! 11-15
       frac_t3,alb_t2,alb_t3,emis_t2,emis_t3,&        ! 16-20
       veg_t2,veg_t3,droot_t1,droot_t2,droot_t3,       &        ! 21-25
       dsoil_t1,dsoil_t2,dsoil_t3,rsmin_t1,rsmin_t2,   &        ! 26-30
       rsmin_t3,alb_veg_t1,alb_veg_t2,alb_veg_t3,texture,        &        ! 31-35
       minlai_t1,minlai_t2,minlai_t3,        &        ! 36-38
       maxlai_t1,maxlai_t2,maxlai_t3,        &        ! 39-41
       frac_lake,soil_carb,        &        ! 42-43
       ahyb,bhyb,hybf,hybh,        &
       emc,csusa,c1er,cml2,fabso3,cadd, &
       dtdt,dqdt,dcwdt,dudt,dvdt,  &
       totcov,cucov,dtdtph,dqdtph,dcdtph,    &
       svar,dsvardt,nsvar,         &
       preta, dtp, dtpr, dtpt, dqp, dqpt, dcp, dcpt,   &
       slwnet,sswnet, slwdn, sswdn, tlwnet, tswnet, tswdn,       &
       slwrsea,sswrsea,sswri,slwri, pblh,    &
       omega_kf,dtdt_kf,dqdt_kf,dqcdt_kf,    &
       div_kf,oldcov,zvarcu,cwcu,       &	!added cwcu cgj300611
       raincv_kf,snowcv_kf,umfb,   &
       shal_cgj,zdx_kf,zdy_kf,     &
       nca,kf_ind,       &
       kf_base,kf_top,afull,bfull,cu2d,shal2d)



    use confys  
    use hybridmod
    use radiation
    use co2mod
    use pblhgtmod               
    use cbr                     
    use comdfb                  
    use qnegatmod               
    use flake
    use condensmod              
    use condsmod                
    use config                  
    use surface
    use decomp
    implicit none               
    integer nhor,nlev,kstart,kstop                                    
    integer yearc,monthc,dayc,hourc,minc,secc                         
    integer jk,jl,nsvdif,jstep,l
    logical lprint,lprint7,lprint8,lprint9                      
    logical lsamedt,lnewpht,ldynvd                                    
    real(kind=realkind) dtime,dtvdif,conacc,dtheta                                   
    real(kind=realkind) dtacc,zalat            
    real(kind=realkind),intent(in)::t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),        &
         uvel(nhor,nlev),vvel(nhor,nlev)
    real(kind=realkind)::omf(nhor,nlev),       &
         dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev),         &
         dudt(nhor,nlev),dvdt(nhor,nlev),      &
         totcov(nhor,nlev),cucov(nhor,nlev),   &
         dtdtph(nhor,nlev),dqdtph(nhor,nlev),dcdtph(nhor,nlev)        
    real(kind=realkind),intent(in):: ps(nhor)
    real(kind=realkind)::frice(nhor),       &
         along(nhor),coslat(nhor),sinlat(nhor),&
         rough(nhor),rousea(nhor),   &
         tsea(nhor),       &
         cov2d(nhor),cwpath(nhor),   &
         sw1opl(nhor),sw1for(nhor),sw2opl(nhor),sw2for(nhor),      &
         sw3opl(nhor),sw3for(nhor),  &
         dsw1opldt(nhor),dsw1fordt(nhor),dsw2opldt(nhor),&
         dsw2fordt(nhor),dsw3opldt(nhor),dsw3fordt(nhor),&
         prcpst(nhor),stsnow(nhor),cusnow(nhor),prcpcu(nhor),      &
         dpsdin(nhor)           
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1),hybf(nlev),hybh(nlev+1),        &
         emc(nlev),csusa(nlev),c1er(nlev),cml2(nlev+1),  &
         fabso3(nlev),cadd(nlev) 
!         stddif(nlev),stdif(nlev),stcon(nlev),stdcon(nlev),        &
!         strad(nlev),stdrad(nlev),stclo(nlev),stdclo(nlev),        &
!         stcov(nlev),stdcov(nlev),stccov(nlev),stscov(nlev),       &
!         stcw(nlev),stcpnt(nlev),sttpnt(nlev),stscal(nlev)            
    real(kind=realkind) ts2ns(nhor),ts2sn(nhor),ts2c(nhor),ts2csn(nhor),&
         dts2nsdt(nhor),dts2sndt(nhor),dts2cdt(nhor),dts2csndt(nhor)  
    integer kstep         
    real(kind=realkind) accsunny(nhor),dtphysh 
    real(kind=realkind) orosigm(nhor),&
         !q2d(nhor),    &
         !accprl(nhor),&
         !accprc(nhor),&
         !sswr(nhor), &
         drolddt(nhor),dsolddt(nhor)                                  
    real(kind=realkind) weg(nhor),wevhv(nhor),wesn(nhor)                             
    real(kind=realkind) wveg(nhor),wlai(nhor)  

    real(kind=realkind) laifor_int(nhor)       

    real(kind=realkind) svar(nhor,nlev,*),dsvardt(nhor,nlev,*)                       
    integer nsvar               
    real(kind=realkind) wsv(nhor,nlev,nsvar),dsvdtpr(nhor,nlev,nsvar)                
    real(kind=realkind) preta(nhor,nlev)       
    real(kind=realkind) dtp(nhor,nlev),dtpr(nhor,nlev),dtpt(nhor,nlev), &
         dqp(nhor,nlev),dqpt(nhor,nlev),       &
         dcp(nhor,nlev),dcpt(nhor,nlev)                               
    real(kind=realkind) slwnet(nhor), sswnet(nhor), slwdn(nhor), sswdn(nhor),     &
         tlwnet(nhor), tswnet(nhor), tswdn(nhor)                      
    real(kind=realkind) slwrsea(nhor),sswrsea(nhor),slwri(nhor),sswri(nhor)          
    real(kind=realkind) pblh(nhor)             
    real(kind=realkind) alat(nhor)             
    real(kind=realkind) tice(nhor),wtsns(nhor),wts(nhor)                             

    !    following are diagnostic variables stored in array svar_surf      
    !    msvars in gemini should correspond to current number of variables in this list.                         
    !    varibles are stored in grib-output by specifying the following grib code in run script:                 
    !    par=250, type=105, lev=nn (in list below) for instantaneous output
    !    par=242, type=105, lev=nn (in list below) for mean value over forecast output interval                  
    !                                
    !    variable   description                                    nn      
    !    --------   -----------------------------------------      --      
    !                                
    !    l    t2         grid averaged t2m                              1  
    !    l    t2ml       land averaged t2m                              2  
    !    l    t2mopsn    open land and snow averaged t2m                3  
    !    l    t2mfor     forest (bare soil and snow) t2m                4  
    !    l    t2ms       lake and/or sea water t2m                      5  
    !    l    t2mi       lake and/or sea ice and snow averaged t2m      6  
    !                                
    !    l    q2         grid averaged q2m                              7  
    !    l    q2ml       land averaged q2m                              8  
    !    l    q2mopsn    open land and snow averaged q2m                9  
    !    l    q2mfor     forest (bare soil and snow) q2m                10 
    !    l    q2ms       lake and/or sea water q2m                      11 
    !    l    q2mi       lake and/or sea ice and snow averaged q2m      12 
    !                                
    !    l    u10        grid averaged u10                              13 
    !    l    u10l       land averaged u10                              not defined                              
    !    l    u10opsn    open land and snow averaged u10                15 
    !    l    u10for     forest (bare soil and snow) u10                not defined                              
    !    l    u10ms      lake and/or sea water u10                      17 
    !    l    u10mi      lake and/or sea ice and snow averaged u10      18 
    !                                
    !    l    v10        grid averaged v10                              19 
    !    l    v10l       land averaged v10                              not defined                              
    !    l    v10opsn    open land and snow averaged v10                21 
    !    l    v10for     forest (bare soil and snow) v10                not defined                              
    !    l    v10ms      lake and/or sea water v10                      23 
    !    l    v10mi      lake and/or sea ice and snow averaged v10      24 
    !                                
    !    l    senf       grid averaged senisble heat flux (h)           25 
    !    l    senfl      land averaged h                                26 
    !    l    senfopsn   open land and snow averaged h                  not defined                              
    !    l    senffor    forest (bare soil and snow) h                  not defined                              
    !    l    senfs      lake and/or sea water h                        29 
    !    l    senfi      lake and/or sea ice and snow averaged h        30 
    !                                
    !    l    latf       grid averaged latent heat flux (le)            31 
    !    l    latfl      land averaged le                               32 
    !    l    latfopsn   open land and snow averaged le                 not defined                              
    !    l    latffor    forest (bare soil and snow) le                 not defined                              
    !    l    latfs      lake and/or sea water le                       35 
    !    l    latfi      lake and/or sea ice and snow averaged le       36 
    !                                
    !    l    momfu      grid averaged u momentum flux              37 
    !    l    momful     land averaged                                not defined                              
    !    l    momfuopsn  open land and snow averaged                  not defined                              
    !    l    momfufor   forest (bare soil and snow)                  not defined                              
    !    l    momfus     lake and/or sea water                        41 
    !    l    momfui     lake and/or sea ice and snow averaged        42 
    !                                
    !    l    momfv      grid averaged v momentum flux (vw)             43 
    !    l    momfvl     land averaged vw                               not defined                              
    !    l    momfvopsn  open land and snow averaged vw                 not defined                              
    !    l    momfvfor   forest (bare soil and snow) vw                 not defined                              
    !    l    momfvs     lake and/or sea water vw                       47 
    !    l    momfvi     lake and/or sea ice and snow averaged vw       48 
    !                                
    !    l    momf       grid averaged momentum flux (rho*ustar**2)     49 
    !    l    ustar      grid averaged ustar                            50 
    !                                
    !    l    frcw       fraction forest of fraction land (0-1)         51 
    !    l    vegopl     vegetation cover for open land (0-1)           52 
    !    l    frsn       fraction snow on open land (0-1)               53 
    !    l    frsnfor    fraction snow in forest (0-1)                  54 
    !    l    frsnice    fraction snow on sea ice (0-1)                 55 
    !                                
    !    l    laiopn_int lai for open land vegetation                   56 
    !    l    lai_conif  lai for coniferous forest                      57 
    !    l    lai_decid  lai for deciduous forest                       58 
    !                                
    !    l    albedo     grid averaged albedo                           59 
    !    l    albsnowl   albedo of open land snow                       60 
    !    l    albsnice   albedo of sea ice snow                         61 
    !    l    albicenl   albedo of sea/lake ice                         62 
    !                                
    !    l    tseff      grid averaged surface temperature              63 
    !    l    tskin      grid averaged radiation surface temp           64 
    !    l    tca        forest canopy air temperature                  65 
    !    l    evap       latf in (mm/time step    )                     66 
    !    l    faopet     potential evapotranspiration  (mm /time step)  67 
    !    l    sn         grid averaged snow water eq.                   68 
    !                                
    !    l    tsopsn1    open land/snow averaged 1 soil temp            69 
    !    l    tsopsn2    open land/snow averaged 2 soil temp            70 
    !    l    tsopsn3    open land/snow averaged 3 soil temp            71 
    !    l    tsopsn4    open land/snow averaged 4 soil temp            72 
    !    l    tsopsn5    open land/snow averaged 5 soil temp            73 
    !                                
    !    l    tsfor1     forest bare soil/snow averaged 1 soil temp     74 
    !    l    tsfor2     forest bare soil/snow averaged 2 soil temp     75 
    !    l    tsfor3     forest bare soil/snow averaged 3 soil temp     76 
    !    l    tsfor4     forest bare soil/snow averaged 4 soil temp     77 
    !    l    tsfor5     forest bare soil/snow averaged 5 soil temp     78 
    !                                
    !    l    swa        soil water availability                        79 
    !    l    lwlai      land-averaged lai                              80 
    !    l    storage_w  water storage terms (for water-bal purpose)    81 
    !    l    flux_w     water flux terms (for water-bal purpose)       82 
    !                                
    !    l    rh2        grid averaged relative humidity 2m (rh2m)	83     
    !    l    rh2ml      land averaged rh2m				84     
    !    l    rh2mopsn   open land and snow averaged rh2m		85     
    !    l    rh2mfor    forest (bare soil and snow) rh2m		86     
    !    l    rh2ms      lake and/or sea water rh2m			87     
    !    l    rh2mi      lake and/or sea ice and snow averaged rh2m	88     
    !                                
    !    l    frsngrid   grid averaged snow cover			89     
    !    l    wevopl     total evapotr from open land veg	        90     
    !    l    cloudbot   altitude of cloud bottom			91     
    !    l    gustest    estimated gust wind speed                      92 
    !    l    gustlow    lower bound of gust wind speed                 93 
    !    l    gustup     upper bound of gust wind speed                 94 
    !                                
    !    l    tsland1    land averaged 1 soil temp                      95 
    !    l    tsland2    land averaged 2 soil temp                      96 
    !    l    tsland3    land averaged 3 soil temp                      97 
    !    l    tsland4    land averaged 4 soil temp                      98 
    !    l    tsland5    land averaged 5 soil temp                      99 
    !                                
    !    l    tice1      averaged tice and ticesn temp                  100
    !    l    tice2      averaged ticed and ticesnd temp                101
    !    l    mv10 	 total 10m wind speed (grid)			102    
    !    l    frdecid 	 fraction decidious forest			103                                  
    !    l    frfor 	 fraction forest			        104                                  
    !    l    fr_rain_hrs number of hours with freezing rain            105
    !    l    emsnowl    emissivity of snow over open land              106
    !                                
    !    l    t2mopsnsi  open land, snow, water and ice averaged t2m    107
    !    l    q2mopsnsi  open land, snow, water and ice averaged q2m    108
    !    l    rh2mopsnsi open land, snow, water and ice averaged rh2m	109                                  
    !    l    u10opsnsi  open land, snow, water and ice averaged u10	110                                  
    !    l    v10opsnsi  open land, snow, water and ice averaged v10	111                                  
    !    l    mv10opsnsi open land, snow, water and ice averaged u10	112                                  
    !                                
    !    l    soilwmm    total soil water content (kg/m2 = mm)		113                                  
    !    l    soilfrwmm  total soil frozen water content (kg/m2 = mm)	114                                  
    !                                
    !    l    dzsnow     grid snow depth (m), not w.e.! (wjb)           115
    !    l    dzsnowopl  open land snow depth (m), not w.e.! (wjb)      116
    !    l    dzsnowfor  forest snow depth (m), not w.e.!    (wjb)      117     
    !    l    dzsnowice  sea ice snow depth (m), not w.e.!              118     
    !                                
    !    l    snowmeltland  snow melt rate over land (mm/timestep)	119                                  
    !                                
    !    l    gustestuc  uncorrected estimated gust wind speed          120 

    integer  msvars             
    real(kind=realkind) t2(nhor),t2ml(nhor),t2mopsn(nhor),t2mfor(nhor), &
         t2ms(nhor),t2mi(nhor),      &
         q2(nhor),q2ml(nhor),q2mopsn(nhor),q2mfor(nhor), &
         q2ms(nhor),q2mi(nhor),      &
         u10(nhor),&!u10l(nhor),&
         u10opsn(nhor),&!u10for(nhor),&
         u10ms(nhor),u10mi(nhor),    &
         v10(nhor),&!v10l(nhor),
         v10opsn(nhor),&!v10for(nhor),&
         v10ms(nhor),v10mi(nhor),    &
         senf(nhor),senfl(nhor),senfopsn(nhor),&
         senffor(nhor),senfs(nhor),senfi(nhor),&
         latf(nhor),latfl(nhor),latfopsn(nhor),&
         latffor(nhor),latfs(nhor),latfi(nhor),&
         momfu(nhor),momful(nhor),momfuopsn(nhor),       &
         momfufor(nhor),momfus(nhor),momfui(nhor),       &
         momfv(nhor),momfvl(nhor),momfvopsn(nhor),       &
         momfvfor(nhor),momfvs(nhor),momfvi(nhor),       &
         momf(nhor),ustar(nhor),     &
         frcw(nhor),vegopl(nhor),frsn(nhor),frsnfor(nhor),         &
         frsnice(nhor),    &
         laiopn_int(nhor),lai_conif(nhor),lai_decid(nhor),         &
         albedo(nhor),albsnowl(nhor),albsnice(nhor),albicenl(nhor),&
         tseff(nhor),tskin(nhor),evap(nhor),faopet(nhor),sn(nhor)     
    real(kind=realkind) tsopsn1(nhor),tsopsn2(nhor),tsopsn3(nhor),tsopsn4(nhor),  &
         tsopsn5(nhor),tsfor1(nhor),tsfor2(nhor),tsfor3(nhor),     &
         tsfor4(nhor),tsfor5(nhor),swa(nhor),lwlai(nhor),&
         storage_w(nhor),flux_w(nhor),         &
         rh2(nhor),rh2ml(nhor),rh2mopsn(nhor),rh2mfor(nhor),       &
         rh2ms(nhor),rh2mi(nhor),frsngrid(nhor),wevopl(nhor),      &
         cloudbot(nhor),gustest(nhor),         &
         gustlow(nhor),gustup(nhor),tsland1(nhor),tsland2(nhor),   &
         tsland3(nhor),tsland4(nhor),tsland5(nhor),tice1(nhor),    &
         tice2(nhor),mv10(nhor),     &
         frdecid(nhor),frfor(nhor),fr_rain_hrs(nhor),    &
         emsnowl(nhor),    &
         t2mopsnsi(nhor),q2mopsnsi(nhor),rh2mopsnsi(nhor),         &
         u10opsnsi(nhor),v10opsnsi(nhor),mv10opsnsi(nhor),         &
         soilwmm(nhor),soilfrwmm(nhor),dzsnow(nhor),dzsnowopl(nhor),  &
         dzsnowfor(nhor),dzsnowice(nhor),snowmeltland(nhor),  &
         gustestuc(nhor)
    real(kind=realkind) weightbs(nhor),weightopv(nhor),weightsn(nhor)                     

    !    ecoclimap physiography      
    !                                
    !    variable  description                                     nn      
    !    --------  -----------------------------------------       --      
    !    lai_t1	  lai open land					 1     
    !    lai_t2	  lai conif forest				 2     
    !    lai_t3	  lai broad-leaf forest				 3     
    !    z0_t1	  roughness length open land			 4     
    !    z0_t2	  roughness length conif forest			 5     
    !    z0_t3	  roughness length broad-leaf forest		 6     
    !    emis_t1	  emissivity open land				 7     
    !    alb_t1	  albedo open land				 8     
    !    veg_t1	  vegetation cover open land			 9     
    !    frland	  fraction land					10     
    !    alb_soil	  albedo of bare soil				11     
    !    clay	  percentage of clay				12     
    !    sand	  percentage of sand				13     
    !    frac_t1	  fraction open land				14     
    !    frac_t2	  fraction conif forest				15     
    !    frac_t3	  fraction broad-leaf forest			16     
    !    alb_t2	  albedo conif forest				17     
    !    alb_t3	  albedo broad-leaf forest			18     
    !    emis_t2	  emissivity conif forest			19     
    !    emis_t3	  emissivity broad-leaf forest			20     
    !    veg_t2	  vegetation cover conif forest			21     
    !    veg_t3	  vegetation cover broad-leaf forest		22     
    !    droot_t1	  root depth open land				23     
    !    droot_t2	  root depth conif forest			24     
    !    droot_t3	  root depth broad-leaf forest			25     
    !    dsoil_t1	  soil depth open land				26     
    !    dsoil_t2	  soil depth conif forest			27     
    !    dsoil_t3	  soil depth broad-leaf forest			28     
    !    rsmin_t1	  minimum surface resistance open land		29     
    !    rsmin_t2	  minimum surface resistance conif forest	30     
    !    rsmin_t3 	  minimum surface resistance			31     
    !    alb_veg_t1  albedo vegetation open land			32     
    !    alb_veg_t2  albedo vegetation conif forest		33             
    !    alb_veg_t3  albedo vegetation				34     
    !    texture	  texture according to texture triangle		35     
    !    minlai_t1	  annual min of lai open land			36     
    !    minlai_t2	  annual min of lai conif forest		37     
    !    minlai_t3	  annual min of lai broad-leaf forest		38     
    !    maxlai_t1	  annual max of lai open land			39     
    !    maxlai_t2	  annual max of lai conif forest		40     
    !    maxlai_t3	  annual max of lai broad-leaf forest		41     
    !    frac_lake	  fraction of lake				42     
    !    soil_carb   soil carbon                                   43 				             

    real(kind=realkind)  lai_t1(nhor),lai_t2(nhor),lai_t3(nhor),z0_t1(nhor),      &
         z0_t2(nhor),z0_t3(nhor),emis_t1(nhor),alb_t1(nhor),       &
         veg_t1(nhor),     &
         frland(nhor),alb_soil(nhor),clay(nhor),sand(nhor),        &
         frac_t1(nhor),    &
         frac_t2(nhor),frac_t3(nhor),alb_t2(nhor),alb_t3(nhor),    &
         emis_t2(nhor),emis_t3(nhor),veg_t2(nhor),veg_t3(nhor),    &
         droot_t1(nhor),droot_t2(nhor),droot_t3(nhor),dsoil_t1(nhor),        &
         dsoil_t2(nhor),dsoil_t3(nhor),rsmin_t1(nhor),rsmin_t2(nhor),        &
         rsmin_t3(nhor),alb_veg_t1(nhor),alb_veg_t2(nhor),         &
         alb_veg_t3(nhor),texture(nhor),       &
         minlai_t1(nhor),minlai_t2(nhor),minlai_t3(nhor),&
         maxlai_t1(nhor),maxlai_t2(nhor),maxlai_t3(nhor),&
         frac_lake(nhor),soil_carb(nhor)                                   


    real(kind=realkind) vsw_mix1(nhor),vsw_mix2(nhor),vcc_mix1(nhor),vcc_mix2(nhor),        &
         vfl_mix1(nhor),vfl_mix2(nhor),psis_mix1(nhor),  &
         psis_mix2(nhor),  &
         bw_mix1(nhor),bw_mix2(nhor),ks_mix1(nhor),ks_mix2(nhor),  &
         zcapdry_mix1(nhor),zcapdry_mix2(nhor),&
         zcaps_mix1(nhor),zcaps_mix2(nhor)
    real(kind=realkind) cs_mix1(nhor),cs_mix2(nhor),cs_min,cs_w,cs_i                      
    real(kind=realkind) zemopl(nhor)

    !    === for soil thermal conductivity output                          
    real(kind=realkind) ztlambda1_forns(nhor),ztlambda2_forns(nhor), & !cj thermal cond layer 1-5 forest no snow            
         ztlambda3_forns(nhor),ztlambda4_forns(nhor),    &
         ztlambda5_forns(nhor)  
    real(kind=realkind) ztlambda1_oplns(nhor),ztlambda2_oplns(nhor), & !cj thermal cond layer 1-5 open land no snow         
         ztlambda3_oplns(nhor),ztlambda4_oplns(nhor),    &
         ztlambda5_oplns(nhor)  
    real(kind=realkind) ztlambda1_forsn(nhor),ztlambda2_forsn(nhor), & !cj thermal cond layer 1-5 forest snow               
         ztlambda3_forsn(nhor),ztlambda4_forsn(nhor),    &
         ztlambda5_forsn(nhor)  
    real(kind=realkind) ztlambda1_oplsn(nhor),ztlambda2_oplsn(nhor),  & !cj thermal cond layer 1-5 open land snow            
         ztlambda3_oplsn(nhor),ztlambda4_oplsn(nhor),    &
         ztlambda5_oplsn(nhor)  
    real(kind=realkind) fo1_out(nhor),fo2_out(nhor) !cj                                   

    real(kind=realkind) vegfor(nhor),     &
         soil3wopl(nhor),soil3wfor(nhor),      &
         frroot1wopl(nhor),frroot1wfor(nhor),  &
         frroot2wopl(nhor),frroot2wfor(nhor)                          
    real(kind=realkind) swaopl(nhor),swafor(nhor),swaopl12(nhor),       &
         swafor12(nhor),swaopl3(nhor),swafor3(nhor),     &
         evapl(nhor),evaps(nhor),evapi(nhor),  &
         accrunoffopl(nhor),accrunofffor(nhor),accrunofflake(nhor),&
         accrunoff(nhor)        
    real(kind=realkind) f_wato(nhor),frac_t1o(nhor),&
         frac_t2o(nhor),frac_t3o(nhor),textureo(nhor)                 

    real(kind=realkind) dtdt_kf(nhor,nlev),         &
         dqdt_kf(nhor,nlev),         &
         dqcdt_kf(nhor,nlev),        &
         omega_kf(nhor,nlev),        &
         oldcov(nhor,nlev),&
         rad_cloud(nhor,nlev),       &
         zvarcu(nhor,nlev),&
         wdzf(nhor,nlev),&
!cgj300611
         relh(nhor,nlev), &
         udr(nhor,nlev), &
         udrf(nhor,nlev), &
         rhcrit(nhor,nlev), &
         cwcu(nhor,nlev), &
         div_kf(nhor,nlev),&
         raincv_kf(nhor),  &
         snowcv_kf(nhor),  &
         oldraincv_kf(nhor),         &
         oldsnowcv_kf(nhor),         &
         umfb(nhor),       &
         zdx_kf(nhor),     &
         zdy_kf(nhor),     &
         afull(nlev),bfull(nlev)
    real(kind=realkind) om_cpbl(nhor),    &
         zfrtop(nhor)           
    integer nca(nhor),     &
         kf_ind(nhor),     &
         kf_base(nhor),    &
         kf_top(nhor),     &
         shal_cgj(nhor),   &
         lsb(nhor)              
    integer zfrpbl(nhor),  &
         zlevtop(nhor)          
    real(kind=realkind) cu2d(nhor),shal2d(nhor)
    real(kind=realkind) zwin(nhor)         

    real(kind=realkind) dqdt_cbr(nhor,nlev)    
    real(kind=realkind) dtdt_cbr(nhor,nlev)    
    real(kind=realkind) dcwdt_cbr(nhor,nlev)   
    real(kind=realkind) snowh(nhor)            
    real(kind=realkind) tsns(nhor),tc(nhor),tsc(nhor),tsnow(nhor),tssn(nhor),     &
         svegopl(nhor),svegfor(nhor),snopl(nhor),snfor(nhor),      &
         swsn(nhor),rhosn(nhor),snmax(nhor),snmaxf(nhor)              

    real(kind=realkind) dtsnsdt(nhor),dtcdt(nhor),dtscdt(nhor),dtsndt(nhor),      &
         dtssndt(nhor),dsvegopldt(nhor),dsvegfordt(nhor),&
         dsnopldt(nhor),   &
         dsnfordt(nhor),dswsndt(nhor),drhosndt(nhor),dsnmaxdt(nhor),         &
         dsnmaxfdt(nhor)        
    real(kind=realkind) dticedt(nhor),icethick(nhor)                                 
    real(kind=realkind),intent(in):: lcounttice(nhor)       
    integer guessid(nhor)
    real(kind=realkind) swrad_net_for(nhor),swrad_net_opl(nhor),        &
         cond_for(nhor),cond_opl(nhor),lai_forunder(nhor)             
    real(kind=realkind) dummy1(nhor),dummy2(nhor)   
    real zducld

    !    work space                       

    real(kind=realkind) acprcu,acprst,acsncu,acsnst,sevapr,sevapc                    
    real(kind=realkind) zslask,zfrlim,zcw,zsnw,zfrop,zsnpland,zsnpice,zsnp           
    real(kind=realkind) zfrsea,zsafe,zfrice    
    real(kind=realkind) zlandon,zlandoff,ziceon,ziceoff                              
    real(kind=realkind) ztdumout               
    real(kind=realkind) fracdt                 
    real(kind=realkind) zralfon,zralfoff,zfrtot
    real(kind=realkind) radf(nhor),draindt(nhor),dsnowdt(nhor)                       
    real(kind=realkind) ph(nhor,nlev+1),pf(nhor,nlev),dph(nhor,nlev+1), &
         dpf(nhor,nlev),gpot(nhor,nlev),virt(nhor,nlev)               

    real(kind=realkind) wt(nhor,nlev),wq(nhor,nlev),wu(nhor,nlev),wv(nhor,nlev),  &
         dtdtpr(nhor,nlev),dqdtpr(nhor,nlev),dudtpr(nhor,nlev),    &
         dvdtpr(nhor,nlev),wcw(nhor,nlev),dcdtpr(nhor,nlev)           
    real(kind=realkind) dtdtin(nhor,nlev), dqdtin(nhor,nlev),dcdtin(nhor,nlev),   &
         wsenf(nhor),wlatf(nhor), wmomf(nhor), &
         wmomfu(nhor),wmomfv(nhor)                                    
    real(kind=realkind) wheatv(nhor),wobkl(nhor)                                     
    integer kht(nhor)           
    real(kind=realkind) prsl(nhor,0:nlev),prsi(nhor,0:nlev),  &
         prcl(nhor,0:nlev),prci(nhor,0:nlev)                          
    real(kind=realkind) wustarg(nhor)          
    real(kind=realkind) wevhvfor(nhor),wehvopl(nhor),         &
         wlatfnsbs(nhor),  &
         tca(nhor),        &
         wustarl(nhor),wustari(nhor),wz0land(nhor),      &
         wustaropsn(nhor),wustarfor(nhor),     &
         wdhdtsns(nhor),wdhdtsn(nhor),wdhcdtc(nhor),wdhcdtsc(nhor),&
         wdhscdtsc(nhor),wdhscdtc(nhor),wustars(nhor)                 
    real(kind=realkind)  wfrop(nhor),wfrsnw(nhor),  &
         wfrsn(nhor),wfrsnfor(nhor), &
         wscos(nhor),emskin(nhor),   &
         wevhopl(nhor),    &
         wsenfc(nhor),wlatfc(nhor),wradfc(nhor),         &
         wsenfsn(nhor),wlatfsn(nhor),wradfsn(nhor),      &
         wsenfsc(nhor),wlatfsc(nhor),wradfsc(nhor),      &
         wradfice(nhor),wradsnice(nhor),       &
         wevhvopl(nhor),wradfsnice(nhor),      &
         wsenfns(nhor),wlatfns(nhor),wradfns(nhor),      &
         wetropl1(nhor),wetropl2(nhor),wetropl3(nhor),   &
         wetrfor1(nhor),wetrfor2(nhor),wetrfor3(nhor)                 
    real(kind=realkind) tsnc(nhor),wdhcdtsnc(nhor),wdhsncdtc(nhor),     &
         wdhsncdtsnc(nhor),&
         wsenfsnc(nhor),wlatfsnc(nhor),wradfsnc(nhor),   &
         tscsn(nhor),swsnc(nhor),rhosnc(nhor), &
         dtsncdt(nhor),dtscsndt(nhor),dswsncdt(nhor),    &
         drhosncdt(nhor)        
    real(kind=realkind) wdsenfnsdtns(nhor),wdlatfnsdtns(nhor),&
         wdsenfsndtsn(nhor),wdlatfsndtsn(nhor),&
         wdsenfcdtc(nhor),wdlatfcdtc(nhor),wdsenfcdtsc(nhor),      &
         wdlatfcdtsc(nhor),wdsenfcdtsnc(nhor),wdlatfcdtsnc(nhor),  &
         wdsenfscdtsc(nhor),wdlatfscdtsc(nhor),wdsenfscdtc(nhor),  &
         wdlatfscdtc(nhor),wdsenfsncdtsnc(nhor),wdlatfsncdtsnc(nhor),        &
         wdsenfsncdtc(nhor),wdlatfsncdtc(nhor)                             

    !    snow interception           
    real(kind=realkind) wlatsnowcan(nhor),snowcan(nhor),dsnowcandt(nhor),         &
         wvegvel(nhor),z0wind(nhor)                                   

    real(kind=realkind)  tsns2(nhor),tsns3(nhor), tsns4(nhor),tsns5(nhor),        &
         tssn2(nhor),tssn3(nhor),tssn4(nhor),tssn5(nhor),&
         tsc2(nhor),tsc3(nhor),tsc4(nhor),tsc5(nhor),    &
         tscsn2(nhor),tscsn3(nhor), tscsn4(nhor),tscsn5(nhor),     &
         dtsns2dt(nhor),dtsns3dt(nhor),dtsns4dt(nhor),dtsns5dt(nhor),        &
         dtssn2dt(nhor),dtssn3dt(nhor), dtssn4dt(nhor),dtssn5dt(nhor),       &
         dtsc2dt(nhor),dtsc3dt(nhor), dtsc4dt(nhor),dtsc5dt(nhor), &
         dtscsn2dt(nhor),dtscsn3dt(nhor),      &
         dtscsn4dt(nhor),dtscsn5dt(nhor)                              
    real(kind=realkind) ticed(nhor),ticesn(nhor),ticesnd(nhor),tsnice(nhor),      &
         snice(nhor),swsnice(nhor),rhosnice(nhor),snmaxice(nhor),  &
         dticeddt(nhor),dticesndt(nhor),dticesnddt(nhor),&
         dtsnicedt(nhor),dsnicedt(nhor),dswsnicedt(nhor),&
         drhosnicedt(nhor),dsnmaxicedt(nhor)                               


    !    variables for the lake model flake                                
    integer, intent(in) ::  lake_types, &! number of lake types nb! must be = nlaketype in rcalakeconst!!! , = 3
         lake_no_prog,       & ! number of lake prognostic variables, = 12                                  
         lake_no_tend,       & ! number of lake tendencies = 12       
         lake_no_diag         ! number of lake diagnostic variables (without parameters) = 13               
    real(kind=realkind), dimension(nhor,lake_types), intent(in) :: frlake, depthlake 
    real(kind=realkind), dimension(nhor,lake_types,lake_no_prog), intent(in) ::   &
         prog_lake              
    real(kind=realkind), dimension(nhor,lake_types,lake_no_tend), intent(out) ::  &
         tend_lake              
    real(kind=realkind), dimension(nhor,lake_types,lake_no_diag), intent(inout) ::&
         diag_lake                   

    real(kind=realkind) lakemax1,lakemax8,lakemax9                                   
    real(kind=realkind) tsclim_years(nhor)     
    real(kind=realkind) soilw_surf_GUESS_opl, soilw_deep_GUESS_opl, &
         soilw_surf_GUESS_for, soilw_deep_GUESS_for, snow_GUESS

    !    variable        type            content                           
    !    --------        ------          --------------------------------  
    !                               
    !*integer input  
    !    l                           
    !    l      nhor            input           dimension length in the horizont                                 
    !    l      nlev            input           number of vertical levels  
    !    l      kstart          input           starting index for horizontal lo                                 
    !    l      kstop           input           ending index for hor. loops
    !    l                           
    !    l      yearc,monthc,dayc,   
    !    l      hourc,minc,secc input           current time               
    !    l      nsvar           input           number of extra scalars    
    !    l                           
    !*real input * 
    !    l                           
    !    l      dtime           input           time step for preliminary  
    !    lforecasts of temperatures and                                    
    !    lhumidities (from t, q, dtdt and                                  
    !    ldqdt) to be used in the    
    !    lcondensation calculations  
    !    l      conacc          input           fraction of precipitation at                                     
    !    lthe surface over the time- 
    !    lperiod dtime to be added   
    !    lto accumulated precipitation                                     
    !    l(e.g. =0.5 for leap-frog step)                                   
    !    l                           
    !    ffield accumulation is done in                                    
    !    fphtask                    
    !*full fields input  
    !    l                           
    !    l      t(nhor,nlev)    input           temperatures               
    !    l      q(nhor,nlev)    input           specific humidities        
    !    l                           
    !    l      u(nhor,nlev)         
    !    l      v(nhor,nlev)    input           winds interpolated to t,q-points                                 
    !    l      svar(nhor,nlev,*) input         extra scalars              
    !    l                          
    !*horizontal fields input  
    !    l                           
    !    l      ps(nhor)        input           surface pressure           
    !    l      tsea(nhor)      input           sea surface temperature    
    !    l                           
    !    l      frice(nhor)	input		fraction ice                   
    !    l                           
    !    l      along,coslat,        
    !    l      sinlat(nhor)    input           geometrical data           
    !    l                           
    !    l      dpsdin(nhor)    input           dps/dt                     
    !    l      rough(nhor)     input           roughness length           
    !    l                          
    !*horizontal field transfer  
    !    l                 (transferred between routines in phys)          
    !    l                           
    !    l      radf(nhor)       transfer       radiation at surface       
    !    l      dhdts(nhor)      transfer       sum of derivatives of sensible                                   
    !    land latent heat fluxes with
    !    lrespect to surface temperature                                   
    !    l                          
    !*horizontal fields input,output * 
    !    uh                          
    !    l      accsunny(nhor)  input,output    accumulated sunshine time  
    !    l                          
    !*horizontal field output  
    !    l                       (to be used outside phys)                 
    !    l                           
    !    l     draindt(nhor)    output          rain intensity             
    !    l     dsnowdt(nhor)    output          snowfall intensity         
    !    l                           
    !    f     prcpst(nhor)   output  stratiform precipitation             
    !    f     stsnow(nhor)   output  stratiform snow                      
    !    f     cusnow(nhor)   output convective snow                       
    !    f     prcpcu(nhor)   output convective precipitation (rain and snow)                                    
    !    f     slnet(nhor)    output     sfc longwave net rad.             
    !    f     ssnet(nhor)    output     sfc shortwave net rad.            
    !    f     sldn(nhor)     output     sfc longwave down                 
    !    f     ssdn(nhor)     output     sfc shortgwave down               
    !    f     tlnet(nhor)    output     toa longwave net rad.             
    !    f     tsnet(nhor)    output     toa shortwave net rad.            
    !    f     tsdn(nhor)     output    toa shortwave down                 
    !    f                           
    !    f      pblh(nhor)       output        boundary layer height      
    !*horizontal fields input,output * 
    !    l                           
    !    l      rousea(nhor)    input,output    modified roughness over sea
    !    linitialized as rough(nhor) 
    !    f                          
    !*1-d vertical arrays input  
    !    l                           
    !    l      ahyb(nlev+1)    input           "a" for hybrid coord.      
    !    l      bhyb(nlev+1)    input           "b" for hybrid coord.      
    !    l      hybf(nlev)      input           hybrid cord. on full levels
    !    l      hybh(nlev+1)    input           hybrid cord. on half levels
    !    l                          
    !*1-d vertical arrays input/output * 
    !    l                           
    !    l      emc,csusa      input/output arrays defined in inirad       
    !    l                           
    !    l      c1er(nlev)      input/output    vertical array             
    !    ldefined "first time" in inicond                                  
    !    l      cml2(nlev)      input           vertical array used in vdifnl                                    
    !    l                           
    !99  0613/pr                     
    !    fabso3(nlev),cadd(nlev)        vertical distribution of o3 sw abs 
    !                                
    !    l      stdif,stddif,stcon,  
    !    l      stdcon,strad,stdrad, 
    !    l      stclo,stdclo,stcov,  
    !    l      stdcov,stccov,stscov,
    !    l      stcw,stcpnt,sttpnt,  
    !    l      stscal(nlev)   input,output     arrays used for diagnostics
    !    l                          
    !*full fields input/output * 
    !    l                           
    !    l      dtdt(nhor,nlev) output          total temperature tendency 
    !    l      dqdt(nhor,nlev) output          total tendency of specific humid                                 
    !    l      dsvardt(nhor,nlev) output       total tendency of extra scalars                                  
    !    l                          
    !*full field transfer  
    !    l                           
    !    l      pf(nhor,nlev)   transfer        pressure on "full" levels  
    !    l      ph(nhor,nlev+1) transfer        pressure on "half" levels  
    !    l      dpf(nhor,nlev)  transfer        pressure diff. on "full" levels                                  
    !    l      dph(nhor,nlev+1 transfer        pressure diff. on "full" levels                                  
    !    l      gpot(nlev)      transfer        geopotential               
    !    l      virt(nlev)      transfer        virtual temperature correction                                   
    !    l                          
    !*full fields output * 
    !    l      dudt(nhor,nlev),     
    !    l      dvdt(nhor,nlev) output          wind tendency in t,q-points due                                  
    !    lphysics                    
    !    l                           
    !    l      preta(nhor,nlev)  output    precip. at model levels        
    !    f      dtp(nhor,nlev)    output    tendency of t by phys(ics)     
    !    f      dtpr(nhor,nlev)   output    tendency of t by radia(tion)   
    !    f      dtpt(nhor,nlev)   output    tendency of t by vdiff(usion)  
    !    f      dqp(nhor,nlev)    output    tendency of q by phys          
    !    f      dqpt(nhor,nlev)   output    tendency of q by vdiff(usion)  
    !    f      dcp(nhor,nlev)    output    tendency of c by phys          
    !    f      dcpt(nhor,nlev)   output    tendency of c by vdiff(usion) 
    !*work variables * 
    !    l   wt,wq,wu,wv(nhor,nlev)       t,q,u,and v for temporary use    
    !    l   dtdtpr,dqdtpr,dudtpr,dcdtpr,dtdt,dqdtq,dudt,and dvdt for temporary                                  
    !    l   dvdtpr(nhor,nlev)            use. note that wt is used for temporar                                 
    !    l                                storage of dcwdt!                
    !    f   dtdtin(nhor,nlev)            dtdt on input                    
    !    f   dqdtin(nhor,nlev)            dqdt on input                    
    !    f   dcdtin(nhor,nlev)            dcdt on input                    
    !    f                                the above three arrays are needed
    !    f                                to facilitate the storage of     
    !    f                                of individual physical tendencies.                                     
    !    f   wsenf(nhor)                 incremental surf s.h flux         
    !    f   wlatf(nhor)                 incremental surf l.h. flux        
    !    f   wmomf(nhor)                 incremental surf scalar mom. flux 
    !    f   wmomfu(nhor)                 incremental surf u-mom. flux     
    !    f   wmomfv(nhor)                 incremental surf v-mom. flux     
    !    f   wsenfl(nhor)                 incremental land surf s.h flux   
    !    f   wlatfl(nhor)                 incremental land surf l.h. flux  
    !    l   wts(nhor)                    temporary values for ts          
    !    l   ustarg                       nhor dim. work arrays in holtslag scheme                               
    !    l                           
    !                                
    !      here follows real scalars               
    !                                
    !    acprcu - accumulated precipitation from convective clouds         
    !    acprst - accumulated precipitation from stratiform clouds         
    !    acsncu - accumulated snow from convective clouds                  
    !    acsnst - accumulated snow from stratiform clouds                  
    !    sevapr - accumulated evaporation from precipitation               
    !    sevapc - accumulated evaporation of cloud water                   
    !                                
    !      here follows 1-d real arrays (nhor)       
    !    these are output (cf 199903)
    !    prcpst - stratiform precipitation                                 
    !    stsnow - stratiform snow    
    !    cusnow - convective snow    
    !    prcpcu - convective precipitation (rain and snow)                 
    !                               
    !*         
    !                                
    !    new variables over land:    
    !                                
    !    l    tsns:      no snow temperature                               
    !    l    tc:        canopy temperature                                
    !    l    tsc:       ground temperature below trees                    
    !    l    tsnow:     snow temperature over land                        
    !    l    tssn:      soil temperature under snow                       
    !    l    svegopl:   water on open land canopy (m)                     
    !    l    svegfor:   water on forest canopy (m)                        
    !    l    snopl:     land average snow depth over land in water equivalent (m)                               
    !    l    snfor:     land average snow depth in forest in water equivalent (m)                               
    !    l    swsn:      actual water amount in snow over land             
    !    l    rhosn:     current density of snow over land                 
    !    l    albsnowl:  current albedo of snow over land                  
    !    l    snmax:     snowmax over open land (hydrological treatment)   
    !    l    snmaxf:    snowmax over forest land (hydrological treatment) 
    !    l output:                   
    !    l    dtsnsdt:   time tendency of tsns                             
    !    l    dtcdt:     time tendency of tc                               
    !    l    dtscdt:    time tendency of tsc                              
    !    l    dtsndt:    time tendency of tsnow                            
    !    l    dtssndt:   time tendency of tssn                             
    !    l    dsvegopldt:time tendency of svegopl                          
    !    l    dsvegfordt:time tendency of svegfor                          
    !    l    dsnopldt:  time tendency of snopl                            
    !    l    dsnfordt:  time tendency of snfor                            
    !    l    dswsndt:   time tendency of swsn                             
    !    l    drhosndt:  time tendency of rhosn                            
    !    l    dsnmaxdt:  time tendency of snmax                            
    !    l    dsnmaxfdt:  time tendency of snmaxf                          
    !                               
    !*         
    !                                
    logical locw              ! cloud water scheme                    
    logical lesat             ! controls calculation of esat in partially cloudy boxes                      
    integer  iut                


    dudtpr = 666.0_realkind
    !                                
    locw = ( noption(2)>0 )       

    !    gj300408 initialise logical lesat just once and pass into respectively                                  
    !    gj300408 cldfrc.f findsp.f pcond.f vcbr.f and vqsatd.f            
    lesat=.false.               
    zfrlim=0.01_realkind                 
    ztdumout=99._realkind                

    do jk=1,nlev                
       do jl=kstart,kstop       
          dtdt_cbr(jl,jk)=0._realkind    
          dcwdt_cbr(jl,jk)=0._realkind   
          dqdt_cbr(jl,jk)=0._realkind    
          dtdtin(jl,jk) = dtdt(jl,jk)                                 
          dqdtin(jl,jk) = dqdt(jl,jk)                                 
          dcdtin(jl,jk) = dcwdt(jl,jk)                                     

          !    gj160600  get total cloud cover for use in radia                       

          rad_cloud(jl,jk)=max(0._realkind,min(1.0_realkind,     &
               (cucov(jl,jk)+totcov(jl,jk)) ))                        
       enddo
    enddo


    !    compute gpot, pressures and pressure diferences                   
    !    using 'ahyb' and 'bhyb' in subroutine 'ahybrid'                        

    call ahybrid (nhor,nlev,kstart,kstop,      &
         t,q,    &
         ps,     &
         gpot,ph,pf,dph,dpf,         &
         ahyb,bhyb,hybf)             

    !     cacculate lowest cloud bottom height (meters) = lowest layer  with cloud fraction>=50%         
    do jl=kstart,kstop          
       cloudbot(jl)=0._realkind         
    enddo
    do jk=1,nlev                
       do jl=kstart,kstop       
          if (rad_cloud(jl,jk) > 0.5_realkind)       &
               cloudbot(jl)=gpot(jl,jk)/gravit                        
       enddo
    enddo

    !    compute physiographic information needed for the land surface subroutines                               
    !    based on ecoclimap input    
!RCA4

    call conv_ecocli(nhor,kstart,kstop,        &
         lai_t1,lai_t2,lai_t3,&
         veg_t1,emis_t1,frland,   &
         frac_t2,   &
         frac_t3,&
         veg_t2,veg_t3,droot_t1,droot_t2,droot_t3,       &
         texture,        &
         minlai_t1,minlai_t2,minlai_t3,        &
         maxlai_t1,maxlai_t2,maxlai_t3,        &
         tsns4,tsc4, &
         sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,      &
         soil_carb,        &
         along,coslat,     &
         vsw_mix1,vsw_mix2,vcc_mix1,vcc_mix2,  &
         vfl_mix1,vfl_mix2,psis_mix1,psis_mix2,&
         bw_mix1,bw_mix2,ks_mix1,ks_mix2,      &
         zcapdry_mix1,zcapdry_mix2,zcaps_mix1,zcaps_mix2,&
         cs_mix1,cs_mix2,cs_min,cs_w,cs_i,     &
         frfor,frdecid,vegopl,       &
         laiopn_int,lai_conif,lai_decid,       &
         vegfor,soil3wopl,soil3wfor, &
         frroot1wopl,frroot1wfor,    &
         frroot2wopl,frroot2wfor,    &
         tsclim_years,zemopl)
!RCA35
!      call conv_ecocli(nhor,kstart,kstop,  &
!          lai_t1,lai_t2,lai_t3,z0_t1,z0_t2,  &
!          z0_t3,emis_t1,alb_t1,veg_t1,frland,  &
!          alb_soil,clay,sand,frac_t1,frac_t2,  &
!          frac_t3,alb_t2,alb_t3,emis_t2,emis_t3,  &
!          veg_t2,veg_t3,droot_t1,droot_t2,droot_t3,  &
!          dsoil_t1,dsoil_t2,dsoil_t3,rsmin_t1,rsmin_t2,  &
!          rsmin_t3,alb_veg_t1,alb_veg_t2,alb_veg_t3,texture,  &
!          minlai_t1,minlai_t2,minlai_t3,  &
!          maxlai_t1,maxlai_t2,maxlai_t3,  &
!          frac_lake,  &
!          frice,tsns4,tsc4,  &
!          sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,  &
!          along,coslat,  &
!          frfor,frdecid,vegopl,  &
!          laiopn_int,lai_conif,lai_decid,  &
!          vegfor,soil3wopl,soil3wfor,  &
!          frroot1wopl,frroot1wfor,  &
!          frroot2wopl,frroot2wfor)



    !    compute effective temperatures used in the radiation scheme (tskin)                                     
    !    and fractions etc:          
!RCA4
    call calctemps(nhor,kstart,kstop,&
         monthc,dayc,hourc,minc,secc,&
         along,coslat,sinlat,        &
         tsns,tc,tsc,tsnow,tsea,tice,tsnice,   &
         snopl,snfor,snice,       &
         snmax,snmaxf,snmaxice,      &
         frice,frland,frfor,         &
         orosigm,lcounttice,albsnowl,&
         albsnice,albicenl,&
         lai_conif,lai_decid,vegopl,cov2d,     &
         frdecid,&
         alb_t1,alb_t2,alb_t3,       &
         zemopl,emis_t2,emis_t3,    &
         texture,&
         tsnc,   &
         prog_lake(:,:,1), frlake,   &
         prog_lake(:,:,8) , prog_lake(:,:,9),  &
         wfrop,frcw,wfrsnw,&
         frsn,frsnfor,frsnice,       &
         albedo,wscos,tskin,emskin)                                   

!RCA35
!      call calctemps(nhor,kstart,kstop,  &
!                monthc,dayc,hourc,minc,secc,  &
!                along,coslat,sinlat,  &
!                tsns,tc,tsc,tsnow,tsea,tice,tsnice,  &
!                ps,snopl,snfor,snice,  &
!                snmax,snmaxf,snmaxice,  &
!                frice,frland,frfor,  &
!                orosigm,lcounttice,albsnowl,  &
!                albsnice,albicenl,  &
!                lai_conif,lai_decid,vegopl,cov2d,  &
!                frdecid,  &
!                alb_t1,alb_t2,alb_t3,  &
!                emis_t1,emis_t2,emis_t3,  &
!                texture,  &
!                tsnc,  &
!                prog_lake(:,:,1), frLake,   &
!                prog_lake(:,:,8) , prog_lake(:,:,9),  &
!                wfrop,frcw,wfrsnw,  &
!                frsn,frsnfor,frsnice,  &
!                albedo,wscos,tskin,emskin)



    do jl=kstart,kstop          
       f_wato(jl)=frland(jl)    
       frac_t1o(jl)=emskin(jl)  
       frac_t2o(jl)=frcw(jl)    
       frac_t3o(jl)=frfor(jl)   
       oldraincv_kf(jl)=raincv_kf(jl)                                 
       oldsnowcv_kf(jl)=snowcv_kf(jl)                                 

    enddo

    !    radiation ( parameter noption(1) )                                     

    if( noption(1)/=0 ) then  
       !    gj   slwdn(jl) is used before it is calculated in radia.f...ie    
       !    gj   previous timestep value is used in correct slwup(sfc) calc.  
       !    gj   thus at timestep 1, slwdn is initialised in gemini to 200wm-2
       do jl=kstart,kstop       
          !                                
          slwnet(jl)=0._realkind         
          sswnet(jl)=0._realkind         
          sswdn(jl)=0._realkind          
          tswnet(jl)=0._realkind         
          tswdn(jl)=0._realkind          
          tlwnet(jl)=0._realkind         
          snowh(jl)=snopl(jl)+snfor(jl)                               

          !    compute grid average snow cover fraction                          
          !    frsngrid=1 means that all land+ice is snow covered.               

          zsafe=0._realkind              
          if((frland(jl)+(1._realkind-frland(jl))*frice(jl))<=0._realkind)zsafe=1._realkind    
          frsngrid(jl)=( frland(jl)*(frcw(jl)*frsnfor(jl)+         &
               (1._realkind-frcw(jl))*frsn(jl))+        &
               (1._realkind-frland(jl))*frice(jl)*frsnice(jl) )   &
               /(frland(jl)+(1._realkind-frland(jl))*frice(jl)+zsafe)          

       enddo

       !    radiation scheme                 


       call aradia( nhor, nlev, kstart, kstop, &
            hourc, minc, secc,  locw, hybf, hybh,  &
            emc, csusa, fabso3, cadd, along, coslat, sinlat,       &
            albedo, tskin, t, q, cw, pf, ph, frland, emskin,       &
            dpf, rad_cloud, dtdtpr, &
            accsunny  , dtphysh  ,   &
            slwdn, sswdn, slwnet, sswnet, tlwnet,tswnet,tswdn,     &
            radf,frice,ps,frsngrid,  &
!cgj300611
            wfrop,vegopl,  &
            texture,gpot,pblh,totcov,cucov,cwcu,zwin)
!cgj300611    

       do jk=1,nlev             
          do jl=kstart,kstop    
             dtdt(jl,jk)=dtdt(jl,jk)+dtdtpr(jl,jk)                    
             dtpr(jl,jk) = dtdtpr(jl,jk)                              
          enddo
       enddo

       !    put variables to zero which are not used                          
       !    slwrs,sswrs will be overwritten in slfluxo_sea_ice                

       do jl=kstart,kstop       

          dticedt(jl)=0._realkind        
          slwrsea(jl)=0._realkind        
          sswrsea(jl)=0._realkind        
          swaopl(jl)=0._realkind         
          swaopl12(jl)=0._realkind       
          swaopl3(jl)=0._realkind        
          swafor(jl)=0._realkind         
          swafor12(jl)=0._realkind       
          swafor3(jl)=0._realkind        
          swrad_net_for(jl)=0._realkind  
          swrad_net_opl(jl)=0._realkind  
          cond_for(jl)=0._realkind       
          cond_opl(jl)=0._realkind       
          lai_forunder(jl)=0._realkind   
          dummy1(jl)=0._realkind         
          dummy2(jl)=0._realkind         
       enddo

    endif

    !    now compute tendencies due to vertcal difusion                    
    !    ( noption(3) different from zero ) of                             
    !    forecast variables in (t,q)-points                                
    if( noption(3)/=0 ) then  
       do l=1,nsvar             
          do jk=1,nlev          
             do jl=kstart,kstop 
                wsv(jl,jk,l) = svar(jl,jk,l)                          
             enddo
          enddo
       enddo
       do jk=1,nlev             
          do jl=kstart,kstop    
             wt(jl,jk) = t(jl,jk)                                     
             wq(jl,jk) = q(jl,jk)                                     
             wu(jl,jk) = uvel(jl,jk)                                     
             wv(jl,jk) = vvel(jl,jk)                                     
             wcw(jl,jk) = cw(jl,jk)                                   
          enddo
       enddo

       do jl=kstart,kstop       
          wtsns(jl) = tsns(jl)  
       enddo

       dtvdif = dtime/real(nsvdif,realkind)                                   
       fracdt = dtvdif/dtime         

       do 90 jstep=1,nsvdif     
          if(ldynvd) then       
             do l=1,nsvar       
                do jk=1,nlev    
                   do jl=kstart,kstop                                 
                      wsv(jl,jk,l)=wsv(jl,jk,l)+         &
                           dtvdif*dsvardt(jl,jk,l)                    
                   enddo
                enddo
             enddo
             
             
             do jk=1,nlev       
                do jl=kstart,kstop                                    
                   wt(jl,jk) = wt(jl,jk) + dtvdif*dtdt(jl,jk)         
                   wq(jl,jk) = wq(jl,jk) + dtvdif*dqdt(jl,jk)         
                   wu(jl,jk) = wu(jl,jk) + dtvdif*dudt(jl,jk)         
                   wv(jl,jk) = wv(jl,jk) + dtvdif*dvdt(jl,jk)         
                   wcw(jl,jk) = wcw(jl,jk) + dtvdif*dcwdt(jl,jk)      
                enddo
             enddo

          endif

          !    surface flux computations   
          !    ( noption(4)=1 )            

          if( noption(4)==1 ) then                                  
             !    call slfluxo for land and sea_ice, respectively, followed by      
             !    slfluxo_average             
!RCA4
             call slfluxo_land(nhor,nlev,kstart,kstop,   &
                  dtvdif,  &
                  wt,wq,wu,wv,wsv(1,1,1),gpot,dph,       &
                  ps,wtsns,tsc,tc,tsnow,tskin, &
                  sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,       &
                  svegfor,svegopl,   &
                  frland,texture,      &
                  drolddt,dsolddt,   &
                  oldraincv_kf,oldsnowcv_kf,   &
                  radf,albedo,wscos,sswdn,emskin,        &
                  wfrop,frcw,wfrsnw,frsnfor,albsnowl,    &
                  tsnc,tscsn,        &
                  snowcan, &
                  tsns2,tsns3,tsns4,tsns5,     &
                  tssn2,tssn3,tssn4,tssn5,     &
                  tsc2,tsc3,tsc4,tsc5,         &
                  tscsn2,tscsn3,tscsn4,tscsn5, &
                  laiopn_int,lai_conif,lai_decid,        &
                  vegopl,cov2d,      &
                  frdecid,vegfor,    &
                  soil3wopl,soil3wfor,         &
                  frroot1wopl,frroot1wfor,     &
                  frroot2wopl,frroot2wfor,     &
                  alb_t1,alb_t2,alb_t3,zemopl,emis_t2,emis_t3,    &
                  z0_t1,z0_t2,z0_t3,rsmin_t1,rsmin_t2,rsmin_t3,    &
                  along,coslat,sinlat,kstep,   &
                  vcc_mix1,vcc_mix2,vfl_mix1,vfl_mix2,   &
                  senffor,latffor,tca,q2mfor,  &
                  wetropl1,wetropl2,wetropl3,wetrfor1,wetrfor2,    &
                  wetrfor3,&
                  wevhvfor,wevhvopl,wlatfnsbs, &
                  wsenfns,wlatfns,wradfns,     &
                  wsenfc,wlatfc,wradfc,        &
                  wsenfsn,wlatfsn,wradfsn,     &
                  wsenfsc,wlatfsc,wradfsc,     &
                  wdhdtsns,wdhdtsn,wdhcdtc,wdhcdtsc,     &
                  wdhscdtsc,wdhscdtc,&
                  tseff,t2ml,q2ml,t2mopsn,q2mopsn,t2mfor,&
                  rh2ml,rh2mopsn,rh2mfor,      &
                  u10opsn,v10opsn,   &
                  senfl,latfl,evapl,wustarl,wustaropsn,wustarfor,  &
                  wz0land, &
                  wsenfsnc,wlatfsnc,wradfsnc,  &
                  wdhcdtsnc,wdhsncdtc,wdhsncdtsnc,       &
                  wdsenfnsdtns,wdlatfnsdtns,   &
                  wdsenfsndtsn,wdlatfsndtsn,   &
                  wdsenfcdtc,wdlatfcdtc,wdsenfcdtsc,wdlatfcdtsc,   &
                  wdsenfcdtsnc,wdlatfcdtsnc,   &
                  wdsenfscdtsc,wdlatfscdtsc,wdsenfscdtc,wdlatfscdtc,         &
                  wdsenfsncdtsnc,wdlatfsncdtsnc,wdsenfsncdtc,      &
                  wdlatfsncdtc,      &
                  swrad_net_opl,swrad_net_for,cond_opl,cond_for,   &
                  wlatsnowcan,wvegvel,wevopl,  &
                  weightbs,weightopv,weightsn,faopet)                             
!RCA35
!      call slfluxo_land(nhor,nlev,kstart,kstop,  &
!      monthc,dayc,  &
!      dtvdif,  &
!      wt,wq,wu,wv,gpot,dph,  &
!      ps,wtsns,tsc,tc,tsnow,tskin,  &
!      sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,  &
!      svegfor,svegopl,  &
!      frland,orosigm,texture,  &
!      drolddt,dsolddt,  &
!      oldRAINCV_KF,oldSNOWCV_KF,  &
!      radf,albedo,wscos,sswdn,emskin,  &
!      wfrop,frcw,wfrsnw,frsnfor,albsnowl,  &
!      tsnc,tscsn,  &
!      snowcan,  &
!      tsns2,tsns3,tsns4,tsns5,  &
!      tssn2,tssn3,tssn4,tssn5,  &
!      tsc2,tsc3,tsc4,tsc5,  &
!      tscsn2,tscsn3,tscsn4,tscsn5,  &
!      laiopn_int,lai_conif,lai_decid,  &
!      vegopl,cov2d,  &
!      frdecid,vegfor,  &
!      soil3wopl,soil3wfor,  &
!      frroot1wopl,frroot1wfor,  &
!      frroot2wopl,frroot2wfor,  &
!      alb_t1,alb_t2,alb_t3,emis_t1,emis_t2,emis_t3,  &
!      z0_t1,z0_t2,z0_t3,rsmin_t1,rsmin_t2,rsmin_t3,  &
!      along,coslat,sinlat,kstep,  &
!      senffor,latffor,tca,q2mfor,  &
!      wetropl1,wetropl2,wetropl3,wetrfor1,wetrfor2,wetrfor3,  &
!      wevhvfor,wevhvopl,wlatfnsbs,  &
!      wsenfns,wlatfns,wradfns,  &
!      wsenfc,wlatfc,wradfc,  &
!      wsenfsn,wlatfsn,wradfsn,  &
!      wsenfsc,wlatfsc,wradfsc,  &
!      wdhdtsns,wdhdtsn,wdhcdtc,wdhcdtsc,  &
!      wdhscdtsc,wdhscdtc,  &
!      tseff,t2ml,q2ml,t2mopsn,q2mopsn,t2mfor,  &
!      rh2ml,rh2mopsn,rh2mfor,  &
!      u10opsn,v10opsn,  &
!      senfl,latfl,evapl,wustarl,wustaropsn,wustarfor,wz0land,  &
!      wsenfsnc,wlatfsnc,wradfsnc,  &
!      wdhcdtsnc,wdhsncdtc,wdhsncdtsnc,  &
!      wdsenfnsdtns,wdlatfnsdtns,  &
!      wdsenfsndtsn,wdlatfsndtsn,  &
!      wdsenfcdtc,wdlatfcdtc,wdsenfcdtsc,wdlatfcdtsc,  &
!      wdsenfcdtsnc,wdlatfcdtsnc,  &
!      wdsenfscdtsc,wdlatfscdtsc,wdsenfscdtc,wdlatfscdtc,  &
!      wdsenfsncdtsnc,wdlatfsncdtsnc,wdsenfsncdtc,wdlatfsncdtc,  &
!      wlatsnowcan,wvegvel,wevopl)


             !    note that tseff from slfluxo_land is average temperature over land
             !    will be updated to grid average in slfluxo_average                

!RCA4
             call surf_land(nhor,kstart,kstop, &
                  dtime,monthc,dayc, &
                  tsns,    &
                  sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,       &
                  snopl,snfor,svegfor,svegopl, &
                  snmax,snmaxf,      &
                  wradfns,wsenfns,wlatfns,wdhdtsns,      &
                  frice,frland,      &
                  wfrop,frcw,wfrsnw,frsnfor,        &
                  drolddt,dsolddt,   &
                  oldraincv_kf,oldsnowcv_kf,   &
                  tc,wdhcdtc,wdhcdtsc,         &
                  wsenfc,wlatfc,wradfc,        &
                  tsnow,wsenfsn,wlatfsn,wradfsn,wdhdtsn, &
                  wetropl1,wetropl2,wetropl3,wetrfor1,wetrfor2,    &
                  wetrfor3,&
                  wevhvopl,wevhvfor,orosigm,   &
                  conacc,texture,    &
                  tssn,rhosn,swsn,albsnowl,    &
                  tsc,wdhscdtc,wdhscdtsc,      &
                  wsenfsc,wlatfsc,wradfsc,wlatfnsbs,     &
                  tsnc,wdhcdtsnc,wdhsncdtc,wdhsncdtsnc,  &
                  wsenfsnc,wlatfsnc,wradfsnc,  &
                  tscsn,swsnc,rhosnc,&
                  vsw_mix1,vsw_mix2,vcc_mix1,vcc_mix2,   &
                  vfl_mix1,vfl_mix2,psis_mix1,psis_mix2, &
                  bw_mix1,bw_mix2,ks_mix1,ks_mix2,       &
                  zcapdry_mix1,zcapdry_mix2,zcaps_mix1,zcaps_mix2, &
                  cs_mix1,cs_mix2,cs_min,cs_w,cs_i,      &
                  ztlambda1_forns,ztlambda2_forns, & !cj thermal cond layer 1-5 forest no snow               
                  ztlambda3_forns,ztlambda4_forns,       &
                  ztlambda5_forns,   &
                  ztlambda1_oplns,ztlambda2_oplns, & !cj thermal cond layer 1-5 open land no snow            
                  ztlambda3_oplns,ztlambda4_oplns,       &
                  ztlambda5_oplns,   &
                  ztlambda1_forsn,ztlambda2_forsn, & !cj thermal cond layer 1-5 forest snow                  
                  ztlambda3_forsn,ztlambda4_forsn,       &
                  ztlambda5_forsn,   &
                  ztlambda1_oplsn,ztlambda2_oplsn, & !cj thermal cond layer 1-5 open land snow               
                  ztlambda3_oplsn,ztlambda4_oplsn,       &
                  ztlambda5_oplsn,   &
                  wlatsnowcan,snowcan,wvegvel, &
                  laiopn_int,lai_conif,lai_decid,        &
                  vegopl,frdecid,vegfor,       &
                  soil3wopl,soil3wfor,         &
                  zemopl,emis_t2,emis_t3,     &
                  along,coslat,sinlat,kstep,   &
                  accrunoffopl,accrunofffor,   &
                  dtsnsdt, &
                  dsw1opldt,dsw1fordt,dsw2opldt,dsw2fordt,dsw3opldt,         &
                  dsw3fordt,         &
                  dsnopldt,dsnfordt, &
                  dsvegfordt,dsvegopldt,dsnmaxdt,dsnmaxfdt,        &
                  dtcdt,dtscdt,dtsndt,         &
                  drhosndt,dtssndt,dswsndt,    &
                  dtsncdt,dtscsndt,dswsncdt,drhosncdt,   &
                  dsnowcandt,        &
                  tsns2,tsns3,tsns4,tsns5,     &
                  tssn2,tssn3,tssn4,tssn5,     &
                  tsc2,tsc3,tsc4,tsc5,         &
                  tscsn2,tscsn3,tscsn4,tscsn5, &
                  dtsns2dt,dtsns3dt,dtsns4dt,dtsns5dt,   &
                  dtssn2dt,dtssn3dt,dtssn4dt,dtssn5dt,   &
                  dtsc2dt,dtsc3dt,dtsc4dt,dtsc5dt,       &
                  dtscsn2dt,dtscsn3dt,dtscsn4dt,dtscsn5dt,         &
                  swaopl,swafor,swa,swaopl12,swafor12,   &
                  swaopl3,swafor3,soilwmm,soilfrwmm,     &
                  lwlai,storage_w,flux_w,fr_rain_hrs,    &
                  emsnowl, snowmeltland, &
                  dzsnowopl,dzsnowfor)                                
!RCA35
!      call surf_land(nhor,kstart,kstop,  &
!      dtime,monthc,dayc,  &
!      tsns,  &
!      sw1opl,sw1for,sw2opl,sw2for,sw3opl,sw3for,  &
!      snopl,snfor,svegfor,svegopl,  &
!      snmax,snmaxf,  &
!      wradfns,wsenfns,wlatfns,wdhdtsns,  &
!      frice,frland,  &
!      wfrop,frcw,wfrsnw,frsn,frsnfor,  &
!      drolddt,dsolddt,  &
!      oldRAINCV_KF,oldSNOWCV_KF,  &
!      tc,wdhcdtc,wdhcdtsc,  &
!      wsenfc,wlatfc,wradfc,  &
!      tsnow,wsenfsn,wlatfsn,wradfsn,wdhdtsn,  &
!      wetropl1,wetropl2,wetropl3,wetrfor1,wetrfor2,wetrfor3,  &
!      wevhvopl,wevhvfor,orosigm,  &
!      conacc,texture,  &
!      tssn,rhosn,swsn,albsnowl,  &
!      tsc,wdhscdtc,wdhscdtsc,  &
!      wsenfsc,wlatfsc,wradfsc,wlatfnsbs,  &
!      tsnc,wdhcdtsnc,wdhsncdtc,wdhsncdtsnc,  &
!      wsenfsnc,wlatfsnc,wradfsnc,  &
!      tscsn,swsnc,rhosnc,  &
!      wlatsnowcan,snowcan,wvegvel,  &
!      laiopn_int,lai_conif,lai_decid,  &
!      vegopl,frdecid,vegfor,  &
!      soil3wopl,soil3wfor,  &
!      emis_t1,emis_t2,emis_t3,  &
!      along,coslat,sinlat,kstep,  &
!      accrunoffopl,accrunofffor,  &
!      dtsnsdt,  &
!      dsw1opldt,dsw1fordt,dsw2opldt,dsw2fordt,dsw3opldt,dsw3fordt,  &
!      dsnopldt,dsnfordt,  &
!      dsvegfordt,dsvegopldt,dsnmaxdt,dsnmaxfdt,  &
!      dtcdt,dtscdt,dtsndt,  &
!      drhosndt,dtssndt,dswsndt,  &
!      dtsncdt,dtscsndt,dswsncdt,drhosncdt,  &
!      dsnowcandt,  &
!      tsns2,tsns3,tsns4,tsns5,  &
!      tssn2,tssn3,tssn4,tssn5,  &
!      tsc2,tsc3,tsc4,tsc5,  &
!      tscsn2,tscsn3,tscsn4,tscsn5,  &
!      dtsns2dt,dtsns3dt,dtsns4dt,dtsns5dt,  &
!      dtssn2dt,dtssn3dt,dtssn4dt,dtssn5dt,  &
!      dtsc2dt,dtsc3dt,dtsc4dt,dtsc5dt,  &
!      dtscsn2dt,dtscsn3dt,dtscsn4dt,dtscsn5dt,  &
!      swaopl,swafor,swa,swaopl12,swafor12,  &
!      swaopl3,swafor3,soilwmm,  &
!      lwlai,storage_w,flux_w,fr_rain_hrs,  &
!      emsnowl)


!ps       call mod_flux_for_implicit_solution(nhor,kstart,kstop,   &
!ps            dtime,             &
!ps! ---- input ----         
!ps!      alphaimp,         
!ps            wfrop,frcw,wfrsnw,frsnfor,                               &
!ps!      time tendencies    
!ps            dtsnsdt,dtcdt,dtscdt,dtsndt,dtsncdt,                     &
!ps!      flux derivatives   
!ps            wdsenfnsdtns,wdlatfnsdtns,                               &
!ps            wdsenfsndtsn,wdlatfsndtsn,                               &
!ps            wdsenfcdtc,wdlatfcdtc,wdsenfcdtsc,wdlatfcdtsc,           &
!ps            wdsenfcdtsnc,wdlatfcdtsnc,                               &
!ps            wdsenfscdtsc,wdlatfscdtsc,wdsenfscdtc,wdlatfscdtc,       &
!ps            wdsenfsncdtsnc,wdlatfsncdtsnc,wdsenfsncdtc,wdlatfsncdtc, &
!ps! ---- input/output ----  
!ps!      fluxes             
!ps            wsenfns,wlatfns,wsenfsn,wlatfsn,wsenfc,wlatfc,           &
!ps            wsenfsc,wlatfsc,wsenfsnc,wlatfsnc,                       &
!ps            senffor,latffor,senfl,latfl                              &
!ps            )                       

             if(use_guess)then  
                !    lai calculated by guess     
                do jl=kstart,kstop                                    
                   alat(jl)=acos(coslat(jl))*180._realkind/3.141592654_realkind         
                   if(sinlat(jl)<0._realkind)alat(jl)=-alat(jl)             
                   lprint7=.false.                                    
                   if(along(jl)>16.9_realkind .and. along(jl)<16.2_realkind )then 
                      if(alat(jl)>59.7_realkind .and. alat(jl)< 59.8_realkind)then
                         lprint7=.true.                               
                      endif
                   endif

                   frfor(jl)=0.5_realkind
                   frdecid(jl)=0.5_realkind                                    
                   vegopl(jl)=0.5_realkind                                     
                   laiopn_int(jl)=2._realkind                                  
                   lai_conif(jl)=2._realkind                                   
                   lai_decid(jl)=2._realkind     

                 ! Weight swa to fit GUESS upper (0-0.5m) and lower (0.5-1.5m) soil layers. For open land (here) and forest (below).
                 ! soilw_surf_GUESS_opl =(swaopl12(jl)*(dz1w+dz2w)    + swaopl3(jl)*(0.5_realkind-(dz1w+dz2w)))      /(dz1w+dz2w   +0.5_realkind-(dz1w+dz2w))
                 ! soilw_deep_GUESS_opl =(swaopl3(jl) *1.0_realkind   + soilw_surf_GUESS_opl*(0.5_realkind-(dz1w+dz2w))) /(1.0_realkind+0.5_realkind-(dz1w+dz2w))
                   soilw_surf_GUESS_opl =(swaopl12(jl)*0.282_realkind + swaopl3(jl)*0.218_realkind)      /0.5_realkind
                   soilw_deep_GUESS_opl =(swaopl3(jl)                 + soilw_surf_GUESS_opl*0.218_realkind) /1.218_realkind
           
                   soilw_surf_GUESS_for =(swafor12(jl)*0.282_realkind + swafor3(jl)*0.218_realkind)      /0.5_realkind
                   soilw_deep_GUESS_for =(swafor3(jl)                 + soilw_surf_GUESS_for*0.218_realkind) /1.218_realkind

                   snow_GUESS = max(0._realkind,dsolddt(jl))

                   lai_forunder(jl)=0.0_realkind

                   call lai_guess(guessid(jl), alat(jl), along(jl),   &
                        yearc,monthc,dayc,hourc,minc,secc,frland(jl), &
                        const_co2, t2mopsn(jl), tca(jl),              &
                        tsopsn1(jl),tsopsn2(jl),tsopsn3(jl),          &
                        tsopsn4(jl),                                  &
                        tsopsn5(jl),                                  &
                        tsfor1(jl),tsfor2(jl),tsfor3(jl),tsfor4(jl),  &
                        tsfor5(jl),                                   &
                        soilw_surf_GUESS_opl, soilw_deep_GUESS_opl,   &
                        soilw_surf_GUESS_for, soilw_deep_GUESS_for,   &
                        swrad_net_opl(jl), swrad_net_for(jl),         &
                        drolddt(jl), snow_GUESS,                      &
                        vegopl(jl),laiopn_int(jl),                    &
                        lai_conif(jl),lai_decid(jl),lai_forunder(jl), &
                        frdecid(jl),frfor(jl))

                   if(along(jl)>1.95_realkind .and. along(jl)<2.05_realkind)then  
                      if(alat(jl)>32.1_realkind.and. alat(jl)<32.2_realkind)then  
                         write(6,*) 'errsearch_phys_before: lon=', &
                              along(jl),       &
                              ' lat=',         &
                              alat(jl),        &
                              ' laiopn_int=',laiopn_int(jl),' ymdhm=',       &
                              yearc,monthc,dayc,hourc,minc,        &
                              ' vegopl=',vegopl(jl),     &
                              ' t2mopsn=',t2mopsn(jl),   &
                              ' tsopsn1=',tsopsn1(jl),   &
                              ' swaopl12=',swaopl12(jl), &
                              ' swaopl3=',swaopl3(jl)                 
                      endif
                   endif

                enddo

             endif

!RCA4
             call slfluxo_surf_sea_ice(nhor,nlev,kstart,kstop,     &
                  dtvdif,  &
                  wt,wq,wu,wv,pf,gpot,dph,     &
                  ps,tsea, &
                  tice,ticed,ticesn,ticesnd,tsnice,      &
                  tskin,snice,swsnice,rhosnice,snmaxice, &
                  frice,frsnice,frland,        &
                  drolddt,dsolddt,   &
                  rousea,  &
                  radf,albedo,albicenl,albsnice,wscos,sswdn,emskin,&
                  lcounttice,icethick,         &
                  along,coslat,kstep,&
                  t2ms,q2ms,t2mi,q2mi,         &
                  rh2ms,rh2mi,       &
                  u10ms,v10ms,u10mi,v10mi,     &
                  momfus,momfvs,momfui,momfvi, &
                  senfs,latfs,senfi,latfi,     &
                  evaps,evapi,accrunofflake,   &
                  sswrsea,slwrsea,sswri,slwri, &
                  wustars,wustari,   &
                  dzsnowice,         &
                  dticedt,dticeddt,  &
                  dticesndt,dticesnddt,dtsnicedt,        &
                  dsnicedt,dswsnicedt,drhosnicedt,dsnmaxicedt)             
!RCA35
!      call slfluxo_surf_sea_ice(nhor,nlev,kstart,kstop,  &
!      dtvdif,  &
!      wt,wq,wu,wv,pf,gpot,dph,  &
!      ps,tsea,  &
!      tice,ticed,ticesn,ticesnd,tsnice,  &
!      tskin,snice,swsnice,rhosnice,snmaxice,  &
!      frice,frsnice,frland,  &
!      drolddt,dsolddt,  &
!      rousea,  &
!      radf,albedo,albicenl,albsnice,wscos,sswdn,emskin,  &
!      lcounttice,icethick,  &
!      along,coslat,kstep,  &
!      t2ms,q2ms,t2mi,q2mi,  &
!      rh2ms,rh2mi,  &
!      u10ms,v10ms,u10mi,v10mi,  &
!      momfus,momfvs,momfui,momfvi,  &
!      senfs,latfs,senfi,latfi,  &
!      evaps,evapi,accrunofflake,  &
!      sswrsea,slwrsea,sswri,slwri,  &
!      wustars,wustari,  &
!      dticedt,dticeddt,  &
!      dticesndt,dticesnddt,dtsnicedt,  &
!      dsnicedt,dswsnicedt,drhosnicedt,dsnmaxicedt)


             call slfluxo_surf_lake_ice(nhor,nlev,kstart,kstop,        &
                  dtvdif, coslat,lcounttice,   &
                  wt,wq,wu,wv,gpot,dph,     &
                  ps,      &
                  prog_lake(:,:,1), prog_lake(:,:,2),prog_lake(:,:,3),       &
                  prog_lake(:,:,4), prog_lake(:,:,5),prog_lake(:,:,6),       &
                  prog_lake(:,:,7), prog_lake(:,:,8),prog_lake(:,:,9),       &
                  prog_lake(:,:,10),prog_lake(:,:,11),   &
                  prog_lake(:,:,12), &
                  frlake,  &
                  depthlake,         &
                  tskin,   &
                  drolddt,dsolddt,   &
                  radf,    &
                  albedo,  &
                  wscos,   &
                  sswdn,   &
                  emskin,  &
                  diag_lake(:,:,1),  &
                  t2ms,t2mi,         &
                  q2ms,q2mi,         &
                  rh2ms,rh2mi,       &
                  u10ms,u10mi,       &
                  v10ms,v10mi,       &
                  senfs,senfi,       &
                  latfs,latfi,       &
                  evaps,evapi,       &
                  accrunofflake,     &
                  momfus,momfui,     &
                  momfvs,momfvi,     &
                  wustars,wustari,   &
                  rousea,  &
                  tsnice,frsnice,    &
                  diag_lake(:,:,2),diag_lake(:,:,3),diag_lake(:,:,13),       &
                  diag_lake(:,:,4), diag_lake(:,:,5),    &
                  diag_lake(:,:,6), diag_lake(:,:,7) ,   &
                  diag_lake(:,:,8), diag_lake(:,:,9),    &
                  diag_lake(:,:,10), diag_lake(:,:,11),  &
                  diag_lake(:,:,12), &
                  tend_lake(:,:,1), tend_lake(:,:,2),tend_lake(:,:,3),       &
                  tend_lake(:,:,4), tend_lake(:,:,5),tend_lake(:,:,6),       &
                  tend_lake(:,:,7), tend_lake(:,:,8),tend_lake(:,:,9),       &
                  tend_lake(:,:,10),tend_lake(:,:,11),   &
                  tend_lake(:,:,12))
                                             
             !    write(*,*)'in phys before slfluxo_average,kstep=',kstep           
                         
             
!RCA4
             call slfluxo_average(nhor,nlev,kstart,kstop,&
                  wt,wq,wu,wv,gpot,dph,        &
                  ps,tsea,tice,frland,frice,   &
                  tsnice,frsnice,    &
                  rough,rousea,wz0land,        &
                  senfl,latfl,evapl, &
                  senfs,latfs,evaps,senfi,latfi,evapi,   &
                  wustarl,wustars,wustari,     &
                  t2ml,t2ms,t2mi,t2mopsn,      &
                  q2ml,q2ms,q2mi,q2mopsn,      &
                  rh2ml,rh2ms,rh2mi,rh2mopsn,  &
                  u10ms,u10mi,u10opsn,         &
                  v10ms,v10mi,v10opsn,         &
                  tseff,along,coslat,kstep,    &
                  momf,momfu,momfv,  &
                  t2,q2,rh2,u10,v10, &
                  t2mopsnsi,q2mopsnsi,rh2mopsnsi,        &
                  u10opsnsi,v10opsnsi,         &
                  senf,latf,evap,ustar,z0wind)                             
!RCA35
!      call slfluxo_average(nhor,nlev,kstart,kstop,  &
!      wt,wq,wu,wv,gpot,dph,  &
!      ps,tsea,tice,frland,frice,  &
!      tsnice,frsnice,  &
!      rough,rousea,wz0land,  &
!      senfl,latfl,evapl,  &
!      senfs,latfs,evaps,senfi,latfi,evapi,  &
!      wustarl,wustars,wustari,  &
!      t2ml,t2ms,t2mi,t2mopsn,  &
!      q2ml,q2ms,q2mi,q2mopsn,  &
!      rh2ml,rh2ms,rh2mi,rh2mopsn,  &
!            u10ms,u10mi,u10opsn,  &
!            v10ms,v10mi,v10opsn,  &
!      tseff,along,coslat,kstep,  &
!      momf,momfu,momfv,  &
!      t2,q2,rh2,u10,v10,  &
!      t2mopsnsi,q2mopsnsi,rh2mopsnsi,  &
!      u10opsnsi,v10opsnsi,  &
!      senf,latf,evap,ustar,z0wind)



             diag_lake(:,1,14)=frlake(:,1)                            
             diag_lake(:,2,14)=frlake(:,2)                            
             diag_lake(:,3,14)=frlake(:,3)                            
             diag_lake(:,1,15)=depthlake(:,1)                         
             diag_lake(:,2,15)=depthlake(:,2)                         
             diag_lake(:,3,15)=depthlake(:,3)                              

             do jl=kstart,kstop 
                mv10(jl)=sqrt(u10(jl)**2_realkind+v10(jl)**2_realkind)                  
                mv10opsnsi(jl)=sqrt(u10opsnsi(jl)**2_realkind+v10opsnsi(jl)**2_realkind)
                gustest(jl)=mv10(jl)                                  
                gustlow(jl)=mv10(jl)                                  
                gustup(jl)=mv10(jl)                                   
                gustestuc(jl)=mv10(jl)                                  

                if( frland(jl)>zfrlim ) then                       
                   zlandon=1.0_realkind  
                   zlandoff=0.0_realkind 
                else            
                   zlandon=0.0_realkind  
                   zlandoff=1.0_realkind 
                endif
                ziceon=1._realkind       
                ziceoff=0._realkind      
                zfrice=frice(jl)*(1._realkind-frland(jl))                      
                if(zfrice<=zfrlim)then                              
                   ziceon=0._realkind    
                   ziceoff=1._realkind   
                endif

                !    compute grid average snow water equivalent                                   
                !    snopl=open land snow average over land area                       
                !    snfor=forest snow average over land area                          
                !    snice=snow on ice average over ice area                           
                !    sn= grid average snow !!    

                zsnpland=snopl(jl)+dtime*dsnopldt(jl)+   &
                     snfor(jl)+dtime*dsnfordt(jl)                     
                zsnpland=frland(jl)*zsnpland + (1._realkind-frland(jl))*    &
                     frice(jl)*(snice(jl)+dtime*dsnicedt(jl))         
                zsafe=0._realkind        
                if((frland(jl)+(1._realkind-frland(jl))*frice(jl))<=0._realkind)then  
                   zsafe=1._realkind     
                endif
                sn(jl)=zsnpland/(frland(jl)+(1._realkind-frland(jl))*       &
                     frice(jl)+zsafe)                                 

                ! Average snow depth
                dzsnow(jl)=0._realkind
                if((wfrsnw(jl)+frsnfor(jl)+frsnice(jl))>0._realkind)then
                  dzsnow(jl)=(dzsnowopl(jl)*wfrsnw(jl)+dzsnowfor(jl)*frsnfor(jl)+dzsnowice(jl)*frsnice(jl))/  &
                             (wfrsnw(jl)+frsnfor(jl)+frsnice(jl))
                endif

                !    compute open land + snow averaged soil temperatures (tsopsn) and  
                !    forest averaged soil temperatures (tsfor)                         

                tsopsn1(jl)=(1._realkind-frsn(jl))*     &
                     (tsns(jl)+dtime*dtsnsdt(jl))+       &
                     frsn(jl)*(tssn(jl)+dtime*dtssndt(jl))            
                tsopsn2(jl)=(1._realkind-frsn(jl))*(tsns2(jl)+    &
                     dtime*dtsns2dt(jl))+      &
                     frsn(jl)*(tssn2(jl)+dtime*dtssn2dt(jl))          
                tsopsn3(jl)=(1._realkind-frsn(jl))*     &
                     (tsns3(jl)+dtime*dtsns3dt(jl))+     &
                     frsn(jl)*(tssn3(jl)+dtime*dtssn3dt(jl))          
                tsopsn4(jl)=(1._realkind-frsn(jl))*(tsns4(jl)+    &
                     dtime*dtsns4dt(jl))+      &
                     frsn(jl)*(tssn4(jl)+dtime*dtssn4dt(jl))          
                tsopsn5(jl)=(1._realkind-frsn(jl))*(tsns5(jl)+    &
                     dtime*dtsns5dt(jl))+      &
                     frsn(jl)*(tssn5(jl)+dtime*dtssn5dt(jl))          

                tsfor1(jl)=(1._realkind-frsnfor(jl))*(tsc(jl)+    &
                     dtime*dtscdt(jl))+        &
                     frsnfor(jl)*(tscsn(jl)+dtime*dtscsndt(jl))       
                tsfor2(jl)=(1._realkind-frsnfor(jl))*(tsc2(jl)+   &
                     dtime*dtsc2dt(jl))+       &
                     frsnfor(jl)*(tscsn2(jl)+dtime*dtscsn2dt(jl))     
                tsfor3(jl)=(1._realkind-frsnfor(jl))*(tsc3(jl)+   &
                     dtime*dtsc3dt(jl))+       &
                     frsnfor(jl)*(tscsn3(jl)+dtime*dtscsn3dt(jl))     
                tsfor4(jl)=(1._realkind-frsnfor(jl))*(tsc4(jl)+   &
                     dtime*dtsc4dt(jl))+       &
                     frsnfor(jl)*(tscsn4(jl)+dtime*dtscsn4dt(jl))     
                tsfor5(jl)=(1._realkind-frsnfor(jl))*(tsc5(jl)+   &
                     dtime*dtsc5dt(jl))+       &
                     frsnfor(jl)*(tscsn5(jl)+dtime*dtscsn5dt(jl))     

                tsland1(jl)=(frcw(jl)*tsfor1(jl)+        &
                     (1._realkind-frcw(jl))*tsopsn1(jl))*         &
                     zlandon + ztdumout*zlandoff                      
                tsland2(jl)=(frcw(jl)*tsfor2(jl)+        &
                     (1._realkind-frcw(jl))*tsopsn2(jl))*         &
                     zlandon + ztdumout*zlandoff                      
                tsland3(jl)=(frcw(jl)*tsfor3(jl)+        &
                     (1._realkind-frcw(jl))*tsopsn3(jl))*         &
                     zlandon + ztdumout*zlandoff                      
                tsland4(jl)=(frcw(jl)*tsfor4(jl)+        &
                     (1._realkind-frcw(jl))*tsopsn4(jl))*         &
                     zlandon + ztdumout*zlandoff                      
                tsland5(jl)=(frcw(jl)*tsfor5(jl)+        &
                     (1._realkind-frcw(jl))*tsopsn5(jl))*         &
                     zlandon + ztdumout*zlandoff                      

                tice1(jl)=(frsnice(jl)*ticesn(jl)+       &
                     (1._realkind-frsnice(jl))*tice(jl))*         &
                     ziceon+ztdumout*ziceoff                          
                tice2(jl)=(frsnice(jl)*ticesnd(jl)+      &
                     (1._realkind-frsnice(jl))*ticed(jl))*        &
                     ziceon+ztdumout*ziceoff                          

                accrunoff(jl)=accrunoffopl(jl)+accrunofffor(jl)+   &
                     accrunofflake(jl)                                

             enddo

          endif

          !    boundary layer height (subroutine 'apblhgt')                      
          !    the computed boudary layer height is purely diagnostic            
          !    for later use as output.    
          !    ( turbulence parameterization needs to be active )                

          call pblhgt( nhor,nlev,kstart,kstop,&
               ustar,ps,tseff,senf,latf,pf,dph,gpot,     &
               wt,wq,wcw,wu,wv,pblh,wheatv,wobkl  )                        

          !    turbulence parameterization 
          !    noptions(3) used for choosing among turbulence schemes.           
          !    noptions(3)=2  means cbr turbulence scheme                        
          !    noptions(3)=1  means holtslag scheme.                             
          if( noption(3)==2 ) then                                  
! RCA4
             call vcbr(nhor,nlev,kstart,kstop,dtvdif,  &
                  wt,wq,wu,wv,wcw,wsv(1,1,1),gpot,       &
                  pf,ph,dpf,dph,     &
                  tseff, ps,momfu,momfv,       &
                  senf,latf,ustar,   &
                  dtdtpr ,dqdtpr, dudtpr ,dvdtpr,        &
                  dcdtpr,dsvdtpr(1,1,1),       &
!                  stddif,stdif,      &
                  lesat,pblh,om_cpbl,zfrpbl,zfrtop,zlevtop,        &
                  z0wind,gustest,gustlow,gustup,gustestuc,         &
                  zvarcu,wdzf)         

! RCA3
!             call vcbr(nhor,nlev,kstart,kstop,dtvdif,  &
!                  lstat,wt,wq,wu,wv,wcw,wsv(1,1,1),gpot,       &
!                  pf,ph,dpf,dph,     &
!                  tseff, ps,momfu,momfv,       &
!                  senf,latf,ustar,   &
!                  dtdtpr ,dqdtpr, dudtpr ,dvdtpr,        &
!                  dcdtpr,dsvdtpr(1,1,1),       &
!                  stddif,stdif,      &
!                  pblh,om_cpbl,zfrtop,zlevtop,        &
!                  z0wind,gustest,gustlow,gustup) 



             do l=1,nsvar       
                do jk=1,nlev    
                   do jl=kstart,kstop                                 
                      wsv(jl,jk,l) = wsv(jl,jk,l) +      &
                           dtvdif*dsvdtpr(jl,jk,l)                    
                      dqdt_cbr(jl,jk)=dqdtpr(jl,jk)                   
                      dtdt_cbr(jl,jk)=dtdtpr(jl,jk)                   
                      dcwdt_cbr(jl,jk)=dcdtpr(jl,jk)                  
                   enddo
                enddo
             enddo
          endif

          if( noption(3)==1 ) then                                  
             !    ( holtslag scheme )         
             call avdifnl(nhor,nlev,kstart,kstop,        &
                  dtvdif,  &
                  wt,wq,wcw,totcov,wu,wv,gpot, &
                  wheatv,wobkl,      &
                  pf,ph,dpf,dph,     &
                  wts,ps,t2,q2,wmomfu,wmomfv,  &
                  wsenf,wlatf,ustar,frland, &
                  wustarg,pblh,rough,&
                  dtdtpr,dqdtpr,dcdtpr,dudtpr,dvdtpr)!,    &
             !stddif,stdif)                        
          endif
          
          do jk=1,nlev          
             do jl=kstart,kstop 
                wt(jl,jk) = wt(jl,jk) + dtvdif*dtdtpr(jl,jk)          
                wq(jl,jk) = wq(jl,jk) + dtvdif*dqdtpr(jl,jk)          
                wu(jl,jk) = wu(jl,jk) + dtvdif*dudtpr(jl,jk)          
                wv(jl,jk) = wv(jl,jk) + dtvdif*dvdtpr(jl,jk)          
                wcw(jl,jk) = wcw(jl,jk) + dtvdif*dcdtpr(jl,jk)        
             enddo
          enddo
90     enddo

       !    sum up tendencies (all in t,q-points)                             

       if(ldynvd) then          
          do l=1,nsvar          
             do jk=1,nlev       
                do jl=kstart,kstop                                    
                   dsvardt(jl,jk,l)= (wsv(jl,jk,l)-svar(jl,jk,l))/dtime
                enddo
             enddo
          enddo

          do jk=1,nlev          
             do jl=kstart,kstop 
                dtdt(jl,jk)=(wt(jl,jk)-t(jl,jk))/dtime                
                dqdt(jl,jk)=(wq(jl,jk)-q(jl,jk))/dtime                
                dudt(jl,jk)=(wu(jl,jk)-uvel(jl,jk))/dtime                
                dvdt(jl,jk)=(wv(jl,jk)-vvel(jl,jk))/dtime                
                dcwdt(jl,jk)=(wcw(jl,jk)-cw(jl,jk))/dtime             
             enddo
          enddo
       else                     
          do l=1,nsvar          
             do jk=1,nlev       
                do jl=kstart,kstop                                    
                   dsvardt(jl,jk,l)=dsvardt(jl,jk,l)+(wsv(jl,jk,l)-svar(jl,jk,l))/dtime            
                enddo
             enddo
          enddo
          do jk=1,nlev          
             do jl=kstart,kstop 
                dtdt(jl,jk)=dtdt(jl,jk)+(wt(jl,jk)-t(jl,jk))/dtime    
                dqdt(jl,jk)=dqdt(jl,jk)+(wq(jl,jk)-q(jl,jk))/dtime    
                dudt(jl,jk)=dudt(jl,jk)+(wu(jl,jk)-uvel(jl,jk))/dtime    
                dvdt(jl,jk)=dvdt(jl,jk)+(wv(jl,jk)-vvel(jl,jk))/dtime    
                dcwdt(jl,jk)=dcwdt(jl,jk)+(wcw(jl,jk)-cw(jl,jk))/dtime 
             enddo
          enddo
       endif
    endif

    if( .not.lsamedt ) then     
       if( lnewpht ) then       
          do 150 jk=1,nlev      
             do jl=kstart,kstop                                   
                dtdtph(jl,jk)=dtdt(jl,jk)                             
                dqdtph(jl,jk)=dqdt(jl,jk)                             
             enddo
             if( locw ) then    
                do jl=kstart,kstop                                
                   dcdtph(jl,jk)=dcwdt(jl,jk)                         
                enddo
             endif
150       enddo
       else                     
          do 170 jk=1,nlev      
             do jl=kstart,kstop                                   
                dtdt(jl,jk)=dtdt(jl,jk)+dtdtph(jl,jk)                 
                dqdt(jl,jk)=dqdt(jl,jk)+dqdtph(jl,jk)                 
             enddo
             if( locw ) then    
                do jl=kstart,kstop                                
                   dcwdt(jl,jk)=dcwdt(jl,jk)+dcdtph(jl,jk)            
                enddo
             endif
170       enddo

          !    correction for negative humidity                                  
          call aqnegat(nhor,nlev,kstart,kstop, &
               dtime,      &
               dqdt,q,dpf, &
               dqdtpr)          

          do  jk=1,nlev      
             do  jl=kstart,kstop                                   
                dqdt(jl,jk)=dqdt(jl,jk)+dqdtpr(jl,jk)                 
             enddo
          enddo
          return                
       endif
    endif

    !    compute length of accumulation time step 'dtacc'                  

    dtacc=conacc*dtime               

    !    put rain and snow intencities to zero                             
    do jl=kstart,kstop      
       draindt(jl)=0._realkind           
       dsnowdt(jl)=0._realkind           
    enddo



    !    condensation ( noption(2) ) 
    !    noption(2)=1  for the straco condensation scheme                  
    !    noption(2)=2  for hirlam (tiedtke) mass flux                      
    !    scheme plus sundqvist stratiform scheme                           
    !    gjkf                        
    !    gjkf  noption(2)=4  kain-fritsch convection + rak ls condensation 
    !    gjkf                        
    if( noption(2)==4 ) then  
       !    here wts is valid for the temperature of everything but sea       
       do jl=kstart,kstop       
          zfrsea=(1._realkind-frland(jl))*(1._realkind-frice(jl))                       
          if(zfrsea>(1._realkind-zfrlim))then                               
             wts(jl)=tseff(jl)  
          else                  
             wts(jl)=(tseff(jl)-zfrsea*tsea(jl))/(1._realkind-zfrsea)          
          endif
       enddo

       call akfrak(nhor,nlev,kstart,kstop,     &
            dtime,conacc,  &
            frland,frice,tseff,ps,     &
            ahyb,bhyb,     &
            uvel,vvel, &
            t,q,cw,svar,pf,dpf,      &
            dcwdt,dtdt,dqdt,         &
            cov2d,         &
            prcpst,stsnow, &
            draindt,dsnowdt,         &
            dtdtpr,dqdtpr,wt,        &
            totcov,cucov,  &
            preta,     &       !stratiform 3d precip field
            orosigm,       &
            raincv_kf,     &
            snowcv_kf,     &
            frsngrid,      &
            dph, &
            zdx_kf,        &
            zdy_kf,        &
            oldcov,        &
            omf,div_kf,    &
            lesat,pblh,    &
            om_cpbl,zfrpbl,zfrtop,zlevtop,     &
            relh,udr,udrf,      &
            zvarcu,        &
            gpot,&
            wdzf,&
            kf_ind,shal_cgj,lsb,     &
!cgj300611
            cwcu,rhcrit)             

       !    gjkf dtdt,dqdt&dcwdt on entry to acondens were named dtdtin,dqdtin & dcwdin                             
       !    gjkf internal to acondens and on entry contained the tends due to dyn,rad                               
       !    gjkf and vdiff. inside acondens the tends due to convection get added in                                
       !    gjkf to these terms. the terms  dtdtpr,dqdtpr,wt are dtdt,dqdt & dcwdt                                  
       !    gjkf internally in acondens. these on exit from condense contain only                                   
       !    gjkf tends due to stratiform condensation.                        
       !    gjkf                        
       !    gjkf  save conv prec elements to cusnow & prcpcu which go out as  
       !    gjkf  zprc znd zsnoc into phtask for accumulation                 
       !    gjkf                        
       do jl=kstart,kstop       
          cusnow(jl)=snowcv_kf(jl)                                    
          prcpcu(jl)=raincv_kf(jl)                                    
       enddo

       do jk=1,nlev             
          do jl=kstart,kstop    
             dtdt(jl,jk)=dtdt(jl,jk)+dtdtpr(jl,jk)                    
             dqdt(jl,jk)=dqdt(jl,jk)+dqdtpr(jl,jk)                    
             dcwdt(jl,jk)=dcwdt(jl,jk)+wt(jl,jk)                      
!wdzf
!!!             zducld=udrf(jl,jk)*5.e04_realkind
!!!
!                   alat(jl)=acos(coslat(jl))*180._realkind/3.141592654_realkind         
!                   if(sinlat(jl)<0._realkind)alat(jl)=-alat(jl)             
!                   if(along(jl)>10.0_realkind .and. along(jl)<12.0_realkind )then 
!                      if(alat(jl)>0.0_realkind .and. alat(jl)< 1.5_realkind)then
!!                       write(*,232)rhcrit(jl,jk),wdzf(jl,jk),pf(jl,jk), &
!!                                   frland(jl),jl,jk
!                       write(*,233)relh(jl,jk),rhcrit(jl,jk),  &
!                                   totcov(jl,jk),cucov(jl,jk),&
!                                   pf(jl,jk),wdzf(jl,jk),frland(jl), &
!                                   alat(jl),along(jl),jk       
!                      endif
!                      if(alat(jl)>10.0_realkind .and. alat(jl)< 11.5_realkind)then
!!                       write(*,232)rhcrit(jl,jk),wdzf(jl,jk),pf(jl,jk), &
!!                                   frland(jl),jl,jk
!                       write(*,233)relh(jl,jk),rhcrit(jl,jk),  &
!                                   totcov(jl,jk),cucov(jl,jk),&
!                                   pf(jl,jk),wdzf(jl,jk),frland(jl), &
!                                   alat(jl),along(jl),jk       
!                      endif
!                      if(alat(jl)>20.0_realkind .and. alat(jl)< 21.5_realkind)then
!!                       write(*,232)rhcrit(jl,jk),wdzf(jl,jk),pf(jl,jk), &
!!                                   frland(jl),jl,jk
!                       write(*,233)relh(jl,jk),rhcrit(jl,jk),  &
!                                   totcov(jl,jk),cucov(jl,jk),&
!                                   pf(jl,jk),wdzf(jl,jk),frland(jl), &
!                                   alat(jl),along(jl),jk       
!                      endif
!!232                   format(5f14.5,3x,1e13.5,3x,2i4)
!233                   format(5f12.4,3x,1e13.5,3x,3f12.4,3x,i4)
!                   endif
          enddo
       enddo
    endif                     !end akfrak                                  

    if( noption(2)==2 ) then  
       !    here wt is temporary used for output of dcwdt                     

       call acondens  (nhor,nlev,kstart,kstop, &
            dtime,         &
            frland,tseff,ps,dpsdin,  &
            bhyb,&
            bfull,         &
            ph,gpot,senf,latf,uvel,vvel,omf,         &
            t,q,cw,pf,dpf, &
            dcwdt,dtdt,dqdt,dudt,dvdt,         &
            cov2d,cwpath,kht,        &
            sevapr,sevapc, &
            prcpst,stsnow,cusnow,prcpcu,       &
            draindt,dsnowdt,         &
            dtdtpr,dqdtpr,wt,dudtpr,dvdtpr,    &
            totcov,cucov,  &
            preta,prsl,prsi,prcl,prci)                                
       sevapr = sevapr*dtacc    
       sevapc = sevapc*dtacc    

       acprcu = 0._realkind              
       acprst = 0._realkind              
       acsncu = 0._realkind              
       acsnst = 0._realkind              

       do jl = kstart,kstop     
          acprcu     = acprcu + dtacc*prcpcu(jl)                      
          acprst     = acprst + dtacc*prcpst(jl)                      
          acsncu     = acsncu + dtacc*cusnow(jl)                      
          acsnst     = acsnst + dtacc*stsnow(jl)                      
       enddo


       !    sum tendencies of t,q and store tendencies in t,q-points          
       !    of u and v accumulate precipitation at model levels.              

       do jk=1,nlev             
          do jl=kstart,kstop    
             dtdt(jl,jk)=dtdt(jl,jk)+dtdtpr(jl,jk)                    
             dqdt(jl,jk)=dqdt(jl,jk)+dqdtpr(jl,jk)                    
             dcwdt(jl,jk)=dcwdt(jl,jk)+wt(jl,jk)                      
          enddo
       enddo

    endif

    !    call straco condensation scheme                                   
    !    (first argument = 1)        
    if( noption(2)==1 ) then  
       call aconds(1,nhor,nlev,kstart,kstop,kstep,       &
            dtime,conacc,dtheta,     &
            t,q,cw,dcwdt,dtdt,dqdt,pf,dpf,     &
            ps,dpsdin,     &
            cov2d,cwpath,draindt,dsnowdt,      &
            preta,prcpst,prcpcu,stsnow,cusnow, &
            ahyb,bhyb,&
            dtdtpr,dqdtpr,wt,totcov,cucov)                            

       do jk=1,nlev             
          do jl=kstart,kstop    
             dtdt(jl,jk)=dtdt(jl,jk)+dtdtpr(jl,jk)                    
             dqdt(jl,jk)=dqdt(jl,jk)+dqdtpr(jl,jk)                    
             dcwdt(jl,jk)=dcwdt(jl,jk)+wt(jl,jk)                      
          enddo
       enddo
    endif

    !    the following loops are not currently used                        
    do jk=1,nlev                
       do jl=kstart,kstop       
          virt(jl,jk)=1._realkind+0.607717_realkind*q(jl,jk) -cw(jl,jk)                 
       enddo
    enddo

    !    compute tendencies of q due to no humidity < 0                    
    call aqnegat(nhor,nlev,kstart,kstop,       &
         dtime,  &
         dqdt,q,dpf,       &
         dqdtpr)

    !     sum up tendencies of q
    do jk=1,nlev
       do jl=kstart,kstop
          dqdt(jl,jk)=dqdt(jl,jk)+dqdtpr(jl,jk)
       enddo
    enddo

    !     store tendencies for use next timestep

    do jl=kstart,kstop
       drolddt(jl)=draindt(jl)
       dsolddt(jl)=dsnowdt(jl)
    enddo

    !     compute tendencies due to all of phys
    do jk=1,nlev
       do jl=kstart,kstop
          dtp(jl,jk) =  (dtdt(jl,jk)-dtdtin(jl,jk))
          dqp(jl,jk) =  (dqdt(jl,jk)-dqdtin(jl,jk))
          dcp(jl,jk) = (dcwdt(jl,jk)-dcdtin(jl,jk))
       enddo
    enddo

    if( lsamedt ) return
    if( lnewpht ) then
       do jk=1,nlev
          do jl=kstart,kstop
             dtdtph(jl,jk)=dtdt(jl,jk)-dtdtph(jl,jk)
             dqdtph(jl,jk)=dqdt(jl,jk)-dqdtph(jl,jk)
          enddo

          if( locw ) then
             do jl=kstart,kstop
                dcdtph(jl,jk) = dcwdt(jl,jk) - dcdtph(jl,jk)
             enddo
          endif
       enddo

    endif
    return
  end subroutine phys


  subroutine phtask(kstep,year,month,day,hour,mnt,sec,            &
       nhorph,nlon,nlat,nlev,kstart,kstop,                         &
       dtime,nsvdif,conacc,dtheta, &
       lpp,            &
       lsamedt,lnewpht,ldynvd,                                     &
       ahalf,bhalf,hybf,hybh,&
       t, q, cw, u, v, omf,  &
       ps,                   &
       tsea, frice,          &
       along, coslat, sinlat, rough, rousea,                       &
       drolddt,dsolddt,      &
       accsunny,dtphysh,     &
       accrunoff,accrunoffopl,accrunofffor,                        &
       accrunofflake,        &
       q2d,                  &
       orogsigm,             &
       accprl,accprc,        &
       dtdt, dqdt, dcwdt, dudt, dvdt,                              &
       totcov,  cucov, cov2d,  cwpath,                             &
       dpsdt,                &
       dtdtph, dqdtph,dcdtph,&
       pblh,                 &
       slwr,sswr,tlwr,tswr,  &
       slwrsea,sswrsea,slwrice,sswrice,                            &
       accsnow,accsnoc,accsnol,                                    &
       slwdn,sswdn,tswdn,    &
       svars, dsvarsdt, nsvars,                                    &
       icethick,             &
       guessid,                                                    &
       lcounttice,           &
       msvars,svar_surf,     &
       eco,meco,             &
       lake_types,lake_no_prog,                                    &
       lake_no_tend,lake_no_diag,                                  &
       frac_lakes,depth_lakes,                                     &
       prog_lakes,tend_lakes,diag_lakes,                           &
       tsclim_years,         &
       svar,dsvardt,nsvar,   &
       sacdg,nsacdg,sacdg2,nsacdg2,                                &
       dtdt_kf      ,        &
       dqdt_kf      ,        &
       dqcdt_kf     ,        &
       oldcov       ,        &
       zvarcu       ,        &
!cgj300611
       cwcu         ,        &
       div_kf       ,        &
       om_rmean     ,        &
       kf_ind       ,        &
       kf_base      ,        &
       kf_top       ,        &
       raincv_kf    ,        &
       snowcv_kf    ,        &
       umfb         ,        &
       shal_cgj     ,        &
       nca          ,        &
       zdx_kf       ,        &
       zdy_kf       ,        &
       afull        ,        &
       bfull,ztropo,ttropo)                


    implicit none                    

    integer kstep,year,month,day,hour,mnt,sec                         
    integer nhorph,nlon,nlat,nlev,kstart,kstop                        

    real(kind=realkind) dtime                  
    integer nsvdif              
    real(kind=realkind) conacc,dtheta          

    logical lpp,         &
         lsamedt,lnewpht,ldynvd 
    real(kind=realkind) ahalf(nlev+1),bhalf(nlev+1),hybf(nlev),hybh(nlev+1)          
    real(kind=realkind) t(nlon*nlat,nlev),q(nlon*nlat,nlev),                        &
         cw(nlon*nlat,nlev),u(nlon*nlat,nlev),                       &
         v(nlon*nlat,nlev),omf(nlon*nlat,nlev),ps(nlon*nlat),        &
         tsea(nlon*nlat),frice(nlon*nlat),                           &
         along(nlon*nlat),     &
         coslat(nlon*nlat),sinlat(nlon*nlat),                        &
         rousea(nlon*nlat),    &
         rough(nlon*nlat)       
    real(kind=realkind) slwr(nlon*nlat),sswr(nlon*nlat),                            &
         tlwr(nlon*nlat),tswr(nlon*nlat)                              
    real(kind=realkind) slwrsea(nlon*nlat),sswrsea(nlon*nlat),                      &
         slwrice(nlon*nlat),sswrice(nlon*nlat)                        
    real(kind=realkind) accsnow(nlon*nlat)     
    real(kind=realkind) accsnol(nlon*nlat)     
    real(kind=realkind) accsnoc(nlon*nlat)     
    real(kind=realkind) slwdn(nlon*nlat),sswdn(nlon*nlat),tswdn(nlon*nlat)           

    real(kind=realkind) dtdt(nlon*nlat,nlev),dqdt(nlon*nlat,nlev),                  &
         dcwdt(nlon*nlat,nlev),&
         dudt(nlon*nlat,nlev),dvdt(nlon*nlat,nlev),                  &
         totcov(nlon*nlat,nlev),cucov(nlon*nlat,nlev),               &
         cov2d(nlon*nlat),cwpath(nlon*nlat),                         &
         dpsdt(nlon*nlat)       
    real(kind=realkind) dtdtph(nlon*nlat,nlev),dqdtph(nlon*nlat,nlev),              &
         dcdtph(nlon*nlat,nlev) 
    real(kind=realkind) svar(nlon*nlat,nlev,*),dsvardt(nlon*nlat,nlev,*)             
    real(kind=realkind) accsunny(nlon*nlat),dtphysh                                  
    real(kind=realkind) accprl(nlon*nlat),accprc(nlon*nlat),                        &
         drolddt(nlon*nlat),dsolddt(nlon*nlat),                      &
         accrunoff(nlon*nlat), &
         accrunoffopl(nlon*nlat),accrunofffor(nlon*nlat),            &
         accrunofflake(nlon*nlat),                                   &
         q2d(nlon*nlat),       &
         orogsigm(nlon*nlat)    

    integer nsvar               

    integer,intent(in):: nsacdg,nsacdg2      
    real(kind=realkind):: sacdg(nlon*nlat,nlev,nsacdg)                                 
    real(kind=realkind):: sacdg2(nlon*nlat,nsacdg2)                                    
    real(kind=realkind):: pblh(nlon*nlat)        

    real(kind=realkind):: dtdt_kf(nlon*nlat,nlev),                                    &
         dqdt_kf(nlon*nlat,nlev),                                    &
         dqcdt_kf(nlon*nlat,nlev),                                   &
         oldcov(nlon*nlat,nlev),                                     &
         zvarcu(nlon*nlat,nlev),                                     &
         cwcu(nlon*nlat,nlev),                                     &
         div_kf(nlon*nlat,nlev),                                     &
         om_rmean(nlon*nlat,nlev),                                   &
         raincv_kf(nlon*nlat), &
         snowcv_kf(nlon*nlat), &
         umfb(nlon*nlat),      &
         zdx_kf(nlon*nlat),    &
         zdy_kf(nlon*nlat),    &
         ttropo(nlon*nlat),    &
         ztropo(nlon*nlat),    &
         afull(nlev),bfull(nlev)
    integer nca(nlon*nlat),    &
         kf_ind(nlon*nlat),    &
         kf_base(nlon*nlat),   &
         kf_top(nlon*nlat),    &
         shal_cgj(nlon*nlat)    
    integer  nsvars             
    real(kind=realkind),intent(in):: lcounttice(nlon*nlat)  
    real(kind=realkind) svars(nlon*nlat,*), dsvarsdt(nlon*nlat,*),                  &
         icethick(nlon*nlat)    
    integer guessid(nlon*nlat)
    integer  msvars             
    real(kind=realkind)   svar_surf(nlon*nlat,*)    

    integer  meco               
    real(kind=realkind)     eco(nlon*nlat,*)   

    integer lake_types,lake_no_prog,lake_no_tend,lake_no_diag         
    real(kind=realkind) frac_lakes(nlon*nlat,lake_types),                           &
         depth_lakes(nlon*nlat,lake_types),                          &
         prog_lakes(nlon*nlat,lake_types,lake_no_prog),              &
         tend_lakes(nlon*nlat,lake_types,lake_no_tend),              &
         diag_lakes(nlon*nlat,lake_types,lake_no_diag)                     

    real(kind=realkind) tsclim_years(nlon*nlat)     

    real(kind=realkind)   ttsk(nhorph,nlev),    qtsk(nhorph,nlev),                  &
         cwtsk(nhorph,nlev),    utsk(nhorph,nlev),                   &
         vtsk(nhorph,nlev),   dttsk(nhorph,nlev),                    &
         dqtsk(nhorph,nlev),  dcwtsk(nhorph,nlev),                   &
         dutsk(nhorph,nlev),   dvtsk(nhorph,nlev),                   &
         totctsk(nhorph,nlev),  cuctsk(nhorph,nlev),                 &
         dtphtsk(nhorph,nlev), dqphtsk(nhorph,nlev),                 &
         dcphtsk(nhorph,nlev),  omftsk(nhorph,nlev)                   
    real(kind=realkind) svtsk(nhorph,nlev,nsvar), dsvtsk(nhorph,nlev,nsvar)          
    real(kind=realkind) zprtsk(nhorph,nlev)    
    real(kind=realkind) zttsk(nhorph,nlev),ztrtsk(nhorph,nlev),ztttsk(nhorph,nlev), &
         zqtsk(nhorph,nlev),zqttsk(nhorph,nlev),                     &
         zctsk(nhorph,nlev),zcttsk(nhorph,nlev)                       
    real(kind=realkind) zprl(nhorph),zprc(nhorph),                                  &
         zsnoc(nhorph),zsnol(nhorph)                                  
    integer kpar,klev,k,jx,jy,l,jc                                    
    real(kind=realkind) dtacc                  
    real(kind=realkind) omega_tsk(nhorph,nlev)   ,            &
         dtdt_kftsk(nhorph,nlev)  ,                                  &
         dqdt_kftsk(nhorph,nlev)  ,                                  &
         dqcdt_kftsk(nhorph,nlev) ,                                  &
         div_tsk(nhorph,nlev)    , &
         oldcov_tsk(nhorph,nlev)  ,                                  &
         zvarcu_tsk(nhorph,nlev) , &                                     
         cwcu_tsk(nhorph,nlev)                                     
    real(kind=realkind) svtsks(nhorph,nsvars), dsvtsks(nhorph,nsvars)                
    real(kind=realkind) svar_tsks(nhorph,msvars+1)  

    real(kind=realkind) eco_tsks(nhorph,meco+1)
    real(kind=realkind) frac_l_tsks(nhorph,lake_types),                             &
         depth_l_tsks(nhorph,lake_types),                            &
         prog_l_tsks(nhorph,lake_types,lake_no_prog),                &
         tend_l_tsks(nhorph,lake_types,lake_no_tend),                &
         diag_l_tsks(nhorph,lake_types,lake_no_diag)                  
    real(kind=realkind) lakemax1,lakemax9  

    do 90 klev=1,nlev           
       do 80 k=kstart,kstop     
          ttsk(k-kstart+1,klev)  =    t(k,klev)                       
          dttsk(k-kstart+1,klev) = dtdt(k,klev)                       
          qtsk(k-kstart+1,klev)  =    q(k,klev)                       
          omftsk(k-kstart+1,klev)  =    omf(k,klev)                   
          dqtsk(k-kstart+1,klev) = dqdt(k,klev)                       
          cwtsk(k-kstart+1,klev)   =     cw(k,klev)                   
          dcwtsk(k-kstart+1,klev)  =  dcwdt(k,klev)                   
          totctsk(k-kstart+1,klev) = totcov(k,klev)                   
          cuctsk(k-kstart+1,klev)  =  cucov(k,klev)                   
          !    cf   initialize 3-d tendencies and fluxes:                        
          zprtsk(k-kstart+1,klev) = 0._realkind                                
          zttsk(k-kstart+1,klev) = 0._realkind                                 
          ztrtsk(k-kstart+1,klev) = 0._realkind                                
          ztttsk(k-kstart+1,klev) = 0._realkind                                
          zqtsk(k-kstart+1,klev) = 0._realkind                                
          zqttsk(k-kstart+1,klev) = 0._realkind                                
          zctsk(k-kstart+1,klev) = 0._realkind                                 
          zcttsk(k-kstart+1,klev) = 0._realkind                                

          omega_tsk(k-kstart+1,klev) = om_rmean(k,klev)               
          dtdt_kftsk(k-kstart+1,klev) = dtdt_kf(k,klev)               
          dqdt_kftsk(k-kstart+1,klev) = dqdt_kf(k,klev)               
          dqcdt_kftsk(k-kstart+1,klev) = dqcdt_kf(k,klev)             
          div_tsk(k-kstart+1,klev) = div_kf(k,klev)                   
          oldcov_tsk(k-kstart+1,klev) = oldcov(k,klev)                
          zvarcu_tsk(k-kstart+1,klev) = zvarcu(k,klev)                
          cwcu_tsk(k-kstart+1,klev) = cwcu(k,klev)                
80     enddo
90  enddo

    do jc=1,nsvars              
       do k=kstart,kstop        
          svtsks(k-kstart+1,jc)  =    svars(k,jc)                     
          dsvtsks(k-kstart+1,jc) = dsvarsdt(k,jc)                     
       enddo
    enddo

    do l=1,nsvar                
       do klev=1,nlev           
          do k=kstart,kstop     
             svtsk(k-kstart+1,klev,l)  =    svar(k,klev,l)            
             dsvtsk(k-kstart+1,klev,l) = dsvardt(k,klev,l)            
          enddo
       enddo
    enddo
    do jc=1,msvars              
       do k=kstart,kstop        
          svar_tsks(k-kstart+1,jc)  =    svar_surf(k,jc)              
       enddo
    enddo

    do jc=1,meco                
       do k=kstart,kstop        
          eco_tsks(k-kstart+1,jc)  =    eco(k,jc)                     
       enddo
    enddo
    do jc=1,lake_types          
       do k=kstart,kstop        
          frac_l_tsks(k-kstart+1,jc)  = frac_lakes(k,jc)              
          depth_l_tsks(k-kstart+1,jc) = depth_lakes(k,jc)             
          do l = 1,lake_no_prog 
             prog_l_tsks(k-kstart+1,jc,l) = prog_lakes(k,jc,l)        
          enddo
          do l = 1,lake_no_tend 
             tend_l_tsks(k-kstart+1,jc,l) = tend_lakes(k,jc,l)        
          enddo
          do l = 1,lake_no_diag 
             diag_l_tsks(k-kstart+1,jc,l) = diag_lakes(k,jc,l)        
          enddo
       enddo
    enddo
    if( .not.lsamedt .and. .not.lnewpht ) then                        
       do 96 klev=1,nlev        
          do 94 k=kstart,kstop  
             dtphtsk(k-kstart+1,klev) = dtdtph(k,klev)                
             dqphtsk(k-kstart+1,klev) = dqdtph(k,klev)                
             dcphtsk(k-kstart+1,klev) = dcdtph(k,klev)                
94        enddo
96     enddo
    endif

    do 180 klev=1,nlev       
       do 160 k=kstart,kstop 
          jy = (k-1)/nlon+1  
          jx = k-(jy-1)*nlon 
          dutsk(k-kstart+1,klev) = dudt(k,klev)                    
          dvtsk(k-kstart+1,klev) = dvdt(k,klev)                    
          if(jx==1) then   
             utsk(k-kstart+1,klev)  =       u(k,klev)              
          else               
             utsk(k-kstart+1,klev)  =   0.5_realkind*(  u(k,klev)          &
                  + u(k-1,klev) )                                  
          endif
          if(jy==1) then   
             vtsk(k-kstart+1,klev)  =       v(k,klev)              
          else               
             vtsk(k-kstart+1,klev) =   0.5_realkind*( v(k,klev)            &
                  + v(k-nlon,klev) )                               
          endif
160    enddo
180 enddo



    call phys(nhorph,nlev,1,kstop-kstart+1,                          &
         kstep,          &
         year,month,day,hour,mnt,sec,                                &
         dtime,nsvdif,conacc,dtheta,                                 &
         lsamedt,lnewpht,ldynvd,                       &
         ttsk,qtsk,cwtsk,utsk,vtsk,omftsk,                           &
         ps(kstart),tsea(kstart),frice(kstart),                      &
         along(kstart),coslat(kstart),sinlat(kstart),                &
         dpsdt(kstart),        &
         cov2d(kstart),cwpath(kstart),                               &
         rough(kstart),rousea(kstart),                               &
         zprl,zprc, zsnoc,zsnol,                                     &
         drolddt(kstart),dsolddt(kstart),                            &
         accsunny(kstart),dtphysh,                                   &
         orogsigm(kstart),                               &
         accrunoff(kstart),accrunoffopl(kstart),                     &
         accrunofffor(kstart),accrunofflake(kstart),                 & 
                                !    following are prognostic variables stored in arrays svarsm/z/p    
                                !    variable description is found below.                              
                                !    ksvars in gemini should correspond to current number of variables 
                                !    in this list.              
         svtsks(1,19),svtsks(1,20),svtsks(1,21),svtsks(1,22),        &
         svtsks(1,23),svtsks(1,24),svtsks(1,25),svtsks(1,26),        &
         svtsks(1,27),svtsks(1,28),svtsks(1,29),svtsks(1,30),        &
         svtsks(1,31),svtsks(1,32),svtsks(1,33),svtsks(1,34),        &
         svtsks(1,35),         &
         svtsks(1,36),svtsks(1,37),svtsks(1,38),svtsks(1,39),        &
         svtsks(1,40),svtsks(1,41),svtsks(1,42),svtsks(1,43),        &
         svtsks(1,44),svtsks(1,45),                                  &
         svtsks(1,46),svtsks(1,47),svtsks(1,48),svtsks(1,49),        &
         svtsks(1,1),svtsks(1,2),svtsks(1,3),                        &
         svtsks(1,4),svtsks(1,5),svtsks(1,6),                        &
         svtsks(1,7),svtsks(1,8),svtsks(1,9),                        &
         svtsks(1,10),svtsks(1,11),svtsks(1,12),                     &
         svtsks(1,13),svtsks(1,14),svtsks(1,15),svtsks(1,16),        &
         svtsks(1,17),svtsks(1,18),                                   &
                                !    following are time tendency variables stored in arrays dsvarsdt   
                                !    variable description is found below.                              
                                !    ksvars in gemini should correspond to current number of variables 
                                !    in this list.              
         dsvtsks(1,1),dsvtsks(1,2),dsvtsks(1,3),                     &
         dsvtsks(1,4),dsvtsks(1,5),dsvtsks(1,6),                     &
         dsvtsks(1,7),dsvtsks(1,8),dsvtsks(1,9),                     &
         dsvtsks(1,10),dsvtsks(1,11),dsvtsks(1,12),                  &
         dsvtsks(1,13),dsvtsks(1,14),dsvtsks(1,15),dsvtsks(1,16),    &
         dsvtsks(1,17),dsvtsks(1,18),                                &
         dsvtsks(1,19),dsvtsks(1,20),dsvtsks(1,21),dsvtsks(1,22),    &
         dsvtsks(1,23),dsvtsks(1,24),dsvtsks(1,25),dsvtsks(1,26),    &
         dsvtsks(1,27),dsvtsks(1,28),dsvtsks(1,29),dsvtsks(1,30),    &
         dsvtsks(1,31),dsvtsks(1,32),dsvtsks(1,33),dsvtsks(1,34),    &
         dsvtsks(1,35),        &
         dsvtsks(1,36),dsvtsks(1,37),dsvtsks(1,38),dsvtsks(1,39),    &
         dsvtsks(1,40),dsvtsks(1,41),dsvtsks(1,42),dsvtsks(1,43),    &
         dsvtsks(1,44),dsvtsks(1,45),                                &
         dsvtsks(1,46),dsvtsks(1,47),dsvtsks(1,48),dsvtsks(1,49),    & 
!    following are diagnostic variables stored in array svar_surf      
!    variable description is found below.                              
!    msvars in gemini should correspond to current number of variables 
!    in this list.              
         svar_tsks(1,1),svar_tsks(1,2),svar_tsks(1,3),               &
         svar_tsks(1,4),       &
         svar_tsks(1,5),svar_tsks(1,6),svar_tsks(1,7),               &
         svar_tsks(1,8),       &
         svar_tsks(1,9),svar_tsks(1,10),svar_tsks(1,11),             &
         svar_tsks(1,12),      &
         svar_tsks(1,13),&
         svar_tsks(1,15),            &
         svar_tsks(1,17),svar_tsks(1,18),svar_tsks(1,19),            &
         svar_tsks(1,21),&
         svar_tsks(1,23),            &
         svar_tsks(1,24),      &
         svar_tsks(1,25),svar_tsks(1,26),svar_tsks(1,27),            &
         svar_tsks(1,28),      &
         svar_tsks(1,29),svar_tsks(1,30),svar_tsks(1,31),            &
         svar_tsks(1,32),      &
         svar_tsks(1,33),svar_tsks(1,34),svar_tsks(1,35),            &
         svar_tsks(1,36),      &
         svar_tsks(1,37),svar_tsks(1,38),svar_tsks(1,39),            &
         svar_tsks(1,40),      &
         svar_tsks(1,41),svar_tsks(1,42),svar_tsks(1,43),            &
         svar_tsks(1,44),      &
         svar_tsks(1,45),svar_tsks(1,46),svar_tsks(1,47),            &
         svar_tsks(1,48),      &
         svar_tsks(1,49),svar_tsks(1,50),svar_tsks(1,51),            &
         svar_tsks(1,52),      &
         svar_tsks(1,53),svar_tsks(1,54),svar_tsks(1,55),            &
         svar_tsks(1,56),      &
         svar_tsks(1,57),svar_tsks(1,58),svar_tsks(1,59),            &
         svar_tsks(1,60),      &
         svar_tsks(1,61),svar_tsks(1,62),svar_tsks(1,63),            &
         svar_tsks(1,64),      &
         svar_tsks(1,65),svar_tsks(1,66),svar_tsks(1,67),            &
         svar_tsks(1,68),      &
         svar_tsks(1,69),svar_tsks(1,70),svar_tsks(1,71),            &
         svar_tsks(1,72),      &
         svar_tsks(1,73),svar_tsks(1,74),svar_tsks(1,75),            &
         svar_tsks(1,76),      &
         svar_tsks(1,77),svar_tsks(1,78),svar_tsks(1,79),            &
         svar_tsks(1,80),      &
         svar_tsks(1,81),svar_tsks(1,82),svar_tsks(1,83),            &
         svar_tsks(1,84),      &
         svar_tsks(1,85),svar_tsks(1,86),svar_tsks(1,87),            &
         svar_tsks(1,88),      &
         svar_tsks(1,89),svar_tsks(1,90),svar_tsks(1,91),            &
         svar_tsks(1,92),      &
         svar_tsks(1,93),svar_tsks(1,94),svar_tsks(1,95),            &
         svar_tsks(1,96),      &
         svar_tsks(1,97),svar_tsks(1,98),svar_tsks(1,99),            &
         svar_tsks(1,100),     &
         svar_tsks(1,101),svar_tsks(1,102),svar_tsks(1,103),         &
         svar_tsks(1,104),svar_tsks(1,105),svar_tsks(1,106),         &
         svar_tsks(1,107),svar_tsks(1,108),svar_tsks(1,109),         &
         svar_tsks(1,110),     &
         svar_tsks(1,111),svar_tsks(1,112),svar_tsks(1,113),         &
         svar_tsks(1,114),svar_tsks(1,115),                          &
         svar_tsks(1,116),svar_tsks(1,117),                          &
         svar_tsks(1,118),svar_tsks(1,119),                          &
         svar_tsks(1,120),                                           &
         lake_types,lake_no_prog,lake_no_tend,lake_no_diag,          &
         frac_l_tsks(1,1),depth_l_tsks(1,1),                         &
         prog_l_tsks(1,1,1),tend_l_tsks(1,1,1),diag_l_tsks(1,1,1),   &
         tsclim_years(kstart), &
         lcounttice(kstart),icethick(kstart),                        &
         guessid(kstart),                                            & 
                                !    following is physiographic information from ecoclimap stored in a rray eco
                                !    variable description is found below.                              
                                !    meco in gemini should correspond to current number of variables   
                                !    in this list.              
         eco_tsks(1,1),eco_tsks(1,2),eco_tsks(1,3),eco_tsks(1,4),    &
         eco_tsks(1,5),eco_tsks(1,6),eco_tsks(1,7),eco_tsks(1,8),    &
         eco_tsks(1,9),eco_tsks(1,10),eco_tsks(1,11),                &
         eco_tsks(1,12),       &
         eco_tsks(1,13),eco_tsks(1,14),eco_tsks(1,15),               &
         eco_tsks(1,16),       &
         eco_tsks(1,17),eco_tsks(1,18),eco_tsks(1,19),               &
         eco_tsks(1,20),       &
         eco_tsks(1,21),eco_tsks(1,22),eco_tsks(1,23),               &
         eco_tsks(1,24),       &
         eco_tsks(1,25),eco_tsks(1,26),eco_tsks(1,27),               &
         eco_tsks(1,28),       &
         eco_tsks(1,29),eco_tsks(1,30),eco_tsks(1,31),               &
         eco_tsks(1,32),       &
         eco_tsks(1,33),eco_tsks(1,34),eco_tsks(1,35),               &
         eco_tsks(1,36),       &
         eco_tsks(1,37),eco_tsks(1,38),eco_tsks(1,39),               &
         eco_tsks(1,40),       &
         eco_tsks(1,41),eco_tsks(1,42),eco_tsks(1,43),               &
         ahalf,bhalf,hybf,hybh,&
         phcoef%emc,phcoef%csusa,phcoef%c1er,phcoef%cml2,phcoef%fabso3,phcoef%cadd,   &
         dttsk,dqtsk,dcwtsk,dutsk,dvtsk,                             &
         totctsk,cuctsk,dtphtsk,dqphtsk,dcphtsk,                     &
         svtsk,dsvtsk,nsvar,   &
         zprtsk, zttsk,ztrtsk,ztttsk, zqtsk, zqttsk, zctsk,zcttsk,   &
         slwr(kstart),sswr(kstart),slwdn(kstart),sswdn(kstart),      &
         tlwr(kstart),tswr(kstart),tswdn(kstart),                    &
         slwrsea(kstart),sswrsea(kstart),sswrice(kstart),            &
         slwrice(kstart),pblh(kstart),                               &
         omega_tsk,dtdt_kftsk,dqdt_kftsk,dqcdt_kftsk,                &
         div_tsk,oldcov_tsk,zvarcu_tsk,cwcu_tsk ,                    &
         raincv_kf(kstart),snowcv_kf(kstart),umfb(kstart),           &
         shal_cgj(kstart),zdx_kf(kstart),zdy_kf(kstart),             &
         nca(kstart),kf_ind(kstart),                                 &
         kf_base(kstart),kf_top(kstart),afull,bfull,                 &
         ztropo(kstart),ttropo(kstart))                       

    !    gjkf  accumulations of precip modified slightly as conv & ls prec 
    !    gjkf  are in their component forms, looks here that sacdg2(k,1)   
    !    gjkf  is total convective precip (both phases) and                
    !    gjkf  sacdg2(k,2) is total ls precip (both phases).        

    if (lnewpht) then           
       dtacc=dtime*conacc       

       if(nsacdg2>= 1)then    
          do k=kstart,kstop     
             !sacdg2(k,1)=sacdg2(k,1)+dtacc*zprc(k-kstart +1)
             sacdg2(k,1)=zprc(k-kstart+1)+zsnoc(k-kstart+1)           
             accprc(k) = accprc(k)+dtacc*(zprc(k-kstart+1)           &
                  +zsnoc(k-kstart+1))                                 
          enddo
       endif
       if(nsacdg2>= 2)then    
          do k=kstart,kstop     
             !sacdg2(k,2)=sacdg2(k,2)+dtacc*zprl(k-kstart+1)  
             !sacdg2(k,2)=sacdg2(k,2)+dtacc*(zprl(k-kstar t+1)
             !+zsnol(k-ksta rt+1))
             sacdg2(k,2)=zprl(k-kstart+1)+zsnol(k-kstart+1)           
             accprl(k) = accprl(k)+dtacc*(zprl(k-kstart+1)           &
                  +zsnol(k-kstart+1))                                 
          enddo
       endif
       if(nsacdg2>= 3)then    
          do k=kstart,kstop     
             !sacdg2(k,3)=sacdg2(k,3)+dtacc*zsnoc(k-kstar t+1)
             sacdg2(k,3)=zsnoc(k-kstart+1)                            
             accsnoc(k) = accsnoc(k)+dtacc*zsnoc(k-kstart+1)          
          enddo
       endif
       if(nsacdg2>= 4)then    
          do k=kstart,kstop     
             !sacdg2(k,4)=sacdg2(k,4)+dtacc*zsnol(k-kstar t+1)
             sacdg2(k,4)=zsnol(k-kstart+1)                            
             accsnow(k) = accsnow(k)+dtacc*(zsnoc(k-kstart+1)        &
                  +zsnol(k-kstart+1))                                 
             accsnol(k) = accsnol(k)+dtacc*zsnol(k-kstart+1)          
          enddo
       endif
       if(nsacdg2>=10)then    
          do k=kstart,kstop     
             !sacdg2(k,10)=sacdg2(k,10)+dtacc*zslnet(k-kstart+1 )
             !sacdg2(k,10)=sacdg2(k,10)+dtacc*slwr(k)     
             sacdg2(k,10)=slwr(k)                                     
             !slwr(k)=zslnet(k-kstart+1)  
          enddo
       endif
       if(nsacdg2>=11)then    
          do k=kstart,kstop     
             !sacdg2(k,11)=sacdg2(k,11)+dtacc*zssnet(k-kstart+1 )
             !sacdg2(k,11)=sacdg2(k,11)+dtacc*sswr(k)     
             sacdg2(k,11)=sswr(k)                                     
             !sswr(k)=zssnet(k-kstart+1)  
          enddo
       endif
       if(nsacdg2>=12)then    
          do k=kstart,kstop     
             !sacdg2(k,12)=sacdg2(k,12)+dtacc*zsldn(k-kstart+1) 
             !sacdg2(k,12)=sacdg2(k,12)+dtacc*slwdn(k)    
             sacdg2(k,12)=slwdn(k)                                    
          enddo
       endif
       if(nsacdg2>=13)then    
          do k=kstart,kstop     
             !sacdg2(k,13)=sacdg2(k,13)+dtacc*zssdn(k-kstart+1) 
             !sacdg2(k,13)=sacdg2(k,13)+dtacc*sswdn(k)    
             sacdg2(k,13)=sswdn(k)                                    
          enddo
       endif
       if(nsacdg2>=14)then    
          do k=kstart,kstop     
             !sacdg2(k,14)=sacdg2(k,14)+dtacc*ztlnet(k-kstart+1 )
             !sacdg2(k,14)=sacdg2(k,14)+dtacc*tlwr(k)     
             sacdg2(k,14)=tlwr(k)                                     
             !tlwr(k)=ztlnet(k-kstart+1)  
          enddo
       endif
       if(nsacdg2>=15)then    
          do k=kstart,kstop     
             !sacdg2(k,15)=sacdg2(k,15)+dtacc*ztsnet(k-kstart+1 )
             !sacdg2(k,15)=sacdg2(k,15)+dtacc*tswr(k)     
             sacdg2(k,15)=tswr(k)                                     
             !tswr(k)=ztsnet(k-kstart+1)  
          enddo
       endif
       if(nsacdg2>=16)then    
          do k=kstart,kstop     
             !sacdg2(k,16)=sacdg2(k,16)+dtacc*ztsdn(k-kstart+1) 
             !sacdg2(k,16)=sacdg2(k,16)+dtacc*tswdn(k)    
             sacdg2(k,16)=tswdn(k)                                    
          enddo
       endif

       do klev=1,nlev           
          
          if(nsacdg>= 2)then  
             do k=kstart,kstop  
                sacdg(k,klev,2)=sacdg(k,klev,2)+                     &
                     dtacc*zttsk(k-kstart+1,klev)                     
             enddo
          endif
          if(nsacdg>= 3)then  
             do k=kstart,kstop  
                sacdg(k,klev,3)=sacdg(k,klev,3)+                     &
                     dtacc*ztrtsk(k-kstart+1,klev)                    
             enddo
          endif
          if(nsacdg>= 4)then  
             do k=kstart,kstop  
                sacdg(k,klev,4)=sacdg(k,klev,4)+                     &
                     dtacc*ztttsk(k-kstart+1,klev)                    
             enddo
          endif
          if(nsacdg>= 5)then  
             do k=kstart,kstop  
                sacdg(k,klev,5)=sacdg(k,klev,5)+                     &
                     dtacc*zqtsk(k-kstart+1,klev)                     
             enddo
          endif
          if(nsacdg>= 6)then  
             do k=kstart,kstop  
                sacdg(k,klev,6)=sacdg(k,klev,6)+                     &
                     dtacc*zqttsk(k-kstart+1,klev)                    
             enddo
          endif
          if(nsacdg>= 7)then  
             do k=kstart,kstop  
                sacdg(k,klev,7)=sacdg(k,klev,7)+                     &
                     dtacc*zctsk(k-kstart+1,klev)                     
             enddo
          endif
          if(nsacdg>= 8)then  
             do k=kstart,kstop  
                sacdg(k,klev,8)=sacdg(k,klev,8)+                     &
                     dtacc*zcttsk(k-kstart+1,klev)
             enddo
          endif
       enddo                  ! end klev-loop
    endif

    do 410 klev=1,nlev
       do 400 k=kstart,kstop
          dtdt(k,klev) = dttsk(k-kstart+1,klev)
          dqdt(k,klev) = dqtsk(k-kstart+1,klev)
          dudt(k,klev) = dutsk(k-kstart+1,klev)
          dvdt(k,klev) = dvtsk(k-kstart+1,klev)
          dcwdt(k,klev)  = dcwtsk(k-kstart+1,klev)
          totcov(k,klev) = totctsk(k-kstart+1,klev)
          cucov(k,klev)  = cuctsk(k-kstart+1,klev)

          dtdt_kf(k,klev) = dtdt_kftsk(k-kstart+1,klev)
          dqdt_kf(k,klev) = dqdt_kftsk(k-kstart+1,klev)
          dqcdt_kf(k,klev) = dqcdt_kftsk(k-kstart+1,klev)

          oldcov(k,klev) = oldcov_tsk(k-kstart+1,klev)
          zvarcu(k,klev) = zvarcu_tsk(k-kstart+1,klev)
          cwcu(k,klev) = cwcu_tsk(k-kstart+1,klev)
400    enddo
410 enddo

    do jc=1,lake_types
       do k=kstart,kstop
          do l = 1,lake_no_prog
             prog_lakes(k,jc,l) = prog_l_tsks(k-kstart+1,jc,l) 
          enddo
          do l = 1,lake_no_tend
             tend_lakes(k,jc,l) = tend_l_tsks(k-kstart+1,jc,l) 
          enddo
          do l = 1,lake_no_diag
             diag_lakes(k,jc,l) = diag_l_tsks(k-kstart+1,jc,l) 
          enddo
       enddo
    enddo
    do jc=1,msvars
       do k=kstart,kstop
          svar_surf(k,jc) = svar_tsks(k-kstart+1,jc)
       enddo
    enddo

    do jc=1,nsvars
       do k=kstart,kstop
          dsvarsdt(k,jc) = dsvtsks(k-kstart+1,jc)
       enddo
    enddo

    do l=1,nsvar
       do klev=1,nlev
          do k=kstart,kstop
             dsvardt(k,klev,l) = dsvtsk(k-kstart+1,klev,l)
          enddo
       enddo
    enddo

    if( lnewpht .and. (.not.lsamedt) ) then
       do 416 klev=1,nlev
          do 414 k=kstart,kstop
             dtdtph(k,klev) = dtphtsk(k-kstart+1,klev)
             dqdtph(k,klev) = dqphtsk(k-kstart+1,klev)
             dcdtph(k,klev) = dcphtsk(k-kstart+1,klev)
414       enddo
416    enddo

    endif
    return
  end subroutine phtask



  subroutine phcall(kstep,        &
       time,                 &
       klon, klat, klev,   ksvar,   &
       dtime, dtdyn, dtvdif, dtphys,                               &
       conacc,               &
       timesu,               &
       lpp,ldynvd,                                    &
       linit,                &
       phiz,&
       omf,  &
       tsea, frice,          &
       along, coslat, sinlat,&
       rough,rousea,         &
       drolddt, dsolddt, accsunny,                                 &
       accrunoff,            &
       accrunoffopl,accrunofffor,accrunofflake,                    &
       q2d, &
       accprl, accprc,       &
       slwr,sswr,tlwr,tswr,  &
       slwrsea,sswrsea,      &
       slwrice,sswrice,      &
       accsnow,accsnoc,accsnol,                                    &
       slwdn,sswdn,tswdn,    &
       svars, dsvarsdt, nsvars,                                    &
       icethick,             &
       guessid,                                                    &
       lcounttice,           &
       dphidt,&
       totcov,  cucov, cov2d,  cwpath,                             &
       dtdtph, dqdtph, dcdtph,                                     &
       pblh,                 &
       msvars,svar_surf,     &
       eco,meco,             &
       lakes,&
       tsclim_years,         &
       sacdg, nsacdg, sacdg2, nsacdg2,                             &
       dtdt_kf,              &
       dqdt_kf,              &
       dqcdt_kf,             &
       oldcov,               &
       zvarcu,               &
       cwcu,                 &
       div_kf,               &
       om_rmean,             &
       kf_ind,               &
       kf_base,              &
       kf_top,               &
       raincv_kf,            &
       snowcv_kf,            &
       umfb,                 &
       shal_cgj,             &
       nca,                  &
       zdx_kf,               &
       zdy_kf)


    use decomp                  
    use radiation               
    use timers
    use condensmod,only:suphec 
    use flake,only:laketype
    use derived_types
    implicit none    
    type(atm),intent(inout)::dphidt
    type(atm),intent(in)::phiz
    type(laketype),intent(inout)::lakes
    integer,intent(in):: kstep 

    integer,intent(in):: klon,klat,klev      

    real(kind=realkind),intent(in):: time,dtime,dtdyn,dtvdif,dtphys,conacc,timesu               

    logical,intent(in):: lpp,ldynvd    
    logical,intent(in):: linit            
    real(kind=realkind)::omf(klon,klat,klev),                  &
         tsea(klon,klat),      frice(klon,klat),                     &
         along(klon,klat),     &
         coslat(klon,klat),     sinlat(klon,klat),                   &
         rousea(klon,klat), rough(klon,klat),                        &
         pblh(klon,klat)        
    real(kind=realkind) drolddt(klon,klat), dsolddt(klon,klat),                     &
         accsunny(klon,klat), dtphysh,                               &
         accrunoff(klon,klat), &
         accrunoffopl(klon,klat),accrunofffor(klon,klat),            &
         accrunofflake(klon,klat),                                   &
         q2d(klon,klat),       &
         accprl(klon,klat), accprc(klon,klat)                              

    real(kind=realkind) slwr(klon,klat), sswr(klon,klat),                           &
         tlwr(klon,klat), tswr(klon,klat)                             
    real(kind=realkind) accsnow(klon,klat)     
    real(kind=realkind) accsnoc(klon,klat)     
    real(kind=realkind) accsnol(klon,klat)     
    real(kind=realkind) slwrsea(klon,klat),sswrsea(klon,klat),                      &
         slwrice(klon,klat),sswrice(klon,klat)                        
    real(kind=realkind) slwdn(klon,klat),sswdn(klon,klat),tswdn(klon,klat)           
    save dtphysh                     

!    real(kind=realkind):: dtdt(klon,klat,klev),  dqdt(klon,klat,klev), &
!         dcwdt(klon,klat,klev),&
!         dudt(klon,klat,klev),  dvdt(klon,klat,klev)                

    real(kind=realkind):: totcov(klon,klat,klev), cucov(klon,klat,klev),              &
         cov2d(klon,klat),     cwpath(klon,klat),                    &
         dpsdt(klon,klat)            

    real(kind=realkind) dtdtph(klon,klat,klev),dqdtph(klon,klat,klev),              &
         dcdtph(klon,klat,klev)      

!    real(kind=realkind) svar(klon,klat,klev,*),dsvardt(klon,klat,klev,*)             
    integer ksvar               
    integer nsacdg,nsacdg2      
    real(kind=realkind) sacdg(klon,klat,klev,nsacdg)                                 
    real(kind=realkind) sacdg2(klon,klat,nsacdg2)                                    
    integer  nsvars             
    real(kind=realkind):: svars(klon,klat,*), dsvarsdt(klon,klat,*),                  &
         icethick(klon,klat)
    real(kind=realkind),intent(in):: lcounttice(klon,klat)  
    integer:: guessid(klon,klat)
    integer::  msvars             
    real(kind=realkind)::   svar_surf(klon,klat,msvars)                                     

    integer::  meco               
    real(kind=realkind)::     eco(klon,klat,meco)     

    real(kind=realkind)::  tsclim_years(klon,klat)    

    integer,allocatable,dimension(:)::kstart,kstop


    logical lsamedt,lnewpht     
    integer jx,jy,kpar,jlev,k,ktaskind,ntaskx,ktaskx               
    real(kind=realkind) conacp                      

    !    de-staggered dynamic wind component tendencies                         

    real(kind=realkind) duds(klon,klat,klev), dvds(klon,klat,klev)                        

    !    local (saved) adjusted time-stepping parameters                        

    integer nsphys            ! Number of dyn. time steps per phys    
    integer nsvdif            ! Number of vert. diff. time steps p    
    real(kind=realkind)    tsuadj              
    save    nsphys, nsvdif, tsuadj   

    integer ilatl,ilonl,ihorph,ilatf,ilonf                                 

    real(kind=realkind) dtdt_kf(klon,klat,klev),                                    &
         dqdt_kf(klon,klat,klev),                                    &
         dqcdt_kf(klon,klat,klev),                                   &
         oldcov(klon,klat,klev),                                     &
         zvarcu(klon,klat,klev),                                     &
         cwcu(klon,klat,klev),                                       &
         div_kf(klon,klat,klev),                                     &
         om_rmean(klon,klat,klev),                                   &
         raincv_kf(klon,klat), &
         snowcv_kf(klon,klat), &
         umfb(klon,klat),      &
         zdx_kf(klon,klat),    &
         zdy_kf(klon,klat),    &
         ztropo(klon,klat),    &
         ttropo(klon,klat)
    integer nca(klon,klat),    &
         kf_ind(klon,klat),    &
         kf_base(klon,klat),   &
         kf_top(klon,klat),    &
         shal_cgj(klon,klat)         
    integer::idum

#ifdef TIMING                     
    call timer(1,'DESTAG')     
#endif                                 


    if(.not.namelist_is_read)then
       call read_namprc()!read nhorph and lserial and noption
       namelist_is_read = .true.
       if(lserial)then
          ntask = 1
          nhorph = klon*klat
       else
          idum = klon/nhorph
          if(idum*nhorph < klon) idum = idum + 1
          ntask  = klat*idum
       endif
    endif

    allocate(kstart(ntask),kstop(ntask))

    if(lserial) then            
       kstart(1) = 1
       kstop(1) = klon*klat
       ihorph = nhorph          
    else                        
       !    Number of rows/columns excluded from calculations                 
       if(atbase) then          
          ilatf = 0             
       else                     
          ilatf = 1             
       endif
       if(attop) then           
          ilatl = 0             
       else                     
          ilatl = 1             
       endif
       if(atleft) then          
          ilonf = 0             
       else                     
          ilonf = 1             
       endif
       if(atright) then         
          ilonl = 0             
       else                     
          ilonl = 1             
       endif

       !    No of tasks in x-direction and total                              

       ntaskx = (klon-ilonf-ilonl)/nhorph                             
       if(ntaskx*nhorph <klon-ilonf-ilonl) ntaskx = ntaskx + 1     
       ntask  = (klat-ilatf-ilatl)*ntaskx                             

       !    Start- and end-points       

       ktaskind = 0                
       ihorph = 0               
       do jy=1+ilatf,klat-ilatl 
          do ktaskx=1,ntaskx    
             ktaskind = ktaskind + 1  
             kstart(ktaskind)=(jy-1)*klon + (ktaskx-1)*nhorph + ilonf +  1
             if( ktaskx==ntaskx ) then                              
                kstop(ktaskind) = jy*klon - ilonl                        
             else               
                kstop(ktaskind) = kstart(ktaskind) + nhorph - 1             
             end if
             ihorph = MAX(ihorph,kstop(ktaskind)-kstart(ktaskind)+1)        
          enddo
       enddo

       if(ktaskind/=ntask) then  
          write(6,*)' KTASK#NTASK ',ktaskind,ntask                       
          call stop_program(' KTASK#NTASK ')                  
       endif
    endif

    if( kstep==0 .or. linit )then                                  
       !    Initialize some coefficients in physics.                          
       call iniphys(klev,RCAdomain%hybk,RCAdomain%hybi,&
            current_time%year,current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec)

       !    constants for ecmwf physic (mass-flux)                            
       !    ( noption(2) = 2 in this case )                                   
       if ( noption(2)==2 ) call suphec(klev,RCAdomain%ahyb,RCAdomain%bhyb)


       accprl=0._realkind          
       accprc=0._realkind          

       accsnow=0._realkind         
       accsnoc=0._realkind         
       accsnol=0._realkind         

       accrunoff=0._realkind       
       accrunoffopl=0._realkind    
       accrunofffor=0._realkind    
       accrunofflake=0._realkind   

       drolddt=0._realkind         
       dsolddt=0._realkind         
       accsunny=0._realkind        
       q2d=0._realkind             
       totcov = 0._realkind
       cucov  = 0._realkind

       dtphysh=dtphys/3600._realkind          

       !    adjust and save local time-stepping param. at time-step 0              

       if( dtphys<dtdyn ) then                                     
          nsphys = 1            
       else                     
          nsphys = nint(dtphys / dtdyn)                               
       endif

       if( dtvdif>dtdyn ) then                                     
          nsvdif = 1            
       else                     
          nsvdif = nint(dtdyn / dtvdif)                               
       endif

       tsuadj = real(nint(timesu/(real(nsphys,realkind)*dtdyn)),realkind) * (real(nsphys,realkind)*dtdyn) - 0.5_realkind         

       !    Update sun elevation Uandrae
    elseif(current_time%hour==00 .and. current_time%min==00 .and. &
         abs(real(current_time%sec,realkind))<dtime/2._realkind) then 
       call inirad(current_time%year, current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec,                          &
            klev,RCAdomain%hybk,RCAdomain%hybi, &
            phcoef%emc,&
            phcoef%csusa,&
            phcoef%fabso3,&
            phcoef%cadd)

       call partly_solarupdate(current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec)
    endif
    if(linit) then              
       call inirad(current_time%year,current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec,    &
            klev,RCAdomain%hybk,RCAdomain%hybi,&
            phcoef%emc,&
            phcoef%csusa,&
            phcoef%fabso3,&
            phcoef%cadd)

       call partly_solarupdate(current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec)
    endif

    !    select physics work for this particular time-step                 

    if( time <= tsuadj .or. nsphys == 1 ) then                    
       lsamedt = .true.         
       lnewpht = .true.         
       conacp = conacc          
    else                        
       lsamedt = .false.        
       if( mod(kstep, nsphys) == 0 ) then                           
          lnewpht = .true.      
          conacp  = conacc*real(nsphys,realkind)                                     
       else                     
          lnewpht = .false.     
       endif
    endif


    !    save dynamic destaggered wind component tendencies                

    do jlev=1,klev          
       do jy=1,klat
          do jx=1,klon
             if(jx==1) then   
                duds(jx,jy,jlev) = dphidt%u(jx,jy,jlev)                        
             else               
                duds(jx,jy,jlev) =0.5_realkind*(dphidt%u(jx,jy,jlev)+dphidt%u(jx-1,jy,jlev) )
             endif
             if(jy==1) then   
                dvds(jx,jy,jlev) =  dphidt%v(jx,jy,jlev)                        
             else               
                dvds(jx,jy,jlev) = 0.5_realkind*( dphidt%v(jx,jy,jlev) + dphidt%v(jx,jy-1,jlev))
             endif
          enddo
       enddo
    enddo

#ifdef TIMING                     
    call timer(2,'DESTAG')     
#endif                            
#ifdef TIMING                     
    call timer(1,'PHTASK')     
#endif                            
    !    loop over all tasks (to be multi-tasked)                          



    do 1000 ktaskind=1,ntask       
       call phtask(kstep,current_time%year,current_time%month,current_time%day,&
            current_time%hour,current_time%min,current_time%sec,                &
            ihorph,klon,klat,klev,kstart(ktaskind),kstop(ktaskind),        &
            dtime,nsvdif,conacp,RCAdomain%dlat, &
            lpp,lsamedt,lnewpht,ldynvd,                        &
            RCAdomain%ahyb,RCAdomain%bhyb,RCAdomain%hybk,RCAdomain%hybi,                                   &
            phiz%t,  phiz%q,  phiz%cw,  phiz%u,  phiz%v,  omf,                                &
            phiz%ps,                &
            tsea, frice,       &
            along, coslat, sinlat, rough, rousea,                    &
            drolddt, dsolddt, accsunny, dtphysh,                     &
            accrunoff,         &
            accrunoffopl,accrunofffor,accrunofflake,                 &
            q2d, RCAdomain%orogsigm,     &
            accprl, accprc,    &
            dphidt%t, dphidt%q, dphidt%cw, duds, dvds,                           &
            totcov, cucov, cov2d, cwpath,                            &
            dpsdt,             &
            dtdtph, dqdtph, dcdtph,                                  &
            pblh,              &
            slwr,sswr,tlwr,tswr,                                     &
            slwrsea,sswrsea,slwrice,sswrice,                         &
            accsnow,accsnoc,accsnol,                                 &
            slwdn,sswdn,tswdn, &
            svars, dsvarsdt, nsvars,                                 &
            icethick,          &
            guessid,                                                 &
            lcounttice,        &
            msvars,svar_surf,  &
            eco,meco,          &
            lakes%lake_types,lakes%lake_no_prog,                                 &
            lakes%lake_no_tend,lakes%lake_no_diag,                               &
            lakes%frac_lakes,lakes%depth_lakes,                                  &
            lakes%prog_lakes,lakes%tend_lakes,lakes%diag_lakes,                        &
            tsclim_years,      &
            phiz%svar,  dphidt%svar, ksvar,                                    &
            sacdg, nsacdg, sacdg2, nsacdg2,                          &
            dtdt_kf,           &
            dqdt_kf,           &
            dqcdt_kf,          &
            oldcov,            &
            zvarcu,            &
            cwcu,              &
            div_kf,            &
            om_rmean,          &
            kf_ind,            &
            kf_base,           &
            kf_top,            &
            raincv_kf,         &
            snowcv_kf,         &
            umfb,              &
            shal_cgj,          &
            nca,               &
            zdx_kf,            &
            zdy_kf,            &
            RCAdomain%afull,             &
            RCAdomain%bfull,ztropo,ttropo)                  

1000 enddo



#ifdef DEBUG                      
    if(mype==0) write(6,*) ' tsfor5 after phtask ',svar_surf(1,1,78) 
#endif                                 

#ifdef TIMING                     
    call timer(2,'PHTASK')     
#endif                                 

    !    stagger physical wind component tendencies                             

#ifdef TIMING                     
    call timer(1,'DUVSTAG')    
#endif                                 

    !call swap(duds,dvds,klon,klat,klev)
    call swap2(duds,dvds,klon,klat,klev)
    call duvstg(klon,klat,klev, dphidt%u, dphidt%v,duds,dvds)   


#ifdef TIMING                     
    call timer(2,'DUVSTAG')    
#endif                                 

    deallocate(kstart,kstop)
    return

  end subroutine phcall



  subroutine duvstg(klon,klat,klev,dudt,dvdt,duds,dvds)
    use decomp
    implicit none
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(inout)::dudt(klon,klat,klev),dvdt(klon,klat,klev)
    real(kind=realkind),intent(in)::duds(klon,klat,klev),dvds(klon,klat,klev)
    integer::i,j,k
    real(kind=realkind)::du(klon,klat)
    real(kind=realkind)::dv(klon,klat)

    do k=1,klev
       !fix du
       do j=1,klat
          du(1,j) = duds(1,j,k) - dudt(1,j,k)
          do i=2,klon
             du(i,j) = duds(i,j,k)-0.5_realkind*(dudt(i,j,k)+dudt(i-1,j,k))
          enddo
       enddo

       !fix dudt
       do j=1,klat
          do i=1,klon-1
             dudt(i,j,k) = dudt(i,j,k)+0.5_realkind*(du(i,j)+du(i+1,j))
          enddo
          dudt(klon,j,k) = dudt(klon,j,k) + du(i,j)
       enddo


       !fix dv
       do i=1,klon
          dv(i,1) = dvds(i,1,k) - dvdt(i,1,k)
       enddo
       do j=2,klat
          do i=1,klon
             dv(i,j) = dvds(i,j,k) - 0.5_realkind*(dvdt(i,j,k)+dvdt(i,j-1,k))
          enddo
       enddo


       !fix dvdt
       do j=1,klat-1
          do i=1,klon
             dvdt(i,j,k) = dvdt(i,j,k) + 0.5_realkind*(dv(i,j)+dv(i,j+1))
          enddo
       enddo
       do i=1,klon
          dvdt(i,klat,k) = dvdt(i,klat,k) + dv(i,klat)
       enddo

    enddo

    !call swap(dudt,dvdt,klon,klat,klev)
    call swap2(dudt,dvdt,klon,klat,klev)

    return
  end subroutine duvstg



  subroutine iniphys(nlev, &
       hybf,hybh,yearc,monthc,dayc,hourc,minc,secc)

    use kflut
    use kuomod
    use radiation
    use comdfb
    use surface
    use condsmod
    use ccons
    implicit none
    !     
    integer nlev
    integer,intent(in)::yearc, monthc,dayc,hourc,minc,secc
    real(kind=realkind)::hybf(nlev),hybh(nlev+1)
    !     initialize physical parameterizations 
    !     noption(1)=0 means no radiation initialized
    !     noption(1)=1 means that the savijaervi scheme is used.
    
    if(.not.phcoef%isallocated)then
       allocate(phcoef%emc(nlev),phcoef%csusa(nlev),phcoef%c1er(nlev),phcoef%cml2(nlev+1),&
            phcoef%fabso3(nlev),phcoef%cadd(nlev))
       phcoef%isallocated = .true.
    endif

    if( noption(1)==0 ) then
       write(6,*) ' radiation processes has not been initialized'
    endif
    if( noption(1)==1 ) then 
       call inirad(yearc, monthc,dayc,hourc,minc,secc, &
            nlev,hybf,hybh,phcoef%emc,phcoef%csusa,phcoef%fabso3,phcoef%cadd  )
    endif
    !     noption(2)=0 means no condensation processes initialized.
    !     noption(2)=1 means that the straco scheme is used
    !     noption(2)=2 means that the hirlam (tiedtke)  mass flux
    !     scheme is used together with sundqvist 
    !     stratiform condensation.
    !     noption(2)=3 means that the hirlam (tiedtke)  mass flux
    !     scheme is used together with rasch-kristjans-
    !     son condensations scheme
    !     gjkf
    !     noption(2)=4 means that the kain-fritsch convection scheme
    !     is used together with the rasch-kristjans-
    !     son condensations scheme
    !     gjkf              n.b. 
    !     gjkf              presently use of kf+rak presupposes a use of
    !     gjkf              a tke mixing scheme, this is bacause shallow
    !     gjkf              convection in kf requires tke as a closure
    !     gjkf              varaible. kf can be run without shallow
    !     gjkf              convection, but results look worse.
    !     gjkf              in the future an alternative closure for
    !     gjkf              shallow convection can be implemented but
    !     gjkf              pbl tke seems the natural choice.
    !     gjkf              
    !     gjkf

    !     
    if( noption(2)==0 ) then
       write(6,*) ' condensation processes not initialized'
    endif
    !     
    if( noption(2)==1 ) then
       call iniconds
    endif
    !     
    if( noption(2)==2 ) then
       call inicons
    endif
    !     gjkf
    if( noption(2)==4 ) then
       call inikf
       call lutab
       !     write(6,*)'kain-fritsch initialised'
    endif
    !     gjkf
    !     
    !     (the infomation in 'inicond' is currently not used )
    !     
    call inicond(nlev,hybf,phcoef%c1er)
    !     
    !     -----------------------------------------------------------------    
    !     
    !     3) 
    !     
    !     turbulence
    !     noption(3)=0 means no condensation processes initialized.
    !     noption(3)=1 means that the holtslag scheme is used
    !     noption(3)=2 means that the cbr turbulent kinetic energy based
    !     scheme is used.
    !     ---------------------------------------------------------------
    !     
    !     
    !     ps020320      if( (noption(3)==1) ) then
    if( (noption(3)==1) .or. (noption(3)==2)) then
       call inidif2(nlev,hybf,hybh,phcoef%cml2)
    endif
    !     
    if( (noption(3)==0) ) then
       write(6,*) 'warning: turbulence processes not initialized'
    endif
    !     gj   lutab writes look up tables needed in kfcumulus
    !     gj   these variables get wqritten into common kflut
    call lutab
    !     ---------------------------------------------------------------
    !     
    !     surface processes
    !     (currently not initialized at this stage, 
    !     isba to be initialized here using noption(4)=2)
    !     ---------------------------------------------------------------

    call inisurf
    return
  end subroutine iniphys



  subroutine akfrak(nhor, nlev, kstart, kstop, dtime, conacc, &
       frland, frice, ts, ps,   ahyb, bhyb, u, v, t, &
       q, cw, tke, pf, dpf, dcwdin, dtdtin, dqdtin,  &
       cov2d,  prcpst, stsnow, draindt, dsnowdt, dtdt, &
       dqdt, dcwdt, totcov, hcucov, preta, orogsigm, raincv_kf, &
       snowcv_kf, snowh, dph, zdx_kf, zdy_kf, oldcov, omega_pa, div_kf, &
       lesat, pblh, om_cpbl, zfrpbl, zfrtop, zlevtop,  &
       relhum, udr, udrf, zvarcu, gpot, wdzf, &
       kf_ind, shal_cgj, lsb, 	    &
!cgj300611
       cwcu,rhcrit)
    !large scale rh cld from t-1
    !     gjmoistcbr
    !     gjmoistcbr
    !
    !     front end to kain fritsch convection & rak condensation
    use confys  
    use bkf_mod  
    use akfcum_mod  
    use condensmod
    use pcondmod  
    use cldfrcmod
    use ccons  
    implicit none  
    integer, intent(in) ::nhor, nlev, kstart, kstop  
    !     gj   i added for array looping
    integer ::i, jl, jk  
    real(kind=realkind) ::dtime, conacc  
    real(kind=realkind) ::frland(nhor), ts(nhor), ps(nhor)
    real(kind=realkind) ::frice(nhor)  
    real(kind=realkind) ::bhyb(nlev+1), ahyb(nlev+1)  
    real(kind=realkind) ::t(nhor, nlev), q(nhor, nlev), cw(nhor, nlev), &
         dcwdin(nhor, nlev)
    real(kind=realkind) ::dtdtin(nhor, nlev), dqdtin(nhor, nlev), pf(nhor,nlev)
    real(kind=realkind) ::dtdtin_ttt(nhor, nlev)
    real(kind=realkind) ::dpf(nhor, nlev)  
    !     gjshallow
    real(kind=realkind) ::tke(nhor, nlev)  
    !     gjshallow
    real(kind=realkind) ::cov2d(nhor), dsnowdt(nhor), draindt(nhor)

    real(kind=realkind) ::dtdt(nhor, nlev), dqdt(nhor, nlev), dcwdt(nhor, nlev) &
         , totcov(nhor, nlev), hcucov(nhor, nlev), shalcld(nhor, nlev)

    real(kind=realkind) :: prcpst(nhor), stsnow(nhor)  
    real(kind=realkind) :: hldcp2(nhor, nlev), dlnpdt(nhor, nlev)  
    real(kind=realkind) :: omega_kf(nhor, nlev),  div_kf(nhor, nlev),  zdx_kf(nhor), &
         zdy_kf(nhor), u(nhor, nlev),  v(nhor, nlev),  &
         dtdt_kf(nhor, nlev),  dqdt_kf(nhor, nlev),  dqcdt_kf(nhor, nlev),&
         raincv_kf(nhor),  snowcv_kf(nhor),  umf(nhor, nlev),&
         udrf(nhor, nlev),                 &
         udr(nhor,nlev),  snowh(nhor),  dph(nhor, nlev+1), rdph(nhor,nlev), &
         oldcov(nhor, nlev),  rad_cloud(nhor, nlev)
    real(kind=realkind) :: orogsigm(nhor)  
    real(kind=realkind) :: qvap(nhor, nlev), qliq(nhor, nlev), qice(nhor, nlev)  
    real(kind=realkind) :: dz_half(nhor, nlev)  
    !     gj   shallow conv vars
    integer :: lsb(nhor), shal_cgj(nhor), oro(nhor)  
    real(kind=realkind) :: relhum(nhor, nlev)  
    real(kind=realkind) :: lstmp(nhor, nlev), hcutmp(nhor, nlev), shaltmp(nhor,nlev)

    real(kind=realkind) :: gpot(nhor,nlev), wdzf(nhor,nlev), pblh(nhor), rhcrit(nhor,nlev)  
    !     gj  vars for rk cond
    real(kind=realkind) :: preta(nhor, nlev), dt, ro_cgj(nhor, nlev), omega_pa(nhor, nlev)

    integer :: nca(nhor), ncldck, kstep,  k, kf_base(nhor), &
         kf_ind(nhor), kf_top(nhor)
    !     gj
    real(kind=realkind) :: wkf1(nhor), wkf2(nhor), wkf3(nhor)  
    real(kind=realkind) :: zfrtop(nhor)  
    integer :: zfrpbl(nhor), zlevtop(nhor)  

    real(kind=realkind) :: om_cpbl(nhor)  
    real(kind=realkind) :: tlow, fice(nhor, nlev)  

    real(kind=realkind) :: zqcu, zqt, zdqdz, zwrk1, zwrk3, zvarcu1, zvarcu(nhor,nlev),&
         clddep(nhor), wzz
    !     gjmoistcbr
!cgj300611
    real(kind=realkind) :: cwcu(nhor, nlev), rhenv  

    logical :: loxu, lesat  
    loxu = .false.  
    tlow = 250.16_realkind  
    dt = conacc * dtime  
    !     gj
    !     gj   dt is 1 timestep in secs everywhere but nstep=0
    !     gj   dtphys=dt
    !
    !     gj   intialise arrays to zero
    do jk = 1, nlev  
       do jl = kstart, kstop  
          hldcp2(jl, jk) = 0._realkind  
          dlnpdt(jl, jk) = 0._realkind  
          dtdt(jl, jk) = 0._realkind  
          dqdt(jl, jk) = 0._realkind  
          dcwdt(jl, jk) = 0._realkind  
          !     gj1dcode
          lstmp(jl, jk) = 0._realkind  
          hcutmp(jl, jk) = 0._realkind  
          shaltmp(jl, jk) = 0._realkind  
          !     gj1dcode
          !     gj
          ro_cgj(jl, jk) = pf(jl, jk) /(rair *(t(jl, jk) *(1._realkind+0.61_realkind *  q(jl, jk) ) ) )
          omega_kf(jl, jk) = -(omega_pa(jl, jk) /(ro_cgj(jl, jk)  * gravit) )
!cgj300611
          cwcu(jl,jk)=0.	!diagnostic conv cloud water exists over 1 timestep
!cgj300611
       enddo

    enddo

    call ficefun(nhor, nlev, kstart, kstop, t, fice)  
    do i = kstart, kstop  
       !     gj
       !     gj    oro used in kfcumulus for sfc inhomeneity calc
       !     gj    if frland>5% then asuume a land point.
       if(frland(i) <0.05_realkind) oro(i) = 0  
       if(frland(i) >=0.05_realkind) oro(i) = 1  


    enddo
    !     gj   get some prelim values
    call prcond(nhor, nlev, kstart, kstop, t, hldcp2)  
    !     gjkf
    ncldck = nint(600._realkind / dtime)  
    !     gj  check value of nca, if zero a new convective cycle is to
    !     gj  begin. set convective tends at these points to zero, otherwise
    !     gj  retain existing tends to complete nca convective period
    !     gj
    do i = kstart, kstop  
       nca(i) = 0  
       shal_cgj(i) = 0  
       kf_ind(i) = 0  
       kf_base(i) = 0  
       kf_top(i) = 0  
       raincv_kf(i) = 0._realkind  
       snowcv_kf(i) = 0._realkind  
       !     gj
       do k = 1, nlev  
          dtdt_kf(i, k) = 0._realkind  
          dqdt_kf(i, k) = 0._realkind  
          dqcdt_kf(i, k) = 0._realkind  
          !     gj1dcode
          qvap(i, k) = 0._realkind  
          qliq(i, k) = 0._realkind  
          qice(i, k) = 0._realkind  
          !     gj1dcode
          umf(i, k) = 0._realkind  
          udr(i, k) = 0._realkind  
          udrf(i, k) = 0._realkind  
          zvarcu(i, k) = 0._realkind  
       enddo
    enddo
    !     gj
    if(lbkf) then  
       !     gjbkf
       !     gjbkf n.b. on exit from bkfcall, the convective tendencies and out
       !     gjbkf variables such as umf, qvap etc are combined deep and shallo
       !     gjbkf convection values. this combination is done in bkfcall.f90 a
       !     gjbkf will be a one or other combination as if deep convection occ
       !     gjbkf shallow convection is not allowed at the same grid point and
       !     gjbkf timestep. shallow convection is only tested for if deep conv
       !     gjbkf is not supported.
       !     gjbkf
       call bkfcall(nhor, nlev, kstart, kstop, nca, dtime, zdx_kf, &
            zdy_kf, orogsigm, div_kf, zfrtop, zlevtop,             &
            pf, t, q, u, v, omega_kf, dtdt_kf, dqdt_kf,            &
            dqcdt_kf, raincv_kf, snowcv_kf, kf_base, kf_top, umf,  &
            udr, udrf, om_cpbl, &
            qvap, qliq, qice, clddep, shal_cgj,frland,gpot )
    else  
       call akfcum(nhor, nlev, kstart, kstop, nca, dtime, zdx_kf, &
            zdy_kf, orogsigm, pf, t, q, u, v, omega_kf, dtdt_kf, dqdt_kf, &
            dqcdt_kf, raincv_kf, snowcv_kf, kf_base, kf_top, umf, om_cpbl, &
            qvap, qliq, qice, clddep, shal_cgj,(ncldck),(ahyb),(bhyb), &
            (ps),(ts),(tke),(div_kf),(kf_ind),(udr),(zfrpbl), &
            (zfrtop),(zlevtop),(dz_half),(lsb),(wkf1),(wkf2),(wkf3) )


    endif
    do i = kstart, kstop  
       nca(i) = 0  

    enddo
    do k = 1, nlev  
       do i = kstart, kstop  
          !nlev+1 is ps-pfull(nk)..inv not neede
          rdph(i,k) = 1._realkind / dph(i,k)  
          dz_half(i,k) = dz_half(i,k)/gravit  
       enddo

    enddo
    call acldfrc(nhor, nlev, kstart, kstop, pf, rdph, t, q, cw, fice, &
         omega_pa, kf_top, kf_base, lstmp, oldcov, dpf, umf, oro, snowh, &
         frland, frice, hcutmp, ts, ps, udrf, cov2d, shal_cgj, kstep, &
         rad_cloud, zdx_kf, zdy_kf, orogsigm, wdzf, rhcrit, lesat, qvap, &
         qliq, qice, gpot, pblh, shaltmp, relhum)
    !     gj1dcode
    !     gj1dcode
    dtdtin_ttt=dtdtin
    do k = 1, nlev  

       do i = kstart, kstop  
!cgj170711
!cgj170711 Calculate clear sky fraction of grid box relative humidity
!          if(rad_cloud(i,k)<1._realkind)
!           rhenv=(relhum(i,k)-rad_cloud(i,k))/                       &
!                       (1._realkind-rad_cloud(i,k))
!          else
!           rhenv=relhum(i,k)
!          endif
!cgj170711
!cgj170711 Evaporation of detraining convective cloud water can follow 3 options: All cloud water
!cgj170711 evaporates in a timestep as a function of (i) 1-rad_cloud (clear fractio of grid box)
!cgj170711 (ii) 1-rhenv (RH of non-cloud fraction) or (iii) 1-relhum
!cgj170711
          if(shaltmp(i, k) >0._realkind) then  
             dtdtin(i, k) = dtdtin(i, k)+dtdt_kf(i, k) -(hldcp2(i, k) &
                  * dqcdt_kf(i, k) *(1._realkind - relhum(i,k)) )
             dqdtin(i, k) = dqdtin(i, k)+dqdt_kf(i, k)+(dqcdt_kf(i, &
                  k) *(1._realkind - relhum(i,k)) )
             dcwdin(i, k) = dcwdin(i, k)+(dqcdt_kf(i, k)*relhum(i,k))
             cwcu(i,k)=dqcdt_kf(i,k)*relhum(i,k)*dtime
          else  
             dtdtin(i, k) = dtdtin(i, k)+dtdt_kf(i, k) -(hldcp2(i, k) &
                  * dqcdt_kf(i, k) *(1._realkind - relhum(i,k) ) )
             dqdtin(i, k) = dqdtin(i, k)+dqdt_kf(i, k)+(dqcdt_kf(i, &
                  k) *(1._realkind - relhum(i,k)) )
             dcwdin(i, k) = dcwdin(i, k)+(dqcdt_kf(i, k)*relhum(i,k))
             cwcu(i,k)=dqcdt_kf(i,k)*relhum(i,k)*dtime
          endif
!cgj170711
          !     gj
          !     gj     calculate a convective scale variance
          !     gj
          if(k>=2) then  
             wzz =(gpot(i, k - 1) - gpot(i, k) ) / gravit  
             if(k>=kf_top(i) .and.k<=kf_base(i) ) then  
                if(k>kf_base(i) ) then  
                   zqcu = qvap(i, k)  
                   zqt = q(i, k)+cw(i, k)  
                   zdqdz =((q(i, k)+cw(i, k) ) -(q(i, k - 1) &
                        +cw(i, k - 1) ) ) / wzz
                   zdqdz = sqrt(zdqdz * zdqdz)  
                else  
                   zqcu = qvap(i, k)+qliq(i, k)+qice(i, k)  
                   zqt = q(i, k)+cw(i, k)  
                   zdqdz =((q(i, k)+cw(i, k) ) -(q(i, k - 1) &
                        +cw(i, k - 1) ) ) / wzz
                endif
                zwrk1 = umf(i, k) *(zqcu - zqt) * clddep(i) * zdqdz  
                zwrk3 = ro_cgj(i, k) * om_cpbl(i)  
                if(om_cpbl(i) >0._realkind.and.zwrk1>0._realkind) then  
                   zvarcu1 = sqrt(zwrk1 / zwrk3)  
                else  
                   zvarcu1 = 0._realkind  
                endif
                zvarcu(i, k) = zvarcu1  
             endif
          endif
          !     gj
       enddo
    enddo
    !     gj
    !     gj110608  passing in lstmp below enables condensation with rh clou
    !     gj110608  passing in totcov, means cloud fraction from statcld in
    !     gj
!RCA4
    call apcond(nhor, nlev, kstart, kstop, dqdtin, dtdtin, omega_pa, &
!cgj040711         dtime, q, cw, t, pf, dpf, fice, lesat, lstmp, totcov, hcutmp, &
         dtime, q, cw, t, pf, dpf, fice, lstmp, totcov, hcutmp, &
         shaltmp, dtdt, dqdt, dcwdt, prcpst, preta, frland, frice, snowh, &
         ps, dz_half, gpot, pblh)

!RCA3
!    call apcond(nhor, nlev, kstart, kstop, dqdtin, dtdtin, omega_pa, &
!         dtime, q, cw, t, pf, dpf, fice, lstmp, totcov, hcutmp, &
!         shaltmp, dtdt, dqdt, dcwdt, prcpst, preta, frland, frice, snowh, &
!         ps, dz_half, gpot, pblh)


    !     gjstcld    +              totcov,oldcov,hcutmp,shaltmp,
    !     gj
    !     gj  calculate xu-randall cloud fraction using the cw & rh values
    !     gj
    if(loxu) then  
       call axucld(nhor, nlev, kstart, kstop, t, q, cw, pf, ps, &
            rhcrit, lesat, lstmp, totcov, hcutmp, shaltmp, shal_cgj, cov2d)
       !     gjshallowcld
       !     gjshallowcld
    endif
    !     gj
    !     gj   distinguish between startiform and convective precip
    do jl = kstart, kstop  
       if(t(jl, nlev - 1) <tmelt.and.t(jl, nlev) <tmelt.and.ts( &
            jl) <tmelt) then
          stsnow(jl) = prcpst(jl)  
          prcpst(jl) = 0._realkind  
       else  
          stsnow(jl) = 0._realkind  
       endif
       draindt(jl) = draindt(jl)+prcpst(jl)+raincv_kf(jl)  
       dsnowdt(jl) = dsnowdt(jl)+stsnow(jl)+snowcv_kf(jl)  
    enddo
    !     gj
    !     gj   set tends of precip and calc accumulations
    do jk = 1, nlev  
       do jl = kstart, kstop  
          !     gj for statcld from moist tke to be active comment next 2 lines
          if(.not.loxu) totcov(jl, jk) = lstmp(jl, jk)  
          oldcov(jl, jk) = lstmp(jl, jk)  
          hcucov(jl, jk) = hcutmp(jl, jk)  
          shalcld(jl, jk) = shaltmp(jl, jk)  
       enddo
    enddo
    !
    return  
  end subroutine akfrak






end module physpar
