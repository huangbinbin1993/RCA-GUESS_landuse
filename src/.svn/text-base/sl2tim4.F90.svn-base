module sl2tim4
  use decomp
  use comhkp
  use referenceParameters
  use physpar
  use decomp
  use sldynm5
  use domainMod
  use rcadomainmod,only:ffmean,hxv,ra,rdlam,hyu,rdth,plong,pclat,pslat,rhxu,rhyv,fpar
  use mod_diffh
  use mod_implicit_solver
  use timers
  use derived_types
  use boundaryRelaxation,only:npbpts
  use lateralBc
#ifdef USING_NETCDF
    use netcdfrestart
#else
    use restart
#endif
  use calendar
  use confys
  implicit none
  private

  integer,save::khalo=10
  integer,save::nslpqi=1     !number of iterations for calculation of displacement in the semi-lagrangian scheme
  integer,save::nslinc=4    !type of interpolation at the midpoint in the semi-lagrangian scheme (not used)
  integer,save::nslind=4    !type of interpolation at the departure point in the semi-lagrangian scheme =1 linear, =2 quadratic =3 mixed linear/cubic 
  integer,save::nslint(100)=4! list of interpolation types for each iteration i=1,nslpqui =1 linear, =2 quadratic =3 cubic  4mixed linear/cubic
  real(kind=realkind),save::epsg=0.2_realkind    ! coefficient for the gravity wave damper in the semi-lagrangian sckeme
  real(kind=realkind),public,save::epsn=0.2_realkind ! coefficient for the gravity wave damper 
  integer,save::nslext=2 !order of time-extrapolation        
  logical,public,save::nlslan=.true.! .true. if semi-lagrangian advection
  logical,save::nlsl3d = .true.! .true. if 3-dim. semi-lagrangian advection

  logical,save::nlitrh = .false. !true for iterative Helmholtz-solver
  real(kind=realkind),save::aerrih= 1.0E-13_realkind ! abs. error tol. for iterative Helmholtz solver
  real(kind=realkind),save:: rerrih = 1.0E-08_realkind ! rel. error tol. for iterative Helmholtz solver
  integer,save::nityph = 2! iteration type for iterative Helmholtz solve
  integer,save::nptyph = 1! preconditioning type for iterative Helmholtz solver
  logical,save::dynamic_halo=.true.
  real(kind=realkind),save::safety_factor=0.7_realkind

  public sl2tim,omcomp,read_namsl

contains

  subroutine read_namsl()

    namelist /namsl/ nlslan,nslpqi,nslinc,nslind,nslint,nlsl3d,epsg,epsn,nslext,&
         nlitrh,aerrih,rerrih,nityph,nptyph,dynamic_halo,khalo,safety_factor 
    open(57,file='namelists.dat',status='old')
    read(57,nml=namsl)
    close(57)
    if(mype==0)then
       write(6,nml=namsl)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namsl)
       close(1)
    endif

  end subroutine read_namsl

  subroutine sl2tim(nstep,klon,klat,klev,ksvar,RCAdomain,&
        phim,phiz,phip,&
        surf_vars,&
        totcov,cucov,cov2d,cwpath,pblh,&
        sacdg,kacdg,sacdg2,kacdg2,&
        nlinit,drolddt,dsolddt,accsunny,accrunoff,&
        accrunoffopl,accrunofffor,accrunofflake,q2d,&
        accprl,accprc,slwr,sswr,tlwr,tswr,slwrsea,sswrsea,&
        slwrice,sswrice,slwdn,sswdn,tswdn,svarsz,svarsp,ksvars, &
        icethick,guessid,lcounttice,&
        accsnow,accsnoc,accsnol,msvars,svar_surf,eco,meco, &
        tsclim_years, &
        dtdtph, dqdtph, dsdtph, &
        nca,dtdt_kf,dqdt_kf,dqcdt_kf,raincv_kf,snowcv_kf,umfb,&
        shal_cgj,kf_ind,kf_base,kf_top,oldcov,zvarcu,cwcu,lakes)

    !     SL2TIM - SEMI-LAGRANGIAN SEMI-IMPLICIT 2-TIME-LEVEL SCHEME
        
    use flake,only:laketype
    use surface,only:surfVariables
    implicit none
    type(laketype),intent(inout)::lakes
    type(domain),intent(in)::RCAdomain
    integer,intent(in)::klon,klat,klev,ksvar,nstep
    type(surfVariables),intent(inout)::surf_vars

    type(atm),intent(inout)::phim,phip
    type(atm),intent(inout)::phiz 

   
    real(kind=realkind)::     totcov(klon,klat,klev), & ! total cloud cover
         cucov(klon,klat,klev), & ! cumulus cloud cover
         cov2d(klon,klat),  &
         cwpath(klon,klat),  &! vertically integrated cloud water content
         pblh(klon,klat),    &! boundary layer height
         dtdtph(klon,klat,klev), &
         dqdtph(klon,klat,klev),  &
         dsdtph(klon,klat,klev)!,  &

    real(kind=realkind)::accrunoff(klon,klat),                          &
         accrunoffopl(klon,klat),accrunofffor(klon,klat),  &
         accrunofflake(klon,klat),                         &
         drolddt(klon,klat), dsolddt(klon,klat),           &
         accsunny(klon,klat),                              &
         q2d(klon,klat),                                   &
         accprl(klon,klat), accprc(klon,klat)

    real(kind=realkind):: slwr(klon,klat), sswr(klon,klat),tlwr(klon,klat), tswr(klon,klat)
    real(kind=realkind):: slwrsea(klon,klat),sswrsea(klon,klat),slwrice(klon,klat), sswrice(klon,klat)
    real(kind=realkind):: accsnow(klon,klat)
    real(kind=realkind):: accsnoc(klon,klat)
    real(kind=realkind):: accsnol(klon,klat)
    real(kind=realkind):: slwdn(klon,klat),sswdn(klon,klat),tswdn(klon,klat)
    logical::  nlinit
    integer::  ksvars
    real(kind=realkind),dimension(klon,klat,ksvars):: svarsz,svarsp
    real(kind=realkind):: icethick(klon,klat)
    real(kind=realkind),intent(in)::lcounttice(klon,klat)
    integer:: guessid(klon,klat)
    integer::  msvars
    real(kind=realkind)::   svar_surf(klon,klat,msvars)
    integer:: meco
    real(kind=realkind),target:: eco(klon,klat,meco)
    real(kind=realkind)::tsclim_years(klon,klat)
    integer::pretyp(klev)
    integer,intent(in)::kacdg, kacdg2

    real(kind=realkind)::sacdg(klon,klat,klev,kacdg)
    real(kind=realkind)::sacdg2(klon,klat,kacdg2)

    integer::jx,jy,jk,js,icall
    real(kind=realkind)::omfz(klon,klat,klev),omhz(klon,klat,klev+1)         
    real(kind=realkind)::zapp(klon,klat,klev)                                 
    real(kind=realkind)::zapspx(klon+1,klat+1)
    real(kind=realkind)::zahalf(klev+1),zbhalf(klev+1)
    real(kind=realkind)::zdt,zrdt,zconacc,ztimeph

    real(kind=realkind)::dtdt_kf(klon,klat,klev) 
    real(kind=realkind)::dqdt_kf(klon,klat,klev)  
    real(kind=realkind)::dqcdt_kf(klon,klat,klev) 
    real(kind=realkind)::oldcov(klon,klat,klev)   
    real(kind=realkind)::zvarcu(klon,klat,klev)   
!cgj300611
    real(kind=realkind)::cwcu(klon,klat,klev)   
    real(kind=realkind)::om_rmean(klon,klat,klev) 
    real(kind=realkind)::div_kf(klon,klat,klev)

    real(kind=realkind)::raincv_kf(klon,klat) 
    real(kind=realkind)::snowcv_kf(klon,klat) 
    real(kind=realkind)::umfb(klon,klat)      
    real(kind=realkind)::zdx_kf(klon,klat)    
    real(kind=realkind)::zdy_kf(klon,klat)    
    real(kind=realkind)::pkf,rokf,omkf

    integer::nca(klon,klat)  
    integer::kf_ind(klon,klat)   
    integer::kf_top(klon,klat)   
    integer::kf_base(klon,klat)  
    integer::shal_cgj(klon,klat)

    real(kind=realkind):: eps
    real(kind=realkind):: hhdiv(klon,klat,klev)              


    real(kind=realkind),allocatable,dimension(:,:):: apsx, alnpsx
    real(kind=realkind),allocatable,dimension(:,:,:)::tx,ux,vx,qx,sx 


    lakes%frac_lakes => eco(:,:,44:46)
    lakes%depth_lakes => eco(:,:,47:49)


    allocate(apsx(klon,klat), alnpsx(klon,klat))
    allocate(tx(klon,klat,klev))
    allocate(ux(klon,klat,klev))
    allocate(vx(klon,klat,klev))  
    allocate(qx(klon,klat,klev))  
    allocate(sx(klon,klat,klev))

    apsx = -666.0_realkind
    alnpsx = -666.0_realkind
    tx = -666.0_realkind
    ux = -666.0_realkind
    vx = -666.0_realkind
    qx = -666.0_realkind
    sx = -666.0_realkind




    if(ksvar/=0.and..not.nlsl3d) then
       write(6,*)'illegal combination of options:'
       write(6,*)'ksvar=',ksvar
       write(6,*)'nlsl3d=',nlsl3d
       write(6,*)'with the present coding structure this option cannot be implemented'
       call stop_program('')
    endif



    zdt = real(nSeconds(ndtime),realkind) 

    icall = nslext
    if(nstep==0) icall = 1
    if(nslext>=3) then
       call stop_program( 'nslext>=3 not supported')
    endif

    !     extrapolation to time n+1/2 in compfx
    call compfx(klon,klat,klev,icall, &
         phim%ps,phim%lnps,phim%u,phim%v,phim%t,phim%q,phim%cw,phim%edot,     &
         phiz%ps,phiz%lnps,phiz%u,phiz%v,phiz%t,phiz%q,phiz%cw,phiz%edot,     &
         apsx,alnpsx,ux,vx,tx,qx,sx)



    !     compute tendencies at time n+1/2 (or time n) in sldyn


    zapspx = -666.0_realkind
    div_kf =  0.0_realkind
    zapspx(1:klon,1:klat) = apsx !this is consistent across all processors since compfx does a swap
    if(atright) then
       zapspx(klon+1,1:klat) = apsx(klon,1:klat)
    endif
    if(attop) then
       zapspx(1:klon,klat+1) = apsx(1:klon,klat)
    endif
    if(attop.and.atright) then
       zapspx(klon+1,klat+1) = apsx(klon,klat)
    endif
    call swap_ps(zapspx,klon,klat)

    call sldyn(klon , klat , klev ,  &
          ffmean,&
          alnpsx,zapspx, ux, vx, tx, qx, sx,               &
          phip%ps, phip%u, phip%v, phip%t, phip%q, phip%cw,          &!these are used as temporaries (d/dt variables)
          zapp,                                     &
          nltvir, nlsl3d,RCAdomain%fis,div_kf,RCAdomain%ahyb,RCAdomain%bhyb)                                  


    call swap2d(phip%ps,klon,klat)
    call swap6(phip%u,phip%v,phip%t,phip%q,phip%cw,zapp,klon,klat,klev)
    call swap(div_kf,klon,klat,klev)

    !     1:st order Euler explicit time stepping

    phip%ps = phiz%ps+zdt*phip%ps
    phip%u = phiz%u+zdt*phip%u
    phip%v = phiz%v+zdt*phip%v
    phip%t = phiz%t+zdt*phip%t
    phip%q = phiz%q+zdt*phip%q
    phip%cw = phiz%cw+zdt*phip%cw

    phip%svar = phiz%svar
    phip%timestamp = phip%timestamp+ndtime
    svarsp = svarsz

    !     4:th ORDER LINEAR HORIZONTAL DIFFUSION SCHEME
    if (nlhdif) then
       call hdiff4(klon,klat,klev,ksvar,nltcrf,zdt,         &
             phiz%lnps,phiz%t,phiz%u,phiz%v,phiz%q,phiz%cw,phiz%svar, &
             phip%t,phip%u,phip%v,phip%q,phip%cw,phip%svar)

       call swap5(phip%u,phip%v,phip%t,phip%q,phip%cw,klon,klat,klev)
       do js = 1,ksvar
          call swap(phip%svar(:,:,:,js),klon,klat,klev)
       enddo

    endif

    !     SEMI-IMPLICIT CALCULATIONS


    !     compute the vertical velocity at time level n
    call comped(klon,klat,klev,phiz%ps,phiz%u,phiz%v,phiz%edot,RCAdomain%ahyb,RCAdomain%bhyb)

    call swapklevp1(phiz%edot,klon,klat,klev+1)

    !     semi-lagrangian calculations in sldynm


    call sldynm(klon,klat,klev,ksvar,khalo,dynamic_halo,safety_factor,    &
         nslpqi,nslint(1),nslind,     &
         zdt,ffmean,RCAdomain%fis,     &
         phim%ps, phim%u,phim%v, phim%edot, &
         phiz%ps,phiz%lnps,phiz%u,phiz%v,phiz%t,phiz%q,phiz%cw,phiz%edot,phiz%svar, &
         phip%ps,phip%lnps,phip%u,phip%v,phip%t,phip%q,phip%cw,phip%edot,phip%svar, &
         zapp,nlsl3d,epsg,RCAdomain%ahyb,RCAdomain%bhyb)

    call swap2d(phip%lnps,klon,klat)
    call swap2d(phip%ps,klon,klat)
    call swap5(phip%u,phip%v,phip%t,phip%q,phip%cw,klon,klat,klev)
    do js = 1,ksvar
       call swap(phip%svar(:,:,:,js),klon,klat,klev)
    enddo
    !     prevent negativ values of cloud water and TKE
    phip%cw = max(phip%cw,0.0_realkind)
    phip%svar = max(phip%svar,0.0_realkind)

    !     SEMI-LAGRANGIAN EXPLICIT ADJUSTMENT IN SLEXPA

    call slexpa(klon  , klat, klev, npbpts, &
          zdt   , ffmean, RCAdomain%fis,             &
          phim%lnps, phim%u, phim%v, phim%t,               &
          phiz%lnps, phiz%u, phiz%v, phiz%t,               &
          phip%lnps, phip%u, phip%v, phip%t,               &
          hhdiv )

    call swap2d(phip%lnps,klon,klat)
    call swap4(phip%u,phip%t,phip%v,hhdiv,klon,klat,klev)

    !     SOLVE HELMHOLTZ EQUATIONS IN HHSOLV
    if ( nlitrh ) then
       if (nptyph > 1 .and. nptyph < klev) then
          pretyp(1:nptyph) = 1
          pretyp(nptyph+1:klev) = 0
       elseif ( nptyph == 1 ) then
          pretyp(1:klev) = 1
       else
          pretyp(1:klev) = 0
       endif
       call hhsolvitr(klon,klat,klev,npbpts,zdt,hhdiv,nlsimp,epsg,nityph,pretyp,aerrih,rerrih)
    else
       call hhsolv(klon,klat,klev,npbpts,zdt,hhdiv,nlsimp,epsg,RCAdomain%ahyb,RCAdomain%bhyb)
    endif
    call swap(hhdiv ,klon ,klat ,klev)



    !     IMPLICIT ADJUSTMENT IN IMPADJ
    call impadj(klon,klat,klev,npbpts,zdt,ffmean,  &
         phip%lnps,phip%u,phip%v,phip%t,hhdiv)


    call swap2d(phip%lnps,klon ,klat)
    call swap3(phip%u,phip%v,phip%t,klon,klat,klev)

    call setLateral_bc(klon,klat,klev,ksvar,phip,nlslan,RCAdomain)


    !     physical parameterization for the lagrangian model in phcall
    if(.not.nlphys)then
       if(mype==0)then
          print *,'not doing call to physics parametrization'
       endif
    endif

    if (nlphys) then
       zrdt = 1._realkind/zdt
       phip%ps = exp(phip%lnps)
       phip%ps = (phip%ps-phiz%ps)*zrdt
       phip%u = (phip%u-phiz%u)*zrdt
       phip%v = (phip%v-phiz%v)*zrdt
       phip%t = (phip%t-phiz%t)*zrdt
       phip%q = (phip%q-phiz%q)*zrdt
       phip%cw = (phip%cw-phiz%cw)*zrdt
       phip%svar = (phip%svar-phiz%svar)*zrdt
       svarsp = (svarsp-svarsz)*zrdt
       zdx_kf=hxv/(ra*rdlam)
       zdy_kf=hyu/(ra*rdth)

       zconacc = 1.0_realkind
       ztimeph = real(Nseconds(current_time-initialTime),realkind) !real(nstep)*Nseconds(ndtime) 

       eps = 1.e-6_realkind
       do jk=1,klev+1
          zahalf(jk) = RCAdomain%ahyb(jk)
          zbhalf(jk) = RCAdomain%bhyb(jk)
       enddo

       if( abs(RCAdomain%afull(1))<=eps .and. abs(RCAdomain%bfull(1))>eps ) then
          zahalf(1) = 0._realkind
          zbhalf(1) = 0.25_realkind*RCAdomain%bhyb(2)
       else if( abs(RCAdomain%afull(1))>eps .and. abs(RCAdomain%bfull(1))<=eps ) then
          zahalf(1) = 0.25_realkind*RCAdomain%ahyb(2)
          zbhalf(1) = 0._realkind
       endif

       !     vertical velocity computations needed for
       !     convection mass flux scheme

       call omcomp(phiz%u,phiz%v,phiz%ps,omfz,omhz,klon,klat,klev,RCAdomain%ahyb,RCAdomain%bhyb,RCAdomain%afull,&
            RCAdomain%bfull)
       call swap(omfz,klon,klat,klev)
       call swapklevp1(omhz,klon,klat,klev+1)

       !     total omega = omfz + omhz
       do jk = 1,klev
          do jy = 1,klat
             do jx = 1,klon
                omfz(jx,jy,jk)=omfz(jx,jy,jk)+omhz(jx,jy,jk)
                !     gjkf   get running mean omega into m/s for kfcumulus
                pkf=RCAdomain%afull(jk)+RCAdomain%bfull(jk)*zapspx(jx,jy)
                rokf=pkf/(287.04_realkind*(phiz%t(jx,jy,jk)*(1._realkind+0.61_realkind*phiz%q(jx,jy,jk))))
                omkf=omfz(jx,jy,jk)
                om_rmean(jx,jy,jk)=-(omkf/(9.81_realkind*rokf))
             enddo
          enddo
       enddo

       call phcall(nstep,&
            ztimeph,&
            klon,klat,klev,ksvar,&
            zdt,&                                            
            real(Nseconds(ndtime),realkind),  & 
            dtvdif,&
            dtphys,&
            zconacc,&
            timesu,&
            nlpost,&
            nldynvd,&
            nlinit,&
            phiz,&
            omfz,&
            surf_vars%tsea,&
            surf_vars%fri,&
            plong,&
            pclat,&
            pslat,&
            surf_vars%roc,&
            surf_vars%rou,&
            drolddt, dsolddt, accsunny,&
            accrunoff,&
            accrunoffopl,accrunofffor,accrunofflake,&
            q2d, &
            accprl, accprc,&
            slwr, sswr, tlwr, tswr,&
            slwrsea,sswrsea,&
            slwrice,sswrice,&
            accsnow,accsnoc,accsnol,&
            slwdn,sswdn,tswdn,&
            svarsz, svarsp, ksvars,&
            icethick, &
            guessid, &
            lcounttice,&
            phip,&
            totcov,cucov,&
            cov2d,cwpath,&
            dtdtph,dqdtph,dsdtph,&
            pblh,&
            msvars,svar_surf,&
            eco,meco,&
            lakes,&
            tsclim_years,&
            sacdg ,kacdg ,sacdg2 ,kacdg2,&
            dtdt_kf,dqdt_kf,dqcdt_kf,&
            oldcov,zvarcu,cwcu,div_kf,om_rmean,&
            kf_ind,kf_base,kf_top,&
            raincv_kf,snowcv_kf,umfb,shal_cgj,&
            nca,zdx_kf,zdy_kf)                       


       call swap2d(phip%ps,klon,klat)
       call swap5(phip%u,phip%v,phip%t,phip%q,phip%cw,klon,klat,klev)


       phip%ps   = phiz%ps+zdt*phip%ps
       phip%lnps = log(phip%ps)
       phip%u = phiz%u+zdt*phip%u
       phip%v = phiz%v+zdt*phip%v
       phip%t = phiz%t+zdt*phip%t
       phip%q = phiz%q+zdt*phip%q
       phip%cw = phiz%cw+zdt*phip%cw
       phip%svar = phiz%svar+zdt*phip%svar

       svarsp = svarsz+zdt*svarsp
       lakes%prog_lakes = lakes%prog_lakes+ zdt*lakes%tend_lakes
    endif



    !     implicit horisontal diffusion in diffh
    call diffh( klon   , klat   , klev,                          &
          phip%u,phip%v,phip%t,phip%q,phip%lnps,phip%cw,       &
          phip%svar,ksvar,zdt,atcref, npbpts,RCAdomain%ahyb,RCAdomain%bhyb  )                                                



    deallocate(apsx, alnpsx)
    deallocate(tx)
    deallocate(ux)
    deallocate(vx)
    deallocate(qx)
    deallocate(sx)
    return
  end subroutine sl2tim


  subroutine compfx(klon   , klat   , klev   , kcall  , &
        ppsm, plnpsm, pum, pvm, ptm, pqm, psm, pedotm,  &
        ppsz, plnpsz, puz, pvz, ptz, pqz, psz, pedotz,  &
        ppsx, plnpsx, pux, pvx, ptx, pqx, psx)
    !     
    !     compfx -  compute fields at timestep n+1/2
    !     purpose.
    !     --------
    !     
    !     this routine extrapolate values in time to timestep n+1/2
    !     in the semi-lagrangian two-time-level advection scheme
    !     
    !     interface.
    !     ----------
    !     
    !     compfx is called from the main time stepping routines
    !     
    !     input parameters:
    !     -----------------
    !     
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels
    !     kcall     switch between different formulas for extrapolation
    !     = 1: f(n+1/2) = f(n)
    !     = 2: f(n+1/2) = f(n,n-1)
    !     = 3: f(n+1/2) = f(n,n-1,n-2)
    !     
    !     ppsx      surface pressure
    !     plnpsx    ln( surface pressure )
    !     pux       velocity component in the x-direction
    !     pvx       velocity component in the y-direction
    !     ptx       temperature
    !     pqx       spesific humidity
    !     psx       passive scalar (liquid water)
    !     
    !     index x = x,m,z refer to timestep n-2,n-1 and n respectively
    !     
    !     output parameters:
    !     ------------------
    !     
    !     psx   
    !     plnpsx
    !     pux   
    !     pvx   
    !     ptx   
    !     pqx   
    !     psx   
    implicit none
    integer,intent(in):: klon,klat,klev,kcall !,khalo
    real(kind=realkind),intent(in)::ppsm(klon,klat),       plnpsm(klon,klat),  &
         pum(klon,klat,klev),     pvm(klon,klat,klev), &
         ptm(klon,klat,klev),     pqm(klon,klat,klev), &
         psm(klon,klat,klev),                          &
         pedotm(klon,klat,klev+1)
    real(kind=realkind),intent(in)::ppsz(klon,klat),      plnpsz(klon,klat),      &
         puz(klon,klat,klev),     pvz(klon,klat,klev), &
         ptz(klon,klat,klev),     pqz(klon,klat,klev), &
         psz(klon,klat,klev),                          &
         pedotz(klon,klat,klev+1)
    real(kind=realkind),intent(out)::ppsx(klon,klat),      plnpsx(klon,klat)

    real(kind=realkind),intent(out)::pux(klon,klat,klev), &
          pvx(klon,klat,klev), &
          ptx(klon,klat,klev), &
          pqx(klon,klat,klev), &
          psx(klon,klat,klev)

    character(len=80):: mes
    real(kind=realkind) zc1,zc2

    !     compute weights of different time levels and do extrapolation
    if (kcall==1) then
       !     copy from timelevel n
       pux = puz
       pvx = pvz
       ptx = ptz
       pqx = pqz
       psx = psz
       ppsx   = ppsz
       plnpsx = plnpsz

    elseif(kcall==2) then
       !     two-time-level extrapolation
       zc1 = 1.5_realkind
       zc2 = 0.5_realkind
       pux = zc1*puz - zc2*pum
       pvx = zc1*pvz - zc2*pvm
       ptx = zc1*ptz - zc2*ptm
       pqx = zc1*pqz - zc2*pqm
       psx = zc1*psz - zc2*psm
       ppsx = zc1*  ppsz - zc2*  ppsm
       plnpsx = zc1*plnpsz - zc2*plnpsm
    else
       write(mes,'(/,1x,''invalid call to compfx, kcall='',i9)')kcall
!       if(mype==0)print *,mes
       call stop_program(mes)
    endif

    return
  end subroutine compfx



  subroutine sldyn (klon   , klat   , klev    ,  &
        pfpara, &
        alnpsx,apspx,ux,vx,tx,qx,sx, &
        pdpsdt,pdudt,pdvdt,pdtdt,pdqdt,pdsdt,pdpdt,  &
        ltvir,lsl3d,fis,div_kf,ahyb,bhyb)


    !     sldyn - explicit dynamical tendency in semil-lagrangian scheme
    !     
    !     purpose.
    !     --------
    !     
    !     calculate the dynamic tendencies.
    !     input parameters:
    !     -----------------
    !     
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels
    !     pfpara    mean-value of the coriolis-parameter
    !     alnpsx    ln( surface pressure )
    !     apspx      surface pressure
    !     ux       velocity-component in x-direction
    !     vx       velocity-component in y-direction
    !     tx       temperature
    !     qx       specific humidity
    !     sx       passive scalar (liquid water)
    !     ltvir     true if virtual temperature corrections
    !     lsl3d     true if 3 dimensional semi-lagrangian advection
    !     
    !     output parameters:
    !     ------------------
    !     
    !     pdpsdt    tendency of surface pressure
    !     pdudt     tendency of velocity-component in x-direction
    !     pdvdt     tendency of velocity-component in y-direction
    !     pdtdt     tendency of temperature
    !     pdqdt     tendency of specific humidity
    !     pdsdt     tendency of passive scalar
    !     pdpdt     nonlinear terms in the continuity equation
    !     
    !     vertical advection of cloud water (s) removed (comments).
    !     horisontal advection of cloud water by upstream scheme.
    !     
    !     ---------------------------------------------------------------
    !     
    !     declaration of global parameters
    !     --------------------------------
    !     
    implicit none

    integer,intent(in)::klon, klat, klev

    real(kind=realkind),intent(in)::alnpsx(klon,klat)
    real(kind=realkind),intent(in)::apspx(klon+1,klat+1)
    real(kind=realkind),intent(in)::ux(klon,klat,klev)
    real(kind=realkind),intent(in)::vx(klon,klat,klev)
    real(kind=realkind),intent(in)::tx(klon,klat,klev)
    real(kind=realkind),intent(in)::qx(klon,klat,klev)
    real(kind=realkind),intent(in)::sx(klon,klat,klev)
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)

    real(kind=realkind),intent(out)::pdpsdt(klon,klat)
    real(kind=realkind),intent(out)::pdudt(klon,klat,klev)
    real(kind=realkind),intent(out)::pdvdt(klon,klat,klev)
    real(kind=realkind),intent(out)::pdtdt(klon,klat,klev)
    real(kind=realkind),intent(out)::pdqdt(klon,klat,klev)
    real(kind=realkind),intent(out)::pdsdt(klon,klat,klev)
    real(kind=realkind),intent(out)::pdpdt(klon,klat,klev)

    logical,intent(in):: ltvir, lsl3d
    real(kind=realkind),intent(in)::   pfpara, fis(klon,klat)
    real(kind=realkind),intent(inout):: div_kf(klon,klat,klev)

    !     declaration of local workspace
    real(kind=realkind)::dp_cgj
    integer:: jx,jy,jk, kp1,km1
    real(kind=realkind)::zrt0  , zffo4 , zalfa0, zbeta0,zrdlo2, zrdla2,                                        &
         zreps , zcpvd1, zdak  , zdbk,      &
         zcappa,  zrgasv, zcpv , zrdloh, zrdlah,   &
         zakm  , zbkm  , zlntwo, zrv , zrt  , zrdlor,           &
         zrdlar, zrgash

    real(kind=realkind)::  zrhxhy(klon,klat), zhxhy(klon,klat), zdpk(klon+1,klat+1),&
         zphi(klon,klat), zdsum(klon,klat), zpkm(klon+1,klat+1),  &
         zlogm(klon,klat), zlogp(klon,klat), zpkp(klon+1,klat+1),  &
         zrdpk(klon,klat), zpp(klon,klat), zdlnpk(klon,klat),      &
         zdivk(klon,klat), zomega(klon,klat), zalfa(klon,klat),      &
         zbeta(klon,klat), ztv(klon,klat), zedpde(klon,klat),      &
         zuu(klon,klat), zvv(klon,klat), zzk   (klon,klat),      &
         zek(klon,klat), zalfa1(klon,klat), zbeta1(klon,klat),      &
         zgradx(klon,klat), zgrady(klon,klat),                         &
         zw2(klon,klat), zlpspr(klon,klat)                          

    


    !     initialisation of physical constants
    zcappa = rair/cpair 
    zrgasv = 461.51_realkind
    zcpv   = 1869.46_realkind

    zreps  = zrgasv/rair
    zcpvd1 = zcpv/cpair-1.0_realkind
    zrt0   = rair*sit0
    zffo4  = pfpara*0.25_realkind

    do jy = 1,klat
       do jx = 1,klon
          zrhxhy(jx,jy) = ra*rhxu(jx,jy)*rhyv(jx,jy)
          zhxhy (jx,jy) = 1.0_realkind/zrhxhy(jx,jy)
       enddo
    enddo

    !     surface pressure tendency for diagnostic computation of sdot
    zuu = 0.0_realkind
    zvv = 0.0_realkind

    do  jk=1,klev
       zdak = ahyb(jk+1)-ahyb(jk)
       zdbk = bhyb(jk+1)-bhyb(jk)

       do jy = 1,klat+1
          do jx = 1,klon+1
             zdpk(jx,jy) = zdak + zdbk*apspx(jx,jy)
          enddo
       enddo

       do jy = 1,klat
          do jx = 1,klon
             zuu(jx,jy) = zuu(jx,jy)-ux(jx,jy,jk)*(zdpk(jx,jy)+zdpk(jx+1,jy))
             zvv(jx,jy) = zvv(jx,jy)-vx(jx,jy,jk)*(zdpk(jx,jy)+zdpk(jx,jy+1))
          enddo
       enddo
    enddo


    zrdloh = rdlam*0.5_realkind
    zrdlah = rdth*0.5_realkind

    pdpsdt(:,1) = 0.0_realkind
    pdpsdt(1,:) = 0.0_realkind

    !pdpsdt = D_(zuu)+D_(zvv)
    do jy = 2,klat
       do jx = 2,klon
          pdpsdt(jx,jy) = zrhxhy(jx,jy)*( &
               (zuu(jx,jy)*hyu(jx,jy)-zuu(jx-1,jy)*hyu(jx-1,jy))*zrdloh &
               +(zvv(jx,jy)*hxv(jx,jy)-zvv(jx  ,jy-1) *hxv(jx,jy-1))*zrdlah)
       enddo
    enddo

    !     store [ln(ps)]' in arraay zlpspr
    do jy = 1,klat
       do jx = 1,klon
          zlpspr(jx,jy) = alnpsx(jx,jy) + fis(jx,jy)/zrt0
       enddo
    enddo

    !     store tendency of lnps in array pdpdt
    do jk = 1,klev
       do jy = 2,klat
          do jx = 2,klon
             pdpdt(jx,jy,jk) = pdpsdt(jx,jy)/apspx(jx,jy)
          enddo
       enddo
    enddo

    !     initialization before vertical integration
    do jy = 1,klat+1
       do jx = 1,klon+1
          zpkm(jx,jy)   = apspx(jx,jy)
       enddo
    enddo

    do jy = 1,klat
       do jx = 1,klon
          zphi(jx,jy)   = -zrt0*alnpsx(jx,jy)      
          zlogm(jx,jy)  = alnpsx(jx,jy)
          zedpde(jx,jy) = 0.0_realkind
          zdsum (jx,jy) = 0.0_realkind
       enddo
    enddo

    !     vertical integration from bottom to top

    do 1000 jk=klev,1,-1
       kp1  = jk+1
       km1  = jk-1
       zakm = ahyb(jk)
       zbkm = bhyb(jk)
       zdbk = bhyb(jk+1)-bhyb(jk)
       zalfa0 = sigam1(jk)/rair
       zbeta0 = sigam2(jk)/rair - zalfa0

       !     pressure variabels
       do jy = 1,klat+1
          do jx = 1,klon+1
             zpkp(jx,jy) = zpkm(jx,jy)
             zpkm(jx,jy) = zakm + zbkm*apspx(jx,jy)
             zdpk(jx,jy) = zpkp(jx,jy) - zpkm(jx,jy)
          enddo
       enddo

       do jy = 1,klat
          do jx = 1,klon
             zrdpk(jx,jy) = 1.0_realkind/zdpk(jx,jy)
             zlogp(jx,jy) = zlogm(jx,jy) 
             pdudt(jx,jy,jk) = 0.0_realkind
             pdvdt(jx,jy,jk) = 0.0_realkind
             pdtdt(jx,jy,jk) = 0.0_realkind
             pdqdt(jx,jy,jk) = 0.0_realkind
             pdsdt(jx,jy,jk) = 0.0_realkind
          enddo
       enddo

       if (jk==1) then
          zlntwo = log(2.0_realkind)
          do jy = 1,klat
             do jx = 1,klon
                zlogm (jx,jy) = 0.0_realkind
                zdlnpk(jx,jy) = 0.0_realkind
                zpp   (jx,jy) = (zpkp(jx,jy)*zlogp(jx,jy) - zpkm(jx,jy)*zlogm(jx,jy) )*zrdpk(jx,jy)
                zalfa (jx,jy) =  zlntwo
                zbeta (jx,jy) = -zlntwo
             enddo
          enddo
          do jy = 1,klat
             do jx = 1,klon
                zalfa1(jx,jy) =  1.0_realkind
                zbeta1(jx,jy) = -1.0_realkind
             enddo
          enddo
       else
          do jy = 1,klat
             do jx = 1,klon
                zlogm(jx,jy) = log( zpkm(jx,jy) )
             enddo
          enddo
          do jy = 1,klat
             do jx = 1,klon
                zdlnpk(jx,jy) = zlogp(jx,jy) - zlogm(jx,jy)
                zpp   (jx,jy) = ( zpkp(jx,jy)*zlogp(jx,jy) - zpkm(jx,jy)*zlogm(jx,jy) )*zrdpk(jx,jy)
                zalfa (jx,jy) = 1.0_realkind - zpkm(jx,jy)*zrdpk(jx,jy) *zdlnpk(jx,jy)
                zbeta (jx,jy) = zdlnpk(jx,jy) - zalfa(jx,jy)
                zalfa1(jx,jy) = zalfa(jx,jy)
                zbeta1(jx,jy) = zbeta(jx,jy)
             enddo
          enddo
       endif

       !     virtual temperature tv and nonlinear part of geopotential phi-p
       if (ltvir) then
          !     take into account the mixing ratios of vapor, liquid and ice
          !     in the computation of the virtual temperature.
          !     at present ice mixing ratio is not known in the dynamics.
          !     tv = t*(1.+r(vapor)/eps)/(1.+r(vapor)+r(liquid)+r(ice))
          do jy = 1,klat
             do jx = 1,klon
                zrv = qx(jx,jy,jk)/(1.0_realkind-qx(jx,jy,jk))
                zrt = sx(jx,jy,jk)/(1.0_realkind-sx(jx,jy,jk))
                ztv(jx,jy) = tx(jx,jy,jk)*(1.0_realkind+zreps*zrv)/(1._realkind+zrv+zrt)
                zphi(jx,jy) = zphi(jx,jy)+zalfa(jx,jy)*rair*ztv(jx,jy)-zalfa0*rair*tx(jx,jy,jk)
             enddo
          enddo
       else
          do jy = 1,klat
             do jx = 1,klon
                ztv(jx,jy) = tx(jx,jy,jk)
                zphi(jx,jy)=zphi(jx,jy)+zalfa(jx,jy)*rair*ztv(jx,jy) - zalfa0  *rair*tx(jx,jy,jk)
             enddo
          enddo
       endif
       !     auxiliary velocities
       do jy = 1,klat
          do jx = 1,klon
             zuu(jx,jy)= 0.5_realkind * (zdpk(jx+1,jy) + zdpk(jx,jy)) *  ux(jx,jy,jk)* hyu(jx,jy)
             zvv(jx,jy) = 0.5_realkind * (zdpk(jx,jy+1) + zdpk(jx,jy))*  vx(jx,jy,jk)* hxv(jx,jy)
          enddo
       enddo


       !     divergence 
       do jy = 2,klat
          do jx = 2,klon
             zdivk(jx,jy) = zrhxhy(jx,jy) * ( ( zuu(jx,jy) - zuu(jx-1,jy  ) ) * rdlam &
                  + ( zvv(jx,jy) - zvv(jx  ,jy-1) ) * rdth )
             dp_cgj = 0.25_realkind*(zdpk(jx+1,jy)+zdpk(jx,jy)+ zdpk(jx,jy+1)+zdpk(jx,jy))
             div_kf(jx,jy,jk) = zdivk(jx,jy)/dp_cgj
             if(div_kf(jx,jy,jk)<1.e-12_realkind .and. div_kf(jx,jy,jk)>-1.e-12_realkind) &
                div_kf(jx,jy,jk)=0._realkind
          enddo
       enddo

       !     vertical advection terms
       !     inflow from below
       if(jk/=klev .and. .not.lsl3d)then
          do jy = 2,klat-1
             do jx = 2,klon-1
                pdtdt(jx,jy,jk) = pdtdt(jx,jy,jk) - 0.5_realkind * zrdpk(jx,jy)* zedpde(jx,jy) &
                     * ( tx(jx,jy,kp1)-tx(jx,jy,jk) )
                pdqdt(jx,jy,jk) = pdqdt(jx,jy,jk) - 0.5_realkind * zrdpk(jx,jy)* zedpde(jx,jy) &
                     *( qx(jx,jy,kp1)-qx(jx,jy,jk) )

!                pdsdt(jx,jy,jk) = pdsdt(jx,jy,jk)-0.5_realkind * zrdpk(jx,jy)* zedpde(jx,jy) &
!                     *( sx(jx,jy,kp1)-sx(jx,jy,jk) )

             enddo
          enddo

          do jy = 2,klat-1
             do jx = 2,klon-1
                pdudt(jx,jy,jk) = pdudt(jx,jy,jk) - 0.5_realkind /( zdpk(jx,jy) +   zdpk(jx+1,jy   ) ) * &
                     (zedpde(jx,jy) + zedpde(jx+1,jy   ) ) * ( ux(jx,jy,kp1) - ux(jx,jy,jk) )
             enddo
          enddo

          do jy = 2,klat-1
             do jx = 2,klon-1
                pdvdt(jx,jy,jk) = pdvdt(jx,jy,jk) -  0.5_realkind /( zdpk(jx,jy) +  zdpk(jx,jy+1) ) * &
                     (zedpde(jx,jy)+ zedpde(jx,jy+1) ) * ( vx(jx,jy,kp1)-vx(jx,jy,jk) )
             enddo
          enddo
       endif

       !     update vertical velocity
       if(jk/=1)then
          do jy = 2,klat
             do jx = 2,klon
                zedpde(jx,jy) = zedpde(jx,jy) +  zdbk*pdpsdt(jx,jy) + zdivk(jx,jy)
             enddo
          enddo

          !     inflow from above
          if(.not.lsl3d)then
             do jy = 2,klat-1
                do jx = 2,klon-1
                   pdtdt(jx,jy,jk)=pdtdt(jx,jy,jk) - 0.5_realkind * zrdpk(jx,jy)*zedpde(jx,jy)* ( tx(jx,jy,jk) -tx(jx,jy,km1) )
                   pdqdt(jx,jy,jk)=pdqdt(jx,jy,jk) - 0.5_realkind * zrdpk(jx,jy)*zedpde(jx,jy)* ( qx(jx,jy,jk) -qx(jx,jy,km1) )


!                   pdsdt(jx,jy,jk)=pdsdt(jx,jy,jk)-0.5_realkind * zrdpk(jx,jy)*zedpde(jx,jy) &
!                        *( sx(jx,jy,jk) -sx(jx,jy,km1) )

                enddo
             enddo

             do jy = 2,klat-1
                do jx = 2,klon-1
                   pdudt(jx,jy,jk) = pdudt(jx,jy,jk)-0.5_realkind/(zdpk(jx,jy)+zdpk(jx+1,jy))*(zedpde(jx,jy)+zedpde(jx+1,jy))* &
                        ( ux(jx,jy,jk)-ux(jx,jy,km1) )
                enddo
             enddo

             do jy = 2,klat-1
                do jx = 2,klon-1
                   pdvdt(jx,jy,jk) = pdvdt(jx,jy,jk) - 0.5_realkind /( zdpk(jx,jy)+  zdpk(jx,jy+1) ) * &
                        (zedpde(jx,jy)+zedpde(jx,jy+1) ) * ( vx(jx,jy,jk) -vx(jx,jy,km1) )
                enddo
             enddo
          endif
       endif

       !     nonlinear part of energy conversion term
       if (ltvir) then
          do jy = 2,klat-1
             do jx = 2,klon-1
                zomega(jx,jy)= zcappa/( 1.0_realkind + zcpvd1*qx(jx,jy,jk) )
             enddo
          enddo
       else
          do jy = 2,klat-1
             do jx = 2,klon-1
                zomega(jx,jy)= zcappa
             enddo
          enddo
       endif

       zrdloh = rdlam*0.5_realkind
       zrdlah = rdth*0.5_realkind

       do jy = 2,klat-1
          do jx = 2,klon-1
             zomega(jx,jy) =  zomega(jx,jy) * zrdpk(jx,jy) * ( ztv(jx,jy) *  &
                  ( zdlnpk(jx,jy) * ( zdsum(jx,jy) + pdpsdt(jx,jy) )          &
                  + zbeta1(jx,jy) *   zdivk(jx,jy)               )            &
                  + 0.25_realkind * zrhxhy(jx,jy) *                                    &
                  ( (   zuu(jx,jy)*                                       &
                  ( ztv(jx+1,jy ) + ztv(jx,jy   ) ) *                         &
                  ( zpp(jx+1,jy ) - zpp(jx,jy   ) )                           &
                  + zuu(jx-1,jy ) *                                           &
                  ( ztv(jx,jy   ) + ztv(jx-1,jy ) ) *                         &
                  ( zpp(jx,jy   ) - zpp(jx-1,jy ) )   ) * rdlam               &
                  + (   zvv(jx,jy   ) *                                       &
                  ( ztv(jx,jy+1 ) + ztv(jx,jy   ) ) *                         &
                  ( zpp(jx,jy+1 ) - zpp(jx,jy   ) )                           &
                  + zvv(jx,jy-1 ) *                                           &
                  ( ztv(jx,jy   ) + ztv(jx,jy-1) ) *                          &
                  ( zpp(jx,jy   ) - zpp(jx,jy-1) )   ) * rdth ) )           

             pdtdt(jx,jy,jk) = pdtdt(jx,jy,jk) + zomega(jx,jy) 

             !     add the advection of [ln(ps)]' at level k

             pdpdt(jx,jy,jk) = pdpdt(jx,jy,jk) +                       &
                  zrdpk(jx,jy) * zrhxhy(jx,jy) *                        &
                  ( zrdloh *                                            &
                  ( zuu(jx  ,jy) * (zlpspr(jx+1,jy) - zlpspr(jx,jy) )   &
                  - zuu(jx-1,jy) * (zlpspr(jx-1,jy) - zlpspr(jx,jy)))   &
                  + zrdlah *                                            &
                  ( zvv(jx,jy  ) * (zlpspr(jx,jy+1) - zlpspr(jx,jy))    &
                  - zvv(jx,jy-1) * (zlpspr(jx,jy-1) - zlpspr(jx,jy))))
          enddo
       enddo
       !     absolute vorticity and energy
       zrdlo2 = 2.0_realkind*rdlam
       zrdla2 = 2.0_realkind*rdth

       do jy = 1,klat-1
          do jx = 1,klon-1
             zzk(jx,jy) = ( fpar(jx,jy) * ( zhxhy(jx  ,jy  )       &
                  + zhxhy(jx+1,jy  )                                &
                  + zhxhy(jx  ,jy+1)                                &
                  + zhxhy(jx+1,jy+1) )                              &
                  + zrdlo2*( vx(jx,jy,jk)     + vx(jx+1,jy,jk) )* &
                  ( 1._realkind/rhyv(jx+1,jy ) - 1._realkind/rhyv(jx,jy)  )           &
                  - zrdla2*( ux(jx,jy,jk)     + ux(jx,jy+1,jk) )* &
                  ( 1._realkind/rhxu(jx,jy+1)  - 1._realkind/rhxu(jx,jy)  )           &
                  ) / ( zhxhy(jx  ,jy  )*zdpk(jx  ,jy  )            &
                  + zhxhy(jx+1,jy  )*zdpk(jx+1,jy  )                &
                  + zhxhy(jx  ,jy+1)*zdpk(jx  ,jy+1)                &
                  + zhxhy(jx+1,jy+1)*zdpk(jx+1,jy+1) )
          enddo
       enddo

       zrdlor = rdlam*ra !/zrear
       zrdlar = rdth*ra !/zrear
       zrgash = rair*0.5_realkind

       !     pressure gradient, vorticity and energy minus linearized coriolis
       do jy = 2,klat-1
          do jx = 2,klon-1
             zgradx(jx,jy) = zrdlor * rhxu(jx,jy) *        &
                  ( zphi(jx+1,jy   )-zphi(jx,jy)            &
                  + zrgash*( ztv(jx,jy)+ztv(jx+1,jy   ) ) * &
                  ( zpp(jx+1,jy   )-zpp(jx,jy) ) )

             pdudt(jx,jy,jk) = pdudt(jx,jy,jk) -                 &
                  zgradx(jx,jy) + rhxu (jx,jy) *                  &
                  0.125_realkind * ( zzk(jx,jy-1)+zzk(jx  ,jy  ) )         &
                  * ( zvv(jx,jy  )+zvv(jx+1,jy  )                 &
                  +zvv(jx,jy-1)+zvv(jx+1,jy-1) )                  &
                  - zffo4 * ( vx(jx,jy  ,jk) + vx(jx+1,jy  ,jk) &
                  + vx(jx,jy-1,jk) + vx(jx+1,jy-1,jk) )
          enddo
       enddo

       do jy = 2,klat-1
          do jx = 2,klon-1
             zgrady(jx,jy) = zrdlar * rhyv (jx,jy) *       &
                  ( zphi(jx,jy+1)-zphi(jx,jy)               &
                  + zrgash*( ztv(jx,jy)+ztv(jx,jy+1) ) *    &
                  ( zpp(jx,jy+1)-zpp(jx,jy) ) )

             pdvdt(jx,jy,jk) = pdvdt(jx,jy,jk) -                 &
                  zgrady(jx,jy) - rhyv (jx,jy) *                  &
                  0.125_realkind * ( zzk(jx-1,jy  )+zzk(jx,jy  ) )         &
                  * ( zuu(jx-1,jy+1)+zuu(jx,jy+1)                 &
                  +zuu(jx-1,jy  )+zuu(jx,jy  ) )                  &
                  + zffo4 * ( ux(jx,jy  ,jk) + ux(jx-1,jy  ,jk) &
                  + ux(jx,jy+1,jk) + ux(jx-1,jy+1,jk) )
          enddo
       enddo

       !     update vertical integrated geopotential (phi-p) and divergence
       if (jk>=2 ) then
          do jy = 1,klat
             do jx = 1,klon
                zphi(jx,jy)=zphi(jx,jy)+zbeta(jx,jy)*rair*ztv(jx,jy) - zbeta0*rair*tx(jx,jy,jk)
             enddo
          enddo
          do jy = 2,klat-1
             do jx = 2,klon-1
                zdsum(jx,jy) = zdsum(jx,jy) + zdivk(jx,jy)
             enddo
          enddo
       endif
1000 enddo

    zw2 = 0.0_realkind

    do jk=1,klev
       do jy = 2,klat-1
          do jx = 2,klon-1
             !     compute the divergence
             zek(jx,jy) = zrhxhy(jx,jy) *                      &
                  ( ( hyu(jx  ,jy  )*ux(jx  ,jy  ,jk) -        &
                  hyu(jx-1,jy  )*ux(jx-1,jy  ,jk) ) * rdlam    &
                  + ( hxv(jx  ,jy  )*vx(jx  ,jy  ,jk) -        &
                  hxv(jx  ,jy-1)*vx(jx  ,jy-1,jk) ) * rdth )

             !     subtract the linear term in the energy conversion term from the
             !     eulerian term to complete the computation of nt
             pdtdt(jx,jy,jk) = pdtdt(jx,jy,jk)+ sitau1(jk)*zw2(jx,jy) + sitau2(jk)*zek(jx,jy)

             !     accumulate the sum of divergence
             zw2(jx,jy) = zw2(jx,jy) + sidpk0(jk)*zek(jx,jy)
          enddo
       enddo
    enddo

    !     add the linear term to pdpdt to complete the computation of np

    do jk=1,klev
       do jy = 2,klat-1
          do jx = 2,klon-1
             pdpdt(jx,jy,jk) = pdpdt(jx,jy,jk) + zw2(jx,jy)/sip0
          enddo
       enddo
    enddo

    return
  end subroutine sldyn



  subroutine slexpa(klon   , klat   , klev   , kpbpts,  &
       ptwodt , pfpara,                                   &
       fis,                                               &
       plnpsm , pum    , pvm    , ptm,                     &
       plnpsz , puz    , pvz    , ptz,                     &
       plnpsp , pup    , pvp    , ptp,                     &
       pdiv )
    ! slexpa - semi-lagrangian explicit adjustment
    !
    !     purpose.
    !     --------
    !
    !     this routine prepares the right hand side of the helmholtz
    !     equation for the divergence in the semi-implicit scheme in
    !     case of semi-lagrangian advection.
    !
    !   interface.
    !     ----------
    !
    !     slexpa is called from the main time stepping routines
    !
    !     input parameters:
    !     -----------------
    !
    !     klon      number of gridpoints in the x-direction
    !     klat      number of gridpoints in the y-direction
    !     klat      number of vertical levels
    !     kpbpts    number of passive boundary lines
    !     ptwodt    current double timestep
    !     pfpara    mean-value of the coriolis-parameter
    !
    !     plnpsx    ln( surface pressure )
    !     pux       velocity component in the x-direction
    !     pvx       velocity component in the y-direction
    !     ptx       temperature
    !
    !     index x = m,z,p refer to previous, current and following timestep
    !
    !     output parameters:
    !     ------------------
    !
    !     pup       adjusted velocity component in the x-direction
    !     pvp       adjusted velocity component in the y-direction
    !     pdiv      right hand side of the helmholtz-equation
    implicit none

    integer,intent(in):: klon,klat,klev,kpbpts

    real(kind=realkind),intent(in):: ptwodt,pfpara

    real(kind=realkind),intent(in)::  plnpsm(klon,klat)     , pum(klon,klat,klev), &
         pvm(klon,klat,klev), ptm(klon,klat,klev),  &
         plnpsz(klon,klat)     , puz(klon,klat,klev),  &
         pvz(klon,klat,klev), ptz(klon,klat,klev)
    real(kind=realkind),intent(inout)::plnpsp(klon,klat)     , pup(klon,klat,klev),  &
         pvp(klon,klat,klev), ptp(klon,klat,klev)
    real(kind=realkind),intent(inout)::pdiv(klon,klat,klev)
    real(kind=realkind),intent(in)::fis(klon,klat)   
    integer:: i,j,k,kpbpts_base,kpbpts_top,kpbpts_left,kpbpts_right

    real(kind=realkind):: zdt,zrt0,zdtrdx,zdtrdy,zrdlam,zrdth,zffcon,zaa,zgam1,zgam2

    real(kind=realkind) zw1(klon,klat), zw2(klon,klat), zw3(klon,klat)


    !    fix various local constants

    zdt    = (1.0_realkind+epsg)*0.5_realkind*ptwodt


    zrt0   = rair*sit0
    zdtrdx = zdt*ra*rdlam
    zdtrdy = zdt*ra*rdth
    zrdlam = rdlam*ra
    zrdth  = rdth *ra
    zffcon = zdt*pfpara*0.25_realkind
    zaa    = zdt*pfpara
    zaa    = 1.0_realkind/(1.0_realkind+zaa*zaa)

    !    fix local passive boundary point counters
    if( attop ) then
       kpbpts_top = kpbpts
    else
       kpbpts_top = 0
    endif

    if( atbase ) then
       kpbpts_base = kpbpts
    else
       kpbpts_base = 0
    endif

    if( atright ) then
       kpbpts_right = kpbpts
    else
       kpbpts_right = 0
    endif

    if( atleft ) then
       kpbpts_left = kpbpts
    else
       kpbpts_left = 0
    endif

    !    give zero values to divergence field and work arrays
    do k=1,klev
       do j=1,klat
          do i=1,klon
             pdiv(i,j,k) = 0.0_realkind
          enddo
       enddo
    enddo

    do j=1,klat
       do i=1,klon
          zw1(i,j) = 0.0_realkind
          zw2(i,j) = 0.0_realkind
          zw3(i,j) = 0.0_realkind
       enddo
    enddo

    !    give boundary values to ln(ps) for pressure gradient term
    if( atbase ) then
       do i=1+kpbpts_left,klon-kpbpts_right
          plnpsp(i,1+kpbpts) = 2.0_realkind*plnpsz(i,1+kpbpts)-plnpsm(i,1+kpbpts)
       enddo
    endif

    if( attop ) then
       do i=1+kpbpts_left,klon-kpbpts_right
          plnpsp(i,klat-1-kpbpts)=2.0_realkind*plnpsz(i,klat-1-kpbpts)-&
               plnpsm(i,klat-1-kpbpts)
       enddo
    endif

    if( atleft ) then
       do j=1+kpbpts_base,klat-kpbpts_top
          plnpsp(1+kpbpts,j)=2.0_realkind*plnpsz(1+kpbpts,j)-plnpsm(1+kpbpts,j)
       enddo
    endif

    if( atright ) then
       do j=1+kpbpts_base,klat-kpbpts_top
          plnpsp(klon-1-kpbpts,j)=2.0_realkind*plnpsz(klon-1-kpbpts,j)- & 
               plnpsm(klon-1-kpbpts,j)
       enddo
    endif

    !    initiate pressure gradient term
    do j=1,klat
       do i=1,klon
          zw1(i,j) = fis(i,j) + zrt0 * plnpsp(i,j)
       enddo
    enddo

    !    verical integration from bottom to top
    do k=klev,1,-1
       !    give boundary values to temperature for pressure gradient term
       if( atbase ) then
          do i=1+kpbpts_left,klon-kpbpts_right
             ptp(i,1+kpbpts,k)=2.0_realkind*ptz(i,1+kpbpts,k)-ptm(i,1+kpbpts,k)
          enddo
       endif

       if( attop ) then
          do i=1+kpbpts_left,klon-kpbpts_right
             ptp(i,klat-1-kpbpts,k) = 2.0_realkind*ptz(i,klat-1-kpbpts,k) -  &
                  ptm(i,klat-1-kpbpts,k)
          enddo
       endif

       if( atleft ) then
          do j=1+kpbpts_base,klat-kpbpts_top
             ptp(1+kpbpts,j,k)=2.0_realkind*ptz(1+kpbpts,j,k)-ptm(1+kpbpts,j,k)
          enddo
       endif

       if( atright ) then
          do j=1+kpbpts_base,klat-kpbpts_top
             ptp(klon-1-kpbpts,j,k)=2.0_realkind*ptz(klon-1-kpbpts,j,k) - &
                  ptm(klon-1-kpbpts,j,k)
          enddo
       endif


       !    pressure gradient term
       zgam1 = sigam1(k)
       zgam2 = sigam2(k)

       do j=1,klat
          do i=1,klon
             zw2(i,j) = zgam1 * ptp(i,j,k) + zw1(i,j)
          enddo
       enddo

       !    update vertical integrated sum
       if (k>1) then
          do j=1,klat
             do i=1,klon
                zw1(i,j) = zw1(i,j) + zgam2 * ptp(i,j,k)
             enddo
          enddo
       endif

       !    give boundary values to velocities for divergence
       if( atbase ) then
          do i=1+kpbpts_left,klon-kpbpts_right
             pup(i,1+kpbpts,k)=2.0_realkind*puz(i,1+kpbpts,k) - &
                  pum(i,1+kpbpts,k)
             pvp(i,1+kpbpts,k) = 2.0_realkind*pvz(i,1+kpbpts,k) -  &
                  pvm(i,1+kpbpts,k)
          enddo
       endif

       if( attop ) then
          do i=1+kpbpts_left,klon-kpbpts_right
             pup(i,klat-1-kpbpts,k ) = 2.0_realkind*puz(i,klat-1-kpbpts,k) - &
                  pum(i,klat-1-kpbpts,k)
             pvp(i,klat-1-kpbpts,k ) = 2.0_realkind*pvz(i,klat-1-kpbpts,k) -  &
                  pvm(i,klat-1-kpbpts,k)
          enddo
       endif

       if( atleft ) then
          do j=1+kpbpts_base,klat-kpbpts_top
             pup(1+kpbpts,j,k) = 2.0_realkind*puz(1+kpbpts,j,k) -  &
                  pum(1+kpbpts,j,k)
             pvp(1+kpbpts,j,k) = 2.0_realkind*pvz(1+kpbpts,j,k) -  &
                  pvm(1+kpbpts,j,k)
          enddo
       endif

       if( atright ) then
          do j=1+kpbpts_base,klat-kpbpts_top
             pup(klon-1-kpbpts,j,k) = 2.0_realkind*puz(klon-1-kpbpts,j,k) - &
                  pum(klon-1-kpbpts,j,k)
             pvp(klon-1-kpbpts,j,k) = 2.0_realkind*pvz(klon-1-kpbpts,j,k) - &
                  pvm(klon-1-kpbpts,j,k)
          enddo
       endif

       !    modify velocities
       do j=1,klat
          do i=1,klon-1
             pup(i,j,k) = pup(i,j,k)-zdtrdx*rhxu(i,j)*(zw2(i+1,j)-zw2(i,j))
          enddo
       enddo

       do j=1,klat-1
          do i=1,klon
             pvp(i,j,k)=pvp(i,j,k)-zdtrdy*rhyv(i,j)*(zw2(i,j+1)-zw2(i,j))
          enddo
       enddo

       !    decopling of velocity-components
       do j=2,klat-1
          do i=1,klon-1
             zw2(i,j) = pup(i,j,k)                     &
                  + zffcon*( pvp(i,j,k) + pvp(i+1,j,k) &
                  + pvp(i,j-1,k) + pvp(i+1,j-1,k) )
          enddo
       enddo

       do j=1,klat-1
          do i=2,klon-1
             zw3(i,j) = pvp(i,j,k) &
                  - zffcon*( pup(i,j,k) + pup(i-1,j,k) &
                  + pup(i,j+1,k) + pup(i-1,j+1,k) )
          enddo
       enddo

       !    divergence
       do j=2,klat-1
          do i=2,klon-1
             pdiv(i,j,k) = zaa * rhxu(i,j)*rhyv(i,j) *        &
                  ( zrdlam*                                      &
                  (zw2(i,j)*hyu(i,j) - zw2(i-1,j)*hyu(i-1,j))  &
                  + zrdth *                                      &
                  (zw3(i,j)*hxv(i,j) - zw3(i,j-1)*hxv(i,j-1)) )
          enddo
       enddo
    enddo
    return
  end subroutine slexpa

  subroutine impadj(klon   , klat   , klev   , kpbpts, &
        ptwodt , pfpara,                               &
        plnpsp , pup    , pvp    , ptp,                &
        pdiv)
    !     
    !     impadj - implicit semi-implicit adjustment
    !     
    !     purpose.
    !     final part of the semi-implicit adjustment
    !     
    !     input parameters:
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels
    !     kpbpts    number of passiv boundarypoints
    !     ptwodt    current double timestep
    !     pfpara    mean-value of the coriolis-parameter
    !     pdiv      the solution of the helmholtz-equation
    !     
    !     output parameters:
    !     plnpsp    ln( surface pressure )
    !     pup       velocity component in the x-direction
    !     pvp       velocity component in the y-direction
    !     ptp       temperature
    implicit none

    integer,intent(in)::klon   , klat   , klev   , kpbpts
    real(kind=realkind),intent(in)::ptwodt,pfpara
    real(kind=realkind),intent(inout):: plnpsp(klon,klat),                          &
         pup(klon,klat,klev),pvp(klon,klat,klev), &
         ptp(klon,klat,klev)
    real(kind=realkind),intent(inout)::pdiv(klon,klat,klev)


    !     declaration of local workspace
    real(kind=realkind)::zw1(klon,klat), zw2(klon,klat), zw3(klon,klat)
    integer::jx,jy,jk,istart,istop,istartp1,istopm1,jstart,jstartp1 
    integer::jstop,jstopm1 !, jkp1
    integer::kpbpts_top, kpbpts_base , kpbpts_right, kpbpts_left   
    real(kind=realkind)::zdt,zdtrdx,zdtrdy         
    real(kind=realkind)::zgam1,zgam2,zrdt,zdtrp0, zdpk0         
    real(kind=realkind)::zffcon,zaa,ztau2,ztau1,zconst !, zbk

    kpbpts_top = 0
    kpbpts_base = 0
    kpbpts_right = 0
    kpbpts_left = 0

    if( attop ) then
       kpbpts_top = kpbpts
    endif
    if( atbase ) then
       kpbpts_base = kpbpts
    endif
    if( atright ) then
       kpbpts_right = kpbpts
    endif
    if( atleft ) then
       kpbpts_left = kpbpts
    endif

    istart   =  kpbpts_left + 1
    istartp1 =  istart + 1
    istop    =  klon -  kpbpts_right
    istopm1  =  istop  - 1
    jstart   =  kpbpts_base + 1
    jstartp1 =  jstart + 1
    jstop    =  klat -  kpbpts_top
    jstopm1  =  jstop  - 1

    !     give zero value to work arrays
    do jy = 1,klat
       do jx = 1,klon
          zw1(jx,jy) = 0.0_realkind
          zw2(jx,jy) = 0.0_realkind
          zw3(jx,jy) = 0.0_realkind
       enddo
    enddo



    zdt    = (1.0_realkind+epsg)*0.5_realkind*ptwodt

    zrdt   = 1.0_realkind/zdt

    zdtrp0 = zdt/sip0
    zdtrdx = zdt*ra*rdlam
    zdtrdy = zdt*ra*rdth
    zffcon = zdt*pfpara*0.25_realkind
    zaa    = zdt*pfpara
    zaa    = 1.0_realkind/(1.0_realkind+zaa*zaa)
    !     integration of continuity equation and adjustment of temperature
    zdpk0 =  sidpk0(1)
    ztau2 = -sitau2(1) * zdt

    do jy = jstart,jstop
       do jx = istart,istop
          zw1(jx,jy)   = zdpk0 * pdiv(jx,jy,1)
          pdiv(jx,jy,1) = ztau2 * pdiv(jx,jy,1)
          ptp(jx,jy,1) = ptp(jx,jy,1) + pdiv(jx,jy,1)
       enddo
    enddo

    do jk=2,klev
       ztau1 = -sitau1(jk) * zdt
       ztau2 = -sitau2(jk) * zdt
       zdpk0 =  sidpk0(jk)
       do jy = jstart,jstop
          do jx = istart,istop
             zw2(jx,jy)     = ztau1 * zw1(jx,jy)
             zw1(jx,jy)     = zw1(jx,jy) + zdpk0 * pdiv(jx,jy,jk)
             pdiv(jx,jy,jk) = zw2(jx,jy) + ztau2 * pdiv(jx,jy,jk)
             ptp(jx,jy,jk)  = ptp(jx,jy,jk) + pdiv(jx,jy,jk)
          enddo
       enddo
    enddo
    !     adjustment of ln(ps). initiate pressure gradient term
    zconst = -zdtrp0 * rair * sit0

    do jy = jstart,jstop
       do jx = istart,istop
          plnpsp(jx,jy) = plnpsp(jx,jy) - zdtrp0 * zw1(jx,jy)
          zw1(jx,jy)    = zconst * zw1(jx,jy)
       enddo
    enddo

    do jk=klev,1,-1
       !     pressure gradient term
       zgam1 = sigam1(jk)
       do jy = jstart,jstop
          do jx = istart,istop
             zw2(jx,jy) = zgam1 * pdiv(jx,jy,jk) + zw1(jx,jy)
          enddo
       enddo

       !     boundary values of pressure gradient term
       if( attop ) then
          do jx = istart,istop
             zw2(jx,klat-1-kpbpts) = 0.0_realkind
          enddo
       endif
       if( atbase ) then
          do jx = istart,istop
             zw2(jx,kpbpts+1)      = 0.0_realkind
          enddo
       endif
       if( atright ) then
          do jy=jstart,jstop
             zw2(klon-1-kpbpts,jy) = 0.0_realkind
          enddo
       endif
       if( atleft ) then
          do jy=jstart,jstop
             zw2(kpbpts+1,jy)      = 0.0_realkind
          enddo
       endif

       !     update vertical sum
       if (jk>1) then
          zgam2 = sigam2(jk)
          do jy = jstart,jstop
             do jx = istart,istop
                zw1(jx,jy) = zw1(jx,jy) + zgam2 * pdiv(jx,jy,jk)
             enddo
          enddo
       endif
       !     adjustment of velocities

       do jy = jstart,jstop
          do jx = istart,istopm1
             pup(jx,jy,jk) = pup(jx,jy,jk) &
                  - zdtrdx * rhxu(jx,jy) * (zw2(jx+1,jy) - zw2(jx,jy))
          enddo
       enddo

       do jy = jstart,jstopm1
          do jx = istart,istop
             pvp(jx,jy,jk) = pvp(jx,jy,jk) &
                  - zdtrdy * rhyv(jx,jy) * (zw2(jx,jy+1) - zw2(jx,jy))
          enddo
       enddo
       !     decopling of velocity-components


       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1
             zw2(jx,jy) = pup(jx,jy,jk) + zffcon *         &
                  ( pvp(jx,jy  ,jk) + pvp(jx+1,jy  ,jk) +   &
                  pvp(jx,jy-1,jk) + pvp(jx+1,jy-1,jk) )      
             zw3(jx,jy) = pvp(jx,jy,jk)  - zffcon *         &
                  ( pup(jx,jy  ,jk) + pup(jx-1,jy  ,jk) +   &
                  pup(jx,jy+1,jk) + pup(jx-1,jy+1,jk) )
          enddo
       enddo

       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1
             pup(jx,jy,jk) = zaa*zw2(jx,jy)
             pvp(jx,jy,jk) = zaa*zw3(jx,jy)
          enddo
       enddo

    enddo

    return
  end subroutine impadj



 


  subroutine omcomp(u,v,ps,om1,om2,klon,klat,klev,ahalf,bhalf,afull,bfull)

    implicit none

    !     omcomp - compute half and full level part of omega

    integer,intent(in):: klon,klat,klev
    real(kind=realkind),intent(in),dimension(klon,klat,klev):: u,v
    real(kind=realkind),intent(in),dimension(klon,klat):: ps
    real(kind=realkind),intent(out):: om1(klon,klat,klev), om2(klon,klat,klev+1)
    real(kind=realkind),intent(in):: ahalf(klev+1),bhalf(klev+1),afull(klev),  bfull(klev)

    integer j,i,k
    real(kind=realkind):: zda,zdb,zpm1,zpm2,zpp2,zpm3,zpp3,zpp1

    real(kind=realkind):: utemp(klon),vtemp(klon,2)
    real(kind=realkind):: dptemp(klon,3),ptemp(klon,3)

    !     dptemp = dp=p(k+.5)-p(k-.5)  at t-latitudes
    !     ptemp  = ln(p)               at t-latitudes
    !     utemp  = u * dp * hyu        at u-latitudes
    !     vtemp  = v * dp * hxv        at v-latitudes

    om2 = 0.0_realkind
    om1 = 0.0_realkind

    do  k=1,klev
       zda = ahalf(k+1) - ahalf(k)
       zdb = bhalf(k+1) - bhalf(k)
       !     initiate dptemp and ptemp at the two southern t-latitudes
       !     and vtemp at southern v-latitude
       do i=1,klon
          dptemp(i,1) = zda + zdb*ps(i,1)
          dptemp(i,2) = zda + zdb*ps(i,2)
          vtemp(i,1) = v(i,1,k)*hxv(i,1)*0.5_realkind*( dptemp(i,1)+dptemp(i,2) )
       enddo

       if ((abs(afull(k))>1.e-14_realkind).and.(abs(bfull(k))>1.e-14_realkind))then  !hybrid surface
          do i=1,klon
             zpm1 = ahalf(k  ) + bhalf(k  )*ps(i,1)
             zpp1 = ahalf(k+1) + bhalf(k+1)*ps(i,1)
             ptemp(i,1) = ( zpp1* log(zpp1) - zpm1* log(zpm1) )/dptemp(i,1)
             zpm2 = ahalf(k  ) + bhalf(k  )*ps(i,2)
             zpp2 = ahalf(k+1) + bhalf(k+1)*ps(i,2)
             ptemp(i,2) = ( zpp2* log(zpp2) - zpm2* log(zpm2) )/dptemp(i,2)
          enddo
       elseif((abs(afull(k))< 1.e-14_realkind).and.(abs(bfull(k))>1.e-14_realkind))then  !     sigma surface
          do  i=1,klon
             ptemp(i,1) =  log(ps(i,1))
             ptemp(i,2) =  log(ps(i,2))
          enddo
       endif

       do  j=2,klat-1
          !     evaluate dptemp and ptemp at next t-latitude,
          !     vtemp at this v-latitude and
          !     utemp at this u-latitude
          do i=1,klon
             dptemp(i,3) = zda + zdb*ps(i,j+1)
             vtemp(i,2) = v(i,j,k)*hxv(i,j)*0.5_realkind*( dptemp(i,2)+dptemp(i,3) )
          enddo

          do i=1,klon-1
             utemp(i) = u(i,j,k)*hyu(i,j)*0.5_realkind*( dptemp(i,2)+dptemp(i+1,2) )
          enddo

          if ((abs(afull(k))>1.e-14_realkind).and.(abs(bfull(k))>1.e-14_realkind)) then             !     hybrid surface
             do  i=1,klon
                zpm3 = ahalf(k  ) + bhalf(k  )*ps(i,j+1)
                zpp3 = ahalf(k+1) + bhalf(k+1)*ps(i,j+1)
                ptemp(i,3) = ( zpp3* log(zpp3)-zpm3*log(zpm3))/dptemp(i,3)
             enddo
          elseif((abs(afull(k))<1.e-14_realkind).and.(abs(bfull(k))>1.e-14_realkind)) then             !     sigma surface
             do i=1,klon
                ptemp(i,3) =  log(ps(i,j+1))
             enddo
          endif

          if (abs(bfull(k))<1.e-14_realkind) then              !     pressure surface
             do i=2,klon-1
                om1(i,j,k) = 0._realkind
             enddo
          else             !     hybrid or sigma surface
             do  i=2,klon-1
                om1(i,j,k) =  (1._realkind/dptemp(i,2))*              &
                     ( afull(k)+bfull(k)*ps(i,j) )*           &
                     ra*rhxu(i,j)*rhyv(i,j)*0.5_realkind*              &
                     ( rdlam*( utemp(i-1)*                          &
                     (ptemp(i  ,2)-ptemp(i-1,2))                  &
                     + utemp(i  )*(ptemp(i+1,2)-ptemp(i,2)))    &
                     + rdth* ( vtemp(i,1)*(ptemp(i  ,2)-          &
                     ptemp(i  ,1))                                  &
                     + vtemp(i,2)*(ptemp(i,3)-ptemp(i  ,2))))
             enddo
          endif

          do  i=2,klon-1
             om2(i,j,k+1) = om2(i,j,k) - ra*rhxu(i,j)*rhyv(i,j)*     &
                  ( rdlam*(utemp(i  )-utemp(i-1))  + rdth *(vtemp(i,2)-vtemp(i,1)) )
          enddo

          if (j<klat-1) then             !     prepare dptemp, ptemp and vtemp for next latitude
             do i=1,klon
                dptemp(i,2) = dptemp(i,3)
                ptemp(i,1) =  ptemp(i,2)
                ptemp(i,2) =  ptemp(i,3)
                vtemp(i,1) =  vtemp(i,2)
             enddo
          endif
       enddo
    enddo

    return
  end subroutine omcomp



end module sl2tim4
