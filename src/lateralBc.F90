module lateralBc

  use asimof
  use timetype
  use calendar
  use derived_types
  use boundaryRelaxation
  use domainmod
  use decomp
  use mod_grib,only:realgribkind
  use util
  implicit none
  private

  type(atm),save::bc(2)
  type(time),save::boundaryTime(2)
  type(deltatime),save::dtBoundary
  logical,allocatable,save,dimension(:,:)::lsv
  logical,save::lcw(2)
  public initLateral_bc,setLateral_bc,setInitialCondition,comped

contains

    subroutine comped(klon,klat,klev,pps,pu,pv,pedot,ahyb,bhyb)

    !     comped - compute vertical velocity
    !     
    !     j.e. haugen         hirlam
    !     k.s. eerola         hirlam 4 (2d arrays introduced)
    !     
    !     purpose.
    !     --------
    !     
    !     compute  pedot = pedot*d(p)/d(eta) / ps
    !     
    !     input parameters:
    !     -----------------
    !     
    !     klon      number of gridpoints in the x-direction
    !     klat      number of gridpoints in the y-direction
    !     klev      number of vertical levels
    !     pps       surface pressure
    !     pu        velocity in the x-direction
    !     pv        velocity in the x-direction
    !  
    !     output parameters:
    !     pedot     vertical velocity variable



    use rcaDomainMod
    use comhkp
    use decomp
    implicit none

    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(in)::pps(klon,klat),pu(klon,klat,klev),pv(klon,klat,klev)
    real(kind=realkind),intent(out):: pedot(klon,klat,klev+1)                              
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)

    !     declaration of local work-space

    integer:: klonm1,klatm1
    integer:: jk,jx,jy,ilevp1
    real(kind=realkind)::   zakm , zrdlah , zrdloh ,zdbk,zbkm, zdak
    real(kind=realkind):: zahxhy(klon,klat), zdpk  (klon,klat), zdpsdt(klon,klat), &
         zpkm  (klon,klat), zpkp  (klon,klat),                     &
         zdivk (klon,klat), zedpde(klon,klat), zuu   (klon,klat),  &
         zvv   (klon,klat)

    klonm1 = klon - 1
    klatm1 = klat - 1

!    zrdlo = rdlam
!    zrdla = rdth
!    zrear = 1.0/ra
    !     computation surface pressure tendency
    do jy = 1,klat
       do jx = 1,klon
          zuu   (jx,jy) = 0.0_realkind
          zvv   (jx,jy) = 0.0_realkind
          zahxhy(jx,jy) = rhxu(jx,jy)*rhyv(jx,jy)*ra
       enddo
    enddo

    do jk=1,klev
       zdak = ahyb(jk+1)-ahyb(jk)
       zdbk = bhyb(jk+1)-bhyb(jk)

       do jy = 1,klat
          do jx = 1,klon
             zdpk(jx,jy) = zdak+zdbk*pps(jx,jy)
          enddo
       enddo

       do jy = 1,klat-1
          do jx = 1,klon-1
             zuu(jx,jy)=zuu(jx,jy)-pu(jx,jy,jk)*(zdpk(jx,jy)+zdpk(jx+1,jy))
             zvv(jx,jy)=zvv(jx,jy)-pv(jx,jy,jk)*(zdpk(jx,jy)+zdpk(jx,jy+1))
          enddo
       enddo
    enddo

    zrdloh = rdlam*0.5_realkind
    zrdlah = rdth*0.5_realkind


    do jy = 2,klat-1
       do jx = 2,klon-1
          zdpsdt(jx,jy) = zahxhy(jx,jy)* &
               (( zuu(jx,jy)*hyu(jx,jy) - zuu(jx-1,jy)*hyu(jx-1,jy) )*zrdloh &
               + ( zvv(jx,jy)*hxv(jx,jy) - zvv(jx,jy-1)*hxv(jx,jy-1) )*zrdlah)
       enddo
    enddo


    !     initiate vertical loop
    ilevp1 = klev+1

    do jy = 1,klat
       do jx = 1,klon
          zedpde(jx,jy) = 0.0_realkind
          zpkm  (jx,jy) = pps(jx,jy)
          pedot(jx,jy,1     ) = 0.0_realkind
          pedot(jx,jy,ilevp1) = 0.0_realkind
       enddo
    enddo


    !     vertical integration of the continuity equation
    do 1000 jk=klev,2,-1
       zakm  = ahyb(jk)
       zbkm  = bhyb(jk)
       zdbk  = bhyb(jk+1)-bhyb(jk)

       do jy = 1,klat
          do jx = 1,klon
             zpkp(jx,jy) = zpkm(jx,jy)
             zpkm(jx,jy) = zakm + zbkm*pps(jx,jy)
             zdpk(jx,jy) = zpkp(jx,jy) - zpkm(jx,jy)
          enddo
       enddo

       do jy = 1,klat
          do jx = 1,klon-1
             zuu(jx,jy)= 0.5_realkind*(zdpk(jx+1,jy)+zdpk(jx,jy))*pu(jx,jy,jk)*hyu(jx,jy)
          enddo
       enddo

       do jy = 1,klat-1
          do jx = 1,klon
             zvv(jx,jy)=0.5_realkind*(zdpk(jx,jy+1)+zdpk(jx,jy))*pv(jx,jy,jk)*hxv(jx,jy)
          enddo
       enddo


       !     compute divergence at level k


       do jy = 2,klat-1
          do jx = 2,klon-1
             zdivk(jx,jy)= zahxhy(jx,jy)*((zuu(jx,jy)-zuu(jx-1,jy))*rdlam &
                  + ( zvv(jx,jy)-zvv(jx  ,jy-1) ) * rdth )
          enddo
       enddo

       do jy = 2,klat-1
          do jx = 2,klon-1
             zedpde(jx,jy)   = zedpde(jx,jy) + zdbk*zdpsdt(jx,jy)+ zdivk(jx,jy) 
             pedot(jx,jy,jk) = zedpde(jx,jy) / pps(jx,jy)
          enddo
       enddo

       do jx=1,klon
          pedot(jx,   1,jk) = 0.0_realkind
          pedot(jx,klat,jk) = 0.0_realkind
       enddo

       do jy=2,klat-1
          pedot(1   ,jy,jk) = 0.0_realkind
          pedot(klon,jy,jk) = 0.0_realkind
       enddo
1000 enddo
    call swapklevp1(pedot,klon,klat,klev+1)
    return
  end subroutine comped

  subroutine setInitialCondition(phi,nlsimp,nlslan)
    implicit none
    type(atm),intent(inout)::phi
    logical,intent(in)::nlsimp,nlslan
    real(kind=realkind)::timrat
    
    if(phi%timestamp>=bc(1)%timestamp.and. &
         phi%timestamp<=bc(2)%timestamp)then
       timrat = real(Nseconds(phi%timestamp-bc(1)%timestamp),realkind)/ &
            real(Nseconds(dtBoundary),realkind)
       if(timrat<0._realkind.or.timrat>1._realkind)then
          call stop_program( 'timrat OutOfBounds lateralBc')
       endif
       if(nlsimp)then
          phi%lnps = bc(1)%lnps+ timrat*(bc(2)%lnps-bc(1)%lnps)
          phi%ps = exp(phi%lnps)
       else
          phi%ps = bc(1)%ps+ timrat*(bc(2)%ps-bc(1)%ps) 
          phi%lnps = log(phi%ps)
       endif
       phi%u = bc(1)%u + timrat*(bc(2)%u -bc(1)%u ) 
       phi%v = bc(1)%v + timrat*(bc(2)%v -bc(1)%v ) 
       phi%t = bc(1)%t + timrat*(bc(2)%t -bc(1)%t ) 
       phi%q = bc(1)%q + timrat*(bc(2)%q -bc(1)%q ) 
       phi%cw = bc(1)%cw + timrat*(bc(2)%cw -bc(1)%cw ) 
       if(nlslan) then
          phi%edot = bc(1)%edot+timrat*(bc(2)%edot-bc(1)%edot)
       endif
       phi%svar = bc(1)%svar + timrat*(bc(2)%svar -bc(1)%svar ) 
    else
       call stop_program( 'Initial data cannot be interpolated from boundary data')
    endif
  end subroutine setInitialCondition

  subroutine setLateral_bc(klon,klat,klev,ksvar,phip,nlslan,RCAdom)

    use comhkp,only:nwmosv,nlhumc,nlsimp
    use physpar,only:noption
    implicit none
    integer,intent(in)::klon,klat,klev,ksvar
    type(atm),intent(inout)::phip
    logical,intent(in)::nlslan


    type(domain),intent(in)::RCAdom

    do while(phip%timestamp>bc(2)%timestamp) !if time is enough read new data
       !a loop statement if the timestep exceeds dtBoundary
       bc(1) = bc(2)
       lcw(1) = lcw(2)
       lsv(:,1) = lsv(:,2)
       boundaryTime(1) = boundaryTime(2)
       boundaryTime(2) = boundaryTime(2) + dtBoundary

       call getlateral_bc(klon,klat,klev,ksvar,nwmosv, &
       boundaryTime(2), &
       bc(2), &
       nlhumc,noption(2)>0,&
       RCAdom,&
       lsv(:,2),lcw(2))
    enddo
    call bndrel(klon,klat,klev,ksvar,phip,bc,nlsimp,nlslan,&
         noption(2)==1,lcw,lsv )
  end subroutine setLateral_bc


  subroutine initLateral_bc(startTime,klon,klat,klev,ksvar,&
       nwmosv,nlhumc,RCAdom)
    use physpar, only:noption,read_namprc
    implicit none
    type(time),intent(in)::startTime
    integer,intent(in)::klon,klat,klev,ksvar
    integer,intent(in)::nwmosv(ksvar)

    logical,intent(in)::nlhumc   
    type(domain),intent(in)::RCAdom
    integer::k,dtsec(4),secmin,kmin
    allocate(lsv(ksvar,2))
    call initBoundary(klon,klat,klev,RCAdom%ahyb,RCAdom%bhyb,startTime)

    dtBoundary = makedeltatime(0,0,0,6,0,0)
    call read_namprc()

    boundaryTime(1) = startTime 
    
    if(boundaryTime(1)%hour/= 0 .and.boundaryTime(1)%hour/= 6 .and. &
         boundaryTime(1)%hour/= 12 .and.  boundaryTime(1)%hour/= 18)then
       do k=1,4
          dtsec(k) = startTime%hour - (k-1)*6
       enddo
       !We need to pick out the smallest positive (or zero) dtsec, since the time interpolation assumes data(1) to be 
       !at startTime or before
       secmin = -99999
       kmin = -1
       do k=1,4
          if(dtsec(k)>=0.and.secmin<dtsec(k))then
             secmin = dtsec(k)
             kmin = k
          endif
       enddo
       if(kmin==-1)call stop_program('Error init lateral_bc kmin==-1')
       boundaryTime(1) = maketime(startTime%year,startTime%month,startTime%day,(kmin-1)*6,0,0)
       
    endif


    boundaryTime(2) = boundaryTime(1) + dtBoundary
    bc(1) = alloc_atm(klon,klat,klev,ksvar,boundaryTime(1))
    bc(2) = alloc_atm(klon,klat,klev,ksvar,boundaryTime(2))

    call getlateral_bc(klon,klat,klev,ksvar,nwmosv, &
       boundaryTime(1),&
       bc(1), &
       nlhumc,noption(2)>0, &
       RCAdom,&
       lsv(:,1),lcw(1))
    
    call getlateral_bc(klon,klat,klev,ksvar,nwmosv, &
       boundaryTime(2), &
       bc(2), &
       nlhumc,noption(2)>0,&
       RCAdom,&
       lsv(:,2),lcw(2))

  end subroutine initLateral_bc



    
  subroutine getlateral_bc(klon,klat,klev,ksvar,kwmosv,&
       ktime, &
       bc, &
       lhumc,lrdcw,RCAdom,&
       lsv,lcw)

    use referenceParameters, only:sip0
    use escom 
    use confys
    use gcm
    use modddr

    implicit none

    integer,intent(in):: klon, klat,  klev, ksvar 
    integer,intent(in)::kwmosv(ksvar)
    type(time),intent(in)::ktime

    type(atm),intent(inout)::bc
    logical,intent(inout)::lsv(ksvar),lcw
    type(domain),intent(in)::RCAdom
    logical,intent(in)::lhumc   
    logical,intent(in)::lrdcw    ! indicator whether to read cloud water
    logical::losv(ksvar)       ! indicator of extra scalars in file
    logical::locw               ! indicator of cloud water in file
    integer::l
    logical,save::debug=.false.
    type(time)::filetime
    type(ddr)::gcmddr

    call pre_getgrb(klon,klat, klev,ksvar,91, kwmosv,     &
          bc%ps, bc%u   , bc%v   , bc%t   , bc%q  , bc%cw   , bc%svar, &
          ktime, & 
          RCAdom, &
          losv, lrdcw, locw,gcmddr)

    !use ddr to check that the time is right

    !Patrick is supisious against the adding of lhadgem here.
    !Should probabaly be checked via calendar...
    if(lhadley.or.lecham2.or.lhadgem)then
       gcmddr%ndtvhl=ktime%year*10000+ktime%month*100+ktime%day
       gcmddr%nscvhl=ktime%hour*10000
    endif

    filetime%year  = gcmddr%ndtvhl/10000
    filetime%month = gcmddr%ndtvhl/100-ktime%year*100
    filetime%day   = gcmddr%ndtvhl-ktime%year*10000-ktime%month*100
    filetime%hour  = gcmddr%nscvhl/10000
    
    if(ktime/=filetime)then
       call stop_program( 'Wrong date in demanding data getlateral_bc')
    endif

    if(mype==0.and.debug)then
       write(6,'(1x,''date/time for boundary data:'',i8.8,a,i6.6)') &
            gcmddr%ndtvhl,'/',gcmddr%nscvhl
    endif
    
    
    !          compute ln(ps)
    bc%lnps = log( bc%ps )
    !          compute vertical velocity in *comped*
    call comped(klon,klat,klev,bc%ps,bc%u,bc%v,bc%edot,RCAdom%ahyb,RCAdom%bhyb)

    !          check humidity for critical values in *crihum*
    if (lhumc) then
       call crihum(klon,klat,klev,RCAdom%ahyb,RCAdom%bhyb,bc%ps,bc%t,bc%q)
    endif

    !     estimate initial cloud condensate in boundary zone
    !     when cloud water is not in boundary file ( relevant
    !     only initially if cloud water is not
    !    boundary relaxed later during forecast )
    if( .not.locw ) then
       call estimateCloudCondensate(klon,klat,klev,RCAdom,bc)
    endif

    !     IS CLOUD-WATER IN THE BOUNDARY DATA FILE?

    lcw=locw

    !     IF SEMI-LAGRANGIAN SCHEME:
    !     1. SET BOUNDARY-VALUES TO ZERO IF NOT IN INPUT FILE.
    !     2. DO BOUNDARY RELAXATION.

    if (.not.lcw) then
       bc%cw = 0._realkind
    endif
    lcw=.true.

    !     ARE THE EXTRA SCALARS IN THE BOUNDARY FILE?
    do l=1,ksvar
       lsv(l)=losv(l)
       if (.not.lsv(l)) then
          bc%svar(:,:,:,l)=0.01_realkind
       endif
       lsv(l)=.true.
    enddo
    bc%timestamp = ktime !set the b.c. to be at time of input/query_time
    return
  end subroutine getlateral_bc


  subroutine pre_getgrb(klon,klat, klev,ksvar,lundir,nwmosv, &
       ps,u,v,t,q,cw,svar,&
       ktime,&
       RCAdom,&
       lasv,lrdcw,lacw,gcmddr)       
#ifdef USING_NETCDF
    use netcdfrestart
#else
    use restart
#endif
    use gcm
    use interpol
    use IO
    use modddr

    implicit none

    type(ddr),intent(out)::gcmddr
    integer,intent(in):: klev,klon,klat
    integer,intent(in):: ksvar

    real(kind=realkind),intent(inout)::ps(klon,klat),    &
         u(klon,klat,klev),   v(klon,klat,klev), &
         t(klon,klat,klev),   q(klon,klat,klev), &
         cw(klon,klat,klev)

    real(kind=realkind),intent(inout)::svar(klon,klat,klev,ksvar)
    integer,intent(in)::lundir
    logical,intent(in)::lrdcw
    logical,intent(inout)::lacw
    type(domain),intent(in)::RCAdom
    type(time),intent(in)::ktime
    logical::lasv(ksvar)
    integer::nwmosv(ksvar)
    real(kind=realkind)::north,east
    logical:: lfis
    integer, parameter::jpbuf=700000
    real(kind=realgribkind):: buf(jpbuf)
    integer:: ilen,iy,im,id,ih,imin,length
    character(len=132):: prefix,filenam 

    north = RCAdom%south + real(klat_global-1,realkind)*RCAdom%dlat
    east  = RCAdom%west  + real(klon_global-1,realkind)*RCAdom%dlon


    call readboundpath(prefix,ilen)

    iy=ktime%year
    im=ktime%month
    id=ktime%day
    ih=ktime%hour
    imin=0
    length=0

!!$    if(lecmwf.or.lccsm)then
!!$       call cre_file_ma(prefix,ilen,iy,im,id,ih,length,filenam)
!!$    elseif(lhadley.or.lecham.or.lecham2.or.lcanesm2.or.lcnrm)then
!!$       call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
!!$    endif

    call createFilename(prefix,ilen,iy,im,id,ih,length,filenam)
    call groploc(lundir,filenam,buf,jpbuf,gcmddr)

    if(lhadley)then
       lfis = .false.
    else
!uh!
       lfis = .true.       
    endif
!uh next not necessary?
    if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl) then
       lfis = .true.       
    endif

    call interpol_bd(klon,klat,klev,lrdcw,lacw,lundir, &
         RCAdom,&
         lfis, &
         ps,u,v,t,q,cw,svar, &
         gcmdomain,&
         ksvar, &
         nwmosv,&
         lasv)

    call grclloc(lundir,mype)

    return
  end subroutine pre_getgrb








  subroutine interpol_bd(klon,klat,klev,lrdcw,lacw,lundir,&
       RCAdom,&
       lfis, &
       ps,u,v,t,q,cw,svar,&
       gcmdomainLoc,&
       ksvar, &
       nwmosv, &
       lasv)

    use util
    use grw1
    use IO
    use ctun, only: tseafr
    use confys
    use gcm
    use routmodvar
    use config
    use referenceParameters, only:sip0
    implicit none

    integer,intent(in)::klon,klat,klev,ksvar,lundir
    type(domain),intent(in)::RCAdom,gcmdomainLoc

    logical,intent(in):: lfis

    real(kind=realkind),intent(inout):: ps(klon,klat),                  &
         u(klon,klat,klev),              &
         v(klon,klat,klev),              &
         t(klon,klat,klev),              &
         q(klon,klat,klev)
    real(kind=realkind):: cw(klon,klat,klev)
    real(kind=realkind) graviti
    real(kind=realkind),allocatable:: psrel(:,:) 

    logical,intent(in):: lrdcw
    logical::lacw
    real(kind=realkind) dlat_gribi,dlon_gribi
    integer::nwmosv(ksvar)
    real(kind=realkind)::svar(klon,klat,klev,ksvar)
    logical:: lasv(ksvar)
    real(kind=realkind),allocatable:: fis_grib(:,:)
    real(kind=realkind),allocatable:: ps_grib(:,:),uu(:,:,:),vv(:,:,:),uu_grib(:,:,:), &
         vu_grib(:,:,:),uv_grib(:,:,:),vv_grib(:,:,:),t_grib(:,:,:),q_grib(:,:,:)
    real(kind=realkind),allocatable::lat_m(:,:),lon_m(:,:),lat_u(:,:),lon_u(:,:),lat_v(:,:),lon_v(:,:)
    real(kind=realkind),allocatable,dimension(:,:)::zlat1_m,zlon1_m,zlat1_u,zlon1_u,zlat1_v,zlon1_v
    real(kind=realkind),allocatable,dimension(:,:)::zlat2_m,zlon2_m,zlat2_u,zlon2_u,zlat2_v,zlon2_v
    
    integer,allocatable:: i_m(:,:),j_m(:,:),i_u(:,:),j_u(:,:),i_v(:,:),j_v(:,:) 
    real(kind=realkind),allocatable:: w_x_m(:,:),w_y_m(:,:),w_x_u(:,:),w_y_u(:,:),w_x_v(:,:),w_y_v(:,:)
    real(kind=realkind),allocatable::  pa(:,:),pb(:,:),pc(:,:),pd(:,:), fi500ec(:,:)         

    
    integer::i,j,itype,ipar,ilev,jlev,kcall
    integer::ierr

    real(kind=realkind)::eps,x_m,y_m,y_u,x_u,y_v,x_v

    logical::rotate
    integer::lev_start
    real(kind=realkind)::tempo(klon,klat)
    real(kind=realkind)::td1echam(klon,klat)
    real(kind=realkind)::td2echam(klon,klat)
    real(kind=realkind)::td3echam(klon,klat)
    real(kind=realkind)::td4echam(klon,klat)
    real(kind=realkind)::dist,w1,w2,dec1,dec2,dec3,dec4,drc3,drc4,drc5

!cau1106
    real(kind=realkind)::plev,rdcp,sav, rmin,rmax

    data     eps/1.e-6_realkind/



    allocate(psrel(klon,klat))
    allocate(fis_grib(klon,klat))
    allocate(ps_grib(klon,klat))
    allocate(uu(klon,klat,gcmdomainLoc%nlev))
    allocate(vv(klon,klat,gcmdomainLoc%nlev))
    allocate(uu_grib(klon,klat,gcmdomainLoc%nlev))
    allocate(vu_grib(klon,klat,gcmdomainLoc%nlev))
    allocate(uv_grib(klon,klat,gcmdomainLoc%nlev))
    allocate(vv_grib(klon,klat,gcmdomainLoc%nlev))
    allocate(t_grib(klon,klat,gcmdomainLoc%nlev))
    allocate(q_grib(klon,klat,gcmdomainLoc%nlev))

    allocate(lat_m(klon,klat), &
         lon_m(klon,klat),      &
         lat_u(klon,klat),      &
         lon_u(klon,klat),      &
         lat_v(klon,klat),      &
         lon_v(klon,klat),      &
         zlat1_m(klon,klat),    &
         zlon1_m(klon,klat),    &
         zlat1_u(klon,klat),    &
         zlon1_u(klon,klat),    &
         zlat1_v(klon,klat),    &
         zlon1_v(klon,klat),    &
         zlat2_m(klon,klat),    &
         zlon2_m(klon,klat),    &
         zlat2_u(klon,klat),    &
         zlon2_u(klon,klat),    &
         zlat2_v(klon,klat),    &
         zlon2_v(klon,klat))

    allocate(i_m(klon,klat),   &
         j_m(klon,klat),        &
         i_u(klon,klat),        &
         j_u(klon,klat),        &
         i_v(klon,klat),        &
         j_v(klon,klat))              

    allocate(w_x_m(klon,klat),  &
         w_y_m(klon,klat),      &
         w_x_u(klon,klat),      &
         w_y_u(klon,klat),      &
         w_x_v(klon,klat),      &
         w_y_v(klon,klat))            

    allocate(pa(klon,klat),     &
         pb(klon,klat),         &
         pc(klon,klat),         &
         pd(klon,klat),         &
         fi500ec(klon,klat))


    !     initiate flags for presence of data

    lacw  =.true.

    do j=1,ksvar
       lasv(j)=.true.
    enddo


    !     1.1   calculate grid-point coordinates of non-staggered output field
    do j=1,klat
       lat_m(:,j) = RCAdom%south + real( jdatastart + j - 2,realkind )*RCAdom%dlat
    enddo
    do i=1,klon
       lon_m(i,:) = RCAdom%west + real( idatastart + i - 2,realkind  )*RCAdom%dlon
    enddo

    !     add staggering to wind component coordinates

    do j=1,klat
       do i=1,klon
          lat_u(i,j) = lat_m(i,j) + RCAdom%stagu_y*RCAdom%dlat
          lon_u(i,j) = lon_m(i,j) + RCAdom%stagu_x*RCAdom%dlon
          lat_v(i,j) = lat_m(i,j) + RCAdom%stagv_y*RCAdom%dlat
          lon_v(i,j) = lon_m(i,j) + RCAdom%stagv_x*RCAdom%dlon
       enddo
    enddo

    !     rotate coordinates if needed


    if(abs(RCAdom%polon-gcmdomainLoc%polon)>eps.or. &
         abs(RCAdom%polat-gcmdomainLoc%polat)>eps)then
       rotate = .true.

       !     de-rotate if needed coordinates from output grid geometry
       if( abs(RCAdom%polat+90._realkind)>eps ) then
          call regrot(zlon1_m,zlat1_m,lon_m,lat_m, &
               klon,klat,&
               RCAdom%polon,RCAdom%polat,-1)                        
          call  regrot(zlon1_u,zlat1_u,lon_u,lat_u,  &
               klon,klat,&
               RCAdom%polon,RCAdom%polat,-1)                        
          call  regrot(zlon1_v,zlat1_v,lon_v,lat_v,  &
               klon,klat,&
               RCAdom%polon,RCAdom%polat,-1)
       else
          do j=1,klat
             do i=1,klon
                zlat1_m(i,j) = lat_m(i,j)
                zlon1_m(i,j) = lon_m(i,j)
                zlat1_u(i,j) = lat_u(i,j)
                zlon1_u(i,j) = lon_u(i,j)
                zlat1_v(i,j) = lat_v(i,j)
                zlon1_v(i,j) = lon_v(i,j)
             enddo
          enddo
       endif

       !     1.3.2  rotate if needed coordinates to input grid geometry 

       if( abs(gcmdomainLoc%polat+90._realkind)>eps ) then
          call  regrot(zlon1_m,zlat1_m,zlon2_m,zlat2_m, &
               klon,klat,&
               gcmdomainLoc%polon,gcmdomainLoc%polat,+1)                  
          call  regrot(zlon1_u,zlat1_u,zlon2_u,zlat2_u,  &
               klon,klat,&
               gcmdomainLoc%polon,gcmdomainLoc%polat,+1)                  
          call  regrot(zlon1_v,zlat1_v,zlon2_v,zlat2_v,  &
               klon,klat,&
               gcmdomainLoc%polon,gcmdomainLoc%polat,+1)
       else
          do j=1,klat
             do i=1,klon
                zlat2_m(i,j) = zlat1_m(i,j)
                zlon2_m(i,j) = zlon1_m(i,j)
                zlat2_u(i,j) = zlat1_u(i,j)
                zlon2_u(i,j) = zlon1_u(i,j)
                zlat2_v(i,j) = zlat1_v(i,j)
                zlon2_v(i,j) = zlon1_v(i,j)
             enddo
          enddo
       endif
       !     no rotation needed, just copy coordinates
    else
       rotate = .false.
       do j=1,klat
          do i=1,klon
             zlat2_m(i,j) = lat_m(i,j)
             zlon2_m(i,j) = lon_m(i,j)
             zlat2_u(i,j) = lat_u(i,j)
             zlon2_u(i,j) = lon_u(i,j)
             zlat2_v(i,j) = lat_v(i,j)
             zlon2_v(i,j) = lon_v(i,j)
          enddo
       enddo
    endif


    ! calculate indeces and weights needed for horizontal interpolation

    dlat_gribi = 1.0_realkind/gcmdomainLoc%dlat
    dlon_gribi = 1.0_realkind/gcmdomainLoc%dlon
    do j=1,klat
       do i=1,klon
          if( zlon2_m(i,j)<gcmdomainLoc%west )then
             zlon2_m(i,j) = zlon2_m(i,j) + 360._realkind      
          endif
          if( zlon2_u(i,j)<gcmdomainLoc%west )then
             zlon2_u(i,j) = zlon2_u(i,j) + 360._realkind     
          endif
          if( zlon2_v(i,j)<gcmdomainLoc%west )then            
             zlon2_v(i,j) = zlon2_v(i,j) + 360._realkind
          endif

          y_m = ( zlat2_m(i,j) - gcmdomainLoc%south )*dlat_gribi + 1.0_realkind
          j_m(i,j)    = int( y_m )
          w_y_m(i,j)  = y_m - real(j_m(i,j),realkind)

          x_m = ( zlon2_m(i,j) - gcmdomainLoc%west )*dlon_gribi + 1.0_realkind
          i_m(i,j)    = int( x_m )

          if ( i_m(i,j) > gcmdomainLoc%nlon ) then
             x_m = x_m - real(gcmdomainLoc%nlon,realkind)
             i_m(i,j) = i_m(i,j) - gcmdomainLoc%nlon
          endif
          w_x_m(i,j)  = x_m - real(i_m(i,j),realkind)

          y_u = ( zlat2_u(i,j) - gcmdomainLoc%south )*dlat_gribi - gcmdomainLoc%stagu_y + 1.0_realkind
          j_u(i,j)    = int( y_u )
          w_y_u(i,j)  = y_u - real(j_u(i,j),realkind)
          x_u = ( zlon2_u(i,j) - gcmdomainLoc%west )*dlon_gribi - gcmdomainLoc%stagu_x + 1.0_realkind
          i_u(i,j)    = int( x_u )
          if ( i_u(i,j) > gcmdomainLoc%nlon ) then
             x_u = x_u - real(gcmdomainLoc%nlon,realkind)
             i_u(i,j) = i_u(i,j) - gcmdomainLoc%nlon
          endif
          w_x_u(i,j)  = x_u - real(i_u(i,j),realkind)

          y_v = ( zlat2_v(i,j) - gcmdomainLoc%south )*dlat_gribi - gcmdomainLoc%stagv_y + 1.0_realkind
          j_v(i,j)    = int( y_v )
          w_y_v(i,j)  = y_v - real(j_v(i,j),realkind)
          x_v = ( zlon2_v(i,j) - gcmdomainLoc%west )*dlon_gribi - gcmdomainLoc%stagv_x + 1.0_realkind
          i_v(i,j)    = int( x_v )
          if ( i_v(i,j) > gcmdomainLoc%nlon ) then
             x_v = x_v - real(gcmdomainLoc%nlon,realkind)
             i_v(i,j) = i_v(i,j) - gcmdomainLoc%nlon
          endif
          w_x_v(i,j)  = x_v - real(i_v(i,j),realkind)
       enddo
    enddo



    !     read fis for the global model e.g. era40
    if(lfis) then
       itype = 105
       ipar  = 006
       ilev  = 000
       if(gcmdomainLoc%dlon<=0.77)itype = 001 !THIS IS A HACK (MARCO)

!cau1106
       if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
          itype = 001
          ipar  = 006
          ilev  = 000
       endif
!       if(lhadley.or.lcanesm2) then
!          do j=1,klat
!             do i=1,klon
!                fis_grib(i,j) = orog_bound(i,j)
!             end do
!          end do
!          ierr=0
!       else
          call grrdloc_hint(lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
            fis_grib,fis_grib,                                 &
            klon,klat,                                         &
            gcmdomainLoc%nlon,gcmdomainLoc%nlat,               &
            i_m,w_x_m,j_m,w_y_m,                               &
            .false.,                                           &
            i_m,w_x_m,j_m,w_y_m,                               &
            ierr)
!       endif

!cau1190623
          if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
             rmin=1.e10
             rmax=0.
             do j=1,klat
                do i=1,klon
                   fis_grib(i,j)=fis_grib(i,j)*9.81
                   rmin=min(rmin,fis_grib(i,j))
                   rmax=max(rmax,fis_grib(i,j))
                enddo
             enddo
!             write(6,*)' orog minmax ',mype,rmin,rmax
          endif
          
          if (ierr/=0) then
             write(6,*) ' lecmwf ',lecmwf
             write(6,*) ' klon,klat = ',klon,klat
             write(6,*) ' nx_grib,ny_grib = ',gcmdomainLoc%nlon,gcmdomainLoc%nlat
             write(6,*)'error with fis first ',ipar,itype,ilev,ierr
             call stop_program('')
          endif
       endif                     !  lfis
    

    !     read ps
    if( lecmwf ) then
       itype = 109
       ipar  = 152
       ilev  = 1
    elseif(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
       itype = 001
       ipar  = 001
       ilev  = 000
    else
       itype = 105
       ipar  = 001
       ilev  = 000
    endif
    if(lccsm)then
       itype = 105
       ipar  = 001
       ilev  = 000
    endif

    call grrdloc_hint(                         &
         lundir,itype,ipar,real(ilev,realgribkind),lecmwf,   &
         ps_grib,ps_grib,                       &
         klon,klat,                             &
         gcmdomainLoc%nlon,gcmdomainLoc%nlat,                       &
         i_m,w_x_m,j_m,w_y_m,                   &
         .false.,                               &
         i_m,w_x_m,j_m,w_y_m,                   &
         ierr)
    if(ierr/=0)call stop_program( 'could not read ps')

    if( lecmwf ) then
       do j=1,klat
          do i=1,klon
             ps_grib(i,j) = exp( ps_grib(i,j) )
          enddo
       enddo
    endif



    !     read model level fields

    itype = 109
    lev_start = 1
    if(lecmwf)then
       lev_start = 13
    endif

    rdcp=287./1004.

    uu_grib = 0._realgribkind
    vv_grib = 0._realgribkind
    uv_grib = 0._realgribkind
    vu_grib = 0._realgribkind


    do ilev=lev_start,gcmdomainLoc%nlev

       !     u-component

       ipar  = 033
       call grrdloc_hint(                        &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf,  &
            uu_grib(:,:,ilev),uv_grib(:,:,ilev),  &
            klon,klat,                            &
            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                      &
            i_u,w_x_u,j_u,w_y_u,                  &
            rotate,                               &
            i_v,w_x_v,j_v,w_y_v,                  &
            ierr)
       if(ierr/=0)call stop_program( 'could not read u')


       !     v-component

       ipar  = 034
       call grrdloc_hint(                             &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf,       &
            vv_grib(:,:,ilev),vu_grib(:,:,ilev),       &
            klon,klat,                                 &
            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                           &
            i_v,w_x_v,j_v,w_y_v,                       &
            rotate,                                    &
            i_u,w_x_u,j_u,w_y_u,                       &
            ierr)                                            
       if(ierr/=0)call stop_program( 'could not read v')

       !    temperature                                         
       ipar  = 011                                      
       call grrdloc_hint(                              &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf,       &
            t_grib(:,:,ilev),t_grib(:,:,ilev),         &
            klon,klat,                                 &
            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                           &
            i_m,w_x_m,j_m,w_y_m,                       &
            .false.,                                   &
            i_m,w_x_m,j_m,w_y_m,                       &
            ierr)
       if(ierr/=0)call stop_program( 'could not read T')
       

       if(lhadley) then
          do j=1,klat
             do i=1,klon
                plev=gcmdomain%afull(ilev)+gcmdomain%bfull(ilev)*ps_grib(i,j)
                t_grib(i,j,ilev)=t_grib(i,j,ilev)*(plev/100000.)**rdcp
             enddo
          enddo
          if(mype.eq.0.and.ilev.eq.1) write(6,*) ' THETA to Temp'
       endif

       !     humidity
       ipar  = 051
       call grrdloc_hint(                       &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
            q_grib(:,:,ilev),q_grib(:,:,ilev),   &
            klon,klat,                           &
            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                     &
            i_m,w_x_m,j_m,w_y_m,                 &
            .false.,                             &
            i_m,w_x_m,j_m,w_y_m,                 &
            ierr)
       if(ierr/=0)call stop_program( 'could not read q')

       !     cloud water
       !     note cw not vertically interpolated 
       !     only used when nlev == gcmdomainLoc%nlev
       if( lrdcw ) then
          if(klev==gcmdomainLoc%nlev.and. .not.lmiroc5)then !a special case??
             if(lecmwf)then
                ipar  = 247
             else
                ipar  = 058
             endif
             call grrdloc_hint(lundir,itype,ipar,real(ilev,realgribkind),lecmwf,     &
                  fi500ec,fi500ec,                         &
                  klon,klat,                               &
                  gcmdomainLoc%nlon,gcmdomainLoc%nlat,                         &
                  i_m,w_x_m,j_m,w_y_m,                     &
                  .false.,                                 &
                  i_m,w_x_m,j_m,w_y_m,                     &
                  ierr)                                     
             if (ierr/=0) then                            
                write(6,*)' cloud error  ',ipar,ierr        
                call stop_program( 'could not read fi500ec')
             endif

             if(lecmwf)then                                 
                ipar  = 246                                 
             else                                           
                ipar  = 076                                 
             endif
             call grrdloc_hint(                            &
                  lundir,itype,ipar,real(ilev,realgribkind),lecmwf,     &
                  cw(:,:,ilev),cw(:,:,ilev),                 &
                  klon,klat,                               &
                  gcmdomainLoc%nlon,gcmdomainLoc%nlat,                         &
                  i_m,w_x_m,j_m,w_y_m,                     &
                  .false.,                                 &
                  i_m,w_x_m,j_m,w_y_m,                     &
                  ierr)
             if(ierr/=0) then
                write(6,*)' cloud error  ',ipar,ierr
                call stop_program( 'could not read cw')
             endif
             do j=1,klat
                do i=1,klon
                   cw(i,j,ilev) = cw(i,j,ilev) + fi500ec(i,j)
                enddo
             enddo
             if(ierr/=0) lacw=.false.
          else
             lacw = .false.
          endif               ! klev==gcmdomainLoc%nlev
       else
          do j=1,klat
             do i=1,klon
                cw(i,j,ilev) = 0._realkind
             enddo
          enddo
       endif                  ! lrdcw





       !     extra scalars
       !     note svar not vertically interpolated 
       !     only used when klev == gcmdomainLoc%nlev
       do j=1,ksvar
          if(lasv(j))then
             if(klev==gcmdomainLoc%nlev .and..not.lmiroc5) then
                ipar = nwmosv(j)
                call grrdloc_hint(                       &
                     lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
                     svar(:,:,ilev,j),                    &
                     svar(:,:,ilev,j),                    &
                     klon,klat,                           &
                     gcmdomainLoc%nlon,gcmdomainLoc%nlat,                     &
                     i_m,w_x_m,j_m,w_y_m,                 &
                     .false.,                             &
                     i_m,w_x_m,j_m,w_y_m,                 &
                     ierr)

                if(ierr/=0 )then
                   call stop_program( 'could not read tke??')
                   lasv(j) = .false.
                endif
             else
                lasv(j) = .false.
             endif            ! klev==gcmdomainLoc%nlev
          endif               ! lasv(j)
       enddo
    enddo                     ! gcmdomainLoc%nlev
!!$    if(.not.lacw.and.mype==0)then
!!$       print *,'SOMEBODY should check this and remove if not needed'
!!$       print *,' failure to find cw   on file ',lundir
!!$    endif
!!$    do j=1,ksvar
!!$       if(.not.lasv(j).and.mype==0)then
!!$          print *,'SOMEBODY should check this and remove if not needed'
!!$          print *,' failure to find TKE on file ',lundir
!!$       endif
!!$    enddo


!       if (lcanesm2) then
!          do ilev = 1,gcmdomainLoc%nlev/2
!             do j=1,klat
!                do i=1,klon
!                   sav=uu_grib(i,j,ilev)
!                   uu_grib(i,j,ilev)=uu_grib(i,j,gcmdomainLoc%nlev-ilev+1)
!                   uu_grib(i,j,gcmdomainLoc%nlev-ilev+1)=sav
!                   sav=vv_grib(i,j,ilev)
!                   vv_grib(i,j,ilev)=vv_grib(i,j,gcmdomainLoc%nlev-ilev+1)
!                   vv_grib(i,j,gcmdomainLoc%nlev-ilev+1)=sav
!                   sav=t_grib(i,j,ilev)
!                   t_grib(i,j,ilev)=t_grib(i,j,gcmdomainLoc%nlev-ilev+1)
!                   t_grib(i,j,gcmdomainLoc%nlev-ilev+1)=sav
!                   sav=q_grib(i,j,ilev)
!                   q_grib(i,j,ilev)=q_grib(i,j,gcmdomainLoc%nlev-ilev+1)
!                   q_grib(i,j,gcmdomainLoc%nlev-ilev+1)=sav
!                enddo
!             enddo
!          enddo
!          if(mype.eq.0) write(6,*) 'canesm2 turn upside down'
!       endif

    !     read tsd
    if(lecham2)then
       if(lecham5)then
          itype = 111
          ipar  = 85
          ilev  = 3
          call grrdloc_hint(                          &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,    &
               td1echam,td1echam,                      &
               klon,klat,                              &
               gcmdomainLoc%nlon,gcmdomainLoc%nlat,                        &
               i_m,w_x_m,j_m,w_y_m,                    &
               .false.,                                &
               i_m,w_x_m,j_m,w_y_m,                    &
               ierr)                                         
          if(ierr/=0)call stop_program( 'could not read td1echam')


          itype = 111                                   
          ipar  = 85                                    
          ilev  = 19                                    
          call grrdloc_hint(                           &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,    &
               td2echam,td2echam,                      &
               klon,klat,                              &
               gcmdomainLoc%nlon,gcmdomainLoc%nlat,                        &
               i_m,w_x_m,j_m,w_y_m,                    &
               .false.,                                &
               i_m,w_x_m,j_m,w_y_m,                    &
               ierr)
          if(ierr/=0)call stop_program( 'could not read td2echam')
       

          itype = 111
          ipar  = 85
          ilev  = 78

          call grrdloc_hint(                         &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,   &
               td3echam,td3echam,                     &
               klon,klat,                             &
               gcmdomainLoc%nlon,gcmdomainLoc%nlat,                       &
               i_m,w_x_m,j_m,w_y_m,                   &
               .false.,                               &
               i_m,w_x_m,j_m,w_y_m,                   &
               ierr)
          if(ierr/=0)call stop_program( 'could not read td3echam')


          itype = 111
          ipar  = 85
          ilev  = 268

          call grrdloc_hint(                       &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
               td4echam,td4echam,                   &
               klon,klat,                           &
               gcmdomainLoc%nlon,gcmdomainLoc%nlat,                     &
               i_m,w_x_m,j_m,w_y_m,                 &
               .false.,                             &
               i_m,w_x_m,j_m,w_y_m,                 &
               ierr)
          if(ierr/=0)call stop_program( 'could not read td4echam')



          !     depths (cm)
          dec1=3._realkind
          dec2=19._realkind
          dec3=78._realkind
          dec4=268._realkind
          drc3=17.7_realkind
          drc4=64.2_realkind
          drc5=194.7_realkind

          dist=1.0_realkind/(dec2-dec1)
          w1=abs((dec2-drc3))*dist
          w2=abs((dec1-drc3))*dist


          dist=1.0_realkind/(dec3-dec2)
          w1=abs((dec3-drc4))*dist
          w2=abs((dec2-drc4))*dist

          dist=1.0_realkind/(dec4-dec3)
          w1=abs((dec4-drc5))*dist
          w2=abs((dec3-drc5))*dist


          graviti = 0.0065_realkind/9.81_realkind


          !     i.e. echam4
          itype = 105
          ipar  = 85
          ilev  = 5

          call grrdloc_hint(                       &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
               tempo,tempo,                         &
               klon,klat,                           &
               gcmdomainLoc%nlon,gcmdomainLoc%nlat,                     &
               i_m,w_x_m,j_m,w_y_m,                 &
               .false.,                             &
               i_m,w_x_m,j_m,w_y_m,                 &
               ierr)                                 
          if(ierr/=0)call stop_program( 'could not read tempo')



       endif

    endif



   


!!$    !     for routing
!!$    if(use_routing) then
!!$
!!$       itype = 001
!!$       ipar  = 219
!!$       ilev  = 000
!!$
!!$       call grrdloc_global(                                        &
!!$            lundir,itype,ipar,real(ilev,realgribkind),                           &
!!$            runoff_bound,                                           &
!!$            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                                        &
!!$            klon,klat,                                              &
!!$            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                                        &
!!$            ierr)                                                    
!!$       if(ierr/=0)call stop_program( 'could not read runoff_bound')
!!$
!!$       rmin=minval(runoff_bound)                                     
!!$       rmax=maxval(runoff_bound)                                     
!!$       rmean=sum(runoff_bound)/(gcmdomainLoc%nlon*gcmdomainLoc%nlat)                          
!!$
!!$       if(mype==0)then                                             
!!$          write(6,661) ktime%hour,gcmdomainLoc%nlon,gcmdomainLoc%nlat,klon_bound,klat_bound,  &
!!$               rmin,rmax,rmean
!!$       endif
!!$
!!$       itype = 001
!!$       ipar  = 160
!!$       ilev  = 000
!!$       call grrdloc_global(                                         &
!!$            lundir,itype,ipar,real(ilev,realgribkind),                            &
!!$            runoff_bound,                                            &
!!$            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                                         &
!!$            klon,klat,                                               &
!!$            gcmdomainLoc%nlon,gcmdomainLoc%nlat,                                         &
!!$            ierr)                                                          
!!$       if(ierr/=0)call stop_program( 'could not read runoff_bound')
!!$       rmin=minval(runoff_bound)                                      
!!$       rmax=maxval(runoff_bound)                                      
!!$       rmean=sum(runoff_bound)/(gcmdomainLoc%nlon*gcmdomainLoc%nlat)                      
!!$
!!$       write(6,662) ktime%hour,gcmdomainLoc%nlon,gcmdomainLoc%nlat,klon_bound,klat_bound,kk,   &
!!$            rmin,rmax,rmean
!!$661    format(' echam5 river  ',6i6,3f10.6)
!!$662    format(' echam5 runoff ',6i6,3f10.6)
!!$    endif
!!$    !     end for routing


    !     rotation of wind components
    if( rotate ) then
       !     rotate from rotated input geometry to non-rotated geometry, if needed
       if( abs(gcmdomainLoc%polat+90._realkind)>eps ) then
          call turnwirot2reg(pa,pb,pc,pd,  &
               zlon1_u,zlat1_u,zlon2_u,zlat2_u,            &
               klon,klat,&
               gcmdomainLoc%polon,gcmdomainLoc%polat,-2)       
          call turnwi1(uu_grib,vu_grib,uu,vv,pa,pb,pc,pd,   &
               klon,klat,gcmdomainLoc%nlev)
          do jlev=1,gcmdomainLoc%nlev
             do j=1,klat
                do i=1,klon
                   uu_grib(i,j,jlev) = uu(i,j,jlev)
                   vu_grib(i,j,jlev) = vv(i,j,jlev)
                enddo
             enddo
          enddo
          call turnwirot2reg(pa,pb,pc,pd, &
               zlon1_v,zlat1_v,zlon2_v,zlat2_v,           &
               klon,klat,&
               gcmdomainLoc%polon,gcmdomainLoc%polat,-2)      

          call turnwi1(uv_grib,vv_grib,uu,vv,pa,pb,pc,pd,  &
               klon,klat,gcmdomainLoc%nlev)

          do jlev=1,gcmdomainLoc%nlev
             do j=1,klat
                do i=1,klon
                   uv_grib(i,j,jlev) = uu(i,j,jlev)
                   vv_grib(i,j,jlev) = vv(i,j,jlev)
                enddo
             enddo
          enddo
       endif

       !     rotate to rotated output geometry if needed
       if( abs(RCAdom%polat+90._realkind)>eps ) then
          call turnwireg2rot(pa,pb,pc,pd,zlon1_u,zlat1_u,lon_u,lat_u, &
               klon,klat,&
               RCAdom%polon,RCAdom%polat,2)                
          call turnwi1(uu_grib,vu_grib,uu,vv,pa,pb,pc,pd,  &
               klon,klat,gcmdomainLoc%nlev)

          do jlev=1,gcmdomainLoc%nlev
             do j=1,klat
                do i=1,klon
                   uu_grib(i,j,jlev) = uu(i,j,jlev)
                   vu_grib(i,j,jlev) = vv(i,j,jlev)
                enddo
             enddo
          enddo


          call turnwireg2rot(pa,pb,pc,pd, &
               zlon1_v,zlat1_v,lon_v,lat_v,               &
               klon,klat,&
               RCAdom%polon,RCAdom%polat,2)

          call turnwi1(uv_grib,vv_grib,uu,vv,pa,pb,pc,pd,  &
               klon,klat,gcmdomainLoc%nlev)

          do jlev=1,gcmdomainLoc%nlev
             do j=1,klat
                do i=1,klon
                   uv_grib(i,j,jlev) = uu(i,j,jlev)
                   vv_grib(i,j,jlev) = vv(i,j,jlev)
                enddo
             enddo
          enddo

       endif
    endif


    !     4.0 vertical interpolation



    do jlev=1,lev_start-1
       do j=1,klat
          do i=1,klon
             uu_grib(i,j,jlev) = uu_grib(i,j,lev_start) 
             vv_grib(i,j,jlev) = vv_grib(i,j,lev_start) 
             t_grib(i,j,jlev) = t_grib(i,j,lev_start) 
             q_grib(i,j,jlev) = q_grib(i,j,lev_start) 
          enddo
       enddo
    enddo
!
!cau110623
!    write(6,*) ' before etaeta ',mype,RCAdom%fis(1,1),fis_grib(1,1)
!    write(6,*) ' before etaeta ',mype,klev,gcmdomainLoc%nlev
!    if(mype.eq.0) then
!    do i =1,35,5
!       write(6,*) ' pss glob ',i,gcmdomainLoc%afull(i)+gcmdomainLoc%bfull(i)*ps_grib(1,1)
 !   enddo
 !   endif

    call etaeta(klon,klat,klev,                             &
         sip0,.true.,                                       &
         gcmdomainLoc%nlev,                                 &
         gcmdomainLoc%afull,gcmdomainLoc%bfull,             &
         gcmdomainLoc%ahyb,gcmdomainLoc%bhyb,               &
         fis_grib,                                          &
         ps_grib,                                           &
         t_grib,uu_grib,vv_grib,q_grib,                     &
         RCAdom%afull,RCAdom%bfull,RCAdom%fis,ps,t,u,v,q,   &
         psrel,fi500ec)
!    if(mype.eq.0) then
!    do i = 1,24,4
!       write(6,*) ' pss rca ',i,RCAdom%afull(i)+RCAdom%bfull(i)*ps(1,1)
!    enddo
!    endif

    deallocate(fis_grib)
    deallocate(ps_grib, uu, vv, uu_grib, vu_grib, uv_grib,     &
         vv_grib, t_grib,q_grib)                                              

    deallocate(lat_m, lon_m, lat_u, lon_u,lat_v,lon_v,zlat1_m,zlon1_m,  &
         zlat1_u,zlon1_u,zlat1_v,zlon1_v,zlat2_m,zlon2_m,zlat2_u,       &
         zlon2_u,zlat2_v,zlon2_v)

    deallocate(i_m,j_m, i_u, j_u, i_v,j_v)

    deallocate(w_x_m, w_y_m, w_x_u, w_y_u,w_x_v,w_y_v)

    deallocate(pa, pb, pc, pd, fi500ec)
    deallocate(psrel)

    return                                                               

  end subroutine interpol_bd



  subroutine etaeta(klon,klat,klev,                  &
       refpre,lpsmod,                                &
       nlev_grib,                                    &
       af_grib,bf_grib,                              &
       ah_grib,bh_grib,                              &
       fis_grib,ps_grib,                             &
       t_grib,u_grib,v_grib,q_grib,                  &
       af,  bf,  fis,  ps,  temperature,  u,  v,  q, &
       psprel,fi500)

    use asimof
    use interpol
    use confys
    use escom
    use gcm
    implicit none

    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(in):: refpre
    logical,intent(in)::lpsmod
    integer,intent(in)::nlev_grib
    real(kind=realkind),intent(in):: af_grib(nlev_grib),bf_grib(nlev_grib)
    real(kind=realkind),intent(in):: ah_grib(nlev_grib+1),bh_grib(nlev_grib+1)
    real(kind=realkind),intent(in)::fis_grib(klon,klat)
    real(kind=realkind),intent(in)::ps_grib(klon,klat)
    real(kind=realkind),intent(in)::t_grib(klon,klat,nlev_grib)
    real(kind=realkind),intent(in)::u_grib(klon,klat,nlev_grib)
    real(kind=realkind),intent(in)::v_grib(klon,klat,nlev_grib)
    real(kind=realkind),intent(in)::q_grib(klon,klat,nlev_grib)

    real(kind=realkind),intent(in)::af(klev),bf(klev)
    real(kind=realkind),intent(in)::fis(klon,klat)
    real(kind=realkind),intent(inout):: ps(klon,klat),temperature(klon,klat,klev),&
         u(klon,klat,klev),v(klon,klat,klev),q(klon,klat,klev),&
         psprel(klon,klat), fi500(klon,klat)

    real(kind=realkind):: graviti
    real(kind=realkind),allocatable,dimension(:)::zah_grib,zbh_grib,eta,ah,bh
    real(kind=realkind),allocatable,dimension(:,:)::psu,psv,psu_grib,psv_grib,dteta,dfi,inttv_grib
    real(kind=realkind),allocatable,dimension(:,:,:)::lnp_grib,tv_grib, p_grib,pu_grib,&
         pv_grib, pt_grib,rh_grib, ph_grib,fi_grib, ptpbl,&
         rhpbl,upbl,vpbl,peta,pueta,pveta,peh,pehl,rh,fi,eta_grib,zh_grib
    integer,allocatable:: jbltxx(:,:),jblt_grib(:,:),jpsc_grib(:,:)
    real(kind=realkind),allocatable:: ps_x(:,:)
    real(kind=realkind),allocatable:: ps_grib_x(:,:)


    integer jk,jblt,jx,jy,jz,&
         jl1,novblt,jsblt,jbltx,jk_grib
    integer ibpw, ibpwio
    real(kind=realkind) qsat,etapbl,etapsc,epsm1i,w1,w2,ztv

    real(kind=realkind) ztvs_grib,ztvs,zps,zp

    integer imdi

    data ibpw /0/, ibpwio /0/

    !     STATURATION VAPOUR PRESSURE TABLES AND FUNCTIONS


    allocate(eta_grib(klon,klat,nlev_grib),zah_grib(nlev_grib+1),              &
         zbh_grib(nlev_grib+1),                                       &
         eta(klev),     ah(klev+1),     bh(klev+1),                  &
         psu(klon,klat),    psv(klon,klat),                          &
         psu_grib(klon,klat),  psv_grib(klon,klat),dteta(klon,klat), &
         dfi(klon,klat),                                            &
         zh_grib(klon,klat,nlev_grib+1),inttv_grib(klon,klat))

    allocate(lnp_grib(klon,klat,nlev_grib),                           &
         tv_grib(klon,klat,nlev_grib),                                 &
         p_grib(klon,klat,nlev_grib), pu_grib(klon,klat,nlev_grib),    &
         pv_grib(klon,klat,nlev_grib), pt_grib(klon,klat,nlev_grib),   &
         rh_grib(klon,klat,nlev_grib), ph_grib(klon,klat,nlev_grib+1), &
         fi_grib(klon,klat,nlev_grib+1), ptpbl(klon,klat,klev),        &
         rhpbl(klon,klat,klev),    upbl(klon,klat,klev),               &
         vpbl(klon,klat,klev),    peta(klon,klat,klev),                &
         pueta(klon,klat,klev),   pveta(klon,klat,klev),               &
         peh(klon,klat,klev+1),  pehl(klon,klat,klev+1),               &
         rh(klon,klat,klev),      fi(klon,klat,klev+1))
    allocate(ps_x(klon+1,klat+1), ps_grib_x(klon+1,klat+1))
    allocate(jbltxx(klon,klat),jblt_grib(klon,klat),jpsc_grib(klon,klat))

    graviti = 1.0_realkind/gravit
    epsm1i = 1._realkind/epsilo - 1._realkind
    etapbl = 0.80_realkind
    novblt = 3
    etapsc = 0.5_realkind              ! only look for layers to calculate hirlam ps if
    ! eta-ecmwf>=etapsc

    !     get the missing data indicator.

    call asimhm(ibpw, ibpwio, imdi)

    !     preset all output fields to 'missing data'.

    ps  = real(imdi,realkind) !zmdi
    temperature =  real(imdi,realkind) !zmdi
    u =  real(imdi,realkind) !zmdi
    v =  real(imdi,realkind) !zmdi
    q =  real(imdi,realkind) !zmdi

    !     compute large scale full eta values,
    !     find indeces to planetary boundary layer top and
    !     highest possible level of small scale orography

    !     find a:s and b:s of large scale eta half levels
    if(lhadgem)then
    ! Calculate height coordinate for HadGEM
      zah_grib=ah_grib
      zbh_grib=bh_grib
      do  jk=1,nlev_grib+1       
         do jy=1,klat
            do jx=1,klon
               zh_grib(jx,jy,jk)=zah_grib(jk)+zbh_grib(jk)*fis_grib(jx,jy)/gravit
            enddo
         enddo
      enddo
    else  ! all other GCMs
      zah_grib(nlev_grib+1) = 0._realkind
      zbh_grib(nlev_grib+1) = 1._realkind
      do  jk=nlev_grib,2,-1
         zah_grib(jk)=2._realkind*af_grib(jk) - zah_grib(jk+1)
         if(af_grib(jk)<af_grib(jk-1))then
            if(zah_grib(jk)>af_grib(jk-1).or.                  &
                 zah_grib(jk)<af_grib(jk))                      &
                 zah_grib(jk)=0.5_realkind*(af_grib(jk-1)+af_grib(jk))     
         else                                                          
            if(zah_grib(jk)<af_grib(jk-1).or.                   &
                 zah_grib(jk)>af_grib(jk))                      &
                 zah_grib(jk)=0.5_realkind*(af_grib(jk-1)+af_grib(jk))     
         endif
         zbh_grib(jk)=2._realkind*bf_grib(jk) - zbh_grib(jk+1)              
         if(bf_grib(jk)<bf_grib(jk-1))then                      
            if(zbh_grib(jk)>bf_grib(jk-1).or.                   &
                 zbh_grib(jk)<bf_grib(jk))                      &
                 zbh_grib(jk)=0.5_realkind*(bf_grib(jk-1)+bf_grib(jk))     
         else                                                          
            if(zbh_grib(jk)<bf_grib(jk-1).or.                   &
                 zbh_grib(jk)>bf_grib(jk))                      &
                 zbh_grib(jk)=0.5_realkind*(bf_grib(jk-1)+bf_grib(jk))
         endif
      enddo
      zah_grib(1)        = 0._realkind
      zbh_grib(1)        = 0._realkind
    endif           ! lhadgem

    !     FIND ETA FULL VALUES OF SMALL SCALE FIELD AND
    !     INDEX TO PLANETARY BOUNDARY LAYER TOP

    do jk=1,klev
       eta(jk) =  af(jk)/refpre + bf(jk)
    enddo
    do jk=klev,1,-1
       if( eta(jk)>etapbl ) jblt = jk
    enddo

    !     find a:s and b:s of small scale eta half levels

    ah(klev+1) = 0._realkind
    bh(klev+1) = 1._realkind
    do jk=klev,1,-1
       ah(jk) = max(2._realkind*af(jk) - ah(jk+1), 0._realkind)
       bh(jk) = max(2._realkind*bf(jk) - bh(jk+1), 0._realkind)
    enddo

    !     COMPUTE PRELIMINARY SURFACE PRESSURE OF SMALL
    !     SCALE T AND Q- POINTS
    !                 DO THE CALCULATIONS FOR klon*klat POINTS AT A TIME
    !                 IN ORDER TO SAVE WORKSPACE
    !     tf 060202 not any longer, do the computation over whole klon*klat
    !     = klon*klat at once in order to be able to do a correct
    !     surface pressure interpolation


    !     FIND LOGARITHS OF PRESSURE AND VIRTUAL
    !     TEMPERATURES NEEDED OR TEMPERATURE INTERPOLATION
    !     AND EXTRAPOLATION

    do  jk=1,nlev_grib
       do jy=1,klat
          do jx=1,klon
             p_grib(jx,jy,jk)  = af_grib(jk) +  bf_grib(jk)*ps_grib(jx,jy)
             lnp_grib(jx,jy,jk)=log(p_grib(jx,jy,jk))
             tv_grib(jx,jy,jk)=(1.0_realkind+epsm1i*q_grib(jx,jy,jk))*t_grib(jx,jy,jk)
          enddo
       enddo
    enddo

    !FIND PRESSURES NEEDED FOR GEOPOTENTIAL
    !                     CALCULATION
    !
    if(lhadgem)then
! Calculate pressure at half levels, ph_grib, by the hypsometric equation:
! p(k-1/2) = ps * exp( -gravit/rair INT_0^z(k-1/2) ( 1/Tv(z) dz ) )
! where the integral INT_0^z(k-1/2) ( 1/Tv(z) dz ) = SUM ( 0.5*( 1/Tv(k-1/2) + 1/Tv(k+1/2) ) * ( z(k-1/2)-z(k+1/2) ) )
! which to a good approximation, considering the values of Tv, can be written as
!                                                    SUM ( 1/Tv(k) * ( z(k-1/2)-z(k+1/2) ) )
      inttv_grib=0._realkind
      do  jk=nlev_grib+1,2,-1
         do jy=1,klat
            do jx=1,klon
               ph_grib(jx,jy,jk)=ps_grib(jx,jy)*exp(-gravit/rair*inttv_grib(jx,jy))
               inttv_grib(jx,jy)=inttv_grib(jx,jy) +  &
                                 1._realkind/tv_grib(jx,jy,jk-1) * ( zh_grib(jx,jy,jk-1) - zh_grib(jx,jy,jk) )
            enddo
         enddo
      enddo
! Pressure at top half level
      do jy=1,klat
         do jx=1,klon
            ph_grib(jx,jy,1)=ps_grib(jx,jy)*exp(-gravit/rair*inttv_grib(jx,jy))
         enddo
      enddo
! Pressure at full levels is calculated as arithmetic mean between half levels
      do  jk=1,nlev_grib
         do jy=1,klat
            do jx=1,klon
               p_grib(jx,jy,jk)=0.5*( ph_grib(jx,jy,jk) + ph_grib(jx,jy,jk+1) )
               lnp_grib(jx,jy,jk)=log( ph_grib(jx,jy,jk) )
            enddo
         enddo
      enddo

    else    ! All other GCMs

      do  jk=2,nlev_grib+1
         do jy=1,klat
            do jx=1,klon
               ph_grib(jx,jy,jk)=zah_grib(jk)+zbh_grib(jk)*ps_grib(jx,jy)
            enddo
         enddo
      enddo
    endif

    if(lhadgem)then
      do jk=1,nlev_grib
         do jy=1,klat
            do jx=1,klon
              eta_grib(jx,jy,jk) = p_grib(jx,jy,jk)/ps_grib(jx,jy)
            enddo
         enddo
      enddo
    else   ! All other GCMs
      do jk=1,nlev_grib
         do jy=1,klat
            do jx=1,klon
! RCA original version
              eta_grib(jx,jy,jk) = af_grib(jk)/refpre + bf_grib(jk)
! RCA new alternative version
!!!!              eta_grib(jx,jy,jk) = p_grib(jx,jy,jk)/ps_grib(jx,jy)
            enddo
         enddo
      enddo
    endif
    do jk=nlev_grib,1,-1
       do jy=1,klat
          do jx=1,klon
             if( eta_grib(jx,jy,jk)>=etapbl ) jblt_grib(jx,jy) = jk
             if( eta_grib(jx,jy,jk)>=etapsc ) jpsc_grib(jx,jy) = jk
          enddo
       enddo
    enddo
    !                      FIND GEOPOTENTIALS NEEDED FOR SURFACE PRESSURE
    !                      CALCULATIONS

    do jy=1,klat
       do jx=1,klon
          fi_grib(jx,jy,nlev_grib+1) = fis_grib(jx,jy)
       enddo
    enddo

    do  jk=nlev_grib,2,-1
       do jy=1,klat
          do jx=1,klon
             fi_grib(jx,jy,jk) = fi_grib(jx,jy,jk+1) +&
                  rair*tv_grib(jx,jy,jk)*log( ph_grib(jx,jy,jk+1)/ph_grib(jx,jy,jk))
          enddo
       enddo
    enddo



    do jy=1,klat
       do jx=1,klon
          if( fis(jx,jy)<=fis_grib(jx,jy) ) then
             ztvs_grib = tv_grib(jx,jy,nlev_grib) +                        &
                  ( log(ps_grib(jx,jy)) - lnp_grib(jx,jy,nlev_grib) )/       &
                  (lnp_grib(jx,jy,nlev_grib-1) -lnp_grib(jx,jy,nlev_grib))*  &
                  ( tv_grib(jx,jy,nlev_grib-1) - tv_grib(jx,jy,nlev_grib) )   
             ztvs   =  ztvs_grib +                                         &
                  ( fis_grib(jx,jy) - fis(jx,jy) )*graviti*0.0065_realkind             
             ps(jx,jy) = ps_grib(jx,jy)*                                     &
                  exp(  ( fis_grib(jx,jy) - fis(jx,jy) )/                    &
                  (  rair*0.5_realkind*(ztvs_grib+ztvs)  )      )
          endif
       enddo
    enddo


    do jy=1,klat
       do jx=1,klon
          if( fis(jx,jy)<=fi_grib(jx,jy,nlev_grib) .and.                   &
               fis(jx,jy)>fi_grib(jx,jy,nlev_grib+1)     ) then             
             zps = ph_grib(jx,jy,nlev_grib+1)*                               &
                  exp( ( -fis(jx,jy) + fi_grib(jx,jy,nlev_grib+1) )/          &
                  (rair*tv_grib(jx,jy,nlev_grib))  )                          
             zp = 0.5_realkind*( zps + ph_grib(jx,jy,nlev_grib+1) )                    
             ztv = tv_grib(jx,jy,nlev_grib) +                                &
                  ( log(zp) - lnp_grib(jx,jy,nlev_grib) )/                   &
                  (lnp_grib(jx,jy,nlev_grib-1)-lnp_grib(jx,jy,nlev_grib) )*   &
                  ( tv_grib(jx,jy,nlev_grib-1) - tv_grib(jx,jy,nlev_grib) )    
             ps(jx,jy) =                                                     &
                  ph_grib(jx,jy,nlev_grib+1)*                                &
                  exp( ( -fis(jx,jy) + fi_grib(jx,jy,nlev_grib+1) )/          &
                  (rair*ztv)  )
          endif
       enddo
    enddo

    do jy=1,klat
       do jx=1,klon
          do jk=nlev_grib-1,max(jpsc_grib(jx,jy)-1,1),-1
             if( fis(jx,jy)<=fi_grib(jx,jy,jk) .and.               &
                  fis(jx,jy)>fi_grib(jx,jy,jk+1)     ) then         
                zps = ph_grib(jx,jy,jk+1)*                           &
                     exp( ( -fis(jx,jy) + fi_grib(jx,jy,jk+1) )/      &
                     (RAIR*tv_grib(jx,jy,jk))  )                      
                zp = 0.5_realkind*( zps + ph_grib(jx,jy,jk+1) )                
                ztv = tv_grib(jx,jy,jk) +                            &
                     ( log(zp) - lnp_grib(jx,jy,jk) )/               &
                     ( lnp_grib(jx,jy,jk+1) - lnp_grib(jx,jy,jk) )* &
                     ( tv_grib(jx,jy,jk+1) - tv_grib(jx,jy,jk) )     
                ps(jx,jy) =  ph_grib(jx,jy,jk+1)*                     &
                     exp( ( -fis(jx,jy) + fi_grib(jx,jy,jk+1) )/      &
                     (rair*ztv)  )
             endif
          enddo
       enddo
    enddo
    !     FIND 500 HPA GEOPOTENTIALS

    if( lpsmod ) then
       do jk=nlev_grib,2,-1
          do jy=1,klat
             do jx=1,klon
                if( ph_grib(jx,jy,jk)<50000._realkind .and.                &
                     ph_grib(jx,jy,jk+1)>=50000._realkind ) then              
                   fi500(jx,jy) =                                       &
                        (  log(ph_grib(jx,jy,jk+1)/50000._realkind)*           &
                        fi_grib(jx,jy,jk) +                           &
                        log(50000._realkind/ph_grib(jx,jy,jk))*                &
                        fi_grib(jx,jy,jk+1) )/                        &
                        log(ph_grib(jx,jy,jk+1)/ph_grib(jx,jy,jk))
                endif
             enddo
          enddo
       enddo
    endif

    psprel = ps
    !            DO THE REST OF THE CALCULATIONS FOR KLON*KLAT POINTS AT A TIME
    !            IN ORDER TO SAVE WORKSPACE


    !     SURFACE PRESSURE VALUES ARE INTERPOLATED TO THE STAGGERED
    !     U- AND V- GRIDPOINTS (IN CASE OF C-GRID)


!hadgemcheck this for hadgem. All CMIP5 GCM data is on an Arakawa A grid, so this may not do anything.
!hadgem It should work directly anyhow as ps_grib is directly available in all files same as for other
!hadgem CMIP5 GCMs
    call  add_dim(klon,klat,ps_x,ps_grib_x,ps,ps_grib)
    do jy=1,klat
       do jx=1,klon
          if (abs(ps_grib_x(jx+1,jy)- real(imdi,realkind))<1.e-14_realkind)then !zmdi) then
             psu(jx,jy) = ps_x(jx,jy)
             psu_grib(jx,jy) = ps_grib_x(jx,jy)
          else
             psu(jx,jy) = 0.5_realkind*(ps_x(jx,jy)+ ps_x(jx+1,jy) )
             psu_grib(jx,jy) = 0.5_realkind*( ps_grib_x(jx,jy) + ps_grib_x(jx+1,jy) )
          endif

          if (abs(ps_grib_x(jx,jy+1) - real(imdi,realkind))<1.e-14_realkind)then !zmdi) then
             psv(jx,jy) =  ps_x(jx,jy)
             psv_grib(jx,jy) = ps_grib_x(jx,jy)
          else
             psv(jx,jy) = 0.5_realkind*(ps_x(jx,jy) + ps_x(jx,jy+1) )
             psv_grib(jx,jy) = 0.5_realkind*( ps_grib_x(jx,jy) + ps_grib_x(jx,jy+1) )

          endif
       enddo
    enddo

    !              COMPUTE INPUT PRESSURES, POTENTIAL TEMPERATURES AND
    !                RELATIVE HUMIDITIES

!hadgem I am not sure if the code above tries to stagger U,V onto the RCA arakawa C grid
!hadgem If yes then we will need to directly calulate the staggering for all pu/pv_grib
!hadgem vertical levels. I think trhsi can work by pu_grib=p_grib*(psu_grib/ps_grib)
!hadgem and same for pv_grib
    if(lhadgem)then

      do  jk=1,nlev_grib
         do jy=1,klat
            do jx=1,klon
               pu_grib(jx,jy,jk) = p_grib(jx,jy,jk) * (psu_grib(jx,jy)/ps_grib(jx,jy))
               pv_grib(jx,jy,jk) = p_grib(jx,jy,jk) * (psv_grib(jx,jy)/ps_grib(jx,jy))
            enddo
         enddo
      enddo

    else    ! All other GCMs

      do  jk=1,nlev_grib
         do jy=1,klat
            do jx=1,klon
               pu_grib(jx,jy,jk) = af_grib(jk) +  bf_grib(jk)*psu_grib(jx,jy)
               pv_grib(jx,jy,jk) = af_grib(jk) +  bf_grib(jk)*psv_grib(jx,jy)
            enddo
         enddo
      enddo

    endif

    do  jk=1,nlev_grib
       do  jy=1,klat
          do jx=1,klon
             pt_grib(jx,jy,jk) = t_grib(jx,jy,jk)*(refpre/p_grib(jx,jy,jk))**(rair/cpair)
          enddo
       enddo
    enddo

    do  jk=1,nlev_grib
       do jy=1,klat
          do jx=1,klon
             qsat = epsilo*esat(t_grib(jx,jy,jk))/p_grib(jx,jy,jk)
             rh_grib(jx,jy,jk) = q_grib(jx,jy,jk)/qsat
          enddo
       enddo
    enddo


    !     interpolate in boundary layer with eta
    !     as vertical coordinate
    jsblt = max(1,jblt-novblt)
    do 590 jk=jsblt,klev
       do jy=1,klat
          do jx=1,klon
             if( eta(jk)<=eta_grib(jx,jy,1) ) then
                jl1 = 1
                w1 = 1._realkind
                w2 = 0._realkind
             else if( eta(jk)>=eta_grib(jx,jy,nlev_grib) ) then
                jl1 = nlev_grib - 1
                w1 = 0._realkind
                w2 = 1._realkind
             else
                do jk_grib=1,nlev_grib-1
                   if( eta(jk)>=eta_grib(jx,jy,jk_grib) .and.&
                        eta(jk)<eta_grib(jx,jy,jk_grib+1) ) then
                      jl1 = jk_grib
                      w2 = log(eta(jk)/eta_grib(jx,jy,jk_grib) )/ &
                           log( eta_grib(jx,jy,jk_grib+1)/eta_grib(jx,jy,jk_grib) )
                      w1 = 1._realkind - w2
                   endif
                enddo
             endif

             ptpbl(jx,jy,jk) = w1*pt_grib(jx,jy,jl1) +   w2*pt_grib(jx,jy,jl1+1)
             rhpbl(jx,jy,jk) = w1*rh_grib(jx,jy,jl1) +   w2*rh_grib(jx,jy,jl1+1)
             upbl(jx,jy,jk)  = w1*u_grib(jx,jy,jl1) +   w2*u_grib(jx,jy,jl1+1)
             vpbl(jx,jy,jk)  = w1*v_grib(jx,jy,jl1) +   w2*v_grib(jx,jy,jl1+1)
          enddo
       enddo
590 enddo


    !     compute pressures for small-scale fields and
    !     interpolate in free atmosphere with pressure
    !     as vertical coordinate

    do jk=1,klev
       do jy=1,klat
          do jx=1,klon
             peta(jx,jy,jk)   = af(jk) + bf(jk)*ps(jx,jy)
             pueta(jx,jy,jk)  = af(jk) + bf(jk)*psu(jx,jy)
             pveta(jx,jy,jk)  = af(jk) + bf(jk)*psv(jx,jy)
             peh(jx,jy,jk+1)  = ah(jk+1) + bh(jk+1)*ps(jx,jy)
             pehl(jx,jy,jk+1) =  log( peh(jx,jy,jk+1) )
          enddo
       enddo
    enddo


    do  jk=1,jblt
       call vint(nlev_grib,klon*klat,t_grib,p_grib, temperature(:,:,jk),peta(:,:,jk))
       call vint(nlev_grib,klon*klat,rh_grib,p_grib,rh(:,:,jk)         ,peta(:,:,jk))
       call vint(nlev_grib,klon*klat,u_grib,pu_grib,u(:,:,jk)          ,pueta(:,:,jk))
       call vint(nlev_grib,klon*klat,v_grib,pv_grib,v(:,:,jk)          ,pveta(:,:,jk))
    enddo


    !     COMPUTE BOUNDARY LAYER TEMPERATURE CORRECTION FACTOR
    !     AND MERGE BOUNDARY LAYER AND FREE ATMOSPHERE
    !     PROFILES

    do  jy=1,klat
       do jx=1,klon
          jbltxx(jx,jy) = jblt
       enddo
    enddo

    do jk=jblt-1,jsblt,-1
       do  jy=1,klat
          do jx=1,klon
             jz=jblt_grib(jx,jy)
             if( peta(jx,jy,jk)>p_grib(jx,jy,jz) )jbltxx(jx,jy) = jk
          enddo
       enddo
    enddo

    do jy=1,klat
       do jx=1,klon
          jbltx = jbltxx(jx,jy)
          dteta(jx,jy) = temperature(jx,jy,jbltx)* &
               (refpre/peta(jx,jy,jbltx))**(rair/cpair) - ptpbl(jx,jy,jbltx)
       enddo
    enddo

    do  jk=jsblt,klev
       do  jy=1,klat
          do jx=1,klon
             if( jk>=jbltxx(jx,jy)  )then
                temperature(jx,jy,jk) = (ptpbl(jx,jy,jk)+ dteta(jx,jy))/&
                     (refpre/peta(jx,jy,jk))**(rair/cpair)

                rh(jx,jy,jk) = rhpbl(jx,jy,jk)
                u(jx,jy,jk) = upbl(jx,jy,jk)
                v(jx,jy,jk) = vpbl(jx,jy,jk)
             endif
          enddo
       enddo
    enddo


    !     compute 500 hpa geopotential from interpolated profile,
    !     compare with large scale 500 hpa geopotential and
    !     correct surface pressure so that large scale  and
    !     small scale geopotentials agree

    if( lpsmod ) then
       do  jy=1,klat
          do jx=1,klon
             fi(jx,jy,klev+1) = fis(jx,jy)
          enddo
       enddo

       do  jk=klev,2,-1
          do  jy=1,klat
             do jx=1,klon
                qsat = epsilo*esat(temperature(jx,jy,jk))/peta(jx,jy,jk)
                q(jx,jy,jk) = rh(jx,jy,jk)*qsat
                fi(jx,jy,jk) = fi(jx,jy,jk+1) + rair*(1._realkind+epsm1i*q(jx,jy,jk))*    &
                     temperature(jx,jy,jk)*(pehl(jx,jy,jk+1)-pehl(jx,jy,jk))
             enddo
          enddo
       enddo

       do jk=klev,2,-1
          do jy=1,klat
             do jx=1,klon
                if(peh(jx,jy,jk)<50000._realkind.and.peh(jx,jy,jk+1)>=50000._realkind)then
                   dfi(jx,jy) = fi500(jx,jy) -                  &
                        (log(peh(jx,jy,jk+1)/50000._realkind)*         &
                        fi(jx,jy,jk) +                        &
                        log(50000._realkind/peh(jx,jy,jk))*            &
                        fi(jx,jy,jk+1) )/                     &
                        log(peh(jx,jy,jk+1)/peh(jx,jy,jk))
                endif
             enddo
          enddo
       enddo

       do jy=1,klat
          do jx=1,klon
             ztv=(1._realkind+epsm1i*q(jx,jy,klev))*temperature(jx,jy,klev)
             ps(jx,jy) = ps(jx,jy)*exp( dfi(jx,jy)/(rair*ztv) )
          enddo
       enddo
    endif

    !     compute corrected pressures and
    !     convert relative humidities to specific humidities
    do  jk=1,klev
       do   jy=1,klat
          do jx=1,klon
             peta(jx,jy,jk) =  af(jk) + bf(jk)*ps(jx,jy)
             qsat = epsilo*esat(temperature(jx,jy,jk))/peta(jx,jy,jk)
             q(jx,jy,jk) = rh(jx,jy,jk)*qsat
          enddo
    enddo
 enddo

    deallocate(eta_grib,   zah_grib,   zbh_grib,    eta,     ah,     bh)
    deallocate(psu,    psv,psu_grib,  psv_grib,    dteta,  dfi, &
         zh_grib,inttv_grib)
    deallocate(lnp_grib,tv_grib,p_grib,pu_grib,  pv_grib,    pt_grib, &
         rh_grib,    ph_grib, fi_grib, ptpbl, rhpbl,    upbl, &
         vpbl,    peta, pueta,   pveta,peh,  pehl, rh,      fi)

    deallocate(jbltxx)
    deallocate(ps_x,ps_grib_x)
    return

  end subroutine etaeta


  subroutine add_dim(klon,klat,ps_x,ps_grib_x,ps,ps_grib)

    implicit none
    integer,intent(in):: klon,klat
    real(kind=realkind):: ps_x(klon+1,klat+1),ps_grib_x(klon+1,klat+1)
    real(kind=realkind),intent(in)::ps(klon,klat),ps_grib(klon,klat)

    !     Copy 'inner' into _x
    ps_x(1:klon,1:klat) = ps
    ps_grib_x(1:klon,1:klat) = ps_grib

    !     the rightmost column of _x is extrapolated using 0:th order extrapolation
    if(atright) then
       ps_x(klon+1,1:klat) =  ps_x(klon,1:klat)
       ps_grib_x(klon+1,1:klat) =  ps_grib_x(klon,1:klat)
    endif
    if(attop) then
       ps_x(1:klon,klat+1) =  ps_x(1:klon,klat)
       ps_grib_x(1:klon,klat+1) =  ps_grib_x(1:klon,klat)
    endif
    if(attop.and.atright)then
       ps_x(klon+1,klat+1) =  ps_x(klon,klat)
       ps_grib_x(klon+1,klat+1) =  ps_grib_x(klon,klat)
    endif

    call swap_ps(ps_x,klon,klat)
    call swap_ps(ps_grib_x,klon,klat)
    return
  end subroutine add_dim


  subroutine  estimateCloudCondensate(klon,klat,klev,RCAdom,bc)
    use confys
    use escom 
    use modddr
    implicit none
    type(domain),intent(in)::RCAdom
    integer,intent(in)::klon,klat,klev
    type(atm),intent(inout)::bc
    real(kind=realkind)::zrhuc(klev),zrelf0,zrelf1,zrelf2
    real(kind=realkind)::zdrelf,zpdps,zcloud,zpi,zearth,zrhu,zmepsi
    real(kind=realkind)::zdist,zrelfh,zqsat
    real(kind=realkind)::zcval1,zcval2
    real(kind=realkind)::zesatd 
    integer:: jx,jy,jk
    zmepsi=1.0_realkind-epsilo
    zpi=pi
    zearth=rearth
    zdist=sqrt( zearth*zpi*RCAdom%dlat/180.0_realkind )
    zrelf0=0.0200_realkind
    zrelf1=0.3000_realkind
    zrelf2=0.0030_realkind
    zdrelf=0.0300_realkind
    zrelfh=zrelf1*(1._realkind -exp( -zrelf2*zdist ))

    do jk=1,klev
       zpdps=min(RCAdom%hybk(jk),1.0_realkind)
       zrhuc(jk)=zrelfh*(zrelf0+zdrelf*(1._realkind-zpdps*zpdps*zpdps))/( zrelf0 +zdrelf )
    enddo

    !     Correct cloud condensate field 'PSZ' on the lateral
    !     boundaries to be consistent with a stratiform condensation
    !     scheme.
    do jk=1,klev
       do jy=1,klat
          do jx=1,klon
             zesatd=esatfun(bc%t(jx,jy,jk)) 
             zqsat=epsilo*zesatd*1._realkind/(RCAdom%afull(jk)+RCAdom%bfull(jk)*bc%ps(jx,jy) - &
                  zmepsi*zesatd )
             if (bc%q(jx,jy,jk)>zqsat)  bc%q(jx,jy,jk) = zqsat
             if (bc%q(jx,jy,jk)<0._realkind    ) bc%q(jx,jy,jk) = 1.e-7_realkind
             zrhu=max( bc%q(jx,jy,jk)/zqsat, zrhuc(jk) )
             zrhu=min( zrhu, 1.0_realkind)
             zcloud=1.0_realkind -sqrt( (1._realkind-zrhu)/(1._realkind-zrhuc(jk)) )

             !     Estimate cloud water on the lateral boundaries.
             !     The estimate should be approximately consistent
             !     with the stratiform condensation scheme in order
             !     to work in an optimal way.
             bc%cw(jx,jy,jk)=0._realkind ! it was not on the file, remember!
             zcval1=(bc%q(jx,jy,jk) +bc%cw(jx,jy,jk))*(1._realkind-zrhuc(jk))
             zcval2=(bc%q(jx,jy,jk) +bc%cw(jx,jy,jk))*(1._realkind+zrhuc(jk))

             if( zcval1>zqsat ) then
                bc%cw(jx,jy,jk)=zqsat*zrhuc(jk)/(1._realkind-zrhuc(jk))
             endif

             if( (zcval1<=zqsat).and.(zcval2>zqsat) ) then
                bc%cw(jx,jy,jk)=zqsat/(1._realkind+zrhuc(jk)-2._realkind*zrhuc(jk)*zcloud)-bc%q(jx,jy,jk)
             endif
          enddo
       enddo
    enddo
  end subroutine estimateCloudCondensate


end module lateralBc
