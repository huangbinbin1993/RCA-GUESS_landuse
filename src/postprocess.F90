module postprocess
  use timetype
  use referenceParameters
  use calendar
  use gcm
  use RCAdomainMod
  use decomp
  use derived_types
  use mod_grib,only:realgribkind
  implicit none
  private

  !    declaration of post-processing parameters
  integer,parameter:: jpfix=2
  integer,parameter:: jsfix=2
  integer,parameter:: maxSize=200


  integer,save::npplon=-666
  integer,save::npplat=-666
  integer,save::iminpp=-666
  integer,save::jminpp=-666

  logical,save::lphys=.true.
  logical,save::lomega=.true.
  logical,save::month_file=.true.

  integer,public,save::nppstr=0 !#output files

  type,public::file
     integer::unitnr=-1 !unit nr.
     character(len=jpfix):: prefix='fc'
     character(len=jsfix):: sufix='  '
     type(time)::write2fileTime
     type(deltatime)::outdt

     integer::nlevml=0 !# ml-fields
     integer::nwmoml=0 !number of components of multi-level variables
     integer::ltypml=-666 !grib-code  type for multi-level variable
     integer::alevml(maxSize)=-666 !grib-code level for multi-level variable
     integer::iwmoml(maxSize)=-666 !grib-code parameter for multi-level variable
     integer::npplon=-666
     integer::npplat=-666
     integer::iminpp=-666
     integer::jminpp=-666
     integer::nsl=0    !# sl-fields
     integer::ltypsl(maxSize)=-666
     integer::iwmosl(maxSize)=-666
     integer::alevsl(maxSize)=-666

     integer::meantime_surf(maxSize)=0
     integer::meantime_3006(maxSize)=0

     integer::scratchtime_acum(maxSize)=0
     integer::scratchtime_4006(maxSize)=0
     integer::scratchtime_extrem(maxSize)=0
     integer::scratchtime_hilo(maxSize)=0

  end type file

  type(file),allocatable,public,save,dimension(:)::files


  public postproc,intvert2,intvert_moist,pslcom,find2points,read_nampos, cloud_cover
  public read_namvar

contains
  subroutine extend(klon,klat,pa)
    !     purpose.
    !     extend variable pa(klon,klat) from inner points
    !     to the boundary points where values are not given.
    use boundaryRelaxation, only:npbpts,nbdpts
    implicit none

    integer,intent(in):: klon  , klat  
    integer::iw , ie , is , in
    real(kind=realkind) pa(klon,klat)
    integer jx,jy,istart,jstart,istop,jstop
    !     start- and endpoints in the local array

    iw = npbpts + 1
    ie = klon_global  - npbpts 
    is = npbpts + 1
    in = klat_global  - npbpts 

    istart = 1
    istop = klon
    jstart = 1
    jstop = klat
    if(idatastart < iw) then
       istart = 1 + iw - idatastart
    endif
    if(idatastart+klon-1 > ie) then
       istop = 1 + ie - idatastart
    endif
    if(jdatastart < is) then
       jstart = 1 + is - jdatastart
    endif
    if(jdatastart+klat-1 > in) then
       jstop = 1 + in - jdatastart
    endif
    
    do jx=1,istart-1
       pa(jx,:) = pa(istart,:)
    enddo
    do jx=istop+1,klon
       pa(jx,:) = pa(istop,:)
    enddo

    do jy=1,jstart-1
       pa(:,jy) = pa(:,jstart)
    enddo
    do jy=jstop+1,klat
       pa(:,jy) = pa(:,jstop)
    enddo
    return
  end subroutine extend


  subroutine init_file(fileI,suff,nfiles,npplon,npplat,iminpp,jminpp)
    use comhkp
#ifdef USING_NETCDF
    use netcdfrestart
#else
    use restart
#endif
    implicit none
    type(file)::fileI
    character(len=jsfix),intent(in):: suff
    integer,intent(in)::nfiles,npplon,npplat,iminpp,jminpp
    integer::it

    integer::nlevmlp,nwmomlp,nslp
    integer::lunppfp
    integer::ltypmlp
    character(len=jpfix)::prefixp
    character(len=jsfix)::sufixp
    integer::timeppp
    integer::alevmlp(maxSize), alevslp(maxSize)
    integer:: iwmomlp(maxSize), ltypslp(maxSize), iwmoslp(maxSize)
    integer::i
    integer::outputInterval

    namelist /namppp/ nlevmlp, nwmomlp, nslp, &
         lunppfp , timeppp, ltypmlp, alevmlp, iwmomlp , ltypslp, iwmoslp, alevslp, &
         prefixp, sufixp   

    if(nfiles>10)call stop_program( 'We currently do not allow more than 10 different output GRIB files')

    nlevmlp=-666
    nwmomlp=-666
    nslp=0
    lunppfp=810
    timeppp=-666
    ltypmlp=-666
    alevmlp=-666
    iwmomlp=-666
    ltypslp=-666
    iwmoslp=-666
    alevslp=-666
    prefixp = 'fc'
    sufixp = '  '

    !    if(suff/='  ')then
    open(57,file='namelists_namppp.dat',status='old')
    read(57,nml=namppp)
    it=1
    do while(suff /= sufixp .and. it<=nfiles) !find the correct namelist
       nlevmlp=-666
       nwmomlp=-666
       nslp=0
       lunppfp=810
       timeppp=-666
       ltypmlp=-666
       alevmlp=-666
       iwmomlp=-666
       ltypslp=-666
       iwmoslp=-666
       alevslp=-666
       prefixp = 'fc'
       sufixp = '  '
       read(57,nml=namppp)
       it = it + 1
    enddo
    if(it>nfiles)then
       print *,__FILE__,__LINE__
       print *,suff, 'not found'
       call stop_program( 'namelist not found: suffix does not match to present namelists')
    endif
    close(57)


    fileI%nwmoml = nwmomlp
    fileI%unitnr = lunppfp
    if(fileI%unitnr<801 .or. fileI%unitnr>810)then
       print *,__FILE__,__LINE__
       print *,fileI%unitnr,fileI%prefix,fileI%sufix
       call stop_program( 'unitNR must be [801-810]')
    endif
    fileI%ltypml = ltypmlp
    fileI%sufix  = sufixp
    fileI%prefix = prefixp
    fileI%npplon = npplon
    fileI%npplat = npplat
    fileI%iminpp = iminpp
    fileI%jminpp = jminpp


    select case(fileI%sufix)
    case('ss')
       outputInterval = (ndtime%hour*60 +ndtime%min)*60 + ndtime%sec !every timestep
    case('dd')
       outputInterval = 86400 !24h
    case('pp')
       outputInterval = 21600 !6h
    case('qq')
       outputInterval = 10800 !3h
    case('hh')
       outputInterval = 3600 !1h
    case('  ')
       outputInterval = 21600 !6h   
    case('gb')
       outputInterval =86400*30 ! makedeltatime(0,1,0,0,0,0) !monthly
    case default
       outputInterval = timeppp !user specified
    end select

    fileI%outdt = makedeltatime(0,0,0,0,0,outputInterval)
    if(fileI%sufix=='gb')then
       fileI%outdt = makedeltatime(0,1,0,0,0,0) !monthly
    endif

    fileI%write2fileTime = startTime 


    if(ltypmlp /= -666)then                      !output levels
       fileI%nlevml=0
       do while(alevmlp(fileI%nlevml+1) /= -666) !.and.iwmomlp(fileI%nlevml+1)/= -666) cj110608
         fileI%nlevml = fileI%nlevml+1
          if(fileI%nlevml>maxSize)then
             print *,__FILE__,__LINE__
             call stop_program( '# ML-fields exeeds maxSize')
          endif
          fileI%alevml(fileI%nlevml) = alevmlp(fileI%nlevml)
       enddo
    endif

    if(ltypmlp /= -666)then                      !output variables
       fileI%nwmoml=0
       do while(iwmomlp(fileI%nwmoml+1)/= -666)
         fileI%nwmoml = fileI%nwmoml+1
         fileI%iwmoml(fileI%nwmoml) = iwmomlp(fileI%nwmoml)
       enddo
    endif


    fileI%nsl=0
    do while(ltypslp(fileI%nsl+1) /= -666 .and. iwmoslp(fileI%nsl+1) /= -666 .and. &
         alevslp(fileI%nsl+1) /= -666)
       fileI%nsl=fileI%nsl+1
       fileI%ltypsl(fileI%nsl) = ltypslp(fileI%nsl)
       fileI%iwmosl(fileI%nsl) = iwmoslp(fileI%nsl)
       fileI%alevsl(fileI%nsl) = alevslp(fileI%nsl)
       if(ltypslp(fileI%nsl)==109.and.(alevslp(fileI%nsl)>klev_global.or.alevslp(fileI%nsl)<1))then
          print *,__FILE__,__LINE__
          print *,'nampp for file ',fileI%prefix,'**',fileI%sufix,alevslp(fileI%nsl),'>',klev_global
          call stop_program( 'trying to require output on a level outside domain')
       endif
    enddo


    !outputInterval is expressed in seconds!
    do i=1,fileI%nsl
       if(fileI%iwmosl(i) == 242) then
          fileI%meantime_surf(i)=outputInterval
       endif
       if(fileI%iwmosl(i) == 245) then
          fileI%scratchtime_acum(i)=outputInterval
       endif
       if(fileI%alevsl(i) == 3006) then
          fileI%meantime_3006(i)=outputInterval
       endif
       if(fileI%alevsl(i) == 4006) then
          fileI%scratchtime_4006(i)=outputInterval
       endif
       if(fileI%iwmosl(i) == 241) then
          fileI%scratchtime_extrem(i)=outputInterval
       endif
       if(fileI%iwmosl(i) == 15 .or.fileI%iwmosl(i) == 16 .or. fileI%iwmosl(i) == 32  .or. fileI%iwmosl(i) == 64) then
          fileI%scratchtime_hilo(i)=outputInterval
       endif
    enddo
  end subroutine init_file

  subroutine read_namvar(ksvar,iacdg,iacdg2)
    implicit none
    !ksvar  = no of extra variables
    !iacdg  = no of 3d fields to be accumulated
    !iacdg2 = no of 2d fields to be accumulated
    integer,intent(out)::ksvar,iacdg,iacdg2
    namelist /namvar/ksvar,iacdg,iacdg2
    ksvar  = 1
    iacdg  = 1   
    iacdg2 = 13 
    open(57,file='namelists.dat',status='old')
    read(57,nml=namvar)
    close(57)
    iacdg  = max(iacdg,1)
    iacdg2 = max(iacdg2,2)

    if(mype==0)then
       write(6,nml=namvar)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namvar)
       close(1)
    endif
  end subroutine read_namvar

  subroutine read_nampos()
    use boundaryRelaxation
    implicit none
    
    integer::j
    logical::linner=.true.
    character(len=jsfix)::suff(10)
    namelist /nampos/lphys,lomega,npplon,npplat,iminpp,jminpp,&
         linner
    namelist /nampp/nppstr,suff,month_file
    
    suff(1)='pp'
    suff(2)='dd'
    suff(3)='qq'
    suff(4)='ss'
    suff(5)='gb'
    suff(6)='hh'
    suff(7:10)='  '
    npplon=-666
    npplat=-666
    iminpp=-666
    jminpp=-666


    open(57,file='namelists_namppp.dat',status='old')
    read(57,nml=nampp)
    close(57)
    open(57,file='namelists.dat',status='old')
    read(57,nml=nampos)
    close(57)
    if(linner)then
       if(npplon==-666.or.npplat==-666.or.iminpp==-666.or.jminpp==-666)then!if not set in namelist
          npplon=klon_global-nbdpts-npbpts
          npplat=klat_global-nbdpts-npbpts
          iminpp=1+nbdpts+npbpts
          jminpp=1+nbdpts+npbpts
       endif
    else
       npplon=klon_global
       npplat=klat_global
       iminpp=1
       jminpp=1
    endif
    if(iminpp<=0.or.jminpp<=0.or.(iminpp+npplon-1)>klon_global.or.(jminpp+npplat-1)>klat_global)then
       print *,__FILE__,__LINE__
       call stop_program( '(i/j)minpp or nppl(at/on) is set outside domain')
    endif

    if(mype==0)then
       write(6,nml=nampos)
       write(6,nml=nampp)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=nampos)
       write(1,nml=nampp)
       close(1)
    endif

    allocate(files(nppstr))
    do j=1,nppstr
       call init_file(files(j),suff(j),nppstr,npplon,npplat,iminpp,jminpp)
    enddo
  end subroutine read_nampos



  subroutine postproc(klon,klat,klev,ksvar,&
       ktime, phim, & 
       surf_vars,&
       totcov   , cucov   , cw2d, acccw2d,                    &
       pblh,                                                 &
       kwmosv,  sacdg  , kacdg   ,nwmodg,  sacdg2 , kacdg2, &
       ci2d, accci2d, rad_cloud,  accqu2d, accqv2d, accvimfc2d, &
       accsunny, ice_thic,&
       accrunoff, accprl, accprc, tsmax,tsmin,t2max,t2min,maxdayprecip, &
       uv10max,accq2d,&
       accslwr,accsswr,acctlwr,acctswr,                  &
       accslwdn,accsswdn,acctswdn,                        &
       accsnow,                                             &
       accsnoc,accsnol,                                    &
       tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
       lcounttice,                                          &
       msvars,svar_surf,ksvars, svarsz,esvars,extrem_surf,  &
       acum_surf,mean_surf, meco, eco,lakes)                             

    use modddr
    use calendar
    use flake,only:laketype 
    use surface,only:surfVariables
    implicit none

    type(laketype),intent(in)::lakes
    integer,intent(in):: ksvar,klon,klat,klev
    type(time),intent(in)::ktime
    type(atm),intent(in)::phim
    type(surfVariables),intent(in)::surf_vars
    real(kind=realkind),pointer,dimension(:,:)::pfrl
    real(kind=realkind),pointer,dimension(:,:)::pfr_lake
    real(kind=realkind),intent(in)::totcov(klon,klat,klev), cucov(klon,klat,klev),        &
         cw2d(klon,klat), acccw2d(klon,klat)

    real(kind=realkind), intent(in)::pblh(klon,klat)
    real(kind=realkind),intent(in)::accsunny(klon,klat) 
    real(kind=realkind),intent(in)::ice_thic(klon,klat)

    real(kind=realkind),intent(in)::accrunoff(klon,klat),accprl(klon,klat) , accprc(klon,klat)
    real(kind=realkind),intent(in)::tsmax(klon,klat), tsmin(klon,klat),    &
         t2max(klon,klat), t2min(klon,klat),      &
         maxdayprecip(klon,klat),      &
         uv10max(klon,klat),                       &
         accq2d(klon,klat),                        &
         accslwr(klon,klat), accsswr(klon,klat),  &
         acctlwr(klon,klat), acctswr(klon,klat)
    real(kind=realkind),intent(in)::tcov(klon,klat),lcov(klon,klat),mcov(klon,klat),hcov(klon,klat),&
         acctcov(klon,klat),acclcov(klon,klat),accmcov(klon,klat),acchcov(klon,klat)
    real(kind=realkind),intent(in)::accslwdn(klon,klat), accsswdn(klon,klat)
    real(kind=realkind),intent(in)::accsnow(klon,klat),acctswdn(klon,klat)
    real(kind=realkind),intent(in)::accsnoc(klon,klat)
    real(kind=realkind),intent(in)::accsnol(klon,klat)

    real(kind=realkind),intent(in)::lcounttice(klon,klat)
    integer,intent(in)::ksvars
    real(kind=realkind),intent(in)::svarsz(klon,klat,ksvars)
    integer,intent(in)::msvars
    real(kind=realkind),intent(in)::svar_surf(klon,klat,msvars)
    integer,intent(in)::esvars
    real(kind=realkind),intent(in)::extrem_surf(klon,klat,esvars)
    real(kind=realkind),intent(in)::acum_surf(klon,klat,msvars)
    real(kind=realkind),intent(in)::mean_surf(klon,klat,msvars)
    integer,intent(in)::meco
    real(kind=realkind),target,intent(in)::eco(klon,klat,meco)
    real(kind=realkind),intent(in)::ci2d(klon,klat), accci2d(klon,klat), rad_cloud(klon,klat,klev)
    real(kind=realkind),intent(in)::accqu2d(klon,klat),accqv2d(klon,klat), &
         accvimfc2d(klon,klat)
    integer,intent(in)::kwmosv(ksvar)
    integer,intent(in)::kacdg ,nwmodg(kacdg)
    integer,intent(in)::kacdg2
    real(kind=realkind),intent(in)::sacdg(klon,klat,klev,kacdg)
    real(kind=realkind),intent(in)::sacdg2(klon,klat,kacdg2)

    integer j

    pfrl => eco(:,:,10)
    pfr_lake => eco(:,:,42)

    do 200 j=1,nppstr
       if(ktime==files(j)%write2fileTime)then
          !     extend fields with dummy values to boundary points in extend
          call extend(klon,klat,svar_surf(:,:,13)) !pu10
          call extend(klon,klat,svar_surf(:,:,19)) !pv10
          call extend(klon,klat,svar_surf(:,:,1))  !pt2m
          call extend(klon,klat,svar_surf(:,:,7))  !pq2m
          call extend(klon,klat,svar_surf(:,:,59)) !palb
          call extend(klon,klat,pblh)              !pblh


          call putdat(files(j),ktime, & 
                klon     ,klat   ,klev,ksvar,                             &
                surf_vars,pfrl,&
                phim, &
                rad_cloud ,cucov ,cw2d,acccw2d,            &
                pblh,                                               &
                kwmosv,                           &
                sacdg    ,kacdg  ,nwmodg,                           &
                sacdg2   ,kacdg2, &
                accsunny, ice_thic,                                &
                accrunoff,                                            &
                accprl, accprc,                                   &
                tsmax,tsmin,t2max,t2min,maxdayprecip,                        &
                uv10max,                                           &
                accq2d,                                            &
                accslwr,accsswr,acctlwr,acctswr,                &
                accslwdn,accsswdn,acctswdn,                      &
                ci2d,accci2d,                                           &
                accqu2d, accqv2d, accvimfc2d,                    &
                accsnow,                                           &
                accsnoc,accsnol,                                  &
                tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
                lcounttice,                                        &
                ksvars, svarsz,                                      &
                msvars,svar_surf,                                   &
                esvars,extrem_surf,                                 &
                acum_surf, mean_surf,                          &
                lakes,&
                meco,eco)

          files(j)%write2fileTime=files(j)%write2fileTime+files(j)%outdt

       endif
200 enddo
    return
  end subroutine postproc




  subroutine postpp(fileI, & 
       klev  , klon   , klat, ksvar  ,                                   &
       intime,                  &
       ah    , bh,                                               &
       phim,&
       totcov ,cucov ,cw2d, cwpath,    &
       pblh,                                                     &
       accsunny, ice_thick,                                      &
       accrunoff,                                                &
       accprl, accprc,                                           &
       tsmax,tsmin,t2max,t2min,maxdayprecip,                                  &
       uv10max,                                                  &
       accq2d,                                                   &
       accslwr,accsswr,acctlwr,acctswr,                          &
       ci2d,accci2d,                                                  &
       accqu2d,accqv2d,accvimfc2d,                               &
       accslwdn,accsswdn,acctswdn,                               &
       accsnow,                                                  &
       accsnoc,accsnol,                                          &
       tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
       lcounttice,                                               &
       ksvars, svars,                                           &
       msvars, svar_surf,                                       &
       esvars, extrem_surf,                                     &
       acum_surf, mean_surf,                               &
       lakes,&
       meco, eco,                                               &
       frland,&
       omf,omh,                                                  &
       pnwmosv,                                   &
       sacdg ,nsacdg ,nwmodg,                                    &
       sacdg2,nsacdg2, &
       etadot,surf_vars)

    use util
    use IO

    use confys
    use calendar
    use escom
    use modddr
    use surface,only:surfVariables
    use flake,only:laketype 
    implicit none
    type(surfVariables),intent(in),target::surf_vars
    type(laketype),intent(in),target::lakes
    type(file),intent(in)::fileI
    type(atm),intent(in),target::phim
    integer,intent(in):: klev,klon,klat

    type(time),intent(in)::intime
    real(kind=realkind),intent(in):: ah(klev+1),bh(klev+1)

    real(kind=realkind)::graviti
    real(kind=realkind),target,intent(in):: totcov(klon,klat,klev),         &
         cucov(klon,klat,klev),                &
         cwpath(klon,klat),  cw2d(klon,klat)

    real(kind=realkind),target,intent(in)::pblh(klon,klat),frland(klon,klat)

    real(kind=realkind),target::td(klon,klat,klev),td2m(klon,klat),rh2m(klon,klat)
    real(kind=realkind),intent(in)::omf(klon,klat,klev)
    real(kind=realkind),intent(in)::omh(klon,klat,klev+1)

    real(kind=realkind),allocatable,target::fstor2(:,:)
    real(kind=realkind),pointer::fstore(:,:)
    real(kind=realkind),allocatable::psl(:,:),fi(:,:,:),rh(:,:,:)

    integer,allocatable::ind(:,:)

    real(kind=realkind),allocatable::pretah(:,:,:),pretaf(:,:,:)
    real(kind=realkind),target,intent(in)::accsunny(klon,klat), ice_thick(klon,klat), &
         accrunoff(klon,klat),                       &
         accprl(klon,klat), accprc(klon,klat)         
    real(kind=realkind),target,intent(in)::tsmax(klon,klat), tsmin(klon,klat),         &
         t2max(klon,klat), t2min(klon,klat),maxdayprecip(klon,klat),         &
         uv10max(klon,klat),                         &
         accq2d(klon,klat),      &
         accslwr(klon,klat), accsswr(klon,klat),     &
         acctlwr(klon,klat), acctswr(klon,klat)
    real(kind=realkind),target,intent(in)::tcov(klon,klat),lcov(klon,klat),mcov(klon,klat),hcov(klon,klat),&
         acctcov(klon,klat),acclcov(klon,klat),accmcov(klon,klat),acchcov(klon,klat)

    real(kind=realkind),target::frsnow(klon,klat)
    real(kind=realkind),target,intent(in)::accsnow(klon,klat)

    !     add transient variables here
    real(kind=realkind),target,intent(in)::accsnoc(klon,klat)
    real(kind=realkind),target,intent(in)::accsnol(klon,klat)
    real(kind=realkind),target,intent(in)::ci2d(klon,klat),accci2d(klon,klat)
    real(kind=realkind),target,intent(in)::accqu2d(klon,klat),accqv2d(klon,klat),accvimfc2d(klon,klat)
    real(kind=realkind),target,intent(in)::accslwdn(klon,klat), accsswdn(klon,klat),acctswdn(klon,klat)

    real(kind=realkind)::covu(klon,klat,klev)
    real(kind=realkind) snopl_out(klon,klat), snfor_out(klon,klat)
    real(kind=realkind),target,intent(in)::lcounttice(klon,klat)
    integer,intent(in):: ksvars
    real(kind=realkind),target,intent(in)::svars(klon,klat,ksvars)
    integer,intent(in)::msvars
    real(kind=realkind),intent(in)::svar_surf(klon,klat,msvars)
    integer,intent(in)::esvars
    real(kind=realkind),intent(in)::extrem_surf(klon,klat,esvars)
    real(kind=realkind),intent(in)::acum_surf(klon,klat,msvars)
    real(kind=realkind),intent(in)::mean_surf(klon,klat,msvars)
    integer,intent(in)::meco
    real(kind=realkind),intent(in)::eco(klon,klat,meco)
    integer ilake

    integer,intent(in)::ksvar,pnwmosv(*)
    integer,intent(in)::nsacdg,nwmodg(*)
    integer,intent(in):: nsacdg2 
!    real(kind=realkind),intent(in)::svar(klon,klat,klev,ksvar)
    real(kind=realkind),target,intent(in):: sacdg(klon,klat,klev,nsacdg)
    real(kind=realkind),target,intent(in):: sacdg2(klon,klat,nsacdg2)
    real(kind=realkind),intent(in):: etadot(klon,klat,klev)

    !     GRIB FILE DDR


    logical relhum,geopot,presl,omegal
    real(kind=realkind) refpre
    integer jk,length,lev,jwmo,jsl,i,j
    integer iklow
    integer ik1(klon,klat),ik2(klon,klat)
    integer ikhigh2d(klon,klat),ikmean2d(klon,klat)
    integer indsac
    real(kind=realkind) zphigh,zpmean,zdhigh,zdmean,zdp
    real(kind=realkind) tdpre
    integer mype_out, inr
    real(kind=realkind) zsea

    logical:: lfirstopen(801:810)
    logical hit

    type(ddr)::myddr

    data lfirstopen/10*.true./
    save lfirstopen

    refpre=sip0
    allocate(fstor2(klon,klat))
    allocate(psl(klon,klat))
    allocate(ind(klon,klat))
    allocate(pretah(klon,klat,klev+1), pretaf(klon,klat,klev))
    allocate(fi(klon,klat,klev+1),rh(klon,klat,klev))
    relhum = .false.
    geopot = .false.
    presl  = .false.
    omegal = lomega
    indsac = 0




    !     DEFINE MODEL FULL LEVELS AND MODEL HALF LEVELS
    !     PRESSURES PRETAH & PRETAF


    do jk=1,klev+1
       do j=1,klat
          do i=1,klon
             if(jk<=klev) then
                pretaf(i,j,jk) = rcaDomain%afull(jk) + phim%ps(i,j)*rcaDomain%bfull(jk)
             endif
             pretah(i,j,jk) =   ah(jk) + phim%ps(i,j)*bh(jk)
          enddo
       enddo
    enddo

    !     Correct output of snopl and snfor
    do j=1,klat
       do i=1,klon
          snopl_out(i,j)=0._realkind
          if(svar_surf(i,j,53)>0._realkind)then
             snopl_out(i,j) = svars(i,j,8)
          endif
          snfor_out(i,j)=0._realkind
          if(svar_surf(i,j,54)>0._realkind)then
             snfor_out(i,j) = svars(i,j,9)
          endif
       enddo
    enddo


    call reset_ddr(myddr)  !     form ddr
    length = 0
    myddr%ndtbhl =(intime%year*100+intime%month)*100+intime%day
    myddr%nscbhl =(intime%hour*100+intime%min)*100+00 ! ss=0
    call adddtg(myddr%ndtbhl, myddr%nscbhl, 0, myddr%ndtvhl, myddr%nscvhl)
    myddr%nflshl = 0

    myddr%nlonhl = fileI%npplon
    myddr%nlathl = fileI%npplat

    myddr%nltphl = fileI%ltypml
    myddr%nlevhl = fileI%nlevml
    myddr%nrflhl = 0
    do lev=1,fileI%nlevml
       myddr%nlpthl(lev) = lev - 1
    enddo
    myddr%nmlfhl= fileI%nwmoml
    do jwmo=1,fileI%nwmoml !loopen over 
       myddr%nwmmhl(jwmo) = fileI%iwmoml(jwmo)!typen
       myddr%nmpthl(jwmo) =(jwmo-1)*fileI%nlevml !vilken level
    enddo

    myddr%nslfhl=fileI%nsl
    do jsl=1,fileI%nsl          !single level fields
       myddr%nwmshl(jsl)=fileI%iwmosl(jsl)
       myddr%nspthl(jsl)=fileI%nwmoml*fileI%nlevml+jsl-1
       myddr%nslthl(jsl)=fileI%ltypsl(jsl)
    enddo
    !geometric information

    myddr%dlonhl=RCAdomain%dlon
    myddr%dlathl=RCAdomain%dlat
    myddr%aweshl=RCAdomain%west+real(iminpp-1,realgribkind)*RCAdomain%dlon
    myddr%aeashl=RCAdomain%west+real(iminpp-1+myddr%nlonhl-1,realgribkind)*RCAdomain%dlon
    myddr%alafhl=RCAdomain%south+real(jminpp-1,realgribkind)*RCAdomain%dlat
    myddr%alalhl=RCAdomain%south+real(jminpp-1+myddr%nlathl-1,realgribkind)*RCAdomain%dlat
    myddr%aplohl=RCAdomain%polon
    myddr%aplahl=RCAdomain%polat

    if(fileI%ltypml==0.or.fileI%ltypml==109)then
       myddr%npplhl=2
       do lev=1,fileI%nlevml
          myddr%alevhl(lev,1)=rcaDomain%afull(fileI%alevml(lev))
          myddr%alevhl(lev,2)=rcaDomain%bfull(fileI%alevml(lev))
       enddo
    else
       myddr%npplhl=1
       do lev=1,fileI%nlevml
          myddr%alevhl(lev,1)=real(fileI%alevml(lev),realgribkind)
       enddo
    endif

    do jsl=1,fileI%nsl
       if(fileI%ltypsl(jsl)==0.or.fileI%ltypsl(jsl)==109)then
          myddr%alevhl(jsl,1) = rcaDomain%afull(fileI%alevsl(jsl))
          myddr%alevhl(jsl,2) = rcaDomain%bfull(fileI%alevsl(jsl))
       else
          myddr%alevhl(jsl,1) = real(fileI%alevsl(jsl),realgribkind)
       endif
    enddo


    !     OPEN/CLOSE THE GRIB FILE
    mype_out = 0
    if(mype == mype_out) then
       if(month_file)then
          
          if( fileI%prefix == 'fc'.or.fileI%prefix == 'cl' ) then
             if( fileI%sufix=='qq' .or. fileI%sufix=='pp' .or.  &
                  fileI%sufix=='hh' .or. fileI%sufix=='dd'.or. &
                  fileI%sufix=='  '.or.fileI%sufix=='ss'.or.fileI%sufix=='gb' ) then

                if(lfirstopen(fileI%unitnr).or. &
                     (intime%day==1.and.intime%hour==0.and.intime%min==0.and.intime%sec==0))then
                   if(.not. lfirstopen(fileI%unitnr)) then
                      call gwclos2(fileI%unitnr)
                   endif

                   call gwopen2(fileI%unitnr, &
                        intime%year,intime%month,intime%day,intime%hour,intime%min,length, &
                        fileI%prefix,fileI%sufix,myddr)
                   lfirstopen(fileI%unitnr)=.false.
                endif
             endif
!          elseif( fileI%prefix == 'cl' ) then
!             write(6,*) 'postpp: call gwopen2 cl..gb',fileI%unitnr
!             call gwopen2(fileI%unitnr,&
!                  intime%year,intime%month,intime%day,intime%hour,intime%min,length, &
!                  fileI%prefix,fileI%sufix,myddr)
!!$          else
!!$             call gwopen2(fileI%unitnr,&
!!$                  intime%year,intime%month,intime%day,intime%hour,intime%min,length, &
!!$                  fileI%prefix,fileI%sufix,myddr)
          endif
!!$       else
!!$          call gwopen2(fileI%unitnr, &
!!$               intime%year,intime%month,intime%day,intime%hour,intime%min,length, &
!!$               fileI%prefix,fileI%sufix,myddr)
       endif
    endif

    !     COMPUTATION OF DEW POINT TEMPERATURE

    !     DEFAULT VALUE T - 100

    do lev=1,klev
       do j=1,klat
          do i=1,klon
             tdpre=rcaDomain%afull(lev)+rcaDomain%bfull(lev)*phim%ps(i,j)
             if(phim%q(i,j,lev)>0._realkind) then
                td(i,j,lev)=es2t(q2e(tdpre,phim%q(i,j,lev)))
                if(td(i,j,lev)<0._realkind) then
                   td(i,j,lev)=phim%t(i,j,lev)-100._realkind
                endif
             else
                td(i,j,lev)=phim%t(i,j,lev)-100._realkind
             endif
          enddo
       enddo
    enddo

    !     write all the requested multilevel fields
    do  jwmo=1,fileI%nwmoml
       do  jk=1,fileI%nlevml
          if( fileI%ltypml==100 ) then 
             call interpolate2pressureLevel(filei%iwmoml(jwmo),real(filei%alevml(jk),realkind),  & 
                  fstor2,klon,klat,klev,            &
                  ah,bh,rcaDomain%afull,rcaDomain%bfull,pretah,pretaf,refpre, &
!ps111115                  RCAdomain%fis,phim%ps,phim%t,td,surf_vars%tsm,phim%u,phim%v,phim%q,phim%cw,totcov,   &
                  RCAdomain%fis,phim%ps,phim%t,td,svar_surf(:,:,63),phim%u,phim%v,phim%q,phim%cw,totcov,   &
                  rh,fi,psl,omf,omh,                &
                  phim%svar,ksvar,pnwmosv,               &
                  lphys, & 
                  relhum,geopot,presl,omegal)
              fstore=>fstor2

          elseif( filei%ltypml==0.or.filei%ltypml==109 ) then
             call xtreta(filei%iwmoml(jwmo),filei%alevml(jk), & !     extraction of eta level fields
                  fstor2,klon,klat,klev,                 &
                  pretah,pretaf,                         &
                  phim%ps,RCAdomain%fis,fi,phim%t,phim%u,phim%v,phim%q,omf,omh,             &
                  rh,phim%cw,totcov,cucov,                    &
                  phim%svar,ksvar,pnwmosv,                    &
                  sacdg,nsacdg,nwmodg,                   &
                  etadot,                                &
                  rcaDomain%afull,rcaDomain%bfull,                                 &
                  relhum,geopot,omegal)
             fstore=>fstor2
          else
             print *,__FILE__,__LINE__
             call stop_program( 'Not allowed ML')
          endif


          call gwwrloc(mype_out,fileI%unitnr,myddr%nltphl,fileI%iwmoml(jwmo),fileI%alevml(jk), &
               fstore,klon,klat,myddr)
       enddo
    enddo


    !     write all the requested single level fields
    do  jsl=1,fileI%nsl

       hit=.false.

       select case (fileI%ltypsl(jsl))
       case(0)!ltypsl(0)
          myddr%npplhl=2
          call xtreta(fileI%iwmosl(jsl),fileI%alevsl(jsl),fstor2,klon,klat, &
               klev,pretah,pretaf,phim%ps,RCAdomain%fis,fi,phim%t,phim%u,phim%v,phim%q,omf,omh,           &
               rh,phim%cw,totcov,cucov,phim%svar,ksvar,pnwmosv,sacdg,nsacdg,     &
               nwmodg,etadot,rcaDomain%afull,rcaDomain%bfull,relhum,geopot,omegal)
          fstore=>fstor2
          hit=.true.

       case(100)!ltypsl(100)
          myddr%npplhl=1
          call interpolate2pressureLevel(fileI%iwmosl(jsl),real(fileI%alevsl(jsl),realkind),fstor2,klon,klat,klev, &
!ps111115               ah,bh,rcaDomain%afull,rcaDomain%bfull,pretah,pretaf,refpre,RCAdomain%fis,phim%ps,phim%t,td,surf_vars%tsm,phim%u,phim%v,phim%q, &
               ah,bh,rcaDomain%afull,rcaDomain%bfull,pretah,pretaf,refpre,RCAdomain%fis,phim%ps,phim%t,td,svar_surf(:,:,63),phim%u,phim%v,phim%q, &
               phim%cw,totcov,rh,fi,psl,omf,omh,phim%svar,ksvar,pnwmosv,        &
               lphys, relhum,geopot,presl,omegal)
          fstore=>fstor2
          hit=.true.

       case(102)!ltypsl(102)
          myddr%npplhl=1

          select case(fileI%iwmosl(jsl))
          case(11)
             if(fileI%alevsl(jsl)==0)then!ltypsl(102),iwmosl(11),alevsl(0)sea surface temperature
                fstore=>surf_vars%tsea
                hit=.true.
             endif

          case(12)
             if(fileI%alevsl(jsl)==0) then!ltypsl(102),iwmosl(12),alevsl(0)
                !     extraction of fields at sea surface : sst + tsm
                !     sea surface temperature + tsm
                do j=1,klat
                   do i=1,klon
                      zsea=(1.0_realkind-frland(i,j))*(1.0_realkind-surf_vars%fri(i,j))
                      if(zsea<0.99_realkind) then 
!ps111115                         fstor2(i,j)=surf_vars%tsm(i,j)
                         fstor2(i,j)=svar_surf(i,j,63)
                      else
                         fstor2(i,j)=surf_vars%tsea(i,j)
                      endif
                   enddo
                enddo
                fstore=>fstor2
                hit=.true.
             endif

          case(83)
             if(fileI%alevsl(jsl)==0)then!ltypsl(102),iwmosl(83),alevsl(0)
                fstore=>surf_vars%rou
                hit=.true.
             endif

          case(91)
             if(fileI%alevsl(jsl)==0)then !ltypsl(102),iwmosl(91),alevsl(0)
                fstore=>surf_vars%fri
                hit=.true.
             endif

          case(92)
             if(fileI%alevsl(jsl)==0)then!ltypsl(102),iwmosl(92),alevsl(0)
                fstore=>ice_thick
                hit=.true.
             endif

          case(239)
             if(fileI%alevsl(jsl)==0)then!ltypsl(102),iwmosl(239),alevsl(0)
                fstore=>lcounttice
                hit=.true.
             endif

          case default
             print *,__FILE__,__LINE__
             write(6,*)' field missing in model data file'
             write(6,*)' wmo=',fileI%iwmosl(jsl),' lev.type=',&
                  fileI%ltypsl(jsl),' level=',fileI%alevsl(jsl)
             call stop_program(' field missing in model data file') 
          end select

       case(103)!ltypsl(103)
          myddr%npplhl=1

          !     extraction of fields at sea surface
          select case(fileI%iwmosl(jsl))
          case(1)!ltypsl(103),iwmosl(1)
             call interpolate2pressureLevel(fileI%iwmosl(jsl),&
                  real(fileI%alevsl(jsl),realkind),fstor2,klon,klat,klev, &
!111115                  ah,bh,rcaDomain%afull,rcaDomain%bfull,pretah,pretaf,refpre,RCAdomain%fis,phim%ps,phim%t,td,surf_vars%tsm,phim%u,phim%v,phim%q, &
                  ah,bh,rcaDomain%afull,rcaDomain%bfull,pretah,pretaf,refpre,RCAdomain%fis,phim%ps,phim%t,td,svar_surf(:,:,63),phim%u,phim%v,phim%q, &
                  phim%cw,totcov,rh,fi,psl,omf,omh,phim%svar,ksvar,pnwmosv,        &
                  lphys, relhum,geopot,presl,omegal)
             fstore=>fstor2
             hit=.true.

          case(11)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(103),iwmosl(11),alevsl(0)
                !     sea surface temperature
                fstore=>surf_vars%tsea
                hit=.true.
             endif

          case(12)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(103),iwmosl(12),alevsl(0)
                !     extraction of fields at sea surface : sst + tsm
                !     sea surface temperature + tsm
                do j=1,klat
                   do i=1,klon
                      zsea=(1.0_realkind-frland(i,j))*(1.0_realkind-surf_vars%fri(i,j))
                      if(zsea<0.99_realkind) then 
!ps111115                         fstor2(i,j)=surf_vars%tsm(i,j)
                         fstor2(i,j)=svar_surf(i,j,63)
                      else
                         fstor2(i,j)=surf_vars%tsea(i,j)
                      endif
                   enddo
                enddo
                fstore=>fstor2
                hit=.true.
             endif

          case(83)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(103),iwmosl(83),alevsl(0)
                !     roughness over sea
                fstore=>surf_vars%rou
                hit=.true.
             endif

          case(91)
             if(fileI%alevsl(jsl)==0)then!ltypsl(103),iwmosl(91),alevsl(0)
                fstore=>surf_vars%fri
                hit=.true.
             endif

          case(92)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(103),iwmosl(92),alevsl(0)
                fstore=>ice_thick
                hit=.true.
             endif

          case(239)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(103),iwmosl(239),alevsl(0)
                fstore=>lcounttice
                hit=.true.
             endif

          case default
             print *,__FILE__,__LINE__
             write(6,*)' field missing in model data file'
             write(6,*)' wmo=',fileI%iwmosl(jsl),' lev.type=',&
                  fileI%ltypsl(jsl),' level=',fileI%alevsl(jsl)
             call stop_program( 'old stop 890')
          end select



       case(105)!ltypsl(105)
          myddr%npplhl=1
          select case(fileI%iwmosl(jsl))
          case(1)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(1),alevsl(0)
                fstore=>phim%ps
                hit=.true.
             endif
          case(2)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(2),alevsl(0)
                do j=1,klat
                   do i=1,klon
                      fstor2(i,j) = real(idatastart-1+i,realkind) !i-index
                   enddo
                enddo
                fstore=>fstor2
                hit=.true.
             endif
          case(3)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(3),alevsl(0)
                do j=1,klat
                   do i=1,klon
                      fstor2(i,j) = real(jdatastart-1+j,realkind) !j-index
                   enddo
                enddo
                fstore=>fstor2
                hit=.true.
             endif
          case(4)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(4),alevsl(0)
                do j=1,klat
                   do i=1,klon
                      fstor2(i,j) = real(mype,realkind) !processor id
                   enddo
                enddo
                fstore=>fstor2
                hit=.true.
             endif
!!$          case(5)
!!$             if(fileI%alevsl(jsl)==0 ) then
!!$                do j=1,klat
!!$                   do i=1,klon
!!$                      fstore(i,j) = real(idatastart)
!!$                   enddo
!!$                enddo
!!$             endif

          case(6,7)
             if(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(6,7),alevsl(0)
                fstore=>RCAdomain%fis
                hit=.true.
                if(fileI%iwmosl(jsl)==7) then
                   graviti = 1.0_realkind/gravit
                   fstor2 = RCAdomain%fis*graviti
                   fstore=>fstor2
                   hit=.true.
                endif
             endif

          case(9)
             if(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(9),alevsl(0)
                fstore=>RCAdomain%orogsigm
                hit=.true.
             endif

          case(11)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(11),alevsl(0)
                fstore=>surf_vars%tsm
                hit=.true.
             elseif(fileI%alevsl(jsl)==999)then!ltypsl(105),iwmosl(11),alevsl(999)
                fstore=>surf_vars%tdm
                hit=.true.
             elseif(fileI%alevsl(jsl)==998)then!ltypsl(105),iwmosl(11),alevsl(998)
                fstore=>surf_vars%tsc
                hit=.true.
             endif

          case(15)
             if(fileI%alevsl(jsl)==2) then!ltypsl(105),iwmosl(15),alevsl(2)
                fstore=>t2max
                hit=.true.
             elseif(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(15),alevsl(0)
                fstore=>tsmax
                hit=.true.
             endif

          case(16)
             if(fileI%alevsl(jsl)==2)then!ltypsl(105),iwmosl(16),alevsl(2)
                fstore=>t2min
                hit=.true.
             elseif(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(16),alevsl(0)
                fstore=>tsmin
                hit=.true.
             endif

          case(17)
             if(fileI%alevsl(jsl)==2)then!ltypsl(105),iwmosl(17),alevsl(2)
                do j=1,klat
                   do i=1,klon
                      if(svar_surf(i,j,7)>0._realkind) then
                         td2m(i,j)=es2t(q2e(phim%ps(i,j),svar_surf(i,j,7)))
                         if(td2m(i,j)<0._realkind) then
                            td2m(i,j)=svar_surf(i,j,1)-100._realkind
                         endif
                      else
                         td2m(i,j)=svar_surf(i,j,1)-100._realkind
                      endif
                   enddo
                enddo
                fstore=>td2m
                hit=.true.
             endif

          case(32)
             if(fileI%alevsl(jsl)==10)then!ltypsl(105),iwmosl(32),alevsl(10)
                fstore=>uv10max
                hit=.true.
             endif

          case(52)
             if(fileI%alevsl(jsl)==2 ) then!ltypsl(105),iwmosl(52),alevsl(2)
                do j=1,klat
                   do i=1,klon
                      if(svar_surf(i,j,7)>0._realkind) then
                         rh2m(i,j)=q2e(phim%ps(i,j),svar_surf(i,j,7))/t2es(svar_surf(i,j,1))
                         if(rh2m(i,j)<0._realkind) then
                            rh2m(i,j)=0.01_realkind
                         endif
                      else
                         rh2m(i,j)=0.01_realkind
                      endif
                      !     remove values above 100%
                      rh2m(i,j)=min(rh2m(i,j),1.0_realkind)
                   enddo
                enddo
                fstore=>rh2m
                hit=.true.
             endif

          case(54)
             if(fileI%alevsl(jsl)==3006)then!ltypsl(105),iwmosl(54),alevsl(3006)
                !     averaged vertically integrated precipitable water
                fstore=>accq2d
                hit=.true.

             elseif(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(54),alevsl(0)
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2 = fstor2 +phim%q(:,:,jk)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore=>fstor2
                hit=.true.

             endif
          case(57)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(57),alevsl(0)
                !     acc. surface evaporation
                indsac=9
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif
          case(58)
             if(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(58),alevsl(0)
                !     cloud ice
                fstore=>ci2d
                hit=.true.
             elseif(fileI%alevsl(jsl)==3006)then!ltypsl(105),iwmosl(58),alevsl(3006)
                !     cloud ice
                fstore=>accci2d
                hit=.true.
             endif
          case(61)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(61),alevsl(0)
                !     total precipitation
                if(2>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = sacdg2(:,:,2) + sacdg2(:,:,1)
                fstore=>fstor2
                hit=.true.

             elseif(fileI%alevsl(jsl)==4006) then!ltypsl(105),iwmosl(61),alevsl(4006)
                !     total precipitation
                fstor2 = accprl + accprc
                fstore=>fstor2
                hit=.true.
             endif
          case(62)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(62),alevsl(0)
                !     large scale precipitation
                indsac=2
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.
             elseif(fileI%alevsl(jsl)==4006) then!ltypsl(105),iwmosl(62),alevsl(4006)
                !     acc. prl
                fstore=>accprl
                hit=.true.
             endif
          case(63)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(63),alevsl(0)
                !     convective precipitation
                indsac=1
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.
             elseif(fileI%alevsl(jsl)==4006) then!ltypsl(105),iwmosl(63),alevsl(4006)
                !     acc prc
                fstore=>accprc
                hit=.true.
             endif

          case(64)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     maximum hourly precipitation rate
                fstore=>maxdayprecip
                hit=.true.

             endif

          case(65)
             if(fileI%alevsl(jsl)==4006) then!ltypsl(105),iwmosl(65),alevsl(4006)
                !     acc snow
                fstore=>accsnow
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(65),alevsl(0)
                !     snow precipitation
                if(4>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = sacdg2(:,:,3)+sacdg2(:,:,4)
                fstore=>fstor2
                hit=.true.
             endif
          case(66)
             if(fileI%alevsl(jsl)==0 )then!ltypsl(105),iwmosl(66),alevsl(0)
                fstore=>surf_vars%snm
                hit=.true.
             endif
          case(67)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(67),alevsl(0)
                !     boundary layer height
                fstore=>pblh
                hit=.true.
             endif
          case(71)
             if(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(71),alevsl(0)
                !     2-dimensional cloud cover
                fstore=>tcov
                hit=.true.
             elseif(fileI%alevsl(jsl)==3006)then
                fstore=>acctcov
                hit=.true.
             elseif(fileI%alevsl(jsl)==2 ) then!ltypsl(105),iwmosl(71),alevsl(2)
                !     fog
                fstor2=totcov(:,:,klev)
                fstore => fstor2
                hit=.true.
             endif
          case(73)
             if(fileI%alevsl(jsl)==0)then
                fstore=>lcov
                hit=.true.
             elseif(fileI%alevsl(jsl)==3006)then
                fstore=>acclcov
                hit=.true.
             endif
          case(74)
             if(fileI%alevsl(jsl)==0)then
                fstore=>mcov
                hit=.true.
             elseif(fileI%alevsl(jsl)==3006)then
                fstore=>accmcov
                hit=.true.
             endif
          case(75)
             if(fileI%alevsl(jsl)==0)then
                fstore=>hcov
                hit=.true.
             elseif(fileI%alevsl(jsl)==3006)then
                fstore=>acchcov
                hit=.true.
             endif
          case(76)
             if(fileI%alevsl(jsl)==3006)then!ltypsl(105),iwmosl(76),alevsl(3006)
                !     vertically integrated cloud water
                fstore=>cwpath
                hit=.true.
             elseif(fileI%alevsl(jsl)==0)then!ltypsl(105),iwmosl(76),alevsl(0)
                !     vertically integrated cloud water
                fstore=>cw2d
                hit=.true.
             endif
          case(78)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(78),alevsl(0)
                !     accumulated convective snowfall
                indsac=3
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             elseif(fileI%alevsl(jsl)==4006 ) then!ltypsl(105),iwmosl(78),alevsl(4006)
                fstore=>accsnoc
                hit=.true.
             endif

          case(79)
             if(fileI%alevsl(jsl)==0 ) then!ltypsl(105),iwmosl(79),alevsl(0)
                !     accumulated stratiform snowfall
                indsac=4
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.
             elseif(fileI%alevsl(jsl)==4006)then!ltypsl(105),iwmosl(79),alevsl(4006)
                fstore=>accsnol
                hit=.true.
             endif

          case(80)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     sea surface temperature
                fstore=>surf_vars%tsea
                hit=.true.

             endif

          case(81)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     fraction of land
                fstore=>frland
                hit=.true.

             endif

          case(83)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     climatological roughness
                fstore=>surf_vars%roc
                hit=.true.

             endif

          case(86)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0)then
                !     soil wetness
                fstore=>surf_vars%swm
                hit=.true.

             elseif(fileI%alevsl(jsl)==999) then
                fstore=>surf_vars%sdm
                hit=.true.

             elseif(fileI%alevsl(jsl)==998 ) then
                !     climatological deep soil wetness
                fstore=>surf_vars%swc
                hit=.true.

             endif

          case(90)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==4006)then
                !     accumulated water runoff
                fstore=>accrunoff
                hit=.true.

             endif

          case(111)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. net surface shortwave radiation
                indsac=11
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             elseif(fileI%alevsl(jsl)==3006)then
                !     averaged surface short-wave radiation
                fstore=>accsswr
                hit=.true.

             endif

          case(112)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     averaged surface long-wave radiation
                fstore=>accslwr
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. net surface longwave radiation
                indsac=10
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif


          case(113)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     averaged top short-wave radiation
                fstore=>acctswr
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. net toa shortwave radiation
                indsac=15
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif
          case(114)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006) then
                !     averaged top long-wave radiation
                fstore=>acctlwr
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. net toa longwave radiation
                indsac=14
                if(indsac>nsacdg2)then
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif
          case(115)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006 ) then
                !     surface down long-wave radiation
                fstore=>accslwdn
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. downward surface longtwave radiation
                indsac=12
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif

          case(116)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. downward surface shortwave radiation
                indsac=13
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             elseif(fileI%alevsl(jsl)==3006)then
                !     surface down short-wave radiation
                fstore=>accsswdn
                hit=.true.

             endif

          case(117)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     toa down short-wave radiation
                fstore=>acctswdn
                hit=.true.

             elseif(fileI%alevsl(jsl)==0 ) then                                       
                !     acc. global toa radiation
                indsac=16
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif

          case(121)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006) then
                !     averaged  surface latent heat flux
                print *,__FILE__,__LINE__
                call stop_program('105,121,3006 not available now')
             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. surface latent heat flux
                indsac=6
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif
          case(122)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006 ) then
                !     averaged surface sensible heat flux
                print *,__FILE__,__LINE__
                call stop_program('105,122,3006 not available now')
             elseif(fileI%alevsl(jsl)==0 ) then
                !     acc. surface sensible heat flux
                indsac=5
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif
          case(124)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. surface u-momentum flux
                indsac=7
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.

             endif

          case(125)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. surface v-momentum flux
                indsac=8
                if(indsac>nsacdg2)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstore=>sacdg2(:,:,indsac)
                hit=.true.
             endif


          case(128)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006) then
                !     averaged  surface momentum flux
                print *,__FILE__,__LINE__
                call stop_program('105 128 3006 not available now')
             endif

          case(129)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0)then
                fstore=>frsnow
                hit=.true.

             endif
          case(151)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     vertically_integrated_zonal_moisture_flux
                fstore=>accqu2d
                hit=.true.

             endif

          case(152)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     vertically_integrated_meridional_moisture_flux
                fstore=>accqv2d
                hit=.true.

             endif

          case(153)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==3006)then
                !     vertically_integrated_moisture_flux_convergence
                fstore=>accvimfc2d
                hit=.true.

             endif


          case(197)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     fraction of forest
                fstore=>surf_vars%frf
                hit=.true.

             endif
          case(240)!ltypsl(105),iwmosl
             if(fileI%alevsl(jsl)==4006)then
                !     accum sunshine hours
                fstore=>accsunny
                hit=.true.

             endif
          case(241)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > esvars) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 241 105 xxx, ',inr,', is larger than allowed: ',esvars
                print *,'Increase esvars in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = extrem_surf(:,:,inr)
                fstore => fstor2
                hit=.true.

             endif

          case(242)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > msvars ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 242 105 xxx, ',inr,', is larger than allowed: ',msvars
                print *,'Increase msvars in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = mean_surf(:,:,inr)
                fstore => fstor2
                hit=.true.
             endif
          case(243)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             ilake =(inr-1) / lakes%lake_no_prog + 1
             inr = inr -(ilake-1)*lakes%lake_no_prog
             if( inr > lakes%lake_types*lakes%lake_no_prog ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 243 105 xxx, ',inr,', is larger than allowed: ',lakes%lake_types*lakes%lake_no_prog
                print *,'Increase lake_no_prog in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2  = lakes%prog_lakes(:,:,ilake,inr)
                fstore => fstor2
                hit=.true.
             endif

          case(244)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             ilake =(inr-1) / lakes%lake_no_diag + 1
             inr = inr -(ilake-1)*lakes%lake_no_diag
             if( inr > lakes%lake_types*lakes%lake_no_diag ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 244 105 xxx, ',inr,', is larger than allowed: ',lakes%lake_types*lakes%lake_no_diag
                print *,'Increase lake_no_diag in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = lakes%diag_lakes(:,:,ilake,inr)
                fstore => fstor2
                hit=.true.
             endif

          case(245)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > msvars ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 245 105 xxx, ',inr,', is larger than allowed: ',msvars
                print *,'Increase msvars in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = acum_surf(:,:,inr)
                fstore => fstor2
                hit=.true.
             endif
          case(246)!ltypsl(105),iwmosl
             !     full ECOCLIMAP field
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > meco) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 246 105 xxx, ',inr,', is larger than allowed: ',meco
                print *,'Increase meco in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = eco(:,:,inr)
                fstore => fstor2
                hit=.true.
             endif

          case(250)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > msvars ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 250 105 xxx, ',inr,', is larger than allowed: ',msvars
                print *,'Increase msvars in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = svar_surf(:,:,inr)
                fstore => fstor2
                hit=.true.
             endif

          case(252)!ltypsl(105),iwmosl
             fstor2 = 0._realkind
             inr = fileI%alevsl(jsl)
             if( inr > ksvars ) then
                print *,'------------------------------'
                print *,'Error in file and line number:'
                print *,__FILE__,__LINE__
                print *,'Number of requested variables for 252 105 xxx, ',inr,', is larger than allowed: ',ksvars
                print *,'Increase ksvars in gemini.F90 or reduce variables in namelists_namppp.dat'
                print *,'------------------------------'
                call stop_program('')
             else
                fstor2 = svars(:,:,inr)
                fstore => fstor2
                hit=.true.
             endif
          case default
             write(6,*)' field missing in model data file'
             write(6,*)' wmo=',fileI%iwmosl(jsl),' lev.type=',&
                  fileI%ltypsl(jsl),' level=',fileI%alevsl(jsl)
             print *,__FILE__,__LINE__
             call stop_program(' field missing in model data file')
          end select

       case(109)!ltypsl(109)
          myddr%npplhl=2
          call xtreta(fileI%iwmosl(jsl),fileI%alevsl(jsl),fstor2,klon,klat, &
               klev,pretah,pretaf,phim%ps,RCAdomain%fis,fi,phim%t,phim%u,phim%v,phim%q,omf,omh,           &
               rh,phim%cw,totcov,cucov,phim%svar,ksvar,pnwmosv,sacdg,nsacdg,     &
               nwmodg,etadot,rcaDomain%afull,rcaDomain%bfull,relhum,geopot,omegal)
          fstore => fstor2
          hit=.true.


       case(200)!ltypsl(200)
          myddr%npplhl=1
          !     * vertically integrated variables
          !     vertically integrated tendencies given by physics
          select case(fileI%iwmosl(jsl))
          case(210)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of temperature
                indsac=2         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2 = fstor2 + sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif

          case(211)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of temperature from v.diffn
                indsac=4         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2=fstor2+sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif

          case(212)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of temperature from radiation
                indsac=3         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2 = fstor2 + sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif


          case(220)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of spec humidity
                indsac=5         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2=fstor2+sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif

          case(221)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of spec. hum. from v.diffn
                indsac=6         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2=fstor2+sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif

          case(230)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of cloud condensate
                indsac=7         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2=fstor2+sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif

          case(231)!ltypsl(200),iwmosl
             if(fileI%alevsl(jsl)==0 ) then
                !     acc. tendency of cloud condensate from v.diffn
                indsac=8         !hard coding also in combig and phtask
                if(indsac>nsacdg)then 
                   print *,__FILE__,__LINE__
                   call stop_program( 'old stop 890')
                endif
                fstor2 = 0._realkind
                do jk=1,klev
                   fstor2=fstor2+sacdg(:,:,jk,indsac)*(pretah(:,:,jk+1)-pretah(:,:,jk))/gravit
                enddo
                fstore => fstor2
                hit=.true.
             endif
          case default
             write(6,*)' field missing from model '
             write(6,*)' wmo=',fileI%iwmosl(jsl),' lev.type=',     &
                  fileI%ltypsl(jsl),' level=',fileI%alevsl(jsl)
             print *,__FILE__,__LINE__
             call stop_program(' field missing from model ')
          end select
       end select    !ltyp

       !     WRITING OF THE FIELD
       if(hit)then
          call gwwrloc(mype_out,fileI%unitnr,fileI%ltypsl(jsl),fileI%iwmosl(jsl), &
               fileI%alevsl(jsl),fstore,klon,klat,myddr)
       else
          if(mype==0)then
             write(6,*)' field missing from model '
             write(6,*)' wmo=',fileI%iwmosl(jsl),' lev.type=',     &
                  fileI%ltypsl(jsl),' level=',fileI%alevsl(jsl)
          endif
          print *,__FILE__,__LINE__
          call stop_program(' field missing from model ')
       endif

    enddo !jsl


    deallocate(fstor2,psl)
    deallocate(ind)
    deallocate(pretah, pretaf)
    deallocate(fi,rh)
    return
  end subroutine postpp


  subroutine putdat(fileI, intime, & 
       klon    ,klat     ,klev,ksvar,&
       surf_vars,&
       frland , &
       phim, &
       totcov,cucov ,cw2d, cwpath,&
       pblh, &
       nwmosv,&
       sacdg    ,kacdg  ,nwmodg,&
       sacdg2   ,kacdg2, &
       accsunny, ice_thic,&
       accrunoff,&
       accprl, accprc,&
       tsmax,tsmin,t2max,t2min,maxdayprecip,&
       uv10max,&
       accq2d,&
       accslwr,accsswr,acctlwr,acctswr,&
       accslwdn,accsswdn,acctswdn,&
       ci2d,accci2d,&
       accqu2d, accqv2d, accvimfc2d,&
       accsnow,&
       accsnoc,accsnol,&
       tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
       lcounttice,&
       ksvars, svars,&
       msvars, svar_surf,&
       esvars, extrem_surf, &
       acum_surf, mean_surf, &
       lakes,&
       meco, eco)

    use  sl2tim4,only:omcomp
    use calendar
    use surface,only:surfVariables
    use flake,only:laketype 
    implicit none
    type(laketype),intent(in)::lakes
    type(surfVariables),intent(in)::surf_vars
    type(file)::fileI
    type(atm),intent(in)::phim
    type(time),intent(in)::intime
    integer,intent(in):: klon,klat,klev

    real(kind=realkind),intent(in)::frland(klon,klat) 
    real(kind=realkind),intent(in)::totcov(klon,klat,klev), cucov(klon,klat,klev), &
         cwpath(klon,klat), cw2d(klon,klat),                        &
         pblh(klon,klat)

    real(kind=realkind),allocatable,dimension(:,:,:)::omf,omh 

    real(kind=realkind):: zahalf(klev+1),zbhalf(klev+1)

    real(kind=realkind),intent(in)::accsunny(klon,klat), ice_thic(klon,klat)
    real(kind=realkind),intent(in)::accrunoff(klon,klat)
    real(kind=realkind),intent(in)::accprl(klon,klat), accprc(klon,klat)
    real(kind=realkind),intent(in)::accsnow(klon,klat)
    real(kind=realkind),intent(in)::accsnoc(klon,klat)
    real(kind=realkind),intent(in)::accsnol(klon,klat)

    real(kind=realkind),intent(in)::tsmax(klon,klat),tsmin(klon,klat),t2max(klon,klat),t2min(klon,klat)
    real(kind=realkind),intent(in)::maxdayprecip(klon,klat)
    real(kind=realkind),intent(in)::uv10max(klon,klat)
    real(kind=realkind),intent(in)::accq2d(klon,klat)
    real(kind=realkind),intent(in)::accslwr(klon,klat),accsswr(klon,klat),acctlwr(klon,klat),acctswr(klon,klat)
    real(kind=realkind),intent(in)::tcov(klon,klat),lcov(klon,klat),mcov(klon,klat),hcov(klon,klat),&
         acctcov(klon,klat),acclcov(klon,klat),accmcov(klon,klat),acchcov(klon,klat)

    real(kind=realkind),intent(in)::accslwdn(klon,klat),accsswdn(klon,klat),acctswdn(klon,klat)
    real(kind=realkind),intent(in)::ci2d(klon,klat),accci2d(klon,klat)
    real(kind=realkind),intent(in)::accqu2d(klon,klat),accqv2d(klon,klat), accvimfc2d(klon,klat)
    real(kind=realkind),intent(in)::lcounttice(klon,klat)

    integer,intent(in)::ksvars
    real(kind=realkind),intent(in)::svars(klon,klat,ksvars)
    integer,intent(in)::msvars
    real(kind=realkind),intent(in)::svar_surf(klon,klat,msvars)
    integer,intent(in)::esvars
    real(kind=realkind),intent(in)::extrem_surf(klon,klat,esvars)
    real(kind=realkind),intent(in)::acum_surf(klon,klat,msvars)
    real(kind=realkind),intent(in)::mean_surf(klon,klat,msvars)
    integer,intent(in)::meco
    real(kind=realkind),intent(in)::eco(klon,klat,meco)

    integer::jlev


    integer,intent(in)::ksvar,nwmosv(ksvar)
    integer,intent(in):: kacdg,nwmodg(kacdg)
    integer,intent(in):: kacdg2 
    real(kind=realkind),intent(in):: sacdg(klon,klat,klev,kacdg)
    real(kind=realkind),intent(in):: sacdg2(klon,klat,kacdg2)
    real(kind=realkind):: etadot(klon,klat,klev)
    real(kind=realkind) eps 

    allocate(omf(klon,klat,klev),omh(klon,klat,(klev+1)))

    !     MODEL GEOMETRY

    if(fileI%ltypml==100.or.fileI%ltypml==0.or.fileI%ltypml==109.or.fileI%ltypml==-1) then

       eps = 1.e-6_realkind

       zahalf = rcaDomain%ahyb
       zbhalf = rcaDomain%bhyb
       if(abs(rcaDomain%afull(1))<=eps .and.abs(rcaDomain%bfull(1))>eps ) then
          zahalf(1) = 0._realkind
          zbhalf(1) = rcaDomain%bhyb(2)*0.25_realkind
       elseif( abs(rcaDomain%afull(1))>eps .and.abs(rcaDomain%bfull(1))<=eps)then
          zahalf(1) = rcaDomain%ahyb(2)*0.25_realkind 
          zbhalf(1) = 0._realkind
       else
          write(6,*)' wrong definition of a:s and b:s!'
       endif

       if(lomega) then        !     compute omegas: omh and omf
          call omcomp(phim%u,phim%v,phim%ps,omf,omh,klon,klat,klev,zahalf,zbhalf,rcaDomain%afull,rcaDomain%bfull )
          call edcomp(klon,klat,klev,phim%ps,phim%u,phim%v,etadot,rcaDomain%ahyb,rcaDomain%bhyb)
       endif

       call postpp(fileI,klev,klon,klat,ksvar,   &
            intime,          &
            zahalf,zbhalf,                            &
            phim,& 
            totcov,cucov,cw2d,cwpath,               &
            pblh,                                    &
            accsunny, ice_thic,                      &
            accrunoff,                               &
            accprl, accprc,                          &
            tsmax,tsmin,t2max,t2min,maxdayprecip,    &
            uv10max,                                 &
            accq2d,                                  &
            accslwr,accsswr,acctlwr,acctswr,         &
            ci2d,accci2d,                                 &
            accqu2d,accqv2d,accvimfc2d,              &
            accslwdn,accsswdn,acctswdn,              &
            accsnow,                                 &
            accsnoc,accsnol,                         &
            tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
            lcounttice,                              &
            ksvars, svars,                           &
            msvars, svar_surf,                       &
            esvars, extrem_surf,                     &
            acum_surf, mean_surf,               &
            lakes,&
            meco, eco,                               &
            frland,&
            omf,omh,                                 &
            nwmosv,                       &
            sacdg,kacdg,nwmodg,                      &
            sacdg2,kacdg2,&
            etadot,surf_vars)
    else
       write(6,*) ' post processing not for this level type !',fileI%ltypml,fileI%sufix
       stop
    endif

    deallocate(omf,omh)
    return
  end subroutine putdat


  subroutine cloud_cover(klon,klat,klev,ps,cucov,totcov, &
       rad_cloud,tcov,lcov,mcov,hcov)

    implicit none
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(in)::totcov(klon,klat,klev), cucov(klon,klat,klev), &
         ps(klon,klat)
    real(kind=realkind),intent(out)::tcov(klon,klat),lcov(klon,klat),mcov(klon,klat),hcov(klon,klat)
    real(kind=realkind),intent(out)::rad_cloud(klon,klat,klev)
    real(kind=realkind)::pfull(klev)
    integer::jx,jy,jk
    integer:: phigh,pmid
    !    real(kind=realkind)::cov
    real(kind=realkind)::covu
    ! Limits for high and mid clouds
    !RCA3    real(kind=realkind),parameter::lbd_high = 45000._realkind
    !RCA3    real(kind=realkind),parameter::lbd_mid  = 80000._realkind
    real(kind=realkind),parameter::lbd_high = 44000._realkind
    real(kind=realkind),parameter::lbd_mid  = 68000._realkind

    !CUW110301beg calculate low, mid and high clouds. OBS Remove corresponding calculations in postprocess
    !new var phigh, pmid, cov, covu, tcov(x,y)
    rad_cloud =  max(0._realkind,min(1.0_realkind,(cucov+totcov)))
    do jy=1,klat
       do jx=1,klon
          phigh = 1
          pmid = 1
          do jk=1,klev
             pfull(jk)=RCAdomain%afull(jk)+RCAdomain%bfull(jk)*ps(jx,jy)
             if(pfull(jk)<=lbd_high)phigh=jk
             if(pfull(jk)<=lbd_mid)pmid=jk
          enddo
          !           cov=1.0_realkind
          covu=1.0_realkind
          do jk=2,phigh
             covu=covu*(1._realkind-max(rad_cloud(jx,jy,jk-1),&
                  rad_cloud(jx,jy,jk)))/&
                  (1._realkind-min(rad_cloud(jx,jy,jk-1),0.99_realkind))
             !             covu=1._realkind-cov  
             !             hcov(jx,jy)=covu
          enddo
          hcov(jx,jy)=1.0_realkind-covu
          !           cov=1.0_realkind
          covu=1.0_realkind
          do jk=phigh,pmid
             covu=covu*(1._realkind-max(rad_cloud(jx,jy,jk-1),&
                  rad_cloud(jx,jy,jk)))/&
                  (1._realkind-min(rad_cloud(jx,jy,jk-1),0.99_realkind))
             !             covu=1._realkind-cov  
             !             mcov(jx,jy)=covu
          enddo
          mcov(jx,jy)=1.0_realkind-covu
          !           cov=1.0_realkind
          covu=1.0_realkind
          do jk=pmid,klev
             covu=covu*(1._realkind-max(rad_cloud(jx,jy,jk-1),&
                  rad_cloud(jx,jy,jk)))/&
                  (1._realkind-min(rad_cloud(jx,jy,jk-1),0.99_realkind))
             !             covu=1._realkind-cov  
             !             lcov(jx,jy)=covu
          enddo
          lcov(jx,jy)=1.0_realkind-covu
          !           cov=1.0_realkind
          covu=1.0_realkind
          do jk=2,klev
             !             cov=covu*(1._realkind-max(rad_cloud(jx,jy,jk-1),&
             covu=covu*(1._realkind-max(rad_cloud(jx,jy,jk-1),&
                  rad_cloud(jx,jy,jk)))/&
                  (1._realkind-min(rad_cloud(jx,jy,jk-1),0.99_realkind))
             !             covu=1._realkind-cov  
             !             tcov(jx,jy)=covu
          enddo
          tcov(jx,jy)=1.0_realkind-covu
       enddo
    enddo
    !CUW110301end

    return
  end subroutine cloud_cover









  subroutine interpolate2pressureLevel(param, pres, fxtr, klon, klat, klev, &
       ah, bh, &
       af, bf, pretah, pretaf, refpre, fis, ps, t, td, ts, u, v, q, cw, &
       totcov, rh, fi, psl, omf, omh, svar, nsvar, nwmosv, lphys, &
       relhum, geopot, presl, omegal)

    use confys
    implicit none  
    !     subroutine interpolate2pressureLevel
    !     -----------------
    !
    !     purpose
    !     -------
    !     organization of the interpolation to pressure levels
    !     during post-processing
    !
    !     input parameters
    !     ----------------
    !     param  =  parameter to be extracted according to
    !     the wmo grib code convention
    !     pres   =  pressure of output pressure level
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) full model levels
    !     ah     =  a-values of model half levels
    !   (dimension klev+1)
    !     bh     =  b-values of model half levels
    !   (dimension klev+1)
    !     af     =  a-values of model full levels
    !   (dimension klev)
    !     bf     =  b-values of model full levels
    !   (dimension klev)
    !     pretah =  pressure of the model half levels
    !   (dimension klon,klat*(klev+1))
    !     pretaf =  pressure of the model full levels
    !   (dimension klon,klat*klev)
    !     refpre =  reference pressure of the model vertical
    !     coordinate system
    !     fis    =  surface geopotential(dimension klon,klat)
    !     ps     =  surface pressure(dimension klon,klat)
    !     t      =  temperature( dimension klon,klat*klev )
    !     td     =  dew point temperature( dimension klon,klat*klev )
    !     ts     =  surface temperature( dimension klon,klat )
    !     u      =  u-comp. of wind( dimension klon,klat*klev )
    !     v      =  v-comp. of wind( dimension klon,klat*klev )
    !     q      =  spec. humidity( dimension klon,klat*klev )
    !     cw     =  cloud water content( dimension klon,klat*klev )
    !     totcov =  total cloud cover( dimension klon,klat*klev )
    !     omf    =  part of omega given at model full levels
    !   (dimension klon,klat*klev)
    !     omh    =  part of omega given at model half levels
    !   (dimension klon,klat*(klev+1))
    !     svar   =  extra scalars
    !     nsvar  =  number of extra scalars
    !     nwmosv =  wmo-code for extra scalars
    !     lphys  =  .t. if physics varaiables are present( e.g. ts)
    !     .f. if physics varaiables are not present
    !     relhum    logical
    !     geopot    logical
    !     presl     logical
    !     omegal    logical
    !
    !     output parameters
    !     -----------------
    !     fxtr   =  output pressure level field
    !   (dimension klat*klon)
    !     rh     =  work field for model level relative humidity
    !   (may also be input, dimension klon,klat*klev)
    !     fi     =  work field for model level geopotential
    !   (may also be input, dimension klon,klat*(klev+1))
    !     psl    =  work field for sea-level pressure
    !   (may also be input, dimension klon,klat)
    integer,intent(in) :: param  
    integer,intent(in) :: klon, klat, klev  
    real(kind=realkind),intent(in) :: pres  
    real(kind=realkind) :: fxtr(klon,klat)  
    real(kind=realkind),intent(in)::fis(klon,klat), ps(klon,klat), t(klon,klat, klev), &
         td(klon,klat, klev), ts(klon,klat), u(klon,klat, klev), &
         v(klon,klat, klev), q(klon,klat, klev), rh(klon,klat,klev),&
         cw(klon,klat, klev), totcov(klon,klat, klev), &
         omf(klon,klat,klev),&
         omh(klon,klat, klev + 1),pretaf(klon,klat, klev), &
         pretah(klon,klat, klev + 1), &
         af(klev),ah(klev+1),bf(klev), bh(klev + 1)
    real(kind=realkind),intent(inout)::psl(klon,klat),fi(klon,klat, klev + 1)
    real(kind=realkind),intent(in):: refpre  
    logical,intent(in) :: lphys  
    logical :: do_cheating  
    integer,intent(in) :: nsvar, nwmosv(*)  
    logical,intent(in) :: omegal  
    real(kind=realkind),intent(in) :: svar(klon,klat, klev,nsvar)  
    logical:: relhum, geopot, presl
    integer ::  l
    do_cheating = .true.  
    !     mean sea-level pressure
    if(param==001) then  
       if(.not.presl) then  
          call pslcom(.false., ps, fis, t, ts, psl, klon, klat, klev, bf)
          presl = .true.  
       endif
       fxtr = psl

       !     geopotential
    elseif(param==006.or.param==007) then  
       if(.not.geopot) then  
          call getGeoPotentialFullLevel(fis,pretah,t,q,fi,klon,klat,klev)  
          geopot = .true.  
       endif
       call fipint(lphys, fi, pretah, ps, fis, ts, t, klon, klat, &
            klev, pres, fxtr, ah, bh, bf, refpre)
       if(param==007) then  
          fxtr = fxtr/gravit ! 9.80665  
       endif
       !     temperature
    elseif(param==011.or.param==013) then  
       call tpint(lphys, t, pretaf, ps, ts, fis, klon, klat, klev, &
            pres, fxtr, af, bf, refpre)
       if(param==013) then  
          fxtr = fxtr*(refpre/pres)**0.286_realkind  
       endif
    elseif(param==017) then  
       !     dew point temperature

       call interpolateModel2pressureLevel(pretaf, td, klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
    elseif(param==051) then  
       !     specific humidity

       call interpolateModel2pressureLevel(pretaf, q, klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
    elseif(param==052) then  
       !     relative humidity
       if(.not.relhum) then  
          !rh is output
          call qtorh(q, rh, t, ps, af, bf, klon, klat, klev)  
          relhum = .true.  
       endif

       call interpolateModel2pressureLevel(pretaf, rh, klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
       if(do_cheating) then  
          !     values between 0. and 1.0_realkind
          fxtr = max(min(fxtr,1.0_realkind),0.0_realkind)  
       endif
    elseif(param==033) then  
       !     wind components

       call interpolateModel2pressureLevel(pretaf, u, klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
    elseif(param==034) then  
       !     wind components

       call interpolateModel2pressureLevel(pretaf, v, klon, klat, klev, pres, fxtr, 1.0_realkind, af,    bf, refpre)
    elseif(param==003) then  
       !     surface pressure tendency
       if(.not.omegal) then  
          print *,__FILE__,__LINE__
          write(6,*) ' omega not available for extraction'  
          write(6,*) ' of dpsdt!'  
          call stop_program( ' omega not available for extraction'  ) 
       endif

       call xtr(omh(:,:,klev+1),klon,klat,fxtr)
    elseif(param==039) then  
       !     vertical velocity
       if(.not.omegal) then  
          print *,__FILE__,__LINE__
          write(6,*) ' omega not available for vertical interpolation!'  
          call stop_program(' omega not available for vertical interpolation!')  
       endif
       call omint(klon, klat, klev, omf, pretaf, omh, pretah, pres,   fxtr, af, bf, ah, bh, refpre)
    elseif(param==076) then  
       !     cloud water
       call interpolateModel2pressureLevel(pretaf, cw, klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
    elseif(param==071) then  
       !     cloud cover
       call interpolateModel2pressureLevel(pretaf, totcov, klon, klat, klev, pres, fxtr, 1.0_realkind,    af, bf, refpre)
    else  
       do l = 1, nsvar  
          if(param==nwmosv(l) )then !extra scalars
             call interpolateModel2pressureLevel(pretaf, svar(:,:,:,l), klon, klat, klev, pres, fxtr, 1.0_realkind, af, bf, refpre)
             return  
          endif


       enddo
       write(6,*) 'pressure-field storing not introduced for param=', &
            & param, ' and level=', pres
       print *,__FILE__,__LINE__
       call stop_program(  'pressure-field storing not introduced for param=')
       return  

    endif
  end subroutine interpolate2pressureLevel






  subroutine interpolateModel2pressureLevel(pretaf, x, klon, klat, klev, pres, xpres, scale,af,bf,refpre)

    implicit none  
    !     subroutine interpolateModel2pressureLevel
    !     ----------------
    !
    !     purpose
    !     -------
    !     vertical interpolation from model levels to
    !     a pressure level
    !
    !     input parameters
    !     ----------------
    !     pretaf =  pressure of the model full levels
    !   (dimension klon,klat*klev)
    !     x      =  input field to be interpolated to a pressure
    !     level(dimension klon,klat*klev)
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) model levels
    !     pres   =  pressure of the output pressure level
    !     numlon =  number of output grid-points in the x-direction
    !     numlat =  number of output grid-points in the y-direction
    !     scale  =  scaling factor to be applied to the output field
    !     af     =  a-values of the model vertical coordinate
    !   (dimension klev)
    !     bf     =  b-values of the model vertical coordinate
    !   (dimension klev)
    !     refpre =  reference pressure of the model vertical
    !     coordinate

    !     ind    =  index field giving references of the gridpoints in
    !     the output field to the gridpoints of the
    !     input field
    !
    !     output parameters
    !     -----------------
    !     xpres  =  output pressure level field
    !   (dimension numlat*numlon)
    !
    !----------------------------------------------------------------------
    !
    integer,intent(in) :: klon, klat, klev  
    real(kind=realkind),intent(in) :: pretaf(klon, klat, klev), x(klon, klat, klev), pres, &
         refpre,  scale
    real(kind=realkind),intent(inout)::xpres(klon, klat)
    real(kind=realkind),intent(in) :: af(klev), bf(klev)  

    !     0.1 local work space and local variables

    real(kind=realkind) :: etaf(klev), aletaf(klev), daleta(klev), zaleta  
    integer :: jx, jy, jk  
    !     extract eta levels information and
    !     compute vertical interpolation factors
    do 200 jk = 1, klev  
       etaf(jk) = af(jk) / refpre+bf(jk)  
       aletaf(jk) = log(etaf(jk) )  
200 enddo
    do 203 jk = 2, klev  
       daleta(jk - 1) = 1.0_realkind /(aletaf(jk) - aletaf(jk - 1) )  
203 enddo
    !     extrapolation above the highest eta level
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres<=pretaf(jx, jy, 1) ) then  
             xpres(jx, jy) = x(jx, jy, 1)  
          endif
       enddo


    enddo
    !     interpolation between two eta levels
    do 410 jk = 2, klev  
       do jy = 1, klat  
          do jx = 1, klon  
             if(pres>pretaf(jx, jy, jk - 1) .and.pres<=pretaf(jx, &
                  jy, jk) ) then

                zaleta = log(etaf(jk - 1) +(pres - pretaf(jx, jy, jk - &
                     1) ) /(pretaf(jx, jy, jk) - pretaf(jx, jy, jk - 1) ) &
                     *(etaf(jk) - etaf(jk - 1) ) )

                xpres(jx, jy) = x(jx, jy, jk) +(aletaf(jk) - zaleta) &
                     * daleta(jk - 1) *(x(jx, jy, jk - 1) - x(jx, jy, jk) )
             endif
          enddo
       enddo
410 end do

    !     extrapolation below the lowest model level
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres>=pretaf(jx, jy, klev) ) then  
             xpres(jx, jy) = x(jx, jy, klev)  
          endif
       enddo
    enddo
    !     carry out required re-scaling of output
    do jy = 1, klat  
       do jx = 1, klon  
          xpres(jx, jy) = xpres(jx, jy) * scale  
       enddo
    enddo
    return  
  end subroutine interpolateModel2pressureLevel




  subroutine xtreta(param, level, fxtr, klon, klat, klev, pretah, &
       pretaf, ps, fis, fi, t, u, v, q, omf, omh, rh, cw, totcov, cucov, &
       svar, nsvar, nwmosv, sacdg, nsacdg, nwmodg, etadot, af, bf, &
       relhum, geopot, omegal)
    use interpol

    implicit none  
    !     subroutine xtreta
    !     purpose
    !     -------
    !     organization of the extraction of model level data
    !     during postprocessing
    !
    !     input parameters
    !     ----------------
    !     param  =  parameter to be extracted according to
    !     the wmo grib code convention
    !     level  =  model level number of the field to be
    !     extracted
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) model levels
    !     ps     =  surface pressure of dimension klon,klat
    !     fis    =  surface geopotential
    !     fi     =  geopotential
    !     t      =  temperature of dimension klon,klat*klev
    !     u      =  u-component of wind(dim. klon,klat*klev)
    !     v      =  v-component of wind(dim. klon,klat*klev)
    !     q      =  specific humidity(dim. klon,klat*klev)
    !     omf    =  contribution of full levels to vertical velocity
    !     omh    =  contribution of half levels to vertical velocity
    !     cw     =  cloud water(dim. klon,klat*klev)
    !     totcov =  total cloud cover(dim. klon,klat*klev)
    !     cucov  =  convective cloud cover(dim. klon,klat*klev)
    !     svar   =  extra scalars
    !     nsvar  =  number of extra scalars
    !     nwmosv =  wmo-code for extra scalars
    !     etadot =  hybrid vertical velocity( dimension klon,klat*klev )
    !     f    sacdg(1,l): accumulated tendencies at model levels
    !     af     =  a-values of model full levels(dimension klev)
    !     bf     =  b-values of model full levels(dimension klev)
    !     ind    =  index field giving references of the gridpoints in
    !     the output field to the gridpoints of the
    !     input field
    !
    !     output parameters
    !     -----------------
    !     fxtr   =  extracted field of dimension klon,klat
    !     rh     =  work filed of dim. klon,klat*klev used for
    !     storing of relative humidity
    !     ok     =  .true. if field was extracted
    !     .false. if field was not present(error in
    !     selection parameters)
    !
    !------------------------------------------------------------------
    !
    integer,intent(in) :: param, level, klon, klat, klev  
    real(kind=realkind),intent(inout) :: fxtr(klon,klat)  
    real(kind=realkind)::pretah(klon,klat,klev+1), pretaf(klon,klat,klev)  
    real(kind=realkind),intent(inout)::fi(klon,klat,klev+1)
    real(kind=realkind),intent(in)::ps(klon,klat), fis(klon,klat), &
         t(klon,klat,klev), u(klon,klat, klev), v(klon,klat,klev),&
         q(klon,klat,klev), omf(klon,klat,klev), &
         omh(klon,klat,klev+1),rh(klon,klat,klev),cw(klon,klat,klev),&
         totcov(klon,klat,klev), cucov(klon,klat, klev) &
         , af(klev), bf(klev)
    integer,intent(in) :: nsvar, nwmosv(*)  
    real(kind=realkind),intent(in) :: svar(klon,klat, klev,nsvar)  
    logical :: do_cheat  
    integer,intent(in):: nsacdg, nwmodg(*)  
    real(kind=realkind),intent(in)::sacdg(klon,klat,klev,nsacdg),etadot(klon,klat,klev)
    real(kind=realkind) :: zomega(klon,klat), zfi(klon,klat)  
    integer :: i,j, l
    logical,intent(in) ::   omegal  
    logical::relhum, geopot
    do_cheat = .true.  
    !     geopotential
    if(param==006) then  
       if(.not.geopot) then  
          call getGeoPotentialFullLevel(fis, pretah, t, q, fi, klon, klat, klev)  
          geopot = .true.  
       endif
       call vint(klev + 1, klon * klat, fi, pretah, zfi,pretaf(:,:,level) )
       call xtr(zfi, klon, klat, fxtr)
       return  
    endif
    !     temperature
    if(param==011.or.param==013) then  
       call xtr(t(:,:,level), klon, klat, fxtr)
       if(param==013) then  
          call xtr(pretaf(:,:, level), klon, klat, zfi)
          do j = 1,klat  
             do i=1,klon
                fxtr(i,j) = fxtr(i,j) *(1.e+5_realkind / zfi(i,j) ) **0.286_realkind  
             enddo
          enddo
       endif
       return  
    endif
    !     u-component of wind
    if(param==033) then  
       call xtr(u(:,:,level), klon, klat, fxtr)
       return  

    endif
    !     v-component of wind
    if(param==034) then  
       call xtr(v(:,:,level), klon, klat, fxtr)
       return  

    endif
    !     vertical velocity
    if(param==039) then  
       call vint(klev + 1, klon * klat, omh, pretah, zomega, pretaf(:,:,level))
       do 150 j = 1,klat 
          do i=1,klon
             zomega(i,j) = zomega(i,j) + omf(i,j, level)  
          enddo
150    enddo
       call xtr(zomega, klon, klat, fxtr)
       return  


    endif
    !     water vapour mixing ration
    if(param==051) then  
       call xtr(q(:,:, level), klon, klat, fxtr)
       return  


    endif
    !     relative humidity
    if(param==052) then  
       if(.not.relhum) then  
          call qtorh(q, rh, t, ps, af, bf, klon, klat, klev)  
          relhum = .true.  

       endif
       call xtr(rh(:,:, level), klon, klat, fxtr)
       if(do_cheat) then  
          !     *        keep values between 0. and 1.
          do j = 1,  klat  
             do i=1,klon
                fxtr(i,j) = min(fxtr(i,j), 1.0_realkind)  
                fxtr(i,j) = max(fxtr(i,j), 0._realkind)  
             enddo
          enddo

       endif
       return  
    endif
    !     cloud water
    if(param==076) then  
       call xtr(cw(:,:, level), klon, klat, fxtr)
       return  

    endif
    !     total cloud cover
    if(param==071) then  
       call xtr(totcov(:,:, level), klon, klat, fxtr)
       return  


    endif
    !     convective cloud cover
    if(param==072) then  
       call xtr(cucov(:,:, level), klon, klat, fxtr)
       return  

    endif
    !     *    extra scalars
    do l = 1, nsvar  
       if(param==nwmosv(l) ) then  
          call xtr(svar(:,:, level, l), klon, klat, fxtr)
          return  
       endif
    enddo
    !     *    "diagnostics" accumulated 3-d fields
    do l = 1, nsacdg  
       if(param==nwmodg(l) ) then  
          call xtr(sacdg(:,:, level, l), klon, klat, fxtr)
          return  
       endif
    enddo
    !     *    hybrid vertical velocity
    if(param==040) then  
       if(.not.omegal) then  
          print *,__FILE__,__LINE__
          write(6,*) 'etadot not available'  
          write(6,*) 'choose lomega=.true. in input namelist'  
          call stop_program( 'etadot not available') 
       endif
       call xtr(etadot(:,:, level), klon, klat, fxtr)
       return  
    endif
    write(6,*) 'model-field storing not introduced for param=', &
         param, ' and level=', level
    print *,__FILE__,__LINE__
    call stop_program( 'model-field storing not introduced for param=') 
  end subroutine xtreta





  subroutine xtr(fin, klon, klat, fout)
    implicit none  
    !     extraction of grid-point values in a sub-area of a
    !     larger area
    !
    !     input parameters
    !     ----------------
    !     klon = number of gridpoints in the x-dir. of input field
    !     klat = number of gridpoints in the y-dir. of input field
    !     fin  = input field of dimension klon,klat
    !     klon = number of gridpoints in the x-dir. of output field
    !     klat = number of gridpoints in the y-dir. of output field
    !     fout = output field of dimension klon,klat
    !     ind  = index(pointer) field of dimension klon,klat
    !            containing indeces of grid-points in output field
    !            to grid-points in the input field
    integer,intent(in) :: klon, klat
    real(kind=realkind),intent(in) ::fin(klon,klat)  
    real(kind=realkind)::fout(klon,klat)  
    fout = fin
    return  
  end subroutine xtr




  subroutine qtorh(q,rh,t,ps,afull,bfull,klon,klat,klev)
    use confys
    use escom
    implicit none

    !     conversion of specific humidity into relative humidity
    !     for model level data

    !     input parameters
    !     ----------------
    !     klon   =  number of gridpoints in the x-direction
    !     klat   =  number of gridpoints in the y-direction
    !     klev   =  number of vertical levels
    !     q      =  specific humidity of dimension klon,klat*klev
    !     t      =  temperature of fimension klon,klat*klev
    !     ps     =  surface pressure of dimension klon,klat
    !     afull  =  a-values(dimension klev) defing model
    !     full levels
    !     bfull  =  b-values(dimension klev) defing model
    !     full levels
    !     
    !     output parameters
    !     -----------------
    !     rh     =  relative humidity of dimension klon,klat

    integer,intent(in):: klon,klat,klev
    real(kind=realkind):: q(klon,klat,klev),rh(klon,klat,klev)
    real(kind=realkind),intent(in)::t(klon,klat,klev),ps(klon,klat)
    real(kind=realkind),intent(in):: afull(klev), bfull(klev)

    integer k,i,j

    do  k=1,klev
       do j=1,klat
          do i=1,klon
             rh(i,j,k) = (afull(k) + ps(i,j)*bfull(k))*q(i,j,k)/(epsilo*esat( t(i,j,k) ))
             if( rh(i,j,k)<0._realkind ) rh(i,j,k) = 0._realkind
          enddo
       enddo
    enddo

    return
  end subroutine qtorh




  subroutine fipint(lphys, fi, pretah, ps, fis, ts, t, klon, klat, &
       klev, pres, fipres, ah, bh, bf, refpre)
    use confys  
    implicit none  
    !
    !---------------------------------------------------------------------
    !
    !     subroutine fipint
    !     -----------------
    !
    !     purpose:
    !     --------
    !     vertical interpolation of geopotential from hybrid model
    !     levels to pressure levels
    !
    !     input parameters
    !     ----------------
    !     lphys  =  .t. if physics variables(ts) are available
    !     fi     =  model half level geopotentials
    !             ( dimension klon*klat*(klev+1) )
    !     pretah =  pressure of the model half levels
    !             (dimension klon,klat*(klev+1))
    !     fis    =  surface geopotential(dimension klon,klat)
    !     ps     =  surface pressure(dimension klon,klat)
    !     ts     =  surface temperature( dimension klon,klat )
    !     t      =  temperature( dimension klon,klat*klev )
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) full model levels
    !     pres   =  pressure of output pressure level
    !     numlon =  number of output grid-points in the x-direction
    !     numlat =  number of output grid-points in the y-direction
    !     ah     =  a-values of model half levels
    !             (dimension klev+1)
    !     bh     =  b-values of model half levels
    !             (dimension klev+1)
    !     bf     =  b-values of model full levels
    !             (dimension klev)
    !     refpre =  reference pressure of the model vertical
    !               coordinate system
    !     ind    =  index field giving references of the gridpoints in
    !               the output field to the gridpoints of the
    !               input field
    !
    !     output parameters
    !     -----------------
    !     fipres =  output pressure level geopotential field
    !             (dimension numlat*numlon)
    !
    !*    0.1 input and output variables
    !
    logical,intent(in)::lphys  
    integer,intent(in)::klon, klat, klev  
    real(kind=realkind),intent(in)::fi(klon, klat, klev + 1), pretah(klon, klat, klev + 1), &
         ps(klon, klat), fis(klon, klat), ts(klon, klat), t(klon, klat, &
         klev), pres, ah(klev + 1), bh(klev + 1), &
         bf(klev), refpre
    real(kind=realkind),intent(inout):: fipres(klon, klat)

    ! work space and local variables
    real(kind=realkind) :: etah(klev + 1), aletah(klev + 1), dalehi(klev)  
    integer :: jx, jy, jk  
    real(kind=realkind) :: aleta, tstar, alfa, tzero,  alsint  

    real(kind=realkind),parameter::tlapse = 0.0065_realkind  
    ! initalize vertical integration weights

    do 200 jk = 1, klev + 1  
       etah(jk) = ah(jk) / refpre+bh(jk)  
       aletah(jk) = log(etah(jk) )  
200 enddo
    !
    do 202 jk = 2, klev + 1  
       dalehi(jk - 1) = 1.0_realkind /(aletah(jk) - aletah(jk - 1) )  
202 enddo
    !      loop over all horizontal gridpoints and do, if needed,
    !         extrapolation above the highest model half level

    do jy = 1, klat  
       do jx = 1, klon  
          if(pres<=pretah(jx, jy, 1) ) then  
             aleta = log(etah(1) +(etah(2) - etah(1) ) *(pretah(jx, &
                  jy, 1) - pres) /(pretah(jx, jy, 1) - pretah(jx, jy, 2) ) )
             fipres(jx, jy) = fi(jx, jy, 1) +(aletah(1) - aleta) &
                  * dalehi(1) *(fi(jx, jy, 1) - fi(jx, jy, 2) )
          endif
       enddo
    enddo
    !      loop over all horizontal gridpoints and do, if needed,
    !         interpolation between two model half levels

    do 410 jk = 2, klev + 1  
       do jy = 1, klat  
          do jx = 1, klon  
             if(pres>pretah(jx, jy, jk - 1) .and.pres<=pretah(jx, &
                  jy, jk) ) then
                aleta = log(etah(jk - 1) +(pres - pretah(jx, jy, jk - 1) &
                     ) *(etah(jk) - etah(jk - 1) ) /(pretah(jx, jy, jk) &
                     - pretah(jx, jy, jk - 1) ) )
                fipres(jx, jy) = fi(jx, jy, jk) +(aletah(jk) - aleta) &
                     * dalehi(jk - 1) *(fi(jx, jy, jk - 1) - fi(jx, jy, jk) )
             endif
          enddo
       enddo
410 end do
    !     loop over all horizontal gridpoints and do, if needed,
    !         extrapolation below the ground
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres>pretah(jx, jy, klev + 1) ) then  
             !*       5.1 extract the  surface temperature
             if(lphys) then  
                tstar = ts(jx, jy)  
             else  
                tstar = t(jx, jy, klev) *(1.0_realkind + tlapse * rair / gravit * &
                     (1.0_realkind / bf(klev) - 1.0_realkind) )
             endif
             !*           5.1.1 determine temperature reduction factor
             if(tstar<255.0_realkind) then  
                !*               5.1.1.1  prevent too high pressures
                !*                        under cold terrain
                tstar = 0.5_realkind *(tstar + 255._realkind)  
                alfa = tlapse * rair / gravit  
             else  
                !*                5.1.1.2 trial extrapolation of
                !*                        of sea-level temperature
                tzero = tstar + tlapse * fis(jx, jy) / gravit  
                if(tzero>290.5_realkind) then  
                   !*                    5.1.1.2.1  prevent too low pressures
                   !*                               under hot terrain
                   if(tstar<=290.5_realkind) then  
                      alfa = rair *(290.5_realkind - tstar) / fis(jx, jy)  
                   else  
                      alfa = 0._realkind  
                      tstar = 0.5_realkind *(290.5_realkind + tstar)  
                   endif
                else  
                   !*                    5.1.1.2.2  normal temp.reduction factor
                   alfa = tlapse * rair / gravit  
                   !
                endif
                !
             endif
             !*       5.2 extrapolation of geopotential below the surface
             alsint = log(pres / ps(jx, jy) )  
             fipres(jx, jy) = fis(jx, jy) - rair * tstar * alsint * &
                  (1._realkind + 0.5_realkind * alfa * alsint + 1._realkind / 6._realkind *(alfa * alsint) **2.0_realkind)
          endif
       enddo
    enddo
    return  
  end subroutine fipint






  subroutine getGeoPotentialFullLevel(fis, pretah, t, q, fi, klon, klat, klev)  
    use confys  
    implicit none  
    !
    !--------------------------------------------------------------------
    !
    !     subroutine getGeoPotentialFullLevel
    !     -----------------
    !
    !     purpose
    !     -------
    !     calculation of geopotential of model full levels
    !
    !     input parameters
    !     ----------------
    !     fis    =  surface geopotential(dim.klon,klat)
    !     pretaf =  pressure of the model full levels
    !             (dimension klon,klat*klev)
    !     t      =  temperature of the model full levels
    !             (dimension klon,klat*klev)
    !     q      =  specific humidity of the model full levels
    !             (dimension klon,klat*klev)
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) model levels
    !     numlon =  number of output grid-points in the x-direction
    !     numlat =  number of output grid-points in the y-direction
    !     ind    =  index field giving references of the gridpoints in
    !               the output field to the gridpoints of the
    !               input field(dimension numlon*numlat)
    !
    !     output parameters
    !     -----------------
    !     fi     =  output geopotential on model half levels
    !             (dimension numlat*numlon*(klev+1))
    !
    !-----------------------------------------------------------------------
    !
    integer,intent(in)::klon, klat, klev  
    real(kind=realkind),intent(in)::fis(klon, klat), pretah(klon, klat, klev + 1),&
         t(klon,klat, klev), q(klon, klat, klev)
    real(kind=realkind),intent(inout)::fi(klon,klat,klev+1)
    integer :: jx, jy, jk  
    real(kind=realkind) :: tv, zepsm1  

    zepsm1 = 1.0_realkind / epsilo - 1.0_realkind  

    !*     2.0 surface value of geopotential
    do jy = 1, klat  
       do jx = 1, klon  
          fi(jx, jy, klev + 1) = fis(jx, jy)  
       enddo
    enddo

    !*      3.0   layer by layer summation
    do 310 jk = klev, 1, - 1  
       do jy = 1, klat  
          do jx = 1, klon  
             tv =(1.0_realkind + zepsm1 * q(jx, jy, jk) ) * t(jx, jy, jk)  
             fi(jx, jy, jk) = fi(jx, jy, jk + 1) + log(pretah(jx, jy, &
                  jk + 1) / pretah(jx, jy, jk) ) * rair * tv
          enddo
       enddo
310 end do

    return  

  end subroutine getGeoPotentialFullLevel









  subroutine omint(klon, klat, klev, omf, pretaf, omh, pretah, &
       pres, ompres, af, bf, ah, bh, refpre)
    implicit none  
    !     vertical interpolation of omega(vertical velocity)
    !
    !     input parameters
    !     ----------------
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) full model levels
    !     omf    =  part of omega given at model full levels
    !             (dimension klon,klat*klev)
    !     pretaf =  pressure of the model full levels
    !             (dimension klon,klat*klev)
    !     omh    =  part of omega given at model half levels
    !             (dimension klon,klat*(klev+1))
    !     pretah =  pressure of the model half levels
    !             (dimension klon,klat*(klev+1))
    !     pres   =  pressure of output pressure level
    !     af     =  a-values of model full levels
    !             (dimension klev)
    !     bf     =  b-values of model full levels
    !             (dimension klev)
    !     ah     =  a-values of model half levels
    !             (dimension klev+1)
    !     bh     =  b-values of model half levels
    !             (dimension klev+1)
    !     refpre =  reference pressure of the model vertical
    !               coordinate
    !
    !     output parameters
    !     -----------------
    !     ompres =  pressure level vertical velocity field
    !             (dimension klat*klon)
    integer,intent(in)::klon, klat, klev  
    real(kind=realkind),intent(in)::omf(klon,klat,klev), pretaf(klon,klat,klev), &
         omh(klon,klat,klev+1), pretah(klon,klat,klev+1)
    real(kind=realkind),intent(inout)::ompres(klon,klat)
    real(kind=realkind),intent(in)::af(klev),bf(klev),ah(klev+1),bh(klev+1), &
         pres, refpre
    real(kind=realkind)::etaf(klev),etah(klev+1),aletaf(klev),aletah(klev+ 1), &
         dalefi(klev-1), dalehi(klev), zaleta

    integer :: jx, jy, jk  

    !     1.0 initialization and generation
    !         of  vertical interpolation functions
    do 100 jk = 1, klev  
       etaf(jk) = af(jk) / refpre+bf(jk)  
       aletaf(jk) = log(etaf(jk) )  
100 end do

    do 150 jk = 2, klev  
       dalefi(jk - 1) = 1.0_realkind /(aletaf(jk) - aletaf(jk - 1) )  
150 end do

    do 160 jk = 1, klev + 1  
       etah(jk) = ah(jk) / refpre+bh(jk)  
       aletah(jk) = log(etah(jk) )  
160 end do

    do 170 jk = 2, klev + 1  
       dalehi(jk - 1) = 1.0_realkind /(aletah(jk) - aletah(jk - 1) )  
170 end do

    !*    2.0  omega = 0. above the highest model full level
    !*         or below the lowest model full level
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres<=pretaf(jx,jy,1).or.pres>=pretaf(jx,jy,klev)) then
             ompres(jx, jy) = 0.0_realkind  
          endif
       enddo
    enddo
200 continue  

    !     3.0 interpolation between two model full levels
    do 310 jk = 2, klev  
       do jy = 1, klat  
          do jx = 1, klon  
             if(pres>pretaf(jx,jy,jk-1).and.pres<=pretaf(jx,jy,jk))then
                zaleta = log(etaf(jk - 1) +(pres - pretaf(jx, jy, jk - &
                     1) ) /(pretaf(jx, jy, jk) - pretaf(jx, jy, jk - 1) ) &
                     *(etaf(jk) - etaf(jk - 1) ) )
                ompres(jx, jy) = omf(jx, jy, jk - 1) +(zaleta - aletaf( &
                     jk - 1) ) * dalefi(jk - 1) *(omf(jx, jy, jk) - omf(jx, &
                     jy, jk - 1) )
             endif
          enddo
       enddo
310 end do
    !    4.0  below the lowest half eta level(below ground)
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres>=pretah(jx, jy, klev + 1) ) then  
             ompres(jx, jy) = omh(jx, jy, klev + 1)  
          endif
       enddo
    enddo

    !*    5.0  between two model half levels
    do 510 jk = 2, klev + 1  
       do jy = 1, klat  
          do jx = 1, klon  
             if(pres>pretah(jx,jy,jk-1).and.pres<=pretah(jx,jy,jk))then
                zaleta = log(etah(jk - 1) +(pres - pretah(jx, jy, jk - &
                     1) ) /(pretah(jx, jy, jk) - pretah(jx, jy, jk - 1) ) &
                     *(etah(jk) - etah(jk - 1) ) )
                ompres(jx,jy) = ompres(jx, jy) + omh(jx, jy, jk - 1) &
                     +(zaleta-aletah(jk-1))*dalehi(jk-1)*(omh(jx,jy,jk) &
                     - omh(jx, jy, jk - 1) )
             endif
          enddo
       enddo
510 end do
    return  
  end subroutine omint





  subroutine tpint(lphys, t, pretaf, ps, ts, fis, klon, klat, klev, &
       pres, tpres, af, bf, refpre)
    use confys  
    implicit none  
    !
    !----------------------------------------------------------------------
    !
    !     subroutine tpint
    !     ----------------
    !
    !     purpose
    !     -------
    !     vertical interpolation of temperature to pressure levels
    !
    !     input parameters
    !     ----------------
    !     lphys  =  .t. if physics variables(ts) are available
    !               .f. if physics variables are not available
    !     t      =  temperature at model full levels
    !             (dimension klon,klat*klev)
    !     pretaf =  pressure of the model full levels
    !             (dimension klon,klat*klev)
    !     ps     =  surface pressure(dim. klon,klat)
    !     ts     =  surface temperature(dim. klon,klat)
    !     fis    =  surface geopotential(dim. klon,klat)
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) model levels
    !     pres   =  pressure of the output pressure level
    !     numlon =  number of output grid-points in the x-direction
    !     numlat =  number of output grid-points in the y-direction
    !     af     =  a-values of the model vertical coordinate
    !             (dimension klev)
    !     bf     =  b-values of the model vertical coordinate
    !             (dimension klev)
    !     refpre =  reference pressure of the model vertical
    !               coordinate
    !     ind    =  index field giving references of the gridpoints in
    !               the output field to the gridpoints of the
    !               input field
    !
    !     output parameters
    !     -----------------
    !     tpres  =  pressure level temperature field
    !             (dimension numlat*numlon)
    !
    !----------------------------------------------------------------------
    !
    logical,intent(in) :: lphys  
    integer,intent(in) :: klon, klat, klev  
    real(kind=realkind),intent(in)::pres, refpre  
    real(kind=realkind),intent(in)::t(klon,klat,klev),pretaf(klon,klat,klev)
    real(kind=realkind),intent(in)::ps(klon,klat),ts(klon,klat),fis(klon,klat)
    real(kind=realkind),intent(inout)::tpres(klon, klat)
    real(kind=realkind),intent(in)::af(klev), bf(klev)

    real(kind=realkind) :: etaf(klev), aletaf(klev), dalefi(klev - 1)  
    integer :: jx, jy, jk  
    real(kind=realkind) :: tstar, alsint, fisdvg, alfa, tzero, tplat, tzprim, zaleta
    real(kind=realkind),parameter::tlapse = 0.0065_realkind  

    !  compute vertical interpolation factors
    do  jk = 1, klev  
       etaf(jk) = af(jk) / refpre+bf(jk)  
       aletaf(jk) = log(etaf(jk) )  
 enddo
    !
    do  jk = 2, klev  
       dalefi(jk - 1) = 1.0_realkind /(aletaf(jk) - aletaf(jk - 1) )  
 enddo
    !  extrapolation above the highest model level
    do jy = 1, klat  
       do jx = 1, klon  
          if(pres<=pretaf(jx, jy, 1) ) then  
             tpres(jx, jy) = t(jx, jy, 1)  
          endif
       enddo
    enddo
    !  extrapolation below the lowest model level
    do jy = 1, klat  
       do jx = 1, klon  
          !*        4.1 extract surface temperature
          if(lphys) then  
             tstar = ts(jx, jy)  
          else  
             tstar = t(jx, jy, klev) *(1.0_realkind + tlapse * rair / gravit * &
                  (1.0_realkind / bf(klev) - 1.0_realkind) )
          endif
          if(pres>=pretaf(jx, jy, klev).and.pres<=ps(jx, jy))then
             ! interpolation between ground and lowest model level
             tpres(jx, jy) =((ps(jx, jy) - pres) * t(jx, jy, klev) &
                  +(pres-pretaf(jx,jy,klev))*tstar)/(ps(jx,jy) &
                  - pretaf(jx, jy, klev) )
          endif
          ! extrapolation below the ground
          if(pres>ps(jx, jy) ) then  
             alsint = log(pres / ps(jx, jy) )  
             fisdvg = fis(jx, jy) / gravit  
             alfa = tlapse * rair / gravit  
             if(fisdvg>2000._realkind) then  
                tzero = tstar + tlapse * fisdvg  
                tplat = min(tzero, 298._realkind)  
                if(fisdvg>2500._realkind) then  
                   tzprim = tplat  
                else  
                   tzprim=0.002_realkind*((2500._realkind-fisdvg)*tzero+(fisdvg-2000._realkind)*tplat)
                endif
                alfa = rair *(tzprim - tstar) / fis(jx, jy)  
             endif
             tpres(jx,jy)=tstar*(1.0_realkind+alfa*alsint+0.5_realkind*(alfa*alsint) **2.0_realkind &
                  + 1._realkind / 6._realkind *(alfa * alsint) **3.0_realkind)
          endif
       enddo
    enddo
    !       interpolation between two model levels
    do  jk = 2, klev  
       do jy = 1, klat  
          do jx = 1, klon  
             if(pres>pretaf(jx,jy,jk-1).and.pres<=pretaf(jx,jy,jk)) then
                zaleta = log(etaf(jk - 1) +(pres - pretaf(jx, jy, jk - &
                     1) ) /(pretaf(jx, jy, jk) - pretaf(jx, jy, jk - 1) ) &
                     *(etaf(jk) - etaf(jk - 1) ) )
                tpres(jx, jy) = t(jx, jy, jk - 1) +(zaleta - aletaf(jk - 1) ) &
                     * dalefi(jk - 1) *(t(jx, jy, jk) - t(jx, jy, jk -1) )
             endif
          enddo
       enddo
 enddo
    return  
  end subroutine tpint






  subroutine edcomp(klon, klat, klev, pps, pu, pv, pedot, ahyb, bhyb)
    ! edcomp - compute vertical velocity
    !
    !     j.e. haugen         hirlam
    !     k.s. eerola         hirlam4(revised)
    !     compute etadot at full levels in the equation
    !
    !     df/dt = etadot*(df/deta)
    !
    !     etadot(k) = 0.5*((etadot*(dp/deta)(k-1/2)
    !     +(etadot*(dp/deta)(k+1/2) )*(deta/dp)(k)
    !
    !     where
    !
    !     deta(k) = dak(k)/p0 + dbk
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
    !     ------------------
    !
    !     pedot     vertical velocity variable
    use RCAdomainMod

    implicit none  
    integer,intent(in) :: klon, klat, klev  
    real(kind=realkind),intent(in)::pps(klon, klat), pu(klon, klat, klev), pv(klon,klat,klev)
    real(kind=realkind),intent(out)::pedot(klon, klat, klev)

    real(kind=realkind),intent(in):: ahyb(klev + 1), bhyb(klev + 1)
    !         declaration of local work-space
    integer :: jx, jy, k  
    real(kind=realkind) :: zdek, zrdlah, zrdloh, zdbk, zdak, zrdlo, zrdla, zrear  

    real(kind=realkind),  allocatable, dimension(:, :) ::zahxhy, zdpk, zdpsdt, &
         zdivk, zedpde, zuu, zvv

    allocate(zahxhy(klon,klat), zdpk(klon,klat), zdpsdt(klon,klat),&
         zdivk(klon,klat), zedpde(klon,klat), zuu(klon,klat), &
         zvv(klon,klat) )

    zrdlo = rdlam  
    zrdla = rdth  

    zrear = 1.0_realkind / ra  
    !         computation surface pressure tendency
    do jy = 1, klat  
       do jx = 1, klon  
          zuu(jx, jy) = 0.0_realkind  
          zvv(jx, jy) = 0.0_realkind  
          zahxhy(jx, jy) = rhxu(jx,jy)*rhyv(jx,jy)*ra
       enddo
    enddo

    do k = 1, klev  
       zdak = ahyb(k + 1) - ahyb(k)  
       zdbk = bhyb(k + 1) - bhyb(k)  
       do jy = 1, klat  
          do jx = 1, klon  
             zdpk(jx, jy) = zdak + zdbk * pps(jx, jy)  
          enddo
       enddo
       do jy = 1, klat - 1  
          do jx = 1, klon - 1  
             zuu(jx,jy) = zuu(jx,jy)-pu(jx,jy,k)*(zdpk(jx,jy)+zdpk(jx+1,jy))
             zvv(jx,jy) = zvv(jx,jy)-pv(jx,jy,k)*(zdpk(jx,jy)+zdpk(jx,jy+1))
          enddo
       enddo
    enddo

    zrdloh = rdlam * 0.5_realkind  
    zrdlah = rdth * 0.5_realkind  

    do jy = 2, klat-1  
       do jx = 2, klon-1  
          zdpsdt(jx,jy)=zahxhy(jx,jy)*((zuu(jx,jy)*hyu(jx, jy) &
               - zuu(jx-1,jy)*hyu(jx-1,jy))*zrdloh + &
               (zvv(jx, jy)*hxv(jx,jy)-zvv(jx,jy-1)*hxv(jx,jy-1))*zrdlah)
       enddo
    enddo
    !         initiate vertical loop
    do jy = 1, klat  
       do jx = 1, klon  
          zedpde(jx,jy) = 0.0_realkind  
       enddo
    enddo
    !         vertical integration of the continuity equation

    do 1000 k = klev,1,-1  
       zdak = ahyb(k + 1) - ahyb(k)  
       zdbk = bhyb(k + 1) - bhyb(k)  

       zdek = zdak /sip0 + zdbk  
       do jy = 1, klat  
          do jx = 1, klon  
             zdpk(jx,jy)=zdak+zdbk*pps(jx,jy)  
          enddo

       enddo
       do jy = 1, klat  
          do jx = 1, klon - 1  
             zuu(jx,jy)=0.5_realkind*(zdpk(jx+1,jy)+zdpk(jx,jy)) &
                  * pu(jx,jy,k) * hyu(jx,jy)
          enddo

       enddo
       do jy = 1, klat - 1  
          do jx = 1, klon  
             zvv(jx,jy) = 0.5_realkind*(zdpk(jx,jy+1)+zdpk(jx,jy)) &
                  * pv(jx,jy,k)*hxv(jx,jy)
          enddo
       enddo
       !        compute divergence at level k
       do jy = 2,klat-1  
          do jx = 2,klon-1  
             zdivk(jx,jy)=zahxhy(jx,jy)*((zuu(jx,jy)-zuu(jx-1,jy))*zrdlo &
                  +(zvv(jx,jy)-zvv(jx,jy-1))*zrdla)
          enddo
       enddo
       !     update vertical velocity
       do jy = 2,klat-1  
          do jx = 2,klon-1  
             pedot(jx,jy,k) = zedpde(jx,jy)  
             zedpde(jx,jy)  = zedpde(jx,jy)+zdbk*zdpsdt(jx,jy)+zdivk(jx, jy)
             pedot(jx,jy,k)=0.5_realkind*(pedot(jx,jy,k)+zedpde(jx,jy))*zdek/zdpk(jx,jy)
          enddo
       enddo
       do jx = 1, klon  
          pedot(jx, 1, k) = 0.0_realkind  
          pedot(jx, klat, k) = 0.0_realkind  
       enddo
       do jy = 2, klat - 1  
          pedot(1, jy, k) = 0.0_realkind  
          pedot(klon, jy, k) = 0.0_realkind  
       enddo

1000 enddo
    call swap(pedot,klon,klat,klev)

    deallocate(zahxhy, zdpk, zdpsdt, zdivk, zedpde, zuu, zvv)  
    return  
  end subroutine edcomp





  subroutine intvert2(klon, klat, klev, q, ps, ahyb, bhyb, q2d)  
    !
    ! version 2
    ! called from GEMINI
    ! version one was called from PHYS
    !

    use confys  
    implicit none  
    !
    !  written by Ulf Hansson, Rossby Centre, 9807
    !
    ! vertical integration of water vapor originally
    ! and then cloud water
    ! or anything
    !  q      - specific humidity
    !  q2d - precipitable water vertically integrated
    !  rgravit - 1/gravit
    !
    integer,intent(in)::klon,klat,klev  
    real(kind=realkind),intent(in) :: q(klon, klat, klev), ps(klon, klat)
    real(kind=realkind)::q2d(klon, klat)  

    real(kind=realkind),intent(in)::ahyb(klev + 1), bhyb(klev + 1)  
    real(kind=realkind) :: dpf(klon, klat, klev)  
    real(kind=realkind) :: phalfp1(klon, klat), phalf(klon, klat)  
    real(kind=realkind) :: rgravit  
    integer :: i, j, k

    do j=1,klat  
       do i=1,klon 
          q2d(i,j)=0.0_realkind 
          phalfp1(i,j) = ahyb(klev+1)+ bhyb(klev+1)*ps(i,j) 
          phalf(i,j)    = ahyb(klev)  + bhyb(klev)*ps(i,j)
          dpf(i,j,klev) = phalfp1(i,j)-phalf(i,j)
          do k=klev-1,1,-1
             phalfp1(i,j) = phalf(i,j)
             phalf(i,j)=ahyb(k)+bhyb(k)*ps(i,j)
             dpf(i,j,k)=phalfp1(i,j)-phalf(i,j)
          enddo
       enddo
    enddo
    !

    rgravit=1.0_realkind/gravit
    do k=1,klev
       do j=1,klat
          do i=1,klon
             q2d(i,j)=q2d(i,j)+dpf(i,j,k)*rgravit*q(i,j,k)
          enddo
       enddo
    enddo
    return
  end subroutine intvert2




  subroutine intvert_moist(klon,klat,klev,q,u,v,ps,qu2d,qv2d,vimfc2d,ahyb,bhyb)

    use confys
    use RCAdomainMod
    implicit none
    real(kind=realkind),intent(in)::q(klon,klat,klev),ps(klon,klat)
    real(kind=realkind),intent(in)::u(klon,klat,klev),v(klon,klat,klev)
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind),intent(out)::qu2d(klon,klat),qv2d(klon,klat)

    real(kind=realkind),intent(out)::vimfc2d(klon,klat)
    real(kind=realkind)::zrdlo,zrdla,zrear
    real(kind=realkind),allocatable::zhxhy(:,:),zrhxhy(:,:)
    real(kind=realkind),allocatable::dpf(:,:,:)
    real(kind=realkind),allocatable::phalfp1(:,:),phalf(:,:)
    real(kind=realkind)::rgravit
    integer::klon,klat,klev


    integer::i,j,k


    allocate(zhxhy(klon,klat),zrhxhy(klon,klat),dpf(klon,klat,klev),&
         phalfp1(klon,klat),phalf(klon,klat))

    do j=1,klat
       do i=1,klon
          qu2d(i,j)=0.0_realkind
          qv2d(i,j)=0.0_realkind
          vimfc2d(i,j)=0.0_realkind
          phalfp1(i,j)=ahyb(klev+1)+bhyb(klev+1)*ps(i,j)
          phalf(i,j)=ahyb(klev)+bhyb(klev)*ps(i,j)
          dpf(i,j,klev)=phalfp1(i,j)-phalf(i,j)
          do k=klev-1,1,-1
             phalfp1(i,j)=phalf(i,j)
             phalf(i,j)=ahyb(k)+bhyb(k)*ps(i,j)
             dpf(i,j,k)=phalfp1(i,j)-phalf(i,j)
          enddo
       enddo
    enddo

    rgravit=1.0_realkind/gravit
    do k=1,klev
       do j=1,klat
          do i=1,klon
             qu2d(i,j)=qu2d(i,j)+dpf(i,j,k)*rgravit*q(i,j,k)*u(i,j,k)
             qv2d(i,j)=qv2d(i,j)+dpf(i,j,k)*rgravit*q(i,j,k)*v(i,j,k)
          enddo
       enddo
    enddo

    zrdlo=rdlam
    zrdla=rdth
    zrear=1.0_realkind/ra

    do j=1,klat
       do i=1,klon
          zhxhy(i,j)=zrear/(rhxu(i,j)*rhyv(i,j))
          zrhxhy(i,j)=1.0_realkind/zhxhy(i,j)
       enddo
    enddo

    do j=2,klat-1
       do i=2,klon-1
          vimfc2d(i,j)=-0.5_realkind*zrhxhy(i,j)*((qu2d(i+1,j)* &
               (hyu(i+1,j)+hyu(i,j))*0.5_realkind-qu2d(i-1,j)*&
               (hyu(i-1,j)+hyu(i,j))*0.5_realkind)*zrdlo+(qv2d(i,j+1)&
               *(hxv(i,j+1)+hxv(i,j))*0.5_realkind-qv2d(i,j-1)*&
               (hxv(i,j)+hxv(i,j-1))*0.5_realkind)*zrdla)
       enddo

    enddo
    call swap2d(vimfc2d,klon,klat)
    deallocate(zhxhy,zrhxhy,dpf,phalfp1,phalf)
    return
  end subroutine intvert_moist




  subroutine pslcom(lphys, ps, fis, t, ts, psl, klon, klat, klev,bf)
    use confys  
    implicit none  
    !     computation of mean sea-level pressure
    !
    !     input parameters
    !     ----------------
    !     lphys  =  .t. if physics variables(ts) are available
    !               .f. if physics variables are not available
    !     ps     =  surface pressure(dim. klon,klat)
    !     fis    =  surface geopotential(dim.klon,klat)
    !     t      =  temperature of the model full levels
    !             (dimension klon,klat*klev)
    !     ts     =  surface temperature(dim. klon,klat)
    !     klon   =  number of input grid-points in the x-direction
    !     klat   =  number of input grid-points in the y-direction
    !     klev   =  number of(input) model levels
    !     numlon =  number of output grid-points in the x-direction
    !     numlat =  number of output grid-points in the y-direction
    !     bf     =  b-values of model full levels(dimension klev)
    !     ind    =  index field giving references of the gridpoints in
    !               the output field to the gridpoints of the
    !               input field(dimension numlon*numlat)
    !
    !     output parameters
    !     -----------------
    !     psl    =  output mean sea-level pressure
    !             (dimension numlat*numlon)
    !
    logical,intent(in) :: lphys  
    integer,intent(in) :: klon, klat, klev  

    real(kind=realkind),intent(in)::ps(klon,klat),fis(klon,klat),t(klon,klat,klev),ts(klon,klat), bf(klev)
    real(kind=realkind),intent(out)::psl(klon, klat)
    real(kind=realkind) ::  tstar, alfa, tzero, arg  
    integer :: jx, jy  
    real(kind=realkind),parameter:: tlapse = 0.0065_realkind  

    do jy = 1, klat  
       do jx = 1, klon  
          !        surface temperature and temp.reduction factor
          if(lphys) then  
             tstar = ts(jx, jy)  
          else  
             tstar = t(jx,jy,klev)*(1.0_realkind+tlapse*rair/gravit*(1.0_realkind/bf(klev)-1.0_realkind))
          endif

          if(tstar<255.0_realkind) then  
             !         prevent too high pressures
             !            under cold terrain
             tstar = 0.5_realkind *(tstar + 255._realkind)  
             alfa = tlapse*rair/gravit  
          else  
             !         trial extra-polation of sea-level
             !             temperature
             tzero = tstar + tlapse * fis(jx, jy) / gravit  
             if(tzero>290.5_realkind) then  
                !        prevent too low pressures under
                !            hot terrain
                if(tstar<=290.5_realkind) then  
                   alfa = rair *(290.5_realkind - tstar) / fis(jx, jy)  
                else  
                   alfa = 0._realkind  
                   tstar = 0.5_realkind *(290.5_realkind + tstar)  
                endif
             else  
                !        normal temp.reduction factor
                alfa = tlapse * rair / gravit  
             endif
          endif
          !        compute sea-level pressure
          !            end loop over gridpoints
          arg = fis(jx, jy) / rair / tstar  
          psl(jx,jy) = ps(jx,jy)*exp(arg*(1.0_realkind-0.5_realkind*alfa*arg+1.0_realkind/3.0_realkind* &
               (alfa*arg)**2.0_realkind) )
       enddo
    enddo
    return  
  end subroutine pslcom





  subroutine find2points( klon, klat, alnglob, cltglob, pslglob,  slpdebilt, &
       slpoksoy, lphys, ps, fis, t, ts, klev, bhyb)
    !
    ! written by ulf hansson, rossby centre - 010904
    ! changed by ulf hansson, rossby centre - 030103
    !                         use  colfld instead of collect2global(using co
    !
    ! find coordinates close to debilt and oksoy
    !
    ! extract sea level pressure in these points
    !
    ! output
    !   slpdebilt
    !   slpoksoy
    !
    ! input
    !   the rest

    use RCAdomainMod

    implicit none  
    integer,intent(in)::klon, klat,klev  
    real(kind=realkind),intent(inout)::alnglob(klon_global, klat_global)  
    real(kind=realkind),intent(inout)::cltglob(klon_global, klat_global)  
    real(kind=realkind),intent(inout)::pslglob(klon_global, klat_global)  
    real(kind=realkind),intent(inout)::slpdebilt, slpoksoy  
    real(kind=realkind),intent(in)::ps(klon, klat), fis(klon, klat), t(klon, klat, klev)  
    real(kind=realkind),intent(in)::ts(klon, klat)  
    real(kind=realkind),intent(in)::bhyb(klev + 1)  

    logical,intent(in)::lphys  
    real(kind=realkind)::bfull(klev)  
    real(kind=realkind) :: psl(klon, klat)  
    real(kind=realkind):: pi, rad2deg  
    real(kind=realkind):: lonI, latI  !override from a use?
    real(kind=realkind):: lon1, lat1  
    real(kind=realkind):: lon2, lat2  
    real(kind=realkind):: dist1, dist2  
    real(kind=realkind):: smalldist1, smalldist2  
    integer :: i, j, i1, j1, i2, j2  
    logical,save ::  lfirst=.true.  

    save i1, j1, i2, j2  
    pi = 2._realkind * asin(1.0_realkind)  

    rad2deg = 180._realkind / pi  
    do j = 1, klev  
       bfull(j) = 0.5_realkind *(bhyb(j) + bhyb(j + 1) )  

    enddo

    call pslcom(lphys, ps, fis, t, ts, psl, klon, klat, klev, bfull)  
    call colfld(0, alnglob, plong, klon, klat)  
    call colfld(0, cltglob, pclat, klon, klat)  

    call colfld(0, pslglob, psl, klon, klat)  

    if(mype/=0) return  

    !
    ! find coordinates close to debilt and oksoy
    ! do it once only and save them
    !

    if(lfirst) then  
       ! 52 deg  6 min ! debilt
       lat1 = 52.1_realkind  
       !  5 deg 11 min ! debilt
       lon1 = 5.1833_realkind  
       ! 58 deg  4 min ! oksoy
       lat2 = 58.0667_realkind  
       !  8 deg  3 min ! oksoy

       lon2 = 8.05_realkind  
       smalldist1 = 999._realkind  
       smalldist2 = 999._realkind  
       i1 = - 9  
       j1 = - 9  
       i2 = - 9  

       j2 = - 9  
       do i = 1, klon_global  
          do j = 1, klat_global  
             cltglob(i, j) = rad2deg * acos(cltglob(i, j) )  
             lonI = alnglob(i, j)  
             latI = cltglob(i, j)  
             ! to avoid crossing greenwichproblem
             if(lonI>180._realkind) lonI = lonI - 360._realkind  
             dist1 = sqrt((latI - lat1) **2.0_realkind +(lonI - lon1) **2.0_realkind)  
             dist2 = sqrt((latI - lat2) **2.0_realkind +(lonI - lon2) **2.0_realkind)  
             if(dist1<smalldist1) then  
                smalldist1 = dist1  
                i1 = i  
                j1 = j  
             endif
             if(dist2<smalldist2) then  
                smalldist2 = dist2  
                i2 = i  
                j2 = j  
             endif
          enddo
       enddo

       lfirst = .false.  
    endif
    slpdebilt = pslglob(i1, j1)  
    slpoksoy = pslglob(i2, j2)  
    return  
  end subroutine find2points




  subroutine gwwrloc(mype_out,lun,type,par,plev,f,klon,klat,myddr)
    !  purpose:
    !      write a field to an asimof file
    !  description:
    !    collects a field on one pe and writes it into an asimof file
    !  author:
    !      kalle eerola   hirlam    1998


    use modddr
    implicit none  
    integer,intent(in)::mype_out,plev, lun, type, par, klon, klat
    type(ddr)::myddr

    real(kind=realkind),intent(in) ::  f(klon, klat)  
    integer :: jx, jy 

    real(kind=realkind) , allocatable, dimension(:, :) ::f_global
    real(kind=realgribkind) , allocatable, dimension(:, :) ::f_out  

    allocate(f_global(klon_global,klat_global),f_out(npplon,npplat))

    call colfld(mype_out, f_global, f, klon, klat)  

    if(mype==mype_out) then  
       do jy=1,npplat  
          do jx=1,npplon  
             f_out(jx,jy) = f_global(jx+iminpp-1,jy+jminpp-1)  
          enddo

       enddo
       call gwrite(lun, type, par, plev, f_out, npplon, npplat,myddr)
    endif
    deallocate(f_global,f_out)  
    return  

  end subroutine gwwrloc


  subroutine gwopen2(ludir,&
       year,month,day,hour,minut,length,prefix,sufix,myddr)
    use asimof
    use modddr
    implicit none
    !     Asimof parameters
    !     Input
    !     N.B. The old parameter IBUF is here used as a real containing
    !     the vertical parameters A and B.

    integer,intent(in)::ludir 
    real(kind=realgribkind),allocatable::pab(:,:) 
    integer,intent(in)::year,month,day,hour,minut,length
    character(len=*),intent(in)::prefix,sufix
    integer iyear
    integer,parameter::klb1=30,nstr=75,knfldx=8000
    integer jprex,jsufx
    integer knfld,knlev,kerr,knlevx
    integer kbks1(klb1,knfldx)
    logical lapab
    character(len=nstr)::filename
    integer k,jstr,ilen
    type(ddr)::myddr

    !     set file name
    iyear=year
    jprex = index(prefix,' ') - 1
    if(jprex<=0)jprex=len(prefix)
    jsufx = index(sufix,' ')  - 1
    if(jsufx<=0)jsufx=len(sufix)
    jstr=jprex+4+4*2+jsufx

    if(jstr>nstr)then
       print *,__FILE__,__LINE__
       write(6,*)'error in gwopen2:  string too long, length=',jstr
       call stop_program('error in gwopen2:  string too long')
    endif

    ilen=nint(real(length)/3600.)
    write(filename,900)prefix(1:jprex),iyear,month,day,hour,minut,sufix(1:jsufx)
900 format(a,i4.4,4i2.2,a)


    lapab = myddr%nlevhl > 0
    knlevx = myddr%nlevhl
    knlev = myddr%nlevhl

    allocate(pab(2,myddr%nlevhl))
    !     Set the A:s and B:s in PAB
    do 200 k=1,myddr%nlevhl
       pab(1,k)=myddr%alevhl(k,1)
       pab(2,k)=myddr%alevhl(k,2)
200 enddo

    !     load file
    call asimho(ludir,filename(1:jstr),-1,kerr)
    if(kerr /= 0) then
       print *,__FILE__,__LINE__
       write(6,*)'gwopen2 error opening ',filename(1:jstr),' code:',kerr
       call stop_program('gwopen2 error opening ')
    endif
    call loadfd(ludir,filename(1:jstr),kbks1,klb1,knfldx,knfld, &
         lapab,pab,knlevx,knlev,-890,kerr)


    if(kerr/=0)then
       print *,__FILE__,__LINE__
       write(6,*)' error in loadfd: kerr=',kerr
       write(6,*)' called by gwopen'
       call stop_program(' error in loadfd called by gwopen')
    endif
    deallocate(pab)
  end subroutine gwopen2

  subroutine gwclos2(ludir)
    use asimof
    implicit none

    integer,intent(in):: ludir
    integer ierr
    call asimhc(ludir,ierr)
  end subroutine gwclos2

  subroutine gwrite(ludir,ktyp,kwmo,alev,field,mx,my,myddr)

    use asimof
    use modddr
    use mod_grib
    implicit none
    integer,parameter::klb1=30
    integer,parameter::knfldx=8000
    integer,parameter::klevx=100
    integer,parameter::klb2=400
    !     old grid input parameters

    integer,intent(in)::alev,ludir,ktyp,kwmo,mx,my
    real(kind=realgribkind),intent(in):: field(mx,my)
    type(ddr)::myddr
    !     interface to asimof
    !     asimof parameters

    integer ibpw,ibpwio

    integer kdev,fieldSize,kerr
    integer kb1(klb1),kb2(klb2)
    real(kind=realgribkind)    pb2(klb2),plev,pacc
    real(kind=realgribkind)    amax,amin

    integer  jtyp,jpar,jlev
    data     ibpw /0/, ibpwio /0/

    call asimhm(ibpw,ibpwio,notdef)
    kdev=ludir
    fieldSize=mx*my

    jpar=kwmo
    jlev=alev
    jtyp=ktyp
    if(ktyp==109.or.ktyp==0) then
       jtyp=109
       plev=myddr%alevhl(jlev,1)+myddr%alevhl(jlev,2)*sip0  !rlevhl(1)
    else if(ktyp==100)then
       plev=real(alev,realgribkind)*100.0_realgribkind
    else
       plev= sip0 !rlevhl(1)
    endif
    pacc=0.0_realgribkind !accu()

    call maxmin(field,fieldSize,amax,amin,real(notdef,realgribkind))
#ifdef DEBUG
    write(6,900) 'writing: unit,par,typ,alev,max,min', &
         kdev,jpar,jtyp,real(alev,realgribkind),amax,amin
900 format(a,3i4,f12.2,2e12.4)
#endif
    !     seting parameters in block 1 and 2

    call block1(jtyp,jpar,real(alev,realgribkind),kb1,myddr)

    ! parameter table version number 138 if jpar>999
    if(jpar>=1000)then
       jpar=jpar-1000
       kb1(4)=138
    endif

    call block2(jtyp,jpar,real(alev,realgribkind),kb2,pb2,myddr)
    call putfd(kdev,kb1,klb1,kb2,klb2,pb2,klb2,field,fieldSize,pacc,kerr)
    if(kerr/=0)then
       print *,__FILE__,__LINE__
       write(6,*)' write error in putfd: unit=',ludir,' kerr=',kerr
       write(6,*)' called by gwrite'
       call stop_program(' write error in putfd:called by gwrite')
    endif
  end subroutine gwrite




  subroutine block1(ktyp,kpar,alev,kdef,myddr)

    !     subroutine block1 packs the field information of block 1.
    !     identification data from the common area comddr plus some
    !     extra parameters are packed into 24 octets(i.e. lenb=8 )
    !     
    !     input :
    !     ktyp : level type
    !     kpar : parameter
    !     alev : level of parameter
    !     
    !     output :
    !     kdef  : integer array of grib - code - block 1
    !     
    use modddr
    implicit none

    type(ddr)::myddr
    integer kpar,iflag,ktyp,lev,ilenfc,issrest
    real(kind=realgribkind) alev
    integer kdef(*)
    logical loaccu            ! for accumulated parameters

    !     accumulated parameters: precip, net surface rad, surface fluxes
    loaccu=                                                    &
         kpar == 57 .or.                                      &
         kpar == 61 .or. kpar == 62 .or.                    &
         kpar == 63 .or. kpar == 65 .or.                    &
         kpar == 78 .or. kpar == 79 .or.                    &
         kpar ==111 .or. kpar ==112 .or.                    &
         kpar ==113 .or. kpar ==114 .or.                    &
         kpar ==115 .or. kpar ==117 .or.                    &
         kpar ==121 .or. kpar ==122 .or.                    &
         kpar ==123 .or. kpar ==124 .or.                    &
         kpar ==125 .or.                                      &
         kpar ==128 .or.                                      &
         kpar ==210 .or. kpar ==211 .or. kpar == 212 .or.  &
         kpar ==220 .or. kpar ==221 .or.                    &
         kpar ==230 .or. kpar ==231 

    iflag = 128

    kdef(1) = 24
    kdef(2) = 0
    kdef(3) = 0
    kdef(4) = 1
    kdef(5) = 96
    kdef(6) = 1
    kdef(7) = 255
    kdef(8) = iflag

    if(ktyp==100)then
       lev=nint(alev*0.01_realgribkind)
    else
       lev=nint(alev)
    endif

    kdef( 9) = kpar
    kdef(10) = ktyp
    kdef(11) = lev
    kdef(12) = 0
    kdef(13) = myddr%ndtbhl / 10000
    kdef(14) = myddr%ndtbhl / 100 - kdef(13) * 100
    kdef(15) = myddr%ndtbhl       - kdef(13) * 10000 - kdef(14) * 100
    kdef(16) = myddr%nscbhl / 10000
    kdef(17) = myddr%nscbhl / 100 - kdef(16) * 100

    !     forecast length in minutes:
    ilenfc = myddr%nflshl/60
    !     can forecast length be expressed in hours?
    issrest  = myddr%nflshl -(myddr%nflshl/3600)*3600

    !     accumulated parameters(in hours, but short non-integer fcsts in min)
    if(loaccu)then
       if(issrest==0.or.ilenfc>255)then
          kdef(18) = 1        ! forecast length in hours
          kdef(20) = nint(real(myddr%nflshl,realgribkind)/3600._realgribkind)
       else
          kdef(18)= 0         ! forecast length in minutes
          kdef(20)=ilenfc
       endif
       kdef(19)=0             ! accumulated from start of forecast
       kdef(21)=4             ! accumulated quantity
    else
       if( issrest == 0) then
          !     integer hour forecast length: forecast length in hours
          kdef(18) = 1        ! forecast length in hours
          kdef(19) = myddr%nflshl/3600
          kdef(20) = 0
!!$          if( nflshl == 0 .and. ndorhl == 1 ) then
!!$             kdef(21) = 1     ! initialised analysis
!!$          else
          kdef(21) = 0     ! forecast or un-initialised analysis
!!$          endif
       else


          !     not in integer hours: forecast length in minutes
          kdef(18) = 0        ! forecast length in minutes
          kdef(19) = ilenfc
          kdef(20) = 0
          kdef(21) = 10       ! forecast length occupies two octets
       endif
    endif

    kdef(22) = 0
    kdef(23) = 0
    kdef(24) = 0

    return
  end subroutine block1

  subroutine block2(ktyp,kpar,alev,kdes,prwb2,myddr)
    !
    !  subroutine block2 packs the field information of block 2.
    !   identification data from the common area comddr plus some
    !   extra parameters are packed(i.e. lenb=8 ).
    !
    !  input :
    !    ktyp : level type
    !    kpar : parameter according to wmo - parameter table
    !    alev : level of parameter
    !
    !  output :
    !    kdes : integer array of grib - code - block 2
    !    prwb2 : real array of grib - code - block 2
    !
    use modddr
    implicit none


    type(ddr)::myddr
    integer kdes(*)
    real(kind=realgribkind) prwb2(*)
    integer lendes,nuusb,j,ktyp,kpar,iscm,klev
    real(kind=realgribkind) alev,vcp

    !   nuusb is the number of unused bits at the end of block 2.
    !   nuusb will always be zero, since lendes is always even.

    lendes = 42 + myddr%npplhl * 4
    nuusb = 0
    do 100 j=1,74
       kdes(j) = 0
100 enddo

    kdes(1) = lendes
    kdes(4) = nuusb
    kdes(5) = 0
    kdes(6) = 10
    kdes(7) = myddr%nlonhl
    kdes(9) = myddr%nlathl
    kdes(11) = nint( myddr%alafhl * 1000._realgribkind )
    kdes(14) = nint( myddr%aweshl * 1000._realgribkind )
    kdes(17) = 136
    kdes(18) = nint( myddr%alalhl * 1000._realgribkind )
    kdes(21) = nint( myddr%aeashl * 1000._realgribkind )
    kdes(24) = nint( abs( myddr%dlonhl * 1000._realgribkind ) )
    kdes(26) = nint( abs( myddr%dlathl * 1000._realgribkind ) )

    !   staggering
    if(ktyp/=105)then
       if(kpar==33)then
          kdes(14) = nint((myddr%aweshl+0.5_realgribkind*myddr%dlonhl) * 1000._realgribkind )
          kdes(21) = nint((myddr%aeashl+0.5_realgribkind*myddr%dlonhl) * 1000._realgribkind )
       endif
       if(kpar==34)then
          kdes(11) = nint((myddr%alafhl+0.5_realgribkind*myddr%dlathl) * 1000._realgribkind )
          kdes(18) = nint((myddr%alalhl+0.5_realgribkind*myddr%dlathl) * 1000._realgribkind )
       endif
    endif

    if( myddr%dlathl < 0._realgribkind ) then
       iscm = 0
    else
       iscm = 64
    endif

    kdes(28) = iscm

    if( kdes(6) == 10 ) then
       !     kdes(33) and kdes(36) contain coordinates of the south pole
       kdes(33) = nint( myddr%aplahl * 1000._realgribkind )
       kdes(36) = nint( myddr%aplohl * 1000._realgribkind )
    else if( kdes(6) == 60 ) then
       !     kdes(33) and kdes(36) contain coordinates of the north pole
       !     with the longitude rotated over 180 degrees
       kdes(33) = nint( myddr%aplahl * 1000._realgribkind )
       kdes(36) = nint( myddr%aplohl * 1000._realgribkind )
    else
       print *,__FILE__,__LINE__
       write(6,*) 'kdes(6) = ' , kdes(6)
       write(6,*) ' !!! think about definition kdes(33,36) !!! '
       call stop_program( ' stop in block2 ')
    end if

    kdes(39) = 1
    kdes(40) = 1

    !    vertical coordinate parameter(s) come(s) next and last
    prwb2(1) = 0._realgribkind
    klev=nint(alev)
    do 200 j=1,myddr%npplhl
       if( ktyp == 109 .or. ktyp == 0) then
          vcp = myddr%alevhl(klev,j)
       else
          vcp = alev
       endif
       kdes(42+j*2-1) = 1 + j
       kdes(42+j*2  ) = 0
       prwb2(1+j) = vcp
200 enddo
    return
  end subroutine block2



end module postprocess
