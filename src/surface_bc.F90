module surface_bc
  use timetype
  use calendar
  use RCAdomainMod
  use domainmod
  use gcm
  use decomp
  use mod_grib,only:realgribkind
  use config,only:use_oasis
  use coupling
  implicit none
  private
  real(kind=realkind),allocatable,dimension(:,:,:)::tsea_bound, frice_bound,tice_bound
  type(time),save::boundaryTime(2)
  type(deltatime),save::dtBoundary
  
  public initSurface_bc,setSurface_bc
contains


  subroutine readoceanpath(ocean_path,ilen)
     use decomp
     implicit none
     character(len=132)::ocean_path
     integer ilen
     logical,save::done=.false.
     character(len=20)namelistfile,inst
     namelist/institute/inst
     namelist /ocean/ocean_path
     ilen=1
     open(66,file='namelists.dat',status='old',form='formatted')
     read(66,nml=institute)
     close(66)
     namelistfile='gcmpaths.'//trim(inst)
     open(61,file=namelistfile,status='old',form='formatted')
     read(61,nml=ocean)
     close(61)
     if(mype.eq.0.and..not.done)then
        write(6,nml=ocean)
        done=.true.
     endif
     ilen=index(ocean_path,' ')-1
   end subroutine readoceanpath



   subroutine setSurface_bc(klon,klat,ktime,tsea,frice,tice,eco,meco,&
        initialization,RCAdom,lcounttice)
     use decomp, only:mype
     implicit none
     integer,intent(in)::klon,klat,meco
     real(kind=realkind),intent(in)::eco(klon,klat,meco)
     real(kind=realkind),intent(in)::lcounttice(klon,klat)
     type(domain),intent(in)::RCAdom
     type(time),intent(in)::ktime
     logical,intent(in)::initialization
     real(kind=realkind),dimension(klon,klat),intent(inout)::tsea,frice,tice
     real(kind=realkind)::timrat
     integer::i,j

     do while(ktime>boundaryTime(2))
        tsea_bound(:,:,1)  = tsea_bound(:,:,2) 
        frice_bound(:,:,1) = frice_bound(:,:,2)
        tice_bound(:,:,1)  = tice_bound(:,:,2) 
        boundaryTime(1) = boundaryTime(2)
        boundaryTime(2) = boundaryTime(2) + dtBoundary

        call getSurface_bc(2,klon,klat, &
             boundaryTime(2),eco(:,:,10),&
             tsea_bound(:,:,2),frice_bound(:,:,2),tice_bound(:,:,2),&
             RCAdom)
     enddo
     if(ktime>=boundaryTime(1).and.ktime<=boundaryTime(2))then
        timrat = real(Nseconds(ktime-boundaryTime(1)),realkind)/real(Nseconds(dtBoundary),realkind)
        if(timrat<0._realkind.or.timrat>1._realkind)then
           print *,'timrat=',timrat
           print *,'ktime=',ktime
           print *,'btime(1)',boundaryTime(1)
           print *,'ktime-btime(1)',ktime-boundaryTime(1)
           print *,Nseconds(ktime-boundaryTime(1)),'seconds'
           print *,Nseconds(dtBoundary),'seconds for dtBoundary'
           call stop_program( 'Error in surface_bc timrat errr')
        endif

        tsea = tsea_bound(:,:,1) + timrat*(tsea_bound(:,:,2)-tsea_bound(:,:,1))
        frice = frice_bound(:,:,1) + timrat*(frice_bound(:,:,2)-frice_bound(:,:,1))

!wangsy
!        if(use_oasis)then
!           where(lcounttice>0.5_realkind) 
!              tsea = sst_fieldm
!              frice = icc_fieldm
!           end where
!        endif
#ifdef DEBUG
        do j=1,klat
           do i=1,klon
              if(frice(i,j)<0_realkind .or.frice(i,j)>1_realkind)then
write(6,*)' i j frice_bound 1 2 frice ',i,j,frice_bound(i,j,1),frice_bound(i,j,2),frice(i,j)
                 call stop_program('frice is not correct')
              endif
           enddo
        enddo
#endif       
        if(initialization)then
           if(lecmwf.or.lhadley)then !we read tice initially from ECMWF also for HADLEY
              tice = tice_bound(:,:,1) + timrat*(tice_bound(:,:,2)-tice_bound(:,:,1))
           endif
        endif
!wangsy
!        if(use_oasis)then
!           where(lcounttice>0.5_realkind) 
!              tice =  ict_fieldm
!           end where
!        endif
     else
        call stop_program( 'Time interpolation error in SFC_BC')
     endif
   end subroutine setSurface_bc


  subroutine initSurface_bc(klon,klat,startTime,eco,meco,RCAdom)

    implicit none
    integer,intent(in)::klon,klat,meco
    type(time),intent(in)::startTime
    type(domain)::RCAdom
    real(kind=realkind),intent(in)::eco(klon,klat,meco)
    integer::dtsec(4)
    integer::k,secmin,kmin

    allocate(tsea_bound(klon,klat,2), frice_bound(klon,klat,2))
    allocate(tice_bound(klon,klat,2))
    
    dtBoundary = makedeltatime(0,0,0,6,0,0)
    boundaryTime(1) = startTime 
    if(boundaryTime(1)%hour/= 0 .and.boundaryTime(1)%hour/= 6 .and. &
         boundaryTime(1)%hour/= 12 .and.  boundaryTime(1)%hour/= 18)then
       do k=1,4
          dtsec(k) = startTime%hour - (k-1)*6
       enddo
       !We need to pick out the positive negative (or zero) dtsec, since the time interpolation assumes data(1) to be 
       !at startTime or before
       secmin = -99999
       kmin = -1
       do k=1,4
          if(dtsec(k)>=0.and.secmin<dtsec(k))then
             secmin = dtsec(k)
             kmin = k
          endif
       enddo
       if(kmin==-1)call stop_program('Error init surface_bc kmin==-1')
       boundaryTime(1) = maketime(startTime%year,startTime%month,startTime%day,(kmin-1)*6,0,0)

    endif

    boundaryTime(2) = boundaryTime(1) + dtBoundary
    
    call getSurface_bc(1,klon,klat,        &
         boundaryTime(1),eco(:,:,10),&
         tsea_bound(:,:,1),frice_bound(:,:,1),tice_bound(:,:,1),&
         RCAdom)
    
    call getSurface_bc(2,klon,klat, &
         boundaryTime(2),eco(:,:,10),&
         tsea_bound(:,:,2),frice_bound(:,:,2),tice_bound(:,:,2),&
         RCAdom)
     
  end subroutine initSurface_bc


  subroutine getSurface_bc(nbdnum,klon,klat,&
       ktime,frl, &
       tsea_bound,frice_bound,tice_bound,&
       RCAdom)

    use decomp

    use modddr
    use calendar

    implicit none
    integer,intent(in)::nbdnum
    integer,intent(in):: klon, klat
    type(domain),intent(in)::RCAdom
    type(time),intent(in)::ktime
    real(kind=realkind),intent(in)::frl(klon,klat)
    real(kind=realkind),intent(inout)::tsea_bound(klon,klat)
    real(kind=realkind),intent(inout)::frice_bound(klon,klat)
    real(kind=realkind),intent(inout)::tice_bound(klon,klat)
    type(time)::filetime
    type(ddr)::gcmddr

    call pre_getgrblow(nbdnum,klon,klat, 16,   frl,   &
          tsea_bound, tice_bound,frice_bound, &
          ktime, & 
          RCAdom, &
          gcmddr)



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
       call stop_program( 'Wrong date in demanding data getdat')
    endif

    return
  end subroutine getSurface_bc



  subroutine pre_getgrblow(nbdnum,klon,klat, lundir,frland, &
       tsea_bound, tice_bound, frice_bound, &
       ktime,&
       RCAdom,&
       gcmddr)       
    use gcm
    use decomp
    use IO
    use climate
    use modddr
    use calendar

    implicit none
    type(domain),intent(in)::RCAdom
    type(ddr),intent(out)::gcmddr
    integer,intent(in):: nbdnum
    integer,intent(in):: klon,klat
    real(kind=realkind),intent(in):: frland(klon,klat)
    real(kind=realkind),intent(inout)::frice_bound(klon,klat),tsea_bound(klon,klat)
    real(kind=realkind),intent(inout)::tice_bound(klon,klat)
    integer,intent(in)::lundir
    type(time),intent(in)::ktime
    integer,save::lunlsm=64
    logical,save::lfirstlsm=.true.

    real(kind=realkind)::north,east

    logical::larou

    logical,save::ltest=.false.

    logical:: ltsea,lfrice
    integer:: jlev
    real(kind=realkind):: stagu_x_grib,stagu_y_grib,stagv_x_grib,stagv_y_grib
    integer:: nbuf
    integer, parameter::jpbuf=700000
    real(kind=realgribkind):: buf(jpbuf)

    integer::nx_grib,ny_grib
    real(kind=realkind)::west_grib,south_grib
    real(kind=realkind)::dlat_grib,dlon_grib
    real(kind=realkind)::polon_grib,polat_grib
    character(len=1)::arakawa_grib

    !in case of hadley there is a grid with other resolution for sst 
    integer:: nx_grib2,ny_grib2
    real(kind=realkind):: dlat_grib2,dlon_grib2



    real(kind=realkind):: zeps
    integer:: jx,jy
    integer:: ilen,istat,iy,im,id,ih,imin,length
    character(len=132):: prefix,filenam,filenam2
    integer:: year1,year2,month1,month2,day2
    real(kind=realkind)::    weigh1,weigh2
    real(kind=realkind)::    tsea2(klon,klat),frice2(klon,klat)
    integer:: lundir2,ilen1,ilen2
    character(len=132):: prefix1,prefix2
    character(len=132):: prefix3,filenam3

    real west_grib_u,south_grib_u,west_grib_v,south_grib_v
    integer::  hour0
    Logical:: Lfis
    real:: west_model, south_model, dlon_model, dlat_model, &
         polon_model, polat_model


    zeps = 1.e-6_realkind

    !     Set indicator for roughness over sea

    larou  = .true.



    !     If required, retrieve surface data in actual model grid

    ltsea   = .true.
    lfrice  = .true.



    if(ltest .and. mype==0.and.lhadley) then
       write(6,*) ' south,west,dlon,dlat,polon,polat = ', &
            RCAdom%south,RCAdom%west,RCAdom%dlon,RCAdom%dlat,RCAdom%polon,RCAdom%polat
    endif


    north = RCAdom%south + real(klat_global-1,realkind)*RCAdom%dlat
    east  = RCAdom%west  + real(klon_global-1,realkind)*RCAdom%dlon

    if( ltest .and. mype==0 ) then
       write(6,*)
       write(6,*)'Output model geometry to be used for'
       write(6,*)'horizontal and vertical interpolations:'
       write(6,*)'south,north=',RCAdom%south,north
       write(6,*)'west,east=',RCAdom%west,east
       write(6,*)'dlon,dlat=',RCAdom%dlon,RCAdom%dlat
       write(6,*)'polon,polat=',RCAdom%polon,RCAdom%polat
       write(6,*)'nlon,nlat=',klon_global,klat_global
       write(6,*)
    endif


    if(lhadley)then

! ***   NOTE  orog at 3.75,2.50 but at global grid

      if(nbdnum.le.2 ) then 
!         nx_grib2     = 97
         nx_grib2     = 96
         ny_grib2     = 73
         west_grib    = 0.
         south_grib   = -90.
         dlon_grib2   = 3.75
         dlat_grib2   = 2.50
         polon_grib   = 0.
         polat_grib   = -90.
         if(ltest .and. mype.eq.0) write(6,*) ' orog grid ', &
          nx_grib2,ny_grib2,dlon_grib2,dlat_grib2,lecmwf

         lfis=.true.

!         call interpol_fis(klon,klat,           &
!                      nbdnum,                   &
!                      lundir,                   &
!                      ktime%year,ktime%month,ktime%day,ktime%hour,      &
!                      RCAdom,    &
!                      lfis,                     &
!                      nx_grib2,ny_grib2,                   &
!                      west_grib,south_grib,dlat_grib2,dlon_grib2,  &
!                      polon_grib,polat_grib,                       &
!cau1106                      lecmwf,gcmddr )
!                      gcmddr )
!                      lprinteta,iprneta,jprneta,lunprneta)

      endif      !  nbdnum
    endif

    if(lhadley)then
       !     NOTE  sst, frice  at finer global grid 
       if(nbdnum<=4 .or. ktime%hour==0 ) then 
          hour0=0
!          nx_grib2     = 193
!          ny_grib2     = 145
          nx_grib2     = 96
          ny_grib2     = 73
          west_grib    = 0._realkind
          south_grib   = -90._realkind
          dlon_grib2   = 1.875_realkind
          dlat_grib2   = 1.25_realkind
          polon_grib   = 0._realkind
          polat_grib   = -90._realkind
          if(ltest .and. mype==0) write(6,*) ' finer grid ', &
               nx_grib2,ny_grib2,dlon_grib2,dlat_grib2,lecmwf

          ltsea=.true.
          lfrice=.true.

          year1 = ktime%year
          year2 = ktime%year
          day2  = 16
          if( ktime%day >= 16 ) then
             month1 = ktime%month
             month2 = ktime%month+1
             if(month2==13) then
                month2=1
                year2 =ktime%year+1
             endif
          else
             month1 = ktime%month-1
             month2 = ktime%month
             if (month1==0) then 
                month1=12
                year1 = ktime%year-1
             endif
          endif

          call interpol_sea(klon,klat,              &
               nbdnum,                                    &
               lundir,                                    &
!cau1106               year1,month1,day2,ktime%hour,                    &
               year1,month1,day2,hour0,                    &
               RCAdom,&
               ltsea,tsea_bound,                                &
               lfrice,frice_bound,                          &
               nx_grib2,ny_grib2,                         &
               west_grib,south_grib,dlat_grib2,dlon_grib2,&
               polon_grib,polat_grib,                     &
               gcmddr)

          call interpol_sea(klon,klat,                &
               nbdnum,                                      &
               lundir,                                      &
!cau1106               year2,month2,day2,ktime%hour,                      &
               year2,month2,day2,hour0,                      &
               RCAdom,&
               ltsea,tsea2,                                 &
               lfrice,frice2,                           &
               nx_grib2,ny_grib2,                           &
               west_grib,south_grib,dlat_grib2,dlon_grib2,  &
               polon_grib,polat_grib,                       &
               gcmddr)

          call int_clim(ktime%day,month1,month2,weigh1,weigh2)
          
          tsea_bound=weigh1*tsea_bound+weigh2*tsea2
          frice_bound=weigh1*frice_bound+weigh2*frice2
          tice_bound=tsea_bound
             
       endif
    endif


    !     Open input data file to find out geometry - needed for dynamical
    !     memory allocation and interpolations

    nbuf = jpbuf

    !get path to global data
    call readboundpath(prefix1,ilen1)
    prefix = prefix1
    ilen=ilen1



    if(lecmwf)then
       lundir2=63
       call readecicepath(prefix2,ilen2)
       prefix2 = prefix2(1:ilen2)//'tice'
       ilen2=ilen2+4
       iy=ktime%year
       im=ktime%month
       id=ktime%day
       ih=00
       imin=0
       length=0
       call cre_filnam(prefix2,ilen2,iy,im,id,ih,length,filenam2)
       call groploc(lundir2,filenam2,buf,nbuf,gcmddr)
       nx_grib = gcmddr%nlonhl
       ny_grib = gcmddr%nlathl
    endif

    if(lecham5)then
       if (lfirstlsm) then
          lunlsm=64
          if (mype==0 ) then
             prefix3 = prefix1(1:ilen1)
             filenam3 = prefix3
          endif
          call groploc(lunlsm,filenam3,buf,nbuf,gcmddr)
          nx_grib = gcmddr%nlonhl
          ny_grib = gcmddr%nlathl

          lfirstlsm=.false.
       endif
    endif

    iy=ktime%year
    im=ktime%month
    id=ktime%day
    ih=ktime%hour
    imin=0
    length=0

!!$    if(lecmwf.or.lccsm)then
!!$       call cre_file_ma(prefix,ilen,iy,im,id,ih,length,filenam)
!!$    endif
!!$
!!$    if(lhadley.or.lecham.or.lecham2.or.lcanesm2.or.lcnrm)then
!!$       call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
!!$    endif
    call createFilename(prefix,ilen,iy,im,id,ih,length,filenam)

    call groploc(lundir,filenam,buf,nbuf,gcmddr)
    nx_grib = gcmddr%nlonhl
    ny_grib = gcmddr%nlathl

if(lecearth.and.lsstpath)then
    lundir2=62
!    prefix2=
    call readsstpath(prefix2,ilen2)
!write(6,'(a,a)')'prefix2 ',prefix2
!write(6,'(a,i3)')'ilen2 ',ilen2
    call createFilename(prefix2,ilen2,iy,im,id,ih,length,filenam2)
!write(6,'(a,a)')'filenam2 ',filenam2
    call groploc(lundir2,filenam2,buf,nbuf,gcmddr)
    if(mype==0) then
      print *,''
      print*,'NOTE,  this is a special version where we read'
      print*,'bias corrected SST from a separate run made by Klaus'
      print*,'The rest of the boundaries are read from the normal run'
      print *,'NOTE'
      print *,''
    endif
endif



    if(lhadley)then
       if(ltest .and. mype==0) write(6,*) ' grid values from script'

       west_model=0.
       south_model=-90.
       dlon_model=3.75
       dlat_model=2.50
       polon_model=0.
       polat_model=-90.

       gcmddr%aweshl=west_model
       gcmddr%alafhl=south_model
       gcmddr%dlonhl=dlon_model
       gcmddr%dlathl=dlat_model
       gcmddr%aplohl=polon_model
       gcmddr%aplahl=polat_model
    endif
    west_grib   = gcmddr%aweshl
    south_grib  = gcmddr%alafhl
    dlat_grib   = gcmddr%dlathl
    dlon_grib   = gcmddr%dlonhl
    polon_grib  = gcmddr%aplohl
    polat_grib  = gcmddr%aplahl

    if( abs(polon_grib)<zeps .and.abs(polat_grib)<zeps ) then
       polon_grib = 0._realkind
       polat_grib = -90._realkind
    endif


    if(lccsm)then
       if ( mype==0 ) then
          write(6,*) ' PRE_GETGRB org south_grib,dlat_grib = ',south_grib,dlat_grib
       endif
       south_grib = - south_grib
       dlat_grib = - dlat_grib
       if ( mype==0 ) then
          write(6,*) ' PRE_GETGRB corrected south_grib,dlat_grib = ',&
               south_grib,dlat_grib
       endif
    endif


    if(lhadley)then
       west_grib_u  = 1.875
       south_grib_u = -88.75
       west_grib_v  = west_grib_u
       south_grib_v = south_grib_u
    endif

    if(lecmwf.or.lecham2.or.lccsm.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl) then
       arakawa_grib = 'a'
    elseif(lhadley)then
       arakawa_grib = 'b'
    elseif(lecearth)then
       arakawa_grib = 'c'
    else
       if(mype==0)then
          write(6,*)'You have not defined the Arakawa grid for the lateral boundaries'
          write(6,*)'Have to stop'
          call stop_program('pre_getgrblow Arakawa?')
       endif
    endif


    call arakawastag(arakawa_grib,stagu_x_grib,stagu_y_grib, &
         stagv_x_grib,stagv_y_grib)


    !     Read remaining fields in another geometry and
    !     carry out necessary interpolations and adjustments

    !     Do not read time dependent climate fields from bd file 

    if(lecmwf.or.lhadley.or.lecham2.or.lecearth.or.lmpiesm) then
       ltsea   = .false.
       lfrice = .false.
    endif
    if(lccsm)then
       ltsea   = .false.
    endif
    if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl)then
       ltsea   = .false.
    endif



    call interpol_bdlow(klon,klat,lundir,lundir2, &
         ktime, & 
         RCAdom,&
         lunlsm,&
         ltsea,tsea_bound,&
         tice_bound,&
         frland,lfrice,frice_bound,                &
         nx_grib,ny_grib,               &
         west_grib,south_grib,dlat_grib,dlon_grib,&
         west_grib_u,south_grib_u,west_grib_v,south_grib_v, &
         polon_grib,polat_grib,                   &
         stagu_x_grib,stagu_y_grib,               &
         stagv_x_grib,stagv_y_grib)


    if(lecmwf.or.lhadley)then
       call grclloc(lundir2,mype)
    endif

    if(lecearth.and.lsstpath)then
       call grclloc(lundir2,mype)
    endif

    call grclloc(lundir,mype)

    return
  end subroutine pre_getgrblow





  subroutine readecicepath(ecice_path,ilen)
    use decomp
    implicit none
    character(len=132),intent(out)::ecice_path
    integer,intent(out)::ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /ecice/ecice_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=ecice)
    close(61)
    if(mype==0.and..not.done)then
       write(6,nml=ecice)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=ecice)
       close(1)
       done = .true.
    endif
    ilen=index(ecice_path,' ')-1
  end subroutine readecicepath

!!$  subroutine readboundpath(bound_path,ilen)
!!$    use decomp
!!$    implicit none
!!$    character(len=132),intent(out)::bound_path
!!$    integer,intent(out)::ilen
!!$    logical,save::done=.false.
!!$    namelist /bound/bound_path
!!$    ilen=1
!!$    open(66,file='namelists.dat',status='old',form='formatted')
!!$    read(66,nml=bound)
!!$    close(66)
!!$    if(mype==0.and..not.done)then
!!$       write(6,nml=bound)
!!$       done = .true.
!!$    endif
!!$    ilen=index(bound_path,' ')-1
!!$  end subroutine readboundpath



  subroutine interpol_bdlow(klon,klat,lundir,lundir2,&
       ktime, & 
       RCAdom,&
       lunlsm,            &
       ltsea,tsea_bound,                         &
       tice_bound,                                      &
       frland,lfrice,frice_bound,               &
       nx_grib,ny_grib,                 &
       west_grib,south_grib,dlat_grib,dlon_grib,  &
       west_grib_u,south_grib_u,west_grib_v,south_grib_v, &
       polon_grib,polat_grib,                     &
       stagu_x_grib,stagu_y_grib,                 &
       stagv_x_grib,stagv_y_grib)
    use decomp
    use util
    use grw1
    use IO
    use ctun, only: tseafr
    use confys
    use gcm
    use routmodvar
    use config
    use calendar

    implicit none
    type(domain),intent(in)::RCAdom
    integer,intent(in)::klon,klat,lundir
    type(time),intent(in)::ktime
    integer,intent(in):: lundir2

    real,intent(in):: west_grib_u,south_grib_u,west_grib_v,south_grib_v

    real(kind=realkind),intent(in):: frland(klon,klat)
!    real(kind=realkind) graviti
    real(kind=realkind),allocatable:: psrel(:,:) 

    integer lunlsm
    logical,intent(in)::ltsea,lfrice
    real(kind=realkind) tice_bound(klon,klat),tsea_bound(klon,klat),frice_bound(klon,klat)
    integer,intent(in)::nx_grib,ny_grib
    real(kind=realkind),intent(in):: west_grib,south_grib,dlat_grib,dlon_grib,polon_grib,polat_grib
    real(kind=realkind) dlat_gribi,dlon_gribi
    real(kind=realkind),intent(in):: stagu_x_grib,stagu_y_grib,  stagv_x_grib,stagv_y_grib



    integer::lunfis

    real(kind=realkind),allocatable::lat_m(:,:),lon_m(:,:),lat_u(:,:),lon_u(:,:),lat_v(:,:),lon_v(:,:)
    real(kind=realkind),allocatable,dimension(:,:)::zlat1_m,zlon1_m,zlat1_u,zlon1_u,zlat1_v,zlon1_v
    real(kind=realkind),allocatable,dimension(:,:)::zlat2_m,zlon2_m,zlat2_u,zlon2_u,zlat2_v,zlon2_v
    
    integer,allocatable:: i_m(:,:),j_m(:,:),i_u(:,:),j_u(:,:),i_v(:,:),j_v(:,:) 
    real(kind=realkind),allocatable:: w_x_m(:,:),w_y_m(:,:),w_x_u(:,:),w_y_u(:,:),w_x_v(:,:),w_y_v(:,:)
    real(kind=realkind),allocatable::  pa(:,:),pb(:,:),pc(:,:),pd(:,:), fi500ec(:,:)         

    integer::i,j,itype,ipar,ilev
!    integer::i,j,itype,ipar,ilev,no_corr,no_corr3
    integer::ierr

    real(kind=realkind)::eps,x_m,y_m,y_u,x_u,y_v,x_v

    logical::rotate
    integer::kk
    real(kind=realkind)::rmin,rmax,rmean

    real(kind=realkind)::lsmecham(nx_grib,ny_grib)
!    real(kind=realkind)::tempo(klon,klat)
!    real(kind=realkind)::td1echam(klon,klat)
!    real(kind=realkind)::td2echam(klon,klat)
!    real(kind=realkind)::td3echam(klon,klat)
!    real(kind=realkind)::td4echam(klon,klat)
!    real(kind=realkind)::dist,w1,w2,dec1,dec2,dec3,dec4,drc3,drc4,drc5
    logical,save::lfirstlsm=.true.

    data     eps/1.e-6_realkind/


    allocate(psrel(klon,klat))


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

    if(lecham2)then
       lfirstlsm=.true.
    endif

    !     1.1   calculate grid-point coordinates of non-staggered output field
    do j=1,klat
       lat_m(:,j) = RCAdom%south + real( jdatastart + j - 2,realkind )*RCAdom%dlat
    enddo
    do i=1,klon
       lon_m(i,:) = RCAdom%west + real( idatastart + i - 2,realkind )*RCAdom%dlon
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


    if(abs(RCAdom%polon-polon_grib)>eps.or. &
         abs(RCAdom%polat-polat_grib)>eps)then
       rotate = .true.

       !     de-rotate if needed coordinates from output grid geometry
       if( abs(RCAdom%polat+90._realkind)>eps ) then
          call  regrot(zlon1_m,zlat1_m,lon_m,lat_m, &
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

       if( abs(polat_grib+90._realkind)>eps ) then
          call  regrot(zlon1_m,zlat1_m,zlon2_m,zlat2_m, &
               klon,klat,&
               polon_grib,polat_grib,+1)                  
          call  regrot(zlon1_u,zlat1_u,zlon2_u,zlat2_u,  &
               klon,klat,&
               polon_grib,polat_grib,+1)                  
          call  regrot(zlon1_v,zlat1_v,zlon2_v,zlat2_v,  &
               klon,klat,&
               polon_grib,polat_grib,+1)
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

    dlat_gribi = 1.0_realkind/dlat_grib
    dlon_gribi = 1.0_realkind/dlon_grib
    do j=1,klat
       do i=1,klon
          if( zlon2_m(i,j)<west_grib )then
             zlon2_m(i,j) = zlon2_m(i,j) + 360._realkind      
          endif
          if( zlon2_u(i,j)<west_grib )then
             zlon2_u(i,j) = zlon2_u(i,j) + 360._realkind      
          endif
          if( zlon2_v(i,j)<west_grib )then            
             zlon2_v(i,j) = zlon2_v(i,j) + 360._realkind
          endif

          y_m = ( zlat2_m(i,j) - south_grib )*dlat_gribi + 1.0_realkind
          j_m(i,j)    = int( y_m )
          w_y_m(i,j)  = y_m - real(j_m(i,j),realkind)

          x_m = ( zlon2_m(i,j) - west_grib )*dlon_gribi + 1.0_realkind
          i_m(i,j)    = int( x_m )

          if ( i_m(i,j) > nx_grib ) then
             x_m = x_m - real(nx_grib,realkind)
             i_m(i,j) = i_m(i,j) - nx_grib
          endif
          w_x_m(i,j)  = x_m - real(i_m(i,j),realkind)

          if (lhadley) then
             y_u = ( zlat2_u(i,j) - south_grib_u )*dlat_gribi - &
                  stagu_y_grib + 1.0
          else
             y_u = ( zlat2_u(i,j) - south_grib )*dlat_gribi - stagu_y_grib + 1.0
          endif
          j_u(i,j)    = int( y_u )
          w_y_u(i,j)  = y_u - real(j_u(i,j),realkind)

          if (lhadley) then
             x_u = ( zlon2_u(i,j) - west_grib_u )*dlon_gribi - &
                  stagu_x_grib + 1.0
          else
             x_u = ( zlon2_u(i,j) - west_grib )*dlon_gribi - stagu_x_grib + 1.0
          endif
          i_u(i,j)    = int( x_u )
          if ( i_u(i,j) > nx_grib ) then
             x_u = x_u - real(nx_grib,realkind)
             i_u(i,j) = i_u(i,j) - nx_grib
          endif
          w_x_u(i,j)  = x_u - real(i_u(i,j),realkind)


          if (lhadley) then
             y_v = ( zlat2_v(i,j) - south_grib_v )*dlat_gribi - &
                  stagv_y_grib + 1.0
          else
             y_v = ( zlat2_v(i,j) - south_grib )*dlat_gribi - stagv_y_grib + 1.0
          endif
          j_v(i,j)    = int( y_v )
          w_y_v(i,j)  = y_v - real(j_v(i,j),realkind)


          if (lhadley) then
             x_v = ( zlon2_v(i,j) - west_grib_v )*dlon_gribi - &
                  stagv_x_grib + 1.0
          else
             x_v = ( zlon2_v(i,j) - west_grib )*dlon_gribi - stagv_x_grib + 1.0
          endif
          i_v(i,j)    = int( x_v )
          if ( i_v(i,j) > nx_grib ) then
             x_v = x_v - real(nx_grib,realkind)
             i_v(i,j) = i_v(i,j) - nx_grib
          endif
          w_x_v(i,j)  = x_v - real(i_v(i,j),realkind)
       enddo
    enddo

!    if(lhadley) then
!       itype = 105
!       ipar  = 006
!       ilev  = 000
!       lunfis = lundir

!       do j=1,klat
!          do i=1,klon
!             fis_grib(i,j) = orog_bound(i,j)
!          end do
!       end do
!       ierr=0
!    endif



    !     read tsd
    if(lecham2)then
       if(lecham5)then
!!$          itype = 111
!!$          ipar  = 85
!!$          ilev  = 3
!!$          call grrdloc_hint(                          &
!!$               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,    &
!!$               td1echam,td1echam,                      &
!!$               klon,klat,                              &
!!$               nx_grib,ny_grib,                        &
!!$               i_m,w_x_m,j_m,w_y_m,                    &
!!$               .false.,                                &
!!$               i_m,w_x_m,j_m,w_y_m,                    &
!!$               ierr)                                         
!!$          if(ierr/=0)call stop_program( 'could not read td1echam')


!!$          itype = 111                                   
!!$          ipar  = 85                                    
!!$          ilev  = 19                                    
!!$          call grrdloc_hint(                           &
!!$               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,    &
!!$               td2echam,td2echam,                      &
!!$               klon,klat,                              &
!!$               nx_grib,ny_grib,                        &
!!$               i_m,w_x_m,j_m,w_y_m,                    &
!!$               .false.,                                &
!!$               i_m,w_x_m,j_m,w_y_m,                    &
!!$               ierr)
!!$          if(ierr/=0)call stop_program( 'could not read td2echam')
       

!!$          itype = 111
!!$          ipar  = 85
!!$          ilev  = 78
!!$
!!$          call grrdloc_hint(                         &
!!$               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,   &
!!$               td3echam,td3echam,                     &
!!$               klon,klat,                             &
!!$               nx_grib,ny_grib,                       &
!!$               i_m,w_x_m,j_m,w_y_m,                   &
!!$               .false.,                               &
!!$               i_m,w_x_m,j_m,w_y_m,                   &
!!$               ierr)
!!$          if(ierr/=0)call stop_program( 'could not read td3echam')


!!$          itype = 111
!!$          ipar  = 85
!!$          ilev  = 268
!!$
!!$          call grrdloc_hint(                       &
!!$               lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
!!$               td4echam,td4echam,                   &
!!$               klon,klat,                           &
!!$               nx_grib,ny_grib,                     &
!!$               i_m,w_x_m,j_m,w_y_m,                 &
!!$               .false.,                             &
!!$               i_m,w_x_m,j_m,w_y_m,                 &
!!$               ierr)
!!$          if(ierr/=0)call stop_program( 'could not read td4echam')



!!$          !     depths (cm)
!!$          dec1=3
!!$          dec2=19
!!$          dec3=78
!!$          dec4=268
!!$          drc3=17.7
!!$          drc4=64.2
!!$          drc5=194.7
!!$
!!$          dist=1.0/(dec2-dec1)
!!$          w1=abs((dec2-drc3))*dist
!!$          w2=abs((dec1-drc3))*dist
!!$
!!$
!!$          dist=1.0/(dec3-dec2)
!!$          w1=abs((dec3-drc4))*dist
!!$          w2=abs((dec2-drc4))*dist
!!$
!!$          dist=1.0/(dec4-dec3)
!!$          w1=abs((dec4-drc5))*dist
!!$          w2=abs((dec3-drc5))*dist
!!$
!!$
!!$          graviti = 0.0065/9.81


          !     i.e. echam4
!!$          itype = 105
!!$          ipar  = 85
!!$          ilev  = 5
!!$
!!$          call grrdloc_hint(                       &
!!$               lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
!!$               tempo,tempo,                         &
!!$               klon,klat,                           &
!!$               nx_grib,ny_grib,                     &
!!$               i_m,w_x_m,j_m,w_y_m,                 &
!!$               .false.,                             &
!!$               i_m,w_x_m,j_m,w_y_m,                 &
!!$               ierr)                                 
!!$          if(ierr/=0)call stop_program( 'could not read tempo')



       endif

    endif



    if( .not.ltsea ) then
       if(lecham2)then
          !     
          if (ktime%hour == 0) then
             if (lfirstlsm .and. mype==0) then
                itype = 105
                ipar  = 081
                ilev  = 000
                if(lecham5)then
                   call gread(lunlsm,itype,ipar,real(ilev,realgribkind),&
                        lsmecham,nx_grib,ny_grib,ierr)
                else
                   call gread(lundir,itype,ipar,real(ilev,realgribkind), &
                        lsmecham,nx_grib,ny_grib,ierr)
                endif
                if(ierr/=0)call stop_program( 'could not read lsmecham')
                lfirstlsm=.false.
             endif
             itype = 102
             ipar  = 011
             ilev  = 000

             call grrdloc_hint_lsmmatch(       &
                  lundir,itype,ipar,real(ilev,realgribkind), &
                  tsea_bound,                         &
                  frland,lsmecham,              &
                  klon,klat,                    &
                  nx_grib,ny_grib,              &
                  i_m,w_x_m,j_m,w_y_m,          &
                  ierr)
             if(ierr/=0)call stop_program( 'could not read tsea')
          endif               ! hour=0
       endif ! lecham2

       if(lecmwf.or.lccsm)then
          if(.not.lecmwf)then
! i.e. ccsm
             itype = 105
             ipar  = 011
             ilev  = 000

             call grrdloc_hint(                      &
                  lundir,itype,ipar,real(ilev,realgribkind),lecmwf,&
                  tsea_bound,tsea_bound,                          &
                  klon,klat,                          &
                  nx_grib,ny_grib,                    &
                  i_m,w_x_m,j_m,w_y_m,                &
                  .false.,                            &
                  i_m,w_x_m,j_m,w_y_m,                &
                  ierr)
             if(ierr/=0)call stop_program( 'could not read tsea')

             !     &**  get ice_conc ( 0 or 1 ) fron tsea !!!!

             do j=1,klat
                do i=1,klon
                   !     if(tsea_bound(i) < 271.36 .and. ! < -1.8c
                   if(tsea_bound(i,j) < tseafr .and. &
                        frland(i,j) < 0.30_realkind ) then ! sea
                      frice_bound(i,j) = 1._realkind
                      tice_bound(i,j)=tsea_bound(i,j)
                   else
                      frice_bound(i,j) = 0._realkind
                      tice_bound(i,j)=tmelt
                   endif
                enddo
             enddo
          else !read from tice_asim/tice_*-files
! lecmwf
             itype = 001
             ipar  = 034
             ilev  = 0
             call grrdloc_hint(                                &
                  lundir2,itype,ipar,real(ilev,realgribkind),lecmwf,         &
                  tsea_bound,tsea_bound,                                    &
                  klon,klat,                                    &
                  nx_grib,ny_grib,                              &
                  i_m,w_x_m,j_m,w_y_m,                          &
                  .false.,                                      &
                  i_m,w_x_m,j_m,w_y_m,                          &
                  ierr)    
             if(ierr/=0)call stop_program( 'could not read tsea_bound')

             !read from tice_asim/tice_*-files  
             ilev  = 000                                         
             ipar  = 035      !kw test                           
             itype = 112      !kw test                           
             call grrdloc_hint(                                 &
                  lundir2,itype,ipar,real(ilev,realgribkind),lecmwf,         &
                  tice_bound,tice_bound,                                    &
                  klon,klat,                                    &
                  nx_grib,ny_grib,                              &
                  i_m,w_x_m,j_m,w_y_m,                          &
                  .false.,                                      &
                  i_m,w_x_m,j_m,w_y_m,                          &
                  ierr)
             if(ierr/=0)call stop_program( 'could not read tice_bound'            )
             
             tice_bound=min(tmelt,tice_bound)
             tsea_bound=max(tseafr,tsea_bound)

          endif ! ecmwf
       endif ! ecmwf or ccsm
       if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
          itype = 001
          ipar  = 235
          ilev  = 0
          call grrdloc_hint(                                &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,         &
               tsea_bound,tsea_bound,                                    &
               klon,klat,                                    &
               nx_grib,ny_grib,                              &
               i_m,w_x_m,j_m,w_y_m,                          &
               .false.,                                      &
               i_m,w_x_m,j_m,w_y_m,                          &
               ierr)    
!          if(mype.eq.0) write(6,*) ' tsea in from canesm2'
          if(ierr/=0)call stop_program( 'could not read tsea_bound')
!!cau110629
!!cau110629   dummy test
!          no_corr=0
!          no_corr3=0
!          do j=1,klat
!             do i=1,klon
!                if(tsea_bound(i,j)<-800._realkind) then
!!                   if(frland(i,j)<0.3_realkind) then
!!                      tsea_bound(i,j)=273._realkind
!!                      no_corr3=no_corr3+1
!!                   else
!!                      tsea_bound(i,j)=273._realkind
!!                      no_corr=no_corr+1
!!                   endif
!                   no_corr=no_corr+1
!                endif
!             enddo
!          enddo
       endif ! canesm or cnrm or noresm or hadgem or miroc5 or ipsl or gfdl
       if(lecearth.or.lmpiesm) then
          itype = 102
          ipar  = 11
          ilev  = 0
if(lecearth.and.lsstpath)then
          call grrdloc_hint(                                &
               lundir2,itype,ipar,real(ilev,realgribkind),lecmwf,         &
               tsea_bound,tsea_bound,                                    &
               klon,klat,                                    &
               nx_grib,ny_grib,                              &
               i_m,w_x_m,j_m,w_y_m,                          &
               .false.,                                      &
               i_m,w_x_m,j_m,w_y_m,                          &
               ierr)    
  if(mype==0) then
    print *,''
    print*,'NOTE,  this is a special version where we read'
    print*,'bias corrected SST from  a separate run made by Klaus'
    print*,'The rest of the boundaries are read from the normal run'
    print *,'NOTE'
    print *,''
  endif
else
          call grrdloc_hint(                                &
               lundir,itype,ipar,real(ilev,realgribkind),lecmwf,         &
               tsea_bound,tsea_bound,                                    &
               klon,klat,                                    &
               nx_grib,ny_grib,                              &
               i_m,w_x_m,j_m,w_y_m,                          &
               .false.,                                      &
               i_m,w_x_m,j_m,w_y_m,                          &
               ierr)    
endif
          if(lecearth)&
          call meanmnmx('tsea_bound from EC-EARTH read',tsea_bound,klon,klat)
          if(lmpiesm)&
          call meanmnmx('tsea_bound from MPIESM read',tsea_bound,klon,klat)
          if(ierr/=0)call stop_program( 'could not read tsea_bound read')

       endif ! EC-EARTH
    endif ! ltsea

    !     read frice_bound
    if(lccsm)then
       itype = 105
       ipar  = 091
       ilev  = 000

       call grrdloc_hint(                         &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf,   &
            frice_bound,frice_bound,                           &
            klon,klat,                             &
            nx_grib,ny_grib,                       &
            i_m,w_x_m,j_m,w_y_m,                   &
            .false.,                               &
            i_m,w_x_m,j_m,w_y_m,                   &
            ierr)
       if(ierr/=0)call stop_program( 'could not read frice_bound')
    endif
    if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl)then
       itype = 001
       ipar  = 091
       ilev  = 000

       call grrdloc_hint(                         &
            lundir,itype,ipar,real(ilev,realgribkind),lecmwf,   &
            frice_bound,frice_bound,                           &
            klon,klat,                             &
            nx_grib,ny_grib,                       &
            i_m,w_x_m,j_m,w_y_m,                   &
            .false.,                               &
            i_m,w_x_m,j_m,w_y_m,                   &
            ierr)
       if(ierr/=0)call stop_program( 'could not read frice_bound')
       if(lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
          do j=1,klat
          do i=1,klon
          if(frice_bound(i,j)<0_realkind) frice_bound(i,j)=0._realkind
          end do
          end do
!
! should do mod_sst! (mod_ice instead?????????
!
          if(mype==0)write(6,*)'changed frice_bound'
          if(mype==0)write(6,*)'changed missing value -999. to 0.'
          if(mype==0)write(6,*)'should do mod_ice in gcm instead???'
       endif
    endif
    if(lecmwf)then
       if( .not.lfrice ) then
          itype  = 1
          ipar  = 031
          ilev  = 000
          call grrdloc_hint(                        &
               lundir2,itype,ipar,real(ilev,realgribkind),lecmwf, &
               frice_bound,frice_bound,                          &
               klon,klat,                            &
               nx_grib,ny_grib,                      &
               i_m,w_x_m,j_m,w_y_m,                  &
               .false.,                              &
               i_m,w_x_m,j_m,w_y_m,                  &
               ierr)
          if(ierr/=0)call stop_program( 'could not read frice_bound')
          frice_bound=min(.9999_realkind,max(.0001_realkind,frice_bound))
       endif
    else
       if( .not.lfrice ) then
          itype = 102
          ipar  = 091
          ilev  = 000
          if(lecham2)then
             if (ktime%hour == 0) then
                if (lfirstlsm .and. mype==0) then
                   itype = 105
                   ipar  = 081
                   ilev  = 000
                   if(lecham5)then
                      call gread(lunlsm,itype,ipar,real(ilev,realgribkind), &
                           lsmecham,nx_grib,ny_grib,ierr)
                   else
                      call gread(lundir,itype,ipar,real(ilev,realgribkind), &
                           lsmecham,nx_grib,ny_grib,ierr)
                   endif
                   if(ierr/=0)call stop_program( 'could not read lsmecham')
                endif
                itype = 102
                ipar  = 091
                ilev  = 000
                call grrdloc_hint_lsmmatch(        &
                     lundir,itype,ipar,real(ilev,realgribkind),  &
                     frice_bound,                         &
                     frland,lsmecham,               &
                     klon,klat,                     &
                     nx_grib,ny_grib,               &
                     i_m,w_x_m,j_m,w_y_m,           &
                     ierr)
                if(ierr/=0)call stop_program( 'could not read frice_bound')
             endif            ! hour=0
          else
! not echam2
             call grrdloc_hint(                       &
                  lundir,itype,ipar,real(ilev,realgribkind),lecmwf, &
                  frice_bound,frice_bound,                         &
                  klon,klat,                           &
                  nx_grib,ny_grib,                     &
                  i_m,w_x_m,j_m,w_y_m,                 &
                  .false.,                             &
                  i_m,w_x_m,j_m,w_y_m,                 &
                  ierr)
             if(ierr/=0)call stop_program( 'could not read frice_bound')
             if(lmpiesm) then
                do j=1,klat
                do i=1,klon
                if(frice_bound(i,j)<0_realkind) frice_bound(i,j)=0._realkind
                end do
                end do
                if(mype==0)write(6,*)'changed frice_bound'
                if(mype==0)write(6,*)'changed missing value -999. to 0.'
             endif
             !     
          endif
       endif ! not lfrice
    endif ! not lecmwf



    !     for routing
    if(use_routing) then
       itype = 001
       ipar  = 219
       ilev  = 000
       call grrdloc_global(                                        &
            lundir,itype,ipar,real(ilev,realgribkind),                           &
            runoff_bound,                                           &
            nx_grib,ny_grib,                                        &
            ierr)                                                    
       if(ierr/=0)call stop_program( 'could not read runoff_bound')

       rmin=minval(runoff_bound)                                     
       rmax=maxval(runoff_bound)                                     
       rmean=sum(runoff_bound)/real(nx_grib*ny_grib,realkind)                          

       if(mype==0)then                                             
          write(6,661) ktime%hour,nx_grib,ny_grib,klon_bound,klat_bound,  &
               rmin,rmax,rmean
       endif

       itype = 001
       ipar  = 160
       ilev  = 000
       call grrdloc_global(                                         &
            lundir,itype,ipar,real(ilev,realgribkind),                            &
            runoff_bound,                                            &
            nx_grib,ny_grib,                                         &
            ierr)                                                          
       if(ierr/=0)call stop_program( 'could not read runoff_bound')
       rmin=minval(runoff_bound)                                      
       rmax=maxval(runoff_bound)                                      
       rmean=sum(runoff_bound)/real(nx_grib*ny_grib,realkind)                      

       write(6,662) ktime%hour,nx_grib,ny_grib,klon_bound,klat_bound,kk,   &
            rmin,rmax,rmean
661    format(' echam5 river  ',6i6,3f10.6)
662    format(' echam5 runoff ',6i6,3f10.6)
    endif
    !     end for routing




    !     4.0 vertical interpolation


    deallocate(lat_m, lon_m, lat_u, lon_u,lat_v,lon_v,zlat1_m,zlon1_m,  &
         zlat1_u,zlon1_u,zlat1_v,zlon1_v,zlat2_m,zlon2_m,zlat2_u,       &
         zlon2_u,zlat2_v,zlon2_v)

    deallocate(i_m,j_m, i_u, j_u, i_v,j_v)

    deallocate(w_x_m, w_y_m, w_x_u, w_y_u,w_x_v,w_y_v)

    deallocate(pa, pb, pc, pd, fi500ec)
    deallocate(psrel)

    return                                                               

  end subroutine interpol_bdlow



  subroutine interpol_sea (klon, klat,  nbdnum, lundir, &
       year, month, day, hour, &
       RCAdom,&
       ltsea, tsea, lfrice, frice,&
       nx_grib, ny_grib, west_grib, &
       south_grib, dlat_grib, dlon_grib, polon_grib, polat_grib, gcmddr)

    use decomp  
    use util
    use IO
    use modddr

    implicit none 
    type(domain),intent(in)::RCAdom
    integer :: klon, klat,  lundir, ierr, year, month, &
         day, hour
    type(ddr),intent(inout)::gcmddr
    integer :: nbdnum  

    logical :: ltsea, lfrice  

    real(kind=realkind) :: tsea(klon,klat), frice(klon,klat)
    integer :: nx_grib, ny_grib  
    real(kind=realkind) :: west_grib, south_grib, dlat_grib, dlon_grib, polon_grib, &
         polat_grib

    real(kind=realkind) :: stagu_x_grib, stagu_y_grib, stagv_x_grib, stagv_y_grib  



    character (len=80) :: mes  
    integer::ilen, istat, iy, im, id, ih, imin, length  
    character (len=132) :: prefix, filenam  
    integer :: jpbuf, lenbuf  
    parameter (jpbuf = 70000)  

    real(kind=realgribkind) :: buf (jpbuf)  

    real(kind=realkind) :: lat_m (klon, klat), lon_m (klon, klat), &
         zlat1_m (klon, klat), zlon1_m (klon, klat), &
         zlat2_m (klon, klat), zlon2_m (klon, klat)

    integer :: i_m (klon, klat), j_m (klon, klat)  

    real(kind=realkind) :: w_x_m (klon, klat), w_y_m (klon, klat)  
    integer:: lunsea, ddc, hhc, nxny  
    real(kind=realkind)::plev, psum  

    real(kind=realkind) :: frice_grib (klon, klat), tsea_grib (klon,klat)

    integer:: i, j, k, itype, ipar, ilev, jlev
    real(kind=realkind)::zlat, zlon, eps, x_m, y_m

    logical::rotate,  lopen  
    integer :: kk  

    real(kind=realkind) ::  rmean  
    data eps / 1.e-6_realkind /  



    lenbuf = jpbuf  
    !     1.1   Calculate grid-point coordinates of non-staggered output fie
    do j = 1, klat  
       zlat = RCAdom%south + real (jdatastart + j - 2,realkind) * RCAdom%dlat  
       do i = 1, klon  
          zlon = RCAdom%west + real (idatastart + i - 2,realkind) * RCAdom%dlon  
          lat_m (i, j) = zlat  
          lon_m (i, j) = zlon  
       enddo
    enddo
    !     1.3   Rotate coordinates if needed

    if (abs (RCAdom%polon - polon_grib) >eps.or.abs (RCAdom%polat - polat_grib) &
         >eps) then
       rotate = .true.  
       !     1.3.1  De-rotate if needed coordinates from output grid geometry


       if (abs (RCAdom%polat + 90._realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, lon_m, lat_m, klon, &
               klat,  RCAdom%polon, RCAdom%polat, - 1)
       else  
          do j = 1, klat  
             do i = 1, klon  
                zlat1_m (i, j) = lat_m (i, j)  
                zlon1_m (i, j) = lon_m (i, j)  
             enddo
          enddo
       endif
       !     1.3.2  Rotate if needed coordinates to input grid geometry

       if (abs (polat_grib + 90._realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, zlon2_m, zlat2_m, klon, &
               klat,  polon_grib, polat_grib,+ 1)
       else  
          do j = 1, klat  
             do i = 1, klon  
                zlat2_m (i, j) = zlat1_m (i, j)  
                zlon2_m (i, j) = zlon1_m (i, j)  
             enddo
          enddo
       endif
       !     1.3.3  No rotation needed, just copy coordinates
    else  
       rotate = .false.  
       do j = 1, klat  
          do i = 1, klon  
             zlat2_m (i, j) = lat_m (i, j)  
             zlon2_m (i, j) = lon_m (i, j)  
          enddo
       enddo
    endif
    !     1.4  calculate indeces and weights needed for horizontal interpola

    do j = 1, klat  
       do i = 1, klon  
          if(zlon2_m (i, j)<west_grib)zlon2_m(i, j)=zlon2_m(i,j)+360._realkind
          y_m = (zlat2_m (i, j) - south_grib) / dlat_grib + 1.0_realkind  
          j_m (i, j) = int(y_m)  
          w_y_m (i, j) = y_m - real (j_m (i, j),realkind )  
          x_m = (zlon2_m (i, j) - west_grib) / dlon_grib + 1.0_realkind  
          i_m (i, j) = int (x_m)  
          w_x_m (i, j) = x_m - real (i_m (i, j),realkind )  
       enddo

    enddo
    !     2.4  read tsea, frice from Hadley  NOTE  half grid !!!

    lunsea = lundir + 5  
    !     always use same unit number,
    !     to avoid problem with reopening file with new number (1. 19, 2. 21
    !     don't know why this happens
    !     and it does not happen in Anders's atmosphere-only version
    lunsea = 19  
    if (ltsea.or.lfrice) then  
       if ( (day==16.and.hour==0) .or.nbdnum==1) then  
          if (mype==0) then  
             call readoceanpath(prefix,ilen)
             prefix = prefix (1:ilen) //'HAD_OCEAN'  
             ilen = ilen + 9  
             iy = year  
             im = month  
             id = day  
             ih = hour  
             imin = 0  
             length = 0  
             ih = 00  
             id = 16  
             call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)

          endif
          call groploc(lunsea,filenam, buf, lenbuf,gcmddr)
          nx_grib = gcmddr%nlonhl
          ny_grib = gcmddr%nlathl

          lopen = .true.  
       else  
          lopen = .false.  
          !  if day hour
       endif
       if (lopen) then  
          if (ltsea) then  
             itype = 105  
             ipar = 091  
             ilev = 000  
             call grrdloc_hint(lunsea, itype, ipar, real(ilev,realgribkind), &
                  lecmwf, frice_grib, frice_grib, klon, klat, &
                  nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, .false., i_m, &
                  w_x_m, j_m, w_y_m, ierr)
             if(ierr/=0)call stop_program( 'could not read frice_grib')
          endif
          if (lfrice) then  
             itype = 105  
             ipar = 011  
             ilev = 000  
             call grrdloc_hint (lunsea, itype, ipar, real(ilev,realgribkind), &
                  lecmwf, tsea_grib, tsea_grib, klon, klat, &
                  nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, .false., i_m, &
                  w_x_m, j_m, w_y_m, ierr)
             if(ierr/=0)call stop_program( 'could not read tsea_grib')
          endif

          psum = 0._realkind 
          rmean = 0._realkind  
          do j = 1, klat  
             do i = 1, klon  
                if (ltsea) then  
                   tsea(i,j) = tsea_grib (i, j)  
                   rmean = rmean + tsea_grib (i, j)  
                endif
                if (lfrice) then  
                   frice(i,j) = frice_grib (i, j)  
                   if (frice(i,j) >0._realkind.and.frice(i,j)<0.005_realkind)then
                      frice(i,j) = 0._realkind  
                   endif
                   if (frice_grib (i, j) >1.0_realkind) then  
                      print *, 'frice skum', mype, i, j, frice_grib (i,j)
                   endif
                   psum = psum + frice(i,j)  
                endif
             enddo
          enddo
          call grclloc (lunsea, mype)  
          
       endif
    endif

    return  
  end subroutine interpol_sea



end module surface_bc
