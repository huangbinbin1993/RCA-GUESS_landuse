module climate
  use timetype
  use gcm
  use comhkp
  use decomp
  use config
  use modddr
  use calendar
  use grw1
  use rcaDomainMod
  use IO
  use util
  use ctun, only:tseafr  
  use confys
  use flake
  use surface
  use mod_grib,only:realgribkind
  implicit none
  private
  integer,save:: year1,month1,year2,month2,ddc,hhc
  real(kind=realkind),save::weight1,weight2

  public setInitialCondPhys,read_clim_year,int_clim
contains

  subroutine setInitialCondPhys(klon,klat,tsclim_years, &
       svar_surf,svarsm,prog_lakes,eco, &
       intime, west_in,south_in,dlon_in,dlat_in,polon_in,polat_in  )

    implicit none

    integer,intent(in):: klon,klat
    real(kind=realkind),intent(out)::tsclim_years(klon,klat)
    real(kind=realkind),allocatable,save::clim(:,:,:)
    real(kind=realkind),intent(in):: west_in,south_in,dlon_in,dlat_in,polon_in,polat_in  
    type(time),intent(in)::intime
    real(kind=realkind)::svarsm(:,:,:),prog_lakes(:,:,:,:),svar_surf(:,:,:)
    real(kind=realkind)::eco(:,:,:)

    !local arrays
    real(kind=realkind)::tsdcli(klon,klat),swdcli(klon,klat), tsea(klon,klat),tsm(klon,klat)
    real(kind=realkind)::frice(klon,klat), tsd(klon,klat), swd(klon,klat),snm(klon,klat),swm(klon,klat)
    integer,parameter::jpbuf=700000
    integer:: lenbuf
    integer::nsvarns,nlakevarns, nlaketypeIn,meco
    type(ddr)::myddr
    logical,save::firstTime=.true.


    nsvarns = ubound(svarsm,3)
    nlakevarns = ubound(prog_lakes,4)
    nlaketypeIn = ubound(prog_lakes,3)
    meco = ubound(eco,3)

    allocate(clim(klon,klat,18)) 



    lenbuf = jpbuf
    ddc = 01
    hhc = 00
    year1 = intime%year
    year2 = intime%year

    if ( intime%day < 16 ) then
       month1=intime%month-1
       month2=intime%month
       if ( month1 == 0 ) then
          month1=12
          year1 = intime%year-1
       endif
    else 	
       month1=intime%month
       month2=intime%month+1
       if ( month2 == 13 ) then
          month2=1
          year2 =intime%year+1
       endif
    endif


    !     Open the time-dependent surface data file

    !     fields needed at initial time and, if option
    !     lclupdate is selected, also when boundary data are
    !     changed

    !     read two months climate data


    call read_clim(1,klon,klat,clim,year1,month1,ddc,hhc,&
         west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)
    call read_clim(2,klon,klat,clim,year2,month2,ddc,hhc,&
         west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)


    !     interpolate between two months
    call int_clim(intime%day,month1,month2,weight1,weight2)
    call mix_clim(mype,weight1,weight2,klon,klat,&
         tsdcli,swdcli,tsea,frice,tsd,swd,snm,tsm,swm,clim)

!!$    if(use_oasis)then
!!$       do jy=1,klat
!!$          do jx=1,klon 
!!$             if(lcounttice(jx,jy)>0.5) then
!!$                if(lecmwf)then
!!$                   svarsm(jx,jy,14) = tice_bound(jx,jy,1)+ timrat*(tice_bound(jx,jy,2)-tice_bound(jx,jy,1))
!!$                endif
!!$             endif
!!$          enddo
!!$       enddo
!!$    else
!!$       if(lecmwf)then
!!$          svarsm(:,:,14)=tice_bound(:,:,1)+ timrat*(tice_bound(:,:,2)-tice_bound(:,:,1))
!!$       endif
!!$    endif

    if(firstTime)then
       call iniphy(klon,klat,intime,nsvarns,nlakevarns,nlaketypeIn,&
            tsm,swm,snm,tsd,swd,tsea,tsdcli,eco(:,:,35),swdcli,eco(:,:,10), &
            frice,eco(:,:,15),eco(:,:,16),&
            svar_surf(:,:,62),svarsm,eco(:,:,47:49),prog_lakes, &
            west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)


       call read_clim_many_years(klon,klat,tsclim_years,myddr )
       firstTime=.false.
    endif
    deallocate(clim)

    return
  end subroutine setInitialCondPhys

  subroutine grrdloc_clim(lun,type,param,alev,f,&
       klon,klat,nx_grib,ny_grib,iind,wghx,jind,wghy,ierr)
    implicit none

    integer:: lun,type,param,klon,klat,nx_grib,ny_grib

    real(kind=realgribkind):: alev
    real(kind=realkind)::f(klon,klat)
    integer:: iind(klon,klat),jind(klon,klat)
    real(kind=realkind)::  wghx(klon,klat),wghy(klon,klat)

    real(kind=realkind),allocatable,dimension(:,:):: w

    integer:: i,j,ierr,ii,jj
    integer:: iip,jjp

#ifdef MPI_SRC
#include"mpif.h"
    integer:: mpierr
#endif
    allocate(w(nx_grib,ny_grib))
    !     Read the field with PE 0 and distribute the error code
    ierr = -1
    if( mype==0 ) then
       call gread(lun,type,param,alev,w,nx_grib,ny_grib,ierr)
    endif
#ifdef MPI_SRC
    call mpi_bcast(ierr,1,MPI_INTEGER,0,localComm,mpierr)
#endif
    if(ierr/=0) return

#ifdef MPI_SRC
    call mpi_bcast(w,nx_grib*ny_grib,REALTYPE,0,localComm,mpierr)
#endif

    !     Do the interpolation 
    do j=1,klat
       do i=1,klon
          ii = iind(i,j)
          jj = jind(i,j)
          iip = 1+mod(ii,nx_grib) !ii+1 taking into account periodicity
          jjp = min(jj+1,ny_grib) !jj+1 taking into account periodicity

          f(i,j) = (1.0_realkind-wghx(i,j))*(1.0_realkind-wghy(i,j))*w(ii,jj)    + &
               (1.0_realkind-wghx(i,j))*wghy(i,j)     *w(ii,jjp)  + &
               wghx(i,j)*(1.0_realkind-wghy(i,j))*w(iip,jj)  + &
               wghx(i,j)*wghy(i,j)     *w(iip,jjp)

       enddo
    enddo
    deallocate(w)
    return

  end subroutine grrdloc_clim



  subroutine int_clim(day,m1,m2,w1,w2)
    !     
    !     971028 Ulf Hansson, Rossby Centre
    !     020409 - change by Ulf Hansson, introduced ECHAM2
    !     
    !     interpolate between 2 sets of monthly climate datasets
    !     do it daily
    !     get a new month at the middle of a month
    !     and throw away the oldest one at the same time
    !     generally 1 corresponds to the earlier month
    !     and       2 to the later month

    !     
    !     nodpm   no of days per month
    !     midmonx middle of month
    !     wx      weight
    !     mx      month
    !     nodtp   no of days in the period (NOD This Period)
    !     lnewm   set true when you enter the middle of month,
    !     set false at beginning of month

    implicit none

    integer:: nodpm(12)
    integer:: day,m1,m2,midmon,midmon1,midmon2,nodtp,iw1,iw2
    real(kind=realkind):: w1,w2

    if(lecham.or.lecham2)then
       if(lecham5)then
          nodpm=(/31,28,31,30,31,30,31,31,30,31,30,31/)
       else
          nodpm = 30
       endif
    elseif(lhadley.or.lhadgem)then
       nodpm = 30
    elseif(lecmwf.or.lccsm.or.lecearth.or.lcanesm2.or.lcnrm.or.lnoresm.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
       nodpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    else
       if (mype == 0) then
          write(6,*) 'calendar not defined for this set of boundary data'
          write(6,*) 'have to stop in int_clim (grdy/intclim.f)'
       endif
       call stop_program('calendar not defined for this set of boundary data')
    endif
    !     check leap year
    !     skip for Hadley
    !     if (month == 2) then
    !     if ( (mod(year,4) == 0 .and. mod(year,100) /= 0) .or.
    !     x    mod(year,400) == 0) then
    !     nodpm(2)=29
    !     else
    !     nodpm(2)=28
    !     endif
    !     endif
    !     temporary return as the  test program always tries 31 days per month
    !     if (day > nodpm(month)) return

    !     calculate  weights
    !     midmon=nodpm(month)/2+1

    midmon=16

    midmon1=nodpm(m1)/2+1 
    midmon2=nodpm(m2)/2+1  

    !     number of days this period
    nodtp=nodpm(m1)-midmon1+1+midmon2-1
    if(day<midmon) then
       iw1=midmon1-day
       iw2=nodtp-iw1
    else
       iw1=nodtp-(day-midmon1)
       iw2=nodtp-iw1
    endif

    w1=real(iw1,realkind)/real(nodtp,realkind)
    w2=real(iw2,realkind)/real(nodtp,realkind)   

    return
  end subroutine int_clim



  subroutine mix_clim( mype, w1, w2, klon, &
       klat, tsdcli, swdcli, tsea, frice, tsd, swd,snm,tsm,swm, clim)

    implicit none  

    integer,intent(in) :: klon, klat,  mype  
    real(kind=realkind),intent(in) :: w1, w2  

    real(kind=realkind) :: tsdcli (klon, klat), swdcli (klon, klat),tsm(klon,klat), &
         tsea (klon, klat), frice (klon, klat),&
         tsd(klon, klat), swd (klon,klat),snm(klon,klat),swm(klon,klat), &
         clim (klon,  klat, 18)
    integer ::  kk
    logical::debug = .false.  

    kk = 9  
    tsdcli = w1*clim(:,:, 1)+w2*clim(:,:,kk+1)  
    swdcli = w1*clim(:,:, 2)+w2*clim(:,:,kk+2)  
    if (.not.lhadley.and..not.lecham) then  
       tsea  = w1*clim(:,:,3)+w2*clim(:,:,kk+3)  
       frice = w1*clim(:,:,4)+w2*clim(:,:,kk+4)  
    endif

    tsd = w1*clim(:,:,5)+w2*clim(:,:,kk+5)  
    swd = w1*clim(:,:,6)+w2*clim(:,:,kk+6) 
    snm = w1*clim(:,:,7)+w2*clim(:,:,kk+7) 
    tsm = w1*clim(:,:,8)+w2*clim(:,:,kk+8) 
    swm = w1*clim(:,:,9)+w2*clim(:,:,kk+9) 

    if (mype==0.and.debug) then  
       write (6, * ) ' mix_clim. interp. : tsdcli, swdcli'  
       write (6, * ) ' mix_clim. interp. ; tsd, swd'  
       if (.not.lhadley.and..not.lecham) then  
          write (6, * ) ' mix_clim. interp. : tsea, frice'  
       endif
    endif
    return  

  end subroutine mix_clim



  subroutine read_clim(kcall,klon,klat,clim,year,month,day,hour,&
       west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)
    implicit none
    integer:: kcall,klon,klat,year,month,day,hour
    real(kind=realkind):: clim(klon,klat,18)

    real(kind=realkind),intent(in)::west_in,south_in,dlon_in,dlat_in,polon_in,polat_in

    !     kcall    tsdcli   swdcli   tsea   frac.ice   tsd    swd  snm  tsm  swm 
    !     1         1         2       3       4        5      6     7    8    9
    !     2         10        11      12     13       14     15     16   17   18


    integer:: ilen,iy,im,id,ih,imin,length
    integer, parameter::jpbuf=700000
    real(kind=realgribkind) buf(jpbuf)
    integer ::lenbuf,nxny
    integer ::lundir
    character(len=132):: prefix,filenam

    integer:: ii_tsclim,ii_tsd,ii_swdclim,ii_swd,ii_tsea,ii_frice,ii_snm,ii_tsm,ii_swm
    integer:: nx_grib,ny_grib
    real(kind=realkind)::    west_grib,south_grib,dlon_grib,dlat_grib,polon_grib,polat_grib



    type(ddr),intent(inout)::myddr

    lenbuf = jpbuf



    !cau1106    if(lhadley.or.lecham)then
    if(lecham)then
       lundir=60
       call stop_program('need to read from a global file instead of a local cl-file')
    endif

    if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
       !     Open climate file from ECMWF ( in lat/long-grid )
       nx_grib     = 360
       ny_grib     = 181
       west_grib   = -180.0_realkind
       south_grib  = -90.0_realkind
       dlon_grib   = 1.0_realkind
       dlat_grib   = 1.0_realkind
       polon_grib  = 0.0_realkind
       polat_grib  = -90.0_realkind

       lundir=63

       if (mype==0) then
          call readgcpath(prefix,ilen)
          prefix = prefix(1:ilen)//'clim'
          ilen=ilen+4
          iy=year
          im=month
          id=1 !day
          ih=0 !hour
          imin=0
          length=0
          call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
       endif

       nxny=nx_grib*ny_grib
       call groploc(lundir,filenam,buf,nxny,myddr)
       nx_grib = myddr%nlonhl
       ny_grib = myddr%nlathl


       ii_tsclim = (kcall-1)*9+1
       ii_swdclim = (kcall-1)*9+2
       ii_tsea = (kcall-1)*9+3
       ii_frice = (kcall-1)*9+4
       ii_tsd  = (kcall-1)*9+5
       ii_swd = (kcall-1)*9+6
       ii_snm = (kcall-1)*9+7
       ii_tsm = (kcall-1)*9+8
       ii_swm = (kcall-1)*9+9

       call interpol_clim(klon,klat, &
            lundir,&
            west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,& 
            clim(:,:,ii_tsclim),clim(:,:,ii_tsd),clim(:,:,ii_swdclim),&
            clim(:,:,ii_swd), clim(:,:,ii_tsea),clim(:,:,ii_frice), &
            clim(:,:,ii_snm),clim(:,:,ii_tsm),clim(:,:,ii_swm),&
            nx_grib,ny_grib,west_grib,south_grib,dlat_grib,dlon_grib,&
            polon_grib,polat_grib)


       call grclloc(lundir,mype)
    endif
    !   ecmwf

    return
  end subroutine read_clim

  subroutine readgcpath(global_clim_path,ilen)

    implicit none
    character(len=132)::global_clim_path
    integer:: ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist/global_clim/global_clim_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=global_clim)
    close(61)
    if(mype==0.and..not.done)then
       write(6,nml=global_clim)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=global_clim)
       close(1)
       done = .true.
    endif
    ilen=index(global_clim_path,' ')-1
  end subroutine readgcpath


  subroutine interpol_clim (klon, klat, lundir,  &
       west, south, dlon, dlat, polon, polat, &
       tsclim, tsd,swdclim, &
       swd,tsea,frice,snm,tsm,swm, &
       nx_grib, ny_grib, west_grib, &
       south_grib, dlat_grib, dlon_grib, polon_grib, polat_grib)


    implicit none  
    integer :: klon, klat, lundir, ierr


    real(kind=realkind) :: west, south, dlat, dlon, polon, polat, r2, r3  
    real(kind=realkind) :: tsea(klon,klat), tsclim(klon,klat)
    real(kind=realkind)::frice(klon,klat), tsd(klon,klat), swdclim(klon,klat)
    real(kind=realkind)::swd(klon,klat),swm(klon,klat)
    real(kind=realkind)::snm(klon,klat),tsm(klon,klat)
    logical :: lwet_scal  
    integer :: nx_grib, ny_grib  
    real(kind=realkind) :: west_grib, south_grib, dlat_grib, dlon_grib  
    real(kind=realkind) :: polon_grib, polat_grib  

    real(kind=realkind), allocatable, dimension (:,:)::lat_m,lon_m,zlat1_m,zlon1_m,zlat2_m,zlon2_m
    integer , allocatable, dimension (:, :) ::i_m, j_m  
    real(kind=realkind), allocatable, dimension (:, :) ::w_x_m, w_y_m  


    integer:: i, j, k, itype, ipar, ilev,  iut  
    real(kind=realkind):: zlat, zlon, eps, x_m, y_m
    logical:: rotate
    real(kind=4) rlev

    logical::debug = .false.  
    data eps / 1.e-6_realkind /  

    real(kind=realkind) :: frland(klon,klat)  

    allocate (lat_m (klon, klat), lon_m (klon, klat), &
         zlat1_m (klon, klat), zlon1_m (klon, klat), &
         zlat2_m (klon, klat), zlon2_m (klon, klat) )
    allocate (i_m (klon, klat), j_m (klon, klat) )  

    allocate (w_x_m (klon, klat), w_y_m (klon, klat) )  
    iut = 6  
    if (mype==0.and.debug) then  
       write (6, * ) ' interpol_clim local grib ', klon, klat, &
            nx_grib, ny_grib

    endif
    !     calculate grid-point coordinates of non-staggered output field
    do j = 1, klat  
       zlat = south + real (jdatastart + j - 2,realkind) * dlat  
       do i = 1, klon  
          zlon = west + real (idatastart + i - 2,realkind) * dlon  
          lat_m (i, j) = zlat  
          lon_m (i, j) = zlon  
       enddo



    enddo
    !     rotate coordinates if needed
    if (abs (polon - polon_grib) >eps.or.abs (polat - polat_grib) >eps) then
       rotate = .true.  
       !     de-rotate if needed coordinates from output grid geometry
       if (abs (polat + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, lon_m, lat_m, klon, &
               klat, polon, polat, - 1)
       else  
          do j = 1, klat  
             do i = 1, klon  
                zlat1_m (i, j) = lat_m (i, j)  
                zlon1_m (i, j) = lon_m (i, j)  
             enddo
          enddo


       endif
       !     rotate if needed coordinates to input grid geometry
       if (abs (polat_grib + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, zlon2_m, zlat2_m, klon, &
               klat, polon_grib, polat_grib, &
               + 1)
       else  
          do j = 1, klat  
             do i = 1, klon  
                zlat2_m (i, j) = zlat1_m (i, j)  
                zlon2_m (i, j) = zlon1_m (i, j)  
             enddo
          enddo

       endif
       !     no rotation needed, just copy coordinates
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
          if(zlon2_m(i, j)<west_grib) then
             zlon2_m(i, j) = zlon2_m(i, j)+360.0_realkind
          endif
          if (zlon2_m(i, j) >=west_grib + 360.0_realkind) then  
             zlon2_m(i, j) = zlon2_m(i, j) - 360.0_realkind  
          endif
          y_m = (zlat2_m(i, j) - south_grib) / dlat_grib + 1.0_realkind  
          j_m(i, j) = int(y_m)  
          if(j_m(i,j)<1 .or. j_m(i,j)>ny_grib)then
             print *,j_m(i,j),ny_grib
             call stop_program('i out of bounds')
          endif

          w_y_m(i, j) = y_m - real(j_m(i, j),realkind )  
          x_m = (zlon2_m(i,j) - west_grib) / dlon_grib + 1.0_realkind  
          i_m(i,j) = int(x_m)  
          if(i_m(i,j)<1 .or. i_m(i,j)>nx_grib)then
             print *,i_m(i,j),nx_grib
             call stop_program( 'i out of bounds')
          endif
          w_x_m(i,j) = x_m - real (i_m(i,j),realkind )  
       enddo



    enddo
    !     read tsea,... from ecmwf climate files in lat/long
    if (mype==0.and.debug) then  
       write (6, * ) ' interpol_clim start read tsclim'  
    endif
    itype = 105  
    ipar = 011  
    ilev = 998  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), tsclim, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)

    itype = 105  
    ipar = 011  
    ilev = 999  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), tsd, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)
    !
    itype = 105  
    ipar = 098  
    ilev = 998  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), swdclim, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)
    !
    itype = 105  
    ipar = 098  
    ilev = 999  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), swd, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)

    itype = 105  
    ipar = 081  
    ilev = 000  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), frland, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)

    itype = 105  
    ipar = 011  
    ilev = 000  
    call grrdloc_clim(lundir, itype, ipar, real(ilev,realgribkind), tsm, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)

    itype = 105  
    ipar = 098  
    ilev = 000  
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), swm, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)
    where (swm<0._realkind) swm = 0.0_realkind
    swm = swm * 0.07_realkind

    if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
       tsea = tsm
    else
       call stop_program('We do not know what to do for other GCM with tsm climate')
    endif


    itype = 105 
    ipar = 66
    ilev = 0
    call grrdloc_clim (lundir, itype, ipar, real(ilev,realgribkind), snm, &
         klon, klat, nx_grib, ny_grib, i_m, w_x_m, j_m, w_y_m, &
         ierr)




    lwet_scal = .true.  
    r2 = 0.0_realkind  
    r3 = 0.0_realkind  
    k = 0  
    do j = 1, klat  
       do i = 1, klon  
          k = k + 1  
          if(tsea(i,j)<tseafr.and.frland(i,j)<0.30_realkind)then  
             frice(i,j) = 1.0_realkind  
          else  
             frice(i,j) = 0.0_realkind  
          endif
          if (lwet_scal) then  
             swd (i,j) = 0.07_realkind * (swd (i,j) + swdclim (i,j) ) * 0.5_realkind  
             swdclim (i,j) = 0.07_realkind * swdclim (i,j)  
          endif
          r2 = max (r2, swd (i,j) )  
          r3 = max (r3, swdclim (i,j) )  
       enddo
    enddo
    if (debug) then  
       write (6, * ) ' interpol_clim max(swd) ', r2, mype  
       write (6, * ) ' interpol_clim max(swdclim) ', r3, mype  
    endif
    deallocate (lat_m, lon_m, zlat1_m, zlon1_m, zlat2_m, zlon2_m)  
    deallocate (i_m, j_m, w_x_m, w_y_m)  
    return  
  end subroutine interpol_clim




  subroutine read_clim_year(year,  day, hour, klon, &
       klat, tsclim_year, west_in, south_in, dlon_in, dlat_in, &
       polon_in, polat_in,myddr)
    implicit none  
    integer , intent (in) ::klon, klat, year, day,  hour
    real(kind=realkind) :: west_in, south_in, dlon_in, dlat_in, polon_in, polat_in  
    real(kind=realkind) :: tsclim_year (klon, klat)  

    real(kind=realkind) , allocatable::tsclim (:, :)  
    integer  ::iy, im, id, ih, imin  
    integer :: ilen, length  
    real(kind=realkind), parameter::otwth = 1.0_realkind / 12.0_realkind  
    integer :: nxny, nx_global, ny_global  
    integer :: lundir  

    integer :: nx_grib, ny_grib  
    character (len=132) :: prefix, filenam  
    real(kind=realkind) :: west_grib2, south_grib2, dlon_grib2, dlat_grib2, &
         polon_grib2, polat_grib2
    integer, parameter::jpbuf = 700000  
    real(kind=realgribkind) :: gribbuf (jpbuf)  
    integer :: i, j  

    real(kind=realkind) :: west, south, dlamda, dtheta, aplon, aplat  
    type(ddr),intent(inout)::myddr

    allocate (tsclim (klon, klat) )  
    nx_global = klon_global  


    ny_global = klat_global  

    tsclim_year = 0.0_realkind  
    nx_grib = 360  
    ny_grib = 181  
    west_grib2 = - 180.0_realkind  
    south_grib2 = - 90.0_realkind  
    dlon_grib2 = 1.0_realkind  
    dlat_grib2 = 1.0_realkind  
    polon_grib2 = 0.0_realkind  
    polat_grib2 = - 90.0_realkind

    lundir = 61  
    west = west_in  
    south = south_in  
    dlamda = dlon_in  
    dtheta = dlat_in  
    aplon = polon_in  

    aplat = polat_in  
    !     Open climate file from ECMWF ( interpolated to finer grid )
    do im = 1, 12  
       call readgcpath(prefix, ilen)  
       prefix = prefix (1:ilen) //'clim'  
       ilen = ilen + 4  
       iy = year  
       id = 1!day  
       ih = 0!hour  
       imin = 0  
       length = 0  
       call cre_filnam (prefix,ilen,iy,im,id,ih,length,filenam)
       nxny = nx_grib * ny_grib  
       call groploc(lundir,filenam,gribbuf,nxny,myddr)
       nx_grib = myddr%nlonhl
       ny_grib = myddr%nlathl

       call interpolate_xtra2 (klon, klat, lundir, west, south, &
            dlamda, dtheta, aplon, aplat, tsclim, nx_grib, ny_grib, &
            west_grib2, south_grib2, dlat_grib2, dlon_grib2, polon_grib2, &
            polat_grib2)


       call grclloc (lundir, mype)  
       !     averaging over a whole year
       do i = 1, klon  
          do j = 1, klat  
             tsclim_year(i,j) = tsclim_year(i,j)+tsclim(i,j)*otwth  
          enddo
       enddo

    enddo
    deallocate (tsclim)  
    return  
  end subroutine read_clim_year


  subroutine read_clim_many_years(klon,klat,tsclim_years,myddr)
    implicit  none
    integer,intent(in)::klon,klat
    integer::year,day,hour
    real(kind=realkind),intent(inout)::tsclim_years(klon,klat)
    type(ddr),intent(inout)::myddr

    integer::ilen,iy,im,id,ih,imin,length
    integer,parameter::jpbuf=70000
    real(kind=realgribkind)::buf(jpbuf)
    integer::lenbuf,itype,ipar,ilev,nxny,ierr,i,j
    integer::lundir
    integer::nx_grib,ny_grib
    character(len=132)::prefix,filenam
    real(kind=realkind)::lat_m(klon,klat), lon_m(klon,klat), zlat1_m(klon,klat), &
         zlon1_m(klon,klat),zlat2_m(klon,klat), zlon2_m(klon,klat)
    integer::i_m(klon,klat), j_m(klon,klat)                                   
    real(kind=realkind)::w_x_m(klon,klat), w_y_m(klon,klat)
    real(kind=realkind)::eps                         
    real(kind=realkind)::x_m,y_m
    logical::rotate
    real(kind=realkind)::west_grib,south_grib,dlon_grib,dlat_grib,polon_grib,polat_grib
    real(kind=4) rlev

    eps=1.e-6_realkind

    nx_grib     = 360
    ny_grib     = 181
    west_grib   = -180.0_realkind
    south_grib  = -90.0_realkind
    dlon_grib   = 1.0_realkind
    dlat_grib   = 1.0_realkind
    polon_grib  = 0.0_realkind
    polat_grib  = -90.0_realkind

    do j=1,klat                                               
       lat_m(:,j) = RCAdomain%south + real( jdatastart + j - 2,realkind )*RCAdomain%dlat
    enddo
    do i=1,klon                                            
       lon_m(i,:) = RCAdomain%west + real( idatastart + i - 2,realkind )*RCAdomain%dlon
    enddo


    !    1.3   Rotate coordinates if needed                            
    if( abs(RCAdomain%polon-polon_grib)>eps .or.                        &
         abs(RCAdomain%polat-polat_grib)>eps        ) then               
       rotate = .true.                                            
       if( abs(RCAdomain%polat+90.0_realkind)>eps ) then                           
          call  regrot(zlon1_m,zlat1_m,lon_m,lat_m,              &
               klon,klat,&
               RCAdomain%polon,RCAdomain%polat,-1) 
       else                                                       
          zlat1_m = lat_m                         
          zlon1_m = lon_m                         
       endif
       if( abs(polat_grib+90.0_realkind)>eps ) then                      
          call  regrot(zlon1_m,zlat1_m,zlon2_m,zlat2_m,          &
               klon,klat,&
               polon_grib,polat_grib,+1)                          
       else                                                       
          zlat2_m = zlat1_m                       
          zlon2_m = zlon1_m                       
       endif
    else                                                          
       rotate = .false.                                           
       zlat2_m = lat_m                            
       zlon2_m = lon_m                            
    endif

    do j=1,klat                                                    
       do i=1,klon                                            
          if( zlon2_m(i,j)<west_grib )then
             zlon2_m(i,j)=zlon2_m(i,j)+360.0_realkind
          endif
          if( zlon2_m(i,j)>=west_grib+360.0_realkind )then
             zlon2_m(i,j)=zlon2_m(i,j)-360.0_realkind
          endif
          y_m = ( zlat2_m(i,j) - south_grib )/dlat_grib + 1.0_realkind     
          j_m(i,j)    = int( y_m )                                
          w_y_m(i,j)  = y_m - real(j_m(i,j),realkind)                     
          x_m = ( zlon2_m(i,j) - west_grib )/dlon_grib + 1.0_realkind      
          i_m(i,j)    = int( x_m )                                
          if ( i_m(i,j) > nx_grib ) then                       
             x_m = x_m - real(nx_grib,realkind)                           
             i_m(i,j) = i_m(i,j) - nx_grib                        
          endif
          w_x_m(i,j)  = x_m - real(i_m(i,j),realkind)                     
       enddo
    enddo

    lenbuf = jpbuf

    year=00
    day=00
    hour=00

    !     open climate file from ecmwf ( interpolated to finer grid )
    lundir=67
    call readclimyearspath(prefix,ilen)

    iy=1960
    im=01
    id=01
    ih=00
    imin=0
    length=0
    call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
    nxny=klon_global*klat_global
    call groploc(lundir,filenam,buf,lenbuf,myddr)
    nx_grib = myddr%nlonhl
    ny_grib = myddr%nlathl

    itype = 105                                                   
    ipar  = 011                                                   
    ilev  = 000                                                   
    call grrdloc_hint(lundir,itype,ipar,real(ilev,realgribkind),.false.,       &
         tsclim_years,tsclim_years,                    &
         klon,klat,&
         nx_grib,ny_grib,                                        &
         i_m,w_x_m,j_m,w_y_m,                                    &
         .false.,                                                &
         i_m,w_x_m,j_m,w_y_m,                                    &
         ierr)

    call grclloc(lundir,mype)
    return
  end subroutine read_clim_many_years

  subroutine readclimyearspath(climyears_path,ilen)

    implicit none
    character(len=132)::climyears_path
    integer:: ilen
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /climyears/climyears_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=climyears)
    close(61)
    if(mype==0)then
       write(6,nml=climyears)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=climyears)
       close(1)
    endif
    ilen=index(climyears_path,' ')-1
  end subroutine readclimyearspath



  subroutine read_soiltype(klon,klat,lfirst,soiltype,myddr)

    ! written by Ulf Hansson, Rossby Centre 020409 -
    !
    ! this is a small part of read_clim
    ! the part that reads soiltype,
    ! used only #ifdef ECHAM2
    ! the reason is that I don't want to call the 
    ! getclim-read_clim-move_clim-mix_clim
    ! for the ECHAM2 case
    ! as I do most of it in iNTERPOL_BD instead
    !
    implicit none
    integer:: klon,klat
    real(kind=realkind):: soiltype(klon,klat)
    logical:: lfirst

    integer:: jpbuf
    integer ::ilen,iy,im,id,ih,imin,length
    parameter (jpbuf=700000)
    real(kind=realgribkind):: buf(jpbuf)
    integer:: lenbuf,itype,ipar,ilev
    integer:: lundir
    integer:: nx_grib,ny_grib
    character(len=132):: prefix,filenam

    integer:: lstat
    logical:: lexist,lopen
    type(ddr),intent(inout)::myddr

    call stop_program( ' read_soiltype is errornous ')

    lenbuf = jpbuf


    !   1.0    Open climate file from HIRLAM

    if(mype==0)write(6,*)'read_soiltype called, lfirst= ', lfirst

    if (lfirst) then
       lundir=60
       if(mype==0) then
          inquire(unit=lundir, exist=lexist, opened=lopen,&
               iostat=lstat,name=filenam)
          if( lexist .and. lopen ) then
             close(lundir)
             inquire(unit=lundir, exist=lexist, opened=lopen,iostat=lstat)
          endif


          IY=0
          ! any month will do for soiltype            IM=month
          IM=1
          ID=0
          IH=0
          IMIN=0
          LENGTH=-1
          call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
       endif

       call groploc(lundir,filenam,buf,lenbuf,myddr)
       nx_grib = myddr%nlonhl
       ny_grib = myddr%nlathl

       if ( mype == 0 ) then
          inquire(UNIT=lundir, EXIST=Lexist, OPENED=Lopen,IOSTAT=Lstat)
       endif
       if ( nx_grib == 0 .or. ny_grib == 0 ) then
          nx_grib=klon_global
          ny_grib=klat_global
       endif



       itype = 105
       ipar  = 195
       ilev  = 000
       !       call grrdloc(lundir,itype,ipar,real(ilev,realgribkind),&
       !            soiltype(1,1),klon_global,klat_global,&
       !            klon,klat,nx_grib,ny_grib,ierr)


       !   1.6 close the orography field data base file

       call grclloc(lundir,mype)
    endif ! lfirst

    return
  end subroutine read_soiltype


  subroutine interpolate_xtra2 (nx, ny, lundir, west, south, dlon, &
       dlat, aplon, aplat, tsclim, nx_grib, ny_grib, west_grib, &
       south_grib, dlat_grib, dlon_grib, polon_grib, polat_grib)
    implicit none  
    integer :: nx, ny, lundir  
    real(kind=realkind) :: west, south, dlon, dlat, aplon, aplat  
    integer :: nx_grib, ny_grib  
    real(kind=realkind) :: west_grib, south_grib, dlat_grib, dlon_grib, polon_grib, &
         polat_grib


    real(kind=realkind) :: tsclim (nx * ny)  
    real(kind=realkind) :: lat_m (nx, ny), lon_m (nx, ny), zlat1_m (nx, ny), zlon1_m &
         (nx, ny), zlat2_m (nx, ny), zlon2_m (nx, ny)
    integer :: i_m (nx, ny), j_m (nx, ny)  
    real(kind=realkind) :: w_x_m (nx, ny), w_y_m (nx, ny)  
    real(kind=realkind) :: eps, x_m, y_m  
    logical :: rotate  
    integer :: i, j, ierr  

    data eps / 1.e-6_realkind /  
    !     Calculate grid-point coordinates of non-staggered output field
    do j = 1, ny  
       lat_m (:, j) = south + real (jdatastart + j - 2,realkind) * dlat  
    enddo
    do i = 1, nx  
       lon_m (i, :) = west + real (idatastart + i - 2,realkind) * dlon  


    enddo
    !     rotate coordinates if needed
    if (abs (aplon - polon_grib) >eps.or.abs (aplat - polat_grib) &
         >eps) then
       rotate = .true.  
       !     De-rotate if needed coordinates from output grid geometry
       if (abs (aplat + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, lon_m, lat_m,  nx, ny, &
               aplon, aplat, - 1)
       else  
          zlat1_m = lat_m  
          zlon1_m = lon_m  
       endif
       !     Rotate if needed coordinates to input grid geometry
       if (abs (polat_grib + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, zlon2_m, zlat2_m,  nx, &
               ny, polon_grib, polat_grib, + 1)
       else  
          zlat2_m = zlat1_m  
          zlon2_m = zlon1_m  
       endif
       !     No rotation needed, just copy coordinates
    else  
       rotate = .false.  
       zlat2_m = lat_m  
       zlon2_m = lon_m  


    endif
    !     calculate indeces and weights needed for horizontal interpolation
    do j = 1, ny  
       do i = 1, nx  
          if (zlon2_m (i, j) <west_grib) zlon2_m (i, j) = zlon2_m (i, j) &
               + 360.0_realkind
          if (zlon2_m (i, j) >=west_grib + 360._realkind) then  
             zlon2_m (i, j) = zlon2_m (i, j) - 360.0_realkind  

          endif
          y_m = (zlat2_m (i, j) - south_grib) / dlat_grib + 1.0_realkind  
          j_m (i, j) = int (y_m)  
          w_y_m (i, j) = y_m - real (j_m (i, j),realkind )  
          x_m = (zlon2_m (i, j) - west_grib) / dlon_grib + 1.0_realkind  
          i_m (i, j) = int (x_m)  
          w_x_m (i, j) = x_m - real (i_m (i, j),realkind )  
       enddo


    enddo
    call grrdloc_clim (lundir, 105, 11, 0.0_realgribkind, tsclim, nx, ny, nx_grib, &
         ny_grib, i_m, w_x_m, j_m, w_y_m, ierr)
    if (ierr/=0) then  
       if (mype==0) print * , 'Could not read TSCLIM'  
       call stop_program( 'interpolate_xtra2'  )


    endif
    return  
  end subroutine interpolate_xtra2


  subroutine interpolate_xtra (nx, ny, lundir, west, south, dlon, &
       dlat, aplon, aplat, ts, sw, sn, nx_grib, ny_grib, west_grib, &
       south_grib, dlat_grib, dlon_grib, polon_grib, polat_grib)
    implicit none  
    integer :: nx, ny, lundir  
    real(kind=realkind) :: west, south, dlon, dlat, aplon, aplat  
    integer :: nx_grib, ny_grib  
    real(kind=realkind) :: west_grib, south_grib, dlat_grib, dlon_grib, polon_grib, &
         polat_grib


    real(kind=realkind) :: ts (nx * ny), sw (nx * ny), sn (nx * ny)  
    real(kind=realkind) :: lat_m (nx, ny), lon_m (nx, ny), zlat1_m (nx, ny), zlon1_m &
         (nx, ny), zlat2_m (nx, ny), zlon2_m (nx, ny)
    integer :: i_m (nx, ny), j_m (nx, ny)  
    real(kind=realkind) :: w_x_m (nx, ny), w_y_m (nx, ny)  
    real(kind=realkind) :: eps, x_m, y_m  
    logical :: rotate  
    integer :: i, j, ierr  

    data eps / 1.e-6_realkind /  
    !     Calculate grid-point coordinates of non-staggered output field
    do j = 1, ny  
       lat_m (:, j) = south + real (jdatastart + j - 2,realkind) * dlat  
    enddo
    do i = 1, nx  
       lon_m (i, :) = west + real (idatastart + i - 2,realkind) * dlon  


    enddo
    !     rotate coordinates if needed
    if(abs(aplon-polon_grib)>eps &
         .or.abs(aplat - polat_grib)>eps) then
       rotate = .true.  
       !     De-rotate if needed coordinates from output grid geometry
       if (abs (aplat + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, lon_m, lat_m, nx, ny, &
               aplon, aplat, - 1)
       else  
          zlat1_m = lat_m  
          zlon1_m = lon_m  
       endif
       !     Rotate if needed coordinates to input grid geometry
       if (abs (polat_grib + 90.0_realkind) >eps) then  
          call regrot (zlon1_m, zlat1_m, zlon2_m, zlat2_m, nx, &
               ny, polon_grib, polat_grib, + 1)
       else  
          zlat2_m = zlat1_m  
          zlon2_m = zlon1_m  
       endif
       !     No rotation needed, just copy coordinates
    else  
       rotate = .false.  
       zlat2_m = lat_m  
       zlon2_m = lon_m  


    endif
    !     calculate indeces and weights needed for horizontal interpolation
    do j = 1, ny  
       do i = 1, nx  
          if (zlon2_m(i, j)<west_grib)then
             zlon2_m(i, j) = zlon2_m(i, j) + 360.0_realkind
          endif
          if (zlon2_m(i, j)>=west_grib + 360._realkind) then  
             zlon2_m(i, j) = zlon2_m(i, j) - 360.0_realkind  
          endif
          y_m = (zlat2_m(i, j) - south_grib) / dlat_grib + 1.0_realkind  
          j_m(i, j) = int(y_m)  
          w_y_m(i, j) = y_m - real (j_m(i, j),realkind )  
          x_m = (zlon2_m(i, j) - west_grib) / dlon_grib + 1.0_realkind  
          i_m(i, j) = int(x_m)  
          w_x_m(i, j) = x_m - real(i_m(i, j),realkind)  
       enddo


    enddo
    call grrdloc_clim (lundir, 105, 11, 0.0_realgribkind, ts, nx, ny, nx_grib, &
         ny_grib, i_m, w_x_m, j_m, w_y_m, ierr)
    if (ierr/=0) then  
       if (mype==0) print * , 'Could not read TS'  
       call stop_program( 'interpolate_xtra'  )

    endif
    call grrdloc_clim (lundir, 105, 98, 0.0_realgribkind, sw, nx, ny, nx_grib, &
         ny_grib, i_m, w_x_m, j_m, w_y_m, ierr)

    where (sw<0._realkind) sw = 0.0_realkind  
    !same scaling as in interpol_clim.F
    sw = sw * 0.07_realkind  
    if (ierr/=0) then  
       if (mype==0) print * , 'Could not read SW'  
       call stop_program( 'interpolate_xtra second place'  )

    endif
    call grrdloc_clim (lundir, 105, 66, 0.0_realgribkind, sn, nx, ny, nx_grib, &
         ny_grib, i_m, w_x_m, j_m, w_y_m, ierr)
    if (ierr/=0) then  
       if (mype==0) print * , 'Could not read SN '  
       call stop_program( 'interpolate xtra snow'  )

    endif
    return  
  end subroutine interpolate_xtra


  subroutine iniphy(klon,klat,intime,nsvarns,nlakevarns,nlaketypeIn,    &
       ptsm,pswm,psnm,ptsd,pswd, &
       ptss,ptsc, &
       soiltype,pswc,&
       pfrl,pfri,frac_t2,frac_t3,&
       palbicenl,psvarns,   &
       depthlake, plakevarns,                                 &
       west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)

    !**** *iniphy* - initiate physics
    !     *klon*      number of gridpoints in the x-direction
    !     *klat*      number of gridpoints in the y-direction
    !     *klev*      number of vertical levels
    !     *ptsm*      surface temperature
    !     *pswm*      soil water content
    !     *psnm*      snow depth
    !     *ptsd*      deep surface temperature
    !     *pswd*      deep soil water content
    !     *prou*      roughness parameter
    !     *ptss*      sea surface temperature
    !     *ptsc*      climatological surface temperature
    !     *pswc*      climatological soil water content



    !     *pfrl*      fraction of land
    !     *pfri*      fraction of ice
    !     *pfrf*      fraction of forest



    !     output parameters:



    !     fraction of land, pfrl, and fraction of ice, pfri, is checked
    !     so that the values are within the interval zero to one.

    !     surface temperature, ptsm, in defined as the temperature at
    !     the lowest model level, ptm(klev), over land, and as the sea
    !     surface temperature, ptss, over sea (i.e. pfrl=0.0).

    !    soil water content is checked so that the minimum value
    !     is not less than the value 0.0005.

    !     roughness parameter is checked so that the minimum value
    !     is not less than the value 0.00001.

    implicit none
    real(kind=realkind):: west_in,south_in,dlon_in,dlat_in,polon_in,polat_in
    integer :: klon , klat 
    type(time),intent(in)::intime
    integer:: nsvarns
    integer:: nlakevarns        
    integer:: nlaketypeIn         

    real(kind=realkind),intent(in)::ptsc(klon,klat)
    real(kind=realkind):: ptsm(klon,klat), pswm(klon,klat), psnm(klon,klat),           &
         ptsd(klon,klat), pswd(klon,klat), &
         ptss(klon,klat),  pswc(klon,klat),                            &
         pfrl(klon,klat), pfri(klon,klat), &
         pfrf(klon,klat),                             &
         soiltype(klon,klat), frac_t2(klon,klat), frac_t3(klon,klat)
    real(kind=realkind),allocatable::ptsclim_year(:,:)
    real(kind=realkind):: ptice
    !    real(kind=realkind) ws_echam(klon,klat), wsmx_echam(klon,klat), ptd5(klon,klat)
    real(kind=realkind) ::palbicenl(klon,klat),psvarns(klon,klat,nsvarns)
    real(kind=realkind), dimension(klon,klat,nlaketypeIn), intent(in) :: depthlake !kk0405 
    real(kind=realkind) , intent(out) :: plakevarns (klon,klat,nlaketypeIn,nlakevarns)

    real(kind=realkind):: zfrlim,zfrclr,zsnw,zfrfor,zfropl

    !     picethick   thickness of ice from oceanographic model
    !     
    !     psvarns are new surface variables, that shall be initialized
    !     where n=1,2,3,....18 and the time derivatives for the
    !     are stored in the model in dsvarnsdt(klon,klat,nsvarns)
    !     
    !     variable                                                  n
    !     
    !     l    tsns:      no snow temperature                            1
    !     l    tc:        canopy temperature                             2
    !     l    tsc:       ground temperature below trees                 3
    !     l    tsnow:     snow temperature over land                     4
    !     l    tssn:      soil temperature under snow                    5
    !     l    svegopl:   water on open land canopy (m)                  6 
    !     l    svegfor:   water on forest canopy (m)                     7
    !     l    snopl:     snow depth over land in water equivalent (m)   8
    !     l    snfor:     snow depth in forest in water equivalent (m)   9
    !     l    swsn:      water amount in snow over land                10
    !     l    rhosn:     current density of snow over land             11 
    !     l    snmax:     snowmax over open land                        12
    !     l    snmaxf:    snowmax for forest                            13
    !     
    !     l    tice:        ice surface temperature                     14
    !     l
    !     l    tsnc:      snow temperature below trees                  15
    !     l    tscsn:     soil temperature under snow in forest         16
    !     l    swsnc:     water amount in forest-snow                   17
    !     l    rhosnc:    current density of forest-snow                18
    !     sg021007
    !     l    tsns2:     second layer no snow temperature              19
    !     l    tsns3:     third layer no snow temperature               20
    !     l    tsns4:     fourth layer no snow temperature              21
    !     l    tsns5:     fifth layer no snow temperature               22
    !     l    tssn2:     second layer soil temp under snow             23
    !     l    tssn3:     third layer soil temp under snow              24
    !     l    tssn4:     fourth layer soil temp under snow             25
    !     l    tssn5:     fifth layer soil temp under snow              26
    !     l    tsc2:      second layer ground temperature below trees   27
    !     l    tsc3:      third layer ground temperature below trees    28
    !     l    tsc4:      fourth layer ground temperature below trees   29
    !     l    tsc5:      fifth layer ground temperature below trees    30
    !     l    tscsn2:    second layer soil temp under snow in forest   31
    !     l    tscsn3:    third layer soil temp under snow in forest    32
    !     l    tscsn4:    fourth layer soil temp under snow in forest   33
    !     l    tscsn5:    fifth layer soil temp under snow in forest    34
    !     sg021007
    !     l    snowcan:   intercepted snow on forest canopy             35
    !     ps040120
    !     l    ticed      deep ice temperature                          36
    !     l    ticesn     ice temperature under snow                    37
    !     l    ticesnd    deep ice temperature under snow               38
    !     l    tsnice     snow temperature on ice                       39
    !     l    snice      snow depth over ice in water equivalent (m)   40
    !     l    swsnice    water amount in snow over ice                 41
    !     l    rhosnice   current density of snow over ice              42
    !     l    snmaxice   snowmax over ice                              43
    !     ps040120
    !     l    sw         top soil water                                44
    !     l    swd        deep soil water                               45
    !     
    !     kk0405
    !     !      depthlake : depth of lake of each type, (m)
    !     plakevarns are lake variables, that shall be initialized
    !     where n=1,2,3,....12 and the time derivatives for the
    !     are stored in the model in dlakevarnsdt(klon,klat,nlakevarns,nlaketypeIn)
    !     variable                                                  n
    !     
    !     tlakesfc  :   lake surface temperature,                           1
    !     including ice and snow (pseudo-prognostic) (i,j)
    !     tlakesnow :  temperature of snow over the lake ice (k)            2
    !     tlakeice  :  temperature of lake ice (k)                          3
    !     tlakemean :  mean water lake temperature (k)                      4
    !     tlakeml   : mixed layer lake temperature (k)                      5
    !     tlakebot : lake bottom temperature (k)                            6
    !     ctlake   : chape-factor                                           7
    !     sndepthlake : snow depth over lake ice (m)                        8
    !     icedepthlake : lake ice depth (m)                                 9
    !     mldepth : mixed layer depth (m)                                  10
    !     t_b1lake : temperature on the extremum on bottom sedoments (k)   11
    !     : when bottom sediments = off, use dummy value
    !     h_b1lake: depth of the bottom sediments temperature extremum (m) 12
    !     : when bottom sediments = off, use dummy value
    integer:: jx,jy,ilan0,ilanf,ilan1,iice0,iicef,iice1,iprpe, i,j


    real(kind=realkind):: rhosnow_int
    integer::  jut,jice
    real(kind=realkind):: z1w,zslask,assat
    integer ::isoil

    integer:: ilaketype         ! lake type loop index
    type(ddr),intent(inout)::myddr
    integer::istart,istop,jstart,jstop
#ifdef MPI_SRC
#include"mpif.h"
    integer::ierr
    integer::ilan0loc,ilan1loc,ilanfloc,iice0loc,iice1loc,iicefloc
#endif

    if(mype==0)print *,' --------- in iniphy --------- '

    if(.not.lecham5)then
       !     create annual mean soil climate temperature field and
       !     a mask field for mountain areas to be used in surf_land
       allocate(ptsclim_year(klon,klat))
       call read_clim_year(intime%year,intime%day,intime%hour,klon,klat,ptsclim_year,&
            west_in,south_in,dlon_in,dlat_in,polon_in,polat_in,myddr)        
    endif

    call inisurf()

    z1w=0.072_realkind
    assat=0.02_realkind

    iprpe = 0

    ilan0 = 0
    ilanf = 0
    ilan1 = 0
    iice0 = 0
    iicef = 0
    iice1 = 0

    jut=6


    jice=0
    do jy = 1,klat
       do jx = 1,klon
          pfrl(jx,jy) = max( pfrl(jx,jy),0.0_realkind )
          pfrl(jx,jy) = min( pfrl(jx,jy),1.0_realkind )
          pfri(jx,jy) = max( pfri(jx,jy),0.0_realkind )
          pfri(jx,jy) = min( pfri(jx,jy),1.0_realkind )

          pfrf(jx,jy) = (frac_t2(jx,jy)+frac_t3(jx,jy))*pfrl(jx,jy)
          if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then

             ptice=psvarns(jx,jy,14)

             !     corr ptice for a few points near greenland
             if ( ptice < 0.01_realkind .and. pfri(jx,jy) > 0.01_realkind ) then
                ptice = ptss(jx,jy)
             endif
          endif
          if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
             ptsm(jx,jy)=ptsm(jx,jy)*pfrl(jx,jy)+(1.0_realkind-pfrl(jx,jy)) &
                  *((1.0_realkind-pfri(jx,jy))*ptss(jx,jy)+pfri(jx,jy)*ptice)
          else

             ptsm(jx,jy)=ptsm(jx,jy)*pfrl(jx,jy)+(1.0_realkind-pfrl(jx,jy)) &
                  *ptss(jx,jy)
          endif

          if ( ptsm(jx,jy) < 1.0_realkind ) then
             write(6,6011) mype,jx,jy,pfrl(jx,jy),ptss(jx,jy), &
                  pfri(jx,jy),ptice,ptsm(jx,jy)
             call stop_program( 'iniphy 305')
          endif
6011      format('iniphy skumt ',3i5,5f10.2)
          pswc(jx,jy) = max( pswc(jx,jy),0.0005_realkind  )

       enddo
    enddo

    call int_snow(intime%month,intime%day,rhosnow_int)

    !     new surface variables

    zfrlim=0.01_realkind
    zfrclr=0.20_realkind

    do j=1,klat
       do i=1,klon
          jice=1
          if(pfri(i,j)<0.01_realkind)jice=0
          !     no snow temperature is put equal to ts

          psvarns(i,j,1) = ptsm(i,j)

          !     canopy temperature is put equal to ts

          psvarns(i,j,2) = ptsm(i,j)

          !     ground temperature below trees is put equal to ts

          psvarns(i,j,3) = ptsm(i,j)

          !     snow temperature is put equal to ts, but not warmer than zero

          psvarns(i,j,4) = min(ptsm(i,j),tmelt)

          !     soil temperature under the snow is put equal to ts

          psvarns(i,j,5) = ptsm(i,j)

          !     both canopy waters are put to zero            

          psvarns(i,j,6) = 0.0_realkind          
          psvarns(i,j,7) = 0.0_realkind          

          !     snow depth over open land is given as zsnw*sn (sn=snmax if sn not = 0) 

          zfropl=0.0_realkind
          zfrfor = min(pfrl(i,j),pfrf(i,j))
          if(pfrl(i,j)>zfrlim) then
             zfrfor=zfrfor/pfrl(i,j) 
             zfropl=1.0_realkind-zfrfor
          endif

          !     frsn=1 since sn=snmax

          zsnw=zfropl+zfrfor*zfrclr

          psvarns(i,j,8) = zsnw*psnm(i,j)

          !     snow depth in the forest is given as (1-zsnw)*sn

          psvarns(i,j,9) = (1.0_realkind-zsnw)*psnm(i,j)

          !     water amount in snow over land is put equal to zero       

          psvarns(i,j,10) = 0.0_realkind                

          !     snow density over land is put equal to climatological value

          psvarns(i,j,11) = rhosnow_int

          !     snowmax for open land  is put equal to zero

          psvarns(i,j,12) = 0.0_realkind              

          !     snowmax for forest  is put equal to sn 

          psvarns(i,j,13) = psvarns(i,j,9)

          !     ice temperature is put equal to ts
          psvarns(i,j,14) = min(ptss(i,j),tseafr)

          if(.not.(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl))then
             psvarns(i,j,14) = real(jice,realkind)*psvarns(i,j,14) +  (1.0_realkind-real(jice,realkind))*tseafr
          endif

          !     forest-snow temperature is put equal to ts but not warmer than zero

          psvarns(i,j,15) = min(ptsm(i,j),tmelt)

          !     top-soil temperature under forest-snow is put equal to ts

          psvarns(i,j,16) = ptsm(i,j)

          !     water amount in forest-snow is put equal to zero

          psvarns(i,j,17) = 0.0_realkind

          !     forest-snow density is put equal to climatological value

          psvarns(i,j,18) = rhosnow_int

          !     put tsd in third layer for tsns3,tssn3,tsc3,tscsn3

          psvarns(i,j,20) = ptsd(i,j)
          psvarns(i,j,24) = ptsd(i,j)
          psvarns(i,j,28) = ptsd(i,j)
          psvarns(i,j,32) = ptsd(i,j)

          !     put tsdcli in fourth layer for tsns4,tssn4,tsc4,tscsn4

          psvarns(i,j,21) = ptsc(i,j)
          psvarns(i,j,25) = ptsc(i,j)
          psvarns(i,j,29) = ptsc(i,j)
          psvarns(i,j,33) = ptsc(i,j)

          !     put tsclim_year in fifth layer for tsns5,tssn5,tsc5,tscsn5

          if(lecham5)then
             call stop_program( 'solve this iniphy')
             !             psvarns(i,j,22) = ptd5(i,j)
             !             psvarns(i,j,26) = ptd5(i,j)
             !             psvarns(i,j,30) = ptd5(i,j)
             !             psvarns(i,j,34) = ptd5(i,j)
          else
             psvarns(i,j,22) = ptsclim_year(i,j)
             psvarns(i,j,26) = ptsclim_year(i,j)
             psvarns(i,j,30) = ptsclim_year(i,j)
             psvarns(i,j,34) = ptsclim_year(i,j)
          endif

          !     second layer ground temperatures, tsns2,tssn2,tsc2,tscsn2, 
          !     interpolated between first and third layers

          psvarns(i,j,19) = 0.9_realkind*ptsm(i,j)+0.1_realkind*ptsd(i,j)
          psvarns(i,j,23) = 0.9_realkind*ptsm(i,j)+0.1_realkind*ptsd(i,j)
          psvarns(i,j,27) = 0.9_realkind*ptsm(i,j)+0.1_realkind*ptsd(i,j)
          psvarns(i,j,31) = 0.9_realkind*ptsm(i,j)+0.1_realkind*ptsd(i,j)

          !     intercepted snow on forest canopy is put equal to zero

          psvarns(i,j,35) = 0.0_realkind

          !     ticed      deep ice temperature
          psvarns(i,j,36) = 0.5_realkind*(psvarns(i,j,14)+tseafr)
          if(.not.(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl))then
             psvarns(i,j,36) = real(jice,realkind)*psvarns(i,j,36) +(1.0_realkind-real(jice,realkind))*tseafr
          endif

          !     ticesn     ice temperature under snow

          psvarns(i,j,37) = psvarns(i,j,36)

          !     ticesnd    deep ice temperature under snow

          psvarns(i,j,38) = psvarns(i,j,36)

          !     tsnice     snow temperature on ice

          psvarns(i,j,39) = psvarns(i,j,14)

          !     snice      snow depth over ice in water equivalent (m)

          psvarns(i,j,40) = real(jice,realkind)*psnm(i,j)

          !     swsnice    water amount in snow over ice

          psvarns(i,j,41) = 0.0_realkind

          !     rhosnice   current density of snow over ice

          psvarns(i,j,42) = real(jice,realkind)*rhosnow_int + (1.0_realkind-real(jice,realkind))*100.0_realkind

          !     snmaxice   snowmax over ice

          psvarns(i,j,43) = 0.0_realkind

          !     sw   top soil water

          isoil=nint(soiltype(i,j))
          if(isoil< 1 .or. isoil> 12) then
             isoil=6
          endif

          if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
             zslask=pswm(i,j)/assat
          elseif(lecham2)then
             call stop_program( 'solve iniphy!!' )
             !             zslask=ws_echam(i,j)/wsmx_echam(i,j)
          else
             zslask=0.0_realkind
          endif
          zslask=max(min(zslask,1.0_realkind),0.0_realkind)
          zslask=zslask*(vcc(isoil)-vfl(isoil))+vfl(isoil)
          psvarns(i,j,44)=zslask*z1w
          psvarns(i,j,45)=zslask*z1w

          !     swd   deep soil water

          if(lecmwf.or.lhadley.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
             zslask=pswd(i,j)/assat
          endif
          if(lecham2)then
             call stop_program( 'solve iniphy')
             !             zslask=ws_echam(i,j)/wsmx_echam(i,j)
          endif
          zslask=max(min(zslask,1.0_realkind),0.0_realkind)
          zslask=vfl(isoil)+(vcc(isoil)-vfl(isoil))*zslask
          psvarns(i,j,46)=zslask*z1w
          psvarns(i,j,47)=zslask*z1w
          psvarns(i,j,48)=zslask*z1w
          psvarns(i,j,49)=zslask*z1w
          if( pfri(i,j) > 0.0_realkind) then
             !     au010412            palbicenl(i,j) = albice
             !     value from inisurf.f       albice not set yet !!!
             palbicenl(i,j) = 0.20_realkind
          else
             palbicenl(i,j) = 0.0_realkind
          endif

          do ilaketype= 1, nlaketypeIn
             !     lake surface temperature (pseudo-prognostic) put to deep soil temperature
             plakevarns(i,j,ilaketype,1) = ptsd(i,j)
             call lakeinit(plakevarns(i,j,ilaketype,1), &
                  depthlake(i,j,ilaketype),              &
                  plakevarns(i,j,ilaketype,2),           &
                  plakevarns(i,j,ilaketype,3),           &
                  plakevarns(i,j,ilaketype,4),           &
                  plakevarns(i,j,ilaketype,5),           &
                  plakevarns(i,j,ilaketype,6),           &
                  plakevarns(i,j,ilaketype,7),           &
                  plakevarns(i,j,ilaketype,8),           &
                  plakevarns(i,j,ilaketype,9),           &
                  plakevarns(i,j,ilaketype,10),          &
                  plakevarns(i,j,ilaketype,11),          &
                  plakevarns(i,j,ilaketype,12))
          enddo
       enddo
    enddo

    if(.not.lecham5)then
       deallocate(ptsclim_year)
    endif

    call jlimits(klon,klat,istart,istop,jstart,jstop)
    do jy = jstart,jstop       
       do jx = istart,istop
          if(abs(pfrl(jx,jy))<1.e-14_realkind) ilan0 = ilan0 + 1
          if(abs(pfrl(jx,jy)-1.0_realkind)<1.e-14_realkind) ilan1 = ilan1 + 1
          if(pfrl(jx,jy)>0.0_realkind .and. pfrl(jx,jy)<1.0_realkind) ilanf = ilanf + 1
          if(abs(pfri(jx,jy))<1.e-14_realkind) iice0 = iice0 + 1
          if(abs(pfri(jx,jy)-1.0_realkind)<1.e-14_realkind) iice1 = iice1 + 1
          if(pfri(jx,jy)>0.0_realkind .and. pfri(jx,jy)<1.0_realkind) iicef = iicef + 1
       enddo
    enddo
#ifdef MPI_SRC
    ilan0loc=ilan0
    ilan1loc=ilan1
    ilanfloc=ilanf
    iice0loc=iice0
    iice1loc=iice1
    iicefloc=iicef
    call MPI_reduce(ilan0loc,ilan0,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
    call MPI_reduce(ilan1loc,ilan1,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
    call MPI_reduce(ilanfloc,ilanf,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
    call MPI_reduce(iice0loc,iice0,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
    call MPI_reduce(iice1loc,iice1,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
    call MPI_reduce(iicefloc,iicef,1,MPI_INTEGER,MPI_SUM,iprpe,LOCALCOMM,ierr)
#endif
    if(mype == iprpe) then
       write(6,'(1x,'' '')')
       write(6,'(1x,''numb. of points with pure   land = '',i9)')ilan1
       write(6,'(1x,''numb. of points with no     land = '',i9)')ilan0
       write(6,'(1x,''numb. of points with partial land= '',i9)')ilanf
       write(6,'(1x,''numb. of points with pure    ice = '',i9)')iice1
       write(6,'(1x,''numb. of points with no      ice = '',i9)')iice0
       write(6,'(1x,''numb. of points with partial ice = '',i9)')iicef
       write(6,'(1x,'' '')')
       write(6,'(1x,"mean-, min- and max-values of surface fields:")')
    endif


    call minmax ('ts      ',ptsm,klon,klat)
    call minmax ('ts-sea  ',ptss,klon,klat)
    call minmax ('ts-deep ',ptsd,klon,klat)
    call minmax ('ts-clim ',ptsc,klon,klat)
    call minmax ('sw      ',pswm,klon,klat)
    call minmax ('sw-deep ',pswd,klon,klat)
    call minmax ('sw-clim ',pswc,klon,klat)
    call minmax ('sn      ',psnm,klon,klat)
    !    call minmax ('sn-clim ',psnc,klon,klat)
    !    call minmax ('z0      ',prou,klon,klat)
    !    call minmax ('z0-clim ',proc,klon,klat)
    !    call minmax ('albedo  ',palb,klon,klat)
    call minmax ('fr-land ',pfrl,klon,klat)
    call minmax ('fr-ice  ',pfri,klon,klat)

    call minmax ('fr-ice  ',pfri,klon,klat)
    if(mype==0)print*,' --------- iniphy done --------- '

    return
  end subroutine iniphy


  ! ++++++++++++++++++++++++++++++++++++++++++++++++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine minmax(message, pa, klon, klat)  
    !     minmax - statistics computations
    !     compute and print
    !
    !     mean-value
    !     minimum -value and position
    !     maximum -value and position
    !
    !   interface.
    !     ----------
    !
    !    minmax is an auxiliary routine
    !
    !     input parameters:
    !     -----------------
    !
    !     message        message to be printed
    !     pa           local array array with data
    !     klon         local  number of gridpoints in x-direction
    !     klat         local  number of gridpoints in y-direction
    !     externals.



    implicit none  
    character (len=*) :: message  
    integer :: klon, klat  
    real(kind=realkind) :: pa (klon, klat)  
    !    real(kind=realkind) :: pa_g(klon_global*klat_global)  
    !    integer :: imin, imax, iymin, ixmin, iymax, ixmax  
    real(kind=realkind) :: zsum, zmin, zmax  
    integer::istart,istop,jstart,jstop
#ifdef MPI_SRC
#include"mpif.h"
    integer::ierr
    real(kind=realkind)::zminloc,zmaxloc,zsumloc
#endif
    call jlimits(klon,klat,istart,istop,jstart,jstop)
    zmin = minval(pa)
    zmax = maxval(pa)
    zsum = sum(pa(istart:istop,jstart:jstop))
#ifdef MPI_SRC
    zminloc=zmin
    zmaxloc=zmax
    zsumloc=zsum
    call MPI_reduce(zminloc,zmin,1,REALTYPE,MPI_MIN,0,LOCALCOMM,ierr)
    call MPI_reduce(zmaxloc,zmax,1,REALTYPE,MPI_MAX,0,LOCALCOMM,ierr)
    call MPI_reduce(zsumloc,zsum,1,REALTYPE,MPI_SUM,0,LOCALCOMM,ierr)
#endif    
    zsum = zsum/real(klon_global*klat_global,realkind)

    if(mype==0)then
       write (6, '(1x,/,1x,(a))') message  
       write (6, '(1x,''mean-value: '',g20.8)') zsum  
       write (6, '(1x,''min -value: '',g20.8)') zmin 
       write (6, '(1x,''max -value: '',g20.8)') zmax
    endif
    return  
  end subroutine minmax



end module climate

