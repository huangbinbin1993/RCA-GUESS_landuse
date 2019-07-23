Module util
  use decomp,only:stop_program,realkind
  implicit none
  private
  public adddtg,regrot,turnwi1,turnwireg2rot,turnwirot2reg,meanmnmx
contains
  subroutine turnwi1(uarg,varg,ures,vres,a,b,c,d,klon,klat,klev)
    implicit none

    integer,intent(in):: klon,klat,klev
    real(kind=realkind),intent(in):: uarg(klon,klat,klev),varg(klon,klat,klev)
    real(kind=realkind),intent(out):: ures(klon,klat,klev),vres(klon,klat,klev)
    real(kind=realkind),intent(in):: a(klon,klat),b(klon,klat),c(klon,klat),d(klon,klat)

    integer:: jy,jx,jz                                              

    !         multiplication between regular and rotated spherical grid
    do jy = 1,klat
       do jx = 1,klon
          do jz = 1,klev
             ures(jx,jy,jz)=a(jx,jy)*uarg(jx,jy,jz)+b(jx,jy)*varg(jx,jy,jz)
             vres(jx,jy,jz)=c(jx,jy)*uarg(jx,jy,jz)+d(jx,jy)*varg(jx,jy,jz)
          enddo
       enddo
    enddo
    
    return
  end subroutine turnwi1




  subroutine turnwireg2rot(a,b,c,d,&
       xreg,yreg,xrot,yrot,klon,klat,        &
       xcen,ycen,kcall)

    !         turn horizontal velocity components between regular and
    !         rotated spherical coordinates.
    implicit none

    integer,intent(in):: klon,klat
    integer,intent(in):: kcall 
    real(kind=realkind),intent(out),dimension(klon,klat)::a,b,c,d
    real(kind=realkind),intent(in),dimension(klon,klat)::xreg,yreg,xrot,yrot
    real(kind=realkind),intent(in):: xcen,ycen                                                 

    integer:: jy,jx                                              
    real(kind=realkind):: pi,zrad,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,  &
         zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot

    if (kcall>0) then! precalculate matrixes from regular to rotated grid
       pi = 4._realkind*atan(1._realkind)
       zrad = pi/180._realkind
       zsyc = sin(zrad*(ycen+90._realkind))
       zcyc = cos(zrad*(ycen+90._realkind))

       do jy = 1,klat
          do jx = 1,klon
             zsxreg = sin(zrad*xreg(jx,jy))
             zcxreg = cos(zrad*xreg(jx,jy))
             zsyreg = sin(zrad*yreg(jx,jy))
             zcyreg = cos(zrad*yreg(jx,jy))

             zxmxc  = zrad*(xreg(jx,jy) - xcen)
             zsxmxc = sin(zxmxc)
             zcxmxc = cos(zxmxc)

             zsxrot = sin(zrad*xrot(jx,jy))
             zcxrot = cos(zrad*xrot(jx,jy))
             zsyrot = sin(zrad*yrot(jx,jy))
             zcyrot = cos(zrad*yrot(jx,jy))

             a(jx,jy) = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
             b(jx,jy) = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - zsxmxc*zsyreg*zcxrot
             c(jx,jy) = zsyc*zsxmxc/zcyrot
             d(jx,jy) = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
          enddo
       enddo
!!$    elseif (kcall<0) then! precalculate matrixes from rotated to regular grid
!!$       pi = 4._realkind*atan(1._realkind)
!!$       zrad = pi/180._realkind
!!$       zsyc = sin(zrad*(pycen+90._realkind))
!!$       zcyc = cos(zrad*(pycen+90._realkind))
!!$
!!$       do jy = 1,klat
!!$          do jx = 1,klon
!!$             zsxreg = sin(zrad*pxreg(jx,jy))
!!$             zcxreg = cos(zrad*pxreg(jx,jy))
!!$             zsyreg = sin(zrad*pyreg(jx,jy))
!!$             zcyreg = cos(zrad*pyreg(jx,jy))
!!$
!!$             zxmxc  = zrad*(pxreg(jx,jy) - pxcen)
!!$             zsxmxc = sin(zxmxc)
!!$             zcxmxc = cos(zxmxc)
!!$
!!$             zsxrot = sin(zrad*pxrot(jx,jy))
!!$             zcxrot = cos(zrad*pxrot(jx,jy))
!!$             zsyrot = sin(zrad*pyrot(jx,jy))
!!$             zcyrot = cos(zrad*pyrot(jx,jy))
!!$
!!$             pa(jx,jy) = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
!!$             pb(jx,jy) = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -  zcxmxc*zsxrot*zsyrot
!!$             pc(jx,jy) = -zsyc*zsxrot/zcyreg
!!$             pd(jx,jy) = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
!!$          enddo
!!$       enddo
    endif

    return
  end subroutine turnwireg2rot


  subroutine turnwirot2reg(a,b,c,d,&
       xreg,yreg,xrot,yrot,klon,klat,        &
       xcen,ycen,kcall)

    !         turn horizontal velocity components between regular and
    !         rotated spherical coordinates.
    !         puarg(klon,klat,klev) : input u components
    !         pvarg(klon,klat,klev) : input v components
    !         pures(klon,klat,klev) : output u components
    !         pvres(klon,klat,klev) : output v components
    !         pa(klon,klat,klev)    : transformation coefficients
    !         pb(klon,klat,klev)    :    -"-
    !         pc(klon,klat,klev)    :    -"-
    !         pd(klon,klat,klev)    :    -"-
    !         pxreg(klon,klat,klev) : regular longitudes
    !         pyreg(klon,klat,klev) : regular latitudes
    !         pxrot(klon,klat,klev) : rotated longitudes
    !         pyrot(klon,klat,klev) : rotated latitudes
    !         klon              : dimension in the x (longitude) direction
    !         klat              : dimension in the y (latitude) direction
    !         klev               : number of vertical layers
    !         kx                 : number of gridpoints in the x direction
    !         ky                 : number of gridpoints in the y direction
    !         pxcen              : regular longitude of the south pole of the
    !                              transformed grid
    !         pycen              : regular latitude of the south pole of the
    !                              transformed grid
    !         kcall=-2 or 2      : preparation of coefficients.
    !         kcall=-1 or 1      : multiplication after preparations.
    !     
    !         kcall < 0          : find wind components in regular coordinates
    !                              from wind components in rotated coordinates
    !         kcall > 0          : find wind components in rotated coordinates
    !                              from wind components in regular coordinates
    !         Note that all coordinates are given in degrees N and degrees E.
    !            (negative values for S and W) 
    !     
    !         j.e. haugen   hirlam   june -92
    !         altered slightly by magnus nov -96   

    implicit none


    integer,intent(in):: klon,klat
    integer,intent(in):: kcall 
    real(kind=realkind),intent(out):: a(klon,klat),   b(klon,klat),  c(klon,klat),   d(klon,klat)
    real(kind=realkind),intent(in)::xreg(klon,klat),yreg(klon,klat), xrot(klon,klat),yrot(klon,klat)                  
    real(kind=realkind),intent(in):: xcen,ycen                                                 

    integer:: jy,jx                                              
    real(kind=realkind):: pi,zrad,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,  &
         zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot

!!$    if (kcall>0) then! precalculate matrixes from regular to rotated grid
!!$       
!!$       pi = 4._realkind*atan(1._realkind)
!!$       zrad = pi/180._realkind
!!$       zsyc = sin(zrad*(pycen+90._realkind))
!!$       zcyc = cos(zrad*(pycen+90._realkind))
!!$
!!$       do jy = 1,klat
!!$          do jx = 1,klon
!!$             zsxreg = sin(zrad*pxreg(jx,jy))
!!$             zcxreg = cos(zrad*pxreg(jx,jy))
!!$             zsyreg = sin(zrad*pyreg(jx,jy))
!!$             zcyreg = cos(zrad*pyreg(jx,jy))
!!$
!!$             zxmxc  = zrad*(pxreg(jx,jy) - pxcen)
!!$             zsxmxc = sin(zxmxc)
!!$             zcxmxc = cos(zxmxc)
!!$
!!$             zsxrot = sin(zrad*pxrot(jx,jy))
!!$             zcxrot = cos(zrad*pxrot(jx,jy))
!!$             zsyrot = sin(zrad*pyrot(jx,jy))
!!$             zcyrot = cos(zrad*pyrot(jx,jy))
!!$
!!$             pa(jx,jy) = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
!!$             pb(jx,jy) = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - zsxmxc*zsyreg*zcxrot
!!$             pc(jx,jy) = zsyc*zsxmxc/zcyrot
!!$             pd(jx,jy) = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
!!$          enddo
!!$       enddo
!!    else
    if (kcall<0) then! precalculate matrixes from rotated to regular grid
       pi = 4._realkind*atan(1._realkind)
       zrad = pi/180._realkind
       zsyc = sin(zrad*(ycen+90._realkind))
       zcyc = cos(zrad*(ycen+90._realkind))

       do jy = 1,klat
          do jx = 1,klon
             zsxreg = sin(zrad*xreg(jx,jy))
             zcxreg = cos(zrad*xreg(jx,jy))
             zsyreg = sin(zrad*yreg(jx,jy))
             zcyreg = cos(zrad*yreg(jx,jy))

             zxmxc  = zrad*(xreg(jx,jy) - xcen)
             zsxmxc = sin(zxmxc)
             zcxmxc = cos(zxmxc)

             zsxrot = sin(zrad*xrot(jx,jy))
             zcxrot = cos(zrad*xrot(jx,jy))
             zsyrot = sin(zrad*yrot(jx,jy))
             zcyrot = cos(zrad*yrot(jx,jy))

             a(jx,jy) = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
             b(jx,jy) = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -  zcxmxc*zsxrot*zsyrot
             c(jx,jy) = -zsyc*zsxrot/zcyreg
             d(jx,jy) = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
          enddo
       enddo
    endif

    return
  end subroutine turnwirot2reg

  subroutine regrot(xreg,yreg,xrot,yrot,klon,klat, &
       xcen,ycen,kcall)

    !     conversion between regular and rotated spherical coordinates.
    !     
    !     pxreg     longitudes of the regular coordinates
    !     pyreg     latitudes of the regular coordinates
    !     xrot     longitudes of the rotated coordinates
    !     yrot     latitudes of the rotated coordinates
    !     all coordinates given in degrees n (negative for s)
    !     and degrees e (negative values for w)
    !     klon     dimension of the gridpoint fields in the x-direction
    !     klat     dimension of the gridpoint fields in the y-direction
    !     kx        number of gridpoint in the x-direction
    !     ky        number of gridpoints in the y-direction
    !     xcen     regular longitude of the south ole of the rotated grid
    !     ycen     regular latitude of the south pole of the rotated grid
    !     
    !     kcall=-1: find regular as functions of rotated coordinates.
    !     kcall= 1: find rotated as functions of regular coordinates.
    !     
    !     j.e. haugen   hirlam   june -92
    implicit none

    integer,intent(in):: klon,klat,kcall
    real(kind=realkind):: xreg(klon,klat),yreg(klon,klat),&
         xrot(klon,klat),yrot(klon,klat)
    real(kind=realkind),intent(in)::xcen,ycen

    real(kind=realkind)::pi,zrad,zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg,&
         zsyrot,zcyrot,zcxrot,zsxrot,zradi
    integer jy,jx
    pi = 4._realkind*atan(1._realkind)
    zrad = pi/180._realkind
    zradi = 1._realkind/zrad
    zsycen = sin(zrad*(ycen+90._realkind))
    zcycen = cos(zrad*(ycen+90._realkind))

    if (kcall==1) then
       do jy = 1,klat
          do jx = 1,klon
             zxmxc  = zrad*(xreg(jx,jy) - xcen)
             zsxmxc = sin(zxmxc)
             zcxmxc = cos(zxmxc)
             zsyreg = sin(zrad*yreg(jx,jy))
             zcyreg = cos(zrad*yreg(jx,jy))
             zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
             zsyrot = max(zsyrot,-1.0_realkind)
             zsyrot = min(zsyrot,+1.0_realkind)
             yrot(jx,jy) = asin(zsyrot)*zradi
             zcyrot = cos(yrot(jx,jy)*zrad)
             zcxrot = (zcycen*zcyreg*zcxmxc + zsycen*zsyreg)/zcyrot
             zcxrot = max(zcxrot,-1.0_realkind)
             zcxrot = min(zcxrot,+1.0_realkind)
             zsxrot = zcyreg*zsxmxc/zcyrot
             xrot(jx,jy) = acos(zcxrot)*zradi
             if (zsxrot<0.0_realkind) xrot(jx,jy) = -xrot(jx,jy)
          enddo
       enddo
    elseif (kcall==-1) then
       do jy = 1,klat
          do jx = 1,klon
             zsxrot = sin(zrad*xrot(jx,jy))
             zcxrot = cos(zrad*xrot(jx,jy))
             zsyrot = sin(zrad*yrot(jx,jy))
             zcyrot = cos(zrad*yrot(jx,jy))
             zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
             zsyreg = max(zsyreg,-1.0_realkind)
             zsyreg = min(zsyreg,+1.0_realkind)
             yreg(jx,jy) = asin(zsyreg)*zradi
             zcyreg = cos(yreg(jx,jy)*zrad)
             zcxmxc = (zcycen*zcyrot*zcxrot -  zsycen*zsyrot)/zcyreg
             zcxmxc = max(zcxmxc,-1.0_realkind)
             zcxmxc = min(zcxmxc,+1.0_realkind)
             zsxmxc = zcyrot*zsxrot/zcyreg
             zxmxc  = acos(zcxmxc)*zradi
             if (zsxmxc<0.0_realkind) zxmxc = -zxmxc
             xreg(jx,jy) = zxmxc + xcen
          enddo
       enddo
    else
       write(6,'(1x,''invalid kcall in regrot'')')
       call stop_program('invalid kcall in regrot')
    endif

    return
  end subroutine regrot


  subroutine adddtg(kdi,kti,kfp,kd,kt)

    !     returns new date yyyymmdd and time hhmmss by adding kfp (seconds,
    !     may be negative) to old values

    implicit none
    integer,intent(in)::kdi
    integer,intent(in)::kti
    integer,intent(in)::kfp
    integer,intent(out)::kd
    integer,intent(out)::kt

    integer,parameter::ispd=86400 !seconds/day
    integer::ic
    integer::id
    integer::is
    integer::ih
    integer::ir
    integer::im

    !     convert to century date
    ic=idat2c(kdi)

    !     split kfp into days and seconds
    id=kfp/ispd
    if(kfp<0)id=id-1
    is=kfp-id*ispd

    !     decode kti into seconds, and add to seconds from kti
    ih=kti/10000
    ir=kti-ih*10000
    im=ir/100
    is=is + ir-im*100 + im*60 + ih*3600

    !     split seconds into days and seconds
    ir=is/ispd
    id=id+ir
    is=is-ir*ispd

    !     add days, convert to yyyymmdd format
    ic=ic+id
    kd=ic2dat(ic)

    !     convert seconds to hhmmss format
    ih=is/3600
    ir=is-ih*3600
    im=ir/60
    is=ir-im*60
    kt=ih*10000 + im*100 + is
  end subroutine adddtg


  subroutine gregor(jd,iy,im,id)

    !     converts julian daynumber into gregorian ( normal ) date
    !     ( year,month,day ) .
    implicit none
    integer  jd, iy, im, id
    integer  l, n

    l  = jd + 68569 + 2415020
    n  = 4*l / 146097
    l  = l - ( 146097*n + 3 ) / 4
    iy = 4000 * ( l+1 ) / 1461001
    l  = l - 1461 * iy / 4 + 31
    im = 80 * l / 2447
    id = l - 2447 * im / 80
    l  = im / 11
    im = im + 2 - 12 * l
    iy = 100 * ( n- 49 ) + iy + l

    return
  end subroutine gregor
  integer function ic2dat(kd)

    !     converts julian daynumber since 19000101 into gregorian yyyymmdd
    implicit none
    integer  kd
    integer  iy, im, id, l, n

    l  = kd + 68569 + 2415020
    n  = 4*l / 146097
    l  = l - ( 146097*n + 3 ) / 4
    iy = 4000 * ( l+1 ) / 1461001
    l  = l - 1461 * iy / 4 + 31
    im = 80 * l / 2447
    id = l - 2447 * im / 80
    l  = im / 11
    im = im + 2 - 12 * l
    iy = 100 * ( n- 49 ) + iy + l
    ic2dat=10000*iy+im*100+id
  end function ic2dat

  integer function ic2ymd ( ic )

    !     conversion of date representation.

    !     input:   ic - integer - day in centuary (1=01/01/1900)

    !     output:       integer - yymmdd
    !     yy = year  (00..99)
    !     mm = month (01..12)
    !     dd = day   (01..31)

    implicit none
    integer  ic, ky, km, kd

    call gregor  ( ic, ky, km, kd )

    ky = mod ( ky , 100 )
    ic2ymd = ky*10000 + km*100 + kd

    return

  end function ic2ymd
  integer function idat2c(kd)

    !     for given kd(yyyymmdd) return day number since 19000101
    implicit none
    integer kd,iy,ir,im,id
    iy=kd/10000
    ir=kd-iy*10000
    if(iy<100)iy=iy+1900
    im=ir/100+1
    id=ir-im*100

    if(im<=3)then
       iy=iy-1
       im=im+12
    endif

    ir=iy/100

    idat2c=365*iy-693923+iy/4-ir+ir/4+int(30.6001_realkind*real(im,realkind))+id

  end function idat2c
  integer function idat2wd(kd)
    !     for given kd(yyyymmdd) return weekday (0=Sunday)
    implicit none
    integer kd
    idat2wd=mod(idat2c(kd)+2483593,7)
  end function idat2wd
  integer function iy2eastr(ky)
    !     for given year(yyyy) return date of easter (mmdd)
    implicit none
    integer ky,ia,ib,ic,id,ie,iff,ig,ih,ii,ik,il,im,in

    ia=mod(ky,19)
    ib=ky/100
    ic=ky-ib*100
    id=ib/4
    ie=ib-id*4
    iff=(ib+8)/25
    ig=(ib-iff+1)/3
    ih=mod(19*ia+ib-id-ig+15,30)
    ii=ic/4
    ik=ic-ii*4
    il=mod(32+2*(ie+ii)-ih-ik,7)
    im=(ia+11*ih+22*il)/451
    in=(ih+il-7*im+114)/31
    iy2eastr=in*100+(ih+il-7*im+114-in*31)+1
  end function iy2eastr
  integer function iymd2c(kymd)
    !        this function returns the century day given an integer yyyymmdd
    !        counting from 1-1-1900
    implicit none
    integer kymd
    iymd2c=idat2c(kymd)
  end function iymd2c
  logical function lcmpdtg(kd1,kt1,yop,kd2,kt2)
    implicit none
    integer:: kd1, kt1, kd2, kt2
    character(len=2):: yop
    integer:: iop
    iop=0

    if(yop(1:1)=='g'.or.yop(1:1)=='G')then
       if(yop(2:2)=='t'.or.yop(2:2)=='T')iop=1  !gt
       if(yop(2:2)=='e'.or.yop(2:2)=='E')iop=2  !ge
    endif
    if(yop(1:1)=='e'.or.yop(1:1)=='E')then
       if(yop(2:2)=='q'.or.yop(2:2)=='Q')iop=3  !eq
    endif
    if(yop(1:1)=='l'.or.yop(1:1)=='L')then
       if(yop(2:2)=='e'.or.yop(2:2)=='E')iop=4  !le
       if(yop(2:2)=='t'.or.yop(2:2)=='T')iop=5  !lt
    endif
    if(yop(1:1)=='n'.or.yop(1:1)=='N')then
       if(yop(2:2)=='e'.or.yop(2:2)=='E')iop=6  !ne
    endif
    if(iop==0)then
       write(0,*)'lcmpdtg called with illegal operand: ',yop
       call stop_program('lcmpdtg called with illegal operand: ')
    endif

    if(kd1>kd2)then
       lcmpdtg=iop==1.or.iop==2.or.iop==6
    elseif(kd1==kd2)then
       if(kt1>kt2)then
          lcmpdtg=iop==1.or.iop==2.or.iop==6
       elseif(kt1==kt2)then
          lcmpdtg=iop==2.or.iop==3.or.iop==4
       else
          lcmpdtg=iop==4.or.iop==5.or.iop==6
       endif
    else
       lcmpdtg=iop==4.or.iop==5.or.iop==6
    endif
  end function lcmpdtg

  subroutine meanmnmx(y,field,klon,klat)
    use decomp
    implicit none
    character(len=*)::y
    integer,intent(in)::klon,klat
    real(kind=realkind),intent(in):: field(klon,klat)
    real(kind=realkind) rmin,rmax,rsum,rminG,rmaxG,rsumG
    integer i,j,is,ie,js,je
#ifdef MPI_SRC  
#include"mpif.h"
    integer::ierr
#endif  
    call jlimits(klon,klat,is,ie,js,je)
    rmin = minval(field(is:ie,js:je))
    rmax = maxval(field(is:ie,js:je))
    rsum = sum(field(is:ie,js:je))
    rminG = rmin
    rmaxG = rmax
    rsumG = rsum
#ifdef MPI_SRC  
    call MPI_reduce(rmin,rminG,1,REALTYPE,MPI_MIN,0,LOCALCOMM,ierr)
    call MPI_reduce(rmax,rmaxG,1,REALTYPE,MPI_MAX,0,LOCALCOMM,ierr)
    call MPI_reduce(rsum,rsumG,1,REALTYPE,MPI_SUM,0,LOCALCOMM,ierr)
#endif
    if(mype==0)then
       write(6,'(a,a,a,3e12.4)')'uhtmp ',y,' minmaxmean ',rminG,rmaxG,&
            rsumG/real(klon_global*klat_global,realkind)
    endif
    return
  end subroutine meanmnmx

end module util
