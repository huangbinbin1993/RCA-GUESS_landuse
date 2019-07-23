module escom
  use confys
  use ctun, only:acrit
  use decomp

  implicit none
  private
  real(kind=realkind),parameter:: xndegr=100.0_realkind
  integer,parameter::nstarti=1316,nstopi=37316
  real(kind=realkind),save:: estab(nstarti:nstopi)
  real(kind=realkind),save:: destab(nstarti:nstopi)
  real(kind=realkind),save:: estabw(nstarti:nstopi)
  real(kind=realkind),save:: destabw(nstarti:nstopi)
  real(kind=realkind),save:: estabi(nstarti:nstopi)
  real(kind=realkind),save:: destabi(nstarti:nstopi)
  logical,save::esatsfun=.false.

  real(kind=realkind),parameter::c1es=610.78_realkind
  real(kind=realkind),parameter::c2es=17.269_realkind
  real(kind=realkind),parameter::c2is=21.875_realkind
  real(kind=realkind),parameter::c3es=273.16_realkind !tmelt from confys
  real(kind=realkind),parameter::c4es=35.86_realkind
  real(kind=realkind),parameter::c4is=7.66_realkind
  real(kind=realkind),parameter::zc5p=10000.0_realkind
  real(kind=realkind),parameter::zc6e=0.378_realkind

  public tabdef,esatw,esati,esat,desdt,desdtw,desdti,pecor,crihum, esatfun
contains 
  real(kind=realkind) function zfracp(pp)
    implicit none
    real(kind=realkind),intent(in)::pp
    zfracp = min( pp, zc5p )/zc5p
  end function zfracp

  real(kind=realkind) function pecor(pp,tt)
    implicit none
    real(kind=realkind),intent(in)::pp,tt
    pecor = pp-zfracp(pp)*zc6e*esat(tt)
  end function pecor

  real(kind=realkind) function esatw(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       esatw=esatwfun(tt)
    else
       esatw=estabw(nint(xndegr*tt))
    endif
  end function esatw

  real(kind=realkind) function esati(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       esati=esatifun(tt)
    else
       esati=estabi(nint(xndegr*tt))
    endif
  end function esati

  real(kind=realkind) function esat(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       esat=esatfun(tt)
    else
       esat=estab(nint(xndegr*tt))
    endif
  end function esat

  real(kind=realkind) function desdt(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       desdt=desdtfun(tt)
    else
       desdt=destab(nint(xndegr*tt))
    endif
  end function desdt

  real(kind=realkind) function desdtw(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       desdtw=desdtwfun(tt)
    else
       desdtw=destabw(nint(xndegr*tt))
    endif
  end function desdtw

  real(kind=realkind) function desdti(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    if(esatsfun)then
       desdti=desdtifun(tt)
    else
       desdti=destabi(nint(xndegr*tt))
    endif
  end function desdti


  subroutine tabdef()

    implicit none
    integer jk
    real(kind=realkind) zdiv,ztemp
    zdiv = 1.0_realkind/xndegr
    do jk=nstarti,nstopi
       ztemp = real(jk,realkind)*zdiv
       estab(jk) = esatfun(ztemp)
       destab(jk) = desdtfun(ztemp)

       if( ztemp<36.0_realkind)then !Avoid division by zero in tt-c4is or tt-c4es :
          estabw(jk) = 0.0_realkind
          destabw(jk) = 0.0_realkind
          estabi(jk) = 0.0_realkind
          destabi(jk) = 0.0_realkind
       else
          estabw(jk) = esatwfun(ztemp)
          destabw(jk) = desdtwfun(ztemp)
          estabi(jk) = esatifun(ztemp)
          destabi(jk) = desdtifun(ztemp)
       endif
       if (ztemp>tmelt) then !Avoid ice calculation above freezing level:
          estabi(jk) = esatwfun(ztemp)
          destabi(jk) = desdtwfun(ztemp)
       endif
    enddo
    return
  end subroutine tabdef

  real(kind=realkind) function esatwfun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,c2es*(tt-c3es)/(tt-c4es))
    esatwfun=c1es*exp(tmp)
  end function esatwfun

  real(kind=realkind) function esatifun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,c2is*(tt-c3es)/(tt-c4is))
    esatifun=c1es*exp(tmp)
  end function esatifun

  real(kind=realkind) function zes(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    zes=zc2(tt)*(tt-c3es)/(tt-zc4(tt))
  end function zes

  real(kind=realkind) function esatfun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,zes(tt))
    esatfun=c1es*exp(tmp) !this gives underflow
  end function esatfun

  real(kind=realkind) function desdtfun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,zes(tt))
    desdtfun=c1es*exp(tmp)*zc2(tt)*(c3es-zc4(tt))/(tt-zc4(tt))**2.0_realkind
  end function desdtfun

  real(kind=realkind) function desdtwfun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,c2es*(tt-c3es)/(tt-c4es))
    desdtwfun = c1es*exp(tmp)*c2es*(c3es-c4es)/(tt-c4es)**2.0_realkind! pure saturation pressure over water
  end function desdtwfun

  real(kind=realkind) function desdtifun(tt)
    implicit none
    real(kind=realkind),intent(in)::tt
    real(kind=realkind)::tmp

    tmp = max(-89._realkind,c2is*(tt-c3es)/(tt-c4is))
    desdtifun =  c1es*exp(tmp)*c2is*(c3es-c4is)/(tt-c4is)**2.0_realkind! pure saturation pressure over ice (Has no meaning for T > 0 C)
  end function desdtifun



  real(kind=realkind) function zc2(tt)

    implicit none
    real(kind=realkind) fice,tt
    call ficefun(tt,fice)
    zc2=c2is*fice+(1.0_realkind-fice)*c2es
    return
  end function zc2

  real(kind=realkind) function zc4(tt)

    implicit none
    real(kind=realkind) tt,fice
    call ficefun(tt,fice)
    zc4=c4is*fice+(1.0_realkind-fice)*c4es
    return
  end function zc4


  subroutine crihum(klon,klat,klev,pahyb,pbhyb,pps,pt,pq) !,prpk )

    !crihum - critical humidity
    !     
    !     j.e. haugen           hirlam
    !     k.s. eerola           hirlam4 (revised)
    !     
    !     purpose.
    !     --------
    !     
    !     check humidity for negative values and
    !     remove moisture if values greater that acrit*saturation
    !     
    !     full level values of pressure is computed as the mean between
    !     half level values, instead of the more consistent

    !     1./ p(k) = ln ( p(k+0.5) / p(k-0.5) ) / ( p(k+0.5) - p(k-0.5) )
    !     crihum is called from gemini
    !     
    !     input parameters:
    !     -----------------
    !     
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of levels
    !     pahyb     vertical a-parameter at half levels
    !     pbhyb     vertical b-parameter at half levels
    !     pcrit     critical humidity parameter
    !     pps       surface pressure
    !     pt        temperature
    !     pq        specific humidity

    implicit none

    integer::  klon, klat, klev
    real(kind=realkind) pahyb(klev+1), pbhyb(klev+1),pps(klon,klat) , &
         pt(klon,klat,klev), pq(klon,klat,klev)
    real(kind=realkind)::prpk(klon,klat)


    integer:: jk,jx,jy
    real(kind=realkind)::  zqsat, zbk , zak

    do  jk=1,klev
       zak = 0.5_realkind*( pahyb(jk) + pahyb(jk+1) )
       zbk = 0.5_realkind*( pbhyb(jk) + pbhyb(jk+1) )
       do jy = 1,klat
          do jx = 1,klon
             prpk(jx,jy) = 1.0_realkind / (zak + zbk*pps(jx,jy))
          enddo
       enddo
       do jy = 1,klat
          do jx = 1,klon
             zqsat=epsilo*esat(pt(jx,jy,jk))*prpk(jx,jy)*acrit
             if (pq(jx,jy,jk)>zqsat) pq(jx,jy,jk) = zqsat
             if (pq(jx,jy,jk)<0.0_realkind   ) pq(jx,jy,jk) = 0.6e-6_realkind
          enddo
       enddo
    enddo

    return

    !  call swap(pq,klon,klat,klev)

  end subroutine crihum



end module escom
