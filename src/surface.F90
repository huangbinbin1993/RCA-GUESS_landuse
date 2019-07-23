module surface
  use decomp,only:realkind
  implicit none
  private
  type,public::surfVariables
     real(kind=realkind),allocatable,dimension(:,:):: tsm  ! surface temperature
     real(kind=realkind),allocatable,dimension(:,:):: tdm  ! deep surface temperature
     real(kind=realkind),allocatable,dimension(:,:):: tsc  ! climatological deep surface temprature
     real(kind=realkind),allocatable,dimension(:,:):: tsea ! sea surface temperature
     real(kind=realkind),allocatable,dimension(:,:):: swm  ! soil water content
     real(kind=realkind),allocatable,dimension(:,:):: sdm  ! deep soil water content
     real(kind=realkind),allocatable,dimension(:,:):: swc  ! climatological deep soil water content
     real(kind=realkind),allocatable,dimension(:,:):: snm  ! snow depth
     real(kind=realkind),allocatable,dimension(:,:):: snc  ! climatological snow depth
     real(kind=realkind),allocatable,dimension(:,:):: rou  ! roughness length
     real(kind=realkind),allocatable,dimension(:,:):: roc  ! climatological roughness length
     real(kind=realkind),allocatable,dimension(:,:):: alb  ! albedo
     real(kind=realkind),allocatable,dimension(:,:):: fri  ! fraction of ice
     real(kind=realkind),allocatable,dimension(:,:):: frf  
  end type surfVariables


  real(kind=realkind),public,save::vsw(12)
  real(kind=realkind),public,save::vcc(12)
  real(kind=realkind),public,save::vfl(12)
  real(kind=realkind),public,save::bw(12)
  real(kind=realkind),public,save::psis(12)
  real(kind=realkind),public,save::ks(12)
  real(kind=realkind),public,save::cquartz(12)
  real(kind=realkind),public,save::sfdist
  logical,public,save::snowint=.false.

  public inisurf,int_snow,conv_ecocli,slfluxo_land,calctemps
  public surf_land,slfluxo_surf_sea_ice,slfluxo_average,surf_land_tiles_splitup

contains

#include "surface_parameterisations/calctemps.incf90"
#include "surface_parameterisations/conv_ecocli.incf90"
#include "surface_parameterisations/slfluxo_average.incf90"
#include "surface_parameterisations/slfluxo_land.incf90"
#include "surface_parameterisations/slfluxo_surf_sea_ice.incf90"
#include "surface_parameterisations/surf_land.incf90"

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine inisurf
    use comrpar
    implicit none
    !
    !
    ! texture
    ! nr      texture
    ! ------- -------
    !   1     silty loam
    !   2     sand
    !   3     silty clay loam
    !   4     loam
    !   5     clay loam
    !   6     sany loam
    !   7     silty clay
    !   8     sandy clay loam
    !   9     loamy sand
    !  10     clay
    !  11     silt
    !  12     sandy clay
    !
    ! np89 m3/m3
    ! vfl vsw for peat dwd=7:mc cumber-pielke tab1;peat vcc:interp. 
    ! psis bw:clapp hornberger  peat: mc cumber pielke
    !
    ! hydraulic conductivity m s-1 * 10**6: same sources
    !
    ! diffusion of temperature
    ! dry volumetr.heat capac.(1.e-3*j m-3 k-1) mc cumber pielke tab1
    !
    !
    ! total porosity (%)
    vsw(1)=0.485_realkind
    vsw(2)=0.395_realkind
    vsw(3)=0.477_realkind
    vsw(4)=0.451_realkind
    vsw(5)=0.476_realkind
    vsw(6)=0.435_realkind
    vsw(7)=0.492_realkind
    vsw(8)=0.420_realkind
    vsw(9)=0.410_realkind
    vsw(10)=0.482_realkind
    vsw(11)=0.485_realkind
    vsw(12)=0.426_realkind
    !
    ! field capacity (m3/m3)
    vcc(1)=0.369_realkind
    vcc(2)=0.174_realkind
!cgj    vcc(2)=0.1_realkind	!value more in midrange of quoted values in literature
    vcc(3)=0.357_realkind
    vcc(4)=0.314_realkind
    vcc(5)=0.391_realkind
    vcc(6)=0.249_realkind	
!cgj    vcc(6)=0.2_realkind	!value more in midrange of quoted values in literature
    vcc(7)=0.409_realkind
    vcc(8)=0.299_realkind
    vcc(9)=0.179_realkind
    vcc(10)=0.400_realkind
    vcc(11)=0.369_realkind
    vcc(12)=0.316_realkind
    !
    ! wilting point (m3/m3)
    vfl(1)=0.179_realkind
    vfl(2)=0.068_realkind
!cgj    vfl(2)=0.045_realkind	!value more in midrange of quoted values in literature
    vfl(3)=0.218_realkind
    vfl(4)=0.155_realkind
    vfl(5)=0.250_realkind
    vfl(6)=0.114_realkind
    vfl(7)=0.283_realkind
    vfl(8)=0.175_realkind
    vfl(9)=0.075_realkind
    vfl(10)=0.286_realkind
    vfl(11)=0.179_realkind
    vfl(12)=0.219_realkind
    !
    ! sat. soil matric pot. (m)
    psis(1)=0.786_realkind
    psis(2)=0.121_realkind
    psis(3)=0.356_realkind
    psis(4)=0.478_realkind
    psis(5)=0.630_realkind
    psis(6)=0.218_realkind
    psis(7)=0.490_realkind
    psis(8)=0.299_realkind
    psis(9)=0.090_realkind
    psis(10)=0.405_realkind
    psis(11)=0.786_realkind
    psis(12)=0.153_realkind
    !
    ! clapp and hornberger exponent, (b-parameter)
    bw(1)=5.30_realkind
    bw(2)=4.05_realkind
    bw(3)=7.75_realkind
    bw(4)=5.39_realkind
    bw(5)=8.52_realkind
    bw(6)=4.90_realkind
    bw(7)=10.40_realkind 
    bw(8)=7.12_realkind
    bw(9)=4.38_realkind
    bw(10)=11.40_realkind 
    bw(11)=5.30_realkind
    bw(12)=10.40_realkind 
    !
    ! dry soil density (kg/m3)
    !
    ! is calculated in surf_land according to
    ! (1-vsw)*2650
    ! where 2650 is density of solid material
    !
    ! sat. hyd. cond. (m/s)
    ks(1)=7.20_realkind
    ks(2)=176.0_realkind
    ks(3)=1.70_realkind
    ks(4)=6.95_realkind
    ks(5)=2.45_realkind
    ks(6)=34.7_realkind
    ks(7)=1.03_realkind
    ks(8)=6.30_realkind
    ks(9)=156.0_realkind
    ks(10)=1.28_realkind
    ks(11)=7.20_realkind
    ks(12)=2.17_realkind
    !
    ! fraction of quartz, see peters-lidard et al. jas 1998, vol 55, nr 7, pp1209-1224
    cquartz(1)=0.25_realkind
    cquartz(2)=0.92_realkind
    cquartz(3)=0.10_realkind
    cquartz(4)=0.40_realkind
    cquartz(5)=0.35_realkind
    cquartz(6)=0.60_realkind
!    cquartz(6)=0.30_realkind	!slu maps say 30% for Sweden
    cquartz(7)=0.10_realkind
    cquartz(8)=0.60_realkind
    cquartz(9)=0.82_realkind
    cquartz(10)=0.25_realkind
    cquartz(11)=0.10_realkind
    cquartz(12)=0.52_realkind
    !
    sfdist=0.6_realkind
    !
    ! Thickness of soil layers w.r.t. temperature
    !
    dz1=1.e-2_realkind
    dz2=6.2e-2_realkind
    dz3=0.21_realkind
    dz4=0.72_realkind
    dz5=1.89_realkind
    !
    ! Thickness of soil layers w.r.t. soil moisture
    !
    dz1w=0.072_realkind
    dz2w=0.21_realkind
    !
    ! albedo specifications:
    albsnow=0.75_realkind
    !cgj300712albsnow=0.70_realkind
    !
    alboplveg=0.19_realkind
    alboplsoil=0.19_realkind
    albforconif=0.10_realkind
    albfordecid=0.15_realkind
    albforfloor=0.15_realkind
    albdesert=0.33_realkind
    !
    albsnfor=0.50_realkind
    albsnlmin=0.50_realkind     !rca4cordex was 0.6
    albsnlmax=0.85_realkind
    !cgj300712albsnlmin=0.40_realkind
    !cgj300712albsnlmax=0.70_realkind
    albice=0.50_realkind   ! good for sea ice
    albwater=0.07_realkind
    !
    ! emissivity specifications:
    !
    emsnow=0.99_realkind
    !
    emoplveg=0.96_realkind
    emoplsoil=0.96_realkind
    emforconif=0.97_realkind
    emfordecid=0.97_realkind
    emforfloor=0.97_realkind
!    emdesert=0.94_realkind
    emdesert=0.94_realkind
    !
    emsnlmin=0.97_realkind
    emsnlmax=0.99_realkind
    emice=0.94_realkind
    emwater=0.94_realkind
    !
    return
  end subroutine inisurf

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine int_z0for ( month,day,z0for_int)
    implicit   none
    integer    month,day
    real(kind=realkind)       z0for_int,z0for(12)

    integer    m1,m2
    real(kind=realkind)       s1,s2

    !gj280705 suggested values of annual cycle of z0 for mixed forest from
    !gj280705 http://ldas.gsfc.nasa.gov/ldas8th/mapped.veg/web.veg.monthly.table.html

    data z0for /0.81_realkind,0.81_realkind,0.88_realkind,1.0_realkind,&
         1.05_realkind,1.06_realkind,1.06_realkind,1.06_realkind,&
         1.06_realkind,1.0_realkind,0.88_realkind,0.81_realkind/

    print *,'This is not true for any month...'
    print *,'use functionality from calendar'
    print *,__FILE__,__LINE__
    if ( day > 15 ) then
       m1 = month
       m2 = m1+1
       if (m2>12) m2=1
       s2 = abs((real(day,realkind)-15._realkind))/30._realkind
       s1 = 1._realkind - s2
    else
       m2 = month
       m1 = m2-1
       if (m1==0) m1=12
       s1 = abs((real(day,realkind)-15._realkind))/30._realkind
       s2 = 1._realkind - s1
    endif
    z0for_int = s1 * z0for(m1) + s2 * z0for(m2)
    return
  end subroutine int_z0for

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine int_snow ( month,day,rhosnow_int)
    implicit   none
    integer,intent(in):: month,day
    real(kind=realkind)       rhosnow_int,rhosnow(12)

    integer    m1,m2
    real(kind=realkind)       s1,s2

    !mtr980323 "typical (?)" monthly density of snow from january-december.
    data rhosnow /220._realkind,230._realkind,240._realkind,280._realkind,&
         320._realkind,320._realkind,320._realkind,320._realkind,&
         100._realkind,160._realkind,180._realkind,210._realkind/
    !
!    print *,'This is not true for any month...'
!    print *,'use functionality from calendar'
!    print *,__FILE__,__LINE__
    if ( day > 15 ) then
       m1 = month
       m2 = m1+1
       if (m2>12) m2=1
       s2 = abs((real(day,realkind)-15._realkind))/30._realkind
       s1 = 1._realkind - s2
    else	
       m2 = month
       m1 = m2-1
       if (m1==0) m1=12
       s1 = abs((real(day,realkind)-15._realkind))/30._realkind
       s2 = 1._realkind - s1
    endif
    rhosnow_int = s1 * rhosnow(m1) + s2 * rhosnow(m2)
    return
  end subroutine int_snow


  subroutine tridag(a,b,c,r,u,n)
    implicit none
    integer,parameter::nmax=100
    real(kind=realkind):: gam(nmax),a(n),b(n),c(n),r(n),u(n)
    real(kind=realkind):: bet
    integer:: j,n
    bet=b(1)
    u(1)=r(1)/bet
    do  j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    enddo
    return
  end subroutine tridag



  subroutine surf_land_tiles_splitup(rzcw,rzsnopl,rzsnoplm,rzsnfor, rzsnform, &
       rlandon,   rsrfdist,  rfrlim,   loldlim, &
       rfrop,     rfrsnw,    rfrsnti,  rfrsnfor, rzsnoplmu, rzsnformu  )
    ! calculation of the tile areas for a given snowload
    ! input
    !             rzcw         forest fraction
    !             rzsnopl      total open land snow load in gridbox
    !             rzsnoplm     seasonal max open land snow load
    !             rzsnfor      total forest snow load in gridbox
    !             rzsnform     seasonal max forest snow load
    !             rlandon      if 1 when land
    !             rsrfdist     orographic distribution factor
    !             rfrlim       mimimum tile area
    !             loldlim      use old rca35 calculations
    ! output
    !             rfrop        fraction of grid box that is snow free
    !             rfrsnw       fraction of grid box that is snow covered  
    !             rfrsnti      fraction of open land tile that is snow covered
    !             rfrsnfor     fraction of forest tile that is snow covered  
    !             rzsnoplmu    updated seasonal max open land snow load
    !             rzsnformu    updated seasonal max forest snow load
    !
    ! version 1.1 wj van de berg, 6 nov 2009
    !
    real(kind=realkind)    rzcw,      rzsnopl,   rzsnoplm, rzsnfor, rzsnform
    real(kind=realkind)    rlandon,   rsrfdist,  rfrlim
    logical loldlim
    real(kind=realkind)    rfrop,     rfrsnw,    rfrsnti,  rfrsnfor
    real(kind=realkind)    rzsnoplmu, rzsnformu  

    real(kind=realkind)    zopa
    real(kind=realkind)    zsnoplr,   zsnforr,   zsnoplg,  zfreeoplg

    if ( rlandon <= 0._realkind ) then
       rfrop     = 0._realkind
       rfrsnw    = 0._realkind
       rfrsnti   = 0._realkind
       rfrsnfor  = 0._realkind 
       rzsnoplmu = rzsnoplm
       rzsnformu = rzsnform
       return
    endif

    zopa = 1._realkind - rzcw

    call surf_snowfraction(zopa,rzsnopl,rzsnoplm,rsrfdist, loldlim,&
         zsnoplr,  rzsnoplmu )

    call surf_snowfraction(rzcw,rzsnfor,rzsnform,rsrfdist,loldlim, &
         zsnforr,  rzsnformu )

    ! apply rules on min/max allowable areas      
    ! and calculate final areas
    if ( loldlim ) then
       ! require that all fractions are larger than rfrlim
       ! and that the snow free open land area is always larger than rfrlim.
       zfreeoplg = zopa * (1._realkind-zsnoplr)
       zsnoplg   = zopa * zsnoplr

       if ( zfreeoplg<rfrlim ) zsnoplg = zopa - rfrlim
       if ( zsnoplg  <rfrlim ) zsnoplg = 0._realkind
       if ( zsnforr*rzcw<rfrlim ) zsnforr = 0._realkind

       rfrop    = ( zopa - zsnoplg ) * rlandon
       rfrsnw   = (        zsnoplg ) * rlandon
       rfrsnti  = ( zsnoplg / zopa ) * rlandon
       rfrsnfor =   zsnforr          * rlandon
    else
       ! only require that the snow cover is more than rfrlim of the tile fraction
       if ( zsnoplr<rfrlim ) zsnoplr = 0._realkind
       if ( zsnforr<rfrlim ) zsnforr = 0._realkind          

       rfrop    = zopa * (1._realkind-zsnoplr) * rlandon
       rfrsnw   = zopa * zsnoplr      * rlandon
       rfrsnti  =        zsnoplr      * rlandon
       rfrsnfor = zsnforr             * rlandon
    endif

    return
  end subroutine surf_land_tiles_splitup


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine surf_snowfraction(rfrac,rsnow,rsnowm, rsrfdist, loldlim,&
       rsnwfrac, rsnowmp )



    ! calculation of the snowfraction for a given snowload
    ! formula works for open land and forest
    ! input
    !             rfrac        total area of tile (snowfree+snowcovered)
    !             rsnow        total snowload in tile on whole gridbox
    !             rsnowm       seasonal max total snowload in tile on whole gridbox 
    !             rsrfdist     orographic distribution factor
    !             loldlim      use old rca35 calculations
    ! output
    !             rsnwfrac     snow covered fraction (1=whole tile)
    !             rsnowmp      updated seasonal max total snowload in tile on whole gridbox
    ! 
    ! version 1.0 wj van de berg, 5 nov 2009
    !
    implicit none
    ! input
    real(kind=realkind)    rfrac,      rsnow,   rsnowm
    real(kind=realkind)    rsrfdist
    logical loldlim
    !output
    real(kind=realkind)    rsnwfrac,   rsnowmp

    ! constants
    real(kind=realkind)    zfrsnasymp, zsnlim

    ! local 
    real(kind=realkind)    zfinv,      zsn0,    zslask      


    ! return to main function if whole tile is nonexistant      
    if ( rfrac < 1.e-6_realkind ) then
       rsnwfrac = 0._realkind
       rsnowmp  = 0._realkind
       return
    endif

    zfrsnasymp = 0.985_realkind
    zsnlim     = 0.0015_realkind

    zfinv      = 1._realkind/rfrac
    rsnowmp    = rsnowm

    if ( loldlim ) then ! do it the old way	      
       if ( rsnowm < 1.e-6_realkind) then
          rsnwfrac = zfrsnasymp*tanh(100._realkind*rsnow)
          rsnowmp  = 0._realkind
          if (rsnwfrac>(zfrsnasymp-0.001_realkind)) rsnowmp = rsnow
       else
          zsn0 = max(0._realkind,rsnow) 
          if ( zsn0>=zsnlim ) then
             rsnwfrac = min(   zsn0/(rsnowm*rsrfdist), zfrsnasymp )
          else
             zslask   = min( zsnlim/(rsnowm*rsrfdist), zfrsnasymp )
             rsnwfrac = zsn0*zslask/zsnlim
          endif
          rsnwfrac = min(max( 0.0_realkind, rsnwfrac), zfrsnasymp ) 
       endif
    else      
       if ( rsnowm < rfrac*1.e-6_realkind) then
          rsnwfrac = zfrsnasymp*tanh(100._realkind*rsnow*zfinv)
          rsnowmp  = 0._realkind
          if ( rsnwfrac>(zfrsnasymp-0.001_realkind) ) rsnowmp = rsnow
       else
          zsn0 = max(0._realkind,rsnow*zfinv) ! snow amount if tile would cover whole gridbox
          if ( zsn0>=zsnlim ) then
             rsnwfrac = min(   zsn0/(rsnowm*rsrfdist), zfrsnasymp )
          else
             zslask   = min( zsnlim/(rsnowm*rsrfdist), zfrsnasymp )
             rsnwfrac = zsn0*zslask/zsnlim
          endif
          rsnwfrac = min(max( 0.0_realkind, rsnwfrac), zfrsnasymp ) 
       endif
    endif
    return
  end subroutine surf_snowfraction

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module surface
