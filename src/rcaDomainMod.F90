module rcaDomainMod
  use filter_orog
  use domainMod
  use decomp
  use util
  implicit none
  private
  real(kind=realkind),public,save::rdlam  ! inverse grid distance in x-direction (radians)
  real(kind=realkind),public,save::rdth   ! inverse grid distance in y-direction (radians)
  real(kind=realkind),public,save::ra     ! inverse radius of the earth (meter)
  real(kind=realkind),public,save::ffmean ! mean value of coriolis parameter
  real(kind=realkind),public,save::rhxumm ! inverse mean of mapfactor in x-direction
  real(kind=realkind),public,save::rhyvmm ! inverse mean of mapfactor in y-direction


  real(kind=realkind),public,save,allocatable,dimension(:,:)::fpar ! coriolis parameter
  real(kind=realkind),public,save,allocatable,dimension(:,:)::hxv ! map-factor in x-direction, v-points
  real(kind=realkind),public,save,allocatable,dimension(:,:)::hyu ! map-factor in y-direction, u-points
  real(kind=realkind),public,save,allocatable,dimension(:,:)::rhxu ! inverse map-factor in x-dirctn, u-pts
  real(kind=realkind),public,save,allocatable,dimension(:,:)::rhyv! inverse map-factor in y-dirctn, v-pts 
  real(kind=realkind),public,save,allocatable,dimension(:,:)::plong! longitude
  real(kind=realkind),public,save,allocatable,dimension(:,:)::pslat! sine(latitude)
  real(kind=realkind),public,save,allocatable,dimension(:,:)::pclat! cosine(latitude)
  real(kind=realkind),public,save,allocatable,dimension(:,:)::lat! (latitude)
  real(kind=realkind),public,save,allocatable,dimension(:,:)::lon! (latitude)

  
  type(domain),public,save,target::RCAdomain

  public fld_mean, fld_mean_dbl, initRca
contains

  subroutine initRca(klon,klat,klev)
    use confys, only:gravit
    use referenceParameters,only:sip0
    implicit none
    external readgtopo30data
    integer,intent(in)::klon,klat,klev
    integer::halo,jk,i,j
    logical::lfilter
    real(kind=realkind)::rlat(klon,klat),rlon(klon,klat)
    halo=1
    lfilter=.true.
    allocate(RCAdomain%ahyb(klev+1),RCAdomain%bhyb(klev+1),RCAdomain%hybi(klev+1),&
         RCAdomain%hybk(klev),RCAdomain%afull(klev),RCAdomain%bfull(klev))
    allocate(RCAdomain%fis(klon,klat))
    allocate(RCAdomain%orogsigm(klon,klat))

    call readgtopo30data(klon,klat,klon_global,klat_global,halo, &
         RCAdomain%south,RCAdomain%west,RCAdomain%dlon,RCAdomain%dlat,&
         RCAdomain%polon,RCAdomain%polat,RCAdomain%fis,RCAdomain%orogsigm, &
         idatastart,jdatastart,gravit,localComm)

    allocate(fpar(klon,klat))
    allocate(hxv(klon,klat))
    allocate(hyu(klon,klat))
    allocate(rhxu(klon,klat))
    allocate(rhyv(klon,klat))
    allocate(plong(klon,klat))
    allocate(pslat(klon,klat))
    allocate(pclat(klon,klat))
    allocate(lat(klon,klat),lon(klon,klat))
    do j=1,klat
       do i=1,klon
          rlat(i,j) = RCAdomain%south + (j+jdatastart-2)*RCAdomain%dlat
          rlon(i,j) = RCAdomain%west  + (i+idatastart-2)*RCAdomain%dlon
       enddo
    enddo
    call regrot(lon,lat,rlon,rlat,klon,klat,RCAdomain%polon,RCAdomain%polat,-1)


    lfilter=.true.
    if (lfilter)then
       call filter_orography(RCAdomain%fis,klon,klat)
       if(mype==0) write(6,*) ' filtering of orography done'
    else
       if(mype==0) write(6,*) ' filtering of orography not done'
    endif
    call mapfac(klon,klat)

    !initiate all pressure/hybrid level coefficients
    call get_a_b(RCAdomain%ahyb,RCAdomain%bhyb,klev)

    do  jk=1,klev+1
       RCAdomain%hybi(jk) = RCAdomain%ahyb(jk)/sip0 + RCAdomain%bhyb(jk)
    enddo
    do jk=1,klev
       RCAdomain%hybk(jk) = (RCAdomain%hybi(jk) + RCAdomain%hybi(jk+1))/2.0_realkind
    enddo
    do jk=1,klev
       RCAdomain%afull(jk) = 0.5_realkind*(RCAdomain%ahyb(jk)+RCAdomain%ahyb(jk+1) )
       RCAdomain%bfull(jk) = 0.5_realkind*(RCAdomain%bhyb(jk)+RCAdomain%bhyb(jk+1) )
    enddo
    RCAdomain%nlev = klev_global
    RCAdomain%nlon = klon_global
    RCAdomain%nlat = klat_global
    call arakawastag('c',RCAdomain%stagu_x,RCAdomain%stagu_y,RCAdomain%stagv_x,RCAdomain%stagv_y)

  end subroutine initRca

  subroutine mapfac(klon, klat)
    !
    ! mapfac - geographically dependent variables
    !
    !     j.e. haugen         hirlam
    !     k.s. eerola         hirlam4(revised)
    !     compute coriolis-parameter and mapping factors
    !             sin(latitude) - pslat,
    !             cos(latitude) - pclat,
    !             longtude      - plong.
    !
    !     input parameters:
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels

    use decomp  
    !    use comhkp  
    use confys
    implicit none  
    integer :: klon, klat  
    integer :: jx, jy
    real(kind=realkind) :: zpi, zpir18,  zsc, zwc,  zdlam, zdth, zlonc, &
         zlatc, zrhxu, zhxv, z18rpi, zomega, zloncx, zlatcx

    real(kind=realkind) :: zlon(klon, klat), zlat(klon, klat), zxlon(klon),zylat(klat)

    if(mype==0)print *, ' -------- in mapfac ------------'

    zpi = 2.0_realkind * asin(1._realkind)  
    zpir18 = zpi / 180.0_realkind  

    zsc = RCAdomain%south * zpir18  
    zwc = RCAdomain%west * zpir18  


    zdlam = RCAdomain%dlon * zpir18  
    zdth = RCAdomain%dlat * zpir18  

    rdlam = 1.0_realkind / zdlam  
    rdth = 1.0_realkind / zdth  
    ra = 1.0_realkind / rearth

    zlonc = RCAdomain%polon * zpir18  
    zlatc = zpi * 0.5_realkind + RCAdomain%polat * zpir18  

    !    metric coefficients used in dynamic routines
    do jy = 1, klat  
       do jx = 1, klon  
          rhyv(jx, jy) = 1.0_realkind  
          hyu(jx, jy) = 1.0_realkind  
       enddo
    enddo
    !
    do jy = 1, klat  
       zrhxu = 1._realkind/cos(zsc+real(jy+jdatastart-2,realkind)*zdth)  
       zhxv = cos(zsc +(real(jy + jdatastart,realkind) - 1.5_realkind) * zdth)  
       do jx = 1, klon  
          rhxu(jx, jy) = zrhxu  
          hxv(jx, jy) = zhxv  
       enddo
    enddo
    !     compute convensional longitude and latitude in lonlat
    do jx = 1, klon  
       zxlon(jx) = zwc + real(jx + idatastart - 2,realkind) * zdlam  
    enddo
    do jy = 1, klat  
       zylat(jy) = zsc + real(jy + jdatastart - 2,realkind) * zdth  
    enddo
    call lonlat(klon,klat,zlon,zlat,zxlon,zylat,zlonc,zlatc)
    !    compute sin(latitude), cos(latitude) and longitude in degrees
    z18rpi = 1.0_realkind / zpir18  
    do jy = 1, klat  
       do jx = 1, klon  
          plong(jx, jy) = zlon(jx,jy)*z18rpi  
          pslat(jx, jy) = sin(zlat(jx,jy))  
          pclat(jx, jy) = cos(zlat(jx,jy))  
       enddo
    enddo
    !    compute conventional latitude in lonlat and coriolis param.
    do jx = 1, klon  
       zxlon(jx) = zwc +(real(jx + idatastart - 1,realkind) - 0.5_realkind) * zdlam  
    enddo
    do jy = 1, klat  
       zylat(jy) = zsc +(real(jy + jdatastart - 1,realkind) - 0.5_realkind) * zdth  
    enddo
    call lonlat(klon,klat,zlon,zlat,zxlon,zylat,zlonc,zlatc)
    zomega = 4.0_realkind * zpi /(23.934_realkind * 60.0_realkind * 60.0_realkind)  
    do jy = 1, klat  
       do jx = 1, klon  
          fpar(jx, jy) = zomega*sin(zlat(jx, jy))  
       enddo
    enddo

    ! compute mean value of coriolis parameter
    ffmean = fld_mean(fpar,klon,klat)
    if(mype==0)print *,' mean value of coriolis, ffmean =', ffmean

    zloncx = zlonc / zpir18  
    zlatcx = zlatc / zpir18  
    if(mype==0)then
       print *,' zlonc ', zloncx  
       print *,' zlatc ', zlatcx
       print *, ' --------  mapfac done ------------'
    endif
    return  
  end subroutine mapfac


  subroutine lonlat(klon,klat,plon,plat,pxlon,pylat,plonc,platc)
    ! latlon - computes latitude and longitude
    !     compute convensional longitude and latitude(plon,plat)
    !     from coordinates pxlon and pylat used as model grid
    !
    !     input parameters:
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     kgrid     type of grid
    !     pxlon     model grid x-coordinate
    !     pylat     model grid y-coordinate
    !     plonc     x-position of rotated equator
    !     platc     y-position of rotated equator
    !
    !     output parameters:
    !     plon      conventional longitude
    !     plat      convensional latitude

    use decomp
    implicit none  

    integer,intent(in)::  klon, klat
    real(kind=realkind),intent(in):: plonc, platc  
    real(kind=realkind) :: plon(klon, klat), plat(klon, klat)
    real(kind=realkind),intent(in)::pxlon(klon),pylat(klat)
    !     declaration of lobal parameters
    integer :: jx, jy  

    real(kind=realkind) :: zcyc, zsyc, zpi, zpi2, zsy, zcy, zcx, zsx, zsyl, zcyl, &
         zcdx1, zsdx1, zdx1
    zpi = 4.0_realkind * atan(1.0_realkind)  
    zpi2 = zpi * 2.0_realkind  

    !    compute rotated longitude and latitude grid
    do  jy = 1, klat  
       do  jx = 1, klon  
          plon(jx, jy) = pxlon(jx)  
          plat(jx, jy) = pylat(jy)  
       enddo
    enddo
    !    compute convensional longitude and latitude
    zcyc = cos(platc)  
    zsyc = sin(platc)  

    do  jy = 1, klat  
       do  jx = 1, klon  
          zsy = sin(plat(jx, jy) )  
          zcy = cos(plat(jx, jy) )  
          zcx = cos(plon(jx, jy) )  
          zsx = sin(plon(jx, jy) )  
          zsyl = zsy * zcyc + zsyc * zcy * zcx  
          zsyl = min(zsyl, 1.0_realkind)  
          zsyl = max(zsyl, - 1.0_realkind)  
          plat(jx, jy) = asin(zsyl)  
          zcyl = cos(plat(jx, jy) )  
          if(zcyl<= 1.e-7_realkind) then  
             plon(jx, jy) = 0.0_realkind  
          else  
             zcdx1 =(zcyc * zcy * zcx - zsyc * zsy) / zcyl  
             zsdx1 = zcy * zsx / zcyl  
             zcdx1 = max(min(zcdx1, 1.0_realkind), - 1.0_realkind)  
             zdx1 = acos(zcdx1)  
             if(zsdx1<0.0_realkind) zdx1 = zpi2 - zdx1  
             plon(jx, jy) = zdx1 + plonc  
             if(plon(jx, jy) >zpi2)then
                plon(jx, jy) = plon(jx, jy) - zpi2
             endif
          endif

       enddo
    enddo
    return  
  end subroutine lonlat



  real(kind=realkind) function fld_mean(fld, klon, klat )
    !     calculate a mean of a field
    !     interface:
    !     description:
    !     collects the field on pe 0, calculates there the global
    !     mean and broadcasts the result to every pe
    use decomp
    implicit none
    integer klon,klat
    real(kind=realkind) fld(klon,klat)
    real(kind=realkind)  zsmean
    integer jmin,jmax,imin,imax,i,j
    real(kind=realkind) wrk
#ifdef MPI_SRC
#include"mpif.h"
    integer::ierr
#endif
    call jlimits(klon,klat,imin,imax,jmin,jmax)
    zsmean = 0.0_realkind
    do j=jmin,jmax
       do i=imin,imax
          zsmean = zsmean + fld(i,j)
       enddo
    enddo
    zsmean = zsmean/real(klon_global*klat_global,realkind)
#ifdef MPI_SRC
    wrk = zsmean
    call mpi_allreduce(wrk,zsmean,1,REALTYPE, mpi_sum,localComm,ierr)
#endif      
    fld_mean = zsmean
    return
  end function fld_mean


  real(kind=8) function fld_mean_dbl(fld,klon,klat)
    use decomp
    implicit none
    integer klon,klat
    real(kind=8) fld(klon,klat)
    real(kind=8) wrk,zsmean
#ifdef MPI_SRC
#include"mpif.h"
    integer ierr
#endif

    integer  i,j,imin,imax,jmin,jmax
    call jlimits(klon,klat,imin,imax,jmin,jmax)
    zsmean = 0d0
    do j=jmin,jmax
       do i=imin,imax
          zsmean = zsmean + fld(i,j)
       enddo
    enddo

    zsmean = zsmean/real(klon_global*klat_global,kind(fld))

#ifdef MPI_SRC
    wrk = zsmean
    call mpi_allreduce(wrk, zsmean,1,MPI_DOUBLE_PRECISION, mpi_sum,  localComm, ierr)
#endif
    fld_mean_dbl = zsmean
    return
  end function fld_mean_dbl


  subroutine get_a_b(ahalf,bhalf,klev)
    implicit none
    integer,intent(in)::klev
    real(kind=realkind):: ahalf(klev+1),bhalf(klev+1)

    select case(klev)!available are 19,24,31,40,50,60,62,91
    case(19)!from ecmwf
       ahalf(1:20)=(/0.0_realkind,2000.0_realkind,4000.0_realkind,6046.110595_realkind,&
            8267.927560_realkind,10609.513232_realkind,12851.100169_realkind,14698.498086_realkind, &
            15861.125180_realkind,16116.236610_realkind, 15356.924115_realkind, 13621.460403_realkind,&
            11101.561987_realkind,  8127.144155_realkind,  5125.141747_realkind, &
            2549.969411_realkind,   783.195032_realkind,     0.0_realkind,     0.0_realkind,     0.0_realkind/)

       bhalf(1:20)=(/0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0003389933_realkind,&
            0.0033571866_realkind,&
            0.0130700434_realkind,0.0340771467_realkind,0.0706498323_realkind,0.1259166826_realkind,&
            0.2011954093_realkind,0.2955196487_realkind,0.4054091989_realkind,&
            0.5249322235_realkind,0.6461079479_realkind,0.7596983769_realkind,0.8564375573_realkind,&
            0.9287469142_realkind,0.9729851852_realkind,0.9922814815_realkind,1.0000000000_realkind/)
    case(24)
       ahalf(1:25)=(/0.0_realkind,      3806.742399768224_realkind,      7653.614367556312_realkind,&
            11253.36349801443_realkind,      14390.64919396709_realkind,      16915.32999628355_realkind, &
            18735.76024534827_realkind,      19812.07464310479_realkind,      20149.48434288351_realkind,&
            19791.56461262767_realkind,      18813.54797638499_realkind,      17315.61247280784_realkind,&
            15416.17779901790_realkind,      13245.18809670882_realkind,      10937.41500946619_realkind,&
            8625.729944404449_realkind,      6434.416341760640_realkind,      4472.449976936451_realkind,&
            2826.782241573884_realkind,      1555.648116570220_realkind,      681.8463580350899_realkind,&
            186.0290483950827_realkind,   0.0_realkind,   0.0_realkind,  0.0_realkind/)
       bhalf(1:25)=(/0.0_realkind,    0.0_realkind, 2.2695333717478114E-09_realkind,&
            3.1713357874686724E-03_realkind, 1.1910524859928885E-02_realkind, 2.7905492032898480E-02_realkind,&
            5.2201797942528856E-02_realkind, 8.5268807095453703E-02_realkind, 0.1270659877502267_realkind,&
            0.1771089674283446_realkind,     0.2345359330504250_realkind,     0.2981737778528530_realkind,&
            0.3666041544631837_realkind,     0.4382300145782184_realkind,     0.5113413834541805_realkind,&
            0.5841821876920823_realkind,     0.6550157897438893_realkind,     0.7221915275815412_realkind,&
            0.7842112033585109_realkind,     0.8397947168279032_realkind,     0.8879466076938281_realkind,&
            0.9250224344760386_realkind,     0.9553045730432707_realkind,     0.9788056889896553_realkind,&
            1.000000000000000_realkind/)
    case(31) !from ecmwf
       ahalf(1:32)=(/0.000000_realkind, 2000.000000_realkind, 4000.000000_realkind, 6000.000000_realkind,&
            8000.000000_realkind, 9976.135361_realkind,11820.539617_realkind,&
            13431.393926_realkind,14736.356909_realkind,15689.207458_realkind,16266.610500_realkind,&
            16465.005734_realkind,16297.619332_realkind,15791.598604_realkind,14985.269630_realkind,&
            13925.517858_realkind,12665.291662_realkind,11261.228878_realkind, 9771.406290_realkind,&
            8253.212096_realkind, 6761.341326_realkind, 5345.914240_realkind, 4050.717678_realkind,&
            2911.569385_realkind, 1954.805296_realkind, 1195.889791_realkind,  638.148911_realkind,&
            271.626545_realkind, 72.063577_realkind, 0.000000_realkind,0.000000_realkind,0.000000_realkind/)

       bhalf(1:32)=(/0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0003908582_realkind,0.0029197006_realkind,&
            0.0091941320_realkind,0.0203191555_realkind,0.0369748598_realkind,0.0594876397_realkind,&
            0.0878949492_realkind,0.1220035886_realkind,0.1614415235_realkind,0.2057032385_realkind,&
            0.2541886223_realkind,0.3062353873_realkind,0.3611450218_realkind,0.4182022749_realkind,&
            0.4766881754_realkind,0.5358865832_realkind,0.5950842740_realkind,0.6535645569_realkind,&
            0.7105944258_realkind,0.7654052430_realkind,0.8171669567_realkind,0.8649558510_realkind,&
            0.9077158297_realkind,0.9442132326_realkind,0.9729851852_realkind,0.9922814815_realkind,&
            1.0000000000_realkind/)
    case(40)
       ahalf(1:41)=(/ 0.00000000_realkind,      2006.05757759_realkind,      3996.76329133_realkind,&
            5923.68296148_realkind,      7748.44259336_realkind,      9441.11750699_realkind,&
            10978.86779645_realkind,     12344.78661268_realkind,     13526.95383523_realkind,&
            14517.65048286_realkind,     15312.73221231_realkind,     15911.12695982_realkind,&
            16314.44423948_realkind,     16526.65906480_realkind,     16553.90675567_realkind,&
            16404.28906623_realkind,     16087.74976371_realkind,     15615.95510868_realkind,&
            15002.23043794_realkind,     14261.41117411_realkind,     13409.81329939_realkind,&
            12465.02671268_realkind,     11445.83987802_realkind,     10372.00083822_realkind,&
            9263.98994152_realkind,      8142.77678679_realkind,      7029.44805180_realkind,&
            5944.80979003_realkind,      4909.01397857_realkind,      3940.95327361_realkind,&
            3057.76346642_realkind,      2274.18680381_realkind,      1601.89182737_realkind,&
            1048.74500687_realkind,       618.01601383_realkind,       307.62477035_realkind,&
            109.20228347_realkind,         7.26729959_realkind,         0.00000000_realkind,&
            0.00000000_realkind,         0.00000000_realkind/)
       bhalf(1:41)=(/ 0.00000000_realkind,         0.00000000_realkind,         0.00000000_realkind, &
            0.00079242_realkind,         0.00304974_realkind,         0.00733081_realkind, &
            0.01408708_realkind,         0.02366891_realkind,         0.03633186_realkind, &
            0.05224293_realkind,         0.07148688_realkind,         0.09407249_realkind, &
            0.11993881_realkind,         0.14896147_realkind,         0.18095902_realkind,&
            0.21569908_realkind,         0.25290469_realkind,         0.29226059_realkind,&
            0.33341944_realkind,         0.37600830_realkind,         0.41963462_realkind,&
            0.46389262_realkind,         0.50836986_realkind,         0.55265285_realkind,&
            0.59633421_realkind,         0.63901802_realkind,         0.68032690_realkind,&
            0.71990793_realkind,         0.75743866_realkind,         0.79263406_realkind,&
            0.82525205_realkind,         0.85510040_realkind,         0.88204274_realkind,&
            0.90600465_realkind,         0.92698038_realkind,         0.94503868_realkind,&
            0.96032944_realkind,         0.97308980_realkind,         0.98365033_realkind,&
            0.99244155_realkind,         1.00000000_realkind/)
    case(50)!from ecmwf
       ahalf(1:51)=(/0.000000_realkind,20.006149_realkind,43.297810_realkind,75.346230_realkind,&
            115.082146_realkind,161.897491_realkind,215.896912_realkind,278.005798_realkind,350.138184_realkind,&
            435.562286_realkind,539.651489_realkind,668.615540_realkind,828.398987_realkind,&
            1026.366943_realkind,1271.644531_realkind,1575.537842_realkind,1952.054443_realkind,2418.549805_realkind,&
            2996.526611_realkind,3712.626221_realkind,4599.856934_realkind,5699.114746_realkind,&
            6998.388184_realkind,8507.411133_realkind,10181.707031_realkind, 11883.089844_realkind, &
            13442.915039_realkind,&
            14736.354492_realkind, 15689.206055_realkind, 16266.609375_realkind, 16465.003906_realkind,&
            16297.620117_realkind, 15791.597656_realkind, 14985.269531_realkind, 13925.519531_realkind,&
            12665.294922_realkind, 11261.230469_realkind,  9771.406250_realkind,  8253.210938_realkind,&
            6761.339844_realkind,  5345.917969_realkind,  4050.718750_realkind,  2911.570313_realkind,&
            1954.804688_realkind,  1195.890625_realkind,   638.148438_realkind,   271.625000_realkind,&
            72.062500_realkind,     0.000000_realkind,     0.000000_realkind,     0.000000_realkind/)

       bhalf(1:51)=(/0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,0.0000000000_realkind,&
            0.0001003604_realkind,0.0006727143_realkind,0.0031633405_realkind,0.0092923399_realkind,&
            0.0203191563_realkind,0.0369748585_realkind,0.0594876409_realkind,0.0878949761_realkind,&
            0.1220036149_realkind,0.1614415050_realkind,0.2057032585_realkind,0.2541885972_realkind,&
            0.3062353730_realkind,0.3611450195_realkind,0.4182022810_realkind,0.4766881466_realkind,&
            0.5358865857_realkind,0.5950842500_realkind,0.6535645723_realkind,0.7105944157_realkind,&
            0.7654052377_realkind,0.8171669841_realkind,0.8649558425_realkind,0.9077158570_realkind,&
            0.9442132115_realkind,0.9729852080_realkind,0.9922814965_realkind,1.0000000000_realkind/)
    case(60)
       ahalf(1:61)=(/0.00000000_realkind,      2000.05135426_realkind,      4033.40907891_realkind,&
            6074.40328124_realkind,      8099.63670588_realkind,     10087.89367010_realkind,    &
            12020.06658606_realkind,     13879.06676239_realkind,     15649.75041085_realkind,    &
            17318.83964105_realkind,     18874.85209658_realkind,     20308.02335589_realkind,&
            21610.23752964_realkind,     22774.95349175_realkind,     23797.14016962_realkind,&
            24673.19673574_realkind,     25400.90458140_realkind,     25979.33203475_realkind,&
            26408.79136269_realkind,     26690.74701491_realkind,     26827.77478534_realkind,&
            26823.47214753_realkind,     26682.39550616_realkind,     26409.99528681_realkind,&
            26012.54977341_realkind,     25497.08526683_realkind,     24871.31302372_realkind,&
            24143.55208429_realkind,     23322.66080848_realkind,     22417.97990382_realkind,&
            21439.22595692_realkind,     20396.44798444_realkind,     19299.92669234_realkind,&
            18160.12861384_realkind,     16987.58390227_realkind,     15792.86806530_realkind,&
            14586.45478044_realkind,     13378.69029188_realkind,     12179.68413384_realkind,&
            10999.25363904_realkind,      9846.77794029_realkind,      8731.20481400_realkind,&
            7660.88648378_realkind,      6643.53270367_realkind,      5686.12140321_realkind,&
            4794.80922574_realkind,      3974.84219208_realkind,      3230.45458673_realkind,&
            2564.80865572_realkind,      1979.88863358_realkind,      1476.41314284_realkind,&
            1053.75138659_realkind,       709.81576221_realkind,       441.03262309_realkind,&
            242.16324537_realkind,       106.30163568_realkind,        24.73731477_realkind,&
            0.00000000_realkind,         0.00000000_realkind,         0.00000000_realkind,&
            0.00000000_realkind/)
       bhalf(1:61)=(/ 0.00000000_realkind,         0.00000000_realkind,         0.00000000_realkind,&
            0.00023670_realkind,         0.00092292_realkind,         0.00224844_realkind,&
            0.00438079_realkind,         0.00746615_realkind,         0.01163013_realkind,&
            0.01697863_realkind,         0.02359864_realkind,         0.03155913_realkind,&
            0.04091181_realkind,         0.05169201_realkind,         0.06391953_realkind,&
            0.07759939_realkind,         0.09272279_realkind,         0.10926775_realkind,&
            0.12720022_realkind,         0.14647459_realkind,         0.16703484_realkind,&
            0.18881506_realkind,         0.21174059_realkind,         0.23572858_realkind,&
            0.26068905_realkind,         0.28652551_realkind,         0.31313601_realkind,&
            0.34041375_realkind,         0.36824815_realkind,         0.39652539_realkind,&
            0.42512962_realkind,         0.45394336_realkind,         0.48284882_realkind,&
            0.51172809_realkind,         0.54046468_realkind,         0.56894389_realkind,&
            0.59705391_realkind,         0.62468618_realkind,         0.65173687_realkind,&
            0.67810731_realkind,         0.70370500_realkind,         0.72844406_realkind,&
            0.75224667_realkind,         0.77504346_realkind,         0.79677456_realkind,&
            0.81739014_realkind,         0.83685163_realkind,         0.85513240_realkind,&
            0.87221834_realkind,         0.88810912_realkind,         0.90281875_realkind,&
            0.91637651_realkind,         0.92882764_realkind,         0.94023440_realkind,&
            0.95067674_realkind,         0.96025322_realkind,         0.96908170_realkind,&
            0.97730036_realkind,         0.98506843_realkind,         0.99256705_realkind,&
            1.000_realkind/)
    case(62)!from ecmwf
       ahalf(1:63)=(/0.000000_realkind,   988.835876_realkind,  1977.676270_realkind,  2966.516602_realkind,&
            3955.356934_realkind,  4944.197266_realkind,  5933.037598_realkind,&
            6921.870117_realkind,  7909.441406_realkind,  8890.707031_realkind,  9860.528320_realkind, &
            10807.783203_realkind, 11722.749023_realkind, 12595.006836_realkind,&
            13419.463867_realkind, 14192.009766_realkind, 14922.685547_realkind, 15638.053711_realkind,&
            16329.560547_realkind, 16990.623047_realkind, 17613.281250_realkind,&
            18191.029297_realkind, 18716.968750_realkind, 19184.544922_realkind, 19587.513672_realkind, &
            19919.796875_realkind, 20175.394531_realkind, 20348.916016_realkind,&
            20434.158203_realkind, 20426.218750_realkind, 20319.011719_realkind, 20107.031250_realkind, &
            19785.357422_realkind, 19348.775391_realkind, 18798.822266_realkind,&
            18141.296875_realkind, 17385.595703_realkind, 16544.585938_realkind, 15633.566406_realkind, &
            14665.645508_realkind, 13653.219727_realkind, 12608.383789_realkind,&
            11543.166992_realkind, 10471.310547_realkind,  9405.222656_realkind,  8356.252930_realkind,  &
            7335.164551_realkind,  6353.920898_realkind,  5422.802734_realkind,&
            4550.215820_realkind,  3743.464355_realkind,  3010.146973_realkind,  2356.202637_realkind, &
            1784.854614_realkind,  1297.656128_realkind,   895.193542_realkind,&
            576.314148_realkind,   336.772369_realkind,   162.043427_realkind,    54.208336_realkind, &
            6.575628_realkind,     0.003160_realkind,     0.000000_realkind/)
       bhalf(1:63)=(/0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000013_realkind,&
            0.000087_realkind,&
            0.000275_realkind,0.000685_realkind,0.001415_realkind,0.002565_realkind,0.004187_realkind,&
            0.006322_realkind,0.009035_realkind,0.012508_realkind,0.016860_realkind,0.022189_realkind,&
            0.028610_realkind,&
            0.036227_realkind,0.045146_realkind,0.055474_realkind,0.067316_realkind,0.080777_realkind,&
            0.095964_realkind,0.112979_realkind,0.131935_realkind,0.152934_realkind,0.176091_realkind,&
            0.201520_realkind,&
            0.229315_realkind,0.259554_realkind,0.291993_realkind,0.326329_realkind,0.362203_realkind,&
            0.399205_realkind,0.436906_realkind,0.475016_realkind,0.513280_realkind,0.551458_realkind,&
            0.589317_realkind,&
            0.626559_realkind,0.662934_realkind,0.698224_realkind,0.732224_realkind,0.764679_realkind,&
            0.795385_realkind,0.824185_realkind,0.850950_realkind,0.875518_realkind,0.897767_realkind,&
            0.917651_realkind,&
            0.935157_realkind,0.950274_realkind,0.963007_realkind,0.973466_realkind,0.982238_realkind,&
            0.989153_realkind,0.994204_realkind,0.997630_realkind,1.000000_realkind/)
    case(91)!from ecmwf
       ahalf(1:92)=(/ 0.000000_realkind,  2.000040_realkind,  3.980832_realkind, 7.387186_realkind, &
            12.908319_realkind, 21.413612_realkind, & 
            33.952858_realkind, 51.746601_realkind,  76.167656_realkind, 108.715561_realkind, &
            150.986023_realkind,204.637451_realkind,271.356506_realkind,& 
            352.824493_realkind, 450.685791_realkind, 566.519226_realkind, 701.813354_realkind, &
            857.945801_realkind, 1036.166504_realkind, 1237.585449_realkind,& 
            1463.163940_realkind, 1713.709595_realkind, 1989.874390_realkind, 2292.155518_realkind, &
            2620.898438_realkind, 2976.302246_realkind, 3358.425781_realkind, &
            3767.196045_realkind, 4202.416504_realkind, 4663.776367_realkind, 5150.859863_realkind, &
            5663.156250_realkind, 6199.839355_realkind, &
            6759.727051_realkind, 7341.469727_realkind, 7942.926270_realkind, 8564.624023_realkind, &
            9208.305664_realkind, 9873.560547_realkind, &
            10558.881836_realkind, 11262.484375_realkind, 11982.662109_realkind, 12713.897461_realkind, &
            13453.225586_realkind, 14192.009766_realkind, &
            14922.685547_realkind, 15638.053711_realkind, 16329.560547_realkind, 16990.623047_realkind,&
            17613.281250_realkind, 18191.029297_realkind,&
            18716.968750_realkind, 19184.544922_realkind, 19587.513672_realkind, 19919.796875_realkind, &
            20175.394531_realkind, 20348.916016_realkind, &
            20434.158203_realkind, 20426.218750_realkind, 20319.011719_realkind, 20107.031250_realkind, &
            19785.357422_realkind, 19348.775391_realkind, 18798.822266_realkind,&
            18141.296875_realkind, 17385.595703_realkind, 16544.585938_realkind, 15633.566406_realkind,&
            14665.645508_realkind, 13653.219727_realkind, 12608.383789_realkind, &
            11543.166992_realkind, 10471.310547_realkind,  9405.222656_realkind,  8356.252930_realkind,&
            7335.164551_realkind,  6353.920898_realkind,  5422.802734_realkind,&
            4550.215820_realkind,  3743.464355_realkind,  3010.146973_realkind,  2356.202637_realkind,&
            1784.854614_realkind,  1297.656128_realkind,   895.193542_realkind,&
            576.314148_realkind,   336.772369_realkind,   162.043427_realkind,    54.208336_realkind,&
            6.575628_realkind,     0.003160_realkind,     0.000000_realkind /)

       bhalf(1:92)=(/0.000000_realkind, 0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,&
            0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000000_realkind,0.000014_realkind,&
            0.000055_realkind,0.000131_realkind,0.000279_realkind,0.000548_realkind,0.001000_realkind,0.001701_realkind,&
            0.002765_realkind,0.004267_realkind,0.006322_realkind,0.009035_realkind,0.012508_realkind,&
            0.016860_realkind,0.022189_realkind,0.028610_realkind,0.036227_realkind,0.045146_realkind,0.055474_realkind,&
            0.067316_realkind,0.080777_realkind,0.095964_realkind,0.112979_realkind,0.131935_realkind,&
            0.152934_realkind,0.176091_realkind,0.201520_realkind,0.229315_realkind,0.259554_realkind,0.291993_realkind,&
            0.326329_realkind,0.362203_realkind,0.399205_realkind,0.436906_realkind,0.475016_realkind,&
            0.513280_realkind,0.551458_realkind,0.589317_realkind,0.626559_realkind,0.662934_realkind,0.698224_realkind,&
            0.732224_realkind,0.764679_realkind,0.795385_realkind,0.824185_realkind,0.850950_realkind,&
            0.875518_realkind,0.897767_realkind,0.917651_realkind,0.935157_realkind,0.950274_realkind,0.963007_realkind,&
            0.973466_realkind,0.982238_realkind,0.989153_realkind,0.994204_realkind,0.997630_realkind,1.000000_realkind/)
    case default
       write(*,*)klev
       call stop_program( 'ahalf,bhalf not available for klev=')
    end select
  end subroutine get_a_b


end module rcaDomainMod
