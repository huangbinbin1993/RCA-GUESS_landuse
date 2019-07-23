module mod_diffh
  use RCAdomainMod
  use decomp
  implicit none
  private

  real(kind=realkind),save::rimpdt = - 99999._realkind 

  integer,save::ndifu=6 ! horizontal velocity component u
  integer,save::ndifv=6 ! horizontal velocity component v
  integer,save::ndift=6 ! temperature
  integer,save::ndifq=6 ! specific humidity
  integer,save::ndifs=6 ! specific cloud water
  integer,save,allocatable::ndifx(:) ! extra scalar prognostic variables

  real(kind=realkind),allocatable,dimension(:),save::cdifu 
  real(kind=realkind),allocatable,dimension(:),save::cdifv                                
  real(kind=realkind),allocatable,dimension(:),save::cdift
  real(kind=realkind),allocatable,dimension(:),save::cdifq
  real(kind=realkind),allocatable,dimension(:),save::cdifs 
  real(kind=realkind),allocatable,dimension(:,:),save::cdifx

  logical,save::ldifu=.true.!  .true. if impl. hor. diffusion for u
  logical,save::ldifv=.true.!  .true. if impl. hor. diffusion for v
  logical,public,save::ldift=.true.!  .true. if impl. hor. diffusion for t
  logical,save::ldifq=.true.!  .true. if impl. hor. diffusion for q
  logical,save::ldifs=.true.!  .true. if impl. hor. diffusion for cloud water
  logical,allocatable,dimension(:),save::ldifx ! .true. if impl. hor. diffusion for extra passive scalars


  logical,public,save::nlhdif=.false.! .true. if explicit hor. diffusion in dynamics
  real(kind=realkind),public,allocatable,dimension(:),save::atcref! list of coefficients for temperature along pseudo pressure levels
  real(kind=realkind),allocatable,dimension(:),save::ak4lev
  real(kind=realkind),save::ak4=1.0e+14_realkind !diffusion coefficient for 2nd order hor. diffu computed statistics

  real(kind=realkind),save::rprers        ! reference pressure ?
  real(kind=realkind),save::rtemrs        ! reference temperature trs
  real(kind=realkind),save::rtemrt        ! reference temperature trt

  !     implicit diffusion coefficient files
  real(kind=realkind),allocatable,dimension(:,:,:)::diffcx,diffcy 
  integer,save::ndiffx,ndiffy
  public diffh,difhini,hdiff4
contains

  subroutine read_namdiffh(klev,ksvar)

    integer,intent(in)::klev,ksvar
    integer,parameter::maxLevels=91,nc=2!!
    real(kind=realkind)::cdifuin(maxLevels)=1.0_realkind
    real(kind=realkind)::cdifvin(maxLevels)=1.0_realkind
    real(kind=realkind)::cdiftin(maxLevels)=1.0_realkind
    real(kind=realkind)::cdifqin(maxLevels)=1.0_realkind
    real(kind=realkind)::cdifsin(maxLevels)=1.0_realkind
    real(kind=realkind)::cdifxin(maxLevels,nc)=1.0_realkind
    real(kind=realkind)::ak4levin(maxLevels*(nc+5))
    integer::ndifxin(nc)=6
    logical::ldifxin(nc)=.true.
    integer i
    namelist/namdiffh/ndifu,ndifv,ndift,ndifq,ndifs,ndifxin,&
         cdifuin,cdifvin,cdiftin,cdifqin,cdifsin,cdifxin,&
         ldifu,ldifv,ldift,ldifq,ldifs,ldifxin,nlhdif,ak4levin,ak4
    do i=1,4
       ak4levin(((i-1)*klev+1):i*klev) = 1.0_realkind
    enddo
    cdifu = 0._realkind
    cdifv = 0._realkind
    cdift = 0._realkind
    cdifq = 0._realkind
    cdifs = 0._realkind
    cdifx = 0._realkind
    ldifx = .false.
    open(57,file='namelists.dat',status='old')
    read(57,nml=namdiffh)
    close(57)
    if(mype==0)then
       write(6,nml=namdiffh)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namdiffh)
       close(1)
    endif
    do i=1,ksvar
       cdifx(:,1) = cdifxin(1:klev,i)
    enddo
    ak4lev = ak4levin(1:klev*(5+ksvar)) 
    ndifx = ndifxin(1:ksvar)
    cdifu = cdifuin(1:klev)
    cdifv = cdifvin(1:klev)
    cdift = cdiftin(1:klev)
    cdifq = cdifqin(1:klev)
    cdifs = cdifsin(1:klev)
  end subroutine read_namdiffh


  subroutine hdiff4(klon,klat,klev,ksvar,nltcrf,ptwodt ,&
       plnpsm,ptm,pum,pvm,pqm,psm,psvarm,&
       ptp,pup,pvp,pqp,psp,psvarp)

    !     hdiff4 - horizontal diffusion
    !     
    !     j.e. haugen               hirlam
    !     k.s. eerola               hirlam4(revised)
    !     hdiff4 applies horizontal diffusion on t,u,v,q,s
    !     and extra scalars
    !     
    !     input parameters:
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels
    !     ptwodt    current double timestep
    !     plnpsm    ln( surface pressure )
    !     ptm       temperature
    !     pum       velocity-component in x-direction
    !     pvm       velocity-component in x-direction
    !     pqm       spesific humidity
    !     psm       passive scalar(liquid water)
    !     psvarm    extra scalars
    !     ksvar     number of extra scalars
    !     
    !     ptp       updated temperature
    !     pup       updated velocity-component in x-direction
    !     pvp       updated velocity-component in x-direction
    !     pqp       updated spesific humidity
    !     psp       updated passive scalar
    !     psvarp    updated extra scalars
    !     
    !     method.
    !     -------
    !     
    !     method is linear fourth order diffusion
    !     diffusion of t is performed along reference pressure surfaces
    !     defined in subroutine difibni
    !     
    !     diffusion of q along reference pressure surfaces added
    !     diffusion of s and extra scalars added

    implicit none

    integer  klon   ,klat   ,klev   ,ksvar
    logical  nltcrf
    real(kind=realkind)     ptwodt !,ak4    
    real(kind=realkind) plnpsm(klon,klat)     ,                   &
         ptm(klon,klat,klev),   pum(klon,klat,klev),&
         pvm(klon,klat,klev),   pqm(klon,klat,klev),&
         psm(klon,klat,klev),                       &
         ptp(klon,klat,klev),   pup(klon,klat,klev),&
         pvp(klon,klat,klev),   pqp(klon,klat,klev),&
         psp(klon,klat,klev)                         &
         ,psvarm(klon,klat,klev,*)                  &
         ,psvarp(klon,klat,klev,*)

    !     declaration of local workspace

    integer jx,jy,jk,l
    integer istart,istop,jstart,jstop,istartp1,jstartp1,istopm1,jstopm1

    real(kind=realkind) zrxtry(klon,klat),zytrx(klon,klat),zxtry(klon,klat),&
         zrxby(klon,klat),zrxbry(klon,klat),                    &
         zxby(klon,klat),zrybx(klon,klat),zybx(klon,klat), &
         zrybrx(klon,klat),ztemp(klon,klat),ztempt(klon,klat), &
         ztempu(klon,klat),ztempv(klon,klat),ztempq(klon,klat), &
         ztemps(klon,klat),zlaps(klon,klat)

    real(kind=realkind) zfactr,zra2,zrdx2,zrdy2,zfactt,zfacuu ,zfacvv,zfacqq &
         ,zfacss ,ztcref,zqcref


    if(atright) then
       istop = klon -1
    else
       istop = klon
    endif
    if(attop) then
       jstop = klat -1
    else
       jstop = klat
    endif
    istart   = 1
    jstart   = 1
    jstartp1 = jstart + 1
    istartp1 = istart + 1
    jstopm1  = jstop  - 1
    istopm1  = istop  - 1

    zfactr = ptwodt*ak4

    !     prepare mapping factors

    zra2  = ra*ra
    zrdx2 = rdlam*rdlam
    zrdy2 = rdth*rdth

    do jy = 1,klat
       do jx = 1,klon
          zrxtry(jx,jy) = rhxu(jx,jy) * rhyv(jx,jy) * zra2
          zytrx(jx,jy) =  hyu(jx,jy) * rhxu(jx,jy) * zrdx2
          zxtry(jx,jy) =  hxv(jx,jy) * rhyv(jx,jy) * zrdy2

          zrxby(jx,jy) = rhxu(jx,jy) /  hyu(jx,jy) * zra2
          zrxbry(jx,jy) = rhxu(jx,jy) / rhyv(jx,jy) * zrdx2
          zxby(jx,jy) =  hxv(jx,jy) /  hyu(jx,jy) * zrdy2

          zrybx(jx,jy) = rhyv(jx,jy) /  hxv(jx,jy) * zra2
          zybx(jx,jy) =  hyu(jx,jy) /  hxv(jx,jy) * zrdx2
          zrybrx(jx,jy) = rhyv(jx,jy) / rhxu(jx,jy) * zrdy2

       enddo
    enddo


    !     correction term for temperature along pressure surfaces
    if(nltcrf) then
       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1
             ztemp(jx,jy) = zrxtry(jx,jy) *                          &
                (zytrx(jx  ,jy) *(plnpsm(jx+1,jy)-plnpsm(jx,jy))   &
                  +  zytrx(jx-1,jy) *(plnpsm(jx-1,jy)-plnpsm(jx,jy)) &
                  +  zxtry(jx,jy  ) *(plnpsm(jx,jy+1)-plnpsm(jx,jy)) &
                  +  zxtry(jx,jy-1) *(plnpsm(jx,jy-1)-plnpsm(jx,jy)))
          enddo
       enddo

       if(atbase) then
          do jx=1,klon
             ztemp(jx,1) = 0.0_realkind
          enddo
       endif

       if(attop) then
          do jx=1,klon
             ztemp(jx,klat)   = 0.0_realkind
             ztemp(jx,klat-1) = 0.0_realkind
          enddo
       endif

       if(atleft) then
          do jy=1,klat
             ztemp(1,jy) = 0.0_realkind
          enddo
       endif
       if(atright) then
          do jy=1,klat
             ztemp(klon  ,jy) = 0.0_realkind
             ztemp(klon-1,jy) = 0.0_realkind
          enddo
       endif

       call swap2d(ztemp,klon,klat)

       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1
             zlaps(jx,jy) = zrxtry(jx,jy) *                         &
                (zytrx(jx  ,jy) *(ztemp(jx+1,jy)-ztemp(jx,jy))    &
                  +  zytrx(jx-1,jy) *(ztemp(jx-1,jy)-ztemp(jx,jy))  &
                  +  zxtry(jx,jy  ) *(ztemp(jx,jy+1)-ztemp(jx,jy))  &
                  +  zxtry(jx,jy-1) *(ztemp(jx,jy-1)-ztemp(jx,jy)))
          enddo
       enddo
       call swap2d(zlaps,klon,klat)
    else
       do jy = 1,klat
          do jx = 1,klon
             zlaps(jx,jy) = 0.0_realkind
          enddo
       enddo
    endif

    !     loop over all levels
    do jk=1,klev
       zfactt = zfactr * ak4lev(       jk)
       zfacuu = zfactr * ak4lev(  klev+jk)
       zfacvv = zfactr * ak4lev(2*klev+jk)
       zfacqq = zfactr * ak4lev(3*klev+jk)
       zfacss = zfactr * ak4lev(4*klev+jk)
       ztcref = atcref(jk)
       zqcref = atcref(jk)*5.5e-2_realkind

       !     compute laplacian

       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1

             ztempt(jx,jy) = zrxtry(jx,jy) *                          &
                ( zytrx(jx  ,jy  )*( ptm(jx+1,jy,jk)-ptm(jx,jy,jk))  &
                  + zytrx(jx-1,jy  )*( ptm(jx-1,jy,jk)-ptm(jx,jy,jk))  &
                  + zxtry(jx  ,jy  )*( ptm(jx,jy+1,jk)-ptm(jx,jy,jk))  &
                  + zxtry(jx  ,jy-1)*( ptm(jx,jy-1,jk)-ptm(jx,jy,jk)))

             ztempq(jx,jy) = zrxtry(jx,jy) *                          &
                (  zytrx(jx  ,jy)*( pqm(jx+1,jy,jk)-pqm(jx,jy,jk))   &
                  +  zytrx(jx-1,jy)*( pqm(jx-1,jy,jk)-pqm(jx,jy,jk))   &
                  +  zxtry(jx,jy  )*( pqm(jx,jy+1,jk)-pqm(jx,jy,jk))   &
                  +  zxtry(jx,jy-1)*( pqm(jx,jy-1,jk)-pqm(jx,jy,jk)))

             ztemps(jx,jy) = zrxtry(jx,jy) *                          &
                (  zytrx(jx  ,jy)*(psm(jx+1,jy,jk)-psm(jx,jy,jk))    &
                  +  zytrx(jx-1,jy)*(psm(jx-1,jy,jk)-psm(jx,jy,jk))    &
                  +  zxtry(jx,jy  )*(psm(jx,jy+1,jk)-psm(jx,jy,jk))    &
                  +  zxtry(jx,jy-1)*(psm(jx,jy-1,jk)-psm(jx,jy,jk)))

             ztempu(jx,jy) = zrxby(jx,jy) *                          &
                ( zrxbry(jx+1,jy) *(pum(jx+1,jy,jk)-pum(jx,jy,jk))  &
                  + zrxbry(jx  ,jy) *(pum(jx-1,jy,jk)-pum(jx,jy,jk))  &
                  +   zxby(jx,jy  ) *(pum(jx,jy+1,jk)-pum(jx,jy,jk))  &
                  +   zxby(jx,jy-1) *(pum(jx,jy-1,jk)-pum(jx,jy,jk)))

             ztempv(jx,jy) = zrybx(jx,jy) *                          &
                (   zybx(jx  ,jy) *(pvm(jx+1,jy,jk)-pvm(jx,jy,jk))  &
                  +   zybx(jx-1,jy) *(pvm(jx-1,jy,jk)-pvm(jx,jy,jk))  &
                  + zrybrx(jx,jy+1) *(pvm(jx,jy+1,jk)-pvm(jx,jy,jk))  &
                  + zrybrx(jx,jy  ) *(pvm(jx,jy-1,jk)-pvm(jx,jy,jk)))

          enddo
       enddo

       call swap2d(ztempt,klon,klat)
       call swap2d(ztempq,klon,klat)
       call swap2d(ztempu,klon,klat)
       call swap2d(ztempv,klon,klat)
       call swap2d(ztemps,klon,klat)
!       call swap5(ztempt,ztempq,ztempu,ztempv,ztemps,klon,klat,1)

       !     boundary values of laplacian assumed equal zero

       if(atbase) then
          do jx = istart,istop
             ztempt(jx,1) = 0.0_realkind
             ztempq(jx,1) = 0.0_realkind
             ztemps(jx,1) = 0.0_realkind
             ztempu(jx,1) = 0.0_realkind
             ztempv(jx,1) = 0.0_realkind
          enddo
       endif

       if(attop) then
          do jx = istart,istop
             ztempt(jx,klat) = 0.0_realkind
             ztempq(jx,klat) = 0.0_realkind
             ztemps(jx,klat) = 0.0_realkind
             ztempu(jx,klat) = 0.0_realkind
             ztempv(jx,klat) = 0.0_realkind

             ztempt(jx,klat-1) = 0.0_realkind
             ztempq(jx,klat-1) = 0.0_realkind
             ztemps(jx,klat-1) = 0.0_realkind
             ztempu(jx,klat-1) = 0.0_realkind
             ztempv(jx,klat-1) = 0.0_realkind
          enddo
       endif

       if(atleft) then
          do  jy= jstart,jstop
             ztempt(1,jy) = 0.0_realkind
             ztempq(1,jy) = 0.0_realkind
             ztemps(1,jy) = 0.0_realkind
             ztempu(1,jy) = 0.0_realkind
             ztempv(1,jy) = 0.0_realkind
          enddo
       endif

       if(atright) then
          do  jy= jstart,jstop
             ztempt(klon,jy) = 0.0_realkind
             ztempq(klon,jy) = 0.0_realkind
             ztemps(klon,jy) = 0.0_realkind
             ztempu(klon,jy) = 0.0_realkind
             ztempv(klon,jy) = 0.0_realkind

             ztempt(klon-1,jy) = 0.0_realkind
             ztempq(klon-1,jy) = 0.0_realkind
             ztemps(klon-1,jy) = 0.0_realkind
             ztempu(klon-1,jy) = 0.0_realkind
             ztempv(klon-1,jy) = 0.0_realkind
          enddo
       endif

       !     compute laplacian of laplacian and add changes


       do jy = jstartp1,jstopm1
          do jx = istartp1,istopm1

             ptp(jx,jy,jk) = ptp(jx,jy,jk) - zfactt *(                 &
                  zrxtry(jx,jy) *                                        &
                (  zytrx(jx  ,jy) *(ztempt(jx+1,jy)-ztempt(jx,jy))    &
                  +  zytrx(jx-1,jy) *(ztempt(jx-1,jy)-ztempt(jx,jy))    &
                  +  zxtry(jx,jy  ) *(ztempt(jx,jy+1)-ztempt(jx,jy))    &
                  +  zxtry(jx,jy-1) *(ztempt(jx,jy-1)-ztempt(jx,jy)))   &
                  - ztcref * zlaps(jx,jy) )                                   

             pqp(jx,jy,jk) = pqp(jx,jy,jk) - zfacqq *(                  &
                  zrxtry(jx,jy) *                                        &
                (  zytrx(jx  ,jy) *(ztempq(jx+1,jy)-ztempq(jx,jy))    &
                  +  zytrx(jx-1,jy) *(ztempq(jx-1,jy)-ztempq(jx,jy))    &
                  +  zxtry(jx,jy  ) *(ztempq(jx,jy+1)-ztempq(jx,jy))    &
                  +  zxtry(jx,jy-1) *(ztempq(jx,jy-1)-ztempq(jx,jy)))   &
                  - zqcref * zlaps(jx,jy) * pqp(jx,jy,jk))                


             psp(jx,jy,jk) = psp(jx,jy,jk) - zfacss *                    &
                  zrxtry(jx,jy) *                                        &
                (  zytrx(jx  ,jy) *(ztemps(jx+1,jy)-ztemps(jx,jy))    &
                  +  zytrx(jx-1,jy) *(ztemps(jx-1,jy)-ztemps(jx,jy))    &
                  +  zxtry(jx,jy  ) *(ztemps(jx,jy+1)-ztemps(jx,jy))    &
                  +  zxtry(jx,jy-1) *(ztemps(jx,jy-1)-ztemps(jx,jy)))        

             pup(jx,jy,jk) = pup(jx,jy,jk) - zfacuu *                    &
                  zrxby(jx,jy) *                                        &
                ( zrxbry(jx+1,jy) *(ztempu(jx+1,jy)-ztempu(jx,jy))    &
                  + zrxbry(jx  ,jy) *(ztempu(jx-1,jy)-ztempu(jx,jy))    &
                  +   zxby(jx,jy  ) *(ztempu(jx,jy+1)-ztempu(jx,jy))    &
                  +   zxby(jx,jy-1) *(ztempu(jx,jy-1)-ztempu(jx,jy)))    

             pvp(jx,jy,jk) = pvp(jx,jy,jk) - zfacvv *                    &
                  zrybx(jx,jy) *                                        &
                (   zybx(jx  ,jy) *(ztempv(jx+1,jy)-ztempv(jx,jy))    &
                  +   zybx(jx-1,jy) *(ztempv(jx-1,jy)-ztempv(jx,jy))    &
                  + zrybrx(jx,jy+1) *(ztempv(jx,jy+1)-ztempv(jx,jy))    &
                  + zrybrx(jx,jy  ) *(ztempv(jx,jy-1)-ztempv(jx,jy)))

          enddo
       enddo
    enddo


    !     diffusion of extra scalars
    do l=1,ksvar
       do jk=1,klev
          zfacss = zfactr*ak4lev((4+l)*klev+jk)

          if(abs(zfacss)>1.e-14_realkind) then
             do jy = jstartp1,jstopm1
                do jx = istartp1,istopm1
                   ztemps(jx,jy) = zrxtry(jx,jy)*                    &
                      ( zytrx(jx  ,jy  )                            &
                        *(psvarm(jx+1,jy  ,jk,l)-psvarm(jx,jy,jk,l)) &
                        + zytrx(jx-1,jy  )                            &
                        *(psvarm(jx-1,jy  ,jk,l)-psvarm(jx,jy,jk,l)) &
                        + zxtry(jx  ,jy  )                            &
                        *(psvarm(jx  ,jy+1,jk,l)-psvarm(jx,jy,jk,l)) &
                        + zxtry(jx  ,jy-1)                            &
                        *(psvarm(jx  ,jy-1,jk,l)-psvarm(jx,jy,jk,l)))
                enddo
             enddo

             call swap2d(ztemps,klon,klat)

             if(atleft)then
                do jx = istart,istop
                   ztemps(jx,  1) = 0.0_realkind
                enddo
             endif

             if(atright)then
                do jx = istart,istop
                   ztemps(jx,klat) = 0.0_realkind
                   ztemps(jx,klat-1) = 0.0_realkind
                enddo
             endif

             if(atbase) then
                do  jy= jstartp1,jstopm1
                   ztemps(   1,jy) = 0.0_realkind
                enddo
             endif
             if(attop) then
                do  jy= jstartp1,jstopm1
                   ztemps(klon,jy) = 0.0_realkind
                   ztemps(klon-1,jy) = 0.0_realkind
                enddo
             endif

             do jy = jstartp1,jstopm1
                do jx = istartp1,istopm1
                   psvarp(jx,jy,jk,l) = psvarp(jx,jy,jk,l)            &
                        - zfacss*zrxtry(jx,jy)*                        &
                      ( zytrx(jx,jy)*(ztemps(jx+1,jy)-ztemps(jx,jy)) &
                        +zytrx(jx-1,jy)*(ztemps(jx-1,jy)-ztemps(jx,jy)) &
                        +zxtry(jx,jy)*(ztemps(jx,jy+1)-ztemps(jx,jy))   &
                        +zxtry(jx,jy-1)*(ztemps(jx,jy-1)-ztemps(jx,jy)))
                enddo
             enddo
          endif
       enddo
    enddo

    return
  end subroutine hdiff4

  subroutine difhini(kmode,klon,klat,klev,dt,kpbpts,ksvar,ahyb,bhyb)
    !     Initialize the implicit horizontal diffusion
    !     DESCRIPTION:
    !     KMODE = 1 : Rescales the horizontal diffusion coefficients
    !     KMODE != 1 : Calculates the coefficients used in the implicit
    !     horizontal diffusion.
    !     KALLE EEROLA   HIRLAM    1998
    use transpose
    implicit none  
    
    integer,intent(in) ::klon,klat,klev,ksvar,kmode,kpbpts 
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind) :: dt 
    real(kind=realkind),allocatable,dimension(:,:):: div_fft,div_tri  
    integer :: jx,jy,jk,ipara,j,i,jrec_fft,jstart,jslab,&
         jrec_tri,jrec,jl
    logical :: ldiffi  
    real(kind=realkind) :: zxu,zxv,zxt,zxq,zxs,zxx,zyu,zyv,zyt,zyq,zys,&
         zyx,rbdeps,zcdif,dtheta
    real(kind=realkind),allocatable,dimension(:,:,:) ::cx,cy 
    integer::status
    !     THE FOLLOWING ARE THE DEFAULT VALUES FOR THE HORIZONTAL DIFFUSION
    !     COEFFICIENTS ON THE 0.5X0.5 DEGREE GRID.
    real(kind=realkind),parameter:: ppdif52 = 5.0e+5_realkind
    real(kind=realkind),parameter:: ppdif54 = 3.5e+14_realkind  
    real(kind=realkind),parameter:: ppdif56 = 1.0e+24_realkind  
    real(kind=realkind),parameter:: ppdth5 = 0.5_realkind  
    !     KMODE=1: RE-SCALE THE HORIZONTAL DIFFUSUION COEFFICIENTS.
    !     K1 = K2 * {(DX1/DX2)**N }.
    !     N = ORDER OF DIFFUSION
    !     THIS SCALES THE E-FOLDING TIME APPROXIMATELY.
    !     THE COEFFICIENTS ARE INCREASED FOR TOP 5 LAYERS.
    integer::status1,status2

    allocate(cx(klon,klat,klev),cy(klon,klat,klev)  )
    ldiffi = .false.  

    if(.not.allocated(cdifu))then
       ndiffx = klat_global*klev/nproc + 1
       ndiffy = klon_global*klev/nproc + 1
       allocate(cdifx(klev,ksvar))
       allocate(ndifx(ksvar))
       allocate(cdifu(klev),cdifv(klev),cdift(klev))
       allocate(cdifq(klev),cdifs(klev))
       allocate(ldifx(ksvar))
       allocate(diffcx(klon_global,ndiffx,5+ksvar),stat=status1)
       allocate(diffcy(klat_global,ndiffy,5+ksvar),stat=status2)
       if(status1+status2/=0)then
          stop 'memory allocation error diff_h'
       endif
       allocate(ak4lev(klev*(5+ksvar)))
       allocate(atcref(klev))
       call read_namdiffh(klev,ksvar)
       if(nlhdif .or. ldift) then
          call difini(klev,ksvar,ahyb,bhyb)
       endif
    endif


    if(kmode==1) then  
       dtheta = 1._realkind /(rdth * 2._realkind * asin(1._realkind) / 180._realkind)  
       if(cdifu(1) >1.e8_realkind) then  
          if(mype==0) print * ,'cdifu(1) is probably too large.'  
          if(mype==0) print  * ,'did you use an obsolete namelist?'  
          if(mype==0) print * ,'cdifu/v/t/q should be order 1'  
          if(mype==0) print  * ,'see release notes for version 4.1.1'  
          call stop_program('')  
       endif
       if(ldifu) then  
          ldiffi = .true.  
          if(ndifu==2) zcdif = ppdif52 *((dtheta / ppdth5) **2.0_realkind)  
          if(ndifu==4) zcdif = ppdif54 *((dtheta / ppdth5) **4.0_realkind)  
          if(ndifu==6) zcdif = ppdif56 *((dtheta / ppdth5) **6.0_realkind)  
          cdifu(1) = cdifu(1) * zcdif * 16.0_realkind  
          cdifu(2) = cdifu(2) * zcdif * 16.0_realkind  
          cdifu(3) = cdifu(3) * zcdif * 8.0_realkind  
          cdifu(4) = cdifu(4) * zcdif * 4.0_realkind  
          cdifu(5) = cdifu(5) * zcdif * 2.0_realkind  
          do jk = 6,klev  
             cdifu(jk) = cdifu(jk) * zcdif  
          enddo
       endif
       !
       if(ldifv) then  
          ldiffi = .true.  
          if(ndifv==2) zcdif = ppdif52 *((dtheta / ppdth5) **2.0_realkind)  
          if(ndifv==4) zcdif = ppdif54 *((dtheta / ppdth5) **4.0_realkind)  
          if(ndifv==6) zcdif = ppdif56 *((dtheta / ppdth5) **6.0_realkind)  
          cdifv(1) = cdifv(1) * zcdif * 16.0_realkind 
          cdifv(2) = cdifv(2) * zcdif * 16.0_realkind  
          cdifv(3) = cdifv(3) * zcdif * 8.0_realkind  
          cdifv(4) = cdifv(4) * zcdif * 4.0_realkind  
          cdifv(5) = cdifv(5) * zcdif * 2.0_realkind  
          do jk = 6,klev  
             cdifv(jk) = cdifv(jk) * zcdif  
          enddo
       endif
       !
       if(ldift) then  
          ldiffi = .true.  
          if(ndift==2) zcdif = ppdif52 *((dtheta / ppdth5) **2.0_realkind)  
          if(ndift==4) zcdif = ppdif54 *((dtheta / ppdth5) **4.0_realkind)  
          if(ndift==6) zcdif = ppdif56 *((dtheta / ppdth5) **6.0_realkind)  
          cdift(1) = cdift(1) * zcdif * 16.0_realkind  
          cdift(2) = cdift(2) * zcdif * 16.0_realkind  
          cdift(3) = cdift(3) * zcdif * 8.0_realkind  
          cdift(4) = cdift(4) * zcdif * 4.0_realkind  
          cdift(5) = cdift(5) * zcdif * 2.0_realkind  
          do jk = 6,klev  
             cdift(jk) = cdift(jk) * zcdif  
          enddo
       endif
       !
       if(ldifq) then  
          ldiffi = .true.  
          if(ndifq==2) zcdif = ppdif52 *((dtheta / ppdth5) **2.0_realkind)  
          if(ndifq==4) zcdif = ppdif54 *((dtheta / ppdth5) **4.0_realkind)  
          if(ndifq==6) zcdif = ppdif56 *((dtheta / ppdth5) **6.0_realkind)  
          cdifq(1) = cdifq(1) * zcdif * 16.0_realkind  
          cdifq(2) = cdifq(2) * zcdif * 16.0_realkind  
          cdifq(3) = cdifq(3) * zcdif * 8.0_realkind  
          cdifq(4) = cdifq(4) * zcdif * 4.0_realkind  
          cdifq(5) = cdifq(5) * zcdif * 2.0_realkind  
          do jk = 6,klev  
             cdifq(jk) = cdifq(jk) * zcdif  
          enddo
       endif
       !
       if(ldifs) then  
          ldiffi = .true.  
          if(ndifs==2) zcdif = ppdif52 *((dtheta / ppdth5) **2.0_realkind)  
          if(ndifs==4) zcdif = ppdif54 *((dtheta / ppdth5) **4.0_realkind)  
          if(ndifs==6) zcdif = ppdif56 *((dtheta / ppdth5) **6.0_realkind)  
          cdifs(1) = cdifs(1) * zcdif * 16.0_realkind  
          cdifs(2) = cdifs(2) * zcdif * 16.0_realkind  
          cdifs(3) = cdifs(3) * zcdif * 8.0_realkind 
          cdifs(4) = cdifs(4) * zcdif * 4.0_realkind  
          cdifs(5) = cdifs(5) * zcdif * 2.0_realkind  
          do jk = 6,klev  
             cdifs(jk) = cdifs(jk) * zcdif  
          enddo
       endif
       do jl = 1,ksvar  
          if(ldifx(jl) ) then  
             ldiffi = .true.  
             if(ndifx(jl) ==2) zcdif = ppdif52 *((dtheta / ppdth5)**2.0_realkind)
             if(ndifx(jl) ==4) zcdif = ppdif54 *((dtheta / ppdth5)**4.0_realkind)
             if(ndifx(jl) ==6) zcdif = ppdif56 *((dtheta / ppdth5)**6.0_realkind)
             cdifx(1,jl) = cdifx(1,jl) * zcdif * 16.0_realkind  
             cdifx(2,jl) = cdifx(2,jl) * zcdif * 16.0_realkind  
             cdifx(3,jl) = cdifx(3,jl) * zcdif * 8.0_realkind  
             cdifx(4,jl) = cdifx(4,jl) * zcdif * 4.0_realkind  
             cdifx(5,jl) = cdifx(5,jl) * zcdif * 2.0_realkind  
             do jk = 6,klev  
                cdifx(jk,jl) = cdifx(jk,jl) * zcdif  
             enddo
          endif
       enddo

       if(mype==0.and.ldiffi) then  
          write(6,* )  
          write(6, * ) ' the values of the horiz diffusion coefficients are'
          write(6,'(a)') ' lev    u         v         t         q    s' &
               //'         tke'
          write(6,'(i3,7e10.3)') 1,cdifu(1) ,cdifv(1) ,cdift( &
               1) ,cdifq(1) ,cdifs(1) ,cdifx(1,1)
          do jk = 2,klev  
             if(abs(cdifu(jk)- cdifu(jk - 1))>1.e-14_realkind     .or. &
                  abs(cdifv(jk) - cdifv(jk - 1))> 1.e-14_realkind .or.&
                  abs(cdift(jk) - cdift(jk - 1))> 1.e-14_realkind.or. &
                  abs(cdifq(jk) - cdifq(jk - 1))> 1.e-14_realkind.or. &
                  abs(cdifs(jk) - cdifs(jk - 1))> 1.e-14_realkind.or. &
                  abs(cdifx(jk,1) - cdifx(jk - 1,1))>1.e-14_realkind )&
                  write(6,'(i3,7e10.3)') jk,cdifu(jk) ,cdifv(jk) &
                  & ,cdift(jk) ,cdifq(jk) ,cdifs(jk) ,cdifx(jk,1)
          enddo
          write(6, * ) '(level not mentioned?: same as the level above)'  
          write(6,* ) 'consult hirlam htlm documentation for '// &
               'conversion of these coefficients to e-folding times'
          !

       endif
    endif
    !
    !     KMODE#1: CALCLULATE THE CONTSNTS USED IN THE DIFFUSION SCHEME
    !     --------------------------------------------------------------
    !     The constants are calculated in the beginning and if the
    !     timestep is changing. Normally in the beginning and after the
    !     first timestep.
    !

    if(kmode/=1.and.abs(dt-rimpdt)>1.e-14_realkind) then  
       allocate(div_fft(klon_global,klat_global*klev_global ))
       allocate(div_tri(klat_global,klon_global*klev_global))

       if(mype==0) print * ,' recalculate the diff coefficients'  
       zxu = dt*(ra*rdlam)**real(ndifu,realkind)  
       zxv = dt*(ra*rdlam)**real(ndifv,realkind)  
       zxt = dt*(ra*rdlam)**real(ndift,realkind)  
       zxq = dt*(ra*rdlam)**real(ndifq ,realkind) 
       zxs = dt*(ra*rdlam)**real(ndifs ,realkind) 
       zxx = dt*(ra*rdlam)**real(ndifx(1),realkind)  

       zyu = dt*(ra*rdth)**real(ndifu ,realkind) 
       zyv = dt*(ra*rdth)**real(ndifv ,realkind) 
       zyt = dt*(ra*rdth)**real(ndift ,realkind) 
       zyq = dt*(ra*rdth)**real(ndifq ,realkind) 
       zys = dt*(ra*rdth)**real(ndifs  ,realkind)
       zyx = dt*(ra*rdth)**real(ndifx(1),realkind)  


       !     diffusion of u.
       if(ldifu) then  
          ipara = 1  
          do jk = 1,klev  
             do jy = 1,klat  
                do jx = 1,klon  
                   cx(jx,jy,jk)=zxu*cdifu(jk)*(1._realkind/hxv(jx,jy))**real(ndifu,realkind)
                   cy(jx,jy,jk)=zyu*cdifu(jk)*(rhyv(jx,jy))**real(ndifu,realkind)
                enddo
             enddo
          enddo
          !     coefficients in y-direction
          call twod_to_fft(cy,klon,klat,klev,div_fft)  

          call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
          jrec_tri = 0  
          jstart = 1  
          do jslab = 1,nslab_tri  
             jstart = jstart + jrec_tri  
             jrec_tri = 0  
             do i = imin_tri(jslab),imax_tri(jslab)  
                jrec_tri = jrec_tri + 1  
             enddo
             do j = 1,jrec_tri  
                do i = 1,klat_global  
                   diffcy(i,jstart+j-1,ipara)=div_tri(i,jstart+j-1)  
                enddo
             enddo
          enddo
          !     coefficients in x-direction

          call twod_to_fft(cx,klon,klat,klev,div_fft)   
          jrec_fft = 0  

          jstart = 1  
          do jslab = 1,nslab_fft  
             jstart = jstart + jrec_fft  
             jrec_fft = 0  
             do jrec = 1,jlen_fft(jslab)  
                jrec_fft = jrec_fft + 1  

             enddo
             do j = 1,jlen_fft(jslab)  
                do i = 1,klon_global  
                   diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                enddo
             enddo
          enddo
       endif
       !     diffusion of v.
       if(ldifv) then  
          ipara = 2  
          do jk = 1,klev  
             do jy = 1,klat  
                do jx = 1,klon  
                   cx(jx,jy,jk) = zxv * cdifv(jk) *(rhxu(jx,jy) )**real(ndifv,realkind)
                   cy(jx,jy,jk) = zyv * cdifv(jk) *(1._realkind / hyu(jx,jy) )**real(ndifv,realkind)
                enddo
             enddo


          enddo
          !     coefficients in y-direction
          call twod_to_fft(cy,klon,klat,klev,div_fft)   

          call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
          jrec_tri = 0  
          jstart = 1  
          do jslab = 1,nslab_tri  
             jstart = jstart + jrec_tri  
             jrec_tri = 0  
             do i = imin_tri(jslab),imax_tri(jslab)  
                jrec_tri = jrec_tri + 1  
             enddo
             do j = 1,jrec_tri  
                do i = 1,klat_global  
                   diffcy(i,jstart + j - 1,ipara) = div_tri(i,jstart + j - 1)  
                enddo
             enddo


          enddo
          !     coefficients in x-direction

          call twod_to_fft(cx,klon,klat,klev,div_fft)   
          jrec_fft = 0  

          jstart = 1  
          do jslab = 1,nslab_fft  
             jstart = jstart + jrec_fft  
             jrec_fft = 0  
             do jrec = 1,jlen_fft(jslab)  
                jrec_fft = jrec_fft + 1  

             enddo
             do j = 1,jlen_fft(jslab)  
                do i = 1,klon_global  
                   diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                enddo
             enddo

          enddo

       endif
       !     diffusion of temperature.
       if(ldift) then  
          ipara = 3  
          do jk = 1,klev  
             do jx = 1,klon  
                do jy = 1,klat  
                   cx(jx,jy,jk) = zxt * cdift(jk) *(rhxu(jx,jy) )**real(ndift,realkind)
                   cy(jx,jy,jk) = zyt * cdift(jk) *(rhyv(jx,jy) )**real(ndift,realkind)
                enddo
             enddo
          enddo
          !     coefficients in y-direction
          call twod_to_fft(cy,klon,klat,klev,div_fft)   

          call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
          jrec_tri = 0  
          jstart = 1  
          do jslab = 1,nslab_tri  
             jstart = jstart + jrec_tri  
             jrec_tri = 0  
             do i = imin_tri(jslab),imax_tri(jslab)  
                jrec_tri = jrec_tri + 1  

             enddo
             do j = 1,jrec_tri  
                do i = 1,klat_global  
                   diffcy(i,jstart + j - 1,ipara) = div_tri(i,jstart + j - 1)  
                enddo
             enddo
          enddo
          !     coefficients in x-direction

          call twod_to_fft(cx,klon,klat,klev,div_fft)   
          jrec_fft = 0  

          jstart = 1  
          do jslab = 1,nslab_fft  
             jstart = jstart + jrec_fft  
             jrec_fft = 0  
             do jrec = 1,jlen_fft(jslab)  
                jrec_fft = jrec_fft + 1  

             enddo
             do j = 1,jlen_fft(jslab)  
                do i = 1,klon_global  
                   diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                enddo
             enddo
          enddo
       endif
       !     diffusion of humidity.

       if(ldifq) then  
          ipara = 4  
          do jk = 1,klev  
             do jy = 1,klat  
                do jx = 1,klon  
                   cx(jx,jy,jk) = zxq * cdifq(jk) *(rhxu(jx,jy) )**real(ndifq,realkind)
                   cy(jx,jy,jk) = zyq * cdifq(jk) *(rhyv(jx,jy) )**real(ndifq,realkind)
                enddo
             enddo
          enddo
          !     coefficients in y-direction
          call twod_to_fft(cy,klon,klat,klev,div_fft)   

          call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
          jrec_tri = 0  
          jstart = 1  
          do jslab = 1,nslab_tri  
             jstart = jstart + jrec_tri  
             jrec_tri = 0  
             do i = imin_tri(jslab),imax_tri(jslab)  
                jrec_tri = jrec_tri + 1  

             enddo
             do j = 1,jrec_tri  
                do i = 1,klat_global  
                   diffcy(i,jstart + j - 1,ipara) = div_tri(i,jstart + j - 1)  
                enddo
             enddo
          enddo
          !
          !     coefficients in x-direction
          !
          call twod_to_fft(cx,klon,klat,klev,div_fft)   
          !
          jrec_fft = 0  

          jstart = 1  
          do jslab = 1,nslab_fft  
             jstart = jstart + jrec_fft  
             jrec_fft = 0  
             do jrec = 1,jlen_fft(jslab)  
                jrec_fft = jrec_fft + 1  

             enddo
             do j = 1,jlen_fft(jslab)  
                do i = 1,klon_global  
                   diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                enddo
             enddo
          enddo
          !
       endif
       !
       !----------------------------------------------------------------------
       !
       !     diffusion of liquid water
       !     -------------------------
       !
       if(ldifs) then  
          !
          ipara = 5  
          do jk = 1,klev  
             do jy = 1,klat  
                do jx = 1,klon  
                   !     mac
                   cx(jx,jy,jk) = zxs * cdifs(jk) *(rhxu(jx,jy) )**real(ndifs,realkind)
                   cy(jx,jy,jk) = zys * cdifs(jk) *(rhyv(jx,jy) )**real(ndifs,realkind)
                   !     mac
                enddo
             enddo
          enddo
          !
          !
          !     coefficients in y-direction
          !
          call twod_to_fft(cy,klon,klat,klev,div_fft)   
          call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
          !
          jrec_tri = 0  
          jstart = 1  
          do jslab = 1,nslab_tri  
             jstart = jstart + jrec_tri  
             jrec_tri = 0  
             do i = imin_tri(jslab),imax_tri(jslab)  
                jrec_tri = jrec_tri + 1  
             enddo
             do j = 1,jrec_tri  
                do i = 1,klat_global  
                   diffcy(i,jstart + j - 1,ipara) = div_tri(i,jstart + j - 1)  
                enddo
             enddo
          enddo
          !
          !     coefficients in x-direction
          !
          call twod_to_fft(cx,klon,klat,klev,div_fft)   
          !
          jrec_fft = 0  

          jstart = 1  
          do jslab = 1,nslab_fft  
             jstart = jstart + jrec_fft  
             jrec_fft = 0  
             do jrec = 1,jlen_fft(jslab)  
                jrec_fft = jrec_fft + 1  

             enddo
             do j = 1,jlen_fft(jslab)  
                do i = 1,klon_global  
                   diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                enddo
             enddo
          enddo
          !
       endif
       !
       !----------------------------------------------------------------------
       !
       !     diffusion of extra scalars
       !     -------------------------
       !

       do jl = 1,ksvar  
          if(ldifx(jl) ) then  
             ipara = 5 + jl  
             do jk = 1,klev  
                do jy = 1,klat  
                   do jx = 1,klon  
                      cx(jx,jy,jk) = zxx * cdifx(jk,jl) *(rhxu(jx,jy) )**real(ndifx(jl),realkind)
                      cy(jx,jy,jk) = zyx * cdifx(jk,jl) *(rhyv(jx,jy) )**real(ndifx(jl),realkind)
                   enddo
                enddo

             enddo
             !     coefficients in y-direction
             call twod_to_fft(cy,klon,klat,klev,div_fft)   

             call fft_to_tri(div_fft,kpbpts,klev,div_tri)  
             jrec_tri = 0  
             jstart = 1  
             do jslab = 1,nslab_tri  
                jstart = jstart + jrec_tri  
                jrec_tri = 0  
                do i = imin_tri(jslab),imax_tri(jslab)  
                   jrec_tri = jrec_tri + 1  

                enddo
                do j = 1,jrec_tri  
                   do i = 1,klat_global  
                      diffcy(i,jstart + j - 1,ipara) = div_tri(i,jstart + j - 1)  
                   enddo
                enddo
             enddo
             !     coefficients in x-direction
             call twod_to_fft(cx,klon,klat,klev,div_fft)   
             !
             jrec_fft = 0  

             jstart = 1  
             do jslab = 1,nslab_fft  
                jstart = jstart + jrec_fft  
                jrec_fft = 0  
                do jrec = 1,jlen_fft(jslab)  
                   jrec_fft = jrec_fft + 1  

                enddo
                do j = 1,jlen_fft(jslab)  
                   do i = 1,klon_global  
                      diffcx(i,jstart + j - 1,ipara) = div_fft(i,jstart + j - 1)  
                   enddo
                enddo
             enddo
          endif
       enddo
       rimpdt = dt  
       deallocate(div_fft,div_tri)
    endif
    deallocate(cx,cy )
    return  
  end subroutine difhini


  subroutine diffh(klon  ,klat  ,klev,               &
        u     ,v     ,t     ,q     ,alnps ,s,     &
        svar  ,ksvar,                                 &
        dt,   & 
        atcref,                                        &
        kpbpts,ahyb,bhyb)

    !     a mcdonald.
    !     modified by j.e. haugen

    !     this is a program to perform implicit horizontal diffusion on
    !     temperature(t),humidity(q),winds(u,v),liquid water(s) and
    !     extra scalars(svar(ksvar))

    !     i am solving the following equations:
    !     
    !     [1 - dt*ku*gradm]u(n+1) = u(n)
    !     
    !     [1 - dt*kv*gradm]v(n+1) = v(n)
    !     
    !     [1 - dt*kt*gradm]{t(n+1) - tc(n)} = t(n) -tc(n)

    !     where tc is defined in(2.1.4.7) of the hirlam black book as
    !     tc = ln(ps)*[ps*(dp/dps)*(dt/dp)]ref

    !     [1 - dt*kq*gradm]q(n+1) = q(n)

    !     [1 - dt*ks*gradm]s(n+1) = s(n)

    !     [1 - dt*ksvar*gradm]svar(n+1) = svar(n)

    !     kx = diffusion coefficient
    !     dt = time step

    !     gradm =(d/dx)**2 +(d/dy)**2 for ndifx = 2
    !     gradm =(d/dx)**4 +(d/dy)**4 for ndifx = 4
    !     gradm =(d/dx)**6 +(d/dy)**6 for ndifx = 6

    !     klon      : number of grid points in the 'x' direction.
    !     klat      : number of grid points in the 'y' direction.
    !     klev      : number of grid points in the 'z' direction.
    !     ldifu     : is .true. if diffusing the u-component of wind.
    !     ldifv     : is .true. if diffusing the v-component of wind.
    !     ldift     : is .true. if diffusing the temperature.
    !     ldifq     : is .true. if diffusing the humidity.
    !     ldifs     : is .true. if diffusing the liquid water.
    !     ldifx     : is .true. if diffusing the extra scalars.
    !     ndifu     : order of diffusion for u-component of wind.
    !     ndifv     : order of diffusion for v-component of wind.
    !     ndift     : order of diffusion for temperature.
    !     ndifq     : order of diffusion for humidity.
    !     ndifs     : order of diffusion for liquid water.
    !     ndifx     : order of diffusion for extra scalars.
    !     cdifx     : level dependent diffusion coefficient
    !     ksvar     : number of extra scalars.

    !     u,v,t,q,alnps,s,svar:
    !     u-velocity,v-velocity,tepmerature,specific
    !     humidity,log(surface pressure),liquid water
    !     and extra scalars.

    !     atcref    : is generated by subroutine difini which must be called
    !     before diffh.

    implicit none

    integer,intent(in):: klon,klat,klev,ksvar,kpbpts             
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind)    dt    

    real(kind=realkind)   u(klon,klat,klev),v(klon,klat,klev),      &
         t(klon,klat,klev),q(klon,klat,klev),         &
         alnps(klon,klat),                            &
         s(klon,klat,klev),svar(klon,klat,klev,ksvar),&
         atcref(klev)

    !     declaration of local workspace
    integer jx,jy,l,jk,ipara
    real(kind=realkind)  w1(klon,klat,klev)
    real(kind=realkind)  rbdeps

    rbdeps=1._realkind/((ra*rdth)**2.0_realkind)

    if(abs(dt - rimpdt)>1.e-14_realkind) then

       call difhini(0,klon,   klat,  klev,    &
            dt    ,&
            kpbpts,ksvar,ahyb,bhyb)

    endif

    if(ldifu) then
       !     diffusion of u.

       ipara = 1
       call impsub(klon,  klat,  klev,ksvar,  &
            ndifu, cdifu, ipara,          &
            u,                              &
            rbdeps,kpbpts)


    endif

    if(ldifv) then
       !     diffusion of v.

       ipara = 2
       call impsub(klon,  klat,  klev,ksvar,    &
            ndifv, cdifv, ipara,            &
            v,                                &
            rbdeps,kpbpts)

    endif
    if(ldift) then
       !     diffusion of temperature.

       !     add modification to take mountains into account.

       ipara = 3
       do jk=1,klev
          do jx=1,klon
             do jy=1,klat
                w1(jx,jy,jk) = t(jx,jy,jk) - atcref(jk)*alnps(jx,jy)
             enddo
          enddo
       enddo


       call impsub(klon,  klat,  klev,ksvar,   &
            ndift,cdift, ipara,            &
            w1,                              &
            rbdeps,kpbpts)

       do jk=1,klev
          do jy=1,klat
             do jx=1,klon
                t(jx,jy,jk) = w1(jx,jy,jk) + atcref(jk)*alnps(jx,jy)
             enddo
          enddo
       enddo

    endif

    if(ldifq) then
       !     diffusion of humidity.

       ipara = 4
       call impsub(klon,  klat,  klev,ksvar,    &
            ndifq, cdifq, ipara,           &
            q,                               &
            rbdeps,kpbpts)


    endif
    if(ldifs) then
       !     diffusion of liquid water.

       ipara = 5
       call impsub(klon,  klat,  klev,ksvar,    &
            ndifs, cdifs, ipara,           &
            s,                               &
            rbdeps,kpbpts)

    endif
    !     diffusion of extra scalars.

    do l=1,ksvar
       if(ldifx(l)) then
          ipara = l+5
          call impsub(klon,  klat,  klev,ksvar,    &
               ndifx(l), cdifx(1,l), ipara,   &
               svar(1,1,1,l),                   &
               rbdeps,kpbpts)

       endif
    enddo
    return
  end subroutine diffh



  subroutine impsub(klon,  klat,  klev, ksvar ,     &
       ndifx, cdifx, kpara,         &
       x,&
       rbdeps,kpbpts)
    !       implicit filtering of a two dimensional field using
    !       a second,fourth,or sixth order implicit sin filter
    !       see raymond,m.w.r.,116,2132-2141,and my notes.
    !
    !
    !  this is a subroutine to perform implicit horizontal diffusion on
    !  a given field 'x'.
    !
    !  klon      : number of grid points in the 'x' direction.
    !  klat      : number of grid points in the 'y' direction.
    !  klev      : number of grid points in the 'z' direction.
    !  ndifx     : order of diffusion for variable x.
    !  cdifx     : level dependent diffusion coefficient
    !  x         : variable to be diffused
    !  cy,cx     : coefficients in x- and y-direction
    !  rbdeps    : a parameter used to get the diffusion constant
    !              correct at the first line in for invlo4 and at the first
    !              two lines in for invlo6.
    !  kpbpts    : number of passive points
    !      kalle eerola   hirlam    1998

    use transpose
    implicit none

    integer,intent(in):: klon, klat, klev, ndifx,kpbpts,kpara, ksvar
    real(kind=realkind)    rbdeps

    real(kind=realkind),intent(in)::  cdifx(klev)
    real(kind=realkind)::x(klon,klat,klev)

    integer  i ,j,k,jstart,jslab,jrec_fft,jrec,jrec_tri,jj
    integer  iactive(klon_global),jactive(klat_global)
    real(kind=realkind),allocatable,dimension(:,:):: div_fft,div_tri 
    logical iskip
    integer::status
    data iskip /.false./

    !     1. transform variable 'x' to x-k-plane
    allocate(div_fft(klon_global,klat_global*klev_global),stat=status)
    if(status/=0)then
       stop 'memory allocation error div_fft'
    endif
    allocate(div_tri(klat_global,klon_global*klev_global),stat=status)
    if(status/=0)then
       stop 'memory allocation error div_tri'
    endif
    call twod_to_fft(x,klon,klat,klev,div_fft) 

    !      filter in x-direction
    jrec_fft = 0
    jstart   = 1

    do jslab = 1,nslab_fft
       jstart   = jstart + jrec_fft
       jrec_fft = 0
       j        = jmin_fft(jslab)
       do jrec=1,jlen_fft(jslab)
          if(j ==kpbpts+1.or.j == klat_global-kpbpts-1.or.&
               j == klat_global-kpbpts) then
             jactive(jrec) = 0
          else
             jactive(jrec) = 1
          endif
          jrec_fft = jrec_fft+1
          j        = j + 1
       enddo

       call lowpass(div_fft(1,jstart)  ,               &
            klon_global-1      ,jrec_fft,              &
            klon_global    ,klat_global*klev_global+1, &
            diffcx(1,jstart,kpara),                        &
            klon_global             ,ndiffx,                &
            ndifx,iskip,rbdeps ,jactive)

    enddo

    !     3. transformation to y-k-plane

    call fft_to_tri(div_fft,KPBPTS,KLEV,div_tri)  

    !     4. filter in y-direction
    jrec_tri = 0
    jstart   = 1

    do jslab=1,nslab_tri

       jstart   = jstart + jrec_tri
       jrec_tri = 0
       do i = imin_tri(jslab),imax_tri(jslab)
          jrec_tri = jrec_tri + 1
          if(i==kpbpts+1.or.i==klon_global-kpbpts-1.or. &
               i == klon_global-kpbpts   ) then
             iactive(jrec_tri)=0
          else
             iactive(jrec_tri)=1
          endif
       enddo

       call lowpass(div_tri(1,jstart)   ,                      &
            klat_global-1       ,jrec_tri,             &
            klat_global     ,klon_global*klev_global+1,&
            diffcy(1,jstart,kpara) ,                       &
            klat_global              ,ndiffy,               &
            ndifx,iskip,rbdeps  ,iactive)

    enddo

    !     6. back to 2-d
    call tri_to_fft(div_tri,kpbpts,klev,div_fft)
    call fft_to_twod(div_fft,klon,klat,klev,x)

    call swap(x,klon,klat,klev)
    deallocate(div_fft,div_tri)
    return

  end subroutine impsub



  SUBROUTINE LOWPASS(XY,N,M,K1,K2,EPS,KE1,KE2,NORDER,&
       ISKIP,RBDEPS,KACTIVE)
    !       SECOND,FOURTH,SIXTH-ORDER LOW-PASS IMPLICIT SINE FILTER
    !     (RAYMOND,MWR,116,2132-2141)
    !
    !       XY     UNFILTERED VALUES ON INPUT
    !              FILTERED VALUES ON OUTPUT.
    !       K1,K2  DIMENSION OF THE ARRAY 'XY'
    !       N      NUMBER OF VALUES.
    !       M      NUMBER OF VALUES.
    !       EPS    FILTER PARAMETER
    !            (DETERMINES CUTOFF)
    !       NORDER = 2 FOR SECOND ORDER DIFFUSION.
    !       NORDER = 4 FOR FOURTH ORDER DIFFUSION.
    !       NORDER = 6 FOR SIXTH  ORDER DIFFUSION.
    !       ISKIP  =.TRUE.  : SKIP INITIALIZATION
    !              =.FALSE. : DO INITIALIZATION
    !       RBDEPS 1/(dy**2). A parameter use d to get the diffusion constant
    !              correct at the first line in for invlo4 and at the first
    !              two lines in for invlo6.


    IMPLICIT NONE  
    INTEGER :: N,M,K1,K2,KE1,KE2,NORDER,KFIRST,KLAST  
    INTEGER :: KACTIVE(K2)  
    REAL(KIND=REALKIND) :: XY(K1,K2),EPS(KE1,KE2),RBDEPS  
    LOGICAL :: ISKIP  
    !
    !
    INTEGER :: I,J,NM1,NM2,NM3,NM4,IND  
    REAL(KIND=REALKIND) :: BB(M,N),XANS(M,N),ZEPS(M,N),ZXY(M,N)  

    NM1 = N - 1  
    NM2 = N - 2  
    NM3 = N - 3  
    NM4 = N - 4  
    !     CHANGE THE ORDER OF ROWS AND COLUMNS
    DO I = 1,N  
       DO J = 1,M  
          ZEPS(J,I) = EPS(I,J)  
          ZXY(J,I) = XY(I,J)  
       ENDDO
    ENDDO
    IF(NORDER==2) THEN  
       DO J = 1,M  
          BB(J,1) = 0._realkind  
          BB(J,N) = 0._realkind  
       ENDDO
       DO I = 2,NM1  
          DO J = 1,M  
             BB(J,I) = ZEPS(J,I) *(ZXY(J,I - 1) - 2._realkind * ZXY(J,I) &
                  + ZXY(J,I + 1) )
          ENDDO
       ENDDO
       !       SOLVE THE LINEAR SYSTEM FOR XANS
       CALL INVLO2(BB,XANS,M,N,ZEPS,ISKIP)  
    ELSEIF(NORDER==4) THEN  
       DO J = 1,M  
          BB(J,1) = 0._realkind  
          BB(J,N) = 0._realkind  
          BB(J,2) = ZEPS(J,2) * RBDEPS *(ZXY(J,1) - 2._realkind * ZXY(J,&
               2) + ZXY(J,3) )
          BB(J,NM1) = ZEPS(J,NM1) * RBDEPS *(ZXY(J,NM2) - 2._realkind * &
               ZXY(J,NM1) + ZXY(J,N) )
       ENDDO
       DO I = 3,NM2  
          DO J = 1,M  
             BB(J,I) = ZEPS(J,I) *( - 1._realkind *(ZXY(J,I - 2) + ZXY(J,I &
                  + 2) ) + 4._realkind *(ZXY(J,I - 1) + ZXY(J,I + 1) ) - 6._realkind * ZXY(J,&
                  I) )
          ENDDO
       ENDDO
       !       SOLVE THE LINEAR SYSTEM FOR XANS
       CALL INVLO4(BB,XANS,M,N,ZEPS,ISKIP,RBDEPS)  
    ELSEIF(NORDER==6) THEN  
       DO J = 1,M  
          BB(J,1) = 0._realkind  
          BB(J,N) = 0._realkind  
          BB(J,2) = ZEPS(J,2) * RBDEPS * RBDEPS *(ZXY(J,1) &
               - 2._realkind * ZXY(J,2) + ZXY(J,3) )
          BB(J,NM1) = ZEPS(J,NM1) * RBDEPS * RBDEPS *(ZXY(J,NM2) &
               - 2._realkind * ZXY(J,NM1) + ZXY(J,N) )
          BB(J,3) = ZEPS(J,3) * RBDEPS *( - 1._realkind *(ZXY(J,1) &
               + ZXY(J,5) ) + 4._realkind *(ZXY(J,2) + ZXY(J,4) ) - 6._realkind * ZXY(J,&
               3) )
          BB(J,NM2) = ZEPS(J,NM2) * RBDEPS *( - 1._realkind *(ZXY(J,N) &
               + ZXY(J,NM4) ) + 4._realkind *(ZXY(J,NM1) + ZXY(J,NM3) ) - 6._realkind * &
               ZXY(J,NM2) )
       ENDDO
       DO I = 4,NM3  
          DO J = 1,M  
             BB(J,I) = ZEPS(J,I) *((ZXY(J,I - 3) + ZXY(J,I + 3) ) &
                  - 6._realkind *(ZXY(J,I - 2) + ZXY(J,I + 2) ) + 15._realkind *(ZXY(J,I - &
                  1) + ZXY(J,I + 1) ) - 20._realkind * ZXY(J,I) )
          ENDDO
       ENDDO
       !       SOLVE THE LINEAR SYSTEM FOR XANS
       CALL INVLO6(BB,XANS,M,N,ZEPS,ISKIP,RBDEPS)  
    ELSE  

       WRITE(6,* ) 'IN *LOWPASS* - NORDER MUST BE 2,4 OR 6 not',norder  
       CALL STOP_PROGRAM('')  
    ENDIF
    !     ADD CORRECTION TO GET FILTERED VALUES.
    IND = 0  
    DO J = 1,M  
       IND = J  
       IF(KACTIVE(J) ==1) THEN  
          DO I = 1,N  
             XY(I,J) = XY(I,J) + XANS(IND,I)  
          ENDDO
       ENDIF
    ENDDO
    RETURN  
  END SUBROUTINE LOWPASS



SUBROUTINE INVLO2(BB,XANS,M,N,EPS,ISKIP)  


  implicit none  
  !     GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
  !     SECOND-ORDER LOW-PASS IMPLICIT SINE FILTER.
  !   (REF: WILLIAM H RAYMOND,MWR,116,2132-2124)
  INTEGER :: M,N  
  REAL(KIND=REALKIND) :: BB(M,N),XANS(M,N),EPS(M,N)  


  LOGICAL :: ISKIP  
  INTEGER :: I,J  


  REAL(KIND=REALKIND) :: B(M,N),C(M,N),D(M,N),DELTA(M,N),BETA(M,N),&
       W(M,N),H(M,N)
  !     SKIP INITIALIZATION OF MATRIX ON REPEAT CALLS.

  IF(ISKIP) GOTO 100  
  !     INITIALIZE THE MATRIX
  DO 10 I = 2,N - 1  
     DO J = 1,M  
        B(J,I) = - EPS(J,I)  
        C(J,I) = 1._realkind + 2._realkind * EPS(J,I)  
        D(J,I) = B(J,I)  
     ENDDO

10 END DO
  DO J = 1,M  
     B(J,1) = 0._realkind  
     C(J,1) = 1._realkind  
     D(J,1) = 0._realkind  
     B(J,N) = 0._realkind  
     C(J,N) = 1._realkind  
     D(J,N) = 0._realkind  

  ENDDO
  !     STEP ONE.
  DO J = 1,M  
     BETA(J,1) = D(J,1) / C(J,1)  
     DELTA(J,2) = B(J,2)  
     W(J,1) = C(J,1)  
     W(J,2) = C(J,2) - DELTA(J,2) * BETA(J,1)  
     BETA(J,2) = D(J,2) / W(J,2)  
     DELTA(J,3) = B(J,3)  
     W(J,3) = C(J,3) - DELTA(J,3) * BETA(J,2)  
     BETA(J,3) = D(J,3) / W(J,3)  

  ENDDO
  !     STEP TWO
  DO 20 I = 4,N  
     DO J = 1,M  
        DELTA(J,I) = B(J,I)  
        W(J,I) = C(J,I) - DELTA(J,I) * BETA(J,I - 1)  
        BETA(J,I) = D(J,I) / W(J,I)  
     ENDDO

20 END DO

100 CONTINUE  
  !     STEP THREE
  DO J = 1,M  
     H(J,1) = BB(J,1) / W(J,1)  
     H(J,2) =(BB(J,2) - DELTA(J,2) * H(J,1) ) / W(J,2)  
     H(J,3) =(BB(J,3) - DELTA(J,3) * H(J,2) ) / W(J,3)  
  ENDDO
  DO 30 I = 4,N  
     DO J = 1,M  
        H(J,I) =(BB(J,I) - DELTA(J,I) * H(J,I - 1) ) / W(J,I)
     ENDDO

30 END DO
  !     STEP FOUR
  DO J = 1,M  
     XANS(J,N) = H(J,N)  
     XANS(J,N - 1) = H(J,N - 1) - BETA(J,N - 1) * XANS(J,N)  
     XANS(J,N - 2) = H(J,N - 2) - BETA(J,N - 2) * XANS(J,N - 1)  
  ENDDO
  DO 40 I = N - 3,1,- 1  
     DO J = 1,M  
        XANS(J,I) = H(J,I) - BETA(J,I) * XANS(J,I + 1)  
     ENDDO

40 END DO
  RETURN  
END subroutine INVLO2


SUBROUTINE INVLO4(BB,XANS,M,N,EPS,ISKIP,RBDEPS)  
  !
  !       GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
  !
  !       FOURTH-ORDER LOW-PASS IMPLICIT SINE FILTER.
  !     (REF: WILLIAM H RAYMOND,MWR,116,2132-2124)
  !
  IMPLICIT NONE  
  !
  INTEGER :: M,N  
  REAL(KIND=REALKIND) :: BB(M,N),XANS(M,N),EPS(M,N)  
  REAL(KIND=REALKIND) :: RBDEPS  
  LOGICAL :: ISKIP  
  !
  !
  INTEGER :: I,J  
  REAL(KIND=REALKIND) :: A(M,N),B(M,N),C(M,N),D(M,N),E(M,N),DELTA( &
       M,N),BETA(M,N),W(M,N),GAM(M,N),H(M,N),AP(M,N)
  !       SKIP INITIALIZATION OF MATRIX ON REPEAT CALLS.
  IF(ISKIP) GOTO 100  
  !       INITIALIZE THE MATRIX
  DO 10 I = 3,N - 2  
     DO J = 1,M  
        A(J,I) = EPS(J,I)  
        B(J,I) = - 4._realkind * EPS(J,I)  
        C(J,I) = 1._realkind + 6._realkind * EPS(J,I)  
        D(J,I) = B(J,I)  
        E(J,I) = A(J,I)  
     ENDDO
10 END DO
  !
  DO J = 1,M  
     A(J,1) = 0._realkind  
     A(J,2) = 0._realkind  
     B(J,1) = 0._realkind  
     B(J,2) = - EPS(J,2) * RBDEPS  
     C(J,1) = 1._realkind  
     C(J,2) = 1._realkind + 2._realkind * EPS(J,2) * RBDEPS  
     D(J,1) = 0._realkind  
     D(J,2) = - EPS(J,2) * RBDEPS  
     E(J,1) = 0.0_realkind  
     E(J,2) = 0.0_realkind  
     A(J,N - 1) = 0.0_realkind  
     A(J,N) = 0.0_realkind  
     B(J,N - 1) = - EPS(J,N - 1) * RBDEPS  
     B(J,N) = 0._realkind  
     C(J,N - 1) = 1._realkind + 2._realkind * EPS(J,N - 1) * RBDEPS  
     C(J,N) = 1._realkind  
     D(J,N - 1) = - EPS(J,N - 1) * RBDEPS  
     D(J,N) = 0.0_realkind  
     E(J,N - 1) = 0.0_realkind  
     E(J,N) = 0.0_realkind  
  ENDDO
  !
  !       STEP ONE.
  DO J = 1,M  
     BETA(J,1) = D(J,1) / C(J,1)  
     DELTA(J,2) = B(J,2)  
     W(J,1) = C(J,1)  
     AP(J,1) = 0._realkind  
     AP(J,2) = 0._realkind  
     AP(J,3) = A(J,3)  
     W(J,2) = C(J,2) - DELTA(J,2) * BETA(J,1)  
     GAM(J,1) = E(J,1) / C(J,1)  
     BETA(J,2) =(D(J,2) - DELTA(J,2) * GAM(J,1) ) / W(J,2)  
     GAM(J,2) =(E(J,2) ) / W(J,2)  
     DELTA(J,3) =(B(J,3) - AP(J,3) * BETA(J,1) )  
     W(J,3) = C(J,3) - DELTA(J,3) * BETA(J,2) - AP(J,3) * GAM(J,1)
     BETA(J,3) =(D(J,3) - DELTA(J,3) * GAM(J,2) ) / W(J,3)  
     GAM(J,3) =(E(J,3) ) / W(J,3)  
  ENDDO
  !

  !       STEP TWO
  DO 20 I = 4,N  
     DO J = 1,M  
        AP(J,I) = A(J,I)  
        DELTA(J,I) = B(J,I) - AP(J,I) * BETA(J,I - 2)  
        W(J,I) = C(J,I) - AP(J,I) * GAM(J,I - 2) - DELTA(J,I) &
             * BETA(J,I - 1)
        BETA(J,I) =(D(J,I) - DELTA(J,I) * GAM(J,I - 1) ) &
             / W(J,I)
        GAM(J,I) =(E(J,I) ) / W(J,I)  
     ENDDO
20 END DO
  !
100 CONTINUE  
  !

  !       STEP THREE
  DO J = 1,M  
     H(J,1) = BB(J,1) / W(J,1)  
     H(J,2) =(BB(J,2) - DELTA(J,2) * H(J,1) ) / W(J,2)  
     H(J,3) =(BB(J,3) - DELTA(J,3) * H(J,2) - AP(J,3) &
          * H(J,1) ) / W(J,3)
  ENDDO
  DO 30 I = 4,N  
     DO J = 1,M  
        H(J,I) =(BB(J,I) - DELTA(J,I) * H(J,I - 1) - AP(J,I) &
             * H(J,I - 2) ) / W(J,I)
     ENDDO
30 END DO
  !       STEP FOUR
  DO J = 1,M  
     XANS(J,N) = H(J,N)  
     XANS(J,N - 1) = H(J,N - 1) - BETA(J,N - 1) * XANS(J,N)  
     XANS(J,N - 2) = H(J,N - 2) - BETA(J,N - 2) * XANS(J,N - 1) &
          - GAM(J,N - 2) * XANS(J,N)
  ENDDO
  DO 40 I = N - 3,1,- 1  
     DO J = 1,M  
        XANS(J,I) = H(J,I) - BETA(J,I) * XANS(J,I + 1) - GAM( &
             J,I) * XANS(J,I + 2)
     ENDDO
40 END DO
  
  RETURN  
END SUBROUTINE INVLO4


SUBROUTINE INVLO6(BB,XANS,M,N,EPS,ISKIP,RBDEPS)  
  !       GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
  !
  !       SIXTH-ORDER LOW-PASS IMPLICIT SINE FILTER.
  !     (REF: WILLIAM H RAYMOND,MWR,116,2132-2124)
  !

  IMPLICIT NONE  
  INTEGER :: M,N  
  REAL(KIND=REALKIND) :: RBDEPS  
  REAL(KIND=REALKIND) :: BB(M,N),XANS(M,N),EPS(M,N)  
  LOGICAL :: ISKIP  
  !
  !
  INTEGER :: J,I  

  REAL(KIND=REALKIND) :: A(M,N),B(M,N),C(M,N),D(M,N),E(M,N),DELTA( &
       M,N),BETA(M,N),W(M,N),GAM(M,N),H(M,N),PI(M,N),&
       AP(M,N),F(M,N),Z(M,N)
  !       SKIP INITIALIZATION OF MATRIX ON REPEAT CALLS.

  IF(ISKIP) GOTO 100  
  !       INITIALIZE THE MATRIX
  !
  DO 10 I = 4,N - 3  
     DO J = 1,M  
        Z(J,I) = - EPS(J,I)  
        A(J,I) = 6._realkind * EPS(J,I)  
        B(J,I) = - 15._realkind * EPS(J,I)  
        C(J,I) = 1._realkind + 20._realkind * EPS(J,I)  
        D(J,I) = B(J,I)  
        E(J,I) = A(J,I)  
        F(J,I) = Z(J,I)  
     ENDDO
10 END DO
  !
  DO J = 1,M  
     Z(J,1) = 0.0_realkind  
     Z(J,2) = 0.0_realkind  
     Z(J,3) = 0.0_realkind  
     A(J,1) = 0.0_realkind  
     A(J,2) = 0.0_realkind  
     A(J,3) = EPS(J,3) * RBDEPS  
     B(J,1) = 0.0_realkind  
     B(J,2) = - EPS(J,2) * RBDEPS * RBDEPS  
     B(J,3) = - 4._realkind * EPS(J,3) * RBDEPS  
     C(J,1) = 1._realkind  
     C(J,2) = 1._realkind + 2._realkind * EPS(J,2) * RBDEPS * RBDEPS  
     C(J,3) = 1._realkind + 6._realkind * EPS(J,3) * RBDEPS  
     D(J,1) = 0.0_realkind  
     D(J,2) = - EPS(J,2) * RBDEPS * RBDEPS  
     D(J,3) = - 4._realkind * EPS(J,3) * RBDEPS  
     E(J,1) = 0.0_realkind  
     E(J,2) = 0.0_realkind  
     E(J,3) = EPS(J,3) * RBDEPS  
     F(J,1) = 0.0_realkind  
     F(J,2) = 0.0_realkind  
     F(J,3) = 0.0_realkind 
     Z(J,N - 2) = 0.0_realkind 
     Z(J,N - 1) = 0.0_realkind  
     Z(J,N) = 0._realkind 
     A(J,N - 2) = EPS(J,N - 2) * RBDEPS  
     A(J,N - 1) = 0.0_realkind  
     A(J,N) = 0._realkind  
     B(J,N - 2) = - 4._realkind * EPS(J,N - 2) * RBDEPS  
     B(J,N - 1) = - EPS(J,N - 1) * RBDEPS * RBDEPS  
     B(J,N) = 0.0_realkind  
     C(J,N - 2) = 1._realkind + 6._realkind * EPS(J,N - 2) * RBDEPS  
     C(J,N - 1) = 1._realkind + 2._realkind * EPS(J,N - 1) * RBDEPS * RBDEPS  
     C(J,N) = 1._realkind  
     D(J,N - 2) = - 4._realkind * EPS(J,N - 2) * RBDEPS  
     D(J,N - 1) = - EPS(J,N - 1) * RBDEPS * RBDEPS  
     D(J,N) = 0.0_realkind  
     E(J,N - 2) = EPS(J,N - 2) * RBDEPS  
     E(J,N - 1) = 0.0_realkind 
     E(J,N) = 0.0_realkind  
     F(J,N - 2) = 0.0_realkind  
     F(J,N - 1) = 0.0_realkind  
     F(J,N) = 0.0_realkind  
  ENDDO

  !       STEP ONE.
  DO J = 1,M  
     BETA(J,1) = D(J,1) / C(J,1)  
     DELTA(J,2) = B(J,2)  
     W(J,1) = C(J,1)  
     PI(J,1) = F(J,1) / W(J,1)  
     AP(J,1) = 0.0_realkind  
     AP(J,2) = 0.0_realkind 
     AP(J,3) = A(J,3)  
     W(J,2) = C(J,2) - DELTA(J,2) * BETA(J,1)  
     GAM(J,1) = E(J,1) / C(J,1)  
     BETA(J,2) =(D(J,2) - DELTA(J,2) * GAM(J,1) ) / W(J,2)  
     GAM(J,2) =(E(J,2) - PI(J,1) * DELTA(J,2) ) / W(J,2)  
     PI(J,2) = F(J,2) / W(J,2)  
     DELTA(J,3) =(B(J,3) - AP(J,3) * BETA(J,1) )  
     W(J,3) = C(J,3) - DELTA(J,3) * BETA(J,2) - AP(J,3) &
          * GAM(J,1)
     BETA(J,3) =(D(J,3) - AP(J,3) * PI(J,1) - DELTA(J,3) &
          * GAM(J,2) ) / W(J,3)
     GAM(J,3) =(E(J,3) - DELTA(J,3) * PI(J,2) ) / W(J,3)  
     PI(J,3) = F(J,3) / W(J,3)  
  ENDDO

  !       STEP TWO
  DO 20 I = 4,N  
     DO J = 1,M  
        AP(J,I) = A(J,I) - Z(J,I) * BETA(J,I - 3)  
        DELTA(J,I) = B(J,I) - AP(J,I) * BETA(J,I - 2) - Z(J,&
             I) * GAM(J,I - 3)
        W(J,I) = C(J,I) - AP(J,I) * GAM(J,I - 2) - DELTA(J,I) &
             * BETA(J,I - 1) - Z(J,I) * PI(J,I - 3)
        BETA(J,I) =(D(J,I) - AP(J,I) * PI(J,I - 2) - DELTA(J,&
             I) * GAM(J,I - 1) ) / W(J,I)
        GAM(J,I) =(E(J,I) - DELTA(J,I) * PI(J,I - 1) ) &
             / W(J,I)
        PI(J,I) = F(J,I) / W(J,I)  
     ENDDO
20 END DO

100 CONTINUE  
  !       STEP THREE
  DO J = 1,M  
     H(J,1) = BB(J,1) / W(J,1)  
     H(J,2) =(BB(J,2) - DELTA(J,2) * H(J,1) ) / W(J,2)  
     H(J,3) =(BB(J,3) - DELTA(J,3) * H(J,2) - AP(J,3) &
          * H(J,1) ) / W(J,3)
  ENDDO
  DO 30 I = 4,N  
     DO J = 1,M  
        H(J,I) =(BB(J,I) - DELTA(J,I) * H(J,I - 1) - AP(J,I) &
             * H(J,I - 2) - Z(J,I) * H(J,I - 3) ) / W(J,I)
     ENDDO
30 END DO

  !       STEP FOUR
  DO J = 1,M  
     XANS(J,N) = H(J,N)  
     XANS(J,N - 1) = H(J,N - 1) - BETA(J,N - 1) * XANS(J,N)  
     XANS(J,N - 2) = H(J,N - 2) - BETA(J,N - 2) * XANS(J,N - 1) &
          - GAM(J,N - 2) * XANS(J,N)
  ENDDO
  DO 40 I = N - 3,1,- 1  
     DO J = 1,M  
        XANS(J,I) = H(J,I) - BETA(J,I) * XANS(J,I + 1) - GAM( &
             J,I) * XANS(J,I + 2) - PI(J,I) * XANS(J,I + 3)
     ENDDO
40 END DO
  !
  RETURN  
END SUBROUTINE INVLO6



subroutine difini(klev,ksvar,ahyb,bhyb)  
  !     difini - diffusion scheme
  !     j.e. haugen        hirlam
  !     k.s. eerola        hirlam4(revised)
  !     define constants for linear explicit 4.order horizontal diffusion
  !     also the constansts for temperature correction due to the
  !     orography is calculated here.

  use comhkp,only:nltcrf  
  use referenceParameters,only:sip0
  implicit none  

  integer,intent(in) :: klev,ksvar  
  real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
  integer :: jk,j1,j2,j,l  
  real(kind=realkind) :: zalpha,zlnprs,zprerk,ztemrc,ztcref  

  character(len=80) :: mes  
  logical::debug
  debug=.false.

  zlnprs = log(sip0)  
  rprers = sip0  
  rtemrs = 288.0_realkind  
  rtemrt = 216.5_realkind  
  zalpha = 1._realkind/5.256_realkind  


  if(mype==0.and.debug)print *,' ----------- in difini ------------'  
  if(.not.nltcrf) then  
     if(mype==0)print *,' no correction term for temp.-field'  
     goto 110  
  endif


  if(mype==0.and.debug)then
     print *,'rprers',rprers
     print *,'rtemrs',rtemrs
     print *,'rtemrt',rtemrt  
  endif


  if(debug)then
     if(mype==0)print *,' '  
     write(mes, * ) 'level           prk        trc      tcref     tc'
     if(mype==0)print *,mes  
     write(mes, * ) '-----------------------------------------------'
     if(mype==0)print *,mes  
  endif

  do jk = 1,klev  
     zprerk = 0.5_realkind*(ahyb(jk+1)+ahyb(jk))+0.5_realkind*(bhyb(jk+1)+bhyb(jk))*rprers
     ztemrc = rtemrs *(zprerk / rprers) **zalpha  
     atcref(jk) = 0._realkind 

     if(ztemrc>rtemrt)then
        atcref(jk) = 0.5_realkind*(bhyb(jk+1)+bhyb(jk))*zalpha*ztemrc*rprers/zprerk
     endif

     ztcref = atcref(jk) * zlnprs  

     if(debug)then
        write(mes,'(1x,i5,f15.5,3f11.5)')jk,zprerk,ztemrc,atcref(jk),ztcref
        if(mype==0)print *,mes  
     endif
  enddo

  if(debug)then
     write(mes, * ) '-----------------------------------------------'
     if(mype==0)print *,mes  
  endif

110 continue  

  if(nlhdif) then  
     if(mype==0.and.debug)then
        print *,' '  
!        print *,'twodt ',twodt
!        print *,'dlamda',dlamda  
        print *,'ak4   ',ak4  
        print *,' '  
        print *,' weight factors multiplied with ak4: '  
        print *,' ----------------------------------- '  
        print *,' '  
        print *,'------------------------------------------- -'
        print *,'level      t       u       v       q        &s'
        print *,mes  
        do jk = 1,klev  
           write(mes,'(1x,i5,5(1x,f8.4))')jk,(ak4lev((j-1)*klev+jk),j=1,5)
           print *,mes  
        enddo
        if(mype==0)print *,' '  

        !        extra scalars
        if(ksvar>0) then  
           do l = 1,ksvar,5  
              j1 = l  
              j2 = min(l + 4,ksvar)  
              print *,' '  
              write(mes,'(1x,''level'',5('' svar('',i3,'')''))')(j,j=j1,j2)
              print *,mes  
              print *,'-------------------------------------------'
              write(mes,'(1x,''level'',5('' svar('',i3,'')''))')(j,j=j1,j2)
              print *,mes  
              print *,'-------------------------------------------'
              do jk = 1,klev  
                 write(mes,'(1x,i5,5(1x,f9.4))') jk,(ak4lev((4+j)*klev+jk),j=j1,j2)
                 print *,mes  
              enddo
              print *,'-------------------------------------------'
           enddo
        endif

     endif
     print *,' '  
     print *,'---------- difini done -----------'  
  endif
  if(mype==0) write(6,*) ' ak4 = ',ak4
end subroutine difini

end module mod_diffh
