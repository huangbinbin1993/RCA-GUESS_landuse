module sldynm5
  use decomp
  use comhkp
  use referenceParameters
  use calendar
  use RCAdomainMod
  use mod_implicit_solver
  use boundaryRelaxation
  use confys

  implicit none
  private
!     combix - work-space variables in the interpolation routines
!     parameter statement for tabulated interpolation coefficient:
!     maxint - number of subinterval between each model level
!     (if the spacing between levels is equally distributed)
!     maxtab - length of tables

  type:: slweights
     integer,allocatable,dimension(:,:)::kp
     integer,allocatable,dimension(:,:)::kq
     integer,allocatable,dimension(:,:)::kr
     real(kind=realkind),allocatable,dimension(:,:,:)::alfa
     real(kind=realkind),allocatable,dimension(:,:,:)::beta
     real(kind=realkind),allocatable,dimension(:,:,:)::gama
  end type slweights

  real(kind=realkind),allocatable,dimension(:),save::etafull
  real(kind=realkind),allocatable,dimension(:),save::etahalf
  real(kind=realkind),allocatable,dimension(:,:),save::tablinf
  real(kind=realkind),allocatable,dimension(:,:),save::tabqdrf
  real(kind=realkind),allocatable,dimension(:,:),save::tabcubf
  real(kind=realkind),allocatable,dimension(:),save::tablinh
  real(kind=realkind),allocatable,dimension(:),save::tabqdrh
  real(kind=realkind),allocatable,dimension(:),save::tabcubh
  integer,save:: maxtab
  integer,allocatable,dimension(:,:),save::tabetaf
  integer,allocatable,dimension(:),save::tabetah
  real(kind=realkind),save::etainc,etamin
  integer,save::iwest=-1
  integer,save::isouth=-1
  integer,save::ieast=-1
  integer,save::inorth=-1  
  integer,save::khalo=-666
#ifdef MPI_SRC
  integer,allocatable,dimension(:),target::sliside,sljside
#endif
  integer,save::maxHalo=20

  public sldynm 
contains  

  type(slweights) function allocWeights(klon,klat)
    integer,intent(in)::klon,klat
    allocate(allocWeights%kp(klon,klat))
    allocate(allocWeights%kq(klon,klat))
    allocate(allocWeights%kr(klon,klat))
    allocate(allocWeights%alfa(klon,klat,4))
    allocate(allocWeights%beta(klon,klat,4))
    allocate(allocWeights%gama(klon,klat,4))
  end function allocWeights

  subroutine deallocWeights(w)
    type(slweights)::w
    deallocate(w%kp)
    deallocate(w%kq)
    deallocate(w%kr)
    deallocate(w%alfa)
    deallocate(w%beta)
    deallocate(w%gama)
  end subroutine deallocWeights

  subroutine allocate_combix_memory(klon,klat,klev)

    implicit none
    integer, intent(in)::klon,klat,klev
    integer,parameter::maxint=1000
    integer::halo
#ifdef MPI_SRC
#include"mpif.h"
    integer::ierr
#endif
    maxtab = maxint*klev + 3
    etainc = 1.0_realkind/real(maxtab-3,realkind)
    etamin = -etainc
    allocate(etafull(klev), etahalf(klev+1),&
         tabetaf(maxtab,4), tabetah(maxtab*4), &
         tablinf(maxtab,2), tabqdrf(maxtab,3), tabcubf(maxtab,4),&
         tablinh(maxtab*2), tabqdrh(maxtab*3), tabcubh(maxtab*4))
    
#ifdef MPI_SRC
    halo = min(klon,klat) !the halo cannot be larger than the smallest innerDomain size
    call MPI_allreduce(halo,maxhalo,1,mpi_integer,mpi_min,localcomm,ierr)
    allocate(sliside(maxHalo),sljside(maxHalo))
    !count, blocklen, stride
    do halo=1,maxHalo
       call MPI_TYPE_VECTOR(klev, & !count
            (klon+2*(halo-1))*halo,   & !blocklen
            (klon+2*(halo-1))*(klat+2*(halo-1)),&!stride
            REALTYPE,sljside(halo),ierr)
       call MPI_TYPE_COMMIT(sljside(halo),ierr)
       call MPI_TYPE_VECTOR((klat+2*(halo-1))*klev, &!count
            halo, &!blocklen 
            klon+2*(halo-1),&!stride 
            REALTYPE,sliside(halo),ierr)
       call MPI_TYPE_COMMIT(sliside(halo),ierr)
    enddo
#endif
  end subroutine allocate_combix_memory

  integer function compute_halo(uz,vz,nslind,safety_factor)
    implicit none
    integer,intent(in)::nslind
    real(kind=realkind),intent(in)::safety_factor
    real(kind=realkind),intent(in),dimension(:,:,:)::uz,vz
    real(kind=realkind)::mspeed,dist,dt

#ifdef MPI_SRC
#include"mpif.h"
    integer:: ierr
    real(kind=realkind) tmp
#endif
    !m  this is an approximate conversion from degrees to meters and an approximation to non-uniform grid
    dist = safety_factor*min(RCAdomain%dlon,RCAdomain%dlat)*111000._realkind   
    mspeed = sqrt(maxval(uz*uz+vz*vz)) !max(maxval(abs(uz)),maxval(abs(vz)))!m/s 
    dt = real(Nseconds(ndtime),realkind)
#ifdef MPI_SRC
    tmp=mspeed
    call MPI_allreduce(tmp,mspeed,1,REALTYPE,mpi_max,localComm,ierr)
#endif

    if(dt<0._realkind)print *,ndtime
    compute_halo = ceiling(dt*mspeed/dist) + nslind
    if(compute_halo<1)then
       print *,dt,mspeed,dist,nslind
       call stop_program( 'compute_halo<1')
    endif
    if(compute_halo>maxHalo)then
       if(mype==0)then
          print *,'dt',dt
          print *,'maxspeed=',mspeed
          print *,'dist=',dist
          print *,'orderOfinterpolation=',nslind
          print *,compute_halo,'>',maxHalo
       endif
       call stop_program( 'halo>maxHalo')
    endif
  end function compute_halo


  

  subroutine sldynm(klon,klat,klev,ksvar,khalo_in,dynamic_halo,safety_factor,&
       kslpqi , kslint ,  kslind, &
       zdt , pfpara, &
       phis,   &             ! surface geopotential (input)
       ppsm, &
       pum,pvm,&
       pedotm, &
       ppsz, plnpsz,puz,pvz,ptz,pqz,psz,pedotz,psvarz, &
       ppsp, plnpsp,pup,pvp,ptp,pqp,psp,pedotp,psvarp, &
       ppp,lsl3d  ,epsg,ahyb,bhyb )
    use referenceParameters,only:sip0
    implicit none
    real(kind=realkind),intent(in)::epsg,safety_factor
    integer,intent(in)::klon,klat,klev,ksvar,kslpqi,kslind,khalo_in
    integer,intent(in):: kslint(kslpqi)
    real(kind=realkind),intent(in):: zdt,pfpara
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind),intent(in):: phis(klon,klat)
    real(kind=realkind),intent(in)::ppsm(klon,klat),&  
         pum(klon,klat,klev)  ,    pvm(klon,klat,klev), &
         pedotm(klon,klat,klev+1) 
    real(kind=realkind),intent(in)::ppsz(klon,klat)       , plnpsz(klon,klat)     , &
         puz(klon,klat,klev)  ,    pvz(klon,klat,klev), &
         ptz(klon,klat,klev)  ,    pqz(klon,klat,klev), &
         psz(klon,klat,klev) , pedotz(klon,klat,klev+1), psvarz(klon,klat,klev,ksvar)
    real(kind=realkind),intent(inout)::ppsp(klon,klat),plnpsp(klon,klat)     , &
         pup(klon,klat,klev),pvp(klon,klat,klev), &
         ptp(klon,klat,klev),pqp(klon,klat,klev),&
         psp(klon,klat,klev),pedotp(klon,klat,klev+1),psvarp(klon,klat,klev,ksvar)
    logical,intent(in)::lsl3d ,dynamic_halo
    real(kind=realkind),intent(in)::ppp(klon,klat,klev)

    real(kind=realkind)::hirtmp(klon,klat,klev)
    real(kind=realkind),allocatable,dimension(:,:,:)::ppx,psvarx
    real(kind=realkind),allocatable,dimension(:,:)::ppsx,plnpsx
    real(kind=realkind),allocatable,dimension(:,:,:)::pux,pvx,ptx,pqx,psx
    
    

    real(kind=realkind),allocatable,dimension(:,:,:)::palfa,pbeta,pgama
    integer:: jx,jy,jk,js 
    logical,save::initialized=.false.

    real(kind=realkind)::zcappa, zrsit0 , zdtrdx , zdtrdy , zffdt1 , zffdt2 , z1meps, &
          z1peps, z1mepsh, z1pepsh, zdtrdxm, zdtrdym, zffdt1m, &
          zdbk   , zgam1 , zgam2
    real(kind=realkind)::zalfh(klon,klat,klev), zbeth(klon,klat,klev)

    real(kind=realkind)::zpsum(klon,klat), &
         zw1(klon,klat),zw2(klon,klat), zw3(klon,klat), &
         zw4(klon,klat),zw5(klon,klat), zrhxy(klon,klat)


    real(kind=realkind)::zwk(klon,klat,klev)
    real(kind=realkind)::log10k
    integer::tmp
!    integer,allocatable,dimension(:,:)::kp,kq,kr
    type(slweights)::slw
    !compute khalo if needed
    if(dynamic_halo)then
       tmp = khalo
       khalo=compute_halo(pup,pvp,kslind,safety_factor)
       if(mype==0.and.tmp/=khalo)print *,'khalo is updated = ',khalo
    else
       khalo = khalo_in
    endif

    if(klon<khalo.or.klat<khalo)then
       print *, 'inner dimension too small sl2tim',klon,klat
       call stop_program('inner dimension too small sl2tim')
    endif


    allocate(palfa(klon,klat,klev))
    allocate(pbeta(klon,klat,klev))
    allocate(pgama(klon,klat,klev))
    allocate(ppx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(psvarx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(ppsx(klon,klat),plnpsx(klon,klat))
    allocate(pux(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev), &
         pvx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev), &
         ptx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev), &
         pqx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev), &
         psx(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))

    ppx = -666.0_realkind
!test

    log10k=log(100000._realkind)

    if(.not.initialized)then
       call allocate_combix_memory(klon,klat,klev)
       do jk = 1,klev+1
          etahalf(jk) = ahyb(jk)/sip0 + bhyb(jk)
       enddo
       do jk = 1,klev
          etafull(jk) = 0.5_realkind*(etahalf(jk)+etahalf(jk+1))
       enddo
       call preint(klev,tabetaf,tablinf,tabqdrf,tabcubf,etafull)
       call preint(klev+1,tabetah,tablinh,tabqdrh,tabcubh,etahalf)
       initialized = .true.
       !
       iwest  = max(npbpts+2-(idatastart-1),2) !since weight==1 at i<=npbpts+1
       isouth = max(npbpts+2-(jdatastart-1),2)
       !since weigth==1 at i,j>= k-npbpts
       ieast  = min(klon_global-npbpts-1-(idatastart-1),klon-1)
       inorth = min(klat_global-npbpts-1-(jdatastart-1),klat-1)
    endif



    zcappa = rair/cpair
    zrsit0 = rair*sit0
    zdtrdx = 0.5_realkind*zdt*rdlam*ra
    zdtrdy = 0.5_realkind*zdt*rdth*ra
    zffdt1 = 0.5_realkind*zdt*pfpara*0.25_realkind
    zffdt2 = zdt*pfpara*0.25_realkind

    z1meps = 1.0_realkind-epsg
    z1peps = 1.0_realkind+epsg

    z1mepsh = z1meps*0.5_realkind
    z1pepsh = z1peps*0.5_realkind
    zdtrdxm = z1meps*zdtrdx
    zdtrdym = z1meps*zdtrdy
    zffdt1m = z1meps*zffdt1

    zrhxy = ra*rhxu*rhyv
    zw1 = 0.0_realkind
    zw2 = 0.0_realkind
    zw3 = 0.0_realkind
    zw4 = 0.0_realkind
    zw5 = 0.0_realkind
    pup = pup - puz
    pvp = pvp - pvz
    ptp = ptp - ptz
    pqp = pqp - pqz
    psp = psp - psz
    psvarp = psvarp - psvarz

    palfa=0.0_realkind
    pbeta=0.0_realkind
    pgama=0.0_realkind


    call calpqr(klon,klat,klev,kslpqi,kslint, &
         khalo, &
         npbpts,&
         pum,pvm,pedotm,ppsm, &
         puz,pvz,pedotz,ppsz,&
         palfa,pbeta,pgama, &
         zalfh,zbeth, &
         lsl3d,kslind,zdt,ahyb,bhyb )

    pqx(1:klon,1:klat,:) = 0.5_realkind*pqp + pqz
    psx(1:klon,1:klat,:) = 0.5_realkind*psp + psz

    call slswap(pqx,klon,klat,klev,khalo)
    call slswap(psx,klon,klat,klev,khalo)


    slw = allocWeights(klon,klat)
    do jk=1,klev
       call prepbixint(klon,klat,klev,jk,kslind,khalo,npbpts, &
            slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk), &
            lsl3d )

       call dobixint(klon,klat,klev,jk,kslind,khalo,&
            pqx,zw3,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk), &
            lsl3d )
       call dobixint(klon,klat,klev,jk,kslind,khalo,&
            psx,zw4,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk),pbeta(:,:,jk), pgama(:,:,jk), &
            lsl3d )

       do jy=isouth,inorth
          do jx=iwest,ieast
             pqp(jx,jy,jk) = 0.5_realkind*pqp(jx,jy,jk) + zw3(jx,jy)
             psp(jx,jy,jk) = 0.5_realkind*psp(jx,jy,jk) + zw4(jx,jy)
          enddo
       enddo
    enddo



    do jy=1,klat
       do jx=1,klon
          plnpsp(jx,jy) = 0.0_realkind
          zw4   (jx,jy) = 0.0_realkind
       enddo
    enddo

    do jk=1,klev
       do jy=2,klat-1
          do jx=2,klon-1
             !    compute the divergence
             zw1(jx,jy) = zrhxy(jx,jy)*(rdlam*(puz(jx,jy,jk)*hyu(jx,jy)-puz(jx-1,jy,jk)*hyu(jx-1,jy)) &
                  + rdth*(pvz(jx,jy,jk)*hxv(jx,jy) - pvz(jx,jy-1,jk)*hxv(jx,jy-1)) )
             !    compute ptx = t(n) + (dt-/2)*{nt(n+1/2) + lt(n)}
             ptx(jx,jy,jk) = ptz(jx,jy,jk)+z1mepsh*ptp(jx,jy,jk)- &
                  z1meps*0.5_realkind*zdt*(sitau1(jk)*zw4(jx,jy)+&
                  sitau2(jk)*zw1(jx,jy))
             !    partial sum of divergence in zw4
             zw4(jx,jy) = zw4(jx,jy) + sidpk0(jk)*zw1(jx,jy)
          enddo
       enddo
    enddo

    !    compute ppx = [ln(ps)]' + (dt-/2)*{np(n+1/2) + lp(n)}

    do jk=1,klev
       do jy=2,klat-1
          do jx=2,klon-1
             ppx(jx,jy,jk) = plnpsz(jx,jy) - log10k + &
                  phis(jx,jy)/zrsit0 + z1meps*0.5_realkind*zdt*ppp(jx,jy,jk) - &
                  z1meps*0.5_realkind*zdt*zw4(jx,jy)/sip0
          enddo
       enddo
    enddo

    !    for the momentum equations
    zpsum = phis + zrsit0*plnpsz

    do jk=klev,1,-1
       zgam1 = sigam1(jk)
       zgam2 = sigam2(jk)
       zwk(:,:,jk) = zpsum + zgam1*ptz(:,:,jk)
       if (jk>1) then
          zpsum = zpsum + zgam2*ptz(:,:,jk)
       endif
    enddo

    call slswap(ptx,klon,klat,klev,khalo)
    call slswap(ppx,klon,klat,klev,khalo)
    do jk=1,klev
       call prepbixint (klon,klat,klev,jk,kslind,khalo,npbpts,&
            slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk), pbeta(:,:,jk), pgama(:,:,jk),&
            lsl3d )
       call dobixint (klon,klat,klev,jk,kslind,khalo,&
            ptx,zw1,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk), pbeta(:,:,jk), pgama(:,:,jk),&
            lsl3d )

       call dobixint ( klon,klat,klev,jk,kslind,khalo,&
            ppx,zw2,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            palfa(:,:,jk), pbeta(:,:,jk), pgama(:,:,jk), &
            lsl3d )

       zdbk  = bhyb(jk+1)-bhyb(jk)

       do jy=isouth,inorth
          do jx=iwest,ieast
             zw2(jx,jy) = zw2(jx,jy) + z1peps*0.5_realkind*zdt*ppp(jx,jy,jk)
             ptp(jx,jy,jk) = z1pepsh*ptp(jx,jy,jk) + zw1(jx,jy)
             hirtmp(jx,jy,jk)=zw2(jx,jy)*zdbk
          enddo
       enddo
    enddo
    
    do jk=1,klev
       do jy=isouth,inorth
          do jx=iwest,ieast
             plnpsp(jx,jy) = plnpsp(jx,jy)+hirtmp(jx,jy,jk)
          enddo
       enddo
    enddo

    do jy=isouth,inorth
       do jx=iwest,ieast
          plnpsp(jx,jy) = plnpsp(jx,jy) + log10k
       enddo
    enddo

    do js=1,ksvar
       do jk=1,klev
          do jy=2,klat-1
             do jx=2,klon-1
                psvarx(jx,jy,jk) = 0.5_realkind*psvarp(jx,jy,jk,js) +  psvarz(jx,jy,jk,js)
             enddo
          enddo
       enddo
       call slswap(psvarx,klon,klat,klev,khalo)
       do jk=1,klev
          call prepbixint( klon,klat,klev,jk,kslind,khalo,npbpts, &
               slw,&
               zalfh(:,:,jk),zbeth(:,:,jk),&
               palfa(:,:,jk), pbeta(:,:,jk), pgama(:,:,jk), &
               lsl3d )
          call dobixint( klon,klat,klev,jk,kslind,khalo,&
               psvarx,zw1,slw,&
               zalfh(:,:,jk),zbeth(:,:,jk),&
               palfa(:,:,jk), pbeta(:,:,jk), pgama(:,:,jk), &
               lsl3d )

          do jy=isouth,inorth
             do jx=iwest,ieast
                psvarp(jx,jy,jk,js) = 0.5_realkind*psvarp(jx,jy,jk,js) + zw1(jx,jy)
             enddo
          enddo
       enddo
    enddo

    do jk=1,klev
       do jy=2,klat-1
          do jx=2,klon-1
             pux(jx,jy,jk) = puz(jx,jy,jk) + z1mepsh*pup(jx,jy,jk) &
                  - zdtrdxm*rhxu(jx,jy) * (zwk(jx+1,jy,jk) - zwk(jx,jy,jk))&
                  + zffdt1m*( pvz(jx,jy  ,jk) + pvz(jx+1,jy  ,jk) &
                  + pvz(jx,jy-1,jk) + pvz(jx+1,jy-1,jk) )
             pvx(jx,jy,jk) = pvz(jx,jy,jk) + z1mepsh*pvp(jx,jy,jk)      &
                  - zdtrdym*rhyv(jx,jy) * (zwk(jx,jy+1,jk) - zwk(jx,jy,jk)) &
                  - zffdt1m*( puz(jx,jy  ,jk) + puz(jx-1,jy  ,jk) &
                  + puz(jx,jy+1,jk) + puz(jx-1,jy+1,jk) )

             pqx(jx,jy,jk) = palfa(jx,jy,jk)
             psx(jx,jy,jk) = pbeta(jx,jy,jk)
             ptx(jx,jy,jk) = pgama(jx,jy,jk)
          enddo
       enddo

    enddo

    call slswap(pux,klon,klat,klev,khalo)
    call slswap(pvx,klon,klat,klev,khalo)

    if (lsl3d) then
       call slswap(pqx,klon,klat,klev,khalo)
       call slswap(psx,klon,klat,klev,khalo)
       call slswap(ptx,klon,klat,klev,khalo)
    else
       call slswap(pqx,klon,klat,klev,khalo)
       call slswap(psx,klon,klat,klev,khalo)
    end if

    do jk=1,klev
       call intpqr(klon,klat,1,khalo,pqx(:,:,jk),zw1 )
       call intpqr(klon,klat,1,khalo,psx(:,:,jk),zw2 )

       if (lsl3d)then
          call intpqr(klon,klat,1,khalo,ptx(:,:,jk),zw3 )
       endif
       call prepbixint(klon,klat,klev,jk,kslind,khalo,npbpts,&
            slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            zw1,zw2,zw3, &
            lsl3d )
       call dobixint(klon,klat,klev,jk,kslind,khalo,&
            pux,zw4,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            zw1,zw2,zw3, &
            lsl3d )

       do jy=isouth,inorth
          do jx=iwest,ieast
             pup(jx,jy,jk) = z1pepsh*pup(jx,jy,jk) + zw4(jx,jy)
          enddo
       enddo

       call intpqr(klon,klat,2,khalo,pqx(:,:,jk),zw1 )
       call intpqr(klon,klat,2,khalo,psx(:,:,jk),zw2 )

       if (lsl3d)then
          call intpqr(klon,klat,2,khalo,ptx(:,:,jk),zw3 )
       endif
       
       call prepbixint(klon,klat,klev,jk,kslind,khalo,npbpts, &
            slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            zw1,zw2,zw3,&
            lsl3d )
       call dobixint(klon,klat,klev,jk,kslind,khalo,&
            pvx,zw4,slw,&
            zalfh(:,:,jk),zbeth(:,:,jk),&
            zw1,zw2,zw3,&
            lsl3d )

       do jy=isouth,inorth
          do jx=iwest,ieast
             pvp(jx,jy,jk) = z1pepsh*pvp(jx,jy,jk) + zw4(jx,jy)
          enddo
       enddo

    enddo
    do jy=isouth,inorth
       do jx=iwest,ieast
          plnpsp(jx,jy) = plnpsp(jx,jy) - phis(jx,jy)/zrsit0
       enddo
    enddo

    call deallocWeights(slw)
    deallocate(palfa)
    deallocate(pbeta)
    deallocate(pgama)
    deallocate(ppx)
    deallocate(psvarx)
    deallocate(ppsx,plnpsx)
    deallocate(pux,pvx,ptx,pqx,psx)

    return
  end subroutine sldynm


  subroutine preint(klev ,tabeta,tablin,tabqdr,tabcub,etalev)

    !     preint - preparation for vertical interpolation

    !     purpose:

    !     compute time independent arrays for vertical interpolation

    !     input parameters:
    !     
    !     klev      number of vertical levels
    !     maxtab    dimension of tabulated arrays
    !     etalev    array with level values
    !     etamin    minimum of eta-values in tables
    !     etainc    increment of eta-vaues in tables
    !     
    !     output parameters:
    !     
    !     tabeta    array of tabulated level indexes
    !     tablin    array of coefficients for linear interpolation
    !     tabqdr    array of coefficients for quadratic interpolation
    !     tabcub    array of coefficients for cubic interpolation

    !     j.e. haugen       hirlam      1992

    implicit none

    integer,intent(in):: klev
    integer,intent(out)::tabeta(maxtab,4)
    real(kind=realkind),intent(out)::tablin(maxtab,2),tabqdr(maxtab,3),tabcub(maxtab,4) 
    real(kind=realkind),intent(in)::etalev(klev)

    integer::jl,ilev,jk 
    real(kind=realkind)::zeta,zetam1,zetak,zmin,zetap1,zetam2,zdiff

    if(mype==0)then
       print *,' preint: klev ',klev
       print *,' maxtab       ',maxtab
       print *,' etamin       ',etamin
       print *,' etainc       ',etainc
       print *,' etamax       ',etamin+real(maxtab-1,realkind)*etainc
    endif

    !     compute table of "nearest level" to the interpolation point
    !     and table of coefficients.

    !     linear interpolation

    do jl = 1,maxtab
       zeta = etamin + real(jl-1,realkind)*etainc

       !     find nearest level jk,where etalev(jk)>=zeta
       ilev = klev
       do jk = 1,klev
          if (etalev(jk)>=zeta) then
             ilev = jk
             goto 100
          endif
       enddo
100    continue

       !     check that ilev is >= 2 and <= klev
       ilev = max(ilev,2)
       ilev = min(ilev,klev)

       !     table of "nearest level"
       tabeta(jl,1) = ilev

       !     prevent extrapolation below surface and above upper boundary
       zeta = max(zeta,0.0_realkind)
       zeta = min(zeta,1.0_realkind)

       !     table of coefficients
       zetam1 = etalev(ilev-1)
       zetak  = etalev(ilev)

       tablin(jl,1) =(zeta  -zetak)/(zetam1-zetak)
       tablin(jl,2) =(zeta  -zetam1)/(zetak -zetam1)
    enddo

    !     quadratic interpolation

    do jl = 1,maxtab
       zeta = etamin + real(jl-1,realkind)*etainc

       !     find nearest level jk,with minimum difference between
       !     zeta and etalev(jk)
       zmin = 1.0_realkind
       ilev = 1
       do jk = 1,klev
          zdiff = abs(zeta-etalev(jk))
          if (zdiff<zmin) then
             zmin = zdiff
             ilev = jk
          endif
       enddo

       !     table of "nearest level"
       tabeta(jl,2) = ilev

       !     check that ilev is >= 2 and <= klev-1
       ilev = max(ilev,2)
       ilev = min(ilev,klev-1)
       !     prevent extrapolation below surface and above upper boundary
       zeta = max(zeta,0.0_realkind)
       zeta = min(zeta,1.0_realkind)

       !     table of coefficients
       zetam1 = etalev(ilev-1)
       zetak  = etalev(ilev)
       zetap1 = etalev(ilev+1)

       tabqdr(jl,1) =((zeta  -zetak)*(zeta  -zetap1))/((zetam1-zetak)*(zetam1-zetap1))
       tabqdr(jl,2) =((zeta  -zetam1)*(zeta  -zetap1))/((zetak -zetam1)*(zetak -zetap1))
       tabqdr(jl,3) =((zeta  -zetam1)*(zeta  -zetak))/((zetap1-zetam1)*(zetap1-zetak))
       !     use linear interpolation where cubic not possible
       if (tabeta(jl,2)<2) then
          tabeta(jl,2) = 2
          tabqdr(jl,1) = tablin(jl,1)
          tabqdr(jl,2) = tablin(jl,2)
          tabqdr(jl,3) = 0.0_realkind
       endif
       if (tabeta(jl,2)>klev-1) then
          tabeta(jl,2) = klev-1
          tabqdr(jl,1) = 0.0_realkind
          tabqdr(jl,2) = tablin(jl,1)
          tabqdr(jl,3) = tablin(jl,2)
       endif
    enddo

    !     cubic interpolation

    do jl = 1,maxtab
       zeta = etamin + real(jl-1,realkind)*etainc

       !     find nearest level jk,where etalev(jk)>=zeta
       ilev = klev
       do jk = 1,klev
          if (etalev(jk)>=zeta) then
             ilev = jk
             goto 300
          endif
       enddo
300    continue

       !     table of "nearest level"
       tabeta(jl,3) = ilev

       !     check that ilev is >= 3 and <= klev-1
       ilev = max(ilev,3)
       ilev = min(ilev,klev-1)

       !     prevent extrapolation below surface and above upper boundary
       zeta = max(zeta,0.0_realkind)
       zeta = min(zeta,1.0_realkind)

       !     table of coefficients
       zetam2 = etalev(ilev-2)
       zetam1 = etalev(ilev-1)
       zetak  = etalev(ilev)
       zetap1 = etalev(ilev+1)

       tabcub(jl,1) = ((zeta  -zetam1)*(zeta  -zetak)*(zeta  -zetap1))/ &
            ((zetam2-zetam1)*(zetam2-zetak)*(zetam2-zetap1))
       tabcub(jl,2) = ((zeta  -zetam2)*(zeta  -zetak)*(zeta  -zetap1))/ &
            ((zetam1-zetam2)*(zetam1-zetak)*(zetam1-zetap1))
       tabcub(jl,3) = ((zeta  -zetam2)*(zeta  -zetam1)*(zeta  -zetap1))/ &
            ((zetak -zetam2)*(zetak -zetam1)*(zetak -zetap1))
       tabcub(jl,4) = ((zeta  -zetam2)*(zeta  -zetam1)*(zeta  -zetak))/ &
            ((zetap1-zetam2)*(zetap1-zetam1)*(zetap1-zetak))

       !     use linear interpolation where cubic not possible
       if (tabeta(jl,3)<3) then
          tabeta(jl,3) = 3
          tabcub(jl,1) = tablin(jl,1)
          tabcub(jl,2) = tablin(jl,2)
          tabcub(jl,3) = 0.0_realkind
          tabcub(jl,4) = 0.0_realkind
       endif
       if (tabeta(jl,3)>klev-1) then
          tabeta(jl,3) = klev-1
          tabcub(jl,1) = 0.0_realkind
          tabcub(jl,2) = 0.0_realkind
          tabcub(jl,3) = tablin(jl,1)
          tabcub(jl,4) = tablin(jl,2)
       endif
    enddo

    !     mixed cubic/linear interpolation

    do jl = 1,maxtab
       tabeta(jl,4) = tabeta(jl,3)
    enddo

    if(mype==0)print *,' preint done'

    return
  end subroutine preint


  subroutine prepbixint(klon,klat,klev,klevel,kint,khalo,kpbpts, & 
       slw, &
       palfh  , pbeth  ,&
       palfa  , pbeta  , pgama,   l3dim   )
    

    !  interpolation in the semi-lagrangian scheme
    !
    !  klon      number of gridpoints in the x-direction
    !  klat      number of gridpoints in the y-direction
    !  klev      number of vertical levels
    !  klevel    current vertical level
    !  kcall     = 1 - preparation of weights in addition to
    !                  the interpolation part
    !            = 2 - interpolation part only
    !  palfh     alfa hat
    !  pbeth     beta hat
    !  palfa     alfa
    !  pbeta     beta
    !  pgama     gama
    !  l3dim     true if 3-dimensional interpolation
    !
    !  output parameters:
    !
    !  externals:
    !
    !  verint    three-dimensional interpolation between gridpoints
    !  horint    two  -dimensional interpolation between gridpoints
    !

    implicit none

    integer,intent(in)::klon,klat,klev,klevel,kint,khalo
    integer,intent(in)::kpbpts

    logical,intent(in)::l3dim


    real(kind=realkind),intent(inout)::palfh(klon,klat)  , pbeth(klon,klat) !!!these could be local variables!!!
    real(kind=realkind),intent(inout)::palfa(klon,klat)  , pbeta(klon,klat)  , pgama(klon,klat) 

    real(kind=realkind)::zsum
    type(slweights),intent(inout)::slw
 
    integer:: jx, jy, ilev,io
    integer:: ipmin, ipmax, iqmin, iqmax
    real(kind=realkind):: zeps, z100 , zmax, zmin , zeta
    real(kind=realkind):: zr2 , z1ma , z1pa, z1mb , z1pb
    real(kind=realkind):: zr6 , z1ma2, z2ma, z1mb2, z2mb

    real(kind=realkind) zfix(4), zlim(4)
    data zfix /100.0_realkind,100.5_realkind,2*100.0_realkind/
    data zlim /1.000005_realkind,1.500005_realkind,2*2.000005_realkind/

    zeps = zlim(kint)

    slw%kp = -666
    slw%kq = -666
    slw%kr = -666
    
    !  check horizontal displacements close to the boundaries
    !  and do trajectory truncation (if kpbpts > 0)
    !  or determine position of interpolation molecule
    !  for which trajectory extrapolation starts (if kpbpts = 0)
    
    if (kpbpts>0)then
       !  in x-direction
       if (atleft)then 
          do jx = 1,klon!iwest,ieast
             zmax =  (real(jx-1,realkind)-zeps)
             do jy = 1,klat!isouth,inorth
                palfa(jx,jy) = min(palfa(jx,jy),zmax) !truncation
             enddo
          enddo
       else
          zsum = 0._realkind
          do jx = 1,klon!iwest,ieast
             zmax =  (real(jx+khalo-1,realkind)-zeps) !this is what is available
             do jy = 1,klat!isouth,inorth
                zsum = zsum  + abs(min(palfa(jx,jy),zmax) - palfa(jx,jy))
             enddo
          enddo
          if (zsum>0._realkind) then
             write(6,*)'truncation in halo zone in bixint',zsum
             call stop_program('truncation in halo zone in bixint')
          endif
       endif
       if (atright) then
          do jx = 1,klon!iwest,ieast
             zmin = -(real(klon-jx-1,realkind)-zeps)
             do jy = 1,klat!isouth,inorth
                palfa(jx,jy) = max(palfa(jx,jy),zmin)!truncation
             enddo
          enddo
       else
          zsum = 0._realkind
          do jx = 1,klon!iwest,ieast
             zmin = -(real(klon+khalo-jx,realkind)-zeps)
             do jy = 1,klat!isouth,inorth
                zsum = zsum  + abs( max(palfa(jx,jy),zmin) - palfa(jx,jy))
             enddo
          enddo
          if (zsum>0._realkind) then
             write(6,*)'truncation in halo zone in bixint',zsum
             call stop_program('truncation in halo zone in bixint')
          endif
       endif
       !  in y-direction
       if (atbase) then
          do jy = 1,klat!isouth,inorth
             zmax =  (real(jy-1,realkind)-zeps)
             do jx = 1,klon!iwest,ieast
                pbeta(jx,jy) = min(pbeta(jx,jy),zmax) !truncation
             enddo
          enddo
       else
          zsum = 0._realkind
          do jy = 1,klat!isouth,inorth
             zmax =  (real(jy+khalo-1,realkind)-zeps)
             do jx = 1,klon!iwest,ieast
                zsum = zsum  + abs(min(pbeta(jx,jy),zmax) - pbeta(jx,jy))
             enddo
          enddo
          if (zsum>0._realkind) then
             write(6,*)'truncation in halo zone in bixint',zsum
             call stop_program('truncation in halo zone in bixint')
          endif
       endif
       if (attop) then
          do jy = 1,klat!isouth,inorth
             zmin = -(real(klat-jy-1,realkind)-zeps)
             do jx = 1,klon!iwest,ieast
                pbeta(jx,jy) = max(pbeta(jx,jy),zmin) !truncation
             enddo
          enddo
       else
          zsum = 0._realkind
          do jy = 1,klat!isouth,inorth
             zmin = -(real(klat+khalo-jy,realkind)-zeps)
             do jx = 1,klon!iwest,ieast
                zsum = zsum  + abs(max(pbeta(jx,jy),zmin) - pbeta(jx,jy))
             enddo
          enddo
          if (zsum>0._realkind) then
             write(6,*)'truncation in halo zone in bixint',zsum
             call stop_program('truncation in halo zone in bixint')
          endif
       endif
    else
       if(atleft) then
          if (kint<=2) then
             ipmin = 3
          else
             ipmin = 4
          endif
       else
          if (kint<=2) then
             ipmin = 3 - khalo
          else
             ipmin = 4 - khalo
          endif
       endif
       if (atright) then
          if (kint<=1) then
             ipmax = klon - 1
          else
             ipmax = klon - 2
          endif
       else
          if (kint<=1) then
             ipmax = klon - 1 + khalo
          else
             ipmax = klon - 2 + khalo
          endif
       endif
       if (atbase) then
          if (kint<=2) then
             iqmin = 3
          else
             iqmin = 4
          endif
       else
          if (kint<=2) then
             iqmin = 3 - khalo
          else
             iqmin = 4 - khalo
          endif
       endif
       if (attop) then
          if (kint<=1) then
             iqmax = klat - 1
          else
             iqmax = klat - 2
          endif
       else
          if (kint<=1) then
             iqmax = klat - 1 + khalo
          else
             iqmax = klat - 2 + khalo
          endif
       endif
    endif
    
    !  preparation of horizontal weights -
    !  do trajectory extrapolation if kpbpts = 0
    
    z100 = zfix(kint)
    
    if (kpbpts>0) then
       do jy = 1,klat
          do jx = 1,klon
             slw%kp(jx,jy)    = int( palfa(jx,jy) + z100 ) - 100
             slw%kq(jx,jy)    = int( pbeta(jx,jy) + z100 ) - 100
             palfh(jx,jy) = palfa(jx,jy) - real(slw%kp(jx,jy),realkind)
             pbeth(jx,jy) = pbeta(jx,jy) - real(slw%kq(jx,jy),realkind)
             slw%kp(jx,jy)    = jx - slw%kp(jx,jy)
             slw%kq(jx,jy)    = jy - slw%kq(jx,jy)
          enddo
       enddo
    else
       do jy = 1,klat
          do jx = 1,klon
             slw%kp(jx,jy)    = jx - int( palfa(jx,jy) + z100 ) + 100
             slw%kq(jx,jy)    = jy - int( pbeta(jx,jy) + z100 ) + 100
             slw%kp(jx,jy)    = max(slw%kp(jx,jy),ipmin)
             slw%kp(jx,jy)    = min(slw%kp(jx,jy),ipmax)
             slw%kq(jx,jy)    = max(slw%kq(jx,jy),iqmin)
             slw%kq(jx,jy)    = min(slw%kq(jx,jy),iqmax)
             palfh(jx,jy) = palfa(jx,jy) - real(jx-slw%kp(jx,jy),realkind)
             pbeth(jx,jy) = pbeta(jx,jy) - real(jy-slw%kq(jx,jy),realkind)
          enddo
       enddo
    endif
    
    
    if (kint==1) then          !  linear interpolation
       do jy = 1,klat
          do jx = 1,klon
             slw%alfa(jx,jy,1) = palfh(jx,jy)
             slw%alfa(jx,jy,2) = 1.0_realkind - palfh(jx,jy)
             slw%beta(jx,jy,1) = pbeth(jx,jy)
             slw%beta(jx,jy,2) = 1.0_realkind - pbeth(jx,jy)
          enddo
       enddo
    elseif (kint==2) then        !  quadratic interpolation
       zr2 = 0.5_realkind
       do jy = 1,klat
          do jx = 1,klon
             z1ma = 1.0_realkind - palfh(jx,jy)
             z1pa = 1.0_realkind + palfh(jx,jy)
             z1mb = 1.0_realkind - pbeth(jx,jy)
             z1pb = 1.0_realkind + pbeth(jx,jy)
             slw%alfa(jx,jy,1) =  zr2*palfh(jx,jy)*z1pa
             slw%alfa(jx,jy,2) =  z1ma*z1pa
             slw%alfa(jx,jy,3) = -zr2*palfh(jx,jy)*z1ma
             slw%beta(jx,jy,1) =  zr2*pbeth(jx,jy)*z1pb
             slw%beta(jx,jy,2) =  z1mb*z1pb
             slw%beta(jx,jy,3) = -zr2*pbeth(jx,jy)*z1mb
          enddo
       enddo
    elseif (kint>=3) then   !  cubic and mixed cubic/linear interpolation
       zr2 = 0.5_realkind
       zr6 = 1.0_realkind/6.0_realkind
       do jy = 1,klat
          do jx = 1,klon
             z1ma2 = 1.0_realkind - palfh(jx,jy)*palfh(jx,jy)
             z2ma  = 2.0_realkind - palfh(jx,jy)
             z1mb2 = 1.0_realkind - pbeth(jx,jy)*pbeth(jx,jy)
             z2mb  = 2.0_realkind - pbeth(jx,jy)
             slw%alfa(jx,jy,1) = -zr6*palfh(jx,jy)*z1ma2
             slw%alfa(jx,jy,2) =  zr2*palfh(jx,jy)*(1.0_realkind + palfh(jx,jy))*z2ma
             slw%alfa(jx,jy,3) =  zr2*z1ma2*z2ma
             slw%alfa(jx,jy,4) = -zr6*palfh(jx,jy)*(1.0_realkind - palfh(jx,jy))*z2ma
             slw%beta(jx,jy,1) = -zr6*pbeth(jx,jy)*z1mb2
             slw%beta(jx,jy,2) =  zr2*pbeth(jx,jy)*(1.0_realkind + pbeth(jx,jy))*z2mb
             slw%beta(jx,jy,3) =  zr2*z1mb2*z2mb
             slw%beta(jx,jy,4) = -zr6*pbeth(jx,jy)*(1.0_realkind - pbeth(jx,jy))*z2mb
          enddo
       enddo
    endif
    
    if (l3dim) then     !  check vertical displacements close to the boundaries
       do jy = 1,klat!isouth,inorth
          do jx = 1,klon!iwest,ieast
             pgama(jx,jy) = max(etafull(1),pgama(jx,jy))
             pgama(jx,jy) = min(1.0_realkind,pgama(jx,jy))
          enddo
       enddo
       
       !  preparation of vertical weights
       if (kint==1) then    !  linear interpolation
          do jy = 1,klat!isouth,inorth
             do jx = 1,klon!iwest,ieast
                zeta = (pgama(jx,jy) - etamin)/etainc
                ilev = int(zeta + 1.5_realkind)
                slw%kr(jx,jy) = tabetaf(ilev,kint)
                slw%gama(jx,jy,1) = tablinf(ilev,1)
                slw%gama(jx,jy,2) = tablinf(ilev,2)
             enddo
          enddo
       elseif (kint==2) then         !  quadratic interpolation
          do jy = 1,klat!isouth,inorth
             do jx = 1,klon!iwest,ieast
                zeta = (pgama(jx,jy) - etamin)/etainc
                ilev = int(zeta + 1.5_realkind)
                slw%kr(jx,jy) = tabetaf(ilev,kint)
                slw%gama(jx,jy,1) = tabqdrf(ilev,1)
                slw%gama(jx,jy,2) = tabqdrf(ilev,2)
                slw%gama(jx,jy,3) = tabqdrf(ilev,3)
             enddo
          enddo
       elseif (kint>=3) then         !  cubic and mixed cubic/linear interpolation
          do jy = 1,klat!isouth,inorth
             do jx = 1,klon!iwest,ieast
                zeta = (pgama(jx,jy) - etamin)/etainc
                ilev = int(zeta + 1.5_realkind)
                slw%kr(jx,jy) = tabetaf(ilev,kint)
                slw%gama(jx,jy,1) = tabcubf(ilev,1)
                slw%gama(jx,jy,2) = tabcubf(ilev,2)
                slw%gama(jx,jy,3) = tabcubf(ilev,3)
                slw%gama(jx,jy,4) = tabcubf(ilev,4)
             enddo
          enddo
       endif
    endif
 

    return
  end subroutine prepbixint



  subroutine dobixint(klon,klat,klev,klevel,kint,khalo, &
       parg_bixint   , pres,slw, &
       palfh  , pbeth  ,&
       palfa  , pbeta  , pgama,   l3dim   )
    
    !  bixint - interpolation from parg_bixint to pres
    !
    !  interpolation in the semi-lagrangian scheme
    !
    !  klon      number of gridpoints in the x-direction
    !  klat      number of gridpoints in the y-direction
    !  klev      number of vertical levels
    !  klevel    current vertical level
    !  parg_bixint      field to be interpolated
    !  palfh     alfa hat
    !  pbeth     beta hat
    !  palfa     alfa
    !  pbeta     beta
    !  pgama     gama
    !  l3dim     true if 3-dimensional interpolation
    !
    !  output parameters:
    !
    !  pres      interpolated field
    !
    !  externals:
    !
    !  verint    three-dimensional interpolation between gridpoints
    !  horint    two  -dimensional interpolation between gridpoints
    !

    implicit none

    integer,intent(in)::klon,klat,klev,klevel,kint,khalo
    real(kind=realkind),intent(in)::parg_bixint(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,*)
    logical,intent(in)::l3dim
    real(kind=realkind),intent(inout)::pres(klon,klat) 
    real(kind=realkind),intent(inout)::palfh(klon,klat),pbeth(klon,klat)
    real(kind=realkind),intent(inout)::palfa(klon,klat),pbeta(klon,klat),pgama(klon,klat) 

    type(slweights),intent(in)::slw


    !  interpolation in verint and horint
    if(l3dim)then   !  3-dimensional case ( horizontal and vertical )
       call verint (klon,klat,klev,kint,khalo,&
            slw%kp,slw%kq,slw%kr, parg_bixint,pres,  palfh,pbeth,slw%alfa,slw%beta,slw%gama )
    else       !  2-dimensional case (horizontal)
       call horint( klon,klat,kint,khalo, slw%kp, slw%kq, &
            parg_bixint(:,:,klevel),pres, palfh,&
            slw%alfa,slw%beta )
    endif

    return
  end subroutine dobixint



  subroutine slswap(a,klon,klat,klev,khalo)

    implicit none
    integer,intent(in)::klon,klat,klev,khalo
    real(kind=realkind),intent(inout)::a(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev)
#ifdef MPI_SRC
#include"mpif.h"
    integer:: status(mpi_status_size)
    integer:: info
    integer,pointer::is,js
#endif
    if(nproc==1)return
    

#ifdef MPI_SRC
    is => sliside(khalo)
    js => sljside(khalo)

    if(pe_top/=-1)then   
       call mpi_recv(a(2-khalo,klat,1),1,js,pe_top,1,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(a(2-khalo,2,1),1,js,pe_base,1,localcomm,info)
       call mpi_recv(a,1,js,pe_base,2,localComm,status,info)
    endif
    if(pe_top/=-1)then
      call mpi_send(a(2-khalo,klat-khalo,1),1,js,pe_top,2,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(a(klon,2-khalo,1),1,is,pe_right,3,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(a(2,2-khalo,1),1,is,pe_left,3,localComm,info)
       call mpi_recv(a,1,is,pe_left,4,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(a(klon-khalo,2-khalo,1),1,is,pe_right,4,localComm,info)
    endif


#endif

    return

  end subroutine slswap

  subroutine intpqr(klon ,klat ,kcall ,khalo,  parg_intpqr ,pres)

    !     
    ! intpqr - interpolation from masspoints to velocitypoints
    !     
    !     j.e. haugen            hirlam
    !     
    !     purpose.
    !     --------
    !     
    !     the displacements parg_intpqr in the semi-lagrangian advection
    !     scheme is interpolated from masspoints to velocitypoints
    !     
    !   interface.
    !     ----------
    !     
    !     intpqr is called from sldynm
    !     
    !     input parameters:
    !     -----------------
    !     
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     kpbpts    number of passive boundary lines
    !     kcall     = 1 : interpolation to u velocity points
    !     = 2 : interpolation to v velocity points
    !     parg_intpqr      displacement
    !     
    !     output parameters:
    !     ------------------
    !     
    !     pres      interpolated field
    !     
    !     ---------------------------------------------------------------
    !     
    !         declaration of global parameters
    !     --------------------------------
    !     

    implicit none
    integer,intent(in):: klon,klat,kcall,khalo
    real(kind=realkind),intent(in):: parg_intpqr(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1)
    real(kind=realkind),intent(inout):: pres(klon,klat)

    !         local declarations
    integer:: i,j,i1,i2,j1,j2
    real(kind=realkind):: zah,za1,za2,za3,za4,zahl,za1l,za2l,za3l,zahr,za1r,za2r,za3r 

    !         weights for cubic interpolation alfa-hat = 0.5
    zah = 0.5_realkind
    za1 = -zah*(1.0_realkind-zah)*(1.0_realkind+zah)/6.0_realkind
    za2 = zah*(1.0_realkind+zah)*(2.0_realkind-zah)/2.0_realkind
    za3 = (1.0_realkind-zah)*(1.0_realkind+zah)*(2.0_realkind-zah)/2.0_realkind
    za4 = -zah*(1.0_realkind-zah)*(2.0_realkind-zah)/6.0_realkind

    !         weights for quadratic interpolation alfa-hat = 0.5
    zahl = 0.5_realkind
    za1l = 0.5_realkind*zahl*(1.0_realkind+zahl)
    za2l = (1.0_realkind-zahl)*(1.0_realkind+zahl)
    za3l = -0.5_realkind*zahl*(1.0_realkind-zahl)
    !         weights for quadratic interpolation alfa-hat = -0.5
    zahr = -0.5_realkind
    za1r = 0.5_realkind*zahr*(1.0_realkind+zahr)
    za2r = (1.0_realkind-zahr)*(1.0_realkind+zahr)
    za3r = -0.5_realkind*zahr*(1.0_realkind-zahr)

    if (kcall==1) then

       !         interpolation from mass points to u points
       if(atleft) then
          do j=isouth,inorth
             pres(iwest ,j) = za1l*parg_intpqr(iwest,j)+ za2l*parg_intpqr(iwest+1,j)+ za3l*parg_intpqr(iwest+2,j)
          enddo
          i1 = iwest + 1
       else
          i1 = iwest
       endif

       if(atright) then
          do j=isouth,inorth
             pres(ieast-1,j) = za1r*parg_intpqr(ieast -2,j)+ za2r*parg_intpqr(ieast -1,j)+ za3r*parg_intpqr(ieast ,j)
             pres(ieast ,j) = 1.5_realkind *parg_intpqr(ieast ,j) - 0.5_realkind *parg_intpqr(ieast -1,j)
          enddo
          i2 = ieast - 2
       else
          i2 = ieast
       endif

       do i=i1,i2               
          do j=isouth,inorth
             pres(i,j) = za1*parg_intpqr(i-1,j)+za2*parg_intpqr(i,j) + za3*parg_intpqr(i+1,j)+za4*parg_intpqr(i+2,j)
          enddo
       enddo

    elseif (kcall==2) then
       !          interpolation from mass points to v points
       if(atbase) then
          do i=iwest,ieast
             pres(i,isouth)=za1l*parg_intpqr(i,isouth)+ za2l*parg_intpqr(i,isouth+1)+za3l*parg_intpqr(i,isouth+2)
          enddo
          j1 = isouth + 1
       else
          j1 = isouth
       endif

       if(attop) then
          do i=iwest,ieast
             pres(i,inorth-1) = za1r*parg_intpqr(i,inorth-2)+ za2r*parg_intpqr(i,inorth -1)+ za3r*parg_intpqr(i,inorth)
             pres(i,inorth) = 1.5_realkind *parg_intpqr(i,inorth) - 0.5_realkind *parg_intpqr(i,inorth -1)
          enddo
          j2 = inorth - 2
       else
          j2 = inorth
       endif

       do j=j1,j2              
          do i=iwest,ieast
             pres(i,j) = za1*parg_intpqr(i,j-1) + za2*parg_intpqr(i,j)+ za3*parg_intpqr(i,j+1) + za4*parg_intpqr(i,j+2)
          enddo
       enddo
    else
       if(mype==0)then
!          print *,'invalid call to intpqr,kcall ',kcall
          call stop_program('invalid call to intpqr,kcall ')
       endif
    endif
    return
  end subroutine intpqr


  subroutine bivint(klon ,klat ,klev ,kint ,khalo, &
    tabeta,                                          &
    parg_bivint ,pres,                                            &
    tablin ,tabqdr ,tabcub,                                 &
    etalev)

    !
    !
    !
    !  bivint - vertical interpolation
    !
    !  purpose:
    !
    !  vertical interpolation from half levels to full levels
    !
    !  input parameters:
    !
    !  klon      number of gridpoints in x-direction
    !  klat      number of gridpoints in y-direction
    !  klev      number of vertical levels
    !  parg_bivint      field to be interpolated
    !
    !  output parameters:
    !
    !  pres      interpolated field
    !
    !  history:
    !
    !  j.e. haugen     hirlam    1992
    !
    !
    !
    implicit none

    integer,intent(in):: klon,klat,klev,kint,khalo 

    integer,intent(in)::tabeta(maxtab,4)

    real(kind=realkind),intent(in)::parg_bivint(klon,klat,klev)
    real(kind=realkind),intent(inout)::pres(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1)
    real(kind=realkind),intent(in)::tablin(maxtab,2) ,tabqdr(maxtab,3) ,tabcub(maxtab,4)

    real(kind=realkind),intent(in):: etalev
    integer::jlev,ilev,ilm1,jx,jy,ilp1,ilm2
    real(kind=realkind)::zeta,zgama1,zgama2,zgama3,zgama4

    if (kint==1) then
       !  linear interpolation
       zeta = (etalev - etamin)/etainc
       jlev = int(zeta + 1.5_realkind)
       ilev = tabeta(jlev,kint)
       ilm1 = ilev - 1
       zgama1 = tablin(jlev,1)
       zgama2 = tablin(jlev,2)

       do jy = 1,klat!isouth,inorth
          do jx = 1,klon!iwest,ieast
             pres(jx,jy) = zgama1*parg_bivint(jx,jy,ilm1)  + zgama2*parg_bivint(jx,jy,ilev)
          enddo
       enddo

    elseif (kint==2) then!  quadratic interpolation
       zeta = (etalev - etamin)/etainc
       jlev = int(zeta + 1.5_realkind)
       ilev = tabeta(jlev,kint)
       ilm1 = ilev - 1
       ilp1 = ilev + 1
       zgama1 = tabqdr(jlev,1)
       zgama2 = tabqdr(jlev,2)
       zgama3 = tabqdr(jlev,3)
       do jy = 1,klat!isouth,inorth
          do jx = 1,klon!iwest,ieast
             pres(jx,jy) = zgama1*parg_bivint(jx,jy,ilm1)+ zgama2*parg_bivint(jx,jy,ilev) + zgama3*parg_bivint(jx,jy,ilp1)
          enddo
       enddo
    elseif (kint>=3) then!  cubic and mixed cubic/linear interpolation
       zeta = (etalev - etamin)/etainc
       jlev = int(zeta + 1.5_realkind)
       ilev = tabeta(jlev,kint)
       ilm2 = ilev - 2
       ilm1 = ilev - 1
       ilp1 = ilev + 1
       zgama1 = tabcub(jlev,1)
       zgama2 = tabcub(jlev,2)
       zgama3 = tabcub(jlev,3)
       zgama4 = tabcub(jlev,4)

       do jy = 1,klat!isouth,inorth
          do jx = 1,klon!iwest,ieast
             pres(jx,jy) = zgama1*parg_bivint(jx,jy,ilm2)+ zgama2*parg_bivint(jx,jy,ilm1)+ &
                  zgama3*parg_bivint(jx,jy,ilev) + zgama4*parg_bivint(jx,jy,ilp1)
          enddo
       enddo
    endif
    return
  end subroutine bivint


  subroutine calpqr( klon   , klat   , klev   , kslpqi , kslint, &
       khalo, &
       kpbpts, &
       pum    , pvm    , pedotm , ppsm, &
       puz    , pvz    , pedotz , ppsz, &
       palfa  , pbeta  , pgama, &
       palfh  , pbeth  ,&
       l3dim  ,  kslind , zdt,ahyb,bhyb  )

    implicit none

    integer,intent(in)::klon,klat,klev,kslpqi
    integer,intent(in)::khalo
    integer,intent(in)::kpbpts
    integer,intent(in)::kslint(kslpqi)
    real(kind=realkind),intent(in)::zdt
    real(kind=realkind),intent(in),dimension(klev+1)::ahyb,bhyb
    logical,intent(in)::l3dim 
    integer,intent(in)::kslind
    real(kind=realkind),intent(in)::ppsm(klon,klat),pum(klon,klat,klev),pvm(klon,klat,klev),pedotm(klon,klat,klev+1)
    real(kind=realkind),intent(in)::ppsz(klon,klat),puz(klon,klat,klev),pvz(klon,klat,klev),pedotz(klon,klat,klev+1)

    real(kind=realkind),dimension(klon,klat,klev),intent(inout)::palfa,pbeta,pgama
    real(kind=realkind),dimension(klon,klat,klev),intent(inout)::palfh,pbeth 


    real(kind=realkind),allocatable,dimension(:,:,:)::parga,pargb,pargg,zarga,zargb,zargg

    real(kind=realkind),allocatable,dimension(:,:)::zw1,zw2,zw3,zw4,zw5,zw6,zwa,zwb,zwg
    real(kind=realkind),allocatable,dimension(:,:)::zza,zzb,zzg
    real(kind=realkind),allocatable,dimension(:,:,:)::zwu,zwv

    type(slweights)::slw
    integer:: jx,jy,jk,jiter
    real(kind=realkind)::zdtrdx,zdtrdy,zdak,zdbk,zdek
    real(kind=realkind)::aa,bb,cc,dd,ee,ff,dtrmnt,undt,unm1dt,vndt,vnm1dt,endt,enm1dt

    allocate(zw1(klon,klat))
    allocate(zw2(klon,klat))
    allocate(zw3(klon,klat))
    allocate(zw4(klon,klat)) 
    allocate(zw5(klon,klat))
    allocate(zw6(klon,klat))
    allocate(zwa(klon,klat))
    allocate(zwb(klon,klat))
    allocate(zwg(klon,klat))
    allocate(zza(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1))
    allocate(zzb(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1))
    allocate(zzg(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1))
    allocate(zwu(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(zwv(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(parga(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(pargb(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(pargg(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(zarga(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(zargb(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev))
    allocate(zargg(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev)) 
    
    if (kslpqi<1) then
       call stop_program('in calpqr: kslpqi must be >= 1')
    endif

    zdtrdx = zdt*rdlam*ra
    zdtrdy = zdt*rdth*ra

!!!!!!!!!!!!!!!!!!!!
    do jk = 1,klev
       call bivint(klon,klat,klev+1,kslind,khalo, &
            tabetah,pedotm,zzg, tablinh,tabqdrh,tabcubh, &
            etafull(jk))
       zdak = ahyb(jk+1)-ahyb(jk)
       zdbk = bhyb(jk+1)-bhyb(jk)
       zdek = zdt*(etahalf(jk+1)-etahalf(jk))
       do jy = 1,klat
          do jx = 1,klon
             zwu(jx,jy,jk) = zdtrdx*pum(jx,jy,jk)*rhxu(jx,jy) !@ time m (displacement)
             zwv(jx,jy,jk) = zdtrdy*pvm(jx,jy,jk)*rhyv(jx,jy)
             zargg(jx,jy,jk) = zdek*ppsm(jx,jy)*zzg(jx,jy)/(zdak+zdbk*ppsm(jx,jy))
          enddo
       enddo
    enddo
    call slswap(zwu,klon,klat,klev,khalo)
    call slswap(zwv,klon,klat,klev,khalo)
    do jk=1,klev
       call destag(klon,klat,1,khalo,zwu(:,:,jk),zarga(:,:,jk))!in grid point
       call destag(klon,klat,2,khalo,zwv(:,:,jk),zargb(:,:,jk))
    enddo
!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!
    do jk = 1,klev
       call bivint(klon,klat,klev+1,kslind,khalo,&
            tabetah, pedotz,zzg, tablinh,tabqdrh,tabcubh, etafull(jk))
       zdak = ahyb(jk+1)-ahyb(jk)
       zdbk = bhyb(jk+1)-bhyb(jk)
       zdek = zdt*(etahalf(jk+1)-etahalf(jk))
       do jy = 1,klat
          do jx = 1,klon
             zwu(jx,jy,jk) = zdtrdx*puz(jx,jy,jk)*rhxu(jx,jy) !@ time z (displacement)
             zwv(jx,jy,jk) = zdtrdy*pvz(jx,jy,jk)*rhyv(jx,jy)
             pargg(jx,jy,jk) = zdek*ppsz(jx,jy)*zzg(jx,jy)/(zdak+zdbk*ppsz(jx,jy))
          enddo
       enddo
    enddo
    call slswap(zwu,klon,klat,klev,khalo)
    call slswap(zwv,klon,klat,klev,khalo)
    do jk=1,klev
       call destag(klon,klat,1,khalo,zwu(:,:,jk),parga(:,:,jk))
       call destag(klon,klat,2,khalo,zwv(:,:,jk),pargb(:,:,jk))
    enddo
!!!!!!!!!!!!!!!!!!!!!


    dd =   0.50_realkind
    ee =   0.25_realkind
    ff =  -1.00_realkind

    aa =   1.0_realkind + dd + ee + 2.0_realkind*ff
    bb =  -0.5_realkind - dd - ff
    cc =   0.5_realkind - dd - 2.0_realkind*ee - 2.0_realkind*ff

    dtrmnt = ff*cc - dd*ee

    if (abs(dtrmnt)<1.0e-5_realkind) then
       call stop_program('abort in calpqr because f*c - d*e = 0.')
    endif

    do jk = 1,klev
       do jy = 1,klat
          do jx = 1,klon
             undt = parga(jx,jy,jk)
             unm1dt = zarga(jx,jy,jk)
             palfa(jx,jy,jk) = undt
             parga(jx,jy,jk) = cc*undt + dd*unm1dt !time weighted
             zarga(jx,jy,jk) = ee*undt + ff*unm1dt
             vndt = pargb(jx,jy,jk)
             vnm1dt = zargb(jx,jy,jk)
             pbeta(jx,jy,jk) = vndt
             pargb(jx,jy,jk) = cc*vndt + dd*vnm1dt
             zargb(jx,jy,jk) = ee*vndt + ff*vnm1dt
             endt = pargg(jx,jy,jk)
             enm1dt = zargg(jx,jy,jk)
             pgama(jx,jy,jk) = etafull(jk) - endt
             pargg(jx,jy,jk) = cc*endt + dd*enm1dt
             zargg(jx,jy,jk) = ee*endt + ff*enm1dt
          enddo
       enddo
    enddo

    !      call slswap6(parga,pargb,pargg,zarga,zargb,zargg,klon,klat,klev,khalo)
    call slswap(parga,klon,klat,klev,khalo)
    call slswap(pargb,klon,klat,klev,khalo)
    call slswap(pargg,klon,klat,klev,khalo)
    call slswap(zarga,klon,klat,klev,khalo)
    call slswap(zargb,klon,klat,klev,khalo)
    call slswap(zargg,klon,klat,klev,khalo)


    slw = allocWeights(klon,klat)
    do jk = 1,klev
       do jy = 1,klat
          do jx = 1,klon
             zza(jx,jy) = ( ( aa*ff-bb*ee)*parga(jx,jy,jk) + (-aa*dd+bb*cc)*zarga(jx,jy,jk) )/dtrmnt
             zzb(jx,jy) = ( ( aa*ff-bb*ee)*pargb(jx,jy,jk) + (-aa*dd+bb*cc)*zargb(jx,jy,jk) )/dtrmnt
             zzg(jx,jy) = ( ( aa*ff-bb*ee)*pargg(jx,jy,jk) + (-aa*dd+bb*cc)*zargg(jx,jy,jk) )/dtrmnt
          enddo
       enddo


       do jiter = 1,kslpqi
          call prepbixint( klon,klat,klev,jk,kslint(jiter),khalo,kpbpts, &
               slw,&
               palfh(:,:,jk),pbeth(:,:,jk), &
               palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk),&
               l3dim )
          call dobixint( klon,klat,klev,jk,kslint(jiter),khalo, &
               parga,zw1,slw,&
               palfh(:,:,jk),pbeth(:,:,jk), &
               palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk),&
               l3dim )

          call dobixint(klon,klat,klev,jk,kslint(jiter),khalo, &
               pargb,zw2,slw,&
               palfh(:,:,jk),pbeth(:,:,jk), &
               palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk),&
               l3dim )

          if (l3dim) then
             call dobixint (klon,klat,klev,jk,kslint(jiter),khalo, &
                  pargg,zw3,slw,&
                  palfh(:,:,jk),pbeth(:,:,jk),&
                  palfa(:,:,jk),pbeta(:,:,jk),pgama(:,:,jk),&
                  l3dim )
          endif
          do jy = 1,klat!isouth,inorth
             do jx = 1,klon!iwest,ieast
                zwa(jx,jy) = 2.0_realkind*palfa(jx,jy,jk)
                zwb(jx,jy) = 2.0_realkind*pbeta(jx,jy,jk)
                zwg(jx,jy) = 2.0_realkind*pgama(jx,jy,jk) - etafull(jk)
             enddo
          enddo

          call prepbixint( klon,klat,klev,jk,kslint(jiter),khalo,kpbpts, &
               slw,&
               palfh(:,:,jk),pbeth(:,:,jk),&
               zwa,zwb,zwg,&
               l3dim )
          call dobixint( klon,klat,klev,jk,kslint(jiter),khalo, &
               zarga,zw4,slw,&
               palfh(:,:,jk),pbeth(:,:,jk),&
               zwa,zwb,zwg,&
               l3dim )

          call dobixint(klon,klat,klev,jk,kslint(jiter),khalo,&
               zargb,zw5,slw,&
               palfh(:,:,jk),pbeth(:,:,jk),&
               zwa,zwb,zwg,&
               l3dim )

          if (l3dim) then
             call dobixint(klon,klat,klev,jk,kslint(jiter),khalo, &
                  zargg,zw6,slw,&
                  palfh(:,:,jk),pbeth(:,:,jk),&
                  zwa,zwb,zwg,  l3dim )
          endif

          do jy = isouth,inorth
             do jx = iwest,ieast
                palfa(jx,jy,jk) = zza(jx,jy) + zw1(jx,jy) + zw4(jx,jy)
                pbeta(jx,jy,jk) = zzb(jx,jy) + zw2(jx,jy) + zw5(jx,jy)
                if(l3dim)then
                   pgama(jx,jy,jk) = etafull(jk) - (zzg(jx,jy) + zw3(jx,jy) + zw6(jx,jy))
                endif
             enddo
          enddo
       enddo

    enddo
    call deallocWeights(slw)
    deallocate(zw1,zw2,zw3,zw4,zw5,zw6,zwa,zwb,zwg)
    deallocate(zza,zzb,zzg)
    deallocate(zwu,zwv)
    deallocate(parga)
    deallocate(pargb)
    deallocate(pargg)
    deallocate(zarga)
    deallocate(zargb)
    deallocate(zargg) 
    return
  end subroutine calpqr

  subroutine destag(klon,klat,kcall,khalo,velo,pres)
    !     destag - interpolation from velocitypoints to masspoints
    !     j.e. haugen            hirlam
    !     purpose.
    !     in the semi-lagrangian advection scheme the velocities
    !     are interpolated to masspoints to create arguments for
    !     the iterative procedure for computation of alfa and beta.
    !
    !     input parameters:
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     kcall     = 1 : interpolation from u-points
    !                 = 2 : interpolation from v-points
    !     velo      u- or v-velocity
    !
    !     output parameters:
    !     pres      interpolated field

    implicit none

    integer,intent(in)::klon,klat,kcall,khalo

    real(kind=realkind),intent(in):: velo(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1)
    real(kind=realkind),intent(inout)::pres(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1)

    integer::i,j,imin,imax,jmin,jmax

    real(kind=realkind):: zah,za1,za2,za3,za4,zahl,za1l,za2l,za3l,zahr,za1r,za2r,za3r

    !    weights for cubic interpolation alfa-hat = 0.5
    zah = 0.5_realkind
    za1 = -zah*(1.0_realkind-zah)*(1.0_realkind+zah)/6.0_realkind
    za2 = zah*(1.0_realkind+zah)*(2.0_realkind-zah)/2.0_realkind
    za3 = (1.0_realkind-zah)*(1.0_realkind+zah)*(2.0_realkind-zah)/2.0_realkind
    za4 = -zah*(1.0_realkind-zah)*(2.0_realkind-zah)/6.0_realkind

    !    weights for quadratic interpolation alfa-hat = 0.5
    zahl = 0.5_realkind
    za1l = 0.5_realkind*zahl*(1.0_realkind+zahl)
    za2l = (1.0_realkind-zahl)*(1.0_realkind+zahl)
    za3l = -0.5_realkind*zahl*(1.0_realkind-zahl)

    !    weights for quadratic interpolation alfa-hat = -0.5
    zahr = -0.5_realkind
    za1r = 0.5_realkind*zahr*(1.0_realkind+zahr)
    za2r = (1.0_realkind-zahr)*(1.0_realkind+zahr)
    za3r = -0.5_realkind*zahr*(1.0_realkind-zahr)

    if (kcall==1) then       !    interpolation of u velocity to mass points
       if(atleft) then
          do j=1,klat
             pres(1 ,j) = 1.5_realkind *velo(1,j) - 0.5_realkind *velo(2,j)
             pres(2 ,j) = za1l*velo(1,j)+ za2l*velo(2,j)   + za3l*velo(3,j)
          end do
          imin = 3
       else
          imin = 1
       end if

       if(atright) then
          do j=1,klat
             pres(klon,j) = za1r*velo(klon-2,j)+ za2r*velo(klon-1,j) + za3r*velo(klon ,j)
          end do
          imax = klon - 1
       else
          imax = klon
       end if

       do i=imin,imax
          do j=1,klat
             pres(i,j) = za1*velo(i-2,j) + za2*velo(i-1,j) &
                  + za3*velo(i ,j) + za4*velo(i+1,j)
          end do
       end do
    elseif (kcall==2) then   !    interpolation of v velocity to mass points
       if(atbase) then
          do i=1,klon
             pres(i,1) = 1.5_realkind *velo(i,1) - 0.5_realkind*velo(i,2)
             pres(i,2) = za1l*velo(i,1) + za2l*velo(i,2) + za3l*velo(i,3)
          end do
          jmin = 3
       else
          jmin = 1
       endif

       if(attop) then
          do i=1,klon
             pres(i,klat) = za1r*velo(i,klat-2) + za2r*velo(i,klat-1) + za3r*velo(i,klat)
          end do
          jmax = klat - 1
       else
          jmax = klat
       endif
       do j=jmin,jmax
          do i=1,klon
             pres(i,j) = za1*velo(i,j-2) + za2*velo(i,j-1) &
                  + za3*velo(i,j) + za4*velo(i,j+1)
          enddo
       enddo
    else
       if(mype==0)print *,'invalid call to destag,kcall=',kcall
       call stop_program('invalid call to destag,kcall=')
    endif
    return
  end subroutine destag


  subroutine verint (klon ,klat ,klev ,kint ,khalo, &
    kp  ,kq  ,kr,                                       &
    parg_verint ,pres,                                              &
    palfh ,pbeth,                                             &
    palfa ,pbeta ,pgama)

    !  verint - three dimensional interpolation
    !  three dimensional interpolation
    !  input parameters:
    !
    !  klon      number of gridpoints in x-direction
    !  klat      number of gridpoints in y-direction
    !  klev      number of vertical levels
    !  kint      type of interpolation
    !            = 1 - linear
    !            = 2 - quadratic
    !            = 3 - cubic
    !            = 4 - mixed cubic/linear
    !  kp        array of indexes for horizontal displacements
    !  kq        array of indexes for horizontal displacements
    !  kr        array of indexes for vertical   displacements
    !  parg_verint      array of arguments
    !  palfh     alfa hat
    !  pbeth     beta hat
    !  palfa     array of weights in x-direction
    !  pbeta     array of weights in y-direction
    !  pgama     array of weights in vertical direction
    !
    !  output parameters:
    !
    !  pres      interpolated field
    !
    implicit none

    integer,intent(in):: klon ,klat ,klev ,kint ,khalo
    integer,intent(in)::kp(klon,klat),kq(klon,klat),kr(klon,klat)
    real(kind=realkind),intent(in)::parg_verint(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1,klev) 
    real(kind=realkind),intent(in)::palfh(klon,klat),pbeth(klon,klat)
    real(kind=realkind),intent(in)::palfa(klon,klat,4),pbeta(klon,klat,4),pgama(klon,klat,4)
    real(kind=realkind),intent(inout)::pres(klon,klat)  

    integer::jx,jy,idx,idy
    integer:: ilev,ilm1,ilm2,ilp1
    real(kind=realkind)::z1mah,z1mbh

    if (kint==1) then

       !  linear interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             ilev = kr(jx,jy)
             ilm1 = ilev - 1

             pres(jx,jy) = pgama(jx,jy,1)*(                           &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-1,idy-1,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy-1,ilm1))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-1,idy ,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy ,ilm1))) &
                  + pgama(jx,jy,2)*(                            &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-1,idy-1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy-1,ilev))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-1,idy ,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy ,ilev)))
          enddo
       enddo

    elseif(kint==2) then
       !  quadratic interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             ilev = kr(jx,jy)
             ilm1 = ilev - 1
             ilp1 = ilev + 1

             pres(jx,jy) = pgama(jx,jy,1)*(                          &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-1,idy-1,ilm1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy-1,ilm1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy-1,ilm1))  &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-1,idy ,ilm1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy ,ilm1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy ,ilm1))  &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-1,idy+1,ilm1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy+1,ilm1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy+1,ilm1)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,2)*(             &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-1,idy-1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy-1,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy-1,ilev))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-1,idy ,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy ,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy ,ilev))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-1,idy+1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy+1,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy+1,ilev)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,3)*(            &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-1,idy-1,ilp1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy-1,ilp1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy-1,ilp1))  &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-1,idy ,ilp1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy ,ilp1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy ,ilp1))  &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-1,idy+1,ilp1)    &
                  + palfa(jx,jy,2)*parg_verint(idx ,idy+1,ilp1)    &
                  + palfa(jx,jy,3)*parg_verint(idx+1,idy+1,ilp1)))
          enddo
       enddo

    elseif (kint==3) then
       !  cubic interpolation
       !
       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             ilev = kr(jx,jy)
             ilm2 = ilev - 2
             ilm1 = ilev - 1
             ilp1 = ilev + 1

             pres(jx,jy) = pgama(jx,jy,1)*(                             &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-2,ilm2)       &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-2,ilm2)       &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-2,ilm2)       &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-2,ilm2))     &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilm2)       &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilm2)       &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilm2)       &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilm2))     &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilm2)       &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilm2)       &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilm2)       &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilm2))     &
                  + pbeta(jx,jy,4)*(palfa(jx,jy,1)*parg_verint(idx-2,idy+1,ilm2)       &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy+1,ilm2)       &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy+1,ilm2)       &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy+1,ilm2)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,2)*(             &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-2,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-2,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-2,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-2,ilm1))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilm1))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilm1))   &
                  + pbeta(jx,jy,4)*(palfa(jx,jy,1)*parg_verint(idx-2,idy+1,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy+1,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy+1,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy+1,ilm1)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,3)*(             &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-2,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-2,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-2,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-2,ilev))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilev))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilev))   &
                  + pbeta(jx,jy,4)*(palfa(jx,jy,1)*parg_verint(idx-2,idy+1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy+1,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy+1,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy+1,ilev)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,4)*(             &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-2,ilp1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-2,ilp1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-2,ilp1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-2,ilp1))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilp1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilp1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilp1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilp1))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilp1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilp1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilp1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilp1))   &
                  + pbeta(jx,jy,4)*(palfa(jx,jy,1)*parg_verint(idx-2,idy+1,ilp1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy+1,ilp1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy+1,ilp1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy+1,ilp1)))
          enddo
       enddo
    elseif (kint==4) then
       !  mixed cubic/linear interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             ilev = kr(jx,jy)
             ilm2 = ilev - 2
             ilm1 = ilev - 1
             ilp1 = ilev + 1

             z1mah = 1.0_realkind - palfh(jx,jy)
             z1mbh = 1.0_realkind - pbeth(jx,jy)

             pres(jx,jy) = pgama(jx,jy,1)*(                           &
                  pbeth(jx,jy)  *(palfh(jx,jy)  *parg_verint(idx-1,idy-1,ilm2)     &
                  + z1mah         *parg_verint(idx ,idy-1,ilm2))   &
                  + z1mbh         *(palfh(jx,jy)  *parg_verint(idx-1,idy ,ilm2)     &
                  + z1mah         *parg_verint(idx ,idy ,ilm2))) &
                  + pgama(jx,jy,4)*(                            &
                  pbeth(jx,jy)  *(palfh(jx,jy)  *parg_verint(idx-1,idy-1,ilp1)     &
                  + z1mah         *parg_verint(idx ,idy-1,ilp1))   &
                  + z1mbh         *(palfh(jx,jy)  *parg_verint(idx-1,idy ,ilp1)     &
                  + z1mah         *parg_verint(idx ,idy ,ilp1)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,2)*(             &
                  pbeta(jx,jy,1)*(palfh(jx,jy)  *parg_verint(idx-1,idy-2,ilm1)     &
                  + z1mah         *parg_verint(idx ,idy-2,ilm1))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilm1))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilm1)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilm1)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilm1)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilm1))   &
                  + pbeta(jx,jy,4)*(palfh(jx,jy)  *parg_verint(idx-1,idy+1,ilm1)     &
                  + z1mah         *parg_verint(idx ,idy+1,ilm1)))

             pres(jx,jy) = pres(jx,jy) + pgama(jx,jy,3)*(             &
                  pbeta(jx,jy,1)*(palfh(jx,jy)  *parg_verint(idx-1,idy-2,ilev)     &
                  + z1mah         *parg_verint(idx ,idy-2,ilev))   &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_verint(idx-2,idy-1,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy-1,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy-1,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy-1,ilev))   &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_verint(idx-2,idy ,ilev)     &
                  + palfa(jx,jy,2)*parg_verint(idx-1,idy ,ilev)     &
                  + palfa(jx,jy,3)*parg_verint(idx ,idy ,ilev)     &
                  + palfa(jx,jy,4)*parg_verint(idx+1,idy ,ilev))   &
                  + pbeta(jx,jy,4)*(palfh(jx,jy)  *parg_verint(idx-1,idy+1,ilev)     &
                  + z1mah         *parg_verint(idx ,idy+1,ilev)))
          enddo
       enddo

    endif

    return
  end subroutine verint

  subroutine horint(klon ,klat ,kint ,khalo, &
       kp  ,kq,                                      &
       parg_horint ,pres,                                    &
       palfh , &
       palfa ,pbeta)
    !
    !
    !
    !  horint - horizontal interpolation
    !
    !  purpose:
    !
    !  two dimensional (horizontal) interpolation
    !
    !  input parameters:
    !
    !  klon      number of gridpoints in x-direction
    !  klat      number of gridpoints in y-direction
    !  kint      type of interpolation
    !            = 1 - linear
    !            = 2 - quadratic
    !            = 3 - cubic
    !            = 4 - mixed cubic/linear
    !  kp        array of indexes for displacements
    !  kq        array of indexes for displacements
    !  parg_horint      array of arguments
    !  palfh     alfa hat
    !  palfa     array of weights in x-direction
    !  pbeta     array of weights in y-direction

    !  output parameters:

    !  pres      interpolated field

    !  j.e. haugen       hirlam      1992


    implicit none
    integer,intent(in)::klon,klat,kint ,khalo
    integer,intent(in)::kp(klon,klat) ,kq(klon,klat)
    real(kind=realkind),intent(in)::parg_horint(2-khalo:klon+khalo-1,2-khalo:klat+khalo-1) 
    real(kind=realkind),intent(inout)::pres(klon,klat) 
    real(kind=realkind),intent(in)::palfh(klon,klat) 
    real(kind=realkind),intent(in)::palfa(klon,klat,*),pbeta(klon,klat,*)

    integer:: jx,jy,idx,idy
    real(kind=realkind):: z1mah

    if (kint==1) then
       !  linear interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             pres(jx,jy) =                                      &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_horint(idx-1,idy-1)   &
                  + palfa(jx,jy,2)*parg_horint(idx ,idy-1)) &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_horint(idx-1,idy)   &
                  + palfa(jx,jy,2)*parg_horint(idx ,idy))
          enddo
       enddo

    elseif (kint==2) then
       !  quadratic interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             pres(jx,jy) =                                       &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_horint(idx-1,idy-1)    &
                  + palfa(jx,jy,2)*parg_horint(idx ,idy-1)    &
                  + palfa(jx,jy,3)*parg_horint(idx+1,idy-1))  &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_horint(idx-1,idy)    &
                  + palfa(jx,jy,2)*parg_horint(idx ,idy)    &
                  + palfa(jx,jy,3)*parg_horint(idx+1,idy))  &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_horint(idx-1,idy+1)    &
                  + palfa(jx,jy,2)*parg_horint(idx ,idy+1)    &
                  + palfa(jx,jy,3)*parg_horint(idx+1,idy+1))
          enddo
       enddo

    elseif (kint==3) then
       !  cubic interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             pres(jx,jy) =                                      &
                  pbeta(jx,jy,1)*(palfa(jx,jy,1)*parg_horint(idx-2,idy-2)   &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy-2)   &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy-2)   &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy-2)) &
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_horint(idx-2,idy-1)   &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy-1)   &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy-1)   &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy-1)) &
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_horint(idx-2,idy)   &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy)   &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy)   &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy)) &
                  + pbeta(jx,jy,4)*(palfa(jx,jy,1)*parg_horint(idx-2,idy+1)   &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy+1)   &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy+1)   &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy+1))
          enddo
       enddo

    elseif (kint==4) then
       !  mixed cubic/linear interpolation

       do jy = isouth,inorth
          do jx = iwest,ieast
             idx = kp(jx,jy)
             idy = kq(jx,jy)
             z1mah = 1.0_realkind - palfh(jx,jy)
             pres(jx,jy) =                                     &
                  pbeta(jx,jy,1)*(palfh(jx,jy)  *parg_horint(idx-1,idy-2)  &
                  + z1mah         *parg_horint(idx ,idy-2))&
                  + pbeta(jx,jy,2)*(palfa(jx,jy,1)*parg_horint(idx-2,idy-1)  &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy-1)  &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy-1)  &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy-1))&
                  + pbeta(jx,jy,3)*(palfa(jx,jy,1)*parg_horint(idx-2,idy)  &
                  + palfa(jx,jy,2)*parg_horint(idx-1,idy)  &
                  + palfa(jx,jy,3)*parg_horint(idx ,idy)  &
                  + palfa(jx,jy,4)*parg_horint(idx+1,idy))&
                  + pbeta(jx,jy,4)*(palfh(jx,jy)  *parg_horint(idx-1,idy+1)  &
                  + z1mah         *parg_horint(idx ,idy+1))
          enddo
       enddo

    endif

    return
  end subroutine horint

end module sldynm5
