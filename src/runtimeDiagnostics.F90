module runtimeDiagnostics
  use timetype
  use comhkp
  use decomp
  use RCAdomainMod
  use calendar
  use util
  implicit none
  private



  public diagnos,hilot,hilorh,higust,hiuv10

contains

  subroutine diagnos(klon,klat,klev,      &
       ptwodt,pps,pu,pv,pt,pq,ps,ppsm,phis,ahyb,bhyb)
 
    implicit none
#ifdef MPI_SRC
#include"mpif.h"
    real(kind=realkind)::slask(6)
    integer::info
    real(kind=realkind)::ztmp(6)
    real(kind=realkind)::zmax(2),ZRECV(2)
    integer::ind
#endif     
    real(kind=realkind)::pressureTrend        ! mean absolute pressure tendency (hpa/3 hours)
    real(kind=realkind)::meanPressure          ! area mean surface pressure
    real(kind=realkind)::meanHumid           ! area mean vertically integrated specific humid
    real(kind=realkind)::meanCw           ! area mean vertically integrated cloud water
    real(kind=realkind)::meanPotEnergy          ! area mean vertically integrated potential ener
    real(kind=realkind)::meanKinteticEnergy          ! area mean vertically integrated kinetic energy
    real(kind=realkind)::meanTotalEnergy          ! area mean vertically integrated total energy


    integer,intent(in):: klon, klat, klev
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind),intent(in)::ptwodt
    real(kind=realkind),intent(in):: pps(klon,klat),pu(klon,klat,klev), pv(klon,klat,klev), &
         pt(klon,klat,klev),  pq(klon,klat,klev),ps(klon,klat,klev), &
         ppsm(klon,klat),    phis(klon,klat) 

    !         declaration of work-space for runtimeDiagnostics
    real(kind=realkind)  zhxy(klon,klat), zdpk(klon,klat), zps(klon,klat),   &
         zpeb(klon,klat), zql (klon,klat),  zpe(klon,klat), &
         zke(klon,klat),  zabsdp(klon,klat), zsl (klon,klat)

    integer  i, j, k,    istart, istop, jstart, jstop 
    real(kind=realkind) zrg, zcappa, zrgasd, zcpd, z3hour,  zpbx,  &
         zdak, zdbk, zkes, zpes, zpsx, zqls, zsls


    real(kind=realkind)::zabsps
    real(kind=realkind)::zhxhy
    real(kind=realkind)::maxwind
    integer::imax(3)!index of max-velocity
      




    zcappa = 0.2857143_realkind                                                
    zrgasd = 287.04_realkind                                                   
    zcpd   = zrgasd/zcappa                                            

    do j = 1,klat
       do i = 1,klon
          zhxy(i,j) = 0._realkind
          zps(i,j) = 0._realkind
          zql(i,j) = 0._realkind !
          zsl(i,j) = 0._realkind !
          zpe(i,j) = 0._realkind !
          zke(i,j) = 0._realkind !
          zdpk(i,j) = 0._realkind
          zabsdp(i,j) = 0._realkind
       enddo
    enddo

    call jlimits(klon,klat,istart,istop,jstart,jstop)
    !         mass, humidity and energi                                         
    do j = 1,klat
       do i = 1,klon 
          zhxy(i,j) = 1._realkind/(rhxu(i,j)*rhyv(i,j))                     
          zps (i,j) = pps(i,j)*zhxy(i,j)                                    
          zpeb(i,j) = zps(i,j)*phis(i,j)
       enddo
    enddo

    zhxhy = sum(zhxy(istart:istop,jstart:jstop))
#ifdef MPI_SRC
    slask(1) = zhxhy
    call mpi_reduce(slask,zhxhy,1,REALTYPE,mpi_sum,0,localcomm,info)
#endif
    zhxhy = 1.0_realkind/zhxhy
    zrg   = 1.0_realkind/9.8066_realkind

    zpsx = sum(zps(istart:istop,jstart:jstop))
    zpbx = sum(zpeb(istart:istop,jstart:jstop))

    zqls = 0._realkind                                                         
    zsls = 0._realkind                                                         
    zpes = 0._realkind                                                         
    zkes = 0._realkind                                                         

    do k=1,klev                                                  
       zdak = ahyb(k+1) - ahyb(k)                                      
       zdbk = bhyb(k+1) - bhyb(k)                                      

       do j = 1,klat 
          do i = 1,klon 
             zdpk(i,j) = ( zdak + zdbk*pps(i,j) )*zhxy(i,j)
             zql(i,j) = pq(i,j,k)*zdpk(i,j)
             zsl(i,j) = ps(i,j,k)*zdpk(i,j)
             zpe(i,j) = pt(i,j,k)*zdpk(i,j)
          enddo
       enddo
       zqls = zqls + sum(abs(zql(istart:istop,jstart:jstop)))
       zsls = zsls + sum(abs(zsl(istart:istop,jstart:jstop)))
       zpes = zpes + sum(abs(zpe(istart:istop,jstart:jstop)))

       zke=0.0_realkind
       do j = 2,klat !jstop
          do i = 2,klon !istop
             zke(i,j) = zdpk(i,j)*( rhyv(i,j)*               &
                  ( (pu(i-1,j,k)*pu(i-1,j,k))*hyu(i-1,j)     &      
                  + (pu(i,  j,k)*pu(i,  j,k))*hyu(i,  j) )   &      
                  + rhxu(i,j)*                                 &
                  ( (pv(i,j-1,k)*pv(i,j-1,k))*hxv(i,j-1)     &      
                  + (pv(i,j  ,k)*pv(i,j  ,k))*hxv(i,j  ) ) )       
          enddo
       enddo
       zkes = zkes + sum(abs(zke(istart:istop,jstart:jstop)))
    enddo

    
#ifdef MPI_SRC
    ztmp(1) = zqls
    ztmp(2) = zsls
    ztmp(3) = zpes
    ztmp(4) = zkes
    ztmp(5) = zpsx
    ztmp(6) = zpbx
    slask=ztmp
    call mpi_reduce(slask,ztmp,6,REALTYPE,mpi_sum,0,localcomm,info)
    zqls = ztmp(1)
    zsls = ztmp(2)
    zpes = ztmp(3)
    zkes = ztmp(4)
    zpsx = ztmp(5)
    zpbx = ztmp(6)
#endif

    

    meanPressure = zrg*zhxhy*zpsx                                                 
    meanHumid  = zrg*zhxhy*zqls                                                 
    meanCw  = zrg*zhxhy*zsls       
    meanPotEnergy = zrg*zhxhy*(zpbx + zcpd*zpes)                                   
    meanKinteticEnergy = zrg*zhxhy*0.25_realkind*zkes                                            
    meanTotalEnergy = meanPotEnergy + meanKinteticEnergy        


    !     tendency of surface pressure                                      
    zabsdp = abs(pps-ppsm)
    
    zabsps = sum(zabsdp(istart:istop,jstart:jstop))     
#ifdef MPI_SRC
    slask(1)=zabsps
    call mpi_reduce(slask,zabsps,1,REALTYPE,mpi_sum,0,localcomm,info)
#endif
    z3hour = 3._realkind*3600._realkind/(ptwodt*0.5_realkind*100._realkind)                               
    pressureTrend = zabsps*z3hour/real(klon_global*klat_global,realkind)     


    maxwind = sqrt(maxval(pu*pu + pv*pv))
    imax = maxloc(pu*pu + pv*pv)
#ifdef MPI_SRC
    ind = imax(1)+(idatastart-1) + (imax(2)-1+(jdatastart-1))*klon_global+ &
         (imax(3)-1)*klon_global*klat_global
    zmax(1) = maxwind
    zmax(2) = real(ind,realkind)
    if(realkind==4)then
       call mpi_reduce(zmax,zrecv,1,mpi_2real,mpi_maxloc,0,localcomm,info)
    elseif(realkind==8)then
       call mpi_reduce(zmax,zrecv,1,mpi_2double_precision,mpi_maxloc,0,localcomm,info)
    endif
    maxwind = zrecv(1)
    ind = int(zrecv(2))
    call i2ijk(ind,klon_global,klat_global,imax(1),imax(2),imax(3))
#endif    

    if(mype==0) then
       write(6,*)
       write(6,*)' -------------------------------'
       write(6,'(1x,''abs(dps)/3h ='',f11.6,'' mb'')') pressureTrend
       write(6,'(1x,''maxwind   ='',f11.6,'' m/s in (x,y,lev)='',3i4)') &
            maxwind, imax(1),imax(2),imax(3)
       write(6,'(1x,''stps        ='',f30.12)') meanPressure
       write(6,'(1x,''stq         ='',f30.12)') meanHumid
       write(6,'(1x,''sts         ='',f30.12)') meanCw
       write(6,'(1x,''stpe        ='',f30.12)') meanPotEnergy
       write(6,'(1x,''stke        ='',f30.12)') meanKinteticEnergy
       write(6,'(1x,''stte        ='',f30.12)') meanTotalEnergy
    endif


    return
  end subroutine diagnos

  subroutine i2ijk(ind,klon,klat,i,j,k)
    implicit none
    integer,intent(in)::ind,klon,klat
    integer,intent(out)::i,j,k
    k = (ind-mod(ind,klon*klat))/(klon*klat)+1
    j = (ind-(k-1)*klon*klat-mod(ind-(k-1)*klon*klat,klon))/klon+1
    i = ind-(k-1)*klon*klat-(j-1)*klon
  end subroutine i2ijk


  subroutine hilot (klon, klat, t, tmax, tmin)  
    !
    !  purpose:  collect max and min temperatures during a period
    !
    ! written by ulf hansson, rossby centre, 980209
    ! rewritten by ulf hansson, rossby centre, 060728
    !
    !

    implicit none  
    integer :: klon, klat  
    real(kind=realkind)::t(klon,klat),tmax(klon,klat),tmin(klon,klat)  
    integer :: i, j, istart, istop, jstart, jstop  


    call jlimits (klon, klat, istart, istop, jstart, jstop)  
    do j = jstart, jstop  
       do i = istart, istop  
          !au040317                              skippa 99 gradk !!!!!!!!!!!!!
          if (t (i, j) >100._realkind) then  
             if (t (i, j) >tmax (i, j) ) tmax (i, j) = t (i, j)  
             if (t (i, j) <tmin (i, j) ) tmin (i, j) = t (i, j)  
          else  
             tmax (i, j) = t (i, j)  
             tmin (i, j) = t (i, j)  
          endif
       enddo
    enddo
    !
    !-----------------------------------------------------------------------
    !
    return  
  end subroutine hilot

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine higust (klon, klat, gustest, gustlo, gusthi, gustucest, &
       gustestmax, gustlomax, gusthimax, gustucestmax)
    !
    !  purpose:  collect max wind gust speed during a period
    !
    ! written by ulf hansson, rossby centre, 050920
    ! modified by ulf hansson, rossby centre, 091016, thanks marco
    !
    !

    implicit none  

    integer :: klon, klat  
    real(kind=realkind)::gustest(klon, klat),gustlo(klon, klat),gusthi(klon,klat),gustucest(klon,klat)
    real(kind=realkind) :: gustestmax (klon, klat),  gustlomax (klon, klat), gustucestmax(klon, klat)
    real(kind=realkind)::gusthimax(klon, klat)
    integer :: i, j, istart, istop, jstart, jstop  

    call jlimits (klon, klat, istart, istop, jstart, jstop)  
    !
    ! initialise
    ! no
    ! done in gemini
    !
    ! the work
    !
    do j = jstart, jstop  
       do i = istart, istop  
          if (gustest (i, j) >gustestmax (i, j) ) then  
             gustestmax (i, j) = gustest (i, j)  
          endif
          if (gustlo (i, j) >gustlomax (i, j) ) then  
             gustlomax (i, j) = gustlo (i, j)  
          endif
          if (gusthi (i, j) >gusthimax (i, j) ) then  
             gusthimax (i, j) = gusthi (i, j)  
          endif
          if (gustucest (i, j) >gustucestmax (i, j) ) then  
             gustucestmax (i, j) = gustucest (i, j)  
          endif
       enddo
    enddo
    !
    !-----------------------------------------------------------------------
    !
    return  
  end subroutine higust

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine hilorh (klon, klat, rh, rhmax, rhmin)  
    !
    !  purpose:  collect max and min rh2m or other during a period
    !
    ! written by ulf hansson, rossby centre, 980209
    ! rewritten by ulf hansson, rossby centre, 060728
    !
    !

    implicit none  
    integer :: klon, klat  
    real(kind=realkind) :: rh (klon, klat), rhmax (klon, klat), rhmin (klon, klat)  
    !
    !-----------------------------------------------------------------------
    !
    integer :: i, j, istart, istop, jstart, jstop  
    !
    !-----------------------------------------------------------------------
    !
    call jlimits (klon, klat, istart, istop, jstart, jstop)  
    !
    ! the work
    !
    do j = jstart, jstop  
       do i = istart, istop  
          if (rh (i, j) >rhmax (i, j) ) rhmax (i, j) = rh (i, j)  
          if (rh (i, j) <rhmin (i, j) ) rhmin (i, j) = rh (i, j)  
       enddo
    enddo
    !
    !
    !-----------------------------------------------------------------------
    !
    return  
  end subroutine hilorh

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine hiuv10(klon,klat,u,v,uvmax,intime)

    !
    !  purpose:  collect max wind speed during a period
    !            dump and reinitialise at the end of the period
    !
    ! written by ulf hansson, rossby centre, 980223
    !
    !

    implicit none  
    integer :: klon, klat  
    type(time),intent(in)::intime
    real(kind=realkind),intent(in):: u(klon, klat), v(klon, klat)
    real(kind=realkind),intent(out)::uvmax(klon, klat)  
    real(kind=realkind) :: sum  
    real(kind=realkind) :: z, zmax, zmin, savzmax  
    real(kind=realkind) :: rnopts, speed, uij, vij  
    real(kind=realkind) , allocatable::uvmax_global (:, :)  
    integer :: i, j, istart, istop, jstart, jstop  
    integer :: nopts  
    type(time),save::savtime
    integer :: savi, savj  

    integer :: putstep  

    logical :: lfirst, ltest  
    data lfirst / .true. /  
    data ltest / .false. /  

    data putstep / 0 /  
    save lfirst, ltest, putstep  

    save savi, savj, savzmax  


    call jlimits (klon, klat, istart, istop, jstart, jstop)  

    if (lfirst) then  
       savzmax = - 999._realkind  
       if (mype==0.and.ltest) then  
          write (6,  * ) '               maximum 10m windspeed '  
          write (6,  * ) '     max       min                 mean'  
       endif
       lfirst = .false.  
    endif
    do j = jstart, jstop  
       do i = istart, istop  
          uij = u (i, j)  
          vij = v (i, j)  
          if (uij<98._realkind.and.vij<98._realkind) then  
             speed = sqrt (uij * uij + vij * vij)  
             if (speed>uvmax (i, j) ) uvmax (i, j) = speed  
          endif
          uvmax (i, j) = max (uvmax (i, j), 0._realkind)  
       enddo
    enddo

    putstep = putstep + 1  

    if(ltest)then  
       putstep = 0  
       allocate(uvmax_global (klon_global, klat_global) )  
       call colfld(0, uvmax_global, uvmax, klon, klat)  
       if (mype==0) then  
          zmin = 1.e32_realkind  
          zmax = - 1.e32_realkind  
          sum = 0._realkind  
          do i = 1, klon_global  
             do j = 1, klat_global  
                z = uvmax_global(i,j)  
                if(z<zmin) zmin = z  
                if(z>zmax)then  
                   zmax = z  
                   savi = i  
                   savj = j  
                   savtime = intime
                endif
                sum = sum + z  
             enddo
          enddo
          deallocate(uvmax_global)  
          if (zmax>savzmax) then  
             write (6, '(a,f10.5,a,i4,a,i4,a,i5,5i3)') 'max= ', &
                  zmax, ' at i=', savi, ', j=', savj, ' time: ',savtime 
             savzmax = zmax  
          endif
          rnopts = 1._realkind / real(klon_global * klat_global,realkind)  
          write (6, '(2f10.5,10x,1f10.5)') zmax, zmin, sum*rnopts
       endif
    endif

    return  
  end subroutine hiuv10





end module runtimeDiagnostics
