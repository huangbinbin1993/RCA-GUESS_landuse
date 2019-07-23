module bkf_mod
  use decomp,only:realkind
  implicit none
  private

  public bkfcall
contains
  subroutine bkfcall(nhor,nlev,kstart,kstop,nca, &
       dtime, zdx_kf,zdy_kf,orogsigm,            &
       pzdiv,pzfrtop,ppbltop,             & 
       pf,t,q,u,v,omega_kf,                      &
       dtdt_kfd,dqdt_kfd,dqcdt_kfd,                     &
       raincv_kf,snowcv_kf,kfd_base,kfd_top,         &
       umfd,udrd,udrdf,om_cpbl,                                     &
       qvap,qliq,qice,                                  &
       clddep,shal_cgj,&
       oro,gpot)

    use confys
    use modd_convpar
    use modd_convpar_shal
    implicit none

    !     variables en appel
    integer:: kstart, kstop, nhor, nlev
    integer:: nca(nhor),kfd_base(nhor),kfd_top(nhor)
    integer:: shal_cgj(nhor)
    integer:: ctrig(nhor),ppbltop(nhor)
    real(kind=realkind)::oro(nhor),om_cpbl(nhor)

    real(kind=realkind):: dtime
    real(kind=realkind):: t(nhor,nlev), q(nhor,nlev)
    real(kind=realkind):: u(nhor,nlev), v(nhor,nlev)
    real(kind=realkind):: omega_kf(nhor,nlev),pzdiv(nhor,nlev)
    real(kind=realkind):: pf(nhor,nlev),gpot(nhor,nlev) 
    real(kind=realkind):: zdx_kf(nhor),zdy_kf(nhor)
    real(kind=realkind):: qvap(nhor,nlev),qliq(nhor,nlev),qice(nhor,nlev)
    real(kind=realkind)::  dtdt_kfd(nhor,nlev), dqdt_kfd(nhor,nlev),          &
         dqcdt_kfd(nhor,nlev)
    real(kind=realkind):: raincv_kf(nhor),snowcv_kf(nhor)                
    real(kind=realkind):: umfd(nhor,nlev)
    real(kind=realkind):: udrd(nhor,nlev),udrdf(nhor,nlev)
    !cgj
    real(kind=realkind):: clddep(nhor)
    real(kind=realkind):: orogsigm(nhor),pzfrtop(nhor)

    integer,dimension(:),allocatable :: kfs_base,kfs_top,indexcv
    real(kind=realkind), dimension(:),  allocatable :: cape
    real(kind=realkind), dimension(:,:),allocatable :: dmfd,umfs,dtdt_kfs,dqdt_kfs,dqcdt_kfs
    real(kind=realkind), dimension(:,:),allocatable :: dqidt_kfs,dqidt_kfd,svap,sliq,sice
    !cgj
    real(kind=realkind), dimension(:,:),allocatable :: sem,sems,sed

    integer:: i, k, nk
    integer:: kbdia, ktdia, kice ,kch
    logical:: orefresh, odown, osettadj,ochconv
    real(kind=realkind),dimension(:),allocatable :: ptimec

    integer,dimension(:),allocatable :: kcount
    real(kind=realkind), dimension(:,:),allocatable :: rmvap,rmliq,rmice
    real(kind=realkind), dimension(:,:),allocatable :: wg,pzz,pprlflx,pprsflx,pch,pchten
    real(kind=realkind), dimension(:),  allocatable :: wrk,pdxdy
    !cgj
    integer,dimension(:),allocatable :: isink,pbltop
    real(kind=realkind), dimension(:),allocatable :: pddep,psdep,woro,convw,zfrtop

    real(kind=realkind):: chlc,chls,cpd
    !cgj
    integer:: l700
    real(kind=realkind):: es,qs,ztv,aliq,bliq,cliq,dliq,p700


    allocate ( kfs_base  (nhor) )
    allocate ( kfs_top   (nhor) )
    allocate ( indexcv   (nhor) )
    allocate ( cape      (nhor) )
    allocate ( dmfd      (nhor,nlev) )
    allocate ( umfs      (nhor,nlev) )
    allocate ( dtdt_kfs  (nhor,nlev) )
    allocate ( dqdt_kfs  (nhor,nlev) )
    allocate ( dqcdt_kfs (nhor,nlev) )
    allocate ( dqidt_kfs (nhor,nlev) )
    allocate ( dqidt_kfd (nhor,nlev) )
    allocate ( svap      (nhor,nlev) )
    allocate ( sliq      (nhor,nlev) )
    allocate ( sice      (nhor,nlev) )
    allocate ( wg        (nhor,nlev) )
    allocate ( pzz       (nhor,nlev) )
    allocate ( pprlflx   (nhor,nlev) )
    allocate ( pprsflx   (nhor,nlev) )
    allocate ( pch       (nhor,nlev) )
    allocate ( pchten    (nhor,nlev) )
    allocate ( rmvap     (nhor,nlev) )
    allocate ( rmliq     (nhor,nlev) )
    allocate ( rmice     (nhor,nlev) )
    !cgj
    allocate ( sem     (nhor,nlev) )
    allocate ( sems    (nhor,nlev) )
    allocate ( sed     (nhor,nlev) )
    !cgj
    allocate ( kcount    (nhor) )
    allocate ( isink     (nhor) )
    !      allocate ( ctrig     (nhor) )
    allocate ( wrk       (nhor) )
    allocate ( pdxdy     (nhor) )
    !cgj
    allocate ( pddep     (nhor) )
    allocate ( psdep     (nhor) )
    allocate ( woro      (nhor) )
    allocate ( convw     (nhor) )
    allocate ( zfrtop    (nhor) )
    allocate ( pbltop    (nhor) )
    allocate(ptimec(nhor))
    kch   = 1
    kice  = 1
    kbdia = 1
    ktdia = 1
    odown    = .true. 
    osettadj = .false. 
    ochconv  = .false. 
    !     constantes
    chlc=2.5008e+6_realkind !j/kg chaleur latente de condensation
    chls=2.8345e+6_realkind !j/kg chaleur latente de sublimation
    cpd =1004.5_realkind    !j/kg/k chaleur specifique air sec
    !cgj
    aliq=613.3_realkind
    bliq=17.502_realkind
    cliq=4780.8_realkind
    dliq=32.19_realkind
    l700=2
    !       
    !        initialiser champs 
    !        kcount = 0 a modifier si osettadj = true
    !        dans ce cas, = 0 au premier pas de temps seulement
    !
    !      do 895 i=1,nhor
    do 895 i=kstart,kstop
       if(nca(i)==0) orefresh = .true.
       pdxdy     (i) = zdx_kf(i) * zdy_kf(i)
       indexcv   (i) = 0 
       kfd_top   (i) = 0 
       kfd_base  (i) = 0 
       kfs_top   (i) = 0 
       kfs_base  (i) = 0 
       kcount    (i) = 0
       raincv_kf (i) = 0._realkind
       snowcv_kf (i) = 0._realkind
       cape      (i) = 0._realkind
       !cgj
       clddep    (i) = 0._realkind
       pddep	   (i) = 0._realkind
       psdep     (i) = 0._realkind
       woro      (i) = orogsigm(i)
       convw     (i) = om_cpbl(i)
       zfrtop    (i) = pzfrtop(i)
       pbltop    (i) = nlev-ppbltop(i)+1 !flip in vertical from nlev--1 to 1--nlev
       !cgj
       isink     (i) = 1
       ctrig     (i) = 0
895 enddo

    do k=1,nlev
       !         do 896 i=1,nhor
       do 896 i=kstart,kstop
          pprlflx   (i,k) = 0._realkind
          pprsflx   (i,k) = 0._realkind
          dqcdt_kfd (i,k) = 0._realkind
          dqidt_kfd (i,k) = 0._realkind
          dtdt_kfd  (i,k) = 0._realkind
          dqdt_kfd  (i,k) = 0._realkind
          umfd      (i,k) = 0._realkind
          dmfd      (i,k) = 0._realkind
          udrd      (i,k) = 0._realkind
          udrdf     (i,k) = 0._realkind
          dqcdt_kfs (i,k) = 0._realkind
          dqidt_kfs (i,k) = 0._realkind
          dtdt_kfs  (i,k) = 0._realkind
          dqdt_kfs  (i,k) = 0._realkind
          umfs      (i,k) = 0._realkind
          svap      (i,k) = 0._realkind
          sliq      (i,k) = 0._realkind
          sice      (i,k) = 0._realkind
          pch       (i,k) = 0._realkind
          pchten    (i,k) = 0._realkind
          !cgj
          !cgj  calculate static energy terms to filter convective points
          !cgj
          es=aliq*exp((bliq*t(i,k)-cliq)/(t(i,k)-dliq))
          qs=epsilo*es/(pf(i,k)-es)
          ztv=t(i,k)*(1.0_realkind+0.608_realkind*q(i,k))
          sed(i,k)=(cpd*ztv)+gpot(i,k)  !cp*tv*gz
          sems(i,k)=sed(i,k)+chlc*qs
          sem(i,k)=sed(i,k)+chlc*q(i,k)
          !cgj
          p700=pf(i,nlev)-7.e4
          if(pf(i,k)<p700 .and. pf(i,k)<25000._realkind)l700=k
          !cgj
          !--------------------------------------------------------------------
          !         ajouter traceurs diffuses d'eau liquide et solide a la temp 
          !         et l'humidite 
          !         **** a enlever si eau disponible dans le modele****
          !
          t(i,k) = t(i,k) - chlc*qliq(i,k)/cpd - chls*qice(i,k)/cpd
          q(i,k) = q(i,k) + qliq(i,k) + qice(i,k)

          qliq(i,k) = 0.0_realkind    !it has been set to zero in akfrak.f
          qice(i,k) = 0.0_realkind    !it has been set to zero in akfrak.f
          sliq(i,k) = 0.0_realkind
          sice(i,k) = 0.0_realkind

          !        calcul du rapport de melange  
          rmvap(i,k) = q(i,k) / (1._realkind - q(i,k))
          rmliq(i,k) = 0.0_realkind
          rmice(i,k) = 0.0_realkind

          !        caculate wg 
          !        wg(i,k) = omega_kf(i,k) / (-1.0*pf(i,k)*gravit/(rair*t(i,k)))
          wg(i,k) = omega_kf(i,k)                             ! jyj omega_kf is wg(m/s) already in akfrak

          !        caculate pzz
          pzz(i,k) = gpot(i,k) / gravit
896    enddo
    enddo
    !
    !cgj
    do 897 k=nlev,l700,-1
       !         do 897 i=1,nhor
       do  i=kstart,kstop
          if(isink(i)==1)then
             if(omega_kf(i,k)>0._realkind)isink(i)=0     !vert vel is upward 
          endif
       enddo
897 enddo
    !cgj
    !cgj   if air sinking everywhere in up to 250hpa. can only have/need
    !cgj   convection if we have a superadiabatic layer.
    !cgj   regions of ascent we require a moist unstable layer somewhere
    !
    do k=nlev,l700,-1
       do nk=k-1,l700,-1
          !          do i=1,nhor
          do i=kstart,kstop
             if(isink(i)==1)then
                if(sed(i,k)>sed(i,nk))ctrig(i)=1
             endif
          enddo
       enddo
    enddo
    !
    do k=nlev,l700,-1
       do nk=k-1,l700,-1
          !          do i=1,nhor
          do i=kstart,kstop
             if(isink(i)==0)then
                if(sem(i,k)>sems(i,nk))ctrig(i)=1
             endif
          enddo
       enddo
    enddo
    !cgj
    ! reverse the order
    call revert (pf       ,wrk,nhor,nlev)
    call revert (pzz      ,wrk,nhor,nlev)
    call revert (pzdiv    ,wrk,nhor,nlev)
    call revert (t        ,wrk,nhor,nlev)
    call revert (rmvap    ,wrk,nhor,nlev)
    call revert (rmliq    ,wrk,nhor,nlev)
    call revert (rmice    ,wrk,nhor,nlev)
    call revert (u        ,wrk,nhor,nlev)
    call revert (v        ,wrk,nhor,nlev)
    call revert (wg       ,wrk,nhor,nlev)
    !c    appel de la convection profonde
    
    ptimec = 1800._realkind   
    call deep_convection (nhor,nlev,kstart,kstop,kbdia,ktdia,        &
         dtime,kice,orefresh,odown,osettadj,                         &
         pf,pzz,pdxdy,ptimec,                                  &  !cgj
         indexcv,                                                    &
         t,rmvap,rmliq,rmice,u,v,wg,woro,oro,                        &
         kcount,ctrig,dtdt_kfd,dqdt_kfd,dqcdt_kfd,dqidt_kfd,         &
         raincv_kf,snowcv_kf,                                        &
         kfd_top,kfd_base,pprlflx,pprsflx,                           &
         umfd,dmfd,udrd,udrdf,cape,                                  &
         qvap,qliq,qice,                                             &
         ochconv,kch,pch,pchten,pddep) 
    !c    appel de la convection peu profonde
    ptimec = 5400._realkind
    call shallow_convection(nhor,nlev,kstart,kstop,kbdia,ktdia,  &
         dtime,kice,osettadj,ptimec(1),                             &
         indexcv,oro,om_cpbl,                                    &
         pf,pzz,t,rmvap,rmliq,rmice,wg,                          &
         dtdt_kfs,dqdt_kfs,dqcdt_kfs,dqidt_kfs,                  &
         kfs_top,kfs_base,umfs,                                  &
         svap,sliq,sice,                                         &
         ochconv,kch,pch,pchten,psdep)
    !------------------------------------------------------------------------
    ! inverse the input variables
    call revert (pf       ,wrk,nhor,nlev)
    call revert (pzz      ,wrk,nhor,nlev)
    call revert (pzdiv      ,wrk,nhor,nlev)
    call revert (t        ,wrk,nhor,nlev)
    call revert (rmvap    ,wrk,nhor,nlev)
    call revert (rmliq    ,wrk,nhor,nlev)
    call revert (rmice    ,wrk,nhor,nlev)
    call revert (u        ,wrk,nhor,nlev)
    call revert (v        ,wrk,nhor,nlev)
    call revert (wg       ,wrk,nhor,nlev)
    ! inverse the output variables
    call revert (dtdt_kfd ,wrk,nhor,nlev)
    call revert (dqdt_kfd ,wrk,nhor,nlev)
    call revert (dqcdt_kfd,wrk,nhor,nlev)
    call revert (dqidt_kfd,wrk,nhor,nlev)
    call revert (umfd     ,wrk,nhor,nlev)
    call revert (dmfd     ,wrk,nhor,nlev)
    call revert (udrd     ,wrk,nhor,nlev)
    call revert (udrdf    ,wrk,nhor,nlev)
    call revert (qvap     ,wrk,nhor,nlev)
    call revert (qliq     ,wrk,nhor,nlev)
    call revert (qice     ,wrk,nhor,nlev)

    call revert (dtdt_kfs ,wrk,nhor,nlev)
    call revert (dqdt_kfs ,wrk,nhor,nlev)
    call revert (dqcdt_kfs,wrk,nhor,nlev)
    call revert (dqidt_kfs,wrk,nhor,nlev)
    call revert (umfs     ,wrk,nhor,nlev)
    call revert (svap     ,wrk,nhor,nlev)
    call revert (sliq     ,wrk,nhor,nlev)
    call revert (sice     ,wrk,nhor,nlev)

    !     units of rain and snow from m/s to mm/s
    !      do 921 i=1,nhor
    do 921 i=kstart,kstop
       raincv_kf(i)=raincv_kf(i)*1000._realkind
       snowcv_kf(i)=snowcv_kf(i)*1000._realkind

       !     if deep convection occurs, shallow convection is not allowed
       !     so kfs_top and kfs_base will only be non-zero if kfd_top and
       !     kfd_base are identically zero.
       !     culbute des niveaux de base et sommet des nuages convectifs
       if(kfd_top (i)/=0.0_realkind) kfd_top (i)=nlev-kfd_top (i)
       if(kfd_base(i)/=0.0_realkind) kfd_base(i)=nlev-kfd_base(i)
       if(kfs_top (i)/=0.0_realkind) kfd_top (i)=nlev-kfs_top (i)   
       if(kfs_base(i)/=0.0_realkind) kfd_base(i)=nlev-kfs_base(i)  
       !cgj
       if(pddep(i)/=0.0_realkind) clddep(i)=pddep(i)

       !cgj
       if(indexcv (i)==1  ) then 
          shal_cgj(i) = 1
          if(clddep(i)==0.0_realkind .and. psdep(i)/=0.0_realkind) clddep(i)=psdep(i)
       endif
       !        shal_cgj(i)=indexcv (i)
921 enddo

    do  k=1,nlev
       do 922 i=kstart,kstop
          ! in this scm, t, q and cw are updated in apcond.f-->pcond.f
          !         t    (i,k) =     t(i,k) + dtdt_kfd(i,k) * dtime
          !         rmvap(i,k) = rmvap(i,k) + dqdt_kfd(i,k) * dtime
          !         q    (i,k) = rmvap(i,k) / (1. + rmvap(i,k))
          !
          dtdt_kfd (i,k) = dtdt_kfd (i,k) + dtdt_kfs(i,k)
          dqdt_kfd (i,k) = dqdt_kfd (i,k) + dqdt_kfs(i,k)
          dqcdt_kfd (i,k) =dqcdt_kfd (i,k) +dqcdt_kfs(i,k)
          dqidt_kfd (i,k) =dqidt_kfd (i,k) +dqidt_kfs(i,k)
          umfd (i,k) =     umfd (i,k) +     umfs(i,k)
          qvap (i,k) =     qvap (i,k) +     svap(i,k)
          qliq (i,k) =     qliq (i,k) +     sliq(i,k)
          qice (i,k) =     qice (i,k) +     sice(i,k)
          !convert mixing ratio into specific humidity
          qvap (i,k) = qvap(i,k) / (1._realkind + qvap(i,k))
          qliq (i,k) = qliq(i,k) / (1._realkind + qliq(i,k))
          qice (i,k) = qice(i,k) / (1._realkind + qice(i,k))
922    enddo
    enddo

    deallocate(ptimec)
    deallocate ( kfs_base  )
    deallocate ( kfs_top   )
    deallocate ( indexcv   )
    deallocate ( cape      )
    deallocate ( dmfd      )
    deallocate ( umfs      )
    deallocate ( dtdt_kfs  )
    deallocate ( dqdt_kfs  )
    deallocate ( dqcdt_kfs )
    deallocate ( dqidt_kfs )
    deallocate ( dqidt_kfd )
    deallocate ( svap      )
    deallocate ( sliq      )
    deallocate ( sice      )
    deallocate ( wg        )
    deallocate ( pzz       )
    deallocate ( rmvap     )
    deallocate ( rmliq     )
    deallocate ( rmice     )
    deallocate ( pprlflx   )
    deallocate ( pprsflx   )
    deallocate ( pch       )
    deallocate ( pchten    )
    deallocate ( kcount    )
    deallocate ( wrk       )
    deallocate ( pdxdy     )
    deallocate ( woro      )
    deallocate ( isink     )
    deallocate ( sed       )
    deallocate ( sem       )
    deallocate ( sems      )
    deallocate ( pddep     )
    deallocate ( psdep     )
    deallocate ( convw     )
    deallocate ( zfrtop    )
    deallocate ( pbltop    )

    return
  end subroutine bkfcall

subroutine revert (fin,wrk,ni,nk)
  implicit none
  integer ni,nk,nkh,k,i
  real(kind=realkind):: fin(ni,nk),wrk(ni)

  nkh=nk/2
  do  k=1,nkh
     do  i=1,ni
        wrk(i)        = fin (i,k)
        fin(i,k)      = fin (i,nk-k+1)
        fin(i,nk-k+1) = wrk (i)
     enddo
  enddo
  return
end subroutine revert



end module bkf_mod
