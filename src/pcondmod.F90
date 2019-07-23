module pcondmod
  use decomp,only:stop_program,realkind
  implicit none
  private

  public apcond
contains

  subroutine apcond (nhor,nlev,kstart,kstop,    &
       dqdtin,dtdtin,omf,dtime,    &
!cgj040711      q,cw,t,pf,dpf,fice,lesat,   &
       q,cw,t,pf,dpf,fice,         &
       totcov,oldcov,cucov,shalcld,&
       dtdt,dqdt,dcwdt,            &
       prcpst,preta,frland,frice,  &
       snowh,ps,dzh,gpot,pblh)

    !
    use confys
    implicit none
    !     
    !     input 0-d
    !     
    logical lesat
    integer nhor,nlev,kstart,kstop
    real(kind=realkind) dtime
    !
    !  input 2-d
    !
    integer jwmask(nhor,nlev)
    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),    &
         totcov(nhor,nlev),oldcov(nhor,nlev),         &
         omf(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev), &
         dqdtin(nhor,nlev),dtdtin(nhor,nlev),         &
         preta(nhor,nlev),                            &
         prc(nhor,nlev),cucov(nhor,nlev)
    !
    !  output 1-d
    !
    real(kind=realkind) prcpst(nhor),frland(nhor),frice(nhor)
    !
    !  output 2-d
    !
    real(kind=realkind) dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev)
    !
    ! work space:
    !
    !  here follows 0-d integer work numbers
    integer jl,jk
    !  here follows 2-d real work arrays
    real(kind=realkind) cme(nhor,nlev)      ! rate of cond-evap within the cloud
    real(kind=realkind) evapr(nhor,nlev)    ! rate of evaporation of falling precipitation (1/s)
    real(kind=realkind) prain(nhor,nlev)    ! rate of conversion of condensate to precipitation (1/s)
    real(kind=realkind) rmelt(nhor,nlev)    ! heating rate due to precip phase change (k/s) (disabled)
    real(kind=realkind) hlat
    !
    real(kind=realkind) ps(nhor)
    real(kind=realkind) snowh(nhor)
    !gjkfarm
    real(kind=realkind) shalcld(nhor,nlev)
    real(kind=realkind) dzh(nhor,nlev)
    real(kind=realkind) gpot(nhor,nlev)
    real(kind=realkind) fice(nhor,nlev)
    real(kind=realkind) pblh(nhor)
    !gjkfarm
    !gj

    call pcond(kstart,kstop,nhor,nlev,dtime, &
         dqdtin,dtdtin,omf, &
         q,cw,t,pf,dpf,fice,lesat, &
         totcov,oldcov,cucov, & 
         cme,evapr,prain,rmelt, &
         prc,prcpst,preta,frland, frice, &
         snowh,ps,dzh,gpot,pblh)
    !
    !gj   add in the tendency terms due to ls condensation
    !gj   these were formally set to zero at the start of acondens
    !gj   they therefore will contain only ls tends on leaving acondens
    do jk=1,nlev
       do jl=kstart,kstop
          !gj
          !gj
          !gj     modify hlat for snow/water single or double transition
          !gj
          if(t(jl,jk)<tmelt)then      !assume evaped precip is snow->water->vapour
             hlat=latvap+latice
          else
             hlat=latvap
          endif
          !
          dqdt(jl,jk)=dqdt(jl,jk)-cme(jl,jk)+evapr(jl,jk)
          dtdt(jl,jk)=dtdt(jl,jk)+hlat/cpair*(cme(jl,jk)-evapr(jl,jk))
          dcwdt(jl,jk)=dcwdt(jl,jk)+cme(jl,jk)-prain(jl,jk)
       enddo
    enddo
    !

    return
  end subroutine apcond



  subroutine pcond (kstart ,kstop  ,plond  ,plev  ,deltat, &
       qtend  ,ttend  ,omga   ,                &
       q      ,cwat   ,t      ,p     ,pdel,    &
       fice   ,lesat,                          &
       cldn   ,cldo   ,cucov  ,                &
       cme    ,evapr  ,prain  , rmelt,         &
       pcflx  ,precab , preta, frland,frice,   &
       snowh,ps,dzh,gpot,pblh)
    ! calculate the prognostic condensate amount and tendencies
    !     p. rasch and j.e. kristjansson, april 1997

    use confys
    use escom
    use condensmod
    use ccons
    implicit none
    integer,intent(in):: kstart,kstop,plond,plev
    real(kind=realkind),intent(in):: deltat                ! 2*dt
    real(kind=realkind),intent(in):: cldn(plond,plev) !new value of cloud fraction(fraction)
    real(kind=realkind),intent(in):: cldo(plond,plev) !old value of cloud fraction(fraction)
    real(kind=realkind),intent(in):: cwat(plond,plev)     ! cloud water (kg/kg)
    real(kind=realkind),intent(in):: omga(plond,plev)     ! vert pressure vel (pa/s)
    real(kind=realkind),intent(in):: p(plond,plev)        ! pressure          (k)
    real(kind=realkind),intent(in):: pdel(plond,plev)     ! pressure thickness (pa)
    real(kind=realkind),intent(in):: q(plond,plev)        ! water vapor (kg/kg)
    real(kind=realkind),intent(in):: qtend(plond,plev)    ! mixing ratio tend (kg/kg/s)
    real(kind=realkind),intent(in):: t(plond,plev)        ! temperature       (k)
    real(kind=realkind),intent(in):: ttend(plond,plev)    ! temp tend         (k/s)
    real(kind=realkind),intent(in):: pcflx(plond,plev)!convective precip level by level(kg/m2/s)

    !  output 1-d

    real(kind=realkind) precab(plond)       ! rate of precipitation (kg / (m**2 * s))

    !  output 2-d

    real(kind=realkind) rmelt(plond,plev)   ! heating rate due to precip phase change (k/s) (disabled)
    real(kind=realkind) prain(plond,plev)   ! rate of conversion of condensate to precipitation (1/s)
    real(kind=realkind) evapr(plond,plev)   ! rate of evaporation of falling precipitation (1/s)
    real(kind=realkind) cme(plond,plev)     ! rate of cond-evap within the cloud
    !vo +++ additional hirlam output
    real(kind=realkind) preta(plond,plev)   ! stratiformous precip level by level
    !gjkfarm
    real(kind=realkind) cucov(plond,plev)   ! stratiformous precip level by level
    !vo ---
    !---------------------------------------------------------------------
    ! internal variables
    !  here follows 0-d integer work numbers
    integer i                   ! work variable
    integer iter                ! no. iterations for precipitation calculation
    integer k                   ! work variable
    integer l                   ! work variable
    !  here follows 0-d real work numbers
    real(kind=realkind) conke                  ! rate of evaporation of precipitation:
    real(kind=realkind) denom                  ! work variable
    real(kind=realkind) dqsdp                  ! change in sat spec. hum. wrt pressure
    real(kind=realkind) dqsdt                  ! change in sat spec. hum. wrt temperature
    real(kind=realkind) grav
    real(kind=realkind) icwc                   ! in-cloud water content (kg/kg)
    real(kind=realkind) mincld                 ! a small cloud fraction to avoid / zero
    real(kind=realkind) omeps                  ! 1 minus epsilon
    real(kind=realkind) omsm                   ! a number just less than unity (for rounding)
    real(kind=realkind) pcme1                  ! work variable
    real(kind=realkind) prtmp                  ! work variable
    real(kind=realkind) qsn(plond), esn               ! work variable
    real(kind=realkind) qtmp, ttmp             ! work variable
    real(kind=realkind) tc                     ! crit temp of transition to ice
    real(kind=realkind) cp                     ! heat capacity for dry air
    real(kind=realkind) hlatv                  ! latent heat of vaporization
    real(kind=realkind) hlatf                  ! latent heat of freezing
    real(kind=realkind) epsqs                  ! epsilon
    !
    !  here follows 1-d real work numbers
    real(kind=realkind) cldm(plond)            ! mean cloud fraction over the time step
    real(kind=realkind) cldmax(plond)          ! max cloud fraction above
    real(kind=realkind) coef(plond)            ! conversion time scale for condensate to rain
    real(kind=realkind) cwm(plond)             ! cwat mixing ratio at midpoint
    real(kind=realkind) cwn(plond)             ! cwat mixing ratio at end
    real(kind=realkind) es(plond)              ! sat. vapor pressure
    real(kind=realkind) fice(plond,plev)            ! fraction of cwat that is ice
    real(kind=realkind) gamma(plond)           ! d qs / dt
    real(kind=realkind) iceab(plond)           ! rate of ice only from above
    real(kind=realkind) pcme(plond)            ! provisional condensation minus evaporation 
    real(kind=realkind) prprov(plond)          ! provisional value of precip at bottom of layer
    real(kind=realkind) qs(plond)              ! spec. hum. of water vapor
    real(kind=realkind) qtl(plond)             ! tendency which would saturate the grid box in deltat
    real(kind=realkind) relhum(plond)          ! relative humidity
    real(kind=realkind) prect(plond)           ! rate of precipitation including convection (kg / (m**2 * s))

    !  here follows 2-d real work numbers
    real(kind=realkind) cldnt(plond,plev)      ! sane new cloud fraction    (temp , fraction)
    real(kind=realkind) cldot(plond,plev)      ! sane old cloud fraction    (temp , fraction)
    real(kind=realkind) convm(plond,plev)      ! moistening rate   (kg/kg/s)
    real(kind=realkind) qn(plond,plev)         ! mixing rat at end of time step ignoring condensate
    real(kind=realkind) qsp(plond,plev)        ! sat pt mixing ratio
    real(kind=realkind) tn(plond,plev)         ! temp at end of time step ignoring condensate
    real(kind=realkind) tsp(plond,plev)        ! sat pt temperature
    !gj
    real(kind=realkind) frland(plond)
    real(kind=realkind) frice(plond)
    real(kind=realkind) ps(plond)
    real(kind=realkind) snowh(plond)
    real(kind=realkind) cldtmp
    !gjecmwfevap
    real(kind=realkind) zwrk1,zwrk2,zwrk3,zwrk4
    !gj
    !gjfire
    real(kind=realkind) ppt_iter(plond)
    real(kind=realkind) dzh(plond,plev)
    real(kind=realkind) gpot(plond,plev)
    real(kind=realkind) pblh(plond)
    real(kind=realkind) wtthick
    real(kind=realkind) wtthick2
    real(kind=realkind) zbeta,tlow
    integer jintt
    logical lesat
    !gjfire


    lesat = .false.
    tlow = 250.16_realkind
    cp = cpair
    epsqs = epsilo
    hlatv = latvap
    hlatf = latice
    omeps = 1.0_realkind - epsqs
    mincld = 1.e-14_realkind

    grav = 9.80664_realkind
#if ( WORDSIZE == 8 )
    omsm = 0.99999999_realkind
#else
    omsm = 0.99999_realkind
#endif
    !
    !     number of times to iterate the precipitation calculation
    iter = 2

    !     constant for computing rate of evaporation of precipitation:
    conke = 1.e-5_realkind 
    !
    ! initialize a few single level fields
    do i = kstart,kstop
       precab(i) = 0._realkind
       prect(i) = 0._realkind
       iceab(i) = 0._realkind		! latent heat of precip above
       cldmax(i) = 0._realkind
    enddo
    !
    ! initialize some multi-level fields
    do k = 1,plev
       do i = kstart,kstop
          qn(i,k) = q(i,k) + deltat*qtend(i,k)
          tn(i,k) = t(i,k) + deltat*ttend(i,k)
          preta(i,k)=0._realkind
       enddo
    enddo

    ! find the saturation point for the provisional t and q without condensation
    !vo      t1 = tsecnd()

    call findsp (plond,plev,kstart,kstop,lesat, &
         qn,tn,cwat,cldn,p,tsp,qsp)
    !vo      timesp(lat) = timesp(lat) + tsecnd()-t1

    !
    do k = 1,plev
       call vqsatd (plond,kstart,kstop,lesat, &
            t(1,k),p(1,k),cwat(1,k),cldn(1,k),es,qs,gamma)
       do i = kstart,kstop
          relhum(i) = q(i,k)/qs(i)

          ! sane limits on cloud fractions
          !gjkfarm
          !gjkfarm            cldtmp=cldn(i,k)+cucov(i,k)
          cldtmp=cldn(i,k)
          !gjkfarm

          cldnt(i,k) = max(0._realkind,min(1._realkind,cldtmp,relhum(i)))
          cldot(i,k) = max(0._realkind,min(1._realkind,cldo(i,k),relhum(i)))

          ! the mean cloud fraction over the time step
          cldm(i) = max((cldnt(i,k) + cldot(i,k))*0.5_realkind,mincld)

          ! the max cloud fraction above this level
          cldmax(i) = max(cldmax(i), cldm(i))
       enddo


       dqsdp = 0.0_realkind
       do i = kstart,kstop

          ! fractions of ice at this level
          tc = t(i,k) - tmelt
          !gjorig            fice(i) = max(0.,min(-tc*0.05,1.0))
          !gj080304          jintt = max(1,1 + int((tmelt - t(i,k))/dttabl))
          !gj080304          jintt = min(jintt, 1750)
          !gj080304          fice(i)  = prbice(jintt)
          !gjrca3          fice(i) = max(0.,min(-tc*0.025,1.0))


          dqsdt = cp*gamma(i)/hlatv

          ! the moistening term
          convm(i,k) = qtend(i,k)- cldm(i)*dqsdt*ttend(i,k) &
               -cldm(i)*dqsdp*omga(i,k)
          !     
          ! amount of water which must be evaporated to change the cloud
          ! from fraction cldo to cldn keeping same incloud amount
          !gj200499            icwc = cwat(i,k)/max(cldot(i,k),0.0001)
          !gj200499            pcme1 = min(
          !gj200499     &                 (cldnt(i,k) - cldot(i,k))*icwc/deltat,
          !gj200499     &                 0.
          !gj200499     &                 )
          ! the above should be changed to allow either evap or cond to take place
          !gjnewcld  attempt to do this here ane include new cw formation if cloud growing
          !gj
          if(cldnt(i,k)<=cldot(i,k))then
             icwc = cwat(i,k)/max(cldot(i,k),0.0001_realkind)
             pcme1 = min((cldnt(i,k) - cldot(i,k))*icwc/deltat,0._realkind)
          else
             icwc = cwat(i,k)/max(cldnt(i,k),0.0001_realkind)
             pcme1 = max((cldnt(i,k) - cldot(i,k))*icwc/deltat,0._realkind)
          endif
          !gjnewcld
          !
          !
          ! first guess on q-e, have not yet checked to see that there
          ! is enough cloud water to support evaporation
          !     
          denom = 1._realkind + cldm(i) * hlatv/cp * dqsdt
          pcme(i) = (cldm(i)*convm(i,k) + pcme1) / denom


          ! calculate the cooling due to a phase change of the rainwater
          ! from above
          if (t(i,k)>=tmelt) then
             !              rmelt(i,k) =  -hlatf/cp*iceab(i)*grav/pdel(i,k) 
             rmelt(i,k) = 0._realkind
             iceab(i) = 0._realkind
          else
             rmelt(i,k) = 0._realkind
          endif

       enddo

       ! put reasonable bounds on the first guess
       do i = kstart,kstop
          qtmp = q(i,k) + ( qtend(i,k) - pcme(i) )*deltat
          ttmp = t(i,k) + deltat*(ttend(i,k)+rmelt(i,k)+hlatv/cp*pcme(i))
          !gj080304
          if(lesat)then
             if(cwat(i,k)<=0._realkind)then
                es(i) = esatw(ttmp)
             else
                es(i)=(cldn(i,k)*esat(ttmp))+((1._realkind-cldn(i,k))*esatw(ttmp))
             endif
          else
             es(i) = esat(ttmp)
          endif
          qs(i) = min(epsqs*es(i)/(p(i,k) - omeps*es(i)),1._realkind)
          !           qtl has the tend required to bring to saturation
          qtl(i) =  max((qsp(i,k) - qtmp)/deltat,0._realkind)

          ! make sure we condense enough to bring it back to saturation
          ! this guards (partially) against inconsistent cloud fractions
          ! if no cloud, then make sure we evaporate all the water
          if (qtmp>qsp(i,k)) then
             pcme(i) = pcme(i) + (qtmp-qsp(i,k))/deltat
          else
             if (cldm(i)>mincld) then
                !                    limit the estimated evaporation of cloud water
                !                    to that available 
                !                    or that which brings the box to saturation
                !                     pcme(i) = max(-cwat(i,k)/deltat,pcme(i),-qtl(i))
             else
                pcme(i) = - min(cwat(i,k)/deltat,qtl(i))
             endif
          endif
          pcme(i) = max(-cwat(i,k)/deltat,pcme(i),-qtl(i))*omsm

       enddo

       do i = kstart,kstop

          !           a provisional value of the cloud water from cond-evap at midpoint of timestep
          cwm(i) = max(cwat(i,k) + pcme(i)*deltat*0.5_realkind,0._realkind)

          !           a provisional value at endpoint 
          cwn(i) = max(cwat(i,k) + pcme(i)*deltat,0._realkind)

          !           a provisional value of precip
          prain(i,k) = 0.0_realkind

          !           move provisional cond-evap into final location
          cme(i,k) = pcme(i)
       enddo

       ! calculate the formation of precip. since this is a highly nonlinear
       ! calculation, we do it iteratively, using values from the midpoint of
       ! the time step
       !gjfire011202
       !gj
       do l = 1,iter

          do i = kstart,kstop
             !gj              prprov(i) = prect(i) + prain(i,k)*pdel(i,k)/grav
             wtthick=max(0._realkind,min(1._realkind,(dzh(i,k)/1000._realkind)**2._realkind))
             !gj               wtthick2=max(0.,min(1.,dzh(i,k)/1000.))
             prprov(i) = prect(i) + (prain(i,k)*pdel(i,k)/grav)*wtthick
          enddo
          !vo            t1 = tsecnd()
          call findmcnew (plond,plev,kstart,kstop,&
               k, prprov, t, p, cwm, cldm, cldmax,&
               fice, coef, frland,frice, snowh,ps,gpot,pblh)

          !vo            timemc(lat) = timemc(lat) + tsecnd()-t1


          !        calculate the precip rate
          do i = kstart,kstop
             if (cldm(i)>0._realkind) then
                !                 first predict the cloud water
                if (coef(i)>0._realkind) then
                   !gjorig                     cwn(i) = max(
                   !gjorig     &                             0.,
                   !gjorig     &                            (cwat(i,k)-cme(i,k)/coef(i))
                   !gjorig     &                              *exp(-coef(i)*deltat)
                   !gjorig     &                             + cme(i,k)/coef(i)
                   !gjorig     &                           )
                   cwn(i) = max(0._realkind,(cwat(i,k)-cldm(i)*cme(i,k)/coef(i)) &
                        *exp(-coef(i)*deltat)+cldm(i)*cme(i,k)/coef(i))
                   !gjnew
                else
                   cwn(i) = max(0._realkind,cwat(i,k) + cme(i,k)*deltat)
                endif
                !                 now back out the tendency
                if(coef(i)>0._realkind)then
                   prain(i,k) = max(0._realkind,cme(i,k) - (cwn(i)-cwat(i,k))/deltat )
                   !gj                   prain(i,k)=0.
                else
                   prain(i,k)=0._realkind
                endif
             else
                prain(i,k) = 0.0_realkind
                cwn(i) = 0._realkind
             endif

             !              update any remaining  provisional values
             cwm(i) = (cwn(i) + cwat(i,k))*0.5_realkind
          enddo              ! end of do i = kstart,kstop

       enddo                 ! end of do l = 1,iter


       ! we now have estimates of the condensation-evaporation terms, and
       ! the precipitation
       ! so given decent estimates of pcme, calculate provisional value
       ! of cloud water for evapr calculation

       do i = kstart,kstop
          qtmp = q(i,k) + (qtend(i,k) - cme(i,k))*deltat
          ttmp = t(i,k) + (ttend(i,k)+rmelt(i,k))*deltat + hlatv/cp*deltat*cme(i,k)
          !gj080304
          if(lesat)then
             if(cwat(i,k)<=0._realkind)then
                esn = esatw(ttmp)
             else
                esn=(cldn(i,k)*esat(ttmp))+ ((1._realkind-cldn(i,k))*esatw(ttmp))
             endif
          else
             esn = esat(ttmp)
          endif
          !gj080304
          qsn(i) = min(epsqs*esn/(p(i,k) - omeps*esn),1._realkind)
          qtl(i) = max((qsn(i) - qtmp)/deltat,0._realkind)
          relhum(i) = qtmp/qsn(i)

       enddo
       !
       ! evaporation of stratiformous rain, evaporation from convective rain
       ! is done in cuflx (no evaporation from snow yet)
       do i = kstart,kstop
          !gjrca3            evapr(i,k) = conke*(1. - cldot(i,k))*sqrt(precab(i))
          !gjrca3     &                        *(1. - min(relhum(i),1.))
          !           limit the evaporation to the amount which is entering the box
          !           or saturates the box
          prtmp = precab(i)*grav/pdel(i,k)
          !gj191007
          !gj191007  to run a very preliminary version of ecmwf precipitation evapoartion
          !gj191007  uncomment the 4 lines below between the comment lines cgjecmwf
          !gj191007  this is based on p87 eqn 6.43 of the ifs cy31r1 documentation.
          !gj191007  it is temporary code only
          !gjecmwf
          zwrk1=(1._realkind-cldot(i,k))*5.44e-4_realkind*(1._realkind-relhum(i))*qsn(i)
          zwrk4=precab(i)/5.9e-3_realkind
          evapr(i,k)=zwrk1*( zwrk4**0.577_realkind )
          !gjecmwf
          evapr(i,k) = min(evapr(i,k), prtmp, qtl(i))*omsm
       enddo
       !
       ! precipitation 
       do i = kstart,kstop
          prtmp = pdel(i,k) / grav *(prain(i,k) - evapr(i,k))
          iceab(i) = iceab(i) + fice(i,k)*prtmp
          precab(i) = precab(i) + prtmp
          !vo+++ additional hirlam output - precipitation level by level
          !gjorig            preta(i,k)= preta(i,k)+ prtmp
          preta(i,k)= precab(i)+ prtmp
          if (abs(preta(i,k))<1.e-14_realkind) then
             preta(i,k) = 0._realkind
          endif
          !vo---
          !au            prect(i) = prect(i) + prtmp + pcflx(i,k+1)
          !cj            prect(i) = prect(i) + prtmp + pcflx(i,k)
          prect(i) = prect(i) + prtmp
          if (abs(precab(i))<1.e-14_realkind) then
             precab(i) = 0._realkind
          endif
          if (abs(prect(i))<1.e-14_realkind) then
             prect(i) = 0._realkind
          endif
          !gjcliwanet
          !gj240302   3d rain rates converted to mixing ratio (g/kg)
          !gj3ddout            preta_liq(i,k)= (precab(i)*(1.-fice(i)))
          !gj280602     &                       *(grav/pdel(i,k))*deltat*1000.
          !gj3ddout            preta_ice(i,k)= (precab(i)*fice(i))
          !gj280602     &                       *(grav/pdel(i,k))*deltat*1000.
          !gjcliwanet
       enddo

    enddo                    ! level loop (k=1,plev)
    !
    return
  end subroutine pcond



  subroutine findsp(plond,plev,kstart,kstop,lesat, &
       q,t,cwat,cloud,p,tsp, qsp)

    !     find the saturation point for a given t and q
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q


    !-----------------------------------------------------------------------
    use confys
    use escom
    implicit none
    !     

    !     
    integer plond             ! number of horizontal gridpoints (nlon*nlat)
    integer plev              ! number of vertical levels
    integer kstart            ! starting gridpoint of the call
    integer kstop             ! final gridpoint of the call


    real(kind=realkind) q(plond,plev)        ! water vapor (kg/kg)
    real(kind=realkind) t(plond,plev)        ! temperature (k)
    real(kind=realkind) p(plond,plev)        ! pressure    (pa)
    real(kind=realkind) cwat(plond,plev)
    real(kind=realkind) cloud(plond,plev)
    !     

    !     
    real(kind=realkind) tsp(plond,plev)      ! saturation temp (k)
    real(kind=realkind) qsp(plond,plev)      ! saturation mixing ratio (kg/kg)
    !     
    !     
    integer i                 ! work variable
    integer k                 ! work variable
    integer iter              ! work variable
    integer l                 ! work variable
    !     

    !     
    real(kind=realkind) omeps                ! 1 minus epsilon
    real(kind=realkind) trinv                ! work variable
    real(kind=realkind) es                   ! sat. vapor pressure
    real(kind=realkind) desdtm               ! change in sat vap pressure wrt temperature
    real(kind=realkind) dqsdt                ! change in sat spec. hum. wrt temperature
    real(kind=realkind) dgdt                 ! work variable
    real(kind=realkind) g                    ! work variable
    real(kind=realkind) hlatsb               ! (sublimation)
    real(kind=realkind) hlatvp               ! (vaporization)
    real(kind=realkind) tterm                ! work var.
    real(kind=realkind) qs                   ! spec. hum. of water vapor
    real(kind=realkind) tc                   ! crit temp of transition to ice
    real(kind=realkind) t1,q1, dt, dq
    real(kind=realkind) dtm, dqm
    !     real qvd, r1, a1, tmp
    real(kind=realkind) qvd, a1, tmp
    !     real const
    !     real rair
    real(kind=realkind) ttrice
    real(kind=realkind) r1b, c1, c2, c3
    real(kind=realkind) cp                   ! heat capacity for dry air
    real(kind=realkind) hlatv                ! latent heat of vaporization
    real(kind=realkind) hlatf                ! latent heat of freezing
    real(kind=realkind) epsqs                ! epsilon
    real(kind=realkind) rgasv                ! gas constant for water vapour
    !     vo replacing the subroutine gestbl.f
    real(kind=realkind) pcf(1:5)             ! diff. sat. vapor press. over water and ice
    !     

    logical lflg              ! work variable
    !     

    !     
    real(kind=realkind) weight(plond)        ! work variable
    !     

    !     
    real(kind=realkind) hltalt(plond,plev)   ! lat. heat. of vap.
    logical lesat
    !     
    !     -------------------------------------------------------------


    !     computations
    !     
    tsp = -666.0_realkind
    epsqs = epsilo
    hlatv = latvap
    hlatf = latice
    !     gj
    !     gj   if we want water phase to be present below -20c we
    !     gj   need to modify other temp dependent variables in pcond
    !     gj   and findmcnew also to be consistent with ttrice=40 here
    !     gj   at present -20c is consistent throughout code.
    !     gj
    !     gjorig      ttrice = 20.
    ttrice=40._realkind
    cp = cpair
    omeps = 1.0_realkind - epsqs
    trinv = 1.0_realkind/ttrice
    a1 = 7.5_realkind*log(10._realkind)
    !     rair =  287.04
    rgasv = 461.50_realkind
    c3 = rair*a1/cp
    !     vo +++
    !     replacing the subroutine gestbl.f
    !     set coefficients for polynomial approximation of
    !     difference between saturation vapor press over water and saturation
    !     pressure over ice for -ttrice < t < 0 (degrees c). note: polynomial
    !     is valid in the range -40 < t < 0 (degrees c).
    !     
    !     --- degree 5 approximation ---
    !     
    pcf(1) =  5.04469588506e-01_realkind
    pcf(2) = -5.47288442819e+00_realkind
    pcf(3) = -3.67471858735e-01_realkind
    pcf(4) = -8.95963532403e-03_realkind
    pcf(5) = -7.78053686625e-05_realkind
    !     vo ---
    !     number of times to iterate the calculation
    iter = 6

    do k = 1,plev
       !     first guess on the saturation point
       do i = kstart, kstop

          !     gj080304
          if(lesat)then
             if(cwat(i,k)<=0._realkind)then
                es = esatw(t(i,k))
             else
                es = (cloud(i,k)*esat(t(i,k)))+&
                     ((1._realkind-cloud(i,k))*esatw(t(i,k)))
             endif
          else
             es = esat(t(i,k))
          endif
          es=min(es,p(i,k))
          qs = epsqs*es/(p(i,k) - omeps*es)
          !     gj080304
          !     saturation specific humidity
          !     
          !     "generalized" analytic expression for t derivative of es
          !     accurate to within 1 percent for 173.16 < t < 373.16

          !     weighting of hlat accounts for transition from water to ice
          !     polynomial expression approximates difference between es over
          !     water and es over ice from 0 to -ttrice (c) (min of ttrice is
          !     -40): required for accurate estimate of es derivative in transition 
          !     range from ice to water also accounting for change of hlatv with t 
          !     above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw

          tc     = t(i,k) - tmelt
          lflg   = (tc>=-ttrice .and. tc<0.0_realkind)
          weight(i) = min(-tc*trinv,1.0_realkind)
          hlatsb = hlatv + weight(i)*hlatf
          hlatvp = hlatv - 2369.0_realkind*tc
          if (t(i,k)<tmelt) then
             hltalt(i,k) = hlatsb
          else
             hltalt(i,k) = hlatvp
          endif

          !     e1 = hltalt(i,k)*q(i,k) + cp*t(i,k)
          tmp =  q(i,k) - qs
          c1 = hltalt(i,k)*c3
          c2 = (t(i,k) + 36._realkind)**2.0_realkind
          r1b    = c2/(c2 + c1*qs)
          qvd   = r1b*tmp
          tsp(i,k) = t(i,k) + ((hltalt(i,k)/cp)*qvd)
          !     gj080304
          if(lesat)then
             if(cwat(i,k)<=0._realkind)then
                es = esatw(tsp(i,k))
             else
                es = (cloud(i,k)*esat(tsp(i,k)))+&
                     ((1._realkind-cloud(i,k))*esatw(tsp(i,k)))
             endif
          else
             es = esat(tsp(i,k))
          endif
          es=min(es,p(i,k))
          qsp(i,k) = epsqs*es/(p(i,k) - omeps*es)
       enddo

       !     now iterate on first guess
       do i = kstart,kstop
          dt=0._realkind
          dq=0._realkind
          do l = 1, iter
             if ( qsp(i,k) < 1.0_realkind ) then
                !     gj080304
                if(lesat)then
                   if(cwat(i,k)<=0._realkind)then
                      es = esatw(tsp(i,k))
                   else
                      es = (cloud(i,k)*esat(tsp(i,k)))+ &
                           ((1._realkind-cloud(i,k))*esatw(tsp(i,k)))
                   endif
                   !     gj080304 
                else
                   es = esat(tsp(i,k))
                endif
                es=min(es,p(i,k))
                qs = epsqs*es/(p(i,k) - omeps*es)
                !     saturation specific humidity
                !     
                !     "generalized" analytic expression for t derivative of es
                !     accurate to within 1 percent for 173.16 < t < 373.16

                !     weighting of hlat accounts for transition from water to ice
                !     polynomial expression approximates difference between es over
                !     water and es over ice from 0 to -ttrice (c) (min of ttrice is
                !     -40): required for accurate estimate of es derivative in transition 
                !     range from ice to water also accounting for change of hlatv with t 
                !     above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw

                tc     = tsp(i,k) - tmelt
                lflg   = (tc>=-ttrice .and. tc<0.0_realkind)
                weight(i) = min(-tc*trinv,1.0_realkind)
                hlatsb = hlatv + weight(i)*hlatf
                hlatvp = hlatv - 2369.0_realkind*tc
                if (tsp(i,k)<tmelt) then
                   hltalt(i,k) = hlatsb
                else
                   hltalt(i,k) = hlatvp
                endif
                if (lflg) then
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)  &
                        + tc*(pcf(4) + tc*pcf(5))))
                else
                   tterm = 0.0_realkind
                endif


                !     the following check is to avoid the generation of negative
                !     values that can occur in the upper stratosphere and mesosphere
                !     kw070116 check removed, unnecessary if es<=p(i,k)

                desdtm = hltalt(i,k)*es &
                     /(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs) &
                     /(p(i,k) - omeps*es)*desdtm
                g = cp*(t(i,k)-tsp(i,k)) + hltalt(i,k)*q(i,k) &
                     - hltalt(i,k)*qsp(i,k)
                dgdt = -(cp + hltalt(i,k)*dqsdt)
                !     e2 = hltalt(i,k)*qsp(i,k) + cp*tsp(i,k)
                t1 = tsp(i,k) - g/dgdt
                dt = abs(t1 - tsp(i,k))/t1
                tsp(i,k) = t1
                !     gj080304
                if(lesat)then
                   if(cwat(i,k)<=0._realkind)then
                      es = esatw(tsp(i,k))
                   else
                      es = (cloud(i,k)*esat(tsp(i,k)))+ &
                           ((1._realkind-cloud(i,k))*esatw(tsp(i,k)))
                   endif
                   !     gj080304 
                else
                   es = esat(tsp(i,k))
                endif
                es=min(es,p(i,k))
                q1 = epsqs*es/(p(i,k) - omeps*es)
                dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_realkind)
                qsp(i,k) = q1

                !     write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq
             endif

             if (dt<1.e-4_realkind.and.dq<1.e-4_realkind) go to 10

          enddo               ! do l = 1,iter
10        continue
       enddo                  ! do i = 1,plon

       !     if (dtm>1.e-4.or.dqm>1.e-4) then
       !     write (6,*) ' findsp not converging ', dtm, dqm, k
       !     call stop_program('')
       !     endif

    enddo                     ! level loop (k=1,plev)

    return
  end subroutine findsp


  subroutine vqsatd(plond,kstart,kstop,lesat, &
       t    ,p  ,cwat   ,cloud  ,    &
       es      ,qs      ,gam)
    !-----------------------------------------------------------------------
    !
    ! utility procedure to look up and return saturation vapor pressure from 
    ! precomputed table, calculate and return saturation specific humidity 
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dt).  the same
    ! function as qsatd, but operates on vectors of temperature and pressure
    !
    !----------------------------code history-------------------------------
    !
    ! original version:  j. hack
    ! standardized:      j. rosinski, june 1992
    !                    t. acker, march 1996
    ! reviewed:          j. hack, august 1992
    !
    use confys
    use escom
    implicit none
    !------------------------------arguments--------------------------------
    !
    ! input arguments
    !
    integer plond       ! vector length
    integer kstart      ! start point of vector
    integer kstop       ! end point of vector
    real(kind=realkind) t(plond)       ! temperature
    real(kind=realkind) p(plond)       ! pressure
    real(kind=realkind) cwat(plond)
    real(kind=realkind) cloud(plond)
    ! 
    ! output arguments
    !
    real(kind=realkind) es(plond)   ! saturation vapor pressure
    real(kind=realkind) qs(plond)   ! saturation specific humidity
    real(kind=realkind) gam(plond)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------local variables------------------------------
    !
    logical lflg   ! true if in temperature transition region
    !
    integer i      ! index for vector calculations
    !
    real(kind=realkind) omeps     ! 1. - epsilo
    real(kind=realkind) trinv     ! reciprocal of ttrice (transition range)
    real(kind=realkind) tc        ! temperature (in degrees c)
    real(kind=realkind) weight    ! weight for es transition from water to ice
    real(kind=realkind) hltalt    ! appropriately modified hlat for t derivatives  
    !
    real(kind=realkind) hlatsb    ! hlat weighted in transition region
    real(kind=realkind) hlatvp    ! hlat modified for t changes above 273.16
    real(kind=realkind) tterm     ! account for d(es)/dt in transition region
    real(kind=realkind) desdtm     ! d(es)/dt
    real(kind=realkind) cp        ! heat capacity for dry air
    real(kind=realkind) hlatv     ! latent heat of vaporization
    real(kind=realkind) hlatf     ! latent heat of freezing
    real(kind=realkind) ttrice    ! temperature treshold for ice
    real(kind=realkind) rgasv     ! gas constant for water vapour
    !vo replacing the routine gestbl.f
    real(kind=realkind) pcf(1:5)  ! diff. sat. vapour press. over water and ice
    logical lesat


    hlatv = latvap
    hlatf = latice
    rgasv = 461.50_realkind
    omeps = 1.0_realkind - epsilo
    !gj
    !gj  see comments in findsp
    !gj
    !gjorig      ttrice = 20.
    ttrice = 40._realkind
    !gj      ttrice = 0.
    cp = cpair
    !vo +++
    ! replacing the routine gestbl.f
    ! set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees c). note: polynomial
    ! is valid in the range -40 < t < 0 (degrees c).
    !
    !                  --- degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_realkind
    pcf(2) = -5.47288442819e+00_realkind
    pcf(3) = -3.67471858735e-01_realkind
    pcf(4) = -8.95963532403e-03_realkind
    pcf(5) = -7.78053686625e-05_realkind
    !vo ---

    do i=kstart,kstop
       !gj080304
       if(lesat)then
          if(cwat(i)<=0._realkind)then
             es(i)=esatw(t(i))
          else
             es(i)=(cloud(i)*esat(t(i)))+((1._realkind-cloud(i))*esatw(t(i)))
          endif
          !gj080304
       else
          es(i) = esat(t(i))
       endif
       !
       ! saturation specific humidity
       !
       qs(i) = epsilo*es(i)/(p(i) - omeps*es(i))
       !
       ! the following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = min(1.0_realkind,qs(i))
       !
       if (qs(i) < 0.0_realkind) then
          qs(i) = 1.0_realkind
          es(i) = p(i)
       end if
    end do
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_realkind
    if(abs(ttrice)<1.e-14_realkind) goto 10
    trinv = 1.0_realkind/ttrice
    do i=kstart,kstop
       !
       ! weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (c) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition 
       ! range from ice to water also accounting for change of hlatv with t 
       ! above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc>=-ttrice .and. tc<0.0_realkind)
       weight = min(-tc*trinv,1.0_realkind)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_realkind*tc
       if (t(i)<tmelt) then
          hltalt = hlatsb
       else
          hltalt = hlatvp
       end if
       if (lflg) then
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       else
          tterm = 0.0_realkind
       end if
       desdtm  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdtm/(cp*es(i)*(p(i) - omeps*es(i)))
       if(abs(qs(i)-1.0_realkind)<1.e-14_realkind) gam(i) = 0.0_realkind
    end do
    return
    !
    ! no icephs or water to ice transition
    !
10  continue
    do i=kstart,kstop
       !
       ! account for change of hlatv with t above 273.16 where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_realkind*(t(i)-tmelt)
       hlatsb = hlatv + hlatf
       if (t(i)<tmelt) then
          hltalt = hlatsb
       else
          hltalt = hlatvp
       end if
       desdtm  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdtm/(cp*es(i)*(p(i) - omeps*es(i)))
       if(abs(qs(i)-1.0_realkind)<1.e-14_realkind) gam(i) = 0.0_realkind
    end do
    !
    return
    !
  end subroutine vqsatd

  subroutine findmcnew (plond,plev,kstart,kstop,           &
       k, precab, t, p, cwm, cldm, cldmax, &
       fice, coef,                         &
       frland,frice,                       &
       snowh,ps,gpot,pblh)

    ! calculate the conversion of condensate to precipitate

    !     written by phil rasch april 1997
    use confys
    use escom

    implicit none
    !

    !
    integer plond             ! horizontal vector length
    integer plev              ! vertical vector length
    integer kstart            ! start point of vector
    integer kstop             ! end point of vector
    integer k                 ! level index 
    !

    !
    real(kind=realkind) cldm(plond)          ! cloud fraction
    real(kind=realkind) cldmax(plond)        ! max cloud fraction above this level
    real(kind=realkind) cwm(plond)           ! condensate mixing ratio (kg/kg)
    real(kind=realkind) fice(plond,plev)          ! fraction of cwat that is ice
    real(kind=realkind) precab(plond)        ! rate of precipitation from above (kg / (m**2 * s))
    !

    !
    real(kind=realkind) t(plond,plev)        ! temperature       (k)
    real(kind=realkind) p(plond,plev)        ! pressure          (pa)
    !

    !
    real(kind=realkind) coef(plond)          ! conversion rate (1/s)
    real(kind=realkind) fwaut(plond)         ! relative importance of liquid autoconversion (a diagnostic)
    real(kind=realkind) fsaut(plond)         ! relative importance of ice autoconversion (a diagnostic)
    real(kind=realkind) fracw(plond)         ! relative  importance of rain accreting liquid (a diagnostic)
    real(kind=realkind) fsacw(plond)         ! relative  importance of snow accreting liquid (a diagnostic)
    real(kind=realkind) fsaci(plond)         ! relative  importance of snow accreting ice (a diagnostic)

    ! work variables

    integer i
    integer ii
    integer nlons

    real(kind=realkind) alpha                ! ratio of 3rd moment radius to 2nd
    !gjfire                         ! constant for autoconversion
    real(kind=realkind) capc
    !gjfire                         ! constant for autoconversion
    real(kind=realkind) capn                 ! local cloud particles / cm3
    real(kind=realkind) capnc                ! cold and oceanic cloud particles / cm3
    real(kind=realkind) capnw                ! warm continental cloud particles / cm3
    real(kind=realkind) ciaut                ! coefficient of autoconversion of ice (1/s)
    real(kind=realkind) ciautb               ! coefficient of autoconversion of ice (1/s
    real(kind=realkind) con1                 ! work constant
    real(kind=realkind) con2                 ! work constant
    real(kind=realkind) convfw               ! constant used for fall velocity calculation
    real(kind=realkind) cracw                ! constant used for rain accreting water
    real(kind=realkind) critpr               ! critical precip rate collection efficiency changes
    real(kind=realkind) csacx                ! constant used for snow accreting liquid or ice
    real(kind=realkind) dtice                ! interval for transition from liquid to ice
    !gjfire
    real(kind=realkind) effc                 ! collection efficiency
    !gjfire
    real(kind=realkind) icrit                ! threshold for autoconversion of ice
    real(kind=realkind) icritc               ! threshold for autoconversion of cold ice
    real(kind=realkind) icritw               ! threshold for autoconversion of warm ice
!!cgj230611
    real(kind=realkind) icritw1              ! threshold for autoconversion of warm ice
    real(kind=realkind) icritw2              ! threshold for autoconversion of warm ice
!!cgj230611
    real(kind=realkind) kconst               ! const for terminal velocity (stokes regime)
    real(kind=realkind) pracw                ! rate of rain accreting water
    real(kind=realkind) psaci                ! rate of collection of ice by snow (lin et al 1983)
    real(kind=realkind) psacw                ! rate of collection of liquid by snow (lin et al 1983)
    real(kind=realkind) psaut                ! rate of autoconversion of ice condensate
    real(kind=realkind) ptot                 ! total rate of conversion
    real(kind=realkind) pwaut                ! rate of autoconversion of liquid condensate
    real(kind=realkind) r3l                  ! volume radius
    real(kind=realkind) r3lcrit              ! critical radius at which autoconversion become efficient
    real(kind=realkind) rat1                 ! work constant
    real(kind=realkind) rat2                 ! work constant
    real(kind=realkind) rdtice               ! recipricol of dtice
    real(kind=realkind) rhocgs               ! density (cgs units)
    real(kind=realkind) snowfr               ! fraction of precipate existing as snow
    real(kind=realkind) vfallw               ! fall speed of precipitate as liquid
    real(kind=realkind) wp                   ! weight factor used in calculating pressure dep of autoconversion
    real(kind=realkind) wt                   ! fraction of ice
    real(kind=realkind) rhonot,t0,cldmin,small,c,d,esi,esw,nos,pir,prhonos, &
         thrpd,gam3pd,gam4pd,rhoi,rhosr,rhow,mcon05,mcon07,mcon08
    !gjtest210700
    real(kind=realkind) mcon01,mcon02,mcon03,mcon04,lamdas,csacw
    real(kind=realkind) snowmr(plond)
    real(kind=realkind) ps(plond)
    real(kind=realkind) pblh(plond),wzz
    real(kind=realkind) gpot(plond,plev)

    !gjtest210700

    !

    !
    integer ind(plond)
    !

    !
    real(kind=realkind) cldloc(plond)        ! non-zero amount of cloud
    real(kind=realkind) cldpr(plond)         ! assumed cloudy volume occupied by rain and cloud
    real(kind=realkind) icemr(plond)         ! in-cloud ice mixing ratio
    real(kind=realkind) liqmr(plond)         ! in-cloud liquid water mixing ratio
    real(kind=realkind) prlloc(plond)        ! local rain flux in mm/day
    real(kind=realkind) prscgs(plond)        ! local snow amount in cgs units
    real(kind=realkind) rainmr(plond)        ! in-cloud rain mixing ratio
    real(kind=realkind) rho(plond)           ! density (mks units)
    real(kind=realkind) totmr(plond)         ! in-cloud total condensate mixing ratio
    !

    !
    !      real hltalt(plond,plev)   ! lat. heat. of vap.
    !gj
    !gjfire
    real(kind=realkind) r3term,r3thresh
    real(kind=realkind) r3lc2
    !gjfire
    real(kind=realkind) frland(plond)
    real(kind=realkind) frice(plond)
    integer landm(plond)	!=0 if frland=0, otherwise=1
    real(kind=realkind) wpterm,wptt,ttbase,ttdiv,efact,efact1
    !gj160404
    real(kind=realkind) icritcl,icritch,wzkm,zwk1,zwrk,zwrk1,zpwr
    real(kind=realkind) zkk1,zkk2
    !gj080605
    real(kind=realkind) capnsi,zfrice,snowh(plond)



    !     giving value to constants
    !
    call inimc (rhonot, t0, cldmin, small,         &
         c, d, nos, pir, prhonos, thrpd,     &
         rhosr, gam3pd, gam4pd, rhow, rhoi,  &
         esi, esw,                           &
         mcon05, mcon07, mcon08,             &
         mcon01, mcon02, mcon03, mcon04)
    !gjtest210700
    !
    ! critical precip rate at which we assume the collector drops can change the
    ! drop size enough to enhance the auto-conversion process (mm/day)
    critpr=1.0_realkind
    !gj140800
    convfw = 1.94_realkind*2.13_realkind*sqrt(rhow*1000._realkind*9.81_realkind*2.7e-4_realkind)
    ! liquid microphysics
    !      cracw = 6                 ! beheng
    cracw = .884_realkind*sqrt(9.81_realkind/(rhow*1000._realkind*2.7e-4_realkind)) ! tripoli and cotton
    ciautb = 1.e-3_realkind	!lin83
!cgj230611    icritw = 1.e-3_realkind	!lin83
!!cgj050411    icritw = 4.e-4	
    icritw1 = 1.e-3_realkind
    icritw2 = 4.e-4_realkind
!cgjtest231011    icritw2 = 5.e-5_realkind
!!cgj230611
    icritc = 5.e-6_realkind
    icritcl=1.e-4_realkind
    icritch=5.e-6_realkind
!cgj030711
    zkk1=2.47_realkind
    zkk2=-1.79_realkind

    dtice = 20._realkind
    rdtice = 1._realkind/dtice

!    capnw = 400._realkind              ! warm continental cloud particles / cm3
    capnw = 250._realkind              ! warm continental cloud particles / cm3
!gultepetest    capnw = 200._realkind              ! warm continental cloud particles / cm3
    capnc =  80._realkind              ! cold and oceanic cloud particles / cm3
    capnc = 150._realkind              ! cold and oceanic cloud particles / cm3
    capnc = 100._realkind              ! cold and oceanic cloud particles / cm3
!gultepetest    capnc = 120._realkind              ! cold and oceanic cloud particles / cm3
    capnsi = 75._realkind              ! cold and oceanic cloud particles / cm3
    !
    ! for europe
    r3lcrit = 8.e-6_realkind           ! 9u crit radius where liq conversion begins
!cgj020711
    r3lcrit = 10.e-6_realkind           ! 9u crit radius where liq conversion begins
    ! for africa: test only: 05/10/11
!    r3lcrit = 7.5e-6           ! 10u crit radius where liq conversion begins

    !gj140800
    !gjwptt
    ttbase=243.15_realkind
    ttdiv=30._realkind
    !gjwptt


    ! find all the points where we need to do the microphysics
    ! and set the output variables to zero
    nlons = 0
    do i = kstart,kstop
       !gj
       if(frland(i)>0._realkind)then
          landm(i)=1
       else
          landm(i)=0
       endif
       !          write(*,*)'in findmcnew, landm=',landm(i)
       !gj
       coef(i) = 0._realkind
       fwaut(i) = 0._realkind
       fsaut(i) = 0._realkind
       fracw(i) = 0._realkind
       fsacw(i) = 0._realkind
       fsaci(i) = 0._realkind
       if (cwm(i)>1.e-20_realkind) then
          nlons = nlons + 1
          ind(nlons) = i
       endif
    end do

    !dir$ ivdep
    do ii = 1,nlons

       i = ind(ii)

       !        the local cloudiness at this level
       cldloc(i) = max(cldmin,cldm(i))

       !        a weighted mean between max cloudiness above, and this layer
       cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5_realkind)

       !        decompose the suspended condensate into 
       !        an incloud liquid and ice phase component
       totmr(i) = cwm(i)/cldloc(i)
       icemr(i) = totmr(i)*fice(i,k)
       liqmr(i) = totmr(i)*(1._realkind-fice(i,k))

       !        density
       rho(i) = p(i,k)/(287._realkind*t(i,k))
       rhocgs = rho(i)*1.e-3_realkind     ! density in cgs units

       !        decompose the precipitate into a liquid and ice phase 
       if (t(i,k)>t0) then
          vfallw = convfw/sqrt(rho(i))
          rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
          snowfr = 0._realkind
       else
          snowfr = 1._realkind
          rainmr(i) = 0._realkind
       endif

       prscgs(i) = precab(i)/cldpr(i)*0.1_realkind*snowfr ! local snow amount in cgs units
       prlloc(i) = precab(i)*86400._realkind/cldpr(i)    ! local rain amount in mm/day

    end do
    !gjfire
    kconst = 1.18e6_realkind           ! const for terminal velocity

    !      effc = 1.                 ! autoconv collection efficiency following boucher 96
    !      effc = .55*0.05           ! autoconv collection efficiency following baker 93
    effc = 0.55_realkind                ! autoconv collection efficiency following tripoli and cotton
    !gjfire
    !       do i=kstart,kstop
    !        if( (prlloc(i)/critpr) <=0.001)then
    !         effc(i)=0.55*0.001
    !        else
    !         effc(i)=0.55*( (prlloc(i)/critpr)**2. )
    !        endif
    !        if(effc(i)>0.55)effc(i)=0.55
    !gjfire
    !
    alpha = 1.1_realkind**4.0_realkind
    capc = pir**(-.333_realkind)*kconst*effc & ! constant for autoconversion
           *(0.75_realkind)**(1.333_realkind)*alpha
    !      enddo


    con1 = 1._realkind/(1.333_realkind*pir)**0.333_realkind * 0.01_realkind ! meters

    !     calculate the conversion terms
    !dir$ ivdep
    do ii = 1,nlons

       i = ind(ii)

       rhocgs = rho(i)*1.e-3_realkind     ! density in cgs units

       !        exponential temperature factor
       !gjpsaut
       efact = exp(0.025_realkind*(t(i,k)-t0))
       !gjpsaut

       wzz=gpot(i,k)/gravit
       !gj160404 make icritc a function of height
       !gj160404
       !gj110605         wzkm=wzz/1000.
       !gj110605         zwk1=wzkm
       !gj110605         icritc=icritcl*exp(-zwk1)
       !gj110605         icritc=max(icritc,icritch)
       !gj160404
       !        some temperature dependent constants
       wt = min(1._realkind,max(0._realkind,(t0-t(i,k))*rdtice))

!!cgj230611 Make efact=efact below 0.75*ps and over land only and 1 elsewhere,efact slows down the autoconversion 
!!cgj230611 (psaut) and collection (psaci) rate of cloud ice to snow processes in the temp range 0C to -20C. N.B.
!!cgj230611 this essentially leads to all cloud ice to snow processes in the tropics being without the efact term
!!cgj230611 as they generally do not occur below 0.75*ps altitudes due to temperatures being too warm. This process
!!cgj230611 mimics what is done for liquid autoconversion, where we assume it occurs more slowly over lower levels
!!cgj230611 over land due to a higher number of CCN, and more rapidly elsewhere (free atmos over land and everywhere
!!cgj230611 over ocean and sea-ice. Modified icritw is basically doing a similar job, making the cloud ice threshold
!!cgj230611 in the 0C to -20C temp range, higher before autoconversion of ice can start over land below 0.75ps and 
!!cgj230611 lower elsewhere.
!!cgj230611
       zpwr = 0.1_realkind
       zwrk  = min(1._realkind,max(0._realkind, &
             (p(i,k)-0.8_realkind*ps(i))/(0.2_realkind*ps(i))))
       zwrk1 = zwrk**zpwr
       icritw = icritw1*zwrk1 + icritw2*(1._realkind-zwrk1)
!cgj       icritc = icritcl*zwrk1 + icritch*(1._realkind-zwrk1)
!cgj060711       if(efact1>0.)then
!cgj060711        efact=efact1*((1.-zwrk1)/efact1 + zwrk1)
!cgj060711       else
!cgj060711        efact=efact1
!cgj060711       endif
!!cgj230611 Force icritw=icritw2 through all the atmosphere over non-land points.
       icritw = icritw + (icritw2-icritw)*min(1._realkind,max &
                (0._realkind,1._realkind-real(landm(i),realkind)))
!!cgj230611
       icrit = icritc*wt + icritw*(1._realkind-wt)
!cgjlin83_240911
!
       !        linear weight factor in pressure (1 near sfc, 0 at .8 of sfc) 
       !gj210800         wp = min(1.,max(0.,(p(i,k)-0.8*p(i,plev))/(0.2*p(i,plev))))
       !gj061003         wp = min(1.,max(0.,(p(i,k)-0.7*ps(i))/(0.3*ps(i))))
       !gjwptt
       !gjwptt  also make the capn term a function of temperature use
       !gj      capn=150 at -30c and below & capn=capn as fn of p only
       !gj      above 0c, smooth transition in between
       !gjwptt
       !gj         wpterm=(t(i,k)-ttbase)/ttdiv 
       !gj         wptt =min(wp,max(0.,wp*wpterm))
       !gjwptt
       !gj080605         if(wzz<=pblh(i))then
       !gj080605          wp=1.
       !gj080605         else
       !gj080605          wp=min(1.,max(0.,(1.-(0.00033*(wzz-pblh(i))))))
       !gj080605         endif
       !        near land near sfc raise the number concentration
       !gj080605
       !gj     modify ccn to be lased on land/sea value with mods for snow
       !gj     covered land and sea-ice fraction.
       !gj
       !gj100605
       zfrice=frice(i)*(1._realkind-frland(i))
       capn = capnw+(capnc-capnw)*min(1._realkind,max(0._realkind, &
            1._realkind-(p(i,k)-0.8_realkind*ps(i))/(0.2_realkind*ps(i))))
       ! ramp between polluted air over land to clean over sea
       capn = capn + (capnc-capn)*min(1._realkind,max(0._realkind,1._realkind-real(landm(i),realkind)))
       !modify for snow depth over land
       !gj          capn = capn + (capnsi-capn)*min(1.,max(0.,snowh(i)))
       capn = capn + (capnsi-capn)*min(1._realkind,max(0._realkind,snowh(i)*10._realkind))
       ! ramp between resulting value and sea-ie value
       capn = capn + (capnsi-capn)*min(1._realkind,max(0._realkind,zfrice))
       !gj100605
       !gj080605
       !gj080605            capn =  landm(i)*(capnw*wp + capnc*(1-wp))
       !gj080605     $          +(1.-landm(i))*capnc
       !gj            capn =  landm(i)*(capnw*wptt + capnc*(1-wptt))
       !gj     $          +(1.-landm(i))*capnc
       !gjwptt
       !gj131299            capn =  capnw*wp + capnc*(1-wp)

       !        useful terms in following calculations
       rat1 = rhocgs/rhow
       rat2 = liqmr(i)/capn
       con2 = (rat1*rat2)**0.333_realkind 

       !        volume radius
       !        r3l = (rhocgs*liqmr(i)/(1.333*pir*capn*rhow))**0.333 * 0.01 ! meters
       r3l = con1*con2 

       !        critical threshold for autoconversion if modified for mixed phase
       !        clouds to mimic a bergeron findeisen process
       r3lc2 = r3lcrit*(1._realkind-0.5_realkind*fice(i,k)*(1._realkind-fice(i,k)))
       r3lc2 = max(r3lc2,1.e-6_realkind)

       !        autoconversion of liquid
       !        cwaut = 2.e-4
       !        cwaut = 1.e-3
       !        lcrit = 2.e-4
       !        lcrit = 5.e-4
       !        pwaut = max(0.,liqmr(i)-lcrit)*cwaut

       ! pwaut is following tripoli and cotton (and many others)
       ! we reduce the autoconversion below critpr, because these are regions where
       ! the drop size distribution is likely to imply much smaller collector drops than
       ! those relevant for a cloud distribution corresponding to the value of effc = 0.55 
       ! suggested by cotton (see austin 1995 jas, baker 1993)

       !gjfire
       r3term=0._realkind
       r3thresh=max(0._realkind,r3l-r3lcrit)
       !gj         r3thresh=max(0.,r3l-r3lc2)
       if(r3thresh>0._realkind)r3term=1._realkind
       !gjfire
       pwaut = capc*liqmr(i)**2.0_realkind*rat1                   &
            *con2                                    &
            *max(0.00000001_realkind,sign(1._realkind,r3l-r3lcrit))    &
            *max(0.05_realkind,min(1._realkind,((prlloc(i)/critpr))))
       !gj
       !gjfire
       !gjrca2         pwaut = capc(i)*liqmr(i)**2*rat1
       !gjrca2     $           *con2
       !gjrca2     $           *r3term
       !gjnewpwaut
       !!pwaut=1350.*(liqmr(i)**zkk1)*(capn**zkk2)
       !gjnewpwaut
       !        autoconversion of ice
       !gjlin83
       !gj110605        ciaut = ciautb*efact
       !gj110605
!!cgj050411       ciaut = ciautb
!!cgj050411 reset ciaut to rca3 equivalent ciaut with efact included will be
!!cgj050411 than when not included in the temp range -20C to 0C
!cgj300911       ciaut = ciautb*efact
       ciaut = ciautb
       !        psaut = capc*totmr(i)**2*rhocgs/rhoi
       !     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
       !   
!!cgj050411 effect of using rca3 icrit will be that ice cloud water conversion
!!cgj050411 will start (later) requiring a higher value of icemr than in rca4
!!cgj050411 the effect of ciaut=ciautb*efact will furthermore mean that when 
!!cgj050411 does start it will also occur more slowly. All this is for the
!!cgj050411 temperature range -20C to 0C.
       psaut = max(0._realkind,icemr(i)-icrit)*ciaut ! autoconversion of ice condensate


       !        collection of liquid by rain 
       !        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
       pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)

       ! the following lines calculate the slope parameter and snow mixing ratio
       ! from the precip rate using the equations found in lin et al 83
       ! in the most natural form, but it is expensive, so after some tedious
       ! algebraic manipulation you can use the cheaper form found below
       !            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
       !     $               *0.01   ! convert from cm/s to m/s
       !            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
       !gjtest210700
       !            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
       !            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
       !            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd) 
       !gjtest210700


       !        coefficient for collection by snow independent of phase
       csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05

       !        collection of liquid by snow (lin et al 1983)
       psacw = csacx*liqmr(i)*esw
       !         psacw = csacw*liqmr(i)*esw

       !        collection of ice by snow (lin et al 1983)
       !ps110706 psaci = csacx*icemr(i)*esi
       !cgj only allow efact effect on ice autoconversion
       psaci = csacx*icemr(i)*esi
!cgj30091       psaci = csacx*icemr(i)*(min(1._realkind,efact))	!rca3 version of code

       !        total conversion of condensate to precipitate 
       ptot = pwaut + psaut + pracw + psacw + psaci

       !        the recipricol of cloud water amnt (or zero if no cloud water)
       !         rcwm =  totmr(i)/(max(totmr(i),small)**2)

       !        turn the tendency back into a loss rate (1/seconds)
       if (totmr(i)>0._realkind) then
          coef(i) = ptot/totmr(i)
       else
          coef(i) = 0._realkind
       endif

       !gj         if(ptot>small)then
       !gj        fwaut(i) = pwaut/ptot
       !gj          fsaut(i) = psaut/ptot
       !gj          fracw(i) = pracw/ptot
       !gj          fsacw(i) = psacw/ptot
       !gj          fsaci(i) = psaci/ptot
       !gj         endif
       fwaut(i) = pwaut/max(ptot,small)
       fsaut(i) = psaut/max(ptot,small)
       fracw(i) = pracw/max(ptot,small)
       fsacw(i) = psacw/max(ptot,small)
       fsaci(i) = psaci/max(ptot,small)
       !gj

    end do

    return
  end subroutine findmcnew

  subroutine inimc (rhonot, t0, cldmin, small, &
       c, d, nos, pir, prhonos, thrpd,          &
       rhosr, gam3pd, gam4pd, rhow, rhoi,       &
       esi, esw,                                &
       mcon05, mcon07, mcon08,                  &
       mcon01, mcon02, mcon03, mcon04)

    use confys, only:tmelt
    implicit none
    !     initialize the common block variables for the prognostic condensate parameterization

    !     phil rasch april 1997      

    real(kind=realkind) rhonot,t0,cldmin,small,c,d,esi,esw,nos,pir,prhonos,&
         thrpd,gam3pd,gam4pd,rhoi,rhosr,rhow,mcon01,mcon02,  &
         mcon03,mcon04,mcon05,mcon06,mcon07,mcon08,rhos
    !     gj
    real(kind=realkind) mcon02_cj,mcon07_cj


    rhonot = 1.275e-3_realkind         ! air density at surface (gm/cm3)

    rhos = 0.1_realkind                ! assumed snow density (gm/cm3)
    rhow = 1._realkind                 ! water density
    rhoi = 1._realkind                 ! ice density

    esi = 1.0_realkind                 ! collection efficient for ice by snow
    !     gj      esi = 0.5                 ! collection efficient for ice by snow
    esw = 0.1_realkind                 ! collection efficient for water by snow

    t0 = tmelt                ! approximate freezing temp

    cldmin = 0.02_realkind             ! assumed minimum cloud amount 

    small = 1.e-22_realkind            ! a small number compared to unity

    c = 152.93_realkind                ! constant for graupel like snow cm**(1-d)/s

    d = 0.25_realkind                  ! constant for graupel like snow

    nos = 3.e-2_realkind               ! particles snow / cm**4

    pir = 4._realkind*atan(1.0_realkind)         

    prhonos = pir*rhos*nos

    thrpd = 3._realkind + d

    if(abs(d-0.25_realkind)<1.e-14_realkind) then
       gam3pd = 2.549256966718531_realkind ! only right for d = 0.25
       gam4pd = 8.285085141835282_realkind
    else
       write (6,*) ' can only use d ne 0.25 on a cray '
       call stop_program(' can only use d ne 0.25 on a cray ')
    endif

    mcon01 = pir*nos*c*gam3pd/4._realkind
    mcon02 = 1._realkind/(c*gam4pd*sqrt(rhonot)/(6._realkind*prhonos**(d/4._realkind)))
    !     mcon02_cj = ( 6*( prhonos**(d+4.) ) )/(c*gam4pd*sqrt(rhonot))
    mcon03 = -(0.5_realkind+d/4._realkind)
    mcon04 = 4._realkind/(4._realkind+d)
    mcon05 = (3._realkind+d)/(4._realkind+d)
    mcon06 = (3._realkind+d)/4._realkind
    mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
    mcon08 = -0.5_realkind/(4._realkind+d)
    !     mcon07_cj = mcon01*sqrt(rhonot)*mcon02_cj**mcon05/
    !     +            ((rhos*nos)**mcon06)

    !     write (6,*) ' inimc complete ' 

    return
  end subroutine inimc

end module pcondmod
