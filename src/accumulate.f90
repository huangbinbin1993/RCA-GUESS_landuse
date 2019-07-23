module accumulate
  use timetype
  use calendar
  use comhkp,only:ndtime
  use postprocess
  use decomp
  implicit none
  private
  public accumulate1,acc_surf,acc4mean_surf,accmean1,accscratch4mean_surf,&
       accmean_surf,accscratch_surf,accscratch1,oasacc,oasaccmean,oasaccscratch
  public scratch_post_process_surface,scratch_post_process
  public accumulate_post_process_surface,accumulate_post_process

contains

!  subroutine mean_and_accumulation(svar_surf,account_surf,mean_surf,acum_surf)
!
!    implicit none
!    real(kind=realkind),dimension(:,:,:),intent(in)   :: svar_surf
!    integer,intent(in)                 :: account_surf
!    real(kind=realkind),dimension(:,:,:),intent(inout):: mean_surf,acum_surf,extrem_surf
!
!    integer:: n
!
!    do n=1,length(account_surf)
!       mean_surf(:,:,n) = 1./(account_surf(n)+1)*svar_surf(:,:,n) + account_surf(n)/(account_surf(n)+1)*mean_surf(:,:,n)
!       account_surf(n) = account_surf(n)+1
!    enddo
!
!    acum_surf = acum_surf + svar_surf 
!
!  end subroutine scratch_mean_and_accumulation

  subroutine accumulate_post_process(nstep,klon,klat,accsswr,accslwr,acctswr,acctlwr,accslwdn,accsswdn,acctswdn,accq2d,accci2d,&
       acccw2d,accqu2d,accqv2d,accvimfc2d,acctcov,acclcov,accmcov,acchcov)
    implicit none
    integer,intent(in)::nstep,klon,klat
    real(kind=realkind),dimension(:,:)::acccw2d,accslwr,accsswr,acctlwr
    real(kind=realkind),dimension(:,:)::acctswr,accsswdn,accslwdn,acctswdn,accq2d,accci2d,accqu2d,accqv2d,accvimfc2d, &
                                        acctcov,acclcov,accmcov,acchcov 
    integer::jy,i,ii

    integer::dtime

    dtime = Nseconds(ndtime) !(ndtime%hour*60 + ndtime%min)*60 + ndtime%sec
    if( nstep > 0 )then
       do jy = 1,nppstr
          do i=1,files(jy)%nsl
             if( files(jy)%meantime_3006(i) /= 0 )then
                if( mod(nstep,files(jy)%meantime_3006(i)/dtime) == 0)then
                   ii=files(jy)%iwmosl(i)
                   if(ii == 111) then
                      call accmean1(klon,klat,accsswr,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 112) then
                      call accmean1(klon,klat,accslwr,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 113) then
                      call accmean1(klon,klat,acctswr,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 114) then
                      call accmean1(klon,klat,acctlwr,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 115) then
                      call accmean1(klon,klat,accslwdn,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 116) then
                      call accmean1(klon,klat,accsswdn,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 117) then
                      call accmean1(klon,klat,acctswdn,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 54) then
                      call accmean1(klon,klat,accq2d,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 58) then
                      call accmean1(klon,klat,accci2d,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 76) then
                      call accmean1(klon,klat,acccw2d,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 151) then
                      call accmean1(klon,klat,accqu2d,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 152) then
                      call accmean1(klon,klat,accqv2d,files(jy)%meantime_3006(i)/dtime)
                   elseif(ii == 153) then
                      call accmean1(klon,klat,accvimfc2d,files(jy)%meantime_3006(i)/dtime)
!CUW110301beg
                    elseif (ii .eq. 71) then
                       call accmean1(klon,klat,acctcov,files(jy)%meantime_3006(i)/dtime)
                    elseif (ii .eq. 73) then
                       call accmean1(klon,klat,acclcov,files(jy)%meantime_3006(i)/dtime)
                    elseif (ii .eq. 74) then
                       call accmean1(klon,klat,accmcov,files(jy)%meantime_3006(i)/dtime)
                    elseif (ii .eq. 75) then
                       call accmean1(klon,klat,acchcov,files(jy)%meantime_3006(i)/dtime)
!CUW110301end
                   endif
                endif
             endif
          enddo
       enddo               ! jy
    endif                  !   nstep>0

  end subroutine accumulate_post_process

  
  subroutine accumulate_post_process_surface(nstep,klon,klat,msvars,account_surf,acum4mean_surf,extrem_surf,utot10ms)
    implicit none
    integer,intent(in)::nstep,klon,klat,msvars
    integer::account_surf(:)
    real(kind=realkind),dimension(:,:,:)::acum4mean_surf
    real(kind=realkind),dimension(:,:,:)::extrem_surf
    real(kind=realkind),dimension(:,:):: utot10ms
    integer::i,jy,ii,jk,jx
    integer::dtime

    dtime =(ndtime%hour*60 + ndtime%min)*60 + ndtime%sec

    if( nstep > 0 )then
       do jy = 1,nppstr
          do i=1,files(jy)%nsl
             if( files(jy)%meantime_surf(i) /= 0 )then
                if( mod(nstep,files(jy)%meantime_surf(i)/dtime) == 0)then
                   ii=files(jy)%alevsl(i)
                   call accmean_surf(klon,klat,account_surf,msvars,ii,acum4mean_surf) !surface
                endif
             endif
          enddo

          !     special case for extrem(18) that is not extreme but mean
          do i=1,files(jy)%nsl
             ii=files(jy)%alevsl(i)
             if(ii == 18) then
                if( files(jy)%scratchtime_extrem(i) /= 0 )then
                   if(mod(nstep,files(jy)%scratchtime_extrem(i)/dtime)==0)then
                      do jk=1,klat
                         do jx=1,klon
                            if( utot10ms(jx,jk) > 0.0_realkind ) then
                               extrem_surf(jx,jk,ii)=utot10ms(jx,jk)/real(files(jy)%scratchtime_extrem(i),realkind)*&
                                    real(dtime,realkind)
                            else
                               extrem_surf(jx,jk,ii) = -1.0_realkind
                            endif
                         enddo
                      enddo
                   endif
                endif
             endif
          enddo
       enddo
    endif
  end subroutine accumulate_post_process_surface



  subroutine scratch_post_process_surface(nstep,klon,klat,msvars,account_surf,acum4mean_surf,acum_surf, &
       extrem_surf,utot10ms )
    implicit none
    integer,intent(in)::nstep,klon,klat,msvars
    integer::account_surf(:)
    real(kind=realkind),dimension(:,:,:)::acum4mean_surf,acum_surf
    real(kind=realkind),dimension(:,:,:)::extrem_surf
    real(kind=realkind),dimension(:,:)::utot10ms
    integer::jy,i,ii

    integer::dtime

    dtime =(ndtime%hour*60 + ndtime%min)*60 + ndtime%sec

    do jy = 1,nppstr

       do i=1,files(jy)%nsl
          if( files(jy)%meantime_surf(i) /= 0 )then
             if( mod(nstep,files(jy)%meantime_surf(i)/dtime)== 0)then
                ii=files(jy)%alevsl(i)
                call accscratch4mean_surf(klon,klat,account_surf,msvars,ii,acum4mean_surf)
             endif
          endif
       enddo

       do i=1,files(jy)%nsl
          if( files(jy)%scratchtime_acum(i) /= 0 )then
             if( mod(nstep,files(jy)%scratchtime_acum(i)/dtime)== 0)then
                ii=files(jy)%alevsl(i)
                call accscratch_surf(klon,klat, msvars,ii,acum_surf)
             endif
          endif
       enddo

       do i=1,files(jy)%nsl
          if( files(jy)%scratchtime_extrem(i) /= 0 )then
             if(mod(nstep,files(jy)%scratchtime_extrem(i)/dtime) ==0)then
                ii=files(jy)%alevsl(i)
                if(ii==1.or.ii==3.or.ii==5.or.ii==7.or. &
                     ii==9.or.ii==11.or.                      &
                     ii==13.or.ii==14.or.ii==15.or.         &
                     ii==16.or.ii==19.or.                     &
                     ii==21.or.ii==22.or.ii==23             &
                     .or.ii==24.or.                             &
                     ii==26.or.ii==27) then
                   !call scratchhi(klon,klat,extrem_surf(:,:,ii))
                   extrem_surf(:,:,ii) = -999.0_realkind
                elseif(ii==2.or.ii==4.or.ii==6.or.    &
                     ii==8.or.                           &
                     ii==10.or.ii==12.or.              &
                     ii==17.or.ii==20.or.ii==25) then
!                   call scratchlo(klon,klat,extrem_surf(:,:,ii))
                   extrem_surf(:,:,ii) = 999.0_realkind
                elseif(ii==18) then
                   utot10ms = 0.0_realkind
                endif
             endif
          endif
       enddo
    enddo !do nppstr jy
  end subroutine scratch_post_process_surface



  subroutine scratch_post_process(nstep,klon,klat,accsswr,accslwr,acctswr,acctlwr,accslwdn,accsswdn,acctswdn,&
       accq2d,accci2d,acccw2d,accqu2d,accqv2d,accvimfc2d,accsunny,accrunoff,accrunoffopl,accrunofffor,accrunofflake,&
       accsnow,accsnol,accsnoc,accprl,lastpr,accprc,tsmax,t2max,tsmin,t2min,uv10max,maxdayprecip,&
       acctcov,acclcov,accmcov,acchcov)
    implicit none
    integer,intent(in)::nstep,klon,klat
    real(kind=realkind),dimension(:,:)::accsswr,accslwr,acctswr,acctlwr,accslwdn,accsswdn,acctswdn,accq2d,accci2d,acccw2d,accqu2d
    real(kind=realkind),dimension(:,:)::accqv2d,accvimfc2d,accsunny,accrunoff,accrunoffopl,accrunofffor,accrunofflake
    real(kind=realkind),dimension(:,:)::accsnow,accsnol,accsnoc,accprl,lastpr,accprc,tsmax,t2max,tsmin,t2min,uv10max
    real(kind=realkind),dimension(:,:)::acctcov,acclcov,accmcov,acchcov
    real(kind=realkind),dimension(:,:)::maxdayprecip
    integer::jy,i,ii

    integer::dtime

    dtime =(ndtime%hour*60 + ndtime%min)*60 + ndtime%sec

    do jy = 1,nppstr
       do i=1,files(jy)%nsl
          if( files(jy)%meantime_3006(i) /= 0 )then
             if( mod(nstep,files(jy)%meantime_3006(i)/dtime)== 0)then
                ii=files(jy)%iwmosl(i)
                if(ii == 111) then
                   call accscratch1(klon,klat,accsswr)
                elseif(ii == 112) then
                   call accscratch1(klon,klat,accslwr)
                elseif(ii == 113) then
                   call accscratch1(klon,klat,acctswr)
                elseif(ii == 114) then
                   call accscratch1(klon,klat,acctlwr)
                elseif(ii == 115) then
                   call accscratch1(klon,klat,accslwdn)
                elseif(ii == 116) then
                   call accscratch1(klon,klat,accsswdn)
                elseif(ii == 117) then
                   call accscratch1(klon,klat,acctswdn)
                elseif(ii == 54) then
                   call accscratch1(klon,klat,accq2d)
                elseif(ii == 58) then
                   call accscratch1(klon,klat,accci2d)
                elseif(ii == 76) then
                   call accscratch1(klon,klat,acccw2d)
                elseif(ii == 151) then
                   call accscratch1(klon,klat,accqu2d)
                elseif(ii == 152) then
                   call accscratch1(klon,klat,accqv2d)
                elseif(ii == 153) then
                   call accscratch1(klon,klat,accvimfc2d)
!CUW110301beg
                elseif (ii .eq. 71) then
                   call accscratch1(klon,klat,acctcov)
                elseif (ii .eq. 73) then
                   call accscratch1(klon,klat,acclcov)
                elseif (ii .eq. 74) then
                   call accscratch1(klon,klat,accmcov)
                elseif (ii .eq. 75) then
                   call accscratch1(klon,klat,acchcov)
!CUW110301end
                endif
             endif
          endif
       enddo

       do i=1,files(jy)%nsl
          if(files(jy)%scratchtime_4006(i) /= 0 )then
             if(mod(nstep,files(jy)%scratchtime_4006(i)/dtime)==0)then
                ii=files(jy)%iwmosl(i)
                if(ii == 240) then
                   call accscratch1(klon,klat,accsunny)
                elseif(ii == 90) then
                   call accscratch1(klon,klat,accrunoff)
                   call accscratch1(klon,klat,accrunoffopl)
                   call accscratch1(klon,klat,accrunofffor)
                   call accscratch1(klon,klat,accrunofflake)
                elseif(ii == 65) then
                   call accscratch1(klon,klat,accsnow)
                elseif(ii == 79) then
                   call accscratch1(klon,klat,accsnol)
                elseif(ii == 78) then
                   call accscratch1(klon,klat,accsnoc)
                elseif(ii == 62) then
                   call accscratch1(klon,klat,accprl)
                   call accscratch1(klon,klat,lastpr)
                elseif(ii == 63) then
                   call accscratch1(klon,klat,accprc)
                   call accscratch1(klon,klat,lastpr)
                endif
             endif
          endif
       enddo

       do i=1,files(jy)%nsl
          if( files(jy)%scratchtime_hilo(i) /= 0 )then
             if(mod(nstep,files(jy)%scratchtime_hilo(i)/dtime)==0)then
                ii=files(jy)%iwmosl(i)
                if(ii == 15) then
                   if(files(jy)%alevsl(i) == 0) then
                      tsmax = -999.0_realkind
                      !call scratchhi(klon,klat,tsmax)
                   elseif(files(jy)%alevsl(i) == 2) then
                      !call scratchhi(klon,klat,t2max)
                      t2max = -999.0_realkind
                   endif
                elseif(ii == 16) then
                   if(files(jy)%alevsl(i) == 0) then
                      !call scratchlo(klon,klat,tsmin)
                      tsmin = 999.0_realkind
                   elseif(files(jy)%alevsl(i) == 2) then
                      !call scratchlo(klon,klat,t2min)
                      t2min = 999.0_realkind
                   endif
                elseif(ii == 32) then
                   !call scratchhi(klon,klat,uv10max)
                   uv10max = -999.0_realkind
                elseif(ii == 64) then
                   maxdayprecip = 0._realkind
                endif
             endif
          endif
       enddo
    enddo
  end subroutine scratch_post_process


  subroutine acc4mean_surf(klon,klat,acccount,msvars,svar_surf,acum4mean_surf)
    !  purpose:  accumulate
    ! written by Ulf Hansson, Rossby centre, 980313
    ! Changes:

    implicit none  

    integer :: klon, klat  
    integer :: msvars  
    integer :: acccount(msvars)  
    real(kind=realkind) :: svar_surf(klon, klat, msvars)  
    real(kind=realkind) :: acum4mean_surf(klon, klat, msvars)  


    !     Accumulate all variables in svar_surf
    acum4mean_surf=acum4mean_surf+svar_surf
    acccount = acccount + 1  

    return  
  end subroutine acc4mean_surf

  subroutine accmean1(klon, klat, accvar, acccount)  
    !  purpose:  make mean from accumulated variable
    ! written by Ulf Hansson, Rossby centre, 050712
    ! Changes:

    implicit none  

    integer :: klon, klat, acccount  
    real(kind=realkind) :: accvar(klon, klat)  

    accvar = accvar/real(acccount,realkind)  

    return  
  end subroutine accmean1


  subroutine accmean_surf(klon, klat, acccount, msvars, kk, acum4mean_surf)
    !     purpose:  make mean from accumulated surface fluxes
    !     written by Patrick Samuelsson, Rossby centre, 040527

    implicit none  
    integer :: klon, klat  
    integer :: msvars, kk  
    integer :: acccount(msvars)  
    real(kind=realkind) :: acum4mean_surf(klon, klat, msvars)  
    acum4mean_surf(:,:, kk) = acum4mean_surf(:,:,kk)/real(acccount(kk),realkind)

    return  
  end subroutine accmean_surf


  subroutine accscratch1(klon, klat, accvar)  
    !  purpose:  reinitialise accumulated variable
    ! written by ulf hansson, rossby centre, 050712
    ! changes:

    implicit none  
    integer :: klon, klat  
    real(kind=realkind) :: accvar(klon, klat)  


    accvar = 0.0_realkind  

    return  
  end subroutine accscratch1

  subroutine accscratch4mean_surf(klon, klat, acccount, msvars, kk, &
       acum4mean_surf)
    !     purpose:  reinitialise accumulated fluxes  and counter
    !     based on accscratch. modified by patrick samuelsson 040527

    implicit none  

    integer,intent(in) :: klon, klat  
    integer,intent(in) :: msvars, kk  
    integer,intent(inout) :: acccount(msvars)  
    real(kind=realkind),intent(inout) :: acum4mean_surf(klon, klat, msvars)  
    acum4mean_surf(:,:,kk) = 0.0_realkind  
    acccount(kk) = 0  
    return  
  end subroutine accscratch4mean_surf


  subroutine accscratch_surf(klon, klat, msvars, kk, acum_surf)  
    !  purpose:  reinitialise accumulated fluxes  and counter
    ! based on accscratch. modified by patrick samuelsson 040527

    implicit none  

    integer :: klon, klat  
    integer :: msvars, kk  
    real(kind=realkind) :: acum_surf(klon, klat, msvars)  

    acum_surf(:, :, kk) = 0.0_realkind  
    return  
  end subroutine accscratch_surf


  subroutine acc_surf(klon, klat, msvars, svar_surf, acum_surf)  
    !     purpose:  accumulate
    !     written by ulf hansson, rossby centre, 980313
    !     changes:

    implicit none  

    integer :: klon, klat  
    integer :: msvars  
    real(kind=realkind) :: svar_surf(klon, klat, msvars)  
    real(kind=realkind) :: acum_surf(klon, klat, msvars)  

    !     accumulate all variables in svar_surf

    acum_surf = acum_surf + svar_surf 

    return  
  end subroutine acc_surf

  subroutine accumulate1(klon, klat, var, accvar)  
    !
    !  purpose:  accumulate variable
    !
    ! written by ulf hansson, rossby centre, 050712
    !
    ! changes:
    !
    !
    !

    implicit none  

    integer :: klon, klat  
    real(kind=realkind) :: var(klon, klat)  
    real(kind=realkind) :: accvar(klon, klat)  

    accvar = accvar + var
    return  
  end subroutine accumulate1



  subroutine oasacc(klon, klat, acccount, t2ms, acct2ms, t2mi, &
       acct2mi, senfs, accsenfs, senfi, accsenfi, latfs, acclatfs, latfi, &
       acclatfi, slwrs, accslwrs,slwri, accslwri, sswrs, accsswrs,sswri, accsswri, slwdn, accslwrdown, &
       sswdn, accsswrdown, q2m, accq2m, ps, accps)
!       sswdn, accsswrdown, q2m, accq2m, cov2d, acccov2d, ps, accps)
    !  purpose:  accumulate varuiables for oasis couupler
    ! written by ulf hansson, rossby centre, 001005
    ! based on accumulate

    implicit none  

    integer :: klon, klat  


    real(kind=realkind) :: t2ms(klon, klat), acct2ms(klon, klat), t2mi(klon, klat) &
         , acct2mi(klon, klat), senfs(klon, klat), accsenfs(klon, klat), &
         senfi(klon, klat), accsenfi(klon, klat), latfs(klon, klat), &
         acclatfs(klon, klat), latfi(klon, klat), acclatfi(klon, klat), &
         slwrs(klon, klat), accslwrs(klon, klat),slwri(klon, klat), accslwri(klon, klat), &
          sswrs(klon, klat),sswri(klon, klat), &
         accsswrs(klon, klat), accsswri(klon, klat),slwdn(klon, klat), accslwrdown(klon, &
         klat), sswdn(klon, klat), accsswrdown(klon, klat), q2m(klon, &
         klat), accq2m(klon, klat), &
         ps(klon, klat), accps(klon, klat)

    integer :: acccount  

    logical,save :: lfirst=.true.


    if(lfirst) then         ! initialise
       acct2ms   = 0.0_realkind  
       acct2mi   = 0.0_realkind  
       accsenfs  = 0.0_realkind  
       accsenfi  = 0.0_realkind  
       acclatfs  = 0.0_realkind  
       acclatfi  = 0.0_realkind  
       accslwrs  = 0.0_realkind  
       accsswrs  = 0.0_realkind  
       accslwri  = 0.0_realkind
       accsswri  = 0.0_realkind

       accslwrdown  = 0.0_realkind  
       accsswrdown  = 0.0_realkind  
       accq2m  = 0.0_realkind  
!      acccov2d  = 0.0_realkind  
       accps  = 0.0_realkind  
       acccount = 0  

       lfirst = .false.  
    endif

    acct2ms  = acct2ms  + t2ms   
    acct2mi  = acct2mi  + t2mi   
    accsenfs  = accsenfs  + senfs   
    accsenfi  = accsenfi  + senfi   
    acclatfs  = acclatfs  + latfs   
    acclatfi  = acclatfi  + latfi   
    accslwrs  = accslwrs  + slwrs   
    accsswrs  = accsswrs  + sswrs   
    accslwri  = accslwri  + slwri
    accsswri  = accsswri  + sswri

    accslwrdown  = accslwrdown  + slwdn   
    accsswrdown  = accsswrdown  + sswdn   
    accq2m  = accq2m  + q2m   
!    acccov2d  = acccov2d  + cov2d   
    accps  = accps  + ps   

    acccount = acccount + 1  
    return  
  end subroutine oasacc


  subroutine oasaccmean(klon,klat,acccount,acct2ms,&
       acct2mi,accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri, &
       accq2m,accps,accslwrdown,accsswrdown )

    !     purpose:  make mean from accumulated oasis variables 
    !     written by Ulf Hansson, Rossby centre, 001005
    !     based on accmean



    implicit none
    integer klon,klat
    real(kind=realkind)::acct2ms(klon,klat),acct2mi(klon,klat) , &
         accsenfs(klon,klat),accsenfi(klon,klat) ,&
         acclatfs(klon,klat),acclatfi(klon,klat) ,&
         accslwrs(klon,klat),accsswrs(klon,klat), &
         accslwri(klon,klat),accsswri(klon,klat), &
         accq2m(klon,klat),  &
         accps(klon,klat)
    real(kind=realkind):: accslwrdown(klon,klat),accsswrdown(klon,klat)
!    real(kind=realkind):: zconacc,dtime
    integer acccount

    acct2ms = acct2ms/real(acccount,realkind)
    acct2mi = acct2mi/real(acccount,realkind)
    accsenfs = accsenfs/real(acccount,realkind)
    accsenfi = accsenfi/real(acccount,realkind)
    acclatfs = acclatfs/real(acccount,realkind)
    acclatfi = acclatfi/real(acccount,realkind)
    accslwrs = accslwrs/real(acccount,realkind)
    accslwri = accslwri/real(acccount,realkind)
    accsswrs = accsswrs/real(acccount,realkind)
    accsswri = accsswri/real(acccount,realkind)
    accq2m = accq2m/real(acccount,realkind)
!    acccov2d = acccov2d/real(acccount,realkind)
    accps = accps/real(acccount,realkind)
    accslwrdown = accslwrdown/real(acccount,realkind)
    accsswrdown = accsswrdown/real(acccount,realkind)


    return
  end subroutine oasaccmean


  subroutine oasaccscratch(klon, klat, acccount, acct2ms, acct2mi, &
       accsenfs, accsenfi, acclatfs, acclatfi, accslwrs,accslwri, accsswrs,accsswri, &
       accslwrdown, accsswrdown, accq2m, accps)
    !  purpose:  reinitialise accumulated oasis variables and counter
    ! written by ulf hansson, rossby centre, 001005
    ! based on accscratch

    implicit none  

    integer :: klon, klat  
    real(kind=realkind)::acct2ms(klon, klat), acct2mi(klon, klat),&
         accsenfs(klon,klat),accsenfi(klon,klat),acclatfs(klon, klat), &
         acclatfi(klon, klat), accslwrs(klon, klat), accsswrs(klon, &
         klat), accslwrdown(klon, klat), accsswrdown(klon, klat), &
         accq2m(klon, klat), accps(klon, klat),  &
         accslwri(klon, klat), accsswri(klon,klat)
    integer :: acccount  


    acct2ms  = 0.0_realkind  
    acct2mi  = 0.0_realkind  
    accsenfs  = 0.0_realkind  
    accsenfi  = 0.0_realkind  
    acclatfs  = 0.0_realkind  
    acclatfi  = 0.0_realkind  
    accslwrs  = 0.0_realkind  
    accslwri  = 0.0_realkind  
    accsswrs  = 0.0_realkind  
    accsswri  = 0.0_realkind  

    accslwrdown  = 0.0_realkind  
    accsswrdown  = 0.0_realkind  
    accq2m  = 0.0_realkind  
!    acccov2d  = 0.0_realkind  
    accps  = 0.0_realkind  
    
    acccount = 0  

    return  
  end subroutine oasaccscratch

end module accumulate
