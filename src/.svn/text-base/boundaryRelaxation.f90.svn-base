module boundaryRelaxation
  use timetype
  use calendar
  use derived_types
  use decomp
  implicit none
  private
  
  integer,public,save::nbdpts=8 ! number of gridpoints in the boundary relaxation zone for the boundary relaxation
  integer,public,save::npbpts=2 ! number of extra (passive) boundary lines  where the boundary relaxation function
  
  logical,public,save::nltanh = .true.! .true. if tanh-shape boundary relaxation function

  real(kind=realkind),allocatable,dimension(:,:):: weightu,weightv,weightm 
  integer:: p450,p750,nbdpts_min,nbdpts_max
  real(kind=realkind):: pdep
  integer::pdel
  real(kind=realkind),allocatable,dimension(:,:,:):: weightuz,weightvz,weightmz 
  real(kind=realkind),allocatable,dimension(:)::pfull
  
  
  public bndrel,bdinit,bdy_swap,initBoundary
contains


  subroutine initBoundary(klon,klat,klev,ahyb_in,bhyb_in,initTime)
    use decomp

    use referenceParameters, only:sip0
    implicit none
    type(time),intent(in)::initTime
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(in),dimension(klev+1)::ahyb_in,bhyb_in
    real(kind=realkind)::afull(klev),bfull(klev)
    integer::jk

    namelist/nambc/npbpts,nbdpts,nltanh 
    
    open(57,file='namelists.dat',status='old')
    read(57,nml=nambc)
    close(57)
    if(mype==0)then
       write(6,nml=nambc)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=nambc)
       close(1)
    endif

    allocate(weightu(klon,klat),weightv(klon,klat),weightm(klon,klat))
    allocate(pfull(klev))                        
    allocate(weightuz(klon,klat,klev),weightvz(klon,klat,klev))
    allocate(weightmz(klon,klat,klev))      
    p450 = 1
    p750=1
    do jk=1,klev
       afull(jk) = 0.5_realkind*( ahyb_in(jk) + ahyb_in(jk+1) )
       bfull(jk) = 0.5_realkind*( bhyb_in(jk) + bhyb_in(jk+1) )
       pfull(jk)=afull(jk)+(bfull(jk)*sip0)
       if(pfull(jk)<=40000.0_realkind)p450=jk
       if(pfull(jk)<=70000.0_realkind)p750=jk
    enddo
    pdel=p750-p450

    nbdpts_min=nbdpts
    nbdpts_max=nbdpts*2

    call bdinit(klon,klat,nbdpts,npbpts,weightu,weightv,weightm,nltanh)

    do jk=1,klev
!!$     if(jk<=p450)nbdpts=nbdpts_max
!!$     if(jk>=p750)nbdpts=nbdpts_min
!!$     if(jk>p450 .and. jk<p750)then
!!$        pdep=jk-p450
!!$        nbdpts=int(nbdpts_max-((pdel/(nbdpts_max-nbdpts_min))*pdep))
!!$        nbdpts=max(nbdpts_min,(min(nbdpts,nbdpts_max)))
!!$     endif
!!$     call bdinit (klon, klat, nbdpts, npbpts, weightu, weightv &
!!$          , weightm,  nltanh , nlpwei   )

       weightuz(:,:,jk)=weightu
       weightvz(:,:,jk)=weightv
       weightmz(:,:,jk)=weightm
    enddo                     !	end jk vertical loop
    !  nbdpts=nbdpts_min



  end subroutine initBoundary

  subroutine bdinit(klon,klat,kbdpts,npass,weightu,weightv,weightm,ltanh)
    use decomp
    implicit none
    integer,intent(in)::klon,klat
    integer,intent(in)::kbdpts,npass
    logical,intent(in)::ltanh
    real(kind=realkind),dimension(klon,klat),intent(out)::weightu,weightv,weightm
    real(kind=realkind)::lonfun(klon_global),latfun(klat_global)
    integer::i,j

    call get_relax_fun(klon_global,lonfun,npass,kbdpts,3.0_realkind,ltanh,.true.)
    call get_relax_fun(klat_global,latfun,npass,kbdpts,3.0_realkind,ltanh,.true.)
    
    do j=1,klat
       weightu(:,j) = lonfun(idatastart:idatastart+klon-1)
    enddo
    do i=1,klon
       weightv(i,:) = latfun(jdatastart:jdatastart+klon-1)
    enddo

    call get_relax_fun(klon_global,lonfun,npass,kbdpts,3.0_realkind,ltanh,.false.)
    call get_relax_fun(klat_global,latfun,npass,kbdpts,3.0_realkind,ltanh,.false.)

    weightm = 1.0_realkind

    do j=1,klat
       do i=1,klon
          weightm(i,j) = max(lonfun(idatastart+i-1),latfun(jdatastart+j-1))
          weightu(i,j) = max(weightu(i,j),latfun(jdatastart+j-1))
          weightv(i,j) = max(weightv(i,j),lonfun(idatastart+i-1))
       enddo
    enddo

  end subroutine bdinit

  subroutine bdy_swap(bcout,bcin)
    !     purpose:
    !     move the boundary fields from one arrays to the others
    !     description:
    !     all the boundary fields are copied from time level 2 to time
    !     level 1 before reading in new boundary fileds
    use derived_types
    implicit none  

    type(atm),intent(in)::bcin
    type(atm),intent(inout)::bcout
    bcout = bcin
    return  
  end subroutine bdy_swap

subroutine bndrel(klon,klat,klev,ksvar,phip,bc,lsimp,lslan,lcond,lcw,lsv)

    use derived_types
    use calendar

    implicit none

    integer,intent(in):: klon, klat, klev, ksvar
    type(atm),intent(in)::bc(2)
    type(atm),intent(inout)::phip
    logical,intent(in):: lsimp , lslan, lcw(2), lsv(ksvar,2),lcond


    integer:: jx,jy,jk,js
    real(kind=realkind):: zt,zw
    integer:: lrelcw
    real(kind=realkind)::timrat !expressed in seconds!

    real(kind=realkind)::dtbc
    if(bc(2)%timestamp<bc(1)%timestamp.or.bc(2)%timestamp==bc(1)%timestamp)then
       print *,bc(1)%timestamp
       print *,bc(2)%timestamp
       if(bc(2)%timestamp<bc(1)%timestamp)then
          call stop_program( 'bc1 is larger than bc2 in TIME!!!')
       endif
       if(bc(2)%timestamp==bc(1)%timestamp)then
          call stop_program( 'bc1 is equal bc2 in TIME!!!')
       endif
    endif
    dtbc = real(Nseconds(bc(2)%timestamp - bc(1)%timestamp),realkind) !bc time frequency
    timrat = real(Nseconds(phip%timestamp-bc(1)%timestamp),realkind)/dtbc
    
    if(timrat<0._realkind.or.timrat>1._realkind)then
       print *,timrat
       call stop_program( 'timraterror relaxation')
    endif
!!$    if(mype==0)then
!!$       print *,'BNDREL',timrat
!!$       print *,bc(1)%timestamp
!!$       print *,bc(2)%timestamp
!!$    endif

    !         interpolate and relax single level fields (ps/lnps)
    if(lsimp) then
       do jy=1,klat
          do jx=1,klon
             zw = bc(1)%lnps(jx,jy) + timrat*(bc(2)%lnps(jx,jy)-bc(1)%lnps(jx,jy))
             phip%lnps(jx,jy)=bsf(phip%lnps(jx,jy),zw,weightmz(jx,jy,klev))
          enddo
       enddo
       phip%ps = exp(phip%lnps)
    else
       do jy=1,klat
          do jx=1,klon
             zw = bc(1)%ps(jx,jy) + timrat*(bc(2)%ps(jx,jy)-bc(1)%ps(jx,jy))
             phip%ps(jx,jy)=bsf(phip%ps(jx,jy),zw,weightmz(jx,jy,klev))
          enddo
       enddo
       do jy=1,klat
          do jx=1,klon
             phip%lnps(jx,jy) = log(phip%ps(jx,jy))
          enddo
       enddo
    endif

    !     interpolate and relax multi-level fields (u,v,t,q,s,edot)

    do jk=1,klev
       do jy=1,klat
          do jx=1,klon
             zw = bc(1)%u(jx,jy,jk) + timrat*(bc(2)%u(jx,jy,jk)-bc(1)%u(jx,jy,jk))
             phip%u(jx,jy,jk)=bsf(phip%u(jx,jy,jk),zw,weightuz(jx,jy,klev))
             zw = bc(1)%v(jx,jy,jk) + timrat*(bc(2)%v(jx,jy,jk)-bc(1)%v(jx,jy,jk))
             phip%v(jx,jy,jk)=bsf(phip%v(jx,jy,jk),zw,weightvz(jx,jy,klev))
             zw = bc(1)%t(jx,jy,jk) + timrat*(bc(2)%t(jx,jy,jk)-bc(1)%t(jx,jy,jk))
             phip%t(jx,jy,jk)=bsf(phip%t(jx,jy,jk),zw,weightmz(jx,jy,klev))
             zw = bc(1)%q(jx,jy,jk) + timrat*(bc(2)%q(jx,jy,jk)-bc(1)%q(jx,jy,jk))
             phip%q(jx,jy,jk)=bsf(phip%q(jx,jy,jk),zw,weightmz(jx,jy,klev))
          enddo
       enddo

       !     gjkf do not relax cloud water when running with era boundary conditions
       !     gjkf later put in an (if.lecmwf.)test
       !     gjkf
       !     do not relax cloud water field in the case of eulerian with
       !     straco condensation scheme due to overspecification
       !     of boundary conditions - currently cloud water is
       !     relaxed for other condensation schemes if both
       !     fields are there.

       lrelcw=1
       if(lrelcw==0)then
          write(*,*)'relaxing cloud water in bndrel!!'
          if ( lcw(1) .and. lcw(2)) then
             if ( lslan .or. .not.lcond ) then
                do jy=1,klat
                   do jx=1,klon
                      zw = bc(1)%cw(jx,jy,jk) + timrat*(bc(2)%cw(jx,jy,jk)-bc(1)%cw(jx,jy,jk))  
                      phip%cw(jx,jy,jk)=bsf(phip%cw(jx,jy,jk),zw,weightmz(jx,jy,klev))
                   enddo
                enddo
             endif            ! end test lslan .or. .not.lcond
          endif               ! end test lcw(1) .and. lcw(2)
          !     gjkf
       endif                  !endif hop around cw relaxation

       if (lslan) then
          do jy=1,klat
             do jx=1,klon
                zw = bc(1)%edot(jx,jy,jk)+timrat*(bc(2)%edot(jx,jy,jk)-bc(1)%edot(jx,jy,jk))
                phip%edot(jx,jy,jk)=bsf(phip%edot(jx,jy,jk),zw,weightmz(jx,jy,klev))
             enddo
          enddo
       endif

       !         extra scalars

       !     parameter is relaxed if both fields are there
       do js = 1,ksvar
          if(lsv(js,1) .and. lsv(js,2)) then
             do jy = 1,klat
                do jx = 1,klon
                   zw = bc(1)%svar(jx,jy,jk,js)+ timrat*(bc(2)%svar(jx,jy,jk,js)-bc(1)%svar(jx,jy,jk,js))
                   phip%svar(jx,jy,jk,js)=bsf(phip%svar(jx,jy,jk,js),zw,weightmz(jx,jy,klev))
                enddo
             enddo
          endif
       enddo
    enddo                     ! end jk loop over levels

    !          upper/lower boundary condition for vertical velocity
    if (lslan) then
       do jy=1,klat
          do jx=1,klon
             phip%edot(jx,jy,1) = 0.0_realkind
             phip%edot(jx,jy,klev+1) = 0.0_realkind
          enddo
       enddo
    endif
  end subroutine bndrel

  real(kind=realkind) function bsf(pa,pb,pw)
    real(kind=realkind),intent(in)::pa,pb,pw
    bsf = (1.0_realkind - pw)*pa + pw*pb
  end function bsf
  real(kind=realkind) function relaxFunction(x,N)
    !The relaxfunction starts at x==1 (1.5) and ends at x==N (N+0.5) 
    implicit none
    real(kind=realkind),intent(in)::x
    integer,intent(in)::N
    real(kind=realkind),parameter::a=0.2675_realkind,tol=0.001_realkind
    real(kind=realkind)::f
    
    if(N>18 .or. N<7)then
       print *,'N=',N
       call stop_program( 'relaxation function has not been tuned for more than 18 points or less than 7 points')
    endif
    f = 0.5_realkind*(1.0_realkind-tanh((2.0_realkind*x-real(N,realkind)-1.5_realkind)/(a*real(N,realkind)))) !this works upto 14 relaxation points
    if(f<tol)then
       f = 0.0_realkind
    endif
    if(abs(f-1.0_realkind)<tol)then
       f=1.0_realkind
    endif
    relaxFunction = f
    return
  end function relaxFunction

  subroutine get_relax_fun(n,fun,npass,nrelpts,shiftC,ltanh,lstaggered)
    implicit none
    integer,intent(in)::n,npass,nrelpts
    real(kind=realkind),intent(in)::shiftC
    logical,intent(in)::ltanh,lstaggered
    real(kind=realkind),intent(inout)::fun(n)
    real(kind=realkind)::acon
    integer::jk
    real(kind=realkind)::x

    if(ltanh)then
       acon=2.0_realkind/(2.0_realkind*(real(nrelpts,realkind)-1.0_realkind)-4.0_realkind)
    else
       acon=4._realkind*atan(1.0_realkind)/(2.0_realkind*real(nrelpts,realkind)-2.0_realkind)
    endif

    fun = 1.0_realkind
    fun(npass+nrelpts+1:n-npass-nrelpts)=0.0_realkind !inner set to 0.0

    if(ltanh)then
       do jk=npass+1,npass+nrelpts+1
          x = real(jk-npass,realkind)
          if(lstaggered) x=x+0.5_realkind
          fun(jk) = relaxFunction(x,nrelpts)
          if(lstaggered)then
             fun(n-jk) = fun(jk)
          else
             fun(n-jk+1) = fun(jk)
          endif
       enddo
    else
       do jk = npass+2,npass+nrelpts 
          !fun(jk)=0.5*(1.+cos(acon*(2.0*(jk-npass)-3.0)))
          fun(jk)=0.5_realkind*(1.0_realkind+cos(acon*(2.0_realkind*real(jk-npass,realkind)-shiftC)))
          fun(n-jk+1) = fun(jk)
       enddo
       fun(npass+1) = 1.0_realkind
       fun(n-npass) = 1.0_realkind
    endif


  end subroutine get_relax_fun



end module boundaryRelaxation
