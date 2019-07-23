module hybridmod
  use confys
  use decomp,only:realkind
  implicit none
  private

  public ahybrid
contains
  subroutine ahybrid (nhor,nlev,kstart,kstop, &
       t,q,                                     &
       ps,                                      &
       gpot,ph,pf,dph,dpf,                      &
       ahyb,bhyb,hybf)
    implicit none

!!$    integer,intent(in):: nhor,nlev,kstart,kstop
!!$    real(kind=realkind),intent(in):: t(nhor,nlev),q(nhor,nlev)
!!$    real(kind=realkind),intent(out)::gpot(nhor,nlev),ph(nhor,nlev+1),&
!!$         pf(nhor,nlev),dph(nhor,nlev+1),dpf(nhor,nlev)
!!$    real(kind=realkind),intent(in):: ps(nhor)
!!$    real(kind=realkind),intent(in):: ahyb(nlev+1),bhyb(nlev+1),hybf(nlev)
!!$
!!$    call hybrid (nhor,nlev,kstart,kstop, &
!!$         t,q,                              &
!!$         ps,                               &
!!$         gpot,ph,pf,dph,dpf,               &
!!$         ahyb,bhyb,hybf)
!!$
!!$    return
!!$  end subroutine ahybrid
!!$
!!$
!!$  subroutine hybrid (nhor,nlev,kstart,kstop, &
!!$       t,q,                                   &
!!$       ps,                                    &
!!$       gpot,ph,pf,dph,dpf,                    &
!!$       ahyb,bhyb,hybf)
!!$    implicit none
    integer,intent(in):: nhor,nlev,kstart,kstop
    integer:: jk,jl,jkp1
    real(kind=realkind):: zcons1,zgpot,zpfll,zphalf,ztmvir,zcrdq,zphnp1
    real(kind=realkind),intent(in):: t(nhor,nlev),q(nhor,nlev)
    real(kind=realkind),intent(out)::gpot(nhor,nlev),ph(nhor,nlev+1),&
         pf(nhor,nlev),dph(nhor,nlev+1),dpf(nhor,nlev)
    real(kind=realkind),intent(in):: ps(nhor)
    real(kind=realkind),intent(in):: ahyb(nlev+1),bhyb(nlev+1),hybf(nlev)
    !     l    purpose:
    !     l    calculate pressure ('pf' and 'ph' ) in the middle of
    !     l    and at the intersections of model layers. also corre-
    !     l    sponding pressure thicknesses ( 'dpf' and 'dph' ) as
    !     l    well as geopotential ( 'gpot' ) in the middle of model
    !     l    layers are determined.
    !     l    method:
    !     l    see 'hirlam documentation manual', 1988.
    !     variable        type            content
    !     --------        ------          --------------------------------
    !************************integer input ********************************
    !     l
    !     l      nhor            input           dimension length in the horizont
    !     l
    !     l      nlev            input           number of vertical levels
    !     l
    !     l      kstart          input           starting index for horizontal lo
    !     l
    !     l      kstop           input           ending index for hor. loops
    !     
    !************************full fields input ****************************
    !     l
    !     l      t(nhor,nlev)    input           temperatures
    !     l
    !            q(nhor,nlev)    input           specific humidities
    !     
    !************************horizontal fields input **********************
    !     l      ps(nhor)        input           surface pressure
    !     l
    !************************1-d vertical arrays input ********************
    !     l
    !     l      ahyb,bhyb,hybf   input           vertical coordinate data
    !     l
    !************************full fields output ***************************
    !     l
    !     l      gpot(nhor,nlev) input           virtual geopotential
    !     l      pf(nhor,nlev)   input           pressure on "full" levels
    !     l      ph(nhor,nlev+1) input           pressure on "half" levels
    !     l      dpf(nhor,nlev)  input           pressure diff. on "full" levels
    !     l      dph(nhor,nlev+1 input           pressure diff. on "full" levels
    !     l
    !***********************************************************************
    !--------------------------------------------------------------------


    !     --------------------------------------------------------
    !     l          1.   determine 'pf','ph','dpf' and 'dph'.

    zcons1 = -rair*log(hybf(nlev))
    zcrdq = 1._realkind /epsilo-1._realkind 

    do 136 jl=kstart,kstop
       zphnp1 = ahyb(nlev+1)+bhyb(nlev+1)*ps(jl)
       zphalf = ahyb(nlev)+bhyb(nlev)*ps(jl)
       zpfll = ( zphalf + zphnp1 )*0.5_realkind 
       dpf(jl,nlev) = zphnp1 - zphalf
       ph(jl,nlev) = zphalf
       ph(jl,nlev+1) = zphnp1
       dph(jl,nlev+1) = ps(jl)-zpfll
       gpot(jl,nlev) = zcons1*t(jl,nlev)*(1._realkind +zcrdq*q(jl,nlev))
       pf(jl,nlev) = zpfll
136 enddo

    do 140 jk=nlev-1,1,-1

       jkp1=jk+1

       do 138 jl=kstart,kstop
          zphalf=ahyb(jk)+bhyb(jk)*ps(jl)
          zpfll=( zphalf + ph(jl,jkp1) )*0.5_realkind 
          dpf(jl,jk)=ph(jl,jkp1) - zphalf
          ph(jl,jk)=zphalf

          ztmvir=0.5_realkind *(t(jl,jk)*(1._realkind +zcrdq*q(jl,jk))+&
               t(jl,jkp1)*(1._realkind +zcrdq*q(jl,jkp1)))
          zgpot=gpot(jl,jkp1)+rair*ztmvir*log(pf(jl,jkp1)/zpfll)
          gpot(jl,jk)=zgpot

          pf(jl,jk)=zpfll
          dph(jl,jkp1)=pf(jl,jkp1)-zpfll

138    enddo

140 enddo

    do 150 jl=kstart,kstop
       dph(jl,1) = pf(jl,1)
150 enddo

    return
  end subroutine ahybrid
!  end subroutine hybrid
end module hybridmod
