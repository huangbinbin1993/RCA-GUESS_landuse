module derived_types
  use timetype
  use calendar
  use decomp,only:realkind
  implicit none
  private


!!$  type,public::ecoclim
!!$     type(time)::timestamp
!!$     real(kind=realkind),allocatable,dimension(:,:)::
!!$     real(kind=realkind),pointer,dimension(:,:)::eco(:)
!!$  end type ecoclim
  
  type,public::atm
     type(time)::timestamp
     real(kind=realkind),allocatable,dimension(:,:)::ps,lnps
     real(kind=realkind),allocatable,dimension(:,:,:)::u,v,t,q,cw,edot
     real(kind=realkind),allocatable,dimension(:,:,:,:)::svar
  end type atm

  interface operator(*)
     module procedure mult_atm_real
  end interface

  public operator(*) 
  public alloc_atm,dealloc_atm

contains
  

  type(atm) function mult_atm_real(z,at)
    type(atm), intent(in)::at
    real(kind=realkind),intent(in)::z
    mult_atm_real%ps = z*at%ps
    mult_atm_real%lnps = z*at%lnps
    mult_atm_real%u = z*at%u
    mult_atm_real%v = z*at%v
    mult_atm_real%t = z*at%t
    mult_atm_real%cw = z*at%cw
    mult_atm_real%edot = z*at%edot
    mult_atm_real%q = z*at%q
    mult_atm_real%svar = z*at%svar
    mult_atm_real%timestamp = at%timestamp
  end function mult_atm_real

  subroutine assign_atm(field1,from)
    implicit none 
    type(atm), intent(in)::from
    type(atm),intent(inout)::field1
    field1%ps = from%ps
    field1%lnps = from%lnps
    field1%u = from%u
    field1%v = from%v
    field1%t = from%t
    field1%cw = from%cw
    field1%edot = from%edot
    field1%q = from%q
    field1%svar = from%svar
    field1%timestamp = from%timestamp
  end subroutine assign_atm

  type(atm) function alloc_atm(klon,klat,klev,ksvar,t)
    implicit none
    type(time),intent(in)::t
    integer,intent(in)::klon,klat,klev,ksvar
    alloc_atm%timestamp = t
    allocate(alloc_atm%ps(klon,klat))
    allocate(alloc_atm%lnps(klon,klat))
    allocate(alloc_atm%u(klon,klat,klev))
    allocate(alloc_atm%v(klon,klat,klev))
    allocate(alloc_atm%t(klon,klat,klev))
    allocate(alloc_atm%q(klon,klat,klev))
    allocate(alloc_atm%cw(klon,klat,klev))
    allocate(alloc_atm%edot(klon,klat,klev+1))
    allocate(alloc_atm%svar(klon,klat,klev,ksvar))

    alloc_atm%ps = -666.0_realkind
    alloc_atm%lnps = -666.0_realkind
    alloc_atm%u = -666.0_realkind
    alloc_atm%v = -666.0_realkind
    alloc_atm%t = -666.0_realkind
    alloc_atm%q = -666.0_realkind
    alloc_atm%cw = -666.0_realkind
    alloc_atm%edot = -666.0_realkind
    alloc_atm%svar = -666.0_realkind
  end function alloc_atm


  subroutine dealloc_atm(phi)
    implicit none
    type(atm),intent(inout)::phi
    phi%timestamp = maketime(-1,-1,-1,-1,-1,-1)
    deallocate(phi%ps)
    deallocate(phi%lnps)
    deallocate(phi%u)
    deallocate(phi%v)
    deallocate(phi%t)
    deallocate(phi%q)
    deallocate(phi%cw)
    deallocate(phi%edot)
    deallocate(phi%svar)

  end subroutine dealloc_atm


  type(atm) function alloc_atmCopy(at)
    type(atm),intent(in)::at
    integer::klon,klat,klev,ksvar
    klon = ubound(at%ps,1)
    klat = ubound(at%ps,2)
    klev = ubound(at%u,3)
    ksvar = ubound(at%svar,4)
    alloc_atmCopy = alloc_atm(klon,klat,klev,ksvar,at%timestamp)
    alloc_atmCopy = at
  end function alloc_atmCopy



    
end module derived_types
