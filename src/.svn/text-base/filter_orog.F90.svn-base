module filter_orog
  use decomp
  implicit none
  private
  public filter_orography
contains
  subroutine filter_orography(fld,klon,klat)

    implicit none
    integer,intent(in)::klon,klat
    real(kind=realkind)::fld(klon,klat)
    real(kind=realkind),allocatable,dimension(:,:)::tmp
    integer:: i,j,rank
    real(kind=realkind)::epsilon
#ifdef MPI_SRC
    real(kind=realkind),allocatable,dimension(:,:)::f_global
#include"mpif.h"
    integer::ierr
#endif
    integer::status
    !real(kind=realkind),parameter::ost=1d0/8d0,four=2d0!gives classic 1 2 1 filter
    real(kind=realkind),parameter::ost=1d0/12d0,four=4d0 !1 4 1
    
    allocate(tmp(klon,klat))
    
    epsilon=1.0_realkind

    tmp = fld
    do j=2,klat-1
       do i=2,klon-1
          tmp(i,j) = ost*(fld(i-1,j)+four*fld(i,j)+fld(i+1,j) )+ &
               ost*(fld(i,j-1)+four*fld(i,j)+fld(i,j+1) )
       enddo
    enddo
    

    !call swap(tmp,klon,klat) 
    call swap2d(tmp,klon,klat)
    fld = tmp
    deallocate(tmp)

!!$    if(nproc==1) then
!!$       call orogfil(fld,epsilon,klon_global,klat_global)
!!$#ifdef MPI_SRC       
!!$    else
!!$       allocate(f_global(klon_global,klat_global),stat=status)
!!$       if(status/=0)then
!!$          stop 'memory allocation error klon_global*klat_global too large?'
!!$       endif
!!$       call colfld(0,f_global,fld,klon,klat)
!!$       if(mype==0)then
!!$          call orogfil(f_global,epsilon,klon_global,klat_global)
!!$       endif
!!$
!!$
!!$       call MPI_Bcast(f_global,klon_global*klat_global,REALTYPE,0, &
!!$            localComm,ierr)
!!$   
!!$
!!$       call myshare(f_global,fld,klon,klat)
!!$       deallocate(f_global)
!!$#endif  
!!$    endif


    return
  end subroutine filter_orography




  subroutine orogfil(f, epsilon,klon_global,klat_global)
    !     implicit filtering of a two dimensional field using
    !     the sixth-order low-pass implicit tangent filter
    !     described in raymond, m.w.r., 116, 2132-2141.
    !     epsilon:         ----  filter parameter to determine cutoff
    !     where filter response is one-half.
    implicit none  
    integer,intent(in)::klon_global,klat_global
    real(kind=realkind),intent(inout)::f(klon_global,klat_global)
    real(kind=realkind),intent(in)::epsilon
    real(kind=8):: edubl  
    real(kind=8),allocatable::fdubl(:)
    integer :: j, i 

    allocate(fdubl(max(klon_global,klat_global)))

    edubl = dble(epsilon)

    !filter the rows in the x-direction.
    do j = 1,klat_global
       fdubl(1:klon_global) = dble(f(:, j))  
       call lowpas(fdubl, klon_global, edubl)  
       f(:, j) = real(fdubl(1:klon_global),realkind)
    enddo

    !filter the rows in the y-direction.
    do i = 1,klon_global
       fdubl(1:klat_global) = dble(f(i, :))
       call lowpas(fdubl, klat_global, edubl)  
       f(i, :) = real(fdubl(1:klat_global),realkind)
    enddo

    return  

  end subroutine orogfil


  subroutine lowpas(xy,n,eps)

    !     SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
    !     (RAYMOND, MWR,116,2132-2141)
    !     
    !*************************************************************
    !***  THIS CODE IS COPIED FROM A LISTING PROVIDED       ***
    !***  BY WILLIAM H RAYMOND. SOME NOTATIONAL CHANGES     ***
    !***  HAVE BEEN MADE IN THE ROUTINE LOWPAS. THE         ***
    !***  ROUTINE INVLOW HAS BEEN COPIED ALMOST VERBATIM.   ***
    !*************************************************************
    !     
    !     XY     UNFILTERED VALUES ON INPUT
    !     FILTERED VALUES ON OUTPUT.
    !     N      NUMBER OF VALUES.
    !     EPS    FILTER PARAMETER 
    !     (DETERMINES CUTOFF)
    implicit none
    integer:: n,nm1,nm2,nm3,nm4,j
    real(kind=8)::eps
    real(kind=8)::xy(n),rhs(n),xdash(n)
    logical,save::iskip=.false.



    nm1 = n-1
    nm2 = n-2
    nm3 = n-3
    nm4 = n-4

    rhs(1) = 0d0
    rhs(2) = eps*(xy(1)-2.d0*xy(2)+xy(3))
    rhs(3) = eps*(-1.d0*(xy(1)+xy(5)) &
         +4.d0*(xy(2)+xy(4)) &
         -6.d0*xy(3))

    do j=4,nm3
       rhs(j) = eps*((xy(j-3)+xy(j+3)) &
            - 6.d0*(xy(j-2)+xy(j+2)) &
            +15.d0*(xy(j-1)+xy(j+1)) &
            -20.d0* xy(  j)         )
    enddo

    rhs(nm2) = eps*(-1.d0*(xy(  n)+xy(nm4)) &
         +4.d0*(xy(nm1)+xy(nm3)) &
         -6.d0* xy(nm2)         )
    rhs(nm1) = eps*(xy(nm2)-2.d0*xy(nm1)+xy(  n)) 
    rhs(n) = 0d0     
    !     solve the linear system for xdash
    call invlow(rhs,n,xdash,eps,iskip)

    !     add correction to get filtered values.
    do j=1,n
       xy(j) = xy(j) + xdash(j)
    enddo

  end subroutine lowpas

  subroutine invlow(bb,n,xans,ep,iskip)

    !     gaussian elimination for low-pass filter.

    !     sixth-order low-pass implicit tangent filter.
    !     (ref: william h raymond, mwr, 116, 2132-2124)

    implicit none
    real(kind=8):: ep
    integer:: n,i
    logical:: iskip
    real(kind=8):: a(n),b(n),c(n),d(n),e(n), &
         delta(n),beta(n),w(n),gam(n), &
         h(n),xans(n),bb(n),pi(n), &
         ap(n),f(n),z(n)

    if (.not.iskip)then      
       !     initialize the matrix
       do 10 i=4,n-3
          z(i) = 1.d0-ep
          a(i) = 6.d0*(1.d0+ep)
          b(i) = 15.d0*(1.d0-ep)
          c(i) = 20.d0*(1.d0+ep)
          d(i) = b(i)
          e(i) = a(i)
          f(i) = z(i)
10     enddo

       z(1) = 0d0
       z(2) = 0d0
       z(3) = 0d0

       a(1) = 0d0
       a(2) = 0d0
       a(3) = 1.d0+ep

       b(1) = 0d0
       b(2) = 1.d0-ep
       b(3) = 4.d0*(1.d0-ep)

       c(1) = 1.d0
       c(2) = 2.d0*(1.d0+ep)
       c(3) = 6.d0*(1.d0+ep)

       d(1) = 0d0
       d(2) = 1.d0-ep
       d(3) = 4.d0*(1.d0-ep)

       e(1) = 0d0
       e(2) = 0d0
       e(3) = 1.d0+ep

       f(1) = 0d0
       f(2) = 0d0
       f(3) = 0d0


       z(n-2) = 0d0 
       z(n-1) = 0d0
       z(n) = 0d0

       a(n-2) = 1.d0+ep
       a(n-1) = 0d0
       a(n) = 0d0

       b(n-2) = 4.d0*(1.d0-ep)
       b(n-1) = 1.d0-ep
       b(n) = 0d0

       c(n-2) = 6.d0*(1.d0+ep)
       c(n-1) = 2.d0*(1.d0+ep)
       c(n) = 1.d0

       d(n-2) = 4.d0*(1.d0-ep)
       d(n-1) = 1.d0-ep
       d(n) = 0d0

       e(n-2) = 1.d0+ep
       e(n-1) = 0d0
       e(n) = 0d0

       f(n-2) = 0d0 
       f(n-1) = 0d0
       f(n) = 0d0

       !     step one.
       beta(1) = d(1)/c(1)
       delta(2) = b(2)
       w(1) = c(1)
       pi(1) = f(1)/w(1)
       ap(1) = 0d0
       ap(2) = 0d0
       ap(3) = a(3)
       w(2) = c(2)-delta(2)*beta(1)
       gam(1) = e(1)/c(1)
       beta(2) = (d(2)-delta(2)*gam(1))/w(2)
       gam(2) = (e(2)-pi(1)*delta(2))/w(2)
       pi(2) = f(2)/w(2)
       delta(3) = (b(3)-ap(3)*beta(1))
       w(3) = c(3)-delta(3)*beta(2)-ap(3)*gam(1)
       beta(3) = (d(3)-ap(3)*pi(1)-delta(3)*gam(2))/w(3)
       gam(3) = (e(3)-delta(3)*pi(2))/w(3)
       pi(3) = f(3)/w(3)

       !     step two
       do 20 i=4,n
          ap(i) = a(i)-z(i)*beta(i-3)
          delta(i) = b(i)-ap(i)*beta(i-2)-z(i)*gam(i-3)
          w(i) = c(i)-ap(i)*gam(i-2)-delta(i)*beta(i-1) &
               -z(i)*pi(i-3)
          beta(i) = (d(i)-ap(i)*pi(i-2)-delta(i)*gam(i-1))/w(i)
          gam(i) = (e(i)-delta(i)*pi(i-1))/w(i)
          pi(i) = f(i)/w(i)
20     enddo
100 endif

    !     step three
    h(1) = bb(1)/w(1)
    h(2) = (bb(2)-delta(2)*h(1))/w(2)
    h(3) = (bb(3)-delta(3)*h(2)-ap(3)*h(1))/w(3)
    do 30 i=4,n
       h(i) = (bb(i)-delta(i)*h(i-1)-ap(i)*h(i-2) &
            -z(i)*h(i-3))/w(i)
30  enddo

    !     step four
    xans(n) = h(n)
    xans(n-1) = h(n-1)-beta(n-1)*xans(n)
    xans(n-2) = h(n-2)-beta(n-2)*xans(n-1)-gam(n-2)*xans(n)
    do 40 i=n-3,1,-1
       xans(i) = h(i)-beta(i)*xans(i+1)-gam(i)*xans(i+2) &
            -pi(i)*xans(i+3)
40  enddo
  end subroutine invlow


end module filter_orog

