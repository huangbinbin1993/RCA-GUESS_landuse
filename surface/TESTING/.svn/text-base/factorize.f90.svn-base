subroutine localfactorize(npr,nig,njg,nkg,np1,np2,np3)
  use decomp
  implicit none
  integer,intent(in)::npr,nig,njg,nkg
  integer,intent(inout)::np1,np2,np3
  double precision nor,norloc
  integer p1,p2,p3

  call factorize(npr,nig,njg,nkg,np1,np2,np3)

!!$  nor = nig+2*njg+nkg
!!$  do p3=1,npr
!!$     do p1=1,npr
!!$        if(mod(npr,p1*p3).eq.0)then
!!$           p2 = npr/(p1*p3)
!!$           norloc = abs(nig/real(p1,kind(norloc)) - &
!!$                njg/real(p2,kind(norloc))  ) + &
!!$                abs(njg/real(p2,kind(norloc)) - &
!!$                nkg/real(p3,kind(norloc)) )
!!$           if(norloc.lt.nor)then
!!$              np1 = p1
!!$              np2 = p2
!!$              np3 = p3
!!$              nor = norloc
!!$           endif
!!$        endif
!!$     enddo
!!$  enddo
end subroutine localfactorize
