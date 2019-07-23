module modddr
  use mod_grib,only:realgribkind

  implicit none
  private
  integer,parameter::jpnslfx=120
  integer,parameter::jpnmlfx=40
!  integer,parameter::jplevx=60
! ecearth has 62 levels
  integer,parameter::jplevx=62
  !  integer,parameter::mddrsz=21 + 5*jpnslfx + 2*jpnmlfx + jplevx

  type,public::ddr
     !time section
     integer::ndtvhl  !kyear*10000+kmonth*100+kday
     integer::nscvhl  !khour*10000
     integer::ndtbhl  !(year*100+month)*100+day
     integer::nscbhl  !(hour*100+minut)*100+00
     integer::nflshl
     !grid section
     integer::nlonhl
     integer::nlathl
     integer::nlevhl
     !level section
     integer::nltphl
     integer::npplhl
     integer::nrflhl
     integer::nlpthl(jplevx)
     !multi level fields
     integer::nmlfhl
     integer::nwmmhl(jpnmlfx) !par ml
     integer::nmpthl(jpnmlfx) !lev
     !sinle level fields
     integer::nslfhl
     integer::nwmshl(jpnslfx) !par sl
     integer::nspthl(jpnslfx) !lev sl
     integer::nslthl(jpnslfx) !typ sl
     !geometric info
     real(kind=realgribkind)::aplohl
     real(kind=realgribkind)::aplahl
     real(kind=realgribkind)::aweshl
     real(kind=realgribkind)::aeashl
     real(kind=realgribkind)::alalhl
     real(kind=realgribkind)::alafhl
     real(kind=realgribkind)::dlonhl
     real(kind=realgribkind)::dlathl
     !level
     real(kind=realgribkind)::alevhl(jpnslfx,2)
  end type ddr
!  type(ddr),public::myddr

  public reset_ddr,commddr

contains
  subroutine commddr(rank,ddr_in)
    use decomp
    implicit none
    integer,intent(in)::rank
    type(ddr),intent(inout)::ddr_in
#ifdef MPI_SRC
#include"mpif.h"
    integer::ibuf(13+jplevx+2*jpnmlfx+3*jpnslfx)
    real(kind=realgribkind)::rbuf(8+2*jpnslfx)
    integer::ierr,msgSize
    msgSize = size(ibuf)

    ibuf(1) = ddr_in%ndtvhl 
    ibuf(2) = ddr_in%nscvhl 
    ibuf(3) = ddr_in%ndtbhl 
    ibuf(4) = ddr_in%nscbhl 
    ibuf(5) = ddr_in%nflshl 
    ibuf(6) = ddr_in%nlonhl  
    ibuf(7) = ddr_in%nlathl 
    ibuf(8) = ddr_in%nlevhl 
    ibuf(9) = ddr_in%nltphl 
    ibuf(10) = ddr_in%npplhl 
    ibuf(11) = ddr_in%nrflhl 
    ibuf(12) = ddr_in%nmlfhl 
    ibuf(13) = ddr_in%nslfhl 
    ibuf(14:14+jplevx-1) = ddr_in%nlpthl 
    ibuf(14+jplevx:14+jplevx+jpnmlfx-1) = ddr_in%nwmmhl 
    ibuf(14+jplevx+jpnmlfx:14+jplevx+2*jpnmlfx-1) = ddr_in%nmpthl 
    ibuf(14+jplevx+2*jpnmlfx:14+jplevx+2*jpnmlfx+jpnslfx-1) = ddr_in%nwmshl 
    ibuf(14+jplevx+2*jpnmlfx+jpnslfx:14+jplevx+2*jpnmlfx+2*jpnslfx-1) = ddr_in%nspthl 
    ibuf(14+jplevx+2*jpnmlfx+2*jpnslfx:14+jplevx+2*jpnmlfx+3*jpnslfx-1) = ddr_in%nslthl 


    call MPI_bcast(ibuf,msgSize,MPI_INTEGER,rank,localComm,ierr)

    ddr_in%ndtvhl = ibuf(1)   
    ddr_in%nscvhl = ibuf(2)   
    ddr_in%ndtbhl = ibuf(3)   
    ddr_in%nscbhl = ibuf(4)   
    ddr_in%nflshl = ibuf(5)   
    ddr_in%nlonhl = ibuf(6)    
    ddr_in%nlathl = ibuf(7)   
    ddr_in%nlevhl = ibuf(8)   
    ddr_in%nltphl = ibuf(9)   
    ddr_in%npplhl = ibuf(10)  
    ddr_in%nrflhl = ibuf(11)  
    ddr_in%nmlfhl = ibuf(12)  
    ddr_in%nslfhl = ibuf(13)  
    ddr_in%nlpthl = ibuf(14:14+jplevx-1) 
    ddr_in%nwmmhl = ibuf(14+jplevx:14+jplevx+jpnmlfx-1) 
    ddr_in%nmpthl = ibuf(14+jplevx+jpnmlfx:14+jplevx+2*jpnmlfx-1) 
    ddr_in%nwmshl = ibuf(14+jplevx+2*jpnmlfx:14+jplevx+2*jpnmlfx+jpnslfx-1) 
    ddr_in%nspthl = ibuf(14+jplevx+2*jpnmlfx+jpnslfx:14+jplevx+2*jpnmlfx+2*jpnslfx-1) 
    ddr_in%nslthl = ibuf(14+jplevx+2*jpnmlfx+2*jpnslfx:14+jplevx+2*jpnmlfx+3*jpnslfx-1)  

    msgSize = size(rbuf)
    rbuf(1)=ddr_in%aplohl
    rbuf(2)=ddr_in%aplahl
    rbuf(3)=ddr_in%aweshl
    rbuf(4)=ddr_in%aeashl
    rbuf(5)=ddr_in%alalhl
    rbuf(6)=ddr_in%alafhl
    rbuf(7)=ddr_in%dlonhl
    rbuf(8)=ddr_in%dlathl
    rbuf(9:9+jpnslfx-1)=ddr_in%alevhl(:,1)
    rbuf(9+jpnslfx:9+2*jpnslfx-1)=ddr_in%alevhl(:,2)

    if(realgribkind==4)then
       call MPI_bcast(rbuf,msgSize,MPI_REAL,rank,localComm,ierr)
    elseif(realgribkind==8)then
       call MPI_bcast(rbuf,msgSize,MPI_DOUBLE_PRECISION,rank,localComm,ierr)
    else
       stop 'moddr'
    endif
    ddr_in%aplohl = rbuf(1)
    ddr_in%aplahl = rbuf(2)
    ddr_in%aweshl = rbuf(3)
    ddr_in%aeashl = rbuf(4)
    ddr_in%alalhl = rbuf(5)
    ddr_in%alafhl = rbuf(6)
    ddr_in%dlonhl = rbuf(7)
    ddr_in%dlathl = rbuf(8)
    ddr_in%alevhl(:,1)=rbuf(9:9+jpnslfx-1)
    ddr_in%alevhl(:,2)=rbuf(9+jpnslfx:9+2*jpnslfx-1)
#endif    
    return
  end subroutine commddr

  subroutine reset_ddr(ddr_in)
    implicit none
    integer,parameter::unini=0
    real(kind=realgribkind),parameter::uninr=0.0_realgribkind
    type(ddr),intent(inout)::ddr_in
    ddr_in%ndtvhl = unini
    ddr_in%nscvhl = unini
    ddr_in%ndtbhl = unini
    ddr_in%nscbhl = unini
    ddr_in%nflshl = unini
    ddr_in%nlonhl = unini 
    ddr_in%nlathl = unini
    ddr_in%nlevhl = unini
    ddr_in%nltphl = unini
    ddr_in%npplhl = unini
    ddr_in%nrflhl = unini
    ddr_in%nlpthl = unini
    ddr_in%alevhl = uninr
    ddr_in%nmlfhl = unini
    ddr_in%nwmmhl = unini
    ddr_in%nmpthl = unini
    ddr_in%nslfhl = unini
    ddr_in%nwmshl = unini
    ddr_in%nspthl = unini
    ddr_in%nslthl = unini
    ddr_in%aplohl = uninr
    ddr_in%aplahl = uninr 
    ddr_in%aweshl = uninr
    ddr_in%aeashl = uninr
    ddr_in%alalhl = uninr 
    ddr_in%alafhl = uninr
    ddr_in%dlonhl = uninr
    ddr_in%dlathl = uninr
  end subroutine reset_ddr

end module modddr
