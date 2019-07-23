module timers
  use decomp
  implicit none
  private
  character(len=10),save::cc(1000)
  integer,save::itim
  integer,save::ntim
  real(kind=realkind),save::t1(1000)
  real(kind=realkind),save::t2(1000)

  public timer,second
contains

  subroutine second(time)
  implicit none
#ifdef MPI_SRC
#include "mpif.h"
#endif

  real(kind=realkind) time
  real(kind=8):: sdtime
  integer first
  save first,sdtime
  if ( first /= 12789 ) then
#ifdef MPI_SRC
     sdtime=mpi_wtime()
#else
     call csecond(sdtime)
#endif
     first=12789
  endif
#ifdef MPI_SRC
  time=mpi_wtime()-sdtime
#else
  time = sdtime
  call csecond(sdtime)
  time = sdtime-time
#endif
  return
end subroutine second


  subroutine timer(ktype,char)
    use decomp
    implicit none
    integer::ktype
    character(len=*)::char
    real(kind=realkind):: sec
    integer:: ktim

    if(ktype==0)then
       itim=0
       ntim=0
       return
    endif

    !     start timing of 'char'
    if(ktype==1)then
       if(ntim>0) then
          do ktim=1,ntim
             if(cc(ktim)==char) then
                itim=ktim
                goto 100
             end if
          end do
       end if

       ntim=ntim+1
       itim=ntim
       t2(itim) = 0._realkind
       cc(itim)=char

100    continue

       call second(sec)
       t1(itim)=sec
print*,'itim t1(itim) ',itim,t1(itim),char
       return
    endif


    !     stop timing  of 'char'

    if(ktype==2)then
       if(ntim>0) then
          do ktim=1,ntim
             if(cc(ktim)==char) then
                itim=ktim
                goto 200
             end if
          end do
       end if

       ntim=ntim+1
       itim=ntim
       t2(itim) = 0._realkind
       cc(itim)=char

200    continue

       call second(sec)
       t2(itim)=t2(itim)+sec-t1(itim)
print*,'itim t2(itim) ',itim,t2(itim),char
       return
    endif
  end subroutine timer
end module timers
