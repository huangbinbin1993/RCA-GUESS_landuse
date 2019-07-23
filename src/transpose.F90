module transpose
  use comtrix
  use decomp

  implicit none
  private

  !variables associated with the fft-decomposition
  integer,public,save::nslab_fft  
  integer,allocatable,dimension(:,:),public,save::pe_fft  !(klat_global, nproc) 
  integer,allocatable,dimension(:,:),public,save::kpro_fft!(klat_global,klev_global)
  integer,allocatable,dimension(:,:),public,save::krec_fft!(klat_global,klev_global
  integer,allocatable,dimension(:),public,save::jmin_fft  !(klev_global) 
  integer,allocatable,dimension(:),public,save::jlen_fft  !(klev_global)  
  integer,allocatable,dimension(:),public,save::jlev_fft  !(klev_global)

  !variables associated with the tri-decomposition
  integer,public,save::nslab_tri  
  integer,allocatable,dimension(:),public,save::imin_tri  !(NSLAB_TRI_MAX) NSLAB_TRI_MAX=KLEV_GLOBAL_MAX+2
  integer,allocatable,dimension(:),public,save::imax_tri  
  integer,allocatable,dimension(:),public,save::jlev_tri  
  integer,allocatable,dimension(:,:),public,save::kpro_tri 
  integer,allocatable,dimension(:,:),public,save::krec_tri 
  integer,allocatable,dimension(:),public,save::pe_tri     
  integer,allocatable,dimension(:),public,save::nproc_tri  
  integer,allocatable,dimension(:,:),public,save::jlen_tri 


  public twod_to_fft,fft_to_tri,tri_to_fft,fft_to_twod,decompose_hh,allocate_transpose
contains

  subroutine allocate_transpose(klon_global, klat_global, klev_global,nproc)
    implicit none
    integer,intent(in)::klon_global, klat_global, klev_global,nproc
    allocate(kpro_fft(klat_global, klev_global))
    allocate(krec_fft(klat_global, klev_global))
    allocate(jmin_fft(klev_global),jlen_fft(klev_global), jlev_fft(klev_global))
    allocate(pe_fft(klat_global, nproc))

    allocate(kpro_tri(klon_global, klev_global)) 
    allocate(krec_tri(klon_global, klev_global))
    allocate(pe_tri(klev_global))
    allocate(nproc_tri(klev_global))
    allocate(jlen_tri(klev_global, nproc))
    allocate(imin_tri(klev_global),imax_tri(klev_global),jlev_tri(klev_global))
  end subroutine allocate_transpose


  subroutine decompose_hh(klev,kpbpts)

    !     PURPOSE:
    !     Initializes the rquired arrays needed when transposing data
    !     PHASE1 --> PHASE2 -->  PHASE3 --> PHASE2 --> PHASE1  
    !     where
    !     PHASE1 = normal xy-decomposition
    !     PHASE2 = phase used during Fourier transformation or
    !     filtering in x-direction 
    !     PHASE3 = phase used during Gaussian elimination or
    !     filtering in y-direction 
    !     DESCRIPTION:
    !     The arrays in DECOMP.inc defines the layout of data in
    !     different phases
    !     AUTHOR:
    !     From the SCIFI-model modified by
    !     KALLE EEROLA   HIRLAM    1998

    implicit none

    integer,intent(in):: klev,kpbpts
    integer:: jpe,jmin,jmax
    integer:: kdist,krem,kproc,krec,maxrec,k,j,i

!    klev_hh = klev
    triwgt  = -1

    !      fft decomposition
    kdist = klat_global*klev/nproc
    krem  = klat_global*klev-nproc*kdist

    kproc = 0
    krec  = 0
    do k=1,klev
       do j=1,klat_global
          krec = krec + 1
          if( krem==0 .or. kproc>=krem ) then
             maxrec = kdist
          else
             maxrec = kdist + 1
          endif
          kpro_fft(j,k) = kproc
          krec_fft(j,k) = krec
          if( krec==maxrec ) then
             kproc = kproc + 1
             krec  = 0
          endif
       enddo
    enddo


    nslab_fft = 0
    do k=1,klev
       do j=1,klat_global
          if( kpro_fft(j,k)==mype ) then
             if( nslab_fft==0 ) then
                nslab_fft = 1
                jmin_fft(nslab_fft) = j
                jlen_fft(nslab_fft) = 1
                jlev_fft(nslab_fft) = k
             else if( jlev_fft(nslab_fft)==k ) then
                jlen_fft(nslab_fft) = jlen_fft(nslab_fft) + 1
             else
                nslab_fft = nslab_fft + 1
                jmin_fft(nslab_fft) = j
                jlen_fft(nslab_fft) = 1
                jlev_fft(nslab_fft) = k
             endif
          endif
       enddo
    enddo

    do jpe=1,nprocy
       if(jpe==1) then
          jmin=jdatastart_ar2d(jpe)
       else
          jmin=jdatastart_ar2d(jpe)+1
       endif
       if(jpe==nprocy) then
          jmax=klat_global
       else
          jmax=jdatastart_ar2d(jpe+1)
       endif
       do j=jmin,jmax
          do k=1,nprocx
             pe_fft(j,k)=(jpe-1)*nprocx+k-1  
          enddo
       enddo
    enddo



    !      tri-diagonal decomposition
    kdist =(klon_global-2*kpbpts-3)*klev/nproc
    krem  =(klon_global-2*kpbpts-3)*klev-nproc*kdist

    kproc = 0
    krec  = 0
    do k=1,klev
       do i=kpbpts+2,klon_global-kpbpts-2
          krec = krec + 1
          if( krem==0 .or. kproc>=krem ) then
             maxrec = kdist
          else
             maxrec = kdist + 1
          endif
          kpro_tri(i,k) = kproc
          krec_tri(i,k) = krec
          if( krec==maxrec ) then
             kproc = kproc + 1
             krec  = 0
          endif
       enddo
    enddo


    nslab_tri = 0
    do k=1,klev
       do i=kpbpts+2,klon_global-kpbpts-2
          if( kpro_tri(i,k)==mype ) then
             if( nslab_tri==0 ) then
                nslab_tri = 1
                imin_tri(nslab_tri) = i
                imax_tri(nslab_tri) = i
                jlev_tri(nslab_tri) = k
             else if( jlev_tri(nslab_tri)==k ) then
                imax_tri(nslab_tri) = i
             else
                nslab_tri = nslab_tri + 1
                imin_tri(nslab_tri) = i
                imax_tri(nslab_tri) = i
                jlev_tri(nslab_tri) = k
             endif
          endif
       enddo
    enddo
    !    print *,mype,'nslab_tri=',nslab_tri,nprocx,nprocy
    !    call stop_program('')
    do k=1,klev
       pe_tri(k)=kpro_fft(1,k)
       nproc_tri(k)=1
       jlen_tri(k,nproc_tri(k))=0
       kproc=pe_tri(k)
       do j=1,klat_global
          if(kpro_fft(j,k)==kproc) then
             jlen_tri(k,nproc_tri(k))=jlen_tri(k,nproc_tri(k)) + 1
          else
             nproc_tri(k)=nproc_tri(k)+1
             kproc=kproc+1
             jlen_tri(k,nproc_tri(k))=1
          endif
       enddo
    enddo
    return
  end subroutine decompose_hh



  subroutine twod_to_fft(div, klon, klat, klev ,div_fft )
    !     AUTHOR:
    !     KALLE EEROLA   HIRLAM    1998


    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(in):: div(klon,klat,klev)
    real(kind=realkind),intent(out)::div_fft(klon_global,*) !(klon_global,klat_global*klev_global)

    !     declaration of local variables
    integer:: k,j,rempro,remrec,jmin,jmax,imin,nvalues,imax
    integer:: error,ipe,fft_tag,iminr,nvaluesr,idatastartr
    integer:: jslab,jlen,jlev
    integer:: np, t_send, t_recv

    call jlimits(klon,klat,imin,imax,jmin,jmax)

    nvalues = klon - 2

    if( atleft ) then
       nvalues = nvalues + 1
    endif

    if( atright ) then
       nvalues = nvalues + 1
    endif

    if (nproc == 1 ) then
       do k=1,klev
          do j=jmin,jmax
             remrec = krec_fft(jdatastart+j-1,k)
             !copy div_fft(1:klon_global,?,?)
             call scopy(nvalues,div(imin,j,k),1,div_fft(idatastart+imin-1,remrec),1)
          enddo
       enddo
       return
    endif

#ifdef MPI_SRC
    call fft2dsend(div,klon,klat,klev,mype,div_fft)
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
    do np=0,nproc-1
       t_send=mod(nproc+np-mype,nproc)
       t_recv=t_send
       if (t_send/=mype) then
          if (mype<t_send) then
             call fft2dsend(div,klon,klat,klev,t_send,div_fft)
             call fft2drecv(div,klon,klat,klev,t_recv,div_fft)
          else
             call fft2drecv(div,klon,klat,klev,t_recv,div_fft)
             call fft2dsend(div,klon,klat,klev,t_send,div_fft)
          endif
       endif
    enddo
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
#endif
    return
  end subroutine twod_to_fft


  subroutine fft_to_tri(div_fft,kpbpts, klev,div_tri)

    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer,intent(in):: kpbpts,klev
    real(kind=realkind),intent(in)::div_fft(klon_global,*) 
    real(kind=realkind),intent(out)::div_tri(klat_global,*)
    integer:: irec,jslab,i,jmin,jlen,jlev,rempro,remrec
    integer:: error,kproc,tri_tag,nrempro
    integer:: islab,imin,imax
    integer:: ncol,isend
    integer:: j, ii,jj,kk, ndx, lenw
    integer:: np, t_send, t_recv
    !    IN CASE OF SHARED MEMORY ONLY COPY IS NEEDED

    if(nproc==1) then
       irec = 1
       do jslab=1,nslab_fft
          jmin = jmin_fft(jslab)
          jlen = jlen_fft(jslab)
          jlev = jlev_fft(jslab)
          do i=kpbpts+2,klon_global-kpbpts-2
             rempro = kpro_tri(i,jlev)
             remrec = krec_tri(i,jlev)
             call scopy(jlen, div_fft(i,irec),klon_global, div_tri(jmin,remrec),1)
          enddo
          irec = irec + jlen
       enddo
       return
    endif

    call fft_tri_send(div_fft,kpbpts,klev,mype,div_tri)
#ifdef MPI_SRC
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
#endif
    do np=0,nproc-1
       t_send=mod(nproc+np-mype,nproc)
       t_recv=t_send
       if (t_send/=mype) then
          if (mype<t_send) then
             call fft_tri_send(div_fft,kpbpts,klev,t_send,div_tri)
             call fft_tri_recv(kpbpts,klev,t_recv,div_tri)
          else
             call fft_tri_recv(kpbpts,klev,t_recv,div_tri)
             call fft_tri_send(div_fft,kpbpts,klev,t_send,div_tri)
          endif
       endif
    enddo
#ifdef MPI_SRC
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
#endif
    return
  end subroutine fft_to_tri


  subroutine tri_to_fft(div_tri,kpbpts, klev,div_fft )

    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    real(kind=realkind),intent(out)::div_fft(klon_global,*) 
    real(kind=realkind),intent(in)::div_tri(klat_global,*)
    integer,intent(in):: kpbpts,klev
    integer irec,jslab,i,jmin,jlen,jlev,rempro,remrec,error
    integer np, t_send, t_recv
    if(nproc == 1) then
       irec = 1
       do jslab=1,nslab_fft
          jmin = jmin_fft(jslab)
          jlen = jlen_fft(jslab)
          jlev = jlev_fft(jslab)
          do i=kpbpts+2,klon_global-kpbpts-2
             rempro = kpro_tri(i,jlev)
             remrec = krec_tri(i,jlev)
             call scopy(jlen,div_tri(jmin,remrec),1,div_fft(i,irec),klon_global)
          enddo
          irec = irec + jlen
       enddo
       return
    endif

    call tri_send(div_tri,kpbpts,klev,mype,div_fft)
#ifdef MPI_SRC
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
#endif
    do np=0,nproc-1
       t_send=mod(nproc+np-mype,nproc)
       t_recv=t_send
       if (t_send/=mype) then
          if (mype<t_send) then
             call tri_send(div_tri,kpbpts,klev,t_send,div_fft)
             call tri_recv(kpbpts,klev,t_recv,div_fft)
          else
             call tri_recv(kpbpts,klev,t_recv,div_fft)
             call tri_send(div_tri,kpbpts,klev,t_send,div_fft)
          endif
       endif
    enddo
#ifdef MPI_SRC
#ifdef USE_BARRIER
    call mpi_barrier(localcomm,error)
#endif
#endif

    return
  end subroutine tri_to_fft


  subroutine fft_to_twod(div_fft, klon, klat, klev,div )


    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    real(kind=realkind),intent(in)::div_fft(klon_global,*)
    integer,intent(in):: klon,klat,klev
    real(kind=realkind),intent(out)::div(klon,klat,klev)

    ! DECLARATION OF LOCAL VARIABLES

    integer k,j,rempro,remrec,jmin,jmax,imin,nvalues,imax
    integer error,ipe,fft_tag,iminr,nvaluesr,idatastartr
    integer np, t_send, t_recv

    call jlimits(klon,klat,imin,imax,jmin,jmax)
    nvalues = klon - 2
    if( atleft ) then
       nvalues = nvalues + 1
    endif
    if( atright ) then
       nvalues = nvalues + 1
    endif

    if(nproc == 1) then
       do k=1,klev
          do j=jmin,jmax
             rempro = kpro_fft(jdatastart+j-1,k)
             remrec = krec_fft(jdatastart+j-1,k)
             call scopy(nvalues,div_fft(idatastart+imin-1,remrec),1,div(imin,j,k),1)
          enddo
       enddo
       return
    endif


    call f2dsend(div,klon,klat,klev,mype,div_fft)
#ifdef MPI_SRC
#ifdef USE_BARRIER
    call MPI_barrier(localComm,error)
#endif
#endif
    do np=0,nproc-1
       t_send=mod(nproc+np-mype,nproc)
       t_recv=t_send
       if (t_send/=mype) then
          if (mype<t_send) then
             call f2dsend(div,klon,klat,klev,t_send,div_fft)
             call f2drecv(div,klon,klat,klev,t_recv)
          else
             call f2drecv(div,klon,klat,klev,t_recv)
             call f2dsend(div,klon,klat,klev,t_send,div_fft)
          endif
       endif
    enddo

#ifdef MPI_SRC
#ifdef USE_BARRIER
    call MPI_barrier(localComm,error)
#endif
#endif
    RETURN
  END subroutine fft_to_twod




  subroutine fft2dgath(div,klon,klat,klev,rb,rblen,ib,iblen,destpro,div_fft)

    IMPLICIT NONE
    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv

    integer kpbpts, rblen, iblen, destpro
    real(kind=realkind) rb(rblen)
    !     The IB is Index Buffer, containing index information
    !     for databuffer RB
    integer ib(3,iblen)
    integer irec,jslab,i,jmin,jlen,jlev,rempro,remrec

    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax,jmax,iminr
    integer ncol,isend,fft_tag
    integer j, k,ii,jj,kk, ndx, lenw,sumlenw
    integer rbc,ibc
    integer nvalues

    sumlenw=0
    call jlimits(klon,klat,imin,imax,jmin,jmax)
    nvalues = klon - 2
    if( atleft ) then
       nvalues = nvalues + 1
    endif
    if( atright ) then
       nvalues = nvalues + 1
    endif
    irec = 1
    rbc = 0
    ibc = 1
    if (destpro==mype) then
       do k=1,klev
          do j=jmin,jmax
             rempro = kpro_fft(jdatastart+j-1,k)
             remrec = krec_fft(jdatastart+j-1,k)
             if (rempro==mype) then
                call scopy(nvalues,div(imin,j,k),1, &
                     div_fft(idatastart+imin-1,remrec),1)
             endif
          enddo
       enddo
    else
       do k=1,klev
          do j=jmin,jmax
             rempro = kpro_fft(jdatastart+j-1,k)
             if (rempro==destpro) then
                fft_tag=1910+idatastart+imin-1+ &
                     ((jdatastart+j-2)+(k-1)*klat_global)*klon_global
                if (rempro/=mype) then
                   do i=imin,imin+nvalues-1
                      rbc = rbc+1
                      rb(rbc)=div(i,j,k)
                   enddo
                endif
             endif
          enddo
       enddo
    endif

    rblen = rbc
    iblen = ibc

    return
  end subroutine fft2dgath


  subroutine fft2dscat(div,klon,klat,klev,rb,rblen,ib,iblen,destpro,div_fft)

    implicit none

    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv
    integer kpbpts, rblen, iblen, destpro

    real(kind=realkind) rb(rblen)
    integer ib(3,iblen)
    integer irec,jslab,jmin,jlen,jlev,rempro,remrec
    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax,idatastartr,iminr
    integer ncol,isend,fft_tag,nvaluesr,ipe
    integer i,j,k, ii,jj,kk, ndx, lenw,sumlenw
    integer rbc,ibc

    rbc = 0
    sumlenw=0
    do jslab=1,nslab_fft
       jmin = jmin_fft(jslab)
       jlen = jlen_fft(jslab)
       jlev = jlev_fft(jslab)
       do j=jmin,jmin+jlen-1
          do ipe=1,nprocx
             idatastartr=idatastart_ar2d(ipe)
             nvaluesr=klon_ar2d(ipe)-2
             if(ipe==1) then
                iminr = 1
                nvaluesr = nvaluesr + 1
             else
                iminr = 2
             endif
             if(ipe==nprocx) then
                nvaluesr = nvaluesr + 1
             endif
             if (pe_fft(j,ipe)==destpro) then
                fft_tag=1910+idatastartr+iminr-1+ &
                     ((j-1)+(jlev-1)*klat_global)*klon_global
                sumlenw = sumlenw+nvaluesr
                if (iblen /= 0) then
                   ibc=1
                   remrec = krec_fft(j,jlev)
                   do i=idatastartr+iminr-1,idatastartr+iminr-1+nvaluesr-1
                      rbc = rbc+1
                      div_fft(i,remrec) = rb(rbc)
                   enddo
                endif
             endif
          enddo
       enddo
    enddo
    rblen = sumlenw

    return
  end subroutine fft2dscat


  subroutine fft2dsend(div,klon,klat,klev,t_send,div_fft)

    implicit none
    integer,intent(in)::klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv
    integer kpbpts,t_send
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer rbuflen, ibuflen
    integer error
    integer,parameter::maxibuf=100000,maxrbuf=3000000
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)

    rbuflen=maxrbuf
    ibuflen=maxibuf-1
    call fft2dgath(div,klon,klat,klev,rbuf,rbuflen,ibuf,ibuflen,t_send,div_fft)
    ibuflen=1
    if (ibuflen>=maxibuf) then
       write(*,*) 'pe',mype,'fatal in fft2dsend(), maxibuf too small',maxibuf,ibuflen
       call stop_program('')
    endif
    if (rbuflen>=maxrbuf) then
       write(*,*) 'pe',mype,'fatal in fft2dsend(), maxrbuf too small',maxrbuf,rbuflen
       call stop_program('')
    endif

    if (t_send/= mype) then
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_send(rbuf,rbuflen,REALTYPE,t_send,2200+t_send,localComm,error)
#endif
       endif
    endif
    return
  end subroutine fft2dsend


  subroutine fft2drecv(div,klon,klat,klev,t_recv,div_fft)

    implicit none
    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv
#ifdef MPI_SRC
#include"mpif.h"
    integer status(mpi_status_size)
#endif
    integer rbuflen, ibuflen
    integer error
    integer,parameter::maxibuf=100000,maxrbuf=3000000
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)
    

    if (t_recv/=mype) then
       call fft2dscat(div,klon,klat,klev,rbuf,rbuflen,ibuf,0,t_recv,div_fft)
       ibuflen=1
       if(rbuflen/=0) then
#ifdef MPI_SRC
          call mpi_recv(rbuf,rbuflen,REALTYPE,t_recv,2200+mype,localComm,status,error)
#endif
          call fft2dscat(div,klon,klat,klev,rbuf,rbuflen,ibuf,ibuflen,t_recv,div_fft)
       endif
    endif
    return
  end subroutine fft2drecv




  subroutine fftgath(div_fft,kpbpts,klev,rb,rblen,ib,iblen,destpro,div_tri)

    implicit none

    real(kind=realkind) div_fft(klon_global,*)
    real(kind=realkind) div_tri(klat_global,*)
    INTEGER KPBPTS,KLEV, rblen, iblen, destpro
    REAL(KIND=REALKIND) RB(RBLEN)
    !     The IB is Index Buffer, containing index information
    !     for databuffer RB
    INTEGER IB(3,IBLEN)
    INTEGER IREC,JSLAB,I,JMIN,JLEN,JLEV,REMPRO,REMREC

    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax
    integer ncol,isend
    integer j, ii,jj,kk, ndx, lenw
    integer rbc,ibc

    irec = 1
    rbc = 0
    ibc = 1
    if (destpro==mype) then
       irec = 1
       do jslab=1,nslab_fft
          jmin = jmin_fft(jslab)
          jlen = jlen_fft(jslab)
          jlev = jlev_fft(jslab)
          do i=kpbpts+2,klon_global-kpbpts-2
             rempro = kpro_tri(i,jlev)
             remrec = krec_tri(i,jlev)
             if (rempro==mype) then
                call scopy(jlen,div_fft(i,irec),klon_global,div_tri(jmin,remrec),1)
             endif
          enddo
          irec = irec + jlen
       enddo
    else

       irec=1
       do jslab=1,nslab_fft
          jmin = jmin_fft(jslab)
          jlen = jlen_fft(jslab)
          jlev = jlev_fft(jslab)
          isend=kpbpts+2
          ncol=1
          do i=kpbpts+2,klon_global-kpbpts-2
             rempro = kpro_tri(i,jlev)
             if (i == klon_global-kpbpts-2) then
                nrempro = -1
             else
                nrempro = kpro_tri(i+1,jlev)
             endif
             if(nrempro /= rempro) then
                if (rempro == destpro) then
                   tri_tag=2200+isend+(jmin-1+(jlev-1)* klat_global)*klev_global
                   if (rempro /= mype) then
                      do kk=irec,irec+jlen-1
                         do ii=isend,i
                            rbc = rbc + 1
                            rb(rbc) = div_fft(ii,kk)
                         enddo
                      enddo
                   endif
                endif
                isend=i+1
                ncol=1
             else
                ncol=ncol+1
             endif
          enddo
          irec = irec + jlen
       enddo
    endif

    rblen = rbc
    iblen = ibc
    return
  end subroutine fftgath


  subroutine fftscat(kpbpts,klev,rb,rblen,ib,iblen,destpro,div_tri)

    implicit none

    real(kind=realkind) div_tri(klat_global,*)

    INTEGER KPBPTS,KLEV, rblen, iblen, destpro
    REAL(KIND=REALKIND) RB(RBLEN)
    INTEGER IB(3,IBLEN)
    INTEGER IREC,JSLAB,I,JMIN,JLEN,JLEV,REMPRO,REMREC

    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax
    integer ncol,isend
    integer j, ii,jj,kk, ndx, lenw,sumlenw
    integer rbc,ibc

    sumlenw=0
    rbc = 0
    do islab=1,nslab_tri
       imin = imin_tri(islab)
       imax = imax_tri(islab)
       jlev = jlev_tri(islab)
       kproc = nproc_tri(jlev)
       remrec=krec_tri(imin,jlev)
       jmin = 1
       do j=1,kproc
          tri_tag=2200+imin+(jmin-1+(jlev-1)*klat_global)*klev_global
          lenw = jlen_tri(jlev,j) * (imax-imin+1)
          if (destpro==pe_tri(jlev)+j-1) then
             sumlenw = sumlenw + lenw
             if (iblen /= 0) then
                ibc=1
                do jj=jmin,jmin+jlen_tri(jlev,j)-1
                   do kk=remrec,remrec+imax-imin
                      rbc = rbc + 1
                      div_tri(jj,kk) = rb(rbc)
                   enddo
                enddo
             endif
          endif
          jmin = jmin + jlen_tri(jlev,j)
       enddo
    enddo
    rblen = sumlenw

    return
  end subroutine fftscat


  subroutine fft_tri_send(div_fft,kpbpts,klev,t_send,div_tri)

    implicit none
    real(kind=realkind) div_fft(klon_global,*)
    real(kind=realkind) div_tri(klat_global,*)
    integer kpbpts,klev,t_send
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)

    rbuflen=maxrbuf
    ibuflen=maxibuf-1
    call fftgath(div_fft,kpbpts,klev,rbuf,rbuflen,ibuf,ibuflen,t_send,div_tri)
    if (ibuflen>=maxibuf) then
       write(*,*) 'pe',mype,'FATAL in FFT_TRI_SEND(), MAXIBUF TOO SMALL', maxibuf,ibuflen
    endif
    if (rbuflen>=maxrbuf) then
       write(*,*) 'pe',mype,'FATAL in FFT_TRI_SEND(), MAXRBUF TOO SMALL ', maxrbuf,rbuflen
    endif
    ibuflen=1
    if (t_send /= mype) then
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_send(rbuf,rbuflen,REALTYPE,t_send,2201+t_send,localComm,error)
#endif
       endif
    endif
    return
  end subroutine fft_tri_send


  subroutine fft_tri_recv(kpbpts,klev,t_recv,div_tri)

    implicit none
    real(kind=realkind) div_tri(klat_global,*)
    integer kpbpts,klev,t_recv
    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)
#ifdef MPI_SRC
#include"mpif.h"
    integer status(mpi_status_size)
#endif

    if (t_recv /= mype) then
       call fftscat(kpbpts,klev, rbuf,rbuflen,ibuf,0,t_recv,div_tri)
       ibuflen=1
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_recv(rbuf,rbuflen,REALTYPE,t_recv,2201+mype,localComm,status,error)
#endif
          call fftscat(kpbpts,klev, rbuf,rbuflen,ibuf,ibuflen,t_recv,div_tri)
       endif
    endif
    return
  end subroutine fft_tri_recv




  subroutine trigath(div_tri,kpbpts,klev,rb,rblen,ib,iblen,destpro,div_fft)

    IMPLICIT NONE

    real(kind=realkind) div_tri(klat_global,*)
    real(kind=realkind) div_fft(klon_global,*)
    INTEGER KPBPTS,KLEV, rblen, iblen, destpro
    REAL(KIND=REALKIND) RB(RBLEN)
    !     The IB is Index Buffer, containing index information
    !     for databuffer RB
    INTEGER IB(3,IBLEN)


    integer irec,jslab,i,imin,jlen,jlev,rempro,remrec

    integer error,kproc,tri_tag,nrempro
    integer islab,imax,jmin
    integer ncol,isend
    integer j, ii,jj,kk, ndx, lenw
    integer rbc,ibc

    rbc = 0
    ibc = 1
    if (destpro==mype) then
       irec = 1
       do jslab=1,nslab_fft
          jmin = jmin_fft(jslab)
          jlen = jlen_fft(jslab)
          jlev = jlev_fft(jslab)
          do i=kpbpts+2,klon_global-kpbpts-2
             rempro = kpro_tri(i,jlev)
             remrec = krec_tri(i,jlev)
             if (rempro==mype) then
                call scopy(jlen,div_tri(jmin,remrec),1,div_fft(i,irec),klon_global)
             endif
          enddo
          irec = irec + jlen
       enddo

    else

       do islab=1,nslab_tri
          imin = imin_tri(islab)
          imax = imax_tri(islab)
          jlev = jlev_tri(islab)
          kproc = nproc_tri(jlev)
          remrec=krec_tri(imin,jlev)
          jmin = 1
          do j=1,kproc
             if (pe_tri(jlev)+j-1 == destpro) then
                tri_tag=2100+imin+(jmin-1+(jlev-1)*klat_global)*klev_global
                if (destpro /= mype) then
                   do jj=jmin,jmin+jlen_tri(jlev,j)-1
                      do kk=remrec,remrec+imax-imin
                         rbc = rbc +1
                         rb(rbc) = div_tri(jj,kk)
                      enddo
                   enddo
                endif
             endif
             jmin = jmin + jlen_tri(jlev,j)
          enddo
       enddo
    endif

    rblen = rbc
    iblen = ibc
    return
  end subroutine trigath

  subroutine triscat(kpbpts,klev,rb,rblen,ib,iblen,destpro,div_fft)

    implicit none
    real(kind=realkind) div_fft(klon_global,*)

    integer kpbpts,klev, rblen, iblen, destpro
    real(kind=realkind) rb(rblen)
    integer ib(3,iblen)
    integer irec,jslab,i,jmin,jlen,jlev,rempro,remrec

    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax
    integer ncol,isend,irecv
    integer j, ii,jj,kk, ndx, lenw,sumlenw
    integer rbc,ibc,jfound

    sumlenw=0
    irec = 1
    rbc = 0
    do jslab=1,nslab_fft
       jmin = jmin_fft(jslab)
       jlen = jlen_fft(jslab)
       jlev = jlev_fft(jslab)
       irecv=kpbpts+2
       ncol=1
       do i=kpbpts+2,klon_global-kpbpts-2
          rempro = kpro_tri(i,jlev)
          if (i == klon_global-kpbpts-2) then
             nrempro = -1
          else
             nrempro = kpro_tri(i+1,jlev)
          endif
          if(nrempro /= rempro) then
             if (rempro==destpro) then
                tri_tag=2100+irecv+(jmin-1+(jlev-1)*klat_global)*klev_global
                lenw = ncol*jlen
                sumlenw = sumlenw + lenw
                if (iblen /= 0) then
                   ibc=1
                   if (destpro /= mype) then
                      do kk=irec,irec+jlen-1
                         do ii=irecv,irecv+ncol-1
                            rbc = rbc + 1
                            div_fft(ii,kk) = rb(rbc)
                         enddo
                      enddo
                   endif
                endif
             endif
             irecv=i+1
             ncol=1
          else
             ncol=ncol+1
          endif
       enddo
       irec = irec + jlen
    enddo
    rblen = sumlenw
    return
  end subroutine triscat

  subroutine tri_send(div_tri,kpbpts,klev,t_send,div_fft)

    implicit none
    real(kind=realkind) div_tri(klat_global,*)
    real(kind=realkind) div_fft(klon_global,*)
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer kpbpts,klev,t_send
    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)
    integer ismax
    data ismax/0/
    save ismax
    integer rsmax
    data rsmax/0/
    save rsmax

    rbuflen=maxrbuf
    ibuflen=maxibuf-1
    call trigath(div_tri,kpbpts,klev,rbuf,rbuflen,ibuf,ibuflen,t_send,div_fft)
    if (ibuflen>=maxibuf) then
       write(*,*) 'pe',mype,'FATAL in TRI_SEND(), MAXIBUF TOO SMALL',maxibuf,ibuflen
    endif
    if (rbuflen>=maxrbuf) then
       write(*,*) 'pe',mype,'FATAL in TRI_SEND(), MAXRBUF TOO SMALL ',maxrbuf,rbuflen
    endif
    if (t_send /= mype) then
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_send(rbuf,rbuflen,REALTYPE,t_send,2200+t_send,localComm,error)
#endif
       endif
    endif
    return
  end subroutine tri_send


  subroutine tri_recv(kpbpts,klev,t_recv,div_fft)

    implicit none
    real(kind=realkind) div_fft(klon_global,*)
    integer kpbpts,klev,t_recv
#ifdef MPI_SRC
#include"mpif.h"
    integer status(mpi_status_size)
#endif
    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)

    if (t_recv /= mype) then
       call triscat(kpbpts,klev, rbuf,rbuflen,ibuf,0,t_recv,div_fft)
       ibuflen=1
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_recv(rbuf,rbuflen,REALTYPE,t_recv,2200+mype,localComm,status,error)
#endif
          call triscat(kpbpts,klev, rbuf,rbuflen,ibuf,ibuflen,t_recv,div_fft)

       endif
    endif
    return
  end subroutine tri_recv


  subroutine f2dgath(div,klon,klat,klev,rb,rblen,ib,iblen,destpro,div_fft)


    implicit none
    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv

    integer kpbpts, rblen, iblen, destpro
    real(kind=realkind) rb(rblen)
    !     The IB is Index Buffer, containing index information
    !     for databuffer RB
    integer ib(3,iblen)
    integer irec,jslab,i,jmin,jlen,jlev,rempro,remrec
    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax,jmax,iminr
    integer ncol,isend,fft_tag
    integer j, k,ii,jj,kk, ndx, lenw
    integer rbc,ibc,ipe
    integer nvalues,nvaluesr,idatastartr

    call jlimits(klon,klat,imin,imax,jmin,jmax)

    nvalues = klon - 2
    if( atleft ) then
       nvalues = nvalues + 1
    endif
    if( atright ) then
       nvalues = nvalues + 1
    endif

    irec = 1
    rbc = 0
    ibc = 1
    if (destpro==mype) then
       do k=1,klev
          do j=jmin,jmax
             rempro = kpro_fft(jdatastart+j-1,k)
             remrec = krec_fft(jdatastart+j-1,k)
             if (rempro==mype) then
                call scopy(nvalues,div_fft(idatastart+imin-1,remrec),1,div(imin,j,k),1)
             endif
          enddo
       enddo
    else
       do k=1,klev
          do j=1,klat_global
             rempro = kpro_fft(j,k)
             if(rempro==mype) then
                do ipe=1,nprocx
                   if(pe_fft(j,ipe)==destpro)then
                      idatastartr=idatastart_ar2d(ipe)
                      nvaluesr=klon_ar2d(ipe)-2
                      iminr = 2
                      if(ipe==1) then
                         iminr = 1
                         nvaluesr = nvaluesr + 1
                      endif
                      if(ipe==nprocx) then
                         nvaluesr = nvaluesr + 1
                      endif
                      fft_tag=2010+idatastartr+iminr-1+((j-1)+(k-1)*klat_global)*klon_global
                      remrec = krec_fft(j,k)
                      if (destpro/=mype) then
                         do i=idatastartr+iminr-1, idatastartr+iminr-1+nvaluesr-1
                            rbc = rbc+1
                            rb(rbc)=div_fft(i,remrec)
                         enddo
                      endif
                   endif
                enddo
             endif
          enddo
       enddo
    endif

    rblen = rbc
    iblen = ibc

    return
  end subroutine f2dgath



  subroutine f2dscat(div,klon,klat,klev, rb,rblen,ib,iblen,destpro)


    implicit none

    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    integer t_recv
    integer kpbpts, rblen, iblen, destpro
    real(kind=realkind) rb(rblen)
    integer ib(3,iblen)
    integer irec,jslab,jmin,jmax,jlen,jlev,rempro,remrec
    integer error,kproc,tri_tag,nrempro
    integer islab,imin,imax,idatastartr,iminr
    integer ncol,isend,fft_tag,nvalues,ipe
    integer i,j,k, ii,jj,kk, ndx, lenw,sumlenw
    integer rbc,ibc

    call jlimits(klon,klat,imin,imax,jmin,jmax)

    nvalues = klon - 2

    if( atleft ) then
       nvalues = nvalues + 1
    endif

    if( atright ) then
       nvalues = nvalues + 1
    endif

    sumlenw = 0
    rbc = 0
    do k=1,klev
       do j=jmin,jmax
          rempro = kpro_fft(jdatastart+j-1,k)
          if (rempro==destpro) then
             fft_tag=2010+idatastart+imin-1+((jdatastart+j-2)+(k-1)*klat_global)*klon_global
             sumlenw = sumlenw + nvalues
             if (iblen /= 0) then
                ibc=1
                do i=imin,imin+nvalues-1
                   rbc=rbc+1
                   div(i,j,k) = rb(rbc)
                enddo
             endif
          endif
       enddo
    enddo
    rblen = sumlenw

    return
  end subroutine f2dscat



  subroutine f2dsend(div,klon,klat,klev,t_send,div_fft)

    implicit none
    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    real(kind=realkind) div_fft(klon_global,*)
    integer t_recv
    integer kpbpts,t_send


    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)
#ifdef MPI_SRC
#include"mpif.h"      
#endif      
    rbuflen=maxrbuf
    ibuflen=maxibuf-1
    call f2dgath(div,klon,klat,klev,rbuf,rbuflen,ibuf,ibuflen,t_send,div_fft)
    ibuflen=1
    if (ibuflen>=maxibuf) then
       write(*,*) 'pe',mype,'FATAL in F2DSEND(), MAXIBUF TOO SMALL', maxibuf,ibuflen
       call stop_program('')
    endif
    if (rbuflen>=maxrbuf) then
       write(*,*) 'pe',mype,'FATAL in F2DSEND(), MAXRBUF TOO SMALL',maxrbuf,rbuflen
       call stop_program('')
    endif


    if (t_send/= mype) then
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call mpi_send(rbuf,rbuflen,REALTYPE,t_send,2200+t_send,LOCALCOMM,error)
#endif  
       endif
    endif
    return
  end subroutine f2dsend

  subroutine f2drecv(div,klon,klat,klev,t_recv)
    implicit none
    integer klon,klat,klev
    real(kind=realkind) div(klon,klat,klev)
    integer t_recv
    integer rbuflen, ibuflen
    integer maxrbuf,maxibuf,error
    parameter(maxibuf=100000,maxrbuf=3000000)
    real(kind=realkind) rbuf(maxrbuf)
    integer ibuf(3,maxibuf)
#ifdef MPI_SRC
#include"mpif.h"  
    integer status(mpi_status_size)   
#endif 
    if (t_recv /= mype) then
       call f2dscat(div,klon,klat,klev,rbuf,rbuflen,ibuf,0,t_recv)
       ibuflen=1
       if (rbuflen /= 0) then
#ifdef MPI_SRC
          call MPI_recv(rbuf,rbuflen,REALTYPE, t_recv,2200+mype,localComm,status,error)
#endif
          call f2dscat(div,klon,klat,klev,rbuf,rbuflen,ibuf,ibuflen,t_recv)
       endif
    endif
    return
  end subroutine f2drecv

  SUBROUTINE SCOPY (N, SX, INCX, SY, INCY)  
    !     COPIES A VECTOR, X, TO A VECTOR, Y.
    !     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
    !     JACK DONGARRA, LINPACK, 3/11/78.
    !     MODIFIED 12/3/93, ARRAY(1) DECLARATIONS CHANGED TO ARRAY(*)
    REAL(KIND=REALKIND) :: SX ( * ), SY ( * )  

    INTEGER :: I, INCX, INCY, IX, IY, M, MP1, N  
    IF (N<=0) RETURN  


    IF (INCX==1.AND.INCY==1) GOTO 20  
    !     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
    !     NOT EQUAL TO 1
    IX = 1  
    IY = 1  
    IF (INCX<0) IX = ( - N + 1) * INCX + 1  
    IF (INCY<0) IY = ( - N + 1) * INCY + 1  
    DO 10 I = 1, N  
       SY (IY) = SX (IX)  
       IX = IX + INCX  
       IY = IY + INCY  
10  END DO




    RETURN  
    !     CODE FOR BOTH INCREMENTS EQUAL TO 1
    !     CLEAN-UP LOOP
20  M = MOD (N, 7)  
    IF (M==0) GOTO 40  
    DO 30 I = 1, M  
       SY (I) = SX (I)  
30  END DO
    IF (N<7) RETURN  
40  MP1 = M + 1  
    DO 50 I = MP1, N, 7  
       SY (I) = SX (I)  
       SY (I + 1) = SX (I + 1)  
       SY (I + 2) = SX (I + 2)  
       SY (I + 3) = SX (I + 3)  
       SY (I + 4) = SX (I + 4)  
       SY (I + 5) = SX (I + 5)  
       SY (I + 6) = SX (I + 6)  
50  END DO
    RETURN  


  END SUBROUTINE SCOPY


end module transpose
