module decomp
  implicit none
  private
#ifndef SINGLE_PRECISION
  integer,public,parameter::realkind=8
#else
  integer,public,parameter::realkind=4
#endif

  integer,public,save::localComm
#ifdef MPI_SRC
#include"mpif.h"
  integer,public,save::iside2p,jside2p
  integer,public,save::isidep,jsidep
  integer,public,save::iside,jside
  integer,public,save::iside_levp1,jside_levp1
  integer,public,save::iside_lev,jside_lev
  integer,public,save::iside_lev2,jside_lev2
  integer,public,save::iside_lev3,jside_lev3
  integer,public,save::iside_lev4,jside_lev4
  integer,public,save::iside_lev5,jside_lev5
  integer,public,save::iside_lev6,jside_lev6
  integer,public,save::REALTYPE
#endif
  integer,public,save::mype                                   !processor id in localComm
  integer,public,save::nproc                  !number of processors in 2d
  integer,public,save::idatastart,jdatastart                 !local to global index mappings
  integer,public,save::pe_top,pe_base,pe_right,pe_left !neighbor processors
  integer,public,save::klon_global,klat_global,klev_global   !global number of points
  logical,public,save::atbase,attop,atright,atleft


  integer,public,save::nprocx,nprocy


  integer,allocatable,dimension(:),public,save::klon_ar,klat_ar,idatastart_ar,jdatastart_ar 
  integer,allocatable,dimension(:),public,save::klon_ar2d,klat_ar2d,idatastart_ar2d,jdatastart_ar2d
  public decompose
  public jlimits
  public swap
  public swap2
  public swap3
  public swap4
  public swap5
  public swap6
  public swap2d
  public swap_ps
  public swapklevp1
  public readDeconam
  public myshare
  public colfld
  public factorize
  public stop_program
  public HeaviSide
contains
  real(kind=realkind) function HeaviSide(x,xlim,a,b,eps)
    !computes the continuous version of the Heaviside function
    !if(x<xlim)then
    !   f=a
    !else
    !   f=b
    !endif
    real(kind=realkind),intent(in)::x,xlim,a,b,eps

    if(abs(x-xlim)<eps)then
       HeaviSide = 1._realkind/eps *( (a-b)*x +a*eps + (b-a)*(xlim+0.5_realkind*eps))
    else
       if(x<xlim)then
          HeaviSide = a
       else
          HeaviSide = b
       endif
    endif
    return

  end function HeaviSide


  subroutine myshare(za_global,za_local,klon,klat) 
    !     picks up from a global array the sub area belonging to the
    !     current node(pe)
    implicit none  

    integer::klon,klat  
    real(kind=realkind)::za_global(klon_global,klat_global),za_local(klon,klat)

    integer::jx,jy  
    do jy = 1,klat  
       do jx = 1,klon  
          za_local(jx,jy)= za_global(jx+idatastart-1,jy+jdatastart-1)
       enddo
    enddo
    return  
  end subroutine myshare

  subroutine jlimits(klon,klat,istart,istop,jstart,jstop) 
    implicit none  
    integer,intent(in)::klon,klat  
    integer,intent(out)::istart,istop,jstart,jstop
    jstart = 2  
    jstop = klat - 1  
    istart = 2  
    istop = klon - 1  

    if(atbase)then  
       jstart = 1  
    endif
    if(attop)then  
       jstop = klat  
    endif
    if(atleft)then  
       istart = 1  
    endif
    if(atright)then  
       istop = klon  
    endif
    return  
  end subroutine jlimits

  subroutine readDeconam()
    implicit none

    real(kind=realkind)::south,west,dlon,dlat,polon,polat !used as dummy
    integer::nprocz
#if defined(MPI_SRC)
#include"mpif.h"
    integer::ierr,ibuf(6)
#endif    
    namelist/namprocessor/nprocx,nprocy
    namelist/domain/south,west,dlon,dlat,polon,polat,klon_global,klat_global,klev_global
    nprocx=-1
    nprocy=-1


    if(mype==0)then

       open(57,file='namelists.dat',status='old')
       read(57,nml=namprocessor)
       close(57)

       open(57,file='namelists.dat',status='old')
       read(57,nml=domain)
       close(57)

       write(6,nml=domain)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=domain)
       close(1)
    endif
#ifdef MPI_SRC
    ibuf = 666
    if(mype==0)then
       ibuf(1) = nproc        !decomp        namprc
       ibuf(2) = nprocx       !decomp        namprc
       ibuf(3) = nprocy       !decomp        namprc
       ibuf(4) = klon_global  !decomp        domain
       ibuf(5) = klat_global  !decomp        domain
       ibuf(6) = klev_global  !decomp        domain
    endif
    call mpi_bcast(ibuf,6,mpi_integer,0,localComm,ierr)
    nproc       = ibuf(1)
    nprocx      = ibuf(2)
    nprocy      = ibuf(3)
    klon_global = ibuf(4)
    klat_global = ibuf(5)
    klev_global = ibuf(6)
#endif

    if(nprocx==-1.or.nprocy==-1)then
       call factorize(nproc,klon_global,klat_global,1,nprocx,nprocy,nprocz)
    endif
    if(nproc/=nprocx*nprocy)then
       print *, __FILE__,__LINE__,nproc,'!=',nprocx,'*',nprocy
       call stop_program('')
    endif

    call allocate_decomp_memory()
    if(mype==0)then
       write(6,nml=namprocessor)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namprocessor)
       close(1)
    endif
  end subroutine readDeconam

  subroutine allocate_decomp_memory()
    implicit none
    allocate(klon_ar(nproc),klat_ar(nproc),idatastart_ar(nproc),jdatastart_ar(nproc))
    allocate(klon_ar2d(nprocx),klat_ar2d(nprocy),idatastart_ar2d(nprocx),jdatastart_ar2d(nprocy))
  end subroutine allocate_decomp_memory


  subroutine decompose(klon,klat,klev) 
    implicit none  

    integer,intent(out)::klon,klat,klev  
    integer::i,j,rem_x,rem_y  
    integer::ierr
    logical::debug = .false.  
    integer,save::i_pe_grid,j_pe_grid                  !processor id in 2d processor grid

    !     2d pe identities
    i_pe_grid = mod(mype,nprocx) ![0...nprocx-1]
    j_pe_grid = mype/nprocx      ![0...nprocy-1]

    !     who are my neighbours?
    pe_base = -1
    pe_top = -1
    pe_right = -1
    pe_left = -1
    if(j_pe_grid>0)then
       pe_base = i_pe_grid+(j_pe_grid-1)*nprocx  
    endif
    if(j_pe_grid<nprocy)then
       pe_top = i_pe_grid+(j_pe_grid+1)*nprocx  
    endif
    if(i_pe_grid<nprocx)then
       pe_right = i_pe_grid+1 + j_pe_grid*nprocx  
    endif
    if(i_pe_grid>0)then
       pe_left = i_pe_grid-1 + j_pe_grid*nprocx  
    endif


    !     what is remainder?
    rem_x = mod((klon_global-2),nprocx) 
    do i = 1,nprocx  
       klon_ar2d(i)= (klon_global-2)/nprocx 
       if(i - 1<rem_x) klon_ar2d(i) = klon_ar2d(i) + 1  !     add extra to low pe's
       klon_ar2d(i) = klon_ar2d(i) + 2  !     add on halo
    enddo

    do i = 1,nprocx  
       if(i-1 < rem_x) then  
          idatastart_ar2d(i) = 1+(i-1)*(klon_ar2d(i) - 2)  
       else  
          idatastart_ar2d(i) = 1+(i-1)*(klon_ar2d(i) - 2) + rem_x  
       endif
    enddo
    klon = klon_ar2d(i_pe_grid+1)  
    idatastart = idatastart_ar2d(i_pe_grid+1)  


    rem_y = mod((klat_global-2),nprocy) 
    do j = 1,nprocy  
       klat_ar2d(j) = (klat_global-2)/nprocy 
       if(j - 1<rem_y) klat_ar2d(j) = klat_ar2d(j) + 1  !     add extra to low pe's
       klat_ar2d(j) = klat_ar2d(j) + 2  !     add on halo
    enddo
    !     position of the start of the local data
    do j = 1,nprocy  
       if(j - 1<rem_y) then  
          jdatastart_ar2d(j) = 1 +(j-1)*(klat_ar2d(j)-2)  
       else  
          jdatastart_ar2d(j) = 1 +(j-1)*(klat_ar2d(j)-2) + rem_y  
       endif
    enddo
    klat = klat_ar2d(j_pe_grid+1)  
    jdatastart = jdatastart_ar2d(j_pe_grid+1)  

    do j = 1,nprocy  
       do i = 1,nprocx  
          klon_ar(i+(j-1)*nprocx) = klon_ar2d(i)
          idatastart_ar(i+(j-1)*nprocx) = idatastart_ar2d(i)
          klat_ar(i+(j-1)*nprocx) = klat_ar2d(j)
          jdatastart_ar(i+(j-1)*nprocx) = jdatastart_ar2d(j)          
       enddo
    enddo
    klev = klev_global


    !     are we at the base,top,right or left boundary of the region?
    atleft = .false.  
    atright = .false.  
    atbase = .false.  
    attop = .false.  

    if(i_pe_grid==0) then  
       atleft = .true.  
       pe_left = - 1  
    endif
    if(i_pe_grid==(nprocx - 1) ) then  
       atright = .true.  
       pe_right = - 1  
    endif
    if(j_pe_grid==0) then  
       atbase = .true.  
       pe_base = - 1  
    endif
    if(j_pe_grid==(nprocy - 1) ) then  
       attop = .true.  
       pe_top = - 1  
    endif

    if(mype==0.and.debug) then  
       write(6,*) ' pe i_pe_grid j_pe_grid klon klat idatastart,jdat astart'
       do j=1,nprocy  
          do i=1,nprocx  
             write(6,*)i,j,klon_ar2d(i),klat_ar2d(j),idatastart_ar2d(i),jdatastart_ar2d(j)
          enddo
       enddo
    endif

    !     helmholz solver: initial values,possibly to changed

#ifdef MPI_SRC
    if(realkind==4)then
       realtype = MPI_REAL
    elseif(realkind==8)then
       realtype = MPI_DOUBLE_PRECISION
    else
       call stop_program('r8 or r4 needs to be specified')
    endif

    call MPI_TYPE_CONTIGUOUS(2*(klon+1),REALTYPE,jside2p,ierr)
    call MPI_TYPE_COMMIT(jside2p,ierr)
    call MPI_TYPE_CONTIGUOUS(klon+1,REALTYPE,jsidep,ierr)
    call MPI_TYPE_COMMIT(jsidep,ierr)
    call MPI_TYPE_VECTOR(klat+1,2,klon+1,REALTYPE,iside2p,ierr)
    call MPI_TYPE_COMMIT(iside2p,ierr)
    call MPI_TYPE_VECTOR(klat+1,1,klon+1,REALTYPE,isidep,ierr)
    call MPI_TYPE_COMMIT(isidep,ierr)
    call MPI_TYPE_CONTIGUOUS(klon,REALTYPE,jside,ierr)
    call MPI_TYPE_COMMIT(jside,ierr)
    call MPI_TYPE_VECTOR(klat,1,klon,REALTYPE,iside,ierr)
    call MPI_TYPE_COMMIT(iside,ierr)

    call MPI_TYPE_VECTOR(klev,klon,klon*klat,REALTYPE,jside_lev,ierr)
    call MPI_TYPE_COMMIT(jside_lev,ierr)
    call MPI_TYPE_VECTOR(klat*klev,1,klon,REALTYPE,iside_lev,ierr)
    call MPI_TYPE_COMMIT(iside_lev,ierr)

    call MPI_TYPE_VECTOR(klev+1,klon,klon*klat,REALTYPE,jside_levp1,ierr)
    call MPI_TYPE_COMMIT(jside_levp1,ierr)
    call MPI_TYPE_VECTOR(klat*(klev+1),1,klon,REALTYPE,iside_levp1,ierr)
    call MPI_TYPE_COMMIT(iside_levp1,ierr)


    call MPI_TYPE_VECTOR(2*klev,klon,klon*klat,REALTYPE,jside_lev2,ierr)
    call MPI_TYPE_COMMIT(jside_lev2,ierr)

    call MPI_TYPE_VECTOR(2*klat*klev,1,klon,REALTYPE,iside_lev2,ierr)
    call MPI_TYPE_COMMIT(iside_lev2,ierr)

    call MPI_TYPE_VECTOR(3*klev,klon,klon*klat,REALTYPE,jside_lev3,ierr)
    call MPI_TYPE_COMMIT(jside_lev3,ierr)

    call MPI_TYPE_VECTOR(3*klat*klev,1,klon,REALTYPE,iside_lev3,ierr)
    call MPI_TYPE_COMMIT(iside_lev3,ierr)

    call MPI_TYPE_VECTOR(4*klev,klon,klon*klat,REALTYPE,jside_lev4,ierr)
    call MPI_TYPE_COMMIT(jside_lev4,ierr)

    call MPI_TYPE_VECTOR(4*klat*klev,1,klon,REALTYPE,iside_lev4,ierr)
    call MPI_TYPE_COMMIT(iside_lev4,ierr)

    call MPI_TYPE_VECTOR(5*klev,klon,klon*klat,REALTYPE,jside_lev5,ierr)
    call MPI_TYPE_COMMIT(jside_lev5,ierr)

    call MPI_TYPE_VECTOR(5*klat*klev,1,klon,REALTYPE,iside_lev5,ierr)
    call MPI_TYPE_COMMIT(iside_lev5,ierr)

    call MPI_TYPE_VECTOR(6*klev,klon,klon*klat,REALTYPE,jside_lev6,ierr)
    call MPI_TYPE_COMMIT(jside_lev6,ierr)

    call MPI_TYPE_VECTOR(6*klat*klev,1,klon,REALTYPE,iside_lev6,ierr)
    call MPI_TYPE_COMMIT(iside_lev6,ierr)

#endif

    return  
  end subroutine decompose



  subroutine factorize(npr,nig,njg,nkg,np1,np2,np3)
    integer,intent(in)::npr,nig,njg,nkg
    integer,intent(inout)::np1,np2,np3
    real(kind=8)::nor,norloc
    integer::p1,p2,p3

    nor = real(nig+2*njg+nkg,kind(nor))
    do p3=1,npr
       do p1=1,npr
          if(mod(npr,p1*p3)==0)then
             p2 = npr/(p1*p3)
             norloc = abs(real(nig,kind(norloc))/real(p1,kind(norloc)) -       &
                  real(njg,kind(norloc))/real(p2,kind(norloc))  ) +   &
                  abs(real(njg,kind(norloc))/real(p2,kind(norloc)) -  &
                  real(nkg,kind(norloc))/real(p3,kind(norloc)) )
             if(norloc<nor)then
                np1 = p1
                np2 = p2
                np3 = p3
                nor = norloc
             endif
          endif
       enddo
    enddo
  end subroutine factorize

  subroutine swap2d(a,klon,klat)
    !  purpose:
    !     updates the halo zone from the neighbouring processors
    !  author:
    !      kalle eerola   hirlam    1998

    implicit none
    integer,intent(in)::klon,klat
    real(kind=realkind),intent(inout)::a(klon,klat)

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif

    if(nproc==1) return
#ifdef MPI_SRC
    if(pe_top/=-1)then   
       call mpi_recv(a(1,klat),1,jside,pe_top,51,localComm,status,info)
    endif
    if(pe_base/=-1) then
       call mpi_send(a(1,2),1,jside,pe_base,51,localcomm,info)
       call mpi_recv(a,1,jside,pe_base,61,localcomm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(a(1,klat-1),1,jside,pe_top,61,localComm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(a(klon,1),1,iside,pe_right,21,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(a(2,1),1,iside,pe_left,21,localComm,info)
       call mpi_recv(a,1,iside,pe_left,31,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(a(klon-1,1),1,iside,pe_right,31,localComm,info)
    endif

#ifdef USE_BARRIER
    call mpi_barrier(localComm,info)
#endif
#endif

    return
  end subroutine swap2d


  subroutine swap(a,klon,klat,klev)
    !  purpose:
    !     updates the halo zone from the neighbouring processors
    !  author:
    !      kalle eerola   hirlam    1998

    implicit none
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(inout)::a(klon,klat,klev)

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif

    if(nproc==1)return

#ifdef MPI_SRC

    if(pe_top/=-1)then   
       call mpi_recv(a(1,klat,1),1,jside_lev,pe_top,11051,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(a(1,2,1),1,jside_lev,pe_base,11051,localcomm,info)
       call mpi_recv(a,1,jside_lev,pe_base,11061,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(a(1,klat-1,1),1,jside_lev,pe_top,11061,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(a(klon,1,1),1,iside_lev,pe_right,11211,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(a(2,1,1),1,iside_lev,pe_left,11211,localComm,info)
       call mpi_recv(a,1,iside_lev,pe_left,11311,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(a(klon-1,1,1),1,iside_lev,pe_right,11311,localComm,info)
    endif


#endif

    return
  end subroutine swap


  subroutine swapklevp1(a,klon,klat,klev)
    !  purpose:
    !     updates the halo zone from the neighbouring processors
    !  author:
    !      kalle eerola   hirlam    1998

    implicit none
    integer,intent(in)::klon,klat,klev
    real(kind=realkind),intent(inout)::a(klon,klat,klev)

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif

    if(nproc==1)return

#ifdef MPI_SRC

    if(pe_top/=-1)then   
       call mpi_recv(a(1,klat,1),1,jside_levp1,pe_top,51,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(a(1,2,1),1,jside_levp1,pe_base,51,localcomm,info)
       call mpi_recv(a,1,jside_levp1,pe_base,61,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(a(1,klat-1,1),1,jside_levp1,pe_top,61,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(a(klon,1,1),1,iside_levp1,pe_right,211,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(a(2,1,1),1,iside_levp1,pe_left,211,localComm,info)
       call mpi_recv(a,1,iside_levp1,pe_left,311,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(a(klon-1,1,1),1,iside_levp1,pe_right,311,localComm,info)
    endif


#endif

    return
  end subroutine swapklevp1


  subroutine swap2(a,b,klon,klat,klev)
    !  PURPOSE:
    !     Multiple variable updates the halo zone from the
    !     neighbouring processors
    !
    !  DESCRIPTION:
    !
    !  MODIFICATIONS:
    !
    !  AUTHOR:
    !     JUSSI MAKI      CSC       1999

    implicit none

    integer::klon,klat,klev
    real(kind=realkind)::a(klon,klat,klev)
    real(kind=realkind)::b(klon,klat,klev)

    ! LEN_X IS THE LENGTH OF THE BUFFER TO SEND IN X-DIRECTION
!!$    INTEGER LEN_NS,LEN_EW
#ifdef MPI_SRC
#include"mpif.h"
    integer::info
    integer::status(mpi_status_size)
#endif
    real(kind=realkind),allocatable,dimension(:,:,:,:)::tmp

    if(nproc==1)return
#ifdef MPI_SRC
    allocate(tmp(klon,klat,klev,2))
    tmp(:,:,:,1) = a
    tmp(:,:,:,2) = b
    if(pe_top/=-1)then   
       call mpi_recv(tmp(1,klat,1,1),1,jside_lev2,pe_top,110212,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(tmp(1,2,1,1),1,jside_lev2,pe_base,110212,localcomm,info)
       call mpi_recv(tmp,1,jside_lev2,pe_base,110612,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(tmp(1,klat-1,1,1),1,jside_lev2,pe_top,110612,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(tmp(klon,1,1,1),1,iside_lev2,pe_right,112112,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(tmp(2,1,1,1),1,iside_lev2,pe_left,112112,localComm,info)
       call mpi_recv(tmp,1,iside_lev2,pe_left,113112,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(tmp(klon-1,1,1,1),1,iside_lev2,pe_right,113112,localComm,info)
    endif
    a = tmp(:,:,:,1) 
    b = tmp(:,:,:,2) 
    deallocate(tmp)
#endif

    !    call swap(a,klon,klat,klev)
    !    call swap(b,klon,klat,klev)


    return

  end subroutine swap2


  subroutine swap3(a,b,c,klon,klat,klev)

    !  PURPOSE:
    !     Multiple variable updates the halo zone from the
    !     neighbouring processors
    !
    !  DESCRIPTION:
    !
    !  MODIFICATIONS:
    !
    !  AUTHOR:
    !     JUSSI MAKI      CSC       1999

    implicit none

    integer::klon,klat,klev
    real(kind=realkind)::a(klon,klat,klev)
    real(kind=realkind)::b(klon,klat,klev)
    real(kind=realkind)::c(klon,klat,klev)
    real(kind=realkind),allocatable,dimension(:,:,:,:)::tmp

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif

    if(nproc==1)return
#ifdef MPI_SRC
    allocate(tmp(klon,klat,klev,3))
    tmp(:,:,:,1) = a
    tmp(:,:,:,2) = b
    tmp(:,:,:,3) = c
    if(pe_top/=-1)then   
       call mpi_recv(tmp(1,klat,1,1),1,jside_lev3,pe_top,110313,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(tmp(1,2,1,1),1,jside_lev3,pe_base,110313,localcomm,info)
       call mpi_recv(tmp,1,jside_lev3,pe_base,110613,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(tmp(1,klat-1,1,1),1,jside_lev3,pe_top,110613,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(tmp(klon,1,1,1),1,iside_lev3,pe_right,112113,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(tmp(2,1,1,1),1,iside_lev3,pe_left,112113,localComm,info)
       call mpi_recv(tmp,1,iside_lev3,pe_left,113113,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(tmp(klon-1,1,1,1),1,iside_lev3,pe_right,113113,localComm,info)
    endif
    a = tmp(:,:,:,1) 
    b = tmp(:,:,:,2) 
    c = tmp(:,:,:,3) 
    deallocate(tmp)
#endif
!!$    call swap(a,klon,klat,klev)
!!$    call swap(b,klon,klat,klev)
!!$    call swap(c,klon,klat,klev)

    return
  end subroutine swap3


  subroutine swap4(a,b,c,d,klon,klat,klev)

    !  PURPOSE:
    !     Multiple variable updates the halo zone from the
    !     neighbouring processors
    !
    !  DESCRIPTION:
    !
    !  MODIFICATIONS:
    !
    !  AUTHOR:
    !     JUSSI MAKI      CSC       1999

    implicit none

    integer::klon,klat,klev
    real(kind=realkind)::a(klon,klat,klev)
    real(kind=realkind)::b(klon,klat,klev)
    real(kind=realkind)::c(klon,klat,klev)
    real(kind=realkind)::d(klon,klat,klev)

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif
    real(kind=realkind),allocatable,dimension(:,:,:,:)::tmp

    if(nproc==1)return

#ifdef MPI_SRC
    allocate(tmp(klon,klat,klev,4))
    tmp(:,:,:,1) = a
    tmp(:,:,:,2) = b
    tmp(:,:,:,3) = c
    tmp(:,:,:,4) = d

    if(pe_top/=-1)then   
       call mpi_recv(tmp(1,klat,1,1),1,jside_lev4,pe_top,110414,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(tmp(1,2,1,1),1,jside_lev4,pe_base,110414,localcomm,info)
       call mpi_recv(tmp,1,jside_lev4,pe_base,110614,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(tmp(1,klat-1,1,1),1,jside_lev4,pe_top,110614,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(tmp(klon,1,1,1),1,iside_lev4,pe_right,112114,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(tmp(2,1,1,1),1,iside_lev4,pe_left,112114,localComm,info)
       call mpi_recv(tmp,1,iside_lev4,pe_left,113114,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(tmp(klon-1,1,1,1),1,iside_lev4,pe_right,113114,localComm,info)
    endif
    a = tmp(:,:,:,1) 
    b = tmp(:,:,:,2) 
    c = tmp(:,:,:,3) 
    d = tmp(:,:,:,4) 

    deallocate(tmp)
#endif
    !    call swap(a,klon,klat,klev)
    !    call swap(b,klon,klat,klev)
    !    call swap(c,klon,klat,klev)
    !    call swap(d,klon,klat,klev)
    return
  end subroutine swap4


  subroutine swap5(a,b,c,d,e,klon,klat,klev)

    !  PURPOSE:
    !     Multiple variable updates the halo zone from the
    !     neighbouring processors
    !
    !  DESCRIPTION:
    !
    !  MODIFICATIONS:
    !
    !  AUTHOR:
    !     JUSSI MAKI      CSC       1999

    implicit none

    integer::klon,klat,klev
    real(kind=realkind),target::a(klon,klat,klev)
    real(kind=realkind),target::b(klon,klat,klev)
    real(kind=realkind),target::c(klon,klat,klev)
    real(kind=realkind),target::d(klon,klat,klev)
    real(kind=realkind),target::e(klon,klat,klev)

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif
    real(kind=realkind),allocatable,dimension(:,:,:,:)::tmp

    if(nproc==1)return

#ifdef MPI_SRC
    allocate(tmp(klon,klat,klev,5))
    tmp(:,:,:,1) = a
    tmp(:,:,:,2) = b
    tmp(:,:,:,3) = c
    tmp(:,:,:,4) = d
    tmp(:,:,:,5) = e
    if(pe_top/=-1)then   
       call mpi_recv(tmp(1,klat,1,1),1,jside_lev5,pe_top,110515,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(tmp(1,2,1,1),1,jside_lev5,pe_base,110515,localcomm,info)
       call mpi_recv(tmp,1,jside_lev5,pe_base,110615,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(tmp(1,klat-1,1,1),1,jside_lev5,pe_top,110615,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(tmp(klon,1,1,1),1,iside_lev5,pe_right,112115,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(tmp(2,1,1,1),1,iside_lev5,pe_left,112115,localComm,info)
       call mpi_recv(tmp,1,iside_lev5,pe_left,113115,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(tmp(klon-1,1,1,1),1,iside_lev5,pe_right,113115,localComm,info)
    endif
    a = tmp(:,:,:,1) 
    b = tmp(:,:,:,2) 
    c = tmp(:,:,:,3) 
    d = tmp(:,:,:,4) 
    e = tmp(:,:,:,5) 
    deallocate(tmp)
#endif
    !    call swap(a,klon,klat,klev)
    !    call swap(b,klon,klat,klev)
    !    call swap(c,klon,klat,klev)
    !    call swap(d,klon,klat,klev)
    !    call swap(e,klon,klat,klev)
    return
  end subroutine swap5



  subroutine swap6(a,b,c,d,e,f,klon,klat,klev)
    implicit none

    integer::klon,klat,klev
    real(kind=realkind)::a(klon,klat,klev)
    real(kind=realkind)::b(klon,klat,klev)
    real(kind=realkind)::c(klon,klat,klev)
    real(kind=realkind)::d(klon,klat,klev)
    real(kind=realkind)::e(klon,klat,klev)
    real(kind=realkind)::f(klon,klat,klev)
    real(kind=realkind),allocatable,dimension(:,:,:,:)::tmp
#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif

    if(nproc==1)return
#ifdef MPI_SRC
    allocate(tmp(klon,klat,klev,6))
    tmp(:,:,:,1) = a
    tmp(:,:,:,2) = b
    tmp(:,:,:,3) = c
    tmp(:,:,:,4) = d
    tmp(:,:,:,5) = e
    tmp(:,:,:,6) = f
    if(pe_top/=-1)then   
       call mpi_recv(tmp(1,klat,1,1),1,jside_lev6,pe_top,110616,localComm,status,info)
    endif
    if(pe_base/=-1)then
       call mpi_send(tmp(1,2,1,1),1,jside_lev6,pe_base,110616,localcomm,info)
       call mpi_recv(tmp,1,jside_lev6,pe_base,110616,localComm,status,info)
    endif
    if(pe_top/=-1)then
       call mpi_send(tmp(1,klat-1,1,1),1,jside_lev6,pe_top,110616,localcomm,info)
    endif

    if(pe_right/=-1)then  
       call mpi_recv(tmp(klon,1,1,1),1,iside_lev6,pe_right,112116,localComm,status,info)
    endif
    if(pe_left/=-1)then
       call mpi_send(tmp(2,1,1,1),1,iside_lev6,pe_left,112116,localComm,info)
       call mpi_recv(tmp,1,iside_lev6,pe_left,113116,localComm,status,info)
    endif
    if(pe_right/=-1)then
       call mpi_send(tmp(klon-1,1,1,1),1,iside_lev6,pe_right,113116,localComm,info)
    endif
    a = tmp(:,:,:,1) 
    b = tmp(:,:,:,2) 
    c = tmp(:,:,:,3) 
    d = tmp(:,:,:,4) 
    e = tmp(:,:,:,5) 
    f = tmp(:,:,:,6) 
    deallocate(tmp)
#endif
!!$    call swap(a,klon,klat,klev)
!!$    call swap(b,klon,klat,klev)
!!$    call swap(c,klon,klat,klev)
!!$    call swap(d,klon,klat,klev)
!!$    call swap(e,klon,klat,klev)
!!$    call swap(f,klon,klat,klev)
    return
  end subroutine swap6



  subroutine colfld(kmype,p_global,p_local,klon,klat)  
    !  purpose:
    !     this subroutine collects a 2-dimensional field from all
    !     processors and creates a global field on a given processor
    !  description:
    !     given the processor number "kmype" and a local field on
    !     arrays "p_local" of dimension "klon,klat" on every
    !     processor this routines collects the field on one processor
    !     using gc-routines. in case of one processor
    !     the input field is only copied to the output array.
    implicit none
    integer,intent(in)::kmype
    integer,intent(in)::klon,klat
    real(kind=realkind),intent(in)::p_local(klon,klat)
    real(kind=realkind)::p_global(klon_global,klat_global)

#ifdef MPI_SRC
#include"mpif.h"
    integer::istat
#endif
    integer::istart,istop,jstart,jstop
    integer::i,j,pe,error,ibuf(8)
    integer,allocatable,dimension(:)::ibuf_all 
    real(kind=realkind),allocatable,dimension(:)::tmp,tmp_all
    integer::n,m

    if(nproc == 1)then
       p_global = p_local
    else
#ifdef MPI_SRC
       call jlimits(klon,klat,istart,istop,jstart,jstop) !local indecies!!
       n = maxval(klon_ar)
       m = maxval(klat_ar)
       allocate(tmp(n*m),tmp_all(n*m*nproc))  
       ibuf=(/jstart,jstop,jdatastart,klat,istart,istop,idatastart,klon/)

       allocate(ibuf_all(nproc*8))
       call MPI_gather(ibuf,8,MPI_INTEGER,ibuf_all,8,MPI_INTEGER,kmype,localComm,error)
       tmp=-666.0_realkind
       do j=1,klat
          do i=1,klon
             tmp(i+(j-1)*klon) = p_local(i,j)
          enddo
       enddo
       call MPI_Gather(tmp,n*m,REALTYPE,tmp_all,n*m,REALTYPE,kmype,localComm,error)

       if(mype==kmype)then
          do pe=0,nproc-1
             do j=1,8
                ibuf(j) = ibuf_all(j+pe*8)
             enddo
             do j=1,n*m
                tmp(j) = tmp_all(j+pe*n*m)
             enddo
             do j=ibuf(1),ibuf(2)
                do i=ibuf(5),ibuf(6)
                   p_global(ibuf(7)+i-1,ibuf(3)+j-1)= tmp(i + (j-1)*ibuf(8))
                enddo
             enddo
          enddo
       endif
       deallocate(tmp,tmp_all,ibuf_all) 
#endif
    endif
    return
  end subroutine colfld



  subroutine swapI(a,klon,klat)
    !  purpose:
    !     updates the halo zone from the neighbouring processors
    !  author:
    !      kalle eerola   hirlam    1998

    implicit none
    integer::klon,klat
    integer::a(klon,klat)
    integer::len_ns,len_ew

#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
    integer::info
#endif
    integer,allocatable,dimension(:)::buf_nss,buf_nsr
    integer,allocatable,dimension(:)::buf_ews,buf_ewr

    if(nproc==1)return
#ifdef MPI_SRC
    !     Alway make sure that send and receive will be pairs
    len_ns = klon
    len_ew = klat
    allocate(buf_nss(klon),buf_ews(klat))
    allocate(buf_nsr(klon),buf_ewr(klat))

    if(pe_top/=-1)then
       call mpi_recv(buf_nsr,len_ns,mpi_integer,pe_top,11051,localComm,status,info)
       a(:,klat)= buf_nsr
    endif
    if(pe_base/=-1)then
       buf_nss = a(:,2)
       call  mpi_send(buf_nss,len_ns,mpi_integer,pe_base,11051,LOCALCOMM,info)
       call mpi_recv(buf_nsr,len_ns,MPI_INTEGER,pe_base,11061,LOCALCOMM,STATUS,info)
       a(:,1)= buf_nsr
    endif
    if(pe_top/=-1)then
       buf_nss = a(:,klat-1)
       call mpi_send(buf_nss,len_ns,mpi_integer,pe_top,11061,localComm,info)
    endif


    if(pe_right/=-1)then
       call mpi_recv(buf_ewr,len_ew,mpi_integer,pe_right,11211,localComm,status,info)
       a(klon,:)= buf_ewr
    endif
    if(pe_left/=-1)then
       buf_ews = a(2,:)
       call mpi_send(buf_ews,len_ew,mpi_integer,pe_left,11211,localComm,info)
       call mpi_recv(buf_ewr,len_ew,mpi_integer,pe_left,11311,localComm,status,info)
       a(1,:)= buf_ewr
    endif
    if(pe_right/=-1)then
       buf_ews = a(klon-1,:)
       call mpi_send(buf_ews,len_ew,mpi_integer,pe_right,11311,localComm,info)
    endif

    deallocate(buf_nss,buf_ews)
    deallocate(buf_nsr,buf_ewr)

#ifdef USE_BARRIER
    call mpi_barrier(localComm,info)
#endif
#endif

    return
  end subroutine swapI


  subroutine swap_ps(a,klon,klat)
    !  purpose:
    !     exchange halo zone of surface pressure.
    !     this routine exchanges two rows and columns, which
    !     is needed for surface pressure because of staggering
    !     in 'dyn' and 'sldyn'.
    !  author:
    !      kalle eerola   hirlam    1998

    implicit none

    integer,intent(in):: klon,klat
    real(kind=realkind):: a(klon+1,klat+1)
    integer:: info 
#ifdef MPI_SRC
#include"mpif.h"
    integer::status(mpi_status_size)
#endif
    if(nproc==1) return

#ifdef MPI_SRC

#ifdef USE_BARRIER
    call mpi_barrier(localComm,info)
#endif


    if (.not.attop) then
       call mpi_recv(a(1,klat),1,jside2p,pe_top,1201,localComm,status,info)
    endif
    if (.not.atbase) then
       call mpi_send(a(1,2),1,jside2p,pe_base,1201,localComm,info)
    endif

    if (.not. atbase) then    
       call mpi_recv(a,1,jsidep,pe_base,1211,localComm,status,info)
    endif
    if (.not. attop) then
       call mpi_send(a(1,klat-1),1,jsidep,pe_top,1211,localComm,info)
    endif


    if (.not. atright) then 
       call mpi_recv(a(klon,1),1,iside2p,pe_right,1221,localComm,status,info)
    endif
    if (.not. atleft) then
       call mpi_send(a(2,1),1,iside2p,pe_left,1221,localComm,info)
    endif

    if (.not. atright) then 
       call mpi_send(a(klon-1,1),1,isidep,pe_right,1231,localComm,info)
    endif
    if (.not. atleft) then
       call mpi_recv(a,1,isidep,pe_left,1231,localComm,status,info)
    endif

#ifdef USE_BARRIER
    call mpi_barrier(localComm,info)
#endif
#endif
    return
  end subroutine swap_ps

  subroutine stop_program(message)
    implicit none
    character(len=*)::message
    integer::ierr
    print *,'stop_program called from mype=',mype
    print *,message
    !#ifdef MPI_SRC
    !    call MPI_Barrier(localComm,ierr)
    !    call MPI_Finalize(ierr)
    !#endif
    stop 
  end subroutine stop_program

end module decomp
