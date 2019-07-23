module transposeJB
  implicit none
  private
  integer, parameter::SEND=1,RECV=2,PNORTH=1,PSOUTH=2,PEAST=3,PWEST=4
  integer, parameter :: POINT_TO_POINT = 1
  integer, parameter :: ALL_TO_ALL     = 2
  integer,save :: hl_comm_x, hl_comm_y, mp_method
  integer,allocatable::size_buf_fft(:),size_div_fft(:)
  integer,allocatable::size_buf_tri(:),size_div_tri(:)
  integer,allocatable::offs_fft_inp(:),offs_fft_self(:)
  integer,allocatable::offs_fft_peer(:),offs_tri_peer(:)
  integer,allocatable::offs_buf_fft(:),offs_buf_tri(:)
  integer,allocatable::offs_tri_inp(:),offs_tri_self(:)
  integer,save::nrec_fft,nbuf_fft,nrec_tri,nbuf_tri
!  integer,save::offs_swap(SEND:RECV,PNORTH:PWEST)
  real,public,pointer,save::div_fft(:,:)
  real,public,pointer,save::div_tri(:,:)
  real,pointer,save::buf_fft(:)
  real,pointer,save::buf_tri(:)
  integer,allocatable,save::nlonx(:), nlony(:)
  integer,allocatable,save::nlaty(:), nlevx(:)
  integer,save::imin_tri,imax_tri,jlev_tri

  public twod_to_fft,fft_to_tri,tri_to_fft,fft_to_twod,decompose_hh,allocate_transpose
contains
  subroutine allocate_transpose(klon_global, klat_global, klev_global,nproc)
    implicit none
    integer,intent(in)::klon_global, klat_global, klev_global,nproc
    allocate(size_buf_fft (nproc),size_div_fft (nproc), &
         size_buf_tri (nproc),size_div_tri (nproc), &
         offs_fft_inp (nproc),offs_fft_self(nproc), &
         offs_fft_peer(nproc),offs_tri_peer(nproc), & 
         offs_buf_fft (nproc),offs_buf_tri (nproc), &
         offs_tri_inp (nproc),offs_tri_self(nproc), & 
         nlonx(nproc), nlony(nproc), &
         nlaty(nproc), nlevx(nproc))

  end subroutine allocate_transpose

  subroutine twod_to_fft(div, klon, klat, klev)
    use decomp
    implicit none
#ifdef MPI_SRC
#include"mpif.h"
    integer :: stat(MPI_STATUS_SIZE,nprocx*2)
#endif    
    integer,intent(in)::klon,klat,klev
    real,intent(in):: div(klon,klat,klev)
    integer :: ierror, ipe, nreq, req(nprocx*2)
    integer :: i,imin,j,jmin,ibuf,jbuf,k
    integer :: nlonx_self, nlaty_self, nlevx_peer


    nlonx_self = nlonx(i_pe_grid+1)
    nlaty_self = nlaty(j_pe_grid+1)
    do ipe = 1, nprocx
       nlevx_peer = nlevx(ipe)
       imin = offs_fft_inp(ipe)
       jmin = offs_buf_fft(ipe)
       do i = 1, nlonx_self
          do k = 1, nlevx_peer
             jbuf = jmin + nlaty_self*( k-1 + (i-1)*nlevx_peer )
             do j = 1, nlaty_self
                buf_fft(jbuf+j) = div(imin+i,j,k)
             end do
          end do
       end do
    end do

#ifdef MPI_SRC
    select case( mp_method )

    case( POINT_TO_POINT )
       nreq = 0
       do ipe = 1, nprocx
          call mpi_isend(buf_fft(offs_buf_fft(ipe)+1),size_buf_fft(ipe),MPI_REAL,&
               ipe-1, 1000, hl_comm_x, req(nreq+1), ierror )
          call mpi_irecv(div_fft(offs_fft_self(ipe)+1,1),size_div_fft(ipe),&
               MPI_REAL,ipe-1,1000, hl_comm_x, req(nreq+2), ierror )
          nreq = nreq + 2
       end do
       call mpi_waitall( nreq, req, stat, ierror )

    case( ALL_TO_ALL )
       call mpi_alltoallv(buf_fft, size_buf_fft, offs_buf_fft,  MPI_REAL, &
            div_fft, size_div_fft, offs_fft_self, MPI_REAL, &
            hl_comm_x, ierror )

    end select
#endif

  end subroutine twod_to_fft



  subroutine decompose_hh(klon, klat, klev,kpbpts)
    use comtrix
    use decomp
    implicit none
    integer,intent(in) :: klon,klat,klev,kpbpts
#ifdef MPI_SRC
#include "mpif.h"
#endif

    integer ::i, imin, irec, ipe, ierr, j, jmin, jrec, jpe, k, rem,&
         n,  nlevx_self, nlevx_peer, nlonx_self, nlonx_peer,&
         nlony_self, nlony_peer,nlaty_self, nlaty_peer,ilonx_self, ilonx_peer,&
         jlaty_self, nlon, nlon_x, nlon_y, nlat, tp, tp1, tp2, tp3, &
         offset_fft, offset_tri, offs_inp,klon_slice, klev_slice,&
         dims(6), nbytes

#ifdef MPI_SRC
    character*80 :: env
    character*20 :: method_name(2)
#endif

    klev_hh = klev
    triwgt  = -1

    nlon_x = klon_global
    nlon_y = klon_global-2*kpbpts-3
    nlat   = klat_global

    ! Compute nlony: number of longitudes for each of the nprocy processors

    klon_slice  = nlon_y / nprocy
    rem = nlon_y - nprocy * klon_slice
    nlony(rem+1:nprocy) = klon_slice
    if( rem > 0 ) nlony(1:rem) = klon_slice + 1

    ! Compute nlevx: number of levels for each of the nprocx processors

    klev_slice = klev / nprocx
    rem = klev - nprocx * klev_slice
    nlevx(rem+1:nprocx) = klev_slice
    if( rem > 0 ) nlevx(1:rem) = klev_slice + 1

    ! Compute nlonx: number of longitudes for each of the nprocx processors

    nlonx(     1) = idatastart_ar(2) - idatastart_ar(1) + 1
    nlonx(nprocx) = klon_global      - idatastart_ar(nprocx)
    do ipe = 2, nprocx-1
       nlonx(ipe) = idatastart_ar(ipe+1) - idatastart_ar(ipe)
    end do

    ! Compute nlaty: number of latitudes for each of the nprocy processors

    nlaty(     1) = jdatastart_ar(2) - jdatastart_ar(1) + 1
    nlaty(nprocy) = klat_global      - jdatastart_ar(nprocy)
    do jpe = 2, nprocy-1
       nlaty(jpe) = jdatastart_ar(jpe+1) - jdatastart_ar(jpe)
    end do

    ! Common counts and offsets for FFT and TRI transpositions

    nlevx_self = nlevx( i_pe_grid+1 )
    nlony_self = nlony( j_pe_grid+1 )
    nlaty_self = nlaty( j_pe_grid+1 )
    nlonx_self = nlonx( i_pe_grid+1 )

    nrec_fft = nlevx_self * nlaty_self
    nrec_tri = nlevx_self * nlony_self

    imin = 1
    jmin = 1
    if( atleft  ) imin = 0
    if( atbase  ) jmin = 0
    ilonx_self = idatastart + imin - 1
    jlaty_self = jdatastart + jmin - 1

    ! Counts and offsets for div <--> fft transposition

    offs_inp = imin + klon * jmin
    nbuf_fft = 0
    do ipe = 1, nprocx
       nlevx_peer = nlevx( ipe )
       nlonx_peer = nlonx( ipe )
       ilonx_peer = idatastart_ar( ipe )
       if( ipe == 1 ) ilonx_peer = ilonx_peer - 1
       size_buf_fft (ipe) = nlonx_self * nlaty_self * nlevx_peer
       size_div_fft (ipe) = nlonx_peer * nlaty_self * nlevx_self
       offs_fft_peer(ipe) = ilonx_self * nlaty_self * nlevx_peer
       offs_fft_self(ipe) = ilonx_peer * nlaty_self * nlevx_self
       offs_fft_inp (ipe) = offs_inp
       offs_buf_fft (ipe) = nbuf_fft
       offs_inp = offs_inp + klon * klat * nlevx_peer
       nbuf_fft = nbuf_fft + size_buf_fft(ipe)
    end do

    ! Counts and offsets for fft <--> tri transposition

    nbuf_tri   = 0
    offset_tri = 0
    offset_fft = (kpbpts+1) * nrec_fft
    do jpe = 1, nprocy
       nlony_peer = nlony( jpe )
       nlaty_peer = nlaty( jpe )
       size_buf_tri (jpe) = nlony_peer * nlaty_self * nlevx_self
       size_div_tri (jpe) = nlony_self * nlaty_peer * nlevx_self
       offs_tri_peer(jpe) = nlony_peer * jlaty_self * nlevx_self
       offs_tri_self(jpe) = offset_tri
       offs_tri_inp (jpe) = offset_fft
       offs_buf_tri (jpe) = nbuf_tri
       nbuf_tri   = nbuf_tri   + size_buf_tri(jpe)
       offset_fft = offset_fft + size_buf_tri(jpe)
       offset_tri = offset_tri + size_div_tri(jpe)
    end do

    !  TRI range of logitudinal points

    imin_tri = kpbpts+2
    do jpe = 1, j_pe_grid
       imin_tri = imin_tri + nlony(jpe)
    end do
    imax_tri = imin_tri + nlony_self - 1

    !  TRI Levels offset

    jlev_tri = 0
    do ipe = 1, i_pe_grid
       jlev_tri = jlev_tri + nlevx(ipe)
    end do

    ! Halo swap buffer size and segment offsets

    nlon = maxval( nlonx(1:nprocx) )
    nlat = maxval( nlaty(1:nprocy) )


#ifdef MPI_SRC
    !   Create communicators for processors handling
    !   subgrids in X and Y direction
    call mpi_comm_split( localComm, j_pe_grid, i_pe_grid,hl_comm_x, ierr )
    call mpi_comm_split( localComm, i_pe_grid, j_pe_grid,hl_comm_y, ierr )

    mp_method = POINT_TO_POINT

    ! Overwrite default with env var setting

    !      call getenv( 'MP_METHOD', env )
!!$      if( env /= ' ' )then
!!$         mp_method = 0
!!$         if(      index(env,'POINT_TO_POINT') > 0 .or.
!!$     +            index(env,'point_to_point') > 0 .or.
!!$     +            index(env,'PTP'            ) > 0 .or.
!!$     +            index(env,'ptp'            ) > 0      )then
!!$
!!$            mp_method = POINT_TO_POINT
!!$
!!$         else if( index(env,'ALL_TO_ALL'    ) > 0 .or.
!!$     +            index(env,'all_to_all'    ) > 0 .or.
!!$     +            index(env,'ATA'            ) > 0 .or.
!!$     +            index(env,'ata'            ) > 0      )then
!!$
!!$            mp_method = ALL_TO_ALL
!!$         else
!!$            if( mype == 0 )then
!!$               print *,'ERROR in DECOMPOSE_HH: invalid value ',
!!$     +                 'of env var MP_METHOD: ', env
!!$               print *,'   Valid MP_METHOD settings are ',
!!$     +                 'POINT_TO_POINT (or PTP), ',
!!$     +                 'ALL_TO_ALL (or ATA), all uppercase or lowercase'
!!$            end if
!!$            call abort
!!$         endif
!!$      end if

    if( mype == 0 )then

       ! Report back the method used

       method_name(1) = 'Point-to-Point'
       method_name(2) = 'All-to-All'
       print '(/,2a,/)','Message passing method: ',trim(method_name(mp_method))
    end if
#endif
    allocate(div_fft (nrec_fft,klon_global),buf_fft (nbuf_fft),&
         div_tri (nrec_tri,klat_global),buf_tri (nbuf_tri))

    ! Preparations for HHSOLV

!    allocate( hhdia(nrec_tri,klat_global) )
!    dtold = 0.

  end subroutine decompose_hh

  SUBROUTINE FFT_TO_TRI( KPBPTS, KLEV )
    use decomp
    IMPLICIT NONE
    INTEGER KPBPTS,KLEV
    integer :: ierror, jpe, nreq, req(nprocy*2)
    integer :: i,imin,ibuf,j,jmin,jfft,k,m
    integer :: nlony_peer, nlaty_self, nlevx_self

#ifdef MPI_SRC
#include"mpif.h"
    integer :: stat(MPI_STATUS_SIZE,nprocy*2)
#endif

    nlaty_self = nlaty(j_pe_grid+1)
    nlevx_self = nlevx(i_pe_grid+1)
    do jpe = 1, nprocy
       imin = offs_buf_tri(jpe)
       jmin = offs_tri_inp(jpe)
       nlony_peer = nlony(jpe)
       do k = 1, nlevx_self
          jfft = jmin + (k-1) * nlaty_self
          do j = 1, nlaty_self
             ibuf = imin + nlony_peer*( k-1 + (j-1)*nlevx_self )
             do i = 1, nlony_peer
                buf_tri(ibuf+i) = div_fft(jfft+j,i)
             end do
          end do
       end do
    end do

#ifdef MPI_SRC
    select case( mp_method )
    case( POINT_TO_POINT )
       nreq = 0
       do jpe = 1, nprocy
          call mpi_isend(buf_tri(offs_buf_tri(jpe)+1),size_buf_tri(jpe),MPI_REAL, &
               jpe-1, 1000, hl_comm_y, req(nreq+1), ierror )

          call mpi_irecv(div_tri(offs_tri_self(jpe)+1,1),size_div_tri(jpe),MPI_REAL, &
               jpe-1,1000, hl_comm_y, req(nreq+2), ierror )
          nreq = nreq + 2
       end do
       call mpi_waitall( nreq, req, stat, ierror )
    case( ALL_TO_ALL )
       call mpi_alltoallv(buf_tri, size_buf_tri, offs_buf_tri,  MPI_REAL, &
            div_tri, size_div_tri, offs_tri_self, MPI_REAL, &
            hl_comm_y, ierror )
    end select
#endif
  end subroutine fft_to_tri

  SUBROUTINE TRI_TO_FFT( KPBPTS, KLEV )
    use decomp
    IMPLICIT NONE
    INTEGER KPBPTS,KLEV
    integer :: ierror, jpe, nreq, req(nprocy*2)
    integer :: i,imin,ibuf,j,jmin,jfft,k,m
    integer :: nlony_peer, nlaty_self, nlevx_self

#ifdef MPI_SRC
#include"mpif.h"
    integer :: stat(MPI_STATUS_SIZE,nprocy*2)
#endif

#ifdef MPI_SRC
    select case( mp_method )
    case( POINT_TO_POINT )
       nreq = 0
       do jpe = 1, nprocy
          call mpi_irecv( buf_tri(offs_buf_tri(jpe)+1),size_buf_tri(jpe),MPI_REAL, &
               jpe-1, 1000, hl_comm_y, req(nreq+1), ierror )

          call mpi_isend(div_tri(offs_tri_self(jpe)+1,1),size_div_tri(jpe),MPI_REAL, &
               jpe-1,1000, hl_comm_y, req(nreq+2), ierror )
          nreq = nreq + 2
       end do
       call mpi_waitall( nreq, req, stat, ierror )
    case( ALL_TO_ALL )
       call mpi_alltoallv(div_tri, size_div_tri, offs_tri_self, MPI_REAL, &
            buf_tri, size_buf_tri, offs_buf_tri,  MPI_REAL, hl_comm_y, ierror )
    end select
#endif

    nlaty_self = nlaty(j_pe_grid+1)
    nlevx_self = nlevx(i_pe_grid+1)
    do jpe = 1, nprocy
       imin = offs_buf_tri(jpe)
       jmin = offs_tri_inp(jpe)
       nlony_peer = nlony(jpe)
       do k = 1, nlevx_self
          jfft = jmin + (k-1) * nlaty_self
          do j = 1, nlaty_self
             ibuf = imin + nlony_peer*( k-1 + (j-1)*nlevx_self )
             do i = 1, nlony_peer
                div_fft(jfft+j,i) = buf_tri(ibuf+i)
             end do
          end do
       end do
    end do

  END SUBROUTINE TRI_TO_FFT

  SUBROUTINE FFT_TO_TWOD(DIV, KLON, KLAT, KLEV)
    use decomp
    IMPLICIT NONE
    INTEGER :: KLON,KLAT,KLEV
    REAL    :: DIV(KLON,KLAT,KLEV)
    integer :: ierror, ipe, nreq, req(nprocx*2)
    integer :: i,imin,j,jmin,ibuf,jbuf,k
    integer :: nlonx_self,nlaty_self,nlevx_peer

#ifdef MPI_SRC
#include"mpif.h"
    integer :: stat(MPI_STATUS_SIZE,nprocx*2)
#endif

    !-------------------------------------------------------------------

#ifdef MPI_SRC
    select case( mp_method )
    case( POINT_TO_POINT )
       nreq = 0
       do ipe = 1, nprocx
          call mpi_irecv(buf_fft(offs_buf_fft(ipe)+1),size_buf_fft(ipe),MPI_REAL, &
               ipe-1, 1000, hl_comm_x, req(nreq+1), ierror )
          call mpi_isend(div_fft(offs_fft_self(ipe)+1,1),size_div_fft(ipe), &
               MPI_REAL,ipe-1,1000, hl_comm_x, req(nreq+2), ierror )
          nreq = nreq + 2
       end do

       call mpi_waitall( nreq, req, stat, ierror )

    case( ALL_TO_ALL )

       call mpi_alltoallv(div_fft, size_div_fft, offs_fft_self, MPI_REAL, &
            buf_fft, size_buf_fft, offs_buf_fft,  MPI_REAL,  hl_comm_x, ierror )

    end select

#endif


    nlonx_self = nlonx(i_pe_grid+1)
    nlaty_self = nlaty(j_pe_grid+1)

    do ipe = 1, nprocx
       nlevx_peer = nlevx(ipe)
       imin = offs_fft_inp(ipe)
       jmin = offs_buf_fft(ipe)
       do i = 1, nlonx_self
          do k = 1, nlevx_peer
             jbuf = jmin + nlaty_self * ( k-1 + (i-1)*nlevx_peer )
             do j = 1, nlaty_self
                div(imin+i,j,k) = buf_fft(jbuf+j)
             end do
          end do
       end do
    end do

  END SUBROUTINE FFT_TO_TWOD


end module transposeJB
