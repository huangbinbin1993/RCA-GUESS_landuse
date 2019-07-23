subroutine hlprog
  !     HLPROG - MAIN PROGRAM

  use comhkp
  use decomp
  use gemini__genmod
  use config
  use coupling
  use transpose
  use sl2tim4
#ifdef USING_NETCDF
  use netcdfrestart
#else
  use restart
#endif
  use boundaryRelaxation
  use gcm

  implicit none
#if defined(MPI_SRC)
#include"mpif.h"
  integer ierr
  integer ver, subver
#endif


  integer::klat               !
  integer::klon               ! no of points in x-dirctn in a subarea
  integer::klev               ! no of points in y-dirctn in a subarea
  integer::MPICOMP


#ifdef TIMING
  real::zt1, zt2,tim,timL
#endif

  mype = 0
  nproc = 1
  call read_config()
  model_name = 'RCA'


  if(use_oasis)then

#ifndef MPI_SRC
     stop 'using oasis without MPI N/A'
#endif

#ifdef MPI_SRC
#ifdef OASIS3
     call couple_prism_init(localcomm)
     write(6,*) 'after prism init',localComm
!     call flush(6)
     call mpi_comm_rank(localcomm,mype,ierr)
     call mpi_comm_size(localcomm,nproc,ierr)
     call mpi_get_version( ver, subver,ierr )
     if( mype==0)then
        write(*,*)"Initializing MPI Version ", ver,subver
     endif
     write (*,*) 'I am the ', trim(model_name), &
          ' comp ', comp_id,                            &
          ' local rank ',mype
#endif
  else
     call mpi_init(ierr)
     if(ierr/=MPI_SUCCESS)stop 'MPI_INIT_ERROR'
     call mpi_comm_dup(mpi_comm_world,localcomm,ierr)
     if(ierr/=MPI_SUCCESS)stop 'MPI_COMM_DUP_ERROR'
     call mpi_comm_rank(localcomm,mype,ierr)
     if(ierr/=MPI_SUCCESS)stop 'MPI_COMM_RANK_ERROR'
     call mpi_comm_size(localcomm,nproc,ierr)
     if(ierr/=MPI_SUCCESS)stop 'MPI_COMM_SIZE_ERROR'
     call mpi_get_version( ver, subver,ierr )
     if(ierr/=MPI_SUCCESS)stop 'MPI_COMM_GET_VERSION_ERROR'

     call MPI_Comm_Compare(localComm,MPI_COMM_WORLD,MPICOMP,ierr)
     if(mype==0)then
        if(MPICOMP==MPI_IDENT)THEN
           print *,mype,'::identical'
        elseif(MPICOMP==MPI_CONGRUENT)then
           print *,mype,'MPI_CONGRUENT'
     call flush(6)
        elseif(MPICOMP==MPI_SIMILAR)then
           print *,mype,'similar'
        elseif(MPICOMP==MPI_UNEQUAL)then
           print *,mype,'unequal'
        else
           stop 'N/A'
        endif
     endif
     if( mype==0)then
        write(*,*)"Initializing MPI Version ", ver,subver
     endif
#endif

  endif


  if(mype==0)then
     open(1,file='printed_namelists.dat',status='replace')
     write(1,*)'Printed namelists:'
     close(1)
  endif

  call readDeconam()

  call init_hkp_memory(klev_global,1)
  call read_namrun(1)

  call read_namsl()
  call read_restart()

  call decompose(klon,klat,klev)
  call allocate_transpose(klon_global, klat_global, klev_global,nproc)
  call decompose_hh(klev_global, npbpts)  
  !  call initRca(klon,klat,klev)
  !cau  call getgcm(startTime,klon,klat,RCAdom) move to gemini



  if(mype==0)print *,'start forecast'
#ifdef TIMING
  call second (zt1)
#endif

  call gemini(klon,klat,klev)!,ksvar,iacdg,iacdg2) 

#ifdef TIMING
  call second (zt2)
  tim = zt2-zt1
#ifdef MPI_SRC
  timL=tim
  call mpi_reduce(rminl,tim,1,MPI_REAL,mpi_max,0,localcomm,ierr)
#endif  
  if(mype==0)print *,'forecast took ',tim,'seconds'
#endif

  !     EXIT PARALLEL SYSTEM
  if(use_oasis)then

#if defined(OASIS3)
        CALL prism_terminate_proto (ierr)
#endif


#ifdef MPI_SRC
  else
     call mpi_finalize(ierr)
#endif
  endif
  !  write(6,*) 'Atmosphere done!'

  return
end subroutine hlprog


