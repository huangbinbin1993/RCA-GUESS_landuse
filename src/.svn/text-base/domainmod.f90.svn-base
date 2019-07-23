module domainmod
  use decomp,only:stop_program,realkind
  implicit none
  private
  type,public::domain
     real(kind=realkind)::polon
     real(kind=realkind)::polat
     real(kind=realkind)::south
     real(kind=realkind)::west
     real(kind=realkind)::dlon
     real(kind=realkind)::dlat
     !global sizes only!
     integer::nlon
     integer::nlat
     integer::nlev
     real(kind=realkind),allocatable,dimension(:,:)::fis
     real(kind=realkind),allocatable,dimension(:,:)::orogsigm
     real(kind=realkind),allocatable,dimension(:)::ahyb
     real(kind=realkind),allocatable,dimension(:)::bhyb
     real(kind=realkind),allocatable,dimension(:)::afull
     real(kind=realkind),allocatable,dimension(:)::bfull
     real(kind=realkind),allocatable,dimension(:)::hybi
     real(kind=realkind),allocatable,dimension(:)::hybk
     real(kind=realkind)::stagu_x
     real(kind=realkind)::stagu_y
     real(kind=realkind)::stagv_x
     real(kind=realkind)::stagv_y
     character(len=1)::arakawa='c'
  end type domain

  public arakawastag
contains

  subroutine arakawastag(arakawa, stag_u_x, stag_u_y, stag_v_x,stag_v_y)
    ! 971030 Ulf Hansson, Rossby Centre
    ! find staggered grid coefficients for deifferent Arakawa grids
    implicit none  
    character(len = 1),intent(in)::arakawa  
    real(kind=realkind),intent(out)::stag_u_x, stag_u_y, stag_v_x, stag_v_y  

    if (arakawa=='a'.or.arakawa=='A') then  
       stag_u_x = 0.0_realkind  
       stag_u_y = 0.0_realkind  
       stag_v_x = 0.0_realkind  
       stag_v_y = 0.0_realkind  
    elseif (arakawa=='b'.or.arakawa=='B') then  
       stag_u_x = 0.5_realkind  
       stag_u_y = 0.5_realkind  
       stag_v_x = 0.5_realkind  
       stag_v_y = 0.5_realkind  
    elseif (arakawa=='c'.or.arakawa=='C') then  
       stag_u_x = 0.5_realkind  
       stag_u_y = 0.0_realkind  
       stag_v_x = 0.0_realkind  
       stag_v_y = 0.5_realkind  
    else  
       write(6, * ) 'arakwastag called with illegal value: ', &
            arakawa
       write(6, * ) 'The only allowed values are a,b & c '
       write(6, * ) 'Have to stop'  
       call stop_program( 'arakawastag'  )

    endif
    return  
  end subroutine arakawastag


end module domainmod
