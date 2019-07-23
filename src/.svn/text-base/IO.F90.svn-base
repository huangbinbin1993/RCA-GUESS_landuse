module IO
  use decomp,only:realkind
  use mod_grib,only:realgribkind
  implicit none
  private

  public grrdloc_global 
contains

  subroutine grrdloc_global(lun,type,param,alev,fld, &  
       nx_grib,ny_grib,ierr)
    use decomp
    use grw1
    implicit none
    integer lun,type,param, nx_grib,ny_grib
    real(kind=realgribkind) alev
    real(kind=realkind)::fld(nx_grib,ny_grib)

    integer ierr

    ierr = 0
    if( mype==0 ) then
       call gread(lun,type,param,alev,fld,nx_grib,ny_grib,ierr)
    endif
    if(ierr/=0) return
    return
  end subroutine grrdloc_global

end module IO
