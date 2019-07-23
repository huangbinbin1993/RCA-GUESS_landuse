module config
  implicit none
  private
  logical,public,save:: use_oasis=.false.
  logical,public,save:: use_oasis_arctic=.false.
  logical,public,save:: use_oasis_netcdf_output=.false.
  logical,public,save:: use_guess=.false.
  logical,public,save:: use_routing=.false.
  logical,public,save:: use_match=.false.!for future use
  logical,public,save:: use_jarvisco2=.false.

  public read_config
contains
  subroutine read_config()
    implicit none
    namelist/namconfig/use_oasis,use_guess,use_routing,use_match,&
    use_oasis_arctic,use_oasis_netcdf_output,use_jarvisco2

    
    open(57,file='namelists.dat',status='old')
    read(57,nml=namconfig)
    close(57)

  end subroutine read_config

end module config
