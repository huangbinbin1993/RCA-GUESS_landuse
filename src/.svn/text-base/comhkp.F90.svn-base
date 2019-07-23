module comhkp
  use decomp,only:realkind
  use calendar
  use timetype
  implicit none
  private

  logical,public,save::nlsimp=.true.          ! .true. if semi-implicit sheme
  logical,public,save::nlphys=.true.          ! .true. if physical paramerization
  logical,public,save::nlstat=.true.          ! .true. if computation and print of statics
  logical,public,save::nltvir =.true.         ! .true. if virtual temperature in dynamics
  logical,public,save::nlhumc =.true.         ! .true. if check of critical humidity for input data
  logical,public,save::nltcrf=.true.          !.true.if correctionfor horizontal diffusion of temperature and humidity along pseud-pressure levels   boundary relaxations constants

  logical,public,save::nldynvd=.false.        ! .true. if dynamic tendency used in the vertical diffusion scheme


  logical,save::nlusug=.false.         !

  integer,public,save,allocatable::nwmosv(:)

  type(time),public,save::current_time !current time
  type(deltatime),public,save::ndtime  !timestep

  real(kind=realkind),public,save   ::dtphys=-1.0_realkind           ! timestep for physics in seconds
  real(kind=realkind),public,save   ::dtvdif=-1.0_realkind           ! timestep for vertical diffusion in seconds
  real(kind=realkind),public,save   ::timesu=-1.0_realkind          ! spinup time in seconds (the physics is called at every dynamic timestep during spinup time)   geographical area etc


  logical,public,save::nlpost=.true.         ! .true. if postprocessing at current time step

  public read_namrun,init_hkp_memory
contains

  subroutine read_namrun(ksvar)
    use decomp
    use calendar
    use referenceParameters

    implicit none
    integer,parameter::nc=2
    integer,intent(in)::ksvar 
    integer::nwmosvin(nc)
    integer::dtime

    namelist / namrun / dtime,nlsimp,nlphys,nlstat,nltvir,nlhumc, &
         nltcrf ,   dtphys, dtvdif, timesu, nldynvd, &
            nwmosvin    , sit0   , sip0 , nlusug
    
    dtime = -360

    open(57,file='namelists.dat',status='old')
    read(57,nml=namrun)
    close(57)
    if(mype==0)then
       write(6,nml=namrun)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namrun)
       close(1)
    endif
    if(dtime<0)call stop_program( 'dtime<0 (error) or not given')
    ndtime = makedeltatime(0,0,0,0,0,dtime)
    if(dtphys<0.0_realkind)dtphys=real(dtime,realkind)
    if(dtvdif<0.0_realkind)dtvdif=real(dtime,realkind)
    if(timesu<0.0_realkind)timesu=real(2*dtime,realkind)
    
    nwmosv = nwmosvin(1:ksvar)


  end subroutine read_namrun


  subroutine init_hkp_memory(klev,ksvar)
    implicit none
    integer,intent(in)::klev,ksvar
    allocate(nwmosv(ksvar))
    
  end subroutine init_hkp_memory





  

end module comhkp
