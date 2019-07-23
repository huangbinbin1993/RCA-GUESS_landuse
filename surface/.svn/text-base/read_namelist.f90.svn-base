subroutine getgtopo30path(gtopo30_path)
  implicit none
  character(len=128)::gtopo30_path
  character(len=20)namelistfile,inst
  namelist/institute/inst
  namelist/gtopo30/gtopo30_path 
  
  open(66,file='namelists.dat',status='old',form='formatted')
  read(66,nml=institute)
  close(66)
  namelistfile='gcmpaths.'//trim(inst)
  open(61,file=namelistfile,status='old',form='formatted')
  read(61,nml=gtopo30)
  close(61)
end subroutine getgtopo30path


subroutine read_domain(south,west,dlon,dlat,polon,polat, &
     klon_global,klat_global,klev_global)
  use decomp, only:mype,realkind
  implicit none
  real(kind=8)::south,west,dlon,dlat,polon,polat
  integer::klon_global,klat_global,klev_global

  namelist/domain/south,west,dlon,dlat,polon,polat,&
       klon_global,klat_global,klev_global


  open(66,file='namelists.dat',status='old',form='formatted')
  read(66,nml=domain)
  close(66)
end subroutine read_domain


subroutine getecopath(ecopath)
  implicit none
  character(len=128)::ecopath
  character(len=20)namelistfile,inst
  namelist/institute/inst
  namelist/ecoclimap/ecopath
  open(66,file='namelists.dat',status='old',form='formatted')
  read(66,nml=institute)
  close(66)
  namelistfile='gcmpaths.'//trim(inst)
  open(61,file=namelistfile,status='old',form='formatted')
  read(61,nml=ecoclimap)
  close(61)
end subroutine getecopath


