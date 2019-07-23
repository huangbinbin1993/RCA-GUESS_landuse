module ctun
  use decomp
  implicit none
  private
  real(kind=realkind),public,save::acrit = 1.00_realkind
  real(kind=realkind),public,save::tseafr= 271.15_realkind

  public contun
contains
  subroutine contun()
    !    subroutine contun:   stefan gollvik feb 1991
    !    give values to physical tuning constant

    implicit none

    namelist/namtun/ acrit,tseafr

    open(57,file='namelists.dat',status='old')
    read(57,nml=namtun)
    close(57)
    if(mype==0) then
       write(6,nml=namtun)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namtun)
       close(1)
    endif

    return
  end subroutine contun


end module ctun
