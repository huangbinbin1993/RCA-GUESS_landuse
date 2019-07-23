module co2mod
  use timetype
  use decomp,only:stop_program,realkind
  implicit none
  private

  real(kind=realkind), public,save::const_co2,const_co2only
  logical,save::lsres=.false.
  logical,save::lsresa1b=.false.
  logical,save::lsresa2 =.false.
  logical,save::lsresb2 =.false.
  logical,save::lsresb1 =.false.
  logical,save::lrcp =.false.
  logical,save::lrcp26 =.false.
  logical,save::lrcp45 =.false.
  logical,save::lrcp85 =.false.


  public gco2,read_scenario
contains
  subroutine read_scenario()
    use decomp
    implicit none
#ifdef MPI_SRC
#include"mpif.h"
    integer ierr,sum
    logical lbuf(9)
#endif
    namelist/scenario/lsres,lsresa1b,lsresa2,lsresb2,lsresb1,lrcp,lrcp26,lrcp45,lrcp85
    if(mype==0)then
       open(57,file='namelists.dat',status='old')
       read(57,nml=scenario)
       close(57)
       write(6,nml=scenario)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=scenario)
       close(1)
    endif
!
! lsres & lrcp are not scenarios
! They indicate what to be used in historical runs
! Old values (SRES) or new values (RCP)
! They are almost identical, so it does not really matter
! Default is RCP
!
    if(mype==0) then
       sum=0
       if(lsresa1b)sum = sum+1
       if(lsresa2)sum = sum+1
       if(lsresb2)sum = sum+1
       if(lrcp26)sum = sum+1
       if(lrcp45)sum = sum+1
       if(lrcp85)sum = sum+1
       if(sum==0)then
          print *,' '
          print *,'You have not choosen a scenario'
          print *,'This only makes sense if it is a historical run'
          print *,' '
       elseif(sum/=1)then
          print *,'Only one scenario can be selected at a time!!!'
          print *,'You have choosen ',sum
          print *,'lrcp26,lrcp45,lrcp85,lsresa1b,lsresa2,lsresb2=',&
                   lrcp26,lrcp45,lrcp85,lsresa1b,lsresa2,lsresb2
          print *,'Have to stop'
          call stop_program('')
       endif
       if(lsres.and.lrcp) then
          print*,'You have chosen both lsres and lrcp, it does not make sense'
          print*,'They decide what GHG values to use in historical runs'
          print*,'New values (RCP) or old values (SRES)'
          print*,'Have to stop'
          call stop_program('')
       endif
       if(.not.lsres.and..not.lrcp) then
          print *,' '
          print *,'You have not choosen lsres or lrcp!'
          print *,'This only makes sense if it is a simulation with true'
          print *,'boundary conditions (ERA40 or ERAInterim) or a scenario run.'
          print *,'For a scenario they decide what GHG values to use in historical'
          print *,'runs. New values (RCP) or old values (SRES)'
          print *,' '
       endif
    endif   ! mype
    
#ifdef MPI_SRC
    lbuf = (/lsres,lsresa1b,lsresa2,lsresb2,lsresb1,lrcp,lrcp45,lrcp85,lrcp26/) 
    call MPI_Bcast(lbuf,9,MPI_LOGICAL,0,LOCALCOMM,ierr)
    lsres=lbuf(1)
    lsresa1b=lbuf(2)
    lsresa2=lbuf(3)
    lsresb2=lbuf(4)
    lsresb1=lbuf(5)
    lrcp=lbuf(6)
    lrcp45=lbuf(7)
    lrcp85=lbuf(8)
    lrcp26=lbuf(9)
#endif      
  end subroutine read_scenario

  subroutine gco2(ntime,nlinit)
    use gcm
    use calendar
    use decomp
    implicit none
    type(time),intent(in)::ntime 
    logical,intent(in)::nlinit 
    logical,save::lfirst=.true.

    if(lfirst)then
       call read_scenario
       lfirst=.false.
    endif

    if(lecmwf)then
       call era_co2(mype,ntime,nlinit)
    elseif(lsres.or.lsresa1b.or.lsresa2.or.lsresb1.or.lsresb2)then
       if(lecham2.or.lhadley)then
          if(mype.eq.0) print *, ' call gcm_co2'
          call gcm_co2(mype,ntime)
          if(mype.eq.0) print *, ' call gcm_co2 ok'
       elseif(lccsm)then
          call ccsm_co2(mype,ntime)
       endif
    elseif(lrcp.or.lrcp26.or.lrcp45.or.lrcp85)then
       call get_co2(mype,ntime,nlinit)
    else
       if(mype==0) then
          print *,'Scenario or historical co2 not specified'
          print *,' '
          if(.not.(lrcp.or.lrcp26.or.lrcp45.or.lrcp85)) then
             print *,'lrcp lrcp26 lrcp45 lrcp85 ',lrcp,lrcp26,lrcp45,lrcp85
          endif
          if(.not.(lsres.or.lsresa1b.or.lsresa2.or.lsresb1.or.lsresb2)) then
             print *,'lsres,lsresa1b,lsresa2,lsresb1,lsresb2 ',lsres,lsresa1b,lsresa2,lsresb1,lsresb2
          endif
          print *,' '
          print *,'Have to stop'
          print *,__FILE__,__LINE__
          print *,' '
          call stop_program('')
       endif
    endif

  end subroutine gco2

  subroutine get_co2(mype,intime,nlinit)
    !
    ! read CO2-equivalent from file
    ! file depends on scenario
    !
    use calendar
    use gcm
    implicit none

    integer,intent(in)::mype
    type(time),intent(in)::intime
    logical,intent(in)::nlinit

    integer:: i,year_co2,ilen,noheaderlines
    integer,save:: lunghg
    logical,save::lfirst=.true.
    character(len=132),save:: filename
    character(len=132):: header,prefix

    real(kind=realkind):: zco2,zdummy,zco2only


    if(.not.nlinit)then
       if(intime%month/=1 .or. intime%day/=1 .or. intime%hour/=0 .or. intime%min/=0) return
    endif

    !cau1106
    const_co2=0.
    const_co2only=0.

    if(lfirst)then
       call readghgpath(prefix,ilen)
       if(lrcp)then
          filename=prefix(1:ilen)//'/PRE2005_MIDYR_CONC.DAT'
       elseif(lrcp26)then
          filename=prefix(1:ilen)//'/RCP3PD_MIDYR_CONC.DAT'
       elseif(lrcp45)then
          filename=prefix(1:ilen)//'/RCP45_MIDYR_CONC.DAT'
       elseif(lrcp85)then
          filename=prefix(1:ilen)//'/RCP85_MIDYR_CONC.DAT'
       else
          if(mype==0) then
             print *,'Scenario not specified'
             print *,'lrcp lrcp26 lrcp45 lrcp85 ',lrcp,lrcp26,lrcp45,lrcp85
             print *,'Have to stop'
             print *,' '
             call stop_program('')
          endif
       endif

       lunghg=65
       open(lunghg,file=filename)
       noheaderlines=39
       do i=1,noheaderlines
          read(lunghg,'(a)') header
          if(mype==0) print *, header
       enddo
       if(mype==0) print *,' '

       year_co2=0
       do while  (year_co2 < intime%year)
          read(lunghg,*) year_co2,zco2,zdummy,zco2only
       enddo

       lfirst=.false.

    else

       read(lunghg,*) year_co2,zco2,zdummy,zco2only
       if(year_co2/=intime%year .and. mype==0)then
          print *,'The year we read ',year_co2,' is not equal to actual year',intime%year
          print *, 'the co2 value was read from file ',filename
          print *,'Have to stop'
          print *,' '
          call stop_program('')
       endif

    endif

    if(mype==0) print *,'get_co2 read ', year_co2,zco2
    const_co2=zco2
    const_co2only=zco2only

    ! interpolate between months?
    ! interpolate between years, not any more
    !          if ( year_true >= iy(i) .and. year_true < iy(i+1)) then
    !             ss = (real(year_true,realkind)-real(iy(i),realkind))/10.0_realkind
    !             const_co2= ss * zco2(i+1) + (1.0_realkind-ss)*zco2(i)
    !          endif

    if (mype == 0) then
       print *, ' '
       print *, 'set co2 equivalent ',intime%year,const_co2
       print *, 'the co2 value was read from file ',filename
       print *, ' '
       if(abs(const_co2)< 1.e-14_realkind) then
          print *,'ERROR in gcm_co2, co2 value is zero'
          print *,'maybe the year is outside the range'
          print *,'have to stop'
          print *,' '
          call stop_program('')
       endif
    endif

    return
  end subroutine get_co2

  subroutine ccsm_co2(mype,intime)
    use calendar
    implicit none
    integer,intent(in):: mype
    type(time),intent(in)::intime
    integer :: iy(15),i,year_true
    real(kind=realkind)::     ss, zco2(15)
    data iy/1960,1970,1980,1990,2000,2010,2020,2030,2040, &
         2050,2060,2070,2080,2090,2100/

    const_co2only=380.0_realkind

    !     A1B prel from CCSM output
    if(lsresa1b)then
       zco2=(/317.0_realkind,325.0_realkind,338.0_realkind,354.0_realkind,371.0_realkind,&
            395.0_realkind,427.0_realkind,462.0_realkind,498.0_realkind,&
            535.0_realkind,570.0_realkind,604.0_realkind,635.0_realkind,664.0_realkind,&
            675.0_realkind/)
    endif
    if(lsresa2)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            403.0_realkind,426.0_realkind,470.0_realkind,532.0_realkind, &
            602.0_realkind,702.0_realkind,823.0_realkind,963.0_realkind,1123.0_realkind,&
            1316.0_realkind/)
    endif
    if(lsresb2)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            409.0_realkind,453.0_realkind,492.0_realkind,536.0_realkind,&
            581.0_realkind,628.0_realkind,678.0_realkind,730.0_realkind,787.0_realkind,&
            847.0_realkind/)
    endif

    year_true = intime%year
    if ( year_true < 1960 ) then
       const_co2= zco2(1)
    else if ( year_true > 2100 ) then
       const_co2= zco2(15)
    else
       do i = 1,14
          if(year_true >= iy(i) .and. year_true < iy(i+1)) then
             ss = real((year_true-iy(i)),realkind)/10.0_realkind
             const_co2= ss * zco2(i+1) + (1.0_realkind-ss)*zco2(i)
          endif
       enddo
    endif

    if (mype == 0) then
       print *, ' '
       print *, 'set co2 equivalent for CCSM ',intime%year,const_co2
       if(lsresa2)then
          print *, 'the co2 value is from the SRES A2 scenario'
       endif
       if(lsresb2)then
          print *, 'the co2 value is from the SRES B2 scenario'
       endif
       if(lsresa1b)then
          print *, 'the co2 value is from the SRES A1B scenario'
       endif
       print *, ' '
       if(abs(const_co2)< 1.e-14_realkind) then
          print *,'ERROR in gcm_co2, co2 value is zero'
          print *,'maybe the year is outside the range'
          print *,'or'
          print *,'maybe you forgot to define scenario in Makefile'
          print *,'have to stop'
          print *,' '
          call stop_program('')
       endif
    endif

    return
  end subroutine ccsm_co2

  subroutine era_co2(mype,intime,nlinit)
    use calendar
    implicit none

    logical,intent(in)::nlinit
    integer,intent(in)::mype
    type(time),intent(in)::intime
 
    integer,save::iyref=1990
    real(kind=realkind),save::zco2ref=353.0_realkind

    
    integer:: iyr1,imn,irfy

    if(.not.nlinit)then
       if(intime%day/=1 .or. intime%hour/=0 .or. intime%min/=0) return
    endif

    if (intime%month<=6) then
       iyr1=intime%year-1
       imn=intime%month+6
    elseif (intime%month>6) then
       iyr1=intime%year
       imn=intime%month-6
    endif

    irfy=iyr1-iyref

    const_co2=zco2ref+real(irfy,realkind)*1.5_realkind+real(imn,realkind)*0.125_realkind
    const_co2only=const_co2

    if (mype == 0) then
       print *, ' '
       print *, 'set co2 equivalent ',intime%year,intime%month,const_co2
       print *, 'the co2 value was set in subroutine era_co2'
       print *, ' '
    endif

    return
  end subroutine era_co2


  subroutine gcm_co2(mype,intime)
    use calendar
    use gcm
    implicit none

    integer,intent(in)::mype
    type(time),intent(in)::intime
    integer:: iy(15),i,year_true
    real(kind=realkind):: ss,zco2(15)
    data iy/1960,1970,1980,1990,2000,2010,2020,2030,2040, &
         2050,2060,2070,2080,2090,2100/
    if(lsresa2)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            403.0_realkind,426.0_realkind,470.0_realkind,532.0_realkind, &
            602.0_realkind,702.0_realkind,823.0_realkind,963.0_realkind,1123.0_realkind,&
            1316.0_realkind/)
    elseif(lsresb2)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            409.0_realkind,453.0_realkind,492.0_realkind,536.0_realkind,&
            581.0_realkind,628.0_realkind,678.0_realkind,730.0_realkind,787.0_realkind,&
            847.0_realkind/)
    elseif(lsresb1)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            402.0_realkind,435.0_realkind,470.0_realkind,504.0_realkind,&
            540.0_realkind,576.0_realkind,606.0_realkind,625.0_realkind,636.0_realkind,&
            637.0_realkind/)
    elseif(lsresa1b)then
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            396.0_realkind,436.0_realkind,495.0_realkind,572.0_realkind,&
            634.0_realkind,713.0_realkind,781.0_realkind,832.0_realkind,871.0_realkind,&
            902.0_realkind/)
    else
       zco2=(/313.0_realkind,314.0_realkind,331.0_realkind,353.0_realkind,373.0_realkind,&
            -1.0_realkind,-1.0_realkind,-1.0_realkind,-1.0_realkind,-1.0_realkind,&
            -1.0_realkind,-1.0_realkind,-1.0_realkind,-1.0_realkind,-1.0_realkind/)
    endif

!cau1106
    const_co2=0.0_realkind
    const_co2only=380.0_realkind

    ! i.e. ECHAM4
    if(lecham2.and..not.lecham5)then
       year_true = intime%year + 1760
    else
       year_true = intime%year
    endif
!cau1106
    if(lhadley) then
       year_true = intime%year
    endif
    if (mype.eq.0) print *, ' gcm_co2 ',year_true

    if ( year_true < 1960 ) then
       const_co2= zco2(1)
    else if ( year_true > 2100 ) then
       const_co2= zco2(15)
    else
       do i = 1,14
          if ( year_true >= iy(i) .and. year_true < iy(i+1)) then
             ss = (real(year_true,realkind)-real(iy(i),realkind))/10.0_realkind
             const_co2= ss * zco2(i+1) + (1.0_realkind-ss)*zco2(i)
          endif
       enddo
    endif

    if (mype == 0) then
       print *, ' '
       print *, 'set co2 equivalent for GCM ',intime%year,const_co2
       if(lsresa2)then
          print *, 'the co2 value is from the SRES A2 scenario'
       endif
       if(lsresb2)then
          print *, 'the co2 value is from the SRES B2 scenario'
       endif
       if(lsresa1b)then
          print *, 'the co2 value is from the SRES A1B scenario'
       endif
       if(lsresb1)then
          print *, 'the co2 value is from the SRES B1 scenario'
       endif
       print *, ' '
       if(abs(const_co2)< 1.e-14_realkind) then
          print *,'ERROR in gcm_co2, co2 value is zero'
          print *,'maybe the year is outside the range'
          print *,'or'
          print *,'maybe you forgot to define scenario in Makefile'
          print *,'have to stop'
          print *,' '
          call stop_program('')
       endif
    endif

    return
  end subroutine gcm_co2

  subroutine readghgpath(ghg_path,ilen)
    use decomp,only:mype
    implicit none
    character(len=132),intent(out)::ghg_path
    integer,intent(out)::ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /ghg/ghg_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=ghg)
    close(61)
    if(mype==0.and..not.done)then
       write(6,nml=ghg)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=ghg)
       close(1)
       done = .true.
    endif
    ilen=index(ghg_path,' ')-1
  end subroutine readghgpath

end module co2mod
