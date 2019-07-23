module calendar
  use decomp
  use timetype
  implicit none
  private


!BCM         real
!CNRM        real
!ECHAM5      real
!MPIESM      real
!EC-EARTH    real
!ERA40       real
!ERA-Interim real     
!NORESM      365 days, i.e. no leap years 31,28,31,30,31,30,31,31,30,31,30,31
!CANESM2     365 days, i.e. no leap years 31,28,31,30,31,30,31,31,30,31,30,31
!CCSM3       365 days, i.e. no leap years 31,28,31,30,31,30,31,31,30,31,30,31
!MIROC5      365 days, i.e. no leap years 31,28,31,30,31,30,31,31,30,31,30,31
!GFDL-ESM2M  365 days, i.e. no leap years 31,28,31,30,31,30,31,31,30,31,30,31
!ECHAM4      360 days = 12*30
!ECHO-G      360 days = 12*30
!Hadley      360 days = 12*30
!IPSL        360 days = 12*30
!HadGEM      360 days = 12*30


  

  interface operator(+)
     module procedure addtime !,addexp
  end interface
  interface operator(-)
     module procedure subtime,subtime2
  end interface
  interface operator(==)
     module procedure sametime
  end interface
  interface operator(/=)
     module procedure notsametime
  end interface
  interface operator(>)
     module procedure gt
  end interface
  interface operator(<)
     module procedure lt
  end interface
  interface operator(>=)
     module procedure ge
  end interface
  interface operator(<=)
     module procedure le
  end interface
  
  public operator(+),operator(-),operator(==),operator(>), &
       operator(<),operator(/=),operator(>=),operator(<=)
  public makedeltatime,maketime
  public secondOfDay,Nseconds,nDaysInMonth,Nseconds2,nDaysInYear

contains
  logical function gt(t1,t2)
    type(time),intent(in)::t1,t2

    gt = .false.
    if(t1%year>t2%year)then
       gt = .true.
       return
    elseif(t1%year==t2%year)then
       if(t1%month>t2%month)then
          gt = .true.
          return
       elseif(t1%month==t2%month)then
          if(t1%day>t2%day)then
             gt = .true.
             return
          elseif(t1%day==t2%day)then
             if(t1%hour>t2%hour)then
                gt = .true.
                return
             elseif(t1%hour==t2%hour)then  
                if(t1%min>t2%min)then
                   gt = .true.
                   return
                elseif(t1%min==t2%min)then
                   if(t1%sec>t2%sec)then
                      gt = .true.
                      return
                   endif
                endif
             endif
          endif
       endif
    endif
    return
  end function gt

  logical function le(t1,t2)
    type(time),intent(in)::t1,t2
    le = t1<t2.or.t1==t2
  end function le

  logical function ge(t1,t2)
    type(time),intent(in)::t1,t2
    ge = t1>t2.or.t1==t2
  end function ge

  logical function lt(t1,t2)
    type(time),intent(in)::t1,t2
    lt = gt(t2,t1)
    return
  end function lt

  integer function secondOfDay(t)
    implicit none
    type(time),intent(in)::t
    secondOfDay = (t%hour*60+t%min)*60+t%sec
  end function SecondOfDay

  integer function Nseconds2(dt,t)
    implicit none
    type(deltatime),intent(in)::dt
    type(time)::t
    integer::nd
    nd = nDaysInMonth(t%year,t%month)

    Nseconds2 = ((( (dt%month)*nd+dt%day)*24 +dt%hour)*60+dt%min)*60+dt%sec

  end function Nseconds2

  integer function Nseconds(dt)
    implicit none
    type(deltatime),intent(in)::dt

    Nseconds = ((dt%day*24 +dt%hour)*60+dt%min)*60+dt%sec
  end function Nseconds

  logical function sametime(t1,t2)
    implicit none  
    type(time),intent(in)::t1,t2
    sametime = t1%year==t2%year.and.t1%month==t2%month.and.t1%day==t2%day.and. &
         t1%hour==t2%hour.and.t1%min==t2%min.and.t1%sec==t2%sec
  end function sametime
  logical function notsametime(t1,t2)
    implicit none  
    type(time),intent(in)::t1,t2
    notsametime = .not.sametime(t1,t2)
  end function notsametime

  type(time) function subtime(t1,dt)
    implicit none  
    type(time),intent(in)::t1
    type(deltatime),intent(in)::dt

    subtime = maketime(t1%year-dt%year,t1%month-dt%month,t1%day-dt%day,&
         t1%hour-dt%hour,t1%min-dt%min,t1%sec-dt%sec)
  end function subtime

  type(deltatime) function dtLeft2endOfMonth(t)
    implicit none
    type(time),intent(in)::t
    dtLeft2endOfMonth%year=0
    dtLeft2endOfMonth%month=0
    dtLeft2endOfMonth%day=nDaysInMonth(t%year,t%month)-t%day
    dtLeft2endOfMonth%hour = 0
    dtLeft2endOfMonth%min = 0
    dtLeft2endOfMonth%sec = 0
    if(t%hour>0)then
       dtLeft2endOfMonth%hour=24-t%hour
    endif
    if(t%min>0)then
       dtLeft2endOfMonth%min=60-t%min
    endif
    if(t%sec>0)then
       dtLeft2endOfMonth%sec=60-t%sec
    endif
  end function dtLeft2endOfMonth

  type(deltatime) function dtLeft2endOfyear(t)
    implicit none
    type(time),intent(in)::t
    dtLeft2endOfyear%year=0
    dtLeft2endOfyear%month=12-t%month
    dtLeft2endOfyear%day=nDaysInMonth(t%year,t%month)-t%day
    dtLeft2endOfyear%hour = 0
    dtLeft2endOfyear%min = 0
    dtLeft2endOfyear%sec = 0
    if(t%hour>0)then
       dtLeft2endOfyear%hour=24-t%hour
    endif
    if(t%min>0)then
       dtLeft2endOfyear%min=60-t%min
    endif
    if(t%sec>0)then
       dtLeft2endOfyear%sec=60-t%sec
    endif
  end function dtLeft2endOfyear

  type(deltatime) function subtime2(t1,t2)
    implicit none  
    type(time),intent(in)::t1
    type(time),intent(in)::t2

    if(t1%month==t2%month.and.t1%year==t2%year)then
       subtime2 = makedeltatime(t1%year-t2%year, &
            t1%month-t2%month,&
            t1%day-t2%day,&
            t1%hour-t2%hour,&
            t1%min-t2%min, &
            t1%sec-t2%sec)
    elseif(t1%year==t2%year)then
       subtime2 = dtLeft2endOfMonth(t2)
       subtime2%day = subtime2%day + t1%day - 1 !days start with 1 NOT 0
       subtime2%hour = subtime2%hour + t1%hour
       subtime2%min = subtime2%min + t1%min
       subtime2%sec = subtime2%sec + t1%sec
    else
       subtime2  = dtLeft2endOfyear(t2)
       subtime2%month = subtime2%month + t1%month -1
       subtime2%day = subtime2%day + t1%day -1
       subtime2%hour = subtime2%hour + t1%hour
       subtime2%min = subtime2%min + t1%min
       subtime2%sec = subtime2%sec + t1%sec
    endif
  end function subtime2

  type(time) function maketime(year,month,day,hour,min,sec)
    implicit none
    integer,intent(in)::year,month,day,hour,min,sec
    integer::iyear,imonth,iday,ihour,imin,isec

    isec = mod(sec,60)
    imin = min + sec/60

    ihour = hour + imin/60
    imin = mod(imin,60)

    iday = day + ihour/24
    ihour = mod(ihour,24)

    iyear = year 
    imonth = month 



    do while(iday<1)
       imonth = imonth - 1
       do while(imonth<=0)
          imonth = 12-imonth
          iyear = iyear -1
       enddo
       iday = nDaysInMonth(iyear,imonth) - abs(iday)
    enddo

    do while(iday>nDaysInMonth(iyear,imonth)) !what if imonth==13?
       iday = iday - nDaysInMonth(iyear,imonth)
       imonth = imonth + 1 
       do while(imonth>12)
          imonth = imonth -12
          iyear = iyear + 1
       enddo
    enddo

    do while(imonth>12)
       iyear = iyear + 1
       imonth = imonth - 12
    enddo
    do while(imonth<1)
       iyear = iyear - 1
       imonth = imonth + 12
    enddo

    maketime%year  = iyear
    maketime%month = imonth
    maketime%day   = iday
    maketime%hour  = ihour 
    maketime%min   = imin
    maketime%sec   = isec
  end function maketime


  type(deltatime) function makedeltatime(dyear,dmonth,dday,dhour,dmin,dsec)
    implicit none
    integer,intent(in)::dyear,dmonth,dday,dhour,dmin,dsec
    integer::dPm
    logical::illegal
    !Restriction to positive deltatimes
    illegal = (dyear<0).or.&
         (dyear==0.and.dmonth<0).or. &
         (dyear==0.and.dmonth==0.and.dday<0).or. &
         (dyear==0.and.dmonth==0.and.dday==0.and.dhour<0).or. &
         (dyear==0.and.dmonth==0.and.dday==0.and.dhour==0.and.dmin<0).or.&
         (dyear==0.and.dmonth==0.and.dday==0.and.dhour==0.and.dmin==0.and.dsec<0)
    if(illegal)then
       if(mype==0)then
          print *,dyear,dmonth,dday,dhour,dmin,dsec
       endif
       call stop_program( 'illegal call to makedeltatime <0')
    endif
    makedeltatime%year  = dyear 
    makedeltatime%month = dmonth
    makedeltatime%day   = dday  
    makedeltatime%hour  = dhour     
    makedeltatime%min   = dmin      
    makedeltatime%sec   = dsec  

    if(makedeltatime%sec>=60)then
       makedeltatime%min = makedeltatime%min + (makedeltatime%sec-mod(makedeltatime%sec,60))/60
       makedeltatime%sec = mod(makedeltatime%sec,60)
    endif
    if(makedeltatime%min>=60)then
       makedeltatime%hour = makedeltatime%hour + (makedeltatime%min-mod(makedeltatime%min,60))/60
       makedeltatime%min = mod(makedeltatime%min,60)
    endif
    if(makedeltatime%hour>=24)then
       makedeltatime%day = makedeltatime%day + (makedeltatime%hour-mod(makedeltatime%hour,24))/24
       makedeltatime%hour = mod(makedeltatime%hour,24)
    endif
    if(dmonth>0)then
       dPm = nDaysInMonth(makedeltatime%year,makedeltatime%month)
       do while(makedeltatime%day>dPm)       
          makedeltatime%month = makedeltatime%month + 1
          makedeltatime%day = makedeltatime%day - dPm
          dPm = nDaysInMonth(makedeltatime%year,makedeltatime%month)
       end do
       do while(makedeltatime%month>12)
          makedeltatime%year = makedeltatime%year + 1
          makedeltatime%month = makedeltatime%month - 12
       enddo
    endif
  end function makedeltatime
    
  type(time) function addtime(t1,dt)
    implicit none  
    type(time),intent(in)::t1
    type(deltatime),intent(in)::dt

    addtime = maketime(t1%year +dt%year,t1%month+dt%month, t1%day+dt%day,&  
         t1%hour+dt%hour, t1%min+dt%min, t1%sec+dt%sec)
    if(isValid(addtime))then
       return
    else
       call stop_program('addtime not a valid time')
       return
    endif
!!$    do while(addtime%sec>=60)
!!$       addtime%sec = addtime%sec -60
!!$       addtime%min = addtime%min + 1
!!$    enddo
!!$    do while(addtime%min>=60)
!!$       addtime%min = addtime%min-60
!!$       addtime%hour = addtime%hour+1
!!$    enddo
!!$    do while(addtime%hour>=24)
!!$       addtime%hour = addtime%hour-24
!!$       addtime%day = addtime%day+1
!!$    enddo
!!$    if(addtime%day > nDaysInMonth(addtime%year,addtime%month))then
!!$       addtime%day = nDaysInMonth(addtime%year,addtime%month) - addtime%day
!!$       addtime%month = addtime%month + 1
!!$    endif
!!$    
!!$    if(addtime%month>12)then
!!$       addtime%year = addtime%year+1
!!$       addtime%month = 1
!!$    endif
  end function addtime

  type(time) function addexp(t1,dt)
    implicit none  
    type(time),intent(in)::t1
    integer,intent(in)::dt(6) !(y,m,d,h,m,s)
    type(deltatime)::tmp
    tmp%year = dt(1)  !y
    tmp%month = dt(2) !mnt
    tmp%day = dt(3) !d
    tmp%hour = dt(4) !h
    tmp%min = dt(5) !m
    tmp%sec = dt(6) !s
    addexp = addtime(t1,tmp)
  end function addexp

  integer function nDaysInYear(year)
    implicit none
    integer,intent(in)::year
    integer::m
    nDaysInYear=0
    do m=1,12
       nDaysInYear = nDaysInYear + nDaysInMonth(year,m)
    enddo
  end function nDaysInYear
    
  integer function nDaysInMonth(year,month)
    use gcm
    implicit none
    integer,intent(in)::year,month
    logical leap

    if(lhadley.or.lhadgem)then
       nDaysInMonth = 30
    elseif(lecham2)then
       nDaysInMonth = 30
       if(lecham5)then
          leap=.true.
          nDaysInMonth = ndayReal(year,month,leap)
       endif
    else
       if(lcanesm2.or.lccsm.or.lnoresm.or.lmiroc5.or.lipsl.or.lgfdl) then
          leap=.false.
       else
          leap=.true.
       endif
       nDaysInMonth = ndayReal(year,month,leap)
    endif
  end function nDaysInMonth


  integer function ndayReal(year,month,leap)
    implicit none
    integer,intent(in)::year,month
    logical,intent(in)::leap
    integer::iyear,imonth
    integer::daysPerMonth(12)
    daysPerMonth = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    iyear = year 
    imonth = month

    iyear = year + month/12
    imonth = max(1,mod(month,13))
    !determine if leapyear
! not for Canesm2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(leap) then
       if(mod(iyear,4)==0.and..not.mod(iyear,100)==0.or.mod(iyear,400)==0)then
          daysPerMonth(2) = 29
       endif
    endif
    ndayReal = daysPerMonth(imonth)
  end function ndayReal

  integer function ndaysPerYear(year)
    implicit none
    integer,intent(in)::year
    integer::m
    ndaysPerYear=0
    do m=1,12
       ndaysPerYear = ndaysPerYear + nDaysInMonth(year,m)
    enddo
  end function ndaysPerYear

  logical function isValid(t)
    type(time),intent(in)::t
    isValid = .true.
    if(t%month>12.or.t%month<1) isValid =.false.
    if(t%day>nDaysInMonth(t%year,t%month).or.t%day<1)isValid =.false.
    if(t%hour>24.or.t%hour<0) isValid =.false.
    if(t%min>60.or.t%min<0) isValid =.false.
    if(t%sec>60.or.t%sec<0) isValid =.false.
    return
  end function isValid

end module calendar
