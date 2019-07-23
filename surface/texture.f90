  !--------------------------------------------------------------------------
  !                            INPOLY
  !   Function to tell if a point is inside a polygon or not.
  !--------------------------------------------------------------------------
  !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
  !
  !   Please feel free to use this source code for any purpose, commercial
  !   or otherwise, as long as you don't restrict anyone else's use of
  !   this source code.  Please give credit where credit is due.
  !
  !   Point-in-polygon algorithm, created especially for World-Wide Web
  !   servers to process image maps with mouse-clickable regions.
  !
  !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
  !
  !                                       6/19/95 - Bob Stein & Craig Yap
  !                                       stein@gcomm.com
  !                                       craig@cse.fau.edu
  !--------------------------------------------------------------------------
  !   Modified by:
  !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
  !                            2.  made it work with real coordinates
  !                            3.  now resolves the case where point falls
  !                                on polygon border.
  !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
  !   Aris Gerakis, july 1999: Now all borderline cases are valid
  !--------------------------------------------------------------------------
  !   Glossary:
  !   function inpoly: true=inside, false=outside (is target point inside
  !                    a 2D polygon?)
  !   poly(*,2):  polygon points, [0]=x, [1]=y
  !   npoints: number of points in polygon
  !   xt: x (horizontal) of target point
  !   yt: y (vertical) of target point
  !--------------------------------------------------------------------------
logical function inpoly(poly, npoints, sand, clay)
  implicit none
  integer,intent(in):: npoints
  real(kind=8),intent(in):: sand, clay
  real(kind=8),intent(in):: poly(7, 2)

  real(kind=8)::xnew, ynew, xold, yold, x1, y1, x2, y2

  integer:: i
  logical:: inside, on_border

  inside = .false.
  on_border = .false.

  if(npoints>7)then
     stop 'npoints too large in inpoly'
  endif
  if(npoints<3)then
     inpoly = .false.
     return
  endif

  xold = poly(npoints,1)
  yold = poly(npoints,2)

  do i = 1 , npoints
     xnew = poly(i,1)
     ynew = poly(i,2)

     if (xnew > xold)  then
        x1 = xold
        x2 = xnew
        y1 = yold
        y2 = ynew
     else
        x1 = xnew
        x2 = xold
        y1 = ynew
        y2 = yold
     end if

     ! The outer IF is the 'straddle' test and the 'vertical border' test.
     ! The inner IF is the 'non-vertical border' test and the 'north' test.  

     ! The first statement checks whether a north pointing vector crosses  
     ! (stradles) the straight segment.  There are two possibilities, depe-
     ! nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
     ! must be "open" at left, which is necessary to keep correct count when 
     ! vector 'licks' a vertix of a polygon.  

     if ((xnew < sand .and. sand <= xold) .or. &
          (.not. xnew < sand .and.  .not. sand <= xold)) then
        ! The test point lies on a non-vertical border:
        if (abs((clay-y1)*(x2-x1) - (y2-y1)*(sand-x1))<1.e-14_8) then
           on_border = .true. 
           ! Check if segment is north of test point.  If yes, reverse the 
           ! value of INSIDE.  The +0.001 was necessary to avoid errors due   
           ! arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
        elseif ((clay-y1)*(x2-x1) < ((y2-y1)*(sand-x1) + DBLE(0.001))) then
           inside = .not.inside ! cross a segment
        endif
        ! This is the rare case when test point falls on vertical border or  
        ! left edge of non-vertical border. The left x-coordinate must be  
        ! common.  The slope requirement must be met, but also point must be
        ! between the lower and upper y-coordinate of border segment.  There 
        ! are two possibilities,  depending on whether ynew < yold or ynew > 
        ! yold:
     elseif ((abs(xnew - sand)<1.e-14_8 .or. &
          abs(xold - sand)<1.e-14_8) .and.&
          abs((clay-y1)*(x2-x1) - (y2-y1)*(sand-x1))<1.e-14_8 .and.&
          ((ynew <= clay .and. clay <= yold) .or.&
          (.not. ynew < clay .and. .not. clay < yold))) then
        on_border = .true. 
     endif

     xold = xnew
     yold = ynew

  enddo

  ! If test point is not on a border, the function result is the last state 
  ! of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
  ! inside the polygon if it falls on any of its borders:

  if (.not. on_border) then
     inpoly = inside
  else
     inpoly = .true.
  endif

  return

end function inpoly


integer function what_texture(sand, clay)
  !* +-----------------------------------------------------------------------
  !* | WHAT TEXTURE?
  !* | Function to classify a soil in the triangle based on sand and clay %
  !* +-----------------------------------------------------------------------
  !* | Created by: aris gerakis, apr. 98
  !* | Modified by: aris gerakis, june 99.  Now check all polygons instead of
  !* | stopping when a right solution is found.  This to cover all borderline 
  !* | cases.
  !* +-----------------------------------------------------------------------

  implicit none

  ! Declare arguments:

  real(kind=8),intent(in):: clay, sand

  ! Declare local variables:

  integer   itex
  interface
     logical function inpoly(poly, npoints, xt, yt)
       integer,intent(in):: npoints
       real(kind=8),intent(in):: xt, yt
       real(kind=8),intent(in):: poly(7, 2)
     end function inpoly
  end interface
  real(kind=8)::silty_loam(7,2), sandy(7,2), silty_clay_loam(7,2)
  real(kind=8)::loam(7,2), clay_loam(7,2), sandy_loam(7,2)
  real(kind=8)::silty_clay(7,2), sandy_clay_loam(7,2), loamy_sand(7,2)
  real(kind=8)::clayey(7,2), silt(7,2), sandy_clay(7,2)

  ! Initalize polygon coordinates:
  silty_loam(:,1) = (/0d0, 0d0, 0.23d0, 0.50d0, 0.20d0, 0.08d0, 0d0/)
  silty_loam(:,2) = (/0.12d0,0.27d0, 0.27d0, 0d0, 0d0, 0.12d0, 0d0/)
  sandy(:,1) = (/0.85d0, 0.90d0, 1.00d0, 0d0, 0d0, 0d0, 0d0/)
  sandy(:,2) = (/0d0,  0.10d0,   0d0, 0d0, 0d0, 0d0, 0d0/)
  silty_clay_loam(:,1) = (/ 0d0,  0d0, 0.20d0, 0.20d0, 0d0, 0d0, 0d0/)
  silty_clay_loam(:,2) = (/0.27d0, 0.40d0, 0.40d0, 0.27d0, 0d0, 0d0, 0d0/)
  loam(:,1) = (/0.43d0, 0.23d0, 0.45d0, 0.52d0, 0.52d0, 0d0, 0d0/)
  loam(:,2) = (/ 0.07d0, 0.27d0, 0.27d0, 0.20d0,  0.07d0, 0d0, 0d0/)
  clay_loam(:,1) = (/0.20d0, 0.20d0, 0.45d0, 0.45d0, 0d0, 0d0, 0d0/)
  clay_loam(:,2) = (/0.27d0, 0.40d0, 0.40d0, 0.27d0, 0d0, 0d0, 0d0/)
  sandy_loam(:,1) = (/0.50d0, 0.43d0, 0.52d0, 0.52d0, 0.80d0, 0.85d0, 0.70d0/)
  sandy_loam(:,2) = (/0d0,   0.07d0,  0.07d0, 0.20d0, 0.20d0, 0.15d0, 0d0/)
  silty_clay(:,1) = (/ 0d0,  0d0, 0.20d0, 0d0, 0d0, 0d0, 0d0/)
  silty_clay(:,2) = (/0.40d0, 0.60d0, 0.40d0, 0d0, 0d0, 0d0, 0d0/)
  sandy_clay_loam(:,1) = (/0.52d0, 0.45d0, 0.45d0, 0.65d0, 0.80d0, 0d0, 0d0/)
  sandy_clay_loam(:,2) = (/0.20d0, 0.27d0, 0.35d0, 0.35d0, 0.20d0, 0d0, 0d0/)
  loamy_sand(:,1) = (/0.70d0, 0.85d0, 0.90d0, 0.85d0, 0d0, 0d0, 0d0/)
  loamy_sand(:,2) = (/ 0d0, 0.15d0, 0.10d0,  0d0, 0d0, 0d0, 0d0/)
  clayey(:,1) = (/0.20d0,  0d0,   0d0, 0.45d0, 0.45d0, 0d0, 0d0/)
  clayey(:,2) = (/0.40d0, 0.60d0, 1.00d0, 0.55d0, 0.40d0, 0d0, 0d0/)
  silt(:,1) = (/0d0,  0d0,  0.08d0, 0.20d0, 0d0, 0d0, 0d0/)
  silt(:,2) = (/0d0, 0.12d0, 0.12d0,  0d0, 0d0, 0d0, 0d0/)
  sandy_clay(:,1) = (/0.45d0, 0.45d0, 0.65d0, 0d0, 0d0, 0d0, 0d0/)
  sandy_clay(:,2) = (/0.35d0, 0.55d0, 0.35d0, 0d0, 0d0, 0d0, 0d0/)
  

  ! Find polygon(s) where the point is.  
  itex = 0

  if (sand > 0d0 .and. clay > 0d0) then

     if (inpoly(silty_loam, 6, sand, clay)) then
        itex=1
        !      texture = trim(texture)//'/silt loam'
        !     endif
     elseif(inpoly(sandy, 3, sand, clay)) then
        itex=2
        !      texture = trim(texture)//'/sand'
        !     endif
     elseif(inpoly(silty_clay_loam, 4, sand, clay)) then
        itex=3
        !      texture = trim(texture)//'/silty clay loam'
        !     endif
     elseif(inpoly(loam, 5, sand, clay)) then
        itex=4
        !      texture = trim(texture)//'/loam'
        !     endif
     elseif(inpoly(clay_loam, 4, sand, clay)) then
        itex=5
        !      texture = trim(texture)//'/clay loam'
        !     endif
     elseif(inpoly(sandy_loam, 7, sand, clay)) then
        itex=6
        !      texture = trim(texture)//'/sandy loam'
        !     endif
     elseif(inpoly(silty_clay, 3, sand, clay)) then
        itex=7
        !      texture = trim(texture)//'/silty clay'
        !     endif
     elseif(inpoly(sandy_clay_loam, 5, sand, clay)) then
        itex=8
        !      texture = trim(texture)//'/sandy clay loam'
        !     endif
     elseif(inpoly(loamy_sand, 4, sand, clay)) then
        itex=9
        !      texture = trim(texture)//'/loamy sand'
        !     endif
     elseif(inpoly(clayey, 5, sand, clay)) then
        itex=10
        !      texture = trim(texture)//'/clay'
        !     endif
     elseif(inpoly(silt, 4, sand, clay))then
        itex=11
        !      texture = trim(texture)//'/silt'
        !     endif
     elseif(inpoly(sandy_clay, 3, sand, clay)) then
        itex=12
        !      texture = trim(texture)//'/sandy clay'
     endif
  endif

  what_texture = itex

  return
end function what_texture






