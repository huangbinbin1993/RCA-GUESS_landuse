subroutine rot2reg(xreg,yreg,xrot,yrot,kx,ky,xcen,ycen)

  implicit none
  integer,intent(in)::kx,ky
  real(kind=8),intent(in):: xcen,ycen
  real(kind=8),intent(in):: xrot(kx,ky),yrot(kx,ky)
  real(kind=8),intent(out):: xreg(kx,ky),yreg(kx,ky)

  real(kind=8):: pi,zrad,zsycen,zcycen
  real(kind=8):: zxmxc,zsxmxc,zcxmxc,zsyreg
  real(kind=8):: zcyreg,zsyrot,zcyrot,zcxrot,zsxrot,zradi
  integer:: jy,jx
  
  pi = 4d0*atan(1d0)
  zrad = pi/180d0
  zradi = 1d0/zrad
  zsycen = sin(zrad*(ycen+90d0))
  zcycen = cos(zrad*(ycen+90d0))
  
  do jy = 1,ky
     do jx = 1,kx
        zsxrot = sin(zrad*xrot(jx,jy))
        zcxrot = cos(zrad*xrot(jx,jy))
        zsyrot = sin(zrad*yrot(jx,jy))
        zcyrot = cos(zrad*yrot(jx,jy))
        zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
        zsyreg = max(zsyreg,-1d0)
        zsyreg = min(zsyreg,+1d0)

        yreg(jx,jy) = asin(zsyreg)*zradi

        zcyreg = cos(yreg(jx,jy)*zrad)
        zcxmxc = (zcycen*zcyrot*zcxrot - zsycen*zsyrot)/zcyreg
        zcxmxc = max(zcxmxc,-1d0)
        zcxmxc = min(zcxmxc,+1d0)
        zsxmxc = zcyrot*zsxrot/zcyreg
        zxmxc  = acos(zcxmxc)*zradi
        if (zsxmxc<0d0) zxmxc = -zxmxc
        xreg(jx,jy) = zxmxc + xcen
        call CC(xreg(jx,jy),yreg(jx,jy))
     enddo
  enddo
  
  
  return
end subroutine rot2reg


subroutine CC(lon,lat)
  !Makes lon/lat into correct coordinates
  implicit none
  real(kind=8),intent(inout)::lon,lat
  integer::iter
  iter = 0
  if(abs(lat)>90d0)then
     if(lat>90d0)then
        lat = 90d0 - (lat-90d0)
     else
        lat = -90d0 - (lat+90d0)
     endif
     lon = lon + 180d0
  endif
  
  do while( abs(lon)>180d0 .and. iter<10) 
     iter = iter+1
     lon = lon - sign(1d0,lon)*360d0
  enddo
  if(iter==10)stop'out of iterations in CC'

end subroutine CC
