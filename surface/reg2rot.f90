subroutine reg2rot(xreg,yreg,xrot,yrot,kx,ky,xcen,ycen)

  implicit none
  integer,intent(in):: kx,ky
  real(kind=8),intent(in):: xcen,ycen
  real(kind=8),intent(in):: xreg(kx,ky),yreg(kx,ky)
  real(kind=8),intent(out):: xrot(kx,ky),yrot(kx,ky)
  real(kind=8):: pi,zrad,zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg
  real(kind=8):: zcyreg, zsyrot,zcyrot,zcxrot,zsxrot,zradi
  integer:: jy,jx
  pi = 4d0*atan(1d0)
  zrad = pi/180d0
  zradi = 1d0/zrad
  zsycen = sin(zrad*(ycen+90d0))
  zcycen = cos(zrad*(ycen+90d0))

  do jy = 1,ky
     do jx = 1,kx
        zxmxc  = zrad*(xreg(jx,jy) - xcen)
        zsxmxc = sin(zxmxc)
        zcxmxc = cos(zxmxc)
        zsyreg = sin(zrad*yreg(jx,jy))
        zcyreg = cos(zrad*yreg(jx,jy))
        zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
        zsyrot = max(zsyrot,-1d0)
        zsyrot = min(zsyrot,+1d0)

        yrot(jx,jy) = asin(zsyrot)*zradi !output

        zcyrot = cos(yrot(jx,jy)*zrad)
        zcxrot = (zcycen*zcyreg*zcxmxc +  zsycen*zsyreg)/zcyrot
        zcxrot = max(zcxrot,-1d0)
        zcxrot = min(zcxrot,+1d0)
        zsxrot = zcyreg*zsxmxc/zcyrot

        xrot(jx,jy) = acos(zcxrot)*zradi !output

        if (zsxrot<0d0) xrot(jx,jy) = -xrot(jx,jy)

        call CC(xrot(jx,jy),yrot(jx,jy))
     enddo
  enddo

  return
end subroutine reg2rot
