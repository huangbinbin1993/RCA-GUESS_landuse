module interpol
  use decomp  
  use util
  
  implicit none
  private

  public vint,turn_wind

contains

  subroutine vint (nlev, nhor, fin, vcin, fout, vcout)  
    implicit none  
    integer,intent(in)::nlev,nhor  
    real(kind=realkind)::vcin(nhor,nlev),fin(nhor,nlev)
    real(kind=realkind),intent(inout)::vcout(nhor),fout(nhor)

    integer::j,k
    do j=1,nhor
       if(vcout(j)<=vcin(j,1))then
          fout(j)=fin(j,1)  
       elseif(vcout(j)>vcin(j,nlev))then  
          fout(j)=fin(j,nlev)  
       endif
    enddo
    do k=1,nlev-1  
       do j=1,nhor  
          if(vcout(j)>vcin(j,k).and.vcout(j)<=vcin(j,k + 1))then
             fout(j)=fin(j,k)+log(vcout(j)/vcin(j,k))/ &
                  log(vcin(j,k+1)/vcin(j,k))*(fin(j,k+1)-fin(j,k))
          endif
       enddo

    enddo
    return  
  end subroutine vint




  subroutine turn_wind (nlon, nlat, west, south, aplon, aplat, &
       dlamda, dtheta, u10, v10, u10_reg, v10_reg)
    implicit none  
    !     Anders Ullerstig, 991119,
    integer :: nlev, nlon, nlat  
    real(kind=realkind) :: west, south, aplon, aplat, dlamda, dtheta  
    real(kind=realkind) :: u10 (nlon, nlat), v10 (nlon, nlat)  


    real(kind=realkind) :: u10_reg (nlon, nlat), v10_reg (nlon, nlat)  
    integer :: i, j  
    real(kind=realkind) :: reg_lon(1,1), reg_lat(1,1), rot_lon(1,1), rot_lat(1,1)

    real(kind=realkind) :: pa (1, 1), pb (1, 1), pc (1, 1), pd (1, 1)  
    do j = 1, nlat  
       rot_lat (1, 1) = south + real (jdatastart + j - 2,realkind) * dtheta  
       do i = 1, nlon  
          rot_lon (1, 1) = west + real (idatastart + i - 2,realkind) * dlamda  
          call regrot(reg_lon, reg_lat, rot_lon, rot_lat, 1, 1,   aplon, aplat, - 1)
          call turnwirot2reg(pa,pb,pc,pd,reg_lon,reg_lat,rot_lon,rot_lat,1,1,&
               aplon, aplat, -2)
          call turnwi1(u10(i,j),v10(i,j),u10_reg(i,j),v10_reg(i, j),pa,pb,pc,pd,1,1,1)
       enddo
    enddo
    return  
  end subroutine turn_wind

!!$  subroutine interpolateGCM2LAM(gcm,LAM,fgcm,gLAM)
!!$    use gcm
!!$    use rcaDomainMod
!!$    implicit none
!!$    integer::i
!!$    
!!$  end subroutine interpolateGCM2LAM


end module interpol
