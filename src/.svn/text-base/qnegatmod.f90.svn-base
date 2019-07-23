module qnegatmod
  use decomp,only:realkind
  implicit none
  private

  public aqnegat
contains
  subroutine aqnegat (nhor,nlev,kstart,kstop,  &
       dtime,                                  &
       dqdtin,q,dpf,                           &
       dqdt)
    implicit none
    integer nhor,nlev,kstart,kstop
    real(kind=realkind) dtime
    real(kind=realkind) dqdtin(nhor,nlev),dqdt(nhor,nlev),q(nhor,nlev),dpf(nhor,nlev)
    !
    ! local workspace

    real(kind=realkind) wqneg(nhor)

    call qnegat (nhor,nlev,kstart,kstop, &
         dtime,                          &
         dqdtin,q,dpf,                   &
         dqdt,                           &
         wqneg)
    !
    return
  end subroutine aqnegat

  subroutine qnegat (nhor,nlev,kstart,kstop, &
       dtime,                                 &
       dqdtin,q,dpf,                          &
       dqdt,                                  &
       wqneg)
    !     ----------------------------------------------------------

    !     l   subroutine qnegat:  nwn and bhs , dmi, dec 1987
    !     l   restructured        sg , smhi,         dec 1990
    !     ----------------------------------------------------------
    !     
    implicit none
    integer nhor,nlev,kstart,kstop
    integer jk,jl,jkm1
    real(kind=realkind) dtime
    real(kind=realkind) zmin,zqp,zqpin
    real(kind=realkind) dqdtin(nhor,nlev),dqdt(nhor,nlev),q(nhor,nlev),dpf(nhor,nlev)
    real(kind=realkind) wqneg(nhor)
    !     l    do vertical mixing of moisture for gridpoints
    !     l    with negative values of specific humidity.
    !     see 'hirlam documentation manual'  1988.
    !************************integer input ********************************
    !     l
    !     l      nhor            input           dimension length in the horizont
    !     l
    !     l      nlev            input           number of vertical levels
    !     l
    !     l      kstart          input           starting index for horizontal lo
    !     l
    !     l      kstop           input           ending index for hor. loops
    !     l
    !     l      dtime           input           time step for preliminary foreca
    !     l
    !************************full fields input ****************************
    !     l
    !     l      dqdtin(nhor,nlev) input        "preliminary" specific humidity
    !     l                                      tendencies
    !     l
    !     l      q(nhor,nlev)    input           specific humidities
    !     l
    !     l      dpf(nhor,nlev)  input           pressure differences at full lev
    !     l
    !************************full field output ****************************
    !     l
    !     l      dqdt(nhor)      output          tendency due to non negative hum
    !     l
    !***********************************************************************

    zmin=1.e-6_realkind

    !     do vertical mixing of moisture
    !     highest level

    do 210 jl=kstart,kstop
       zqpin = q(jl,1) + dtime*dqdtin(jl,1)
       wqneg(jl) = min(zqpin-zmin, 0._realkind)
       zqp = max(zqpin,zmin)
       dqdt(jl,1) = (zqp-zqpin)/dtime
210 enddo


    !     levels below the highest level

    do 250 jk=2,nlev
       jkm1 = jk - 1
       do 230 jl=kstart,kstop
          zqpin = q(jl,jk) + dtime*dqdtin(jl,jk)
          zqp = zqpin + wqneg(jl)*dpf(jl,jkm1)/dpf(jl,jk)
          wqneg(jl) = min(zqp-zmin, 0._realkind )
          zqp = max(zqp,zmin)
          dqdt(jl,jk)=(zqp-zqpin)/dtime
230    enddo
250 enddo

    return
  end subroutine qnegat
end module qnegatmod
