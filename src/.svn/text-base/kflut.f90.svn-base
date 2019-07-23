module kflut
  use confys, only:epsilo
  use decomp,only:realkind
  implicit none
  private
  integer,public,parameter::kfnt=250,kfnp=220
  real(kind=realkind),public,save::ttab(kfnt,kfnp)
  real(kind=realkind),public,save::qstab(kfnt,kfnp)
  real(kind=realkind),public,save::the0k(kfnp)
  real(kind=realkind),public,save::alu(200)
  real(kind=realkind),public,save::rdpr
  real(kind=realkind),public,save::rdthk
  real(kind=realkind),public,save::ptop
  public lutab

contains
  subroutine lutab

    !     this subroutine is a lookup table.
    !     given a series of series of saturation equivalent potential 
    !     temperatures, the temperature is calculated.

    implicit none

    real(kind=realkind) dth,tmin,toler,pbot,aliq,bliq,cliq,dliq,dpr,temp,p
    integer kp,it,itcnt,i

    real(kind=realkind) es,qs,pi,thes,thgues,f0,t1,t0,thtgs,f1,tgues,dt,astrt,ainc,a1
    !     gjesat
    !     gjesat-----------------------------------------------------------
    !     gjesat
    !     gjesat  use the hirlam reference esat calcs in kfcumulus to get a more
    !     gjesat  accurate value at low temps that takes ice into account


    !     equivalent potential temperature increment
    data dth/1._realkind/
    !     minimum starting temp 
    data tmin/150._realkind/
    !     top pressure (pascals)
    !     data ptop/100.0/
    !     bottom pressure (pascals)
    !     data pbot/110000.0/
    !     tolerance for accuracy of temperature 
    data toler/0.001_realkind/
    !     top pressure (pascals)
    ptop=100.0_realkind
    !     bottom pressure (pascals)
    pbot=110000.0_realkind

    !     
    !...  define constants for calculation of saturation vapor pressure
    !...  according to buck (j. appl. meteo., december, 1981)...
    !     
    aliq = 613.3_realkind
    bliq = 17.502_realkind
    cliq = 4780.8_realkind
    dliq = 32.19_realkind
    !     
    !     compute parameters
    !     
    !     1./(sat. equiv. theta increment)
    rdthk=1._realkind/dth
    !     pressure increment
    dpr=(pbot-ptop)/real(kfnp-1,realkind)
    !     1./(pressure increment)
    rdpr=1._realkind/dpr
    !     compute the spread of thes
    !     thespd=dth*(kfnt-1)
    !     
    !     calculate the starting sat. equiv. theta
    !     
    temp=tmin 
    p=ptop-dpr
    do 100 kp=1,kfnp
       p=p+dpr
       es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
       !     gjesat      es=esat(temp)
       qs=epsilo*es/(p-es)
       pi=(1.e5_realkind/p)**(0.2854_realkind*(1._realkind-0.28_realkind*qs))
       the0k(kp)=temp*pi*exp((3374.6525_realkind/temp-2.5403_realkind)*qs*(1._realkind+0.81_realkind*qs))
100 enddo

    !     
    !     compute temperatures for each sat. equiv. potential temp.
    !     
    p=ptop-dpr
    do 110 kp=1,kfnp
       thes=the0k(kp)-dth
       p=p+dpr
       do 120 it=1,kfnt
          !     define sat. equiv. pot. temp.
          thes=thes+dth
          !     iterate to find temperature
          !     find initial guess
          if(it==1) then
             tgues=tmin
          else
             tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          !     gjesat       es=esat(tgues)
          qs=epsilo*es/(p-es)
          pi=(1.e5_realkind/p)**(0.2854_realkind*(1._realkind-0.28_realkind*qs))
          thgues=tgues*pi*exp((3374.6525_realkind/tgues-2.5403_realkind)*qs*(1._realkind+0.81_realkind*qs))
          f0=thgues-thes
          t1=tgues-0.5_realkind*f0
          t0=tgues
          itcnt=0
          !     iteration loop
130       continue
          es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
          !     gjesat       es=esat(t1)
          qs=epsilo*es/(p-es)
          pi=(1.e5_realkind/p)**(0.2854_realkind*(1._realkind-0.28_realkind*qs))
          thtgs=t1*pi*exp((3374.6525_realkind/t1-2.5403_realkind)*qs*(1._realkind+0.81_realkind*qs))
          f1=thtgs-thes
          if(abs(f1)<toler)goto 140
          itcnt=itcnt+1
          if(itcnt>10) then
             print *,' itcnt > 10',' it=',it,' p=',p,' t1=',t1,' thes=',thes
             goto 140
          endif
          dt=f1*(t1-t0)/(f1-f0)
          t0=t1
          f0=f1
          t1=t1-dt
          goto 130
140       continue
          ttab(it,kp)=t1 
          qstab(it,kp)=qs
120    enddo
110 enddo
    !     
    !     lookup table for tlog(emix/aliq)
    !     
    !     set up intial values for lookup tables
    !     
    astrt=1.e-3_realkind
    ainc=0.075_realkind
    !     
    a1=astrt-ainc
    do 200 i=1,200
       a1=a1+ainc
       alu(i)=log(a1)
200 enddo
    !     
    return
  end subroutine lutab


end module kflut
