module pblhgtmod
  use confys
  use decomp,only:realkind
  implicit none
  private
  
  public pblhgt
contains

  subroutine pblhgt(nhor,nlev,kstart,kstop,ustar,ps,ts,senf,latf, &
       pf,dph,gpot,t,q,cw,u,v, &
       pblh,     &           ! boundary layer height
       heatv,    &           ! surface virtual heat flux
       obklen)              ! obukhov length
    !     compute the boundary layer height, (pblh),
    !     surface virtual heat flux  (heatv),
    !     obukhov length             (obklen)
    !     
    !     
    !     the scheme is described in
    !     troen,i. and l. mahrt (1986), a simple model of the atmospheric
    !     boundary layer; sensitivity to surface evaporation.
    !     boundary-layer meteorol., 37, 129-148.
    !     the original bulk formulation (holtslag and boville, 1993)
    !     has been upgraded (niels woetmannn nielsen on 23/5 1996)
    !     with a new formula suggested by vogelzang and holtslag, 1996.
    !     in the new formula the virtual tmp. and vind vel. differences
    !     are taken between nlev-k and nlev. in addition the term
    !     b*ustar*ustar (b=100.) is added to the norm of the wind
    !     velocity shear.
 
    implicit none
    integer,intent(in):: nhor,nlev,kstart,kstop
    integer::nlevp,kpbl,ktopfl
         
    !      arrays received from the calling routine (phys)

    real(kind=realkind),intent(inout)::ustar(nhor)
    real(kind=realkind),intent(in)::ps(nhor),ts(nhor),senf(nhor),latf(nhor)
    !     
    real(kind=realkind),intent(in)::pf(nhor,nlev),gpot(nhor,nlev),dph(nhor,nlev+1), &
         t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev), &
         u(nhor,nlev),v(nhor,nlev)
         
    !     rrays produced in pblhgt for use in other parts of physics

    real(kind=realkind),intent(inout)::pblh(nhor),heatv(nhor),obklen(nhor)
         
    !     workspace declarations for the 'pblhgt'-routine.

    integer jwork(nhor)
    real(kind=realkind) wthm(nhor,nlev),wzm(nhor,nlev),wsenf(nhor),wcflx(nhor),wtseff(nhor),&
         wtav(nhor),wtlv(nhor),wrino(nhor)
     
    integer i,k,km,jck
    !     
    real(kind=realkind) tiny,vvk,fmt,wsc,tkv,zrino,ztherm,umun,vmvn
    real(kind=realkind) zunson,zunsoff
    real(kind=realkind) zslask
    real(kind=realkind) zlatw,zlati,zlatin,ztmp,zginv,zkappa,zp0,zrho
    real(kind=realkind) zbetam,zfak,zonet,zricr
    
    tiny = 1.e-9_realkind              ! bound wind**2 to guard divide by 0
    nlevp=nlev+1
    jck=0
    kpbl=int(real(nlev,realkind)*0.5_realkind)+1
    ktopfl=1
    !     
    zginv=1._realkind/gravit
    zkappa=rair/cpair
    zp0=1.e5_realkind
    zlatw=1._realkind/latvap
    zlati=1._realkind/(latvap +latice)
    zbetam=15.0_realkind
    zfak=8.5_realkind
    zonet=1._realkind/3._realkind
    zricr=0.25_realkind
    
    
    do 100 i=kstart,kstop
       
       jwork(i)=0
       wtseff(i)=ts(i)
       if( wtseff(i)<tmelt ) then
          zlatin=zlati
       else
          zlatin=zlatw
       endif
            
       zrho=gpot(i,nlev)/dph(i,nlev+1)
       wsenf(i)=-senf(i)*zrho/cpair
       wcflx(i)=-latf(i)*zrho*zlatin
            
       wzm(i,nlev)=gpot(i,nlev)*zginv
       wthm(i,nlev)= t(i,nlev)*(zp0/pf(i,nlev))**zkappa
       
100 enddo
        
    do k=nlev-1,ktopfl,-1
       
       do i=kstart,kstop
          
          !     compute relative height and the
          !     potential temperature.
          
          wzm(i,k)=gpot(i,k)*zginv
          wthm(i,k)= t(i,k)*(zp0/pf(i,k))**zkappa
       enddo
    enddo
    !     !  set up local arrays: compute bottom level virtual temperature
    !     !  and heat flux

    do 120 i=kstart,kstop
       if( t(i,nlev)>=tmelt ) then
          zlatin=zlatw
       else
          zlatin=zlati
       endif
       ztmp=1._realkind -cw(i,nlev)/(zlatin*cpair*t(i,nlev))
       wtav(i)=wthm(i,nlev)*ztmp*(1.0_realkind +0.61_realkind*q(i,nlev) -cw(i,nlev))
       heatv(i) = wsenf(i) + 0.61_realkind*wthm(i,nlev)*wcflx(i)
       ustar(i) = max(ustar(i),0.01_realkind)
       obklen(i) = -wtav(i)*ustar(i)**3.0_realkind/ &
            (gravit*carman*(heatv(i) + sign(1.e-10_realkind,heatv(i))))
120 enddo
         
    !  calculation of boundary layer height
         
    do 501 i=kstart,kstop
       jwork(i) = 1
       pblh(i) = 0._realkind
       wrino(i) = 0._realkind
       !     wtlv(i) is temporarily used to store 100.*ustar(i)*ustar(i).
       !     in the convective pbl (loops 507 and 512) wtlv(i) is used to
       !     store near surface temperature excess ( wtav(i)+ztherm).
       wtlv(i) = 100._realkind*ustar(i)*ustar(i)
       !     
501 enddo
         
    do 505 k=2,kpbl
       km = k - 1
       do 504 i=kstart,kstop
          umun = u(i,nlevp-k)-u(i,nlev)
          vmvn = v(i,nlevp-k)-v(i,nlev)
          vvk  = umun*umun + vmvn*vmvn + wtlv(i)
          vvk = max(vvk,tiny)
               
          if( t(i,nlevp-k)>=tmelt ) then
             zlatin=zlatw
          else
             zlatin=zlati
          endif
               
          ztmp=1._realkind -cw(i,nlevp-k)/(zlatin*cpair*t(i,nlevp-k))
          tkv=wthm(i,nlevp-k)*ztmp*(1._realkind+0.61_realkind*q(i,nlevp-k)-cw(i,nlevp-k))
          zrino = gravit*(tkv-wtav(i))*(wzm(i,nlevp-k)-wzm(i,nlev))/ &
               (wtav(i)*vvk)
          zrino=real(jwork(i),realkind)*zrino
               
          zslask=1.0_realkind
          if( zrino>=zricr ) then
             zslask=(zrino-wrino(i))/(wzm(i,nlevp-k)-wzm(i,nlevp-km))
             pblh(i)=wzm(i,nlevp-km) + (zricr-wrino(i))/zslask
             jwork(i)=0
          endif
          !     if pblh is found already avoid dividing by zero next level
          !     by a faked value of wrino(i) (jwork=0 already)
          
          zrino=real(jwork(i),realkind)*zrino
          wrino(i)=zrino+real(1-jwork(i),realkind)*0.1_realkind
               
504    enddo
505 enddo
         
    do 507 i=kstart,kstop
       pblh(i) = real(1-jwork(i),realkind)*pblh(i)+wzm(i,nlevp-kpbl)*real(jwork(i),realkind)
            
       if( heatv(i)>0.0_realkind ) then
          zunson=1.0_realkind
          zunsoff=0.0_realkind
       else
          zunson=0.0_realkind
          zunsoff=1.0_realkind
       endif
            
       !      unstable case; compute first estimate of velocity scale,
       !      and thermal temperature excess
            
       fmt = (zunson*(1._realkind -zbetam*wzm(i,nlev)/obklen(i)) + zunsoff)**zonet
       wsc = ustar(i)*fmt
       ztherm = heatv(i)*zfak/wsc
            
       !     improve phbl estimate under convective conditions using
       !     convective temperature excess (ztherm)
            
       wrino(i) = zunsoff*0.1_realkind
       wtlv(i) = zunson*( wtav(i) + ztherm)
       jwork(i) = int(zunson) 
            
507 enddo
    do 512 k=2,kpbl
            
       km = k - 1
            
       do 510 i=kstart,kstop
               
          umun = u(i,nlevp-k) - u(i,nlev)
          vmvn = v(i,nlevp-k) - v(i,nlev)
          vvk  = umun*umun + vmvn*vmvn + 100._realkind*ustar(i)*ustar(i)
               
          if( t(i,nlevp-k)>=tmelt ) then
             zlatin=zlatw
          else
             zlatin=zlati
          endif
               
          ztmp=1._realkind -cw(i,nlevp-k)/(zlatin*cpair*t(i,nlevp-k))
          tkv=wthm(i,nlevp-k)*ztmp*(1._realkind+0.61_realkind*q(i,nlevp-k)-cw(i,nlevp-k))
          vvk=max(vvk,tiny)
          zrino = gravit*(tkv-wtlv(i))*(wzm(i,nlevp-k)-wzm(i,nlev))/(vvk*wtav(i))
          zrino=zrino*real(jwork(i),realkind)
          !     
          zslask=1._realkind
          if( zrino>=zricr ) then
             zslask=(zrino-wrino(i))/(wzm(i,nlevp-k)-wzm(i,nlevp-km))
             pblh(i)=wzm(i,nlevp-km) + (zricr-wrino(i))/zslask
             jwork(i)=0
          endif
          !     if pblh is found already avoid dividing by zero next level
          !     by a faked value of wrino(i) (jwork=0 already)
          zrino=zrino*real(jwork(i),realkind)
          wrino(i)=zrino + real(1-jwork(i),realkind)*0.1_realkind
               
510    enddo
512 enddo

    do 514 i=kstart,kstop
       pblh(i) = real(1-jwork(i),realkind)*pblh(i)+real(jwork(i),realkind)*wzm(i,nlevp-kpbl)
514 enddo
    return
  end subroutine pblhgt
end module pblhgtmod
