module kuomod
  use confys
  use escom
  use decomp,only:realkind
  implicit none
  private

  public akuo,inicond
contains
  subroutine inicond(nlev, hybf,c1er)

    !  define constants (first time only) used in subroutine cond and kuo

    implicit none
    integer,intent(in):: nlev
    real(kind=realkind),intent(in):: hybf(nlev)
    real(kind=realkind),intent(out)::c1er(nlev)
    integer::jk
    do 250 jk=1,nlev
       c1er(jk)=1.93e-6_realkind*282._realkind*sqrt( 1.e3_realkind/7.35_realkind*sqrt(hybf(jk)))
250 enddo

    return
  end subroutine inicond

  subroutine akuo(nhor,nlev,kstart,kstop,      &
       dtime,                         &
       nlstat,                        &
       t,q,virt,dtdtin,dqdtin,pf,dpf, &
       ps,dpsdin,ts,dtsdin,           &
       draindt,dsnowdt,               &
       prerat,                        &
       c1er,stcov,stdcov,             &
       dtdt,dqdt)

    implicit none
    integer nhor,nlev,kstart,kstop
    real(kind=realkind) dtime
    logical nlstat
    real(kind=realkind)     t(nhor,nlev),      q(nhor,nlev), virt(nhor,nlev), &
         dtdtin(nhor,nlev), dqdtin(nhor,nlev),   pf(nhor,nlev),  &
         dpf(nhor,nlev),   dtdt(nhor,nlev), dqdt(nhor,nlev),  &
         ps(nhor), dpsdin(nhor),     ts(nhor), dtsdin(nhor), &
         prerat(nhor), draindt(nhor),dsnowdt(nhor),              &
         c1er(nlev),stcov(nlev),stdcov(nlev)

    !     work space
    integer label(nhor,nlev)
    real(kind=realkind) wqp(nhor,nlev),wtp(nhor,nlev),xqp,xtp
    real(kind=realkind) wcldcp(nhor,nlev),wqsat1(nhor,nlev),xtvirt(nhor,nlev+1)
    real(kind=realkind) wcpqm1(nhor,nlev),wbfrac(nhor,nlev),wparcl(nhor,nlev)
    real(kind=realkind) wqac(nhor,nlev),wqc(nhor,nlev),wtc(nhor,nlev)
    real(kind=realkind) wgpot(nhor,nlev),wdqtoq(nhor,nlev),wdqtot(nhor,nlev)
    real(kind=realkind) wpsp(nhor),wtsp(nhor)
    real(kind=realkind) wcucov(nhor),wcpqn(nhor),wsdsig(nhor)
    real(kind=realkind) wcuprt(nhor),wqacpb(nhor)

    call kuo(nhor,nlev,kstart,kstop,                   &
         dtime,                               &
         nlstat,                              &
         t,q,virt,dtdtin,dqdtin,pf,dpf,       &
         ps,dpsdin,ts,dtsdin,                 &
         draindt,dsnowdt,                     &
         prerat,                              &
         c1er,stcov,stdcov,                   &
         dtdt,dqdt,                           &
         label,                               &
         xtvirt,                              &
         wtp,wqp,wqac,wparcl,wtc,             &
         wqc,wcldcp,wbfrac,wqsat1,wcpqm1,   &
         wgpot,wdqtoq,wdqtot,               &
         wpsp,wtsp,wcucov,wcpqn,wsdsig,       &
         wcuprt,wqacpb)
    !
    return
  end subroutine akuo

  subroutine kuo(nhor,nlev,kstart,kstop, & !dim,loop contr.
       dtime,              & !input coeff.
       nlstat,             & !input switch
       t,q,virt,dtdtin,dqdtin,pf,dpf, & !2-d input
       ps,dpsdin,ts,dtsdin, & !1-d input
       draindt,dsnowdt,    & !1-d output
       prerat,             & !1-d input/outpu
       c1er,stcov,stdcov,  & !1-d vert. input
       dtdt,dqdt,          & !2-d output
       label,              & !2-d work (int)
       xtvirt,             & !2-d(+1)work(flo
       wtp,wqp,wqac,wparcl,wtc, &!2-d work (flo)
       wqc,wcldcp,wbfrac,wqsat1,wcpqm1, &
       wgpot,wdqtoq,wdqtot, &
       wpsp,wtsp,wcucov,wcpqn,wsdsig, & !1-d work (flo)
       wcuprt,wqacpb)
    !     
    !====================================================================
    !     
    !     ----------------------------------------------------------

    !     l    subroutine kuo:   nwn , bhs ,  dmi ,  december 1987.
    !     l    restructured            pk,   smhi ,  december 1990.
    !     l    optimised         djs , meon , cray , august 1991.
    !     ----------------------------------------------------------
    !     --------
    !     l    purpose:
    !     --------
    !     
    !     l    routine to determine the vertical exchange of heat and
    !     l    moisture through cumulus convection. also precipitation
    !     l    fluxes 'rain' and 'snow' at the ground are determined
    !     l    as well as convective precipitation rate,'prerat'.
    !     
    !     ---------------------------------------------------------
    !     ----------
    !     l    interface:
    !     ----------
    !     
    !     l    subroutine 'kuo' is called from subroutine 'phys'.
    !     
    !     variable        type            content
    !     --------        ------          -------------------------------
    !************************integer input ********************************
    !     l
    !     l      nhor            input           horizontal dimension
    !     l      nlev            input           number of vertical levels
    !     l      kstart          input           startindex for horizontal loops
    !     l      kstop           input           endindex for hor. loops
    !     l
    !************************real input ***********************************
    !     l
    !     l      dtime           input           time step for preliminary
    !     l                                      forecasts of temperatures and
    !     l                                      humidities (from t, q, dtdt and
    !     l                                      dqdt) to be used in the
    !     l                                      condensation calculations
    !************************logical input ********************************
    !     l
    !     l      nlstat          input           switch for statistics
    !     l
    !************************full fields input ****************************
    !     l      t(nhor,nlev)    input           temperatures
    !     l      q(nhor,nlev)    input           specific humidities
    !     l      virt(nhor,nlev)  input          virtual temperature correction
    !     l                                      (made in vdiff)
    !     l      dtdtin(nhor,nlev)  input        "preliminary" temperature tend.
    !     l      dqdtin(nhor,nlev)  input        "preliminary" specific humidity
    !     l                                      tendencies
    !     l      pf(nhor,nlev)   input           pressure of model full levels
    !     l      dpf(nhor,nlev)  input           pressure thicknesses between
    !     l                                      model half levels
    !************************horizontal fields input **********************
    !     l      ps(nhor)        input           ps at tau-1
    !     l      dpsdin(nhor)    input           dps/dt
    !     l      ts(nhor)        input           ts at tau-1
    !     l      dtsdin(nhor)    input           dts/dt
    !     l
    !************************horizontal fields output *********************
    !     l      draindt(nhor)    output         convective rainfall intensity
    !     l
    !     l      dsnowdt(nhor)    output         convective snowfall intensity
    !     l
    !************************horizontal fields input/output ***************
    !     l      prerat(nhor)    input/output    conv. precipitation rate
    !************************1-d vertical arrays input ********************
    !     l      c1er(nlev)      input           defined "first time" in cond
    !************************1-d vertical arrays input output *************
    !     l      stcov(nlev)      input/output    statistics
    !     l      stdcov(nlev)     input/output    statistics
    !     l
    !************************full fields output ***************************
    !     l      dtdt(nhor,nlev) output          temperature tendency as result
    !     l                                      of convective condensation
    !     l      dqdt(nhor,nlev) output          specific humidity tendency as
    !     
    !--------------------------------------------------------------------
    !     
    !     ------------------------------------------------------------
    !     -------
    !     l    method:
    !     -------
    !     
    !     l    see 'hirlam documentation manual ' , 1988.
    !     
    !--------------------------------------------------------------------
    implicit none

    integer nhor,nlev,kstart,kstop

    real(kind=realkind) dtime

    logical nlstat
    real(kind=realkind)     t(nhor,nlev),      q(nhor,nlev), virt(nhor,nlev), &
         dtdtin(nhor,nlev), dqdtin(nhor,nlev),   pf(nhor,nlev),&
         dpf(nhor,nlev),   dtdt(nhor,nlev), dqdt(nhor,nlev),   &
         ps(nhor), dpsdin(nhor),     ts(nhor), dtsdin(nhor),   &
         prerat(nhor), draindt(nhor),dsnowdt(nhor),            &
         c1er(nlev),stcov(nlev),stdcov(nlev)
    !     
    !     -------------------------------------------------------------

    !     local work variables
    !     
    integer jk,jl
    integer label(nhor,nlev)
    real(kind=realkind) wqp(nhor,nlev),wtp(nhor,nlev),xqp,xtp
    real(kind=realkind) wcldcp(nhor,nlev),wqsat1(nhor,nlev),xtvirt(nhor,nlev+1)
    real(kind=realkind) wcpqm1(nhor,nlev),wbfrac(nhor,nlev),wparcl(nhor,nlev)
    real(kind=realkind) wqac(nhor,nlev),wqc(nhor,nlev),wtc(nhor,nlev)
    real(kind=realkind) wgpot(nhor,nlev),wdqtoq(nhor,nlev),wdqtot(nhor,nlev)
    real(kind=realkind) wpsp(nhor),wtsp(nhor)
    real(kind=realkind) wcucov(nhor),wcpqn(nhor),wsdsig(nhor)
    real(kind=realkind) wcuprt(nhor),wqacpb(nhor)
    real(kind=realkind) zdgpot,zicevap,zlvap,zlice,zcpqc
    real(kind=realkind) zqsat,zt,zqcd,ztcm,zdtneg,zdqdtq,zbfrac
    real(kind=realkind) zcuprq,zcpdls,ztpvir,ztcvir,zrevap,zdqdtt
    !     local constants.
    real(kind=realkind) zccov,zcons1,zcrit,zucrit,zeprsi,zdt,zrgrav,zrcpa,zcrdq, z1r1mu,zdti

    real(kind=realkind),parameter::csecur= 1.e-10_realkind

    zccov=600._realkind
    zcons1=cpair*ccpq
    zcrdq=1.0_realkind/epsilo -1.0_realkind
    zcrit=1.0_realkind
    zucrit=0.0_realkind
    z1r1mu=1._realkind/(1._realkind-zucrit)
    zeprsi=epsilo*zcrit
    zdt=dtime
    zdti=1._realkind/zdt
    zrgrav=1._realkind/gravit
    zrcpa=1._realkind/cpair
    !     
    !     ---------------------------------------------------------------

    !     l         2.1    calculate moisture convergence at nlev (loop 215)
    !     l      calculate moisture convergence at remaining levels (loop 218)
    !     l      calculate moisture convergence for whole pbl (loop 218)-
    !     l      if negative assume zero moisture convergence at every layer
    !     l      for the computation of convection  (loop 224)
    !     --------------------------------------------------------------
    !     
    do 215 jl=kstart,kstop
       !     
       !     temporary forward timestepping
       wqp(jl,nlev) = q(jl,nlev) + zdt*dqdtin(jl,nlev)
       wtp(jl,nlev) = t(jl,nlev) + zdt*dtdtin(jl,nlev)
       wpsp(jl) = ps(jl) + zdt*dpsdin(jl)
       wtsp(jl) = ts(jl) + zdt*dtsdin(jl)
       wqac(jl,nlev)=(wqp(jl,nlev)-q(jl,nlev))*dpf(jl,nlev)
       wqacpb(jl)=wqac(jl,nlev)
       wcpqm1(jl,nlev)=cpair + zcons1*q(jl,nlev)
       wcpqn(jl)=wcpqm1(jl,nlev)
       xtvirt(jl,nlev+1)=(1._realkind+zcrdq*q(jl,nlev))*wtsp(jl)
       xtvirt(jl,nlev)=virt(jl,nlev)*wtp(jl,nlev)
       !     
       wgpot(jl,nlev)=log(wpsp(jl)/pf(jl,nlev))
       wgpot(jl,nlev)=rair*xtvirt(jl,nlev)*wgpot(jl,nlev)
       wparcl(jl,nlev)=( wtsp(jl)*wcpqn(jl) - wgpot(jl,nlev) )/wcpqn(jl)
215 enddo
    !     
    do 218 jk=nlev-1,1,-1
       do  jl=kstart,kstop
          wqp(jl,jk) = q(jl,jk) + zdt*dqdtin(jl,jk)
          wtp(jl,jk) = t(jl,jk) + zdt*dtdtin(jl,jk)
          wqac(jl,jk)=(wqp(jl,jk)-q(jl,jk))*dpf(jl,jk)
          wcpqm1(jl,jk)=cpair + zcons1*q(jl,jk)
          xtvirt(jl,jk)=virt(jl,jk)*wtp(jl,jk)
          wgpot(jl,jk)=log(pf(jl,jk+1)/pf(jl,jk))
          !     
          wgpot(jl,jk)=wgpot(jl,jk+1) + rair * wgpot(jl,jk) * &
               0.5_realkind * (xtvirt(jl,jk+1) + xtvirt(jl,jk))
          wparcl(jl,jk)=( wtsp(jl)*wcpqn(jl)-wgpot(jl,jk) )/wcpqm1(jl,jk)
          if( wtp(jl,jk)<=wparcl(jl,jk) ) wqacpb(jl)=wqacpb(jl) + wqac(jl,jk)
       enddo
218 enddo
    !     
    do 224 jk=1,nlev
       do  jl=kstart,kstop
          if(wtp(jl,jk)<=wparcl(jl,jk).and.wqacpb(jl)<=0._realkind) wqac(jl,jk)=0._realkind
       enddo
224 enddo
    !     
    !     --------------------------------------------------------------

    !     l    2.3  specify tc and qc at the surface for cumulus convection
    !     
    do 235 jl=kstart,kstop
       !     
       !     -----------------------------------------------------------
       !     l    23/10/86    initial values for cumulus ascent defined as
       !     l                  environmental values at nlev
       !     
       !     
       wcpqm1(jl,nlev)=wcpqm1(jl,nlev-1)
       if( wqacpb(jl)<=0._realkind) then
          wtc(jl,nlev)=wtsp(jl) - 20._realkind
          wqc(jl,nlev)=0._realkind
       else
          wtc(jl,nlev)=wtp(jl,nlev)
          wqc(jl,nlev)=q(jl,nlev)
       endif
       !     
235 enddo
    !     
    !     -----------------------------------------------------------

    !     l    2.5  calculate tc and qc at next level jk by dry adiabatic
    !     l         lifting and considering latent heat release
    !     l         specify parameter label 0 for stable layers,
    !     l         1 for unstable cloudfree layers, 2 for cloud layers.
    !     -----------------------------------------------------------
    !     
    do 270 jk=nlev-1,1,-1
       do 250 jl=kstart,kstop
          !     
          wqc(jl,jk)=wqc(jl,jk+1)
          zdgpot =wgpot(jl,jk+1) - wgpot(jl,jk)
          wtc(jl,jk)=wtc(jl,jk+1) +  wtc(jl,jk+1)*(1._realkind+zcrdq*wqc(jl,jk+1))* &
               zdgpot /(xtvirt(jl,jk+1) * wcpqm1(jl,jk+1) )
          if( wtc(jl,jk)>wtp(jl,jk) )   then
             label(jl,jk)=1
          else
             label(jl,jk)=0
          endif
          !     
          if (wtc(jl,jk)<tmelt) then
             zicevap = 1.0_realkind
          else
             zicevap = 0.0_realkind
          endif
          !     
          zlvap =latvap

          !     zlvap =latvap + zcons1*wtc(jl,jk)
          zlice =zicevap*latice
          zcpqc=1._realkind/(cpair + zcons1*wqc(jl,jk))
          wcldcp(jl,jk)=( zlvap + zlice) * zcpqc
          !     
          zt=wtc(jl,jk)
          !     itt=zt*100.0
          !     zdif=zdiftabk(itt)
          !     zqsat = zesattabk(itt) / pf(jl,jk)
          zqsat = zeprsi * esat(zt) / pf(jl,jk)
          !     zqcd =( wqc(jl,jk)-zqsat)/(1.+ wcldcp(jl,jk)*zdif/pf(jl,jk))
          zqcd =( wqc(jl,jk)-zqsat)/(1._realkind+ wcldcp(jl,jk)*zeprsi*desdt(zt)/pf(jl,jk))
          zqcd =max(zqcd ,0._realkind)
          wqc(jl,jk)=wqc(jl,jk)-zqcd
          ztcm =wtc(jl,jk)+zqcd * wcldcp(jl,jk)
          !     
          wtc(jl,jk)= ( (wtc(jl,jk+1)*(1._realkind + zcrdq*wqc(jl,jk+1)) +          &
               ztcm *(1._realkind + zcrdq*wqc(jl,jk)) )*                      &
               zdgpot /(xtvirt(jl,jk+1) + xtvirt(jl,jk)) +           &
               (cpair + zcons1*wqc(jl,jk+1))*wtc(jl,jk+1) )*zcpqc +  &
               zqcd *wcldcp(jl,jk)
          !     
          if( wtc(jl,jk)>tmelt) wcldcp(jl,jk)= zlvap *zcpqc
          !     
          if(zqcd >0.0_realkind)  then
             zt=wtc(jl,jk)

             !     itt=zt*100.0
             !     zdif=zdiftabk(itt)
             !     zesat=zesattabk(itt)
             !     
             zqsat =zeprsi*esat(zt)/pf(jl,jk)
             !     zqcd  =  (wqc(jl,jk)-zqsat )/
             !     &         (1.+wcldcp(jl,jk)*zdif/pf(jl,jk))
             zqcd  =  (wqc(jl,jk)-zqsat )/(1._realkind+wcldcp(jl,jk)*zeprsi*desdt(zt)/pf(jl,jk))
             wqc(jl,jk)=wqc(jl,jk) - zqcd
             wtc(jl,jk)=wtc(jl,jk) + zqcd *wcldcp(jl,jk)
          else
             zqcd =0.0_realkind
          endif
          !     
          if ( abs(zqcd) >1.e-14_realkind .and. wtc(jl,jk)>wtp(jl,jk)  .and. &
               wqc(jl,jk)>q(jl,jk))          label(jl,jk)=2
          !     
          if((wtc(jl,jk)-tmelt)>=0.0_realkind ) wcldcp(jl,jk)= &
               zlvap /(cpair + zcons1*wqc(jl,jk))
          !     
          !     -----------------------------------------------------------

          !     l    2.6  check for new lifting level,i.e.
          !     l         specify tc and qc at level jk in case that the layer
          !     l         below was stable and there exists moisture convergence
          !     -----------------------------------------------------------
          !     
          if(label(jl,jk)==0) then
             if( wqac(jl,jk)<=0.0_realkind)  then
                wtc(jl,jk)=wtp(jl,jk)-20.0_realkind
                wqc(jl,jk)=0.0_realkind
             else
                if(jk/=1) then
                   wtc(jl,jk)=wtp(jl,jk)
                   wqc(jl,jk)=q(jl,jk)
                endif
             endif
          endif
250    enddo
270 enddo
    !     
    !     -------------------------------------------------------------

    !     l    3.   set label=0 for dry unstable layers if no cloud is above
    !     l         set label=3 for lifting level if layer above is unstable
    !     
    do 320 jk=2,nlev
       do  jl=kstart,kstop
          if( label(jl,jk)==1.and.label(jl,jk-1)==0) label(jl,jk)=0
       enddo
320 enddo

    do 340 jk=nlev,2,-1
       do  jl=kstart,kstop
          if( label(jl,jk)==0 .and.label(jl,jk-1)/=0) label(jl,jk)=3
       enddo
340 enddo

    !     l    3.5  calculate total moisture accession for unstable layers
    !     l         and moisture/temperature surplus of clouds.

    do 360 jl=kstart,kstop
       if( wtp(jl,nlev)<=wparcl(jl,nlev) ) then
          label(jl,nlev) = 1
       else
          label(jl,nlev) = 0
       endif
       zdtneg  = (tmelt-wtc(jl,nlev))/15._realkind
       if( zdtneg>0.0_realkind) then
          zlice  = latice*min(zdtneg*zdtneg*zdtneg,1.0_realkind)
       else
          zlice  = 0.0_realkind
       endif
       wcldcp(jl,nlev)=(latvap+ zlice)*zrcpa
       wqsat1(jl,nlev)=zeprsi*esat(wtp(jl,nlev))/pf(jl,nlev)
       wdqtot(jl,nlev)=0._realkind
       wdqtoq(jl,nlev)=0._realkind
       wbfrac(jl,nlev)=0._realkind
       wsdsig(jl)=0._realkind
       if( label(jl,nlev)==0)  wqac(jl,nlev) = 0.0_realkind
360 enddo

    do 365 jk=nlev-1,1,-1
       do  jl=kstart,kstop
          if( (label(jl,jk)/=3).and. (label(jl,jk+1)>0) ) &
               wqac(jl,jk)=wqac(jl,jk+1) + wqac(jl,jk)
          if( label(jl,jk)==0 ) wqac(jl,jk)=0._realkind
          wqsat1(jl,jk)=zeprsi*esat(wtp(jl,jk))/pf(jl,jk)
          if(label(jl,jk)==2)then
             zdqdtt=  ( wtc(jl,jk)*(1._realkind+zcrdq*wqc(jl,jk)) - &
                  xtvirt(jl,jk) )*  dpf(jl,jk)/wcldcp(jl,jk)
             wdqtot(jl,jk)=wdqtot(jl,jk+1) + zdqdtt
             zdqdtq=( wqsat1(jl,jk)-q(jl,jk) )*dpf(jl,jk)
             wdqtoq(jl,jk)=wdqtoq(jl,jk+1) + zdqdtq
             zbfrac=(1._realkind-q(jl,jk)/wqsat1(jl,jk))*z1r1mu
             wbfrac(jl,jk)=( wbfrac(jl,jk+1)*wsdsig(jl) + &
                  zbfrac*dpf(jl,jk) )/( wsdsig(jl)+dpf(jl,jk) )
             wsdsig(jl)=wsdsig(jl) + dpf(jl,jk)
          else
             wdqtot(jl,jk)=0._realkind
             wdqtoq(jl,jk)=0._realkind
             wbfrac(jl,jk)=0._realkind
             wsdsig(jl)=0._realkind
          endif
       enddo
365 enddo

    !     l       3.8  replace the moisture accession at cloud layers by the
    !     l            total moisture accession of the whole unstable layer
    !     l            do the same for the moisture + temperature surplus
    !     l            set label= 0 if total moisture accesion is le.0.
    !     -------------------------------------------------------------
    !     
    do 382 jl=kstart,kstop
       if( (wqac(jl,1)<=0._realkind).or.(wdqtot(jl,1)<=0._realkind).or.&
            (wdqtoq(jl,1)<=0._realkind) )   label(jl,1)=0
       prerat(jl)=0._realkind
       wcucov(jl)=0._realkind
382 enddo

    do 388 jk=2,nlev
       do jl=kstart,kstop
          if( (label(jl,jk)==2).and.(label(jl,jk-1)==2) ) then
             wqac(jl,jk)=wqac(jl,jk-1)
             wdqtot(jl,jk)=wdqtot(jl,jk-1)
             wdqtoq(jl,jk)=wdqtoq(jl,jk-1)
             wbfrac(jl,jk)=wbfrac(jl,jk-1)
          endif

          if((wdqtot(jl,jk)<=0._realkind.or.wdqtoq(jl,jk)<=0._realkind &
               .or.wqac(jl,jk)<=0._realkind).and.(label(jl,jk)==2)) &
               label(jl,jk)=0

          if((label(jl,jk)/=2).and.(label(jl,jk-1)==0))label(jl,jk)=0
       enddo
388 enddo

    !     l       4.   do mixing of the cumulus air with the environment
    do 450 jk=1,nlev
       do  jl=kstart,kstop
          if( label(jl,jk)==2.and. abs(wqac(jl,jk))>1.e-14_realkind) then
             !     l       the moistening parameter wbfrac(jl,jk) is decreased
             !     l       by computing wbfrac(jl,jk)**3 . ( ecmwf 1985 ).
             wbfrac(jl,jk)=wbfrac(jl,jk)*wbfrac(jl,jk)*wbfrac(jl,jk)
             if( wdqtot(jl,jk)<=0._realkind )   wdqtot(jl,jk)=-1.0_realkind
             if( wdqtoq(jl,jk)<=0._realkind )   wdqtoq(jl,jk)=-1.0_realkind
             wcuprt(jl)=(1._realkind-wbfrac(jl,jk)) * wqac(jl,jk)/wdqtot(jl,jk)
             zcuprq=wbfrac(jl,jk) * wqac(jl,jk)/wdqtoq(jl,jk)
             zcpdls=1._realkind/wcldcp(jl,jk)
             ztpvir=wtp(jl,jk)*( 1._realkind+zcrdq*q(jl,jk) )
             ztcvir=wtc(jl,jk)*( 1._realkind+zcrdq*wqc(jl,jk) )
             prerat(jl)=prerat(jl) + &
                  wcuprt(jl)*(ztcvir -ztpvir)*dpf(jl,jk)*zrgrav*zcpdls
             xtp=wtp(jl,jk)+ wcuprt(jl)*( ztcvir-ztpvir )
          else
             wcuprt(jl)=0._realkind
             xtp=wtp(jl,jk)
             zcuprq=0._realkind
          endif
          if( label(jl,jk)>0 )then
             xqp=( q(jl,jk)*ps(jl)/wpsp(jl) ) &
                  + zcuprq*( wqsat1(jl,jk)-q(jl,jk) )
          else
             xqp=wqp(jl,jk) + zcuprq*( wqsat1(jl,jk)-q(jl,jk) )
          endif

          !     l       4.2  calculate evaporation of rain
          if(prerat(jl)>0.0_realkind) then
             zqsat =zeprsi*esat(xtp)/pf(jl,jk)
             wcucov(jl)= max( wcucov(jl),wcuprt(jl)*zccov + csecur )

             if ( (label(jl,jk)/=2).and.((zqsat -wqp(jl,jk))>0._realkind))then
                zrevap =(zqsat -wqp(jl,jk))* &
                     wcucov(jl)*c1er(jk)*sqrt(prerat(jl)/wcucov(jl))
             else
                zrevap =0._realkind
             endif

             zqcd =prerat(jl)*gravit/dpf(jl,jk)
             zrevap =min(zqcd , zrevap )
             xqp=xqp +  zrevap
             xtp=xtp -  zrevap *wcldcp(jl,jk)

             prerat(jl)=prerat(jl) - zrevap *dpf(jl,jk)*zrgrav
          else
             prerat(jl)=0.0_realkind
          endif
          dqdt(jl,jk)=( xqp-wqp(jl,jk) ) / zdt
          dtdt(jl,jk)=( xtp-wtp(jl,jk) ) / zdt
       enddo
450 enddo

    !     l       5.1  calculate precipitation and precipitation fluxes
    !     l            'rain', 'snow' .
    do 515 jl=kstart,kstop
       prerat(jl)=prerat(jl)*zdti
       if( ts(jl)>=tmelt )   then
          draindt(jl)=draindt(jl) + prerat(jl)
       else
          dsnowdt(jl)=dsnowdt(jl) + prerat(jl)
       endif
515 enddo

    !     l    update 'stcov(jk)' and 'stdcov(jk)'
    if(nlstat)then
       do 2500 jk=1,nlev
          do 2400 jl=kstart,kstop
             stcov(jk)=stcov(jk) + dtdt(jl,jk)
             stdcov(jk)=stdcov(jk) + dtdt(jl,jk)*dtdt(jl,jk)
2400      enddo
2500   enddo
    endif
    return
  end subroutine kuo


end module kuomod
