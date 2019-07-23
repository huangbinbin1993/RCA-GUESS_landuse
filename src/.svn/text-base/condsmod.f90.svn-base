module condsmod
  use confys
  use decomp
  use escom
  use referenceParameters,only:sip0

  implicit none
  private
  real(kind=realkind),public,save::he273= 611.0_realkind
  real(kind=realkind),public,save::tvirtc= 0.61_realkind
  real(kind=realkind),public,save::u00max = 0.975_realkind
  real(kind=realkind),public,save::hu00 = 0.85_realkind
  real(kind=realkind),public,save::aecon
  real(kind=realkind),public,save::conae= 0.15_realkind
  real(kind=realkind),public,save::coales= 100.0_realkind
  real(kind=realkind),public,save::hccu= 2.5e-4_realkind
  real(kind=realkind),public,save::hcunrm = 3.e4_realkind
  real(kind=realkind),public,save::hmrcu = 1.5e-3_realkind
  real(kind=realkind),public,save::hmrst = 5.0e-4_realkind
  real(kind=realkind),public,save::htaucu = 3600.0_realkind
  real(kind=realkind),public,save::hp0 = 1.e5_realkind
  real(kind=realkind),public,save::cbfeff = 4.0_realkind
  real(kind=realkind),public,save::hvterm= 5.0_realkind
  real(kind=realkind),public,save::hvsnow= 1.0_realkind
  real(kind=realkind),public,save::hkmelt= 5.e-5_realkind
  real(kind=realkind),public,save::tanvil= 253.0_realkind
  real(kind=realkind),public,save::hcst = 1.e-4_realkind
  real(kind=realkind),public,save::pmoist = 3.0_realkind
  real(kind=realkind),public,save::prbice(1750)
  real(kind=realkind),public,save::bfeff(1750)
  real(kind=realkind),public,save::hdewi(1750)
  real(kind=realkind),public,save::hmroft(1750)
  real(kind=realkind),public,save::dttabl = 0.1_realkind
  real(kind=realkind),public,save::hdl=0.334e6_realkind
  real(kind=realkind),public,save::ht273= 273.0_realkind
  real(kind=realkind),public,save::stpevp = 5.e-5_realkind
  real(kind=realkind),public,save::cfreez = 0.12_realkind
  real(kind=realkind),public,save::hkevap = 5.0E-05_realkind
  real(kind=realkind),public,save::tcir= 235.0_realkind
  real(kind=realkind),public,save::hkap   !     hkap   - rair/cpair
  real(kind=realkind),public,save::hdldcp !     hdldcp - latice/cpair
  real(kind=realkind),public,save::hecdr  !     hecdr  - epsilo / hkap
  real(kind=realkind),public,save::hldcp  !     hldcp  - latvap/cpair
  real(kind=realkind),public,save::coalcu
  real(kind=realkind),public,save::coalst
  real(kind=realkind),public,save::sqvsno
  real(kind=realkind),public,save::heldr
  real(kind=realkind),public,save::elotci !epsilo*latice/rair/tcir
  real(kind=realkind),public,save::asnow = 9.0_realkind
  real(kind=realkind),public,save::bsnow = 4.0_realkind
  real(kind=realkind),public,save::snoref = 0.05_realkind

  public aconds, iniconds,inicons
contains
  subroutine iniconds()
    implicit none

    integer i10,it,i
    real(kind=realkind) hedldr,tci,topeq0,todpmx,tscale,demax,apri,test, &
         tphalf,podpmx,fxt,xti,x,x232,bfmax,xpt1mp

    !     processor no of the current PE and function to get it


    !     Initialize variables used in common block 'COMCONDS'
    !     used in connection with the STRACO cloud scheme.


    heldr  = epsilo*latvap/rair
    hedldr = epsilo*hdl/rair

    tci    = 232.0_realkind
    topeq0 = 273.0_realkind
    todpmx = 299.0_realkind
    tscale = (todpmx - tci)*sqrt(2.0_realkind)

    !     the table look up is done by finding the integer
    !     it = 1 + nint((t- h273)/dttabl)
    !     
    dttabl = 0.1_realkind

    !     calculate saturation vapor pressure to be used as table
    !     
    !     calculate difference of saturation pressure between
    !     water and ice, dewi, for t< t0
    !     calculate a probability function, prbice, for existence of
    !     ice crystals
    !     calculate a temperature function, hmroft,
    !     to multiply hmrcu and hmrst
    !     
    demax = 0.0_realkind
    apri = 1.0_realkind/(1.0_realkind-exp(-((topeq0-tci)/tscale)**2.0_realkind))
    test = 2.0_realkind*apri/(2.0_realkind*apri-1.0_realkind)
    if(test <= 0.0_realkind) write(6,*)'warning: test is < 0 !!!!'
    tphalf = log(2.0_realkind*apri/(2.0_realkind*apri-1.0_realkind))
    tphalf = tci + tscale*sqrt(tphalf)
    if(tphalf<0.0_realkind) write(6,*) 'before do 25,  tphalf=',tphalf
    podpmx = apri*(exp(-0.5_realkind)-1.0_realkind)+1.0_realkind

    do 25 i = 1,1750
       fxt = ht273 - (real(i,realkind)-1.0_realkind) * dttabl

       xti = 1.0_realkind/ht273 - 1.0_realkind/fxt
       hdewi(i)  = he273/fxt * exp(heldr*xti) * (1.0_realkind - exp(hedldr*xti))
       demax = max(demax,hdewi(i))

       prbice(i) = apri * (exp(-(((fxt - tci)/tscale)**2.0_realkind)) - 1.0_realkind) &
            + 1.0_realkind
       prbice(i) = max( prbice(i) , 0.0_realkind)
       if(fxt < tci) prbice(i) = 1.0_realkind
       prbice(i) = min(prbice(i) , 1.0_realkind)
       if(fxt > 250.0_realkind) hmroft(i) = 1.33_realkind*&
            exp(-((fxt-ht273)*0.066_realkind)**2.0_realkind)

       if(fxt <= 250.0_realkind) then
          x = abs (fxt - 232.0_realkind) / 18.0_realkind

          x = x * (1.0_realkind + x * (1.0_realkind+ 1.333_realkind * x))
          x232 = 1.0_realkind
          if (fxt < 232.0_realkind) x232 = -1.0_realkind

          x = x / (1.0_realkind + x) * x232

          hmroft(i) = 0.5_realkind*0.15_realkind*(1.07_realkind+x)
       endif
       hmroft(i) = min( hmroft(i) , 1.0_realkind)
       hmroft(i) = (1.0_realkind-prbice(i))**2.0_realkind + prbice(i) * hmroft(i)
       hmroft(i)   = max(hmroft(i) , 3.e-2_realkind)
25  enddo
    prbice(1) = 0.0_realkind
    hmroft(1) = 1.0_realkind

    !     normalize with max difference, demax
    !     calculate the product of dewi and prbice, called bfeff

    bfmax = 0._realkind

    do 30 i = 1,1750

       hdewi(i)  = hdewi(i)/demax
       bfeff( i) = prbice(i) * (1.0_realkind-prbice(i)) * hdewi(i)
       bfmax = max(bfmax,bfeff(i))
30  enddo

    do 32 i = 1,1750
       bfeff( i) = bfeff( i) / bfmax
32  enddo

    if(mype==0) then
       write(6,600)
       write(6,605) topeq0,tphalf,podpmx,todpmx
       write(6,615)
       do 60 i10 = 1,47,3
          it = 1 - i10
          i = 1 + (i10-1) * 10
          xpt1mp = prbice(i)*(1.0_realkind-prbice(i))
          write(6,617) it, hdewi(i),prbice(i),bfeff(i), hmroft(i),xpt1mp
60     enddo
    endif

600 format(//,' from iniconds:',/)


605 format('   PROBABILITY TABLES BETWEEN 273 AND 99 K FOR:',//,        &
         '     DIFFERENCE of saturation pressure over WATER and ICE;',  &
         /,                                                             &
         '     TEMPERATURE FUNCTION to multiply hmrst and hmrcu;',//,   &
         '     ICE CRYSTAL PROBABILITY P = 0 for T = ',f5.1,            &
         '  and P = 0.50                                                &
         for T = ',f5.1,/,20x,'     and P and T of dP/dT-max are P ='   &
         ,f5.2,'  and T = ',f5.1,/)

615 format(9x,'t',5x,'de',6x,'prob-ice',5x,'bf-f',4x,'hmroft', &
         5x,'p*(1-p)')
617 format(i10,f9.4,4x,f7.4,4x,f7.4,2x,f7.4,4x,f7.4)

    !     parameter values in si units
    conae = 0.15_realkind
    aecon = exp(conae)

    cfreez = 0.12_realkind
    coales = 100.0_realkind
    hccu = 2.5e-4_realkind
    hcunrm = 3.e4_realkind
    hmrcu =1.5e-3_realkind
    hmrst =5.0e-4_realkind
    htaucu = 3600.0_realkind
    tanvil = 253.0_realkind
    hcst = 1.e-4_realkind

    cbfeff = 4.0_realkind
    hvterm = 5.0_realkind
    hvsnow = 1.0_realkind
    hkmelt = 5.e-5_realkind

    stpevp = 5.e-5_realkind
    u00max = 0.975_realkind
    hu00 = 0.85_realkind
    pmoist = 3.0_realkind

    !     lr for stamic et al.:

    hkevap= 5.0E-05_realkind
    tcir= 235.0_realkind

    !     asnow  - a parameter together with bsnow and snoref to govern the
    !     rate of ice precip before the Bergeron-Findeisen mechanism
    !     becomes effective, even if the precip from above is pure
    !     ice;
    !     the factor cbfsno = asnow*(snowrate/snoref)**bsnow
    !     and cbfsno = cbfsno/(1+cbfsno)
    !     which multiplies the modified ice probability
    !     bsnow  - see asnow
    !     snoref - see asnow

    asnow = 9.0_realkind
    bsnow = 4.0_realkind
    snoref = 0.05_realkind

    !     hkmelt - coefficient for melting of ice in precipitation
    !     tanvil - temperature below which convective anvil is allowed and
    !     treated as stratiform
    !     hcst   - conversion rate of cloud to precipitation drops in stratifor
    !     clouds
    !     pmoist - exponent for computations of hmoist (eq. 3 in OSZ94)
    !     
    hkmelt = 2.5e-4_realkind
    tanvil = 253.0_realkind
    hcst = 1.e-4_realkind
    pmoist = 3.0_realkind
    !     

    !     
    hkap   = rair/cpair
    hdldcp = latice/cpair
    hecdr  = epsilo / hkap
    hldcp  = latvap/cpair
    !     
    !     coalcu - coalescence factor for convection clouds (c0.mr/2g)
    !     coalst - coalescence factor for convection clouds (c0.mr/2g)
    !     sqvsno - square root of the snow terminal velocity
    !     elotci - epsilo*latice/rair/tcir
    !     
    coalcu = 0.5_realkind* hccu* hmrcu / gravit
    coalst = 0.5_realkind* hcst* hmrst / gravit
    sqvsno = sqrt(hvsnow)
    elotci = epsilo*latice/rair/tcir
    !     
    !     print some of the parameters
    !     
    if(mype==0) then
       !     
       write(6,3) tanvil,hu00,htaucu,hcunrm,hccu,hmrcu,cfreez,coales,hkmelt,pmoist
3      format(1h ,//,'  with ice phase,'                                 &
            ,/,'  convection scheme                                     &
            for T(cutop) <=',f5.0,' K',//,' other parameters are',/,  &
            '   u0   taucu    cunrm     cuc0' //                        &
            '   cumr    cfreez  coales   hkm                            &
            elt    pmoist',                                             &
            /,f6.2,f7.0,f10.0,e10.2,e10.2,f7.2,f8.0,e11.2,f7.1,//)
       !     
       write(6,4) hu00,stpevp,cbfeff,hcst,hmrst,cfreez,coales,hkmelt,pmoist
4      format(1h ,//,'ice phase  stratiform cloud of sub-grid scale, i.e., even for u < 1.',&
            //,' other parameters are   u0     kevp    cbfeff    stc0      stmr    cfreez  coales   hkmelt    pmoist'              &
            , /,f6.2,e10.2,f6.1,2x,e9.2,e10.2,f7.2,f8.0,e11.2,f7.1,//)
       !     
    endif
    !     
    return
  end subroutine iniconds


  subroutine aconds(kcall,nhor,nlev,              &
       kstart,kstop,kstep, &
       dtime,conacc,         &                                   
       dtheta, &
       t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf, &
       ps,dpsdin,                           &           
       cov2d,cwpath,draindt,dsnowdt,  &
       preta,prcpst,prcpcu,stsnow,cusnow, &            
       ahyb,bhyb, &
       dtdt,dqdt,dcwdt,totcov,cucov)                          
    implicit none                                                           
    integer kcall,nhor,nlev,kstart,kstop,kstep
    real(kind=realkind) dtime,conacc,dtheta

    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),dcwdin(nhor,nlev), &        
         dtdtin(nhor,nlev),dqdtin(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev)       
    real(kind=realkind) ps(nhor),dpsdin(nhor),    &                                
         cov2d(nhor),cwpath(nhor),dsnowdt(nhor),  &
         draindt(nhor), prcpst(nhor),prcpcu(nhor),stsnow(nhor),cusnow(nhor),&
         preta(nhor,nlev)
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1) !, stcon(nlev),stdcon(nlev), &
    !stcov(nlev),stdcov(nlev),stccov(nlev),          &                      
    !stscov(nlev),stcw(nlev),stcpnt(nlev),sttpnt(nlev),stscal(nlev)        
    real(kind=realkind) dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev),  &                
         totcov(nhor,nlev),cucov(nhor,nlev)                                   
    !                                                                             
    !

    !     
    !     The routine is called from suboutine 'PHYS'
    ! work space:                                                                 
    !                                                                             

    !  here follows 1-D integer work arrays                                       
    !
    !  here follows 2-D integer work array                                        
    integer klabcv(nhor,nlev),klabc(nhor,nlev),klab(nhor,nlev)
    !
    !  here follows 1-D real work arrays                                          
    real(kind=realkind) wrhuc(nhor,nlev), wprst(nhor),wprcu(nhor), &
         wcusno(nhor),wstsno(nhor), wmr(nhor),wpcoef(nhor), &
         wdcusn(nhor),wdstsn(nhor),wpsp(nhor),wqsat(nhor),wdqsdt(nhor), &
         wqmax(nhor),wrk1(nhor),wrk2(nhor) 
    !  here follows 2-D real work arrays                                          
    real(kind=realkind) wldcp(nhor,nlev),wqp(nhor,nlev),wtp(nhor,nlev),wcwp(nhor,nlev),&
         wtc(nhor,nlev),wqc(nhor,nlev),wrhu(nhor,nlev),&
         wqliq(nhor,nlev),wtcdif(nhor,nlev),wqcdif(nhor,nlev),wqac(nhor,nlev),&
         wdptot(nhor,nlev),wdpcld(nhor,nlev),wbfrac(nhor,nlev),&
         wdift(nhor,nlev),wdifq(nhor,nlev)
    real(kind=realkind) wbu(nhor,nlev),wgpot(nhor,nlev)  
    real(kind=realkind) wscdif(nhor,nlev),wdifs(nhor,nlev)

    !     Call the main routine 'conds' of the STRACO
    !     cloud-and condensation scheme 
    call conds (kcall,nhor,nlev,kstart,kstop,kstep,&
         dtime,conacc,                      &                     
         dtheta,&
         t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,    &
         ps,dpsdin,                              &       
         cov2d,cwpath,               &
         draindt,dsnowdt, &
         preta,prcpst,prcpcu,stsnow,cusnow,&
         ahyb,bhyb, &
         !stcon,stdcon,stcov,stdcov,stccov,  &
         !stscov,stcw,stcpnt,sttpnt,stscal,     &          
         dtdt,dqdt,dcwdt,totcov,cucov,       &
         wrhuc,&
         wdcusn,wdstsn,&
         wmr,wpcoef,wpsp,wrk1,wrk2, &
         klabcv,klabc,klab, &
         wldcp,wqp,wtp,wcwp, &
         wtc,wqc,wqsat,wdqsdt,wrhu, &
         wtcdif,wqcdif,wqac,wdptot,wdpcld, &
         wbfrac,wdift,wdifq,wbu,wgpot, &
         wscdif,wdifs )
    return                                                                  
  end subroutine aconds

  subroutine conds (kcall,nhor,nlev,               &
       kstart,kstop,kstep,                          &
       dtime,conacc,                                &
       dtheta,                                      &
       t,q,cw,dcwdin,dtdtin,dqdtin,pf,dpf,          &
       ps,dpsdin,                                   &
       cov2d,cwpath,                                &
       draindt,dsnowdt,                             &
       preta,prcpst,prcpcu,stsnow,cusnow,           &
       ahyb,bhyb, &
       !stcon,stdcon,stcov,stdcov,stccov,  &
       !stscov,stcw,stcpnt,sttpnt,stscal,            &
       dtdt,dqdt,dcwdt,totcov,cucov,                &
       wrhuc,                                       &
       wdcusn,wdstsn,                               &
       wmr,wpcoef,wpsp,wrk1,wrk2,                   &
       klabcv,klabc,klab,                           &
       wldcp,wqp,wtp,wcwp,                          &
       wtc,wqc,wqsat,wdqsdt,wrhu,                   &
       wtcdif,wqcdif,wqac,wdptot,wdpcld,            &
       wbfrac,wdift,wdifq,wbu,wgpot,                &
       wscdif,wdifs )
    !     


    !     Purpose:     
    !     --------   
    !     
    !     Parameterize the effect of condensation and precipitation processes.
    !     A semi-implicit numerical treatment of the precipitation
    !     related processes has been implemented.
    !     The scheme is called 'STRACO' (Soft TRAnsition COndensation)
    !     
    !     i) 
    !     First, the moist convective processes  are treated ( subroutine
    !     'condcv' )
    !     As a first step, preparations for the convective parame-
    !     terization are made by assigning labels and parameters used
    !     in the convection scheme which is an elaborate 'KUO' type 
    !     of scheme. Convective 'entities' can be diagnosed from 
    !     any level in the atmosphere, provided that the necessary 
    !     closure conditions are fulfilled. 
    !     Convective updates of temperature, specific humidity and
    !     cloud water are carried out. The preliminary updates due to 
    !     dynamics, vertical diffusion and radiation are 
    !     taken into account. 
    !     
    !     ii)
    !     Cloud cover is computed in any model layer. A distinction 
    !     is made between 'current' and 'equlibrium' cloud cover at 
    !     a given instant.  Also a distinction is made between 
    !     convective and stratiform cloud cover. 
    !     Drastic changes in cloud cover associated 
    !     with oscillations between stratiform and convective regime are
    !     prevented by formulating gradual transitions between the 
    !     two regimes.
    !     
    !     iii) 
    !     Stratiform condensation step ( sub-routine 'condst' )
    !     Sub-grid scale condensation starts below a 100 percent grid box
    !     relative humidity. Stratiform condensation is done always if
    !     the grid box relative humidity exceeds 100 % after other model
    !     processes have been active. Otherwise, the stratiform condensation
    !     is not done inside the convective entities.
    !     
    !     iv)                                     
    !     precipitation release and evaporation of precipitation
    !     is considered next, based on the updated values from previous 
    !     processes (subroutine 'prevap'). 
    !     The original 'microphysics' parameterizations as described 
    !     in the HIRLAM 2.5 Documentation Manual (June 1996, p 2.44 -2.49) 
    !     has been extended to include the work of Sundqvist (1993)
    !     (Inclusion of ice phase of hydrometeors in cloud parameterization
    !     for mesoscale and large scale models (Beitr. Phys. Atmosph., 66,
    !     137-147 )
    !     
    !     v) 
    !     Finally, some diagnostics are computed related to codensation 
    !     and precipitation processes. 
    !     
    !     Interface:    subroutine 'conds' called from subroutine 'aconds' 
    !     ----------
    !     
    !     Subroutines:  
    !     ------------ 
    !     'condcv'     ( computes heating and moistening due to 
    !     convection )
    !     'cloudcv'    ( computes cloud cover in all model layers )
    !     
    !     'condst'     ( computes heating and moistening due to
    !     stratiform condensation )
    !     'prevap'     ( computes precipitation release and evaporation
    !     of precipitation )
    !     
    !     Author:        Bent H. Sass, (DMI), based on the research during the
    !     -------        the HIRLAM I-III projects. 
    !     
    !     Literature:    HIRLAM Newsletter No. 29, November 1997, p 37-45   
    !     ------------   ( Available from Met Eireann, Dublin, Ireland ) 
    !     II) (Beitr. Phys. Atmosph., 66, 137-147 )
    !     

    !     1. Declarations     
    !     ---------------
    !     
    implicit none                                                           

    integer kcall,nhor,nlev,kstart,kstop,kstep
    real(kind=realkind) dtime,conacc,dtheta              
    integer jl,jk
    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),dcwdin(nhor,nlev), &        
         dtdtin(nhor,nlev),dqdtin(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev) 
    real(kind=realkind) ps(nhor),dpsdin(nhor),                               &
         cov2d(nhor),cwpath(nhor),draindt(nhor),dsnowdt(nhor), &     
         prcpst(nhor),prcpcu(nhor),stsnow(nhor),cusnow(nhor),  &
         preta(nhor,nlev)
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1) !,stcon(nlev),stdcon(nlev),  &   
    !stcov(nlev),stdcov(nlev),stccov(nlev),                &
    !stscov(nlev),stcw(nlev),stcpnt(nlev),sttpnt(nlev),stscal(nlev)
    real(kind=realkind) dtdt(nhor,nlev),dqdt(nhor,nlev),dcwdt(nhor,nlev),  &               
         totcov(nhor,nlev),cucov(nhor,nlev)                         
    !     
    !     1.1  work space:
    !     ----------------
    !     here follows 1-D integer work arrays
    !     
    !     here follows 2-D integer work array
    !     
    integer klabcv(nhor,nlev),klabc(nhor,nlev),klab(nhor,nlev)
    !     
    !     here follows 1-D real work arrays
    !     
    real(kind=realkind) wmr(nhor),wpcoef(nhor), wdcusn(nhor),wdstsn(nhor), &
         wpsp(nhor),wqsat(nhor),wdqsdt(nhor), wqmax(nhor),&
         wrk1(nhor),wrk2(nhor)  
    !     
    !     here follows 2-D real work arrays
    !     
    real(kind=realkind) wldcp(nhor,nlev),wrhuc(nhor,nlev),                   &
         wqp(nhor,nlev),wtp(nhor,nlev),wcwp(nhor,nlev),        &
         wtc(nhor,nlev),wqc(nhor,nlev),wrhu(nhor,nlev),        &
         wqliq(nhor,nlev),                                     &
         wtcdif(nhor,nlev),wqcdif(nhor,nlev),wqac(nhor,nlev),  &
         wdptot(nhor,nlev),wdpcld(nhor,nlev),wbfrac(nhor,nlev),&
         wdift(nhor,nlev),wdifq(nhor,nlev) 
    real(kind=realkind) wbu(nhor,nlev),wgpot(nhor,nlev)
    real(kind=realkind) wscdif(nhor,nlev),wdifs(nhor,nlev)
    !     
    !     local variables 
    !     
    real(kind=realkind) zwrk1,zwrk2,zwrk3,zwrk4,zwrk5,zwrk6, zwrk7,zwrk8
    !-----------------------------------------------------------------------

    !     1.2)  Nomenclature 
    !     ------------------------------------------------------------      
    !     
    !     1.2.1 Arrays. 
    !     -------------
    !     
    !     bfeff(1750) bergeron-findeisen effect from (dewi*prbice)              
    !     
    !     cov2d(jl)    total cloud cover (2-dimensional) 
    !     
    !     cw(jl,jk)   cloud water content at (t-dt)                             
    !     
    !     dpf(jk)      delta-p around level jk                                  
    !     
    !     dqdt(jk)      tendency of humidity due to condensation processes
    !     dqdtin(jl,jk) tendency of humidity due to effects other than          
    !     condensation                                             
    !     dtdt(jl,jk)   tendency of temperature due to convection                
    !     dtdtin(jl,jk) tendency of temperature due to effects other than       
    !     condensation                                             
    !     cucov(jk)     contains cumulus cloud cover                             
    !     hdewi(1750)   difference in saturation vapour pressure over            
    !     water and ice                                            
    !     wldcp(jl,jk)  (1/cpair)*latent heat as a function of                  
    !     temperature by weighting in latent heat of               
    !     freezing times probability of ice crystals               
    !     hmroft(1750)  temperature function to multiply mr with                 
    !     for t<273                                                
    !     wqsat(jl,jk)  saturation mixing ratio with respect to                  
    !     gridpoint temperature                                    
    !     wdqsdt(jl,jk) derivative of saturation water vapour pressure 
    !     with respect to temperature 
    !     pf(jk)        pressure at full level
    !     prbice(1750)  probability for ice crystals as a function of            
    !     temperature                                              
    !     wprcu(jl)     rate of convective precipitation at level jk,            
    !     wprst(jl)     rate of stratiform precipitation at level jk,            
    !     stcon(jk)     mean temperature tendency due to condensation for        
    !     each model layer                                         
    !     stdcon(jk)    mean variance of temperature tendency due to             
    !     condensation for each model layer                        
    !     stccov(jk)    mean convective cloud cover for each model layer         
    !     stscov(jk)    mean   total    cloud cover for each model layer         
    !     stcpnt(jk)    number of points with convective condensation            
    !     for each model layer                                     
    !     sttpnt(jk)    total number of points with condensation                 
    !     for each model layer                                     
    !     stscal(1)     volume mean cloud cover from assumption of               
    !     maximum/random overlapping                               
    !     stscal(2)     volume mean vertically integrated cloud water            
    !     stscal(3)     average accumulated convective precipitation             
    !     stscal(4)     average accumulated stratiform precipitation             
    !     stscal(5)     average accumulated amount of water from                 
    !     evaporation of precipitation                             
    !     stscal(6)     average accumulated amount of water from                 
    !     evaporation of clouds                                    
    !     stscal(7)     total number of convective columns                       
    !     stscal(8)     total number of points with condensation in the          
    !     model volume                                             
    !     stscal(9)     volume mean convective cloud cover from                  
    !     assumption of maximum overlapping                        
    !     stscal(10)    average convective precipitation rate                    
    !     stscal(11)    average stratiform precipitation rate                    
    !     
    !     klabcv       label to identify convective cloud layers
    !     klabc        label to identify convective cloud layers (1 or 0)
    !     klab         label to identify bottom vertical level of new 
    !     convective 'entity', It contains the number 
    !     of the convective entity. 
    !     wstsno       stratiform precipitation rate as snow
    !     wcusno       convective precipitation rate as snow 
    !     wpsp         provisional value of surface pressure 
    !     wqp          provisional value of specific humidity 
    !     wtp          provisional value of temperature 
    !     wcwp         provisional value of specific cloud water
    !     wtc          cloud parcel temperature 
    !     wqc          cloud parcel specific humidity
    !     wqsat        saturation specific humidity 
    !     wdqsdt       derivative of saturation specific humidity with
    !     respect to temperature 
    !     wrhu         relative humidity 
    !     wrhuc        relative humidity threshold 
    !     wdift, 
    !     wtcdif       weighting parameters for temperature change
    !     in the convection scheme 
    !     wdifq,
    !     wqcdif       weighting parameters for humidity change 
    !     in the convection scheme 
    !     wdifs,
    !     wscdif       weighting parameters for cloud condensate
    !     change in the convection scheme
    !     wqac         vertically integrated moisture convergence
    !     convective 'entity'.
    !     wdptot       total pressure thickness of convective 'entity'
    !     wdpcld       total pressure thickness of convective 
    !     'cloud entity'
    !     wbfrac       moisture partitioning parameter in the
    !     convection scheme
    !     wbu          averarage buoyancy '(tc -t)/t)' in the               
    !     convective 'entity' 
    !     
    !     1.2.2) scalars                                                               
    !     --------------                                                               
    !     cfreez      increases conversion rate for temp low than 273 k         
    !     coales      increases conversion rate due to precipitation            
    !     coming in from above                                      
    !     dttabl      increment in temperature dependent tables                 
    !     cpair       specific heat of dry air at constant pressure           
    !     hccu        conversion rate of cloud to precip drops in               
    !     convective cloud                                          
    !     hcst        conversion rate of cloud to precip drops in               
    !     stratiform cloud                                          
    !     hdl         latent heat of melting                                    
    !     heldr       =epsilo*latvap/rair
    !     he273       saturation vapour pressure at t=273k                      
    !     hkmelt      coefficient for melting of ice in precipitation           
    !     hmrcu       cloud water mixing ratio at which conversion becomes      
    !     efficient in convective cloud                             
    !     hmrst       cloud water mixing ratio at which conversion becomes      
    !     efficient in stratiform cloud                             
    !     hp0         reference pressure (=100 kpa)                             
    !     htaucu      characteristic time (not used at present)
    !     ht273       = 273 K                                                   
    !     hu00        threshold relative humidity for stratiform                
    !     condensation                                              
    !     hvsnow      terminal velocity of ice/snow precipitation               
    !     (not used at present)
    !     hvterm      terminal velocity of precipitation                        
    !     (not used at present)
    !     nlev        number of model levels. the level nearest to              
    !     to the ground                                             
    !-----------------------------------------------------------------------
    !     

    !     1.3)  Common block variables.
    !     


    !     2. Provisional updates due to 
    !     dynamics, vertical diffusion and radiation 

    !     
    !     
    do jk=1,nlev
       do jl=kstart,kstop
          !     
          wtp(jl,jk)=t(jl,jk) +dtime*dtdtin(jl,jk)
          wqp(jl,jk)=q(jl,jk) +dtime*dqdtin(jl,jk)
          wcwp(jl,jk)=cw(jl,jk) +dtime*dcwdin(jl,jk)
          dtdt(jl,jk)=( wtp(jl,jk) -t(jl,jk) )/dtime 
          dqdt(jl,jk)=( wqp(jl,jk) -q(jl,jk) )/dtime 
          dcwdt(jl,jk)=( wcwp(jl,jk) -cw(jl,jk) )/dtime
          !     
       enddo
    enddo
    !     
    do jl=kstart,kstop
       !     
       wpsp(jl)=ps(jl) +dtime*dpsdin(jl) 
       !     
    enddo
    !     


    !     3.  Determine label arrays defining convective layers,
    !     and subsequently convective heating and moistening 
    !     (subroutine 'condcv')

    !     
    call condcv( kcall,nhor,nlev,kstart,kstop, &
         dtime,conacc,dtheta,ahyb,bhyb,         &
         t,q,cw,pf,dpf,                         &
         klabcv,klabc,klab,                     &
         wpsp,wldcp,wtc,wtcdif,wqc,wqcdif,      &
         wqp,wtp,wcwp,wqsat,wdqsdt,wrhu,        &
         wqmax,wqliq,wqac,wbfrac,wdift,wdifq,   &
         wscdif,wdifs,                          &
         wbu,wdpcld,wdptot,wgpot,               &
         cucov,totcov )
    !     


    !     4.  Determine preliminary cloud cover needed primarily in the 
    !     precipitation release routine 'prevap'. 
    !     Currently, the radiation scheme applies during the first time
    !     step a cloud cover field equal to zero.
    !     As a refinement, cloud cover could be initialized to a 
    !     more realistic value to be applied at the first time step.
    !     In that case the parameter 'kcall' should be changed from 1 to 0.

    !     
    call cloudcv( nhor,nlev,kstart,kstop,kstep,  &
         dtime,conacc,dtheta,                     &
         ahyb,bhyb,ps,wrk1,                       &
         t,q,cw,totcov,cucov,pf,dpf,              &
         klabcv,klab,                             &
         wpsp,wldcp,wtc,wqc,                      &
         wqp,wtp,wcwp,wqsat,wdqsdt,wrhuc )
    !     
    if( kcall/=0 )then
       !     


       !     5. Determine stratiform condensation and moistening,
       !     and compute final tendencies due to condensation processes,
       !     by calling subroutine 'condst' 

       !     
       call condst( 1,nhor,nlev,kstart,kstop,   &
            dtime,conacc,                        &
            ahyb,bhyb,                           &
            t,q,cw,pf,dpf,                       &
            klabcv,klab,                         &
            wpsp,wldcp,                          &
            wqp,wtp,wcwp,wqsat,wdqsdt,wrhuc )
       !     


       !     6. Determine precipitation and associated evaporation
       !     by calling subroutine 'prevap'
       !     The computations are based on forecast variables
       !     updated due to condensation processes during the
       !     current time step. 

       !     
       call prevap( nhor,nlev,kstart,kstop,         &
            dtime,conacc,                            &
            ahyb,bhyb,dpf,                           &
            cucov,totcov,draindt,dsnowdt,            &
            preta,                                   &
            klabcv,                                  &
            prcpcu,prcpst,cusnow,stsnow,wpcoef,wmr,  &
            wdcusn,wdstsn,wrk1,wrk2,                 &
            wpsp,wldcp,wqp,wtp,wcwp,wqsat,wdqsdt )
       !     


       !     7. compute final tendencies of temperature, specific humidity and 
       !     cloud water due to condensation and precipitation processes

       !     
       do jk=1,nlev
          do jl=kstart,kstop
             dtdt(jl,jk)=(wtp(jl,jk) -t(jl,jk))/dtime  -dtdtin(jl,jk) 
             dqdt(jl,jk)=(wqp(jl,jk) -q(jl,jk))/dtime  -dqdtin(jl,jk) 
             dcwdt(jl,jk)=(wcwp(jl,jk) -cw(jl,jk))/dtime -dcwdin(jl,jk)  
          enddo
       enddo
       !     


       !     8. Diagnostic calculations  

       !     

       !     8.1  Initialize arrays and local constants.
       !     -------------------------------------------------------------
       !     'wrk1'  : work array
       !     'wrk2'  : work array
       !     'cov2d' : horizontal total cloud cover field (dimensionless)
       !     'cwpath': vericallly integrated cloud water (kg water/m2)
       !     -------------------------------------------------------------
       !     _realkind
       do jl = kstart,kstop                                               
          wrk1(jl) = 0.0_realkind                                                          
          wrk2(jl) = 0.0_realkind                                                          
          cov2d(jl)= 1.0_realkind 
          cwpath(jl)=0.0_realkind     
       enddo
       !     

       !     8.2 compute total horizontal cloud cover in two steps
       !     compute  vertically integrated cloud water.
       !     ---------------------------------------------------------
       do jk=1,nlev
          do jl=kstart,kstop
             if(jk > 1) then
                cov2d (jl) = cov2d(jl)*(1.0_realkind - max(totcov(jl,jk-1),&
                     totcov(jl,jk))) /(1.0_realkind- min(totcov(jl,jk-1),0.99_realkind))
             endif

             cwpath(jl)=cwpath(jl) +cw(jl,jk)*dpf(jl,jk)/gravit 
          enddo
       enddo

       !     final horizontal field 'cov2d' for total cloud cover

       do jl=kstart,kstop
          cov2d(jl) = 1.0_realkind - cov2d(jl)
       enddo

    endif
    return                                                                  
  end subroutine conds

  subroutine cloudcv(nhor,nlev,kstart,kstop,kstep, &
       dtime,conacc,dtheta,           &
       ahyb,bhyb,ps,wrk1,             &
       t,q,cw,totcov,cucov,pf,dpf,    &
       klabcv,klab,                   &
       wpsp,wldcp,wtc,wqc,            &
       wqp,wtp,wcwp,wqsat,wdqsdt,wrhuc )
    !
    !     --------------------------------------------------------------------------
    !

    !     Purpose:      Compute cloud cover in every model layer taking into account 
    !     --------      both stratiform and convective condensation processes.
    !                    
    !                   To counteract unrealistic 'on/off' switches a relaxation  
    !                   of both convective and stratiform cloud cover is carried 
    !                   out towards 'equilibrium conditions' determined from 
    !                   specified statistical frequency distributions of total
    !                   specific humidity ( vapour plus cloud condensate ).
    !  
    !                   The relaxation of cloud cover X is of the following type: 
    !                   ( Xeq is equilibrium value for the given atmospheric 
    !                     state and 'zlifex' is the time scale for adjustment )
    !                   
    !                     X = X +(dt/zlifex)*(Xeq -X) 
    !                         
    !     --------------------------------------------------------------------------
    !     Interface  :  subroutine 'condst' is called from subroutine 'conds'
    !     ------------
    !
    !     Input      :  ( from argument list or common blocks )
    !     ------------
    !       t            temperature at old time step 
    !       q            specific humidity at old time step           
    !       cw           cloud water content at old time step.
    !       dpf          delta-p around 'full' level
    !       pf           pressure at 'full' levels
    !       wpsp         provisional value of surface pressure
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !       klab         identifies potential convective 'entities', that is,
    !                    a cluster of consecutive layers where air parcel
    !                    lifting is buoyant.
    !       klabcv       is a convective label. Layers with the same
    !                    value of 'klab(jl,jk)' can only be accepted as
    !                    convective. Convective cloud layers are first
    !                    assigned a value of klabcv(jl,jk)=2, and a final
    !                    value of 3. Convective sub-cloud layers are
    !                    assigned a value of 1. New lifting levels
    !                    however, are assigned a value of 4.
    !                    convective 'entity'.
    !
    !     Output     :  
    !     ------------
    !     'cucov'        convective cloud cover (3-D) field 
    !     'totcov'       total cloud cover (3-D) field 
    !
    !                                    
    !     Work arrays:            
    !     ------------
    !     'wrk1' 
    !
    !     Subroutines:  none
    !     ------------     
    !
    !     Author:       Bent H. Sass, (DMI), based on the research during the
    !     -------                     the HIRLAM I-III projects.
    !

    !     1. Declarations. 

    !
    implicit none                                                           


    integer nhor,nlev,kstart,kstop,kstep, jk,jl,jkm1,jkp1,jintt
    integer klabcv(nhor,nlev),klab(nhor,nlev)
    real(kind=realkind) dtime,conacc,dtheta 
    real(kind=realkind) ztceff,zweight
    real(kind=realkind) zlimhu,zlimhuc
    real(kind=realkind) zwqp,zsec,zseva,zcpi,                     &
         zlifec,zlifes,zadc,zads,                   &
         zpf,zpdps,zqctot,zfr,zqmax,zqmin,          &
         zs1,zs2,zs3,zs4,zs5,zs6,zs7,               &
         zeqcov,zeqcos,zeqcot,                      &
         zss1,zss2,zss3,zss4,zss5,zss6,zss7,zcwlim 
    real(kind=realkind) zmepsi,zresat,zrelf0,zrelf1,zrelf2,&
         zdrelf,zrelfh,zdist,                &
         zak,zbk,zp0  
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1)
    real(kind=realkind) wrhuc(nhor,nlev) 
    real(kind=realkind) ps(nhor),wrk1(nhor),wpsp(nhor), wqsat(nhor),wdqsdt(nhor)
    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),                  &
         cw(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev),  &
         wtc(nhor,nlev),wqc(nhor,nlev),               &
         wldcp(nhor,nlev),wqp(nhor,nlev),             &
         wtp(nhor,nlev),wcwp(nhor,nlev),              &
         totcov(nhor,nlev),cucov(nhor,nlev)            

    !     2. Initialize various constants and arrays including labels

    !
    !     2.1 local constants   
    !     -------------------
    !     'zlifec'   is a relaxation time (s) for adjusting
    !               cloud cover in convective conditions
    !               towards current 'equilibrium conditions'.
    !
    !     'zlifes'   is a relaxation time (s) for adjusting
    !               cloud cover in stratiform conditions
    !               towards current 'equilibrium conditions'.
    !
    !     'zrelf1'  relative humidity threshold
    !               for stratiform condensation
    !
    !     'zrelf2'  coefficient connected to threshold
    !               for stratiform condensation
    !
    !     'zcwlim'  minimum 'cutoff' cloud condensate for clouds 
    !               to exist. 
    !    
    !     'zlimhuc' minimum relative humidity increment for security
    !
    !     'ztceff'  effective convective cloud parcel temperature
    !               (initialized to dummy value =300 K)
    !
    !     'zweight' weight factor to determine effective cloud parcel
    !               temperature
    !
    !     'zsec'    is a security constant.
    !     ------------------------------------------------------------
    !
    zlifec=900.0_realkind
    zlifes=900.0_realkind 
    zadc=min( 0.5_realkind*dtime/zlifec, 1.0_realkind )
    zads=min( 0.5_realkind*dtime/zlifes, 1.0_realkind )
    if(real(kstep,realkind)*dtime<=7200._realkind)then
       zads=min(1.0_realkind,(7200._realkind-real(kstep,realkind)*dtime)/&
            7200.0_realkind+  real(kstep,realkind)*dtime/7200._realkind*zads)
       zadc=min(1.0_realkind,(7200._realkind-real(kstep,realkind)*dtime)/&
            7200.0_realkind+  real(kstep,realkind)*dtime/7200._realkind*zadc)
    endif
    zcwlim=1.E-7_realkind
    zlimhuc=0.005_realkind
    zsec=0.001_realkind
    zcpi=1.0_realkind/cpair 
    zdist=sqrt( rearth*pi*dtheta/180.0_realkind )
    zrelf0=0.0200_realkind
    zrelf1=0.3000_realkind
    zrelf2=0.0030_realkind
    zdrelf=0.0300_realkind
    zrelfh=zrelf1*(1.0_realkind -exp( -zrelf2*zdist ))
    zp0=sip0
    zmepsi=1.0_realkind-epsilo
    ztceff=300.0_realkind
    zweight=0.50_realkind
    !


    !     3.  Compute convective cloud cover 'cucov'
    !         convective conditions are fulfilled ( 'klabcv =3').
    !
    !         The computations are based on the variables:
    !         'wqp':   preliminary grid box value of specific humidity.    
    !         'wcwp':  preliminary grid box value of specific 
    !                  cloud water.
    !         'wqsat': the saturation specific humidity at the
    !                  convective cloud temperature. 
    !         'zqctot': total specific humidity, that is, water             
    !                   vapour plus cloud water. 
    !         'wrhuc': the threshold value of relative fraction for
    !                  total scpecific humidity. 
    !         
    !        3.0  Compute threshold relative fraction  'wrhuc'
    !             of total specific humidity involved in cloud cover and
    !             liquid water computations.

    !
    do jl=kstart,kstop
       cucov(jl,1)=0.0_realkind
       totcov(jl,1)=0.0_realkind 
       wrhuc(jl,1)=zlimhuc
    enddo
    !
    do jk=2,nlev
       jkm1=jk-1
       jkp1=jk+1
       zak=0.5_realkind*( ahyb(jk) +ahyb(jkp1) )
       zbk=0.5_realkind*( bhyb(jk) +bhyb(jkp1) )
       zpdps=(zak +zbk*zp0)/zp0
       !
       do jl=kstart,kstop
          !       ------------------
          !        
          zpf=0.50_realkind*( ahyb(jk) +ahyb(jkp1) +  (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )
          !        
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk))
          wqsat(jl)=epsilo*esat(wtp(jl,jk))*zresat
          zqctot=wqp(jl,jk) +wcwp(jl,jk)
          zqctot=max( zqctot,zcwlim )
          zlimhu=1.0_realkind + (wqsat(jl) -2.0_realkind*wqp(jl,jk) )/zqctot
          wrhuc(jl,jk)=zrelfh*   (zrelf0 +zdrelf*(1.0_realkind -zpdps*zpdps*zpdps))/ &
               (zrelf0 +zdrelf)
          wrhuc(jl,jk)=max( min( wrhuc(jl,jk),zlimhu ), zlimhuc )
          !

          !        3.1 Compute convective 'equilibrium' cloud cover 'zeqcov'
          !            (first initialization)
          !
          zeqcov=0.0_realkind         
          !
          !        -----------------------------------------------------------

          !         Ajust convective cloud by relaxation towards 'equilibrium'
          !        -----------------------------------------------------------
          !
          if( klabcv(jl,jk)==3 ) then

             !      
             ztceff=zweight*wtc(jl,jk) +(1.0_realkind -zweight)*wtp(jl,jk)
             zresat=1.0_realkind/pecor(zpf,ztceff)
             wqsat(jl)=epsilo*esat(ztceff)*zresat
             !
             !          ----------------------------------------------------------

             !          3.2 First convective case, that is,
             !              zqctot <= wqsat 
             !          -----------------------------------------------------------
             if( zqctot<=wqsat(jl) ) then 
                zfr=max( wqsat(jl) -wqp(jl,jk), 0.0_realkind )/(max( wcwp(jl,jk), zcwlim))
                zeqcov=1.0_realkind/( 1.0_realkind + sqrt( zfr ) )
             endif
             !          ----------------------------------------------------------

             !          3.3 Second convective case, that is,
             !              zqctot*(1-wrhuc) <= wqsat < zqctot.
             !          ----------------------------------------------------------
             if( (zqctot*(1.0_realkind-wrhuc(jl,jk))<=wqsat(jl)) .and.(wqsat(jl)<zqctot) ) then 
                zeqcov= 0.50_realkind + ( zqctot  -wqsat(jl) )/  &
                     ( 2.0_realkind*wrhuc(jl,jk)*zqctot +zcwlim )
             endif
             !          ----------------------------------------------------------

             !          3.4 Special convective case, that is,
             !              wqsat <= zqctot*(1-wrhuc)
             !          ----------------------------------------------------------
             if( wqsat(jl)<=zqctot*(1.0_realkind -wrhuc(jl,jk)) ) then
                zeqcov=1.0_realkind
             endif
             !          ----------------------------------------------------------
             !          Ajust convective cloud by relaxation towards 'equilibrium'
             !          ----------------------------------------------------------
             !
             cucov(jl,jk)=cucov(jl,jk) +zadc*( zeqcov -cucov(jl,jk) )
             !          ----------------------------------------------------------
             !
          else  
             cucov(jl,jk)=cucov(jl,jk)*(1.0_realkind -zadc)
          endif

          !      

          !        4. Treat stratiform cloud cover
          !           equilibrium value 'zeqcos' to be determined 
          !        ------------------------------------------------------------       
          !
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk))
          wqsat(jl)=epsilo*esat(wtp(jl,jk))*zresat 
          zqmax=zqctot*(1.0_realkind +wrhuc(jl,jk))
          zqmin=zqctot*(1.0_realkind -wrhuc(jl,jk))
          !
          !    

          !        4.1 First case ( wqsat>=zqmax )
          !        -----------------------------------------------------------
          zeqcos=0.0_realkind
          !

          !        4.2 Second case ( intermediate situation )
          !        -----------------------------------------------------------
          if( (zqmin<=wqsat(jl)).and. (wqsat(jl)<zqmax) ) then  
             zeqcos=( zqmax -wqsat(jl))/  (2.0_realkind*zqctot*wrhuc(jl,jk))
          endif
          !        -----------------------------------------------------------

          !        5.0 Ajust total cloud cover towards 'equilibrium'
          !            being maximum of 'zeqcov' and 'zeqcos'
          !        -----------------------------------------------------------
          !
          zeqcot=totcov(jl,jk)
          !
          if( klabcv(jl,jk)==3 ) then 
             zeqcot=zeqcov
             if( klabcv(jl,jkm1)/=3 ) zeqcot=max( zeqcov, zeqcos )
          endif
          !
          if( klabcv(jl,jk)/=3 ) zeqcot=zeqcos
          !
          totcov(jl,jk)=totcov(jl,jk) +zads*(zeqcot -totcov(jl,jk))
          !
          !        -----------------------------------------------------------

          !        5.1 saturated case ( wqsat<zqmin )
          !        -------------------------------------
          if( wqsat(jl)<zqmin ) then
             totcov(jl,jk)=1.0_realkind
          endif
          !

       enddo
    enddo

    !


    !     6. Final cloud cover check. It is set to zero if
    !        if the minimum value 'zcwlim' of cloud condensate
    !        has been reached. This is consistent with the 
    !        adjustment in stratiform condensation which 
    !        completely evaporates this small amount of 
    !        cloud condensate to enable zero-values.
    !        For security, totcov >= cucov is current imposed
    !        (although this should not be needed)
    !     --------------------------------------------------------------
    !
    do jk=1,nlev
       do jl=kstart,kstop
          if( wcwp(jl,jk)<zcwlim ) then
             cucov(jl,jk)=0.00_realkind
             totcov(jl,jk)=0.00_realkind
          endif

          totcov(jl,jk)=max( cucov(jl,jk),totcov(jl,jk) )

       enddo
    enddo

    return                                                                  
  end subroutine cloudcv

  subroutine   condcv( kcall,nhor,nlev,kstart,kstop,  &
       dtime,conacc,dtheta,ahyb,bhyb,          &
       t,q,cw,pf,dpf,                          &
       klabcv,klabc,klab,                      &
       wpsp,wldcp,wtc,wtcdif,wqc,wqcdif,       &
       wqp,wtp,wcwp,wqsat,wdqsdt,wrhu,         &
       wqmax,wqliq,wqac,wbfrac,wdift,wdifq,    &
       wscdif,wdifs,                           &
       wbu,wdpcld,wdptot,wgpot,                &
       cucov,totcov )
    !
    !     -------------------------------------------------------------------------

    !     --------   
    !     Purpose       Compute changes of temperature, specific humidity and
    !     --------      specific cloud water due to convective condensation
    !                   processes. An advanced  'KUO'- type of closure is used. 
    !                   Several convective 'entities' can be treated.
    !                   First, convective labels ('klab' ,'klabc' and 'klabcv')
    !                   are determined.
    !
    !     -------------------------------------------------------------------------
    !     Interface  :  subroutine 'condcv' is called from subroutine 'conds'
    !     ------------
    !
    !     Input      :  ( from argument list or common blocks )
    !     ------------
    !       cw           cloud water content at (t-dt)
    !       dpf          delta-p around 'full' level
    !       pf           pressure at 'full' levels
    !       wpsp         provisional value of surface pressure
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !       klab         identifies potential convective 'entities', that is,
    !                    a cluster of consecutive layers where air parcel
    !                    lifting is buoyant.
    !       klabc        is a label assigned a value of 1 in all convective layers
    !                    of the convective entities that contain at least one 
    !                    convective cloud layer, otherwise a label value of 0. 
    !       klabcv       is a convective label. Layers with the same
    !                    value of 'klab(jl,jk)' can only be accepted as
    !                    convective. Convective cloud layers are first
    !                    assigned a value of klabcv(jl,jk)=2, and a final
    !                    value of 3. Convective sub-cloud layers are
    !                    assigned a value of 1. New lifting levels
    !                    however, are assigned a value of 4.
    !       wdift        weighting parameters for temperature change 
    !       wtcdif       in the convection scheme
    !       wdifq        weighting parameters for humidity change  
    !       wqcdif       in the convection scheme
    !       wdifs        weighting parameters for cloud condensate
    !       wscdif       change in the convection scheme
    !       wqac         vertically integrated moisture convergence
    !                    convective 'entity'.
    !       wdptot       total pressure thickness of convective 'entity'
    !       wdpcld       total pressure thickness of convective
    !                    'cloud entity'
    !       wbfrac       moisture partitioning parameter in the
    !                    convection scheme
    !       wbu          averarage buoyancy '(tc -t)/t)' in the
    !                    convective 'entity'
    !     Output     :  ( new output or modified input) 
    !     ------------
    !       wldcp(jl,jk)  (1/cpair)*latent heat as a function of
    !                     temperature by weighting in latent heat of
    !                     freezing times probability of ice crystals
    !       wqsat(jl,jk)  saturation mixing ratio with respect to
    !                     gridpoint temperature
    !       wdqsdt(jl,jk) derivative of saturation water vapour pressure
    !                     with respect to temperature
    !       prbice(1750) probability for ice crystals as a function of
    !                    temperature
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !       wtc          cloud parcel temperature
    !       wqc          cloud parcel specific humidity
    !
    !     Subroutines:  none
    !     ------------     
    !
    !     Author:       Bent H. Sass, (DMI), based on the research during the
    !     -------                     the HIRLAM I-III projects.
    !

    !     1. Declarations. 

    !
    implicit none                                                           


    integer nhor,nlev,kstart,kstop,knlevm,knlevp,kcall,jk,jl,jkm1,jkp1,jkp2,kklab,kkla,jintt
    integer klabcv(nhor,nlev),klabc(nhor,nlev),klab(nhor,nlev)
    real(kind=realkind) dtime,conacc,dtheta  
    real(kind=realkind) zdthecr
    real(kind=realkind) zpf,zpfp1,zsec,zpcr,ztvir,ztvirp,ztvirc,zqsat,     &
         zceva,zdift,zdifq,zdifq1,zt,zdifp,zqcd,zqcd1,zdgpot,&
         zcev,zbucr,zdpf,zcpi,zc,zc1,zc2,zc3,zc4,zc5,zc6,    &  
         zc7,ztunc,zcc1,zcc2,zcc3,zcc4,zcc5,zcc6,zcc7,       &
         zfrc1,zpardt,zmepsi,zresat  
    real(kind=realkind) zkaccr,zqpar,ztpar,zdelta,ztmp  
    real(kind=realkind) zspar,zqcmax 
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1)
    real(kind=realkind) wqmax(nhor),wpsp(nhor),wqsat(nhor),wdqsdt(nhor)
    real(kind=realkind) wgpot(nhor,nlev),t(nhor,nlev),q(nhor,nlev),        &
         cw(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev),         &
         wldcp(nhor,nlev),wtc(nhor,nlev),wtcdif(nhor,nlev),  &
         wqc(nhor,nlev),wqcdif(nhor,nlev),wqp(nhor,nlev),    &
         wtp(nhor,nlev),wcwp(nhor,nlev),                     &
         wrhu(nhor,nlev),wqliq(nhor,nlev),wdift(nhor,nlev),  &
         wqac(nhor,nlev),wbfrac(nhor,nlev),wdifq(nhor,nlev), &
         wbu(nhor,nlev),wdpcld(nhor,nlev),wdptot(nhor,nlev), &
         totcov(nhor,nlev),cucov(nhor,nlev)            
    real(kind=realkind) wscdif(nhor,nlev),wdifs(nhor,nlev)
    !


    !     2. Initialize various constants and arrays including labels

    !
    !     2.1 local constants   
    !     -------------------
    !     'zceva'   is a coefficient for cloud water evaporation
    !               reflecting the invese life time of convective 
    !               clouds in the absence of a convective moisture 
    !               source.
    !     'zpardt'  air parcel initial excess temperature increment.
    !     'zbucr'  is a threshold buoyancy that must be reached 
    !               as a mean vertical buoyancy in the convective
    !               'entities', in order for convective activity 
    !               to be fully developed. 
    !     'zkaccr'  tuning constant to control maximum convective
    !               cloud depth as a function of moisture                
    !               availability. 
    !     'zdthecr' scaling threshold (tunable) related to resolution.
    !               (as default active for hor. resolution below 0.125 deg.)
    !     'zqpar'   moistening parameter not available in the ECMWF
    !               type of KUO convection scheme.
    !     'ztpar'   heating parameter not available in the ECMWF
    !               type of KUO convection scheme.
    !     'zspar'   cloud condensate parameter similar to 'zqpar'
    !               for specific humidity.
    !     'zqcmax'  maximum value for cloud condensate in cloud          
    !               ascent (tuning parameter)
    !     'zpcr'    is a threshold value of the pressure thickness
    !               of convective 'entities' to exceed.
    !     'zcpi'    inverse of specific heat capacity
    !     'zsec'    is a security constant.
    !     'ztunc'   tuning constant to control numerical stability
    !               in computing evaporation of cloud water.
    !     'zfrc1'   absolute minimim specific cloud water 
    !               allowed in the atmosphere.

    !


    !
    zceva=0.18_realkind
    zpardt=1.00_realkind
    ztunc=0.50_realkind
    zbucr=0.0030_realkind
    zkaccr=2.2E8_realkind/dtime
    zdthecr=0.125_realkind
    zqpar=0.0015_realkind
    ztpar=0.0_realkind
    zspar=0.0_realkind
    zqcmax=0.0030_realkind
    zsec=0.00001_realkind 
    zfrc1=1.0E-7_realkind  
    zpcr=10000.0_realkind
    zcpi=1.0_realkind/cpair 
    zmepsi=1.0_realkind -epsilo
    knlevm=nlev-1  
    knlevp=nlev+1
    !     ---------------------------------------------------------
    !     2.2 arrays 
    !     ---------------------------------------------------------
    do jk=1,nlev
       do jl=kstart,kstop
          !        
          klabcv(jl,jk)=0
          klabc(jl,jk)=0
          klab(jl,jk)=0
          wqac(jl,jk)=0.0_realkind
          wdptot(jl,jk)=0.0_realkind
          wdpcld(jl,jk)=0.0_realkind
          wbfrac(jl,jk)=0.0_realkind
          wbu(jl,jk)=0.0_realkind
          wgpot(jl,jk)=0.0_realkind
          wtc(jl,jk)=0.0_realkind
          wtcdif(jl,jk)=0.0_realkind
          wdift(jl,jk)=0.0_realkind
          wdifq(jl,jk)=0.0_realkind
          wscdif(jl,jk)=0.0_realkind
          wdifs(jl,jk)=0.0_realkind
          wqc(jl,jk)=0.0_realkind
          wqcdif(jl,jk)=0.0_realkind
          wqliq(jl,jk)=0.0_realkind
          !
          !         latent heat as a function of temperature.
          !         -----------------------------------------
          jintt = max(1, 1 + nint((tmelt - wtp(jl,jk))/dttabl))
          jintt = min(jintt, 1750)
          wldcp(jl,jk) = zcpi * (latvap + hdl * prbice(jintt))

       enddo
    enddo

    do jl=kstart,kstop
       wqsat(jl)=0.0_realkind
       wdqsdt(jl)=0.0_realkind
       wqmax(jl)=0.0_realkind
    enddo
    !


    !     3.  First step towards determining labels and arrays used in the
    !         convective parameterization. This involves air parcel lifting.
    !
    !     'klab(jl,jk)'   identifies potential convective 'entities', that 
    !                     is, a cluster of consecutive layers where air 
    !                     parcel lifting is buoyant.
    !     'klabcv(jl,jk)' is a convective label. Layers with the same
    !                     value of 'klab(jl,jk)' can only be accepted as
    !                     convective. Convective cloud layers are first
    !                     assigned a value of klabcv(jl,jk)=2, and a final
    !                     value of 3. Convective sub-cloud layers are
    !                     assigned a value of 1. However, new lifting
    !                     levels are assigned a value of 4.
    !     'klabc(jl,jk)'  is a label assigned a value of 1 in all convec-
    !                     tive layers of the convective entities that
    !                     contain at least one convective cloud layer,
    !                     otherwise a label value of 0.

    !
    do jl=kstart,kstop
       !
       wtc(jl,nlev)=wtp(jl,nlev) +zpardt 
       wqc(jl,nlev)=wqp(jl,nlev)
       wqmax(jl)=wqp(jl,nlev)
       zdifp=wpsp(jl)*(1.0_realkind -bhyb(nlev))  -ahyb(nlev)
       zpf=0.5_realkind*( ahyb(nlev) +ahyb(knlevp) +  (bhyb(nlev) +bhyb(knlevp))*wpsp(jl) )
       ztvir=wtp(jl,nlev)*(1.0_realkind +tvirtc*wqp(jl,nlev) -wcwp(jl,nlev))
       wgpot(jl,nlev)=log(wpsp(jl)/(zpf +zsec))
       wgpot(jl,nlev)=rair*ztvir*wgpot(jl,nlev)
       !
       klabcv(jl,nlev)=4
       klab(jl,nlev)=1
       !
       wdptot(jl,1)=2.0_realkind*(wpsp(jl) -zpf)
       wqac(jl,1)=zdifp*wqp(jl,nlev) -dpf(jl,nlev)*q(jl,nlev)
       !
    enddo
    !
    do jk=knlevm,1,-1
       !     --------------------
       jkp1=jk+1
       jkp2=jk+2
       do jl=kstart,kstop
          zpf=0.5_realkind*( ahyb(jk) +ahyb(jkp1) +   (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )
          zpfp1=0.5_realkind*( ahyb(jkp1) +ahyb(jkp2) + (bhyb(jkp1) +bhyb(jkp2))*wpsp(jl) )
          zdifp=ahyb(jkp1) -ahyb(jk) + ( bhyb(jkp1) -bhyb(jk))*wpsp(jl)
          !
          ztvir=wtp(jl,jk)*(1.0_realkind +tvirtc*wqp(jl,jk) -wcwp(jl,jk))
          ztvirp=wtp(jl,jkp1)* (1.0_realkind +tvirtc*wqp(jl,jkp1) -wcwp(jl,jkp1))
          wgpot(jl,jk)=wgpot(jl,jkp1) + rair*log(zpfp1/zpf)* 0.5_realkind * (ztvirp + ztvir)
          !
          jintt = max(1, 1 + nint((tmelt - wtc(jl,jkp1))/dttabl))
          jintt = min(jintt, 1750)
          wldcp(jl,jk) =(latvap + latice* prbice(jintt))/cpair
          !
          zdgpot =wgpot(jl,jkp1) - wgpot(jl,jk)
          !
          wqliq(jl,jkp1)=wqmax(jl)-wqc(jl,jkp1)
          !
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk)) 
          zqsat=epsilo*esat(wtp(jl,jk))*zresat 
          wrhu(jl,jk)=wqp(jl,jk)/zqsat
          !
          !         do dry adiabatic lifting
          !
          wtc(jl,jk)=wtc(jl,jkp1) + wtc(jl,jkp1)*(1.0_realkind +tvirtc*wqc(jl,jkp1)&
               -wqliq(jl,jkp1))* zdgpot /(ztvirp*cpair )
          !
          !         release latent heat during parcel lifting
          !
          zresat=1.0_realkind/pecor(zpf,wtc(jl,jk))
          wqsat(jl)=epsilo*esat(wtc(jl,jk))*zresat
          wqc(jl,jk)=wqc(jl,jkp1)
          zdifq=wqc(jl,jk)-wqsat(jl)
          zt=wtc(jl,jk)
          wdqsdt(jl)=epsilo*zpf*desdt(zt)* zresat*zresat 
          zqcd =zdifq/(1.0_realkind +wldcp(jl,jk)*wdqsdt(jl))
          zqcd =max(zqcd ,0.0_realkind)
          !
          !         new cloud parcel humidity and temperature
          !
          wqc(jl,jk)=wqc(jl,jk) -zqcd
          wtc(jl,jk)=wtc(jl,jk) +zqcd * wldcp(jl,jk)
          !
          !         adjust a second time towards cloud parcel saturation.
          !
          if( zqcd>0.0_realkind ) then 
             zresat=1.0_realkind/pecor(zpf,wtc(jl,jk)) 
             wqsat(jl)=epsilo*esat(wtc(jl,jk))*zresat 
             zdifq1=wqc(jl,jk)-wqsat(jl)
             zqcd1=zdifq1/(1.0_realkind +wldcp(jl,jk)*wdqsdt(jl))
             wqc(jl,jk)=wqc(jl,jk) -zqcd1
             wtc(jl,jk)=wtc(jl,jk) +zqcd1 * wldcp(jl,jk)
          endif
          !
          !         buoyancy check 
          !
          wqliq(jl,jk)=wqmax(jl) -wqc(jl,jk)
          wqliq(jl,jk)=min( wqliq(jl,jk), zqcmax )
          ztvirc=wtc(jl,jk)*(1.0_realkind +tvirtc*wqc(jl,jk) -wqliq(jl,jk))
          zdift=ztvirc -ztvir
          !
          !         -------------------------------------------------------

          !         Check if the convection can continue, i.e. the cloud
          !         parcel should be positively buoyant.
          !         It is also assumed that there is a maximum depth of a
          !         convective entity depending on the available moisture
          !         convergence (relevant for small moisture supplies only)
          !         -------------------------------------------------------
          !
          kklab=klab(jl,jkp1)
          !
          ztmp=zkaccr*wqac(jl,kklab) -wdptot(jl,kklab) -zdifp
          !
          if((zdift>0.0_realkind).and.(ztmp>0.0_realkind)) then
             !         -------------------------------------------------------
             wqac(jl,kklab)=wqac(jl,kklab) +  zdifp*wqp(jl,jk) -dpf(jl,jk)*q(jl,jk)
             wbu(jl,kklab)=wbu(jl,kklab) +zdifp*zdift/ztvir
             wdptot(jl,kklab)=wdptot(jl,kklab) +zdifp
             !
             klab(jl,jk)=klab(jl,jkp1)
             klabcv(jl,jk)=1
             !
             if( zdifq>0.0_realkind ) then
                !           ----------------------
                klabcv(jl,jk)=2
                wdift(jl,jk)=ztvirc -ztvir +ztpar 
                wdifq(jl,jk)=wqc(jl,jk) -wqp(jl,jk) +zqpar 
                wdifq(jl,jk)=max( wdifq(jl,jk),0.0_realkind )
                wdifs(jl,jk)=wqliq(jl,jk) +zspar 
                wdifs(jl,jk)=max( wdifs(jl,jk),0.0_realkind )
                wtcdif(jl,kklab)=wtcdif(jl,kklab) +zdifp*wdift(jl,jk)    
                wqcdif(jl,kklab)=wqcdif(jl,kklab) +zdifp*wdifq(jl,jk)
                wscdif(jl,kklab)=wscdif(jl,kklab) +zdifp*wdifs(jl,jk)
                wbfrac(jl,kklab)=wbfrac(jl,kklab) +zdifp*wrhu(jl,jk)
                wdpcld(jl,kklab)=wdpcld(jl,kklab) +zdifp
             endif
             !         ------------------------------------------------------------
          else
             !         ------------------------------------------------------------

             !         prepare for finding new convective entities
             !
             wtc(jl,jk)=wtp(jl,jk) +zpardt 
             wqc(jl,jk)=wqp(jl,jk)
             wqmax(jl)=wqp(jl,jk)
             klabcv(jl,jk)=4
             kkla=klab(jl,jkp1)+1
             klab(jl,jk)=kkla
             wtcdif(jl,kkla)=0.0_realkind
             wqcdif(jl,kkla)=0.0_realkind
             wscdif(jl,kkla)=0.0_realkind
             wqac(jl,kkla)=zdifp*wqp(jl,jk) -dpf(jl,jk)*q(jl,jk)
             wbfrac(jl,kkla)=0.0_realkind
             wbu(jl,kkla)=0.0_realkind
             wdptot(jl,kkla)=zdifp
             wdpcld(jl,kkla)=0.0_realkind
             !         ----------------------------------------------------------------
          endif
          !         ----------------------------------------------------------------
       enddo
       !       -----
    enddo
    !     -----


    !     4.  Last steps towards determining labels and arrays used in the
    !         convective parameterization.
    !
    !     (*) The final value of a convective cloud label 'klabcv' is set to 3.
    !         The model layers are investigated 'from top to bottom'
    !         of the atmosphere searching first the top cloud layer 
    !         of a convective entity. This layer is assigned a value of 3
    !         provided that 
    !
    !            i) 'klabcv' has already been assigned a value of 2 
    !         
    !           ii) the vertical integral of moisture accession over 
    !               the entire convective entity is positive 
    !          
    !         If the top cloud layer convective label 'klabcv' is set to 3 
    !         the convective indicator 'klabc' is set to 1, which means 
    !         that the convective entity is finally classified as being active.
    ! 
    !     (*) Secondly it is assured that all layers of the convective entity 
    !         have the same value of 'klabc' as the uppermost one.
    !
    !     (*) Thirdly, the label 'klabcv' is adjusted to a value of 3 in other
    !         layers below the cloud top provided that the appropriate 
    !         conditions are met as mentioned above for the top layer. 

    !
    do jk=2,nlev
       jkm1=jk-1
       do jl=kstart,kstop
          !       -------------------------------------------
          kklab=klab(jl,jk)
          if( (klab(jl,jk)/=klab(jl,jkm1)).and. (klabcv(jl,jk)==2).and. &
               (wqac(jl,kklab)>0.0_realkind) ) then 
             !         -----------------------------------------
             klabcv(jl,jk)=3
             klabc(jl,jk)=1 
             !         -----------------------------------------
          endif
          !         -----
          !
          !         -----------------------------------------
          if( (klab(jl,jk)==klab(jl,jkm1)).and.  (klabc(jl,jkm1)==1) ) then
             klabc(jl,jk)=1
             !         -----
          endif
          !         -----
          !
          !         -----------------------------------------
          if( (klab(jl,jk)==klab(jl,jkm1)).and.(klabcv(jl,jk)==2).and. &
               (wqac(jl,kklab)>0.0_realkind) ) then 
             !         -----------------------------------------
             klabcv(jl,jk)=3
             !         -----
          endif
          !         -----
          !
       enddo
    enddo
    !
    !     -----------------------------------------------------------
    do jk=2,nlev
       jkm1=jk-1
       do jl=kstart,kstop
          if( (klab(jl,jk)/=klab(jl,jkm1)).and.(klabcv(jl,jk)==3) ) then
             !         -------------------------------------------------------

             !           do final computation of convective work arrays
             !
             !           i) buoyancy 'wbu' is averaged over all layers 
             !              including convective sub-cloud layers)
             !           ii) other parameters are to be averaged  
             !               over the total convective cloud depth. 
             !
             !            Final computation of moistening parameter 
             !            ('beta'-parameter in KUO-type of schemes). 
             !            It is noted that this parameter is difficult 
             !            to parameterize, i.e., there is uncertainty involved.
             !            Currently, the formulation according to
             !            (1 - 'vertically averaged relative humidity')**2
             !            is used.  In the code below this parameter is 
             !            called 'wbfrac'. 
             !            There are indications that small precipitation 
             !            amounts are forecast too often. Some improvement 
             !            in this respect may be achieved by using 
             !            ( 1 - 'vertically averaged relative humidity')
             !            instead, that is a power of one instead of a 
             !            power of two. It seems that an improved formulation  
             !            should allow for effective moistening in relation to 
             !            relatively shallow convective clouds.
             !           -----------------------------------------------------
             !
             kklab=klab(jl,jk)
             wbu(jl,kklab)=wbu(jl,kklab)/(wdptot(jl,kklab)+zsec)
             wqac(jl,kklab)=wqac(jl,kklab)/(wdpcld(jl,kklab)+zsec)
             wtcdif(jl,kklab)=wtcdif(jl,kklab)/(wdpcld(jl,kklab) +zsec)
             wqcdif(jl,kklab)=wqcdif(jl,kklab)/(wdpcld(jl,kklab) +zsec)
             wscdif(jl,kklab)=wscdif(jl,kklab)/(wdpcld(jl,kklab) +zsec)
             wbfrac(jl,kklab)=wbfrac(jl,kklab)/(wdpcld(jl,kklab) +zsec)
             wbfrac(jl,kklab)=(1.0_realkind -wbfrac(jl,kklab))*(1.0_realkind -wbfrac(jl,kklab))
             ! 
          endif
          !         -----------------------------------------------------
       enddo
       !       -----
    enddo
    !     -----
    !


    !     5.  Convective updates. 
    !         First  adjust values of specific humidity in 
    !         convective layers, that is, update to old time
    !         values must be done before redistributing moisture 
    !         due to convection in the convective 'entities'.

    !
    if( kcall==1 ) then 

       do jk=1,nlev
          jkp1=jk+1
          do jl=kstart,kstop
             !       ------------------
             !

             !         5.1  Define 'help' coefficients 
             !         -------------------------------------------------------
             !         'zdpf' is pressure thickness between 'half' levels at
             !         new time step
             !         Do only convective updates if 'kcall=1'
             !         --------------------------------------------------------
             !
             zdpf=ahyb(jkp1) -ahyb(jk) + wpsp(jl)*(bhyb(jkp1) -bhyb(jk))    
             zpf=0.5_realkind*( ahyb(jk) +ahyb(jkp1) + (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )
             zresat=1.0_realkind/pecor(zpf,t(jl,jk))
             wqsat(jl)=epsilo*esat(t(jl,jk))*zresat
             wdqsdt(jl)=epsilo*zpf*desdt(t(jl,jk))* zresat*zresat 
             !
             kklab=klab(jl,jk)
             !

             !         5.2  Adjust variables in convective entities
             !         --------------------------------------------------------
             !         For shallow convective entities there
             !         is a linear weighting (paramameter 'zdelta') between
             !         updates due to dynamics and vertical diffusion (as 
             !         an average over the convective entity) on one
             !         hand and updates due to convection on the other. 
             !         The idea is to describe a gradual transition (linear)
             !         between 'small scale vertical mixing', that is ,
             !         vertical diffusion dominating for shallow phenomena 
             !         ( small 'wdptot' )  and deep convective phenomena 
             !         ( large 'wdptot' ). Also it is assumed that a sufficient 
             !         buoyancy is needed to activate fully the deep convection. 
             !         This feature also tends to prevent that full convective
             !         activity switches 'on and off'  between alternating
             !         time steps. Finally a linear term depending on  
             !         horizontal resolution ('dtheta') has been introduced
             !         in the expression for 'zdelta'. The idea is that convec-
             !         tion parameterization should gradually 'switch off'
             !         as the model resolution increases, that is, the
             !         convection should not be active in the limit where 
             !         horizontal resolution is very high because the dynamics
             !         should then resolve all vertical motions responsible 
             !         for convective precipitation. This approach is 
             !         described further in the HIRLAM scientific documentation.
             !         ---------------------------------------------------------
             !
             zdelta=min( wbu(jl,kklab)/zbucr, 1.0_realkind )*&
                  min( wdptot(jl,kklab)/zpcr, 1.0_realkind )* min( dtheta/zdthecr, 1.0_realkind )
             !         ---------------------------------------------------------

             !         Adjustment towards old values (taking into account mass
             !         changes), before adding convective updates during the
             !         current time step.
             !
             if( klabc(jl,jk)==1 ) then
                zc1=(1.0_realkind-zdelta)*wqp(jl,jk)
                wqp(jl,jk)=zdelta*q(jl,jk)*dpf(jl,jk)/(zdpf+zsec) +zc1
             endif
             !
             !         --------------------------------------------------------

             !         5.3 add convective updates (only if 'kcall =1' )
             !         numerical control of cloud water evaporation 
             !         by parameters 'ztunc' and 'zcev'.
             !         ----------------------------------------------------------
             !
             if( klabcv(jl,jk)==3 ) then
                !         ----------------------------------------------------------
                zc1=wqsat(jl)
                zc2=wdqsdt(jl)
                zc3=zdelta*wqac(jl,kklab)*(1.0_realkind -wbfrac(jl,kklab))* &
                     wdifs(jl,jk)/(wscdif(jl,kklab) +zsec)
                zc4=zdelta*wqac(jl,kklab)*wbfrac(jl,kklab)*  wdifq(jl,jk)/(wqcdif(jl,kklab) +zsec)
                zcev=ztunc/(dtime*max(wqsat(jl) -q(jl,jk),zsec))
                zcev=min( zceva, zcev )
                zc5=wldcp(jl,jk)*zdelta*wqac(jl,kklab)*(1.0_realkind -&
                     wbfrac(jl,kklab))*wdift(jl,jk)/(wtcdif(jl,kklab) +zsec)
                zc6=zcev*dtime*wcwp(jl,jk)
                zc7=wldcp(jl,jk)*zc6
                !
                zcc1=( wtp(jl,jk) +zc5 -zc1*zc7 +zc2*zc7*t(jl,jk) )/(1.0_realkind +zc2*zc7)
                zcc2=zc7/(1.0_realkind +zc2*zc7)
                zcc3=( wqp(jl,jk) +zc4 +zc1*zc6 -zc2*zc6*t(jl,jk) )/(1.0_realkind +zc6)
                zcc4=zc2*zc6/(1.0_realkind +zc6) 
                zcc5=wcwp(jl,jk) +zc3 -zc1*zc6 +zc2*zc6*t(jl,jk) 
                zcc6=-zc2*zc6 
                ! 
                wtp(jl,jk)=(zcc1 +zcc2*zcc3)/(1.0_realkind -zcc2*zcc4)
                wqp(jl,jk)=(zcc3 +zcc1*zcc4)/(1.0_realkind -zcc2*zcc4)
                wcwp(jl,jk)=zcc5 +zcc6*wtp(jl,jk) +zc6*wqp(jl,jk)
                !
             endif
             !
             !         -------------------------------------------------------

             !         5.4  Allow cloud water to possibly evaporate completely
             !              Make associated correction of specific humidity
             !              and temperature 
             !         -------------------------------------------------------
             if( wcwp(jl,jk)<zfrc1 ) then 
                wqp(jl,jk)=wqp(jl,jk) +wcwp(jl,jk)
                wtp(jl,jk)=wtp(jl,jk) -wldcp(jl,jk)*wcwp(jl,jk)
                wcwp(jl,jk)=0.0_realkind 
             endif
          enddo
       enddo
    endif
    return                                                                  
  end subroutine condcv


  subroutine condst( kscall,nhor,nlev,kstart,kstop,&
       dtime,conacc,                  &
       ahyb,bhyb,                     &
       t,q,cw,pf,dpf,                 &
       klabcv,klab,                   &
       wpsp,wldcp,                    &
       wqp,wtp,wcwp,wqsat,wdqsdt,wrhuc )
    !
    !     ---------------------------------------------------------------
    !

    !
    !     Purpose: (*)  Compute changes of temperature, specific humidity
    !     --------      and specific cloud water due to stratiform 
    !                   condensation processes.
    !
    !              (*)  A necessary condition for condensation
    !                   is that the relative humidity exceeds a threshold 
    !                   value that currently depends on the the vertical 
    !                   coordinate. 
    !        
    !              (*)  For lower relative humidity the existing 
    !                   cloud water decays exponentially towards zero. 
    !                   within a suitable time scale.
    !
    !              (*)  If convection has been active stratiform 
    !                   condensation is not active unless a given layer
    !                   has become supersaturated. In the latter situa-
    !                   tion the stratiform scheme takes care of the 
    !                   supersaturation.
    !
    !     ----------------------------------------------------------------
    !     Interface  :  subroutine 'condst' is called 
    !     ------------  from subroutine 'conds'        
    !
    !         Nomenclature:
    !         -------------
    !
    !         'wqsat'  saturation specific humidity at preliminary 
    !                  input temperature.
    !         'wdqsdt' partial derivative of saturation vapour pressure with
    !                  respect to temperature.
    !         'wtp'    preliminary updated temperature due to all other
    !                  processes
    !                  than stratiform condensation
    !         'wqp'    preliminary updated specific humidity due to all
    !                  other processes
    !                  than stratiform condensation
    !         'wcwp'   preliminary updated temperature due to all other
    !                  processes than stratiform condensation
    !         'wrhuc'  threshold fraction of total specific humidity 'zqtot'
    !         'zqtot'  sum of preliminary cloud water and specific humidity
    !         'zs1'    intermediate variable
    !         'zs2'    intermediate variable
    !         'zs3'    intermediate variable
    !         'zs4'    intermediate variable
    !         'zss1'   coefficient in final update of temperature and
    !                  humidity due to evaporation of cloud water
    !         'zss2'   coefficient in final update of temperature and 
    !                  humidity due to evaporation of cloud water
    !         'zss3'   coefficient in final update of temperature and 
    !                  humidity due to evaporation of cloud water
    !         'zss4'   coefficient in final update of temperature and
    !                  humidity due to evaporation of cloud water
    !         'zss5'   coefficient in final update of cloud water
    !                  due to evaporation of cloud water.
    !         'zss6'   coefficient in final update of cloud water
    !                  due to evaporation of cloud water.
    !         'zss7'   coefficient in final update of cloud water
    !                  due to evaporation of cloud water.
    !
    !     Input      :  ( from argument list or common blocks )
    !     ------------
    !       t            temperature at old time step 
    !       q            specific humidity at old time step           
    !       cw           cloud water content at old time step.
    !       dpf          delta-p around 'full' level
    !       pf           pressure at 'full' levels
    !       dqdt         total provisional humidity tendency due to 
    !                    all processes including dynamics. 
    !       wpsp         provisional value of surface pressure
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !       klab         identifies potential convective 'entities',
    !                    that is a cluster of consecutive layers where
    !                    air parcel lifting is buoyant.
    !       klabcv       is a convective label. Layers with the same
    !                    value of 'klab(jl,jk)' can only be accepted as
    !                    convective. Convective cloud layers are first
    !                    assigned a value of klabcv(jl,jk)=2, and a final
    !                    value of 3. Convective sub-cloud layers are
    !                    assigned a value of 1. New lifting levels
    !                    however, are assigned a value of 4.
    !                    convective 'entity'.
    !
    !     Output     :   (new output or modified input) 
    !     ------------
    !       wldcp(jl,jk)  (1/cpair)*latent heat as a function of
    !                     temperature by weighting in latent heat of
    !                     freezing times probability of ice crystals
    !       wqsat(jl,jk)  saturation mixing ratio with respect to
    !                     gridpoint temperature
    !       wdqsdt(jl,jk) derivative of saturation water vapour pressure
    !                     with respect to temperature
    !       prbice(1750)  probability for ice crystals as a function of
    !                     temperature
    !       wqp           provisional value of specific humidity
    !       wtp           provisional value of temperature
    !       wcwp          provisional value of specific cloud water
    !
    !     Subroutines:  none
    !     ------------     
    !
    !     Author:       Bent H. Sass, (DMI), based on the research 
    !     -------                     during the HIRLAM I-III projects.
    !

    !     1. Declarations. 

    !
    implicit none                                                           


    integer kscall,nhor,nlev,kstart,kstop,jkp1,jk,jl,kit,jj,jintt
    integer klabcv(nhor,nlev),klab(nhor,nlev)
    real(kind=realkind) dtime,conacc 
    real(kind=realkind) zcwmin,zsec,zseva,zqmin,zqmax,zcpi,z1,&
         zsev,zstun,zdqliq,zqtot,zqceq,zpf,     &
         zads,zlife,zs1,zs2,zs3,zs4,            &
         zss1,zss2,zss3,zss4,zss5,zss6,zss7,    &
         zmepsi,zresat  
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1)
    real(kind=realkind) wpsp(nhor),wqsat(nhor),wdqsdt(nhor)
    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),                 &
         cw(nhor,nlev),pf(nhor,nlev),dpf(nhor,nlev), &
         wldcp(nhor,nlev),wqp(nhor,nlev),            &
         wtp(nhor,nlev),wcwp(nhor,nlev),             &
         wrhuc(nhor,nlev)



    !


    !     2. Initialize various constants and arrays including labels

    !
    !     2.1 local constants   
    !     -------------------
    !     'zseva'   is a coefficient for cloud water evaporation
    !     'zlife'   adjustment time scale
    !     'zads'    adjustment factor
    !     'zcwmin'  minimum scecific humidity or cloud water content
    !     'zsec'    is a security constant.
    !     'zstun'   numerical tuning constant
    !               for evaporation of cloud water.
    !     'kit'     number of iterations in stratiform calculations.
    !     ----------------------------------------------------------
    !
    zseva=0.20_realkind
    zlife=900.0_realkind
    zads=min( 0.5_realkind*dtime/zlife , 1.0_realkind )
    zsec=1.E-7_realkind
    zcwmin=1.0E-7_realkind 
    zstun=0.5_realkind 
    zcpi=1.0_realkind/cpair 
    zmepsi=1.0_realkind -epsilo
    kit=2
    !     ------------------------------------------------------------
    !     2.2 arrays 
    !     ------------------------------------------------------------
    do jk=1,nlev
       do jl=kstart,kstop
          !        
          !         latent heat as a function of temperature.
          !
          jintt = max(1, 1 + nint((ht273 - wtp(jl,jk))/dttabl))
          jintt = min(jintt, 1750)
          wldcp(jl,jk) = zcpi * (latvap + hdl * prbice(jintt))
          !
       enddo
    enddo
    !


    !         3.  Do stratiform condensation computations.
    !         The previously updates of temperature, specific humidity 
    !         and cloud water are used, and also a preliminary 
    !         'supersaturation' specific humidity inside the 
    !         stratiform specific humidity inside the stratiform 
    !         clouds. 
    !
    !         First case: Cloud water evaporates if neither convective
    !         conditions are fulfilled ('klabcv' is different from 3')
    !         neither stratiform condensation can occur.
    !         --------------------------------------------------------
    !
    do jk=1,nlev
       jkp1=jk+1
       do jl=kstart,kstop
          !       ------------------
          !
          zpf=0.5_realkind*( ahyb(jk) +ahyb(jkp1) + (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk))
          wqsat(jl)=epsilo*esat(wtp(jl,jk))*zresat
          wdqsdt(jl)=epsilo*zpf*desdt(wtp(jl,jk))* zresat*zresat
          zqtot=wqp(jl,jk) +wcwp(jl,jk) 
          zqtot=max( zqtot, zcwmin )
          zqmax=zqtot*(1.0_realkind +wrhuc(jl,jk))
          !
          !         ----------------------------------------------------------
          if( kscall==1 ) then

             !
             if( (klabcv(jl,jk)/=3).and.(zqmax<wqsat(jl)) ) then 
                !           --------------------------------------------------------
                ! 
                zsev=zstun/( dtime*max( wqsat(jl)-wqp(jl,jk), zsec )) 
                zsev=min( zseva, zsev ) 
                zs1=dtime*zsev*wldcp(jl,jk)*wcwp(jl,jk) 
                zs2=wqsat(jl)                                       
                zs3=wdqsdt(jl)
                zs4=zs1/wldcp(jl,jk)
                !          
                zss1=(wtp(jl,jk) -zs1*zs2 +zs1*zs3*wtp(jl,jk))/ (1.0_realkind +zs1*zs3)
                zss2=zs1/(1.0_realkind +zs1*zs3)
                zss3=(wqp(jl,jk) +zs2*zs4 -zs3*zs4*wtp(jl,jk))/(1.0_realkind +zs4)
                zss4=zs3*zs4/(1.0_realkind +zs4) 
                zss5=wcwp(jl,jk) -zs2*zs4 +zs3*zs4*wtp(jl,jk)
                zss6=-zs3*zs4 
                zss7=zs4   
                wtp(jl,jk)=zss1 + zss2*(zss3 +zss1*zss4)/ (1.0_realkind -zss2*zss4)
                wqp(jl,jk)=(zss3 +zss1*zss4)/(1.0_realkind -zss2*zss4) 
                wcwp(jl,jk)=zss5 +zss6*wtp(jl,jk) +zss7*wqp(jl,jk)
                !
                !           -----
             endif
             !           -----
             !
             !         -------------------------------------------------------
          endif
          !       -------
       enddo
    enddo
    !     -----------------------------------------------------------
    !
    !     set integer counter 'jj' in iteration loop.
    jj=1 
    !


    !     Allow for 'kit' number of condensation iterations.
    !     ----------------------------------------------------------
1000 continue                        

    !

    !         3.2 Second case: Normal stratiform conditions 
    !         ---------------------------------------------
    !
    do jk=1,nlev
       jkp1=jk+1
       do jl=kstart,kstop
          !       ------------------
          !
          zpf=0.5_realkind*( ahyb(jk) +ahyb(jkp1) +   (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk))
          wqsat(jl)=epsilo*esat(wtp(jl,jk))*zresat
          wdqsdt(jl)=epsilo*zpf*desdt(wtp(jl,jk))* zresat*zresat 
          zqtot=wqp(jl,jk) +wcwp(jl,jk)
          zqtot=max( zqtot, zcwmin ) 
          zqmin=zqtot*(1.0_realkind -wrhuc(jl,jk))
          zqmax=zqtot*(1.0_realkind +wrhuc(jl,jk))
          !
          if( klabcv(jl,jk)/=3 ) then
             if( (zqmin<=wqsat(jl)).and.(wqsat(jl)<zqmax).and.(kscall==1) ) then
                zqceq=( wqsat(jl)*wqsat(jl)/zqtot +  (1.0_realkind +wrhuc(jl,jk))*&
                     (1.0_realkind +wrhuc(jl,jk))*zqtot - &
                     2.0_realkind*wqsat(jl)*(1.0_realkind +wrhuc(jl,jk)) )/&
                     (4.0_realkind*wrhuc(jl,jk))
                z1=(zqmax -wqsat(jl))/(2.0_realkind*wrhuc(jl,jk)*zqtot)
                zdqliq=(zqceq -wcwp(jl,jk))/ (1.0_realkind +wldcp(jl,jk)*wdqsdt(jl)*&
                     min( z1,1.0_realkind))

                wcwp(jl,jk)=wcwp(jl,jk) +zads*zdqliq
                wqp(jl,jk)=wqp(jl,jk) -zads*zdqliq
                wtp(jl,jk)=wtp(jl,jk) +wldcp(jl,jk)*zads*zdqliq
             endif
          endif
          if( wqsat(jl)<=wqp(jl,jk) ) then
             zqceq=zqtot -wqsat(jl)
             zdqliq=(zqceq -wcwp(jl,jk))/(1.0_realkind +wldcp(jl,jk)*wdqsdt(jl) )    

             wcwp(jl,jk)=wcwp(jl,jk) +zdqliq
             wqp(jl,jk)=wqp(jl,jk) -zdqliq
             wtp(jl,jk)=wtp(jl,jk) +wldcp(jl,jk)*zdqliq
          endif
       enddo
    enddo
    !     adjust iteration counter variable 'jj'
    jj=jj+1 
    if( jj<=kit) goto 1000

    !     end of stratiform iterations 

    !
    !     ---------------------------------------------------------

    !     4.  Allow cloud water to possibly evaporate completely
    !         Make associated correction of specific humidity.
    !         and temperature.
    !     ---------------------------------------------------------
    !
    do jk=1,nlev 
       do jl=kstart,kstop 
          if( wcwp(jl,jk)<zcwmin ) then
             wqp(jl,jk)=wqp(jl,jk) +wcwp(jl,jk)
             wtp(jl,jk)=wtp(jl,jk) -wldcp(jl,jk)*wcwp(jl,jk)
             wcwp(jl,jk)=0.0_realkind
          endif
       enddo
    enddo
    !
    !     ------
    return                                                                  
  end subroutine condst


  subroutine prevap( nhor,nlev,kstart,kstop,          &
       dtime,conacc,                           &
       ahyb,bhyb,dpf,                          &
       cucov,totcov,                           &
       draindt,dsnowdt,preta,                  &
       klabcv,                                 &
       wprcu,wprst,wcusno,wstsno,wpcoef,wmr,   &
       wdcusn,wdstsn,wrk1,wrk2,                &
       wpsp,wldcp,wqp,wtp,wcwp,wqsat,wdqsdt )
    !
    !     -------------------------------------------------------------------------

    !     --------   
    !     Purpose       Compute precipitation release and the evaporation of
    !     --------      precipitation. Also compute the associated modifications 
    !                   of cloud water, specific humidity and temperature.
    !
    !     -------------------------------------------------------------------------
    !     Interface  :  subroutine 'prevap' is called from subroutine 'condens'
    !     ------------
    !
    !     Input      :  ( from argument list or common blocks )
    !     ------------
    !       cw           cloud water content at (t-dt)
    !       dpf          delta-p around 'full' level
    !       pf           pressure at 'full' levels
    !       wpsp         provisional value of surface pressure
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !                    lifting is buoyant.
    !       klabcv       is a convective label. Layers with the same
    !                    value of 'klab(jl,jk)' can only be accepted as
    !                    convective. Convective cloud layers are first
    !                    assigned a value of klabcv(jl,jk)=2, and a final
    !                    value of 3. Convective sub-cloud layers are
    !                    assigned a value of 1. New lifting levels
    !                    however, are assigned a value of 4.
    !
    !     Output     :  ( new output or modified input) 
    !     ------------
    !       wldcp(jl,jk)  (1/cpair)*latent heat as a function of
    !                     temperature by weighting in latent heat of
    !                     freezing times probability of ice crystals
    !       wqsat(jl,jk)  saturation mixing ratio with respect to
    !                     gridpoint temperature
    !       wdqsdt(jl,jk) derivative of saturation water vapour pressure
    !                     with respect to temperature
    !       wqp          provisional value of specific humidity
    !       wtp          provisional value of temperature
    !       wcwp         provisional value of specific cloud water
    !
    !lr
    !lr     preta(jl,jk) precipitation intensity at model levels
    !lr
    !     Subroutines:  none
    !     ------------     
    !
    !     Author:       Bent H. Sass, (DMI), based on the research during the
    !     -------                     the HIRLAM I-III projects.
    !

    !     1. Declarations. 

    !
    implicit none                                                           

    integer nhor,nlev,kstart,kstop,jk,jl,jkp1,jintt
    integer klabcv(nhor,nlev)
    real(kind=realkind) dtime,conacc 
    real(kind=realkind) zpf,zdpf,zsec,zsec1,zecoef,zcpi,zdpdg,&
         zdqlim,                                &
         zpodif,zsqrt,zevapc,zevaps,            &
         zfbf,zprbi,zfcocv,zfcost,zftco,        &
         z1,zp1,zmepsi,zresat,                  &
         zprtc,zprts,                           &
         zcwcu,zcwst,zphase
    real(kind=realkind) ahyb(nlev+1),bhyb(nlev+1), accprc(nhor),accprl(nhor),&
         draindt(nhor),dsnowdt(nhor)
    real(kind=realkind) preta(nhor,nlev)
    real(kind=realkind) wpsp(nhor),wqsat(nhor),wdqsdt(nhor),                 &
         wprst(nhor),wprcu(nhor),wstsno(nhor),wcusno(nhor),    &
         wmr(nhor),wpcoef(nhor),                               &
         wdcusn(nhor),wdstsn(nhor),wrk1(nhor),wrk2(nhor)        
    real(kind=realkind) dpf(nhor,nlev),                                       &
         wldcp(nhor,nlev),wqp(nhor,nlev),                      &
         wtp(nhor,nlev),wcwp(nhor,nlev),                       &
         totcov(nhor,nlev),cucov(nhor,nlev)
    !


    !     2. Initialize various constants and arrays.

    !
    !     2.1 local constants   
    !     -------------------
    !     'zcpi'    inverse of specific heat capacity
    !     'zsec'    is a security constant.
    !     'zecoef'  is an evaporation coefficient applied 
    !               in computing evaporation of precipitation.
    !     'zphase'  melting coefficient to describe melting 
    !               of snow. 
    !     ---------------------------------------------------------
    !
    zsec=1.E-7_realkind
    zsec1=1.E-3_realkind
    zecoef=1.0E-3_realkind
    zcpi=1.0_realkind/cpair 
    zphase=4.0E-4_realkind
    zmepsi=1.0_realkind -epsilo 

    do jl=kstart,kstop
       !     initialize precipitation intensity arrays to zero.
       wprcu(jl)=0.0_realkind
       wprst(jl)=0.0_realkind
       wcusno(jl)=0.0_realkind
       wstsno(jl)=0.0_realkind      
       wrk1(jl)=0.0_realkind
       wrk2(jl)=0.0_realkind
       do jk=1,nlev
          preta(jl,jk)=0.0_realkind
       enddo
    enddo

    !         Start main loop to compute precipitation release and evaporation 
    !         of precipitation.
    do jk=1,nlev
       jkp1=jk+1
       do jl=kstart,kstop
          !         3. Evaporation of precipitation
          !         3.1 compute help parameters such as
          !             pressure in middle of a layer (zpf) 
          !             pressure thickness (zdpf), 
          !             thickness in terms of mass (zdpdg)
          !             latent heat divided with heat capacity (wldcp)
          !             saturation specific humidity (wqsat)
          !             the derivative of 'wqsat' with respect 
          !             to temperature (wdqsdt)  
          !         ---------------------------------------------------------
          !
          zpf=0.5_realkind*( ahyb(jk) +ahyb(jkp1) + (bhyb(jk) +bhyb(jkp1))*wpsp(jl) )

          zdpf=ahyb(jkp1) -ahyb(jk) +wpsp(jl)*(bhyb(jkp1) -bhyb(jk))
          zdpdg=zdpf/gravit    

          jintt = max(1, 1 + nint((tmelt - wtp(jl,jk))/dttabl))
          jintt = min(jintt, 1750)
          wldcp(jl,jk) =(latvap + latice*prbice(jintt))/cpair
          zresat=1.0_realkind/pecor(zpf,wtp(jl,jk))
          wqsat(jl)=epsilo*esat(wtp(jl,jk))*zresat
          wdqsdt(jl)=epsilo*zpf*desdt(wtp(jl,jk))*zresat*zresat

          !         3.2 Convective evaporation and stratiform evaporation are
          !             treated separately. Two constraints should apply:
          !
          !         i)  First,  the evaporation in a layer must not exceed
          !             the precipitation passing through the layer in a
          !             given time step.
          !         ii) Second, the evaporation must not exceed the satu-
          !             ration deficit ( wqsat -wqp )
          !         ---------------------------------------------------------
          !             Update specific humidity and temperature
          !             as a result of evaporation of convective
          !             precipitation.
          !         ---------------------------------------------------------
          !
          zpodif=max( wqsat(jl) -wqp(jl,jk), zsec )/(1.0_realkind +wldcp(jl,jk)*wdqsdt(jl) )
          zsqrt=sqrt( wprcu(jl) +zsec )
          zevapc=dtime*zecoef*zsqrt*zpodif
          zdqlim=wprcu(jl)*dtime*gravit/dpf(jl,jk)
          zevapc=min( zevapc,zdqlim )
          zevapc=min( zevapc,zpodif )

          wqp(jl,jk)=wqp(jl,jk) +zevapc
          wtp(jl,jk)=wtp(jl,jk) -zevapc*zcpi*(latvap +latice*wcusno(jl)/(wprcu(jl) +zsec))

          !             Update specific humidity and temperature
          !             as a result of evaporation of stratiform
          !             precipitation.
          !         --------------------------------------------------------
          !
          zpodif=max( wqsat(jl) -wqp(jl,jk), zsec )/(1.0_realkind +wldcp(jl,jk)*wdqsdt(jl) )
          zsqrt=sqrt( wprst(jl) +zsec )
          zevaps=dtime*zecoef*zsqrt*zpodif
          zdqlim=wprst(jl)*dtime*gravit/dpf(jl,jk)
          zevaps=min( zevaps,zdqlim )
          zevaps=min( zevaps,zpodif )
          !
          wqp(jl,jk)=wqp(jl,jk) +zevaps
          wtp(jl,jk)=wtp(jl,jk) -zevaps*zcpi*(latvap +latice*wstsno(jl)/(wprst(jl) +zsec))

          !         3.3 Update convective precipitation due to evaporation

          wprcu(jl)=wprcu(jl) -zdpdg*zevapc/dtime
          wprcu(jl)=max( wprcu(jl),0.00_realkind )

          !         3.4 Update stratiform precipitation due to evaporation
          !         -------------------------------------------------------
          !
          wprst(jl)=wprst(jl) -zdpdg*zevaps/dtime
          wprst(jl)=max( wprst(jl),0.00_realkind )
          !
          !


          !         4. Precipitation release 
          !         -------------------------------------------------------
          !         4.1 modified ice probability 'zprbi'
          !             (eq.2.2.3.20) in HIRLAM-2.5 documentation )
          !         --------------------------------------------------------- 
          !
          zprbi=prbice(jintt) +(1.0_realkind -prbice(jintt))*(wcusno(jl) +&
               wstsno(jl))/(wprcu(jl) +wprst(jl) +zsec)
          zprbi=max( 0.0_realkind, min( zprbi,1.0_realkind ) )
          zfbf=zprbi*(1.0_realkind -prbice(jintt))*hdewi(jintt)
          !
          !         ---------------------------------------------------------

          !         4.2 extra temperature dependency term 'zftco'
          !             (eq. 2.2.3.25)
          !         ---------------------------------------------------------
          !
          zftco=1.0_realkind
          if( wtp(jl,jk)<=238.0_realkind ) zftco=1.0_realkind +0.5_realkind*&
               (238.0_realkind -wtp(jl,jk))
          if( wtp(jl,jk)<=230.0_realkind ) zftco=5.0_realkind
          !
          zcwcu=wcwp(jl,jk)
          zcwst=wcwp(jl,jk)
          !
          !         ---------------------------------------------------------
          !
          if( klabcv(jl,jk)==3) then

             !

             !         4.3 take in to account collection effect and
             !             Bergeron-Findeisen effect in parameter 'zfcocv'
             !             (eq. 2.2.3.15) in HIRLAM-2.5 documentation )
             !         ---------------------------------------------------------
             !
             zfcocv=1.0_realkind +coales*sqrt((wprcu(jl)+zsec)/(totcov(jl,jk)+zsec1))+cbfeff*zfbf
             !

             !         4.4 compute cloud water humidity threshold 'wmr' 
             !             for convective precipitation release (eq. 2.2.3.22)
             !         ---------------------------------------------------------
             !
             wmr(jl)=hmrcu*( hmroft(jintt)/(zfcocv +zsec) )
             !
             !         ---------------------------------------------------------

             !         4.5 precipitation coefficient 'wpcoef' (eq. 2.2.3.14)     
             !         ---------------------------------------------------------
             !            
             wpcoef(jl)=hccu*zfcocv*zftco
             !        
             !         --------------------------------------------------------

             !         4.6 release of convective precipitation in the layer.
             !         -------------------------------------------------------
             !
             z1=wcwp(jl,jk)/((zsec +cucov(jl,jk))*wmr(jl))
             zp1=dtime*wpcoef(jl)*(1.0_realkind -exp( -z1*z1 ))
             !
             wcwp(jl,jk)=wcwp(jl,jk)/(1.0_realkind +zp1)
             !
             zprtc=zcwcu -wcwp(jl,jk)
             wprcu(jl)=wprcu(jl) + zprtc*zdpdg/dtime 
             !       
             !         --------------------------------------------------------

             !         4.7 Compute snow fraction
             !         -------------------------------------------------------
             !
             wcusno(jl)=wcusno(jl) + zprbi*zprtc*zdpdg/dtime
             !
             !         -------------------------------------------------------

             !         4.8 Release heat due to precipitation induced freezing
             !             of convective precipitation
             !         ------------------------------------------------------
             wtp(jl,jk)=wtp(jl,jk) + latice*zcpi*(zprbi -prbice(jintt))*zprtc
             !         -------------------------------------------------------
             !
          else 

             !

             !
             !         4.90 take in to account collection effect and
             !              Bergeron-Findeisen effect in parameter 'zfcost'
             !              (eq. 2.2.3.15) in HIRLAM-2.5 documentation )
             !         ---------------------------------------------------------
             !
             zfcost=1.0_realkind+coales*sqrt((wprst(jl)+zsec)/(totcov(jl,jk)+zsec1))+cbfeff*zfbf
             !
             !         4.9 Compute cloud water humidity threshold 'wmr'
             !             for stratiform precipitation release (eq. 2.2.3.22)
             !         ---------------------------------------------------------
             !
             wmr(jl)=hmrst*( hmroft(jintt)/(zfcost +zsec) )
             !
             !         ---------------------------------------------------------

             !         4.10 Precipitation coefficient 'wpcoef' (eq. 2.2.3.14)
             !         ---------------------------------------------------------
             wpcoef(jl)=hcst*zfcost*zftco
             !
             !         --------------------------------------------------------

             !         4.11 Release of convective precipitation in the layer.
             !         -------------------------------------------------------
             !
             z1=wcwp(jl,jk)/((zsec +totcov(jl,jk))*wmr(jl))
             zp1=dtime*wpcoef(jl)*(1.0_realkind -exp( -z1*z1 ))
             !
             wcwp(jl,jk)=wcwp(jl,jk)/(1.0_realkind +zp1)
             !
             zprts=zcwst -wcwp(jl,jk)
             wprst(jl)=wprst(jl) + zprts*zdpdg/dtime
             !
             !         --------------------------------------------------------

             !         4.12 Partition between rain and snow and release energy
             !              in association with freezing and melting of
             !              precipitation.
             !         -------------------------------------------------------
             !
             wstsno(jl)=wstsno(jl) + zprbi*zprts*zdpdg/dtime
             !
             !         -------------------------------------------------------

             !         4.13 Release heat due to precipitation induced freezing
             !              of stratiform precipitation
             !         -------------------------------------------------------
             wtp(jl,jk)=wtp(jl,jk) + latice*zcpi*(zprbi -prbice(jintt))*zprts
          endif


          !         5.  Adjust snow flux due to melting of precipitation 
          !             and release associated latent heat 
          !             Note that the case of re-freezing rain 
          !             is not considered at present. 
          wdcusn(jl)=zphase*cpair*zdpdg* max( wtp(jl,jk) -tmelt, 0.0_realkind )/latice
          wdcusn(jl)=-min( wdcusn(jl), wcusno(jl) ) 
          wcusno(jl)=wcusno(jl) +wdcusn(jl) 
          wtp(jl,jk)=wtp(jl,jk) +dtime*latice*zcpi*wdcusn(jl)/zdpdg 
          !
          wdstsn(jl)=zphase*cpair*zdpdg*max( wtp(jl,jk) -tmelt, 0.0_realkind )/latice
          wdstsn(jl)=-min( wdstsn(jl), wstsno(jl) )
          wstsno(jl)=wstsno(jl) +wdstsn(jl)
          wtp(jl,jk)=wtp(jl,jk) +dtime*latice*zcpi*wdstsn(jl)/zdpdg

          !lr     precipitation flux at a model level, for output

          preta(jl,jk)=wprst(jl)+wprcu(jl)
       enddo
    enddo

    !     6. Update rain intensity 'draindt', snow intensity 'dsnowdt',
    !        and accumulated precipitations 'accprc' and 'accprl'. 


    do jl=kstart,kstop
       draindt(jl)=wprst(jl) +wprcu(jl) -wstsno(jl) -wcusno(jl)
       dsnowdt(jl)=wstsno(jl) +wcusno(jl)
    enddo
    return
  end subroutine prevap



  subroutine inicons()

    !     inicons - setup of Sundqvist scheme parameters

    implicit none
    !     
    !     REFERENCES:
    !     
    !     Sundqvist, H. (1993).-"Inclusion of Ice Phase of Hydrometeors in Clou
    !     Parameterization for Mesoscale and Largescale Models". Beit. Phys.
    !     Atmosph., 66,137-147. Later will be S93.
    !     
    !     Olofsson, P-O.; Sundqvist, H.; Zurovac-Jevtic, D. (1994).-"Documentat
    !     of the routine condens". MISU. Later will be OSZ94.
    !     
    !     some index for loops
    !     
    integer i10,it,i
    logical debug

    !     local constants to build the tables

    real(kind=realkind) hedldr,tci,todpmx,tscale,demax,apri,test,tphalf,podpmx,&
         fxt,xti,x,x232,bfmax,xpt1mp

    debug=.false.
    hdl   =0.334e6_realkind                                                         

    heldr  = epsilo*latvap/rair                                             
    hedldr = epsilo*hdl/rair                                                
    !     
    !     
    !     he273  - saturation vapor pressure for t=273K
    !     tvirtc - constant needded to compute virtual temperature
    !     hkevap - coefficient for evaporation from stratiform precipitation
    !     u00max - maximum allowable value of modified hu00
    !     
    ht273  = 273.0_realkind                                                           
    he273  = 611.0_realkind
    tvirtc = 0.61_realkind
    hkevap = 5.e-5_realkind
    stpevp = 5.e-5_realkind                                                          
    u00max = 0.975_realkind
    !     

    !     hu00   - threshold relative humidity for stratiform condensation over
    !     tcir   - temperature below which cirrus (pure ice crystal) clouds are
    !     considered
    !     conae  - constant for the development in series of the equivalent
    !     potential temperature
    !     aecon  - exp(conae)
    !     
    hu00 = 0.85_realkind
    tcir  = 235.0_realkind
    conae = 0.15_realkind
    aecon = exp(conae)
    !     
    !     coales - parameter of the coalescence factor
    !     hccu   - parameter c00 for cumulus, that is the conversion rate of
    !     cloud to precipitation drops in convective clouds
    !     hcunrm - that is cloud cover depens on cloud depth compared to hcunrm
    !     hmrcu  - mr for convective clouds, that is cloud water mixing ratio
    !     at which conversion becomes efficient in convective clouds
    !     hmrst  - mr for stratiform clouds, that is cloud water mixing ratio
    !     at which conversion becomes efficient in stratiform clouds
    !     htaucu - characteristic time used in convective cloud cover scheme
    !     
    cfreez = 0.12_realkind                                                           
    coales = 100.0_realkind
    hccu = 2.e-4_realkind
    hcunrm = 3.e4_realkind
    hmrcu = 8.e-4_realkind
    hmrst = 4.e-4_realkind
    htaucu = 3600.0_realkind

    !     hp0    _realkind- reference pressure (=100 kpa)
    !     cbfeff - parameter for the Bergeron-Findeisen effect
    !     hvterm - terminal velocity of precipitation
    !     hvsnow - terminal velocity of snow
    !     
    hp0    = 1.e5_realkind
    cbfeff = 4.0_realkind
    hvterm = 5.0_realkind
    hvsnow = 1.0_realkind
    !     
    !     asnow  - a parameter together with bsnow and snoref to govern the
    !     rate of ice precip before the Bergeron-Findeisen mechanism
    !     becomes effective, even if the precip from above is pure
    !     ice;
    !     the factor cbfsno = asnow*(snowrate/snoref)**bsnow
    !     and cbfsno = cbfsno/(1+cbfsno)
    !     which multiplies the modified ice probability
    !     bsnow  - see asnow
    !     snoref - see asnow
    !     
    asnow = 9.0_realkind
    bsnow = 4.0_realkind
    snoref = 0.05_realkind

    !     hkmelt - coefficient for melting of ice in precipitation
    !     tanvil - temperature below which convective anvil is allowed and
    !     treated as stratiform
    !     hcst   - conversion rate of cloud to precipitation drops in stratifor
    !     clouds
    !     pmoist - exponent for computations of hmoist (eq. 3 in OSZ94)
    !     
    hkmelt = 2.5e-4_realkind
    tanvil = 253.0_realkind
    hcst = 1.e-4_realkind
    pmoist = 3.0_realkind
    !     
    !     hkap   - rair/cpair
    !     hdldcp - latice/cpair
    !     hecdr  - epsilo / hkap
    !     hldcp  - latvap/cpair
    !     
    hkap   = rair/cpair
    hdldcp = latice/cpair
    hecdr  = epsilo / hkap
    hldcp  = latvap/cpair

    !     coalcu - coalescence factor for convection clouds (c0.mr/2g)
    !     coalst - coalescence factor for convection clouds (c0.mr/2g)
    !     sqvsno - square root of the snow terminal velocity
    !     heldr  - epsilo*latvap/rair
    !     hedldr - epsilo*latice/rair
    !     elotci - epsilo*latice/rair/tcir
    !     
    coalcu = 0.5_realkind* hccu* hmrcu / gravit
    coalst = 0.5_realkind* hcst* hmrst / gravit
    sqvsno = sqrt(hvsnow)
    heldr  = epsilo*latvap/rair
    hedldr = epsilo*latice/rair
    elotci = epsilo*latice/rair/tcir

    !     dttabl - temperature increment to build the table
    !     prbice - porbability of ice (eq. 1 in S93)
    !     bfeff  - first value of the Bergeron-Findeisen effect (without take
    !     into account the ice precipitation) (eq. 13 in S93)
    !     hdewi  - difference between both saturation vapor pressure over water
    !     over ice
    !     hmroft - temperature function to multiply mr for T<273 (eq. 25 in OSZ
    !     
    !     
    !     Parameter values in si units
    !     ----------------------------
    !     
    tci    = 232.0_realkind
    todpmx = 299.0_realkind
    tscale = (todpmx - tci)*sqrt(2.0_realkind)
    !     
    !     the table look up is done by finding the integer
    !     it = 1 + nint((t- h273)/dttabl)
    !     
    dttabl = 0.1_realkind
    !     
    !     calculate saturation vapor pressure to be used as table
    !     
    !     calculate difference of saturation pressure between
    !     water and ice, dewi, for t< t0
    !     calculate a probability function, prbice, for existence of
    !     ice crystals
    !     calculate a temperature function, hmroft,
    !     to multiply hmrcu and hmrst
    !     
    !     demax  - maximum value of hdewi
    !     apri   - parameter A in S93
    !     test   - 2A/(2A-1)
    !     tphalf - temperature with a probalility of ice equal to 0.5 (T3 in S9
    !     podpmx - ice probability for todpmx
    !     
    demax = 0.0_realkind
    apri = 1.0_realkind/(1.0_realkind-exp(-((tmelt-tci)/tscale)**2.0_realkind))
    test = 2.0_realkind*apri/(2.0_realkind*apri-1.0_realkind)
    if(test <= 0.0_realkind)  then
       write(6,*)'warning: test is < 0 !!!! in calc of table for'
       write(6,*)'  ice crystal probability'
    endif
    tphalf = log(2.0_realkind*apri/(2.0_realkind*apri-1.0_realkind))
    tphalf = tci + tscale*sqrt(tphalf)
    if(tphalf<0.0_realkind) then
       write(6,*) 'before do 25  in calc of table for ice'
       write(6,*) '   crystal probability       tphalf=',tphalf
    endif
    podpmx = apri*(exp(-0.5_realkind)-1.0_realkind)+1.0_realkind
    !     

    !     table from T=273 up to T=99K
    !     
    !     fxt    - temperature
    !     xti    - (1/T0-1/T)
    !     hdewi  - difference between both saturation vapor pressure over water
    !     over ice
    !     prbice - porbability of ice (eq. 1 in S93)
    !     hmroft - temperature function to multiply mr for T<273 (eq. 25 in OSZ
    !     
    do 25 i = 1,1750
       !     
       fxt = tmelt - (real(i,realkind)-1.0_realkind) * dttabl
       xti = 1._realkind/tmelt - 1.0_realkind/fxt
       hdewi(i)  = he273/fxt * exp(heldr*xti) * (1.0_realkind- exp(hedldr*xti))
       demax = max(demax,hdewi(i))
       !     
       prbice(i) = apri * (exp(-(((fxt - tci)/tscale)**2.0_realkind)) - 1.0_realkind) + 1.0_realkind
       prbice(i) = max( prbice(i) , 0.0_realkind)
       if(fxt < tci) prbice(i) = 1.0_realkind
       prbice(i) = min(prbice(i) , 1.0_realkind)
       if(fxt > 250.0_realkind) hmroft(i) = 1.33_realkind*exp(-((fxt-tmelt)*&
            0.066_realkind)**2.0_realkind)
       !     
       if(fxt <= 250.0_realkind) then
          x = abs (fxt - 232.0_realkind) / 18.0_realkind
          x = x * (1.0_realkind + x * (1.0_realkind + 1.333_realkind * x))
          x232 = 1.0_realkind
          if (fxt < 232.0_realkind) x232 = -1.0_realkind
          x = x / (1.0_realkind + x) * x232
          hmroft(i) = 0.5_realkind*0.15_realkind*(1.07_realkind+x)
       endif
       hmroft(i) = min( hmroft(i) , 1.0_realkind)

       !     hmroft - becomes fmr (eq. 24 in OSZ94)
       !     
       hmroft(i) = (1.0_realkind-prbice(i))**2.0_realkind + prbice(i) * hmroft(i)
       hmroft(i)   = max(hmroft(i) , 3.e-2_realkind)
25  enddo
    !     
    !     setting values for i=1
    !     
    prbice(1) = 0.0_realkind
    hmroft(1) = 1.0_realkind
    !     
    !     normalize with max difference, demax
    !     calculate the product of dewi and prbice, called bfeff
    !     
    !     bfmax  - maximum value of bfeff
    !     
    bfmax = 0._realkind

    !     bfeff  - first value of the Bergeron-Findeisen effect (without take
    !     into account the ice precipitation) (eq. 13 in S93)
    !     
    do 30 i = 1,1750
       !     
       hdewi(i)  = hdewi(i)/demax
       bfeff( i) = prbice(i) * (1.0_realkind-prbice(i)) * hdewi(i)
       bfmax = max(bfmax,bfeff(i))
       !     
30  enddo
    !     
    !     normalizing bfeff with hir maximum value
    !     
    do 32 i = 1,1750
       bfeff( i) = bfeff( i) / bfmax
32  enddo

    !     print a selection of the computed values
    !     

    if(mype==0.and.debug) then
       write(6,605) tmelt,tphalf,podpmx,todpmx
       write(6,615)
       do 60 i10 = 1,47,3

          it = 1 - i10
          i  = 1 + (i10-1) * 10
          xpt1mp = prbice(i)*(1._realkind-prbice(i))
          write(6,617)it,hdewi(i),prbice(i),bfeff(i),hmroft(i),xpt1mp

60     enddo
    endif

600 format(//,' from inicons:',/)
605 format('   PROBABILITY TABLES BETWEEN 273 AND 99 K FOR:',//,                &
         '     DIFFERENCE of saturation pressure over WATER and ICE;',/,         &
         '     TEMPERATURE FUNCTION to multiply hmrst and hmrcu;',//,            &
         '     ICE CRYSTAL PROBABILITY P = 0 for T = ',f5.1,'  and P = 0.50      &
         for T = ',f5.1,/,20x,'     and P and T of dP/dT-max are P ='            &
         ,f5.2,'  and T = ',f5.1,/)

615 format(9x,'t',5x,'de',6x,'prob-ice',5x,'bf-f',4x,'hmroft',5x,'p*(1-p)')
617 format(i10,f9.4,4x,f7.4,4x,f7.4,2x,f7.4,4x,f7.4)

    !     Print-out of some of the parameters
    if (mype==0) then

       write(6,3) tanvil,hu00,htaucu,hcunrm,hccu,hmrcu, coales,hkmelt,pmoist

3      format(1h ,//,' FROM CONDENS with ice phase, part CONVEC:'         &
            ,/,'  convection scheme including anvil cloud calculation      &
            for T(cutop) <=',f5.0,' K',//,' other parameters are',/,     &
            '   u0   taucu    cunrm     cuc0      cumr    coales   hkm     &
            elt    pmoist',                                                &
            /,f6.2,f7.0,f10.0,e10.2,e10.2,f7.2,e11.2,f7.1,//)
       !     
       write(6,4) hu00,hkevap,cbfeff,hcst,hmrst,coales,hkmelt,pmoist,asnow,bsnow,snoref
4      format(1h ,//,' FROM CONDENS with ice phase, part STCOND:'              &
            ,/,'  stratiform cloud of sub-grid scale, i.e., even for u < 1.',// &
            ,' other parameters are',&
            '   u0     kevp    cbfeff         stc0      stmr    coales   hkmelt    pmoist' &
       , /,f6.2,e10.2,f6.1,2x,e9.2,e10.2,f7.2,f8.0,e11.2,f7.1 &
            ,//,'  asnow',3x,'bsnow',3x,'snoref mm/h',/,3x,f4.1,3x,f4.1  ,5x,f6.3,//)

       write(6,*) '  '
       write(6,*) '  '

    endif

    snoref = snoref/3600.0_realkind


    return
  end subroutine inicons

end module condsmod
