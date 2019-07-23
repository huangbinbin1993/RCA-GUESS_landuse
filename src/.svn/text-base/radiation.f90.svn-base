module radiation
  use confys
  use co2mod
  use decomp
  use gcm
  use ccons
  implicit none
  private
  integer,parameter::mtabm=5000
  integer,save:: ntabm
  real(kind=realkind),save:: cqpmin=2.e-5_realkind !     minimum water path length
  real(kind=realkind),save::cpat !computed in inirad
  real(kind=realkind),save::c10e=0.4343_realkind !     conversion factor between logarithmic functions.
  real(kind=realkind),save::cem1=0.60_realkind !     emissivity constants
  real(kind=realkind),save::cem2=0.17_realkind !     emissivity constants
  real(kind=realkind),save::cem3=0.0082_realkind !     emissivity constants
  real(kind=realkind),save::cem4=0.0045_realkind !     emissivity constants
  real(kind=realkind),save::ckilgr=1000._realkind !conversion factor from kilo to gram.
  real(kind=realkind),save::csur1=35._realkind !emission constants for surface flux calculations
  real(kind=realkind),save::csur2=3000._realkind !emission constants for surface flux calculations
  real(kind=realkind),save::csur3=1.00_realkind !emission constants for surface flux calculations
  real(kind=realkind),save::cosmin=0.00001_realkind
  real(kind=realkind),save::cdegr=0.0174533_realkind
  real(kind=realkind),save::cosd !computed in inirad
  real(kind=realkind),save::sind !computed in inirad
  real(kind=realkind),save::cpatsh=0.05_realkind
  real(kind=realkind),save::caadry=0.05_realkind
  real(kind=realkind),save::caahum=0.03_realkind
  real(kind=realkind),save::cbbdry=0.63_realkind
  real(kind=realkind),save::cbbhum=0.81_realkind
  real(kind=realkind),save::sa   !computed in inirad
  real(kind=realkind),save::csac=1.e-9_realkind
  real(kind=realkind),save::ctrmin=0.10_realkind
  real(kind=realkind),save::calbsn
  real(kind=realkind),save::csndep
  real(kind=realkind),save::cemc1=0.25_realkind !     emissivity constants in connection with clouds
  real(kind=realkind),save::cemc2=0.20_realkind !     emissivity constants in connection with clouds
  real(kind=realkind),save::cfrac1=1.0_realkind
  real(kind=realkind),save::cfrac2=2.0e4_realkind
  real(kind=realkind),save::cetyp=11.5_realkind
  real(kind=realkind),save::caak=1.20_realkind
  real(kind=realkind),save::cask=1.25_realkind
  real(kind=realkind),save::chybcr=0.05_realkind
  real(kind=realkind),save::cpf=2500._realkind
  real(kind=realkind),save::cpc=7.0e4_realkind
  real(kind=realkind),save::cpsref=100000._realkind !should this be the same as sip0??
  real(kind=realkind),save::cdpsur=2500._realkind
  real(kind=realkind),save::crefrad=1.66_realkind
  real(kind=realkind),save::shtab(mtabm) !computed in inirad
  !END COMRAD

  real(kind=realkind),save::cqco2 
  real(kind=realkind),save::ceml1 = 6.54860_realkind 
  real(kind=realkind),save::ceml2 = -11.9669_realkind
  real(kind=realkind),save::ceml3 = 7.82396_realkind
  real(kind=realkind),save::ceml4 = -1.79180_realkind 
  real(kind=realkind),save::cpate != 10.**(-2.5)non-constant initialization expression at

  real(kind=realkind),save::ceme1a = 31.811_realkind
  real(kind=realkind),save::ceme2a = -1463.1_realkind
  real(kind=realkind),save::ceme1b = -8.2253e-2_realkind
  real(kind=realkind),save::ceme2b = 3.1258_realkind
  real(kind=realkind),save::ceme3b = -2.1966_realkind
  real(kind=realkind),save::ceme4b = -3.6719_realkind
  real(kind=realkind),save::ceme5b = 4.0475_realkind
  real(kind=realkind),save::cpatc  = 0.001_realkind
  real(kind=realkind),save::cemc1a = 4.5297_realkind
  real(kind=realkind),save::cemc2a = -109.83_realkind
  real(kind=realkind),save::cemc3a = 1264.4_realkind
  real(kind=realkind),save::cemc1b = 0.21332_realkind
  real(kind=realkind),save::cemc2b = 2.8000e-2_realkind
  real(kind=realkind),save::cemc3b = 2.2341e-3_realkind
  real(kind=realkind),save::cemc4b = 1.6113e-4_realkind
  real(kind=realkind),save::ca1 = 1.8062_realkind
  real(kind=realkind),save::ca2 = -8.0828_realkind
  real(kind=realkind),save::ca3 = 10.029_realkind
  real(kind=realkind),save::ca4 = -4.948_realkind
  real(kind=realkind),save::ca5 = 0.88512_realkind
  real(kind=realkind),save::cb1 = 0.11293_realkind
  real(kind=realkind),save::cb2 !computed in inirad
  real(kind=realkind),save::cb3 !computed in inirad
  real(kind=realkind),save::cb4 !computed in inirad
  real(kind=realkind),save::cb5 !computed in inirad
  real(kind=realkind),save::cemadd = 0.12_realkind
  !END COMRAD2
  real(kind=realkind),save::cmmy =1.0e+06_realkind
  real(kind=realkind),save::remn =4._realkind
!cgj040711  real(kind=realkind),save::remn =7._realkind
  real(kind=realkind),save::remx =20._realkind
!  real(kind=realkind),save::creland = 9.528e-5_realkind
  real(kind=realkind),save::creland = 1.114e-4_realkind
  real(kind=realkind),save::cresea = 1.434e-4_realkind
  real(kind=realkind),save::cresice =1.678e-4_realkind
  real(kind=realkind),save::b10a  = 1.55e-4_realkind
  real(kind=realkind),save::b10b = 8.18e-3_realkind
  real(kind=realkind),save::b11  = 1.29_realkind
  real(kind=realkind),save::b12 = 0.545_realkind
  real(kind=realkind),save::b13a = 7.00_realkind
  real(kind=realkind),save::b13b = -4.75_realkind
  real(kind=realkind),save::b14 = 8.30e-2_realkind
  real(kind=realkind),save::a1 =0.522_realkind
  real(kind=realkind),save::a2 =-4.551_realkind
  real(kind=realkind),save::a3 =4.115_realkind
  !END COMREFF

  public inirad,aradia,partly_solarupdate
contains

  subroutine inirad(yearc, monthc,dayc,hourc,minc,secc,&
       nlev,hybf,hybh,emc,csusa,fabso3,cadd)
    use calendar
    use referenceParameters, only:sip0
    implicit none

    integer,intent(in):: nlev
    integer,intent(in)::yearc, monthc,dayc,hourc,minc,secc
    real(kind=realkind),intent(in):: hybf(nlev),hybh(nlev+1)
    real(kind=realkind):: emc(nlev),csusa(nlev)
    real(kind=realkind):: fabso3(nlev),cadd(nlev)
    integer:: j,jk,jkp1
    real(kind=realkind):: zgravi,zpref,zday,zdec,zincr
    real(kind=realkind):: x,cabs1,cabs2,zco2
    integer::ndaysYear

    !l    nomenclature

    !l    input:

    !l    parameter       content

    !l    monthc      :   current month  (utc)
    !l    dayc        :   current day    (utc)
    !l    hourc       :   current hour   (utc)
    !l    minc        :   current minut  (utc)
    !l    secc        :   current second (utc)
    !l    nlev        :   number of vertical model layers.
    !l    hybf,hybh   :   hybrid coordinates in full and half
    !l                    levels

    !l    output:

    !l    common block "comrad"
    !l    with all variables and arrays initialized.

    !l    author: bent h.sass ,
    !l    the danish meteorological institute ,      november 1992
    !l    modified by laura rontu and bent h.sass    february 1994
    !lrcf modified by carl fortelius and laura rontu august   1998
    !lrkw modified by laura rontu and klaus wyser    december 1998
    !lrpr modified by laura rontu and petri räisänen october  2000
    !gjlr modified by colin jones and laura rontu    june     2002

    !l    purpose:

    !l    compute constants applied in the radiation scheme of
    !l    hannu savijaervi.
    !l    -------------------------------------------------------

    !l    interface:

    !l    subroutine "inirad" is called from subroutine "iniphys"
    !l    which is called from subroutine "phcall".
    !l    -------------------------------------------------------
    !l    constants used in longwave radiation.

    !
    !     minimum water path length "cqpmin"
    !     -------------------------------------------------------
    cqpmin=2.e-5_realkind
    !     -------------------------------------------------------
    !     scaling factor "cpat" involved in calculating water
    !     path in cm of precipitable water.
    !     -------------------------------------------------------
    zgravi=gravit
    zpref= sip0 
    cpat=1._realkind/( zgravi*zpref*10._realkind )
    !     -------------------------------------------------------
    !     conversion factor "c10e" between logarithmic functions.
    !     -------------------------------------------------------
    c10e=0.4343_realkind
    !     -------------------------------------------------------
    !     emissivity constants "cem1,cem2,cem3,cem4"
    !     -------------------------------------------------------
    cem1=0.60_realkind
    cem2=0.17_realkind
    cem3=0.0082_realkind
    cem4=0.0045_realkind
    !     -------------------------------------------------------

    !     emissivity constants "cemc1,cemc2" in connection
    !     with clouds.
    !     emissivity constants
    !     -------------------------------------------------------
    cemc1=0.25_realkind
    cemc2=0.20_realkind
    cfrac1=1.0_realkind
    cfrac2=2.0e4_realkind
    cetyp=11.5_realkind
    !lrpr   cadd=2.3e-6  ! 'cadd' defined as a function of height, see below
    !     -------------------------------------------------------
    do jk=1,nlev
       !
       jkp1=jk+1
       !       -----------------------------------------------------

       !       emissivity factor "emc(jk)" used in the determination
       !       of cloud emissivity.
       !       -----------------------------------------------------
       emc(jk)=0.05_realkind
       if( hybf(jk)>cemc1 ) then
          emc(jk)=0.05_realkind + cemc2*(hybf(jk)-cemc1)/(1._realkind-cemc1)
       endif
       !       -----------------------------------------------------
       !       constant array "csusa(jk)" involved in estimating
       !       layer integrated cloud water, when cloud water is not
       !       an explicit forecast variable.
       !       -----------------------------------------------------
       csusa(jk)=0.10_realkind -0.05_realkind*hybf(jk)
       !
    enddo
    !     -------------------------------------------------------

    !     emission constants "csur1,csur2,csur3" for surface flux
    !     calculations. values of 'csur1' and 'csur2' have been
    !     tuned by b.h.sass 30/3 1993 to produce improved downward
    !     radiative flux in the atmosphere.
    !     -------------------------------------------------------
    csur1=35._realkind
    !bhs removed by lrpr     tuning action associated with longwave radiation
    csur2=3000._realkind
    !bhs removed by lrpr     csur2=750.
    csur3=1.00_realkind
    !     -------------------------------------------------------
    !     conversion factor "ckilgr" from kilo to gram.
    !     -------------------------------------------------------
    ckilgr=1000._realkind
    !

    !     constants and arrays used in shortwave radiation.

    !     solar radiative flux "sa" on top of the atmosphere
    !     as a funtion of day of the year.
    !     also determine cosine "cosd" and the sine "sind" of
    !     the sun's declination.
    !     formulae are from paltridge and platt 1976;
    !     solar constant is 1365 w/m2.
    !     -------------------------------------------------------
    !
    !

    ndaysYear = nDaysInYear(yearc)
    zday=(real(monthc-1,realkind)*30.4375_realkind + real(dayc,realkind) &
         -1.9_realkind)*2._realkind*pi/real(ndaysYear,realkind)
!         -1.9_realkind)*2._realkind*pi/365._realkind
    zdec= 0.006918_realkind - 0.399912_realkind*cos(zday) + 0.070257_realkind*sin(zday) &
         -0.006758_realkind*cos(2._realkind*zday) + 0.000907_realkind*sin(2._realkind*zday) &
         -0.002697_realkind*cos(3._realkind*zday) + 0.001480_realkind*sin(3._realkind*zday)
    !
    !
    !
    cosd=cos(zdec)
    sind=sin(zdec)
    sa=solar*(1._realkind +.034221_realkind*cos(zday)+.00128_realkind*sin(zday)+ &
         .000719_realkind*cos(2._realkind*zday))
    csac=1.e-9_realkind
    !
    cosmin=0.00001_realkind
    cpatsh=0.05_realkind
    ctrmin=0.10_realkind
    !
    cdegr=0.0174533_realkind
    !lrb  snow variables are not used
    !lr         csndep=0.015
    !lr         calbsn=0.75
    !lre  and were removed from common block
    caadry=0.05_realkind
    caahum=0.03_realkind
    cbbdry=0.63_realkind
    cbbhum=0.81_realkind
    !     ---------------------------------------------------------

    !     table "shtab(j)" applied in humid clear air absorbtion
    !     ---------------------------------------------------------
    !
    ntabm=mtabm
    zincr=0.01_realkind
    !
    do j=1,ntabm
       shtab(j)=( real(j,realkind)*0.5_realkind*zincr )**(-cbbhum)
    enddo
    !
    caak=1.20_realkind
    !gjrca3      caak=1.30
    cask=1.25_realkind
    !gjrca3      cask=1.35
    chybcr=0.05_realkind
    cpf=2500._realkind
    !bhs990521-b
    !     correct coefficient associated with cooling
    !     rate below clouds, to get better agreement
    !     with reference comparisons for thick cirrus clouds.
    !     cpc=4.0e4
    !
    cpc=7.0e4_realkind
    !bhs990521-e
    cpsref= 100000._realkind
    cdpsur=2500._realkind
    crefrad=1.66_realkind
    !     ---------------------------------------------------------
    !lrpr990328-begin constants and coefficients for the modified 
    !lrpr             computation of longwave clear-sky emissivities

    !    --------------------------------------------
    !    modify the values of csur1 and chybcr
    !    --------------------------------------------
    !      csur1 = 13.
    !mtr991203 (pr)
!!!!!!!!!!!!!!!!!!!OVERRIDES PREVIOUS DEFINITION
    csur1 = 15._realkind
    csur1 = 13._realkind	!cgj280811
    chybcr= 0._realkind

    !    ------------------------------------------- 
    !    carbon dioxide volume and mass mixing ratio
    !    -------------------------------------------
    !      zco2 = 882.5e-6   ! co2 hadley fakt=2.5  
    !      zco2 = 706.e-6 ! co2 echam, era  fakt=2.0
    !      zco2 = 353.e-6 ! cntr hadley, echam, era

    zco2 = const_co2 * 1.e-6_realkind 
    cqco2 = zco2 *44._realkind/28.9644_realkind                         
    !      print *,'co2:',zco2

    !     ----------------------------------------------------------
    !     water vapor line (&p-type continuum) emissivity parameters
    !     ----------------------------------------------------------
    ceml1 = 6.54860_realkind
    ceml2 = -11.9669_realkind
    ceml3 = 7.82396_realkind
    ceml4 = -1.79180_realkind  

    !     ----------------------------------------------------------
    !     water vapor e-type continuum emissivity parameters
    !     ----------------------------------------------------------

    cpate  = 10._realkind**(-2.5_realkind)
    ceme1a = 31.811_realkind
    ceme2a = -1463.1_realkind
    ceme1b = -8.2253e-2_realkind
    ceme2b = 3.1258_realkind
    ceme3b = -2.1966_realkind
    ceme4b = -3.6719_realkind
    ceme5b = 4.0475_realkind

    !     ----------------------------------------------------------
    !     carbon dioxide emissivity parameters
    !     ----------------------------------------------------------

    cpatc  = 0.001_realkind
    cemc1a = 4.5297_realkind
    cemc2a = -109.83_realkind
    cemc3a = 1264.4_realkind
    cemc1b = 0.21332_realkind
    cemc2b = 2.8000e-2_realkind
    cemc3b = 2.2341e-3_realkind
    cemc4b = 1.6113e-4_realkind

    !     -----------------------------------------------------------    
    !     coefficients for water vapour line and continuum absorption
    !     overlap factor
    !     -----------------------------------------------------------

    ca1 = 1.8062_realkind
    ca2 = -8.0828_realkind
    ca3 = 10.029_realkind
    ca4 = -4.948_realkind
    ca5 = 0.88512_realkind

    !     ---------------------------------------------------------------    
    !     coefficients for water vapour and carbon dioxide overlap factor
    !     ---------------------------------------------------------------
    x = 2.1449_realkind*((zco2/500.e-6_realkind)**(1._realkind/3._realkind)-1._realkind)
    cb1 = 0.11293_realkind
    cb2 = 4.4780_realkind +1.4774e-2_realkind*x
    cb3 = -9.777_realkind -8.5438e-2_realkind*x
    cb4 = 6.313_realkind +0.19318_realkind*x
    cb5 = -1.317_realkind -8.1497e-2_realkind*x  

    !     ---------------------------------------------------------------    
    !     typical contribution of other factors than h2o line absorption
    !     to single-layer clear-sky emissivity:
    !     --------------------------------------------------------------- 
    cemadd = 0.12_realkind

    !     ----------------------------------------------------------------
    !     prescribe the vertical distribution of ozone shortwave absorption
    !     ----------------------------------------------------------------    

    cabs1 = 0._realkind
    do jk =2,nlev+1
       x = log(1._realkind+1000._realkind*hybh(jk))
       cabs2 = 1._realkind-exp(x*(-6.52978e-3_realkind+x*(-0.178512_realkind+x*(7.13869e-2_realkind &
            +x*(-1.45003e-2_realkind+x*7.36790e-4_realkind)))))
       fabso3(jk-1) = cabs2-cabs1
       cabs1 = cabs2
    enddo

    !     ------------------------------------------------------------------
    !     tuning factor for longwave cooling (now vertically non-constant)
    !     ------------------------------------------------------------------

    do jk = 1,nlev
       !pr/mtr, 991012>
       !        cadd(jk)=1.16e-6*(6.*hybf(jk)-2.)
       !
       cadd(jk)=1.5e-6_realkind*(3._realkind*hybf(jk)-1._realkind)
    enddo

    !lrpr990328-end
    ! lrkw
    !     ---------------------------------------------------------
    !  definitions needed for effective radius parameterization
    !  filling common block comrad
    !
    cmmy=1.0e+06_realkind
    remn=4._realkind
!cgj040711    remn=7._realkind
    !gj        remn=5.5
    remx=20._realkind
    !gj090605
    !gj090605 define creland etc to be ((3/4*pi*density of water)/(kn))**0.333
    !gj090605
!    creland=9.528e-5_realkind
    creland=1.114e-4_realkind
    cresea= 1.434e-4_realkind
    cresice=1.678e-4_realkind
    !gj090605
    b10a = 1.55e-4_realkind
    b10b = 8.18e-3_realkind
    b11 = 1.29_realkind
    b12 = 0.545_realkind
    b13a = 7.00_realkind
    b13b = -4.75_realkind
    b14 = 8.30e-2_realkind
    a1=0.522_realkind
    a2=-4.551_realkind
    a3=4.115_realkind

    return
  end subroutine inirad

  subroutine partly_solarupdate(monthc,dayc,hourc,mnt,sec)
    implicit none

    integer,intent(in)::monthc,dayc,hourc,mnt,sec
    real(kind=realkind)::zday,zdec

    logical lfirst
    data lfirst/.true./
    save lfirst

    if (lfirst) then
       if(mype==0)then
          write(6,*) 'solar constant set to ',solar
       endif
       lfirst=.false.
    endif
!!$    if(mype==0)then
!!$       print *,'call 2 solarupdate',monthc,dayc,hourc,mnt,sec
!!$    endif
    
    if(lhadley.or.lecham)then
       zday=(real(monthc-1,realkind)*30._realkind + real(dayc,realkind) )*2._realkind*pi/360._realkind
    else
       zday=(real(monthc-1,realkind)*30.4375_realkind +  real(dayc,realkind) - 1.9_realkind)*2._realkind*pi/365.25_realkind
    endif

    zdec= 0.006918_realkind - 0.399912_realkind*cos(zday) + 0.070257_realkind*sin(zday) &
         -0.006758_realkind*cos(2._realkind*zday) + 0.000907_realkind*sin(2._realkind*zday) &
         -0.002697_realkind*cos(3._realkind*zday) + 0.001480_realkind*sin(3._realkind*zday)

    cosd=cos(zdec)
    sind=sin(zdec)
    sa=solar*(1._realkind +.034221_realkind*cos(zday)+.00128_realkind*sin(zday)+ &
         .000719_realkind*cos(2._realkind*zday))


    return
  end subroutine partly_solarupdate


  subroutine aradia( nhor, nlev, kstart, kstop,             &
       hourc,   minc,   secc,     locw,   hybf,   hybh, &
       emc,  csusa, fabso3,   cadd,  along, coslat, sinlat,  &
       albedo, tsk,  t,  q, cw, pf,     ph, frland,emskin, &
       dpf, totcov,   dtdt, &!  strad, stdrad,  stclo, stdclo, &
       accsunny,dtphysh, &
       slwdn ,  sswdn, wslwnet, wsswnet, tlwnet, tswnet, tswdn, &
       radf,frice,ps,snowh, &
!cgj300611
       wfrop,vegopl, &
       texture,gpot,pblh,lscov,cucov,cwcu,zwin)
!cgj300611

    implicit none
    !
    integer,intent(in):: nhor,nlev,kstart,kstop,hourc,minc,secc
    logical locw
    !
    real(kind=realkind) hybf(nlev),hybh(nlev+1),emc(nlev),csusa(nlev)
    !strad(nlev),stdrad(nlev),stclo(nlev),stdclo(nlev)
    real(kind=realkind) fabso3(nlev), cadd(nlev), scos(nhor)
    real(kind=realkind) frland(nhor),along(nhor), &
         coslat(nhor),sinlat(nhor),albedo(nhor),tsk(nhor)
    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),    &
         pf(nhor,nlev),ph(nhor,nlev+1),dpf(nhor,nlev), &
         totcov(nhor,nlev),dtdt(nhor,nlev), &
         slwdn(nhor) ,  sswdn(nhor) , &
         tlwnet(nhor), tswnet(nhor), tswdn(nhor) ,  salbcor(nhor)
!cgj300611
    real(kind=realkind) lscov(nhor,nlev),cucov(nhor,nlev), &
                        cwcu(nhor,nlev),gpot(nhor,nlev)
!cgj300611
    !gj
    real(kind=realkind) radf(nhor)
    !gj
    !gj090605
    real(kind=realkind) frice(nhor)
    real(kind=realkind) ps(nhor)
    real(kind=realkind) snowh(nhor),pblh(nhor),texture(nhor), &
    wfrop(nhor),vegopl(nhor)
    !gj090605
    real(kind=realkind) accsunny(nhor)
    real(kind=realkind) dtphysh


    !     workspace to subroutine 'radia' is allocated here.
    !lr   this might not be necessary in present-day computer systems...

    real(kind=realkind) slwup(nhor),tmpcov(nhor),fstopc(nhor), &
         trans(nhor),covmax(nhor),sumemi(nhor), &
         sumcwi(nhor),scos03(nhor), &
         scosre(nhor),wslwnet(nhor), wsswnet(nhor), wsalb(nhor)
    real(kind=realkind) cloabs(nhor)
    real(kind=realkind) emc1(nhor,nlev),fice(nhor,nlev)
    real(kind=realkind) re(nhor,nlev),resw(nhor,nlev),reice(nhor,nlev), &
         req(nhor,nlev),rewat(nhor,nlev)
    integer ibot(nhor)
    real(kind=realkind) bb(nhor,nlev+1),qpath(nhor,nlev+1),    &
         qbpath(nhor,nlev+1),emiqu(nhor,nlev+1), &
         emiqb(nhor,nlev+1),emtotu(nhor,nlev), &
         emtotb(nhor,nlev),cwin(nhor,nlev), &
         efcov(nhor,nlev),efcovu(nhor,nlev+1), &
         dfrac(nhor,nlev), &
         qepath(nhor,nlev+1),qcpath(nhor,nlev+1),  &
         econt(nhor),eco2(nhor),lwdncl(nhor)
    integer ima(nhor),imau(nhor), icovu(nhor,nlev+1)
    real(kind=realkind) dtsw(nhor,nlev+1),dtlw(nhor,nlev+1)
    real(kind=realkind) swnet(nhor,nlev+1),lwnet(nhor,nlev+1)
    real(kind=realkind) emskin(nhor),zwin(nhor)

    !     call subroutine 'radia'

    call radia(  nhor  ,  nlev , kstart,  kstop , hourc ,            &
         minc  ,  secc  ,   locw  ,                              &
         hybf  ,  hybh  , emc   , csusa , fabso3,  cadd  ,             &
         frland,  along , coslat,  sinlat,                             &
         albedo,  tsk   ,  slwup , tmpcov,  scos  , scos03,            &
         scosre,  salbcor, wsalb,  fstopc,  cloabs,                    &
         trans , covmax, ima   ,   imau  , sumemi,  sumcwi,            &
         t     , q     , cw    ,   pf    , ph    ,  dpf   , totcov,    &
         bb    , qpath , qbpath,   emiqu , emiqb ,  emtotu, emtotb,    &
         cwin  , efcov , efcovu,   icovu , dfrac ,  dtdt  ,            &
!         strad , stdrad, stclo ,   stdclo,                             &
         accsunny, dtphysh,                                            &
         emc1,fice,re,resw,reice,req,rewat,ibot,                       &
         wslwnet,wsswnet, slwdn ,   sswdn , tlwnet, tswnet, tswdn,     &
         lwnet  ,  swnet,  dtlw ,   dtsw  ,                            &
         qepath, qcpath, econt ,  eco2  , lwdncl,                      &
         emskin,radf,frice,ps,snowh,				       &
!cgj300611
         wfrop,vegopl,  &
         texture,gpot,pblh,lscov,cucov,cwcu,zwin)
!cgj300611

    return
  end subroutine aradia

  subroutine radia(   nhor   ,  nlev  , kstart,  kstop , hourc , &
       minc   ,  secc  ,   locw   ,                         &
       hybf  , hybh  , emc    ,  csusa , fabso3,  cadd  ,         &
       frland, along , coslat ,  sinlat,                          &
       albedo, tsk    ,  slwup , tmpcov , scos  ,  scos03,        &
       scosre, salbcor, wsalb ,  fstopc,  cloabs,                 &
       trans , covmax, ima    ,  imau  ,  sumemi, sumcwi,         &
       t     , q     ,  cw    ,  pf    ,  ph    , dpf   ,  totcov,&
       bb    , qpath ,  qbpath,  emiqu , emiqb  , emtotu,  emtotb,&
       cwin  , efcov ,  efcovu,  icovu , dfrac  , dtdt  ,         &
!       strad , stdrad,  stclo ,  stdclo,                          &
       accsunny, dtphysh,                                         &
       emc1,fice,re,resw,reice,req,rewat,ibot,                    &
       slwnet, sswnet, slwdn ,  sswdn , tlwnet, tswnet, tswdn,    &
       lwnet,  swnet,  dtlw ,   dtsw  ,                           &
       qepath, qcpath, econt , eco2   , lwdncl,                   &
       emskin,radf,frice,ps,snowh,				  &
!cgj300611
       wfrop,vegopl,  &
       texture,gpot,pblh,lscov,cucov,cwcu,csur4)
!cgj300611

 
    
    implicit none
    !     
    integer,intent(in):: nhor,nlev,kstart,kstop,hourc,minc,secc
    integer::nlevp1,nlevm1,jl,jk,jkp1,jkm1,itab,itabr
    !     
    logical,intent(in):: locw
    !     
    real(kind=realkind) hybh(nlev+1),hybf(nlev), &!strad(nlev), &
         !stdrad(nlev),stclo(nlev),stdclo(nlev),&
         emc(nlev),csusa(nlev),scos03(nhor)


    !     
    real(kind=realkind) fabso3(nlev),cadd(nlev), salbcor(nhor),wsalb(nhor)
    real(kind=realkind) zgrrcp,zrgrav,zcdpr,z10log,zemc,zeps,                   &
         zutc,zhrang,zcos,zpath,dtcle,ztrc,zuzc,zemi,             &
         zusq,zrads,ztrad,ztdrad,ztclo,ztdclo,zdtdt,zcl,zt,zabc,  &
         zefcov,zw1,zw2,zw3,zw4,zw5,zw6,zw7,zfls,zq,zabs1,zabs2,  &
         zq1,zq2,zq3,zq4,zf1,zf2,zf3,zf4,zdif1,zdif3,zdif4,       &
         zdtcle,ztrs,zrad,zsur,ztran                               
    real(kind=realkind) zrpath,zrbpath,zw8,zw9,zetyp,zetmax,                     &
         eline,eh2o,a,b,x,x0,cx,t3,u,ue,uc,                       &
         zalbr,ztrf,zabf,zapr,zfl,zarfl,zden1,zden2                
    real(kind=realkind) frland(nhor),along(nhor),                                &
         coslat(nhor),sinlat(nhor),albedo(nhor),tsk(nhor),        &
         slwup(nhor),tmpcov(nhor),scos(nhor),                     &
         covmax(nhor),sumemi(nhor),sumcwi(nhor),                  &
         trans(nhor),fstopc(nhor)
    real(kind=realkind) cloabs(nhor)

    real(kind=realkind) t(nhor,nlev),q(nhor,nlev),cw(nhor,nlev),               &
         pf(nhor,nlev),ph(nhor,nlev+1),dpf(nhor,nlev),           &
         totcov(nhor,nlev),dtdt(nhor,nlev),bb(nhor,nlev+1),      &
         qpath(nhor,nlev+1),qbpath(nhor,nlev+1),                 &
         emiqu(nhor,nlev+1),emiqb(nhor,nlev+1),                  &
         cwin(nhor,nlev),dfrac(nhor,nlev),                       &
         efcov(nhor,nlev),efcovu(nhor,nlev+1),                   &
         emtotu(nhor,nlev),emtotb(nhor,nlev),                    &
         qepath(nhor,nlev+1),qcpath(nhor,nlev+1),                &
         econt(nhor),eco2(nhor),scosre(nhor),lwdncl(nhor)         
    integer ima(nhor),imau(nhor),icovu(nhor,nlev+1)               
    real(kind=realkind) slwnet(nhor),sswnet(nhor),slwdn(nhor),sswdn(nhor),      &
         tlwnet(nhor),tswnet(nhor),tswdn(nhor)                    
    real(kind=realkind) dtsw(nhor,nlev+1),dtlw(nhor,nlev+1)                      
    real(kind=realkind) swnet(nhor,nlev+1), lwnet(nhor,nlev+1)                   
    real(kind=realkind) emc1(nhor,nlev),fice(nhor,nlev)                          
    real(kind=realkind) re(nhor,nlev),resw(nhor,nlev),reice(nhor,nlev),         &
         req(nhor,nlev),rewat(nhor,nlev)
!cgj300611
    real(kind=realkind) lscov(nhor,nlev),cucov(nhor,nlev), &
                        cwcu(nhor,nlev),gpot(nhor,nlev)
!cgj300611
    integer ibot(nhor)
    !     gj
    real(kind=realkind) radf(nhor)
    real(kind=realkind) zucons
    real(kind=realkind) emskin(nhor),texture(nhor), &
    wfrop(nhor),vegopl(nhor)
    !     gj
    real(kind=realkind) accsunny(nhor),pblh(nhor),zinhm
    real(kind=realkind) dtphysh

    real(kind=realkind) uwcov2d(nhor),uwzcovu(nhor,nlev+1), uwzcovb(nhor,nlev+1)

    !     gj210404
    real(kind=realkind) swo3,swh2o,swsc,swaer,zhum1,zhum2,caahum2
    !     gj210404
    !     gj090605
    real(kind=realkind) ps(nhor),frice(nhor),snowh(nhor),zfrice,cretot
    !     gj090605
    !     gjcsur1
    real(kind=realkind) csur4(nhor),zwin(nhor),ztem(nhor),zwt,zwgas
    !     ---------------------------------------------------------
    real(kind=realkind)  emice,emwat,frmod,tc
    integer  jintt
    real(kind=realkind)      wei, weisum, reqsum, rewsum
    real(kind=realkind)  tau75,sumtau(1:nhor,0:nlev)
!cgj300611
    real(kind=realkind) zcweps,zcwls,zcwcu,zwrk1,zwrk2,zwrk3,zcweff
!cgj040711
    real(kind=realkind) remnice ,demnice,deice,zinhmfac,zzm

    !     interface:
    !     
    !     subroutine "radia" is called from subroutine "aradia" called
    !     from "phys".

    !     variable        type         content

    !     
    !     l    nhor            input        number of horizontal grid points pro-
    !     l                                 cessed at a time.
    !     l    nlev            input        number of vertical levels
    !     l    kstart          input        starting index for horizontal loops
    !     l    kstop           input        ending index for horizontal loops
    !     l    hourc           input        current hour (utc)
    !     l    minc            input        current minut (utc)
    !     l    secc            input        current second (utc)
    !     l    locw            input        switch for cloud water availabilty
    !     l                                 condensation.
    !     l    hybf            input        hybrid levels (full)
    !     l    hybh            input        hybrid levels (half)
    !     l    frland          input        fraction of land in grid square
    !     l    along           input        geographical longitude.
    !     l    coslat          input        cosine of geographical latitude.
    !     l    sinlat          input        sine of geographical latitude.
    !     l    albedo          input        grid-square albedo without solar height 
    !     l                                 correction.
    !     l    tsk             input        grid-square averaged surface temperature (k)
    !     l    t               input        layer temperature (k)
    !     l    q               input        layer specific humidity (kg/kg)
    !     l    cw              input        layer specific cloud water.(kg/kg)
    !     l    pf              input        pressure (pa) at model full levels.
    !     l    ph              input        pressure (pa) at model half levels.
    !     l    dpf             input        pressure difference between model
    !     l                                 half levels.
    !     l    totcov          input        fractional cloud cover for each
    !     l                                 model layer
    !     l    emc             input        cloud emissivity factor
    !     l    csusa           input        array used to estimate layer
    !     l                                 integrated cloud water
    !     pr   fabso3          input        vertical distribution of o3 solar absorption
    !     pr   cadd            input        tuning factor for longwave cooling (now
    !     pr                                defined as a function of height)
    !     l    ----------------------------------------------------------------
    !     l    dtdt            output       temperature tendency (k/s) as a
    !     l                                 result of radiative processes.
    !     l    strad           output       sum of radiative temperature ten-
    !     l                                 dencies for each layer.
    !     l                                 ( for statistics )
    !     l    stdrad          output       sum of squares of radiative tempe-
    !     l                                 rature tendencies for each level.
    !     l                                 ( for statistics )
    !     l    stclo           output       sum of cloud covers for every layer.
    !     l                                 ( for statistics )
    !     l    stdclo          output       sum of squares of cloud covers
    !     l                                 for each model layer.
    !     l                                 ( for statistics )
    !     l    sswnet          output       net surface shortwave radiation
    !     l    slwnet          output       net surface longwave radiation
    !     l    sswdn           output       downwelling sfc shortwave radiation
    !     l    slwdn           output       downwelling sfc longwave radiation
    !     l    tswnet          output       net t.o.a shortwave radiation
    !     l    tlwnet          output       net t.o.a longwave radiation
    !     l    tswdn           output       downwelling t.o.a shortwave radiation
    !     l    salbcor         output       solar height correction for albedo
    !     l    ----------------------------------------------------------------
    !     l    swnet           work array   net radiation flux density at
    !     l                                 model half levels. (w/m2)
    !     l    lwnet           work array   net radiation flux density at
    !     l                                 model half levels. (w/m2)
    !     l    slwup           work array   grid average upward 
    !     l                                 longwave  radiation from the ground
    !     l    tmpcov          work array   temporary cloud cover.
    !     l    scos            work array   cosine of solar zenith angle.
    !     l    scos03          work array   scos**0.3
    !     pr   scosre          work array   true cosine of solar zenith angle
    !     pr                                (allowed to be smaller than 'cosmin')
    !     l    fstopc          work array   shortwave downward flux at the top
    !     l                                 of the uppermost cloud layer.
    !     l    trans           work array   transmission (dimensionless)
    !     l                                 for solar radiation.
    !     l    cloabs          work array   cloud absorption (dimensionless)
    !     l                                 for solar radiation.
    !     l    covmax          work array   maximum fractional cloud cover
    !     l                                 in a vertical column.
    !     l                                 half levels.
    !     l    bb              work array   black body radiation at model half
    !     l                                 levels
    !     l    qpath           work array   scaled water vapour path from top
    !     l                                 of the atmosphere to the top of
    !     l                                 layer "jk"
    !     l    qbpath          work array   scaled water vapour from bottom of
    !     l                                 the atmosphere to the top of layer
    !     l                                 "jk"
    !     l    emiqu           work array   clear air emissivity from top of
    !     l                                 the atmosphere to the top of layer
    !     l                                 "jk"
    !     l    emiqb           work array   clear air emissivity corresponding
    !     l                                 to water path "qbpath"
    !     l    emtotu          work array   cloud emissivity over several model
    !     l                                 layers above layer "jk"
    !     l    emtotb          work array   cloud emissivity over several model
    !     l                                 layers below layer "jk"
    !     l    cwin            work array   vertical integral of cloud water
    !     l                                 ( g/m2 ) in every layer.
    !     l    dfrac          work array    dimensionless "fractional emissi-
    !     l                                 vity" for every model layer.
    !     l    efcov           work array   effective cloud cover for every
    !     l                                 layer.
    !     l    efcovu          work array   effective total cloud cover above
    !     l                                 layer "jk"
    !     l    sumemi          work array   vertical 'integral' of grid box
    !     l                                 cloud water amount.
    !     l    sumcwi          work array   cloud water (g/m2) in layer scaled
    !     l                                 by factor "totcov(jl,jk)/covmax(jl)"
    !     l                                 to be valid for cloud cover "covmax"
    !     l    ima             work array   layer number with maximum effective
    !     l                                 cloud cover "efcov(jl,ima(jl))"
    !     l                                 below layer number "jk".
    !     l    imau            work array   layer number with maximum effective
    !     l                                 cloud cover "efcovu(jl,imau(jl))"
    !     l                                 above layer number "jk"
    !     lrkw emc1,fice,re,resw,reice,req,rewat,ibot
    !     lrkw                 work arrays  for effective radius calculation
    !     pr   qepath          work array   scaled path length for h2o e-type continuum
    !     pr                                absorption from toa to the top of layer "jk"
    !     pr   qcpath          work array   scaled path length for co2 absorption
    !     pr                                from toa to the top of layer "jk"
    !     pr   econt           work array   emissivity for h2o e-type continuum 
    !     pr   eco2            work array   emissivity for carbon dioxide
    !     gjlr lwdncl          work array   lw flux below cloud
    !     l
    !     l    -----------------------------------------------------------------
    !     l    the contents of common block "comrad" with constants and arrays
    !     l    used specifically in the radiation scheme are defined in subrou-
    !     l    tine "inirad".
    !     l

    !     l    authors of subroutine "radia":

    !     l    hannu savijaervi, university of helsinki, october 1991.
    !     l    modified for operational use by
    !     l    bent h.sass,      the danish meteorological institute,
    !     l                      november 1992 and january 1993
    !     l    bent h. sass and laura rontu (fmi) ,
    !     l                      revision february 1994.
    !     l    carl fortelius (fmi), klaus wyser (misu), petri räisänen (dmuh),
    !     l         laura rontu (fmi), colin jones (smhi)
    !     l         revisions december 1998, january 2000, october 2000, june 2002


    !     l    purpose:

    !     l
    !     l       1)   compute temperature tendencies "dtdt(jl,jk)"
    !     l            due to radiative processes (in k/s).
    !     l       2)   compute radiative fluxes at
    !     l            model half levels ( w/m2 ) including sfc and toa
    !     l


    !     l    method:

    !     l            see  j.appl.meteor. june 1990 vol. 29 p437-p447.
    !     l            for effective radius: contr atm phys, 1999,72,205-218
    !     l            for longwave modifications: report 49, 
    !     l                department of meteorology, university of helsinki


    !     l    definition of local constants.
    zeps=0.001_realkind
    zgrrcp=gravit/cpair
    zrgrav=1._realkind/gravit
    zcdpr=cask/cpsref
    zalbr=0.70_realkind
    ztrf=40._realkind
    zabf=0.013_realkind
    zapr=0.8_realkind
    zetmax=3.47e-5_realkind
    nlevp1=nlev+1
    nlevm1=nlev-1
!cgj300611
!    zinhmfac=0.3_realkind/20._realkind
    zinhm=0.6
    zcweps=1.e-10_realkind
    zwrk1=0._realkind
    zwrk2=0._realkind
!cgj300611
!cgj030711
!cgj040711    remnice=10._realkind
    remnice=4._realkind
    if(maximum_random)then

       do jl=kstart,kstop
          uwzcovu(jl,1)=totcov(jl,1)
          uwzcovb(jl,nlev+1)=totcov(jl,nlev)
          uwzcovb(jl,nlev)=totcov(jl,nlev)
       enddo
       do jl=kstart,kstop
          uwcov2d(jl)=1._realkind
       enddo
       do jk=2,nlev
          do jl=kstart,kstop
             uwcov2d(jl)=uwcov2d(jl)* &
                  (1._realkind-max(totcov(jl,jk-1),totcov(jl,jk)) )/ &
                  (1._realkind-min(totcov(jl,jk-1),.99_realkind) )
             uwzcovu(jl,jk)=max(0._realkind,min(1._realkind,1._realkind-uwcov2d(jl)))
          enddo
       enddo
       do jl=kstart,kstop
          uwzcovu(jl,nlev+1)=uwzcovu(jl,nlev)
       enddo
       !     gj
       do jl=kstart,kstop
          uwcov2d(jl)=1._realkind
       enddo
       do jk=nlev-1,1,-1
          do jl=kstart,kstop
             uwcov2d(jl)=uwcov2d(jl)* &
                  (1._realkind-max(totcov(jl,jk+1),totcov(jl,jk)) )/ &
                  (1._realkind-min(totcov(jl,jk+1),.99_realkind) )
             uwzcovb(jl,jk)=max(0._realkind,min(1._realkind,1._realkind-uwcov2d(jl)))
          enddo
       enddo
    else
       do jl=kstart,kstop
          uwzcovu(jl,1)=totcov(jl,1)
          uwzcovb(jl,nlev)=totcov(jl,nlev)
       enddo
       do jk=2,nlev
          do jl=kstart,kstop
             uwzcovu(jl,jk)=max(totcov(jl,jk),uwzcovu(jl,jk-1))
          enddo
       enddo
       uwzcovu(:,nlev+1)=uwzcovu(:,nlev)
       do jk=nlev-1,1,-1
          do jl=kstart,kstop
             uwzcovb(jl,jk)=max(totcov(jl,jk),uwzcovb(jl,jk+1))
          enddo
       enddo
       uwzcovb(:,nlev+1)=uwzcovb(:,nlev)
    endif
    !     
    !     longwave radiation

    !     
    !     -----------------------------------------------------------------
    !     "slwup(jl)" is the area averaged infrared
    !     radiative flux from surface. blackbody radiation at model
    !     half levels is "bb(jl,jk)". at top jk=1, at bottom jk=nlevp1.
    !     -----------------------------------------------------------------
    !     
    do jl=kstart,kstop
       !     

       !     lrtest this is to test consistent formulation of surface emissivity
       !     lrtest in radia and surface scheme. slwup now contains upward
       !     lrtest lw flux emitted by the surface plus reflected downwelling lw
       !     lrtest flux. downwelling lw flux is taken from previous time step 
       !     lrtest (zero at kstep 1, initialized in phys). 
       zt=tsk(jl)*tsk(jl)*tsk(jl)*tsk(jl)
       slwup(jl)=stebol*emskin(jl)*zt+(1._realkind-emskin(jl))*slwdn(jl)
       !     gj        slwup(jl)=stebol*emsurf*zt+(1.-emsurf)*slwdn(jl)
       !     lrtest        slwup(jl)=stebol*zt
       bb(jl,1)=stebol*t(jl,1)*t(jl,1)*t(jl,1)*t(jl,1)

       !     bb(jl,nlevp1) is used in this scheme as the effective emission temperature 
       !     of the lowest layer. in order not to overestimate "heating from
       !     ground" it should be closer to the surface temperature than the layer
       !     mid-point temperature is. the semi-empirical approximation below is based
       !     on an assumption of a logarithmical temperature profile in the lowest half-
       !     layer with a roughness length of a few cm.
       !     
       zt= (2._realkind*t(jl,nlev) + tsk(jl))/3._realkind
       !     gjcsur1
       ztem(jl)=zt
       !     gj050405        zt= (t(jl,nlev) + 2.*tsk(jl))/3.
       bb(jl,nlevp1)=stebol*zt*zt*zt*zt
    enddo
    !     
    !     -----------------------------------------------------------------
    !     security check for negative humidity.
    !     -----------------------------------------------------------------
    !     
    do jk=1,nlev
       do jl=kstart,kstop
          !     
          q(jl,jk)=max(q(jl,jk),0._realkind)
          !     
       enddo
    enddo
    !     -----------------------------------------------------------------

    !     atmospheric blackbody radiation "bb(jl,jk)" at model half levels.
    !     -----------------------------------------------------------------
    do jk=2,nlev
       jkm1=jk-1
       do jl=kstart,kstop
          zt=t(jl,jkm1) +(t(jl,jk)-t(jl,jkm1))* &
               dpf(jl,jkm1)/(dpf(jl,jk)+dpf(jl,jkm1)+zeps)
          bb(jl,jk)=stebol*zt*zt*zt*zt
          !     
       enddo
    enddo
    !     ------------------------------------------------------------------

    !     pr  the modification of lw clear-sky emissivities starts here (pr 990329)


    !     pr initializations
    do jl=kstart,kstop
       qpath(jl,1) = 0._realkind
       qepath(jl,1) = 0._realkind
       qcpath(jl,1) = 0._realkind
       qbpath(jl,nlevp1)=0._realkind
       emiqu(jl,1) = 0._realkind
       emiqb(jl,nlevp1) =0._realkind
    enddo

    do jk =2,nlevp1
       jkm1=jk-1
       do jl=kstart,kstop

          !     pr compute effective absorber amounts from the toa downwards
          !     pr qpath(jl,jk), qepath(jl,jk), and qcpath(jl.jk) are the path lengths for 
          !     pr h2o line absorption, h2o e-type continuum absorption and co2 absorption
          !     pr from the top of the atmosphere to the top of layer 'jk'.

          cx = pf(jl,jkm1)*dpf(jl,jkm1)*cpat
          t3 = t(jl,jkm1)*t(jl,jkm1)*t(jl,jkm1)
          zpath = q(jl,jkm1)*cx

          qpath(jl,jk) = qpath(jl,jkm1) + zpath
          qepath(jl,jk)= qepath(jl,jkm1)+ zpath &
               *1.6078_realkind*q(jl,jkm1)/(1._realkind+0.6078_realkind*q(jl,jkm1)) &
               *6.0202e+14_realkind/(t3*t3+zeps)
          qcpath(jl,jk)= qcpath(jl,jkm1)+ cqco2* cx *t3/1.3824e+7_realkind

          !     pr computation of individual layer emissivity for the the calculation of
          !     pr 'dfrac(jl,jkm1)'. 

          x = sqrt(sqrt(sqrt(zpath)))
          zemc=(x*x*x*x*(ceml1+x*(ceml2+x*(ceml3+x*ceml4))) &
               +cemadd) / (dpf(jl,jkm1)+zeps)

          dfrac(jl,jkm1) = cfrac1 / (cfrac1 +cfrac2*zemc) 
       enddo


       !     pr computation of upward clear air emissivities (emiqu) from level jk to toa


       !     pr h2o e-type continuum 
       do jl=kstart,kstop
          ue= qepath(jl,jk)
          if (ue<cpate) then
             econt(jl) = ue*(ceme1a+ue*ceme2a)
          else
             x = sqrt(ue)
             econt(jl) = ceme1b+x* &
                  (ceme2b+x*(ceme3b+x*(ceme4b+x*ceme5b)))
          end if
       enddo
       !     pr co2
       do jl=kstart,kstop
          uc= qcpath(jl,jk)
          if (uc<cpatc) then
             x = sqrt(uc)
             eco2(jl) = x*(cemc1a+x*(cemc2a+x*cemc3a))
          else
             x = log(uc)
             eco2(jl) = cemc1b+x*(cemc2b+x*(cemc3b+x*cemc4b))
          end if
       enddo

       do jl=kstart,kstop
          !     pr h2o lines and p-type continuum
          u = qpath(jl,jk)
          x0 = sqrt(sqrt(u))
          x = sqrt(x0)
          eline=x0*x0*(ceml1+x*(ceml2+x*(ceml3+x*ceml4)))
          !     pr a and b are the "overlap factors" between water vapor line and continuum
          !     pr absorption, and between water vapor and co2 absorption, respectively 
          a=1._realkind+x0*(ca1+x0*(ca2+x0*(ca3+x0*(ca4+x0*ca5))))
          b=1._realkind+x0*(cb1+x0*(cb2+x0*(cb3+x0*(cb4+x0*cb5))))
          !     pr total emissivity for water vapor
          eh2o= eline +a*(1._realkind-eline)*econt(jl)
          !     pr total clear air emissivity
          emiqu(jl,jk) = eh2o+b*(1._realkind-eh2o)*eco2(jl)
       enddo
    enddo


    !     pr computation of "downward" clear air emissivities (emiqb) between level jk
    !     pr and the surface


    do jk=1,nlev
       !     pr h2o e-type continuum
       do jl=kstart,kstop
          ue= qepath(jl,nlevp1)-qepath(jl,jk)
          if (ue<cpate) then
             econt(jl) = ue*(ceme1a+ue*ceme2a)
          else
             x = sqrt(ue)
             econt(jl) = ceme1b+x*(ceme2b+x*(ceme3b+x*(ceme4b+x*ceme5b)))
          end if
       enddo
       !     pr co2
       do jl=kstart,kstop
          uc= qcpath(jl,nlevp1)-qcpath(jl,jk)
          if (uc<cpatc) then
             x = sqrt(uc)
             eco2(jl) = x*(cemc1a+x*(cemc2a+x*cemc3a))
          else
             x = log(uc)
             eco2(jl) = cemc1b+x*(cemc2b+x*(cemc3b+x*cemc4b))
          end if
       enddo
       !     pr h2o lines and p-type continuum
       do jl = kstart,kstop
          u = qpath(jl,nlevp1)-qpath(jl,jk)
          !     pr qbpath (the water vapor path length from the surface upwards) is still
          !     pr needed for shortwave calculations
          qbpath(jl,jk) = u                 
          x0= sqrt(sqrt(u))
          x = sqrt(x0)
          eline=x0*x0*(ceml1+x*(ceml2+x*(ceml3+x*ceml4)))
          !     pr a and b are the "overlap factors" between water vapor line and continuum
          !     pr absorption, and between water vapor and co2 absorption, respectively 
          a=1._realkind+x0*(ca1+x0*(ca2+x0*(ca3+x0*(ca4+x0*ca5))))
          b=1._realkind+x0*(cb1+x0*(cb2+x0*(cb3+x0*(cb4+x0*cb5))))
          !     pr total emissivity for water vapor
          eh2o= eline +a*(1._realkind-eline)*econt(jl)
          !     pr total clear air emissivity
          emiqb(jl,jk) = eh2o+b*(1._realkind-eh2o)*eco2(jl)
          !     gjcsur
          zwin(jl)=b*(1._realkind-eh2o)
       enddo
    enddo
    !     gjcsur
    zwt=0.004_realkind
    do jl=kstart,kstop
!!       if(abs(texture(jl)-2.0)<1.e-14_realkind.and.wfrop(jl)>0.7_realkind.and.vegopl(jl)<0.2_realkind)then
       if( (abs(texture(jl)-2.0_realkind)<1.e-14_realkind) .or. &
           (abs(texture(jl)-6.0_realkind)<1.e-14_realkind) .or.&
           (abs(texture(jl)-9.0_realkind)<1.e-14_realkind) )then
        csur4(jl)=csur1+11.0_realkind
       else
        csur4(jl)=csur1
       endif
       zw1=1._realkind-((305._realkind-ztem(jl))*zwt)
       zw2=max(0.625_realkind,min(1._realkind,zw1))
       zwgas=zwin(jl)*zwin(jl)*1.9_realkind*zw2
       zwgas=max(1._realkind,min(1.65_realkind,zwgas))
       csur4(jl)=csur4(jl)*zwgas
    enddo
    !     gjcsur

    !     pr  the modification of lw clear-sky emissivities ends here (pr 990329)

    !     ------------------------------------------------------------------

    !     estimate vertically integrated cloud water content "cwin(jl,jk)"
    !     ( g/m**2 ), the integral inside clouds, not grid box value. .
    !     compute emissivity "zemc" of clouds and a reduced effective
    !     fractional cloud cover "efcov(jl,jk)"
    !     ------------------------------------------------------------------
    !     
    !     initialize for layer 1, assuming there are no clouds
    !     
    do jl=kstart,kstop
       re(jl,1)=0._realkind
       cwin(jl,1)=0._realkind
       efcov(jl,1)=0._realkind
    enddo
    !     
    call ficefun(nhor,nlev,kstart,kstop,t,fice)

    do jk=1,nlev
       do jl=kstart,kstop
          if( locw ) then
!             if (totcov(jl,jk)>zeps) then
             if (totcov(jl,jk)>zeps .and. cw(jl,jk)>zcweps) then
!cgj300611
!cgj300611 Average diagnosed cwcu (convective cloud water) and ls cloud water as in
!cgj300611 Tiedtke 1996
!cgj030911 zinhm is reduced in the cold/stable boundary layer.....less inhomogeneity...test with tke later
!cgj300611
!                zinhm(jl)=0.6_realkind
!                zzm=gpot(jl,jk)*zrgrav
!                if(zzm<=pblh(jl))then
!                  if(tsk(jl)<273.16_realkind .and. zzm<1750._realkind)then
!                    zwrk1=max(0._realkind,min(0.3_realkind,   &
!                             ((273.16_realkind-tsk(jl))*zinhmfac)))
!                    zinhm(jl)=zinhm(jl)+zwrk1
!                  endif
!                else
!                   if(pf(jl,jk)<20000._realkind .and.      &
!                                t(jl,jk)<230._realkind)then
!                     zwrk1=max(0._realkind,min(0.3_realkind,   &
!                             ((230._realkind-t(jl,jk))*zinhmfac)))
!                     zinhm(jl)=zinhm(jl)+zwrk1
!                   endif
!                endif
!
                zcwls=max(0._realkind,min(cw(jl,jk),(cw(jl,jk)-cwcu(jl,jk))))
                zcwcu=cw(jl,jk)-zcwls
                if(zcwcu>zcweps)zwrk1=log10(zinhm*zcwcu)
                if(zcwls>zcweps)zwrk2=log10(zinhm*zcwls)
                zwrk3=(cucov(jl,jk)/totcov(jl,jk))*zwrk1 + &
                      (lscov(jl,jk)/totcov(jl,jk))*zwrk2
                zcweff=10._realkind**zwrk3
!cgj300611
                zw1=zcweff/(totcov(jl,jk)+zeps)
!                zw1=cw(jl,jk)/(totcov(jl,jk)+zeps)
             else
                zw1=0._realkind
             endif
          else
             zw1=csusa(jk)*q(jl,jk)
          endif
          cwin(jl,jk)=zrgrav*ckilgr*dpf(jl,jk)*zw1
          !     gj040405
          !     gj040405  apply the subgrid scale tau variability to both sw and lw 
!cgj300611 subgrid scale term done with Tiedtke avergaing above
!cgj300611          cwin(jl,jk)=0.8_realkind*cwin(jl,jk)
!cgj300611 Increase cloud inhomogeneity factor
!          cwin(jl,jk)=0.65_realkind*cwin(jl,jk)
          !     gj040405
!cgj300611
!cgj300611 Average diagnosed cwcu (convective cloud water) and ls cloud water as in
!cgj300611 Tiedtke 1996
!cgj300611
          
          !     
          !     calculate the effective radius for water and ice
          !     relax the droplet number towards the background value at
          !     sigma=0.7, the background value is identical with the sea value
          !--------------------------------------------
          !     
          !     gj090605
          zfrice=frice(jl)*(1._realkind-frland(jl))
          cretot=creland+(cresea-creland)*min(1._realkind,max(0._realkind,&
               1._realkind-(pf(jl,jk)-0.8_realkind*ps(jl))/(0.2_realkind*ps(jl))))
!cgjrca4               1._realkind-(pf(jl,jk)-0.7_realkind*ps(jl))/(0.3_realkind*ps(jl))))
          !ramp between polluted air over land and clean over sea
          cretot=cretot+(cresea-cretot)*min(1._realkind,max(0._realkind,1._realkind-frland(jl)))
          !modify for snow depth over land
!cgj050911         increase in cretot and rewat acts to reduce cloud albedo by around 2-4% over snow at low zenith angles
!cgj050911          cretot=cretot+(cresice-cretot)*min(1._realkind,max(0._realkind,snowh(jl)))
          !ramp between resulting value and sea-ice value
          cretot=cretot+(cresice-cretot)*min(1._realkind,max(0._realkind,zfrice))
          !     gj090605
          frmod=frland(jl)*max(0._realkind,(hybf(jk)-.7_realkind)/.3_realkind)
          !     gj        rewat(jl,jk)=max(remn, cmmy*(creland*frmod+cresea*(1.-frmod))*
          rewat(jl,jk)=max(remn, cmmy*cretot* &
               (zw1*pf(jl,jk)/(rair*t(jl,jk)))**(1.0_realkind/3.0_realkind))
          !     gj090605
          !     
          !     re for ice is function of t, after ou & liou,1993
          tc=min(-20._realkind,t(jl,jk)-273._realkind)
!cgj030711          reice(jl,jk)=max(remn,0.5_realkind*(((.0012_realkind*tc+.197_realkind)*tc+12.42_realkind)*tc+326.3_realkind))
          !cgj150811reice(jl,jk)=max(remnice,0.5_realkind*(((.0012_realkind*tc+.197_realkind)*tc+12.42_realkind)*tc+326.3_realkind))
!cgj150811
          demnice=2.*remnice
          deice=max(demnice,((.0012_realkind*tc+.197_realkind)*tc+12.42_realkind)*tc+326.3_realkind)
          reice(jl,jk)=max(remnice,(((1.2601e-7_realkind*deice- &
               3.0954e-5_realkind)*deice+5.6416e-3_realkind)*&
               deice+0.56383_realkind)*deice-2.2054_realkind)
          reice(jl,jk)=min(125._realkind,reice(jl,jk))
!cgj150811
          !     gj200705
          !     gj210705        reice(jl,jk)=1.5*(((.0012*tc+.197)*tc+12.42)*tc+326.3)
          !     gj210705        reice(jl,jk)=min(120.,reice(jl,jk))
          !     gj200705
          !     gj
          !     gj        jintt=max(1,1+nint((t(jl,jk)-200.)/dttabl))
          !     gj        jintt=min(jintt,740)
          !     gj        reice(jl,jk)=rehm(jintt)
          !     
          !     emissivity is a function of effective radius, different for
          !     water and ice.
          !     
          emwat=.0255_realkind+.2855_realkind*exp(-.0890_realkind*rewat(jl,jk))
          emice=.0202_realkind+.2059_realkind*exp(-.0676_realkind*reice(jl,jk))
          emc1(jl,jk)=(1._realkind-fice(jl,jk))*emwat+fice(jl,jk)*emice
          zw1=emc1(jl,jk)*cwin(jl,jk)
          !     
          !     gjorig          if (zw1<0.5) then
          !     gjorig             zemc=zw1

          !     !1.-max(1.-zw1,0.)
          !     gjorig          else if (zw1<5.) then
          if(zw1<5._realkind)then
             zemc=1._realkind-exp(-zw1)
          else
             zemc=1._realkind
          endif
          !     
          efcov(jl,jk)=zemc*totcov(jl,jk)
          !     
       enddo
    enddo
    !     
    !     ---------------------------------------------------------------

    !     determine maximum effective cloud cover "efcovu(jl,jk)"
    !     valid for layers above "jk" , and the associated layer
    !     number "icovu(jl,jk)".
    !     ---------------------------------------------------------------
    !     
    do jl=kstart,kstop
       !     
       imau(jl)=1
       sumemi(jl)=0._realkind
       tmpcov(jl)=0._realkind
       !     
    enddo
    !     
    do jk=1,nlev
       do jl=kstart,kstop
          !     
          !sumemi(jl)=sumemi(jl) + totcov(jl,jk)*cwin(jl,jk)
          !     
          !zw1=emc1(jl,jk)*sumemi(jl)
          !     kw080506 reduce emissivity from above by absorption in layer jk
          !     kw080506 and add emissivity from layer jk
          zw1=max(0._realkind,sumemi(jl)*(1._realkind-emc1(jl,jk))+ &
               emc1(jl,jk)*totcov(jl,jk)*cwin(jl,jk))
          sumemi(jl)=sumemi(jl)+emc1(jl,jk)*totcov(jl,jk)*cwin(jl,jk)

          !     gjorig          if (zw1<0.5) then
          !     gjorig             emtotu(jl,jk)=zw1

          !     1.-max(1.-zw1,0.)
          !     gjorig          else if (zw1<5.) then
          if(zw1<5._realkind)then
             emtotu(jl,jk)=1._realkind-exp(-zw1)
          else
             emtotu(jl,jk)=1._realkind
          endif
          !     
          efcovu(jl,jk)=efcov(jl,1)
          icovu(jl,jk)=imau(jl)
          !     
       enddo
    enddo
    !     
    do jk=2,nlev
       jkm1=jk-1
       do jl=kstart,kstop
          !     
          efcovu(jl,jk)=tmpcov(jl)
          icovu(jl,jk)=imau(jl)
          if(.not.maximum_random)then
             zemi=max( totcov(jl,jkm1)*emtotu(jl,jkm1),efcov(jl,jkm1) )
          else
             zemi=uwzcovu(jl,jkm1)*emtotu(jl,jkm1)
          endif
          if( zemi>=tmpcov(jl) ) then
             efcovu(jl,jk)=zemi
             icovu(jl,jk)=jkm1
             tmpcov(jl)=zemi
             imau(jl)=jkm1
          endif
          !     
       enddo
    enddo
    !     
    do jl=kstart,kstop
       !     
       efcovu(jl,nlevp1)=tmpcov(jl)
       icovu(jl,nlevp1)=imau(jl)
       if(.not.maximum_random)then
          zemi=max( totcov(jl,nlev)*emtotu(jl,nlev),efcov(jl,nlev) )
       else
          zemi=uwzcovu(jl,nlev)*emtotu(jl,nlev)
       endif
       if( zemi>=tmpcov(jl) ) then
          efcovu(jl,nlevp1)=zemi
          icovu(jl,nlevp1)=nlev
       endif
       !     
    enddo
    !     
    !     --------------------------------------------------------------

    !     estimate a total cloud emissivity "emtotb" from the bottom of
    !     the atmosphere to the top of layer "jk"
    !     ---------------------------------------------------------------
    !     
    do jl=kstart,kstop
       !     
       sumemi(jl)=0._realkind
       ima(jl)=nlev
       tmpcov(jl)=0._realkind
       !     
    enddo
    !     
    do jk=nlev,1,-1
       do jl=kstart,kstop
          !     
          !sumemi(jl)=sumemi(jl) + totcov(jl,jk)*cwin(jl,jk)
          !zw1=emc1(jl,jk)*sumemi(jl)
          !     kw080506 reduce emissivity from below by absorption in layer jk
          !     kw080506 and add emissivity from layer jk
          zw1=max(0._realkind,sumemi(jl)*(1._realkind-emc1(jl,jk))+&
               emc1(jl,jk)*totcov(jl,jk)*cwin(jl,jk))
          sumemi(jl)=sumemi(jl)+emc1(jl,jk)*totcov(jl,jk)*cwin(jl,jk)
          !     gjorig          if (zw1<0.5) then
          !     gjorig             emtotb(jl,jk)=zw1

          !     1.-max(1.-zw1,0.)
          !     gjorig          else if (zw1<5.) then
          if(zw1<5._realkind)then
             emtotb(jl,jk)=1._realkind-exp(-zw1)
          else
             emtotb(jl,jk)=1._realkind
          endif
          !     
       enddo
    enddo
    !     
    !     initialize longwave fluxes
    !---------------------------
    do jl=kstart,kstop
       slwdn(jl)=0._realkind
       slwnet(jl)=0._realkind
       tlwnet(jl)=0._realkind
    enddo

    do jk=1,nlev+1
       do jl=kstart,kstop
          dtlw(jl,jk)=0._realkind
          lwnet(jl,jk)=0._realkind
       enddo
    enddo

    !     compute infrared heating effect of all layers except the upper-
    !     most and the lowest one. the heating rate consists of the fol-
    !     lowing terms "zq1,zq2,zq3,zq4" described below.

    do jk=nlevm1,2,-1
       jkm1=jk-1
       jkp1=jk+1
       do jl=kstart,kstop
          zefcov=tmpcov(jl)
          if(.not.maximum_random)then
             zemi=max( totcov(jl,jkp1)*emtotb(jl,jkp1),efcov(jl,jkp1) )
          else
             zemi=uwzcovb(jl,jkp1)*emtotb(jl,jkp1)
          endif
          if( zemi>=tmpcov(jl) ) then
             ima(jl)=jkp1
             zefcov=zemi
             tmpcov(jl)=zemi
          endif
          !     
          !     ------------------------------------------------------------

          !     "zq1" is the flux difference across the clear fraction of
          !     the layer with no overlying clouds.
          !     ------------------------------------------------------------
          zw1=1._realkind-max( efcov(jl,jk),efcovu(jl,jk) )
          zt=0.5_realkind*(bb(jl,jk) +bb(jl,jkp1))
          zfls=(slwup(jl)-zt)*(emiqb(jl,jk)-emiqb(jl,jkp1))
          !     pr990329: the empirical correction for e-type continuum is no longer needed!
          !     pr         zetyp=cetyp*q(jl,jk)*q(jl,jk)*q(jl,jk)
          !     pr         zetyp=min( zetyp, zetmax )
          !     pr         zq=zt*(emiqu(jl,jk)-emiqu(jl,jkp1))
          !     pr  +      -(zetyp +cadd)*dpf(jl,jk)/zgrrcp
          zq=zt*(emiqu(jl,jk)-emiqu(jl,jkp1))-cadd(jk)*dpf(jl,jk)/zgrrcp
          zq1=zw1*( zq +zfls )
          !     ------------------------------------------------------------

          !     "zq2" is the flux difference across the clear fraction of
          !     the layer with overlying clouds.
          !     ------------------------------------------------------------
          zw2=max( efcovu(jl,jk) -efcov(jl,jk),0._realkind )
          zw3=min( max( zefcov-efcov(jl,jk), 0._realkind), zw2 )
          zf2=min( (ph(jl,jk) -ph(jl,icovu(jl,jk)+1))/cpc, 1._realkind )
          zq2=zw2*zq*zf2 + (zw2-zw3)*zfls
          !     ------------------------------------------------------------
          !     "zq3" is the flux difference across the cloudy fraction of
          !     the layer with no overlying clouds.
          !     ------------------------------------------------------------
          zw4=max( efcov(jl,jk)-efcovu(jl,jk), 0._realkind )
          zw5=min( efcov(jl,jkp1),efcov(jl,jk) )
          zw6=max( zw5-efcovu(jl,jk), 0._realkind )
          zw7=min( zefcov,efcov(jl,jk) )
          zw8=max( zw7-efcovu(jl,jk),0._realkind )
          !     
          zf3=min( cpf/(ph(jl,ima(jl))-ph(jl,jkp1)+zeps),1._realkind )
          zdif3=bb(jl,ima(jl)) -bb(jl,jkp1)
          zf4=min( cpf/(ph(jl,nlevp1)-ph(jl,jkp1)+zeps),1._realkind )
          zdif4=slwup(jl) -bb(jl,jkp1)
          !     
          !     pr990329: the empirical correction for e-type continuum is no longer needed!
          !     pr            zq3=zw4*( -bb(jl,jk) +0.5*(bb(jl,jk) +bb(jl,jkm1))*
          !     pr     &          emiqu(jl,jk) +csur1 +zcsur4*sqrt( q(jl,jkm1) ) +
          !     pr     &                               csur2*q(jl,jkm1) -
          !     pr     &          zcs1*min( (1.-hybh(jk))/(1.-chybcr) ,1. ) )
          zq3=zw4*( -bb(jl,jk) +0.5_realkind*(bb(jl,jk) +bb(jl,jkm1))* &
          !cgjcsur4     emiqu(jl,jk) +csur1*hybh(jk))
               emiqu(jl,jk) +csur4(jl)*hybh(jk))
          zq3=zq3 + (zw8-zw6)*zdif3*(zf3+dfrac(jl,jkp1)*(1._realkind-zf3))
          zq3=zq3 + (zw4-zw8)*zdif4*(zf4+dfrac(jl,jkp1)*(1._realkind-zf4))
          !     ------------------------------------------------------------

          !     "zq4" is the flux difference across the cloudy fraction of
          !     the layer with overlying clouds.
          !     ------------------------------------------------------------
          zw6=min( efcovu(jl,jk), efcov(jl,jk) )
          zw7=min(min( zefcov,efcov(jl,jk) ),efcovu(jl,jk) )
          zw8=max( zw6 -efcov(jl,jkm1) ,0._realkind)
          zw9=max( zw7 -efcov(jl,jkp1), 0._realkind )
          zf1=min( cpf/(ph(jl,jk)-ph(jl,icovu(jl,jk)+1)+zeps),1._realkind )
          zdif1=bb(jl,icovu(jl,jk)+1) -bb(jl,jk)
          !     
          zq4=zw8*zdif1*( zf1 +dfrac(jl,jkm1)*(1._realkind-zf1) )
          zq4=zq4 + zw9*zdif3*( zf3 +dfrac(jl,jkp1)*(1._realkind-zf3) )
          zq4=zq4 + (zw6-zw7)*zdif4*( zf4+dfrac(jl,jkp1)*(1._realkind-zf4) )
          !     
          !     -----------------------------------------------------------

          !     infrared heating rate "dtlw(jl,jk) in layer "jk" due to
          !     all terms zq1, zq2, zq3 and zq4
          !     ---------------------------------------------------------
          !     
          dtlw(jl,jk)=zgrrcp*( zq1 +zq2 + zq3 + zq4 )/dpf(jl,jk)
          !     
       enddo
    enddo
    !     
    !     ---------------------------------------------------------------
    !     temperature tendencies of uppermost layer and of lowest model
    !     layer.
    !     ---------------------------------------------------------------
    !     
    do jl=kstart,kstop
       !     
       !     pr990329 assume pure "cooling to space" in the uppermost layer, instead of
       !     pr        setting the temperature tendency to zero
       !     pr     dtdt(jl,1)=0.
       dtlw(jl,1)=0.5_realkind*(bb(jl,1)+bb(jl,2))*(emiqu(jl,1)-emiqu(jl,2)) &
            *zgrrcp/(dpf(jl,1)+zeps)
       !     -------------------------------------------------------------
       zw1=1._realkind-max( efcov(jl,nlev),efcovu(jl,nlev) )
       zt=bb(jl,nlevp1)
       zfls=(slwup(jl)-zt)*(emiqb(jl,nlev)-emiqb(jl,nlevp1))
       !     pr990329: the empirical correction for e-type continuum is no longer needed!
       !     pr        zetyp=cetyp*q(jl,nlev)*q(jl,nlev)*q(jl,nlev)
       !     pr        zetyp=min( zetyp, zetmax )
       !     pr        zq=zt*(emiqu(jl,nlev)-emiqu(jl,nlevp1))
       !     pr     &  -(zetyp + cadd)*dpf(jl,nlev)/zgrrcp
       zq=zt*(emiqu(jl,nlev)-emiqu(jl,nlevp1)) &
            -cadd(nlev)*dpf(jl,nlev)/zgrrcp
       zq1=zw1*( zq +zfls )
       !     -------------------------------------------------------------
       zw2=max( efcovu(jl,nlev) -efcov(jl,nlev), 0._realkind )
       zf2=min( (ph(jl,nlev) -ph(jl,icovu(jl,nlev)+1))/cpc, 1._realkind )
       zq2=zw2*( zq*zf2 + zfls )
       !     -------------------------------------------------------------
       zw4=max( efcov(jl,nlev) -efcovu(jl,nlev), 0._realkind )
       !     
       !     pr990329: the empirical correction for e-type continuum is no longer needed!
       zq3=zw4*( -bb(jl,nlev) +0.5_realkind*( bb(jl,nlev) +bb(jl,nlevm1))* &
       !cgj     emiqu(jl,nlev) +csur1*hybh(nlev) -bb(jl,nlevp1) +slwup(jl))
            emiqu(jl,nlev) +csur4(jl)*hybh(nlev) - &
            bb(jl,nlevp1) +slwup(jl))
       !     
       !     -------------------------------------------------------------
       zw6=min( efcovu(jl,nlev), efcov(jl,nlev) )
       zw9=max( zw6 -efcov(jl,nlevm1), 0._realkind )
       zf1=min( cpf/(ph(jl,nlev)-ph(jl,icovu(jl,nlev)+1)+zeps),1._realkind )
       zdif1=bb(jl,icovu(jl,nlev)+1) -bb(jl,nlev)
       !     
       zq4=zw9*zdif1*(zf1 +dfrac(jl,nlevm1)*(1._realkind-zf1))
       zq4= zq4 + zw6*( slwup(jl) -bb(jl,nlevp1) )
       !     
       !     -------------------------------------------------------------
       !     
       dtlw(jl,nlev)=zgrrcp*( zq1 + zq2 + zq3 + zq4 )/dpf(jl,nlev)
       !     
    enddo
    !     
    !     ---------------------------------------------------------------

    !     determination of surface downwelling longwave radiation "slwdn(jl)"
    !     at first a cloud free contribution is calculated
    !     in a loop over all model layers.the final net radiation is
    !     determined by adding a cloudy contribution  with weight
    !     "efcovu(jl,nlevp1)".
    !     ---------------------------------------------------------------
    !     
    do jl=kstart,kstop
       lwdncl(jl)=0._realkind
    enddo
    !     
    !     pr990410-begin
    !     to be consistent with the computation of "heating from ground"
    !     in the lowest layer, compute the contribution of the lowest
    !     layer to the downward flux at the surface separately

    do jl=kstart,kstop
       slwdn(jl) = bb(jl,nlevp1)*(emiqb(jl,nlev)-emiqb(jl,nlevp1))
       lwdncl(jl) = slwdn(jl)
    enddo
    !     pr990410-end
    !     pr990410 'nlev' changed to 'nlev-1'
    do jk=nlev-1,1,-1
       jkp1=jk+1
       do jl=kstart,kstop
          !     
          slwdn(jl)=slwdn(jl) + 0.5_realkind*( bb(jl,jk)+bb(jl,jkp1) )* &
               ( emiqb(jl,jk) - emiqb(jl,jkp1) )
          !     
          !     gjlr------------------------------
          !     gjlr calculate clear sky emission below cloud base
          !     gjlr
          if(jk>=(icovu(jl,nlevp1)))then
             lwdncl(jl)=lwdncl(jl)+0.5_realkind*( bb(jl,jk)+bb(jl,jkp1) )* &
                  ( emiqb(jl,jk) - emiqb(jl,jkp1) )
          endif
          !     gjlr------------------------------
       enddo
    enddo
    !     
    do jl=kstart,kstop
       !     
       !     pr990329: the empirical correction for e-type continuum is no longer needed!
       !     pr         slwdn(jl)=slwdn(jl) + csur1 + zcsur4*sqrt( q(jl,nlev) ) +
       !     pr     &                       csur2*q(jl,nlev)
       !     gj300305
       !     gj300305        slwdn(jl)=slwdn(jl) + csur1
       !     gj300305
       !     
       !     gjlr------------------------------
       !     gjlr caluclate clear sky contribution to lw flux at surface
       !     gjlr it is a combination of clear sky flux form
       !     gjlr depth of atmosphere in clear portion of grid box plus
       !     gjlr clear contribution from below cloud base, multiplied
       !     gjlr by cloud fraction plus csur1 (missing gas)
       !     gjlr contribution, which is assumed to be a full depth of atmos
       !     gjlr contribution (for want of knowing how better to partition it)
       !     gjlr
       !     gjcsur1        slwdn(jl)=( (1.-efcovu(jl,nlevp1))*(slwdn(jl)+csur1) )+
       slwdn(jl)=( (1._realkind-efcovu(jl,nlevp1))*(slwdn(jl)+csur4(jl)) )+ &
       !cgjcsur4     ( efcovu(jl,nlevp1)*(lwdncl(jl)+csur1) )
            ( efcovu(jl,nlevp1)*(lwdncl(jl)+csur4(jl)) )
       !     gjlr------------------------------
       zsur=t(jl,nlev) -cdpsur*(t(jl,nlev) -t(jl,nlevm1))/&
            (pf(jl,nlev) -pf(jl,nlevm1)+zeps)
       zrad=stebol*zsur*zsur*zsur*zsur
       !     gj150404
       zrad=bb(jl,nlevp1)
       !     gj150404
       !     gj300305
       !     gj300305        zemi=min( slwdn(jl)/zrad, 1. )
       !     gj300305
       !     gjlr------------------------------
       !     gjlr zemi normalised cloudy flux contribution as a function
       !     gjlr clear sky emissivity. this is now calculated from cloud base
       !     gjlr to surface in lwdncl array
       !     gjlr
       !     gj200705
       zemi=min(lwdncl(jl)/zrad, 1._realkind )
       !     gj200705
       !     gjorig        zemi=min( lwdncl(jl)/zrad, 1. )
       !     gjlr------------------------------
       slwdn(jl)=slwdn(jl) +csur3*(1._realkind-zemi)*efcovu(jl,nlevp1)* &
            bb(jl,icovu(jl,nlevp1)+1)
       !     gj110605     &           (0.5*(bb(jl,icovu(jl,nlevp1)+1)+
       !     gj110605     &                 bb(jl,icovu(jl,nlevp1))))
       !     gj110605
       slwnet(jl)=slwdn(jl)-slwup(jl)
       !     lrtest        slwnet(jl)=emsurf(jl)*(slwdn(jl)-slwup(jl))
       !     lrtest use the form with emsurf if slwup=sigma*tskin**4 
       !     lrtest the form without emsurf is identical to this in case
       !     lrtest slwup(jl)=stebol*emsurf(jl)*zt+(1.-emsurf(jl))*slwdn(jl)
       !     lrtest both forms correspond the definitions of subroutine rad2surf
       !     
    enddo

    !     


    !     shortwave radiation

    !     
    !     calculate local solar zenith angle:
    !     -----------------------------------
    !15 degrees/hour 15/60 = 0.25 degrees/min 
    zutc = real(hourc,realkind)*15._realkind+ &
         real(minc,realkind)*0.25_realkind + real(secc,realkind)*15._realkind/3600._realkind - &
         180._realkind

    do jl=kstart,kstop
       zhrang=(zutc+along(jl))*cdegr
       zcos=coslat(jl)*cosd*cos(zhrang)+sinlat(jl)*sind
       
       scosre(jl) = zcos
       zcos=max( zcos, cosmin )
       scos(jl)=zcos
       scos03(jl)=zcos**0.3_realkind
            
       fstopc(jl)=0._realkind
       sswdn(jl)=0._realkind
       sswnet(jl)=0._realkind
       tswnet(jl)=0._realkind
       tswdn(jl)=0._realkind
       do jk=1,nlev+1
          swnet(jl,jk)=0._realkind
          dtsw(jl,jk)=0._realkind
       enddo

    enddo

    !     
    !     initialize horizontal arrays
    !     ----------------------------
    !     
    do jl=kstart,kstop
       !     
       covmax(jl)=0._realkind
       sumcwi(jl)=0._realkind
       sumtau(jl,0)=0._realkind
       trans(jl)=1._realkind
       cloabs(jl)=0._realkind
       ima(jl)=nlev+1
       ibot(jl)=nlev
       !     
    enddo
    !     
    !     solar angle correction to surface albedo (slightly modified by pr 990327)
    !     -------------------------------------------------------------------------
    do jl = kstart,kstop
       salbcor(jl)=0.2_realkind/(1._realkind+scos(jl)) - 0.12_realkind
       wsalb(jl)=albedo(jl)+salbcor(jl)
    enddo
    !     ------------------------------------------------------

    !     determine maximum cloud cover "covmax" in the vertical
    !     column.
    !     ------------------------------------------------------
    !     
    !     kw080507 use uwzcovu that has been defined for maximum and maximum/random
    !     kw080507 overlap at the beginning of radia
    do jl=kstart,kstop
       covmax(jl)=uwzcovu(jl,nlev)

       !increase 2-d cloud cover at very high solar angles
       covmax(jl)=min(1._realkind,covmax(jl)*(1._realkind+.03_realkind*(1._realkind/scos(jl)-1._realkind)))
    enddo

    do jk=1,nlev
       do jl=kstart,kstop
          if(totcov(jl,jk)>zeps) then
             if (ima(jl)>nlev) ima(jl)=jk
             ibot(jl)=jk
          endif
          !     
          !add tau, assume cw in g/m3
          !use overlap assumption as in sw
          sumtau(jl,jk)=sumtau(jl,jk-1)+                          &
               1.5e3_realkind*cwin(jl,jk)*totcov(jl,jk)/(covmax(jl)+zeps)*  &
               ((1._realkind-fice(jl,jk))/(1000._realkind*rewat(jl,jk))+             &
               fice(jl,jk)/(917._realkind*reice(jl,jk)))
       enddo
    enddo
    !     
    !     

    !     compute downward solar flux "fstopc(jl)" at the top of
    !     the uppermost cloud layer.
    !     ------------------------------------------------------
    !     
    do jl=kstart,kstop
       !     
       if( ima(jl)<=nlev .and. ima(jl)>1 ) then
          !     
          !     lrkw define zalbr from cloud properties
          !     lrkw calculation copied here from trans/abs calculations below
          !     lrkw
          tau75=sumtau(jl,nlev)**.75_realkind
          ztrc=5.0_realkind+7.0_realkind*scos(jl)
          ztrc=ztrc/(ztrc+3.0_realkind*tau75)*exp(-.017_realkind*tau75)
          zalbr=(1._realkind-ztrc)*zapr


          zuzc=sqrt( qpath(jl,ima(jl))/scos(jl))
          zusq=sqrt( zuzc )
          !     gj030204
          if(qpath(jl,ima(jl))<0.01_realkind)then
             zucons=0.155_realkind
          elseif(qpath(jl,ima(jl))>=0.01_realkind .and.qpath(jl,ima(jl))<=1._realkind)then
             zucons=0.155_realkind-(0.02525_realkind*(qpath(jl,ima(jl))-0.01_realkind))
          else
             zucons=0.13_realkind
          endif
          !     gj030204
          !     gj210404          fstopc(jl)=sa*scos(jl)*(1.-0.024/sqrt(scos(jl)) -
          !     gj210404     &               0.11*caak*zusq -zcdpr*pf(jl,ima(jl)-1)*
          !     gj210404     &               ( 0.28/(1.+6.43*scos(jl)) -0.07*zalbr ) )
          !     gj210404
          !     gj210404   follow more exactly lacis &  hansen (jas 1974) upon 
          !     gj210404   which savijärvi(1990) the sw sfc flux calculation
          !     gj210404
!cgj060911 Increase aerosols over desert region.....mimics sand in air
!ps111117 Included open land fraction (wfrop) and vegetation cover (vegopl) in test
!cgj          if(abs(texture(jl)-2.0)<1.e-14_realkind.and.wfrop(jl)>0.7_realkind.and.vegopl(jl)<0.2_realkind)then
          if( (abs(texture(jl)-2.0_realkind)<1.e-14_realkind) .or. &
              (abs(texture(jl)-6.0_realkind)<1.e-14_realkind) .or. &
              (abs(texture(jl)-9.0_realkind)<1.e-14_realkind) )then
           zcdpr=(cask+0.05)/cpsref
           caak=1.30
          else
           zcdpr=cask/cpsref
           caak=1.25
          endif
!cgj060911
          swo3=0.024_realkind/sqrt(scos(jl))
          !     gjorig          swh2o=0.11*caak*zusq
          swh2o=zucons*caak*zusq
          swsc=0.28_realkind/(1._realkind+6.43_realkind*scos(jl))
          swaer=zcdpr*pf(jl,ima(jl)-1)
          !     gj210404
          !     kw040524          fstopc=sa*scos(jl)*
          fstopc(jl)=sa*scos(jl)*( (0.353_realkind-swh2o) + &
               ( (0.647_realkind-swo3-swaer*swsc)/(1._realkind-0.07_realkind*wsalb(jl)) ) )
          !     gj210404
          fstopc(jl)=max(0._realkind,fstopc(jl))
       else
          fstopc(jl)=sa*scosre(jl)
          fstopc(jl)=max(0._realkind,fstopc(jl))
       endif
       !     
    enddo
    !     
    !     solar radiative heating by h2o-absorption
    !     (gas+cloud droplets);
    !     start main loops:
    !     ------------------------------------------------------
    !     
    !     $dir no_peel
    do jk=1,nlev
       jkp1=jk+1
       do jl=kstart,kstop
          !     gj220404
          !     gj220404  let caahum2 decrease from 0.05 to 0.03(caahum) between pf(nlev)
          !     gj220404  and 0.8*pf(nlev)
          zhum1=1._realkind-(pf(jl,jk)/pf(jl,nlev))
          zhum2=min(0.2_realkind,zhum1)
          caahum2=0.04_realkind-(0.05_realkind*zhum2)
          !     gj220404

          !     
          !     ----------------------------------------------------------

          !     clear sky part as in savijaervi (j. appl. meteor. 1990)
          !     ( correction by bhs,i.e. multiplication of "dtcle"
          !     by scos(jl) to produce correct effect per unit area.
          !     ----------------------------------------------------------
          zpath=0.5_realkind*( qpath(jl,jkp1) + qpath(jl,jk) )
          !     lr
          zrbpath=qpath(jl,nlevp1)/scos(jl)
          zrpath=0.5_realkind*crefrad*( qbpath(jl,jkp1) + qbpath(jl,jk))+ zrbpath
          !     lr
          !     
          if( zpath<cpatsh ) then
             !     --------------------------
             !     lr
             !     lr in dry case exponents for incoming radiation and 'shtab' for
             !     lr reflected radiation
             !     lr
             !     au0610111/k-i               zw1=caadry*(zpath/scos(jl))**(-cbbdry)
             if ( abs(zpath)< 1.e-14_realkind ) then
                zw1=0._realkind
             else
                zw1=caadry*(zpath/scos(jl))**(-cbbdry)
             endif
             itabr=nint( 200._realkind*zrpath )
             itabr=max( min(itabr,ntabm),1 )
             dtcle=sa*csac*q(jl,jk)*pf(jl,jk)* &
                  (zw1+caahum*albedo(jl)*crefrad*scos(jl)*shtab(itabr)) &
                  +1.7e-06_realkind*scos03(jl)
          else

             !     lr
             !     lr in humid case 'shtab' for incoming and reflected radiation
             !     lr
             itab=nint( 200._realkind*zpath/scos(jl) )
             itab=max(min( itab,ntabm ),1 )
             itabr=nint( 200._realkind*zrpath )
             itabr=max( min(itabr,ntabm),1 )
             dtcle=sa*csac*q(jl,jk)*pf(jl,jk)*caahum2* &
                  (shtab(itab)+albedo(jl)*crefrad*scos(jl)*shtab(itabr)) &
                  +1.7e-06_realkind*scos03(jl)
          endif
          !     
          !     pr990329 add stratospheric ozone absorption to clear-sky sw heating
          dtcle = dtcle+0.024_realkind*sa*sqrt(scos(jl))*zgrrcp &
               *fabso3(jk)/(dpf(jl,jk)+zeps)
          !     
          !     ----------------------------------------------------------

          !     determine shortwave radiative temperature tendency "dtsw". 
          !     cloud transmission and absorbtion from slingo, liou, et al.
          !     ----------------------------------------------------------
          !     
          tau75=sumtau(jl,jk-1)**.75_realkind
          ztrc=1._realkind+2._realkind*scos(jl)
          cloabs(jl)=(.15_realkind+.08_realkind*scos(jl))*(1._realkind-ztrc/(tau75+ztrc)) &
               *(1._realkind-exp(-.18_realkind*tau75))

          tau75=sumtau(jl,jk)**.75_realkind
          zabs2=(.15_realkind+.08_realkind*scos(jl))*(1._realkind-ztrc/(tau75+ztrc)) &
               *(1._realkind-exp(-.18_realkind*tau75))

          zw1=max(0._realkind,zabs2-cloabs(jl))

          !     lr----------------------------------------------------------------------

          zdtcle=dtcle*covmax(jl)*trans(jl)
          !     
          if (scosre(jl)<=0._realkind) then
             dtsw(jl,jk)=0._realkind
          else
             dtsw(jl,jk)=zgrrcp/(dpf(jl,jk)+zeps)* &
                  covmax(jl)*fstopc(jl)*zw1*trans(jl)+ &
                  zdtcle + dtcle*(1._realkind-covmax(jl))
          endif
          !     
          ztran=trans(jl)
          ztrc=5.0_realkind+7.0_realkind*scos(jl)
          trans(jl)=ztrc/(ztrc+3.0_realkind*tau75)*exp(-.017_realkind*tau75)
          zabc=ztran -trans(jl)
          zabc=max( zabc, ztran*dtcle*dpf(jl,jk)/ &
               (zgrrcp*fstopc(jl)+zeps) )
          trans(jl)=min( 1.0_realkind, max( ztran-zabc, 0.0_realkind ) )
          trans(jl)=min(trans(jl),ztran)
          !     
          !     ----------------------------------------------------------

          !     net surface shortwave radiative flux density at the ground
          !     ----------------------------------------------------------
          if( jk==nlev ) then
             zuzc=sqrt( qpath(jl,nlevp1)/scos(jl) )
             zusq=sqrt( zuzc )
             !     gj030204
             if(qpath(jl,nlevp1)<0.01_realkind)then
                zucons=0.155_realkind
             elseif(qpath(jl,nlevp1)>=0.01_realkind .and.qpath(jl,nlevp1)<=1._realkind)then
                zucons=0.155_realkind-(0.02525_realkind*(qpath(jl,nlevp1)-0.01_realkind))
             else
                zucons=0.13_realkind
             endif
             !     gj061103
             !     gj210404               ztrs=sa*scos(jl)*(1. -0.024/sqrt(scos(jl)) -
             !     gj210404     &             zucons*caak*zusq -zcdpr*pf(jl,nlev)*
             !     gj210404     &             ( 0.28/(1. +6.43*scos(jl)) -0.07*wsalb(jl) ) )
             !     gj210404
             !     gj210404   follow more exactly lacis &  hansen (jas 1974) upon 
             !     gj210404   which savijärvi(1990) the sw sfc flux calculation
             !     gj210404
!cgj060911 Increase aerosols over desert region.....mimics sand in air
!ps111117 Included open land fraction (wfrop) and vegetation cover (vegopl) in test
!cgj          if(abs(texture(jl)-2.0)<1.e-14_realkind.and.wfrop(jl)>0.7_realkind.and.vegopl(jl)<0.2_realkind)then
          if( (abs(texture(jl)-2.0_realkind)<1.e-14_realkind) .or. &
              (abs(texture(jl)-6.0_realkind)<1.e-14_realkind) .or. &
              (abs(texture(jl)-9.0_realkind)<1.e-14_realkind) )then
           zcdpr=(cask+0.05)/cpsref
           caak=1.30
          else
           zcdpr=cask/cpsref
           caak=1.25
          endif
!cgj060911
             swo3=0.024_realkind/sqrt(scos(jl))
             swh2o=zucons*caak*zusq
             swsc=0.28_realkind/(1._realkind+6.43_realkind*scos(jl))
             swaer=zcdpr*pf(jl,nlev)
             !     gj210404
             ztrs=sa*scos(jl)*( (0.353_realkind-swh2o) + &
                  ( (0.647_realkind-swo3-swaer*swsc)/(1._realkind-0.07_realkind*wsalb(jl)) ) )
             !     gj210404
             !     gj061103
             !     gj020304
             !     gj130404              ztrs=sa*scos(jl)*(1. -0.024/sqrt(scos(jl)) -
             !     gj130404     &             0.11*caak*zusq -zcdpr*pf(jl,nlev)*
             !     gj130404     &             ( 0.28/(1. +6.43*scos(jl)) -0.07*wsalb(jl) ) )
             !     
             zrads=ztrs*(1._realkind-covmax(jl))+fstopc(jl)*covmax(jl)*trans(jl)/ &
                  (1._realkind-wsalb(jl)*(1._realkind-trans(jl))*zapr)
             !     pr940411 for computational reasons, a minimum value 'scosmin' has been assumed
             !     above for the cosine of solar zenith angle. however, if the sun is
             !     really below the horizon, 'zrads' should be set to zero
             !     if(scosre(jl)<=0.) zrads=0.
             zrads=max(0._realkind,zrads)
             sswdn(jl)=zrads
             sswnet(jl)=(1._realkind-wsalb(jl))*sswdn(jl)
             !     
             !     gj
             radf(jl)=slwnet(jl)+sswnet(jl)
             !     gj
             !     
             !     save sunshine hours
             !     
             !     it is sunny if lowest level short-wave radiation
             !     is larger than 120 w / squaremeter
             !     
             if (zrads > 120._realkind) accsunny(jl)= accsunny(jl)+dtphysh
          endif
          !     
       enddo
    enddo
    !     ----------------------------------------------------------------

    !     bhs940411-begin
    !     include absorption of radiation from reflected beams
    !     in cloudy part of grid box.
    !     ----------------------------------------------------------------
    do jk=1,nlev
       do jl=kstart,kstop
          !     
          zfl=fstopc(jl)*covmax(jl)*trans(jl)
          zden1=wsalb(jl)*(1._realkind-trans(jl))
          zden2=zden1*zapr
          if( jk>ibot(jl) ) then
             zarfl=zgrrcp*zfl*( 1._realkind/(1._realkind-zden1) -1._realkind/(1._realkind-zden2) )/ &
                  (ph(jl,nlevp1)-ph(jl,ibot(jl)+1)+zeps)
          else
             zarfl=0.0_realkind
          endif

          !     pr 990327 the absorption of reflected radiation should also be accounted for
          !     pr 990327 in shortwave heating. 
          !     pr940411 for computational reasons, a minimum value 'scosmin' has been assumed
          !     above for the cosine of solar zenith angle. however, if the sun is
          !     really below the horizon, 'zarfl' should be set to zero

          if (scosre(jl)<=0._realkind) zarfl=0._realkind
          dtsw(jl,jk)=dtsw(jl,jk) + zarfl

          !     lr ---------------------------------------------------
          !     lr total radiative heating as a sum of lw and sw parts
          !     lr ---------------------------------------------------

          dtdt(jl,jk)=dtlw(jl,jk)+dtsw(jl,jk)

       enddo
    enddo
    !     
    !     ----------------------------------------------------------
    !     predicted temperature changes in the uppermost stra-
    !     tosphere ( for "hybf(j) < hybcr" ) are neglected.
    !     lrpr now hybcr = 0, thus these lines do not have any effect
    !     ----------------------------------------------------------
    !     
    !     lr      do jk=1,nlev
    !     lr        do jl=kstart,kstop
    !     lr
    !     lr          if( hybf(jk)<chybcr ) then
    !     lr             dtdt(jl,jk)=0.0
    !     lr             dtsw(jl,jk)=0.0
    !     lr          endif
    !     lr
    !     lr        enddo
    !     lr      enddo
    !     
    !     ----------------------------------------------------------------

    !     determine sw and lw net radiative flux densities (w/m2) 
    !     between layers.
    !     ----------------------------------------------------------------
    !     
    do jl=kstart,kstop
       swnet(jl,nlevp1)=sswnet(jl)
       lwnet(jl,nlevp1)=slwnet(jl)
    enddo

    do jk=nlev,1,-1
       jkp1=jk+1
       do jl=kstart,kstop
          swnet(jl,jk)=swnet(jl,jkp1) + dtsw(jl,jk)*dpf(jl,jk)/zgrrcp
          lwnet(jl,jk)=lwnet(jl,jkp1) + dtlw(jl,jk)*dpf(jl,jk)/zgrrcp
       enddo
    enddo
    !     
    do jl=kstart,kstop
       tswnet(jl)=swnet(jl,1)
       tswdn(jl)=max(0._realkind,sa*scosre(jl))
       tlwnet(jl)=lwnet(jl,1)
    enddo

    return
  end subroutine radia






end module radiation
