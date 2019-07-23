module ccons
  use confys
  use decomp,only:realkind
  implicit none
  private
  real(kind=realkind),public,save::he273= 611.0_realkind !saturation vapor pressure for t=273K
  real(kind=realkind),public,save::tvirtc = 0.61_realkind !constant needded to compute virtual temperature
  real(kind=realkind),public,save::hkevap= 5.e-5_realkind !coefficient for evaporation from stratiform precipitation
  real(kind=realkind),public,save::u00max= 0.975_realkind !maximum allowable value of modified hu00
  real(kind=realkind),public,save::hu00 = 0.85_realkind  !threshold relative humidity for stratiform condensation over
  real(kind=realkind),public,save::tcir = 235.0_realkind  !temperature below which cirrus (pure ice crystal) clouds are considered
  real(kind=realkind),public,save::aecon  !exp(conae)
  real(kind=realkind),public,save::conae = 0.15_realkind  !constant for the development in series of the equivalent potential temperature
  real(kind=realkind),public,save::coales = 100.0_realkind !parameter of the coalescence factor
  real(kind=realkind),public,save::hccu  = 2.e-4_realkind !parameter c00 for cumulus, that is the conversion rate of cloud to precipitation drops in convective clouds
  real(kind=realkind),public,save::hcunrm = 3.e4_realkind!that is cloud cover depens on cloud depth compared to hcunrm
  real(kind=realkind),public,save::hmrcu = 8.e-4_realkind !mr for convective clouds, that is cloud water mixing ratio at which conversion becomes efficient in convective clouds
  real(kind=realkind),public,save::hmrst = 4.e-4_realkind !mr for stratiform clouds, that is cloud water mixing ratio at which conversion becomes efficient in stratiform clouds
  real(kind=realkind),public,save::htaucu  = 3600.0_realkind!characteristic time used in convective cloud cover scheme
  real(kind=realkind),public,save::hp0    = 1.e5_realkind !reference pressure (=100 kpa)
  real(kind=realkind),public,save::cbfeff = 4.0_realkind!parameter for the Bergeron-Findeisen effect
  real(kind=realkind),public,save::hvterm = 5.0_realkind!terminal velocity of precipitation
  real(kind=realkind),public,save::hvsnow = 1.0_realkind!terminal velocity of snow
  real(kind=realkind),public,save::hkmelt = 2.5e-4_realkind!coefficient for melting of ice in precipitation
  real(kind=realkind),public,save::tanvil  = 253.0_realkind!temperature below which convective anvil is allowed and treated as stratiform
  real(kind=realkind),public,save::hcst = 1.e-4_realkind  !conversion rate of cloud to precipitation drops in stratifor clouds
  real(kind=realkind),public,save::pmoist = 3.0_realkind! exponent for computations of hmoist (eq. 3 in OSZ94)
  real(kind=realkind),public,save::prbice(1750) !porbability of ice (eq. 1 in S93)
  real(kind=realkind),public,save::bfeff(1750)  ! first value of the Bergeron-Findeisen effect (without take into account the ice precipitation) (eq. 13 in S93)
  real(kind=realkind),public,save::hdewi(1750)  !difference between both saturation vapor pressure over water over ice
  real(kind=realkind),public,save::hmroft(1750) !temperature function to multiply mr for T<273 (eq. 25 in OSZ
  real(kind=realkind),public,save::dttabl !temperature increment to build the table
  real(kind=realkind),public,save::hkap   !rair/cpair
  real(kind=realkind),public,save::hdldcp !latice/cpair
  real(kind=realkind),public,save::hecdr  !epsilo/hkap
  real(kind=realkind),public,save::hldcp  !latvap/cpair
  real(kind=realkind),public,save::coalcu !coalescence factor for convection clouds (c0.mr/2g)
  real(kind=realkind),public,save::coalst !coalescence factor for convection clouds (c0.mr/2g)
  real(kind=realkind),public,save::sqvsno !square root of the snow terminal velocity
  real(kind=realkind),public,save::heldr  !epsilo*latvap/rair
  real(kind=realkind),public,save::elotci !epsilo*latice/rair/tcir
  real(kind=realkind),public,save::asnow= 9._realkind !a parameter together with bsnow and snoref to govern the
  !           rate of ice precip before the Bergeron-Findeisen mechanism
  !           becomes effective, even if the precip from above is pure
  !           ice;
  !           the factor cbfsno = asnow*(snowrate/snoref)**bsnow
  !           and cbfsno = cbfsno/(1+cbfsno)
  !           which multiplies the modified ice probability
  real(kind=realkind),public,save::bsnow = 4.0_realkind!see asnow
  real(kind=realkind),public,save::snoref = 0.05_realkind!see asnow
  real(kind=realkind),public,save::rehm(740)

  public inikf
contains
  subroutine inikf

    ! inicons - setup of Sundqvist scheme parameters


    !  REFERENCES:
    !
    !  Sundqvist, H. (1993).-"Inclusion of Ice Phase of Hydrometeors in Clou
    !     Parameterization for Mesoscale and Largescale Models". Beit. Phys.
    !     Atmosph., 66,137-147. Later will be S93.
    !
    !  Olofsson, P-O.; Sundqvist, H.; Zurovac-Jevtic, D. (1994).-"Documentat
    !     of the routine condens". MISU. Later will be OSZ94.
    !
    !  some index for loops
    !

    !  local constants to build the tables
    !

    implicit none
    real(kind=realkind)::hedldr,tci,todpmx,tscale,demax,apri,test,tphalf,podpmx, &
         fxt,xti,x,x232,bfmax,xpt1mp

    integer i10,it,i

    !  he273  - saturation vapor pressure for t=273K
    !  tvirtc - constant needded to compute virtual temperature
    !  hkevap - coefficient for evaporation from stratiform precipitation
    !  u00max - maximum allowable value of modified hu00
    !
    he273  = 611.0_realkind
    tvirtc = 0.61_realkind
    hkevap = 5.e-5_realkind
    u00max = 0.975_realkind
    !
    !  hu00   - threshold relative humidity for stratiform condensation over
    !  tcir   - temperature below which cirrus (pure ice crystal) clouds are
    !           considered
    !  conae  - constant for the development in series of the equivalent
    !           potential temperature
    !  aecon  - exp(conae)
    !
    hu00 = 0.85_realkind
    tcir  = 235.0_realkind
    conae = 0.15_realkind
    aecon = exp(conae)
    !
    !  coales - parameter of the coalescence factor
    !  hccu   - parameter c00 for cumulus, that is the conversion rate of
    !           cloud to precipitation drops in convective clouds
    !  hcunrm - that is cloud cover depens on cloud depth compared to hcunrm
    !  hmrcu  - mr for convective clouds, that is cloud water mixing ratio
    !           at which conversion becomes efficient in convective clouds
    !  hmrst  - mr for stratiform clouds, that is cloud water mixing ratio
    !           at which conversion becomes efficient in stratiform clouds
    !  htaucu - characteristic time used in convective cloud cover scheme
    !
    coales = 100.0_realkind
    hccu = 2.e-4_realkind
    hcunrm = 3.e4_realkind
    hmrcu = 8.e-4_realkind
    hmrst = 4.e-4_realkind
    htaucu = 3600.0_realkind
    !
    !  hp0    - reference pressure (=100 kpa)
    !  cbfeff - parameter for the Bergeron-Findeisen effect
    !  hvterm - terminal velocity of precipitation
    !  hvsnow - terminal velocity of snow
    !
    hp0    = 1.e5_realkind
    cbfeff = 4.0_realkind
    hvterm = 5.0_realkind
    hvsnow = 1.0_realkind
    !
    !  asnow  - a parameter together with bsnow and snoref to govern the
    !           rate of ice precip before the Bergeron-Findeisen mechanism
    !           becomes effective, even if the precip from above is pure
    !           ice;
    !           the factor cbfsno = asnow*(snowrate/snoref)**bsnow
    !           and cbfsno = cbfsno/(1+cbfsno)
    !           which multiplies the modified ice probability
    !  bsnow  - see asnow
    !  snoref - see asnow
    !
    asnow = 9.0_realkind
    bsnow = 4.0_realkind
    snoref = 0.05_realkind
    !
    !  hkmelt - coefficient for melting of ice in precipitation
    !  tanvil - temperature below which convective anvil is allowed and
    !           treated as stratiform
    !  hcst   - conversion rate of cloud to precipitation drops in stratifor
    !           clouds
    !  pmoist - exponent for computations of hmoist (eq. 3 in OSZ94)
    !
    hkmelt = 2.5e-4_realkind
    tanvil = 253.0_realkind
    hcst = 1.e-4_realkind
    pmoist = 3.0_realkind
    !
    !  hkap   - rair/cpair
    !  hdldcp - latice/cpair
    !  hecdr  - epsilo / hkap
    !  hldcp  - latvap/cpair
    !
    hkap   = rair/cpair
    hdldcp = latice/cpair
    hecdr  = epsilo / hkap
    hldcp  = latvap/cpair
    !
    !  coalcu - coalescence factor for convection clouds (c0.mr/2g)
    !  coalst - coalescence factor for convection clouds (c0.mr/2g)
    !  sqvsno - square root of the snow terminal velocity
    !  heldr  - epsilo*latvap/rair
    !  hedldr - epsilo*latice/rair
    !  elotci - epsilo*latice/rair/tcir
    !
    coalcu = 0.5_realkind* hccu* hmrcu / gravit
    coalst = 0.5_realkind* hcst* hmrst / gravit
    sqvsno = sqrt(hvsnow)
    heldr  = epsilo*latvap/rair
    hedldr = epsilo*latice/rair
    elotci = epsilo*latice/rair/tcir
    !
    !  dttabl - temperature increment to build the table
    !  prbice - porbability of ice (eq. 1 in S93)
    !  bfeff  - first value of the Bergeron-Findeisen effect (without take
    !           into account the ice precipitation) (eq. 13 in S93)
    !  hdewi  - difference between both saturation vapor pressure over water
    !           over ice
    !  hmroft - temperature function to multiply mr for T<273 (eq. 25 in OSZ
    !
    !
    !             Parameter values in si units
    !             ----------------------------
    !
    tci    = 232.0_realkind
    todpmx = 299.0_realkind
    tscale = (todpmx - tci)*sqrt(2.0_realkind)
    !
    !               the table look up is done by finding the integer
    !               it = 1 + int((t- h273)/dttabl)
    !
    dttabl = 0.1_realkind
    !
    !             calculate saturation vapor pressure to be used as table
    !
    !             calculate difference of saturation pressure between
    !                       water and ice, dewi, for t< t0
    !             calculate a probability function, prbice, for existence of
    !                       ice crystals
    !             calculate a temperature function, hmroft,
    !                       to multiply hmrcu and hmrst
    !
    !  demax  - maximum value of hdewi
    !  apri   - parameter A in S93
    !  test   - 2A/(2A-1)
    !  tphalf - temperature with a probalility of ice equal to 0.5 (T3 in S9
    !  podpmx - ice probability for todpmx
    !

    demax = 0.0_realkind
    apri = 1.0_realkind/(1.0_realkind-exp(-((tmelt-tci)/tscale)**2))
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
    !
    !  table from T=273 up to T=99K
    !
    !  fxt    - temperature
    !  xti    - (1/T0-1/T)
    !  hdewi  - difference between both saturation vapor pressure over water
    !           over ice
    !  prbice - porbability of ice (eq. 1 in S93)
    !  hmroft - temperature function to multiply mr for T<273 (eq. 25 in OSZ
    !
    do 25 i = 1,1750
       fxt = tmelt - (real(i,realkind)-1.0_realkind) * dttabl
       xti = 1.0_realkind/tmelt - 1.0_realkind/fxt
       hdewi(i)  = he273/fxt * exp(heldr*xti) * (1.0_realkind- exp(hedldr*xti))
       demax = max(demax,hdewi(i))

       prbice(i) = apri * (exp(-(((fxt - tci)/tscale)**2)) - 1.0_realkind) + 1.0_realkind
       prbice(i) = max( prbice(i) , 0.0_realkind)
       if(fxt < tci) prbice(i) = 1.0_realkind
       prbice(i) = min(prbice(i) , 1.0_realkind)
       if(fxt > 250.0_realkind) hmroft(i) = 1.33_realkind*exp(-((fxt-tmelt)*0.066_realkind)**2)

       if(fxt <= 250.0_realkind) then
          x = abs (fxt - 232.0_realkind) / 18.0_realkind
          x = x * (1.0_realkind + x * (1.0_realkind + 1.333_realkind * x))
          x232 = 1.0_realkind
          if (fxt < 232.0_realkind) x232 = -1.0_realkind
          x = x / (1.0_realkind + x) * x232
          hmroft(i) = 0.5_realkind*0.15_realkind*(1.07_realkind+x)
       endif
       hmroft(i) = min( hmroft(i) , 1.0_realkind)
       !        hmroft - becomes fmr (eq. 24 in OSZ94)
       hmroft(i) = (1.0_realkind-prbice(i))**2 + prbice(i) * hmroft(i)
       hmroft(i)   = max(hmroft(i) , 3.e-2_realkind)
25  enddo

    !  setting values for i=1

    prbice(1) = 0.0_realkind
    hmroft(1) = 1.0_realkind

    !             normalize with max difference, demax
    !             calculate the product of dewi and prbice, called bfeff

    !  bfmax  - maximum value of bfeff

    bfmax = 0.0_realkind
    !
    !  bfeff  - first value of the Bergeron-Findeisen effect (without take
    !           into account the ice precipitation) (eq. 13 in S93)
    !
    do 30 i = 1,1750
       hdewi(i)  = hdewi(i)/demax
       bfeff( i) = prbice(i) * (1.0_realkind-prbice(i)) * hdewi(i)
       bfmax = max(bfmax,bfeff(i))
30  enddo
    !
    !  normalizing bfeff with hir maximum value
    !
    do 32 i = 1,1750
       bfeff( i) = bfeff( i) / bfmax
32  enddo
    !
    !  print a selection of the computed values
    !
    !MtR

    do 60 i10 = 1,47,3
       it = 1 - i10
       i = 1 + (i10-1) * 10
       xpt1mp = prbice(i)*(1.0_realkind-prbice(i))
60  enddo



    snoref = snoref/3600.0_realkind
    return
  end subroutine inikf


end module ccons
