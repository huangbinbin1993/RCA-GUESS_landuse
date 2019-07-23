module vari
  use mod_grib,only:realgribkind
  implicit none
  private 

  public gb2r2v,gb2gcc,gb2gca
  
contains
  
  subroutine gb2gca(pxold,pyold,pxnew,pynew,pdiri,k, &
       kb2i,pb2i,kb2o,pb2o)
    !
    !     subroutine gb2gca: transform coordinates as per gribblocks 2
    !
    ! input:
    !     pxold,pyold: input coordinates (metres, degrees)
    !     pdiri:       initial rotation angle values
    !     k:           number of coordinate pairs to be transformed
    !     kb2i,pb2i:   integer resp. real values in block 2
    !                  describing input grid (cf grib code fm 92)
    !     kb2o,pb2o:   idem, describing output grid
    !
    ! output:
    !     pxnew,pynew: transformde coordinates
    !     pdiri:       as on input, but incremented by
    !                  rotation angle (degrees) btwn coordinate systems
    !
    ! notes:

    !     1.           pxold may be allocated the same memory words
    !                  as pxnew; this routine then overwrites pxold
    !                  similarly for pyold and pynew.
    !     2.           pdiri is incremented (see notes in gb2rla)
    !                  in degrees, but for efficiency purposes pdiri
    !                  is not reduced to within -180...180. most
    !                  often, pdiri will be used as an argument to
    !                  trigonometric functions, so that the reduction
    !                  is not required at all
    !     3.           the real parameters in the fm-92 grib block 2
    !                  definition are expected in the arrays pb2i/o
    !                  the corresponding integer elements of kb2i/o
    !                  contain pointers to the relevant element of
    !                  pb2i/o.  e.g., for rotated  and stretched
    !                  regular lat/lon grids, we could have:
    !                  kb2(39)=1, to indicate that
    !                  pb2(1) is the angle of rotation, and
    !                  kb2(49)=2, i.e. pb2(kb2(49)) is the
    !                  stretching factor.
    !     4.           implemented are:
    !                  a)regular lat/lon projections,
    !                  perhaps stretched, perhaps rotated
    !                  (coordinates in degrees wrt 0 degrees lat,
    !                  0 degrees lon)
    !                  b)polar stereographic projections
    !                  (coordinates are in metres wrt the grid
    !                  lower left corner, which itself is
    !                  given in 0.001 degrees in the greographical
    !                  lat/lon grid, as words 11 and 14 of grib
    !                  block 2)
    !
    ! method:
    !     1. set a number of parameters to define in -and output grids
    !     2. convert polar-stereographic input to either polar-stereo
    !        graphic output or geographical lat/lon (i.e. pole at 90)
    !     3. destretch stretched input to regular lat/lon with pole
    !        at the pole of stretching
    !     4. de-gauss gaussian input (not implemented)
    !        or de-spectralise spectral grid (not implemented)
    !     5. convert regular lat/lon to regular lat/lon
    !        (for stretched output: pole at pole of stretching)
    !     6. interpolate to gaussian output grid (not implemented)
    !        or spectralise to spectral grid (not implemented)
    !     7. stretch stretched output grid
    !     8. convert regular lat/lon to polar-stereographic output
    !
    !     the meaning of the local variables is as follows:
    !     names ending in 1 or i refer to the input grid,
    !     those ending in 2 or o to the output grid.
    !     for lat/lon grids,
    !     (zlam,zthe) give longitude and latitude of the pole,
    !     zlov the angle of rotation, for polar-stereograpihic grids
    !     zthe is the latitude of the projection plane, and zlov the
    !     grid orientation.
    !     names ending in 1 and 2 define the grids involved in the
    !     central transformation (lat/lon to lat/lon); for lat/lon
    !     grids these are the in- or output grids, and for the
    !     following cases these grids are the geographical
    !     lat/lon projections:
    !     a. polar-stereographic projections,
    !     b. gaussian grids,
    !     c. spectral grids,
    !     while for (d.) stretched grids these grids are the
    !     lat/lon projections, with the northpole rotated to the
    !     pole of stretching.
    !     names ending in i or o define the grids from resp. to which
    !     the grids defined by names ending
    !     in 1 or 2 have to be transformed to get the
    !     in- resp. output grids in resp. from lat/lon.
    !     schematically: i to 1 (1 is regular lat/lon)
    !                    1 to 2 (2 is regular lat/lon)
    !                    2 to o
    !     short-cuts are taken for specific combinations of
    !     in- and output grids
    !
    !     u/v may be specified wrt the grid (called 'matching') or
    !     wrt to true north (called 'true').  the following combinations
    !     are possible, iuvtm is used to indicate which is the case
    !     iuvtm:
    !      0    both in- and output grid match
    !      1    input grid matches, output grid true
    !      2    output grid matches, input grid true
    !      3    both in- and output grid true
    !     these cases are treated as follows:
    !      0: use direction changing routines gb2xxa
    !      1: first transformation to be done is to geographical lat/lon,
    !         using direction changing routines gb2xxa, then continue
    !         as if the input grid is geographical lat/lon and iuvtm=3
    !      2: first transformation to be done is to geographical lat/lon,
    !         using direction non-changing routines gb2xxc, then continue
    !         as if the input grid is geographical lat/lon and iuvtm=0
    !      3: use direction non-changing routines gb2xxc
    !
    !
    !                                 g.j.cats, 05 jan 89.
    !
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    integer kb2i(*),kb2o(*)
    real(kind=realgribkind) pb2i(*),pb2o(*)
    real(kind=realgribkind) zlono(1),zlato(1),zlati(1),zloni(1)
    logical losame
    real(kind=realgribkind) p1,p2,zsame,zthe1,zlam1,zlov1,zthei,zlovi,zlonoi,zlatoi,zlami
    real(kind=realgribkind) zstri
    integer iuvtm,j,k
    losame(p1,p2)=abs(p1-p2)<zsame
    data zsame/0.01_realgribkind/
    logical lopsi,lolali, logeoi,logaui,loroti,lostri,lonew,louvmi
    logical lolalo,logeoo,logauo,lopso,loroto,lostro,locoin,louvmo
    real(kind=realgribkind) zthe2,zlam2,zlov2,ztheo,zlovo,zlonoo,zlatoo,zlamo,zstro




    lonew=.false.

    !     1.2 prepare interpolation parameters for input grid

    lopsi=.false.
    lolali=.false.
    logeoi=.false.
    logaui=.false.
    loroti=.false.
    lostri=.false.
    !     louvmi indicates whether to bother about u/v representation

    louvmi=.true.
    iuvtm=(1-mod(kb2i(17)/8,2))*2

    !     1.2.0 geographical lat/lon
    if(kb2i(6)==0)then
       lolali=.true.
       logeoi=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       louvmi=.false.

       !     1.2.2 polar-stereographic (defined by g.cats, knmi/dm-88-03)
    elseif(kb2i(6)==2)then
       lopsi=.true.
       zthei=real(kb2i(29),realgribkind)
       zlovi=0.01_realgribkind*real(kb2i(31),realgribkind)
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       zlonoi=0._realgribkind
       zlatoi=0._realgribkind

       !     1.2.4 gaussian geographical lat/lon
    elseif(kb2i(6)==4)then
       lolali=.true.
       logeoi=.true.
       logaui=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       louvmi=.false.

       !     1.2.5 polar-stereographic
    elseif(kb2i(6)==5)then
       lopsi=.true.
       zthei=60._realgribkind
       if(kb2i(27)==1)zthei=-zthei
       zlovi=0.001_realgribkind*real(kb2i(18),realgribkind)
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       zlati(1)=0.001_realgribkind*real(kb2i(11),realgribkind)
       zloni(1)=0.001_realgribkind*real(kb2i(14),realgribkind)
       call gb2lpc(zloni,zlati,zlono,zlato,1,zthei,zlovi)
       zlonoi=zlono(1)
       zlatoi=zlato(1)

       !     1.2.10 rotated lat/lon
    elseif(kb2i(6)==10)then
       lolali=.true.
       loroti=.true.
       zthe1=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlov1=pb2i(kb2i(39))

       !     1.2.14 rotated gaussian lat/lon
    elseif(kb2i(6)==14)then
       lolali=.true.
       loroti=.true.
       logaui=.true.
       zthe1=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlov1=pb2i(kb2i(39))

       !     1.2.20 stretched lat/lon
    elseif(kb2i(6)==20)then
       lolali=.true.
       lostri=.true.
       zthei=-90._realgribkind
       zlami=0._realgribkind
       zlovi=0._realgribkind
       zthe1=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(39))
       louvmi=.false.

       !     1.2.24 stretched gaussian lat/lon
    elseif(kb2i(6)==24)then
       lolali=.true.
       lostri=.true.
       logaui=.true.
       zthei=-90._realgribkind
       zlami=0._realgribkind
       zlovi=0._realgribkind
       zthe1=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(39))
       louvmi=.false.

       !     1.2.30 rotated and stretched lat/lon (get south pole of
       !            stretching in geographical lat/lon)
    elseif(kb2i(6)==30)then
       lolali=.true.
       loroti=.true.
       lostri=.true.
       zthei=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlami=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlovi=pb2i(kb2i(39))
       zlati(1)=-0.001_realgribkind*real(kb2i(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2i(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1, &
            zlami,zthei,zlovi,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe1=zlato(1)
       zlam1=zlono(1)
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(49))

       !     1.2.34 rotated and stretched gaussian lat/lon
    elseif(kb2i(6)==34)then
       lolali=.true.
       loroti=.true.
       lostri=.true.
       logaui=.true.
       zthei=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlami=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlovi=pb2i(kb2i(39))
       zlati(1)=-0.001_realgribkind*real(kb2i(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2i(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1, &
            zlami,zthei,zlovi,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe1=zlato(1)
       zlam1=zlono(1)
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(49))

       !     1.2.60 rotated lat/lon (hirlam definition)
    elseif(kb2i(6)==60)then
       lolali=.true.
       loroti=.true.
       zthe1=real(-kb2i(31),realgribkind)
       zlam1=real(kb2i(29),realgribkind)
       zlov1=0._realgribkind

       !     1.2.256 undefined
    else
       !        write(6,*)'undefined grid encountered, kb2i(6)=',kb2i(6)
       return
    endif

    !     1.3 idem for output grid, check if coinciding with input grid
    !
    lolalo=.false.
    logeoo=.false.
    logauo=.false.
    lopso=.false.
    loroto=.false.
    lostro=.false.
    locoin=.false.

    !     louvmo indicates whether to bother about u/v representation
    louvmo=.true.
    iuvtm=iuvtm+(1-mod(kb2o(17)/8,2))
    !
    !     1.3.0 geographical lat/lon
    if(kb2o(6)==0)then
       lolalo=.true.
       logeoo=.true.
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and. losame(zlam1,zlam2).and. &
            losame(zlov1,zlov2)
       louvmo=.false.

       !     1.3.2 polar-stereographic (defined by g.cats, knmi/dm-88-03)
    elseif(kb2o(6)==2)then
       lopso=.true.
       ztheo=real(kb2o(29),realgribkind)
       zlovo=0.01_realgribkind*real(kb2o(31),realgribkind)
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lopsi.and.  &
            losame(zlonoo,zlonoi).and. &
            losame(zlatoo,zlatoi).and. &
            losame(zthe1,zthe2).and. &
            losame(zlam1,zlam2).and. &
            losame(zlov1,zlov2).and. &
            losame(zthei,ztheo).and. &
            losame(zlovi,zlovo)

       !     1.3.4 gaussian geographical lat/lon
    elseif(kb2o(6)==4)then
       lolalo=.true.
       logeoo=.true.
       logauo=.true.
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lolali.and.logaui.and..not.lostri.and.&
            losame(zthe1,zthe2).and. &
            losame(zlam1,zlam2).and. &
            losame(zlov1,zlov2) 
       louvmo=.false.

       !     1.3.5 polar-stereographic
    elseif(kb2o(6)==5)then
       lopso=.true.
       ztheo=60._realgribkind
       if(kb2o(27)==1)ztheo=-ztheo
       zlovo=0.001_realgribkind*real(kb2o(18),realgribkind)
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       zlati(1)=0.001_realgribkind*real(kb2o(11),realgribkind)
       zloni(1)=0.001_realgribkind*real(kb2o(14),realgribkind)
       call gb2lpc(zloni,zlati,zlono,zlato,1,ztheo,zlovo)
       zlonoo=zlono(1)
       zlatoo=zlato(1)
       locoin=lopsi.and.          &
            losame(zlonoo,zlonoi).and.  &
            losame(zlatoo,zlatoi).and.  &
            losame(zthe1,zthe2).and.    &
            losame(zlam1,zlam2).and.    &
            losame(zlov1,zlov2).and.    &
            losame(zthei,ztheo).and.    &
            losame(zlovi,zlovo)

       !     1.3.10 rotated lat/lon
    elseif(kb2o(6)==10)then
       lolalo=.true.
       loroto=.true.
       zthe2=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlam2=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlov2=pb2o(kb2o(39))
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                            &
            losame(zlam1,zlam2).and.                            &
            losame(zlov1,zlov2)

       !     1.3.14 rotated gaussian lat/lon
    elseif(kb2o(6)==14)then
       lolalo=.true.
       loroto=.true.
       logauo=.true.
       zthe2=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlam2=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlov2=pb2o(kb2o(39))
       locoin=lolali.and.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2)

       !     1.3.20 stretched lat/lon
    elseif(kb2o(6)==20)then
       lolalo=.true.
       lostro=.true.
       ztheo=-90._realgribkind
       zlamo=0._realgribkind
       zlovo=0._realgribkind
       zthe2=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam2= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(39))
       locoin=lolali.and..not.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2).and.                       &
            losame(zthei,ztheo).and.                       &
            losame(zlami,zlamo).and.                       &
            losame(zlovi,zlovo).and.                       &
            losame(zstri,zstro)
       louvmo=.false.

       !     1.3.24 stretched gaussian lat/lon
    elseif(kb2o(6)==24)then
       lolalo=.true.
       lostro=.true.
       logauo=.true.
       ztheo=-90._realgribkind
       zlamo=0._realgribkind
       zlovo=0._realgribkind
       zthe2=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam2= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(39))
       locoin=lolali.and.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                  &
            losame(zlam1,zlam2).and.                  &
            losame(zlov1,zlov2).and.                  &
            losame(zthei,ztheo).and.                  &
            losame(zlami,zlamo).and.                  &
            losame(zlovi,zlovo).and.                  &
            losame(zstri,zstro)
       louvmo=.false.

       !     1.3.30 rotated and stretched lat/lon (get pole of
       !            stretching in geographical lat/lon)
    elseif(kb2o(6)==30)then
       lolalo=.true.
       loroto=.true.
       lostro=.true.
       ztheo=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlamo=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlovo=pb2o(kb2o(39))
       zlati(1)=-0.001_realgribkind*real(kb2o(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2o(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,zlamo,ztheo,zlovo,&
            0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe2=zlato(1)
       zlam2=zlono(1)
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(49))
       locoin=lolali.and..not.logaui.and.lostri.and.&
            losame(zthe1,zthe2).and.                      &
            losame(zlam1,zlam2).and.                      &
            losame(zlov1,zlov2).and.                      &
            losame(zthei,ztheo).and.                      &
            losame(zlami,zlamo).and.                      &
            losame(zlovi,zlovo).and.                      &
            losame(zstri,zstro)

       !     1.3.34 rotated and stretched gaussian lat/lon
    elseif(kb2o(6)==34)then
       lolalo=.true.
       loroto=.true.
       lostro=.true.
       logauo=.true.
       ztheo=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlamo=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlovo=pb2o(kb2o(39))
       zlati(1)=-0.001_realgribkind*real(kb2o(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2o(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,zlamo,ztheo,zlovo,&
            0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe2=zlato(1)
       zlam2=zlono(1)
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(49))
       locoin=lolali.and.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                  &
            losame(zlam1,zlam2).and.                  &
            losame(zlov1,zlov2).and.                  &
            losame(zthei,ztheo).and.                  &
            losame(zlami,zlamo).and.                  &
            losame(zlovi,zlovo).and.                  &
            losame(zstri,zstro)

       !     1.3.60 rotated lat/lon (hirlam definition)
    elseif(kb2o(6)==60)then
       lolalo=.true.
       loroto=.true.
       zthe2=real(-kb2o(31),realgribkind)
       zlam2=real(kb2o(29),realgribkind)
       zlov2=0._realgribkind
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                            &
            losame(zlam1,zlam2).and.                            &
            losame(zlov1,zlov2)

       !     1.3.256 undefined
    else
       !        write(6,*)'undefined grid encountered, kb2o(6)=',kb2o(6)
       return
    endif

    !     1.4 find cheapest combination of locoin, iuvtm, louvmi, louvmo
    !
    if(louvmi)then
       if(louvmo)then
          if(iuvtm==1.or.iuvtm==2)locoin=.false.
       else
          if(iuvtm==1)then
             iuvtm=0
          elseif(iuvtm==2)then
             iuvtm=3
          endif
       endif
    else
       if(louvmo)then
          if(iuvtm==1)then
             iuvtm=3
          elseif(iuvtm==2)then
             iuvtm=0
          endif
       else
          iuvtm=3
       endif
    endif

    !     1.5 treat coinciding grids
    !
    if(locoin)then
       do 1401 j=1,k
          pxnew(j)=pxold(j)
          pynew(j)=pyold(j)
1401   enddo
       return
    endif

    !     1.6 treat features not implemented
    !
    if(logaui.or.logauo)then
       !        write(6,*)'gaussian grid not implemented'
       return
    endif

    !     2.  input grid is polar-stereographic
    !
    if(lopsi)then
       if(kb2i(6)==2)then
          do 2002 j=1,k
             pxnew(j)=pxold(j)*1000._realgribkind
             pynew(j)=pyold(j)*1000._realgribkind
2002      enddo
       else
          do 2005 j=1,k
             pxnew(j)=pxold(j)+zlonoi
             pynew(j)=pyold(j)+zlatoi
2005      enddo
       endif
       lonew=.true.

       !     2.1 output grid is also polar-stereographic
       !
       if(lopso.and.(iuvtm==0.or.iuvtm==3))then
          if(iuvtm==0)then
             call gb2ppa(pxnew,pynew,pxnew,pynew,pdiri,k,&
                  zthei,zlovi,ztheo,zlovo)
          else
             call gb2ppc(pxnew,pynew,pxnew,pynew,k,    &
                  zthei,zlovi,ztheo,zlovo)
          endif
          if(kb2o(6)==2)then
             do 2112 j=1,k
                pxnew(j)=0.001_realgribkind*pxnew(j)
                pynew(j)=0.001_realgribkind*pynew(j)
2112         enddo
          elseif(kb2o(6)==5)then
             do 2115 j=1,k
                pxnew(j)=pxnew(j)-zlonoo
                pynew(j)=pynew(j)-zlatoo
2115         enddo
          endif
          return
       endif

       !     2.2 else convert to regular lat/lon
       !
       if(iuvtm==0.or.iuvtm==1)then
          call gb2pla(pxnew,pynew,pxnew,pynew,pdiri,k,zthei,zlovi)
          if(iuvtm==1)iuvtm=3
       else
          call gb2plc(pxnew,pynew,pxnew,pynew,k,zthei,zlovi)
          if(iuvtm==2)iuvtm=0
       endif

       !     2.5 if output grid is regular lat/lon return inmediately
       !
       locoin=lolalo.and..not.logauo.and..not.lostro.and. &
            zthe1==zthe2.and.zlam1==zlam2.and.zlov1==zlov2
       if(locoin)return

       !     2.9  end polar-sterographic input grid
       !
    endif
    !
    !     3.  treat stretched input grid
    !
    if(lostri)then
       !
       !     3.1 output grid is also stretched, with same pole-of-stretching
       !
       !
       !     3.2 convert from user pole to pole of stretching, no iuvtm change
       !
       if(lonew)then
          if(iuvtm==0.or.iuvtm==1)then
             call gb2lla(pxnew,pynew,pxnew,pynew,pdiri,k, &
                  zlami,zthei,zlovi,zlam1,zthe1,zlov1)
          else
             call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
                  zlami,zthei,zlovi,zlam1,zthe1,zlov1)
          endif
       else
          if(iuvtm==0.or.iuvtm==1)then
             call gb2lla(pxold,pyold,pxnew,pynew,pdiri,k,     & 
                  zlami,zthei,zlovi,zlam1,zthe1,zlov1)
          else
             call gb2llc(pxold,pyold,pxnew,pynew,k,     &
                  zlami,zthei,zlovi,zlam1,zthe1,zlov1)
          endif
          lonew=.true.
       endif

       !     3.3 de-stretch
       !
       call gb2tlc(pynew,pynew,k,1._realgribkind/zstri)
       !
       !     3.9 end stretched input grid
       !
    endif
    !
    !     5.  rotate regular lat/lon grid, but first adjust iuvtm
    !
    !     5.1 adjust iuvtm by rotating to geographical lat/lon
    !
    if(iuvtm==1)then
       if(lonew)then
          call gb2lla(pxnew,pynew,pxnew,pynew,pdiri,k,     &
               zlam1,zthe1,zlov1,0._realgribkind,-90._realgribkind,0._realgribkind)
       else
          call gb2lla(pxold,pyold,pxnew,pynew,pdiri,k,     &
               zlam1,zthe1,zlov1,0._realgribkind,-90._realgribkind,0._realgribkind)
          lonew=.true.
       endif
       lolali=.true.
       logeoi=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       iuvtm=3
    elseif(iuvtm==2)then
       if(lonew)then
          call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
               zlam1,zthe1,zlov1,0._realgribkind,-90._realgribkind,0._realgribkind)
       else
          call gb2llc(pxold,pyold,pxnew,pynew,k,     &
               zlam1,zthe1,zlov1,0._realgribkind,-90._realgribkind,0._realgribkind)
          lonew=.true.
       endif
       lolali=.true.
       logeoi=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       iuvtm=0
    endif

    !     5.2 rotation
    !
    if(lonew)then
       if(iuvtm==0)then
          call gb2lla(pxnew,pynew,pxnew,pynew,pdiri,k,     &
               zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
       else
          call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
               zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
       endif
    else
       if(iuvtm==0)then
          call gb2lla(pxold,pyold,pxnew,pynew,pdiri,k,     &
               zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
       else
          call gb2llc(pxold,pyold,pxnew,pynew,k,     & 
               zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
       endif
       lonew=.true.
    endif

    !     7.  treat stretched output grid
    !
    if(lostro)then
       !
       !     7.1 stretch
       !
       call gb2tlc(pynew,pynew,k,zstro)
       !
       !     7.2 convert from pole of stretching to user pole
       !
       if(iuvtm==0)then
          call gb2lla(pxnew,pynew,pxnew,pynew,pdiri,k,     &
               zlam2,zthe2,zlov2,zlamo,ztheo,zlovo)
       else
          call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
               zlam2,zthe2,zlov2,zlamo,ztheo,zlovo)
       endif

       !     7.9 end stretched output grid
       !
    endif
    !
    !     8.  convert to polar-stereoghraphic output grid
    !
    if(lopso)then
       if(iuvtm==0)then
          call gb2lpa(pxnew,pynew,pxnew,pynew,pdiri,k,     &
               ztheo,zlovo)
       else
          call gb2lpc(pxnew,pynew,pxnew,pynew,k,     &
               ztheo,zlovo)
       endif
       if(kb2o(6)==2)then
          do 8112 j=1,k
             pxnew(j)=0.001_realgribkind*pxnew(j)
             pynew(j)=0.001_realgribkind*pynew(j)
8112      enddo
       elseif(kb2o(6)==5)then
          do 8115 j=1,k
             pxnew(j)=pxnew(j)-zlonoo
             pynew(j)=pynew(j)-zlatoo
8115      enddo
       endif
    endif

    return
  end subroutine gb2gca

  subroutine gb2gcc(pxold,pyold,pxnew,pynew,k,kb2i,pb2i,kb2o,pb2o)
    !
    !     subroutine gb2gcc: transform coordinates as per gribblocks 2
    !
    ! input:
    !     pxold,pyold: input coordinates (metres, degrees)
    !     k:           number of coordinate pairs to be transformed
    !     kb2i,pb2i:   integer resp. real values in block 2
    !                  describing input grid (cf grib code fm 92)
    !     kb2o,pb2o:   idem, describing output grid
    !
    ! output:
    !     pxnew,pynew: transformde coordinates
    !
    ! notes:

    !     1.           pxold may be allocated the same memory words
    !                  as pxnew; this routine then overwrites pxold
    !                  similarly for pyold and pynew.
    !     3.           the real parameters in the fm-92 grib block 2
    !                  definition are expected in the arrays pb2i/o
    !                  the corresponding integer elements of kb2i/o
    !                  contain pointers to the relevant element of
    !                  pb2i/o.  e.g., for rotated  and stretched
    !                  regular lat/lon grids, we could have:
    !                  kb2(39)=1, to indicate that
    !                  pb2(1) is the angle of rotation, and
    !                  kb2(49)=2, i.e. pb2(kb2(49)) is the
    !                  stretching factor.
    !     4.           implemented are:
    !                  a)regular lat/lon projections,
    !                  perhaps stretched, perhaps rotated
    !                  (coordinates in degrees wrt 0 degrees lat,
    !                  0 degrees lon)
    !                  b)polar stereographic projections
    !                  (coordinates are in metres wrt the grid
    !                  lower left corner, which itself is
    !                  given in 0.001 degrees in the geographical
    !                  lat/lon grid, as words 11 and 14 of grib
    !                  block 2)
    !
    ! method:
    !     1. set a number of parameters to define in- and output grids
    !     2. convert polar-stereographic input to either polar-stereo
    !        graphic output or geographical lat/lon (i.e. pole at 90)
    !     3. destretch stretched input to regular lat/lon with pole
    !        at the pole of stretching
    !     4. de-gauss gaussian input (not implemented)
    !        or de-spectralise spectral grid (not implemented)
    !     5. convert regular lat/lon to regular lat/lon
    !        (for stretched output: pole at pole of stretching)
    !     6. interpolate to gaussian output grid (not implemented)
    !        or spectralise to spectral grid (not implemented)
    !     7. stretch stretched output grid
    !     8. convert regular lat/lon to polar-stereographic output
    !
    !     the meaning of the local variables is as follows:
    !     names ending in 1 or i refer to the input grid,
    !     those ending in 2 or o to the output grid.
    !     for lat/lon grids,
    !     (zlam,zthe) give longitude and latitude of the pole,
    !     zlov the angle of rotation, for polar-stereograpihic grids
    !     zthe is the latitude of the projection plane, and zlov the
    !     grid orientation.
    !     names ending in 1 and 2 define the grids involved in the
    !     central transformation (lat/lon to lat/lon); for lat/lon
    !     grids these are the in- or output grids, and for the
    !     following cases these grids are the geographical
    !     lat/lon projections:
    !     a. polar-stereographic projections,
    !     b. gaussian grids,
    !     c. spectral grids,
    !     while for (d.) stretched grids these grids are the
    !     lat/lon projections, with the northpole rotated to the
    !     pole of stretching.
    !     names ending in i or o define the grids from resp. to which
    !     the grids defined by names ending
    !     in 1 or 2 have to be transformed to get the
    !     in- resp. output grids in resp. from lat/lon.
    !     schematically: i to 1 (1 is regular lat/lon)
    !                    1 to 2 (2 is regular lat/lon)
    !                    2 to o
    !     short-cuts are taken for specific combinations of
    !     in- and output grids
    !
    !                                 g.j.cats, 05 jan 89.
    !
!!!      implicit logical(l)
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*)
    integer kb2i(*),kb2o(*)
    real(kind=realgribkind) pb2i(*),pb2o(*)
    real(kind=realgribkind) zlono(1),zlato(1),zlati(1),zloni(1)
    integer j,k
    logical losame
    real(kind=realgribkind) p1,p2,zsame,zstro

    losame(p1,p2)=abs(p1-p2)<zsame
    data zsame/0.01_realgribkind/
    logical lonew,lopsi,lolali,logeoi,logaui,loroti,lostri,lolalo
    logical logeoo,logauo,lopso,loroto,lostro,locoin
    real(kind=realgribkind) zthe1,zlam1,zlov1,zthei,zlovi,zlonoi,zlatoi,zlami
    real(kind=realgribkind) zstri,zthe2,zlam2,zlov2,ztheo,zlovo,zlonoo,zlatoo,zlamo


    lonew=.false.

    lopsi=.false.
    lolali=.false.
    logeoi=.false.
    logaui=.false.
    loroti=.false.
    lostri=.false.

    !     1.2.0 geographical lat/lon
    if(kb2i(6)==0)then
       lolali=.true.
       logeoi=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind

       !     1.2.2 polar-stereographic (defined by g.cats, knmi/dm-88-03)
    elseif(kb2i(6)==2)then
       lopsi=.true.
       zthei=real(kb2i(29),realgribkind)
       zlovi=0.01_realgribkind*real(kb2i(31),realgribkind)
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       zlonoi=0._realgribkind
       zlatoi=0._realgribkind

       !     1.2.4 gaussian geographical lat/lon
    elseif(kb2i(6)==4)then
       lolali=.true.
       logeoi=.true.
       logaui=.true.
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind

       !     1.2.5 polar-stereographic
    elseif(kb2i(6)==5)then
       lopsi=.true.
       zthei=60._realgribkind
       if(kb2i(27)==1)zthei=-zthei
       zlovi=0.001_realgribkind*real(kb2i(18),realgribkind)
       zthe1=-90._realgribkind
       zlam1=0._realgribkind
       zlov1=0._realgribkind
       zlati(1)=0.001_realgribkind*real(kb2i(11),realgribkind)
       zloni(1)=0.001_realgribkind*real(kb2i(14),realgribkind)
       call gb2lpc(zloni,zlati,zlono,zlato,1,zthei,zlovi)
       zlonoi=zlono(1)
       zlatoi=zlato(1)

       !     1.2.10 rotated lat/lon
    elseif(kb2i(6)==10)then
       lolali=.true.
       loroti=.true.
       zthe1=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlov1=pb2i(kb2i(39))

       !     1.2.14 rotated gaussian lat/lon
    elseif(kb2i(6)==14)then
       lolali=.true.
       loroti=.true.
       logaui=.true.
       zthe1=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlov1=pb2i(kb2i(39))

       !     1.2.20 stretched lat/lon
    elseif(kb2i(6)==20)then
       lolali=.true.
       lostri=.true.
       zthei=-90._realgribkind
       zlami=0._realgribkind
       zlovi=0._realgribkind
       zthe1=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(39))

       !     1.2.24 stretched gaussian lat/lon
    elseif(kb2i(6)==24)then
       lolali=.true.
       lostri=.true.
       logaui=.true.
       zthei=-90._realgribkind
       zlami=0._realgribkind
       zlovi=0._realgribkind
       zthe1=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam1= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(39))

       !     1.2.30 rotated and stretched lat/lon (get south pole of
       !            stretching in geographical lat/lon)
    elseif(kb2i(6)==30)then
       lolali=.true.
       loroti=.true.
       lostri=.true.
       zthei=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlami=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlovi=pb2i(kb2i(39))
       zlati(1)=-0.001_realgribkind*real(kb2i(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2i(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,     &
            zlami,zthei,zlovi,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe1=zlato(1)
       zlam1=zlono(1)
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(49))

       !     1.2.34 rotated and stretched gaussian lat/lon
    elseif(kb2i(6)==34)then
       lolali=.true.
       loroti=.true.
       lostri=.true.
       logaui=.true.
       zthei=0.001_realgribkind*real(kb2i(33),realgribkind)
       zlami=0.001_realgribkind*real(kb2i(36),realgribkind)
       zlovi=pb2i(kb2i(39))
       zlati(1)=-0.001_realgribkind*real(kb2i(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2i(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,     &
            zlami,zthei,zlovi,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe1=zlato(1)
       zlam1=zlono(1)
       zlov1=0._realgribkind
       zstri=pb2i(kb2i(49))

       !     1.2.60 rotated lat/lon (hirlam definition)
    elseif(kb2i(6)==60)then
       lolali=.true.
       loroti=.true.
       zthe1=real(-kb2i(31),realgribkind)
       zlam1=real(kb2i(29),realgribkind)
       zlov1=0._realgribkind

       !     1.2.256 undefined
    else
       !        write(6,*)'undefined grid encountered, kb2i(6)=',kb2i(6)
       return
    endif

    !     1.3 idem for output grid, check if coinciding with input grid
    !
    lolalo=.false.
    logeoo=.false.
    logauo=.false.
    lopso=.false.
    loroto=.false.
    lostro=.false.
    locoin=.false.

    !     1.3.0 geographical lat/lon
    if(kb2o(6)==0)then
       lolalo=.true.
       logeoo=.true.
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                            &
            losame(zlam1,zlam2).and.                            &
            losame(zlov1,zlov2)

       !     1.3.2 polar-stereographic (defined by g.cats, knmi/dm-88-03)
    elseif(kb2o(6)==2)then
       lopso=.true.
       ztheo=real(kb2o(29),realgribkind)
       zlovo=0.01_realgribkind*real(kb2o(31),realgribkind)
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lopsi.and.         &
            losame(zlonoo,zlonoi).and. &
            losame(zlatoo,zlatoi).and. &
            losame(zthe1,zthe2).and.   &
            losame(zlam1,zlam2).and.   &
            losame(zlov1,zlov2).and.   &
            losame(zthei,ztheo).and.   &
            losame(zlovi,zlovo)

       !     1.3.4 gaussian geographical lat/lon
    elseif(kb2o(6)==4)then
       lolalo=.true.
       logeoo=.true.
       logauo=.true.
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       locoin=lolali.and.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2)

       !     1.3.5 polar-stereographic
    elseif(kb2o(6)==5)then
       lopso=.true.
       ztheo=60._realgribkind
       if(kb2o(27)==1)ztheo=-ztheo
       zlovo=0.001_realgribkind*real(kb2o(18),realgribkind)
       zthe2=-90._realgribkind
       zlam2=0._realgribkind
       zlov2=0._realgribkind
       zlati(1)=0.001_realgribkind*real(kb2o(11),realgribkind)
       zloni(1)=0.001_realgribkind*real(kb2o(14),realgribkind)
       call gb2lpc(zloni,zlati,zlono,zlato,1,ztheo,zlovo)
       zlonoo=zlono(1)
       zlatoo=zlato(1)
       locoin=lopsi.and.         &
            losame(zlonoo,zlonoi).and. &
            losame(zlatoo,zlatoi).and. &
            losame(zthe1,zthe2).and.   &
            losame(zlam1,zlam2).and.   &
            losame(zlov1,zlov2).and.   &
            losame(zthei,ztheo).and.   &
            losame(zlovi,zlovo)

       !     1.3.10 rotated lat/lon
    elseif(kb2o(6)==10)then
       lolalo=.true.
       loroto=.true.
       zthe2=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlam2=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlov2=pb2o(kb2o(39))
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                            &
            losame(zlam1,zlam2).and.                            &
            losame(zlov1,zlov2)

       !     1.3.14 rotated gaussian lat/lon
    elseif(kb2o(6)==14)then
       lolalo=.true.
       loroto=.true.
       logauo=.true.
       zthe2=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlam2=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlov2=pb2o(kb2o(39))
       locoin=lolali.and.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2)

       !     1.3.20 stretched lat/lon
    elseif(kb2o(6)==20)then
       lolalo=.true.
       lostro=.true.
       ztheo=-90._realgribkind
       zlamo=0._realgribkind
       zlovo=0._realgribkind
       zthe2=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam2= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(39))
       locoin=lolali.and..not.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2).and.                       &
            losame(zthei,ztheo).and.                       &
            losame(zlami,zlamo).and.                       &
            losame(zlovi,zlovo).and.                       &
            losame(zstri,zstro)

       !     1.3.24 stretched gaussian lat/lon
    elseif(kb2o(6)==24)then
       lolalo=.true.
       lostro=.true.
       logauo=.true.
       ztheo=-90._realgribkind
       zlamo=0._realgribkind
       zlovo=0._realgribkind
       zthe2=-0.001_realgribkind*real(kb2i(33),realgribkind)
       zlam2= 0.001_realgribkind*real(kb2i(36),realgribkind)+180._realgribkind
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(39))
       locoin=lolali.and.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                  &
            losame(zlam1,zlam2).and.                  &
            losame(zlov1,zlov2).and.                  &
            losame(zthei,ztheo).and.                  &
            losame(zlami,zlamo).and.                  &
            losame(zlovi,zlovo).and.                  &
            losame(zstri,zstro)

       !     1.3.30 rotated and stretched lat/lon (get pole of
       !            stretching in geographical lat/lon)
    elseif(kb2o(6)==30)then
       lolalo=.true.
       loroto=.true.
       lostro=.true.
       ztheo=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlamo=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlovo=pb2o(kb2o(39))
       zlati(1)=-0.001_realgribkind*real(kb2o(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2o(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,    &
            zlamo,ztheo,zlovo,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe2=zlato(1)
       zlam2=zlono(1)
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(49))
       locoin=lolali.and..not.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                       &
            losame(zlam1,zlam2).and.                       &
            losame(zlov1,zlov2).and.                       &
            losame(zthei,ztheo).and.                       &
            losame(zlami,zlamo).and.                       &
            losame(zlovi,zlovo).and.                       &
            losame(zstri,zstro)

       !     1.3.34 rotated and stretched gaussian lat/lon
    elseif(kb2o(6)==34)then
       lolalo=.true.
       loroto=.true.
       lostro=.true.
       logauo=.true.
       ztheo=0.001_realgribkind*real(kb2o(33),realgribkind)
       zlamo=0.001_realgribkind*real(kb2o(36),realgribkind)
       zlovo=pb2o(kb2o(39))
       zlati(1)=-0.001_realgribkind*real(kb2o(43),realgribkind)
       zloni(1)= 0.001_realgribkind*real(kb2o(46),realgribkind)+180._realgribkind
       call gb2llc(zloni,zlati,zlono,zlato,1,     &  
            zlamo,ztheo,zlovo,0._realgribkind,-90._realgribkind,0._realgribkind)
       zthe2=zlato(1)
       zlam2=zlono(1)
       zlov2=0._realgribkind
       zstro=pb2o(kb2o(49))
       locoin=lolali.and.logaui.and.lostri.and. &
            losame(zthe1,zthe2).and.                  &
            losame(zlam1,zlam2).and.                  &
            losame(zlov1,zlov2).and.                  &
            losame(zthei,ztheo).and.                  &
            losame(zlami,zlamo).and.                  &
            losame(zlovi,zlovo).and.                  &
            losame(zstri,zstro)

       !     1.3.60 rotated lat/lon (hirlam definition)
    elseif(kb2o(6)==60)then
       lolalo=.true.
       loroto=.true.
       zthe2=real(-kb2o(31),realgribkind)
       zlam2=real(kb2o(29),realgribkind)
       zlov2=0._realgribkind
       locoin=lolali.and..not.logaui.and..not.lostri.and. &
            losame(zthe1,zthe2).and.                            &
            losame(zlam1,zlam2).and.                            &
            losame(zlov1,zlov2)

       !     1.3.256 undefined
    else
       !        write(6,*)'undefined grid encountered, kb2o(6)=',kb2o(6)
       return
    endif

    !     1.4 treat coinciding grids
    !
    if(locoin)then
       do 1401 j=1,k
          pxnew(j)=pxold(j)
          pynew(j)=pyold(j)
1401   enddo
       return
    endif

    !     1.5 treat features not implemented
    !
    if(logaui.or.logauo)then
       !        write(6,*)'gaussian grid not implemented'
       return
    endif

    !     2.  input grid is polar-stereographic
    !
    if(lopsi)then
       if(kb2i(6)==2)then
          do 2002 j=1,k
             pxnew(j)=pxold(j)*1000._realgribkind
             pynew(j)=pyold(j)*1000._realgribkind
2002      enddo
       else
          do 2005 j=1,k
             pxnew(j)=pxold(j)+zlonoi
             pynew(j)=pyold(j)+zlatoi
2005      enddo
       endif
       lonew=.true.

       !     2.1 output grid is also polar-stereographic
       !
       if(lopso)then
          call gb2ppc(pxnew,pynew,pxnew,pynew,k,     &
               zthei,zlovi,ztheo,zlovo)
          if(kb2o(6)==2)then
             do 2112 j=1,k
                pxnew(j)=0.001_realgribkind*pxnew(j)
                pynew(j)=0.001_realgribkind*pynew(j)
2112         enddo
          elseif(kb2o(6)==5)then
             do 2115 j=1,k
                pxnew(j)=pxnew(j)-zlonoo
                pynew(j)=pynew(j)-zlatoo
2115         enddo
          endif
          return
       endif

       !     2.2 else convert to regular lat/lon
       !
       call gb2plc(pxnew,pynew,pxnew,pynew,k,zthei,zlovi)
       !
       !     2.5 if output grid is regular lat/lon return inmediately
       !
       locoin=lolalo.and..not.logauo.and..not.lostro.and.     & 
            zthe1==zthe2.and.zlam1==zlam2.and.zlov1==zlov2
       if(locoin)return

       !     2.9  end polar-sterographic input grid
       !
    endif
    !
    !     3.  treat stretched input grid
    !
    if(lostri)then
       !
       !     3.1 output grid is also stretched, with same pole-of-stretching
       !
       !     3.2 convert from user pole to pole of stretching
       !
       if(lonew)then
          call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
               zlami,zthei,zlovi,zlam1,zthe1,zlov1)
       else
          call gb2llc(pxold,pyold,pxnew,pynew,k,     &
               zlami,zthei,zlovi,zlam1,zthe1,zlov1)
          lonew=.true.
       endif

       !     3.3 de-stretch
       !
       call gb2tlc(pynew,pynew,k,1._realgribkind/zstri)
       !
       !     3.9 end stretched input grid
       !
    endif
    !
    !     5.  rotate regular lat/lon grid
    !
    if(lonew)then
       call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
            zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
    else
       call gb2llc(pxold,pyold,pxnew,pynew,k,     &
            zlam1,zthe1,zlov1,zlam2,zthe2,zlov2)
       lonew=.true.
    endif

    !     7.  treat stretched output grid
    !
    if(lostro)then
       !
       !     7.1 stretch
       !
       call gb2tlc(pynew,pynew,k,zstro)
       !
       !     7.2 convert from pole of stretching to user pole
       !
       call gb2llc(pxnew,pynew,pxnew,pynew,k,     &
            zlam2,zthe2,zlov2,zlamo,ztheo,zlovo)

       !     7.9 end stretched output grid
       !
    endif
    !
    !     8.  convert to polar-stereoghraphic output grid
    !
    if(lopso)then
       call gb2lpc(pxnew,pynew,pxnew,pynew,k,     & 
            ztheo,zlovo)
       if(kb2o(6)==2)then
          do 8112 j=1,k
             pxnew(j)=0.001_realgribkind*pxnew(j)
             pynew(j)=0.001_realgribkind*pynew(j)
8112      enddo
       elseif(kb2o(6)==5)then
          do 8115 j=1,k
             pxnew(j)=pxnew(j)-zlonoo
             pynew(j)=pynew(j)-zlatoo
8115      enddo
       endif
    endif
    return
  end subroutine gb2gcc
  subroutine gb2lla(pxold,pyold,pxnew,pynew,pdiri,k,     &
       plam1,pthe1,phi1,plam2,pthe2,phi2)
    !
    !     subroutine gb2lla: as gb2llc, but also calculate angle of rotation
    !
    !     input:
    !     pxold,pyold:     arrays of longitudes,latitudes to be transformed
    !     pdiri:           array of initial orientation angles
    !     k:               number of gridpoints to be transformed
    !     plam1,pthe1,phi1:parameters defining input grid
    !     plam2,pthe2,phi2:parameters defining output grid
    !                      plam: longitude of south pole
    !                      pthe: latitude of south pole
    !                      phi:  angle of rotation as defined by grib
    !                      note: coordinates of south poles of in- and
    !                            output grids must be given in the same
    !                            coordinate system, e.g. geographical
    !
    !     output:
    !     pxnew,pynew: arrays of transformed longitudes,latitudes
    !     pdiri:           array of incremented orientation angles
    !
    !     all parameters to be in degrees,
    !     output is between -180(excl) and 180. (incl)
    !
    !     in- and output arrays of coordinates may coincide
    !
    !     all comments concerning the parameter pdiri in subroutine
    !     gb2sla are applicable here - in fact, gb2lla is simply a
    !     copy of gb2llc with gb2slc replaced by gb2sla
    !
    !     it is anticipated that this routine will mainly be used
    !     to transform from or to the geographical system, therefore
    !     the cases with pthe1=-90 or pthe2=-90 are treated separately
    !
    !                                                 g.j.cats  16 dec 88.
    !
    !      implicit logical(l)
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    real(kind=realgribkind) zlasp(1),zlosp(1)
    real(kind=realgribkind) zsame,pthe1,zlam1,plam1,pthe2,zlam2,plam2,phi2,phi1,zrot
    integer k,j
    data zsame/0.01_realgribkind/

    !     1.  preliminaries
    !
    !     1.1 eliminate ambiguities near coordinate singularities
    !
    if(pthe1<-89.9_realgribkind)then
       zlam1=0._realgribkind
    else
       zlam1=plam1
    endif
    if(pthe2<-89.9_realgribkind)then
       zlam2=0._realgribkind
    else
       zlam2=plam2
    endif

    !     2.  special cases:
    !
    !     2.1 south poles coincide
    !
    if(abs(pthe1-pthe2)<zsame.and.     &
         (abs(zlam1-zlam2)<zsame.or.   & 
         abs(zlam1-zlam2)>360._realgribkind-zsame))then
       call gb2rlc(pxold,pxnew,k,phi2-phi1)
       do 2101 j=1,k
          pynew(j)=pyold(j)
2101   enddo

       !     2.2 input grid is geographical
       !
    elseif(pthe1<-89.9_realgribkind)then
       call gb2rlc(pxold,pxnew,k,zlam2-phi1)
       call gb2sla(pxnew,pyold,pxnew,pynew,pdiri,k,pthe2)
       call gb2rlc(pxnew,pxnew,k,phi2)

       !     2.3 output grid is geographical
       !
    elseif(pthe2<-89.9_realgribkind)then
       call gb2rlc(pxold,pxnew,k,180._realgribkind-phi1)
       call gb2sla(pxnew,pyold,pxnew,pynew,pdiri,k,pthe1)
       zrot=phi2+180._realgribkind-zlam1
       if(zrot> 180._realgribkind)zrot=zrot-360._realgribkind
       call gb2rlc(pxnew,pxnew,k,zrot)

       !     3.  general case
       !
    else
       !
       !     3.1 coordinates of output grid south pole in input grid
       !
       zlasp(1)=pthe2
       zlosp(1)=zlam2-zlam1
       call gb2slc(zlosp,zlasp,zlosp,zlasp,1,pthe1)

       !     3.2 coordinate transformation from one pole to the other
       !
       call gb2rlc(pxold,pxnew,k,zlosp(1)-phi1)
       call gb2sla(pxnew,pyold,pxnew,pynew,pdiri,k,zlasp(1))
       !
       !     3.3 pole1 is on 180w longitude of pole2 grid; rotate back
       !
       zlosp(1)=zlam1-zlam2
       zlasp(1)=pthe1
       call gb2slc(zlosp,zlasp,zlosp,zlasp,1,pthe2)
       zrot=phi2-zlosp(1)+180._realgribkind
       if(zrot> 180._realgribkind)zrot=zrot-360._realgribkind
       call gb2rlc(pxnew,pxnew,k,zrot)

    endif
    return
  end subroutine gb2lla
  subroutine gb2llc(pxold,pyold,pxnew,pynew,k,     &
       plam1,pthe1,phi1,plam2,pthe2,phi2)
    !
    !     subroutine gb2llc: transformation from one latlon grid to another
    !
    !     input:
    !     pxold,pyold:     arrays of longitudes,latitudes to be transformed
    !     k:               number of gridpoints to be transformed
    !     plam1,pthe1,phi1:parameters defining input grid
    !     plam2,pthe2,phi2:parameters defining output grid
    !                      plam: longitude of south pole
    !                      pthe: latitude of south pole
    !                      phi:  angle of rotation as defined by grib
    !                      note: coordinates of south poles of in- and
    !                            output grids must be given in the same
    !                            coordinate system, e.g. geographical
    !
    !     output:
    !     pxnew,pynew: arrays of transformed longitudes,latitudes
    !
    !     all parameters to be in degrees,
    !     output is between -180(excl) and 180. (incl)
    !
    !     in- and output arrays of coordinates may coincide
    !
    !     it is anticipated that this routine will mainly be used
    !     to transform from or to the geographical system, therefore
    !     the cases with pthe1=-90 or pthe2=-90 are treated separately
    !
    !                                                 g.j.cats  09 nov 88.
    !
    !      implicit logical(l)
    implicit none
    real(kind=realgribkind)  pxold(*),pyold(*),pxnew(*),pynew(*)
    real(kind=realgribkind) zlasp(1),zlosp(1)
    real(kind=realgribkind) pthe1,zlam1,plam1,pthe2,zlam2,plam2,phi2,phi1,zrot
    integer k,j
    !     1.1 eliminate ambiguities near coordinate singularities
    !
    if(pthe1<-89.9_realgribkind)then
       zlam1=0._realgribkind
    else
       zlam1=plam1
    endif
    if(pthe2<-89.9_realgribkind)then
       zlam2=0._realgribkind
    else
       zlam2=plam2
    endif

    !     2.  special cases:
    !
    !     2.1 south poles coincide
    !
    if(abs(pthe1-pthe2)<1.e-6_realgribkind.and.     &
         (abs(zlam1-zlam2)<1.e-6_realgribkind.or.     & 
         abs(zlam1-zlam2)>360._realgribkind-1.e-4_realgribkind))then
       call gb2rlc(pxold,pxnew,k,phi2-phi1)
       do 2101 j=1,k
          pynew(j)=pyold(j)
2101   enddo

       !     2.2 input grid is geographical
       !
    elseif(pthe1<-89.9_realgribkind)then
       call gb2rlc(pxold,pxnew,k,zlam2-phi1)
       call gb2slc(pxnew,pyold,pxnew,pynew,k,pthe2)
       call gb2rlc(pxnew,pxnew,k,phi2)

       !     2.3 output grid is geographical
       !
    elseif(pthe2<-89.9_realgribkind)then
       call gb2rlc(pxold,pxnew,k,180._realgribkind-phi1)
       call gb2slc(pxnew,pyold,pxnew,pynew,k,pthe1)
       zrot=phi2+180._realgribkind-zlam1
       if(zrot> 180._realgribkind)zrot=zrot-360._realgribkind
       call gb2rlc(pxnew,pxnew,k,zrot)

       !     3.  general case
       !
    else
       !
       !     3.1 coordinates of output grid south pole in input grid
       !
       zlasp(1)=pthe2
       zlosp(1)=zlam2-zlam1
       call gb2slc(zlosp,zlasp,zlosp,zlasp,1,pthe1)

       !     3.2 coordinate transformation from one pole to the other
       !
       call gb2rlc(pxold,pxnew,k,zlosp(1)-phi1)
       call gb2slc(pxnew,pyold,pxnew,pynew,k,zlasp(1))
       !
       !     3.3 pole1 is on 180w longitude of pole2 grid; rotate back
       !
       zlosp(1)=zlam1-zlam2
       zlasp(1)=pthe1
       call gb2slc(zlosp,zlasp,zlosp,zlasp,1,pthe2)
       zrot=phi2-zlosp(1)+180._realgribkind
       if(zrot> 180._realgribkind)zrot=zrot-360._realgribkind
       call gb2rlc(pxnew,pxnew,k,zrot)

    endif
    return
  end subroutine gb2llc
  subroutine gb2lpa(pxold,pyold,pxnew,pynew,pdiri,k,pthe,plov)
    !     
    !     subroutine gb2lpa: transform from geographical grid to polar-sterg
    !     
    !     input:
    !     pxold,pyold: arrays of longitudes,latitudes to be transformed
    !     pdiri:       array of initial orientation angles
    !     k:           number of gridpoints to be transformed
    !     pthe:        latitude of projection plane
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe)
    !     plov:        longitude of positive y-axis (grid orientation)
    !     
    !     output:
    !     pxnew,pynew: arrays of x and y coordinates
    !     pdiri:       array of incremented orientation angles
    !     
    !     all input parameters to be in degrees,
    !     output coordinates are in metres, pdiri is in degrees
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     refer to gb2sla for comments concerning the use of pdiri
    !     and the choice between gb2lpa and gb2lpc
    !     
    !     g.j.cats  02 jan 89.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    real(kind=realgribkind) zrho,zlac,pthe,zpol,plov,zlat,zr,zlam
    integer j,k


    data zrho/ 0.0174532925_realgribkind/
    zlac=6371000._realgribkind*(1._realgribkind+abs(sin(pthe*zrho)))
    zpol=sign(1._realgribkind,pthe)
    do 100 j=1,k
       pdiri(j)=pdiri(j)-(pxold(j)-plov)*zpol
       zlat=pyold(j)*zrho
       zr=zlac/(1._realgribkind+zpol*sin(zlat))*cos(zlat)
       zlam=(pxold(j)-90._realgribkind+zpol*sin(zlat))*cos(zlat)
       zlam=(pxold(j)-90._realgribkind-plov)*zrho
       pxnew(j)=zr*cos(zlam)
       pynew(j)=zr*sin(zlam)*zpol
100 enddo
    return
  end subroutine gb2lpa
  subroutine gb2lpc(pxold,pyold,pxnew,pynew,k,pthe,plov)
    !     
    !     subroutine gb2lpc: transform from geographical grid to polar-sterg
    !     
    !     input:
    !     pxold,pyold: arrays of longitudes,latitudes to be transformed
    !     k:           number of gridpoints to be transformed
    !     pthe:        latitude of projection plane
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe)
    !     plov:        longitude of positive y-axis (grid orientation)
    !     
    !     output:
    !     pxnew,pynew: arrays of x and y coordinates
    !     
    !     all input parameters to be in degrees,
    !     output is in metres
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     g.j.cats  02 jan 89.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*)
    real(kind=realgribkind) zrho,zlac,zpol,pthe,zlat,zr,zlam,plov
    integer j,k

    data zrho/ 0.0174532925_realgribkind/
    zlac=6371000._realgribkind*(1._realgribkind+abs(sin(pthe*zrho)))
    zpol=sign(1._realgribkind,pthe)
    do 100 j=1,k
       zlat=pyold(j)*zrho
       zr=zlac/(1._realgribkind+zpol*sin(zlat))*cos(zlat)
       zlam=(pxold(j)-90._realgribkind-plov)*zrho
       pxnew(j)=zr*cos(zlam)
       pynew(j)=zr*sin(zlam)*zpol
100 enddo
    return
  end subroutine gb2lpc
  subroutine gb2pla(pxold,pyold,pxnew,pynew,pdiri,k,pthe,plov)
    !     
    !     subroutine gb2pla: transform from polar-stereographic to lat/lon
    !     
    !     input:
    !     pxold,pyold: arrays of x and y coordinates to be transformed
    !     k:           number of gridpoints to be transformed
    !     pdiri:       array of initial orientation angles
    !     pthe:        latitude of projection plane
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe)
    !     plov:        longitude of positive y-axis (grid orientation)
    !     
    !     output:
    !     pxnew,pynew: arrays of latitudes and longitudes
    !     pdiri:       array of incremented orientation angles
    !     
    !     pxold and pyold in metres, pthe in degrees
    !     output is in degrees
    !     
    !     refer to gb2sla for comments concerning the use of pdiri
    !     and the choice between gb2plc and gb2pla
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     g.j.cats  03 jan 89.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    real(kind=realgribkind) zrad,zrho,zlac,zpol,zshf,pthe,plov,zr
    integer j,k

    data zrad/57.2957795_realgribkind   /
    data zrho/ 0.0174532925_realgribkind/
    zlac=6371000._realgribkind*(1._realgribkind+abs(sin(pthe*zrho)))
    zpol=sign(1._realgribkind,pthe)
    zshf=plov+90._realgribkind
    do 100 j=1,k
       zr=sqrt(pxold(j)**2._realgribkind+pyold(j)**2._realgribkind)
       if(abs(pxold(j))<1.e-14_realgribkind.and.abs(pyold(j))<1.e-14_realgribkind)then
          pxnew(j)=zshf
       else
          pxnew(j)=atan2(zpol*pyold(j),pxold(j))*zrad+zshf
       endif
       if(pxnew(j)> 180._realgribkind)pxnew(j)=pxnew(j)-360._realgribkind
       if(pxnew(j)<=-180._realgribkind)pxnew(j)=pxnew(j)+360._realgribkind
       pynew(j)=zpol*(90._realgribkind-2._realgribkind*atan(zr/zlac)*zrad)
       pdiri(j)=pdiri(j)+(pxnew(j)-plov)*zpol
100 enddo
    return
  end subroutine gb2pla
  subroutine gb2plc(pxold,pyold,pxnew,pynew,k,pthe,plov)
    !
    !     subroutine gb2plc: transform from polar-stereographic to lat/lon
    !
    !     input:
    !     pxold,pyold: arrays of x and y coordinates to be transformed
    !     k:           number of gridpoints to be transformed
    !     pthe:        latitude of projection plane
    !                  (it is assumed that the projection centre is the pole
    !                   of which the latitude has opposite sign to pthe)
    !     plov:        longitude of positive y-axis (grid orientation)
    !
    !     output:
    !     pxnew,pynew: arrays of latitudes and longitudes
    !
    !     pxold and pyold in metres, pthe in degrees
    !     output is in degrees
    !
    !     in- and output arrays of coordinates may coincide
    !
    !                                                 g.j.cats  03 jan 89.
    !
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*)
    real(kind=realgribkind) zrad,zrho,zlac,zpol,pthe,zshf,plov,zr
    integer j,k

    data zrad/57.2957795_realgribkind   /
    data zrho/ 0.0174532925_realgribkind/
    zlac=6371000._realgribkind*(1._realgribkind+abs(sin(pthe*zrho)))
    zpol=sign(1._realgribkind,pthe)
    zshf=plov+90._realgribkind
    do 100 j=1,k
       zr=sqrt(pxold(j)**2._realgribkind+pyold(j)**2._realgribkind)
       if(abs(pxold(j))<1.e-14_realgribkind.and.abs(pyold(j))<1.e-14_realgribkind)then
          pxnew(j)=zshf
       else
          pxnew(j)=atan2(zpol*pyold(j),pxold(j))*zrad+zshf
       endif
       if(pxnew(j)> 180._realgribkind)pxnew(j)=pxnew(j)-360._realgribkind
       if(pxnew(j)<=-180._realgribkind)pxnew(j)=pxnew(j)+360._realgribkind
       pynew(j)=zpol*(90._realgribkind-2._realgribkind*atan(zr/zlac)*zrad)
100 enddo
    return
  end subroutine gb2plc
  subroutine gb2ppa(pxold,pyold,pxnew,pynew,pdiri,k,pthe1,plov1,     &
       pthe2,plov2)
    !     
    !     subroutine gb2ppa: transform from polar-stereographic to polar-s
    !     
    !     input:
    !     pxold,pyold: arrays of x and y coordinates to be transformed
    !     k:           number of gridpoints to be transformed
    !     pdiri:       array of initial orientation angles
    !     pthe1:       latitude of projection plane of input grid
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe1)
    !     plov1:       orientation of input grid, as defined by gribcode
    !     pthe2:       latitude of projection plane of output grid
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe2)
    !     plov2:       orientation of output grid, as defined by gribcode
    !     
    !     output:
    !     pxnew,pynew: arrays of  x and y coordinates in output grid
    !     pdiri:       array of incremented orientation angles
    !     
    !     coordinates in metres, pthe and plov in degrees
    !     output is in degrees
    !     
    !     refer to gb2sla for comments concerning the use of pdiri
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     g.j.cats  03 jan 89.
    !     
    implicit none

    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    logical lsamhs
    real(kind=realgribkind) zrad,zrho,pthe1,pthe2,zlac1,zlac2,zrot,plov2,plov1,zcos,zsin
    real(kind=realgribkind) zscale,zx,zy,zr
    integer j,k

    data zrad/57.2957795_realgribkind   /
    data zrho/ 0.0174532925_realgribkind/
    lsamhs=pthe1*pthe2>0._realgribkind.or.&
         (abs(pthe1)<1.e-14_realgribkind.and.abs(pthe2)<1.e-14_realgribkind)
    zlac1=1._realgribkind+abs(sin(pthe1*zrho))
    zlac2=1._realgribkind+abs(sin(pthe2*zrho))
    if(lsamhs)then
       zrot=plov2-plov1
       if(pthe1<0._realgribkind)zrot=-zrot
       zcos=cos(zrot*zrho)
       zsin=sin(zrot*zrho)
       zscale=zlac2/zlac1
       zcos=zcos*zscale
       zsin=zsin*zscale
       do 100 j=1,k
          zx=pxold(j)
          zy=pyold(j)
          pxnew(j)=zx*zcos+zy*zsin
          pynew(j)=zy*zcos-zx*zsin
          pdiri(j)=pdiri(j)+zrot
100    enddo
    else
       zscale=-zlac2*zlac1*6371000._realgribkind**2._realgribkind
       zrot=plov2-plov1+180._realgribkind
       if(pthe1<0._realgribkind)zrot=-zrot
       zcos=cos(zrot*zrho)
       zsin=sin(zrot*zrho)
       do 101 j=1,k
          zx=pxold(j)
          zy=pyold(j)
          zr=zscale/(zx*zx+zy*zy)
          pdiri(j)=pdiri(j)+2._realgribkind*atan2(zy,zx)*zrad-zrot
          pxnew(j)=zr*(zx*zcos+zy*zsin)
          pynew(j)=zr*(zx*zsin-zy*zcos)
101    enddo
    endif
    return
  end subroutine gb2ppa
  subroutine gb2ppc(pxold,pyold,pxnew,pynew,k,pthe1,plov1,pthe2,plov2)
    !     
    !     subroutine gb2ppc: transform from polar-stereographic to polar-s
    !     
    !     input:
    !     pxold,pyold: arrays of x and y coordinates to be transformed
    !     k:           number of gridpoints to be transformed
    !     pthe1:       latitude of projection plane of input grid
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe1)
    !     plov1:       orientation of input grid, as defined by gribcode
    !     pthe2:       latitude of projection plane of output grid
    !     (it is assumed that the projection centre is the pole
    !     of which the latitude has opposite sign to pthe2)
    !     plov2:       orientation of output grid, as defined by gribcode
    !     
    !     output:
    !     pxnew,pynew: arrays of  x and y coordinates in output grid
    !     
    !     coordinates in metres, pthe and plov in degrees
    !     output is in degrees
    !     
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     g.j.cats  03 jan 89.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*)
    logical lsamhs
    real(kind=realgribkind) zrho,pthe1,pthe2,zlac1,zlac2,zrot,plov2,plov1
    real(kind=realgribkind) zcos,zsin,zscale,zx,zy,zr
    integer j,k


    data zrho/ 0.0174532925_realgribkind/
    lsamhs=pthe1*pthe2>0._realgribkind.or. &
         (abs(pthe1)<1.e-14_realgribkind.and.abs(pthe2)<1.e-14_realgribkind)
    zlac1=1._realgribkind+abs(sin(pthe1*zrho))
    zlac2=1._realgribkind+abs(sin(pthe2*zrho))
    if(lsamhs)then
       zrot=plov2-plov1
       if(pthe1<0._realgribkind)zrot=-zrot
       zcos=cos(zrot*zrho)
       zsin=sin(zrot*zrho)
       zscale=zlac2/zlac1
       zcos=zcos*zscale
       zsin=zsin*zscale
       do 100 j=1,k
          zx=pxold(j)
          zy=pyold(j)
          pxnew(j)=zx*zcos+zy*zsin
          pynew(j)=zy*zcos-zx*zsin
100    enddo
    else
       zscale=-zlac2*zlac1*6371000._realgribkind**2._realgribkind
       zrot=plov2-plov1+180._realgribkind
       if(pthe1<0._realgribkind)zrot=-zrot
       zcos=cos(zrot*zrho)
       zsin=sin(zrot*zrho)
       do 101 j=1,k
          zx=pxold(j)
          zy=pyold(j)
          zr=zscale/(zx*zx+zy*zy)
          pxnew(j)=zr*(zx*zcos+zy*zsin)
          pynew(j)=zr*(zx*zsin-zy*zcos)
101    enddo
    endif
    return
  end subroutine gb2ppc
  subroutine gb2r2v(pu,pv,pdiri)
    !     
    !     subroutine gb2r2v: change vector components by rotation over pdiri
    !     
    !     input:
    !     pu,pv: arrays of vector components to be changed
    !     pdiri: array of rotation angles (see gb2sla for sign)
    !     k:     number of gridpoints to be transformed
    !     
    !     output:
    !     pu,pv: changed vector components
    !     
    !     pdiri in degrees, pu and pv dimensions are immaterial
    !     
    !     g.j.cats  12 jan 89.
    !     
    implicit none
    real(kind=realgribkind) pu,pv,pdiri
    real(kind=realgribkind) zrho,zs,zc,z

    data zrho/ 0.0174532925_realgribkind/
    zs=sin(pdiri*zrho)
    zc=cos(pdiri*zrho)
    z=pu
    pu= zc*z+zs*pv
    pv=-zs*z+zc*pv
    return
  end subroutine gb2r2v
  subroutine gb2rlc(pxold,pxnew,k,plam)
    !
    !     subroutine gb2rlc: rotate regular lat/lon coordinates about axis
    !
    !     input:
    !     pxold:       arrays of longitudes to be transformed
    !     k:           number of gridpoints to be transformed
    !     plam:        rotation angle
    !
    !     output:
    !     pxnew:       arrays of transformed longitudes
    !
    !     all parameters to be in degrees,
    !     output is between -180(excl) and 180. (incl) (provided that
    !     pxold was between these limits and -360<=plam<=360.).
    !     in and output arrays may be the same
    !
    !                                                 g.j.cats  08 nov 88.
    !
    implicit none
    real(kind=realgribkind) pxold(*),pxnew(*)
    real(kind=realgribkind) plam
    integer j,k

    if(plam>0._realgribkind)then
       do 100 j=1,k
          pxnew(j)=pxold(j)-plam
          if(pxnew(j)<=-180._realgribkind)pxnew(j)=pxnew(j)+360._realgribkind
100    enddo
    else
       do 101 j=1,k
          pxnew(j)=pxold(j)-plam
          if(pxnew(j)> 180._realgribkind)pxnew(j)=pxnew(j)-360._realgribkind
101    enddo
    endif
    return
  end subroutine gb2rlc
  subroutine gb2sla(pxold,pyold,pxnew,pynew,pdiri,k,pthe)
    !     
    !     subroutine gb2sla: as gb2slc, but also calculate angle of rotation
    !     
    !     input:
    !     pxold,pyold: arrays of longitudes,latitudes to be transformed
    !     pdiri:       array of initial orientation angles
    !     k:           number of gridpoints to be transformed
    !     pthe:        latitude of new south pole, measured in old grid
    !     
    !     output:
    !     pxnew,pynew: arrays of transformed longitudes,latitudes
    !     pdiri:       array of incremented orientation angles
    !     
    !     all parameters to be in degrees,
    !     output is between -180(excl) and 180. (incl)
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     this subroutine does the same as gb2slc, but also
    !     calculates changes to the angles pdiri due to the
    !     angle between the coordinate axes.  use of gb2slc
    !     if this angle is not required, is slightly more
    !     efficient, but if the angle is required, use
    !     gb2sla also for the coordinate transformation itself
    !     as use of gb2sla is much more efficient than use
    !     of gb2slc followed by use of gb2sla.
    !     
    !     note that this routine calculates the changes to
    !     the angles as follows:
    !     pdiri on output=pdiri on input + angle between
    !     coordinate axes
    !     
    !     if you are interested in the angle between the
    !     coordinate axes only, ensure that pdiri is zero
    !     on input.
    !     
    !     the sign of the increment of pdiri is such that it
    !     denotes the angle of rotation to be applied to
    !     the old coordinate system to let its axes
    !     coincide with the new one, where the angle is
    !     measured in the old coordinate system, so pdiri is positive
    !     for anti-clockwise rotation from old x-axis to new one. as a
    !     consequence, if the initial value of pdiri is
    !     a wind direction dd (measured conform meteorological use)
    !     in the old coordinate system, the output value of
    !     pdiri will be dd in the new coordinate system
    !     (modulo 360). note that the meteorological use of dd
    !     has effectively an opposite sign to the mathematically
    !     usual way of defining angles.
    !     
    !     g.j.cats  15 dec 88.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*),pdiri(*)
    real(kind=realgribkind) zrad,zrho,zst,zct,pthe,zlat,zlon,zsinla,zcosla
    real(kind=realgribkind) zcoslo,zsinlo,zarg,zalph,zdiv,zxn,zdx
    integer j,k

    data zrad/57.2957795_realgribkind   /
    data zrho/ 0.0174532925_realgribkind/
    zst=-sin(pthe*zrho)
    zct=-cos(pthe*zrho)
    do 100 j=1,k
       zlat=pyold(j)*zrho
       zlon=pxold(j)
       if(zlon>180._realgribkind)zlon=zlon-360._realgribkind
       zsinla=sin(zlat)
       zcosla=cos(zlat)
       zcoslo=cos(zlon*zrho)
       zsinlo=sin(zlon*zrho)
       zarg=zst*zsinla+zct*zcosla*zcoslo
       zarg=max(min(zarg,1._realgribkind),-1._realgribkind)
       zalph=asin(zarg)
       zdiv=cos(zalph)
       if(abs(zdiv)<1.e-14_realgribkind)zdiv=1._realgribkind
       zdiv=1._realgribkind/zdiv
       pynew(j)=zalph*zrad
       zarg=(zst*zcosla*zcoslo-zct*zsinla)*zdiv
       zarg=max(min(zarg,1._realgribkind),-1._realgribkind)
       zxn=sign(acos(zarg),zlon*zcosla)
       pxnew(j)=zxn*zrad
       zarg=zst*zsinlo*sin(zxn)+zcoslo*zarg
       zarg=max(min(zarg,1._realgribkind),-1._realgribkind)
       zdx=sign(acos(zarg),-zlon*zct)
       pdiri(j)=pdiri(j)-zdx*zrad
100 enddo
    return
  end subroutine gb2sla
  subroutine gb2slc(pxold,pyold,pxnew,pynew,k,pthe)
    !     
    !     subroutine gb2slc: shift pole of regular lat/lon coordinate system
    !     
    !     input:
    !     pxold,pyold: arrays of longitudes,latitudes to be transformed
    !     k:           number of gridpoints to be transformed
    !     pthe:        latitude of new south pole, measured in old grid
    !     
    !     output:
    !     pxnew,pynew: arrays of transformed longitudes,latitudes
    !     
    !     all parameters to be in degrees,
    !     output is between -180(excl) and 180. (incl)
    !     
    !     in- and output arrays of coordinates may coincide
    !     
    !     g.j.cats  07 nov 88.
    !     
    implicit none
    real(kind=realgribkind) pxold(*),pyold(*),pxnew(*),pynew(*)
    real(kind=realgribkind) zrad,zrho,zst,zct,zlat,zlon,pthe,zsinla,zcosla
    real(kind=realgribkind) zcoslo,zarg,zalph,zdiv
    integer j,k

    data zrad/57.2957795_realgribkind   /
    data zrho/ 0.0174532925_realgribkind/
    zst=-sin(pthe*zrho)
    zct=-cos(pthe*zrho)
    do 100 j=1,k
       zlat=pyold(j)*zrho
       zlon=pxold(j)
       if(zlon>180._realgribkind)zlon=zlon-360._realgribkind
       zsinla=sin(zlat)
       zcosla=cos(zlat)
       zcoslo=cos(zlon*zrho)
       zarg=zst*zsinla+zct*zcosla*zcoslo
       zarg=max(min(zarg,1._realgribkind),-1._realgribkind)
       zalph=asin(zarg)
       zdiv=cos(zalph)
       if(abs(zdiv)<1.e-14_realgribkind)zdiv=1._realgribkind
       zdiv=1._realgribkind/zdiv
       pynew(j)=zalph*zrad
       zarg=(zst*zcosla*zcoslo-zct*zsinla)*zdiv
       zarg=max(min(zarg,1._realgribkind),-1._realgribkind)
       pxnew(j)=sign(acos(zarg)*zrad,zlon*zcosla)
100 enddo
    return
  end subroutine gb2slc
  subroutine gb2tlc(pyold,pynew,k,ps)
    !     
    !     subroutine gb2tlc: stretch regular lat/lon coords wrt. north pole
    !     
    !     input:
    !     pyold:       array of latitudes to be stretched
    !     k:           number of gridpoints to be transformed
    !     ps:          stretching factor
    !     
    !     output:
    !     pynew:       array of transformed latitudes
    !     
    !     all parameters to be in degrees,
    !     output is between -90(incl) and 90. (incl)
    !     in and output arrays may be the same
    !     note that de-stretching is achieved by using 1/ps as input
    !     
    !     g.j.cats  11 nov 88.
    !     
    implicit none

    real(kind=realgribkind) pyold(*),pynew(*)
    real(kind=realgribkind) zrad,zrho,zm,zp,ps,zs
    integer j,k

    data zrad/57.2957795_realgribkind  /
    data zrho/ 0.0174532925_realgribkind/
    zm=1._realgribkind-ps*ps
    zp=1._realgribkind+ps*ps
    do 100 j=1,k
       zs=sin(zrho*pyold(j))
       pynew(j)=zrad*asin((zm+zp*zs)/(zp+zm*zs))
100 enddo
    return
  end subroutine gb2tlc
end module vari
