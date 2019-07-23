module grw1
  use decomp,only:realkind
  use mod_grib,only:realgribkind
  implicit none
  private
  integer,parameter::klb1=30
  integer,parameter::knfldx=8000
  integer,parameter::klevx=100
  integer,parameter::klb2=400
  
  public gread !,priddr

contains

  subroutine gread(ludir,ktyp,kwmo,alev,field,mx,my,kerr)
    !     INTERFACE TO AS2DDR AND ASIMOF
    use asimof
#ifdef DEBUG
    use mod_grib
#endif

    implicit none

    integer,intent(in):: ludir,ktyp,kwmo,mx,my
    integer kerr
    real(kind=realgribkind)    alev
    real(kind=realkind)::field(mx*my)
    real(kind=realgribkind)::field4(mx*my)
    integer,parameter::klpsup=2
    integer kdev,fieldSize,kpf
    integer kb1(klb1),kb2(klb2)
    real(kind=realgribkind)    pb2(klb2),psup(klpsup)
#ifdef DEBUG
    real(kind=realgribkind):: zmax,zmin
#endif
    integer ibpw,ibpwio

    integer i,jpar,jtyp,jlev,readSize

    data    ibpw /0/, ibpwio /0/

    kdev = ludir
    fieldSize = mx*my
    kpf = 1
    call asimhm(ibpw,ibpwio,notdef)

    !     SET KEYS

    do 100 i=1,klb1
       kb1(i)=notdef
100 enddo

    jpar=kwmo


    !     parameter table version number 138 if jpar>999
    if(jpar>=1000)then
       jpar=jpar-1000
       kb1(4)=138
    endif
    if(ktyp==109.or.ktyp==0)then
       jlev=nint(alev)
       jtyp=109
    else if(ktyp==100)then
       jlev=nint(alev*0.01_realgribkind)
       jtyp=ktyp
    else
       jlev=nint(alev)
       jtyp=ktyp
    endif
    kb1(9)=jpar
    kb1(10)=jtyp
    kb1(11)=jlev

    field4 = real(field,realgribkind)

    call getfd(kdev,kb1,klb1,kb2,klb2,pb2,klb2,psup,klpsup,&
         field4,fieldSize,kpf,kerr,readSize)
    field = real(field4,realkind)

    if(kerr>=0)then
#ifdef DEBUG
       call maxmin(field4,readSize,zmax,zmin,real(notdef,realgribkind))
       write(6,*)'kerr=',kerr
       write(6,900)'read:par,typ,lev,max,min',jpar,jtyp,jlev,zmax,zmin
#endif
    else
       print *,'ERROR READING kerr=',kerr,jpar,jtyp,jlev,ludir, &
            __FILE__,__LINE__
    endif

900 format(a,3i4,2e12.4)
  end subroutine gread




  subroutine gwclos(ludir)
    use asimof
    implicit none
    integer ludir
    integer kerr
    call asimhc(ludir,kerr)
    return
  end subroutine gwclos




end module grw1
