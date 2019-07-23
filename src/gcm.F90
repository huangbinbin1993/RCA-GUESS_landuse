module gcm
  use domainmod
  use timetype
  use decomp
  use util
  use mod_grib,only:realgribkind

  implicit none
  private

  logical,public,save:: lecmwf    = .false. !era era-interim
  logical,public,save:: lhadley   = .false.        
  logical,public,save:: lipsl     = .false.          
  logical,public,save:: lecham    = .false.          
  logical,public,save:: lecham2   = .false.          
  logical,public,save:: lecham5   = .false.          
  logical,public,save:: lccsm     = .false.   !ccsm3
  logical,public,save:: lecearth  = .false.          
  logical,public,save:: lcanesm2  = .false.          
  logical,public,save:: lcnrm     = .false.          
  logical,public,save:: lnoresm   = .false.          
  logical,public,save:: lhadgem   = .false.          
  logical,public,save:: lmiroc5   = .false.          
  logical,public,save:: lmpiesm   = .false.          
  logical,public,save:: lgfdl     = .false.          
  logical,public,save:: lsstpath  = .false.          

!!$  type,public::gcmfile
!!$     character(len=128)::path
!!$     integer::nfields=0
!!$     integer,allocatable::typ(:)
!!$     integer,allocatable::lev(:) !levels, 1-60
!!$     integer,allocatable::par(:) 
!!$     integer,allocatable::nlev(:) !number of levels
!!$     integer,allocatable::special(:) !code for special treatment
!!$     integer::typeOfFile=-666 !many definitions to be made ask Ulf -> vilken subrutin ska läsa?
!!$  end type gcmfile
!!$
!!$  type,public::gcmfiles
!!$     integer::nfiles=0
!!$     type(gcmfile),allocatable::gcmfile(:)
!!$  end type gcmfiles
!!!  type(gcmfiles),public,save::filesGCM


  type(domain),public,save::gcmdomain
  
  public getgcm
  public readboundpath
  public readsstpath
  public grrdloc_hint
  public groploc
  public grrdloc_hint_lsmmatch
  public grclloc
  !public cre_file_ma
  public cre_filnam

  public createFileName

contains

  subroutine createFileName(prefix,ilen,iy,im,id,ih,length,filenam)
    implicit none
    integer :: ilen, iy, im, id, ih,  length  
    character (len=132) :: prefix, filenam  
    if(lecmwf.or.lccsm)then
       call cre_file_ma(prefix,ilen,iy,im,id,ih,length,filenam)
    elseif(lhadley.or.lecham.or.lecham2.or.lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl)then
       call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
    endif
  end subroutine createFileName
  
  subroutine cre_filnam(pref,ilen,iy,im,id,ih,length,filnam)
    !climate
    !gcm
    !lateral_bc
    !surface_bc
    implicit none  
    integer :: ilen, iy, im, id, ih,  length  
    character (len=132) :: pref, filnam  
    filnam = pref(1:ilen)  
    if(length==0)then  
       filnam (ilen + 1:ilen + 1) = '_'  
       write (filnam (ilen + 2:ilen + 5) , '(I4.4)') iy  
       write (filnam (ilen + 6:ilen + 7) , '(I2.2)') im  
       write (filnam (ilen + 8:ilen + 9) , '(I2.2)') id  
       write (filnam (ilen + 10:ilen + 11) , '(I2.2)') ih  
       filnam (ilen + 12:ilen + 21) = '00+000H00M'  
       !   climate
    else	  
       write (filnam (ilen + 1:ilen + 2) , '(I2.2)') iy  
       write (filnam (ilen + 3:ilen + 4) , '(I2.2)') im  
       write (filnam (ilen + 5:ilen + 6) , '(I2.2)') id  
       write (filnam (ilen + 7:ilen + 8) , '(I2.2)') ih  
    endif
    return  
  end subroutine cre_filnam


  subroutine cre_file_ma(pref,ilen,iy,im,id,ih,length,filnam)
    !gcm
    !lateral_bc
    !surface_bc
    implicit none  
    integer :: ilen, iy, im, id, ih, length  
    character (len=132) :: pref, filnam  
    filnam = pref (1:ilen)  
    if (length==0) then  
       write (filnam (ilen + 1:ilen + 4) , '(I4.4)') iy  
       write (filnam (ilen + 5:ilen + 6) , '(I2.2)') im  
       write (filnam (ilen + 7:ilen + 8) , '(I2.2)') id  
       write (filnam (ilen + 9:ilen + 10) , '(I2.2)') ih  
       filnam (ilen + 11:ilen + 13) = '00m'  
       !   climate
    else	  
       write (filnam (ilen + 1:ilen + 2) , '(I2.2)') iy  
       write (filnam (ilen + 3:ilen + 4) , '(I2.2)') im  
       write (filnam (ilen + 5:ilen + 6) , '(I2.2)') id  
       write (filnam (ilen + 7:ilen + 8) , '(I2.2)') ih  
    endif
    return  
  end subroutine cre_file_ma


  subroutine grclloc(lun, mype)  
    use grw1
    implicit none  
    integer :: lun, mype  
    if (mype==0) then  
       call grclos(lun)  
    endif
    return  
  end subroutine grclloc

  

  subroutine init_gcm(ktime,klon,klat,RCAdom)
    use referenceParameters, only:sip0
    use modddr
    use confys
    implicit none
    type(time),intent(in)::ktime
    integer,intent(in)::klon,klat
    type(domain),intent(in)::RCAdom
    integer::iy,im,id,ih,imin,length,ilen,lundir
    character(len=132):: prefix,filenam 
    integer::ierr
    integer, parameter::jpbuf=700000
    real(kind=realgribkind):: buf(jpbuf)
    type(ddr)::gcmddr
    real(kind=realkind):: zeps=1.0e-6_realkind
    logical::lfis
    integer::i,j, ilev
    real(kind=realkind)::dlon_gribi,dlat_gribi
    integer,allocatable:: iind(:,:),jind(:,:)
    real(kind=realkind),allocatable:: wx(:,:),wy(:,:)
    real(kind=realkind),allocatable,dimension(:,:)::zlon2,zlat2,lat,lon,zlon1,zlat1,w_x,w_y
    real(kind=realkind)::x,y,eps, rmin,rmax,sav
    logical::rotate

    integer::lun_orog, lenbuf, itype

    eps = 1.0e-6_realkind

    lfis=.true.
    if(lecmwf)then
       gcmdomain%arakawa = 'a'
    elseif(lhadley)then
       gcmdomain%arakawa = 'b'
       lfis=.false.
    elseif(lipsl)then
       gcmdomain%arakawa = 'a'
    elseif(lecham2)then
       gcmdomain%arakawa = 'a'
    elseif(lecham5)then
       gcmdomain%arakawa = 'c'
    elseif(lccsm)then
       gcmdomain%arakawa = 'a' 
    elseif(lecearth)then
       gcmdomain%arakawa = 'c'
    elseif(lcanesm2)then
       gcmdomain%arakawa = 'a'
    elseif(lcnrm)then
       gcmdomain%arakawa = 'a'
    elseif(lnoresm)then
       gcmdomain%arakawa = 'a'
    elseif(lhadgem)then
       gcmdomain%arakawa = 'a'
    elseif(lmiroc5)then
       gcmdomain%arakawa = 'a'
    elseif(lmpiesm)then
       gcmdomain%arakawa = 'a'
    elseif(lgfdl)then
       gcmdomain%arakawa = 'a'
    else
       if(mype==0)then
          write(6,*)'You have not defined Arakawa grdi for lateral boundaries'
          write(6,*)'Have to stop'
          call stop_program('init_gcm Arakawa?')
       endif
    endif
    call arakawastag(gcmdomain%arakawa,gcmdomain%stagu_x,gcmdomain%stagu_y, &
         gcmdomain%stagv_x,gcmdomain%stagv_y)
    call readboundpath(prefix,ilen)
    iy=ktime%year
    im=ktime%month
    id=ktime%day
    ih=0 !ktime%hour
    imin=0
    length=0


    call createFileName(prefix,ilen,iy,im,id,ih,length,filenam)

    lundir = 91
    call groploc(lundir,filenam,buf,jpbuf,gcmddr)
    gcmdomain%nlon = gcmddr%nlonhl
    gcmdomain%nlat = gcmddr%nlathl
    gcmdomain%nlev = gcmddr%nlevhl

    if(lhadley)then
       gcmddr%aweshl = 0._realkind    
       gcmddr%alafhl = -90._realkind  
       gcmddr%dlonhl = 3.75_realkind  
       gcmddr%dlathl = 2.50_realkind  
       gcmddr%aplohl = 0._realkind    
       gcmddr%aplahl = -90._realkind  
       gcmdomain%nlev=19
    endif
    
    gcmdomain%west   = real(gcmddr%aweshl,realkind)
    gcmdomain%south  = real(gcmddr%alafhl,realkind)
    gcmdomain%dlat   = real(gcmddr%dlathl,realkind)
    gcmdomain%dlon   = real(gcmddr%dlonhl,realkind)
    gcmdomain%polon  = real(gcmddr%aplohl,realkind)
    gcmdomain%polat  = real(gcmddr%aplahl,realkind)

    if( abs(gcmdomain%polon)<zeps .and.abs(gcmdomain%polat)<zeps ) then
       gcmdomain%polon = 0._realkind
       gcmdomain%polat = -90._realkind
    endif

    if(lccsm)then
       gcmdomain%south = -gcmdomain%south
       gcmdomain%dlat = -gcmdomain%dlat
    endif

    allocate(gcmdomain%afull(gcmdomain%nlev))
    allocate(gcmdomain%bfull(gcmdomain%nlev))
    allocate(gcmdomain%ahyb(gcmdomain%nlev+1))
    allocate(gcmdomain%bhyb(gcmdomain%nlev+1))

    if(lhadley)then
       if(gcmdomain%nlev/=19)call stop_program( 'number of levels is wrong for hadley')
       gcmdomain%afull = (/460.588_realkind,1479.717_realkind,2959.433_realkind,5529.418_realkind, &
            8861.730_realkind,11801.445_realkind,13660.086_realkind,14577.727_realkind,          &
            14688.062_realkind,14006.977_realkind,12323.531_realkind,9469.934_realkind,          &
            5809.324_realkind,2408.037_realkind,472.518_realkind, 0._realkind,0._realkind,0._realkind,0._realkind/)
       gcmdomain%bfull = (/0._realkind,0._realkind,0._realkind,0.002_realkind,0.011_realkind,0.031_realkind, &
            0.063_realkind,0.104_realkind,0.153_realkind,0.215_realkind,0.299_realkind,0.410_realkind,    &
            0.541_realkind,0.675_realkind,0.788_realkind,0.870_realkind,0.930_realkind,0.975_realkind,0.997_realkind/)
    elseif(lccsm)then
       if(gcmdomain%nlev/=26)call stop_program( 'number of levels is wrong for ccsm')
       gcmdomain%afull = (/354.4638_realkind, 738.88135_realkind, 1396.7214_realkind, 2394.4625_realkind,     &
            3723.029_realkind, 5311.4605_realkind, 7005.915_realkind, 7791.257_realkind, 7660.701_realkind, 7507.1085_realkind,&
            7326.415_realkind, 7113.8385_realkind, 6863.7535_realkind, 6569.5415_realkind, 6223.4155_realkind,        &
            5816.2165_realkind, 5337.168_realkind, 4773.5925_realkind, 4110.5755_realkind, 3330.57_realkind,          &
            2496.844_realkind, 1709.591_realkind, 1021.471_realkind, 480.3175_realkind, 126.068_realkind, 0._realkind /)
       gcmdomain%bfull = (/0._realkind,0._realkind,0._realkind,0._realkind,0._realkind,0._realkind,0._realkind,&
            0.00752654500000002_realkind, &
            0.023907685_realkind,      &
            0.04317925_realkind,0.0658512450000003_realkind, 0.0925236850000004_realkind,        &
            0.1239024_realkind, 0.16081785_realkind, 0.204247_realkind, 0.2553391_realkind,               &
            0.315446300000001_realkind, 0.386159300000001_realkind, 0.469349500000002_realkind,  &
            0.567218500000003_realkind, 0.671827850000003_realkind, 0.770606150000003_realkind,  &
            0.856946050000001_realkind, 0.924845700000002_realkind, 0.969294150000001_realkind,  &
            0.9925561_realkind/)
       elseif(lcanesm2) then
          if(gcmdomain%nlev/=35)call stop_program( 'number of levels is wrong for canesm2')
       gcmdomain%afull = (/25272.3577727806_realkind,70281.864878971_realkind,139062.283468519_realkind, &
          206330.475534117_realkind,272075.42524483_realkind,336285.872464643_realkind,398950.303586331_realkind, &
          460056.941876517_realkind,532146.203349759_realkind,613963.547365183_realkind,707822.401086833_realkind, &
          811280.602202879_realkind,924712.62365136_realkind,1043698.05111908_realkind,1165514.53124781_realkind, &
          1282507.53599184_realkind,1385435.2331673_realkind,1462196.76078162_realkind,1501144.89552048_realkind, &
          1491369.55106738_realkind,1424310.24627159_realkind,1296282.12390124_realkind,1145805.20507116_realkind, &
          991270.658409403_realkind,844603.04655789_realkind,716076.828212889_realkind,594882.021789457_realkind, &
          492459.933162568_realkind,408338.997197718_realkind,296913.100872_realkind,198800.309802059_realkind, &
          117534.533316727_realkind,57728.5093401557_realkind,25834.3812810661_realkind,10212.1412543364_realkind/)
       gcmdomain%bfull = (/0.992505689126255_realkind, 0.979063376936541_realkind, 0.958274942413293_realkind, & 
          0.937635760409187_realkind, 0.91714691815586_realkind, 0.896809526997173_realkind, & 
          0.876624723293888_realkind, 0.856593669376578_realkind, 0.832478661335397_realkind, & 
          0.804403518815122_realkind, 0.771139913039199_realkind, 0.732928878582424_realkind, & 
          0.68873345601546_realkind, 0.638989927840596_realkind, 0.582966982703533_realkind, & 
          0.521420101066735_realkind, 0.455261425861893_realkind, 0.386685278248952_realkind, & 
          0.317841206521863_realkind, 0.25080600561909_realkind, 0.188424571035177_realkind, & 
          0.133060587850253_realkind, 0.0939122379519188_realkind, 0.0661643645470388_realkind, & 
          0.0466400467274092_realkind, 0.0333252242190201_realkind, 0.0232868119039225_realkind, & 
          0.0163955849622416_realkind, 0.0116980855509555_realkind, 0.00669550919147256_realkind, & 
          0.00337896665988366_realkind, 0.00139967101098229_realkind, 0.000422357941161105_realkind, &
          0.000100218981339704_realkind, 1.2090282832963e-05_realkind/)
       elseif(lcnrm) then
          if(gcmdomain%nlev/=31)call stop_program( 'number of levels is wrong for cnrm')
       gcmdomain%afull = (/0.0_realkind, 0.0_realkind, 0.000355693896219267_realkind, &
          0.00169639707906107_realkind, 0.00449050044730997_realkind, &
          0.00905251144303739_realkind, 0.0155513083334619_realkind, 0.0240196175306947_realkind, & 
          0.0343646945181392_realkind, 0.0463802163669176_realkind, 0.0597593973654985_realkind, &
          0.0741093325826752_realkind, 0.0889665211623889_realkind, 0.103813606706565_realkind, &
          0.118097351616486_realkind, 0.131247840572557_realkind, 0.1426988601308_realkind, &
          0.151909512277887_realkind, 0.158387061944719_realkind, 0.161710883514314_realkind, &
          0.161557814813672_realkind, 0.157728610099951_realkind, 0.150175530602172_realkind, &
          0.13903134254689_realkind, 0.124639356336377_realkind, 0.107584776807749_realkind, &
          0.0887272296088351_realkind, 0.069101678183613_realkind, 0.0493583415597236_realkind, &
          0.0296150049358342_realkind, 0.00987166831194472_realkind/) 

       gcmdomain%bfull = (/0.996140748262405_realkind, 0.982633352279663_realkind, &
          0.958599209785461_realkind, 0.925964564085007_realkind, 0.886335849761963_realkind, &
          0.841061383485794_realkind, 0.791286081075668_realkind, 0.737999796867371_realkind, &
          0.682079493999481_realkind, 0.624324381351471_realkind, 0.565485388040543_realkind, &
          0.506287336349487_realkind, 0.447445198893547_realkind, 0.38967365026474_realkind, &
          0.333690196275711_realkind, 0.280212000012398_realkind, 0.229945950210094_realkind, &
          0.18357240408659_realkind, 0.141722552478313_realkind, 0.104949299246073_realkind, &
          0.0736912991851568_realkind, 0.0482312496751547_realkind, 0.0286470493301749_realkind, &
          0.0147566497325897_realkind, 0.00605690013617277_realkind, 0.00165530000231229_realkind, &
          0.000195450003957376_realkind, 0.0_realkind, 0.0_realkind, 0.0_realkind, 0.0_realkind/) 
       elseif(lnoresm) then
          if(gcmdomain%nlev/=26)call stop_program( 'number of levels is wrong for noresm')
       gcmdomain%afull = (/0.0_realkind, 0.00126068_realkind, 0.00480317500000001_realkind, &
          0.01021471_realkind, 0.01709591_realkind, 0.02496844_realkind, &
          0.0333057_realkind, 0.041105755_realkind, 0.0477359250000001_realkind, & 
          0.0533716800000001_realkind, 0.0581621650000002_realkind, 0.0622341550000001_realkind, & 
          0.065695415_realkind, 0.0686375349999999_realkind, 0.071138385_realkind, &
          0.0732641500000002_realkind, 0.0750710850000003_realkind, 0.0766070100000003_realkind, &
          0.0779125700000003_realkind, 0.0700591500000003_realkind, 0.0531146050000002_realkind, &
          0.0372302900000001_realkind, 0.023944625_realkind, 0.013967214_realkind, &
          0.00738881350000001_realkind, 0.00354463800000001_realkind/)
       gcmdomain%bfull = (/0.9925561_realkind, 0.969294150000001_realkind, 0.924845700000002_realkind, &
          0.856946050000001_realkind, 0.770606150000003_realkind, 0.671827850000003_realkind, &
          0.567218500000003_realkind, 0.469349500000002_realkind, 0.386159300000001_realkind, &
          0.315446300000001_realkind, 0.2553391_realkind, 0.204247_realkind, &
          0.16081785_realkind, 0.1239024_realkind, 0.0925236850000004_realkind, &
          0.0658512450000003_realkind, 0.04317925_realkind, 0.023907685_realkind, &
          0.00752654500000002_realkind, 0.0_realkind, 0.0_realkind, &
          0.0_realkind, 0.0_realkind, 0.0_realkind, &
          0.0_realkind, 0.0_realkind/)
       elseif(lhadgem) then
!      Note thet HadGEM is a non-hydrostatic model which means that a and b coefficients are given with
!      resepct to height coordinate instead of pressure, i.e. z(k)=a(k)+b(k)*fis/gravit.
!      For this reason there is special treatment of HadGEM in subroutine etaeta.
          if(gcmdomain%nlev/=38)call stop_program( 'number of levels is wrong for HadGEM')
       gcmdomain%afull = (/39254.83203125_realkind, 32908.69140625_realkind, 29219.080078125_realkind, 26583.640625_realkind, &
          24458.28515625_realkind, 22626.08203125_realkind, 20990.1875_realkind, 19503.009765625_realkind, &
          18138.626953125_realkind, 16875.310546875_realkind, 15694.6396484375_realkind, 14580.7998046875_realkind, &
          13520.0009765625_realkind, 12499.9990234375_realkind, 11519.998046875_realkind, 10579.998046875_realkind, &
          9679.9990234375_realkind, 8820._realkind, 8000.00146484375_realkind, 7220._realkind, &
          6479.99951171875_realkind, 5779.99951171875_realkind, 5120._realkind, 4500.00146484375_realkind, &
          3919.99951171875_realkind, 3379.99829101562_realkind, 2880.00146484375_realkind, 2420.00170898438_realkind, &
          1999.99841308594_realkind, 1619.99987792969_realkind, 1279.998046875_realkind, 980.000854492188_realkind, &
          720.000366210938_realkind, 500.000579833984_realkind, 320.00146484375_realkind, 179.999114990234_realkind, &
          80.001350402832_realkind, 20.000337600708_realkind/)
       gcmdomain%bfull = (/ 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, &
          0._realkind, 0._realkind, 0.00130179093685001_realkind, 0.0107164792716503_realkind, &
          0.0279368180781603_realkind, 0.0518637150526047_realkind, 0.0817952379584312_realkind, 0.116947874426842_realkind, &
          0.156554222106934_realkind, 0.199878215789795_realkind, 0.24621507525444_realkind, 0.294891387224197_realkind, &
          0.34526526927948_realkind, 0.39672577381134_realkind, 0.44869339466095_realkind, 0.500619947910309_realkind, &
          0.55198872089386_realkind, 0.602314412593842_realkind, 0.651142716407776_realkind, 0.698050200939178_realkind, &
          0.742646217346191_realkind, 0.784570515155792_realkind, 0.823493480682373_realkind, 0.859118342399597_realkind, &
          0.891178011894226_realkind, 0.919438362121582_realkind, 0.943695485591888_realkind, 0.9637770652771_realkind, &
          0.979542553424835_realkind, 0.990881502628326_realkind, 0.99771648645401_realkind/)
       gcmdomain%ahyb = (/42427.90234375_realkind, 36081.76171875_realkind, 31063.890625_realkind, 27901.359375_realkind, &
          25520.9609375_realkind, 23542.18359375_realkind, 21808.13671875_realkind, 20246.6015625_realkind, 18820.8203125_realkind, &
          17506.96875_realkind, 16284.97265625_realkind, 15137.71875_realkind, 14050.3984375_realkind, 13010._realkind, &
          12010._realkind, 11050._realkind, 10130._realkind, 9250._realkind, 8410._realkind, 7610._realkind, 6850._realkind, &
          6130._realkind, 5450._realkind, 4810._realkind, 4210._realkind, 3650.00073242188_realkind, 3129.99975585938_realkind, &
          2649.99951171875_realkind, 2210._realkind, 1810.0009765625_realkind, 1449.9990234375_realkind, 1130.00146484375_realkind, &
          850.00048828125_realkind, 610.00048828125_realkind, 410.0009765625_realkind, 249.998336791992_realkind, 130.000228881836_realkind, &
          49.9988861083984_realkind, 0._realkind/)
       gcmdomain%bhyb = (/0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, 0._realkind, &
          0._realkind, 0._realkind, 0.00487210974097252_realkind, 0.0183146893978119_realkind, 0.0389823913574219_realkind, &
          0.0659807920455933_realkind, 0.0985881090164185_realkind, 0.136030256748199_realkind, 0.177555441856384_realkind, &
          0.222443282604218_realkind, 0.270004868507385_realkind, 0.319582045078278_realkind, 0.370548844337463_realkind, &
          0.422309935092926_realkind, 0.474301397800446_realkind, 0.525990784168243_realkind, 0.576877355575562_realkind, 0.626490533351898_realkind, &
          0.674392521381378_realkind, 0.720175802707672_realkind, 0.763464510440826_realkind, 0.80391401052475_realkind, &
          0.84121161699295_realkind, 0.875074565410614_realkind, 0.905253052711487_realkind, 0.931527435779572_realkind, &
          0.953709840774536_realkind, 0.971644043922424_realkind, 0.985203862190247_realkind, 0.994296252727509_realkind, &
          1._realkind/)

       elseif(lmiroc5) then
!
! from /home/sm_grini/Scripts/CMIP5/lbc/MIROC5/vert_coor/MIROC5_a_and_b.txt
!
          if(gcmdomain%nlev/=40)call stop_program( 'number of levels is wrong for miroc5')
       gcmdomain%afull = (/ 0.00111274564266205_realkind, 0.00378219056129456_realkind, 0.0075644998550415_realkind, &
          0.012458571434021_realkind, 0.0184659862518311_realkind, 0.0258085975646973_realkind, &
          0.0344865684509277_realkind, 0.0442758674621582_realkind, 0.0551801223754883_realkind, &
          0.0671956329345703_realkind, 0.0801021499633789_realkind, 0.094804084777832_realkind, &
          0.112640205383301_realkind, 0.134059173583984_realkind, 0.159080047607422_realkind, &
          0.189771148681641_realkind, 0.222394271850586_realkind, 0.251527801513672_realkind, &
          0.277208648681641_realkind, 0.303096038818359_realkind, 0.28772509765625_realkind, &
          0.250411437988281_realkind, 0.217940170288086_realkind, 0.189681274414063_realkind, &
          0.165080352783203_realkind, 0.143672698974609_realkind, 0.125043327331543_realkind, &
          0.108827621459961_realkind, 0.0947157287597656_realkind, 0.0824342575073242_realkind, &
          0.0717448806762695_realkind, 0.0624412078857422_realkind, 0.0542777328491211_realkind, &
          0.0469373054504395_realkind, 0.0401577377319336_realkind, 0.0336126861572266_realkind, &
          0.0265931816101074_realkind, 0.0186132392883301_realkind, 0.0106025438308716_realkind, &
          0.00290465235710144_realkind  /)
       gcmdomain%bfull = (/ 0.996386528015137_realkind, 0.987716317176819_realkind, 0.975432455539703_realkind, &
          0.959537029266357_realkind, 0.940026998519897_realkind, 0.916181147098541_realkind, &
          0.88799923658371_realkind, 0.856206655502319_realkind, 0.82079690694809_realkind, 0.78177684545517_realkind, &
          0.739865183830261_realkind, 0.692146897315979_realkind, 0.634282648563385_realkind, &
          0.564825654029846_realkind, 0.483753234148026_realkind, 0.384700953960419_realkind, &
          0.279083371162415_realkind, 0.184922024607658_realkind, 0.102643676102161_realkind, &
          0.0274994932115078_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, &
          0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, &
          0_realkind, 0_realkind, 0_realkind, 0_realkind, &
          0.0_realkind, 0.0_realkind/)

       elseif(lipsl) then
!
! from /home/sm_grini/Scripts/CMIP5/lbc/IPSL-CMA5-LR/r1i1p1/vert_coor/IPSL-LR_a_and_b.txt
!

          if(gcmdomain%nlev/=39)call stop_program( 'number of levels is wrong for miroc5')
       gcmdomain%afull = (/  140.853424072266_realkind, 459.041107177734_realkind, 886.3955078125_realkind, 1492.93151855469_realkind, &
    2343.68768310547_realkind, 3497.71838378906_realkind, 5004.88647460938_realkind, 6898.5458984375_realkind, &
    9182.28002929688_realkind, 11809.6362304688_realkind, 14658.189453125_realkind, 17504.6948242188_realkind, &
    20017.6689453125_realkind, 21794.8330078125_realkind, 22472.3193359375_realkind, 21893.60546875_realkind, &
    20237.9482421875_realkind, 17956.2314453125_realkind, 15511.3125_realkind, 13171.06640625_realkind, &
    11026.7978515625_realkind, 9101.6884765625_realkind, 7402.44189453125_realkind, 5927.81982421875_realkind, &
    4670.0166015625_realkind, 3615.9208984375_realkind, 2748.46923828125_realkind, 2047.99407958984_realkind, &
    1493.48248291016_realkind, 1063.68096923828_realkind, 738.006469726562_realkind, 497.243316650391_realkind, &
    324.027206420898_realkind, 203.132019042969_realkind, 121.585121154785_realkind, 68.6437549591064_realkind, &
    35.6667098999023_realkind, 15.9143776893616_realkind, 4.30683565139771_realkind /)

       gcmdomain%bfull = (/ 0.994381487369537_realkind, 0.981732308864594_realkind, 0.964865922927856_realkind, &
    0.941147953271866_realkind, 0.908268451690674_realkind, 0.864292740821838_realkind, &
    0.807745635509491_realkind, 0.737765997648239_realkind, 0.654367953538895_realkind, &
    0.558830052614212_realkind, 0.454186111688614_realkind, 0.34568452835083_realkind, &
    0.240891009569168_realkind, 0.148888442665339_realkind, 0.0780392717570066_realkind, &
    0.0325538357719779_realkind, 0.00987227889709175_realkind, 0.00190454708354082_realkind, &
    0.000188917819286871_realkind, 6.67837648293812e-06_realkind, 4.47209505816115e-08_realkind, &
    1.96520117896045e-11_realkind, 9.28085928930852e-17_realkind, 1.90709342508476e-25_realkind, &
    4.81193280962963e-40_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, &
    0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind /)

       elseif(lgfdl) then
!
! from  /home/sm_grini/Scripts/CMIP5/lbc/GFDL-ESM2M/vert_coor/GFDL_ESM2M_a_and_b.txt
!

          if(gcmdomain%nlev/=24)call stop_program( 'number of levels is wrong for gfdl')
       gcmdomain%afull = (/    &
    0.0024054256263204_realkind, 0.0077853397319501_realkind, 0.0144052052491981_realkind,   &
    0.022457359421262_realkind, 0.0321076594497903_realkind, 0.0434867883358006_realkind,   &
    0.0566793646295954_realkind, 0.0717101240593388_realkind, 0.0885263038181594_realkind,   &
    0.106974733222304_realkind, 0.126771132108932_realkind, 0.147458219605848_realkind,   &
    0.168345901184308_realkind, 0.188422445565014_realkind, 0.206214937237849_realkind,   &
    0.219553263015051_realkind, 0.225317010547742_realkind, 0.22077037842339_realkind,   &
    0.202152957685665_realkind, 0.15740741811621_realkind, 0.10013717655132_realkind,   &
    0.0541838125308414_realkind, 0.0216049381823649_realkind, 0.00495162335143335_realkind   &
     /)

       gcmdomain%bfull = (/   &
    0.99367618560791_realkind, 0.979488730430603_realkind, 0.961919605731964_realkind,   &
    0.940375208854675_realkind, 0.914288520812988_realkind, 0.883127689361572_realkind,   &
    0.846406996250153_realkind, 0.80370044708252_realkind, 0.754659652709961_realkind,   &
    0.699036598205566_realkind, 0.636714458465576_realkind, 0.567749261856079_realkind,   &
    0.492430031299591_realkind, 0.411367237567902_realkind, 0.325632661581039_realkind,   &
    0.236995249986649_realkind, 0.151226207613945_realkind, 0.0768977031111717_realkind,   &
    0.0217839498072863_realkind,  &
    0_realkind, 0_realkind, 0_realkind, 0_realkind, 0_realkind /)



! RCA GRIB code provides mean values of half levels which is wrong       gcmdomain%afull = real(gcmddr%alevhl(1:gcmdomain%nlev,1),realkind)
! RCA GRIB code provides mean values of half levels which is wrong       gcmdomain%bfull = real(gcmddr%alevhl(1:gcmdomain%nlev,2),realkind)
    else
       gcmdomain%afull = real(gcmddr%alevhl(1:gcmdomain%nlev,1),realkind)
       gcmdomain%bfull = real(gcmddr%alevhl(1:gcmdomain%nlev,2),realkind)
    endif

    if(lnoresm.or.lmiroc5) then
       do ilev = 1,gcmdomain%nlev
          gcmdomain%afull(ilev)=gcmdomain%afull(ilev)*100000._realkind
       enddo
    endif
    if(lcnrm) then
       do ilev = 1,gcmdomain%nlev
          gcmdomain%afull(ilev)=gcmdomain%afull(ilev)*101300._realkind
       enddo
    endif
    if(lgfdl) then
       do ilev = 1,gcmdomain%nlev
          gcmdomain%afull(ilev)=gcmdomain%afull(ilev)*101325._realkind
       enddo
    endif
    if(lcanesm2) then
       do ilev = 1,gcmdomain%nlev
          gcmdomain%afull(ilev)=gcmdomain%afull(ilev)*.01_realkind
       enddo
    endif
    if(lcanesm2.or.lcnrm.or.lnoresm.or.lmiroc5.or.lipsl.or.lgfdl) then
!cau       upside down
       do ilev = 1,gcmdomain%nlev/2
          sav=gcmdomain%afull(ilev)
          gcmdomain%afull(ilev)=gcmdomain%afull(gcmdomain%nlev-ilev+1)
          gcmdomain%afull(gcmdomain%nlev-ilev+1)=sav
          sav=gcmdomain%bfull(ilev)
          gcmdomain%bfull(ilev)=gcmdomain%bfull(gcmdomain%nlev-ilev+1)
          gcmdomain%bfull(gcmdomain%nlev-ilev+1)=sav
       enddo
    endif
!!$
!!$    allocate(gcmdomain%fis(klon,klat))
!!$    !     read fis for the global model 
!!$    if(lfis) then
!!$       allocate(iind(klon,klat),jind(klon,klat))
!!$       allocate(wx(klon,klat),wy(klon,klat))
!!$       allocate(zlat2(klon,klat),zlon2(klon,klat))
!!$       allocate(lat(klon,klat),lon(klon,klat))
!!$       allocate(zlat1(klon,klat),zlon1(klon,klat))
!!$       allocate(w_x(klon,klat),w_y(klon,klat))
!!$       dlat_gribi = 1.0/gcmdomain%dlat
!!$       dlon_gribi = 1.0/gcmdomain%dlon
!!$
!!$       do j=1,klat
!!$          lat(:,j) = RCAdomain%south + real( jdatastart + j - 2 )*RCAdomain%dlat
!!$       enddo
!!$       do i=1,klon
!!$          lon(i,:) = RCAdomain%west + real( idatastart + i - 2 )*RCAdomain%dlon
!!$       enddo
!!$
!!$
!!$       if(abs(RCAdomain%polon-gcmdomain%polon)>eps.or. &
!!$            abs(RCAdomain%polat-gcmdomain%polat)>eps)then
!!$          rotate = .true.
!!$
!!$          !     de-rotate if needed coordinates from output grid geometry
!!$          if( abs(RCAdomain%polat+90.)>eps ) then
!!$             call regrot(zlon1,zlat1,lon,lat, &
!!$                  klon,klat,klon,klat,                  &
!!$                  RCAdomain%polon,RCAdomain%polat,-1)                        
!!$          else
!!$             do j=1,klat
!!$                do i=1,klon
!!$                   zlat1(i,j) = lat(i,j)
!!$                   zlon1(i,j) = lon(i,j)
!!$                enddo
!!$             enddo
!!$          endif
!!$
!!$          !     1.3.2  rotate if needed coordinates to input grid geometry 
!!$
!!$          if( abs(gcmdomain%polat+90.)>eps ) then
!!$             call  regrot(zlon1,zlat1,zlon2,zlat2, &
!!$                  klon,klat,klon,klat,                      &
!!$                  gcmdomain%polon,gcmdomain%polat,+1)                  
!!$          else
!!$             do j=1,klatruntime_canesm2.log
!!$                do i=1,klon
!!$                   zlat2(i,j) = zlat1(i,j)
!!$                   zlon2(i,j) = zlon1(i,j)
!!$                enddo
!!$             enddo
!!$          endif
!!$          !     no rotation needed, just copy coordinates
!!$       else
!!$          rotate = .false.
!!$          do j=1,klat
!!$             do i=1,klon
!!$                zlat2(i,j) = lat(i,j)
!!$                zlon2(i,j) = lon(i,j)
!!$             enddo
!!$          enddo
!!$       endif
!!$
!!$
!!$       ! calculate indeces and weights needed for horizontal interpolation
!!$
!!$       do j=1,klat
!!$          do i=1,klon
!!$             if( zlon2(i,j)<gcmdomain%west )then
!!$                zlon2(i,j) = zlon2(i,j) + 360.      
!!$             endif
!!$
!!$             y = ( zlat2(i,j) - gcmdomain%south )*dlat_gribi + 1.0
!!$             jind(i,j)    = int( y )
!!$             w_y(i,j)  = y - real(jind(i,j))
!!$
!!$             x = ( zlon2(i,j) - gcmdomain%west )*dlon_gribi + 1.0
!!$             iind(i,j)    = int( x )
!!$
!!$             if ( iind(i,j) > gcmdomain%nlon ) then
!!$                x = x - real(gcmdomain%nlon)
!!$                iind(i,j) = iind(i,j) - gcmdomain%nlon
!!$             endif
!!$             w_x(i,j)  = x - real(iind(i,j))
!!$
!!$          enddo
!!$       enddo
!!$
!!$       call grrdloc_hint(lundir,105,6,0.0,lecmwf, &
!!$            gcmdomain%fis,gcmdomain%fis,&
!!$            klon,klat,              &
!!$            gcmdomain%nlon,gcmdomain%nlat,                                   &
!!$            iind,wx,jind,wy,                               &
!!$            .false.,                                           &
!!$            iind,wx,jind,wy,                               &
!!$            ierr)
!!$
!!$       if (ierr/=0) then
!!$          print *, ' lecmwf ',lecmwf
!!$          print *, ' klon,klat = ',klon,klat
!!$          print *, ' nx_grib,ny_grib = ',gcmdomain%nlon,gcmdomain%nlat
!!$          print *,'error with fis GCM'
!!$          call stop_program('')
!!$       endif
!!$    endif

!cau110623    if(lhadley.or.lcanesm2) then
    if(lhadley) then

       allocate(gcmdomain%fis(klon,klat))
    !     read fis for the global model 
       allocate(iind(klon,klat),jind(klon,klat))
       allocate(wx(klon,klat),wy(klon,klat))
       allocate(zlat2(klon,klat),zlon2(klon,klat))
       allocate(lat(klon,klat),lon(klon,klat))
       allocate(zlat1(klon,klat),zlon1(klon,klat))
       allocate(w_x(klon,klat),w_y(klon,klat))
       dlat_gribi = 1.0_realkind/gcmdomain%dlat
       dlon_gribi = 1.0_realkind/gcmdomain%dlon


       do j=1,klat
          lat(:,j) = RCAdom%south + real( jdatastart + j - 2,realkind )*RCAdom%dlat
       enddo
       do i=1,klon
          lon(i,:) = RCAdom%west + real( idatastart + i - 2,realkind )*RCAdom%dlon
       enddo


       if(abs(RCAdom%polon-gcmdomain%polon)>eps.or. &
            abs(RCAdom%polat-gcmdomain%polat)>eps)then
          rotate = .true.

          !     de-rotate if needed coordinates from output grid geometry
          if( abs(RCAdom%polat+90._realkind)>eps ) then
             call regrot(zlon1,zlat1,lon,lat, &
                  klon,klat,&
                  RCAdom%polon,RCAdom%polat,-1)                        
          else
             do j=1,klat
                do i=1,klon
                   zlat1(i,j) = lat(i,j)
                   zlon1(i,j) = lon(i,j)
                enddo
             enddo
          endif

          !     1.3.2  rotate if needed coordinates to input grid geometry 

          if( abs(gcmdomain%polat+90._realkind)>eps ) then
             call  regrot(zlon1,zlat1,zlon2,zlat2, &
                  klon,klat,&
                  gcmdomain%polon,gcmdomain%polat,+1)                  
          else
             do j=1,klat
                do i=1,klon
                   zlat2(i,j) = zlat1(i,j)
                   zlon2(i,j) = zlon1(i,j)
                enddo
             enddo
          endif
          !     no rotation needed, just copy coordinates
       else
          rotate = .false.
          do j=1,klat
             do i=1,klon
                zlat2(i,j) = lat(i,j)
                zlon2(i,j) = lon(i,j)
             enddo
          enddo
       endif


       ! calculate indeces and weights needed for horizontal interpolation

       do j=1,klat
          do i=1,klon
             if( zlon2(i,j)<gcmdomain%west )then
                zlon2(i,j) = zlon2(i,j) + 360._realkind      
             endif

             y = ( zlat2(i,j) - gcmdomain%south )*dlat_gribi + 1.0_realkind
             jind(i,j)    = int( y )
             w_y(i,j)  = y - real(jind(i,j),realkind)

             x = ( zlon2(i,j) - gcmdomain%west )*dlon_gribi + 1.0_realkind
             iind(i,j)    = int( x )

             if ( iind(i,j) > gcmdomain%nlon ) then
                x = x - real(gcmdomain%nlon,realkind)
                iind(i,j) = iind(i,j) - gcmdomain%nlon
             endif
             w_x(i,j)  = x - real(iind(i,j),realkind)

          enddo
       enddo

       lun_orog = lundir + 4  
       if (mype.eq.0) then  
          call readorogpath(prefix,ilen)
          print *, ' orog prefix ',prefix
          if(lhadley) then
             prefix = prefix (1:ilen) //'HAD_OROG2'  
             ilen = ilen + 9

             iy = 0  
             im = 0  
             id = 0  
             ih = 0  
             length = 0
!!$             !!call cre_filnam(prefix,ilen,iy,im,id,ih,length,filenam)
             call createFileName(prefix,ilen,iy,im,id,ih,length,filenam)
             itype=105
          endif
          if (lcanesm2) then
             filenam = prefix (1:ilen) //'CanESM2_196001010600+000H00M'
             ilen = ilen + 28
             itype=1
          endif
       endif

       call groploc (lun_orog,filenam, buf, lenbuf, gcmddr)

       ilev=0

       call grrdloc_hint(lun_orog,itype,6,real(ilev,realgribkind),lecmwf, &
            gcmdomain%fis,gcmdomain%fis,&
            klon,klat,              &
            gcmdomain%nlon,gcmdomain%nlat,                 &
            iind,wx,jind,wy,                               &
            .false.,                                           &
            iind,wx,jind,wy,                               &
            ierr)

       rmin=1.e10
       rmax=0.
       do j=1,klat
          do i=1,klon
             rmin=min(rmin,gcmdomain%fis(i,j))
             rmax=max(rmax,gcmdomain%fis(i,j))
          enddo
       enddo

       if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl) then
          do j=1,klat
             do i=1,klon
                gcmdomain%fis(i,j) = gcmdomain%fis(i,j) * gravit
             enddo
          enddo
       endif

       if (ierr.ne.0) then
          print *, ' lecmwf ',lecmwf
          print *, ' klon,klat = ',klon,klat
          print *, ' nx_grib,ny_grib = ',gcmdomain%nlon,gcmdomain%nlat
          print *,'error with fis GCM'
          call stop_program('')
       endif
    endif

    call grclloc(lundir,mype)

  end subroutine init_gcm


  subroutine getgcm(ktime,klon,klat,RCAdom)
    implicit none
    type(time),intent(in)::ktime
    integer,intent(in)::klon,klat
    type(domain),intent(in)::RCAdom
    integer sum
    namelist/namgcm/lecmwf,lhadley,lipsl,lecham2,lecham5,lccsm,lecearth,lcanesm2,lcnrm,lnoresm,lhadgem,lmiroc5,lmpiesm,lgfdl,lsstpath

    open(57,file='namelists.dat',status='old')
    read(57,nml=namgcm)
    close(57)
    if(mype==0)then
       write(6,nml=namgcm)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namgcm)
       close(1)
    endif
    sum=0
    if(lecmwf)sum = sum+1
    if(lhadley)sum = sum+1
    if(lipsl)sum = sum+1
    if(lecham2)sum = sum+1
    if(lecham5)sum = sum+1
    if(lccsm)sum = sum+1
    if(lecearth)sum = sum+1
    if(lcanesm2)sum = sum+1
    if(lcnrm)sum = sum+1
    if(lnoresm)sum = sum+1
    if(lhadgem)sum = sum+1
    if(lmiroc5)sum = sum+1
    if(lmpiesm)sum = sum+1
    if(lgfdl)sum = sum+1
    if(sum/=1)then
       print *, 'Only one and exactly one GCM can be selected at a time!!!'
       print *, 'You chose ',sum
       print *,'lecmwf,lhadley,lipsl,lecham2,lecham5,lccsm,lecearth,lcanesm2,lcnrm,lnoresm,lhadgem,lmiroc5,lmpiesm,lgfdl=',&
            lecmwf,lhadley,lipsl,lecham2,lecham5,lccsm,lecearth,lcanesm2,lcnrm,lnoresm,lhadgem,lmiroc5,lmpiesm,lgfdl
       print *, 'Have to stop'
       call stop_program('')
    endif
    
    call init_gcm(ktime,klon,klat,RCAdom)
    


  end subroutine getgcm

!!$  subroutine init_gcm(gcmddrI)
!!$    implicit none
!!$    type(gcmddr)::gcmddrI
!!$    gcmddrI%nfields = numberOfFieldsIGCM()
!!$    allocate(gcmddrI%typ(nfields),&
!!$         gcmddrI%lev(nfields),&
!!$         gcmddrI%par(nfields),&
!!$         gcmddrI%filename(nfields))
!!$
!!$    if(lecmwf)then
!!$       gcmddrI%nfiles = 
!!$       gcmddrI%nlev = 60 !future? 
!!$       gcmddrI%typ = (\ \)
!!$       gcmddrI%lev = (\ \)
!!$       gcmddrI%par = (\ \)
!!$    elseif(lhadley)then
!!$       gcmddrI%typ = (\ \)
!!$       gcmddrI%lev = (\ \)
!!$       gcmddrI%par = (\ \)
!!$    elseif(lipsl)then
!!$       gcmddrI%typ = (\ \)
!!$       gcmddrI%lev = (\ \)
!!$       gcmddrI%par = (\ \)
!!$    elseif(lecham)then
!!$       if(lecham2)then
!!$          gcmddrI%typ = (\ \)
!!$          gcmddrI%lev = (\ \)
!!$          gcmddrI%par = (\ \)
!!$       elseif(lecham5)then
!!$          gcmddrI%typ = (\ \)
!!$          gcmddrI%lev = (\ \)
!!$          gcmddrI%par = (\ \)
!!$       else
!!$          gcmddrI%typ = (\ \)
!!$          gcmddrI%lev = (\ \)
!!$          gcmddrI%par = (\ \)
!!$       endif
!!$    elseif(lccsm)then
!!$       gcmddrI%typ = (\ \)
!!$       gcmddrI%lev = (\ \)
!!$       gcmddrI%par = (\ \)
!!$    elseif(lecearth)then
!!$       gcmddrI%typ = (\ \)
!!$       gcmddrI%lev = (\ \)
!!$       gcmddrI%par = (\ \)
!!$   endif
!!$  end subroutine init_gcm


!!$  integer function numberOfFieldsInGCM()
!!$    if(lecmwf)then
!!$       numberOfFieldsInGCM = 
!!$    elseif(lhadley)then
!!$       numberOfFieldsInGCM = 
!!$    elseif(lipsl)then
!!$       numberOfFieldsInGCM = 
!!$    elseif(lecham)then
!!$       if(lecham2)then
!!$          numberOfFieldsInGCM = 
!!$       elseif(lecham5)then
!!$          numberOfFieldsInGCM = 
!!$       else
!!$          numberOfFieldsInGCM =    
!!$       endif
!!$    elseif(lccsm)then
!!$       numberOfFieldsInGCM = 
!!$    elseif(lecearth)then
!!$       numberOfFieldsInGCM = 
!!$    endif
!!$  end function numberOfFieldsInGCM


  subroutine readboundpath(bound_path,ilen)
    implicit none
    character(len=132),intent(out)::bound_path
    integer,intent(out)::ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /bound/bound_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=bound)
    close(61)
    if(mype==0.and..not.done)then
       write(6,nml=bound)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=bound)
       close(1)
       done = .true.
    endif
    ilen=index(bound_path,' ')-1
  end subroutine readboundpath

  subroutine readsstpath(sst_path,ilen)
    implicit none
    character(len=132),intent(out)::sst_path
    integer,intent(out)::ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /sst/sst_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=sst)
    close(61)
    if(mype==0.and..not.done)then
       write(6,nml=sst)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=sst)
       close(1)
       done = .true.
    endif
    ilen=index(sst_path,' ')-1
  end subroutine readsstpath


  subroutine grrdloc_hint(lun,type,param,alev,lecmwfIn, &
       fFIELD,gFIELD, nx_local,ny_local, &
       nx_grib,ny_grib,  i_f,w_x_f,j_f,w_y_f, &
       lgfield,i_g,w_x_g,j_g,w_y_g,    ierr)
    use decomp
    use confys, only:tmelt
    use grw1

    implicit none

    integer lun,type,param,nx_local,ny_local,nx_grib,ny_grib
    logical lecmwfIn

    real(kind=realgribkind) alev
    real(kind=realkind) ::fFIELD(nx_local,ny_local),gFIELD(nx_local,ny_local)
    integer i_f(nx_local,ny_local),j_f(nx_local,ny_local)
    real(kind=realkind)  w_x_f(nx_local,ny_local),w_y_f(nx_local,ny_local)

    logical,intent(in):: lgfield
    integer i_g(nx_local,ny_local),j_g(nx_local,ny_local)
    real(kind=realkind)::  w_x_g(nx_local,ny_local),w_y_g(nx_local,ny_local)
    real(kind=realkind):: w(nx_grib,ny_grib)
    real(kind=realkind)::w_ec(nx_grib+1,ny_grib)

    real(kind=realkind):: w1(nx_grib,ny_grib)

    real(kind=realkind) wnew2(nx_grib,ny_grib)
    real(kind=realkind) w_left(nx_grib,ny_grib),w_right(nx_grib,ny_grib)
    integer wcnt(nx_grib,ny_grib),wcnt2(nx_grib,ny_grib)
    integer wcnt_left(nx_grib,ny_grib),wcnt_right(nx_grib,ny_grib)
    integer i,j,ierr,ii,jj,null,iip,jjp
    real(kind=realkind)  rmin,rmax,rmean
#ifdef MPI_SRC
#include"mpif.h"
    integer::error
#endif
    !     1. Read the field with PE 0 and distribute the error code

    ierr = 0

    if( mype==0 ) then
       call gread(lun,type,param,alev,w,nx_grib,ny_grib,ierr)

       if(lecmwfIn)then
          if(param==31 .and. type==1) then
             wcnt=1
             where(w<0._realkind .or. w>1._realkind)
                w=0._realkind
                wcnt=0
             end where
             w1=w
             w_left=cshift(w,1,1)
             w_right=cshift(w,-1,1)
             wnew2=w_left+cshift(w_left,1,2)+cshift(w_left,-1,2) &
                  +w_right+cshift(w_right,1,2)+cshift(w_right,-1,2) &
                  +cshift(w,1,2)+cshift(w,-1,2)
             wcnt_left=cshift(wcnt,1,1)
             wcnt_right=cshift(wcnt,-1,1)
             wcnt2=wcnt_left+cshift(wcnt_left,1,2)+             &
                  cshift(wcnt_left,-1,2)               &
                  +wcnt_right+cshift(wcnt_right,1,2)   &
                  +cshift(wcnt_right,-1,2)             &
                  +cshift(wcnt,1,2)+cshift(wcnt,-1,2)
             where(wcnt==0 .and. wcnt2>0) 
                w=wnew2/real(wcnt2,realkind)
             end where
             w(:,1)=w1(:,1)
             w(:,ny_grib)=w1(:,ny_grib)

             !kw 060410 read deep soil temperature (layer3) as a proxy for icefrac
             call gread(lun,112,183,28._realgribkind,w1,nx_grib,ny_grib,ierr)
             where(wcnt2==0 .and. w1<tmelt) 
                w=1._realkind
             end where
          endif

          if(param==34 .and. type==1) then
             call gread(lun,105,11,0._realgribkind,w1,nx_grib,ny_grib,ierr)
             where(w<200._realkind .or. w>400._realkind) w=w1
          endif

          if((param==80.and.type==105.and. alev< 0.5_realgribkind).or. &
               (param==85.and.type==105.and.alev<0.5_realgribkind).or. &
               param==91.and.type==105.and.alev< 0.5_realgribkind)then
             print *, ' GRRDLOC_HINT no swapns for ', &
                  param,type,alev,w(1,1)
          else
             call swapns(w,nx_grib,ny_grib)
          endif
       endif


       if(lhadley)then
          if(param==11.and.type==105.and.nx_grib==193)then
             call mod_sst(nx_grib,ny_grib,w)
             rmin=1.e10_realkind
             rmax=-1.e10_realkind
             rmean=0._realkind
             null=0
             do i = 1,nx_grib
                do j = 1,ny_grib
                   rmin=min(rmin,w(i,j))
                   rmax=max(rmax,w(i,j))
                   rmean=rmean+w(i,j)
                   if(w(i,j)<0._realkind ) null=null+1
                enddo
             enddo
             print *,'MODSST',null,rmin,rmax,rmean/real(nx_grib*ny_grib,realkind)
          endif
       endif

       if(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lmiroc5.or.lipsl.or.lgfdl)then
          if(param==235.and.type==1)then
             call mod_sst(nx_grib,ny_grib,w)
             rmin=1.e10_realkind
             rmax=-1.e10_realkind
             rmean=0._realkind
             null=0
             do i = 1,nx_grib
                do j = 1,ny_grib
                   rmin=min(rmin,w(i,j))
                   rmax=max(rmax,w(i,j))
                   rmean=rmean+w(i,j)
                   if(w(i,j)<0._realkind ) null=null+1
                enddo
             enddo
           print *,'MODSST',null,rmin,rmax,rmean/real(nx_grib*ny_grib,realkind)
          endif
       elseif(lecearth.or.lmpiesm)then
          if(param==11.and.type==102)then
             call mod_sst(nx_grib,ny_grib,w)
             rmin=1.e10_realkind
             rmax=-1.e10_realkind
             rmean=0._realkind
             null=0
             do i = 1,nx_grib
                do j = 1,ny_grib
                   rmin=min(rmin,w(i,j))
                   rmax=max(rmax,w(i,j))
                   rmean=rmean+w(i,j)
                   if(w(i,j)<0._realkind ) null=null+1
                enddo
             enddo
             print *,'MODSST',null,rmin,rmax,rmean/real(nx_grib*ny_grib,realkind)
          endif
       endif
    endif                  !  mype


#ifdef MPI_SRC
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,localComm,error)
#endif
    if(ierr/=0) then
       if( mype==0 ) print *, ' grrdloc_hint broadcast err', &
            param,type,alev,ierr
       return
    endif

#ifdef MPI_SRC
    call MPI_Bcast(w,nx_grib*ny_grib,REALTYPE,0,localComm,error)
#endif


    !     Do the interpolation 

    if( lecmwfIn ) then
       do j=1,ny_grib
          do i=1,nx_grib
             w_ec(i,j) = w(i,j)
          enddo
          w_ec(nx_grib+1,j) = w_ec(1,j)
       enddo

       do j=1,ny_local
          do i=1,nx_local
             ii = i_f(i,j)
             jj = j_f(i,j)
             iip = 1+mod(ii,nx_grib)
             jjp = min(jj+1,ny_grib)

             fFIELD(i,j)=(1.0_realkind-w_x_f(i,j))*(1.0_realkind-w_y_f(i,j))*w_ec(ii,jj) + &
                  (1.0_realkind-w_x_f(i,j))*w_y_f(i,j)*w_ec(ii,jjp)  + &
                  w_x_f(i,j)*(1.0_realkind-w_y_f(i,j))*w_ec(iip,jj)  + &
                  w_x_f(i,j)*w_y_f(i,j)*w_ec(iip,jjp)
          enddo
       enddo

       if( lgfield ) then
          do j=1,ny_local
             do i=1,nx_local
                ii = i_g(i,j)
                jj = j_g(i,j)
                iip = 1+mod(ii,nx_grib)
                jjp = min(jj+1,ny_grib)
                gFIELD(i,j)=(1._realkind-w_x_g(i,j))*(1._realkind-w_y_g(i,j))*w_ec(ii,jj)+ &
                     (1._realkind-w_x_g(i,j))*w_y_g(i,j)     *w_ec(ii,jjp)+ &
                     w_x_g(i,j)*(1._realkind-w_y_g(i,j))*w_ec(iip,jj)  + &
                     w_x_g(i,j)*w_y_g(i,j)     *w_ec(iip,jjp)
             enddo
          enddo
       endif               ! lgfield
    else
       do j=1,ny_local
          do i=1,nx_local
             ii = i_f(i,j)
             jj = j_f(i,j)
             iip = 1+mod(ii,nx_grib)
             jjp = min(jj+1,ny_grib)
! ii may be == 0 !
             ii=max(1,ii)
             fFIELD(i,j) = (1._realkind-w_x_f(i,j))*(1._realkind-w_y_f(i,j))*w(ii,jj)    + &
                  (1._realkind-w_x_f(i,j))*w_y_f(i,j)     *w(ii,jjp)  + &
                  w_x_f(i,j)*(1._realkind-w_y_f(i,j))*w(iip,jj)  + &
                  w_x_f(i,j)*w_y_f(i,j)     *w(iip,jjp)
          enddo
       enddo
       if( lgfield ) then
          do j=1,ny_local
             do i=1,nx_local
                ii = i_g(i,j)
                jj = j_g(i,j)
                iip = 1+mod(ii,nx_grib)
                jjp = min(jj+1,ny_grib)
! ii may be == 0 !
                ii=max(1,ii)
                gFIELD(i,j) = (1._realkind-w_x_g(i,j))*(1._realkind-w_y_g(i,j))*w(ii,jj)+ &
                     (1._realkind-w_x_g(i,j))*w_y_g(i,j)     *w(ii,jjp)  + &
                     w_x_g(i,j)*(1._realkind-w_y_g(i,j))*w(iip,jj)  + &
                     w_x_g(i,j)*w_y_g(i,j)     *w(iip,jjp)
             enddo
          enddo
       endif               ! lgfield
    endif                  ! lecmwfIn
    return

  end subroutine grrdloc_hint


  subroutine groploc(lun,filename,buf,nbuf,myddr)
    use decomp
    use modddr
    implicit none
    integer,intent(in)::nbuf
    integer::lun
    real(kind=realgribkind)::buf(nbuf)
    character(len=*)::filename
    type(ddr),intent(inout)::myddr
    ! Open the field data base file on PE 0
    if( mype==0 ) then
       call as2ddr(lun,filename,buf,nbuf,myddr)
    endif

    !     broadcast the ddr to all pes
    call commddr(0,myddr)

    return
  end subroutine groploc

  subroutine grrdloc_hint_lsmmatch(lun,type,param,alev,f,&
       lsmrca,lsmecham,nx_local,ny_local,nx_grib,ny_grib, &
       i_f,w_x_f,j_f,w_y_f,ierr)

    !     Special method (manual) for the Baltic Sea
    !     because we found some strange results
    !     when nearest point to Northern Baltic
    !     was found north of Scandinavia
    !     This happens as the ECHAM/OPYC land/sea mask
    !     is very far from reality
    !     
    !     The method is to only look south, east and west in the Baltic case
    !     
    !     method:
    !     
    !     the original method is to use the 4 surrounding points
    !     from the driving model (here ECHAM/OPYC)
    !     to produce one RCA-point
    !     using weights for each point's relative importance.
    !     this method is still used as long as all points have matching
    !     landseamask value.
    !     if one or two or three points have nonmatching values
    !     this/these point(s) is/are excluded from the calculation
    !     and the weights are adjusted
    !     
    !     if all 4 points have nonmatching values we have to look in a larger area.
    !     from now on we are looking for ONE point with matching landsea masks
    !     
    !     if all 4 points have nonmatching values
    !     look in the larger square, 4*4=16 points
    !     look in a special order, see below
    !     look until you find a match
    !     
    !     if all 16 points have nonmatching values
    !     look in a larger square, 6*6=36 points
    !     look in a special order, see below
    !     look until you find a match
    !     
    !     if all 36 points have nonmatching values
    !     look in a larger square, 8*8=64 points
    !     look until you find a match
    !     
    !     if all 64 points have nonmatching values
    !     look in a larger square, 10*10=100 points
    !     look in a special order, see below
    !     look until you find a match
    !     
    !     !omment:

    !     it might be surprising that we have to look in such a large area,
    !     10*10 is approximately 28 by 28 degrees
    !     and the point you pick in such a case is VERY far from
    !     the RCA-point.
    !     the reason is that ECHAM has no lakes whereas RCA does.
    !     the (few) most extreme cases are (in the grid choosen 9903)
    !     near the Caspian Sea and near Lake Onega (northeast of Lake Ladoga)
    !     i.e. lakes far from sea.
    !     
    use grw1
    use decomp
    use confys, only:tmelt
!    use gcm
    implicit none

    integer,intent(in):: lun,type,param,nx_local,ny_local,nx_grib,ny_grib

    real(kind=realgribkind),intent(in):: alev
    real(kind=realkind),intent(in)::lsmrca(nx_local,ny_local)
    real(kind=realkind),intent(out)::f(nx_local,ny_local)
    integer,intent(in):: i_f(nx_local,ny_local), j_f(nx_local,ny_local)
    real(kind=realkind),intent(in)::  w_x_f(nx_local,ny_local), w_y_f(nx_local,ny_local)

    integer nomismatch
    real(kind=realkind) w1,w2,w3,w4,sumw,ourlsm
    logical lmismatch1,lmismatch2,lmismatch3,lmismatch4



    real(kind=realkind):: w(nx_grib,ny_grib)
    real(kind=realkind):: lsmecham(nx_grib,ny_grib)
    logical lbothnia(nx_grib,ny_grib)


    real(kind=realkind)  rmin,rmax,rmean
    integer i,j,ierr,ii,jj,iut
    integer ii1,ii2,ii3,ii4,jj1,jj2,jj3,jj4
    integer ii0,ii5,jj0,jj5
    integer ilon,jlat
    logical ltest,lfirst
    logical ltest2
    data ltest,lfirst,ltest2/.false.,.true.,.false./
    save ltest,lfirst,ltest2
#ifdef MPI_SRC
#include"mpif.h"
    integer::error
#endif
    !     1. read the field with pe 0 and distribute the error code


    if(.not.lecham5)then
       if (mype == 0 .and. lfirst) then
          write(6,'(a)') ' '
          write(6,'(a)') 'warning warning warning warning warning'
          write(6,'(a)') ' '
          write(6,'(a)') 'warning from grrdloc_hint_lsmmatch'
          write(6,'(a)') 'reading ice cover and sst from echam/opyc'
          write(6,'(a)') ' '
          write(6,'(a)') 'special method (manual) for the baltic sea'
          write(6,'(a)') 'because we found some strange results'
          write(6,'(a)') 'in the northern baltic'
          write(6,'(a)') 'when nearest echam sea-point'
          write(6,'(a)') 'was found north of scandinavia'
          write(6,'(a)') 'and not in the baltic.'
          write(6,'(a)') 'this happensastheecham/opycland/sea mask'
          write(6,'(a)') 'is very far from reality.'
          write(6,'(a)') 'the method is to only look south, east and '
          write(6,'(a)') 'west in the bothnia sea case.'
          write(6,'(a,a)') 'i.e. do not look north when you look ', &
               'for a sea-point'
          write(6,'(a)') ' '
          write(6,'(a)') 'the method is tailored for a specific grid'
          write(6,'(a)') 'in this case the 106*102 european grid'
          write(6,'(a)') 'used for the first transient runs in 2004.'
          write(6,'(a)') 'if you change the grid you will probably'
          write(6,'(a)') 'find some funny ice cover and sst'
          write(6,'(a)') 'in some part of the world.'
          write(6,'(a)') ' '
          write(6,'(a)') 'warning warning warning warning warning'
          write(6,'(a)') ' '
       endif
    endif

    ierr = 0
    if (mype==0) then
       call gread(lun,type,param,alev,w,nx_grib,ny_grib,ierr)

       rmin=1.e10_realkind 
       rmax=-1.e10_realkind 
       rmean=0._realkind 
       do i = 1,nx_grib
          do j = 1,ny_grib
             rmin=min(rmin,w(i,j))
             rmax=max(rmax,w(i,j))
             rmean=rmean+w(i,j)
          enddo
       enddo
       if (param == 11) then
          print *, ' tsea echam ',rmin,rmax,rmean/real(nx_grib*ny_grib,realkind)
       elseif (param == 91) then
          print *, ' frice echam ',rmin,rmax,rmean/real(nx_grib*ny_grib,realkind)
       endif
    endif                     !  mype


#ifdef MPI_SRC
    call MPI_Bcast(ierr,1,MPI_INTEGER,0,localComm,error)
#endif
    if(ierr/=0) return

    !   2. send to all processors 
#ifdef MPI_SRC
    call MPI_Bcast(w,nx_grib*ny_grib,REALTYPE,0,localComm,error)
    call MPI_Bcast(lsmecham,nx_grib*ny_grib,REALTYPE,0,localComm,error)
#endif     

    if (ltest) then
       rmin=1.e10_realkind 
       rmax=-1.e10_realkind 
       rmean=0._realkind 
       do i = 1,nx_local
          do j = 1,ny_local
             rmin=min(rmin,lsmrca(i,j))
             rmax=max(rmax,lsmrca(i,j))
             rmean=rmean+lsmrca(i,j)
          enddo
       enddo
       print *, ' lsm rca ',mype,rmin,rmax,rmean/real(nx_local*ny_local,realkind)
    endif


    !    define Bothnia = northern part Baltic Sea


    lbothnia=.false.

    if(.not.lecham5)then
       do i=40,42
          lbothnia(i,17)=.true.
       enddo
       do i=41,42
          lbothnia(i,18)=.true.
       enddo


       if(lfirst .and. mype==0) then
          print *, 'Bothnic area '
          do j=18,14,-1
             print *,(lbothnia(i,j),i=36,45)
          enddo
          rmin=0._realkind
          rmax=0._realkind
          do i = 1,nx_grib
             do j = 1,ny_grib
                if (lbothnia(i,j)) then
                   rmin=rmin+1._realkind
                else
                   rmax=rmax+1._realkind
                endif
             enddo
          enddo
          print *, 'Bothnia points ',rmin
          print *, 'non-Bothnia points ',rmax
       endif                  ! mype
    endif


    if (ltest .and. mype == 10) then
       print *, 'lsm RCA Baltic area '
       do j=27,7,-1
          write(6,'(17f5.2)')(lsmrca(i,j),i=12,28)
       enddo

    endif

    iut=6

    do j=1,ny_local
       do i=1,nx_local
          ii = i_f(i,j)
          jj = j_f(i,j)
          if (lbothnia(ii,jj)) goto 2

          nomismatch=0
          lmismatch1=.false.
          lmismatch2=.false.
          lmismatch3=.false.
          lmismatch4=.false.
          !     uh the boundary file area might be too small
          if ( ii < 1 .or. ii > nx_grib-1 .or. &
               jj < 1 .or. jj > ny_grib-1 ) then
             print *, ' autmp ',i,j,ii,jj,mype
          endif

          !    lsmecham have the values 0 and 1
          !     lsmrca vary between 0 and 1

          if (abs(lsmecham(ii,jj))>  1.e14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch1=.true.
          endif
          if (abs(lsmecham(ii,jj+1)) > 1.e-14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch2=.true.
          endif
          if ((abs(lsmecham(ii+1,jj)))>  1.e14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch3=.true.
          endif
          if(abs(lsmecham(ii+1,jj+1))>  1.e14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch4=.true.
          endif



          if (nomismatch == 0) then
             f(i,j) = (1._realkind -w_x_f(i,j))*(1._realkind -w_y_f(i,j))*w(ii,jj)    +  &
                  (1._realkind -w_x_f(i,j))*w_y_f(i,j)     *w(ii,jj+1)  + &
                  w_x_f(i,j)*(1._realkind -w_y_f(i,j))*w(ii+1,jj)  + &
                  w_x_f(i,j)*w_y_f(i,j)     *w(ii+1,jj+1)
          elseif(nomismatch>=1 .and. nomismatch<=3) then ! 1-3 mismatch
             sumw=0._realkind 
             w1 = (1._realkind -w_x_f(i,j))*(1._realkind -w_y_f(i,j))
             w2 = (1._realkind -w_x_f(i,j))*w_y_f(i,j)
             w3 =      w_x_f(i,j)*(1._realkind -w_y_f(i,j))
             w4 =      w_x_f(i,j)*w_y_f(i,j)
             if (.not. lmismatch1) sumw=sumw+w1
             if (.not. lmismatch2) sumw=sumw+w2
             if (.not. lmismatch3) sumw=sumw+w3
             if (.not. lmismatch4) sumw=sumw+w4
             f(i,j) = 0._realkind 
             if(.not.lmismatch1) f(i,j)=f(i,j) + w1/sumw*w(ii,jj)
             if(.not.lmismatch2) f(i,j)=f(i,j) + w2/sumw*w(ii,jj+1)
             if(.not.lmismatch3) f(i,j)=f(i,j) + w3/sumw*w(ii+1,jj)
             if(.not.lmismatch4) f(i,j)=f(i,j) + w4/sumw*w(ii+1,jj+1)
             !     if(ltest)then
             !     if(nomismatch==1) print *,mype,' use 3 ',i,j,ii,jj
             !     if(nomismatch==2) print *,mype,' use 2 ',i,j,ii,jj
             !     if(nomismatch==3) print *,mype,' use 1 ',i,j,ii,jj
             !     endif
          elseif (nomismatch == 4) then
             if(ltest)then
                print *,mype,'4 mismatchs, have to do something'
                call stop_program('')
             endif
             !      
             !     so far we searched in a square of 4 points
             !     now look in the next larger square, 16 points (- the 4 non-usable points)
             !     
             !     first idea:
             !===========
             !     first find real latitudes & longitudes
             !     De-rotate coordinates
             !     or get local lat lon
             !     then calculate distance
             !     then find smallest distance
             !     and pick the point, if masks agree
             !     
             !     simpler
             !==========
             !     find a point in outer square where masks agree
             !     first look east, north, west and south   of inner square
             !     then look NE, NW, SW and SE of inner square
             !     
             ii1=ii-1
             ii2=ii
             ii3=ii+1
             ii4=ii+2
             jj1=jj-1
             jj2=jj
             jj3=jj+1
             jj4=jj+2

             if ( ii1 < 1 .or. ii4 > nx_grib .or. &
                  jj1 < 1 .or. jj4 > ny_grib ) then
                print *, 'grrdlocnewecham ',i,j,ii1,ii4,jj1,jj4,mype
                print *, ' larger square (4*4) outside border'
                print *, ' have to improve programme, should abort?'
                print *, ' should abort?'
                print *, ' mype param ',mype,param
                goto 3
                call stop_program('')
             endif

             ourlsm=lsmrca(i,j)

             !     east
             if (abs(lsmecham(ii4,jj2)- 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj2)
             elseif(abs(lsmecham(ii4,jj3) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj3)

                !     north
             elseif (abs(lsmecham(ii3,jj4) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii3,jj4)
             elseif (abs(lsmecham(ii2,jj4) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii2,jj4)

                !     west
             elseif (abs(lsmecham(ii1,jj3) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj3)
             elseif (abs(lsmecham(ii1,jj2) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj2)

                !     south
             elseif (abs(lsmecham(ii2,jj1) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii2,jj1)
             elseif (abs(lsmecham(ii3,jj1) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii3,jj1)

                !     NE
             elseif (abs(lsmecham(ii4,jj4) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj4)
                !     NW
             elseif (abs(lsmecham(ii1,jj4) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj4)
                !     SW
             elseif (abs(lsmecham(ii1,jj1) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj1)
                !     SE
             elseif (abs(lsmecham(ii4,jj1) - 1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj1)

             else             !4*4=16 not enough, check 6*6=36
                !     
                !     so far we searched in a square of 16 points
                !     now look in the next larger square, 36 points (- the 16 non-usable points)

                !     find a point in outer square where masks agree
                !     look from SE corner , then anti-clockwise
                !     
                ii0=ii-2
                ii5=ii+3
                jj0=jj-2
                jj5=jj+3
                if( ltest)then 
                   print *, ' 16 mismatchs, have to do something '
                   print *,'mism=16 mype i j ',mype,i,j
                   print *,'ii jj ',ii,jj
                endif

                if ( ii0 < 1 .or. ii5 > nx_grib .or. &
                     jj0 < 1 .or. jj5 > ny_grib ) then
                   print *,  i,j,ii0,ii5,jj0,jj5,mype
                   print *,' larger square (6*6) outside border'
                   print *,'have to improve program, should abort?'
                   print *,' mype param ',mype,param
                   goto 3

                endif


                !     first eastern side

                if (abs(lsmecham(ii5,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj0)
                   if (ltest)print *, ' match in ',17,' try',mype
                elseif (abs(lsmecham(ii5,jj1) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj1)
                   if (ltest)print *, ' match in ',18,' try',mype
                elseif (abs(lsmecham(ii5,jj2) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj2)
                   if (ltest)print *, ' match in ',19,' try',mype
                elseif (abs(lsmecham(ii5,jj3) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj3)
                   if(ltest)print *, ' match in ',20,' try',mype
                elseif (abs(lsmecham(ii5,jj4) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj4)
                   if(ltest)print *, ' match in ',21,' try',mype
                elseif (abs(lsmecham(ii5,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj5)
                   if(ltest)print *, ' match in ',22,' try',mype
                   !     northern side

                elseif (abs(lsmecham(ii4,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii4,jj5)
                   if(ltest)print *, ' match in ',23,' try',mype
                elseif (abs(lsmecham(ii3,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii3,jj5)
                   if(ltest)print *, ' match in ',24,' try',mype
                elseif (abs(lsmecham(ii2,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii2,jj5)
                   if(ltest)print *, ' match in ',25,' try',mype
                elseif (abs(lsmecham(ii1,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii1,jj5)
                   if(ltest)print *, ' match in ',26,' try',mype
                elseif (abs(lsmecham(ii0,jj5) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj5)
                   if(ltest)print *, ' match in ',27,' try',mype

                   !     western side

                elseif (abs(lsmecham(ii0,jj4) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj4)
                   if(ltest)print *, ' match in ',28,' try',mype
                elseif (abs(lsmecham(ii0,jj3) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj3)
                   if(ltest)print *, ' match in ',29,' try',mype
                elseif (abs(lsmecham(ii0,jj2) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj2)
                   if(ltest)print *, ' match in ',30,' try',mype
                elseif (abs(lsmecham(ii0,jj1) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj1)
                   if(ltest)print *, ' match in ',31,' try',mype
                elseif (abs(lsmecham(ii0,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj0)
                   if(ltest)print *, ' match in ',32,' try',mype

                   !     southern side

                elseif (abs(lsmecham(ii1,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii1,jj0)
                   if(ltest)print *, ' match in ',33,' try',mype
                elseif (abs(lsmecham(ii2,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii2,jj0)
                   if(ltest)print *, ' match in ',34,' try',mype
                elseif (abs(lsmecham(ii3,jj0))>1.e-14_realkind ) then
                   f(i,j)=w(ii3,jj0)
                   if(ltest)print *, ' match in ',35,' try',mype
                elseif (abs(lsmecham(ii4,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii4,jj0)
                   if(ltest)print *, ' match in ',36,' try',mype


                else          !6*6=36 not enough, check 8*8=64

                   !     problem no fit among 36 nearest points
                   !     not very likely
                   !     but it did happen
                   if(ltest)then 
                      print *,'36 mismatchs, have to do something'
                      print *,'mism=36 mype i j ',mype,i,j
                      print *,'ii jj ',ii,jj
                   endif

                   !     go through whole square although I  know that inner points
                   !     won't give a fit
                   !     - simpler programming
                   if ( ii-3 < 1 .or. ii+4 > nx_grib .or.   &
                        jj-3 < 1 .or. jj+4 > ny_grib ) then
                      print *,i,j,ii1,ii4,jj1,jj4,   mype
                      print *, 'larger square (8*8) outside border'
                      print *,'have to improve program,abort?'
                      !     print *, 'have to abort'
                      print *, ' mype param ',mype,param
                      goto 3
                      !     call stop_program('')
                   endif

                   !     south to north
                   jlat=jj-3
                   f(i,j)=-999._realkind

                   do while(f(i,j) < 0._realkind  .and. jlat <= jj+4)
                      !     west to east
                      ilon=ii-3

                      do while(f(i,j) < 0._realkind  .and. ilon <= ii+4) 
                         if (abs(lsmecham(ilon,jlat) -1._realkind)>1.e-14_realkind ) then
                            f(i,j)=w(ilon,jlat)
                            if(ltest)then
                               print *, 'match in 37-64 try',mype
                               print *, 'lon lat ',ilon,jlat
                            endif
                         endif
                         ilon=ilon+1
                      enddo   ! inner while

                      jlat=jlat+1
                   enddo      ! outer while

                   if (f(i,j) < 0._realkind )then
                      !     problem no fit among 64 nearest points
                      !     not very likely
                      !     but it did happen
                      if(ltest)then 
                         print *,'64 mismatchs, do something'
                         print *,'mism=64 mype i j ',mype,i,j
                         print *,'ii jj ',ii,jj
                      endif

                      if ( ii-4 < 1 .or. ii+5 > nx_grib .or. &
                           jj-4 < 1 .or. jj+5 > ny_grib ) then
                         print *, i,j,ii1,ii4,jj1,jj4, mype
                         print *,'larger square outside border'
                         print *,'have to improve program, abort?'
                         print *, ' mype param ',mype,param
                         goto 3
                      endif

                      !     south to north
                      jlat=jj-4

                      do while(f(i,j) < 0._realkind  .and. jlat <= jj+5)
                         !     west to east
                         ilon=ii-4

                         do while(f(i,j) < 0._realkind  .and. ilon <= ii+5) 
                            if (abs(lsmecham(ilon,jlat) -1._realkind)>1.e-14_realkind ) then
                               f(i,j)=w(ilon,jlat)
                               if(ltest)then
                                  print *,'match in 85-100 try',mype
                                  print *, ' lon lat ',ilon,jlat
                               endif
                            endif
                            ilon=ilon+1
                         enddo ! inner while

                         jlat=jlat+1
                      enddo   ! outer while

                      if (f(i,j) < 0._realkind )then
                         !     problem no fit among 100 nearest points
                         !     not very likely
                         print *,' 100 mismatchs, abort? '
                         print *, 'Bothnia mype param ',mype,param
                         write(6,'(a,5i3,a)')                &
                              '100 mism i j ii jj mype ',i,j,ilon,jlat,mype
                         goto 3

                      endif   !inner f<0
                   endif      !outer f<0
                endif         ! 36,64,100
             endif            ! 16
          endif               ! 4

          !     cycle INNER
          goto 3

2         continue

          nomismatch=0
          lmismatch1=.false.
          lmismatch2=.false.
          lmismatch3=.false.
          lmismatch4=.false.
          if ( ii < 1 .or. ii > nx_grib-1 .or.           & 
               jj < 1 .or. jj > ny_grib-1 ) then
             print *, ' autmp2 ',i,j,ii,jj,mype

          endif
          !     Which one are we missing and count the number of misses
          if (abs(lsmecham(ii,jj))>1.e-14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch1=.true.
          endif
          if (abs(lsmecham(ii,jj+1))>1.e-14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch2=.true.
          endif
          if (abs(lsmecham(ii+1,jj))>1.e-14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch3=.true.
          endif
          if(abs(lsmecham(ii+1,jj+1))>1.e-14_realkind ) then
             nomismatch = nomismatch+1
             lmismatch4=.true.
          endif

          if(ltest2)then
             print *,'i j ii jj mype ',i,j,ii,jj,mype
             print *,'LSM'
             print *,lsmecham(ii,jj+1),lsmecham(ii+1,jj+1)
             print *,lsmecham(ii,jj),lsmecham(ii+1,jj)

             print *,'Values'
             print *,w(ii,jj+1),w(ii+1,jj+1)
             print *,w(ii,jj),w(ii+1,jj)

             print *,'mismatchs ',nomismatch

          endif

          if(nomismatch>=0 .and. nomismatch<=3) then ! 0-3 mismatch
             !compute weights
             w1 = (1._realkind -w_x_f(i,j))*(1._realkind -w_y_f(i,j))
             w2 = (1._realkind -w_x_f(i,j))*w_y_f(i,j)
             w3 =  w_x_f(i,j)*(1._realkind -w_y_f(i,j))
             w4 =  w_x_f(i,j)*w_y_f(i,j)

             if(ltest2)then
                print *,'mismatchs ',nomismatch
                print *,lmismatch2,lmismatch4
                print *,lmismatch1,lmismatch3

             endif

             !     test with sumw
             sumw=0._realkind 
             if (.not. lmismatch1) sumw=sumw+w1
             if (.not. lmismatch2 .and. lbothnia(ii,jj+1))sumw=sumw+w2
             if (.not. lmismatch3) sumw=sumw+w3
             if (.not. lmismatch4 .and. lbothnia(ii+1,jj+1))sumw=sumw+w4

             if(ltest2)then
                print *,'weights'
                print *,w2,w4
                print *,w1,w3

                print *,'weights/sum ',sumw
                print *,w2/sumw,w4/sumw
                print *,w1/sumw,w3/sumw

                print *,'Bothnia'
                print *,lbothnia(ii,jj+1),lbothnia(ii+1,jj+1)
                print *,lbothnia(ii,jj),lbothnia(ii+1,jj)

             endif

             f(i,j) = 0._realkind 
             if (.not. lmismatch1)then
                f(i,j) = f(i,j) + w1/sumw*w(ii,jj)
             endif
             if (.not. lmismatch2 .and. lbothnia(ii,jj+1)) &
                  f(i,j) = f(i,j) + w2/sumw*w(ii,jj+1)
             if (.not. lmismatch3) f(i,j)=f(i,j) + w3/sumw*w(ii+1,jj)
             if (.not. lmismatch4 .and. lbothnia(ii+1,jj+1))  &
                  f(i,j) = f(i,j) + w4/sumw*w(ii+1,jj+1)
             if(ltest2 .and. abs(sumw)>1.e-14_realkind )then
                write(6,'(a,f9.2,4i3,i6)')'result i j ii jj mype ', &
                     f(i,j),i,j,ii,jj,mype
             endif

             if(ltest2)then
                if(nomismatch==0) print *,mype, &
                     'Bothnia  use 4 ',i,j,ii,jj
                if(nomismatch==1) print *,mype, &
                     'Bothnia  use 3 ',i,j,ii,jj
                if(nomismatch==2) print *,mype, &
                     'Bothnia  use 2 ',i,j,ii,jj
                if(nomismatch==3) print *,mype, &
                     'Bothnia  use 1 ',i,j,ii,jj
                if(.not.lbothnia(ii,jj+1) .and. .not. lmismatch2)  &
                     print *,mype,'Bothnia  MINUS 1 '
                if(.not.lbothnia(ii+1,jj+1) .and. .not. lmismatch4)  &
                     print *,mype,'Bothnia  MINUS 1 '
             endif            !end ltest2

          elseif(nomismatch==4)then
             if(ltest2)then
                print *, mype,'iBhnia  4 mismatchs, have to do something'
             endif

             ii1=ii-1
             ii2=ii
             ii3=ii+1
             ii4=ii+2
             jj1=jj-1
             jj2=jj
             jj3=jj+1
             jj4=jj+2

             if ( ii1 < 1 .or. ii4 > nx_grib .or. &
                  jj1 < 1 .or. jj4 > ny_grib ) then
                print *, 'grrdlocnewecham ',i,j,ii1,ii4,jj1,jj4,mype
                print *, ' larger square (4*4) outside border'
                print *, ' have to improve programme, should abort?'
             endif

             ourlsm=lsmrca(i,j)

             !     east
             if (abs(lsmecham(ii4,jj2) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj2)
                if (ltest2)then
                   print *,'Bothnia  match in ',5,' try',mype,f(i,j)
                endif

                !     west
             elseif (abs(lsmecham(ii1,jj2) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj2)
                if (ltest2)then
                   print *,'Bothnia  match in ',6,' try',mype,f(i,j)
                endif

                !     south
             elseif (abs(lsmecham(ii2,jj1) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii2,jj1)
                if (ltest2)then
                   print *,'Bothnia  match in ',7,' try',mype,f(i,j)
                endif
             elseif (abs(lsmecham(ii3,jj1) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii3,jj1)
                if (ltest2)then
                   print *,'Bothnia  match in ',8,' try',mype,f(i,j)
                endif

                !     SW
             elseif (abs(lsmecham(ii1,jj1) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii1,jj1)
                if (ltest2)then
                   print *,'Bothnia  match in ',9,' try',mype,f(i,j)
                endif
                !     SE
             elseif (abs(lsmecham(ii4,jj1) -1._realkind)>1.e-14_realkind ) then
                f(i,j)=w(ii4,jj1)
                if (ltest2)then
                   print *,'Bothnia  match in ',10,'try',mype,f(i,j)
                endif

             else             !4*4=16 not enough, check 6*6=36

                ii0=ii-2
                ii5=ii+3
                jj0=jj-2
                jj5=jj+3
                if(ltest2)then 
                   print *,'Bothnia  16 mismatchs,do something'
                   print *,'mism=16 mype i j ',mype,i,j
                   print *,'ii jj ',ii,jj
                endif

                if ( ii0 < 1 .or. ii5 > nx_grib .or.                 &
                     jj0 < 1 .or. jj5 > ny_grib ) then
                   print *,   i,j,ii1,ii4,jj1,jj4,mype
                   print *,'larger square (6*6) outside border'
                   print *,'have to improve program, abort?'

                   print *, 'Bothnia mype param ',mype,param
                endif


                !     first eastern side

                if (abs(lsmecham(ii5,jj0) -1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj0)
                   if (ltest2)then
                      print *,'Bothni match in',11,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii5,jj1)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj1)
                   if (ltest2)then
                      print *,'Bothni match in',12,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii5,jj2)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii5,jj2)
                   if (ltest2)then
                      print *,'Bothni match in',13,'try',mype,f(i,j)
                   endif

                   !     western side

                elseif (abs(lsmecham(ii0,jj2)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj2)
                   if(ltest2)then
                      print *,'Bothni match in',14,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii0,jj1)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj1)
                   if(ltest2)then
                      print *,'Bothni match in',15,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii0,jj0)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii0,jj0)
                   if(ltest2)then
                      print *,'Bothni match in',16,'try',mype,f(i,j)
                   endif

                   !     southern side

                elseif (abs(lsmecham(ii1,jj0)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii1,jj0)
                   if(ltest2)then
                      print *,'Bothni match in',17,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii2,jj0)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii2,jj0)
                   if(ltest2)then
                      print *,'Bothni match in',18,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii3,jj0)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii3,jj0)
                   if(ltest2)then
                      print *,'Bothni match in',19,'try',mype,f(i,j)
                   endif
                elseif (abs(lsmecham(ii4,jj0)-1._realkind)>1.e-14_realkind ) then
                   f(i,j)=w(ii4,jj0)
                   if(ltest2)then
                      print *,'Bothni match in',20,'try',mype,f(i,j)
                   endif


                else          !6*6=36 not enough, check 8*8=64

                   if(ltest2)then 
                      print *,'Bothnia 36 mismatchs, do something'
                      print *,'mism=36 mype i j ',mype,i,j
                      print *,'ii jj ',ii,jj
                   endif

                   if ( ii-3 < 1 .or. ii+4 > nx_grib .or. &
                        jj-3 < 1 .or. jj+4 > ny_grib ) then
                      print *,'grrdlocnewecham',i,j,ii1,ii4,jj1,jj4,mype
                      print *, 'larger square (8*8) outside border'
                      print *, 'have to improve programme,abort?'

                      print *, 'Bothnia mype param ',mype,param

                   endif

                   !     south to north
                   jlat=jj-3
                   f(i,j)=-999._realkind 

                   do while(f(i,j) < 0._realkind  .and. jlat <= jj)
                      !     west to east
                      ilon=ii-3

                      do while(f(i,j) < 0._realkind  .and. ilon <= ii+4) 
                         if (abs(lsmecham(ilon,jlat)-1._realkind)>1.e-14_realkind ) then
                            f(i,j)=w(ilon,jlat)
                            if(ltest2)then
                               print *,'Bothnia match in 37-64 try',mype,f(i,j)
                               print *, 'lon lat ',ilon,jlat
                            endif
                         endif
                         ilon=ilon+1
                      enddo   ! inner while

                      jlat=jlat+1
                   enddo      ! outer while

                   if (f(i,j) < 0._realkind )then
                      if(ltest2)then 
                         print *,'Bothnia 64 mismatchs, do something'
                         print *,'mism=64 mype i j ',mype,i,j
                         print *,'ii jj ',ii,jj
                      endif

                      if ( ii-4 < 1 .or. ii+5 > nx_grib .or. &
                           jj-4 < 1 .or. jj+5 > ny_grib ) then
                         print *,'grrdlocnewecham',i,j,ii1,ii4,jj1,jj4,mype
                         print *,'largr square (10*10)outsid border'
                         print *, 'have to improve programme,abort?'
                         print *, 'Bothnia mype param ',mype,param
                      endif

                      !     south to north
                      jlat=jj-4

                      do while(f(i,j) < 0._realkind  .and. jlat <= jj)
                         !     west to east
                         ilon=ii-4

                         do while(f(i,j) < 0._realkind  .and. ilon <= ii+5) 
                            if (abs(lsmecham(ilon,jlat)-1._realkind)>1.e-14_realkind ) then
                               f(i,j)=w(ilon,jlat)
                               if(ltest2)then
                                  print *,'Bothnia match in 85-100 try',mype,f(i,j)
                                  print *, ' lon lat ',ilon,jlat
                               endif
                            endif
                            ilon=ilon+1
                         enddo ! inner while

                         jlat=jlat+1
                      enddo   ! outer while

                      if (f(i,j) < 0._realkind )then
                         print *,' 100 Bothnia mismatchs,abort?'
                         print *, 'Bothnia mype param ',mype,param
                         write(6,'(a,5i3,a)')'100 mism i j ii jj mype' &
                              ,i,j,ilon,jlat,mype
                         goto 3
                         call stop_program('')
                      endif   !inner f<0
                   endif      !outer f<0
                endif         ! 36,64,100
             endif            ! 16
          endif               ! 4

3         continue
          if(param==11 .and. abs(f(i,j))<1.e-14_realkind ) write(6,'(a,5i3,l2)') &
               'sst=0K i j ii jj mype ',i,j,ii,jj,mype,lbothnia(ii,jj)
          !     end do INNER !i
       enddo                 !i
    enddo                    !j

    if (ltest .and. mype == 10) then
       if(param==11) then
          print *, 'sst Baltic area'
          write(66,*) 'sst Baltic area'
          do j=27,7,-1
             write(6,'(17f5.0)')(f(i,j)-tmelt,i=12,28)
             write(66,'(17f5.0)')(f(i,j)-tmelt,i=12,28)
          enddo
       else
          print *, 'ice Baltic area'
          write(66,*) 'ice Baltic area'
          do j=27,7,-1
             write(6,'(17f5.2)')(f(i,j),i=12,28)
             write(66,'(17f5.2)')(f(i,j),i=12,28)
          enddo
       endif
    endif

    if (ltest .and. mype == 14) then
       if(param==11) then
          print *, 'sst Baltic area 14'
          write(67,*) 'sst Baltic area 14'
          do j=21,1,-1
             write(6,'(17f5.0)')(f(i,j)-tmelt,i=12,28)
             write(67,'(17f5.0)')(f(i,j)-tmelt,i=12,28)
          enddo
       else
          print *, 'ice Baltic area 14'
          write(67,*) 'ice Baltic area 14'
          do j=21,1,-1
             write(6,'(17f5.2)')(f(i,j),i=12,28)
             write(67,'(17f5.2)')(f(i,j),i=12,28)
          enddo
       endif


    endif

    lfirst=.false.

    return
  end subroutine grrdloc_hint_lsmmatch
subroutine grclos(ludir)
    use asimof
    implicit none
    integer,intent(in):: ludir
    integer:: kfld(1),klenf,kb1(30),kmode,kerr

    klenf=-1
    call asimhw(ludir,kfld,klenf,kb1,kmode,kerr)
    close(ludir)
  end subroutine grclos

  subroutine swapns(f,klon,klat)
    implicit none
    integer,intent(in):: klon,klat
    real(kind=realkind),intent(inout):: f(klon,klat)

    integer i,j
    real(kind=realkind) zf(klat)

    do i=1,klon
       do j=1,klat
          zf(j) = f(i,klat-j+1)
       enddo
       do j=1,klat
          f(i,j) = zf(j)
       enddo
    enddo

    return
  end subroutine swapns

  subroutine mod_sst (nx, ny, a)  
    implicit none  
    integer :: nx, ny  
    real(kind=realkind) :: a(nx, ny)

    real(kind=realkind):: b(nx,ny)
  
    integer :: i, j, ii, jj, k, mymodl, mymodr, spec

    mymodl = 0  
    mymodr = 0  
    spec = 0  
    !***  search for sea-point in longitud direction only
    do j = 1, ny  
       jj = j  
       do i = 1, nx  
          !     miss val = -1.
          if (a (i, j) <0._realkind ) then  
!             do k = 1, 50
             do k = 1, 100  
                ii = i-k
!                if (ii<1) ii=nx+ii
                if (ii>=1 ) then  
                if (a(ii,jj) >0._realkind ) then  
                   b (i, j) = a (ii, jj)  
                   mymodl = mymodl+1  
                   goto 100
                endif
                endif
                ii = i+k  
!                if (ii>nx) ii = ii-nx
                if (ii<=nx) then
                if (a(ii,jj) >0._realkind ) then
                   b (i,j) = a(ii, jj)  
                   mymodr = mymodr+1  
                   goto 100
                endif
                endif
             enddo  ! k
             !     au  special
             b(i,j) = 273.16_realkind   
             spec = spec + 1  
             if(.not.(lcanesm2.or.lcnrm.or.lnoresm.or.lhadgem.or.lecearth.or.lmiroc5.or.lipsl.or.lmpiesm.or.lgfdl))print *, ' mod sst spec ',i,j
          else
             b(i,j)=a(i,j)
          endif
100       continue  
       enddo  ! i
    enddo     ! j
    write (6,  * ) ' mod_sst in x=  ', mymodl, mymodr, spec  

    do j = 1, ny  
       do i = 1, nx  
          a(i,j)=b(i,j)
       enddo
    enddo
    return  
  end subroutine mod_sst

  subroutine as2ddr(kdev,filename,buf,lenrec,myddr)
    use asimof
    use modddr
    implicit none
    type(ddr),intent(inout)::myddr
    integer,intent(in)::kdev
    character(len=*),intent(in):: filename
    integer,intent(in)::lenrec
    real(kind=realgribkind),intent(in)::buf(lenrec)
    integer,parameter::klb1=30
    integer,parameter::knfldx=8000
    integer,parameter::klevx=100
    integer,parameter::kppb2=400
    integer,parameter::kpkb2=kppb2
    integer:: kbs1(klb1,knfldx),kb2(kpkb2)
    real(kind=realgribkind):: pab(2,klevx)
    real(kind=realgribkind)::pb2(kppb2)
    logical::logab
    integer:: kdtg,knfld,klev,knmf,knsf,ksptr
    integer:: ierr,j,ipf
    integer ivdtg,ivsec,ibdtg,ibsec,iflen,islen
    real(kind=realgribkind):: zsup(1)
    integer,parameter::jplb1=40
    integer ib1(jplb1)
    logical::debug=.false.
    integer::retlen

    !     purpose: open an asimof file and construct comddr from it

    !     interface:
    !     kdev   i  unit number of asimof file            |
    !     kbs1   o  array of keys in the file             |
    !     klb1   i  first dimension of kbs1               |
    !     knfldx i  second dimension of kbs1              | as used by
    !     pab    o  list of eta parameters of the levels  |
    !     klevx  i  second dimension of pab (first is 2)  |
    !     kb2    o  block 2 of the orography              |
    !     kpkb2  i  dimension of kb2                      |
    !     pb2    o  block 2 of the orography              |
    !     kppb2  i  dimension of pb2                      | as used by getfd
    !     buf      o  orography                             |
    !     lenrec     io on input:  dimension of p             |
    !     on output: number of field values     |
    !     
    !     method:
    !     loadfd is used to construct a list of fields and levels.
    !     getfd is used to extract orography, which field is used to
    !     define the area.  also, orography is returned in p, because
    !     it has been decoded already by this routine.
    !     the other fields are found in the list kbs1 at positions
    !     kmprt and ksprt resp.  e.g. orography is found with keys
    !     kbs1(*,ksptr) - orography is always the first in the list
    !     of fields. kmptr(1) points to the first multi-level field on
    !     the first level, kmptr(2) points to the first multi-level
    !     field on the second level, kmptr(klev+1) is the second
    !     multi-field on the first level etc.
    !     all fields valid at kdtg are selected.  if kdtg=0, all
    !     fields valid at that of the first orography field are used.
    !     
    !     the ddr is constructed, assuming it is the ddr of a grib file
    !     (ntyphl=10, pointers are to 'record' numbers - but the
    !     record lengths lrechl are not set.

    !     asimof parameters

    data logab /.true./

    call reset_ddr(myddr)

    !     load asimof input file
    call loadfd(kdev,filename,kbs1,klb1,knfldx,knfld,logab,pab,klevx,klev,0,ierr)

    if(ierr/=0)then
       if(ierr==20)then
          print *,'error opening ',filename
       elseif(ierr==30)then
          print *,'error reading ',filename
       elseif(ierr==40)then
          print *,'file ',filename,'is empty'
       elseif(ierr==31)then
          print *,'more than knfldx fields in file. return only the first knfldx',knfldx
       elseif(ierr==51)then
          print *,'klbk1 < 30',filename
       elseif(ierr==70)then
          print *,'decode error',filename
       elseif(ierr==101)then
          print *,'more than knlevx levels in file. return only the first klevx',klevx,filename
       elseif(ierr==102)then
          print *,'not enough vertical level parameters to'// &
               'define the level of some field in the file', &
               filename
       elseif(ierr>200)then
          print *,'error reported by getfd or putfd',filename
       endif
       print *,'error code ',ierr,'returned by loadfd as2ddr'
       print *,filename
       call stop_program('')
    endif

    !     create list of fields and extract area from orography

    !     list of fields

    call setlof(kbs1,klb1,knfld,ksptr,kdtg,knmf,klev,knsf,myddr%nwmmhl,myddr%nwmshl,&
         myddr%nslthl,myddr%alevhl)

    !     if orography was not found for this time, take any one
    if(ksptr<=0)then
       print *,'orography not found - abort'
       call stop_program('')
    endif

    !     define area
    ipf=1
    do 231 j=1,min(klb1,jplb1)
       ib1(j)=kbs1(j,ksptr)
231 enddo
    call getfd(kdev,ib1,min(jplb1,klb1),kb2,kpkb2,pb2,kppb2,zsup,0,&
         buf,lenrec,ipf,ierr,retlen)
    if(ierr/=0)then
       print *,'error code ',ierr,'returned by getfd as2ddr'
       print *,filename
       call stop_program('')
    endif

    !     define dates and times
    call getbvt(kbs1,ivdtg,ivsec,ibdtg,ibsec,iflen,islen)

    !     complete ddr

    !     preliminaries : file type

    ipf=1
    myddr%nlonhl=kb2(7)
    ipf=myddr%nlonhl

    if(debug)print *,' as2ddr, kdtg:',kdtg
    myddr%ndtvhl=ivdtg
    myddr%nscvhl=ivsec
    myddr%ndtbhl=ibdtg
    myddr%nscbhl=ibsec
    myddr%nflshl=islen + iflen*3600

    !     horizontal grid description section
    myddr%nlathl = kb2(9)

    !     vertical grid description section

    myddr%nlevhl=klev
    myddr%npplhl=2
    do j=2,klev
       myddr%nlpthl(j)=myddr%nlpthl(j-1)+ipf
    enddo

    !     multi-level fields description section
    myddr%nmlfhl=knmf
    myddr%nmpthl(1)=2*ipf
    do j=2,knmf
       myddr%nmpthl(j)=myddr%nmpthl(j-1)+ipf*klev
    enddo

    !     single-level fields description section
    myddr%nslfhl=knsf
    myddr%nspthl(2)=ipf
    if (knmf == 0) then
       myddr%nspthl(3)=myddr%nspthl(2)+ipf
    else
       myddr%nspthl(3)=myddr%nmpthl(knmf)+ipf*klev
    endif
    do j=4,knsf
       myddr%nspthl(j)=myddr%nspthl(j-1)+ipf
    enddo

    !     horizontal grid description section
    myddr%aplohl=real(kb2(36),realgribkind)*.001_realgribkind
    myddr%aplahl=real(kb2(33),realgribkind)*.001_realgribkind
    myddr%aweshl=real(kb2(14),realgribkind)*.001_realgribkind
    myddr%aeashl=real(kb2(21),realgribkind)*.001_realgribkind

    !     swap order if kb2(28)==64
    if(kb2(28)==64)then
       myddr%alafhl=real(kb2(11),realgribkind)*.001_realgribkind
       myddr%alalhl=real(kb2(18),realgribkind)*.001_realgribkind
    else
       myddr%alalhl=real(kb2(11),realgribkind)*.001_realgribkind
       myddr%alafhl=real(kb2(18),realgribkind)*.001_realgribkind
    endif
    myddr%dlonhl=real(real(myddr%aeashl-myddr%aweshl,realkind)/ &
         real(myddr%nlonhl-1,realkind),realgribkind)
    myddr%dlathl=real(real(myddr%alalhl-myddr%alafhl,realkind)/ &
         real(myddr%nlathl-1,realkind),realgribkind)

    !     vertical grid description section
    do j=1,klev
       myddr%alevhl(j,1)=pab(1,j)
       myddr%alevhl(j,2)=pab(2,j)
    enddo

    return

  end subroutine as2ddr


  subroutine getbvt(kbs1,kvdtg,kvsec,kbdtg,kbsec,klenhh,klensc)
    use util
    implicit none
    integer kvdtg,kvsec,kbdtg,kbsec,klenhh,klensc
    integer kbs1(*)

    !     getbvt:  get base and verification time of kbs1
    !     kvdtg=-1 if unknown time range indicator
    !     kvdtg=-2 if unknown unit of time range indicator
    !     
    !     gerard cats   knmi  4 november 1991

    integer iadd,isec

    !       accumulation time

    if(kbs1(21)==4)then
       iadd=kbs1(20)

       !     fail if not accumulation or analysis or forecast
    elseif (kbs1(21) > 2 .and. kbs1(21) /= 10) then
       kvdtg=-1
       return
       !       analysis/forecast
    else
       iadd=kbs1(19)
    endif

    !       convert time to add to seconds

    isec=0
    if(kbs1(18)==0)then

       isec=iadd*60
    elseif(kbs1(18)==1)then
       isec=iadd*3600
    elseif(kbs1(18)==2)then
       isec=iadd*24*3600
    elseif(kbs1(18)==254)then
       isec=iadd
    else
       kvdtg=-2
       return
    endif
    klenhh=0
    klensc=isec

    !     base time

    kbdtg=(kbs1(13)*100+kbs1(14))*100+kbs1(15)
    kbsec=(kbs1(16)*100+kbs1(17))*100

    !     verification time
    call adddtg(kbdtg, kbsec, klensc, kvdtg, kvsec)
  end subroutine getbvt


  subroutine setlof(kbs1,kl1,knfld,ksptr, &
       kdtg,kmfld,knlev,ksfld,kwmmhl,kwmshl,kslthl,plevhl)
    implicit none

    !     subroutine setlof
    !     purpose:
    !     create a list of fields from blocks 1, valid for time kdtg
    !     knlev is the number of levels, input
    !     ksptr is a pointer to the orography field.
    !     pass the other parameters to create comddr

    !     gerard cats   knmi  4 november 1991

    integer kl1,knfld,kdtg,kmfld,knlev,ksfld,ksptr,idum
    integer kbs1(kl1,knfld)
    integer kwmmhl(*),kwmshl(*),kslthl(*)
    real(kind=realgribkind)::plevhl(*)

    integer jfld,jm,ilev,iptr,ivs,ibd,ibs,idtg

    integer junit,ii
    logical::debug=.false.
!cau    logical::debug=.true.
    junit=6

    kmfld=0
    ksfld=2
    ksptr=0

    !     if dtg is not set, get it from the first orography field

    if(debug)print *, ' setlof knfld ',knfld
    if(debug) then
       do jfld=1,2
          write(6,6123) jfld,(kbs1(ii,jfld),ii=1,20)
       enddo
    endif
6123 format(i4,20I5)

    do 121 jfld=1,knfld
       if((kbs1( 9,jfld)== 6 .or. kbs1( 9,jfld)==81 .or.          &
            kbs1( 9,jfld)==86 .or. kbs1( 9,jfld)==11 ) .and.&
            (kbs1(11,jfld)==0 .or. kbs1(11,jfld)==1) .and.  &
            (kbs1(10,jfld)==102.or.kbs1(10,jfld)==105 .or.  &
            kbs1(10,jfld)==1) ) then
          call getbvt(kbs1(1,jfld),kdtg,ivs,ibd,ibs,idum,idum)
          if(debug)print *, ' setlof kdtg ',jfld,kdtg

          if ( kdtg == -1 ) then
             kdtg=kbs1(13,jfld)*10000+kbs1(14,jfld)*100+kbs1(15,jfld)
          endif
          goto 122
       endif
121 enddo

    !     no orography, try second field

    if(debug)print *, ' setlof second field'

    call getbvt(kbs1(1,min(2,knfld)),kdtg,ivs,ibd,ibs,idum,idum)
    if(debug)print *, ' kdtg ',kdtg

    ksptr=2
122 continue
    if(kdtg<=0)then
       print *,'no dtg found - abort ',kdtg
       call stop_program('')
    endif
    !     loop over fields

    do 701 jfld=1,knfld

       !     check time

       call getbvt(kbs1(1,jfld),idtg,ivs,ibd,ibs,idum,idum)

       !     a multi-level field

       if(kbs1(10,jfld)==109.or.kbs1(10,jfld)==0)then
          !     find its level number

          ilev=kbs1(11,jfld)
          if(ilev>knlev)goto 700

          !     check if found before

          do 421 jm=1,kmfld
             if(kwmmhl(jm)==kbs1(9,jfld))then
                goto 700
             endif
421       enddo
          !     new multi-level field

          kmfld=kmfld+1
          kwmmhl(kmfld)=kbs1(9,jfld)

          !     a single-level field

       else
          !     orography to start
          if((kbs1( 9,jfld)== 6 .or. kbs1( 9,jfld)==81 .or. &
               kbs1( 9,jfld) == 86 .or. kbs1( 9,jfld)==11 ) .and. &
               (kbs1(11,jfld)==0 .or. kbs1(11,jfld)==1 ).and. &
               (kbs1(10,jfld)==102.or.kbs1(10,jfld)==105 .or. &
               kbs1(10,jfld)==1) )then
             iptr=1
             ksptr =jfld
             !     surface pressure to follow
          elseif(kbs1( 9,jfld)== 1.and.kbs1(11,jfld)==0.and. &
               (kbs1(10,jfld)==102.or.kbs1(10,jfld)==105))then
             iptr=2
             !     other fields
          else
             ksfld=ksfld+1
             iptr=ksfld
          endif
          !     add

          kwmshl(iptr)=kbs1( 9,jfld)
          kslthl(iptr)=kbs1(10,jfld)
          plevhl(iptr)=real(kbs1(11,jfld),realgribkind)
          if(kslthl(iptr)==100)then
             plevhl(iptr)=plevhl(iptr)*100._realgribkind
          endif
       endif

700    continue
701 enddo
  end subroutine setlof

 subroutine readorogpath(orog_path,ilen)
    use decomp
    implicit none
    character(len=132)::orog_path
    integer ilen
    logical,save::done=.false.
    character(len=20)namelistfile,inst
    namelist/institute/inst
    namelist /orog/orog_path
    ilen=1
    open(66,file='namelists.dat',status='old',form='formatted')
    read(66,nml=institute)
    close(66)
    namelistfile='gcmpaths.'//trim(inst)
    open(61,file=namelistfile,status='old',form='formatted')
    read(61,nml=orog)
    close(61)
    if(mype.eq.0.and..not.done)then
       write(6,nml=orog)
       done = .true.
    endif
    ilen=index(orog_path,' ')-1
  end subroutine readorogpath


end module gcm
