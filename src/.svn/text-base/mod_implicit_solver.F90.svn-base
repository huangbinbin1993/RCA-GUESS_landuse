module mod_implicit_solver
  use RCAdomainMod
  use decomp
  implicit none
  private
  real(kind=realkind),save,allocatable,dimension(:)::cosu,cosv,sk,trigs

  real(kind=realkind),allocatable,dimension(:),save:: edepth     ! list of depths for vertical normal modes
  integer,save::nfax(10)      ! list of factors for the fourier transformation

  real(kind=realkind),allocatable,dimension(:,:),save::eiginv! eigenvector matrix for vertical coupling
  real(kind=realkind),allocatable,dimension(:,:),save::eig   ! eigenvector matrix for vertical decoupling
  real(kind=realkind),allocatable,dimension(:),save::c2      ! eigenvalues for normal modes
  real(kind=realkind),allocatable,dimension(:,:),save::gmat  ! vertical structure function g
  real(kind=realkind),allocatable,dimension(:,:),save::gmati ! inverse  vertical structure function
  real(kind=realkind),allocatable,dimension(:,:),save::pmat  ! vertical coupling matrix  gamma for the  auxillary variable p
  real(kind=realkind),allocatable,dimension(:,:),save::pmati ! inverse of  vertical coupling matrix

  !the following varaibles must be initialized before any call to sl2tim!
  real(kind=realkind),allocatable,dimension(:),public,save::sitau1 ! list of level dependent coefficients in the linear term in the thermodynamic equation
  real(kind=realkind),allocatable,dimension(:),public,save::sitau2 ! list of level dependent coefficients in the linear term in the thermodynamic equation
  real(kind=realkind),allocatable,dimension(:),public,save::sigam1 ! list of level dependent coefficients in the linear term in the hydrostatic equation
  real(kind=realkind),allocatable,dimension(:),public,save::sigam2 ! list of level dependent coefficients in the linear term in the hydrostatic equation
  real(kind=realkind),allocatable,dimension(:),public,save::sidpk0 ! list of level dependent coefficients in the linear term in the continuity equation

  real(kind=realkind),save,allocatable,dimension(:,:)::zhhdia 
  real(kind=realkind),save:: dtold=0._realkind

  public impini,hhsolv,hhsolvitr

contains
  subroutine hhsolv(klon, klat, klev, kpbpts, ptwodt, pdiv,  nlsimpx, epsg,&
       ahyb,bhyb)

    use fft_helmholtz
    use comtrix
    use transpose
    implicit none
    !      hhsolv  - subroutine for solution of helmholtz-equations
    !                in the semi-implicit timestepping scheme
    !
    !     j.e. haugen            hirlam
    !
    !     re-written to fit mpp-usage 12 july 1995 by
    !
    !     nils gustafsson        hirlam
    !     deborah salmond        cray resarch
    !
    !     purpose.
    !     --------
    !
    !     this routine solves the three-dimensional helmholtz-equation
    !     for the divergence by vertical decoupling into klev horizontal
    !     helmholtz-equations.
    !
    !     interface.
    !     hhsolv is called from the main time stepping routines

    !     input parameters.
    !     klon      number of gridpoints in x-direction
    !     klat      number of gridpoints in y-direction
    !     klev      number of vertical levels
    !     kpbpts    number of passive boundarylines
    !     ptwodt    current double timestep
    !     pdiv      right hand side of the hh-equations
    !
    !     output parameters:
    !     ------------------
    !
    !      pdiv       solution of hh-equations


    !     declaration of global parameters
    integer,intent(in)::klon,klat,klev,kpbpts
    real(kind=realkind),intent(in)::ptwodt
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    real(kind=realkind):: pdiv(klon,klat,klev)

    real(kind=realkind),allocatable,dimension(:,:):: div_fft,div_tri  
    real(kind=realkind),intent(in)::epsg 
    logical::nlsimpx
    integer:: i,j,k,jrec_fft,jrec_tri,jrec,kslab,jrec_min,jrec_max
    real(kind=realkind)::zdiv(klon,klat,klev)
    real(kind=realkind):: zwork(klon_global,klat_global*klev_global+1)
    real(kind=realkind):: z1peps,za,zaa,zcosv,zdt,c2(klev)
    integer::status
    
    allocate(div_fft(klon_global,klat_global*klev_global),stat=status)
    if(status/=0)then
       stop 'memory allocation error div_fft2'
    endif
    allocate(div_tri(klat_global,klon_global*klev_global),stat=status)
    if(status/=0)then
       stop 'memory allocation error div_tri2'
    endif
    if(.not.allocated(zhhdia))then
       allocate(zhhdia(klat_global,klon_global*klev_global+1),stat=status)
       if(status/=0)then
       stop 'memory allocation error zhhdia'
    endif
    endif
    !     pre-calculated tri-diagonal weights

    

    if(.not.allocated(cosu).and.nlsimpx)then
       call impini(klon, klat, klev, kpbpts,ahyb,bhyb)
    endif
    !    redo helmholz solver cecomposition if needed
!    if( klev/=klev_hh ) then 
!       call decompose_hh(klon,klat,klev,kpbpts)
!    endif

    !    vertical decoupling in sgemm
    if( klev>1 ) then
       call sgemm('N','T',klon*klat,klev,klev,1.0_realkind,pdiv,klon*klat, & 
            eiginv,klev, 0.0_realkind,zdiv,klon*klat)
    else
       do j=1,klat
          do i=1,klon
             zdiv(i,j,1) = pdiv(i,j,1)
          enddo
       enddo
    endif

    !    data transposition for fft's
    
    call twod_to_fft(zdiv,klon,klat,klev,div_fft)

    !    zero dtt divergence boundary condition
    jrec_fft = 0
    do kslab=1,nslab_fft
       j = jmin_fft(kslab)
       do jrec=1,jlen_fft(kslab)
          jrec_fft = jrec_fft + 1
          if(j==(kpbpts+1).or.j==(klat_global-kpbpts-1) .or. &
               j==(klat_global-kpbpts))then
             do i=kpbpts+1,klon_global-kpbpts-2
                div_fft(i,jrec_fft) = 0.0_realkind
             enddo
          endif
          if(j>=(kpbpts+2).and.j<=(klat_global-kpbpts-2))then
             div_fft(kpbpts+1,jrec_fft) = 0.0_realkind
          endif
          j = j + 1
       enddo
    enddo

    !    sin-transform of right hand side

    call fft44(div_fft(kpbpts+1,1), zwork(1,1),1,klon_global,klon_global-2*kpbpts-2,&
         jrec_fft,-1,nfax,trigs,klon_global)

    !    vertical part of helmholz term
    zdt    = 0.5_realkind*ptwodt
    z1peps = 1.0_realkind+epsg
    zdt    = z1peps*zdt
    zaa    = zdt*ffmean
    zaa    = 1.0_realkind+zaa*zaa


    if (nlsimpx) then
       do k=1,klev
          c2(k) = zaa/(ra*ra*rdth*rdth*edepth(k)*(zdt)**2.0_realkind)
       enddo
    else
       do k=1,klev
          c2(k) = ffmean*ffmean/(ra*ra*rdth*rdth*edepth(k))
       enddo
    endif

    !    data transposition for tri-diagonal solver
    call fft_to_tri(div_fft,kpbpts,klev,div_tri)

    !     loop over vertical slabs for tri-diagonal solver
    jrec_tri = 0

    do kslab=1,nslab_tri
       k        = jlev_tri(kslab)
       jrec_min = jrec_tri + 1
       jrec_max = jrec_tri + imax_tri(kslab) - imin_tri(kslab) + 1
       !     diagonals
       if( abs(dtold-ptwodt)>0.000001_realkind ) triwgt = -1
       if( triwgt<0 .or. klev==1 ) then
          jrec = jrec_tri
          do i=imin_tri(kslab),imax_tri(kslab)
             jrec = jrec + 1
             do j=kpbpts+2,klat_global-kpbpts-2
                zhhdia(j,jrec) = cosu(j)*c2(k) + ((2.0_realkind*rdlam)/(sk(i-kpbpts)*rdth))**2.0_realkind/cosu(j) &
                     + cosv(j) + cosv(j-1)
             enddo
             do j=kpbpts+3,klat_global-kpbpts-2
                zhhdia(j,jrec) = zhhdia(j,jrec) - cosv(j-1)*cosv(j-1)/zhhdia(j-1,jrec)
             enddo
             do j=kpbpts+2,klat_global-kpbpts-2
                zhhdia(j,jrec) = 1.0_realkind/zhhdia(j,jrec)
             enddo
          enddo
       endif


       !     backward elimination with order of loops to permit vectorization
       do jrec=jrec_min,jrec_max
          div_tri(kpbpts+1,jrec) = 0.0_realkind
       enddo

       za    = c2(k)*cosu(kpbpts+2)

       do jrec=jrec_min,jrec_max
          div_tri(kpbpts+2,jrec) = div_tri(kpbpts+2,jrec)*za
       enddo

       do j=kpbpts+3,klat_global-kpbpts-2
          za    = c2(k)*cosu(j)
          zcosv = cosv(j-1)
          do jrec=jrec_min,jrec_max
             div_tri(j,jrec) = div_tri(j,jrec)*za
             div_tri(j,jrec) = div_tri(j,jrec) + zcosv*zhhdia(j-1,jrec)*div_tri(j-1,jrec)
          enddo
       enddo

       do jrec=jrec_min,jrec_max
          div_tri(klat_global-kpbpts-1,jrec) = 0._realkind
       enddo


       !     forward substitution with order of loops to permit vectorization
       do j=klat_global-kpbpts-2,kpbpts+2,-1
          zcosv = cosv(j)
          do jrec=jrec_min,jrec_max
             div_tri(j,jrec) = (div_tri(j+1,jrec)*zcosv + div_tri(j,jrec))*zhhdia(j,jrec)
          enddo
       enddo
       do jrec=jrec_min,jrec_max
          div_tri(kpbpts+1,jrec) = 0.0_realkind
       enddo
       !    endloop over vertical slabs
       jrec_tri = jrec_max
    enddo

    !    data transposition for fft

    call tri_to_fft(div_tri,kpbpts,klev,div_fft)

    !    inverse fourier transform

    call fft44(div_fft(kpbpts+1,1),zwork(1,1),1,klon_global,klon_global-2*kpbpts-2, &
         jrec_fft,+1,nfax,trigs,klon_global)

    !    boundary conditions
    jrec_fft = 0

    do kslab=1,nslab_fft
       j = jmin_fft(kslab)
       do jrec=1,jlen_fft(kslab)
          jrec_fft = jrec_fft + 1
          if(j==(kpbpts+1).or.j==(klat_global-kpbpts-1) .or. j==(klat_global-kpbpts))then
             do i=kpbpts+1,klon_global-kpbpts-2
                div_fft(i,jrec_fft) = 0.0_realkind
             enddo
          endif

          if( j>=(kpbpts+1).and. j<=(klat_global-kpbpts))then
             div_fft(kpbpts+1,jrec_fft)             = 0.0_realkind
             div_fft(klon_global-kpbpts-1,jrec_fft) = 0.0_realkind
             div_fft(klon_global-kpbpts,jrec_fft)   = 0.0_realkind
          endif
          j = j + 1
       enddo
    enddo

    !    data transposition back to 2d

    call fft_to_twod(div_fft,klon,klat,klev,zdiv)

    !    vertical coupling
    if( klev>1 ) then
       call sgemm('N','T',klon*klat,klev,klev, 1.0_realkind,zdiv,klon*klat,eig,klev,0.0_realkind,pdiv,klon*klat)
    else
       do j=1,klat
          do i=1,klon
             pdiv(i,j,1) = zdiv(i,j,1)
          enddo
       enddo
    endif

    triwgt = +1
    dtold  = ptwodt

    return

  end subroutine hhsolv



  subroutine impini(klon,klat,klev,kpbpts,ahyb,bhyb)
    !     
    ! IMPINI - INITIALIZATION OF IMPLICIT SCHEMES
    !     
    !     J.E. HAUGEN             HIRLAM
    !     
    !     PURPOSE.
    !     --------
    !     
    !     SET UP TIME INDEPENDENT CONSTANTS USED IN THE SEMI-IMPLICIT
    !     TIME STEPPING SCHEME AND IN THE IMPLICIT NONLINEAR NORMAL
    !     MODE INITIALIZATION SCHEME. THE CONSTANTS ARE STORED IN
    !     COMMON BLOKCS COMIMP, COMNMI AND COMFFT.
    !     
    !   INTERFACE.
    !     ----------
    !     
    !     IMPINI IS CALLED FROM GEMINI
    !     
    !     INPUT PARAMETERS:
    !     -----------------
    !     
    !     KLON      NUMBER OF GRIDPOINTS IN X-DIRECTION
    !     KLAT      NUMBER OF GRIDPOINTS IN Y-DIRECTION
    !     KLEV      NUMBER OF VERTICAL LEVELS
    !     KPBPTS    NUMBER OF PASSIVE BOUNDARY LINES
    !     KLAM      CHOISE OF LAM-FORMULATION
    !     KLAM=1 : HIRLAM HYBRID COORDINATES, LEVEL 1
    !     KLAM=2 : HIRLAM HYBRID COORDINATES, LEVEL 2
    !     PTWODT    DOUBLE TIME-STEP 2.0DT
    !     
    !     EXTERNALS.
    !     ----------
    !     
    !     MATIN1    MATRIX INVERSION ROUTINE
    !     GETEIG    EIGEN VALUES AND EIGEN VECTOR ROUTINE
    !     FAX       PRIME FACTORS FOR FFT-ROUTINES
    !     FFTRIG    TRIGONOMETRIC FUNCTIONS FOR FFT-ROUTINES

    use comhkp
    use referenceParameters

    implicit none

    integer,intent(in):: klon, klat , klev , kpbpts 
    real(kind=realkind),intent(in)::ahyb(klev+1),bhyb(klev+1)
    
    real(kind=8),allocatable::wdiv(:),wtemp(:),wpot(:),wlnps(:)
    real(kind=8),allocatable::wsum(:) 
    real(kind=8):: zrgasd,zcappa,zcapt0,zp0m,zp0p,zdlnp0, zalfa0 
    real(kind=8):: zalfa1 ,ztau1 , ztau2  , zdpk0, zrp0  , zrt0 
    real(kind=8):: zgam1  ,zgam2 , zdet   ,zpi   , zpir18   
    real(kind=8):: zsc   , zdth   ,zhxumm, zhyvmm ,zpin

    real(kind=realkind),allocatable,dimension(:,:):: aa,g
    integer,allocatable::index(:) 
    integer jk,ilev,isqr,idum,ierr,jl,jlev,jy,jx,ilnm2

    if(.not.allocated(cosu))then
       allocate(cosu(klat_global), cosv(klat_global) )
       allocate(sk(klon_global)  , trigs(klon_global/2*5))
       allocate(edepth(klev))
       allocate(eiginv(klev,klev))
       allocate(eig(klev,klev))
       allocate(c2(klev))
       allocate(gmat(klev,klev),gmati(klev,klev),pmat(klev,klev))
       allocate(pmati(klev,klev))
       allocate(sitau1(klev),sitau2(klev),sigam1(klev))
       allocate(sigam2(klev),sidpk0(klev))
    endif

    allocate(aa(klev,klev),   g(klev,klev))
    allocate(wdiv(klev*klev), wtemp(klev*klev)  ,&
         wpot(klev*klev), wlnps(klev)  ,  wsum(klev))
    allocate(index(klev))



    !         CHECK DIMENSION OF WORK SPACE ARRAYS
    if(mype==0) print *,'---------- in impini ----------'


    !         GENERATE CONSTANTS FOR REFERENCE STATE
    zrgasd = 287.04_realkind
    zcappa = 0.2857143_realkind
    zcapt0 = zcappa * sit0

    zp0m = ahyb(1) + bhyb(1)*sip0
    do  jk=1,klev
       zp0p       = ahyb(jk+1) + bhyb(jk+1)*sip0
       sidpk0(jk) = zp0p - zp0m
       if (jk==1) then
          zdlnp0 = 0.0_realkind
          zalfa0 = log(2.0_realkind)
          zalfa1 = 1.0_realkind
       else
          zdlnp0 = log( zp0p/zp0m )
          zalfa0 = 1._realkind - zp0m*zdlnp0/sidpk0(jk)
          zalfa1 = zalfa0
       endif
       sitau1(jk) = zcapt0 * zdlnp0 / sidpk0(jk) !used in sl2tim
       sitau2(jk) = zcapt0 * zalfa1              !used in sl2tim
       sigam1(jk) = zrgasd * zalfa0              !used in sl2tim
       sigam2(jk) = zrgasd * zdlnp0              !used in sl2tim
       zp0m       = zp0p
    enddo


    !         COMPUTE VERTICAL STRUCTURE MATRIX G
    !         SET MATRIX OF DIVERGENCE FIELDS
    do  jk =1,klev
       ilev = (jk-1)*klev
       do  jl=1,klev
          wdiv(jl+ilev) = 0.0_realkind
       enddo
       wdiv(jk+ilev) = 1.0_realkind
    enddo

    !         TEMPERATURE AND SURFACE PRESSURE TERMS
    do  jl=1,klev
       wsum(jl) = 0._realkind
    enddo

    do jk=1,klev
       ilev  = (jk-1)*klev
       ztau1 = sitau1(jk)
       ztau2 = sitau2(jk)
       zdpk0 = sidpk0(jk)
       do  jl=1,klev
          wtemp(jl+ilev) = ztau1*wsum(jl) + ztau2*wdiv(jl+ilev)
          wsum(jl)      =       wsum(jl) + zdpk0*wdiv(jl+ilev)
       enddo
    enddo

    zrp0 = 1._realkind/sip0
    do  jl=1,klev
       wlnps(jl) = zrp0*wsum(jl)
    enddo

    !         pressure gradient term
    zrt0 = zrgasd*sit0

    do  jl=1,klev
       wsum(jl) = zrt0 * wlnps(jl)
    enddo

    do  jk=klev,1,-1
       zgam1 = sigam1(jk)
       zgam2 = sigam2(jk)
       do  jl=1,klev
          wpot(jl+(jk-1)*klev) = wsum(jl) + zgam1*wtemp(jl+(jk-1)*klev)
          wsum(jl)      = wsum(jl) + zgam2*wtemp(jl+(jk-1)*klev)
       enddo
    enddo

    do jk=1,klev
       do  jl=1,klev
          g(jl,jk) = wpot(jk+(jl-1)*klev)
       enddo
    enddo

    !         COMPUTE INVERSE OF G-MATRIX IN *MATIN1
    gmat = g
    aa   = g

    isqr = klev*klev

    call matin1(aa,klev,klev,idum,0,index,ierr,zdet,isqr)

    gmati = aa

    !         COMPUTE MATRIX FOR PRESSURE GRADIENT TERM

    do  jk=1,klev
       ilev = (jk-1)*klev
       do  jl=1,klev
          if (jl<jk) wpot(jl+ilev) = 0.0_realkind
          if (jl==jk) wpot(jl+ilev) = sigam1(jl)
          if (jl>jk) wpot(jl+ilev) = sigam2(jl)
       enddo
    enddo

    do  jk=1,klev
       !       ilev = (jk-1)*klev
       do  jl=1,klev
          jlev = (jl-1)*klev
          !        pmat(jl+(jk-1)*klev) = wpot(jk+jlev)
          pmat(jl,jk) = wpot(jk+jlev)
       enddo
    enddo

    !         compute inverse of pmat in matin1
    aa = pmat
    isqr = klev*klev

    call matin1(aa,klev,klev,idum,0,index,ierr,zdet,isqr)

    pmati = aa

    !         eigen values, eigen vectors of g and invers eigen vector matrix
    call geteig(g,eig,eiginv,c2,aa,klev)

    do  jk=1,klev
       edepth(jk) = c2(jk)
    enddo
    !         cosu and cosv contain cos(theta) in reversed order
    zpi    = 2._realkind*asin(1._realkind)
    zpir18 = zpi/180.0_realkind

    zsc    = RCAdomain%south*zpir18
    zdth   = RCAdomain%dlat*zpir18


    do  jy=1,klat_global
       cosu(jy) = cos( zsc +  real(jy-1,realkind)   *zdth )
       cosv(jy) = cos( zsc + (real(jy-1,realkind)+0.5_realkind)*zdth )
    enddo


    !         initialise sine-coefficients
    ilnm2 = klon_global - 2 - 2*kpbpts
    zpin  = asin(1._realkind)/real(ilnm2,realkind)

    sk(1) = 0.0_realkind
    do  jl=2,ilnm2
       sk(jl) = 1.0_realkind/(sin(real(jl-1,realkind)*zpin))
    enddo

    !         initialize  comfft
    call fax(nfax, ilnm2, 4)
    call fftrig(trigs, ilnm2, 4)

    if(mype==0) print *,' nfax =',nfax


    if(mype==0)print *,'---------- impini done ----------'
    deallocate(aa,   g)
    deallocate(wdiv, wtemp,wpot, wlnps,wsum)
    deallocate(index)
    return
  end subroutine impini

  subroutine fftrig (trigs, n, mode)  

    implicit none  
    real(kind=realkind) :: trigs ( * )  
    real(kind=realkind) :: pi, del, angle  


    integer :: imode, mode, nn, n, l, i, nh, la  
    pi = 2.0_realkind * asin (1.0_realkind)  
    imode = abs (mode)  
    nn = n  
    if (imode>1.and.imode<6) nn = n / 2  
    del = (pi + pi) / real (nn,realkind)  

    l = nn + nn  
    do i = 1, l, 2  
       angle = 0.5_realkind * real (i - 1,realkind) * del  
       trigs (i) = cos (angle)  
       trigs (i + 1) = sin (angle)  

    enddo
    if (imode==1) return  

    if (imode==8) return  
    del = 0.5_realkind * del  
    nh = (nn + 1) / 2  
    l = nh + nh  

    la = nn + nn  
    do i = 1, l, 2  
       angle = 0.5_realkind * real (i - 1,realkind) * del  
       trigs (la + i) = cos (angle)  
       trigs (la + i + 1) = sin (angle)  

    enddo

    if (imode<=3) return  
    del = 0.5_realkind * del  

    la = la + nn  
    if (mode/=5) then  
       do i = 2, nn  
          angle = real(i - 1,realkind) * del  
          trigs (la + i) = 2.0_realkind * sin (angle)  
       enddo
       return  


    endif

    del = 0.5_realkind * del  
    do 50 i = 2, n  
       angle = real (i - 1,realkind) * del  
       trigs (la + i) = sin (angle)  
50  enddo
    return  
  end subroutine fftrig

  subroutine fax(ifax, n, mode)  

    implicit none  
    integer :: n, mode, nn, k, l, inc, nfaxL, ii, istop, i, item  

    integer :: ifax (10)  
    nn = n  
    if (.not. (abs (mode) ==1.or.abs (mode) ==8) ) then  
       nn = n / 2  
       if (.not. ( (nn + nn) ==n) ) then  
          ifax (1) = - 99  
          return  
       endif

    endif
    k = 1  
    l = 5  

    inc = 2  

    do while (nn/=1)  
       !test for factors of 4
       do while (mod (nn, 4) ==0.and.nn/=1)  
          k = k + 1  
          ifax (k) = 4  
          nn = nn / 4  

       enddo
       !test for factor of 2
       do while (mod (nn, 2) ==0.and.nn/=1)  
          k = k + 1  
          ifax (k) = 2  
          nn = nn / 2  

       enddo
       !test for factors of 3
       do while (mod (nn, 3) ==0.and.nn/=1)  
          k = k + 1  
          ifax (k) = 3  
          nn = nn / 3  

       enddo
       !5,7,9,11,13,15,
       do while (mod (nn, l) ==0.and.nn/=1)  
          k = k + 1  
          ifax (k) = l  
          nn = nn / l  
       enddo
       l = l + inc  
       inc = 6 - inc  


    enddo

    ifax (1) = k - 1  
    !     ifax(1) contains number of factors
    nfaxL = ifax (1)  
    !     sort factors into ascending order
    if (nfaxL/=1) then  
       do ii = 2, nfaxL  
          istop = nfaxL + 2 - ii  
          do i = 2, istop  
             if (ifax (i + 1) <ifax (i) ) then  
                item = ifax (i)  
                ifax (i) = ifax (i + 1)  
                ifax (i + 1) = item  
             endif
          enddo
       enddo


    endif
    return  
  end subroutine fax

  subroutine geteig(am,evects,einv,eigval,atest,nmax)

    !     FIND THE EIGENVALUES AND EIGENVECTORS
    !     WHICH DESCRIBES THE VERTICAL STRUCTURE OF THE MODEL ATMOSPHERE.

    !  use eisrg1__genmod
    !  use matin1__genmod
    implicit none

    integer,intent(in)::nmax
    real(kind=realkind)::am(:,:), evects(:,:), einv(:,:)
    real(kind=realkind)::eigval(:) , atest (:,:)

    !    DECLARATION OF LOCAL WORKSPACE
    integer:: il,jl,j,isqr,idum,ierr,jk
    real(kind=8):: zdet
    real(kind=realkind)::  evalsr(nmax), evalsi(nmax) !automatic allocation
    integer:: index(nmax)                             !automatic allocation

    if(ubound(am,1)/=nmax)then
       print *,'geteig bug',ubound(am,1),nmax
       call stop_program('')
    endif
    if(ubound(am,2)/=nmax)then
       print *,'geteig bug',ubound(am,2),nmax
       call stop_program('')
    endif

    !     WRITE(6,'(/,1X,''IN GETEIG, NMAX='',I9)') NMAX
    do il=1,nmax
       do jl=1,nmax
          atest(il,jl) = am(il,jl)
       enddo
    enddo

    !    CALL THE EISPAC ROUTINE FOR REAL MATRICES.
    call eisrg1(nmax,atest,evalsr,evalsi,evects,ierr)

    do j=1,nmax
       eigval(j) = evalsr(j)
    enddo

    if (ierr/=0) then
       write(6,'(1x,''error in eisrg1'')')
       call stop_program('')
    endif

    !    PRINT EIGEN VALUES AND EIGEN-VECTORS
    !  write(6,'(/,1x,''eigval='',10(6e10.3,/,8x))')(eigval(j),j=1,nmax)

    !  do jk=1,nmax
    !     write(6,'(/,1x,''eigvec='',10(6e10.3,/,8x))')(evects(j,jk),j=1,nmax)
    !  enddo

    !    CALCULATE THE INVERSE OF EVECTS.

    do  il=1,nmax
       do  jl=1,nmax
          atest(il,jl) = evects(il,jl)
       enddo
    enddo

    isqr=nmax*nmax

    call matin1(atest,nmax,nmax,idum,0,index,ierr,zdet,isqr)

    do  il=1,nmax
       do  jl=1,nmax
          einv(il,jl) = atest(il,jl)
       enddo
    enddo
    
#ifdef DEBUG    
    write(6,'(/,1x,''geteig done'')')
#endif
    return
  end subroutine geteig

  subroutine matin1(aa, idim1, n1, idim2, n2, index, nerror,determ, nsqr)
    implicit none  
    integer :: nsqr, i, n, n1, iemat, n2, idim, nmin1, IPIVC, main  
    integer :: ipivc1, ipivc2, i1, lpiv, icol, i3, jcol  


    integer :: i2, nerr, nerror, idim1, idim2  
    !     MATRIX INVERSION WITH ACCOMPANYING SOLUTION OF LINEAR EQUATIONS
    !     CERN LIBRARY
    real(kind=8):: deter, pivot, swap, determ  
    real(kind=8):: a (nsqr)  
    real(kind=realkind) :: aa (nsqr)  
    integer :: index (idim1)  
    do i = 1, nsqr  
       a (i) = dble (aa (i) )  
    enddo
    deter = 1.0d0  
    n = n1  
    iemat = n + n2  
    idim = idim1  

    nmin1 = n - 1  
    !     THE ROUTINE DOES ITS OWN EVALUATION FOR DOUBLE SUBSCRIPTING OF
    !     ARRAY A.
    ipivc = 1 - idim  
    !     MAIN LOOP TO INVERT THE MATRIX
    do 11 main = 1, n  
       pivot = 0.0d0  
       ipivc = ipivc + idim  
       !     SEARCH FOR NEXT PIVOT IN COLUMN MAIN.
       ipivc1 = ipivc + main - 1  
       ipivc2 = ipivc + nmin1  
       do  i1 = ipivc1, ipivc2  
          if ( (dabs (dble (a (i1) ) ) - dabs (pivot) ) <=0._realkind) then  
             goto 2  
          else  
             goto 1  
          endif
1         pivot = a (i1)  
          lpiv = i1  
2         continue
       enddo
       !     IS PIVOT DIFFERENT FROM ZERO
       if (abs(pivot)<0.d-14) then  
          goto 15  
       else  
          goto 3  
       endif
       !     GET THE PIVOT-LINE INDICATOR AND SWAP LINES IF NECESSARY
3      icol = lpiv - ipivc + 1  
       index (main) = icol  
       if ( (icol - main) <=0) then  
          goto 6  
       else  
          goto 4  
       endif
       !     COMPLEMENT THE DETERMINANT
4      deter = - deter  
       !     POINTER TO LINE PIVOT FOUND

       icol = icol - idim  
       !     POINTER TO EXACT PIVOT LINE
       i3 = main - idim  
       do 5 i = 1, iemat  
          icol = icol + idim  
          i3 = i3 + idim  
          swap = a (i3)  
          a (i3) = a (icol)  

          a (icol) = swap  
5      enddo
       !     COMPUTE DETERMINANT
6      deter = deter * pivot  
       pivot = 1d0 / pivot  
       !     TRANSFORM PIVOT COLUMN
       i3 = ipivc + nmin1  
       do 7 i = ipivc, i3  
          a (i) = - a (i) * pivot  
7      enddo
       a (ipivc1) = pivot  
       !     PIVOT ELEMENT TRANSFORMED
       !
       !     NOW CONVERT REST OF THE MATRIX
       i1 = main - idim  
       !     POINTER TO PIVOT LINE ELEMENTS
       icol = 1 - idim  
       !     GENERAL COLUMN POINTER
       do i = 1, iemat  
          icol = icol + idim  
          i1 = i1 + idim  
          !     POINTERS MOVED
          if ( (i - main) ==0) then  
             goto 10  
          else  
             goto 8  
          endif
          !     PIVOT COLUMN EXCLUDED
8         jcol = icol + nmin1  
          swap = a (i1)  
          i3 = ipivc - 1  
          do 9 i2 = icol, jcol  
             i3 = i3 + 1  
             a (i2) = a (i2) + swap * a (i3)  
9         enddo
          a (i1) = swap * pivot  
10        continue
       end do

11  end do
    !     NOW REARRANGE THE MATRIX TO GET RIGHT INVERS
    do  i1 = 1, n  
       main = n + 1 - i1  
       lpiv = index (main)  
       if ( (lpiv - main) ==0) then  
          goto 14  
       else  
          goto 12  
       endif
12     icol = (lpiv - 1) * idim + 1  
       jcol = icol + nmin1  
       ipivc = (main - 1) * idim + 1 - icol  
       do 13 i2 = icol, jcol  
          i3 = i2 + ipivc  
          swap = a (i2)  
          a (i2) = a (i3)  
          a (i3) = swap  
13     enddo
14     continue
    end do
    nerr = 0  
    goto 16  
15  nerror = main  
    goto 16  
16  continue  
    if (abs (deter) <1.0d-38) then  
       determ = sign (1.0d-38, deter)  
    elseif (abs (deter) >1.0d+38) then  
       determ = sign (1.0d+38, deter)  
    else  
       determ = deter  
    endif
    do i = 1, nsqr  
       aa (i) = real (a (i),realkind )  
    enddo
  END subroutine matin1

  subroutine elmhes(nm, n, low, igh, a, intvar)  
    implicit none  
    integer :: i, j, m, n, la, nm, igh, kp1, low, mm1, mp1  
    real(kind=realkind) :: a (nm, n)  
    real(kind=realkind) :: x, y  
    integer :: intvar (igh)  
    !
    !     this subroutine is a translation of the algol procedure elmhes,
    !     num. math. 12, 349-368(1968) by martin and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
    !
    !     given a real general matrix, this subroutine
    !     reduces a submatrix situated in rows and columns
    !     low through igh to upper hessenberg form by
    !     stabilized elementary similarity transformations.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrix.
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been use d,
    !          set low=1, igh=n.
    !
    !        a contains the input matrix.
    !
    !     on output
    !
    !        a contains the hessenberg matrix.  the multipliers
    !          which were use d in the reduction are stored in the
    !          remaining triangle under the hessenberg matrix.
    !
    !        intvar contains information on the rows and columns
    !          interchanged in the reduction.
    !          only elements low through igh are use d.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    la = igh - 1  
    kp1 = low + 1  
    if (la<kp1) goto 200  
    !
    do m = kp1, la  
       mm1 = m - 1  
       x = 0.0_realkind  
       i = m  
       !
       do j = m, igh  
          if (abs (a (j, mm1) ) <=abs (x) ) goto 100  
          x = a (j, mm1)  
          i = j  
100       continue
       end do
       !
       intvar (m) = i  
       if (i==m) goto 130  
       !     .......... interchange rows and columns of a ..........
       do 110 j = mm1, n  
          y = a (i, j)  
          a (i, j) = a (m, j)  
          a (m, j) = y  
110    end do
       !
       do 120 j = 1, igh  
          y = a (j, i)  
          a (j, i) = a (j, m)  
          a (j, m) = y  
120    end do
       !     .......... end interchange ..........
130    if (abs(x)<1.e-14_realkind) goto 180  

       mp1 = m + 1  
       do i = mp1, igh  
          y = a (i, mm1)  
          if (abs(y)<1.e-14_realkind) goto 160  
          y = y / x  

          a (i, mm1) = y  
          do 140 j = m, n  
             a (i, j) = a (i, j) - y * a (m, j)  
140       enddo
          do 150 j = 1, igh  
             a (j, m) = a (j, m) + y * a (j, i)  
150       enddo
160       continue
       end do
180    continue
    end do
200 return  
  end subroutine elmhes

  subroutine eltran (nm, n, low, igh, a, intvar, z)  
    implicit none  
    integer :: i, j, n, kl, mm, mp, nm, igh, low, mp1  
    real(kind=realkind) :: a (nm, igh), z(nm, n)  
    integer :: intvar(igh)  
    !
    !     this subroutine is a translation of the algol procedure elmtrans,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine accumulates the stabilized elementary
    !     similarity transformations use d in the reduction of a
    !     real general matrix to upper hessenberg form by  elmhes.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrix.
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been use d,
    !          set low=1, igh=n.
    !
    !        a contains the multipliers which were use d in the
    !          reduction by  elmhes  in its lower triangle
    !          below the subdiagonal.
    !
    !        intvar contains information on the rows and columns
    !          interchanged in the reduction by  elmhes.
    !          only elements low through igh are use d.
    !
    !     on output
    !
    !        z contains the transformation matrix produced in the
    !          reduction by  elmhes.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    !     .......... initialize z to identity matrix ..........
    do 80 j = 1, n  
       do 60 i = 1, n  
          z (i, j) = 0.0_realkind  
60     enddo
       z (j, j) = 1.0_realkind  
80  end do

    kl = igh - low - 1  
    if (kl<1) goto 200  
    !     .......... for mp=igh-1 step -1 until low+1 do -- ..........
    do mm = 1, kl  
       mp = igh - mm  
       mp1 = mp + 1  
       !
       do 100 i = mp1, igh  
          z (i, mp) = a (i, mp - 1)  
100    enddo
       !
       i = intvar(mp)  
       if (i==mp) goto 140  
       !
       do 130 j = mp, igh  
          z (mp, j) = z (i, j)  
          z (i, j) = 0.0_realkind  
130    end do
       !
       z (i, mp) = 1.0_realkind  
140    continue
    end do
    !
200 return  
  end subroutine eltran

  subroutine hqr2 (nm, n, low, igh, h, wr, wi, z, ierr)  
    implicit none  
    integer :: i, j, k, l, m, n, en, ii, jj, ll, mm, na, nm, nn, igh, &
         itn, its, low, mp2, enm2, ierr
    real(kind=realkind) :: h (nm, n), wr (n), wi (n), z (nm, n)  
    real(kind=realkind) :: p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, tst1, &
         tst2
    logical :: notlas  
    !
    !     this subroutine is a translation of the algol procedure hqr2,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine finds the eigenvalues and eigenvectors
    !     of a real upper hessenberg matrix by the qr method.  the
    !     eigenvectors of a real general matrix can also be found
    !     if  elmhes  and  eltran  or  orthes  and  ortran  have
    !     been use d to reduce this general matrix to hessenberg form
    !     and to accumulate the similarity transformations.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrix.
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been use d,
    !          set low=1, igh=n.
    !
    !        h contains the upper hessenberg matrix.
    !
    !        z contains the transformation matrix produced by  eltran
    !          after the reduction by  elmhes, or by  ortran  after the
    !          reduction by  orthes, if performed.  if the eigenvectors
    !          of the hessenberg matrix are desired, z must contain the
    !          identity matrix.
    !
    !     on output
    !
    !        h has been destroyed.
    !
    !        wr and wi contain the real and imaginary parts,
    !          respectively, of the eigenvalues.  the eigenvalues
    !          are unordered except that complex conjugate pairs
    !          of values appear consecutively with the eigenvalue
    !          having the positive imaginary part first.  if an
    !          error exit is made, the eigenvalues should be correct
    !          for indices ierr+1,...,n.
    !
    !        z contains the real and imaginary parts of the eigenvectors.
    !          if the i-th eigenvalue is real, the i-th column of z
    !          contains its eigenvector.  if the i-th eigenvalue is complex
    !          with positive imaginary part, the i-th and (i+1)-th
    !          columns of z contain the real and imaginary parts of its
    !          eigenvector.  the eigenvectors are unnormalized.  if an
    !          error exit is made, none of the eigenvectors has been found.
    !
    !        ierr is set to
    !          zero       for normal return,
    !          j          if the limit of 30*n iterations is exhausted
    !                     while the j-th eigenvalue is being sought.
    !
    !     calls cdiv for complex division.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    ierr = 0  
    norm = 0.0_realkind  
    k = 1  
    !     .......... store roots isolated by balanc
    !                and compute matrix norm ..........
    do  i = 1, n  
       !
       do 40 j = k, n  
          norm = norm + abs (h (i, j) )  
40     enddo
       !
       k = i  
       if (i>=low.and.i<=igh) goto 50  
       wr (i) = h (i, i)  
       wi (i) = 0.0_realkind  
50     continue
    enddo
    !
    en = igh  
    t = 0.0_realkind  
    itn = 30 * n  
    !     .......... search for next eigenvalues ..........
60  if (en<low) goto 340  
    its = 0  
    na = en - 1  
    enm2 = na - 1  
    !     .......... look for single small sub-diagonal element
    !                for l=en step -1 until low do -- ..........
70  do 80 ll = low, en  
       l = en + low - ll  
       if (l==low) goto 100  
       s = abs (h (l - 1, l - 1) ) + abs (h (l, l) )  
       if (abs(s)<1.e-14_realkind) s = norm  
       tst1 = s  
       tst2 = tst1 + abs (h (l, l - 1) )  
       if (abs(tst2-tst1)<1.e-14_realkind) goto 100  
80  enddo
    !     .......... form shift ..........
100 x = h (en, en)  
    if (l==en) goto 270  
    y = h (na, na)  
    w = h (en, na) * h (na, en)  
    if (l==na) goto 280  
    if (itn==0) goto 1000  
    if (its/=10.and.its/=20) goto 130  
    !     .......... form exceptional shift ..........
    t = t + x  
    !
    do 120 i = low, en  
       h (i, i) = h (i, i) - x  
120 enddo
    !
    s = abs (h (en, na) ) + abs (h (na, enm2) )  
    x = 0.75_realkind * s  
    y = x  
    w = - 0.4375_realkind * s * s  
130 its = its + 1  
    itn = itn - 1  
    !     .......... look for two consecutive small
    !                sub-diagonal elements.
    !                for m=en-2 step -1 until l do -- ..........
    do 140 mm = l, enm2  
       m = enm2 + l - mm  
       zz = h (m, m)  
       r = x - zz  
       s = y - zz  
       p = (r * s - w) / h (m + 1, m) + h (m, m + 1)  
       q = h (m + 1, m + 1) - zz - r - s  
       r = h (m + 2, m + 1)  
       s = abs (p) + abs (q) + abs (r)  
       p = p / s  
       q = q / s  
       r = r / s  
       if (m==l) goto 150  
       tst1 = abs (p) * (abs (h (m - 1, m - 1) ) + abs (zz) + abs (h ( &
            m + 1, m + 1) ) )
       tst2 = tst1 + abs (h (m, m - 1) ) * (abs (q) + abs (r) )  
       if (abs(tst2-tst1)<1.e-15_realkind) goto 150  
140 enddo
    !
150 mp2 = m + 2  
    !
    do  i = mp2, en  
       h (i, i - 2) = 0.0_realkind  
       if (i==mp2) goto 160  
       h (i, i - 3) = 0.0_realkind  
160    continue
    enddo
    !     .......... double qr step involving rows l to en and
    !                columns m to en ..........
    do  k = m, na  
       notlas = k/=na  
       if (k==m) goto 170  
       p = h (k, k - 1)  
       q = h (k + 1, k - 1)  
       r = 0.0_realkind  
       if (notlas) r = h (k + 2, k - 1)  
       x = abs (p) + abs (q) + abs (r)  
       if (abs(x)<1.e-15_realkind) goto 260  
       p = p / x  
       q = q / x  
       r = r / x  
170    s = sign (sqrt (p * p + q * q + r * r), p)  
       if (k==m) goto 180  
       h (k, k - 1) = - s * x  
       goto 190  
180    if (l/=m) h (k, k - 1) = - h (k, k - 1)  
190    p = p + s  
       x = p / s  
       y = q / s  
       zz = r / s  
       q = q / p  
       r = r / p  
       if (notlas) goto 225  
       !     .......... row modification ..........
       do 200 j = k, n  
          p = h (k, j) + q * h (k + 1, j)  
          h (k, j) = h (k, j) - p * x  
          h (k + 1, j) = h (k + 1, j) - p * y  
200    enddo
       !
       j = min (en, k + 3)  
       !     .......... column modification ..........
       do 210 i = 1, j  
          p = x * h (i, k) + y * h (i, k + 1)  
          h (i, k) = h (i, k) - p  
          h (i, k + 1) = h (i, k + 1) - p * q  
210    enddo
       !     .......... accumulate transformations ..........
       do 220 i = low, igh  
          p = x * z (i, k) + y * z (i, k + 1)  
          z (i, k) = z (i, k) - p  
          z (i, k + 1) = z (i, k + 1) - p * q  
220    enddo
       goto 255  
225    continue  
       !     .......... row modification ..........
       do 230 j = k, n  
          p = h (k, j) + q * h (k + 1, j) + r * h (k + 2, j)  
          h (k, j) = h (k, j) - p * x  
          h (k + 1, j) = h (k + 1, j) - p * y  
          h (k + 2, j) = h (k + 2, j) - p * zz  
230    enddo
       !
       j = min (en, k + 3)  
       !     .......... column modification ..........
       do 240 i = 1, j  
          p = x * h (i, k) + y * h (i, k + 1) + zz * h (i, k + 2)  
          h (i, k) = h (i, k) - p  
          h (i, k + 1) = h (i, k + 1) - p * q  
          h (i, k + 2) = h (i, k + 2) - p * r  
240    enddo
       !     .......... accumulate transformations ..........
       do 250 i = low, igh  
          p = x * z (i, k) + y * z (i, k + 1) + zz * z (i, k + 2)  
          z (i, k) = z (i, k) - p  
          z (i, k + 1) = z (i, k + 1) - p * q  
          z (i, k + 2) = z (i, k + 2) - p * r  
250    enddo
255    continue  
       !
260    continue
    enddo
    !
    goto 70  
    !     .......... one root found ..........
270 h (en, en) = x + t  
    wr (en) = h (en, en)  
    wi (en) = 0.0_realkind  
    en = na  
    goto 60  
    !     .......... two roots found ..........
280 p = (y - x) / 2.0_realkind  
    q = p * p + w  
    zz = sqrt (abs (q) )  
    h (en, en) = x + t  
    x = h (en, en)  
    h (na, na) = y + t  
    if (q<0.0_realkind) goto 320  
    !     .......... real pair ..........
    zz = p + sign (zz, p)  
    wr (na) = x + zz  
    wr (en) = wr (na)  
    if (abs(zz)>1.e-14_realkind) wr (en) = x - w / zz  
    wi (na) = 0.0_realkind  
    wi (en) = 0.0_realkind  
    x = h (en, na)  
    s = abs (x) + abs (zz)  
    p = x / s  
    q = zz / s  
    r = sqrt (p * p + q * q)  
    p = p / r  
    q = q / r  
    !     .......... row modification ..........
    do 290 j = na, n  
       zz = h (na, j)  
       h (na, j) = q * zz + p * h (en, j)  
       h (en, j) = q * h (en, j) - p * zz  
290 enddo
    !     .......... column modification ..........
    do 300 i = 1, en  
       zz = h (i, na)  
       h (i, na) = q * zz + p * h (i, en)  
       h (i, en) = q * h (i, en) - p * zz  
300 enddo
    !     .......... accumulate transformations ..........
    do 310 i = low, igh  
       zz = z (i, na)  
       z (i, na) = q * zz + p * z (i, en)  
       z (i, en) = q * z (i, en) - p * zz  
310 enddo
    !
    goto 330  
    !     .......... complex pair ..........
320 wr (na) = x + p  
    wr (en) = x + p  
    wi (na) = zz  
    wi (en) = - zz  
330 en = enm2  
    goto 60  
    !     .......... all roots found.  backsubstitute to find
    !                vectors of upper triangular form ..........
340 if (abs(norm)<1.e-14_realkind) goto 1001  
    !     .......... for en=n step -1 until 1 do -- ..........
    do  nn = 1, n  
       en = n + 1 - nn  
       p = wr (en)  
       q = wi (en)  
       na = en - 1  
       if(q<0.0_realkind)then
          goto 710
       elseif(q>0.0_realkind)then
          goto 800 
       else
          goto 600
       endif
       !     .......... real vector ..........
600    m = en  
       h (en, en) = 1.0_realkind  
       if (na==0) goto 800  
       !     .......... for i=en-1 step -1 until 1 do -- ..........
       do ii = 1, na  
          i = en - ii  
          w = h (i, i) - p  
          r = 0.0_realkind  
          !
          do 610 j = m, en  
             r = r + h (i, j) * h (j, en)  
610       enddo
          !
          if (wi (i) >=0.0_realkind) goto 630  
          zz = w  
          s = r  
          goto 700  
630       m = i  
          if (abs(wi(i))>1.e-14_realkind) goto 640  
          t = w  
          if (abs(t)>1.e-14_realkind) goto 635  
          tst1 = norm  
          t = tst1  
632       t = 0.01_realkind * t  
          tst2 = norm + t  
          if (tst2>tst1) goto 632  
635       h (i, en) = - r / t  
          goto 680  
          !     .......... solve real equations ..........
640       x = h (i, i + 1)  
          y = h (i + 1, i)  
          q = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i)  
          t = (x * s - zz * r) / q  
          h (i, en) = t  
          if (abs (x) <=abs (zz) ) goto 650  
          h (i + 1, en) = ( - r - w * t) / x  
          goto 680  
650       h (i + 1, en) = ( - s - y * t) / zz  
          !
          !     .......... overflow control ..........
680       t = abs (h (i, en) )  
          if (abs(t)<1.e-14_realkind) goto 700  
          tst1 = t  
          tst2 = tst1 + 1.0_realkind / tst1  
          if (tst2>tst1) goto 700  
          do 690 j = i, en  
             h (j, en) = h (j, en) / t  
690       enddo
          !
700       continue
       enddo
       !     .......... end real vector ..........
       goto 800  
       !     .......... complex vector ..........
710    m = na  
       !     .......... last vector component chosen imaginary so that
       !                eigenvector matrix is triangular ..........
       if (abs (h (en, na) ) <=abs (h (na, en) ) ) goto 720  
       h (na, na) = q / h (en, na)  
       h (na, en) = - (h (en, en) - p) / h (en, na)  
       goto 730  
720    call cdiv (0.0_realkind, - h (na, en), h (na, na) - p, q, h (na, na), &
            h (na, en) )
730    h (en, na) = 0.0_realkind  
       h (en, en) = 1.0_realkind  
       enm2 = na - 1  
       if (enm2==0) goto 800  
       !     .......... for i=en-2 step -1 until 1 do -- ..........
       do ii = 1, enm2  
          i = na - ii  
          w = h (i, i) - p  
          ra = 0.0_realkind  
          sa = 0.0_realkind  
          !
          do 760 j = m, en  
             ra = ra + h (i, j) * h (j, na)  
             sa = sa + h (i, j) * h (j, en)  
760       enddo
          !
          if (wi (i) >=0.0_realkind) goto 770  
          zz = w  
          r = ra  
          s = sa  
          goto 795  
770       m = i  
          if (abs(wi(i))> 1.e-14_realkind) goto 780  
          call cdiv ( - ra, - sa, w, q, h (i, na), h (i, en) )  
          goto 790  
          !     .......... solve complex equations ..........
780       x = h (i, i + 1)  
          y = h (i + 1, i)  
          vr = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i) - q * q  
          vi = (wr (i) - p) * 2.0_realkind * q  
          if (abs(vr)>1.e-14_realkind.or.abs(vi)>1.e-14_realkind) goto 784  
          tst1 = norm * (abs (w) + abs (q) + abs (x) + abs (y)  + abs (zz) )
          vr = tst1  
783       vr = 0.01_realkind * vr  
          tst2 = tst1 + vr  
          if (tst2>tst1) goto 783  
784       call cdiv (x * r - zz * ra + q * sa, x * s - zz * sa - q * &
               ra, vr, vi, h (i, na), h (i, en) )
          if (abs (x) <=abs (zz) + abs (q) ) goto 785  
          h (i + 1, na) = ( - ra - w * h (i, na) + q * h (i, en) ) &
               / x
          h (i + 1, en) = ( - sa - w * h (i, en) - q * h (i, na) ) &
               / x
          goto 790  
785       call cdiv ( - r - y * h (i, na), - s - y * h (i, en), &
               zz, q, h (i + 1, na), h (i + 1, en) )
          !
          !     .......... overflow control ..........
790       t = max (abs (h (i, na) ), abs (h (i, en) ) )  
          if (abs(t)<1.e-14_realkind) goto 795  
          tst1 = t  
          tst2 = tst1 + 1.0_realkind / tst1  
          if (tst2>tst1) goto 795  
          do 792 j = i, en  
             h (j, na) = h (j, na) / t  
             h (j, en) = h (j, en) / t  
792       enddo
          !
795       continue
       enddo
       !     .......... end complex vector ..........
800 continue
    enddo
    !     .......... end back substitution.
    !                vectors of isolated roots ..........
    do i = 1, n  
       if (i>=low.and.i<=igh) goto 840  
       !
       do 820 j = i, n  
          z (i, j) = h (i, j)  
820    enddo
840    continue
    enddo
    !     .......... multiply by transformation matrix to give
    !                vectors of original full matrix.
    !                for j=n step -1 until low do -- ..........
    do 881 jj = low, n  
       j = n + low - jj  
       m = min (j, igh)  
       !
       do 880 i = low, igh  
          zz = 0.0_realkind  
          !
          do 860 k = low, m  
             zz = zz + z (i, k) * h (k, j)  
860       enddo
          !
          z (i, j) = zz  
880    enddo
881 enddo
    !
    goto 1001  
    !     .......... set error -- all eigenvalues have not
    !                converged after 30*n iterations ..........
1000 ierr = en  
1001 return  
  end subroutine hqr2



  subroutine cdiv (ar, ai, br, bi, cr, ci)  
    implicit none  
    real(kind=realkind) :: ar, ai, br, bi, cr, ci  
    real(kind=realkind) :: s, ars, ais, brs, bis  
    s = abs (br) + abs (bi)  
    ars = ar / s  
    ais = ai / s  
    brs = br / s  
    bis = bi / s  
    s = brs**2.0_realkind + bis**2.0_realkind  
    cr = (ars * brs + ais * bis) / s  
    ci = (ais * brs - ars * bis) / s  
    return  
  end subroutine cdiv



  subroutine balbak (nm, n, low, igh, scale, m, z)  
    implicit none  
    integer :: i, j, k, m, n, ii, nm, igh, low  
    real(kind=realkind) :: scale (n), z (nm, m)  
    real(kind=realkind) :: s  
    !
    !     this subroutine is a translation of the algol procedure balbak,
    !     num. math. 13, 293-304(1969) by parlett and reinsch.
    !     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
    !
    !     this subroutine forms the eigenvectors of a real general
    !     matrix by back transforming those of the corresponding
    !     balanced matrix determined by  balanc.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrix.
    !
    !        low and igh are integers determined by  balanc.
    !
    !        scale contains information determining the permutations
    !          and scaling factors use d by  balanc.
    !
    !        m is the number of columns of z to be back transformed.
    !
    !        z contains the real and imaginary parts of the eigen-
    !          vectors to be back transformed in its first m columns.
    !
    !     on output
    !
    !        z contains the real and imaginary parts of the
    !          transformed eigenvectors in its first m columns.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    if (m==0) goto 200  
    if (igh==low) goto 120  
    !
    do 110 i = low, igh  
       s = scale (i)  
       !     .......... left hand eigenvectors are back transformed
       !                if the foregoing statement is replaced by
       !                s=1.0_realkind/scale(i). ..........
       do 100 j = 1, m  
          z (i, j) = z (i, j) * s  
100    enddo
       !
110 end do
    !     ......... for i=low-1 step -1 until 1,
    !               igh+1 step 1 until n do -- ..........
120 do  ii = 1, n  
       i = ii  
       if (i>=low.and.i<=igh) goto 140  
       if (i<low) i = low - ii  
       k = int(scale(i)) !this used to be scale(i) !marco!  
       if (k==i) goto 140  
       !
       do 130 j = 1, m  
          s = z (i, j)  
          z (i, j) = z (k, j)  
          z (k, j) = s  
130    end do
       !
140    continue
    end do
    !
    
200 return  
  end subroutine balbak

  subroutine eisrg1 (nm, ar, wr, wi, zr, ierr)  
    !     CERN EIGEN-VALUE PACKAGE (WITH MODIFICATIONS)
    !     EXTRACT FROM CERN LIBRARY (PETERLYNCH 19-OCT-1982)
    !
    !     COMPUTE ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS
    !     OF A GENERAL MATRIX
    implicit none  
    integer :: nm,  ierr,  low, igh 

    real(kind=realkind) :: AR(NM, NM), WR(NM), WI(NM), ZR(NM, NM), WORK(NM)
    integer::work2(nm)
    !scaling, general matrix
    CALL BALANC (NM, NM, AR, LOW, IGH, WORK)  
    !reduction, real matrix
    CALL ELMHES (NM, NM, LOW, IGH, AR, WORK2 )  
    !reduction, real matrix
    CALL ELTRAN (NM, NM, LOW, IGH, AR, WORK2, ZR)  
    !all eigenvalues and eig

    CALL HQR2 (NM, NM, LOW, IGH, AR, WR, WI, ZR, IERR)  
    IF (IERR/=0) RETURN  
    !back scaling, general m

    CALL BALBAK (NM, NM, LOW, IGH, WORK, NM, ZR)  
    RETURN  

  END subroutine EISRG1

  subroutine balanc (nm, n, a, low, igh, scale)  
    implicit none  
    integer :: i, j, k, l, m, n, jj, nm, igh, low, iexc  
    real(kind=realkind) :: a (nm, n), scale (n)  
    real(kind=realkind) :: c, f, g, r, s, b2, radix  
    logical :: noconv  
    !
    !     this subroutine is a translation of the algol procedure balance,
    !     num. math. 13, 293-304(1969) by parlett and reinsch.
    !     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
    !
    !     this subroutine balances a real matrix and isolates
    !     eigenvalues whenever possible.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrix.
    !
    !        a contains the input matrix to be balanced.
    !
    !     on output
    !
    !        a contains the balanced matrix.
    !
    !        low and igh are two integers such that a(i,j)
    !          is equal to zero if
    !           (1) i is greater than j and
    !           (2) j=1,...,low-1 or i=igh+1,...,n.
    !
    !        scale contains information determining the
    !           permutations and scaling factors use d.
    !
    !     suppose that the principal submatrix in rows low through igh
    !     has been balanced, that p(j) denotes the index interchanged
    !     with j during the permutation step, and that the elements
    !     of the diagonal matrix use d are denoted by d(i,j).  then
    !        scale(j) = p(j),    for j = 1,...,low-1
    !                 = d(j,j),      j = low,...,igh
    !                 = p(j)         j = igh+1,...,n.
    !     the order in which the interchanges are made is n to igh+1,
    !     then 1 to low-1.
    !
    !     note that 1 is returned for igh if igh is zero formally.
    !
    !     the algol procedure exc contained in balance appears in
    !     balanc  in line.  (note that the algol roles of identifiers
    !     k,l have been reversed.)
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------
    !
    radix = 16.0_realkind  
    !
    b2 = radix * radix  
    k = 1  
    l = n  
    goto 100  
    !     .......... in-line procedure for row and
    !                column exchange ..........
20  scale(m) = real(j,realkind)  
    if (j==m) goto 50  
    !
    do 30 i = 1, l  
       f = a (i, j)  
       a (i, j) = a (i, m)  
       a (i, m) = f  
30  end do
    !
    do 40 i = k, n  
       f = a (j, i)  
       a (j, i) = a (m, i)  
       a (m, i) = f  
40  end do
    !
50  goto (80, 130), iexc  
    !     .......... search for rows isolating an eigenvalue
    !                and push them down ..........
80  if (l==1) goto 280  
    l = l - 1  
    !     .......... for j=l step -1 until 1 do -- ..........
100 do  jj = 1, l  
       j = l + 1 - jj  
       !
       do  i = 1, l  
          if (i==j) goto 110  
          if (abs(a(j, i))> 1.e-14_realkind) goto 120  
110       continue 
       end do
       !
       m = l  
       iexc = 1  
       goto 20 
120    continue
    end do
    !
    goto 140  
    !     .......... search for columns isolating an eigenvalue
    !                and push them left ..........
130 k = k + 1  
    !
140 do j = k, l  
       !
       do  i = k, l  
          if (i==j) goto 150  
          if (abs(a(i,j))> 1.e-12_realkind) goto 170  
150       continue
       end do
       !
       m = k  
       iexc = 2  
       goto 20  
170    continue
    end do
    !     .......... now balance the submatrix in rows k to l ..........
    do 180 i = k, l  
       scale (i) = 1.0_realkind  
180 enddo
    !     .......... iterative loop for norm reduction ..........
190 noconv = .false.  
    !
    do i = k, l  
       c = 0.0_realkind  
       r = 0.0_realkind  
       !
       do j = k, l  
          if (j==i) goto 200  
          c = c + abs (a (j, i) )  
          r = r + abs (a (i, j) )  
200       continue
       end do
       !     .......... guard against zero c or r due to underflow ..........
       if (abs(c)<1.e-14_realkind.or.abs(r)<1.e-14_realkind) goto 270  
       g = r / radix  
       f = 1.0_realkind  
       s = c + r  
210    if (c>=g) goto 220  
       f = f * radix  
       c = c * b2  
       goto 210  
220    g = r * radix  
230    if (c<g) goto 240  
       f = f / radix  
       c = c / b2  
       goto 230  
       !     .......... now balance ..........
240    if ( (c + r) / f>=0.95_realkind * s) goto 270  
       g = 1.0_realkind / f  
       scale (i) = scale (i) * f  
       noconv = .true.  
       !
       do 250 j = k, n  
          a (i, j) = a (i, j) * g  
250    enddo
       !
       do 260 j = 1, l  
          a (j, i) = a (j, i) * f  
260    enddo
270    continue       !
    end do
    !
    if (noconv) goto 190  
    !
280 low = k  
    igh = l  
    return  
  end subroutine balanc



  SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K,KALPHA, A, LDA, B, LDB, &
       BETA, C, LDC )
    !    .. Scalar Arguments ..
    CHARACTER*1        TRANSA, TRANSB
    INTEGER            M, N, K, LDA, LDB, LDC
    REAL(KIND=REALKIND)               KALPHA, BETA
    !    .. Array Arguments ..
    REAL(KIND=REALKIND)               A( LDA, * ), B( LDB, * ), C( LDC, * )

    ! Purpose


    ! SGEMM  performs one of the matrix-matrix operations

    !    C := alpha*op( A )*op( B ) + beta*C,

    ! where  op( X ) is one of

    !    op( X ) = X   or   op( X ) = X',

    ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
    ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

    ! Parameters


    ! TRANSA - CHARACTER*1.
    !          On entry, TRANSA specifies the form of op( A ) to be used in
    !          the matrix multiplication as follows:

    !             TRANSA = 'N' or 'n',  op( A ) = A.

    !             TRANSA = 'T' or 't',  op( A ) = A'.

    !             TRANSA = 'C' or 'c',  op( A ) = A'.

    !          Unchanged on exit.

    ! TRANSB - CHARACTER*1.
    !          On entry, TRANSB specifies the form of op( B ) to be used in
    !          the matrix multiplication as follows:

    !             TRANSB = 'N' or 'n',  op( B ) = B.

    !             TRANSB = 'T' or 't',  op( B ) = B'.

    !             TRANSB = 'C' or 'c',  op( B ) = B'.

    !          Unchanged on exit.

    ! M      - INTEGER.
    !          On entry,  M  specifies  the number  of rows  of the  matrix
    !          op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !          Unchanged on exit.

    ! N      - INTEGER.
    !          On entry,  N  specifies the number  of columns of the matrix
    !          op( B ) and the number of columns of the matrix C. N must be
    !          at least zero.
    !          Unchanged on exit.

    ! K      - INTEGER.
    !          On entry,  K  specifies  the number of columns of the matrix
    !          op( A ) and the number of rows of the matrix op( B ). K must
    !          be at least  zero.
    !          Unchanged on exit.

    ! KALPHA - REAL            .
    !          On entry, KALPHA specifies the scalar alpha.
    !          Unchanged on exit.

    ! A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
    !          k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !          Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !          part of the array  A  must contain the matrix  A,  otherwise
    !          the leading  k by m  part of the array  A  must contain  the
    !          matrix A.
    !          Unchanged on exit.

    ! LDA    - INTEGER.
    !          On entry, LDA specifies the first dimension of A as declared
    !          in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !          LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !          least  max( 1, k ).
    !          Unchanged on exit.

    ! B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
    !          n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !          Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !          part of the array  B  must contain the matrix  B,  otherwise
    !          the leading  n by k  part of the array  B  must contain  the
    !          matrix B.
    !          Unchanged on exit.

    ! LDB    - INTEGER.
    !          On entry, LDB specifies the first dimension of B as declared
    !          in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !          LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !          least  max( 1, n ).
    !          Unchanged on exit.

    ! BETA   - REAL            .
    !          On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !          supplied as zero then C need not be set on input.
    !          Unchanged on exit.

    ! C      - REAL             array of DIMENSION ( LDC, n ).
    !          Before entry, the leading  m by n  part of the array  C must
    !          contain the matrix  C,  except when  beta  is zero, in which
    !          case C need not be set on entry.
    !          On exit, the array  C  is overwritten by the  m by n  matrix
    !          ( alpha*op( A )*op( B ) + beta*C ).

    ! LDC    - INTEGER.
    !          On entry, LDC specifies the first dimension of C as declared
    !          in  the  calling  (sub)  program.   LDC  must  be  at  least
    !          max( 1, m ).
    !          Unchanged on exit.


    ! Level 3 Blas routine.

    ! -- Written on 8-February-1989.
    !    Jack Dongarra, Argonne National Laboratory.
    !    Iain Duff, AERE Harwell.
    !    Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !    Sven Hammarling, Numerical Algorithms Group Ltd.


    INTRINSIC          MAX
    !    .. Local Scalars ..
    LOGICAL            NOTA, NOTB
    INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
    REAL(KIND=REALKIND)               TEMP
    !    .. Parameters ..
    REAL(KIND=REALKIND)               ONE         , ZERO
    PARAMETER        ( ONE = 1.0E+0_realkind, ZERO = 0.0E+0_realkind )
    !    ..
    !    .. Executable Statements ..

    !    Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    !    transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    !    and  columns of  A  and the  number of  rows  of  B  respectively.

    NOTA  = LSAME( TRANSA, 'N' )
    NOTB  = LSAME( TRANSB, 'N' )
    IF( NOTA )THEN
       NROWA = M
       NCOLA = K
    ELSE
       NROWA = K
       NCOLA = M
    END IF
    IF( NOTB )THEN
       NROWB = K
    ELSE
       NROWB = N
    END IF

    !    Test the input parameters.

    INFO = 0
    IF(      ( .NOT.NOTA                 ).AND. &
         ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
       INFO = 1
    ELSE IF( ( .NOT.NOTB                 ).AND. &
         ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
       INFO = 2
    ELSE IF( M  <0               )THEN
       INFO = 3
    ELSE IF( N  <0               )THEN
       INFO = 4
    ELSE IF( K  <0               )THEN
       INFO = 5
    ELSE IF( LDA<MAX( 1, NROWA ) )THEN
       INFO = 8
    ELSE IF( LDB<MAX( 1, NROWB ) )THEN
       INFO = 10
    ELSE IF( LDC<MAX( 1, M     ) )THEN
       INFO = 13
    END IF
    IF( INFO/=0 )THEN
       CALL XERBLA( 'SGEMM ', INFO )
       RETURN
    END IF

    !    Quick return if possible.

    IF( ( M==0 ).OR.( N==0 ).OR. &
         ( ( ( abs(KALPHA-ZERO)<1.e-15_realkind ).OR.( K==0 ) ).AND.&
         ( abs(BETA-ONE)<1.e-15_realkind ) ) ) &
         RETURN

    !    And if  alpha==zero.

    IF( abs(KALPHA-ZERO)<1.e-15_realkind  )THEN
       IF( abs(BETA-ZERO)<1.e-15_realkind )THEN
          DO 20, J = 1, N
             DO 10, I = 1, M
                C( I, J ) = ZERO
10           ENDDO
20        ENDDO
       ELSE
          DO 40, J = 1, N
             DO 30, I = 1, M
                C( I, J ) = BETA*C( I, J )
30           ENDDO
40        ENDDO
       END IF
       RETURN
    END IF

    !    Start the operations.

    IF( NOTB )THEN
       IF( NOTA )THEN

          !          Form  C := alpha*A*B + beta*C.

          DO 90, J = 1, N
             IF( abs(BETA-ZERO)<1.e-15_realkind )THEN
                DO 50, I = 1, M
                   C( I, J ) = ZERO
50              ENDDO
             ELSE IF( abs(BETA-ONE)>1.e-15_realkind )THEN
                DO 60, I = 1, M
                   C( I, J ) = BETA*C( I, J )
60              ENDDO
             END IF
             DO 80, L = 1, K
                IF( abs(B( L, J )-ZERO)>1.e-15_realkind )THEN
                   TEMP = KALPHA*B( L, J )
                   DO 70, I = 1, M
                      C( I, J ) = C( I, J ) + TEMP*A( I, L )
70                 ENDDO
                END IF
80           ENDDO
90        ENDDO
       ELSE

          !          Form  C := alpha*A'*B + beta*C

          DO 120, J = 1, N
             DO 110, I = 1, M
                TEMP = ZERO
                DO 100, L = 1, K
                   TEMP = TEMP + A( L, I )*B( L, J )
100             ENDDO
                IF( abs(BETA-ZERO)<1.e-15_realkind )THEN
                   C( I, J ) = KALPHA*TEMP
                ELSE
                   C( I, J ) = KALPHA*TEMP + BETA*C( I, J )
                END IF
110          ENDDO
120       ENDDO
       END IF
    ELSE
       IF( NOTA )THEN

          !          Form  C := alpha*A*B' + beta*C

          DO 170, J = 1, N
             IF( abs(BETA-ZERO)<1.e-15_realkind )THEN
                DO 130, I = 1, M
                   C( I, J ) = ZERO
130             ENDDO
             ELSE IF( abs(BETA-ONE)>1.e-15_realkind )THEN
                DO 140, I = 1, M
                   C( I, J ) = BETA*C( I, J )
140             ENDDO
             END IF
             DO 160, L = 1, K
                IF( abs(B( J, L )-ZERO)>1.e-15_realkind )THEN
                   TEMP = KALPHA*B( J, L )
                   DO 150, I = 1, M
                      C( I, J ) = C( I, J ) + TEMP*A( I, L )
150                ENDDO
                END IF
160          ENDDO
170       ENDDO
       ELSE

          !          Form  C := alpha*A'*B' + beta*C

          DO 200, J = 1, N
             DO 190, I = 1, M
                TEMP = ZERO
                DO 180, L = 1, K
                   TEMP = TEMP + A( L, I )*B( J, L )
180             ENDDO
                IF( abs(BETA-ZERO)<1.e-15_realkind )THEN
                   C( I, J ) = KALPHA*TEMP
                ELSE
                   C( I, J ) = KALPHA*TEMP + BETA*C( I, J )
                END IF
190          ENDDO
200       ENDDO
       END IF
    END IF

    RETURN

    !    End of SGEMM .

  end SUBROUTINE SGEMM


  logical function lsame( ca, cb )
    implicit none
    character          ca, cb
    intrinsic          ichar
    integer            inta, intb, zcode

    lsame = ca==cb
    if( lsame )  return

    zcode=ichar( 'Z' )

    inta = ichar( ca )
    intb = ichar( cb )


    if( zcode==90 .or. zcode==122 ) then

       !        ascii is assumed - zcode is the ascii code of either lower or
       !        upper case 'z'.

       if( inta>=97 .and. inta<=122 ) inta = inta - 32
       if( intb>=97 .and. intb<=122 ) intb = intb - 32

    else if( zcode==233 .or. zcode==169 ) then

       !        ebcdic is assumed - zcode is the ebcdic code of either lower or
       !        upper case 'z'.

       if( inta>=129 .and. inta<=137 .or.      &
            inta>=145 .and. inta<=153 .or.   &
            inta>=162 .and. inta<=169 ) inta = inta + 64
       if( intb>=129 .and. intb<=137 .or.      &
            intb>=145 .and. intb<=153 .or.   &
            intb>=162 .and. intb<=169 ) intb = intb + 64

    else if( zcode==218 .or. zcode==250 ) then

       !        ascii is assumed, on prime machines - zcode is the ascii code
       !        plus 128 of either lower or upper case 'z'.

       if( inta>=225 .and. inta<=250 ) inta = inta - 32
       if( intb>=225 .and. intb<=250 ) intb = intb - 32
    end if
    lsame = inta==intb

  end function lsame



  subroutine xerbla( srname, info )
    implicit none
    character(len=6)::srname
    integer:: info
    write( *, fmt = 9999 )srname, info
9999 format( ' ** on entry to ', a6, ' parameter number ', i2, ' had ', &
         'an illegal value' )
  end subroutine xerbla







  subroutine hhsolvitr (klon, klat,  klev , kpbpts &
       , ptwodt , pdiv , nlsimpx  , epsg  , itrtyp, pretyp  , abserr, relerr )
    !                                                                       
    ! HHSOLVITR - MASTER ROUTINE FOR SOLUTION OF HELMHOLTZ-EQUATIONS     
    !                   IN THE SEMI-IMPLICIT TIMESTEPPING SCHEME               
    !                                                                       
    !     J.E. HAUGEN            HIRLAM
    !
    !     PARALLEL ITERATIVE VERSION:
    !
    !     D. BJOERGE             DNMI           1994-12-20
    !     D. BJOERGE             DNMI           1995-04-22, UPDATE
    !     D. BJOERGE             DNMI           1999-04-28, HIRLAM
    !                                                                       
    !                                                                       
    !     PURPOSE.                                                          
    !     --------                                                          
    !                                                                       
    !     THIS ROUTINE SOLVES THE THREE-DIMENSIONAL HELMHOLTZ-EQUATION      
    !     FOR THE DIVERGENCE BY VERTICAL DECOPLING INTO KLEV HORIZONTAL     
    !     HELMHOLTZ-EQUATIONS. THE EQUATIONS ARE SOLVED BY A CHOICE OF
    !     ITERATIVE METHODS.                                            
    !                                                                       
    !   INTERFACE.                                                        
    !     ----------                                                        
    !                                                                       
    !     *HHSOLVITR* IS CALLED FROM THE MAIN TIME STEPPING ROUTINES           
    !                                                                       
    !     INPUT PARAMETERS.                                                 
    !     -----------------                                                 
    !                                                                       
    !     *KLON*      NUMBER OF GRIDPOINTS IN X-DIRECTION                   
    !     *KLAT*      NUMBER OF GRIDPOINTS IN Y-DIRECTION                   
    !     *KLEV*      NUMBER OF VERTICAL LEVELS                             
    !     *KPBPTS*    NUMBER OF PASSIVE BOUNDARYLINES                       
    !     *PTWODT*    CURRENT DOUBLE TIMESTEP                               
    !     *FFMEAN*    MEAN-VALUE OF CORIOLIS-PARAMETER                      
    !     *PDIV*      RIGHT HAND SIDE OF THE HH-EQUATIONS                   
    !                                                                       
    !     OUTPUT PARAMETERS:                                                
    !     ------------------                                                
    !                                                                       
    !     *PDIV*     SOLUTION OF HH-EQUATIONS                              
    !                                                                       
    !     EXTERNALS.                                                        
    !     ----------                                                        
    !                                                                       
    !     *SGEMM*     MATRIX MULTIPLICATION (MATRIX OUTER PRODUCT)
    !     *HHPRE*     PREPARATIONS FOR ITERATIVE SOLVER
    !     *CGSPRE*    ITERATIVE HELMHOLTZ SOLVERS
    !     *CGPRE*
    !     *RICHPRE*
    !
    !---------------------------------------------------------------   

    IMPLICIT NONE

    !    DECLARATION OF GLOBAL PARAMETERS

    INTEGER KLON,KLAT,KLEV,KPBPTS
    REAL(KIND=REALKIND)    PTWODT,PDIV(KLON,KLAT,KLEV)
    !
    REAL(KIND=REALKIND)  EPSG 
    LOGICAL NLSIMPX

    INTEGER ITRTYP, PRETYP(KLEV)
    REAL(KIND=REALKIND) ABSERR, RELERR
    !
    ! LOCAL VARIABLES
    !
    REAL(KIND=REALKIND) ZDIV(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) WORK1(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) WORK2(KLON,KLAT,KLEV)                     
    REAL(KIND=REALKIND) AMAT(KLON,KLAT,4,KLEV)
    INTEGER IERR, INFO
    INTEGER ISTART,ISTOP,JSTART,JSTOP
    INTEGER I, J

    !db      SAVE AMAT
    ! DEFINE THE MATRIX TO BE USED IN 'AMULT' AND 'PRECOND'

    !db      IF (NSTEP <= 1) THEN
    CALL HHPRE(AMAT,PTWODT,KLON,KLAT,KLEV, ISTART,ISTOP,JSTART,JSTOP,epsg)
    !db      ENDIF

    !    VERTICAL DECOUPLING
    IF( KLEV>1 ) THEN
       CALL SGEMM(              &
            'N','T', &
            KLON*KLAT,KLEV,KLEV, &
            1.0_realkind,PDIV, &
            KLON*KLAT, &
            EIGINV, &
            KLEV, &
            0.0_realkind,ZDIV, &
            KLON*KLAT)
    ELSE
       DO J=1,KLAT
          DO I=1,KLON
             ZDIV(I,J,1) = PDIV(I,J,1)
          ENDDO
       ENDDO
    ENDIF

    ! THE COUPLED MATRIX PDIV MAY NOW BE USED AS A WORK-ARRAY

    ! CHECK THE CHOICE OF ITERATIVE SOLVER

    IF (ITRTYP <1 .OR. ITRTYP >3) THEN
       IF (MYPE == 0) THEN
          WRITE(*,*) ' YOU HAVE REQUESTED WRONG ITERATION TYPE = ',ITRTYP
          WRITE(*,*) ' ITRTYP IS SET TO 2, CONJUGATE GRADIENT SQUARED'
       ENDIF
       ITRTYP = 2
    ENDIF

    ! CALL THE ITERATIVE SOLVER

    IF (ITRTYP == 1) THEN
       CALL CGPRE(KLON,KLAT,KLEV,KLEV,      &
            AMAT, ZDIV, &
            PDIV, WORK1, WORK2, &
            ISTART,ISTOP,JSTART,JSTOP, &
            PRETYP, ABSERR, RELERR,  IERR)
    ELSEIF (ITRTYP == 2) THEN
       CALL CGSPRE(KLON,KLAT,KLEV,KLEV,  &
            AMAT, ZDIV, &
            PDIV, WORK1, WORK2, &
            ISTART,ISTOP,JSTART,JSTOP, &
            PRETYP, ABSERR, RELERR, IERR)
    ELSEIF (ITRTYP == 3) THEN
       CALL RICHPRE(KLON,KLAT,KLEV,KLEV, &
            AMAT, ZDIV, &
            PDIV, WORK1, WORK2, &
            ISTART,ISTOP,JSTART,JSTOP, &
            PRETYP, ABSERR, RELERR, IERR)
    ENDIF

    IF (IERR < 0) THEN
       CALL STOP_PROGRAM( 'HHSOLVITR DID NOT CONVERGE')
    ENDIF

    !    VERTICAL COUPLING
    IF( KLEV>1 ) THEN

       CALL SGEMM( 'N','T', KLON*KLAT,KLEV,KLEV, 1.0_realkind,ZDIV, &
            KLON*KLAT, EIG, KLEV, 0.0_realkind,PDIV,   KLON*KLAT)
    ELSE
       DO J=1,KLAT
          DO I=1,KLON
             PDIV(I,J,1) = ZDIV(I,J,1)
          ENDDO
       ENDDO
    ENDIF

    RETURN                                                            
  END subroutine hhsolvitr

  SUBROUTINE HHPRE( AMAT  , PTWODT,  KLON,KLAT,KLEV, &
       ISTART, ISTOP , JSTART, JSTOP,epsg)
    !
    !    DEFINE MATRIX FOR ITERATIVE SOLVER AND PRECONDITIONER
    !
    !
    !     D. BJOERGE      DNMI    1994-12-20
    !     D. BJOERGE      DNMI    1997-10-23, polar-stereographic
    !     D. BJOERGE      DNMI    1999-04-28, HIRLAM
    !
    !
    ! called by:
    ! ----------
    !
    !         hhsolv
    !         delnmi
    !
    !
    ! input parameters:
    ! -----------------
    !
    !     ptwodt   double time step length in seconds
    !     ffmean   mean value of Coriolis parameter
    !
    !
    ! output parameters:
    ! ------------------
    !
    !     amat   matrix defining the left hand side of the helmholtz
    !            equation.
    !
    ! This is a new version where the original matrix 'amat' 
    ! is scaled with the original diagonal elements. 
    ! The diagonal element is redefined to be equal to the scaling
    ! factor, and is used to scale the RHS in the iterative solver.
    !
    !
    !                   AMAT(JX,JY,3)
    !                    |
    !  AMAT(JX,JY,2) -- AMAT(JX,JY,1) -- AMAT(JX,JY,2)
    !                    |
    !                   AMAT(JX,JY,4)
    !

    use comhkp

    use boundaryRelaxation
    IMPLICIT NONE
    ! GLOBAL VARIABLES

    INTEGER KLON,KLAT,KLEV
    REAL(KIND=REALKIND) AMAT(KLON,KLAT,4,klev)
    REAL(KIND=REALKIND) PTWODT
    real(kind=realkind),intent(in)::epsg

    INTEGER ISTART,ISTOP,JSTART,JSTOP

    ! LOCAL VARIABLES

    INTEGER JX, JY, JK, JJ, JYG
    REAL(KIND=REALKIND) RDLAM2, RDTH2, CCONST, SCALE
    real(kind=realkind) a1, a2, a3, a4
    REAL(KIND=REALKIND) ZDT, ZFF, Z1PEPS, ZAA

    do jk=1,klev
       do jj=1,4
          do jy=1,klat
             do jx=1,klon
                AMAT(jx,jy,jj,jk) = 1._realkind
             enddo
          enddo
       enddo
    enddo

    if(atleft)then
       istart = NPBPTS + 2
    else
       istart = 2
    endif

    if(atright)then
       istop = klon - 2 - NPBPTS
    else
       istop = klon - 1
    endif

    if(atbase)then
       jstart = NPBPTS + 2
    else
       jstart = 2
    endif

    if(attop)then
       jstop = klat - 2 - NPBPTS
    else
       jstop = klat - 1
    endif

    DO JK=1,klev
       ! for the semi-implicit scheme
       ZDT = 0.5_realkind*PTWODT                                                  
       ZFF = FFMEAN                                                      

       Z1PEPS = 1.0_realkind+EPSG
       ZDT    = Z1PEPS*ZDT

       ZAA = ZDT*ZFF                                                     
       ZAA = 1.0_realkind+ZAA*ZAA 

       CCONST = ( ZDT*ZDT * EDEPTH(JK) )/ZAA

       CCONST = CCONST * RA * RA
       RDLAM2 = RDLAM* RDLAM * CCONST
       RDTH2  = RDTH * RDTH  * CCONST

       DO JY=JSTART,JSTOP
          JYG = JDATASTART - 1 + JY
          a2 = -RDLAM2 / (cosu(JYG)  * cosu(JYG))
          a3 = -RDTH2  * (cosv(JYG)  / cosu(JYG))
          a4 = -RDTH2  * (cosv(JYG-1)/ cosu(JYG))
          a1 = 1._realkind - 2._realkind*a2 - a3 - a4

          SCALE = 1._realkind / a1
          DO JX=ISTART,ISTOP
             AMAT(JX,JY,1,JK) = SCALE
             AMAT(JX,JY,2,JK) = SCALE * a2
             AMAT(JX,JY,3,JK) = SCALE * a3
             AMAT(JX,JY,4,JK) = SCALE * a4
          ENDDO
       ENDDO

    ENDDO

    RETURN
  END SUBROUTINE HHPRE

  SUBROUTINE CGPRE(KLON,KLAT,KLEV,KKLEV,           &
       AMATP,RHS, &
       P, Q, R, &
       ISTART,ISTOP,JSTART,JSTOP, &
       PRETYP, ABSERR, RELERR, &
       IERR)
    !
    ! ITERATIVE HELMHOLTZ SOLVER
    ! PRECONDITIONED CONJUGATE GRADIENT METHOD
    !
    !
    !     D. BJOERGE      DNMI     1994-12-22, original version
    !                              1999-06-23, modified for HIRLAM MPP
    !                                                                       
    !
    ! CALLED BY:
    ! ----------
    !
    !     HHSOLV
    !     DELNMI
    !
    !
    ! INPUT PARAMETERS:
    ! -----------------
    !
    !     KKLEV   NUMBER OF EQUATIONS TO SOLVE
    !     AMATP   MATRIX OF THE LEFT HAND SIDE
    !     RHS     RIGHT HAND SIDE OF THE HELMHOLTZ EQUATION
    !
    !     P,Q,R   WORK ARRAYS
    !     ABSERR   ABSOLUTE ERROR LIMIT FOR THE ITERATION
    !     RELERR   RELATIVE ERROR LIMIT FOR THE ITERATION
    !
    !
    ! OUTPUT PARAMETERS:
    ! ------------------
    !
    !     IERR
    !          = -1    DID NOT CONVERGE
    !          >= 0    NUMBER OF ITERATION NEEDED
    !     RHS  THE SOLUTION
    !
    !
    !
    ! EXTERNALS
    ! ---------
    !
    !          AMULT   -- LEFT HAND SIDE OF HELMHOLTZ EQUATION
    !          HHGV    -- GLOBAL OPERATIONS
    !          PRECOND -- PRECONDITIONER
    !
    !-------------------------------------------------------------------
    !

    IMPLICIT NONE
    ! GLOBAL VARIABLES

    INTEGER KLON, KLAT, KLEV
    REAL(KIND=REALKIND)    AMATP(KLON,KLAT,4,KLEV)
    INTEGER KKLEV
    REAL(KIND=REALKIND)    RHS(KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    P  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    Q  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    R  (KLON,KLAT,KKLEV)
    INTEGER ISTART,ISTOP,JSTART,JSTOP
    INTEGER PRETYP(KLEV)
    REAL(KIND=REALKIND) ABSERR, RELERR
    INTEGER IERR

    ! LOCAL VARIABLES

    REAL(KIND=REALKIND) SDOT(KLEV), ROM1(KLEV), ROM2(KLEV)
    REAL(KIND=REALKIND) ALPHA, BETA
    REAL(KIND=REALKIND) ZRR0(KLEV), ZRR(KLEV), TYPICAL(KLEV), RMS(KLEV)
    REAL(KIND=REALKIND) SCALE, DUMMY(KLON,KLAT,KKLEV)
    INTEGER ITER, LAT, LON, LEV, KLEVE

    ! SET HOW OFTEN TO CHECK FOR CONVERGENCE

    INTEGER    TSTITR
    PARAMETER (TSTITR = 3)

    ! SET SECURITY TRESHOLD FOR ITERATION

    INTEGER    MAXITR
    PARAMETER (MAXITR = 100)

    KLEVE = KKLEV

    ! SET BOUNDARIES TO ZERO

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          RHS(ISTART-1,LAT,LEV) = 0.0_realkind
          RHS(ISTOP +1,LAT,LEV) = 0.0_realkind
          P  (ISTART-1,LAT,LEV) = 0.0_realkind
          P  (ISTOP +1,LAT,LEV) = 0.0_realkind
       ENDDO
       DO LON = ISTART,ISTOP
          RHS(LON,JSTART-1,LEV) = 0.0_realkind
          RHS(LON,JSTOP +1,LEV) = 0.0_realkind
          P  (LON,JSTART-1,LEV) = 0.0_realkind
          P  (LON,JSTOP +1,LEV) = 0.0_realkind
       ENDDO
    ENDDO

    ! INITIALIZE ITERATIONS.
    ! SET FIRST GUESS EQUAL TO THE SOLUTION OF LAST TIME STEP.
    ! ALSO THE RIGHT HAND SIDE MUST BE SCALED IN THE SAME
    ! WAY AS THE AMATP-MATRIX WAS SCALED IN HHPRE.
    !
    DO LEV=1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             SCALE = AMATP(LON,LAT,1,LEV)
             RHS(LON,LAT,LEV) =  RHS(LON,LAT,LEV) * SCALE
          ENDDO
       ENDDO
    ENDDO


    ! START ITERATION 

    ! FIRST ITERATION IS HARDCODED TO AVOID IF-TESTS

    ITER = 1

    CALL AMULT(KLON,KLAT,KLEV,  &
         AMATP, RHS,R, KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP)

    DO LEV=1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             R(LON,LAT,LEV) =  RHS(LON,LAT,LEV) - R(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO


    ! CHECK FOR CONVERGENCE

    CALL HHGV(KLON,KLAT,KLEV, &
         ISTART,ISTOP,JSTART,JSTOP, &
         2,DUMMY,R,RHS,KLEVE,DUMMY,RMS,TYPICAL)

    DO LEV=1,KLEVE
       IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
          KLEVE = LEV-1
          GOTO 10
       ENDIF
    ENDDO

10  CONTINUE
    IF (KLEVE == 0) THEN
       IERR = 0
       RETURN
    ENDIF

    ! END OF CHECK FOR CONVERGENCE

    CALL PRECOND(KLON,KLAT,KLEV, &
         AMATP,R,Q,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP,PRETYP)

    CALL HHGV(KLON,KLAT,KLEV, &
         ISTART,ISTOP,JSTART,JSTOP, &
         1,R,Q,DUMMY,KLEVE,ROM1,DUMMY,DUMMY)

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             P(LON,LAT,LEV) = Q(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    CALL AMULT(KLON,KLAT,KLEV, &
         AMATP,P,Q,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP)

    CALL HHGV(KLON,KLAT,KLEV, &
         ISTART,ISTOP,JSTART,JSTOP, &
         1,P,Q,DUMMY,KLEVE,SDOT,DUMMY,DUMMY)     

    DO LEV=1,KLEVE
       ALPHA     = ROM1(LEV)/SDOT(LEV)
       ROM2(LEV) = ROM1(LEV)
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             RHS(LON,LAT,LEV)=RHS(LON,LAT,LEV) +  ALPHA * P(LON,LAT,LEV)
             R(LON,LAT,LEV)=    R(LON,LAT,LEV) -  ALPHA * Q(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    DO 100 ITER = 2, MAXITR

       CALL PRECOND(KLON,KLAT,KLEV, AMATP,R,Q,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP,PRETYP)


       ! POSSIBLY CHECK FOR STOPPING THE ITERATION.

       IF (MOD(ITER,TSTITR) == 0) THEN

          CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
               3,R,Q,RHS,KLEVE,ROM1,RMS,TYPICAL)
          ! CHECK FOR CONVERGENCE

          DO LEV=1,KLEVE
             IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
                KLEVE = LEV-1
                GOTO 20
             ENDIF
          ENDDO

20        CONTINUE
          IF (KLEVE == 0) THEN
             IERR = ITER -1
             RETURN
          ENDIF

          ! END OF CHECK FOR CONVERGENCE

       ELSE
          CALL HHGV(KLON,KLAT,KLEV,  ISTART,ISTOP,JSTART,JSTOP, &
               1,R,Q,DUMMY,KLEVE,ROM1,DUMMY,DUMMY)
       ENDIF

       DO LEV = 1,KLEVE
          BETA = ROM1(LEV)/ROM2(LEV)
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                P(LON,LAT,LEV) = Q(LON,LAT,LEV) + BETA * P(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

       CALL AMULT(KLON,KLAT,KLEV,AMATP,P,Q,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP)

       CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
            1,P,Q,DUMMY,KLEVE,SDOT,DUMMY,DUMMY)

       DO LEV=1,KLEVE
          ALPHA     = ROM1(LEV)/SDOT(LEV)
          ROM2(LEV) = ROM1(LEV)
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                RHS(LON,LAT,LEV)=RHS(LON,LAT,LEV) + ALPHA* P(LON,LAT,LEV)
                R  (LON,LAT,LEV)=R  (LON,LAT,LEV) - ALPHA* Q(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

100 enddo

    IERR = -1
    RETURN
  END SUBROUTINE CGPRE

  SUBROUTINE CGSPRE(KLON,KLAT,KLEV,KKLEV, &
       AMATP,RHS, &
       P, Q, R, &
       ISTART,ISTOP,JSTART,JSTOP, &
       PRETYP, ABSERR, RELERR,   IERR)
    !
    ! ITERATIVE HELMHOLTZ SOLVER
    ! PRECODITIONED CONUGATE SQUARED METHOD
    !
    !
    !     D. BJOERGE            DNMI           1994-12-22
    !
    !
    ! CALLED BY:
    ! ----------
    !
    !     HHSOLV
    !     DELNMI
    !
    !
    ! INPUT PARAMETERS:
    ! -----------------
    !
    !     KKLEV   NUMBER OF EQUATIONS TO SOLVE
    !     AMATP   MATRIX OF THE LEFT HAND SIDE
    !     RHS     RIGHT HAND SIDE OF THE HELMHOLTZ EQUATION
    !
    !     P,Q,R   WORK ARRAYS
    !     ABSERR   ABSOLUTE ERROR LIMIT FOR THE ITERATION
    !     RELERR   RELATIVE ERROR LIMIT FOR THE ITERATION
    !
    ! OUTPUT PARAMETERS:
    ! ------------------
    !
    !     IERR
    !          = -1    DID NOT CONVERGE
    !          >= 0    NUMBER OF ITERATION NEEDED
    !     RHS  THE SOLUTION
    !
    !
    !
    ! EXTERNALS
    ! ---------
    !
    !          AMULT   -- LEFT HAND SIDE OF HELMHOLTZ EQUATION
    !          HHGV    -- GLOBAL OPERATIONS
    !          PRECOND -- PRECONDITIONER
    !
    !------------------------------------------------------------------
    !

    IMPLICIT NONE
    ! GLOBAL VARIABLES

    INTEGER KLON, KLAT, KLEV
    REAL(KIND=REALKIND)    AMATP(KLON,KLAT,4,KLEV)
    INTEGER KKLEV
    REAL(KIND=REALKIND)    RHS(KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    P  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    Q  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    R  (KLON,KLAT,KKLEV)
    INTEGER ISTART,ISTOP,JSTART,JSTOP
    INTEGER PRETYP(KLEV)
    REAL(KIND=REALKIND) ABSERR, RELERR
    INTEGER IERR

    ! LOCAL VARIABLES

    REAL(KIND=REALKIND) SDOT(KLEV), ROM1(KLEV),ROM2(KLEV)
    REAL(KIND=REALKIND) RORO(KLEV), ALFA(KLEV)
    REAL(KIND=REALKIND) TYPICAL(KLEV), RMS(KLEV)
    REAL(KIND=REALKIND) SCALE, BETAM, DUMMY(KLON,KLAT,KKLEV)
    INTEGER ITER, LAT, LON, LEV, KLEVE
    REAL(KIND=REALKIND)    U(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) PHAT(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) VHAT(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) RHAT(KLON,KLAT,KLEV)

    ! SET HOW OFTEN TO CHECK FOR CONVERGENCE

    INTEGER    TSTITR
    PARAMETER (TSTITR = 1)

    ! SET SECURITY TRESHOLD FOR ITERATION

    INTEGER    MAXITR
    PARAMETER (MAXITR = 50)

    KLEVE = KKLEV

    ! SET BOUNDARIES TO ZERO

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          RHS (ISTART-1,LAT,LEV) = 0.0_realkind
          RHS (ISTOP +1,LAT,LEV) = 0.0_realkind
          PHAT(ISTART-1,LAT,LEV) = 0.0_realkind
          PHAT(ISTOP +1,LAT,LEV) = 0.0_realkind
       ENDDO
       DO LON = ISTART,ISTOP
          RHS (LON,JSTART-1,LEV) = 0.0_realkind
          RHS (LON,JSTOP +1,LEV) = 0.0_realkind
          PHAT(LON,JSTART-1,LEV) = 0.0_realkind
          PHAT(LON,JSTOP +1,LEV) = 0.0_realkind
       ENDDO
    ENDDO
    !     
    ! INITIALIZE ITERATIONS.
    ! SET FIRST GUESS EQUAL TO THE SOLUTION OF LAST TIME STEP.
    ! ALSO THE RIGHT HAND SIDE MUST BE SCALED IN THE SAME
    ! WAY AS THE AIJ-MATRIX WAS SCALED IN HHPRE.

    DO LEV=1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             SCALE = AMATP(LON,LAT,1,LEV)
             RHS(LON,LAT,LEV) =  RHS(LON,LAT,LEV) * SCALE
          ENDDO
       ENDDO
    ENDDO

    CALL AMULT(KLON,KLAT,KLEV,AMATP, RHS,R, KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP)

    DO LEV=1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             R(LON,LAT,LEV) =  RHS(LON,LAT,LEV) - R(LON,LAT,LEV)
             RHAT(LON,LAT,LEV) = R(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    ! START ITERATION 

    ! FIRST ITERATION IS HARDCODED TO AVOID IF-TESTS

    ITER = 1
    CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
         3,RHAT,R,RHS,KLEVE,ROM1,RMS,TYPICAL)

    ! CHECK FOR CONVERGENCE

    DO LEV=1,KLEVE
       IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
          KLEVE=LEV-1
          GOTO 10
       ENDIF
    ENDDO

10  CONTINUE
    IF (KLEVE == 0) THEN
       IERR = 0
       RETURN
    ENDIF

    ! END OF CHECK FOR CONVERGENCE

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             U(LON,LAT,LEV) = R(LON,LAT,LEV)
             P(LON,LAT,LEV) = U(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    CALL PRECOND(KLON,KLAT,KLEV, AMATP,P,PHAT,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP,PRETYP)

    CALL AMULT(KLON,KLAT,KLEV,AMATP,PHAT,VHAT,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP)

    CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
         1,RHAT,VHAT,DUMMY,KLEVE,SDOT,DUMMY,DUMMY)     

    DO LEV = 1,KLEVE
       ALFA(LEV) = ROM1(LEV)/SDOT(LEV)
       ROM2(LEV) = ROM1(LEV)
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             Q(LON,LAT,LEV) = U(LON,LAT,LEV)- ALFA(LEV)*VHAT(LON,LAT,LEV)
             U(LON,LAT,LEV) = U(LON,LAT,LEV) + Q(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    CALL PRECOND(KLON,KLAT,KLEV,AMATP,U,PHAT,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP,PRETYP)

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             RHS(LON,LAT,LEV) = RHS(LON,LAT,LEV) + ALFA(LEV)*PHAT(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    CALL AMULT(KLON,KLAT,KLEV, AMATP,PHAT,VHAT,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP)

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             R(LON,LAT,LEV) = R(LON,LAT,LEV) -ALFA(LEV)*VHAT(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    DO 100 ITER = 2, MAXITR


       ! POSSIBLY CHECK FOR STOPPING THE ITERATION.

       IF (MOD(ITER,TSTITR) == 0) THEN

          CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
               3,RHAT,R,RHS,KLEVE,ROM1,RMS,TYPICAL)

          ! CHECK FOR CONVERGENCE

          DO LEV=1,KLEVE
             IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
                KLEVE=LEV-1
                GOTO 20
             ENDIF
          ENDDO

20        CONTINUE
          IF (KLEVE == 0) THEN
             IERR = ITER -1
             RETURN
          ENDIF

          ! END OF CHECK FOR CONVERGENCE

       ELSE
          CALL HHGV(KLON,KLAT,KLEV, ISTART,ISTOP,JSTART,JSTOP, &
               1,RHAT,R,DUMMY,KLEVE,ROM1,DUMMY,DUMMY)
       ENDIF

       DO LEV = 1,KLEVE
          BETAM = ROM1(LEV)/ROM2(LEV)
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                U(LON,LAT,LEV) = R(LON,LAT,LEV) + BETAM * Q(LON,LAT,LEV)
                P(LON,LAT,LEV) = U(LON,LAT,LEV) + BETAM * (Q(LON,LAT,LEV) &
                     + BETAM *  P(LON,LAT,LEV))
             ENDDO
          ENDDO
       ENDDO

       CALL PRECOND(KLON,KLAT,KLEV,AMATP,P,PHAT,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP,PRETYP)

       CALL AMULT(KLON,KLAT,KLEV, AMATP,PHAT,VHAT,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP)

       CALL HHGV(KLON,KLAT,KLEV, ISTART,ISTOP,JSTART,JSTOP, &
            1,RHAT,VHAT,DUMMY,KLEVE,SDOT,DUMMY,DUMMY)

       DO LEV = 1,KLEVE
          ALFA(LEV) = ROM1(LEV)/SDOT(LEV)
          ROM2(LEV) = ROM1(LEV)
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                Q(LON,LAT,LEV) = U(LON,LAT,LEV) - ALFA(LEV)*VHAT(LON,LAT,LEV)
                U(LON,LAT,LEV) = U(LON,LAT,LEV) + Q(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

       CALL PRECOND(KLON,KLAT,KLEV,AMATP,U,PHAT,KLEVE,&
            ISTART,ISTOP,JSTART,JSTOP,PRETYP)

       DO LEV = 1,KLEVE
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                RHS(LON,LAT,LEV) = RHS(LON,LAT,LEV) +  ALFA(LEV)*PHAT(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

       CALL AMULT(KLON,KLAT,KLEV, AMATP,PHAT,VHAT,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP)

       DO LEV = 1,KLEVE
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                R(LON,LAT,LEV) = R(LON,LAT,LEV) -ALFA(LEV)*VHAT(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

100 enddo

    IERR = -1
    RETURN
  END SUBROUTINE CGSPRE

  SUBROUTINE RICHPRE(KLON,KLAT,KLEV,KKLEV, &
       AMATP,RHS, &
       P, Q, R, &
       ISTART,ISTOP,JSTART,JSTOP, &
       PRETYP, ABSERR, RELERR,IERR)
    !
    ! ITERATIVE HELMHOLTZ SOLVER
    ! PRECONDITIONED RICHARDSON METHOD
    !
    !
    !     D. BJOERGE      DNMI     1994-12-22, original version
    !                              1999-06-23, modified for HIRLAM MPP
    !                                                                       
    ! CALLED BY:
    ! ----------
    !
    !     HHSOLV
    !     DELNMI
    !
    !
    ! INPUT PARAMETERS:
    ! -----------------
    !
    !     KKLEV   NUMBER OF EQUATIONS TO SOLVE
    !     AMATP   MATRIX OF THE LEFT HAND SIDE
    !     RHS     RIGHT HAND SIDE OF THE HELMHOLTZ EQUATION
    !
    !     P,Q,R   WORK ARRAYS
    !     ABSERR   ABSOLUTE ERROR LIMIT FOR THE ITERATION
    !     RELERR   RELATIVE ERROR LIMIT FOR THE ITERATION
    !
    ! OUTPUT PARAMETERS:
    ! ------------------
    !
    !     IERR
    !           = -1    DID NOT CONVERGE
    !          >=  0    NUMBER OF ITERATIONS NEEDED
    !
    !     RHS  THE SOLUTION
    !
    !
    ! EXTERNALS
    ! ---------
    !
    !          AMULT   -- CALCULATES LEFT HAND SIDE OF THE
    !                     HELMHOLTZ EQUATION AX=B
    !          PRECOND -- PRECONDITIONING (SSOR OR NO PRECONDITIONER)
    !          HHGV    -- GLOBAL VECTOR OPERATIONS
    !
    !-------------------------------------------------------------------
    !

    IMPLICIT NONE

    ! GLOBAL VARIABLES

    INTEGER KLON, KLAT, KLEV
    REAL(KIND=REALKIND)    AMATP(KLON,KLAT,4,KLEV)
    INTEGER KKLEV
    REAL(KIND=REALKIND)    RHS(KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    P  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    Q  (KLON,KLAT,KKLEV)
    REAL(KIND=REALKIND)    R  (KLON,KLAT,KKLEV)
    INTEGER ISTART,ISTOP,JSTART,JSTOP
    INTEGER PRETYP(KLEV)
    REAL(KIND=REALKIND) ABSERR, RELERR
    INTEGER IERR

    ! LOCAL VARIABLES

    REAL(KIND=REALKIND) PQ(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND) ZRR0(KLEV), ZRR(KLEV), TYPICAL(KLEV), RMS(KLEV)
    REAL(KIND=REALKIND) SCALE
    REAL(KIND=REALKIND) DUMMY(KLON,KLAT,KKLEV)
    INTEGER ITER, IITER, LAT, LON, LEV, KLEVE

    ! SET HOW OFTEN TO CHECK FOR CONVERGENCE

    INTEGER    TSTITR
    PARAMETER (TSTITR = 3)

    ! SET SECURITY TRESHOLD FOR ITERATION

    INTEGER    MAXITR
    PARAMETER (MAXITR = 100)

    ! SET RELAXATION FACTOR

    REAL(KIND=REALKIND)       ALPHA
    PARAMETER (ALPHA = 0.8_realkind)

    KLEVE = KKLEV

    ! SET BOUNDARIES TO ZERO

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          RHS (ISTART-1,LAT,LEV) = 0.0_realkind
          RHS (ISTOP +1,LAT,LEV) = 0.0_realkind
       ENDDO
       DO LON = ISTART,ISTOP
          RHS (LON,JSTART-1,LEV) = 0.0_realkind
          RHS (LON,JSTOP +1,LEV) = 0.0_realkind
       ENDDO
    ENDDO

    ! INITIALIZE ITERATIONS.
    ! SET FIRST GUESS EQUAL TO THE SOLUTION OF LAST TIME STEP.
    ! ALSO THE RIGHT HAND SIDE MUST BE SCALED IN THE SAME
    ! WAY AS THE AMATP-MATRIX WAS SCALED IN HHPRE.

    DO LEV = 1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             SCALE = AMATP(LON,LAT,1,LEV)
             RHS(LON,LAT,LEV) = RHS(LON,LAT,LEV)*SCALE
             R(LON,LAT,LEV) = RHS(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    CALL AMULT(KLON,KLAT,KLEV, AMATP,RHS,Q,KLEVE,&
         ISTART,ISTOP,JSTART,JSTOP)

    DO LEV = 1,KLEVE
       DO LAT=JSTART,JSTOP
          DO LON=ISTART,ISTOP
             Q(LON,LAT,LEV) =   Q(LON,LAT,LEV) - R(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    ! START ITERATION 

    ! FIRST ITERATION IS HARDCODED TO AVOID IF-TESTS

    ITER = 1

    CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
         2,DUMMY,Q,RHS,KLEVE,DUMMY,RMS,TYPICAL)

    ! CHECK FOR CONVERGENCE

    DO LEV=1,KLEVE
       IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
          KLEVE = LEV-1
          GOTO 10
       ENDIF
    ENDDO

10  CONTINUE
    IF (KLEVE == 0) THEN
       IERR = 0
       RETURN
    ENDIF

    ! END OF CHECK FOR CONVERGENCE

    CALL PRECOND(KLON,KLAT,KLEV, AMATP,Q,PQ,KLEVE, &
         ISTART,ISTOP,JSTART,JSTOP,PRETYP)

    DO LEV = 1,KLEVE
       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             RHS(LON,LAT,LEV)=RHS(LON,LAT,LEV) - ALPHA*PQ(LON,LAT,LEV)
          ENDDO
       ENDDO
    ENDDO

    DO ITER=2,MAXITR

       CALL AMULT(KLON,KLAT,KLEV,AMATP,RHS,P,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP)

       DO LEV = 1,KLEVE
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                Q(LON,LAT,LEV) = P(LON,LAT,LEV) - R(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

       CALL PRECOND(KLON,KLAT,KLEV,AMATP,Q,PQ,KLEVE, &
            ISTART,ISTOP,JSTART,JSTOP,PRETYP)

       DO LEV = 1,KLEVE
          DO LAT = JSTART,JSTOP
             DO LON = ISTART,ISTOP
                RHS(LON,LAT,LEV)=RHS(LON,LAT,LEV) - ALPHA*PQ(LON,LAT,LEV)
             ENDDO
          ENDDO
       ENDDO

       ! POSSIBLY CHECK FOR STOPPING THE ITERATION.

       IF (MOD(ITER,TSTITR) == 0) THEN

          CALL HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
               2,DUMMY,PQ,RHS,KLEVE,DUMMY,RMS,TYPICAL)

          ! CHECK FOR CONVERGENCE

          DO LEV=1,KLEVE
             IF (RMS(LEV) < (ABSERR + RELERR*TYPICAL(LEV))) THEN
                KLEVE = LEV-1
                GOTO 20
             ENDIF
          ENDDO

20        CONTINUE
          IF (KLEVE == 0) THEN
             IERR = ITER - 1
             RETURN
          ENDIF

       ENDIF

       ! END OF CHECK FOR CONVERGENCE

    ENDDO

    ! REACHED 'MAXITR', REPORT ERROR.

    IERR = -1
    RETURN
  END SUBROUTINE RICHPRE

  SUBROUTINE PRECOND(KLON,KLAT,KLEV,AMATP,INN,UT,KLEVE, &
       ISTART,ISTOP,JSTART,JSTOP,PRETYP)
    !
    ! PRECONDITIONING FOR THE ITERATIVE HELMHOLTZ SOLVER
    !
    !
    !     D. BJOERGE            DNMI           1994-12-20
    !                                                                       
    !

    IMPLICIT NONE


    ! GLOBAL VARIABLES

    INTEGER KLON, KLAT, KLEV

    REAL(KIND=REALKIND) AMATP(KLON,KLAT,4,KLEV)

    REAL(KIND=REALKIND)   INN(KLON,KLAT,KLEV)
    REAL(KIND=REALKIND)    UT(KLON,KLAT,KLEV)

    INTEGER ISTART,ISTOP,JSTART,JSTOP
    INTEGER KLEVE
    INTEGER PRETYP(KLEV)

    ! LOCAL VARIABLES

    INTEGER I, J, LEV
    DO LEV=1,KLEVE

       ! SSOR PRECONDITION

       IF (PRETYP(LEV) == 1) THEN
          CALL SSOR(KLON,KLAT,INN(1,1,LEV),UT(1,1,LEV),AMATP(1,1,1,LEV), &
               ISTART,ISTOP,JSTART,JSTOP)

          ! NO PRECONDITION (I.E. IMPLICIT JACOBI PRECONDITION)

       ELSE
          DO J=JSTART,JSTOP
             DO I=ISTART,ISTOP
                UT(I,J,LEV) = INN(I,J,LEV)
             ENDDO
          ENDDO
       ENDIF

    ENDDO

    RETURN
  END SUBROUTINE PRECOND

  SUBROUTINE SSOR(KLON,KLAT,VIN,VOUT,AMAT, &
       ISTART,ISTOP,JSTART,JSTOP)
    !
    ! SSOR PRECONDITIONING FOR THE ITERATIVE HELMHOLTZ SOLVER
    !
    !
    !  DAG BJOERGE, DNMI           1994-12-20
    !
    !

    IMPLICIT NONE

    ! GLOBAL VARIABLES

    INTEGER KLON,KLAT
    REAL(KIND=REALKIND)   VIN(KLON,KLAT)
    REAL(KIND=REALKIND)  VOUT(KLON,KLAT)
    REAL(KIND=REALKIND)  AMAT(KLON,KLAT,4)
    INTEGER ISTART,ISTOP,JSTART,JSTOP

    ! LOCAL VARIABLES

    REAL(KIND=REALKIND) AIPMJ, AIJP, AIJM
    INTEGER I, J, INFO, LAT

    ! OVERRELAXATION FACTOR

    REAL(KIND=REALKIND)       OMEGA
    PARAMETER( OMEGA=1.2_realkind)

    ! SET UP BOUNDARIES

    DO I=ISTART,ISTOP
       VOUT(I,JSTART-1) = 0.0_realkind
    ENDDO

    DO J=JSTART,JSTOP
       VOUT(ISTART-1,J) = 0.0_realkind
    ENDDO

    DO J=JSTART,JSTOP
       DO I=ISTART,ISTOP
          AIPMJ=AMAT(I,J,2)
          AIJM =AMAT(I,J,4)
          VOUT(I,J) = OMEGA*(VIN(I,J) -(AIPMJ*VOUT(I-1,J) + AIJM*VOUT(I,J-1)))
       ENDDO
    ENDDO

    IF ( ATTOP ) THEN
       DO I=ISTART,ISTOP
          VOUT(I,JSTOP+1) = 0.0_realkind
       ENDDO
    ELSE
       DO I=ISTART,ISTOP
          VOUT(I,JSTOP+1) = VOUT(I,JSTOP)
       ENDDO
    ENDIF

    IF ( ATRIGHT ) THEN
       DO J=JSTART,JSTOP
          VOUT(ISTOP+1,J) = 0.0_realkind
       ENDDO
    ELSE
       DO J=JSTART,JSTOP
          VOUT(ISTOP+1,J) = VOUT(ISTOP,J)
       ENDDO
    ENDIF

    DO J=JSTOP,JSTART,-1
       DO I=ISTOP,ISTART,-1
          AIPMJ=AMAT(I,J,2)
          AIJP =AMAT(I,J,3)
          AIJM =AMAT(I,J,4)
          VOUT(I,J) = OMEGA*(VIN(I,J) - &
               ( AIPMJ*VOUT(I-1,J) + AIJM*VOUT(I,J-1) + &
               AIPMJ*VOUT(I+1,J) + AIJP*VOUT(I,J+1))) + &
               (1.0_realkind-OMEGA)*VOUT(I,J)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE SSOR

  SUBROUTINE AMULT(KLON,KLAT,KLEV, AMAT, INN,UT, KLEVE, &
       ISTART,ISTOP,JSTART,JSTOP)
    !
    ! MATRIX MULTIPLICATION FOR THE PARALLEL ITERATIVE HELMHOLTZ SOLVER
    !
    ! UT=AMAT*INN
    !                                                                       
    !     D. BJOERGE            DNMI     1994-12-20
    !     D. BJOERGE            DNMI     1995-11-07, UPDATE
    !
    !
    !

    IMPLICIT NONE

    ! GLOBAL VARIABLES

    INTEGER KLON,KLAT,KLEV,KLEVE
    REAL(KIND=REALKIND)  AMAT(KLON,KLAT,4,KLEV)
    REAL(KIND=REALKIND)   INN(KLON,KLAT,KLEVE)
    REAL(KIND=REALKIND)    UT(KLON,KLAT,KLEVE)
    INTEGER ISTART,ISTOP,JSTART,JSTOP

    ! LOCAL VARIABLES

    INTEGER LON, LAT, LEV
    REAL(KIND=REALKIND) AIPMJ, AIJP, AIJM, DUMMY

    ! EXCHANGE BOUNDARIES

    CALL SWAP(INN, KLON, KLAT, KLEVE)

    ! MATRIX MULTIPLY

    DO LEV=1,KLEVE

       DO LAT = JSTART,JSTOP
          DO LON = ISTART,ISTOP
             AIPMJ=AMAT(LON,LAT,2,LEV)
             AIJP =AMAT(LON,LAT,3,LEV)
             AIJM =AMAT(LON,LAT,4,LEV)
             UT(LON,LAT,LEV) = INN(LON,  LAT,  LEV) &
                  + AIPMJ * INN(LON-1,LAT,  LEV) &
                  + AIPMJ * INN(LON+1,LAT,  LEV) &
                  + AIJP  * INN(LON,  LAT+1,LEV) &
                  + AIJM  * INN(LON,  LAT-1,LEV)
          ENDDO
       ENDDO

    ENDDO

    RETURN
  END SUBROUTINE AMULT

  SUBROUTINE HHGV(KLON,KLAT,KLEV,ISTART,ISTOP,JSTART,JSTOP, &
       ITYPE,VEC1,VEC2,VEC3,KLEVE,prod,rms,typical)


    !     D. BJOERGE   DNMI   1994-12-22, original version
    !                         1999-06-23, modified for HIRLAM MPP
    !
    !
    ! purpose:
    ! --------
    !
    !       Perform global operations necessary in the iterative solvers.
    !
    !
    ! called by:
    ! ----------
    !
    !       GCPRE     --- iterative solvers ---
    !       CGSPRE
    !       RICHPRE
    !
    !
    ! input parameters:
    ! -----------------
    !
    !	vec1        the vectors involved
    !       vec2
    !       vec3
    !
    !	itype
    !	      = 1   Calculate vector-product of 'VEC1' and 'VEC2',
    !                   and return result in 'prod'.
    !
    !             = 2   Calculate rms of 'VEC2', and return result in 'rms'.
    !                   Calculate typical values of 'VEC3',
    !                   and return result in 'typical'.
    !
    !             = 3   A combination of itype=1 and itype=2.
    !
    !     	kleve       Third dimension of 'VECx' = number of
    !                   vector products, rms values and typical values
    !
    ! output parameters:
    ! ------------------
    !
    !	prod     The vector-products of 'VEC1' and 'VEC2'
    !       rms      Root mean square values of 'VEC2'
    !       typical  Maximum values of 'ABS(VEC3)'
    !
    !	Output parameters are undefined if not requested by itype.
    !
    !-----------------------------------------------------------------

    IMPLICIT NONE
    INTEGER KLON,KLAT,KLEV, ITYPE, KLEVE
    INTEGER ISTART,ISTOP,JSTART,JSTOP
    REAL(KIND=REALKIND) VEC1(KLON,KLAT,KLEVE)
    REAL(KIND=REALKIND) VEC2(KLON,KLAT,KLEVE)
    REAL(KIND=REALKIND) VEC3(KLON,KLAT,KLEVE)
    REAL(KIND=REALKIND) prod(KLEV), rms(KLEV), typical(KLEV)

    ! LOCAL VARIABLES

    REAL(KIND=REALKIND) VV, VVMAX, SUM
    INTEGER LAT, LON, LEV, INFO
    integer, parameter :: klev_global_max = 100
    REAL(KIND=REALKIND) SSUM(2*KLEV_GLOBAL_MAX), TTYP(KLEV_GLOBAL_MAX)

#ifdef MPI_SRC
#include"mpif.h"
    REAL(KIND=REALKIND) WORK(2*KLEV_GLOBAL_MAX)
#endif

    IF (ITYPE == 1 .OR. ITYPE == 3) THEN
       DO LEV=1,KLEVE	
          SUM = 0.0_realkind	
          DO LAT=JSTART,JSTOP
             DO LON=ISTART,ISTOP
                SUM = SUM + VEC1(LON,LAT,LEV) * VEC2(LON,LAT,LEV)
             ENDDO
          ENDDO
          SSUM(LEV) = SUM
       ENDDO
    ENDIF
    !
    IF (ITYPE == 2 .OR. ITYPE == 3) THEN
       DO LEV=1,KLEVE
          SUM   = 0.0_realkind
          VVMAX = 0.0_realkind
          DO LAT=JSTART,JSTOP
             DO LON=ISTART,ISTOP
                SUM = SUM + VEC2(LON,LAT,LEV) * VEC2(LON,LAT,LEV)
                VV  = ABS(VEC3(LON,LAT,LEV))
                IF (VV > VVMAX) VVMAX = VV
             ENDDO
          ENDDO
          SSUM(KLEVE+LEV) = SUM
          TTYP(LEV)      = VVMAX
       ENDDO

    ENDIF
#ifdef MPI_SRC
    IF (ITYPE == 1) THEN
       WORK(1:KLEVE) = SSUM(1:KLEVE)
       CALL MPI_ALLREDUCE(WORK, SSUM(1), KLEVE, REALTYPE, MPI_SUM,&
            localComm, INFO)
    ELSEIF (ITYPE == 2) THEN
       WORK(1:KLEVE) = SSUM(KLEVE+1:2*KLEVE)
       CALL MPI_ALLREDUCE(WORK, SSUM(KLEVE+1), KLEVE, REALTYPE, &
            MPI_SUM, localComm, INFO)
       WORK(1:KLEVE) = TTYP(1:KLEVE)
       CALL MPI_ALLREDUCE(WORK, TTYP, KLEVE, REALTYPE, &
            MPI_MAX, localComm, INFO)
    ELSEIF (ITYPE == 3) THEN
       WORK(1:2*KLEVE) = SSUM(1:2*KLEVE)
       CALL MPI_ALLREDUCE(WORK, SSUM(1), 2*KLEVE, REALTYPE, &
            MPI_SUM, localComm, INFO)
       WORK(1:KLEVE) = TTYP(1:KLEVE)
       CALL MPI_ALLREDUCE(WORK, TTYP, KLEVE, REALTYPE, &
            MPI_MAX, localComm, INFO)
    ELSE
       if(MYPE==0)write(*,*)' YOU CAN NOT CALL HHGV WITH ITYPE= ',ITYPE
       call stop_program('')
    ENDIF
#endif

    IF (ITYPE == 1 .OR. ITYPE == 3) THEN
       DO LEV = 1,KLEVE
          prod(LEV) = SSUM(LEV)
       ENDDO
    ENDIF

    IF (ITYPE == 2 .OR. ITYPE == 3) THEN
       DO LEV = 1,KLEVE
          rms(LEV) = sqrt(SSUM(KLEVE + LEV)) / real(KLON_GLOBAL*KLAT_GLOBAL,realkind)
          typical(LEV) = TTYP(LEV)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE HHGV


end module mod_implicit_solver
