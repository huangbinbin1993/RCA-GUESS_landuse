module coupling
  use decomp

#if defined(OASIS3)
   USE mod_kinds_model
   USE mod_prism_proto              ! OASIS3 prism module
   USE mod_prism_def_partition_proto! OASIS3 prism module for partitioning
   USE mod_prism_grids_writing
   USE mod_prism_put_proto          ! OASIS3 prism module for snding
   USE mod_prism_get_proto          ! OASIS3 prism module for receiving
   USE mod_comprism_proto           ! OASIS3 prism module to get coupling frequency
#endif

  use timetype
  use domainmod
!  use surface_bc,only:tsea
  implicit none

  private
  character(len=128),public,save :: model_name
  integer,public,save::comp_id
  integer,public,save::nnstep=99999
  type(time)::exchangeTime
  type(deltatime)::exchangeInterval
  logical,public,save:: lfirst_oasis=.true.

#if defined(OASIS3)
!    real(kind=8),public :: seconds_since_begin
    integer,public :: seconds_since_begin
#endif


!#if defined(OASIS3)
integer, allocatable,public :: mask (:,:)
integer,public         :: NSENDID(34),NRECVID(4)
character(len=8),public::  CDRECVID(4),CDSENDID(34)

!#endif

  real(kind=8), allocatable,public :: coords_lon (:, :)
  real(kind=8), allocatable,public :: coords_lat (:, :)
  real(kind=8), allocatable,public :: corners_lon (:,:,:)
  real(kind=8), allocatable,public :: corners_lat (:,:,:)


!#if defined(OASIS3)
  public read_namcoupling
  public couple_prism_init
  public couple_prism_define
  public couple_prism_recv
  public couple_prism_send 
!#endif


  public set_dqndt

  real(kind=realkind),public,allocatable,dimension(:,:)::sst_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::ict_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::icc_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::alb_fieldm !left to use
  real(kind=realkind),public,allocatable,dimension(:,:)::momuo_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::momvo_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::momui_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::momvi_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::netheato_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::netheati_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::solaro_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::solari_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:) ::dqndt_fieldm
!  real(kind=realkind),allocatable,dimension(:,:) ::dqndt_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::oemp_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::rain_local

  real(kind=realkind),public,allocatable,dimension(:,:)::t2ms_fieldm,t2mi_fieldm &
  ,evas_fieldm,shfs_fieldm,shfi_fieldm,lhfs_fieldm,lhfi_fieldm &
  ,swd_fieldm,swns_fieldm,lwd_fieldm,lwns_fieldm &
  ,u10s_fieldm,u10i_fieldm,v10s_fieldm,v10i_fieldm &
  ,rai_fieldm,sno_fieldm,rqa_fieldm,clo_fieldm,slp_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::u10s_fieldm1,u10i_fieldm1,v10s_fieldm1,v10i_fieldm1
 
  real(kind=realkind),public,allocatable,dimension(:,:)::amus_fieldm,amvs_fieldm,amui_fieldm,amvi_fieldm
  real(kind=realkind),public,allocatable,dimension(:,:)::tsea_save,tice_save,alb_save,frice_save


contains

!#if defined(OASIS3)



  subroutine read_namcoupling(inTime,deltaTimeStep)
    use decomp
    implicit none
    type(time),intent(in)::inTime
    type(deltatime),intent(in)::deltaTimeStep
#ifdef MPI_SRC
#include"mpif.h"
    integer ierr
#endif
    namelist /namcoupling/ nnstep !number of timesteps of whole simulation period 
    if(mype==0)then
       open(57,file='namelists.dat',status='old')
       read (57,nml=namcoupling) !coupling frequency
       close(57)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namcoupling)
       close(1)
       write(6,nml=namcoupling)
    endif

#ifdef MPI_SRC
    call mpi_bcast(nnstep,1,mpi_integer,0,localcomm,ierr)
#endif
    exchangeTime = inTime
    !exchangeInterval = nexchange*deltaTimeStep
  end subroutine read_namcoupling


!
! the follwoing part is for OASIS3 coupler
!

!#if defined(OASIS3)

  subroutine couple_prism_init(localcomm)
    implicit none
    integer, intent(out) :: localcomm
    integer              :: ierr
#if defined(OASIS3)

    model_name = 'rcaoa4'

      CALL prism_init_comp_proto (comp_id, 'rcaoa4', ierr)
!wangsy      if ( ierr .ne. 0) CALL abort
       write(6,*)'rca check after prism_init'
      CALL prism_get_localcomm_proto (localcomm,ierr)
       write(6,*)'check localComm  =',localComm,ierr

      IF (ierr .NE. 0) THEN
        WRITE (6,*) ' atm : pb prism_init_comp_proto'
         CALL prism_abort_proto(comp_id, 'atm.F90','abort2')
      ELSE
         WRITE(6,*) 'atm : prism_init_comp_proto ok '
      ENDIF

    write(6,*) 'check comp_id  =',comp_id

#else
    !not to get compilation errror
    localComm=999
#endif
  end subroutine couple_prism_init

  subroutine couple_prism_define(ndtime,klon,klat,meco, RCAdom, lcounttice,eco)
    use config
    use decomp
    use util
    implicit none
    type(domain)::RCAdom
    integer,intent(in)   :: klon,klat,meco
    real(kind=realkind), dimension(klon,klat), intent(inout) :: lcounttice
    real(kind=realkind),dimension(klon,klat,meco), target,intent(inout) :: eco
    real(kind=realkind),dimension(:,:),allocatable::frlandglobal,frlakeglobal
    real(kind=realkind),dimension(:,:), pointer:: fr_lake,fr_land
    integer:: ierror,i,j,k,lunut
    integer:: ndtime
    real(kind=realkind):: rmin,rmax,aplon,aplat, south,west,dlamda,dtheta
    real(kind=realkind),dimension(1,1)::rotlat,rotlon,reglon,reglat
    real(kind=realkind)::xxmean,newtime

     real(kind=realkind),dimension(1:klon,1:klat) :: zarea
!     integer,dimension(1:klon,1:klat) :: mask
!     integer         :: INODIM(2),IPSHAPE(2,2),IPARTID,NSENDID(10),NRECVID(4)
!     character(len=8)::  CDRECVID(4),CDSENDID(10)
    integer         :: INODIM(2),IPSHAPE(2,2),IPARTID
    integer::status,ierr,igparal(5)

    integer      :: masktmp(klon_global, klat_global)
    integer      :: iextent,jextent
    real(kind=8) :: masknew(klon,klat)
    character(len=11) :: form191
    integer      :: iline(1)

#if defined(OASIS3)

#ifdef MPI_SRC
#include "mpif.h"
#endif
    allocate(frlandglobal(klon_global,klat_global),&
         frlakeglobal(klon_global,klat_global),stat=status)
    if(status/=0)then
       stop 'memory allocation errror klon_global*klat_global too large?'
    endif

    allocate(sst_fieldm(klon,klat),ict_fieldm(klon,klat),&
         icc_fieldm(klon,klat),alb_fieldm(klon,klat))
    allocate(dqndt_fieldm(klon,klat))
    allocate(mask(klon,klat))

    allocate(momuo_fieldm(klon,klat),momvo_fieldm(klon,klat),momui_fieldm(klon,klat),  &
            momvi_fieldm(klon,klat),netheato_fieldm(klon,klat),netheati_fieldm(klon,klat),  &
            solaro_fieldm(klon,klat),solari_fieldm(klon,klat),oemp_fieldm(klon,klat))

    if(use_oasis_arctic) then
    allocate(t2ms_fieldm(klon,klat),t2mi_fieldm(klon,klat) &
            ,evas_fieldm(klon,klat),shfs_fieldm(klon,klat),shfi_fieldm(klon,klat),lhfs_fieldm(klon,klat),lhfi_fieldm(klon,klat) &
            ,swd_fieldm(klon,klat),swns_fieldm(klon,klat),lwd_fieldm(klon,klat),lwns_fieldm(klon,klat) &
            ,u10s_fieldm(klon,klat),u10i_fieldm(klon,klat),v10s_fieldm(klon,klat),v10i_fieldm(klon,klat) &
            ,rai_fieldm(klon,klat),sno_fieldm(klon,klat),rqa_fieldm(klon,klat),clo_fieldm(klon,klat),slp_fieldm(klon,klat))
    allocate(u10s_fieldm1(klon,klat),u10i_fieldm1(klon,klat),  &
             v10s_fieldm1(klon,klat),v10i_fieldm1(klon,klat))
    allocate(amus_fieldm(klon,klat),amvs_fieldm(klon,klat),amui_fieldm(klon,klat),amvi_fieldm(klon,klat))
    allocate(tsea_save(klon,klat),tice_save(klon,klat),alb_save(klon,klat),frice_save(klon,klat))
    endif

    allocate (coords_lon(klon,klat),coords_lat(klon,klat),corners_lon(klon,klat,4),corners_lat(klon,klat,4))


    fr_lake => eco(:,:,42)
    fr_land => eco(:,:,10)

    lunut=6
    if(use_oasis)then
       if (lfirst_oasis) then
!wangsy          where(lcounttice<0.5) fr_land=0._realkind

!wangsy          call colfld(0,frlandglobal,fr_land,klon,klat)
!wangsy          call colfld(0,frlakeglobal,fr_lake,klon,klat)

            if(use_oasis_arctic) then
    
            if ( mype == 0) then    
               masktmp=0
               write(6,*) 'lsm_ read common mask,mype',mype
               open(191,file="mask4atm_arctic_glimpse_1.asc",form="formatted")
!	       write(form191,'(a4,i3,a3)') '(4x,',klon_global,'i3)'
	       write(form191,'(a4,i3,a3)') '(i4,',klon_global,'i3)'
	       print*,'form191 ',form191,klon_global,klat_global
	       do j=klat_global,1,-1
		 write(6,*) 'j ',j,1+(klon_global*(j-1)),klon_global
		 read(191,fmt=form191) iline,(masktmp(i,j), i=1,klon_global)
		 write(6,*) 'lsm ',iline,(masktmp(i,j), i=1,10)
		 call flush(6)
	       enddo
	       close(191)
	       print*,'masktmp(1,1)  ',masktmp(1,1),masktmp(1,2)
	       print*,'a-masktmp(60,50)  ',masktmp(60,50)
	       
	       !convert
               do j=1,klat_global
               do i=1,klon_global
                  masktmp(i,j)=abs(masktmp(i,j)-1)
               enddo
               enddo
	       print*,'b-masktmp(60,50)  ',masktmp(60,50)
	    endif
	    	    
            call mpi_bcast(masktmp,klon_global*klat_global,mpi_integer,0,localComm,ierr)

            do j=1,klat
            do i=1,klon
              masknew(i,j)=real(masktmp(i+idatastart-1,j+jdatastart-1),realkind)
            enddo
            enddo
            print*,'masknew(1,1)  ',masknew(1,1),masknew(1,2)
	    print*,'c-masktmp(60,50)  ',masktmp(60,50)
    
            else   !use_oasis_arctic
    
!
!==================================
!wangsy
!  read in new mask
             if ( mype == 0) then
              open(23,file='rca_cou_mask.txt',form='formatted')

              do j=klat_global,1,-1
!              read(23,'(456I1)')(masktmp(i,j),i=1,klon_global)
              read(23,'(186I1)')(masktmp(i,j),i=1,klon_global)
              enddo
    
              close(23)
             endif

             call mpi_bcast(masktmp,klon_global*klat_global,mpi_integer,0,localComm,ierr)

             do j=1,klat
             do i=1,klon
                masknew(i,j)=real(masktmp(i+idatastart-1,j+jdatastart-1),realkind)

!                if(masknew(i,j) > 0) lcounttice(i,j)=0.

             enddo
             enddo
	    endif   !use_oasis_arctic
!====================================


#ifdef notusednow
             do j=1,klat
                rotlat(1,1)=RCAdom%south + (jdatastart+j-1-1)*RCAdom%dlat
                do i=1,klon
                   rotlon(1,1)=RCAdom%west + (idatastart+i-1-1)*RCAdom%dlon
                   call regrot(reglon,reglat,rotlon,rotlat,1,1,RCAdom%polon,RCAdom%polat,-1)
                   coords_lon(i,j)=reglon(1,1)
                   coords_lat(i,j)=reglat(1,1)
                   if(use_oasis_arctic) then
                     coords_lon(i,j)=rotlon(1,1)
                     coords_lat(i,j)=rotlat(1,1)
                   endif
                   write(97,'(i4,2x,3f14.6)')j,reglon,reglat
                   rmin=min(rmin,coords_lon(i,j))
                   rmax=max(rmax,coords_lon(i,j))
                enddo
             enddo

            zarea(i,j)=1.0

            do j=1,klat
                rotlat(1,1)=RCAdom%south + (jdatastart+j-1-1.5)*RCAdom%dlat
                do i=1,klon
                   rotlon(1,1)=RCAdom%west + (idatastart+i-1-1.5)*RCAdom%dlon

                   call regrot(reglon,reglat,rotlon,rotlat,1,1,RCAdom%polon,RCAdom%polat,-1)
                   corners_lon(i,j,1)=reglon(1,1)
                   corners_lat(i,j,1)=reglat(1,1)
                   if(use_oasis_arctic) then
                     corners_lon(i,j,1)=rotlon(1,1)
                     corners_lat(i,j,1)=rotlat(1,1)
                   endif
                enddo
             enddo

             do j=1,klat
                rotlat(1,1)=RCAdom%south + (jdatastart+j-1-1.5)*RCAdom%dlat
                do i=1,klon
                   rotlon(1,1)=RCAdom%west + (idatastart+i-1-0.5)*RCAdom%dlon

                   call regrot(reglon,reglat,rotlon,rotlat,1,1,RCAdom%polon,RCAdom%polat,-1)
                   corners_lon(i,j,2)=reglon(1,1)
                   corners_lat(i,j,2)=reglat(1,1)
                   if(use_oasis_arctic) then
                     corners_lon(i,j,2)=rotlon(1,1)
                     corners_lat(i,j,2)=rotlat(1,1)
                   endif
                enddo
             enddo

             do j=1,klat
                rotlat(1,1)=RCAdom%south + (jdatastart+j-1-0.5)*RCAdom%dlat
                do i=1,klon
                   rotlon(1,1)=RCAdom%west + (idatastart+i-1-0.5)*RCAdom%dlon
                   call regrot(reglon,reglat,rotlon,rotlat,1,1,RCAdom%polon,RCAdom%polat,-1)
                   corners_lon(i,j,3)=reglon(1,1)
                   corners_lat(i,j,3)=reglat(1,1)
                   if(use_oasis_arctic) then
                     corners_lon(i,j,3)=rotlon(1,1)
                     corners_lat(i,j,3)=rotlat(1,1)
                   endif
                enddo
             enddo
             !     upper left
             do j=1,klat
                rotlat(1,1)=RCAdom%south + (jdatastart+j-1-0.5)*RCAdom%dlat
                do i=1,klon
                   rotlon(1,1)=RCAdom%west + (idatastart+i-1-1.5)*RCAdom%dlon

                   call regrot(reglon,reglat,rotlon,rotlat,1,1,RCAdom%polon,RCAdom%polat,-1)
                   corners_lon(i,j,4)=reglon(1,1)
                   corners_lat(i,j,4)=reglat(1,1)
                   if(use_oasis_arctic) then
                     corners_lon(i,j,4)=rotlon(1,1)
                     corners_lat(i,j,4)=rotlat(1,1)
                   endif
                enddo
             enddo
#endif  ! notusednow

             do j=1,klat
                do i=1,klon
!                   if (eco(i,j,10) < 0.5) then
                   if (masknew(i,j) > 0.5) then
!wangsy                      mask_array(i,j,1)=.true.
                      mask(i,j)=1
                      lcounttice(i,j)=0.
                   else
                      mask(i,j)=0
                   endif
                enddo
             enddo

          write(6,*)'before start grid wiriting'
!          CALL PRISM_START_GRIDS_WRITING(IERR)
!          write(6,*)'check flag grid write ierr  =',ierr
!          CALL PRISM_WRITE_GRID("atmo",klon,klat,coords_lon,coords_lat)
!          CALL PRISM_WRITE_CORNER('atmo',klon,klat,4,corners_lon,corners_lat)

!          CALL PRISM_WRITE_MASK('atmo',klon,klat,mask)

!          CALL PRISM_WRITE_AREA('atmo',klon,klat,zarea)

!          CALL PRISM_TERMINATE_GRIDS_WRITING
          write(6,*)'after grids writing'

         igparal(1) = 2              ! box partitioning
         if(atbase) then
         igparal(2) = idatastart-1
         else
!         igparal(2) = klon_global*(jdatastart-2)+idatastart-1              ! lower left corner global offset
         igparal(2) = klon_global*(jdatastart-1)+idatastart-1              ! lower left corner global offset
         endif
         igparal(3) = klon         ! local extent in i
         igparal(4) = klat         ! local extent in j

!          igparal(2) = klon_global*(jdatastart-1)+idatastart
!          if(atleft) then
!          igparal(2) = klon_global*(jdatastart-1)+idatastart-1
!          endif
!          iextent=klon-2
!          jextent=klat-2
!          if(atbase .or. attop)then
!          jextent = klat-1
!          endif
!          if(atleft .or. atright)then
!          iextent = klon-1
!          endif
!         igparal(3) =iextent
!         igparal(4) =jextent

         igparal(5) = klon_global         ! global extent in x



          write(6,*) 'rca before partition definition'
          call PRISM_DEF_PARTITION_PROTO(IPARTID,IGPARAL,IERR)
          if (ierr .ne.0) then
             write(6,*)'OASIS3 partition definition failed',mype
             write(6,'(A,I8)')'Return code : ',IERR
             call flush(6)
             call abort
           else
            write(6,*)'rca partition define is ok'
          endif

          INODIM(1)=2
          INODIM(2)=0

          IPSHAPE(1,1) =1
          IPSHAPE(2,1) =klon
          IPSHAPE(1,2) =1
          IPSHAPE(2,2) =klat

          CDSENDID(1)='taux_oce'
          CDSENDID(2)='tauy_oce'
          CDSENDID(3)='taux_ice'
          CDSENDID(4)='tauy_ice'
          CDSENDID(5)='qsol_oce'
          CDSENDID(6)='qsol_ice'
          CDSENDID(7)='qnet_oce'
          CDSENDID(8)='qnet_ice'
          CDSENDID(9)='oemp_oce'
          CDSENDID(10)='qndt_oce'

       if(.not.use_oasis_arctic) then
         do i=1,10
          call PRISM_DEF_VAR_PROTO(NSENDID(I),CDSENDID(I),IPARTID,  &
          INODIM,PRISM_Out,IPSHAPE,PRISM_Real,IERR)
         enddo
       endif

       if(use_oasis_arctic) then
          CDSENDID(11)='t2ms_a'
          CDSENDID(12)='t2mi_a'
          CDSENDID(13)='evas_a'
          CDSENDID(14)='shfs_a'
          CDSENDID(15)='shfi_a'
          CDSENDID(16)='lhfs_a'
          CDSENDID(17)='lhfi_a'
          CDSENDID(18)='swd_a'
          CDSENDID(19)='swns_a'
          CDSENDID(20)='lwd_a'
          CDSENDID(21)='lwns_a'
          CDSENDID(22)='u10s_a'
          CDSENDID(23)='u10i_a'
          CDSENDID(24)='v10s_a'
          CDSENDID(25)='v10i_a'
          CDSENDID(26)='rai_a'
          CDSENDID(27)='sno_a'
          CDSENDID(28)='rqa_a'
          CDSENDID(29)='clo_a'
          CDSENDID(30)='slp_a'
          CDSENDID(31)='amus_a'
          CDSENDID(32)='amvs_a'
          CDSENDID(33)='amui_a'
          CDSENDID(34)='amvi_a'
         do i=11,34
          call PRISM_DEF_VAR_PROTO(NSENDID(I),CDSENDID(I),IPARTID,  &
          INODIM,PRISM_Out,IPSHAPE,PRISM_Real,IERR)
         enddo
       endif

 
          CDRECVID(1)='SST__ATM'
          CDRECVID(2)='ICT__ATM'
          CDRECVID(3)='ALB__ATM'
          CDRECVID(4)='ICC__ATM'
          do i=1,4
          CALL PRISM_DEF_VAR_PROTO(NRECVID(i),CDRECVID(i),IPARTID,  &
          INODIM,PRISM_In,IPSHAPE,PRISM_Real,IERR)
          enddo


           CALL PRISM_ENDDEF_PROTO(IERR)

          write(6,*)'after definition',mype,IERR

          lfirst_oasis=.false.
       endif                     ! lfirst_oasis
    endif !use_oasis
    deallocate(frlandglobal,frlakeglobal)
    sst_fieldm = 5.0_realkind + 273.155_realkind
    ict_fieldm = 0.0_realkind + 273.155_realkind
    alb_fieldm = 0.1_realkind
    icc_fieldm = 0.0_realkind

     momuo_fieldm = 0.1_realkind
     momvo_fieldm = 0.1_realkind
     momui_fieldm = 0.1_realkind
     momvi_fieldm = 0.1_realkind
     netheato_fieldm = -10.0_realkind
     netheati_fieldm = -10.0_realkind
     solaro_fieldm = 20.0_realkind
     solari_fieldm = 20.0_realkind
     dqndt_fieldm = -10.0_realkind
     oemp_fieldm = 0.0_realkind

     if(use_oasis_arctic) then
       t2ms_fieldm(:,:)=273.155_realkind
       t2mi_fieldm = 273.155_realkind
       evas_fieldm = 1.0_realkind
       shfs_fieldm = 1.0_realkind
       shfi_fieldm = 1.0_realkind
       lhfs_fieldm = 1.0_realkind
       lhfi_fieldm = 1.0_realkind
       swd_fieldm = 1.0_realkind
       swns_fieldm = 1.0_realkind
       lwd_fieldm = 1.0_realkind
       lwns_fieldm = 1.0_realkind
       u10s_fieldm = 1.0_realkind
       u10i_fieldm = 1.0_realkind
       v10s_fieldm = 1.0_realkind
       v10i_fieldm = 1.0_realkind
       rai_fieldm = 1.0_realkind
       sno_fieldm = 1.0_realkind
       rqa_fieldm = 1.0_realkind
       clo_fieldm = 1.0_realkind
       slp_fieldm = 100000.0_realkind     
       amus_fieldm = 1.0_realkind
       amvs_fieldm = 1.0_realkind
       amui_fieldm = 1.0_realkind
       amvi_fieldm = 1.0_realkind
     endif   ! use_oasis_arctic
#endif   ! OASIS3


  end subroutine couple_prism_define

  subroutine couple_prism_recv(klon,klat,tsea,frice,tice,alb,ndtime)
    implicit none
    integer   :: i,j
    integer   :: klon,klat,ndtime
    real(kind=8),allocatable,dimension(:,:)::oas_r8_fieldm
    real(kind=realkind),dimension(klon,klat)::tsea,alb
    real(kind=realkind),dimension(klon,klat)::frice,tice
    integer   :: lunout_oas,info,ierror
    
#ifdef OASIS3    
    
    lunout_oas =6
    write(lunout_oas,*) 'atmos: read from coupler',seconds_since_begin
    call flush(lunout_oas)

    allocate(oas_r8_fieldm(klon,klat))
    !  get the data from the ocean (sst)
    ! sst ocean -> atm
!oas_r8_fieldm=-1.0*mype
oas_r8_fieldm=tsea

   call prism_get_proto(NRECVID(1),seconds_since_begin,oas_r8_fieldm,info)

!wsy    call prism_get (prism_var_id(1), model_time,model_time_bounds, oas_r8_fieldm, info, ierror)
    write(lunout_oas,*) 'info error ',info,ierror
!    sst_fieldm=real(oas_r8_fieldm,realkind)
!    tsea=real(oas_r8_fieldm,realkind)

    if(info==3) then   ! if a new field comes in
      sst_fieldm=real(oas_r8_fieldm,realkind)

      do j=1,klat
      do i=1,klon 
      if(mask(i,j)==1) then
      tsea(i,j)=real(oas_r8_fieldm(i,j),realkind)
      endif
      enddo
      enddo
    else
      do j=1,klat
      do i=1,klon
        if(mask(i,j)==1) then
          tsea(i,j)=sst_fieldm(i,j)
        endif
      enddo
      enddo
    endif

    ! ict ocean -> atm
oas_r8_fieldm=tice
   call prism_get_proto (NRECVID(2),seconds_since_begin,oas_r8_fieldm,info)
!wsy    call prism_get (prism_var_id(2), model_time,model_time_bounds, oas_r8_fieldm, info, ierror)
    write(lunout_oas,*) 'info error ',info,ierror
!    ict_fieldm=real(oas_r8_fieldm,realkind)
!    tice=real(oas_r8_fieldm,realkind)
    if(info==3) then   ! if a new field comes in
      ict_fieldm=real(oas_r8_fieldm,realkind)

      do j=1,klat
      do i=1,klon
        if(mask(i,j)==1) then
        tice(i,j)=real(oas_r8_fieldm(i,j),realkind)
        endif
      enddo
      enddo
    else
      do j=1,klat
      do i=1,klon
!        if(mask_array(i,j,1)) then
        if(mask(i,j)==1) then
          tice(i,j)=ict_fieldm(i,j)
        endif
      enddo
      enddo
    endif


    ! alb ocean -> atm
 oas_r8_fieldm=alb
   call prism_get_proto (NRECVID(3),seconds_since_begin,oas_r8_fieldm,info)
!wsy    call prism_get (prism_var_id(3), model_time,model_time_bounds, oas_r8_fieldm, info, ierror)
    write(lunout_oas,*) 'info error ',info,ierror
!    alb_fieldm=real(oas_r8_fieldm,realkind)
    if(info==3) then   ! if a new field comes in
      alb_fieldm=real(oas_r8_fieldm,realkind)
      do j=1,klat
      do i=1,klon
        if(mask(i,j)==1) then
        alb(i,j)=real(oas_r8_fieldm(i,j),realkind)
        endif
      enddo
      enddo
    else
      do j=1,klat
      do i=1,klon
!        if(mask_array(i,j,1)) then
        if(mask(i,j)==1) then
          alb(i,j)=alb_fieldm(i,j)
        endif
      enddo
      enddo
    endif

    ! icc ocean -> atm
 oas_r8_fieldm=frice
   call prism_get_proto (NRECVID(4),seconds_since_begin,oas_r8_fieldm,info)
!wsy    call prism_get (prism_var_id(4), model_time,model_time_bounds, oas_r8_fieldm, info, ierror)
    write(lunout_oas,*) 'info error ',info,ierror
!    icc_fieldm=real(oas_r8_fieldm,realkind)
!    frice=real(oas_r8_fieldm,realkind)
    if(info==3) then   ! if a new field comes in
      icc_fieldm=real(oas_r8_fieldm,realkind)
      do j=1,klat
      do i=1,klon
      if(mask(i,j)==1) then
      frice(i,j)=real(oas_r8_fieldm(i,j),realkind)
        if(i.eq.144.and.j.eq.114) then
          write(6,'("-a-i,j,mask,frice,mask",2(1x,i4),1x,f8.2,1x,i4)') i,j,mask,frice,info
        endif
      endif
      enddo
      enddo
    else
      do j=1,klat
      do i=1,klon
!        if(mask_array(i,j,1)) then
        if(mask(i,j)==1) then
          frice(i,j)=icc_fieldm(i,j)
          if(i.eq.144.and.j.eq.114) then
            write(6,'("-b-i,j,mask,frice,mask",2(1x,i4),1x,f8.2,1x,i4)') i,j,mask,frice,info
          endif
        endif
      enddo
      enddo
    endif
    deallocate(oas_r8_fieldm)

#endif   ! OASIS3

  end subroutine couple_prism_recv

  subroutine couple_prism_send(klon,klat,ndtime)

    use decomp
    use config,only:use_oasis_arctic,use_oasis_netcdf_output
    use netcdf

    implicit none
    integer,intent(in)::ndtime
    integer  :: timestep,nstart1,info,ierror
    integer  :: i,j,k
    integer,intent(in):: klon,klat
!EM Variable for netcdf coupling fields output
    integer  :: ier, status, il_file_id, dimids(2), il_var_id(34), ib
    character(len=8)::  csendname(34)

    real(kind=realkind):: t2ms_fieldmin,t2ms_fieldmax
    real(kind=8)::oas_r8_fieldm(klon,klat)

#ifdef OASIS3

    write(6,*) 'atmos: send to coupler',seconds_since_begin
    write(6,*) 'mask: ', mask(1,1)




   write(6,*) 'atmos: send to coupler',seconds_since_begin

    ! write to coupler

    !C          write the field at the location agreed so that it can be
    !C          read by the coupler at a further phase. Use of the oasis
    !C          locwrite routine is made.
    !C
    !c --- t2m atm -> ocean ---

!EM Coupling field output





!    IF (mype == 0) then
!       status = nf90_create("rca_cpl.nc", NF90_CLOBBER, il_file_id)
!       status = nf90_def_dim(il_file_id, "x", klon_global, dimids(1))
!       status = nf90_def_dim(il_file_id, "y", klat_global, dimids(2))
!       status = nf90_def_var(il_file_id, "qsol_oce", NF90_DOUBLE, dimids, il_var_id(1))
!       status = nf90_enddef(il_file_id)
!       status = nf90_close(il_file_id)
!    ENDIF

!    DO ib = 0, nproc - 1

!       call MPI_Barrier(localcomm, ierror)

!       IF (mype == ib) THEN
!          status = nf90_open("rca_cpl.nc", NF90_WRITE, il_file_id)
!          status = nf90_inq_varid(il_file_id, "qsol_oce" , il_var_id(1))
!          oas_r8_fieldm=DBLE(solaro_fieldm)
!          oas_r8_fieldm=20.*mask
!          status = nf90_put_var(il_file_id, il_var_id(1), oas_r8_fieldm, &
!                   start = (/ idatastart , jdatastart /), &
!                   count = (/ klon, klat /) )
!          status = nf90_close(il_file_id)
!       ENDIF
!    ENDDO

!EM Coupling field output

    !  Send the data to the ocean

    if(.not.use_oasis_arctic) then
    oas_r8_fieldm=DBLE(momuo_fieldm)
!    oas_r8_fieldm=0.1*mask
!    oas_r8_fieldm=0.1
     CALL PRISM_PUT_PROTO(NSENDID(1),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(momvo_fieldm)
!    oas_r8_fieldm=0.1*mask
     CALL PRISM_PUT_PROTO(NSENDID(2),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(momui_fieldm)
!    oas_r8_fieldm=0.1*mask
     CALL PRISM_PUT_PROTO(NSENDID(3),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(momvi_fieldm)
!    oas_r8_fieldm=0.1*mask
     CALL PRISM_PUT_PROTO(NSENDID(4),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(solaro_fieldm)
!    oas_r8_fieldm=20.0*mask
     CALL PRISM_PUT_PROTO(NSENDID(5),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(solari_fieldm)
!    oas_r8_fieldm=20.0*mask
     CALL PRISM_PUT_PROTO(NSENDID(6),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(netheato_fieldm)
!    oas_r8_fieldm=-10.0*mask
     CALL PRISM_PUT_PROTO(NSENDID(7),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(netheati_fieldm)
!    oas_r8_fieldm=-10.0*mask
     CALL PRISM_PUT_PROTO(NSENDID(8),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(oemp_fieldm)
!    oas_r8_fieldm=0.0
     CALL PRISM_PUT_PROTO(NSENDID(9),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(dqndt_fieldm)
!    oas_r8_fieldm=-10.0*mask
     CALL PRISM_PUT_PROTO(NSENDID(10),seconds_since_begin,oas_r8_fieldm,INFO)
    endif


    if(use_oasis_arctic) then
    oas_r8_fieldm=DBLE(t2ms_fieldm)
    write(6,*)' before prism_put t2ms',seconds_since_begin
    CALL PRISM_PUT_PROTO(NSENDID(11),seconds_since_begin,oas_r8_fieldm,INFO) 
    write(6,*)' rca after prism_put t2ms',t2ms_fieldm(10,10),info
    write(6,*)' rca after prism_put t2ms',t2ms_fieldm(10,10),seconds_since_begin

    oas_r8_fieldm=DBLE(t2mi_fieldm)
    do j=1,klat
    do i=1,klon
      if(t2mi_fieldm(i,j)<100.0) oas_r8_fieldm(i,j)=DBLE(t2ms_fieldm(i,j))    ! to prevent 99 dummy values to be send
    enddo
    enddo
    CALL PRISM_PUT_PROTO(NSENDID(12),seconds_since_begin,oas_r8_fieldm,INFO) 

    write(6,*)' rca after prism_put t2mi',t2mi_fieldm(10,10),info
    write(6,*)' rca after prism_put t2mi',t2mi_fieldm(10,10),seconds_since_begin
    
    oas_r8_fieldm=DBLE(evas_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(13),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(shfs_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(14),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(shfi_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(15),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(lhfs_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(16),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(lhfi_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(17),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(swd_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(18),seconds_since_begin,oas_r8_fieldm,INFO)

    oas_r8_fieldm=DBLE(swns_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(19),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(lwd_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(20),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(lwns_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(21),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(u10s_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(22),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(u10i_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(23),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(v10s_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(24),seconds_since_begin,oas_r8_fieldm,INFO)
 
    oas_r8_fieldm=DBLE(v10i_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(25),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(rai_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(26),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(sno_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(27),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(rqa_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(28),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(clo_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(29),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(slp_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(30),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(amus_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(31),seconds_since_begin,oas_r8_fieldm,INFO)
 
    oas_r8_fieldm=DBLE(amvs_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(32),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(amui_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(33),seconds_since_begin,oas_r8_fieldm,INFO)
    
    oas_r8_fieldm=DBLE(amvi_fieldm)
    CALL PRISM_PUT_PROTO(NSENDID(34),seconds_since_begin,oas_r8_fieldm,INFO)


    if(use_oasis_netcdf_output) then
    write(6,*) 'atmos: send to coupler',seconds_since_begin,INFO
    if(seconds_since_begin==10800) then
    if(INFO==4) then

       csendname(11)='t2ms_fieldm'
       csendname(12)='t2mi_fieldm'
       csendname(13)='evas_fieldm'
       csendname(14)='shfs_fieldm'
       csendname(15)='shfi_fieldm'
       csendname(16)='lhfs_fieldm'
       csendname(17)='lhfi_fieldm'
       csendname(18)='swd_fieldm'
       csendname(19)='swns_fieldm'
       csendname(20)='lwd_fieldm'
       csendname(21)='lwns_fieldm'
       csendname(22)='u10s_fieldm'
       csendname(23)='u10i_fieldm'
       csendname(24)='v10s_fieldm'
       csendname(25)='v10i_fieldm'
       csendname(26)='rai_fieldm'
       csendname(27)='sno_fieldm'
       csendname(28)='rqa_fieldm'
       csendname(29)='clo_fieldm'
       csendname(30)='slp_fieldm'
       csendname(31)='amus_fieldm'
       csendname(32)='amvs_fieldm'
       csendname(33)='amui_fieldm'
       csendname(34)='amvi_fieldm'

    IF (mype == 0) then
       status = nf90_create("rca_cpl_out.nc", NF90_CLOBBER, il_file_id)
       status = nf90_def_dim(il_file_id, "x", klon_global, dimids(1))
       status = nf90_def_dim(il_file_id, "y", klat_global, dimids(2))
       do i=11,34
         status = nf90_def_var(il_file_id, csendname(i), NF90_DOUBLE, dimids, il_var_id(i))
       enddo
       status = nf90_enddef(il_file_id)
       status = nf90_close(il_file_id)
    ENDIF

    DO ib = 0, nproc - 1

       call MPI_Barrier(localcomm, ierror)

       IF (mype == ib) THEN
	  status = nf90_open("rca_cpl_out.nc", NF90_WRITE, il_file_id)
	  
       do i=11,34
	  status = nf90_inq_varid(il_file_id, csendname(i) , il_var_id(i))
	  if(i==11) oas_r8_fieldm=DBLE(t2ms_fieldm)
	  if(i==12) oas_r8_fieldm=DBLE(t2mi_fieldm)
	  if(i==13) oas_r8_fieldm=DBLE(evas_fieldm)
	  if(i==14) oas_r8_fieldm=DBLE(shfs_fieldm)
	  if(i==15) oas_r8_fieldm=DBLE(shfi_fieldm)
	  if(i==16) oas_r8_fieldm=DBLE(lhfs_fieldm)
	  if(i==17) oas_r8_fieldm=DBLE(lhfi_fieldm)
	  if(i==18) oas_r8_fieldm=DBLE(swd_fieldm)
	  if(i==19) oas_r8_fieldm=DBLE(swns_fieldm)
	  if(i==20) oas_r8_fieldm=DBLE(lwd_fieldm)
	  if(i==21) oas_r8_fieldm=DBLE(lwns_fieldm)
	  if(i==22) oas_r8_fieldm=DBLE(u10s_fieldm)
!	  if(i==23) oas_r8_fieldm=DBLE(u10i_fieldm)
	  if(i==23) oas_r8_fieldm=DBLE(u10s_fieldm1)
	  if(i==24) oas_r8_fieldm=DBLE(v10s_fieldm)
!	  if(i==25) oas_r8_fieldm=DBLE(v10i_fieldm)
	  if(i==25) oas_r8_fieldm=DBLE(v10s_fieldm1)
	  if(i==26) oas_r8_fieldm=DBLE(rai_fieldm)
	  if(i==27) oas_r8_fieldm=DBLE(sno_fieldm)
	  if(i==28) oas_r8_fieldm=DBLE(rqa_fieldm)
	  if(i==29) oas_r8_fieldm=DBLE(clo_fieldm)
	  if(i==30) oas_r8_fieldm=DBLE(slp_fieldm)
	  if(i==31) oas_r8_fieldm=DBLE(amus_fieldm)
	  if(i==32) oas_r8_fieldm=DBLE(amvs_fieldm)
	  if(i==33) oas_r8_fieldm=DBLE(amui_fieldm)
	  if(i==34) oas_r8_fieldm=DBLE(amvi_fieldm)
	  !oas_r8_fieldm=20.*mask
	  !oas_r8_fieldm=DBLE(mype)
	  status = nf90_put_var(il_file_id, il_var_id(i), oas_r8_fieldm, &
		   start = (/ idatastart , jdatastart /), &
		   count = (/ klon, klat /) )
       enddo
		   
!	  status = nf90_inq_varid(il_file_id, csendname(12) , il_var_id(12))
!	  oas_r8_fieldm=DBLE(t2mi_fieldm)
!	  status = nf90_put_var(il_file_id, il_var_id(12), oas_r8_fieldm, &
!		   start = (/ idatastart , jdatastart /), &
!		   count = (/ klon, klat /) )
		   
	  status = nf90_close(il_file_id)
       ENDIF
    ENDDO
    
    endif
    endif
    endif    !use_oasis_netcdf_output



    endif    ! use_oasis_arctic
#endif   ! OASIS3


  end subroutine couple_prism_send


!#endif

!#endif    ! OASIS3

  subroutine set_dqndt(svar_surf,t2m)
    implicit none
    integer::i,j,klon,klat
    real(kind=realkind),intent(in),target::svar_surf(:,:,:)
    real(kind=realkind),intent(in)::t2m(:,:)
!wangsy    real(kind=realkind)::t2m2(klon,klat)
    real(kind=realkind),pointer,dimension(:,:)::u10ms,u10mi,v10ms,v10mi
     real(kind=realkind)::zst2,zst3,zu10
    klon = ubound(t2m,1)
    klat = ubound(t2m,2)
!wangsy    t2m2(:,:)=svar_surf(:,:,2)
    u10ms=>svar_surf(:,:,13)
!    u10me=>svar_surf(:,:,17)
    u10mi=>svar_surf(:,:,18)
    v10ms=>svar_surf(:,:,19)
!    v10ms=>svar_surf(:,:,23)
    v10mi=>svar_surf(:,:,24)

    !calculate no-solar heat flux derivative
    !Taken from NEMO bulk core formula
    do i=1,klon
       do j=1,klat
          zst2=t2m(i,j)*t2m(i,j)
          zst3=zst2*t2m(i,j)
          zu10=sqrt(u10mi(i,j)*u10mi(i,j)+v10mi(i,j)*v10mi(i,j))
!wangsy          zu10=sqrt(u10ms(i,j)*u10ms(i,j)+v10ms(i,j)*v10ms(i,j))
          dqndt_fieldm(i,j)=-4.0_realkind*0.95_realkind*5.67e-8_realkind*zst3-1.22_realkind* &
               1000.5_realkind*1.63e-3_realkind*zu10  &
               +2.83e6_realkind*1.63e-3_realkind*11637800_realkind*zu10/zst2*exp(-5897.8_realkind/t2m(i,j))

       enddo
    enddo

  end subroutine set_dqndt


end module coupling
