module ecoclimap
  use timetype
  use calendar
  use decomp
  use domainmod
  use flake
  implicit none
  private
  real(kind=realkind),allocatable,dimension(:,:,:,:),save::ecoData
  type(time),save::dataTime(2)
  type(deltatime),save::dtdata
  integer::halo=1
  

  public initEcoclimap,setInitialEcoclimap,setEcoclimap
contains

  subroutine initEcoclimap(startTime,klon,klat,meco,RCAdom,lcounttice)
    implicit none
    
    type(domain),intent(in)::RCAdom
    type(time),intent(in)::startTime
    real(kind=realkind),intent(inout)::lcounttice(klon,klat)
    integer,intent(in)::klon,klat,meco
    integer::k,ierr,i,j,c
    external readecoclimapdata
    
    dtdata = makedeltatime(0,1,0,0,0,0)
    dataTime(1) = startTime
    if(startTime%day < nDaysInMonth(startTime%year,startTime%month)/2+1 )then
       dataTime(1)%month = startTime%month - 1
       if(dataTime(1)%month==0)then
          dataTime(1)%year=dataTime(1)%year-1
          dataTime(1)%month=12
       endif
    endif
    dataTime(1)%day = nDaysInMonth(dataTime(1)%year,dataTime(1)%month)/2+1
    dataTime(1)%hour = 0
    dataTime(1)%min = 0
    dataTime(1)%sec = 0
       
  
    dataTime(2) = dataTime(1) + dtdata
    dataTime(2)%day = nDaysInMonth(dataTime(2)%year,dataTime(2)%month)/2+1
    dataTime(2)%hour = 0
    dataTime(2)%min = 0
    dataTime(2)%sec = 0
    
    
    allocate(ecoData(klon,klat,meco,2))

    do k=1,2
       call readecoclimapdata(klon,klat,klon_global,klat_global,halo,  &
            RCAdom%south,RCAdom%west,RCAdom%dlon,RCAdom%dlat,  &
            RCAdom%polon,RCAdom%polat,&
            ecoData(:,:,:,k),dataTime(k)%month,&
            idatastart,jdatastart,localComm)
       where ( ecoData(:,:,47:49,k)>40.0_realkind ) ecoData(:,:,47:49,k)=40.0_realkind 
       call lakemasks(klon, klat, meco,ecoData(:,:,:,k),lcounttice)  

       
       if(k==1)then
          ecoData(:,:,:,2) = ecoData(:,:,:,1)
#ifdef MPI_SRC
          call MPI_BARRIER(localComm,ierr)
#endif          
       endif
    enddo
  

  end subroutine initEcoclimap




  subroutine setInitialEcoclimap(ktime,klon,klat,meco,eco)
    implicit none
    type(time),intent(in)::ktime
    integer,intent(in)::klon,klat,meco
    real(kind=realkind),intent(inout)::eco(klon,klat,meco)
    real(kind=realkind)::timrat,tim,zdt
    integer::i,j,k
    if(ktime>=dataTime(1).and.ktime<=dataTime(2))then
       tim = real(Nseconds(ktime-dataTime(1)),realkind) 
       zdt = real(Nseconds2(dtdata,dataTime(1)),realkind)
       timrat = tim/zdt
       if(timrat<0._realkind.or.timrat>1._realkind)then
          call stop_program( 'timrat OutOfBounds ecoclimapInit')
       endif
       do k=1,meco
          eco(:,:,k) = (1.0_realkind-timrat)*ecoData(:,:,k,1) + timrat*ecoData(:,:,k,2)

          do j=1,klat
             do i=1,klon
                if(eco(i,j,k)<-665.0_realkind)then
                   call stop_program('ecoclimap error read')
                endif
             enddo
          enddo

       enddo
    else
       call stop_program &
            ('Initial data(ecoclimap) cannot be interpolated from data')
    endif
  end subroutine setInitialEcoclimap

  

  subroutine setEcoclimap(ktime,klon,klat,meco,eco,RCAdom,lcounttice)
    implicit none
    type(time),intent(in)::ktime
    type(domain),intent(in)::RCAdom
    integer,intent(in)::klon,klat,meco
    real(kind=realkind),intent(inout)::eco(klon,klat,meco)
    real(kind=realkind),intent(inout)::lcounttice(klon,klat)
    real(kind=realkind)::timrat
    integer::k
    external readecoclimapdata
    
    do while(ktime>dataTime(2))
       ecoData(:,:,:,1) = ecoData(:,:,:,2)
       dataTime(1) = dataTime(2)
       dataTime(2) = dataTime(1) + dtdata

       dataTime(2)%day = nDaysInMonth(dataTime(2)%year,dataTime(2)%month)/2+1
       dataTime(2)%hour = 0
       dataTime(2)%min = 0
       dataTime(2)%sec = 0
       call readecoclimapdata(klon,klat,klon_global,klat_global,halo,  &
            RCAdom%south,RCAdom%west,RCAdom%dlon,RCAdom%dlat,  &
            RCAdom%polon,RCAdom%polat,&
            ecoData(:,:,:,2),dataTime(2)%month,idatastart,jdatastart,localComm)
       where ( ecoData(:,:,47:49,2)>40.0_realkind ) ecoData(:,:,47:49,2)=40.0_realkind 
       call lakemasks(klon, klat, meco,ecoData(:,:,:,2),lcounttice)   
    enddo

        if(ktime>=dataTime(1).and.ktime<=dataTime(2))then
           timrat = real(Nseconds(ktime-dataTime(1)),realkind)/ &
                real(Nseconds2(dtData,dataTime(1)),realkind)
           if(timrat<0.0_realkind.or.timrat>1.0_realkind)then
              print *,timrat
              print *,ktime
              print *,dataTime(1)
              print *,ktime-dataTime(1)
              print *,Nseconds(ktime-dataTime(1)),'seconds'
              print *,Nseconds(dtData),'seconds for dtData'
              call stop_program( 'Error in ecoclimap timrat errr')
           endif
           do k=1,meco
              eco(:,:,k)=(1.0_realkind-timrat)*ecoData(:,:,k,1)+timrat*ecoData(:,:,k,2)
           enddo
           
        endif

  end subroutine setEcoclimap

  
end module ecoclimap
