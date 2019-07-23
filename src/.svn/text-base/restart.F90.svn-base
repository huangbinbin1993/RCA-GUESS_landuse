module restart
  use timetype
  use calendar
  use gcm
  use decomp
  use comhkp

  implicit none
  private

  logical,public,save:: doRestart=.false.
  type(time),public,save::initialTime
  type(time),public,save::startTime
  type(time),public,save::stopTime
  type(time),public,save::dumpTime
  type(deltatime),save::dumpInterval
  logical,save::monthly=.false. !else yearly
  public read_restart,writeDump,readDump
  interface myWrite !polymorphism
     module procedure myWrite2 
     module procedure myWrite3
     module procedure myWrite4
     module procedure myWrite1I
  end interface
  interface myRead !polymorphism
     module procedure myRead2 
     module procedure myRead3
     module procedure myRead4
     module procedure myRead1I
  end interface
  interface myReadSPECIAL
     module procedure myReadSPECIAL2,myReadSPECIAL3
  end interface
  interface myWriteSPECIAL
     module procedure myWriteSPECIAL2,myWriteSPECIAL3
  end interface

contains

  subroutine read_restart()

    implicit none
    integer::reyear=1990
    integer::remonth=1
    integer::reday=1
    integer::rehour=0
    integer::remin=0
    integer::resec=0

    integer::stopYear=1990
    integer::stopMonth=1
    integer::stopDay=1
    integer::stopHour=0
    integer::stopMin=0
    integer::stopSec=0
    integer:: ntimesteps=-1  !number of timesteps to take before dump

    namelist /namrestart/ntimesteps,&
         stopYear,stopMonth,stopDay,stopHour,stopMin,stopSec,&
         reyear,  remonth,  reday,  rehour,  remin  ,resec  , &
         doRestart,monthly


    open(57,file='namelists.dat',status='old')
    read(57,nml=namrestart)
    close(57)
    if(mype==0)then
       write(6,nml=namrestart)
       open(1,file='printed_namelists.dat',status='old',position='append')
       write(1,nml=namrestart)
       close(1)
    endif

    startTime = maketime(reyear,  remonth,  reday,  rehour,  remin  ,resec)
    stopTime = maketime(stopYear,stopMonth,stopDay,stopHour,stopMin,stopSec)

    if(.not.doRestart)then
       initialTime = startTime !else it is read from dumpFile
    endif

    if(monthly)then
       dumpInterval = makedeltatime(0,1,0,0,0,0)
    else
       dumpInterval = makedeltatime(1,0,0,0,0,0)
    endif
    if(ntimesteps/=-1)then !if user defined to dump after ntimesteps
       dumpInterval = makedeltatime(0,0,0,0,0,ntimesteps*Nseconds(ndtime))
    endif


    if(monthly)then
       dumpTime = maketime(startTime%year,startTime%month+1,1,0,0,0)
    elseif(ntimesteps==-1)then !yearly
       dumpTime = maketime(startTime%year+1,1,1,0,0,0)
    else !userdefined Ntimesteps
       dumpTime = startTime + dumpInterval
    endif

  end subroutine read_restart


  subroutine writeDump(ktime, &
       apsm,alnpsm,tm,um,vm,qm,sm,edotm,svarm, &
       apsz,alnpsz,tz,uz,vz,qz,sz,edotz,svarz, &
       apsp,alnpsp,tp,up,vp,qp,sp,edotp,svarp, &
       dtdtph,dqdtph,dsdtph,drolddt,dsolddt,accsunny,&
       accrunoff,q2d,accrunoffopl,accrunofffor,accrunofflake,&
       accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr,slwr,tlwr,&
       accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri,&
       accq2d,acccw2d,accslwr,accsswr,acctlwr,acctswr,accsnow,&
       sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn,slwrsea,sswrsea,slwrice,sswrice,&
       svarsm,svarsz,svarsp,svar_surf,acum_surf,account_surf,mean_surf,&
       slwrs,sswrs,eco,tdm,tsc,tss,swm,sdm,swc,snm,snc,rou,roc,alb,fri,frf, &
       totcov,cucov,cov2d,cwpath,pblh,dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,&
       tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
       lcounttice,oldcov,zvarcu,cwcu,raincv_kf,snowcv_kf,umfb,rad_cloud,ci2d,accci2d,&
       cwice,cwliq,qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d,&
       u10_reg,v10_reg,u10s_reg,v10s_reg,u10i_reg,v10i_reg,&
       ice_thic,extrem_surf,sacdg,sacdg2,accpr,lastpr,&
       prog_lakes,tend_lakes,&
       accsnol,accsnoc,utot10ms)

    implicit none
    type(time),intent(in)::ktime

    integer,intent(in),dimension(:)::account_surf

    real(kind=realkind),intent(in),dimension(:,:)::apsm,alnpsm,apsz,alnpsz,apsp,alnpsp
    real(kind=realkind),intent(in),dimension(:,:)::utot10ms
    real(kind=realkind),intent(in),dimension(:,:)::ice_thic
    real(kind=realkind),intent(in),dimension(:,:)::ci2d,accci2d,u10_reg,v10_reg,accsnol,accsnoc
    real(kind=realkind),intent(in),dimension(:,:)::u10s_reg,v10s_reg,u10i_reg,v10i_reg,lastpr,accpr
    real(kind=realkind),intent(in),dimension(:,:)::drolddt,dsolddt,accsunny,accrunoff 
    real(kind=realkind),intent(in),dimension(:,:)::q2d,accrunoffopl,accrunofffor,accrunofflake
    real(kind=realkind),intent(in),dimension(:,:)::accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr
    real(kind=realkind),intent(in),dimension(:,:)::accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri
    real(kind=realkind),intent(in),dimension(:,:)::slwr,tlwr,accq2d,acccw2d
    real(kind=realkind),intent(in),dimension(:,:)::accslwr,accsswr,acctlwr,acctswr,accsnow
    real(kind=realkind),intent(in),dimension(:,:)::sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn
    real(kind=realkind),intent(in),dimension(:,:)::slwrsea,sswrsea,slwrice,sswrice,slwrs,sswrs
    real(kind=realkind),intent(in),dimension(:,:)::tdm,tsc,tss,swm,sdm
    real(kind=realkind),intent(in),dimension(:,:)::swc,snm,snc,rou,roc,alb,fri,frf
    real(kind=realkind),intent(in),dimension(:,:)::cov2d,cwpath
    real(kind=realkind),intent(in),dimension(:,:)::tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov
    real(kind=realkind),intent(in),dimension(:,:)::pblh,lcounttice
    real(kind=realkind),intent(in),dimension(:,:)::raincv_kf,snowcv_kf,umfb
    real(kind=realkind),intent(in),dimension(:,:)::qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d

    real(kind=realkind),intent(in),dimension(:,:,:)::tm,um,vm,qm,sm,edotm,tz,uz,vz,qz,sz,edotz
    real(kind=realkind),intent(in),dimension(:,:,:)::tp,up,vp,qp,sp,edotp
    real(kind=realkind),intent(in),dimension(:,:,:)::dtdtph,dqdtph
    real(kind=realkind),intent(in),dimension(:,:,:)::dsdtph,extrem_surf,sacdg2
    real(kind=realkind),intent(in),dimension(:,:,:)::svarsm,svarsz,svarsp
    real(kind=realkind),intent(in),dimension(:,:,:)::svar_surf
    real(kind=realkind),intent(in),dimension(:,:,:)::acum_surf
    real(kind=realkind),intent(in),dimension(:,:,:)::mean_surf
    real(kind=realkind),intent(in),dimension(:,:,:)::eco
    real(kind=realkind),intent(in),dimension(:,:,:)::totcov,cucov,rad_cloud
    real(kind=realkind),intent(in),dimension(:,:,:)::dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,oldcov,zvarcu,cwcu
    real(kind=realkind),intent(in),dimension(:,:,:)::cwice,cwliq
    real(kind=realkind),intent(in),dimension(:,:,:,:)::svarm,svarz,svarp,sacdg
    real(kind=realkind),intent(in),dimension(:,:,:,:)::prog_lakes,tend_lakes
    integer::iostd
    character(len=80)::fildump
    integer::root=0

    if(ktime==dumpTime)then

       if(mype==root) then
          fildump='dump_'
          if ( dumpTime%year > 100 ) then
             write(fildump,'(''dump_'',i4.4,3i2.2)') &
                  dumpTime%year,dumpTime%month,dumpTime%day,dumpTime%hour
          else
             write(fildump,'(''dump_'',i4.4,3i2.2)') &
                  dumpTime%year+1900,dumpTime%month,dumpTime%day,&
                  dumpTime%hour
          endif

          print *,'RESTART file write ',fildump
          
          open(95,form='unformatted',file=fildump,status='unknown',iostat=iostd)
          write(95) current_time%year, &
               current_time%month,&
               current_time%day,&   
               current_time%hour,&  
               current_time%min, &
               current_time%sec

          write(95) initialtime%year,&
               initialtime%month,&
               initialtime%day,&   
               initialtime%hour,&  
               initialtime%min, &
               initialtime%sec 
       endif

       dumpTime = dumpTime + dumpInterval


       call mywrite(root,95,apsm)
       call mywrite(root,95,alnpsm)
       call mywrite(root,95,tm)
       call mywrite(root,95,um)
       call mywrite(root,95,vm)
       call mywrite(root,95,qm)
       call mywrite(root,95,sm)
       call mywrite(root,95,edotm)
       call mywrite(root,95,svarm)
       call mywrite(root,95,apsz)
       call mywrite(root,95,alnpsz)
       call mywrite(root,95,tz)
       call mywrite(root,95,uz)
       call mywrite(root,95,vz)
       call mywrite(root,95,qz)
       call mywrite(root,95,sz)
       call mywrite(root,95,edotz)
       call mywrite(root,95,svarz)
       call mywrite(root,95,apsp)
       call mywrite(root,95,alnpsp)
       call mywrite(root,95,tp)
       call mywrite(root,95,up)
       call mywrite(root,95,vp)
       call mywrite(root,95,qp)
       call mywrite(root,95,sp)
       call mywrite(root,95,edotp)
       call mywrite(root,95,svarp)
       call mywrite(root,95,dtdtph)
       call mywrite(root,95,dqdtph)
       call mywrite(root,95,dsdtph)
       call mywrite(root,95,drolddt)
       call mywrite(root,95,dsolddt)
       call mywrite(root,95,accsunny)
       call mywrite(root,95,accrunoff)
       call mywrite(root,95,q2d)
       call mywrite(root,95,accrunoffopl)
       call mywrite(root,95,accrunofffor)
       call mywrite(root,95,accrunofflake)
       call mywrite(root,95,accprl)
       call mywrite(root,95,accprc)
       call mywrite(root,95,sswr)
       call mywrite(root,95,tsmax)
       call mywrite(root,95,tsmin)
       call mywrite(root,95,t2max)
       call mywrite(root,95,t2min)
       call mywrite(root,95,uv10max)
       call mywrite(root,95,cw2d)
       call mywrite(root,95,tswr)
       call mywrite(root,95,slwr)
       call mywrite(root,95,tlwr)
       call mywrite(root,95,accsenfs)
       call mywrite(root,95,accsenfi)
       call mywrite(root,95,acclatfs)
       call mywrite(root,95,acclatfi)
       call mywrite(root,95,accslwrs)
       call mywrite(root,95,accslwri)
       call mywrite(root,95,accsswrs)
       call mywrite(root,95,accsswri)
       call mywrite(root,95,accq2d)
       call mywrite(root,95,acccw2d)
       call mywrite(root,95,accslwr)
       call mywrite(root,95,accsswr)
       call mywrite(root,95,acctlwr)
       call mywrite(root,95,acctswr)
       call mywrite(root,95,accsnow)
       call mywrite(root,95,sswdn)
       call mywrite(root,95,slwdn)
       call mywrite(root,95,accsswdn)
       call mywrite(root,95,accslwdn)
       call mywrite(root,95,tswdn)
       call mywrite(root,95,acctswdn)
       call mywrite(root,95,slwrsea)
       call mywrite(root,95,sswrsea)
       call mywrite(root,95,slwrice)
       call mywrite(root,95,sswrice)
       call mywrite(root,95,svarsm)
       call mywrite(root,95,svarsz)
       call mywrite(root,95,svarsp)
       call mywrite(root,95,svar_surf)
       call mywrite(root,95,acum_surf)
       call mywrite(root,95,account_surf)
       call mywrite(root,95,mean_surf)
       call mywrite(root,95,slwrs)
       call mywrite(root,95,sswrs)
       call mywrite(root,95,eco(:,:,42))
       call mywrite(root,95,tdm)
       call mywrite(root,95,tsc)
       call mywrite(root,95,tss)
       call mywrite(root,95,swm)
       call mywrite(root,95,sdm)
       call mywrite(root,95,swc)
       call mywrite(root,95,snm)
       call mywrite(root,95,snc)
       call mywrite(root,95,rou)
       call mywrite(root,95,roc)
       call mywrite(root,95,alb)
       call mywrite(root,95,eco(:,:,10))
       call mywrite(root,95,fri)
       call mywrite(root,95,frf)
       call mywrite(root,95,totcov)
       call mywrite(root,95,cucov)
       call mywrite(root,95,cov2d)
       call mywrite(root,95,cwpath)
       call mywrite(root,95,pblh)
       call mywrite(root,95,dtdt_kf)
       call mywrite(root,95,dqdt_kf)
       call mywrite(root,95,dqcdt_kf)
       call mywrite(root,95,om_t1)
       call mywrite(root,95,tcov)
       call mywrite(root,95,lcov)
       call mywrite(root,95,mcov)
       call mywrite(root,95,hcov)
       call mywrite(root,95,acctcov)
       call mywrite(root,95,acclcov)
       call mywrite(root,95,accmcov)
       call mywrite(root,95,acchcov)
       call mywrite(root,95,lcounttice)
       call mywrite(root,95,om_t2)
       call mywrite(root,95,om_t3)
       call mywrite(root,95,oldcov)
       call mywrite(root,95,zvarcu)
       call mywrite(root,95,cwcu)
       call mywrite(root,95,raincv_kf)
       call mywrite(root,95,snowcv_kf)
       call mywrite(root,95,umfb)
       call mywrite(root,95,rad_cloud)
       call mywrite(root,95,ci2d)
       call mywrite(root,95,accci2d)
       call mywrite(root,95,cwice)
       call mywrite(root,95,cwliq)
       call mywrite(root,95,qu2d)
       call mywrite(root,95,accqu2d)
       call mywrite(root,95,qv2d)
       call mywrite(root,95,accqv2d)
       call mywrite(root,95,vimfc2d)
       call mywrite(root,95,accvimfc2d)
       call mywrite(root,95,u10_reg)
       call mywrite(root,95,v10_reg)
       call mywrite(root,95,u10s_reg)
       call mywrite(root,95,v10s_reg)
       call mywrite(root,95,u10i_reg)
       call mywrite(root,95,v10i_reg)
       call mywrite(root,95,ice_thic)
       call mywrite(root,95,tsmax)
       call mywrite(root,95,tsmin)
       call mywrite(root,95,t2max)
       call mywrite(root,95,t2min)
       call mywrite(root,95,uv10max)
       call mywrite(root,95,extrem_surf)
       call mywrite(root,95,sacdg)
       call mywrite(root,95,sacdg2)
       call mywrite(root,95,accpr)
       call mywrite(root,95,lastpr)
       call mywrite(root,95,accsnol)
       call mywrite(root,95,accsnoc)
       call mywrite(root,95,utot10ms)
       call mywrite(root,95,prog_lakes)
       call mywrite(root,95,tend_lakes)

       if(mype==root)then
          close(unit=95)
       endif
    endif
  end subroutine writeDump

  subroutine readDump(apsm,alnpsm,tm,um,vm,qm,sm,edotm,svarm, &
       apsz,alnpsz,tz,uz,vz,qz,sz,edotz,svarz, &
       apsp,alnpsp,tp,up,vp,qp,sp,edotp,svarp, &
       dtdtph,dqdtph,dsdtph,drolddt,dsolddt,accsunny,&
       accrunoff,q2d,accrunoffopl,accrunofffor,accrunofflake,&
       accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr,slwr,tlwr,&
       accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri,&
       accq2d,acccw2d,accslwr,accsswr,acctlwr,acctswr,accsnow,&
       sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn,slwrsea,sswrsea,slwrice,sswrice,&
       svarsm,svarsz,svarsp,svar_surf,acum_surf,account_surf,mean_surf,&
       slwrs,sswrs,eco,tdm,tsc,tss,swm,sdm,swc,snm,snc,rou,roc,alb,fri,frf, &
       totcov,cucov,cov2d,cwpath,pblh,dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,&
       tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov,&
       lcounttice,oldcov,zvarcu,cwcu,raincv_kf,snowcv_kf,umfb,rad_cloud,ci2d,accci2d,&
       cwice,cwliq,qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d,&
       u10_reg,v10_reg,u10s_reg,v10s_reg,u10i_reg,v10i_reg,&
       ice_thic,extrem_surf,sacdg,sacdg2,accpr,lastpr,&
       prog_lakes,tend_lakes,&
       accsnol,accsnoc,utot10ms)

    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer::ierr
    real(kind=realkind),intent(out),dimension(:,:)::apsm,alnpsm,apsz,alnpsz,apsp,alnpsp
    real(kind=realkind),intent(out),dimension(:,:,:)::tm,um,vm,qm,sm,edotm,tz,uz,vz,qz,sz,edotz,tp,up,vp,qp,sp,edotp
    real(kind=realkind),intent(out),dimension(:,:,:,:)::svarm,svarz,svarp,sacdg
    real(kind=realkind),intent(out),dimension(:,:)::utot10ms
    real(kind=realkind),intent(out),dimension(:,:)::ice_thic
    real(kind=realkind),intent(out),dimension(:,:)::ci2d,accci2d,u10_reg,v10_reg,accsnol,accsnoc
    real(kind=realkind),intent(out),dimension(:,:)::u10s_reg,v10s_reg,u10i_reg,v10i_reg,lastpr,accpr
    real(kind=realkind),intent(out),dimension(:,:,:)::dtdtph,dqdtph
    real(kind=realkind),intent(out),dimension(:,:,:)::dsdtph,extrem_surf,sacdg2
    real(kind=realkind),intent(out),dimension(:,:)::drolddt,dsolddt,accsunny,accrunoff 
    real(kind=realkind),intent(out),dimension(:,:)::q2d,accrunoffopl,accrunofffor,accrunofflake
    real(kind=realkind),intent(out),dimension(:,:)::accprl,accprc,sswr,tsmax,tsmin,t2max,t2min,uv10max,cw2d,tswr
    real(kind=realkind),intent(out),dimension(:,:)::accsenfs,accsenfi,acclatfs,acclatfi,accslwrs,accslwri,accsswrs,accsswri
    real(kind=realkind),intent(out),dimension(:,:)::slwr,tlwr,accq2d,acccw2d
    real(kind=realkind),intent(out),dimension(:,:)::accslwr,accsswr,acctlwr,acctswr,accsnow
    real(kind=realkind),intent(out),dimension(:,:)::sswdn,slwdn,accsswdn,accslwdn,tswdn,acctswdn
    real(kind=realkind),intent(out),dimension(:,:)::slwrsea,sswrsea,slwrice,sswrice,slwrs,sswrs
    real(kind=realkind),intent(out),dimension(:,:,:)::svarsm,svarsz,svarsp
    integer,intent(out)::account_surf(:)
    real(kind=realkind),intent(out),dimension(:,:,:)::svar_surf
    real(kind=realkind),intent(out),dimension(:,:,:)::acum_surf
    real(kind=realkind),intent(out),dimension(:,:,:)::mean_surf
    real(kind=realkind),intent(out):: eco(:,:,:)
    real(kind=realkind),intent(out),dimension(:,:)::tdm,tsc,tss,swm,sdm
    real(kind=realkind),intent(out),dimension(:,:)::swc,snm,snc,rou,roc,alb,fri,frf
    real(kind=realkind),intent(out),dimension(:,:,:)::totcov,cucov,rad_cloud
    real(kind=realkind),intent(out),dimension(:,:)::cov2d,cwpath
    real(kind=realkind),intent(out),dimension(:,:)::tcov,lcov,mcov,hcov,acctcov,acclcov,accmcov,acchcov
    real(kind=realkind),intent(out),dimension(:,:)::pblh,lcounttice
    real(kind=realkind),intent(out),dimension(:,:,:)::dtdt_kf,dqdt_kf,dqcdt_kf,om_t1,om_t2,om_t3,oldcov,zvarcu,cwcu
    real(kind=realkind),intent(out),dimension(:,:)::raincv_kf,snowcv_kf,umfb
    real(kind=realkind),intent(out),dimension(:,:,:)::cwice,cwliq
    real(kind=realkind),intent(out),dimension(:,:)::qu2d,accqu2d,qv2d,accqv2d,vimfc2d,accvimfc2d
    real(kind=realkind),intent(out),dimension(:,:,:,:)::prog_lakes,tend_lakes
    integer::iostd,ibuf(12)
    character(len=80)::fildump
    integer::root=0,error

    if(mype==root) then
       fildump='dump_'
       if ( startTime%year > 100 ) then
          write(fildump,'(''dump_'',i4.4,3i2.2)') &
               startTime%year,startTime%month,startTime%day,startTime%hour
       else
          write(fildump,'(''dump_'',i4.4,3i2.2)') &
               startTime%year+1900,startTime%month,startTime%day,&
               startTime%hour
       endif
       print *,'READ RESTART FILE',fildump
       open(95,form='unformatted',file=fildump,status='unknown',iostat=iostd)
    endif
#ifdef MPI_SRC
    call MPI_BCAST(iostd,1,MPI_INTEGER,root,localComm,error)
#endif
    if(iostd/=0)then
       stop 'Could not open restart file!'
    endif
    if(mype==root) then
       read(95)current_time%year, &
            current_time%month,&   
            current_time%day,&  
            current_time%hour, &
            current_time%min,&
            current_time%sec
       read(95) initialtime%year,&
            initialtime%month,&   
            initialtime%day,&  
            initialtime%hour, &
            initialtime%min,&
            initialtime%sec 
    endif

    ibuf(1:6) = (/current_time%year, current_time%month, current_time%day, current_time%hour, current_time%min,current_time%sec/)
    ibuf(7:12) = (/initialtime%year, initialtime%month, initialtime%day, initialtime%hour,  initialtime%min,   initialtime%sec/)
#ifdef MPI_SRC
    call MPI_Bcast(ibuf,12,MPI_INTEGER,root,localComm,error)
#endif
    current_time%year = ibuf(1)
    current_time%month= ibuf(2)
    current_time%day  = ibuf(3)
    current_time%hour = ibuf(4)
    current_time%min  = ibuf(5)
    current_time%sec  = ibuf(6)
    initialtime%year = ibuf(7)
    initialtime%month= ibuf(8)
    initialtime%day  = ibuf(9)
    initialtime%hour = ibuf(10)
    initialtime%min  = ibuf(11)
    initialtime%sec  = ibuf(12)

    call myread(root,95,apsm)
    call myread(root,95,alnpsm)
    call myread(root,95,tm)
    call myread(root,95,um)
    call myread(root,95,vm)
    call myread(root,95,qm)
    call myread(root,95,sm)
    call myread(root,95,edotm)
    call myread(root,95,svarm)
    call myread(root,95,apsz)
    call myread(root,95,alnpsz)
    call myread(root,95,tz)
    call myread(root,95,uz)
    call myread(root,95,vz)
    call myread(root,95,qz)
    call myread(root,95,sz)
    call myread(root,95,edotz)
    call myread(root,95,svarz)
    call myread(root,95,apsp)
    call myread(root,95,alnpsp)
    call myread(root,95,tp)
    call myread(root,95,up)
    call myread(root,95,vp)
    call myread(root,95,qp)
    call myread(root,95,sp)
    call myread(root,95,edotp)
    call myread(root,95,svarp)
    call myread(root,95,dtdtph)
    call myread(root,95,dqdtph)
    call myread(root,95,dsdtph)
    call myread(root,95,drolddt)
    call myread(root,95,dsolddt)
    call myread(root,95,accsunny)
    call myread(root,95,accrunoff)
    call myread(root,95,q2d)
    call myread(root,95,accrunoffopl)
    call myread(root,95,accrunofffor)
    call myread(root,95,accrunofflake)
    call myread(root,95,accprl)
    call myread(root,95,accprc)
    call myread(root,95,sswr)
    call myread(root,95,tsmax)
    call myread(root,95,tsmin)
    call myread(root,95,t2max)
    call myread(root,95,t2min)
    call myread(root,95,uv10max)
    call myread(root,95,cw2d)
    call myread(root,95,tswr)
    call myread(root,95,slwr)
    call myread(root,95,tlwr)
    call myread(root,95,accsenfs)
    call myread(root,95,accsenfi)
    call myread(root,95,acclatfs)
    call myread(root,95,acclatfi)
    call myread(root,95,accslwrs)
    call myread(root,95,accslwri)
    call myread(root,95,accsswrs)
    call myread(root,95,accsswri)
    call myread(root,95,accq2d)
    call myread(root,95,acccw2d)
    call myread(root,95,accslwr)
    call myread(root,95,accsswr)
    call myread(root,95,acctlwr)
    call myread(root,95,acctswr)
    call myread(root,95,accsnow)
    call myread(root,95,sswdn)
    call myread(root,95,slwdn)
    call myread(root,95,accsswdn)
    call myread(root,95,accslwdn)
    call myread(root,95,tswdn)
    call myread(root,95,acctswdn)
    call myread(root,95,slwrsea)
    call myread(root,95,sswrsea)
    call myread(root,95,slwrice)
    call myread(root,95,sswrice)
    call myread(root,95,svarsm)
    call myread(root,95,svarsz)
    call myread(root,95,svarsp)
    call myread(root,95,svar_surf)
    call myread(root,95,acum_surf)
    call myread(root,95,account_surf)
    call myread(root,95,mean_surf)
    call myread(root,95,slwrs)
    call myread(root,95,sswrs)
    call myread(root,95,eco(:,:,42))
    call myread(root,95,tdm)
    call myread(root,95,tsc)
    call myread(root,95,tss)
    call myread(root,95,swm)
    call myread(root,95,sdm)
    call myread(root,95,swc)
    call myread(root,95,snm)
    call myread(root,95,snc)
    call myread(root,95,rou)
    call myread(root,95,roc)
    call myread(root,95,alb)
    call myread(root,95,eco(:,:,10))
    call myread(root,95,fri)
    call myread(root,95,frf)
    call myread(root,95,totcov)
    call myread(root,95,cucov)
    call myread(root,95,cov2d)
    call myread(root,95,cwpath)
    call myread(root,95,pblh)
    call myread(root,95,dtdt_kf)
    call myread(root,95,dqdt_kf)
    call myread(root,95,dqcdt_kf)
    call myread(root,95,om_t1)
    call myread(root,95,tcov)
    call myread(root,95,lcov)
    call myread(root,95,mcov)
    call myread(root,95,hcov)
    call myread(root,95,acctcov)
    call myread(root,95,acclcov)
    call myread(root,95,accmcov)
    call myread(root,95,acchcov)
    call myread(root,95,lcounttice)
    call myread(root,95,om_t2)
    call myread(root,95,om_t3)
    call myread(root,95,oldcov)
    call myread(root,95,zvarcu)
    call myread(root,95,cwcu)
    call myread(root,95,raincv_kf)
    call myread(root,95,snowcv_kf)
    call myread(root,95,umfb)
    call myread(root,95,rad_cloud)
    call myread(root,95,ci2d)
    call myread(root,95,accci2d)
    call myread(root,95,cwice)
    call myread(root,95,cwliq)
    call myread(root,95,qu2d)
    call myread(root,95,accqu2d)
    call myread(root,95,qv2d)
    call myread(root,95,accqv2d)
    call myread(root,95,vimfc2d)
    call myread(root,95,accvimfc2d)
    call myread(root,95,u10_reg)
    call myread(root,95,v10_reg)
    call myread(root,95,u10s_reg)
    call myread(root,95,v10s_reg)
    call myread(root,95,u10i_reg)
    call myread(root,95,v10i_reg)
    call myread(root,95,ice_thic)
    call myread(root,95,tsmax)
    call myread(root,95,tsmin)
    call myread(root,95,t2max)
    call myread(root,95,t2min)
    call myread(root,95,uv10max)
    call myread(root,95,extrem_surf)
    call myread(root,95,sacdg)
    call myread(root,95,sacdg2)
    call myread(root,95,accpr)
    call myread(root,95,lastpr)
    call myread(root,95,accsnol)
    call myread(root,95,accsnoc)
    call myread(root,95,utot10ms)
    call myread(root,95,prog_lakes)
    call myread(root,95,tend_lakes)

    if(mype==root)then
       close(unit=95)
    endif

  end subroutine readDump


  subroutine myWriteSPECIAL2(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:)

    if(mype==root)then
       write(unit)p
    endif

  end subroutine myWriteSPECIAL2
 subroutine myWriteSPECIAL3(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:,:)
    if(mype==root)then
       write(unit)p
    endif
  end subroutine myWriteSPECIAL3
  subroutine myReadSPECIAL2(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(inout)::p(:,:)
    integer::ni,nj,nsize,error
    if(mype==root)then
       read(unit)p
    endif
#ifdef MPI_SRC
    ni = ubound(p,1)
    nj = ubound(p,2)
    nsize = ni
    if(nj>1)nsize = nsize*nj
    call MPI_Bcast(p,nsize,realtype,root,localComm,error)
#endif
  end subroutine myReadSPECIAL2

  subroutine myReadSPECIAL3(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(inout)::p(:,:,:)
    integer::ni,nj,nk,nsize,error
    if(mype==root)then
       read(unit)p
    endif
#ifdef MPI_SRC
    ni = ubound(p,1)
    nj = ubound(p,2)
    nk = ubound(p,3)
    nsize = ni
    if(nj>1)nsize = nsize*nj
    if(nk>1)nsize = nsize*nk

    call MPI_Bcast(p,nsize,realtype,root,localComm,error)
#endif
  end subroutine myReadSPECIAL3

  subroutine myWrite2(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat
    allocate(p_global(klon_global,klat_global))
    klon = ubound(p,1)
    klat = ubound(p,2)
    call colfld(root, p_global, p, klon, klat)  

    if(mype==root)then
       write(unit)p_global
    endif
    deallocate(p_global)
  end subroutine myWrite2

  subroutine myWrite3(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:)
    integer::klon,klat,klev,k
    klon = ubound(p,1)
    klat = ubound(p,2)
    klev = ubound(p,3)
    allocate(p_global(klon_global,klat_global,klev))
    
    do k=1,klev
       call colfld(root, p_global(:,:,k), p(:,:,k), klon, klat)  
    enddo
    if(mype==root)then
       write(unit)p_global
    endif
    deallocate(p_global)
  end subroutine myWrite3


  subroutine myWrite4(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:,:)
    integer::klon,klat,klev,ksvar,k,j
    klon = ubound(p,1)
    klat = ubound(p,2)
    klev = ubound(p,3)
    ksvar = ubound(p,4)
    allocate(p_global(klon_global,klat_global,klev,ksvar))
    
    do j=1,ksvar
       do k=1,klev
          call colfld(root, p_global(:,:,k,j), p(:,:,k,j), klon, klat)  
       enddo
    enddo
    if(mype==root)then
       write(unit)p_global
    endif
    deallocate(p_global)

  end subroutine myWrite4

  subroutine myWrite1I(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    integer,intent(in)::p(:)

    if(mype==root)then
       write(unit)p
    endif
  end subroutine myWrite1I





  subroutine myRead2(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(inout)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat,i,j,error
    allocate(p_global(klon_global,klat_global))

    if(mype==root)then
       read(unit)p_global
    endif
#ifdef MPI_SRC
    call MPI_Bcast(p_global,klon_global*klat_global,realtype,root,localComm,error)
#endif    
    klon = ubound(p,1)
    klat = ubound(p,2)
    do j=1,klat
       do i=1,klon
          p(i,j) = p_global(i+idatastart-1,j+jdatastart-1)
       enddo
    enddo

    deallocate(p_global)
  end subroutine myRead2

  subroutine myRead3(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(inout)::p(:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:)
    integer::klon,klat,klev,k,error,i,j
    klon = ubound(p,1)
    klat = ubound(p,2)
    klev = ubound(p,3)
    allocate(p_global(klon_global,klat_global,klev))
    if(mype==root)then
       read(unit)p_global
    endif
#ifdef MPI_SRC
    call MPI_Bcast(p_global,klon_global*klat_global*klev,realtype,root,localComm,error)
#endif  

    
    do k=1,klev
       do j=1,klat
          do i=1,klon
             p(i,j,k) = p_global(i+idatastart-1,j+jdatastart-1,k)
          enddo
       enddo
    enddo
    
    deallocate(p_global)
  end subroutine myRead3


  subroutine myRead4(root,unit,p)
    implicit none
    integer,intent(in)::unit,root
    real(kind=realkind),intent(inout)::p(:,:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:,:)
    integer::klon,klat,klev,ksvar,k,l,i,j,error
    klon = ubound(p,1)
    klat = ubound(p,2)
    klev = ubound(p,3)
    ksvar = ubound(p,4)
    allocate(p_global(klon_global,klat_global,klev,ksvar))
    if(mype==root)then
       read(unit)p_global
    endif
#ifdef MPI_SRC
    call MPI_Bcast(p_global,klon_global*klat_global*klev*ksvar,realtype,root,localComm,error)
#endif  
    do l=1,ksvar
       do k=1,klev
          do j=1,klat
             do i=1,klon
                p(i,j,k,l) = p_global(i+idatastart-1,j+jdatastart-1,k,l)
             enddo
          enddo
       enddo
    enddo
    deallocate(p_global)
  end subroutine myRead4

  subroutine myRead1I(root,unit,p)
    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer,intent(in)::unit,root
    integer,intent(inout)::p(:)
    integer::l,error
    l = ubound(p,1)
    if(mype==root)then
       read(unit)p
    endif
#ifdef MPI_SRC
    call MPI_Bcast(p,l,MPI_INTEGER,root,localComm,error)
#endif      

  end subroutine myRead1I





!!$  subroutine comm_restart()
!!$    implicit none
!!$#ifdef MPI_SRC
!!$#include"mpif.h"
!!$    integer::ierr
!!$    integer::ibuf(10)
!!$    ibuf(1) = nstop1
!!$    ibuf(2) = nstart1
!!$    ibuf(3) = nut1
!!$    ibuf(4) = startTime%year
!!$    ibuf(5) = startTime%month
!!$    ibuf(6) = startTime%day
!!$    ibuf(7) = startTime%hour
!!$    ibuf(8) = startTime%min
!!$    ibuf(9) = startTime%sec
!!$    if(doRestart)then
!!$       ibuf(10)=0
!!$    else
!!$       ibuf(10)=1
!!$    endif
!!$    call MPI_Bcast(ibuf,10,MPI_INTEGER,0,LOCALCOMM,ierr)
!!$    nstop1     =   ibuf(1)
!!$    nstart1    =   ibuf(2)
!!$    nut1       =   ibuf(3)
!!$    startTime%year     =   ibuf(4)
!!$    startTime%month    =   ibuf(5)
!!$    startTime%day     =   ibuf(6)
!!$    startTime%hour    =   ibuf(7)
!!$    startTime%min     =   ibuf(8)
!!$    startTime%sec    =   ibuf(9)
!!$    doRestart = (ibuf(10)==0)
!!$#endif    
!!$  end subroutine comm_restart
end module restart
