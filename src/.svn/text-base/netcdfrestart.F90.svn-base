module netcdfrestart
#ifdef USING_NETCDF
  use rcadomainmod
  use timetype
  use calendar
  use gcm
  use decomp
  use comhkp
  use netcdf
  use util
  implicit none
  private

#include"netcdf.inc"

#ifdef MPI_SRC
#include"mpif.h"
#endif

#define COUNT_FIELDS 1
#define HEADER_INFO 2
#define FIELDRW 3

  character(len=*),parameter:: UNITS = "units"
  character(len=*),parameter:: LVL_NAME = "level"
  character(len=*),parameter:: LAT_NAME = "latitude"
  character(len=*),parameter:: LON_NAME = "longitude"
  character(len=*),parameter:: KSVAR_NAME = "component"


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
     module procedure myWrite2t
     module procedure myWrite3t
     module procedure mywriteMean
  end interface
  interface myRead !polymorphism
     module procedure myRead2t
     module procedure myRead3t
     module procedure myreadMean
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
       accsnol,accsnoc,utot10ms,tsclim_years)

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
    real(kind=realkind),intent(in),dimension(:,:)::slwr,tlwr,accq2d,acccw2d,tsclim_years
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

    integer::iostd,field,ind
    integer::id
    character(len=80)::fildump,fieldname

    integer::root=0,fid=95
    integer::operations(3)=(/COUNT_FIELDS,HEADER_INFO,FIELDRW/),operation
    integer::varid(1000),dimid(5)
    logical,parameter::timeDep=.true.
    varid=-1

    if(ktime==dumpTime)then

       if(mype==root) then
          fildump='dump_'
          if ( dumpTime%year > 100 ) then
             write(fildump,'(''dump_'',i4.4,3i2.2,''.nc'')') &
                  dumpTime%year,dumpTime%month,dumpTime%day,dumpTime%hour
          else
             write(fildump,'(''dump_'',i4.4,3i2.2,''.nc'')') &
                  dumpTime%year+1900,dumpTime%month,dumpTime%day,&
                  dumpTime%hour
          endif
          call check(nf90_create(path=fildump, &
               cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid=fid),'')

          call check(nf90_def_dim(fid,"lon",klon_global,dimId(1)),'')
          call check(nf90_def_dim(fid,"lat",klat_global,dimId(2)),'')
          call check(nf90_def_dim(fid,"height",klev_global,dimId(3)),'')
          call check(nf90_def_dim(fid,"height+1",klev_global+1,dimId(4)),'')
          call check(nf90_def_dim(fid,"time",nf90_unlimited,dimId(5)),'')

          call check(nf90_put_att(fid,nf90_global,"current year",current_time%year),'')
          call check(nf90_put_att(fid,nf90_global,"current month",current_time%month),'')
          call check(nf90_put_att(fid,nf90_global,"current day",current_time%day),'')
          call check(nf90_put_att(fid,nf90_global,"current hour",current_time%hour),'')
          call check(nf90_put_att(fid,nf90_global,"current min",current_time%min),'')
          call check(nf90_put_att(fid,nf90_global,"current sec",current_time%sec),'')

          call check(nf90_put_att(fid,nf90_global,"initial year",initialtime%year),'')
          call check(nf90_put_att(fid,nf90_global,"initial month",initialtime%month),'')
          call check(nf90_put_att(fid,nf90_global,"initial day",initialtime%day),'')
          call check(nf90_put_att(fid,nf90_global,"initial hour",initialtime%hour),'')
          call check(nf90_put_att(fid,nf90_global,"initial min",initialtime%min),'')
          call check(nf90_put_att(fid,nf90_global,"initial sec",initialtime%sec),'')
       endif

       dumpTime = dumpTime + dumpInterval



       !VARIABLES
       do operation=operations(1),operations(3)
          id = 1
          call writeHeaderInfo(root,fid,varid,id,operation,dimId) !id is updated inside
          call mywrite(root,fid,lon,"lon",varid(id),operation,dimId)
          id = id + 1
          call mywrite(root,fid,lat,"lat",varid(id),operation,dimId)
          id = id + 1
          call mywrite(root,fid,apsm,"apsm",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,alnpsm,"log(surface pressure) m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tm,"temperature m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,um,"u-velocity m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,vm,"v-velocity m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,qm,"q m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sm,"s m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,edotm,"edot m-1",varid(id),operation,dimId,timeDep)
          id = id + 1
          do field=1,ubound(svarm,4)
             write(fieldname,'(''svar m-1'',i3.3)')field
             call mywrite(root,fid,svarm(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          call mywrite(root,fid,apsz,"surface pressure m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,alnpsz,"log(surface pressure) m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tz,"temperature m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,uz,"u-velocity m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,vz,"v-velocity",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,qz,"q m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sz,"s m",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,edotz,"edot m",varid(id),operation,dimId,timeDep)
          id = id + 1 !klev+1
          do field=1,ubound(svarz,4)
             write(fieldname,'(''svar m'',i3.3)')field
             call mywrite(root,fid,svarz(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          call mywrite(root,fid,apsp,"surface pressure m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,alnpsp,"log(surface pressure) m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tp,"temperature m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,up,"u-velocity m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,vp,"v-velocity m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,qp,"q m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sp,"s m+1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,edotp,"edot m+1",varid(id),operation,dimId,timeDep)
          id = id + 1!klev+1
          do field=1,ubound(svarp,4)
             write(fieldname,'(''svar m+1'',i3.3)')field
             call mywrite(root,fid,svarp(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          call mywrite(root,fid,dtdtph,"dtdtph",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dqdtph,"dqdtph",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dsdtph,"dsdtph",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,drolddt,"drolddt",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dsolddt,"dsolddt",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsunny,"accsunny",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accrunoff,"accrunoff",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,q2d,"q2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accrunoffopl,"accrunoffopl",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accrunofffor,"accrunofffor",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accrunofflake,"accrunofflake",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accprl,"accprl",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accprc,"accprc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sswr,"sswr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cw2d,"cw2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tswr,"tswr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,slwr,"slwr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tlwr,"tlwr",varid(id),operation,dimId,timeDep)
          id = id + 1

          call mywrite(root,fid,accsenfs,"accsenfs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsenfi,"accsenfi",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acclatfs,"acclatfs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acclatfi,"acclatfi",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accslwrs,"accslwrs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accslwri,"accslwri",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsswrs,"accsswrs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsswri,"accsswri",varid(id),operation,dimId,timeDep)
          id = id + 1

          call mywrite(root,fid,accq2d,"accq2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acccw2d,"acccw2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accslwr,"accslwr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsswr,"accsswr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acctlwr,"acctlwr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acctswr,"acctswr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsnow,"accsnow",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sswdn,"sswdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,slwdn,"slwdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsswdn,"accsswdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accslwdn,"accslwdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tswdn,"tswdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acctswdn,"acctswdn",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,slwrsea,"slwrsea",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sswrsea,"sswrsea",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,slwrice,"slwrice",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sswrice,"sswrice",varid(id),operation,dimId,timeDep)
          id = id + 1
          do field=1,ubound(svarsm,3)
             write(fieldname,'(''svarsm'',i3.3)')field
             call mywrite(root,fid,svarsm(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(svarsz,3)
             write(fieldname,'(''svarsz'',i3.3)')field
             call mywrite(root,fid,svarsz(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(svarsp,3)
             write(fieldname,'(''svarsp'',i3.3)')field
             call mywrite(root,fid,svarsp(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(svar_surf,3)
             write(fieldname,'(''svar_surf'',i3.3)')field
             call mywrite(root,fid,svar_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(acum_surf,3)
             write(fieldname,'(''acum_surf'',i3.3)')field
             call mywrite(root,fid,acum_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(mean_surf,3)
             write(fieldname,'(''mean_surf'',i3.3)')field
             call mywrite(root,fid,mean_surf(:,:,field),account_surf(field),fieldname,&
                  varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          call mywrite(root,fid,slwrs,"slwrs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sswrs,"sswrs",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,eco(:,:,42),"eco(:,:,42)",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tdm,"tdm",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tsc,"tsc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tss,"tss",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,swm,"swm",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,sdm,"sdm",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,swc,"swc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,snm,"snm",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,snc,"snc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,rou,"rou",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,roc,"roc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,alb,"alb",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,eco(:,:,10),"eco(:,:,10)",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,fri,"fri",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,frf,"frf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,totcov,"totcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cucov,"cucov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cov2d,"cov2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cwpath,"cwpath",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,pblh,"pblh",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dtdt_kf,"dtdt_kf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dqdt_kf,"dqdt_kf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,dqcdt_kf,"dqcdt_kf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,om_t1,"om_t1",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tcov,"tcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,lcov,"lcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,mcov,"mcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,hcov,"hcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acctcov,"acctcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acclcov,"acclcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accmcov,"accmcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,acchcov,"acchcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,lcounttice,"lcounttice",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,om_t2,"om_t2",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,om_t3,"om_t3",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,oldcov,"oldcov",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,zvarcu,"zvarcu",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cwcu,"cwcu",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,raincv_kf,"raincv_kf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,snowcv_kf,"snowcv_kf",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,umfb,"umfb",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,rad_cloud,"rad_cloud",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,ci2d,"ci2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accci2d,"accci2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cwice,"cwice",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,cwliq,"cwliq",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,qu2d,"qu2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accqu2d,"accqu2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,qv2d,"qv2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accqv2d,"accqv2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,vimfc2d,"vimfc2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accvimfc2d,"accvimfc2d",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,u10_reg,"u10_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,v10_reg,"v10_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,u10s_reg,"u10s_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,v10s_reg,"v10s_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,u10i_reg,"u10i_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,v10i_reg,"v10i_reg",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,ice_thic,"ice_thic",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tsmax,"tsmax",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,tsmin,"tsmin",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,t2max,"t2max",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,t2min,"t2min",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,uv10max,"uv10max",varid(id),operation,dimId,timeDep)
          id = id + 1
          do field=1,ubound(extrem_surf,3)
             write(fieldname,'(''extrem_surf'',i3.3)')field
             call mywrite(root,fid,extrem_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(sacdg,4)
             write(fieldname,'(''sacdg'',i3.3)')field
             call mywrite(root,fid,sacdg(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          do field=1,ubound(sacdg2,3)
             write(fieldname,'(''sacdg2'',i3.3)')field
             call mywrite(root,fid,sacdg2(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
          call mywrite(root,fid,accpr,"accpr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,lastpr,"lastpr",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsnol,"accsnol",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,accsnoc,"accsnoc",varid(id),operation,dimId,timeDep)
          id = id + 1
          call mywrite(root,fid,utot10ms,"utot10ms",varid(id),operation,dimId,timeDep)
          id = id + 1
          do field =1,ubound(prog_lakes,4)
             do ind=1,ubound(prog_lakes,3)
                write(fieldname,'(''prog_lakes'',i3.3,i3.3)')ind,field
                call mywrite(root,fid,prog_lakes(:,:,ind,field),fieldname,varid(id),operation,dimId,timeDep)
                id = id + 1
             enddo
          enddo

          do field =1,ubound(prog_lakes,4)
             do ind=1,ubound(prog_lakes,3)
                write(fieldname,'(''tend_lakes'',i3.3,i3.3)')ind,field
                call mywrite(root,fid,tend_lakes(:,:,ind,field),fieldname,varid(id),operation,dimId,timeDep)
                id = id + 1
             enddo
          enddo
          call mywrite(root,fid,tsclim_years,"tsclim_years",varid(id),operation,dimId,timeDep)
          id = id + 1

          id=id-1
          if(mype==root)then
             if(operation==COUNT_FIELDS)then
                print *,id,'fields are being printed...'
             endif
             if(operation==HEADER_INFO)then
                print *,'definitions are made'
                call check( nf90_enddef(fid),'enddef' )
             endif
             if(operation==FIELDRW) call check( nf90_close(fid) ,'close')
          endif
       enddo !end operation
    endif
  end subroutine writeDump



  subroutine writeHeaderInfo(root,unit,varids,id,operation,dimids)
    implicit none
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation

    integer,intent(inout)::varids(:),id
    integer,intent(in)::unit,root
    real(kind=realkind),allocatable::lon(:),lat(:)
    character(len=80)::timename
    real(kind=realkind)::time(1)=0._realkind
    integer::i

    if(operation==COUNT_FIELDS)then
       id = id + 4
       return
    endif
    if(operation==HEADER_INFO)then
       if(mype==root)then
          call check( nf90_def_var(unit, "rlon", nf90_double, dimids(1), varids(id)),"rlon" )
          call check( nf90_put_att(unit,varids(id),units,"degrees"),"rlon")
          call check( nf90_put_att(unit,varids(id),"standard_name","grid_longitude"),"grid_longitude")
          call check( nf90_put_att(unit,varids(id),"long_name","longitude in rotated pole grid"),"grid_lae")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_def_var(unit, "rlat", nf90_double, dimids(2), varids(id)),"rlat" )
          call check( nf90_put_att(unit,varids(id),units,"degrees"),"rlan")
          call check( nf90_put_att(unit,varids(id),"standard_name","grid_latitude"),"grid_latitude")
          call check( nf90_put_att(unit,varids(id),"long_name","latitude in rotated pole grid"),"grid_latitude")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_def_var(unit, "time", nf90_double, dimids(5), varids(id)),"time" )
          write(timename,'(''hours since '',i4.4,''-'',i2.2,''-'',i2.2,'' '',i2.2,'':'',i2.2,'':'',i2.2)')&
               current_time%year,current_time%month,current_time%day,current_time%hour,&
               current_time%min,current_time%sec
          call check( nf90_put_att(unit,varids(id),units,timename),"hours sice X")
          call check( nf90_put_att(unit,varids(id),"calendar","proleptic_gregorian"),"time")
          call check( nf90_put_att(unit,varids(id),"standard_name","time"),"time")
          call check( nf90_put_att(unit,varids(id),"long_name","time"),"time")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_def_var(unit, "rotated_pole", nf90_char,  varids(id)),"rotated_pole def var" )
          call check( nf90_put_att(unit,varids(id),"grid_mapping_name","rotated_latitude_longitude"),"mapping")
          call check( nf90_put_att(unit,varids(id),"grid_north_pole_latitude",rcaDomain%polat-180.0_realkind),"rotlat")
          call check( nf90_put_att(unit,varids(id),"grid_north_pole_longitude",rcaDomain%polon),"rotlon")
       endif
       id = id + 1
    endif
    if(operation==FIELDRW)then   
       allocate(lon(klon_global),lat(klat_global))
       do i=1,klon_global
          lon(i) = rcaDomain%west + (i-1)*rcaDomain%dlon
       enddo
       do i=1,klat_global
          lat(i) = rcaDomain%south + (i-1)*rcaDomain%dlat
       enddo
       if(mype==root)then
          call check( nf90_put_var(unit,varids(id),lon),"rlon")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_put_var(unit,varids(id),lat),"rlat")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_put_var(unit,varids(id),time),"time")
       endif
       id = id + 1
       if(mype==root)then
          call check( nf90_put_var(unit,varids(id),""),"rot char")
       endif
       id = id + 1
       deallocate(lon,lat)

    endif

    return

  end subroutine writeHeaderInfo




  subroutine myWrite2t(root,unit,p,name,id,operation,dimids,timedep)
    implicit none
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation
    character(len=*),intent(in)::name
    logical,intent(in)::timedep

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    if(operation==HEADER_INFO .and. mype==0)then
       call check( nf90_def_var(unit, name, netcdfreal, (/dimids(1:2), dimids(5)/), id),name )
       call check( nf90_put_att(unit, id, units, "NOT IMPLEMENTED"),name ) 
       call check( nf90_put_att(unit, id, "long_name", name),name ) 
       call check( nf90_put_att(unit, id, "standard_name", name),name ) 
       call check( nf90_put_att(unit, id, "coordinates", "lon lat"),name ) 
       call check( nf90_put_att(unit, id, "grid_mapping", "rotated_pole"),name ) 
    elseif(operation==FIELDRW)then   
       allocate(p_global(klon_global,klat_global))
       klon = ubound(p,1)
       klat = ubound(p,2)
       call colfld(root, p_global, p, klon, klat)  

       if(mype==root)then
          call check( nf90_put_var(unit,id,p_global),name)
       endif
       deallocate(p_global)
    endif
  end subroutine myWrite2t

  subroutine myWrite2(root,unit,p,name,id,operation,dimids)
    implicit none
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation
    character(len=*),intent(in)::name
    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    if(operation==HEADER_INFO .and. mype==0)then
       call check( nf90_def_var(unit, name, netcdfreal, dimids(1:2), id),name )
       if(name=='lon'.or.name=='lat')then
          if(name=='lon')then
             call check( nf90_put_att(unit, id, units, "degrees_east"),name ) 
          else
             call check( nf90_put_att(unit, id, units, "degrees_north"),name ) 
          endif
       else
          call check( nf90_put_att(unit, id, "coordinates", "lon lat"),name ) 
          call check( nf90_put_att(unit, id, units, "NOT IMPLEMENTED"),name ) 
          call check( nf90_put_att(unit, id, "grid_mapping", "rotated_pole"),name ) 
       endif
       call check( nf90_put_att(unit, id, "long_name", name),name ) 
       call check( nf90_put_att(unit, id, "standard_name", name),name ) 

    elseif(operation==FIELDRW)then   
       allocate(p_global(klon_global,klat_global))
       klon = ubound(p,1)
       klat = ubound(p,2)
       call colfld(root, p_global, p, klon, klat)  

       if(mype==root)then
          call check( nf90_put_var(unit,id,p_global),name)
       endif
       deallocate(p_global)
    endif
  end subroutine myWrite2



  subroutine myWrite3t(root,unit,p,name,id,operation,dimids,timedep)
    implicit none
    integer,intent(in)::dimids(:)
    logical,intent(in)::timedep
    integer,intent(in)::operation
    character(len=*),intent(in)::name

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:)
    integer::klon,klat,klev,k
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    klon = ubound(p,1)
    klat = ubound(p,2)
    klev = ubound(p,3)

    if(operation==HEADER_INFO.and. mype==0)then
       if(klev==klev_global)then
          call check( nf90_def_var(unit, name, netcdfreal, (/dimids(1:3),dimids(5)/), id),name )
       elseif(klev==klev_global+1)then
          call check( nf90_def_var(unit, name, netcdfreal, (/dimids(1:2),dimids(4:5)/), id),name )
       endif

       call check( nf90_put_att(unit, id, units, "hPa"),name ) 
       call check( nf90_put_att(unit, id, "long_name", name),name ) 
       call check( nf90_put_att(unit, id, "standard_name", name),name ) 
       call check( nf90_put_att(unit, id, "coordinates", "lon lat"),name ) 
       call check( nf90_put_att(unit, id, "grid_mapping", "rotated_pole"),name )
    elseif(operation==FIELDRW)then
       allocate(p_global(klon_global,klat_global,klev))

       do k=1,klev
          call colfld(root, p_global(:,:,k), p(:,:,k), klon, klat)  
       enddo
       if(mype==root)then
          call check( nf90_put_var(unit,id,p_global),name) 
       endif
       deallocate(p_global)
    endif

  end subroutine myWrite3t

  subroutine mywriteMean(root,unit,p,n,name,id,operation,dimids,timedep)
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation,n
    character(len=*),intent(in)::name
    logical,intent(in)::timedep

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(in)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    if(operation==HEADER_INFO .and. mype==0)then
       call check( nf90_def_var(unit, name, netcdfreal, dimids(1:2), id),name )
       call check( nf90_put_att(unit, id, units, "hPa"),name ) 
       call check( nf90_put_att(unit, id, "long_name", name),name ) 
       call check( nf90_put_att(unit, id, "standard_name", name),name ) 
       call check( nf90_put_att(unit, id, "coordinates", "lon lat"),name ) 
       call check( nf90_put_att(unit, id, "grid_mapping", "rotated_pole"),name )
       call check( nf90_put_att(unit, id, "averaging_steps", n),name ) 
    elseif(operation==FIELDRW)then   
       allocate(p_global(klon_global,klat_global))
       klon = ubound(p,1)
       klat = ubound(p,2)
       call colfld(root, p_global, p, klon, klat)  

       if(mype==root)then
          call check( nf90_put_var(unit,id,p_global),name)
       endif
       deallocate(p_global)
    endif
  end subroutine mywriteMean






  subroutine check(status,msg)
    integer,intent(in)::status
    character(len=*)::msg
    if(status/=nf90_noerr)then
       print *, trim(nf90_strerror(status)),'  ::',msg
       stop 2
    endif
  end subroutine check





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
       accsnol,accsnoc,utot10ms,tsclim_years)

    implicit none
#ifdef MPI_SRC
#include"mpif.h"
#endif
    integer::ierr
    real(kind=realkind),intent(out),dimension(:,:)::apsm,alnpsm,apsz,alnpsz,apsp,alnpsp
    real(kind=realkind),intent(out),dimension(:,:,:)::tm,um,vm,qm,sm,edotm,tz,uz,vz,qz,sz,edotz,tp,up,vp,qp,sp,edotp
    real(kind=realkind),intent(out),dimension(:,:,:,:)::svarm,svarz,svarp,sacdg
    real(kind=realkind),intent(out),dimension(:,:)::utot10ms
    real(kind=realkind),intent(out),dimension(:,:)::ice_thic,tsclim_years
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
    integer::iostd,ibuf(12),field
    integer::id,ind
    character(len=80)::fildump,fieldname
    integer::root=0,error,fid=95
    integer::operations(3)=(/COUNT_FIELDS,HEADER_INFO,FIELDRW/),operation
    integer::varid(1000),dimid(5)
    logical,parameter::timeDep=.true.
    varid=-1


    if(mype==root) then
       fildump='dump_'
       if ( startTime%year > 100 ) then
          write(fildump,'(''dump_'',i4.4,3i2.2,''.nc'')') &
               startTime%year,startTime%month,startTime%day,startTime%hour
       else
          write(fildump,'(''dump_'',i4.4,3i2.2,''.nc'')') &
               startTime%year+1900,startTime%month,startTime%day,&
               startTime%hour
       endif
       call check(nf90_open(fildump,NF90_NOWRITE,fid),"nf90open")
       call check(nf90_inq_dimid(fid,"lon",dimid(1)),"inqdimid klon")
       call check(nf90_inq_dimid(fid,"lat",dimid(2)),"inqdimid klat")
       call check(nf90_inq_dimid(fid,"height",dimid(3)),"inqdimid klev")
!!$       if(dimid(1)/=klon_global .or.dimid(2)/=klat_global.or.&
!!$            dimid(3)/=klev_global)then
!!$          print *,'klon_global',dimid(1),klon_global
!!$          print *,'klat_global',dimid(2),klat_global
!!$          print *,'klev_global',dimid(3),klev_global
!!$          stop 'dimensions of restart field do not match with specified namelist'
!!$       endif
       call check(nf90_get_att(fid,nf90_global,"current year",current_time%year),'')
       call check(nf90_get_att(fid,nf90_global,"current month",current_time%month),'')
       call check(nf90_get_att(fid,nf90_global,"current day",current_time%day),'')
       call check(nf90_get_att(fid,nf90_global,"current hour",current_time%hour),'')
       call check(nf90_get_att(fid,nf90_global,"current min",current_time%min),'')
       call check(nf90_get_att(fid,nf90_global,"current sec",current_time%sec),'')

       call check(nf90_get_att(fid,nf90_global,"initial year",initialtime%year),'')
       call check(nf90_get_att(fid,nf90_global,"initial month",initialtime%month),'')
       call check(nf90_get_att(fid,nf90_global,"initial day",initialtime%day),'')
       call check(nf90_get_att(fid,nf90_global,"initial hour",initialtime%hour),'')
       call check(nf90_get_att(fid,nf90_global,"initial min",initialtime%min),'')
       call check(nf90_get_att(fid,nf90_global,"initial sec",initialtime%sec),'')
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

    do operation=operations(1),operations(3)
       id = 1
       call myread(root,fid,apsm,"apsm",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,alnpsm,"log(surface pressure) m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tm,"temperature m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,um,"u-velocity m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,vm,"v-velocity m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,qm,"q m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sm,"s m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,edotm,"edot m-1",varid(id),operation,dimId,timeDep)
       id = id + 1
       do field=1,ubound(svarm,4)
          write(fieldname,'(''svar m-1'',i3.3)')field
          call myread(root,fid,svarm(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       call myread(root,fid,apsz,"surface pressure m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,alnpsz,"log(surface pressure) m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tz,"temperature m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,uz,"u-velocity m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,vz,"v-velocity",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,qz,"q m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sz,"s m",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,edotz,"edot m",varid(id),operation,dimId,timeDep)
       id = id + 1 !klev+1
       do field=1,ubound(svarz,4)
          write(fieldname,'(''svar m'',i3.3)')field
          call myread(root,fid,svarz(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       call myread(root,fid,apsp,"surface pressure m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,alnpsp,"log(surface pressure) m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tp,"temperature m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,up,"u-velocity m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,vp,"v-velocity m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,qp,"q m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sp,"s m+1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,edotp,"edot m+1",varid(id),operation,dimId,timeDep)
       id = id + 1!klev+1
       do field=1,ubound(svarp,4)
          write(fieldname,'(''svar m+1'',i3.3)')field
          call myread(root,fid,svarp(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       call myread(root,fid,dtdtph,"dtdtph",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dqdtph,"dqdtph",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dsdtph,"dsdtph",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,drolddt,"drolddt",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dsolddt,"dsolddt",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsunny,"accsunny",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accrunoff,"accrunoff",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,q2d,"q2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accrunoffopl,"accrunoffopl",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accrunofffor,"accrunofffor",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accrunofflake,"accrunofflake",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accprl,"accprl",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accprc,"accprc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sswr,"sswr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cw2d,"cw2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tswr,"tswr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,slwr,"slwr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tlwr,"tlwr",varid(id),operation,dimId,timeDep)
       id = id + 1

       call myread(root,fid,accsenfs,"accsenfs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsenfi,"accsenfi",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acclatfs,"acclatfs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acclatfi,"acclatfi",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accslwrs,"accslwrs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accslwri,"accslwri",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsswrs,"accsswrs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsswri,"accsswri",varid(id),operation,dimId,timeDep)
       id = id + 1

       call myread(root,fid,accq2d,"accq2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acccw2d,"acccw2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accslwr,"accslwr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsswr,"accsswr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acctlwr,"acctlwr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acctswr,"acctswr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsnow,"accsnow",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sswdn,"sswdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,slwdn,"slwdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsswdn,"accsswdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accslwdn,"accslwdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tswdn,"tswdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acctswdn,"acctswdn",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,slwrsea,"slwrsea",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sswrsea,"sswrsea",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,slwrice,"slwrice",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sswrice,"sswrice",varid(id),operation,dimId,timeDep)
       id = id + 1
       do field=1,ubound(svarsm,3)
          write(fieldname,'(''svarsm'',i3.3)')field
          call myread(root,fid,svarsm(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(svarsz,3)
          write(fieldname,'(''svarsz'',i3.3)')field
          call myread(root,fid,svarsz(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(svarsp,3)
          write(fieldname,'(''svarsp'',i3.3)')field
          call myread(root,fid,svarsp(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(svar_surf,3)
          write(fieldname,'(''svar_surf'',i3.3)')field
          call myread(root,fid,svar_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(acum_surf,3)
          write(fieldname,'(''acum_surf'',i3.3)')field
          call myread(root,fid,acum_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(mean_surf,3)
          write(fieldname,'(''mean_surf'',i3.3)')field
          call myread(root,fid,mean_surf(:,:,field),account_surf(field),fieldname,&
               varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       call myread(root,fid,slwrs,"slwrs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sswrs,"sswrs",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,eco(:,:,42),"eco(:,:,42)",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tdm,"tdm",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tsc,"tsc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tss,"tss",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,swm,"swm",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,sdm,"sdm",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,swc,"swc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,snm,"snm",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,snc,"snc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,rou,"rou",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,roc,"roc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,alb,"alb",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,eco(:,:,10),"eco(:,:,10)",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,fri,"fri",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,frf,"frf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,totcov,"totcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cucov,"cucov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cov2d,"cov2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cwpath,"cwpath",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,pblh,"pblh",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dtdt_kf,"dtdt_kf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dqdt_kf,"dqdt_kf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,dqcdt_kf,"dqcdt_kf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,om_t1,"om_t1",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tcov,"tcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,lcov,"lcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,mcov,"mcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,hcov,"hcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acctcov,"acctcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acclcov,"acclcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accmcov,"accmcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,acchcov,"acchcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,lcounttice,"lcounttice",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,om_t2,"om_t2",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,om_t3,"om_t3",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,oldcov,"oldcov",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,zvarcu,"zvarcu",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cwcu,"cwcu",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,raincv_kf,"raincv_kf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,snowcv_kf,"snowcv_kf",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,umfb,"umfb",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,rad_cloud,"rad_cloud",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,ci2d,"ci2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accci2d,"accci2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cwice,"cwice",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,cwliq,"cwliq",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,qu2d,"qu2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accqu2d,"accqu2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,qv2d,"qv2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accqv2d,"accqv2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,vimfc2d,"vimfc2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accvimfc2d,"accvimfc2d",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,u10_reg,"u10_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,v10_reg,"v10_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,u10s_reg,"u10s_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,v10s_reg,"v10s_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,u10i_reg,"u10i_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,v10i_reg,"v10i_reg",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,ice_thic,"ice_thic",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tsmax,"tsmax",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,tsmin,"tsmin",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,t2max,"t2max",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,t2min,"t2min",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,uv10max,"uv10max",varid(id),operation,dimId,timeDep)
       id = id + 1
       do field=1,ubound(extrem_surf,3)
          write(fieldname,'(''extrem_surf'',i3.3)')field
          call myread(root,fid,extrem_surf(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(sacdg,4)
          write(fieldname,'(''sacdg'',i3.3)')field
          call myread(root,fid,sacdg(:,:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       do field=1,ubound(sacdg2,3)
          write(fieldname,'(''sacdg2'',i3.3)')field
          call myread(root,fid,sacdg2(:,:,field),fieldname,varid(id),operation,dimId,timeDep)
          id = id + 1
       enddo
       call myread(root,fid,accpr,"accpr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,lastpr,"lastpr",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsnol,"accsnol",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,accsnoc,"accsnoc",varid(id),operation,dimId,timeDep)
       id = id + 1
       call myread(root,fid,utot10ms,"utot10ms",varid(id),operation,dimId,timeDep)
       id = id + 1
       do field =1,ubound(prog_lakes,4)
          do ind=1,ubound(prog_lakes,3)
             write(fieldname,'(''prog_lakes'',i3.3,i3.3)')ind,field
             call myread(root,fid,prog_lakes(:,:,ind,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
       enddo

       do field =1,ubound(prog_lakes,4)
          do ind=1,ubound(prog_lakes,3)
             write(fieldname,'(''tend_lakes'',i3.3,i3.3)')ind,field
             call myread(root,fid,tend_lakes(:,:,ind,field),fieldname,varid(id),operation,dimId,timeDep)
             id = id + 1
          enddo
       enddo
       call myread(root,fid,tsclim_years,"tsclim_years",varid(id),operation,dimId,timeDep)
       id = id + 1


       id=id-1
       if(mype==root)then
          if(operation==COUNT_FIELDS)then
             print *,id,'fields are being read...'
          endif
          if(operation==FIELDRW) call check( nf90_close(fid) ,'close')
       endif
    enddo !end operation

  end subroutine readDump


  subroutine myRead2t(root,unit,p,name,id,operation,dimids,timedep)
    implicit none
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation
    character(len=*),intent(in)::name
    logical,intent(in)::timedep

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(out)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat,error,i,j
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    if(operation==HEADER_INFO .and. mype==0)then

       call check( nf90_inq_varid(unit, name,  id),name )

    elseif(operation==FIELDRW)then   
       allocate(p_global(klon_global,klat_global))

       if(mype==root)then

          call check( nf90_get_var(unit,id,p_global),name)
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
    endif
  end subroutine myRead2t

  subroutine myRead3t(root,unit,p,name,id,operation,dimids,timedep)
    implicit none
    integer,intent(in)::dimids(:)
    logical,intent(in)::timedep
    integer,intent(in)::operation
    character(len=*),intent(in)::name

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(out)::p(:,:,:)
    real(kind=realkind),allocatable::p_global(:,:,:)
    integer::klon,klat,klev,k,i,j,error
    integer::NETCDFREAL

    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return


    if(operation==HEADER_INFO.and. mype==0)then
       call check( nf90_inq_varid(unit, name,  id),name )
    elseif(operation==FIELDRW)then
       klon = ubound(p,1)
       klat = ubound(p,2)
       klev = ubound(p,3)
       allocate(p_global(klon_global,klat_global,klev))

       if(mype==root)then
          call check( nf90_get_var(unit,id,p_global),name) 
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
    endif

  end subroutine myRead3t

  subroutine myreadMean(root,unit,p,n,name,id,operation,dimids,timedep)
    integer,intent(in)::dimids(:)
    integer,intent(in)::operation
    integer,intent(out)::n
    character(len=*),intent(in)::name
    logical,intent(in)::timedep

    integer,intent(inout)::id
    integer,intent(in)::unit,root
    real(kind=realkind),intent(out)::p(:,:)
    real(kind=realkind),allocatable::p_global(:,:)
    integer::klon,klat,i,j,error
    integer::NETCDFREAL


    if(realkind==8)then
       NETCDFREAL = nf90_double
    else
       NETCDFREAL = nf90_real
    endif

    if(operation==COUNT_FIELDS)return

    if(operation==HEADER_INFO .and. mype==0)then
       call check( nf90_inq_varid(unit, name, id),name )
       call check( nf90_get_att(unit, id, "averaging_steps", n),name ) 
    elseif(operation==FIELDRW)then   
       allocate(p_global(klon_global,klat_global))
       klon = ubound(p,1)
       klat = ubound(p,2)

       if(mype==root)then
          call check( nf90_get_var(unit,id,p_global),name)
       endif
#ifdef MPI_SRC
       call MPI_Bcast(p_global,klon_global*klat_global,realtype,root,localComm,error)
       call MPI_Bcast(n,1,mpi_integer,root,localComm,error)
#endif 
       do j=1,klat
          do i=1,klon
             p(i,j) = p_global(i+idatastart-1,j+jdatastart-1)
          enddo
       enddo
       deallocate(p_global)
    endif
  end subroutine myreadMean

#endif
end module netcdfrestart
