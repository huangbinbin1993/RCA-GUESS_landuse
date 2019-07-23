module akfcum_mod
  use decomp,only:stop_program,realkind
  implicit none
  private
  public akfcum
contains

  subroutine akfcum(nhor,nlev,kstart,kstop,nca, &
       dtime,zdx_kf,zdy_kf,orogsigm,&
       pf,t,q,u,v,omega_kf,&
       dtdt,dqdt,dqcdt,&
       raincv,snowcv,kf_base,kf_top, &
       umf,om_cpbl,&
       qvap,qliq,qice,clddep,&
       shal_cgj, &
       ncldck,ahyb,bhyb,ps,ts,tke,div_kf,kf_ind,&
       udr,zfrpbl,zfrtop,zlevtop,dz_half, &
       lsb,wkf1,wkf2,wkf3)
    use confys  
    implicit none  

    integer :: nhor,nlev,kf_ind(nhor),kf_base(nhor),kf_top(nhor)
    real(kind=realkind) :: div_kf(nhor,nlev),omega_kf(nhor,nlev)  
    real(kind=realkind) :: ahyb(nlev + 1),bhyb(nlev + 1)
    real(kind=realkind) :: zdx_kf(nhor),zdy_kf(nhor),dtime  

    real(kind=realkind) :: u(nhor,nlev),v(nhor,nlev),t(nhor,nlev),q(nhor,nlev),&
         pf(nhor,nlev),tke(nhor,nlev),ps(nhor),ts(nhor)

    real(kind=realkind) :: p200,p400,p500,zpkm,zpkp,dphalf,rdphalf,zlogp,&
         zlogm,zdlnpk,zalfa,zbeta,ztv,es,qs,zphi,&
         sed(nlev),sem(nlev),&
         sems(nlev),zphi_full(nhor,nlev),zphi_half(nhor,nlev + 1),&
         dz_full(nhor,nlev),dz_half(nhor,nlev),dp_half(nhor,nlev)

    !gj   output fields from kfcumulus
    real(kind=realkind) :: dtdt(nhor,nlev),dqdt(nhor,nlev),dqcdt(nhor,nlev),&
         raincv(nhor),snowcv(nhor),umf(nhor,nlev),udr(nhor,nlev),&
         qvap(nhor,nlev),qliq(nhor,nlev),qice(nhor,nlev)
    real(kind=realkind) :: orogsigm(nhor),om_cpbl(nhor)  

    !gj   shallow conv vars
    integer :: shal_cgj(nhor)  
    real(kind=realkind) :: ahyb_kf(nlev + 1)  

    !gjfree conv
    integer::k,nk,l200,l400,l500,isink,lsb(nhor)  
    integer::ncuyes,icuyes(nhor),nca(nhor),ncldck,nc
    integer::kstop,kstart,i,nshall

    real(kind=realkind)::clddep(nhor)  
    real(kind=realkind)::wkf1(nhor),wkf2(nhor),wkf3(nhor)  
    real(kind=realkind)::zfrtop(nhor)  
    integer::zfrpbl(nhor),zlevtop(nhor)  

    real(kind=realkind), parameter::cp = 1005.7_realkind
    real(kind=realkind), parameter::xlv = 2.5e6_realkind
    real(kind=realkind), parameter::grav = 9.8_realkind
    real(kind=realkind), parameter::zrgasd = 287.04_realkind
    real(kind=realkind), parameter::aliq = 613.3_realkind
    real(kind=realkind), parameter::bliq = 17.502_realkind
    real(kind=realkind), parameter::cliq = 4780.8_realkind
    real(kind=realkind), parameter::dliq = 32.19_realkind

    integer :: l800  
    real(kind=realkind) :: p800  

    nshall = 0  
    ahyb_kf(1) = 1.0_realkind  
    do k = 2,nlev + 1  
       ahyb_kf(k) = ahyb(k)  
    enddo

    do i = kstart,kstop  
       icuyes(i) = 0  
       lsb(i) = 0  
    enddo
    !  call to the convective part of the scheme
    ncuyes = 0  
    do i = kstart,kstop  
       if(nca(i)>ncldck) goto 60
       p200 = ps(i)-2.e4_realkind  
       p400 = ps(i)-4.e4_realkind  
       p500 = ps(i)-5.e4_realkind  
       p800 = ps(i)-8.e4_realkind  

       !gj     calculate delp between model half levels(dphalf) and the
       !gj     full model level geopotential(zphi)(following 2.1.2.25 in
       !gj     hirlam 2.5 doc). use these to calculate static energy(sed)
       !gj     moist stat energy(sem) and sat. sem(sems)
       zphi = 0.0_realkind  
       zpkm = ps(i)  
       zlogm = log(zpkm)  

       do k = nlev,1,-1  
          nk = k + 1  
          zphi_half(i,nk) = zphi  
          zpkp = zpkm  
          zpkm = ahyb_kf(k) + bhyb(k)*ps(i)  
          dphalf = zpkp-zpkm
          dp_half(i,k) = dphalf  
          rdphalf = 1.0_realkind/dphalf  
          zlogp = zlogm  
          zlogm = log(zpkm)  
          zdlnpk = zlogp-zlogm  
          zalfa = 1.0_realkind-zpkm*rdphalf*zdlnpk  
          zbeta = zdlnpk-zalfa  

          ztv = t(i,k)*(1.0_realkind + 0.608_realkind*q(i,k))  
          !gj       calculate the full(level k) geopotential height
          !gj       inegrating upward from the k+1/2 level directly
          !gj       below. before looping around to k+1 full level
          !gj       we will later in this loop integrate from k to k-1/2
          !gj       to update zphi before this calc to be that representative
          !gj       of the k-1/2 level
          zphi = zphi + zalfa*zrgasd*ztv  
          zphi_full(i,k) = zphi  
          es = aliq*exp((bliq*t(i,k)-cliq)/(t(i,k)-dliq))

          qs = epsilo*es/(pf(i,k)-es)  
          sed(k) =(cp*ztv) + zphi
          sems(k) = sed(k) + xlv*qs  
          sem(k) = sed(k) + xlv*q(i,k)  
          !gj     find model full levels 200 & 400mb above sfcp
          if(pf(i,k) >p200) l200 = k  
          if(pf(i,k) >p400) l400 = k  
          if(pf(i,k) >p500) l500 = k  

          l800=nlev
          if(pf(i,k) >p800.and.pf(i,k) >20000.0_realkind) l800 = k  
          !gj
          !gj   finished using full level zphi update zphi to be k-1/2 level
          !gj   value as input in next portion of loop to get full level k-1
          !gj   geopotential value
          zphi = zphi + zbeta*zrgasd*ztv  
          if(k==1) zphi_half(i,k) = zphi  
          dz_half(i,k) = zphi-zphi_half(i,nk)  
       enddo
       dz_full(i,1) = 0.0_realkind  
       do k = nlev,2,-1  
          dz_full(i,k) = zphi_full(i,k-1)-zphi_full(i,k)  
       enddo
       !
       !gj     check to see if vertical velocity upward at some level
       !gj     in lowest 400mb
       isink = 0  

       do  k = nlev,l800,-1  
          !vert vel is
          if(omega_kf(i,k) >0.0_realkind.or.tke(i,k) >0.3_realkind) goto 25
       enddo
       isink = 1  


25     continue  
       !gj   if air sinking everywhere in lowest 400mb. can only have/need
       !gj   convection if we have a superadiabatic layer. check lowest
       !gj   400mb for one. if we find a superadiabatic layer ncuyes
       !gj   is incremented and icuyes location of potential convecting
       !gj   point is savec in reduced convection array icuyes.
       !subsidence from sfc to 400mb

       if(isink==1) then
          do k = nlev,l800,-1  
             do nk = k-1,l800,-1  
                if(sed(k) >sed(nk)) then  
                   ncuyes = ncuyes + 1  
                   icuyes(ncuyes) = i  
                   lsb(i) = nlev-k + 1  
                   goto 60  
                endif
             enddo
          enddo
          !motion upward
       else
          !
          !gj  if there is upward motion,require conditional instability for
          !gj  a parcel originating in the lowest 400mb. as above if one
          !gj  found then save location to icuyes and jump to testing next
          !gj  i point in do 60 loop.
          !gj  lsb(i) lowest model full level at which convection can be
          !gj  supported by below test. it is inverted in the vertical as
          !gj  kain-fritsch scheme runs sfc(k=1) to top(k=nlev).
          do  k = nlev,l800,-1  
             do  nk = k-1,1,-1  
                if(sem(k) >sems(nk)) then  
                   ncuyes = ncuyes + 1  
                   icuyes(ncuyes) = i  
                   lsb(i) = nlev-k + 1  
                   goto 60  
                endif
             enddo

          enddo
          !gj260404

       endif
60     continue
    enddo

    if(ncuyes>0)then
       call kfcumulus(nhor,nlev,nca,       &
            ncldck,dtime,ncuyes,icuyes,lsb,&
            zdx_kf,zdy_kf,ps,pf,t,q,            &
            u,v,tke,orogsigm,                   &
            zphi_full,dz_full,        &
            dz_half,dp_half,omega_kf,div_kf,    &
            dtdt,dqdt,dqcdt,raincv,             &
            kf_ind,kf_base,kf_top,              &
            umf,udr,om_cpbl,                    &
            zfrpbl,zfrtop,zlevtop,              &
            qvap,qliq,qice,clddep,              &
            shal_cgj,wkf1,wkf2,wkf3)             
    endif
    !gj  need to invert in the vertical tendencies from kfcumulus
    !gj  as these run from 1(lowest model level) to nlev(top model level).
    !gj  also invert base of convection for use in stamic

    do nc=1,ncuyes
       i=icuyes(nc)
       !gjkf convert some of convective rain into convective snow
       !gjkf do this as a function of lowest model temp. this will
       !gjkf be changed but is done just to get the diags in stamic
       !gjkf working.

       if(t(i,nlev-1)<tmelt .and. t(i,nlev)<tmelt.and.ts(i)<tmelt)then
          snowcv(i)=raincv(i)
          raincv(i)=0.0_realkind
       else
          snowcv(i)=0.0_realkind
       endif
    enddo
    return
  end subroutine akfcum


  subroutine kfcumulus(nhor,nlev,nca,ncldck, dtphys,&
       ncuyes,icuyes,lsb,zdx,zdy,ps,pf,t,q,u,v,cbr_tke,&
       orogsigm,zphi_full,dz_full,dz_half,dp_half,&
       omega_kf,div_kf,dtdt,dqdt,dqcdt,rainrate,kf_ind,kf_base,&
       kf_top,umf_out,udr_out,om_cpbl,zfrpbl,zfrtop,zlevtop,&
       qvap_cj,qliq_cj,qice_cj,clddep,fbfrc_out,wkf1,wkf2,wkf3)

    use confys,only:epsilo,tmelt
    use kflut
    use escom
    implicit none  

    integer :: i,k,nhor,nlev, nca(nhor),ncuyes,icuyes(nhor),&
         nc,lsb(nhor),kf_ind(nhor),kf_base(nhor),kf_top(nhor),&
         ncldck,kftrig(nhor)

    integer :: nccnt,nchm,ishall,ml,kl,klm,nk,l5, llfc,&
         kmix,low,lc,nlayrs,kpbl,lcl,let,klcl,kstart

    real(kind=realkind) :: zdx(nhor),zdy(nhor),ps(nhor),pf(nhor,nlev),&
         t(nhor,nlev),q(nhor,nlev),u(nhor,nlev),v(nhor,nlev),&
         cbr_tke(nhor,nlev),zphi_full(nhor,nlev),&
         dz_full(nhor,nlev),dp_half(nhor,nlev),dz_half(nhor,nlev),&
         omega_kf(nhor,nlev),dtphys,div_kf(nhor,nlev),&
         umf_out(nhor,nlev),dmf_out(nhor,nlev),udr_out(nhor,nlev),&
         uer_out(nhor,nlev),ddr_out(nhor,nlev)
    !     gj  look for snow_fb and rain_fb as these are related to
    !     gj  these vars and give a rain and snow flux that can be
    !     gj  related to convective snow and rain rates subsequently

    real(kind=realkind) :: p00,t00,g,cp,rlf,ttfrz,tbfrz,&
         rate,rad,r,aliq,bliq,cliq,dliq,aice,bice,cice,dice,&
         xlv0,xlv1,rovg,gdry,dt2
    real(kind=realkind) :: fbfrc,dmffrc,dxsq,dx,p200,p300,p500,pm15,p600  

    real(kind=realkind) :: p0(nlev),t0(nlev),q0(nlev),es,qes(nlev),rh(nlev) &
         ,ql0(nlev),qi0(nlev),qr0(nlev),qs0(nlev),u0(nlev),&
         v0(nlev),tv0(nlev),rhoe(nlev),dzq(nlev),dp(nlev),&
         w0(nlev),z00(nlev),dza(nlev),cldhgt(nlev),theteu(nlev),&
         thtes(nlev),div0(nlev)

    real(kind=realkind) :: dpthmx,qef,quer(nlev),tmix,qmix,zmix,pmix,emix,&
         tdpt,tlcl,tvlcl,zlcl,dlp,tenv,qenv,tven,wklcl,wsigne,&
         dtlcl,wabs,u00,qslcl,rhlcl,dqsdt,dtrh,thmix,tvavg,gdt,&
         plcl,qese,wtw,rholcl,wkl,wlcl,rg
    real(kind=realkind) :: astrt,ainc,a1,tp,value,aintrp,tlog  


    integer :: indlu  
    !     gj     variables used in updraft calculations.

    integer :: iflag,kfrz,ltop,nj,ltop1,ltopm1,lvf,nk1,nic,lfc
    real(kind=realkind) :: wu(nlev),tu(nlev),tvu(nlev),qu(nlev),qliq(nlev),&
         qice(nlev),qlqout(nlev),qicout(nlev),detlq(nlev),detic(nlev),&
         pptliq(nlev),pptice(nlev),umf(nlev),ratio2(nlev),&
         uer(nlev),eqfrc(nlev),thetee(nlev),udr(nlev),tvqu(nlev),&
         qdt(nlev),ems(nlev),emsd(nlev),tg(nlev),qg(nlev),&
         qlg(nlev),qig(nlev),qrg(nlev),qsg(nlev),omg(nlev + 1),&
         exn(nlev),thtau(nlev),thta0(nlev)


    real(kind=realkind) :: au0,vmflcl,upold,upnew,abe,trppt,&
         ttemp,pmix0,qfrz,qnewic,qnewlq,be,boterm,enterm,dzz,&
         wsq, rl,rei,ee2,ud2,ttmp,f1,f2,thttmp,qtmp,tmpliq,&
         tmpice,tu95,tu10,chmin,dptt,dumfdp,tsat,thta,p150,usr,&
         wspd(nlev),vconv,timec,tadvec,shsign,vws,pef,cbh,rcbh,&
         peff,peff2,ee1,ud1,ee,pefcbh,frc1
    !     gj     variables used in downdraft calculations.

    integer ::  klfs,lfs,nd,nd1,ldt,ldb,lmax,ncount,noitr,istop,ndk


    real(kind=realkind) :: tder,dppp,qd(nlev),tz(nlev),qss,&
         thtad(nlev),tvd(nlev),rdd,dmf(nlev),der(nlev),ddr(nlev),&
         pptmlt,dtmelt,qsrh,t1rh,pptflx,cpr,cndtnf,updinc,wd(nlev) &
         ,aincm2,ddinc,aincmx,tder2,pptfl2,detic2(nlev),udr2(nlev) &
         ,uer2(nlev),ddr2(nlev),der2(nlev),umf2(nlev),dmf2(nlev),&
         fabe,stab,aincm1,detlq2(nlev),dssdt,rhh,dtmp,qsd(nlev),&
         theted(nlev),dpdd
    !     gj     variables used in calculation of compensational subsidence
    integer :: nstep,ntc  


    real(kind=realkind) :: dtt,dtt1,dtime,tma,tmb,tmm,bcoeff,acoeff,qvdiff,&
         topomg,domgdp(nlev),thpa(nlev),qpa(nlev),fxm(nlev),&
         thfxin(nlev),thfxout(nlev),qfxin(nlev),qfxout(nlev),&
         thtag(nlev),tvg(nlev)
    !     gj     variables used in new cloud and buoyant energy calcs

    real(kind=realkind) :: rocpq,cpm,dq,cporq,thtesg(nlev),abeg,dabe,dfda,&
         aincold,fabeold,frc2,qlpa(nlev),qipa(nlev),qrpa(nlev),&
         qspa(nlev),rainfb(nlev),snowfb(nlev),qlfxin(nlev),qlfxout &
         (nlev),qifxin(nlev),qifxout(nlev),qrfxin(nlev),qrfxout(nlev),&
         qsfxin(nlev),qsfxout(nlev),qinit,qfnl,dpt,err2,&
         relerr,rnc
    real(kind=realkind) :: dtdt(nhor,nlev),dqdt(nhor,nlev),dqcdt(nhor,nlev),&
         raincv(nhor),rainrate(nhor)
    !     gjkfarm

    real(kind=realkind) :: qvap_cj(nhor,nlev),qliq_cj(nhor,nlev),qice_cj(nhor,nlev)
    !     gjkfarm
    !     gjkfvars
    real(kind=realkind) ::rad_kf(nhor)  
    !     gj
    !     gj  extra shallow conv variables
    !     gj
    real(kind=realkind) :: tke(nlev),dilfrc(nlev),pclb(nhor),umfb(nhor),&
         psrc(nhor),emst,evac,tkemax,rhbar,plcl0,dpmin,chmax
    integer::nloop
    real(kind=realkind) :: shal_cgj(nhor)  
    !     gj
    !     gjfree conv
    real(kind=realkind) :: orogsigm(nhor)  


    !     gj3dppt
    integer ::    cldbase,ltop_cj  

    integer :: fbfrc_out(nhor),fbfrc_shal  

    real(kind=realkind) :: dptt_cj(nlev),uer_shal  

    integer :: cldlow  
    !     gj250404
    !     gjomcpbl
    real(kind=realkind) :: wkf1(nhor),wkf2(nhor),wkf3(nhor)  
    real(kind=realkind) :: oro(nhor)  
    real(kind=realkind) :: tvmix,themix,exnmix,tvneg,dpneg,theneg,exneg,&
         dz_free,tvbase,tvtop,thedz,exncj(nlev),thecj(nlev),zsource,zrf1,zrf2
    integer ::  zitop,zlevtop(nhor),zfrpbl(nhor)  
    !     gjomcpbl
    real(kind=realkind) :: clddep(nhor)  
    real(kind=realkind) :: tvcgj(nlev),exncgj(nlev),thecgj(nlev)  
    real(kind=realkind) :: zfrtop(nhor),om_cpbl(nhor)  
    real(kind=realkind) :: dlp_sce,tenv_sce,qenv_sce,tven_sce,dz_sce,rhoe_sce,psce,tvp_sce
    real(kind=realkind) :: divmix,pert1,pert2,zcgj1,zcgj2,zcgj3,zcgj4,om_pert,&
         exndry(nlev),tvdry(nlev),wneg,cin_neg,cin_rf,zltop_kf,pert1_base
    integer :: lsce  
    !     gjnewent
    !     gjnewent cape now uses an entrainment influenced dilute ascent
    !     gjnewent profile rather than the original undilute ascent cape
    !     gjnewent
    real(kind=realkind) :: dilbe,ddilfrc(nlev),tgu(nlev),qgu(nlev),thteeg(nlev)
    real(kind=realkind) :: dilout(nhor,nlev),dilfrc2(nlev)  
    !     gjnewterminal detrainment
    integer :: ldetst  

    real(kind=realkind) :: dumfdz,dztt  
    !     gjnewent
    real(kind=realkind) :: thesmin  
    integer :: thelev  



    logical :: iprnt  
    !     gj   to activate the new trigger function remove the comments
    !     gj   cgjnewtrig from all lines where it is present.
    !     gjnewtrig
    !     gjesat-----------------------------------------------------------
    !     gjesat
    !     gjesat  use the hirlam reference esat calcs in kfcumulus to get a
    !     gjesat  accurate value at low temps that takes ice into account
    !     gjesat



    !  define constants needed for kfpara subroutine
    data p00,t00,g,cp/1.e5_realkind,273.16_realkind,9.81_realkind,1004.6_realkind/
    data rlf/3.339e5_realkind/
    data ttfrz,tbfrz/268.16_realkind,248.16_realkind/
    data rate,r/0.03_realkind,287.04_realkind/
    data aliq,bliq,cliq,dliq/613.3_realkind,17.502_realkind,4780.8_realkind,32.19_realkind/ 
    data aice,bice,cice,dice/613.2_realkind,22.452_realkind,6133.0_realkind,0.61_realkind/ 
    data xlv0,xlv1/3.147e6_realkind,2369.0_realkind/ 


    t00 = tmelt  
    rovg = r/g  
    gdry =-g/cp  
    rg = 1.0_realkind/g  
    dt2 = 2.0_realkind*dtphys  
    !
    !
    !     ! ppt fb mods
    !  option to feed convectively generated rainwater    ! ppt fb mods
    !  into grid-resolved rainwater(or snow/graupel)     ! ppt fb mods
    !  field.  'fbfrc' is the fraction of available       ! ppt fb mods
    !  precipitation to be fed back(0.0-1.0)        ! ppt fb mods
    ! ppt fb mods
    fbfrc = 0.0_realkind  
    fbfrc_shal = 0  
    !*
    !  loop over all potential convective grid points
    !*
    !
    nccnt = 0  
    iprnt = .false.  
    !     gjomcpbl
    zrf1 = 5.e-7_realkind  
    zrf2 = 3.e4_realkind  
    !     gjomcpbl
    !  specify downdraft mass flux as fraction of updraft mass flux(dmff
    !
    dmffrc = 0.9_realkind  
    !
    do nc = 1,ncuyes  
       !  mods to allow shallow convection
       nchm = 0  
       ishall = 0  
       rad = 1500.0_realkind  
       fbfrc = 0.0_realkind  
       dpmin = 5.e3_realkind  
       fbfrc_shal = 0  
       !
       i = icuyes(nc)  
       !     gjomcpbl
       wkf1(i) = 0.0_realkind  
       wkf2(i) = 0.0_realkind  
       wkf3(i) = 0.0_realkind  
       !     gjomcpbl
       dxsq = zdx(i)*zdy(i)  

       dx = sqrt(dxsq)  
       p200 = ps(i)-2.0e4_realkind  
       p300 = max((ps(i)-3.0e4_realkind),3.0e4_realkind)  
       p600 = max((ps(i)-6.0e4_realkind),4.0e4_realkind)  
       !     gj150905
       !halfway between psfc and 100hp
       p500 = ps(i)-((ps(i)-1.0e4_realkind)/2.0_realkind)  
       !     gj150905
       ml = 0  
       thelev = 0  
       !     gj150905
       llfc = 0  
       !     gj150905
       !  allow printed output for grid points within the(expanded) meaprs
       !  domain
       !
       !  define number of layers above ground level-call this kl
       !
       !     gjomcpbl
       zitop = nlev-zlevtop(i) + 1  
       zitop = min(nlev,zitop)  
       zltop_kf = real(nlev,realkind)-real(zlevtop(i),realkind) + 1._realkind  
       !     gjomcpbl
       kl = nlev  
       klm = kl-1  
       !
       !
       !     input a vertical sounding
       !
       !     djs  do loop 15 is crucial.  the kf scheme it set up to have level
       !     djs  starting with 1 at the lowest model level and going up,which
       !     djs  inverted from what sigma and eta coordinate models tend to us
       !     djs  thus,the first step is to switch the order of the variables
       !     djs  within the kf scheme.  note that the scheme needs the followi
       !     djs  variables:
       !     input:   temperture(t0,k) ;    specific humidity(q0,kg/kg) ; c
       !     horizontal wind speed(u0 and v0,m/s) ;                c
       !     pressure(p0,pascal) ;  height(z0,m);                c
       !     vertical motion(w0,m/s).                              c
       !     djs
       do  k = 1,kl  
          nk = nlev-k + 1  
          !     gj   pressure here must be that relevant to the colocated
          !     gj   temp and q as it is used below for rh etc calcs.
          !     gj   these must be done on full model levels ie
          !     gj   we can pass pf in from acondens and let it get
          !     gj   get set to p0 here.
          p0(k) = pf(i,nk)  
          t0(k) = t(i,nk)  
          q0(k) = q(i,nk)  
          tke(k) = cbr_tke(i,nk)  
          div0(k) = div_kf(i,nk)  
          !     gj
          exncgj(k) =(p00/p0(k))**(0.2854_realkind*(1.0_realkind-0.2854_realkind*q0(k)))
          tvcgj(k) = t0(k)*(1.0_realkind-0.608_realkind*q0(k))  
          thecgj(k) = t0(k)*exncgj(k)  
          !     gj
          !
          !  saturation vapor pressure(es) is calculated using formula given b
          !  buck(1981).  if q0 is above saturation value using this method,
          !  reduce it to saturation level
          !
          !     gjesat   call esat tables for a more accurate calculation of es
          !     gjesat   and calculate qes more correctly
          !     gjesat           es=esat(t0(k))
          !     gjesat           qes(k)=0.622*es/(p0(k)-(0.378*es))
          !     gjesat
          es = aliq*exp((bliq*t0(k)-cliq)/(t0(k)-dliq))  
          qes(k) = epsilo*es/(p0(k)-es)  
          q0(k) = min(qes(k),q0(k))  
          if(q0(k) <=0.1e-8_realkind) then  
             q0(k) = 0.1e-8_realkind  
          endif
          rh(k) = q0(k)/qes(k)  
          !
          !  set hydrometeor concentrations to zero initially.  although some
          !  hydrometeors may actually be present,this is done to disallow
          !  entrainment of these into convective updrafts and downdrafts.
          !  necessary because convective scheme operates over multiple time
          !  steps assuming steady-state environmental conditions and there
          !  is no way of guaranteeing that liquid water will continue to
          !  remain available for entrainment even if it is there initially
          !
          ql0(k) = 0.0_realkind  
          qi0(k) = 0.0_realkind  
          qr0(k) = 0.0_realkind  
          qs0(k) = 0.0_realkind  
          !     gjrca3           dilfrc(k) = 0.0_realkind
          dilfrc(k) = 1.0_realkind 
          dilfrc2(k) = 1.0_realkind  
          !     gjout
          dilout(i,k) =-1.0_realkind  
          !     gjout
          !
          !  calculate wind at h-point as average of 4 surrounding v-point valu
          !
          u0(k) = u(i,nk)  
          v0(k) = v(i,nk)  
          tv0(k) = t0(k)*(1.0_realkind + 0.608_realkind*q0(k))  
          rhoe(k) = p0(k)/(r*tv0(k))  
          !     gjomcpbl
          exncj(k) =(p00/p0(k))**(0.2854_realkind*(1.0_realkind-0.2854_realkind*q0(k)))
          thecj(k) = t0(k)*exncj(k)  
          !     gjomcpbl
          !
          dzq(k) = dz_half(i,nk)*rg  
          dp(k) = dp_half(i,nk)  
          w0(k) = omega_kf(i,nk)  
          !
          !  dzq is dz between eta surfaces,dza is dz between model half level
          !     dp is the pressure interval between full eta levels
          if(p0(k) >=p500) l5 = k  
          !test ps-600hpa or to 500pa
          if(p0(k) >=p600) llfc = k  
          if(t0(k) >t00) ml = k  
          cldhgt(k) = 0.0_realkind  
       enddo
       !ps is so high we only test for convecti
       if(llfc==0) llfc = 1  
       !     1st model level
       do  k = 1,kl  
          nk = nlev-k + 1  
          z00(k) = zphi_full(i,nk)*rg  
          dza(k) = dz_full(i,nk)*rg  

       enddo

       nloop = 0
       kmix = lsb(i)  
       !     gj  we come back to here after shallow convection tests
25     low = kmix  
       !     gj  we come back to here after shallow convection tests
       nloop = nloop + 1  
       if(nloop>100) then  
          goto 265  
       endif
       !
       !  if parcel originates from 300 mb above surface or higher then sche
       !     is not used-check next grid point
       !
       !  mods for shallow convectionllfc is 300mb above sfc pressure
       !     if(low>llfc)goto 325
       !     gj   this implies we keep looking for a deep convective sourcs lay
       !     gj   up yo llfc only if we haven't found one do we then allow
       !     gj   for the possible testing of potential shallow convective poin
       !     gj   for these points ishall will have been set to 1 and a cldhgt
       !     gj   will be available.loop 527 finds the source level that
       !     gj   supports the deepest cloud in this potential shallow convecti
       !     gj   column.the source layer for shallow convection is the nk
       !     gj   layer that supports the deepest cloud between sfc & llfc
       !     gj   nchm gets set to this vertical level and chmax to the cldhgt
       !     gj   from this level
       !     gj   first time we find a shallow convective point and we have che
       !     gj   up to level llfc without finding a deep convective point,we
       !     gj   have low>llfc & ishall==1 so we find the deepest cloud a
       !     gj   nchm equal to the source level for that cloud and set kmix to
       !     gj   we now go back to 25 low=kmix. presumably now low<llfc and
       !     gj   skip the below if test and set lc=low=kmix=nchmwe now have
       !     gj   ishall==1 & lc==nchm so we set a rad=1500 for shallow clo
       !     gj
       !     gjkfarm presently set to find source layer(lc) between 1 to llfc
       !     gjkfarm that supports the deepest cloud. change this to find the
       !     gjkfarm first layer from the surface that can support a cloud of
       !     gjkfarm any depthtested in arm case no difference
       !     gjkfarm
       if(low>llfc) then  
          if(ishall==1) then  
             !     kmix = lsb(i)
             chmax = 0.0_realkind  
             nchm = 0  
             do nk = 1,llfc  
                if(cldhgt(nk) >chmax) then  
                   nchm = nk  
                   chmax = cldhgt(nk)  
                endif
             enddo
             !set kmix to the diagnosed source of sha
             kmix = nchm  
             ! and loop back to 25 continue to do ful
             goto 25  
          endif
          goto 325  
       endif
       lc = low  
       if(ishall==1.and.lc==nchm) then  
          !     gjshallow give shallow clds a small rad
          rad = 300.0_realkind  
          !     gj
          fbfrc_shal = 1  
          fbfrc = 0.0_realkind  
          !     gjkfarm
          fbfrc = 1.0_realkind  
          !     gjkfarm
          !  3_24_99comment out statemnet below so that shallow cloud draws
          !  from same layer as deep cloud-to prevent shallow cloud
          !  from becoming too deep!
          !
          !     gjkfarm
          dpmin = 2.5e3_realkind  
          !     gjkfarm
       endif
       !
       !  assume that in order to support a deep updraft you need a layer of
       !  unstable air 50 to 100 mb deepto approximate this,isolate a
       !  group of adjacent individual model layers,with the base at level
       !  lc,such that the combined depth of these layers is at least 60 mb
       !
       !
       !cc   for eta model,introduce a moisture perturbation in updraft sou
       !cc   layers if rh > 75%note that you should consider changing thi
       !cc   approach if the grid-scale condensation threshold is modified..
       nlayrs = 0  
       dpthmx = 0.0_realkind  
       do nk = lc,kl  
          dpthmx = dpthmx + dp(nk)  
          nlayrs = nlayrs + 1  
          qef = 0.0_realkind  
          quer(nk) = q0(nk)*(1.0_realkind + qef)  
          if(dpthmx>dpmin) goto 64  
       enddo
       goto 325  
64     kpbl = lc + nlayrs-1  
       !
       !  go ahead and determine what level to start with for the
       !  next mixture in case the current mixture,with base at
       !  level lc,is not buoyant
       !  instead of checking mixtures using every single layer,
       !  move up in increments of at least 15 mb
       pm15 = p0(lc)-15.e2_realkind  
       do  nk = lc + 1,kl  
          if(p0(nk) <pm15) then  
             kmix = nk  
             goto 267  
          endif
       enddo
       goto 325  

267    continue  
       !  for computational simplicity without much loss in accuracy,
       !  mix temperature instead of theta for evaluating convective
       !  initiation(triggering) potential
       tmix = 0.0_realkind  
       qmix = 0.0_realkind  
       zmix = 0.0_realkind  



       pmix = 0.0_realkind  
       !  find the thermodynamic characteristics of the layer by
       !  mass-weighting the characteristics of the individual model
       !  layers
       tkemax = 0.0_realkind  
       tvmix = 0.0_realkind  
       themix = 0.0_realkind  
       exnmix = 0.0_realkind  
       divmix = 0.0_realkind
       do nk = lc,kpbl  
          tmix = tmix + dp(nk)*t0(nk)  
          qmix = qmix + dp(nk)*quer(nk)  
          zmix = zmix + dp(nk)*z00(nk)  
          pmix = pmix + dp(nk)*p0(nk)  
          tkemax = max(tkemax,tke(nk))  
          tvmix = tvmix + dp(nk)*tv0(nk)  
          themix = themix + dp(nk)*thecj(nk)  
          divmix = divmix + dp(nk)*div0(nk)  
       enddo
       tmix = tmix/dpthmx  
       qmix = qmix/dpthmx  
       zmix = zmix/dpthmx  
       pmix = pmix/dpthmx  
       emix = qmix*pmix/(epsilo + qmix)  
       !     gjomcpbl
       tvmix = tvmix/dpthmx  
       themix = themix/dpthmx  
       divmix = divmix/dpthmx  
       !     gjrftrig
       zsource = zmix  
       do nk = 1,kl  
          if(z00(nk) <zsource) lsce = nk  
       enddo
       !     gj     get tvenv(zource) & tv(parcel) at zsource
       dlp_sce =(zsource-z00(lsce))/(z00(lsce+1)-z00(lsce))  
       tenv_sce = t0(lsce) +(t0(lsce+1)-t0(lsce))*dlp_sce  
       qenv_sce = q0(lsce) +(q0(lsce+1)-q0(lsce))*dlp_sce  
       tven_sce = tenv_sce*(1.0_realkind + 0.608_realkind*qenv_sce)  
       !     gj
       !     gj    get delta-z between zsource and lsce in order to get
       !     gj    pressure ate zsource,namely psce and then exner function
       !     gj    and tv of parcel at zsource.
       !     gj
       dz_sce = zsource-z00(lsce)  
       rhoe_sce =(rhoe(lsce) + rhoe(lsce+1))/2.0_realkind  
       psce = pmix  
       tvp_sce = tvmix  
       !     gj
       !     gj     now we have the source height we can caluclate the
       !     gj     new vertical velocity perturbations as in rogers & fritsch
       !     gj     we need dz_free which is the vertical distance between the
       !     gj     top of the free convective pbl and zsource and the gross
       !     gj     stability of this layer. we also need the mean virtual
       !     gj     potential temp of this layer.
       !     gj
       !     gjkfnewtrigger
       if(oro(i) >=0.05_realkind) then  
          pert1_base = 0.2_realkind*dx/10000.0_realkind*(1.0_realkind-(zsource/5000.0_realkind))  
          if(orogsigm(i) <=100.0_realkind) then  
             pert1 = pert1_base  
          else  
             pert1 = pert1_base*(((orogsigm(i)-100.0_realkind)/300.0_realkind) + 1.0_realkind)
          endif
       else  
          pert1 = 0.0_realkind  
       endif
       !     gj
       !free conv pbl exists
       if(zfrpbl(i) ==1) then  
          !conv source is inside free pbl
          if(real(lsce,realkind)<=zltop_kf) then  
             !om pert is free conv vert. vel.
             pert2 = om_cpbl(i)  
          else  
             !     gj
             !     gj     calculate contribution of free conv vert vel. to convective
             !     gj     source above free conv. pbl.
             !     gj
             exneg = 0.0_realkind  
             tvneg = 0.0_realkind  
             dpneg = 0.0_realkind  
             theneg = 0.0_realkind  
             dz_free = zsource-zfrtop(i)  
             !     gj
             do nk = nint(zltop_kf),lsce   
                theneg = theneg + dp(nk)*thecgj(nk)  
                dpneg = dpneg + dp(nk)  
                if(nk==lsce) then  
                   tvbase = t0(nk)*exncgj(nk)  
                endif
                if(nk==nint(zltop_kf)) then  
                   tvtop = t0(nk)*exncgj(nk)  
                endif
             enddo
             !     gj
             theneg = theneg/dpneg  
             thedz =(tvtop-tvbase)/dz_free  
             zcgj1 = 5.e-7_realkind*(dz_free**2.0_realkind)  
             zcgj2 = 3.e4_realkind*thedz/theneg  
             zcgj3 = zcgj1 + zcgj2  
             zcgj4 = exp(-zcgj3)  
             pert2 = om_cpbl(i)*zcgj4  
             !end test for source loc relative to zfr
          endif
          !no free conv pbl
       else  
          pert2 = 0.0_realkind  
       endif
       !     gj
       !divergence
       if(divmix>=0.0_realkind) then  
          om_pert =(pert1 + pert2)*(1.0_realkind-10.0_realkind*(divmix**(1.0_realkind/3.0_realkind)))!0.333_realkind))  
       else  
          !convergence
          divmix = abs(divmix)  
          om_pert =(pert1 + pert2)*(1.0_realkind + 10.0_realkind*(divmix**(1.0_realkind/3.0_realkind)))!0.333_realkind))  
       endif
       !     cgj
       !     gj
       !     gj     later cin1 is calculated as the cine that exists between
       !     gj     zsource and zlfc for convection to occur we must have
       !     gj     om_pert > sqrt(cin1).
       !     gj     this is the new trigger function
       !     gj
       !
       !  find the temperature of the mixture at its lcl
       !
       !     tlog=log(emix/aliq)
       astrt = 1.e-3_realkind  
       ainc = 0.075_realkind  
       a1 = emix/aliq  
       tp =(a1-astrt)/ainc  
       indlu = nint(tp) + 1  
       value =real(indlu-1,realkind)*ainc + astrt  
       aintrp =(a1-value)/ainc  
       tlog = aintrp*alu(indlu + 1) +(1._realkind-aintrp)*alu(indlu)  
       tdpt =(cliq-dliq*tlog)/(bliq-tlog)  
       tlcl = tdpt-(0.212_realkind + 1.571e-3_realkind*(tdpt-t00)-4.36e-4_realkind*(tmix-t00))*(tmix-tdpt)
       tlcl = min(tlcl,tmix)  
       tvlcl = tlcl*(1.0_realkind + 0.608_realkind*qmix)  
       zlcl = zmix +(tlcl-tmix)/gdry  
       do nk = lc,kl  
          klcl = nk  
          if(zlcl<=z00(nk)) goto 35  
       enddo
       goto 325  
35     k = klcl-1  
       !
       !  calculate dlp using z instead og log(p)
       dlp =(zlcl-z00(k))/(z00(klcl)-z00(k))  
       !
       !  estimate environmental temperature and mixing ratio at the lcl
       !
       tenv = t0(k) +(t0(klcl)-t0(k))*dlp  
       qenv = q0(k) +(q0(klcl)-q0(k))*dlp  
       tven = tenv*(1.0_realkind + 0.608_realkind*qenv)  
       !
       !  check to see if cloud is buoyant using fritsch-chappell trigger
       !  function described in kain and fritsch(1992)w0 is an
       !  aproximate value for the running-mean grid-scale vertical
       !  velocity,which gives smoother fields of convective initiation
       !  than the instantaneous valueformula relating temperature
       !  perturbation to vertical velocity has been used with the most
       !  success at grid lengths near 25 km.  for different grid-lengths,
       !  adjust vertical velocity to equivalent value for 25 km grid
       !  length,assuming linear dependence of w on grid length
       !
       if(zlcl<2.e3_realkind) then  
          wklcl = 0.02_realkind*zlcl/2.e3_realkind  
       else  
          wklcl = 0.02_realkind  
       endif
       wkl =(w0(k) +(w0(klcl)-w0(k))*dlp)*dx/25.e3_realkind-wklcl
       wabs = abs(wkl)  
       if(abs(wabs)<1.e-14_realkind) then  
          wsigne = 1.0_realkind  
          dtlcl = 0.0_realkind  
          goto 26  
       endif
       wsigne = wkl/wabs  
       !     gjomcpbl
       !     ps059005        wkf1(i)=wkl
       !     ps059005        wkf2(i)=wkf
       !     ps059005        wkf3(i)=div0(lsource)
       !     gj        if(wkl>0.0_realkind)then
       !     gj         wknett=max(wkl,wkf)
       !     gj        else
       !     gj         wknett=wkf+wkl
       !     gj        endif
       !     gjomcpbl
       dtlcl = 4.64_realkind*wsigne*wabs**(1.0_realkind/3.0_realkind)!0.33_realkind  
       !     gj110605        if(wkf>0.0_realkind)then
       !     gj110605         dtlcl=4.64_realkind*wkf*0.33_realkind
       !     gj110605        else
       !     gj110605         dtlcl=0.0_realkind
       !     gj110605        endif
       dtlcl = max(dtlcl,0.0_realkind)  
26     continue  
       !
       !  for eta model,give parcel an extra temperature perturbation based
       !  the threshold rh for condensation
       !     gjfire140102
       !     gjfire rationale for applying the rh pertuurbation is that there
       !     gjfire is a sug-grid scale distribution of rh values about the gri
       !     gjfire box mean value. the distribution will be much wider as the
       !     gjfire surface forcing varies therefore the dtrh term can justifia
       !     gjfire be made a function of orogsigm. presently here we use orogs
       !     gjfire to indicate seas(orogsigm=0) and land(orogsigm>1) and dtrh
       !     gjfire assumed zero over sea where variance of lcl rh(due to surf
       !     gjfire variability forcing) is zero.
       !
       !  for now,just assume u00=0.75_realkind
       u00 = 0.75_realkind  
       if(u00<1.0_realkind) then  
          qslcl = qes(k) +(qes(klcl)-qes(k))*dlp  
          rhlcl = qenv/qslcl  
          dqsdt = qmix*(cliq-bliq*dliq)/((tlcl-dliq)*(tlcl-dliq))
          !     gjfire
          if(orogsigm(i) >0.0_realkind) then  
             !     gjrca3          if(rhlcl>=0.75 .and. rhlcl<=0.95)then
             if(rhlcl>=0.6_realkind.and.rhlcl<=0.95_realkind) then  
                !     gjrca3            dtrh = 0.25*(rhlcl-0.75)*qmix/dqsdt
                dtrh = 0.25_realkind*(rhlcl-0.6_realkind)*qmix/dqsdt  
             elseif(rhlcl>0.95_realkind) then  
                dtrh =(1.0_realkind/rhlcl-1.0_realkind)*qmix/dqsdt  
             else  
                dtrh = 0.0_realkind  
             endif
          else  
             !     gj          if(zlcl>500. .and. tke(1)>0.2)then
             if(zlcl>500.0_realkind.and.tke(1) >0.5_realkind) then  
                !     gjrca3            if(rhlcl>=0.75 .and. rhlcl<=0.95)then
                if(rhlcl>=0.6_realkind.and.rhlcl<=0.95_realkind) then  
                   !     gjrca3              dtrh = 0.25*(rhlcl-0.75)*qmix/dqsdt
                   dtrh = 0.25_realkind*(rhlcl-0.6_realkind)*qmix/dqsdt  
                elseif(rhlcl>0.95_realkind) then  
                   dtrh =(1.0_realkind/rhlcl-1.0_realkind)*qmix/dqsdt  
                else  
                   dtrh = 0.0_realkind  
                endif
             else  
                dtrh = 0.0_realkind  
             endif
             !     gjnodtrhocean          dtrh=0.0_realkind

          endif
       endif
       !     gj
       !     gj---------------------------------------------------------------
       !     gj
       !     gj     lfc has been defined higher up. old code here calculated
       !     gj     cin1 between zlcl-zlfc but did not consider the cin1
       !     gj     that may exist between the convective source region
       !     gj     and zlcl itself. add this in here and apply the new trigger
       !     gj     function test.
       !     gj
       !     gj     in the loop below we calculate cine from kpbl to
       !     gj     klcl-1up to this level we are sure processes are
       !     gj     dry. zlcl lies between klcl-1 & klcl.
       !     gk     we assume we calc cine from one level below kpbl
       !     gj     where kpbl is the top of the 60mb level being checked
       !     gj     bottom of this level is lc(ie 60mb layer goes lc-kpbl)
       !     gj     how sensitive is cine to where in layer we consider
       !     gj     relevant point to start calc. nb..rogers & fritsch use
       !     gj     the mid-height of the source layer.
       !     gj
       cin_rf = 0.0_realkind  
       !src layer below zlcl for cine to exist
       if(lsce<klcl-1) then  
          do nk = lsce+1,klcl-1  
             !     gj
             !     gj     assume thetav_mix is conserved during ascent to klcl-
             !     gj     during ascent calculate tvcgj(nk) as a function of
             !     gj     thetav_mix and the local pressure(po(nk))
             !     gj
             exndry(nk) = 1.0_realkind/exncgj(nk)  
             tvdry(nk) =(exndry(nk)*themix)*(1.0_realkind + 0.608_realkind*qmix)  
             !     gj
             if(nk==lsce+1) then  
                cin_rf = cin_rf + g*(z00(nk)-zsource)*((tvdry( &
                     nk) + tvp_sce)/(tv0(nk) + tven_sce)-1.0_realkind)
             else  
                cin_rf = cin_rf + g*(z00(nk)-z00(nk-1))*&
                     ((tvdry(nk) + tvdry(nk-1))/(tv0(nk) + tv0(nk-1))-1.0_realkind)
             endif
             !cin_rf contains cine top of source laye
          enddo
          !     gj
          !     gj    add in cine between klcl-1 & zlcl(this is dry!)
          !     gj    we can do this because below we use kf original and they
          !     gj    assumed cine to zlcl was zero and began calc from that point
          !     gj    as done below
          !     gj
          cin_rf = cin_rf + g*(zlcl-z00(klcl-1))*((tvlcl + &
               tvdry(klcl-1))/(tven + tv0(klcl-1))-1.0_realkind)
       else  
          cin_rf = 0.0_realkind  
       endif
       !     gj
       !     gj     we now have cine between lsource and zlcl,we now need
       !     gj     to compare the vertical velocity pertubations the source
       !     gj     parcel experiences with that needed to overcome this cine
       !     gj     to determine if convection is to be allowed. if conv allowe
       !     gj     we goto 45 as with the old trigger.
       !     gj     cin_rf at this point is negative if it is convective inhibi
       !     gj     to get inhibition vertical motion to overcome,multiply by
       !     gj
       if(cin_rf<0.0_realkind) then  
          cin_neg =-cin_rf  
          wneg = sqrt(2.0_realkind*cin_neg)  
       else  
          cin_neg = 0.0_realkind  
          wneg = 0.0_realkind  
       endif
       !     gj
       if(om_pert>=wneg.or.cin_neg==0.0_realkind) goto 45  
       !     gj
       !ps110701 if(tlcl+dtlcl+dtrh>tenv)goto 45
       !     gj
       if(ishall==1.and.lc==nchm) goto 325  
       if(kmix<=llfc) goto 25  
       !  shallow convection mods
       if(ishall==1) then  
          chmax = 0.0_realkind  
          nchm = 0  
          do nk = 1,llfc  
             if(cldhgt(nk) >chmax) then  
                nchm = nk  
                chmax = cldhgt(nk)  
             endif
          enddo
          kmix = nchm  
          goto 25  
       endif
       goto 325  
       !  convective triggering criteria has been satisfiedcompute
       !  equivalent potential temperature
       !(theteu) and vertical velocity of the rising parcel at the lcl
45     continue  
       !     gj
       kftrig(i) = 1  
       !     gj
       !     gjfree conv
       thmix = tmix*(1.e5_realkind/pmix)**(0.2854_realkind*(1.0_realkind-0.28_realkind*qmix))  
       theteu(k) = thmix*exp((3374.6525_realkind/tlcl-2.5403_realkind)*qmix*(1.0_realkind + 0.81_realkind*qmix))
       !     gjesat
       !     gjesat
       !     gjesat           es=esat(tenv)
       es = aliq*exp((tenv*bliq-cliq)/(tenv-dliq))  
       tvavg = 0.5_realkind*(tv0(klcl) + tenv*(1.0_realkind + 0.608_realkind*qenv))  
       !  modify calculation of initial parcel vertical velocityjsk 11/26
       !     gdt=g*dtlcl*(zlcl-z00(lc))/(tv0(lc)+tven)
       !     wlcl=1.+.5*wsigne*sqrt(abs(gdt)+1.e-14)
       gdt = 2.0_realkind*g*(dtlcl + dtrh)*500.0_realkind/tven  
       !     wlcl=1.+0.5*wsigne*sqrt(abs(gdt)+1.e-14)
       !     gjkfarm
       if(ishall==1.and.lc==nchm) then  
          wlcl = 1.0_realkind + 0.25_realkind*sqrt(abs(gdt) + 1.e-14_realkind)  
       else  
          wlcl = 1.0_realkind + 0.5_realkind*sqrt(abs(gdt) + 1.e-14_realkind)  
       endif
       !     gj210206
       !     if(orogsigm(i)<=0.0_realkind)then
       !     wlcl = 0.5*sqrt(abs(gdt)+1.e-14)
       !     endif
       !     gj210206
       !     gjkfarm
       wlcl = min(wlcl,3.0_realkind)  
       plcl = p0(k) +(p0(klcl)-p0(k))*dlp  
       !     gjesat
       !     gjesat        qese=0.622*es/(plcl-(0.378*es))
       !     gjesat
       qese = epsilo*es/(plcl-es)  
       thtes(k) = tenv*(1.e5_realkind/plcl)**(0.2854_realkind*(1.0_realkind-0.28_realkind* &
            qese))*exp((3374.6525_realkind/tenv-2.5403_realkind)*qese*(1.0_realkind + 0.81_realkind*qese))
       wtw = wlcl*wlcl  
       if(wlcl<0.0_realkind) goto 25  
       tvlcl = tlcl*(1.0_realkind + 0.608_realkind*qmix)  
       rholcl = plcl/(r*tvlcl)  
       plcl0 = plcl  
       !
       lcl = klcl  
       let = lcl  
       !     gj
       !     gj  this is to make rad a function of w at lcl.my
       !     gj  alternative follows within the cgjomega comments
       !     gjshallow
       !     gj  this if test was commented previously 31/5/99 if(lc/=nchm)th
       !     gj  this sets cloud radius for deep convection only
       !     gjjapplmet
       if(lc/=nchm) then  
          if(wkl<0.0_realkind) then  
             rad = 1000.0_realkind  
          elseif(wkl>0.1_realkind) then  
             rad = 2000.0_realkind  
          else  
             rad = 1000.0_realkind + 1000.0_realkind*wkl/0.1_realkind  
          endif
          !     gjjapplmet
          !     gjnewcode
          !     gjnewcode          if(wkl<-0.1)then
          !     gjnewcode            rad = 1000.
          !     gjnewcode          elseif(wkl>=-0.1 .and. wkl<=0.0)then
          !     gjnewcode            rad=1500.
          !     gjnewcode
          !     gjnewcode          elseif(wkl>0.1)then
          !     gjnewcode            rad = 2000.
          !     gjnewcode
          !     gjnewcode            rad = 1500.+500*wkl/0.1
       endif
       rad_kf(i) = rad  
       !
       !*
       !    *
       !     compute updraft properties                      *
       !    *
       !*
       !
       !
       !
       !  estimate initial updraft mass flux(umf(k))
       !
       wu(k) = wlcl  
       !     gj        au0=pie*rad*rad
       !spec for wrf
       au0 = 0.01_realkind*dxsq  
       umf(k) = rholcl*au0  
       vmflcl = umf(k)  
       upold = vmflcl  
       upnew = upold  
       !
       !  ratio2 is the degree of glaciation in the cloud(0 to 1),
       !  uer is the envir entrainment rate,abe is available
       !  buoyant energy,trppt is the total rate of precipitation
       !  production
       !
       ratio2(k) = 0.0_realkind  
       uer(k) = 0.0_realkind  
       abe = 0.0_realkind  
       trppt = 0.0_realkind  
       tu(k) = tlcl  
       tvu(k) = tvlcl  
       qu(k) = qmix  
       eqfrc(k) = 1.0_realkind  
       qliq(k) = 0.0_realkind  
       qice(k) = 0.0_realkind  
       qlqout(k) = 0.0_realkind  
       qicout(k) = 0.0_realkind  
       detlq(k) = 0.0_realkind  
       detic(k) = 0.0_realkind  
       pptliq(k) = 0.0_realkind  
       pptice(k) = 0.0_realkind  
       iflag = 0  
       kfrz = lc  
       !
       !  the amount of conv avail pot energy(cape) is calculated with
       !  respect to undilute parcel ascent; eq pot temp of undilute
       !  parcel is thtudl,undilute temperature is given by tudl
       !
       !     gjnewent
       !     gjnewent this is now modified so cape uses dilute parcel ascent
       !     gjnewent using calculated ascent profile including entrainment
       !     gjnewent effects
       !     gjnewent        thtudl=theteu(k)
       !     gjnewent        tudl=tlcl
       !
       !  ttemp is used during calculation of the linear glaciation
       !  process; it is initially set to the temperature at which
       !  freezing is specified to begin.  within the glaciation
       !  interval,it is set equal to the updraft temp at the
       !  previous model level
       !
       ttemp = ttfrz  
       !
       !  enter the loop for updraft calculationscalculate updraft temp,
       !  mixing ratio,vertical mass flux,lateral detrainment of mass and
       !  moisture,precipitation rates at each model level
       !
       !
       !     gj
       thesmin = 600.0_realkind  
       !     gj
       !
       !  !!! diagnostic stuff !!!!
       plcl0 = plcl  
       lfc = 0  
       pmix0 = pmix  
       rei = 0.0_realkind  
       !loop goes klcl-1-->kl-1
       do nk = k,klm  
          nk1 = nk + 1  
          ratio2(nk1) = ratio2(nk)  
          frc1 = 0.0_realkind  
          tu(nk1) = t0(nk1)  
          theteu(nk1) = theteu(nk)  
          qu(nk1) = qu(nk)  
          qliq(nk1) = qliq(nk)  
          qice(nk1) = qice(nk)  
          call tpmix2(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),&
               qliq(nk1),qice(nk1),qnewlq,qnewic,ratio2(nk1),&
               xlv1,xlv0)
          !
          !  check to see if updraft temp is above the temperature at which
          !  glaciation is assumed to initiate; if it is,calculate the
          !  fraction of remaining liquid water to freezettfrz is the
          !  temp at which freezing begins,tbfrz the temp below which all
          !  liquid water is frozen at each level
          !
          if(tu(nk1) <=ttfrz) then  
             if(tu(nk1) >tbfrz) then  
                if(ttemp>ttfrz) ttemp = ttfrz  
                frc1 =(ttemp-tu(nk1))/(ttemp-tbfrz)  
             else  
                frc1 = 1.0_realkind  
                iflag = 1  
             endif
             ttemp = tu(nk1)  
             !
             !     determine the effects of liquid water freezing when temperature
             !  is below ttfrz
             !
             qfrz =(qliq(nk1) + qnewlq)*frc1  
             qnewic = qnewic + qnewlq*frc1  
             qnewlq = qnewlq-qnewlq*frc1  
             qice(nk1) = qice(nk1) + qliq(nk1)*frc1  
             qliq(nk1) = qliq(nk1)-qliq(nk1)*frc1  
             call dtfrznew(tu(nk1),p0(nk1),theteu(nk1),qu(nk1) &
                  ,qfrz,qice(nk1),aliq,bliq,cliq,dliq)
          endif
          tvu(nk1) = tu(nk1)*(1.0_realkind + 0.608_realkind*qu(nk1))  
          !
          !     calculate updraft vertical velocity and precipitation fallout
          !
          if(nk==k) then  
             be =(tvlcl + tvu(nk1))/(tven + tv0(nk1))-1.0_realkind  
             boterm = 2.0_realkind*(z00(nk1)-zlcl)*g*be/1.5_realkind  
             !     gjnewenterm
             enterm = 0.0_realkind  
             dzz = z00(nk1)-zlcl  
          else  
             be =(tvu(nk) + tvu(nk1))/(tv0(nk) + tv0(nk1)) -1.0_realkind
             boterm = 2.0_realkind*dza(nk)*g*be/1.5_realkind  
             !     gjnewenterm
             enterm = 2.0_realkind*uer(nk)*wtw/upold  
             dzz = dza(nk)  
          endif
          !     gjtestasia         enterm=2.*rei*wtw/upold
          !  diagnostics
          if(tvu(nk1) >tv0(nk1)) then  
             if(tvu(nk) <tv0(nk) .or.nk1==klcl) lfc = nk1  
          endif
          wsq = wtw  
          !     gj
          if(ishall==1) then  
             rate = 0.0_realkind  
          else  
             rate = 0.03_realkind  
          endif
          !     gj1dcode
          !     gj111103         rate_phase=rate-(0.001*(ttfrz-tu(nk1)))
          !     gj         rate_phase=rate-(0.0003*(ttfrz-tu(nk1)))
          !     gj111103         rate=max(0.01,min(0.03,rate_phase))
          !     gj
          !     gj         if(oro(i)<=0.)rate=0.001
          !     gj         rate=max(0.015,min(0.03,rate_phase))
          !     gjtest
          !     gj1dcode
          !     gj         if(oro(i)<=0.0_realkind)rate=0.01
          call condload(qliq(nk1),qice(nk1),wtw,dzz,boterm,&
               enterm,rate,qnewlq,qnewic,qlqout(nk1),qicout(nk1))
          wabs = sqrt(abs(wtw))  
          wu(nk1) = wtw/wabs  
          !
          !  if vert velocity is less than zero,exit the updraft loop and,
          !  if cloud is tall enough,finalize updraft calculations
          !
          if(wu(nk1) <0.0_realkind) goto 65  
          !
          !     update the abe for undilute ascent
          !
          thtes(nk1) = t0(nk1)*(1.e5_realkind/p0(nk1))**(0.2854_realkind* &
               (1.0_realkind-0.28_realkind*qes(nk1)))*exp((3374.6525_realkind/t0(nk1)-2.5403_realkind) &
               *qes(nk1)*(1.0_realkind + 0.81_realkind*qes(nk1)))
          !     gj
          !     gj   save minimum in env theta-es for later use as source
          !     gj   level for downdraft
          if(thtes(nk1) <thesmin) then  
             thesmin = thtes(nk1)  
             thelev = nk1  
          endif
          !
          !  call subroutine to calculate environmental equivalent potential
          !  tempwithin glaciation interval,thetae must be calculated
          !  with respect to the same degree of glaciation for all entraining
          !  air
          !
          !  for lookup table version,calculate thetae with respect to
          !  liquid water at all levels
          !     call envirtht(p0(nk1),t0(nk1),q0(nk1),thetee(nk1),ratio2(nk1),
          !    *                rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          if(nk1<=kpbl) then  
             call envirtht(p0(nk1),t0(nk1),quer(nk1),thetee( &
                  nk1),0.0_realkind,rl,aliq,bliq,cliq,dliq,aice,bice,cice,&
                  dice)
          else  
             call envirtht(p0(nk1),t0(nk1),q0(nk1),thetee(nk1) &
                  ,0.0_realkind,rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          endif
          !
          !  rei is the rate of environmental inflow
          !
          rei = vmflcl*dp(nk1)*0.03_realkind/rad  
          tvqu(nk1) = tu(nk1)*(1.0_realkind + 0.608_realkind*qu(nk1)-qliq(nk1)-qice(nk1))
          !     gjnewent
          if(nk==k) then  
             dilbe =((tvlcl + tvqu(nk1))/(tven + tv0(nk1)) -1.0_realkind)*dzz
          else  
             dilbe =((tvqu(nk) + tvqu(nk1))/(tv0(nk) + tv0(nk1))-1.0_realkind)*dzz
          endif
          if(dilbe>0.0_realkind) abe = abe+dilbe*g  
          !     gjnewent
          !
          !  if cloud parcels are virtually colder than the environment,no
          !     entrainment is allowed at this level
          !
          if(tvqu(nk1) <=tv0(nk1)) then  
             !     uer(nk1)=0.0
             !     udr(nk1)=rei
             uer(nk1) = 0.5_realkind*rei  
             udr(nk1) = 1.5_realkind*rei  
             !     ee2=0.0_realkind
             !     ud2=1.
             ee2 = 0.5_realkind  
             ud2 = 1.0_realkind  
             eqfrc(nk1) = 0.0_realkind  
             goto 55  
          endif
          let = nk1  
          ttmp = tvqu(nk1)  
          !
          !  determine the critical mixed fraction of updraft and environmental
          !
          f1 = 0.95_realkind  
          f2 = 1.0_realkind-f1  
          thttmp = f1*thetee(nk1) + f2*theteu(nk1)  
          qtmp = f1*q0(nk1) + f2*qu(nk1)  
          tmpliq = f2*qliq(nk1)  
          tmpice = f2*qice(nk1)  
          call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,&
               qnewlq,qnewic,ratio2(nk1),xlv1,xlv0)
          tu95 = ttmp*(1.0_realkind + 0.608_realkind*qtmp-tmpliq-tmpice)  
          if(tu95>tv0(nk1)) then  
             ee2 = 1.0_realkind  
             ud2 = 0.0_realkind  
             eqfrc(nk1) = 1.0_realkind  
             goto 50  
          endif
          f1 = 0.10_realkind  
          f2 = 1.0_realkind-f1  
          thttmp = f1*thetee(nk1) + f2*theteu(nk1)  
          qtmp = f1*q0(nk1) + f2*qu(nk1)  
          tmpliq = f2*qliq(nk1)  
          tmpice = f2*qice(nk1)  
          call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,&
               qnewlq,qnewic,ratio2(nk1),xlv1,xlv0)
          tu10 = ttmp*(1.0_realkind + 0.608_realkind*qtmp-tmpliq-tmpice)  
          if(tu10==tvqu(nk1)) then  
             ee2 = 1.0_realkind  
             ud2 = 0.0_realkind  
             eqfrc(nk1) = 1.0_realkind  
             goto 50  
          endif
          eqfrc(nk1) =(tv0(nk1)-tvqu(nk1))*f1/(tu10-tvqu (nk1))
          eqfrc(nk1) = max(0.0_realkind,eqfrc(nk1))  
          eqfrc(nk1) = min(1.0_realkind,eqfrc(nk1))  
          if(abs(eqfrc(nk1)-1._realkind)<1.e-14_realkind) then  
             ee2 = 1.0_realkind  
             ud2 = 0.0_realkind  
             goto 50  
          elseif(abs(eqfrc(nk1))< 1.e-14_realkind) then  
             ee2 = 0.0_realkind  
             ud2 = 1.0_realkind  
             goto 50  
          else  
             !
             !  subroutine prof5 integrates over the gaussian dist to determine th
             !     fractional entrainment and detrainment rates
             !
             call prof5(eqfrc(nk1),ee2,ud2)  
          endif
          !
          !base of cloud ent=1,det=0
50        if(nk==k) then  
             ee1 = 1.0_realkind  
             ud1 = 0.0_realkind  
          endif
          !
          !  net entrainment and detrainment rates are given by the average fra
          !     values in the layer
          !
          !
          !     gjjapplmetv43
          !     gj         ee2 = max(ee2,0.5)
          !     gjreducedent
          ee2 = max(ee2,0.3_realkind)  
          !     gjreducedent
          ud2 = 1.5_realkind*ud2  
          !     gjbkf
          ee2 = 0.5_realkind  
          ud2 = 0.5_realkind  
          !     gjbkf
          uer(nk1) = 0.5_realkind*(umf(nk)*0.03_realkind*dp(nk1)/rad) *(ee1 + ee2)
          udr(nk1) = 0.5_realkind*(umf(nk)*0.03_realkind*dp(nk1)/rad) *(ud1 + ud2)
          !     gjbkf
          !     gjbkf         uer(nk1)=0.5*rei*(ee1+ee2)
          !     gjbkf         udr(nk1)=0.5*rei*(ud1+ud2)
          !
          !  if the calculated updraft detrainment rate is greater than the tot
          !     updraft mass flux,all cloud mass detrains,exit updraft calculati
          !
55        if(umf(nk)-udr(nk1) <10.0_realkind) then  
             !
             !  if the calculated detrained mass flux is greater than the total up
             !     flux,impose total detrainment of updraft mass at the previous mod
             !
             if(dilbe>0.0_realkind) abe = abe-dilbe*g  
             !     gjnewent
             let = nk  
             !     write(98,1015)p0(nk1)/100.0_realkind
             goto 65  
          endif
          ee1 = ee2  
          ud1 = ud2  
          upold = umf(nk)-udr(nk1)  
          upnew = upold+uer(nk1)  
          umf(nk1) = upnew  
          dilfrc(nk1) = upnew/upold  
          dilfrc2(nk1) =(umf(nk) + uer(nk1))/umf(nk)  
          !     gj
          dilout(i,nk1) = dilfrc2(nk1)  
          !     gj
          !
          !  detlq and detic are the rates of detrainment of liquid and
          !  ice in the detraining updraft mass
          !
          detlq(nk1) = qliq(nk1)*udr(nk1)  
          detic(nk1) = qice(nk1)*udr(nk1)  
          qdt(nk1) = qu(nk1)  
          if(nk1<=kpbl) then  
             qu(nk1) =(upold*qu(nk1) + uer(nk1)*quer(nk1)) /upnew
          else  
             qu(nk1) =(upold*qu(nk1) + uer(nk1)*q0(nk1)) /upnew
          endif
          theteu(nk1) =(theteu(nk1)*upold+thetee(nk1)*uer(nk1))/upnew
          qliq(nk1) = qliq(nk1)*upold/upnew  
          qice(nk1) = qice(nk1)*upold/upnew  
          !
          !  kfrz is the highest model level at which liquid condensate
          !  is generatedpptliq is the rate of generation(fallout) of
          !  liquid precip at a given model lvl,pptice the same for ice,
          !  trppt is the total rate of production of precip up to the
          !  current model level
          !
          if(abs(ratio2(nk1)-1.0_realkind) >1.e-6_realkind) kfrz = nk1  
          !  reverse the mod that allows feedback of rain/snow that originates
          !  detrained air12/9/98jsk
          !     pptliq(nk1)=qlqout(nk1)*(umf(nk)-udr(nk1))
          !     pptice(nk1)=qicout(nk1)*(umf(nk)-udr(nk1))
          pptliq(nk1) = qlqout(nk1)*umf(nk)  
          pptice(nk1) = qicout(nk1)*umf(nk)  
          !
          trppt = trppt + pptliq(nk1) + pptice(nk1)  
          if(nk1<=kpbl) uer(nk1) = uer(nk1) + vmflcl*dp(nk1)/dpthmx

       enddo
       !  check cloud depthif cloud is tall enough,estimate the equilibr
       !     temperature level(let) and adjust mass flux profile at cloud top
       !     that mass flux decreases to zero as a linear function of pressure
       !     the let and cloud top
       !
       !  ltop is the model level just below the level at which vertical vel
       !     first becomes negative
       !
65     ltop = nk  
       cldlow = 0  
       cldhgt(lc) = z00(ltop)-zlcl  
       clddep(i) = cldhgt(lc)  
       if(ltop==klcl) cldlow = 1  
       !     gj3dppt
       !     gj3dppt cldbase is saved here as teh first level at which
       !     gj3dppt a 3d precip field is generated in an updraft.this
       !     gj3dppt is at level k+1 in loop 60(ie nk1=nk+1)
       !     gj3dppt
       cldbase = k + 1  
       ltop_cj = ltop  
       !     gj3dppt
       !     gj     ltop is the top level of convection save this to kf_top
       !     gj     cldbase is saved as base of convection.
       !     gj     n.b. these have been inverted from kfcumulus vertical grid
       !     gj     to hirlam grid
       !     gj
       kf_top(i) = nlev-ltop + 1  
       kf_base(i) = nlev-cldbase+1  
       !     gj
       !
       !  if cloud top height is less than the specified minimum for deep
       !  convection,save value to consider this level as source for
       !  shallow convection,go back up to check next level
       !
       !  try specifying minimum cloud depth as a function of tlcl
       !
       !     gjjapplmet
       if(tlcl>293.0_realkind) then  
          !     gjrca3        chmin = 4.e3
          chmin = 3.e3_realkind  
       elseif(tlcl<=293.0_realkind.and.tlcl>=273.0_realkind) then  
          !     gjrca3        chmin = 2.e3 + 100.0_realkind*(tlcl-273.)
          chmin = 2.e3_realkind + 50.0_realkind*(tlcl-273.0_realkind)  
       elseif(tlcl<273.0_realkind) then  
          chmin = 2.e3_realkind  
       endif
       !     gjjapplmet
       !     print*,'chmin =',chmin
       !     chmin = 3.e3
       !     gjasiatest
       !     gj      chmin = 1.5e3
       !     gjasiatest
       kstart = max0(kpbl,klcl)  
       !
       !  do not allow any cloud from this layer if:
       !
       !  1.) if there is no cape,or
       !  2.) cloud top is at model level just above lcl,or
       !  3.) cloud top is within updraft source layer,or
       !  4.) cloud-top detrainment layer begins within
       !  updraft source layer.
       !
       if(ltop<=klcl.or.ltop<=kpbl.or.let + 1<=kpbl) then  
          cldhgt(lc) = 0.0_realkind  
          !  if this is selected shallow source and still does not make a
          !  significant cloud,go on to next grid point.
          if(lc==nchm) goto 325  
          !  if all layers in the specified portion of lower atmosph have been
          !  checked,then
          if(kmix>llfc) then  
             !  if no possible shallow cloud layers were found,goto next grid
             !  point
             if(ishall==0) then  
                goto 325  
             else  
                !  if some potential shallow cloud source layers were found,
                !  find the one that gives the tallest cloud
                goto 68  
             endif
          else  
             !  if there are more layers to check,reset cloud characteristics,
             !  go to next level
             goto 67  
          endif
          !  if this layer has been selected as a shallow convective source,
          !  allow shallow convection(even if the shallow-cloud parameters
          !  allow cldhgt to exceed chmin)
          !     gj all tests for shallow convection done and we have a shallow clo
          !     gj continue with calcsgoto 69
       elseif(lc==nchm) then  
          goto 69  
          !  if cloud depth is greater than minimum depth criterion,and
          !  this layer has not been marked as a shallow convective source laye
          !  allow deep convection
       elseif(cldhgt(lc) >chmin.and.abe>1._realkind) then  
          ishall = 0  
          goto 69  
       endif
       !  at this point,we have a parcel that is able to maintain upward mo
       !  for at least a short distance,but the cloud from this source laye
       !  is not deep enough for "deep"(precipitating) convection;
       !  set ishall=1 to save lc as a possible for source
       !  layer for shallow convection
       !
       !  to disallow shallow convection,comment out next line !!!!!!!!
       ishall = 1  
       !     gj
       !     gj   points where shallow conv occruing saved as value 1
       shal_cgj(i) = 1.0_realkind  
       !     gj
       if(kmix>llfc) then  
          goto 68  
       else  
          goto 67  
       endif
68     continue  
       if(ishall==0) then  
          goto 325  
       else  
          chmax = 0.0_realkind  
          !     gj        chmax = 300.0_realkind
          nchm = 0  
          do  nk = 1,llfc  
             if(cldhgt(nk) >chmax) then  
                nchm = nk  
                chmax = cldhgt(nk)  
             endif
          enddo
          kmix = nchm  
       endif
67     continue  
       !loop goes from klcl-1-->top of convecti
       do nk = k,ltop  
          umf(nk) = 0.0_realkind  
          udr(nk) = 0.0_realkind  
          uer(nk) = 0.0_realkind  
          detlq(nk) = 0.0_realkind  
          detic(nk) = 0.0_realkind  
          pptliq(nk) = 0.0_realkind  
          pptice(nk) = 0.0_realkind  
       enddo
       goto 25  
69     continue  
       !     gj      if(ishall==1)let=klcl
       !     gjasiatest
       kf_ind(i) = 1  
       !     gjasiatest
       if(ishall==1) then  
          kstart = max0(kpbl,klcl)  
          !     gjkfarm preventing let=kstart means shallow convective clouds
          !     gjkfarm retain the entrainment-detrainment terms calculated above
          !     gjkfarm and now does not get treated as a shedding plume from
          !     gjkfarm kpbl to ltop as would occur in loop 75 below.
          !     gjkfarm
          !     gjeurocssheeding plume
          let = kstart  
          !     gjeurocssheeding plume
       endif
       !     gj new terminal detrainment
       !     gj
       if(let<-1) then  
          !     gjreal       if(let==ltop)then      !only deep convection
          dztt = 0.0_realkind  
          dptt = 0.0_realkind  
          !
          ldetst = 0
          if(thelev==0.and.ml==0) then  
             ldetst = let  
          elseif(thelev==0.and.ml>0) then  
             ldetst = min(let,ml)  
          elseif(ml==0.and.thelev>0) then  
             ldetst = min(let,thelev)  
          else  
             ldetst = min(let,min(thelev,ml))  
          endif
          !
          ldetst = max(klcl,ldetst)  
          !     ldetst=max(klcl,10)
          !
          do nj = ldetst,ltop  
             dztt = dztt + dza(nj)  
             dptt = dptt + dp(nj)  
          enddo
          dumfdz = umf(ldetst)/dztt  
          dumfdp = umf(ldetst)/dptt  
          !     gj new terminal detrainment
          !     gj
          !loop goes klcl+1(k+2) to ltop
          do nk = ldetst,ltop  
             if(nk==ltop) then  
                udr(nk) = umf(nk-1)  
                uer(nk) = 0.0_realkind  
                umf(nk) = umf(nk-1)-udr(nk)  
                detlq(nk) = udr(nk)*qliq(nk)*upnew/upold  
                detic(nk) = udr(nk)*qice(nk)*upnew/upold  
                goto 85  
             else  
                udr(nk) = dza(nk)*dumfdz  
                !     gj          udr(nk)=dp(nk)*dumfdp
                umf(nk) = umf(nk-1)-udr(nk)  
                uer(nk) = umf(nk)*(dilfrc(nk)-1.0_realkind)  
                umf(nk) = umf(nk) + uer(nk)  
                detlq(nk) = udr(nk)*qliq(nk)*upnew/upold  
                detic(nk) = udr(nk)*qice(nk)*upnew/upold  
             endif
          enddo
       endif
       !     gjnewdet
       !
       !  if the let and ltop are the same,detrain all of the updraft mass
       !     this level
       !
       if(let==ltop) then  
          udr(ltop) = umf(ltop) + udr(ltop)-uer(ltop)  
          detlq(ltop) = qliq(ltop)*udr(ltop)*upnew/upold  

          detic(ltop) = qice(ltop)*udr(ltop)*upnew/upold  
          !  reverse the mod that allows feedback of rain/snow that originates
          !  detrained air12/9/98jsk
          !     trppt=trppt-(pptliq(ltop)+pptice(ltop))
          !
          !     pptliq(ltop)=0.0_realkind
          !     pptice(ltop)=0.0_realkind
          !
          !
          uer(ltop) = 0.0_realkind  
          umf(ltop) = 0.0_realkind  
          goto 85  
       endif
       !
       !     begin total detrainment at the level above the let
       !     gj
       !     gj  this is where jkain refers too shallow convective cloud as
       !     gj   being a adetraining plume between klcl(cloud base) to
       !     gj   ltop. ltop is calculated as for a normal 1500m radius
       !     gj   cloudwon't this overestimate cloud top for shallow
       !     gj   clouds and place the detrainemnt over too deep a layer.
       !     gj   if we reduce ltop via rad for shallow clouds,we will
       !     gj   have a shallower cloud layer,detraining over a ashallower
       !     gj   depth
       !
       !     gj   do 71 loop goes from 1 level above level of neutral buoyancy
       !     gj   of updraft to diagnosed convective top
       dptt = 0.0_realkind  
       !     gj290402
       uer_shal = 0.0_realkind  
       !     gjoriginal        dptt_cj(ltop+1)=0.0_realkind!level above ltop
       !     gjoriginal        do nj=ltop,let+1,-1
       !     gjoriginal         dptt_cj(nj)=dptt_cj(nj+1)+dp(nj)
       !     gjoriginal        enddo
       !level above ltop
       dptt_cj(ltop + 1) = 0.0_realkind  
       do nj = ltop,let + 2,-1  
          dptt_cj(nj) = dptt_cj(nj + 1) + dp(nj)  

       enddo
       do  nj = let + 1,ltop  
          dptt = dptt + dp(nj)  
       enddo

       dumfdp = umf(let)/dptt  
       !  adjust mass flux profiles,detrainment rates,and precipitation fa
       !     rates to reflect the linear decrease in mass flx between the let a
       !
       do  nk = let + 1,ltop  
          !
          !  entrainment is allowed at every level except for ltop,so disallow
          !  entrainment at ltop and adjust entrainment rates between let and l
          !  so the the dilution factor due to entyrianment is not changed but
          !  the actual entrainment rate will change due due forced total
          !  detrainment in this layer
          !
          !     gjstdcode          if(nk==ltop)then
          !     gjstdcode            udr(nk) = umf(nk-1)
          !     gjstdcode            uer(nk) = 0.0_realkind
          !     gjstdcode            detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
          !     gjstdcode            detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
          !     gjstdcode          else
          !     gjstdcode            udr(nk)=dp(nk)*dumfdp
          !     gjstdcode            umf(nk)=umf(nk-1)-udr(nk)
          !     gjstdcode            uer(nk)=umf(nk)*(dilfrc(nk)-1.)
          !     gjstdcode            umf(nk)=umf(nk)+uer(nk)
          !     gjstdcode            detlq(nk)=udr(nk)*qliq(nk)*dilfrc(nk)
          !     gjstdcode            detic(nk)=udr(nk)*qice(nk)*dilfrc(nk)
          !     gjstdcode          endif
          if(nk==ltop) then  
             udr(nk) = umf(nk-1)  
             !     gj290402
             umf(nk) = umf(nk-1)-udr(nk)  
             !     gj290402
             uer(nk) = 0.0_realkind  
             detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)  
             detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)  
          else  
             !     gj180602
             if(nk>=let + 2) then  
                udr(nk) = dp(nk)*dumfdp +(uer_shal*dp(nk)/dptt_cj(nk))
                uer_shal = uer_shal-(uer_shal*dp(nk)/dptt_cj(nk))
             else  
                udr(nk) = dp(nk)*dumfdp  
             endif
             !     gj180602
             umf(nk) = umf(nk-1)-udr(nk)  
             uer(nk) = umf(nk)*(dilfrc(nk)-1.0_realkind)  
             !     gj180602
             uer_shal = uer_shal + uer(nk)  
             !     gj180602
             umf(nk) = umf(nk) + uer(nk)  
             detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)  
             detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)  
          endif
          !  reverse the mod that allows feedback of rain/snow that originates
          !  detrained air12/9/98jsk
          !     trppt=trppt-pptliq(nk)-pptice(nk)
          !     pptliq(nk)=(umf(nk-1)-udr(nk))*qlqout(nk)
          !     pptice(nk)=(umf(nk-1)-udr(nk))*qicout(nk)
          !     trppt=trppt+pptliq(nk)+pptice(nk)
          if(nk>=let + 2) then  
             trppt = trppt-pptliq(nk)-pptice(nk)  
             pptliq(nk) = umf(nk-1)*qlqout(nk)  
             pptice(nk) = umf(nk-1)*qicout(nk)  
             trppt = trppt + pptliq(nk) + pptice(nk)  

          endif
       enddo
       !
       !  send updraft characteristics to output files
       !
85     continue  
       !  when moist perturbation is added to updraft and klcl<=kpbl,
       !  reset thetee so that moisture perturbation only affects updrafts..
       if(klcl<=kpbl) then  
          do  nk1 = klcl,kpbl  
             call envirtht(p0(nk1),t0(nk1),q0(nk1),thetee(nk1) &
                  ,0.0_realkind,rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          enddo
       endif
       !
       !  extend the updraft mass flux profile down to the source layer for
       !     updraft airalso,define thetae for levels below the lcl
       !
       do  nk = 1,k  
          if(nk>=lc) then  
             if(nk==lc) then  
                umf(nk) = vmflcl*dp(nk)/dpthmx  
                uer(nk) = vmflcl*dp(nk)/dpthmx  
             elseif(nk<=kpbl) then  
                uer(nk) = vmflcl*dp(nk)/dpthmx  
                umf(nk) = umf(nk-1) + uer(nk)  
             else  
                umf(nk) = vmflcl  
                uer(nk) = 0.0_realkind  
             endif
             tu(nk) = tmix +(z00(nk)-zmix)*gdry  
             qu(nk) = qmix  
             wu(nk) = wlcl  
          else  
             tu(nk) = 0.0_realkind  
             qu(nk) = 0.0_realkind  
             !     gjextend updraft diags below source layer lc.
             !     gj            tu(nk)=t0(nk)
             !     gj            qu(nk)=q0(nk)
             !     gjextend updraft diags below source layer lc.
             umf(nk) = 0.0_realkind  
             wu(nk) = 0.0_realkind  
             uer(nk) = 0.0_realkind  
          endif
          udr(nk) = 0.0_realkind  
          qdt(nk) = 0.0_realkind  
          qliq(nk) = 0.0_realkind  
          qice(nk) = 0.0_realkind  
          qlqout(nk) = 0.0_realkind  
          qicout(nk) = 0.0_realkind  
          pptliq(nk) = 0.0_realkind  
          pptice(nk) = 0.0_realkind  
          detlq(nk) = 0.0_realkind  
          detic(nk) = 0.0_realkind  
          ratio2(nk) = 0.0_realkind  
          ee = q0(nk)*p0(nk)/(epsilo + q0(nk))  
          tlog = log(ee/aliq)  
          tdpt =(cliq-dliq*tlog)/(bliq-tlog)  
          tsat = tdpt-(0.212_realkind + 1.571e-3_realkind*(tdpt-t00)-4.36e-4_realkind*(t0(nk)-t00))*(t0(nk)-tdpt)
          thta = t0(nk)*(1.e5_realkind/p0(nk))**(0.2854_realkind*(1.0_realkind-0.28_realkind*q0(nk)))
          thetee(nk) = thta*exp((3374.6525_realkind/tsat-2.5403_realkind)*q0(nk)*(1.0_realkind + 0.81_realkind*q0(nk)))
          thtes(nk) = thta*exp((3374.6525_realkind/t0(nk)-2.5403_realkind)*qes(nk)*(1.0_realkind + 0.81_realkind*qes(nk)))
          eqfrc(nk) = 1.0_realkind  

       enddo
       ltop1 = ltop + 1  
       ltopm1 = ltop-1  
       !
       !  define variables above cloud top
       !
       do  nk = ltop1,kl  
          umf(nk) = 0.0_realkind  
          udr(nk) = 0.0_realkind  
          uer(nk) = 0.0_realkind  
          qdt(nk) = 0.0_realkind  
          qliq(nk) = 0.0_realkind  
          qice(nk) = 0.0_realkind  
          qlqout(nk) = 0.0_realkind  
          qicout(nk) = 0.0_realkind  
          detlq(nk) = 0.0_realkind  
          detic(nk) = 0.0_realkind  
          pptliq(nk) = 0.0_realkind  
          pptice(nk) = 0.0_realkind  
          if(nk>ltop1) then  
             tu(nk) = 0.0_realkind  
             qu(nk) = 0.0_realkind  
             wu(nk) = 0.0_realkind  
          endif
          thta0(nk) = 0.0_realkind  
          thtau(nk) = 0.0_realkind  
          ems(nk) = 0.0_realkind  
          emsd(nk) = 0.0_realkind  
          tg(nk) = t0(nk)  
          qg(nk) = q0(nk)  
          qlg(nk) = 0.0_realkind  
          qig(nk) = 0.0_realkind  
          qrg(nk) = 0.0_realkind  
          qsg(nk) = 0.0_realkind  
          omg(nk) = 0.0_realkind  
       enddo
       omg(kl + 1) = 0.0_realkind  
       p150 = p0(klcl)-1.50e4_realkind  
       do nk = 1,ltop  
          ems(nk) = dp(nk)*dxsq/g  
          emsd(nk) = 1.0_realkind/ems(nk)  
          !
          !  initialize some variables to be used later in the vert advection s
          !
          exn(nk) =(p00/p0(nk))**(0.2854_realkind*(1.0_realkind-0.28_realkind*qdt(nk)))
          thtau(nk) = tu(nk)*exn(nk)  
          exn(nk) =(p00/p0(nk))**(0.2854_realkind*(1.0_realkind-0.28_realkind*q0(nk)))
          thta0(nk) = t0(nk)*exn(nk)  
          !     gjentnew
          ddilfrc(nk) = 1.0_realkind/dilfrc2(nk)  
          !     gjentnew
          !
          !  lvf is the level at which moisture flux is estimated as the basis
          !  precipitation efficiency calculations
          !
          if(p0(nk) >p150) lvf = nk  
          omg(nk) = 0.0_realkind  
       enddo
       ! jsk mods
       lvf = min0(lvf,let)  
       usr = umf(lvf + 1)*(qu(lvf + 1) + qliq(lvf + 1) + qice(lvf + 1))
       usr = min(usr,trppt)  
       !
       !  compute convective time scale(timec). the mean wind at the lcl
       !  and midtroposphere is used.
       !
       wspd(klcl) = sqrt(u0(klcl)*u0(klcl) + v0(klcl)*v0(klcl))
       wspd(l5) = sqrt(u0(l5)*u0(l5) + v0(l5)*v0(l5))  
       wspd(ltop) = sqrt(u0(ltop)*u0(ltop) + v0(ltop)*v0(ltop))
       vconv = 0.5_realkind*(wspd(klcl) + wspd(l5))  
       timec = dx/vconv  
       tadvec = timec  
       timec = max(1800.0_realkind,timec)  
       timec = min(3600.0_realkind,timec)  

       if(ishall==1) timec = 2400.0_realkind  
       nic = nint(timec/(0.5_realkind*dt2))  
       timec = real(nic,realkind)*0.5_realkind*dt2  
       !
       !  compute wind shear and precipitation efficiency.
       !
       !     shsign = cvmgt(1.,-1.,wspd(ltop)>wspd(klcl))
       if(wspd(ltop) >wspd(klcl)) then  
          shsign = 1.0_realkind  
       else  
          shsign =-1.0_realkind  
       endif
       vws =(u0(ltop)-u0(klcl))*(u0(ltop)-u0(klcl)) &
            +(v0(ltop)-v0(klcl))*(v0(ltop)-v0(klcl))
       vws = 1.e3_realkind*shsign*sqrt(vws)/(z00(ltop)-z00(lcl))  
       pef = 1.591_realkind + vws*(-0.639_realkind + vws*(9.53e-2_realkind-vws*4.96e-3_realkind) )
       pef = max(pef,0.2_realkind)  
       pef = min(pef,0.9_realkind)  
       !
       !  precipitation efficiency is a function of the height of cloud base
       !
       cbh =(zlcl-z00(1))*3.281e-3_realkind  
       if(cbh<3.0_realkind) then  
          rcbh = 0.02_realkind  
       else  
          rcbh = 0.96729352_realkind + cbh*(-0.70034167_realkind + cbh*(0.162179896_realkind + &
               cbh*(-1.2569798e-2_realkind + cbh*(4.2772e-4_realkind-cbh*5.44e-6_realkind)) &
               ))
       endif
       if(cbh>25._realkind) rcbh = 2.4_realkind  
       pefcbh = 1.0_realkind/(1.0_realkind + rcbh)  
       pefcbh = min(pefcbh,0.9_realkind)  
       !
       !  mean pef. is used to compute rainfall.
       !
       peff = 0.5_realkind*(pef + pefcbh)  
       ! jsk mods
       peff2 = peff  
       !
       !    
       !     compute downdraft properties                 
       !    
       !
       !
       !
       !     gj180100  no downdraft for shallow convection,but for
       !     gk180100  safety zero out downdraft arrays here.
       tder = 0.0_realkind  
       if(ishall==1) then  
          lfs = 1  
          do  nk = 1,kl  
             dmf(nk) = 0.0_realkind  
             der(nk) = 0.0_realkind  
             ddr(nk) = 0.0_realkind  
             wd(nk) = 0.0_realkind  
             tz(nk) = 0.0_realkind  
             qd(nk) = 0.0_realkind  
             thtad(nk) = 0.0_realkind  
          enddo
          !skip downdraft calcs
          goto 141  
       endif
       !     gj180100
       !
       !  start downdraft about 150 mb above cloud base
       !
       !     kstart=max0(kpbl,klcl)
       kstart = kpbl  
       do  nk = kstart + 1,kl  
          dppp = p0(kstart)-p0(nk)  
          !     if(dppp>200.e2)then
          if(dppp>150.0e2_realkind) then  
             klfs = nk  
             goto 684  
          endif
       enddo
684    continue  
       klfs = min0(klfs,let-1)  
       lfs = klfs  
       !     gjdd
       !     gjnewlfs
       klfs = min0(thelev,let-1)  
       lfs = klfs  
       !     gjnewlfs
       !     gjdd
       !     print*,'lfs =',lfs
       !
       !  if lfs is not at least 50 mb above cloud base(implying that the
       !  level of equil temp,let,is just above cloud base) do not allow a
       !  downdraft
       !
       if((p0(kstart)-p0(lfs)) <50.0e2_realkind) then  
          !     if((p0(klcl)-p0(lfs))<50.e2)then
          tder = 0.0_realkind  
          !     print*,'p0(kstart),p0(lfs) =',p0(kstart),p0(lfs)
          !     write(98,*)'downdraft not allowedtoo shallow'
          goto 141  

       endif

       theted(lfs) = thetee(lfs)  

       qd(lfs) = q0(lfs)  
       !     gjmosierr
       !
       !  call tpmix2dd to find wet-bulb temp,qv
       !
       call tpmix2dd(p0(lfs),theted(lfs),tz(lfs),qss)  
       thtad(lfs) = tz(lfs)*(p00/p0(lfs))**(0.2854_realkind*(1.0_realkind-0.28_realkind*qss))
       !
       !  take a first guess at the initial downdraft mass flux
       !
       tvd(lfs) = tz(lfs)*(1.0_realkind + 0.608_realkind*qss)  
       rdd = p0(lfs)/(r*tvd(lfs))  
       a1 =(1.0_realkind-peff)*au0  
       dmf(lfs) =-a1*rdd  
       der(lfs) = dmf(lfs)  
       ddr(lfs) = 0.0_realkind  
       rhbar = rh(lfs)*dp(lfs)  
       dptt = dp(lfs)  
       do  nd = lfs-1,kstart,-1  
          nd1 = nd+1  
          der(nd) = der(lfs)*ems(nd)/ems(lfs)  
          ddr(nd) = 0.0_realkind  
          dmf(nd) = dmf(nd1) + der(nd)  
          theted(nd) =(theted(nd1)*dmf(nd1) + thetee(nd)*der(nd))/dmf(nd)
          qd(nd) =(qd(nd1)*dmf(nd1) + q0(nd)*der(nd))/dmf(nd)
          dptt = dptt + dp(nd)  
          rhbar = rhbar + rh(nd)*dp(nd)  
       enddo
       rhbar = rhbar/dptt  
       !     print*,' rhbar =',rhbar
       dmffrc = 2.0_realkind*(1.0_realkind-rhbar)  
       dpdd = 0.0_realkind  
       !     ldt = klcl-1
       !     do 138 nd = klcl-1,1,-1
       !
       !  3/2/98:  calculate melting effect
       !  first,compute total frozen precipitation generated
       !
       pptmlt = 0.0_realkind  
       do  nk = klcl,ltop  
          pptmlt = pptmlt + pptice(nk)  
       enddo
       if(lc<ml) then  
          !     dtmelt = rlf*pptmlt/(cp*dmffrc*umf(klcl))
          !
          !  for now,calculate melting effect as if dmf =-umf at klcl,i.e.,
          !  if dmffrc=1.  otherwise,for small dmffrc,dtmelt gets too large!
          !  12/14/98 jsk
          dtmelt = rlf*pptmlt/(cp*umf(klcl))  
       else  
          dtmelt = 0.0_realkind  
       endif
       !     print*,'dtmelt =',dtmelt
       !     ldt = min0(lfs-1,klcl-1)
       ldt = min0(lfs-1,kstart-1)  
       !     call tpmix2dd(p0(klcl),theted(klcl),tz(klcl),qss)
       call tpmix2dd(p0(kstart),theted(kstart),tz(kstart),qss)
       !     print*,'p(klcl),tz before melt =',p0(klcl)/100.,tz(klcl)
       !     tz(klcl) = tz(klcl)-dtmelt
       tz(kstart) = tz(kstart)-dtmelt  
       !     print*,'p(klcl),tz after melt =',p0(klcl)/100.0_realkind,tz(klcl)
       !     es=aliq*exp((bliq*tz(klcl)-cliq)/(tz(klcl)-dliq))
       !     qss=0.622*es/(p0(klcl)-es)
       !     gjesat
       !     gjesat          es=esat(tz(kstart))
       !     gjesat          qss=0.622*es/(p0(kstart)-(0.378*es))
       es = aliq*exp((bliq*tz(kstart)-cliq)/(tz(kstart) -dliq))
       qss = epsilo*es/(p0(kstart)-es)  
       !     qd(klcl)=qs
       !     theted(klcl)=tz(klcl)*(1.e5/p0(klcl))**(0.2854*(1.-0.28*qss))*
       !     +           exp((3374.6525/tz(klcl)-2.5403)*qss*(1.+0.81*qss))
       theted(kstart) = tz(kstart)*(1.e5_realkind/p0(kstart))**( &
            0.2854_realkind*(1.0_realkind-0.28_realkind*qss))*exp((3374.6525_realkind/tz(kstart) &
            -2.5403_realkind)*qss*(1.0_realkind + 0.81_realkind*qss))
       !.
       !     gj3dppt
       !     gj3dppt top of the part of downdraft that actually evaporates conv
       !     gj3dppt
       !     ldt = min0(lfs-1,klcl-1)
       ldt = min0(lfs-1,kstart-1)  
       do  nd = ldt,1,-1  
          dpdd = dpdd+dp(nd)  
          !     theted(nd) = theted(klcl)
          !     qd(nd)     = qd(klcl)
          theted(nd) = theted(kstart)  
          qd(nd) = qd(kstart)  
          !
          !  call tpmix2dd to find wet bulb temp,saturation mixing ratio
          !
          call tpmix2dd(p0(nd),theted(nd),tz(nd),qss)  
          qsd(nd) = qss  
          !
          !  specify rh decrease of 10%/km in downdraft
          !
          !     rhh = 1.-0.1/1000.*(z00(klcl)-z00(nd))
          !     rhh = 1.-0.1/1000.*(z00(kstart)-z00(nd))
          rhh = 1.0_realkind-0.2_realkind/1000.0_realkind*(z00(kstart)-z00(nd))  
          !     gjtest
          !     gj         rhh = 1.0-0.2/1000.*(z00(kstart)-z00(nd))
          !     gj         rhh = 0.8-0.2/1000.*(z00(kstart)-z00(nd))
          rhh = max(0.0_realkind,rhh)  
          !
          !  adjust downdraft temp,q to specified rh:
          !
          if(rhh<1.0_realkind) then  
             dssdt =(cliq-bliq*dliq)/((tz(nd)-dliq) *(tz(nd)-dliq))
             rl = xlv0-xlv1*tz(nd)  
             dtmp = rl*qss*(1.0_realkind-rhh)/(cp + rl*rhh*qss*dssdt)
             t1rh = tz(nd) + dtmp  
             !     gjesat
             !     gjesat           es=esat(t1rh)
             !     gjesat           qsrh=0.622*es/(p0(nd)-(0.378*es))
             es = rhh*aliq*exp((bliq*t1rh-cliq)/(t1rh-dliq))
             qsrh = epsilo*es/(p0(nd)-es)  
             !
             !  check to see if mixing ratio at specified rh is less than actual
             !  mixing ratioif so,adjust to give zero evaporation
             !
             if(qsrh<qd(nd)) then  
                qsrh = qd(nd)  
                t1rh = tz(nd) +(qss-qsrh)*rl/cp  
                !     t1rh=tz(nd)
             endif
             tz(nd) = t1rh  
             qss = qsrh  
             qsd(nd) = qss  
          endif
          tvd(nd) = tz(nd)*(1.0_realkind + 0.608_realkind*qsd(nd))  
          if(tvd(nd) >tv0(nd) .or.nd==1) then  
             !     if(tvd(nd)>tv0(nd).or.nd==lc)then
             ldb = nd  
             goto 139  
          endif
       enddo
139    continue  
       !     print*,'lfs,ldb =',p0(lfs),p0(ldb)
       ! no downdraft allowed!
       if((p0(ldb)-p0(lfs)) <50.0e2_realkind) then  

          tder = 0.0_realkind  
          !     write(98,*)'downdraft not allowedtoo shallow'
          !     write(6,*)'downdraft not allowedtoo shallow'
          !     print*,'p0(ldb),p0(lfs) =',p0(ldb),p0(lfs)
          goto 141  
       endif
       tder = 0.0_realkind  
       !
       !  calculate an evaporation rate for given mass flux
       !
       !     print*,'lfs,ldb,ldt,i,j =',p0(lfs),p0(ldb),p0(ldt)
       !    *,i,j
       !     if(i==96 .and. j==44)then
       !     write(15,4433) istid,rlat,rlon,idate(3)-1900,
       !     & idate(1),idate(2),ihrst,ifhr
       !     write(15,4422) lmh+1
       !     write(96,4422) i,j,kl
       !     do 550 nl=1,kl
       !     l = kl-nl+1
       !     write(96,4455) p0(l)*0.01,t0(l)-273.16,q0(l)*1000.,
       !    *                 u0(l),v0(l),w0avg(i,j,nl),dp(l),q2(i,j,nl)
       !     550   continue
       !     endif
       !44   33  format(i6,2f8.2,i5,4i3)
       !44   22  format(i6)
       !44   55  format(8f12.3)
       do  nd = ldt,ldb,-1  
          nd1 = nd+1  
          ddr(nd) =-dmf(kstart)*dp(nd)/dpdd  
          der(nd) = 0.0_realkind  
          dmf(nd) = dmf(nd1) + ddr(nd)  
          tder = tder +(qsd(nd)-qd(nd))*ddr(nd)  
          qd(nd) = qsd(nd)  
          thtad(nd) = tz(nd)*(p00/p0(nd))**(0.2854_realkind*&
               (1.0_realkind-0.28_realkind*qd(nd)))
          !
       enddo
       !     gjnoddevap
       !     gj          tder=0.0_realkind
       !     gjnoddevap
       !
       !  if downdraft does not evaporate any water for specified relative
       !  humidity,no downdraft is allowed
       !
141    if(tder<1.0_realkind) then  
          pptflx = trppt  
          cpr = trppt  
          tder = 0.0_realkind  
          cndtnf = 0.0_realkind  
          updinc = 1.0_realkind  
          ldb = lfs  

          ldt = lfs  
          do ndk = 1,ltop  
             dmf(ndk) = 0.0_realkind  
             der(ndk) = 0.0_realkind  
             ddr(ndk) = 0.0_realkind  
             thtad(ndk) = 0.0_realkind  
             wd(ndk) = 0.0_realkind  
             tz(ndk) = 0.0_realkind  
             qd(ndk) = 0.0_realkind  
          enddo
          aincm2 = 100.0_realkind  
          goto 165  

       endif
       !
       !  adjust downdraft mass flux so that evaporation rate in downdraft i
       !  consistent with precipitation efficiency relationship
       !
       !  ppr is the total amount of precipitation that falls  out of the up
       !  from cloud base to the lfsupdraft mass flux will be increased u
       !  the lfs to account for updraft air mixing with environmental air t
       !  the updraft,so ppr will increase proportionately
       !
       !  cndtnf is the amount of condensate transferred along with updraft
       !  the downdraft at the lfs
       !
       !  ddinc is the factor by which to increase the first-guess downdraft
       !  flux to satisfy the precip efficiency relationship,updinc is the
       !  which to increase the updraft mass flux below the lfs to account f
       !  transfer of mass from updraft to downdraft
       ddinc =-dmffrc*umf(klcl)/dmf(kstart)  
       updinc = 1.0_realkind  
       if(tder*ddinc>trppt) then  
          ddinc = trppt/tder  
       endif
       tder = tder*ddinc  
       do  nk = ldb,lfs  
          dmf(nk) = dmf(nk)*ddinc  
          der(nk) = der(nk)*ddinc  
          ddr(nk) = ddr(nk)*ddinc  
       enddo
       !     cpr=trppt+ppr*(updinc-1.)
       cpr = trppt  
       pptflx = trppt-tder  
       !     gjnotder        pptflx=trppt
       !     gjeurocstest
       !     print*,'pptflx = ',pptflx
       !     pptflx=pptflx+peff*ppr*(updinc-1.)
       !     peff=peff2                                ! jsk mods

       peff = pptflx/trppt  
       !     write(*,*)'precip efficiency =',peff
       !     if(iprnt)then
       !     write(98,*)'precip efficiency =',peff
       !     endif
       !     tder=tder*ddinc
       !
       !  adjust updraft mass flux,mass detrainment rate,and liquid water
       !     detrainment rates to be consistent with the transfer of the estima
       !     from the updraft to the downdraft at the lfs
       !
       do nk = lc,lfs  
          umf(nk) = umf(nk)*updinc  
          udr(nk) = udr(nk)*updinc  
          uer(nk) = uer(nk)*updinc  
          pptliq(nk) = pptliq(nk)*updinc  
          pptice(nk) = pptice(nk)*updinc  
          detlq(nk) = detlq(nk)*updinc  
          detic(nk) = detic(nk)*updinc  
       enddo
       !     gj071299
165    continue  
       !     gj071299
       !
       !  zero out the arrays for downdraft data at levels above and below t
       !  downdraft
       !
       !     gj180100  don't try this for shallow convective points downdraft
       !     gj180100  arrays are zeroed out higher up
       if(ishall/=1) then  
          if(ldb>1) then  
             do nk = 1,ldb-1  
                dmf(nk) = 0.0_realkind  
                der(nk) = 0.0_realkind  
                ddr(nk) = 0.0_realkind  
                wd(nk) = 0.0_realkind  
                tz(nk) = 0.0_realkind  
                qd(nk) = 0.0_realkind  
                thtad(nk) = 0.0_realkind  
             enddo
          endif
          do nk = lfs + 1,kl  
             dmf(nk) = 0.0_realkind  
             der(nk) = 0.0_realkind  
             ddr(nk) = 0.0_realkind  
             wd(nk) = 0.0_realkind  
             tz(nk) = 0.0_realkind  
             qd(nk) = 0.0_realkind  
             thtad(nk) = 0.0_realkind  
          enddo
          do  nk = ldt + 1,lfs-1  
             tz(nk) = 0.0_realkind  
             qd(nk) = 0.0_realkind  
             ! jsk mods
             thtad(nk) = 0.0_realkind  
          enddo

       endif
       !  set limits on the updraft and downdraft mass fluxes so that the in
       !     into convective drafts from a given layer is no more than is avail
       !     in that layer initially
       !
       !     gj071299  165   aincmx=1000.0_realkind
       aincmx = 1000.0_realkind  
       lmax = max0(klcl,lfs)  
       do  nk = lc,lmax  
          if((uer(nk)-der(nk)) >0.0_realkind)then
             aincm1 = ems(nk)/((uer(nk)-der(nk))*timec)
          endif
          aincmx = min(aincmx,aincm1)  
       enddo
       ainc = 1.0_realkind  
       if(aincmx<ainc) ainc = aincmx  
       !
       !  save the relevent variables for a unit updrft and downdrftthey
       !  iteratively adjusted by the factor ainc to satisfy the stabilizati
       !  closure
       !
       ncount = 0  
       tder2 = tder  
       pptfl2 = pptflx  
       do nk = 1,ltop  
          detlq2(nk) = detlq(nk)  
          detic2(nk) = detic(nk)  
          udr2(nk) = udr(nk)  
          uer2(nk) = uer(nk)  
          ddr2(nk) = ddr(nk)  
          der2(nk) = der(nk)  
          umf2(nk) = umf(nk)  
          dmf2(nk) = dmf(nk)  
       enddo
       fabe = 1.0_realkind  
       stab = 0.95_realkind  
       !     if(xnin<5)stab=stab*(xnin-1.)*0.25
       noitr = 0  
       istop = 0  
       if(ainc/aincmx>0.999_realkind) then  
          ncount = 0  
          goto 255  
       endif
       !  shallow convection mods
       if(ishall==1) then  
          !
          !  find the maximum tke value between lc and klcl
          !     gj and set up the ainc term used to detrmine convergence of
          !     gj closure requrirements for deep convection
          !     gj
          !     gj   ainc=[1/20.*tkemax*(delp of layer)*(delx*delx)]/
          !     gj                   [(cloud base mass flux)*grav*2400secs)]
          !     gj
          tkemax = 0.0_realkind  
          do  k = lc,klcl  
             nk = kl-k + 1  
             tkemax = max(tkemax,cbr_tke(i,nk))  
          enddo
          !
          tkemax = min(tkemax,10.0_realkind)  
          tkemax = max(tkemax,5.0_realkind)  
          !     tkemax = 10.0_realkind
          !  3_24_99dpmin was changed for shallow convection so that it is t
          !  the same as for deep convection(5.e3).  since this doubles
          !(roughly) the value of dpthmx,add a factor of 0.5 to calcu-
          !  lation of evac
          !     gj06921         evac  = tkemax*0.1
          evac = 0.5_realkind*tkemax*0.1_realkind  
          !     ainc = 0.1*dpthmx*dxij*dxij/(vmflcl*g*timec)
          ainc = evac*dpthmx*dx*dx/(vmflcl*g*timec)  
          goto 255  
          !     gj
          !     gj   we have ainc now go to 255 continue make the test for
          !     gj   ainc convergence and where necesaary go back up to 175 contin
          !     gj   for another go at the ainc calculation
       endif
       !  end of shallow mods
175    ncount = ncount + 1  

       !     compute properties for compensational subsidence    
       !
       !  determine omega value necessary at top and bottom of each layer to
       !  satisfy mass continuity
       !

       dtt = timec  
       do nk = 1,ltop  
          domgdp(nk) =-(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
          if(nk>1) then  
             omg(nk) = omg(nk-1)-dp(nk-1)*domgdp(nk-1)  
             dtt1 = 0.75_realkind*dp(nk-1)/(abs(omg(nk)) + 1.0e-10_realkind)  
             dtt = min(dtt,dtt1)  
          endif
       enddo
       do  nk = 1,ltop  
          thpa(nk) = thta0(nk)  
          qpa(nk) = q0(nk)  
          nstep = nint(timec/dtt + 1._realkind)  
          dtime = timec/real(nstep,realkind)  
          fxm(nk) = omg(nk)*dxsq/g  
       enddo
       !
       !  do an upstream/forward-in-time advection of theta,qv
       !
       do ntc = 1,nstep  
          !
          !  assign theta and q values at the top and bottom of each layer base
          !  sign of omega
          !
          do  nk = 1,ltop  
             thfxin(nk) = 0.0_realkind  
             thfxout(nk) = 0.0_realkind  
             qfxin(nk) = 0.0_realkind  
             qfxout(nk) = 0.0_realkind  
          enddo
          do nk = 2,ltop  
             if(omg(nk) <=0.0_realkind) then  
                thfxin(nk) =-fxm(nk)*thpa(nk-1)  
                qfxin(nk) =-fxm(nk)*qpa(nk-1)  
                thfxout(nk-1) = thfxout(nk-1) + thfxin(nk)  
                qfxout(nk-1) = qfxout(nk-1) + qfxin(nk)  
             else  
                thfxout(nk) = fxm(nk)*thpa(nk)  
                qfxout(nk) = fxm(nk)*qpa(nk)  
                thfxin(nk-1) = thfxin(nk-1) + thfxout(nk)  
                qfxin(nk-1) = qfxin(nk-1) + qfxout(nk)  
             endif
          enddo
          !
          !  update the theta and qv values at each level
          !
          do nk = 1,ltop  
             thpa(nk)=thpa(nk)+(thfxin(nk)+udr(nk)*thtau(nk)+ddr(nk)*&
                  thtad(nk)-thfxout(nk)-(uer(nk)-der(nk))*thta0(nk))*dtime*emsd(nk)
             if(nk>=lc.and.nk<=kpbl) then
                qpa(nk)=qpa(nk)+(qfxin(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk) - &
                     qfxout(nk)-uer(nk)*quer(nk)+der(nk)*q0(nk))*dtime*emsd(nk) 
             else  
                qpa(nk)=qpa(nk)+(qfxin(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk)-&
                     qfxout(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)
             endif
          enddo
       enddo
       !     gj
       !     gj   extract the individual terms making up the convective tendenc
       !     gj   presently do it just for temperature
       !     gj
       do nk = 1,ltop  
          thtag(nk) = thpa(nk)  
          qg(nk) = qpa(nk)  
       enddo
       !
       !  check to see if mixing ratio dips below zero anywhere;  if so,bor
       !  moisture from adjacent layers to bring it back up above zero
       !
       do nk = 1,ltop  
          if(qg(nk) <0.0_realkind) then  
             ! jsk mods
             if(nk==1) then  
                ! jsk mods
                print*,'!!!!! problem with kf scheme:  '  
                ! jsk mods
                print*,'qg = 0 at the surface!!!!!!!'  
             endif
             nk1 = nk + 1  
             if(nk==ltop) nk1 = klcl  
             tma = qg(nk1)*ems(nk1)  
             tmb = qg(nk-1)*ems(nk-1)  
             tmm =(qg(nk)-1.0e-9_realkind)*ems(nk)  
             bcoeff =-tmm/((tma*tma)/tmb + tmb)  
             acoeff = bcoeff*tma/tmb  
             tmb = tmb*(1.0_realkind-bcoeff)  
             tma = tma*(1.0_realkind-acoeff)  
             if(nk==ltop) then  
                qvdiff =(qg(nk1)-tma*emsd(nk1))*100.0_realkind/qg(nk1)
                if(abs(qvdiff) >1.0_realkind) then  
                   print*,'!!!warning!!! cloud base water vapor chang'
                endif
             endif
             qg(nk) = 1.0e-9_realkind
             qg(nk1) = tma*emsd(nk1)  
             qg(nk-1) = tmb*emsd(nk-1)  
          endif
       enddo
       topomg =(udr(ltop)-uer(ltop))*dp(ltop)*emsd(ltop)  
       if(abs(topomg-omg(ltop)) >1.0e-3_realkind) then  
          print *, 'error:  mass does not balance in kf scheme; &
               & topomg,omg =',topomg,omg(ltop)
          print *, 'error:  mass does not balance in kf scheme;      t &
               &opomg,omg =',topomg,omg(ltop)
          istop = 1  
          iprnt = .true.  
          goto 265  
       endif
       
       !  convert theta to t
       
       do  nk = 1,ltop  
          exn(nk) =(p00/p0(nk))**(0.2854_realkind*(1.0_realkind-0.28_realkind*qg(nk)))
          tg(nk) = thtag(nk)/exn(nk)  
          tvg(nk) = tg(nk)*(1.0_realkind + 0.608_realkind*qg(nk))  
       enddo
       !  shallow mods
       if(ishall==1) then  
          goto 265  
       endif
       !
       !
       !    
       !  compute new cloud and change in available buoyant energy.   *
       !    
       !
       !
       !  the following computations are similar to that for updraft
       !
       thmix = 0.0_realkind  
       qmix = 0.0_realkind  
       pmix = 0.0_realkind  
       do  nk = lc,kpbl  
          rocpq = 0.2854_realkind*(1.0_realkind-0.28_realkind*qg(nk))  
          thmix = thmix + dp(nk)*tg(nk)*(p00/p0(nk))**rocpq
          qmix = qmix + dp(nk)*qg(nk)  
          pmix = pmix + dp(nk)*p0(nk)  
       enddo
       thmix = thmix/dpthmx  
       qmix = qmix/dpthmx  
       pmix = pmix/dpthmx  
       rocpq = 0.2854_realkind*(1.0_realkind-0.28_realkind*qmix)  
       tmix = thmix*(pmix/p00)**rocpq  
       !     gjesat
       !     gjesat        es=esat(tmix)
       !     gjesat        qss=0.622*es/(pmix-(0.378*es))
       es = aliq*exp((tmix*bliq-cliq)/(tmix-dliq))  
       qss = epsilo*es/(pmix-es)  
       !
       !  remove supersaturation for diagnostic purposes,if necessary
       !
       if(qmix>qss) then  
          rl = xlv0-xlv1*tmix  
          cpm = cp*(1.0_realkind + 0.887_realkind*qmix)  
          dssdt = qss*(cliq-bliq*dliq)/((tmix-dliq)*(tmix-dliq))
          dq =(qmix-qss)/(1.0_realkind + rl*dssdt/cpm)  
          tmix = tmix + rl/cp*dq  
          qmix = qmix-dq  
          rocpq = 0.2854_realkind*(1.0_realkind-0.28_realkind*qmix)  
          thmix = tmix*(p00/pmix)**rocpq  
          tlcl = tmix  
          plcl = pmix  
       else  
          qmix = max(qmix,0.0_realkind)  
          emix = qmix*pmix/(epsilo + qmix)  
          tlog = log(emix/aliq)  
          tdpt =(cliq-dliq*tlog)/(bliq-tlog)  
          tlcl = tdpt-(0.212_realkind + 1.571e-3_realkind*(tdpt-t00)-4.36e-4_realkind*(tmix-t00))*(tmix-tdpt)
          tlcl = min(tlcl,tmix)  
          cporq = 1.0_realkind/rocpq  
          plcl = p00*(tlcl/thmix)**cporq  
       endif
       tvlcl = tlcl*(1.0_realkind + 0.608_realkind*qmix)  
       do  nk = lc,kl  
          klcl = nk  
          if(plcl>=p0(nk)) goto 240  
       enddo
240    k = klcl-1  
       dlp = log(plcl/p0(k))/log(p0(klcl)/p0(k))  
       !
       !  estimate environmental temperature and mixing ratio at the lcl
       !
       tenv = tg(k) +(tg(klcl)-tg(k))*dlp  
       qenv = qg(k) +(qg(klcl)-qg(k))*dlp  
       tven = tenv*(1.0_realkind + 0.608_realkind*qenv)  
       !     tvbar=0.5*(tvg(k)+tven)
       !     zlcl=z00(k)+r*tvbar*log(p0(k)/plcl)/g
       zlcl = z00(k) +(z00(klcl)-z00(k))*dlp  
       tvavg = 0.5_realkind*(tven + tg(klcl)*(1.0_realkind + 0.608_realkind*qg(klcl)))  
       plcl = p0(klcl)*exp(g/(r*tvavg)*(z00(klcl)-zlcl))  
       theteu(k) = tmix*(1.0e5_realkind/pmix)**(0.2854_realkind*(1.0_realkind-0.28_realkind*qmix))*&
            exp((3374.6525_realkind/tlcl-2.5403_realkind)*qmix*(1.0_realkind + 0.81_realkind*qmix))
       !     gjesat
       !     gjesat          es=esat(tenv)
       !     gjesat          qese=0.622*es/(plcl-(0.378*es))
       es = aliq*exp((tenv*bliq-cliq)/(tenv-dliq))  
       qese = epsilo*es/(plcl-es)  
       thtesg(k) = tenv*(1.0e5_realkind/plcl)**(0.2854_realkind*(1.0_realkind-0.28_realkind*qese))* &
            exp((3374.6525_realkind/tenv-2.5403_realkind)*qese*(1.0_realkind +  0.81_realkind*qese))
       !
       !  compute adjusted abe(abeg).
       !
       !     gjnewent
       abeg = 0.0_realkind  
       do nk = k,ltopm1  
          nk1 = nk + 1  
          theteu(nk1) = theteu(nk)  
          !
          call tpmix2dd(p0(nk1),theteu(nk1),tgu(nk1),qgu(nk1))  
          !
          tvqu(nk1) = tgu(nk1)*(1.0_realkind + 0.608_realkind*qgu(nk1)-qliq(nk1) &
               -qice(nk1))
          if(nk==k) then  
             dzz = z00(klcl)-zlcl  
             dilbe =((tvlcl + tvqu(nk1))/(tven + tvg(nk1))-1.0_realkind)**dzz
          else  
             dzz = dza(nk)  
             dilbe =((tvqu(nk) + tvqu(nk1))/(tvg(nk) + tvg(nk1))-1.0_realkind)*dzz
          endif
          if(dilbe>0.0_realkind) abeg = abeg + dilbe*g  
          !
          !  dilute by entrainment by the rate as original updraft
          !
          call envirtht(p0(nk1),tg(nk1),qg(nk1),thteeg(nk1),&
               0.0_realkind,rl,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          !
          theteu(nk1) = theteu(nk1)*ddilfrc(nk1) + thteeg(nk1)*(1.0_realkind-ddilfrc(nk1))
       enddo
       !     gjnewent
       !
       !
       !  assume at least 90% of cape(abe) is removed by convection during
       !  the period timec
       !
       if(noitr==1) then  
          !     write(98,1060)fabe
          goto 265  
       endif
       dabe = max(abe-abeg,0.1_realkind*abe)  
       fabe = abeg/(abe+1.0e-8_realkind)  
       if(fabe>1.0_realkind.and.ishall==0) then  
          !     write(98,*)'updraft/downdraft couplet increases cape at this
          !    *grid point; no convection allowed!'
          goto 325  
       endif
       if(ncount/=1) then  
          dfda =(fabe-fabeold)/(ainc-aincold)  
          if(dfda>0.0_realkind) then  
             noitr = 1  
             ainc = aincold  
             goto 255  
          endif
       endif
       aincold = ainc  
       fabeold = fabe  
       if(ainc/aincmx>0.999_realkind.and.fabe>1.05_realkind-stab) then  
          !     write(98,*)' '
          !     write(98,*)'tau,i,j,=',ntsd,i,j
          !     write(98,1055)fabe
          goto 265  
       endif
       if(fabe<=1.05_realkind-stab.and.fabe>=0.95_realkind-stab) goto 265  
       if(ncount>10) then  
          !     write(98,*)' '
          !     write(98,*)'tau,i,j,=',ntsd,i,j
          !     write(98,1060)fabe
          goto 265  
       endif
       !
       !  if more than 10% of the original cape remains,increase the convec
       !  mass flux by the factor ainc:
       !
       if(abs(fabe)<1.e-14_realkind) then  
          ainc = ainc*0.5_realkind  
       else  
          ainc = ainc*stab*abe/(dabe+1.0e-8_realkind)  
       endif
255    ainc = min(aincmx,ainc)  
       !  if ainc becomes very small,effects of convection ! jsk mods
       !  will be minimal so just ignore it              ! jsk mods
       !     if(ainc<0.05)goto 325                        ! jsk mods
       if(ainc<0.05_realkind) then  
          ! jsk mods
          goto 325  
       endif
       !     ainc=max(ainc,0.05)                        ! jsk mods
       tder = tder2*ainc  
       pptflx = pptfl2*ainc  
       !     if(xtime<10.)then
       !     write(98,1080)lfs,ldb,ldt,timec,tadvec,nstep,ncount,
       !    *             fabeold,aincold
       !     endif
       !     gj
       !     gj   it is here that umf etc get converted into units of kg/s and
       !     gj   correctly scaled by ainc to satisfy the closure criterion tha
       !     gj   90% of cape is used up
       !     gj   this leads to correct raincv units later of cm/timestep
       !     gj
       !     gj   up to this point udr etc and udr2 etc have had units of
       !     gj   kg2/m4s2(derived from rei units). scaling is donr at this po
       !     gj   using ainc=(ems/(uer-der)*timec) with units of 1/kgm4s1 the r
       !     gj   uer etc now all being in units of kg/s
       !     gj
       do  nk = 1,ltop  
          umf(nk) = umf2(nk)*ainc  
          dmf(nk) = dmf2(nk)*ainc  
          detlq(nk) = detlq2(nk)*ainc  
          detic(nk) = detic2(nk)*ainc  
          udr(nk) = udr2(nk)*ainc  
          uer(nk) = uer2(nk)*ainc  
          der(nk) = der2(nk)*ainc  
          ddr(nk) = ddr2(nk)*ainc  
       enddo
       !
       !  go back up for another iteration
       !
       goto 175  
265    continue  
       !     gj
       do k = 1,kl-1  
          nk = kl-k + 1  
          if(k<=ltop) then  
             umf_out(i,nk) = umf(k)/dxsq  
             dmf_out(i,nk) = dmf(k)/dxsq  
             !     gj           umf_out(i,nk)=umf(k)/au0
             !     gj           dmf_out(i,nk)=dmf(k)/au0
             !     gjkfarm
             !     gjkfarm diagnostics
             if(umf_out(i,nk) >0.0_realkind) then  
                udr_out(i,nk) = udr(k)/(dxsq*rhoe(k))  
                uer_out(i,nk) = uer(k)/(dxsq*rhoe(k))  
                ddr_out(i,nk) = ddr(k)/(dxsq*rhoe(k))  
             else  
                udr_out(i,nk) = 0.0_realkind  
                ddr_out(i,nk) = 0.0_realkind  
                uer_out(i,nk) = 0.0_realkind  
             endif
             !     gjkfarm
             !k is higher than ltop
          else  
             umf_out(i,nk) = 0.0_realkind  
             dmf_out(i,nk) = 0.0_realkind  
             udr_out(i,nk) = 0.0_realkind  
             ddr_out(i,nk) = 0.0_realkind  
             uer_out(i,nk) = 0.0_realkind  
          endif
       enddo
       !     
       !     ccompute hydrometeor tendencies as is done for t,qv
       !
       !  frc2 is the fraction of total condensate      !  ppt fb mods
       !  generated that goes into precipitiation       !  ppt fb mods
       !     if(i==95 .and. j==82)then
       !     write(98,*)'pptflx,cpr,ainc =',pptflx,cpr,ainc
       !     endif
       if(cpr>0.0_realkind) then  
          !  ppt fb mods
          frc2 = pptflx/(cpr*ainc)  
       else  
          frc2 = 0.0_realkind  
       endif
       do nk = 1,ltop  
          qlpa(nk) = ql0(nk)  
          qipa(nk) = qi0(nk)  
          qrpa(nk) = qr0(nk)  
          qspa(nk) = qs0(nk)  
          !  ppt fb mods
          rainfb(nk) = pptliq(nk)*ainc*fbfrc*frc2  
          !  ppt fb mods
          snowfb(nk) = pptice(nk)*ainc*fbfrc*frc2  
       enddo
       do  ntc = 1,nstep  
          !
          !  assign hydrometeors concentrations at the top and bottom of each l
          !  based on the sign of omega
          !
          do  nk = 1,ltop  
             qlfxin(nk) = 0.0_realkind  
             qlfxout(nk) = 0.0_realkind  
             qifxin(nk) = 0.0_realkind  
             qifxout(nk) = 0.0_realkind  
             qrfxin(nk) = 0.0_realkind  
             qrfxout(nk) = 0.0_realkind  
             qsfxin(nk) = 0.0_realkind  
             qsfxout(nk) = 0.0_realkind  
          enddo
          do nk = 2,ltop  
             if(omg(nk) <=0.0_realkind) then  
                qlfxin(nk) =-fxm(nk)*qlpa(nk-1)  
                qifxin(nk) =-fxm(nk)*qipa(nk-1)  
                qrfxin(nk) =-fxm(nk)*qrpa(nk-1)  
                qsfxin(nk) =-fxm(nk)*qspa(nk-1)  
                qlfxout(nk-1) = qlfxout(nk-1) + qlfxin(nk)  
                qifxout(nk-1) = qifxout(nk-1) + qifxin(nk)  
                qrfxout(nk-1) = qrfxout(nk-1) + qrfxin(nk)  
                qsfxout(nk-1) = qsfxout(nk-1) + qsfxin(nk)  
             else  
                qlfxout(nk) = fxm(nk)*qlpa(nk)  
                qifxout(nk) = fxm(nk)*qipa(nk)  
                qrfxout(nk) = fxm(nk)*qrpa(nk)  
                qsfxout(nk) = fxm(nk)*qspa(nk)  
                qlfxin(nk-1) = qlfxin(nk-1) + qlfxout(nk)  
                qifxin(nk-1) = qifxin(nk-1) + qifxout(nk)  
                qrfxin(nk-1) = qrfxin(nk-1) + qrfxout(nk)  
                qsfxin(nk-1) = qsfxin(nk-1) + qsfxout(nk)  
             endif
          enddo
          !
          !  update the hydrometeor concentration values at each level
          !
          do  nk = 1,ltop  
             qlpa(nk) = qlpa(nk) +(qlfxin(nk) + detlq(nk) -qlfxout(nk))*dtime*emsd(nk)
             qipa(nk) = qipa(nk) +(qifxin(nk) + detic(nk) -qifxout(nk))*dtime*emsd(nk)
             !  reverse the mod that allows feedback of rain/snow that originates
             !  detrained air12/9/98jsk
             !     qrpa(nk)=qrpa(nk)+(qrfxin(nk)+qlqout(nk)*udr(nk)-qrfxout(nk)
             qrpa(nk) = qrpa(nk) +(qrfxin(nk)-qrfxout(nk) + rainfb(nk))*dtime*emsd(nk)
             !  ppt fb mods
             !     &              )*dtime*emsd(nk)                   !  ppt fb mods
             !  reverse the mod that allows feedback of rain/snow that originates
             !  detrained air12/9/98jsk
             !     qspa(nk)=qspa(nk)+(qsfxin(nk)+qicout(nk)*udr(nk)-qsfxout(nk)
             qspa(nk) = qspa(nk) +(qsfxin(nk)-qsfxout(nk) + snowfb(nk))*dtime*emsd(nk)
             !  ppt fb mods
             !     &              )*dtime*emsd(nk)                   !  ppt fb mods
             !
          enddo
       enddo
       do  nk = 1,ltop  
          qlg(nk) = qlpa(nk)  
          qig(nk) = qipa(nk)  
          qrg(nk) = qrpa(nk)  
          qsg(nk) = qspa(nk)  
       enddo
       if(iprnt) then  
          if(istop==1) then  
             if(istop==1) then  
                call stop_program('kain-fritsch')
             endif
          endif
       endif

       cndtnf =(1.0_realkind-eqfrc(lfs))*(qliq(lfs) + qice(lfs))*dmf(lfs)
       !     evaluate moisture budget
       qinit = 0.0_realkind  
       qfnl = 0.0_realkind  
       dpt = 0.0_realkind  
       do  nk = 1,ltop  
          dpt = dpt + dp(nk)  
          qinit = qinit + q0(nk)*ems(nk)  
          qfnl = qfnl + qg(nk)*ems(nk)  
          qfnl = qfnl +(qlg(nk) + qig(nk) + qrg(nk) + qsg(nk))*ems(nk)
       enddo
       !  ppt fb mods
       qfnl = qfnl + pptflx*timec*(1.0_realkind-fbfrc)  
       err2 =(qfnl-qinit)*100.0_realkind/qinit  
       if(abs(err2) >0.05_realkind.and.istop==0) then  
          print *, '!!!!!!!! moisture budget error in kfpara !!!'  
          !          write(*,1110) qinit,qfnl,err2  
          iprnt = .true.  
          istop = 1  
       endif



       relerr = err2*qinit/(pptflx*timec + 1.0e-10_realkind)  
       !  feedback to resolvable scale tendencies.
       !  if the advective time period(tadvec) is less than specified minim
       !  timec,allow feedback to occur only during tadvec
       if(tadvec<timec) nic = nint(tadvec/(0.5_realkind*dt2))  
       nca(i) = nic  
       if(ishall==1) then  
          timec = 2400.0_realkind  
          nca(i) = ncldck  
       endif
       do  k = 1,nlev  
          nk = nlev-k + 1  
          !  for eta model,evaporate or sublimate all hydrometeors so that the
          !  feedbacks are temperature and water vapor tendenciesthis may
          !  create supersaturated values of tg,but these will be removed by
          !  normal supersaturation-removal mechanisms
          !
          !     if(imoist(inest)/=2)then
          !
          !  if hydrometeors are not allowed,they must be evaporated or sublim
          !  and fed back as vapor,along with associated changes in temperatur
          !  note:  this will introduce changes in the convective temperature a
          !  water vapor feedback tendencies and may lead to supersaturated val
          !  of qg
          !
          !
          !  activate the following 2 statements to disallow hydrometeor
          !  feedbacks for eta model
          !
          !     dqcdt(i,k)=0.0_realkind
          !     dqrdt(i,k)=0.0_realkind
          !     dqldt(i,nk)=0.0_realkind
          !     dqidt(i,nk)=0.0_realkind
          !     dqrdt(i,nk)=0.0_realkind
          !     dqsdt(i,nk)=0.0_realkind
          !     else
          !     if(iexice/=1 .and. iice/=1) then
          !     if(imphys(inest)==3)then
          !
          !  if ice phase is not allowed,melt all frozen hydrometeors
          !
          !     cpm=cp*(1.+0.887*qg(k))
          !     tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
          !     dqldt(i,nk)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
          !     dqidt(i,nk)=0.
          !     dqrdt(i,nk)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
          !     dqsdt(i,nk)=0.0_realkind
          !     elseif(iexice==1 .and. iice==0)then
          !     elseif(imphys(inest)==4)then
          !
          !  if ice phase is allowed,but mixed phase is not,melt frozen hydro
          !  below the melting level,freeze liquid water above the melting lev
          !
          !     cpm=cp*(1.+0.887*qg(k))
          !     if(k<=ml)then
          !     tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
          !     elseif(k>ml)then
          !     tg(k)=tg(k)+(qlg(k)+qrg(k))*rlf/cpm
          !     endif
          !     dqldt(i,nk)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
          !     dqidt(i,nk)=0.0_realkind
          !     dqrdt(i,nk)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
          !     dqsdt(i,nk)=0.0_realkind
          !     elseif(iice==1 .and. iexice==0)then
          !     elseif(imphys(inest)>=5)then
          !
          !  if mixed phase hydrometeors are allowed,feed back convective tend
          !  of hydrometeors directly
          !
          !     dqldt(i,nk)=(qlg(k)-ql0(k))/timec
          !     dqidt(i,nk)=(qig(k)-qi0(k))/timec
          !     dqrdt(i,nk)=(qrg(k)-qr0(k))/timec
          !     dqsdt(i,nk)=(qsg(k)-qs0(k))/timec
          !     else
          !     print*,'this combination of imoist,iexice,iice not allowed!'
          !     stop 'kain-fritsch'
          !     endif
          !     endif
          !     dtdt(i,nk)=(tg(k)-t0(k))/timec
          !     dqdt(i,nk)=(qg(k)-q0(k))/timec
          !     qlg(k) = qlg(k)+qig(k)+qrg(k)+qsg(k)
          !     dqcdt(i,nk)=(qlg(k)-ql0(k))/timec
          !     if(ishall==1)then
          !     dtdt(i,nk)=0.0_realkind
          !     dqdt(i,nk)=0.0_realkind
          !     dqcdt(i,nk)=0.0_realkind
          !     else
          !     gj
          !     gj   here it is assumed all liquid water produced by
          !     gj   shallow conv is evaporated..it would also be
          !     gj   possiblr to do this in acondens using oldcov for
          !     gj   determining how much cw to evap.
          !     gj
          qlg(k) = qlg(k) + qig(k) + qrg(k) + qsg(k)  
          !     gj
          !     gj   evap shallow conv detraining cld water in acondens as
          !     gj   a function of oldcovshould allow for shallow conv
          !     gj   cldsthese should participate in radiation calcs
          !     gj   but not really in condensation calcs as that has been
          !     gj   done here already
          !     gj
          !     gj          if(ishall==1)then
          !     gj             rl=xlv0-xlv1*tg(k)
          !     gj             qg(k) = qg(k)+qlg(k)
          !     gj             tg(k) = tg(k)-rl*qlg(k)/cp
          !     gj             qlg(k) = 0.0_realkind
          !     gj          endif
          !     gj
          dtdt(i,nk) =(tg(k)-t0(k))/timec  
          dqdt(i,nk) =(qg(k)-q0(k))/timec  

          dqcdt(i,nk) =(qlg(k)-ql0(k))/timec  
          !below conv source layer
          if(k<lc) then  
             qliq_cj(i,nk) = 0.0_realkind  
             qice_cj(i,nk) = 0.0_realkind  
             qvap_cj(i,nk) = q0(k)  
          elseif(k>=lc.and.k<cldbase) then  
             qliq_cj(i,nk) = 0.0_realkind  
             qice_cj(i,nk) = 0.0_realkind  
             qvap_cj(i,nk) = qmix  
          else  
             qliq_cj(i,nk) = qliq(k) +(rainfb(k)*fbfrc*dtime*emsd(k))
             qice_cj(i,nk) = qice(k) +(snowfb(k)*fbfrc*dtime*emsd(k))
             qvap_cj(i,nk) = qdt(k)  
          endif
          if(ishall==1) then  
             psrc(i) = pmix0/100.0_realkind  
             pclb(i) = plcl0/100.0_realkind  
             emst = dpthmx*dxsq/g  
             umfb(i) = 100.0_realkind*vmflcl*timec*ainc/emst  
          endif
       enddo
       fbfrc_out(i) = fbfrc_shal  
       !  ppt fb mods
       raincv(i) = 0.1_realkind*0.5_realkind*dt2*pptflx*(1.0_realkind-fbfrc)/dxsq  
       rainrate(i) =(raincv(i)*10.0_realkind)/(0.5_realkind*dt2)  
       rnc = raincv(i)*real(nic,realkind)  
       nccnt = nccnt + 1  
325    continue 
    enddo
    return  
  end subroutine kfcumulus




  subroutine tpmix2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic,&
       ratio2,xlv1,xlv0)
    use kflut
    implicit none  
    real(kind=realkind) :: p,thes,tu,qu,qliq,qice,qnewlq,qnewic,ratio2,xlv1,xlv0
    real(kind=realkind) :: tp,qq,bth,pp,t00,t10,t01,t11,q00,q10,q01,q11,&
         tth,temp,qs,qnew,dq,qtot,rll,cp
    integer :: iptb,ithtb  
    !lookup table variables*

    !     scaling pressure & tt table index
    tp =(p-ptop)*rdpr  
    qq = tp-aint(tp)  


    iptb = int(tp) + 1  
    !              base and scaling factor for the
    !  scaling the & tt table index
    bth =(the0k(iptb + 1)-the0k(iptb))*qq + the0k(iptb)  
    tth =(thes-bth)*rdthk  
    pp = tth-aint(tth)  
    ithtb = int(tth) + 1  
    
    t00 = ttab(ithtb,iptb)  
    t10 = ttab(ithtb + 1,iptb)  
    t01 = ttab(ithtb,iptb + 1)  
    t11 = ttab(ithtb + 1,iptb + 1)  
    
    q00 = qstab(ithtb,iptb)  
    q10 = qstab(ithtb + 1,iptb)  
    q01 = qstab(ithtb,iptb + 1)  
    q11 = qstab(ithtb + 1,iptb + 1)  
    
    
    !              parcel temperature
    
    
    temp =(t00+(t10-t00)*pp+(t01-t00)*qq +(t00-t10-t01 + t11)*pp*qq)
    
    qs =(q00+(q10-q00)*pp+(q01-q00)*qq +(q00-q10-q01 + q11)*pp*qq)
    
    if(qs<=qu) then  
       qnew = qu-qs  
       qu = qs  
       goto 96  
    endif
    !
    !   if the parcel is subsaturated,temperature and mixing ratio must be
    !   adjustedif liquid water is present,it is allowed to evaporate
    !
    qnew = 0.0_realkind  
    dq = qs-qu  
    qtot = qliq + qice  
    !
    !   if there is enough liquid or ice to saturate the parcel,temp stays
    !   wet bulb value,vapor mixing ratio is at saturated level,and the mi
    !   ratios of liquid and ice are adjusted to make up the original satura
    !   deficit otherwise,any available liq or ice vaporizes and appropr
    !   adjustments to parcel temp; vapor,liquid,and ice mixing ratios are
    !
    !note that the liq and ice may be present in proportions slightly dif
    !   than suggested by the value of ratio2check to make sure that liq
    !   ice concentrations are not reduced to below zero when evaporation/
    !   sublimation occurs
    !
    !subsaturated values only occur in calculations involving various mix
    !updraft and environmental air for estimation of entrainment and detr
    !for these purposes,assume that reasonable estimates can be given us
    !liquid water saturation calculations only-i.e.,ignore the effect
    !ice phase in this process only
    !
    if(qtot>=dq) then  
       qliq = qliq-dq*qliq/qtot  
       qice = qice-dq*qice/qtot  
       qu = qs  
       goto 96  
    else  
       rll = xlv0-xlv1*temp  
       cp = 1005.7_realkind*(1.0_realkind + 0.89_realkind*qu)  
       if(qtot<1.0e-10_realkind) then  
          !if no liquid water or ice is available,temperature is given by:
          temp = temp + rll*(dq/(1.0_realkind + dq))/cp  
          goto 96  
       else  
          !if some liq water/ice is available,but not enough to achieve satura
          !   the temperature is given by:
          temp = temp + rll*((dq-qtot)/(1._realkind + dq-qtot))/cp
          qu = qu + qtot  
          qtot = 0.0_realkind  
       endif
       qliq = 0.0_realkind  
       qice = 0.0_realkind  
    endif
96  tu = temp  
    qnewlq = qnew  

    qnewic = 0.0_realkind  
    return  
  end subroutine tpmix2

  !  this subroutine integrates the area under the curve in the gaussian
  !  distributionthe numerical approximation to the integral is taken f
  !  "handbook of mathematical functions with formulas,graphs and mathema
  !  tables" ed. by abramowitz and stegun,nat'l bureau of standards appli
  !  mathematics series.  june,1964.,may,1968.
  !                                     jack kain
  !                                     7/6/89
  !    gaussian type mixing profile

  subroutine prof5(eq,ee,ud)  

    implicit none  

    real(kind=realkind) :: eq,ee,ud,sqrt2p,a1,a2,a3,p,sigma,fe,x,y,ey,&
         e45,t2,t1,c1,c2
    data sqrt2p,a1,a2,a3,p,sigma,fe/2.506628_realkind,0.4361836_realkind,&
         -0.1201676_realkind,0.9372980_realkind,0.33267_realkind,0.166666667_realkind,0.202765151_realkind/
    x =(eq-0.5_realkind)/sigma  
    y = 6.0_realkind*eq-3.0_realkind  
    ey = exp(y*y/(-2.0_realkind))  
    e45 = exp(-4.5_realkind)  
    t2 = 1.0_realkind/(1.0_realkind + p*abs(y))  
    t1 = 0.500498_realkind  
    c1 = a1*t1 + a2*t1*t1 + a3*t1*t1*t1  
    c2 = a1*t2 + a2*t2*t2 + a3*t2*t2*t2  
    if(y>=0.0_realkind) then  
       ee = sigma*(0.5_realkind*(sqrt2p-e45*c1-ey*c2) + sigma*&
            (e45-ey))-e45*eq*eq/2.0_realkind
       ud = sigma*(0.5_realkind*(ey*c2-e45*c1) + sigma*(e45-ey)) &
            -e45*(0.5_realkind + eq*eq/2.0_realkind-eq)
    else  
       ee = sigma*(0.5_realkind*(ey*c2-e45*c1) + sigma*(e45-ey)) -e45*eq*eq/2.0_realkind
       ud = sigma*(0.5_realkind*(sqrt2p-e45*c1-ey*c2) + sigma*&
            (e45-ey))-e45*(0.5_realkind + eq*eq/2.0_realkind-eq)
    endif
    ee = ee/fe  
    ud = ud/fe  
    return  
  end subroutine prof5



  subroutine tpmix2dd(p,thes,ts,qs)  
    use kflut
    implicit none  
    integer :: iptb,ithtb  



    real(kind=realkind) :: p,thes,ts,qs,tp,qq,bth,tth,pp,t00,t10,t01,t11,&
         q00,q10,q01,q11
    !lookup table variables*

    !     scaling pressure & tt table index
    tp =(p-ptop)*rdpr  
    qq = tp-aint(tp)  

    iptb = int(tp) + 1  
    !              base and scaling factor for the
    !  scaling the & tt table index
    bth =(the0k(iptb + 1)-the0k(iptb))*qq + the0k(iptb)  
    tth =(thes-bth)*rdthk  
    pp = tth-aint(tth)  
    ithtb = int(tth) + 1  
    !
    t00 = ttab(ithtb,iptb)  
    t10 = ttab(ithtb + 1,iptb)  
    t01 = ttab(ithtb,iptb + 1)  
    t11 = ttab(ithtb + 1,iptb + 1)  
    !
    q00 = qstab(ithtb,iptb)  
    q10 = qstab(ithtb + 1,iptb)  
    q01 = qstab(ithtb,iptb + 1)  
    q11 = qstab(ithtb + 1,iptb + 1)  

    !              parcel temperature
    ts =(t00 +(t10-t00)*pp +(t01-t00)*qq +(t00-t10-t01 + t11)*pp*qq)
    qs =(q00 +(q10-q00)*pp +(q01-q00)*qq +(q00-q10-q01 + q11)*pp*qq)
    return  
  end subroutine tpmix2dd


  subroutine envirtht(p1,t1,q1,tht1,r1,rl,aliq,bliq,cliq,&
       dliq,aice,bice,cice,dice)
    use confys,only:epsilo,tmelt  

    implicit none  

    real(kind=realkind) :: p1,t1,q1,tht1,r1,rl,aliq,bliq,cliq,dliq,aice,&
         bice,cice,dice,t00,p00,c1,c2,c3,c4,c5,ee,tlog,tdpt,&
         tsat,tht,tlogic,tfpt,tsatlq,tsatic

    data p00,c1,c2,c3,c4,c5/1.e5_realkind,3374.6525_realkind,2.5403_realkind,3114.834_realkind,&
         0.278296_realkind,1.0723e-3_realkind/
    t00 = tmelt  
    !
    !  calculate environmental equivalent potential temperature
    !
    if(r1<1.0e-6_realkind) then  
       ee = q1*p1/(epsilo + q1)  
       tlog = log(ee/aliq)  
       tdpt =(cliq-dliq*tlog)/(bliq-tlog)  
       tsat = tdpt-(0.212_realkind + 1.571e-3_realkind*(tdpt-t00)-4.36e-4_realkind*(t1-t00))*(t1-tdpt)
       tht = t1*(p00/p1)**(0.2854_realkind*(1.0_realkind-0.28_realkind*q1))  
       tht1 = tht*exp((c1/tsat-c2)*q1*(1.0_realkind + 0.81_realkind*q1))  
    elseif(abs(r1-1.0_realkind) <1.0e-6_realkind) then  
       ee = q1*p1/(epsilo + q1)  
       tlog = log(ee/aice)  
       tfpt =(cice-dice*tlog)/(bice-tlog)  
       tht = t1*(p00/p1)**(0.2854_realkind*(1.0_realkind-0.28_realkind*q1))  
       tsat = tfpt-(0.182_realkind + 1.13e-3_realkind*(tfpt-t00)-3.58e-4_realkind*(t1-t00))*(t1-tfpt)
       tht1 = tht*exp((c3/tsat-c4)*q1*(1.0_realkind + 0.81_realkind*q1))  
    else  
       ee = q1*p1/(epsilo + q1)  
       tlog = log(ee/aliq)  
       tdpt =(cliq-dliq*tlog)/(bliq-tlog)  
       tlogic = log(ee/aice)  
       tfpt =(cice-dice*tlogic)/(bice-tlogic)  
       tht = t1*(p00/p1)**(0.2854_realkind*(1.0_realkind-0.28_realkind*q1))  
       tsatlq = tdpt-(0.212_realkind + 1.571e-3_realkind*(tdpt-t00)-4.36e-4_realkind*(t1-t00))*(t1-tdpt)
       tsatic = tfpt-(0.182_realkind + 1.13e-3_realkind*(tfpt-t00)-3.58e-4_realkind*(t1-t00))*(t1-tfpt)
       tsat = r1*tsatic +(1.0_realkind-r1)*tsatlq  
       tht1 = tht*exp(rl*q1*c5/tsat*(1.0_realkind + 0.81_realkind*q1))  
    endif
    return  
  end subroutine envirtht



  subroutine condload(qliq,qice,wtw,dz,boterm,enterm,rate,&
       qnewlq,qnewic,qlqout,qicout)

    implicit none  

    real(kind=realkind) :: qliq,qice,wtw,dz,boterm,enterm,rate,qnewlq,qnewic,&
         qlqout,qicout,g,qnew,qtot,qest,g1,wavg,conv,ratio3,oldq,&
         ratio4,dq,pptdrg
    real(kind=realkind) :: eps
    !  9/18/88this precipitation fallout scheme is based on the scheme us
    !  by ogura and cho(1973).  liquid water fallout from a parcel is cal-
    !  culated using the equation dq=-rate*q*dt,but to simulate a quasi-
    !  continuous process,and to eliminate a dependency on vertical
    !  resolution this is expressed as q=q*exp(-rate*dz).
    data g/9.81_realkind/ 
    qtot = qliq + qice  
    qnew = qnewlq + qnewic  
    !
    !  estimate the vertical velocity so that an average vertical velocity c
    !  be calculated to estimate the time required for ascent between model
    !  levels
    !
    qest = 0.5_realkind*(qtot + qnew)  
    g1 = wtw + boterm-enterm-2.0_realkind*g*dz*qest/1.5_realkind  
    if(g1<0.0_realkind) g1 = 0.0_realkind  
    wavg =(sqrt(wtw) + sqrt(g1))/2.0_realkind  
    eps=1.0e-6_realkind
    conv = rate*dz/(wavg+eps)  
    !gj
    !gj     ffrac determines the fraction of frozen condensate that
    !gj     should stay in the updraft and not precipitate out.
    !gj     beyond that already calculated. this term activates if the
    !gj     updraft is greater than 1m/s. at updraft velocity 3m/s
    !gj     or greater all frozen condensate stays in the updraft.
    !gj
    !gj        ffrac=1.-((wavg-1.0_realkind)/3.)
    !gj        ffrac=max(0.0_realkind,min(1.,ffrac))
    !gj
    !
    !  ratio3 is the fraction of liquid water in fresh condensate,ratio4 is
    !  the fraction of liquid water in the total amount of condensate involv
    !  in the precipitation process-note that only 75% of the fresh conden
    !  sate is is allowed to participate in the conversion process
    !gj mod to 75% newcond
    !
    ratio3 = qnewlq/(qnew + 1.0e-10_realkind)  
    !     oldq=qtot
    qtot = qtot + 0.60_realkind*qnew  
    oldq = qtot  
    ratio4 =(0.60_realkind*qnewlq + qliq)/(qtot + 1.0e-10_realkind)  
    qtot = qtot*exp(-conv)  
    !
    !  determine the amount of precipitation that falls out of the updraft
    !  parcel at this level
    !
    dq = oldq-qtot  
    qlqout = ratio4*dq  
    qicout =(1.0_realkind-ratio4)*dq  
    !gj
    !gj   decrease qicout by the assumed fractional amount of
    !gj   ice that should stay in the updraft. icesusp contains
    !gj   the extra ice that is now considered to stay in the updraft.
    !gj   qtot_tmp contains the total condensate(irrespective of phase)
    !gj   after precip removal. icesusp gets added directly back into
    !gj   qice array after pptdrg calculation.
    !gj   n.b. on output qicout & qlqout contain solid & liquid
    !gj   precipi fluxes assumed to have fallen out of this layer
    !gj   qice & qliq contains the phase seperated condensate remaining
    !gj   in and ascending with the updraft
    !gj
    !new      qicout2=qicout*ffrac
    !new      icesusp=qicout-qicout2
    !new      qtot_tmp=qtot+icesusp
    !new      qicout=qicout2
    !gj
    !
    !  estimate the mean load of condensate on the updraft in the layer,cal
    !  late vertical velocity
    !
    !new      pptdrg=0.5*(oldq+qtot_tmp-0.2*qnew)
    pptdrg = 0.5_realkind*(oldq + qtot-0.2_realkind*qnew)  
    wtw = wtw + boterm-enterm-2.0_realkind*g*dz*pptdrg/1.5_realkind  
    !
    !  determine the new liquid water and ice concentrations including losse
    !  due to precipitation and gains from condensation
    !
    qliq = ratio4*qtot + ratio3*0.4_realkind*qnew  
    qice =(1.0_realkind-ratio4)*qtot +(1.0_realkind-ratio3)*0.4_realkind*qnew  
    !new      qice=((1.-ratio4)*qtot)+((1.-ratio3)*0.4*qnew)+icesusp
    !gj
    qnewlq = 0.0_realkind  
    qnewic = 0.0_realkind  
    return  
  end subroutine condload


  subroutine dtfrznew(tu,p,thteu,qu,qfrz,qice,aliq,bliq,cliq,dliq)
    use confys,only:epsilo,tmelt  
    implicit none  

    real(kind=realkind) :: rlc,rls,rlf,cp,a,dtfrz,tu,es,qs,dqevap,qice,qu,&
         aliq,bliq,cliq,dliq,p,pi,thteu,qfrz
    !
    !allow the freezing of liquid water in the updraft to proceed as an
    !approximately linear function of temperature in the temperature rang
    !ttfrz to tbfrz
    !for colder termperatures,freeze all liquid water
    !thermodynamic properties are still calculated with respect to liquid
    !to allow the use of lookup table to extract tmp from thetae
    !
    rlc = 2.5e6_realkind-2369.276_realkind*(tu-tmelt)  
    rls = 2833922.0_realkind-259.532_realkind*(tu-tmelt)  
    rlf = rls-rlc  
    cp = 1005.7_realkind*(1.0_realkind + 0.89_realkind*qu)  
    !
    !  a = d(es)/dt is that calculated from buck's(1981) emperical formulas
    !  for saturation vapor pressure
    !
    a =(cliq-bliq*dliq)/((tu-dliq)*(tu-dliq))  
    dtfrz = rlf*qfrz/(cp + rls*qu*a)  

    tu = tu + dtfrz  
    es = aliq*exp((bliq*tu-cliq)/(tu-dliq))  
    qs = es*epsilo/(p-es)  
    !
    !freezing warms the air and it becomes unsaturatedassume that some
    !liquid water that is available for freezing evaporates to maintain s
    !tionsince this water has already been transferred to the ice cate
    !subtract it from ice concentration,then set updraft mixing ratio at
    !temperature to the saturation value
    !
    dqevap = qs-qu  
    qice = qice-dqevap  
    qu = qu + dqevap  
    pi =(1.0e5_realkind/p)**(0.2854_realkind*(1.0_realkind-0.28_realkind*qu))  
    thteu = tu*pi*exp((3374.6525_realkind/tu-2.5403_realkind)*qu*(1.0_realkind + 0.81_realkind*qu))
    return  
  end subroutine dtfrznew


end module akfcum_mod
