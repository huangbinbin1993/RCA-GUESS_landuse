module mod_grib
  use decomp,only:stop_program
  implicit none
  private
  integer,public,parameter::realgribkind=4

  integer,public,save::ngbmod
  integer,public,save::nfref
  integer,public,save::nrnd
  integer,public,save::ndbg
  integer,public,save::nvck
  integer,public,save::nuser
  real(kind=realgribkind),public,save::fref

  public gsbite,maxmin,engrib,degrib,grbmod
contains

subroutine grbmod(kgbmod,lanew)

    !     grbmod:  set grib software mode of operation
    !              but only on first call or if lanew is true
    !              (so lanew forces a new value to be written)

    implicit none

    integer:: kgbmod
    logical lofirs,lanew
    save lofirs
    data lofirs/.true./
    if(.not.lofirs.and..not.lanew)return
    lofirs=.false.
    ngbmod=kgbmod
  end subroutine grbmod

  subroutine grchk1(ksec1,kret)
    implicit none
    !     
    ! grchk1 - check parameters for section 1 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     check parameters for section 1 of grib code against
    !     valid values for grib code edition 1.
    !     
    !   interface.
    !     ----------
    !     
    !     call grchk1 (ksec1,kret)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec1      - array containing parameters for section

    !     1 of grib code.
    !     
    !     output parameters.
    !     ------------------
    !     
    !     kret       - return code.
    !     0   , no error encountered.
    !     1   , error in grib code parameter.
    !     
    !     method.
    !     -------
    !     
    !     values checked against current code/flag tables
    !     and against maximum or minimum permitted values.
    !     they are also checked against the current status
    !     of the implementation of gribex and ecmwf usage.
    !     
    !     externals.
    !     ----------
    !     
    !     none.
    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes for grib code.
    !     
    !     comments.
    !     ---------
    !     
    !     routine contains sections 0 to 2 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     j. hennessy      ecmwf      18.06.91
    !     
    !     modifications.
    !     --------------
    !     
    !     j. hennessy      ecmwf      30.08.91
    !     checks for bit-map present removed..
    !     check for negative value for century removed.
    !     
    !     j. hennessy      ecmwf      02.12.91
    !     changes to table 2 and parameter number checks.
    !     
    !     ----------------------------------------------------------------
    !     
    integer:: iac
    integer:: iret


    integer:: j208
    integer:: j210
    integer:: j220
    integer:: j230
    integer:: j232

    integer:: kret

    integer,parameter::jpavac=9
    integer,parameter::jp1=4
    integer,parameter::jp3=26
    integer,parameter::jp4=9
    integer,parameter::jp5=14

    integer:: ksec1(*)

    integer:: itab1(jp1)
    integer:: itab3(jp3)
    integer:: itab4(jp4)
    integer:: itab5(jp5)
    integer:: iavac(jpavac)

    logical lpsbufr

    !     valid values given in code table 1.
    data itab1  /0,128,64,192/
    !     valid values given in code table 3.
    data itab3 /1,2,3,4,5,6,7,8,100,101,102,103,104,105,106,107,108,&
         109,110,111,112,121,128,141,160,200/
    !     valid values given in code table 4.
    data itab4 /0,1,2,3,4,5,6,7,254/
    !     valid values given in code table 5.
    data itab5 /0,1,2,3,4,5,10,113,114,115,116,117,123,124/
    !     valid values given in code table 5, for averages and
    !     accumulations.
    data iavac /-1,-1,113,114,115,116,117,123,124/
    kret = 0
    !         section 2 . check values against code tables and extreme values.
    !         check parameter table version number.
    if(ksec1(1)<1.or.ksec1(1)>254)then
       write(*,9001) ksec1(1)
       kret = 1
    endif

    lpsbufr=ksec1(1)==137

    !         check identification of centre. code table 0.
    !     currently only values 1 to 98 inclusive are used.
    !     missing value indicator(255) is considered an error
    !     when coding data.

    if(ksec1(2)<1.or.ksec1(2)>254)then
       write(*,9002) ksec1(2)
       kret = 1
    endif

    !         check generating process identification number.

    if(ksec1(3)<1.or.ksec1(3)>255) then
       write(*,9003) ksec1(3)
       kret = 1
    endif
    !         check grid definition.

    if(ksec1(4)<1.or.ksec1(4)>255)then
       write(*,9004) ksec1(4)
       kret = 1
    endif

    !         check flag. code table 1.

    do 208 j208=1,jp1
       if(ksec1(5)==itab1(j208)) go to 209
208 enddo

    write(*,9005) ksec1(5)
    kret = 1

209 continue
    !         cross check that, if uncatalogued grid is specified, section
    !     2, grid description section, is included.

    if(ksec1(4)==255.and.(ksec1(5)==0.or.ksec1(5)==64))then
       kret = 1
       write(*,9014)
    endif

    !         check parameter indicator. code table 2.
    !     missing value indicator(255) is considered an error
    !     when coding data.

    if(ksec1(6)<1.or.ksec1(6)>254)then
       if(lpsbufr.and.(ksec1(6)>=0.and.ksec1(6)<=255))then
       else
          write(*,9006) ksec1(6)
          kret = 1

       endif
    endif
    !         check ecmwf parameter table number and parameter indicators.


    !         check international table useage.

    if(ksec1(1)<127.and.ksec1(6)>127)then
#ifdef debug
       write(*,9024) ksec1(6) , ksec1(1)
#endif
    endif
    !         check indicator of type of level. code table 3.

    iret = 0


    if(lpsbufr.and.(ksec1(7)>=0.and.ksec1(7)<=255))goto 211


    do 210 j210=1,jp3
       if(ksec1(7)==itab3(j210)) go to 211
210 enddo

    iret = 1

211 continue
    !     ecmwf uses 200 for pseudo-levels.
    !     ecmwf uses 120 for experimental space view perspective.

    if((ksec1(7)==200.or.ksec1(7)==120).and.ksec1(2)==98)then
       iret = 0
    endif


    !     allow level type 0 to mean 'meaningless'

    if(ksec1(7)==0)then
       iret=0
       write(*,9015) ksec1(7)
    endif

    if(iret==1)then
       write(*,9015) ksec1(7)
       kret = 1
    endif
    !     for certain level types no description is necessary and
    !     those fields should be 0.

    !         check year of century.

    if(ksec1(10)<1.or.ksec1(10)>100)then
       write(*,9008) ksec1(10)
       kret = 1
    endif

    !         month check.
    !     the following avoids failures for yearly climatology yymmddhh=00000000
    if(ksec1(11)<0.or.ksec1(11)>12)then

       write(*,9009) ksec1(11)
       kret = 1
    endif

    !         day check.

    !     the following avoids failures for climatology yymmddhh=00mm0000
    if(ksec1(12)<0.or.ksec1(12)>31)then

       write(*,9010) ksec1(12)
       kret = 1
    endif

    !         hour check.

    if(ksec1(13)<0.or.ksec1(13)>23)then
       write(*,9011) ksec1(13)
       kret = 1
    endif
    !         minute check.

    if(ksec1(14)<0.or.ksec1(14)>59)then
       write(*,9012) ksec1(14)
       kret = 1
    endif

    !         indicator of unit of time check. code table 4.

    do 220 j220=1,jp4
       if(ksec1(15)==itab4(j220)) go to 221
220 enddo

    write(*,9013) ksec1(15)
    kret = 1

221 continue
    !         time range indicator check. code table 5.

    do 230 j230=1,jp5
       if(ksec1(18)==itab5(j230)) go to 231
230 enddo

    write(*,9019) ksec1(18)
    kret = 1

231 continue
    !         cross check time range indicator and number averaged or
    !     accumulated.

    iac = 0
    do 232 j232=1,jpavac
       if(ksec1(18)==iavac(j232)) iac = 1
232 enddo

    !     if average or accumulation check for valid numbers
    !     of included and missing values.

    if(iac==1)then
       if(ksec1(19)<2)then
          write(*,9016) ksec1(18) , ksec1(19)
          kret = 1
       endif
       if(ksec1(20)<0.or.ksec1(20)>=ksec1(19))then
          write(*,9020) ksec1(18) , ksec1(20)
          kret = 1
       endif
    endif
    !         century check.

    !     ecmwf data start in 20th century(1900 for some climate fields).

    if(ksec1(21)<19.and.ksec1(2)==98)then
       write(*,9021) ksec1(21)
       kret = 1
    endif

    !         decimal scale factor check.

    !     at ecmwf the scale factor is always 0.

    if(ksec1(23)/=0.and.ksec1(2)==98)then
       write(*,9022) ksec1(23)
       kret = 1
    endif
    !         section 9 . return to calling routine. format statements.
900 continue

9001 format(1h ,'grchk1 : invalid parameter table version number - ',&
         i5)
9002 format(1h ,'grchk1 : invalid identification of centre - ',i5)
9003 format(1h ,'grchk1 : invalid generating process - ',i5)
9004 format(1h ,'grchk1 : invalid grid definition - ',i5)
9005 format(1h ,'grchk1 : invalid flag field - ',i5,' decimal.')
9006 format(1h ,'grchk1 : invalid indicator of parameter - ',i5)
9007 format(1h ,'grchk1 : inconsistent version number ecmwf code ',&
         'table 2 - ',i5,' and indicator of parameter - ',i5)
9008 format(1h ,'grchk1 : invalid year of century - ',i5)
9009 format(1h ,'grchk1 : invalid month - ',i5)
9010 format(1h ,'grchk1 : invalid day - ',i5)
9011 format(1h ,'grchk1 : invalid hour - ',i5)
9012 format(1h ,'grchk1 : invalid minute - ',i5)
9013 format(1h ,'grchk1 : invalid indicator of unit of time - ',i5)
9014 format(1h ,'grchk1 : uncatalogued grid and no section 2.')

9015 format(1h ,'grchk1 : invalid indicator of type of level - ',i5)

9016 format(1h ,'grchk1 : inconsistent time range indicator',&
         ' - ',i5,' and number included in averages - ',i5)
9019 format(1h ,'grchk1 : invalid time range indicator - ',i5)
9020 format(1h ,'grchk1 : inconsistent time range indicator', &
         ' - ',i5,' and number missing from averages - ',i5)
9021 format(1h ,'grchk1 : invalid century of reference time - ',i5)
9022 format(1h ,'grchk1 : invalid decimal scale factor - ',i5)
9023 format(1h ,'grchk1 : for level type ',i3,' descriptions are',&
         'invalid - ',i5,3x,i5)
9024 format(1h ,'grchk1 : ** warning ** parameter number ',i3,&
         ' is not defined in international table number ',i3,'.')

    return

  end subroutine grchk1

  subroutine confp3(pval,kexp,kmant,kbits,kround)
    implicit none
    !     
    ! confp3 - convert floating point number to grib representation.
    !     
    !     purpose.
    !     --------
    !     
    !     convert floating point number from machine
    !     representation to grib representation.
    !     
    !   interface.
    !     ----------
    !     
    !     call confp3(pval,kexp,kmant,kbits,kround)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     pval    - floating point number to be converted.
    !     kbits   - number of bits in computer word.
    !     kround  - conversion type.

    !     0 , closest number in grib format less than
    !     original number.
    !     1 , closest number in grib format to the
    !     original number(equal to, greater than or
    !     less than original number).
    !     10 , as for 0 but with debug printout.
    !     11 , as for 1 but with debug printout.
    !     
    !     output parameters.
    !     -----------------
    !     
    !     kexp    - 8 bit signed exponent.
    !     
    !     kmant   - 24 bit mantissa.
    !     
    !     method.
    !     -------
    !     
    !     floating point number represented as 8 bit signed
    !     exponent and 24 bit mantissa in integer values.
    !     
    !     externals.
    !     ----------
    !     
    !     decfp2

    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes re grib representation.
    !     
    !     comments.
    !     ---------
    !     
    !     routine aborts if an invalid conversion type
    !     parameter is used or if a 24 bit mantissa is not
    !     produced.
    !     
    !     routine contains sections 0 to 2 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     john hennessy   ecmwf   18.06.91
    !     
    !     modifications.
    !     --------------
    !     
    !     john hennessy   ecmwf   24.09.91
    !     corrections made to descriptions of input parameter kround.
    !     changes to comments and format statements.
    !     
    !     -----------------------------------------------------------------
    !     
    integer:: iexp
    integer:: ipr
    integer:: iround
    integer:: isign

    integer:: kbits
    integer:: kexp
    integer:: kmant
    integer:: kround

    real(kind=realgribkind) pval

    real(kind=realgribkind) zeps
    real(kind=realgribkind) zref
    real(kind=realgribkind) zval

    !     debug print switch.

    if(kround>=10)then
       ipr    = 1
       iround = kround - 10
    else
       ipr    = 0
       iround = kround
    endif
    !     check conversion type parameter.

    if(iround/=0.and.iround/=1)then
       write(*,9004) kround
       call stop_program( 'confp3 stop 1')
    endif
    !          convert value of 0.0.


    if(abs(pval)<1.e-14_realgribkind)then
       kexp  = 0
       kmant = 0
       iexp  = 0
       isign = 0
       go to 900
    endif

    zeps = 1.0e-12_realgribkind
    if(kbits==32) zeps = 1.0e-8_realgribkind
    zref = pval

    !     sign of value.

    isign = 0
    if(zref<0._realgribkind)then
       isign = 128
       zref  = - zref
    endif
    !     exponent.

    iexp = int(log(zref)*(1.0_realgribkind/log(16.0_realgribkind))+64.0_realgribkind+1.0_realgribkind+zeps)

    if(iexp<0  ) iexp = 0
    if(kbits==32) iexp = max(39,iexp)
    if(iexp>127) iexp = 127
    if(kbits==32) iexp = min(101,iexp)

    !     mantissa.

    if(iround==0)then
       !     closest number in grib format less than original number.

       if(isign==0)then
          !     truncate for positive numbers.
          kmant = int(zref/16.0_realgribkind**real(iexp-70,realgribkind))
       else
          !     round up for negative numbers.
          kmant = nint(zref/16.0_realgribkind**real(iexp-70,realgribkind)+0.5_realgribkind)
       endif
    else
       !     closest number in grib format to the original number
       !    (equal to, greater than or less than original number).
       kmant = nint(zref/16.0_realgribkind**real(iexp-70,realgribkind))
    endif
    !     check that mantissa value does not exceed 24 bits.
    !     16777215 = 2**24 - 1

    if(kmant>16777215)then
       iexp = iexp + 1
       if(iround==0)then
          !     closest number in grib format less than original
          !     number.

          if(isign==0)then
             !     truncate for positive numbers.
             kmant = int(zref/16.0_realgribkind**real(iexp-70,realgribkind))
          else
             !     round up for negative numbers.
             kmant = nint(zref/16.0_realgribkind**real(iexp-70,realgribkind)+0.5_realgribkind)
          endif
       else
          !     closest number in grib format to the original number
          !    (equal to, greater or less than original number).
          kmant = nint(zref/16.0_realgribkind**real(iexp-70,realgribkind))
       endif

       if(kmant>16777215)then
          write(*,9001)
          write(*,9002) pval
          write(*,9003) isign , iexp , kmant

          call stop_program( 'confp3 stop 2')
       else
          write(*,9005)
       endif
    endif
    !     add sign bit to exponent.

    kexp = iexp + isign

900 continue

    if(ipr==1)then
       write(*,9006) kround
       write(*,9002) pval
       call decfp2(zval,kexp,kmant)
       write(*,9007) zval
       write(*,9003) isign , iexp , kmant
    endif

    return

9001 format(1h ,'confp3 : mantissa overflow fatal error.')

9002 format(1h ,'confp3 : original number = ',f30.20)

9003 format(1h ,'confp3 : sign = ',i3,', exponent = ',i3,', mantissa = ',i12)

9004 format(1h ,'confp3 : invalid conversion type parameter = ',i4)

9005 format(1h ,'confp3 : mantissa overflow recoverable error.')

9006 format(1h ,'confp3 : conversion type parameter = ',i4)

9007 format(1h ,'confp3 : converted to      ',f30.20)

  end subroutine confp3

  subroutine decfp2(pval,kexp,kmant)
    implicit none
    ! decfp2 - grib representation to floating point representation.
    !           convert grib representation of a floating point
    !           number to machine representation.
    !               input parameters.
    !               kexp    - 8 bit signed exponent.
    !               kmant   - 24 bit mantissa.
    !
    !               output parameters.
    !               pval    - floating point number represented
    !                         by kexp and kmant.
    !
    !     method.
    !           floating point number represented as 8 bit exponent
    !           and 24 bit mantissa in integer values converted to
    !           machine floating point format.
    !     reference.
    !           wmo manual on codes re grib representation.
    !
    !           rewritten from decfp, to conform to programming standards.
    !           sign bit on 0 value now ignored, if present.
    !
    !           routine contains sections 0 to 2 and section 9.
    !
    !     author.
    !           john hennessy   ecmwf   18.06.91
    integer::iexp
    integer::isign
    integer::kexp
    integer::kmant

    real(kind=realgribkind)::pval
    if((kexp==128.or.kexp==0).and.kmant==0)then
       pval = 0.0_realgribkind
    else
       !     go to 900
       !  endif
       !     sign of value.

       iexp  = kexp
       isign = 1

       if(iexp>=128)then
          iexp  = iexp - 128
          isign = -1
       end if
       !     decode value.
       pval = real(isign*kmant,realgribkind)*16._realgribkind**real(iexp-64-6,realgribkind)
    endif

    return

  end subroutine decfp2
  subroutine degrib(kgrib,kleng,field,fieldsize,kb1,kb2,pw2,kb4, &
       pmdi,kerr,returnSize)

    implicit none
    !     
    !     subroutine degrib is the new interface between the old routines
    !     and the new gribex and replaces the old routine degrib. therefore
    !     it has the same parameter list(some of them aren't of much use).
    !     
    !     input :

    !     kgrib : bit-stream
    !     kleng : length of bit-stream array
    !     
    !     output :
    !     field : field with grib - code
    !     fieldsize : length of field
    !     kerr  : error index
    !     kb1, kb2, kb3, kb4 : arrays with integer values of blocks 1 - 4
    !     klb1, klb2, klb3, klb4 : length of arrays kb1 - kb4
    !     pw2 : array with the real values of block 2(=vertical levels)
    !     
    !     rest :
    !     ksec0, ksec1, ksec2, ksec3, ksec4 : arrays with the
    !     integer values of block 0 - 4
    !     nsec0, nsec1, nsec2, nsec3, nsec4 : length of ksec0 - ksec4
    !     psec2, psec3 : arrays with real parts of block 2 and block 3
    !     npsec2, npsec3 : length of arrays psec2 and psec3
    !     hfunc : d to decode, i to decode only sections 1 and 2
    !     

    integer,intent(in)::fieldsize
    integer,intent(out)::returnSize
    real(kind=realgribkind)::  field(fieldsize) , pw2(*)
!    integer:: kb1(*) , kb2(*) ,  kb4(*),kgrib(128000) 
    integer:: kb1(*) , kb2(*) ,  kb4(*),kgrib(600000) 
    integer:: ksec0(2) , ksec1(25) , ksec2(50) , ksec3(2) , ksec4(42)
    real(kind=realgribkind):: psec2(400) , psec3(2)
    character(len=1)::hfunc
    integer:: kret,kleng,kword,kerr,ipt,j,ierr
    real(kind=realgribkind):: pmdi,znorm


    returnSize=fieldSize

    !     set variables

    kret = 0
    if( fieldsize == 0 ) then
       hfunc = 'I'
    else
       hfunc = 'D'
    endif

    !     decode
    ksec3(2)= int(pmdi)!imdi
    psec3(2)= pmdi !zmdi

    call gribex(ksec0,ksec1,ksec2,psec2,ksec3,psec3,ksec4, &
         field,fieldsize,kgrib,kleng,kword,hfunc,kret)
    if(kret<0)kret=0
    kerr = kret
    if(fieldsize>0)then
       !fieldsize = min(ksec2(2)*ksec2(3),fieldsize)
       returnsize = min(ksec2(2)*ksec2(3),fieldsize)
    endif
    !     in case, it's an older version, shift parameter into new wmo - code

    if(mod(ngbmod,2)==0)then
       if( ksec0(2) /= 1 .or. ksec1(1) /= 1 ) then
          !     normalise ukmo
          if(ksec1(1)==0.and.ksec1(2)==74)then
             if(fieldsize/=0)then
                !                  ipt=min(ksec2(2)*ksec2(3),fieldsize)
                ipt=min(ksec2(2)*ksec2(3),returnsize)
                if(ksec1(6)==1)then
                   !     pressure in hpa becometh pa
                   znorm=100._realgribkind
                else if( ksec1(6)==2)then
                   !     geopotential in dam becometh m2/s2
                   znorm=98.065_realgribkind
                else
                   znorm=1._realgribkind
                endif
                do j=1,ipt
                   field(j)=field(j)*znorm
                enddo
             endif
          endif
          call gb1tb2(ksec0,ksec1,ierr)
          if(ierr/=0)then
             write(6,*) 'warning: table 2 parameter numbers not converted'
          endif
       endif
    endif

    !     shift values

    call gbe2h1(ksec0,ksec1,kb1)
    call gbe2h2(ksec2,kb2,psec2,pw2)
    call gbe2h4(ksec2,ksec4,kb4)
    return
  end subroutine degrib
  subroutine engrib(kgrib,kleng,field,fieldSize,kb1,kb2,pw2,kb4, &
       pacc,kb,pmdi,kerr)
    implicit none
    !     
    !     subroutine engrib is the new interface between the old routines
    !     and the new gribex and replaces the old routine engrib. therefore
    !     it has the same parameter list(some of them aren't of much use).
    !     
    !     input :
    !     field : field with grib - code
    !     fieldSize : length of field
    !     kb1, kb2, kb3, kb4 : arrays with integer values of blocks 1 - 4
    !     klb1, klb2, klb3, klb4 : length of arrays kb1 - kb4
    !     pw2 : array with the real values of block 2(=vertical levels)
    !     pacc : accuracy
    !     kb : number of bits per packed value
    !     
    !     output :
    !     kgrib : bit-stream
    !     kleng : length of bit-stream array
    !     kerr  : error index
    !     
    !     rest :
    !     ksec0, ksec1, ksec2, ksec3, ksec4 : arrays with the
    !     integer values of block 0 - 4
    !     psec2, psec3 : arrays with real parts of block 2 and block 3
    !     hfunc : c to code
    !     
    !     overbodig :
    !     kbpw  : dummy value, only here because of the parameter list
    !     kmask : dummy value, only here because of the parameter list
    !     
    !     pmdi  : missing data indicator
    !     
    !     variables
    !     
    character(len=1):: hfunc
    integer,intent(in)::fieldSize
    integer::kret,kb,mdat,indat,jdat,kleng,kword,kerr
    real(kind=realgribkind)::pacc,zmax,zmin,pmdi
!    integer:: kb1(*) , kb2(*) , kb4(*),kgrib(128000) 
    integer:: kb1(*) , kb2(*) , kb4(*),kgrib(600000) 
    real(kind=realgribkind)::  field(fieldSize) , pw2(*)
    integer::ksec0(2), ksec1(25) , ksec2(50) , ksec3(2) , ksec4(42)
    real(kind=realgribkind):: psec2(400) , psec3(2)


    !     set variables

    kret = 0
    hfunc = 'C'
    !     
    !     set accuracy
    !     pacc = 0  ==>  kb = kb(mostly 12, for hirlam 15)
    !     pacc > 0  ==>  kb = log((max-min)/pacc)/log(2)
    !     pacc < 0  ==>  kb = -pacc
    !     
    if( abs(pacc) < 1.e-14_realgribkind ) then
       kb4(11) = kb
    else
       if( pacc < 0._realgribkind ) then
          kb4(11) = int(abs( pacc ))
       else
          call maxmin(field,fieldSize,zmax,zmin,pmdi)
          if(zmax-zmin<pacc)then
             kb4(11)=0
          else
             kb4(11)=int(log((zmax-zmin)/pacc)/log(2._realgribkind)+1._realgribkind)
          endif
       endif
    endif

    !     limit the number of bits per word to some practical value
    !     
    kb4(11)=min(kb4(11),24)
    mdat = fieldSize
    !     
    !     check if a bit-map is required
    indat = kb2(7) * kb2(9)
    kb1(8)=128
    do 210 jdat=1,indat
       if( abs(field(jdat)-pmdi)<1.e-14_realgribkind)then 
          kb1(8)=192
       endif
210 enddo

    !     set data of section 0, 1, 2 and 4 (section 3 = 0)
    !     
    !      call gbset0(ksec0)
    ksec0(1) = 9999
    ksec0(2) = 1
    call gbh2e1(kb1,ksec1)
    call gbh2e2(kb2,ksec2,pw2,psec2)
    call gbh2e4(kb4,ksec4,mdat)

    !     code
    ksec3(1)=0

    ksec3(2)= int(pmdi) !imdi

    psec3(2)=pmdi !zmdi

    !     
    call gribex(ksec0,ksec1,ksec2,psec2,ksec3,psec3,ksec4, &
         field,fieldSize,kgrib,kleng,kword,hfunc,kret)
    if(kret<0)kret=0
    kerr = kret
    kleng = kword
    !     
    return
  end subroutine engrib


  subroutine exscal(pdata,klen,pref,pscale)
    implicit none
    integer,intent(in):: klen
    real(kind=realgribkind),intent(in)::pref
    real(kind=realgribkind),intent(in)::pscale
    real(kind=realgribkind):: pdata(klen)

    if(realgribkind==8)then
       call exscalc(pdata,klen,pref,pscale)
    else
       call exscalcf(pdata,klen,pref,pscale)
    endif
    return
  end subroutine exscal
  subroutine gb1tb2(ksec0,ksec1,kerr)
    !     
    !     This subroutine transforms the old parameter number, used
    !     in HIRLAM and ECMWF files, to the new parameter number as defined
    !     by the WMO.
    !     
    !     Input :
    !     KSEC0 : array with GRIB  section 0
    !     KSEC1 : array with GRIB block number 1.
    !     
    !     Output :
    !     KSEC1 : parameter, level and height of level changed to WMO-code
    !     KERR : error index

    !     0 = no errors
    !     1 = wrong parameter
    !     2 = no WMO equivalent
    !     3 = wrong number for the centre-identifier
    !     4 = trying to convert GRIB version 2 or higher - not coded yet
    !     
    !     HIRLAM      -->        WMO
    !     ------------------------------------------------
    !     Pressure                     101        -->         1
    !     Geopotential                 102        -->         6
    !     Temperature                  104        -->        11
    !     Dew-point temperature        110        -->        17
    !     Specific humidity            112        -->        51
    !     Relative humidity            113        -->        52
    !     U-component of wind          123        -->        33
    !     V-component of wind          124        -->        34
    !     Relative vorticity           130        -->        43
    !     Horizontal divergence        134        -->
    !     Vertical velocity            140        -->        39
    !     Soil wetness                 147        -->        86
    !     Large-scale precipitation    150        -->        62
    !     Total precipitation          150 lev 14 -->        61
    !     Snow depth                   151        -->        66 + level corr.
    !     Convective precipitation     153        -->        63
    !     Sea-surface temperature      161        -->        11 + level corr.
    !     Cloud cover                  179        -->        71
    !     Surface roughness            180        -->        83 + level corr.
    !     Fraction of land             181        -->        81
    !     Urbanization                 182        -->
    !     Fraction of ice              183        -->        91 + level corr.
    !     Albedo                       184        -->        84
    !     Fraction of vegetation       185        -->        87
    !     X-component of stress        187        -->
    !     Y-component of stress        188        -->
    !     Mixed layer depth            217        -->        67
    !     
    !     
    !     ECMWF  -->    WMO
    !     --------------------------------------------------
    !     atmosph. tide                127   -->
    !     budget values                128   -->
    !     geopotential                 129   -->     6,7 + level corr.
    !     temperature                  130   -->    11
    !     u-component wind             131   -->    33
    !     v-component wind             132   -->    34
    !     specific humidity            133   -->    51
    !     surface pressure             134   -->     1 + level corr.
    !     vertical velocity            135   -->    39
    !     vorticity                    138   -->
    !     surface temperature          139   -->    11 + level corr.
    !     surface soil wetness         140   -->
    !     snow depth                   141   -->    66
    !     large scale precipit.        142   -->    62
    !     convective precipit.         143   -->    63
    !     snowfall                     144   -->    64
    !     boundary layer dissipat.     145   -->   123
    !     surface sensible heat flux   146   -->   122 + level corr.
    !     surface latent heat flux     147   -->   121 + level corr.
    !     surface stress               148   -->
    !     surface net radiation        149   -->
    !     top net radiation            150   -->
    !     mean sea level pressure      151   -->     1 + level corr.
    !     ln surface pressure          152   -->   152( better : 1 + corr. )
    !     divergence                   155   -->
    !     height(geopotential)        156   -->     7
    !     relative humidity            157   -->    52
    !     cloud cover(total)          164   -->    71
    !     10 m wind(u)                165   -->    33 + level corr.
    !     10 m wind(v)                166   -->    34 + level corr.
    !     2 m temperature              167   -->    11 + level corr.
    !     2 m dewpoint                 168   -->    17 + level corr.
    !     deep soil temperature        170   -->
    !     deep soil wetness            171   -->
    !     land sea mask                172   -->    81
    !     surface roughness            173   -->    83
    !     albedo                       174   -->    84
    !     surface solar radiation      176   -->
    !     surface thermal radiation    177   -->
    !     top solar radiation          178   -->
    !     top thermal radiation        179   -->
    !     E-W surface stress           180   -->
    !     N-S surface stress           181   -->
    !     evaporation                  182   -->    57
    !     clim. deep soil temperature  183   -->
    !     clim. deep soil wetness      184   -->
    !     convective cloud cover       185   -->    72
    !     low cloud cover              186   -->    73
    !     medium cloud cover           187   -->    74
    !     high cloud cover             188   -->    75
    !     E-W orographic variance      190   -->
    !     N-S orographic variance      191   -->
    !     NW-SE orographic variance    192   -->
    !     NE-SW orographic variance    193   -->
    !     latit. gravity wave stress   195   -->
    !     merid. gravity wave stress   196   -->
    !     gravity wave dissipation     197   -->
    !     skin reservoir content       198   -->
    !     percentage of vegetation     199   -->    87
    !     var. of sub-grid scale oro.  200   -->
    !     max. temperature at 2 m      201   -->    15 + level corr.
    !     min. temperature at 2 m      202   -->    16 + level corr.
    !     precipit. analysis weight    204   -->
    !     runoff                       205   -->
    !     total precipitation          228   -->    61
    !     cloud liquid water content   246   -->    76
    !     cloud ice water content      247   -->    58
    !     
    !     
    !     variables
    !     
    implicit none
    integer:: ipar,jpar,kerr,j
    parameter( ipar = 26 , jpar = 27 )
    integer::   ksec0(*)
    integer::   ksec1(*)
    integer:: ihirl(ipar) , iwmo(ipar)
    integer:: jecmw(jpar) , jwmo(jpar)
    data ihirl /101,102,104,110,112,113,123,124,130,134,140,147, &
         150,151,&
         153,161,179,180,181,182,183,184,185,187,188,217/
    data iwmo /   1,  6, 11, 17, 51, 52, 33, 34, 43,999, 39, 86, &
         62, 66, &
         63, 11, 71, 83, 81,182, 91, 84, 87,999,999,067/
    data jecmw /129,130,131,132,133,135,141,142,143,144,145,152,&
         156,157,164,172,173,174,182,185,186,187,188,199,228,&
         246,247/
    data jwmo /   6, 11, 33, 34, 51, 39, 66, 62, 63, 64,123,152, &
         7, 52, 71, 81, 83, 84, 57, 72, 73, 74, 75, 87, 61, &
         76, 58/

    !     set initial start values
    !     
    kerr = 0
    !     only for non standard codes
    if( ksec1(1) == 1 ) return
    !     only for grib version 0 or 1
    if(ksec0(2)>1)then
       kerr=4
       return
    endif

    !     transform values

    if( ksec1(2) == 96 .or. ksec1(2) == 86 ) then
       !     sea water temperature
       if( ksec1(6) == 161 ) then
          ksec1(6) = 104
          ksec1(7) = 102
          ksec1(8) = 0
          ksec1(9) = 0
       endif

       if(ksec1(6)==150.and.ksec1(8)==14)then
          !     total precipitation was sometimes indicated by level=14
          ksec1(6)=61
          ksec1(8)=0
          ksec1(1)=1
          goto 101
       endif
       do 100 j=1,ipar
          if( ksec1(6) == ihirl(j) ) then
             ksec1(6) = iwmo(j)
             if(ksec1(6)==91)then

                !     fraction of ice over sea
                ksec1(7)=102
                ksec1(8)=0
                ksec1(9)=0
             elseif(ksec1(6)==83.or.ksec1(6)==66)then
                !     roughness over sea or snow over ice?
                if(ksec1(7)==105.and.ksec1(8)==999)then
                   ksec1(7)=102
                   ksec1(8)=0
                   ksec1(9)=0
                endif
             endif
             ksec1(1)=1
             goto 101
          endif
100    enddo

       !     the following would result in an error for undefined transformations:
       !     ksec1(6) = 998
    else if( ksec1(2)==74)then
       !     if ukmo
       if(ksec1(6)==2)then
          ksec1(6)=6
       else if(ksec1(6)==123)then
          ksec1(6)=33
       else if(ksec1(6)==124)then
          ksec1(6)=34
       else if(ksec1(6)==150)then
          ksec1(6)=61
          ksec1(7)=105
          ksec1(8)=0
       endif
       ksec1(1)=1
       goto 101
    else if( ksec1(2) == 98 .and. ksec1(1)==128 ) then

       !     if ecmwf-parameter
       do 200 j=1,jpar
          if( ksec1(6) == jecmw(j) ) then
             ksec1(6) = jwmo(j)
             ksec1(1)=1
             if(ksec1(6)==6.and.ksec1(8)==1)then
                ksec1(7)=105
                ksec1(8)=0
             endif
             goto 101
          endif
200    enddo
       if(ksec1(6)==139)then
          ksec1(6)=11
          ksec1(7)=105
          ksec1(8)=0
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       elseif(ksec1(6)==151)then
          ksec1(6)=1
          ksec1(7)=103
          ksec1(8)=0
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       elseif(ksec1(6)==165)then
          ksec1(6)=33
          ksec1(7)=105
          ksec1(8)=10
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       elseif(ksec1(6)==166)then
          ksec1(6)=34
          ksec1(7)=105
          ksec1(8)=10
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       elseif(ksec1(6)==167)then
          ksec1(6)=11
          ksec1(7)=105
          ksec1(8)=2
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       elseif(ksec1(6)==168)then
          ksec1(6)=17
          ksec1(7)=105
          ksec1(8)=2
          ksec1(9)=0
          ksec1(1)=1
          goto 101
       endif

       !     the following would result in an error for undefined transformations:
       !     ksec1(6) = 998
    else
       kerr = 3
    endif
101 continue

    if( ksec1(6) == 998 ) then
       kerr = 1
    else if( ksec1(6) == 999 ) then
       kerr = 2
    endif

    return
  end subroutine gb1tb2

  subroutine gbctb2(kb1,yaname,kdir,kerr)
    !
    ! GBCTB2 - Convert indicator of parameter to string or vice versa
    !
    !     INPUT:
    !        KB1 :    Block 1 of GRIB message
    !        YANAME : Name of indicator of parameter
    !        KDIR :   Direction of translation:

    !                 0: Return name given indicator of parameter in KB1
    !                 1: Set indicator of parameter in KB1 given name
    !     OUTPUT:
    !        KB1 :    Block 1 of GRIB message
    !        YANAME : Name of indicator of parameter
    !        KERR :   Error indicator:
    !                 0: No error
    !                 1: KDIR not 0 or 1
    !                 2: No translation
    !
    implicit none
    integer:: jtbl,i
    parameter(jtbl=255)
    character(len=*):: yaname
    character(len=50):: yotbl(jtbl)
    integer::   kb1(*), kdir, kerr
    data yotbl(  1) /'pressure'/
    data yotbl(  2) /'pressure reduced to msl'/
    data yotbl(  3) /'pressure tendency'/
    data yotbl(  4)            /  ' '/
    data yotbl(  5)            /  ' '/
    data yotbl(  6) /'geopotential'/
    data yotbl(  7) /'geopotential height'/
    data yotbl(  8) /'geometric height'/
    data yotbl(  9)            /  ' '/
    data yotbl( 10)            /  ' '/
    data yotbl( 11) /'temperature'/
    data yotbl( 12) /'virtual temperature'/
    data yotbl( 13) /'potential temperature'/
    data yotbl( 14) /'pseudo-adiabatic potential temperature'/
    data yotbl( 15) /'maximum temperature'/
    data yotbl( 16) /'minimum temperature'/
    data yotbl( 17) /'dew point temperature'/
    data yotbl( 18) /'dew point depression(or deficit)'/
    data yotbl( 19) /'lapse rate'/
    data yotbl( 20)            /  ' '/
    data yotbl( 21) /'radar spectra(1)'/
    data yotbl( 22) /'radar spectra(2)'/
    data yotbl( 23) /'radar spectra(3)'/
    data yotbl( 24)            /  ' '/
    data yotbl( 25) /'temperature anomaly'/
    data yotbl( 26) /'pressure anomaly'/
    data yotbl( 27) /'geopotential height anomaly'/
    data yotbl( 28) /'wave spectra(1)'/
    data yotbl( 29) /'wave spectra(2)'/
    data yotbl( 30) /'wave spectra(3)'/
    data yotbl( 31) /'wind direction'/
    data yotbl( 32) /'wind speed'/
    data yotbl( 33) /'u-component of wind'/
    data yotbl( 34) /'v-component of wind'/
    data yotbl( 35) /'stream function'/
    data yotbl( 36) /'velocity potential'/
    data yotbl( 37)            /  ' '/
    data yotbl( 38) /'sigma coord. vertical velocity'/
    data yotbl( 39) /'vertical velocity'/
    data yotbl( 40) /'vertical velocity'/
    data yotbl( 41) /'absolute vorticity'/
    data yotbl( 42) /'absolute divergence'/
    data yotbl( 43) /'relative vorticity'/
    data yotbl( 44) /'relative divergence'/
    data yotbl( 45) /'vertical u-component shear'/
    data yotbl( 46) /'vertical v-component shear'/
    data yotbl( 47) /'direction of current'/
    data yotbl( 48) /'speed of current'/
    data yotbl( 49) /'u-component of current'/
    data yotbl( 50) /'v-component of current'/
    data yotbl( 51) /'specific humidity'/
    data yotbl( 52) /'relative humidity'/
    data yotbl( 53) /'humidity mixing ratio'/
    data yotbl( 54) /'precipitable water'/
    data yotbl( 55) /'vapor pressure'/
    data yotbl( 56) /'saturation deficit'/
    data yotbl( 57) /'evaporation'/
    data yotbl( 58)            /  ' '/
    data yotbl( 59) /'precipitation rate'/
    data yotbl( 60) /'thunderstorm probability'/
    data yotbl( 61) /'total precipitation'/
    data yotbl( 62) /'large scale precipitation'/
    data yotbl( 63) /'convective precipitation'/
    data yotbl( 64) /'snowfall rate water equivalent'/
    data yotbl( 65) /'water equivalent of accumulated snow depth'/
    data yotbl( 66) /'snow depth'/
    data yotbl( 67) /'depth of the pbl'/
    data yotbl( 68) /'transient thermocline depth'/
    data yotbl( 69) /'main thermocline depth'/
    data yotbl( 70) /'main thermocline anomaly'/
    data yotbl( 71) /'total cloud cover'/
    data yotbl( 72) /'convective cloud cover'/
    data yotbl( 73) /'low cloud cover'/
    data yotbl( 74) /'medium cloud cover'/
    data yotbl( 75) /'high cloud cover'/
    data yotbl( 76) /'cloud water'/
    data yotbl( 77)            /  ' '/
    data yotbl( 78) /'convective snowfall'/
    data yotbl( 79) /'large scale snowfall'/
    data yotbl( 80)            /  ' '/
    data yotbl( 81) /'fraction of land(1=land; 0=sea)'/
    data yotbl( 82) /'deviation of sea level from mean'/
    data yotbl( 83) /'surface roughness'/
    data yotbl( 84) /'albedo'/
    data yotbl( 85) /'soil temperature'/
    data yotbl( 86) /'soil moisture content'/
    data yotbl( 87) /'vegetation'/
    data yotbl( 88) /'salinity'/
    data yotbl( 89) /'density'/
    data yotbl( 90)            /  ' '/
    data yotbl( 91) /'ice concentration(ice=1; no ice=0)'/
    data yotbl( 92) /'ice thickness'/
    data yotbl( 93) /'direction of ice drift'/
    data yotbl( 94) /'speed of ice drift'/
    data yotbl( 95) /'u-component of ice drift'/
    data yotbl( 96) /'v-component of ice drift'/
    data yotbl( 97) /'ice growth'/
    data yotbl( 98) /'ice divergence'/
    data yotbl( 99)            /  ' '/
    data yotbl(100) /'sign. height of comb. wind waves and swell'/
    data yotbl(101) /'direction of wind waves'/
    data yotbl(102) /'significant height of wind waves'/
    data yotbl(103) /'mean period of wind waves'/
    data yotbl(104) /'direction of swell waves'/
    data yotbl(105) /'significant height of swell waves'/
    data yotbl(106) /'mean period of swell waves'/
    data yotbl(107) /'primary wave direction'/
    data yotbl(108) /'primary wave mean period'/
    data yotbl(109) /'secondary wave direction'/
    data yotbl(110) /'secondary wave mean period'/
    data yotbl(111) /'net shortwave radiation(surface)'/
    data yotbl(112) /'net longwave radiation(surface)'/
    data yotbl(113) /'net shortwave radiation(top of atmos.)'/
    data yotbl(114) /'net longwave radiation(top of atmos.)'/
    data yotbl(115) /'downwelling long wave radiation'/
    data yotbl(116) /'short wave radiation'/
    data yotbl(117) /'downwelling global sw radiation'/
    data yotbl(118)            /  ' '/
    data yotbl(119)            /  ' '/
    data yotbl(120)            /  ' '/
    data yotbl(121) /'latent heat flux'/
    data yotbl(122) /'sensible heat flux'/
    data yotbl(123) /'boundary layer dissipation'/
    data yotbl(124) /'momentum flux u-component '/
    data yotbl(125) /'momentum flux v-component '/
    data yotbl(126)            /  ' '/
    data yotbl(127) /'image data'/
    data( yotbl(i), i=128,169) /42*' '/
    data yotbl(170) /'forest: clearing; fraction'/
    data yotbl(171) /'forest: coniferious, dense; fraction'/
    data yotbl(172) /'forest: coniferious, sparse; fraction'/
    data yotbl(173) /'forest: decidious, dense; fraction'/
    data yotbl(174) /'forest: decidious, sparse; fraction'/
    data yotbl(175) /'forest: mixed(= conif. & decid.); fraction'/
    data yotbl(176) /'forest: bush land; fraction'/
    data yotbl(177)            /  ' '/
    data yotbl(178)            /  ' '/
    data yotbl(179) /'forest: undefined forest; fraction'/
    data yotbl(180)/'open land: agriclt/vegtatd, non-irrigated; frctn'/
    data yotbl(181) /'open land: bare mountain; fraction'/
    data yotbl(182) /'open land: bare soil; fraction'/
    data yotbl(183) /'open land: marsh, wet; fraction'/
    data yotbl(184) /'open land: marsh, dry; fraction'/
    data yotbl(185) /'open land: ice/snow(permanent); fraction'/
    data yotbl(186) /'open land: agriclt/vegtatd, irrigated; fraction'/
    data yotbl(187) /'open land: grass land; fraction'/
    data yotbl(188) /'open land: urban area; fraction'/
    data yotbl(189) /'open land: undefined open land; fraction'/
    data yotbl(190)            /  ' '/
    data yotbl(191)            /  ' '/
    data yotbl(192)            /  ' '/
    data yotbl(193)            /  ' '/
    data yotbl(194)            /  ' '/
    data yotbl(195) /'soil type(code number : 1,2 .... 9)'/
    data yotbl(196) /'lake; fraction'/
    data yotbl(197) /'forest; fraction'/
    data yotbl(198) /'open land; fraction'/
    data yotbl(199)            /  ' '/
    data yotbl(200) /'turbulence kinetic energy'/
    data yotbl(201)            /  ' '/
    data yotbl(202)            /  ' '/
    data yotbl(203)            /  ' '/
    data yotbl(204)            /  ' '/
    data yotbl(205)            /  ' '/
    data yotbl(206)            /  ' '/
    data yotbl(207)            /  ' '/
    data yotbl(208)            /  ' '/
    data yotbl(209)            /  ' '/
    data yotbl(210) /'dtdt due to all parametrizations '/
    data yotbl(211) /'dtdt due to turbulence '/
    data yotbl(212) /'dtdt due to radiation '/
    data yotbl(213)            /  ' '/
    data yotbl(214)            /  ' '/
    data yotbl(215)            /  ' '/
    data yotbl(216)            /  ' '/
    data yotbl(217)            /  ' '/
    data yotbl(218)            /  ' '/
    data yotbl(219)            /  ' '/
    data yotbl(220) /'dqdt due to all parametrizations '/
    data yotbl(221) /'dqdt due to turbulence '/
    data yotbl(222)            /  ' '/
    data yotbl(223)            /  ' '/
    data yotbl(224)            /  ' '/
    data yotbl(225)            /  ' '/
    data yotbl(226)            /  ' '/
    data yotbl(227)            /  ' '/
    data yotbl(228)            /  ' '/
    data yotbl(229)            /  ' '/
    data yotbl(230) /'dcwdt due to all parametrizations '/
    data yotbl(231) /'dcwdt due to turbulence '/
    data yotbl(232)            /  ' '/
    data yotbl(233)            /  ' '/
    data yotbl(234)            /  ' '/
    data yotbl(235)            /  ' '/
    data yotbl(236)            /  ' '/
    data yotbl(237)            /  ' '/
    data yotbl(238)            /  ' '/
    data yotbl(239)            /  ' '/
    !
    data( yotbl(i), i=240,255) /16*' '/

    !     Preset return values
    !
    kerr = 0
    if(kdir == 0) yaname = ' '
    if(kdir == 1) kb1(9) = 0
    !
    !     return either name or indicator of parameter.
    !
    if(kdir == 0 .and. kb1(9) > 0 .and. kb1(9) <= jtbl) then
       if(kb1(4)==1.or.kb1(4)==2)then
          yaname = yotbl(kb1(9))
       elseif(kb1(4)==254)then       ! some hirlam extensions
          if(kb1(9)==226)then
             yaname='longitudes(table 2 vn 254)'
          elseif(kb1(9)==225)then
             yaname='latitudes(table 2 vn 254)'
          elseif(kb1(9)==1)then
             yaname='dates(table 2 vn 254)'
          elseif(kb1(9)==2)then
             yaname='times(table 2 vn 254)'
          elseif(kb1(9)==11)then
             yaname='obs subtypes(table 2 vn 254)'
          else
             yaname='?? - unknown entry in tbl2 vn 254'
          endif
       else
          yaname='?? - unknown table 2 version'
       endif
    else if(kdir == 1) then
       do 100 i = 1, jtbl
          if(yotbl(i) == yaname) then
             kb1(9) = i
             goto 101
          endif
100    enddo
101    continue
    else
       kerr = 1
    endif

    !     2. error checks.
    !
    if(kdir == 0 .and. yaname == ' ') kerr = 2
    if(kdir == 1 .and. kb1(9) == 0)   kerr = 2
    !
    return
  end subroutine gbctb2



  subroutine gbe2h1(ksec0,ksec1,kb1)
    !     
    !     subroutine gbe2h1 swaps the values of ksec1 back to kb1.
    !     
    !     input :
    !     ksec0 : array with the values of block 0
    !     ksec1 : array with the values of block 1
    !     
    !     output :
    !     kb1 : array with the shifted values of block 1 of a grib-code-file
    !     
    !     
    !     variables
    !     
    implicit none
    integer:: ksec0(*) , ksec1(*) , kb1(*)
    !     
    !     setting of kb1
    !     
    if( ksec0(2) == 0 ) then
       kb1(1) = 24
    else
       kb1(1) = 28
    endif
    kb1(4) = ksec1(1)
    kb1(5) = ksec1(2)
    kb1(6) = ksec1(3)
    kb1(7) = ksec1(4)
    kb1(8) = ksec1(5)
    kb1(9) = ksec1(6)
    kb1(10) = ksec1(7)
    kb1(11) = ksec1(8)
    kb1(12) = ksec1(9)
    kb1(13) = ksec1(10) + ksec1(21) * 100 - 100
    kb1(14) = ksec1(11)
    kb1(15) = ksec1(12)
    kb1(16) = ksec1(13)
    kb1(17) = ksec1(14)
    kb1(18) = ksec1(15)
    kb1(19) = ksec1(16)
    kb1(20) = ksec1(17)
    kb1(21) = ksec1(18)
    kb1(22) = ksec1(19)

    !     this is the end.
    !     
    return
  end subroutine gbe2h1
  subroutine gbe2h2(ksec2,kb2,psec2,pw2)
    !     
    !     subroutine gbe2h2 swaps the values from ksec2 back to kb2
    !     and from psec2 back to pw2.
    !     
    !     input :
    !     ksec2 : array with integer values of block 2
    !     psec2 : array with real values of block 2
    !     
    !     output :
    !     kb2 : array with the shifted values of block 2 of grib-code-file
    !     pw2 : array with the real levels of block2(=vertical levels)
    !     
    !     
    !     variables
    !     
    implicit none
    integer:: inc,i
    integer:: ksec2(*) , kb2(*) 
    real(kind=realgribkind) psec2(*) , pw2(*)

    !     setting of kb2

    inc = 0

    if(ksec2(1)/=5)then
       if(ksec2(1)==10.or.ksec2(1)==14.or.ksec2(1)==60) inc = 10
       if(ksec2(1)==20.or.ksec2(1)==24.or.ksec2(1)==70) inc = 10
       if(ksec2(1)==30.or.ksec2(1)==34.or.ksec2(1)==80) inc = 20
       kb2(1) = 32 + inc +(ksec2(12)*4) +(ksec2(17)*2*ksec2(3))
       kb2(4) = 0
       kb2(6) = ksec2(1)
       kb2(7) = ksec2(2)
       kb2(9) = ksec2(3)
       kb2(11) = ksec2(4)
       kb2(14) = ksec2(5)
       kb2(17) = ksec2(6)+ksec2(18)+ksec2(19)
       kb2(18) = ksec2(7)
       kb2(21) = ksec2(8)
       kb2(24) = ksec2(9)
       kb2(26) = ksec2(10)
       kb2(28) = ksec2(11)
       kb2(33) = ksec2(13)
       kb2(36) = ksec2(14)
       if(inc == 10 ) kb2(39) = 1
       if(inc == 20 ) kb2(49) = 2
       kb2(40) = inc / 10
       !     goto 200
    else
       !      endif
       ! 105  continue
       kb2(1) = 32 + inc +(ksec2(12)*4) +(ksec2(17)*2*ksec2(3))
       kb2(4) = 0
       kb2(6) = ksec2(1)
       kb2(7) = ksec2(2)
       kb2(9) = ksec2(3)
       kb2(11) = ksec2(4)
       kb2(14) = ksec2(5)
       kb2(17) = ksec2(6)+ksec2(18)+ksec2(19)
       kb2(18) = ksec2(7)
       kb2(21) = ksec2(9)
       kb2(24) = ksec2(10)
       kb2(27) = ksec2(13)
       kb2(28) = ksec2(11)
       kb2(40) = inc / 10
    endif
    !      goto 200
    !200 continue

    !     vertical levels need shifting too(ksec2(12)= number of vertical leve
    !     
    do i=1,ksec2(12)
       pw2(i+kb2(40)) = psec2(i+10)
    enddo
    do  i=1,kb2(40)
       pw2(i) = psec2(i)
    enddo
    !     
    !     this is the end.
    !     
    return
  end subroutine gbe2h2


  subroutine gbe2h4(ksec2,ksec4,kb4)
    !     
    !     subroutine gbe2h4 swaps the values of ksec4 back to kb4.
    !     
    !     input :
    !     ksec2 : array with values of block2
    !     ksec4 : array with integer values of block4
    !     
    !     output :
    !     kb4 : array with the values of block4 of a grib-code-file
    implicit none
    integer:: irep
    real(kind=realgribkind) hkb
    integer:: ksec2(*) , ksec4(*) , kb4(*)
    if( ksec2(1)==50.or.ksec2(1)==60.or. &
         ksec2(1)==70.or.ksec2(1)==80) then
       kb4(4) = 128
       irep = 1
    else
       kb4(4) = 0
       irep = 0
    endif
    hkb = real( ksec4(2) *( ksec4(1) - irep ),realgribkind ) / 8._realgribkind + real(irep * 4 + 11,realgribkind)
    kb4(1) = int(hkb + 0.5_realgribkind)
    kb4(5) = ksec4(33)
    kb4(7) = 0
    kb4(11) = ksec4(2)

    !     this is the end.
    !     
    return
  end subroutine gbe2h4


  subroutine gbg2fc(kb2,k,plon,plat,px,py,kxy,kinlon,kinlat)
    !
    !     gbg2fc: obtain either coordinates of all gridpoints(k=0) or
    !             gridpoints of all coordinates(k>0)
    !
    ! input:
    !     kb2:    gribcode block2 to define grid
    !     k:      if -1,calculate kinlon and kinlat only
    !             if 0, calculate plon,plat for all gridpoints
    !             if >0, calculate px,py,kxy for all points plon,plat
    !            (1....k)
    !     plon,plat:(only if k>0) coordinates in grid units of points

    !             1....k.
    !              nb. gridunits = degrees for lat/lon, km for polar/stereo
    !
    ! output:
    !     k:     (only if 0 on input): nr of gridpoints
    !     plon,plat(if k==0 on input): coordinates of gridpoints
    !     px,py,kxy(if k>0 on input): positions of points(plon,plat)
    !               within grid: nearest gridpoint is kxy, distance to it
    !               is px,py for x,y resp. 'nearest' is: nearest to the
    !               west when scanning from w to e, else nearest to the
    !               east. similarly, it is nearest to the south, when
    !               scanning from s to n, else nearest to the n
    !     kinlon: distance between consecutive longitude points in the data
    !             array, positive when scanning w to e, else negative
    !     kinlat: distance between consecutive latitude points in the data
    !             array, positive when scanning s to n, else negative
    !
    ! notes:
    !     1.  coordinates are defined as follows:
    !         a) for lat/lon grids:
    !            longitude(x) and latitude(y) in degrees
    !         b) for rotated and/or stretched lat/lon grids:
    !            longitude(x) and latitude(y) as if plain lat/lon grids
    !         c) for polar-stereographic projections(projection type 5):
    !            x and y in metres along the user axes(the y-axis
    !            being parallel to lov, the origin being at the
    !            point lo1,la1, where lov, lo1 and la1 are as in the
    !            block 2 octets 18-20,14-16 and 11-13 resp.)
    !            x and y are to true scale at +/- 60 degrees
    !         d) for polar-stereographic projections(projection type 2):
    !           (defined by g.j.cats in knmi memorandum dm-88-03)
    !            x and y in km along the user axes,(the y-axis
    !            being parallel to labda y-(octets 31-32) and the
    !            origin being at the pole on the projection plane)
    !            x and y are to true scale at latitude phi p(octets 29-34)
    !
    !                                                     g.j.cats 23 dec 87
    !
    !     original version(burroughs/cray/convex)        g.j.cats 23 dec 87
    !     revision to include proj.type 5                 g.j.cats 17 jan 89
    !

    implicit none
    logical lops2,loback,lops5,loncon,lops
    integer:: kb2(*),kxy(*)
    real(kind=realgribkind) plon(*),plat(*),px(*),py(*)
    real(kind=realgribkind) zlate,zlato,zlone,zlono,zscale,zlat,zx,zy,zdlat,zdlon,zlon
    real(kind=realgribkind) zs2n,zw2e
    integer:: i0lon,jlon,jlat,i0lat,k,ix,jgp,iy,kinlat,kinlon
    integer:: inlon,inlat,ilat,i,ilon
    !     determine if polar stereographic

    lops2=kb2(6)==2
    lops5=kb2(6)==5
    lops=lops2.or.lops5
    !     obtain scanning mode and number of gridpoints

    ilon=kb2(7)
    ilat=kb2(9)
    zw2e=1._realgribkind
    zs2n=-1._realgribkind
    loncon=.true.
    i=kb2(28)

    !     correct possible counting errors by a shift over 5 bits:
    if(i<8)i=i*32
    if(ngbgrb(i,1)==1)zw2e=-1._realgribkind
    if(ngbgrb(i,2)==1)zs2n=1._realgribkind
    if(ngbgrb(i,3)==1)loncon=.false.
    if(loncon)then
       inlon=1
       inlat=ilon
    else
       inlon=ilat
       inlat=1
    endif
    kinlon=sign(inlon,nint(zw2e))
    kinlat=sign(inlat,nint(zs2n))
    if(k==-1)return

    !     1.3 coordinates of origin, extreme point

    zscale=0.001_realgribkind
    if(lops5)then
       zscale=1._realgribkind
       zlono=0._realgribkind
       zlone=zlono+zscale*real(kb2(21)*(kb2(7)-1),realgribkind)*zw2e
       zlato=0._realgribkind
       zlate=zlato+zscale*real(kb2(24)*(kb2(9)-1),realgribkind)*zs2n
    else
       if(lops2)zscale=0.1_realgribkind
       zlono=zscale*real(kb2(14),realgribkind)
       zlato=zscale*real(kb2(11),realgribkind)
       zlone=zscale*real(kb2(21),realgribkind)
       zlate=zscale*real(kb2(18),realgribkind)
    endif

    !     1.4 lon/lat increments in grid-points and -units(lat/lon at 180w?

    loback=.false.
    if(.not.lops)then
       if(zlono> 180.001_realgribkind)zlono=zlono-360._realgribkind
       zlono=min(zlono,180._realgribkind)
       if(zlone> 180.001_realgribkind)zlone=zlone-360._realgribkind
       zlone=min(zlone,180._realgribkind)
    endif
    if(.not.lops.and.(zlono-zlone)*zw2e>0._realgribkind)then
       zlone=zlone+360._realgribkind*zw2e
       loback=.true.
    endif
    if(ilon>1)then
       zdlon=(zlone-zlono)/real(ilon-1,realgribkind)
    else
       zdlon=1._realgribkind
    endif
    if(ilat>1)then
       zdlat=(zlate-zlato)/real(ilat-1,realgribkind)
    else
       zdlat=1._realgribkind
    endif

    !     determine lat lon for all gridpoints(i.e., if k==0)

    if(k==0)then

       !     initial longitudes in gridpoints and coordinates

       zlon=zlono
       i0lon=0
       !     loop over longitudes, up to section 2.8

       do 280 jlon=1,ilon

          !     initial latitude

          zlat=zlato
          i0lat=1
          !     loop over latitudes

          do 241 jlat=1,ilat
             plat(i0lat+i0lon)=zlat
             plon(i0lat+i0lon)=zlon
             zlat=zlat+zdlat
             i0lat=i0lat+inlat
241       enddo

          !     update longitudes

          zlon=zlon+zdlon
          if(.not.lops)then
             if(zlon> 180._realgribkind)zlon=zlon-360._realgribkind
             if(zlon<=-180._realgribkind)zlon=zlon+360._realgribkind
          endif
          i0lon=i0lon+inlon

          !     end loop over lontiudes

280    enddo

       !     update k

       k=ilat*ilon

       !       calculate positions of plon,plat in grid, i.e. k>0

    else

       !      loop over gridpoints

       do 390 jgp=1,k
          zlon=plon(jgp)

          !     allow longitude across 180w if lat/lon grid

          if(loback.and.(zlon-zlono)*zw2e<0._realgribkind)zlon=zlon+360._realgribkind*zw2e

          !     integer and fractional positions

          zx=(zlon-zlono)/zdlon
          ix=max(min(int(zx),ilon-2),0)
          px(jgp)=zx-real(ix,realgribkind)
          zy=(plat(jgp)-zlato)/zdlat
          iy=max(min(int(zy),ilat-2),0)
          py(jgp)=zy-real(iy,realgribkind)
          kxy(jgp)=ix*inlon+iy*inlat+1


390    enddo

    endif
    return
  end subroutine gbg2fc



  subroutine gbh2e1(kb1,ksec1)
    !     
    !     subroutine gbh2e1 swaps the values of kb1 to ksec1.
    !     
    !     input :
    !     kb1 : array with the values of block 1 of a grib - code - file
    !     
    !     output :
    !     ksec1 : array with 25 words
    !     
    !     
    !     ksec1(1)  : version number of code table 2(=1 for hirlam-files)
    !     ksec1(2)  : identification of centre(=96 for hirlam-files)
    !     ksec1(3)  : identification number
    !     ksec1(4)  : grid definition
    !     ksec1(5)  : flag indication
    !     ksec1(6)  : indicator of parameter
    !     ksec1(7)  : indicator of type of level
    !     ksec1(8)  : height, pressure, etc. of level
    !     ksec1(9)  : height, pressure, etc. of level
    !     ksec1(10) : year
    !     ksec1(11) : month
    !     ksec1(12) : day
    !     ksec1(13) : hour
    !     ksec1(14) : minute
    !     ksec1(15) : indicator of unit of time
    !     ksec1(16) : period of time
    !     ksec1(17) : period of time
    !     ksec1(18) : time range indicator
    !     ksec1(19) : number included in average
    !     ksec1(20) : number missing from average(=0)
    !     ksec1(21) : century of reference time of data
    !     ksec1(22) : set to zero(=0)
    !     ksec1(23) : decimal scale factor(=0)
    !     ksec1(24) : set to zero(=0)
    !     ksec1(25) : set to zero(=0)
    !     
    !     variables
    !     
    implicit none
    integer:: kb1(*) , ksec1(*)
    !     setting of ksec1
    ksec1(1) = kb1(4)
    ksec1(2) = kb1(5)
    ksec1(3) = kb1(6)
    ksec1(4) = kb1(7)
    ksec1(5) = kb1(8)
    ksec1(6) = kb1(9)
    ksec1(7) = kb1(10)
    ksec1(8) = kb1(11)
    ksec1(9) = kb1(12)
    ksec1(10) = mod(kb1(13), 100)
    ksec1(11) = kb1(14)
    ksec1(12) = kb1(15)
    ksec1(13) = kb1(16)
    ksec1(14) = kb1(17)
    ksec1(15) = kb1(18)
    ksec1(16) = kb1(19)
    ksec1(17) = kb1(20)
    ksec1(18) = kb1(21)
    ksec1(19) = kb1(22)
    ksec1(20) = 0
    ksec1(21) = kb1(13) / 100
    !     the following lines because 1999 is 20th, 2000 is 20th, 2001 is 21st
    if(  ksec1(10) == 0 ) then
       ksec1(10) = 100
    else
       ksec1(21) = ksec1(21) + 1
    endif
    !     upward compatibility: if before ad 100, assume 20th century
    if( kb1(13) < 100 ) ksec1(21) = 20
    ksec1(22) = 0
    ksec1(23) = 0
    ksec1(24) = 0
    ksec1(25) = 0

    return
  end subroutine gbh2e1


  subroutine gbh2e2(kb2,ksec2,pw2,psec2)
    !
    !  subroutine gbh2e2 swaps the values of kb2 to ksec2
    !   and from pw2 to psec2
    !
    !  input :
    !     kb2 : array with the integer values of block 2 of a grib-code-file
    !     pw2 : array with the real values of block 2(=the vertical levels)
    !
    !  output :
    !     ksec2 : array with 22(+n) words
    !     psec2 : array with(ksec2(12) + 10) words
    !
    !
    !  ksec2(1)  : data representation type(=10 for hirlam-files)
    !  ksec2(2)  : number of points along a parallel
    !  ksec2(3)  : number of points along a meridian
    !  ksec2(4)  : latitude of first grid point
    !  ksec2(5)  : longitude of first grid point
    !  ksec2(6)  : increments flag
    !  ksec2(7)  : latitude of last grid point
    !  ksec2(8)  : longitude of last grid point
    !  ksec2(9)  : i - direction increment
    !  ksec2(10) : j - direction increment
    !  ksec2(11) : scanning mode flags
    !  ksec2(12) : number of vertical coordinate parameters
    !  ksec2(13) : latitude of southern pole of rotation
    !  ksec2(14) : longitude of southern pole of rotation
    !  ksec2(15) : latitude of pole of stretching
    !  ksec2(16) : longitude of pole of stretching
    !  ksec2(17) : regular of quasi-regular grid(=0)
    !  ksec2(18) : earth flag(=0)
    !  ksec2(19) : components flag(=0)
    !  ksec2(20) : set to zero(=0)
    !  ksec2(21) : set to zero(=0)
    !  ksec2(22) : set to zero(=0)
    !
    !  variables
    !
    implicit none
    integer:: inc,i
    integer:: kb2(*) , ksec2(*) 
    real(kind=realgribkind) pw2(*) , psec2(*)
    !
    !  setting of ksec2
    !
    ksec2(1) = kb2(6)
    inc=0
    !        if(ksec2(1)==5)goto 105
    if(ksec2(1)/=5)then
       if(ksec2(1)==10.or.ksec2(1)==14.or.ksec2(1)==60) inc = 1
       if(ksec2(1)==20.or.ksec2(1)==24.or.ksec2(1)==70) inc = 1
       if(ksec2(1)==30.or.ksec2(1)==34.or.ksec2(1)==80) inc = 2
       ksec2(2) = kb2(7)
       ksec2(3) = kb2(9)
       ksec2(4) = kb2(11)
       ksec2(5) = kb2(14)
       ksec2(6) =(kb2(17)/128)*128
       ksec2(7) = kb2(18)
       ksec2(8) = kb2(21)
       ksec2(9) = kb2(24)
       ksec2(10) = kb2(26)
       ksec2(11) = kb2(28)
       ksec2(12) =( kb2(1) - 32 - inc * 10 ) / 4
       ksec2(13) = kb2(33)
       ksec2(14) = kb2(36)
       ksec2(15) = kb2( 33 +( inc - 1 ) * 10 )
       ksec2(16) = kb2( 36 +( inc - 1 ) * 10 )
       ksec2(17) = 0
       ksec2(18) = 0
       ksec2(19) =  mod((kb2(17)/8),2)*8
       ksec2(20) = 0
       ksec2(21) = 0
       ksec2(22) = 0
    else
       ksec2(2) = kb2(7)
       ksec2(3) = kb2(9)
       ksec2(4) = kb2(11)
       ksec2(5) = kb2(14)
       ksec2(6) =(kb2(17)/128)*128
       ksec2(7) = kb2(18)
       ksec2(8) = 0
       ksec2(9) = kb2(21)
       ksec2(10) = kb2(24)
       ksec2(11) = kb2(28)
       ksec2(12) =( kb2(1) - 32 - inc * 10 ) / 4
       ksec2(13) = kb2(27)
       ksec2(14) = 0
       ksec2(15) = 0
       ksec2(16) = 0
       ksec2(17) = 0
       ksec2(18) = 0
       ksec2(19) =  mod((kb2(17)/8),2)*8
       ksec2(20) = 0
       ksec2(21) = 0
       ksec2(22) = 0
    endif
    !  vertical levels are shifted too(ksec2(12)=number of vertical levels)
    !
    do  i=1,ksec2(12)
       psec2(i+10) = pw2(i+inc)
    enddo
    do  i=1,inc
       psec2(i) = pw2(i)
    enddo
    !
    !  this is the end.
    !
    return
  end subroutine gbh2e2


  subroutine gbh2e4(kb4,ksec4,mdat)
    !     
    !     subroutine gbh2e4 swaps the values of kb4 to ksec4.
    !     
    !     input :
    !     kb4 : array with the integer values of block 4 of a grib-code-file
    !     mdat : number of data values to be packed
    !     
    !     output :
    !     ksec4 : array with 42 words
    !     
    !     
    !     ksec4(1) : number of data values in psec4
    !     ksec4(2) : number of bits used for each packed value
    !     ksec4(3) : type of packing(=0 for hirlam-files)
    !     ksec4(4) : type of data(=0 for hirlam-files)
    !     ksec4(5) : type of grid points(=0 for hirlam-files)
    !     ksec4(6) : secondary bit-maps(=0 for hirlam-files)
    !     ksec4(7) : second order values(=0 for hirlam-files)
    !     ksec4(8) : number of bits for second order values(=0 for hirlam-file
    !     ksec4(9) - ksec4(42) : set to zero(=0)
    !     
    !     variables
    !     
    implicit none
    integer:: i,mdat
    integer:: kb4(*) , ksec4(*)
    !     
    !     setting of ksec4
    !     
    ksec4(1) = mdat
    ksec4(2) = kb4(11)
    ksec4(3) = 0
    ksec4(4) = 0
    ksec4(5) = 0
    ksec4(6) = 0
    ksec4(7) = 0
    ksec4(8) = 0

    do  i=9,42
       ksec4(i) = 0
    enddo
    return
  end subroutine gbh2e4



  subroutine gbset0(ksec0)
    !     
    !     subroutine gbset0 sets the values of section 0.
    !     
    !     input : none
    !     
    !     output :
    !     ksec0 : array with 2 words
    !     
    !     
    !     ksec0(1) : number of octets in grib message
    !    (will be finally set by gribex)
    !     ksec0(2) : grib edition number(always 1 in coding)
    !     
    !     
    !     declaration of the variables
    !     
    implicit none
    !      real(kind=realgribkind) ksec0(*)
    integer::ksec0(2)

    ksec0(1) = 9999
    ksec0(2) = 1

    return
  end subroutine gbset0

  subroutine grchk2(ksec1,ksec2,kret)
    implicit none
    !     
    !**** grchk2 - check parameters for section 2 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     check parameters for section 2 of grib code against
    !     valid values for grib edition 1.
    !     
    !**   interface.
    !     ----------
    !     
    !     call grchk2(ksec1,ksec2,kret)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec1      - array containing integer parameters for
    !     section 1 of grib code.
    !     
    !     ksec2      - array containing integer parameters for
    !     section 2 of grib code.
    !     
    !     
    !     output parameters.
    !     ------------------
    !     
    !     kret       - return code.

    !     0   , no error encountered.
    !     1   , error in grib code parameter.
    !     
    !     method.
    !     -------
    !     
    !     values checked against current code/flag tables
    !     and against maximum or minimum permitted values.
    !     they are also checked against the current status
    !     of the implementation of gribex and ecmwf usage.
    !     
    !     externals.
    !     ----------
    !     
    !     none.
    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes for grib code.
    !     
    integer:: j201
    integer:: j203
    integer:: j301
    integer:: j401

    integer:: kret

    integer,parameter::jp6=22
    integer,parameter::jp6x=14
    integer,parameter::jp8=8

    integer:: ksec1(*)
    integer:: ksec2(*)

    integer:: itab6(jp6)
    integer:: itab6x(jp6x)
    integer:: itab8(jp8)


    !     valid values given in code table 6.
    data itab6 /0,1,2,3,4,5,6,7,8,9,10,13,14,20,24,30,34,50,60,70,80,90/

    !     code table 6 values currently(04.10.91) supported by gribex.
    data itab6x /0,4,5,10,14,20,24,30,34,50,60,70,80,90/

    !     valid values given in code table 8.
    data itab8 /0,128,64,192,32,160,96,224/

    kret = 0


    ksec2(9)=1
    ksec2(10)=1
    !         number of vertical coordinate parameters.

    if(ksec2(12)<0.or.ksec2(12)>255)then
       kret = 1
       write(*,9019) ksec2(12)
    endif
    !         check data representation type.

    do 201 j201=1,jp6
       if(ksec2(1)==itab6(j201)) go to 202
201 enddo

    kret = 1
    write(*,9001) ksec2(1)

202 continue

    !         check data representation type currently supported.
    do 203 j203=1,jp6x
       if(ksec2(1)==itab6x(j203)) go to 204
203 enddo

    kret = 1
    write(*,9002) ksec2(1)
    go to 900

204 continue

    !         earth flag.
    if(ksec2(18)/=0.and.ksec2(18)/=64)then
       kret = 1
       write(*,9005) ksec2(18)
    endif

    !         check ecmwf usage.(0 except for space view perspective)

    if(ksec2(18)/=0.and.ksec1(2)==98.and.ksec2(1)/=90)then
       kret = 1
       write(*,9007)
    endif
    !          components flag.

    if(ksec2(19)/=0.and.ksec2(19)/=8)then
       kret = 1
       write(*,9006) ksec2(19)
    endif

    !         section 3. checks on latitude/longitude grids.


    if(ksec2(1)==0.or.ksec2(1)==10.or.ksec2(1)==20.or.ksec2(1)==30)then
       
       if(ksec2(2)<1.or.ksec2(2)>65535)then ! number of points along a parallel.
          kret = 1
          write(*,9022) ksec2(2)
       endif
       if(ksec2(3)<1.or.ksec2(3)>65535)then! number of points along a meridian.
          kret = 1
          write(*,9023) ksec2(3)
       endif
       if(ksec2(4)<-90000.or.ksec2(4)>90000)then! latitude of first grid point.
          kret = 1
          write(*,9015) ksec2(4)
       endif
       if(ksec2(5)<-360000.or.ksec2(5)>360000)then!longitude of first grid point.
          kret = 1
          write(*,9016) ksec2(5)
       endif
       if(ksec2(6)/=0.and.ksec2(6)/=128)then!resolution flag.
          kret = 1
          write(*,9003) ksec2(6)
       endif
       if(ksec2(7)<-90000.or.ksec2(7)>90000)then! latitude of last grid point.
          kret = 1
          write(*,9020) ksec2(7)
       endif
       if(ksec2(8)<-360000.or.ksec2(8)>360000)then! longitude of last grid point.
          kret = 1
          write(*,9021) ksec2(8)
       endif
       if(ksec2(6)==128)then ! direction increments, if included.
          if(ksec2(9)<1.or.ksec2(9)>65535)then
             kret = 1
             write(*,9024) ksec2(9)
          endif
          if(ksec2(10)<1.or.ksec2(10)>65535)then
             kret = 1
             write(*,9025) ksec2(10)
          endif
       endif
       do 301 j301=1,jp8 !scanning mode flag.
          if(ksec2(11)==itab8(j301)) go to 302
301    enddo
       kret = 1
       write(*,9004) ksec2(11)

302    continue
       
       if(ksec2(17)/=0.and.ksec2(17)/=1)then!   regular / quasi-regular grid check.
          kret = 1
          write(*,9009) ksec2(17)
       endif
       if(ksec2(17)/=0)then!     check currently supported values.
          kret = 1
          write(*,9010)
       endif

       go to 900
    endif
    


    if(ksec2(1)==4.or.ksec2(1)==14.or.ksec2(1)==24 .or.ksec2(1)==34)then !  checks on gaussian grids.
       if(ksec2(4)<-90000.or.ksec2(4)>90000)then ! latitude of first grid point.
          kret = 1
          write(*,9015) ksec2(4)
       endif
       if(ksec2(5)<-360000.or.ksec2(5)>360000)then !  longitude of first grid point.
          kret = 1
          write(*,9016) ksec2(5)
       endif
       if(ksec2(7)<-90000.or.ksec2(7)>90000)then ! latitude of last grid point.
          kret = 1
          write(*,9020) ksec2(7)
       endif
       if(ksec2(8)<-360000.or.ksec2(8)>360000)then!   longitude of last grid point.
          kret = 1
          write(*,9021) ksec2(8)
       endif
       if(ksec2(6)==128)then ! i-direction increment, if included.
          if(ksec2(9)<1.or.ksec2(9)>65535)then
             kret = 1
             write(*,9024) ksec2(9)
          endif
       endif
       if(ksec2(10)<1.or.ksec2(10)>65535)then !     number of parallels beween pole and equator.
          kret = 1
          write(*,9026) ksec2(10)
       endif
       if(ksec2(6)/=0.and.ksec2(6)/=128)then ! increment flag.
          kret = 1
          write(*,9003) ksec2(6)
       endif
       do 401 j401=1,jp8 !     scanning mode flag.
          if(ksec2(11)==itab8(j401)) go to 402
401    enddo
       kret = 1
       write(*,9004) ksec2(11)

402    continue
       
       if(ksec2(17)/=0.and.ksec2(17)/=1)then !    regular / quasi-regular grid check.
          kret = 1
          write(*,9009) ksec2(17)
       endif
       if(ksec2(17)==1.and.ksec2(6)==128)then !    cross-check increments flag and quasi-regular indicator.
          kret = 1
          write(*,9011)
       endif

       go to 900
    endif

    


    if(ksec2(1)==5)then !   checks on polar stereographic data.
       if(ksec2(2)<1.or.ksec2(2)>65535)then!  number of points along x-axis.
          kret = 1
          write(*,9027) ksec2(2)
       endif
       if(ksec2(3)<1.or.ksec2(3)>65535)then !   number of points along y-axis.
          kret = 1
          write(*,9028) ksec2(3)
       endif
       if(ksec2(4)<-90000.or.ksec2(4)>90000)then ! latitude of first grid point.
          kret = 1
          write(*,9015) ksec2(4)
       endif
       if(ksec2(5)<-360000.or.ksec2(5)>360000)then!    longitude of first grid point.
          kret = 1
          write(*,9016) ksec2(5)
       endif
       if(ksec2(7)<-360000.or.ksec2(7)>360000)then !   orientation of the grid.
          kret = 1
          write(*,9017) ksec2(7)
       endif
       if(ksec2(9)<1.or.ksec2(9)>16777215)then !    grid lengths.
          kret = 1
          write(*,9029) ksec2(9)
       endif
       if(ksec2(10)<1.or.ksec2(10)>16777215)then
          kret = 1
          write(*,9030) ksec2(10)
       endif
       !      projection centre. the use of 1 in this octet is
       !     inconsistent with lambert conformal et al representation
       !     where bit 1 is set 1 to indicate north pole. !!!!!!!!!

       if(ksec2(13)/=0.and.ksec2(13)/=1)then
          kret = 1
          write(*,9018) ksec2(13)
       endif

       go to 900
    endif
    


    if(ksec2(1)==50.or.ksec2(1)==60.or.ksec2(1)==70 .or.ksec2(1)==80)then!checks on spherical harmonic data.
       if(ksec2(5)/=1)then!spectral data representation type.
          kret = 1
          write(*,9012) ksec2(5)
       endif
       if(ksec2(6)/=1.and.ksec2(6)/=2)then!  spectral data representation mode.
          kret = 1
          write(*,9013) ksec2(6)
       endif
       if(ksec2(6)/=1) write(*,9014)
       go to 900
    endif

    if(ksec2(1)==90)then! checks on space view perspective.
       go to 900
    endif

900 continue

9001 format(1h ,'grchk2 : invalid data representation type - ',i3)
9002 format(1h ,'grchk2 : unsupported data representation type - ',i3)
9003 format(1h ,'grchk2 : invalid increments flag - ',i3)
9004 format(1h ,'grchk2 : invalid scanning mode flag - ',i3)
9005 format(1h ,'grchk2 : invalid earth flag - ',i3)
9006 format(1h ,'grchk2 : invalid components flag - ',i3)
9007 format(1h ,'grchk2 : earth flag - ecmwf usage is 0.')
9008 format(1h ,'grchk2 : components flag - ecmwf usage is 0.')
9009 format(1h ,'grchk2 : invalid quasi / regular indicator - ',i3)
9010 format(1h ,'grchk2 : quasi-regular latitude/longitude grids',  &
         ' not catered for.')
9011 format(1h ,'grchk2 : quasi-regular gaussian grid cannot have',&
         ' direction increments included.')
9012 format(1h ,'grchk2 : invalid spectral representation type - ',i3)
9013 format(1h ,'grchk2 : invalid spectral representation mode - ',i3)
9014 format(1h ,'grchk2 : complex spectral representation mode ', &
         'not catered for.')
9015 format(1h ,'grchk2 : invalid latitude of first grid point - ',i10)
9016 format(1h ,'grchk2 : invalid longitude of first grid point - ',i10)
9017 format(1h ,'grchk2 : invalid orientation of the grid - ', i10)
9018 format(1h ,'grchk2 : invalid projection centre flag - ',i3)
9019 format(1h ,'grchk2 : invalid number of vertical coordinate ', &
         'parameters - ',i8)
9020 format(1h ,'grchk2 : invalid latitude of last grid point - ', i10)
9021 format(1h ,'grchk2 : invalid longitude of last grid point - ', i10)
9022 format(1h ,'grchk2 : invalid number of points along a parallel', i10)
9023 format(1h ,'grchk2 : invalid number of points along a meridian', i10)
9024 format(1h ,'grchk2 : invalid i-direction increment - ',i10)
9025 format(1h ,'grchk2 : invalid j-direction increment - ',i10)
9026 format(1h ,'grchk2 : invalid number of parallels - ',i10)
9027 format(1h ,'grchk2 : invalid number of points along x-axis', i10)
9028 format(1h ,'grchk2 : invalid number of points along y-axis', i10)
9029 format(1h ,'grchk2 : invalid x-direction grid length - ',i10)
9030 format(1h ,'grchk2 : invalid y-direction grid length - ',i10)

    return

  end subroutine grchk2
  subroutine grchk3(ksec1,ksec3,kret)
    implicit none
    !     
    ! grchk3 - check parameters for section 3 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     check parameters for section 3 of grib code against
    !     valid values for grib edition 1.
    !     
    !   interface.
    !     ----------
    !     
    !     call grchk3(ksec1,ksec3,kret)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec1      - array containing integer parameters for
    !     section 1 of grib code.
    !     
    !     ksec3      - array containing integer parameters for
    !     section 3 of grib code.
    !     
    !     
    !     output parameters.
    !     ------------------
    !     
    !     kret       - return code.

    !     0   , no error encountered.
    !     1   , error in grib code parameter.
    !     
    !     method.
    !     -------
    !     
    !     values checked against current code/flag tables
    !     and against maximum or minimum permitted values.
    !     they are also checked against the current status
    !     of the implementation of gribex and ecmwf usage.
    !     
    !     externals.
    !     ----------
    !     
    !     none.
    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes for grib code.
    !     
    !     comments.
    !     ---------
    !     
    !     routine contains sections 0 to 5 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     j. hennessy      ecmwf      16.09.91
    !     
    !     modifications.
    !     --------------
    !     
    !     j. hennessy      ecmwf      01.10.91

    integer:: kret
    integer:: ksec1(*)
    integer:: ksec3(*)


    kret = 0
    

    if(ksec3(1)<0.or.ksec3(1)>65535)then! check bit-map table reference field.
       kret = 1
       write(*,9001) ksec3(1)
    endif
    
    if(ksec3(1)/=0.and.ksec1(2)==98)then!  ecmwf usage is to include the bit-map.
       kret = 1
       write(*,9002)
    endif
9001 format(1h ,'grchk3 : invalid bit-map table reference - ',i9)

9002 format(1h ,'grchk3 : ecmwf usage is to include bit-map.')

    return

  end subroutine grchk3
  subroutine grchk4(ksec4,kret)
    implicit none
    !     
    ! grchk4 - check parameters for section 4 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     check parameters for section 4 of grib code against
    !     valid values for grib edition 1.
    !     
    !   interface.
    !     ----------
    !     
    !     call grchk4(ksec4,kret)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec4      - array containing integer parameters for
    !     section 4 of grib code.
    !     
    !     
    !     output parameters.
    !     ------------------
    !     
    !     kret       - return code.

    !     0   , no error encountered.
    !     1   , error in grib code parameter.
    !     
    !     method.
    !     -------
    !     
    !     values checked against current code/flag tables
    !     and against maximum or minimum permitted values.
    !     they are also checked against the current status
    !     of the implementation of gribex.
    !     
    !     externals.
    !     ----------
    !     
    !     none.
    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes for grib code.
    !     
    !     comments.
    !     ---------
    !     
    !     routine contains sections 0 to 2 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     j. hennessy      ecmwf      18.06.91
    !     
    !     modifications.
    !     --------------
    !     
    integer:: kret

    integer:: ksec4(*)



    kret = 0
    !         section 2 . check values against code tables and extreme values.


    !         check number of values to be encoded.

    if(ksec4(1)==0)then
       kret = 1
       write(*,9001)
    endif
    !         check number of bits per data value.
    !     the following allows zero bits per datum(for constant fields)
    if(ksec4(2)<0.or.ksec4(2)>32)then
       kret = 1
       write(*,9002) ksec4(2)
    endif
    !         check on type of data(grid or spherical harmonics).

    if(ksec4(3)/=0.and.ksec4(3)/=128)then
       write(*,9003)
    endif

    !         check type of packing.

    if(ksec4(4)/=0)then
       kret = 1
       write(*,9004)
    endif
    !         check data representation.

    if(ksec4(5)/=0.and.ksec4(5)/=32)then
       kret = 1
       write(*,9005)
    endif

    !         additional flags not supported.

    if(ksec4(6)/=0.and.ksec4(7)/=0.and.ksec4(8)/=0.and. &
         ksec4(9)/=0.and.ksec4(10)/=0.and.ksec4(11)/=0)then
       kret = 1
       write(*,9006)
    endif
    !         section 9 . return to calling routine. format statements.

9001 format(1h ,'grchk4:number of data values to be encoded is 0.')
9002 format(1h ,'grchk4 : invalid number of bits for packed data', &
         ' values - ',i3)
9003 format(1h ,'grchk4 : warning - illegal type of data.')
9004 format(1h ,'grchk4 : only simple packing supported.')
9005 format(1h ,'grchk4 : illegal data representation.')
9006 format(1h ,'grchk4 : additional flags not supported.')

    return

  end subroutine grchk4



  subroutine gribex(ksec0,ksec1,ksec2,psec2,ksec3,psec3,ksec4,&
       field,fieldSize,kgrib,kleng,kword,hoper,kret)
    !     
    ! GRIBEX - Coding and decoding of GRIB format data.
    !     
    !     Purpose.
    !     --------

    !     1) Code data in FM-92 GRIB code, Edition 1.
    !     2) Decode data from FM-92 GRIB code.
    !     3) Decode only identification sections of GRIB
    !     coded data ie Sections 0, 1 and 2.
    !     4) Return length of GRIB message, in bytes, and GRIB
    !     Edition number only.
    !     
    !     A number of options exist when coding or decoding -
    !     see values allowed for requested function, HOPER, below.
    !     
    !     Decoding functions work on Experimental Edition,
    !     Edition 0 and Edition 1 of GRIB code. Decoded values
    !     for Sections 0 to 2 are always in Edition 1 format.
    !     
    !     
    !     Input Parameters for all functions.
    !     -----------------------------------
    !     
    !     HOPER      - Requested function.
    !     
    !     'C' To code data in GRIB code, with or
    !     without bit-maps.
    !     
    !     'D' To decode data from GRIB code. If
    !     ECMWF pseudo-Grib data is encountered,
    !     only sections 0 and 1 are decoded and
    !     the return code is set to -6.
    !     
    !     'I' To decode only identification
    !     sections 0, 1 and 2 of GRIB or
    !     pseudo-Grib data.
    !     
    !     'L' Return length of GRIB message, in
    !     bytes, and GRIB Edition number only.
    !     Length does not include any bytes
    !     added to round message length to a
    !     multiple of 120 bytes. Works also for
    !     pseudo-Grib data.
    !     
    !     'M' To code data in GRIB code and, if a
    !     bit-map is encountered, make GRIB
    !     message full length ie the same length
    !     as if all data values were given.
    !     
    !     'R' To decode data from GRIB code, and if
    !     a quasi-regular Gaussian grid is
    !     encountered, convert it to regular.
    !     
    !     'S' To decode initialised analysis data
    !     from GRIB code, and if data is in the
    !     Experimental Edition of GRIB, set the
    !     Time Range Indicator flag. In the
    !     Experimental Edition there was no
    !     distinction between initialised and
    !     uninitialised analyses.
    !     
    !     'X' To extract data values for up to 4
    !     points from a GRIB coded Gaussian or
    !     Latitude/longitude field, without
    !     unpacking the data at other points.
    !     
    !     'Z' To decode data from GRIB code.
    !     If a bit-map is encountered,
    !     only sections 0,1 and 2 are decoded and
    !     the return code is set to -5.
    !     
    !     FIELDSIZE      - Length of array FIELD.
    !     
    !     KLENG      - Length of array KGRIB.
    !     
    !     KRET       - Response to error indicator.
    !     0         , Abort if error encountered.
    !     Negative return codes are
    !     informative and do not cause
    !     an abort.
    !     Non- zero , Return to calling routine
    !     even if error encountered.
    !     
    !     Input parameters for coding function.
    !     Output Parameters for decoding functions.
    !     -----------------------------------------
    !     
    !     KSEC1      - Integer parameters of Section 1(Product
    !     Definition Section) of GRIB code.
    !     Integer array of at least 25 words.
    !     
    !     Word   Contents.
    !     ----   ---------
    !     1    Version number of Code Table 2.
    !     2    Identification of centre(Code Table 0).
    !     3    Generating process identification number
    !    ( allocated by originating centre ).
    !     4    Grid definition(NNN -  Catalogue number
    !     of grid used by originating centre. See
    !     Volume B of publication WMO - No.9).
    !     5    Flag indication relative to Section 2
    !    (Grid Description Section) and Section
    !     3(Bit Map Section). Code Table 1.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Sections 2 and 3 omitted.
    !     128      Section 2 included, Section 3
    !     omitted.
    !     64      Section 2 omitted, Section 3
    !     included.
    !     192      Sections 2 and 3 included.
    !     
    !     6    Indicator of parameter(Code Table 2).
    !     7    Indicator of type of level(Code Table 3).
    !     8    Height, pressure etc of level(Code Table 3).
    !     Single level or top of layer.
    !     9    Height, pressure etc of level(Code Table 3).
    !     Bottom of layer, if word 6 indicates a layer.
    !     10    Year of century  }
    !     11    Month            } Reference time of data -
    !     12    Day              } Date and time of start of
    !     13    Hour             } averaging or accumulation
    !     14    Minute           } period.
    !     15    Indicator of unit of time(Code Table 4).
    !     16    P1 - Period of time(number of time units)
    !    (0 for analyses or initialised analyses).
    !     17    P2 - Period of time(number of time units);
    !     or time interval between successive
    !     analyses, initialised analyses or forecasts
    !     undergoing averaging or accumulation;
    !     otherwise set to zero.
    !     18    Time range indicator(Code Table 5).
    !     19    Number included in average, when time range
    !     indicator indicates an average or
    !     accumulation; otherwise set to zero.
    !     20    Number missing from average, when time range
    !     indicator indicates an average or
    !     accumulation; otherwise set to zero.
    !     21    Century of reference time of data.
    !     22    Reserved. set to 0.
    !     23    Decimal scale factor.
    !     24-25    Reserved. Set to 0.
    !     
    !     KSEC2      - Integer parameters of Section 2(Grid
    !     Description Section) of GRIB code.
    !     Integer array of at least 22 + n words,
    !     where n is the number of parallels or
    !     meridians in a quasi-regular(reduced)
    !     Gaussian or latitude/longitude grid.
    !     
    !     Notes:- 1) Latitudes, longitudes are in
    !     millidegrees.
    !     2) Latitude values in the range 0-90000.
    !     3) Longitude values in the range 0-360000.
    !     4) Southern latitudes and western
    !     longitudes are negative.
    !     
    !     Word   Contents for latitude/longitude grids or
    !     equidistant cylindrical or Plate Carree.
    !     ----   ----------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Ni - Number of points along a parallel.
    !     3    Nj - Number of points along a meridian.
    !     4    La1 - Latitude of first grid point.
    !     5    Lo1 - Longitude of first grid point.
    !     6    Resolution flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Direction increments not given.
    !     Used for quasi-regular grids, but
    !     can also be used for regular grids.
    !     128      Direction increments given.
    !     Grids must be regular.
    !     
    !     7    La2 - Latitude of last grid point.
    !     8    Lo2 - Longitude of last grid point.
    !     9    Di - i direction increment.
    !     10    Dj - j direction increment.
    !     11    Scanning mode flags(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    Latitude of the southern pole of rotation.
    !     14    Longitude of the southern pole of rotation.
    !     15    Latitude of the the pole of stretching.
    !     16    Longitude of the the pole of stretching.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     At the moment quasi-regular latitude/
    !     longitude grids are not properly defined.
    !     The Resolution flag field indicates both
    !     direction increments are given or not.
    !     One increment needs to be given. Grids
    !     can be irregular in one direction only.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     
    !     20-22    Reserved. Set to 0.
    !     23-nn    Number of points along each parallel
    !     in a Quasi-regular grid. Number of parallels
    !     is given by Nj above.
    !     or
    !     Number of points along each meridian
    !     in a Quasi-regular grid. Number of  meridians
    !     is given by Ni above.
    !     
    !     Scanning mode flags(Code Table 8) indicate
    !     whether points are consecutive on a meridian
    !     or a parallel.
    !     
    !     Notes:- 1) Increments are in millidegrees.
    !     
    !     Word   Contents for Gaussian grids .
    !     ----   ---------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Ni - Number of points along a parallel.
    !     Cannot be used for quasi-regular grids.
    !     3    Nj - Number of points along a meridian.
    !     4    La1 - Latitude of first grid point.
    !     5    Lo1 - Longitude of first grid point.
    !     6    Resolution flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Direction increments not given.
    !     Used for quasi-regular grids, but
    !     can also be used for regular grids.
    !     128      Direction increments given.
    !     Grids must be regular.
    !     
    !     7    La2 - Latitude of last grid point.
    !     8    Lo2 - Longitude of last grid point.
    !     9    Di - i direction increment.
    !     Cannot be used for quasi-regular grids.
    !     10    N - Number of parallels between a Pole and
    !     the Equator.
    !     11    Scanning mode flags(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    Latitude of the southern pole of rotation.
    !     14    Longitude of the southern pole of rotation.
    !     15    Latitude of the the pole of stretching.
    !     16    Longitude of the the pole of stretching.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     
    !     20-22    Reserved. Set to 0.
    !     23-nn    Number of points along each parallel
    !     in a Quasi-regular grid. Number of parallels
    !     is given by Nj above.
    !     
    !     Notes:- 1) Increments are in millidegrees.
    !     
    !     Word   Contents for Spherical Harmonic Coefficients.
    !     ----   --------------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    J - Pentagonal resolution parameter.
    !     3    K - Pentagonal resolution parameter.
    !     4    M - Pentagonal resolution parameter.
    !     5    Representation type( Code Table 9 ).
    !     6    Representation mode( Code Table 10 ).
    !     7-11    Reserved. Set to 0.
    !     12    Number of vertical coordinate parameters.
    !     13    Latitude of the southern pole of rotation.
    !     14    Longitude of the southern pole of rotation.
    !     15    Latitude of the the pole of stretching.
    !     16    Longitude of the the pole of stretching.
    !     17-22    Reserved. Set to 0.
    !     
    !     Word   Contents for Polar Stereographic.
    !     ----   --------------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Nx - Number of points along X-axis.
    !     3    Ny - Number of points along Y-axis.
    !     4    La1 - Latitude of first grid point.
    !     5    Lo1 - Longitude of first grid point.
    !     6    Reserved. Set to 0. Resolution flag is
    !     not applicable to Polar stereographic.
    !     7    LoV - Orientation of the grid ie the
    !     longitude of the meridian which is parallel
    !     to the Y-axis along which latitude increases
    !     as the Y-coordinate increases.
    !     8    Reserved. Set to 0.
    !     9    Dx - X-direction grid length.
    !     10    Dy - Y-direction grid length.
    !     11    Scanning mode flag(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    Projection centre flag.
    !     0 , North pole is on projection plane.
    !     1 , South pole is on projection plane. ??????
    !     128 , South pole is on projection plane. ????
    !     14-16    Reserved. Set to 0.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     20-22    Reserved. Set to 0.
    !     
    !     Notes   1) Grid lengths are in metres, at the 60-
    !     degree parallel nearest to the pole on
    !     the projection plane.
    !     
    !     Word   Contents for Mercator.
    !     ----   ---------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Ni - Number of points along a parallel.
    !     3    Nj - Number of points along a meridian.
    !     4    La1 - Latitude of first grid point.
    !     5    Lo1 - Longitude of first grid point.
    !     6    Resolution flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Direction increments not given.
    !     128      Direction increments given.
    !     
    !     7    La2 - Latitude of last grid point.
    !     8    Lo2 - Longitude of last grid point.
    !     9    Latin - latitude at which the Mercator
    !     projection cylinder intersects the earth.
    !     10    Reserved. set to 0.
    !     11    Scanning mode flags(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    Di - i direction grid length.
    !     14    Dj - j direction grid length.
    !     15-16    Reserved. Set to 0.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     
    !     20-22    Reserved. Set to 0.
    !     
    !     Notes   1) Grid lengths are in units of metres,
    !     at the parallel specified by Latin.
    !     
    !     Word   Contents for Lambert conformal, secant or
    !     tangent, conical or bi-polar(normal or
    !     oblique) or
    !     Albers equal-area, secant or tangent,
    !     conical or bi-polar(normal or oblique).
    !     ----   --------------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Nx - Number of points along X-axis.
    !     3    Ny - Number of points along Y-axis.
    !     4    La1 - Latitude of first grid point.
    !     5    Lo1 - Longitude of first grid point.
    !     6    Resolution flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Direction increments not given.
    !     128      Direction increments given.
    !     
    !     7    LoV - Orientation of the grid ie the  East
    !     longitude of the meridian which is parallel
    !     to the Y-axis along which latitude increases
    !     as the Y-coordinate increases.
    !     8    Reserved. Set to 0.
    !     9    Dx - X-direction grid length.
    !     10    Dy - Y-direction grid length.
    !     11    Scanning mode flag(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    Projection centre flag.
    !     0 , North pole is on projection plane.
    !     Only one projection centre is used.
    !     128 , South pole is on projection plane.
    !     Only one projection centre is used.
    !     64 , North pole is on projection plane.
    !     Projection is bi-polar and symmetric.
    !     192 , South pole is on projection plane.
    !     Projection is bi-polar and symmetric.
    !     14    Latin 1 - First latitude from the pole at
    !     which the secant cone cuts the sphere.
    !     15    Latin 2 - Second latitude from the pole at
    !     which the secant cone cuts the sphere.
    !     16    Reserved. Set to 0.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     
    !     20    Latitude of the southern pole.
    !     21    Longitude of the southern pole.
    !     22    Reserved. Set to 0.
    !     
    !     Notes   1) Grid lengths are in metres, at the 60-
    !     degree parallel nearest to the pole on
    !     the projection plane.
    !     
    !     Word   Contents for Space view perspective
    !     or orthographic.
    !     ----   ---------------------------------------
    !     1    Data representation type(Code Table 6).
    !     2    Nx - Number of points along x-axis.
    !     3    Ny - Number of points along y-axis.
    !     4    Lap - Latitude of sub-satellite point.
    !     5    Lop - Longitude of sub-satellite point.
    !     6    Resolution flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Direction increments not given.
    !     128      Direction increments given.
    !     
    !     7    dx - Apparent diameter of the earth in
    !     grid lengths in the x direction.
    !     8    dy - Apparent diameter of the earth in
    !     grid lengths in the y direction.
    !     9    Xp X-coordinate of sub-satellite point
    !     10    Yp Y-coordinate of sub-satellite point
    !     11    Scanning mode flag(Code Table 8).
    !     12    Number of vertical coordinate parameters.
    !     13    The orientation of the grid.
    !     14    nr - the altitude of the camera from the
    !     earth's centre.
    !     15-16    Reserved. Set to 0.
    !     17    0 , Regular grid.
    !     1 , Quasi-regular(reduced) grid.
    !     18    Earth flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Earth assumed spherical with
    !     radius of 6367.47 km.
    !     64      Earth assumed oblate spheroidal
    !     with size as determined by IAU in
    !     1965 :
    !    (6378.160km,6356.775km,f=1/297.0)
    !     
    !     19    Components flag.
    !     Valid values are :-
    !     
    !     Decimal
    !     value     Meaning
    !     -----     -------
    !     0      Resolved u and v components of
    !     vector quantities relative to
    !     easterly and northerly directions.
    !     
    !     8      Resolved u and v components of
    !     vector quantities relative to the
    !     defined grid in the direction of
    !     increasing x and y(or i and j)
    !     coordinates respectively.
    !     
    !     20-22    Reserved. Set to 0.
    !     
    !     PSEC2      - Real parameters for Section 2(Grid
    !     Definition Section) of GRIB Code.
    !     Real array of at least 10 + nn words, where
    !     nn is the number of vertical coordinate
    !     parameters.
    !     
    !     Word   Contents.
    !     ----   --------------------------------------------
    !     1      Angle of rotation.
    !     2      Stretching factor.
    !     3-10     Reserved. Set to 0.
    !     11-nn     Vertical coordinate parameters.
    !     Number given in KSEC2(12)
    !     
    !     KSEC3      - Integer parameters for Section 3(Bit Map
    !     Section) of GRIB code.
    !     Integer array of at least 2 words.
    !     
    !     Word   Contents.
    !     ----   --------------------------------------------
    !     1      0 , Bit map included in the GRIB message.
    !     Binary data array(FIELD) contains the
    !     missing data indicator at the points
    !     where no data is given.
    !     Non-zero, Number of predetermined bit-map.
    !     Bit map is not included in the message.
    !     Binary data array contains only valid
    !     data values.
    !     
    !     2      The value used to indicate missing data in
    !     an integer binary data array is indicated
    !     here.
    !     This value is user supplied for both
    !     coding and decoding.
    !     
    !     PSEC3      - Real parameters for Section 3(Bit Map
    !     Section) of GRIB code.
    !     Real array of at least 2 words.
    !     
    !     Word   Contents.
    !     ----   --------------------------------------------
    !     1      Not used.
    !     
    !     2      The value used to indicate missing data in
    !     a real binary data array is indicated here.
    !     This value is user supplied for both
    !     coding and decoding.
    !     
    !     KSEC4      - Integer parameters for Section 4(Binary
    !     Data Section) of GRIB code.
    !     Integer array of at least 42 words.
    !     
    !     Word   Contents.
    !     ----   --------------------------------------------
    !     1    Number of data values in array FIELD to be
    !     packed in GRIB code or which have been
    !     unpacked from GRIB code.
    !     
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !                                              !
    !     ! If this number is NEGATIVE, it indicates     !
    !     ! ENTIRE FIELD IS MISSING. All values in FIELD !
    !     ! are 0. This is an ECMWF convention - coded   !
    !     ! data has scale factor, exponent and mantissa !
    !     ! of reference value with all bits set to 1.   !
    !     !                                              !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     
    !     2    Number of bits used for each packed value.
    !     3    Type of data. Used only if Section 2 is
    !     not included when coding data.
    !     0 - Grid point data.
    !     128 - Spherical harmonic coefficients.
    !     4    Type of packing.
    !     0 - Simple packing.
    !     64 - Complex or second order packing.
    !    (Not implemented)
    !     5    Data representation.
    !     0 - Floating point data.
    !     32 - Integer data.
    !     6    Additional flags indicator.
    !     0 - No additional flags.
    !     16 - additional flags.
    !    (Not implemented)
    !     7    Reserved. Set to 0.
    !     8    Number of values indicator.
    !     0 - Single datum at each grid point.
    !     64 - Matrix of values at each grid point.
    !    (Not implemented)
    !     9    Secondary bit maps indicator.
    !     0 - No secondary bit maps.
    !     32 - Secondary bit maps present.
    !    (Not implemented)
    !     10    Values width indicator.
    !     0 - Second order values constant width.
    !     16 - Second order values different widths.
    !    (Not implemented)
    !     11    Number of bits for second order values,
    !     when of constant width.
    !    (Not implemented)
    !     12-32    Reserved. Set to 0.
    !     33    Scale factor used for encoding
    !     
    !     Words 34 to 42 are used only for the 'X'
    !     function. Scanning mode must be from West
    !     to East and from North to South.
    !     34    Number of points(maximum 4) from which
    !     data is to be unpacked.
    !     35    Number of latitude row of first value.
    !     36    Number of longitude point of first value.
    !     37    Number of latitude row of second value.
    !     38    Number of longitude point of second value.
    !     39    Number of latitude row of third value.
    !     40    Number of longitude point of third value.
    !     41    Number of latitude row of fourth value.
    !     42    Number of longitude point of fourth value.
    !     
    !     FIELD      - Array of data values to be packed in GRIB
    !     code or which have been unpacked. Where a
    !     bit-map is included in the GRIB message
    !     this array contains missing data indicator
    !    ( value supplied by the user in PSEC3(2)
    !     or KSEC3(2) ) at the appropriate places.
    !     
    !     Although declared as real in GRIBEX this
    !     can be an array of integer data. The value
    !     in KSEC4(5) indicates whether data is in
    !     integer or floating point format.
    !     
    !     Output parameters for coding function.
    !     Input Parameters for decoding functions.
    !     -----------------------------------------
    !     
    !     KGRIB      - Array containing GRIB coded data.
    !     
    !     KWORD      - Number of words of KGRIB occupied by
    !     coded data. Output parameter for coding
    !     function only. Not required as input
    !     for decoding.
    !     
    !     Output Parameters for all functions.
    !     -----------------------------------
    !     
    !     KSEC0      - Word 1 contains number of octets in
    !     GRIB message(not including padding to
    !     a word boundary or rounding to a multiple
    !     of 120 octets).
    !     - Word 2 contains GRIB edition number.
    !     
    !     KRET       - Return code.
    !     
    !     Informative codes for decoding functions.
    !     
    !     -2  , Bit-map encountered with all bits
    !     set to 1. Array FIELD contains all
    !     real data values.
    !     -3  , Predetermined bit-map encountered.
    !     Data has not been fully decoded ie
    !     array FIELD contains only real data
    !     values. The user must use this data
    !     in conjunction with the defined
    !     bit-map.
    !     -4  , Bit-map encountered. The data has
    !     been fully decoded ie array FIELD
    !     contains real values and missing
    !     data indicators where appropriate.
    !     -5  , Bit-map encountered. The data has
    !     not been decoded. This return code
    !     is set only by the 'Z' function.
    !     -6  , ECMWF pseudo-grib data encountered.
    !     
    !     Error codes.
    !     
    !     0   , No error encountered.
    !     201 , Invalid function requested.
    !     202 , Number of bits per data value exceeds
    !     word length.
    !     203 , Missing data indicated and data field
    !     contains non-zero values.
    !     
    !     301 , Error in inserting letters GRIB.
    !     302 , Error extracting length of GRIB
    !     message.
    !     303 , Error inserting/extracting GRIB
    !     Edition Number.
    !     304 , Error extracting octets 22 and 23
    !     Experimental Edition check.
    !     401 , Error inserting/extracting length of
    !     Section 1.
    !     402 , Error inserting/extracting Parameter
    !     Version Number.
    !     403 , Error inserting/extracting six fields
    !     from Identification of Centre to
    !     Indicator of type of level.
    !     404 , Error inserting/extracting Height,
    !     pressure, etc of levels.
    !     405 , Error inserting/extracting six fields
    !     from Year of century to Indicator
    !     of unit of time range.
    !     406 , Error inserting/extracting Period of
    !     time.
    !     407 , Error inserting/extracting time range
    !     indicator.
    !     408 , Error inserting/extracting number
    !     averaged.
    !     409 , Error inserting/extracting number
    !     missing from averages etc.
    !     410 , Error inserting/extracting century of
    !     data or reserved field.
    !     411 , Error inserting/extracting units
    !     decimal scale factor.
    !     499 , Error found when checking values for
    !     Section 1 against valid GRIB values.
    !     
    !     501 , Error inserting/extracting length of
    !     Section 2.
    !     502 , Error inserting/extracting number of
    !     Vertical coordinate parameters.
    !     503 , Error inserting/extracting location
    !     of List of vertical coordinate
    !     parameters or List of numbers of
    !     points.
    !     504 , Error inserting/extracting data
    !     representation type.
    !     505 , Error inserting/extracting number of
    !     points along a parallel or meridian.
    !     506 , Error inserting/extracting latitude
    !     or longitude of first grid point.
    !     507 , Error inserting/extracting components
    !     flag.
    !     508 , Error inserting/extracting latitude
    !     or longitude of last grid point.
    !     509 , Error inserting/extracting i
    !     direction increment.
    !     510 , Error inserting/extracting number of
    !     parallels between pole and Equator.
    !     511 , Error inserting/extracting scanning
    !     mode flags.
    !     513 , Error inserting/extracting j
    !     direction increment.
    !     514 , Error inserting/extracting J,K,M
    !     pentagonal resolution parameters.
    !     515 , Error inserting/extracting
    !     representation type or mode.
    !     517 , Error inserting/extracting latitude
    !     or longitude of southern pole.
    !     518 , Error inserting/extracting angle
    !     of rotation.
    !     519 , Error inserting/extracting latitude
    !     or of pole of stretching.
    !     520 , Error inserting/extracting
    !     stretching factor.
    !     521 , Error inserting/extracting
    !     vertical coordinate parameters.
    !     522 , Error inserting/extracting list of
    !     numbers of points.
    !     523 , Error inserting/extracting number of
    !     points along X or Y axis.
    !     524 , Error inserting/extracting X or Y
    !     axis grid lengths.
    !     525 , Error inserting/extracting Projection
    !     centre flag.
    !     526 ,  Error inserting/extracting latitude
    !     or  longitude of sub-satellite point.
    !     527 , Error inserting/extracting diameter
    !     of the earth in x or y direction.
    !     528 , Error inserting/extracting X or Y
    !     coordinate of sub-satellite point.
    !     529 , Error inserting/extracting orientatio
    !     of the grid or camera angle.
    !     598 , Representation type not catered for.
    !     599 , Error found when checking values for
    !     Section 2 against valid GRIB values.
    !     601 , Error inserting/extracting length of
    !     Section 3.
    !     602 , Error inserting/extracting number of
    !     unused bits at end of section 3.
    !     603 , Error inserting/extracting bit-map
    !     reference table.
    !     604 , Error inserting/extracting bit-map.
    !     605 , Cannot convert Quasi-regular
    !     Gaussian grid with a bit-map.
    !     699 , Error found when checking values for
    !     Section 3 against valid GRIB values.
    !     701 , Error inserting/extracting length of
    !     Section 4.
    !     705 , Only simple packing catered for.
    !     706 , Error in extracting section 4 flag
    !     field.
    !     707 , Error inserting/extracting scale
    !     factor.
    !     708 , Error inserting/extracting reference
    !     value.
    !     709 , Error inserting/extracting number of
    !     bits per data value.
    !     710 , Output array too small.
    !     711 , Error inserting/extracting real
    !     coefficient.
    !     712 , Error inserting/extracting data
    !     values.
    !     713 , Error inserting/extracting flag
    !     and unused bit field.
    !     714 , Function is 'X' and number of
    !     values is illegal.
    !     715 , Function is 'X' and scanning mode is
    !     not North to South and West to East.
    !     716 , Function is 'X' and field is not
    !     Gaussian or Latitude/longitude grid.
    !     717 , Function is 'X' and a bit-map is
    !     included.
    !     799 , Error found when checking values for
    !     Section 4 against valid GRIB values.
    !     801 , Error inserting/extracting 7777 group
    !     802 , Error inserting/extracting length of
    !     GRIB message.
    !     805 , End of message 7777 group not found.
    !     
    !     Method.
    !     -------
    !     
    !     Input data packed in GRIB code in accordance with
    !     parameters given or set by user or fully or partially
    !     unpacked, depending on function used.
    !     
    !     Externals.
    !     ----------
    !     
    !     CONFP3
    !     DECFP2
    !     EXSCAL
    !     GRCHK1
    !     GRCHK2
    !     GRCHK3
    !     GRCHK4
    !     GRPRS0
    !     GRPRS1
    !     GRPRS2
    !     GRPRS3
    !     GRPRS4
    !     INSCAL
    !     INXBIT
    !     INXMAP
    !     MAXMIN
    !     QU2REG
    !     RORINT
    !     SETPAR
    !     
    !     Reference.
    !     ----------
    !     
    !     WMO Manual on Codes for GRIB definition.
    !     WMO Publication No. 9, Volume B, for grid catalogue numbers.
    !     
    !     Comments.
    !     ---------
    !     
    !     All machine dependent code is in 3 low level routines.
    !     Versions of these exist for the VAX, CYBER, IBM
    !     and SUN workstation, as well as the CRAY.
    !     INXBIT - contains calls to the routines GBYTE(S)
    !     and SBYTE(S) or their equivalents.
    !     SETPAR - to set number of bits in the computer word
    !     and largest negative number.
    !     
    !     This is not a full implementation of GRIB Code
    !     Edition 1. This routine codes/decodes regular
    !     latitude/longitude grids, regular and quasi-regular
    !     Gaussian grids, spherical harmonic , Space view
    !     perspective or orthographic and Polar
    !     Stereographic data. Grids may be rotated, stretched
    !     or rotated and stretched, with or without bit-maps.
    !     Only simple packing of data(real or integer) with
    !     no additional flag field is allowed.
    !     
    !     Apart from the values which can be passed to  this
    !     routine, other values are held in a common area and
    !     are used by default, unless changed by calls to the
    !     appropriate routines before calling GRIBEX. The
    !     following defaults are used. They have been selected
    !     to facilitate the most frequent usage at ECMWF and
    !     to ease the transition to the next version of the
    !     GRIB code.
    !     
    !     1) By default debug printout is switched off.
    !     CALL GRSDBG(I) where
    !     I = Non-zero to switch on debug printout.
    !     0, to switch off debug printout.
    !     
    !     2) By default the reference value used is the minimum
    !     of the data values supplied.
    !     CALL GRSREF(ZREF) to change the value, where ZREF
    !     is real and the required value.
    !     
    !     3) By default GRIB messages are rounded to a
    !     multiple of 120 octets.
    !     CALL GRSRND(I) where
    !     I = 0, to switch off rounding.
    !     Non-zero to switch on rounding.
    !     
    !     4) By default, the values given  are checked for
    !     consistency with GRIB code values as currently
    !     defined, when coding data. Data values are never
    !     checked.
    !     CALL GRSVCK(I) where
    !     I = 0, to switch off checking.
    !     Non-zero to switch on checking.



    implicit none

    character(len=*)::hoper

    character(len=1)::yfunc
    character(len=1)::ytemp

    integer,allocatable,dimension(:)::mytmpv
    integer::mytmp(1)

    integer:: i
    integer:: iblen
    integer:: ibits
    integer:: ibmap
    integer:: icount
    integer:: iexp(1)
    integer:: iflag(1)
    integer:: igrib(4)
    integer:: il
    integer:: ilalo(2)
    integer:: ilen
    integer:: ilenf
    integer:: ilenfm
    integer:: ilen1(1)
    integer:: ilen2(1)
    integer:: ilen3(1)
    integer:: ilen4(1)
    integer:: imant(1)
    integer:: imisng
    integer:: imiss
    integer:: imoday(1)
    integer:: imsval
    integer:: inc
    integer:: inil
    integer:: inital
    integer:: inolat
    integer:: inolng
    integer:: inspt
    integer:: inub(1)
    integer:: inum
    integer:: ioff
    integer:: iparm(4)
    integer:: ipl
    integer:: iplen
    integer:: ipseud
    integer:: ipvpl(1)
    integer:: irep
    integer:: iresol(1)
    integer:: iret
    integer:: isbmap
    integer:: iscale(1)
    integer:: isign
    integer:: iskale
    integer:: iskip
    integer:: itemp(1)
    integer:: itrnd
    integer:: ivals
    integer:: i7777(4)

    integer,parameter:: jpedno=1
    integer,parameter:: jplen1=28

    integer:: j202
    integer:: j204
    integer:: j205
    integer:: j206
    integer:: j207
    integer:: j501
    integer:: j530
    integer:: j711
    integer:: j712
    integer:: j713

    integer:: j716
    integer:: j802
    integer:: kleng

    integer:: kgrib(kleng)
    integer,intent(in):: fieldsize
    integer:: kret
    integer:: ksec0(2)
    integer:: ksec1(25)
    integer:: ksec2(50)
    integer:: ksec3(2)
    integer:: ksec4(42)
    integer:: kword
    integer:: ist,ihe
    real(kind=realgribkind)::    psec2(400)
    real(kind=realgribkind)::    psec3(2)
    real(kind=realgribkind)::    field(fieldsize)

    real(kind=realgribkind)::    zmax
    real(kind=realgribkind)::    zmin
    !      real(kind=realgribkind)::    zmisng
    real(kind=realgribkind)::    zmsval
    real(kind=realgribkind)::    zref
    real(kind=realgribkind)::    zs
    real(kind=realgribkind)::    zscale
    real(kind=realgribkind)::    ztemp
    integer,parameter::jps4=10240
    integer:: itmpd(jps4)


    save igrib
    save i7777
    save inital
    save ibits
    save imisng

    !     characters grib and 7777 in ascii for use in sections 0 and 5
    !     of grib code.

    data igrib /71,82,73,66/
    data i7777 /55,55,55,55/

    !     first call flag.

    data inital /0/

    !         Set number of bits per computer word and missing data
    !     indicator, if first time through.

    if(inital==0)then
       call setpar(ibits,imisng,ndbg)
       inital = 1

       !     initialise mode of operatiom if not done before
       call grbmod(0,.false.)
    endif

    !         Set default values for parameters in common area, if values
    !     not already set , either on previous call or by user via calls
    !     to the GRS--- routines.

    if(nuser/=11041967)then

       !     Set all default values.

       !     User supplied reference value.

       fref   = 0._realgribkind
       !     reference value supplied by user flag. set to off.
       nfref  = 0
       !     set rounding to 120 bytes on.
       nrnd   = 1
       !     set debug print off.
       ndbg   = 0
       !     set grib value checking on.
       nvck   = 1
       !     mark common area values set.
       nuser  = 11041967
    endif
         
    !         when coding, print input parameters, if required.
         
    if(ndbg==1)then
       write(*,*) ' gribex : section 1.'
       if(hoper=='C'.or.hoper=='M'.or. & 
            hoper=='c'.or.hoper=='m' )then
          write(*,*) '          input values used -'
          ksec0(2) = jpedno
          call grprs1(ksec0,ksec1)
          call grprs2(ksec0,ksec2,psec2)
          call grprs3(ksec3,psec3)
          call grprs4(ksec4,field)
       endif
    endif

    !     Reset return code to 0, retaining input value to decide
    !     on abort / no abort, if error encountered later.
         
    iret   = kret
    kret   = 0
         
    !     IPSEUD is used to indicate pseudo-GRIB data encountered,
    !     when decoding.
    !     ISBMAP is the bit-map section flag and indicates what decoding
    !     has been done on bit-maps and data.
    !     See informative return codes for KRET when decoding.
         
    ipseud = 0
    isbmap = 0
         
    !     Reset bit-pointer to 0.
         
    inspt = 0
         
    ! Check input parameters.
    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 2.'

    !         Check that valid function has been requested.
         
    yfunc = hoper
    IF(YFUNC/='C'.AND.YFUNC/='D'.AND.YFUNC/='I'.AND.  &
         YFUNC/='L'.AND.YFUNC/='R'.AND.YFUNC/='S'.AND.  &
         YFUNC/='X'.AND.YFUNC/='Z'.AND.YFUNC/='M'.and.  &
         yfunc/='c'.and.yfunc/='d'.and.yfunc/='i'.and.  &
         yfunc/='l'.and.yfunc/='r'.and.yfunc/='s'.and.  &
         yfunc/='x'.and.yfunc/='z'.and.yfunc/='m') THEN
       kret = 201
       write(*,9201) hoper , kret
       go to 900
    endif

    !         function 'l' returns the length of the grib message and
    !     grib edition number only, so no array initialisation is
    !     necessary.
    !     
    IF(YFUNC=='L'.or.yfunc=='l') go to 300
    !     
    !         Function 'M' is for coding data, and if a bit map is encountered
    !     GRIB messages are made a fixed length. HOPER is passed to the
    !     bit-map handling routine.
    !     
    IF(HOPER=='M'.or.hoper=='m') yfunc = 'C'
    !     
    !         Function 'I' is for decoding of sections 0, 1 and 2
    !     of GRIB code only. Value of HOPER is checked at the start of
    !     decoding section 3.
    !     
    IF(HOPER=='I'.or.hoper=='i') yfunc = 'D'
    !     
    !         Function 'R' is the same as 'D', but if a quasi-regular
    !     Gaussian is encountered, it is converted to a regular one.
    !     Value of HOPER is checked near end of section 8.
    !     
    IF(HOPER=='R'.or.hoper=='r') yfunc = 'D'
    !     
    !         Function 'S' is the same as 'D', but if analysis data in
    !     GRIB Experimental Edition is encountered, the time range
    !     indicator flag is set to indicate initialised analysis.
    !     Value of HOPER is checked  when time range indicator has
    !     been extracted.
    !     
    IF(HOPER=='S'.or.hoper=='s') yfunc = 'D'
    !     
    !         Function 'X' is the same as 'D', but only the data
    !     at the requested points is unpacked.
    !     Value of HOPER is checked prior to unpacking data values.
    !     
    IF(HOPER=='X'.or.hoper=='x') yfunc = 'D'

    !         Function 'Z' is for decoding only the information which could
    !     be handled by the the old decoding routine DECOGB(eg no bit
    !     maps) and is used by the new DECOGB interface.
    !     
    IF(HOPER=='Z'.or.hoper=='z') yfunc = 'D'
    !     
    !         Preset some arrays to 0.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')then
       !     Set array to receive coded GRIB data to 0.
            
       do 202 j202 = 1,kleng
          kgrib(j202) = 0
202    enddo

       !                Check number of bits per data field.

       if(ibits<ksec4(2))then
          kret = 202
          write(*,9202) ksec4(2) , ibits , kret
          go to 900
       endif

       !                If missing field is indicated, check that data contains
       !     all zero values.
       !     
       if(ksec4(1)<0)then
          imiss = 1
          ilenf = - ksec4(1)
          do 204 j204 = 1 , ilenf
             if(abs(field(j204))>1.e-14_realgribkind)then
                kret = 203
                write(*,9203) kret
                go to 900
             endif
204       enddo
       else
          imiss = 0
          ilenf = ksec4(1)
       endif
    else

       !     Preset arrays to receive section header information to 0.
       !     Routine GSBITE resets data array to 0.
       !     
       do 205 j205=1,25
          ksec1(j205) = 0
205    enddo
       !     
       do 206 j206=1,22
          ksec2(j206) = 0
206    enddo
       !     
       do 207 j207=1,12
          ksec4(j207) = 0
207    enddo
       !     
    endif
    !     
    !     ----------------------------------------------------------------

    !         Section 3 . Indicator Section(Section 0) of GRIB code.
    !     ----------------------------------------------------------------
    !     
300 CONTINUE
    !     
    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 3.'
    !     
    !         Octets 1 - 4 : The letters G R I B.
    !     Four 8 bit fields.
    !     
    IF(YFUNC=='C'.or.YFUNC=='c')THEN
       !     
       !     Insert fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IGRIB(1),4,IBITS, &
            8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 301
          WRITE(*,9301) KRET
          GO TO 900
       ENDIF
    ELSE

       !     When decoding data, skip field. No check made on
       !     whether letters are actually GRIB.
       !     Update bit-pointer.
       !     
       INSPT = INSPT + 32
    ENDIF
    !     
    !         Octets 5 - 7 : Length of message.
    !     One 24 bit field.
    !     
    IF(YFUNC=='C'.or.YFUNC=='c')THEN
       !     
       !     When coding data, skip field. Length is inserted
       !     later, when known.
       !     Update bit-pointer.
       !     
       INSPT = INSPT + 24
    ELSE

       !     Extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC0(1),1,IBITS, &
            24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 302
          WRITE(*,9302) KRET
          GO TO 900
       ENDIF
    ENDIF

    !         Octet 8 : GRIB Edition Number.
    !     One 8 bit field.
    !     
    IF(YFUNC=='C'.or.YFUNC=='c')THEN
       !     
       !     Set value, if coding data.
       !     
       KSEC0(2) = JPEDNO
       !     
       !     Insert field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC0(2),1,IBITS, &
            8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 303
          WRITE(*,9303) KRET
          GO TO 900
       ENDIF

    ELSE
       !     Extract field.

       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC0(2),1,IBITS, &
            8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 303
          WRITE(*,9303) KRET
          GO TO 900
       ENDIF

       !                When decoding or calculating length, previous editions
       !     of the GRIB code must be taken into account.
       !     
       !     In the table below, covering sections 0 and 1 of the GRIB
       !     code, octet numbering is from the beginning of the GRIB
       !     message;
       !      indicates that the value is not available in the code
       !     edition;
       !     R indicates reserved, should be set to 0;
       !     Experimental edition is considered as edition -1.
       !     
       !     GRIB code edition -1 has fixed length of 20 octets for
       !     section 1, the length not included in the message.
       !     GRIB code edition 0 has fixed length of 24 octets for
       !     section 1, the length being included in the message.
       !     GRIB code edition 1 can have different lengths for section
       !     1, the minimum being 28 octets, length being included in
       !     the message.
       !     
       !     Octet numbers for code
       !     editions
       !     
       !     Contents.                   -1      0      1
       !     ---------                ----------------------
       !     Letters GRIB                          1-4    1-4    1-4
       !     Total length of GRIB message.                     5-7
       !     GRIB code edition number                           8
       !     Length of Section 1.                        5-7    9-11
       !     Reserved octet(R).                          8(R)   
       !     Version no. of Code Table 2.                      12
       !     Identification of centre.              5      9     13
       !     Generating process.                    6     10     14
       !     Grid definition .                      7     11     15
       !     Flag(Code Table 1).                   8     12     16
       !     Indicator of parameter.                9     13     17
       !     Indicator of type of level.           10     14     18
       !     Height, pressure etc of levels.      11-12  15-16  19-20
       !     Year of century.                      13     17     21
       !     Month.                                14     18     22
       !     Day.                                  15     19     23
       !     Hour.                                 16     20     24
       !     Minute.                               17     21     25
       !     Indicator of unit of time.            18     22     26
       !     P1 - Period of time.                  19     23     27
       !     P2 - Period of time                  20(R)   24     28
       !     or reserved octet(R).
       !     Time range indicator.                21(R)   25     29
       !     or reserved octet(R).
       !     Number included in average.       22-23(R)  26-27  30-31
       !     or reserved octet(R).
       !     Number missing from average.         24(R)  28(R)   32
       !     or reserved octet(R).
       !     Century of data.                                  33
       !     Reserved. set to 0.                               34
       !     Decimal scale factor.                            35-36
       !     Reserved. Set to 0.                              37-48
       !    (Need not be present)
       !     For originating centre use only.                 49-nn
       !    (Need not be present)
       !     
       !                Identify which GRIB code edition is being decoded.
       !     
       !     In GRIB edition 1, the edition number is in octet 8.
       !     In GRIB edition 0, octet 8 is reserved and set to 0.
       !     In GRIB edition -1, octet 8 is a flag field and can have a
       !     a valid value of 0, 1, 2 or 3.
       !     
       !     However, GRIB edition number 0 has a fixed
       !     length of 24, included in the message, for section 1, so
       !     if the value extracted from octets 5-7 is 24 and that from
       !     octet 8 is 0, it is safe to assume edition 0 of the code.
       !     
       IF(KSEC0(1)==24.AND.KSEC0(2)==0) THEN
          !     
          !     Set bit-pointer back by 32 bits(4 octets).
          !     
          INSPT = INSPT - 32
          !     
          !     Set length of GRIB message to missing data
          !     value.
          !     
          KSEC0(1) = IMISNG
          !     
          GO TO 400
       ENDIF

       !     In GRIB Edition -1, octets 22 and 23 are reserved and set
       !     to 0. These octets in Edition 1 are the month and the day,
       !     and must be non-zero.
       !     
       ITEMP(1) = 168
       CALL INXBIT(KGRIB,KLENG,ITEMP(1),IMODAY,1,IBITS, &
            16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 304
          WRITE(*,9304) KRET
          GO TO 900
       ENDIF
       !     
       IF(IMODAY(1)==0)THEN

          !     Set bit-pointer back by 32 bits(4 octets).
          !     
          INSPT = INSPT - 32
          !     
          KSEC0(2) = -1
          !     
          !     Set length of GRIB message to missing data
          !     value.
          !     
          KSEC0(1) = IMISNG
          !     
          !     Set length of section 1 of GRIB code to 20 octets.
          !     
          ILEN1(1) = 20
          !     
          !     Skip next 4 octets, as they do not exist in
          !     the Experimental Edition of the code. ie
          !     length of Section and Table 2 Version Number.
          !     
          GO TO 401
       ENDIF
       !     
    ENDIF

    !         If Grib Edition 1 and only length is required, go to section 9.
    !     
    IF(YFUNC=='L'.or.yfunc=='l') GO TO 900
    !     
    !     ----------------------------------------------------------------
    !     
    !         Section 4 . Product Definition Section(Section 1) of GRIB code.
    !     -----------------------------------------------------------------
    !     
400 CONTINUE

    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 4.'
    !     
    !         Octets 1 - 3 : Length of Section.
    !     One 24 bit field.
    !     
    !     Set value, if coding data.
    !     
    IF(YFUNC=='C'.or.YFUNC=='c') ILEN1(1) = JPLEN1
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,ILEN1,1,IBITS, &
         24,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 401
       WRITE(*,9401) KRET
       GO TO 900
    ENDIF

    !         Octet 4  : Version Number of Table 2.
    !     One 8 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(1),1,IBITS, &
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 402
       WRITE(*,9402) KRET
       GO TO 900
    ENDIF
    !     
401 CONTINUE

    !         Print length of Section 1, if required.
    !     
    IF(NDBG==1) WRITE(*,9101) ILEN1
    !     
    !         Octet 5  : Identification of centre.
    !     Octet 6  : Generating process identification.
    !     Octet 7  : Grid definition.
    !     Octet 8  : Flag.
    !     Octet 9  : Indicator of parameter.
    !     Octet 10 : Indicator of type of level.
    !     Six 8 bit fields.
    !     
    !     Insert / extract fields.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(2),6,IBITS, &
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 403
       WRITE(*,9403) KRET
       GO TO 900
    ENDIF

    !         Unless coding, fix-up for Experimental Edition and Edition 0
    !     of GRIB code.
    !     
    IF(YFUNC/='C'.and.yfunc/='c')THEN
       !     
       !     In GRIB Experimental Edition and Edition 0
       !     the International Table Version Number in use was 0.
       !     ECMWF has always used its own local table. It is the same
       !     for Experimental Edition, Edition 0 and Edition 1 and is
       !     local table number 128.
       !     
       IF(KSEC0(2)<1)THEN
          IF(KSEC1(2)==98)THEN
             !     ECMWF data. Local table number.
             KSEC1(1) = 128
          ELSE
             !     
             !     International table number.
             !     
             !     KSEC1(1) = 0 is already preset.
          ENDIF
       ENDIF
       !     
       !                Fix-up for flag field, which was different in Experimental
       !     Edition.
       !     
       !     Experimental          Editions 0 and 1
       !     Edition
       !     
       !     Sections     Binary    Decimal     Binary    Decimal
       !     included     value     value       value     value
       !     
       !     none       00000000    0         00000000     0

       !     2         00000001    1         10000000   128
       !     3         00000010    2         01000000    64
       !     2 and 3      00000011    3         11000000   192
       !     
       IF(KSEC0(2)==-1)THEN
          IF(KSEC1(5)==1) KSEC1(5) = 128
          IF(KSEC1(5)==2) KSEC1(5) = 64
          IF(KSEC1(5)==3) KSEC1(5) = 192
       ENDIF
    ENDIF

    !         Once the flag field has been extracted, no further fields
    !     from section 1 of the GRIB code are required, when length
    !     of GRIB or pseudo-Grib message only is required.
    !     
    IF(YFUNC=='L'.or.YFUNC=='l')THEN
       !     
       !     Length of section 0 + section 1 is 28 octets(224 bits)
       !     for GRIB Edition 0 and 24 octets(192 bits) for
       !     Experimental edition.
       !     Set bit-pointer and jump
       !     to extraction of length of section 3.
       !     
       INSPT = 224
       IF(KSEC0(2)==-1) INSPT = 192
       GO TO 500
    ENDIF

    !         Octets 11 - 12 : Height, pressure etc of levels.
    !     One 16 bit field or two 8 bit fields.
    !     
    !     For certain levels, no description is necessary, and
    !     when decoding the fields are already set to 0.
    ! not for pseudo-bufr
    IF(YFUNC/='C'.and.YFUNC/='c'.AND.KSEC1(2)/=137 &
         .AND.(KSEC1(7)<8.OR.KSEC1(7)==102))then


       !     Update bit-pointer and skip extraction.
       !     KSEC1(8) and KSEC1(9) already set to 0.

       INSPT = INSPT + 16
       GO TO 402
    ENDIF

    !     Certain level types require that the description occupies
    !     both octets.
    !     
    IF(KSEC1(7)==100.OR.KSEC1(7)==103.OR. &
         KSEC1(7)==105.OR.KSEC1(7)==107.OR. &
         KSEC1(7)==109.OR.KSEC1(7)==111.OR. &
         KSEC1(7)==160 .OR.KSEC1(1)==137)THEN

       !     One 16 bit field.
       !     
       INUM     = 1
       IBLEN    = 16
       KSEC1(9) = 0
    ELSE
       !     
       !     Two 8 bit fields.
       !     
       INUM  = 2
       IBLEN = 8
    ENDIF

    !     Insert / extract fields.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(8),INUM,IBITS,IBLEN,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 404
       WRITE(*,9404) KRET
       GO TO 900
    ENDIF

    !         Fix-up for ECMWF upper-air data incorrectly coded in Experimental
    !     Edition.
    !     
    IF(KSEC0(2)==-1.AND.KSEC1(2)==98)THEN
       ITEMP(1) = INSPT - 16
       INUM  = 2
       IBLEN = 8
       CALL INXBIT(KGRIB,KLENG,ITEMP(1),KSEC1(8),INUM,IBITS,IBLEN,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 404
          WRITE(*,9404) KRET
          GO TO 900
       ENDIF
       KSEC1(8) = KSEC1(8) * 32 + KSEC1(9)
       KSEC1(9) = 0
    ENDIF

402 CONTINUE

    !         Octet 13 : Year of century.
    !     Octet 14 : Month.
    !     Octet 15 : Day.
    !     Octet 16 : Hour.
    !     Octet 17 : Minute.
    !     Octet 18 : Indicator of unit of time range..
    !     Six 8 bit fields.
    !     
    !     Insert / extract fields.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(10),6,IBITS,8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 405
       WRITE(*,9405) KRET
       GO TO 900
    ENDIF

    !         Fix-up for unit of time, which was different in Experimental
    !     Edition.
    !     
    !     Experimental          Editions 0 and 1
    !     Edition
    !     
    !     Meaning      Decimal               Decimal
    !     value                 value
    !     
    !     Minute       0 or 30                 0
    !     Hour         1 or 40                 1
    !     Day          2 or 50                 2
    !     Month        3 or 60                 3
    !     Year         4 or 70                 4
    !     Decade       5 or 80                 5
    !     Normal         6                     6
    !     Century      7 or 90                 7
    !     Second         -                   254
    !     
    IF(KSEC0(2)==-1)THEN
       IF(KSEC1(15)==90) KSEC1(15) = 7
       IF(KSEC1(15)>10) KSEC1(15) =(KSEC1(15) / 10) - 3
    ENDIF
    !     
    !         Octets 19 - 20 : Period of time.
    !     One 16 bit field or two 8 bit fields.
    !     
    IF((YFUNC=='C'.or.yfunc=='c').AND.KSEC1(18)==10)THEN
       !     One 16 bit field.
       INUM  = 1
       IBLEN = 16
    ELSE

       !     Two 8 bit fields.

       INUM  = 2
       IBLEN = 8
    ENDIF
    !     Insert / extract fields.

    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(16),INUM,IBITS, &
         IBLEN,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 406
       WRITE(*,9406) KRET
       GO TO 900
    ENDIF

    !         Octet 21 : Time range indicator.
    !     One 8 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(18),1,IBITS, 8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 407
       WRITE(*,9407) KRET
       GO TO 900
    ENDIF

    !         When decoding, period of time field and time range
    !     indicator may need modification.
    !     
    IF(YFUNC/='C'.and.yfunc/='c')THEN
       !     
       !     When decoding, length of period of time field is known
       !     only at this time. If a 16 bit field is indicated, put
       !     the two extracted 8-bit fields together.
       !     
       IF(KSEC1(18)==10)THEN
          !     
          !     One 16 bit field.
          !     
          KSEC1(16) = KSEC1(16) * 256 + KSEC1(17)
          KSEC1(17) = 0
       ENDIF

       !     If data is known to be initialised analysis and GRIB is
       !     Experimental Edition, set time range indicator flag.
       !     
       IF(KSEC1(16)==0.AND.(HOPER=='S'.or.hoper=='s') &
            .AND.KSEC0(2)==-1)  KSEC1(18) = 1
       !     
    ENDIF
    !     
    !         Octet 22 - 23 : Number averaged.
    !     One 16 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(19),1,IBITS,16,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 408
       WRITE(*,9408) KRET
       GO TO 900
    ENDIF

    !         Octet 24 : Number missing from averages etc.
    !     One 8 bit field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC1(20),1,IBITS,8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 409
       WRITE(*,9409) KRET
       GO TO 900
    ENDIF

    !         This is the end of Section 1 , if Edition 0 or -1 of GRIB code.
    !     Set other fields to be compatible with Edition 1, where possible.
    !     
    IF(KSEC0(2)<1)THEN
       !     
       !     Century of data.
       !     
       IF(KSEC1(2)==98)THEN
          !     
          !     All ECMWF data in Edition 0 or -1 is 20th century.
          !     
          KSEC1(21) = 20
       ELSE

          !     Otherwise set century to missing data value.
          !     
          KSEC1(21) = IMISNG
       ENDIF
       !     
       !     Reserved field and decimal scale factor field(which
       !     was always 0).
       !     
       !     KSEC1(22) and KSEC1(23) already set to 0.
       !     
       GO TO 499
    ENDIF
    !     
    !         Octet 25 : Century of data.
    !     One 8 bit field.
    !     
    !     When coding, set sign bit if value is negative.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       IF(KSEC1(21)<0)THEN
          ITEMP(1) = - KSEC1(21)
          ITEMP(1) = ITEMP(1) + 32768
       ELSE
          ITEMP(1) = KSEC1(21)
       ENDIF
    ENDIF

    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,ITEMP,1,IBITS,8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 410
       WRITE(*,9410) KRET
       GO TO 900
    ENDIF

    !     When decoding, set sign bit if value is negative.
    !     
    IF(YFUNC=='D'.or.yfunc=='d') THEN
       IF(ITEMP(1)>32768)   THEN
          ITEMP(1) = ITEMP(1) - 32768
          KSEC1(21) = - ITEMP(1)
       ELSE
          KSEC1(21) = ITEMP(1)
       ENDIF
    ENDIF

    !         Octet 26 : Reserved field.(set to 0)
    !     One 8 bit field.
    !     
    !     When coding data, array KGRIB is already set to 0.
    !     When decoding, KSEC1(22) is already set to 0.
    !     Update pointer only.
    !     
    INSPT = INSPT + 8
    !     
    !         Octets 27 - 28 : Units decimal scale factor.
    !     One 16 bit field.
    !     
    !     When coding, set sign bit if value is negative.
    !     
    IF(YFUNC=='C'.or.yfunc=='c') THEN
       IF(KSEC1(23)<0)   THEN
          ITEMP(1) = - KSEC1(23)
          ITEMP(1) = ITEMP(1) + 32768
       ELSE
          ITEMP(1) = KSEC1(23)
       ENDIF
    ENDIF

    !     Insert / extract field.

    CALL INXBIT(KGRIB,KLENG,INSPT,ITEMP,1,IBITS, 16,YFUNC,KRET)
    IF(KRET/=0) THEN
       KRET = 411
       WRITE(*,9411) KRET
       GO TO 900
    ENDIF

    !     When decoding, set sign bit if value is negative.
    !     
    IF(YFUNC=='D'.or.yfunc=='d') THEN
       IF(ITEMP(1)>32768)    THEN
          ITEMP(1) = ITEMP(1) - 32768
          KSEC1(23) = - ITEMP(1)
       ELSE
          KSEC1(23) = ITEMP(1)
       ENDIF
    ENDIF

    !         Reserved octets 29-40 need not be present, and 41-nn are
    !     for local use, so when decoding, skip these octets, if
    !     presence is indicated by a length of section > 28 octets.
    !     
    IF(ILEN1(1)>28) INSPT = INSPT +(ILEN1(1) - 28) * 8

499 CONTINUE
    !     
    !         Check for ECMWF pseudo-grib data. This saves calling GRIBEX
    !     with function 'I' to check if the data is GRIB data, and another
    !     call with function 'D' when GRIB data is found.
    !     
    IF(YFUNC/='C'.and.yfunc/='c' &
         .AND.(KSEC1(6)==127.OR.KSEC1(6)==128))THEN
       IF(KSEC1(1)==128.AND.KSEC1(2)==98) THEN

          IPSEUD = -6

          !     Change function to 'L' so that section 0 is
          !     fully decoded.

          YFUNC = 'L'
          GO TO 500

       ENDIF
    ENDIF

    !         Check consistency of values given, with GRIB code, if required.

    IF(NVCK==1.AND.(YFUNC=='C'.or.yfunc=='c')) THEN
       CALL GRCHK1(KSEC1,KRET)
       IF(KRET/=0)THEN
          KRET = 499
          WRITE(*,9499) KRET
          GO TO 900
       ENDIF
    ENDIF
    !     
    !     ----------------------------------------------------------------

    !         Section 5 . Grid Description Section(Section 2) of GRIB code.
    !     ----------------------------------------------------------------
    !     
500 CONTINUE

    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 5.'

    !         Go to section 6, if no grid description included.

    IF(KSEC1(5)==0.OR.KSEC1(5)==64) THEN

       !     Set section 2 values to missing data indicator value,
       !     if decoding data.

       IF(YFUNC=='D'.or.yfunc=='d') THEN
          DO 501 J501=1,22
             KSEC2(J501) = IMISNG
501       ENDDO
       ENDIF
       GO TO 600
    ENDIF

    !         Octets 1 - 3 : Length of section.
    !     One 24 bit field.

    !     Calculate length of section, if coding data.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN

       !     parameters + vertical coordinate parameters + list of
       !     numbers of points.
       !    (Lambert conformal and Mercator are 42 octets in length,
       !     while Space view appears to be 40 !!!!!!!!!!!!!!!!!!)

       !     Ordinary Grid.

       INC = 0

       !     Space view perspective. !!!!!!!!!!! +10 ????

       IF(KSEC2(1)==90) INC = 8
       !     
       !     Rotated grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==10.OR.KSEC2(1)==14.OR.KSEC2(1)==60) INC = 10

       !     Stretched grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==20.OR.KSEC2(1)==24.OR.KSEC2(1)==70) INC = 10
       !     
       !     Stretched and rotated grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==30.OR.KSEC2(1)==34.OR.KSEC2(1)==80) INC = 20
       !     
       ILEN2(1) = 32 + INC +(KSEC2(12)*4) +(KSEC2(17)*2*KSEC2(3))
       !     
    ENDIF

    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,ILEN2,1,IBITS, 24,YFUNC,KRET)
    IF(KRET/=0) THEN
       KRET = 501
       WRITE(*,9501) KRET
       GO TO 900
    ENDIF

    !         Print length of Section 2, if required.

    IF(NDBG==1) WRITE(*,9102) ILEN2(1)

    !         If only length is required, update bit-pointer and jump
    !     to extraction of length of section 3.

    IF(YFUNC=='L'.or.yfunc=='l')THEN
       INSPT = INSPT -24 + ILEN2(1) * 8
       GO TO 600
    ENDIF

    !         Octet 4 : NV - number of vertical coordinate parameters.
    !     One 8 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(12),1,IBITS,8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 502
       WRITE(*,9502) KRET
       GO TO 900
    ENDIF

    !         Fixup for Editions -1 and 0 of GRIB code, where number
    !     of Vertical Coordinate Parameters must be calculated,
    !     as this octet contained the number of unused bits at the
    !     end of the section, which by definition of the section
    !     always had to be 0.
    !     
    !     The number of vertical parameters can only be calculated after
    !     the projection type is known.  Set it to 0 for the moment.
    IF(KSEC0(2)<1) THEN
       KSEC2(12) = 0
    ENDIF
    !     
    !         Octet 5 : PV - location of list of vertical coordinate parameters,
    !     if any,
    !     or
    !     PL - location of list of numbers of points, if no PV,
    !     or

    !     255 - no PV or PL.
    !     One 8 bit field.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       !     Set value, if coding data.
       !     
       !     Neither present is default.
       !     
       IPVPL(1) = 255
       !     
       !     Vertical coordinate parameters present.
       !     
       IF(KSEC2(12)/=0) IPVPL(1) = 32 + INC + 1
       !     
       !     List of number of points present, if no vertical
       !     coordinate parameters present and if quasi-regular grid.
       !     
       IF(KSEC2(17)==1.AND.KSEC2(12)==0) IPVPL(1) = 32 + INC + 1
       !     
       !     Insert field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IPVPL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 503
          WRITE(*,9503) KRET
          GO TO 900
       ENDIF

    ELSE

       !     If decoding data.
       !     
       !     Extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IPVPL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 503
          WRITE(*,9503) KRET
          GO TO 900
       ENDIF

       !     Experimental space view perspective data  received
       !     at ECMWF has all 0 bits.
       !     0 is illegal for all data types, so change it.
       !     
       IF(IPVPL(1)==0) IPVPL(1) = 255
       !     
       !     Neither present, so set regular grid indicator.
       !     
       IF(IPVPL(1)==255) KSEC2(17) = 0
       !     
       !     Vertical coordinate parameters present.
       !     If the length of section is greater than the
       !     end of the vertical coordinate parameters, then
       !     there is a list of numbers of points following, so
       !     set quasi-regular grid indicator.
       !     
       IF(KSEC2(12)/=0) THEN
          IPL = 4 * KSEC2(12) + IPVPL(1) - 1
          IF(IPL<ILEN2(1)) KSEC2(17) = 1
       ENDIF

       !     List of number of points present, no vertical
       !     coordinate parameters present, so set quasi-regular
       !     grid indicator.
       !     
       IF(KSEC2(12)==0.AND.IPVPL(1)/=255) KSEC2(17) = 1
       !     
       !     Fixup for Editions -1 and 0 of GRIB code, where
       !     all grids were regular.
       !     
       IF(KSEC0(2)<1)  THEN
          KSEC2(17) = 0
       ENDIF
       !     
    ENDIF

    !         Octet 6 : Data representation type.
    !     One 8 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(1),1,IBITS, 8,YFUNC,KRET)
    IF(KRET/=0) THEN
       KRET = 504
       WRITE(*,9504) KRET
       GO TO 900
    ENDIF


    !     Calculate length of section when decoding
    !     
    IF(YFUNC=='D'.or.yfunc=='d') THEN
       !     
       !     parameters + vertical coordinate parameters + list of
       !     numbers of points.
       !    (Lambert conformal and Mercator are 42 octets in length,
       !     while Space view appears to be 40 !!!!!!!!!!!!!!!!!!)
       !     
       !     Ordinary Grid.
       !     
       INC = 0
       !     
       !     Space view perspective. !!!!!!!!!!! +10 ????
       !     
       IF(KSEC2(1)==90) INC = 8
       !     
       !     Rotated grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==10.OR.KSEC2(1)==14.OR.KSEC2(1)==60) INC = 10

       !     Stretched grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==20.OR.KSEC2(1)==24.OR.KSEC2(1)==70) INC = 10
       !     
       !     Stretched and rotated grid.
       !    (Gaussian, Latitude/longitude or Spherical Harmonics)
       !     
       IF(KSEC2(1)==30.OR.KSEC2(1)==34.OR.KSEC2(1)==80) INC = 20
       !     
       !     Now calculate number of vertical coordinate parameters
       !     
       IF(KSEC0(2)<1) THEN
          KSEC2(12) =( ILEN2(1) - 32 - INC ) / 4
       ENDIF
    ENDIF

    !         Gaussian grid definition.
    !     
    IF(KSEC2(1)==4.OR.KSEC2(1)==14.OR.KSEC2(1)==24 &
         .OR.KSEC2(1)==34) THEN
       !     
       !                Octets 7 - 8  : Ni - number of points along a parallel.
       !     Octets 9 - 10 : Nj - number of points along a meridian.
       !     Two 16 bit fields.
       !     
       !     For quasi-regular grids Ni is set to all 1 bits, as
       !     the number of points is different on different parallels.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN

          !     When coding, set to all 1 bits.

          IF(KSEC2(17)==1) KSEC2(2) = 65535
       ENDIF
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(2),2,IBITS,16,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 505
          WRITE(*,9505) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     When decoding, change 1 bits to missing data value.

          IF(KSEC2(2)==65535) KSEC2(2) = IMISNG
       ENDIF

       !                Octets 11 - 13 : La1 - latitude of first grid point.
       !     Octets 14 - 16 : Lo1 - longitude of first grid point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.

       IF(YFUNC=='C'.or.yfunc=='c')then
          ILALO(1) = KSEC2(4)
          ILALO(2) = KSEC2(5)
          IF(KSEC2(4)<0) ILALO(1) = -(KSEC2(4)) + 8388608
          IF(KSEC2(5)<0) ILALO(2) = -(KSEC2(5)) + 8388608
       ENDIF

       !     Insert / extract fields.

       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 506
          WRITE(*,9506) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.

       IF(YFUNC=='D'.or.yfunc=='d')then
          KSEC2(4) = ILALO(1)
          KSEC2(5) = ILALO(2)
          IF(KSEC2(4)>8388608)KSEC2(4) = -(KSEC2(4)-8388608)
          IF(KSEC2(5)>8388608)KSEC2(5) = -(KSEC2(5)-8388608)
       ENDIF

       !                Octet 17 : Resolution and components flag.
       !     One 8 bit field.

       IF(YFUNC=='C'.or.yfunc=='c')then
          IRESOL(1) = KSEC2(6)+KSEC2(18)+KSEC2(19)
       endif
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IRESOL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 507
          WRITE(*,9507) KRET
          GO TO 900
       ENDIF

       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     All flag fields are already set to 0, so

          IF(IRESOL(1)==0) GO TO 510

          !     Fix up for flag which was different in Experimental
          !     edition.

          IF(KSEC0(2)==-1.AND.(IRESOL(1)==1.OR.IRESOL(1)== &
               3)) IRESOL(1) = 128
          !     
          !     Set Resolution flag.
          !     
          IF(IRESOL(1)>=128)THEN
             KSEC2(6) = 128
             IRESOL(1)   = IRESOL(1) - 128
          ENDIF

          !     Set earth flag.
          !     
          IF(IRESOL(1)>=64)THEN
             KSEC2(18) = 64
             IRESOL(1)    = IRESOL(1) - 64
          ENDIF
          !     
          !     Set components flag.
          !     
          KSEC2(19) = IRESOL(1)
          !     
       ENDIF
       !     
510    CONTINUE

       !                Octets 18 - 20 : La2 - latitude of last grid point.
       !     Octets 21 - 23 : Lo2 - longitude of last grid point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN
          ILALO(1) = KSEC2(7)
          ILALO(2) = KSEC2(8)
          IF(KSEC2(7)<0) ILALO(1) = -(KSEC2(7)) + 8388608
          IF(KSEC2(8)<0) ILALO(2) = -(KSEC2(8)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS, 24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 508
          WRITE(*,9508) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(7) = ILALO(1)
          KSEC2(8) = ILALO(2)
          IF(KSEC2(7)>8388608)KSEC2(7) = -(KSEC2(7)-8388608)
          IF(KSEC2(8)>8388608)KSEC2(8) = -(KSEC2(8)-8388608)
       ENDIF

       !                Octets 24 - 25 : Di - i direction increment.
       !     One 16 bit field.
       !     
       !     For quasi-regular grids all Di bits are set to 1, as
       !     the increment is different on different parallels.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     When coding, set to all 1 bits.
          !     
          IF(KSEC2(17)==1) KSEC2(9) = 65535
          !     
          !     If increments not given, set all bits to 1.
          !     
          IF(KSEC2(6)==0) KSEC2(9) = 65535
       ENDIF
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(9),1,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 509
          WRITE(*,9509) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     When decoding, change 1 bits to missing data value.
          !     
          IF(KSEC2(9)==65535) KSEC2(9) = IMISNG
       ENDIF
       !     
       !                Octets 26 - 27 : N- number of parallels between a Pole
       !     and the Equator.
       !     One 16 bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(10),1,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 510
          WRITE(*,9510) KRET
          GO TO 900
       ENDIF

       !                Octet 28 : Scanning mode flags.
       !     One 8 bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(11),1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 511
          WRITE(*,9511) KRET
          GO TO 900
       ENDIF

       !     Fix-up for flag which was different in Experimental
       !     Edition.
       !     
       IF(KSEC0(2)==-1.AND.KSEC2(11)==1) KSEC2(11) = 0
       !     
       !                Octets 29 - 32 : Reserved.
       !     Two 16 bit fields.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')  THEN
          !     
          !     All bits already set to 0.
          !     No insertion, only update bit pointer.
          !     
          INSPT = INSPT + 32
          !     
       ELSE
          !     
          !     No extraction, only update bit pointer.
          !     
          INSPT = INSPT + 32
          !     
       ENDIF
       !     
       GO TO 520
    ENDIF

    !         Latitude/longitude grid definition,
    !     Equidistant Cylindrical or Plate Carree.
    !     
    IF(KSEC2(1)==0.OR.KSEC2(1)==10.OR.KSEC2(1)==20 &
         .OR.KSEC2(1)==30)THEN
       !     
       !                Octets 7 - 8  : Ni - number of points along a parallel.
       !     Octets 9 - 10 : Nj - number of points along a meridian.
       !     Quasi-regular grids not catered for.
       !     Two 16 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(2),2,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 505
          WRITE(*,9505) KRET
          GO TO 900
       ENDIF

       !                Octets 11 - 13 : La1 - latitude of first grid point.
       !     Octets 14 - 16 : Lo1 - longitude of first grid point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          ILALO(1) = KSEC2(4)
          ILALO(2) = KSEC2(5)
          IF(KSEC2(4)<0) ILALO(1) = -(KSEC2(4)) + 8388608
          IF(KSEC2(5)<0) ILALO(2) = -(KSEC2(5)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 506
          WRITE(*,9506) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN
          KSEC2(4) = ILALO(1)
          KSEC2(5) = ILALO(2)
          IF(KSEC2(4)>8388608)KSEC2(4) = -(KSEC2(4)-8388608)
          IF(KSEC2(5)>8388608)KSEC2(5) = -(KSEC2(5)-8388608)
       ENDIF

       !                Octet 17 : Resolution and components flag.
       !     One 8 bit field.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')then
          IRESOL(1) = KSEC2(6)+KSEC2(18)+KSEC2(19)
       endif
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IRESOL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 507
          WRITE(*,9507) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     All flag fields are already set to 0, so
          !     
          IF(IRESOL(1)==0) GO TO 511
          !     
          !     Fix up for flag which was different in Experimental
          !     edition.
          !     
          IF(KSEC0(2)==-1.AND.(IRESOL(1)==1.OR.IRESOL(1)== &
               3)) IRESOL(1) = 128
          !     
          !     Set Resolution flag.
          !     
          IF(IRESOL(1)>=128)THEN
             KSEC2(6) = 128
             IRESOL(1)   = IRESOL(1) - 128
          ENDIF

          !     Set earth flag.
          !     
          IF(IRESOL(1)>=64)THEN
             KSEC2(18) = 64
             IRESOL(1)    = IRESOL(1) - 64
          ENDIF
          !     
          !     Set components flag.
          !     
          KSEC2(19) = IRESOL(1)
          !     
       ENDIF
       !     
511    CONTINUE

       !                Octets 18 - 20 : La2 - latitude of last grid point.
       !     Octets 21 - 23 : Lo2 - longitude of last grid point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          ILALO(1) = KSEC2(7)
          ILALO(2) = KSEC2(8)
          IF(KSEC2(7)<0) ILALO(1) = -(KSEC2(7)) + 8388608
          IF(KSEC2(8)<0) ILALO(2) = -(KSEC2(8)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 508
          WRITE(*,9508) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(7) = ILALO(1)
          KSEC2(8) = ILALO(2)
          IF(KSEC2(7)>8388608)KSEC2(7) = -(KSEC2(7)-8388608)
          IF(KSEC2(8)>8388608)KSEC2(8) = -(KSEC2(8)-8388608)
       ENDIF

       !                Octets 24 - 25 : Di - i direction increment.
       !     One 16 bit field.
       !     
       !     Quasi-regular grids not catered for.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     If field not given, set to all bits to 1.
          !     
          IF(KSEC2(6)==0) KSEC2(9) = 65535
       ENDIF
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(9),1,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 509
          WRITE(*,9509) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN

          !     When decoding, change 1 bits to missing data value.
          !     
          IF(KSEC2(9)==65535) KSEC2(9) = IMISNG
       ENDIF
       !     
       !                Octets 26 - 27 : Dj - j direction increment.
       !     One 16 bit field.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     If field not given, set to all bits to 1.
          !     
          IF(KSEC2(6)==0) KSEC2(10) = 65535
       ENDIF

       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(10),1,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 513
          WRITE(*,9513) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN

          !     When decoding, change 1 bits to missing data value.
          !     
          IF(KSEC2(10)==65535) KSEC2(10) = IMISNG
       ENDIF
       !     
       !                Octet 28 : Scanning mode flags.
       !     One 8 bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(11),1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 511
          WRITE(*,9511) KRET
          GO TO 900
       ENDIF

       !     Fix-up for flag which was different in Experimental
       !     Edition.
       !     
       IF(KSEC0(2)==-1.AND.KSEC2(11)==1) KSEC2(11) = 0
       !     
       !                Octets 29 - 32 : Reserved.
       !     Two 16 bit fields.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     All bits already set to 0.
          !     No insertion, only update bit pointer.
          !     
          INSPT = INSPT + 32
          !     
       ELSE
          !     
          !     No extraction, only update bit pointer.
          !     
          INSPT = INSPT + 32
       ENDIF
       !     
       GO TO 520
    ENDIF

    !         Spherical Harmonic format.
    !     
    IF(KSEC2(1)==50.OR.KSEC2(1)==60.OR.KSEC2(1)==70 &
         .OR.KSEC2(1)==80)THEN
       !     
       !                Octets 7 - 8   : J pentagonal resolution parameter.
       !     Octets 9 - 10  : K pentagonal resolution parameter.
       !     Octets 11 - 12 : M pentagonal resolution parameter.
       !     Three 16 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(2),3,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 514
          WRITE(*,9514) KRET
          GO TO 900
       ENDIF

       !                Octet 13 : Representation type.
       !     Octet 14 : Representation mode.
       !     Two 8 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(5),2,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 515
          WRITE(*,9515) KRET
          GO TO 900
       ENDIF

       !                Octets 15 - 32 : Reserved.
       !     Nine 16 bit fields.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')  THEN
          !     
          !     All bits already set to 0.
          !     No insertion, only update bit pointer.
          !     
          INSPT = INSPT + 144
          !     
       ELSE
          !     
          !     No extraction, only update bit pointer.
          !     KSEC2(7) to KSEC2(11) already set to 0.
          !     
          INSPT = INSPT + 144
       ENDIF
       !     
       GO TO 520
    ENDIF

    !         Polar Stereographic.
    !     
    IF(KSEC2(1)==5) THEN
       !     
       !                Octets 7 - 8  : Ni - number of points along X-axis.
       !     Octets 9 - 10 : Nj - number of points along Y-axis.
       !     Two 16 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(2),2,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 523
          WRITE(*,9523) KRET
          GO TO 900
       ENDIF

       !                Octets 11 - 13 : La1 - latitude of first grid point.
       !     Octets 14 - 16 : Lo1 - longitude of first grid point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          ILALO(1) = KSEC2(4)
          ILALO(2) = KSEC2(5)
          IF(KSEC2(4)<0) ILALO(1) = -(KSEC2(4)) + 8388608
          IF(KSEC2(5)<0) ILALO(2) = -(KSEC2(5)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 506
          WRITE(*,9506) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(4) = ILALO(1)
          KSEC2(5) = ILALO(2)
          IF(KSEC2(4)>8388608)KSEC2(4) = -(KSEC2(4)-8388608)
          IF(KSEC2(5)>8388608)KSEC2(5) = -(KSEC2(5)-8388608)
       ENDIF

       !                Octet 17 : Resolution and components flag.
       !     One 8 bit field.
       !     
       !     Resolution flag( KSEC2(6) ) is not applicable.
       !     
       KSEC2(6) = 0
       IF(YFUNC=='C'.or.yfunc=='c')IRESOL(1)=KSEC2(18)+KSEC2(19)
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IRESOL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 507
          WRITE(*,9507) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     All flag fields are already set to 0, so
          !     
          IF(IRESOL(1)==0) GO TO 513
          !     
          !     Fix up for flag which was different in Experimental
          !     edition.
          !     
          IF(KSEC0(2)==-1.AND.(IRESOL(1)==1.OR.IRESOL(1)== &
               3)) IRESOL(1) = 128
          !     
          !     Resolution flag is not applicable.
          !     
          IF(IRESOL(1)>=128) IRESOL(1) = IRESOL(1) - 128
          !     
          !     Set earth flag.
          !     
          IF(IRESOL(1)>=64)THEN
             KSEC2(18) = 64
             IRESOL(1)    = IRESOL(1) - 64
          ENDIF

          !     Set components flag.
          !     
          KSEC2(19) = IRESOL(1)
          !     
       ENDIF
       !     
513    CONTINUE

       !                Octets 18 - 20 : LoV - orientation of the grid.
       !     One 24 bit field.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          ILALO(1) = KSEC2(7)
          IF(KSEC2(7)<0) ILALO(1) = -(KSEC2(7)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),1,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 508
          WRITE(*,9508) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(7) = ILALO(1)
          IF(KSEC2(7)>8388608)KSEC2(7) = -(KSEC2(7)-8388608)
       ENDIF
       !     
       !                Octets 21 - 23 : Dx - X direction grid length.
       !                Octets 24 - 26 : Dy - Y direction grid length.
       !     Two 24 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(9),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 524
          WRITE(*,9524) KRET
          GO TO 900
       ENDIF

       !                Octet 27 : Projection centre flag.
       !     One 8-bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(13),1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 525
          WRITE(*,9525) KRET
          GO TO 900
       ENDIF

       !                Octet 28 : Scanning mode flags.
       !     One 8 bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(11),1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 511
          WRITE(*,9511) KRET
          GO TO 900
       ENDIF

       !     Fix-up for flag which was different in Experimental
       !     Edition.
       !     
       IF(KSEC0(2)==-1.AND.KSEC2(11)==1) KSEC2(11) = 0
       !     
       !                Octets 29 - 32 : Reserved.
       !     Two 16 bit fields.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     All bits already set to 0.
          !     No insertion, only update bit pointer.
          !     
          INSPT = INSPT + 32
          !     
       ELSE
          !     
          !     No extraction, only update bit pointer.
          !     
          INSPT = INSPT + 32
          !     
       ENDIF
       !     
       GO TO 520
    ENDIF

    !         Space view perspective or orthographic.
    !     
    IF(KSEC2(1)==90)THEN
       !     
       !                Octets 7 - 8  : Nx - number of points along X-axis.
       !     Octets 9 - 10 : Ny - number of points along Y-axis.
       !     Two 16 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(2),2,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 523
          WRITE(*,9523) KRET
          GO TO 900
       ENDIF

       !                Octets 11 - 13 : Lap - latitude of sub-satellite point.
       !     Octets 14 - 16 : Lop - longitude of sub-satellite point.
       !     Two 24 bit fields.
       !     
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN
          ILALO(1) = KSEC2(4)
          ILALO(2) = KSEC2(5)
          IF(KSEC2(4)<0) ILALO(1) = -(KSEC2(4)) + 8388608
          IF(KSEC2(5)<0) ILALO(2) = -(KSEC2(5)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 526
          WRITE(*,9526) KRET
          GO TO 900
       ENDIF

       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(4) = ILALO(1)
          KSEC2(5) = ILALO(2)
          IF(KSEC2(4)>8388608)KSEC2(4) = -(KSEC2(4)-8388608)
          IF(KSEC2(5)>8388608)KSEC2(5) = -(KSEC2(5)-8388608)
       ENDIF

       !                Octet 17 : Resolution and components flag.
       !     One 8 bit field.
       !     
       !     Resolution flag( KSEC2(6) ) is not applicable.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')IRESOL(1)=KSEC2(18)+KSEC2(19)
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IRESOL,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 507
          WRITE(*,9507) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN

          !     All flag fields are already set to 0, so
          !     
          IF(IRESOL(1)==0) GO TO 514
          !     
          !     Resolution flag is not applicable.
          !     
          IF(IRESOL(1)>=128) IRESOL(1) = IRESOL(1) - 128
          !     
          !     Set earth flag.
          !     
          IF(IRESOL(1)>=64)THEN
             KSEC2(18) = 64
             IRESOL(1)    = IRESOL(1) - 64
          ENDIF

          !     Set components flag.
          !     
          KSEC2(19) = IRESOL(1)
       ENDIF
       !     
514    CONTINUE
       !     
       !                Octets 18 - 20 : dx Apparent diameter of earth in grid
       !     lengths in x direction.
       !     Octets 21 - 23 : dy Apparent diameter of earth in grid
       !     lengths in y direction.
       !     Two 24 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(7),2,IBITS, 24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 527
          WRITE(*,9527) KRET
          GO TO 900
       ENDIF

       !                Octets 24 - 25 : Xp X-coordinate of sub-satellite point.
       !     Octets 26 - 27 : Yp Y-coordinate of sub-satellite point.
       !     Two 16 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(9),2,IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 528
          WRITE(*,9528) KRET
          GO TO 900
       ENDIF

       !                Octet 28 : Scanning mode flags.
       !     One 8 bit field.
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(11),1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 511
          WRITE(*,9511) KRET
          GO TO 900
       ENDIF

       !                Octets 29 - 31 : The orientation of the grid.
       !     Octets 32 - 34 : nr the altitude of the camera.
       !     Two 24 bit fields.
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(13),2,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 529
          WRITE(*,9529) KRET
          GO TO 900
       ENDIF

       !                Octets 35 - 40 : Reserved.   !!!!!!! 42 ?????????
       !     Three 16 bit fields.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          !     
          !     All bits already set to 0.
          !     No insertion, only update bit pointer.
          !     
          INSPT = INSPT + 48
          !     
       ELSE
          !     
          !     No extraction, only update bit pointer.
          !     
          INSPT = INSPT + 48
          !     
       ENDIF
       !     
       GO TO 520
       !     
    ENDIF

    !         Other representation types not yet catered for.
    !     
    KRET = 598
    WRITE(*,9598) KSEC2(1) , KRET
    GO TO 900
    !     
    !         Rotation parameters for rotated or stretched and rotated grids.
    !    (Gaussian, Latitude/longitude or Spherical Harmonics)
    !     
520 CONTINUE
    !     
    IF(KSEC2(1)==10.OR.KSEC2(1)==30.OR.&
         KSEC2(1)==14.OR.KSEC2(1)==34.OR.&
         KSEC2(1)==60.OR.KSEC2(1)==80)THEN

       !                Octets 33 - 35 : Latitude of the southern pole.
       !     Octets 36 - 38 : Longitude of the southern pole.
       !     Two 24 bit fields.
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c') THEN
          ILALO(1) = KSEC2(13)
          ILALO(2) = KSEC2(14)
          IF(KSEC2(13)<0) ILALO(1) = -(KSEC2(13)) + 8388608
          IF(KSEC2(14)<0) ILALO(2) = -(KSEC2(14)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,   &
            24,YFUNC,KRET)
       !     
       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN
          KSEC2(13) = ILALO(1)
          KSEC2(14) = ILALO(2)
          IF(KSEC2(13)>8388608)KSEC2(13) = -(KSEC2(13)-8388608)
          IF(KSEC2(14)>8388608)KSEC2(14) = -(KSEC2(14)-8388608)
       ENDIF

       IF(KRET/=0)THEN
          KRET = 517
          WRITE(*,9517) KRET
          GO TO 900
       ENDIF

       !                Octets 39 - 42 : Angle of rotation.
       !     One 8 bit and one 24 bit field.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN
          !     
          !     Convert floating point to GRIB representation.
          !     
          ITRND = 1
          CALL CONFP3(PSEC2(1),IEXP(1),IMANT(1),IBITS,ITRND)
          !     
          !     Insert fields.
          !     
          CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,8,YFUNC,KRET)
          IF(KRET/=0)THEN
             KRET = 518
             WRITE(*,9518) KRET
             GO TO 900
          ENDIF
          CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,24,YFUNC,KRET)
          IF(KRET/=0)THEN
             KRET = 518
             WRITE(*,9518) KRET
             GO TO 900
          ENDIF
          !     
       ELSE

          !     Extract fields.
          !     
          CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,8,YFUNC,KRET)
          IF(KRET/=0)THEN
             KRET = 518
             WRITE(*,9518) KRET
             GO TO 900
          ENDIF
          CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,24,YFUNC,KRET)
          IF(KRET/=0)THEN
             KRET = 518
             WRITE(*,9518) KRET
             GO TO 900
          ENDIF

          !     Convert GRIB representation to floating point.
          !     
          CALL DECFP2(PSEC2(1),IEXP(1),IMANT(1))
          !     
       ENDIF
    ENDIF
    !     
    !         Stretching parameters for stretched grids.
    !    (Gaussian, Latitude/longitude or Spherical Harmonics)
    !     
    IF(KSEC2(1)==20.OR.KSEC2(1)==24.OR.KSEC2(1)==70)THEN

       !                Octets 33 - 35 : Latitude of pole of stretching.
       !     Octets 36 - 38 : Longitude of pole of stretching.
       !     Two 24 bit fields.
       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN
          ILALO(1) = KSEC2(13)
          ILALO(2) = KSEC2(14)
          IF(KSEC2(13)<0) ILALO(1) = -(KSEC2(13)) + 8388608
          IF(KSEC2(14)<0) ILALO(2) = -(KSEC2(14)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       !     
       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          KSEC2(13) = ILALO(1)
          KSEC2(14) = ILALO(2)
          IF(KSEC2(13)>8388608)KSEC2(13) = -(KSEC2(13)-8388608)
          IF(KSEC2(14)>8388608)KSEC2(14) = -(KSEC2(14)-8388608)
       ENDIF

       IF(KRET/=0)THEN
          KRET = 519
          WRITE(*,9519) KRET
          GO TO 900
       ENDIF

       !                Octets 39 - 42 : Stretching factor.
       !     One 8 bit and one 24 bit field.
       !     
       ITRND = 1
       IF(YFUNC=='C'.or.yfunc=='c')then
          !     
          !     Convert floating point to GRIB representation.
          !     
          CALL CONFP3(PSEC2(2),IEXP(1),IMANT(1),IBITS,ITRND)
       endif
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 520
          WRITE(*,9520) KRET
          GO TO 900
       ENDIF
       CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 520
          WRITE(*,9520) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')then
          !     Convert GRIB representation to floating point.
          CALL DECFP2(PSEC2(2),IEXP(1),IMANT(1))
       endif
       !     
    ENDIF
    !     
    !         Stretching parameters for stretched and rotated grids.
    !    (Gaussian, Latitude/longitude or Spherical Harmonics)
    !     
    IF(KSEC2(1)==30.OR.KSEC2(1)==34.OR.KSEC2(1)==80)THEN
       !     
       !                Octets 43 - 45 : Latitude of pole of stretching.
       !     Octets 46 - 48 : Longitude of pole of stretching.
       !     Two 24 bit fields.

       !     When coding data, set sign bit to 1, if value is
       !     negative.
       !     
       IF(YFUNC=='C'.or.yfunc=='c')THEN
          ILALO(1) = KSEC2(13)
          ILALO(2) = KSEC2(14)
          IF(KSEC2(13)<0) ILALO(1) = -(KSEC2(13)) + 8388608
          IF(KSEC2(14)<0) ILALO(2) = -(KSEC2(14)) + 8388608
       ENDIF

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILALO(1),2,IBITS,24,YFUNC,KRET)
       !     
       !     When decoding data, if sign bit is 1, value is
       !     negative.
       !     
       IF(YFUNC=='D'.or.yfunc=='d') THEN
          KSEC2(13) = ILALO(1)
          KSEC2(14) = ILALO(2)
          IF(KSEC2(13)>8388608)KSEC2(13) = -(KSEC2(13)-8388608)
          IF(KSEC2(14)>8388608)KSEC2(14) = -(KSEC2(14)-8388608)
       ENDIF
       IF(KRET/=0)THEN
          KRET = 519
          WRITE(*,9519) KRET
          GO TO 900
       ENDIF

       !                Octets 49 - 52 : Stretching factor.
       !     One 8 bit and one 24 bit field.
       !     
       ITRND = 1
       IF(YFUNC=='C'.or.yfunc=='c')then
          !     Convert floating point to GRIB representation.
          CALL CONFP3(PSEC2(2),IEXP(1),IMANT(1),IBITS,ITRND)
       endif
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 520
          WRITE(*,9520) KRET
          GO TO 900
       ENDIF
       CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 520
          WRITE(*,9520) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')then
          !     Convert GRIB representation to floating point.
          CALL DECFP2(PSEC2(2),IEXP(1),IMANT(1))
       endif

    ENDIF
    !     
    !         Vertical coordinate parameters, if any.
    !     
    IF(KSEC2(12)/=0)THEN
       ITRND = 1
       DO 530 J530 = 1 , KSEC2(12)

          !     One 8 bit and one 24 bit field.
          !     
          IF(YFUNC=='C'.or.yfunc=='c')then
             !     Convert floating point to GRIB representation.
             CALL CONFP3(PSEC2(J530+10),IEXP(1),IMANT(1),IBITS,ITRND)
          endif
          !     
          !     Insert / extract fields.
          !     
          CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,8,YFUNC,KRET)
          IF(KRET/=0) THEN
             KRET = 521
             WRITE(*,9521) KRET
             GO TO 900
          ENDIF
          CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,24,YFUNC,KRET)
          IF(KRET/=0)THEN
             KRET = 521
             WRITE(*,9521) KRET
             GO TO 900
          ENDIF
          !     
          IF(YFUNC=='D'.or.yfunc=='d')then
             !     Convert GRIB representation to floating point.
             CALL DECFP2(PSEC2(J530+10),IEXP(1),IMANT(1))
          endif
          !     
530    enddo
    ENDIF
    !     
    !         List of number of points, if any.
    !     Number of 16 bit fields.
    !     
    IF(KSEC2(17)==1)THEN
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,KSEC2(23),KSEC2(3),IBITS,16,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 522
          WRITE(*,9522) KRET
          GO TO 900
       ENDIF
    ENDIF

    !         Check consistency of values given, with GRIB code, if required.
    !     
    IF(NVCK==1.AND.(YFUNC=='C'.or.yfunc=='c'))THEN
       CALL GRCHK2(KSEC1,KSEC2,KRET)
       IF(KRET/=0)THEN
          KRET = 599
          WRITE(*,9599) KRET
          GO TO 900
       ENDIF
    ENDIF
    !     
    !     ----------------------------------------------------------------

    !         Section 6 . Bit Map Section(section 3) of GRIB code.
    !     ----------------------------------------------------------------
    !     
600 CONTINUE
    !     
    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 6.'
    !     
    !         Go to section 9, if decoding of identification sections only and
    !     GRIB Code Edition is higher than 0. If Edition is lower the
    !     length of the GRIB message needs to be calculated, so change
    !     function to 'L' to complete decoding of section 0.
    !     Number of data values decoded( KSEC4(1) ) already set to 0.
    !     
    IF(HOPER(1:1)=='I'.or.HOPER(1:1)=='i')THEN
       IF(KSEC0(2)>0) THEN
          GO TO 900
       ELSE
          YFUNC = 'L'
       ENDIF
    ENDIF

    !         Go to section 7, if no bit map required.
    !     
    IF(KSEC1(5)==0.OR.KSEC1(5)==128) GO TO 700
    !     
    !         Set bit-map flag and attempt no decoding of bit-map, if
    !     routine has been called by the DECOGB interface routine,
    !     which is provided for upward compatibility with old software.
    !     
    IF(HOPER=='Z'.or.HOPER=='z')THEN
       ISBMAP = -5
       WRITE(*,9606) ISBMAP
       GO TO 900
    ENDIF

    !         When coding data, calculate the length of section and number
    !     of unused bits.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       IF(KSEC3(1)/=0)THEN
          !     
          !     Predetermined bit-map table included.
          !     Length of section is 6 octets, number of unused
          !     bits is 0.
          !     
          ILEN3(1) = 6
          INUB(1)  = 0
          !     
       ELSE

          !     Bit-map included in section 3.
          !     Length of section - 6 octets of header + length of
          !     bit-map, rounded to a multiple of 2 octets.
          !     
          ITEMP(1) = 48 + ILENF
          ILEN3(1) =( ITEMP(1) + 15 ) / 16
          ILEN3(1) = ILEN3(1) * 2
          !     
          !     Number of unused bits.
          !     
          INUB(1) = ILEN3(1) * 8 - ITEMP(1)
          !     
       ENDIF
       !     
    ENDIF

    !         Octets 1 - 3 : Length of section.
    !     One 24 bit field.
    !     
    !     Insert/extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,ILEN3,1,IBITS,24,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 601
       WRITE(*,9601) KRET
       GO TO 900
    ENDIF

    !         Print length of Section 3, if required.
    !     
    IF(NDBG==1) WRITE(*,9103) ILEN3(1)
    !     
    !     If only length is required, update bit-pointer and jump
    !     to extraction of length of section 4.
    !     
    IF(YFUNC=='L'.or.YFUNC=='l')THEN
       INSPT = INSPT -24 + ILEN3(1) * 8
       GO TO 700
    ENDIF

    !         Octet 4 : Number of unused bits at end of section.
    !     One 8 bit field.
    !     
    !     Insert/extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,INUB,1,IBITS,&
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 602
       WRITE(*,9602) KRET
       GO TO 900
    ENDIF

    !         Octets 5-6 : Bit-map table reference.
    !     One 16 bit field.
    !     
    !     Insert/extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC3(1),1,IBITS,     &
         16,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 603
       WRITE(*,9603) KRET
       GO TO 900
    ENDIF

    !         Finished if a predetermined bit-map table is given.

    IF(KSEC3(1)/=0) GO TO 699

    !         Bit-map definition included.

    IF(YFUNC=='C'.or.yfunc=='c')THEN

       !     IVALS is the number of bits in the bit-map. It is the same
       !     as the number of data values(including missing data
       !     values) to be coded.

       IVALS = ILENF

       !     Set integer or real missing data value.
       !     ZMSVAL and IMSVAL are euqivalenced.

       IF(KSEC4(5)==0)THEN
          ZMSVAL = PSEC3(2)
       ELSE
          IMSVAL = KSEC3(2)
       ENDIF

       !     Insert bit-map. Set function for fixed length messages,
       !     if required.

       YTEMP = YFUNC
       IF(HOPER=='M'.or.hoper=='m') YTEMP = 'M'

       IBMAP = INSPT

       CALL INXMAP(KGRIB,KLENG,INSPT,FIELD,ILENF,IBITS,     &        
            ISBMAP,ZMSVAL,YTEMP,NDBG,KRET)
       IF(KRET/=0)THEN
          KRET = 604
          WRITE(*,9604) KRET
          GO TO 900
       ENDIF

       !     Number of data values(not including missing data values)
       !     is now in ILENF, which is used when finding maximum and
       !     minimum values etc.
       !     
       !                Unused bits at end of section.
       !     These bits are already set to 0, so update bit-pointer
       !     only.
       !     
       INSPT = INSPT + INUB(1)
       !     
    ELSE
       !     
       !     Retain pointer to bit-map location.
       !     
       IBMAP = INSPT
       !     
       !     IVALS is the number of bits in the bit-map. It is the same
       !     as the number of data values(including missing data
       !     values to be decoded.
       !     
       IVALS =(ILEN3(1) - 6) * 8 - INUB(1)
       !     
       !     Update bit-pointer to start of section 4 of Grib message.
       !     
       INSPT = INSPT - 48 + ILEN3(1) * 8
    ENDIF

    !         Check consistency of values given, with GRIB code, if required.
    !     
699 CONTINUE
    !     
    IF(NVCK==1.AND.(YFUNC=='C'.or.yfunc=='c'))THEN
       CALL GRCHK3(KSEC1,KSEC3,KRET)
       IF(KRET/=0)THEN
          KRET = 699
          WRITE(*,9699) KRET
          GO TO 900
       ENDIF
    ENDIF
    !     
    !     ----------------------------------------------------------------

    !         Section 7 . Binary Data Section(section 4) of GRIB code.
    !     ----------------------------------------------------------------
    !     
700 CONTINUE
    !     
    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 7.'
    !     
    !         Octets 1 - 3 : Length of section.
    !     One 24 bit field.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       !     Retain pointer to location, for later
       !     insertion, if coding data.
       !     Increment pointer.
       !     
       IPLEN = INSPT
       INSPT = INSPT + 24
    ELSE

       !     Extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ILEN4,1,IBITS,     &        
            24,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 701
          WRITE(*,9701) KRET
          GO TO 900
       ENDIF

       !     Print length of Section 4, if required.
       !     
       IF(NDBG==1) WRITE(*,9104) ILEN4(1)
       !     
       !     If only length is required, update bit-pointer and add
       !     length of section 5(32 bits). Return length in bytes.
       !     
       IF(YFUNC=='L'.or.YFUNC=='l')THEN
          INSPT    = INSPT - 24 + ILEN4(1) * 8 + 32
          KSEC0(1) = INSPT / 8
          GO TO 900
       ENDIF
       !     
    ENDIF

    !         Octet 4 : 4 bit flag field and 4 bit unused bit count field.
    !     One 8 bit field for insertion/extraction purposes.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       !     Only simple packing of data(real or integer) with no
       !     additional flag field supported.
       !     
       IF(KSEC4(4)/=0.OR.KSEC4(6)/=0)THEN
          KRET = 705
          WRITE(*,9705) KRET
          GO TO 900
       ENDIF

       !     Type of data(spherical harmonic coefficients or grid
       !     point) is taken from KSEC4(3) only if no Section 2 is
       !     included. This allows coding of data without the use
       !     of Section 2.
       !     
       IREP = 0
       IF(KSEC1(5)==0.OR.KSEC1(5)==64)THEN
          IFLAG(1) = KSEC4(3)
       ELSE
          IFLAG(1) = 0
          IF( KSEC2(1)==50.OR.KSEC2(1)==60.OR.     &           
               KSEC2(1)==70.OR.KSEC2(1)==80) IFLAG(1) = 128
       ENDIF
       !     
       IF(IFLAG(1)/=0) IREP = 1
       IFLAG(1) = IFLAG(1) + KSEC4(4) + KSEC4(5) + KSEC4(6)

       !     When coding data, field is inserted later, when
       !     number of unused bits is known and added to it.
       !     Increment pointer.
       !     
       INSPT = INSPT + 8
       !     
    ELSE
       !     
       !     Extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IFLAG,1,IBITS,     &        
            8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 706
          WRITE(*,9706) KRET
          GO TO 900
       ENDIF

       !     All flags already preset to 0.
       !     
       IF(KSEC0(2)==-1)THEN
          !     
          !     In the Experimental Edition flag field was
          !     0000 for grid point data.
          !     0001 for spherical harmonic data.
          !     
          IREP = IFLAG(1) / 16
          INIL = IFLAG(1) - IREP * 16
          GO TO 710
       ENDIF
       !     
       IF(KSEC0(2)==0)THEN

          !     In Edition 0 flag field was
          !     0000 for grid point data.
          !     1000 for spherical harmonic data.
          !     
          IREP = IFLAG(1) / 128
          INIL = IFLAG(1) - IREP * 128
          GO TO 710
       ENDIF
       !     
       IF(KSEC0(2)==1) THEN

          !     In Edition 1 flag field all 4 bits have
          !     significance.
          !     
          !     0--- for grid point data.
          !     1--- for spherical harmonic data.
          !     -0-- for simple packing.
          !     -1-- for complex or second order packing.
          !     --0- for floating point values.
          !     --1- for integer values.
          !     ---0 for no additional flags at Octet 14.
          !     ---1 for additional flags at Octet 14.
          !     
          IF(IFLAG(1)>=128) THEN
             IREP     = 1
             KSEC4(3) = 128
             IFLAG(1)    = IFLAG(1) - 128
          ELSE
             IREP     = 0
          ENDIF
          !     
          IF(IFLAG(1)>=64)THEN
             KSEC4(4) = 64
             IFLAG(1)    = IFLAG(1) - 64
          ENDIF
          !     
          IF(IFLAG(1)>=32)THEN
             KSEC4(5) = 32
             IFLAG(1)    = IFLAG(1) - 32
          ENDIF
          !     
          IF(IFLAG(1)>=16)THEN
             KSEC4(6) = 16
             IFLAG(1)    = IFLAG(1) - 16
          ENDIF
          !     
          INIL  = IFLAG(1)
       ENDIF
       !     
710    CONTINUE

       !     Print number of unused bits, if required.
       !     
       IF(NDBG==1) WRITE(*,9107) INIL
       !     
       !     Only simple packing of data(real or integer) with no
       !     additional flag field supported.
       !     
       IF(KSEC4(4)/=0.OR.KSEC4(6)/=0)THEN
          KRET = 705
          WRITE(*,9705) KRET
          GO TO 900
       ENDIF
       !     
    ENDIF

    !         Octets 5 - 6 : Scale factor.
    !     One 16 bit field.
    !     
    !     Calculate scale factor, if coding data.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       !     If input data is integer, change it to real.
       !     
       IF(KSEC4(5)==32)THEN
          !CALL RORINT(FIELD,FIELD,ILENF,'R')
          if(realgribkind==8)then
             call rorintr(field,ilenf)
          else
             call rorintrf(field,ilenf)
          endif
       ENDIF
       !     
       !     Change units of data values , if required.
       !     
       IF(KSEC1(23)/=0)THEN
          DO 711 J711 = 1 , ILENF
             FIELD(J711) = FIELD(J711) * real( 10**KSEC1(23),realgribkind)
711       ENDDO
       ENDIF

       !     Find maximum and minimum values in data array, ignoring
       !     any missing-data value. For data in spherical harmonic
       !     format the first word contains the real(0,0) coefficient,
       !     which is treated separately. IREP is 1 for spherical
       !     harmonics, 0 for other data.
       !     
       ILEN = ILENF - IREP
       CALL MAXMIN(FIELD(IREP+1),ILEN,ZMAX,ZMIN,PSEC3(2))

       !     
       IF(NDBG==1) WRITE(*,9106) ZMAX , ZMIN
       !     
       !     Calculate and pack scale factor.
       !     If user has supplied a reference value, use it
       !     unless it exceeds the minimum value. Otherwise
       !     use the minimum value.
       !     
       IF(NFREF==1)THEN
          ZREF = FREF

          !     If integer data being packed, ensure that
          !     reference value represents an integer.
          !     
          IF(KSEC4(5)==32)THEN
             ITEMP(1) = NINT(ZREF)
             ZREF  = real(ITEMP(1),realgribkind)
          ENDIF
          !     
          IF(ZREF<ZMIN)THEN
             WRITE(*,9718) ZREF , ZMIN
             ZREF = ZMIN
          ENDIF
       ELSE
          ZREF = ZMIN
       ENDIF
       !     
       ZS =(ZMAX-ZREF) / real(2**(KSEC4(2)+1)-1,realgribkind)
       IF(abs(ZS)>1.e-14_realgribkind)then
          ZS=LOG(ZS)/LOG(2._realgribkind)+2._realgribkind
       endif
       ISCALE(1) = MIN(INT(ZS),INT(ZS+SIGN(1._realgribkind,ZS)))
       ZSCALE = 2._realgribkind**ISCALE(1)

       !     Set sign bit.
       !     
       IF(ISCALE(1)<0)THEN
          ISCALE(1) = -ISCALE(1)
          ISIGN  = 32768
          ISCALE(1) = ISCALE(1) + ISIGN
       ENDIF

       !     Scale factor has all bits set to 1 for
       !     missing fields.(ECMWF convention only).
       !     
       IF(IMISS==1) ISCALE(1) = 65535
       !     
    ENDIF
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,ISCALE,1,IBITS,16,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 707
       WRITE(*,9707) KRET
       GO TO 900
    ENDIF

    !     If decoding, set scale factor.
    !     
    IF(YFUNC=='D'.or.yfunc=='d')THEN
       ISKALE = ISCALE(1)
       IF(ISKALE>=32768)THEN
          ISCALE(1) = ISCALE(1) - 32768
          ISCALE(1) = - ISCALE(1)
       ENDIF

       KSEC4(33)=ISCALE(1)

       ZSCALE = 2._realgribkind**ISCALE(1)
    ENDIF

    !         Octets 7 - 10 : Reference value.
    !     One 8 bit and one 24 bit field.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       IF(IMISS==1)THEN
          !     
          !     For missing data fields are set to all 1 bits.
          !     
          IEXP(1)   = 255
          IMANT(1)  = 16777215
       ELSE
          !     Convert floating point to GRIB representation.
               
          itrnd = 1
          if(ndbg==1) itrnd = 11
          ztemp = zref
          call confp3(zref,iexp(1),imant(1),ibits,itrnd)
          !     Set reference value to that actually stored
          !     in the GRIB code.
          call decfp2(zref,iexp(1),imant(1))
              
          !     If the nearest number which can be represented in
          !     GRIB format is greater than the reference value,
          !     find the nearest number in GRIB format lower
          !     than the reference value.
               
          if(ztemp<zref.or.ztemp>zref)then
             !     convert floating point to grib representation
             !     using truncation to ensure that the converted
             !     number is smaller than the original one.
             itrnd = itrnd - 1
             zref  = ztemp
             call confp3(zref,iexp(1),imant(1),ibits,itrnd)
             !     set reference value to that actually stored
             !     in the grib code.
             call decfp2(zref,iexp(1),imant(1))
                  
             !if(abs(ztemp-zref)>1.e-16_realgribkind)then
             if(ztemp<zref)then
                WRITE(*,*) 'Reference value error.',ZTEMP,ZREF
                WRITE(*,*) 'Notify J. Hennessy.'
                ZREF = ZTEMP
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    !     Insert / extract fields.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,     &     
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 708
       WRITE(*,9708) KRET
       GO TO 900
    ENDIF
    CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,     &     
         24,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 708
       WRITE(*,9708) KRET
       GO TO 900
    ENDIF

    !     Conversion from GRIB format, if decoding.
    !     
    IF(YFUNC=='D'.or.yfunc=='d')THEN
       !     
       !     Set IMISS to 1 if entire field is missing ie scale
       !     factor, exponent and mantissa with all bits set to 1.
       !     
       IMISS = 0
       IF(ISKALE==65535.AND.IEXP(1)==255.AND.     &        
            IMANT(1)==16777215) IMISS = 1
       !     
       !     Convert GRIB representation to floating point.
       !     
       IF(IMISS==0) THEN
          
          !     Field is present.
          !     
          CALL DECFP2(ZREF,IEXP(1),IMANT(1))
       ELSE
          !     
          !     Field is missing. Print warning message and
          !     field identification sections of Grib code,
          !     forcing field data values to 0.
          !     
          WRITE(*,9719)
          CALL GRPRS0(KSEC0)
          CALL GRPRS1(KSEC0,KSEC1)
          ZREF   = 0._realgribkind
          ZSCALE = 0._realgribkind
       ENDIF
       IF(NDBG==1) WRITE(*,9105) ZREF
    ENDIF

    !         Octet 11 : Number of bits containing each packed value.
    !     One 8 bit field.
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,INSPT,KSEC4(2),1,IBITS,     &     
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 709
       WRITE(*,9709) KRET
       GO TO 900
    ENDIF

    !     If decoding , calculate number of data values.
    !     
    IF(YFUNC=='D'.or.yfunc=='d')THEN

       !     The following allows zero bits per datum(for constant fields)
       IF(KSEC4(2)==0)THEN
          ILENF=KSEC2(2)*KSEC2(3)
       ELSE

          ILENF = ILEN4(1) - 11 - IREP * 4
          ILENF =(ILENF * 8 - INIL ) / KSEC4(2) + IREP

       ENDIF

       KSEC4(1) = ILENF

       !     Check length of output array.
       !     
       IF(ILENF>FIELDSIZE)THEN
          KRET = 710
          WRITE(*,9710) FIELDSIZE , ILENF , KRET
          GO TO 900
       ENDIF
       !     
    ENDIF

    !         Octets 12 - 15 : If spherical harmonic data, real(0,0)
    !     coefficient is in floating point representation.
    !     One 8 bit and one 24 bit field.
    !     
    IF(IREP==1) THEN
       ITRND = 1
       IF(YFUNC=='C'.or.yfunc=='c')then

          !     Convert floating point to GRIB representation.

          CALL CONFP3(FIELD(1),IEXP(1),IMANT(1),IBITS,ITRND)
       endif
       !     
       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IEXP,1,IBITS,     &        
            8,YFUNC,KRET)
       IF(KRET/=0)THEN
          KRET = 711
          WRITE(*,9711) KRET
          GO TO 900
       ENDIF
       CALL INXBIT(KGRIB,KLENG,INSPT,IMANT,1,IBITS,     &        
            24,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 711
          WRITE(*,9711) KRET
          GO TO 900
       ENDIF
       !     
       IF(YFUNC=='D'.or.yfunc=='d')THEN

          !     Convert GRIB representation to floating point.
          !     
          IF(IMISS==1)THEN
             FIELD(1) = 0.0_realgribkind
          ELSE
             CALL DECFP2(FIELD(1),IEXP(1),IMANT(1))
          ENDIF
       ENDIF
    ENDIF

    !         Scale and store, or extract and scale data values.
    !     
    !     Only a few points to be unpacked.
    !     
    IF(HOPER=='X'.or.HOPER=='x')THEN
       !     
       !     Check that no bit-map is included.
       !     
       IF(KSEC1(5)==64.OR.KSEC1(5)==192)THEN
          KRET = 717
          WRITE(*,9717) KRET
          GO TO 900
       ENDIF

       !     Check that field is Gaussian or latitude/longitude grid.
       !     
       IF(KSEC2(1)/=0.AND.KSEC2(1)/=4.AND.   &
            KSEC2(1)/=10.AND.KSEC2(1)/=14.AND. &
            KSEC2(1)/=20.AND.KSEC2(1)/=24.AND. &
            KSEC2(1)/=30.AND.KSEC2(1)/=34) THEN
          KRET = 716
          WRITE(*,9716) KRET
          GO TO 900
       ENDIF

       !     Check that scanning mode is West to East and North to
       !     South.
       !     
       IF(KSEC2(11)/=0) THEN
          KRET = 715
          WRITE(*,9715) KRET
          GO TO 900
       ENDIF

       !     Check that number of points required does not exceed
       !     maximum or minimum allowed.
       !     
       IF(KSEC4(34)>4.OR.KSEC4(34)<1) THEN
          KRET = 714
          WRITE(*,9714) KSEC4(34) , KRET
          GO TO 900
       ENDIF
       !     
       ITEMP(1) = 1
       !     
       DO 713 J713=1,KSEC4(34)

          !     Skip down latitude rows.
          !     
          IF(KSEC2(17)==0) THEN
             !     
             !     Regular grid.
             !     
             ISKIP =(KSEC4(34+ITEMP(1))-1) * KSEC2(2)
          ELSE
             !     
             !     Quasi-regular grid.
             !     
             ISKIP = 0
             DO 712 J712=1,KSEC4(34+ITEMP(1)) - 1
                ISKIP = ISKIP + KSEC2(22+J712)
712          ENDDO
          ENDIF

          !     Skip any points not required on this latitude row.
          !     
          ISKIP = ISKIP + KSEC4(34+ITEMP(1)+1) - 1
          !     
          !     Calculate number of bits in these values and add
          !     to current value of bit-pointer.
          !     
          ISKIP = ISKIP * KSEC4(2) + INSPT
          !     
          !     Extract value from 1 point.
          !convert field to integer-bitstream
          if(realgribkind==8)then
             call r2ibts(field(j713),mytmp)
          else
             call r2ibtsf(field(j713),mytmp)
          endif

          CALL INXBIT(KGRIB,KLENG,ISKIP,mytmp,1, &
               IBITS,KSEC4(2),YFUNC,KRET)
          !deconvert field from integer-bitstream to real
          if(realgribkind==8)then
             call i2rbts(field(j713),mytmp)
          else
             call i2rbtsf(field(j713),mytmp)
          endif

          IF(KRET/=0) THEN
             KRET = 712
             WRITE(*,9712) KRET
             GO TO 900
          ENDIF
          ITEMP(1) = ITEMP(1) + 2
713    ENDDO
       !     
       ILENF    = KSEC4(34)
       ILENFM   = KSEC4(34)
       KSEC4(1) = KSEC4(34)
       GO TO 715
    ENDIF

    !     All data to be unpacked or packed.
    !     
    ILENFM = ILENF - IREP
    !     Avoid overwriting FIELD by strip-mining into array ITMPD
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       IST=IREP+1
702    CONTINUE
       IHE=MIN(JPS4,ILENFM)
       CALL INSCAL(FIELD(IST),ITMPD,IHE,ZREF,ZSCALE)

       !     Insert / extract fields.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,ITMPD,IHE,IBITS, &
            KSEC4(2),YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 712
          WRITE(*,9712) KRET
          GO TO 900
       ENDIF
       IST=IST+IHE
       ILENFM=ILENFM-IHE
       IF(ILENFM>0)GOTO 702
       ILENFM=ILENF-IREP
    ELSE

       !     Insert / extract fields.
       allocate(mytmpv(ilenfm))
       if(realgribkind==8)then
          call r2ibtsv(field(irep+1),mytmpv,ilenfm)
       else
          call r2ibtsvf(field(irep+1),mytmpv,ilenfm)
       endif

       CALL INXBIT(KGRIB,KLENG,INSPT,MYTMPV,ILENFM,IBITS, &
            KSEC4(2),YFUNC,KRET)
       if(realgribkind==8)then
          call i2rbtsv(field(irep+1),mytmpv,ilenfm)
       else
          call i2rbtsvf(field(irep+1),mytmpv,ilenfm)
       endif
       deallocate(mytmpv)
       IF(KRET/=0)THEN
          KRET = 712
          WRITE(*,9712) KRET
          GO TO 900
       ENDIF
    ENDIF
    !     
715 CONTINUE
    !     
    IF(YFUNC=='D'.or.yfunc=='d')THEN
       !     
       !         do j716=1,ilenfm
       !            field(irep+1+j716) = zref + field(irep+1+j716)*zscale 
       !         enddo
       !         if(mype==0)then
       !            do j716=1,ilenfm
       !               print *, field(irep+j716),'gribex'
       !            enddo
       !         endif
       CALL EXSCAL(FIELD(IREP+1),ILENFM,ZREF,ZSCALE)



       !     Change units of data values, if required.
       !     
       IF(KSEC1(23)/=0)THEN
          DO 716 J716 = 1 , ILENF
             FIELD(J716) = FIELD(J716) / real(10**KSEC1(23),realgribkind)
716       ENDDO
       ENDIF

       !     Convert to integer if original data was integer.
       !     
       IF(KSEC4(5)==32)THEN
          !CALL RORINT(FIELD,FIELD,ILENF,'I')
          if(realgribkind==8)then
             call rorinti(field,ilenf)
          else
             call rorintif(field,ilenf)
          endif
       ENDIF
       !     
       !     Finish, if only a few points extracted.
       !     
       IF(HOPER=='X'.or.HOPER=='x') GO TO 900
    ENDIF
    !     
    !         Enter length of binary data section, ensuring that the
    !     length is an even number of octets, padding with binary
    !     zeroes as required.
    !     One 24 bit field.
    !     
    IF(YFUNC/='C'.and.yfunc/='c') GO TO 799
    !     
    !     Length of section 4, in bits.
    !     
    ILEN4(1) = INSPT - IPLEN
    IL    = ILEN4(1) / 16
    IL    = ILEN4(1) -( IL * 16 )
    INIL  = 0
    IF(IL/=0) THEN
       INIL = 16 - IL
       INSPT = INSPT + INIL
       ILEN4(1) = ILEN4(1) + INIL
    ENDIF
    !     
    ILEN4(1) = ILEN4(1) / 8

    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,IPLEN,ILEN4,1,IBITS, &
         24,YFUNC,KRET)
    IF(KRET/=0) THEN
       KRET = 701
       WRITE(*,9701) KRET
       GO TO 900
    ENDIF

    !         Enter flag / unused bits field.
    !     One 8 bit field.
    !     Two 4 bit fields.
    !     
    IFLAG(1) = IFLAG(1) + INIL
    !     
    !     Print number of unused bits, if required.
    !     
    IF(NDBG==1) WRITE(*,9107) INIL
    !     
    !     Insert / extract field.
    !     
    CALL INXBIT(KGRIB,KLENG,IPLEN,IFLAG,1,IBITS, &
         8,YFUNC,KRET)
    IF(KRET/=0)THEN
       KRET = 713
       WRITE(*,9713) KRET
       GO TO 900
    ENDIF

    !         Check consistency of values given, with GRIB code, if required.
    !     
799 CONTINUE
    !     
    IF(NVCK==1.AND.(YFUNC=='C'.or.yfunc=='c'))THEN
       CALL GRCHK4(KSEC4,KRET)
       IF(KRET/=0) THEN
          KRET = 799
          WRITE(*,9799) KRET
          GO TO 900
       ENDIF
    ENDIF
    !     
    !     ----------------------------------------------------------------

    !         Section 8 . Code/decode End Section(Section 5) of GRIB code.
    !     ----------------------------------------------------------------
    !     
800 CONTINUE
    !     
    IF(NDBG==1) WRITE(*,*) ' GRIBEX : Section 8.'
    !     
    !         Ascii 7 7 7 7 at end of coded data.
    !     Four 8 bit fields.
    !     
    IF(YFUNC=='C'.or.yfunc=='c')THEN
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,I7777(1),4,IBITS, &
            8,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 801
          WRITE(*,9801) KRET
          GO TO 900
       ENDIF

       !                Length of GRIB message.
       !     
       KSEC0(1) = INSPT / 8
       ITEMP(1) = 32
       !     
       !     Insert / extract field.
       !     
       CALL INXBIT(KGRIB,KLENG,ITEMP(1),KSEC0(1),1,IBITS, &
            24,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 802
          WRITE(*,9802) KRET
          GO TO 900
       ENDIF
       !     
    ELSE

       !     Skip padding.
       !     
       INSPT = INSPT + INIL
       !     
       CALL INXBIT(KGRIB,KLENG,INSPT,IPARM(1),4,IBITS, &
            8,YFUNC,KRET)
       IF(KRET/=0) THEN
          KRET = 801
          WRITE(*,9801) KRET
          GO TO 900
       ENDIF

       !     Check that 7777 group is found where expected.
       !     
       ICOUNT = 0
       DO 802 J802 = 1 , 4
          IF(IPARM(J802)/=55) ICOUNT = ICOUNT + 1
802    ENDDO
       IF(ICOUNT/=0)THEN
          KRET = 805
          WRITE(*,9805) KRET
       ENDIF
    ENDIF

    !                Final handling when bit-map section included.
    !     
    IF(KSEC1(5)==64.OR.KSEC1(5)==192)THEN
       !     
       IF(KSEC3(1)==0)THEN
          !     
          !     Bit-map included in GRIB message.
          !     
          !     Set integer or real missing data value.
          !     ZMSVAL and IMSVAL are euqivalenced.
          !     
          IF(KSEC4(5)==0)THEN
             ZMSVAL = PSEC3(2)
          ELSE
             IMSVAL = KSEC3(2)
          ENDIF

          !     Set number of values decoded to include
          !     missing data values.
          !     
          KSEC4(1) = IVALS
          !     
          CALL INXMAP(KGRIB,KLENG,IBMAP,FIELD,IVALS, &
               IBITS,ISBMAP,ZMSVAL,'D',NDBG,KRET)

          IF(KRET/=0)THEN
             KRET = 604
             WRITE(*,9604) KRET
             GO TO 900
          ENDIF
       ELSE

          !     Predetermined bit-map reference only.
          !     
          ISBMAP = -3
          !     
       ENDIF
    ENDIF
    !     
    !                If required, convert quasi-regular Gaussian grid to
    !     regular.
    !     
    IF((HOPER=='R'.or.HOPER=='r') &
         .AND.KSEC2(1)==4.AND.KSEC2(17)==1)THEN
       !     
       !     Cannot handle data with bit-map.
       !     
       IF(KSEC1(5)==64.OR.KSEC1(5)==192)THEN
          KRET = 605
          WRITE(*,9605) KRET
          GO TO 900
       ENDIF
       !     
       INOLAT    = KSEC2(10) * 2
       INOLNG    = KSEC2(10) * 4
       CALL QU2REG(FIELD,KSEC2(23),INOLAT,INOLNG,1)
       KSEC4(1)  = INOLAT * INOLNG
       KSEC2(2)  = 4 * KSEC2(10)
       KSEC2(17) = 0
    ENDIF

    !     Set number of values decoded negative, if missing data.
    !     
    IF(IMISS==1) KSEC4(1) = - KSEC4(1)
    !     
    !     If GRIB Edition number is -1 or 0, set GRIB message
    !     length for return to user.
    !     
    IF(KSEC0(2)==-1.OR.KSEC0(2)==0) &
         KSEC0(1) = INSPT / 8
    IF(YFUNC/='C'.and.yfunc/='c')GO TO 900


    !         Any unused part of last word is already set to binary zeroes.
    !     Increment pointers as necessary.
    !     
    KWORD    = INSPT / IBITS
    ITEMP(1)    = KWORD * IBITS
    IOFF     = INSPT - ITEMP(1)
    IF(IOFF/=0) THEN
       INSPT = INSPT + IBITS - IOFF
       KWORD = KWORD + 1
    ENDIF

    !         Round length to a multiple of 120 octets, if required,
    !     any additional words are already set to 0.
    !     
    IF(NRND==1) THEN
       I = INSPT / 960
       I = I * 960
       I = INSPT - I
       IF(I/=0) I =(960 - I) / IBITS
       KWORD = KWORD + I
    ENDIF
    !     
    !     ----------------------------------------------------------------

    !         Section 9 . Abort/return to calling routine. Format statements.
    !     ----------------------------------------------------------------
    !     
900 CONTINUE
    !     
    IF(NDBG==1) THEN
       WRITE(*,*) ' GRIBEX : Section 9.'
       WRITE(*,*) '          Output values set -'
       IF(YFUNC=='D'.or.yfunc=='d')THEN
          CALL GRPRS0(KSEC0)
          CALL GRPRS1(KSEC0,KSEC1)
          CALL GRPRS2(KSEC0,KSEC2,PSEC2)
          CALL GRPRS3(KSEC3,PSEC3)
          CALL GRPRS4(KSEC4,FIELD)
       ENDIF
    ENDIF

    !         If no error has been encountered, set return code to informative
    !     value, if required.
    !     
    !     Set pseudo-GRIB data encountered.
    !     
    IF(KRET==0.AND.IPSEUD/=0) KRET = IPSEUD
    !     
    !     Set data with bit-map encountered.
    !     
    IF(KRET==0.AND.ISBMAP/=0) KRET = ISBMAP
    !     
    !         Abort if an error has been encountered and user has requested
    !     an abort. Informative values are negative and do not cause an
    !     abort.
    !     
    IF(IRET==0.AND.KRET>0)THEN
       CALL STOP_PROGRAM('')
    ELSE
       RETURN
    ENDIF

9101 FORMAT(1H ,'GRIBEX : Length of Section 1 of GRIB code is ',I4,' octets.')
9102 FORMAT(1H ,'GRIBEX : Length of Section 2 of GRIB code is ',I4,' octets.')
9103 FORMAT(1H ,'GRIBEX : Length of Section 3 of GRIB code is ',I8,' octets.')
9104 FORMAT(1H ,'GRIBEX : Length of Section 4 of GRIB code is ',I8,' octets.')
9105 FORMAT(1H ,'GRIBEX : Decoded reference value = ',F30.20)
9106 FORMAT(1H ,'GRIBEX : Maximum and minimum values = ',F30.20,&
         2X,F30.20)
9107 FORMAT(1H ,'GRIBEX : Number of unused bits is ',I2,'.')
9201 FORMAT(1H ,'GRIBEX : Invalid function requested - ',A,     & 
         '. Return code = ',I3)
9202 FORMAT(1H ,'GRIBEX : Number of bits per data value, ',     &
         I3,' exceeds word length ',I3,'. Return code = ',I3)
9203 FORMAT(1H ,'GRIBEX : Non-zero value in missing data field',     &
         '. Return code = ',I3)
9301 FORMAT(1H ,'GRIBEX : Error in inserting letters GRIB', &
         '. Return code = ',I3)
9302 FORMAT(1H ,'GRIBEX : Error extracting length of GRIB message', &
         '. Return code = ',I3)
9303 FORMAT(1H ,'GRIBEX : Error inserting/extracting GRIB Edition',&
         ' Number. Return code = ',I3)
9304 FORMAT(1H ,'GRIBEX : Error extracting octets 22 and 23 for ',&
         'Experimental Edition check. Return code = ',I3)
9401 FORMAT(1H ,'GRIBEX : Error inserting/extracting length of ',&
         ' Section 1. Return code = ',I3)
9402 FORMAT(1H ,'GRIBEX : Error inserting/extracting'// &
         ' Parameter Table Version Number. Return code = ',I3)
9403 FORMAT(1H ,'GRIBEX : Error inserting/extracting six fields ',&
         'from Identification of Centre to'//                      &
         ' Indicator of type of level.',                           &
         ' Return code = ',I3)
9404 FORMAT(1H ,'GRIBEX : Error inserting/extracting Height, ',   &
         ' pressure, etc of levels. Return code = ',I3)
9405 FORMAT(1H ,'GRIBEX : Error inserting/extracting six fields ',  &
         'from Year of century to Indicator of unit of time range.', &
         ' Return code = ',I3)
9406 FORMAT(1H ,'GRIBEX : Error inserting/extracting Period of time.', &
         ' Return code = ',I3)
9407 FORMAT(1H ,'GRIBEX : Error inserting/extracting time range ', &
         'indicator. Return code = ',I3)
9408 FORMAT(1H ,'GRIBEX : Error inserting/extracting number ', &
         'averaged. Return code = ',I3)
9409 FORMAT(1H ,'GRIBEX : Error inserting/extracting number ', &
         'missing from averages etc. Return code = ',I3)
9410 FORMAT(1H ,'GRIBEX : Error inserting/extracting century of ', &
         'data or reserved field. Return code = ',I3)
9411 FORMAT(1H ,'GRIBEX : Error inserting/extracting units ', &
         'decimal scale factor. Return code = ',I3)
9499 FORMAT(1H ,'GRIBEX : Error found when checking values for ', &
         'Section 1 against valid GRIB values. Return code = ',I3)
9501 FORMAT(1H ,'GRIBEX : Error inserting/extracting length of ', &
         ' Section 2. Return code = ',I3)
9502 FORMAT(1H ,'GRIBEX : Error inserting/extracting number of ', &
         ' Vertical coordinate parameters. Return code = ',I3)
9503 FORMAT(1H ,'GRIBEX : Error inserting/extracting location of ',&
         'List of vertical coordinate parameters',/, &
         '         or List of numbers of points. Return code = ',I3)
9504 FORMAT(1H ,'GRIBEX : Error inserting/extracting data ', &
         'representation type. Return code = ',I3)
9505 FORMAT(1H ,'GRIBEX : Error inserting/extracting number of ', &
         'points along a parallel or meridian. Return code = ',I3)
9506 FORMAT(1H ,'GRIBEX : Error inserting/extracting latitude or ', &
         'longitude of first grid point. Return code = ',I3)
9507 FORMAT(1H ,'GRIBEX : Error inserting/extracting components', &
         ' flag. Return code = ',I3)
9508 FORMAT(1H ,'GRIBEX : Error inserting/extracting latitude or ', &
         'longitude of last grid point. Return code = ',I3)
9509 FORMAT(1H ,'GRIBEX : Error inserting/extracting i direction', &
         ' increment. Return code = ',I3)
9510 FORMAT(1H ,'GRIBEX : Error inserting/extracting number of ', &
         'parallels between pole and Equator. Return code = ',I3)
9511 FORMAT(1H ,'GRIBEX : Error inserting/extracting scanning ', &
         'mode flags. Return code = ',I3)
9513 FORMAT(1H ,'GRIBEX : Error inserting/extracting j direction', &
         ' increment. Return code = ',I3)
9514 FORMAT(1H ,'GRIBEX : Error inserting/extracting J,K,M ', &
         'pentagonal resolution parameters. Return code = ',I3)
9515 FORMAT(1H ,'GRIBEX : Error inserting/extracting representation', &
         ' type or mode. Return code = ',I3)
9517 FORMAT(1H ,'GRIBEX : Error inserting/extracting latitude or ', &
         'longitude of southern pole. Return code = ',I3)
9518 FORMAT(1H ,'GRIBEX : Error inserting/extracting angle ', &
         'of rotation. Return code = ',I3)
9519 FORMAT(1H ,'GRIBEX : Error inserting/extracting latitude or ', &
         'of pole of stretching. Return code = ',I3)
9520 FORMAT(1H ,'GRIBEX : Error inserting/extracting ', &
         'stretching factor. Return code = ',I3)
9521 FORMAT(1H ,'GRIBEX : Error inserting/extracting ', &
         'vertical coordinate parameters. Return code = ',I3)
9522 FORMAT(1H ,'GRIBEX : Error inserting/extracting list of ', &
         'numbers of points. Return code = ',I3)
9523 FORMAT(1H ,'GRIBEX : Error inserting/extracting number of ',&
         'points along X or Y axis. Return code = ',I3)
9524 FORMAT(1H ,'GRIBEX : Error inserting/extracting X or Y axis',&
         ' grid length. Return code = ',I3)
9525 FORMAT(1H ,'GRIBEX : Error inserting/extracting Projection',&
         ' centre flag. Return code = ',I3)
9526 FORMAT(1H ,'GRIBEX : Error inserting/extracting latitude',&
         ' or longitude of sub-satellite point. Return code = ',I3)
9527 FORMAT(1H ,'GRIBEX : Error inserting/extracting diameter',&
         ' of the earth in x or y direction. Return code = ',I3)
9528 FORMAT(1H ,'GRIBEX : Error inserting/extracting X or Y',&
         ' coordinate of sub-satellite point. Return code = ',I3)
9529 FORMAT(1H ,'GRIBEX : Error inserting/extracting orientation ', &
         'of the grid or camera angle. Return code = ',I3)
9598 FORMAT(1H ,'GRIBEX : Representation type not catered for - ',I3, &
         '. Return code = ',I3)
9599 FORMAT(1H ,'GRIBEX : Error found when checking values for ',&
         'Section 2 against valid GRIB values. Return code = ',I3)
9601 FORMAT(1H ,'GRIBEX : Error inserting/extracting length of ',&
         ' Section 3. Return code = ',I3)
9602 FORMAT(1H ,'GRIBEX : Error inserting/extracting number',&
         ' of unused bits at the end of Section 3. Return code = ',I3)
9603 FORMAT(1H ,'GRIBEX : Error inserting/extracting bit-map',&
         ' reference table. Return code = ',I3)
9604 FORMAT(1H ,'GRIBEX : Error inserting/extracting bit-map.',&
         ' Return code = ',I3)
9605 FORMAT(1H ,'GRIBEX : Cannot convert Quasi-regular Gaussian',&
         ' grid with bit-map. Return code = ',I3)
9606 FORMAT(1H ,'GRIBEX : Bit-map found. No data decoded.',&
         ' Return code = ',I3)
9699 FORMAT(1H ,'GRIBEX : Error found when checking values for ',&
         'Section 3 against valid GRIB values. Return code = ',I3)
9701 FORMAT(1H ,'GRIBEX : Error inserting/extracting length of ', &
         ' Section 4. Return code = ',I3)
9705 FORMAT(1H ,'GRIBEX : Only simple packing of data(real or ',&
         'integer) with no additional flag allowed. Return code = ',&
         I3)
9706 FORMAT(1H ,'GRIBEX : Error in extracting section 4 flag field.',&
         ' Return code = ',I3)
9707 FORMAT(1H ,'GRIBEX : Error inserting/extracting scale factor.',&
         ' Return code = ',I3)
9708 FORMAT(1H ,'GRIBEX : Error inserting/extracting reference',&
         ' value. Return code = ',I3)
9709 FORMAT(1H ,'GRIBEX : Error inserting/extracting number of bits',&
         ' per data value. Return code = ',I3)
9710 FORMAT(1H ,'GRIBEX : Output array too small. Length is ',I8,&
         ' words. Number of values is ',I8,' Return code = ',I3)
9711 FORMAT(1H ,'GRIBEX : Error inserting/extracting real ',&
         'coefficient. Return code = ',I3)
9712 FORMAT(1H ,'GRIBEX : Error inserting/extracting data ',&
         'values. Return code = ',I3)
9713 FORMAT(1H ,'GRIBEX : Error inserting/extracting flag ',&
         'and unused bit field. Return code = ',I3)
9714 FORMAT(1H ,'GRIBEX : Function is X and ',&
         'number of values is ',I3,'. Return code = ',I3)
9715 FORMAT(1H ,'GRIBEX : Function is X and scanning mode is not',&
         'North to South and West to East. Return code = ',I3)
9716 FORMAT(1H ,'GRIBEX : Function is X and field is not ',&
         'Gaussian or Latitude/longitude field. Return code = ',I3)
9717 FORMAT(1H ,'GRIBEX : Function is X and a bit-map is included.',&
         ' Return code = ',I3)
9718 FORMAT(1H ,'GRIBEX : User supplied reference value ',F20.10,&
         ' exceeds minimum value ',F20.10,'. Minimum value used.')
9719 FORMAT(1H ,'GRIBEX : Following field is missing ')
9799 FORMAT(1H ,'GRIBEX : Error found when checking values for ',&
         'Section 4 against valid GRIB values. Return code = ',I3)
9801 FORMAT(1H ,'GRIBEX : Error inserting/extracting 7777 group.',&
         ' Return code = ',I3)
9802 FORMAT(1H ,'GRIBEX : Error inserting/extracting length of GRIB',&
         ' message. Return code = ',I3)
9805 FORMAT(1H ,'GRIBEX : End of message 7777 group not found.',&
         ' Return code = ',I3)

  END subroutine gribex
  subroutine grprs0(ksec0)
    implicit none
    !
    !**** grprs0 - print information from section 0 of grib code.
    !
    !     purpose.
    !     --------
    !
    !           print the information in the indicator
    !           section(section 0) of decoded grib data.
    !
    !**   interface.
    !     ----------
    !
    !           call grprs0(ksec0)
    !
    !               input parameters.
    !               -----------------
    !
    !               ksec0 - array of decoded parameters from section 0.
    !
    !               output parameters.
    !               ------------------
    !
    !               none.
    !
    !     method.
    !     -------
    !
    !           fields are printed as integers.
    !
    !     externals.
    !     ----------
    !
    !           none.
    !
    !     reference.
    !     ----------
    !
    !           wmo manual on codes re grib code.
    !           see also routine gribex.
    !
    !     comments.
    !     ---------
    !
    !           routine contains sections 0 to 1 and section 9.
    !
    !     author.
    !     -------
    !
    !           j. hennessy    ecmwf 18.06.91
    !
    !     modifications.
    !     --------------
    !
    !           j. hennessy    ecmwf 30.08.91
    !           changes to some comments only.
    !
    !     -----------------------------------------------------------------
    !
    integer:: ksec0(*)

    write(*,9000)
    write(*,9001)
    write(*,9002)
    write(*,9003) ksec0(1)
    write(*,9004) ksec0(2)
9000 format(1h )
9001 format(1h ,'section 0 - indicator section.       ')
9002 format(1h ,'-------------------------------------')
9003 format(1h ,'length of grib message(octets).     ',i9)
9004 format(1h ,'grib edition number.                 ',i9)

    return

  end subroutine grprs0
  subroutine grprs1(ksec0,ksec1)
    implicit none
    !     
    !**** grprs1 - print information from section 1 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     print the information in the product definition
    !     section(section 1) of decoded grib data.
    !     
    !**   interface.
    !     ----------
    !     
    !     call grprs1(ksec0,ksec1)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec0 - array of decoded parameters from section 0.
    !     
    !     ksec1 - array of decoded parameters from section 1.
    !     
    !     output parameters.
    !     ------------------
    !     
    !     none.
    !     
    !     method.
    !     -------
    !     
    !     flag fields are printed in binary representation.
    !     other fields as integers.
    !     fields printed depend on grib edition.
    !     
    !     externals.
    !     ----------
    !     
    !     prtbin
    !     
    !     reference.
    !     ----------
    !     
    !     wmo manual on codes re grib code.
    !     see also routine gribex.
    !     
    !     comments.
    !     ---------
    !     
    !     when decoding data from experimental edition or edition 0,
    !     routine gribex adds the additional fields available in
    !     edition 1.
    !     
    !     routine contains sections 0 to 1 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     j. hennessy    ecmwf 18.06.91
    !     
    !     modifications.
    !     --------------
    !     
    !     j. hennessy    ecmwf 30.08.91
    !     different print for grib editions up to edition 1
    !     removed.
    !     
    !     j. hennessy    ecmwf 07.01.92
    !     different print for grib editions up to edition 1
    !     added for centres other than ecmwf.
    integer:: ibit
    integer:: ierr
    integer:: iout

    integer:: ksec0(*),ksec1(*)

    write(*,9000)
    write(*,9001)
    write(*,9002)

    write(*,9100) ksec1(1)
    write(*,9101) ksec1(2)
    write(*,9102) ksec1(3)
    write(*,9103) ksec1(4)
    ibit = 8
    call prtbin(ksec1(5),ibit,iout,ierr)
    write(*,9104) iout
    write(*,9105) ksec1(6)
    write(*,9106) ksec1(7)
    write(*,9107) ksec1(8)
    write(*,9108) ksec1(9)
    write(*,9109) ksec1(10)
    write(*,9110) ksec1(11)
    write(*,9111) ksec1(12)
    write(*,9112) ksec1(13)
    write(*,9113) ksec1(14)
    write(*,9114) ksec1(15)
    write(*,9115) ksec1(16)
    write(*,9116) ksec1(17)
    write(*,9117) ksec1(18)
    write(*,9118) ksec1(19)
    write(*,9119) ksec1(20)
    !     all ecmwf data in grib editions before edition 1 is decoded
    !     as 20th century data. other centres are decoded as missing.

    if(ksec0(2)<1.and.ksec1(2)/=98)then
       write(*,9122)
    else
       write(*,9120) ksec1(21)
    endif

    write(*,9121) ksec1(23)


9000 format(1h )
9001 format(1h ,'section 1 - product definition section.')
9002 format(1h ,'---------------------------------------')
9100 format(1h ,'code table 2 version number.         ',i9)
9101 format(1h ,'originating centre identifier.       ',i9)
9102 format(1h ,'model identification.                ',i9)
9103 format(1h ,'grid definition.                     ',i9)
9104 format(1h ,'flag(code table 1)                   ',i8.8)
9105 format(1h ,'parameter identifier(code table 2). ',i9)
9106 format(1h ,'type of level(code table 3).        ',i9)
9107 format(1h ,'value 1 of level(code table 3).     ',i9)
9108 format(1h ,'value 2 of level(code table 3).     ',i9)
9109 format(1h ,'year of data.                        ',i9)
9110 format(1h ,'month of data.                       ',i9)
9111 format(1h ,'day of data.                         ',i9)
9112 format(1h ,'hour of data.                        ',i9)
9113 format(1h ,'minute of data.                      ',i9)
9114 format(1h ,'time unit(code table 4).            ',i9)
9115 format(1h ,'time range one.                      ',i9)
9116 format(1h ,'time range two.                      ',i9)
9117 format(1h ,'time range indicator(code table 5)  ',i9)
9118 format(1h ,'number averaged.                     ',i9)
9119 format(1h ,'number missing from average.         ',i9)
9120 format(1h ,'century of data.                     ',i9)
9121 format(1h ,'units decimal scaling factor.        ',i9)
9122 format(1h ,'century of data.                     not given')
    !     
    return
    !     
  end subroutine grprs1
  subroutine grprs2(ksec0,ksec2,psec2)
    implicit none
    !     
    !**** grprs2 - print information from section 2 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     print the information in the grid description
    !     section(section 2) of decoded grib data.
    !     
    !**   interface.
    !     ----------
    !     
    !     call grprs2(ksec0,ksec2,psec2)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec0 - array of decoded integers from section 0.
    !     
    !     ksec2 - array of decoded integers from section 2.
    !     
    !     psec2 - array of decoded reals from section 2.
    !     
    !     output parameters.
    !     ------------------
    !     
    !     none.
    !     
    !     method.
    !     -------
    !     
    !     flag fields are printed in binary representation
    !     other fields as integers or reals, as appropriate.
    !     fields printed depend on grib edition.
    !     
    character(len=10):: yout

    integer:: i
    integer:: ibit
    integer:: iedit
    integer:: ierr
    integer:: iout
    integer:: ij
    integer:: ik
    integer:: iresol

    integer:: j103

    integer:: ksec0(*)
    integer:: ksec2(*)

    real(kind=realgribkind) psec2(*)


    !     grib edition number.

    iedit = ksec0(2)

    write(*,9000)
    write(*,9001)
    write(*,9002)
    !     *    spherical harmonic data.

    if(ksec2(1)==50.or.ksec2(1)==60.or.ksec2(1)==70 &
         .or.ksec2(1)==80)then
       write(*,9101) ksec2(1)
       write(*,9102) ksec2(2)
       write(*,9103) ksec2(3)
       write(*,9104) ksec2(4)
       write(*,9105) ksec2(5)
       write(*,9106) ksec2(6)
       write(*,9107)(ksec2(i),i=7,11)
       write(*,9212) ksec2(12)
       go to 190
    endif
    !     *    gaussian grid.

    if(ksec2(1)==4.or.ksec2(1)==14.or.ksec2(1)==24 &
         .or.ksec2(1)==34)then
       write(*,9200)
       write(*,9101) ksec2(1)
       !     quasi-regular grids introduced in edition 1.

       if(ksec2(17)==0.or.iedit<1)then
          write(*,9201) ksec2(2)
       else
          write(*,9233)
          write(*,9234)
          ij = 0
          yout = ' '
          do 103 j103 =1,ksec2(3)
             ij = ij + 1
             write(yout(1:3),'(i3)') ij
             if(ij>ksec2(3)) go to 104
             if(ij==ksec2(3))then
                write(*,9235) ksec2(ij+22) , yout
                go to 104
             endif
             ik = 0
102          continue
             ik = ik + 1
             if(ksec2(ij+22+1)==ksec2(ij+22))then
                ij = ij + 1
                go to 102
             endif
             if(ik>1)then
                yout(4:) = ' to '
                write(yout(8:10),'(i3)') ij
             endif
             write(*,9235) ksec2(ij+22) , yout
             yout = ' '
103       enddo
104       continue
       endif
       write(*,9202) ksec2(3)
       write(*,9203) ksec2(4)
       write(*,9204) ksec2(5)
       ibit = 8
       iresol = ksec2(6) + ksec2(18) + ksec2(19)
       call prtbin(iresol,ibit,iout,ierr)
       write(*,9205) iout
       write(*,9206) ksec2(7)
       write(*,9207) ksec2(8)

       !     print increment if given.

       if(ksec2(6)==128)then
          write(*,9208) ksec2(9)
       else
          write(*,9236)
       endif
       write(*,9210) ksec2(10)
       ibit = 8
       call prtbin(ksec2(11),ibit,iout,ierr)
       write(*,9211) iout
       write(*,9212) ksec2(12)
       go to 190

    endif
    !     *    latitude / longitude grids.

    if(ksec2(1)==0.or.ksec2(1)==10.or.ksec2(1)==20     &     
         .or.ksec2(1)==30)then
       write(*,9200)
       write(*,9101) ksec2(1)
       write(*,9201) ksec2(2)
       write(*,9202) ksec2(3)
       write(*,9203) ksec2(4)
       write(*,9204) ksec2(5)
       ibit = 8
       iresol = ksec2(6) + ksec2(18) + ksec2(19)
       call prtbin(iresol,ibit,iout,ierr)
       write(*,9205) iout
       write(*,9206) ksec2(7)
       write(*,9207) ksec2(8)
       !     print increment if given.

       if(ksec2(6)==128)then
          write(*,9208) ksec2(9)
          write(*,9209) ksec2(10)
       else
          write(*,9236)
          write(*,9237)
       endif
       ibit = 8
       call prtbin(ksec2(11),ibit,iout,ierr)
       write(*,9211) iout
       write(*,9212) ksec2(12)
       go to 190
    endif
    !     *    polar stereographic.

    if(ksec2(1)==5)then
       write(*,9200)
       write(*,9101) ksec2(1)
       write(*,9301) ksec2(2)
       write(*,9302) ksec2(3)
       write(*,9203) ksec2(4)
       write(*,9204) ksec2(5)
       ibit = 8
       iresol = ksec2(18) + ksec2(19)
       call prtbin(iresol,ibit,iout,ierr)
       write(*,9205) iout
       write(*,9303) ksec2(7)
       write(*,9304) ksec2(9)
       write(*,9305) ksec2(10)
       ibit = 8
       call prtbin(ksec2(11),ibit,iout,ierr)
       write(*,9211) iout
       write(*,9212) ksec2(12)
       write(*,9306) ksec2(13)
       go to 190
    endif
    !     *    space view perspective or orthographic.

    if(ksec2(1)==90)then
       write(*,9200)
       write(*,9101) ksec2(1)
       write(*,9301) ksec2(2)
       write(*,9302) ksec2(3)
       write(*,9310) ksec2(4)
       write(*,9311) ksec2(5)
       ibit = 8
       iresol = ksec2(18) + ksec2(19)
       call prtbin(iresol,ibit,iout,ierr)
       write(*,9205) iout
       write(*,9312) ksec2(7)
       write(*,9313) ksec2(8)
       write(*,9314) ksec2(9)
       write(*,9315) ksec2(10)
       ibit = 8
       call prtbin(ksec2(11),ibit,iout,ierr)
       write(*,9211) iout
       write(*,9212) ksec2(12)
       write(*,9303) ksec2(13)
       write(*,9316) ksec2(14)
       go to 190
    endif
    !     *    representation type not catered for.

    write(*,9500) ksec2(1)
    go to 900

    !     *    print vertical coordinate parameters, if any.

190 continue

    if(ksec2(12)/=0)then
       write(*,9000)
       write(*,9400)
       write(*,9401)
       write(*,'(4x,f20.12,4x,f20.12)')(psec2(i),i=11,ksec2(12)+10)
    endif
    !     rotated and stretched grids introduced in edition 1.

    if(iedit<1) go to 900

    !     *    rotated grid information, if any.

    if(ksec2(1)==10.or.ksec2(1)==30.or.ksec2(1)==14 &
         .or.ksec2(1)==34.or.ksec2(1)==60.or.ksec2(1)==80 &
         .or.ksec2(1)==30)then
       write(*,9000)
       write(*,9220) ksec2(13)
       write(*,9221) ksec2(14)
       write(*,9222) psec2(1)
    endif
    !     *    stretched grid information, if any.

    if(ksec2(1)==20.or.ksec2(1)==30.or.ksec2(1)==24     &     
         .or.ksec2(1)==34.or.ksec2(1)==70.or.ksec2(1)==80)then
       write(*,9000)
       write(*,9230) ksec2(15)
       write(*,9231) ksec2(16)
       write(*,9232) psec2(2)
    endif
900 continue
9000 format(1h )
9001 format(1h ,'section 2 - grid description section.')
9002 format(1h ,'-------------------------------------')
9101 format(1h ,'data representation type(table 6).          ',i9)
9102 format(1h ,'j - pentagonal resolution parameter.         ',i9)
9103 format(1h ,'k - pentagonal resolution parameter.         ',i9)
9104 format(1h ,'m - pentagonal resolution parameter.         ',i9)
9105 format(1h ,'representation type(table 9)                ',i9)
9106 format(1h ,'representation mode(table 10).              ',i9)
9107 format(1h ,'not used.                                    ',i9)

9200 format(1h ,'(southern latitudes and western longitudes are negative.)')
9201 format(1h ,'number of points along a parallel.           ',i9)
9202 format(1h ,'number of points along a meridian.           ',i9)
9203 format(1h ,'latitude of first grid point.                ',i9)
9204 format(1h ,'longitude of first grid point.               ',i9)
9205 format(1h ,'resolution and components flag.               ',i8.8)
9206 format(1h ,'latitude of last grid point.                 ',i9)
9207 format(1h ,'longitude of last grid point.                ',i9)
9208 format(1h ,'i direction(east-west) increment.           ',i9)
9209 format(1h ,'j direction(north-south) increment.         ',i9)
9210 format(1h ,'number of parallels between pole and equator.',i9)
9211 format(1h ,'scanning mode flags(code table 8)            ',i8.8)
9212 format(1h ,'number of vertical coordinate parameters.    ',i9)

9220 format(1h ,'latitude of southern pole of rotation.       ',i9)
9221 format(1h ,'longitude of southern pole of rotation.      ',i9)
9222 format(1h ,'angle of rotation.                     ',f20.10)

9230 format(1h ,'latitude of pole of stretching.              ',i9)
9231 format(1h ,'longitude of pole of stretching.             ',i9)
9232 format(1h ,'stretching factor.                     ',f20.10)
9233 format(1h ,'number of points along a parallel varies.')
9234 format(1h ,'number of points.   parallel.(numbered from north to south)')
9235 format(1h , i5,16x,a10)
9236 format(1h ,'i direction(east-west) increment not given.')
9237 format(1h ,'j direction(north-south) increment not given.')

9301 format(1h ,'number of points along x axis.               ',i9)
9302 format(1h ,'number of points along y axis.               ',i9)
9303 format(1h ,'orientation of the grid.                     ',i9)
9304 format(1h ,'x direction increment.                       ',i9)
9305 format(1h ,'y direction increment.                       ',i9)
9306 format(1h ,'projection centre flag.                      ',i9)

9310 format(1h ,'latitude of sub-satellite point.             ',i9)
9311 format(1h ,'longitude of sub-satellite point.            ',i9)
9312 format(1h ,'diameter of the earth in x direction.        ',i9)
9313 format(1h ,'diameter of the earth in y direction.        ',i9)
9314 format(1h ,'x coordinate of sub-satellite point.         ',i9)
9315 format(1h ,'y coordinate of sub-satellite point.         ',i9)
9316 format(1h ,'altitude of the camera.                      ',i9)

9400 format(1h ,'vertical coordinate parameters.')
9401 format(1h ,'-------------------------------')

9500 format(1h ,'grprs2 :data representation type not catered for -',i4)

    return

  end subroutine grprs2
  subroutine grprs3(ksec3,psec3)
    implicit none
    !
    !**** grprs3 - print information from section 3 of grib code.
    !
    !     purpose.
    !     --------
    !
    !           print the information in the bit-map section
    !           section(section 3) of decoded grib data.
    !
    !**   interface.
    !     ----------
    !
    !           call grprs3(ksec0,ksec3,psec3)
    !
    !               input parameters.
    !               -----------------
    !
    !               ksec0 - array of decoded integers from section 0.
    !
    !               ksec3 - array of decoded integers from section 3.
    !
    !               psec3 - array of decoded reals from section 3.
    !
    !               output parameters.
    !               ------------------
    !
    !               none.
    !
    integer::  ksec3(*)

    real(kind=realgribkind) psec3(*)

    write(*,9000)
    write(*,9001)
    write(*,9002)

    if(ksec3(1)/=0)then
       write(*,9003) ksec3(1)
    else
       write(*,9004)
    endif
    write(*,9005) ksec3(2)

    write(*,9006) psec3(2)

9000 format(1h )
9001 format(1h ,'section 3 - bit-map section.')
9002 format(1h ,'-------------------------------------')
9003 format(1h ,'predetermined bit-map number.                ',i9)
9004 format(1h ,'no predetermined bit-map.')
9005 format(1h ,'missing data value for integer data.         ',i9)
9006 format(1h ,'missing data value for real data.        ',f20.6)

    return

  end subroutine grprs3
  subroutine grprs4(ksec4,psec4)
    implicit none
    !     
    !**** grprs4 - print information from section 4 of grib code.
    !     
    !     purpose.
    !     --------
    !     
    !     print the information in the binary data section
    !     section(section 4) of decoded grib data.
    !     
    !**   interface.
    !     ----------
    !     
    !     call grprs4(ksec4,psec4)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     ksec4 - array of decoded integers from section 4.
    !     
    !     psec4 - array of decoded reals from section 4.
    !     
    !     output parameters.
    !     ------------------
    !     
    !     none.
    !     
    !     method.
    !     -------
    !     
    !     fields printed as integers or reals.
    !     
    integer:: inum
    integer:: j210

    integer:: ksec4(*)
    real(kind=realgribkind) psec4(*)

    !     *    section 1 . print integer information from ksec4.


    write(*,9000)
    write(*,9001)
    write(*,9002)

    write(*,9003) ksec4(1)
    write(*,9004) ksec4(2)
    write(*,9005) ksec4(3)
    write(*,9006) ksec4(4)
    write(*,9007) ksec4(5)
    write(*,9008) ksec4(6)
    write(*,9009) ksec4(7)
    write(*,9010) ksec4(8)

    !     *    section 2. print real values from psec4.


    write(*,9000)

    inum = ksec4(1)
    if(inum<0)  inum = - inum
    if(inum>20) inum = 20
    !     print first inum values.

    write(*,9011) inum
    do 210 j210=1,inum
       write(*,9012) psec4(j210)
210 enddo
    !     *    section 9 . format statements. return to calling routine.
900 continue

9000 format(1h )
9001 format(1h ,'section 4 - binary data  section.')
9002 format(1h ,'-------------------------------------')
9003 format(1h ,'number of data values coded/decoded.         ',i9)
9004 format(1h ,'number of bits per data value.               ',i9)
9005 format(1h ,'type of data indicator.                      ',i9)
9006 format(1h ,'type of packing indicator.                   ',i9)
9007 format(1h ,'type of data representation.                 ',i9)
9008 format(1h ,'additional flags indicator.                  ',i9)
9009 format(1h ,'reserved.                                    ',i9)
9010 format(1h ,'number of values indicator.                  ',i9)
9011 format(1h ,'first ',i4,' data values.')
9012 format(1h ,f30.15)

    return

  end subroutine grprs4
  subroutine gsbite(ks,kd,kskst,ksize,kskbtw,k,kbpw,kmask,yadir)
    !
    !     gsbite: vectorising extraction/insertion of bits from/to bitstream
    !
    ! input:
    !     ks:     if yadir='d', input bit stream, else output bit stream
    !     kd:     if yadir='d', output words, else input words
    !     kskst:  number of bits skipped at beginning of ks
    !     ksize:  number of bits to be extracted to one word of kd
    !     kskbtw: number of bits skipped between two words to be extracted
    !     k:      number of words to be extracted into kd(if <=0, only
    !             calculate kbpw and kmask
    !     kbpw:   number of bits per word in ks, calculated if 0
    !     kmask:  masks for bit patterns, calculated if kmask(2)==0
    !     yadir:  direction of conversion: 'd' for decoding, i.e.
    !             extract words kd(1...k)  from bits ks(kskst+1....)
    !             if not 'd', encode, i.e. pack words kd(1....k) into bits
    !             ks(kskst+1.....kskst+k*(ksize+kskbtw))
    !
    ! output:
    !     ks,kd:  see above
    !     kskst:  updated to nr of bits used, i.e. to kskst+k*(ksize+kskbtw)
    !     kbpw:  (if 0 on input): number of bits in each word of ks
    !     kmask: (if(kmask(2) was 0 on input): bit pattern masks
    !
    !                                                     g.j.cats 08 dec 87
    !
    implicit none
    integer:: ks(*),kd(*),kmask(*)
    character(len=1):: yadir
    logical     lodecd
    integer:: kbpw,is,j,k,ksize,istep
    integer:: kskbtw,iskws,kskst,ill,ibdl,istd,ists,ilcf,ishf,jbd
    integer:: ienbs,iskw,ista,iskb,imask,ish,ibs,iend,ios,iod,ji
    !     use is made of the us military bit manipulation functions
    !     ibset, not, ishft, iand and ior.

    !     1.  complete kbpw and kmask, return if 0 words are to be extracted

    if(kbpw==0)then
       is=ks(1)
       ks(1)=1
1101   continue
       if(ks(1)/=0)then
          kbpw=kbpw+1
          ks(1)=ishft(ks(1),1)
          goto 1101
       endif
       ks(1)=is
    endif
    if(kmask(2)==0)then
       kmask(kbpw+1)=0
       do 1110 j=kbpw,1,-1
          kmask(j)=ibset(kmask(j+1),kbpw-j)
1110   enddo
    endif
    if(k<=0)return

    !     2.  preset kd to 0 if kd is output i.e. when decoding

    lodecd = yadir == 'D' .or.yadir == 'd'
    if(lodecd) then
       do 2101 j=1,k
          kd(j)=0
2101   enddo
    endif
    if(ksize<=0)return

    !     3.  calculate several parameters for looping(for efficiency, the
    !         code of sections 3.3 and 3.4 for k=1 is separated into 3.2)

    !     3.1 number of bits used per word, initial nr of skipped bits

    istep=ksize+kskbtw
    iskws=kskst
    !
    !     3.2 vector loop length and step size in kd if k=1;ks step irrelvnt
    !
    if(k==1)then
       ill=1
       ibdl=2
       istd=1
       ists=1
    else

       !     3.3 step sizes in ks,kd: inverse of largest factor of istep,kbpw
       !
       ilcf=kbpw
       ishf=istep
331    continue
       if(ilcf==ishf)goto 332
       if(ilcf==1)goto 332
       if(ilcf>ishf)then
          ilcf=ilcf-ishf
       else
          ishf=ishf-ilcf
       endif
       goto 331
332    continue
       istd=kbpw/ilcf
       ists=istep/ilcf
       !     3.4 vector loop length and switch-over point for smaller loop
       !
       ill=(k-1)/istd+1
       ibdl=k-(ill-1)*istd
    endif

    !     4.  loop over first istd words of kd(trails the vector loop)
    !
    do 790 jbd=1,istd
       !
       !     4.1 last bit in ks to be treated
       !
       ienbs=iskws+ksize
       !
       !     4.2 nr of words of ks to be skipped, nr of bits in those and this
       !
       iskw=iskws/kbpw
       ista=iskw*kbpw
       iskb=iskws-ista

       !     4.3 mask and left shift for the remaining bits
       !
       imask=kmask(iskb+1)
       ish=ksize+iskb
       !
       !     4.4 position of current word of ks
       !
       ibs=iskw+1

       !     5.  loop over words of ks contributing to one word of kd
       !
500    continue
       !
       !     5.1 update shift and last bit in current word
       !
       ish=ish-kbpw
       iend=ista+kbpw

       !     5.2 is last bit of current word outside range to be extracted
       !
       if(iend>ienbs)then
          ish=ienbs-iend
          imask=iand(imask,not(kmask(kbpw+ish+1)))
       endif

       !     5.3 initial offsets for vector elements in vector loop
       !
       ios=0
       iod=0
       !
       !     6.  vector loop is over repeatedly occurring bitpatterns/masks
       !
       if(lodecd) then
          do 611 ji=1,ill
             kd(jbd+iod)=ior(kd(jbd+iod),ishft(iand(imask,ks(ibs+ios)),ish))
             iod=iod+istd
             ios=ios+ists
611       enddo
       else
          do 612 ji=1,ill
             ks(ibs+ios)=ior(  &
                  iand(      ks(ibs+ios),      not(imask)), &
                  iand(ishft(kd(jbd+iod),-ish),    imask ))
             iod=iod+istd
             ios=ios+ists
612       enddo
       endif

       !     7.  end loops
       !
       !     7.1 prepare for end of loop over words of ks witihn one kd word
       !
       ista=ista+kbpw
       !
       !     7.2 next word of kd if extraction not completed
       !
       if(ista<ienbs)then
          imask=kmask(1)
          ibs=ibs+1
          goto 500
       endif

       !     7.8 prepare for end of loop over first words of kd
       !
       if(jbd==ibdl)ill=ill-1
       iskws=iskws+istep
       !
       !     7.9 end loop over first words of kd
       !
790 enddo

    !     8.  finished: update kskst and return
    !
    kskst=kskst+k*istep
    !
    !     8.5 swap bytes on vax
    return
  end subroutine gsbite
  !-----------------------------------------------------------------------
  subroutine inscal(pdata,kdata,klen,pref,pscale)
    implicit none
    !
    !**** INSCAL - Vectorise calculation of increments.
    !
    !     Purpose.
    !     --------
    !
    !           Vectorise calculation of increments.
    !
    !**   Interface.
    !     ----------
    !
    !           CALL INSCAL(PDATA,KDATA,KLEN,PREF,PSCALE)
    !
    !               Input Parameters.
    !               -----------------
    !
    !               PDATA      - Array of floating point values.
    !               KLEN       - Number of values to be converted.
    !               PREF       - Reference value.
    !               PSCALE     - Scale factor.
    !
    !               Output Parameters.
    !               -----------------
    !
    !               KDATA      - Array of integer increments
    !
    !     Method.
    !     -------
    !
    !           The reference value is subtracted from each value,
    !           and the result is then divided by the scale factor.
    !
    !     Externals.
    !     ----------
    !
    !           None.
    !
    !     Reference.
    !     ----------
    !
    !           WMO Manual on Codes re GRIB representation.
    !
    !     Comments.
    !     --------
    !
    !           PDATA and KDATA are really the same array. This routine
    !           is just a device to force vectorisation on the Cray,
    !           without the necessity of using another array.
    !
    !           Routine contains section 0 , 1 and 9.
    !
    !     Author.
    !     -------
    !
    !           J. Hennessy     ECMWF     25.06.91
    !
    !     Modifications.
    !     _____________
    !
    !           None.
    !
    !     -----------------------------------------------------------------
    !
    integer:: klen
    integer:: j101

    real(kind=realgribkind) pref
    real(kind=realgribkind) pscale

    real(kind=realgribkind) pdata(klen)
    integer:: kdata(klen)
    !*     calculation of increments.
    do 101 j101 = 1,klen
       kdata(j101) = nint((pdata(j101) - pref)/pscale )
101 enddo
    return
  end subroutine inscal
  subroutine inxbit(kgrib,kleng,knspt,kparm,knum,kbit, &
       kblen,hfunc,kret)
    !     
    !**** inxbit - insert/extract bits consecutively in/from a given array
    !     
    !     purpose.
    !     --------
    !     
    !     take rightmost kblen bits from knum words of kparm
    !     and insert them consecutively in kgrib, starting at
    !     bit after knspt or vice versa.
    !     
    !**   interface.
    !     ----------
    !     
    !     call inxbit(kgrib,kleng,knspt,kparm,knum,kbit,
    !     c                   kblen,kret)
    !     
    !     input parameters.
    !     -----------------
    !     
    !     kgrib      - array containing bitstream.
    !     kleng      - length(words) of this array.
    !     knspt      - bit number after which insertion or
    !     extraction starts.
    !     kparm      - array from which bits are taken for
    !     insertion in the bitstream or to which
    !     bits are extracted from the bitstream.
    !     kbit       - number of bits in computer word.
    !     knum       - number of bit fields inserted/extracted.
    !     kblen      - number of bits per bit field.
    !     hfunc      - requested function.
    !     'c' to insert bits in bitstream,
    !     'd' to extract bits from bitstream.
    !     
    !     output parameters.
    !     ------------------
    !     
    !     knspt      - bit number of last bit inserted/extracted.
    !     
    !     kret       - return code.

    !     0 , no error encountered.
    !     1 , insertion/extraction exceeded
    !     array boundary.
    !     
    !     method.
    !     -------
    !     
    !     word and offset pointer calculated before calling
    !     insertion/extraction routines.
    !     
    !     externals.
    !     ----------
    !     gsbite
    !     
    !     reference.
    !     ----------
    !     
    !     comments.
    !     ---------
    !     
    !     cray version of routine.
    !     this routine should only be used on the cray as it
    !     contains a call to gsbite, a vectorising version of
    !     gbyte(s) and sbyte(s).
    !     for hirlam purposes, we assume gbytes and sbytes not to
    !     be available.  because gsbite handles everything, albeit
    !     less efficiently, calls to gbytes and sbytes have been
    !     eliminated.
    !     routine contains sections 0 to 3 and section 9.
    !     
    !     author.
    !     -------
    !     
    !     j. hennessy      ecmwf      18.06.91
    !     
    !     modifications.
    !     --------------
    !     
    !     j. hennessy      ecmwf      08.11.91
    !     parameter kmach removed from list of input parameters.
    !     
    !     ----------------------------------------------------------------
    !     
    implicit none
    !     
    integer:: imask(65)
    integer:: ind
    integer:: inum

    integer:: ipr

    !     
    integer:: kbit
    integer:: kblen
    integer:: kgrib(kleng)
    integer:: kleng
    integer:: knspt
    integer:: knum
    integer:: kparm(*)
    integer:: kret

    integer:: j901

    character(len=1):: hfunc

    !     values in imask are set in the first call to routine gsbite, and
    !     are used in subsequent calls.
    !     
    save imask

    !     force routine gsbite to calculate bit-masks first time through.

    data imask(2) /0/
    !     
    !     debug print switch.
    !     
    data ipr /0/
    !     
    !     ----------------------------------------------------------------
    !     
    !     *    section 1 . set initial values.
    !     ----------------------------------------------------------------
    !     

    !     
    if(ipr==1)then
       write(*,*) 'inxbit : section 1.'
       write(*,*) '         input values used -'
       write(*,9002) knspt
       write(*,9004) kbit
       write(*,9005) hfunc
    endif
    !     
    kret = 0
    !     
    !     ----------------------------------------------------------------

    !     *    section 2 . bit insertion/extraction.
    !     ----------------------------------------------------------------
    !     

    !     
    if(ipr==1) write(*,*) 'inxbit : section 2.'
    call gsbite(kgrib,kparm,knspt,kblen,0,knum, kbit,imask,hfunc)


    !     *    section 3 . check out of range.


    if(ipr==1) write(*,*) 'inxbit : section 3.'

    ind = knspt / kbit
    if(ind>kleng) then
       kret = 1
       write(*,9001) ind , kleng
    endif

    !     *    section 9 . return to calling routine. format statements.


    if(ipr==1)then
       inum = knum
       if(inum>360)then
          inum = 360
          write(*,9007) inum
       endif
       do 901 j901=1,inum
          if(hfunc=='C'.or.hfunc=='c') then
             write(*,9006) kparm(j901)
          else
             write(*,9008) kparm(j901)
          endif
901    enddo
       write(*,*) 'inxbit : section 9.'
       write(*,*) '         output values set -'
       write(*,9002) knspt
    endif

9001 format(1h ,'inxbit : word ',i8,' is outside array bounds ',i8)

9002 format(1h ,'         knspt = ',i8)

9003 format(1h ,'inxbit : word is',i8,', bit offset is ',i2)

9004 format(1h ,'         kbit  = ',i8)

9005 format(1h ,'         hfunc = ',a)

9006 format(1h ,'         inserted value = ',i20)

9007 format(1h ,'         first ',i9,' values.')

9008 format(1h ,'         extracted value = ',i20)

    return

  end subroutine inxbit
  subroutine inxmap(kgrib,kleng,kmapt,psec4,ksize,kbits,ksbmap, &
       pmiss,hfunc,kpr,kret)
    implicit none
    !     
    !**** INXMAP - Bit map handling for routine GRIBEX.
    !     
    !     Purpose.
    !     --------

    !     1) Extract a bit-map from an array of GRIB coded data and
    !     insert missing data value in appropriate places in the
    !     array of already decoded data values.
    !     
    !     2) Generate a bit-map and insert in array of GRIB coded
    !     data and remove missing data values from the array of
    !     values being encoded in GRIB code.
    !     
    !**   Interface.
    !     ----------
    !     
    !     CALL INXMAP(KGRIB,KLENG,KMAPT,PSEC4,KSIZE,KBITS,KSBMAP,
    !     C                     PMISS,HFUNC,KPR,KRET)
    !     
    !     Input Parameters.
    !     -----------------
    !     
    !     KGRIB      - Array of GRIB data being coded/decoded.
    !     
    !     KLENG      - Length of this array.
    !     
    !     KMAPT      - Bit-pointer to start of bit-map in
    !     array KGRIB.
    !     
    !     PSEC4      - Array of data values decoded or to be
    !     coded in GRIB.
    !     
    !     KSIZE      - Size of bit-map ie number of values,
    !     including missing data values, in array
    !     PSEC4.
    !     
    !     KBITS      - Number of bits in computer word.
    !     
    !     PMISS      - Value indicating missing data in array
    !     PSEC4.
    !     
    !     HFUNC      - 'C' , GRIB data being coded.
    !     'M' , GRIB data being coded in fixed length
    !     messages.
    !     'D' , GRIB data being decoded.
    !     
    !     KPR        - Debug print switch.
    !     0 , No printout.
    !     1 , Debug printout.
    !     
    !     Output Parameters.
    !     ------------------
    !     
    !     KSBMAP     - Bit-map flag.
    !     Used only when decoding data.
    !     -2 , All bits in the bit-map set to 1.
    !     There is no missing data.
    !     -4 , Some points have no data. User
    !     supplied value for missing data
    !     indicator in appropriate places in
    !     the array PSEC4.
    !     
    !     KSIZE      - When coding data, the number of real data
    !     values in array PSEC4. Not changed when
    !     decoding.
    !     
    !     KRET       - Return code.
    !     0 , No error encountered.
    !     1 , Error in routine INXBIT.
    !     2 , Bit-map size exceeds maximum.
    !     3 , Invalid function requested.
    !     
    !     Method.
    !     -------
    !     
    !     The bit-map contains 1 where valid data exists and 0
    !     where data is missing. The corresponding data array
    !     contains valid data and the missing data indicator value.
    !     
    !     Externals.
    !     ----------
    !     
    !     INXBIT
    !     
    !     Reference.
    !     ----------
    !     
    !     See routine GRIBEX.
    CHARACTER(len=*)::HFUNC
    CHARACTER(len=1)::YFUNC
    INTEGER,parameter::JPMAP=116000 !41472

    INTEGER:: IMAP(JPMAP)
    INTEGER:: INEXT
    INTEGER:: J320
    INTEGER:: J330
    INTEGER:: J420
    INTEGER:: J430
    INTEGER:: J431
    INTEGER:: KBITS
    INTEGER:: KLENG
    INTEGER:: KGRIB(kleng)
    INTEGER:: KMAPT
    INTEGER:: KPR
    INTEGER:: KRET
    INTEGER:: KSBMAP
    INTEGER:: KSIZE

    REAL(KIND=REALGRIBKIND) PMISS
    REAL(KIND=REALGRIBKIND) PSEC4(*)

    REAL(KIND=REALGRIBKIND) ZSEC4(JPMAP)





    IF(KPR==1)THEN
       WRITE(*,*) 'INXMAP : Section 1.'
       WRITE(*,*) '         Input values used -'
       WRITE(*,9004) HFUNC
       WRITE(*,9005) KSIZE
    ENDIF
    !     Reset return code to

    KRET  = 0
    INEXT = 0
    YFUNC = HFUNC

    !     *    Section 2 .  Check input parameters.


    IF(KPR==1) WRITE(*,*) 'INXMAP : Section 2.'

    !     *    Check that bit-map size does not exceed maximum permitted.

    IF(JPMAP<KSIZE)THEN
       WRITE(*,9001) KSIZE , JPMAP
       KRET = 2
       GO TO 900
    ENDIF

    !     *    Section 3 . Decoding of bit-map and data.



    IF(KPR==1) WRITE(*,*) 'INXMAP : Section 3.'

    IF(YFUNC=='D'.or.YFUNC=='d')THEN

       !     *           Extract bit-map from GRIB coded data.

       call inxbit(kgrib,kleng,kmapt,imap,ksize,kbits,1,yfunc,kret)

       IF(KRET/=0)THEN
          WRITE(*,9003)
          KRET = 1
          GO TO 900
       ENDIF
       !     *           Copy data to temporary array and insert missing data
       !     indicator in temporary array, in accordance with the
       !     bit map values.

       DO 320 J320=1,KSIZE
          IF(IMAP(J320)==0)THEN
             ZSEC4(J320) = PMISS
          ELSE
             INEXT       = INEXT + 1
             ZSEC4(J320) = PSEC4(INEXT)
          ENDIF
320    enddo
       !     *           Set bit-map flag, and if missing data is indicated
       !     transfer data to original array.

       IF(INEXT/=KSIZE)THEN
          KSBMAP = -4
          DO 330 J330=1,KSIZE
             PSEC4(J330) = ZSEC4(J330)
330       enddo
       ELSE
          KSBMAP = -2
       ENDIF

       GO TO 900

    ENDIF
    !     *    Section 4 . Generation of bit-map.


    IF(KPR==1) WRITE(*,*) 'INXMAP : Section 4.'

    IF(YFUNC=='C'.OR.YFUNC=='M'.or. &
         yfunc=='c'.or.yfunc=='m' )THEN

       !     *           Copy data to temporary array and remove missing data
       !     indicator in temporary array, generating the bit-map
       !     in accordance with the missing data values.

       DO 420 J420=1,KSIZE
          IF(abs(PSEC4(J420)-PMISS)<1.e-14_realgribkind)THEN
             IMAP(J420)   = 0
          ELSE
             IMAP(J420)   = 1
             INEXT        = INEXT + 1
             ZSEC4(INEXT) = PSEC4(J420)
          ENDIF
420    enddo
       !     *           Insert bit-map in GRIB coded data.

       call inxbit(kgrib,kleng,kmapt,imap,ksize,kbits,1,yfunc,kret)

       IF(KRET/=0)THEN
          WRITE(*,9003)
          KRET = 1
          GO TO 900
       ENDIF
       !     *           If missing data is indicated transfer data to original
       !     array.

       IF(INEXT/=KSIZE)THEN
          DO 430 J430=1,INEXT
             PSEC4(J430) = ZSEC4(J430)
430       enddo

          IF(YFUNC=='M'.or.yfunc=='m')THEN
             !     Fixed length messages required, even though
             !     a bit map is used. The otherwise unused part
             !     of the array is set to a genuine data value
             !     so that extraction of minimum and maximum
             !     values remain correct. Number of data values
             !     includes these padding values.
             !     
             DO 431 J431=INEXT+1,KSIZE
                PSEC4(J431) = PSEC4(1)
431          enddo
          ELSE

             !     Return number of real values(excluding
             !     missing data values).

             KSIZE = INEXT
          ENDIF
       ENDIF

       GO TO 900

    ENDIF
    !     *    Invalid function requested.

    WRITE(*,9002) YFUNC
    KRET = 3

    !     *    Section 9 . Return to calling routine. Format statements.
900 CONTINUE

    IF(KPR==1)THEN
       WRITE(*,*) 'INXMAP : Section 9.'
       WRITE(*,*) '         Output values set -'
       WRITE(*,9005) KSIZE
       WRITE(*,9006) KSBMAP
    ENDIF

9001 FORMAT(1H ,'INXMAP : Bit-map size is ',I6,', maximum allowed', &
         ' is ',I6,'.')

9002 FORMAT(1H ,'INXMAP : Invalid function requested - ',A1)

9003 FORMAT(1H ,'INXMAP : Error reported by routine INXBIT.')

9004 FORMAT(1H ,'         HFUNC  = ',A1)

9005 FORMAT(1H ,'         KSIZE  = ',I6)

9006 FORMAT(1H ,'         KSBMAP = ',I6)

    RETURN

  END subroutine inxmap
  subroutine maxmin(parray,klen,pmax,pmin,pmdi)
    implicit none

    !**** maxmin - get maximum and minimum values.
    !     get maximum and minimum values from an array of
    !     floating point numbers.
    !     call maxmin(parray,klen,pmax,pmin,pmdi)
    !     input parameters.
    !     -----------------
    !     
    !     parray     - array of numbers.
    !     klen       - last word of this array.
    !     pmdi       - missing data indicator.
    !     output parameters.
    !     pmax       - maximum value.
    !     pmin       - minimum value.
    !     intrinsic functions max and min are used.
    !     j. hennessy      ecmwf      18:06:91
    integer:: j110

    integer,intent(in)::klen

    real(kind=realgribkind) pmax
    real(kind=realgribkind) pmin
    real(kind=realgribkind),intent(in)::pmdi

    real(kind=realgribkind),intent(in)::parray(klen)


    pmax = parray(1)
    pmin = parray(1)

    do 101 j110 = 2 , klen
       if(abs(parray(j110)-pmdi)>1.e-14_realgribkind)then
          pmax=parray(j110)
          pmin=parray(j110)
       endif
101 enddo
    if(abs(pmax-pmdi)<1.e-14_realgribkind)then
       pmax=0._realgribkind
       pmin=0._realgribkind
       return
    endif

    !     extract maximum and minimum values.
    do 110 j110 = 1 , klen
       if(abs(parray(j110)-pmdi)>1.e-14_realgribkind)then
          pmax = max(pmax,parray(j110))
          pmin = min(pmin,parray(j110))
       endif
110 enddo
    return
  end subroutine maxmin

  integer function ngbgrb(kw,kb)
    implicit none
    !     ngbgrb: read bit kb in word kw(kb=1 for leftmost bit)
    integer:: kw,kb,ib
    ib=256/2**kb
    ngbgrb=mod(kw/ib,2)
    return
  end function ngbgrb

  subroutine prtbin(kin,knbit,kout,kerr)
    implicit none
    !**** prtbin - binary to decimal conversion.
    !     produces a decimal number with ones and zeroes
    !     corresponding to the ones and zeroes of the input
    !     binary number.
    !     eg input number 1011 binary, output number 1011 decimal.
    !     kin   - integer variable containing binary number.
    !     
    !     knbit - number of bits in binary number.
    !     
    !     output parameters.
    !     -----------------
    !     
    !     kout  - integer variable containing decimal value
    !     with ones and zeroes corresponding to those of
    !     the input binary number.
    !     
    !     kerr  - 0, if no error.
    !     1, number of bits in binary number exceeds
    !     maximum allowed or is less than 1.
    !     odd numbers have a binary representation ending in 1, even
    !     numbers end in 0.

    integer:: idec
    integer:: ik
    integer:: itemp
    integer:: j102
    integer:: kerr
    integer:: kin
    integer:: knbit
    integer:: kout

    !     check length of binary number to ensure decimal number
    !     generated will fit in the computer word - in this case will
    !     it fit in a cray 48 bit integer?

    if(knbit<1.or.knbit>14)then
       kerr = 1
       write(*,9000) knbit
       go to 900
    else
       kerr = 0
    endif

    !     *    section 1. generate required number.


    kout = 0
    ik   = kin
    idec = 1

    do 102 j102=1,knbit
       itemp = ik -((ik/2)*2 )
       kout  = kout + itemp * idec
       ik    = ik / 2
       idec  = idec * 10
102 enddo

    !     *    section 9. format statements. return to calling routine.

900 continue

9000 format(1h ,'prtbin : error in binary number length - ',i3, &
         ' bits.')

    return
  end subroutine prtbin
  subroutine qu2reg(pfield,kpoint,klat,klon,kcode)
    implicit none
    !     convert quasi-regular grid data to regular,
    !     using either a linear or cubic interpolation.
    !     pfield     - array containing quasi-regular grid
    !     data.
    !     
    !     kpoint     - array containing list of the number of
    !     points on each latitude(or longitude) of
    !     the quasi-regular grid.
    !     
    !     klat       - number of latitude lines
    !     
    !     klon       - number of longitude lines
    !     
    !     kcode      - interpolation required.
    !     1 , linear - data quasi-regular on
    !     latitude lines.
    !     3 , cubic - data quasi-regular on
    !     latitude lines.
    !     11, linear - data quasi-regular on
    !     longitude lines.
    !     13, cubic - data quasi-regular on
    !     longitude lines.
    !     
    !     output parameters.
    !     pfield     - array containing regular grid data.
    !     data is interpolated and expanded into a temporary array,
    !     which is then copied back into the user's array.
    !     routine aborts if an invalid interpolation is requested or
    !     field size exceeds array dimensions.
    integer:: icode
    integer:: ilii
    integer:: ilio
    integer:: iquano
    integer:: iregno


    integer:: j210
    integer:: j220
    integer:: j225
    integer:: j230
    integer:: j240

    integer:: kcode
    integer:: klat
    integer:: klon
    integer:: kpoint(*)

    !     maximum number of latitudes(or longitudes), for which arrays
    !     are dimensioned.

    integer,parameter::jpmax=320

    real(kind=realgribkind) pfield(*)
    real(kind=realgribkind) ztemp(jpmax*jpmax*2)
    real(kind=realgribkind) zline(jpmax*2)
    real(kind=realgribkind) zwork(0:jpmax*2+2,3)


    !     *    section 1. set initial values.

    !     check input parameters.

    if(kcode/=1.and.kcode/=3.and.kcode/=11.and.kcode/=13)then
       write(*,9001) kcode
       call stop_program( 'qu2reg stop 1')
    endif

    if(klat>jpmax)then
       write(*,9002) klat , jpmax
       call stop_program( 'qu2reg stop 2')
    endif

    if(klon>jpmax*2)then
       write(*,9003) klat , jpmax*2
       call stop_program( 'qu2reg stop 3')
    endif
    !     set array indices to 0.

    ilii  = 0
    ilio  = 0

    !     establish values of loop parameters.

    if(kcode>10)then
       !     quasi-regular along longitude lines.
       iquano = klon
       iregno = klat
       icode  = kcode - 10
    else
       !     quasi-regular along latitude lines.
       iquano = klat
       iregno = klon
       icode  = kcode
    endif

    !     *    section 2. interpolate field from quasi to regular grid.


    do 230 j230=1,iquano
       if(iregno/=kpoint(j230))then
          !     line contains less values than required,so
          !     extract quasi-regular grid values for a line

          do 210 j210=1,kpoint(j230)
             ilii        = ilii+1
             zline(j210) = pfield(ilii)
210       enddo

          !     and interpolate this line.

          call rowina(zline,iregno,kpoint(j230),zwork,icode)

          !     add regular grid values for this line to the temporary
          !     array.

          do 220 j220=1,iregno
             ilio        = ilio+1
             ztemp(ilio) = zline(j220)
220       enddo

       else
          !     line contains the required number of values, so add
          !     this line to the temporary array.

          do 225 j225=1,iregno
             ilio        = ilio+1
             ilii        = ilii+1
             ztemp(ilio) = pfield(ilii)
225       enddo

       endif

230 enddo
    !     copy temporary array to user array.

    do 240 j240=1,klon*klat
       pfield(j240) = ztemp(j240)
240 enddo

    !     *        section 9. return to calling routine. format statements.

9001 format(1h ,'qu2reg : invalid interpolation type code = ',i3)

9002 format(1h ,'qu2reg : number of latitudes is ',i4,', maximum ',     &
         'allowed is ',i3,'.')

9003 format(1h ,'qu2reg : number of longitudes is ',i4,', maximum ',     &
         'allowed is ',i3,'.')

    return

  end subroutine qu2reg
  subroutine rorint(pdata,kdata,klen,hdir)
    implicit none
    !           Converts real arrays to integer and vice versa.
    !               KDATA      - Array of integer increments
    !                            Input for 'R' function.
    !               PDATA      - Array of floating point values.
    !                            Input for 'I' function.
    !               KLEN       - Number of values to be converted.
    !               HDIR       - 'R', convert integer to real.
    !                            'I', convert real to integer.
    !               KDATA      - Array of integer increments
    !                            Output for 'I' function.
    !               PDATA      - Array of floating point values.
    !                            Output for 'R' function.

    character(len=1):: hdir
    integer::klen
    integer::j

    real(kind=realgribkind) pdata(klen)
    integer:: kdata(klen)

    if(hdir=='i'.or.hdir=='I')then
       do j = 1,klen
          kdata(j) = nint(pdata(j))
       enddo
    else
       do j = 1,klen
          pdata(j) = real(kdata(j),realgribkind)
       enddo
    endif

    return
  end subroutine rorint
  subroutine rowina(p,ko,ki,pw,kcode)
    implicit none
    !
    !**** rowina - interpolation of row of values.
    !
    !     purpose.
    !     --------
    !
    !           interpolate a row of values.
    !
    !**   interface.
    !     ----------
    !
    !           call rowina(p,ko,ki,pw,kcode)
    !
    !               input parameters.
    !               -----------------
    !
    !               p     - row of values to be interpolated.
    !                       dimension must be at least ko.
    !
    !               ko    - number of values required.
    !
    !               ki    - number of values in p on input.
    !
    !               pw    - working array.
    !                       dimension must be at least(0:ko+2,3).
    !
    !               kcode - interpolation required.

    !                       1 , linear.
    !                       3 , cubic.
    !
    !               output parameters.
    !               ------------------
    !
    !               p     - now contains ko values.
    !
    !     method.
    !     -------
    !
    !           linear or cubic interpolation performed as required.
    !
    !     externals.
    !     ----------
    !
    !           scm0
    !
    !     reference.
    !     ----------
    !
    !           none.
    !
    !     comments.
    !     ---------
    !
    !           this is a version of rowint which conforms to ansi
    !           standards, achieved by passing the work array as a
    !           parameter and changing lower case letters to upper case.
    !
    !     author.
    !     -------
    !
    !           j. hennessy     ecmwf     09.10.91
    !
    !     modifications.
    !     --------------
    !
    !           j. hennessy     ecmwf     07.01.92

    integer:: ko,kcode,jl,ki,ip
    real(kind=realgribkind) zrdi,zdo,zpos,zwt,zwt1
    real(kind=realgribkind) p(ko),pw(0:ko+2,3)
    !
    if(kcode==1) then
       do 102 jl=1,ki
          pw(jl,1)=p(jl)
102    enddo
       pw(ki+1,1)=p(1)
       zrdi=real(ki,realgribkind)
       zdo=1._realgribkind/real(ko,realgribkind)
       !
       do 105 jl=1,ko
          zpos=real(jl-1,realgribkind)*zdo
          zwt=zpos*zrdi
          ip=int(zwt)
          zwt=zwt-real(ip,realgribkind)
          p(jl)=(1._realgribkind-zwt)*pw(ip+1,1)+zwt*pw(ip+2,1)
105    enddo
       !
    elseif(kcode==3) then
       do 302 jl=1,ki
          pw(jl,1)=p(jl)
302    enddo
       pw(0,1)=p(ki)
       pw(ki+1,1)=p(1)
       pw(ki+2,1)=p(2)
       do 305 jl=1,ki
          pw(jl,2)= - pw(jl-1,1)/3._realgribkind - 0.5_realgribkind*pw(jl,1)&
               + pw(jl+1,1)    - pw(jl+2,1)/6._realgribkind
          pw(jl+1,3)=   pw(jl-1,1)/6._realgribkind - pw(jl,1)     &
               + 0.5_realgribkind*pw(jl+1,1) + pw(jl+2,1)/3._realgribkind
305    enddo
       call scm0(pw(1,2),pw(2,3),pw(1,1),pw(2,1),ki)
       zrdi=real(ki,realgribkind)
       zdo=1._realgribkind/real(ko,realgribkind)
       do 310 jl=1,ko
          zpos=real(jl-1,realgribkind)*zdo
          zwt=zpos*zrdi
          ip=int(zwt+1._realgribkind)
          zwt=zwt+1._realgribkind-real(ip,realgribkind)
          zwt1 = 1._realgribkind - zwt
          p(jl)=((3._realgribkind-2._realgribkind*zwt1)*pw(ip,1) + &
               zwt*pw(ip,2))*zwt1*zwt1 +     &
              ((3._realgribkind-2._realgribkind*zwt) *pw(ip+1,1) - &
               zwt1*pw(ip+1,3))*zwt*zwt
310    enddo
    else
       write(*,9001) kcode
       call stop_program('') 
    endif

    return

9001 format(1h ,'rowina : invalid interpolation code = ',i4)

  end subroutine rowina
  subroutine scm0(pdl,pdr,pfl,pfr,klg)
    !     
    !**** scm0   - apply scm0 limiter to derivative estimates.
    !     
    !     m. hortal    ecmwf february 1991  closely following d. williamson
    !     
    !     apply scm0 limiter to derivative estimates.
    !     
    !     output:
    !     pdl   = the limited derivative at the left edge of the interval
    !     pdr   = the limited derivative at the right edge of the interval
    !     
    !     inputs
    !     pdl   = the original derivative at the left edge
    !     pdr   = the original derivative at the right edge
    !     pfl   = function value at the left edge of the interval
    !     pfr   = function value at the right edge of the interval
    !     klg  = number of intervals where the derivatives are limited
    !     
    implicit none
    integer:: klg
    real(kind=realgribkind):: pdl(klg),pdr(klg),pfl(klg),pfr(klg)
    real(kind=realgribkind):: zeps,zfac,zalpha,zbeta
    integer:: jl

    !     define constants

    zeps=1.e-12_realgribkind
    zfac=3._realgribkind*(1._realgribkind-zeps)

    do 200 jl=1,klg
       if(abs(pfr(jl)-pfl(jl))>zeps) then
          zalpha=pdl(jl)/(pfr(jl)-pfl(jl))
          zbeta =pdr(jl)/(pfr(jl)-pfl(jl))
          if(zalpha<=0._realgribkind) pdl(jl)=0._realgribkind
          if(zbeta <=0._realgribkind) pdr(jl)=0._realgribkind
          if(zalpha>zfac) pdl(jl)=zfac*(pfr(jl)-pfl(jl))
          if(zbeta >zfac) pdr(jl)=zfac*(pfr(jl)-pfl(jl))
       else
          pdl(jl)=0._realgribkind
          pdr(jl)=0._realgribkind
       endif
200 enddo
  end subroutine scm0
  subroutine setpar(kbit,kneg,kpr)
    !
    !**** setpar - set number of bits in word. set maximum negative integer.
    !
    !     purpose.
    !     --------
    !
    !           set number of bits in word. set maximum negative integer.
    !
    !**   interface.
    !     ----------
    !
    !           call setpar(kbit,kneg,kpr)

    !               input parameters.
    !               -----------------
    !
    !               kpr        - debug print switch.

    !                            1 , print out.
    !                            0 , no print out.
    !
    !               output parameters.
    !               ------------------
    !
    !               kbit       - number of bits in computer word.
    !
    !               kneg       - maximum negative integer.
    !
    !     author.
    !     -------
    !
    !           g. cats, knmi/hirlam, 12 feb 1997
    !
    !     ------------------------------------------------------------------
    !
    implicit none
    integer:: kbit
    integer:: kneg
    integer:: kpr


    kbit = 32

    kneg = -2147483647

    if(kpr==1)then
       write(*,*) ' setpar output values: kbit=',kbit,', kneg=',kneg
    endif
  end subroutine setpar
end module mod_grib
