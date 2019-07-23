module asimof
  use decomp,only:stop_program
  use mod_grib
  implicit none
  private

  integer,parameter::jpnfx=8         !     max number of files opened 
  integer,save::nlrec(jpnfx)         !     direct-access record length in computer words
  integer,save::nfld(jpnfx)          !     number of user fields in the file
  integer,save::munits(jpnfx)        !     open units
  integer,save::nbpwf(jpnfx)         !     number of bits in units of pointers and lengths
  integer,save::mopen(jpnfx)         !     mode of file opening


  integer,parameter::jpnkey=30       !     jpnkey = length of the key

  integer,parameter::jpnfdx=50000    !     jpnfdx = max number of records 25000
  integer,save::blocks(jpnkey,jpnfdx,jpnfx)

  integer,parameter::jpnadx=500      !     jpnadx = maximum number of administration records 250
  integer,save::admin(jpnkey,jpnadx,jpnfx)

  integer,parameter::jplrec=3000     !     jplrec = maximum record length in direct access file (must not be less than the number of fields in any one )
!  integer,parameter::jpnrx=50
  integer,parameter::jpnrx=250
  integer,save::mrec(jplrec*jpnrx)   !     mrec   = array to store one record from the file
  logical,save:: ltest(jplrec)!     ltest  = array with test logicals 



!  integer,parameter::jplg=128000      !     jplg   = maximum length of a grib message
!  integer,parameter::jplg=150000      !     jplg   = maximum length of a grib message
!  integer,parameter::jplg=256000      !     jplg   = maximum length of a grib message
  integer,parameter::jplg=600000      !     jplg   = maximum length of a grib message
  integer,save::mg(jplg)             !     mg     = array to store grib message



  integer,parameter::jpbpwx=64 !jpbpwx  must be >=nbpw
  integer,public,save::notdef
  integer,save::nbpw
  integer,save::nbpwio
  integer,save::mask(jpbpwx+1)

  public asimhw,asimhm,asimhc,asimho,asimhr,asimhk,loadfd,getfd,putfd

contains

  subroutine vdate(ydate)
    implicit  none
    character(len=8):: ydate
    ydate='19980101'
  end subroutine vdate


  subroutine vclock(ytime)
    implicit  none
    character(len=6):: ytime
    ytime='120000'
  end subroutine vclock


  subroutine swab(mrec, length)
    implicit none
    integer i, length
    integer mrec(length)
    do i = 1, length
       mrec(i) = ior(           ishft(mrec(i), -24),     &
            ior(ishft(iand(ishft(mrec(i), -16), 255),  8), &
            ior(ishft(iand(ishft(mrec(i),  -8), 255), 16), &
            ishft(iand(      mrec(i),       255), 24))))
    enddo
  end subroutine swab


  subroutine putfd(kdev,kb1,klb1,kb2,klb2,pb2,klpb2,field,fieldSize,pacc,kerr)

    implicit none
    integer,intent(in):: kdev
    integer,intent(in)::klb1,klb2,klpb2,fieldSize
    integer::kerr
    real(kind=realgribkind):: pacc
    integer kb1(klb1),kb2(klb2)
    real(kind=realgribkind) pb2(klpb2),field(fieldSize)
    !     
    !     putfd:  code a field into grib and place it in an asimof database
    !     
    !     parameters, usage (i,o,io)
    !     
    !     kdev:  i unit number of the database. must have been opened by loadfd
    !     kb1:   i to be used as block 1 and as database key
    !     klb1:  i length of kbk1 (at least 30)
    !     kb2:   i block 2 integer data
    !     klb2:  i length of kbk2
    !     pb2:   i block 2 real data
    !     klpb2: i length of pbk2
    !     field:     i field data
    !     fieldSize:     i number of data to be encoded
    !     fieldSize=-1:  release file from list of files opened by loadfd
    !     fieldSize=-2:  mark fields with these keys as deletable
    !     pacc:  i required accuracy for grib encoding
    !     >0.  : accuracy in units of the data
    !     <0.  : negative of number of bits to be used
    !     =0.  : default accuracy (12 bits)
    !     so pacc allows to prescribe absolute accuracy (pacc>0)
    !     or relative accuracy (pacc<0 or pacc=0).
    !     it is the user's responsibility to ensure that the
    !     computer words are big enough to encode the requested
    !     accuracy.
    !     kerr:  o error return code
    !     as returned by engrib or
    !     as returned by asimhw (incremented by 1000)

    !     0:  no errors
    !     1000: no errors....
    !     
    !     method:
    !     engrib is called for grib encoding
    !     asimhw is called for database access
    !     the facilities to delete and/or release files as offered by asimhw
    !     are available to putfd as well: use fieldSize <0

    integer,parameter::jplb3=0,jplb4=11
    integer ib4(jplb4)
    integer kleng,i12,inb,imode
    real(kind=realgribkind) zacc

    !     engrib if the number of field data is positive
    if(fieldSize>0)then
       if(pacc>0._realgribkind)then
          zacc=pacc
          inb=0
       else
          inb=-nint(pacc)
          if(abs(pacc)<1.e-12_realgribkind)inb=12
          zacc=0._realgribkind
       endif
       kleng=jplg
       call engrib(mg,kleng,field,fieldSize, &
            kb1,kb2,pb2,ib4, &
            zacc,inb,real(notdef,realgribkind),kerr)
       if(kerr==210)then

          call stop_program('mg too small in putfd - recompile with bigger jplg')
       endif
       if(kerr/=0)then
          print *,__FILE__,__LINE__,'kerr = ',kerr
          return
       endif
    else
       kleng=fieldSize
    endif

    !     asimhw (putfd stores grib files - mark kb1(2) such)
    i12=kb1(2)
    kb1(2)=0
    imode=0
    call asimhw(kdev,mg,kleng,kb1,imode,kerr)
    kb1(2)=i12
    if(kerr/=0)then
       kerr=kerr+1000
       print *,__FILE__,__LINE__,'kerr = ',kerr
    endif
  end subroutine putfd


  subroutine loadfd(kdev,filename,kbks1,klbk1,knfldx,knfld, &
       lapab,pab,klevx,klev,kmode,kerr)


    !     function:  loadfd opens a direct-access file as grib (wmo fm-92)
    !     database.  if the file exists, tables of contents will
    !     be extracted; else the file will be opened with an
    !     empty table of contents
    !     up to 5 files may be used at any one time
    !     if the file is to be created, an auxiliary array
    !     containing a grib message with nonsense geometry,
    !     and with the complete list of vertical level
    !     parameters can optionally be added to the database
    !     loadfd will not mention this auxilliary data in the
    !     table of contents
    !     
    !     description of parameters (name, usage (i/o/io), description):
    !     
    !     file definition parameters:
    !     
    !     kdev:    i fortran unit number
    !     filename:  i file name
    !     if filename=' ', the unit has been connected already
    !     
    !     table of available fields (for existing files only):
    !     
    !     kbks1:   o data base keys - usually corresponding to the first

    !     30 words of block1 of the grib message
    !     a single key occupies min(30,klblk1) words, and
    !     kbks1(i,j) is the i'th word of the key for the j'th
    !     field
    !     klblk1:  i first dimension of kbks1 (must be at least 30)
    !     knfldx:  i second dimension of kbks1, maximum nr of fields to be
    !     listed in kbks1
    !     knfld:   o number of fields found in the data base
    !     
    !     table of vertical level parameters (existing and new files):
    !     
    !     lapab:   i if true: for new files, write table of parameters
    !     for existing files, reconstruct the table
    !     (the overhead is significant, unless the
    !     table was written during file construction)
    !     pab:    io table of parameters, to be written to new files,
    !     to be reconstructed from existing files
    !     klevx:  i second dimension of pab (first dimension is 2)
    !     klev:  io number of levels
    !     
    !     error reporting:
    !     
    !     kerr:            0 = no errors
    !     <0   as returned by asimho
    !     20 = open inputfile error
    !     30 = read inputfile error
    !     40 = file is empty
    !     31 = more than knfldx fields in file. return only
    !     the first knfldx.
    !     51 = klbk1 < 30
    !     70 = decode error
    !     101 = more than klevx levels in file.  return only
    !     the first klevx
    !     102 = not enough vertical level parameters to
    !     define the level of some field in the file
    !     >200 : error reported by getfd or putfd, incremented
    !     by 10000
    !     
    !     the routine aborts if internal arrays have too small sizes
    !     
    !     
    !     k.n.m.i. 20/03/91
    !     g. cats
    !     r.m. van westrhenen                update 08/07/91
    !     
    !     i,j    = loop variable
    !     inlev  = number of levels in the current field
    !     jlev   = loop index, levels
    !     iu     = cardinal number of unit number
    !     ierr   = error code (0 = no error)
    !     zsup   = supplementary data from getfd

    implicit none
    integer,intent(in)::kdev
    integer::klbk1,knfldx,knfld,klevx,klev,kerr
    integer kmode
    integer kbks1(klbk1,knfldx)
    real(kind=realgribkind) pab(2,klevx)
    logical lapab
    character(len=*) filename
    integer   i,j,iu
    integer   inlev,jlev
    !     ib1..ib2 = definition blocks.
    !     jplb1..jplb2 = lenght of ib1..ib2
    integer,parameter::jplb1=40, jplb2=400,jplsup=2
    integer  ierr
    integer  ib1(jplb1), ib2(jplb2)
    real(kind=realgribkind)     zduma(4),zrwb2(jplb2),zlev,zsup(jplsup)
    integer::dummi
    data zduma/4*0._realgribkind/
    logical found

    if(klbk1<jpnkey)then
       kerr=51
       return
    endif
    knfld =    0
    kerr=0

    !     open file and treat error messages
    call asimho(kdev,filename,kmode,kerr)

    if(kerr/=0)then
       if(kerr==-4.or.kerr==-5)then
          write(6,'(a,i2,a,i3)')' loadfd: error ', &
               kerr,' from routine asimho for unit',kdev
          if(kerr==-4)then
             write(6,'(a,i5,a)')'jpnadx=',jpnadx, &
                  ' in common casihl is too small'
          else
             write(6,'(a,i5,a)')'jpnfdx=',jpnfdx, &
                  ' in common casihl is too small'
          endif
          call stop_program('error in loadfd after call to asimho')
       endif
       if(kerr>0)kerr=20
       return
    endif

    !     find unit in list of open files
    found = .false.
    j=1
    do while(.not.found.and.j<=jpnfx)
       if( munits(j)==kdev )then 
          iu=j
          found=.true.
       endif
       j=j+1
    enddo

    !     do j=1,jpnfx
    !        if( munits(j)==kdev )then 
    !           iu=j
    !           goto 420
    !        endif
    !     enddo

    if(.not.found)then
       write(6,'('' asimof system internal error 4010 from loadfd'', &
            '', contact system support'')')
       call stop_program('asimof system internal error 4010 from loadfd')
       return
    endif


    !     copy blocks without administration blocks to kbks1.

    knfld=0
    do j=1,nfld(iu)
       if(blocks(9,j,iu)/=254.or.blocks(10,j,iu)/=109 &
            .or.blocks(11,j,iu)/=65535)then
          if(knfld>=knfldx)then
             kerr=31
             goto 603
          endif
          knfld=knfld+1
          do i=1,min(klbk1,jpnkey)
             kbks1(i,knfld)=blocks(i,j,iu)
          enddo
       endif
    enddo
603 continue

    !     if lapab then fill array pab
    if (.not.lapab) then
       return
    endif

    !     new files:  only add if the number of levels is positive

    if(nfld(iu)==0.and.klev<=0)goto 760

    !     create/read-if-present administration field

    do  j=1,jplb1
       ib1(j)=0
    enddo

    ib1(1)=24
    ib1(4)=1
    ib1(5)=96
    ib1(6)=1
    ib1(7)=255
    ib1(8)=128
    ib1(9)=254
    ib1(10)=109
    ib1(11)=65535
    ib1(13)=1919
    ib1(14)=1
    ib1(15)=1
    ib1(16)=1
    ib1(17)=1
    ib1(18)=1
    ib1(19)=1
    do  j=1,jplb2
       ib2(j)=0
    enddo
    ib2(1)=32+klev*8+ib2(40)*10
    ib2(7)=2
    ib2(9)=2
    ib2(11)= 1
    ib2(14)=-1
    ib2(17)=128
    ib2(24)= 1
    ib2(26)= 1
    ib2(28)=64
    if(nfld(iu)==0)then
    else
       ib1(8)=notdef
       call getfd(kdev,ib1,jplb1,ib2,jplb2,pab,2*klevx,zsup,jplsup, &
            zduma,0,1,ierr,dummi)
       if(ierr==-1)goto 720
       if(ierr/=0)then
          kerr=10000+ierr
          return
       endif
       klev=(ib2(1)-32-ib2(40)*10)/8
    endif
    goto 760

    !     initialise

720 continue
    i=0
    klev=0

    !     scan through the file, accepting all keys

730 continue
    do j=1,jplb1
       ib1(j)=notdef
    enddo
    i=i+1
    !     next field - blocks 1 and 2 only

    ib2(40)=0
    call getfd(kdev,ib1,jplb1,ib2,jplb2,zrwb2,jplb2,zsup,jplsup, &
         zduma,0,i,ierr,dummi)
    if(ierr==-1)goto 760

    !     if getfd failed, try next field

    if(i<nfld(iu).and.ierr<0)goto 730

    !     skip administration blocks

    if(ib1(10)==109.and.ib1(11)==65535)goto 730
    inlev=(ib2(1)-32-ib2(40)*10)/8
    if(inlev>1)then

       !     this field contains the complete level list (ahalf, bhalf)
       !     -- i only understand hybrid levels here --
       if(ib1(10)/=109)goto 750
       klev=min(klevx,inlev-1)
       if(klev<inlev-1)kerr=101
       do  jlev=1,klev
          pab(1,jlev)=0.5_realgribkind*(zrwb2(jlev+ib2(40))+zrwb2(jlev+ib2(40)+1))
          pab(2,jlev)=0.5_realgribkind*(zrwb2(jlev+ib2(40)+inlev)+ &
               zrwb2(jlev+ib2(40)+inlev+1))
       enddo
       goto 760
    else

       !     i only understand pressure, sigma or hybrid levels here --
       if(ib1(10)==100.or.ib1(10)==107)then
          if(ib1(10)==100)then
             zlev=real(ib1(11),realgribkind)*100._realgribkind
          else
             zlev=real(ib1(11),realgribkind)/10000._realgribkind
          endif
          do  jlev=1,klev
             if (ib1(10) == 100 .and. abs(pab(1, jlev)- zlev)<1.e-14_realgribkind &
                  .and. abs(pab(2, jlev))< 1.e-14_realgribkind)  goto 750
             if (ib1(10) /= 100 .and. abs(pab(1, jlev))<1.e-14_realgribkind &
                  .and. abs(pab(2, jlev) - zlev)<1.e-14_realgribkind) goto 750
          enddo
          if(klev==klevx)then
             kerr=101
             goto 760
          endif
          klev=klev+1
          if (ib1(10) == 100) then
             pab(1, klev) = zlev
             pab(2, klev) = 0._realgribkind
          else
             pab(1, klev) = 0._realgribkind
             pab(2, klev) = zlev
          endif
       elseif(ib1(10)==109.or.ib1(10)==0)then
          if(ib1(11)==0)then
             do  jlev=1,klev
                if( abs(pab(1,jlev)-zrwb2(ib2(40)+1) )<1.e-14_realgribkind.and.&
                     abs(pab(2,jlev)-zrwb2(ib2(40)+2))<1.e-14_realgribkind) goto 750
             enddo
             if(klev>=klevx)then
                kerr=101
                goto 760
             endif
             klev=klev+1
             pab(1,klev)=zrwb2(ib2(40)+1)
             pab(2,klev)=zrwb2(ib2(40)+2)
          else
             inlev=ib1(11)
             if(inlev>klevx)then
                kerr=101
             else
                klev=max(inlev,klev)
                pab(1,inlev)=zrwb2(ib2(40)+1)
                pab(2,inlev)=zrwb2(ib2(40)+2)
             endif
          endif
       endif
    endif

    !     next field

750 continue
    if(i<nfld(iu).and.ierr>0)goto 730

760 continue
  end subroutine loadfd


  subroutine int2char(y,k,kc)
    !     int2char: read kc characters from k into y
    implicit none
    character(len=*)::y
    integer::k,kc,i,j

    i=k
    do j=kc,1,-1
       y(j:j)=char(iand(255,i))
       i=ishft(i,-8)
    enddo
  end subroutine int2char


  subroutine getfd(kdev,kb1,klb1,kb2,klb2,pb2,klpb2,psup,klpsup, &
       buf,lenrec,kpf,kerr,retlen)

    !     getfd:  code a field into grib and place it in an asimof database
    !     
    !     parameters, usage (i,o,io)

    !     kdev:  i unit number of the database. must have been opened by loadfd
    !     kb1:  io to be used as block 1 and as database key
    !     on input:  key to be matched
    !     on output: matching key
    !     see asimhr
    !     klb1:  i length of kbk1 (at least 30)
    !     kb2:   o block 2 integer data
    !     klb2:  i length of kbk2
    !     pb2:   o block 2 real data
    !     klpb2: i length of pbk2
    !     psup:  o supplementary data
    !     psup(1):  real(number of bits used for encoding)
    !     psup(2):  accuracy used for encoding
    !     psup(3):  real(mode of gribcode operation)
    !     klpsup:i length of psup
    !     p:     o field data
    !     k:    io number of data to be decoded
    !     on input:  length of array
    !     k=0 on input: decode for blocks 1 and 2 only
    !     on output: number of data used
    !     kpf:   i deliver kpf'th matching field
    !     kerr:  o error return code
    !     as returned by asimhr or
    !     so: >=0: number of remaining matching fields
    !     -2000: no errors....
    !     
    !     method:
    !     asimhr is called for database access

    !     
    !     gerard cats  18 july 1991


 

    implicit none
    integer,intent(in)::lenrec
    integer,intent(inout)::retlen
    integer::kdev,klb1,klb2,klpb2,kpf,kerr,klpsup
    integer kb1(klb1),kb2(klb2)
    real(kind=realgribkind) pb2(klpb2),buf(lenrec),psup(klpsup)

    integer, parameter::jplb3=0,jplb4=11
    integer ib4(jplb4)
    integer kleng,i12
    !  asimhr (getfd requests grib files - mark kb1(2) such)

    i12=kb1(2)
    kb1(2)=0
    kerr=kpf
    kleng=jplg
    call asimhr(kdev,mg,kleng,kb1,kerr)
    if(kerr==-4)then
       call stop_program('mg too small in getfd - recompile with bigger jplg')
    endif
    if(kerr<0)then
       kb1(2)=i12
       return
    endif


    if(klpsup>=3)call grbmod(nint(psup(3)),.true.)
    call degrib(mg,kleng,buf,lenrec,kb1,kb2,pb2,ib4,real(notdef,realgribkind),i12,retlen)
    if(i12/=0)then
       kerr=i12-2000
       return
    endif
    if(lenrec>0)then
       if(klpsup>=1)psup(1)=real(ib4(11),realgribkind)
       if(klpsup>=2)psup(2)=2._realgribkind**ib4(5)
    endif
  end subroutine getfd


  subroutine char2int(y,k,kc)
    !     char2int: read kc characters from y into k
    implicit none
    character(len=*)::y
    integer::k,kc,j,i

    k=0
    i=0
    do j=kc,1,-1
       k=ior(k,ishft(ichar(y(j:j)),i))
       i=i+8
    enddo
  end subroutine char2int

  subroutine asimhw( kdev, kfld, klenf, kb1, kmode, kretco)

    implicit none
    integer kdev,klenf,kmode,kretco
    integer kfld(klenf),kb1(30)

    !     #!extract:on 3
    !     asimhw:  write an array into the database
    !     
    !     kdev:    i unit number of the file
    !     kfld:    i array to be written
    !     klenf:   i number of words to be written
    !     -1:  reelase unit number
    !     -2:  delete info for key kb1 (all matching info will
    !     be deleted - so beware of non-unique key
    !     specification).  asimhw does not provide
    !     the facility to reuse deleted space; deleted
    !     info can therefore always be retrieved, use
    !     kb1(1)=9999 to access deleted info
    !     kb1:     i key
    !     kmode:   i overwrite mode

    !     1:  allow the old array associated with the first
    !     key matching kb1 to be overwritten by the new
    !     data array when its length does not exceed
    !     that of the old array, or mark the old array
    !     as deleted and add the new array
    !     other:  if klenf/=-2 a key matching kb1 will cause
    !     an error and the data array will not be
    !     written to file
    !     kretco:  o error return code
    !     0 no errors
    !     1 kdev not opened yet by loadfd
    !     2 writing to this file disabled
    !     3 key is not unique -storage failed (klenf>=0)
    !     or key not found (klenf=-2)
    !     4 error when rewriting old adminstration record
    !     5 no space to add administration-increase jpnadx
    !     6 no space to add array to admin-increase jpnfdx
    !     7 error when writing data record
    !     8 error when writing updated administration rec
    !     9 error when reading part of last record written
    !     
    !     #!extract:off 3
    !     #!extract:on 2
    !     description of an asimof file:
    !     nb in the following comment, all (record)lengths, (word) pointers
    !     and 'words' are in units of nbpwio (currently 32) bits.
    !     an asimof file is a data base structure, which is defined in a
    !     computer independent format.  it consists of administration blocks
    !     which contain positive integers (last bit least significant), and
    !     user data, without any record structure.  on most computers,
    !     the file is easily accessed as a fortran random-access file, with
    !     record lengths of jplrec (currently 3000) units.  on big-endian
    !     machines, it is advised to swap bytes before (when reading) or
    !     after (when writing) access to the file. in the following, a
    !     '(physical) record' refers to a group of jplrec units.
    !     the above applies to version 1.01.  in an earlier
    !     definition, the units were computer words, and hence the format
    !     depended on the computer word length.  the asimof routines to
    !     read a file are upward compatible with earlier formats, but files
    !     can only be written in the new format.  the write routines allow
    !     addition to an existing file, but only if the existing file is also
    !     in the new format.
    !     the first record is an administration block.  more administration
    !     blocks may be present in the file.  the data are stored in 'arrays'
    !     (or 'logical' records), that may extend over several physical records.
    !     an administration block consists of 30 words general administration
    !     (array 'admin') and 99 repetitions of 30 words, each relevant to one
    !     array or 'logical' record (array 'blocks').  the following are
    !     the contents of admin and blocks:
    !     admin(1) :  number of arrays administrated in this block (0...99)
    !     admin(3) :  last record referenced by this block
    !     admin(4) :  last word used in last record in this block
    !     admin(5) :  jplrec of producing process           |
    !     admin(12):  file version nr. currently 1.02 (=258)|
    !     admin(13): | admin(13) and admin(14) contain      | only in 1st block
    !     admin(14): | experiment code (format 2a4)         |
    !     admin(15):  1 if this block was changed, else 0
    !     admin(16):  creation date (format yyyymmdd) (0=missing) |   "        "
    !     admin(17):  creation time (format hhmmss) (0=missing if admin(17)=0)|"
    !     admin(30):  record number of next administration block
    !     blocks(1 until 27): data base key of the array
    !     blocks(1) : length of key actually used
    !     9999 if record is to be considered 'deleted'
    !     blocks(2) : code type:  0 gribcode
    !     blocks(2)....blocks(blocks(1)) must define record uniquely
    !     
    !     blocks(28): record number of first word of the array
    !     blocks(29): start word of the array in the first record
    !     blocks(30): length of the array
    !     #!extract:off 2
    !     
    !     gerard cats  17 july 1991
    !     modifications:
    !     gerard cats, 7 september 1994:  introduce version 1.01
    !     
    !     chiel van zijl, 26 augustus 1998: clear unused part of
    !     record mrec(ifr..ito). this makes it possible to use the routine
    !     pbgrib (reads grib messages from unformatted files).
    !     
    logical lofirs,lotst,lover
    integer irlen,iu,icount,ikrt,inj,iblock,inadm,inr,irlast,ilenf
    integer iptr,ib,ibs,i4,i5,ilx,ix,ibpwf
    integer j,jk
    logical found
    !     1.  initialise
    !     
    !     1.1 several counters
    !     
    ikrt=0
    inr=1
    icount=0
    iblock=0

    !     1.2 occurrence of matching key
    !     
    lofirs=.false.
    lotst=.false.
    lover=.false.

    !     1.9 return error code
    !     
    kretco=0
    !     
    !     2.  find unit in list of open files, release unit if requested
    !     
    found = .false.
    j=1
    do while(.not.found.and.j<=jpnfx)
       if( munits(j)==kdev )then
          irlen = nlrec(j)
          ibpwf=nbpwf(j)
          inadm= jplrec/jpnkey-1
          iu=j
          found = .true.
       endif
       j=j+1
    enddo

    if(.not.found)then
       kretco = 1
       return
    endif

    !     unit found - release it if requested
    if(klenf==-1)then
       munits(iu)=-1
       nfld(iu)=0
       return
    endif

    !     check file consistency

    if(ibpwf/=nbpwio)then
       print *,'asimhw: cannot write to this file:'
       write(0,'(a,i3,a,i2,a,i2)')'unit=',kdev,'; bits per word (', &
            ibpwf,') is not ',nbpwio
       call stop_program('cannot write file')
    endif

    !     check write access

    if(mopen(iu)==0)then
       kretco=2
       return
    endif

    !       loop over administration blocks

    !      length of array in io word count

    if(klenf>=0)ilenf=(klenf*nbpw+nbpwio-1)/nbpwio

    !      next block

300 continue

    !     look for matching keys

    inj=admin(1,icount+1,iu)
    if(inj==0)goto 390
    do j=1,inj
       ltest(j)=.true.
    enddo
    if(kb1(1)/=notdef)then
       ilx=min(kb1(1),27)
    else
       ilx=27
    endif
    do jk=1,ilx
       if (kb1(jk)/=notdef)then
          do  j=1,inj
             ltest(j)=ltest(j).and. &
                  blocks(jk,j+iblock,iu)==kb1(jk).and. &
                  blocks(1,j+iblock,iu)/=9999
          enddo
       endif
    enddo

    !     3.2 if info is to be deleted, mark key and block, rewrite block
    !     
    if(klenf==-2.or.kmode==1)then
       do j=1,inj
          if(ltest(j))then
             lover=ilenf<=blocks(30,j+iblock,iu).and.klenf>=0
             if(.not.lover)then
                blocks(1,j+iblock,iu)=9999
             else
                if(ilenf<blocks(30,j+iblock,iu).and. &
                     mopen(iu)/=-890) admin(15,icount,iu)=1
                blocks(30,j+iblock,iu)=ilenf
                ibs=blocks(28,j+iblock,iu)
                i4=blocks(29,j+iblock,iu)-1
             endif
             lotst=.true.
          endif
       enddo
       if(lotst)then
          lofirs=.true.
          if(mopen(iu)==-890)then
             call asimhwa(kdev,iu,inr,icount,inadm, kretco)
          endif
          if(kretco/=0)then
             kretco=4
             return
          endif
          if(lover)then
             goto 521
          endif
       endif
    else

       !      may info be added ? - is key unique?
       !     
       do j=1,inj
          if(ltest(j)) ikrt=ikrt+1
       enddo
       if (ikrt/=0) then
          kretco=3
          return
       endif
    endif

    !      next administration block (if present) - stop if delete fnctn
    !     
390 continue
    irlast=inr
    icount=icount+1
    iblock=iblock+inadm
    inr=admin(30,icount,iu)
    if(inr/=0)goto 300
    !     
    if(klenf==-2)then
       if(.not.lofirs)kretco=3
       return
    endif

    !       actions if last administration block is full
    !     
    !      check if full
    !     
    if (admin(1,icount,iu)<inadm)goto 500
    !     
    !     4.2 rewrite administration with pointer to newly added block
    !     
    iptr=admin(3,icount,iu)+1
    admin(30,icount,iu)=iptr
    if(mopen(iu)==-890)then
       call asimhwa(kdev,iu,irlast,icount-1,inadm, kretco)
    else
       admin(15,icount,iu)=1
    endif
    if(kretco/=0)then
       kretco=4
       return
    endif

    !     initiate new administration record

    icount=icount+1
    if(icount>jpnadx)then
       print *,'icount',icount,'>',jpnadx,'jpnadx',__FILE__,__LINE__
       kretco=5
       return
    endif
    do j=1,jpnkey
       admin(j,icount,iu)=0
    enddo
    admin(3,icount,iu)=iptr+1
    admin(4,icount,iu)=0
    irlast=iptr

    !       write array

500 continue

    !      if there is place in the administration blocks, add array

    if(nfld(iu)>=jpnfdx)then
       print *, nfld(iu) , 'nfld(iu) > ',jpnfdx, 'jpnfdx' ,__FILE__,__LINE__
       kretco=6
       return
    endif
    nfld(iu)=nfld(iu)+1
    do  j=1,min(27,kb1(1))
       blocks(j,nfld(iu),iu)=kb1(j)
    enddo

    !      determine place where to store the array
    !     
    ibs=admin(3,icount,iu)
    i4=admin(4,icount,iu)
521 continue
    ib=ibs

    !      construct mrec
    !     
    !     ix: number of bits this record
    !     ilx: total number of bits to be read
    !     i5: last word read into mrec
    !     
    ix=irlen*nbpw-i4*nbpwio
    ilx=klenf*nbpw
    i5=0
    !     
    !     initialise first record to be overwritten
    if(lover.or.i4/=0)then
       read(kdev,rec=ib,iostat=kretco)(mrec(j),j=1,irlen)
       if(kretco /= 0) then
          kretco=9
          return
       endif
#ifdef LITTLE_ENDIAN
       call swab(mrec(1), irlen)
#endif
    else
       do j=1,irlen
          mrec(j)=0
       enddo
    endif

    !     count records to be (over)written
    !     
531 continue
    if(ix>=ilx)goto 532
    i5=i5+irlen
    if(i5+irlen>jplrec*jpnrx)then
       call stop_program('mrec too small to contain record increase jpnrx in common casihl in asimof')
    endif
    ib=ib+1
    ix=ix+irlen*nbpw
    goto 531
532 continue

    !     when overwriting, also initialise last record
    !     
    if(ib>ibs)then
       if(lover)then
          read(kdev,rec=ib,iostat=kretco)(mrec(j),j=i5+1,i5+irlen)
          if(kretco /= 0) then
             kretco=9
             return
          endif
#ifdef LITTLE_ENDIAN
          call swab(mrec(i5+1), irlen)
#endif
       else
          do j=1,irlen
             mrec(i5+j)=0
          enddo
       endif
    endif

    !     transfer data from kfld into mrec
    !     
    ix=i4*nbpwio
    call gsbite(mrec,kfld,ix,nbpw,0,klenf,nbpw,mask,'e')
    !     
    !     5.4 (re)write records
    !     
    i5=0
    do jk=ibs,ib
#ifdef LITTLE_ENDIAN
       call swab(mrec(i5+1), irlen)
#endif
       write(kdev,rec=jk,iostat=kretco)(mrec(j+i5),j=1,irlen)
#ifdef LITTLE_ENDIAN
       call swab(mrec(i5+1), irlen)
#endif
       if(kretco /= 0) then
          kretco=7
          return
       endif
       i5=i5+irlen
    enddo

    !     5.8 set start record and word in the administration, next record
    !     
    if(lover)return
    blocks(28,nfld(iu),iu)=ibs
    blocks(29,nfld(iu),iu)=admin(4,icount,iu)+1


    !     
    !     6.  complete adminstration block and write it out
    !     
    ix=ix+(ibs-ib)*nbpwio*jplrec
    admin(4,icount,iu)=(ix+nbpwio-1)/nbpwio
    admin(3,icount,iu)=ib
    admin(1,icount,iu)=admin(1,icount,iu)+1
    blocks(30,nfld(iu),iu)=ilenf
    if(mopen(iu)==-890)then
       call asimhwa(kdev,iu,irlast,icount-1,inadm, kretco)
    else
       admin(15,icount,iu)=1
    endif
    if (kretco /= 0) then
       kretco=8
       return
    endif
  end subroutine asimhw


  subroutine asimhwa(kdev,ku,krecnr,kcount,kadm, kioerr)

    !     asihmwa:  (re)write administration record

    !     input:
    !     kdev:   unit nr
    !     ku:     unit nr within ASIMOF system
    !     krecnr: record number
    !     kcount: administration record to be (re)written
    !     kadm:   nr of blocks witihn one administtartion record

    !     output:
    !     kioerr: iostat error indicator

    !     Gerard Cats, KNMI, 6 September 1994


    implicit none
    integer  kdev, ku, krecnr, kcount, kadm, kioerr


    integer  j, jk, iskst, iu, icount, iadm
    iu     = ku
    icount = kcount
    iadm   = kadm
    if(nbpw>nbpwio)then

       !     pack data before writing
       iskst=0
       call gsbite(mrec,admin(1,icount+1,iu),iskst,nbpwio,0,&
            jpnkey,nbpw,mask,'E')
       call gsbite(mrec,blocks(1,1+icount*iadm,iu),iskst,nbpwio,0,&
            jpnkey*iadm,nbpw,mask,'E')
       iskst=(iskst+nbpw-1)/nbpw

#ifdef LITTLE_ENDIAN
       call swab(mrec(1), iskst)
#endif
       write(kdev,rec=krecnr,iostat=kioerr)(mrec(j),j=1,iskst)

#ifdef LITTLE_ENDIAN
       call swab(mrec(1), iskst)
#endif

    else

#ifdef LITTLE_ENDIAN
       call swab(admin(1,icount+1,iu), jpnkey)
       do j = 1+icount*iadm, (icount+1)*iadm
          call swab(blocks(1,j,iu), jpnkey)
       enddo
#endif

       write(kdev,rec=krecnr,iostat=kioerr)     &
            (admin(jk,icount+1,iu),jk=1,jpnkey),&
            ((blocks(jk,j,iu),jk=1,jpnkey),&
            j=1+icount*iadm,(icount+1)*iadm)

#ifdef LITTLE_ENDIAN
       call swab(admin(1,icount+1,iu), jpnkey)
       do j = 1+icount*iadm, (icount+1)*iadm
          call swab(blocks(1,j,iu), jpnkey)
       enddo
#endif

    endif
  end subroutine asimhwa


  subroutine asimhr( kdev, kfld, klenf, kb1, kerr )

    implicit none
    integer kdev,klenf,kerr
    integer kfld(klenf),kb1(30)

    !     asimhr:  write an array into the database
    !     
    !     kdev:    i unit number of the file
    !     kfld:    i array to be filled
    !     klenf:  io on input: length of kfld
    !     on output:actual number of data transferred
    !     kb1:    io key
    !     on input: key to be matched
    !     on output:key found
    !     kb1 may contain missing data indicator 'notdef'
    !     on input - these will always match.  this routine
    !     replaces those missing data by the actual values
    !     kerr: io
    !     on input:  number of matching field to be extracted
    !     <=0: first matching field
    !     on output: error return code
    !     -1 no matching key found
    !     -2 kdev not opened yet by loadfd
    !     -3 error when reading the data
    !     -4 klenf too small to contain the array
    !     -5 file does not contain a single field
    !     >=0: number of fields still matching
    !     so 0 means last matching field succesfully read
    !     
    !     gerard cats  17 july 1991
    !-----------------------------------------------------------------------


    integer irlen,iu,icount,ikrt,inj,ilenf,ierr
    integer ib,i4,i5,ilx,ix
    integer ib1
    integer j,jk
    logical found

    ikrt = 0
    icount = 0

    !     find unit in list of open files
    found = .false.
    j=1
    do while(.not.found.and.j<=jpnfx)
       if( munits(j)==kdev )then
          irlen = nlrec(j)
          iu=j
          found = .true.
       endif
       j = j+1
    enddo

    if(.not.found)then
       kerr=-2
       return
    endif

    !     find matching fields

    !     length of array

    if(klenf>=0)then
       ilenf=klenf
    else
       kerr=-4
       return
    endif

    !     look for matching keys

    inj=nfld(iu)
    if(inj==0)then
       kerr=-5
       return
    endif
    if(inj>jplrec)then
       print *,'jplrec less than the number of fields'
       print *,'jplrec:',jplrec,',nfld:',nfld(iu)
       call stop_program('jplrec less than the number of fields')
    endif

    do j=1,inj
       ltest(j)=.true.
    enddo

    if(kb1(1)/=notdef)then
       ilx=min(kb1(1),27)
    else
       ilx=27
    endif

    do jk=1,ilx
       if(kb1(jk)/=notdef)then 
          ib1=kb1(jk)
          if(jk==13)then
             if(admin(12,1,iu)<258)then
                if(ib1<100)ib1=ib1+1900
             endif
          endif
          do j=1,inj
             ltest(j)=ltest(j).and.blocks(jk,j,iu)==ib1
          enddo
       endif
    enddo

    !     reject match on deleted field, unless explicitly requested

    if(kb1(1)==notdef)then
       do j=1,inj
          ltest(j)=ltest(j).and.(blocks(1,j,iu)/=9999)
       enddo
    endif

    !     3.3 number of matching keys, and index of matching key

    if(kerr<=0)kerr=1
    do j=1,inj
       if (ltest(j))then
          ikrt=ikrt+1
          if(ikrt==kerr)icount=j
       endif
    enddo
    if(ikrt<kerr)then
       kerr=-1
       return
    endif
    kerr=ikrt-kerr

    !     complete key

    do jk=1,27
       kb1(jk)=blocks(jk,icount,iu)
    enddo

    !     read

    !     determine place where to find the array

    ib=blocks(28,icount,iu)
    i4=blocks(29,icount,iu)
    ilx=blocks(30,icount,iu)

    if(ilx*nbpwio>ilenf*nbpw)then
       kerr=-4
       ilx=ilenf*nbpw/nbpwio
    endif
    klenf=0

    !     read records

    !     ix: number of bits this record
    !     ilx: total number of bits to be read
    !     i5: last word read into mrec

    ix=irlen*nbpw-(i4-1)*nbpwio
    ilx=ilx*nbpwio
    i5=0
531 continue
    read(kdev,rec=ib,iostat=ierr)(mrec(j),j=i5+1,i5+irlen)

#ifdef LITTLE_ENDIAN
    call swab(mrec(i5+1), irlen)
#endif
    if(ix<ilx)then
       i5=i5+irlen
       if(i5+irlen>jplrec*jpnrx)then
          print *,'mrec too small to contain record,', &
               ' increase jpnrx in common casihl in asimof'
          call stop_program('mrec too small to contain record, increase jpnrx in common casihl in asimof')
       endif
       ib=ib+1
       ix=ix+irlen*nbpw
       goto 531
    endif
    !     extract the data from the bitstream

    ix = (i4-1)*nbpwio
    ilenf = (ilx+nbpw-1)/nbpw
    call gsbite(mrec,kfld,ix,nbpw,0,ilenf,nbpw,mask,'D')
    klenf = ilenf
  end subroutine asimhr


  subroutine asimho(kdev,filename,kmode,kerr)

    implicit none
    integer,intent(in)::kdev
    integer::kmode,kerr
    character(len=*)::filename

    !     function:  asimho opens an asimof file
    !     up to 5 files may be used at any one time
    !     
    !     description of parameters (name, usage (i/o/io), description):
    !     
    !     file definition parameters:
    !     
    !     kdev:    i fortran unit number
    !     filename:  i file name
    !     if filename=' ', the unit has been connected already
    !     kmode:   i mode of operation:
    !     -890: open by loadfd - file may be open already
    !     -1: don't care whether file exists or not
    !     0: file must exist, else error; write disabled
    !     1: file must not exist, else error
    !     2: if file exists, clear it
    !     3: file must exist, else error; write enabled
    !     
    !     error reporting:
    !     
    !     kerr:            0 = no errors
    !     -1 = no unit number specified
    !     -2 = unit already in use for databases
    !     -3 = maximum number of files already open
    !     -4 = too many administration blocks in file
    !     -5 = too many records in file
    !     postive: io error
    !     
    !     
    !     the routine aborts if internal arrays have too small sizes
    !     
    !     
    !     k.n.m.i. 20/09/94
    !     g. cats

    logical   exist,lofirs
    save      lofirs
    integer   i, j, iu, icount, inadm
    save inadm
    character(len=8)::yostat



    integer  ibpw,ibpwio
    integer  irec(1),iskst,irec1(1)
    integer  ivn
    character(len=8)::ydate
    character(len=6)::ytime

    character(len=256)::cexp
    integer  iexp


    data lofirs/.true./

    !     initialise file contents, error

    kerr = 0

    !     initialise bits per word, mask, etc.

    if(lofirs)then
       lofirs=.false.
       ibpw=0
       ibpwio=0
       call asimhm(ibpw,ibpwio,notdef)
       do 231 j=1,jpnfx
          munits(j)=-1
          mopen(j)=-891
231    enddo
       inadm=jplrec/jpnkey
       if(inadm*jpnkey/=jplrec)then
          print *,'record length not multple of block length'
          call stop_program('record length not multple of block length')
       endif
       inadm=inadm-1
    endif

    !     check if unit exist, and not opened yet, allocate number
    inquire(kdev,exist=exist)
    if (.not. exist) then
       kerr = -1
       return
    endif
    do 301 j=1,jpnfx
       if(munits(j)==kdev)then
          if(kmode==-890)then
             kerr=0
          else
             kerr=-2
          endif
          return
       endif
301 enddo
    do 302 j=1,jpnfx
       if(munits(j)==-1)then
          munits(j)=kdev
          iu=j
          nfld(iu)=0
          nbpwf(iu)=nbpwio
          goto 303
       endif
302 enddo
    kerr=-3
    return
303 continue

    !     if file exist then open file, else make a new file

    !     file name might be 'fort.'unitnr
    if (filename==' ') then
       print *,'fort.* names are not allowed in rca'
       call stop_program('fort.* names are not allowed in rca')
    endif

    !     open new/existing file
    if (kmode==0.or.kmode==3)then
       yostat='old'
    else if (kmode==1) then
       yostat='new'
    else
       yostat='unknown'
    endif
    !      print *,filename,yostat,kmode

    nlrec(iu) = jplrec
#ifdef WORDRECLEN
    open(unit      = kdev,         &
         file      = filename,      &
         status    = yostat,        &
         form      = 'unformatted', &
         access    = 'direct', &
         recl      = nlrec(iu)*nbpwio/nbpw,&
         iostat    = kerr )
#else
    open(unit      = kdev,         &
         file      = filename,      &
         status    = yostat,        &
         form      = 'unformatted', &
         access    = 'direct', &
         recl      = nlrec(iu)*nbpwio/(8*4),&
         iostat    = kerr )
#endif

    if(kerr/=0)return

    !     check if the file exists and is non-empty

    yostat = 'new'
    irec=0
    irec1=0
    if(kmode/=2)then
       read(kdev, rec=1, iostat=kerr) irec,irec1
#ifdef LITTLE_ENDIAN
       call swab(irec,  1)
       call swab(irec1, 1)
#endif
    endif

    !     the file is new; preset administration blocks

    if(kerr/=0.or.(irec(1)==0.and.irec1(1)==0))then
       kerr=0
       do  j=1,jpnkey
          admin(j,1,iu)=0
       enddo
       admin(3,1,iu)=2
       admin(4,1,iu)=0
       admin(5,1,iu)=jplrec

       !     version number 1.2
       admin(12,1,iu)=258

       cexp='rca'

       call char2int(cexp, iexp, 4)
       if( cexp==' ' .or. iexp==0)then
          call char2int('n-a ', admin(13,1,iu), 4)
          call char2int('    ', admin(14,1,iu), 4)
       else
          call char2int(cexp(1:4), admin(13,1,iu), 4)
          call char2int(cexp(5:8), admin(14,1,iu), 4)
       endif
       call vdate(ydate)
       read(ydate,'(i8)') admin(16,1,iu)
       call vclock (ytime)
       read(ytime,'(i6)') admin(17,1,iu)

       nlrec(iu)=nlrec(iu)*nbpwf(iu)/nbpw
    else

       !     open old file

       yostat = 'old'

       !     bits file needs reopening

       if(irec(1)>0.and.nbpw==nbpwio)goto 520
       if(irec(1)>10000)goto 520
       nbpwf(iu)=nbpwf(iu)*nbpw/nbpwio
       close(kdev)
#ifdef WORDRECLEN
       open(unit      = kdev,        &
            file      = filename,     &
            status    = yostat,       &
            form      = 'unformatted',&
            access    = 'direct', &
            recl      = nlrec(iu)*nbpwf(iu)/nbpw, &
            iostat    = kerr )
#else
       open(unit      = kdev,        &
            file      = filename,     &
            status    = yostat,       &
            form      = 'unformatted',&
            access    = 'direct', &
            recl      = nlrec(iu)*nbpwf(iu)/(8*4),&
            iostat    = kerr )
#endif
       if(kerr/=0)then
          kerr=21
          return
       endif

       !     read administration blocks and product definition blocks

520    continue
       nlrec(iu)=nlrec(iu)*nbpwf(iu)/nbpw
       irec(1)   = 1
       icount = 0

521    continue
       if(icount>=jpnadx)then
          kerr=-4
          return
       endif
       if((icount+1)*inadm>jpnfdx)then
          kerr=-5
          return
       endif

       !     treat 64 bits files,

       !     32/32 or 64/64 (file/machine) combinations

       if(nbpwf(iu)==nbpw)then
          read(kdev, rec=irec(1), iostat=kerr) &
               (admin(i,icount+1,iu),i=1,jpnkey), &
               ((blocks(i,j,iu),i=1,jpnkey),j=1+icount*inadm, &
               (icount+1)*inadm)
#ifdef LITTLE_ENDIAN
          call swab(admin(1,icount+1,iu), jpnkey)
          do j = 1+icount*inadm, (icount+1)*inadm
             call swab(blocks(1,j,iu), jpnkey)
          enddo
#endif

          if (kerr/=0) return

          !     32 bits machine, 64 bits file or 64 bits machine, 32 bits file

       else
          read(kdev,rec=irec(1),iostat=kerr)(mrec(j),j=1,nlrec(iu))
#ifdef LITTLE_ENDIAN
          call swab(mrec(1), nlrec(iu))
#endif

          if (kerr/=0) return
          iskst=64-nbpw
          call gsbite(mrec,admin(1,icount+1,iu),iskst,nbpwio,64-nbpw, &
               jpnkey,nbpw,mask,'d')
          call gsbite(mrec,blocks(1,1+icount*inadm,iu),iskst,nbpwio, &
               64-nbpw,jpnkey*inadm,nbpw,mask,'d')
       endif

       !     convert 64 bits pointers and lengths to 32 bits units

       if(nbpwf(iu)/=nbpwio)then
          admin(4,icount+1,iu)=admin(4,icount+1,iu)*nbpwf(iu)/nbpwio
          do j=1+icount*inadm,(icount+1)*inadm
             blocks(29,j,iu)=(blocks(29,j,iu)-1)*nbpwf(iu)/nbpwio+1
             blocks(30,j,iu)=blocks(30,j,iu)*nbpwf(iu)/nbpwio
          enddo
       endif

       !     next adminstration block

       icount = icount+1
       nfld(iu)   = nfld(iu) + admin(1,icount,iu)
       irec(1)   = admin(30,icount,iu)
       if(irec(1)/=0)goto 521
    endif

    !     convert old century blocks to new one
    !     for version 1.01 only; and then it is the same as 1.02

    ivn=admin(12,1,iu)
    if(ivn==257)then
       do j=1,nfld(iu)
          if(blocks(13,j,iu)>100)ivn=258 ! it was already 1.02
       enddo
    endif
    if(ivn==258)admin(12,1,iu)=258
    if(ivn<258)then
       do j=1,nfld(iu)
          if(blocks(13,j,iu)<100)then
             blocks(13,j,iu)=blocks(13,j,iu)+1900
          endif
       enddo
    endif

    mopen(iu)=kmode
  end subroutine asimho


  subroutine asimhm(kbpw,kbpwio,kmdi)


    implicit none
    integer,intent(out)::kbpw,kbpwio
    integer,intent(out)::kmdi
    !     function:  asimhm returns or sets machine constants for asimof
    !     
    !     description of parameters (name, usage (i/o/io), description):
    !     
    !     kbpw:   io number of bits in a computer word
    !     kbpwio: io number of bits in a computer word for io
    !     kmdi:    o missing data indicator
    !     
    !     kbpw and kbpwio will be used to set the values in common cmachl,
    !     but if they were zero on entry, the current values of the
    !     common will be returned, or calculated
    !     kmdi will always be returned, it cannot be set by this routine
    !     
    !     k.n.m.i. 20/03/91
    !     g. cats



    logical lofirs
    integer idum(1)
    integer ilarge
!    real zLarge

    save lofirs
    data lofirs/.true./

    data iLarge/100000000/
    !     in first entry, set or calculate common variables

    if(lofirs)then
       lofirs=.false.

       !     preset nbpw

       if(kbpw>0.and.kbpw<2096)then
          nbpw=kbpw
       else
          nbpw=0
       endif
       !     mask is always calculated, nbpw rides on its back (if 0)

       mask(1)=0
       mask(2)=0
       call gsbite(idum,idum,0,0,0,0,nbpw,mask,'e')
       if(nbpw>jpbpwx)then
          print *,'jpbpwx too small in asimhm - must be >=',nbpw
          call stop_program('jpbpwx too small in asimhm - must be >=')
       endif
       !     bits per word for io

       nbpwio=32

       !     missing data indicator

       !notdef= nint(zLarge)
       notdef = ilarge

       !     end of first execution path

    endif
    !     return values

    kbpw=nbpw
    kbpwio=nbpwio
    kmdi=notdef
  end subroutine asimhm


  subroutine asimhk( kdev, keys, klk1, knkeys, kb1, kretco )

    implicit none
    integer kdev,klk1,knkeys,kretco
    integer keys(klk1,knkeys),kb1(30)

    !     asimhk:  find at most knkeys keys matching kb1

    !     kdev:    i unit number of the file
    !     keys:    i array of matching keys
    !     klk1:    i first dimension of keys
    !     knkeys:  io on input: second dimension of keys
    !     also the maximum nr of keys to be returned
    !     on output: actual number of keys matched
    !     kb1:    i  key to be matched
    !     kb1 may contain missing data indicator 'notdef'
    !     on input - these will always match.
    !     kretco: io
    !     on input:  number of field to be matched
    !     <=0: first matching field
    !     on output: error return code
    !     -1 no matching key found
    !     -2 kdev not opened yet by loadfd/asimho
    !     -3 knkeys <= 0
    !     -5 file does not contain a single field
    !     -6 klk1<27
    !     >=0: number of fields still matching
    !     so 0 means knkeys matching fields
    !     
    !     gerard cats  17 march 1996



    integer iu,icount,ikrt,inj
    integer ilx,im
    integer ib1
    integer j,jk
    logical found

    !     input checks

    if(klk1<27)then
       kretco=-6
       return
    endif
    if(knkeys<0)then
       kretco=-4
       return
    endif


    !     several counters

    ikrt=0
    im=0
    icount=0

    !     find unit in list of open files
    found=.false.
    j=1
    do while(.not.found.and.j<=jpnfx)
       if( munits(j)==kdev )then
          iu=j
          found = .true.
       endif
       j=j+1
    enddo


    if(.not.found)then
       kretco=-2
       return
    endif

    !     find matching fields
    !     length of array
    !     look for matching keys

    inj=nfld(iu)
    if(inj==0)then
       kretco=-5
       return
    endif
    if(inj>jplrec)then
       print *,'jplrec less than the number of fields'
       print *,'jplrec:',jplrec,',nfld:',nfld(iu)
       call stop_program('jplrec less than the number of fields')
    endif
    do  j=1,inj
       ltest(j)=.true.
    enddo
    if(kb1(1)/=notdef)then
       ilx=min(kb1(1),27)
    else
       ilx=27
    endif
    do  jk=1,ilx
       if(kb1(jk)/=notdef)then
          ib1=kb1(jk)
          if(jk==13)then
             if(admin(12,1,iu)<258)then
                if(ib1<100)ib1=ib1+1900
             endif
          endif
          do j=1,inj
             ltest(j)=ltest(j).and.blocks(jk,j,iu)==ib1
          enddo
       endif
    enddo

    !      reject match on deleted field, unless explicitly requested

    if(kb1(1)==notdef)then
       do  j=1,inj
          ltest(j)=ltest(j).and.blocks(1,j,iu)/=9999
       enddo
    endif

    !      number of matching keys

    if(kretco<=0)kretco=1
    do j=1,inj
       if(ltest(j))then
          im=im+1
          if(im>=kretco.and.ikrt<knkeys)then
             ikrt=ikrt+1
             do jk=1,27
                keys(jk,ikrt)=blocks(jk,j,iu)
             enddo
          endif
       endif
    enddo
    knkeys=ikrt
    if(ikrt<=0)then
       kretco=-1
       return
    endif
    kretco=max(im-kretco-ikrt,0)
  end subroutine asimhk


  subroutine asimhc(kdev, kretco)

    implicit none

    !     asimof - close file

    !     arguments:

    !     kdev     i    unit number of file to close
    !     kretco   o    error return code

    !     0: ok.
    !     1: unit not found

    integer,intent(in):: kdev
    integer::kretco

    integer j, iu, inadm, icount, irec
    !     find unit in list of open files.
    inadm= jplrec/jpnkey-1
    iu = -1
    do j = 1, jpnfx
       if(munits(j)==kdev) iu = j
    enddo

    kretco = 1
    if (iu < 0) return

    !     write out modified adminstration blocks

    if(mopen(iu)==0)   goto 290
    if(mopen(iu)==-890)goto 290
    icount=0
    irec=1
210 continue
    if(admin(15,icount+1,iu)==1)then
       admin(15,icount+1,iu)=0
       call asimhwa(kdev,iu,irec,icount,inadm, kretco)
       if(kretco/=0)return
    endif
    icount=icount+1
    irec=admin(30,icount,iu)
    if(irec/=0)goto 210
290 continue

    !     unit found - release it.

    kretco = 0
    close(kdev)
    munits(iu) = -1
    mopen(iu)=-891
    nfld(iu)   =  0

    return
  end subroutine asimhc


  subroutine asimha(kdev, kmode, yadm, kadm, kladm, kerr)

    implicit none
    integer kdev, kmode, kladm, kerr
    integer kadm(kladm)
    character*(*) yadm

    integer ilyadm            ! length of yadm

    integer jpladm
    parameter(jpladm=3)
    integer iadm(jpladm)      ! temporary administration array
    integer iu, j             ! file locator and counter




    kerr=0

    ilyadm = min(8,len(yadm))

    if(kladm<=0)then
       kerr=-1
       return
    endif

    !     find unit in list of open files
    kerr=1
    j=1
    do while(j<jpnfx .and. kerr==1)
       if( munits(j)==kdev )then
          kerr = 0
          iu=j
       endif
       j=j+1
    enddo

    if(kerr==1)then
       return
    endif

    if(kmode==1)then
       if(mopen(iu)==0)then
          kerr=2
          return
       endif
       if(ilyadm>=4)call char2int(yadm(1:4),admin(13,1,iu),4)
       if(ilyadm>=8)call char2int(yadm(5:8),admin(14,1,iu),4)
       if(kladm>= 2)then
          if(kadm(2)>19980101)then
             admin(16,1,iu)=kadm( 2)
             if(kladm>= 3)admin(17,1,iu)=kadm( 3)
          endif
       endif
       return
    elseif(kmode==0)then
       kadm(1) = admin(12,1,iu)
       if(jpladm>= 2)iadm( 2)=admin(16,1,iu)
       if(jpladm>= 3)iadm( 3)=admin(17,1,iu)
       if(admin(12,1,iu)<258)then
          if(admin(13,1,iu)==0.and.admin(14,1,iu)==0)then
             yadm='n-a'
          else
             if(ilyadm>=4)call int2char(yadm(1:4),admin(13,1,iu),4)
             if(ilyadm>=8)call int2char(yadm(5:8),admin(14,1,iu),4)
          endif
          if(iadm(2)==0)then
             iadm(2)=19980101
             if(iadm(3)==0)iadm( 3)=120000
          endif
       else
          if(ilyadm>=4)call int2char(yadm(1:4),admin(13,1,iu),4)
          if(ilyadm>=8)call int2char(yadm(5:8),admin(14,1,iu),4)
       endif
       do j=2, min(jpladm, kladm)
          kadm(j)=iadm(j)
       enddo
       return
    else
       kerr=-1
       return
    endif
  end subroutine asimha
end module asimof
