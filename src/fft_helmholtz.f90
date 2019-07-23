module fft_helmholtz
  use decomp,only:stop_program,realkind
  implicit none
  private

  public fft44
contains
  subroutine fft44a(a,b,c,d,kinc1,kinc2,kinc3,kinc4,klot,kn,trigs)

    !     subroutine fft44a - first stage of postprocessing for
    !     multiple sine transform
    implicit none
    real(kind=realkind) a(*),b(*),c(*),d(*),trigs(*)
    integer kinc1,kinc2,kinc3,kinc4,klot,kn
    integer inc1,inc2,inc3,inc4,lot,n
    integer ia,ja,j,nh,ib,jb,k,jbbase,iabase,ibbase,jabase

    real(kind=realkind) si,co

    inc1=kinc1
    inc2=kinc2
    inc3=kinc3
    inc4=kinc4
    lot=klot
    n=kn
    ia=1-inc3
    ja=1-inc4
    do  j=1,lot
       ia=ia+inc3
       ja=ja+inc4
       c(ja)=a(ia)-b(ia)
       d(ja)=a(ia)+b(ia)
    enddo
    nh=n/2

    if(nh-(nh/2)*2==0)then  !      if(nh-(nh/2)*2/=0) go to 30
       ib=(nh*inc1)/2+1-inc3
       jb=(nh*inc2)/2+1-inc4
       do j=1,lot
          ib=ib+inc3
          jb=jb+inc4
          c(jb)=2.0_realkind*b(ib)
          d(jb)=2.0_realkind*a(ib)
       enddo
    endif

    iabase=inc1+1
    ibbase=(nh-1)*inc1+1
    jabase=inc2+1
    jbbase=(nh-1)*inc2+1
    do  k=3,nh,2
       ia=iabase-inc3
       ib=ibbase-inc3
       ja=jabase-inc4
       jb=jbbase-inc4
       co=trigs(n+k)
       si=trigs(n+k+1)
       do  j=1,lot
          ia=ia+inc3
          ib=ib+inc3
          ja=ja+inc4
          jb=jb+inc4
          c(ja)=(si*(b(ia)+b(ib))-co*(a(ia)-a(ib)))+(b(ia)-b(ib))
          c(jb)=(si*(b(ia)+b(ib))-co*(a(ia)-a(ib)))-(b(ia)-b(ib))
          d(ja)=(a(ia)+a(ib))+(si*(a(ia)-a(ib))+co*(b(ia)+b(ib)))
          d(jb)=(a(ia)+a(ib))-(si*(a(ia)-a(ib))+co*(b(ia)+b(ib)))
       enddo
       iabase=iabase+inc1
       ibbase=ibbase-inc1
       jabase=jabase+inc2
       jbbase=jbbase-inc2
    enddo

    return
  end subroutine fft44a
  subroutine fft44b(a,b,kinc1,kinc2,kinc3,kinc4,klot,kn,ksep)

    !     subroutine fft44b - second stage of postprocessing for
    !     multiple sine transform   (isep=0 if a,b are the same array)
    implicit none
    real(kind=realkind) a(*),b(*)
    integer kinc1,kinc2,kinc3,kinc4,ksep,kn,klot
    integer inc1,inc2,inc3,inc4,lot,n
    integer isep,ink1,ink2,i,l
    integer ia,ja,j,nh,ib,jb,k,ibase,jbase


    inc1=kinc1
    inc2=kinc2
    inc3=kinc3
    inc4=kinc4
    lot=klot
    n=kn
    isep=ksep

    ia=1-inc3
    ib=(n-2)*inc1+1-inc3
    ja=1-inc4
    jb=(n-2)*inc2+1-inc4
    do j=1,lot
       ia=ia+inc3
       ib=ib+inc3
       ja=ja+inc4
       jb=jb+inc4
       b(jb+inc2)=-a(ia)
       b(ja+inc2)=a(ia+inc1)
       b(ja)=0.0_realkind
    enddo

    ibase=3*inc1+1
    jbase=3*inc2+1
    ink1=inc1+inc1
    ink2=inc2+inc2
    nh=n/2
    do k=3,nh
       i=ibase-inc3
       j=jbase-inc4
       do l=1,lot
          i=i+inc3
          j=j+inc4
          b(j)=b(j-ink2)+a(i)
       enddo
       ibase=ibase+ink1
       jbase=jbase+ink2
    enddo

    if(isep/=0)then
       ibase=2*inc1+1
       jbase=2*inc2+1
       do k=2,nh
          i=ibase-inc3
          j=jbase-inc4
          do l=1,lot
             i=i+inc3
             j=j+inc4
             b(j)=a(i)
          enddo
          ibase=ibase+ink1
          jbase=jbase+ink2
       enddo
    endif

    return
  end subroutine fft44b
  subroutine fft44(a,work,kinc,kjump,kn,klot,ksign,nfax,trigs,klon_global)

    !         subroutine "fft44" - multiple fast sine transform
    !         procedure used to convert to convert to half-length complex
    !         transforms is the inverse of that given by cooley, lewis ' welch
    !         (j. sound vib., vol. 12 (1970), 315-337)
    !     
    !         a is the array containing input ' output data
    !         work is an area of size (n+1)lot
    !         trigs is a previously prepared list of trig function values
    !         nfax is a previously prepared list of factors of n/2
    !         inc is the increment within each data "vector"
    !             (e.g. inc=1 for consecutively stored data)
    !         jump is the increment between the start of each data vector
    !         n is the length of the data vectors (including a(1)=0.0)
    !         lot is the number of data vectors
    !         isign = +1 for transform from spectral to gridpoint
    !               = -1 for transform from gridpoint to spectral
    !     
    !         vectorization is achieved on cray by doing the transforms in
    !             parallel
    !     
    !          n.b. n is assumed to be an even number
    !     
    !         definition of transforms:
    !         -------------------------
    !     
    !         isign=+1: x(j)=sum(k=0,1,...,n-1)(b(k)*sin(j*k*pi/n)),
    !                                                j=0,1,...,n-1
    !     
    !         isign=-1: b(k)=(2/n)*sum(j=0,1,...,n-1)(x(j)*sin(j*k*pi/n)))
    !                                                k=0,1,...,n-1
    !     


    implicit none
    integer kinc1,kinc2,kinc3,kinc4,klot,kn
    integer inc1,inc2,inc3,inc4,lot,n
    integer igo,k,la,infax,j,ib,ja,jb,i,iabase,jbbase
    integer nh,jump,jabase,ibbase,ia,nx,nh3,ksign,ink
    integer isign,kjump,kinc,inc
    integer   klon_global
    real(kind=realkind) a(*), work(*),sink,scale
    integer   nfax(10)
    real(kind=realkind)      trigs(klon_global/2*5)

    inc=kinc
    jump=kjump
    n=kn
    lot=klot
    isign=ksign
    nh=n/2
    nh3=3*n/2
    nx=n+1
    ink=inc+inc

    !     preprocessing

    scale=0.25_realkind
    if (isign==-1) scale=0.5_realkind/real(n,realkind)

    iabase=inc+1
    ibbase=(n-1)*inc+1
    jabase=2
    jbbase=n

    do  i=2,nh
       ia=iabase-jump
       ib=ibbase-jump
       ja=jabase-nx
       jb=jbbase-nx
       sink=trigs(nh3+i)
       do  j=1,lot
          ia=ia+jump
          ib=ib+jump
          ja=ja+nx
          jb=jb+nx
          work(ja)=scale*(sink*(a(ia)+a(ib))+(a(ia) -a(ib)))
          work(jb)=scale*(sink*(a(ia)+a(ib))-(a(ia) -a(ib)))
       enddo
       iabase=iabase+inc
       ibbase=ibbase-inc
       jabase=jabase+1
       jbbase=jbbase-1
    enddo
    ia=iabase-jump
    ja=1-nx
    jb=nh+1-nx
    do j=1,lot
       ia=ia+jump
       ja=ja+nx
       jb=jb+nx
       work(ja)=0.0_realkind
       work(jb)=4.0_realkind*scale*a(ia)
    enddo

    !     complex fft
    infax=nfax(1)
    la=1
    igo=40
    do k=1,infax
       if (igo==50)then 
          call vpassm(a(1),a(1+inc),work(1),work(2),ink,2, &
               jump,nx,lot, nh,nfax(k+1),trigs(1),la)
          igo=40
       else
          call vpassm(work(1),work(2),a(1),a(1+inc),2, &
               ink,nx,jump,lot,nh,nfax(k+1),trigs(1),la)
          igo=50
       endif
       la=la*nfax(k+1)
    enddo

    !     postprocessing - first stage
    if (infax-(infax/2)*2==1)then
       call fft44a(a(1),a(1+inc),work(1),work(2),ink,2,jump,nx,  lot,n,trigs(1))
    else
       call fft44a(work(1),work(2),a(1),a(1+inc),2,ink,nx,jump, lot,n,trigs(1))
    endif

    !     postprocessing - second stage
    if(infax-(infax/2)*2==1)then
       call fft44b(work,a(1),1,inc,nx,jump,lot,n,1)
    else
       call fft44b(a(1),a(1),inc,inc,jump,jump,lot,n,0)
    endif
    return
  end subroutine fft44

  subroutine vpassm(a,b,c,d,kinc1,kinc2,kinc3,kinc4,klot,kn,kfac,trigs,kla)

    !     subroutine 'vpassm' - multiple version of 'vpassa'
    !     performs one pass through data
    !     as part of multiple complex fft routine
    !     a is first real input vector
    !     b is first imaginary input vector
    !     c is first real output vector
    !     d is first imaginary output vector
    !     trigs is precalculated table of sines ' cosines
    !     kinc1 is addressing increment for a and b
    !     kinc2 is addressing increment for c and d
    !     kinc3 is addressing increment between a's + b's
    !     kinc4 is addressing increment between c's + d's
    !     klot is the number of vectors
    !     kn is length of vectors
    !     kfac is current factor of kn
    !     array of trigonometric functions
    !     kla is product of previous factors

    implicit none
    real(kind=realkind) sin36,cos36,sin72,cos72,sin60
    real(kind=realkind) a(*),b(*),c(*),d(*),trigs(*)
    integer inc1,kinc1,inc2,kinc2,inc3,kinc3,inc4,kinc4
    integer lot,klot,n,kn,ifac,kfac,la,kla,m,iink,jink,jump
    integer ibase,jbase,igo,ia,ja,ib,jb,l,i,j,ijk
    integer la1,k,kb,ic,jc,kc,id,jd,kd,ie,je,ke
    real(kind=realkind) c1,s1,c2,s2,c3,s3,c4,s4
    data sin36/0.587785252292473_realkind/,cos36/0.809016994374947_realkind/, &
         sin72/0.951056516295154_realkind/,cos72/0.309016994374947_realkind/, &
         sin60/0.866025403784437_realkind/

    inc1=kinc1
    inc2=kinc2
    inc3=kinc3
    inc4=kinc4
    lot=klot
    n=kn
    ifac=kfac
    la=kla
    m=n/ifac
    iink=m*inc1
    jink=la*inc2
    jump=(ifac-1)*jink
    ibase=0
    jbase=0
    igo=ifac-1
    if (igo>4) call stop_program( 'vpassm - factor in fast fourier transform > 5')
    if(igo==1)then

       !     coding for factor 2

       ia=1
       ja=1
       ib=ia+iink
       jb=ja+jink
       do 20 l=1,la
          i=ibase-inc3
          j=jbase-inc4

          do 15 ijk=1,lot
             i=i+inc3
             j=j+inc4
             c(ja+j)=a(ia+i)+a(ib+i)
             d(ja+j)=b(ia+i)+b(ib+i)
             c(jb+j)=a(ia+i)-a(ib+i)
             d(jb+j)=b(ia+i)-b(ib+i)
15        enddo
          ibase=ibase+inc1
          jbase=jbase+inc2
20     enddo
       if (la==m) return
       la1=la+1
       jbase=jbase+jump
       do 40 k=la1,m,la
          kb=k+k-2
          c1=trigs(kb+1)
          s1=trigs(kb+2)
          do 30 l=1,la
             i=ibase-inc3
             j=jbase-inc4
             do 25 ijk=1,lot
                i=i+inc3
                j=j+inc4
                c(ja+j)=a(ia+i)+a(ib+i)
                d(ja+j)=b(ia+i)+b(ib+i)
                c(jb+j)=c1*(a(ia+i)-a(ib+i))- s1*(b(ia+i)-b(ib+i))
                d(jb+j)=s1*(a(ia+i)-a(ib+i))+  c1*(b(ia+i)-b(ib+i))
25           enddo
             ibase=ibase+inc1
             jbase=jbase+inc2
30        enddo
          jbase=jbase+jump
40     enddo
       return
    endif
    if(igo==2)then

       !     coding for factor 3

       ia=1
       ja=1
       ib=ia+iink
       jb=ja+jink
       ic=ib+iink
       jc=jb+jink
       do 60 l=1,la
          i=ibase-inc3
          j=jbase-inc4
          do 55 ijk=1,lot
             i=i+inc3
             j=j+inc4
             c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
             d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
             c(jb+j)=(a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))- &
                  (sin60*(b(ib+i)-b(ic+i)))
             c(jc+j)=(a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))+ &
                  (sin60*(b(ib+i)-b(ic+i)))
             d(jb+j)=(b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))+ &
                  (sin60*(a(ib+i)-a(ic+i)))
             d(jc+j)=(b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))- &
                  (sin60*(a(ib+i)-a(ic+i)))
55        enddo
          ibase=ibase+inc1
          jbase=jbase+inc2
60     enddo
       if (la==m) return
       la1=la+1
       jbase=jbase+jump
       do 80 k=la1,m,la
          kb=k+k-2
          kc=2*kb
          c1=trigs(kb+1)
          s1=trigs(kb+2)
          c2=trigs(kc+1)
          s2=trigs(kc+2)
          do 70 l=1,la
             i=ibase-inc3
             j=jbase-inc4
             do 65 ijk=1,lot
                i=i+inc3
                j=j+inc4
                c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                c(jb+j)=                                   &
                     c1*((a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))-   &
                     (sin60*(b(ib+i)-b(ic+i))))             &
                     -s1*((b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))+  &
                     (sin60*(a(ib+i)-a(ic+i))))              
                d(jb+j)=                                    &
                     s1*((a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))-   &
                     (sin60*(b(ib+i)-b(ic+i))))             &
                     +c1*((b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))+  &
                     (sin60*(a(ib+i)-a(ic+i))))              
                c(jc+j)=                                    &
                     c2*((a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))+   &
                     (sin60*(b(ib+i)-b(ic+i))))             &
                     -s2*((b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))-  &
                     (sin60*(a(ib+i)-a(ic+i))))              
                d(jc+j)=                                    &
                     s2*((a(ia+i)-0.5_realkind*(a(ib+i)+a(ic+i)))+   &
                     (sin60*(b(ib+i)-b(ic+i))))             &
                     +c2*((b(ia+i)-0.5_realkind*(b(ib+i)+b(ic+i)))-  &
                     (sin60*(a(ib+i)-a(ic+i))))
65           enddo
             ibase=ibase+inc1
             jbase=jbase+inc2
70        enddo
          jbase=jbase+jump
80     enddo
       return
    endif
    if(igo==3)then
       !     coding for factor 4

       ia=1
       ja=1
       ib=ia+iink
       jb=ja+jink
       ic=ib+iink
       jc=jb+jink
       id=ic+iink
       jd=jc+jink
       do 100 l=1,la
          i=ibase-inc3
          j=jbase-inc4
          do 95 ijk=1,lot
             i=i+inc3
             j=j+inc4
             c(ja+j)=(a(ia+i)+a(ic+i))+  &
                  (a(ib+i)+a(id+i))        
             c(jc+j)=(a(ia+i)+a(ic+i))-   &
                  (a(ib+i)+a(id+i))        
             d(ja+j)=(b(ia+i)+b(ic+i))+   &
                  (b(ib+i)+b(id+i))        
             d(jc+j)=(b(ia+i)+b(ic+i))-   &
                  (b(ib+i)+b(id+i))        
             c(jb+j)=(a(ia+i)-a(ic+i))-   &
                  (b(ib+i)-b(id+i))        
             c(jd+j)=(a(ia+i)-a(ic+i))+   &
                  (b(ib+i)-b(id+i))        
             d(jb+j)=(b(ia+i)-b(ic+i))+   &
                  (a(ib+i)-a(id+i))        
             d(jd+j)=(b(ia+i)-b(ic+i))-   &
                  (a(ib+i)-a(id+i))
95        enddo
          ibase=ibase+inc1
          jbase=jbase+inc2
100    enddo
       if (la==m) return
       la1=la+1
       jbase=jbase+jump
       do 120 k=la1,m,la
          kb=k+k-2
          kc=2*kb
          kd=kc+kb
          c1=trigs(kb+1)
          s1=trigs(kb+2)
          c2=trigs(kc+1)
          s2=trigs(kc+2)
          c3=trigs(kd+1)
          s3=trigs(kd+2)
          do 110 l=1,la
             i=ibase-inc3
             j=jbase-inc4
             do 105 ijk=1,lot
                i=i+inc3
                j=j+inc4
                c(ja+j)=(a(ia+i)+a(ic+i))+        &
                     (a(ib+i)+a(id+i))              
                d(ja+j)=(b(ia+i)+b(ic+i))+         &
                     (b(ib+i)+b(id+i))              
                c(jc+j)=                           &
                     c2*((a(ia+i)+a(ic+i))-        &
                     (a(ib+i)+a(id+i)))            &
                     -s2*((b(ia+i)+b(ic+i))-       &
                     (b(ib+i)+b(id+i)))             
                d(jc+j)=                           &
                     s2*((a(ia+i)+a(ic+i))-        &
                     (a(ib+i)+a(id+i)))            &
                     +c2*((b(ia+i)+b(ic+i))-       &
                     (b(ib+i)+b(id+i)))             
                c(jb+j)=                           &
                     c1*((a(ia+i)-a(ic+i))-        &
                     (b(ib+i)-b(id+i)))            &
                     -s1*((b(ia+i)-b(ic+i))+       &
                     (a(ib+i)-a(id+i)))             
                d(jb+j)=                           &
                     s1*((a(ia+i)-a(ic+i))-        &
                     (b(ib+i)-b(id+i)))            &
                     +c1*((b(ia+i)-b(ic+i))+       &
                     (a(ib+i)-a(id+i)))             
                c(jd+j)=                           &
                     c3*((a(ia+i)-a(ic+i))+        &
                     (b(ib+i)-b(id+i)))            &
                     -s3*((b(ia+i)-b(ic+i))-       &
                     (a(ib+i)-a(id+i)))             
                d(jd+j)=                           &
                     s3*((a(ia+i)-a(ic+i))+        &
                     (b(ib+i)-b(id+i)))            &
                     +c3*((b(ia+i)-b(ic+i))-       &
                     (a(ib+i)-a(id+i)))
105          enddo
             ibase=ibase+inc1
             jbase=jbase+inc2
110       enddo
          jbase=jbase+jump
120    enddo
       return
    endif
    if(igo==4)then
       !     coding for factor 5

       ia=1
       ja=1
       ib=ia+iink
       jb=ja+jink
       ic=ib+iink
       jc=jb+jink
       id=ic+iink
       jd=jc+jink
       ie=id+iink
       je=jd+jink
       do 140 l=1,la
          i=ibase-inc3
          j=jbase-inc4
          do 135 ijk=1,lot
             i=i+inc3
             j=j+inc4
             c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+          &
                  (a(ic+i)+a(id+i))                        
             d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+           &
                  (b(ic+i)+b(id+i))                        
             c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-    &
                  cos36*(a(ic+i)+a(id+i)))                &
                  -(sin72*(b(ib+i)-b(ie+i))+              &
                  sin36*(b(ic+i)-b(id+i)))                 
             c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-    &
                  cos36*(a(ic+i)+a(id+i)))                &
                  +(sin72*(b(ib+i)-b(ie+i))+              &
                  sin36*(b(ic+i)-b(id+i)))                 
             d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                  cos36*(b(ic+i)+b(id+i)))                &
                  +(sin72*(a(ib+i)-a(ie+i))+              &
                  sin36*(a(ic+i)-a(id+i)))                 
             d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                  cos36*(b(ic+i)+b(id+i)))                &
                  -(sin72*(a(ib+i)-a(ie+i))+              &
                  sin36*(a(ic+i)-a(id+i)))                 
             c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+    &
                  cos72*(a(ic+i)+a(id+i)))                &
                  -(sin36*(b(ib+i)-b(ie+i))-              &
                  sin72*(b(ic+i)-b(id+i)))                 
             c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+    &
                  cos72*(a(ic+i)+a(id+i)))                &
                  +(sin36*(b(ib+i)-b(ie+i))-              &
                  sin72*(b(ic+i)-b(id+i)))                 
             d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                  cos72*(b(ic+i)+b(id+i)))                &
                  +(sin36*(a(ib+i)-a(ie+i))-              &
                  sin72*(a(ic+i)-a(id+i)))                 
             d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                  cos72*(b(ic+i)+b(id+i)))                &
                  -(sin36*(a(ib+i)-a(ie+i))-              &
                  sin72*(a(ic+i)-a(id+i)))
135       enddo
          ibase=ibase+inc1
          jbase=jbase+inc2
140    enddo
       if (la==m) return
       la1=la+1
       jbase=jbase+jump
       do 160 k=la1,m,la
          kb=k+k-2
          kc=2*kb
          kd=kc+kb
          ke=kd+kb
          c1=trigs(kb+1)
          s1=trigs(kb+2)
          c2=trigs(kc+1)
          s2=trigs(kc+2)
          c3=trigs(kd+1)
          s3=trigs(kd+2)
          c4=trigs(ke+1)
          s4=trigs(ke+2)
          do 150 l=1,la
             i=ibase-inc3
             j=jbase-inc4
             do 145 ijk=1,lot
                i=i+inc3
                j=j+inc4
                c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+            &
                     (a(ic+i)+a(id+i))                          
                d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+             &
                     (b(ic+i)+b(id+i))                          
                c(jb+j)=                                       &
                     c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-     &
                     cos36*(a(ic+i)+a(id+i)))                  &
                     -(sin72*(b(ib+i)-b(ie+i))+                &
                     sin36*(b(ic+i)-b(id+i))))                 &
                     -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                     cos36*(b(ic+i)+b(id+i)))                  &
                     +(sin72*(a(ib+i)-a(ie+i))+                &
                     sin36*(a(ic+i)-a(id+i))))                  
                d(jb+j)=                                       &
                     s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-     &
                     cos36*(a(ic+i)+a(id+i)))                  &
                     -(sin72*(b(ib+i)-b(ie+i))+                &
                     sin36*(b(ic+i)-b(id+i))))                 &
                     +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                     cos36*(b(ic+i)+b(id+i)))                  &
                     +(sin72*(a(ib+i)-a(ie+i))+                &
                     sin36*(a(ic+i)-a(id+i))))                  
                c(je+j)=                                       &
                     c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-     &
                     cos36*(a(ic+i)+a(id+i)))                  &
                     +(sin72*(b(ib+i)-b(ie+i))+                &
                     sin36*(b(ic+i)-b(id+i))))                 &
                     -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                     cos36*(b(ic+i)+b(id+i)))                  &
                     -(sin72*(a(ib+i)-a(ie+i))+                &
                     sin36*(a(ic+i)-a(id+i))))                  
                d(je+j)=                                       &
                     s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-     &
                     cos36*(a(ic+i)+a(id+i)))                  &
                     +(sin72*(b(ib+i)-b(ie+i))+                &
                     sin36*(b(ic+i)-b(id+i))))                 &
                     +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-    &
                     cos36*(b(ic+i)+b(id+i)))                  &
                     -(sin72*(a(ib+i)-a(ie+i))+                &
                     sin36*(a(ic+i)-a(id+i))))                  
                c(jc+j)=                                       &
                     c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+     &
                     cos72*(a(ic+i)+a(id+i)))                  &
                     -(sin36*(b(ib+i)-b(ie+i))-                &
                     sin72*(b(ic+i)-b(id+i))))                 &
                     -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                     cos72*(b(ic+i)+b(id+i)))                  &
                     +(sin36*(a(ib+i)-a(ie+i))-                &
                     sin72*(a(ic+i)-a(id+i))))                  
                d(jc+j)=                                       &
                     s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+     &
                     cos72*(a(ic+i)+a(id+i)))                  &
                     -(sin36*(b(ib+i)-b(ie+i))-                &
                     sin72*(b(ic+i)-b(id+i))))                 &
                     +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                     cos72*(b(ic+i)+b(id+i)))                  &
                     +(sin36*(a(ib+i)-a(ie+i))-                &
                     sin72*(a(ic+i)-a(id+i))))                  
                c(jd+j)=                                       &
                     c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+     &
                     cos72*(a(ic+i)+a(id+i)))                  &
                     +(sin36*(b(ib+i)-b(ie+i))-                &
                     sin72*(b(ic+i)-b(id+i))))                 &
                     -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                     cos72*(b(ic+i)+b(id+i)))                  &
                     -(sin36*(a(ib+i)-a(ie+i))-                &
                     sin72*(a(ic+i)-a(id+i))))                  
                d(jd+j)=                                       &
                     s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+     &
                     cos72*(a(ic+i)+a(id+i)))                  &
                     +(sin36*(b(ib+i)-b(ie+i))-                &
                     sin72*(b(ic+i)-b(id+i))))                 &
                     +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+    &
                     cos72*(b(ic+i)+b(id+i)))                  &
                     -(sin36*(a(ib+i)-a(ie+i))-                &
                     sin72*(a(ic+i)-a(id+i))))
145          enddo
             ibase=ibase+inc1
             jbase=jbase+inc2
150       enddo
          jbase=jbase+jump
160    enddo
       return
    endif
  end subroutine vpassm
end module fft_helmholtz
