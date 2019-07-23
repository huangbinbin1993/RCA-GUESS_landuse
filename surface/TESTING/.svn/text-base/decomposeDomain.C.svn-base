#include"froutines.h"

using namespace std;

void decomp1d(int N, int nproc, int halo,int myid, int &n, int &is){
  n = (N-2*halo)/nproc;
  int rem = (N-2*halo)%nproc;

  int *n_ar = new int[nproc];
  for(int i=0;i<nproc;i++){
    n_ar[i]=n;
    if(i<rem)
      n_ar[i]++;
    n_ar[i]+= 2*halo;
  }
  n = n_ar[myid];

  int *is_ar = new int[nproc];
  for(int i=0;i<nproc;i++){
    if(i<rem)
      is_ar[i] = 1 + i*(n_ar[i]-2*halo);
    else
      is_ar[i] = 1 + i*(n_ar[i]-2*halo)+rem;
  }
  is = is_ar[myid];
  delete[] is_ar;
  delete[] n_ar;
}


void decompose_domain(int klon_global,int klat_global,int &klon,int &klat,int &idatastart,int &jdatastart,int halo,int myid,int nproc){
  int one=1;
  int npx=1,npy=1,npz=1;
  factorize_(&nproc,&klon_global,&klat_global,&one,&npx,&npy,&npz);
  
  int myI = myid%npx;
  int myJ = myid/npx;
  
  decomp1d(klon_global,npx,halo,myI,klon,idatastart);
  decomp1d(klat_global,npy,halo,myJ,klat,jdatastart);

}
