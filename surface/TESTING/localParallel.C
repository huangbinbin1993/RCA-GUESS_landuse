#ifdef MPI_SRC
#include"mpi.h"
#endif
#include<stdlib.h>
#include<iostream>
#include"localParallel.h"
using namespace std;

void initialize_parallel( int argc, char** argv, int& myid, int& nproc, MPI_Comm& localComm)
{
  myid = 0;
  nproc = 1;
  
#ifdef MPI_SRC
  MPI_Init( &argc, &argv );
  MPI_Comm_dup(MPI_COMM_WORLD,&localComm);
  MPI_Comm_rank( localComm, &myid );
  MPI_Comm_size( localComm, &nproc);

  if( myid == 0 ){
      int ver, subver;
      MPI_Get_version( &ver, &subver );
      cout << "Initializing MPI Version " << ver << "." << subver << endl;
    }
#else
  localComm = 0;
#endif    

#ifdef DEBUG_MPI
  char hnam[80];
  gethostname(hnam,80);
  cout << " Proc no " << myid << " running on " << hnam << endl;
  
  struct rlimit rlp;
  getrlimit( RLIMIT_CPU, &rlp );
  cout << "Cpu limits: current " << rlp.rlim_cur << " maximum: " << rlp.rlim_max << endl;
  
  if( rlp.rlim_cur < rlp.rlim_max ){
    rlp.rlim_cur = rlp.rlim_max;
    setrlimit( RLIMIT_CPU, &rlp );
  }
  cout << "Changed Cpu limits: current " << rlp.rlim_cur << " maximum: " << rlp.rlim_max << endl;
#endif
}
void close_parallel(MPI_Comm localComm )
{
#ifdef MPI_SRC
  int myid=0;
  MPI_Comm_rank(localComm, &myid );
  cout << myid << "  ";
  MPI_Barrier(localComm);
  for(int k=0;k<10000;k++)
    {}
  if(myid==0)
    cout <<" has reached the end of the program" << endl;
  
  MPI_Finalize( );
#endif    
}

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
  localfactorize_(&nproc,&klon_global,&klat_global,&one,&npx,&npy,&npz);
  
  int myI = myid%npx;
  int myJ = myid/npx;
  
  decomp1d(klon_global,npx,halo,myI,klon,idatastart);
  decomp1d(klat_global,npy,halo,myJ,klat,jdatastart);

}
