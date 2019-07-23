#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif

void initialize_parallel( int argc, char** argv, int& myid, int& nproc,  MPI_Comm& localComm);
void close_parallel(MPI_Comm localComm );
void decompose_domain(int klon_global,int klat_global,int &klon,int &klat,int &idatastart,int &jdatastart,int halo,int myid,int nproc);
extern "C" { 
  void localfactorize_(int* npr,int* nig,int* njg,int* nkg,int* np1,int* np2,int* np3);
}
