#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif

#include"EcoclimapDataSetReader.h"
#include"froutines.h"
#include<iostream>
#include<math.h>
#include"Domain.h"
#include"myRoutines.h"
#include"localParallel.h"

using namespace std;


extern "C" {
  
  void readgtopo30data_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* halo, 
			float* southF, float* westF, float* dlonF, float* dlatF,
			float* polonF, float* polatF,
			float *orographyF, float *stddevF,
			int* istart,int* jstart, float* gravit, MPI_Comm *localCommIn);
  void readecoclimapdata_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* halo, 
			  float* southF, float* westF, float* dlonF, float* dlatF,
			  float* polonF, float* polatF,
			  float* ecoclimap,int* monTh, int* dAy,
			  int* istart,int* jstart, MPI_Comm *localCommIn);
}

//tests ecoclimapreader and gtopo30reader

int main( int argc, char** argv ){
   int myid,nproc;
   MPI_Comm localComm;
   initialize_parallel(argc,argv,myid,nproc,localComm);

   int klon_global,  klat_global,    klev_global;
   double south,west,dlon,dlat,polon,polat;

   int month=3, day=1;
   int idatastart=1,jdatastart=1;
   int halo=1;
   read_domain_(&south, &west, &dlon, &dlat, &polon, &polat,&klon_global,&klat_global,&klev_global);
   int klon = klon_global, klat=klat_global;
   decompose_domain(klon_global,klat_global,klon,klat,idatastart,jdatastart,halo,myid,nproc);

   float *oro = new float[klon*klat];
   float *stdDev = new float[klon*klat];
   float *ecoclimap = new float[klon*klat*(43+6)];
   float gravit = 9.81;
   float polonF=(float)polon, polatF=(float)polat;
   float southF =(float)south, westF=(float)west;
   float dlonF=(float)dlon, dlatF=(float)dlat;
   int klonF=klon,klatF=klat,klonGlobal=klon_global,klatGlobal=klat_global;

   readecoclimapdata_( &klonF,&klatF,&klonGlobal, &klatGlobal, &halo, 
		       &southF, &westF, &dlonF, &dlatF,
		       &polonF, &polatF,
		       ecoclimap,&month,&day,
		       &idatastart,&jdatastart,&localComm);


   delete[] oro;
   delete[] stdDev;
   delete[] ecoclimap;
   close_parallel(localComm);
   return 0;
}




















