#include"Gtopo30DataSetReader.h"
#include"froutines.h"
#include"myRoutines.h"
#include"localParallel.h"
#include<iostream>
#include<math.h>
using namespace std;
#define NOTDEFINED -666

using namespace std;
extern "C" {
  
  void readgtopo30data_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* halo, 
			float* southF, float* westF, float* dlonF, float* dlatF,
			float* polonF, float* polatF,
			float *orographyF, float *stddevF,
			int* istart,int* jstart, float* gravit, MPI_Comm *localCommIn);
}


//Tests the gtopo30reader
int main( int argc, char** argv ){
   int myid,nproc;
   MPI_Comm localComm;
   initialize_parallel(argc,argv,myid,nproc,localComm);

   double south,west,dlon,	dlat,	polon,	polat;
   int klon_global,  klat_global,    klev_global;
   int idatastart,jdatastart;
   int halo=1;
   read_domain_(&south, &west, &dlon, &dlat, &polon, &polat,&klon_global,&klat_global,&klev_global);

   int klon = klon_global, klat=klat_global;
   decompose_domain(klon_global,klat_global,klon,klat,idatastart,jdatastart,halo,myid,nproc);
   float *oro = new float[klon*klat];
   float *stdDev = new float[klon*klat];
   float gravit = 9.81;
   float polonF=(float)polon, polatF=(float)polat;
   float southF =(float)south, westF=(float)west;
   float dlonF=(float)dlon, dlatF=(float)dlat;
   int klonF=klon,klatF=klat,klonGlobal=klon_global,klatGlobal=klat_global;


   

   readgtopo30data_( &klonF,&klatF,&klonGlobal, &klatGlobal, &halo, 
		     &southF, &westF, &dlonF, &dlatF,
		     &polonF, &polatF,
		     oro, stdDev,
		     &idatastart,&jdatastart, &gravit,&localComm);

   delete[] oro;
   delete[] stdDev;
   close_parallel(localComm);
   return 0;
}




















