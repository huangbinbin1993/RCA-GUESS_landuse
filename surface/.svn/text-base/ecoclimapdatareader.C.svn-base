
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Fint int
#endif

#include"EcoclimapDataSetReader.h"
#include"froutines.h"
#include<iostream>
#include<string.h>
#include<math.h>
#include"myRoutines.h"
using namespace std;

extern "C" {
#ifndef SINGLE_PRECISION
  void readecoclimapdata_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* haloIn, 
			  double* southF, double* westF, double* dlonF, double* dlatF,
			  double* polonF, double* polatF,
			  double* ecoclimap,int* monTh, 
			  int* istart,int* jstart, MPI_Fint *localCommF)
#else
   void readecoclimapdata_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* haloIn, 
 			  float* southF, float* westF, float* dlonF, float* dlatF,
 			  float* polonF, float* polatF,
 			  float* ecoclimap,int* monTh, 
 			  int* istart,int* jstart, MPI_Fint *localCommF)
#endif
    {
#ifdef MPI_SRC
    MPI_Comm localComm=MPI_Comm_f2c(*localCommF);
#else
    int localComm = *localCommF;//in this case an integer
#endif
    double south,west,dlon,dlat,polon,polat;
    int kloG,klaG,kleG;
    read_domain_(&south,&west,&dlon,&dlat,&polon,&polat,&kloG,&klaG,&kleG);
    
    //output to RCA
#ifndef SINGLE_PRECISION
    *southF = (float)south;
    *westF = (float)west;
    *dlonF = (float)dlon;
    *dlatF = (float)dlat;
    *polonF = (float)polon;
    *polatF = (float)polat;
#else
    *southF = (double)south;
    *westF = (double)west;
    *dlonF = (double)dlon;
    *dlatF = (double)dlat;
    *polonF = (double)polon;
    *polatF = (double)polat;
#endif

    int month = *monTh;
    int klon = *klonF;
    int klat = *klatF;
    int klon_global=*klonGlobal;
    int klat_global=*klatGlobal;
    int idatastart= *istart;
    int jdatastart= *jstart;
    int halo = *haloIn;
    bool global = false;
    Domain *modelDomain = new Domain("model", klon_global,   klat_global,  klon,   klat,    
				     west,   south,   dlon,   dlat, localComm, global, idatastart,   jdatastart, 
				     polon,   polat,  halo); //This is used to initiate ed, and also deleted from inside ed

     char fpath[128];
     getecopath_(fpath);
     int ilen=0;
     while(fpath[ilen]!=' ')ilen++;
     fpath[ilen]='\0';
     char* path = new char[ilen+1];
     strcpy(path,fpath);

     
     EcoclimapDataSetReader *ed = new EcoclimapDataSetReader(path,month,modelDomain,localComm);
#ifndef SINGLE_PRECISION
     double undefi = 999.0;
     double eps = 0.001;
#else
     float undefi = 999.0;
     float eps = 0.001;
#endif

     for(int c=1;c<=43+6 ;c++){
       if(c==10){//F_wat
	 for(int j=1;j<=klat;j++){
	   for(int i=1;i<=klon;i++){
#ifndef SINGLE_PRECISION
	     ecoclimap[(i-1)+(j-1)*klon+(c-1)*klon*klat] = 1.0-ed->model->data(c,i,j); //since RCA  uses real*8
#else
	     ecoclimap[(i-1)+(j-1)*klon+(c-1)*klon*klat] = 1.0-(float)ed->model->data(c,i,j); //since RCA  uses real*4
#endif
	   }
	 }
       }
       else{
	 for(int j=1;j<=klat;j++){
	   for(int i=1;i<=klon;i++){
#ifndef SINGLE_PRECISION
	     ecoclimap[(i-1)+(j-1)*klon+(c-1)*klon*klat] = ed->model->data(c,i,j); //since RCA  uses real*8
#else
	     ecoclimap[(i-1)+(j-1)*klon+(c-1)*klon*klat] = (float)ed->model->data(c,i,j); //since RCA  uses real*4
#endif
	   }
	 }
       }
     }
     delete ed;
     delete[] path;
  }
}
