#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Fint int
#endif
#include"Gtopo30DataSetReader.h"
#include"myRoutines.h"
#include"froutines.h"
#include<string.h>
#include<iostream>
#include<math.h>
using namespace std;
#define NOTDEFINED -666

extern "C" {
#ifndef SINGLE_PRECISION
  void readgtopo30data_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* halo, 
			double* southF, double* westF, double* dlonF, double* dlatF,
			double* polonF, double* polatF,
			double *orographyF, double *stddevF,
			int* istart,int* jstart, double* gravit, MPI_Fint *localCommF)

#else
  void readgtopo30data_(int* klonF,int* klatF,int* klonGlobal, int* klatGlobal, int* halo, 
			float* southF, float* westF, float* dlonF, float* dlatF,
			float* polonF, float* polatF,
			float *orographyF, float *stddevF,
			int* istart,int* jstart, float* gravit, MPI_Fint *localCommF)
#endif
{
    int klon = *klonF;
    int klat = *klatF;
    int klon_global=*klonGlobal;
    int klat_global=*klatGlobal;
    int idatastart= *istart;
    int jdatastart= *jstart;
    int khalo=*halo;
    int klev_global=0;
#ifdef MPI_SRC
    MPI_Comm localComm=MPI_Comm_f2c(*localCommF);
#else
    int localComm = *localCommF;//in this case an integer
#endif
    char fpath[128];
    getgtopo30path_(fpath);
    double south,west,dlon,dlat,polon,polat;
    read_domain_(&south, &west, &dlon, &dlat, &polon, &polat,
		 &klon_global,&klat_global,&klev_global);

    if(dlon<=0 || dlat<=0){
      cout<<"dlon="<<dlon<<"  dlat="<<dlat<<endl;
      cout <<"must be positive"<<endl;
      terminate_program(9);
    }
    //These coordinates are on the rotated grid, which must be on the equator:
    double north = south + (klat_global-1)*dlat;
    if(!(north>0 && south <0)){
      cout << "south="<<south<<" north="<<north<<endl;
      cout <<"Rotated grid must be over the equator!"<<endl;
      cout <<__FILE__<<","<<__LINE__<<endl;
      terminate_program(9);
    }

    //output to RCA
#ifndef SINGLE_PRECISION
    *southF = south;
    *westF = west;
    *dlonF = dlon;
    *dlatF = dlat;
    *polonF = polon;
    *polatF = polat;
#else
    *southF = (float)south;
    *westF = (float)west;
    *dlonF = (float)dlon;
    *dlatF = (float)dlat;
    *polonF = (float)polon;
    *polatF = (float)polat;
#endif
    int ilen=0;
    while(fpath[ilen]!=' ')ilen++;
    fpath[ilen]='\0';
    char* path = new char[ilen+1];
    strcpy(path,fpath);

    bool global = false;
    
    Domain *modelDomain = new Domain("model", klon_global,   klat_global,  
				     klon,   klat,    
				     west,   south,   dlon,   dlat, localComm,
				     global, idatastart,   jdatastart, 
				     polon,   polat,  halo);

    Gtopo30DataSetReader *g30 = new Gtopo30DataSetReader(path,modelDomain,localComm);    

    for(int j=1;j<=klat;j++){
      for(int i=1;i<=klon;i++){
#ifndef SINGLE_PRECISION
 	orographyF[i-1+(j-1)*klon] = (*gravit)*g30->model->data(i,j); //since RCA  uses real*4
 	stddevF[i-1+(j-1)*klon] = sqrt(g30->model->var(i,j));
#else
 	orographyF[i-1+(j-1)*klon] = (*gravit)*(float)g30->model->data(i,j); //since RCA  uses real*4
 	stddevF[i-1+(j-1)*klon] = (float)(sqrt(g30->model->var(i,j)));
#endif
      }
    }
    delete[] path;
    delete g30;//should be called
    delete modelDomain;
  }
}
