#ifndef DATA_
#define DATA_
#include<iostream>
#include<stdlib.h>
#include"Domain.h"
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif  

using namespace std;
class Data
{
 private:
  bool allocated;

 protected:
  double *datA;
  int nc;

  int readHeader(char* filename,int myid,MPI_Comm localComm);//returns the header size in number of bits 
  int writeHeader(char* filename,int myid,int p=0,int N=-1,int M=-1);//returns the header size in number of bits 
  void writeData(char* f,double *da, int ni, int nj, int headerStride,
		 MPI_Comm localComm,
		 int c=1,
		 bool respectHalo=true);

 public:
  void getBoundingBox(double lonMx[2],double latMy[2]);
  void getBoundingBoxRot(double lonMx[2],double latMy[2]);
  Domain *domain;
  Data(){};
  Data(Domain *d, int nc=1);

  virtual ~Data(){}

  inline double& data(int i, int j){return datA[(i-1) + domain->klon*(j-1)];}

  inline double& data(int c,int i, int j){return datA[(i-1)+domain->klon*(j-1) + domain->klat*domain->klon*(c-1)];}
  void printDomain(int stride=1);
  inline double* getv(){return datA;}

  void writeFile(char* filename,MPI_Comm localComm,int c=1);
  void writeFileSolo(char* filename,int c=1);
  bool readFile(char* filename,MPI_Comm localComm,int c=1);
  double polon(){return domain->polon;}
  double polat(){return domain->polat;}
  double south(){return domain->south;}
  double west(){return domain->west;}
  double dlon(){return domain->dlon;};
  double dlat(){return domain->dlat;};
  double lon(double i, double j){return domain->lon(i,j);};
  double lat(double i, double j){return domain->lat(i,j);};
  double lonRot(double i, double j){return domain->lonRot_(i,j);};
  double latRot(double i, double j){return domain->latRot_(i,j);};
  double lonindex(double x, double y){return domain->lonindex(x,y);};
  double latindex(double x, double y){return domain->latindex(x,y);};
  double lonindexR(double x, double y){return domain->lonindexR(x,y);};
  double latindexR(double x, double y){return domain->latindexR(x,y);};
  int klon(){return domain->klon;};
  int klat(){return domain->klat;};
  int klon_global(){return domain->klon_global;};
  int klat_global(){return domain->klat_global;};
  void findMinMax(double& minlon, double& maxlon, double& minlat, double& maxlat);
  int idatastart(){return domain->idatastart;}
  int jdatastart(){return domain->jdatastart;}
  void getCell(double i,double j,double lonx[4], double laty[4]){
    domain->getCell(i,j,lonx,laty);};
  void getCellR(double i,double j,double lonx[4], double laty[4]){
    domain->getCellR(i,j,lonx,laty);};
  void correctCoordinates(double& lon, double& lat){
    domain->correctCoordinates(lon,lat);};

};
#endif
