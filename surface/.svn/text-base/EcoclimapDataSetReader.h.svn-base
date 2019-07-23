#ifndef ECOCLIMAPDATASETREADER_
#define ECOCLIMAPDATASETREADER_
#include"ModelData.h"
#include"EcoclimapDataSet.h"
#include"Domain.h"
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif 
class ModelData;
class EcoclimapDataSet;
class EcoclimapDataSetReader{
 private:

  int nlon();     //size of the whole file
  int nlat();     //size of the whole file
  double dlonI();
  double dlatI();

  bool initialized;
  int numberOfTimeDepVars;
  Domain* defineBoundingDomain(char* name,MPI_Comm comm);
  int *var;
  int nc;
  int iOrder;
  bool *readFromFile;
  bool *timeDep;
 public:
  ModelData *model;
  ModelData *annualContainer;
  EcoclimapDataSet *ecodataSet;
  EcoclimapDataSetReader(char* path,int month,Domain *domain,MPI_Comm localComm);  
  ~EcoclimapDataSetReader();
  void getEcoclimapData(int month,MPI_Comm localComm);
  double* getData2RCA(){return model->getv();}
  int ncomp(){return nc;}
};
#endif
