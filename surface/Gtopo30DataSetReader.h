#ifndef GTOPO30DATASETREADER_
#define GTOPO30DATASETREADER_
#include"ModelData.h"
#include"Gtopo30DataSet.h"

#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif

class ModelData;
class Gtopo30DataSet;
class Gtopo30DataSetReader{
 private:
  void checkWhichMapUsed();
  bool initialized;
  
 public:
  ModelData *model;
  Gtopo30DataSet *G30dataSet;
  Gtopo30DataSetReader(char* path, Domain *modelDomain,MPI_Comm localComm);
  ~Gtopo30DataSetReader();
  double* getData2RCA(){return model->getv();}
};
#endif
