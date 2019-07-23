#ifndef LAKEDATA_
#define LAKEDATA_
#include"Data.h"
#ifdef MPI_SRC
#include"mpi.h"
#endif 
class LakeData:public Data
{
  int resolution;
  bool initialized;
  int N;
  int M;
  int ijump,jjump;
  char* filename; //current fileName
  void readData();

 public:
  LakeData(int sizeN,int sizeM, Domain *d);
  void readLakeFile(char* fileName);
  ~LakeData();
};
#endif
