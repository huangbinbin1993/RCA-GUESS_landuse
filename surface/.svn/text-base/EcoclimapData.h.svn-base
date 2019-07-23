#ifndef ECOCLIMAPDATA_
#define ECOCLIMAPDATA_
#include"Data.h"
#ifdef MPI_SRC
#include"mpi.h"
#endif 
class EcoclimapData:public Data
{
  bool initialized;
  int N;
  int M;
  int ijump,jjump;
  char* filename; //current fileName
  void readData();

 public:
  EcoclimapData(int sizeN,int sizeM, Domain *d);
  void readEcoFile(char* fileName);
  ~EcoclimapData();
};
#endif
