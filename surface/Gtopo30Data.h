#ifndef G30TOPODATA_
#define G30TOPODATA_
#include"Data.h"
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif

class Gtopo30Data : public Data
{
 private:
#ifdef LITTLE_ENDIAN
  union tal{
    short tal;
    char bytes[sizeof(short)];
  };
#endif
  
  int N,M;
  char* filename;
  void readData();
  int map_;
 public:
  int map(){return map_;}
  Gtopo30Data(Domain *dy, int N, int M,int map,MPI_Comm localComm);
  ~Gtopo30Data();
  
  char* fileName(){return filename;}

  void readFile(char* fileName);
};
#endif
