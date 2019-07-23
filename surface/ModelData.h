#ifndef MODELDATA_
#define MODELDATA_
#include"Data.h"
#include"Domain.h"
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif  
class Gtopo30DataSet;

class EcoclimapData;
class LakeData;
class ModelData: public Data
{

 protected:
  double *nPoints; //number of points used in the upscaling
  bool useVar;

 public:
  ~ModelData();

  ModelData(Domain *d, int nc_in=1);
  ModelData(Domain *d, int nc_in, bool use_variance);
  inline double& var(int i, int j){       return data(2,i,j);}
  inline double& var(int c, int i, int j){return data(nc+c,i,j);} //datA[(i-1)+(j-1)*nc*2 + 2*nc*domain->klon*(c-1)];}
  void printDomain(int stride=1);
  void check();
  //  void loadData(Gtopo30DataSet *g30data); 
  void loadData(int c,LakeData *globaldata, double uninitialized);
  void loadData(int i,Data *ecoclimap,double uninitialized); 
  void loadData_z0(int i,EcoclimapData *ecoclimap,double UNINIECO); 
  void loadData(int c,EcoclimapData *globaldata,double uninitialized,EcoclimapData *F_wat);
  void interpolate(int c,ModelData *m1,int* compo, double* w,double UNINIECO);
  void loadMinData(int c,ModelData *an,int cstart,int cend,double UNINIECO);
  void loadMaxData(int c,ModelData *an,int cstart,int cend,double UNINIECO);
  double echoMax(int c);
  bool readFile(int c,char* file,MPI_Comm localComm);
  void writeFile(int c,char* file,MPI_Comm localComm);

};


#endif
