#ifndef G30TOPODATASET_
#define G30TOPODATASET_
#include"Gtopo30Data.h"
#include"ModelData.h"
#include"Domain.h"
#include"Data.h"
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif

class Gtopo30DataSet{
 private:
  void makeRectangularMaps();
  int *imap;
  int *jmap;

  Domain **domains;
  //Gtopo30Data **g30topoData;
  char *path;
  int numberOfSets;//TOTAL number of sets 
  int nDomains;//the actual number of sets USED
  char **filenames;
  void connectFileName(int mapNumb);
  void connectFileNames();
  void defineDomains(ModelData *model,MPI_Comm localComm);
  bool *mapUsed;
  void setImapJmap();
  void readHeaderFile(char* actFile,double& g30_dlon,double& g30_dlat,
		      int& ncols,int& nrows,double& g30_lon,double& g30_lat);
  Domain* glueTilesIntoOne(MPI_Comm localComm);
  Domain* defineTiledDomain(char* name,ModelData *model,Domain *oneDomain,MPI_Comm localComm);


  //Domain *LAMDomain;
  void loadTiledData(MPI_Comm localComm);
  Domain* defineBoundingDomain(char* name,ModelData *model,MPI_Comm comm);
  int nlon();
  int nlat();
  double dlonI();
  double dlatI();

 public:
  char* getFileName(int map){return filenames[map];}
  int getNumberOfSets(){return numberOfSets;}
  int getMapNumber(double lon,double lat);
  int nSets(){return nDomains;}

  Data *gtopo;

  Gtopo30DataSet(char* path,ModelData *model,MPI_Comm Comm);

  ~Gtopo30DataSet();

  bool useMap(int i){return mapUsed[i];}
};
#endif
