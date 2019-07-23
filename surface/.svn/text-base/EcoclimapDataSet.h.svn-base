#ifndef ECOCLIMAPDATASET
#define ECOCLIMAPDATASET
#include"EcoclimapData.h"
#include"LakeData.h"
class EcoclimapDataSet{
 private:

  char *path;//ptr
  int numberOfVariables;
  char** months; 
  char** tiles;

  char **filenames;
  char **pureFilenames;
  bool *timeDep;
  bool *readFromfile;

  void connectFilenames();
  bool initialized;
  
 public:
  int getVariableNr(int nn,int month);
  int getnn(int variable);
  char* getFileName(int i){return filenames[i-1];}
  char* getPureFileName(int i){return pureFilenames[i-1];}
  EcoclimapDataSet(char* path, Domain *eco, int N, int M);
  ~EcoclimapDataSet();
  EcoclimapData *ecoclimapData;
  LakeData *lakeData;

  EcoclimapData* getData(){return ecoclimapData;}
  int getNumberOfTimeDepVars();
  bool isTimeDep(int i);
  bool readFromFile(int i);
  int getNumberOfVariables(){return numberOfVariables;}
  void setTimeDep(int i,bool dep);
};

#endif
