#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"Gtopo30DataSetReader.h"
#include"froutines.h"
#include<string.h>
#include<iostream>
#include<math.h>


using namespace std;

Gtopo30DataSetReader::Gtopo30DataSetReader(char* path, Domain *modelDomain,MPI_Comm localComm){
  bool use_variance=true;
  model = new ModelData(modelDomain,1,use_variance);

  initialized = false;
  if(!model->readFile(1,"orography",localComm)){
    G30dataSet = new Gtopo30DataSet(path,model,localComm);
    initialized = true;
    
    model->loadData(1,G30dataSet->gtopo,-9999.0);//-9999 is the definition of a seapoint
    for(int j=1;j<=model->klat();j++)
      for(int i=1;i<=model->klon();i++)
	if(fabs(model->data(i,j)+9999.0)<0.01){
	  model->data(i,j) = 0.0;
	  model->var(i,j) = 0.0;
	}

    model->writeFile(1,"orography",localComm);
  }

}


Gtopo30DataSetReader::~Gtopo30DataSetReader(){
#ifdef DEBUG
  cout << "delete Gtopo30DataSetReader called"<<endl;
#endif
  if(initialized)
    delete G30dataSet;
}


