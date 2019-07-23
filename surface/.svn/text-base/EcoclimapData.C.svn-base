#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"EcoclimapData.h"
#include"myRoutines.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include<string.h>
#include<iostream>
#include <errno.h>

//extern int errno;

using namespace std;

//input should be the domain of interest i.e. over what domain to upscale
//need to know klon,klat,south,west,dlon,dlat before
EcoclimapData::~EcoclimapData(){
#ifdef DEBUG
  cout <<"~EcoclimapData() is called"<<endl;
#endif
  delete[] datA;
  delete domain;
}

EcoclimapData::EcoclimapData(int sizN, int sizM,Domain *d):
  Data(d),N(sizN),M(sizM),ijump(d->idatastart),jjump(d->jdatastart){  
  initialized = false;
#ifdef DEBUG
  cout << " size of EcoclimapData domain is"<<N<<" x "<<M<<endl;
#endif
}


void EcoclimapData::readEcoFile(char* fileName){
#ifdef DEBUG
  cout <<"EcoclimapData::readEcoFile "<<fileName<<endl;
#endif
  filename = fileName;
  readData();
}


void EcoclimapData::readData(){
#ifdef DEBUG
  cout <<"EcoclimapData::readData("<<endl;
#endif
  int fp;
  int Klon = domain->klon;
  int Klat = domain->klat;

  float *ecodata = new float[Klon*Klat];
#ifdef DEBUG
  printf(" reading file = %s\n",filename);
#endif
  if ((fp = open(filename, O_RDONLY)) == -1 ) {
    fprintf (stderr, "ecoread:  Could not open file >%s< sterror=%s\n",filename,strerror(errno));
    cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(2);
  }

  //jump to the start of the data
  long int pos = ((ijump-1)+(jjump-1)*N)*sizeof(float);
  long int n=lseek(fp,pos,SEEK_SET);
  long int nr=0;
  for(int k=0;k<Klat;k++){
    long int nread = read(fp,ecodata+k*Klon,Klon*sizeof(float));
    pos += N*sizeof(float);
    n = lseek(fp,pos,SEEK_SET);
    if(nread>0){
      nr += nread;
    }
    else{
      cout << "Error reading file "<<filename<<endl;
      cout <<"nread="<<nread<<" k="<<k<<endl;
      cout <<__FILE__<<"   "<<__LINE__<<endl;
      terminate_program(2);
    }
  }
  close(fp);

  //if need convert little_endian -> Big endian or vice versa

  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++)
      data(i,j) = (double)(ecodata[i-1+(j-1)*Klon]);//here the indexing is NOT reversed!

  delete[] ecodata;
  
#ifdef MPI_DEBUGIO
  int myid=0;
#ifdef MPI_SRC
  MPI_Comm_rank(localComm,&myid);
#endif

  char filena[128];
  strcpy(filena,"TEST");
  domain->print();

  char tmp[sizeof(int)*8+1];  
  sprintf( tmp, "%d", myid );
  char procnumber[5];
  strcpy(procnumber,".");
  if(myid<10)
    strcat(procnumber,"00");
  else if(myid<100)
    strcat(procnumber,"0");
  strcat(procnumber,tmp);// e.g. ".001"
  strcat(filena,procnumber);
  strcat(filena,".bin");
  writeFileSolo(filena,1);
#endif
}


