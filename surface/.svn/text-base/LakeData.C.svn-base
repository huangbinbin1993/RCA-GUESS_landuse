#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"LakeData.h"
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

LakeData::~LakeData(){
#ifdef DEBUG
  cout <<"~LakeData() is called"<<endl;
#endif
  delete[] datA;
  delete domain;
}

LakeData::LakeData(int sizN, int sizM,Domain *d):
  Data(d),N(sizN),M(sizM),ijump(d->idatastart),jjump(d->jdatastart){  
  initialized = false;
#ifdef DEBUG
  cout << " size of LakeData domain is"<<N<<" x "<<M<<endl;
#endif
}


void LakeData::readLakeFile(char* fileName){
#ifdef DEBUG
  cout <<"LakeData::readLakeFile"<<endl;
#endif
  filename = fileName;
  readData();
}


void LakeData::readData(){
#ifdef DEBUG
  cout <<"LakeData::readData("<<endl;
#endif
  //domain->jpos = false;
  int fp;
  int Klon = domain->klon;
  int Klat = domain->klat;
  int16_t *lakedata = new int16_t[Klon*Klat];
#ifdef DEBUG
  printf(" reading file = %s\n",filename);
#endif
  if ((fp = open(filename, O_RDONLY)) == -1 ) {
    fprintf (stderr, "ecoread:  Could not open file >%s< sterror=%s\n",filename,strerror(errno));
    cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(2);
  }

  //jump to the start of the data
  long int pos = ((ijump-1)+(jjump-1)*N)*sizeof(int16_t);
  long int n=lseek(fp,pos,SEEK_SET);
  long int nr=0;
  for(int k=0;k<Klat;k++){
    long int nread = read(fp,lakedata+k*Klon,Klon*sizeof(int16_t));
    pos += N*sizeof(int16_t);
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
      data(i,j) = (double)(lakedata[i-1+(j-1)*Klon]/10.0);

  delete[] lakedata;
  
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


