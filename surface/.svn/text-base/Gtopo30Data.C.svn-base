#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"Gtopo30Data.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include<string.h>
#include<iostream>
#include"myRoutines.h"
#include <errno.h>

using namespace std;

Gtopo30Data::Gtopo30Data(Domain *domain,int N_, int M_, int ma,MPI_Comm comm):
  Data(domain,2),N(N_),M(M_),map_(ma){
}


Gtopo30Data::~Gtopo30Data(){
#ifdef DEBUG
  cout <<"~Gtopo30Data() called"<<endl;
#endif
  delete domain;
  delete[] datA;
}


void Gtopo30Data::readFile(char* fileName){

  filename = fileName;
  readData();
  //delete[] filename;
}

void Gtopo30Data::readData(){
  int Klon = domain->klon;
  int Klat = domain->klat;
  int idatastart = domain->idatastart;
  int jdatastart = domain->jdatastart;
//   int ijump = domain->iLam;
//   int jjump = domain->jLam;
  //  int N = domain->klon_global;

  int fp;
  short *g30data_short = new short[Klon*Klat];
  
  if ((fp = open (filename, O_RDONLY)) == -1 ) {
    fprintf (stderr, "g30read:  Could not open file %s %s\n",filename,strerror(errno));
    cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(2);
  }
  
  //jump to the start of the data
  long int pos = ((idatastart-1)+(jdatastart-1)*N)*sizeof(short);
  long int n = lseek(fp,pos,SEEK_SET);

  long int nr=0;
  for(int k=0;k<Klat;k++){
    long int nread = read(fp,g30data_short+k*Klon,Klon*sizeof(short));
    pos += N*sizeof(short);
    n = lseek(fp,pos,SEEK_SET);
    if(nread>0){
      nr += nread;
    }
    else{
      cout << "Error reading file "<<filename<<endl;
      cout <<"nread="<<nread<<" k="<<k<<"  Klat="<<Klat<<endl;
      fprintf (stderr, "gtoporead   >%s< sterror=%s\n",filename,strerror(errno));
      cout <<__FILE__<<"   "<<__LINE__<<endl;
      terminate_program(2);
    }
  }
  close(fp);

#ifdef LITTLE_ENDIAN
  union tal convert;
  char tmp[sizeof(short)];
  for(int i=0;i<Klon*Klat;i++){
    convert.tal = g30data_short[i];
    for(int d=0;d<sizeof(short);d++){
      tmp[d]=convert.bytes[sizeof(short)-d-1];
    }
    for(int d=0;d<sizeof(short);d++){
      convert.bytes[d]=tmp[d];
    }
    g30data_short[i] = convert.tal;
  }
#endif
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      //datA[i+j*Klon] = (double)(g30data_short[i+(Klat-1-j)*Klon]);
      data(i,j) = (double)(g30data_short[i-1+(j-1)*Klon]);
    }
  delete[] g30data_short;

#ifdef DEBUGIO
  int myid=0;
#ifdef MPI_SRC
  MPI_Comm_rank(localComm,&myid);
#endif

  char filena[128];
  strcpy(filena,domain->Name());
  //domain->print();

  char tmp2[sizeof(int)*8+1];  
  sprintf( tmp2, "%d", myid );
  char procnumber[5];
  strcpy(procnumber,".");
  if(myid<10)
    strcat(procnumber,"00");
  else if(myid<100)
    strcat(procnumber,"0");
  strcat(procnumber,tmp2);// e.g. ".001"
  strcat(filena,procnumber);
  strcat(filena,".bin");
  writeFileSolo(filena,1);
#endif
}


