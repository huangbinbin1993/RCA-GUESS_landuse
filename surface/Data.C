#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"Data.h"
#include<iostream>
#include<stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include<string.h>
#include<math.h>
#include"myRoutines.h"
#define UNINI -666.0
using namespace std;

void Data::getBoundingBox(double lonMx[2],double latMy[2]){
#ifdef DEBUG
  cout <<" Data::getBoundingBox"<<endl;
#endif
  double minLon = 181.0;
  double minLat = 91.0;
  double maxLon = -181.0;
  double maxLat = -91.0;
  double lonx[4],laty[4];
  int I[2] = {1,klon()};
  int J[2] = {1,klat()};
  for(int j=1;j<=klat();j++){
    for(int i=0;i<2;i++){
      getCell(I[i],j,lonx,laty);
      for(int k=0;k<4;k++){
	minLon = MIN(minLon,lonx[k]);
	maxLon = MAX(maxLon,lonx[k]);
	minLat = MIN(minLat,laty[k]);
	maxLat = MAX(maxLat,laty[k]);
      }
    }
  }
  for(int j=0;j<2;j++){
    for(int i=1;i<=klon();i++){
      getCell(i,J[j],lonx,laty);
      for(int k=0;k<4;k++){
	minLon = MIN(minLon,lonx[k]);
	maxLon = MAX(maxLon,lonx[k]);
	minLat = MIN(minLat,laty[k]);
	maxLat = MAX(maxLat,laty[k]);
      }
    }
  }

//   //use index from 0 to N+1 in order to get some safety margin
//   for(int j=0;j<=klat()+1;j++){
//     for(int i=0;i<=klon()+1;i++){   
//       getCell(i,j,lonx,laty);
//       for(int k=0;k<4;k++){
// 	if(lonx[k]<minLon) minLon = lonx[k];
// 	if(lonx[k]>maxLon) maxLon = lonx[k];
// 	if(laty[k]<minLat) minLat = laty[k];
// 	if(laty[k]>maxLat) maxLat = laty[k];
//       }
//     }
//   }
  if(minLon<-180.0)minLon = -180.0;
  if(maxLon>180.0)maxLon = 180.0;
  if(minLat<-90.0)minLat = -90.0;
  if(maxLat>90.0)maxLat = 90.0;

  lonMx[0]=minLon; lonMx[1]=maxLon;
  latMy[0]=minLat; latMy[1]=maxLat;
}


void Data::getBoundingBoxRot(double *lonMx,double *latMy){
#ifdef DEBUG
  cout <<" Data::getBoundingBoxRot"<<endl;
#endif
  double minLon = 181.0;
  double minLat = 91.0;
  double maxLon = -181.0;
  double maxLat = -91.0;
  double lonx[4],laty[4];
  int I[2] = {0,klon()+1};
  int J[2] = {0,klat()+1};
  for(int j=0;j<=klat()+1;j++){
    for(int i=0;i<2;i++){
      getCellR(I[i],j,lonx,laty);
      for(int k=0;k<4;k++){
	minLon = MIN(minLon,lonx[k]);
	maxLon = MAX(maxLon,lonx[k]);
	minLat = MIN(minLat,laty[k]);
	maxLat = MAX(maxLat,laty[k]);
      }
    }
  }
  for(int j=0;j<2;j++){
    for(int i=0;i<=klon()+1;i++){
      getCell(i,J[j],lonx,laty);
      for(int k=0;k<4;k++){
	minLon = MIN(minLon,lonx[k]);
	maxLon = MAX(maxLon,lonx[k]);
	minLat = MIN(minLat,laty[k]);
	maxLat = MAX(maxLat,laty[k]);
      }
    }
  }

  if(minLon<-180.0){terminate_program(2);}//minLon = -180.0;
  if(maxLon>180.0){terminate_program(2);}//maxLon = 180.0;
  if(minLat<-90.0){terminate_program(2);}//minLat = -90.0;
  if(maxLat>90.0){terminate_program(2);}//maxLat = 90.0;

  lonMx[0]=minLon; lonMx[1]=maxLon;
  latMy[0]=minLat; latMy[1]=maxLat;
}



Data::Data(Domain *d,int nc_in): 
  nc(nc_in),  domain(d){
  int Klon=domain->klon;
  int Klat=domain->klat;
  if(Klon<=0 || Klat<=0 || nc<=0)
    terminate_program(9);
  datA = new double[Klon*Klat*nc];
  for(int j=0;j<Klat*Klon*nc;j++)
    datA[j]=UNINI;
  allocated=true;
}

void Data::printDomain(int stride){
  cout <<"%klon,klat="<<domain->klon<<" "<<domain->klat<<endl;
  cout <<"lon=[";
//   int Klon = domain->klon;
//   int Klat = domain->klat;
//   for(int j=1;j<=Klat;j+=stride){
//     for(int i=1;i<=Klon;i+=stride)
//       cout<<domain->lon(i,j)<<" ";
//     cout << endl;
//   }
//   //cout << "];"<<endl;
//   //cout << "lat=[";
//   for(int j=1;j<=Klat;j+=stride){
//     for(int i=1;i<=Klon;i+=stride)
//       //cout<<domain->lat(i,j)<<" ";
//     //cout << endl;
//   }
//   //cout << "];"<<endl;

//   if(nc==1){
//     //cout << "data=[";
//     for(int j=1;j<=Klat;j+=stride){
//       for(int i=1;i<=Klon;i+=stride)
// 	//cout<<datA[(i-1)+(j-1)*Klon]<<" ";
//       //cout << endl;
//     }
//     //cout << "];"<<endl;
//   }
}






int Data::writeHeader(char* filename,int myid,int p, int N, int M){
  int nr=0;
  if(myid==p){//p is default=0
    int fd=-1;
    cout <<"writing "<<filename<<endl;
    fd  = open( filename,  O_CREAT | O_TRUNC | O_WRONLY, 0660);
    if(N==-1 || M==-1){
      nr += write( fd, &(domain->klon_global), sizeof(int) );
      nr += write( fd, &(domain->klat_global), sizeof(int) );
    }
    else{
      nr += write( fd, &N, sizeof(int) );
      nr += write( fd, &M, sizeof(int) );
    }

    double *mmap = new double[6];
    mmap[0]=domain->south;
    mmap[1]=domain->west;
    mmap[2]=domain->dlon;
    mmap[3]=domain->dlat;
    mmap[4]=domain->polon;
    mmap[5]=domain->polat;
    nr += write( fd, mmap, 6*sizeof(double) );
    close(fd);
  }
  
  return 6*sizeof(double)+2*sizeof(int);
}


void Data::writeFile(char* f,MPI_Comm localComm,int c){

  int ni=klon();
  int nj=klat();
  double *da = new double[ni*nj];
  for(int j=1;j<=nj;j++)
    for(int i=1;i<=ni;i++){
      da[(i-1)+(j-1)*ni] = data(c,i,j);
    }
#ifdef DEBUG
  //cout <<"writing "<<f<<endl;
#endif
  int myid=0,nproc=1;
#ifdef MPI_SRC
  MPI_Comm_rank(localComm,&myid);
  MPI_Comm_size(localComm,&nproc);
  int dum=0,tag=999;
#endif  
  
  int headerStride = writeHeader(f,myid);

  if(myid==0){
    writeData(f,da,ni,nj,headerStride,localComm,c);
  }
#ifdef MPI_SRC
  else
    {
      MPI_Status status;
      MPI_Recv( &dum, 1, MPI_INT, myid-1, tag, localComm, &status );
      writeData(f,da,ni,nj,headerStride,localComm,c);
     }
  if( myid < nproc-1 )
    MPI_Send( &dum, 1, MPI_INT, myid+1, tag, localComm );
#endif
  if(myid==nproc-1)cout <<"Done writing "<<f<<endl;


  if(myid==nproc-1){
    int fd;
    fd  = open( f, O_RDONLY);
    long int pos = lseek(fd,0,SEEK_END);
    pos -= headerStride;
    pos = pos/sizeof(double);
    if(pos != domain->klon_global*domain->klat_global && 
       pos != 3*domain->klon_global*domain->klat_global){
      cout <<pos<<" != "<<domain->klon_global*domain->klat_global<<
	__FILE__<<" "<<__LINE__;
      cout <<"   while checking if it has correct size "<<f<<endl;
      cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(99);
    }
    close(fd);
  }
}

    
  

void Data::writeData(char* f,double *da,  int ni,int nj,int headerStride,
		     MPI_Comm localComm,
		     int c,
		     bool respectHalo){
  long int nr=0;
  int imap=0;
  int jmap=0;
  int myid=0;
  int N = domain->klon_global;
#ifdef MPI_SRC
    MPI_Comm_rank(localComm,&myid);
#endif
      
  if(domain->useHalo && respectHalo){
    imap = domain->idatastart-1;
    jmap = domain->jdatastart-1;
  }
  
  int fd;
  fd  = open( f, O_RDWR );
  //jump to the start of data
  int pos = headerStride + (imap + jmap*N)*sizeof(double);
  lseek(fd,pos, SEEK_SET);

  for(int j=0;j<nj;j++){
    nr += write(fd,da+j*ni,ni*sizeof(double));
    pos += N*sizeof(double);
    lseek(fd,pos, SEEK_SET);  
  }
  
  if(nr/sizeof(double)!=ni*nj){
    cout <<"Something wrong in writing "<< f <<" nr="<<nr/sizeof(double)<<" ni="<<ni<<" nj="<<nj<<" c="<<c<<endl;
    cout <<"imap="<<imap<<" jmap="<<jmap<<endl;
    cout <<__FILE__<<" "<<__LINE__<<endl;
    terminate_program(2);
  }
  close(fd);

  delete[] da;
}



int Data::readHeader(char* filename,int myid,MPI_Comm localComm){
  double *mmap = new double[6];
  int Klon_file=0,Klat_file=0;
  int nr=0;
  int fd=-1;
  if(myid==0)
    fd  = open( filename, O_RDONLY);
    
#ifdef MPI_SRC
  MPI_Bcast( &fd, 1, MPI_INT, 0, localComm );
#endif 
   
  if(fd==-1){
//     if(myid==0)
//       cout<<"file "<<filename <<" not found!"<<endl; 
    return -1;
  }


  if(myid==0){
    nr += read( fd, &Klon_file, sizeof(int) );
    nr += read( fd, &Klat_file, sizeof(int) );
  }

#ifdef MPI_SRC
  int dims[2] = {Klon_file,Klat_file};
  MPI_Bcast( dims, 2, MPI_INT, 0, localComm );
  Klon_file = dims[0];
  Klat_file = dims[1];
#endif 
    
  if(Klon_file != domain->klon_global || Klat_file != domain->klat_global){
    if(myid==0){
      cout <<"Klon or Klat do not match!!! file "<<filename<<endl;
      cout <<Klon_file<<" !="<< domain->klon_global<<
	"   "<< Klat_file<<" != "<<domain->klat_global<<endl;
      cout <<"Will have to recreate this file ("<<filename<<")"  <<endl;
    }
    return -1;
  }
  if(myid==0){
    nr += read( fd, mmap, 6*sizeof(double) );

    //test to see if the file is full
    long int pos = lseek(fd,0,SEEK_END);
    pos -= nr;
    pos = pos/sizeof(double);
    if(pos != domain->klon_global*domain->klat_global && pos != 3*domain->klon_global*domain->klat_global){
      cout <<pos<<" != "<<domain->klon_global*domain->klat_global<<
	__FILE__<<" "<<__LINE__;
      cout <<"   while probing "<<filename<<endl;
      nr = -1;
    }
    close(fd);
  }
  
#ifdef MPI_SRC
  double tol=0.00000000001;
  MPI_Bcast( mmap, 6, MPI_DOUBLE, 0, localComm );
  if(fabs(mmap[0]-domain->south)>tol || fabs(mmap[1]-domain->west)>tol || 
     fabs(mmap[2]-domain->dlon)>tol || fabs(mmap[3]-domain->dlat)>tol || 
     fabs(mmap[4]-domain->polon)>tol || fabs(mmap[5]-domain->polat)>tol){
    if(myid==0){
      cout <<"Geometric data does not match! file "<<filename<<endl;
      cout <<mmap[0]<<" !="<< domain->south<<endl;
      cout <<mmap[1]<<" !="<< domain->west<<endl;
      cout <<mmap[2]<<" !="<< domain->dlon<<endl;
      cout <<mmap[3]<<" !="<< domain->dlat<<endl;
      cout <<mmap[4]<<" !="<< domain->polon<<endl;
      cout<< mmap[5]<<" !="<< domain->polat<<endl;

      cout <<"Will have to recreate this file ("<<filename<<")"<<endl;
    }
    return -1;
  }

  domain->south = mmap[0];
  domain->west  = mmap[1];
  domain->dlon  = mmap[2];
  domain->dlat  = mmap[3];
  domain->polon = mmap[4];
  domain->polat = mmap[5];
  MPI_Bcast(&nr,1,MPI_INT,0,localComm);
#endif  
  delete[] mmap;
  return nr;
}



bool Data::readFile(char* f,MPI_Comm localComm,int c){
  
  double *da = new double[domain->klon*domain->klat];

  int myid=0,nproc=1;
  
#ifdef MPI_SRC
  MPI_Comm_rank(localComm, &myid);
  MPI_Comm_size(localComm,&nproc);
#endif

  int headerStride;
  if((headerStride = readHeader(f,myid,localComm))==-1){
    //cout <<"Header info does not match"<<endl;
    return false;
  }
  if(myid==0)cout<<"reading "<<f<<endl;
  int fd  = open( f, O_RDONLY);
  int pos= headerStride + (domain->idatastart-1 + (domain->jdatastart-1)*domain->klon_global)*sizeof(double);
  lseek(fd,pos,SEEK_SET);
  
  //Make a parallell read here...
  long int nr=0;
  for(int j=0;j<domain->klat;j++){
    nr += read( fd, da+j*domain->klon, domain->klon*sizeof(double) );
    pos += domain->klon_global*sizeof(double);
    lseek(fd,pos,SEEK_SET);
  }
  
  if(nr/sizeof(double) != domain->klon*domain->klat){
    cout << nr/sizeof(double)<< "!="<< domain->klon*domain->klat<<__FILE__<<__LINE__<<endl;
    cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(2);
  }
  close(fd);
  
  for(int j=1;j<=domain->klat;j++)
    for(int i=1;i<=domain->klon;i++)
      data(c,i,j) = da[(i-1)+(j-1)*domain->klon];


  delete[] da;
  return true;
}


void Data::findMinMax(double& minlon, double& maxlon, double& minlat, double& maxlat){
  minlat=666,maxlat=-666;
  minlon=666,maxlon=-666;

  //for this it is sufficient ? to sweep the boundaries of the domain
  for(int j=1;j<=klat();j++){
    for(int i=1;i<=klon();i++){
      if(lon(i,j)<minlon) 	minlon=lon(i,j);
      if(lon(i,j)>maxlon) 	maxlon=lon(i,j);
      if(lat(i,j)<minlat)       minlat=lat(i,j);
      if(lat(i,j)>maxlat)       maxlat=lat(i,j);
    }
  }
}


void Data::writeFileSolo(char* f,int c){

  double *da = new double[domain->klon*domain->klat];
  double *ln = new double[domain->klon*domain->klat];
  double *lt = new double[domain->klon*domain->klat];

  for(int j=1;j<=domain->klat;j++)
    for(int i=1;i<=domain->klon;i++){
      da[(i-1)+(j-1)*domain->klon] = data(c,i,j);
      ln[(i-1)+(j-1)*domain->klon] = lon(i,j);
      lt[(i-1)+(j-1)*domain->klon] = lat(i,j);
    }

  cout <<"writing "<<f<<endl;

  int fd=-1;
  long int nr=0;
  fd  = open( f,  O_CREAT | O_TRUNC | O_WRONLY, 0660);
  
  nr += write( fd, &(domain->klon), sizeof(int) );
  nr += write( fd, &(domain->klat), sizeof(int) );
  

  double *mmap = new double[6];
  mmap[0]=domain->south;
  mmap[1]=domain->west;
  mmap[2]=domain->dlon;
  mmap[3]=domain->dlat;
  mmap[4]=domain->polon;
  mmap[5]=domain->polat;
  nr += write( fd, mmap, 6*sizeof(double) );

nr += write(fd,da,domain->klon*domain->klat*sizeof(double));
nr += write(fd,ln,domain->klon*domain->klat*sizeof(double));
nr += write(fd,lt,domain->klon*domain->klat*sizeof(double));

 close(fd);
//CONTINUE TO WRITE A WHOLE NEW WRITE!!!
//   int myid=0;  
// #ifdef MPI_SRC
//   MPI_Comm_rank(localComm,&myid);
// #endif
//   int headerStride = writeHeader(f,myid,myid,domain->klon,domain->klat);

//   writeData(f,da,domain->klon,domain->klat,headerStride,c,false);
//   delete[] da;
}
