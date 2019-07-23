//#ifdef MPI_SRC
//#include"mpi.h"
//#endif
#include<string.h>
#include"Domain.h"
#include"froutines.h"
#include"myRoutines.h"
#include<iostream>
#include<math.h>
#define UNINI -666.0
using namespace std;

Domain::~Domain(){
  //  cout <<"deleting domain "<<name<<endl;
  delete[] name;
}

Domain:: Domain(Domain *d):
  klon_global(d->klon_global),
  klat_global(d->klat_global),
  klon(d->klon),
  klat(d->klat),
  west(d->west),
  south(d->south),
  north(d->north),
  east(d->east),
  dlon(d->dlon),
  dlat(d->dlat),
  idatastart(d->idatastart),
  jdatastart(d->jdatastart),
  polon(d->polon),
  polat(d->polat),
  global(d->global),
  ipos(d->ipos),
  jpos(d->jpos),
  halo(d->halo),
  useHalo(d->useHalo){
  name = new char[(int)strlen(d->name)+1];
  strcpy(name,d->name);
  //north = d->south + (d->klat_global-1)*d->dlat;
  //east = d->west + (d->klon_global-1)*d->dlon;
}

Domain:: Domain(char* n,int klon_global_,  int klat_global_, int klon_,  int klat_,    
		double west_,  double south_,  double dlon_,  double dlat_,MPI_Comm comm, bool global_,
		int idatastart_,int jdatastart_, double polon_,  double polat_,bool ipos_,
		bool jpos_,int halo_):
  klon_global(klon_global_),
  klat_global(klat_global_),
  klon(klon_),
  klat(klat_),
  west(west_),
  south(south_),
  dlon(dlon_),
  dlat(dlat_),
  idatastart(idatastart_),
  jdatastart(jdatastart_),
  polon(polon_),
  polat(polat_),
  global(global_),
  ipos(ipos_),
  jpos(jpos_),
  halo(halo_),
  useHalo(halo_>0){
  name = new char[(int)strlen(n)+1];
  strcpy(name,n);
  north = south + (klat_global-1)*dlon;
  east = west + (klon_global-1)*dlat;
}


// void Domain::print(char* msg,int myid){
//   cout <<msg<<endl;
//   cout <<"Domain " << name<<" :: ("<<myid<<")"<<endl;
//   cout<<"klon       ="<<klon       <<endl;
//   cout<<"klat       ="<<klat       <<endl;
//   if(useHalo)
//     cout<<"halo       ="<<halo       <<endl;
//   cout<<"idatastart ="<<idatastart <<endl;
//   cout<<"jdatastart ="<<jdatastart <<endl;
//   cout<<"klon_global="<<klon_global<<endl;
//   cout<<"klat_global="<<klat_global<<endl;
//   cout<<"west       ="<<west     <<  "  => east="<<west+(klon_global-1)*dlon<<endl;
//   cout<<"south      ="<<south    << "  => north="<<south+(klat_global-1)*dlat<<endl;
//   cout<<"lon(1,1) = "<<lon(1,1)<<"   lon(klon,klat) = "<<lon(klon,klat)<<endl;
//   cout<<"lat(1,1) = "<<lat(1,1)<<"   lat(klon,klat) = "<<lat(klon,klat)<<endl;
//   cout<<"lonR(1,1) = "<<lonRot_(1,1)<<"   lonR(klon,klat) = "<<lonRot_(klon,klat)<<endl;
//   cout<<"latR(1,1) = "<<latRot_(1,1)<<"   latR(klon,klat) = "<<latRot_(klon,klat)<<endl;
//   cout<<"dlon       ="<<dlon       <<endl;
//   cout<<"dlat       ="<<dlat       <<endl;
//   cout<<"polon      ="<<polon      <<endl;
//   cout<<"polat      ="<<polat      <<endl;
//   cout<<"global    ="<< (global  ? "T" : "F")    <<endl;
//   cout<<"ipos    ="<< (ipos  ? "T" : "F")    <<endl;
//   cout<<"jpos    ="<< (jpos  ? "T" : "F")    <<endl;

//   if(!global){
//     double lonxRot[]={west,
// 		      west+(klon_global-1)*dlon,
// 		      west+(klon_global-1)*dlon,
// 		      west};
//     double latyRot[]={south,
// 		      south,
// 		      south+(klat_global-1)*dlat,
// 		      south+(klat_global-1)*dlat};
//     double lonx[4],  laty[4];
//     int four=4;
//     int one=1;
    
//     rot2reg_(lonx,laty,lonxRot,latyRot,&four,&one,&polon,&polat);//fortran routine

//   }
  
// }

// void Domain::print(char* msg,MPI_Comm localComm){
//   int myid=0;
//   int nproc=1;
// #ifdef MPI_SRC
//   MPI_Comm_rank(localComm,&myid);
//   MPI_Comm_size(localComm,&nproc);
//   int dum=0,tag=999;
// #endif
//   if(myid==0)
//     print(msg,myid);
// #ifdef MPI_SRC
//   else
//     {
//       MPI_Status status;
//       MPI_Recv( &dum, 1, MPI_INT, myid-1, tag, localComm, &status );
//       print(msg,myid);
//      }
//    if( myid < nproc-1 )
//      MPI_Send( &dum, 1, MPI_INT, myid+1, tag, localComm );
//    MPI_Barrier(localComm );
// #endif
// }


void Domain::getCell(double i,double j,double lonx[4], double laty[4]){
  double x = lonRot_(i,j);
  double y = latRot_(i,j);
  double dx = 0.5*dlon;
  double dy = 0.5*dlat;
  double lonxRot[]={x+dx,x+dx,x-dx,x-dx};
  double latyRot[]={y+dy,y-dy,y-dy,y+dy};
  int four=4;
  int one=1;
  if(global)
    for(int i=0;i<4;i++){
      lonx[i] = lonxRot[i]; 
      laty[i] = latyRot[i]; 
      correctCoordinates(lonx[i],laty[i]);
    }
  else{
    rot2reg_(lonx,laty,lonxRot,latyRot,&four,&one,&polon,&polat);//fortran routine
  }
}  



void Domain::getCellR(double i,double j,double lonx[4], double laty[4]){
  double x = lonRot_(i,j);
  double y = latRot_(i,j);
  double dx = 0.5*dlon;
  double dy = 0.5*dlat;
  double lonxRot[]={x+dx,x+dx,x-dx,x-dx};
  double latyRot[]={y+dy,y-dy,y-dy,y+dy};
  for(int i=0;i<4;i++){
    lonx[i] = lonxRot[i]; 
    laty[i] = latyRot[i]; 
    correctCoordinates(lonx[i],laty[i]);
  }
}  






// void Domain::setUpCoordinates(){
//   cout <<"This routine is under development !!"<<__FILE__<<__LINE__<<endl;
//   terminate_program(9);
//   double *Rloncoords = new double[klon*klat];
//   double *Rlatcoords = new double[klon*klat];

//   for(int j=1;j<=klat;j++){
//     for(int i=1;i<=klon;i++){
//       Rlatcoords[(i-1)+(j-1)*klon] = south + (j-1+jdatastart-1)*dlat;//Global
//       Rloncoords[(i-1)+(j-1)*klon] = west  + (i-1+idatastart-1)*dlon;
//     }
//   }
//   if(!global)
//     rot2reg_(loncoords,latcoords,Rloncoords,Rlatcoords,&klon,&klat,&polon,&polat);//fortran routine
//   else{
//     for(int j=0;j<klat*klon;j++){
//       loncoords[j] = Rloncoords[j];
//       latcoords[j] = Rlatcoords[j];
//     }
//   }
//   for(int j=1;j<=klat;j++)
//     for(int i=1;i<=klon;i++)
//       correctCoordinates(loncoords[(i-1)+(j-1)*klon],latcoords[(i-1)+(j-1)*klon]);
  
//   delete[] Rloncoords;
//   delete[] Rlatcoords;
// }


double Domain::modulo(double a, double b){
  cout <<"SHOULD NOT BE CALLED"<<endl;
  terminate_program(2);
  int result = static_cast<int>( a / b );
  return a - static_cast<double>( result ) * b;
}
void Domain::correctCoordinates(double& lon, double& lat){
  if(fabs(lat)>90.0){
    if(lat>90.0)
      lat = 90.0 - (lat-90);
    else
      lat = -90.0 - (lat+90);
    lon += 180.0;
  }
  while(fabs(lon)>180){
    lon -= SIGN(lon)*360.0;
  }
}
  //   if(lon<-180.0)lon+=360.0;
  //   if(lon>180.0)lon-=360.0;
  //   if(lat<-90.0){
  //     lat=-90-(lat+90);
  //     lon+=180.0;
  //     if(lon<-180.0)lon+=360.0;
  //     if(lon>180.0)lon-=360.0;
  //   }
  //   if(lat>90.0){
  //     lat=90-(lat-90);
  //     lon += 180.0;
  //     if(lon<-180.0)lon+=360.0;
  //     if(lon>180.0)lon-=360.0;
  //   }






double Domain::polygonArea(double *X, double *Y, int points) {

  double  area=0. ;
  int     i, j=0  ;

  for (i=0; i<points; i++) {
    j++; if (j==points) j=0;
    area+=(X[i]+X[j])*(Y[i]-Y[j]); }

  return area*.5; 
} 

bool Domain::isInside(int ii, int jj,double vertx[4], double verty[4]){
  //Does point (ii,jj) lie inside cell {vertx,verty}?Z
  //vertx = [-180, 180]
  //verty = [-90, 90]
  double x = lon(ii,jj);
  double y = lat(ii,jj);
  return isInside(x,y, vertx,verty);
  
}


bool Domain::isInside(double x,double y,double vertx[4], double verty[4]){
  //Does point (ii,jj) lie inside cell {vertx,verty}?
  //vertx = [-180, 180]
  //verty = [-90, 90]
  bool inside=false;
  
  int i, j;
  for (i = 0, j = 3; i < 4; j = i++) {//j=3,0,1,2  i=0,1,2,3
    if ( ((verty[i]>y) != (verty[j]>y)) &&  //same ((true!=true) or (false!=false)) if the line connecting i to j is on the same 'side' of y
	 (x <= (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      inside =! inside;
  }

  return inside;
}


// bool Domain::isPointInside(int i,int j,double x, double y)
// {
//   //Answers: does (x,y) lie inside the grid cell 'surrounding' grind point with index (i,j)
//   //corners of the mapped square
//   double vertx[]={UNINI,UNINI,UNINI,UNINI};
//   double verty[]={UNINI,UNINI,UNINI,UNINI}; 

//   getCell(i,j,vertx,verty);
//   return isInside(x,y, vertx,verty);
// }


double Domain::lon(double i, double j){ //regular space
  double x=UNINI;
  double y=UNINI;
  if(global){
    x = lonRot_(i,j);
    y = latRot_(i,j);
    correctCoordinates(x,y);
  }
  else{
    double xR = lonRot_(i,j);
    double yR = latRot_(i,j);
    int one=1;
    rot2reg_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
  }
  return x;
}


double Domain:: latRotated(double i, double j){
  double x = lonRot_(i,j);
  double y = latRot_(i,j);
  int one=1;
  double xR,yR;
  reg2rot_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
  //correctCoordinates(xR,yR);
  return xR;
}

double Domain::lonRotated(double i, double j){
  double x = lonRot_(i,j);
  double y = latRot_(i,j);
  int one=1;
  double xR,yR;
  reg2rot_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
  //correctCoordinates(xR,yR);
  return yR;
}

void Domain:: rotated(double i, double j,double &xR, double &yR){
  double x = lonRot_(i,j);
  double y = latRot_(i,j);
  int one=1;
  reg2rot_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
  //correctCoordinates(xR,yR);
}





double Domain::lonRot_(double i, double j){
  if(ipos)
    return west + (i-1+idatastart-1)*dlon;
  else{
    //double east = west + (klon_global-1)*dlon;
    return east - (i-1+idatastart-1)*dlon;
  }
}





/////////////////////////////////////////



double Domain::lat(double i, double j){
  double x = UNINI;
  double y = UNINI;
  if(global){
    x = lonRot_(i,j);
    y = latRot_(i,j);
    correctCoordinates(x,y);
  }
  else{
    double xR = lonRot_(i,j);
    double yR = latRot_(j,j);
    int one=1;
    rot2reg_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
  }
  return y;
}


double Domain::latRot_(double i, double j){
  if(jpos)
    return south + (j-1+jdatastart-1)*dlat;
  else{
    return north - (j-1+jdatastart-1)*dlat;
  }
}

//////////////////////////////////////////




double Domain::lonindex(double x, double y){
  double index;
  if(global){
    index=lonRotindex_(x,y);
  }
  else{
    //project x and y to rotated and get index there
    int one=1;
    double xR,yR;
    reg2rot_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
    index=lonRotindex_(xR,yR);
  }
  return index;
}

double Domain::lonindexR(double x, double y){
  double xR,yR;
  int one=1;
  rot2reg_(&xR,&yR,&x,&y,&one,&one,&polon,&polat);//fortran routine
  //correctCoordinates(xR,yR);
  return lonRotindex_(xR,yR);
}

double Domain::lonRotindex_(double x, double y){
  if(ipos){
    return (x - west)/dlon + (1-idatastart+1);
  }
  else{
    return (west-x)/dlon + klon_global-idatastart+1;
  }
}




double Domain::latindex(double x, double y){
  double index;
  if(global){
    index=latRotindex_(x,y);
  }
  else{
    //project x and y to rotated and get index there
    int one=1;
    double xR,yR;
    reg2rot_(&x,&y,&xR,&yR,&one,&one,&polon,&polat);//fortran routine
    index=latRotindex_(xR,yR);
  }
  return index;
}


double Domain::latindexR(double x, double y){
  double xR,yR;
  int one=1;
  rot2reg_(&xR,&yR,&x,&y,&one,&one,&polon,&polat);//fortran routine
  //correctCoordinates(xR,yR);
 
  return latRotindex_(xR,yR);
}


double Domain::latRotindex_(double x, double y){
  if(jpos){
    return (y - south)/dlat + (1-jdatastart+1);
  }
  else{
    return (north-y)/dlat + 1-jdatastart+1;
  }
}


