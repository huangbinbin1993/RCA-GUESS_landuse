#include"LakeData.h"
#include"ModelData.h"
#include"froutines.h"
#include"myRoutines.h"
#include"localParallel.h"
#include<iostream>
#include<math.h>
using namespace std;
#define NOTDEFINED -666

using namespace std;

int nlon(){ //number of points in e-w direction
  return 43200;
}
int nlat(){ //number of points in n-s direction
  return 21600;
}
double dlonI(){//grid size in e-w direction
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -180.0+halfdeltai;
  double b = 180.0-halfdeltai;
  return (b-a)/(lonPoints-1);
}
double dlatI(){//grid size in n-s direction
  double latPoints = (double)nlat();
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -90.0+halfdeltai;
  double b = 90.0-halfdeltai;
  return (b-a)/(latPoints-1);
}

Domain* defineBoundingDomain(char* name,Domain *modelDomain){
  double halfdeltai = 180.0/((double)nlon());
  double w = -180.0+halfdeltai, s = -90.0+halfdeltai;
  double east = w + (nlon()-1)*dlonI();
  double north = s + (nlat()-1)*dlatI();


  bool global = true;
  Domain *GTOPO = new Domain("EcoclimapOriginal",nlon(),nlat(),nlon(),nlat(),w,s,
 			    dlonI(),dlatI(),global,1,1,modelDomain->polon,
			     modelDomain->polat,
			    true,false,-1); //the 'full' global domain

  int isOnPole=0;
  double polon = modelDomain->polon, polat = modelDomain->polat;
  double lonPole = 0.0;
  double nPole = 90.0;
  double sPole = -90.0;
  double rotatedNPole[2],rotatedSPole[2];
  int one=1;
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  reg2rot_(&lonPole,&sPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);
  double lonxRot[4],latyRot[4];
  for(int j=1;j<modelDomain->klat && isOnPole==0;j++)
    for(int i=1;i<modelDomain->klon && isOnPole==0;i++){
      modelDomain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space
      if(modelDomain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot)){
	isOnPole = 1;
      }
      if(modelDomain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot)){
	isOnPole = -1;
      }
    }
  
  double box[4]={666.0,-666.0,666.0,-666.0};//minLon,maxLon,minLat,maxLat
  for(int j=0;j<=modelDomain->klat+1;j++){
    for(int i=0;i<=modelDomain->klon+1;i++){
      box[0] = MIN(box[0],modelDomain->lon(i,j));
      box[1] = MAX(box[1],modelDomain->lon(i,j));
      box[2] = MIN(box[2],modelDomain->lat(i,j));
      box[3] = MAX(box[3],modelDomain->lat(i,j));
    }
  }
  if(isOnPole!=0){
    box[0] = -180.0;
    box[1] = 180.0;
    if(isOnPole==1)
      box[3] =90.0; 
    else
      box[2] =-90.0;
  }


  int ibounds[2]={nlon(),1},jbounds[2]={1,nlat()}; //index in GTOPO

  while(GTOPO->lon(ibounds[0],1) > box[0] && ibounds[0]>1)
    ibounds[0]--;
  while(GTOPO->lon(ibounds[1],1) < box[1] && ibounds[1]<nlon())
    ibounds[1]++;
  while(GTOPO->lat(1,jbounds[0]) > box[2] && jbounds[0]<nlat())  //lat decreases with growing j
    jbounds[0]++;
  while(GTOPO->lat(1,jbounds[1]) < box[3]  && jbounds[1]>1)
    jbounds[1]--;


  if(isOnPole!=0){
    ibounds[0]=1;
    ibounds[1]=nlon();
    if(isOnPole==1)
      jbounds[1] = 1;
    else
      jbounds[0] = nlat();
  }
  int tmp = jbounds[0];
  jbounds[0] = jbounds[1];
  jbounds[1] = tmp;

  for(int k=0;k<2;k++){
    ibounds[k] = MIN(MAX(1,ibounds[k]),nlon());
    jbounds[k] = MIN(MAX(1,jbounds[k]),nlat());
  }
 
 
  double west  = GTOPO->west;
  double south = GTOPO->south;
  
 int klon = ibounds[1]-ibounds[0]+1;
 int klon_global = klon;
  
 int klat = jbounds[1]-jbounds[0]+1;
 int klat_global = klat;
 
 int halo = -1;//To mark that we are NOT using any halo!
 bool ipos=true;
 bool jpos=false;
 int idatastart = ibounds[0];  
 int jdatastart = jbounds[0];

 Domain *eco = new Domain(name,klon_global,klat_global,klon,klat,west,south,dlonI(),dlatI(),global,idatastart,
			  jdatastart,modelDomain->polon,modelDomain->polat,ipos,jpos,halo);
 eco->north = GTOPO->north;
 eco->south = eco->north - (jdatastart-1+klat_global-1)*dlatI();

 delete GTOPO;
 return eco;
}


//Tests the gtopo30reader
int main( int argc, char** argv ){
   int myid,nproc;
   initialize_parallel(argc,argv,myid,nproc);

   double south,west,dlon,	dlat,	polon,	polat;
   int klon_global,  klat_global,    klev_global;
   int idatastart,jdatastart;
   int halo=1;
   read_domain_(&south, &west, &dlon, &dlat, &polon, &polat,&klon_global,&klat_global,&klev_global);

   int klon = klon_global, klat=klat_global;
   decompose_domain(klon_global,klat_global,klon,klat,idatastart,jdatastart,halo,myid,nproc);
   float polonF=(float)polon, polatF=(float)polat;
   float southF =(float)south, westF=(float)west;
   float dlonF=(float)dlon, dlatF=(float)dlat;
   int klonF=klon,klatF=klat,klonGlobal=klon_global,klatGlobal=klat_global;


   Domain *modelDomain = new Domain("modelLake",klon,klat,klon,klat,
				    west,south,dlon,dlat,false,1,1,polon,polat);
   Domain *ecoDomain = defineBoundingDomain("ecoclimapOriginalLAM",modelDomain);

   LakeData *lake = new LakeData(nlon(),nlat(),ecoDomain);
   lake->readLakeFile("/nobackup/rossby14/sm_psamu/ECOCLIMAP/RUN/GLOBAL_001/GlobalLakeDepth.001");
   //lake->writeFile("LAM_lakeDepth_001.bin");
   ModelData *model = new ModelData(modelDomain,6);
   model->loadData(1,lake,0.0);//loads 6 componenst
   model->writeFile(1,"deepDepth");
   model->writeFile(2,"mediumDepth");
   model->writeFile(3,"shallowDepth");
   model->writeFile(4,"deepFrac");
   model->writeFile(5,"mediumFrac");
   model->writeFile(6,"shallowFrac");
   close_parallel();
   return 0;
}




















