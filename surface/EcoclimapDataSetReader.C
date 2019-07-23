#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"EcoclimapDataSetReader.h"
#include"froutines.h"
#include"texture.h"
#include"Domain.h"
#include"myRoutines.h"
#include<string.h>
#include<iostream>
#include<math.h>
#define NOT_USED -1
using namespace std;


EcoclimapDataSetReader::~EcoclimapDataSetReader(){
#ifdef DEBUG
  cout << "delete EcoclimapDataSetReader called"<<endl;
#endif
  delete annualContainer;
  delete ecodataSet;
  delete model;//new for now
  delete[] timeDep;
  delete[] var;
  delete[] readFromFile;
}

EcoclimapDataSetReader::EcoclimapDataSetReader(char* path, int month, Domain *domain,MPI_Comm localComm){
  nc = 44; //the lake database generates 6 different model fields

  model = new ModelData(domain,43+6); //is deleted

#ifdef DEBUG
  double minlat=666,maxlat=-666;
  double minlon=666,maxlon=-666;
  model->findMinMax(minlon,maxlon,minlat,maxlat);
  cout <<"Ecoclimap According to program in regular space::"<<endl;
  cout<<"lon =["<<minlon<<" , "<<maxlon<<"]"<<endl;
  cout<<"lat =["<<minlat<<" , "<<maxlat<<"]"<<endl;
#endif

  Domain *ecoDomain = defineBoundingDomain("ecoclimapOriginalLAM",localComm);  //This defines the ecoclimap domain

  ecodataSet = new EcoclimapDataSet(path,ecoDomain,nlon(),nlat());//is deleted
  

  Domain *domainA = new Domain(domain);//used to init annualContainer
  annualContainer = new ModelData(domainA,12*3);//3=number of tiles, 12=number of months is deleted
  
  var = new int[nc];//is deleted
  readFromFile = new bool[nc];//is deleted
  timeDep = new bool[nc]; //is deleted
  

  initialized = false;
#ifdef DEBUG
  cout <<"getting ecoclimapdata..."<<endl;
#endif  
  getEcoclimapData(month,localComm);
}



Domain* EcoclimapDataSetReader::defineBoundingDomain(char* name,MPI_Comm localCommIn){
  double halfdeltai = 180.0/((double)nlon());
  double w = -180.0+halfdeltai, s = -90.0+halfdeltai;
  double east = w + (nlon()-1)*dlonI();
  double north = s + (nlat()-1)*dlatI();


  bool global = true;
  Domain *GTOPO = new Domain("EcoclimapOriginal",nlon(),nlat(),nlon(),nlat(),w,s,
			     dlonI(),dlatI(),localCommIn,global,1,1,model->polon(),model->polat(),
			    true,false,-1); //the 'full' global domain

  int isOnPole=0;
  double polon = model->polon(), polat = model->polat();
  double lonPole = 0.0;
  double nPole = 90.0;
  double sPole = -90.0;
  double rotatedNPole[2],rotatedSPole[2];
  int one=1;
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  reg2rot_(&lonPole,&sPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);
  double lonxRot[4],latyRot[4];
  for(int j=1;j<model->klat() && isOnPole==0;j++)
    for(int i=1;i<model->klon() && isOnPole==0;i++){
      model->domain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space
      if(model->domain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot)){
	isOnPole = 1;
      }
      if(model->domain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot)){
	isOnPole = -1;
      }
    }
  
  double box[4]={666.0,-666.0,666.0,-666.0};//minLon,maxLon,minLat,maxLat
  for(int j=0;j<=model->klat()+1;j++){
    for(int i=0;i<=model->klon()+1;i++){
      box[0] = MIN(box[0],model->lon(i,j));
      box[1] = MAX(box[1],model->lon(i,j));
      box[2] = MIN(box[2],model->lat(i,j));
      box[3] = MAX(box[3],model->lat(i,j));
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

 Domain *eco = new Domain(name,klon_global,klat_global,klon,klat,west,south,dlonI(),dlatI(),localCommIn,global,idatastart,
			  jdatastart,model->polon(),model->polat(),ipos,jpos,halo);
 eco->north = GTOPO->north;
 eco->south = eco->north - (jdatastart-1+klat_global-1)*dlatI();

 delete GTOPO;
 return eco;
}




void EcoclimapDataSetReader::getEcoclimapData(int month,MPI_Comm localComm){
  int sandC = 13,clayC=12,soil_carbC=43;

  for(int c=1;c<=nc;c++){//there are a few ops to be saved in this loop
    int cEco = ecodataSet->getVariableNr(c,month);
    var[c-1] = cEco;
    timeDep[c-1] = ecodataSet->isTimeDep(cEco);
    readFromFile[c-1] = ecodataSet->readFromFile(cEco);
  }
    

  bool re_read=false;
  
  if(!initialized)re_read=true;

  for(int c=1;c<=nc;c++){
    if(readFromFile[c-1]){
      if(!timeDep[c-1] && !initialized){  //Read constant data
	if(c==44){ //lakedatabase
	  if(!model->readFile(c+3,"deepDepthLake",localComm) || 
	     !model->readFile(c+4,"mediumDepthLake",localComm) ||
	     !model->readFile(c+5,"shallowDepthLake",localComm) || 
	     !model->readFile(c,"deepFracLake",localComm) || 
	     !model->readFile(c+1,"mediumFracLake",localComm) ||
	     !model->readFile(c+2,"shallowFracLake",localComm) ){
	    ecodataSet->lakeData->readLakeFile(ecodataSet->getFileName(var[c-1]));
	    model->loadData(c,ecodataSet->lakeData,0.0);
	    model->writeFile(c,"deepFracLake",localComm); 
	    model->writeFile(c+1,"mediumFracLake",localComm);
	    model->writeFile(c+2,"shallowFracLake",localComm);
            model->writeFile(c+3,"deepDepthLake",localComm);
	    model->writeFile(c+4,"mediumDepthLake",localComm);
	    model->writeFile(c+5,"shallowDepthLake",localComm);

	  }
	}
	else{
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
	    ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(var[c-1]));
	    if(c==sandC || c==clayC || c==soil_carbC){
	      Domain *tmpDomain = new Domain(ecodataSet->ecoclimapData->domain);
	      EcoclimapData *Fwat = new EcoclimapData(nlon(),nlat(),tmpDomain);
	      
	      Fwat->readEcoFile(ecodataSet->getFileName(var[9]));
	      model->loadData(c,ecodataSet->ecoclimapData,999.0,Fwat);
	      model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	      delete Fwat;//tmpDomain is deleted along with Fwat
	    }
	    else{
	      model->loadData(c,ecodataSet->ecoclimapData,999.0);
	      model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	    }
	  }
	}
      }
      else{//read time-dependent data
	//only if a re-read is needed
	if(re_read){
	  int cEco1 = ecodataSet->getVariableNr(c,month);
	  if(!model->readFile(c,ecodataSet->getPureFileName(cEco1),localComm)){
	    ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
	    if(c==4 ||c==5 || c==6) //z0_t{1,2,3}
	      model->loadData_z0(c,ecodataSet->ecoclimapData,999.0);
	    else
	      model->loadData(c,ecodataSet->ecoclimapData,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(cEco1),localComm);
	  }
	}
      }
    }
  }
  if(!initialized){
    //The derived fields require all/some other fields to be initiated?
    for(int c=1;c<=nc;c++)
      if(!readFromFile[c-1]){ //this means that it cannot be read from file, but rather needs to be computed from the other fields
	switch(c){
	case 35://texture
	  {
	    if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){ 
	      //This assumes that sand and clay are already read....
	      for(int j=1;j<=model->klat();j++)
		for(int i=1;i<=model->klon();i++){
		  double sand = model->data(sandC,i,j);
		  double clay =model->data(clayC,i,j);
		  model->data(c,i,j)=(double)what_texture_(&sand,&clay);
		}
	      model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	    }
	  }
	  break;
	case 36://minlai_t1 lai_t1=annualContainer(1...12)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){ 
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif		    
	    for(int mon=1;mon<=12;mon++){
	      int clai_t1 = 1;
	      int cEco1 = ecodataSet->getVariableNr(clai_t1,mon);
	      if(!annualContainer->readFile(mon,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	      //integrate datasets into model, then
	      //load the minimum into model->data
	    }
	    model->loadMinData(c,annualContainer,1,12,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif	    
	  
	  break;
	case 37://minlai_t2 lai_t2=annualContainer(13...24)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif	
	    for(int mon=1;mon<=12;mon++){
	      int clai_t2 = 2;
	      int cEco1 = ecodataSet->getVariableNr(clai_t2,mon);
	      if(!annualContainer->readFile(mon+12,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon+12,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon+12,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	      //integrate datasets into model, then
	      //load the minimum into model->data
	    }
	    model->loadMinData(c,annualContainer,13,24,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	    MPI_Barrier(localComm);
#endif	    
	  break;
	case 38://minlai_t3 lai_t3=annualContainer(25...36)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif	
	    for(int mon=1;mon<=12;mon++){
	      int clai_t3 = 3;
	      int cEco1 = ecodataSet->getVariableNr(clai_t3,mon);
	      if(!annualContainer->readFile(mon+24,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon+24,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon+24,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	      //integrate datasets into model, then
	      //load the minimum into model->data
	    }
	    model->loadMinData(c,annualContainer,25,36,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	    MPI_Barrier(localComm);
#endif
	  break;
	case 39://maxlai_t1 lai_t1=annualContainer(1...12)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif	
	    for(int mon=1;mon<=12;mon++){
	      int clai_t1 = 1;
	      int cEco1 = ecodataSet->getVariableNr(clai_t1,mon);
	      if(!annualContainer->readFile(mon,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	      //integrate datasets into model, then
	      //load the minimum into model->data
	    }
	    model->loadMaxData(c,annualContainer,1,12,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	    MPI_Barrier(localComm);
#endif
	  
	  break;
	case 40://maxlai_t2 lai_t2=annualContainer(13...24)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
	    for(int mon=1;mon<=12;mon++){
#ifdef MPI_SRC
	      MPI_Bcast( &mon, 1, MPI_INT, 0, localComm );

	      //MPI_Barrier(localComm);//to be sure that that the files that have already begun to written are actually ready
#endif	
	      int clai_t2 = 2;
	      int cEco1 = ecodataSet->getVariableNr(clai_t2,mon);
	      if(!annualContainer->readFile(mon+12,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon+12,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon+12,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	      //integrate datasets into model, then
	      //load the minimum into model->data
	    }
	    model->loadMaxData(c,annualContainer,13,24,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	    MPI_Barrier(localComm);
#endif
	  
	  break;
	case 41://maxlai_t3 lai_t3=annualContainer(25...36)
	  if(!model->readFile(c,ecodataSet->getPureFileName(var[c-1]),localComm)){
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif	
	    for(int mon=1;mon<=12;mon++){
#ifdef MPI_SRC
	      MPI_Bcast( &mon, 1, MPI_INT, 0, localComm );
#endif
	      int clai_t3 = 3;
	      int cEco1 = ecodataSet->getVariableNr(clai_t3,mon);
	      if(!annualContainer->readFile(mon+24,ecodataSet->getPureFileName(cEco1),localComm)){
		ecodataSet->ecoclimapData->readEcoFile(ecodataSet->getFileName(cEco1));
		annualContainer->loadData(mon+24,ecodataSet->ecoclimapData,999.0);
		annualContainer->writeFile(mon+24,ecodataSet->getPureFileName(cEco1),localComm);
	      }
	    }
	    model->loadMaxData(c,annualContainer,25,36,999.0);
	    model->writeFile(c,ecodataSet->getPureFileName(var[c-1]),localComm);
	  }
#ifdef MPI_SRC
	  MPI_Barrier(localComm);
#endif
	  break;
	default:
	  cout <<"this cannot be a derived field "<<c<<endl;
	  break;
	  
	}
      }
  }
  //  initialized = true;
}










//number of points in the global dataset given resolution
int EcoclimapDataSetReader::nlon(){ //number of points in e-w direction
  return 43200;
}
int EcoclimapDataSetReader::nlat(){ //number of points in n-s direction
  return 21600;
}


double EcoclimapDataSetReader::dlonI(){//grid size in e-w direction
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -180.0+halfdeltai;
  double b = 180.0-halfdeltai;
  return (b-a)/(lonPoints-1);
}
double EcoclimapDataSetReader::dlatI(){//grid size in n-s direction
  double latPoints = (double)nlat();
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -90.0+halfdeltai;
  double b = 90.0-halfdeltai;
  return (b-a)/(latPoints-1);
}


