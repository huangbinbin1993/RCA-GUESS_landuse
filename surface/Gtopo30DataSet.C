#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif
#include"Gtopo30DataSet.h"
#include"myRoutines.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<math.h>
#include"froutines.h"
using namespace std;


Gtopo30DataSet::Gtopo30DataSet(char* path_in,ModelData *model,MPI_Comm localComm):
  numberOfSets(33){
  mapUsed = new bool[numberOfSets];
  filenames = new char*[numberOfSets];
  imap = new int[numberOfSets];
  jmap = new int[numberOfSets];

  //path should be given as an input parameter to Gtopo30DataSet constructor
  int len = static_cast<int>(strlen(path_in));
  path = new char[len+1];
  strcpy(path,path_in);
  connectFileNames();

  defineDomains(model,localComm);  //this sets domains to all subdomains and also defines which are used
  setImapJmap();         //imap,jmap for all subdomains

  Domain *boundingDomain = defineBoundingDomain("gtopoOriginalLAM",model,localComm);
  for(int i=0;i<numberOfSets;i++)
    mapUsed[i]=false;
  
  for(int j=1;j<=boundingDomain->klat;j++)
    for(int i=1;i<=boundingDomain->klon;i++){
      int map = getMapNumber(boundingDomain->lon(i,j), boundingDomain->lat(i,j));
      mapUsed[map] = true;
    }
  //makeRectangularMaps(); //??

  gtopo = new Data(boundingDomain); //allocates and initializes gtopo with -666
  loadTiledData(localComm);
}

int Gtopo30DataSet::nlon(){ //number of points in e-w direction
  return 43200;
}
int Gtopo30DataSet::nlat(){ //number of points in n-s direction
  return 21600;
}


double Gtopo30DataSet::dlonI(){//grid size in e-w direction
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -180.0+halfdeltai;
  double b = 180.0-halfdeltai;
  return (b-a)/(lonPoints-1);
}
double Gtopo30DataSet::dlatI(){//grid size in n-s direction
  double latPoints = (double)nlat();
  double lonPoints = (double)nlon();
  double halfdeltai = 180.0/lonPoints;
  double a = -90.0+halfdeltai;
  double b = 90.0-halfdeltai;
  return (b-a)/(latPoints-1);
}



Domain* Gtopo30DataSet::defineBoundingDomain(char* name,ModelData *model,MPI_Comm localComm){
  double halfdeltai = 180.0/((double)nlon());
  double w = -180.0+halfdeltai, s = -90.0+halfdeltai;
  double east = w + (nlon()-1)*dlonI();
  double north = s + (nlat()-1)*dlatI();


  bool global = true;
  Domain *GTOPO = new Domain("GtopoOriginal",nlon(),nlat(),nlon(),nlat(),w,s,
			     dlonI(),dlatI(),localComm,global,1,1,model->polon(),model->polat(),
			     true,false,-1); //the 'full' global domain
  double polon = model->polon(), polat = model->polat();
  double lonPole = 0.0;
  double nPole = 90.0;
  double sPole = -90.0;
  double rotatedNPole[2],rotatedSPole[2];
  int one=1;
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  reg2rot_(&lonPole,&sPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);
  double lonxRot[4],latyRot[4];
  int isOnPole=0; 
  for(int j=1;j<=model->klat() && isOnPole==0;j++)
    for(int i=1;i<=model->klon() && isOnPole==0;i++){
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
//   for(int j=0;j<=model->klat()+1;j++){
//     for(int i=0;i<=model->klon()+1;i+=model->klon()+1){
//       box[0] = MIN(box[0],model->lon(i,j));
//       box[1] = MAX(box[1],model->lon(i,j));
//       box[2] = MIN(box[2],model->lat(i,j));
//       box[3] = MAX(box[3],model->lat(i,j));
//     }
//   }

  if(isOnPole!=0){
    box[0] = -180.0;
    box[1] = 180.0;
    if(isOnPole==1)
      box[3] =90.0; 
    else
      box[2] =-90.0;
  }

//   box[0] = MIN(MAX(GTOPO->lon(1,1),box[0]),GTOPO->lon(nlon(),1));
//   box[1] = MIN(MAX(GTOPO->lon(1,1),box[1]),GTOPO->lon(nlon(),1));
//   box[2] = MIN(MAX(GTOPO->lat(1,1),box[2]),GTOPO->lat(1,nlat()));
//   box[3] = MIN(MAX(GTOPO->lat(1,1),box[3]),GTOPO->lat(1,nlat()));


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

  Domain *eco = new Domain(name,klon_global,klat_global,klon,klat,west,south,dlonI(),dlatI(),localComm,global,idatastart,
			   jdatastart,model->polon(),model->polat(),ipos,jpos,halo);
  eco->north = GTOPO->north;
  eco->south = eco->north - (jdatastart-1+klat_global-1)*dlatI();


  delete GTOPO;
  return eco;
}

void Gtopo30DataSet::defineDomains(ModelData *model,MPI_Comm localComm){

  nDomains=numberOfSets;
  
  domains = new Domain*[numberOfSets];
  //g30topoData = new Gtopo30Data*[numberOfSets];

  bool global = true;
  bool ipos=true;
  bool jpos=false;
  for(int map=0;map<numberOfSets;map++){
    double dLon,dLat,West,South;
    int Klon,Klat;
    int len = static_cast<int>(strlen(getFileName(map)));
    int lenPath = static_cast<int>(strlen(path));
    char* name = new char[len-lenPath+1];
    strcpy(name,getFileName(map)+lenPath);
    int lenName = static_cast<int>(strlen(name));
    name[lenName-3] = 'O';
    name[lenName-2] = 'R';
    name[lenName-1] = 'G';
    
    readHeaderFile(getFileName(map),dLon,dLat,Klon,Klat,West,South);
    //int iLam=1,jLam=1;
    int halo=-1;
    int idatastart=1,jdatastart=1;
    domains[map] = new Domain(name,Klon,Klat,Klon,Klat, //The domain for the whole subdomain
			      West,South,dLon,dLat,localComm,
			      global,idatastart,jdatastart,
			      model->polon(),model->polat(),ipos,jpos,halo);
    delete[] name;
  }
}




void Gtopo30DataSet::makeRectangularMaps(){
   
  
  //                                         if last row is not used 9cols       otherwise 3cols       full dataset
  //make the tile of maps into a rectangle... [0  1  2  3  4  5  6  7  8 ]      [0 -2  3 - 5  6-8 ]  [0  1  2  3  4  5  6  7  8 ]
  //                                          [9  10 11 12 13 14 15 16 17]	[9 -11 12-14 15-17]  [9  10 11 12 13 14 15 16 17]
  //                                          [18 19 20 21 22 23 24 25 26]	[18-20 21-23 24-26]  [18 19 20 21 22 23 24 25 26]
  //                                                                    	[27-28 29-30 31-32]  [27    28 29    30 31    32]
  // s.t. the gtopo can be indexed with (i,j)
  bool rowUsed[4]={false,false,false,false};

  for(int c=0;c<numberOfSets;c++)
    if(mapUsed[c]){
      rowUsed[c/9]=true;
    }
  int nRows=0;
  for(int c=0;c<4;c++)
    if(rowUsed[c])nRows++;

  
  int nc = rowUsed[3] ? 3:9;
  bool *colUsed =  new bool[nc];
  for(int c=0;c<nc;c++)
    colUsed[c]=false;

  int imin=nc,imax=0;
  for(int c=0;c<nc;c++)
    if(nc==9){
      if(mapUsed[c] || mapUsed[c+9] || mapUsed[c+18]){
	imin = MIN(imin,c);
	imax = MAX(imax,c);
      }
    }
    else{
      if(mapUsed[3*c] || mapUsed[3*c+1] || mapUsed[3*c+2] ||
	 mapUsed[3*c+9] || mapUsed[3*c+10] || mapUsed[3*c+11] ||
	 mapUsed[3*c+18] || mapUsed[3*c+19] || mapUsed[3*c+20] ||
	 mapUsed[2*c+27] || mapUsed[2*c+28] ){
	imin = MIN(imin,c);
	imax = MAX(imax,c);
      }
    }

  int jmin=nRows,jmax=0;
  for(int row=0;row<4;row++){
    if(row!=3){
      for(int k=0;k<9;k++)
	if(mapUsed[9*row+k]){
	  jmin = MIN(jmin,row);
	  jmax = MAX(jmax,row);
	}
    }
    else
      for(int k=0;k<6;k++)
	if(mapUsed[9*row+k]){
	  jmin = MIN(jmin,row);
	  jmax = MAX(jmax,row);
	}
  }
  
  int j=jmin;
  while( j<=jmax){
    int i=imin;
    while(i<=imax){
      int ne = j==3 ? 6:9;
      mapUsed[i+j*ne]=true;
      i++;
    }
    j++;
  }
  
#ifdef DEBUG
  //    cout <<"After manipulations"<<endl;
  //    cout <<imin<<" "<<imax<<endl;
  //    cout <<jmin<<" "<<jmax<<endl;
  //    for(int c=0;c<numberOfSets;c++){
  //      cout <<mapUsed[c]<<"  ";
  //      if((c+1)%9==0)
  //        cout <<endl;
  //    } 
  //    cout <<endl;
#endif
}

void Gtopo30DataSet::loadTiledData(MPI_Comm localComm){
  //determine which dataset to read from and store in maps
  int rank=0;
#ifdef MPI_SRC
  MPI_Comm_rank(localComm, &rank );
#endif
  int N=0;
  bool terminate=false;

  for(int map=0;map<numberOfSets;map++){
    if(mapUsed[map]){
      Gtopo30Data *g30topoData  = new Gtopo30Data(domains[map],domains[map]->klon_global,
						  domains[map]->klat_global,map,localComm);
      g30topoData->readFile(getFileName(map));//read the whole map

      int im = imap[map];	
      int jm = jmap[map];
      int il = gtopo->idatastart();
      int jl = gtopo->jdatastart();
	
      for(int j=1;j<=domains[map]->klat;j++)
	for(int i=1;i<=domains[map]->klon;i++){
	  //i+im är på den globala kartan
	  int ilam = i+im-il+1;
	  int jlam = j+jm-jl+1;
 
	  if(ilam>0 && ilam <= gtopo->klon() &&   jlam>0 && jlam <= gtopo->klat()){
	    if( fabs(gtopo->lon(ilam,jlam)-g30topoData->lon(i,j))>0.1*g30topoData->dlon() ||
		fabs(gtopo->lat(ilam,jlam)-g30topoData->lat(i,j))>0.1*g30topoData->dlon()){
	      cout<<gtopo->lon(ilam,jlam)<<"   "<<g30topoData->lon(i,j)<<endl;
	      terminate = true;
	      //terminate_program(3);
	    }
	    //double dat = g30topoData[map]->data(i,j);
	    double dat = g30topoData->data(i,j);
	    if(fabs(dat+9999.0)<0.01)
	      dat = 0.0; //seapoints are given the value of 0
	    gtopo->data(ilam, jlam) = dat;
	    N++;
	  }
	}
      delete g30topoData;
    }
  }
  //Are all points read in?
  if(N!=(gtopo->klon()*gtopo->klat()) ||terminate){
    if(terminate)cout <<"terminating due to index mismatch"<<endl;
    cout <<"N="<<N<< "!="<<gtopo->klon()*gtopo->klat()<<"  "<<numberOfSets<<endl;

    //gtopo->domain->print("Are all point read in Gtopo30DataSet.C 350",localComm);
    int *MMAP = new int[numberOfSets];
    for(int map=0;map<numberOfSets;map++){
      if(mapUsed[map])cout << map <<"  ";
      MMAP[map] = -1;
    }
    cout <<endl;
    cout <<"i    j      klon()   klat()"<<__FILE__<<__LINE__<<endl;
    for(int j=1;j<=gtopo->klat();j++)
      for(int i=1;i<=gtopo->klon();i++){
 	if(fabs(gtopo->data(i,j)+666.0)<0.01)
	  MMAP[getMapNumber(gtopo->lon(i,j),gtopo->lat(i,j))] = 1;
	//cout <<"map "<<getMapNumber(gtopo->lon(i,j),gtopo->lat(i,j))<<"  "; 
	// 	  cout <<i<<"  "<<j<<"  "<<gtopo->klon()<<"  "<<gtopo->klat()<<endl;
      }
    cout <<"Data from the following maps is not loaded correctly::";
    for(int map=0;map<numberOfSets;map++){
      if(MMAP[map]==1)
	cout <<" ("<< map<<","<<(mapUsed[map] ? "T" : "F")<<")  ," ;
    }
    cout <<endl;
    terminate_program(2);
  }
}

void Gtopo30DataSet::setImapJmap(){
  imap[0]=0;
  for(int c=1;c<9;c++)
    imap[c] = domains[c-1]->klon + imap[c-1];
  for(int r=1;r<3;r++)
    for(int c=0;c<9;c++)
      imap[c+r*9] = imap[c]; 
  imap[27]=0;
  for(int c=28;c<33;c++)
    imap[c] = domains[c-1]->klon + imap[c-1];


  jmap[0]=0;
  for(int r=1;r<4;r++)
    jmap[r*9] = domains[(r-1)*9]->klat + jmap[(r-1)*9];
  for(int r=0;r<3;r++)
    for(int c=1;c<9;c++)
      jmap[c+r*9] = jmap[r*9];
  jmap[28]=jmap[29]=jmap[30]=jmap[31]=jmap[32]=jmap[27];

}

Domain* Gtopo30DataSet::defineTiledDomain(char* name,ModelData *model,Domain *oneDomain,MPI_Comm localComm){
  int N = model->klon(), M = model->klat();
  int ibounds[2] = {oneDomain->klon,1}; //these are indecies in the Tiled Domain
  int jbounds[2] = {oneDomain->klat,1};

  for(int j=0;j<=N+1;j++)
    for(int i=0;i<=M+1;i++){
      double lon = model->lon(i,j);
      double lat = model->lat(i,j);
      int map = getMapNumber(lon, lat);
      ibounds[0] = MIN(ibounds[0],imap[map]+1);
      jbounds[0] = MIN(jbounds[0],jmap[map]+1);
      ibounds[1] = MAX(ibounds[1],domains[map]->klon+imap[map]);
      jbounds[1] = MAX(jbounds[1],domains[map]->klat+jmap[map]);
    }
   
  double west  = oneDomain->lon(ibounds[0],jbounds[0]);
  double south = oneDomain->lat(ibounds[0],jbounds[0]);
   
  for(int j=0;j<2;j++)
    for(int i=0;i<2;i++){
      west  = MIN(west, oneDomain->lon(ibounds[i],jbounds[j]));
      south = MIN(south, oneDomain->lat(ibounds[i],jbounds[j]));
    }
   
   
  int klon = ibounds[1]-ibounds[0]+1;
  int klon_global = klon;
   
  int klat = jbounds[1]-jbounds[0]+1;
  int klat_global = klat;
   
  int halo = -1;//To mark that we are NOT using any halo!
  bool ipos=true;
  bool jpos=false;
  //    int iLAMstart = ibounds[0];  
  //    int jLAMstart = jbounds[0];
  //    int idatastart = 1;
  //    int jdatastart = 1;
  int idatastart = ibounds[0];  
  int jdatastart = jbounds[0];
  bool global=true;
  Domain *eco = new Domain(name,klon_global,klat_global,klon,klat,
			   west,south,oneDomain->dlon,oneDomain->dlat,
			   localComm,global,idatastart,
			   jdatastart,model->polon(),model->polat(),
			   ipos,jpos,halo);
   
  return eco;
   
}

Domain* Gtopo30DataSet::glueTilesIntoOne(MPI_Comm localComm){
  //the following information can be taken from any domain
  double dlon = domains[0]->dlon;
  double dlat = domains[0]->dlat;
  double polon = domains[0]->polon;
  double polat = domains[0]->polat;
  bool ipos = domains[0]->ipos;
  bool jpos = domains[0]->jpos;

  double south = 999.0;//domains[0]->south;
  double west = 999.0;//domains[0]->west;
  
  int id_south;
  int id_west;
  for(int i=0;i<numberOfSets;i++)
    if(mapUsed[i]){
      if(domains[i]->south < south){
	south =domains[i]->south;
	id_south = i;
      }
      if(domains[i]->west < west){
	west =domains[i]->west;
	id_west = i;
      }
    }
  

  double eps=0.0001;
  int klon=0,klon_global=0;
  for(int i=0;i<numberOfSets;i++)
    if(mapUsed[i] &&    fabs(domains[id_south]->south - domains[i]->south)<eps){
      klon += domains[i]->klon;
      klon_global += domains[i]->klon_global;
    }

  int klat=0,klat_global=0;
  for(int i=0;i<numberOfSets;i++)
    if(mapUsed[i] &&fabs(domains[id_west]->west - domains[i]->west)<eps){
      klat += domains[i]->klat;
      klat_global += domains[i]->klat_global;
    }
  
  bool global=true;
  int halo=-1;
  
  Domain *one = new Domain("Gtopo30",klon_global,klat_global,klon,klat,west,south,dlon,dlat,localComm,global,1,1,
			   polon,polat,ipos,jpos,halo);
  return one;
}


void Gtopo30DataSet::connectFileNames(){
  for(int i=0;i<numberOfSets;i++)
    connectFileName(i);
}


void Gtopo30DataSet::connectFileName(int mapNumb){
  //This order of files is MANDATORY!!!
  const int fileLen=11;//19
  //char filename[19];
  char filename[fileLen+1];
  switch(mapNumb){
  case 0:
    //strcpy(filename,"w180n90/W180N90.DEM");
    strcpy(filename,"W180N90.DEM");
    break;
  case 1:
    //strcpy(filename,"w140n90/W140N90.DEM");
    strcpy(filename,"W140N90.DEM");
    break;
  case 2:
    //strcpy(filename,"w100n90/W100N90.DEM");
    strcpy(filename,"W100N90.DEM");
    break;
  case 3:
    //strcpy(filename,"w060n90/W060N90.DEM");
    strcpy(filename,"W060N90.DEM");
    break;
  case 4:
    //strcpy(filename,"w020n90/W020N90.DEM");
    strcpy(filename,"W020N90.DEM");
    break;
  case 5:
    //strcpy(filename,"e020n90/E020N90.DEM");
    strcpy(filename,"E020N90.DEM");
    break;
  case 6:
    //strcpy(filename,"e060n90/E060N90.DEM");
    strcpy(filename,"E060N90.DEM");
    break;
  case 7:
    //strcpy(filename,"e100n90/E100N90.DEM");
    strcpy(filename,"E100N90.DEM");
    break;
  case 8:
    //strcpy(filename,"e140n90/E140N90.DEM");
    strcpy(filename,"E140N90.DEM");
    break;
  case 9:
    //strcpy(filename,"w180n40/W180N40.DEM");
    strcpy(filename,"W180N40.DEM");
    break;
  case 10:
    //strcpy(filename,"w140n40/W140N40.DEM");
    strcpy(filename,"W140N40.DEM");
    break;
  case 11:
    //strcpy(filename,"w100n40/W100N40.DEM");
    strcpy(filename,"W100N40.DEM");
    break;
  case 12:
    //strcpy(filename,"w060n40/W060N40.DEM");
    strcpy(filename,"W060N40.DEM");
    break;
  case 13:
    //strcpy(filename,"w020n40/W020N40.DEM");
    strcpy(filename,"W020N40.DEM");
    break;
  case 14:
    //strcpy(filename,"e020n40/E020N40.DEM");
    strcpy(filename,"E020N40.DEM");
    break;
  case 15:
    //strcpy(filename,"e060n40/E060N40.DEM");
    strcpy(filename,"E060N40.DEM");
    break;
  case 16:
    //strcpy(filename,"e100n40/E100N40.DEM");
    strcpy(filename,"E100N40.DEM");
    break;
  case 17:
    //strcpy(filename,"e140n40/E140N40.DEM");
    strcpy(filename,"E140N40.DEM");
    break;
  case 18:
    //strcpy(filename,"w180s10/W180S10.DEM");
    strcpy(filename,"W180S10.DEM");
    break;
  case 19:
    //strcpy(filename,"w140s10/W140S10.DEM");
    strcpy(filename,"W140S10.DEM");
    break;
  case 20:
    //strcpy(filename,"w100s10/W100S10.DEM");
    strcpy(filename,"W100S10.DEM");
    break;
  case 21:
    //strcpy(filename,"w060s10/W060S10.DEM");
    strcpy(filename,"W060S10.DEM");
    break;
  case 22:
    //strcpy(filename,"w020s10/W020S10.DEM");
    strcpy(filename,"W020S10.DEM");
    break;
  case 23:
    //strcpy(filename,"e020s10/E020S10.DEM");
    strcpy(filename,"E020S10.DEM");
    break;
  case 24:
    //strcpy(filename,"e060s10/E060S10.DEM");
    strcpy(filename,"E060S10.DEM");
    break;
  case 25:
    //strcpy(filename,"e100s10/E100S10.DEM");
    strcpy(filename,"E100S10.DEM");
    break;
  case 26:
    //strcpy(filename,"e140s10/E140S10.DEM");
    strcpy(filename,"E140S10.DEM");
    break;
  case 27:
    //strcpy(filename,"w180s60/W180S60.DEM");
    strcpy(filename,"W180S60.DEM");
    break;
  case 28:
    //strcpy(filename,"w120s60/W120S60.DEM");
    strcpy(filename,"W120S60.DEM");
    break;
  case 29:
    //strcpy(filename,"w060s60/W060S60.DEM");
    strcpy(filename,"W060S60.DEM");
    break;
  case 30:
    //strcpy(filename,"w000s60/W000S60.DEM");
    strcpy(filename,"W000S60.DEM");
    break;    
  case 31:
    //strcpy(filename,"e060s60/E060S60.DEM");
    strcpy(filename,"E060S60.DEM");
    break;
  case 32:
    //strcpy(filename,"e120s60/E120S60.DEM");
    strcpy(filename,"E120S60.DEM");
    break;
  default:
    printf("No map like this %i\n",mapNumb);
  }
  int len = static_cast<int>(strlen(path));
  filenames[mapNumb] = new char[len+fileLen+1];
  strcpy(filenames[mapNumb],path);
  strcat(filenames[mapNumb],filename);
  //cout <<"filenames["<<mapNumb<<"]="<<filenames[mapNumb]<<endl;
}


int Gtopo30DataSet::getMapNumber(double lon,double lat){
  //make sure lon=[-180,180], lat=[-90,90]

  if(lat >=40.0 && lat < 90.0){
    if(lon >=-180.0 && lon < -140.0)//W180N90
      return 0;
    if(lon >=-140.0 && lon < -100.0)//W140N90            
      return 1;
    if(lon >=-100.0 && lon < -60.0)//W100N90            
      return 2;
    if(lon >= -60.0  && lon < -20.0)//W060N90            
      return 3;
    if(lon >=-20.0 && lon <  20.0)//W020N90             
      return 4;
    if( lon >=20.0 && lon < 60.0)//E020N90
      return 5;
    if(lon >=60.0 && lon < 100.0)//E060N90              
      return 6;
    if(lon >=100.0 && lon < 140.0)//E100N90             
      return 7;
    if(lon >=140.0 && lon <= 180.0)//E140N90
      return 8;
  }
  if(lat >=-10.0 && lat < 40.0){
    if(lon>=-180.0 && lon<-140.0) //W180N40  
      return 9;
    if(lon>=-140.0 && lon< -100.0)//W140N40  
      return 10;
    if(lon>=-100.0 && lon< -60.0) //W100N40  
      return 11;
    if(lon>=-60.0  && lon< -20.0) //W060N40  
      return 12;
    if(lon>=-20.0  && lon<  20.0) //W020N40  
      return 13;
    if(lon>=20.0   && lon<60.0)  //E020N40  
      return 14;
    if(lon>=60.0   && lon<100.0) //E060N40  
      return 15;
    if(lon>=100.0  && lon<140.0) //E100N40  
      return 16;
    if(lon>=140.0  && lon<=180.0) //E140N40  
      return 17;
  }
  if(lat>=-60.0 && lat<-10.0){
    if(lon>=-180.0 && lon<-140.0)//W180S10             
      return 18;		              
    if(lon>=-140.0 && lon< -100.0)//W140S10             
      return 19;		              
    if(lon>=-100.0 && lon< -60.0)//W100S10             
      return 20;		              
    if(lon>=-60.0  && lon< -20.0)//W060S10             
      return 21;		              
    if(lon>=-20.0  && lon<  20.0)//W020S10             
      return 22;		   
    if(lon>=20.0   && lon<60.0)//E020S10    
      return 23;		   
    if(lon>=60.0   && lon<100.0)//E060S10   
      return 24;		   
    if(lon>=100.0  && lon<140.0)//E100S10   
      return 25;		   
    if(lon>=140.0  && lon<=180.0)//E140S10   
      return 26;
  }
   
  if(lat>=-90.0 && lat<-60.0){
    if(lon>=-180.0 && lon<-120.0)//W180S60  
      return 27;
    if(lon>=-120.0&& lon<-60.0)//W120S60  
      return 28;
    if(lon>=-60.0 && lon<0.0)//W060S60  
      return 29;
    if(lon>=0.0 && lon<60.0)//W000S60  
      return 30;
    if(lon>=60.0 && lon<120.0)//E060S60  
      return 31;
    if(lon>=120.0 && lon<=180.0 )//E120S60  
      return 32;
  }
  return -1;
  //  -90      -60       -180     180   ANTARCPS 
}



void Gtopo30DataSet::readHeaderFile(char* fileName,double& g30_dlon,double& g30_dlat,
				    int& ncols,int& nrows,double& g30_lon,double& g30_lat){
  int len = static_cast<int>(strlen(fileName));
  char* headerFile = new char[len+1];

  strcpy(headerFile,fileName);
  headerFile[len-3] = 'H';
  headerFile[len-2] = 'D';
  headerFile[len-1] = 'R';
  
  FILE	*finp;
  if ((finp = fopen (headerFile, "r")) == NULL ) {
    fprintf (stderr, "g30read:  Could not open file %s \n",headerFile);
    cout <<__FILE__<<"   "<<__LINE__<<endl;
    terminate_program(2);
  }

  //JUMP OVER 2 ROWS
  int cnt=0, safety_limit=10000;
  for(int i=0;i<2;i++){
    cnt=0;
    while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  }
  
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%i",&nrows)){
    cout <<"Error reading nrows"<<endl;
  }

  cnt=0;
  while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%i",&ncols)){
    cout <<"Error reading ncols"<<endl;
  }

  cnt=0;
  while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};

  //jump over 6 rows
  for(int i=0;i<6;i++){
    cnt=0;
    while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  }
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%lf",&g30_lon)){
    cout<<"Error reading g30_lon"<<endl;
  }
  cnt=0;
  while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%lf",&g30_lat)){
    cout<<"Error reading g30_lat"<<endl;
  }
  cnt=0;
  while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%lf",&g30_dlon)){
    cout<<"Error reading g30_dlon"<<endl;
  }
  cnt=0;
  while( getc(finp) != '\n' && cnt < safety_limit ) {cnt++;};
  cnt=0;
  while( getc(finp) != ' ' && cnt < safety_limit ) {cnt++;};
  if(!fscanf(finp,"%lf\n",&g30_dlat)){
    cout << "Error reading g30dlat"<<endl;
  }
  fclose(finp);
  //adjusting for the weird header file....
  g30_lat -= (nrows-1)*g30_dlat;

  delete[] headerFile;
}


Gtopo30DataSet::~Gtopo30DataSet(){
#ifdef DEBUG
  cout << "delete Gtopo30DataSet called"<<endl;
#endif
  //   for(int k=0;k<numberOfSets;k++){
  //     if(mapUsed[k])
  //       delete g30topoData[k];
  //   }
  for(int k=0;k<numberOfSets;k++){
    delete[] filenames[k];
  }
  delete[] filenames;
  delete[] domains;
  //  delete[] g30topoData;
  delete[] mapUsed;
  delete[] path;

  delete[] imap;
  delete[] jmap;

  //delete LAMDomain;
  delete gtopo;
}
