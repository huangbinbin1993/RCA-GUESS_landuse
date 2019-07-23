#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"ModelData.h"
#include"myRoutines.h"
#include"froutines.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include<string.h>
#include<iostream>
#include<stdlib.h>
#include <math.h>



#define UNINI -666.0
#include"Gtopo30Data.h"
#include"EcoclimapData.h"
#include"Gtopo30DataSet.h"
#include"LakeData.h"
using namespace std;



void ModelData::loadData(int c,Data *globaldata,double uninitialized){
  double eps = 0.005;
  int one=1;

  int Klon = domain->klon;
  int Klat = domain->klat;

  double polon = domain->polon, polat = domain->polat;
  double lonPole = 0.0;
  double nPole = 90.0;
  double rotatedNPole[2],rotatedSPole[2];
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  nPole *= -1.0;
  reg2rot_(&lonPole,&nPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);

  
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      double vertx[4], verty[4];
      double lonxRot[4],latyRot[4];
      
      domain->getCell(i,j,vertx,verty);//cell around model->(i,j) in regular space
      domain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space
      int isOnPole=0; 
      if(domain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot))
	isOnPole = 1;
      
      if(domain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot))
	isOnPole = -1;
      
            
      bool crossesDateLine = false;
      int n,m;
      for(n=0,m=3;n<4;m=n++)
	if(vertx[n]*vertx[m]<0.0)
	  if((fabs(vertx[n]) + fabs(vertx[m])) > 350.0)
	    crossesDateLine = true;
      
      //Stuff in regular space
      double minLonModel = MIN(lon(i,j),vertx[0]);
      double maxLonModel = MAX(lon(i,j),vertx[0]);
      
      double minLatModel = MIN(lat(i,j),verty[0]);
      double maxLatModel = MAX(lat(i,j),verty[0]);
      for(int l=1;l<4;l++){
	minLonModel = MIN(vertx[l],minLonModel);
	maxLonModel = MAX(vertx[l],maxLonModel);
	minLatModel = MIN(verty[l],minLatModel);
	maxLatModel = MAX(verty[l],maxLatModel);
      }
	
      double iisD = globaldata->lonindex(minLonModel,minLatModel);      
      double iieD = globaldata->lonindex(maxLonModel,maxLatModel);
	
      double jjsD = globaldata->latindex(minLonModel,minLatModel);
      double jjeD = globaldata->latindex(maxLonModel,maxLatModel);
	
	
      //this is if the index2data is inverse
      if(iisD>iieD){
	double tmp=iisD;
	iisD = iieD;
	iieD = tmp;
      }
      if(jjsD>jjeD){
	double tmp=jjsD;
	jjsD = jjeD;
	jjeD = tmp;
      }

      int iis = MIN(globaldata->klon(),MAX(1,floor(iisD)));
      int iie = MIN(globaldata->klon(),MAX(1,ceil(iieD)));

      int jjs = MIN(globaldata->klat(),MAX(1,floor(jjsD)));
      int jje = MIN(globaldata->klat(),MAX(1,ceil(jjeD)));
	
      if(isOnPole!=0){
	double safe = 0.6*dlon()/globaldata->dlon()+2;
	iisD -= safe;
	jjsD -= safe;
	iieD += safe;
	jjeD += safe;
	if(iisD<1 || iieD>globaldata->klon()){
	  iis = 1;
	  iie = globaldata->klon();
	}
	if(jjsD<1 || jjeD>globaldata->klat()){
	  jjs = 1;
	  jje = globaldata->klat();
	}
      }
	
      if(crossesDateLine){
	iis = 1;
	iie = globaldata->klon();
      }

      bool treated=false;
      double NP=0.0;

    
      for(int jj=jjs;jj<=jje;jj++)
	for(int ii=iis;ii<=iie;ii++){
	  bool inside;
	  double xreg = globaldata->lon(ii,jj);
	  double yreg = globaldata->lat(ii,jj);
	  
	  if(isOnPole!=0 || crossesDateLine){
	    double xrot,yrot;
	    reg2rot_(&xreg,&yreg,&xrot,&yrot,&one,&one,&polon,&polat);
	    inside = domain->isInside(xrot,yrot,lonxRot,latyRot);
	  }
	  else
	    inside = domain->isInside(xreg,yreg,vertx,verty);
	  
	  if(inside){
	    double dat = globaldata->data(ii,jj);
	    treated = true;
	    if(fabs(dat-uninitialized)>eps){ 
	      NP += 1.0;
	      if(useVar){//variance
		if(NP>1)
		  var(c,i,j) = (NP-2.0)/(NP-1.0)*var(c,i,j) + 
		    1.0/(NP)*pow((data(c,i,j)-dat),2);
		else
		  var(c,i,j) = 0.0;	  
	      }
	      data(c,i,j) = ((NP-1.0)*data(c,i,j)+dat)/NP; //mean (upscaling)
	    }
	  }
	}
	
      if(!treated){
	cout <<"trying to load component c="<<c<<endl;
	cout <<"There was an untreated point (i,j)=("<<i<<","<<j<<") NP="<< NP<<endl;
	cout <<__FILE__<<"   "<<__LINE__<<endl;
	cout <<"Data..."<<endl;
	terminate_program(2);
      }
      if(NP<1.0){
	data(c,i,j) = 0.0;//uninitialized;//UNINIECO;//This is correct for ecoclimap and gtopo
      }
      nPoints[(i-1)+(j-1)*Klon] = NP;
    }
}



void ModelData::loadData(int c,LakeData *globaldata, double uninitialized){
  double eps = 0.005;
  int one=1;

  int Klon = domain->klon;
  int Klat = domain->klat;

  double polon = domain->polon, polat = domain->polat;
  double lonPole = 0.0;
  double nPole = 90.0;
  double rotatedNPole[2],rotatedSPole[2];
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  nPole *= -1.0;
  reg2rot_(&lonPole,&nPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);

  
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      double vertx[4], verty[4];
      double lonxRot[4],latyRot[4];
      
      domain->getCell(i,j,vertx,verty);//cell around model->(i,j) in regular space
      domain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space
      int isOnPole=0; 
      if(domain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot))
	isOnPole = 1;
      
      if(domain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot))
	isOnPole = -1;
      
            
      bool crossesDateLine = false;
      int n,m;
      for(n=0,m=3;n<4;m=n++)
	if(vertx[n]*vertx[m]<0.0)
	  if((fabs(vertx[n]) + fabs(vertx[m])) > 350.0)
	    crossesDateLine = true;
      
      //Stuff in regular space
      double minLonModel = MIN(lon(i,j),vertx[0]);
      double maxLonModel = MAX(lon(i,j),vertx[0]);
      
      double minLatModel = MIN(lat(i,j),verty[0]);
      double maxLatModel = MAX(lat(i,j),verty[0]);
      for(int l=1;l<4;l++){
	minLonModel = MIN(vertx[l],minLonModel);
	maxLonModel = MAX(vertx[l],maxLonModel);
	minLatModel = MIN(verty[l],minLatModel);
	maxLatModel = MAX(verty[l],maxLatModel);
      }
	
      double iisD = globaldata->lonindex(minLonModel,minLatModel);      
      double iieD = globaldata->lonindex(maxLonModel,maxLatModel);
	
      double jjsD = globaldata->latindex(minLonModel,minLatModel);
      double jjeD = globaldata->latindex(maxLonModel,maxLatModel);
	
	
      //this is if the index2data is inverse
      if(iisD>iieD){
	double tmp=iisD;
	iisD = iieD;
	iieD = tmp;
      }
      if(jjsD>jjeD){
	double tmp=jjsD;
	jjsD = jjeD;
	jjeD = tmp;
      }

      int iis = MIN(globaldata->klon(),MAX(1,floor(iisD)));
      int iie = MIN(globaldata->klon(),MAX(1,ceil(iieD)));

      int jjs = MIN(globaldata->klat(),MAX(1,floor(jjsD)));
      int jje = MIN(globaldata->klat(),MAX(1,ceil(jjeD)));
	
      if(isOnPole!=0){
	double safe = 0.6*dlon()/globaldata->dlon()+2;
	iisD -= safe;
	jjsD -= safe;
	iieD += safe;
	jjeD += safe;
	if(iisD<1 || iieD>globaldata->klon()){
	  iis = 1;
	  iie = globaldata->klon();
	}
	if(jjsD<1 || jjeD>globaldata->klat()){
	  jjs = 1;
	  jje = globaldata->klat();
	}
      }
	
      if(crossesDateLine){
	iis = 1;
	iie = globaldata->klon();
      }

      bool treated=false;
      double NP[3]={0.0,0.0,0.0};
      double Npo= 0.0;
    
      for(int jj=jjs;jj<=jje;jj++){
	for(int ii=iis;ii<=iie;ii++){
	  bool inside;
	  double xreg = globaldata->lon(ii,jj);
	  double yreg = globaldata->lat(ii,jj);
	  
	  if(isOnPole!=0 || crossesDateLine){
	    double xrot,yrot;
	    reg2rot_(&xreg,&yreg,&xrot,&yrot,&one,&one,&polon,&polat);
	    inside = domain->isInside(xrot,yrot,lonxRot,latyRot);
	  }
	  else
	    inside = domain->isInside(xreg,yreg,vertx,verty);
	  
	  if(inside){
	    Npo += 1.0; // total number of points
	    double dat = globaldata->data(ii,jj);
	    treated = true;
	    if(fabs(dat-uninitialized)>eps){ 
	      double deep=8.0,medium=3.0,shallow=0.0;
	      int ptr = dat>deep ? 0 : (dat>medium ? 1 : 2);
	      NP[ptr] += 1.0;

	      if(useVar){//variance
		if(NP[ptr]>1)
		  var(c+ptr,i,j) = (NP[ptr]-2.0)/(NP[ptr]-1.0)*var(c+ptr,i,j)+ 
		    1.0/(NP[ptr])*pow((data(c+ptr,i,j)-dat),2);
		else
		  var(c+ptr,i,j) = 0.0;	  
	      }
	      data(c+ptr+3,i,j) = ((NP[ptr]-1.0)*data(c+ptr+3,i,j)+dat)/NP[ptr]; //mean (upscaling)
      
	    }
	    for(int cc=0;cc<3;cc++){
	      data(c+cc,i,j) = NP[cc]/Npo; //fraction deep,medium,shallow
	      if(useVar)
		var(c+cc,i,j) = 0.0;
	    }
	  }
	}
      }
	
      if(!treated){
	//globaldata->domain->print("Not treted lake",111);
	if(isOnPole!=0 || crossesDateLine)cout <<"isOnpole or crosses date line"<<endl;
	else cout <<"ordinary point"<<endl;
	cout <<"trying to load component c="<<c<<endl;
	cout <<"There was an untreated point (i,j)=("<<i<<","<<j<<") NP="<< NP<<endl;
	cout <<"lon = "<<lon(i,j)<<" lat="<<lat(i,j)<<endl;
	cout <<"jj="<<jjs<<" , "<<jje <<"jjsD="<<jjsD<<" jjeD="<<jjeD<<endl;
	cout <<"ii="<<iis<<" , "<<iie <<endl;
	cout <<"glob-klat="<<globaldata->klat()<<"glob-klon="<<globaldata->klon()<<endl;
	cout <<"lakedata....";
	cout <<"polon="<< polon<<"  polat= "<<polat<<endl;
	cout <<"Npo = "<<Npo<<endl;
	cout <<__FILE__<<"   "<<__LINE__<<endl;
	terminate_program(2);
      }
      for(int t=0;t<3;t++)
      if(NP[t]<1.0){
	data(c+t,i,j) = 0.0;
	data(c+3+t,i,j) =0.0;
      }

      if(Npo<1.0){
	for(int cc=0;cc<6;cc++)
	  data(c+cc,i,j) = 0.0;//uninitialized;//UNINIECO;//This is correct for ecoclimap and gtopo
      }
      nPoints[(i-1)+(j-1)*Klon] = Npo;
    }
}


void ModelData::loadData_z0(int c,EcoclimapData *ecodata,double UNINIECO){

  double eps = 0.005;
  int one=1;

  int Klon = domain->klon;
  int Klat = domain->klat;
  double polon = domain->polon, polat = domain->polat;
  double lonPole = 0.0;
  double nPole = 90.0;
  double rotatedNPole[2],rotatedSPole[2];
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  nPole *= -1.0;
  reg2rot_(&lonPole,&nPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);


  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      double vertx[4], verty[4];
      double lonxRot[4],latyRot[4];
      domain->getCell(i,j,vertx,verty);//cell around model->(i,j) in regular space
      domain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space

      int isOnPole=0; 
      if(domain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot))
	isOnPole = 1;
      
      if(domain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot))
	isOnPole = -1;
      
            
      bool crossesDateLine = false;
      int n,m;
      for(n=0,m=3;n<4;m=n++)
	if(vertx[n]*vertx[m]<0.0)
	  if((fabs(vertx[n]) + fabs(vertx[m])) > 350.0)
	    crossesDateLine = true;

      double minLonModel = MIN(lon(i,j),vertx[0]);
      double maxLonModel = MAX(lon(i,j),vertx[0]);
      
      double minLatModel = MIN(lat(i,j),verty[0]);
      double maxLatModel = MAX(lat(i,j),verty[0]);
      for(int l=1;l<4;l++){
	minLonModel = MIN(vertx[l],minLonModel);
	maxLonModel = MAX(vertx[l],maxLonModel);
	minLatModel = MIN(verty[l],minLatModel);
	maxLatModel = MAX(verty[l],maxLatModel);
      }
      double iisD = ecodata->lonindex(minLonModel,minLatModel);      
      double iieD = ecodata->lonindex(maxLonModel,maxLatModel);

      double jjsD = ecodata->latindex(minLonModel,minLatModel);
      double jjeD = ecodata->latindex(maxLonModel,maxLatModel);


      //this is if the index2data is inverse
      if(iisD>iieD){
 	double tmp=iisD;
 	iisD = iieD;
 	iieD = tmp;
      }
      if(jjsD>jjeD){
	double tmp=jjsD;
	jjsD = jjeD;
	jjeD = tmp;
      }

      int iis = MIN(ecodata->klon(),MAX(1,floor(iisD)));
      int jjs = MIN(ecodata->klat(),MAX(1,floor(jjsD)));
      int iie = MIN(ecodata->klon(),MAX(1,ceil(iieD)));
      int jje = MIN(ecodata->klat(),MAX(1,ceil(jjeD)));

      if(isOnPole!=0){
	double safe = 0.6*dlon()/ecodata->dlon()+2;
	iisD -= safe;
	jjsD -= safe;
	iieD += safe;
	jjeD += safe;
	if(iisD<1 || iieD>ecodata->klon()){
	  iis = 1;
	  iie = ecodata->klon();
	}
	if(jjsD<1 || jjeD>ecodata->klat()){
	  jjs = 1;
	  jje = ecodata->klat();
	}
      }
	
      if(crossesDateLine){
	iis = 1;
	iie = ecodata->klon();
      }

      bool treated=false;
      double NP=0.0;
      double zlb=40;
      double theta=0.0;
      double logz0sum=0.0;
      double sum2 = 0.0;

      
      for(int jj=jjs;jj<=jje;jj++)
	for(int ii=iis;ii<=iie;ii++){
	  bool inside;
	  double xreg = ecodata->lon(ii,jj);
	  double yreg = ecodata->lat(ii,jj);
	  if(isOnPole!=0 || crossesDateLine){
	    double xrot,yrot;
	    reg2rot_(&xreg,&yreg,&xrot,&yrot,&one,&one,&polon,&polat);
	    inside = domain->isInside(xrot,yrot,lonxRot,latyRot);
	  }
	  else
	    inside = domain->isInside(xreg,yreg,vertx,verty);

	  //	  if(domain->isInside(xreg,yreg,vertx,verty)){//regular space
	  if(inside){
	    double z0 = ecodata->data(ii,jj);
	    treated = true;
	    if(fabs(z0-UNINIECO)>eps){ //dat==-999
	      NP += 1.0;
	      theta = log(zlb/z0);
	      logz0sum += log(z0)/theta;
	      sum2 += 1.0/theta;
	    }
	  }
	}
      if(sum2>0)
	data(c,i,j) = exp(logz0sum/sum2); //z0 upscaling   
      
      if(!treated){
	cout <<"trying to load component c="<<c<<endl;
	cout <<"There was an untreated point (i,j)=("<<i<<","<<j<<") NP="<< NP<<endl;
	cout <<__FILE__<<"   "<<__LINE__<<endl;
	terminate_program(2);
      }
      if(NP<1.0)data(c,i,j) = 0.0;//UNINIECO;
      nPoints[(i-1)+(j-1)*Klon] = NP;
    }
}



void ModelData::loadData(int c,EcoclimapData *globaldata,double uninitialized,EcoclimapData *F_wat){
  double eps = 0.005;
  int one=1;

  int Klon = domain->klon;
  int Klat = domain->klat;

  double polon = domain->polon, polat = domain->polat;
  double lonPole = 0.0;
  double nPole = 90.0;
  double rotatedNPole[2],rotatedSPole[2];
  reg2rot_(&lonPole,&nPole,&(rotatedNPole[0]),&(rotatedNPole[1]),&one,&one,&polon,&polat);
  nPole *= -1.0;
  reg2rot_(&lonPole,&nPole,&(rotatedSPole[0]),&(rotatedSPole[1]),&one,&one,&polon,&polat);

  
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      double vertx[4], verty[4];
      double lonxRot[4],latyRot[4];
      
      domain->getCell(i,j,vertx,verty);//cell around model->(i,j) in regular space
      domain->getCellR(i,j,lonxRot,latyRot);//cell in rotated space
      int isOnPole=0; 
      if(domain->isInside(rotatedNPole[0],rotatedNPole[1],lonxRot,latyRot))
	isOnPole = 1;
      
      if(domain->isInside(rotatedSPole[0],rotatedSPole[1],lonxRot,latyRot))
	isOnPole = -1;
      
            
      bool crossesDateLine = false;
      int n,m;
      for(n=0,m=3;n<4;m=n++)
	if(vertx[n]*vertx[m]<0.0)
	  if((fabs(vertx[n]) + fabs(vertx[m])) > 350.0)
	    crossesDateLine = true;
      
      //Stuff in regular space
      double minLonModel = MIN(lon(i,j),vertx[0]);
      double maxLonModel = MAX(lon(i,j),vertx[0]);
      
      double minLatModel = MIN(lat(i,j),verty[0]);
      double maxLatModel = MAX(lat(i,j),verty[0]);
      for(int l=1;l<4;l++){
	minLonModel = MIN(vertx[l],minLonModel);
	maxLonModel = MAX(vertx[l],maxLonModel);
	minLatModel = MIN(verty[l],minLatModel);
	maxLatModel = MAX(verty[l],maxLatModel);
      }
	
      double iisD = globaldata->lonindex(minLonModel,minLatModel);      
      double iieD = globaldata->lonindex(maxLonModel,maxLatModel);
	
      double jjsD = globaldata->latindex(minLonModel,minLatModel);
      double jjeD = globaldata->latindex(maxLonModel,maxLatModel);
	
	
      //this is if the index2data is inverse
      if(iisD>iieD){
	double tmp=iisD;
	iisD = iieD;
	iieD = tmp;
      }
      if(jjsD>jjeD){
	double tmp=jjsD;
	jjsD = jjeD;
	jjeD = tmp;
      }

      int iis = MIN(globaldata->klon(),MAX(1,floor(iisD)));
      int iie = MIN(globaldata->klon(),MAX(1,ceil(iieD)));

      int jjs = MIN(globaldata->klat(),MAX(1,floor(jjsD)));
      int jje = MIN(globaldata->klat(),MAX(1,ceil(jjeD)));
	
      if(isOnPole!=0){
	double safe = 0.6*dlon()/globaldata->dlon()+2;
	iisD -= safe;
	jjsD -= safe;
	iieD += safe;
	jjeD += safe;
	if(iisD<1 || iieD>globaldata->klon()){
	  iis = 1;
	  iie = globaldata->klon();
	}
	if(jjsD<1 || jjeD>globaldata->klat()){
	  jjs = 1;
	  jje = globaldata->klat();
	}
      }
	
      if(crossesDateLine){
	iis = 1;
	iie = globaldata->klon();
      }

      bool treated=false;
      double NP=0.0;

    
      for(int jj=jjs;jj<=jje;jj++)
	for(int ii=iis;ii<=iie;ii++){
	  bool inside;
	  double xreg = globaldata->lon(ii,jj);
	  double yreg = globaldata->lat(ii,jj);
	  
	  if(isOnPole!=0 || crossesDateLine){
	    double xrot,yrot;
	    reg2rot_(&xreg,&yreg,&xrot,&yrot,&one,&one,&polon,&polat);
	    inside = domain->isInside(xrot,yrot,lonxRot,latyRot);
	  }
	  else
	    inside = domain->isInside(xreg,yreg,vertx,verty);
	  
	  if(inside){
	    double dat = globaldata->data(ii,jj);
	    treated = true;
	    if(fabs(dat-uninitialized)>eps && F_wat->data(ii,jj)<1.0 ){ 
	      NP += 1.0;
	      if(useVar){//variance
		if(NP>1)
		  var(c,i,j) = (NP-2.0)/(NP-1.0)*var(c,i,j) + 
		    1.0/(NP)*pow((data(c,i,j)-dat),2);
		else
		  var(c,i,j) = 0.0;	  
	      }
	      data(c,i,j) = ((NP-1.0)*data(c,i,j)+dat)/NP; //mean (upscaling)
	    }
	  }
	}

      

	
      if(!treated){
	cout <<"trying to load component c="<<c<<endl;
	cout <<"There was an untreated point (i,j)=("<<i<<","<<j<<") NP="<< NP<<endl;
	cout <<__FILE__<<"   "<<__LINE__<<endl;
	terminate_program(2);
      }
      if(NP<1.0){
	data(c,i,j) = 0.0;//uninitialized;//UNINIECO;//This is correct for ecoclimap and gtopo
      }
      nPoints[(i-1)+(j-1)*Klon] = NP;
    }
}




ModelData::ModelData(Domain *domain, int nc_in, bool use_variance):
  Data(domain,2*nc_in){

  useVar = true;
  nc /= 2;
  int Klon=domain->klon;
  int Klat=domain->klat;
  nPoints = new double[Klon*Klat];

  for(int j=1;j<=Klat;j++){
    for(int i=1;i<=Klon;i++){
      nPoints[(i-1)+(j-1)*Klon]=0.0;
      for(int c=1;c<=nc;c++){
	data(c,i,j)=UNINI;
	var(c,i,j)=UNINI;
      }
    }
  }
}


ModelData::ModelData(Domain *domain,  int nc_in):
  Data(domain,nc_in)
{
  useVar = false;
  int Klon=domain->klon;
  int Klat=domain->klat;
  nPoints = new double[Klon*Klat];
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++){
      nPoints[(i-1)+(j-1)*Klon]=0.0;
    }
}


void ModelData::check(){
  //   for(int j=1;j<=Klat;j++)
  //     for(int i=1;i<=Klon;i++)
  //       if(nPoints[(i-1)+(j-1)*Klon]<1.0){
  // // 	cout << "point = ("<<loncoords[(i-1)+(j-1)*Klon]<<","<<latcoords[(i-1)+(j-1)*Klon]<<") = "<<data(i,j)<<endl;
  // 	cout <<" is surrounded by cell"<<endl;
  // 	double vx[4],vy[4];
  // 	getCell( i, j,vx, vy);
  // 	for(int l=0;l<4;l++)
  // 	  correctCoordinates(vx[l],vy[l]);
  // 	double minLonModel = MIN(lon(i,j),vertx[0]);
  // 	double maxLonModel = MAX(lon(i,j),vertx[0]);
	
  // 	double minLatModel = MIN(lat(i,j),verty[0]);
  // 	double maxLatModel = MAX(lat(i,j),verty[0]);
  // 	for(int l=1;l<4;l++){
  // 	  minLonModel = MIN(vertx[l],minLonModel);
  // 	  maxLonModel = MAX(vertx[l],maxLonModel);
	  
  // 	  minLatModel = MIN(verty[l],minLatModel);
  // 	  maxLatModel = MAX(verty[l],maxLatModel);
  // 	}
	
  // 	int iis = MIN(g30data->klon(),MAX(1,floor((minLonModel-west30)/dlon30))); 
  // 	int jjs = MIN(g30data->klat(),MAX(1,floor((minLatModel-south30)/dlat30))); 
  // 	int iie = MAX(1,MIN(g30data->klon(),ceil( (maxLonModel-west30)/dlon30  )));
  // 	int jje = MAX(1,MIN(g30data->klat(),ceil( (maxLatModel-south30)/dlat30 )));
	

  // 	cout <<"lonx=[";
  // 	for(int k=0;k<4;k++)
  // 	  cout <<vx[k]<<" ";
  // 	cout<<vx[0] <<"];"<<endl;
  // 	cout <<"laty=[";
  // 	for(int k=0;k<4;k++)
  // 	  cout <<vy[k]<<" ";
  // 	cout<<vy[0] <<"];"<<endl;
  //       }
}


ModelData::~ModelData(){
  delete[] nPoints;
  delete[] datA;
  delete domain;
} 

void ModelData::printDomain(int stride){
  Data::printDomain(stride);
  int Klon=domain->klon;
  int Klat=domain->klat;

  if(useVar)
    for(int c=1;c<=nc;c++){
      cout << "var=[";
      for(int j=1;j<=Klat;j+=stride){
	for(int i=1;i<=Klon;i+=stride)
	  cout<<var(c,i,j)<<" ";
	cout << endl;
      }
      cout << "];"<<endl;
    }
}

// void ModelData::interpolate(int c,ModelData *m1,int* compo, double* w,double UNINIECO){
//   double eps=0.005;
//   int Klon=domain->klon;
//   int Klat=domain->klat;

//   for(int j=1;j<=Klat;j++)
//     for(int i=1;i<=Klon;i++){
//       double d0 = m1->data(compo[0],i,j);
//       double d1 = m1->data(compo[0],i,j);
//       double v0=0,v1=0;
//       if(useVar){
// 	v0 = m1->var(compo[0],i,j);
// 	v1 = m1->var(compo[1],i,j);
//       }
//       if(fabs(d0-UNINIECO)<eps && fabs(d1-UNINIECO)<eps){
// 	data(c,i,j) = UNINIECO;
// 	if(useVar)
// 	  var(c,i,j) = UNINIECO;
//       }
//       else if(fabs(d0-UNINIECO)<eps){
// 	data(c,i,j) = d1;
// 	if(useVar)
// 	  var(c,i,j) = v1;
//       }
//       else if(fabs(d1-UNINIECO)<eps){
// 	data(c,i,j) = d0;
// 	if(useVar)
// 	  var(c,i,j) = v0;
//       }
//       else{
// 	data(c,i,j) = w[0]*d0 + w[1]*d1;
// 	if(useVar)
// 	  var(c,i,j) = w[0]*v0 + w[1]*v1;
//       }
      
//     }
  
// }


void ModelData::loadMinData(int c,ModelData *an,int cstart,int cend,double UNINIECO){

  int Klon=domain->klon;
  int Klat=domain->klat;

  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++)
      data(c,i,j) = UNINIECO;
  
  for(int k=cstart;k<=cend;k++)
    for(int j=1;j<=Klat;j++)
      for(int i=1;i<=Klon;i++)
	if(data(c,i,j)> an->data(k,i,j))
	  data(c,i,j) =  an->data(k,i,j);
}

void ModelData::loadMaxData(int c,ModelData *an,int cstart,int cend,double UNINIECO){
  int Klon=domain->klon;
  int Klat=domain->klat;
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++)
      data(c,i,j) = -UNINIECO;

  for(int k=cstart;k<=cend;k++)
    for(int j=1;j<=Klat;j++)
      for(int i=1;i<=Klon;i++)
	if(data(c,i,j) < an->data(k,i,j))
	  data(c,i,j) =  an->data(k,i,j);
}


double ModelData::echoMax(int c)
{

  int Klon=domain->klon;
  int Klat=domain->klat;
  double mes = -9999;
  for(int j=1;j<=Klat;j++)
    for(int i=1;i<=Klon;i++)
      if(data(c,i,j) > mes)
	mes = data(c,i,j);
  return mes;
  
}


bool ModelData::readFile(int c,char* filename,MPI_Comm localComm){

  int Klon=domain->klon;
  int Klat=domain->klat;

  int len = strlen(filename);
  char *f=new char[len+11];
  strcpy(f,"model_");
  strcat(f,filename);
  strcat(f,".bin");

  char *fv = new char[len+15];
  strcpy(fv,"model_");
  strcat(fv,filename);
  strcat(fv,"_var.bin");

  //initialize and decompose domain in readFile?
  bool isRead = Data::readFile(f,localComm,c);
  
  if(useVar){
    isRead =Data::readFile(fv,localComm,nc+c);
  }

  
  if(isRead){
    for(int j=1;j<=Klat;j++){
      for(int i=1;i<=Klon;i++){
	nPoints[(i-1)+(j-1)*Klon]=0.0;
      }
    }
  }

  delete[] f;
  delete[] fv;

  return isRead;
}


void ModelData::writeFile(int c,char* filename,MPI_Comm localComm){

  int len = strlen(filename);

  char *f = new char[len+11];
  strcpy(f,"model_");
  strcat(f,filename);
  strcat(f,".bin");

  char *fv = new char[len+15];
  strcpy(fv,"model_");
  strcat(fv,filename);
  strcat(fv,"_var.bin");

  Data::writeFile(f,localComm,c);
  
  if(useVar){
    Data::writeFile(fv,localComm,nc+c);
  }

  delete[] f;
  delete[] fv;
}

