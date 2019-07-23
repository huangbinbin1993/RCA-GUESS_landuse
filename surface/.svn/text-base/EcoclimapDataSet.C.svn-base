#ifdef MPI_SRC
#include"mpi.h"
#endif
#include"EcoclimapDataSet.h"
#include"LakeData.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<stdlib.h>

using namespace std;

void EcoclimapDataSet::setTimeDep(int i,bool dep){
  timeDep[i-1] = dep; 
}

EcoclimapDataSet::~EcoclimapDataSet(){
#ifdef DEBUG
  cout <<"~EcoclimapDataSet() is called"<<endl;
#endif
  for(int i=0;i<12;i++)
    delete[] months[i];
  delete[] months;
  for(int i=0;i<3;i++)
    delete[] tiles[i];
  delete[] tiles;
  for(int i=0;i<numberOfVariables;i++){
    delete[] filenames[i];
    delete[] pureFilenames[i];
  }
  delete[] filenames;
  delete[] pureFilenames;
  delete[] timeDep;
  delete[] readFromfile;
  delete ecoclimapData;
  delete lakeData;
}

EcoclimapDataSet::EcoclimapDataSet(char* path_in,Domain *eco,int N, int M):
  path(path_in){
  int nMonths = 12;
  months = new char*[nMonths];
  
  for(int i=0;i<nMonths;i++){
    months[i] = new char[4];
    char tmp[sizeof(int)*8+1];  
    sprintf( tmp, "%d", i+1 );
    strcpy(months[i],".");
    if(i+1<10)
      strcat(months[i],"0");
    strcat(months[i],tmp);// e.g. ".01" ".11"
  }

  int nTiles = 3; 
  tiles = new char*[nTiles]; //{".t01",".t02",".t03"};
  for(int i=0;i<nTiles;i++){
    tiles[i] = new char[5];
    char tmp[sizeof(int)*8+1];  
    sprintf( tmp, "%d", i+1 );
    strcpy(tiles[i],".t0");
    strcat(tiles[i],tmp);// e.g. ".t01"
  }
  
  //lai_t{1,2,3}_{1-12},
  //z0_t{1,2,3}_{1-12},
  //emis_t{1,2,3}_{1-12}
  //alb_t{1,2,3}_{1-12}
  //veg_t{1,2,3}_{1-12}
  numberOfVariables = 5*nMonths*nTiles;
  
  numberOfVariables += 5*nTiles; //droot_t{1,2,3},
                              //dsoil_t{1,2,3},
                              //rsmin_t{1,2,3}
                              //alb_veg_t{1,2,3}
                          //frac_t{1,2,3},  

  numberOfVariables += 6;//f_wat,
                         //alb_soil,
                         //clay,
                         //sand,
                         //frac_lake,
                         //soil_carb
  numberOfVariables += 2*nTiles +1;

  numberOfVariables++;//this is the lake database
  
  filenames = new char*[numberOfVariables];
  pureFilenames = new char*[numberOfVariables];

  timeDep = new bool[numberOfVariables];
  readFromfile = new bool[numberOfVariables];

  connectFilenames();

  ecoclimapData = new EcoclimapData(N,M,eco);
  Domain *lakeDomain = new Domain(eco);
  lakeData = new LakeData(N,M,lakeDomain);
}

int EcoclimapDataSet::getNumberOfTimeDepVars(){
  int n=0;
  for(int i=0;i<numberOfVariables;i++)
    if(timeDep[i])n++;
  return n/12;
}

bool EcoclimapDataSet::isTimeDep(int i){
  return timeDep[i-1];
}

bool EcoclimapDataSet::readFromFile(int i){
  return readFromfile[i-1];
}


void EcoclimapDataSet::connectFilenames(){
  char* name;
  int len;
  int variable = 0;

  for(int i=0;i<numberOfVariables;i++){
    timeDep[i] = false;
    readFromfile[i]=true;
  }
  int lenTile =  4;
  int lenMonth =  3;
  int lEnd = 4;
  
  char tmp[sizeof(int)*8+1];  
  int resolution =1;
  sprintf( tmp, "%d", resolution );
  char res[5];
  strcpy(res,".");
  if(resolution<10)
    strcat(res,"00");
  else if(resolution<100)
    strcat(res,"0");
  strcat(res,tmp);// e.g. ".001"


  int lenPath = static_cast<int>(strlen(path));
   
  len = 3+lenTile+lenMonth+lEnd;
  //strcpy(name,"lai");//.month.t{tile}.001
  for(int t=0;t<3;t++){
    for(int m=0;m<12;m++){
      timeDep[variable]=true;
      filenames[variable] = new char[lenPath+len+2];
      strcpy(filenames[variable],path);
      strcat(filenames[variable],"/");

      strcat(filenames[variable],"lai");

      strcat(filenames[variable],months[m]);

      strcat(filenames[variable],tiles[t]);
      strcat(filenames[variable],res);
      
      pureFilenames[variable] = new char[3+lenTile+lenMonth+1];
      strcpy(pureFilenames[variable],"lai");
      strcat(pureFilenames[variable],months[m]);
      strcat(pureFilenames[variable],tiles[t]);
      variable++;
    }
  }
  
  len = 2+lenTile+lenMonth+lEnd;
  for(int t=0;t<3;t++){
    for(int m=0;m<12;m++){
      timeDep[variable]=true;
      filenames[variable] = new char[lenPath+len+2];
      strcpy(filenames[variable],path);
      strcat(filenames[variable],"/");
      strcat(filenames[variable],"z0");//.month.t{tile}.001
      strcat(filenames[variable],months[m]);
      strcat(filenames[variable],tiles[t]);
      strcat(filenames[variable],res);
      
      pureFilenames[variable] = new char[2+lenTile+lenMonth+1];
      strcpy(pureFilenames[variable],"z0");
      strcat(pureFilenames[variable],months[m]);
      strcat(pureFilenames[variable],tiles[t]);
      variable++;
    }
  }

  len = 4+lenTile+lenMonth+lEnd;
  for(int t=0;t<3;t++){
    for(int m=0;m<12;m++){
      timeDep[variable]=true;
      filenames[variable] = new char[lenPath+len+2];
      strcpy(filenames[variable],path);
      strcat(filenames[variable],"/");
      strcat(filenames[variable],"emis");//.month.t{tile}.001
      strcat(filenames[variable],months[m]);
      strcat(filenames[variable],tiles[t]);

      strcat(filenames[variable],res);
      
      pureFilenames[variable] = new char[4+lenTile+lenMonth+1];
      strcpy(pureFilenames[variable],"emis");
      strcat(pureFilenames[variable],months[m]);
      strcat(pureFilenames[variable],tiles[t]);
      variable++;
    }
  }

  len = 3+lenTile+lenMonth+lEnd;
  for(int t=0;t<3;t++){
    for(int m=0;m<12;m++){
      timeDep[variable]=true;
      filenames[variable] = new char[lenPath+len+2];
      strcpy(filenames[variable],path);
      strcat(filenames[variable],"/");
      strcat(filenames[variable],"alb");//.month.t{tile}.001
      strcat(filenames[variable],months[m]);
      strcat(filenames[variable],tiles[t]);

      strcat(filenames[variable],res);

      
      pureFilenames[variable] = new char[3+lenTile+lenMonth+1];
      strcpy(pureFilenames[variable],"alb");
      strcat(pureFilenames[variable],months[m]);
      strcat(pureFilenames[variable],tiles[t]);
      variable++;
    }
  }

  len = 3+lenTile+lenMonth+lEnd;
  for(int t=0;t<3;t++){
    for(int m=0;m<12;m++){
      timeDep[variable]=true;
      filenames[variable] = new char[lenPath+len+2];
      strcpy(filenames[variable],path);
      strcat(filenames[variable],"/");
      strcat(filenames[variable],"veg");//.month.t{tile}.001
      strcat(filenames[variable],months[m]);
      strcat(filenames[variable],tiles[t]);

      strcat(filenames[variable],res);

      
      pureFilenames[variable] = new char[3+lenTile+lenMonth+1];
      strcpy(pureFilenames[variable],"veg");
      strcat(pureFilenames[variable],months[m]);
      strcat(pureFilenames[variable],tiles[t]);
      variable++;
    }
  }

  len = 5+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"F_wat");//.001

  strcat(filenames[variable],res);


  pureFilenames[variable] = new char[5+1];
  strcpy(pureFilenames[variable],"F_wat");
  variable++;

  len = 11+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"albedo_soil");//.001

  strcat(filenames[variable],res);

  
  pureFilenames[variable] = new char[11+1];
  strcpy(pureFilenames[variable],"albedo_soil");
  variable++;

  len = 4+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"clay");//.001

  strcat(filenames[variable],res);


  pureFilenames[variable] = new char[4+1];
  strcpy(pureFilenames[variable],"clay");
  variable++;

  len = 4+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"sand");//.001

  strcat(filenames[variable],res);

  
  pureFilenames[variable] = new char[4+1];
  strcpy(pureFilenames[variable],"sand");
  variable++;

  len = 11+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"frac_tile01");//.001

  strcat(filenames[variable],res);

	
  pureFilenames[variable] = new char[11+1];
  strcpy(pureFilenames[variable],"frac_tile01");
  variable++;

  len = 11+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"frac_tile02");//.001

  strcat(filenames[variable],res);

  
  pureFilenames[variable] = new char[11+1];
  strcpy(pureFilenames[variable],"frac_tile02");
  variable++;

  len = 11+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"frac_tile03");//.001

  strcat(filenames[variable],res);

  
  pureFilenames[variable] = new char[11+1];
  strcpy(pureFilenames[variable],"frac_tile03");
  variable++;

  len = 6+lenTile+lEnd;
  for(int t=0;t<3;t++){
    filenames[variable] = new char[lenPath+len+2];
    strcpy(filenames[variable],path);
    strcat(filenames[variable],"/");
    strcat(filenames[variable],"d_root");//.t{tile}.001
    strcat(filenames[variable],tiles[t]);

    strcat(filenames[variable],res);


    
    pureFilenames[variable] = new char[6+lenTile+1];
    strcpy(pureFilenames[variable],"d_root");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  len = 6+lenTile+lEnd;
  for(int t=0;t<3;t++){
    filenames[variable] = new char[lenPath+len+2];
    strcpy(filenames[variable],path);
    strcat(filenames[variable],"/");
    strcat(filenames[variable],"d_soil");//.t{tile}.001
    strcat(filenames[variable],tiles[t]);

    strcat(filenames[variable],res);

    
    pureFilenames[variable] = new char[6+lenTile+1];
    strcpy(pureFilenames[variable],"d_soil");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  len = 5+lenTile+lEnd;
  for(int t=0;t<3;t++){
    filenames[variable] = new char[lenPath+len+2];
    strcpy(filenames[variable],path);
    strcat(filenames[variable],"/");
    strcat(filenames[variable],"rsmin");//.t{tile}.001
    strcat(filenames[variable],tiles[t]);

    strcat(filenames[variable],res);

    
    pureFilenames[variable] = new char[5+lenTile+1];
    strcpy(pureFilenames[variable],"rsmin");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  len = 10+lenTile+lEnd;
  for(int t=0;t<3;t++){
    filenames[variable] = new char[lenPath+len+2];
    strcpy(filenames[variable],path);
    strcat(filenames[variable],"/");
    strcat(filenames[variable],"albedo_veg");//.t{tile}.001
    strcat(filenames[variable],tiles[t]);

    strcat(filenames[variable],res);

    
    pureFilenames[variable] = new char[10+lenTile+1];
    strcpy(pureFilenames[variable],"albedo_veg");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  len = 6+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"F_lake");//.001

  strcat(filenames[variable],res);


  
  pureFilenames[variable] = new char[6+1];
  strcpy(pureFilenames[variable],"F_lake");
  variable++;

  len = 9+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"soil_carb");//.001

  strcat(filenames[variable],res);

  
  pureFilenames[variable] = new char[9+1];
  strcpy(pureFilenames[variable],"soil_carb");
  variable++;
    
  //Below are the derived qtns.
  len =7+lEnd;
  readFromfile[variable] = false;
  filenames[variable] = new char[2];
  strcpy(filenames[variable],"1");
  pureFilenames[variable] = new char[7+1];
  strcpy(pureFilenames[variable],"texture");
  variable++;

  len =6+lEnd+lenTile;
  for(int t=0;t<3;t++){
    readFromfile[variable] = false;
    filenames[variable] = new char[2];
    strcpy(filenames[variable],"2");
    pureFilenames[variable] = new char[6+lenTile+1];
    strcpy(pureFilenames[variable],"minlai");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  len =6+lEnd+lenTile;
  for(int t=0;t<3;t++){
    readFromfile[variable] = false;
    filenames[variable] = new char[2];
    strcpy(filenames[variable],"3");
    pureFilenames[variable] = new char[6+lenTile+1];
    strcpy(pureFilenames[variable],"maxlai");
    strcat(pureFilenames[variable],tiles[t]);
    variable++;
  }

  //lakedatabase GlobalLakeDepth.001
  //if(variable!=210)terminate_program(9);
  //cout << "This is supposed to be lakedatabase"<<variable <<endl;
  //variable = 209;
  readFromfile[variable] = true;
  len = 15+lEnd;
  filenames[variable] = new char[lenPath+len+2];
  strcpy(filenames[variable],path);
  strcat(filenames[variable],"/");
  strcat(filenames[variable],"GlobalLakeDepth");//.001
  strcat(filenames[variable],res);
  pureFilenames[variable] = new char[len+1+lEnd];
  strcpy(pureFilenames[variable],"GlobalLakeDepth");
  strcat(pureFilenames[variable],res);
  variable++;
}

int EcoclimapDataSet::getVariableNr(int nn,int month){
  
  switch(nn){
  case 1: //lai_t1
    return month; //1...12
    break;
  case 2://lai_t2
    return 12+month; //13...24
    break;
  case 3://lai_t3
    return 2*12+month; //25...36
    break;
  case 4://z0_t1
    return 3*12+month; //37...48
    break;
  case 5://z0_t2
    return 4*12+month; //49...60
    break;
  case 6://z0_t3
    return 5*12+month;//61...72
    break;
  case 7://emis_t1
    return 6*12+month;//73...84
    break;
  case 8://alb_t1
    return 9*12+month;//109...120
    break;
  case 9://veg_t1
    return 12*12+month;//145...156
    break;
  case 10://frland
    return 15*12+1;//181
    break;
  case 11://alb_soil
    return 15*12+2;//182
    break;
  case 12://clay
    return 15*12+3; //183
    break;
  case 13://sand
    return 15*12+4; //184
    break;
  case 14://frac_t1
    return 15*12+5; //185
    break;
  case 15://frac_t2
    return 15*12+6; //185
    break;
  case 16://frac_t3
    return 15*12+7; //187
    break;
  case 17://alb_t2
    return 10*12+month;//121...132
    break;
  case 18://alb_t3
    return 11*12+month; //133...144
    break;
  case 19://emis_t2
    return 7*12+month;//85...96
    break;
  case 20://emis_t3
    return 8*12+month;//97...108
    break;
  case 21://veg_t2
    return 13*12+month;//157...168
    break;
  case 22://veg_t3
    return 14*12+month;//169...180
    break;
  case 23://droot_t1
    return 15*12+8;//188
    break;
  case 24://droot_t2
    return 15*12+9;//189
    break;
  case 25://droot_t3
    return 15*12+10;//190
    break;
  case 26://dsoil_t1
    return 15*12+11;//191
    break;
  case 27://dsoil_t2
    return 15*12+12;//192
    break;
  case 28://dsoil_t3
    return 15*12+13;//193
    break;
  case 29://rsmin_t1
    return 15*12+14;//194
    break;
  case 30://rsmin_t2
    return 15*12+15;//195
    break;
  case 31://rsmin_t3
    return 15*12+16;//196
    break;
  case 32://alb_veg_t1
    return 15*12+17;//197
    break;
  case 33://alb_veg_t2
    return 15*12+18;//198
    break;
  case 34://alb_veg_t3
    return 15*12+19;//199
    break;
  case 35://texture
    return 15*12+22;//202

    //return -1;
    break;
  case 36://minlai_t1
    return 15*12+23;//203
    //     cout <<"A derived qtn";//<<endl;
    //return -1;
    break;
  case 37:
    return 15*12+24;//204
    //     cout <<"A derived qtn";//<<endl;
    //return -1;
    break;
  case 38:
    return 15*12+25;//205
    //     cout <<"A derived qtn";//<<endl;
    //return -1;
    break;
  case 39:
    return 15*12+26;//206
    //     cout <<"A derived qtn";//<<endl;
    //return -1;
    break;
  case 40:
    return 15*12+27;//208
    //     cout <<"A derived qtn";//<<endl;
    //return -1;
    break;
  case 41:
    return 15*12+28;//208
    //     cout <<"A derived qtn";//<<endl;
    return -1;
    break;
  case 42://frac_lake
    return 15*12+20;//200
    break;
  case 43://soil_carb
    return 15*12+21;//201
    break;
  case 44: //lake database
    return 209;
    break;

  default:
    cout<<"nn="<<nn<<" is not defined"<<endl;
    return -666;
  }
}


int EcoclimapDataSet::getnn(int variable){
  //return nn corresponding to Patricks RCA list given a variable number,
  //i.e. inverse datamap
  if(variable>=0 && variable < 12) return 1;//lai_t1
  if(variable>=12 && variable < 2*12) return 2;//lai_t2
  if(variable>=2*12 && variable < 3*12) return 3;//lai_t3
  
  if(variable>=3*12 && variable < 4*12) return 4;//z0_t1
  if(variable>=4*12 && variable < 5*12) return 5;//z0_t2
  if(variable>=5*12 && variable < 6*12) return 6;//z0_t3

  if(variable>=6*12 && variable < 7*12) return 7;//emis_t1
  if(variable>=7*12 && variable < 8*12) return 19;//emis_t2
  if(variable>=8*12 && variable < 9*12) return 20;//emis_t3
  
  if(variable>=9*12 && variable < 10*12) return 8;//alb_t1
  if(variable>=10*12 && variable < 11*12) return 17;//alb_t2
  if(variable>=11*12 && variable < 12*12) return 18;//alb_t3

  if(variable>=12*12 && variable < 13*12) return 9;//veg_t1
  if(variable>=13*12 && variable < 14*12) return 21;//veg_t2
  if(variable>=14*12 && variable < 15*12) return 22;//veg_t3

  if(variable==15*12) return 10;//frland
  if(variable==15*12+1) return 11;//alb_soil
  if(variable==15*12+2) return 12;//clay
  if(variable==15*12+3) return 13;//sand

  if(variable==15*12+4) return 14;//frac_t1
  if(variable==15*12+5) return 15;//frac_t2
  if(variable==15*12+6) return 16;//frac_t3

  if(variable==15*12+7) return 23;//d_root_t1
  if(variable==15*12+8) return 24;//d_root_t2
  if(variable==15*12+9) return 25;//d_root_t3

  if(variable==15*12+10) return 26;//d_soil_t1
  if(variable==15*12+11) return 27;//d_soil_t2
  if(variable==15*12+12) return 28;//d_soil_t3

  if(variable==15*12+13) return 29;//rsmin_t1
  if(variable==15*12+14) return 30;//rsmin_t2
  if(variable==15*12+15) return 31;//rsmin_t3

  if(variable==15*12+16) return 32;//alb_veg_t1
  if(variable==15*12+17) return 33;//alb_veg_t2
  if(variable==15*12+18) return 34;//alb_veg_t3
  
  if(variable==15*12+19) return 42;//frac_lake
  if(variable==15*12+20) return 43;//soil_carb
  if(variable==210) return 44;//lakedatabase
  
  if(variable<0 ||variable>210){
    cout <<variable<<" either too large or too small and is not in ecoclimap read files (maybe a derived qnt?)"<<endl;
    return -1;
  }
  return -1;
}



