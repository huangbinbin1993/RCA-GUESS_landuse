#include<math.h>
#include<iostream>
union tal{
  float real;
  int integer;
};
union talF{
  float real;
  int integer;
};
using namespace std;
extern "C" {
  void r2ibts_(float *field, int *ifield){
    union tal convert;
    convert.real = *field;
    *ifield = convert.integer;
  }

  void i2rbts_(float *field, int *ifield){
    union tal convert;
    convert.integer = *ifield; 
    *field = convert.real;
  }

  void r2ibtsv_(float *field, int *ifield,int *len){
    for(int i=0;i<*len;i++)
      r2ibts_(&field[i], &ifield[i]); 
  }

  void i2rbtsv_(float *field, int *ifield,int *len){
    for(int i=0;i<*len;i++)
      i2rbts_( &field[i], &ifield[i]); 
  }

  void exscalc_(float *pdata, int *klen, float *pref, float *pscale){
  
    union tal convert;

    for(int i=0;i<*klen;i++){
      //r2ibts_(&pdata[i],&pdata[i]);
      convert.real = pdata[i];
      pdata[i] = convert.integer;
    }

    for(int i=0;i<*klen;i++){
      pdata[i] = *pref + pdata[i]*(*pscale);
    }
  }
  
  void rorinti_(float *pdata,int *klen){
    for(int i=0;i<*klen;i++){
      pdata[i] = roundf(pdata[i]);
    }
  }


  void rorintr_(float *pdata,int *klen){
    union tal convert;
    for(int i=0;i<*klen;i++){
      convert.real = pdata[i];
      pdata[i] = (float)convert.integer;
    }
  }

  void r2ibtsf_(float *field, int *ifield){
    union talF convert;
    convert.real = *field;
    *ifield = convert.integer;
  }

  void i2rbtsf_(float *field, int *ifield){
    union talF convert;
    convert.integer = *ifield; 
    *field = convert.real;
  }

  void r2ibtsvf_(float *field, int *ifield,int *len){
    for(int i=0;i<*len;i++)
      r2ibtsf_(&field[i], &ifield[i]); 
  }

  void i2rbtsvf_(float *field, int *ifield,int *len){
    for(int i=0;i<*len;i++)
      i2rbtsf_( &field[i], &ifield[i]); 
  }

  void exscalcf_(float *pdata, int *klen, float *pref, float *pscale){
  
    union talF convert;

    for(int i=0;i<*klen;i++){
      convert.real = pdata[i];
      pdata[i] = convert.integer;
    }

    for(int i=0;i<*klen;i++){
      pdata[i] = *pref + pdata[i]*(*pscale);
    }
  }
  
  void rorintif_(float *pdata,int *klen){
    for(int i=0;i<*klen;i++){
      pdata[i] = roundf(pdata[i]);
    }
  }


  void rorintrf_(float *pdata,int *klen){
    union talF convert;
    for(int i=0;i<*klen;i++){
      convert.real = pdata[i];
      pdata[i] = (float)convert.integer;
    }
  }


}
