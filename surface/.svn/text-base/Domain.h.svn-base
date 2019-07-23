#ifndef _DOMAIN
#define _DOMAIN
#ifdef MPI_SRC
#include"mpi.h"
#else
#define MPI_Comm int
#endif  

class Domain{
 private:
  double latRotindex_(double lon,double lat);
  double lonRotindex_(double lon,double lat);


 protected:
  char* name;
  double modulo(double a, double b);
  double polygonArea(double *X, double *Y, int points);

 public:

  double latRot_(double i, double j);//{return south + (j-1+jdatastart-1)*dlat;}
  double lonRot_(double i, double j);//{return west + (i-1+idatastart-1)*dlon;}
  char* Name(){return name;}
  void getCell(double i,double j,double lonx[4], double laty[4]);
  void getCellR(double i,double j,double lonx[4], double laty[4]);
  void getCellHardWay(double i,double j,double lonx[4], double laty[4]);

  bool ipos;//is true if x0 + (i-1)*dx, negative if xN - (i-1)*dx
  bool jpos;//is true if y0 + (j-1)*dy, negative if yN - (j-1)*dy
  int klon;
  int klat;
  int halo;
  int idatastart;
  int jdatastart;
  int klon_global;
  int klat_global;
  double west;
  double south;
  double east;
  double north;
  double dlon;
  double dlat;
  double polon;
  double polat;
  bool global;
  bool useHalo;
  ~Domain();

  Domain(Domain *d);
  Domain(char* n, int klon_global,  int klat_global, int klon,  int klat,    
	 double west,  double south,  double dlon,  double dlat,MPI_Comm localCommIn,
	 bool global=false, int idatastart=1,  int jdatastart=1, 
	 double polon=-180.0,  double polat=-90.0, bool ipos=true,bool jpos=true,int halo=1);
  /* void print(char* msg,MPI_Comm localComm); */
  /* void print(char* msg,int myid); */
  void correctCoordinates(double& lon, double& lat);

  double lat(double i, double j);
  double lon(double i, double j);
  double latRotated(double i, double j);
  double lonRotated(double i, double j);
  void rotated(double i, double j,double &xR, double &yR);

  double lonindex(double lon,double lat);
  double latindex(double lon,double lat);
  double lonindexR(double lon,double lat);
  double latindexR(double lon,double lat);

  bool isInside(int i, int j,double vertx[4], double verty[4]);
  bool isInside(double x,double y,double vertx[4], double verty[4]);

};
#endif
