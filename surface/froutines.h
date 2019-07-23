extern "C" { 
  void getgtopo30path_(char* gtopo30_path);
  void read_dimensions_(int* klon, int* klat);
  void read_domain_(double* south,double* west,double* dlon,double* dlat,double* polon,double* polat,
		    int* klon_global, int* klat_global, int* klev_global);
  void reg2rot_(double* xreg,double* yreg,double* xrot,double* yrot,int* kx, 
		int* ky,double* xcen,double* ycen);
  void rot2reg_(double* xreg,double* yreg,double* xrot,double* yrot,int* kx, 
		int* ky,double* xcen,double* ycen);
  void getecopath_(char* ecopath);
}
