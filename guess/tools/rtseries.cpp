#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gutil.h>
#include <math.h>

// Creates summary time series based on output files from RCA-GUESS

// Format example:
//   Lon   Lat    Id  Year      NE     TBS     IBS      BE   Grass   Total
//  -4.0  48.5     1     0  0.0000  0.0000  0.0000  0.0000  1.7162  1.7162
//  -4.0  48.5     1     1  0.0000  0.0000  0.0000  0.0000  5.7173  5.7173
//  -4.0  48.5     1     2  0.0000  0.0000  0.0000  0.0000  6.0762  6.0762

const int MAXITEM=26;
	// Maximum number of items in a record (row) of an output file
	// (not including longitude,latitude,id and year)
const int MAXGRID=1500;
	// Maximum number of tiles (=unique ids) per process
const int MAXPROCESS=24;
	// Maximum number of processes
const int MAXYEAR=2000;
	// Maximum number of individual years (need not be consecutive) in time series

class Record {

	// Single record in output file

public:
	float lon;
	float lat;
	int id;
	int year;
	int nrec;
	float val[MAXITEM];
	double area;
	
	Record() {
	
		int i;

		id=0;
		lon=lat=0.0;
		year=0;
		nrec=0;
		area=0.0;
		
		for (i=0;i<MAXITEM;i++) val[i]=0.0;
	}
	
	void add_record(Record& rec) {
	
		// Function to add another record to this record
		// Pixel area (for weighting) should be set in rec
		
		int i;
		
		if (year!=rec.year) {
			printf("add_record: year mismatch\n");
			exit(99);
		}
		
		for (i=0;i<MAXITEM;i++) val[i]+=rec.val[i]*rec.area;
		nrec+=rec.nrec;
		area+=rec.area;
	}
	
	void average() {
	
		// Calculates area-based average over number of added records
		
		int i;
		if (area)
			for (i=0;i<MAXITEM;i++) val[i]/=area;
	}
	
};

// Global matrix to store area-averaged data for different years
// 0=open land; 1=forest
Record data[2][MAXYEAR];

// Year number as given in output files from RCA-GUESS
int guessyear[MAXYEAR];

// Number of unique years in output files
int nyear;

double pixelsize(double longpos,double latpos,double longsize,double latsize,int postype) {

	// Returns area in square km of a pixel of a given size at a given point
	// on the world.  The formula applied is the surface area of a segment of
	// a hemisphere of radius r from the equator to a parallel (circular)
	// plane h vertical units towards the pole: S=2*pi*r*h

	// longpos   longitude position (see postype)
	// latpos    latitude position (see postype)
	// longsize  longitude range in degrees
	// latsize   latitude range in degrees
	// postype   declares which part of the pixel longpos and latpos
	//           refer to:
	//           0 = centre
	//           1 = NW corner
	//           2 = NE corner
	//           3 = SW corner
	//           4 = SE corner
      
      double pi,r,h1,h2,lattop,latbot,s;
      
      pi=3.1415926536;
      r=6367.425;   // mean radius of the earth

      lattop=latpos;
      if (postype==0) lattop=latpos+latsize*0.5;
      if (postype==3 || postype==4) lattop=latpos+latsize;
      if (lattop<0.0) lattop=-lattop+latsize;
      latbot=lattop-latsize;
      h1=r*sin(lattop*pi/180.0);
      h2=r*sin(latbot*pi/180.0);
      s=2.0*pi*r*(h1-h2);  //for this latitude band
      return s*longsize/360.0;  //for this pixel
}

bool readheader(FILE*& in,xtring* label,int& ncol) {

	// Reads header row of an RCA-GUESS output file
	// label = array of header labels
	// ncol  = number of columns (labels)
	// Returns false if too many items in file
	
	xtring line;
	int pos;
	
	readfor(in,"a#",&line);
	
	ncol=0;
	pos=line.findnotoneof(" ");
	while (pos!=-1) {
		line=line.mid(pos);
		pos=line.find(' ');
		if (ncol>MAXITEM+2) return false;
		if (pos>0) {
			label[ncol++]=line.left(pos);
			line=line.mid(pos);
			pos=line.findnotoneof(" ");
		}
		else label[ncol++]=line;
	}
	return true;
}

bool readrecord(FILE*& in,Record& rec,int nitem) {

	// Reads one record (row) in output file
	// Returns false on end of file
	// nitem = number of items (not including lon,lat,id and year
	
	xtring fmt;
	double dlon,dlat,dval[MAXITEM];
	int i;
	fmt.printf("f,f,i,i,%df",nitem);
	if (!readfor(in,fmt,&dlon,&dlat,&rec.id,&rec.year,dval)) return false;
	rec.lon=dlon;
	rec.lat=dlat;
	for (i=0;i<nitem;i++) rec.val[i]=dval[i];
	rec.nrec=1;
	return true;
}

int year_number(int year) {

	// Returns index of specified guess year in array guessyear
	// -1 if that year doesn't exist
	
	int y;
	
	for (y=0;y<nyear;y++)
		if (guessyear[y]==year) return y;
		
	return -1;
}

void readdata(xtring dir,xtring prefix,xtring var,xtring label[MAXITEM+4],int& ncol,
	double north,double south,double east,double west) {

	// Reads data from a set of output files for a particular variable
	// dir= directory path where files are found (e.g. /home/ben/)
	// prefix= prefix part of file name (e.g. test)
	// var= variable part of file name (e.g. lai)
	// north,south,east,west = boundaries of window to average over
	// Filenames expected to be called something like "/home/ben/test_lai_0.out"
	
	int nfile,i,j,y,index;
	xtring filename,line;
	xtring lcopy[MAXITEM+4];
	int n;
	Record rec;
	FILE* in;
	
	// Add trailing "/" to directory if necessary
	
	if (dir!="") {
		if (dir.mid(dir.len()-1)!="/") dir+="/";
	}
	
	nyear=0;
	
	nfile=0;
	for (i=0;i<MAXPROCESS;i++) {
		filename.printf("%s%s_%s_%d.out",(char*)dir,(char*)prefix,(char*)var,i);
		in=fopen(filename,"r");
		
		// Check in case empty
		if (in) {
			readfor(in,"a#",&line);
			if (line=="") {
				printf("Warning: %s contains no data!\n",(char*)filename);
				fclose(in);
				in=NULL;
			}
			else {
				fclose(in);
				in=fopen(filename,"rt");
			}
		}
		
		if (in) {
		
			printf("Reading from %s ...\n",(char*)filename);
		
			if (!nfile) {
				
				// Read header row
				
				if (!readheader(in,label,ncol)) {
					printf("Too many items in %s\n",(char*)filename);
					exit(99);
				}
			}
			else {
			
				// Read header row and check that format matches
			
				if (!readheader(in,lcopy,n)) {
					printf("Too many items in %s\n",(char*)filename);
					exit(99);
				}
				if (n!=ncol) {
					printf("Format mismatch in %s - expected %d not %d items\n",
						(char*)filename,ncol,n);
					exit(99);
				}
				for (j=0;j<ncol;j++) {
					if (lcopy[j]!=label[j]) {
						printf("Format mismatch in %s - expected \"%s\" instead of \"%s\"\n",
							(char*)filename,(char*)label[j],(char*)lcopy[j]);
						exit(99);
					}
				}
			}
		
			nfile++;
			
			while (!feof(in)) {
			
				// Read next record in file
			
				if (readrecord(in,rec,ncol-4)) {
				
					if (rec.lon>=west && rec.lon<=east && rec.lat>=south &&
						rec.lat<=north) { // within window
					
						rec.area=pixelsize(rec.lon,rec.lat,0.5,0.5,0);
						
						index=year_number(rec.year);
						if (index==-1) {
							index=nyear;
							nyear++;
							if (nyear>MAXYEAR) {
								printf("Too many years in output files\n");
								exit(99);
							}
							data[0][index].year=rec.year;
							data[1][index].year=rec.year;
							guessyear[index]=rec.year;
						}
						
						if (rec.id>MAXGRID) {
							printf("Too many tiles (%d) in %s\n",rec.id,(char*)filename);
							exit(99);
						}
						
						if (rec.id%2) // odd id - open land
							data[0][index].add_record(rec);
						else  // even id - forest
							data[1][index].add_record(rec);
					}
				}
			}
			
			fclose(in);
		}
	}

	// Now calculate averages for each year
	
	for (i=0;i<nyear;i++) {
		data[0][i].average();
		data[1][i].average();
	}
}

void writedata(xtring prefix,xtring var,int ncol,xtring labels[MAXITEM+4]) {

	// Outputs data files containing averages for each tile across specified
	// time slice. Separate files for open land (odd ids) and forest (even ids)
	// Note: column "Year" will show the NUMBER of years averaged over
	
	xtring fileop,filefo;
	int i,j,k;
	
	labels[1]="nrec";
	labels[2]="area";
	
	// Open open land file and print header row
	
	fileop.printf("%s_%s_open.txt",(char*)prefix,(char*)var);
	FILE* out_op=fopen(fileop,"wt");
	if (!out_op) {
		printf("Could not open %s for output\n",(char*)fileop);
		exit(99);
	}
	for (i=1;i<ncol;i++) {
		fprintf(out_op,"%s",(char*)labels[i]);
		if (i==ncol-1) fprintf(out_op,"\n");
		else fprintf(out_op,"\t");
	}
	
	// Open forest file and print header row
	
	filefo.printf("%s_%s_for.txt",(char*)prefix,(char*)var);
	FILE* out_fo=fopen(filefo,"wt");
	if (!out_fo) {
		printf("Could not open %s for output\n",(char*)filefo);
		exit(99);
	}
	for (i=1;i<ncol;i++) {
		fprintf(out_fo,"%s",(char*)labels[i]);
		if (i==ncol-1) fprintf(out_fo,"\n");
		else fprintf(out_fo,"\t");
	}
	
	// Print data
	
	for (i=0;i<nyear;i++) {
		
		// Open land tile
		fprintf(out_op,"%d\t%g\t%d",data[0][i].nrec,data[0][i].area,data[0][i].year);
		for (j=0;j<ncol-4;j++) fprintf(out_op,"\t%g",data[0][i].val[j]);
		fprintf(out_op,"\n");
			
		// Forest tile
		fprintf(out_fo,"%d\t%g\t%d",data[1][i].nrec,data[1][i].area,data[1][i].year);
		for (j=0;j<ncol-4;j++) fprintf(out_fo,"\t%g",data[1][i].val[j]);
		fprintf(out_fo,"\n");
	}
	
	fclose(out_op);
	fclose(out_fo);
	
	printf("Output is in files:\n");
	printf("Open land: %s\n",(char*)fileop);
	printf("Forest: %s\n",(char*)filefo);
}


int main(int argc,char* argv[]) {

	// Usage:
	// tslice <directory> <prefix> <variable> <west> <south> <east> <north> (pixels within a window)
	// tslice <directory> <prefix> <variable> <lon> <lat>  (one pixel)
	// tslice <directory> <prefix> <variable>  (all pixels)

	bool argsok=true;
	int minyear,maxyear;
	xtring snorth,ssouth,seast,swest,dir,prefix,var;
	xtring label[MAXITEM+4];
	int ncol;
	double north,south,east,west;
	
	if (!(argc==4 || argc==6 || argc==8)) argsok=false;
	
	if (argsok) {
	
		if (argc==4) {
			west=-180.0;
			south=-90.0;
			east=180.0;
			north=90.0;
			printf("Deriving area-averaged time-series for entire modelled domain\n");
		}
		else if (argc==6) {
			swest=argv[4];
			if (!swest.isnum()) argsok=false;
			else west=east=swest.num();			
			ssouth=argv[5];
			if (!ssouth.isnum()) argsok=false;
			else south=north=ssouth.num();
			printf("Deriving area-averaged time-series for grid cell (%g,%g)\n",
				west,south);
		}
		else if (argc==8) {
			swest=argv[4];
			if (!swest.isnum()) argsok=false;
			else west=swest.num();			
			ssouth=argv[5];
			if (!ssouth.isnum()) argsok=false;
			else south=ssouth.num();
			seast=argv[6];
			if (!seast.isnum()) argsok=false;
			else east=seast.num();			
			snorth=argv[7];
			if (!snorth.isnum()) argsok=false;
			else north=snorth.num();
			printf("Deriving area-averaged time-series for window bounded by (%g,%g) and (%g,%g)\n",
				west,south,east,north);
		}
		
		if (argsok) {
			if (east<west || north<south) argsok=false;
		}
	}
	
	if (!argsok) {
		printf("Usage: %s <directory> <prefix> <variable>\n",
			argv[0]);
		printf("       %s <directory> <prefix> <variable> <lon> <lat>  (single pixel)\n",
			argv[0]);
		printf("       %s <directory> <prefix> <variable> <west> <south> <east> <north>  (window)\n",
			argv[0]);

		exit(99);
	}
	
	dir=argv[1];
	prefix=argv[2];
	var=argv[3];
	
	readdata(dir,prefix,var,label,ncol,north,south,east,west);
	writedata(prefix,var,ncol,label);
}
