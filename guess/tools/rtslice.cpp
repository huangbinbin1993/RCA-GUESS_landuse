#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gutil.h>

// Reads and merges a set of output files from RCA-GUESS

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

class Record {

public:
	float lon;
	float lat;
	int id;
	int year;
	int nrec;
	float val[MAXITEM];
	
	Record() {
	
		int i;
		
		nrec=0;
		id=0;
		lon=lat=0.0;
		year=0;
		
		for (i=0;i<MAXITEM;i++) val[i]=0.0;
	}
	
	void add_record(Record& rec) {
	
		// Function to add another record to this record
		
		int i;
		
		if (lon!=rec.lon || lat!=rec.lat || id!=rec.id) {
			printf("add_record: lon, lat or id mismatch\n");
			exit(99);
		}
		
		for (i=0;i<MAXITEM;i++) val[i]+=rec.val[i];
		nrec+=rec.nrec;
	}
	
	void average() {
	
		// Calculates average over number of added records
		
		int i;
		if (nrec)
			for (i=0;i<MAXITEM;i++) val[i]/=(double)nrec;
	}
};

// Global matrix to store data from files
Record data[MAXPROCESS][MAXGRID];

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

void readdata(xtring dir,xtring prefix,xtring var,int minyear,int maxyear,
	xtring label[MAXITEM+4],int& ncol) {

	// Reads data from a set of output files for a particular variable
	// dir= directory path where files are found (e.g. /home/ben/)
	// prefix= prefix part of file name (e.g. test)
	// var= variable part of file name (e.g. lai)
	// minyear, maxyear= timeslice to average over
	// Filenames expected to be called something like "/home/ben/test_lai_0.out"
	
	int nfile,i,j;
	xtring filename,line;
	xtring lcopy[MAXITEM+4];
	int n;
	Record rec;
	FILE* in;
	
	// Add trailing "/" to directory if necessary
	
	if (dir!="") {
		if (dir.mid(dir.len()-1)!="/") dir+="/";
	}
	
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
				
					if (rec.year>=minyear && rec.year<=maxyear) { // within timeslice
					
						if (rec.id>MAXGRID) {
							printf("Too many tiles (%d) in %s\n",rec.id,(char*)filename);
							exit(99);
						}
						
						if (data[i][rec.id-1].id==0) {
							data[i][rec.id-1].id=rec.id;
							data[i][rec.id-1].lon=rec.lon;
							data[i][rec.id-1].lat=rec.lat;
						}
							
						data[i][rec.id-1].add_record(rec);
					}
				}
			}
			
			fclose(in);	
		}
	}

	// Now calculate averages for each timeslice
	
	for (i=0;i<MAXPROCESS;i++)
		for (j=0;j<MAXGRID;j++)
			data[i][j].average();
}

void writedata(xtring prefix,xtring var,int ncol,xtring labels[MAXITEM+4],int minyear,int maxyear) {

	// Outputs data files containing averages for each tile across specified
	// time slice. Separate files for open land (odd ids) and forest (even ids)
	// Note: column "Year" will show the NUMBER of years averaged over
	
	xtring fileop,filefo;
	int i,j,k;
	
	// Open open land file and print header row
	
	fileop.printf("%s_%s_%d-%d_open.txt",(char*)prefix,(char*)var,minyear,maxyear);
	FILE* out_op=fopen(fileop,"wt");
	if (!out_op) {
		printf("Could not open %s for output\n",(char*)fileop);
		exit(99);
	}
	for (i=0;i<ncol;i++) {
		fprintf(out_op,"%s",(char*)labels[i]);
		if (i==ncol-1) fprintf(out_op,"\n");
		else fprintf(out_op,"\t");
	}
	
	// Open forest file and print header row
	
	filefo.printf("%s_%s_%d-%d_for.txt",(char*)prefix,(char*)var,minyear,maxyear);
	FILE* out_fo=fopen(filefo,"wt");
	if (!out_fo) {
		printf("Could not open %s for output\n",(char*)filefo);
		exit(99);
	}
	for (i=0;i<ncol;i++) {
		fprintf(out_fo,"%s",(char*)labels[i]);
		if (i==ncol-1) fprintf(out_fo,"\n");
		else fprintf(out_fo,"\t");
	}
	
	// Print data
	
	for (i=0;i<MAXPROCESS;i++) {
		for (k=0;k<MAXGRID;k++) {
		
			if (data[i][k].id) {
				if (data[i][k].id%2) { // odd id? then open land tile
					fprintf(out_op,"%g\t%g\t%d\t%d",
						data[i][k].lon,data[i][k].lat,data[i][k].id,data[i][k].nrec);
					for (j=0;j<ncol-4;j++) fprintf(out_op,"\t%g",data[i][k].val[j]);
					fprintf(out_op,"\n");
				}
				else { // even id = forest tile
					fprintf(out_fo,"%g\t%g\t%d\t%d",
						data[i][k].lon,data[i][k].lat,data[i][k].id,data[i][k].nrec);
					for (j=0;j<ncol-4;j++) fprintf(out_fo,"\t%g",data[i][k].val[j]);
					fprintf(out_fo,"\n");
				}
			}
		}
	}
	
	fclose(out_op);
	fclose(out_fo);
	
	printf("Output is in files:\n");
	printf("Open land: %s\n",(char*)fileop);
	printf("Forest: %s\n",(char*)filefo);
}


int main(int argc,char* argv[]) {

	// Usage:
	// tslice <directory> <prefix> <variable> <minyear> <maxyear>

	bool argsok=true;
	int minyear,maxyear;
	xtring syear,dir,prefix,var;
	xtring label[MAXITEM+4];
	int ncol;
	
	if (argc!=6) argsok=false;
	
	if (argsok) {
		syear=argv[4];
		if (!syear.isnum()) argsok=false;
		else minyear=syear.num();
		
		syear=argv[5];
		if (!syear.isnum()) argsok=false;
		else maxyear=syear.num();
		
		if (argsok) {
			if (minyear>maxyear) argsok=false;
		}
	}
	
	if (!argsok) {
		printf("Usage: %s <directory> <prefix> <variable> <minyear> <maxyear>\n",
			argv[0]);
		exit(99);
	}
	
	dir=argv[1];
	prefix=argv[2];
	var=argv[3];
	
	readdata(dir,prefix,var,minyear,maxyear,label,ncol);
	writedata(prefix,var,ncol,label,minyear,maxyear);
}
