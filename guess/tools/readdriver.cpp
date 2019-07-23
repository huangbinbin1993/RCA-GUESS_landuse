#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gutil.h>
#include "fastarchive.h"

// Reads control period driver output files created by RCA-GUESS and produces
// fast archive with average values for use in non-feedbacks scenario run

const int MAXGRID=2000;
const int FIRSTCONTROLYEAR=1961;
const int LASTCONTROLYEAR=1990;
const int NYEAR=LASTCONTROLYEAR-FIRSTCONTROLYEAR+1;
const int NARCHIVE=2;
char dir[][256]={
	"./run0","./re1"
};
const int NPROCESS=24;
char prefix[]="echa1b";

double lon[MAXGRID],lat[MAXGRID];
int id[MAXGRID];
float laiphen_grass_opl[MAXGRID][366];
float laiphen_grass_for[MAXGRID][366];
float laiphen_conifer[MAXGRID][366];
float laiphen_broadleaf[MAXGRID][366];
float laimax_conifer[MAXGRID];
float laimax_broadleaf[MAXGRID];
bool yearfound[MAXGRID][LASTCONTROLYEAR-FIRSTCONTROLYEAR+1];
int ngrid;

void readdata(int process) {

	xtring filename;
	int a,n,year,lid,y,d;
	bool found;
	double llon,llat;
	float fval[366],laimax;
	FILE* in;

	ngrid=0;

	for (a=0;a<NARCHIVE;a++) {

		filename.printf("%s/%s_driver_%d.out",dir[a],prefix,process);
		in=fopen(filename,"rb");
		if (!in) {
			printf("Could not open %s for input\n",(char*)filename);
			exit(99);
		}

		printf("Reading data from %s ...\n",(char*)filename);
	
		while (!feof(in)) {

			fread(&llon,sizeof(double),1,in);
			fread(&llat,sizeof(double),1,in);
			fread(&lid,sizeof(int),1,in);
			//printf("lon=%g lat=%g id=%d\n",llon,llat,lid);

			if (!feof(in)) {

				n=0;
				found=false;
				while (n<ngrid && !found) {
					if (lon[n]==llon && lat[n]==llat && lid==id[n]) {
						found=true;
						//printf("Already have (%g,%g,%d)\n",llon,llat,lid);
					}
					else n++;
				}

				if (!found) {

					if (n<MAXGRID) {
						lon[n]=llon;
						lat[n]=llat;
						id[n]=lid;
						laimax_conifer[n]=0.0;
						laimax_broadleaf[n]=0.0;
						for (y=FIRSTCONTROLYEAR;y<=LASTCONTROLYEAR;y++) yearfound[n][y-FIRSTCONTROLYEAR]=false;
						for (d=0;d<366;d++) {
							laiphen_grass_opl[n][d]=0.0;
							laiphen_grass_for[n][d]=0.0;
							laiphen_conifer[n][d]=0.0;
							laiphen_broadleaf[n][d]=0.0;
						}	
						ngrid++;
						//printf("%d: (%g,%g) %d\n",n,lon[n],lat[n],id[n]);
					}
					else {
						printf("Too many grid cells reading %s\n",(char*)filename);
						exit(99);
					}
				}

				fread(&year,sizeof(int),1,in);
				
				if (year<FIRSTCONTROLYEAR || year>LASTCONTROLYEAR) {
					printf("Invalid year (%d) in %s\n",year,(char*)filename);
					exit(99);
				}

				if (yearfound[n][year-FIRSTCONTROLYEAR]) {
					if (year!=1963 || a!=1) printf("WARNING: Duplicate year %d for (%g,%g,%d) in %s\n",year,llon,llat,lid,(char*)filename);
					fread(fval,sizeof(float),366,in);
					fread(fval,sizeof(float),366,in);
					fread(fval,sizeof(float),366,in);
					fread(fval,sizeof(float),366,in);
					fread(fval,sizeof(float),2,in);
					//exit(99);
				}
				else {
					yearfound[n][year-FIRSTCONTROLYEAR]=true;
				
					fread(fval,sizeof(float),366,in);
					for (d=0;d<366;d++) laiphen_grass_opl[n][d]+=fval[d];
					fread(fval,sizeof(float),366,in);
					for (d=0;d<366;d++) laiphen_grass_for[n][d]+=fval[d];
					fread(fval,sizeof(float),366,in);
					for (d=0;d<366;d++) laiphen_conifer[n][d]+=fval[d];
					fread(fval,sizeof(float),366,in);
					for (d=0;d<366;d++) laiphen_broadleaf[n][d]+=fval[d];
					fread(&laimax,sizeof(float),1,in);
					laimax_conifer[n]+=laimax;
					fread(&laimax,sizeof(float),1,in);
					laimax_broadleaf[n]+=laimax;
				}
			}
		}
		
		fclose(in);
	}
	
	for (n=0;n<ngrid;n++) {
		for (y=FIRSTCONTROLYEAR;y<=LASTCONTROLYEAR;y++) {
			if (!yearfound[n][y-FIRSTCONTROLYEAR]) {
				printf("No data for year %d for (%g,%g,%d)\n",y,lon[n],lat[n],id[n]);
				exit(99);
			}
		}
		for (d=0;d<366;d++) {
			laiphen_grass_opl[n][d]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
			laiphen_grass_for[n][d]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
			laiphen_conifer[n][d]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
			laiphen_broadleaf[n][d]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
		}
		laimax_conifer[n]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
		laimax_broadleaf[n]/=(double)(LASTCONTROLYEAR-FIRSTCONTROLYEAR+1);
	}

	printf("Read in data for %d tiles for process %d\n",ngrid,process);
}

void writecheck(int process) {
	
	int n,d;
	xtring filename;
	
	printf("Writing check files for process %d\n",process);
	
	filename.printf("checkmap%d.txt",process);
	FILE* out=fopen(filename,"wt");
	if (!out) {
		printf("Error opening %s for output\n",(char*)filename);
		exit(99);
	}
	
	fprintf(out,"lon\tlat\tid\tlaiphen_grass_opl_0\tlaiphen_grass_opl_180\tlai_phen_grass_opl_364\tlaiphen_grass_for_0\tlaiphen_grass_for_180\tlai_phen_grass_for_364\tlaiphen_conifer_0\tlaiphen_conifer_180\tlai_phen_conifer_364\tlaiphen_broadleaf_0\tlaiphen_broadleaf_180\tlai_phen_broadleaf_364\tlaimax_conifer\tlaimax_broadleaf\n");
	for (n=0;n<ngrid;n++) {
		fprintf(out,"%g\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			lon[n],lat[n],id[n],laiphen_grass_opl[n][0],laiphen_grass_opl[n][180],laiphen_grass_opl[n][364],
			laiphen_grass_for[n][0],laiphen_grass_for[n][180],laiphen_grass_for[n][364],
   			laiphen_conifer[n][0],laiphen_conifer[n][180],laiphen_conifer[n][364],
	  		laiphen_broadleaf[n][0],laiphen_broadleaf[n][180],laiphen_broadleaf[n][364],
	 		laimax_conifer[n],laimax_broadleaf[n]);
	}
	
	fclose(out);
	
	filename.printf("checktile%d.txt",process);
	out=fopen(filename,"wt");
	if (!out) {
		printf("Error opening %s for output\n",(char*)filename);
		exit(99);
	}
	
	fprintf(out,"day\tlaiphen_grass_open\tlaiphen_grass_for\tlaiphen_broadleaf\tlaiphen_conifer\n");
	
	for (d=0;d<366;d++) {
		
		fprintf(out,"%d\t%g\t%g\t%g\t%g\n",
			d,laiphen_grass_opl[0][d],laiphen_grass_for[0][d],laiphen_broadleaf[0][d],
   			laiphen_conifer[0][d]);
	}
	
	fclose(out);
}


int main() {

	int p,n,d;
	double dval,dvec[366];
	
	CFastArchive ark;
	ark.newarchive("GuessVeg");
	ark.defineindex("process_id",0,31,1);
	ark.defineindex("gridcell_id",0,2047,1);
	ark.defineitems("laiphen_grass_opl",366,0,12,0.01);
	ark.defineitems("laiphen_grass_for",366,0,12,0.01);
	ark.defineitems("laiphen_conifer",366,0,12,0.01);
	ark.defineitems("laiphen_broadleaf",366,0,12,0.01);
	ark.defineitems("laimax_conifer",1,0,12,0.01);
	ark.defineitems("laimax_broadleaf",1,0,12,0.01);
	
	for (p=0;p<NPROCESS;p++) {
		readdata(p);
		writecheck(p);
		
		// Write data for this process to archive
		
		for (n=0;n<ngrid;n++) {
			
			ark.setindex("process_id",p);
			ark.setindex("gridcell_id",id[n]);
			for (d=0;d<365;d++) dvec[d]=laiphen_grass_opl[n][d];
			laiphen_grass_opl[n][365]=laiphen_grass_opl[n][364]; // Repeat second last value for leap yuears
			ark.setitem("laiphen_grass_opl",dvec);
			for (d=0;d<365;d++) dvec[d]=laiphen_grass_for[n][d];
			laiphen_grass_for[n][365]=laiphen_grass_for[n][364]; // Repeat second last value for leap yuears
			ark.setitem("laiphen_grass_for",dvec);
			for (d=0;d<365;d++) dvec[d]=laiphen_conifer[n][d];
			laiphen_conifer[n][365]=laiphen_conifer[n][364]; // Repeat second last value for leap yuears
			ark.setitem("laiphen_conifer",dvec);
			for (d=0;d<365;d++) dvec[d]=laiphen_broadleaf[n][d];
			laiphen_broadleaf[n][365]=laiphen_broadleaf[n][364]; // Repeat second last value for leap yuears
			ark.setitem("laiphen_broadleaf",dvec);
			dval=laimax_conifer[n];
			ark.setitem("laimax_conifer",&dval);
			dval=laimax_broadleaf[n];
			ark.setitem("laimax_broadleaf",&dval);
			ark.storerecord();
		}
		
	}

	ark.closearchive();
	
	return 0;
}
