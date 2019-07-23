#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gutil.h>
#include <math.h>
#include "fastarchive.h"

// Program to create a spinup climate driver file for RCA-GUESS based on
// output files for control period, i.e.
// *_temp_*.out             - monthly air temperature
// *_temp_soil_*.out        - monthly soil temperature
// *_mwcont_upper_extra_*.out     - monthly upper soil water
// *_mwcont_lower_extra_*.out     - monthly lower soil water
// *_par_*.out              - monthly PAR
// *_prec_*.out             - monthly precipitation

const int FIRSTYEAR=1981;
const int LASTYEAR=2010;
const int NYEAR=LASTYEAR-FIRSTYEAR+1;
const int MEANFITYEARS=10;
const int MAXGRID=20000;

struct Mdata {
	double lon;
	double lat;
	double temp[12];
	double temp_soil[12];
	double wcont_upper[12];
	double wcont_lower[12];
	double par[12];
	double prec[12];
};


void regress(double* x,double* y,int n,double& a,double& b) {

	// Performs a linear regression of array y on array x (n values)
	// returning parameters a and b in the fitted model: y=a+bx
	// (Used by function soiltemp)
	// Source: Press et al 1986, Sect 14.2

	int i;
	double sx,sy,sxx,sxy,delta;

	sx=0.0;
	sy=0.0;
	sxx=0.0;
	sxy=0.0;
	for (i=0;i<n;i++) {
		sx+=x[i];
		sy+=y[i];
		sxx+=x[i]*x[i];
		sxy+=x[i]*y[i];
	}
	delta=(double)n*sxx-sx*sx;
	a=(sxx*sy-sx*sxy)/delta;
	b=((double)n*sxy-sx*sy)/delta;
}

class Spinup_data {

	// Class for management of climate data for spinup
	// (derived from first few years of historical climate data)

private:
	int nyear;
	int thisyear;
	double* data;
	bool havedata;

public:
	Spinup_data(int nyear_loc) {
		nyear=nyear_loc;
		havedata=false;
		data=new double[nyear*12];
		if (!data) {
			printf("Spinup_data::Spinup_data: out of memory\n");
			exit(99);
		}
		thisyear=0;
		havedata=true;
	}

	~Spinup_data() {
		if (havedata) delete[] data;
	}

	double& operator[](int month) {

		return data[thisyear*12+month];
	}

	void nextyear() {
		if (thisyear==nyear-1) thisyear=0;
		else thisyear++;
	}

	void firstyear() {
		thisyear=0;
	}

	void get_data_from(double source[][12]) {
		
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m]=source[y][m];
			}
		}
	}

	void detrend_data(double set_mean[12]) {

		// Now detrends by month (not year) and fits to specified mean
		
		int y,m;
		double a,b,anomaly;
		double year_number[NYEAR];
		double mval[NYEAR],mean;

		for (m=0;m<12;m++) {
			
			mean=0.0;
			for (y=0;y<NYEAR;y++) {
				year_number[y]=y;
				mval[y]=data[y*12+m];
				mean+=mval[y]/(double)NYEAR;
			}
			
	
			regress(year_number,mval,NYEAR,a,b);
			
			//if (!m) printf("mean=%g set_mean=%g a=%g b=%g\n",mean,set_mean[m],a,b);
	
			for (y=0;y<NYEAR;y++) {
				anomaly=a+b*(double)y;
				//data[y*12+m]-=anomaly+set_mean[m]-mean;
				
			//	if (!m) printf("Year=%d data=%g anomaly=%g new_data=%g\n",
			//		y,data[y*12+m],anomaly,data[y*12+m]+set_mean[m]-anomaly);
				
				data[y*12+m]+=set_mean[m]-anomaly;
			}
		}
	}
};

Mdata forest[MAXGRID][NYEAR];
Mdata open[MAXGRID][NYEAR];
int ngrid_forest,ngrid_open;


void readmonthly(FILE*& in,double& lon,double& lat,int& id,int& year,double* dval1,double* dval2) {

	// Reads one record in a file containing monthly data
	// The expected format is:
	//   <rca-lon> <rca-lat> <id> <year> <jan> ... <dec>
	// or if dval2 is not NULL:
	//   <rca-lon> <rca-lat> <id> <year> <jan1> <jan2> ... <dec1> <dec2>
	
	double dval[24];
	int i;
	
	if (dval2) {
		readfor(in,"f,f,i,i,24f",&lon,&lat,&year,&id,dval);
		for (i=0;i<12;i++) {
			dval1[i]=dval[i*2];
			dval2[i]=dval[i*2+1];
		}
	}
	else {
		readfor(in,"f,f,i,i,12f",&lon,&lat,&year,&id,dval1);
	}
}

int findgrid(Mdata data[MAXGRID][NYEAR],double lon,double lat,int year,int& ngrid) {

	int i;
	for (i=0;i<ngrid;i++) {
		if (fabs(data[i][0].lon-lon)<0.001 && fabs(data[i][0].lat-lat)<0.001) {
			data[i][0].lon=lon;
			data[i][0].lat=lat;
			return i;
		}
	}
	
	if (ngrid<MAXGRID) {
		data[ngrid][0].lon=lon;
		data[ngrid][0].lat=lat;
		ngrid++;
		return ngrid-1;
	}
	else {
		printf("Too many gridpoints\n");
		exit(99);
	}
	
	return -1;
}

void transfer_mdata(Mdata& dest,Mdata& source) {

	int m;
	for (m=0;m<12;m++) {
		dest.temp[m]=source.temp[m];
		dest.temp_soil[m]=source.temp_soil[m];
		dest.wcont_upper[m]=source.wcont_upper[m];
		dest.wcont_lower[m]=source.wcont_lower[m];
		dest.par[m]=source.par[m];
		dest.prec[m]=source.prec[m];
	}
}


void readdata(xtring prefix, int ninst) {

	xtring file_temp,file_temp_soil,file_wcont_upper,file_wcont_lower,file_par,file_prec;
	FILE *in_temp,*in_temp_soil,*in_wcont_upper,*in_wcont_lower,*in_par,*in_prec;
	int i,id,year,index,yearcopy,idcopy;
	bool end_temp_soil,end_wcont_upper,end_wcont_lower,end_par,end_prec;
	double loncopy,latcopy;
	Mdata dat;
	
	ngrid_forest=ngrid_open=0;
	
	for (i=0;i<ninst;i++) {
		file_temp.printf("%s_temp_%d.out",(char*)prefix,i);
		file_temp_soil.printf("%s_temp_soil_%d.out",(char*)prefix,i);
		file_wcont_upper.printf("%s_mwcont_upper_extra_%d.out",(char*)prefix,i);
		file_wcont_lower.printf("%s_mwcont_lower_extra_%d.out",(char*)prefix,i);
		file_par.printf("%s_par_%d.out",(char*)prefix,i);
		file_prec.printf("%s_prec_%d.out",(char*)prefix,i);
		
		in_temp=fopen(file_temp,"rt");
		if (!in_temp) {
			printf("Could not open %s for input\n",(char*)file_temp);
			exit(99);
		}
		readfor(in_temp,"");
		
		in_temp_soil=fopen(file_temp_soil,"rt");
		if (!in_temp_soil) {
			printf("Could not open %s for input\n",(char*)file_temp_soil);
			exit(99);
		}
		readfor(in_temp_soil,"");
		
		in_wcont_upper=fopen(file_wcont_upper,"rt");
		if (!in_wcont_upper) {
			printf("Could not open %s for input\n",(char*)file_wcont_upper);
			exit(99);
		}
		readfor(in_wcont_upper,"");

		in_wcont_lower=fopen(file_wcont_lower,"rt");
		if (!in_wcont_lower) {
			printf("Could not open %s for input\n",(char*)file_wcont_lower);
			exit(99);
		}
		readfor(in_wcont_lower,"");
		
		in_par=fopen(file_par,"rt");
		if (!in_par) {
			printf("Could not open %s for input\n",(char*)file_par);
			exit(99);
		}
		readfor(in_par,"");

		in_prec=fopen(file_prec,"rt");
		if (!in_prec) {
			printf("Could not open %s for input\n",(char*)file_prec);
			exit(99);
		}
		readfor(in_prec,"");
				
		printf("Instance %d\n",i);
		
		end_temp_soil=end_wcont_upper=end_wcont_lower=end_par=end_prec=false;
		
		while (!feof(in_temp)) {
			readmonthly(in_temp,dat.lon,dat.lat,id,year,dat.temp,NULL);

			if (!feof(in_temp)) {
				if (id%2) index=findgrid(open,dat.lon,dat.lat,year,ngrid_open);
				else index=findgrid(forest,dat.lon,dat.lat,year,ngrid_forest);

				if (!end_temp_soil) {
					readmonthly(in_temp_soil,loncopy,latcopy,idcopy,yearcopy,dat.temp_soil,NULL);
					if (feof(in_temp_soil)) {
						printf("Unexpected end-of-file in %s\n",(char*)file_temp_soil);
						end_temp_soil=true;
						//exit(99);
					}
					else if (loncopy!=dat.lon || latcopy!=dat.lat || idcopy!=id || yearcopy!=year) {
						printf("mismatched record in %s: %g/%g %g/%g %d/%d %d/%d\n",
							(char*)file_temp_soil,loncopy,dat.lon,latcopy,dat.lat,idcopy,id,
							yearcopy,year);
						exit(99);
					}
				}
				
				if (!end_wcont_upper) {
					readmonthly(in_wcont_upper,loncopy,latcopy,idcopy,yearcopy,dat.wcont_upper,NULL);
					if (feof(in_wcont_upper)) {
						printf("Unexpected end-of-file in %s\n",(char*)file_wcont_upper);
						end_wcont_upper=true;
						//exit(99);
					}
					else if (loncopy!=dat.lon || latcopy!=dat.lat || idcopy!=id || yearcopy!=year) {
						printf("mismatched record in %s: %g/%g %g/%g %d/%d %d/%d\n",
							(char*)file_wcont_upper,loncopy,dat.lon,latcopy,dat.lat,idcopy,id,
							yearcopy,year);
						exit(99);
					}
				}

				if (!end_wcont_lower) {
					readmonthly(in_wcont_lower,loncopy,latcopy,idcopy,yearcopy,dat.wcont_lower,NULL);
					if (feof(in_wcont_lower)) {
						printf("Unexpected end-of-file in %s\n",(char*)file_wcont_lower);
						end_wcont_lower=true;
						//exit(99);
					}
					else if (loncopy!=dat.lon || latcopy!=dat.lat || idcopy!=id || yearcopy!=year) {
						printf("mismatched record in %s: %g/%g %g/%g %d/%d %d/%d\n",
							(char*)file_wcont_lower,loncopy,dat.lon,latcopy,dat.lat,idcopy,id,
							yearcopy,year);
						exit(99);
					}
				}
				
				if (!end_par) {
					readmonthly(in_par,loncopy,latcopy,idcopy,yearcopy,dat.par,NULL);
					if (feof(in_par)) {
						printf("Unexpected end-of-file in %s\n",(char*)file_par);
						end_par=true;
						//exit(99);
					}
					else if (loncopy!=dat.lon || latcopy!=dat.lat || idcopy!=id || yearcopy!=year) {
						printf("mismatched record in %s: %g/%g %g/%g %d/%d %d/%d\n",
							(char*)file_par,loncopy,dat.lon,latcopy,dat.lat,idcopy,id,
							yearcopy,year);
						exit(99);
					}
				}

				if (!end_prec) {
					readmonthly(in_prec,loncopy,latcopy,idcopy,yearcopy,dat.prec,NULL);
					if (feof(in_prec)) {
						printf("Unexpected end-of-file in %s\n",(char*)file_prec);
						end_prec=true;
						//exit(99);
					}
					else if (loncopy!=dat.lon || latcopy!=dat.lat || idcopy!=id || yearcopy!=year) {
						printf("mismatched record in %s: %g/%g %g/%g %d/%d %d/%d\n",
							(char*)file_prec,loncopy,dat.lon,latcopy,dat.lat,idcopy,id,
							yearcopy,year);
						exit(99);
					}
				}
									
				if (year>=FIRSTYEAR && year<=LASTYEAR) {
					//printf("lon=%g lat=%g year=%d index=%d\n",dat.lon,dat.lat,year,index);
					if (id%2) transfer_mdata(open[index][year-FIRSTYEAR],dat);
					else transfer_mdata(forest[index][year-FIRSTYEAR],dat);
				}
			}
		}
		
		fclose(in_temp);
		fclose(in_temp_soil);
		fclose(in_wcont_upper);
		fclose(in_wcont_lower);
		fclose(in_par);
		fclose(in_prec);
		
		printf("ngrid=%d (%d)\n",ngrid_forest,ngrid_open);
	}
}


void transform(Mdata dat[NYEAR]) {

	// Detrends and rescales data to mean for first 10 years of data set
	
	// temp, temp_soil - detrended arithmetically
	// wcont, par, prec - not detrended
	
	// Ben 2007-12-10: tested detrending - works!
	
	Spinup_data mdata(NYEAR);
	double mean[12];
	double data[NYEAR][12];
	int m,y;
	
	// Temperature
	
	for (m=0;m<12;m++)
		mean[m]=0.0;
	
	for (m=0;m<12;m++) {
		for (y=0;y<MEANFITYEARS;y++)
			mean[m]+=dat[y].temp[m]/(double)MEANFITYEARS;
	}
	
	for (y=0;y<NYEAR;y++)
		for (m=0;m<12;m++)
			data[y][m]=dat[y].temp[m];
	
	//printf("mean=%g\n",mean[0]);
	mdata.get_data_from(data);
	mdata.detrend_data(mean);
	
	mdata.firstyear();
	for (y=0;y<NYEAR;y++) {
		for (m=0;m<12;m++)
			dat[y].temp[m]=mdata[m];
		mdata.nextyear();
	}
	
	// Soil temperature
	
	for (m=0;m<12;m++)
		mean[m]=0.0;
	
	for (m=0;m<12;m++) {
		for (y=0;y<MEANFITYEARS;y++)
			mean[m]+=dat[y].temp_soil[m]/(double)MEANFITYEARS;
	}
	
	for (y=0;y<NYEAR;y++)
		for (m=0;m<12;m++)
			data[y][m]=dat[y].temp_soil[m];
	
	mdata.get_data_from(data);
	mdata.detrend_data(mean);
	
	mdata.firstyear();
	for (y=0;y<NYEAR;y++) {
		for (m=0;m<12;m++)
			dat[y].temp_soil[m]=mdata[m];
		mdata.nextyear();
	}
}

void writedata() {

	int i,y,m;
	Mdata dat;
	
	FILE* out=fopen("result.txt","wt");
	FILE* out_july=fopen("result_july.txt","wt");
	
	if (!out) {
		printf("Could not open output file\n");
		exit(99);
	}
	
	fprintf(out,"Lon\tLat");
	fprintf(out_july,"Lon\tLat");
	
	for (y=FIRSTYEAR;y<=LASTYEAR;y++) {
		fprintf(out,"\tJan%d",y);
		fprintf(out,"\tFeb%d",y);
		fprintf(out,"\tMar%d",y);
		fprintf(out,"\tApr%d",y);
		fprintf(out,"\tMay%d",y);
		fprintf(out,"\tJun%d",y);
		fprintf(out,"\tJul%d",y);
		fprintf(out,"\tAug%d",y);
		fprintf(out,"\tSep%d",y);
		fprintf(out,"\tOct%d",y);
		fprintf(out,"\tNov%d",y);
		fprintf(out,"\tDec%d",y);
		fprintf(out_july,"\tJul%d",y);
	}
	fprintf(out,"\n");
	fprintf(out_july,"\n");
		
	for (i=0;i<ngrid_forest;i++) {
		fprintf(out,"%8.3f\t%8.3f",forest[i][0].lon,forest[i][0].lat);
		fprintf(out_july,"%8.3f\t%8.3f",forest[i][0].lon,forest[i][0].lat);
		
		for (y=0;y<NYEAR;y++) {
			for (m=0;m<12;m++)
				fprintf(out,"\t%7.1f",forest[i][y].prec[m]);
			fprintf(out_july,"\t%0.10f",forest[i][y].prec[6]);
		}
		fprintf(out,"\n");
		fprintf(out_july,"\n");
		
		transform(forest[i]);
		
		fprintf(out,"%8.3f\t%8.3f",forest[i][0].lon,forest[i][0].lat);
		fprintf(out_july,"%8.3f\t%8.3f",forest[i][0].lon,forest[i][0].lat);
		
		for (y=0;y<NYEAR;y++) {
			for (m=0;m<12;m++)
				fprintf(out,"\t%7.1f",forest[i][y].prec[m]);
			fprintf(out_july,"\t%0.10f",forest[i][y].prec[6]);		
		}
		fprintf(out,"\n");
		fprintf(out_july,"\n");
	}
	
	fclose(out);
}

void writebin() {

	// Writes output data to FastArchive binary archive (similar to CRU data file)
	
	bool desert;
	
	if (ngrid_forest!=ngrid_open) {
		printf("Error: ngrid_forest (%d) must equal ngrid_open (%d)\n",
			ngrid_forest,ngrid_open);
		exit(99);
	}
	
	int i,m,y,ct;
	double dval_for[12*NYEAR],dval_open[12*NYEAR];
	CFastArchive ark;
	ark.newarchive("RCAData");
	
	ark.defineindex("lon",-180.0,180.0,0.25);
	ark.defineindex("lat",-90.0,90.0,0.25);
	ark.defineitems("mtemp_for",12*NYEAR,-100.0,100.0,0.1);
	ark.defineitems("mtemp_soil_for",12*NYEAR,-100.0,100.0,0.1);
	ark.defineitems("mwcont_upper_for",12*NYEAR,0,1,0.01);
	ark.defineitems("mwcont_lower_for",12*NYEAR,0,1,0.01);
	ark.defineitems("mpar_for",12*NYEAR,0,20000,10);
	ark.defineitems("mprec_for",12*NYEAR,0,10000,0.1);
	ark.defineitems("mtemp_open",12*NYEAR,-100.0,100.0,0.1);
	ark.defineitems("mtemp_soil_open",12*NYEAR,-100.0,100.0,0.1);
	ark.defineitems("mwcont_upper_open",12*NYEAR,0,1,0.01);
	ark.defineitems("mwcont_lower_open",12*NYEAR,0,1,0.01);
	ark.defineitems("mpar_open",12*NYEAR,0,20000,10);
	ark.defineitems("mprec_open",12*NYEAR,0,10000,0.1);
	
	ct=0;
	
	for (i=0;i<ngrid_forest;i++) {
		ark.setindex("lon",forest[i][0].lon);
		ark.setindex("lat",forest[i][0].lat);
		
		// Detrend temperature and soil temperature
		transform(forest[i]);
		transform(open[i]);
		
		// Temperature
		
		for (y=0;y<NYEAR;y++)
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].temp[m];
				dval_open[y*12+m]=open[i][y].temp[m];
			}
		ark.setitem("mtemp_for",dval_for);
		ark.setitem("mtemp_open",dval_open);
		
		// Soil temperature
		
		for (y=0;y<NYEAR;y++)
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].temp_soil[m];
				dval_open[y*12+m]=open[i][y].temp_soil[m];
			}
		ark.setitem("mtemp_soil_for",dval_for);
		ark.setitem("mtemp_soil_open",dval_open);
		
		// Soil water
		// Added quick fix for 'undefined' soil water in Sahara
		// If sum of upper and lower layer AWC for all months >23.9
		// then set to 0.01 instead
		
		desert=false;
		for (y=0;y<NYEAR;y++) {
			double sum=0;
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].wcont_upper[m];
				dval_open[y*12+m]=open[i][y].wcont_upper[m];
				sum+=open[i][y].wcont_upper[m]+open[i][y].wcont_lower[m];
			}
			if (sum>23.9 && forest[i][0].lat<38.0) desert=true;
		}
		
		if (desert) {
			for (y=0;y<NYEAR;y++)
				for (m=0;m<12;m++)
					dval_for[y*12+m]=dval_open[y*12+m]=0.01;
		}
		
		ark.setitem("mwcont_upper_for",dval_for);
		ark.setitem("mwcont_upper_open",dval_open);
		
		desert=false;
		for (y=0;y<NYEAR;y++) {
			double sum=0;
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].wcont_lower[m];
				dval_open[y*12+m]=open[i][y].wcont_lower[m];
				sum+=open[i][y].wcont_upper[m]+open[i][y].wcont_lower[m];
			}
			if (sum>23.9 && forest[i][0].lat<38.0) desert=true;
		}
		
		if (desert) {
			for (y=0;y<NYEAR;y++)
				for (m=0;m<12;m++)
					dval_for[y*12+m]=dval_open[y*12+m]=0.01;
		}
		
		ark.setitem("mwcont_lower_for",dval_for);
		ark.setitem("mwcont_lower_open",dval_open);
		
		// PAR
		
		for (y=0;y<NYEAR;y++)
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].par[m];
				dval_open[y*12+m]=open[i][y].par[m];
			}
		ark.setitem("mpar_for",dval_for);
		ark.setitem("mpar_open",dval_open);

		// PREC
		
		for (y=0;y<NYEAR;y++)
			for (m=0;m<12;m++) {
				dval_for[y*12+m]=forest[i][y].prec[m];
				dval_open[y*12+m]=open[i][y].prec[m];
			}
		ark.setitem("mprec_for",dval_for);
		ark.setitem("mprec_open",dval_open);
		
		ark.storerecord();
		
		ct++;
		if (!(ct%1000)) printf("Written record %d/%d\n",ct,ngrid_forest);
	}
	
	ark.closearchive();
}

void test() {

	Mdata dat[NYEAR];
	double data[NYEAR]={15.5,15.9,17.3,17.1,15,15.4,15.6,15.4,15.9,
		12.5,15.1,13.2,15.8,14.7,15.8,17.8,13.4,14.4,15.1,14.1};
	
	int y,m;
	for (y=0;y<NYEAR;y++)
		dat[y].temp[0]=data[y];
	
	transform(dat);
	
	for (y=0;y<NYEAR;y++)
		printf("%d\t%g\n",y,dat[y].temp[0]);
	
}

int main(int argc, char** argv) {

	//test();
	//exit(99);

	if (argc != 3) {
		fprintf(stderr, "Usage: calc <prefix> <number-of-instances>\n");
		exit(99);
	}

	int ninst = atoi(argv[2]);

	if (ninst <= 0) {
		fprintf(stderr, "Please specify a positive integer as second argument\n");
		exit(99);
	}

	readdata(argv[1], ninst);
	
	if (ngrid_open!=ngrid_forest) {
		printf("Mismatch in number of grid cells for open land (%d) and forest (%d)\n",
			ngrid_open,ngrid_forest);
		exit(99);
	}
	
	writebin();
	//writedata();
	
	return 0;
}
