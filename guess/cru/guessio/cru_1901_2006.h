//////////////////////////////////////////////////////////////////////////////////////
// CRU_1901_2006.H
// Header file for input from a fast data archive
// Created automatically by FastArchive on Tue Jun 09 11:59:28 2009
//
// Version dated:         2010-11-22
//
//
// The following #includes should appear in your source code file:
//
//   #include <stdio.h>
//   #include <stdlib.h>
//   #include <string.h>
//   #include "E://CRU_TS_3p0//cru_1901_2006.h"
//
// Functionality to retrieve data from the archive is provided by class Cru_1901_2006Archive.
// The following public functions are provided:
//
// bool open(char* filename)
//   Attempts to open the specified file as a fast data archive. The format must be
//   exactly compatible with this version of cru_1901_2006.h (normally the archive and
//   header file should have been produced together by the same program using class
//   CFastArchive). Returns false if the file could not be opened or had format
//   errors. open() with no argument is equivalent to open("cru_1901_2006.bin").
//
// void close()
//   Closes the archive (if open).
//
// bool rewind()
//   Sets the file pointer to the first record in the archive file. Returns false if
//   no archive file is currently open.
//
// bool getnext(Cru_1901_2006& obj)
//   Retrieves the next record in the archive file and advances the file pointer to
//   the next record. Data are written to the member variables of obj. Returns false if
//   no archive file is currently open or if the file pointer is beyond the last
//   record. Use rewind() and getnext() to retrieve data sequentially from the archive.
//
// bool getindex(Cru_1901_2006& obj)
//   Searches the archive for a record matching the values specified for the index
//   items (lon and lat) in obj. If a matching record is found, the data are
//   written to the member variables of obj. Returns true if the archive was open and
//   a matching record was found, otherwise false. The search is iterative and fast.
//
// Sample program:
//
//   Cru_1901_2006Archive ark;
//   Cru_1901_2006 data;
//   bool success,flag;
//
//   // Retrieve all records in sequence and print values of lon and lat:
//
//   success=ark.open("cru_1901_2006.bin");
//   if (success) {
//      flag=ark.rewind();
//      while (flag) {
//         flag=ark.getnext(data);
//         if (flag)
//            printf("Loaded record: lon=%g, lat=%g\n",data.lon,data.lat);
//      }
//   }
//   
//   // Look for a record with lon=-1800, lat=-900:
//
//   data.lon=-1800;
//   data.lat=-900;
//   success=ark.getindex(data);
//   if (success) printf("Found it!\n");
//   else printf("Not found\n");
//
//   ark.close();


struct Cru_1901_2006 {

	// Index part

	double lon;
	double lat;

	// Data part

	double soilcode[1];
	double mtemp[1272];
	double mprec[1272];
	double msun[1272];
};


const long CRU_1901_2006_NRECORD=59191;
const int CRU_1901_2006_DATA_LENGTH=6043;
const int CRU_1901_2006_INDEX_LENGTH=7;
const int CRU_1901_2006_HEADERSIZE=657;
const unsigned char CRU_1901_2006_HEADER[CRU_1901_2006_HEADERSIZE-4]={
	0x01,0x02,0x91,0x00,0x00,0x00,0x07,0x00,0x02,0x04,0x6C,0x6F,0x6E,0x00,0x2D,0x31,0x38,0x30,0x30,0x00,
	0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0x90,0x07,0x35,0x00,0xE4,0x25,
	0x40,0x00,0x40,0x5D,0x42,0x00,0x31,0x38,0x30,0x30,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,
	0x00,0x00,0x04,0x00,0x00,0x00,0x90,0x07,0x35,0x00,0xE4,0x25,0x40,0x00,0x40,0x5D,0x42,0x00,0x31,0x00,
	0x30,0x30,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0x90,0x07,
	0x35,0x00,0xE4,0x25,0x40,0x00,0x40,0x5D,0x42,0x00,0x0C,0x04,0x6C,0x61,0x74,0x00,0x2D,0x39,0x30,0x30,
	0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0xA0,0x11,0x35,0x00,
	0xE4,0x25,0x40,0x00,0x3C,0x5D,0x42,0x00,0x39,0x30,0x30,0x00,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,
	0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,0xA0,0x11,0x35,0x00,0xE4,0x25,0x40,0x00,0x3C,0x5D,0x42,0x00,
	0x31,0x00,0x30,0x00,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x04,0x00,0x00,0x00,
	0xA0,0x11,0x35,0x00,0xE4,0x25,0x40,0x00,0x3C,0x5D,0x42,0x00,0x0B,0x00,0x00,0x17,0x9B,0x00,0x04,0x09,
	0x73,0x6F,0x69,0x6C,0x63,0x6F,0x64,0x65,0x00,0x30,0x00,0x00,0x00,0x09,0x00,0x00,0x00,0x68,0xF6,0x12,
	0x00,0x00,0x00,0x00,0x00,0x09,0x00,0x00,0x00,0x28,0x30,0x35,0x00,0xFF,0x26,0x40,0x00,0x35,0x5D,0x42,
	0x00,0x31,0x32,0x00,0x00,0x09,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x09,0x00,0x00,
	0x00,0x28,0x30,0x35,0x00,0xFF,0x26,0x40,0x00,0x35,0x5D,0x42,0x00,0x31,0x00,0x00,0x00,0x09,0x00,0x00,
	0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x09,0x00,0x00,0x00,0x28,0x30,0x35,0x00,0xFF,0x26,0x40,
	0x00,0x35,0x5D,0x42,0x00,0x04,0x00,0x00,0x00,0x01,0x06,0x6D,0x74,0x65,0x6D,0x70,0x00,0x2D,0x31,0x30,
	0x30,0x30,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xD0,0x30,0x35,
	0x00,0xFF,0x26,0x40,0x00,0x2A,0x5D,0x42,0x00,0x31,0x30,0x30,0x30,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,
	0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xD0,0x30,0x35,0x00,0xFF,0x26,0x40,0x00,0x2A,0x5D,0x42,
	0x00,0x31,0x00,0x30,0x30,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,
	0x00,0xD0,0x30,0x35,0x00,0xFF,0x26,0x40,0x00,0x2A,0x5D,0x42,0x00,0x0B,0x00,0x00,0x04,0xF8,0x06,0x6D,
	0x70,0x72,0x65,0x63,0x00,0x2D,0x31,0x30,0x30,0x30,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,
	0x00,0x06,0x00,0x00,0x00,0xE8,0x4B,0x35,0x00,0xFF,0x26,0x40,0x00,0x22,0x5D,0x42,0x00,0x39,0x39,0x39,
	0x39,0x39,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xE8,0x4B,0x35,
	0x00,0xFF,0x26,0x40,0x00,0x22,0x5D,0x42,0x00,0x31,0x00,0x39,0x39,0x39,0x00,0x00,0x00,0x68,0xF6,0x12,
	0x00,0x00,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0xE8,0x4B,0x35,0x00,0xFF,0x26,0x40,0x00,0x22,0x5D,0x42,
	0x00,0x11,0x00,0x00,0x04,0xF8,0x05,0x6D,0x73,0x75,0x6E,0x00,0x2D,0x31,0x00,0x00,0x05,0x00,0x00,0x00,
	0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x05,0x00,0x00,0x00,0x00,0x45,0x35,0x00,0xFF,0x26,0x40,0x00,
	0x19,0x5D,0x42,0x00,0x31,0x30,0x30,0x30,0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,
	0x05,0x00,0x00,0x00,0x00,0x45,0x35,0x00,0xFF,0x26,0x40,0x00,0x19,0x5D,0x42,0x00,0x31,0x00,0x30,0x30,
	0x00,0x00,0x00,0x00,0x68,0xF6,0x12,0x00,0x00,0x00,0x00,0x00,0x05,0x00,0x00,0x00,0x00,0x45,0x35,0x00,
	0xFF,0x26,0x40,0x00,0x19,0x5D,0x42,0x00,0x0A,0x00,0x00,0x04,0xF8};


class Cru_1901_2006Archive {

private:

	FILE* pfile;
	long recno;
	long datano;
	unsigned char pindex[CRU_1901_2006_INDEX_LENGTH];
	unsigned char pdata[CRU_1901_2006_DATA_LENGTH];
	bool iseof;

	long readbin(int nbyte) {

		unsigned char buf[4];
		long mult[4]={0x1000000,0x10000,0x100,1},val;
		int i;

		fread(buf,nbyte,1,pfile);

		val=0;
		for (i=0;i<nbyte;i++) {
			val+=buf[i]*mult[4-nbyte+i];
		}

		return val;
	}

	void getindex(long n) {

		fseek(pfile,CRU_1901_2006_INDEX_LENGTH*(n-recno),SEEK_CUR);
		fread(pindex,CRU_1901_2006_INDEX_LENGTH,1,pfile);
		datano=pindex[CRU_1901_2006_INDEX_LENGTH-4]*0x1000000+pindex[CRU_1901_2006_INDEX_LENGTH-3]*0x10000+
			pindex[CRU_1901_2006_INDEX_LENGTH-2]*0x100+pindex[CRU_1901_2006_INDEX_LENGTH-1];
		recno=n+1;
		iseof=(recno==CRU_1901_2006_NRECORD);
	}

	void getdata() {

		fseek(pfile,CRU_1901_2006_INDEX_LENGTH*-recno+(datano-CRU_1901_2006_NRECORD)*CRU_1901_2006_DATA_LENGTH,SEEK_CUR);
		fread(pdata,CRU_1901_2006_DATA_LENGTH,1,pfile);
		fseek(pfile,CRU_1901_2006_DATA_LENGTH*(CRU_1901_2006_NRECORD-datano-1)+CRU_1901_2006_INDEX_LENGTH*recno,SEEK_CUR);
	}

	double popreal(unsigned char* bits,int nbyte,int nbit,double scalar,double offset) {

		unsigned char buf;
		int nb=nbit/8,i;
		double rval=0.0;
		long mult[4]={1,0x100,0x10000,0x1000000};

		for (i=0;i<4;i++) {
			if (i<nb) rval+=bits[nbyte-i-1]*mult[i];
			else if (i==nb) {
				buf=bits[nbyte-i-1]<<(8-nbit%8);
				buf>>=8-nbit%8;
				rval+=buf*mult[i];
			}
		}

		for (i=nbyte-1;i>=0;i--) {
			if (i>=nb)
				bits[i]=bits[i-nb];
			else
				bits[i]=0;
		}

		nb=nbit%8;

		for (i=nbyte-1;i>=0;i--) {
			bits[i]>>=nb;
			if (i>0) {
				buf=bits[i-1];
				buf<<=8-nb;
				bits[i]|=buf;
			}
		}

		rval=rval*scalar+offset;

		return rval;
	}

	void bitify(unsigned char buf[4],double fval,double offset,double scalar) {

		long ival = (long)((fval-offset)/scalar + 0.5);
		buf[0]=(unsigned char)(ival/0x1000000);
		ival-=buf[0]*0x1000000;
		buf[1]=(unsigned char)(ival/0x10000);
		ival-=buf[1]*0x10000;
		buf[2]=(unsigned char)(ival/0x100);
		ival-=buf[2]*0x100;
		buf[3]=(unsigned char)(ival);
	}

	void merge(unsigned char ptarget[3],unsigned char buf[4],int bits) {

		int nb=bits/8;
		int i;
		unsigned char nib;
		for (i=0;i<3;i++) {

			if (i<3-nb)
				ptarget[i]=ptarget[i+nb];
			else
				ptarget[i]=0;
		}
		nb=bits%8;
		for (i=0;i<3;i++) {
			ptarget[i]<<=nb;
			if (i<3-1) {
				nib=ptarget[i+1]>>(8-nb);
				ptarget[i]|=nib;
			}
		}

		nb=bits/8;
		if (bits%8) nb++;
		for (i=1;i<=nb;i++)
			ptarget[3-i]|=buf[4-i];
	}

	int compare_index(unsigned char* a,unsigned char* b) {

		int i;
		for (i=0;i<3;i++) {
			if (a[i]<b[i]) return -1;
			else if (a[i]>b[i]) return +1;
		}

		return 0;
	}

	bool initialise(const char* filename) {

		int i;
		unsigned char* pheader;

		if (pfile) fclose(pfile);
		pfile=fopen(filename,"rb");
		if (!pfile) {
			printf("Could not open %s for input\n",filename);
			return false;
		}

		pheader=new unsigned char[CRU_1901_2006_HEADERSIZE-4];
		if (!pheader) {
			printf("Out of memory\n");
			fclose(pfile);
			pfile=NULL;
			return false;
		}
		::rewind(pfile);
		fread(pheader,CRU_1901_2006_HEADERSIZE-4,1,pfile);
		for (i=0;i<CRU_1901_2006_HEADERSIZE-4;i++) {
			if (pheader[i]!=CRU_1901_2006_HEADER[i]) {
				printf("Format of %s incompatible with this version of cru_1901_2006.h\n",filename);
				fclose(pfile);
				pfile=NULL;
				delete pheader;
				return false;
			}
		}
		delete pheader;

		::rewind(pfile);
		fseek(pfile,CRU_1901_2006_HEADERSIZE+CRU_1901_2006_DATA_LENGTH*CRU_1901_2006_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

public:

	Cru_1901_2006Archive() {
		pfile=NULL;
	}

	~Cru_1901_2006Archive() {
		if (pfile) fclose(pfile);
	}

	bool open(const char* filename) {
		return initialise(filename);
	}

	bool open() {
		return open("cru_1901_2006.bin");
	}

	void close() {
		if (pfile) {
			fclose(pfile);
			pfile=NULL;
		}
	}

	bool rewind() {

		if (!pfile) return false;

		::rewind(pfile);
		fseek(pfile,CRU_1901_2006_HEADERSIZE+CRU_1901_2006_DATA_LENGTH*CRU_1901_2006_NRECORD,SEEK_CUR);
		recno=0;
		iseof=false;

		return true;
	}

	bool getnext(Cru_1901_2006& obj) {

		if (!pfile || iseof) return false;

		int i;

		getindex(recno);
		getdata();

		obj.lat=popreal(pindex,3,11,1,-900);
		obj.lon=popreal(pindex,3,12,1,-1800);

		for (i=1271;i>=0;i--) obj.msun[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,10,1,-1);
		for (i=1271;i>=0;i--) obj.mprec[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,17,1,-1000);
		for (i=1271;i>=0;i--) obj.mtemp[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,11,1,-1000);
		for (i=0;i>=0;i--) obj.soilcode[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,4,1,0);

		return true;
	}

	bool getindex(Cru_1901_2006& obj) {

		if (!CRU_1901_2006_NRECORD || !pfile) return false;

		// else

		unsigned char ptarget[3]={0,0,0};
		unsigned char buf[4];
		bitify(buf,obj.lon,-1800,1);
		merge(ptarget,buf,12);
		bitify(buf,obj.lat,-900,1);
		merge(ptarget,buf,11);

		int i,c;
		bool not_done=true;
		long direction=1;
		long offset=CRU_1901_2006_NRECORD/2+CRU_1901_2006_NRECORD%2;
		long thisrecord=0-CRU_1901_2006_NRECORD%2;

		while (not_done) {
			thisrecord+=offset*direction;
			if (thisrecord>=CRU_1901_2006_NRECORD) thisrecord=CRU_1901_2006_NRECORD-1;
			else if (thisrecord<0) thisrecord=0;
			if (offset==1) not_done=false;

			getindex(thisrecord);
			getdata();

			c=compare_index(pindex,ptarget);
			if (c<0) direction=1;
			else if (c>0) direction=-1;
			else { // found

				for (i=1271;i>=0;i--) obj.msun[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,10,1,-1);
				for (i=1271;i>=0;i--) obj.mprec[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,17,1,-1000);
				for (i=1271;i>=0;i--) obj.mtemp[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,11,1,-1000);
				for (i=0;i>=0;i--) obj.soilcode[i]=popreal(pdata,CRU_1901_2006_DATA_LENGTH,4,1,0);

				return true;
			}
			offset=offset/2+offset%2;
		}

		return false;
	}
};
