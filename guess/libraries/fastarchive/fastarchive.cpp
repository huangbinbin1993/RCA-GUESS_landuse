///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    FASTARCHIVE                                    //
//                                                                                   //
//                                Written by Ben Smith                               //
//                                 University of Lund                                //
//                                                                                   //
//  Provides functionality for the storage of large numerical data sets in           //
//  platform-independent compressed binary files. Also produces a C++ header file    //
//  containing code for the rapid retrieval of data from the archive via a fast      //
//  iterative search algorithm. Data to be stored must consist of numerical data     //
//  only. Each record must be identifiable by a unique index value or combination    //
//  of values.                                                                       //
//                                                                                   //
//  IMPORTANT - PLEASE NOTE                                                          //
//  Required headers: The following list of #includes for standard C/C++ libraries   //
//    must appear before the #include for the present header file in any file using  //
//    functionality from GUTIL                                                       //
//                                                                                   //
//    #include <stdio.h>                                                             //
//    #include <stdlib.h>                                                            //
//                                                                                   //
//  This version dated 22 July 2006                                                  //
//                                                                                   //
//  Enquiries to: Ben Smith, Lund University: ben.smith@nateko.lu.se                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "fastarchive.h"
#include <string.h>
#include <time.h>

void fail_fastarchive() {

	::printf("Error in FastArchive library: out of memory\n");
	fprintf(stderr,"Error in FastArchive library: out of memory\n");
	exit(99);
}


bool CFastArchive::alphanumeric(const char* name) {

	if (name[0]>='A' && name[0]<='Z' || name[0]>='a' && name[0]<='z') return true;
	return false;
}

Element& CFastArchive::finditem(ItemList& list,const char*& name) {

	unsigned int i;
	for (i=0;i<list.nobj;i++) {
		if (!strcmp(list[i].name,name)) return list[i];
	}

	printf("CFastArchive::finditem: attempt to reference undeclared item \"%s\"\n",(char*)name);
	exit(99);

	return list[0];
}

void CFastArchive::newarchive(const char* name) {

	// Creates a new archive
	// Data in the archive will be saved in <name>.bin
	// Code for retrieving data by index from the archive in <name>.cpp and <name>.h

	char filename[MAXFILE];
	int c,i;
	bool found;

	if (strlen(name)>MAXFILE-5) {
		printf("CFastArchive::newarchive: name \"%s\" longer than permitted length of %d characters\n",
			name,MAXFILE-5);
		exit(99);		
	}

	// Extract directory part (if any)
	
	strcpy(dir,name);
	c=strlen(dir);
	i=c-1;
	found=false;

	while (i>=0) {
		if (dir[i]=='/' || dir[i]=='\\') {
			strcpy(filebase,name+i+1);
			dir[i+1]='\0';
			i=-1;
			found=true;
		}
		i--;
	}

	if (!found) {
		dir[0]='\0';
		strcpy(filebase,name);
	}

	if (strlen(filebase)==0) {
		printf("CFastArchive::newarchive: \"%s\" does not include a file/class name part\n",name);
		exit(99);
	}
	if (strlen(filebase)>MAXNAME) {
		printf("CFastArchive::newarchive: name \"%s\" longer than permitted length of %d characters\n",
			filebase,MAXNAME);
		exit(99);		
	}
	if (!alphanumeric(filebase)) {
		printf("CFastArchive::newarchive: name \"%s\" must begin with an alphabetic character\n",
			name);
		exit(99);		
	}
	if (pfile) closearchive();
	indexlist.killall();
	itemlist.killall();
	datalist.killall();
	strcpy(filename,dir);
	strcat(filename,filebase);
	strcat(filename,".bin");
	nrecord=0;
	pfile=fopen(filename,"wb");
	if (!pfile) {
		printf("CFastArchive::newarchive: could not open %s for output\n",(char*)filename);
		exit(99);
	}
	writing=false;
}

short CFastArchive::defineindex(const char* name,double minval,double maxval,double precision) {

	// Adds an index item to the list for the archive, returning its reference number
	// (always 0 for the first, 1 for the second etc)

	if (strlen(name)>MAXNAME) {
		printf("CFastArchive::defineindex: name \"%s\" longer than permitted length of %d characters\n",
			name,MAXNAME);
		exit(99);
	}
	if (!alphanumeric(name)) {
		printf("CFastArchive::defineindex: index name \"%s\" must begin with an alphabetic character\n",
			name);
		exit(99);		
	}
	if (writing) {
		printf("CFastArchive::defineindex: \"%s\" illegal define after writing commenced\n",name);
		exit(99);
	}
	Element& e=indexlist.createobj();
	e.name=name;
	long nval=(long)((maxval-minval)/precision+1); // number of possible values
	int i=0;
	long v=1;
	while (v<nval) {
		v*=2;
		i++;
	}
	if (i>32) {
		printf("CFastArchive::defineindex: \"%s\" requires %d bits at specified precision\nMaximum allowed is 32 bits\n",
			name,i);
		exit(99);
	}
	e.nbit=i;
	e.nval=1;
	e.scalar=precision;
	e.offset=minval;
	e.maxval=maxval;
	e.firstbit=indexlist.nbit;
	e.reserve(1);
	indexlist.resize();
	return e.id;
}

short CFastArchive::defineitems(const char* name,int nitem,double minval,double maxval,double precision) {

	// Adds a data item to the list for the archive, returning its reference number
	// (always 0 for the first, 1 for the second etc)

	if (strlen(name)>MAXNAME) {
		printf("CFastArchive::defineitems: name \"%s\" longer than permitted length of %d characters\n",
			name,MAXNAME);
		exit(99);
	}
	if (!alphanumeric(name)) {
		printf("CFastArchive::defineitems: item name \"%s\" must begin with an alphabetic character\n",
			name);
		exit(99);
	}
	if (writing) {
		printf("CFastArchive::defineitems: \"%s\" illegal define after writing commenced\n",name);
		exit(99);
	}
	Element& e=itemlist.createobj();
	e.name=name;
	long nval=(long)((maxval-minval)/precision+3);
		// number of possible values
		// 0=minval; nval-3=maxval; nval-2=underflow; nval-1=overflow
	int i=0;
	long v=1;
	while (v<nval) {
		v*=2;
		i++;
	}
	if (i>32) {
		printf("CFastArchive::defineitems: \"%s\" requires %d bits at specified precision\nMaximum allowed is 32 bits\n",
			name,i);
		exit(99);
	}
	e.nbit=i;
	e.nval=nitem;
	e.scalar=precision;
	e.offset=minval;
	e.maxval=maxval;
	e.firstbit=itemlist.nbit;
//	printf("Here#1\n");
	e.reserve(nitem);
//	printf("Here#2\n");
	itemlist.resize();
//	printf("Here#3\n");
	return e.id;
}

void CFastArchive::setindex(const char* name,double val) {

	Element& e=finditem(indexlist,name);
	if (val<e.offset || val>e.maxval) {
		printf("CFastArchive::setindex: %g out of range for index \"%s\"\n",val,name);
		exit(99);
	}
	e.pbuf[0].bitify(val,e.offset,e.scalar);
}

void CFastArchive::setitem(const char* name,double* val) {

	int i;
	Element& e=finditem(itemlist,name);
	for (i=0;i<e.nval;i++) {
		if (val[i]<e.offset || val[i]>e.maxval) {
			printf("CFastArchive::setitem: %g out of range for item \"%s\"\n",val[i],name);
			exit(99);
		}
		e.pbuf[i].bitify(val[i],e.offset,e.scalar);
	}
}

int CFastArchive::compare(unsigned char* a,unsigned char* b) {

	// Compares data item of two DataIndex objects
	// Returns:  0  if a=b
	//           -1 if a<b
	//           +1 if a>b

	int i;
	for (i=0;i<indexlist.nbyte;i++) {
		if (a[i]<b[i]) return -1;
		else if (a[i]>b[i]) return +1;
	}

	return 0;
}

void CFastArchive::storerecord() {

	if (!pfile) {
		printf("CFastArchive::storerecord: illegal call before archive opened for writing");
		exit(99);
	}

	if (!writing) {
		writeheader();
		writing=true;
	}

	indexlist.bitmap();
	itemlist.bitmap();

	long nrecord = static_cast<long>(datalist.nobj);
	long begin = 0;
	long end = nrecord;

	while (begin < end) {
		long middle = (begin+end)/2;

		unsigned char* thisindex = datalist[middle].pindex;

		int c = compare(indexlist.pdata, thisindex);

		if (c > 0) {
			begin = middle+1;
		}
		else if (c < 0) {
			end = middle;
		}
		else {
			printf("CFastArchive::storerecord: duplicate index");
			exit(99);
		}
	}

	DataIndex& ind=datalist.insertobj(begin);
	ind.initialise(indexlist.nbyte,indexlist.pdata,nrecord);

	fwrite(itemlist.pdata,itemlist.nbyte,1,pfile);
}

void CFastArchive::closearchive() {

	if (!pfile) {
		printf("CFastArchive::closearchive: illegal call before archive opened for writing");
		exit(99);
	}

	if (!writing) {
		writeheader();
	}

	datalist.firstobj();
	while (datalist.isobj) {
		DataIndex& ind=datalist.getobj();

		fwrite(ind.pindex,ind.nbyte,1,pfile);
		datalist.nextobj();
	}
	fclose(pfile);

	makesource();

	pfile=NULL;
	writing=false;
	delete pheader;
}

void CFastArchive::writebin(long val,int nbyte) {

	// Writes value to file and header buffer (pheader)
	// as specified number of bytes (max 4), high order first

	unsigned char buf[4];
	long mult[4]={0x1000000,0x10000,0x100,1};
	int i;

	if (!pfile) return;

	for (i=0;i<nbyte;i++) {
		buf[i]=(unsigned char)(val/mult[4-nbyte+i]);
		val-=buf[i]*mult[4-nbyte+i];
		*hptr++=buf[i];
	}

	fwrite(buf,nbyte,1,pfile);
}

void CFastArchive::writenumstr(double dval) {

	// Writes value to file as a string

	char numstr[NUMSTRING];
	int i;

	if (!pfile) return;

	sprintf(numstr,"%g",dval);

	for (i=0;i<NUMSTRING;i++) {
		*hptr++=numstr[i];
	}

	fwrite(numstr,NUMSTRING,1,pfile);
}

void CFastArchive::writeheader() {

	// Writes header to binary file
	// Current format:
	// <byte>     version key
	// <2 bytes>  total length of header in bytes
	// <4 bytes>  number of bytes in index record (including record number, 4 bytes)
	// <2 bytes>  number of index items
	//    (for each index item:)
	//    <byte>     length of item name in bytes (including trailing \0)
	//    <bytes>    item name as null-terminated string
	//    <32 bytes> offset (as string)
	//    <32 bytes> maxval (as string)
	//    <32 bytes> scalar (as string)
	//    <byte>     number of bits
	// <4 bytes>  number of bytes in data record
	// <2 bytes>  number of data items
	//    (for each data item:)
	//    <byte>     length of item name in bytes (including trailing \0)
	//    <bytes>    item name as null-terminated string
	//    <32 bytes> offset (as string)
	//    <32 bytes> maxval (as string)
	//    <32 bytes> scalar (as string)
	//    <byte>     number of bits in each value
	//    <4 bytes>  number of values
	// <4 bytes>  number of records in archive

	if (!pfile) return;

	long n,i;

	// Calculate header size (see above)

	hsize=1+2+4+2;

	indexlist.firstobj();
	while (indexlist.isobj) {
		Element& e=indexlist.getobj();
		hsize+=strlen(e.name)+1;
		indexlist.nextobj();
	}

	hsize+=(1+32+32+32+1)*indexlist.nobj;

	itemlist.firstobj();
	while (itemlist.isobj) {
		Element& e=itemlist.getobj();
		hsize+=strlen(e.name)+1;
		itemlist.nextobj();
	}

	hsize+=4+2+(1+32+32+32+1+4)*itemlist.nobj+4;

	pheader=new unsigned char[hsize];
	if (!pheader) {
		printf("CFastArchive::writeheader: out of memory\n");
		exit(99);
	}
	hptr=pheader;

	writebin(VERSION_KEY,1);
	writebin(hsize,2);
	writebin(indexlist.nbyte+4,4);
	writebin(indexlist.nobj,2);

	indexlist.firstobj();
	while (indexlist.isobj) {
		Element& e=indexlist.getobj();
		n=strlen(e.name)+1;
		writebin(n,1);
		for (i=0;i<n;i++) *hptr++=e.name[i];
		fwrite(e.name,n,1,pfile);
		writenumstr(e.offset);
		writenumstr(e.maxval);
		writenumstr(e.scalar);
		writebin(e.nbit,1);
		indexlist.nextobj();
	}

	writebin(itemlist.nbyte,4);
	writebin(itemlist.nobj,2);

	itemlist.firstobj();
	while (itemlist.isobj) {
		Element& e=itemlist.getobj();
		n=strlen(e.name)+1;
		writebin(n,1);
		for (i=0;i<n;i++) *hptr++=e.name[i];
		fwrite(e.name,n,1,pfile);
		writenumstr(e.offset);
		writenumstr(e.maxval);
		writenumstr(e.scalar);
		writebin(e.nbit,1);
		writebin(e.nval,4);
		itemlist.nextobj();
	}
	
	writebin(datalist.nobj,4);
}

void CFastArchive::makesource() {

	char filename[MAXFILE];
	char classname[MAXNAME+1];
	char allcaps[MAXNAME+1];
	int i,j;

	strcpy(filename,dir);
	strcat(filename,filebase);
	strcat(filename,".h");

	FILE* out=fopen(filename,"wt");
	if (!out) {
		printf("CFastArchive::makesource: could not open %s for output\n",(char*)filename);
		exit(99);
	}

	strcpy(classname,filebase);
	if (classname[0]>='a' && classname[0]<='z') classname[0]-='a'-'A';
	strcpy(allcaps,filebase);
	for (i=0;i<strlen(allcaps);i++)
		if (allcaps[i]>='a' && allcaps[i]<='z') allcaps[i]-='a'-'A';

	time_t t;
	time(&t);

	fprintf(out,"//////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(out,"// %s.H\n",allcaps);
	fprintf(out,"// Header file for input from a fast data archive\n");
	fprintf(out,"// Created automatically by FastArchive on %s//\n",ctime(&t));
	fprintf(out,"// The following #includes should appear in your source code file:\n//\n");
	fprintf(out,"//   #include <stdio.h>\n");
	fprintf(out,"//   #include <stdlib.h>\n");
	fprintf(out,"//   #include <string.h>\n");
	fprintf(out,"//   #include \"%s\"\n",filename);
	fprintf(out,"//\n");
	fprintf(out,"// Functionality to retrieve data from the archive is provided by class %sArchive.\n",classname);
	fprintf(out,"// The following public functions are provided:\n");
	fprintf(out,"//\n");
	fprintf(out,"// bool open(char* filename)\n");
	fprintf(out,"//   Attempts to open the specified file as a fast data archive. The format must be\n");
	fprintf(out,"//   exactly compatible with this version of %s.h (normally the archive and\n",filebase);
	fprintf(out,"//   header file should have been produced together by the same program using class\n");
	fprintf(out,"//   CFastArchive). Returns false if the file could not be opened or had format\n");
	fprintf(out,"//   errors. open() with no argument is equivalent to open(\"%s.bin\").\n",filebase);
	fprintf(out,"//\n");
	fprintf(out,"// void close()\n");
	fprintf(out,"//   Closes the archive (if open).\n");
	fprintf(out,"//\n");
	fprintf(out,"// bool rewind()\n");
	fprintf(out,"//   Sets the file pointer to the first record in the archive file. Returns false if\n");
	fprintf(out,"//   no archive file is currently open.\n");
	fprintf(out,"//\n");
	fprintf(out,"// bool getnext(%s& obj)\n",classname);
	fprintf(out,"//   Retrieves the next record in the archive file and advances the file pointer to\n");
	fprintf(out,"//   the next record. Data are written to the member variables of obj. Returns false if\n");
	fprintf(out,"//   no archive file is currently open or if the file pointer is beyond the last\n");
	fprintf(out,"//   record. Use rewind() and getnext() to retrieve data sequentially from the archive.\n");
	fprintf(out,"//\n");
	fprintf(out,"// bool getindex(%s& obj)\n",classname);
	fprintf(out,"//   Searches the archive for a record matching the values specified for the index\n");
	fprintf(out,"//   item");
	if (indexlist.nobj>1) fprintf(out, "s");
	fprintf(out," (");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,"%s",e.name);
		if (i==indexlist.nobj-2) fprintf(out," and ");
		else if (i<indexlist.nobj-2) fprintf(out,", ");
	}

	fprintf(out,") in obj. If a matching record is found, the data are\n");
	fprintf(out,"//   written to the member variables of obj. Returns true if the archive was open and\n");
	fprintf(out,"//   a matching record was found, otherwise false. The search is iterative and fast.\n");
	fprintf(out,"//\n");
	fprintf(out,"// Sample program:\n");
	fprintf(out,"//\n");
	fprintf(out,"//   %sArchive ark;\n",classname);
	fprintf(out,"//   %s data;\n",classname);
	fprintf(out,"//   bool success,flag;\n");
	fprintf(out,"//\n");
	fprintf(out,"//   // Retrieve all records in sequence and print values of ");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,"%s",e.name);
		if (i==indexlist.nobj-2) fprintf(out," and ");
		else if (i<indexlist.nobj-2) fprintf(out,", ");
	}
	
	fprintf(out,":\n");
	fprintf(out,"//\n");
	fprintf(out,"//   success=ark.open(\"%s.bin\");\n",filebase);
	fprintf(out,"//   if (success) {\n");
	fprintf(out,"//      flag=ark.rewind();\n");
	fprintf(out,"//      while (flag) {\n");
	fprintf(out,"//         flag=ark.getnext(data);\n");
	fprintf(out,"//         if (flag)\n");
	fprintf(out,"//            printf(\"Loaded record: ");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,"%s=%%g",e.name);
		if (i<indexlist.nobj-1) fprintf(out,", ");
	}

	fprintf(out,"\\n\"");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,",data.%s",e.name);
	}

	fprintf(out,");\n");
	fprintf(out,"//      }\n");
	fprintf(out,"//   }\n");
	fprintf(out,"//   \n");
	fprintf(out,"//   // Look for a record with ");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,"%s=%g",e.name,e.offset);
		if (i<indexlist.nobj-1) fprintf(out,", ");
	}

	fprintf(out,":\n");
	fprintf(out,"//\n");

	for (i=0;i<indexlist.nobj;i++) {
		Element& e=indexlist[i];
		fprintf(out,"//   data.%s=%g;\n",e.name,e.offset);
	}
	
	fprintf(out,"//   success=ark.getindex(data);\n");
	fprintf(out,"//   if (success) printf(\"Found it!\\n\");\n");
	fprintf(out,"//   else printf(\"Not found\\n\");\n");
	fprintf(out,"//\n");
	fprintf(out,"//   ark.close();\n\n\n");

	fprintf(out,"struct %s {\n\n",classname);
	fprintf(out,"	// Index part\n\n");

	indexlist.firstobj();
	while (indexlist.isobj) {
		Element& e=indexlist.getobj();
		fprintf(out,"	double %s;\n",e.name);
		indexlist.nextobj();
	}

	fprintf(out,"\n	// Data part\n\n");

	itemlist.firstobj();
	while (itemlist.isobj) {
		Element& e=itemlist.getobj();
		fprintf(out,"	double %s[%d];\n",e.name,e.nval);
		itemlist.nextobj();
	}

	fprintf(out,"};\n\n\n");

	fprintf(out,"const long %s_NRECORD=%d;\n",allcaps,datalist.nobj);
	fprintf(out,"const int %s_DATA_LENGTH=%d;\n",allcaps,itemlist.nbyte);
	fprintf(out,"const int %s_INDEX_LENGTH=%d;\n",allcaps,indexlist.nbyte+4);
	fprintf(out,"const int %s_HEADERSIZE=%ld;\n",allcaps,hsize);
	fprintf(out,"const unsigned char %s_HEADER[%s_HEADERSIZE-4]={",allcaps,allcaps);

	for (i=0;i<hsize-4;i++) {
		if (!(i%20)) fprintf(out,"\n\t");
		fprintf(out,"0x%02X",(unsigned int)pheader[i]);
		if (i<hsize-5) fprintf(out,",");
	}
	fprintf(out,"};\n\n\n");

	fprintf(out,"class %sArchive {\n\n",classname);
	fprintf(out,"private:\n\n");
	fprintf(out,"	FILE* pfile;\n	long recno;\n	long datano;\n");
	fprintf(out,"	unsigned char pindex[%s_INDEX_LENGTH];\n",allcaps);
	fprintf(out,"	unsigned char pdata[%s_DATA_LENGTH];\n",allcaps);
	fprintf(out,"	bool iseof;\n\n");

	fprintf(out,"	long readbin(int nbyte) {\n\n");

	fprintf(out,"		unsigned char buf[4];\n");
	fprintf(out,"		long mult[4]={0x1000000,0x10000,0x100,1},val;\n");
	fprintf(out,"		int i;\n\n");

	fprintf(out,"		fread(buf,nbyte,1,pfile);\n\n");

	fprintf(out,"		val=0;\n");
	fprintf(out,"		for (i=0;i<nbyte;i++) {\n");
	fprintf(out,"			val+=buf[i]*mult[4-nbyte+i];\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		return val;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	void getindex(long n) {\n\n");
	fprintf(out,"		fseek(pfile,%s_INDEX_LENGTH*(n-recno),SEEK_CUR);\n",allcaps);
	fprintf(out,"		fread(pindex,%s_INDEX_LENGTH,1,pfile);\n",allcaps);
	fprintf(out,"		datano=pindex[%s_INDEX_LENGTH-4]*0x1000000+pindex[%s_INDEX_LENGTH-3]*0x10000+\n",allcaps,allcaps);
	fprintf(out,"			pindex[%s_INDEX_LENGTH-2]*0x100+pindex[%s_INDEX_LENGTH-1];\n",allcaps,allcaps);
	fprintf(out,"		recno=n+1;\n");
	fprintf(out,"		iseof=(recno==%s_NRECORD);\n	}\n\n",allcaps);

	fprintf(out,"	void getdata() {\n\n");
			
	fprintf(out,"		fseek(pfile,%s_INDEX_LENGTH*-recno+(datano-%s_NRECORD)*%s_DATA_LENGTH,SEEK_CUR);\n",allcaps,allcaps,allcaps);
	fprintf(out,"		fread(pdata,%s_DATA_LENGTH,1,pfile);\n",allcaps);
	fprintf(out,"		fseek(pfile,%s_DATA_LENGTH*(%s_NRECORD-datano-1)+%s_INDEX_LENGTH*recno,SEEK_CUR);\n",allcaps,allcaps,allcaps);
	fprintf(out,"	}\n\n");

	fprintf(out,"	double popreal(unsigned char* bits,int nbyte,int nbit,double scalar,double offset) {\n\n");
		
	fprintf(out,"		unsigned char buf;\n");
	fprintf(out,"		int nb=nbit/8,i;\n");
	fprintf(out,"		double rval=0.0;\n");
	fprintf(out,"		long mult[4]={1,0x100,0x10000,0x1000000};\n\n");

	fprintf(out,"		for (i=0;i<4;i++) {\n");
	fprintf(out,"			if (i<nb) rval+=bits[nbyte-i-1]*mult[i];\n");
	fprintf(out,"			else if (i==nb) {\n");
	fprintf(out,"				buf=bits[nbyte-i-1]<<(8-nbit%%8);\n");
	fprintf(out,"				buf>>=8-nbit%%8;\n");
	fprintf(out,"				rval+=buf*mult[i];\n");
	fprintf(out,"			}\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		for (i=nbyte-1;i>=0;i--) {\n");
	fprintf(out,"			if (i>=nb)\n");
	fprintf(out,"				bits[i]=bits[i-nb];\n");
	fprintf(out,"			else\n");
	fprintf(out,"				bits[i]=0;\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		nb=nbit%%8;\n\n");

	fprintf(out,"		for (i=nbyte-1;i>=0;i--) {\n");
	fprintf(out,"			bits[i]>>=nb;\n");
	fprintf(out,"			if (i>0) {\n");
	fprintf(out,"				buf=bits[i-1];\n");
	fprintf(out,"				buf<<=8-nb;\n");
	fprintf(out,"				bits[i]|=buf;\n");
	fprintf(out,"			}\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		rval=rval*scalar+offset;\n\n");

	fprintf(out,"		return rval;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	void bitify(unsigned char buf[4],double fval,double offset,double scalar) {\n\n");

	fprintf(out,"		long ival = (long)((fval-offset)/scalar + 0.5);\n");
	fprintf(out,"		buf[0]=(unsigned char)(ival/0x1000000);\n");
	fprintf(out,"		ival-=buf[0]*0x1000000;\n");
	fprintf(out,"		buf[1]=(unsigned char)(ival/0x10000);\n");
	fprintf(out,"		ival-=buf[1]*0x10000;\n");
	fprintf(out,"		buf[2]=(unsigned char)(ival/0x100);\n");
	fprintf(out,"		ival-=buf[2]*0x100;\n");
	fprintf(out,"		buf[3]=(unsigned char)(ival);\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	void merge(unsigned char ptarget[%d],unsigned char buf[4],int bits) {\n\n",indexlist.nbyte);

	fprintf(out,"		int nb=bits/8;\n");
	fprintf(out,"		int i;\n");
	fprintf(out,"		unsigned char nib;\n");
	fprintf(out,"		for (i=0;i<%d;i++) {\n\n",indexlist.nbyte);
	fprintf(out,"			if (i<%d-nb)\n",indexlist.nbyte);
	fprintf(out,"				ptarget[i]=ptarget[i+nb];\n");
	fprintf(out,"			else\n");
	fprintf(out,"				ptarget[i]=0;\n");
	fprintf(out,"		}\n");
	fprintf(out,"		nb=bits%%8;\n");
	fprintf(out,"		for (i=0;i<%d;i++) {\n",indexlist.nbyte);
	fprintf(out,"			ptarget[i]<<=nb;\n");
	fprintf(out,"			if (i<%d-1) {\n",indexlist.nbyte);
	fprintf(out,"				nib=ptarget[i+1]>>(8-nb);\n");
	fprintf(out,"				ptarget[i]|=nib;\n");
	fprintf(out,"			}\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		nb=bits/8;\n");
	fprintf(out,"		if (bits%%8) nb++;\n");
	fprintf(out,"		for (i=1;i<=nb;i++)\n");
	fprintf(out,"			ptarget[%d-i]|=buf[4-i];\n",indexlist.nbyte);
	fprintf(out,"	}\n\n");

	fprintf(out,"	int compare_index(unsigned char* a,unsigned char* b) {\n\n");

	fprintf(out,"		int i;\n");
	fprintf(out,"		for (i=0;i<%d;i++) {\n",indexlist.nbyte);
	fprintf(out,"			if (a[i]<b[i]) return -1;\n");
	fprintf(out,"			else if (a[i]>b[i]) return +1;\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		return 0;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	bool initialise(const char* filename) {\n\n");

	fprintf(out,"		int i;\n");
	fprintf(out,"		unsigned char* pheader;\n\n");
	
	fprintf(out,"		if (pfile) fclose(pfile);\n");
	fprintf(out,"		pfile=fopen(filename,\"rb\");\n");
	fprintf(out,"		if (!pfile) {\n");
	fprintf(out,"			printf(\"Could not open %%s for input\\n\",filename);\n");
	fprintf(out,"			return false;\n");
	fprintf(out,"		}\n\n");
	
	fprintf(out,"		pheader=new unsigned char[%s_HEADERSIZE-4];\n",allcaps);
	fprintf(out,"		if (!pheader) {\n");
	fprintf(out,"			printf(\"Out of memory\\n\");\n");
	fprintf(out,"			fclose(pfile);\n");
	fprintf(out,"			pfile=NULL;\n");
	fprintf(out,"			return false;\n");
	fprintf(out,"		}\n");
	fprintf(out,"		::rewind(pfile);\n");
	fprintf(out,"		fread(pheader,%s_HEADERSIZE-4,1,pfile);\n",allcaps);
	fprintf(out,"		for (i=0;i<%s_HEADERSIZE-4;i++) {\n",allcaps);
	fprintf(out,"			if (pheader[i]!=%s_HEADER[i]) {\n",allcaps);
	fprintf(out,"				printf(\"Format of %%s incompatible with this version of %s.h\\n\",filename);\n",filebase);
	fprintf(out,"				fclose(pfile);\n");
	fprintf(out,"				pfile=NULL;\n");
	fprintf(out,"				delete pheader;\n");
	fprintf(out,"				return false;\n");
	fprintf(out,"			}\n");
	fprintf(out,"		}\n");
	fprintf(out,"		delete[] pheader;\n\n");

	fprintf(out,"		::rewind(pfile);\n");
	fprintf(out,"		fseek(pfile,%s_HEADERSIZE+%s_DATA_LENGTH*%s_NRECORD,SEEK_CUR);\n",allcaps,allcaps,allcaps);
	fprintf(out,"		recno=0;\n");
	fprintf(out,"		iseof=false;\n\n");
	fprintf(out,"		return true;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"public:\n\n");

	fprintf(out,"	%sArchive() {\n",classname);
	fprintf(out,"		pfile=NULL;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	~%sArchive() {\n",classname);
	fprintf(out,"		if (pfile) fclose(pfile);\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	bool open(const char* filename) {\n");
	fprintf(out,"		return initialise(filename);\n");
	fprintf(out,"	}\n\n");
	fprintf(out,"	bool open() {\n");
	fprintf(out,"		return open(\"%s.bin\");\n",filebase);
	fprintf(out,"	}\n\n");

	fprintf(out,"	void close() {\n");
	fprintf(out,"		if (pfile) {\n");
	fprintf(out,"			fclose(pfile);\n");
	fprintf(out,"			pfile=NULL;\n");
	fprintf(out,"		}\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	bool rewind() {\n\n");

	fprintf(out,"		if (!pfile) return false;\n\n");

	fprintf(out,"		::rewind(pfile);\n");
	fprintf(out,"		fseek(pfile,%s_HEADERSIZE+%s_DATA_LENGTH*%s_NRECORD,SEEK_CUR);\n",allcaps,allcaps,allcaps);
	fprintf(out,"		recno=0;\n");
	fprintf(out,"		iseof=false;\n\n");

	fprintf(out,"		return true;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	bool getnext(%s& obj) {\n\n",classname);

	fprintf(out,"		if (!pfile || iseof) return false;\n\n");

	fprintf(out,"		int i;\n\n");

	fprintf(out,"		getindex(recno);\n");
	fprintf(out,"		getdata();\n\n");

	i=indexlist.nobj;
	for (j=i-1;j>=0;j--) {
		Element& e=indexlist[j];
		fprintf(out,"		obj.%s=popreal(pindex,%d,%d,%g,%g);\n",
			e.name,indexlist.nbyte,e.nbit,e.scalar,e.offset);
	}

	fprintf(out,"\n");

	i=itemlist.nobj;
	for (j=i-1;j>=0;j--) {
		Element& e=itemlist[j];
		fprintf(out,"		for (i=%d;i>=0;i--) obj.%s[i]=popreal(pdata,%s_DATA_LENGTH,%d,%g,%g);\n",
			e.nval-1,e.name,allcaps,e.nbit,e.scalar,e.offset);
	}

	fprintf(out,"\n");

	fprintf(out,"		return true;\n");
	fprintf(out,"	}\n\n");

	fprintf(out,"	bool getindex(%s& obj) {\n\n",classname);

	fprintf(out,"		if (!%s_NRECORD || !pfile) return false;\n\n",allcaps);

	fprintf(out,"		// else\n\n");

	fprintf(out,"		unsigned char ptarget[%d]={",indexlist.nbyte);

	for (i=0;i<indexlist.nbyte;i++) {
		fprintf(out,"0");
		if (i<indexlist.nbyte-1) fprintf(out,",");
	}

	fprintf(out,"};\n");

	fprintf(out,"		unsigned char buf[4];\n");

	indexlist.firstobj();
	while (indexlist.isobj) {
		Element& e=indexlist.getobj();
		fprintf(out,"		bitify(buf,obj.%s,%g,%g);\n",e.name,e.offset,e.scalar);
		fprintf(out,"		merge(ptarget,buf,%d);\n",e.nbit);
		indexlist.nextobj();
	}

	fprintf(out,"\n");

	fprintf(out,"		long begin = 0;\n");
	fprintf(out,"		long end = %s_NRECORD;\n\n",allcaps);

	fprintf(out,"		while (begin < end) {\n");
	fprintf(out,"			long middle = (begin+end)/2;\n\n");

	fprintf(out,"			getindex(middle);\n");
	fprintf(out,"			getdata();\n\n");

	fprintf(out,"			int c = compare_index(pindex, ptarget);\n\n");

	fprintf(out,"			if (c < 0) {\n");
	fprintf(out,"				begin = middle + 1;\n");
	fprintf(out,"			}\n");
	fprintf(out,"			else if (c > 0) {\n");
	fprintf(out,"				end = middle;\n");
	fprintf(out,"			}\n");
	fprintf(out,"			else {\n\n");

	i=itemlist.nobj;
	for (j=i-1;j>=0;j--) {
		Element& e=itemlist[j];
		fprintf(out,"				for (int i=%d;i>=0;i--) obj.%s[i]=popreal(pdata,%s_DATA_LENGTH,%d,%g,%g);\n",
			e.nval-1,e.name,allcaps,e.nbit,e.scalar,e.offset);
	}

	fprintf(out,"\n");
					
	fprintf(out,"				return true;\n");
	fprintf(out,"			}\n");
	fprintf(out,"		}\n\n");

	fprintf(out,"		return false;\n");
	fprintf(out,"	}\n");
	fprintf(out,"};\n");

	fclose(out);
}


