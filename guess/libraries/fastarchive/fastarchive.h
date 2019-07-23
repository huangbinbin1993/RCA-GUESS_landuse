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


///////////////////////////////////////////////////////////////////////////////////////
// Class CFastArchive
// Provides functionality for WRITING to a fast archive.
// Class for reading from a fast archive is generated in conjunction with the archive
// itself as a C++ header file. Documentation is provided in the generated header file.
//
// PUBLIC FUNCTIONS
// 
// void newarchive(char* name);
//   Create a new archive file with specified name. A directory part may be included.
//   The specified name should not include an extension but .bin will be added to the
//   name of the archive file created, and .h to the name of the customised C++ header
//   file for the archive.
//
// void closearchive();
//   Close the archive file opened with newarchive(). This function MUST be called,
//   otherwise the archive file will not be readable.
//
// short defineindex(char* name,double minval,double maxval,double precision);
//   Add a new index item (i.e. variable) to the data structure of the archive. Values
//   will (and must) be in the range minval-maxval and will be stored at the specified
//   precision. An index item is used in the application program as part of the key
//   to retrieve a specified record from the archive.
//
// short defineitems(char* name,int nitem,double minval,double maxval,double precision)
//   Add a new item to the data structure of the archive. Values will (and must) be in
//   the range minval-maxval and will be stored at the specified precision.
//
// void setindex(char* name,double val)
//   Set the value for the specified index item in the current record. The index item
//   must have been previously declared by a call to defineindex().
//
// void setitem(char* name,double* val)
//   Set the value for the specified item in the current record. The item must have been
//   previously declared by a call to defineitems().
//
// void storerecord();
//   Store the current record in the archive file. The values of all the index item(s)
//   should have been previously set with setindex(). The values of all the item(s)
//   should have been previously set with setitem().
	
const int MAXFILE=1024;
const int MAXNAME=64;
const int NUMSTRING=32;
const int SEGSIZE_FASTARCHIVE=16;
const unsigned char VERSION_KEY=1; // Change for updated versions

// Flags
const int FA_ALLOW_OVERFLOW=0x0001;
const int FA_ALLOW_DUPLICATES=0x0002;
const int FA_VERBOSE=0x0004;

void fail_fastarchive();

///////////////////////////////////////////////////////////////////////////////////////
// LISTARRAY_ID
// List Array of objects with id member and no reference members
// Adapted from version in GUTIL library

template<class tdata> class ListArray_id_fastarchive {

	class Item {

	public:
		Item* pnext;
		Item* pprev;
		tdata object;
	};


private:
	Item** array;
	Item* pfirstitem;
	Item* plastitem;
	Item* pthisitem;
	unsigned int nseg;
	unsigned int id;
	unsigned int thisobj;


public:
	bool isobj;
	unsigned int nobj;

public:
	ListArray_id_fastarchive() {
		pfirstitem=NULL;
		plastitem=NULL;
		pthisitem=NULL;
		isobj=false;
		nobj=0;
		nseg=0;
		id=0;
	}

	void killall() {
		Item* pitem=pfirstitem;
		Item* pnext;
		while (pitem) {
			pnext=pitem->pnext;
			delete pitem;
			pitem=pnext;
		}
		if (nseg) delete[] array;
		pfirstitem=NULL;
		plastitem=NULL;
		pthisitem=NULL;
		isobj=false;
		nobj=0;
		nseg=0;
		id=0;
	}

	~ListArray_id_fastarchive() {
		killall();
	}

	void initarray(unsigned int nitem) {
		unsigned int i;
		killall();
		if (nitem<=0) return;
		nseg=nitem%SEGSIZE_FASTARCHIVE;
		if (nseg*SEGSIZE_FASTARCHIVE<nitem) nseg++;
		array=new Item*[nseg*SEGSIZE_FASTARCHIVE];
		if (!array) fail_fastarchive();
		for (i=0;i<nitem;i++) createobj();
		pthisitem=0;
		thisobj=0;
	}

	tdata& createobj() {
		Item* pitem=new Item;
		if (!pitem) fail_fastarchive();
		pitem->object.id=id++;
		pitem->pprev=plastitem;
		pitem->pnext=NULL;
		if (plastitem) plastitem->pnext=pitem;
		else pfirstitem=pitem;
		plastitem=pitem;
		pthisitem=pitem;
		isobj=true;
		thisobj=nobj++;
		if (nobj>nseg*SEGSIZE_FASTARCHIVE) {
			Item** newarray=new Item*[(nseg+1)*SEGSIZE_FASTARCHIVE];
			if (!newarray) fail_fastarchive();
			if (nseg++) {
				unsigned int i;
				for (i=0;i<nobj;i++)
					newarray[i]=array[i];
				delete[] array;
			}
			array=newarray;
		}
		array[nobj-1]=pitem;
		return pitem->object;
	}

	tdata& insertobj(unsigned int i) {
		unsigned int j;
		if (i>=id) {
			return createobj();
		}
		Item* pitem=new Item;
		if (!pitem) fail_fastarchive();
		Item* polditem=array[i];
		pitem->object.id=i;
		pitem->pprev=polditem->pprev;
		pitem->pnext=polditem;
		polditem->pprev=pitem;
		if (pitem->pprev) pitem->pprev->pnext=pitem;
		if (!i) pfirstitem=pitem;
		id++;
		nobj++;
		if (nobj>nseg*SEGSIZE_FASTARCHIVE) {
			Item** newarray=new Item*[(nseg+1)*SEGSIZE_FASTARCHIVE];
			if (!newarray) fail_fastarchive();
			if (nseg++) {
				for (j=0;j<nobj;j++)
					newarray[j]=array[j];
				delete[] array;
			}
			array=newarray;
		}
		for (j=nobj-1;j>i;j--) {
			array[j]=array[j-1];
			array[j]->object.id++;
		}
		array[i]=pitem;
		pthisitem=pitem;
		isobj=true;
		thisobj=i;
		return pitem->object;
	}
	
	bool firstobj() {
		if (pfirstitem) {
			pthisitem=pfirstitem;
			isobj=true;
			thisobj=0;
			return true;
		}
		return false;
	}

	bool nextobj() {
		if (pthisitem) {
			pthisitem=pthisitem->pnext;
			thisobj++;
		}
		if (pthisitem) return true;
		isobj=false;
		return false;
	}

	tdata& getobj() {
		return pthisitem->object;
	}

	tdata& operator[](unsigned int i) {
		return array[i]->object;
	}

	void killobj() {
		unsigned int i;
		if (!pthisitem) return;
		if (pthisitem==pfirstitem) pfirstitem=pthisitem->pnext;
		if (pthisitem==plastitem) plastitem=pthisitem->pprev;
		Item* pnextitem=pthisitem->pnext;
		if (pthisitem->pprev) pthisitem->pprev->pnext=pthisitem->pnext;
		if (pthisitem->pnext) pthisitem->pnext->pprev=pthisitem->pprev;
		delete pthisitem;
		pthisitem=pnextitem;
		if (!pthisitem) isobj=false;
		for (i=thisobj+1;i<nobj;i++)
			array[i-1]=array[i];
		nobj--;
	}
};

class Bitbuf {

public:
	unsigned char buf[4]; // bit buffer with up to 4*8 bits, highest order first

	void bitify(double fval,double offset,double scalar) {

		// Transforms val to binary and stores in bitbuf

		long ival=(fval-offset)/scalar+0.5;
		buf[0]=ival/0x1000000;
		ival-=buf[0]*0x1000000;
		buf[1]=ival/0x10000;
		ival-=buf[1]*0x10000;
		buf[2]=ival/0x100;
		ival-=buf[2]*0x100;
		buf[3]=ival;
	}
};


class Element {

public:
	unsigned int id;
	const char* name;
	int firstbit; // bit number for first bit in data buffer
	int nbit; // number of bits to store value
	short nval; // number of values in element (1 for single value, >1 for an array)
	double scalar;
	double offset; // Actual value restored by VAL*scalar+offset;
	double maxval;
	Bitbuf* pbuf;

public:
	Element() {

		nval=1;
		pbuf=new Bitbuf[1];
		if (!pbuf) {
			printf("Element::Element: out of memory\n");
			exit(99);
		}
	};

	~Element() {

		delete pbuf;
	};

	void reserve(int nval) {
		
		delete pbuf;
		pbuf=new Bitbuf[nval];
		if (!pbuf) {
			printf("Element::reserve: out of memory\n");
			exit(99);
		}
		int i;
		for (i=0;i<nval;i++)
			pbuf[i].bitify(0,offset,scalar);
	}
};


class ItemList : public ListArray_id_fastarchive<Element> {
	
	// List of index or data items in record (e.g. mtemp, soiltype)

public:
	int nbit;
	int nbyte;
	unsigned char* pdata;

public:
	ItemList() {
	
		nbit=0;
		nbyte=1;
		pdata=new unsigned char[nbyte];
		if (!pdata) {
			printf("ItemList::ItemList: out of memory\n");
			exit(99);
		}
	};

	~ItemList() {

		delete pdata;
	};

	void resize() {

		delete pdata;
		nbit=0;
		firstobj();
		while (isobj) {
			Element& e=getobj();
			nbit+=e.nbit*e.nval;
			nextobj();
		}
		nbyte=nbit/8;
		if (nbit%8) nbyte++;
		printf("Creating new pdata at %d bytes\n",nbyte);
		pdata=new unsigned char[nbyte];
		printf("Done\n");
		if (!pdata) {
			printf("ItemList::resize: out of memory\n");
			exit(99);
		}
		int i;
		for (i=0;i<nbyte;i++)
			pdata[i]='\0';
	};

	void shoveup(int bits) {

		// Shoves up bit sequence left by bits

		int nb=bits/8;
		int i,j;
		unsigned char nib;
		for (i=0;i<nbyte;i++) {
			if (i<nbyte-nb)
				pdata[i]=pdata[i+nb];
			else
				pdata[i]=0;
		}
		nb=bits%8;
		for (i=0;i<nbyte;i++) {
			pdata[i]<<=nb;
			if (i<nbyte-1) {
				nib=pdata[i+1]>>(8-nb);
				pdata[i]|=nib;
			}
		}
	};

	void bitmap() {

		// Transfers data from all elements to bit buffer (pdata)

		int i,j,nb;
		for (i=0;i<nbyte;i++)
			pdata[i]='\0';
		firstobj();
		while (isobj) {
			Element& e=getobj();
			nb=e.nbit/8;
			if (e.nbit%8) nb++;
			for (i=0;i<e.nval;i++) {
				shoveup(e.nbit);
				for (j=1;j<=nb;j++)
					pdata[nbyte-j]|=e.pbuf[i].buf[4-j];
			}
			nextobj();			
		}
	};
};


class DataIndex {

	// Stores index as bit sequence to a single actual data record

public:
	int id;
	unsigned char* pindex;
	int nbyte;

public:
	void initialise(int nb,unsigned char* sequence,long recno) {
		int i;
		nbyte=nb+4; // reserve four bytes for record number
		pindex=new unsigned char[nbyte];
		if (!pindex) {
			printf("DataIndex::initialise: out of memory\n");
			exit(99);
		}
		for (i=0;i<nbyte;i++)
			pindex[i]=sequence[i];
		pindex[nb]=recno/0x1000000;
		recno-=pindex[nb]*0x1000000;
		pindex[nb+1]=recno/0x10000;
		recno-=pindex[nb+1]*0x10000;
		pindex[nb+2]=recno/0x100;
		recno-=pindex[nb+2]*0x100;
		pindex[nb+3]=recno;
	}

	DataIndex() {
		nbyte=1;
		pindex=new unsigned char[nbyte];
		if (!pindex) {
			printf("DataIndex::DataIndex: out of memory\n");
			exit(99);
		}
	};
	
	~DataIndex() {
	
		delete pindex;
	}
};

typedef ListArray_id_fastarchive<DataIndex> DataList;

class CFastArchive {

private:
	char filebase[MAXFILE];
	char dir[MAXFILE];
	ItemList indexlist;
	ItemList itemlist;
	DataList datalist;
	int nrecord;
	FILE* pfile;
	bool writing;
	unsigned char* pheader;
	unsigned char* hptr;
	long hsize;

	Element& finditem(ItemList& list,const char*& name);
	int compare(unsigned char* a,unsigned char* b);
	bool alphanumeric(const char* name);
	void makesource();
	void writebin(long val,int nbyte);
	void writenumstr(double dval);
	void writeheader();

public:
	CFastArchive() {
		pfile=NULL;
	}
	void closearchive();
	~CFastArchive() {
		if (pfile) closearchive();
	}
	void newarchive(const char* name);
	short defineindex(const char* name,double minval,double maxval,double precision);
	short defineitems(const char* name,int nitem,double minval,double maxval,double precision);
	void setindex(const char* name,double val);
	void setitem(const char* name,double* val);
	void storerecord();
};
