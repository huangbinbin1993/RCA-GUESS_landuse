#MPI_SPECIAL := -Nmixrpath
MPI_SPECIAL := -Nmpi
ifeq (MPI_SRC,$(findstring MPI_SRC,$(DEFINITIONS)))
	MPI_SPECIAL += -Nmpi
endif
OPT :=      -O3 -unroll-aggressive              -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 
#OPT :=      -O3 -unroll-aggressive              -ip -ipo  -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 
#OPT :=      -O3 -unroll-aggressive -xCORE-AVX2 -ip -ipo  -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 
#OPT :=      -O3 -unroll-aggressive -xAVX -ip -ipo  -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 
#OPT :=      -O3 -unroll-aggressive -xSSE4.2 -axSSE4.2 -ip -ipo  -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 
#OPT := -pg -O3 -unroll-aggressive -xSSE4.2 -axSSE4.2 -ip -ipo  -ansi-alias -fast-transcendentals -finline -inline-level=2  -scalar-rep -no-prec-div 


#OPT := -O1
#DBG :=
#ADVANCED := -funroll-loops
#CODEGEN :=
#COMPATIBILITY :=
#COMPONENT :=
#DATA :=
#DEPRECATED :=
#DIAGNOSTICS :=
#FLOAT :=
#INLINING :=
#IPO :=
#LANGUAGE :=
#LINK :=
#MISC :=
#OPT :=
#OUTPUT :=
#PGO :=
#PREPROC :=
#OPT_REPORT :=
#OPENMP := 


TRACE := #-lVT -I$VT_ROOT/include -L$VT_LIB_DIR $VT_ADD_LIBS #assumes module add itac

CC      := icc 
CCFLAGS := $(OPT) $(DBG) $(DEFINITIONS) $(MPI_SPECIAL) $(TRACE) -strict-ansi

CXX     := icc 
CXXFLAGS := $(OPT) $(DBG) $(DEFINITIONS) $(MPI_SPECIAL) $(TRACE) -strict-ansi

FC := ifort 


FDBG := $(DBG) -traceback #-CB -warn errors -implicitnone -debug full -traceback -check bounds -CU -warn declarations -check all -fp-stack-check 

FHACKS := -heap-arrays 1 -fp-model precise -fp-speculation=safe
#FHACKS := -heap-arrays 1 -fp-model strict 
#FHACKS := -heap-arrays 1 -fp-speculation=safe -fp-model strict 
#FHACKS := -heap-arrays 1 -fp-speculation=safe -fp-model precise 

FINCLUDE := -I$(ROOTDIR)/$(ARCH)/rcamod  -I$(ROOTDIR)/include/ 

ifeq (USING_NETCDF,$(findstring USING_NETCDF,$(DEFINITIONS)))
	FINCLUDE += -I$(NETCDF_DIR)/include/
	NCDF := -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -Wl,-rpath,$(NETCDF_DIR)/lib
endif

FCFLAGS := $(OPT) $(FDBG) $(FHACKS) $(DEFINITIONS) $(MPI_SPECIAL) $(FINCLUDE) $(TRACE)

# xlf (IBM/AIX), ifort (Intel/Linux) (and others) always produce a lower-case module name
MODNAME = $(shell echo $(*F) | tr "[:upper:]" "[:lower:]")
MODEXT := mod

LD           := $(FC) $(DBG)
#LD          := $(FC) $(DBG) -pg
LD_MPP      :=  -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -lstdc++ $(MPI_SPECIAL) $(NCDF)

AR      := xiar
ARFLAGS := rv
MV      := mv
RM      := rm -f
MKDIR   := mkdir
RMDIR   := rmdir
