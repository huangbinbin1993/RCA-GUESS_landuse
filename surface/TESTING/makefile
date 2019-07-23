#OBJ = reg2rot.o rot2reg.o Data.o  ModelData.o  Gtopo30Data.o main.o Gtopo30DataSet.o read_namelist.o Gtopo30DataSetReader.o

OBJ = decomp.o ecoclimapdatareader.o LakeData.o myRoutines.o EcoclimapDataSet.o Domain.o factorize.o ModelData.o Data.o reg2rot.o rot2reg.o EcoclimapDataSetReader.o EcoclimapData.o texture.o Gtopo30DataSetReader.o Gtopo30Data.o Gtopo30DataSet.o gtopo30reader.o read_namelist.o factorize.o  localParallel.o

DEFINITIONS=  -DUSE_LITTLE_ENDIAN  -DMPI_SRC #-DDEBUG

SRC = ..

INCLUDES = -I. -I$(SRC) 
DBG = -g -traceback #-DDEBUG
OPT = -O0 #-O3 -xH -ip -ipo 

CXX = icc 
FC = ifort 
LD = icc 
LIBS=   -Nmixrpath -Nverbose -lstdc++ -lifcore #-Nverbose

ifneq (,$(findstring DMPI_SRC, $(DEFINITIONS)))
	LIBS += -Nmpi -L/opt/scali/lib64 -lfmpi 
	OPT += -Nmpi
endif

COPT=  $(OPT) -Nmixrpath -Nverbose 
FOPT = $(OPT) -Nmixrpath -Nverbose 

#CXX=gcc 
#FC=gfortran 
#LD=gfortran

CFLAGS = $(DEFINITIONS) $(COPT) $(DBG)
FFLAGS = $(DEFINITIONS) $(FOPT) $(DBG) 
#CXX = icpc

default: main mainE 

decomp.o:../../src/decomp.F90 makefile 
	$(FC) $(FFLAGS) $(INCLUDES)  -c $<

main:$(OBJ) main.o makefile 
	$(LD)  $(PROFIL) -o $@ main.o $(OBJ)  $(LIBS)

mainE:$(OBJ) mainE.o makefile
	$(LD)  $(PROFIL) -o $@ mainE.o $(OBJ)  $(LIBS)


localParallel.o:localParallel.C makefile
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<	

main.o:main.C makefile
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

mainE.o:mainE.C makefile
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<


%.o:$(SRC)/%.C makefile
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

%.o:$(SRC)/%.f90 makefile decomp.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

factorize.o:factorize.f90 makefile decomp.o decomp.mod
	$(FC) $(FFLAGS) $(INCLUDES) -c $<


clean:
	rm *.o main mainE *.*~
