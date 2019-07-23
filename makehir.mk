.DEFAULT:
.SUFFIXES:
.PHONY: depend

# SRCDIR should be CURDIR without the ARCH part
SRCDIR := $(subst $(ARCH)/,,$(CURDIR))

# ROOTDIR is where the first make is started
ROOTDIR := $(PWD)

VPATH := .:$(SRCDIR)
VPATH += $(ROOTDIR)/include

include $(ROOTDIR)/config/config.$(ARCH)

#  GC precision ===========================================
ifeq (PREC32,$(findstring PREC32,$(CCFLAGS)))
   CCFLAGS += -DFLP_32B
endif
ifeq (PREC64,$(findstring PREC64,$(CCFLAGS)))
   CCFLAGS += -DFLP_64B
endif


SRCCXX   := $(notdir $(wildcard $(SRCDIR)/*.cpp)) 
SRCCX :=     $(notdir $(wildcard $(SRCDIR)/*.C))
SRCSC   := $(notdir $(wildcard $(SRCDIR)/*.c))
SRCSF   := $(notdir $(wildcard $(SRCDIR)/*.F))
SRCSFMINOR   := $(notdir $(wildcard $(SRCDIR)/*.f))
SRCSF90 := $(notdir $(wildcard $(SRCDIR)/*.F90))
SRCSF90MINOR := $(notdir $(wildcard $(SRCDIR)/*.f90))
#INKLUDE := $(notdir $(wildcard $(SRCDIR)/*.inc))

OBJSC        := $(SRCSC:.c=.o) 
OBJSCXX      := $(SRCCXX:.cpp=.o)  
OBJSCX      := $(SRCCX:.C=.o)
OBJSF        := $(SRCSF:.F=.o) 
OBJSFMINOR   := $(SRCSFMINOR:.f=.o) 
OBJSF90      := $(SRCSF90:.F90=.o) 
OBJSF90MINOR := $(SRCSF90MINOR:.f90=.o) 



OBJS := $(OBJSF90) $(OBJSF) $(OBJSF90MINOR) $(OBJSFMINOR) $(OBJSC) $(OBJSCXX) $(OBJSCX) 
SRCS := $(SRCSF90) $(SRCSF) $(SRCSF90MINOR) $(SRCSFMINOR) $(SRCSC) $(SRCCXX) $(SRCCX) $(INKLUDE)

$(MAKECMDGOALS): $(OBJS) $(INKLUDE) #The MAKECMDGOALS variable contains a list of all the targets specified on the command line...
	$(AR) $(ARFLAGS) $@ $(OBJS)


DEPSCOND := $(ROOTDIR)/config/definitions.global $(ROOTDIR)/config/config.$(ARCH) 

%.o: %.c $(DEPSCOND)  #PREPROCESSSED
	$(CC) $(CCFLAGS) -c $<

%.o: %.cpp   $(DEPSCOND) #PREPROCESSSED
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.C   $(DEPSCOND) #PREPROCESSSED
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.f90  $(DEPSCOND) #not preprocessed
	$(FC)  $(FCFLAGS) -c $<
	(test -f $(MODNAME).$(MODEXT) && $(MV) $(MODNAME).$(MODEXT) $(ROOTDIR)/$(ARCH)/rcamod)  || echo

%.o: %.f $(DEPSCOND) #not preprocessed
	$(FC) $(FCFLAGS) -c $<
	(test -f $(MODNAME).$(MODEXT) && $(MV) $(MODNAME).$(MODEXT) $(ROOTDIR)/$(ARCH)/rcamod) || echo

%.o : %.F $(DEPSCOND) #PREPROCESSSED
	$(FC)  $(FCFLAGS) -c $<
	(test -f $(MODNAME).$(MODEXT) && $(MV) $(MODNAME).$(MODEXT) $(ROOTDIR)/$(ARCH)/rcamod) || echo	

%.o : %.F90 $(DEPSCOND) #PREPROCESSSED
	$(FC)  $(FCFLAGS) -c $<
	(test -f $(MODNAME).$(MODEXT) && $(MV) $(MODNAME).$(MODEXT) $(ROOTDIR)/$(ARCH)/rcamod) || echo


depend dependencies.make: $(SRCS)  
	rm -f $(SRCDIR)/dependencies.make
	touch $(SRCDIR)/dependencies.make
	echo $(SRCS)
	cd $(SRCDIR) ;  $(if $(filter %.f %.F %.f90 %.F90 , $(SRCS)),$(ROOTDIR)/tools/bin/makedepf90 -W $(DEFINITIONS) -I.:$(ROOTDIR)/include:$(ROOTDIR)/modules:$(ROOTDIR)/$(ARCH)/rcamod:$(ROOTDIR)/interfaces:$(ROOTDIR)/$(ARCH)/interfaces $(filter %.f %.F %.f90 %.F90 %.mod, $(SRCS)) >> $(SRCDIR)/dependencies.make;) $(if $(filter %.c %.cpp %.C, $(SRCS)), $(CXX) $(MPI_SPECIAL) -MM $(filter %.c %.cpp %.C, $(SRCS))  >> $(SRCDIR)/dependencies.make;)
	$(MV) $(SRCDIR)/dependencies.make .

-include dependencies.make

-include $(SRCDIR)/Makefile.qrk

