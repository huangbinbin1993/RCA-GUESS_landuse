
.SUFFIXES: .F .F90 .f .f90
TOROOT=..
 
include $(ROOTDIR)/config/config.$(ARCH) #to get the compiler-rules right

VPATH=.:$(TOROOT)/../mainsrc


TAGS =  -Wl,--defsym,__NSCROSSBY_TAGVERSION_1=0
TAGS += -Wl,--defsym,__NSCROSSBY_CODE_RCA=0
TAGS += -Wl,--defsym,__NSCROSSBY_COMPUTER_$(shell uname -n)=0
TAGS += -Wl,--defsym,__NSCROSSBY_ARCH_$(ARCH)=0
TAGS += $(foreach STUFF,$(FCFLAGS), -Wl,--defsym,__NSCROSSBY_FORTRANFLAGS_$(shell echo $(STUFF)|sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0)
TAGS += $(foreach STUFF,$(CXXFLAGS), -Wl,--defsym,__NSCROSSBY_CXXFLAGS_$(shell echo $(STUFF)|sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0)
TAGS += $(foreach STUFF,$(CCFLAGS), -Wl,--defsym,__NSCROSSBY_CFLAGS_$(shell echo $(STUFF)|sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0)
TAGS += -Wl,--defsym,__NSCROSSBY_AR_$(AR)=0
TAGS += $(foreach STUFF,$(LD_MPP), -Wl,--defsym,__NSCROSSBY_LDMPP_$(shell echo $(STUFF)|sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0)
TAGS += -Wl,--defsym,__NSCROSSBY_SVNREVISION_$(shell sed -n '11,11p' $(ROOTDIR)/.svn/entries)=0


TAGS += $(foreach CHFILE,$(shell cat $(ROOTDIR)/localChanges), -Wl,--defsym,__NSCROSSBY_LOCALLY_CHANGED_FILE_$(shell echo $(CHFILE)|sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0)



TAGS += -Wl,--defsym,__NSCROSSBY_FORTRANCOMPILER_$(shell $(FC) --version | sed -n '1,1p' |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0
TAGS += -Wl,--defsym,__NSCROSSBY_CCOMPILER_$(shell $(CC) --version | sed -n '1,1p' |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0
TAGS += -Wl,--defsym,__NSCROSSBY_CXXCOMPILER_$(shell $(CXX) --version | sed -n '1,1p' |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0
TAGS += -Wl,--defsym,__NSCROSSBY_MPIVER1_$(shell mpirun --version | sed -n '1,1p' |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0
TAGS += -Wl,--defsym,__NSCROSSBY_MPIVER2_$(shell mpirun --version | sed -n '2,2p' |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0
TAGS += -Wl,--defsym,__NSCROSSBY_BULDDATE_$(shell date |sed -r 's/[^a-zA-Z0-9_]{1}/_/g')=0


$(ROOTDIR)/$(ARCH)/bin/rca.x: rca.o $(RCALIBS) $(GUESS)
	$(LD) rca.o $(RCALIBS) $(GUESS) $(LIBPSMILE) $(LIBNETCDF) $(TAGS) -o $@ $(LD_MPP)

rca.o: rca.f90
	$(FC) $(FCFLAGS) -c $< 

clean:
	$(RM) *.o
	-$(RM) *.d work.pc
