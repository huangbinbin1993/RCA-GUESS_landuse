
.DELETE_ON_ERROR:

ROOTDIR := $(PWD)

-include $(ROOTDIR)/config/definitions.global 
export DEFINITIONS

ifeq (,$(findstring $(ARCH), $(SUPPORTED_ARCH)))
	ARCH := UNSUPPORTED
endif

include $(ROOTDIR)/config/config.$(ARCH) 

RCADIRS :=  src surface
RCALL :=   src surface
RCALIBS := $(patsubst %,$(ROOTDIR)/$(ARCH)/lib/%.a,$(RCADIRS))
RCALLIBS := $(patsubst %,$(ROOTDIR)/$(ARCH)/lib/%.a,$(RCALL)) 

GUESS_SRC := $(ROOTDIR)/guess
GUESS_BUILD_PATH := $(ROOTDIR)/$(ARCH)/lib/guess_build
GUESS := $(GUESS_BUILD_PATH)/libguess.a

export RCALIBS
export GUESS

default: rca.x  

$(ARCH)/localChanges:
	svn st |grep ^M |sed 's/^M       /MODIFIED_FILE_/g' >& $(ARCH)/localChanges
	svn st |grep ^D |sed 's/^D       /DELETED_FILE_/g' >> $(ARCH)/localChanges
	svn st |grep ^C |sed 's/^C       /CONFLICTING_FILE_/g' >> $(ARCH)/localChanges
	svn st |grep ^A |sed 's/^A       /ADDED_FILE_/g' >> $(ARCH)/localChanges
	svn st |grep ^R |sed 's/^R       /REPLACED_FILE_/g' >> $(ARCH)/localChanges
	svn st |grep ^I |sed 's/^I       /IGNORED_FILE_/g' >> $(ARCH)/localChanges



rca.x:  $(RCALLIBS) $(ARCH)/bin $(ARCH)/mainsrc $(ROOTDIR)/config/definitions.global $(ROOTDIR)/config/config.$(ARCH) $(GUESS) $(ARCH)/localChanges
	$(MAKE) -C $(ARCH)/mainsrc ROOTDIR=$(ROOTDIR) ARCH=$(ARCH) $(ROOTDIR)/$(ARCH)/bin/rca.x $(TAGS)

$(ARCH)/mainsrc: ./$(ARCH)
	test -d $@ || $(MKDIR) $@
	-ln -sf ../../mainsrc/Makefile ./$(ARCH)/mainsrc

clean: 
	-$(RM) -r ./$(ARCH)	
	$(MAKE) -C ./tools/makedepf90-2.8.8 clean


#make all libraries

$(RCALLIBS): ./config/config.$(ARCH)  ./$(ARCH)/lib ./$(ARCH)/rcamod  ./tools/bin/makedepf90
	test -d $(ARCH)/$(patsubst $(ROOTDIR)/$(ARCH)/lib/%.a,%,$@) || $(MKDIR) $(ARCH)/$(patsubst $(ROOTDIR)/$(ARCH)/lib/%.a,%,$@)
	$(MAKE) -C $(ARCH)/$(patsubst $(ROOTDIR)/$(ARCH)/lib/%.a,%,$@) -f $(ROOTDIR)/makehir.mk ARCH=$(ARCH) TOROOT=.. $@

$(RCALL): % : $(ROOTDIR)/$(ARCH)/lib/%.a

$(GUESS) : $(GUESS_BUILD_PATH)/Makefile
	$(MAKE) --no-print-directory -C $(GUESS_BUILD_PATH)

$(GUESS_BUILD_PATH)/Makefile: | $(GUESS_BUILD_PATH)
	cd $(GUESS_BUILD_PATH) && CXX=$(CXX) cmake -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS)" $(GUESS_SRC)

./tools/bin/makedepf90:
	$(MAKE) -C ./tools/makedepf90-2.8.8 
	$(MAKE) -C ./tools/makedepf90-2.8.8 install


# MISC tasks
./$(ARCH):
	test -d $@ || $(MKDIR) $@

./$(ARCH)/bin: ./$(ARCH)
	test -d $@ || $(MKDIR) $@

./$(ARCH)/lib: ./$(ARCH)
	test -d $@ || $(MKDIR) $@

./$(ARCH)/rcamod: ./$(ARCH)
	test -d $@ || $(MKDIR) $@

./$(ARCH)/config: ./$(ARCH)
	test -d $@ || $(MKDIR) $@

$(GUESS_BUILD_PATH): | ./$(ARCH)
	test -d $@ || $(MKDIR) $@

.PHONY: ./$(ARCH)/config/config.h $(GUESS)

./$(ARCH)/config/config.h: ./$(ARCH)/config $(ROOTDIR)/config/config.$(ARCH) Makefile
	echo "ROOTDIR=" $(ROOTDIR) > $(ARCH)/config/config.h
	echo "OBJDIR=" $(ARCH) >> $(ARCH)/config/config.h
	echo "RCALIBS=" $(RCALIBS) >> $(ARCH)/config/config.h
	cat $(ROOTDIR)/config/config.$(ARCH) >> $(ARCH)/config/config.h



