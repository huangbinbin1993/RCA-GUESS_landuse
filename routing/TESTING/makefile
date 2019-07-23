TL = ../routing$(NBITS).a

# Assume variables:
# NBITS   = 32 or 64 as required
# WRAPUP  = wrapup sequences, including ranlib where required

fc = $(FC) -c $(FFLAGS) $*.f90 || echo $< 

# List dependencies on source and included files:

# Fortran objects in this library:

OBJECTS = modvar.o

OJECTS +=  modvar.mod

OBJECTS += river_routing.o routing.o

# C objects in this library:

# Targets:

$(TL) : $(TL)($(OBJECTS))
	$(WRAPUP)

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(TL)

