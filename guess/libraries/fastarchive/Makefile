######################################################################
# Linux/Unix make file for FastArchive library
# Written by Ben Smith, 22 July 2006

# Name of object file
fastarchiveo= fastarchive.o

# Include file location
includepath= .

# Compiler options
options= -w -O3 -I$(includepath)

# Compiler
CC= g++

# Nothing to change past this point
######################################################################

fastarchivec= fastarchive.cpp
fastarchiveh= fastarchive.h

$(fastarchiveo): $(fastarchivec) $(includepath)/$(fastarchiveh)
	$(CC) $(options) -c -o $(fastarchiveo) $(fastarchivec)
