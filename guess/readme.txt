                        LPJ-GUESS Version 2.1
                        =====================

                        PLEASE READ CAREFULLY
 
This directory contains source code files and other files necessary to
run a *demonstation version* of LPJ-GUESS:
     - in population mode (i.e. as LPJ-DGVM) or cohort mode (i.e. as GUESS)
     - on the Unix or Windows operating system
     - in Windows, as either a command-line executable or within the
       LPJ-GUESS Windows Shell

Note that the concept of "plug and play" does not apply to LPJ-GUESS.
It is the responsibility of each individual user to:
     - provide a C++ compiler. For installation under Windows a software
       development "platform" is also strongly recommended. Documentation
       is available for Microsoft Visual C++.
     - optain the appropriate climatic, atmospheric CO2 and soil data
       to drive the model for their particular study
     - provide source code to read their particular input data and
       produce their particular output
     - produce an instruction script (ins file) with run-time parameters
       (data file names, model configuration settings, PFT parameter values
       etc) for the model runs
Some technical advice regarding these matters can be found in the draft
documentation (reference/guessdoc.pdf). Additional help is available as
commenting in model source code and header files.

The model is compiled with the CMake build system. If it is not available
on your system you can download it for free from www.cmake.org. You can
find further information on how to compile the model with CMake in the
draft documentation (reference/guessdoc.pdf).


IMPORTANT: FTP-TRANSFERS FROM UNIX TO WINDOWS:
The following files should be transferred in "ASCII" mode:
     - Files with the extensions .cpp, .h, .grd, .txt, .dat and .ins
     - Makefile
The following files should be transferred in "binary" mode:
     - Files with the extensions .lib, .exe, .pdf, .bin and .gz
     

Structure of this directory:

./reference
     - Current draft of LPJ-GUESS technical manual

./framework
     - Model framework source code and header files (for explanation see
       technical manual)

./modules
     - Source code and header files for model modules (except "Main"
       module). Input/output module (guessio.cpp/guessio.h) is a
       demonstration version compatible with the input data supplied in
       directory ../data. Most users will be able to write an input/output
       module customised to their own study by modifying the demonstration
       version supplied. Further explanations in the technical manual.

./libraries
     - Source code and header files for the custom libraries gutil and plib, 
       required by LPJ-GUESS.

./command_line_version
     - Files required specifically to install LPJ-GUESS as a command-line
       executable on Unix or Windows.

./windows_version
     - Files required specifically to install LPJ-GUESS as a dynamic link
       library (DLL) to run under the LPJ-GUESS Windows Shell. This is the
       recommended configuration for running the model under Windows.

./data
     - Input data files for the demonstration version of LPJ-GUESS.
       Subdirectories contain:
       ./env
            - historical climate data, soil data. See "readme" file.
       ./ins
            - sample instruction script (ins) files for cohort ("guess")
              and population ("lpj") modes, compatible with demonstration
              input/output module. See commenting in files.
       ./gridlist
            - sample grid cell coordinate list files compatible with
              demonstration input/output module and instruction script
              files under directory ../ins

./cru
     - Input/output module version for reading in CRU historical climate data
       for 1901-2006 in custom binary format used by LPJ-GUESS. The data 
       files themselves (GZIP-compressed) are located in directory cru at the 
       same level as this version of LPJ-GUESS was downloaded from.
       
       NOTE: The CRU data set is not in the public domain and should not be
       downloaded without specific permission. Note also that the file is
       very large and may take many hours to transfer.
       Direct enquiries to Joe Lindstrom (e-mail below).
       

Joe Lindstrom
joe.lindstrom@nateko.lu.se
