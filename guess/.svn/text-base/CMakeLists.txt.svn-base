################################################################################
# CMake configuration file for building LPJ-GUESS
#
# To build LPJ-GUESS with this build system cmake needs to be installed.
# If it's not installed it can be downloaded for free from www.cmake.org.
#

cmake_minimum_required(VERSION 2.8)
project(guess)

# Compiler flags for building with Microsoft Visual C++
if (MSVC)
  # Disable warnings about using secure functions like sprintf_s instead of
  # regular sprintf etc.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _CRT_SECURE_NO_WARNINGS")
endif()

# We get compiler flags from the RCA build system,
# if these contain ipo (interprocedural optimization)
# we need to make sure we link libraries with xiar.
# In CMake this is done automatically if the
# INTERPROCEDURAL_OPTIMIZATION property is switched on.
if (CMAKE_CXX_FLAGS MATCHES "-ipo")
  set_directory_properties(PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
endif()

# The following are configuration variables for the build. 
# Rather than editing them here, you should either edit them
# in the graphical cmake tools (cmake-gui on Windows, ccmake
# on Unix) or pass them in as parameters to cmake.

# The IO variable decides which IO module to use. 
set(IO "cru" CACHE STRING "The IO module to use (cru/demo)")

# A variable controlling whether or not to include unit tests
# Unit tests are disabled in old VS 6.0, since CATCH doesn't compile in such 
# an old compiler.
if (NOT MSVC_VERSION EQUAL "1200")
  set(UNIT_TESTS "OFF" CACHE BOOL "Whether to include unit tests")
endif()

if (UNIX)
  # Setup the SYSTEM variable, currently only used to choose which 
  # submit.sh to generate (for submitting to job queue)


  # Figure out what value it should have initially, based on the
  # environment variable ARCH if it's set.
  if (NOT $ENV{ARCH} STREQUAL "")
	set(DEFAULT_SYSTEM $ENV{ARCH})
  else()
	set(DEFAULT_SYSTEM "")
  endif()

  set(SYSTEM ${DEFAULT_SYSTEM} CACHE STRING "System to build for (empty (=simba), gimle, platon, alarik or multicore)")
endif (UNIX)

# Where to search for cmake modules
# (used by cmake when the include() command is used in a cmake file)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# NetCDF - look for libraries and include files, and use them if found
find_package(NetCDF QUIET)

if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${NETCDF_LIBRARIES})
  add_definitions(-DHAVE_NETCDF)
endif()

# MPI - used if found (not needed on Windows)
# RCA: We let the RCA build system configure MPI
#if (NOT CMAKE_HOST_WIN32)
#  find_package(MPI QUIET)
#endif()

# These are deprecated according to documentation in the FindMPI module,
# but for some reason not hidden. Let's not display them for the typical
# LPJ-GUESS user who hasn't got MPI installed.
#mark_as_advanced(MPI_LIBRARY MPI_EXTRA_LIBRARY)

# RCA's build system defines MPI_SRC if MPI is used
# (it will then also determine include directories etc.)
if (CMAKE_CXX_FLAGS MATCHES "MPI_SRC")
#  include_directories(${MPI_INCLUDE_PATH})
#  set(LIBS ${LIBS} ${MPI_LIBRARIES})
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
  add_definitions(-DHAVE_MPI)
endif()

# Where the compiler should search for header files
include_directories(${guess_SOURCE_DIR}/framework ${guess_SOURCE_DIR}/libraries/gutil ${guess_SOURCE_DIR}/libraries/plib ${guess_SOURCE_DIR}/modules ${guess_SOURCE_DIR}/cru/guessio ${guess_SOURCE_DIR}/libraries/fastarchive)

# The following directories contain source code and 
# additional CMakeLists.txt files
add_subdirectory(framework)
add_subdirectory(modules)
add_subdirectory(cru)
add_subdirectory(libraries)

if (UNIT_TESTS)
  add_subdirectory(tests)
endif()


add_library(guess STATIC ${guess_sources})
add_library(fastarchive STATIC libraries/fastarchive/fastarchive.cpp)
add_library(gutil STATIC libraries/gutil/gutil.cpp)

add_executable(physghost tools/physghost.cpp)
target_link_libraries(physghost guess ${LIBS})

add_executable(calc tools/calc.cpp)
target_link_libraries(calc fastarchive gutil)

# Rule for building the unit test binary
if (UNIT_TESTS)
  add_executable(runtests ${guess_sources} ${test_sources})
endif()


add_executable(readdriver tools/readdriver.cpp)
target_link_libraries(readdriver fastarchive gutil)

add_executable(rtseries tools/rtseries.cpp)
target_link_libraries(rtseries gutil)

add_executable(rtslice tools/rtslice)
target_link_libraries(rtslice gutil)


# Rule for running unit tests automatically
if (UNIT_TESTS)
  add_custom_command(TARGET runtests 
    POST_BUILD 
    COMMAND runtests 
    COMMENT "Running tests")
endif()

string(TOUPPER ${IO} IO)
add_definitions(-DUSE_${IO}_IO)

if (UNIX)
   # pgCC 6 doesn't seem to recognize -rdynamic, so remove it
   # (we shouldn't need it anyway)
   # It seems the CMake developers have fixed this in newer versions
   # (sometime after 2.8)
   SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
   SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
endif(UNIX)

if (UNIX)
  # Set default build type to Release on Unix
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
endif(UNIX)