///////////////////////////////////////////////////////////////////////////////////////
/// \file savestate.cpp
/// \brief Help functions for serialization
///
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "savestate.h"
#ifdef _MSC_VER
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <sstream>

namespace {

// Creates a directory (the parent directory must exist)
void make_directory(const std::string& path) {
	// Check if the directory already exists, if so just exit
#ifdef _MSC_VER
	struct _stat buf;
	if (!_stat(path.c_str(), &buf) &&
		buf.st_mode & _S_IFDIR) {
#else
	struct stat buf;
	if (!stat(path.c_str(), &buf) &&
	    S_ISDIR(buf.st_mode)) {
#endif
		return;
	}

#ifdef _MSC_VER
	if (_mkdir(path.c_str())) {
#else
	if (mkdir(path.c_str(), 0777)) {
#endif
		if (errno != EEXIST) {
			fail("Failed to create directory %s\n%s\n", 
			     path.c_str(), strerror(errno));
		}
	}
}

// Converts an integer to a string
std::string to_string(int d) {
	std::ostringstream os;
	os << d;
	return os.str();
}

}

GuessSerializer* create_serializer(xtring state_dir,
                                   xtring state_name,
                                   int calendar_year,
                                   int month,
                                   int dayofmonth,
                                   int instance,
                                   int num_processes) {

	// first create the directory and its parents
	
	std::string full_path((char*)state_dir);

	full_path += std::string("/") + (char*)state_name + "_state_" + to_string(calendar_year);
	make_directory(full_path);
	
	full_path += "/" + to_string(month+1);
	make_directory(full_path);
	
	full_path += "/" + to_string(dayofmonth+1);
	make_directory(full_path);
	
	// now that we have a directory, we can create the serializer
	return new GuessSerializer(full_path.c_str(), instance, num_processes);
}

GuessDeserializer* create_deserializer(xtring state_dir,
                                       xtring state_name,
                                       int calendar_year,
                                       int month,
                                       int dayofmonth) {

	// verify that the directory exists, otherwise we'll return a null pointer

	std::string full_path((char*)state_dir);
	full_path += std::string("/") + (char*)state_name + "_state_" + to_string(calendar_year);

	full_path += "/" + to_string(month+1);

	full_path += "/" + to_string(dayofmonth+1);

#ifdef _MSC_VER
	struct _stat buf;
	if (!_stat(full_path.c_str(), &buf) &&
		buf.st_mode & _S_IFDIR) {
#else
	struct stat buf;
	if (!stat(full_path.c_str(), &buf) &&
	    S_ISDIR(buf.st_mode)) {
#endif
		return new GuessDeserializer(full_path.c_str());
	}
	else {
		return 0;
	}
}
