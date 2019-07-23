////////////////////////////////////////////////////////////////////////////////
/// \file partitionedmapserializer.cpp
/// \brief Implementation file for PartitionedMapSerializer/Deserializer
///
/// Since these classes are templates, most of their implementation is in the
/// header.
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
////////////////////////////////////////////////////////////////////////////////

#include "partitionedmapserializer.h"
#include <sstream>

std::string create_path(const char* directory, int rank) {
	std::ostringstream os;
	os << directory << "/" << rank << ".state";
	return os.str();
}
