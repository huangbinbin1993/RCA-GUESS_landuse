///////////////////////////////////////////////////////////////////////////////////////
/// \file emdi.h
/// \brief Extra code used by the EMDI benchmarks
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_EMDI_H
#define LPJ_GUESS_EMDI_H

#include <gutil.h>
#include <map>
#include <utility>
#include <sstream>

namespace {
std::map<std::pair<double, double>, double> pawcPerGridCell;
}

void rememberPAWC(double dlon, double dlat, xtring pawc) {
	std::istringstream is((char*)pawc);
	double first_number_in_string;
	is >> first_number_in_string;

	pawcPerGridCell[std::make_pair(dlon, dlat)] = first_number_in_string;
}

void overrideAWC(double lon, double lat, Soiltype& soiltype) {
	double pawc = pawcPerGridCell[std::make_pair(lon, lat)];
	if (pawc > 0) {
		// Assume pawc applies to the whole 1.5m, so replace soiltype.awc with a scaled pawc
		double scaled = pawc / 300.0; // as pawc applies to the upper 30cm
		soiltype.awc[0]=SOILDEPTH_UPPER*scaled;
		soiltype.awc[1]=SOILDEPTH_LOWER*scaled;
	}
}

#endif // LPJ_GUESS_EMDI_H
