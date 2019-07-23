///////////////////////////////////////////////////////////////////////////////////////
/// \file euroflux.h
/// \brief Extra code used by the Euroflux benchmarks
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_EUROFLUX_H
#define LPJ_GUESS_EUROFLUX_H

#include <gutil.h>

// guess2008 - euroflux 
/// The value used for missing data in the EUROFLUX files.
const double MISSING_DATA = -9999.0;

/// The number of EUROFLUX years, currently 1996-2006, inclusive.
const int NFLUXYEARS=11;

/// Type for storing EUROFLUX grid cell information
struct EurofluxData {

	xtring desc; 

	int plantation_year;
	int num_dominant_species;
	xtring dom_species[5];
	int dom_species_density[5];
	
	// Flux data for a site
	double fluxNEE[NFLUXYEARS][12];
	double fluxAET[NFLUXYEARS][12];
	double fluxGPP[NFLUXYEARS][12];
	double fluxSWC[NFLUXYEARS][12];

	// Modelled flux data for the same site
	double modelNEE[NFLUXYEARS][12];
	double modelAET[NFLUXYEARS][12];
	double modelGPP[NFLUXYEARS][12];

	EurofluxData() {
		// initialise EUROFLUX arrays with missing values; 
		for (int yr = 0; yr < NFLUXYEARS; yr++) {
			for (int mth = 0; mth < 12; mth++) {
				fluxNEE[yr][mth] = MISSING_DATA;
				fluxAET[yr][mth] = MISSING_DATA;
				fluxGPP[yr][mth] = MISSING_DATA;
				fluxSWC[yr][mth] = MISSING_DATA;
			}
		}

		// New initialisation
		plantation_year = -1;
		num_dominant_species = 0;
		desc = ""; 

	}
};

#endif // LPJ_GUESS_EUROFLUX_H
