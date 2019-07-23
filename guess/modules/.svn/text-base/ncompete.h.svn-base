///////////////////////////////////////////////////////////////////////////////////////
/// \file ncompete.h
/// \brief Distribution of N among individuals according to supply, demand and
///        the individuals' uptake strength
///
/// \author David Wårlind
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_NCOMPETE_H
#define LPJ_GUESS_NCOMPETE_H

#include <vector>

/// Represents an individual competing for nitrogen uptake
/** Contains what the ncompete function below needs to know about
 *  an individual in order to do the distribution of N.
 */
struct NCompetingIndividual {
	/// This individual's N demand
	double ndemand;

	/// A meassure of this individual's uptake strength
	double strength;

	/// Output from ncompete - fraction of the demand satisfied by the distribution
	double fnuptake;
};

/// Distributes N among individuals according to supply, demand and uptake strength
/** Grasses should get at least 5% and no
 *  individual should get more than 100% of its nitrogen demand. 
 */
void ncompete(std::vector<NCompetingIndividual>& individuals, 
              double nmass_avail);

#endif // LPJ_GUESS_NCOMPETE_H
