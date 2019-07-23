///////////////////////////////////////////////////////////////////////////////////////
/// \file rca_driver.h
/// \brief Misc functions for the RCA coupling
///
///
/// \author Joe Siltberg
///
/// $Date: 2012-01-18 11:23:57 +0100 (Wed, 18 Jan 2012) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_RCA_DRIVER_H
#define LPJ_GUESS_RCA_DRIVER_H

#include "lamarquendep.h"

class Stand;

/// To be called every RCA short time step for every vegetated tile
/** Includes per-timestep accounting */
void transfer_drivers(Stand& stand, int& warnct, double lon_rca, double lat_rca, double temp, double netswrad,
                      double wcont_surf_vol, double wcont_deep_vol,
                      double rain, double snow);

/// To be called before we start running a daily time step
/** Should be called once we've received forcing data for all short timesteps.
 *  This function does the averaging of the forcing data to daily values.
 * \see transfer_drivers
 */
void finalize_drivers(Stand& stand);

/// Calculates various LAI terms required by RCA
/** To be called every daily time step for every vegetated tile if vegetation
 *  feedback is switched on.
 */
void compute_tile_lai(Stand& stand);

/// Specify pathname and timeseries for ndep data
void init_RCA_NDep(const char* path, Lamarque::timeseriestype timeseries);

/// Gets an NDepData object for given CRU coordinates
/** Uses an internal cache, so the data isn't read from file every time. */
Lamarque::NDepData& getNDepData(double lon_cru, double lat_cru);

#endif // LPJ_GUESS_RCA_DRIVER_H
