///////////////////////////////////////////////////////////////////////////////////////
/// \file driver.h
/// \brief Environmental driver calculation/transformation
///
/// \author Ben Smith
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_DRIVER_H
#define LPJ_GUESS_DRIVER_H

#include "guess.h"

double randfrac(long& seed);
void soilparameters(Soiltype& soiltype,int soilcode);
void interp_monthly_means(double mvals[12], double dvals[366]);
void interp_monthly_totals(double mvals[12], double dvals[366]);
void interp_monthly_means_conserve(const double* mvals, double* dvals);
void interp_monthly_totals_conserve(const double* mvals, double* dvals);
void distribute_ndep(const double* mndry, const double* mnwet,
                     const double* dprec, double* dndep);
void prdaily(double mval_prec[12],double dval_prec[366],double mval_wet[12], long seed);
void dailyaccounting_gridcell(Gridcell& gridcell);
void dailyaccounting_climate(Climate& climate);
void dailyaccounting_stand(Stand& stand);
void dailyaccounting_patch(Patch& patch);
void respiration_temperature_response(double temp,double& gtemp);
void daylengthinsoleet(Climate& climate);
void soiltemp(Climate& climate,Soil& soil);

#endif // LPJ_GUESS_DRIVER_H
