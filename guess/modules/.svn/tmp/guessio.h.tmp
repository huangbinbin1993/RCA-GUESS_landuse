///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio.h
/// \brief LPJ-GUESS input/output module with input from instruction script
///
/// \author Ben Smith
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GUESSIO_H
#define LPJ_GUESS_GUESSIO_H

#include "guess.h"

bool find_nearest_cru05(float lon,float lat,float& lon_cru,float& lat_cru,int& soilcode);
bool during_rca_spinup(const Date& date);
xtring get_spinup_type();

void initio(const xtring& insfilename, int inst);

bool getgridcell(Gridcell& gridcell);
bool getclimate(Gridcell& gridcell);
void load_rca_spinup(Gridcell& gridcell);
void getlandcover(Gridcell& gridcell);
void outannual(Gridcell& gridcell);
void termio();
void printhelp();

#endif // LPJ_GUESS_GUESSIO_H
