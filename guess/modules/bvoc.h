///////////////////////////////////////////////////////////////////////////////////////
/// \file bvoc.h
/// \brief The BVOC module header file
///
/// Calculation of VOC production and emission by vegetation.
///
/// \author Guy Schurgers (using Almut's previous attempts)
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions
// defined in the module that are to be accessible to the calling framework or
// to other modules.

#ifndef LPJ_GUESS_BVOC_H
#define LPJ_GUESS_BVOC_H

#include "guess.h"

void bvoc(double temp, double hours, double rad, Climate& climate, Patch& patch,
		Individual& indiv, const Pft& pft, const PhotosynthesisResult& phot,
		double adtmm, const Day& day);
void initbvoc();

#endif // LPJ_GUESS_BVOC_H
