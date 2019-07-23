///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H

#include "guess.h"

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell);

#endif // LPJ_GUESS_LANDCOVER_H
