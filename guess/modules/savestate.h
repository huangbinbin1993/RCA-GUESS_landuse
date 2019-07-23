///////////////////////////////////////////////////////////////////////////////////////
/// \file savestate.h
/// \brief Help functions for serialization
///
/// Contains a couple of functions for creating the serializer and deserializer 
/// objects. Takes care of for instance creating/finding a directory for the state
/// files, corresponding to a base directory and a date.
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_SAVESTATE_H
#define LPJ_GUESS_SAVESTATE_H

#include "guess.h"
#include "guessserializer.h"

/// Creates a serializer object for a given base directory and date
/** It is up to the receiver to make sure that the object gets deleted. */
GuessSerializer* create_serializer(xtring state_dir,
                                   xtring state_name,
                                   int calendar_year,
                                   int month,
                                   int dayofmonth,
                                   int instance,
                                   int num_processes);

/// Creates a deserializer object for a given base directory and date
/** Returns a null pointer if the directory for the state files doesn't
 *  exist. 
 *  It is up to the receiver to make sure that the object gets deleted. */
GuessDeserializer* create_deserializer(xtring state_dir,
                                       xtring state_name,
                                       int calendar_year,
                                       int month,
                                       int dayofmonth);

#endif // LPJ_GUESS_SAVESTATE_H
