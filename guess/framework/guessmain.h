///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK SOURCE CODE FILE
//
// Framework:             LPJ-GUESS Combined Modular Framework
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
//                        Version adapted for linkage to RCA at SMHI
// Header file name:      guessmain.h
// Source code file name: guessmain.cpp
// Written by:            Ben Smith
// Version dated:         2005-01-09

#include "flink.h"

void initlog(int inst);

FROUTINE(lai_guess) (
                     FINTEGER guessid,
                     FREAL alat,
                     FREAL along,
                     FINTEGER yearc,
                     FINTEGER monthc,
                     FINTEGER dayc,
                     FINTEGER hourc,
                     FINTEGER minutec,
                     FINTEGER secondc,
                     FREAL frland,
                     FREAL /*co2*/,
                     FREAL temp_2m_opl_inst,
                     FREAL temp_can_for_inst,
                     FREAL temp_soil_opl_1,
                     FREAL temp_soil_opl_2,
                     FREAL temp_soil_opl_3,
                     FREAL temp_soil_opl_4,
                     FREAL temp_soil_opl_5,
                     FREAL temp_soil_for_1,
                     FREAL temp_soil_for_2,
                     FREAL temp_soil_for_3,
                     FREAL temp_soil_for_4,
                     FREAL temp_soil_for_5,
                     FREAL soilw_surf_opl,
                     FREAL soilw_deep_opl,
                     FREAL soilw_surf_for,
                     FREAL soilw_deep_for,
                     FREAL swrad_net_opl_inst,
                     FREAL swrad_net_for_inst,
                     FREAL rain, 
                     FREAL snow,
                     FREAL veg_opl,
                     FREAL lai_opl,
                     FREAL lai_conif,
                     FREAL lai_broad,
                     FREAL lai_under,
                     FREAL frdecid,
                     FREAL frfor
                     );
