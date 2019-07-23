///////////////////////////////////////////////////////////////////////////////////////
// FRAMEWORK SOURCE CODE FILE
//
// Framework:             LPJ-GUESS Combined Modular Framework
//                        Includes modified code compatible with "fast" cohort/
//                        individual mode - see canexch.cpp
//                        Version adapted for linkage to RCA at SMHI
//                        --> Framework called by RCA every fast timestep
//                        --> GUESS called by framework at end of every day
//                        Updated 2006-11-28:
//                        For coupling to ECHAM5: returns with guessid=0
//                        (no initialisation) until time is 00:00:00
//                        Updated 2006-12-13:
//                        Grid cells with <= 1% land according to ECOCLIMAP
//                        treated as non-land grid cells (frland=0)
// Header file name:      guessmain.h
// Source code file name: guessmain.cpp
// Written by:            Ben Smith
// Version dated:         2007-12-13

#include "config.h"
#include "guess.h"
#include "guessmain.h"
#include "canexch.h"
#include <stdarg.h>


void initlog(int inst) {

	printf("Attempting to open log file for instance number %d\n",inst);

	const char* file_log_prefix = "guess";

	xtring filename;
	filename.printf("%s%d.log", file_log_prefix, inst);

	set_shell(new CommandLineShell((char*)filename));
}

double fpc_from_lai(double lai) {

	// Returns projective cover from LAI according to Beer's law formulation
	
	if (lai<=0.0) return 0.0;
	return 1.0-lambertbeer(lai);
}

/// Creates RCA results based on what GUESS knows about LAIs and land cover fractions
/** This function takes care of switching from GUESS perspective to RCA perspective,
 *  converting things from stand basis to tile or grid cell basis.
 *
 *  Note that even though RCA has forest and open land tiles, and GUESS has stands
 *  with and without trees allowed, these concepts are not exactly the same.
 *  For instance, a GUESS "forest" stand with no trees (due to climatic reasons)
 *  will be seen by RCA as open land.
 *
 *  Basically, GUESS stands are about what is prescribed (such as "no trees allowed"),
 *  and RCA tiles are about what actually is.
 */
template<typename real>
void create_RCA_results(double laiphen_grass_opl,
                        double laiphen_grass_for,
                        double laiphen_conifer,
                        double laiphen_broadleaf,
                        double laimax_conifer,
                        double laimax_broadleaf,
                        double fgrid_land,
                        double fland_crop,
                        double fland_bare,
                        real& veg_opl, 
                        real& lai_opl, 
                        real& lai_conif,
                        real& lai_broad,
                        real& lai_under,
                        real& frdecid,
                        real& frfor) {

	lai_opl = laiphen_grass_opl;

	if (!negligible(fland_crop+fland_bare)) {
		veg_opl = fland_crop/(fland_crop+fland_bare);

		if (lai_opl < 1.0) {
					
			// Very sparse vegetation - reduce vegetated fraction
			// to FPC
					
			double frac = fpc_from_lai(lai_opl);
			if (!negligible(frac)) {
				lai_opl /= frac;
			}
			veg_opl*=frac;
		}
	}
	else
		veg_opl=0.0;

	// Fractional cover of land areas by forest
	// taking account of maximum growing season FPC
	// Areas with tree FPC<0.4 (tree lai<1) are taken to be stunted tree growth, e.g.
	// alpine areas, and count towards open land vegetation
			
	double fpcmax_forest = fpc_from_lai(laimax_conifer+laimax_broadleaf);

	if (fpcmax_forest < 0.4) {
				
		// Sparse woody vegetation (scrub) - transfer "trees" to open land tile
				
		// Reduce vegetation fraction of scrub area to FPC
				
		double lai_scrub = laiphen_broadleaf+laiphen_conifer+laiphen_grass_for;

		double frac = fpc_from_lai(lai_scrub);
		if (!negligible(frac)) {
			lai_scrub /= frac; 
		}
				
		frfor = 0.0;

		// New open land LAI is weighted average of original (open land)
		// LAI and additional ("forest") LAI
				
		lai_opl = lai_opl*veg_opl*(fland_crop+fland_bare)+
			lai_scrub*frac*(1.0-fland_crop-fland_bare);
				
		// New veg_opl is weighted average of vegetated fraction of 
		// crops+bare ground area and vegetated fraction of "forest" area
				
		veg_opl = frac*(1.0-fland_crop-fland_bare)+veg_opl*
			(fland_crop+fland_bare);
					
		// Convert lai_opl to vegetated area basis
				
		if (!negligible(veg_opl)) {
			lai_opl /= veg_opl;
		}
		// else keep current lai_opl (even though it won't be used)
			
		lai_conif = 0.0;
		lai_broad = 0.0;
		lai_under = 0.0;
		frdecid = 0.0;
	}
	else {
			
		// Calculate forest cover (grid cell basis)
				
		double frac = (1.0-fland_crop-fland_bare);
		if (frac < 0.0) {
			frac = 0.0;
		}

		frfor = fgrid_land*frac;
				
		// Calculate broadleaf forest cover (fraction of forest tile)
		// Divide up the forest tile in proportion to the FPCs of
		// broadleaved and coniferous trees
				
		double fpcmax_broadleaf = fpc_from_lai(laimax_broadleaf);
		double fpcmax_conifer = fpc_from_lai(laimax_conifer);
		double fpcsum = fpcmax_broadleaf+fpcmax_conifer;
				
		// Note that RCA's frdecid is really fraction of
		// broadleaved, not deciduous, trees
		if (!negligible(fpcsum)) {
			frdecid = fpcmax_broadleaf/fpcsum;
		}
		else {
			frdecid=0.0;
		}
					
		// Convert tree LAIs from stand basis to deciduous/evergreen fraction basis
		// BLARP! temporary: include understorey LAI in tree LAI
				
		if (!negligible(frdecid)) {
			lai_broad = laiphen_broadleaf/frdecid+laiphen_grass_for;
		}
		else {
			lai_broad=0.0;
		}
					
		if (!negligible(1.0-frdecid)) {
			lai_conif = laiphen_conifer/(1.0-frdecid)+laiphen_grass_for;
		}
		else {
			lai_conif = 0.0;
		}
				
		//lai_under = laiphen_grass_for;
		lai_under = 0.0; // BLARP! temporary
	}

	// Validate/constrain return values

	if (frfor < 0.0) frfor = 0.0;
	else if (frfor > 0.98) frfor = 0.98;

	if (veg_opl < 0.1) veg_opl = 0.1;
	else if (veg_opl > 1.0) veg_opl = 1.0;

	if (lai_opl < 0.1) lai_opl = 0.1;

	if (lai_conif < 0) lai_conif = 0.0;
	if (lai_broad < 0) lai_broad = 0.0;
	if (lai_under < 0) lai_under = 0.0;
	if (lai_conif+lai_broad+lai_under < 0.1) {
		lai_conif = 0.1;
		lai_broad = 0.1;
		lai_under = 0.1;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// INTERFACE RCA-GUESS
// RCA calls lai_guess for each grid cell / pixel and minimum
// time step. Guess recognises whether a particular grid cell is being called for the
// first time and initialises if necessary

FROUTINE(lai_guess) (
                     // [IN/OUT] zero-based id for tile in GUESS World (see guess.cpp)
                     FINTEGER guessid,

                     // [IN] latitude for centre of tile (+=N, -=S)
                     FREAL alat,

                     // [IN] longitude for centre of tile (deg E of greenwich)
                     FREAL along,

                     // [IN] year (from RCA), arbitrary base value
                     FINTEGER yearc,
                     
                     // [IN] month (from RCA), 1-12
                     FINTEGER monthc,

                     // [IN] day of month (from RCA), 1-31
                     FINTEGER dayc,
                     
                     // [IN] hour of day (from RCA), 0-23
                     FINTEGER hourc,
                     
                     // [IN] minute of hour (from RCA), 0-59
                     FINTEGER minutec,

                     // [IN] second of minute (from RCA), 0-59
                     FINTEGER secondc,

                     // [IN/OUT] fraction of grid cell / pixel that is land
                     FREAL frland,
                     
                     // [IN] ambient CO2 concentration (ppm)
                     //      NB: RCA sends in CO2 equivalents, which we shouldn't use
                     FREAL /*co2*/,
                     
                     // [IN] instantaneous temperature 2m above ground (K)
                     //      (for open land tile)
                     FREAL temp_2m_opl_inst,

                     // [IN] instantaneous temperature forest canopy (K)
                     //      (for forest tile)
                     FREAL temp_can_for_inst,

                     // [IN] soil temperatures for soil layers 1-5 (K)
                     //      (for open land tile)
                     FREAL temp_soil_opl_1,
                     FREAL temp_soil_opl_2,
                     FREAL temp_soil_opl_3,
                     FREAL temp_soil_opl_4,
                     FREAL temp_soil_opl_5,

                     // [IN] soil temperatures for soil layers 1-5 (K)
                     //      (for forest tile)
                     FREAL temp_soil_for_1,
                     FREAL temp_soil_for_2,
                     FREAL temp_soil_for_3,
                     FREAL temp_soil_for_4,
                     FREAL temp_soil_for_5,

                     // [IN] soil water content AWC (0-1) - surface layer 0-50 cm
                     //      (for open land tile)
                     FREAL soilw_surf_opl,
                     // [IN] soil water content AWC (0-1) - deep layer 50-150 cm
                     //      (for open land tile)
                     FREAL soilw_deep_opl,

                     // [IN] soil water content AWC (0-1) - surface layer 0-50 cm
                     //      (for forest tile)
                     FREAL soilw_surf_for,
                     // [IN] soil water content AWC (0-1) - deep layer 50-150 cm
                     //      (for forest tile)
                     FREAL soilw_deep_for,

                     // [IN] instantaneous net downward shortwave radiation (W/m2)
                     //      (for open land tile)
                     FREAL swrad_net_opl_inst,

                     // [IN] instantaneous net downward shortwave radiation (W/m2)
                     //      (for forest tile)
                     FREAL swrad_net_for_inst,

                     // [IN] rain (kg/(m2 s) = mm/s)
                     FREAL rain, 
                     
                     // [IN] snow (kg/(m2 s))
                     FREAL snow,
                     
                     // [OUT] fraction of open land tile that is vegetated
                     FREAL veg_opl,

                     // [OUT] leaf area index for vegetated part of tile
                     FREAL lai_opl,

                     // [OUT] leaf area index for conifers today
                     FREAL lai_conif,

                     // [OUT] leaf area index for broadleaved trees today
                     FREAL lai_broad,

                     // [OUT] leaf area index for grassy understorey today
                     FREAL lai_under,

                     // [OUT] fraction of forest cover that is broadleaved (deciduous or evergreen)
                     FREAL frdecid,

                     // [OUT] fraction of GRID CELL (not land!!) that is forest (not open)
                     FREAL frfor
                     ) {

	// Return with guessid=0 if tile not yet initialised and time is not 00:00:00
	// This was needed for ECHAM5 (and perhaps other GCMs?), we simply ignore
	// the grid cell until midnight

	//if (!(FGET(along)>-59 && FGET(along)<-58.5 && FGET(alat)>-28.5 && FGET(alat)<-28)){
	// testing for sepecific point, MCW
	//	return;
	//}	

	if (!FGET(guessid) && !(FGET(hourc)==0 && FGET(minutec)==0 && FGET(secondc)==0)) {
	//if (!FGET(guessid) && !(FGET(hourc)==0 && FGET(minutec)==0 && FGET(secondc)==0) && !(FGET(along)>-59 && FGET(along)<-58.5 && FGET(alat)>-28.5 && FGET(alat)<-28)) {
		return;
	}

	
	// air temperature (degC)
	double temp_opl = FGET(temp_2m_opl_inst)-K2degC;  // conversion from K to degC
	double temp_for = FGET(temp_can_for_inst)-K2degC; // conversion from K to degC

	// soil temperature (25 cm depth) (degC)
	double temp_soil_opl = FGET(temp_soil_opl_3)-K2degC; // use layer 3 soil temperature for now - is this correct? BLARP! 
	double temp_soil_for = FGET(temp_soil_for_3)-K2degC; // use layer 3 soil temperature for now - is this correct? BLARP! 


	// Convert longitude to GUESS standard
	double lon; // longitude (deg +=E, -=W)

	if (FGET(along)>180.0) {
		lon = FGET(along)-360.0;
	}
	else {
		lon = FGET(along);
	}
	
	// Set fraction of grid cell that is land to fraction provided by
	// RCA by default (in case we are not reading from a land cover file)
	double fgrid_land = FGET(frland);
	double fland_crop,fland_bare;

	double laiphen_grass_opl;
	double laiphen_grass_for;
	double laiphen_conifer;
	double laiphen_broadleaf;
	double laimax_conifer;
	double laimax_broadleaf;

	guess_coupled(FGET(guessid), lon, FGET(alat), FGET(yearc),
		FGET(monthc)-1, FGET(dayc)-1,
		FGET(hourc), FGET(minutec), FGET(secondc),
		temp_opl, FGET(swrad_net_opl_inst),
		temp_soil_opl,FGET(soilw_surf_opl),FGET(soilw_deep_opl),
		temp_for, FGET(swrad_net_for_inst),
		temp_soil_for, FGET(soilw_surf_for), FGET(soilw_deep_for),
		FGET(rain), FGET(snow),
		laiphen_grass_opl, laiphen_grass_for, laiphen_conifer,
		laiphen_broadleaf, laimax_conifer, laimax_broadleaf,
		fgrid_land, fland_crop, fland_bare);

	// Fraction of grid cell that is land (from ECOCLIMAP database, or PNV)
	// Ben 2007-12-13: Grid cells with <= 1% land treated as ocean
		
	if (fgrid_land>0.01) {
		
		FGET(frland)=fgrid_land;

		// Return forest cover, vegetated open land fraction and stand-level LAIs for
		// open land, broadleaved forest and conifer forest
		

		// Transfer LAIs from call to GUESS to arguments used by RCA
		// (applies whether or not vegetation feedbacks are enabled --
		// if no feedbacks, stored LAIs from control period are used)

		// Use prognostic or saved LAIs and ECOCLIMAP cover fractions
			
		// Vegetated fraction of open land tile
		// = crop area (unless LAI<1)

		create_RCA_results(laiphen_grass_opl,
		                   laiphen_grass_for,
		                   laiphen_conifer,
		                   laiphen_broadleaf,
		                   laimax_conifer,
		                   laimax_broadleaf,
		                   fgrid_land,
		                   fland_crop,
		                   fland_bare,
		                   FGET(veg_opl),
		                   FGET(lai_opl),
		                   FGET(lai_conif),
		                   FGET(lai_broad),
		                   FGET(lai_under),
		                   FGET(frdecid),
		                   FGET(frfor));

	}
	else {
	
		// Not land, so return null data
		// Ben 2007-12-13: Now also for grid cells with <= 1% land
			
		FGET(frland)=0.0;
		FGET(lai_conif)=0.0;
		FGET(lai_broad)=0.0;
		FGET(lai_under)=0.0;
		FGET(lai_opl)=0.0;
		FGET(frfor)=0.0;
		FGET(frdecid)=0.0;
		FGET(veg_opl)=0.0;
	}

}
