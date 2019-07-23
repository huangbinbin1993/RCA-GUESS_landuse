///////////////////////////////////////////////////////////////////////////////////////
/// \file framework.cpp
/// \brief Implementation of the framework() function
///
/// \author Ben Smith
/// $Date: 2014-06-16 12:30:49 +0200 (Mon, 16 Jun 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "framework.h"
#include "commandlinearguments.h"
#include "guessserializer.h"
#include "parallel.h"

#include "guessio.h"
#include "driver.h"
#include "canexch.h"
#include "soilwater.h"
#include "somdynam.h"
#include "growth.h"
#include "vegdynam.h"
#include "landcover.h"
#include "bvoc.h"

#include <memory>

#include "rca_driver.h"
#include "guessmain.h"
#include "savestate.h"
#include "GuessVeg.h"


/// Verify a coming date, to make sure RCA's and GUESS' calendars agree
/** Checks if a coming date agrees with what the 
 *  Date class thinks the next date should be.
 *
 *  \param date        Current date
 *  \param year        year of the next day according to RCA
 *  \param month       month of the next day according to RCA
 *  \param dayofmonth  day number of the next day according to RCA
 */
bool correct_expected_next_day(const Date& date, 
                               int year, 
                               int month, 
                               int dayofmonth) {
	 Date next_day = date;
	 next_day.next_day();

	 return next_day.year == year &&
		  next_day.month == month &&
		  next_day.dayofmonth == dayofmonth;
}

/// Help function to check if it's time to save state file.
/** If state_interval == "spinup", this function always returns false.
 */
bool time_to_save_state(int year, int month, int dayofmonth, int hour, int minute, int second) {
	// Only save at beginning of a day
	if (hour || minute || second) 
		return false;

	// Don't save if the user has chosen to save only after spinup
	if (state_interval == "spinup")
		return false;

	return (state_interval == "yearly" && month == 0 && dayofmonth == 0) ||
	       (state_interval == "monthly" && dayofmonth == 0);
}

///////////////////////////////////////////////////////////////////////////////////////
// THE GUESS WORLD
// In this version they have to be global to allow sharing with fortran RCA

// The one and only world (vector of gridcells, one forest and one vegetated open land stand
// for each grid cell)
typedef std::vector<Gridcell*> World;
World world;

// Initialisation flag (whether GUESS has been initialised yet true/false)
bool ifinit=false;

// FastArchive containing GUESS-generated vegetation from control period
GuessVegArchive gark;

int inst=0;

// Objects for serializing and deserializing Gridcells
static std::auto_ptr<GuessSerializer> serializer;
static std::auto_ptr<GuessDeserializer> deserializer;

///////////////////////////////////////////////////////////////////////////////////////
// INITGUESS
// Called at start of run to read settings from ins file
// For now ins file name is prescribed in the function body

void initguess() {
	
	// Hard-wired ins file name
	xtring insfilename="guess.ins";
	
	// Get instance id
	inst = GuessParallel::get_rank();

	// Open log file
	initlog(inst);

	// Call input/output module to obtain PFT static parameters and simulation
	// settings and initialise input/output
	initio(insfilename, inst);

	// Nitrogen limitation
	if (ifnlim && !ifcentury) {
		fail("\n\nIf nitrogen limitation is switched on then century soil module also needs to be switched on!");
	}

	// bvoc
	if (ifbvoc) {
	  initbvoc();
	}

	// Open FastArchive containing stored vegetation from control period
	// (needed for non-feedback scenario run)

	if (!ifvegfeedback) {
		bool success=gark.open(file_veg);
		dprintf("file_veg=%s\n",(char*)file_veg);
		if (!success) fail("Could not open %s for input",(char*)file_veg);
		dprintf("Vegetation provided to RCA will be read from %s\n",(char*)file_veg);
	}

	// Done it!
	ifinit=true;
}

/// Simulate one day for a given Stand
/**
 * The climate object in the stand needs to be set up with
 * the day's forcing data before calling this function.
 *
 * \param stand            The stand to simulate
 * \param prescribe_wcont  Whether to use prescribed water levels or run GUESS'
 *                         own hydrology.
 */
void simulate_day(
	 Gridcell& gridcell, 
	 bool prescribe_wcont) {

	 int p;

	 // Update daily climate drivers etc
	 dailyaccounting_gridcell(gridcell);
	 dailyaccounting_climate(gridcell[0].climate);
	 dailyaccounting_climate(gridcell[1].climate);
		
	 if(run_landcover && date.day == 0 && date.year >= nyear_spinup) {
		 // Update dynamic landcover and crop fraction data during historical
		 // period and create/kill stands.
		 landcover_dynamics(gridcell);
	 }

	 gridcell.firstobj();
	 while (gridcell.isobj) {

		 // START OF LOOP THROUGH STANDS

		 Stand& stand=gridcell.getobj();

		 // Calculate daylength, insolation and potential evapotranspiration
		 daylengthinsoleet(stand.climate);

		 dailyaccounting_stand(stand);

		 stand.firstobj();
		 while (stand.isobj) {
			 // START OF LOOP THROUGH PATCHES

			 // Get reference to this patch
			 Patch& patch=stand.getobj();
			 // Update daily soil drivers including soil temperature
			 dailyaccounting_patch(patch);
			 // Leaf phenology for PFTs and individuals
			 leaf_phenology(patch,stand.climate);
			 // Interception
			 interception(patch, stand.climate);

			 if (!prescribe_wcont) {
				 initial_infiltration(patch, stand.climate);
			 }

			 // Photosynthesis, respiration, evapotranspiration
			 canopy_exchange(patch, stand.climate);

			 if (prescribe_wcont) {
				 // During RCA-based spinup, or when ifprescribedwcont is set
				 patch.soil.wcont[0] = stand.wcont_prescribed[0];
				 patch.soil.wcont[1] = stand.wcont_prescribed[1];
			 }
			 else {
				 // Soil water accounting, snow pack accounting
				 soilwater(patch, stand.climate);

				 // NB: It's important not to run the GUESS hydrology
				 // during an RCA-based spinup, since we don't have any
				 // precipitation then. If soilwater is allowed to run,
				 // the evaporative soil layer will eventually become NaN
			 }

			 // Soil organic matter and litter dynamics
			 som_dynamics(patch);

			 if (date.islastday && date.islastmonth) {

				 // LAST DAY OF YEAR
				 // Tissue turnover, allocation to new biomass and reproduction,
				 // updated allometry
				 growth(stand,patch);
			 }
			 stand.nextobj();
		 }// End of loop through patches

		 if (date.islastday && date.islastmonth) {
			 // LAST DAY OF YEAR
			stand.firstobj();
			while (stand.isobj) {

				 // For each patch ...
				 Patch& patch=stand.getobj();
				 // Establishment, mortality and disturbance by fire
				 vegetation_dynamics(stand,patch);
				 stand.nextobj();
			 }
		 }

		 gridcell.nextobj();
	 }	// End of loop through stands

	 if (date.islastday && date.islastmonth) {
		 // LAST DAY OF YEAR
		 // Call input/output module to output results for end of year
		 // or end of simulation for this stand
		 outannual(gridcell);
	 }
}

///////////////////////////////////////////////////////////////////////////////////////
// SPINUP
// Spin up stand to given calender year, month and day using CRU data from nearest
// 0.5 x 0.5 degree land pixel, or archived RCA data if doing an RCA spinup.

const int FIRSTHISTYEAR=1901;
	// first year of historical climate data
	// NB: must match value given in guessio.cpp

void spinup(Gridcell& gridcell,int year,int month,int dayofmonth) {

	// Initialise global variable date
	// (argument nyear not used in this implementation)
	date.init(1);

	// We use simulation years during spinup (starting from 0),
	// so to get the calendar year during spinup we must add
	// the following:
	calendar_year_offset = FIRSTHISTYEAR-nyear_spinup;

	bool spinon=true;

	dprintf("Performing spinup to %04d-%02d-%02d ...\n",year,month+1,dayofmonth+1);

	while (spinon) {

		// START OF LOOP THROUGH SIMULATION DAYS

		// Call input/output to obtain climate, insolation and CO2 for this
		// day of the simulation
		getclimate(gridcell);

		// Do the actual simulation for this day
		simulate_day(gridcell, 
		             // prescribe wcont if this is an RCA-based spinup
					 get_spinup_type()=="cru" ? false : ifprescribedwcont  // no soil water data in CRU dataset, therefore
					  													   // we can only use precp. driven mode during CRU spinup
		             );

		// Advance timer to next simulation day
		date.next_day();

		// Exit loop if RCA day reached
		if (date.year+calendar_year_offset==year && date.month==month && date.dayofmonth==dayofmonth) {

			// In case this is an RCA spinup, load spinup data ready for main routine
			getclimate(gridcell);
			spinon=false;
		}

		//dprintf("cyear=%d date.month=%d date.dayofmonth=%d spinon=%d stand.climate.temp=%g\n",
		//	cyear,date.month,date.dayofmonth,spinon,stand.climate.temp);
		// End of loop through simulation days
	}

	// After spinup, the date object represents real calendar dates
	calendar_year_offset = 0;

	dprintf("Spinup complete\n");
}




///////////////////////////////////////////////////////////////////////////////////////
// GUESS_COUPLED
// To be called every RCA short time step for every pixel / grid cell. 
// Called from guessmain.cpp.

void guess_coupled(int& id,double lon,double lat,
	int year,int month,int dayofmonth,int hour,int minute,int second,
	double temp_opl,double netswrad_opl,
	double temp_soil_opl,double soilw_surf_opl,double soilw_deep_opl,
	double temp_for, double netswrad_for,
	double temp_soil_for, double soilw_surf_for, double soilw_deep_for,
	double rain, double snow,
	double& laiphen_grass_opl, double& laiphen_grass_for,double& laiphen_conifer,double& laiphen_broadleaf,
	double& laimax_conifer,double& laimax_broadleaf,
	double& fgrid_land,double& fland_crop,double& fland_bare) {

	// ARGUMENTS:
	// id                = id from Gridcell object, 0 if initialising
	// lon               = longitude for grid cell
	// lat               = latitude for grid cell
	// year              = calender! year (e.g. 1964)
	// month             = month (0-11)
	// dayofmonth        = what it says (0-30)
	// hour              = hour (0-23)
	// minute            = minute (0-59)
	// second            = second (0-59)
	// temp_opl          = instantaneous temperate 2 m above ground (deg C)
	// netswrad_opl      = net downward shortwave radiation (J/m2/s)
	// temp_soil_opl     = soil temperature (deg C)
	// soilw_surf_opl    = soil water content AWC (0-1) - surface layer 0-50 cm
	// soilw_deep_opl    = soil water content AWC (0-1) - deep layer 50-150 cm
	// temp_for          = instantaneous temperate 2 m above ground (deg C)
	// netswrad_for      = net downward shortwave radiation (J/m2/s)
	// temp_soil_for     = soil temperature (deg C)
	// soilw_surf_for    = soil water content AWC (0-1) - surface layer 0-50 cm
	// soilw_deep_for    = soil water content AWC (0-1) - deep layer 50-150 cm
	// rain              = rain (kg/(m2 s) = mm/s)
	// snow              = snow (kg/(m2 s))
	// laiphen_grass_opl = returned grass LAI in open land tile today on stand basis
	// laiphen_grass_for = returned grass LAI in forest tile today on stand basis
	// laiphen_conifer   = returned conifer LAI today on stand basis
	// laiphen_broadleaf = returned broadleaved LAI today on stand basis
	// laimax_conifer    = returned conifer LAI on stand basis
	// laimax_broadleaf  = returned broadleaved tree summer max LAI on stand basis
	// fgrid_land        = fraction of this grid cell that is land according to ECOCLIMAP database
	//                     NB: now returned to RCA by GUESS (Ben, 2007-04-18)
	// fland_crop        = fraction of LAND on this grid cell that is non-forest
	//                     vegetation due to land use (crops, swamp herbaceous and gardens)
	// fland_bare        = fraction of LAND on this grid cell that is non-vegetated
	//                     (rock, permanent snow or bare soil)

	// NB: LAI returned once per day only

	float lon_cru,lat_cru;
	double wcont_surf,wcont_deep; // soil water content (fAWC)
	bool initialising=false;
	int d,p,soilcode;
	GuessVeg rec;

	// Initialise GUESS if not yet done
	if (!ifinit) {
		initguess();
		deserializer.reset(create_deserializer(state_dir, state_name, year, month, dayofmonth));
	}

	if (!id) {
		
		// Create new Gridcell object and initialise

		// On first call, time is assumed to be 00:00 (Ben 2007-11-22)
		if (hour || minute || second) {
			dprintf("Date is %04d-%02d-%02d %02d:%02d:%02d\n",year,month+1,dayofmonth+1,hour,minute,second);
			fail("guess_coupled: first call of GUESS must be at first daily timestep");
		}
 
		initialising=true;

		world.push_back(new Gridcell());
		Gridcell& gridcell = *world.back();
		gridcell.id = world.size();
		gridcell.date.init(1);
		date.init(1); // date.year needs to be zero when landcover_init/getlandcover is called
			
		// Return id to RCA
		id = gridcell.id;

		// Find nearest half-degree land pixel
		
		dprintf("\nLooking for CRU pixel near (%0.3f,%0.3f) ...\n",lon,lat);
		if (!find_nearest_cru05(lon, lat, lon_cru, lat_cru, soilcode)) {
			dprintf("No suitable CRU pixel found - setting to non-land\n");
			gridcell.island = false;
		}
		else {
			dprintf("Nearest is (%0.1f,%0.1f)\n", lon_cru, lat_cru);
			gridcell.island = true;
		}

		gridcell.lon_cru = lon_cru;
		gridcell.lat_cru = lat_cru;

		gridcell.set_coordinates(lon, lat);

		if (gridcell.island) {
			// Get information about land cover fractions from the chosen LandCoverStrategy
			gridcell.landcoverinfo.reset(landcoverstrategy->get_landcoverinfo(gridcell.lon_cru,
			                                                                  gridcell.lat_cru,
			                                                                  fgrid_land));
			
			if (gridcell.landcoverinfo.get() == 0) {
				// No land here according to the LandCoverStrategy
				gridcell.island = false;
			}
		}

		if (gridcell.island) {

			// Initialise certain climate and soil drivers
			soilparameters(gridcell.soiltype,soilcode);

			// SPINUP
			// If state file with previously saved state exists, read data from there

			if(run_landcover) {
				//Read static landcover and cft fraction data from ins-file and/or from data files for the spinup peroid and create stands.
				landcover_init(gridcell);
			}

			// Ben 2007-12-17
			// Read in spinup data (whether needed or not)
			getgridcell(gridcell);

			// Ben 2007-12-19
			// Now initialising to RCA latitude (even if this is a CRU spinup)
			gridcell[0].climate.initdrivers(gridcell.get_lat());
			gridcell[1].climate.initdrivers(gridcell.get_lat());

			if (deserializer.get()) {
				deserializer->deserialize_gridcell(gridcell);
				dprintf("Loaded state\n");

				// Calculate phenology (not available in state file)
				// Ben 2007-11-22
	
				date.set(year,month,dayofmonth,hour,minute,second);

				for (int i = 0; i < 2; i++) {
					Stand& stand = gridcell[i];

					stand.firstobj();
					while (stand.isobj) {

						// For each patch ...
						Patch& patch = stand.getobj();
						leaf_phenology(patch,stand.climate);

						stand.nextobj();
					}
				}

				// Get RCA spinup data for this day if still spinning up
				if (during_rca_spinup(date)) {

					// Retrieve RCA spinup data for this year and assign to stand climate
					load_rca_spinup(gridcell);
					getclimate(gridcell);
				}
			}
			else {
				// Perform spin up to given year and day
				spinup(gridcell,year,month,dayofmonth);
				
				// Save the state so we don't need a spinup next time
				if (!serializer.get()) {
					serializer.reset(create_serializer(state_dir, 
					                                   state_name, 
					                                   year, 
					                                   month, 
					                                   dayofmonth, 
					                                   inst,
					                                   GuessParallel::get_num_processes()));
				}
				serializer->serialize_gridcell(gridcell);
			}

			gridcell.date.set(year,month,dayofmonth,hour,minute,second);
			gridcell.date.timestep_rca=0;

			dprintf("Simulation time is %04d-%02d-%02d %02d:%02d:%02d\n",
				year,month+1,dayofmonth+1,hour,minute,second);

			if (during_rca_spinup(Date(year, month, dayofmonth, hour, minute, second)))
				dprintf("RCA spinup data will be used up until control period\n");

			// Compute daily LAI terms for return to RCA
			if (ifvegfeedback) {
				compute_tile_lai(gridcell[0]);
				compute_tile_lai(gridcell[1]);
			}
			else {

				// Use stored data from control period in binary archive
			
				dprintf("Driving RCA with data from %s (no vegetation feedback)\n",
					(char*)file_veg);
	
				rec.process_id = inst;
				rec.gridcell_id = gridcell.id;

				if (!gark.getindex(rec)) 
					fail("No data for (%g,%g) available in %s",
					     gridcell.get_lon(), gridcell.get_lat(), (char*)file_veg);

				for (d=0;d<366;d++) {
					gridcell[0].laiphen_grass[d]=rec.laiphen_grass_opl[d];
					gridcell[1].laiphen_grass[d]=rec.laiphen_grass_for[d];
					gridcell[1].laiphen_conifer[d]=rec.laiphen_conifer[d];
					gridcell[1].laiphen_broadleaf[d]=rec.laiphen_broadleaf[d];
				}

				gridcell[1].laimax_conifer=rec.laimax_conifer[0];
				gridcell[1].laimax_broadleaf=rec.laimax_broadleaf[0];		
			}
		}

		date = gridcell.date;

		// Transfer environmental drivers from RCA to climate and soil objects for each stand
		// (not if continuing RCA spinup)

		if (!during_rca_spinup(Date(year, month, dayofmonth, hour, minute, second)) && gridcell.island) {
			transfer_drivers(gridcell[0], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_opl, netswrad_opl,
			                 soilw_surf_opl, soilw_deep_opl, rain, snow);
			transfer_drivers(gridcell[1], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_for, netswrad_for,
			                 soilw_surf_for, soilw_deep_for, rain, snow);
		}
		// End of initialisation
	}

	// Retrieve stand object for this tile
	Gridcell& gridcell = *world[id-1];

	// Return with null data for non-land tiles

	if (!gridcell.island) {

		laiphen_grass_opl = laiphen_grass_for = laiphen_conifer = laiphen_broadleaf = 0.0;
		fgrid_land = fland_crop = fland_bare = 0.0;
		return;
	}
		
	// Update RCA time step
	// Doing it here and this way avoids problems with repeated time steps during restart
	// Ben 2007-11-26

	if (gridcell.date<Date(year,month,dayofmonth,hour,minute,second)) {

		gridcell.date.next_timestep();

		date = gridcell.date;
		
		// Transfer environmental drivers from RCA to climate and soil objects for this stand
		// (unless continuing RCA spinup)

		if (!during_rca_spinup(Date(year, month, dayofmonth, hour, minute, second))) {
			transfer_drivers(gridcell[0], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_opl, netswrad_opl,
			                 soilw_surf_opl, soilw_deep_opl, rain, snow);
			transfer_drivers(gridcell[1], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_for, netswrad_for,
			                 soilw_surf_for, soilw_deep_for, rain, snow);
		}
	}

	// Transfer stand's own internal date to global date
	date = gridcell.date;

	// Do we have the same idea as RCA about things like leap years etc.?
	if (date<Date(year,month,dayofmonth,0,0,0) &&
		 !correct_expected_next_day(date, year, month, dayofmonth)) {
		 fail("RCA sent in unexpected date, check calendar settings in ins file!");
	}
	
	if (date<Date(year,month,dayofmonth,0,0,0) || initialising) {

		// NEW DAY COMING (but GUESS doesn't know yet)

		// Calculate averages of drivers for short timesteps, 
		// to get daily forcing before we simulate this day
		finalize_drivers(gridcell[0]);
		finalize_drivers(gridcell[1]);

		// Do the actual simulation for this day
		simulate_day(gridcell, ifprescribedwcont);
		
		// Advance to next day (resets timestep counter)
		gridcell.date.set(year,month,dayofmonth,hour,minute,second);
		date = gridcell.date;

		// Get RCA spinup data for this day if still spinning up
		if (during_rca_spinup(date)) {

			load_rca_spinup(gridcell);
			getclimate(gridcell);
		}
				
		// Transfer environmental drivers from RCA to climate and soil objects for this stand
		// (the first values for today were lost when dailyaccounting_stand_coupled() was called)
		// For now, calculate LAIs for transfer to RCA here also (change by Ben 060321)
		// (unless still doing an RCA spinup - Ben 2007-12-19)

		if (!during_rca_spinup(date)) {
			transfer_drivers(gridcell[0], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_opl, netswrad_opl,
			                 soilw_surf_opl, soilw_deep_opl, rain, snow);
			transfer_drivers(gridcell[1], gridcell.warnct, gridcell.get_lon(), gridcell.get_lat(), temp_for, netswrad_for,
			                 soilw_surf_for, soilw_deep_for, rain, snow);
		}

		// Compute daily LAI terms for return to RCA
		if (ifvegfeedback) {
			compute_tile_lai(gridcell[0]);
			compute_tile_lai(gridcell[1]);
		}
		else {

			// Use stored data from control period in binary archive

			if (date.day==0) {

				rec.process_id = inst;
				rec.gridcell_id = gridcell.id;

				if (!gark.getindex(rec))
					fail("No data for (%g,%g) available in %s",
					     gridcell.get_lon(),gridcell.get_lat(),(char*)file_veg);

					if (!(date.year%4)) {  // should be replaced by is_leap() ? mcw

						// Leap year -- recycle 28/2 value for 29/2

						for (d=0;d<366;d++) {
							if (d<59) { // 1 Jan-28 Feb
								gridcell[0].laiphen_grass[d]=rec.laiphen_grass_opl[d];
								gridcell[1].laiphen_grass[d]=rec.laiphen_grass_for[d];
								gridcell[1].laiphen_conifer[d]=rec.laiphen_conifer[d];
								gridcell[1].laiphen_broadleaf[d]=rec.laiphen_broadleaf[d];
							}
							else { // 29 Feb-31 Dec
								gridcell[0].laiphen_grass[d]=rec.laiphen_grass_opl[d-1];
								gridcell[1].laiphen_grass[d]=rec.laiphen_grass_for[d-1];
								gridcell[1].laiphen_conifer[d]=rec.laiphen_conifer[d-1];
								gridcell[1].laiphen_broadleaf[d]=rec.laiphen_broadleaf[d-1];
							}
						}					
					}
					else { // Non leap year
						for (d=0;d<365;d++) {
							gridcell[0].laiphen_grass[d]=rec.laiphen_grass_opl[d];
							gridcell[1].laiphen_grass[d]=rec.laiphen_grass_for[d];
							gridcell[1].laiphen_conifer[d]=rec.laiphen_conifer[d];
							gridcell[1].laiphen_broadleaf[d]=rec.laiphen_broadleaf[d];
						}
					}
					gridcell[1].laimax_conifer=rec.laimax_conifer[0];
					gridcell[1].laimax_broadleaf=rec.laimax_broadleaf[0];
			}
			if ((date.month==1 && date.dayofmonth==28 || date.day==365) && date.year%4 && (calendarmode != FLAT_30))
				fail("29 February in non-leap year %d (date.month=%d date.dayofmonth=%d date.day=%d)\n",
					date.year,date.month,date.dayofmonth,date.day);
		}
 
	}

	// Save state file?
	if (date.timestep_rca && gridcell.island && !initialising) {

		// Ben 2007-11-22 - not on first call from RCA
		
		if (time_to_save_state(year, month, dayofmonth, hour, minute, second)) {
			
			if (!serializer.get()) {
				serializer.reset(create_serializer(state_dir, state_name, year, month, dayofmonth, inst, GuessParallel::get_num_processes()));
			}
			serializer->serialize_gridcell(gridcell);
		}
		else {
			// When it's no longer time to save state we'll get rid of the serializer
			// object so that the files get flushed. The next time it's time to
			// save state, we'll create a new serializer object with a new date.
			serializer.reset();
		}
	}
	
	// Return values
	
	laiphen_grass_opl = gridcell[0].laiphen_grass[date.day];
	laiphen_grass_for = gridcell[1].laiphen_grass[date.day];
	laiphen_conifer = gridcell[1].laiphen_conifer[date.day];
	laiphen_broadleaf = gridcell[1].laiphen_broadleaf[date.day];
	laimax_conifer = gridcell[1].laimax_conifer;
	laimax_broadleaf = gridcell[1].laimax_broadleaf;
	
	fgrid_land=gridcell.landcoverinfo->get_fgrid_land();
	fland_crop=gridcell.landcoverinfo->get_fland_crop(year);
	fland_bare=gridcell.landcoverinfo->get_fland_bare(year);
}
