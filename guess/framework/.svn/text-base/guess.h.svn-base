///////////////////////////////////////////////////////////////////////////////////////
/// \file guess.h
/// \brief Framework header file, LPJ-GUESS Combined Modular Framework
///
/// This header file contains:
///  (1) definitions of all main classes used by the framework and modules. Modules may
///      require classes to contain certain member variables and functions (see module
///      source files for details).
///  (2) other type, constant and function definitions to be accessible throughout the
///      model code.
///  (3) a forward declaration of the framework function if this is not the main
///      function.
///
/// \author Ben Smith
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_GUESS_H
#define LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// #INCLUDES FOR LIBRARY HEADER FILES
// C/C++ libraries required for member functions of classes defined in this file.
// These libraries will also be available globally (so omit these #includes from source
// files). In addition to various standard C/C++ runtime libraries, the framework
// requires the following libraries (individual modules may use additional libraries)
//
// GUTIL
//   Includes class xtring, providing functionality for pointer-free dynamic handling
//   of character strings; wherever possible in LPJ-GUESS, strings are represented as
//   objects of type xtring rather than simple arrays of type char. GUTIL also provides
//   templates for dynamic collection classes (list arrays of various types), argument
//   processing for printf-style functions, timing functions and other utilities.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "gutil.h"
#include "config.h"
#include "landcoverstrategies.h"
#include <vector>
#include "shell.h"
#include "guessmath.h"
#include "archive.h"

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL ENUMERATED TYPE DEFINITIONS

/// Life form class for PFTs (trees, grasses)
typedef enum {NOLIFEFORM, TREE, GRASS} lifeformtype;

/// Phenology class for PFTs
typedef enum {NOPHENOLOGY, EVERGREEN, RAINGREEN, SUMMERGREEN, ANY} phenologytype;

/// Biochemical pathway for photosynthesis (C3 or C4)
typedef enum {NOPATHWAY, C3, C4} pathwaytype;

/// Leaf physiognomy types for PFTs
typedef enum {NOLEAFTYPE, NEEDLELEAF, BROADLEAF} leafphysiognomytype;

/// Units for insolation driving data
/** Insolation can be expressed as:
 *
 *  - Percentage sunshine
 *  - Net instantaneous downward shortwave radiation flux (W/m2)
 *  - Total (i.e. with no correction for surface albedo) instantaneous downward 
 *    shortwave radiation flux (W/m2)
 *
 *  Radiation flux can be interpreted as W/m2 during daylight hours, or averaged
 *  over the whole time step which it represents (24 hours in daily mode). For
 *  this reason there are two enumerators for these insolation types (e.g. SWRAD
 *  and SWRAD_TS).
 */
typedef enum {
	/// No insolation type chosen
	NOINSOL,
	/// Percentage sunshine
	SUNSHINE,
	/// Net shortwave radiation flux during daylight hours (W/m2)
	NETSWRAD,
	/// Total shortwave radiation flux during daylight hours (W/m2)
	SWRAD,
	/// Net shortwave radiation flux during whole time step (W/m2)
	NETSWRAD_TS,
	/// Total shortwave radiation flux during whole time step (W/m2)
	SWRAD_TS,
	/// Photoactive radiation
	PAR
} insoltype;

/// Vegetation 'mode', i.e. what each Individual object represents
/** Can be one of: 
 *  1. The average characteristics of all individuals comprising a PFT
 *     population over the modelled area (standard LPJ mode)
 *  2. A cohort of individuals of a PFT that are roughly the same age
 *  3. An individual plant
 */
typedef enum {NOVEGMODE, INDIVIDUAL, COHORT, POPULATION} vegmodetype;

/// CENTURY pool names, NSOMPOOL number of SOM pools
typedef enum {SURFSTRUCT, SOILSTRUCT, SOILMICRO, SURFHUMUS, SURFMICRO, SURFMETA, SURFFWD, SURFCWD,
	SOILMETA, SLOWSOM, PASSIVESOM, LEACHED, NSOMPOOL} pooltype;	

/// Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
typedef enum {URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, NLANDCOVERTYPES} landcovertype;

/// Water uptake parameterisations
/** \see water_uptake in canexch.cpp
  */
typedef enum {WR_WCONT, WR_ROOTDIST, WR_SMART, WR_SPECIESSPECIFIC} wateruptaketype;

/// Decides how many days there are per month
/**
 *  STANDARD  is like a real non-leap year, every year.
 *  LEAPYEARS is the same, but february has 29 days in leap years 
 *            (years divisible by 4).
 *  FLAT_30   means all months have 30 days.
 */
typedef enum {STANDARD, LEAPYEARS, FLAT_30} calendartype;


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL CONSTANTS

/// number  of soil layers modelled
const int NSOILLAYER = 2;

// SOIL DEPTH VALUES

/// soil upper layer depth (mm)
const double SOILDEPTH_UPPER = 500.0;
/// soil lower layer depth (mm)
const double SOILDEPTH_LOWER = 1000.0;

/// Year at which to calculate equilibrium soil carbon
const int SOLVESOM_END=400;

/// Year at which to begin documenting means for calculation of equilibrium soil carbon
const int SOLVESOM_BEGIN = 350;

/// Number of years to average growth efficiency over in function mortality
const int NYEARGREFF = 5; 

/// Coldest day in N hemisphere (January 15)
/** Used to decide when to start counting GDD's and leaf-on days 
 *  for summergreen phenology.
 */
const int COLDEST_DAY_NHEMISPHERE = 14;

/// Coldest day in S hemisphere (July 15)
/** Used to decide when to start counting GDD's and leaf-on days 
 *  for summergreen phenology.
 */
const int COLDEST_DAY_SHEMISPHERE = 195;

/// number of years to average aaet over in function soilnadd
const int NYEARAAET = 5;

/// Maximum number of age classes in age structure plots produced by function outannual
/** not used in RCA-GUESS, set to 1 to reduce memory usage */
const int OUTPUT_MAXAGECLASS = 1;

/// Priestley-Taylor coefficient (conversion factor from equilibrium evapotranspiration to PET)
const double PRIESTLEY_TAYLOR = 1.32;

// Solving Century SOM pools 

/// fraction of nyear_spinup minus freenyears at which to begin documenting for calculation of Century equilibrium
const double SOLVESOMCENT_SPINBEGIN  = 0.1;
/// fraction of nyear_spinup minus freenyears at which to end documentation and start calculation of Century equilibrium
const double SOLVESOMCENT_SPINEND    = 0.3;

/// Kelvin to deg c conversion
const double K2degC = 273.15;

/// Conversion factor for CO2 from ppmv to mole fraction
const double CO2_CONV = 1.0e-6;


///////////////////////////////////////////////////////////////////////////////////////
// FORWARD DECLARATIONS OF CLASSES DEFINED IN THIS FILE
// Forward declarations of classes used as types (e.g. for reference variables in some
// classes) before they are actually defined

class Date;
class Stand;
class Patch;
class Vegetation;
class Gridcell;
class Patchpft;

///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH EXTERNAL LINKAGE
// These variables are defined in the framework source code file, and are accessible
// throughout the code

/// Object describing timing stage of simulation
extern Date date;

/// Vegetation mode (population, cohort or individual)
extern vegmodetype vegmode;

/// Number of patches in each stand (should always be 1 in population mode)
extern int npatch;

/// Patch area (m2) (individual and cohort mode only)
extern double patcharea;


/// Whether background establishment enabled (individual, cohort mode)
extern bool ifbgestab;

/// Whether spatial mass effect enabled for establishment (individual, cohort mode)
extern bool ifsme;

/// Whether establishment stochastic (individual, cohort mode)
extern bool ifstochestab;

/// Whether mortality stochastic (individual, cohort mode)
extern bool ifstochmort;

/// Whether fire enabled
extern bool iffire;

/// Whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
extern bool ifdisturb;

/// Generic patch-destroying disturbance interval (individual, cohort mode)
extern double distinterval;
extern xtring state_dir;
	// path to directory in which state will be saved or loaded from
	// (with or without trailing '/')
extern xtring state_name;
	// identifier string for this particular state file
extern xtring state_interval;
	// when to generate state files
	// "spinup" to generate state file only at end of spinup
	// "yearly" to generate state file at the beginning of each year
	// "monthly" to generate state file at the beginning of each month

/// Whether SLA calculated from leaf longevity (alt: prescribed)
extern bool ifcalcsla;

/// Whether leaf C:N ratio minimum calculated from leaf longevity (alt: prescribed)
extern bool ifcalccton;

/// Establishment interval in cohort mode (years)
extern int estinterval;

/// Number of possible PFTs
extern int npft;

/// Whether C debt (storage between years) permitted
extern bool ifcdebt;

/// Water uptake parameterisation
extern wateruptaketype wateruptake;

/// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
extern bool ifcentury;
/// whether plant growth limited by available N	
extern bool ifnlim;
/// number of years to allow spinup without nitrogen limitation	
extern int freenyears;
/// fraction of nitrogen relocated by plants from roots and leaves	
extern double nrelocfrac;
/// first term in nitrogen fixation eqn (Cleveland et al 1999)	
extern double nfix_a;
/// second term in nitrogen fixation eqn (Cleveland et al 1999)	
extern double nfix_b;

/// Whether other landcovers than natural vegetation are simulated.
extern bool run_landcover;

/// Whether a specific landcover type is simulated (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND).
extern bool run[NLANDCOVERTYPES];

/// Whether landcover fractions are read from ins-file.
extern bool lcfrac_fixed;

/// Set to false by initio( ) if fraction input files have yearly data.
extern bool all_fracs_const;

extern bool ifslowharvestpool; 	// If a slow harvested product pool is included in patchpft.
extern int nyear_spinup; // number of spinup years (ML)	Moved to guess.cpp to be accessed globally.
/// silt/sand fractions per soiltype
extern bool textured_soil;

// For coupling to RCA
extern bool ifvegfeedback;
	// whether vegetation changes feed back to RCA
extern bool ifprescribedwcont;
	// whether to use water content sent in from RCA, or use GUESS' own
	// hydrology driven by precipitation

class GlobalCO2File;
extern GlobalCO2File global_co2;
	// Stores co2 values for each year, read from text file

extern xtring file_veg;

extern int first_control_year;
	// first control year of the simulation

extern int calendar_year_offset;
	// offset to add to simulation year to get calendar year

extern calendartype calendarmode;
	// global calendartype used by the model throughout the simulation
	// (varies depending on boundary data set for RCA)

///////////////////////////////////////////////////////////////////////////////////////
// Settings controlling the saving and loading from state files

/// Location of state files
extern xtring state_path;

/// Whether to restart from state files
extern bool restart;

/// Whether to save state files
extern bool save_state;

/// Save/restart year
extern int state_year;


///////////////////////////////////////////////////////////////////////////////////////
// guess2008 - new input variables, from the .ins file
extern bool ifsmoothgreffmort;
	// whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard
	// step-function (0)
extern bool ifdroughtlimitedestab;
	// whether establishment is limited by growing season drought
extern bool ifrainonwetdaysonly;
	// rain on wet days only (1, true), or a little every day (0, false);
// bvoc
extern bool ifbvoc;
        // whether BVOC calculations are included

///////////////////////////////////////////////////////////////////////////////////////
// Global varibles for coupling
extern xtring couple_src;
	// "RCA": online couple with RCA
	// "RCA_SP": Driven by RCA spinup data (used in phyghost)
	// "CRU": Driven by CRU spinup data (used in phyghost)

// For RCA version: function to be called every short time step
// Defined in guess.cpp, called from guessmain.cpp

void guess_coupled(int& id, double lon, double lat,
	int year, int month, int dayofmonth, int hour, int minute, int second,
	double temp_opl, double netswrad_opl,
	double temp_soil_opl, double soilw_surf_opl, double soilw_deep_opl,
	double temp_for, double netswrad_for,
	double temp_soil_for, double soilw_surf_for, double soilw_deep_for,
	double rain, double snow,
	double& laiphen_grass_opl, double& laiphen_grass_for, double& laiphen_conifer, double& laiphen_broadleaf,
	double& laimax_conifer, double& laimax_broadleaf,
	double& fgrid_land, double& fland_crop, double& fland_bare);


/// General purpose object for handling simulation timing.
/** In general, frameworks should use a single Date object for all simulation
 *  timing.
 *
 *  Member variables of the class (see below) provide various kinds of calender
 *  and timing information, assuming init has been called to initialise the
 *  object, and next() has been called at the end of each simulation day.
 */
class Date {

	// MEMBER VARIABLES

public:

	/// Maximum number of days in an LPJ-GUESS simulation year

	// To support leap year, we need to extend the vector size to 366
	// Without this, the variable following a size of 365 vector in a class would be overwritten 
	// in daily update operation, e.g. mwcontupper following dwcontupper[Date::MAX_YEAR_LENGTH] in class Soil
	// in dailyaccounting_patch - mcw100716
	static const int MAX_YEAR_LENGTH = 366;

	// calender year corresponding to first year in CRU climate data set
	// should be set in the ins file - mcw
	static const int FIRSTHISTYEAR = 1901;

	/// number of days in each month (0=January - 11=December)
	int ndaymonth[12];

	/// julian day of year (0-364; 0=Jan 1)
	int day;

	/// day of current month (0=first day)
	int dayofmonth;

	/// month number (0=January - 11=December)
	int month;

	/// year since start of simulation (0=first simulation year)
	int year;

	/// number of subdaily periods in a day (to be set in IO module)
	int subdaily;

	/// julian day for middle day of each month
	int timestep;

	/// RCA timestep (not reset every day)
	int timestep_rca;

	/// 0-based hour from RCA
	int hour;
	/// 0-based minute from RCA
	int minute;
	/// 0-based second from RCA
	int second;
		
	int middaymonth[12];

	/// true if last year of simulation, false otherwise
	bool islastyear;

	/// true if last month of year, false otherwise
	bool islastmonth;

	/// true if last day of month, false otherwise
	bool islastday;

	/// true if middle day of month, false otherwise
	bool ismidday;

	/// number of days in year
	int ndayyear;

	// spinup_type
	xtring spinup_type;

	// calendar year
	int calendar_year;

private:

	int nyear;

	// MEMBER FUNCTIONS

	 /// Updates fields that depend on which calendartype we're using and current year
	 /**
	  * Updates ndaymonth, middaymonth and ndayyear.
	  * This function should be called whenever the year changes.
	  */
	 void update_calendar() {
		  // set up ndaymonth

		  static const int standard_months[] = {31,28,31,30,31,30,31,31,30,31,30,31};

		  // Before RCA day reached, standard calendar have to be used (mcw220216)
		  // This applies to both CRU and RCA spinup forcing and ensure ndaymonth 
		  // is not influenced by calendarmode setting for the coupling period. 
		  // As during this period we don't have leap year concept
		  // BLARP! year is "raw" year and has not been updated yet during spinup
		  //if (Date::calendar_year == 1979) 
			//  fprintf(stdout, "Hello, I am here");
		  if (calendar_year < first_control_year) {
			  for (int i = 0; i < 12; i++) {
				  ndaymonth[i] = standard_months[i];
			  }
		  }
		  else {
			  switch (calendarmode) {
			  case STANDARD:
				  for (int i = 0; i < 12; i++) {
					  ndaymonth[i] = standard_months[i];
				  }
				  break;
			  case LEAPYEARS:
				  for (int i = 0; i < 12; i++) {
					  ndaymonth[i] = standard_months[i];
				  }
				  if (year % 4 == 0) {
					  ndaymonth[1] = 29;
				  }
				  break;
			  case FLAT_30:
				  for (int i = 0; i < 12; i++) {
					  ndaymonth[i] = 30;
				  }
				  break;
			  default:
				  fail("Date::update_calendar: unknown calendartype");
			  }
		  }
		  
		  // set up middaymonth
		  int dayct = 0;
		  for (int m = 0; m < 12; m++) {
				middaymonth[m] = dayct + ndaymonth[m]/2;
				dayct += ndaymonth[m];
		  }

		  ndayyear = dayct;
	 }

public:

	/// Constructor function called automatically when Date object is created
	/** Do not call explicitly. Initialises some member variables. */
	Date() {

		  // Default constructor with initialisation to all zero
		  init(1);
		  set(0,0,0,0,0,0);
	 }

	 Date(int y,int m,int d,int h,int n,int s) {

		  // Constructor with initialisation to given date and time
		  init(1);
		  set(y,m,d,h,n,s);

	 }

	/// Initialises date to day 0 of year 0 and sets intended number of simulation years
	/** Intended number of simulation years is only used to set islastyear flag,
	 *  actual simulation may be longer or shorter.
	 *
	 *  \param nyearsim  Intended number of simulation years
	 */
	 void init(int nyearsim) {
		  nyear = nyearsim;
		  day = month = year = dayofmonth = 0;
		  // Would be a problem if called outside the programe framework (e.g. in Physghost) - mcw
		  calendar_year = FIRSTHISTYEAR - nyear_spinup;
		  update_calendar();
		  hour = minute = second = 0;
		  islastmonth = islastday = ismidday = false;
		  if (nyear == 1) islastyear = true;
		  else islastyear = false;
		  timestep_rca = 0;
		  timestep = 0;
		  subdaily = 1;
	 }

	 void set(int y,int m,int d,int h,int n,int s) {
		  int i;
		  if (y != year || m != month || d != dayofmonth) {
				timestep = 0;
		  }
		  year = y;
		  // assume date.set is called during RCA-GUESS coupling period
		  // and the input y should be calendar year -  mcw
		  calendar_year = y;
		  update_calendar();
		  month = m;
		  dayofmonth = d;
		  day = d;
		  hour = h;
		  minute = n;
		  second = s;
		  for (i = 0; i < m; i++) {
				day += ndaymonth[i];
		  }
		  if (dayofmonth == ndaymonth[m] - 1) islastday = true;
		  else islastday = false;
		  if (m == 11) islastmonth = true;
		  else islastmonth = false;
		  if (y == nyear - 1) islastyear = true;
		  else islastyear = false;
		  if (dayofmonth == ndaymonth[m] / 2) ismidday = true;
		  else ismidday = false;
	 }

	 /// Call at end of every simulation day to update member variables.
	 void next_day() {
	
		  int length[12];
		  int m;
		  for (m = 0; m < 12; m++) {
				length[m] = ndaymonth[m];
		  }
		  if (islastday) {
				if (islastmonth) {
					 dayofmonth = 0;
					 day = 0;
					 month = 0;
					 year++;
					 calendar_year++;
					 update_calendar();
					 if (year == nyear-1) islastyear = true;
					 islastmonth = false;
				}
				else {
					 day++;
					 dayofmonth = 0;
					 month++;
					 if (month == 11) islastmonth = true;
				}
				islastday = false;
		  }
		  else {
				day++;
				dayofmonth++;
				if (dayofmonth == length[month] / 2) ismidday = true;
				else {
					 ismidday = false;
					 if (dayofmonth == length[month] - 1) islastday = true;
				}
		  }
		  timestep = 0;
		  hour = minute = second = 0;
	 }

	 bool istoday(int y,int m,int d) {
	
		  if (y == year && m == month && d == dayofmonth) return true;
		  return false;
	 }

	 void next_timestep() {
		  timestep++;
		  timestep_rca++;
	 }

	// \returns index (0-11) of previous month (11 if currently month 0).
	int prevmonth() {
		if (month > 0) return month - 1;
		return 11;
	}

	/// \returns index of next month (0 if currently month 11)
	int nextmonth() {
		if (month < 11) return month+1;
		return 0;
	}

	 /// \returns date and time as a scalar value
	 double scalar() {

		  return (double)second+(double)minute*60.0+(double)hour*3600.0+
				(double)dayofmonth*86400.0+(double)month*2678400.0+
				(double)year*32140800.0;
	 }

	 bool operator>(Date d) {

		  // Return true if currently stored date and time is later than
		  // the one specified as the arguments to this function

		  if (scalar()>d.scalar()) return true;

		  return false;	 
	 }

	 bool operator<(Date d) {

		  // Return true if currently stored date and time is earlier than
		  // the one specified as the arguments to this function

		  if (scalar()<d.scalar()) return true;

		  return false;
	 }


	/// Check if the year is leap
	/** \param year        Calendar year
	*   The algorith is as follows: only year that are divisible by 4 could
	*   potentially be leap (e.g., 1904), however, not if they're divisble by
	*   100 (e.g., 1900 is not leap), unless they're divisble by 400 (e.g., 2000
	*   is still leap).
	*/
	static bool is_leap(int year) {
		return (!(year % 4) && (year % 100 | !(year % 400)));
	}
	
	/// Whether the current mode is diurnal
	bool diurnal() const { return subdaily > 1; }

	/// \returns the number of days in the current simulation year
	/** For this function to work properly in simulations with varying number
	 *  of days per year, the set_first_calendar_year must have been called first.
	 *
     */
	int year_length() const {
		return is_leap(calendar_year) ? 366 : 365;
// We need to uniform the leap year algo
// now year%4 is used in most of the program - mcw
	}
};

/// Object describing sub-daily periods
class Day {
public:
	/// Whether sub-daily period first/last within day (both true in daily mode)
	bool isstart, isend;
	
	/// Ordinal number of the sub-daily period [0, date.subdaily)
	int period;

	/// Constructs beginning of the day period (the only one in daily mode)
	Day() {
		isstart = true;
		isend = !date.diurnal();
		period = 0;
	}

	/// Advances to the next sub-daily period
	void next() {
		period++;
		isstart = false;
		isend = period == date.subdaily - 1;
	}
};

/// This struct contains the result of a photosynthesis calculation.
/** \see photosynthesis */  
struct PhotosynthesisResult : public Serializable {
	/// Constructs an empty result
	PhotosynthesisResult() {
		clear();
	}

	/// Clears all members
	/** This is returned by the photosynthesis function when no photosynthesis
	 *  takes place.
	 */
	void clear() {
		agd_g       = 0;
		adtmm       = 0;
		rd_g        = 0;
		vm          = 0;
		je          = 0;
		nactive_opt = 0.0;
		vmaxnlim    = 1.0;
	}

	/// RuBisCO capacity (gC/m2/day)
	double vm;

	/// gross daily photosynthesis (gC/m2/day)
	double agd_g;

	/// leaf-level net daytime photosynthesis
	/** expressed in CO2 diffusion units (mm/m2/day) */
    double adtmm;

	/// leaf respiration (gC/m2/day)
	double rd_g;

	/// PAR-limited photosynthesis rate (gC/m2/h)
    double je;

	/// optimal leaf nitrogen associated with photosynthesis (kgN/m2)
	double nactive_opt;

	/// nitrogen limitation on vm
	double vmaxnlim;

	/// net C-assimilation (gross photosynthesis minus leaf respiration) (kgC/m2/day)
    double net_assimilation() const {
		return (agd_g - rd_g) * 1e-3;
    }

	void serialize(ArchiveStream& arch);
};


/// The Climate for a stand
/** Stores all static and variable data relating to climate parameters, as well as
 *  latitude, atmospheric CO2 concentration and daylength for a grid cell. Includes
 *  a reference to the parent Gridcell object (defined below). Initialised by a
 *  call to initdrivers.
 */
class Climate : public Serializable {

	// MEMBER VARIABLES

public:
	/// mean air temperature today (deg C)
	double temp;

	/// total daily net downward shortwave solar radiation today (J/m2/day)
	double rad;

	/// total daily photosynthetically-active radiation today (J/m2/day)
	double par;

	/// precipitation today (mm)
	double prec;

	/// day length today (h)
	double daylength;

	/// atmospheric ambient CO2 concentration today (ppmv)
	double co2;

	/// latitude (degrees; +=north, -=south)
	double lat;
		// latitude (degrees; +=north, -=south)

	/// Insolation today, see also instype
	double insol;

	/// Type of insolation
	/** This decides how to interpret the variable insol,
	 *  see also documentation for the insoltype enum.
	 */
	insoltype instype;

	/// equilibrium evapotranspiration today (mm/day)
	double eet;

	/// mean temperature for the last 31 days (deg C)
	double mtemp;

	/// mean of lowest mean monthly temperature for the last 20 years (deg C)
	double mtemp_min20;

	/// mean of highest mean monthly temperature for the last 20 years (deg C)
	double mtemp_max20;

	/// highest mean monthly temperature for the last 12 months (deg C)
	double mtemp_max;

	/// accumulated growing degree day sum on 5 degree base 
	/** reset when temperatures fall below 5 deg C */
	double gdd5;

	/// total gdd5 (accumulated) for this year (reset 1 January)
	double agdd5;

	/// number of days with temperatures <5 deg C 
	/** reset when temperatures fall below 5 deg C; 
	 *  maximum value is number of days in the year */
	int chilldays;

	/// true if chill day count may be reset by temperature fall below 5 deg C
	bool ifsensechill;

	/** Respiration response to today's air temperature incorporating damping of Q10
	 *  due to temperature acclimation (Lloyd & Taylor 1994)
	 */
	double gtemp;

	/// daily temperatures for the last 31 days (deg C)
	Historic<double, 31> dtemp_31;

	/// minimum monthly temperatures for the last 20 years (deg C)
	double mtemp_min_20[20];

	/// maximum monthly temperatures for the last 20 years (deg C)
	double mtemp_max_20[20];

	/// minimum monthly temperature for the last 12 months (deg C)
	double mtemp_min;
		
	/// mean of monthly temperatures for the last 12 months (deg C)
	double atemp_mean;

	/// annual nitrogen deposition (kgN/m2/year)
	double andep;
	/// daily nitrogen deposition (kgN/m2)
	double dndep;

	/// annual nitrogen fertilization (kgN/m2/year)
	double anfert;
	/// daily nitrogen fertilization (kgN/m2/year)
	double dnfert;

	// Monthly sums (converted to means) used by canopy exchange module

	double temp_mean;
		// accumulated mean temperature for this month (deg C)
	double par_mean;
		// accumulated mean daily net PAR sum (J/m2/day) for this month
	double co2_mean;
		// accumulated mean CO2 for this month (ppmv)
	double daylength_mean;
		// accumulated mean daylength for this month (h)

	// Saved parameters used by function daylengthinsoleet

	double sinelat;
	double cosinelat;
	double qo[Date::MAX_YEAR_LENGTH];
	double u[Date::MAX_YEAR_LENGTH];
	double v[Date::MAX_YEAR_LENGTH];
	double hh[Date::MAX_YEAR_LENGTH];
	double sinehh[Date::MAX_YEAR_LENGTH];
	double daylength_save[Date::MAX_YEAR_LENGTH];
	/// indicates whether saved values exist for this day
	bool doneday[Date::MAX_YEAR_LENGTH];
		
	/// diurnal temperature range, used in daily/monthly BVOC (deg C)
	double dtr;

	std::vector<double> temps, insols, pars, rads, gtemps;
		// containers for sub-daily values of temperature, short-wave downward
		// radiation, par, rad and gtemp (equivalent to temp, insol, par, rad and gtemp)
		// NB: units of these variable are the same as their daily counterparts,
		// i.e. representing daily averages (e.g. pars [J/m2/day])

	// Parameters involved in linkage to RCA
	
	// Cumulative sums of RCA fast-timestep climate parameters
	// averaged at end of day in finalize_drivers()
	
	double temp_sum; // air temperature
	double insol_sum; // net downward shortwave radiation
	double prec_sum; // rain + snow

	double temp_max; // max air temperature from RCA today
	double temp_min; // min air temperature from RCA today
	
	// Climate parameters stored during year for output at end of year (RCA)
	
	double mtemp_out[12];
	double mpar_out[12];
	double mco2_out[12];
	double mprec_out[12];

public:
	Climate() {};
		// constructor function: initialises gridcell member

	void initdrivers(double latitude) {

		// Initialises certain member variables
		// Should be called before Climate object is applied to a new grid cell

		int day, year;

		for (year=0; year<20; year++) {
			mtemp_min_20[year] = 0.0;
			mtemp_max_20[year] = 0.0;
		}
		mtemp = 0.0;
		gdd5 = 0.0;
		chilldays = 0;
		ifsensechill = true; //  guess2008 - CHILLDAYS
		atemp_mean = 0.0;

		lat = latitude;
		for (day = 0; day < 366; day++) doneday[day] = false;
		sinelat = sin(lat * DEGTORAD);
		cosinelat = cos(lat * DEGTORAD);
		
		temp_sum=insol_sum=prec_sum=0.0;
		temp_min = 1000;
		temp_max = -1000;
	}

	void serialize(ArchiveStream& arch);
};


/// Stores accumulated monthly and annual fluxes.
/** This class handles the storage and accounting of fluxes for a single patch.
 *  Different fluxes can be stored in different ways, depending on what kind of
 *  flux it is and what kind of output we want. The details of whether fluxes
 *  are stored per PFT or just as a patch total, or per day, month or only a
 *  yearly sum, is hidden from the 'scientific' code, which merely reports the
 *  fluxes generated.
 */
class Fluxes : public Serializable {

public:

	/// Fluxes stored as totals for the whole patch
	enum PerPatchFluxType { 
		/// Carbon flux to atmosphere from burnt vegetation and litter (kgC/m2)
		FIREC,
		/// Carbon flux to atmosphere from soil respiration (kgC/m2)
		SOILC,
		/// Flux from atmosphere to vegetation associated with establishment (kgC/m2)
		ESTC,
		/// Flux to atmosphere from consumed harvested products (kgC/m2)
		HARVESTC,
		/// Flux to atmosphere from consumed harvested products (kgN/m2)
		HARVESTN,
		/// NH3 flux to atmosphere from fire
		NH3_FIRE,
		/// NO flux to atmosphere from fire	
		NO_FIRE,
		/// NO2 flux to atmosphere from fire
		NO2_FIRE,
		/// N2O flux to atmosphere from fire	
		N2O_FIRE,
		/// N2 flux to atmosphere from fire	
		N2_FIRE,
		/// N flux from soil
		N_SOIL,
		/// Reproduction costs
		REPRC,
		/// Number of types, must be last
		NPERPATCHFLUXTYPES
	};

	/// Fluxes stored per pft
	enum PerPFTFluxType { 
		/// NPP (kgC/m2)
		NPP,
		/// GPP (kgC/m2)
		GPP,
		/// Autotrophic respiration (kgC/m2)
		RA,
		/// Isoprene (mgC/m2)
		ISO,
		/// Monoterpene (mgC/m2)
		MON,
		/// Number of types, must be last
		NPERPFTFLUXTYPES
	};

	// emission ratios from fire (NH3, NO, NO2, N2O) Delmas et al. 1995
	// values in .cpp file

	static const double NH3_FIRERATIO;
	static const double NO_FIRERATIO;
	static const double NO2_FIRERATIO;
	static const double N2O_FIRERATIO;
	static const double N2_FIRERATIO;

	/// Reference to patch to which this Fluxes object belongs
	Patch& patch;

	// MEMBER FUNCTIONS

public:
	/// constructor: initialises members
	Fluxes(Patch& p);

	/// Sets all fluxes to zero (call at the beginning of each year)
	void reset();

	void serialize(ArchiveStream& arch);

	/// Report flux for a certain flux type
	void report_flux(PerPFTFluxType flux_type, int pft_id, double value);

	/// Report flux for a certain flux type
	void report_flux(PerPatchFluxType flux_type, double value);

	/// \returns flux for a given month and flux type (for all PFTs)
	double get_monthly_flux(PerPFTFluxType flux_type, int month) const;

	/// \returns flux for a given month and flux type
	double get_monthly_flux(PerPatchFluxType flux_type, int month) const;

	/// \returns annual flux for a given PFT and flux type
	double get_annual_flux(PerPFTFluxType flux_type, int pft_id) const;

	/// \returns annual flux for a given flux type (for all PFTs)
	double get_annual_flux(PerPFTFluxType flux_type) const;

	/// \returns annual flux for a given flux type
	double get_annual_flux(PerPatchFluxType flux_type) const;

private:

	/// Stores one flux value per PFT and flux type
	std::vector<std::vector<double> > annual_fluxes_per_pft;

	/// Stores one flux value per month and flux type
	/** For the fluxes only stored as totals for the whole patch */
	double monthly_fluxes_patch[12][NPERPATCHFLUXTYPES];

	/// Stores one flux value per month and flux type
	/** For the fluxes stored per pft for annual values */
	double monthly_fluxes_pft[12][NPERPFTFLUXTYPES];
};


/// Holds static functional parameters for a plant functional type (PFT). 
/** There should be one Pft object for each potentially occurring PFT. The same Pft object 
 *  may be referenced (via the pft member of the Individual object; see below) by different
 *  average individuals. Member functions are included for initialising SLA given leaf
 *  longevity, and for initialising sapling/regen characteristics (required for
 * population mode).
 */
class Pft {

	// MEMBER VARIABLES

public:
	/// id code (should be zero based and sequential, 0...npft-1)
	int id;
	/// name of PFT
	xtring name;
	/// life form (tree or grass)
	lifeformtype lifeform;
	/// leaf phenology (raingreen, summergreen, evergreen, rain+summergreen)
	phenologytype phenology;
	/// leaf physiognomy (needleleaf, broadleaf)
	leafphysiognomytype leafphysiognomy;
	/// growing degree sum on 5 degree base required for full leaf cover
	double phengdd5ramp;
	/// water stress threshold for leaf abscission (range 0-1; raingreen PFTs)
	double wscal_min;
	/// biochemical pathway for photosynthesis (C3 or C4)
	pathwaytype pathway;
	/// approximate low temperature limit for photosynthesis (deg C)
	double pstemp_min;
	/// approximate lower range of temperature optimum for photosynthesis (deg C)
	double pstemp_low;
	/// approximate upper range of temperature optimum for photosynthesis (deg C)
	double pstemp_high;
	/// maximum temperature limit for photosynthesis (deg C)
	double pstemp_max;
	/// non-water-stressed ratio of intercellular to ambient CO2 partial pressure
	double lambda_max;
	/// vegetation root profile 
	/** array containing fraction of roots in each soil layer, [0=upper layer] */
	double rootdist[NSOILLAYER];
	/// canopy conductance component not associated with photosynthesis (mm/s)
	double gmin;
	/// maximum evapotranspiration rate (mm/day)
	double emax;
	/// maintenance respiration coefficient (0-1)
	double respcoeff;
	
	/// minimum leaf C:N mass ratio allowed when nitrogen demand is determined
	double cton_leaf_min;
	/// maximum leaf C:N mass ratio	allowed when nitrogen demand is determined
	double cton_leaf_max;
	/// average leaf C:N mass ratio (between min and max)
	double cton_leaf_avr;
	/// average fine root C:N mass ratio (connected cton_leaf_avr)
	double cton_root_avr;
	/// maximum fine root C:N mass ratio (used when mass is negligible) 	
	double cton_root_max;
	/// average sapwood C:N mass ratio (connected cton_leaf_avr)
	double cton_sap_avr;
	/// maximum sapwood C:N mass ratio (used when mass is negligible) 
	double cton_sap_max;
	/// reference fine root C:N mass ratio	
	double cton_root;
	/// reference sapwood C:N mass ratio	
	double cton_sap;
	/// Maximum nitrogen (NH4+ and NO3- seperatly) uptake per fine root [kgN kgC-1 day-1]
	double nuptoroot;
	/// coefficient to compensate for vertical distribution of fine root on nitrogen uptake
	double nupscoeff;
	/// fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen longterm storage scalar
	double fnstorage;

	/// Michaelis-Menten kinetic parameters 
	/** Half saturation concentration for N uptake [kgN l-1] (Rothstein 2000) */
	double km_volume;
		
	/// fraction of NPP allocated to reproduction
	double reprfrac;
	/// annual leaf turnover as a proportion of leaf C biomass
	double turnover_leaf;
	/// annual fine root turnover as a proportion of fine root C biomass
	double turnover_root;
	/// annual sapwood turnover as a proportion of sapwood C biomass
	double turnover_sap;
	/// sapwood and heartwood density (kgC/m3)
	double wooddens;
	/// maximum tree crown area (m2)
	double crownarea_max;
	/// constant in allometry equations
	double k_allom1;
	/// constant in allometry equations
	double k_allom2;
	/// constant in allometry equations
	double k_allom3;
	/// constant in allometry equations
	double k_rp;
	/// tree leaf to sapwood area ratio
	double k_latosa;
	/// specific leaf area (m2/kgC)
	double sla;
	/// leaf longevity (years)
	double leaflong;
	/// leaf to root mass ratio under non-water-stressed conditions
	double ltor_max;
	/// litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// fire resistance (0-1)
	double fireresist;
	/// minimum forest-floor PAR level for growth (grasses) or establishment (trees)
	/** J/m2/day, individual and cohort modes */
	double parff_min;
	/** parameter capturing non-linearity in recruitment rate relative to
	 *  understorey growing conditions for trees (Fulton 1991) (individual and
	 *  cohort modes)
	 */
	double alphar;
	/// maximum sapling establishment rate (saplings/m2/year) (individual and cohort modes)
	double est_max;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect enabled (individual and cohort modes)
	 */
	double kest_repr;
	/// constant affecting amount of background establishment
	/** \see ifbgestab */
	double kest_bg;
	/** constant used in calculation of sapling establishment rate when spatial
	 *  mass effect disabled (individual and cohort modes)
	 */
	double kest_pres;
	/// expected longevity under non-stressed conditions (individual and cohort modes)
	double longevity;
	/// threshold growth efficiency for imposition of growth suppression mortality
	/** kgC/m2 leaf/year, individual and cohort modes */
	double greff_min;

	// Bioclimatic limits (all temperatures deg C)

	/// minimum 20-year coldest month mean temperature for survival
	double tcmin_surv;
	/// maximum 20-year coldest month mean temperature for establishment
	double tcmax_est;
	/// minimum degree day sum on 5 deg C base for establishment
	double gdd5min_est;
	/// minimum 20-year coldest month mean temperature for establishment
	double tcmin_est;
	/// minimum warmest month mean temperature for establishment
	double twmin_est;
	/// continentality parameter for boreal summergreen trees
	double twminusc;
	/// constant in equation for budburst chilling time requirement (Sykes et al 1996)
	double k_chilla;
	/// coefficient in equation for budburst chilling time requirement
	double k_chillb;
	/// exponent in equation for budburst chilling time requirement
	double k_chillk;
	/// array containing values for GDD0(c) given c=number of chill days
	/** Sykes et al 1996, Eqn 1
	 *  gdd0 has one element for each possible value for number of chill days
	 */
	double gdd0[Date::MAX_YEAR_LENGTH + 1];
	/// interception coefficient (unitless)
	double intc;

	/// Drought tolerance level (0 = very -> 1 = not at all) (unitless)
	/** Used to implement drought-limited establishment */
	double drought_tolerance;

	// bvoc
	
	/// aerodynamic conductance (m s-1)
	double ga;
	/// isoprene emission capacity (ug C g-1 h-1)
	double eps_iso;
	/// whether (1) or not (1) isoprene emissions show a seasonality
	bool seas_iso;
	/// monoterpene emission capacity (ug C g-1 h-1)
	double eps_mon;
	/// fraction of monoterpene production that goes into storage pool (-)
	double storfrac_mon;


	/// Sapling/regeneration characteristics (used only in population mode)
	/** For trees, on sapling individual basis (kgC); for grasses, on stand area basis,
	 *  kgC/m2 */
	struct {
		/// leaf C biomass
		double cmass_leaf;
		/// fine root C biomass
		double cmass_root;
		/// sapwood C biomass
		double cmass_sap;
		/// heartwood C biomass
		double cmass_heart;
	} regen;

	/// specifies type of landcover
	/** \see landcovertype */
	landcovertype landcover;

	/// Fraction of residue outtake at harvest.
	double res_outtake;
	/// Harvest efficiency.
	double harv_eff;
	/// Fraction of harvested products that goes into patchpft.harvested_products_slow
	double harvest_slow_frac;
	/// Yearly turnover fraction of patchpft.harvested_products_slow (goes to fluxes.acflux_harvest).
	double turnover_harv_prod;

	// MEMBER FUNCTIONS

public:

	Pft() {

		// Constructor (initialises array gdd0)

		int y;
		for (y=0; y<366; y++)
			gdd0[y] = -1.0; // value<0 signifies "unknown"; see function phenology()

		// guess2008 - DLE
		drought_tolerance = 0.0; // Default, means that the PFT will never be limited by drought.

		res_outtake = 0.0;
		harv_eff = 0.0;
		turnover_harv_prod = 1.0;	// default 1 year turnover time
	}

	/// Calculates SLA given leaf longevity
	void initsla() {

		// Reich et al 1992, Table 1 (includes conversion x2.0 from m2/kg_dry_weight to
		// m2/kgC)

		if (leafphysiognomy == BROADLEAF) {
			sla = 0.2 * pow(10.0, 2.41 - 0.38 * log10(12.0 * leaflong));
		}
		else if (leafphysiognomy == NEEDLELEAF) {
			sla = 0.2 * pow(10.0, 2.29 - 0.4 * log10(12.0 * leaflong));
		}
	}

	/// Calculates minimum leaf C:N ratio given leaf longevity
	void init_cton_min() {
		
		// Reich et al 1992, Table 1 (includes conversion x500 from mg/g_dry_weight to
		// kgN/kgC)

		if (leafphysiognomy == BROADLEAF)
			cton_leaf_min = 500.0 / pow(10.0, 1.75 - 0.33 * log10(12.0 * leaflong));
		else if (leafphysiognomy == NEEDLELEAF)
			cton_leaf_min = 500.0 / pow(10.0, 1.52 - 0.26 * log10(12.0 * leaflong));
	}

	void init_cton_limits() {

		// Fraction between min and max C:N ratio White et al. 2000
		double frac_mintomax = 2.78;

		// Fraction between leaf and root C:N ratio
		double frac_leaftoroot = 1.16; // Friend et al. 1997
		
		// Fraction between leaf and sap wood C:N ratio
		double frac_leaftosap = 6.9;   // Friend et al. 1997

		// Max leaf C:N ratio
		cton_leaf_max = cton_leaf_min * frac_mintomax;
		
		// Average leaf C:N ratio
		cton_leaf_avr = 1.0 / ((1.0 / cton_leaf_min + 1.0 / cton_leaf_max) / 2.0);

		// Average fine root C:N ratio
		cton_root_avr = cton_leaf_avr * frac_leaftoroot;

		// Maximum fine root C:N ratio
		cton_root_max = cton_leaf_min * frac_leaftoroot * frac_mintomax;

		// Average sap C:N ratio
		cton_sap_avr  = cton_leaf_avr * frac_leaftosap;

		// Maximum sap C:N ratio
		cton_sap_max  = cton_leaf_min * frac_leaftosap * frac_mintomax;

		if (lifeform == GRASS)
			respcoeff /= 2.0 * cton_root / (cton_root_avr + cton_leaf_min * frac_leaftoroot);
		else
			respcoeff /= cton_root / (cton_root_avr + cton_leaf_min * frac_leaftoroot) +
			             cton_sap  / (cton_sap_avr  + cton_leaf_min * frac_leaftosap);
	}

	/// Calculates coefficient to compensate for different vertical distribution of fine root on nitrogen uptake
	void init_nupscoeff() {
		
		// Fraction fine root in upper soil layer should have higher possibility for mineralized nitrogen uptake
		// Soil nitrogen profile is considered to have a exponential decline (Franzluebbers et al. 2009) giving 
		// an approximate advantage of 2 of having more roots in the upper soil layer
		const double upper_adv = 2.0;

		nupscoeff = rootdist[0] * upper_adv + rootdist[1];

	}

	/// Initialises sapling/regen characteristics in population mode following LPJF formulation
	void initregen() {

		// see function allometry in growth module.

		// Note: primary PFT parameters, including SLA, must be set before this
		//       function is called
	
		const double PI = 3.14159265;
		const double REGENLAI_TREE = 1.5;
		const double REGENLAI_GRASS = 0.001;
		const double SAPLINGHW = 0.2;

		if (lifeform == TREE) {

			// Tree sapling characteristics

			regen.cmass_leaf = pow(REGENLAI_TREE * k_allom1 * pow(1.0 + SAPLINGHW, k_rp) *
				pow(4.0 * sla / PI / k_latosa, k_rp * 0.5) / sla, 2.0 / (2.0 - k_rp));

			regen.cmass_sap = wooddens * k_allom2 * pow((1.0 + SAPLINGHW) *
				sqrt(4.0 * regen.cmass_leaf * sla / PI / k_latosa), k_allom3) *
				regen.cmass_leaf * sla / k_latosa;

			regen.cmass_heart = SAPLINGHW * regen.cmass_sap;
		}
		else if (lifeform == GRASS) {

			// Grass regeneration characteristics

			regen.cmass_leaf = REGENLAI_GRASS / sla;
		}

		regen.cmass_root = 1.0 / ltor_max * regen.cmass_leaf;
	}
};


/// A list of PFTs
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Pft objects. In general, there should be a single Pftlist object containing
 *  a single list of PFTs. Pft objects within the list are then referenced by the pft
 *  member of each Individual object.
 *
 * Functionality is inherited from the ListArray_id template type in the GUTIL
 * Library. Sequential Pft objects can be referenced as array elements by id:
 *
 *   Pftlist pftlist;
 *   ...
 *   for (i=0; i<npft; i++) {
 *     Pft& thispft=pftlist[i];
 *     // query or modify object thispft here
 *   }
 *
 * or by iteration through the linked list:
 *
 *   pftlist.firstobj();
 *   while (pftlist.isobj) {
 *     Pft& thispft=pftlist.getobj();
 *     // query or modify object thispft here
 *     pftlist.nextobj();
 *   }
 */
class Pftlist : public ListArray_id<Pft> {};

/// The one and only linked list of Pft objects	
extern Pftlist pftlist;


/// A vegetation individual.
/** In population mode this is the average individual of a PFT population; 
 *  in cohort mode: the average individual of a cohort;
 *  in individual mode: an individual plant. Each grass PFT is represented as a single
 *  individual in all modes. Individual objects are collected within list arrays of
 *  class Vegetation (defined below), of which there is one for each patch, and include
 *  a reference to their 'parent' Vegetation object. Use the createobj member function
 *  of class Vegetation to add new individuals.
 */
class Individual : public Serializable {

public:
	/// reference to Pft object containing static parameters for this individual
	Pft& pft;
	/// reference to Vegetation object to which this Individual belongs
	Vegetation& vegetation;
	/// id code (0-based, sequential)
	int id;
	/// leaf C biomass on modelled area basis (kgC/m2)
	double cmass_leaf;
	/// fine root C biomass on modelled area basis (kgC/m2)
	double cmass_root;
	/// sapwood C biomass on modelled area basis (kgC/m2)
	double cmass_sap;
	/// heartwood C biomass on modelled area basis (kgC/m2)
	double cmass_heart;
	/// C "debt" (retrospective storage) (kgC/m2)
	double cmass_debt;

	/// leaf N biomass on modelled area basis (kgC/m2)
	double nmass_leaf;
	/// root N biomass on modelled area basis (kgC/m2)	
	double nmass_root;
	/// sap N biomass on modelled area basis (kgC/m2)	
	double nmass_sap;
	/// heart N biomass on modelled area basis (kgC/m2)
	double nmass_heart;	

	/// foliar projective cover (FPC) under full leaf cover as fraction of modelled area
	double fpc;
	/// fraction of PAR absorbed by foliage over projective area today, taking account of leaf phenological state
	double fpar;
	/// average density of individuals over patch (indiv/m2)
	double densindiv;
	/// vegetation phenological state (fraction of potential leaf cover)
	double phen;
	/// annual sum of daily fractional leaf cover
	/** Equivalent number of days with full leaf cover
	 *  (population mode only; reset on expected coldest day of year)
	 */
	double aphen;
	/// annual number of days with full leaf cover) (raingreen PFTs only; reset on 1 January)
	int aphen_raingreen;

	/// Photosynthesis values for this individual under non-water-stress conditions
	PhotosynthesisResult photosynthesis;

	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;
		
	/// accumulated NPP over modelled area (kgC/m2/year); 
	/** annual NPP following call to growth module on last day of simulation year */
	double anpp;
	/// actual evapotranspiration over projected area (mm/day)
	double aet;
		// actual evapotranspiration over projected area (mm/day)
	/// leaf to root mass ratio
	double ltor;
	/// plant height (m)
	double height;
	/// plant crown area (m2)
	double crownarea;
	/// increment in fpc since last simulation year
	double deltafpc;
	// running sum (converted to annual mean) for wscal
	double wscal_mean;
	/// bole height, i.e. height above ground of bottom of crown cylinder (m)
	/** (individual and cohort modes only) */
	double boleht;
	/// patch-level lai for this individual or cohort (function fpar)
	double lai;
	/// patch-level lai for cohort in current vertical layer (function fpar)
	double lai_layer;
	/// individual leaf area index (individual and cohort modes only)
	double lai_indiv;
	/// growth efficiency (NPP/leaf area) for each of the last five simulation years (kgC/m2/yr)
	Historic<double, NYEARGREFF> greff_5;
	/// individual/cohort age (years)
	double age;
		// individual/cohort age (years)
	/*
	double mlai[12];
		// monthly LAI (including phenology component)
	*/
	/// FPAR assuming full leaf cover for all vegetation
	double fpar_leafon;
	/// LAI for current layer in canopy (cohort/individual mode; see function fpar)
	double lai_leafon_layer;
	/// non-water-stressed canopy conductance on FPC basis (mm/s)
	double gpterm;
	/// sub-daily version of the above variable (mm/s)
	std::vector<double> gpterms;
	/// interception associated with this individual today (patch basis)
	double intercep;

	/// accumulated mean fraction of potential leaf cover
	double phen_mean;
		
	/// whether individual subject to water stress
	bool wstress;

	/// leaf nitrogen that is photosyntetic active
	double nactive;
	/// Nitrogen extinction scalar
	/** Scalar to account for leaf nitrogen not following the optimal light 
	  * extinction, but is shallower.
	  */
	double nextin;
	/// long-term storage of labile nitrogen
	double nstore_longterm;
	/// storage of labile nitrogen
	double nstore_labile;
	/// daily total nitrogen demand
	double ndemand;
	/// fraction of individual nitrogen demand available for uptake
	double fnuptake;
	/// annual nitrogen uptake
	double anuptake;
	/// maximum size of nitrogen storage
	double max_n_storage;
	/// scales annual npp to maximum nitrogen storage
	double scale_n_storage;
	/// annual nitrogen limitation on vmax
	double avmaxnlim;
	/// annual optimal leaf C:N ratio
	double cton_leaf_aopt;
	/// annual average leaf C:N ratio
	double cton_leaf_aavr;
	/// plant mobile nitrogen status
	double cton_status;
	/// total carbon in compartments before growth
	double cmass_veg;
	/// total nitrogen in compartments before growth
	double nmass_veg;
	/// whether individual subject to nitrogen stress
	bool nstress;
	/// daily leaf nitrogen demand calculated from Vmax (kgN/m2)
	double leafndemand;
	/// daily root nitrogen demand based on leafndemand
	double rootndemand;
	/// daily sap wood nitrogen demand based on leafndemand
	double sapndemand;
	/// daily labile nitrogen demand based on npp
	double storendemand;
	/// leaf fraction of total nitrogen demand
	double leaffndemand;
	/// root fraction of total nitrogen demand
	double rootfndemand;
	/// sap fraction of total nitrogen demand
	double sapfndemand;
	/// store fraction of total nitrogen demand
	double storefndemand;
	/// daily leaf nitrogen demand over possible uptake (storage demand)
	double leafndemand_store;
	/// daily root nitrogen demand over possible uptake (storage demand)
	double rootndemand_store;
		
	/// Number of days with non-negligible phenology this month
	int nday_leafon;	
	// Whether this individual is truly alive. 
	/** Set to false for first year after the Individual object is created, then true. */
	bool alive; 


	// bvoc

	/// isoprene production (mg C m-2 d-1)
	double iso;
	/// monoterpene production (mg C m-2 d-1)
	double mon;
	/// monoterpene storage pool (mg C m-2)
	double monstor;
	/// isoprene seasonality factor (-)
	double fvocseas;

	// MEMBER FUNCTIONS

public:

	// Constructor function for objects of class Individual
	// Initialisation of certain member variables

	Individual(int i,Pft& p,Vegetation& v);

	void serialize(ArchiveStream& arch);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPFTFluxType flux_type, double value);

	/// Report a flux associated with this Individual
	/** Fluxes from 'new' Individuals (alive == false) will not be reported */
	void report_flux(Fluxes::PerPatchFluxType flux_type, double value);

	/// Reduce current biomass due to mortality and/or fire
	/** The removed biomass is put into litter pools and/or goes to fire fluxes.
	 *
	 *  \param mortality      fraction of Individual's biomass killed due to
	 *                        mortality (including fire)
	 *  \param mortality_fire fraction of Individual's biomass killed due to
	 *                        fire only
	 */
	void reduce_biomass(double mortality, double mortality_fire);

	/// Total storage of nitrogen
	double nstore() const {
		return nstore_longterm + nstore_labile;
	}

	/// Total carbon wood biomass
	double cmass_wood() const {
		return cmass_sap + cmass_heart - cmass_debt;
	}

	/// Total nitrogen wood biomass
	double nmass_wood() const {
		return nmass_sap + nmass_heart;
	}

	/// Current leaf C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_leaf(bool use_phen = true) const;

	/// Current fine root C:N ratio
	/**
	 *  \param use_phen Set to false if indiv.phen shouldn't be considered
	 *                  when calculating C:N ratio
	 */
	double cton_root(bool use_phen = true) const;

	/// Current sap C:N ratio
	double cton_sap() const;

	/// Gets the individual's Patchpft
	Patchpft& patchpft() const;

	/// Transfers the individual's biomass (C and N) to litter and harvest pools/fluxes
	/** 
	 *  \param harvest Set to true if some of the biomass should be harvested,
	 *                 harvest will be done according to the PFT's harvest efficiency
	 *                 and residue outtake.
	 */
	void kill(bool harvest = false);

	/// Gets the individual's daily cmass_leaf value
	double cmass_leaf_today() const;
	/// Gets the individual's daily cmass_root value
	double cmass_root_today() const;

	/// Gets the individual'sIndividual daily LAI value (patch-level)
	/** Based on total leaf area for whatever the individual represents
	 *  (individual, cohort, population), over the whole patch.
	 */
	double lai_today() const;

	/// Gets the individual's daily LAI value (individual-level)
	/** Based on the leaf area for the average individual and
	 *  the average individual's crown area.
	 */
	double lai_indiv_today() const;

	/// Gets the individual's daily fpc value
	double fpc_today() const;

	double ndemand_storage(double cton_leaf_opt);
};


/// The vegetation in a patch - a list of individuals
/** Functionality for building, maintaining, referencing and destroying a list array of
 *  Individual objects. A single Vegetation object is defined for each patch. A
 *  reference to the parent Patch object (defined below) is included as a member
 *  variable.
 *
 *  Functionality is inherited from the ListArray_idin1 template type in the GUTIL
 *  Library. Sequential Individual objects can be referenced as array elements by id,
 *  or by iteration through the linked list:
 *
 *    Vegetation vegetation
 *    ...
 *    vegetation.firstobj();
 *    while (vegetation.isobj) {
 *      Individual& thisindiv=vegetation.getobj();
 *      // query or modify object thisindiv here
 *      vegetation.nextobj();
 *    }
 */
class Vegetation : public ListArray_idin2<Individual,Pft,Vegetation>, public Serializable {

public:
	// MEMBER VARIABLES

	/// reference to parent Patch object
	Patch& patch;

	// MEMBER FUNCTIONS

	/// constructor (initialises member variable patch)
	Vegetation(Patch& p):patch(p) {};

	void serialize(ArchiveStream& arch);
};


/// Soiltype stores static parameters for soils and the snow pack.
/** One Soiltype object is defined for each Gridcell. State variables for soils
 *  are held by objects of class Soil, of which there is one for each patch
 *  (see below).
 */
class Soiltype {

	// MEMBER VARIABLES

public:

	/// available water holding capacity as fraction of soil volume
	double awc_frac;
	/// available water holding capacity of soil layers [0=upper layer] (mm)
	double awc[NSOILLAYER];

	/// coefficient in percolation calculation (K in Eqn 31, Haxeltine & Prentice 1996)		
	double perc_base;
	/// exponent in percolation calculation (=4 in Eqn 31, Haxeltine & Prentice 1996)
	double perc_exp;

	/// thermal diffusivity at 0% WHC (mm2/s)
	double thermdiff_0;
	/// thermal diffusivity at 15% WHC (mm2/s)
	double thermdiff_15;
	/// thermal diffusivity at 100% WHC (mm2/s)
	double thermdiff_100;

	/// wilting point of soil layers [0=upper layer] (mm) Cosby et al 1984
	double wp[NSOILLAYER];
	/// saturation point. Cosby et al 1984
	double wsats[NSOILLAYER];

	/// year at which to calculate equilibrium soil carbon
	int solvesom_end;
	/// year at which to begin documenting means for calculation of equilibrium soil carbon		
	int solvesom_begin;

	/// water holding capacity plus wilting point for whole soil volume
	double wtot; 

	// For CENTURY ...
	/// fraction of soil that is sand
	double sand_frac;
	/// fraction of soil that is clay
	double clay_frac;
	/// fraction of soil that is silt plus clay	
	double silt_frac;

	// MEMBER FUNCTIONS

public:

	/// Constructor: initialises certain member variables
	Soiltype() {

		solvesom_end = SOLVESOM_END;
		solvesom_begin = SOLVESOM_BEGIN;

	}

	/// Override the default SOM years with 70-80% of the spin-up period length
	void updateSolveSOMvalues(const int& nyrspinup) {
		
		solvesom_end = static_cast<int>(0.8 * nyrspinup);
		solvesom_begin = static_cast<int>(0.7 * nyrspinup);

	}
};

/// CENTURY SOIL POOL
class Sompool : public Serializable {

public:

	/// Constructor
	Sompool() {
		
		// Initialise pool
		
		cmass = 0.0;
		nmass = 0.0;
		cdec = 0.0;
		ndec = 0.0;
		ligcfrac = 0.0;
		delta_cmass = 0.0;
		delta_nmass = 0.0;
		fracremain = 0.0;
		litterme = 0.0;
		fireresist = 0.0;

		for (int m = 0; m < 12; m++) {
			mfracremain_mean[m] = 0.0;
		}
	}

	/// C mass in pool kgC/m2
	double cmass;
	/// Nitrogen mass in pool kgN/m2
	double nmass;
	/// (potential) decrease in C following decomposition today (kgC/m2)
	double cdec; 
	/// (potential) decrease in nitrogen following decomposition today (kgN/m2)
	double ndec; 
	/// daily change in carbon and nitrogen
	double delta_cmass,delta_nmass;
	/// lignin fractions
	double ligcfrac;
	/// fraction of pool remaining after decomposition
	double fracremain;
	/// nitrogen to carbon ratio
	double ntoc;

	// Fire
	/// soil litter moisture flammability threshold (fraction of AWC)
	double litterme;
	/// soil litter fire resistance (0-1)
	double fireresist;

	// Fast SOM spinup variables

	/// monthly mean fraction of carbon pool remaining after decomposition
	double mfracremain_mean[12];

	void serialize(ArchiveStream& arch);
};

/// This struct contains litter for solving Century SOM pools.
/** \see equilsom() */  
struct LitterSolveSOM : public Serializable {
	/// Constructs an empty result
	LitterSolveSOM() {
		clear();
	}

	/// Clears all members
	void clear() {
		for (int p = 0; p < NSOMPOOL; p++) {
			clitter[p] = 0.0;
			nlitter[p] = 0.0;
		}
	}

	/// Add litter
    void add_litter(double cvalue, double nvalue, int pool) {
		clitter[pool] += cvalue;
		nlitter[pool] += nvalue;
    }

	double get_clitter(int pool) {
		return clitter[pool];
	}
	double get_nlitter(int pool) {
		return nlitter[pool];
	}

	void serialize(ArchiveStream& arch);

private:
	/// Carbon litter
	double clitter[NSOMPOOL];
	
	/// Nitrogen litter
	double nlitter[NSOMPOOL];
};

/// Soil stores state variables for soils and the snow pack. 
/** Initialised by a call to initdrivers. One Soil object is defined for each patch. 
 *  A reference to the parent Patch object (defined below) is included as a member 
 *  variable. Soil static parameters are stored as objects of class Soiltype, of which 
 *  there is one for each grid cell. A reference to the Soiltype object holding the 
 *  static parameters for this soil is included as a member variable.
 */
class Soil : public Serializable {

	// MEMBER VARIABLES

public:
	/// reference to parent Patch object
	Patch& patch;
	/// reference to Soiltype object holding static parameters for this soil
	Soiltype& soiltype;
	/// water content of soil layers [0=upper layer] as fraction of available water holding capacity;
	double wcont[NSOILLAYER];
	/// DLE - the average wcont over the growing season, for each soil layer
	double awcont[NSOILLAYER];
	/// water content of sublayer of upper soil layer for which evaporation from the bare soil surface is possible 
	/** fraction of available water holding capacity */
	double wcont_evap;
	/// daily water content in upper soil layer for each day of year
	double dwcontupper[Date::MAX_YEAR_LENGTH];
	/// mean water content in upper soil layer for last month
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mwcontupper;
	/// stored snow as average over modelled area (mm rainfall equivalents)
	double snowpack;
	/// total runoff today (mm/day)
	double runoff;
	/// soil temperature today at 0.25 m depth (deg C)
	double temp;
	/// daily temperatures for the last month (deg C)
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double dtemp[31];
	/// mean soil temperature for the last month (deg C)
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mtemp;
	/** respiration response to today's soil temperature at 0.25 m depth
	 *  incorporating damping of Q10 due to temperature acclimation (Lloyd & Taylor 1994)
	 */
	double gtemp;
	/// soil organic matter (SOM) pool with c. 1000 yr turnover (kgC/m2)
	double cpool_slow;
	/// soil organic matter (SOM) pool with c. 33 yr turnover (kgC/m2)
	double cpool_fast;

	// Running sums (converted to long term means) maintained by SOM dynamics module

	/// mean annual litter decomposition (kgC/m2/yr)
	double decomp_litter_mean;
	/// mean value of decay constant for fast SOM fraction
	double k_soilfast_mean;
	/// mean value of decay constant for slow SOM fraction
	double k_soilslow_mean;
		

	// Parameters used by function soiltemp and updated monthly

	double alag, exp_alag;
		
	// Soil parameters stored during year for output at end of year (RCA)
	
	float mtemp_out[12];

	/// water content of soil layers [0=upper layer] as fraction of available water holding capacity
	double mwcont[12][NSOILLAYER];
	/// daily water content in lower soil layer for each day of year
	double dwcontlower[Date::MAX_YEAR_LENGTH];
	/// mean water content in lower soil layer for last month
	/** (valid only on last day of month following call to daily_accounting_patch) */
	double mwcontlower;

	/// rainfall and snowmelt today (mm)
	double rain_melt;
	/// upper limit for percolation (mm)
	double max_rain_melt;
	/// whether to percolate today
	bool percolate;

//////////////////////////////////////////////////////////////////////////////////
// CENTURY SOM pools and other variables

	Sompool sompool[NSOMPOOL];

	/// daily percolation (mm)
	double dperc;
	/// fraction of decayed organic nitrogen leached each day;
	double orgleachfrac;
	/// soil mineral nitrogen pool (kgN/m2)
	double nmass_avail;		
	/// soil nitrogen input (kgN/m2)
	double ninput;
	/// annual sum of nitrogen mineralisation
	double anmin;			
	/// annual sum of nitrogen immobilisation
	double animmob;			
	/// annual leaching from available nitrogen pool
	double aminleach;		
	/// annual leaching of organics from active nitrogen pool
	double aorgleach;		
	/// total annual nitrogen fixation 
	double anfix;
	/// calculated annual mean nitrogen fixation
	double anfix_calc;
	
	// Variables for fast spinup of SOM pools

	/// monthly fraction of available mineral nitrogen taken up
	double fnuptake_mean[12];
	/// monthly fraction of organic carbon/nitrogen leached
	double morgleach_mean[12];
	/// monthly fraction of available mineral nitrogen leached
	double mminleach_mean[12];
	/// annual nitrogen fixation
	double anfix_mean;

	// Solving Century SOM pools 

	/// years at which to begin documenting for calculation of Century equilibrium
	int solvesomcent_beginyr;
	/// years at which to end documentation and start calculation of Century equilibrium
	int solvesomcent_endyr;

	/// Cumulative litter pools for one year.
	LitterSolveSOM litterSolveSOM;
	std::vector<LitterSolveSOM> solvesom;

	/// stored nitrogen deposition in snowpack
	double snowpack_nmass;

	// MEMBER FUNCTIONS

public:
	/// constructor (initialises member variable patch)
	Soil(Patch& p,Soiltype& s):patch(p),soiltype(s) {
			initdrivers();
	}

	void initdrivers() {

		// Initialises certain member variables

		alag = 0.0;
		exp_alag = 1.0;
		cpool_slow = 0.0;
		cpool_fast = 0.0;
		decomp_litter_mean = 0.0;
		k_soilfast_mean = 0.0;
		k_soilslow_mean = 0.0;
		wcont[0] = 0.0;
		wcont[1] = 0.0;
		wcont_evap = 0.0;
		snowpack = 0.0;
		orgleachfrac = 0.0;


		// guess2008 - extra initialisation
		mwcontupper = 0.0;
		mwcontlower = 0.0;
		for (int mth=0; mth<12; mth++) {
			mwcont[mth][0] = 0.0;
			mwcont[mth][1] = 0.0;
			fnuptake_mean[mth] = 0.0;
			morgleach_mean[mth] = 0.0;
			mminleach_mean[mth] = 0.0;
		}

		for (int d = 0; d<Date::MAX_YEAR_LENGTH; d++) {
			dwcontupper[d] = 0.0;
			//			dwcontlower[d] = 0.0;
		}

		/////////////////////////////////////////////////////
		// Initialise CENTURY pools

		// Set initial CENTURY pool N:C ratios 
		// Parton et al 1993, Fig 4

		sompool[SOILMICRO].ntoc = 1.0 / 15.0;
		sompool[SURFHUMUS].ntoc = 1.0 / 15.0;
		sompool[SLOWSOM].ntoc = 1.0 / 20.0;
		sompool[SURFMICRO].ntoc = 1.0 / 20.0;

		// passive has a fixed value
		sompool[PASSIVESOM].ntoc = 1.0 / 9.0;

		nmass_avail = 0.0;
		ninput = 0.0;
		anmin = 0.0;			
		animmob = 0.0;		
		aminleach = 0.0;
		aorgleach = 0.0;
		anfix = 0.0;
		anfix_calc = 0.0;
		anfix_mean = 0.0;
		snowpack_nmass = 0.0;
		dperc = 0.0;

		solvesomcent_beginyr = (int)(SOLVESOMCENT_SPINBEGIN * (nyear_spinup - freenyears) + freenyears);
		solvesomcent_endyr   = (int)(SOLVESOMCENT_SPINEND   * (nyear_spinup - freenyears) + freenyears);
	}

	void serialize(ArchiveStream& arch);
};


/// State variables common to all individuals of a particular PFT in a particular patch
/** Used in individual and cohort modes only. */
class Patchpft : public Serializable {

	// MEMBER VARIABLES:

public:

	/// id code (equal to value of member variable id in corresponding Pft object)
	int id;
	/// reference to corresponding Pft object in PFT list
	Pft& pft;
	/// potential annual net assimilation (leaf-level net photosynthesis) at forest floor (kgC/m2/year)
	double anetps_ff;
	/// water stress parameter (0-1 range; 1=minimum stress)
	double wscal;
	/// running sum (converted to annual mean) for wscal
	double wscal_mean;
	/// potential annual net assimilation at forest floor averaged over establishment interval (kgC/m2/year)
	double anetps_ff_est;
	/// first-year value of anetps_ff_est
	double anetps_ff_est_initial;
	/// annual mean wscal averaged over establishment interval
	double wscal_mean_est;
	/// vegetation phenological state (fraction of potential leaf cover), updated daily		
	double phen;
	/// annual sum of daily fractional leaf cover 
	/** equivalent number of days with full leaf cover
	 *  (reset on expected coldest day of year)
	 */
	double aphen;
	/// whether PFT can establish in this patch under current conditions
	bool establish;
	/// running total for number of saplings of this PFT to establish (cohort mode)
	double nsapling;
	/// leaf-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_leaf;
	/// fine root-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_root;
	/// sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_sap;
	/// year's sapwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_sap_year;
	/// heartwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_heart;
	/// year's heartwood-derived litter for PFT on modelled area basis (kgC/m2)
	double litter_heart_year;
	/// litter derived from allocation to reproduction for PFT on modelled area basis (kgC/m2)
	double litter_repr;
	
	/// leaf-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_leaf;
	/// root-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_root;
	/// sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_sap;
	/// year's sapwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_sap_year;
	/// heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_heart;
	/// year's heartwood-derived nitrogen litter for PFT on modelled area basis (kgN/m2)
	double nmass_litter_heart_year;

	/// non-FPC-weighted canopy conductance value for PFT under water-stress conditions (mm/s)
	double gcbase;
		// non-FPC-weighted canopy conductance value for PFT under water-stress
		// conditions (mm/s)
	float gcbase_gross;
		// non-FPC-weighted canopy conductance value for PFT under water-stress
		// conditions including component not associated with photosynthesis (mm/s)	
	double gcbase_day;				// daily value of the above variable (mm/s)

	/// evapotranspirational "supply" function for this PFT today (mm/day)
	double wsupply;
	double wsupply_leafon;
	/// fractional uptake of water from each soil layer today
	double fwuptake[NSOILLAYER];
		
	/// whether water-stress conditions for this PFT
	bool wstress;
	/// daily version of the above variable
	bool wstress_day;

	/// carbon depository for long-lived products like wood
	double harvested_products_slow;	
	/// nitrogen depository for long-lived products like wood
	double harvested_products_slow_nmass; 

	// MEMBER FUNCTIONS:

	/// Constructor: initialises id, pft and data members
	Patchpft(int i,Pft& p):id(i),pft(p) {

		litter_leaf  = 0.0;
		litter_root  = 0.0;
		litter_sap   = 0.0;
		litter_sap_year = 0.0;
		litter_heart = 0.0;
		litter_heart_year = 0.0;
		litter_repr  = 0.0;

		nmass_litter_leaf  = 0.0;
		nmass_litter_root  = 0.0;
		nmass_litter_sap   = 0.0;
		nmass_litter_sap_year   = 0.0;
		nmass_litter_heart = 0.0;
		nmass_litter_heart_year = 0.0;

		wscal = 1.0;
		wscal_mean = 0.0;
		anetps_ff = 0.0;
		aphen = 0.0;

		harvested_products_slow = 0.0;
		harvested_products_slow_nmass = 0.0;
	}

	void serialize(ArchiveStream& arch);
};


/// Stores data for a patch. 
/** In cohort and individual modes, replicate patches are
 *  required in each stand to accomodate stochastic variation; in population mode there
 *  should be just one Patch object, representing average conditions for the entire
 *  stand. A reference to the parent Stand object (defined below) is included as a
 *  member variable.
 */
class Patch : public Serializable {

public:

	// MEMBER VARIABLES
	
	/// id code in range 0-npatch for patch
	int id;
	/// reference to parent Stand object
	Stand& stand;
	/// list array [0...npft-1] of Patchpft objects (initialised in constructor)
	ListArray_idin1<Patchpft,Pft> pft;
	/// vegetation for this patch
	Vegetation vegetation;
	/// soil for this patch
	Soil soil;
	/// fluxes for this patch
	Fluxes fluxes;
	/// FPAR at top of grass canopy today
	double fpar_grass;
	/// FPAR at soil surface today
	double fpar_ff;
	/// mean growing season PAR at top of grass canopy (J/m2/day)
	double par_grass_mean;
	/// number of days in growing season, estimated from mean vegetation leaf-on fraction 
	/** \see function fpar in canopy exchange module */
	int nday_growingseason;
	/// total patch FPC
	double fpc_total;
	/// whether patch was disturbed last year
	bool disturbed;
	/// patch age (years since last disturbance)
	int age;
	/// probability of fire this year
	double fireprob;
		
	/// DLE - the number of days over which wcont is averaged for this patch
	/** i.e. those days for which daily temp > 5.0 degC */
	int growingseasondays;


	// Variables used by new hydrology (Dieter Gerten 2002-07)
	
	/// interception by vegetation today on patch basis (mm)
	double intercep;
	/// annual sum of AET (mm/year)
	double aaet;
	/// annual sum of AET (mm/year) for each of the last five simulation years
	Historic<double, NYEARAAET> aaet_5;
	/// annual sum of soil evaporation (mm/year)
	double aevap;
	/// annual sum of interception (mm/year)
	double aintercep;
	/// annual sum of runoff (mm/year)
	double asurfrunoff;
	/// annual sum of runoff (mm/year)
	double adrainrunoff;
	/// annual sum of runoff (mm/year)
	double abaserunoff;
	/// annual sum of runoff (mm/year)
	double arunoff;
	/// annual sum of potential evapotranspiration (mm/year)
	double apet;
		
	/// equilibrium evapotranspiration today, deducting interception (mm)
	double eet_net_veg;
	
	/// transpirative demand for patch, patch vegetative area basis (mm/day)
	double wdemand;
	/// daily average of the above variable (mm/day)
	double wdemand_day;			
	/// transpirative demand for patch assuming full leaf cover today 
	/** mm/day, patch vegetative area basis	*/
	double wdemand_leafon;
	/// rescaling factor to account for spatial overlap between individuals/cohorts populations
	double fpc_rescale;
	
	/// monthly AET (mm/month)
	double maet[12];
	/// monthly soil evaporation (mm/month)
	double mevap[12];
	/// monthly interception (mm/month)
	double mintercep[12];
	/// monthly runoff (mm/month)
	double mrunoff[12];
	/// monthly PET (mm/month)
	double mpet[12];
		// monthly PET (mm/month)

	// For RCA coupling:
	float gc; // mean canopy conductance on patch basis (mm/s)

	/// daily nitrogen demand
	double ndemand;

	// MEMBER FUNCTIONS

	Patch(int i,Stand& s,Soiltype& st):
		id(i),stand(s),vegetation(*this),soil(*this,st),fluxes(*this) {

		// Constructor: initialises various members and builds list array
		// of Patchpft objects.

		pftlist.firstobj();
		while (pftlist.isobj) {
			pft.createobj(pftlist.getobj());
			pftlist.nextobj();
		}

		age = 0;
		disturbed = false;
		
		// guess2008 - initialise
		growingseasondays = 0;

		fireprob = 0.0;
		ndemand = 0.0;
	}

	void serialize(ArchiveStream& arch);
};

/// Container for variables common to individuals of a particular PFT in a stand.
/** Used in individual and cohort modes only
 */
class Standpft : public Serializable {

public:

	// MEMBER VARIABLES

	int id;
	Pft& pft;
	/// net C allocated to reproduction for this PFT in all patches of this stand this year (kgC/m2)
	double cmass_repr;
	/// maximum value of Patchpft::anetps_ff for this PFT in this stand so far in the simulation (kgC/m2/year)
	double anetps_ff_max;
	// non-FPAR-weighted value for canopy conductance component associated with
	// photosynthesis for PFT under non-water-stress conditions (mm/s)
	double gpterm;
	// sub-daily version of the above variable (mm/s)
	std::vector<double> gpterms;		
	// FPC sum for this PFT as average for stand (used by some versions of guessio.cpp)
	double fpc_total;
	/// Photosynthesis values for this PFT under non-water-stress conditions
	PhotosynthesisResult photosynthesis;
	/// sub-daily version of the above variable (NB: daily units)
	std::vector<PhotosynthesisResult> phots;

	/// Is this PFT allowed to grow in this stand?
	bool active;

	// Variables used by RCA-GUESS (Ben S 2007-06-14)
	
	float dlaiphen[366]; // Daily actual LAI for PFT, average for stand

	// MEMBER FUNCTIONS

	/// Constructor: initialises various data members
	Standpft(int i,Pft& p):id(i),pft(p) {

		anetps_ff_max = 0.0;
		active = !run_landcover;
	}

	void serialize(ArchiveStream& arch);
};


/// The stand class corresponds to a modelled area of a specific landcover type in a grid cell.
/** There may be several stands of the same landcover type (but with different settings).
 */
class Stand : public ListArray_idin2<Patch,Stand,Soiltype>, public Serializable {

public:

	// MEMBER VARIABLES

	/// list array [0...npft-1] of Standpft (initialised in constructor)
	ListArray_idin1<Standpft,Pft> pft;

	/// A number identifying this Stand within the grid cell
	int id;

	/// Seed for generating random numbers within this Stand
	/** The reason why Stand has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when using
	 *  different land cover types.
	 *
	 *  Randomness not associated with a specific stand, but rather a whole
	 *  grid cell should instead use the seed in the Gridcell class.
	 *
	 *  \see randfrac()
	 */
	long seed;

	/// reference to parent object
	Gridcell& gridcell;

	/// climate, insolation and CO2 for this stand
	Climate climate;
	
	/// type of landcover
	/** \see landcovertype
	 *  initialised in constructor
	 */
	landcovertype landcover;

	/// The year when this stand was created.
	/** Will typically be year zero unless running with dynamic
	 *  land cover.
	 *
	 *  Needed to set patchpft.anetps_ff_est_initial
	 */
	int first_year;

	/// soil water content 
	/** (sum of values from RCA short time step) */
	float wcont_sum[NSOILLAYER];

	/// prescribed wcont for this stand
	float wcont_prescribed[NSOILLAYER];

	float laiphen_grass[366]; // grass LAI on stand basis
	float laiphen_conifer[366]; // needleleaved tree LAI on stand basis
	float laiphen_broadleaf[366]; // broadleaved tree LAI on stand basis
	float laimax_conifer; // conifer LAI on stand basis
	float laimax_broadleaf; // broadleaved tree LAI on stand basis

	float gc_out[12]; // Daytime canopy conductance (mean among patches) by month (mm/s)

	// MEMBER FUNCTIONS

	/// Constructs a Stand
	/** \param i         The id for the stand within the grid cell
	 *  \param gc        The parent grid cell
	 *  \param landcover The type of landcover to use for this stand
	 */
	Stand(int i, Gridcell& gc,landcovertype landcover); 

	/// Gives the fraction of this Stand relative to the whole grid cell
	double get_gridcell_fraction() const;

	/// Gives the fraction of this Stand relative to its land cover type
	double get_landcover_fraction() const;

	/// Set the fraction of this Stand relative to its land cover type
	void set_landcover_fraction(double fraction);

	/// Returns the number of patches in this Stand
	unsigned int npatch() const { return nobj; }

	void serialize(ArchiveStream& arch);

private:

	/// Fraction of this stand relative to its landcover
	/** used by crop stands; initialized in constructor to 1,
	 *  set in landcover_init()
	 */
	double frac;
};



/// State variables common to all individuals of a particular PFT in a GRIDCELL.
class Gridcellpft : public Serializable {

public:

	// MEMBER VARIABLES

	/// A number identifying this object within its list array
	int id;

	/// A reference to the Pft object for this Gridcellpft
	Pft& pft;

	/// annual degree day sum above threshold damaging temperature
	/** used in calculation of heat stess mortality; Sitch et al 2000, Eqn 55
	 */
	double addtw;

	/// Michaelis-Menten kinetic parameters 
	/** Half saturation concentration for N uptake (Rothstein 2000, Macduff 2002) 
	 */
	double Km;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcellpft object
	/** \param i   The id for this object
	 *  \param p   A reference to the Pft for this Gridcellpft
	 */
	Gridcellpft(int i,Pft& p):id(i),pft(p) {
		addtw = 0.0;
		Km = 0.0;
	}

	void serialize(ArchiveStream& arch);
};


/// The Gridcell class corresponds to a modelled locality or grid cell.
/** Member variables include an object of type Climate (holding climate, insolation and
 *  CO2 data), a object of type Soiltype (holding soil static parameters) and a list
 *  array of Stand objects. Soil objects (holding soil state variables) are associated
 *  with patches, not gridcells. A separate Gridcell object must be declared for each modelled
 *  locality or grid cell.
 */
class Gridcell : public ListArray_idin2<Stand,Gridcell,landcovertype>, public Serializable {

public:

	// MEMBER VARIABLES

    /// soil static parameters for this grid cell
	Soiltype soiltype;

	/// The fractions of the different land cover types.
	/** landcoverfrac is read in from land cover input file or from
	 *  instruction file in getlandcover().
	 */
	double landcoverfrac[NLANDCOVERTYPES];

	// Required for linkage to RCA:
	double lon_cru,lat_cru; // CRU pixel (standard lon, lat below is RCA pixel)

	int id;
	bool island;

	/// Information about the land cover fractions
	std::auto_ptr<LandCoverInfo> landcoverinfo;
	
	// For coupling to RCA:
	Date date; // each gridcell keeps track of its own date in RCA-GUESS
	int warnct;

	// MEMBER FUNCTIONS

	/// The land cover fractions from the previous year
	/** Used to keep track of the changes when running with dynamic
	 *  land cover.
	 */
	double landcoverfrac_old[NLANDCOVERTYPES];

	/// Whether the land cover fractions changed for this grid cell this year
	/** \see landcover_dynamics
	 */
	bool LC_updated;

	/// list array [0...npft-1] of Gridcellpft (initialised in constructor)
	ListArray_idin1<Gridcellpft,Pft> pft;

	/// Seed for generating random numbers within this Gridcell
	/** The reason why Gridcell has its own seed, rather than using for instance
	 *  a single global seed is to make it easier to compare results when for
	 *  instance changing the order in which the simulation proceeds. It also
	 *  gets serialized together with the rest of the Gridcell state to make it
	 *  possible to get exactly identical results after a restart.
	 *
	 *  \see randfrac()
	 */
	long seed;

	// MEMBER FUNCTIONS

	/// Constructs a Gridcell object
	Gridcell() {
		landcovertype landcover;
		LC_updated = false;

		for(unsigned int p=0; p<pftlist.nobj; p++) {
			pft.createobj(pftlist[p]);
		}		

		memset(landcoverfrac, 0, sizeof(double) * NLANDCOVERTYPES);
		memset(landcoverfrac_old, 0, sizeof(double) * NLANDCOVERTYPES);

		if(!run_landcover) {
			landcover = NATURAL;
			createobj(*this,landcover);
			landcoverfrac[NATURAL] = 1.0;
		}

		island=false;
		warnct=0;

		seed = 12345678;
	}

	/// Longitude for this grid cell
	double get_lon() const;

	/// Latitude for this grid cell
	double get_lat() const;

	/// Set longitude and latitude for this grid cell
	void set_coordinates(double longitude, double latitude);

	void serialize(ArchiveStream& arch);

private:

	/// Longitude for this grid cell
	double lon;

	/// Latitude for this grid cell
	double lat;

};


#endif // LPJ_GUESS_GUESS_H

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// LPJF refers to the original FORTRAN implementation of LPJ as described by Sitch
//   et al 2000
// Delmas, R., Lacaux, J.P., Menaut, J.C., Abbadie, L., Le Roux, X., Helaa, G., Lobert, J., 1995. 
//   Nitrogen compound emission from biomass burning in tropical African Savanna FOS/DECAFE 1991 
//   experiment. Journal of Atmospheric Chemistry 22, 175-193.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical 
//   exploration of the relationships of soil moisture characteristic to the 
//   physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Franzlubbers, AJ & Stuedemann, JA 2009 Soil-profile organic carbon and total 
//   nitrogen during 12 years of pasture management in the Southern Piedmont USA. 
//   Agriculture Ecosystems & Environment, 129, 28-36.
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A 
//   process-based, terrestrial biosphere model of ecosystem dynamics 
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Fulton, MR 1991 Adult recruitment rate as a function of juvenile growth in size-
//   structured plant populations. Oikos 61: 102-105.
// Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
//   model based on ecophysiological constraints, resource availability, and
//   competition among plant functional types. Global Biogeochemical Cycles 10:
//   693-709
// Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
//   Functional Ecology 8: 315-323
// Macduff, JH, Humphreys, MO & Thomas, H 2002. Effects of a stay-green mutation on
//   plant nitrogen relations in Lolium perenne during N starvation and after 
//   defoliation. Annals of Botany, 89, 11-21.
// Monsi M & Saeki T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften und
//   seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany 14: 22-52
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W. 
//   & Kelly, R. 2010. ForCent model development and testing using the Enriched 
//   Background Isotope Study experiment. Journal of Geophysical 
//   Research-Biogeosciences, 115.
// Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
//   effects of climate change on forest landscapes. Ecological Modelling 65: 51-70.
// Reich, PB, Walters MB & Ellsworth DS 1992 Leaf Life-Span in Relation to Leaf,
//   Plant, and Stand Characteristics among Diverse Ecosystems. 
//   Ecological Monographs 62: 365-392.
// Sitch, S, Prentice IC, Smith, B & Other LPJ Consortium Members (2000) LPJ - a
//   coupled model of vegetation dynamics and the terrestrial carbon cycle. In:
//   Sitch, S. The Role of Vegetation Dynamics in the Control of Atmospheric CO2
//   Content, PhD Thesis, Lund University, Lund, Sweden.
// Sykes, MT, Prentice IC & Cramer W 1996 A bioclimatic model for the potential
//   distributions of north European tree species under present and future climates.
//   Journal of Biogeography 23: 209-233.
// White, M A, Thornton, P E, Running, S. & Nemani, R 2000 Parameterization and 
//   Sensitivity Analysis of the BIOME-BGC Terrestrial Ecosystem Model: Net Primary 
//   Production Controls. Earth Interactions, 4, 1-55.
