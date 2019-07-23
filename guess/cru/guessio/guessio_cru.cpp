///////////////////////////////////////////////////////////////////////////////////////
/// \file guessio_cru.cpp
/// \brief LPJ-GUESS input/output module with input from instruction script
///
/// This I/O module reads in CRU climate data in a customised binary format.
/// The binary files contain CRU half-degree global historical climate data
/// for 1901-2006.
///
/// \author Ben Smith
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"

#ifdef USE_CRU_IO

#include "guessio.h"

#include "driver.h"
#include "outputchannel.h"
#include "plib.h"
#include <stdio.h>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include "globalco2file.h"
#include "landcoverstrategies.h"
#include "RCAData.h"
#include "lamarquendep.h"
#include "rca_driver.h"


// guess2008 - header file for the CRU TS 3.0 data archives
#include "cru_1901_2006.h"
#include "cru_1901_2006misc.h"

///////////////////////////////////////////////////////////////////////////////////////
//
//                      SECTION: INPUT FROM INSTRUCTION SCRIPT
//
//  - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY - DO NOT MODIFY -
//
// The first section of this module is concerned with reading simulation settings and
// PFT parameters from an instruction script using functionality from the PLIB library.
// In general model users should not modify this section of the input/output module.
// New instructions (PLIB keywords) may be added (this would require addition of a
// declareitem call in function plib_declarations, and possibly some additional code in
// function plib_callback). However, it is probably preferable to use the "param"
// keyword feature, as this does not require any changes to this section of the module.
//
// Custom keywords may be included in the instruction script using syntax similar to
// the following examples:
//
//   param "co2" (num 340)
//   param "file_gridlist" (str "gridlist.txt")
//
// To retrieve the values associated with the "param" strings in the above examples,
// use the following function calls (may appear anywhere in this file; instruction
// script must have been read in first):
//
//   param["co2"].num
//   param["file_gridlist"].str
//
// Each "param" item can store EITHER a number (int or double) OR a string, but not
// both types of data. Function fail is called to terminate output if a "param" item
// with the specified identifier was not read in.
//
///////////////////////////////////////////////////////////////////////////////////////

/// Represents one custom "param" item
struct Paramtype {
	xtring name;
	xtring str;
	double num;
};

/// List for the custom parameters
/** Functionality for storing and retrieving custom "param" items from the instruction
 *  script.
 */
class Paramlist : public ListArray<Paramtype> {

public:
	/// Adds a parameter with a numeric value, overwriting if it already existed
	void addparam(xtring name,xtring value) {
		Paramtype* p = find(name);
		if (p == 0) {
			p = &createobj();
		}
		p->name=name.lower();
		p->str=value;
	}

	/// Adds a parameter with a string value, overwriting if it already existed
	void addparam(xtring name,double value) {
		Paramtype* p = find(name);
		if (p == 0) {
			p = &createobj();
		}
		p->name=name.lower();
		p->num=value;
	}

	/// Fetches a parameter from the list, aborts the program if it didn't exist
	Paramtype& operator[](xtring name) {
		Paramtype* param = find(name);
		
		if (param == 0) {
			fail("Paramlist::operator[]: parameter \"%s\" not found",(char*)name);
		}
		
		return *param;
	}

private:
	/// Tries to find the parameter in the list
	/** \returns 0 if it wasn't there. */
	Paramtype* find(xtring name) {
		name = name.lower();
		firstobj();
		while (isobj) {
			Paramtype& p=getobj();
			if (p.name==name) return &p;
			nextobj();
		}
		// nothing found
		return 0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////
// ENUM DECLARATIONS OF INTEGER CONSTANTS FOR PLIB INTERFACE

enum {BLOCK_GLOBAL,BLOCK_PFT,BLOCK_PARAM};
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_LANDCOVER,CB_PHENOLOGY,CB_LEAFPHYSIOGNOMY,
	CB_PATHWAY,	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM,CB_WATERUPTAKE};


///////////////////////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES WITH FILE SCOPE

Paramlist param;

xtring title; // Title for this run

int first_control_year,last_control_year;
	// first and last year to save GUESS->RCA driver data for a specified control period
bool ifcontrolperiod;
	// whether to save GUESS->RCA driver data for a specified control period

/// Nitrogen deposition time series to use (historic,rcp26,...)
xtring ndep_timeseries("historic");

double searchradius; // search radius to use when finding CRU data

/// Landcover fractions read from ins-file (% area).
int lc_fixed_frac[NLANDCOVERTYPES]={0};

/// Whether gridcell is divided into equal active landcover fractions.
bool equal_landcover_area;

Pft* ppft; // pointer to Pft object currently being assigned to

xtring paramname;
xtring strparam;
double numparam;
bool ifhelp=false;
bool includepft;


// guess2008 - Now declare the output file xtrings here
// Output file names ...
xtring outputdirectory;
xtring file_cmass,file_anpp,file_dens,file_lai,file_cflux,file_cpool,file_runoff;
xtring file_cmass_leaf, file_cmass_sap, file_cmass_heart, file_cmass_root, file_cmass_debt;
xtring file_mnpp,file_mlai,file_mgpp,file_mra,file_maet,file_mpet,file_mevap,file_mrunoff,file_mintercep,file_mrh;
xtring file_mnee,file_mwcont_upper,file_mwcont_lower;
xtring file_firert,file_speciesheights;

// Monthly output by stand
xtring file_mnpp_extra, file_mlai_extra, file_mgpp_extra, file_mra_extra, file_maet_extra, file_mpet_extra;
xtring file_mevap_extra, file_mrunoff_extra, file_mintercep_extra, file_mrh_extra, file_mnee_extra;
xtring file_miso_extra, file_mmon_extra;

// bvoc
xtring file_aiso,file_miso,file_amon,file_mmon;

// New outputs for RCA-GUESS
xtring file_dlai;
xtring file_driver;
xtring file_gc;
xtring file_landuse;
xtring file_par;
xtring file_prec;
xtring file_temp;
xtring file_temp_soil;
// Temporary files until the RCA spinup is forced by precipitation
// Unlike the regular mwcont_* files these are on stand basis so we
// can create the spinup archive based on them.
xtring file_mwcont_upper_extra;
xtring file_mwcont_lower_extra;

xtring file_cton_leaf, file_cton_veg, file_nsources, file_npool, file_nuptake, file_vmaxnlim, file_nflux, file_ngases;

void initsettings() {

	// Initialises global settings
	// Parameters not initialised here must be set in instruction script

	iffire=true;
	ifcalcsla=true;
	ifdisturb=false;
	ifcalcsla=false;
	ifcalccton=true;
	ifcdebt=false;
	distinterval=1.0e10;
	npatch=1;
	vegmode=COHORT;
	searchradius = 0;
	run_landcover = false;

	// guess2008 - initialise filenames here
	outputdirectory = "";
	file_cmass=file_anpp=file_lai=file_cflux=file_dens=file_runoff="";
	file_cmass_leaf=file_cmass_sap=file_cmass_heart=file_cmass_root=file_cmass_debt="";
	file_mnpp=file_mlai=file_maet=file_mpet=file_mevap=file_mrunoff=file_mintercep=file_mrh="";
	file_mgpp=file_mra=file_mnee=file_mwcont_upper=file_mwcont_lower="";
	file_cpool=file_firert=file_speciesheights="";

	// monthly output by stand
	file_mnpp_extra = file_mlai_extra = file_maet_extra = file_mpet_extra = file_mevap_extra = file_mrunoff_extra = ""; 
	file_mintercep_extra = file_mrh_extra = file_mgpp_extra = file_mra_extra = file_mnee_extra = "";
	file_miso_extra = file_mmon_extra = "";

	// bvoc
	file_aiso=file_miso=file_amon=file_mmon="";

	save_state = false;
	restart = false;

	file_cton_leaf=file_cton_veg=file_nsources=file_npool=file_nuptake=file_vmaxnlim=file_nflux=file_ngases="";
}

void initpft(Pft& pft,xtring& setname) {

	// Initialises a PFT object
	// Parameters not initialised here must be set in instruction script

	pft.name=setname;
	pft.lifeform=NOLIFEFORM;
	pft.phenology=NOPHENOLOGY;

	// Set bioclimatic limits so that PFT can establish and survive under all
	// conditions (may be overridden by settings in instruction script)

	pft.tcmin_surv=-1000.0;
	pft.tcmin_est=-1000.0;
	pft.tcmax_est=1000.0;
	pft.twmin_est=-1000.0;
	pft.gdd5min_est=1000.0;
	pft.twminusc=0.0;

	// Set chilling parameters so that no chilling period required for budburst

	pft.k_chilla=0.0;
	pft.k_chillb=0.0;
	pft.k_chillk=0.0;
}


///////////////////////////////////////////////////////////////////////////////////////
// INPUT FROM INSTRUCTION SCRIPT FILE
// The following code uses functionality from the PLIB library to process an
// instruction script (ins) file containing simulation settings and PFT parameters.
// Function readins() is called by the framework to initiate parsing of the script.
// Function printhelp() is called if GUESS is run with '-help' instead of an ins file
// name as a command line argument. Functions plib_declarations, plib_callback and
// plib_receivemessage comprise part of the interface to PLIB.

void plib_declarations(int id,xtring setname) {

	switch (id) {

	case BLOCK_GLOBAL:

		declareitem("title",&title,80,CB_NONE,"Title for run");
		declareitem("nyear_spinup",&nyear_spinup,1,10000,1,CB_NONE,"Number of simulation years to spinup for");
		declareitem("vegmode",&strparam,16,CB_VEGMODE,
			"Vegetation mode (\"INDIVIDUAL\", \"COHORT\", \"POPULATION\")");
		declareitem("ndep_timeseries", &ndep_timeseries, 10, CB_NONE, "Nitrogen deposition time series to use (historic, rcp26, rcp45, rcp60 or rcp85");
		declareitem("ifbgestab",&ifbgestab,1,CB_NONE,
			"Whether background establishment enabled (0,1)");
		declareitem("ifsme",&ifsme,1,CB_NONE,
			"Whether spatial mass effect enabled for establishment (0,1)");
		declareitem("ifstochmort",&ifstochmort,1,CB_NONE,
			"Whether mortality stochastic (0,1)");
		declareitem("ifstochestab",&ifstochestab,1,CB_NONE,
			"Whether establishment stochastic (0,1)");
		declareitem("estinterval",&estinterval,1,10,1,CB_NONE,
			"Interval for establishment of new cohorts (years)");
		declareitem("distinterval",&distinterval,1.0,1.0e10,1,CB_NONE,
			"Generic patch-destroying disturbance interval (years)");
		declareitem("iffire",&iffire,1,CB_NONE,
			"Whether fire enabled (0,1)");
		declareitem("ifdisturb",&ifdisturb,1,CB_NONE,
			"Whether generic patch-destroying disturbance enabled (0,1)");
		declareitem("ifcalcsla",&ifcalcsla,1,CB_NONE,
			"Whether SLA calculated from leaf longevity");
		declareitem("ifcalccton",&ifcalccton,1,CB_NONE,
			"Whether leaf C:N min calculated from leaf longevity");
		declareitem("ifcdebt",&ifcdebt,1,CB_NONE,
			"Whether to allow C storage");
		declareitem("npatch",&npatch,1,1000,1,CB_NONE,
			"Number of patches simulated");
		declareitem("patcharea",&patcharea,1.0,1.0e4,1,CB_NONE,
			"Patch area (m2)");
		declareitem("wateruptake", &strparam, 20, CB_WATERUPTAKE, 
			"Water uptake mode (\"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");

		declareitem("nrelocfrac",&nrelocfrac,0.0,1.0,1,CB_NONE,
			"Fractional nitrogen relocation from shed leaves & roots");
		declareitem("nfix_a",&nfix_a,0.0,0.4,1,CB_NONE,
			"first term in nitrogen fixation eqn");
		declareitem("nfix_b",&nfix_b,-10.0,10.,1,CB_NONE,
			"second term in nitrogen fixation eqn");

		declareitem("ifcentury",&ifcentury,1,CB_NONE,
			"Whether to use CENTURY SOM dynamics (default standard LPJ)");
		declareitem("ifnlim",&ifnlim,1,CB_NONE,
			"Whether plant growth limited by available nitrogen");
		declareitem("freenyears",&freenyears,0,1000,1,CB_NONE,
			"Number of years to spinup without nitrogen limitation");

		// Annual output variables
		declareitem("outputdirectory",&outputdirectory,300,CB_NONE,"Directory for the output files");
		declareitem("file_cmass",&file_cmass,300,CB_NONE,"C biomass output file");
		declareitem("file_cmass_leaf", &file_cmass_leaf, 300, CB_NONE, "C biomass output file");
		declareitem("file_cmass_sap", &file_cmass_sap, 300, CB_NONE, "C biomass output file");
		declareitem("file_cmass_heart", &file_cmass_heart, 300, CB_NONE, "C biomass output file");
		declareitem("file_cmass_root", &file_cmass_root, 300, CB_NONE, "C biomass output file");
		declareitem("file_cmass_debt", &file_cmass_debt, 300, CB_NONE, "C biomass output file");

		declareitem("file_anpp",&file_anpp,300,CB_NONE,"Annual NPP output file");
		declareitem("file_lai",&file_lai,300,CB_NONE,"LAI output file");
		declareitem("file_cflux",&file_cflux,300,CB_NONE,"C fluxes output file");
		declareitem("file_dens",&file_dens,300,CB_NONE,"Tree density output file");
		declareitem("file_cpool",&file_cpool,300,CB_NONE,"Soil C output file");
		declareitem("file_runoff",&file_runoff,300,CB_NONE,"Runoff output file");
		declareitem("file_firert",&file_firert,300,CB_NONE,"Fire retrun time output file");
		
		declareitem("file_cton_leaf",&file_cton_leaf,300,CB_NONE,"Mean leaf C:N output file");
		declareitem("file_cton_veg",&file_cton_veg,300,CB_NONE,"Mean vegetation C:N output file");
		declareitem("file_nsources",&file_nsources,300,CB_NONE,"Annual nitrogen sources output file");
		declareitem("file_npool",&file_npool,300,CB_NONE,"Soil nitrogen output file");
		declareitem("file_nuptake",&file_nuptake,300,CB_NONE,"Annual nitrogen uptake output file");
		declareitem("file_vmaxnlim",&file_vmaxnlim,300,CB_NONE,"Annual nitrogen limitation on vm output file");
		declareitem("file_nflux",&file_nflux,300,CB_NONE,"Annual nitrogen fluxes output file");
		declareitem("file_ngases",&file_ngases,300,CB_NONE,"Annual nitrogen gases output file");
		
		declareitem("file_speciesheights",&file_speciesheights,300,CB_NONE,"Mean species heights");

		// Monthly output variables
		declareitem("file_mnpp",&file_mnpp,300,CB_NONE,"Monthly NPP output file");
		declareitem("file_mlai",&file_mlai,300,CB_NONE,"Monthly LAI output file");
		declareitem("file_mgpp",&file_mgpp,300,CB_NONE,"Monthly GPP-LeafResp output file");
		declareitem("file_mra",&file_mra,300,CB_NONE,"Monthly autotrophic respiration output file");
		declareitem("file_maet",&file_maet,300,CB_NONE,"Monthly AET output file");
		declareitem("file_mpet",&file_mpet,300,CB_NONE,"Monthly PET output file");
		declareitem("file_mevap",&file_mevap,300,CB_NONE,"Monthly Evap output file");
		declareitem("file_mrunoff",&file_mrunoff,300,CB_NONE,"Monthly runoff output file");
		declareitem("file_mintercep",&file_mintercep,300,CB_NONE,"Monthly intercep output file");
		declareitem("file_mrh",&file_mrh,300,CB_NONE,"Monthly heterotrphic respiration output file");
		declareitem("file_mnee",&file_mnee,300,CB_NONE,"Monthly NEE output file");
		declareitem("file_mwcont_upper",&file_mwcont_upper,300,CB_NONE,"Monthly wcont_upper output file");
		declareitem("file_mwcont_lower",&file_mwcont_lower,300,CB_NONE,"Monthly wcont_lower output file");

		// Monthly output variables by stand
		declareitem("file_mnpp_extra", &file_mnpp_extra, 300, CB_NONE, "Monthly NPP output file");
		declareitem("file_mlai_extra", &file_mlai_extra, 300, CB_NONE, "Monthly LAI output file");
		declareitem("file_mgpp_extra", &file_mgpp_extra, 300, CB_NONE, "Monthly GPP-LeafResp output file");
		declareitem("file_mra_extra", &file_mra_extra, 300, CB_NONE, "Monthly autotrophic respiration output file");
		declareitem("file_maet_extra", &file_maet_extra, 300, CB_NONE, "Monthly AET output file");
		declareitem("file_mpet_extra", &file_mpet_extra, 300, CB_NONE, "Monthly PET output file");
		declareitem("file_mevap_extra", &file_mevap_extra, 300, CB_NONE, "Monthly Evap output file");
		declareitem("file_mrunoff_extra", &file_mrunoff_extra, 300, CB_NONE, "Monthly runoff output file");
		declareitem("file_mintercep_extra", &file_mintercep_extra, 300, CB_NONE, "Monthly intercep output file");
		declareitem("file_mrh_extra", &file_mrh_extra, 300, CB_NONE, "Monthly heterotrphic respiration output file");
		declareitem("file_mnee_extra", &file_mnee_extra, 300, CB_NONE, "Monthly NEE output file");
		declareitem("file_miso_extra", &file_miso_extra, 300, CB_NONE, "monthly isoprene flux output file");
		declareitem("file_mmon_extra", &file_mmon_extra, 300, CB_NONE, "monthly monoterpene flux output file");

		// bvoc
		declareitem("file_aiso",&file_aiso,300,CB_NONE,"annual isoprene flux output file");
		declareitem("file_miso",&file_miso,300,CB_NONE,"monthly isoprene flux output file");
		declareitem("file_amon",&file_amon,300,CB_NONE,"annual monoterpene flux output file");
		declareitem("file_mmon",&file_mmon,300,CB_NONE,"monthly monoterpene flux output file");

		declareitem("file_dlai", &file_dlai,300,CB_NONE,"Daily LAI output file");
		declareitem("file_driver", &file_driver,300,CB_NONE,"Driver data GUESS->RCA from control period");
		declareitem("file_gc", &file_gc,300,CB_NONE,"Monthly canopy conductance output file");
		declareitem("file_landuse", &file_landuse,300,CB_NONE,"Landuse output file");
		declareitem("file_par", &file_par,300,CB_NONE,"PAR output file");
		declareitem("file_prec", &file_prec,300,CB_NONE,"Precipitation output file");
		declareitem("file_temp", &file_temp,300,CB_NONE,"Temperature output file");
		declareitem("file_temp_soil", &file_temp_soil,300,CB_NONE,"Soil temperature output file");
		declareitem("file_mwcont_upper_extra",&file_mwcont_upper_extra,300,CB_NONE,"Monthly wcont_upper_extra output file");
		declareitem("file_mwcont_lower_extra",&file_mwcont_lower_extra,300,CB_NONE,"Monthly wcont_lower_extra output file");


		// guess2008 - new options
		declareitem("ifsmoothgreffmort",&ifsmoothgreffmort,1,CB_NONE,
			"Whether to vary mort_greff smoothly with growth efficiency (0,1)");
		declareitem("ifdroughtlimitedestab",&ifdroughtlimitedestab,1,CB_NONE,
			"Whether establishment drought limited (0,1)");
		declareitem("ifrainonwetdaysonly",&ifrainonwetdaysonly,1,CB_NONE,
			"Whether it rains on wet days only (1), or a little every day (0);");
		declareitem("searchradius", &searchradius, 0, 100, 1, CB_NONE,
			"If specified, CRU data will be searched for in a circle");

		// bvoc 
		declareitem("ifbvoc",&ifbvoc,1,CB_NONE,
			"Whether or not BVOC calculations are performed (0,1)");
		declareitem("run_landcover",&run_landcover,1,CB_NONE,"Landcover version");
		declareitem("run_urban",&run[URBAN],1,CB_NONE,"Whether urban land is to be simulated");
		declareitem("run_crop",&run[CROPLAND],1,CB_NONE,"Whether crop-land is to be simulated");
		declareitem("run_pasture",&run[PASTURE],1,CB_NONE,"Whether pasture is to be simulated");
		declareitem("run_forest",&run[FOREST],1,CB_NONE,"Whether managed forest is to be simulated");
		declareitem("run_natural",&run[NATURAL],1,CB_NONE,"Whether natural vegetation is to be simulated");
		declareitem("run_peatland",&run[PEATLAND],1,CB_NONE,"Whether peatland is to be simulated");
		declareitem("ifslowharvestpool",&ifslowharvestpool,1,CB_NONE,"If a slow harvested product pool is included in patchpft.");
		declareitem("lcfrac_fixed",&lcfrac_fixed,1,CB_NONE,"Whether static landcover fractions are set in the ins-file (0,1)");
		declareitem("equal_landcover_area",&equal_landcover_area,1,CB_NONE,"Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
		declareitem("lc_fixed_urban",&lc_fixed_frac[URBAN],0,100,1,CB_NONE,"% lc_fixed_urban");
		declareitem("lc_fixed_cropland",&lc_fixed_frac[CROPLAND],0,100,1,CB_NONE,"% lc_fixed_cropland");
		declareitem("lc_fixed_pasture",&lc_fixed_frac[PASTURE],0,100,1,CB_NONE,"% lc_fixed_pasture");
		declareitem("lc_fixed_forest",&lc_fixed_frac[FOREST],0,100,1,CB_NONE,"% lc_fixed_forest");
		declareitem("lc_fixed_natural",&lc_fixed_frac[NATURAL],0,100,1,CB_NONE,"% lc_fixed_natural");
		declareitem("lc_fixed_peatland",&lc_fixed_frac[PEATLAND],0,100,1,CB_NONE,"% lc_fixed_peatland");

		declareitem("state_path", &state_path, 300, CB_NONE, "State files directory (for restarting from, or saving state files)");
		declareitem("restart", &restart, 1, CB_NONE, "Whether to restart from state files");
		declareitem("save_state", &save_state, 1, CB_NONE, "Whether to save new state files");
		declareitem("state_year", &state_year, 1, 20000, 1, CB_NONE, "Save/restart year. Unspecified means just after spinup");

		declareitem("pft",BLOCK_PFT,CB_NONE,"Header for block defining PFT");
		declareitem("param",BLOCK_PARAM,CB_NONE,"Header for custom parameter block");
		callwhendone(CB_CHECKGLOBAL);


		break;
	
	case BLOCK_PFT:

		if (!ifhelp) {

			// Create and initialise a new Pft object and obtain a reference to it
			
			ppft=&pftlist.createobj();
			initpft(*ppft,setname);
			includepft=true;
		}

		declareitem("include",&includepft,1,CB_NONE,"Include PFT in analysis");
		declareitem("lifeform",&strparam,16,CB_LIFEFORM,
			"Lifeform (\"TREE\" or \"GRASS\")");
		declareitem("landcover",&strparam,16,CB_LANDCOVER,
			"Landcovertype (\"URBAN\", \"CROP\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
		declareitem("phenology",&strparam,16,CB_PHENOLOGY,
			"Phenology (\"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
		declareitem("leafphysiognomy",&strparam,16,CB_LEAFPHYSIOGNOMY,
			"Leaf physiognomy (\"NEEDLELEAF\" or \"BROADLEAF\")");
		declareitem("phengdd5ramp",&ppft->phengdd5ramp,0.0,1000.0,1,CB_NONE,
			"GDD on 5 deg C base to attain full leaf cover");
		declareitem("wscal_min",&ppft->wscal_min,0.0,1.0,1,CB_NONE,
			"Water stress threshold for leaf abscission (raingreen PFTs)");
		declareitem("pathway",&strparam,16,CB_PATHWAY,
			"Biochemical pathway (\"C3\" or \"C4\")");
		declareitem("pstemp_min",&ppft->pstemp_min,-50.0,50.0,1,CB_NONE,
			"Approximate low temp limit for photosynthesis (deg C)");
		declareitem("pstemp_low",&ppft->pstemp_low,-50.0,50.0,1,CB_NONE,
			"Approx lower range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_high",&ppft->pstemp_high,0.0,60.0,1,CB_NONE,
			"Approx higher range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_max",&ppft->pstemp_max,0.0,60.0,1,CB_NONE,
			"Maximum temperature limit for photosynthesis (deg C)");
		declareitem("lambda_max",&ppft->lambda_max,0.1,0.99,1,CB_NONE,
			"Non-water-stressed ratio of intercellular to ambient CO2 pp");
		declareitem("rootdist",ppft->rootdist,0.0,1.0,NSOILLAYER,CB_ROOTDIST,
			"Fraction of roots in each soil layer (first value=upper layer)");
		declareitem("gmin",&ppft->gmin,0.0,1.0,1,CB_NONE,
			"Canopy conductance not assoc with photosynthesis (mm/s)");
		declareitem("emax",&ppft->emax,0.0,50.0,1,CB_NONE,
			"Maximum evapotranspiration rate (mm/day)");
		// guess2008 - increased the upper limit to possible respcoeff values (was 1.2)
		declareitem("respcoeff",&ppft->respcoeff,0.0,3,1,CB_NONE,
			"Respiration coefficient (0-1)");

		declareitem("cton_root",&ppft->cton_root,1.0,1.0e4,1,CB_NONE,
			"Reference Fine root C:N mass ratio");
		declareitem("cton_sap",&ppft->cton_sap,1.0,1.0e4,1,CB_NONE,
			"Reference Sapwood C:N mass ratio");
		declareitem("nuptoroot",&ppft->nuptoroot,0.0,1.0,1,CB_NONE,
			"Maximum nitrogen uptake per fine root");
		declareitem("km_volume",&ppft->km_volume,0.0,10.0,1,CB_NONE,
			"Michaelis-Menten kinetic parameters for nitrogen uptake");
		declareitem("fnstorage",&ppft->fnstorage,0.0,10.0,1,CB_NONE,
			"fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen storage scalar");

		declareitem("reprfrac",&ppft->reprfrac,0.0,1.0,1,CB_NONE,
			"Fraction of NPP allocated to reproduction");
		declareitem("turnover_leaf",&ppft->turnover_leaf,0.0,1.0,1,CB_NONE,
			"Leaf turnover (fraction/year)");
		declareitem("turnover_root",&ppft->turnover_root,0.0,1.0,1,CB_NONE,
			"Fine root turnover (fraction/year)");
		declareitem("turnover_sap",&ppft->turnover_sap,0.0,1.0,1,CB_NONE,
			"Sapwood turnover (fraction/year)");
		declareitem("wooddens",&ppft->wooddens,10.0,1000.0,1,CB_NONE,
			"Sapwood and heartwood density (kgC/m3)");
		declareitem("crownarea_max",&ppft->crownarea_max,1.0,1000.0,1,CB_NONE,
			"Maximum tree crown area (m2)");
		declareitem("k_allom1",&ppft->k_allom1,10.0,1000.0,1,CB_NONE,
			"Constant in allometry equations");
		// guess2008 - changed lower limit for k_allom2 to 1 from 10. This is needed
		// for the shrub allometries.
		declareitem("k_allom2",&ppft->k_allom2,1.0,1.0e4,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_allom3",&ppft->k_allom3,0.1,1.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_rp",&ppft->k_rp,1.0,2.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_latosa",&ppft->k_latosa,100.0,1.0e5,1,CB_NONE,
			"Tree leaf to sapwood xs area ratio");
		declareitem("sla",&ppft->sla,1.0,1000.0,1,CB_NONE,
			"Specific leaf area (m2/kgC)");
		declareitem("cton_leaf_min",&ppft->cton_leaf_min,1.0,1.0e4,1,CB_NONE,
			"Minimum leaf C:N mass ratio");
		declareitem("ltor_max",&ppft->ltor_max,0.1,10.0,1,CB_NONE,
			"Non-water-stressed leaf:fine root mass ratio");
		declareitem("litterme",&ppft->litterme,0.0,1.0,1,CB_NONE,
			"Litter moisture flammability threshold (fraction of AWC)");
		declareitem("fireresist",&ppft->fireresist,0.0,1.0,1,CB_NONE,
			"Fire resistance (0-1)");
		declareitem("tcmin_surv",&ppft->tcmin_surv,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for survival (deg C)");
		declareitem("tcmin_est",&ppft->tcmin_est,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for establishment (deg C)");
		declareitem("tcmax_est",&ppft->tcmax_est,-50.0,1000.0,1,CB_NONE,
			"Max 20-year coldest month mean temp for establishment (deg C)");
		declareitem("twmin_est",&ppft->twmin_est,-1000.0,50.0,1,CB_NONE,
			"Min warmest month mean temp for establishment (deg C)");
		declareitem("twminusc",&ppft->twminusc,0,100,1,CB_NONE,
			"Stupid larch parameter");
		declareitem("gdd5min_est",&ppft->gdd5min_est,0.0,5000.0,1,CB_NONE,
			"Min GDD on 5 deg C base for establishment");
		declareitem("k_chilla",&ppft->k_chilla,0.0,5000.0,1,CB_NONE,
			"Constant in equation for budburst chilling time requirement");
		declareitem("k_chillb",&ppft->k_chillb,0.0,5000.0,1,CB_NONE,
			"Coefficient in equation for budburst chilling time requirement");
		declareitem("k_chillk",&ppft->k_chillk,0.0,1.0,1,CB_NONE,
			"Exponent in equation for budburst chilling time requirement");
		declareitem("parff_min",&ppft->parff_min,0.0,1.0e7,1,CB_NONE,
			"Min forest floor PAR for grass growth/tree estab (J/m2/day)");
		declareitem("alphar",&ppft->alphar,0.01,100.0,1,CB_NONE,
			"Shape parameter for recruitment-juv growth rate relationship");
		declareitem("est_max",&ppft->est_max,1.0e-4,1.0,1,CB_NONE,
			"Max sapling establishment rate (indiv/m2/year)");
		declareitem("kest_repr",&ppft->kest_repr,1.0,1000.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_bg",&ppft->kest_bg,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_pres",&ppft->kest_pres,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("longevity",&ppft->longevity,0.0,3000.0,1,CB_NONE,
			"Expected longevity under lifetime non-stressed conditions (yr)");
		declareitem("greff_min",&ppft->greff_min,0.0,1.0,1,CB_NONE,
			"Threshold for growth suppression mortality (kgC/m2 leaf/yr)");
		declareitem("leaflong",&ppft->leaflong,0.1,100.0,1,CB_NONE,
			"Leaf longevity (years)");
		declareitem("intc",&ppft->intc,0.0,1.0,1,CB_NONE,"Interception coefficient");
		
		// guess2008 - DLE
		declareitem("drought_tolerance",&ppft->drought_tolerance,0.0,1.0,1,CB_NONE,
			"Drought tolerance level (0 = very -> 1 = not at all) (unitless)");
		
		// bvoc
		declareitem("ga",&ppft->ga,0.0,1.0,1,CB_NONE,
			"aerodynamic conductance (m/s)");
		declareitem("eps_iso",&ppft->eps_iso,0.,100.,1,CB_NONE,
			"isoprene emission capacity (ug C g-1 h-1)");
		declareitem("seas_iso",&ppft->seas_iso,1,CB_NONE,
			"whether (1) or not (0) isoprene emissions show seasonality");
		declareitem("eps_mon",&ppft->eps_mon,0.,100.,1,CB_NONE,
			"monoterpene emission capacity (ug C g-1 h-1)");
		declareitem("storfrac_mon",&ppft->storfrac_mon,0.,1.,1,CB_NONE,
			"fraction of monoterpene production that goes into storage pool (-)");
		
		declareitem("harv_eff",&ppft->harv_eff,0.0,1.0,1,CB_NONE,"Harvest efficiency");
		declareitem("harvest_slow_frac",&ppft->harvest_slow_frac,0.0,1.0,1,CB_NONE,
			"Fraction of harvested products that goes into carbon depository for long-lived products like wood");
		declareitem("turnover_harv_prod",&ppft->turnover_harv_prod,0.0,1.0,1,CB_NONE,"Harvested products turnover (fraction/year)");
		declareitem("res_outtake",&ppft->res_outtake,0.0,1.0,1,CB_NONE,"Fraction of residue outtake at harvest");

		callwhendone(CB_CHECKPFT);
		
		break;

	case BLOCK_PARAM:

		paramname=setname;
		declareitem("str",&strparam,300,CB_STRPARAM,
			"String value for custom parameter");
		declareitem("num",&numparam,-1.0e38,1.0e38,1,CB_NUMPARAM,
			"Numerical value for custom parameter");
		
		break;
	}
}

void badins(xtring missing) {

	xtring message=(xtring)"Missing mandatory setting: "+missing;
	sendmessage("Error",message);
	plibabort();
}

void plib_callback(int callback) {

	xtring message;
	int i;
	double numval;

	switch (callback) {

	case CB_VEGMODE:
		if (strparam.upper()=="INDIVIDUAL") vegmode=INDIVIDUAL;
		else if (strparam.upper()=="COHORT") vegmode=COHORT;
		else if (strparam.upper()=="POPULATION") vegmode=POPULATION;
		else {
			sendmessage("Error",
				"Unknown vegetation mode (valid types: \"INDIVIDUAL\",\"COHORT\", \"POPULATION\")");
			plibabort();
		}
		break;
	case CB_WATERUPTAKE:
		if (strparam.upper() == "WCONT") wateruptake = WR_WCONT;
		else if (strparam.upper() == "ROOTDIST") wateruptake = WR_ROOTDIST;
		else if (strparam.upper() == "SMART") wateruptake = WR_SMART;
		else if (strparam.upper() == "SPECIESSPECIFIC") wateruptake = WR_SPECIESSPECIFIC;
		else {
			sendmessage("Error",
				"Unknown water uptake mode (valid types: \"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");
		}
		break;
	case CB_LIFEFORM:
		if (strparam.upper()=="TREE") ppft->lifeform=TREE;
		else if (strparam.upper()=="GRASS") ppft->lifeform=GRASS;
		else {
			sendmessage("Error",
				"Unknown lifeform type (valid types: \"TREE\", \"GRASS\")");
			plibabort();
		}
		break;
	case CB_LANDCOVER:
		if (strparam.upper()=="NATURAL") ppft->landcover=NATURAL;
		else if (strparam.upper()=="URBAN") ppft->landcover=URBAN;
		else if (strparam.upper()=="CROPLAND") ppft->landcover=CROPLAND;
		else if (strparam.upper()=="PASTURE") ppft->landcover=PASTURE;
		else if (strparam.upper()=="FOREST") ppft->landcover=FOREST;			
		else if (strparam.upper()=="PEATLAND") ppft->landcover=PEATLAND;
		else {
			sendmessage("Error",
				"Unknown landcover type (valid types: \"URBAN\", \"CROPLAND\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
			plibabort();
		}
		break;
	case CB_PHENOLOGY:
		if (strparam.upper()=="SUMMERGREEN") ppft->phenology=SUMMERGREEN;
		else if (strparam.upper()=="RAINGREEN") ppft->phenology=RAINGREEN;
		else if (strparam.upper()=="EVERGREEN") ppft->phenology=EVERGREEN;
		else if (strparam.upper()=="ANY") ppft->phenology=ANY;
		else {
			sendmessage("Error",
				"Unknown phenology type\n  (valid types: \"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
			plibabort();
		}
		break;
	case CB_LEAFPHYSIOGNOMY:
		if (strparam.upper()=="NEEDLELEAF") ppft->leafphysiognomy=NEEDLELEAF;
		else if (strparam.upper()=="BROADLEAF") ppft->leafphysiognomy=BROADLEAF;
		else {
			sendmessage("Error",
				"Unknown leaf physiognomy (valid types: \"NEEDLELEAF\", \"BROADLEAF\")");
			plibabort();
		}
		break;
	case CB_PATHWAY:
		if (strparam.upper()=="C3") ppft->pathway=C3;
		else if (strparam.upper()=="C4") ppft->pathway=C4;
		else {
			sendmessage("Error",
				"Unknown pathway type\n  (valid types: \"C3\" or \"C4\")");
			plibabort();
		}
		break;
	case CB_ROOTDIST:
		numval=0.0;
		for (i=0;i<NSOILLAYER;i++) numval+=ppft->rootdist[i];
		if (numval<0.99 || numval>1.01) {
			sendmessage("Error","Specified root fractions do not sum to 1.0");
			plibabort();
		}
		ppft->rootdist[NSOILLAYER-1]+=1.0-numval;
		break;
	case CB_STRPARAM:
		param.addparam(paramname,strparam);
		break;
	case CB_NUMPARAM:
		param.addparam(paramname,numparam);
		break;
	case CB_CHECKGLOBAL:
		if (!itemparsed("title")) badins("title");
		if (!itemparsed("nyear_spinup")) badins("nyear_spinup");
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcalccton")) badins("ifcalccton");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");
		if (!itemparsed("wateruptake")) badins("wateruptake");

		if (!itemparsed("nrelocfrac")) badins("nrelocfrac");
		if (!itemparsed("nfix_a")) badins("nfix_a");
		if (!itemparsed("nfix_b")) badins("nfix_b");

		if (!itemparsed("ifcentury")) badins("ifcentury");
		if (!itemparsed("ifnlim")) badins("ifnlim");
		if (!itemparsed("freenyears")) badins("freenyears");

		if (nyear_spinup <= freenyears) {
			sendmessage("Error", "freenyears must be smaller than nyear_spinup");
			plibabort();
		}

		if (!itemparsed("outputdirectory")) badins("outputdirectory");
		if (!itemparsed("ifsmoothgreffmort")) badins("ifsmoothgreffmort");
		if (!itemparsed("ifdroughtlimitedestab")) badins("ifdroughtlimitedestab");
		if (!itemparsed("ifrainonwetdaysonly")) badins("ifrainonwetdaysonly");
		// bvoc
		if (!itemparsed("ifbvoc")) badins("ifbvoc");

		if (!itemparsed("run_landcover")) badins("run_landcover");
		if (run_landcover) {
			if (!itemparsed("lcfrac_fixed")) badins("lcfrac_fixed");
			if (!itemparsed("equal_landcover_area")) badins("equal_landcover_area");
			if (!itemparsed("lc_fixed_urban")) badins("lc_fixed_urban");
			if (!itemparsed("lc_fixed_cropland")) badins("lc_fixed_cropland");
			if (!itemparsed("lc_fixed_pasture")) badins("lc_fixed_pasture");
			if (!itemparsed("lc_fixed_forest")) badins("lc_fixed_forest");
			if (!itemparsed("lc_fixed_natural")) badins("lc_fixed_natural");
			if (!itemparsed("lc_fixed_peatland")) badins("lc_fixed_peatland");
			if (!itemparsed("run_natural")) badins("run_natural");
			if (!itemparsed("run_crop")) badins("run_crop");
			if (!itemparsed("run_forest")) badins("run_forest");
			if (!itemparsed("run_urban")) badins("run_urban");
			if (!itemparsed("run_pasture")) badins("run_pasture");
			if (!itemparsed("ifslowharvestpool")) badins("ifslowharvestpool");
		}

		if (!itemparsed("pft")) badins("pft");
		if (vegmode==COHORT || vegmode==INDIVIDUAL) {
			if (!itemparsed("ifbgestab")) badins("ifbgestab");
			if (!itemparsed("ifsme")) badins("ifsme");
			if (!itemparsed("ifstochmort")) badins("ifstochmort");
			if (!itemparsed("ifstochestab")) badins("ifstochestab");
			if (itemparsed("ifdisturb") && !itemparsed("distinterval"))
				badins("distinterval");
			if (!itemparsed("npatch")) badins("npatch");
			if (!itemparsed("patcharea")) badins("patcharea");
			if (!itemparsed("estinterval")) badins("estinterval");
		}
		else if (vegmode==POPULATION && npatch!=1) {
			sendmessage("Information",
				"Value specified for npatch ignored in population mode");
			npatch=1;
		}

		if (save_state && restart) {
			sendmessage("Error",
			            "Can't save state and restart at the same time");
			plibabort();
		}

		if (!itemparsed("state_year")) {
			state_year = nyear_spinup;
		}

		if (state_path == "" && (save_state || restart)) {
			badins("state_path");
		}

		break;
	case CB_CHECKPFT:
		if (!itemparsed("lifeform")) badins("lifeform");
		if (!itemparsed("phenology")) badins("phenology");
		if (ppft->phenology==SUMMERGREEN || ppft->phenology==ANY)
			if (!itemparsed("phengdd5ramp")) badins("phengdd5ramp");
		if (ppft->phenology==RAINGREEN || ppft->phenology==ANY)
			if (!itemparsed("wscal_min")) badins("wscal_min");
		if (!itemparsed("pathway")) badins("pathway");
		if (!itemparsed("pstemp_min")) badins("pstemp_min");
		if (!itemparsed("pstemp_low")) badins("pstemp_low");
		if (!itemparsed("pstemp_high")) badins("pstemp_high");
		if (!itemparsed("pstemp_max")) badins("pstemp_max");
		if (!itemparsed("lambda_max")) badins("lambda_max");
		if (!itemparsed("rootdist")) badins("rootdist");
		if (!itemparsed("gmin")) badins("gmin");
		if (!itemparsed("emax")) badins("emax");
		if (!itemparsed("respcoeff")) badins("respcoeff");
		if (!itemparsed("sla") && !ifcalcsla) badins("sla");
		if (!itemparsed("cton_leaf_min") && !ifcalccton) badins("cton_leaf_min");

		if (!itemparsed("cton_root")) badins("cton_root");
		if (!itemparsed("nuptoroot")) badins("nuptoroot");
		if (!itemparsed("km_volume")) badins("km_volume");
		if (!itemparsed("fnstorage")) badins("fnstorage");

		if (!itemparsed("reprfrac")) badins("reprfrac");
		if (!itemparsed("turnover_leaf")) badins("turnover_leaf");
		if (!itemparsed("turnover_root")) badins("turnover_root");
		if (!itemparsed("ltor_max")) badins("ltor_max");
		if (!itemparsed("intc")) badins("intc");
		if (!itemparsed("leafphysiognomy")) badins("leafphysiognomy");

		if (run_landcover)
		{
			if (!itemparsed("landcover")) badins("landcover");
			if (!itemparsed("turnover_harv_prod")) badins("turnover_harv_prod");
			if (!itemparsed("harvest_slow_frac")) badins("harvest_slow_frac");
			if (!itemparsed("harv_eff")) badins("harv_eff");
			if (!itemparsed("res_outtake")) badins("res_outtake");
		}

		// guess2008 - DLE
		if (!itemparsed("drought_tolerance")) badins("drought_tolerance");

		// bvoc
		if(ifbvoc){
		  if (!itemparsed("ga")) badins("ga");
		  if (!itemparsed("eps_iso")) badins("eps_iso");
		  if (!itemparsed("seas_iso")) badins("seas_iso");
		  if (!itemparsed("eps_mon")) badins("eps_mon");
		  if (!itemparsed("storfrac_mon")) badins("storfrac_mon");
		}

		if (ppft->lifeform==TREE) {
			if (!itemparsed("cton_sap")) badins("cton_sap");
			if (!itemparsed("turnover_sap")) badins("turnover_sap");
			if (!itemparsed("wooddens")) badins("wooddens");
			if (!itemparsed("crownarea_max")) badins("crownarea_max");
			if (!itemparsed("k_allom1")) badins("k_allom1");
			if (!itemparsed("k_allom2")) badins("k_allom2");
			if (!itemparsed("k_allom3")) badins("k_allom3");
			if (!itemparsed("k_rp")) badins("k_rp");
			if (!itemparsed("k_latosa")) badins("k_latosa");
			if (vegmode==COHORT || vegmode==INDIVIDUAL) {
				if (!itemparsed("kest_repr")) badins("kest_repr");
				if (!itemparsed("kest_bg")) badins("kest_bg");
				if (!itemparsed("kest_pres")) badins("kest_pres");
				if (!itemparsed("longevity")) badins("longevity");
				if (!itemparsed("greff_min")) badins("greff_min");		
				if (!itemparsed("alphar")) badins("alphar");
				if (!itemparsed("est_max")) badins("est_max");
			}
		}
		if (iffire) {
			if (!itemparsed("litterme")) badins("litterme");
			if (!itemparsed("fireresist")) badins("fireresist");
		}
		if (ifcalcsla) {
			if (!itemparsed("leaflong")) {
				sendmessage("Error",
					"Value required for leaflong when ifcalcsla enabled");
				plibabort();
			}
			if (itemparsed("sla"))
				sendmessage("Warning",
				"Specified sla value not used when ifcalcsla enabled");

			// Calculate SLA
			ppft->initsla();
		}
		if (vegmode==COHORT || vegmode==INDIVIDUAL) {
			if (!itemparsed("parff_min")) badins("parff_min");	
		}

		if (ifcalccton) {
			if (!itemparsed("leaflong")) {
				sendmessage("Error",
					"Value required for leaflong when ifcalccton enabled");
				plibabort();
			}
			if (itemparsed("cton_leaf_min"))
				sendmessage("Warning",
				"Specified cton_leaf_min value not used when ifcalccton enabled");
		
			// Calculate leaf C:N ratio minimum
			ppft->init_cton_min();
		}

		// Calculate C:N ratio limits
		ppft->init_cton_limits();

		// Calculate nitrogen uptake strength dependency on root distribution
		ppft->init_nupscoeff();

		// Calculate regeneration characteristics for population mode
		ppft->initregen();

		ppft->id=npft++;
			// VERY IMPORTANT (cannot rely on internal id counter of collection class)

		//	delete unused pft:s from pftlist

		if (ppft->landcover!=NATURAL) {
			if (!run_landcover || !run[ppft->landcover])
				includepft=0;
		}
		else if (run_landcover && !run[NATURAL]) {
			if (ppft->landcover==NATURAL)
				includepft=0;
		}

		// If "include 0", remove this PFT from list, and set id to correct value

		if (!includepft) {
			pftlist.killobj();
			npft--;
		}

		break;
	}
}

void plib_receivemessage(xtring text) {

	// Output of messages to user sent by PLIB

	dprintf((char*)text);
}

bool readins(xtring filename) {

	// DESCRIPTION
	// Uses PLIB library functions to read instructions from file specified by
	// 'filename', returning true if file could be successfully opened and read, and
	// no errors were encountered.

	// OUTPUT PARAMETERS
	// pftlist  = initialised list array of PFT parameters

	// Initialise PFT count
	npft=0;

	// Initialise certain parameters
	initsettings();
	param.killall();

	// Call PLIB
	return plib(filename);
}

void printhelp() {

	// Calls PLIB to output help text

	ifhelp=true;
	plibhelp();
	ifhelp=false;
}


///////////////////////////////////////////////////////////////////////////////////////
//
//             SECTION: INPUT OF ENVIRONMENTAL DRIVING DATA FOR SIMULATION
//                            OUTPUT OF SIMULATION RESULTS
//
// In general it is the responsibility of the user of the model to provide code for
// this section of the input/output module. The following functions are called by the
// framework at various stages of the simulation and should contain appropriate code:
//
// void initio(const xtring& insfilename)
//   Initialises input/output (e.g. opening files), sets values for the global
//   simulation parameter variables (currently vegmode, npatch, patcharea,
//   ifbgestab, ifsme, ifstochestab, ifstochmort, iffire, estinterval,
//   npft), initialises pftlist (the one and only list of PFTs and their static
//   parameters for this run of the model). Normally all of the above parameters,
//   and possibly others, are read from the ins file (see above). Function readins
//   should be called to input settings from the ins file.
//
// bool getgridcell(Gridcell& gridcell)
//   Obtains coordinates and soil static parameters for the next grid cell to
//   simulate. The function should return false if no grid cells remain to be simulated,
//   otherwise true. Currently the following member variables of gridcell should be
//   initialised: longitude, latitude and climate.instype; the following members of
//   member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
//   thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
//   code (Sitch et al 2000) by a call to function soilparameters in the driver
//   module (driver.cpp):
//
//   soilparameters(gridcell.soiltype,soilcode);
//
//   If the model is to be driven by quasi-daily values of the climate variables
//   derived from monthly means, this function may be the appropriate place to
//   perform the required interpolations. The utility functions interp_monthly_means
//   and interp_monthly_totals in driver.cpp may be called for this purpose.
//
// bool getclimate(Gridcell& gridcell)
//   Obtains climate data (including atmospheric CO2 and insolation) for this day.
//   The function should return false if the simulation is complete for this grid cell,
//   otherwise true. This will normally require querying the year and day member
//   variables of the global class object date:
//
//   if (date.day==0 && date.year==nyear_spinup) return false;
//   // else
//   return true;
//
//   Currently the following member variables of the climate member of gridcell must be
//   initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
//   values of the climate variables derived from monthly means, this day's values
//   will presumably be extracted from arrays containing the interpolated daily
//   values (see function getgridcell):
//
//   gridcell.climate.temp=dtemp[date.day];
//   gridcell.climate.prec=dprec[date.day];
//   gridcell.climate.insol=dsun[date.day];
//
//   Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
//   BVOC:
//   gridcell.climate.dtr=ddtr[date.day]; 
//
//   If model is run in diurnal mode, which requires appropriate climate forcing data, 
//   additional members of the climate must be initialised: temps, insols. Both of the
//   variables must be of type std::vector. The length of these vectors should be equal
//   to value of date.subdaily which also needs to be set either in getclimate or 
//   getgridcell functions. date.subdaily is a number of sub-daily period in a single 
//   day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
//   diurnal mode.
//
// void outannual(Gridcell& gridcell)
//   Called at the end of the last day of each simulation year to permit output of
//   model results.
//
// termio()
//   Called after simulation is complete for all gridcells to allow memory deallocation,
//   closing of files or other cleanup functions.
//
// printhelp()
//   Prints out information about all available ins file parameters. Is typically
//   called when the user starts the program with the -help option.
//
///////////////////////////////////////////////////////////////////////////////////////

struct Coord {

	// Type for storing grid cell longitude, latitude and description text

	int id;
	double lon;
	double lat;
	xtring descrip;

};


ListArray_id<Coord> gridlist;
	// Will maintain a list of Coord objectsc ontaining coordinates
	// of the grid cells to simulate

int ngridcell; // the number of grid cells to simulate

class Spinup_data {

	// Class for management of climate data for spinup
	// (derived from first few years of historical climate data)

private:
	int nyear;
	int thisyear;
	double* data;
	bool havedata;

	// guess2008 - this array holds the climatology for the spinup period
	double dataclim[12];


public:
	Spinup_data(int nyear_loc) {
		nyear=nyear_loc;
		havedata=false;
		data=new double[nyear*12];
		if (!data) fail("Spinup_data::Spinup_data: out of memory");
		thisyear=0;
		havedata=true;
		reset_clim(); // guess2008
	}

	~Spinup_data() {
		if (havedata) delete[] data;
	}

	double& operator[](int month) {

		return data[thisyear*12+month];
	}

	void nextyear() {
		if (thisyear==nyear-1) thisyear=0;
		else thisyear++;
	}

	void firstyear() {
		thisyear=0;
	}

	void get_data_from(double source[][12]) {
		
		int y,m;
		thisyear=0; // guess2008 - ML bugfix
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m]=source[y][m];
			}
		}
	}

	void get_data_from(double* source) {

		int m;
		thisyear=0; // guess2008 - ML bugfix
		for (m=0;m<nyear*12;m++)
			data[m]=source[m];
	}


	// guess2008 - NEW METHODS 

	void reset_clim() {
		for (int ii = 0; ii < 12; ii++) dataclim[ii] = 0.0;
	}


	void make_clim() {
		
		reset_clim(); // Always reset before calculating

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				dataclim[m] += data[y*12+m] / (double)nyear;
			}
		}
	}


	bool extract_data(double source[][12], const int& startyear, const int& endyear) {
		
		// Populate data with data from the middle of source. 
		// Condition: endyear - startyear + 1 == nyear
		// if startyear == 1 and endyear == 30 then this function is identical to get_data_from above.

		if (endyear < startyear) return false;
		if (endyear - startyear + 1 == nyear) {

			int y,m;
			for (y=startyear-1;y<endyear;y++) {
				for (m=0;m<12;m++) {
					data[(y-(startyear-1))*12+m]=source[y][m];
				}
			}

		} else return false;

		return true;
	}


	void adjust_data(double anom[12], bool additive) {
		
		// Adjust the spinup data to the conditions prevailing at a particular time, as given by 
		// the (additive or multiplicative) anomalies in anom 
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (additive)	
					data[y*12+m] += anom[m];
				else
					data[y*12+m] *= anom[m];
			}
		}

	}
	
	
	// Replace interannual data with the period's climatology.
	void use_clim_data() {
	
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m];
			}
		}
	}


	// Alter variability about the mean climatology
	void adjust_data_variability(const double& factor) {
	
		// factor == 0 gives us the climatology (i.e. generalises use_clim_data above)
		// factor == 1 leaves everything unchanged
		// Remember to check the for negative precip or cloudiness values etc. 
		// after calling this method.

		if (factor == 1.0) return;

		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				data[y*12+m] = dataclim[m] + (data[y*12+m] - dataclim[m]) * factor;
			}
		}
	}


	void limit_data(double minval, double maxval) {

		// Limit data to a range
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < minval) data[y*12+m] = minval;
				if (data[y*12+m] > maxval) data[y*12+m] = maxval;
			}
		}

	}
	
	
	void set_min_val(const double& oldval, const double& newval) {

		// Change values < oldval to newval
		int y,m;
		for (y=0;y<nyear;y++) {
			for (m=0;m<12;m++) {
				if (data[y*12+m] < oldval) data[y*12+m] = newval;
			}
		}

	}

	// guess2008 - END OF NEW METHODS


	void detrend_data() {

		int y,m;
		double a,b,anomaly;
		double* annual_mean=new double[nyear];
		double* year_number=new double[nyear];

		if (!annual_mean || !year_number)
			fail("Spinup_driver::detrend_data: out of memory");

		for (y=0;y<nyear;y++) {
			annual_mean[y]=0.0;
			for (m=0;m<12;m++) annual_mean[y]+=data[y*12+m];
			annual_mean[y]/=12.0;
			year_number[y]=y;
		}

		regress(year_number,annual_mean,nyear,a,b);

		for (y=0;y<nyear;y++) {
			anomaly=b*(double)y;
			for (m=0;m<12;m++)
				data[y*12+m]-=anomaly;
		}
		
		// guess2008 - added [] - Clean up
		delete[] annual_mean;
		delete[] year_number;
	}
};

// Constants associated with historical climate data set

// number of years of historical climate
// CRU TS 3.0 has 106 years of data (1901-2006)
const int NYEAR_HIST=106;

// calender year corresponding to first year in CRU climate data set
const int FIRSTHISTYEAR=1901;

// calender year corresponding to first year nitrogen deposition
const int FIRSTHISTYEARNDEP=1850;

// number of years of historical nitrogen deposition 
const int NYEAR_HISTNDEP=16;

// number of years to use for temperature-detrended spinup data set
// (not to be confused with the number of years to spinup model for, which
// is read from the ins file)	
const int NYEAR_SPINUP_DATA=30;

// Stream pointer to binary CRU historical climate data file (read from ins file)
FILE *in_cru;

// Full pathname of ASCII file containing annual CO2 values (read from ins file)
xtring file_co2;


using namespace GuessOutput;

/// The output channel through which all output is sent
OutputChannel* output_channel;

// Output tables
Table out_cmass, out_anpp, out_dens, out_lai, out_cflux, out_cpool, out_firert, out_runoff, out_speciesheights;
Table out_cmass_leaf, out_cmass_sap, out_cmass_heart, out_cmass_root, out_cmass_debt;

Table out_mnpp, out_mlai, out_mgpp, out_mra, out_maet, out_mpet, out_mevap, out_mrunoff, out_mintercep;
Table out_mrh, out_mnee, out_mwcont_upper, out_mwcont_lower;

// Output table for monthly data by stand
Table out_mnpp_extra, out_mlai_extra, out_mgpp_extra, out_mra_extra, out_maet_extra, out_mpet_extra, out_mevap_extra, out_mrunoff_extra, out_mintercep_extra;
Table out_mrh_extra, out_mnee_extra, out_miso_extra, out_mmon_extra;

// bvoc
Table out_aiso, out_miso, out_amon, out_mmon;

FILE* out_dlai;
FILE* out_driver;
Table out_gc;
Table out_landuse;
Table out_par;
Table out_prec;
Table out_temp;
Table out_temp_soil;
// Temporary files until the RCA spinup is forced by precipitation
// Unlike the regular mwcont_* files these are on stand basis so we
// can create the spinup archive based on them.
Table out_mwcont_upper_extra;
Table out_mwcont_lower_extra;

// Spinup type = CRU, RCA (Ben 2007-12-13)
xtring spinup_type;

/// prefix for output file names
xtring out_prefix;

/// Path to RCA climate archive for spinup
xtring file_rcaclim;

Table out_cton_leaf, out_cton_veg, out_nsources, out_npool, out_nuptake, out_vmaxnlim, out_nflux, out_ngases;

// Timers for keeping track of progress through the simulation
Timer tprogress,tmute;
const int MUTESEC=20; // minimum number of sec to wait between progress messages

/// Yearly CO2 data read from file
/**
 * This object is indexed with calendar years, so to get co2 value for
 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
 * more information.
 */
GlobalCO2File global_co2;

// Monthly temperature, precipitation and sunshine data for current grid cell
// and historical period
double hist_mtemp[NYEAR_HIST][12];
double hist_mprec[NYEAR_HIST][12];
double hist_msun[NYEAR_HIST][12];

// Monthly frost days, precipitation days and DTR data for current grid cell
// and historical period
double hist_mfrs[NYEAR_HIST][12];
double hist_mwet[NYEAR_HIST][12];
double hist_mdtr[NYEAR_HIST][12];

/// Monthly data on daily dry NHx deposition (kgN/m2/day)
double NHxDryDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily wet NHx deposition (kgN/m2/day)
double NHxWetDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily dry NOy deposition (kgN/m2/day)
double NOyDryDep[NYEAR_HISTNDEP][12];
/// Monthly data on daily wet NOy deposition (kgN/m2/day)
double NOyWetDep[NYEAR_HISTNDEP][12];

// Spinup data sets for current grid cell
Spinup_data spinup_mtemp(NYEAR_SPINUP_DATA);
Spinup_data spinup_mtemp_soil(NYEAR_SPINUP_DATA);
Spinup_data spinup_mprec(NYEAR_SPINUP_DATA);
Spinup_data spinup_msun(NYEAR_SPINUP_DATA);
Spinup_data spinup_mpar(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_upper(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_lower(NYEAR_SPINUP_DATA);

Spinup_data spinup_mtemp_for(NYEAR_SPINUP_DATA);
Spinup_data spinup_mtemp_soil_for(NYEAR_SPINUP_DATA);
Spinup_data spinup_mpar_for(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_upper_for(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_lower_for(NYEAR_SPINUP_DATA);
Spinup_data spinup_mprec_for(NYEAR_SPINUP_DATA);

Spinup_data spinup_mtemp_opl(NYEAR_SPINUP_DATA);
Spinup_data spinup_mtemp_soil_opl(NYEAR_SPINUP_DATA);
Spinup_data spinup_mpar_opl(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_upper_opl(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwcont_lower_opl(NYEAR_SPINUP_DATA);
Spinup_data spinup_mprec_opl(NYEAR_SPINUP_DATA);

// guess2008
// Spinup data sets for monthly frost days, precipitation days and DTR data for 
// current grid cell
Spinup_data spinup_mfrs(NYEAR_SPINUP_DATA);
Spinup_data spinup_mwet(NYEAR_SPINUP_DATA);
Spinup_data spinup_mdtr(NYEAR_SPINUP_DATA);


// Object for retrieving RCA spinup data from archive file (Ben 2007-12-13)
RCADataArchive ark;

// bvoc
// Daily diurnal temperature range for one year
double ddtr[366];

// Daily N deposition for one year
double dndep[366];

// guess2008 - make file_cru and file_cru_misc global variables
xtring file_cru;
xtring file_cru_misc;

/// Interpolates monthly data to quasi-daily values.
void interp_climate(double mtemp[12], double mprec[12], double msun[12], double mdtr[12],
					double dtemp[366], double dprec[366], double dsun[366], double ddtr[366]) {
	interp_monthly_means(mtemp, dtemp);
	interp_monthly_totals(mprec, dprec);
	interp_monthly_means(msun, dsun);
	interp_monthly_means(mdtr, ddtr);
}

//Landuse:

//#define DYNAMIC_LANDCOVER_INPUT
#if defined DYNAMIC_LANDCOVER_INPUT
//TimeDataD input code may be put here
TimeDataD LUdata(LOCAL_YEARLY);
TimeDataD Peatdata;
#endif
xtring file_lu, file_peat;
const int NYEAR_LU=103;	//only used to get LU data after historical period (after 2003) : only used in AR4-runs, but causes no harm otherwise

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 
// guess2008 - new functions for reading CRU TS 3.0 binary files.
//
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU
// Determine temp, precip, sunshine & soilcode
 
bool searchcru(char* cruark,double dlon,double dlat,int& soilcode,
	double mtemp[NYEAR_HIST][12],double mprec[NYEAR_HIST][12],
	double msun[NYEAR_HIST][12]) {

	// !!!! NEW VERSION OF THIS FUNCTION - guess2008 - NEW VERSION OF THIS FUNCTION !!!!
	// Please note the new function signature. 

	// Archive object. Definition in new header file, cru.h
	Cru_1901_2006Archive ark;

	int target_ilon=(int)(dlon*10.0);
	int target_ilat=(int)(dlat*10.0);

	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006 data; // struct to hold the data

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays. 
		soilcode=(int)data.soilcode[0];


		for (y=0;y<NYEAR_HIST;y++) {
			for (m=0;m<12;m++) {
				mtemp[y][m] = data.mtemp[y*12+m]*0.1; // now degC
				mprec[y][m] = data.mprec[y*12+m]*0.1; // mm (sum over month)
				
				// Limit very low precip amounts because negligible precipitation causes problems 
				// in the prdaily function (infinite loops). 
				if (mprec[y][m] <= 1.0) mprec[y][m] = 0.0;
				
				msun[y][m]  = data.msun[y*12+m]*0.1;   // % sun 

			}
		}


		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}




///////////////////////////////////////////////////////////////////////////////////////
// SEARCHCRU_MISC
// Determine elevation, frs frq, wet frq & DTR

bool searchcru_misc(char* cruark,double dlon,double dlat,int& elevation,
	double mfrs[NYEAR_HIST][12],double mwet[NYEAR_HIST][12],
	double mdtr[NYEAR_HIST][12]) {
	
	// Please note the new function signature. 

	// Archive object
	Cru_1901_2006miscArchive ark; 
	int y,m;

	// Try block to catch any unexpected errors
	try {

		Cru_1901_2006misc data;

		bool success = ark.open(cruark);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) { 
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;


		// The CRU archive index hold lons & lats as whole doubles * 10
		data.lon = dlon * 10.0;
		data.lat = dlat * 10.0;

		// Read the CRU data into the data struct
		success =ark.getindex(data);
		if (!success) {
			ark.close();
			return false;
		}

		// Transfer the data from the data struct to the arrays.
		// Note that the multipliers are NOT the same as in searchcru above!
		elevation=(int)data.elv[0]; // km * 1000

		for (y=0;y<NYEAR_HIST;y++) { 
			for (m=0;m<12;m++) {

				// guess2008 - catch rounding errors 
				mfrs[y][m] = data.mfrs[y*12+m]*0.01; // days
				if (mfrs[y][m] < 0.1) 
					mfrs[y][m] = 0.0; // Catches rounding errors

				mwet[y][m] = data.mwet[y*12+m]*0.01; // days
				if (mwet[y][m] <= 0.1) 
					mwet[y][m] = 0.0; // Catches rounding errors

				mdtr[y][m] = data.mdtr[y*12+m]*0.1;  // degC

				/*
				If vapour pressure is needed:
				mvap[y][m] = data.mvap[y*12+m]*0.01;
				*/
			}
		}

		// Close the archive
		ark.close();

		return true;
	
	}
	catch(...) {
		// Unknown error.
		return false;
	}
}


// Utility function that returns the CRU data from the nearest cell to (lon,lat) within
// a given search radius
bool findnearestCRUdata(double searchradius, char* cruark, double& lon, double& lat, 
                        int& scode, double hist_mtemp1[NYEAR_HIST][12], 
                        double hist_mprec1[NYEAR_HIST][12], 
                        double hist_msun1[NYEAR_HIST][12]) {

	// First try the exact coordinate
	if (searchcru(cruark, lon, lat, scode, hist_mtemp1, hist_mprec1, hist_msun1)) {
		return true;
	}
	
	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(lon*2)/2;
	double center_lat = floor(lat*2)/2;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.5;
	const double EPS = 1e-15;

	for (double y = center_lon-searchradius; y <= center_lon+searchradius+EPS; y += STEP) {
		for (double x = center_lat-searchradius; x <= center_lat+searchradius+EPS; x += STEP) {
			double xdist = x-center_lat;
			double ydist = y-center_lon;
			double dist = sqrt(xdist*xdist + ydist*ydist);
			
			if (dist <= searchradius + EPS) {
				search_points.push_back(make_pair(dist, make_pair(y, x)));
			}
		}
	}

	// Sort by increasing distance
	std::sort(search_points.begin(), search_points.end());

	// Find closest coordinate which can be found in CRU
	for (unsigned int i = 0; i < search_points.size(); i++) {
		point search_point = search_points[i].second;
		double search_lon = search_point.first;
		double search_lat = search_point.second;

		if (searchcru(cruark, search_lon, search_lat, scode, 
		              hist_mtemp1, hist_mprec1, hist_msun1)) {
			lon = search_lon;
			lat = search_lat;
			return true;
		}
	}

	return false;
}


bool find_nearest_cru05(float lon,float lat,float& lon_cru,float& lat_cru,int& soilcode) {
	 double dummy[NYEAR_HIST][12];

	 double dlon = lon;
	 double dlat = lat;

	 bool found = findnearestCRUdata(1, file_cru, dlon, dlat, soilcode, dummy, dummy, dummy);

	 if (found) {
		  lon_cru = dlon;
		  lat_cru = dlat;
	 }

	 return found;
}

/// Checks if the specified date is during an RCA spinup
/** When doing an RCA spinup, we will continue using the forcing
 *  data from the spinup archive up until the first control year
 *  (so not just up until we start getting forcing online from 
 *  RCA). This has to do with the fact that we don't really
 *  trust the forcing we get from RCA in the beginning.
 */
bool during_rca_spinup(const Date& date) {
	return ( spinup_type == "rca" && date.year < first_control_year );
}

xtring get_spinup_type() {
	return spinup_type;
}



bool checkexists(xtring filename) {

	// Ben 2007-11-22

	FILE* in=fopen(filename,"rt");
	if (in) {
		fclose(in);
		dprintf("Warning: %s already exists\n",(char*)filename);
		return true;
	}

	return false;
}


calendartype xtring_to_calendartype(xtring str) {
	 if (str == "STANDARD") {
		  return STANDARD;
	 }
	 else if (str == "LEAPYEARS") {
		  return LEAPYEARS;
	 }
	 else if (str == "FLAT_30") {
		  return FLAT_30;
	 }
	 else {
		  fail("xtring_to_calendartype: Unknown calendartype: %s\n", (char*)str);
		  return STANDARD; // To avoid compiler warning
	 }
}


/// Adds prefix and instance number to a filename
void prepare_filename(xtring& filename, const char* out_prefix, int inst) {

	 if (filename == "") {
		  return;
	 }

	 // For now our tools assume the RCA-GUESS outfiles end in .out,
	 // so we'll make sure they do and add the instance number just before .out

	 long pos = filename.find(".out");

	 if (pos == -1) {
		  fail("Outfiles are required to end in .out in RCA-GUESS");
	 }

	 xtring basename = filename.left(pos);
	 
	 filename.printf("%s_%s_%d.out", out_prefix, (char*)basename, inst);

	 // If one of the output files exist we fail to avoid the risk of
	 // overwriting valuable data from a previous run.
	 // Perhaps a bit ugly to do it here, but all output filenames
	 // goes through here...
	 if (checkexists(filename)) {
		  fail("Output file(s) already exist - please move or delete before running model");
	 }
}

/// Help function to define_output_tables, creates one output table
void create_output_table(Table& table, const char* file, const ColumnDescriptors& columns) {
	 table = output_channel->create_table(TableDescriptor(file, columns));
}

/// Defines all output tables
/** This function specifies all columns in all output tables, their names,
 *  column widths and precision.
 *
 *  For each table a TableDescriptor object is created which is then sent to
 *  the output channel to create the table.
 */
void define_output_tables() {
	// create a vector with the pft names
	std::vector<std::string> pfts;

	pftlist.firstobj();
	while (pftlist.isobj) {
		 Pft& pft=pftlist.getobj();

		 pfts.push_back((char*)pft.name);

		 pftlist.nextobj();
	}

	// create a vector with the landcover column titles
	std::vector<std::string> landcovers;

	if (run_landcover) {
		 const char* landcover_string[]={"Urban_sum", "Crop_sum", "Pasture_sum", "Forest_sum", "Natural_sum", "Peatland_sum"};
		 for (int i=0; i<NLANDCOVERTYPES; i++) {
			  if(run[i]) {
					landcovers.push_back(landcover_string[i]);
			  }
		 }
	}

	// Create the month columns
	ColumnDescriptors month_columns;
	ColumnDescriptors month_columns_wide;
	std::string months_array[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
	std::vector<std::string> months(months_array, months_array+12);
	month_columns = ColumnDescriptors(months, 8, 3);
	month_columns_wide = ColumnDescriptors(months, 10, 3);

	// Create the columns for each output file

	// CMASS
	ColumnDescriptors cmass_columns;
	cmass_columns += ColumnDescriptors(pfts,               8, 3);
	cmass_columns += ColumnDescriptor("Total",             8, 3);
	cmass_columns += ColumnDescriptors(landcovers,        13, 3);

	// ANPP
	ColumnDescriptors anpp_columns = cmass_columns;

	// DENS
	ColumnDescriptors dens_columns;
	dens_columns += ColumnDescriptors(pfts,                8, 4);
	dens_columns += ColumnDescriptor("Total",              8, 4);
	dens_columns += ColumnDescriptors(landcovers,         13, 4);

	// LAI
	ColumnDescriptors lai_columns = dens_columns;

	// CFLUX
	ColumnDescriptors cflux_columns;
	cflux_columns += ColumnDescriptor("Veg",               8, 3);
	cflux_columns += ColumnDescriptor("Soil",              8, 3);
	cflux_columns += ColumnDescriptor("Fire",              8, 3);
	cflux_columns += ColumnDescriptor("Est",               8, 3);
	if (run_landcover) {
		 cflux_columns += ColumnDescriptor("Harvest",      9, 3);
	}
	cflux_columns += ColumnDescriptor("NEE",              10, 5);

	// CPOOL
	ColumnDescriptors cpool_columns;
	cpool_columns += ColumnDescriptor("VegC",              8, 3);

	if (!ifcentury) {
		cpool_columns += ColumnDescriptor("LittC",         8, 3);
		cpool_columns += ColumnDescriptor("SoilfC",        8, 3);
		cpool_columns += ColumnDescriptor("SoilsC",        8, 3);
	}
	else {
		cpool_columns += ColumnDescriptor("LittVC",        8, 3);
		cpool_columns += ColumnDescriptor("LittSC",        8, 3);
		cpool_columns += ColumnDescriptor("CwdC",          8, 3);
		cpool_columns += ColumnDescriptor("SoilC",         8, 3);
	}
	if (run_landcover && ifslowharvestpool) {
		 cpool_columns += ColumnDescriptor("HarvSlowC",   10, 3);
	}
	cpool_columns += ColumnDescriptor("Total",            10, 3);

	// FIRERT
	ColumnDescriptors firert_columns;
	firert_columns += ColumnDescriptor("FireRT",           8, 1);

	// RUNOFF
	ColumnDescriptors runoff_columns;
	runoff_columns += ColumnDescriptor("Surf",             8, 1);
	runoff_columns += ColumnDescriptor("Drain",            8, 1);
	runoff_columns += ColumnDescriptor("Base",             8, 1);
	runoff_columns += ColumnDescriptor("Total",            8, 1);

	// SPECIESHEIGHTS
	ColumnDescriptors speciesheights_columns;
	speciesheights_columns += ColumnDescriptors(pfts,      8, 2);

	// AISO
	ColumnDescriptors aiso_columns;
	aiso_columns += ColumnDescriptors(pfts,               10, 3);
	aiso_columns += ColumnDescriptor("Total",             10, 3);
	aiso_columns += ColumnDescriptors(landcovers,         13, 3);

	// AMON
	ColumnDescriptors amon_columns = aiso_columns;

	// RCA

	ColumnDescriptors landuse_columns;
	landuse_columns += ColumnDescriptor("land", 8, 4);
	landuse_columns += ColumnDescriptor("pnv", 8, 4);
	landuse_columns += ColumnDescriptor("crop", 8, 4);
	landuse_columns += ColumnDescriptor("bare", 8, 4);

	ColumnDescriptors gc_columns;
	gc_columns += ColumnDescriptor("Id", 8, 0);
	gc_columns += month_columns;

	ColumnDescriptors par_columns;
	par_columns += ColumnDescriptor("Id", 8, 0);
	par_columns += ColumnDescriptors(months, 8, 0);

	ColumnDescriptors prec_columns;
	prec_columns += ColumnDescriptor("Id", 8, 0);
	prec_columns += ColumnDescriptors(months, 8, 1);

	ColumnDescriptors temp_columns = prec_columns;

	ColumnDescriptors temp_soil_columns = prec_columns;

	ColumnDescriptors mwcont_upper_extra_columns = gc_columns;
	ColumnDescriptors mwcont_lower_extra_columns = gc_columns;

	ColumnDescriptors month_stand_columns = gc_columns;
	

	//TODO Fix these for landcover

	// CTON
	ColumnDescriptors cton_columns;
	cton_columns += ColumnDescriptors(pfts,                8, 1);
	cton_columns += ColumnDescriptor("Total",              8, 1);
	cton_columns += ColumnDescriptors(landcovers,         12, 1);

	// NSOURCES
	ColumnDescriptors nsources_columns;
	nsources_columns += ColumnDescriptor("dep",            8, 2);
	nsources_columns += ColumnDescriptor("fix",            8, 2);
	nsources_columns += ColumnDescriptor("fert",           8, 2);
	nsources_columns += ColumnDescriptor("input",          8, 2);
	nsources_columns += ColumnDescriptor("min",            7, 2);
	nsources_columns += ColumnDescriptor("imm",            7, 2);
	nsources_columns += ColumnDescriptor("netmin",         7, 2);
	nsources_columns += ColumnDescriptor("Total",          7, 2);

	// NPOOL
	ColumnDescriptors npool_columns;
	npool_columns += ColumnDescriptor("VegN",              9, 5);
	npool_columns += ColumnDescriptor("LittVN",            9, 5);
	npool_columns += ColumnDescriptor("LittSN",            9, 5);
	npool_columns += ColumnDescriptor("CwdN",              9, 5);
	npool_columns += ColumnDescriptor("SoilN",             9, 5);

	if (run_landcover && ifslowharvestpool) {
		npool_columns += ColumnDescriptor("HarvSlowN",     9, 5);
	}

	npool_columns += ColumnDescriptor("Total",            10, 5);

	// NUPTAKE
	ColumnDescriptors nuptake_columns;
	nuptake_columns += ColumnDescriptors(pfts,             7, 2);
	nuptake_columns += ColumnDescriptor("Total",           7, 2);
	nuptake_columns += ColumnDescriptors(landcovers,      10, 2);

	// VMAXNLIM
	ColumnDescriptors vmaxnlim_columns;
	vmaxnlim_columns += ColumnDescriptors(pfts,            8, 2);
	vmaxnlim_columns += ColumnDescriptor("Total",          8, 2);
	vmaxnlim_columns += ColumnDescriptors(landcovers,     13, 2);

	// NFLUX
	ColumnDescriptors nflux_columns;
	nflux_columns += ColumnDescriptor("dep",               8, 2);
	nflux_columns += ColumnDescriptor("fix",               8, 2);
	nflux_columns += ColumnDescriptor("fert",              8, 2);
	nflux_columns += ColumnDescriptor("flux",              8, 2);
	nflux_columns += ColumnDescriptor("leach",             8, 2);
	if (run_landcover) {
		nflux_columns += ColumnDescriptor("harvest",       8, 2);
	}
	nflux_columns += ColumnDescriptor("NEE",               8, 2);

	// NGASES
	ColumnDescriptors ngases_columns;
	ngases_columns += ColumnDescriptor("NH3",              9, 3);
	ngases_columns += ColumnDescriptor("NO",               9, 3);
	ngases_columns += ColumnDescriptor("NO2",              9, 3);
	ngases_columns += ColumnDescriptor("N2O",              9, 3);
	ngases_columns += ColumnDescriptor("NSoil",            9, 3);
	ngases_columns += ColumnDescriptor("Total",            9, 3);

	// *** ANNUAL OUTPUT VARIABLES ***

	create_output_table(out_cmass,          file_cmass,          cmass_columns);
	create_output_table(out_cmass_leaf, file_cmass_leaf, cmass_columns);
	create_output_table(out_cmass_sap, file_cmass_sap, cmass_columns);
	create_output_table(out_cmass_heart, file_cmass_heart, cmass_columns);
	create_output_table(out_cmass_root, file_cmass_root, cmass_columns);
	create_output_table(out_cmass_debt, file_cmass_debt, cmass_columns);
	create_output_table(out_anpp,           file_anpp,           anpp_columns);
	create_output_table(out_dens,           file_dens,           dens_columns);
	create_output_table(out_lai,            file_lai,            lai_columns);
	create_output_table(out_cflux,          file_cflux,          cflux_columns);
	create_output_table(out_cpool,          file_cpool,          cpool_columns);
	create_output_table(out_firert,         file_firert,         firert_columns);
	create_output_table(out_runoff,         file_runoff,         runoff_columns);
	create_output_table(out_speciesheights, file_speciesheights, speciesheights_columns);
	create_output_table(out_aiso,           file_aiso,           aiso_columns);
	create_output_table(out_amon,           file_amon,           amon_columns);

	create_output_table(out_cton_leaf,      file_cton_leaf,      cton_columns);
	create_output_table(out_cton_veg,       file_cton_veg,       cton_columns);
	create_output_table(out_nsources,       file_nsources,       nsources_columns);
	create_output_table(out_npool,          file_npool,          npool_columns);
	create_output_table(out_nuptake,        file_nuptake,        nuptake_columns);
	create_output_table(out_vmaxnlim,       file_vmaxnlim,       vmaxnlim_columns);
	create_output_table(out_nflux,          file_nflux,          nflux_columns);
	create_output_table(out_ngases,         file_ngases,         ngases_columns);

	// *** MONTHLY OUTPUT VARIABLES ***

	create_output_table(out_mnpp,           file_mnpp,           month_columns);
	create_output_table(out_mlai,           file_mlai,           month_columns);
	create_output_table(out_mgpp,           file_mgpp,           month_columns);
	create_output_table(out_mra,            file_mra,            month_columns);
	create_output_table(out_maet,           file_maet,           month_columns);
	create_output_table(out_mpet,           file_mpet,           month_columns);
	create_output_table(out_mevap,          file_mevap,          month_columns);
	create_output_table(out_mrunoff,        file_mrunoff,        month_columns_wide);
	create_output_table(out_mintercep,      file_mintercep,      month_columns);
	create_output_table(out_mrh,            file_mrh,            month_columns);
	create_output_table(out_mnee,           file_mnee,           month_columns);
	create_output_table(out_mwcont_upper,   file_mwcont_upper,   month_columns);
	create_output_table(out_mwcont_lower,   file_mwcont_lower,   month_columns);
	create_output_table(out_miso,           file_miso,           month_columns_wide);
	create_output_table(out_mmon,           file_mmon,           month_columns_wide);

	// *** MONTHLY OUTPUT VARIABLES AT STAND LEVEL ***

	create_output_table(out_mnpp_extra,      file_mnpp_extra,      month_stand_columns);
	create_output_table(out_mlai_extra,      file_mlai_extra,      month_stand_columns);
	create_output_table(out_mgpp_extra,      file_mgpp_extra,      month_stand_columns);
	create_output_table(out_mra_extra,       file_mra_extra,       month_stand_columns);
	create_output_table(out_maet_extra,      file_maet_extra,      month_stand_columns);
	create_output_table(out_mpet_extra,      file_mpet_extra,      month_stand_columns);
	create_output_table(out_mevap_extra,     file_mevap_extra,     month_stand_columns);
	create_output_table(out_mrunoff_extra,   file_mrunoff_extra,   month_stand_columns);
	create_output_table(out_mintercep_extra, file_mintercep_extra, month_stand_columns);
	create_output_table(out_mrh_extra,       file_mrh_extra,       month_stand_columns);
	create_output_table(out_mnee_extra,      file_mnee_extra,      month_stand_columns);
	create_output_table(out_miso_extra,      file_miso_extra,      month_stand_columns);
	create_output_table(out_mmon_extra,      file_mmon_extra,      month_stand_columns);

	// RCA

	create_output_table(out_gc,                 file_gc,                 gc_columns);
	create_output_table(out_landuse,            file_landuse,            landuse_columns);
	create_output_table(out_par,                file_par,                par_columns);
	create_output_table(out_prec,               file_prec,               prec_columns);
	create_output_table(out_temp,               file_temp,               temp_columns);
	create_output_table(out_temp_soil,          file_temp_soil,          temp_soil_columns);
	create_output_table(out_mwcont_upper_extra, file_mwcont_upper_extra, mwcont_upper_extra_columns);
	create_output_table(out_mwcont_lower_extra, file_mwcont_lower_extra, mwcont_lower_extra_columns);
}

///////////////////////////////////////////////////////////////////////////////////////
// INITIO
// Called by the framework at the start of the model run

void initio(const xtring& insfilename, int inst) {

	// DESCRIPTION
	// Initialises input/output (e.g. opening files), sets values for the global
	// simulation parameter variables (currently vegmode, npatch, patcharea,
	// ifbgestab, ifsme, ifstochestab, ifstochmort, iffire,
	// estinterval, npft), initialises pftlist (the one and only list of PFTs and their
	// static parameters for this run of the model). Normally all of the above
	// parameters, and possibly others, are read from the ins file (see above).
	// Function readins should be called to input settings from the ins file.

	// inst = instance number for parallel runs

	///////////////////////////////////////////////////////////////////////////////////
	// GENERIC SECTION - DO NOT MODIFY

	xtring header;
	xtring landcover_type;
	xtring file_landcover;
 

	unixtime(header);
	header=(xtring)"[LPJ-GUESS  "+header+"]\n\n";
	dprintf((char*)header);

	if (!fileexists(insfilename)) {
		fail("Error: could not open %s for input",(const char*)insfilename);
	}

	// Initialise simulation settings and PFT parameters from instruction script
	// Call to readins() returns false if file could not be opened for reading
	// or contained errors (including missing parameters)

	if (!readins(insfilename)) {
		fail("Bad instruction file!");
	}

	// Print the title of this run
	dprintf("\n\n-----------------------------------------------\n%s\n-----------------------------------------------\n",(char*)title);

	///////////////////////////////////////////////////////////////////////////////////
	// USER-SPECIFIC SECTION (Modify as necessary or supply own code)
	//
	// Reads list of grid cells and (optional) description text from grid list file
	// This file should consist of any number of one-line records in the format:
	//   <longitude> <latitude> [<description>]

	/* 
	 *RCA-GUESS doesn't read in a gridlist since the gridcells are supplied 
	 * by RCA.
	 */
	/*
	double dlon,dlat;
	bool eof=false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);
	*/

	file_cru=param["file_cru"].str;
	file_cru_misc=param["file_cru_misc"].str;

	/*
	ngridcell=0;
	while (!eof) {
		
		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
			ngridcell++;
		}
	}


	fclose(in_grid);
	*/

	// Read CO2 data from file
	global_co2.load_file(param["file_co2"].str);

	// Let rca_driver know about the pathname and timeseries for ndep
	// so we can get ndep data also during the coupled period (not just from input module)
	init_RCA_NDep(param["file_ndep"].str, Lamarque::parse_timeseries((char*)ndep_timeseries));

	// Retrieve full pathname to RCA climate spinup archive file and open it
	// Empty string should be specified if no file is available
	// Ben 2007-12-13
	file_rcaclim=param["file_rcaclim"].str;

	// Choose the strategy to use for land cover, and get path name to
	// land cover database.
	landcover_type = param["landcover_type"].str;
	file_landcover = param["file_landcover"].str;

	if (landcover_type.lower() == "pnv") {
		 landcoverstrategy.reset(new PNVStrategy());
	}
	else if (landcover_type.lower() == "ecoclimap") {
		 landcoverstrategy.reset(new EcoclimapStrategy((char*)file_landcover));
	}
	else if (landcover_type.lower() == "yearlytext") {
		 landcoverstrategy.reset(new YearlyStrategy((char*)file_landcover));
	}
	else {
		 fail("Unknown landcover_type: %s", (char*)landcover_type);
	}

	// Spinup type - should be set to "CRU" or "RCA" (Ben 2007-12-13)
	// NB: Spinup will NOT be performed if statefile exists
	spinup_type=param["spinup_type"].str.lower();
	
	if (!(spinup_type=="cru" || spinup_type=="rca"))
		fail("initio: error in %s: spinup_type must be CRU or RCA",(const char*)insfilename);

	// Open RCA spinup data archive
	if (spinup_type=="rca") {

		 if (file_rcaclim=="")
			  fail("initio: file_rcaclim must be specified if spinup type is RCA");

		 if (!ark.open(file_rcaclim))
			  fail("initio: could not open %s for input",(char*)file_rcaclim);
	}

	// Retrieve state file directory and identifier
	state_dir=param["state_dir"].str;
	state_name=param["state_name"].str;
	state_interval=param["state_interval"].str;

	if (state_interval != "spinup" &&
		state_interval != "yearly" &&
		state_interval != "monthly") {
		fail("Please specify 'spinup', 'yearly' or 'monthly' as state_interval");
	}
	
	// Whether vegetation changes feed back to RCA
	ifvegfeedback=param["ifvegfeedback"].num;
	
	// What kind of calendar to use (depends on RCA boundary data)
	calendarmode = xtring_to_calendartype(param["calendarmode"].str);

	// Whether to save GUESS->RCA driver data for control period
	// first and last year to save data for control period

	ifcontrolperiod=param["ifcontrolperiod"].num;
	first_control_year=param["first_control_year"].num;
	last_control_year=param["last_control_year"].num;

	// FastArchive containing GUESS-generated vegetation for control period
	if (!ifvegfeedback)
		file_veg=param["file_veg"].str;

	// Whether to use water content sent in from RCA, or use GUESS' own
	// hydrology driven by precipitation
	ifprescribedwcont=param["ifprescribedwcont"].num;

	// Retrieve output file prefix from ins file and construct output file names
	// e.g. prefix_cmass_0.out

	out_prefix=param["out_prefix"].str;

	// Add the out prefix and instance number to the filenames
	prepare_filename(file_cmass, out_prefix, inst);

	prepare_filename(file_cmass_leaf, out_prefix, inst);
	prepare_filename(file_cmass_sap, out_prefix, inst);
	prepare_filename(file_cmass_heart, out_prefix, inst);
	prepare_filename(file_cmass_root, out_prefix, inst);
	prepare_filename(file_cmass_debt, out_prefix, inst);

	prepare_filename(file_anpp, out_prefix, inst);
	prepare_filename(file_dens, out_prefix, inst);
	prepare_filename(file_lai, out_prefix, inst);
	prepare_filename(file_cflux, out_prefix, inst);
	prepare_filename(file_cpool, out_prefix, inst);
	prepare_filename(file_firert, out_prefix, inst);
	prepare_filename(file_speciesheights, out_prefix, inst);
	prepare_filename(file_runoff, out_prefix, inst);
	prepare_filename(file_mnpp, out_prefix, inst);
	prepare_filename(file_mlai, out_prefix, inst);
	prepare_filename(file_mgpp, out_prefix, inst);
	prepare_filename(file_mra, out_prefix, inst);
	prepare_filename(file_maet, out_prefix, inst);
	prepare_filename(file_mpet, out_prefix, inst);
	prepare_filename(file_mevap, out_prefix, inst);
	prepare_filename(file_mrunoff, out_prefix, inst);
	prepare_filename(file_mintercep, out_prefix, inst);
	prepare_filename(file_mrh, out_prefix, inst);
	prepare_filename(file_mnee, out_prefix, inst);
	prepare_filename(file_mwcont_upper, out_prefix, inst);
	prepare_filename(file_mwcont_lower, out_prefix, inst);
	prepare_filename(file_aiso, out_prefix, inst);
	prepare_filename(file_miso, out_prefix, inst);
	prepare_filename(file_amon, out_prefix, inst);
	prepare_filename(file_mmon, out_prefix, inst);

	prepare_filename(file_mnpp_extra, out_prefix, inst);
	prepare_filename(file_mgpp_extra, out_prefix, inst);
	prepare_filename(file_mra_extra, out_prefix, inst);
	prepare_filename(file_maet_extra, out_prefix, inst);
	prepare_filename(file_mpet_extra, out_prefix, inst);
	prepare_filename(file_mevap_extra, out_prefix, inst);
	prepare_filename(file_mrunoff_extra, out_prefix, inst);
	prepare_filename(file_mintercep_extra, out_prefix, inst);
	prepare_filename(file_mrh_extra, out_prefix, inst);
	prepare_filename(file_mnee_extra, out_prefix, inst);
	prepare_filename(file_miso_extra, out_prefix, inst);
	prepare_filename(file_mmon_extra, out_prefix, inst);

	prepare_filename(file_dlai, out_prefix, inst);
	prepare_filename(file_driver, out_prefix, inst);
	prepare_filename(file_gc, out_prefix, inst);
	prepare_filename(file_landuse, out_prefix, inst);
	prepare_filename(file_par, out_prefix, inst);
	prepare_filename(file_prec, out_prefix, inst);
	prepare_filename(file_temp, out_prefix, inst);
	prepare_filename(file_temp_soil, out_prefix, inst);
	prepare_filename(file_mwcont_upper_extra, out_prefix, inst);
	prepare_filename(file_mwcont_lower_extra, out_prefix, inst);

	prepare_filename(file_cton_leaf, out_prefix, inst);
	prepare_filename(file_cton_veg, out_prefix, inst);
	prepare_filename(file_nsources, out_prefix, inst);
	prepare_filename(file_npool, out_prefix, inst);
	prepare_filename(file_nuptake, out_prefix, inst);
	prepare_filename(file_vmaxnlim, out_prefix, inst);
	prepare_filename(file_nflux, out_prefix, inst);
	prepare_filename(file_ngases, out_prefix, inst);

	/*
	if (run_landcover) {
		all_fracs_const=true;	//If any of the opened files have yearly data, all_fracs_const will be set to false and landcover_dynamics will call get_landcover() each year

		//Retrieve file names for landcover files and open them if static values from ins-file are not used !
		if (!lcfrac_fixed) {	//This version does not support dynamic landcover fraction data

			if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
				file_lu=param["file_lu"].str;
#if defined DYNAMIC_LANDCOVER_INPUT
				if(!LUdata.Open(file_lu))				//Open Bondeau area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_lu);
				else if(LUdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}

			if (run[PEATLAND]) {	//special case for peatland: separate fraction file
				file_peat=param["file_peat"].str;
#if defined DYNAMIC_LANDCOVER_INPUT				
				if(!Peatdata.Open(file_peat))			//Open peatland area fraction file, returned false if problem
					fail("initio: could not open %s for input",(char*)file_peat);
				else if(Peatdata.format==LOCAL_YEARLY)
					all_fracs_const=false;				//Set all_fracs_const to false if yearly data
#endif
			}

		}
	}
	*/
	// We MUST have an output directory
	if (outputdirectory=="") {
		fail("No output directory given in the .ins file!");
	}

	// Create the output channel
	const int COORDINATES_PRECISION = 3; // decimal places for coords in output (3 for RCA coordinates)
	output_channel = new FileOutputChannel((char*)outputdirectory,
														COORDINATES_PRECISION);

	// Define all output tables and their formats
	define_output_tables();

	if (file_dlai!="") {
		file_dlai = outputdirectory + file_dlai;
		out_dlai=fopen(file_dlai,"wb");
		if (!out_dlai) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_dlai);
	}
	else out_dlai=NULL;

	if (file_driver!="" && ifcontrolperiod) {
		file_driver = outputdirectory + file_driver;
		out_driver=fopen(file_driver,"wb");
		if (!out_driver) fail("Could not open %s for output\nClose the file if it is open in another application",(char*)file_driver);
	}
	else out_driver=NULL;

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

///	Loads landcover area fraction data from file(s) for a gridcell.
/** Called from getgridcell() if run_landcover is true. 
  */
bool loadlandcover(Gridcell& gridcell, Coord c)	{
	bool LUerror=false;

	if (!lcfrac_fixed) {
		// Landcover fraction data: read from land use fraction file; dynamic, so data for all years are loaded to LUdata object and 
		// transferred to gridcell.landcoverfrac each year in getlandcover()

		if (run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST]) {
#if defined DYNAMIC_LANDCOVER_INPUT					
			if (!LUdata.Load(c))		//Load area fraction data from Bondeau input file to data object
			{
				dprintf("Problems with landcover fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;		// skip this stand
			}
#endif
		}

		if (run[PEATLAND] && !LUerror) {
#if defined DYNAMIC_LANDCOVER_INPUT
			if(!Peatdata.Load(c))	//special case for peatland: separate fraction file
			{
				dprintf("Problems with natural fractions input file. EXCLUDING STAND at %.3f,%.3f from simulation.\n\n",c.lon,c.lat);
				LUerror=true;	// skip this stand						
			}
#endif
		}
	}

	return LUerror;
}

/// Retrieves nitrogen deposition for a particular gridcell
/** The values are either taken from a binary archive file or when it's not
 *  provided default to pre-industrial level of 2 kgN/ha/year.
 *
 *  The binary archive files have nitrogen deposition in gN/m2 on a monthly timestep
 *  for 16 years with 10 year interval starting from 1850 (Lamarque et. al., 2011).
 *
 *  \param  lon         Longitude
 *  \param  lat         Latitude
 */
/* Not used by RCA, uses ndep functionality from rca_driver.h instead
void getndep(double lon, double lat) {
	
	const double convert = 1e-7;				// converting from gN ha-1 to kgN m-2

	xtring file_ndep = param["file_ndep"].str;

	if (file_ndep == "") {

		// pre-industrial total nitrogen depostion set to 2 kgN/ha/year [kgN m-2]
		double dailyndep = 2000.0 / (4 * 365) * convert;

		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = dailyndep;
				NHxWetDep[y][m] = dailyndep;
				NOyDryDep[y][m] = dailyndep;
				NOyWetDep[y][m] = dailyndep;
			}
		}
	}
	else {
		GlobalNitrogenDepositionArchive ark;
		if (!ark.open(file_ndep)) {
			fail("Could not open %s for input", (char*)file_ndep);
		}

		GlobalNitrogenDeposition rec;
		rec.longitude = lon;
		rec.latitude = lat;

		if (!ark.getindex(rec)) {
			ark.close();
			fail("Grid cell not found in %s", (char*)file_ndep);
		}

		// Found the record, get the values
		for (int y=0; y<NYEAR_HISTNDEP; y++) {
			for (int m=0; m<12; m++) {
				NHxDryDep[y][m] = rec.NHxDry[y*12+m] * convert;
				NHxWetDep[y][m] = rec.NHxWet[y*12+m] * convert;
				NOyDryDep[y][m] = rec.NOyDry[y*12+m] * convert;
				NOyWetDep[y][m] = rec.NOyWet[y*12+m] * convert;
			}
		}
		ark.close();
	}
}
*/

/// Called by the framework at the start of the simulation for a particular grid cell
bool getgridcell(Gridcell& gridcell) {

	// DESCRIPTION
	// Obtains coordinates and soil static parameters for the next grid cell to
	// simulate. The function should return false if no grid cells remain to be simulated,
	// otherwise true. Currently the following member variables of Gridcell should be
	// initialised: longitude, latitude and climate.instype; the following members of
	// member soiltype: awc[0], awc[1], perc_base, perc_exp, thermdiff_0, thermdiff_15,
	// thermdiff_100. The soil parameters can be set indirectly based on an lpj soil
	// code (Sitch et al 2000) by a call to function soilparameters in the driver
	// module (driver.cpp):
	//
	// soilparameters(gridcell.soiltype,soilcode);
	//
	// If the model is to be driven by quasi-daily values of the climate variables
	// derived from monthly means, this function may be the appropriate place to
	// perform the required interpolations. The utility functions interp_monthly_means
	// and interp_monthly_totals in driver.cpp may be called for this purpose.

	// Select coordinates for next grid cell in linked list
	
	int soilcode;
	// guess2008 - elevation
	int elevation;

	// Read spinup data from CRU or RCA archive, depending on spinup_type setting
	// in ins file (Ben 2007-12-13)
	bool LUerror=false;

	if (spinup_type=="cru") {

		// Load environmental data for this grid cell from CRU file
		// NB: soilcode is read in and set in guess.cpp before call to getstand()

		// guess2008 - New searchcru functions takee the CRU filenames as their first 
		// argument, i.e. cru_1901_2002.bin and cru_1901_2002_misc.bin

		// New code:

		double lon = gridcell.lon_cru;
		double lat = gridcell.lat_cru;
		bool gridfound = false;
		gridfound = findnearestCRUdata(searchradius, file_cru, lon, lat, soilcode, 
		                               hist_mtemp, hist_mprec, hist_msun);

		if (gridfound) // Get more historical CRU data for this grid cell
			gridfound = searchcru_misc(file_cru_misc, lon, lat, elevation, 
			                           hist_mfrs, hist_mwet, hist_mdtr);
		/*
		if (run_landcover) {
			Coord& c=gridlist.getobj();
			LUerror=loadlandcover(gridcell, c);
		}
		if (LUerror)
			gridfound=false;
		*/
		if (!gridfound) {
			fail("\nError: could not find gridcell at (%g,%g) in CRU data file\n",
				gridcell.lon_cru, gridcell.lat_cru);
		}

		// Build spinup data sets
		spinup_mtemp.get_data_from(hist_mtemp);
		spinup_mprec.get_data_from(hist_mprec);
		spinup_msun.get_data_from(hist_msun);

		// Detrend spinup temperature data
		spinup_mtemp.detrend_data();

		// guess2008 - new spinup data sets
		spinup_mfrs.get_data_from(hist_mfrs);
		spinup_mwet.get_data_from(hist_mwet);
		spinup_mdtr.get_data_from(hist_mdtr);
		spinup_mdtr.detrend_data();


		dprintf("\nCommencing simulation for gridcell %d at (%g,%g) using CRU data for spinup\n",
				  gridcell.id, gridcell.get_lon(), gridcell.get_lat());
				
		// Get nitrogen deposition data
		// RCA: we're getting ndep with help from rca_driver.h instead
		// getndep(lon, lat);

		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine
		
		gridcell[0].climate.instype=SUNSHINE;
		gridcell[1].climate.instype=SUNSHINE;

	}
	else { // (spinup_type=="rca")

		// IMPORTANT: period covered in RCA climate data archive should
		// correspond exactly to control period specified in ins file

		if (NYEAR_SPINUP_DATA!=last_control_year-first_control_year+1)
			fail("getstand: control period specified in ins file must span %d (NYEAR_SPINUP_DATA) years",NYEAR_SPINUP_DATA);

		// Retrieve data for this grid cell

		RCAData rec; // record for one grid cell in archive	

		rec.lon = gridcell.get_lon();
		if (rec.lon)
			rec.lon=(int)(fabs(rec.lon)*1000.0+0.5)/1000.0*fabs(rec.lon)/rec.lon;

		rec.lat = gridcell.get_lat();
		if (rec.lat)
			rec.lat=(int)(fabs(rec.lat)*1000.0+0.5)/1000.0*fabs(rec.lat)/rec.lat;

		if (!ark.getindex(rec))
			fail("getstand: record for gridcell %d at (%0.3f,%0.3f) not found in %s",
				gridcell.id,rec.lon,rec.lat,(char*)file_rcaclim);
		/*
		if (run_landcover) {
			Coord& c=gridlist.getobj();
			LUerror=loadlandcover(gridcell, c);
		}
		if (LUerror)
			gridfound=false;
		*/
		// Convert PAR from kJ to J

		for (int y=0;y<NYEAR_SPINUP_DATA;y++)
			for (int m=0;m<12;m++) {
				rec.mpar_for[y*12+m]*=1.0e3;
				rec.mpar_open[y*12+m]*=1.0e3;
			}

		// Transfer to global Spinup_data objects for recycling throughout spinup period

		spinup_mtemp_for.get_data_from(rec.mtemp_for);
		spinup_mtemp_soil_for.get_data_from(rec.mtemp_soil_for);
		spinup_mwcont_upper_for.get_data_from(rec.mwcont_upper_for);
		spinup_mwcont_lower_for.get_data_from(rec.mwcont_lower_for);
		spinup_mpar_for.get_data_from(rec.mpar_for);
		spinup_mprec_for.get_data_from(rec.mprec_for);

		spinup_mtemp_opl.get_data_from(rec.mtemp_open);
		spinup_mtemp_soil_opl.get_data_from(rec.mtemp_soil_open);
		spinup_mwcont_upper_opl.get_data_from(rec.mwcont_upper_open);
		spinup_mwcont_lower_opl.get_data_from(rec.mwcont_lower_open);
		spinup_mpar_opl.get_data_from(rec.mpar_open);
		spinup_mprec_opl.get_data_from(rec.mprec_open);

		dprintf("\nCommencing simulation for gridcell %d at (%g,%g) using RCA data for spinup\n",
				  gridcell.id, gridcell.get_lon(), gridcell.get_lat());

		// Get nitrogen deposition data
		// RCA: we're getting ndep with help from rca_driver.h instead
		// getndep(gridcell.lon_cru, gridcell.lat_cru);

		// RCA insolation data is photosynthetically active radiation
		gridcell[0].climate.instype=PAR;
		gridcell[1].climate.instype=PAR;
	}

	// For Windows shell - clear graphical output
	// (ignored on other platforms)
	
	clear_all_graphs();
	
	return true; // simulate this stand
}

/// Synchronise spinup data series with current RCA date (Ben 2007-12-17)
/**
 * The data have to be loaded anew every day during overlap period between RCA and
 * spinup forcing. getclimate() must be called afterwards to transfer data from the
 * global Spinup_data objects to the stand data.
 */
void load_rca_spinup(Gridcell& gridcell) {
	int nyear_wind;
	int m,y;
	double mtemp[12],mtemp_soil[12],mwcont_upper[12],mwcont_lower[12],mpar[12];
	double dtemp[366],dtemp_soil[366],dwcont_upper[366],dwcont_lower[366],dpar[366];
	RCAData *ptr_rec = 0; // record for one grid cell in archive

	if (spinup_type=="rca") {
		
		// Retrieve data for this grid cell from archive

		 static std::map<std::pair<double, double>, RCAData*> memo;

		 if (memo.find(std::make_pair(gridcell.get_lon(), gridcell.get_lat())) != memo.end()) {
			  // get from memo
			  ptr_rec = memo[std::make_pair(gridcell.get_lon(), gridcell.get_lat())];
		 }
		 else {
			  ptr_rec = new RCAData;
			  RCAData& rec(*ptr_rec);
			  
			  rec.lon = gridcell.get_lon();
			  if (rec.lon)
					rec.lon=(int)(fabs(rec.lon)*1000.0+0.5)/1000.0*fabs(rec.lon)/rec.lon;
			  
			  rec.lat = gridcell.get_lat();
			  if (rec.lat)
					rec.lat=(int)(fabs(rec.lat)*1000.0+0.5)/1000.0*fabs(rec.lat)/rec.lat;
			  
			  if (!ark.getindex(rec))
					fail("load_rca_spinup: record for gridcell %d at (%0.3f,%0.3f) not found in %s",
						  gridcell.id,rec.lon,rec.lat,(char*)file_rcaclim);

			  // Convert PAR from kJ to J

			  for (y=0;y<NYEAR_SPINUP_DATA;y++)
					for (m=0;m<12;m++) {
						 rec.mpar_for[y*12+m]*=1.0e3;
						 rec.mpar_open[y*12+m]*=1.0e3;
					}
			  
			  // insert into memo
			  memo[std::make_pair(gridcell.get_lon(), gridcell.get_lat())] = ptr_rec;
		 }
		 RCAData& rec(*ptr_rec);

                // Transfer to global Spinup_data objects for recycling throughout spinup period

		spinup_mtemp_for.get_data_from(rec.mtemp_for);
		spinup_mtemp_soil_for.get_data_from(rec.mtemp_soil_for);
		spinup_mwcont_upper_for.get_data_from(rec.mwcont_upper_for);
		spinup_mwcont_lower_for.get_data_from(rec.mwcont_lower_for);
		spinup_mpar_for.get_data_from(rec.mpar_for);
		spinup_mprec_for.get_data_from(rec.mprec_for);

		spinup_mtemp_opl.get_data_from(rec.mtemp_open);
		spinup_mtemp_soil_opl.get_data_from(rec.mtemp_soil_open);
		spinup_mwcont_upper_opl.get_data_from(rec.mwcont_upper_open);
		spinup_mwcont_lower_opl.get_data_from(rec.mwcont_lower_open);
		spinup_mpar_opl.get_data_from(rec.mpar_open);
		spinup_mprec_opl.get_data_from(rec.mprec_open);

		// No it is not the number of years with moving air but
		// The number of years to go (wind) forward in the spinup data set	
		nyear_wind=(nyear_spinup+date.year-FIRSTHISTYEAR)%NYEAR_SPINUP_DATA;

		spinup_mtemp_for.firstyear();
		spinup_mtemp_soil_for.firstyear();
		spinup_mwcont_upper_for.firstyear();
		spinup_mwcont_lower_for.firstyear();
		spinup_mpar_for.firstyear();
		spinup_mprec_for.firstyear();

		spinup_mtemp_opl.firstyear();
		spinup_mtemp_soil_opl.firstyear();
		spinup_mwcont_upper_opl.firstyear();
		spinup_mwcont_lower_opl.firstyear();
		spinup_mpar_opl.firstyear();
		spinup_mprec_opl.firstyear();

		for (y=0;y<nyear_wind;y++) {
			spinup_mtemp_for.nextyear();
			spinup_mtemp_soil_for.nextyear();
			spinup_mwcont_upper_for.nextyear();
			spinup_mwcont_lower_for.nextyear();
			spinup_mpar_for.nextyear();
			spinup_mprec_for.nextyear();

			spinup_mtemp_opl.nextyear();
			spinup_mtemp_soil_opl.nextyear();
			spinup_mwcont_upper_opl.nextyear();
			spinup_mwcont_lower_opl.nextyear();
			spinup_mpar_opl.nextyear();
			spinup_mprec_opl.nextyear();
		}
	}
}

///	Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
void getlandcover(Gridcell& gridcell) {

	// In RCA we simply get the fractions from the gridcell's LandCoverInfo object.
	for (int i = 0; i < NLANDCOVERTYPES; i++) {
		gridcell.landcoverfrac[i] = 0;
	}

	LandCoverInfo* lci = gridcell.landcoverinfo.get();
	
	int calendar_year = date.year+calendar_year_offset;

	gridcell.landcoverfrac[NATURAL] = (1 - lci->get_fland_crop(calendar_year) - lci->get_fland_bare(calendar_year)) * lci->get_fgrid_land();
	gridcell.landcoverfrac[PASTURE] = lci->get_fland_crop(calendar_year) * lci->get_fgrid_land();

	return;

	int i, year;
	double sum=0.0, sum_tot=0.0, sum_active=0.0;

	if(date.year<nyear_spinup)					//Use values for first historic year during spinup period !
		year=0;
	else if(date.year>=nyear_spinup+NYEAR_LU)	//AR4 adaptation
	{
//		dprintf("setting LU data for scenario period\n");
		year=NYEAR_LU-1;
	}
	else
		year=date.year-nyear_spinup;

	if(lcfrac_fixed)	// If area fractions are set in the ins-file.
	{
		if(date.year==0) // called by landcover_init
		{
			int nactive_landcovertypes=0;

			if(equal_landcover_area)
			{
				for(i=0;i<NLANDCOVERTYPES;i++)
				{
					if(run[i])
						nactive_landcovertypes++;
				}
			}

			for(i=0;i<NLANDCOVERTYPES;i++)
			{
				if(equal_landcover_area)
				{
					sum_active+=gridcell.landcoverfrac[i]=1.0*run[i]/(double)nactive_landcovertypes;	// only set fractions that are active !
					sum_tot=sum_active;
				}
				else
				{
					sum_tot+=gridcell.landcoverfrac[i]=(double)lc_fixed_frac[i]/100.0;					//count sum of all fractions (should be 1.0)

					if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)					//discard unreasonable values
					{
						if(date.year==0)
							dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
						sum_tot-=gridcell.landcoverfrac[i];
						gridcell.landcoverfrac[i]=0.0;
					}

					sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];				//only set fractions that are active !
				}
			}
			
			if(sum_tot<0.99 || sum_tot>1.01)	// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions
				if(date.year==0)
					dprintf("WARNING ! landcover fixed fraction sum is %4.2f, rescaling landcover fractions !\n", sum_tot);

				for(i=0;i<NLANDCOVERTYPES;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}

			//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
			if(sum_active<0.99)	//if landcover types are turned off in the ini-file, always <=1.0 here
			{
				if(date.year==0)
					dprintf("WARNING ! landcover active fraction sum is %4.2f.\n", sum_active);

				if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
				{
					if(date.year==0)
						dprintf("Inactive fractions (%4.2f) transferred to NATURAL fraction.\n", 1.0-sum_active);

					gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference 1.0-(sum of active landcover fractions) are added to the natural fraction
				}
				else
				{
/*					if(date.year==0)
						dprintf("Rescaling landcover fractions !\n");
					for(i=0;i<NLANDCOVERTYPES;i++)
						gridcell.landcoverfrac[i]/=sum_active;			// if NATURAL not simulated, rescale active fractions to 1.0
*/					if(date.year==0)
						dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
				}
																	
			}
		}
	}
	else	//area fractions are read from input file(s);
	{
		if(run[URBAN] || run[CROPLAND] || run[PASTURE] || run[FOREST])
		{	

			for(i=0;i<PEATLAND;i++)		//peatland fraction data is not in this file, otherwise i<NLANDCOVERTYPES.
			{	
#if defined DYNAMIC_LANDCOVER_INPUT
				sum_tot+=gridcell.landcoverfrac[i]=LUdata.Get(year,i);					//count sum of all fractions (should be 1.0)
#endif
				if(gridcell.landcoverfrac[i]<0.0 || gridcell.landcoverfrac[i]>1.0)			//discard unreasonable values
				{		
					if(date.year==0)
						dprintf("WARNING ! landcover fraction size out of limits, set to 0.0\n");
					sum_tot-=gridcell.landcoverfrac[i];
					gridcell.landcoverfrac[i]=0.0;
				}

				sum_active+=gridcell.landcoverfrac[i]=run[i]*gridcell.landcoverfrac[i];
			}

			if(sum_tot!=1.0)		// Check input data, rescale if sum !=1.0
			{
				sum_active=0.0;		//reset sum of active landcover fractions

				if(sum_tot<0.99 || sum_tot>1.01)
				{
					if(date.year==0)
					{
						dprintf("WARNING ! landcover fraction sum is %4.2f for year %d\n", sum_tot, year+FIRSTHISTYEAR);
						dprintf("Rescaling landcover fractions year %d ! (sum is beyond 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);
					}
				}
				else				//added scaling to sum=1.0 (sum often !=1.0)
					dprintf("Rescaling landcover fractions year %d ! (sum is within 0.99-1.01)\n", date.year-nyear_spinup+FIRSTHISTYEAR);

				for(i=0;i<PEATLAND;i++)
					sum_active+=gridcell.landcoverfrac[i]/=sum_tot;
			}
		}
		else
			gridcell.landcoverfrac[NATURAL]=0.0;

		if(run[PEATLAND])
		{
#if defined DYNAMIC_LANDCOVER_INPUT
			sum_active+=gridcell.landcoverfrac[PEATLAND]=Peatdata.Get(year,"PEATLAND");			//peatland fraction data is currently in a separate file !
#endif
		}

		//NB. These calculations are based on the assumption that the NATURAL type area is what is left after the other types are summed. 
		if(sum_active!=1.0)		//if landcover types are turned off in the ini-file, or if more landcover types are added in other input files, can be either less or more than 1.0
		{
			if(date.year==0)
				dprintf("Landcover fraction sum not 1.0 !\n");

			if(run[NATURAL])	//Transfer landcover areas not simulated to NATURAL fraction, if simulated.
			{
				if(date.year==0)
				{
					if(sum_active<1.0)
						dprintf("Inactive fractions (%4.3f) transferred to NATURAL fraction.\n", 1.0-sum_active);
					else
						dprintf("New landcover type fraction (%4.3f) subtracted from NATURAL fraction (%4.3f).\n", sum_active-1.0, gridcell.landcoverfrac[NATURAL]);
				}

				gridcell.landcoverfrac[NATURAL]+=1.0-sum_active;	// difference (can be negative) 1.0-(sum of active landcover fractions) are added to the natural fraction
				
				if(date.year==0)
					dprintf("New NATURAL fraction is %4.3f.\n", gridcell.landcoverfrac[NATURAL]);

				sum_active=1.0;		//sum_active should now be 1.0

				if(gridcell.landcoverfrac[NATURAL]<0.0)	//If new landcover type fraction is bigger than the natural fraction (something wrong in the distribution of input file area fractions)
				{										
					if(date.year==0)
						dprintf("New landcover type fraction is bigger than NATURAL fraction, rescaling landcover fractions !.\n");

					sum_active-=gridcell.landcoverfrac[NATURAL];	//fraction not possible to transfer moved back to sum_active, which will now be >1.0 again
					gridcell.landcoverfrac[NATURAL]=0.0;

					for(i=0;i<NLANDCOVERTYPES;i++)
					{
						gridcell.landcoverfrac[i]/=sum_active;		//fraction rescaled to unity sum
						if(run[i])
							if(date.year==0)
								dprintf("Landcover type %d fraction is %4.3f\n", i, gridcell.landcoverfrac[i]);
					}
				}
			}
			else
			{
//				if(date.year==0)
//					dprintf("Rescaling landcover fractions !\n");
//				for(i=0;i<NLANDCOVERTYPES;i++)
//					gridcell.landcoverfrac[i]/=sum_active;						// if NATURAL not simulated, rescale active fractions to 1.0
				if(date.year==0)
					dprintf("Non-unity fraction sum retained.\n");				// OR let sum remain non-unity
			}
		}
	}
}


/// Called by the framework each simulation day before any process modelling is performed for this day
/** Obtains climate data (including atmospheric CO2 and insolation) for this day. */
bool getclimate(Gridcell& gridcell) {

	// DESCRIPTION
	// The function should returns false if the simulation is complete for this grid cell,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of gridcell must be
	// initialised: co2, temp, prec, insol. If the model is to be driven by quasi-daily
	// values of the climate variables derived from monthly means, this day's values
	// will presumably be extracted from arrays containing the interpolated daily
	// values (see function getgridcell):
	//
	// gridcell.climate.temp=dtemp[date.day];
	// gridcell.climate.prec=dprec[date.day];
	// gridcell.climate.insol=dsun[date.day];
	// 
	// Diurnal temperature range (dtr) added for calculation of leaf temperatures in 
	// BVOC:
	// gridcell.climate.dtr=ddtr[date.day]; 
	//
	// If model is run in diurnal mode, which requires appropriate climate forcing data, 
	// additional members of the climate must be initialised: temps, insols. Both of the
	// variables must be of type std::vector. The length of these vectors should be equal
	// to value of date.subdaily which also needs to be set either in getclimate or 
	// getgridcell functions. date.subdaily is a number of sub-daily period in a single 
	// day. Irrespective of the BVOC settings, climate.dtr variable is not required in 
	// diurnal mode.
	//
	// Note that interpolation to quasi-daily values is done every day in RCA-GUESS
	// because the values are lost when multiple stands are called in succession for
	// the same timestep

	// NB: only called during spinup in RCA-GUESS

	double progress;
	double mtemp[12],mprec[12],msun[12];
	double dtemp[366],dprec[366],dsun[366];
	double mtemp_soil[12],mwcont_upper[12],mwcont_lower[12],mpar[12];
	double dtemp_soil[366],dwcont_upper[366],dwcont_lower[366],dpar[366];

	// guess2008 - changed name from mwet to mwet_all
	// double mwet_all[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // number of rain days per month

	if (spinup_type=="cru") {

		// NORMAL (CRU+HISTORICAL) SPINUP

		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down

		double mndrydep[12], mnwetdep[12];

		Lamarque::NDepData& ndep = getNDepData(gridcell.lon_cru, gridcell.lat_cru);
		ndep.get_one_calendar_year(date.year+calendar_year_offset, mndrydep, mnwetdep);
		
		if (date.year < nyear_spinup) {

			// During spinup period
			// (note conversion of mprec to mean daily precipitation)

			int m;
			double mfrs[12],mwet[12],mdtr[12];

			for (m=0;m<12;m++) {
				mtemp[m] = spinup_mtemp[m];
				mprec[m] = spinup_mprec[m];
				msun[m]	 = spinup_msun[m];

				// guess2008
				mfrs[m] = spinup_mfrs[m];
				mwet[m] = spinup_mwet[m];
				mdtr[m] = spinup_mdtr[m];
			}

			// Interpolate monthly spinup data to quasi-daily values
			interp_climate(mtemp,mprec,msun,mdtr,dtemp,dprec,dsun,ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if rainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(mprec, dprec, mwet, gridcell.seed);
			}
			
			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);

			if (date.islastmonth && date.islastday) {
				spinup_mtemp.nextyear();
				spinup_mprec.nextyear();
				spinup_msun.nextyear();
				
				// guess2008
				spinup_mfrs.nextyear();
				spinup_mwet.nextyear();
				spinup_mdtr.nextyear();
			}

		}
		else if (date.year < nyear_spinup + NYEAR_HIST) {

			// Historical period

			// Interpolate this year's monthly data to quasi-daily values
			interp_climate(hist_mtemp[date.year-nyear_spinup], 
								hist_mprec[date.year-nyear_spinup],
								hist_msun[date.year-nyear_spinup],
								hist_mdtr[date.year-nyear_spinup],
								dtemp, dprec, dsun, ddtr);

			// guess2008 - only recalculate precipitation values using weather generator
			// if ifrainonwetdaysonly is true. Otherwise we assume that it rains a little every day.
			if (ifrainonwetdaysonly) { 
				// (from Dieter Gerten 021121)
				prdaily(hist_mprec[date.year-nyear_spinup], dprec, hist_mwet[date.year-nyear_spinup], gridcell.seed);
			}

			// Distribute N deposition
			distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
		}
		else {
			// Return false if last year was the last for the simulation
			return false;
		}

		gridcell[0].climate.temp=dtemp[date.day];
		gridcell[0].climate.prec=dprec[date.day];
		gridcell[0].climate.insol=dsun[date.day];

		gridcell[1].climate.temp=dtemp[date.day];
		gridcell[1].climate.prec=dprec[date.day];
		gridcell[1].climate.insol=dsun[date.day];
	}
	else  { // (spinup_type=="rca")

		// SPINUP USING ARCHIVED RCA OUTPUT FROM PREVIOUS RUN
		// Use spinup data up to start of control period

		// Extract N deposition to use for this year,
		// monthly means to be distributed into daily values further down
		double mndrydep[12], mnwetdep[12];

		Lamarque::NDepData& ndep = getNDepData(gridcell.lon_cru, gridcell.lat_cru);
		ndep.get_one_calendar_year(date.year+calendar_year_offset, mndrydep, mnwetdep);

		double mtemp_for[12], mtemp_opl[12];
		double dtemp_for[366],dtemp_opl[366];
		double mtemp_soil_for[12], mtemp_soil_opl[12], mwcont_upper_for[12], mwcont_upper_opl[12], mwcont_lower_for[12], mwcont_lower_opl[12], mpar_for[12], mpar_opl[12], mprec_for[12], mprec_opl[12];
		double dtemp_soil_for[366], dtemp_soil_opl[366], dwcont_upper_for[366], dwcont_upper_opl[366], dwcont_lower_for[366], dwcont_lower_opl[366], dpar_for[366], dpar_opl[366], dprec_for[366], dprec_opl[366];

		for (int m=0;m<12;m++) {
			mtemp_for[m]=spinup_mtemp_for[m];
			mtemp_soil_for[m]=spinup_mtemp_soil_for[m];
			mwcont_upper_for[m]=spinup_mwcont_upper_for[m];
			mwcont_lower_for[m]=spinup_mwcont_lower_for[m];
			mpar_for[m]=spinup_mpar_for[m];
			mprec_for[m]=spinup_mprec_for[m];

			mtemp_opl[m]=spinup_mtemp_opl[m];
			mtemp_soil_opl[m]=spinup_mtemp_soil_opl[m];
			mwcont_upper_opl[m]=spinup_mwcont_upper_opl[m];
			mwcont_lower_opl[m]=spinup_mwcont_lower_opl[m];
			mpar_opl[m]=spinup_mpar_opl[m];
			mprec_opl[m]=spinup_mprec_opl[m];
		}

		// Interpolate to quasi-daily values	
		interp_monthly_means(mtemp_for,dtemp_for);
		interp_monthly_means(mtemp_soil_for,dtemp_soil_for);
		interp_monthly_means(mwcont_upper_for,dwcont_upper_for);
		interp_monthly_means(mwcont_lower_for,dwcont_lower_for);
		interp_monthly_means(mpar_for,dpar_for);
		interp_monthly_totals(mprec_for,dprec_for);

		// Distribute N deposition
		// TODO: should perhaps have a separate dndep for forest and open land,
		//       however for now we don't even have number of wet days in the
		//       RCA spinup, so this is probably good enough
		distribute_ndep(mndrydep, mnwetdep, dprec_for, dndep);

		interp_monthly_means(mtemp_opl,dtemp_opl);
		interp_monthly_means(mtemp_soil_opl,dtemp_soil_opl);
		interp_monthly_means(mwcont_upper_opl,dwcont_upper_opl);
		interp_monthly_means(mwcont_lower_opl,dwcont_lower_opl);
		interp_monthly_means(mpar_opl,dpar_opl);
		interp_monthly_totals(mprec_opl,dprec_opl);

		// We currently don't have dtr in the RCA spinup archive
		for (int d = 0; d < 366; d++) {
			ddtr[d] = 0;
		}

		// On last day of year, advance to next year in spinup time series

		if (date.islastmonth && date.islastday) {
			spinup_mtemp_for.nextyear();
			spinup_mtemp_soil_for.nextyear();
			spinup_mwcont_upper_for.nextyear();
			spinup_mwcont_lower_for.nextyear();
			spinup_mpar_for.nextyear();
			spinup_mprec_for.nextyear();

			spinup_mtemp_opl.nextyear();
			spinup_mtemp_soil_opl.nextyear();
			spinup_mwcont_upper_opl.nextyear();
			spinup_mwcont_lower_opl.nextyear();
			spinup_mpar_opl.nextyear();
			spinup_mprec_opl.nextyear();
		}

		// Every day
		// Save today's values in rca_spinup_data member of stand

		gridcell[0].climate.temp = dtemp_opl[date.day];
		gridcell[0].wcont_prescribed[0] = dwcont_upper_opl[date.day];
		gridcell[0].wcont_prescribed[1] = dwcont_lower_opl[date.day];
		gridcell[0].climate.insol = dpar_opl[date.day];
		gridcell[0].climate.prec = dprec_opl[date.day];

		gridcell[1].climate.temp = dtemp_for[date.day];
		gridcell[1].wcont_prescribed[0] = dwcont_upper_for[date.day];
		gridcell[1].wcont_prescribed[1] = dwcont_lower_for[date.day];
		gridcell[1].climate.insol = dpar_for[date.day];
		gridcell[1].climate.prec = dprec_for[date.day];
	}

	// Send environmental values for today to framework

	gridcell[0].climate.co2 = global_co2[date.year+calendar_year_offset];
	gridcell[1].climate.co2 = global_co2[date.year+calendar_year_offset];

	// Nitrogen deposition
	gridcell[0].climate.dndep = dndep[date.day];
	gridcell[1].climate.dndep = dndep[date.day];

	// Nitrogen fertilization
	gridcell[0].climate.dnfert = 0.0;
	gridcell[1].climate.dnfert = 0.0;
	
	// bvoc
	if(ifbvoc){
	  gridcell[0].climate.dtr = ddtr[date.day];
	  gridcell[1].climate.dtr = ddtr[date.day];
	}

	// First day of year only ...

	/* Don't print progress in RCA-GUESS
	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {
			progress=(double)(gridlist.getobj().id*(nyear_spinup+NYEAR_HIST)
				+date.year)/(double)(ngridcell*(nyear_spinup+NYEAR_HIST));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}
	*/
	return true;
}

/// Prints out climate forcing data to be used in the RCA spinup (and also gc)
/** These files all print out one line per stand, with values on stand basis,
 *  so that the RCA spinup archive can have different forcing for different
 *  tiles.
 */
void outclimate(Gridcell& gridcell) {

	// We always print calendar years (!= date.year during spinup)
	int cyear = date.year+calendar_year_offset;

	for (int i = 0; i < 2; i++) {

		OutputRows out(output_channel, gridcell.get_lon(), gridcell.get_lat(), cyear);

		Stand& stand = gridcell[i];

		std::vector<double> mwcont_upper(12,0);
		std::vector<double> mwcont_lower(12,0);

		for (int p = 0; p < stand.npatch(); p++) {
			Patch& patch = stand[p];
			for (int m = 0; m < 12; m++) {
				mwcont_upper[m] += patch.soil.mwcont[m][0]/stand.npatch();
				mwcont_lower[m] += patch.soil.mwcont[m][1]/stand.npatch();
			}
		}

		// The tool which uses the climate out files to create a binary archive
		// assumes that lines with an odd id are open land. Stands don't have
		// their own unique ids anymore since the gridcell class was introduced,
		// but we'll still print out odd and even ids here just to keep that tool
		// working.
		if ((i == 0 && stand.landcover != PASTURE) ||
		    (i == 1 && stand.landcover != NATURAL)) {
			fail("Failed assumption in outclimate: assumed two stands, PASTURE followed by NATURAL");
		}

		out.add_value(out_temp, i+1);
		out.add_value(out_temp_soil, i+1);
		out.add_value(out_par, i+1);
		out.add_value(out_gc, i+1);
		out.add_value(out_prec, i+1);
		out.add_value(out_mwcont_upper_extra, i+1);
		out.add_value(out_mwcont_lower_extra, i+1);

		for (int m = 0; m < 12; m++) {
			 
			 out.add_value(out_temp, stand.climate.mtemp_out[m]);
			 out.add_value(out_par, stand.climate.mpar_out[m]*1.0e-3); // J -> KJ
			 out.add_value(out_gc, stand.gc_out[m]);
			 out.add_value(out_prec, stand.climate.mprec_out[m]);
			 out.add_value(out_mwcont_upper_extra, mwcont_upper[m]);
			 out.add_value(out_mwcont_lower_extra, mwcont_lower[m]);

			 double temp_soil_sum=0.0;

			 for (int p = 0; p < stand.npatch(); p++) {
				  Patch& patch = stand[p];
				  
				  temp_soil_sum += patch.soil.mtemp_out[m]/(double)stand.npatch();
			 }
		
			 out.add_value(out_temp_soil, temp_soil_sum);
		}
	}
}


/// Help function to prepare C:N values for output
/** Avoids division by zero and limits the results to a maximum
 *  value to avoid inf or values large enough to ruin the alignment
 *  in the output.
 *
 *  If both cmass and nmass is 0, the function returns 0.
 */ 
double limited_cton(double cmass, double nmass) {
	const double MAX_CTON = 1000;

	if (nmass > 0.0) {
		return min(MAX_CTON, cmass / nmass);
	}
	else if (cmass > 0.0) {
		return MAX_CTON;
	}
	else {
		return 0.0;
	}
}

/// Prints out monthly ecosystem out at stand level for the analysis at stand level -mcw151016
/** The previous monthly output has been aggregated with landuse change information
* P.S.: PFT information in annual output are already at stand level without aggregating landuse
*/
void outmonth_stand(Gridcell& gridcell) {

	// We always print calendar years (!= date.year during spinup)
	int cyear = date.year + calendar_year_offset;

	for (int i = 0; i < 2; i++) {

		OutputRows out(output_channel, gridcell.get_lon(), gridcell.get_lat(), cyear);

		Stand& stand = gridcell[i];

		// Assuming that lines with an odd id are open land.
		if ((i == 0 && stand.landcover != PASTURE) ||
			(i == 1 && stand.landcover != NATURAL)) {
			fail("Failed assumption in outclimate: assumed two stands, PASTURE followed by NATURAL");
		}

		/* Not necessary to story monthly LAI here
		Since these can be found from RCA (daily data there)  -mcw */
		std::vector<double> mnpp(12, 0);
		std::vector<double> mgpp(12, 0);
		//std::vector<double> mlai(12, 0);
		std::vector<double> maet(12, 0);
		std::vector<double> mpet(12, 0);
		std::vector<double> mevap(12, 0);
		std::vector<double> mintercep(12, 0);
		std::vector<double> mrunoff(12, 0);
		std::vector<double> mrh(12, 0);
		std::vector<double> mra(12, 0);
		std::vector<double> mnee(12, 0);
		std::vector<double> miso(12, 0);
		std::vector<double> mmon(12, 0);

		// Average across patches
		for (int p = 0; p < stand.npatch(); p++) {
			Patch& patch = stand[p];
			for (int m = 0; m < 12; m++) {

				mgpp[m] += patch.fluxes.get_monthly_flux(Fluxes::GPP, m) / (double)stand.npatch();
				mra[m] += patch.fluxes.get_monthly_flux(Fluxes::RA, m) / (double)stand.npatch();
				maet[m] += patch.maet[m] / (double)stand.npatch();
				mpet[m] += patch.mpet[m] / (double)stand.npatch();
				mevap[m] += patch.mevap[m] / (double)stand.npatch();
				mrunoff[m] += patch.mrunoff[m] / (double)stand.npatch();
				mintercep[m] += patch.mintercep[m] / (double)stand.npatch();
				mrh[m] += patch.fluxes.get_monthly_flux(Fluxes::SOILC, m) / (double)stand.npatch();

				miso[m] += patch.fluxes.get_monthly_flux(Fluxes::ISO, m) / (double)stand.npatch();
				mmon[m] += patch.fluxes.get_monthly_flux(Fluxes::MON, m) / (double)stand.npatch();
			}
		}
		// In contrast to annual NEE, monthly NEE does not include fire
		// or establishment fluxes
		for (int m = 0; m<12; m++) {
			mnpp[m] = mgpp[m] - mra[m];
			mnee[m] = mnpp[m] - mrh[m];
		}

		// Write IDs
		out.add_value(out_mnpp_extra, i + 1);
		//out.add_value(out_mlai_extra, i + 1);
		out.add_value(out_mgpp_extra, i + 1);
		out.add_value(out_mra_extra, i + 1);
		out.add_value(out_maet_extra, i + 1);
		out.add_value(out_mpet_extra, i + 1);
		out.add_value(out_mevap_extra, i + 1);
		out.add_value(out_mrunoff_extra, i + 1);
		out.add_value(out_mintercep_extra, i + 1);
		out.add_value(out_mrh_extra, i + 1);
		out.add_value(out_mnee_extra, i + 1);
		out.add_value(out_miso_extra, i + 1);
		out.add_value(out_mmon_extra, i + 1);

		for (int m = 0; m < 12; m++) {
			// Print monthly output variables

			out.add_value(out_mnpp_extra, mnpp[m]);
			//out.add_value(out_mlai_extra, mlai[m]);
			out.add_value(out_mgpp_extra, mgpp[m]);
			out.add_value(out_mra_extra, mra[m]);
			out.add_value(out_maet_extra, maet[m]);
			out.add_value(out_mpet_extra, mpet[m]);
			out.add_value(out_mevap_extra, mevap[m]);
			out.add_value(out_mrunoff_extra, mrunoff[m]);
			out.add_value(out_mintercep_extra, mintercep[m]);
			out.add_value(out_mrh_extra, mrh[m]);
			out.add_value(out_mnee_extra, mnee[m]);
			out.add_value(out_miso_extra, miso[m]);
			out.add_value(out_mmon_extra, mmon[m]);
		}
	}
}

/// Called by the framework at the end of the last day of each simulation year
void outannual(Gridcell& gridcell) {

	// DESCRIPTION
	// Output of simulation results at the end of each year, or for specific years in
	// the simulation of each stand or grid cell. This function does not have to
	// provide any information to the framework.

	int c, m, nclass;
	double flux_veg, flux_soil, flux_fire, flux_est, flux_charvest;
	double c_litter, c_fast, c_slow, c_harv_slow; 

	double surfsoillitterc,surfsoillittern,cwdc,cwdn,centuryc,centuryn,n_litter,n_harv_slow,availn;
	double flux_nh3,flux_no,flux_no2,flux_n2o,flux_nsoil,flux_ntot,flux_nharvest;

	// Nitrogen output is in kgN/ha instead of kgC/m2 as for carbon 
	double m2toha = 10000.0;

	// hold the monthly average across patches
	double mnpp[12];
	double mgpp[12];
	double mlai[12];
	double maet[12];
	double mpet[12];
	double mevap[12];
	double mintercep[12];
	double mrunoff[12];
	double mrh[12];
	double mra[12];
	double mnee[12];
	double mwcont_upper[12];
	double mwcont_lower[12];
	// bvoc
	double miso[12];
	double mmon[12];

	float last=-1;

	if (vegmode == COHORT)
		nclass = min(date.year / estinterval + 1, OUTPUT_MAXAGECLASS);

	outclimate(gridcell);
	outmonth_stand(gridcell);

	
	// guess2008 - yearly output after spinup
		
	// If only yearly output between, say 1961 and 1990 is requred, use: 
	//	if (date.year>=nyear_spinup+60 && date.year<nyear_spinup+90) {

	if (true) { // Output results from spinup as well in RCA-GUESS

		double lon = gridcell.get_lon();
		double lat = gridcell.get_lat();

		// We always print calendar years (!= date.year during spinup)
		int cyear = date.year+calendar_year_offset;

		// The OutputRows object manages the next row of output for each
		// output table
		OutputRows out(output_channel, lon, lat, cyear);

		// Write leader for new record in binary daily LAI file
		// Record format:
		// <lon_rca> <lat_rca> <gridcell.id> <date.year> <ndayyear> <npft>
		// <pft.id> 0 <lai> <day> <lai> <day> <lai> ... <day> <lai> <ndayyear> <lai>
		// <pft.id> 0 <lai> <day> <lai> <day> <lai> ... <day> <lai> <ndayyear> <lai>
		//     .    .   .     .     .     .     .         .     .        .       .
		//     .    .   .     .     .     .     .         .     .        .       .
		//     .    .   .     .     .     .     .         .     .        .       .
		// <pft.id> 0 <lai> <day> <lai> <day> <lai> ... <day> <lai> <ndayyear> <lai>
		
		if (out_dlai) fwrite(&lon,sizeof(lon),1,out_dlai);
		if (out_dlai) fwrite(&lat,sizeof(lat),1,out_dlai);
		if (out_dlai) fwrite(&gridcell.id,sizeof(gridcell.id),1,out_dlai);
		if (out_dlai) fwrite(&date.year,sizeof(date.year),1,out_dlai);
		if (out_dlai) fwrite(&date.ndayyear,sizeof(date.ndayyear),1,out_dlai);
		if (out_dlai) fwrite(&npft,sizeof(npft),1,out_dlai);


		// guess2008 - reset monthly average across patches each year
		for (m=0;m<12;m++)
			mnpp[m]=mlai[m]=mgpp[m]=mra[m]=maet[m]=mpet[m]=mevap[m]=mintercep[m]=mrunoff[m]=mrh[m]=mnee[m]=mwcont_upper[m]=mwcont_lower[m]=miso[m]=mmon[m]=0.0;



		double landcover_cmass[NLANDCOVERTYPES]={0.0};
		double landcover_cmass_leaf[NLANDCOVERTYPES]={0.0};
		double landcover_nmass_leaf[NLANDCOVERTYPES]={0.0};

		double landcover_cmass_sap[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_heart[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_root[NLANDCOVERTYPES] = { 0.0 };
		double landcover_cmass_debt[NLANDCOVERTYPES] = { 0.0 };

		double landcover_cmass_veg[NLANDCOVERTYPES]={0.0};
		double landcover_nmass_veg[NLANDCOVERTYPES]={0.0};
		double landcover_anpp[NLANDCOVERTYPES]={0.0};
		double landcover_lai[NLANDCOVERTYPES]={0.0};
		double landcover_densindiv_total[NLANDCOVERTYPES]={0.0};
		double landcover_aiso[NLANDCOVERTYPES]={0.0};
		double landcover_amon[NLANDCOVERTYPES]={0.0};
		double landcover_nuptake[NLANDCOVERTYPES]={0.0};
		double landcover_vmaxnlim[NLANDCOVERTYPES]={0.0};

		double gcpft_cmass=0.0;
		double gcpft_nmass=0.0;
		double gcpft_cmass_leaf=0.0;
		double gcpft_nmass_leaf=0.0;

		double gcpft_cmass_sap = 0.0;
		double gcpft_cmass_heart = 0.0;
		double gcpft_cmass_root = 0.0;
		double gcpft_cmass_debt = 0.0;

		double gcpft_cmass_veg=0.0;
		double gcpft_nmass_veg=0.0;
		double gcpft_anpp=0.0;
		double gcpft_lai=0.0;
		double gcpft_densindiv_total=0.0;
		double gcpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double gcpft_aiso=0.0;
		double gcpft_amon=0.0;
		double gcpft_nuptake=0.0;
		double gcpft_vmaxnlim=0.0;

		double cmass_gridcell=0.0;
		double nmass_gridcell= 0.0;
		double cmass_leaf_gridcell=0.0;
		double nmass_leaf_gridcell=0.0;

		double cmass_sap_gridcell = 0.0;
		double cmass_heart_gridcell = 0.0;
		double cmass_root_gridcell = 0.0;
		double cmass_debt_gridcell = 0.0;

		double cmass_veg_gridcell=0.0;
		double nmass_veg_gridcell=0.0;
		double anpp_gridcell=0.0;
		double lai_gridcell=0.0;
		double surfrunoff_gridcell=0.0;
		double drainrunoff_gridcell=0.0;
		double baserunoff_gridcell=0.0;
		double runoff_gridcell=0.0;
		double dens_gridcell=0.0;
		double firert_gridcell=0.0;
		double fireprob_gridcell = 0.0;
		double aiso_gridcell=0.0;
		double amon_gridcell=0.0;
		double nuptake_gridcell=0.0;
		double vmaxnlim_gridcell=0.0;

		double andep_gridcell=0.0;
		double anfert_gridcell=0.0;
		double anmin_gridcell=0.0;
		double animm_gridcell=0.0;
		double anfix_gridcell=0.0;
		double n_min_leach_gridcell=0.0;
		double n_org_leach_gridcell=0.0;

		double standpft_cmass=0.0;
		double standpft_nmass=0.0;
		double standpft_cmass_leaf=0.0;
		double standpft_nmass_leaf=0.0;

		double standpft_cmass_sap = 0.0;
		double standpft_cmass_heart = 0.0;
		double standpft_cmass_root = 0.0;
		double standpft_cmass_debt = 0.0;

		double standpft_cmass_veg=0.0;
		double standpft_nmass_veg=0.0;
		double standpft_anpp=0.0;
		double standpft_lai=0.0;
		double standpft_densindiv_total=0.0;
		double standpft_densindiv_ageclass[OUTPUT_MAXAGECLASS]={0.0};
		double standpft_aiso=0.0;
		double standpft_amon=0.0;
		double standpft_nuptake=0.0;
		double standpft_vmaxnlim=0.0;


		// *** Loop through PFTs ***

		pftlist.firstobj();
		while (pftlist.isobj) {
			
			Pft& pft=pftlist.getobj();
			Gridcellpft& gridcellpft=gridcell.pft[pft.id];

			// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs		
			gcpft_cmass=0.0;
			gcpft_nmass=0.0;
			gcpft_cmass_leaf=0.0;
			gcpft_nmass_leaf=0.0;

			gcpft_cmass_sap = 0.0;
			gcpft_cmass_heart = 0.0;
			gcpft_cmass_root = 0.0;
			gcpft_cmass_debt = 0.0;

			gcpft_cmass_veg=0.0;
			gcpft_nmass_veg=0.0;
			gcpft_anpp=0.0;
			gcpft_lai=0.0;
			gcpft_densindiv_total=0.0;		
			gcpft_aiso=0.0;
			gcpft_amon=0.0;
			gcpft_nuptake=0.0;
			gcpft_vmaxnlim=0.0;

			double heightindiv_total = 0.0;

			gridcell.firstobj();

			// Loop through Stands
			while (gridcell.isobj) {
				Stand& stand=gridcell.getobj();

				Standpft& standpft=stand.pft[pft.id];
				// Sum C biomass, NPP, LAI and BVOC fluxes across patches and PFTs
				standpft_cmass=0.0;
				standpft_nmass=0.0;
				standpft_cmass_leaf=0.0;
				standpft_nmass_leaf=0.0;

				standpft_cmass_sap = 0.0;
				standpft_cmass_heart = 0.0;
				standpft_cmass_root = 0.0;
				standpft_cmass_debt = 0.0;

				standpft_cmass_veg=0.0;
				standpft_nmass_veg=0.0;
				standpft_anpp=0.0;
				standpft_lai=0.0;
				standpft_densindiv_total = 0.0;
				standpft_aiso=0.0;
				standpft_amon=0.0;
				standpft_nuptake=0.0;
				standpft_vmaxnlim=0.0;

				// Initialise age structure array

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++){
						standpft_densindiv_ageclass[c] = 0.0;
					}

				stand.firstobj();

				// Loop through Patches
				while (stand.isobj) {
					Patch& patch = stand.getobj();

					standpft_anpp += patch.fluxes.get_annual_flux(Fluxes::NPP, pft.id);
					standpft_aiso += patch.fluxes.get_annual_flux(Fluxes::ISO, pft.id);
					standpft_amon += patch.fluxes.get_annual_flux(Fluxes::MON, pft.id);

					Vegetation& vegetation = patch.vegetation;

					vegetation.firstobj();
					while (vegetation.isobj) {
						Individual& indiv=vegetation.getobj();
							
						// guess2008 - alive check added
						if (indiv.id!=-1 && indiv.alive) { 
							
							if (indiv.pft.id==pft.id) {
								standpft_cmass += indiv.cmass_leaf + indiv.cmass_root + 
								                  indiv.cmass_wood();
								standpft_nmass += indiv.nmass_leaf + indiv.nmass_root + 
								                  indiv.nmass_wood() + indiv.nstore();
								standpft_cmass_leaf += indiv.cmass_leaf;
								standpft_nmass_leaf += indiv.cmass_leaf / indiv.cton_leaf_aavr;

								standpft_cmass_sap += indiv.cmass_sap;
								standpft_cmass_heart += indiv.cmass_heart;
								standpft_cmass_root += indiv.cmass_root;
								standpft_cmass_debt += indiv.cmass_debt;

								standpft_cmass_veg += indiv.cmass_veg;
								standpft_nmass_veg += indiv.nmass_veg;
								standpft_lai += indiv.lai;
								standpft_vmaxnlim += indiv.avmaxnlim * indiv.cmass_leaf;
								standpft_nuptake += indiv.anuptake;

								if (vegmode==COHORT || vegmode==INDIVIDUAL) {
									
									// Age structure
									
									c=(int)(indiv.age/estinterval); // guess2008
									if (c<OUTPUT_MAXAGECLASS)
										standpft_densindiv_ageclass[c]+=indiv.densindiv;

									// guess2008 - only count trees with a trunk above a certain diameter  
									if (pft.lifeform==TREE && indiv.age>0) {
										double diam=pow(indiv.height/indiv.pft.k_allom2,1.0/indiv.pft.k_allom3);
										if (diam>0.03) {
											standpft_densindiv_total+=indiv.densindiv; // indiv/m2

											heightindiv_total+=indiv.height * indiv.densindiv;
										}
									}
								}
							
							}

						} // alive?
						vegetation.nextobj();
					}

					stand.nextobj();
				} // end of patch loop

				standpft_cmass/=(double)stand.npatch();
				standpft_nmass/=(double)stand.npatch();
				standpft_cmass_leaf/=(double)stand.npatch();
				standpft_nmass_leaf/=(double)stand.npatch();

				standpft_cmass_sap /= (double)stand.npatch();
				standpft_cmass_heart /= (double)stand.npatch();
				standpft_cmass_root /= (double)stand.npatch();
				standpft_cmass_debt /= (double)stand.npatch();

				standpft_cmass_veg/=(double)stand.npatch();
				standpft_nmass_veg/=(double)stand.npatch();
				standpft_anpp/=(double)stand.npatch();
				standpft_lai/=(double)stand.npatch();
				standpft_densindiv_total/=(double)stand.npatch();
				standpft_aiso/=(double)stand.npatch();
				standpft_amon/=(double)stand.npatch();
				standpft_nuptake/=(double)stand.npatch();
				standpft_vmaxnlim/=(double)stand.npatch();
				heightindiv_total/=(double)stand.npatch();

				if (!negligible(standpft_cmass_leaf))
					standpft_vmaxnlim /= standpft_cmass_leaf;

				//Update landcover totals
				landcover_cmass[stand.landcover]+=standpft_cmass*stand.get_landcover_fraction();
				landcover_cmass_leaf[stand.landcover]+=standpft_cmass_leaf*stand.get_landcover_fraction();
				landcover_nmass_leaf[stand.landcover]+=standpft_nmass_leaf*stand.get_landcover_fraction();

				landcover_cmass_sap[stand.landcover] += standpft_cmass_sap*stand.get_landcover_fraction();
				landcover_cmass_heart[stand.landcover] += standpft_cmass_heart*stand.get_landcover_fraction();
				landcover_cmass_root[stand.landcover] += standpft_cmass_root*stand.get_landcover_fraction();
				landcover_cmass_debt[stand.landcover] += standpft_cmass_debt*stand.get_landcover_fraction();

				landcover_cmass_veg[stand.landcover]+=standpft_cmass_veg*stand.get_landcover_fraction();
				landcover_nmass_veg[stand.landcover]+=standpft_nmass_veg*stand.get_landcover_fraction();
				landcover_anpp[stand.landcover]+=standpft_anpp*stand.get_landcover_fraction();
				landcover_lai[stand.landcover]+=standpft_lai*stand.get_landcover_fraction();
				landcover_densindiv_total[stand.landcover]+=standpft_densindiv_total*stand.get_landcover_fraction();
				landcover_aiso[stand.landcover]+=standpft_aiso*stand.get_landcover_fraction();
				landcover_amon[stand.landcover]+=standpft_amon*stand.get_landcover_fraction();
				landcover_nuptake[stand.landcover]+=standpft_nuptake*stand.get_landcover_fraction();
				landcover_vmaxnlim[stand.landcover]+=standpft_vmaxnlim*stand.get_landcover_fraction();

				//Update pft totals
				gcpft_cmass+=standpft_cmass;
				gcpft_nmass+=standpft_nmass;
				gcpft_cmass_leaf+=standpft_cmass_leaf;
				gcpft_nmass_leaf+=standpft_nmass_leaf;

				gcpft_cmass_sap += standpft_cmass_sap;
				gcpft_cmass_heart += standpft_cmass_heart;
				gcpft_cmass_root += standpft_cmass_root;
				gcpft_cmass_debt += standpft_cmass_debt;

				gcpft_cmass_veg+=standpft_cmass_veg;
				gcpft_nmass_veg+=standpft_nmass_veg;
				gcpft_anpp+=standpft_anpp;
				gcpft_lai+=standpft_lai;
				gcpft_densindiv_total+=standpft_densindiv_total;
				gcpft_aiso+=standpft_aiso;
				gcpft_amon+=standpft_amon;
				gcpft_nuptake+=standpft_nuptake;
				gcpft_vmaxnlim+=standpft_vmaxnlim;

				if (vegmode==COHORT || vegmode==INDIVIDUAL)
					for (c=0;c<nclass;c++)
						gcpft_densindiv_ageclass[c] += standpft_densindiv_ageclass[c];

				// Update gridcell totals
				double fraction_of_gridcell = stand.get_gridcell_fraction();
				
				cmass_gridcell+=standpft_cmass*fraction_of_gridcell;
				nmass_gridcell+=standpft_nmass*fraction_of_gridcell;

				cmass_leaf_gridcell+=standpft_cmass_leaf*fraction_of_gridcell;
				nmass_leaf_gridcell+=standpft_nmass_leaf*fraction_of_gridcell;

				cmass_sap_gridcell += standpft_cmass_sap*fraction_of_gridcell;
				cmass_heart_gridcell += standpft_cmass_heart*fraction_of_gridcell;
				cmass_root_gridcell += standpft_cmass_root*fraction_of_gridcell;
				cmass_debt_gridcell += standpft_cmass_debt*fraction_of_gridcell;

				cmass_veg_gridcell+=standpft_cmass_veg*fraction_of_gridcell;
				nmass_veg_gridcell+=standpft_nmass_veg*fraction_of_gridcell;
				anpp_gridcell+=standpft_anpp*fraction_of_gridcell;
				lai_gridcell+=standpft_lai*fraction_of_gridcell;
				dens_gridcell+=standpft_densindiv_total*fraction_of_gridcell;
				aiso_gridcell+=standpft_aiso*fraction_of_gridcell;
				amon_gridcell+=standpft_amon*fraction_of_gridcell;
				nuptake_gridcell+=standpft_nuptake*fraction_of_gridcell;
				vmaxnlim_gridcell+=standpft_vmaxnlim*standpft_cmass_leaf*fraction_of_gridcell;

				// Graphical output every 10 years
				// (Windows shell only - "plot" statements have no effect otherwise)
				if (!(date.year%10)) {
					plot("C mass [kg C/m2]",pft.name,date.year,gcpft_cmass);
					plot("NPP [kg C/m2/yr]",pft.name,date.year,gcpft_anpp);
					plot("LAI [m2/m2]",pft.name,date.year,gcpft_lai);
					plot("dens [indiv/ha]",pft.name,date.year,gcpft_densindiv_total*m2toha);
					if (gcpft_cmass_leaf > 0.0 && ifnlim) {
						plot("vmax nitrogen lim [dimless]",pft.name,date.year,gcpft_vmaxnlim);
						plot("leaf C:N ratio [kg C/kg N]",pft.name,date.year,gcpft_cmass_leaf/gcpft_nmass_leaf);
					}
				}
				gridcell.nextobj();
			}//End of loop through stands

			// Print PFT sums to files

			double gcpft_cton_leaf = limited_cton(gcpft_cmass_leaf, gcpft_nmass_leaf);
			double gcpft_cton_veg = limited_cton(gcpft_cmass_veg, gcpft_nmass_veg);
			
			out.add_value(out_cmass,     gcpft_cmass);
			out.add_value(out_cmass_leaf, gcpft_cmass_leaf);
			out.add_value(out_cmass_sap, gcpft_cmass_sap);
			out.add_value(out_cmass_heart, gcpft_cmass_heart);
			out.add_value(out_cmass_root, gcpft_cmass_root);
			out.add_value(out_cmass_debt, gcpft_cmass_debt);
			out.add_value(out_anpp,      gcpft_anpp);
			out.add_value(out_dens,	     gcpft_densindiv_total);
			out.add_value(out_lai,       gcpft_lai);
			out.add_value(out_aiso,      gcpft_aiso);
			out.add_value(out_amon,      gcpft_amon);
			out.add_value(out_cton_leaf, gcpft_cton_leaf);
			out.add_value(out_cton_veg,  gcpft_cton_veg);
			out.add_value(out_vmaxnlim,  gcpft_vmaxnlim);
			out.add_value(out_nuptake,   gcpft_nuptake * m2toha);	

			// print species heights
			double height = 0.0;
			if (gcpft_densindiv_total > 0.0)
				height = heightindiv_total / gcpft_densindiv_total;
			
			out.add_value(out_speciesheights, height);

			// Write next part of record in binary daily LAI file
			// <pft.id> 0 <lai> <day> <lai> <day> <lai> ... <day> <lai> <ndayyear-1> <lai>
		
			if (out_dlai) fwrite(&pft.id,sizeof(pft.id),1,out_dlai);
			
			for (int d = 0; d < date.ndayyear; d++) {

				 // Sort of an ugly way of dealing with this now, but we get the dlai value
				 // (on stand basis) from _the_ stand where the PFT is active.
				 double dlai;
				 if (gridcell[0].pft[pft.id].active) {
					 dlai = gridcell[0].pft[pft.id].dlaiphen[d];
					 if (gridcell[1].pft[pft.id].active) {
						 fail("Failed assumption in outannual: pft %s is active in two stands\n", (char*)pft.name);
					 }
				 }
				 else {
					 if (!gridcell[1].pft[pft.id].active) {
						fail("Failed assumption in outannual: pft %s isn't active in any stand\n", (char*)pft.name);
					 }
					 dlai = gridcell[1].pft[pft.id].dlaiphen[d];
				 }

				 if (d == 0 || dlai != last || d == date.ndayyear-1) {
					  if (out_dlai) fwrite(&d, sizeof(d), 1, out_dlai);
					  last = dlai;
					  if (out_dlai) fwrite(&last, sizeof(last), 1, out_dlai);
				 }
			}


			pftlist.nextobj();
		
		} // *** End of PFT loop ***

		flux_veg = flux_soil = flux_fire = flux_est = flux_charvest = 0.0;

		// guess2008 - carbon pools
		c_litter = c_fast = c_slow = c_harv_slow = 0.0;

		surfsoillitterc = surfsoillittern = cwdc = cwdn = centuryc = centuryn = n_litter = n_harv_slow = availn = 0.0;
		andep_gridcell = anfert_gridcell = anmin_gridcell = animm_gridcell = anfix_gridcell = 0.0;
		n_org_leach_gridcell = n_min_leach_gridcell = 0.0;
		flux_nh3 = flux_no = flux_no2 = flux_n2o = flux_nsoil = flux_ntot = flux_nharvest = 0.0;

		// Sum C fluxes, dead C pools and runoff across patches

		gridcell.firstobj();

		// Loop through Stands
		while (gridcell.isobj) {
			Stand& stand = gridcell.getobj();
			stand.firstobj();

			//Loop through Patches
			while (stand.isobj) {
				Patch& patch = stand.getobj();

				double to_gridcell_average = stand.get_gridcell_fraction() / (double)stand.npatch();

				flux_veg+=-patch.fluxes.get_annual_flux(Fluxes::NPP)*to_gridcell_average;
				flux_soil+=patch.fluxes.get_annual_flux(Fluxes::SOILC)*to_gridcell_average;
				flux_fire+=patch.fluxes.get_annual_flux(Fluxes::FIREC)*to_gridcell_average;
				flux_est+=patch.fluxes.get_annual_flux(Fluxes::ESTC)*to_gridcell_average;
				flux_charvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTC)*to_gridcell_average;
				flux_nharvest+=patch.fluxes.get_annual_flux(Fluxes::HARVESTN)*to_gridcell_average;
				flux_nh3+=patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE)*to_gridcell_average;
				flux_no+=patch.fluxes.get_annual_flux(Fluxes::NO_FIRE)*to_gridcell_average;
				flux_no2+=patch.fluxes.get_annual_flux(Fluxes::NO2_FIRE)*to_gridcell_average;
				flux_n2o+=patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE)*to_gridcell_average;
				flux_nsoil+=patch.fluxes.get_annual_flux(Fluxes::N_SOIL)*to_gridcell_average;	
				flux_ntot+=(patch.fluxes.get_annual_flux(Fluxes::NH3_FIRE) + 
				           patch.fluxes.get_annual_flux(Fluxes::NO_FIRE) + 
				           patch.fluxes.get_annual_flux(Fluxes::NO2_FIRE) +
				           patch.fluxes.get_annual_flux(Fluxes::N2O_FIRE) +
				           patch.fluxes.get_annual_flux(Fluxes::N_SOIL)) * to_gridcell_average;
				
				c_fast+=patch.soil.cpool_fast*to_gridcell_average;
				c_slow+=patch.soil.cpool_slow*to_gridcell_average;

				// Sum all litter
				for (int q=0;q<npft;q++) {
					Patchpft& patchpft = patch.pft[q];
					c_litter += (patchpft.litter_leaf + patchpft.litter_root + patchpft.litter_sap + patchpft.litter_heart + patchpft.litter_repr)  * to_gridcell_average;
					n_litter += (patchpft.nmass_litter_leaf + patchpft.nmass_litter_root + patchpft.nmass_litter_sap + patchpft.nmass_litter_heart) * to_gridcell_average;
				}

				//Sum slow pools of harvested products
				if(run_landcover && ifslowharvestpool)
				{
					for (int q=0;q<npft;q++) 
					{
						Patchpft& patchpft=patch.pft[q];
						c_harv_slow+=patchpft.harvested_products_slow*to_gridcell_average;
						n_harv_slow+=patchpft.harvested_products_slow_nmass*to_gridcell_average;
					}
				}

				surfrunoff_gridcell+=patch.asurfrunoff*to_gridcell_average;
				drainrunoff_gridcell+=patch.adrainrunoff*to_gridcell_average;
				baserunoff_gridcell+=patch.abaserunoff*to_gridcell_average;
				runoff_gridcell+=patch.arunoff*to_gridcell_average;
	
				// Fix output error here: previous firert is sumed directly through patch which 
				// is not correct. Instead, it should be the fire probability summing up here, and then
				// up-scale to gridcell level - mcw
				// Fire probability: not distinguish stand here, assume fire happen at one stand
				// can spread to other stand, fires at stand are summed up
				if (!iffire || patch.fireprob < 0.001)
					fireprob_gridcell = 0.0; // Set a limit of 1000 years
				else
					fireprob_gridcell += patch.fireprob * to_gridcell_average;


				andep_gridcell += stand.climate.andep / (double)stand.npatch();
				anfert_gridcell += stand.climate.anfert / (double)stand.npatch();
				anmin_gridcell += patch.soil.anmin / (double)stand.npatch();
				animm_gridcell += patch.soil.animmob / (double)stand.npatch();
				anfix_gridcell += patch.soil.anfix / (double)stand.npatch();
				n_min_leach_gridcell += patch.soil.aminleach / (double)stand.npatch();
				n_org_leach_gridcell += patch.soil.aorgleach / (double)stand.npatch();
				availn += (patch.soil.nmass_avail + patch.soil.snowpack_nmass)
				          / (double)stand.npatch();

				for (int r = 0; r < NSOMPOOL-1; r++) {
					if (patch.soil.sompool[r].nmass > 0.0) {
						if(r == SURFMETA || r == SURFSTRUCT || r == SOILMETA || r == SOILSTRUCT){
							surfsoillitterc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							surfsoillittern += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
						else if (r == SURFFWD || r == SURFCWD) {
							cwdc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							cwdn += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
						else {	
							centuryc += patch.soil.sompool[r].cmass / (double)stand.npatch();
							centuryn += patch.soil.sompool[r].nmass / (double)stand.npatch();
						}
					}
				}

				// Monthly output variables

				for (m=0;m<12;m++) {
					maet[m] += patch.maet[m]*to_gridcell_average;
					mpet[m] += patch.mpet[m]*to_gridcell_average;
					mevap[m] += patch.mevap[m]*to_gridcell_average;
					mintercep[m] += patch.mintercep[m]*to_gridcell_average;
					mrunoff[m] += patch.mrunoff[m]*to_gridcell_average;
					mrh[m] += patch.fluxes.get_monthly_flux(Fluxes::SOILC, m)*to_gridcell_average;
					mwcont_upper[m] += patch.soil.mwcont[m][0]*to_gridcell_average;
					mwcont_lower[m] += patch.soil.mwcont[m][1]*to_gridcell_average;
					mgpp[m] += patch.fluxes.get_monthly_flux(Fluxes::GPP, m)*to_gridcell_average;
					mra[m] += patch.fluxes.get_monthly_flux(Fluxes::RA, m)*to_gridcell_average;

					miso[m]+=patch.fluxes.get_monthly_flux(Fluxes::ISO, m)*to_gridcell_average;
					mmon[m]+=patch.fluxes.get_monthly_flux(Fluxes::MON, m)*to_gridcell_average;
				}

				// Calculate monthly NPP and LAI
				/* Not necessary to story monthly LAI here
				   Since these can be found from RCA (daily data there)  -mcw
				Vegetation& vegetation = patch.vegetation;

				vegetation.firstobj();
				while (vegetation.isobj) {
					Individual& indiv = vegetation.getobj();

					// guess2008 - alive check added
					if (indiv.id != -1 && indiv.alive) {

						for (m=0;m<12;m++) {
							mlai[m] += indiv.mlai[m] * to_gridcell_average;
						}

					} // alive?

					vegetation.nextobj();

				} // while/vegetation loop
				*/
				stand.nextobj();
			} // patch loop
			gridcell.nextobj();
		} // stand loop



		// Fire return time
		if (!iffire || fireprob_gridcell < 0.001)
			firert_gridcell = 1000.0; // Set a limit of 1000 years
		else
			firert_gridcell = 1.0 / fireprob_gridcell;

		// In contrast to annual NEE, monthly NEE does not include fire
		// or establishment fluxes
		for (m=0;m<12;m++) {
			mnpp[m] = mgpp[m] - mra[m];
			mnee[m] = mnpp[m] - mrh[m];
		}

		// Print gridcell totals to files

		// Determine total leaf C:N ratio
		double cton_leaf_gridcell = limited_cton(cmass_leaf_gridcell, nmass_leaf_gridcell);
		
		// Determine total vegetation C:N ratio
		double cton_veg_gridcell = limited_cton(cmass_veg_gridcell, nmass_veg_gridcell);

		// Determine total vmax nitrogen limitation
		if (cmass_leaf_gridcell > 0.0) {
			vmaxnlim_gridcell /= cmass_leaf_gridcell;
		}

		out.add_value(out_cmass,  cmass_gridcell);
		out.add_value(out_cmass_leaf, cmass_leaf_gridcell);
		out.add_value(out_cmass_sap, cmass_sap_gridcell);
		out.add_value(out_cmass_heart, cmass_heart_gridcell);
		out.add_value(out_cmass_root, cmass_root_gridcell);
		out.add_value(out_cmass_debt, cmass_debt_gridcell);

		out.add_value(out_anpp,   anpp_gridcell);
		out.add_value(out_dens,   dens_gridcell);
		out.add_value(out_lai,    lai_gridcell);
		out.add_value(out_firert, firert_gridcell);
		out.add_value(out_runoff, surfrunoff_gridcell);
		out.add_value(out_runoff, drainrunoff_gridcell);
		out.add_value(out_runoff, baserunoff_gridcell);
		out.add_value(out_runoff, runoff_gridcell);
		out.add_value(out_aiso,   aiso_gridcell);
		out.add_value(out_amon,   amon_gridcell);

		out.add_value(out_cton_leaf, cton_leaf_gridcell);
		out.add_value(out_vmaxnlim,  vmaxnlim_gridcell);
		out.add_value(out_cton_veg,  cton_veg_gridcell);
		out.add_value(out_nuptake,   nuptake_gridcell * m2toha);

		out.add_value(out_nsources, andep_gridcell * m2toha);
		out.add_value(out_nsources, anfix_gridcell * m2toha);
		out.add_value(out_nsources, anfert_gridcell * m2toha);
		out.add_value(out_nsources, (andep_gridcell + anfix_gridcell + anfert_gridcell) * m2toha);
		out.add_value(out_nsources, anmin_gridcell * m2toha);
		out.add_value(out_nsources, animm_gridcell * m2toha);
		out.add_value(out_nsources, (anmin_gridcell - animm_gridcell) * m2toha);
		out.add_value(out_nsources, (anmin_gridcell - animm_gridcell + andep_gridcell + anfix_gridcell + anfert_gridcell) * m2toha);

		if (run_landcover) {
			for(int i=0;i<NLANDCOVERTYPES;i++) {
				if(run[i]) {
					out.add_value(out_cmass, landcover_cmass[i]);
					out.add_value(out_cmass_leaf, landcover_cmass_leaf[i]);
					out.add_value(out_cmass_sap, landcover_cmass_sap[i]);
					out.add_value(out_cmass_heart, landcover_cmass_heart[i]);
					out.add_value(out_cmass_root, landcover_cmass_root[i]);
					out.add_value(out_cmass_debt, landcover_cmass_debt[i]);

					out.add_value(out_anpp,  landcover_anpp[i]);
					out.add_value(out_dens,  landcover_densindiv_total[i]);
					out.add_value(out_lai,   landcover_lai[i]);
					out.add_value(out_aiso,  landcover_aiso[i]);
					out.add_value(out_amon,  landcover_amon[i]);

					double landcover_cton_leaf = limited_cton(landcover_cmass_leaf[i], landcover_nmass_leaf[i]);
					double landcover_cton_veg = limited_cton(landcover_cmass_veg[i], landcover_nmass_veg[i]);

					if (landcover_cmass_leaf[i] > 0.0) {
						landcover_vmaxnlim[i] /= landcover_cmass_leaf[i];
					}

					out.add_value(out_cton_leaf, landcover_cton_leaf);
					out.add_value(out_cton_veg,  landcover_cton_veg);
					out.add_value(out_vmaxnlim,  landcover_vmaxnlim[i]);
					out.add_value(out_nuptake,   landcover_nuptake[i] * m2toha);
				}
			}
		}

		// Print monthly output variables
		for (m=0;m<12;m++) {
			 out.add_value(out_mnpp,         mnpp[m]);
			 out.add_value(out_mlai,         mlai[m]);
			 out.add_value(out_mgpp,         mgpp[m]);
			 out.add_value(out_mra,          mra[m]);
			 out.add_value(out_maet,         maet[m]);
			 out.add_value(out_mpet,         mpet[m]);
			 out.add_value(out_mevap,        mevap[m]);
			 out.add_value(out_mrunoff,      mrunoff[m]);
			 out.add_value(out_mintercep,    mintercep[m]);
			 out.add_value(out_mrh,          mrh[m]);
			 out.add_value(out_mnee,         mnee[m]);
			 out.add_value(out_mwcont_upper, mwcont_upper[m]);
			 out.add_value(out_mwcont_lower, mwcont_lower[m]);
			 out.add_value(out_miso,         miso[m]);
			 out.add_value(out_mmon,         mmon[m]);
		}


		// Print land use to file
	
		const LandCoverInfo* lci = gridcell.landcoverinfo.get();
		double fgrid_land = lci->get_fgrid_land();
		double fland_crop = lci->get_fland_crop(cyear);
		double fland_bare = lci->get_fland_bare(cyear);

		out.add_value(out_landuse, fgrid_land);
		out.add_value(out_landuse, 1.0-fland_crop-fland_bare);
		out.add_value(out_landuse, fland_crop);
		out.add_value(out_landuse, fland_bare);
	
		// Write driver data during control period

		if (out_driver && ifcontrolperiod && date.year>=first_control_year && date.year<=last_control_year) {

			 double lon = gridcell.get_lon();
			 double lat = gridcell.get_lat();
			 fwrite(&lon,sizeof(double),1,out_driver);
			 fwrite(&lat,sizeof(double),1,out_driver);
			 fwrite(&gridcell.id,sizeof(int),1,out_driver);	
			 fwrite(&date.year,sizeof(int),1,out_driver);
			 fwrite(gridcell[0].laiphen_grass,sizeof(float),366,out_driver);
			 fwrite(gridcell[1].laiphen_grass,sizeof(float),366,out_driver);
			 fwrite(gridcell[1].laiphen_conifer,sizeof(float),366,out_driver);
			 fwrite(gridcell[1].laiphen_broadleaf,sizeof(float),366,out_driver);
			 fwrite(&gridcell[1].laimax_conifer,sizeof(float),1,out_driver);
			 fwrite(&gridcell[1].laimax_broadleaf,sizeof(float),1,out_driver);	
		}


		// Graphical output every 10 years
		// (Windows shell only - no effect otherwise)

		if (!(date.year%10)) {
			gridcell.firstobj();
			if(gridcell.isobj)	//Fixed bug here if no stands were present.
			{
				Stand& stand=gridcell.getobj();
				plot("C flux [kg C/m2/yr]","flux_veg",  date.year, flux_veg);
				plot("C flux [kg C/m2/yr]","flux_soil", date.year, flux_soil);
				plot("C flux [kg C/m2/yr]","flux_fire", date.year, flux_fire);
				plot("C flux [kg C/m2/yr]","flux_est",  date.year, flux_est);
				plot("C flux [kg C/m2/yr]","NEE",       date.year, flux_veg + flux_soil + flux_fire + flux_est);

				if (!ifcentury) {
					plot("Soil C [kg C/m2]","slow", date.year, stand[0].soil.cpool_slow);
					plot("Soil C [kg C/m2]","fast", date.year, stand[0].soil.cpool_fast);
				}
				else {
					plot("N flux (kg N/ha/yr)","Fix",   date.year, -anfix_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Dep",   date.year, -andep_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Fert",  date.year, -anfert_gridcell * m2toha);
					plot("N flux (kg N/ha/yr)","Leach", date.year, (n_min_leach_gridcell + n_org_leach_gridcell) * m2toha);
					plot("N flux (kg N/ha/yr)","Flux",  date.year, flux_ntot * m2toha);

					plot("N mineralization [kg N/ha/yr]","N", date.year, (anmin_gridcell - animm_gridcell) * m2toha);

					plot("Soil C [kg C/m2]","fine litter",   date.year, surfsoillitterc);
					plot("Soil C [kg C/m2]","coarse litter", date.year, cwdc);
					plot("Soil C [kg C/m2]","soil",          date.year, centuryc); 
					plot("Soil C [kg C/m2]","total",         date.year, surfsoillitterc + cwdc + centuryc); 

					plot("Soil N [kg N/ha]","fine litter",   date.year, surfsoillittern);
					plot("Soil N [kg N/ha]","coarse litter", date.year, cwdn);
					plot("Soil N [kg N/ha]","soil",          date.year, centuryn); 
					plot("Soil N [kg N/ha]","total",         date.year, surfsoillittern + cwdn + centuryn); 
				}
			}
		}

		// Write fluxes to file

		out.add_value(out_cflux, flux_veg);
		out.add_value(out_cflux, flux_soil);
		out.add_value(out_cflux, flux_fire);
		out.add_value(out_cflux, flux_est);
		if (run_landcover) {
			 out.add_value(out_cflux, flux_charvest);
		}
		out.add_value(out_cflux, flux_veg + flux_soil + flux_fire + flux_est + flux_charvest);

		out.add_value(out_nflux, -andep_gridcell * m2toha);
		out.add_value(out_nflux, -anfix_gridcell * m2toha);
		out.add_value(out_nflux, -anfert_gridcell * m2toha);
		out.add_value(out_nflux, flux_ntot * m2toha);
		out.add_value(out_nflux, (n_min_leach_gridcell + n_org_leach_gridcell) * m2toha);
		if (run_landcover) {
			 out.add_value(out_nflux, flux_nharvest);
		}
		out.add_value(out_nflux, (flux_nharvest + flux_ntot + n_min_leach_gridcell + n_org_leach_gridcell - (andep_gridcell + anfix_gridcell + anfert_gridcell)) * m2toha);

		out.add_value(out_cpool, cmass_gridcell);
		out.add_value(out_cpool, c_litter);
		if (!ifcentury) {
			out.add_value(out_cpool, c_fast);
			out.add_value(out_cpool, c_slow);
		}
		else {
			out.add_value(out_cpool, surfsoillitterc);
			out.add_value(out_cpool, cwdc);
			out.add_value(out_cpool, centuryc);
		}
		
		if (run_landcover && ifslowharvestpool) {
			out.add_value(out_cpool, c_harv_slow);
		}

		// Calculate total cpool, starting with cmass and litter...
		double cpool_total = cmass_gridcell + c_litter;

		// Add SOM pools
		if (!ifcentury) {
			cpool_total += c_fast + c_slow;
		}
		else {
			cpool_total += centuryc + surfsoillitterc + cwdc;
		}

		// Add slow harvest pool if needed
		if (run_landcover && ifslowharvestpool) {
			cpool_total += c_harv_slow;
		}

		out.add_value(out_cpool, cpool_total);

		if (ifcentury) {
			out.add_value(out_npool, nmass_gridcell);
			out.add_value(out_npool, n_litter);
			out.add_value(out_npool, surfsoillittern);
			out.add_value(out_npool, cwdn);
			out.add_value(out_npool, centuryn + availn);

			if(run_landcover && ifslowharvestpool) {
				out.add_value(out_npool, n_harv_slow);
				out.add_value(out_npool, (nmass_gridcell + n_litter + surfsoillittern + cwdn + centuryn + availn + n_harv_slow));
			}
			else {
				out.add_value(out_npool, (nmass_gridcell + n_litter + surfsoillittern + cwdn + centuryn + availn));
			}
		}

		out.add_value(out_ngases, flux_nh3   * m2toha);
		out.add_value(out_ngases, flux_no    * m2toha);
		out.add_value(out_ngases, flux_no2   * m2toha);
		out.add_value(out_ngases, flux_n2o   * m2toha);
		out.add_value(out_ngases, flux_nsoil * m2toha);
		out.add_value(out_ngases, flux_ntot  * m2toha);

		// Output of age structure (Windows shell only - no effect otherwise)

		if (vegmode==COHORT || vegmode==INDIVIDUAL) {

			if (!(date.year%20) && date.year<2000) {
			
				resetwindow("Age structure [yr]");

				pftlist.firstobj();
				while (pftlist.isobj) {
					Pft& pft=pftlist.getobj();

					if (pft.lifeform==TREE) {

						Gridcellpft& gridcellpft=gridcell.pft[pft.id];

						for (c=0;c<nclass;c++)
							plot("Age structure [yr]",pft.name,
								c * estinterval + estinterval / 2,
								gcpft_densindiv_ageclass[c] / (double)npatch);
					}
					
					pftlist.nextobj();
				}
			}
		}
		// RCA: Flush these since guess isn't terminated properly (termio doesn't get called)
		fflush(out_dlai);
		fflush(out_driver);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
// TERMIO
// Called at end of model run (i.e. following simulation of all stands)

void termio() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.
	delete output_channel;

	// Clean up
	gridlist.killall();
}

#endif // USE_CRU_IO


///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
// Lamarque, J.-F., Kyle, G. P., Meinshausen, M., Riahi, K., Smith, S. J., Van Vuuren, 
//   D. P., Conley, A. J. & Vitt, F. 2011. Global and regional evolution of short-lived
//   radiatively-active gases and aerosols in the Representative Concentration Pathways. 
//   Climatic Change, 109, 191-212.
// Nakai, T., Sumida, A., Kodama, Y., Hara, T., Ohta, T. (2010). A comparison between
//   various definitions of forest stand height and aerodynamic canopy height.
//   Agricultural and Forest Meteorology, 150(9), 1225-1233
