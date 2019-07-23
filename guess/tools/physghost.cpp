/* 
 * A simulator for running GUESS without RCA.
 *
 * Gets climate data from CRU binaries and calls GUESS 
 * through the same interface as RCA uses.
 *
 */

#include "guessmain.h"
#include "guess.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <map>
#include <utility>
#include <vector>
#include <fstream>

#include "cru_1901_2006.h"
#include "RCAData.h"
#include "driver.h"

// The precision of the floating point numbers used in the RCA/GUESS interface
// is either single or double, depending on whether RCA is compiled with
// the SINGLE_PRECISION flag set or not. So the simulator needs to take this
// into account as well when calling GUESS.
#ifdef SINGLE_PRECISION
typedef float real;
#else
typedef double real;
#endif

#define CRU_PATH "D:\\DATA\\CRU_TS3.0\\cru_1901_2006.bin"
#define RCA_PATH "HadData.bin"


/* 
* Data definition 
*/
// Represents one grid cell
struct Coordinate {

	 Coordinate(double latitude, double longitude)
				: id(0), lat_d(-999), lon_d(-999), avail_d(true),
				  lat(latitude),
				  lon(longitude) {
	 }

	 int id;

	 double lat;
	 double lon;
	 double lat_d;
	 double lon_d;
	 bool avail_d;
};

// TODO: To avoid duplicated code, better use the same definition of Spinup_data 
// in guessio_cru.cpp
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
		nyear = nyear_loc;
		havedata = false;
		data = new double[nyear * 12];
		if (!data) fail("Spinup_data::Spinup_data: out of memory");
		thisyear = 0;
		havedata = true;
		reset_clim(); // guess2008
	}

	~Spinup_data() {
		if (havedata) delete[] data;
	}

	double& operator[](int month) {

		return data[thisyear * 12 + month];
	}

	// new function for phyghost
	double* get_clim() {

		return dataclim;
	}

	void nextyear() {
		if (thisyear == nyear - 1) thisyear = 0;
		else thisyear++;
	}

	void firstyear() {
		thisyear = 0;
	}

	void get_data_from(double source[][12]) {

		int y, m;
		thisyear = 0; // guess2008 - ML bugfix
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				data[y * 12 + m] = source[y][m];
			}
		}
	}

	void get_data_from(double* source) {

		int m;
		thisyear = 0; // guess2008 - ML bugfix
		for (m = 0; m<nyear * 12; m++)
			data[m] = source[m];
	}


	// guess2008 - NEW METHODS 

	void reset_clim() {
		for (int ii = 0; ii < 12; ii++) dataclim[ii] = 0.0;
	}


	void make_clim() {

		reset_clim(); // Always reset before calculating

		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				dataclim[m] += data[y * 12 + m] / (double)nyear;
			}
		}
	}


	bool extract_data(double source[][12], const int& startyear, const int& endyear) {

		// Populate data with data from the middle of source. 
		// Condition: endyear - startyear + 1 == nyear
		// if startyear == 1 and endyear == 30 then this function is identical to get_data_from above.

		if (endyear < startyear) return false;
		if (endyear - startyear + 1 == nyear) {

			int y, m;
			for (y = startyear - 1; y<endyear; y++) {
				for (m = 0; m<12; m++) {
					data[(y - (startyear - 1)) * 12 + m] = source[y][m];
				}
			}

		}
		else return false;

		return true;
	}


	void adjust_data(double anom[12], bool additive) {

		// Adjust the spinup data to the conditions prevailing at a particular time, as given by 
		// the (additive or multiplicative) anomalies in anom 
		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				if (additive)
					data[y * 12 + m] += anom[m];
				else
					data[y * 12 + m] *= anom[m];
			}
		}

	}


	// Replace interannual data with the period's climatology.
	void use_clim_data() {

		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				data[y * 12 + m] = dataclim[m];
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

		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				data[y * 12 + m] = dataclim[m] + (data[y * 12 + m] - dataclim[m]) * factor;
			}
		}
	}


	void limit_data(double minval, double maxval) {

		// Limit data to a range
		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				if (data[y * 12 + m] < minval) data[y * 12 + m] = minval;
				if (data[y * 12 + m] > maxval) data[y * 12 + m] = maxval;
			}
		}

	}


	void set_min_val(const double& oldval, const double& newval) {

		// Change values < oldval to newval
		int y, m;
		for (y = 0; y<nyear; y++) {
			for (m = 0; m<12; m++) {
				if (data[y * 12 + m] < oldval) data[y * 12 + m] = newval;
			}
		}

	}

	// guess2008 - END OF NEW METHODS


	void detrend_data() {

		int y, m;
		double a, b, anomaly;
		double* annual_mean = new double[nyear];
		double* year_number = new double[nyear];

		if (!annual_mean || !year_number)

			fail("Spinup_driver::detrend_data: out of memory");

		for (y = 0; y<nyear; y++) {
			annual_mean[y] = 0.0;
			for (m = 0; m<12; m++) annual_mean[y] += data[y * 12 + m];
			annual_mean[y] /= 12.0;
			year_number[y] = y;
		}

		regress(year_number, annual_mean, nyear, a, b);

		for (y = 0; y<nyear; y++) {
			anomaly = b*(double)y;
			for (m = 0; m<12; m++)
				data[y * 12 + m] -= anomaly;
		}

		// guess2008 - added [] - Clean up
		delete[] annual_mean;
		delete[] year_number;
	}
};

// number of year for the spin-up period
const int NYEAR_SPINUP_DATA = 30;
// Spinup data sets for one  grid cell
Spinup_data ghost_msun(NYEAR_SPINUP_DATA);
Spinup_data ghost_mpar(NYEAR_SPINUP_DATA);
Spinup_data ghost_mprec(NYEAR_SPINUP_DATA);

Spinup_data ghost_mtemp_for(NYEAR_SPINUP_DATA);
Spinup_data ghost_mtemp_soil_for(NYEAR_SPINUP_DATA);
Spinup_data ghost_mwcont_upper_for(NYEAR_SPINUP_DATA);
Spinup_data ghost_mwcont_lower_for(NYEAR_SPINUP_DATA);

Spinup_data ghost_mtemp_opl(NYEAR_SPINUP_DATA);
Spinup_data ghost_mtemp_soil_opl(NYEAR_SPINUP_DATA);
Spinup_data ghost_mwcont_upper_opl(NYEAR_SPINUP_DATA);
Spinup_data ghost_mwcont_lower_opl(NYEAR_SPINUP_DATA);

struct climate_ghost {

	int id;

	double lat;
	double lon;

	double dsun[366];
	double dpar[366];
	double dprec[366];

	double dtemp_for[366];
	double dtemp_soil_for[366];
	double dwcont_upper_for[366];
	double dwcont_lower_for[366];

	double dtemp_opl[366];
	double dtemp_soil_opl[366];
	double dwcont_upper_opl[366];
	double dwcont_lower_opl[366];

};


/* 
* Routine and processes 
*/
// Calculates shortwave radiation
double swrad_(double lat, int julian_day, double pday, double sun) {

	// julian_day (0-364)
	// pday       = time of day (0-1 range from midnight to midnight)
	// sun        = proportion of bright sunshine

	const double QOO = 1360.0;
	const double PI = 3.1415927;
	const double BETA = 0.17;
	const double A = 107.0;
	const double B = 0.2;
	const double C = 0.25;
	const double D = 0.5;
	const double DEGTORAD = 0.01745329;
	const double FRADPAR = 0.5;

	double qo = QOO*(1.0 + 2.0*0.01675*cos(2.0*PI*(julian_day + 0.5) / 365.0));  // needs to consider leap year later : mcw
	double delta = -23.4*DEGTORAD*cos(2.0*PI*(julian_day + 10.5) / 365.0);
	double h = pday - 0.5;
	if (h<0.0) h += 1.0;

	double swrad = (C + D*sun)*(1.0 - BETA)*qo*(sin(lat*DEGTORAD)*sin(delta) +
		cos(lat*DEGTORAD)*cos(delta)*cos(h*2.0*PI));

	if (swrad<0.0) swrad = 0.0;

	return swrad;
}

/*
 * Reads in a text file with lon, lat pairs.
 *
 * Each line in the text file should contain two doubles,
 * longitude followed by latitude.
 */
std::vector<Coordinate> read_gridlist(const char* path) {
	 std::ifstream in(path, std::ifstream::in);
	 
	 std::vector<Coordinate> result;

	 double lon, lat;
	 while (in >> lon >> lat) {
		  result.push_back(Coordinate(lat, lon));
	 }
	 
	 return result;
}

bool searchrca(double dlon, double dlat, RCAData& rec) {

	// Object for retrieving RCA spinup data from archive file
	RCADataArchive ark;

	rec.lon = dlon;
	if (rec.lon)
		rec.lon = (int)(fabs(rec.lon)*1000.0 + 0.5) / 1000.0*fabs(rec.lon) / rec.lon;

	rec.lat = dlat;
	if (rec.lat)
		rec.lat = (int)(fabs(rec.lat)*1000.0 + 0.5) / 1000.0*fabs(rec.lat) / rec.lat;

	// Try block to catch any unexpected errors
	try {

		bool success = ark.open(RCA_PATH);

		if (success) {
			bool flag = ark.rewind();
			if (!flag) {
				ark.close(); // I.e. we opened it but we couldn't rewind
				return false;
			}
		}
		else
			return false;

		// Read the RCA data into the data struct
		if (!ark.getindex(rec)) {
			ark.close();
			return false;
		}

		// Convert PAR from kJ to J
		for (int y = 0; y < NYEAR_SPINUP_DATA; y++)
			for (int m = 0; m < 12; m++) {
				rec.mpar_for[y * 12 + m] *= 1.0e3;
				rec.mpar_open[y * 12 + m] *= 1.0e3;
			}

		// Close the archive
		ark.close();

		return true;

	}
	catch (...) {
		// Unknown error.
		//dprintf("Error in reading data from %s \n", RCA_PATH);
		fprintf(stderr,"Error in reading data from %s \n", RCA_PATH);
		return false;
	}

}


bool findnearstRCAData(Coordinate& coor, RCAData& rec) {

	double searchradius = 0.25;

	// First try the exact coordinate
	if (searchrca(coor.lon, coor.lat, rec)) {
		coor.lon_d = coor.lon;
		coor.lat_d = coor.lat;
		coor.avail_d = true;
		//dprintf("Exact record for gridcell (%0.3f,%0.3f) in RCA data!", coor.lon_d, coor.lat_d);
		fprintf(stdout,"\n * Exact record for gridcell (%0.3f,%0.3f) in RCA data!\n", coor.lon_d, coor.lat_d);
		return true;
	}

	if (searchradius == 0) {
		// Don't try to search
		return false;
	}

	// Search all coordinates in a square around (lon, lat), but first go down to
	// multiple of 0.5
	double center_lon = floor(coor.lon * 2) / 2;
	double center_lat = floor(coor.lat * 2) / 2;

	// Enumerate all coordinates within the square, place them in a vector of
	// pairs where the first element is distance from center to allow easy 
	// sorting.
	using std::pair;
	using std::make_pair;
	typedef pair<double, double> point;
	std::vector<pair<double, point> > search_points;

	const double STEP = 0.05; // for RCA
	const double EPS = 1e-15;

	for (double y = center_lon - searchradius; y <= center_lon + searchradius + EPS; y += STEP) {
		for (double x = center_lat - searchradius; x <= center_lat + searchradius + EPS; x += STEP) {
			double xdist = x - center_lat;
			double ydist = y - center_lon;
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

		if (searchrca(search_lon, search_lat, rec)) {
			coor.lon_d = search_lon;
			coor.lat_d = search_lat;
			coor.avail_d = true;
			return true;
			//dprintf("Nearest record (%0.3f,%0.3f) in RCA data for gridcell (%0.3f,%0.3f)!", 
				//     coor.lon_d, coor.lat_d, coor.lon, coor.lat);
			fprintf(stdout, "\n * Nearest record (%0.3f,%0.3f) in RCA data for gridcell (%0.3f,%0.3f)!\n", 
			     coor.lon_d, coor.lat_d, coor.lon, coor.lat);
		}
	}

	return false;

}
	// Make one set of climate data (- grid point level)
climate_ghost makeRCAclimate(RCAData rec, int num) {

	climate_ghost aset;

	aset.lon = rec.lon;
	aset.lat = rec.lat;

	aset.id = num;

	// Transfer to global Spinup_data objects for recycling throughout spinup period
	// RCAData should have already detrend and usually for the period 1961-1990
	ghost_mpar.get_data_from(rec.mpar_for); // mpar_for[360], 30 year, 12 month
	ghost_mprec.get_data_from(rec.mprec_for);

	ghost_mtemp_for.get_data_from(rec.mtemp_for);
	ghost_mtemp_soil_for.get_data_from(rec.mtemp_soil_for);
	ghost_mwcont_upper_for.get_data_from(rec.mwcont_upper_for);
	ghost_mwcont_lower_for.get_data_from(rec.mwcont_lower_for);

	ghost_mtemp_opl.get_data_from(rec.mtemp_open);
	ghost_mtemp_soil_opl.get_data_from(rec.mtemp_soil_open);
	ghost_mwcont_upper_opl.get_data_from(rec.mwcont_upper_open);
	ghost_mwcont_lower_opl.get_data_from(rec.mwcont_lower_open);

	// Make climatology
	ghost_mpar.make_clim();
	ghost_mprec.make_clim();

	ghost_mtemp_for.make_clim();
	ghost_mtemp_soil_for.make_clim();
	ghost_mwcont_upper_for.make_clim();
	ghost_mwcont_lower_for.make_clim();

	ghost_mtemp_opl.make_clim();
	ghost_mtemp_soil_opl.make_clim();
	ghost_mwcont_upper_opl.make_clim();
	ghost_mwcont_lower_opl.make_clim();

	// Interpolate to quasi-daily values
	interp_monthly_means(ghost_mpar.get_clim(), aset.dpar);
	interp_monthly_totals(ghost_mprec.get_clim(), aset.dprec);

	interp_monthly_means(ghost_mtemp_for.get_clim(), aset.dtemp_for);
	interp_monthly_means(ghost_mtemp_soil_for.get_clim(), aset.dtemp_soil_for);
	interp_monthly_means(ghost_mwcont_upper_for.get_clim(), aset.dwcont_upper_for);
	interp_monthly_means(ghost_mwcont_lower_for.get_clim(), aset.dwcont_lower_for);

	interp_monthly_means(ghost_mtemp_opl.get_clim(), aset.dtemp_opl);
	interp_monthly_means(ghost_mtemp_soil_opl.get_clim(), aset.dtemp_soil_opl);
	interp_monthly_means(ghost_mwcont_upper_opl.get_clim(), aset.dwcont_upper_opl);
	interp_monthly_means(ghost_mwcont_lower_opl.get_clim(), aset.dwcont_lower_opl);

	// to improve: interp_monthly_means doesn't handle leap year and assume 365 days for all years : mcw
	aset.dpar[365] = aset.dpar[364];
	aset.dprec[365] = aset.dprec[364];
	aset.dtemp_for[365] = aset.dtemp_for[364];
	aset.dtemp_soil_for[365] = aset.dtemp_soil_for[364];
	aset.dwcont_upper_for[365] = aset.dwcont_upper_for[364];
	aset.dwcont_lower_for[365] = aset.dwcont_lower_for[364];
	aset.dtemp_opl[365] = aset.dtemp_opl[364];
	aset.dtemp_soil_opl[365] = aset.dtemp_soil_opl[364];
	aset.dwcont_upper_opl[365] = aset.dwcont_upper_opl[364];
	aset.dwcont_lower_opl[365] = aset.dwcont_lower_opl[364];

	return aset;
}

// Make climate data vecotor
std::vector<climate_ghost> make_ghost(xtring path, xtring gtype) {

	std::vector<Coordinate> gridlist = read_gridlist(path);

	// Size of cdata may not equal to gridlist
	std::vector<climate_ghost> cdata;
	
	fprintf(stdout, "\n==> Building up RCA ghost climate data  ... \n ");
	int num = 0;
	for (size_t i=0; i < gridlist.size(); i++) {
		if (gtype.lower() == "cru") {

			// TODO:
			//cdata = make_ghost_cru(gridlist);

		}
		else if (gtype.lower() == "rca") {

			RCAData rec;

			// First, to check data avaliablity
			if (findnearstRCAData(gridlist[i], rec)) {
				//num++;
				// Data found and make daily climate forcing from monthly archive
				cdata.push_back(makeRCAclimate(rec,num));
				// fprintf(stdout, "\nRCA ghost climate data found for grid point (%3.4f, %3.4f) ... \n ", rec.lon, rec.lat);

			}
			else {

				gridlist[i].avail_d = false;
				//dprintf("Cannot find RCA ghost climate data for grid point (%0.3f,%0.3f) \n ", rec.lon, rec.lat);
				fprintf(stderr, "\nCannot find RCA ghost climate data for grid point (%0.3f,%0.3f) \n ", rec.lon, rec.lat);

			}

		}
	}

	//dprintf("Ghost climate data is maded ... \n");
	fprintf(stdout, "\nGhost climate data is maded! \n\n");
	return cdata;
}

int main(int argc, char** argv) {

	if (argc != 7) {
		fprintf(stderr, "Usage: physghost <gridlist> <start-year> <end-year> <ghost_type> <calendar_mode> <p_scalar>\n");
		exit(99);
	}

	xtring grid_file = argv[1];
	int start_year = atoi(argv[2]);
	int end_year = atoi(argv[3]);
	xtring ghost_type = argv[4];
	xtring caldr_mode = argv[5];     // calendar type: STANDARD, LEAPYEARS, FLAT_30
	double p_scalar = atof(argv[6]);  // mw130216: for p reduction adhoc test

	if (p_scalar < 0 || p_scalar > 1) p_scalar = 1;


	// Prepare ghost climate data
	std::vector<climate_ghost> climate_ghost = make_ghost(grid_file, ghost_type);

	FILE* out_guess = fopen("guess_feedback.out", "w");

	if (out_guess == 0) {

		//dprintf("Failed to open guess feedback file\n");
		fail("\nFailed to open guess feedback file\n");

	}
	else {

		fprintf(out_guess, "%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n",
			"Lon", "Lat", "Year", "Month", "frland", "frdecid", "frfor", "l_conif", "l_broad", "l_under",
			"lai_opl", "veg_opl");

	}

	// forcing data to GUESS
	int minutec = 0;       // simulator runs hourly timesteps
	int secondc = 0;       // dito...

	static const int standard_months[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int ndaymonth[12];

	real frland = 1.0;     // simulator always has 100% land
	real co2 = 360;		   // ignored by GUESS anyway (GUESS reads co2 from file)

	// forcing data from climate data set
	real temp_2m, temp_soil, insol, rain, snow;
	real soilw_surf, soilw_deep;

	for (int yearc = start_year; yearc <= end_year; yearc++) {
		int julian_day = 0;

		if (caldr_mode == "STANDARD") {
			for (int i = 0; i < 12; i++) {
				ndaymonth[i] = standard_months[i];
			}
		}
		else if (caldr_mode == "LEAPYEARS") {
			for (int i = 0; i < 12; i++) {
				ndaymonth[i] = standard_months[i];
			}
			if (Date::is_leap(yearc)) {
				ndaymonth[1] = 29;
			}
		}
		else if (caldr_mode == "FLAT_30") {
			for (int i = 0; i < 12; i++) {
				ndaymonth[i] = 30;
			}
		}
		else
			fail("Date::update_calendar: unknown calendartype");

		// Indeed, monthc and dayc start from 1 here, and will substract by 1
		// in the interface function FROUTINE(lai_guess), can improve this? - mcw
		for (int monthc = 1; monthc <= 12; monthc++) {

			for (int dayc = 1; dayc <= ndaymonth[monthc - 1]; dayc++) {

				for (int hourc = 0; hourc < 23; hourc++) {

					for (size_t pixel = 0; pixel < climate_ghost.size(); pixel++) {

						if (ghost_type.lower() == "cru") {

							couple_src = "CRU";

							// TODO: now it is not a good way to convert monthly data to hourly. 
							// need to reconsider the calculation here
							// only get data, do not read from file

							// CRU use sunshine
							// CRU data set doesn't contain soil water content, 

						}
						else if (ghost_type.lower() == "rca") {
							
							couple_src = "RCA_SP";   // this will set climate.instype as "PAR"

							// temp from RCAData is degC, need to convert to Kelvin for the interface
							temp_2m = climate_ghost[pixel].dtemp_for[julian_day] + K2degC; 

							// best not to use prescribed soil temperature with the simulator
							temp_soil = temp_2m;

							insol = climate_ghost[pixel].dpar[julian_day];
							rain = climate_ghost[pixel].dprec[julian_day] / (3600 * 24); // Convertto mm/s
							snow = 0;
							soilw_surf = climate_ghost[pixel].dwcont_upper_for[julian_day];
							soilw_deep = climate_ghost[pixel].dwcont_lower_for[julian_day];

						}

						rain*=p_scalar;
						snow*=p_scalar;

						// results from GUESS
						real frdecid = 0, frfor = 0;
						real lai_conif = 0, lai_broad = 0, lai_under = 0;
						real lai_opl = 0, veg_opl = 0;

						lai_guess_(&climate_ghost[pixel].id, &climate_ghost[pixel].lat, &climate_ghost[pixel].lon,
							&yearc, &monthc, &dayc, &hourc, &minutec, &secondc,
							&frland, &co2, 
							&temp_2m, &temp_2m,         // temperature for openland tile and forest tile (K)
							&temp_soil, &temp_soil, &temp_soil, &temp_soil, &temp_soil, // for openland tile (not used)
							&temp_soil, &temp_soil, &temp_soil, &temp_soil, &temp_soil, // for forest tile (not used)
							&soilw_surf, &soilw_deep,   // for openland tile (not used)
							&soilw_surf, &soilw_deep,   // for forest tile (not used)
							&insol,						// insolation for openland tile
							&insol,						// insolation for forest tile
							&rain, &snow,               // snow will add to total precip and split later according to air temp. 
							&veg_opl, &lai_opl,
							&lai_conif, &lai_broad, &lai_under,
							&frdecid, &frfor);

						if (dayc == 1) {
							fprintf(out_guess, "%8.2f%8.2f%8d%8d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
								climate_ghost[pixel].lon, climate_ghost[pixel].lat, yearc, monthc,
								frland, frdecid, frfor, lai_conif, lai_broad, lai_under,
								lai_opl, veg_opl);
						}
					} // gridlist loop
				} // hour loop
				julian_day++;
			} // day loop
		}
	}

	fclose(out_guess);
	return 0;
}
