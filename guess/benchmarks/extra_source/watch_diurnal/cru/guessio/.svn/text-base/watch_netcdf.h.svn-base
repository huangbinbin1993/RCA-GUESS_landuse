///////////////////////////////////////////////////////////////////////////////////////
/// \file watch_netcdf.h
/// \brief Code for reading in WATCH diurnal data from NetCDF files
///
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_WATCH_NETCDF_H
#define LPJ_GUESS_WATCH_NETCDF_H

#ifndef HAVE_NETCDF
#error "NetCDF library not found by build system!"
#endif

#include <netcdf.h>
#include <assert.h>

const int SUBDAILY = 8;
const int FIRST_WATCH_YEAR = 1901;

/// A lon,lat pair, used by load_gridlist function below
typedef std::pair<double, double> landpoint;

/// Help function to deal with status codes from NetCDF library
void handle_error(int status, const char* message) {
	if (status != NC_NOERR) {
		fail("NetCDF error: %s: %s\n", nc_strerror(status), message);
	}
}

/// Opens a NetCDF file and returns its id
int open_ncdf(const char* fname) {
	int netcdf_id;
	int status = nc_open(fname, NC_NOWRITE, &netcdf_id);
	handle_error(status, (std::string("Cannot open NetCDF file: ") + fname).c_str());
	return netcdf_id;
}

/// Reads in all the WATCH coordinates into a vector.
/** The index of a coordinate in the landpoints vector is the grid cell's id,
 *  which is used to find the NetCDF file for the grid cell.
 */
void load_gridlist(const char* dir_name, std::vector<landpoint>& landpoints) {
	using std::vector;
	using std::string;

	string gridlist_file = string(dir_name) + "/WFD-land-lat-long-z.nc";
	int ncid = open_ncdf(gridlist_file.c_str());

	// Get the land dimension and figure out how many grid cells there are

	int landid;
	int status = nc_inq_dimid(ncid, "land", &landid);
	handle_error(status, "land dimension");

	size_t num_gridcells;
	status = nc_inq_dimlen(ncid, landid, &num_gridcells);
	handle_error(status, "land dimension length");
	
	// Read in all longitudes and latitudes into two vectors

	vector<double> lons(num_gridcells), lats(num_gridcells);

	int lonid, latid;
	status = nc_inq_varid(ncid, "Longitude", &lonid);
	handle_error(status, "Longitude variable id");

	status = nc_inq_varid(ncid, "Latitude", &latid);
	handle_error(status, "Latitude variable id");
	
	nc_get_var_double(ncid, lonid, &lons.front());
	nc_get_var_double(ncid, latid, &lats.front());

	// Transfer the coordinates to the landpoints vector

	landpoints.resize(num_gridcells);

	for (size_t i = 0; i < landpoints.size(); i++) {
		landpoints[i].first = lons[i];
		landpoints[i].second = lats[i];
	}

	nc_close(ncid);
}

/// Returns the grid cell id for a coordinate
int get_cell_id(double lon, double lat, const std::vector<landpoint>& landpoints) {
	double c_lon = lon + .25;	// offset between a centre of the cell (WATCH) and
	double c_lat = lat + .25;	// lower-left corner (LPJ-GUESS)
	for (size_t i = 0; i < landpoints.size(); i++) {
		if (landpoints[i].first == c_lon && landpoints[i].second == c_lat) {
			return i;
		}
	}
	fail("Cell (%g,%g) wasn't found in the NetCDF grid\n", lon, lat);
	return -1;
}

/// Read in all data for a single variable
/** Data for leap years is skipped.
 *
 * \param dir_name Path to WATCH directory with all NetCDF files
 * \param var_name Which variable to read in
 * \param cell_id  id number for the grid cell to read in
 * \param diurnal  Whether the variable has subdaily or daily data */
void load_watch_data(const char* dir_name,
                     const char* var_name,
                     int cell_id, 
                     std::vector<double>& data,
                     bool diurnal) {
	using std::vector;
	
	// Open the NetCDF file for this variable and grid cell

	xtring fname;
	fname.printf("%s/%s/gc_%i.nc", dir_name, var_name, cell_id);

	int ncid = open_ncdf(fname);

	// Figure out the number of time steps

	int tstep_id;
	int status = nc_inq_dimid(ncid, "tstep", &tstep_id);
	handle_error(status, "tstep dimension");

	size_t num_timesteps;
	status = nc_inq_dimlen(ncid, tstep_id, &num_timesteps);
	handle_error(status, "tstep dimension length");

	// Read in all the data from the variable

	int var_id;
	status = nc_inq_varid(ncid, var_name, &var_id);
	handle_error(status, var_name);

	vector<double> raw_data(num_timesteps);
	status = nc_get_var_double(ncid, var_id, &raw_data.front());
	handle_error(status, var_name);

	// Transfer data from 'raw_data' to 'data', skipping leap days

	data.clear();
	data.reserve(num_timesteps);

	vector<double>::const_iterator raw_data_itr = raw_data.begin();

	const size_t steps_per_day = diurnal ? SUBDAILY : 1;

	int calendar_year = FIRST_WATCH_YEAR;

	// loop through the whole raw_data vector
	while (raw_data_itr != raw_data.end()) {

		// copy one year in each iteration
		for (int i = 0; i < 365; ++i) {
			if (i == 59 && Date::is_leap(calendar_year)) {
				// leap day, skip it
				raw_data_itr += steps_per_day;
			}
			data.insert(data.end(), raw_data_itr, raw_data_itr + steps_per_day);
			raw_data_itr += steps_per_day;
		}
		++calendar_year;
	}

	// Make sure we agree on leap days etc.
	assert(raw_data_itr == raw_data.end());

	nc_close(ncid);
}

#endif // LPJ_GUESS_WATCH_NETCDF_H
