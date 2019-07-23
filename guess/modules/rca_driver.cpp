///////////////////////////////////////////////////////////////////////////////////////
/// \file rca_driver.cpp
/// \brief Misc functions for the RCA coupling
///
/// \author Joe Lindström
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "rca_driver.h"
#include "guess.h"
#include "globalco2file.h"
#include "guessio.h"
#include <map>

xtring couple_src;
void transfer_drivers(Stand& stand, int& warnct, double lon_rca, double lat_rca, double temp, double netswrad,
                      double soilw_surf, double soilw_deep,
                      double rain, double snow) {

	// warnct         = warning counter (incremented for each warning, 
	//                  stops printing warnings when reaching a limit)
	// lon_rca        = longitude (for warning texts)
	// lat_rca        = latitude (for warning texts)
	// temp           = instantaneous temperate 2 m above ground (deg C)
	// netswrad       = net downward shortwave radiation (J/m2/s)
	// soilw_surf     = soil water content AWC (0-1) - surface layer 0-50 cm
	// soilw_deep     = soil water content AWC (0-1) - deep layer 50-150 cm
	// rain           = rain (kg/(m2 s) = mm/s)
	// snow           = snow (kg/(m2 s)

	const int MAXWARNINGS=100;

	Climate& climate=stand.climate;

	if (couple_src.lower() == "rca_sp")
		//PAR is stored in RCAData.bin
		climate.instype = PAR;
	else
		climate.instype = NETSWRAD_TS;

	// Initialize sum variable before accumulation
	// We assume Date.next_timestep() is called before calling this function
	if (date.timestep == 1){
		climate.temp_sum = climate.insol_sum = climate.prec_sum = 0.0;
		stand.wcont_sum[0] = stand.wcont_sum[1] = 0.0;
		climate.temp_min = 1000;
		climate.temp_max = -1000;
	}
		
	// Validate climate data

	if (temp<-100.0 || temp>100.0) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid temperature %gK received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0\n",
				temp+K2degC,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
		temp=0.0;
	}
	
	// TODO: changes the name netswrad to insol to make more sense
	if ((netswrad<0.0 || netswrad>1366.0) && (climate.instype == NETSWRAD_TS)) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid netswrad %g W/m2 received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0\n",
				netswrad,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
		netswrad=0.0;
	}

	if (soilw_surf<0.0) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid surface layer water content %g received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0.01 AWC\n",
				soilw_surf,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
		soilw_surf=0.01;
	}
 	else if (soilw_surf>1.0) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid surface layer water content %g received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 1.0 AWC\n",
				soilw_surf,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
		soilw_surf=1.0;
	}

	if (soilw_deep<0.0) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid deep layer water content %g received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0.01 AWC\n",
				soilw_deep,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
                soilw_deep=0.01;
	}
	else if (soilw_deep>1.0) {
		if (warnct++<MAXWARNINGS)
			dprintf("WARNING: invalid deep layer water content %g received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 1.0 AWC\n",
					soilw_deep,lon_rca,lat_rca,date.year,date.month,date.dayofmonth,date.timestep);
		soilw_deep=1.0;
	}

	if (rain < 0 || snow < 0) {
		 if (warnct++<MAXWARNINGS)
			  dprintf("WARNING: invalid rain/snow (%g/%g) received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0.0\n",
						 rain, snow, lon_rca, lat_rca, date.year, date.month, date.dayofmonth, date.timestep);
		 
		 rain = max(0.0, rain);
		 snow = max(0.0, snow);
	}

	if (rain > 10 || snow > 10) {
		 if (warnct++<MAXWARNINGS)
			  dprintf("WARNING: invalid rain/snow (%g/%g) received from RCA for stand (%0.3f,%0.3f) on %04d-%02d-%02d timestep %d - setting to 0.0\n",
						 rain, snow, lon_rca, lat_rca, date.year, date.month, date.dayofmonth, date.timestep);
		 
		 rain = rain > 10 ? 0 : rain;
		 snow = snow > 10 ? 0 : snow;
	}

	// Update cumulative sums
	climate.temp_sum+=temp;
	climate.insol_sum+=netswrad;
	// Simply add rain and snow, soilwater module will split
	// according to air temperature.
	climate.prec_sum += (rain + snow);

	climate.temp_max = max(temp, climate.temp_max);
	climate.temp_min = min(temp, climate.temp_min);

	climate.co2 = global_co2[date.year+calendar_year_offset];

	stand.wcont_sum[0] += soilw_surf;
	stand.wcont_sum[1] += soilw_deep;
}


void finalize_drivers(Stand& stand) {
	 Climate& climate = stand.climate;

	// Convert cumulative sums to daily average and transfer to Climate object
	// (unless still doing an RCA spinup - Ben 2007-12-19)

	if (!during_rca_spinup(date)) {
		climate.temp=climate.temp_sum/(double)(date.timestep+1);
		climate.insol=climate.insol_sum/(double)(date.timestep+1);
		climate.prec = (climate.prec_sum/(double)(date.timestep+1))*3600*24;
		climate.dtr = climate.temp_max - climate.temp_min;

		if (climate.dtr < 0) {
			fail("Failed to calculate dtr, no forcing data from RCA before call to finalize_drivers?");
		}

		if (ifprescribedwcont) {
			stand.wcont_prescribed[0] = stand.wcont_sum[0]/(double)(date.timestep+1);
			stand.wcont_prescribed[1] = stand.wcont_sum[1]/(double)(date.timestep+1);
		}
	}

	// Ndep
	Lamarque::NDepData& ndep = getNDepData(stand.gridcell.lon_cru, stand.gridcell.lat_cru);
	double mndrydep[12];
	double mnwetdep[12];
	ndep.get_one_calendar_year(date.year + calendar_year_offset, mndrydep, mnwetdep);

	climate.dndep = mndrydep[date.month] + mnwetdep[date.month];
	climate.dnfert = 0.0;

	// Reset counters for next day
	climate.temp_sum=climate.insol_sum=climate.prec_sum=0.0;
	stand.wcont_sum[0]=stand.wcont_sum[1]=0.0;
	climate.temp_min = 1000;
	climate.temp_max = -1000;
}


void compute_tile_lai(Stand& stand) {

	// Calculates various LAI terms required by RCA

	stand.laiphen_grass[date.day]=0.0;
	stand.laiphen_conifer[date.day]=0.0;
	stand.laiphen_broadleaf[date.day]=0.0;
	stand.laimax_conifer=0.0;
	stand.laimax_broadleaf=0.0;

	for (int p = 0; p < stand.nobj; p++) {
		Patch& patch=stand[p];
		Vegetation& vegetation=patch.vegetation;

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			if (indiv.id>=0 && indiv.alive) {

				// Only count individuals that have lived for
				// a year or more (Ben 2007-11-28)

				if (indiv.pft.lifeform==GRASS)
					stand.laiphen_grass[date.day]+=indiv.lai*indiv.phen;
				else if (indiv.pft.leafphysiognomy==NEEDLELEAF) {
					stand.laiphen_conifer[date.day]+=indiv.lai*indiv.phen;
					stand.laimax_conifer+=indiv.lai;
				}
				else {
					stand.laiphen_broadleaf[date.day]+=indiv.lai*indiv.phen;
					stand.laimax_broadleaf+=indiv.lai;
				}
			}

			vegetation.nextobj();
		}
	}

	stand.laiphen_grass[date.day]/=(float)stand.nobj;
	stand.laiphen_conifer[date.day]/=(float)stand.nobj;
	stand.laiphen_broadleaf[date.day]/=(float)stand.nobj;
	stand.laimax_conifer/=(float)stand.nobj;
	stand.laimax_broadleaf/=(float)stand.nobj;
}

namespace {
std::string ndep_path;
Lamarque::timeseriestype ndep_timeseries;
}

void init_RCA_NDep(const char* path, Lamarque::timeseriestype timeseries) {
	ndep_path = path;
	ndep_timeseries = timeseries;
}

Lamarque::NDepData& getNDepData(double lon_cru, double lat_cru) {

	using Lamarque::NDepData;
	using std::make_pair;

	static std::map<std::pair<double, double>, NDepData*> cache;

	std::map<std::pair<double, double>, NDepData*>::iterator itr = cache.find(make_pair(lon_cru, lat_cru));

	if (itr != cache.end()) {
		return *itr->second;
	}
	else {
		NDepData* ndep = new NDepData();
		ndep->getndep(ndep_path.c_str(), lon_cru, lat_cru, ndep_timeseries);
		cache.insert(make_pair(make_pair(lon_cru, lat_cru), ndep));
		return *ndep;
	}
}
