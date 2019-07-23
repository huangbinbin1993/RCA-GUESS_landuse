///////////////////////////////////////////////////////////////////////////////////////
/// \file landcoverstrategies.cpp
/// \brief Implementations of the various land cover strategies
///
/// \author Joe Lindström
///
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "landcoverstrategies.h"
#include "guess.h"

std::auto_ptr<LandCoverStrategy> landcoverstrategy(0);

/// Trivial LandCoverInfo implementation for PNV
/** 
 *  The only thing we store is the original land cover
 *  fraction received from RCA.
 */
class PNVCoverInfo : public LandCoverInfo {
public:
  PNVCoverInfo(double rca_fgrid_land) 
    : fgrid_land(rca_fgrid_land) {
  }

  double get_fgrid_land() const {
    return fgrid_land;
  }

  double get_fland_crop(int year) const {
    return 0;
  }

  double get_fland_bare(int year) const {
    return 0;
  }

private:
  double fgrid_land;
};

/// Simple LandCoverInfo implementation for static fractions
/**
 *  Here we simply store a double for each fraction since they
 *  remain the same throughout the simulation.
 *
 *  Used by the EcoclimapStrategy.
 */
class FixedCoverInfo : public LandCoverInfo {
public:
  FixedCoverInfo(double fgrid_land, double fland_crop, double fland_bare) {
    this->fgrid_land = fgrid_land;
    this->fland_crop = fland_crop;
    this->fland_bare = fland_bare;
  }

  double get_fgrid_land() const {
    return fgrid_land;
  }

  double get_fland_crop(int year) const {
    return fland_crop;
  }

  double get_fland_bare(int year) const {
    return fland_bare;
  }

private:
  double fgrid_land;
  double fland_crop;
  double fland_bare;
};


////////////////////////////////////////////////////////////////////////////////
// Implementation of PNVStrategy
////////////////////////////////////////////////////////////////////////////////

PNVStrategy::PNVStrategy() {
}

LandCoverInfo* PNVStrategy::get_landcoverinfo(double lon, 
					      double lat, 
					      double rca_fgrid_land) {
  if (rca_fgrid_land > 0)
    return new PNVCoverInfo(rca_fgrid_land);
  else
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of EcoclimapStrategy
////////////////////////////////////////////////////////////////////////////////

EcoclimapStrategy::EcoclimapStrategy(const char* path) {
  file = fopen(path, "r");
	 
  if (!file) {
    fail("EcoclimapStrategy: could not open %s for input", path);
  }
}

LandCoverInfo* EcoclimapStrategy::get_landcoverinfo(double lon,
						    double lat,
						    double rca_fgrid_land) {
  // Classification in ECOCLIMAP:
  // tile 1: bare soil
  // tile 2: rocks
  // tile 3: permanent snow
  // tile 4: deciduous broadleaf trees
  // tile 5: coniferous trees
  // tile 6: evergreen broadleaf trees
  // tile 7: C3 crops
  // tile 8: C4 crops
  // tile 9:  IRRIGATED crops
  // tile 10: natural herbaceous NOT irrigated (temperate and polar)
  // tile 11: natural herbaceous NOT irrigated (tropical and equatorial)
  // tile 12: swamp herbaceous and gardens

  const int NCLASS=12;
  double EPS=0.01; // tolerance (degrees) for pixel lon/lat
  double dlon,dlat,frland,frlc[NCLASS];
  bool found=false;
  double frac;

  rewind(file);
  readfor(file,""); // past header row
	
  while (!feof(file) && !found) {
    readfor(file,"f,f,f,12f",&dlon,&dlat,&frland,frlc);
    if (!feof(file)) {
      if (fabs(dlon-lon)<EPS && fabs(dlat-lat)<EPS)
	found=true;
    }
  }

	
  if (!found) {
    return 0;
  }

  // else

  double fgrid_land, fland_crop, fland_bare;
  fgrid_land=frland;
  fland_crop=frlc[6]+frlc[7]+frlc[8]+frlc[11]; // C3 crops, C4 crops, irrigated crops, wetlands
  fland_bare=frlc[1]+frlc[2]; // rocks, permanent snow
	
  if (frland && fland_crop+fland_bare>1.0) {
	
    // Correct minor (rounding) errors
    if (fland_crop+fland_bare<1.1) {
      frac=1.0/(fland_crop+fland_bare);
      fland_crop*=frac;
      fland_bare*=frac;
    }
    else {
      fail("EcoclimapStrategy::get_landcoverinfo: fcrops (%g) + fbare (%g) = %g (should be max 1.0)\n",
	   fland_crop,fland_bare,fland_crop+fland_bare);
    }
  }
	
  if (!frland) {
    return 0;
  }
	
  // else	
  return new FixedCoverInfo(fgrid_land, fland_crop, fland_bare);
}


////////////////////////////////////////////////////////////////////////////////
// Implementation of YearlyStrategy
////////////////////////////////////////////////////////////////////////////////

YearlyStrategy::YearlyCoverInfo::YearlyCoverInfo(
	 double rca_fgrid_land,
	 const std::vector<Record>& records) {

	 fgrid_land = rca_fgrid_land;
	 this->records = records;
}

const YearlyStrategy::Record& 
YearlyStrategy::YearlyCoverInfo::get_record_for_year(int year) const {
	 if (year <= records.front().year) {
		  return records.front();
	 }
	 else if (year >= records.back().year) {
		  return records.back();
	 }
	 else {
		  return records[year-records.front().year];
	 }	 
}

double YearlyStrategy::YearlyCoverInfo::get_fland_crop(int year) const {
	 return get_record_for_year(year).open_veg;
}

double YearlyStrategy::YearlyCoverInfo::get_fland_bare(int year) const {
	 return get_record_for_year(year).non_veg;
}

YearlyStrategy::YearlyStrategy(const char* path) {

	 // opens the file and reads in all records to memory

	 FILE* file = fopen(path, "r");

	 if (!file) {
		  fail("YearlyStrategy: could not open %s for input", path);
	 }

	 readfor(file,""); // past header row
	
	 while (!feof(file)) {
		  double dlon, dlat, non_veg, open_veg, for_veg;
		  int year;
		  readfor(file,"f,f,i,f,f,f",
					 &dlon, &dlat, &year, &non_veg, &open_veg, &for_veg);
		  if (!feof(file)) {

				// If all fractions are set to zero, that means non-land,
				// so we simply skip the record (coordinates for which we
				// have no records are considered to be non-land)
				if (non_veg == 0.0 && open_veg == 0.0 && for_veg == 0.0) {
					 continue;
				}

				// Otherwise the fractions should add up to 1
				if (fabs( 1 - ( non_veg + open_veg + for_veg ) ) > 0.1) {
					 fail("YearlyStrategy: fractions don't add to 1, lon=%.2f, lat=%.2f, year=%d",
							dlon, dlat, year);
				}

				Record record;
				record.lon = dlon;
				record.lat = dlat;
				record.year = year;
				record.non_veg = non_veg;
				record.open_veg = open_veg;
				records.push_back(record);
		  }
	 }
	
	 fclose(file);
}

LandCoverInfo* YearlyStrategy::get_landcoverinfo(
								double lon, // CRU coordiantes
								double lat, // CRU coordiantes
								double rca_fgrid_land) {

	double lon_c = lon + 0.25; // CRU coordiantes at the centre of the grid
	double lat_c = lat + 0.25; // CRU coordiantes at the centre of the grid
    bool found = false;

	// Only if both CRU and RCA regard this land should go ahead
	if (rca_fgrid_land <= 0)
		return 0;

	 // get the records for the coordinate
	 std::vector<Record> gridcell_records;

	// First attempt
	double EPS=0.01; // tolerance (degrees) for pixel lon/lat
	for (size_t i = 0; i < records.size(); i++) {

		if ( fabs(lon_c - records[i].lon) < EPS &&
			 fabs(lat_c - records[i].lat) < EPS) {

			found = true;

			// to improve efficency, break the loop when step to the next coordinate after record is found
			if ( found && gridcell_records.size() >= 1 &&
				!negligible(gridcell_records.front().lon - records[i].lon) &&
				!negligible(gridcell_records.front().lat - records[i].lat) ) {
				break;
			}

			/* 
			// make sure the records are sorted and that there are no gaps
			if (gridcell_records.size() >= 1 &&
				records[i].year != gridcell_records.back().year+1) {
				fail("YearlyStrategy: Not all years in land cover file are consecutive!, missing year: %d, %d,%d", gridcell_records.size(), records[i].year, gridcell_records.back().year + 1);
			} */
				
			// Fast Archive may help improve efficiency -mcw
			gridcell_records.push_back(records[i]);
		}
	}

	if (!gridcell_records.empty()) {
		dprintf("landcover record for this gridcell: %d \n", gridcell_records.size());
		return new YearlyCoverInfo(rca_fgrid_land, gridcell_records);
	}
	else {
		// Second attemp: increase the tolerance degree
		EPS = 0.5;
		found = false;
	 	for (size_t i = 0; i < records.size(); i++) {

			if ( fabs(lon_c - records[i].lon) < EPS &&
				 fabs(lat_c - records[i].lat) < EPS) {

				found = true;

				// why still more than 106 records found? -mcw
	 	     	// to improve efficency, break the loop when step to the next coordinate after record is found
	 	     	if ( found && gridcell_records.size() >= 1 &&
	 	   	  		 !negligible(gridcell_records.front().lon - records[i].lon) &&
					 !negligible(gridcell_records.front().lat - records[i].lat) ) {
	 	   	 	 	break;
			 	}

				/*
				// make sure the records are sorted and that there are no gaps
				if (gridcell_records.size() >= 1 &&
					records[i].year != gridcell_records.back().year+1) {
					fail("YearlyStrategy: Not all years in land cover file are consecutive!, missing year: %d, %d,%d", gridcell_records.size(), records[i].year, gridcell_records.back().year + 1);
				} */ 

	 	   		// Fast Archive may help improve efficiency -mcw
				gridcell_records.push_back(records[i]);
			}
	 	}

	 	if (!gridcell_records.empty()) {
	 	     dprintf("landcover record for this gridcell (tolerance inceased): %d \n", gridcell_records.size());
	 	     return new YearlyCoverInfo(rca_fgrid_land, gridcell_records);
	 	}
	 	else {
	 	   	// Otherwise: regard this as natural vegetation
			// dprintf causes seg. fault, do not use here! 
	 	   	// dprintf("Cannot find landcover record for this gridcell (%6.4f, %6.4f), 100% natural veg is used\n", lon, lat);

			Record rec0;
			rec0.lon = lon_c;
	 	   	rec0.lat = lat_c;
	 	   	rec0.year = 1901; // only one year is needed, get_record_for_year() will hand for other years
	 	   	rec0.non_veg = 0.0;
	 	   	rec0.open_veg = 0.0;
	 	   	gridcell_records.push_back(rec0);

	 	   	return new YearlyCoverInfo(rca_fgrid_land, gridcell_records);
	 	}
	}
}
