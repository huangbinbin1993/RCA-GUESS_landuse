///////////////////////////////////////////////////////////////////////////////////////
/// \file landcoverstrategies.h
/// \brief Classes for handling the different land cover schemes we support.
///
/// We use the abstract factory design pattern here to give a flexible system for 
/// dealing with different land cover schemes. There's a global LandCoverStrategy
/// (which can be for instance a PNV strategy if we simulate Potential Natural 
/// Vegetation, or an ECOCLIMAP strategy). This global object is used to create
/// objects that keep track of the land cover fractions for each gridcell.
///
/// To support new land cover schemes in the future, implement a new LandCoverStrategy
/// sub-class.
///
/// The I/O module is responsible for setting the current LandCoverStrategy, based
/// on what's chosen in the ins file.
///
/// \author Joe Lindström
///
/// $Date: 2012-01-18 11:23:57 +0100 (Wed, 18 Jan 2012) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVERSTRATEGIES_H
#define LPJ_GUESS_LANDCOVERSTRATEGIES_H

#include <memory>
#include <stdio.h>
#include <vector>

/// Abstract base class for keeping track of land cover fractions for a single grid cell
/** 
 *  This class defines the interface for classes that, in different ways, handle
 *  the land cover fractions for a single grid cell.
 *
 *  The interface makes available the following variables:
 *
 *  fgrid_land  - Fraction of this grid cell that is land (static)
 *  fland_crop  - Fraction of LAND that is non-forest vegetation due to land use
 *  fland_bare  - Fraction of LAND that is non-vegetated
 *
 *  The latter two can vary with the calendar year.
 *
 *  This makes it possible to implement various schemes, from simple PNV
 *  (Potential Natural Vegetation), where fland_crop and fland_bare are always
 *  0, and fgrid_land is the value we get from RCA, to more sophisticated
 *  implementations where the fractions are read from file and change each year.
 *
 *  Objects implementing this interface should not be created directly, but rather
 *  by calling get_landcoverinfo() for the currently enabled LandCoverStrategy
 *  (see below).
 */
class LandCoverInfo {
 public:
  /// Fraction of this grid cell that is land
  virtual double get_fgrid_land() const = 0;

  /// Fraction of LAND on this grid cell that is non-forest vegetation due to land use
  /** (crops, swamp herbaceous and gardens) 
   *  
   *  \param year  Calendar year for which to get the fraction
   */
  virtual double get_fland_crop(int year) const = 0;

  /// Fraction of LAND on this grid cell that is non-vegetated
  /** (rock, permanent snow or bare soil) 
   *  
   *  \param year  Calendar year for which to get the fraction
   */
  virtual double get_fland_bare(int year) const = 0;
};

/// A factory class for creating LandCoverInfo objects
/**
 *  Sub-classes of LandCoverStrategy can create LandCoverInfo object of a 
 *  certain type.
 *  
 *  For instance, the PNVStrategy class creates simple LandCoverInfo objects
 *  for which crop and bare is always 0, regardless of the gridcell's 
 *  coordinates.
 */
class LandCoverStrategy {
 public:
  /// Creates a LandCoverInfo object for a given coordinate
  /**
   *  \param lon       Longitude
   *  \param lat       Latitude
   *  \rca_fgrid_land  Fraction of gridcell that is land according to RCA,
   *                   this is used when we run with PNV, but ignored by
   *                   the EcoclimapStrategy which reads land fractions from
   *                   file instead.
   *
   *  \return          A LandCoverInfo object, or 0 if there is no land for
   *                   this gridcell.
   */
  virtual LandCoverInfo* get_landcoverinfo(double lon, 
					   double lat, 
					   double rca_fgrid_land) = 0;
};

/// Our global LandCoverStrategy
/**
 *  This should typically be created by the I/O module once
 *  we've parsed the ins file and know what kind of land cover
 *  scheme to use.
 */
extern std::auto_ptr<LandCoverStrategy> landcoverstrategy;


/// To be used when running with Potential Natural Vegetation
class PNVStrategy : public LandCoverStrategy {
 public:
  PNVStrategy();

  LandCoverInfo* get_landcoverinfo(double lon, 
				   double lat, 
				   double rca_fgrid_land);
};

/// To be used when running with an ECOCLIMAP text file
class EcoclimapStrategy : public LandCoverStrategy {
 public:
  EcoclimapStrategy(const char* path);

  LandCoverInfo* get_landcoverinfo(double lon, 
				   double lat, 
				   double rca_fgrid_land);
	 
 private:
  /// The opened text file
  FILE* file;
};

/// To be used when using a text file with one line per coordinate and year
class YearlyStrategy : public LandCoverStrategy {
public:
	 YearlyStrategy(const char* path);

	 LandCoverInfo* get_landcoverinfo(double lon, double lat, double rca_fgrid_land);

private:

	 struct Record {
		  double lon;
		  double lat;
		  int year;
		  double non_veg;
		  double open_veg;
	 };

	 class YearlyCoverInfo : public LandCoverInfo {
	 public:
		  YearlyCoverInfo(double rca_fgrid_land, const std::vector<Record>& records);
		  
		  double get_fgrid_land() const { return fgrid_land; }
		  double get_fland_crop(int year) const;
		  double get_fland_bare(int year) const;

	 private:

		  /// Figures out which record to use for a given year
		  /**
			* Years before the first record gets the first record,
			* years after the last record gets the last.
			*/
		  const Record& get_record_for_year(int year) const;

		  double fgrid_land;

		  /// All records for this coordinate
		  std::vector<Record> records;
	 };

	 /// All the records from the file, cached
	 std::vector<Record> records;
};

#endif // LPJ_GUESS_LANDCOVERSTRATEGIES_H
