///////////////////////////////////////////////////////////////////////////////////////
/// \file climate_test.cpp
/// \brief Unit tests for functions processing climate data
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "driver.h"
#include <algorithm>

namespace {

/// Convenience function for testing prdaily
/** Tests prdaily for given monthly conditions, the same conditions
  * are used throughout the year, makes it easy to write test cases
  * without specifying 24 values.
  */
bool verify_prdaily_single_month(double prec, double wetdays) {

	double monthly_prec[12];
	std::fill(monthly_prec, monthly_prec+12, prec);

	double monthly_wetdays[12];
	std::fill(monthly_wetdays, monthly_wetdays+12, wetdays);

	double days[365];

	prdaily(monthly_prec, days, monthly_wetdays, 12345678);

	// Verify monthly sums and number of wet days
	const double SUM_TOLERANCE = 0.1;

	Date date;
	int current_day = 0;
	for (int m = 0; m < 12; m++) {
		double sum = 0;
		int wetcount = 0;

		for (int d = 0; d < date.ndaymonth[m]; d++) {
			sum += days[current_day];
			if (days[current_day] > 0) {
				wetcount++;
			}
			current_day++;
		}

		// Check that the sum for this month isn't too far from
		// the prescribed monthly precipitation
		if (fabs(sum-prec) > SUM_TOLERANCE) {
			return false;
		}

		// Verify wetcount?
	}

	return true;
}
}

TEST_CASE("climate/prdaily", "Tests for the prdaily function") {
	// Test no water
	REQUIRE(verify_prdaily_single_month(0, 0));

	// Regression test: 
	// Very small number of wet days used to cause infinite loop
	REQUIRE(verify_prdaily_single_month(1, 0.001));

	// Regression test: 
	// Very little precipitation and many wet days used to cause infinite loop
	REQUIRE(verify_prdaily_single_month(0.1, 30));
}
