///////////////////////////////////////////////////////////////////////////////////////
/// \file date_test.cpp
/// \brief Unit tests for the Date class
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "guess.h"

namespace {

/// Help function that simply calls Date::next n times
void take_n_steps(Date& d, int n) {
	 for (int i = 0; i < n; i++) {
		  d.next();
	 }
}

}


TEST_CASE("date/construction", "Some basic tests of constructing dates") {
	 Date d;
	 d.init(1);

	 REQUIRE(d.day == 0);
	 REQUIRE(d.dayofmonth == 0);
	 REQUIRE(d.month == 0);
	 REQUIRE(d.year == 0);
	 REQUIRE(d.islastyear);
	 REQUIRE(!d.islastmonth);
	 REQUIRE(!d.islastday);
	 REQUIRE(!d.ismidday);
}


TEST_CASE("date/stepping", "Tests Date::next()") {
	 Date d;
	 d.init(1);

	 // Given that we take these number of steps...
	 int steps[] = { 1, 30, 28, 306 };
	 // ...we expect to end up at these combinations of year, month and day
	 int expectations[][3] = { {0, 0, 1}, { 0, 1, 0 }, {0, 2, 0}, {1, 0, 0} };

	 for (int i = 0; i < sizeof(steps)/sizeof(steps[0]); i++) {

		  take_n_steps(d, steps[i]);

		  REQUIRE(d.year == expectations[i][0]);
		  REQUIRE(d.month == expectations[i][1]);
		  REQUIRE(d.dayofmonth == expectations[i][2]);
	 }
}


TEST_CASE("date/months", "Tests prevmonth and nextmonth") {
	 Date d;
	 d.init(1);

	 // We start in January
	 REQUIRE(d.nextmonth() == 1);
	 REQUIRE(d.prevmonth() == 11);

	 // Go to February
	 take_n_steps(d, 31);
	 
	 REQUIRE(d.nextmonth() == 2);
	 REQUIRE(d.prevmonth() == 0);

	 // Go to December
	 take_n_steps(d, 330);

	 REQUIRE(d.nextmonth() == 0);
	 REQUIRE(d.prevmonth() == 10);

	 // Go to January the following year
	 take_n_steps(d, 10);

	 REQUIRE(d.nextmonth() == 1);
	 REQUIRE(d.prevmonth() == 11);
}
