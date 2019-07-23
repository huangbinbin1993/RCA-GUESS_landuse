///////////////////////////////////////////////////////////////////////////////////////
/// \file math_test.cpp
/// \brief Unit tests functionality in guessmath.h
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "guessmath.h"

TEST_CASE("Historic/add", "Some basic tests of adding values to a Historic") {
	Historic<double, 3> history;

	REQUIRE(history.size() == 0);

	history.add(1);

	REQUIRE(history.size() == 1);
	REQUIRE(history.sum() == Approx(1));
	REQUIRE(history[0] == 1);
	REQUIRE(history.mean() == Approx(1));

	history.add(2);

	REQUIRE(history.size() == 2);
	REQUIRE(history.sum() == Approx(3));
	REQUIRE(history[0] == 1);
	REQUIRE(history[1] == 2);
	REQUIRE(history.mean() == Approx(1.5));

	history.add(3);

	REQUIRE(history.size() == 3);
	REQUIRE(history.sum() == Approx(6));
	REQUIRE(history[0] == 1);
	REQUIRE(history[1] == 2);
	REQUIRE(history[2] == 3);
	REQUIRE(history.mean() == Approx(2));

	history.add(4);

	REQUIRE(history.size() == 3);
	REQUIRE(history.sum() == Approx(9));
	REQUIRE(history[0] == 2);
	REQUIRE(history[1] == 3);
	REQUIRE(history[2] == 4);
	REQUIRE(history.mean() == Approx(3));
}
