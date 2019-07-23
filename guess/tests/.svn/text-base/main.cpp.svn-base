///////////////////////////////////////////////////////////////////////////////////////
/// \file main.cpp
/// \brief Main function for the unit tests
///
/// \author Joe Lindström
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "shell.h"

int main(int argc, char** argv) {

	// Set a shell so we have working dprintf, fail etc.
	set_shell(new CommandLineShell("tests.log"));

	// Let CATCH do the rest
	int result = Catch::Main(argc, argv);

	return result;
}
