///////////////////////////////////////////////////////////////////////////////////////
/// \file dllmain.cpp
/// \brief Main module for interface to Windows shell
///
/// The shell should call function dll_main, passing a GuessParam object containing run
/// time argument(s) for GUESS and pointers to the executable's own callback functions.
///
/// \author Ben Smith
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_MAIN_H
#define LPJ_GUESS_MAIN_H

#include "guess.h"

///////////////////////////////////////////////////////////////////////////////////////
// GUESSPARAM
// One object of this class should be passed to function dll_main when called from 
// the shell. It contains the run time argument for guess (corresponding to the single
// command line argument when GUESS is run as a command line/console application), and
// pointers to callback functions in the executable file.

struct PlotArgs {
	
	// Data buffer object to store arguments to "plot" message

	xtring window_name;
	xtring series_name;
	double x;
	double y;
	bool rescale;
};

// Type definitions for callback functions

typedef void(MessagePrintString)(xtring* string);
typedef void(MessagePlot)(PlotArgs* plotargs);
typedef void(MessageResetwindow)(xtring* string);
typedef void(MessageFinished)();
typedef void(MessageClearGraphs)();

struct GuessParam {

	// ARGUMENTS TO BE PASSED ON TO GUESS

	int argc; // number of run time arguments in argv
	char** argv; // run time arguments (a string, normally the name of an ins file)

	// POINTERS TO DATA BUFFER OBJECTS MAINTAINED BY EXECUTABLE

	xtring* poutput;
	PlotArgs* pplotargs;
	
	// POINTERS TO CALLBACK FUNCTIONS

	MessagePrintString* message_print_string;
		// When called, sends a message to shell to print string in *poutput
	MessagePlot* message_plot;
		// When called, sends a message to shell to plot data using information
		// in *pplotargs
	MessageFinished* message_finished;
		// When called, sends a message to shell that the run is finished
	MessageClearGraphs* message_clear_graphs;
		// When called, sends a message to shell to clear data in graph windows
	MessageResetwindow* message_resetwindow;
		// When called, sends a message to shell to clear data in a particular graph
		// window
};

///////////////////////////////////////////////////////////////////////////////////////
// DECLARATIONS FOR THE EXPORTED FUNCTIONS

__declspec(dllexport) int dll_main(GuessParam param);
__declspec(dllexport) void cleanup_print_string(xtring* pxtring);
__declspec(dllexport) void cleanup_plot(PlotArgs* pplotargs);
__declspec(dllexport) void abort_run();

#endif // LPJ_GUESS_MAIN_H
