#!/bin/bash

# This script produces a scatter plot for comparison of monthly values.
# It is meant to be used with the monthly output from euroflux
# benchmarks, so the following applies:
#
# If a value is -9999.0, the data point is ignored.


if [ $# -lt 4 ]; then
    echo "Usage: $0 <file_observed> <file_modeled> <filename>.png <title>"
    echo
    echo "For instance:"
    echo "$0 eurofluxmonthly_nee.out mnee.out euroflux_nee.png \"NEE\""
    exit 1
fi

source `dirname $0`"/scatter_plot.sh"

# Get the arguments
DATA_FILE_OBS=$1
DATA_FILE_MOD=$2
OUT_FILE=$3
TITLE=$4

# Create temporary files for only the data to plot
OBS_DATA=$(mktemp)
MOD_DATA=$(mktemp)
TEMP_FILTERED=$(mktemp)

# awk script for getting only the values from a monthly output file, one per line
GET_ONLY_DATA='NR>1 { for (i=4; i <= NF; i++) print $i }'

# Get the data points into two separate files
awk "$GET_ONLY_DATA" $DATA_FILE_OBS > $OBS_DATA
awk "$GET_ONLY_DATA" $DATA_FILE_MOD > $MOD_DATA

# Combine into one file with two columns for gnuplot, remove missing data
paste $OBS_DATA $MOD_DATA | grep -v "\-9999.000" > $TEMP_FILTERED

XLABEL=$TITLE"_obs"
YLABEL=$TITLE"_mod"

# Call the scatter_plot function to do the actual plotting
scatter_plot "${TITLE}" "${XLABEL}" "${YLABEL}" ${TEMP_FILTERED} ${OUT_FILE}

# Clean up temp files
rm $TEMP_FILTERED
rm $OBS_DATA
rm $MOD_DATA
