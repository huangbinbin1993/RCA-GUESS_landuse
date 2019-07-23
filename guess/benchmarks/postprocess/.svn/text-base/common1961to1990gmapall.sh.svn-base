#!/bin/bash
#
# Runs gmapall (running gmap on all species/PFTs) generating maps
# for LAI and CMASS for 1961 to 1990 (see also common1961to1990.sh)
#
# All parameters to this script is passed along to gmap.

# Render the maps in scalar mode, giving each file an lai or
# cmass prefix.
gmapall lai1961to1990.txt -P lai_ -s $@
describe_images "LAI For All PFTs/Species (1961-90 average). Units: m2 m-2"  lai_*.jpg

gmapall cmass1961to1990.txt -P cmass_ -s $@
describe_images "CMASS For All PFTs/Species (1961-90 average). Units: kgC m-2" cmass_*.jpg

gmapall cton_leaf1961to1990.txt -P cton_leaf_ -legend common/legend_cton.txt
describe_images "Leaf C:N Ratio For All PFTs/Species (1961-90 average). Units: kgC kgN-1" cton_leaf_*.jpg
