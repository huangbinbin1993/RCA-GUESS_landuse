#!/bin/bash
describe_benchmark "LPJ-GUESS - Global EMDI Benchmarks (using global PFTs)"

common1961to1990.sh
common1961to1990gmapall.sh -portrait -pixsize 5 5

emdi_npp_compare.sh gridlist.txt anpp1961to1990.txt anpp_comparison.png "Global NPP"
describe_image anpp_comparison.png "Modelled annual NPP (1961-90 average) versus EMDI annual NPP observations. /
EMDI Global Class A sites only. Units: kgC m-2." embed