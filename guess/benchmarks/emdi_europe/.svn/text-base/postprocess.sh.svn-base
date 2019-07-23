#!/bin/bash
describe_benchmark "LPJ-GUESS - European EMDI Benchmarks (using European species)"

common1961to1990.sh
common1961to1990gmapall.sh -pixsize 0.5 0.5

emdi_npp_compare.sh gridlist.txt anpp1961to1990.txt anpp_comparison.png "Europe NPP"
describe_image anpp_comparison.png "Modelled annual NPP (1961-90 average) versus EMDI annual NPP observations. /
EMDI European Class A and B sites combined. Units: kgC m-2." embed
