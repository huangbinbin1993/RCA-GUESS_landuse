#!/bin/bash
describe_benchmark "LPJ-GUESS - Pristine European Site Benchmarks"

common1961to1990.sh
common1961to1990gmapall.sh

dominance dens1961to1990.txt dens1961to1990max.txt
describe_textfile dens1961to1990max.txt "Dominant PFT/Species (according to 1961-90 average tree density). Units: trees m-2"
