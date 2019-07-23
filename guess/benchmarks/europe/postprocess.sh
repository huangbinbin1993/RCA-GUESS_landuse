#!/bin/bash
describe_benchmark "LPJ-GUESS - European Benchmarks"

common1961to1990.sh
common1961to1990gmapall.sh

gmap lai1961to1990max.txt -t 'Dominant Species (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxLAIeurope.jpg
describe_image maxLAIeurope.jpg "Species With the Highest LAI in Each Gridcell (1961-90 average)"

gmap cmass1961to1990max.txt -t 'Dominant Species (greatest cmass)' -lon 1 -lat 2 -i 3 -legend legend_europe.txt -o maxCMASSeurope.jpg
describe_image maxCMASSeurope.jpg "Species With the Highest CMASS in Each Gridcell (1961-90 average)"

tslice aiso.out -o aiso1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3 
tslice amon.out -o amon1961to1990.txt -f 560 -t 589 -lon 1 -lat 2 -y 3 
gmap aiso1961to1990.txt -t 'Isoprene flux (mg C/m2/y)' -legend legend_bvoc_europe.txt -lon 1 -lat 2 -i 'Total' -o aiso_tot.jpg
gmap amon1961to1990.txt -t 'Monoterpene flux (mg C/m2/y)' -legend legend_bvoc_europe.txt -lon 1 -lat 2 -i 'Total' -o amon_tot.jpg
describe_image aiso_tot.jpg "Annual isoprene flux (1961-90 average)"
describe_image amon_tot.jpg "Annual monoterpene flux (1961-90 average)"

compute cpool.out -n -o cpool_total.out -i Lon Lat Year Total
compute cflux.out -n -o cflux_nee.out -i Lon Lat Year NEE
balance -pool cpool_total.out -flux cflux_nee.out -start 500 -end 605 -matter C
describe_textfile Cbalance_totalerror_GtC.txt "European Terrestrial Carbon Uptake, 1901 to 2006. /
Determined using C pools (pool_GtC), Cumulative C fluxes (flux_GtC), and their absolute difference (absdiff_GtC)"

compute npool.out -n -o npool_total.out -i Lon Lat Year Total
compute nflux.out -n -o nflux_nee.out -i Lon Lat Year 'nee_m2=NEE/10000'
balance -pool npool_total.out -flux nflux_nee.out -start 500 -end 605 -matter N
describe_textfile Nbalance_totalerror_GtN.txt "European Terrestrial Nitrogen Uptake, 1901 to 2006. /
Determined using N pools (pool_GtN), Cumulative N fluxes (flux_GtN), and their absolute difference (absdiff_GtN)"
