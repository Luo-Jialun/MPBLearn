#!/bin/bash

resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn
# cd $resultPath

epsSet=4.84
referenceFile=wvg_no-cavity_eps-4.84_sep-0_excite_fc-0.4_bw-0.4_flux_fc-0.4_df-0.6_flux.csv
testFiles=wvg_with_cavity-1_eps-4.84_sep-1.000_excite_fc-0.4_bw-0.4_flux_fc-0.4_df-0.6_flux.csv

# test2=wvg_with_cavity-5_sep-1.000_excite_fc-0.35_bw-0.9_flux_fc-0.435_df-0.6_flux.csv


python $scriptPath/PostProcessingUtils.py -p $referenceFile  $testFiles 

# for fluxDataFile in *.csv; do
#   echo $fluxDataFile
# done
