#!/bin/bash

resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn
# cd $resultPath

epsSet=4.84
# referenceFile=no-cavity_sep-0_excite_fc-0.400_bw-0.600_flux_fc-0.400_df-0.600_flux.csv
# testFiles=with_cavity-5_sep-2.000_excite_fc-0.400_bw-0.300_flux_fc-0.400_df-0.600_flux.csv

# # test2=wvg_with_cavity-5_sep-1.000_excite_fc-0.35_bw-0.9_flux_fc-0.435_df-0.6_flux.csv


# python $scriptPath/PostProcessingUtils.py -p $referenceFile  $testFiles 

# for fluxDataFile in *.csv; do
#   echo $fluxDataFile
# done



referenceFile=no-cavity_r-0.380_NRow-5_sep-0_excite_fc-0.400_bw-0.700_flux_fc-0.400_df-0.600_flux.csv

testFiles=with_cavity-*_r-0.380_NRow-5_sep-1.560_excite_fc-0.400_bw-0.700_flux_fc-0.400_df-0.600_flux.csv


python $scriptPath/PostProcessingUtils.py -p $referenceFile  $testFiles