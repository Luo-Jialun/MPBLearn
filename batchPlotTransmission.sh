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



referenceFile=no-cavity_r-0.360_NRow-1_sep-0_excite_fc-0.250_bw-0.200_flux_fc-0.250_df-0.200_flux.csv

testFiles=with_cavity-6_r-0.360_NRow-1_sep-1.200_excite_fc-0.250_bw-0.200_flux_fc-0.250_df-0.200_flux.csv

# testFiles=with_cavity-4_r-0.370_sep-1.250_excite_fc-0.300_bw-0.400_flux_fc-0.400_df-0.600_flux.csv
# referenceFile=$testFiles

python $scriptPath/PostProcessingUtils.py -p $referenceFile  $testFiles