#!/bin/bash
pwd
resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn
cd $resultPath
pwd

referenceFile=wvg_with_no_cavity_exciationParam_fcen-0.4319_bw-0.05_fluxParam_fcen-0.435_df-0.1_flux.csv

for fluxDataFile in *with_cavity*exciationParam_fcen-0.4319_bw-0.05_fluxParam_fcen-0.435_df-0.1_flux.csv; do
  echo $fluxDataFile
done

python $scriptPath/PostProcessingUtils.py

# for fluxDataFile in *.csv; do
#   echo $fluxDataFile
# done