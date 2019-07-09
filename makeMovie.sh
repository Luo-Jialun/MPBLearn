#!/bin/bash

resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn

ezFile=brd_wvg_cvt_diff_spc-with_cavity-3_r-0.370_NRow-5_sep-1.990_excite_fc-0.064_bw-0.030_flux_fc-0.400_df-0.600_ez.h5
epsFile=brd_wvg_cvt_diff_spc-with_cavity-3_r-0.370_NRow-5_sep-1.990_excite_fc-0.064_bw-0.030_flux_fc-0.400_df-0.600_eps.h5

python PostProcessingUtils.py --make-movie $ezFile $epsFile



# ezFile=brd_wvg_cvt_diff_spc-with_cavity-5_sep-2.000_excite_fc-0.400_bw-0.300_flux_fc-0.400_df-0.600_ez.h5
# epsFile=brd_wvg_cvt_diff_spc-with_cavity-5_sep-2.000_excite_fc-0.400_bw-0.300_flux_fc-0.400_df-0.600_eps.h5

# python PostProcessingUtils.py --make-movie $ezFile $epsFile