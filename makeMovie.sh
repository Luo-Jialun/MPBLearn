#!/bin/bash

resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn

ezFile=brd_wvg_cvt_diff_spc-with_cavity-3_r-0.380_NRow-5_sep-1.560_excite_fc-0.400_bw-0.700_flux_fc-0.400_df-0.600_field.h5
epsFile=brd_wvg_cvt_diff_spc-with_cavity-3_r-0.380_NRow-5_sep-1.560_excite_fc-0.400_bw-0.700_flux_fc-0.400_df-0.600_eps.h5

python PostProcessingUtils.py --make-movie $ezFile $epsFile



# ezFile=brd_wvg_cvt_diff_spc-with_cavity-5_sep-2.000_excite_fc-0.400_bw-0.300_flux_fc-0.400_df-0.600_ez.h5
# epsFile=brd_wvg_cvt_diff_spc-with_cavity-5_sep-2.000_excite_fc-0.400_bw-0.300_flux_fc-0.400_df-0.600_eps.h5

# python PostProcessingUtils.py --make-movie $ezFile $epsFile