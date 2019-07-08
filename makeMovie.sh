#!/bin/bash

resultPath=/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole
scriptPath=/home/mumaxbaby/Documents/jialun/MPBLearn

ezFile=brd_wvg_cvt_diff_spc-wvg_with_cavity-1_eps-4.84_sep-1.000_excite_fc-0.4_bw-0.4_flux_fc-0.4_df-0.6_ez.h5
epsFile=brd_wvg_cvt_diff_spc-wvg_with_cavity-1_eps-4.84_sep-1.000_excite_fc-0.4_bw-0.4_flux_fc-0.4_df-0.6_eps.h5

python PostProcessingUtils.py --make-movie $ezFile $epsFile