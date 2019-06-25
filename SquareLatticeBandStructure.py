import time
import subprocess
import sys
import math
import meep as mp
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt

def calculate(rIn, epsIn, k_interpIn = 8, resolutionIn= 32, num_bandsIn = 8):

    r = rIn
    eps = epsIn
    k_interp = k_interpIn # number of points between points specified when using interpolate
    resolution = resolutionIn
    num_bands = num_bandsIn

    Dielectrics = mp.Medium(epsilon = eps)

    geometry_lattice = mp.Lattice(  size = mp.Vector3(1,1))

    geometry = [mp.Cylinder(r, material=Dielectrics)]
    Gamma = mp.Vector3()
    X = mp.Vector3(0.5)
    M = mp.Vector3(0.5, 0.5)
    # k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])
    k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])


    

    ms = mpb.ModeSolver(
      geometry_lattice = geometry_lattice,
      geometry = geometry,
      k_points = k_points,
      resolution = resolution,
      num_bands = num_bands
      )

    t0 = time.time()

    ms.run_te()
    print(f"total ime for both TE and TM bands: {time.time() - t0} seconds")
    # ms.run_tm()
    bandDiagram = plt.figure()
    
    NKs = len(ms.all_freqs[:,1])
    for i in range(num_bands):
        plt.plot(np.arange(NKs), ms.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'Basis: {ms.geometry_lattice.basis1}, {ms.geometry_lattice.basis2}\nr={r:.3f}; eps={eps:.5f}')
    plt.xlabel("k index")
    plt.ylabel("Frequency [2pi/a ]")
    bandDiagram.savefig(f'./plots/dispersion_squareLattice_r-{r:.3f}_eps-{eps}.svg')
    plt.show()

if __name__ == "__main__":
    eps0=4.5

    calculate(0.3, eps0, num_bandsIn = 4)

    # for r in np.arange(0.2, 0.3, 0.1):
    #     calculate(r, eps0, num_bandsIn = 4)
