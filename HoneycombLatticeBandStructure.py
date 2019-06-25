import time
import sys
import os
import meep as mp
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize_scalar ### used for minimizing some function
from scipy.optimize import ridder ### used for finding a root of a function

## My scripts
import LJLUtilities
import MPBGraphingUtilities


# For debugging
# def trace(frame, event, arg):
#     print(f'{event}, {frame.f_code.co_filename}:{frame.f_lineno}')
#     return trace

# sys.settrace(trace)
#####################

""" TODO: save the parameters when you run stuff"""



def calculate(rIn, epsIn, geometryIsHole=False, k_interpIn = 8, resolutionIn = 64, num_bandsIn = 8,
                mesh_sizeIn = 3, runTEorTMorAll = 'all'):

    r = rIn
    eps = epsIn
    k_interp = k_interpIn # number of points between points specified when using interpolate
    resolution = resolutionIn
    num_bands = num_bandsIn
    mesh_size = mesh_sizeIn

    dielectricMaterial = mp.Medium(epsilon = eps)

    """ Specify single geometry pattern """
    if (geometryIsHole):
        """Then we drill holes in a materials (e.g. the geometry we specify will be air columns)"""
        backgroundMedium = dielectricMaterial
        geometry = [mp.Cylinder(r, material=mp.air), 
                    mp.Cylinder(r, center=mp.Vector3(1/3, 1/3, 0), material=mp.air)]

    else:
        backgroundMedium = mp.air
        geometry = [mp.Cylinder(r, material=dielectricMaterial), 
                    mp.Cylinder(r, center=mp.Vector3(math.sqrt(3)/3, 0, 0), material=dielectricMaterial)]

    # basis vector and how many n1 a1 + n2 a2, (n1, n2) is the size parameter
    geometry_lattice = mp.Lattice(  size = mp.Vector3(1,1),
                                    basis1 = mp.Vector3(math.sqrt(3)/2, 0.5),
                                    basis2 = mp.Vector3(math.sqrt(3)/2, -0.5))


    """ Repeat the pattern if needed"""
    geometry = mp.geometric_objects_lattice_duplicates(geometry_lattice, geometry)

    """ the k-points are by default in the basis of the lattice vectors"""
    Gamma = mp.Vector3()
    X = mp.Vector3(y=0.5)
    M = mp.Vector3(-1/3, 1/3)
    k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])
    # k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])

    ms = mpb.ModeSolver(
      geometry_lattice = geometry_lattice,
      geometry = geometry,
      k_points = k_points,
      resolution = resolution,
      num_bands = num_bands,
      mesh_size = mesh_size,
      default_material = backgroundMedium
      )

    # t0 = time.time()

    if (runTEorTMorAll.lower() == 'all'):
        ms.run()
    elif (runTEorTMorAll.lower() == 'run_tm'):
        ms.run_tm()
    elif (runTEorTMorAll.lower() == 'run_te'):
        ms.run_te()
    else:
        print("DB: not a valid run mode; do ms.run() instead")
        ms.run()

    """ Plot file parameters"""
    tentativeBasename = f'DispersionTrigLattice_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}'
    subfolderName = 'plots'
    graphSavefileExtension = 'svg'
    
    dielectricMapImageFilenameNoExtension = f'DielectricMap_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}'

    graphSavefileBasenameWithExtension, goodRunID = LJLUtilities.GetUnoccupiedFilename(tentativeBasename, graphSavefileExtension, subfolderName=subfolderName)
    dielectricMapImageFilename = f'./{subfolderName}/{dielectricMapImageFilenameNoExtension}_run-{goodRunID:03}.svg'



    """ Plot the dielectric map of the direct lattice """
    dielectricMap = ms.get_epsilon()
    dielectricMapFig, dielectricMapFigAxes = plt.subplots(1, 1)
    md = mpb.MPBData(rectify = True, periods=4, resolution = resolution)
    dielectricMap = md.convert(dielectricMap)
    dielectricImage = dielectricMapFigAxes.imshow(dielectricMap, cmap='gray_r')
    dielectricMapFig.colorbar(dielectricImage)
    dielectricMapFig.savefig(dielectricMapImageFilename)


    # print(f"total ime for both TE and TM bands: {time.time() - t0} seconds")
    bandDiagram = plt.figure()
    NKs = len(ms.all_freqs[:,1])
    for i in range(num_bands):
        plt.plot(np.arange(NKs), ms.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'Triangular basis:\n<{ms.geometry_lattice.basis1.x:.3f},{ms.geometry_lattice.basis1.y:.3f}>, <{ms.geometry_lattice.basis2.x:.3f},{ms.geometry_lattice.basis2.y:.3f}>\nr= {r:.3f}; eps= {eps:.4f}; Resolution= {resolution}; Mesh= {mesh_size}')
    plt.xlabel("k index")
    plt.ylabel("Frequency [2pi / a]")
    
    bandDiagram.savefig(f'./{graphSavefileBasenameWithExtension}')
    print(f'DB: Dispersion graph saved to ./{graphSavefileBasenameWithExtension}')
    print(f'DB: Dielectric map image saved to {dielectricMapImageFilename}')
    # plt.show()

    return (-1) * ms.retrieve_gap(1)

# def PlotDispersion(modeSolver, saveFigure = True):
#     """ Simply isolating the code out, might want to make this easier to use / more flexible"""

#     bandDiagram = plt.figure()
#     kpointCount = len(modeSolver.all_freqs[:,1])
#     for i in range(modeSolver.num_bands):
#         plt.plot(np.arange(kpointCount), modeSolver.all_freqs[:,i], marker='o', markersize = 3)
#     plt.suptitle(f'Triangular basis:\n{modeSolver.geometry_lattice.basis1:.3f}, {modeSolver.geometry_lattice.basis2:.3f}\nr={r:.3f}; eps={eps:.4f}')
#     plt.xlabel("k index")
#     plt.ylabel("Frequency [2pi / a]")

#     if(saveFigure):
#         bandDiagram.savefig(f'./plots/dispersion_TriangleLattice_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:0}.svg')
#     else:
#         bandDiagram.show()
#         pass

def SaveDispersionData(cylinderParam, frequencies):
    kPointsArchive = np.stack(cylinderParam.kpoints)

    resultFolderName = 'SimulationResults'
    kPointArxFilenameNoExtension = f'k-points_r-{cylinderParam.cylinderRadius}_eps-{cylinderParam.epsilon}'
    npyExtension = 'npy'
    # Initialize the kpoint file for existing file checking
    kPointArxFilenameWithExt, runID = LJLUtilities.GetUnoccupiedFilename(kPointArxFilenameNoExtension, npyExtension, subfolderName=resultFolderName)

    FreqsArxFilename = f'{resultFolderName}/Freqs_r-{cylinderParam.cylinderRadius}_eps-{cylinderParam.epsilon}.{npyExtension}'

    np.save(kPointArxFilenameWithExt, kPointsArchive)
    np.save(FreqsArxFilename, frequencies)


def calculateFirstBandGap(r, runTEorTM='run_te',  eps=1):
    bandgap = calculate(r, eps, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 4, resolutionIn = 64, runTEorTMorAll=runTEorTM)
    return bandgap


# def timeResolved:




if __name__ == "__main__":
    eps0 = 4.84
    # sys.settrace(trace)

    runTEorTMSet = 'run_te'
    result = minimize_scalar(lambda rrr: calculateFirstBandGap(rrr, runTEorTM=runTEorTMSet, eps=eps0) , method='bounded', bounds=[0.1, 0.2], options={'xatol': 0.0001})
    print(f'Optimal radius: {result.x}')
    

    # rtest=0.2
    # calculate(rtest, eps0, geometryIsHole = True, k_interpIn = 1, num_bandsIn = 8, mesh_sizeIn = 3, resolutionIn = 64)

    
    # for r in np.arange(0.30, 0.45, 0.4):
    # #     # calculate(r, eps0, num_bandsIn = 4, mesh_sizeIn = 3, resolutionIn=4)
    #     calculate(r, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 4, resolutionIn = 64)
