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



def calculate(rIn, epsIn, geometryIsHole=False, k_interpIn = 8, resolutionIn = 64, num_bandsIn = 8, mesh_sizeIn = 3, runTEorTMorAll = 'all', suppressInfo = False, isDebugging = True, centerCylinderRadius = 0.2, removeMiddleCylinder=False, removeMiddleRow=False, latticeType='triangle', superCellSize = 3):

    r = rIn
    eps = epsIn
    k_interp = k_interpIn # number of points between points specified when using interpolate
    resolution = resolutionIn
    num_bands = num_bandsIn
    mesh_size = mesh_sizeIn

    dielectricMaterial = mp.Medium(epsilon = eps)

    """ Specify single geometry pattern """
    if (geometryIsHole):
        """Then we drill holes in a materials (e.g. the cylinder we specify will be air columns)"""
        backgroundMedium = dielectricMaterial
        cylinder = [mp.Cylinder(r, material=mp.air)]
    else:
        backgroundMedium = mp.air
        cylinder = [mp.Cylinder(r, material=dielectricMaterial)]


    # basis vector and how many n1 a1 + n2 a2, (n1, n2) is the size parameter
    if(latticeType == 'triangle'):
        geometry_lattice = mp.Lattice(  size = mp.Vector3(superCellSize, superCellSize),
                                    basis1 = mp.Vector3(math.sqrt(3)/2, 0.5),
                                    basis2 = mp.Vector3(math.sqrt(3)/2, -0.5))

        """ the k-points are by default in the basis of the lattice vectors"""
        Gamma = mp.Vector3()
        X = mp.Vector3(y=0.5)
        M = mp.Vector3(-1/3, 1/3)

    elif(latticeType == 'square'):
        geometry_lattice = mp.Lattice(  
                                    size = mp.Vector3(superCellSize, superCellSize),
                                    basis1 = mp.Vector3(1, 0),
                                    basis2 = mp.Vector3(0, 1))
        Gamma = mp.Vector3()
        X = mp.Vector3(y=0.5)
        M = mp.Vector3(0.5, 0.5)

    k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])


    """ Repeat the pattern if needed"""
    geometryAssembly = mp.geometric_objects_lattice_duplicates(geometry_lattice, cylinder)
    if(isDebugging) :
        for gobject in geometryAssembly:
            print(f'{gobject} at {gobject.center} with material={gobject.material}')

    """ Create defect geometry"""
    if(removeMiddleCylinder):
        if(geometryIsHole):
            geometryAssembly.append(mp.Cylinder(r, material=dielectricMaterial))
        else:
            geometryAssembly.append(mp.Cylinder(r, material=mp.air))
    elif(removeMiddleRow):
        for i in range(superCellSize):
            shift = math.ceil(superCellSize/2)-1
            if(isDebugging): 
                print(f'removing {i}th cylinder at {mp.Vector3(x=1) * (i-shift)}')
                print(f'{i-shift}')
            if(geometryIsHole):
                geometryAssembly.append(mp.Cylinder(r, material=dielectricMaterial, center=mp.Vector3(x=1)* (i-shift)))
            else:
                geometryAssembly.append(mp.Cylinder(r, material=mp.air, center=mp.Vector3(x=1)  * (i-shift)))

    if(isDebugging) :
        for gobject in geometryAssembly:
            print(f'{gobject} at {gobject.center} with material={gobject.material}')
   
    # return
    # k_points = mp.interpolate(k_interp, [Gamma, X, M, Gamma])

    # t0 = time.time()

    """ save file parameters"""
    geometryType = f'{latticeType}LatCylAirHole'
    saveFolderBasename = f'{geometryType}_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}'

    # subfolderName = 'plots'
    resultFolderName = 'results'

    saveFolderBasename = f'./{resultFolderName}/{saveFolderBasename}'

    dielectricMapImageFilenameNoExtension = f'EpsMap_{geometryType}_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}'
    dispersionImageFilenameNoExtension = f'f-vs-k_{geometryType}_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}'
    graphSavefileExtension = 'svg'

    # graphSavefileBasenameWithExtension, goodRunID = LJLUtilities.GetUnoccupiedFilename(tentativeBasename, graphSavefileExtension, subfolderName=subfolderName)

    saveFolder, saveID = LJLUtilities.GetUnoccupiedFoldername(saveFolderBasename)
    os.mkdir(saveFolder)

    if(isDebugging): print(f'Saving results to {saveFolder}')


    dielectricMapImageFilename = f'{saveFolder}/{dielectricMapImageFilenameNoExtension}_run-{saveID:03}.{graphSavefileExtension}'
    dispersionImageFilename = f'{saveFolder}/{dispersionImageFilenameNoExtension}_run-{saveID:03}.{graphSavefileExtension}'

    dielectricMapNPFilename = f'{saveFolder}/{dielectricMapImageFilenameNoExtension}_run-{saveID:03}.npy'
    
    eigenfreqsFilename = f'{saveFolder}/freqs_{geometryType}_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}_run-{saveID:03}.npy'

    kpointFilename = f'{saveFolder}/k_{geometryType}_r-{r:.3f}_eps-{eps:.4f}_Nbands-{num_bands:03}_run-{saveID:03}.npy'

    print(f'Saving stdout and results to {saveFolder}')

    """ Save stdout into a log file"""
    originalStdOut = sys.stdout
    originalStdErr = sys.stderr
    sys.stdout = open(f'{saveFolder}/run.log', 'w')
    # sys.stderr = open(f'{saveFolder}/error.log', 'w')



    ms = mpb.ModeSolver(
      geometry_lattice = geometry_lattice,
      geometry = geometryAssembly,
      k_points = k_points,
      resolution = resolution,
      num_bands = num_bands,
      mesh_size = mesh_size,
      default_material = backgroundMedium
      )

    """ Choose to run modeSolver for tm/te/ all"""
    if (runTEorTMorAll.lower() == 'all'):
        ms.run()
    elif (runTEorTMorAll.lower() == 'run_tm'):
        ms.run_tm()
    elif (runTEorTMorAll.lower() == 'run_te'):
        ms.run_te()
    else:
        print("DB: not a valid run mode; do ms.run() instead")
        ms.run()


    """ Plot the dielectric map of the direct lattice """
    dielectricMap = ms.get_epsilon()
    dielectricMapFig, dielectricMapFigAxes = plt.subplots(1, 1)
    md = mpb.MPBData(rectify = True, resolution = resolution)
    dielectricMap = md.convert(dielectricMap)
    dielectricImage = dielectricMapFigAxes.imshow(dielectricMap, cmap='binary')
    dielectricMapFig.colorbar(dielectricImage)
    dielectricMapFig.savefig(dielectricMapImageFilename)

    """ Save some npy data"""
    np.save(eigenfreqsFilename, ms.all_freqs)
    np.save(dielectricMapNPFilename, np.asarray(dielectricMap))
    np.save(kpointFilename, np.stack(ms.k_points))




    # print(f"total ime for both TE and TM bands: {time.time() - t0} seconds")
    bandDiagram = plt.figure()
    NKs = len(ms.all_freqs[:,1])
    for i in range(num_bands):
        plt.plot(np.arange(NKs), ms.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'{latticeType} basis:\n<{ms.geometry_lattice.basis1.x:.3f},{ms.geometry_lattice.basis1.y:.3f}>, <{ms.geometry_lattice.basis2.x:.3f}, {ms.geometry_lattice.basis2.y:.3f}>\nr= {r:.3f}; eps= {eps:.4f}; Resolution= {resolution}; Mesh= {mesh_size}')
    plt.xlabel("k index")
    plt.ylabel("Frequency [2pi / a]")
    
    bandDiagram.savefig(f'{dispersionImageFilename}')
    

    """ Redirect messages from now on back to stdout etc."""
    sys.stdout = originalStdOut
    # sys.stderr = originalStdErr
    
    if (not suppressInfo):
        print(f'INFO: Dispersion graph saved to {dispersionImageFilename}')
        print(f'INFO: Dielectric map image saved to {dielectricMapImageFilename}')
    # plt.show()



    return ms.retrieve_gap(1)





def PlotDispersion(modeSolver, r, eps, saveFigure = True):
    """ Simply isolating the code out, might want to make this easier to use / more flexible"""
    """ Currently have to manually put in r and eps, maybe think of some struct to hold them?"""
    bandDiagram = plt.figure()
    kpointCount = len(modeSolver.all_freqs[:,1])
    for i in range(modeSolver.num_bands):
        plt.plot(np.arange(kpointCount), modeSolver.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'Lattice basis:\n<{modeSolver.geometry_lattice.basis1.x:.3f},{modeSolver.geometry_lattice.basis1.y:.3f}>, <{modeSolver.geometry_lattice.basis2.x:.3f},{modeSolver.geometry_lattice.basis2.y:.3f}>\nr= {r:.3f}; eps= {eps:.4f}; Resolution= {modeSolver.resolution}; Mesh= {modeSolver.mesh_size}')
    plt.xlabel("k index")
    plt.ylabel("Frequency [2pi / a]")

    if(saveFigure):
        bandDiagram.savefig(f'./plots/dispersion_TriangleLattice_r-{r:.3f}_eps-{eps:.4f}_Nbands-{modeSolver.num_bands:0}.svg')
    else:
        bandDiagram.show()
        pass

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


def calculateFirstBandGap(r, runTEorTM='run_te',  eps=1, bandgapMultipler = -1):
    bandgap = calculate(r, eps, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 2, resolutionIn = 32, mesh_sizeIn = 3,runTEorTMorAll=runTEorTM)
    return bandgapMultipler * bandgap 
    # Multliply by -1 for using minimizer


if __name__ == "__main__":
    eps0 = 4.84
    r0 = 0.38190

    # calculate(r0, eps0, geometryIsHole = True, k_interpIn = 8, num_bandsIn = 8, resolutionIn = 32, runTEorTMorAll='run_tm', superCellSize=1, latticeType='square')
    
    # sys.settrace(trace) """ In case of a seg fault; you can use this to track"""
  
    # for epsi in np.arange(3.24, 4.84, 0.01):
    #     print(f'Running TE band gap optimization for dielectric constant epsilon_r = {epsi}')
    #     runTEorTMSet = 'run_te'
    #     result = minimize_scalar(lambda rrr: calculateFirstBandGap(rrr, runTEorTM=runTEorTMSet, eps=epsi) , method='bounded', bounds=[0.1, 0.5], options=           {'xatol': 0.0001})  
    #     print(f'Optimal radius: {result.x}; gap = {(-1) * result.fun}')
    #     print(f'Finished running TE band optimization for dielectric constant epsilon_r = {epsi}')
    #     break
        

    # rtest=0.2
    # calculate(rtest, eps0, geometryIsHole = True, k_interpIn = 1, num_bandsIn = 8, mesh_sizeIn = 3, resolutionIn = 64)

    
    # for r in np.arange(0.30, 0.45, 0.4):
    # #     # calculate(r, eps0, num_bandsIn = 4, mesh_sizeIn = 3, resolutionIn=4)
    #     calculate(r, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 4, resolutionIn = 64)
    

    # for resolution in [32, 64, 128, 256, 512]:
    #     calculate(r0, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 2, resolutionIn = resolution, runTEorTMorAll='run_te')


    """ For our geometry, mesh doesn't seem to matter too much"""
    # for meshSize in [0, 3, 7, 12]:
    #     calculate(r0, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 2, resolutionIn = 32, mesh_sizeIn=meshSize, runTEorTMorAll='run_te')        

    # calculate(r0, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 50, resolutionIn = 16, mesh_sizeIn=0, runTEorTMorAll='run_te')          
# took 5 minutes for 16 kpoints per boundary
    # calculate(r0, eps0, geometryIsHole = True, k_interpIn = 4, num_bandsIn = 50, resolutionIn = 16, mesh_sizeIn=0, runTEorTMorAll='run_te')          

    """ Without any defect, but compute the bands with 25 cylinders in a unit cell """
    # calculate(r0, eps0, geometryIsHole = True, k_interpIn = 4, num_bandsIn = 75, resolutionIn = 16, mesh_sizeIn=0, runTEorTMorAll='run_te')          
    """ This one take 161 seconds to complete"""



    """ With point defect at the center, let's just take out the air hole in the middle"""
    # calculate(r0, eps0, geometryIsHole = True, k_interpIn = 8, num_bandsIn = 36, resolutionIn = 16, mesh_sizeIn=0, runTEorTMorAll='run_te', removeMiddleCylinder=True)          
    """ This one also took 162 seconds to complete"""


    """ With point defect at the center, take out the air hole in the middle row"""
    calculate(r0, eps0, geometryIsHole = True, k_interpIn = 4, num_bandsIn = 50, resolutionIn = 16, mesh_sizeIn=0, runTEorTMorAll='run_te', removeMiddleRow=True, superCellSize=5)          
    """ This one also took   seconds to complete"""