#%%
import meep as mp
from meep import mpb
# help(mpb.MPBData)
help(mp.Lattice)


#%%




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
                mesh_sizeIn = 3, runTEorTMorAll = 'all', suppressInfo = False):

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
        geometry = [mp.Cylinder(r, material=mp.air)]
    else:
        backgroundMedium = mp.air
        geometry = [mp.Cylinder(r, material=dielectricMaterial)]


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

    """ save file parameters"""
    geometryType = 'TrigLatCylAirHole'
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

    dielectricMapImageFilename = f'./{saveFolder}/{dielectricMapImageFilenameNoExtension}_run-{saveID:03}.svg'
    dispersionImageFilename = f'./{saveFolder}/{dispersionImageFilenameNoExtension}_run-{saveID:03}.svg'
    """ Save some npy data"""
    # with dielectricMapImageFilenameNoExtension
    # np.save()


    """ Plot the dielectric map of the direct lattice """
    dielectricMap = ms.get_epsilon()
    dielectricMapFig, dielectricMapFigAxes = plt.subplots(1, 1)
    md = mpb.MPBData(rectify = True, periods=5, resolution = resolution)
    dielectricMap = md.convert(dielectricMap)
    dielectricImage = dielectricMapFigAxes.imshow(dielectricMap, cmap='gray_r')
    dielectricMapFig.colorbar(dielectricImage)
    dielectricMapFig.savefig(dielectricMapImageFilename)


    # print(f"total ime for both TE and TM bands: {time.time() - t0} seconds")
    bandDiagram = plt.figure()
    NKs = len(ms.all_freqs[:,1])
    for i in range(num_bands):
        plt.plot(np.arange(NKs), ms.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'Triangular basis:\n<{ms.geometry_lattice.basis1.x:.3f},{ms.geometry_lattice.basis1.y:.3f}>, <{ms.geometry_lattice.basis2.x:.3f}, {ms.geometry_lattice.basis2.y:.3f}>\nr= {r:.3f}; eps= {eps:.4f}; Resolution= {resolution}; Mesh= {mesh_size}')
    plt.xlabel("k index")
    plt.ylabel("Frequency [2pi / a]")
    
    bandDiagram.savefig(f'{dispersionImageFilename}')
    if (not suppressInfo):
        print(f'INFO: Dispersion graph saved to {dispersionImageFilename}')
        print(f'INFO: Dielectric map image saved to {dielectricMapImageFilename}')
    # plt.show()

    return ms
    # return ms.retrieve_gap(1)





def PlotDispersion(modeSolver, r, eps, saveFigure = True):
    """ Simply isolating the code out, might want to make this easier to use / more flexible"""
    """ Currently have to manually put in r and eps, maybe think of some struct to hold them?"""
    bandDiagram = plt.figure()
    kpointCount = len(modeSolver.all_freqs[:,1])
    for i in range(modeSolver.num_bands):
        plt.plot(np.arange(kpointCount), modeSolver.all_freqs[:,i], marker='o', markersize = 3)
    plt.suptitle(f'Triangular basis:\n<{modeSolver.geometry_lattice.basis1.x:.3f},{modeSolver.geometry_lattice.basis1.y:.3f}>, <{modeSolver.geometry_lattice.basis2.x:.3f},{modeSolver.geometry_lattice.basis2.y:.3f}>\nr= {r:.3f}; eps= {eps:.4f}; Resolution= {modeSolver.resolution}; Mesh= {modeSolver.mesh_size}')
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
    bandgap = calculate(r, eps, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 2, resolutionIn = 64, runTEorTMorAll=runTEorTM)
    return bandgapMultipler * bandgap 
    # Multliply by -1 for using minimizer


# if __name__ == "__main__":
#%%  
    
eps0 = 4.84
myCalculation=calculate(0.382, eps0, geometryIsHole = True, k_interpIn = 16, num_bandsIn = 2, resolutionIn = 64, runTEorTMorAll='run_te')

# # sys.settrace(trace)
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


#%%
import numpy as np
import meep as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import math


#%%
eps0 = 4.84
r0 = 0.382
f0 = 0.402
df = 0.2
computationCell = mp.Vector3(16,16)

backgroundMedium = mp.Medium(epsilon=eps0)
dielectricMaterial = backgroundMedium
airCylinder = mp.Cylinder(r0, material=mp.air)
hBNCylinder = mp.Cylinder(r0, material=dielectricMaterial)

basis1 = mp.Vector3(math.sqrt(3)/2, 0.5)
basis2 = mp.Vector3(math.sqrt(3)/2, -0.5)

geometryLattice = mp.Lattice(  size = mp.Vector3(15, 15),
                                basis1 = basis1,
                                basis2 = basis2)
# geometryLattice = mp.Lattice(  size = mp.Vector3(5, 5),
#                                 basis1 = mp.Vector3(1, 0),
#                                 basis2 = mp.Vector3(0, 1))


geometryAssembly = mp.geometric_objects_lattice_duplicates(geometryLattice, [airCylinder])

for geometricObject in geometryAssembly:
    geometricObject.center = mp.lattice_to_cartesian(geometricObject.center, geometryLattice)

""" make a defect point at the center"""
geometryAssembly.append(hBNCylinder)

# excitationSource = [mp.Source(mp.ContinuousSource(frequency=f0, width=20),
#                      component=mp.Hz,
#                      center=mp.Vector3(0, -5),
#                      size=mp.Vector3(2,0))]


""" Use a Gaussian source to excite """
excitationSource = [mp.Source(mp.GaussianSource(frequency=f0,fwidth=df),
                     component=mp.Ez,
                     center=mp.Vector3(0,0))]

pml_layers = [mp.PML(1.0)]

resolution = 16
sim = mp.Simulation(cell_size=computationCell,
                    boundary_layers=pml_layers,
                    geometry=geometryAssembly,
                    default_material=mp.Medium(epsilon=eps0),
                    sources=excitationSource,
                    resolution=resolution)



#%%
# sim.run(until=200)

# originalStdout = sys.stdout
# sys.stdout = open('./testOutput/run.log', 'w')

fcenterAnal = f0
framerate = 8
# moviePeriod = 2
fBandwidth=0.8 * f0
sim.reset_meep()
sim.use_output_directory('testOutput')
sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0, df)),
    mp.at_beginning(mp.output_epsilon),
    mp.to_appended("ez", mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
    until_after_sources = 150)
# sys.stdout.close()
# sys.stdout = originalStdout


#%%
framerate = 4
moviePeriod = 2

sim.reset_meep()
sim.use_output_directory('testOutput')
# sim.run(mp.in_volume(mp.Volume(center=mp.Vector3(), size=(computationCell))) ,
#     mp.at_beginning(mp.output_epsilon), 
#     mp.at_every(1/f0 / 10, mp.output_png(mp.Ez, "-Zc dkbluered")), until = 300)


sim.run(mp.in_volume(mp.Volume(center=mp.Vector3(), size=(computationCell)) ,
    mp.at_beginning(mp.output_epsilon),
    mp.to_appended("ez", mp.at_every(1 / f0 / framerate, mp.output_efield_z))), until = 50)
#%%
# commands to convert ez.h5 to png files for making a movie
# h5topng -t 0:479   -R  -S 2 -Zc dkbluered  -a yarg -A eps-000000.00.h5 ez.h5
# h5ls testOutput/
# convert ez*png ezMovie2.gif


#%%
eps_data = sim.get_array(center=mp.Vector3(), size=computationCell, component=mp.Dielectric)
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.show()


#%% Plot the fields 
ez_data = sim.get_array(center=mp.Vector3(), size=computationCell, component=mp.Hz)
MPLfig, MPLax = plt.subplots()
MPLax.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')

eAbsMax = np.max([np.abs(np.min(ez_data)),np.abs(np.max(ez_data))])
ezNorm = mpl.colors.Normalize(vmin=-eAbsMax, vmax=eAbsMax)

ezImage = MPLax.imshow(ez_data.transpose(), norm=ezNorm, interpolation='spline36', cmap='RdBu', alpha=0.9)
MPLfig.colorbar(ezImage)
# plt.axis('off')
plt.show()#%%


#%% Process HDF5 files in python? Use matplotlib to output
import h5py
workingDirectory = './results/meepTrigLatCylAirHole'
h5Filename = 'meepPointDefect-wvg_no_cavity_flux_1_ez.h5'

assembleh5Path = f'{workingDirectory}/{h5Filename}'

myfile = h5py.File(assembleh5Path, 'r')




#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    xdata.append(frame)
    ydata.append(np.sin(frame))
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 128),
                    init_func=init, blit=True)
plt.show()

#%%
