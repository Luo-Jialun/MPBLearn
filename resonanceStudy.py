#!/opt/anaconda3/bin/python

"""
Author: Jialun Luo
Calculate time resolved field propagation of some photonic crystal structure

Note: on a different machine, check the #! statement at the beginning of the file
Parameters:
sidebankThickness,
separation - the distance between the centers of air cylinders of the center row
"""


import numpy as np
import meep as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import argparse
import PostProcessingUtils as PPU
import sys
import subprocess


def setupSimulaion(eps=1, r=0.2, fcen=0.4, df=0.2, unitCellCountX=20, unitCellCountY=5, computeCellSizeX=20, computeCellSizeY=10, doFlux = True, geometryLattice=None, makeCavity=False, cavityUnitCellCount=2, pointSourceLocation=None, PMLThickness=1.0, sidebankThickness = 1.0, bridgeWidth = 1.0, separation = 2.0):
  computationCell = mp.Vector3(computeCellSizeX, computeCellSizeY)

  materialHBN = mp.Medium(epsilon=eps)
  dielectricMaterial = materialHBN
  airCylinder = mp.Cylinder(r, material=mp.air)
  hBNCylinder = mp.Cylinder(r, material=dielectricMaterial)

  if(geometryLattice is not None):
    print('No lattice provided, setup triangle lattice...')
    basis1 = mp.Vector3(math.sqrt(3)/2, 0.5)
    basis2 = mp.Vector3(math.sqrt(3)/2, -0.5)

    geometryLattice = mp.Lattice(  size = mp.Vector3(unitCellCountX, unitCellCountY),
                                    basis1 = basis1,
                                    basis2 = basis2)


  # hBNSidebankLeft = mp.Block(mp.Vector3(PMLThickness + sidebankThickness, computeCellSizeY), 
  #         material = materialHBN, 
  #         center = mp.Vector3((PMLThickness + sidebankThickness)/2-computeCellSizeX/2, 0))
  # hBNSidebankRight = mp.Block(mp.Vector3(PMLThickness + sidebankThickness, computeCellSizeY), 
  #         material = materialHBN, 
  #         center = mp.Vector3(- (PMLThickness + sidebankThickness)/2 + computeCellSizeX/2, 0))
  hBNBridge = mp.Block(mp.Vector3(mp.inf, bridgeWidth, mp.inf), material = materialHBN)

  geometryAssembly = [hBNBridge]

  """ """
  airHoles = mp.geometric_objects_lattice_duplicates(geometryLattice, [airCylinder])
  for hole in airHoles:
    geometryAssembly.append(hole)
  
  """ make a defect line at y=0 """
  for i in range(unitCellCountX):
      shift = math.ceil(unitCellCountX/2)-1
      geometryAssembly.append(mp.Cylinder(r, material=dielectricMaterial, center=mp.Vector3(x=1)* (i-shift) ))

  """ make a cavity in the middle with N air cylinders in either direction"""
  if(makeCavity):
    shift = separation/2
    for i in range(cavityUnitCellCount):        
      geometryAssembly.append(mp.Cylinder(r, material=mp.air, center=mp.Vector3(x=1) * (i + shift)))
      geometryAssembly.append(mp.Cylinder(r, material=mp.air, center= -1 * mp.Vector3(x=1) * (i + shift) ))

  """ change the center into cartesian coordinates for meep"""
  # for geometricObject in geometryAssembly:
  #     geometricObject.center = mp.lattice_to_cartesian(geometricObject.center, geometryLattice)

  if (pointSourceLocation is None):
      """ if the source location is not specified"""
      pointSourceLocation = mp.Vector3(0, 0)

  """ Use a Gaussian source to excite """
  excitationSource = [mp.Source(mp.GaussianSource(frequency=fcen,fwidth=df),
                      component=mp.Ey,
                      center=pointSourceLocation, 
                      size=mp.Vector3(0, bridgeWidth))]

  pml_layers = [mp.PML(PMLThickness)]

  resolution = 20
  sim = mp.Simulation(cell_size=computationCell,
                      boundary_layers=pml_layers,
                      geometry=geometryAssembly,
                      default_material=mp.air,
                      sources=excitationSource,
                      resolution=resolution)

  return sim

""" Debug tool"""

def trace(frame, event, arg):
  print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
  return trace

if __name__ == '__main__':
  # sys.settrace(trace)
  """ Set up argparser here"""
  parser = argparse.ArgumentParser(description = 'Configure and run meep on certain geometry')

  pythonScriptName = 'brd_wvg_cvt_diff_spc'

  PMLThickness = 1.0
  
  eps0 = 4.84
  r0 = 0.38
  # f0 = 0.344086 # center frequency of the source
  framerate = 8
  unitCellCountX = 20
  unitCellCountY = 1

  simDomainSizeX = 14
  simDomainSizeY = 6

  bridgeWidth = 1.2

  """ Analysis parameters """
  nfreq = 4000 # number of frequencies at which to compute flux
  # nfreq = 500
  fluxFcen = 0.5
  fluxDF = 0.8
  
  harminvF0 = 0.25
  harminvDf = 0.2

  """ setup the geometry lattice """
  basis1 = mp.Vector3(1, 0)
  basis2 = mp.Vector3(0.5, math.sqrt(3)/2)

  geometryLattice = mp.Lattice(size = mp.Vector3(unitCellCountX, unitCellCountY),
                                  basis1 = basis1,
                                  basis2 = basis2)

  """ run with flux calculation """
  # print(f'{mp.lattice_to_cartesian(mp.Vector3(0,1,0), geometryLattice)}')
  fluxCutline = mp.FluxRegion(center=mp.Vector3( simDomainSizeX/2 - 1 * PMLThickness - 0.5, 0), size=mp.Vector3(0, 2* bridgeWidth))#, direction = mp.X)


  """ Parameter sweeps """
  cavityUnitCellCountQuery  = np.arange(1, 3, 1)
  isMakingCavityQuery       = [True, False]
  separationQuery           = np.arange(1.5, 1.6, 0.01)

  """ Setups for printing stdout into a file"""
  # originalStdout = sys.stdout


  """ Uncomment the following lines to make just a line defect """
  # cavityUnitCellCountQuery  = [3, 4, 5, 6]
  cavityUnitCellCountQuery = [ 3]
  # separationQuery = [1.64] 
  # separationQuery           = [1, 1.2 , 1.5, 2]
  
  # epsQuery = [5, 7, 9, 11, 13]
  # epsQuery  = [13]
  refIsCalculated=False 
  
  # exciteF0Query = np.arange(0.390, 0.490, 0.01)
  exciteF0Query = [0.4]
  df = 0.6 # bandwidth of the source (Gaussian frequency profile, 1 sigma frequency)

  
  harminvF0 = 0.3
  harminvDf = harminvF0
  ptSourceLocation = mp.Vector3(- (simDomainSizeX/2- 1 * PMLThickness) , 0)
  
  # ptSourceLocation = mp.Vector3(- (simDomainSizeX/2 - 1.5 * PMLThickness), 0)
  
  defaultResultFolder = '/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole'
  
  sim = None
  for f0 in exciteF0Query:
    
    for isMakingCavity in isMakingCavityQuery:
      """ End the current loop after one run of without holes (when isMakingCavity == False)""" 
      if(refIsCalculated):
        # refIsCalculated = False
        break
      
      for cavityUnitCellCount in cavityUnitCellCountQuery:
        for separation in separationQuery:
          
          if(isMakingCavity):
            runDescription = f'with_cavity-{cavityUnitCellCount}_r-{r0:.3f}_NRow-{unitCellCountY}_sep-{separation:.3f}_excite_fc-{f0:.3f}_bw-{df:.3f}_flux_fc-{fluxFcen:.3f}_df-{fluxDF:.3f}'
          else:
            refIsCalculated = True
            runDescription = f'no-cavity_r-{r0:.3f}_NRow-{unitCellCountY}_sep-0_excite_fc-{f0:.3f}_bw-{df:.3f}_flux_fc-{fluxFcen:.3f}_df-{fluxDF:.3f}'

          
          fieldFileBasename = f'{runDescription}_ez'
          epsMapFileBasename = f'{runDescription}_eps'
          fluxFileBasename = f'{runDescription}_flux'

          epsMapH5Filename = f'{defaultResultFolder}/{pythonScriptName}-{epsMapFileBasename}.h5'
          fieldH5Filename = f'{defaultResultFolder}/{pythonScriptName}-{fieldFileBasename}.h5'
          runLogFilename = f'{defaultResultFolder}/{runDescription}.log'
          initLogFilename = f'{defaultResultFolder}/{runDescription}.initialization.log'
          fluxDataFilename = f'{defaultResultFolder}/{fluxFileBasename}.csv'
                    
          sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df, unitCellCountX = unitCellCountX, unitCellCountY = unitCellCountY, geometryLattice=geometryLattice, computeCellSizeX=simDomainSizeX, computeCellSizeY=simDomainSizeY, makeCavity=isMakingCavity, cavityUnitCellCount=cavityUnitCellCount, bridgeWidth = bridgeWidth, separation = separation, pointSourceLocation=ptSourceLocation)
          # sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df, unitCellCountX = unitCellCountX, unitCellCountY = unitCellCountY, geometryLattice=geometryLattice, computeCellSizeX=simDomainSizeX, computeCellSizeY=simDomainSizeY, makeCavity=isMakingCavity, cavityUnitCellCount=cavityUnitCellCount, bridgeWidth = bridgeWidth, separation = separation)
          # sim.init_sim()
          sim.use_output_directory(defaultResultFolder)

          """ add_flux for calculate flux """
          trans = sim.add_flux(fluxFcen, fluxDF, nfreq, fluxCutline)

          if (isMakingCavity) :
            sim.run(
              # mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), harminvF0, harminvDf)),
              mp.at_beginning(mp.to_appended(epsMapFileBasename, mp.output_epsilon)),
              # mp.to_appended(fieldFileBasename, mp.at_every(1 / f0 / framerate, mp.output_efield_z)),
              # until=200)
              # mp.during_sources(mp.in_volume(vol, mp.to_appended(f'{runDescription}_ez-slice', mp.at_every(0.4, mp.output_efield_z)))),
              # until_after_sources = 500)
              until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(simDomainSizeX/2 - PMLThickness - 0.5, 0), 1e-3))
          else:
            sim.run(
              # mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), harminvF0, harminvDf)),
              mp.at_beginning(mp.to_appended(epsMapFileBasename, mp.output_epsilon)),
              # mp.to_appended(fieldFileBasename, mp.at_every(1 / f0 / framerate, mp.output_efield_z)),
              # until=50)
              # mp.during_sources(mp.in_volume(vol, mp.to_appended(f'{runDescription}_ez-slice', mp.at_every(0.4, mp.output_efield_z)))),
              until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(simDomainSizeX/2 - PMLThickness - 0.5, 0), 1e-3))

          # print(f'Run description: {runDescription}')
      
          with open(fluxDataFilename, 'w') as csvrecord:
            print(f'point stdout to {fluxDataFilename}')
            sys.stdout = csvrecord
            sim.display_fluxes(trans)  # print out the flux spectrum
            sys.stdout = open("/dev/stdout", "w")
            # print('print back to the real stdout')
            print('Flux data saved at')
            print(fluxDataFilename)
            print('or')
            print(f'{fluxFileBasename}.csv')
            
          """ closing log files """
          # initLog.close()
          # runLog.close()
          
          """ Convert the eps map h5 file into a png file"""
          PPU.PlotDielectricMap(epsMapH5Filename)