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

  # def shiftedGObject(geom, shift):
  #     geom.center=geom.center+shift
  #     return geom

  if(geometryLattice==None):
      print('No lattice provided, setup one')
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

  """ banks """
  # geometryAssembly.append(hBNSidebankLeft) 
  # geometryAssembly.append(hBNSidebankRight)
  # geometryAssembly.append(hBNBridge)

  """ """
  airHoles = mp.geometric_objects_lattice_duplicates(geometryLattice, [airCylinder])
  for hole in airHoles:
    geometryAssembly.append(hole)
  
  # for gobject in geometryAssembly:
  #     print(f'{gobject} at {gobject.center} with material={gobject.material}')

  """ make a defect line at y=0 """
  for i in range(unitCellCountX):
      shift = math.ceil(unitCellCountX/2)-1
      
      geometryAssembly.append(mp.Cylinder(r, material=dielectricMaterial, center=mp.Vector3(x=1)* (i-shift) ))


  """ make a cavity in the middle with N air cylinders in either direction"""
  if(makeCavity):
      for i in range(cavityUnitCellCount):
          # shift = math.ceil(cavityUnitCellCount/2)-1        
          
          shift = separation/2
          geometryAssembly.append(mp.Cylinder(r, material=mp.air, center=mp.Vector3(x=1)* (i + shift)))
          geometryAssembly.append(mp.Cylinder(r, material=mp.air, center= -1 * mp.Vector3(x=1) * (i + shift) ))

  # """ make a defect point at the center"""
  # geometryAssembly.append(hBNCylinder)            


  # for gobject in geometryAssembly:
      # print(f'{gobject} at {gobject.center} with material={gobject.material}')

  """ change the center into cartesian coordinates for meep"""
  for geometricObject in geometryAssembly:
      geometricObject.center = mp.lattice_to_cartesian(geometricObject.center, geometryLattice)


  if (pointSourceLocation is None):
      """ if the source location is not specified"""
      pointSourceLocation = mp.Vector3(0, 0)

  """ Use a Gaussian source to excite """

  excitationSource = [mp.Source(mp.GaussianSource(frequency=fcen,fwidth=df),
                      component=mp.Ey,
                      center=pointSourceLocation, 
                      size=mp.Vector3(0, bridgeWidth))]

  # """ Use a continuous source to excite """
  # excitationSource = [mp.Source(mp.ContinuousSource(frequency=fcen, width=20),
  #                     component=mp.Ey,
  #                     center=mp.lattice_to_cartesian(pointSourceLocation, geometryLattice))]



  pml_layers = [mp.PML(PMLThickness)]

  resolution = 20
  sim = mp.Simulation(cell_size=computationCell,
                      boundary_layers=pml_layers,
                      geometry=geometryAssembly,
                      default_material=mp.air,
                      sources=excitationSource,
                      resolution=resolution)

  print(f'Done a preliminary setup')
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
  
  eps0 = 13
  r0 = 0.36
  f0 = 0.344086 # center frequency of the source
  framerate = 8
  unitCellCountX = 20
  unitCellCountY = 1

  simDomainSizeX = 12.4
  simDomainSizeY = 6

  bridgeWidth = 1.2

  """ Analysis parameters """
  nfreq = 4000 # number of frequencies at which to compute flux
  nfreq = 500
  fluxFcen = 0.25
  fluxDF = 0.2
  
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
  fluxCutline = mp.FluxRegion(center=mp.Vector3( simDomainSizeX/2 - 1 * PMLThickness - 0.5, 0), size=mp.Vector3(0, bridgeWidth))#, direction = mp.X)


  """ Parameter sweeps """
  cavityUnitCellCountQuery  = np.arange(1, 3, 1)
  isMakingCavityQuery       = [True, False]
  separationQuery           = np.arange(1.8, 2.1, 0.01)

  """ Setups for printing stdout into a file"""
  # originalStdout = sys.stdout


  """ Uncomment the following lines to make just a line defect """
  # cavityUnitCellCountQuery  = [3, 4, 5, 6]
  cavityUnitCellCountQuery = [6]
  separationQuery = [1.2] 
  # separationQuery           = [1, 1.2 , 1.5, 2]
  
  # epsQuery = [5, 7, 9, 11, 13]
  # epsQuery  = [13]
  refIsCalculated=False 
  
  # exciteF0Query = np.arange(0.390, 0.490, 0.01)
  exciteF0Query = [0.25]
  df = 0.2 # bandwidth of the source (Gaussian frequency profile, 1 sigma frequency)

  
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
        refIsCalculated = False
        
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
          sim.init_sim()
          sim.use_output_directory(defaultResultFolder)
          

          """ add_flux for calculate flux """
          trans = sim.add_flux(fluxFcen, fluxDF, nfreq, fluxCutline)

          if (isMakingCavity) :
            sim.run(
              mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), harminvF0, harminvDf)),
              mp.at_beginning(mp.to_appended(epsMapFileBasename, mp.output_epsilon)),
              # mp.to_appended(fieldFileBasename, mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
              # until=200)
              # mp.during_sources(mp.in_volume(vol, mp.to_appended(f'{runDescription}_ez-slice', mp.at_every(0.4, mp.output_efield_z)))),
              # until_after_sources = 500)
              until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(0, 0), 1e-5))
          else:
            sim.run(
              mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), harminvF0, harminvDf)),
              mp.at_beginning(mp.to_appended(epsMapFileBasename, mp.output_epsilon)),
              # mp.to_appended(fieldFileBasename, mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
              # until=50)
              # mp.during_sources(mp.in_volume(vol, mp.to_appended(f'{runDescription}_ez-slice', mp.at_every(0.4, mp.output_efield_z)))),
              until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(0, 0), 1e-5))

          print(f'Run description: {runDescription}')
      
          with open(fluxDataFilename, 'w') as csvrecord:
            print(f'point stdout to {fluxDataFilename}')
            sys.stdout = csvrecord
            sim.display_fluxes(trans)  # print out the flux spectrum
            sys.stdout = open("/dev/stdout", "w")
            print('print back to the real stdout')
            print('Flux data saved at')
            print(fluxDataFilename)
            print('or')
            print(f'{fluxFileBasename}.csv')
            
          """ closing log files """
          # initLog.close()
          # runLog.close()
          
          """ Convert the eps map h5 file into a png file"""
          PPU.PlotDielectricMap(epsMapH5Filename)
          

    # PPU.HDF2DImageTimeSeriesToMovie(fieldH5Filename, overlayh5Filename=epsMapH5Filename)


    # sim.run(mp.at_every(0.6 , mp.output_png(mp.Ey, "-Zc dkbluered")), until=200)

    # fcenterAnal = f0
    # moviePeriod = 2
    # fBandwidth=0.8 * f0
    # sim.reset_meep()
    # sim.run(mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), f0, df)),
    #     mp.at_beginning(mp.output_epsilon),
    #     mp.to_appended("ez", mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
    #     until_after_sources = 300)
    # sys.stdout.close()
    # sys.stdout = originalStdout

    # sim.run(mp.at_beginning(mp.output_epsilon), until_after_sources=1)


    """ Effect of the number of lattice cells """
    # for latticeCellCount in np.arange(5, 11, 100, dtype=int):
    #     sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df)
    #     sim.use_output_directory(defaultResultFolder)
    #     sim.init_sim()
    #     sim.run(mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), f0, df)),
    #             mp.at_beginning(mp.output_epsilon),
    #             mp.to_appended(f'ez_latCellCount-{latticeCellCount}', mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
    #             until_after_sources = 400)
    
    """ Excite with sources outside of the lattice """
    # sim.run(mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), f0, df)), until_after_sources = 300)
    """ see the effect of using different excitation """
    # for df in np.arange(0.2 * f0, 0.4 * f0, f0/20):
    #     sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df)
    #     sim.init_sim()

    #     sim.run(mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0), f0, df)), until_after_sources = 500)
    #     print(f'Using a Gaussian source with f_center={f0} and width of {df}; Check harminv')
        