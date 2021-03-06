#!/opt/anaconda3/bin/python

"""
Author: Jialun Luo
Calculate time resolved field propagation of some photonic crystal structure

Note: on a different machine, check the #! statement at the beginning of the file
Parameters:
sidebankThickness,
separation, the distance between the centers of air cylinders of the center row
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


    hBNSidebankLeft = mp.Block(mp.Vector3(PMLThickness + sidebankThickness, computeCellSizeY), 
            material = materialHBN, 
            center = mp.Vector3((PMLThickness + sidebankThickness)/2-computeCellSizeX/2, 0))
    hBNSidebankRight = mp.Block(mp.Vector3(PMLThickness + sidebankThickness, computeCellSizeY), 
            material = materialHBN, 
            center = mp.Vector3(- (PMLThickness + sidebankThickness)/2 + computeCellSizeX/2, 0))
    hBNBridge = mp.Block(mp.Vector3(computeCellSizeX, bridgeWidth), material = materialHBN)

    geometryAssembly = []

    """ banks """
    geometryAssembly.append(hBNSidebankLeft) 
    geometryAssembly.append(hBNSidebankRight)
    geometryAssembly.append(hBNBridge)

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
            
            shift = separation
            
            geometryAssembly.append(mp.Cylinder(r, material=mp.air, center=mp.Vector3(x=1)* (i + shift)))
            geometryAssembly.append(mp.Cylinder(r, material=mp.air, center= -1 * mp.Vector3(x=1) * (i) ))

    # """ make a defect point at the center"""
    # geometryAssembly.append(hBNCylinder)            


    # for gobject in geometryAssembly:
        # print(f'{gobject} at {gobject.center} with material={gobject.material}')

    """ change the center into cartesian coordinates for meep"""
    for geometricObject in geometryAssembly:
        geometricObject.center = mp.lattice_to_cartesian(geometricObject.center, geometryLattice)


    if (pointSourceLocation== None):
        """ if the source location is not specified"""
        pointSourceLocation = mp.Vector3(- (computeCellSizeX/2 - PMLThickness * 1.1),0)


    """ Use a Gaussian source to excite """

    excitationSource = [mp.Source(mp.GaussianSource(frequency=fcen,fwidth=df),
                        component=mp.Ez,
                        # center=mp.Vector3())]
                        center=mp.lattice_to_cartesian(pointSourceLocation, geometryLattice))]

    # """ Use a continuous source to excite """
    # excitationSource = [mp.Source(mp.ContinuousSource(frequency=fcen, width=20),
    #                     component=mp.Ez,
    #                     center=mp.lattice_to_cartesian(pointSourceLocation, geometryLattice))]



    pml_layers = [mp.PML(PMLThickness)]

    resolution = 16
    sim = mp.Simulation(cell_size=computationCell,
                        boundary_layers=pml_layers,
                        geometry=geometryAssembly,
                        default_material=mp.air,
                        sources=excitationSource,
                        resolution=resolution)

    print(f'Done a preliminary setup')
    return sim

if __name__ == '__main__':

  """ Set up argparser here"""
  parser = argparse.ArgumentParser(description = 'Configure and run meep on certain geometry')

  pythonScriptName = 'brd_wvg_cvt_diff_spc'

  PMLThickness = 1.0
  eps0 = 4.84
  r0 = 0.3
  f0 = 0.4
  df = 0.6
  framerate = 8
  unitCellCountX = 40
  unitCellCountY = 5

  simDomainSizeX = 20
  simDomainSizeY = 10

  bridgeWidth = unitCellCountY
  # cavityUnitCellCount = 1

  # isMakingCavity = True
  # isMakingCavity = False

  """ Analysis parameters """
  nfreq = 4000 # number of frequencies at which to compute flux
  fluxDF = 0.6
  fluxFcen = 0.4
  
  harminvDf = 0.3
  harminvF0 = 0.4

  """ setup the geometry lattice """
  basis1 = mp.Vector3(1, 0)
  basis2 = mp.Vector3(0.5, math.sqrt(3)/2)

  geometryLattice = mp.Lattice(size = mp.Vector3(unitCellCountX, unitCellCountY),
                                  basis1 = basis1,
                                  basis2 = basis2)

  """ run with flux calculation """
  # print(f'{mp.lattice_to_cartesian(mp.Vector3(0,1,0), geometryLattice)}')
  fluxCutline = mp.FluxRegion(center=mp.Vector3(8,0), size=mp.Vector3(0, 2 * bridgeWidth), direction = mp.X)


  """ Parameter sweeps """
  cavityUnitCellCountQuery  = np.arange(1, 3, 1)
  isMakingCavityQuery       = [True]
  separationQuery           = np.arange(1, 2.05, 0.25)


  """ Uncomment the following lines to make just a line defect """
  # cavityUnitCellCountQuery  = [1]
  # isMakingCavityQuery       = [False]
  # separationQuery = [1] 
  # separationQuery           = [1, 1.2 , 1.5, 2]
  
  # epsQuery = [5, 7, 9, 11, 13]
  epsQuery  = [4.84]
  endProgramFlag=False 
  for eps1 in epsQuery:
    
    for isMakingCavity in isMakingCavityQuery:
      for cavityUnitCellCount in cavityUnitCellCountQuery:
        for separation in separationQuery:
          """ I don't seem to be able to calculate the flux in the commented out direction"""
          if(isMakingCavity):
            runDescription = f'wvg_with_cavity-{cavityUnitCellCount}_eps-{eps1}_sep-{separation:.3f}_excite_fc-{f0}_bw-{df}_flux_fc-{fluxFcen}_df-{fluxDF}'
          else:
            # endProgramFlag = True
            runDescription = f'wvg_no-cavity_eps-{eps1}_sep-0_excite_fc-{f0}_bw-{df}_flux_fc-{fluxFcen}_df-{fluxDF}'

          defaultResultFolder = '/home/mumaxbaby/Documents/jialun/MPBLearn/results/meepTrigLatCylAirHole'
          
          fieldFileBasename = f'{runDescription}_ez'
          epsMapFileBasename = f'{runDescription}_eps'
          fluxFileBasename = f'{runDescription}_flux'

          epsMapH5Filename = f'{defaultResultFolder}/{pythonScriptName}-{epsMapFileBasename}.h5'
          fieldH5Filename = f'{defaultResultFolder}/{pythonScriptName}-{fieldFileBasename}.h5'
          runLogFilename = f'{defaultResultFolder}/{runDescription}.log'
          initLogFilename = f'{defaultResultFolder}/{runDescription}.initialization.log'
          fluxDataFilename = f'{defaultResultFolder}/{fluxFileBasename}.csv'
          
          sim = setupSimulaion(eps = eps1, r = r0, fcen = f0, df = df, unitCellCountX = unitCellCountX, unitCellCountY = unitCellCountY, geometryLattice=geometryLattice, computeCellSizeX=simDomainSizeX, computeCellSizeY=simDomainSizeY, makeCavity=isMakingCavity, cavityUnitCellCount=cavityUnitCellCount, bridgeWidth = bridgeWidth, separation = separation)
          # sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df, unitCellCountX = unitCellCountX, unitCellCountY = unitCellCountY, geometryLattice=geometryLattice, computeCellSizeX=simDomainSizeX, computeCellSizeY=simDomainSizeY, makeCavity=isMakingCavity, cavityUnitCellCount=cavityUnitCellCount, bridgeWidth = bridgeWidth, separation = separation)
          sim.init_sim()
          sim.use_output_directory(defaultResultFolder)

          """ Try a narrow source outside of the patterned region"""
          # sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(), f0, df)), until_after_sources = 400)

          """ add_flux for calculate flux """
          trans = sim.add_flux(fluxFcen, fluxDF, nfreq, fluxCutline)

          with open(runLogFilename, 'w') as runLog:
            originalStdout = sys.stdout
            print(f'point stdout to {runLogFilename}')
            sys.stdout = runLog

            sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(simDomainSizeX/2 - PMLThickness * 1.1, 0), f0, harminvDf)),
                        mp.at_beginning(mp.to_appended(epsMapFileBasename, mp.output_epsilon)),
                        mp.to_appended(fieldFileBasename, mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
                        # until=500)
                        # mp.during_sources(mp.in_volume(vol, mp.to_appended(f'{runDescription}_ez-slice', mp.at_every(0.4, mp.output_efield_z)))),
                        until_after_sources = mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(0, 0), 1e-3))
                      
            sys.stdout = originalStdout
            print('print back to the real stdout')

          
          
          print(f'Run description: {runDescription}')
      
          with open(fluxDataFilename, 'w') as csvrecord:
            originalStdout = sys.stdout
            print('point stdout to a file')
            sys.stdout = csvrecord
            
            sim.display_fluxes(trans)  # print out the flux spectrum
            sys.stdout = originalStdout
            print('print back to the real stdout')
            print('Flux data saved at')
            print(fluxDataFilename)
            print('or')
            print(f'{fluxFileBasename}.csv')

          """ Convert the eps map h5 file into a png file"""
          PPU.PlotDielectricMap(epsMapH5Filename)
          
          """ End the program after one run of without holes (when isMakingCavity == False)""" 
          if(endProgramFlag):
            exit()
    
          
          
          
          
          


    # PPU.HDF2DImageTimeSeriesToMovie(fieldH5Filename, overlayh5Filename=epsMapH5Filename)


    # sim.run(mp.at_every(0.6 , mp.output_png(mp.Ez, "-Zc dkbluered")), until=200)

    # fcenterAnal = f0
    # moviePeriod = 2
    # fBandwidth=0.8 * f0
    # sim.reset_meep()
    # sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0, df)),
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
    #     sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0, df)),
    #             mp.at_beginning(mp.output_epsilon),
    #             mp.to_appended(f'ez_latCellCount-{latticeCellCount}', mp.at_every(1 / f0 / framerate, mp.output_efield_z)), 
    #             until_after_sources = 400)
    
    """ Excite with sources outside of the lattice """
    # sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0, df)), until_after_sources = 300)
    """ see the effect of using different excitation """
    # for df in np.arange(0.2 * f0, 0.4 * f0, f0/20):
    #     sim = setupSimulaion(eps = eps0, r = r0, fcen = f0, df = df)
    #     sim.init_sim()

    #     sim.run(mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0, df)), until_after_sources = 500)
    #     print(f'Using a Gaussian source with f_center={f0} and width of {df}; Check harminv')
        