#!/opt/anaconda3/bin/python

"""
Author: Jialun Luo
"""

import h5py
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import time
import subprocess
import argparse
import glob

def PlotTransmission(refTrFile, measureTrFiles, workingDirectory = '.', skiprows=0, legends=None, plotDescription=None, saveFigure = False, saveName = 'UntitledSpectra'):
    referenceData = np.loadtxt(f'{workingDirectory}/{refTrFile}', delimiter=',', skiprows=skiprows, converters={0 : lambda s: np.nan})
    
    referenceData = referenceData.transpose()
    refFlux = referenceData[2]
    refFreq = referenceData[1]
    
    nPoints = len(refFlux)
    nCurves = len(measureTrFiles)
    print(f'We have {nPoints} points per curve')
    
    allTestFlux = np.zeros((nCurves, nPoints))
    allTestFreq = np.zeros((nCurves, nPoints))
    allRatios = np.zeros((nCurves, nPoints))
    
    """ 2nd column (index 1) is the frequency, 3rd column (index 2) is the field"""

    for idx, filename in enumerate(measureTrFiles):
        testSetupData = np.loadtxt(f'{workingDirectory}/{filename}', delimiter=',', skiprows=skiprows, converters={0 : lambda s: np.nan})
        testSetupData = testSetupData.transpose()
        allTestFlux[idx] = testSetupData[2]
        allTestFreq[idx] = testSetupData[1]
        allRatios[idx] = testSetupData[2]/refFlux

    plotCount = 2
    legendWidth = 4
    
    fig, axes = plt.subplots(plotCount, 1, figsize=(7 , 4 * plotCount))



    if(plotDescription is not None):
        fig.suptitle(plotDescription)
        # fig.text(0.05, 0.03, plotDescription, ha='left')

    if(legends is not None):
        """ TODO: clean up the legend code"""
        
    else:
        legends = ['With cavity', 'Without']
        axes[0].legend(legends)

    
    idx = 0    
    for ratio in allRatios:
        axes[0].plot(refFreq, ratio, label=legends[idx])
        idx += 1

    axes[0].set_ylabel('Flux ratio (to reference level)')
    axes[0].set_xlabel('Frequency [c/a]')
    
    axes[0].set_ylim([0, 1.1])
    
    axes[1].plot(refFreq, refFlux, '--', label='Without cavities')
    
    
    for i in range(nCurves):
        axes[1].plot(allTestFreq[i], allTestFlux[i], label=legends[i], linewidth = 0.5)

    axes[1].set_ylabel('Poynting Vector Flux')
    axes[1].set_xlabel('Frequency [c/a]')

    plt.tight_layout()
    # plt.tight_layout(pad=4, w_pad=0.5, h_pad=0.5)
    


    trSpectraFig, trSpectraFigAx = plt.subplots(1, 1, figsize=(8, 6))
    transmissionColorScale = mpl.colors.Normalize(vmin=0, vmax=1.)
    transmissionGraph = trSpectraFigAx.imshow(allRatios, aspect='auto', norm=transmissionColorScale,
                                              cmap='inferno',  origin='lower')

    trSpectraFig.colorbar(transmissionGraph)

    xtickCount = 10
    # xtickLabels = np.arange(refFreq[0], refFreq[-1], xtickCount)
    # xtickLabelPositions = np.arange()

    # plt.xticks(np.arange(min(refFreq), max(refFreq), xtickCount))


    ytickLabels = np.asarray(legends).astype(np.float)
    ytickLabelPositions = np.arange(ytickLabels.shape[0])
    plt.yticks(ytickLabelPositions, ytickLabels)

    if (saveFigure):
      # fig.savefig('TODOFIX_MY_NAME.svg')
      trSpectraFig.savefig(f'{saveName}.svg')
    plt.show()


def Plot2DParameterSweep():
    pass




def HDF2DImageTimeSeriesToMovie(h5filename, fps = 20, suppressInfo= False, overlayh5Filename=None, isDebugging=False):
    subprocess.run(['mkdir', '-p', 'tmp'], shell=True, capture_output=True)

    """ Read the time series image """
    data = None
    with h5py.File(h5filename, 'r') as h5file:
        for key in h5file.keys():
            print(f'Making movie with these data: {key}')
            data = np.asarray(h5file[key])
            if(isDebugging): 
                print(data)
                # print(data.shape[2])

    """ Read overlay image """
    overlayData = None
    if (overlayh5Filename != None):
        with h5py.File(overlayh5Filename, 'r') as overlayh5File:
            for key in overlayh5File.keys():
                print(f'Overlay image: Using {key}')
                overlayData = np.asarray(overlayh5File[key])
                # print(overlayData[:, :, 0])

    """ We are assuming time series of 2D images """
    frameCount = data.shape[2]
    if(isDebugging): print(f'We have {frameCount} frames to make')

    fieldAbsMax = np.max([np.abs(np.min(data)),np.abs(np.max(data))])
    fieldNorm = mpl.colors.Normalize(vmin=-fieldAbsMax, vmax=fieldAbsMax)
    encodeStartTime = time.time()

    images = []
    fig = plt.figure()
    fig.suptitle(h5filename)
  
    for i in range(frameCount):
        if(overlayh5Filename is not None):
            image2 = plt.imshow(overlayData[:, :, 0].transpose(),  cmap='gray_r', animated=True, alpha=0.5)
        image = plt.imshow(data[:, :, i].transpose(), norm=fieldNorm, animated=True, cmap='RdBu', alpha=0.7)
    
        images.append([image2, image])

        # plt.show()

    ani = animation.ArtistAnimation(fig, images, interval=1000/fps, blit=True, repeat_delay=1000)
    encodeEndTime = time.time()
    if(not suppressInfo): print(f'Animation took {(encodeEndTime - encodeStartTime):.5f} seconds to finish')

    FFwriter = animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=['-pix_fmt', 'yuv420p', '-crf', '0'])     
    #, 
    encodeStartTime = time.time()
    ani.save(f'{h5filename}.mp4', writer = FFwriter )
    encodeEndTime = time.time()
    if(not suppressInfo): print(f'Encoding took {(encodeEndTime - encodeStartTime):.5f} seconds')
    print(f'Video exported to\n{h5filename}.mp4')
    
    return 0

def PlotDielectricMap(epsH5filename, suppressInfo= False, isDebugging=False):
    
    with h5py.File(epsH5filename, 'r') as h5file:
        for key in h5file.keys(): # I am assuming there is only one key...
            epsMap = np.asarray(h5file[key])
            epsMap = epsMap[:,:,0].transpose()
    subprocess.run(['mkdir', '-p', 'tmp'], shell=True, capture_output=True)
    
    fig = plt.figure()
    epsImage = plt.imshow(epsMap, cmap='binary')
    fig.colorbar(epsImage)

    svgFilename = f'{epsH5filename[:-2]}svg' # CAUTION: -2 because we assume h5filename ended with .h5

    plt.savefig(svgFilename)
    print(f'Use the following command to open:')
    print(f'eog {svgFilename}&')


def GetParamValue(filename, paramName):
    explodedFilename = filename.split('_')
    return [s2.split('-')[1] for s2 in explodedFilename if 'sep' in s2][0]

if __name__ == '__main__':
    defaultWorkingDirectory = './results/meepTrigLatCylAirHole'
    defaultPythonScriptName = 'bridge_wvg_cavity'  
    
    parser = argparse.ArgumentParser(description="Utilities to process h5 outputs and flux tables from the meep simulation")
    parser.add_argument('-p', "--plot-transmission", nargs='+', help="Plot transmission given a reference flux (first file) and the test setup flux data (second file)")
    parser.add_argument("--make-movie", nargs=2, help="Make time domain simulation movie given the field evolution h5 file (must be the first file) (and optionally) an overlay dielectric map (second file)")
    parser.add_argument('-s', "--sim-file-prefix", help="Name of the python script which runs the meep simulation", default=defaultPythonScriptName)
    parser.add_argument('--set-wd', help="change the directory where files reside", default=defaultWorkingDirectory)

    args = parser.parse_args()

    print(f'Arguments are: {args}')

    resultDirectory = args.set_wd
    
    if args.plot_transmission:
        print(f'Here I should plot transmission and save it')
        print(f'the reference is {args.plot_transmission[0]}')
        print(f'The following are test setup data')

        paramName = 'sep'

        testDataFilenames = []
        legends = []
        testDataFilenameGlobbed  = False
        for s in args.plot_transmission[1:]:
            if( "*" in s ):
                print("* detected, using glob")

                allFilenames = glob.glob(f'{resultDirectory}/{s}')
                allFilenames.sort(key=lambda item: GetParamValue(item, paramName))
                for filename in allFilenames:
                    localFilename =filename.split('/')[-1]  #split and take the last element to get just the filename
                    testDataFilenames.append(localFilename)
                testDataFilenameGlobbed = True

        if ( not testDataFilenameGlobbed ):
            testDataFilenames = args.plot_transmission[1:]

        print(f'Plotting the following data:')
        for filename in testDataFilenames:
            print(filename)
            legends.append(GetParamValue(filename, paramName))
        print(f'Output contains parameter(s): {paramName} in list:')
        print(legends)
        PlotTransmission(args.plot_transmission[0], testDataFilenames, workingDirectory=defaultWorkingDirectory, legends = legends, plotDescription=args.plot_transmission[0])

    if args.make_movie:
        ezLocalName  = args.make_movie[0]
        epsLocalName = args.make_movie[1]

        h5Filename = f'{resultDirectory}/{ezLocalName}'
        overlayh5Filename = f'{resultDirectory}/{epsLocalName}'
        HDF2DImageTimeSeriesToMovie(h5Filename, overlayh5Filename=overlayh5Filename)