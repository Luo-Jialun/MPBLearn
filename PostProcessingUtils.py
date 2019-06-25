"""
Author: Jialun Luo
"""

import h5py
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import subprocess

def PlotTransmission(refTrFile, measureTrFile, workingDirectory = '.', skiprows=0, legends=None):
    """ TODO: rename 'active' to something human readable... """
    referenceData = np.loadtxt(f'{workingDirectory}/{refTrFile}', delimiter=',', skiprows=skiprows, converters={0 : lambda s: np.nan})
    activeData = np.loadtxt(f'{workingDirectory}/{measureTrFile}', delimiter=',', skiprows=skiprows, converters={0 : lambda s: np.nan})
    
    """ 2nd column (index 1) is the frequency, 3rd column (index 2) is the field"""
    frequency = referenceData[:,1]
    activeDataField = activeData[:,2]
    refDataField = referenceData[:,2]
    print(f'Number of data: {np.size(refDataField)}')
    ratio = activeDataField/refDataField
    plotCount = 2
    fig, axes = plt.subplots(1, plotCount, figsize=(4 * plotCount, 4))
    
    axes[0].plot(frequency, ratio, 'r.-')
    axes[0].set_ylabel('Transmission')
    axes[0].set_xlabel('Frequency [c/a]')
    
    axes[1].plot(frequency, activeDataField, 'b.-', frequency, refDataField, 'g.-')
    if(legends!=None):
        axes[1].legend()
    else:
        axes[1].legend(['With cavity', 'Without'])

    axes[1].set_ylabel('Field strength')
    axes[1].set_xlabel('Frequency [c/a]')
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.show()


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
  
    for i in range(frameCount):
        if(overlayh5Filename!=None): 
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
    print(f'eog {svgFilename}')



if __name__ == '__main__':

    workingDirectory = './results/meepTrigLatCylAirHole'

    """ Test PlotDielectricMap """
    h5EpsilonFilename = 'meepPointDefect-wvg_with_no_cavity_exciationParam_fcen-0.43569_bw-0.05_fluxParam_fcen-0.435_df-0.1_eps.h5'

    PlotDielectricMap(f'{workingDirectory}/{h5EpsilonFilename}')



    # measurementTrFile = 'flux_testFlux_cavity_N-3.csv' 
    # refFile = 'flux_testFlux_no_cavity.csv'
    # plotTransmission(refFile, measurementTrFile, workingDirectory=workingDirectory)

    # print(np.divide(np.array([5, 4, 3, 2, 1]), np.array([1, 2, 3, 4, 5])))

    # fRef = 'wvg_with_no_cavity_1_fluxParam_fcen-0.435_df-0.05_flux.csv'
    # fActive = 'wvg_with_cavity_1_fluxParam_fcen-0.435_df-0.05_flux.csv'

    # fActive = 'wvg_with_cavity-1_fluxParam_fcen-0.435_df-0.1_flux.csv'
    # fActive = 'wvg_with_cavity-2_fluxParam_fcen-0.435_df-0.1_flux.csv'
    # fActive = 'wvg_with_cavity-3_exciationParam_fcen-0.43569_bw-0.05_fluxParam_fcen-0.435_df-0.1_flux.csv'
    # fRef = 'meepPointDefect-wvg_with_no_cavity_exciationParam_fcen-0.43569_bw-0.05_fluxParam_fcen-0.435_df-0.1_ez.h5'
    # PlotTransmission(fRef, fActive, workingDirectory=workingDirectory)

    
    # h5Filename = f'{workingDirectory}/meepPointDefect-wvg_with_no_cavity_exciationParam_fcen-0.43569_bw-0.05_fluxParam_fcen-0.435_df-0.1_ez.h5'
    # overlayh5Filename = f'{workingDirectory}/meepPointDefect-wvg_with_no_cavity_exciationParam_fcen-0.43569_bw-0.05_fluxParam_fcen-0.435_df-0.1_eps.h5'
    # HDF2DImageTimeSeriesToMovie(h5Filename, overlayh5Filename=overlayh5Filename, isDebugging = True)