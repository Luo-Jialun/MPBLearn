import sys
import numpy as np


# An object for holding various parameters used in the simulation;
class MPB_ParamsCylinderRodGeometry():
    epsilon = 1
    cylinderRadius = 0.5 # in units of unit cell size (a)
    kpointCountPerSide = 8
    resolution = 64
    num_bands = 8
    mesh_size = 3

    # kPoints = np.array([])

    

    def __ini__(self):
        pass



def ImportDispersion(filename, params ):
    pass