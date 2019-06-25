########################################################################
"""
Author: Jialun Luo
emailto:jl3562@cornell.edu
Last updated: 20190617

Description: This is a collection of file io utilities like (1) getting an unoccupied filename for saving data; (2)

"""
########################################################################
import os
import errno
import numpy as np # Only need this for testing the main
import matplotlib.pyplot as plt
import csv


# subfolderName should be a string without slashes or backslashes
def GetUnoccupiedFilename(checkFileBasename, extensionName, subfolderName=None, runID = 1, runIDDigitCount = 3, isDebugging = True):
    """ Given (or use a default) run ID, get an unoccupied filename for creating a new file    """
    checkingFilename = f'{checkFileBasename}_run-{runID:0{runIDDigitCount}}.{extensionName}'

    # print("DB: checkingFilename is " + checkingFilename)
    
    # If user wants to use a subfolder
    ### TODO: finish implementing the subfolder check. Currently not checking subfolders correctly
    if subfolderName:
        checkingFilename = PrependDirectory(subfolderName, checkingFilename)

    while True:
        checkingFile = os.path.isfile(checkingFilename)
        if checkingFile:
            # print(f'DB: {checkingFilename} + exists! Incrementing runID')
            runID += 1
            checkingFilename = f'{checkFileBasename}_run-{runID:0{runIDDigitCount}}.{extensionName}' # can change the conversion convention if you want more digits
            if(subfolderName):
                checkingFilename = PrependDirectory(subfolderName, checkingFilename)
            # print(f'DB: Try new filename: {checkingFilename}')
        else:
            goodID = runID
            goodFilenameWithExtension = checkingFilename
            break # exit the loop once a good filename is found
    return goodFilenameWithExtension, goodID

def GetUnoccupiedFoldername(checkFolderBasename, runID = 1, runIDDigitCount = 3, isDebugging = False, suppressInfo = False):
    """ Get an unused foldername for saving """
    checkFoldername = f'{checkFolderBasename}_run-{runID:0{runIDDigitCount}}'
    if (isDebugging): print(f'Check this name: {checkFoldername}')

    while True:
        checkingFolder = os.path.isdir(checkFoldername)
        if checkingFolder:
            runID += 1
            checkFoldername = f'{checkFolderBasename}_run-{runID:0{runIDDigitCount}}'
            if (isDebugging): print(f'DB: try new name: {checkFoldername}')
        else:
            goodID = runID
            goodFoldername = checkFoldername
            break # exit the loop once a good foldername is found

    return goodFoldername, goodID

        

def PrependDirectory(directoryName, filename):
    """Assemble a path name for a filename by prepending a directory"""
    return f'{directoryName}/{filename}'

def PavePath(filename):
    """creates folder(s) to a path if not already existing"""
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise




### TODO: implement a function to save np array data into some text file
def saveNPArray():
    pass


def csv2npArray(filename, skipColumn = 1):
    """ skipColumn is implemented as row = row[skipColumn:]"""
    with open(filename, 'r') as csvfile:
        fileReader = csv.reader(csvfile, delimiter = ',')

        for row in fileReader:
            row=row[skipColumn:]
            print(row)
            print(type(row))
            




            break





if __name__ == "__main__":
    """Just for testing"""

    """ Test GetUnoccupiedFoldername"""
    resultFolder = 'results'
    testFolderBasename = 'testFoldername'
    runID = 14
    foldername  = f'./{resultFolder}/{testFolderBasename}'
    goodName, goodID = GetUnoccupiedFoldername(foldername, runID, isDebugging=True)
    os.mkdir(goodName)


    """ Test csv2npArray(filename) """
    # csv2npArray('sampleOutput.txt')

    
    # print(testArray)

    """ test FileExistenceCheck() """ 
    # testFilename = "IAmAFileName"
    # txtExtension = 'txt'
    # goodFilenameWithExtension, goodID = GetUnoccupiedFilename(testFilename, txtExtension, "TestFolder", runID = 2)
    # np_Data = np.array([1, 2, 3 , 4 , 5])

    # print(f'Good filename is found: {goodFilenameWithExtension}')
    
    # PavePath(goodFilenameWithExtension)

    # with open(goodFilenameWithExtension, 'w') as f:
    #     np.savetxt(f, np_Data)
    #     pass
        # f.write('Test\n')

    ####################################



    # Random tests for learning python
    # str1 = "I am a string!"
    # int2 = 20
    # strAssem = f'{str1}_{int2:4}'
    # print(strAssem)
    # print(os.getcwd())
    # print(type(os.getcwd())) # string



    ####################################
