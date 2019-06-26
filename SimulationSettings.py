"""
Author: jl3562@cornell.edu
Provides utilities to parse .ini files that contain simulation parameters
"""

import configparser

def ReadSimulationParam(iniFilename):
  config = configparser.ConfigParser()
  config.read(iniFilename)
  return config
  
if __name__ == '__main__':
  iniFilename = f'exampleSimulationSetting.ini'
  simParam = ReadSimulationParam(iniFilename)
  print(f'Read {iniFilename}')
  print(simParam.sections())
  print(type(simParam))
  print(simParam['Material'])
  for section in simParam.sections():
    print(f'Printing each section:')
    print(section)
    print(f'This is the type of a section: {type(section)}')

    print(f'Type of config[section] is: {type(simParam[section])}')
    for paramName in simParam[section]:
      print(paramName)
      print(f'Param: {paramName}; Value: {simParam[section][paramName]}')
      

