

import configparser
from datetime import datetime
import sys
import os

outputPath = 'PSEUGA/output/'

def getSettingsFromConf():
    config = configparser.ConfigParser()
    configLocation = 'PSEUGA/data/config.ini'
    #print(f"Getting config from {configLocation}")
    config.read(configLocation)
    populationSize = int(config['GA']['population'])
    numGenerations = int(config['GA']['generations'])
    limbDarkeningType = config['ANALYSIS']['limbdarkening']
    timestepsToSkip = int(config['ANALYSIS']['stepstoskip'])
    numChildProcesses = int(config['ANALYSIS']['processes'])
    #print(config.sections())
    #print(f"Pop: {populationSize} which is a {type(populationSize)}")
    #print(f"Pop: {numGenerations} which is a {type(numGenerations)}")
    #print(f"Pop: {limbDarkeningType} which is a {type(limbDarkeningType)}")
    #print(f"Pop: {timestepsToSkip} which is a {type(timestepsToSkip)}")
    #print(f"Pop: {numChildProcesses} which is a {type(numChildProcesses)}")
    return populationSize, numGenerations, limbDarkeningType, timestepsToSkip, numChildProcesses

def getIOFromInput(args):
    #print(f"args: {args}, len: {len(args)}")
    now = datetime.now()
    if len(args) == 2:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath = getPathsFromName(pathToInput, outputPath, inputFileName)
        #print(f"InputPath:  {inputFilePath}")
        #print(f"OutputFits: {fitsOutputPath}")
        #print(f"Population: {populationOutputPath}")
        #print(f"Run data:   {runDataOutputPath}")
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath
    if len(args) == 3:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath = getPathsFromName(pathToInput, outputPath, args[2])
        #print(f"InputPath:  {inputFilePath}")
        #print(f"OutputFits: {fitsOutputPath}")
        #print(f"Population: {populationOutputPath}")
        #print(f"Run data:   {runDataOutputPath}")
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath
    print("Usage: python -m PSEUGA <input.fits> (optional)<outputsprefix>")
    sys.exit(-1)
    
def getPathsFromName(pathToInput, pathToOutput, name):
        inputFilePath = pathToInput + name + '.fits'
        fitsOutputPath = pathToOutput + name + "OUT.fits"
        populationOutputPath = pathToOutput + name + "POP.csv"
        runDataOutputPath = pathToOutput + name + "RUNDATA.txt"
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath

def getInputInfoFromArg(pathArg):
        pathToFile, fileNamePlusExt = os.path.split(pathArg)
        pathToFile = pathToFile + '/'
        fileName, extension = os.path.splitext(fileNamePlusExt)
        if(extension.lower() != '.fits'):
            print("Usage: input file must be in .fits format")
            sys.exit(-1)
        #print(f"pathToFile: {pathToFile}\nfileNamePlusExt: {fileNamePlusExt}\nfileName: {fileName}\nextension: {extension}")
        return pathToFile, fileName
