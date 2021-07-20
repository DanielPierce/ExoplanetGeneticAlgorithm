

import configparser
from datetime import datetime
import sys
import os

outputPath = 'PSEUGA/output/'

def getSettingsFromConf():
    config = configparser.ConfigParser()
    configLocation = 'PSEUGA/data/config.ini'
    #print(f"Config at:  {configLocation}")
    config.read(configLocation)
    populationSize = int(config['GA']['population'])
    numGenerations = int(config['GA']['generations'])
    limbDarkeningType = config['ANALYSIS']['limbdarkening']
    timestepsToSkip = int(config['ANALYSIS']['stepstoskip'])
    numChildProcesses = int(config['ANALYSIS']['processes'])
    debug = bool(config['ANALYSIS']['debug'])
    #printConfig(populationSize, numGenerations, limbDarkeningType, timestepsToSkip, numChildProcesses, debug)
    return populationSize, numGenerations, limbDarkeningType, timestepsToSkip, numChildProcesses, debug

def getIOFromInput(args, debug):
    #print(f"args: {args}, len: {len(args)}")
    now = datetime.now()
    if len(args) == 1 and debug:
        pathToInput, inputFileName = getInputInfoFromArg('DefaultFileTIC307210830C.fits')
        inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath = getPathsFromName(pathToInput, outputPath, inputFileName, 'default')
        print('In debug, running default config')
        #printPaths(inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath)
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath
    if len(args) == 2:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath = getPathsFromName(pathToInput, outputPath, inputFileName, inputFileName)
        #printPaths(inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath)
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath
    if len(args) == 3:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath = getPathsFromName(pathToInput, outputPath, inputFileName, args[2])
        #printPaths(inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath)
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath
    print("Usage: python -m PSEUGA <input.fits> (optional)<outputsprefix>")
    sys.exit(-1)
    
def getPathsFromName(pathToInput, pathToOutput, inputName, outputName):
        inputFilePath = pathToInput + inputName + '.fits'
        fitsOutputPath = pathToOutput + outputName + "_OUT.fits"
        populationOutputPath = pathToOutput + outputName + "_POP.csv"
        runDataOutputPath = pathToOutput + outputName + "_RUNDATA.txt"
        bestIndivOutputPath = pathToOutput + outputName + "_BEST.json"
        return inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath

def getInputInfoFromArg(pathArg):
        pathToFile, fileNamePlusExt = os.path.split(pathArg)
        if(pathToFile != ""):
            pathToFile = pathToFile + '/'
        fileName, extension = os.path.splitext(fileNamePlusExt)
        if(extension.lower() != '.fits'):
            print("ERROR: input file must be in .fits format\nEXITING")
            sys.exit(-1)
        pathToFile = 'PSEUGA/data/' + pathToFile
        #print(f"pathToFile: {pathToFile}\nfileNamePlusExt: {fileNamePlusExt}\nfileName: {fileName}\nextension: {extension}")
        return pathToFile, fileName

def printPaths(inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath, bestIndivOutputPath):
    print(f"InputPath:  {inputFilePath}")
    print(f"OutputFits: {fitsOutputPath}")
    print(f"Population: {populationOutputPath}")
    print(f"Run data:   {runDataOutputPath}")
    print(f"Best Indiv: {bestIndivOutputPath}")

def printConfig(populationSize, numGenerations, limbDarkeningType, timestepsToSkip, numChildProcesses, debug):
    print(f"Pop  : {populationSize} which is a {type(populationSize)}")
    print(f"Gens : {numGenerations} which is a {type(numGenerations)}")
    print(f"LimbD: {limbDarkeningType} which is a {type(limbDarkeningType)}")
    print(f"Skip : {timestepsToSkip} which is a {type(timestepsToSkip)}")
    print(f"Child: {numChildProcesses} which is a {type(numChildProcesses)}")
    print(f"Debug : {debug} which is a {type(debug)}")
