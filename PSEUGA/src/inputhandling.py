

import configparser
from datetime import datetime
import shutil
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
    runSettings = {'populationSize':populationSize, 'numGenerations':numGenerations, 'limbDarkeningType':limbDarkeningType, 'timestepsToSkip':timestepsToSkip, 'numChildProcesses':numChildProcesses, 'debugMode':debug}
    printConfig(runSettings)
    return runSettings

def getIOFromInput(args, debug):
    #print(f"args: {args}, len: {len(args)}")
    now = datetime.now()
    if len(args) == 1 and debug:
        pathToInput, inputFileName = getInputInfoFromArg('DefaultFileTIC307210830C.fits')
        paths = getPathsFromName(pathToInput, outputPath, inputFileName, 'default')
        print('In debug with no input, running with inputs set to default')
        #printPaths(paths)
        return paths
    if len(args) == 2:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        paths = getPathsFromName(pathToInput, outputPath, inputFileName, inputFileName)
        #printPaths(paths)
        return paths
    if len(args) == 3:
        pathToInput, inputFileName = getInputInfoFromArg(args[1])
        paths = getPathsFromName(pathToInput, outputPath, inputFileName, args[2])
        #printPaths(paths)
        return paths
    print("Usage: python -m PSEUGA <input.fits> (optional)<outputsprefix>")
    sys.exit(-1)
    
def getPathsFromName(pathToInput, pathToOutput, inputName, outputName):
        outputFolder = pathToOutput + outputName
        outputPrefix = outputFolder + '/' + outputName
        if os.path.isdir(outputFolder):
            shutil.rmtree(outputFolder)
        os.mkdir(pathToOutput + outputName)

        inputFilePath = pathToInput + inputName + '.fits'
        fitsOutputPath = outputPrefix + "_OUT.csv"
        populationOutputPath = outputPrefix + "_POP.csv"
        runDataOutputPath = outputPrefix + "_RUNDATA.txt"
        bestIndivOutputPath = outputPrefix + "_BEST.json"
        popHistoryOutputPath = outputPrefix + "_PHISTORY.json"
        timeHistoryOutputPath = outputPrefix + "_THISTORY.json"

        paths = {'inputFilePath':inputFilePath, 'fitsOutputPath':fitsOutputPath, 'populationOutputPath':populationOutputPath, 'runDataOutputPath':runDataOutputPath, 
                 'bestIndivOutputPath':bestIndivOutputPath, 'popHistoryOutputPath':popHistoryOutputPath, 'timeHistoryOutputPath':timeHistoryOutputPath, 
                 'outputFolderPath':outputFolder + '/'}
        return paths

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

def printPaths(paths):
    print(f"InputPath:    {paths['inputFilePath']}")
    print(f"OutputFits:   {paths['fitsOutputPath']}")
    print(f"Population:   {paths['populationOutputPath']}")
    print(f"Run data:     {paths['runDataOutputPath']}")
    print(f"Best indiv:   {paths['bestIndivOutputPath']}")
    print(f"Pop history:  {paths['runDataOutputPath']}")
    print(f"Time history: {paths['bestIndivOutputPath']}")

def printConfig(settings):
    print('--+-- CONF SETTINGS --+--')
    print(f"Pop  : {settings['populationSize']}")
    print(f"Gens : {settings['numGenerations']}")
    print(f"LimbD: {settings['limbDarkeningType']}")
    print(f"Skip : {settings['timestepsToSkip']}")
    print(f"Child: {settings['numChildProcesses']}")
    print(f"Debug: {settings['debugMode']}")
