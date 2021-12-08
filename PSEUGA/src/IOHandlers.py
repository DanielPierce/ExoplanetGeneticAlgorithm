

import configparser
from datetime import datetime
import shutil
import sys
import os
from deap.tools import init

import jsonpickle
from PSEUGA.common.PlanetarySystem import PlanetarySystem
import PSEUGA.src.Thesis as helpers
import numpy as np

import copy

from multiprocessing import current_process

import PSEUGA.vizualization.visualizer as viz
from PSEUGA.common.CustomLightcurve import CustomLightcurve

class InputHandler:
    __instance = None
    outputPath = 'PSEUGA/output/'
    paths = {}
    runSettings = {}
    runName = ''

    @staticmethod
    def getInstance():
        if InputHandler.__instance == None:
            raise Exception("Input handler must be initialized before an instance can be retrieved!")
        return InputHandler.__instance

    def __init__  (self):
        if InputHandler.__instance != None:
            raise Exception("Input handler is a singleton!")
        InputHandler.__instance = self
        
    #Get run settings from config file (store settings, return dict)
    #Get IO from input - get path to input file or set to default
    #Create list of paths (store, return dict)
    #Print paths
    #Print run settings

    def getSettingsFromConf(self):
        config = configparser.ConfigParser()
        configLocation = 'PSEUGA/data/config.ini'
        #print(f"Config at:  {configLocation}")
        config.read(configLocation)
        populationSize = int(config['GA']['population'])
        numGenerations = int(config['GA']['generations'])
        islandSwapGens = int(config['GA']['islandswap'])
        islandNumSwaps = int(config['GA']['islandnumswaps'])
        limbDarkeningType = config['ANALYSIS']['limbdarkening']
        timestepsToSkip = int(config['ANALYSIS']['stepstoskip'])
        numChildProcesses = int(config['ANALYSIS']['processes'])
        debug = config.getboolean('ANALYSIS', 'debug')
        outputGens = int(config['ANALYSIS']['outputgens'])
        self.runSettings.update({'populationSize':populationSize, 'numGenerations':numGenerations, 'islandSwapGens':islandSwapGens, 'islandNumSwaps':islandNumSwaps, 'limbDarkeningType':limbDarkeningType, 'timestepsToSkip':timestepsToSkip, 'numChildProcesses':numChildProcesses, 'debugMode':debug, 'outputGens':outputGens})
        #printConfig(runSettings)
        return self.runSettings

    def GetIOFromInput(self, args):
        if len(args) == 1 and self.runSettings['debugMode']:
            pathToInput, inputFileName = self.getInputInfoFromArg('DefaultFileTIC307210830C.fits')
            self.runName = 'default'
            self.addInputPathFromName(pathToInput, inputFileName)
            print('In debug with no input, running with inputs set to default')
            #printPaths(paths)
        if len(args) == 2:
            pathToInput, inputFileName = self.getInputInfoFromArg(args[1])
            self.addInputPathFromName(pathToInput, inputFileName)
            self.runName = args[1]
        if len(args) == 3:
            pathToInput, inputFileName = self.getInputInfoFromArg(args[1])
            self.addInputPathFromName(pathToInput, inputFileName)
            self.runName = args[2]
        if(self.runName == ''):
            print("Usage: python -m PSEUGA <input.fits> (optional)<outputsprefix>")
            sys.exit(-1)
        self.addOutputPathsFromName(self.outputPath)
    
    def getInputInfoFromArg(self, pathArg):
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

        
    def getInputInfoFromArg(self, pathArg):
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

    def addInputPathFromName(self, pathToInput, inputName):
        inputFilePath = pathToInput + inputName + '.fits'
        self.paths.update({'inputFilePath':inputFilePath})

    def addOutputPathsFromName(self, pathToOutput):
        outputFolder = pathToOutput + self.runName
        outputPrefix = outputFolder + '/' + self.runName
        if os.path.isdir(outputFolder) and current_process().name == 'MainProcess':
            shutil.rmtree(outputFolder)
            os.makedirs(pathToOutput + self.runName)

        fitsOutputPath = outputPrefix + "_OUT.fits"
        populationOutputPath = outputPrefix + "_POP.csv"
        runDataOutputPath = outputPrefix + "_RUNDATA.txt"
        bestIndivOutputPath = outputPrefix + "_BEST.json"
        popHistoryOutputPath = outputPrefix + "_PHISTORY.json"
        timeHistoryOutputPath = outputPrefix + "_THISTORY.json"

        newPaths = {'fitsOutputPath':fitsOutputPath, 'populationOutputPath':populationOutputPath, 'runDataOutputPath':runDataOutputPath, 
                 'bestIndivOutputPath':bestIndivOutputPath, 'popHistoryOutputPath':popHistoryOutputPath, 'timeHistoryOutputPath':timeHistoryOutputPath, 
                 'outputFolderPath':outputFolder + '/'}
        self.paths.update(newPaths)


class OutputHandler:
    __instance = None
    outputPath = 'PSEUGA/output/'
    paths = {}
    runSettings = {}
    runName = ''

    @staticmethod
    def getInstance():
        if OutputHandler.__instance == None:
            raise Exception("Output handler must be initialized before an instance can be retrieved!")
        return OutputHandler.__instance

    def __init__  (self, input):
        if OutputHandler.__instance != None:
            raise Exception("Output handler is a singleton!")
        self.runSettings = input.runSettings
        self.paths = input.paths
        self.runName = input.runName
        OutputHandler.__instance = self 

    def saveBestIndividual(self, individual):
        self.saveBestIndividualAt(individual, self.paths['bestIndivOutputPath'])

    def saveBestIndividualAt(self, individual, path):
        jsonSystem = jsonpickle.encode(individual)
        
        dataFile = open(path, 'w')
        dataFile.write(jsonSystem)
        dataFile.close()

    def saveRunData(self, runData, settings, timings):
        dataFile = open(self.paths['runDataOutputPath'], 'w')
        dataFile.write(str(settings['populationSize']) + "\n")
        dataFile.write(str(settings['numGenerations']) + "\n")
        dataFile.write(str(settings['limbDarkeningType']) + "\n")
        dataFile.write(str(settings['timestepsToSkip']) + "\n")
        dataFile.write(str(settings['numChildProcesses']) + "\n")
        dataFile.write(str(settings['debugMode']) + "\n")
        dataFile.write("\n\n")
        dataFile.write(str(runData['CXPB']) + "\n")
        dataFile.write(str(runData['MUTPB']) + "\n")
        dataFile.write("\n\n")
        dataFile.write(str(runData['stats']['minFitness']) + "\n")
        dataFile.write(str(runData['stats']['maxFitness']) + "\n")
        dataFile.write(str(runData['stats']['avgFitness']) + "\n")
        dataFile.write(str(runData['stats']['stdFitness']) + "\n")
        dataFile.write("\n\n")
        for i in range(len(timings)):
            dataFile.write(str(timings[i]) + "\n")
        dataFile.close()

    def saveLightcurve(self, individual):
        self.saveLightcurveAt(individual, self.paths['fitsOutputPath'])

    def saveLightcurveAt(self, individual, path):
        #lc = helpers.uniformSourceLightcurveAlgorithm(individual)
        lc = individual.lc
        #need to turn lc into a lightkurve object rather than custom lightcurve object, then remove .csv postfix
        #hdu = lc.to_fits(path, overwrite=True,TELESCOP='SIMULATION')
        time, flux = lc.toXY()
        arr = np.array([time,flux])
        #arr.tofile(path, sep = ',')
        np.savetxt(path, arr, delimiter = ',')


    def savePopulation(self, pop):
        self.savePopulationAt(pop, self.paths['populationOutputPath'])
    
    def savePopulationAt(self, pop, path):
        popList = [indiv.ps.ToList() for indiv in pop]
        #print(f"poplist: {type(popList)}, in poplist: {type(popList[0])}, values: {type(pop[0].fitness)}, len of vals: {len(pop[0].fitness)}")
        for i in range(len(popList)):
            try:
                popList[i].append(0)
                numValues = len(pop[i].fitness.values)
                #print(f"len fitness {numValues}, type: {type(pop[i].fitness.values)}")
                for j in range(numValues):
                    popList[i].append(pop[i].fitness.values[j])
            except Exception as e:
                print(f"save population error: {e}")
            try:
                popList[i].append(0)
                popList[i].append(pop[i].created)
            except Exception as e:
                print(f"save population error 2: {e}")
        popArray = np.array(popList)
        np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
        np.savetxt(path, popArray, delimiter=',')

    def savePopHistory(self, pop):
        jsonPop = jsonpickle.encode(pop)
        
        dataFile = open(self.paths['popHistoryOutputPath'], 'w')
        dataFile.write(jsonPop)
        dataFile.close()

    def saveTimeHistory(self, timingData):
        jsonTimingData = jsonpickle.encode(timingData)
        
        dataFile = open(self.paths['timeHistoryOutputPath'], 'w')
        dataFile.write(jsonTimingData)
        dataFile.close()

    def saveGenerationData(self, genNum, pop, hof):
        #rewrite this to not just save to default
        input = InputHandler.getInstance()
        numDigits = len(str(input.runSettings['numGenerations']))
        zerodGenNum = str(genNum).zfill(numDigits)
        filepath = 'PSEUGA/output/' + input.runName + '/historicaldata/gen' + zerodGenNum + '/'
        os.makedirs(filepath)
        self.savePopulationAt(pop, filepath + 'population.csv')
        try:
            customTarget = CustomLightcurve(helpers.targetCurve)
            viz.createLightcurvePlot(hof[0].lc, filepath+'GeneratedPlot.png')
            viz.createLightcurvePlot(customTarget, filepath+'TargetPlot.png')
            viz.createComparisonPlot(customTarget, hof[0].lc, filepath+'CompPlot.png') 
            self.saveLightcurveAt(hof[0], filepath + 'bestCurve.csv')
            self.saveBestIndividualAt(hof[0], filepath + 'bestIndiv.json')
        except Exception as e:
            print(f"save generation error: {e}")
