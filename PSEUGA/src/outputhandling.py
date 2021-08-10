
import configparser
import sys
import jsonpickle
from PSEUGA.common.PlanetarySystem import PlanetarySystem
import PSEUGA.src.Thesis as helpers
import numpy as np

def saveBestIndividual(individual, filepath):
    system = PlanetarySystem(individual)
    jsonSystem = jsonpickle.encode(system)
    
    dataFile = open(filepath, 'w')
    dataFile.write(jsonSystem)
    dataFile.close()

def saveRunData(runData, settings, timings, filepath):
    dataFile = open(filepath, 'w')
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

def saveLightcurve(individual, filepath):
    lc = helpers.uniformSourceLightcurveAlgorithm(individual)
    #hdu = lc.to_fits(filepath, overwrite=True,TELESCOP='SIMULATION')

def savePopulation(pop, filepath):
    popArray = np.array(pop)
    np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
    np.savetxt(filepath, popArray, delimiter=',')

def savePopHistory(pop, filepath):
    jsonPop = jsonpickle.encode(pop)
    
    dataFile = open(filepath, 'w')
    dataFile.write(jsonPop)
    dataFile.close()

def saveTimeHistory(timingData, filepath):
    jsonTimingData = jsonpickle.encode(timingData)
    
    dataFile = open(filepath, 'w')
    dataFile.write(jsonTimingData)
    dataFile.close()
