
import configparser
import sys
import os
import jsonpickle
from PSEUGA.common.PlanetarySystem import PlanetarySystem
from PSEUGA.common.CustomLightcurve import CustomLightcurve
import PSEUGA.src.Thesis as helpers
import numpy as np

#rewrite this as a class which can store the filepaths and things on its own
#input too, maybe together for sharing?

def saveGenerationData(genNum, runName, pop, hof):
    #rewrite this to not just save to default
    filepath = 'PSEUGA/output/default/historicaldata/gen' + str(genNum) + '/'
    os.makedirs(filepath)
    saveLightcurve(hof[0], filepath + 'bestCurve.csv')
    saveBestIndividual(hof[0], filepath + 'bestIndiv.json')
    savePopulation(pop, filepath + 'population.csv')



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
    times = []
    fluxes = []
    
    for step in lc.timeSteps:
        times.append(step.secondsFromEpoch)
        fluxes.append(step.flux)
    print("%s times in times"% len(times))
    np.savetxt(filepath,(times,fluxes), delimiter=',')
    pass

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
