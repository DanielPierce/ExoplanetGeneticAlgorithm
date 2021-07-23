
import configparser
import sys
import jsonpickle
from PSEUGA.common.PlanetarySystem import PlanetarySystem
import PSEUGA.src.Thesis as helpers
import numpy as np

def saveBestIndividual(individual, filepath):
    system = PlanetarySystem(individual)
    jsonsystem = jsonpickle.encode(system)
    
    dataFile = open(filepath, 'w')
    dataFile.write(jsonsystem)
    dataFile.close()

def saveRunData(runData, settings, timings, filepath):
    dataFile = open(filepath, 'w')
    dataFile.write(str(settings[0]) + "\n")
    dataFile.write(str(settings[1]) + "\n")
    dataFile.write(str(settings[2]) + "\n")
    dataFile.write(str(settings[3]) + "\n")
    dataFile.write(str(settings[4]) + "\n")
    dataFile.write(str(settings[5]) + "\n")
    dataFile.write("\n\n")
    dataFile.write(str(runData[0]) + "\n")
    dataFile.write(str(runData[1]) + "\n")
    dataFile.write("\n\n")
    dataFile.write(str(runData[2]) + "\n")
    dataFile.write(str(runData[3]) + "\n")
    dataFile.write(str(runData[4]) + "\n")
    dataFile.write(str(runData[5]) + "\n")
    dataFile.write("\n\n")
    for i in range(len(timings)):
        dataFile.write(str(timings[i]) + "\n")
    dataFile.close()

def saveLightcurve(individual, filepath):
    lc = helpers.generateLightcurve(individual)
    hdu = lc.to_fits(filepath, overwrite=True,TELESCOP='SIMULATION')

def savePopulation(pop, filepath):
    popArray = np.array(pop)
    np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
    np.savetxt(filepath, popArray, delimiter=',')