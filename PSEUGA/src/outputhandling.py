
import configparser
import sys
import jsonpickle
from PSEUGA.common.PlanetarySystem import PlanetarySystem

def saveBestIndividual(individual, filepath):
    system = PlanetarySystem(individual)
    jsonsystem = jsonpickle.encode(system)
    
    dataFile = open(filepath, 'w')
    dataFile.write(jsonsystem)
    dataFile.close()
    print(f"Storing {system.numActivePlanets} active planets")