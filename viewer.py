
#%%
#%matplotlib inline
from traitlets.traitlets import Integer
import Thesis as helpers
import Constants as const

import lightkurve as lk
import numpy as np
import random

import sys, getopt
import os
import time

import matplotlib


def getInputs():
    global inputFile
    global curveOutputFile
    global dataOutputFile

    #try:
    #    opts, args = getopt.getopt(sys.argv[1:], "i:c:d:")
    #except getopt.GetoptError:
    #    print('viewer.py -i <inputfile> -c <lightcurveoutputfile> -d <dataoutputfile>')
    #for opt, arg in opts:
    #    if opt == '-i':
    #        inputFile = arg
    #    elif opt == '-c':
    #        curveOutputFile = arg
    #    elif opt == '-d':
    #        dataOutputFile = arg 
    
    inputFile = "DefaultFileTIC307210830C.fits"
    curveOutputFile = "randomized.fits"
    dataOutputFile = "randomized.txt"
    
    print('-------------------------------------------------------')
    print("Input: %s" %inputFile)
    print("Curve: %s" %curveOutputFile)
    print("Data : %s" %dataOutputFile)
    const.skippedTimesteps = 1


def openFiles():
    helpers.targetCurve = lk.read(inputFile)
    #result = lk.read(curveOutputFile)
    with open(dataOutputFile) as f:
        individualStrings = f.read().splitlines()
    individual = [i for i in range(5 + const.ATTRPERPLANET * const.MAXPLANETS)]
    for i in range(const.MAXPLANETS):
        individual[i * const.ATTRPERPLANET + const.RADIUS] = float(individualStrings[i * const.ATTRPERPLANET + const.RADIUS])
        individual[i * const.ATTRPERPLANET + const.ECC] = float(individualStrings[i * const.ATTRPERPLANET + const.ECC])
        individual[i * const.ATTRPERPLANET + const.SMA] = float(individualStrings[i * const.ATTRPERPLANET + const.SMA])
        individual[i * const.ATTRPERPLANET + const.INC] = float(individualStrings[i * const.ATTRPERPLANET + const.INC])
        individual[i * const.ATTRPERPLANET + const.LOAN] = float(individualStrings[i * const.ATTRPERPLANET + const.LOAN])
        individual[i * const.ATTRPERPLANET + const.AOP] = float(individualStrings[i * const.ATTRPERPLANET + const.AOP])
        individual[i * const.ATTRPERPLANET + const.MA] = float(individualStrings[i * const.ATTRPERPLANET + const.MA])
        
    individual[const.NUMPLANETS] = int(individualStrings[const.NUMPLANETS])
    individual[const.STARRADIUS] = float(individualStrings[const.STARRADIUS])
    individual[const.STARMASS] = float(individualStrings[const.STARMASS])
    individual[const.STARBASEFLUX] = float(individualStrings[const.STARBASEFLUX])
    individual[const.DISTANCE] = float(individualStrings[const.DISTANCE])


    print('-------------------------------------------------------')
    inputVars = const.DISTANCE + 3
    print("Number of threads:................... %s" %individualStrings[inputVars])
    print("Skipped timesteps:................... %s" %individualStrings[inputVars + 1])
    print("Number of generations:............... %s" %individualStrings[inputVars + 2])
    print("Number of individuals:............... %s" %individualStrings[inputVars + 3])

    print('-------------------------------------------------------')
    recomboVars = inputVars + 6
    print("Crossover prob:...................... %s" %individualStrings[recomboVars])
    print("Mutation prob:....................... %s" %individualStrings[recomboVars + 1])
    
    print('-------------------------------------------------------')
    statsVars = recomboVars + 4
    print("Mean fitness:........................ %s" %individualStrings[statsVars])
    print("Fitness std:......................... %s" %individualStrings[statsVars+1])
    print("Min fitness:......................... %s" %individualStrings[statsVars+2])
    print("Max fitness:......................... %s" %individualStrings[statsVars+3])

    result = helpers.generateLightcurve(individual)
    result.plot()
    helpers.targetCurve.plot()
    flattened = helpers.targetCurve.flatten()
    flattened.plot()
    plotAvg = result[1]
    print(plotAvg)
    


def main():
    getInputs()
    openFiles()

if __name__ == "__main__":
    main()
# %%
