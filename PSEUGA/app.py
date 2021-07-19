

import PSEUGA.src.Thesis as helpers
import PSEUGA.common.Constants as const
import PSEUGA.src.GeneticAlgorithm as ga

import lightkurve as lk
import numpy as np


import sys, getopt
import os

import PSEUGA.src.inputhandling as input
import multiprocessing as mp


inputFile = ""
curveOutputFile = ""
dataOutputFile = ""
populationOutputFile = ""
numThreads = 0
numIndividuals = 0
numGenerations = 0
timings = []



def printResults(pop):
    fits = [ind.fitness.values[0] for ind in pop]
    bestIndiv = pop[fits.index(max(fits))]
    print("-----------BEST INDIVIDUAL'S STATS-----------")
    print("Fitness: %s" %helpers.evalOneMax(bestIndiv))
    for num in range(bestIndiv[const.NUMPLANETS]):
        print(f"  Planet {num}\n   Radius: {bestIndiv[num * const.ATTRPERPLANET + const.RADIUS]} | ECC: {bestIndiv[num * const.ATTRPERPLANET + const.ECC]} | SMA: {bestIndiv[num * const.ATTRPERPLANET + const.SMA]} | INC: {bestIndiv[num * const.ATTRPERPLANET + const.INC]} | LOAN: {bestIndiv[num * const.ATTRPERPLANET + const.LOAN]} | AOP: {bestIndiv[num * const.ATTRPERPLANET + const.AOP]} | MA: {bestIndiv[num * const.ATTRPERPLANET + const.MA]}")
    print(f"Num Planets: {bestIndiv[const.NUMPLANETS]}")
    print(f"Star Radius: {bestIndiv[const.STARRADIUS]} km")
    print(f"Star Mass  : {bestIndiv[const.STARMASS]} kg")
    print(f"Base Flux  : {bestIndiv[const.STARBASEFLUX]} e/s")
    print(f"Distance   : {bestIndiv[const.DISTANCE]} km")
    print("---------------BEST INDIVIDUAL---------------")
    print(bestIndiv)

def saveResults(pop):
    global inputFile
    global curveOutputFile
    global dataOutputFile
    

    if os.path.exists(curveOutputFile):
        os.remove(curveOutputFile)
    if os.path.exists(dataOutputFile):
        os.remove(dataOutputFile)

    fits = [ind.fitness.values[0] for ind in pop]
    bestIndiv = pop[fits.index(max(fits))]

    lc = helpers.generateLightcurve(bestIndiv)
    zero = lc.flux[0]
    last = lc.flux[-1]
    print(f"first: {zero}")
    print(f"last: {last}")
    hdu = lc.to_fits(curveOutputFile, overwrite=True,TELESCOP='SIMULATION')
    #hdu2 = helpers.targetCurve.to_fits('testoutput.txt', overwrite=True)
    #print("------------------------------------------------------------------------")
    #print(hdu[0].header)
    #rint("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    #print(hdu[1].header)
    #print("========================================================================")
    #print("------------------------------------------------------------------------")
    #print(hdu2[0].header)
    #print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    #print(hdu2[1].header)
    #print("========================================================================")

    dataFile = open(dataOutputFile, 'w')
    for characteristic in bestIndiv:
        dataFile.write(str(characteristic) + "\n")
    dataFile.write("\n\n")
    dataFile.write(str(numThreads) + "\n")
    dataFile.write(str(const.skippedTimesteps) + "\n")
    dataFile.write(str(numGenerations) + "\n")
    dataFile.write(str(numIndividuals) + "\n")
    dataFile.write("\n\n")
    dataFile.write(str(CXPB) + "\n")
    dataFile.write(str(MUTPB) + "\n")
    dataFile.write("\n\n")
    dataFile.write(str(meanFitness) + "\n")
    dataFile.write(str(fitnessSTD) + "\n")
    dataFile.write(str(minFitness) + "\n")
    dataFile.write(str(maxFitness) + "\n")
    dataFile.close()
    #reread = lk.read('testoutput.txt')
    #reread2 = lk.read(curveOutputFile)
    popArray = np.array(pop)
    print(len(pop))
    np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
    np.savetxt(populationOutputFile, popArray, delimiter=',')

def getLightCurve():
    global inputFile
    global curveOutputFile
    global dataOutputFile
    global populationOutputFile
    global numThreads
    global numIndividuals
    global numGenerations

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:t:s:p:g:")
    except getopt.GetoptError:
        print('GeneticAlgorithm.py -i <inputfile> -o <outputfilesprefix> -t <number of threads> -s <skipped timesteps> -p <population size> -g <number of generations>')
    for opt, arg in opts:
        if opt == '-i':
            inputFile = arg
        elif opt == '-o':
            curveOutputFile = arg + "curve.txt"
            dataOutputFile = arg + "data.txt" 
            populationOutputFile = arg + "pop.csv"
        elif opt == '-t':
            numThreads = int(arg)
        elif opt == '-s':
            const.skippedTimesteps = int(arg)
        elif opt == '-p':
            numIndividuals = int(arg)
        elif opt == '-g':
            numGenerations = int(arg)
    
    print("Input: %s" %inputFile)
    print("Curve: %s" %curveOutputFile)
    print("Data : %s" %dataOutputFile)
    print("Threads : %s" %numThreads)
    print("skippedTimesteps : %s" %const.skippedTimesteps)
    print("Population : %s" %numIndividuals)
    print("Generations : %s" %numGenerations)

    helpers.targetCurve = lk.read(inputFile)
    #helpers.targetCustom = CustomLightcurve(helpers.targetCurve)
    #helpers.targetCustomSorted = helpers.targetCustom.sortByFlux()


def initializeChildProcesses(target, timesteps):
    helpers.targetCurve = target
    const.skippedTimesteps = timesteps

def main():
    #allToldStart = time.perf_counter()
    #print(f"Start time: {datetime.now()}")
    #getLightCurve()
    #pop = runGA()
    #printResults(pop)
    #saveResults(pop)
    #allToldEnd = time.perf_counter()
    #print(f"All told, ran to completion in {allToldEnd - allToldStart:0.4f} seconds")
    populationSize, numGenerations, limbDarkeningType, timestepsToSkip, numChildProcesses, debug = input.getSettingsFromConf()
    inputFilePath, fitsOutputPath, populationOutputPath, runDataOutputPath = input.getIOFromInput(sys.argv, debug)
    
    try:
        helpers.targetCurve = lk.read(inputFilePath)
    except FileNotFoundError:
        print(f"ERROR: File {inputFilePath} does not exist\nEXITING")
        sys.exit(-1)
    const.skippedTimesteps = timestepsToSkip

    threadPool = mp.Pool(numChildProcesses, initializeChildProcesses, [helpers.targetCurve, timestepsToSkip])
    finalPopulation, generationTimings = ga.runGA(threadPool, populationSize, numGenerations)
    print("finished ga!!!!")

    

if __name__ == "__main__":
    main()
