

import PSEUGA.src.Thesis as helpers
import PSEUGA.common.Constants as const
import PSEUGA.src.GeneticAlgorithm as ga

import lightkurve as lk
import numpy as np


import sys, getopt
import os

import PSEUGA.src.inputhandling as input
import multiprocessing as mp
import PSEUGA.src.outputhandling as output

import PSEUGA.vizualization.viewer as view

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

def initializeChildProcesses(target, timesteps):
    helpers.targetCurve = target
    const.skippedTimesteps = timesteps

def main():
    runSettings = input.getSettingsFromConf()
    paths = input.getIOFromInput(sys.argv, runSettings['debugMode'])

    try:
        helpers.targetCurve = lk.read(paths['inputFilePath'])
    except FileNotFoundError:
        print(f"ERROR: File {paths['inputFilePath']} does not exist\nEXITING")
        sys.exit(-1)
    const.skippedTimesteps = runSettings['timestepsToSkip']

    threadPool = mp.Pool(runSettings['numChildProcesses'], initializeChildProcesses, [helpers.targetCurve, runSettings['timestepsToSkip']])
    finalPopulation, generationTimings, runInfo, hof = ga.runGA(threadPool, runSettings['populationSize'], runSettings['numGenerations'])
    print("Done with GA, starting save")
    bestIndividual = hof[0]
    output.saveBestIndividual(bestIndividual, paths['bestIndivOutputPath'])
    print("Saved best individual")
    output.saveRunData(runInfo, runSettings, generationTimings, paths['runDataOutputPath'])
    print("Saved run data")
    output.savePopulation(finalPopulation, paths['populationOutputPath'])
    print("Saved population")
    output.saveLightcurve(bestIndividual, paths['fitsOutputPath'])
    print("Saved .fits file")
    print("Done saving!\nExiting properly...")
    

if __name__ == "__main__":
    main()
