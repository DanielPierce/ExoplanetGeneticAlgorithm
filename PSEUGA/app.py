

import datetime
import time
import PSEUGA.src.Thesis as helpers
import PSEUGA.common.Constants as const
import PSEUGA.src.GeneticAlgorithm as ga
import PSEUGA.vizualization.visualizer as viz
from PSEUGA.common.CustomLightcurve import CustomLightcurve

import lightkurve as lk
import numpy as np
import copy


import sys, getopt
import os

import multiprocessing as mp

import PSEUGA.vizualization.viewer as view

from PSEUGA.src.IOHandlers import InputHandler as Input, OutputHandler as Output

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

def initializeChildProcesses(target):
    helpers.targetCurve = target
    try:
        input = Input()
    except:
        input = Input.getInstance()
    input.getSettingsFromConf()
    input.GetIOFromInput(sys.argv)
    try:
        output = Output(input)
    except:
        output = Output.getInstance()

def saveResults(hof, runInfo, generationTimings, finalPopulation):
    input = Input.getInstance()
    runSettings = input.runSettings
    output = Output.getInstance()
    print("--+-- SAVING --+--")
    bestIndividual = hof[0]
    startTime = time.perf_counter()
    output.saveBestIndividual(bestIndividual)
    bestIndivTime = time.perf_counter()
    print(f"Save time for Best Individual: {bestIndivTime - startTime:2.4f}s")
    output.saveRunData(runInfo, runSettings, generationTimings)
    runDataTime = time.perf_counter()
    print(f"Save time for Run Data:        {runDataTime - bestIndivTime:2.4f}s")
    output.savePopulation(finalPopulation)
    popTime = time.perf_counter()
    print(f"Save time for Population:      {popTime - bestIndivTime:2.4f}s")
    output.saveLightcurve(bestIndividual)
    lightcurveTime = time.perf_counter()
    print(f"Save time for Lightcurve:      {lightcurveTime - popTime:2.4f}s")
    output.savePopHistory(runInfo['popHistory'])
    popHistTime = time.perf_counter()
    print(f"Save time for Pop History:     {popHistTime - lightcurveTime:2.4f}s")
    output.saveTimeHistory(runInfo['timeHistory'])
    timeHistTime = time.perf_counter()
    print(f"Save time for Time History:    {timeHistTime - popHistTime:2.4f}s")
    bestCurve = helpers.uniformSourceLightcurveAlgorithm(bestIndividual)
    customTarget = CustomLightcurve(helpers.targetCurve)
    viz.createLightcurvePlot(bestCurve, output.paths['outputFolderPath']+'GeneratedPlot.png')
    viz.createLightcurvePlot(customTarget, output.paths['outputFolderPath']+'TargetPlot.png')
    viz.createComparisonPlot(customTarget, bestCurve, output.paths['outputFolderPath']+'CompPlot.png')
    plotTime = time.perf_counter()
    print(f"Save time for plotting:        {plotTime - timeHistTime:2.4f}s")


def main():
    startupTime = time.perf_counter()
    try:
        input = Input()
    except:
        input = Input.getInstance()
    input.getSettingsFromConf()
    ioFromInputTime = time.perf_counter()
    print(f'Conf settings read in {ioFromInputTime - startupTime:2.4f} seconds')
    
    print('\n--+-- SETUP --+--')
    input.GetIOFromInput(sys.argv)
    setupTargetTime = time.perf_counter()
    print(f'IO in {setupTargetTime - ioFromInputTime:2.4f} seconds')
    try:
        output = Output(input)
    except:
        output = Output.getInstance()
    outputSetupTime = time.perf_counter()
    print(f'Output setup in {outputSetupTime - setupTargetTime:2.4f} seconds')

    try:
        helpers.targetCurve = lk.read(input.paths['inputFilePath'])
    except FileNotFoundError:
        print(f"ERROR: File {input.paths['inputFilePath']} does not exist\nEXITING")
        sys.exit(-1)
    poolTime = time.perf_counter()
    print(f'Set up target curve in {poolTime - setupTargetTime:2.4f} seconds')
    sys.stdout.flush()

    processPool = mp.Pool(input.runSettings['numChildProcesses'], initializeChildProcesses, [helpers.targetCurve])
    startGATime = time.perf_counter()
    print(f'Set up process pool in {startGATime - poolTime:2.4f} seconds')
    sys.stdout.flush()

    print('\n--+-- GA RUN --+--')

    print(f'Starting GA at {datetime.datetime.now()}')
    sys.stdout.flush()
    pop = ga.initGA( processPool)
    initGATime = time.perf_counter()
    finalPopulation, generationTimings, runInfo, hof = ga.runGA(processPool, input.runSettings['numGenerations'], pop)
    print("--+-- RUN COMPLETE --+--")
    print(f'Setup complete in {initGATime - startGATime:0.4f}')
    print(f"{input.runSettings['numGenerations']} generations complete in {sum(generationTimings):0.4f} seconds\n")

    saveResults(hof, runInfo, generationTimings, finalPopulation)
    print("\n--+-- ALL DONE --+--")
    

if __name__ == "__main__":
    main()
