

from datetime import datetime, timedelta
from deap.tools.support import Statistics
import Thesis as helpers
import Constants as const

import lightkurve as lk
import numpy as np
import random
import statistics

from deap import base
from deap import creator
from deap import tools

import sys, getopt
import os
import time

import multiprocessing as mp

inputFile = ""
curveOutputFile = ""
dataOutputFile = ""
populationOutputFile = ""
numThreads = 0
numIndividuals = 0
numGenerations = 0
timings = []

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
# Attribute generator 
toolbox.register("attr_int", random.randint, 0, const.MAXPLANETS)
# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, 
    toolbox.attr_int, 5 + const.ATTRPERPLANET * const.MAXPLANETS)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


toolbox.register("evaluate", helpers.evalOneMaxDist)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", helpers.mutation)
toolbox.register("select", tools.selTournament, tournsize=3)

def initializeThreads(target, timesteps):
    helpers.targetCurve = target
    const.skippedTimesteps = timesteps

def runGA():
    veryBeginning = time.perf_counter()
    threadPool = mp.Pool(numThreads, initializeThreads, [helpers.targetCurve, const.skippedTimesteps])


    pop = toolbox.population(n=numIndividuals)
    tic = time.perf_counter()
    tempPop = list(threadPool.map(helpers.randomizeIndividual, pop))
    pop = tempPop
    toc = time.perf_counter()
    print(f"Randomized individuals in {toc - tic:0.4f} seconds")

    # Evaluate the entire populationprint("first fitness")
    tic = time.perf_counter()
    fitnesses = list(threadPool.map(toolbox.evaluate, pop))
    toc = time.perf_counter()
    print(f"Evaluated individuals in {toc - tic:0.4f} seconds")

    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    global CXPB, MUTPB
    CXPB, MUTPB = 0.6, 0.4
    # Extracting all the fitnesses of 
    fits = [ind.fitness.values[0] for ind in pop]
        # Variable keeping track of the number of generations
    g = 0
    
    # Begin the evolution
    while g < numGenerations:
        tic = time.perf_counter()
        # A new generation
        g = g + 1
        print(f"-- Generation {g}/{numGenerations} --")
        # Select the next generation individuals
        # add some sort of elitism or hall of fame
        startSelect = time.perf_counter()
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))       
        endSelect = time.perf_counter()
        startCX = time.perf_counter()
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        endCX = time.perf_counter()
        startMutate = time.perf_counter()
        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        endMutate = time.perf_counter()
        # Evaluate the individuals with an invalid fitness
        startEval = time.perf_counter()
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = threadPool.map(toolbox.evaluate, invalid_ind)
        endEval = time.perf_counter()
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        # shuffle positions of population bc of positional crossover
        pop[:] = offspring

        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        global meanFitness
        meanFitness = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        global fitnessSTD
        fitnessSTD = abs(sum2 / length - meanFitness**2)**0.5
        global minFitness
        minFitness = min(fits)
        global maxFitness
        maxFitness = max(fits)
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % meanFitness)
        print("  Std %s" % fitnessSTD)
        toc = time.perf_counter()
        thisGen = toc - tic
        timings.append(thisGen)
        avg = statistics.mean(timings)
        gensLeft = numGenerations - g
        remainingSeconds = gensLeft * avg
        estimate = datetime.now() + timedelta(seconds=remainingSeconds)
        select = endSelect-startSelect
        cx = endCX-startCX
        mut = endMutate-startMutate
        evalu = endEval-startEval
        print(f"Generation {g} complete in {thisGen:0.4f} seconds, averaging {avg:0.1f} seconds per")
        print(f"Selection: {select}s, CX: {cx}s, Mut: {mut}s, Eval: {evalu}s")
        print(f"Estimate {remainingSeconds} seconds remaining, done at {estimate}")
        
    veryEnd = time.perf_counter()
    print(f"{g} generations complete in {veryEnd - veryBeginning:0.4f} seconds")
    return pop

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

    helpers.targetCurve = lk.read(inputFile).flatten()

def main():
    allToldStart = time.perf_counter()
    print(f"Start time: {datetime.now()}")
    getLightCurve()
    pop = runGA()
    printResults(pop)
    saveResults(pop)
    allToldEnd = time.perf_counter()
    print(f"All told, ran to completion in {allToldEnd - allToldStart:0.4f} seconds")

if __name__ == "__main__":
    main()
