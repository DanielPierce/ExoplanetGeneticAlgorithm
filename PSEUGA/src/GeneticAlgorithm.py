
from math import ceil
import math
from deap.tools.support import Statistics
from deap import base
from deap import creator
from deap import tools
import PSEUGA.src.Thesis as helpers
import PSEUGA.common.Constants as const

from datetime import datetime, timedelta
import time
import random
import statistics

import sys

from PSEUGA.src.IOHandlers import InputHandler as Input, OutputHandler as Output

from PSEUGA.common.CustomLightcurve import CustomLightcurve
from PSEUGA.common.PlanetarySystem import PlanetarySystem
import PSEUGA.vizualization.visualizer as viz

import numpy as np

import copy

creator.create("FitnessMin", base.Fitness, weights=(-1.0))
creator.create("MSEwXCORR", base.Fitness, weights=(-1.0, 1.0))
creator.create("XCORRcenter", base.Fitness, weights=(10.0,-100.0, -1.0))
creator.create("NormalizedXCORRcenter", base.Fitness, weights=(1.0, -1.0, -1.0))
creator.create("Individual", object, ps=PlanetarySystem, lc=CustomLightcurve, fitness=creator.NormalizedXCORRcenter, created=int)

toolbox = base.Toolbox()

toolbox.register("attr_int", random.randint, 0, const.MAXPLANETS)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, 1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


toolbox.register("evaluate", helpers.calculateTargetedTimesteps)
toolbox.register("mate", helpers.mate)
toolbox.register("mutate", helpers.mutation)
toolbox.register("select", tools.selNSGA2)



def runGA(processPool, numGenerations, pop):

    timings = []
    popStatsOverTime = []
    timeStatsOverTime = []
    input = Input.getInstance()
    output = Output.getInstance()
    # Variable keeping track of the number of generations
    g = 0
    # Begin the evolution
    while g < numGenerations:
        # A new generation
        g = g + 1
        print(f"--+-- GEN {g}/{numGenerations} --+--")

        tic = time.perf_counter()
        pop = runGeneration(pop, processPool, g)

        #if(input.runSettings['debugMode']):
            #output.saveGenerationData(g, pop, hof)
        #elif(g % input.runSettings['outputGens'] == 0):
            #output.saveGenerationData(g, pop, hof)

        toc = time.perf_counter()
        thisGenTime = toc - tic
        timings.append(thisGenTime)
        
        popStats = calculatePopStatistics(pop)
        popStatsOverTime.append(popStats)
        timeStats = calculateTimeStatistics(timings, numGenerations, g)
        timeStatsOverTime.append(timeStats)
        printGenerationData(popStats, timeStats, g)
        sys.stdout.flush()

    popStats = calculatePopStatistics(pop)
    runInfo = {'CXPB':CXPB, 'MUTPB':MUTPB, 'stats':popStats, 'popHistory':popStatsOverTime, 'timeHistory':timeStatsOverTime}
    return pop, timings, runInfo, hof

def initGA(processPool):
    global hof
    global CXPB, MUTPB

    input = Input.getInstance()
    output = Output.getInstance()
    hof = tools.HallOfFame(3)
    pop = toolbox.population(n=input.runSettings['populationSize'])
    tic = time.perf_counter()
    #output.saveGenerationData(-1, pop, hof)
    tempPop = list(processPool.map(helpers.randomizeIndividual, pop))
    #helpers.setToKep8b(tempPop[0])
    #helpers.setToKep8b(tempPop[1], -0.57 * math.pi) # -0.57 pi
    #helpers.setToKep8b(tempPop[2], 0.051 * math.pi)
    pop = tempPop
    #output.saveGenerationData(0, pop, hof)
    toc = time.perf_counter()
    print(f"Randomized individuals in {toc - tic:0.4f} seconds")

    # Evaluate the entire population
    tic = time.perf_counter()
    #fitnesses = list(processPool.map(toolbox.evaluate, (indiv.ps for indiv in pop)))
    fitnesses = []
    curves = []
    for fit, curve in processPool.map(toolbox.evaluate, (indiv.ps for indiv in pop)):
        fitnesses.append(fit)
        curves.append(curve)
    toc = time.perf_counter()
    print(f"Evaluated individuals in {toc - tic:0.4f} seconds")

    for ind, fit, curve in zip(pop, fitnesses, curves):
        ind.fitness.values = fit
        ind.lc = curve
        ind.created = -1

    # CXPB  is the probability with which two individuals are crossed
    # MUTPB is the probability for mutating an individual
    CXPB, MUTPB = 0.6, 0.4
    return pop

def runIslandGeneration(pop, processPool, genNum):
    # Select the next generation individuals
    # add some sort of elitism or hall of fame
    setNumIslands = 4
    parentIslands = []
    childIslands = []
    for i in range(setNumIslands):
        beginIndex = int(i/4 * len(pop))
        endIndex = int((i+1)/4 * len(pop))
        parentIslands.append(copy.deepcopy(pop[beginIndex:endIndex]))
        childIslands.append(pop[beginIndex:endIndex])
    offspring = []
    invalid_ind = []
    #deep copy, are originals still there?
    for i in range(len(parentIslands)):
    #mu/lambeda alg, generate the same number of kids (from best), then mutate/cx kids, then select new population from best of old pop + kids
        numToSelect = max(round(len(parentIslands[i])/4),2)
        thisIslandKids = toolbox.select(parentIslands[i], numToSelect)
        print(f"Tried to select {numToSelect}, actually selected {len(thisIslandKids)}")
        offspring.append(thisIslandKids)
        # Clone the selected individuals
        offspring[i] = list(map(toolbox.clone, offspring[i]))       
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[i][::2], offspring[i][1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
                child1.lc = None
                child2.lc = None
                child1.created = genNum
                child2.created = genNum
            else:
                toolbox.mutate(child1)
                toolbox.mutate(child2)
                del child1.fitness.values
                del child2.fitness.values
                child1.lc = None
                child2.lc = None
                child1.created = genNum
                child2.created = genNum
        # Evaluate the individuals with an invalid fitness
        invalid_ind.append([ind for ind in offspring[i] if not ind.fitness.valid])
    numEvald = 0
    allInvalid = []
    for i in range(len(invalid_ind)):
        numEvald = numEvald + len(invalid_ind[i])
        allInvalid = allInvalid + invalid_ind[i]
    fitnesses = []
    curves = []
    numEvald = numEvald + len(invalid_ind[i])

    for fit, curve in processPool.map(toolbox.evaluate, (indiv.ps for indiv in allInvalid)):
        fitnesses.append(fit)
        curves.append(curve)
    
    for ind, fit, curve in zip(allInvalid, fitnesses, curves):
        ind.fitness.values = fit
        ind.lc = curve

    islandSize = len(parentIslands[0])
    numIslands = len(parentIslands)
    islandPops = []
    input = Input.getInstance()
    for i in range(numIslands):
        islandPop = []
        popToSelectFrom = childIslands[i] + offspring[i]
        if genNum % input.runSettings['islandSwapGens'] == 0:
            #print("swapping!")
            for k in range(numIslands):
                if(i != k):
                    numToSwap = input.runSettings['islandNumSwaps']
                    popToSelectFrom = popToSelectFrom + toolbox.select(childIslands[k] + offspring[k], numToSwap)
        islandPop[:] = toolbox.select(popToSelectFrom, islandSize)
        islandPops = islandPops + islandPop
    # shuffle positions of population bc of positional crossover
    pop = islandPops
    hof.update(pop)
    # Gather all the fitnesses in one list and print the stats
    return pop

def runGeneration(pop, processPool, genNum):
    offspring = toolbox.select(pop, ceil(len(pop)/ 10))
    offspring = list(map(toolbox.clone, offspring))

    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values
            child1.lc = None
            child2.lc = None
            child1.created = genNum
            child2.created = genNum
        else:
            toolbox.mutate(child1)
            toolbox.mutate(child2)
            del child1.fitness.values
            del child2.fitness.values
            child1.lc = None
            child2.lc = None
            child1.created = genNum
            child2.created = genNum

    fitnesses = []
    curves = []

    for fit, curve in processPool.map(toolbox.evaluate, (indiv.ps for indiv in offspring)):
        fitnesses.append(fit)
        curves.append(curve)
    
    for ind, fit, curve in zip(offspring, fitnesses, curves):
        ind.fitness.values = fit
        ind.lc = curve

    pop = toolbox.select(pop + offspring, len(pop))

    hof.update(pop)
    print(f"pop len: {len(pop)}")
    return pop
        
def printGenerationData(popStats, timeStats, g):
    print(f"Min: {popStats['minFitness']}")
    print(f"Max: {popStats['maxFitness']}")
    print(f"Avg: {popStats['avgFitness']}")
    print(f"Std: {popStats['stdFitness']}")
    print(f"Generation {g} complete in {timeStats['currentTime']:0.2f} seconds, averaging {timeStats['avg']:0.2f} seconds per")
    print(f"Estimate {timeStats['timeRemaining']:0.2f} seconds remaining, done at {timeStats['estCompletionTime']}")

def calculateTimeStatistics(timings, numGenerations, g):
    avg = statistics.mean(timings)
    gensLeft = numGenerations - g
    timeRemaining = gensLeft * avg
    estCompletionTime = datetime.now() + timedelta(seconds=timeRemaining)
    timeStats = {'avg':avg, 'timeRemaining':timeRemaining, 'estCompletionTime':estCompletionTime, 'currentTime':timings[-1]}
    return timeStats

def calculatePopStatistics(pop):
    length = len(pop)
    fits = [ind.fitness.values[0] for ind in pop]
    avgFitness = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    stdFitness = abs(sum2 / length - avgFitness**2)**0.5
    minFitness = min(fits)
    maxFitness = max(fits)
    popStats = {'minFitness':minFitness, 'maxFitness':maxFitness, 'avgFitness':avgFitness, 'stdFitness':stdFitness}
    return popStats
