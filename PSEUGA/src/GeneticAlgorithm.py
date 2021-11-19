
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

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
#creator.create("Individual", list, fitness=creator.FitnessMin)
creator.create("Individual", object, ps=PlanetarySystem, lc=CustomLightcurve, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator 
toolbox.register("attr_int", random.randint, 0, const.MAXPLANETS)

# Structure initializers
#toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, 5 + const.ATTRPERPLANET * const.MAXPLANETS)
#toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.genome, toolbox.attr_int,toolbox.attr_int,toolbox.attr_int,toolbox.attr_int, toolbox.attr_int,toolbox.attr_int), n=1)

toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, 1)

toolbox.register("population", tools.initRepeat, list, toolbox.individual)


toolbox.register("evaluate", helpers.evalOneMaxDist)
toolbox.register("mate", helpers.mate)
toolbox.register("mutate", helpers.mutation)
toolbox.register("select", tools.selBest)



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
        runGeneration(pop, processPool)

        if(input.runSettings['debugMode']):
            output.saveGenerationData(g, pop, hof)
        elif(g % input.runSettings['outputGens'] == 0):
            output.saveGenerationData(g, pop, hof)

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
    output.saveGenerationData(-1, pop, hof)
    tempPop = list(processPool.map(helpers.randomizeIndividual, pop))
    pop = tempPop
    output.saveGenerationData(0, pop, hof)
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

    # CXPB  is the probability with which two individuals are crossed
    # MUTPB is the probability for mutating an individual
    CXPB, MUTPB = 0.6, 0.4
    return pop

def runGeneration(pop, processPool):
    # Select the next generation individuals
    # add some sort of elitism or hall of fame

    #mu/lambeda alg, generate the same number of kids (from best), then mutate/cx kids, then select new population from best of old pop + kids
    offspring = toolbox.select(pop, max(len(pop),2))
    # Clone the selected individuals
    offspring = list(map(toolbox.clone, offspring))       
    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values
            child1.lc = None
            child2.lc = None
        else:
            toolbox.mutate(child1)
            toolbox.mutate(child2)
            del child1.fitness.values
            del child2.fitness.values
            child1.lc = None
            child2.lc = None
    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = []
    curves = []
    for fit, curve in processPool.map(toolbox.evaluate, (indiv.ps for indiv in invalid_ind)):
        fitnesses.append(fit)
        curves.append(curve)
    #fitnesses = map(toolbox.evaluate, (indiv.ps for indiv in invalid_ind))
    
    for ind, fit, curve in zip(invalid_ind, fitnesses, curves):
        ind.fitness.values = fit
        ind.lc = curve

    popSize = len(pop)
    # shuffle positions of population bc of positional crossover
    pop[:] = toolbox.select(offspring + pop, popSize)
    hof.update(pop)
    # Gather all the fitnesses in one list and print the stats
        
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
