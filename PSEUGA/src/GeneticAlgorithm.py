
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


def runGA(processPool, numGenerations, pop):

    timings = []
    popStatsOverTime = []
    timeStatsOverTime = []
    # Variable keeping track of the number of generations
    g = 0
    # Begin the evolution
    while g < numGenerations:
        # A new generation
        g = g + 1
        print(f"-- Generation {g}/{numGenerations} --")

        tic = time.perf_counter()
        runGeneration(pop, processPool)
        toc = time.perf_counter()
        thisGenTime = toc - tic
        timings.append(thisGenTime)
        
        popStats = calculatePopStatistics(pop)
        popStatsOverTime.append(popStats)
        timeStats = calculateTimeStatistics(timings, numGenerations, g)
        timeStatsOverTime.append(timeStats)
        printGenerationData(popStats, timeStats, g)

    popStats = calculatePopStatistics(pop)
    runInfo = {'CXPB':CXPB, 'MUTPB':MUTPB, 'stats':popStats, 'popHistory':popStatsOverTime, 'timeHistory':timeStatsOverTime}
    return pop, timings, runInfo, hof

def initGA(populationSize, processPool):
    global hof
    global CXPB, MUTPB

    hof = tools.HallOfFame(3)
    pop = toolbox.population(n=populationSize)
    tic = time.perf_counter()
    tempPop = list(processPool.map(helpers.randomizeIndividual, pop))
    pop = tempPop
    toc = time.perf_counter()
    print(f"Randomized individuals in {toc - tic:0.4f} seconds")

    # Evaluate the entire population
    tic = time.perf_counter()
    fitnesses = list(processPool.map(toolbox.evaluate, pop))
    toc = time.perf_counter()
    print(f"Evaluated individuals in {toc - tic:0.4f} seconds")

    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    # CXPB  is the probability with which two individuals are crossed
    # MUTPB is the probability for mutating an individual
    CXPB, MUTPB = 0.9, 0.4
    return pop

def runGeneration(pop, processPool):
    # Select the next generation individuals
    # add some sort of elitism or hall of fame
    offspring = toolbox.select(pop, len(pop))
    # Clone the selected individuals
    offspring = list(map(toolbox.clone, offspring))       
    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values
    for mutant in offspring:
        if random.random() < MUTPB:
            toolbox.mutate(mutant)
            del mutant.fitness.values
    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = processPool.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    # shuffle positions of population bc of positional crossover
    pop[:] = offspring
    hof.update(pop)
    # Gather all the fitnesses in one list and print the stats
        
def printGenerationData(popStats, timeStats, g):
    print(f"  Min: {popStats['minFitness']}")
    print(f"  Max: {popStats['maxFitness']}")
    print(f"  Avg: {popStats['avgFitness']}")
    print(f"  Std: {popStats['stdFitness']}")
    print(f"Generation {g} complete in {timeStats['currentTime']:0.4f} seconds, averaging {timeStats['avg']:0.1f} seconds per")
    print(f"Estimate {timeStats['timeRemaining']} seconds remaining, done at {timeStats['estCompletionTime']}")

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