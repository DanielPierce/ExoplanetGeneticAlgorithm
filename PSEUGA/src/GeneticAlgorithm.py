
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


def runGA(processPool, populationSize, numGenerations):
    print('starting GA')
    veryBeginning = time.perf_counter()

    timings = []

    pop = toolbox.population(n=populationSize)
    tic = time.perf_counter()
    tempPop = list(processPool.map(helpers.randomizeIndividual, pop))
    pop = tempPop
    toc = time.perf_counter()
    print(f"Randomized individuals in {toc - tic:0.4f} seconds")

    # Evaluate the entire populationprint("first fitness")
    tic = time.perf_counter()
    fitnesses = list(processPool.map(toolbox.evaluate, pop))
    toc = time.perf_counter()
    print(f"Evaluated individuals in {toc - tic:0.4f} seconds")

    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    global CXPB, MUTPB
    #increase prob of crossover? try 0.9
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
        fitnesses = processPool.map(toolbox.evaluate, invalid_ind)
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
    return pop, timings, [CXPB, MUTPB, meanFitness, fitnessSTD, minFitness, maxFitness]