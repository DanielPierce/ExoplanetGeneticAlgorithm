

import Thesis as helpers
import Constants as const

import lightkurve as lk
import numpy as np
import random

from deap import base
from deap import creator
from deap import tools

import sys, getopt
import os

inputFile = ""
curveOutputFile = ""
dataOutputFile = ""

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
# Attribute generator 
toolbox.register("attr_int", random.randint, 0, const.MAXPLANETS)
# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, 
    toolbox.attr_int, 2 + const.ATTRPERPLANET * const.MAXPLANETS)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


toolbox.register("evaluate", helpers.evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", helpers.mutation, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

def runGA():
    pop = toolbox.population(n=5)
    # Evaluate the entire populationprint("first fitness")
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    CXPB, MUTPB = 0.5, 0.2
    # Extracting all the fitnesses of 
    fits = [ind.fitness.values[0] for ind in pop]
        # Variable keeping track of the number of generations
    g = 0
    
    # Begin the evolution
    while g < 1:
        # A new generation
        g = g + 1
        print("-- Generation %i --" % g)
        # Select the next generation individuals
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
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        pop[:] = offspring

        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
    return pop

def printResults(pop):
    fits = [ind.fitness.values[0] for ind in pop]
    bestIndiv = pop[fits.index(max(fits))]
    print("Num planets: %s" %bestIndiv[const.NUMPLANETS])
    print("Fitness: %s" %helpers.evalOneMax(bestIndiv))
    for num in range(bestIndiv[const.NUMPLANETS]):
        print("  Planet %s: %f" %(num, bestIndiv[num * const.ATTRPERPLANET]))
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

    lc = helpers.generateLightCurve(bestIndiv)
    lc.to_fits(curveOutputFile, overwrite=True)

    dataFile = open(dataOutputFile, 'w')
    for characteristic in bestIndiv:
        dataFile.write(str(characteristic) + "\n")
    dataFile.close()

def getLightCurve():
    global inputFile
    global curveOutputFile
    global dataOutputFile
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:c:d:")
    except getopt.GetoptError:
        print('GeneticAlgorithm.py -i <inputfile> -c <lightcurveoutput> -d <dataoutput>')
    for opt, arg in opts:
        if opt == '-i':
            inputFile = arg
        elif opt == '-c':
            curveOutputFile = arg
        elif opt == '-d':
            dataOutputFile = arg 
    
    print("Input: %s" %inputFile)
    print("Curve: %s" %curveOutputFile)
    print("Data : %s" %dataOutputFile)

    helpers.targetCurve = lk.read(inputFile)

def main():
    getLightCurve()
    pop = runGA()
    printResults(pop)
    saveResults(pop)

if __name__ == "__main__":
    main()
    print ('finished')
