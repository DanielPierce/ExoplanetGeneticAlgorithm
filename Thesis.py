import Constants as const
from Constants import CONSTANTS

import random
#print(CONSTANTS[const.ECC][const.MAX])




def evalOneMax(individual):
    numPlanets = individual[const.NUMPLANETS]
    if(numPlanets == 0):
        return 0,
    sum = 0
    for num in range(numPlanets):
        sum += individual[num * const.ATTRPERPLANET]
    return ((sum/numPlanets) + numPlanets),

def lerp(a,b,c):
    return (c * b) + ((1-c) * a)

def randomizeAttr(currentValue, mutFactor):
    return currentValue + lerp(-1 * mutFactor, mutFactor, random.random())

def mutation(individual, indpb):
    indexToMutate = random.randint(0,len(individual)-1)
    attrIndex = indexToMutate % const.ATTRPERPLANET
    if(indexToMutate == const.NUMPLANETS):
        individual[indexToMutate] = random.randint(0, const.MAXPLANETS)
    elif (indexToMutate == const.STARRADIUS):
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARRADIUSMUTFACTOR)
        if(individual[indexToMutate] < 0):
            individual[indexToMutate] = 0
    else:
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], CONSTANTS[attrIndex][const.MUTFACTOR])
        if(individual[indexToMutate] > CONSTANTS[attrIndex][const.MAX]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MAX]
        if(individual[indexToMutate] < CONSTANTS[attrIndex][const.MIN]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MIN]
