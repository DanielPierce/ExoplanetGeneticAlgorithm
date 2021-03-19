import Constants as const
from Constants import CONSTANTS

import lightkurve as lk
import numpy as np
import random
import statistics


targetStar = 'TIC 307210830 c'
targetCurve = lk.LightCurve()


def generateLightCurve(individual):
    myTimes = targetCurve.time
    myFlux = [random.randrange(21200,21600,1)  * targetCurve.flux.unit for i in range(len(targetCurve.time))]
    myErr = [0 * targetCurve.flux_err.unit for i in range(len(myTimes))]
    return lk.LightCurve(time=myTimes, flux=myFlux, flux_err=myErr)

def evalOneMax(individual):
    myLightCurve = generateLightCurve(individual)

    diffs = [starFlux - calcFlux for starFlux, calcFlux in zip(targetCurve.flux.value, myLightCurve.flux.value)]
    absDiff = [abs(i) for i in diffs]
    return [statistics.mean(absDiff)]

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
