from PlanetarySystem import PlanetarySystem
from Planet import Planet
import Constants as const
from Constants import ATTRSPERMUTATION, CONSTANTS

import lightkurve as lk
import numpy as np
import random
import statistics
import math

import dateutil.parser as dateparser
import time

import mpmath as mpm

targetStar = 'TIC 307210830 c'
targetCurve = lk.LightCurve()


def generateRandomLightCurve(individual):
    myTimes = targetCurve.time
    myFlux = [random.randrange(21200,21600,1)  * targetCurve.flux.unit for i in range(len(targetCurve.time))]
    myErr = [0 * targetCurve.flux_err.unit for i in range(len(myTimes))]
    return lk.LightCurve(time=myTimes, flux=myFlux, flux_err=myErr)

def uniformSourceResultAlgorithm(d, rp, rstar, z, z2, p, p2, f):
    if(1 + p < z):
        return 0
    elif(abs(1 - p) < z and z <= 1 + p):
        result1 = (1 - p2 + z2) / (2 * z)
        k1 = math.acos(result1)
        result2 = (p2 + z2 - 1) / (2 * p * z)
        k0 = math.acos(result2)
        return f * (1 / math.pi) * ( p2* k0 + k1 - math.sqrt((4 * z2 - math.pow(1 + z2 - p2,2)) / 4))
    elif(z <= 1-p):
        return math.pow(p,2)
    elif(z <= p-1):
        return 1
    else:
        txt = "uniformSourceResultAlgorithm found problematic function return with:\nd:     {dval}\nrp:   {rpval}\netc"
        raise ValueError(txt.format(dval = d, rpval = rp))

def uniformSourceLightcurveAlgorithm(individual):
    baseFlux = individual[const.STARBASEFLUX]
    overallFlux = [baseFlux for i in range(len(targetCurve.time))]
    notInRangeSkips = 0
    inRange = 0
    thisSystem = PlanetarySystem(individual)
    for planetIndex in range(thisSystem.numActivePlanets):
        
        thisPlanet = thisSystem.GetPlanet(planetIndex)
       
        rstar = thisSystem.CalculateStarAngularSize()
                                                                                                              
        thisSystem.CalculatePlanetaryPeriod(planetIndex)
        
        zeroTime = dateparser.parse(targetCurve.time.iso[0])
        myFlux = [thisSystem.star.flux for i in range(len(targetCurve.time))]
        baseStepIndices = list(range(0, len(targetCurve.time), const.skippedTimesteps))
        followUp = []

        for i in baseStepIndices:
            myFlux[i] = calculateTimestep(i, zeroTime, thisPlanet, thisSystem.star)
            # if not baseflux, add all timeslots between this and previous to followup
            # also add timeslots between this and next
            # once done with this for, remove duplicates from follup and then calculate all those
            if(myFlux[i] != baseFlux):
                inRange = inRange + 1
                for x in range(i - const.skippedTimesteps, i):
                    followUp.append(x)
                for x in range(i, i + const.skippedTimesteps):
                    followUp.append(x)
            else:
                notInRangeSkips = notInRangeSkips + 1
        followUpUnique = list(dict.fromkeys(followUp))
        for i in followUpUnique:
            if(i >= len(targetCurve.time)):
                continue
            myFlux[i] = calculateTimestep(i, zeroTime, thisPlanet, thisSystem.star)
        for i in range(len(targetCurve.time)):
            overallFlux[i] = min(overallFlux[i], myFlux[i])
    return overallFlux

def calculateTimestep(i, zeroTime, thisPlanet, thisStar):
    if(i >= len(targetCurve.time)):
        return thisStar.flux
    timeIndex = targetCurve.time.iso[i]
    currentTime = dateparser.parse(timeIndex)
    seconds = (currentTime - zeroTime).total_seconds()

    trueAnomaly = thisPlanet.CalculateCurrentTrueAnomaly(seconds)

    currentRadius = thisPlanet.CalculateAngularSize(trueAnomaly)

    cartesianPosition = thisPlanet.CalculateCartesianPosition(trueAnomaly)
    if(cartesianPosition[1] < 0):
        return thisStar.flux

    rstar = thisStar.radius            
    rp = getAngularSizeFromSizeAndDist(thisPlanet.radius, 2 * thisStar.distanceTo + cartesianPosition[1])
            
    p = rp / rstar # size ratio, km/km = scalar
    p2 = math.pow(p,2) # scalar squared = scalar

    collapser = [1,0,1]
    collapsedPosition = [a * b for a,b in zip(cartesianPosition, collapser)]

    # center to center distance between star and planet
    cartesianDist = math.dist([0,0,0], collapsedPosition)

    d = getAngularSizeFromSizeAndDist(cartesianDist, thisStar.distanceTo)

    if(d > rp + rstar):
        return thisStar.flux

    z = d / rstar
    z2 = math.pow(z,2)

    return (uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2, thisStar.flux))

def getAngularSizeFromSizeAndDist(size, distance):
    angularSize = mpm.mpf('0')
    trueSize = mpm.mpf(size)
    trueDist = mpm.mpf(distance)
    trueRatio = mpm.mpf(trueSize / trueDist)

    return mpm.atan(trueRatio)

def generateLightcurve(individual):
    myTimes = targetCurve.time
    #myFlux = [random.randrange(21200,21600,1)  * targetCurve.flux.unit for i in range(len(targetCurve.time))]
    myFlux = uniformSourceLightcurveAlgorithm(individual)
    if(myFlux is None):
        myFlux = [0 for i in range(len(myTimes))]
    myErr = [0 for i in range(len(myTimes))]
    return lk.LightCurve(time=myTimes, flux=myFlux, flux_err=myErr)

def evalOneMax(individual):
    myLightCurve = generateLightcurve(individual)
    if(myLightCurve is None):
        return 0
    diffs = [starFlux - calcFlux for starFlux, calcFlux in zip(targetCurve.flux.value, myLightCurve.flux.value)]
    absDiff = [abs(i) for i in diffs]
    numCounted = 0
    for i in range(len(myLightCurve.flux)):
        if(myLightCurve.flux[i] == -1):
            absDiff[i] = 10000
            numCounted += 1
    # instead of checking only in a specific timestep, order timesteps by brightness and then do a distance function to see how far in time + brightness the dimmest timestep is, then next, etc, where fitness is sum of distances and trying to minimize that distance
    print(f"In eval, counted {numCounted}")
    return [statistics.mean(absDiff)]

def sortLightcurves(myLightCurve):
    targetFlux = targetCurve.flux.value.tolist()
    targetTimes = targetCurve.time.iso.tolist()
    currentFlux = myLightCurve.flux.value.tolist()
    currentTimes = myLightCurve.time.iso.tolist()

    targetSort = sorted(zip(targetFlux,targetTimes), key=lambda i: i[0], reverse=True)
    currentSort = sorted(zip(currentFlux,currentTimes), key=lambda i: i[0], reverse=True)
    return targetSort,currentSort
    

def evalOneMaxDist(individual):
    myLightCurve = generateLightcurve(individual)
    if(myLightCurve is None):
        return 0
    targetSort, currentSort = sortLightcurves(myLightCurve)
    sumOfDists = 0
    for i in range(len(targetSort)):
        targetFlux, targetTime = targetSort[i]
        currentFlux, currentTime = currentSort[i]
        fluxDist = targetFlux - currentFlux
        timeDist = dateparser.parse(targetTime) - dateparser.parse(currentTime)
        sumOfDists = sumOfDists + math.sqrt(fluxDist * fluxDist + timeDist.total_seconds() * timeDist.total_seconds())
    return [sumOfDists]

def lerp(a,b,c):
    return (c * b) + ((1-c) * a)

def randomizeAttr(currentValue, mutFactor):
    return currentValue + lerp(-1 * mutFactor, mutFactor, random.random())

def mutateAttr(individual, indexToMutate):
    attrIndex = indexToMutate % const.ATTRPERPLANET
    if(indexToMutate == const.NUMPLANETS):
        #print("mutating num planets")
        individual[indexToMutate] = random.randint(0, const.MAXPLANETS)
    elif (indexToMutate == const.STARRADIUS):
        #print("mutating star radius")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARRADIUSMUTFACTOR)
        if(individual[indexToMutate] <= 10000):
            individual[indexToMutate] = 10000
    elif (indexToMutate == const.STARMASS):
        #print("mutating star mass")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARMASSMUTFACTOR)
        if(individual[indexToMutate] <= const.STARMASSMIN):
            individual[indexToMutate] = const.STARMASSMIN
    elif (indexToMutate == const.STARBASEFLUX):
        #print("mutating base flux")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARBASEFLUXMUTFACTOR)
        if(individual[indexToMutate] < 0):
            individual[indexToMutate] = 0
    else:
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], CONSTANTS[attrIndex][const.MUTFACTOR])
        #print(f"mutating attribute {indexToMutate}")
        if(individual[indexToMutate] > CONSTANTS[attrIndex][const.MAX]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MAX]
        if(individual[indexToMutate] < CONSTANTS[attrIndex][const.MIN]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MIN]


def mutation(individual):
    # have probability to mutate each attribute, with prob such that 2 are mutate per
    mutationThreshold = 1 - (const.ATTRSPERMUTATION / len(individual))
    for i in range(len(individual)):
        if(random.random() > mutationThreshold):
            mutateAttr(individual, i)

def randomizeIndividual(individual):
    for index in range(len(individual)):
        attrIndex = index % const.ATTRPERPLANET
        if (index == const.NUMPLANETS):
            individual[index] = random.randint(0, const.MAXPLANETS)
        elif (index == const.STARRADIUS):
            individual[index] = random.uniform(const.STARRADIUSMIN, const.STARRADIUSMIN * 100)
        elif (index == const.STARMASS):
            individual[index] = random.uniform(const.STARMASSMIN, const.STARMASSMIN * 100)
        elif (index == const.STARBASEFLUX):
            individual[index] = random.uniform(const.STARBASEFLUXMIN, const.STARBASEFLUXMIN * 100)
        elif (index == const.DISTANCE):
            individual[index] = random.uniform(const.DISTANCEMIN, const.DISTANCEMIN * 100)
        elif (attrIndex == const.RADIUS):
            individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN], CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.SMA):
            individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.ECC):
            individual[index] = random.uniform(0.1, 0.4)
        elif (attrIndex == const.INC):
            #individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
            individual[index] = 0
        elif (attrIndex == const.LOAN):
            individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.AOP):
            individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.MA):
            individual[index] = random.uniform(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
    return individual

