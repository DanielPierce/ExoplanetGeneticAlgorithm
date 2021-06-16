from Planet import Planet
import Constants as const
from Constants import CONSTANTS

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
    tic = time.perf_counter()
    baseFlux = individual[const.STARBASEFLUX]
    overallFlux = [baseFlux for i in range(len(targetCurve.time))]
    notInRangeSkips = 0
    inRange = 0
    for planetIndex in range(individual[const.NUMPLANETS]):
        print("on planet %s" %planetIndex)
        
        thisPlanet = Planet(
            individual[const.ATTRPERPLANET * planetIndex + const.RADIUS],
            individual[const.ATTRPERPLANET * planetIndex + const.ECC],
            individual[const.ATTRPERPLANET * planetIndex + const.SMA],
            individual[const.ATTRPERPLANET * planetIndex + const.INC],
            individual[const.ATTRPERPLANET * planetIndex + const.LOAN],
            individual[const.ATTRPERPLANET * planetIndex + const.AOP],
            individual[const.ATTRPERPLANET * planetIndex + const.MA])
        #rstar = individual[const.STARRADIUS] # r*, stellar radius in km
        #rp = individual[const.ATTRPERPLANET * planetIndex + const.RADIUS] # rp, planetary radius in km
        rstar = mpm.mpf('0')
        starSize = mpm.mpf(individual[const.STARRADIUS])
        starDist = mpm.mpf(2 * individual[const.DISTANCE])
        starRatio = mpm.mpf(starSize / starDist)

        rstar = mpm.atan(starRatio)
                                                                                                              
        a = thisPlanet.sma * 1000 # semimajor axis in km to m
        mu = const.GRAVITATIONALCONSTANT * individual[const.STARMASS] # m^3 kg^-1 s^-2 * kg = m^3 s^-2
        thisPlanet.CalculatePeriod(mu)
        
        zeroTime = dateparser.parse(targetCurve.time.iso[0])

        myFlux = [const.STARBASEFLUX for i in range(len(targetCurve.time))]
        timelen = len(targetCurve.time)
        baseSteps = range(0, len(targetCurve.time),const.skippedTimesteps)
        baseStepsLen = len(baseSteps)
        baseTimeSteps = [i for i in range(baseStepsLen)]
        followUp = []

        for i in baseTimeSteps:
            myFlux[i] = calculateTimestep(i, zeroTime, thisPlanet, rstar, baseFlux, individual[const.DISTANCE])
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
            #print(f"Did timestep in {toc - tic:0.4f} seconds")
        followUpUnique = list(dict.fromkeys(followUp))
        for i in followUpUnique:
            myFlux[i] = calculateTimestep(i, zeroTime, thisPlanet, rstar, baseFlux, individual[const.DISTANCE])
        followUps = len(followUpUnique)
        leng = len(baseTimeSteps)
        print(f"needed to follow up on {followUps} timesteps, had {notInRangeSkips} skips and {inRange} in range out of {leng} timesteps from {timelen} curvelength w stepslen {baseStepsLen} skipping {const.skippedTimesteps} timesteps per")
        for i in range(len(targetCurve.time)):
            overallFlux[i] = min(overallFlux[i], myFlux[i])
    numRejects = overallFlux.count(-1)
    steps = len(overallFlux)
    toc = time.perf_counter()
    print(f"Lightcurve had {numRejects} convergence errors out of {steps} timesteps in {toc - tic:0.4f} seconds")
    return overallFlux

def calculateTimestep(i, zeroTime, thisPlanet, rstar, baseFlux, distance):
    timeIndex = targetCurve.time.iso[i]
    currentTime = dateparser.parse(timeIndex)
    seconds = (currentTime - zeroTime).total_seconds()

    trueAnomaly = thisPlanet.CalculateCurrentTrueAnomaly(seconds)

    currentRadius = thisPlanet.CalculateAngularSize(trueAnomaly)

    cartesianPosition = thisPlanet.CalculateCartesianPosition(trueAnomaly)
    if(cartesianPosition[1] < 0):
        return baseFlux
            
    rp = getAngularSizeFromSizeAndDist(thisPlanet.radius, 2 * distance + cartesianPosition[1])
            
    p = rp / rstar # size ratio, km/km = scalar
    p2 = math.pow(p,2) # scalar squared = scalar

    collapser = [1,0,1]
    collapsedPosition = [a * b for a,b in zip(cartesianPosition, collapser)]

    # center to center distance between star and planet
    cartesianDist = math.dist([0,0,0], collapsedPosition)

    d = getAngularSizeFromSizeAndDist(cartesianDist, distance)

    if(d > rp + rstar):
        return baseFlux

    z = d / rstar
    z2 = math.pow(z,2)

    return (uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2,baseFlux))

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
    return sumOfDists

def lerp(a,b,c):
    return (c * b) + ((1-c) * a)

def randomizeAttr(currentValue, mutFactor):
    return currentValue + lerp(-1 * mutFactor, mutFactor, random.random())

def mutation(individual, indpb):
    indexToMutate = random.randint(0,len(individual)-1)
    attrIndex = indexToMutate % const.ATTRPERPLANET
    if(indexToMutate == const.NUMPLANETS):
        print("mutating num planets")
        individual[indexToMutate] = random.randint(0, const.MAXPLANETS)
    elif (indexToMutate == const.STARRADIUS):
        print("mutating star radius")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARRADIUSMUTFACTOR)
        if(individual[indexToMutate] <= 10000):
            individual[indexToMutate] = 10000
    elif (indexToMutate == const.STARMASS):
        print("mutating star mass")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARMASSMUTFACTOR)
        if(individual[indexToMutate] <= const.STARMASSMIN):
            individual[indexToMutate] = const.STARMASSMIN
    elif (indexToMutate == const.STARBASEFLUX):
        print("mutating base flux")
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARBASEFLUXMUTFACTOR)
        if(individual[indexToMutate] < 0):
            individual[indexToMutate] = 0
    else:
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], CONSTANTS[attrIndex][const.MUTFACTOR])
        print(f"mutating attribute {indexToMutate}")
        if(individual[indexToMutate] > CONSTANTS[attrIndex][const.MAX]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MAX]
        if(individual[indexToMutate] < CONSTANTS[attrIndex][const.MIN]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MIN]

def validateIndividual(individual):
    for index in range(len(individual)):
        attrIndex = index % const.ATTRPERPLANET
        if (index == const.NUMPLANETS):
            #individual[index] = random.randint(0, const.MAXPLANETS)
            individual[index] = 3
        elif (index == const.STARRADIUS):
            individual[index] = const.STARRADIUSMIN * 10000000
        elif (index == const.STARMASS):
            if(individual[index] <= const.STARMASSMIN):
                individual[index] = const.STARMASSMIN
        elif (index == const.STARBASEFLUX):
            if(individual[index] < const.STARBASEFLUXMIN):
                individual[index] = const.STARBASEFLUXMIN
        elif (index == const.DISTANCE):
            individual[index] = 9460730000000 * 10 # 10 lightyears
        elif (attrIndex == const.RADIUS):
            individual[index] = CONSTANTS[attrIndex][const.MAX]
        elif (attrIndex == const.SMA):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.ECC):
            #individual[index] = random.uniform(0.1, 0.4)
            individual[index] = 0
        elif (attrIndex == const.INC):
            individual[index] = 0
        elif (attrIndex == const.LOAN):
            if(individual[index] >= CONSTANTS[attrIndex][const.MAX]):
                individual[index] = individual[index] % CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.AOP):
            if(individual[index] >= CONSTANTS[attrIndex][const.MAX]):
                individual[index] = individual[index] % CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.MA):
            if(individual[index] >= CONSTANTS[attrIndex][const.MAX]):
                individual[index] = individual[index] % CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]    
    return individual

def randomizeIndividual(individual):
    for index in range(len(individual)):
        attrIndex = index % const.ATTRPERPLANET
        if (index == const.NUMPLANETS):
            individual[index] = random.randint(0, const.MAXPLANETS)
            #individual[index] = 3
        elif (index == const.STARRADIUS):
            individual[index] = random.randint(const.STARRADIUSMIN, const.STARRADIUSMIN * 100)
        elif (index == const.STARMASS):
            individual[index] = random.randint(const.STARMASSMIN, const.STARMASSMIN * 100)
        elif (index == const.STARBASEFLUX):
            individual[index] = random.randint(const.STARBASEFLUXMIN, const.STARBASEFLUXMIN * 100)
        elif (attrIndex == const.RADIUS):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN], CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.SMA):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.ECC):
            individual[index] = random.uniform(0.1, 0.4)
        elif (attrIndex == const.INC):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.LOAN):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.AOP):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
        elif (attrIndex == const.MA):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN],CONSTANTS[attrIndex][const.MAX])
    return individual

