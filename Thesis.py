import Constants as const
from Constants import CONSTANTS

import lightkurve as lk
import numpy as np
import random
import statistics
import math
import orbital

import dateutil.parser as dateparser
import time

targetStar = 'TIC 307210830 c'
targetCurve = lk.LightCurve()


def generateRandomLightCurve(individual):
    myTimes = targetCurve.time
    myFlux = [random.randrange(21200,21600,1)  * targetCurve.flux.unit for i in range(len(targetCurve.time))]
    myErr = [0 * targetCurve.flux_err.unit for i in range(len(myTimes))]
    return lk.LightCurve(time=myTimes, flux=myFlux, flux_err=myErr)

def uniformSourceResultAlgorithm(d, rp, rstar, z, z2, p, p2, f, k0, k1):
    if(1 + p < z):
        return 0
    elif(abs(1 - p) < z and z <= 1 + p):
        return (1 / math.pi) * ( p2* k0 + k1 - math.sqrt((4 * z2 - math.pow(1 + z2 - p2,2)) / 4))
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
    for planetIndex in range(individual[const.NUMPLANETS]):
        print("on planet %s" %planetIndex)
        rstar = individual[const.STARRADIUS] # r*, stellar radius in km
        rp = individual[const.ATTRPERPLANET * planetIndex + const.RADIUS] # rp, planetary radius
        p = rp / rstar # size ratio
                                                                                                              # semimajor axis in km to m
        a = individual[const.ATTRPERPLANET * planetIndex + const.SMA] * 1000
        mu = const.GRAVITATIONALCONSTANT * individual[const.STARMASS]
        period = 2 * math.pi * math.sqrt(math.pow(a, 3) / mu)
        
        zeroTime = dateparser.parse(targetCurve.time.iso[0])

        myFlux = []
        for timeIndex in targetCurve.time.iso:
            sma = individual[const.ATTRPERPLANET * planetIndex + const.SMA]
            ecc = individual[const.ATTRPERPLANET * planetIndex + const.ECC]
            ecc2 = math.pow(ecc,2)

            currentTime = dateparser.parse(timeIndex)
            epochOffset = (currentTime - zeroTime).total_seconds() / period * 360

            try:
                eccentricAnomaly = orbital.utilities.eccentric_anomaly_from_mean(ecc, individual[const.ATTRPERPLANET * planetIndex + const.MA] + epochOffset)
            except:
                toc = time.perf_counter()
                #print(f"Planet {planetIndex} eccentric anomaly could not converge at timestep {timeIndex} in {toc - tic:0.4f} seconds")
                myFlux.append(-1)
                continue
                
            #print(f"Planet {planetIndex} eccentric anomaly did converge at timestep {timeIndex}")
            trueAnomaly = 2 * math.atan(math.sqrt( (1 + ecc)/(1 - ecc) ) * math.tan(eccentricAnomaly / 2))

            #d = (sma * (1 - ecc2))/(1-ecc * math.cos(deg * math.radians(trueAnomaly))) # center to center distance between star and planet
            d = (sma * (1 - ecc2))/(1-ecc * math.cos(math.radians(trueAnomaly))) # center to center distance between star and planet

            z = d / rstar # normalized seperation of centers
            p2 = math.pow(p,2)
            z2 = math.pow(z,2)
            result1 = (1 - p2 + z2) / (2 * z) # wants radians, all planet characteristics are in degrees
            k1 = math.acos(math.radians(result1))
            result2 = (p2 + z2 - 1) / (2 * p * z)
            k0 = math.acos(math.radians(result2)) # wants radians, all planet characteristics are in degrees

            myFlux.append(uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2,baseFlux,k0,k1))
            #print(f"Did timestep in {toc - tic:0.4f} seconds")
        for i in range(len(targetCurve.time)):
            overallFlux[i] = min(overallFlux[i], myFlux[i])
    numRejects = overallFlux.count(-1)
    steps = len(overallFlux)
    toc = time.perf_counter()
    print(f"Lightcurve had {numRejects} convergence errors out of {steps} timesteps in {toc - tic:0.4f} seconds")
    return overallFlux


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
    print(f"In eval, counted {numCounted}")
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
        if(individual[indexToMutate] <= 10000):
            individual[indexToMutate] = 10000
    elif (indexToMutate == const.STARMASS):
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARMASSMUTFACTOR)
        if(individual[indexToMutate] <= const.STARMASSMIN):
            individual[indexToMutate] = const.STARMASSMIN
    elif (indexToMutate == const.STARBASEFLUX):
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], const.STARBASEFLUXMUTFACTOR)
        if(individual[indexToMutate] < 0):
            individual[indexToMutate] = 0
    else:
        individual[indexToMutate] = randomizeAttr(individual[indexToMutate], CONSTANTS[attrIndex][const.MUTFACTOR])
        if(individual[indexToMutate] > CONSTANTS[attrIndex][const.MAX]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MAX]
        if(individual[indexToMutate] < CONSTANTS[attrIndex][const.MIN]):
            individual[indexToMutate] = CONSTANTS[attrIndex][const.MIN]

def validateIndividual(individual):
    for index in range(len(individual)):
        attrIndex = index % const.ATTRPERPLANET
        if (index == const.NUMPLANETS):
            individual[index] = random.randint(0, const.MAXPLANETS)
        elif (index == const.STARRADIUS):
            if(individual[index] <= const.STARRADIUSMIN):
                individual[index] = const.STARRADIUSMIN
        elif (index == const.STARMASS):
            if(individual[index] <= const.STARMASSMIN):
                individual[index] = const.STARMASSMIN
        elif (index == const.STARBASEFLUX):
            if(individual[index] < const.STARBASEFLUXMIN):
                individual[index] = const.STARBASEFLUXMIN
        elif (attrIndex == const.RADIUS):
            individual[index] = random.randint(CONSTANTS[attrIndex][const.MIN], CONSTANTS[attrIndex][const.MIN])
        elif (attrIndex == const.SMA):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.ECC):
            individual[index] = random.uniform(0.1, 0.4)
        elif (attrIndex == const.INC):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.LOAN):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.AOP):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
        elif (attrIndex == const.MA):
            if(individual[index] > CONSTANTS[attrIndex][const.MAX]):
                individual[index] = CONSTANTS[attrIndex][const.MAX]
            if(individual[index] < CONSTANTS[attrIndex][const.MIN]):
                individual[index] = CONSTANTS[attrIndex][const.MIN]
            
        
    return individual

