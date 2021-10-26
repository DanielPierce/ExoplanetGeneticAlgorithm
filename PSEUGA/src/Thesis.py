from PSEUGA.common.TimeStep import TimeStep
from PSEUGA.common.CustomLightcurve import CustomLightcurve
from PSEUGA.common.PlanetarySystem import PlanetarySystem
from PSEUGA.common.Planet import Planet
import PSEUGA.common.Constants as const
from PSEUGA.common.Constants import ATTRSPERMUTATION, CONSTANTS

import lightkurve as lk
import numpy as np
import random
import statistics
import math

import dateutil.parser as dateparser
import time

import mpmath as mpm
import copy

from PSEUGA.src.IOHandlers import InputHandler
from PSEUGA.common.Functs import Clamp

targetStar = 'TIC 307210830 c'
targetCurve = lk.LightCurve()


def generateRandomLightCurve(individual):
    myTimes = targetCurve.time
    myFlux = [random.randrange(21200,21600,1)  * targetCurve.flux.unit for i in range(len(targetCurve.time))]
    myErr = [0 * targetCurve.flux_err.unit for i in range(len(myTimes))]
    return lk.LightCurve(time=myTimes, flux=myFlux, flux_err=myErr)

def uniformSourceResultAlgorithm(d, rp, rstar, z, z2, p, p2):
    if(1 + p < z):
        return 0
    elif(abs(1 - p) < z and z <= 1 + p):
        result1 = (1 - p2 + z2) / (2 * z)
        k1 = math.acos(result1)
        result2 = (p2 + z2 - 1) / (2 * p * z)
        k0 = math.acos(result2)
        return (1 / math.pi) * ( p2* k0 + k1 - math.sqrt((4 * z2 - math.pow(1 + z2 - p2,2)) / 4))
    elif(z <= 1-p):
        return p2
    elif(z <= p-1):
        return 1
    else:
        txt = "uniformSourceResultAlgorithm found problematic function return with:\nd:     {dval}\nrp:   {rpval}\netc"
        raise ValueError(txt.format(dval = d, rpval = rp))

def uniformSourceLightcurveAlgorithm(thisSystem):
    #overallFlux = [baseFlux for i in range(len(targetCurve.time))]
    notInRangeSkips = 0
    inRange = 0
    planetCurves = []
    input = InputHandler.getInstance()
    for planetIndex in range(thisSystem.numActivePlanets):
        
        thisPlanet = thisSystem.GetPlanet(planetIndex)
        customCurve = CustomLightcurve(targetCurve)
        customCurve.setBaseFlux(thisSystem.star.flux)
        #rstar = thisSystem.CalculateStarAngularSize()
                                                                                                              
        thisSystem.CalculatePlanetaryPeriod(planetIndex)
        
        #zeroTime = customCurve.epochTime
        #myFlux = [thisSystem.star.flux for i in range(len(targetCurve.time))]
        #baseStepIndices = list(range(0, len(targetCurve.time), const.skippedTimesteps))
        followUp = []
        baseCustomSteps = customCurve.getUnskippedTimesteps(input.runSettings['timestepsToSkip'])
        for currentStepIndex in baseCustomSteps:
            currentStep = baseCustomSteps[currentStepIndex]
            #myFlux[i] = calculateTimestepFromTime(i, zeroTime, thisPlanet, thisSystem.star)
            currentStep.flux = calculateTimestep(currentStep, thisPlanet, thisSystem.star)
            # if not baseflux, add all timeslots between this and previous to followup
            # also add timeslots between this and next
            # once done with this for, remove duplicates from follup and then calculate all those
            if(currentStep.flux != thisSystem.star.flux):
                inRange = inRange + 1
                for x in range(currentStepIndex - input.runSettings['timestepsToSkip'], currentStepIndex):
                    followUp.append(customCurve.timeSteps[x])
                for x in range(currentStepIndex,currentStepIndex + input.runSettings['timestepsToSkip']):
                    followUp.append(customCurve.timeSteps[x])
            else:
                notInRangeSkips = notInRangeSkips + 1
        followUpUnique = list(dict.fromkeys(followUp))
        for currentStep in followUpUnique:
            currentStep.flux = calculateTimestep(currentStep, thisPlanet, thisSystem.star)
        planetCurves.append(customCurve)
    
    targetCustom = CustomLightcurve(targetCurve)
    finalCurve = CustomLightcurve()
    finalCurve.epochTime = targetCustom.epochTime
    for i in range(len(targetCustom.timeSteps)):
        if(thisSystem.numActivePlanets == 0):
            continue
        try:
            if(i >= len(planetCurves[0].timeSteps)):
                continue
        except IndexError:
            print("Hit index error, dumping current data")
            print(f"Size of system's planets array: {len(thisSystem.planets)}")
            print(f"Size of planetCurves array: {len(planetCurves)}")
            thisSystem.PrettyPrint()
            exit()
        #print(f"{i} {len(planetCurves)} {len(planetCurves[0].timeSteps)}")
        minFlux = planetCurves[0].timeSteps[i].flux
        for j in range(1, len(planetCurves)):
            if(planetCurves[j].timeSteps[i].flux < minFlux):
                minFlux = planetCurves[j].timeSteps[i].flux
        #minFlux = min(planetCurves, lambda c: c.timeSteps[i].flux)
        newStep = TimeStep(targetCustom.timeSteps[i].secondsFromEpoch, minFlux)
        finalCurve.timeSteps.append(newStep)
    return finalCurve


def calculateTimestep(timestep, thisPlanet, thisStar):
    deltaTime = timestep.secondsFromEpoch
    trueAnomaly = thisPlanet.CalculateCurrentTrueAnomaly(deltaTime)

    currentRadius = thisPlanet.CalculateAngularSize(trueAnomaly)

    cartesianPosition = thisPlanet.CalculateCartesianPosition(trueAnomaly)
    if(cartesianPosition[1] < 0):
        return thisStar.flux

    rstar = thisStar.radius            
    rp = getAngularSizeFromSizeAndDist(currentRadius, 2 * thisStar.distanceTo + cartesianPosition[1])
            
    p = rp / rstar
    p2 = math.pow(p,2)

    collapser = [1,0,1]
    collapsedPosition = [a * b for a,b in zip(cartesianPosition, collapser)]

    # center to center distance between star and planet
    cartesianDist = math.dist([0,0,0], collapsedPosition)

    d = getAngularSizeFromSizeAndDist(cartesianDist, thisStar.distanceTo)

    if(d > rp + rstar):
        return thisStar.flux

    z = d / rstar
    z2 = math.pow(z,2)

    return thisStar.flux * (1 - uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2))

def getAngularSizeFromSizeAndDist(size, distance):
    angularSize = mpm.mpf('0')
    trueSize = mpm.mpf(size)
    trueDist = mpm.mpf(distance)
    trueRatio = mpm.mpf(trueSize / trueDist)

    return mpm.atan(trueRatio)

def evalOneMaxDist(individual):
    generatedCustomCurve = uniformSourceLightcurveAlgorithm(individual.ps)
    targetCustom = CustomLightcurve(targetCurve)
    targetSorted = targetCustom.sortByFlux()
    generatedSorted = generatedCustomCurve.sortByFlux()
    sumOfDists = 0
    for i in range(len(generatedSorted.timeSteps)):
        sumOfDists += targetSorted.timeSteps[i].distanceToTimestep(generatedSorted.timeSteps[i])
    individual.lc = generatedCustomCurve
    return [sumOfDists]



def mutation(individual):
    mutationThreshold = 0.9

    for index in range(const.MAXPLANETS + 2):
        if(random.random() > mutationThreshold):
            if index == 21:
                individual.ps.numActivePlanets += random.randint(-1 * const.NUMPLANETSMUTFACTOR, const.NUMPLANETSMUTFACTOR)
                individual.ps.numActivePlanets = Clamp(individual.ps.numActivePlanets, 0, const.MAXPLANETS)
            elif index == 20:
                individual.ps.star.Mutate()
            else:
                individual.ps.planets[index].Mutate()

def randomizeIndividual(individual):
    index = random.randint(0,21)
    if index == 21:
        individual.ps.numActivePlanets = random.randint(0,20)
    elif index == 20:
        individual.ps.star.Randomize()
    else:
        individual.ps.planets[index].Randomize()
    return individual

