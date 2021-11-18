
import random
import PSEUGA.common.Constants as const
from PSEUGA.common.Planet import Planet
from PSEUGA.common.Star import Star
import mpmath as mpm

class PlanetarySystem:
    def __init__(self, individual=None):
        self.planets = []
        if individual is not None:
            for i in range(const.MAXPLANETS):
                self.planets.append(Planet(
                individual[const.ATTRPERPLANET * i + const.RADIUS],
                individual[const.ATTRPERPLANET * i + const.ECC],
                individual[const.ATTRPERPLANET * i + const.SMA],
                individual[const.ATTRPERPLANET * i + const.INC],
                individual[const.ATTRPERPLANET * i + const.LOAN],
                individual[const.ATTRPERPLANET * i + const.AOP],
                individual[const.ATTRPERPLANET * i + const.MA]))
            self.distanceTo = individual[const.DISTANCE]
            self.numActivePlanets = individual[const.NUMPLANETS]
            self.star = Star(individual[const.STARRADIUS], individual[const.STARMASS], individual[const.STARBASEFLUX], self.distanceTo)
        else:
            for i in range(const.MAXPLANETS):
                self.planets.append(Planet())
            self.distanceTo = 0
            self.numActivePlanets = 0
            self.star = Star(0,0,0,0)    

    def GetPlanet(self, index):
        return self.planets[index]

    def CalculateStarAngularSize(self):
        rstar = mpm.mpf('0')
        starSize = mpm.mpf(self.star.radius)
        starDist = mpm.mpf(2 * self.distanceTo)
        starRatio = mpm.mpf(starSize / starDist)

        return mpm.atan(starRatio)

    def CalculatePlanetaryPeriod(self, planetIndex):
        thisPlanet = self.planets[planetIndex]
        a = thisPlanet.sma * 1000 # semimajor axis in km to m
        mu = const.GRAVITATIONALCONSTANT * self.star.mass # m^3 kg^-1 s^-2 * kg = m^3 s^-2
        thisPlanet.CalculatePeriod(mu)

    def PrettyPrint(self):
        if self.numActivePlanets > 0:
            for i in range(self.numActivePlanets):
                print(f"--- Planet {i} ---")
                self.planets[i].PrettyPrint()
        else:
            print("No active planets")
        print("--- Star ---")
        self.star.PrettyPrint()

    def ToList(self):
        retList = []
        for index in range(const.MAXPLANETS):
            retList.append(self.planets[index].radius)
            retList.append(self.planets[index].sma)
            retList.append(self.planets[index].ecc)
            retList.append(self.planets[index].inc)
            retList.append(self.planets[index].loan)
            retList.append(self.planets[index].aop)
            retList.append(self.planets[index].ma)
        retList.append(self.star.radius)
        retList.append(self.star.mass)
        retList.append(self.star.flux)
        retList.append(self.star.distanceTo)
        retList.append(self.numActivePlanets)
        return retList

