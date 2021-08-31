
import PSEUGA.common.Constants as const
from PSEUGA.common.Planet import Planet
from PSEUGA.common.Star import Star
import mpmath as mpm

class PlanetarySystem:
    def __init__(self, individual):
        self.planets = []
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
        

