

from cmath import isnan
import math
from numpy import NaN
import orbital
import random
import PSEUGA.common.Constants as const
from PSEUGA.common.Functs import Clamp, Lerp


class Planet:
    def __init__(self, rad=0, ec=0, sm=0, inclination=0, lo=0, ao=0, meanAnomaly=0):
        self.radius = rad # Planet's radius
        self.ecc = ec # Orbital eccentricity
        self.ecc2 = ec * ec # eccentricity squared, used for some calculations
        self.sma = sm # Orbital semi-major axis
        self.inc = inclination # orbital inclination relative to xy plane
        self.loan = lo # Orbital longitude of ascending node in xy plane relative to reference direction
        self.aop = ao # Orbital argument of periapse, ie angle from ascending node
        self.ma = meanAnomaly # Planetary mean anomaly, ie starting point at beginning of simulation where 0 is (periapse? test this -- 80% sure)
        self.period = NaN
        self.muConst = NaN


    def CalculatePeriod(self,mu):
        smaInMeters = math.pow(self.sma * 1000, 3)
        self.period = 2 * math.pi * math.sqrt(smaInMeters / mu) # m^3 / (m^3 kg^2 s^-2) = s^2, root(s^2) = s
        self.muConst = math.sqrt(mu / smaInMeters)


    def CalculateCartesianCoordsAtTime(self, secondsFromEpoch):
        unNomralMeanAnomalyAtEpoch = self.ma + secondsFromEpoch * self.muConst
        meanAnomalyAtEpoch = unNomralMeanAnomalyAtEpoch % (2 * math.pi)

        eccentricAnomalyAtEpoch = orbital.utilities.eccentric_anomaly_from_mean(self.ecc, meanAnomalyAtEpoch)
        trueAnomalyAtEpoch = orbital.utilities.true_anomaly_from_eccentric(self.ecc, eccentricAnomalyAtEpoch)

        if math.isnan(self.period):
            raise Exception("period has not yet been calculated")

        distanceFromCentralBody = self.sma * (1 - self.ecc * math.cos(eccentricAnomalyAtEpoch))

        posX = distanceFromCentralBody * math.cos(trueAnomalyAtEpoch)
        posY = distanceFromCentralBody * math.sin(trueAnomalyAtEpoch)

        coscos = math.cos(self.aop) * math.cos(self.loan)
        cossin = math.cos(self.aop) * math.sin(self.loan)
        sincos = math.sin(self.aop) * math.cos(self.loan)
        sinsin = math.sin(self.aop) * math.sin(self.loan)

        posXDot = posX * (coscos - sinsin * math.cos(self.inc)) - posY * (sincos + cossin * math.cos(self.inc))
        posYDot = posX * (cossin + sincos * math.cos(self.inc)) + posY * (coscos * math.cos(self.inc) - sinsin)
        posZDot = posX * math.sin(self.aop) * math.sin(self.inc) + posY * math.cos(self.aop) * math.sin(self.inc)

        return [posXDot, posYDot, posZDot]

    def findIntersectionWithXYPlaneIteratively(self, exp = 5):
        smallestPositive = math.inf
        smallestPositiveAngle = 0
        biggestNegative = -math.inf
        biggestNegativeAngle = 0

        smallestRecordedAbsDiff = math.inf
        smallestRecordedAbsDiffAngle = None
        smallestRecordedAbsDiffPos = None
        smallestRecordedAbsDiffEpoch = None
        #gets me the quadrant the intersection is in
        allXPos = []
        for i in range(4):
            angleToSearch = (i * math.pi) / 2
            epoch = (angleToSearch / (2 * math.pi)) * self.period
            positionAtEpoch = self.CalculateCartesianCoordsAtTime(epoch)
            xcoord = positionAtEpoch[0]

            if xcoord > 0 and xcoord < smallestPositive:
                smallestPositive = xcoord
                smallestPositiveAngle = angleToSearch
            if xcoord < 0 and xcoord > biggestNegative:
                biggestNegative = xcoord
                biggestNegativeAngle = angleToSearch
            if abs(xcoord) < smallestRecordedAbsDiff:
                smallestRecordedAbsDiff = abs(xcoord) 
                smallestRecordedAbsDiffAngle = angleToSearch
                smallestRecordedAbsDiffPos = positionAtEpoch
                smallestRecordedAbsDiffEpoch = epoch
            allXPos = allXPos + [xcoord]

        #gets me the 1/2^expth of circle, subtract two since already have quadrant (1/2^2) from above
        #possibly want to redo this so that it goes until distance is less than radius of planet?
        for i in range(exp-2):
            angleToSearch = (smallestPositiveAngle + biggestNegativeAngle) / 2
            epoch = (angleToSearch / (2 * math.pi)) * self.period
            positionAtEpoch = self.CalculateCartesianCoordsAtTime(epoch)
            xcoord = positionAtEpoch[0]
            #if abs(xcoord)  < 0.0000000001:
            #    return epoch
            if xcoord > 0 and xcoord < smallestPositive:
                smallestPositive = xcoord
                smallestPositiveAngle = angleToSearch
            if xcoord < 0 and xcoord > biggestNegative:
                biggestNegative = xcoord
                biggestNegativeAngle = angleToSearch
            if abs(xcoord) < smallestRecordedAbsDiff:
                smallestRecordedAbsDiff = abs(xcoord) 
                smallestRecordedAbsDiffAngle = angleToSearch
                smallestRecordedAbsDiffPos = positionAtEpoch
                smallestRecordedAbsDiffEpoch = epoch

        return smallestRecordedAbsDiffEpoch, smallestRecordedAbsDiffPos


    def CalculateCartesianPosition(self, trueAnomaly, eccentricAnomalyAtPosition):
        distanceFromCentralBody = self.sma * (1 - self.ecc * math.cos(eccentricAnomalyAtPosition))

        posX = distanceFromCentralBody * math.cos(trueAnomaly)
        posY = distanceFromCentralBody * math.sin(trueAnomaly)

        coscos = math.cos(self.aop) * math.cos(self.loan)
        cossin = math.cos(self.aop) * math.sin(self.loan)
        sincos = math.sin(self.aop) * math.cos(self.loan)
        sinsin = math.sin(self.aop) * math.sin(self.loan)

        posXDot = posX * (coscos - sinsin * math.cos(self.inc)) - posY * (sincos + cossin * math.cos(self.inc))
        posYDot = posX * (cossin + sincos * math.cos(self.inc)) + posY * (coscos * math.cos(self.inc) - sinsin)
        posZDot = posX * math.sin(self.aop) * math.sin(self.inc) + posY * math.cos(self.aop) * math.sin(self.inc)

        return [posXDot, posYDot, posZDot]

    def CalculateAngularSize(self, trueAnomaly):
        # how the heck does this have anything to do with angular size? There is nothing related to distance??
        return self.sma * (1 - self.ecc2) / (1 + self.ecc * math.cos(trueAnomaly))     # km?   takes degrees from periapse
    
    def CalculateCurrentTrueAnomaly(self, secondsFromEpoch):
        eccentricAnomaly = self.CalculateCurrentEccentricAnomaly(secondsFromEpoch)
        return self.CalculateTrueAnomalyFromEccentric(eccentricAnomaly)

    def CalculateCurrentEccentricAnomaly(self, secondsFromEpoch):
        epochOffsetRadians = secondsFromEpoch / self.period * 2 * math.pi # s / s * radians = radians

        currentMeanAnomaly = (self.ma + epochOffsetRadians + 2 * math.pi) % (2 * math.pi)

        try:
            eccentricAnomaly = orbital.utilities.eccentric_anomaly_from_mean(self.ecc, currentMeanAnomaly) # radians
        except:
            return -1
                    
        eccentricAnomaly = (eccentricAnomaly + 2 * math.pi) % (2 * math.pi)
        return eccentricAnomaly

    def CalculateTrueAnomalyFromEccentric(self, eccentricAnomaly):        
        trueAnomaly = 2 * math.atan(math.sqrt( (1 + self.ecc)/(1 - self.ecc) ) * math.tan(eccentricAnomaly / 2)) # radians
        trueAnomaly = (trueAnomaly + 2 * math.pi) % (2 * math.pi)
        return trueAnomaly

    def PrettyPrint(self):
        print(f"Radius: {self.radius} km")
        print(f"Eccentricity: {self.ecc}")
        print(f"Semimajor Axis: {self.sma} km")
        print(f"Inclination: {self.inc} rads")
        print(f"Longitude of Ascending Node: {self.loan} rads")
        print(f"Argument of Periapse: {self.aop} rads")
        print(f"Mean Anomaly: {self.ma} rads")
        if(self.period != 0):
            print(f"Period: {self.period} s")

    def Randomize(self):
        #This feels like it could be cleaned up
        self.radius = random.uniform(const.CONSTANTS[const.RADIUS][const.MIN], const.CONSTANTS[const.RADIUS][const.MAX])
        self.sma = random.uniform(const.CONSTANTS[const.SMA][const.MIN], const.CONSTANTS[const.SMA][const.MAX])
        #self.ecc = random.uniform(const.CONSTANTS[const.ECC][const.MIN], const.CONSTANTS[const.ECC][const.MAX])
        self.ecc = random.uniform(0.1, 0.4)
        #self.inc = random.uniform(const.CONSTANTS[const.INC][const.MIN], const.CONSTANTS[const.INC][const.MAX])
        self.inc = 0
        self.loan = random.uniform(const.CONSTANTS[const.LOAN][const.MIN], const.CONSTANTS[const.LOAN][const.MAX])
        self.aop = random.uniform(const.CONSTANTS[const.AOP][const.MIN], const.CONSTANTS[const.AOP][const.MAX])
        self.ma = random.uniform(const.CONSTANTS[const.MA][const.MIN], const.CONSTANTS[const.MA][const.MAX])

    def Mutate(self):
        numTraitsToMutate = 1
        mutationThreshold = 1 - (numTraitsToMutate / 7)

        if random.random() > mutationThreshold:
            self.radius = self.MutateAttribute(self.radius, const.CONSTANTS[const.RADIUS][const.MUTFACTOR], const.CONSTANTS[const.RADIUS][const.MIN], const.CONSTANTS[const.RADIUS][const.MAX])
        if random.random() > mutationThreshold:
            self.sma = self.MutateAttribute(self.sma, const.CONSTANTS[const.SMA][const.MUTFACTOR], const.CONSTANTS[const.SMA][const.MIN], const.CONSTANTS[const.SMA][const.MAX])
        if random.random() > mutationThreshold:
            self.ecc = self.MutateAttribute(self.ecc, const.CONSTANTS[const.ECC][const.MUTFACTOR], const.CONSTANTS[const.ECC][const.MIN], const.CONSTANTS[const.ECC][const.MAX])
        if random.random() > mutationThreshold:
            self.inc = self.MutateAttribute(self.inc, const.CONSTANTS[const.INC][const.MUTFACTOR], const.CONSTANTS[const.INC][const.MIN], const.CONSTANTS[const.INC][const.MAX])
        if random.random() > mutationThreshold:
            self.loan = self.MutateAttribute(self.loan, const.CONSTANTS[const.LOAN][const.MUTFACTOR], const.CONSTANTS[const.LOAN][const.MIN], const.CONSTANTS[const.LOAN][const.MAX])
        if random.random() > mutationThreshold:
            self.aop = self.MutateAttribute(self.aop, const.CONSTANTS[const.AOP][const.MUTFACTOR], const.CONSTANTS[const.AOP][const.MIN], const.CONSTANTS[const.AOP][const.MAX])
        if random.random() > mutationThreshold:
            self.ma = self.MutateAttribute(self.ma, const.CONSTANTS[const.MA][const.MUTFACTOR], const.CONSTANTS[const.MA][const.MIN], const.CONSTANTS[const.MA][const.MAX])

    def MutateAttribute(self, attr, mutFactor, lowerBound=-float('inf'), upperBound=float('inf')):
            attr += Lerp(random.random(), -1 * mutFactor, mutFactor)
            attr = Clamp(self.radius, lowerBound, upperBound)
            return attr

