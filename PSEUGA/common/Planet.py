

import math
import orbital
import random
import PSEUGA.common.Constants as const


class Planet:
    def __init__(self, rad=0, ec=0, sm=0, inclination=0, lo=0, ao=0, meanAnomaly=0):
        self.radius = rad
        self.ecc = ec
        self.ecc2 = ec * ec
        self.sma = sm
        self.inc = inclination
        self.loan = lo
        self.aop = ao
        self.ma = meanAnomaly
        self.period = 0


    def CalculatePeriod(self,mu):
        if self.period != 0:
            return
        self.period = 2 * math.pi * math.sqrt(math.pow(self.sma * 1000, 3) / mu) # m^3 / (m^3 kg^2 s^-2) = s^2, root(s^2) = s

    def CalculateCartesianPosition(self, trueAnomaly):
        posX = self.radius * math.cos(trueAnomaly)
        posY = self.radius * math.sin(trueAnomaly)
        posZ = 0

        coscos = math.cos(self.aop) * math.cos(self.loan)
        cossin = math.cos(self.aop) * math.sin(self.loan)
        sincos = math.sin(self.aop) * math.cos(self.loan)
        sinsin = math.sin(self.aop) * math.sin(self.loan)

        posXDot = posX * coscos - sinsin * math.cos(self.inc) - posY * sincos + cossin * math.cos(self.inc)
        posYDot = posX * cossin + sincos * math.cos(self.inc) + posY * coscos * math.cos(self.inc) - sinsin
        posZDot = posX * math.sin(self.aop) * math.sin(self.inc) + posY * math.cos(self.aop) * math.sin(self.inc)

        # intertial reference frame, but relative to what? orbital frame has z axis perpendicular to orbital plane and x axis pointing to periapsis of orbit
        # May have to rotate to meet expected reference frame
        return [posXDot, posYDot, posZDot]

    def CalculateAngularSize(self, trueAnomaly):
        return self.sma * (1 - self.ecc2) / (1 + self.ecc * math.cos(trueAnomaly))     # km?   takes degrees from periapse

    def CalculateCurrentTrueAnomaly(self, secondsFromEpoch):
        epochOffsetRadians = secondsFromEpoch / self.period * 2 * math.pi # s / s * radians = radians

        currentMeanAnomaly = (self.ma + epochOffsetRadians + 2 * math.pi) % (2 * math.pi)

        try:
            eccentricAnomaly = orbital.utilities.eccentric_anomaly_from_mean(self.ecc, currentMeanAnomaly) # radians???
        except:
            #print(f"Planet {planetIndex} eccentric anomaly could not converge at timestep {timeIndex} in {toc - tic:0.4f} seconds")
            return -1
                    
        eccentricAnomaly = (eccentricAnomaly + 2 * math.pi) % (2 * math.pi)
                
        trueAnomaly = 2 * math.atan(math.sqrt( (1 + self.ecc)/(1 - self.ecc) ) * math.tan(eccentricAnomaly / 2)) # radians???
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
            print(f"Period: {self.period} km")

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
