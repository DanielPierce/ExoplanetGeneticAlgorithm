

import math
import orbital


class Planet:
    def __init__(self, rad, ec, sm, inclination, lo, ao, meanAnomaly):
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

