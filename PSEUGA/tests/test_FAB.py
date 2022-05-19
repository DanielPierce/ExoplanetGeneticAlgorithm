
from cmath import inf
import unittest
from PSEUGA.common.Planet import Planet
from orbital.utilities import *
import math
import csv
import numpy as np

class TestFAB(unittest.TestCase):

    #@unittest.skip('drawing stuff first')
    def testAngleToZYPlane(self):
        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 1,000,000,000kg
        testMu = G * 1000000000000000000000000000
        

        ecc = 0.4  # range 0-1   yes  e no
        inc = 40  # range 0-360 not w loan e no
        loan = 50 # range 0-360 yes e probably
        aop = 30  # range 0-360 yes e maybe
        ma = 20   # range 0-360 yes e maybe
        p1 = Planet(1, ecc, 1000000000000, math.radians(inc), math.radians(loan), math.radians(aop), math.radians(ma))
        p1.CalculatePeriod(testMu)

        print(f"Period: {p1.period:2.0f} s = {(p1.period/60/60/24):2.3f} days\nAt t=0s:")

        #observerAngle = 270
        observerAngleRads = (3 * math.pi) / 2
        t0 = 0 # seconds
        
        #eccentricAnomalyAtEpoch = p1.CalculateCurrentEccentricAnomaly(t0)
        #trueAnomalyAtEpoch = p1.CalculateTrueAnomalyFromEccentric(eccentricAnomalyAtEpoch)
        truePos = p1.CalculateCartesianCoordsAtTime(t0)
        
        #TAinDegrees = math.degrees(trueAnomalyAtEpoch)
        #change observer angle based on loan/aop/inclination??
        #angleDistance =(observerAngle- TAinDegrees - math.degrees(p1.loan)-math.degrees(p1.aop)) % 360 # not used any more
        #loanTerm = 1 / math.cos(p1.loan)
        #incTerm = p1.inc
        #angleDistanceRads = ((observerAngleRads - p1.loan - p1.aop * math.cos(p1.inc))) % (2 * math.pi)
        #angleDistance = math.degrees(angleDistanceRads)
        
        #angleToTravelTo = trueAnomalyAtEpoch + angleDistanceRads
        #angleToTravelTo = angleDistanceRads

        experimentalMeanAnomalyAngle = (observerAngleRads - p1.ma * math.sin(p1.ecc) - p1.loan - p1.aop * math.cos(p1.inc)) % (2 * math.pi)
        experimentalTimeToPropogateTo = experimentalMeanAnomalyAngle / (2 * math.pi) * p1.period

        #print(f"Planet TA = {TAinDegrees} with observer angle of {observerAngle} gives angle distance of {angleDistance}, in rads is {angleDistanceRads}")
        print(f"Experimental MA angle: {math.degrees(experimentalMeanAnomalyAngle)} gives time of {experimentalTimeToPropogateTo} out of {p1.period}s")
        #predictedMean = mean_anomaly_from_true(p1.ecc,angleToTravelTo)
        #predictedEccentric = eccentric_anomaly_from_true(p1.ecc, angleToTravelTo)
        #predictedTrue = angleToTravelTo

        experimentalPosition = p1.CalculateCartesianCoordsAtTime(experimentalTimeToPropogateTo)
        iterativeEpoch, iterativePos = self.findIntersectionWithXYPlaneIteratively(p1)

        #print(f"Mean: {math.degrees(predictedMean)}, Eccentric: {math.degrees(predictedEccentric)}, true: {math.degrees(predictedTrue)}")
        #print(f"Position at origin time   :    {p1.CalculateCartesianPosition(trueAnomalyAtEpoch, eccentricAnomalyAtEpoch)}")
        print(f"Position at t=0   :    {truePos}")
        #pos = p1.CalculateCartesianPosition(predictedTrue, predictedEccentric)
        #isCorrect = abs(pos[0]) < 0.0000000001
        experimentCorrect = abs(experimentalPosition[0]) < 0.0000000001
        #print(f"Position at predicted time: {'ok' if isCorrect else 'no'} {pos}")
        print(f"Position at test t: {'ok' if experimentCorrect else 'no'} {experimentalPosition}")
        experimentCorrect = abs(iterativePos[0]) < 0.0000000001
        print(f"Position at iter t: {'ok' if experimentCorrect else 'no'} {iterativePos}")
        #print(f"found quadrant to be {self.findIntersectionWithXYPlaneIteratively(p1)}")

    def findIntersectionWithXYPlaneIteratively(self, planet, exp = 50):
        smallestPositive = math.inf
        smallestPositiveAngle = 0
        biggestNegative = -math.inf
        biggestNegativeAngle = 0

        smallestRecordedAbsDiff = math.inf
        smallestRecordedAbsDiffAngle = None
        smallestRecordedAbsDiffPos = None
        #gets me the quadrant the intersection is in
        for i in range(4):
            angleToSearch = (i * math.pi) / 2
            epoch = (angleToSearch / (2 * math.pi)) * planet.period
            positionAtEpoch = planet.CalculateCartesianCoordsAtTime(epoch)
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
            #print(f"testing angle {angleToSearch:<18}, time {epoch:<18}, ended with xcoord {xcoord:<22}, sPA: {smallestPositiveAngle:<18}, bNA: {biggestNegativeAngle:<18}, avgang: {(smallestPositiveAngle + biggestNegativeAngle) / 2:<18} quads")
        #print("")

        #gets me the 1/2^expth of circle, subtract two since already have quadrant from above
        #possibly want to redo this so that it goes until distance is less than radius of planet?
        for i in range(exp-2):
            angleToSearch = (smallestPositiveAngle + biggestNegativeAngle) / 2
            epoch = (angleToSearch / (2 * math.pi)) * planet.period
            positionAtEpoch = planet.CalculateCartesianCoordsAtTime(epoch)
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

            #print(f"testing angle {angleToSearch:<18}, time {epoch:<18}, ended with xcoord {xcoord:<22}, sPA: {smallestPositiveAngle:<18}, bNA: {biggestNegativeAngle:<18}, avgang: {(smallestPositiveAngle + biggestNegativeAngle) / 2:<18}  1/{math.pow(2,i+3)}")

        return smallestRecordedAbsDiffAngle, smallestRecordedAbsDiffPos


    #@unittest.skip("holup")
    def testDrawingOrbits(self):
        testMeanAngles = self.getTestPointsRadians(16)

        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 1,000,000,000kg
        testMu = G * 1000000000

        ecc = 0.4  # range 0-1   yes
        inc = 0  # range 0-360 not w loan
        loan = 0 # range 0-360 yes
        aop = 0  # range 0-360 yes
        ma = 0   # range 0-360 yes
        p1 = Planet(1, ecc, 1, math.radians(inc), math.radians(loan), math.radians(aop), math.radians(ma))
        p1.CalculatePeriod(testMu)

        resultPoints = []
        p1points = self.getResultPointsFromTimes(p1, 1)
        resultPoints = p1points

        
        p2 = Planet(1, ecc, 1, math.radians(inc), math.radians(loan), math.radians(aop + 90), math.radians(ma))
        p2.CalculatePeriod(testMu)

        p2points = self.getResultPointsFromTimes(p2, 3)
        resultPoints.extend(p2points)


        p3 = Planet(1, ecc, 1, math.radians(inc), math.radians(loan), math.radians(aop), math.radians(ma + 10))
        p3.CalculatePeriod(testMu)

        p3points = self.getResultPointsFromTimes(p3, 5)
        resultPoints.extend(p3points)
        
        resultPoints.append([0,0,0,6])

        fields = ['x', 'y', 'z', 'group']

        with open('orbit.csv','w') as f:
            write = csv.writer(f)
            write.writerow(fields)
            write.writerows(resultPoints)
        
        unitp1 = p1points[0][0:2] / np.linalg.norm(p1points[0][0:2])
        unitp2 = p2points[0][0:2] / np.linalg.norm(p2points[0][0:2])
        unitp3 = p3points[0][0:2] / np.linalg.norm(p3points[0][0:2])


        angle1 = np.arccos(np.clip(np.dot(unitp1, unitp2), -1.0, 1.0))
        angle2 = np.arccos(np.clip(np.dot(unitp2, unitp3), -1.0, 1.0))


        #print("\n")
        #print(f"p1: {p1points[0]}")
        #print(f"p2: {p2points[0]} angle1: {math.degrees(angle1)}")
        #print(f"p3: {p3points[0]} angle2: {math.degrees(angle2)}")


    def getResultPointsFromTimes(self, planet, group):
        resultPoints = []
        numTimes = 16
        for i in range(numTimes):
            time = i / numTimes * planet.period
            truePos = planet.CalculateCartesianCoordsAtTime(time)
            truePos.append(group-1 if time == 0 else group)
            resultPoints.append(truePos)
        return resultPoints



    def getResultPointsFromPlanetsPoints(self, planet, testPoints, group):
        resultPoints = []
        for meanAngle in testPoints:
            eccentricAnomalyAtEpoch = eccentric_anomaly_from_mean(planet.ecc, meanAngle)
            trueAngle = planet.CalculateTrueAnomalyFromEccentric(eccentricAnomalyAtEpoch)
            position = planet.CalculateCartesianPosition(trueAngle, eccentricAnomalyAtEpoch)

            position.append(group-1 if meanAngle == testPoints[0] else group)
            resultPoints.append(position)

        return resultPoints


    def getTestPointsRadians(self, numPoints):
        testPoints = []
        currentAngle = 0
        while currentAngle < 360:
            testPoints.append(math.radians(currentAngle))
            currentAngle += 360/numPoints
        return testPoints