from cgi import test
import unittest
from PSEUGA.common.Planet import Planet
import math
import orbital

class TestPlanet(unittest.TestCase):
    def test_period(self):
        # Orbit semi-major axis is 10km, all else are not relevent
        p1 = Planet(1, 0, 10)
        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 10,000kg
        testMu = G * 10000
        p1.CalculatePeriod(testMu)
        self.assertAlmostEqual(p1.period, 7.69107 * math.pow(10,9), delta=100)

    # Expected results generated with the help of
    # http://orbitsimulator.com/sheela/kepler.htm
    # for going from MA to EA, and
    # https://www.vcalc.com/wiki/MichaelBartmess/True+Anomaly%2C+%60nu%60%5BE%5D
    # for going from EA to TA
    def test_calculate_true_anomaly(self):
        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 10,000kg
        testMu = G * 10000
        
        # Test with eccentricity of 0
        p1 = Planet(1, 0, 10)
        p1.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180] # angles (in degrees) of mean anomaly to test at
        expectedResults = [0, 60, 120, 180] # expected true anomaly (in degrees) of each test 
        self.calculatePlanetTrueAnomalyAtTestPoints(p1,testPoints, expectedResults, 'ecc = 0')

        # Test with eccentricity of 0.1
        p2 = Planet(1,0.1,10)
        p2.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180]
        expectedResults = [0, 70.523686162, 129.298289029, 180]
        self.calculatePlanetTrueAnomalyAtTestPoints(p2,testPoints, expectedResults, 'ecc = 0.1')

        # Test with eccentricity of 0.5
        p3 = Planet(1, 0.5, 10)
        p3.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180]
        expectedResults = [0, 118.815000927, 155.544022571, 180]
        self.calculatePlanetTrueAnomalyAtTestPoints(p3,testPoints, expectedResults, 'ecc = 0.5')

        # Test with mean anomaly at t=0 of p/3 radians, 60 degrees
        p4 = Planet(1, 0, 10, 0, 0, 0, (math.pi / 3))
        p4.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180]
        expectedResults = [60, 120, 180, 240]
        self.calculatePlanetTrueAnomalyAtTestPoints(p4,testPoints, expectedResults, 'ma = pi/3')

        
        # Test with mean anomaly of p/3 radians and ecc = 0.1
        p5 = Planet(1, 0.1, 10, 0, 0, 0, (math.pi / 3))
        p5.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180]
        expectedResults = [70.523686379, 129.298289029, 180, 230.701710971]
        self.calculatePlanetTrueAnomalyAtTestPoints(p5,testPoints, expectedResults, 'ma = pi/3, ecc = 0.1')

    def calculatePlanetTrueAnomalyAtTestPoints(self, planet, testPoints, expectedResults, message):
        for pointIndex in range(len(testPoints)):
            with self.subTest(msg=message, pointIndex=pointIndex):
                secondsAtEpoch = (testPoints[pointIndex] / 360) * planet.period
                trueAnomalyAtEpoch = math.degrees(planet.CalculateCurrentTrueAnomaly(secondsAtEpoch))
                self.assertAlmostEqual(trueAnomalyAtEpoch, expectedResults[pointIndex], places=6)

    @unittest.skip("Need to calculate correct expected results")
    def test_calculate_cartesian_position(self):
        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 10,000kg
        testMu = G * 10000
        
        # Test with eccentricity of 0
        p1 = Planet(1, 0, 10)
        p1.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180, 240, 300] # angles (in degrees) of true anomaly to test at
        expectedResults = [[10,0,0], [5, 8.660254037,0], [-5, 8.660254037, 0], [-10, 0, 0], [-5, -8.660254037, 0], [5, -8.660254037,0]] # expected cartesian position
        self.calculatePlanetCartesianPositionAtTestPoints(p1,testPoints, expectedResults, 'ecc = 0')

        # Test with eccentricity of 0.1
        p2 = Planet(1, 0.1, 10)
        p2.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180, 240, 300] # angles (in degrees) of true anomaly to test at
        # v recalculate these v
        expectedResults = [[10,0,0], [5, 8.660254037,0], [-5, 8.660254037, 0], [-10, 0, 0], [-5, -8.660254037, 0], [5, -8.660254037,0]] # expected cartesian position
        self.calculatePlanetCartesianPositionAtTestPoints(p2,testPoints, expectedResults, 'ecc = 0.1')
        
        # Test with eccentricity of 0.5
        p3 = Planet(1, 0.5, 10)
        p3.CalculatePeriod(testMu)
        testPoints = [0, 60, 120, 180, 240, 300] # angles (in degrees) of true anomaly to test at
        # v recalculate these v
        expectedResults = [[10,0,0], [5, 8.660254037,0], [-5, 8.660254037, 0], [-10, 0, 0], [-5, -8.660254037, 0], [5, -8.660254037,0]] # expected cartesian position
        self.calculatePlanetCartesianPositionAtTestPoints(p3,testPoints, expectedResults, 'ecc = 0.5')

    def calculatePlanetCartesianPositionAtTestPoints(self, planet, testPoints, expectedResults, message):
        for pointIndex in range(len(testPoints)):
            with self.subTest(msg=message, pointIndex=pointIndex):
                eccentricAnomalyAtPosition = planet.CalculateCurrentEccentricAnomaly(math.radians(testPoints[pointIndex]))
                cartesianCoordinates = planet.CalculateCartesianPosition(math.radians(testPoints[pointIndex]), eccentricAnomalyAtPosition)

                self.assertAlmostEqual(cartesianCoordinates[0], expectedResults[pointIndex][0])
                self.assertAlmostEqual(cartesianCoordinates[1], expectedResults[pointIndex][1])
                self.assertAlmostEqual(cartesianCoordinates[2], expectedResults[pointIndex][2])
