import unittest
from common.Planet import Planet
import math

class TestPlanet(unittest.TestCase):
    def test_period(self):
        # Orbit semi-major axis is 10km, all else are not relevent
        p1 = Planet(1, 0, 10, 0, 0, 0, 0)
        G = 6.674 * math.pow(10, -11)
        # Stellar mass is 10,000kg
        testMu = G * 10000
        p1.CalculatePeriod(testMu)
        self.assertAlmostEqual(p1.period, 7.69107 * math.pow(10,9), delta=100)