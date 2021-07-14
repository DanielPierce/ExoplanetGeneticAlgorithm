

import unittest
from Thesis import uniformSourceResultAlgorithm
import math

class TestUniformSourceAlgorithm(unittest.TestCase):

    def calcPP2ZZ2(self, rstar, rp, d):
        p = rp / rstar
        p2 = math.pow(p,2)

        z = d / rstar
        z2 = math.pow(z,2)
        return p, p2, z, z2

    def test_case_1(self):
        # Planet is not obscuring star at all, should be no loss in flux
        rstar = 100
        rp = 5
        d = 110
        
        p, p2, z, z2 = self.calcPP2ZZ2(rstar, rp, d)
        
        flux = 1000

        self.assertLess(1 + p, z)
        self.assertEqual(uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2, flux), 0)
    
    def test_case_2(self):
        # Part of planet covers disk of star
        rstar = 100
        rp = 5
        d = 100

        p, p2, z, z2 = self.calcPP2ZZ2(rstar, rp, d)

        flux = 1000
        
        self.assertLess(abs(1 - p), z)
        self.assertLessEqual(z, 1 + p)
        self.assertAlmostEqual(uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2, flux), 1.23673625)

    def test_case_3(self):
        # Entirety of planetary disk covers only part of star
        rstar = 100
        rp = 5
        d = 90

        p, p2, z, z2 = self.calcPP2ZZ2(rstar, rp, d)

        flux = 1000
        
        self.assertLessEqual(z, 1 - p)
        self.assertAlmostEqual(uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2, flux), p2)
        
    def test_case_4(self):
        # Entirety of stellar disk is obscured
        rstar = 100
        rp = 500
        d = 90

        p, p2, z, z2 = self.calcPP2ZZ2(rstar, rp, d)

        flux = 1000
        
        self.assertLessEqual(z, p - 1)
        self.assertEqual(uniformSourceResultAlgorithm(d,rp,rstar,z,z2,p,p2, flux), 1)

