

from TimeStep import TimeStep
from datetime import datetime
import unittest
from Thesis import uniformSourceResultAlgorithm
import lightkurve as lk
import CustomLightcurve
import math
import dateutil.parser as dateparser

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

class TestCustomLightcurves(unittest.TestCase):
    def test_creation(self):
        times = [9,8,7,6,5,4,3,2,1,0]
        epoch = datetime.now()
        test = CustomLightcurve.CustomLightcurve(times,epoch)
        for i in range(len(test.timeSteps)):
            test.timeSteps[i].flux = (i * 37) % 7
        sorted = test.sortByFlux()
        self.assertEqual(sorted.timeSteps[0].distanceToTimestep(sorted.timeSteps[1]), 7)
    
    def test_distance(self):
        times = [9,8,7,6,5,4,3,2,1,0]
        epoch = datetime.now()
        test = CustomLightcurve.CustomLightcurve(times,epoch)
        for i in range(len(test.timeSteps)):
            test.timeSteps[i].flux = (i * 37) % 7
        sorted = test.sortByFlux()

        self.assertAlmostEqual(sorted.timeSteps[0].distanceToTimestep(sorted.timeSteps[3]), 2.2, places=1)
        sorted.timeSteps[0].error = 2
        self.assertEqual(sorted.timeSteps[0].distanceToTimestep(sorted.timeSteps[3]), 1)
    
    def test_from_lightkurve(self):
        # generally dont include these tests because they take a relatively huge amount of time and cause warnings (from the library)
        include = False
        if(include):
            kurve = lk.read("DefaultFileTIC307210830C.fits")
            custom = CustomLightcurve.CustomLightcurve(kurve)
            self.assertEqual(kurve.time.iso[0], custom.epochTime)
            self.assertEqual(kurve.flux[0], custom.timeSteps[0].flux)
            self.assertEqual(kurve.flux_err[0], custom.timeSteps[0].error)