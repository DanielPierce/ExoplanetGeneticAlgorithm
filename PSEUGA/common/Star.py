
import PSEUGA.common.Constants as const
import random
from PSEUGA.common.Functs import Clamp, Lerp


class Star:
    
    
    def __init__(self, radius, mass, flux, distance):
        self.radius = radius
        self.mass = mass
        self.flux = flux
        self.distanceTo = distance

    def PrettyPrint(self):
        print(f"Radius: {self.radius} km")
        print(f"Mass: {self.mass} kg")
        print(f"Base Flux: {self.flux} e/s")
        print(f"Distance: {self.distanceTo} km")

    def Randomize(self):
        self.radius = random.uniform(const.STARRADIUSMIN, const.STARRADIUSMIN * 100)
        self.mass = random.uniform(const.STARMASSMIN, const.STARMASSMIN * 100)
        self.flux = random.uniform(const.STARBASEFLUXMIN, const.STARBASEFLUXMIN * 100)
        self.distanceTo = random.uniform(const.DISTANCEMIN, const.DISTANCEMIN * 100)

    def Mutate(self):
        numTraitsToMutate = 1
        mutationThreshold = 1 - (numTraitsToMutate / 4)
        if random.random() > mutationThreshold:
            self.radius = self.MutateAttribute(self.radius, const.STARRADIUSMUTFACTOR, const.STARRADIUSMIN)
        if random.random() > mutationThreshold:
            self.mass = self.MutateAttribute(self.radius, const.STARMASSMUTFACTOR, const.STARMASSMIN)
        if random.random() > mutationThreshold:
            self.flux = self.MutateAttribute(self.flux, const.STARBASEFLUXMUTFACTOR, const.STARBASEFLUXMIN)
        if random.random() > mutationThreshold:
            self.distanceTo = self.MutateAttribute(self.distanceTo, const.DISTANCEMUTFACTOR, const.DISTANCEMIN)

    def MutateAttribute(self, attr, mutFactor, lowerBound=-float('inf'), upperBound=float('inf')):
        attr += Lerp(random.random(), -1 * mutFactor, mutFactor)
        attr = Clamp(self.radius, lowerBound, upperBound)
        return attr

    