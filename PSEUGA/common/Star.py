
import PSEUGA.common.Constants as const
import random


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
    