

from datetime import datetime, timedelta
import math


class TimeStep:
    def __init__(self, *args):
        if not type(args[0]) is int:
            raise TypeError("Seconds from epoch must be an int") 
        self.secondsFromEpoch = args[0]
        self.flux = None
        self.error = 0
        if(len(args) >= 2):
            self.flux = args[1]
        if(len(args) == 3):
            self.error = args[2]

    def printTimestep(self):
        print(f"{self.secondsFromEpoch} seconds from epoch, reading {self.flux} e/s with an error of +/-{self.error} e/s")

    def distanceToTimestep(self, targetTimestep):
        timeDistance = self.secondsFromEpoch - targetTimestep.secondsFromEpoch
        fluxDistance = self.flux - targetTimestep.flux
        totalError = self.error + targetTimestep.error
        if(abs(fluxDistance) <= totalError):
            fluxDistance = 0
        return math.sqrt(timeDistance * timeDistance + fluxDistance * fluxDistance)