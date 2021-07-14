
from datetime import datetime, timedelta
import dateutil.parser as dateparser
import TimeStep
import lightkurve
import copy

class CustomLightcurve:
    def __init__(self, *args):
        if(len(args) == 2):
            self.timeSteps = list(map(TimeStep.TimeStep, args[0]))
            self.epochTime = args[1]
        elif(len(args) == 1):
            kurve = args[0]
            self.epochTime = kurve.time.iso[0]
            self.timeSteps = []
            for i in range(len(kurve.time.iso)):
                secondsFromEpoch = dateparser.parse(kurve.time.iso[i]) - dateparser.parse(kurve.time.iso[0])
                self.timeSteps.append(TimeStep.TimeStep(secondsFromEpoch, kurve.flux[i], kurve.flux_err[i]))

    @staticmethod
    def createFromCopy(targetCLC):
        newCurve = CustomLightcurve()
        newCurve.epochTime = targetCLC.epochTime
        newCurve.timeSteps = copy.deepcopy(targetCLC.timeSteps)
        for ts in newCurve.timeSteps:
            ts.flux = 0
            ts.error = 0
        return newCurve

    def printTimesteps(self):
        print(f"Epoch time: {self.epochTime}")
        for i in range(len(self.timeSteps)):
            self.timeSteps[i].printTimestep()
    
    def getTimes(self):
        times = []
        for i in range(len(self.timeSteps)):
            times.append(self.epochTime + timedelta(seconds=self.timeSteps[i].secondsFromEpoch))
        return times

    def sortByFlux(self):
        sortedCurve = CustomLightcurve(self.getTimes(), self.epochTime)
        sortedCurve.timeSteps = sorted(self.timeSteps, key= lambda x: (x.flux, x.secondsFromEpoch))
        return sortedCurve