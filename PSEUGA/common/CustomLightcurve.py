
from datetime import datetime, timedelta
import dateutil.parser as dateparser
import PSEUGA.common.TimeStep as TimeStep
import lightkurve
import copy

class CustomLightcurve:
    def __init__(self, *args):
        if(len(args) == 2):
            epochTime = args[1]
            stepTimings = args[0]
            #stepTimings.sort()
            secondsFromEpochPerStep = list(map(lambda x: (x - epochTime).seconds, stepTimings))
            self.timeSteps = list(map(TimeStep.TimeStep, secondsFromEpochPerStep))
            self.epochTime = args[1]
            if not type(self.epochTime) is datetime:
                print(self.epochTime)
                raise TypeError(f"Epoch time must be a datetime, currently is {type(self.epochTime)}")
            self.sortByDate()
        elif(len(args) == 1):
            kurve = args[0]
            self.epochTime = dateparser.parse(kurve.time.iso[0])
            self.timeSteps = []
            for i in range(len(kurve.time.iso)):
                secondsFromEpoch = dateparser.parse(kurve.time.iso[i]) - dateparser.parse(kurve.time.iso[0])
                self.timeSteps.append(TimeStep.TimeStep(secondsFromEpoch.seconds, kurve.flux[i], kurve.flux_err[i]))
            if not type(self.epochTime) is datetime:
                print(self.epochTime)
                raise TypeError(f"Epoch time must be a datetime, currently is {type(self.epochTime)}")
        elif(len(args) == 0):
            self.epochTime = None
            self.timeSteps = []

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

    def sortByDate(self):
        self.timeSteps = sorted(self.timeSteps, key= lambda x: (x.secondsFromEpoch))

    def setBaseFlux(self, baseFlux):
        for i in range(len(self.timeSteps)):
            self.timeSteps[i].flux = baseFlux

    def getUnskippedTimesteps(self, timestepsToSkip):
        unskipped = []
        for i in range(len(self.timeSteps), timestepsToSkip):
            unskipped.append(self.timeSteps[i])
        return unskipped

    def toXY(self):
        x = []
        y = []
        for step in self.timeSteps:
            x.append(self.epochTime + timedelta(seconds=step.secondsFromEpoch))
            y.append(step.flux)
        return x,y