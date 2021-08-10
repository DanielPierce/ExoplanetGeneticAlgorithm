
import matplotlib
import matplotlib.pyplot as plt
from PSEUGA.common.CustomLightcurve import CustomLightcurve

def createLightcurvePlot(lightcurve, path):
    print(type(lightcurve))
    x,y = lightcurve.toXY()
    #plt.plot(x,y)
    #plt.xlabel('Seconds from epoch')
    #plt.ylabel('Flux')
    #plt.savefig(path)