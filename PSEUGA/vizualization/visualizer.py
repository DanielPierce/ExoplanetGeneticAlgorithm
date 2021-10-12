
import matplotlib
import matplotlib.pyplot as plt
from PSEUGA.common.CustomLightcurve import CustomLightcurve

def createLightcurvePlot(lightcurve, path):
    lightcurve.sortByDate()
    x,y = lightcurve.toXY()

    f = plt.figure()
    f.set_figwidth(25)
    f.set_figheight(5)
    
    plt.plot(x,y)
    plt.xlabel('Seconds from Epoch')
    plt.ylabel('Flux')
    plt.savefig(path, bbox_inches='tight')
    #print(f'saved to {path}')

def createComparisonPlot(targetCurve, generatedCurve, path):
    
    targetCurve.sortByDate()
    tX, tY = targetCurve.toXY()

    generatedCurve.sortByDate()
    gX, gY = generatedCurve.toXY()

    f = plt.figure()
    f.set_figwidth(25)
    f.set_figheight(5)
    
    plt.plot(tX, tY)
    plt.plot(gX, gY, 'r')
    plt.xlabel('Seconds from Epoch')
    plt.ylabel('Flux')
    plt.savefig(path, bbox_inches='tight')
    #print(f'saved to {path}')
