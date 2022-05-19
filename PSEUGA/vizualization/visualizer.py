
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
    plt.ylabel('PDCSAP Flux')
    plt.savefig(path, bbox_inches='tight')
    plt.close()

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
    plt.close()

def createXCorrPlot(correlationMap, path):
    x = [i for i in range(len(correlationMap))]
    y = correlationMap

    f = plt.figure()
    f.set_figwidth(25)
    f.set_figheight(5)
    
    plt.plot(x,y)
    plt.xlabel('index')
    plt.ylabel('match')
    plt.savefig(path, bbox_inches='tight')
    plt.close()

def createNSGAPosPlot(pop, path):
    fitnesses = [indiv.fitness.values for indiv in pop]
    x  = [fit[0] for fit in fitnesses]
    y  = [fit[1] for fit in fitnesses]
    z  = [fit[2] for fit in fitnesses]
    
    fig = plt.figure()

    ax = fig.add_subplot(projection = '3d')

    ax.scatter(x,y,z)

    ax.set_xlabel('Cross Correlation Match')
    ax.set_ylabel('XC Distance from Center')
    ax.set_zlabel('Minimum Planetary Distance')
    
    plt.savefig(path,bbox_inches='tight')
    plt.close()