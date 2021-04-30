import math
import numpy as np
import Constants as const




def calculatePosition(individual, planet, radius, trueAnomaly):
    posX = radius * math.cos(trueAnomaly)
    posY = radius * math.sin(trueAnomaly)
    posZ = 0

    coscos = math.cos(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.cos(individual[planet * const.ATTRPERPLANET + const.LOAN])
    cossin = math.cos(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.sin(individual[planet * const.ATTRPERPLANET + const.LOAN])
    sincos = math.sin(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.cos(individual[planet * const.ATTRPERPLANET + const.LOAN])
    sinsin = math.sin(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.sin(individual[planet * const.ATTRPERPLANET + const.LOAN])

    posXDot = posX * coscos - sinsin * math.cos(individual[planet * const.ATTRPERPLANET + const.INC]) - posY * sincos + cossin * math.cos(individual[planet * const.ATTRPERPLANET + const.INC])
    posYDot = posX * cossin + sincos * math.cos(individual[planet * const.ATTRPERPLANET + const.INC]) + posY * coscos * math.cos(individual[planet * const.ATTRPERPLANET + const.INC]) - sinsin
    posZDot = posX * math.sin(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.sin(individual[planet * const.ATTRPERPLANET + const.INC]) + posY * math.cos(individual[planet * const.ATTRPERPLANET + const.AOP]) * math.sin(individual[planet * const.ATTRPERPLANET + const.INC])

    # intertial reference frame, but relative to what? orbital frame has z axis perpendicular to orbital plane and x axis pointing to periapsis of orbit
    return [posXDot, posYDot, posZDot]
