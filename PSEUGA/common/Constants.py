import math

ATTRPERPLANET = 7
MAXPLANETS = 5
NUMPLANETS = MAXPLANETS * ATTRPERPLANET
STARRADIUS = NUMPLANETS + 1 # radius of star in km
STARMASS = STARRADIUS + 1 # mass of star in kg
STARBASEFLUX = STARMASS + 1 # Base flux of the star
DISTANCE = STARBASEFLUX + 1 # Distance to the star, in km

RADIUS = 0  # in km
ECC = 1     # eccentricity, unitless 1>x>0
SMA = 2     # semi-major axis in km
INC = 3     # inclination in degress, 180>x>=-180
LOAN = 4    # longitude of ascending node in radians, 2pi>x>=0
AOP = 5     # argument of periapsis in radians, 2pi>x>=0
MA = 6      # mean anomaly at epoch in radians, 2pi>x>=0 

MUTFACTOR = 0
MAX = 1
MIN = 2

CONSTANTS = [[12500, 500000, 1000],
             [0.25, 0.999, 0],
             [2*math.pow(10,9), 6*math.pow(10,9), 100000], #min is 1/570th the way to mercury, max is just past pluto
             [90, 180, -180],
             [0.5 * math.pi, 2 * math.pi, 0],
             [0.5 * math.pi, 2 * math.pi, 0],
             [0.5 * math.pi, 2 * math.pi, 0]]

NUMPLANETSMUTFACTOR = 2

GRAVITATIONALCONSTANT = 6.674 * math.pow(10, -11) # m^3 kg^-1 s^-2

STARMASSMUTFACTOR = math.pow(10,27)
STARMASSMIN = 1.898 * math.pow(10,27) * 100 # hypothetical minimum mass of star is 75x the mass of jupiter, went with 100 to exclude brown dwarfs
STARMASSMAX = 6.26 * math.pow(10,31) # mass of most massive star discovered, R136a1 at 315 solar masses per https://www.space.com/41313-most-massive-star.html

STARRADIUSMUTFACTOR = math.pow(10,6)
STARRADIUSMIN = math.pow(10,6)

STARBASEFLUXMUTFACTOR = 1000
STARBASEFLUXMIN = 1

DISTANCEMUTFACTOR = 9460730000000 * 3 # 3 lightyear
DISTANCEMIN = 9460730000000 * 5 # 5 lightyears

ATTRSPERMUTATION = 3