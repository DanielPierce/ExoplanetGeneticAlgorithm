import math

ATTRPERPLANET = 7
MAXPLANETS = 20
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

CONSTANTS = [[10, 500000, 1000],
             [0.001, 0.999, 0],
             [1000, 70000000000, 100000],
             [0.1, 180, -180],
             [0.001, 2 * math.pi, 0],
             [0.001, 2 * math.pi, 0],
             [0.001, 2 * math.pi, 0]]

NUMPLANETSMUTFACTOR = 1
STARRADIUSMUTFACTOR = 1000
STARMASSMUTFACTOR = math.pow(10,25)
STARBASEFLUXMUTFACTOR = 1000
DISTANCEMUTFACTOR = 9460730000000 # 1 lightyear

GRAVITATIONALCONSTANT = 6.674 * math.pow(10, -11) # m^3 kg^-1 s^-2

STARMASSMIN = 1.898 * math.pow(10,27) * 50 # hypothetical minimum mass of star is 75x the mass of jupiter, went with 50x to include brown dwarfs
STARRADIUSMIN = STARRADIUSMUTFACTOR * 1000
STARBASEFLUXMIN = 1000 
DISTANCEMIN = 9460730000000 * 5 # 5 lightyears

skippedTimesteps = 0
ATTRSPERMUTATION = 3