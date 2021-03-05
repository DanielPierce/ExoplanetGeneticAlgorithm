
ATTRPERPLANET = 7
MAXPLANETS = 20
NUMPLANETS = MAXPLANETS * ATTRPERPLANET
STARRADIUS = NUMPLANETS + 1

RADIUS = 0  # in km
ECC = 1     # eccentricity, unitless 1>x>0
SMA = 2     # semi-major axis in km
INC = 3     # inclination in degress, 180>x>=-180
LOAN = 4    # longitude of ascending node in degrees, 360>x>=0
AOP = 5     # argument of periapsis in degrees, 360>x>=0
TA = 6      # true anomaly in degrees, 360>x>=0 

MUTFACTOR = 0
MAX = 1
MIN = 2

CONSTANTS = [[10, 500000, 1000],
             [0.001, 1, 0],
             [1000, 70000000000, 10000],
             [0.1, 180, -180],
             [0.1, 360, 0],
             [0.1, 360, 0],
             [0.1, 360, 0]]

NUMPLANETSMUTFACTOR = 1
STARRADIUSMUTFACTOR = 10000
