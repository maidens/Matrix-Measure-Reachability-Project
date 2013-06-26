from numpy import *
import pylab as p
from scipy import integrate


###############################
##    SETUP AND SOLVE ODE    ##
###############################



## Model parameters
C = 1.0
L = 1.0
R = 0.2
Vin = 0.3

# Polynomial describing diode characteristic and its derivative
pol    = [83.72*32,   -226.31*16,   229.62*8,   -103.79*4,   17.76*2, 0]
pprime = [83.72*32*5, -226.31*16*4, 229.62*8*3, -103.79*4*2, 17.76*2 ]


# System dynamics
def deriv(x,t=0):
    V = x[0]
    I = x[1]
    return array([(-polyval(pol,V) + I)/C, (-V - R*I + Vin)/L ])

# Horizon for solution
tmax = 14.5
num_points = 1000
t = linspace(0, tmax,  num_points)

# Initial condition
X0 = array([0.5, 0.5])

# Integrate system dynamics
X, infodict = integrate.odeint(deriv, X0, t, full_output=True)
V, I = X.T


###################################
##    CREATE PHASE PLANE PLOT    ##
###################################

# Create phase plane plot 
f1 = p.figure()
p.plot(V, I)
p.grid()
p.xlabel('V')
p.ylabel('I')
p.title('Tunnel diode oscillator')
f1.savefig('tunnel_diode_1.png')