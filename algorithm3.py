from numpy import *
import pylab as p
from scipy.integrate import ode
from scipy.linalg import eig
from matplotlib.patches import Ellipse
import pickle 




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
pol    = [0.3*ent for ent in pol]
pprime = [0.3*ent for ent in pprime]

# System dynamics
def f(t,x):
    V = x[0]
    I = x[1]
    return array([(-polyval(pol,V) + I)/C, (-V - R*I + Vin)/L ])

# System Jacobian
def Jac(t,x):
    V = x[0]
    I = x[1]
    return array([ [-1/C*polyval(pprime,V), 1/C], [ -1/L, -R/L] ])

# Initial condition
X0 = array([0.5, 0.1])

# Integrate system dynamics
tmax = 10.01
num_points = 1e06
#num_points = 10000
tolerance = 1e-08
r = ode(f, Jac).set_integrator('vode', method='bdf', with_jacobian=True, atol=tolerance)
r.set_initial_value(X0, 0)
dt = tmax/num_points
t = 0
X = [X0]
V = [X0[0]]
I = [X0[1]]
T = [0]
while r.successful() and r.t +dt < tmax:
    r.integrate(r.t+dt)
    t += 1
    X.append(r.y)
    V.append(r.y[0])
    I.append(r.y[1])
    T.append(r.t)
    print r.t, r.y

# Save data V, I and T in pickle format

output = open('data_odesolve_short.pkl', 'wb')
pickle.dump([V,I,T], output)
output.close()



###############################
##    PERFORM ALGORITHM 3    ##
###############################

# Global bound M on vector field magnitude
M = 3.5

# Diameter of ball about initial condition
e0 = 2e-04

# Compute global maximum expansion rate
Vtest = linspace(0,0.5,1000)
muplot = []
for Vi in Vtest:
    J = Jac(0,[Vi,0])
    muplot.append( max(real(eig(J+J.T)[0])) )
index = argmax(muplot)
mustar = muplot[index]
Vstar = Vtest[index]

d = [e0]
c = []
for i in range(len(T)-1):
    print 'time', T[i]  
    # compute maximal expansion rate c_i in a neigbourhood of V[i] using global vector field bound M
    if abs(V[i] - Vstar) <= d[i] + M*(T[i+1]-T[i]):
        c.append(mustar)
    elif V[i] + d[i] + M*(T[i+1]-T[i]) < Vstar:
        J = Jac(0,[V[i] + d[i] + M*(T[i+1]-T[i]),0])
        c.append(max(real(eig(J+J.T)[0])))
    else:
        J = Jac(0,[V[i] - d[i] - M*(T[i+1]-T[i]),0])
        c.append(max(real(eig(J+J.T)[0])))
    print 'c = ', c[i]


    # compute diameter of ball based on bound on expansion rate in neighbourhood of current state
    d.append(exp(c[i]*(T[i+1]-T[i]))*d[i])


# Save ball diameter data d in pickle format 

output = open('data_alg3_short.pkl', 'wb')
pickle.dump(d, output) 
output.close()



f2 = p.figure()
p.plot(d)
f2.savefig('expansion_rate.png')



