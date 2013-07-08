from numpy import *
import pylab as p
from scipy.integrate import ode
from scipy.linalg import eig
from matplotlib.patches import Ellipse
import json 




#############################
##    SYSTEM PARAMETERS    ##
#############################

# Model parameters
C = 1.0
L = 1.0
R = 0.2
Vin = 0.3

# Polynomial describing diode characteristic and its derivative
pol    = [83.72*32,   -226.31*16,   229.62*8,   -103.79*4,   17.76*2, 0]
pprime = [83.72*32*5, -226.31*16*4, 229.62*8*3, -103.79*4*2, 17.76*2 ]
pol    = [0.3*entry for entry in pol]
pprime = [0.3*entry for entry in pprime]

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
# X0 = array([0.45, 0.1])

# ODE Integration parameters
tmax = 9.01
num_points = 400
tolerance = 1e-08

# Initial ball size
e0 = 1e-04

# Bound on vector field
M = 0




##############################
##     ReachTrace CLASS     ##
##############################

# Class to contain ellipsoids covering a reach trace

class ReachTrace:
    
    X = []  # Centers of ellipsoids
    V = []  # V coordinate of center
    I = []  # I coordinate of center
    T = []  # Time
    d = []  # Diameter of ellipsoid (for circle -- change later)
    Jac = [] # Jacobian of system dynamics
    tolerance = 0 # Tolerance for ODE solution
    
    # Solve ODE when object is initialized 
    def __init__(self, f, Jac, X0, tmax, param):
        self.Jac = Jac
        num_points, self.tolerance = param
        r = ode(f, Jac).set_integrator('vode', method='bdf', with_jacobian=True, atol=tolerance)
        r.set_initial_value(X0, 0)
        dt = tmax/num_points
        t = 0
        self.X = [X0]
        self.V = [X0[0]]
        self.I = [X0[1]]
        self.T = [0]
        while r.successful() and r.t +dt < tmax:
            r.integrate(r.t+dt)
            t += 1
            self.X.append(r.y)
            self.V.append(r.y[0])
            self.I.append(r.y[1])
            self.T.append(r.t)

    # Perform Algorithm 3
    def algorithm3(self, e0, M):

        # Compute global maximum expansion rate
        Vtest = linspace(0,0.5,1000)
        muplot = []
        for Vi in Vtest:
            J = self.Jac(0,[Vi,0])
            muplot.append( max(real(eig(J+J.T)[0])) )
        index = argmax(muplot)
        mustar = muplot[index]
        Vstar = Vtest[index]
        
        self.d = [e0]
        c = []
        for i in range(len(self.T)-1):
            #  print 'time', T[i]
            # compute maximal expansion rate c_i in a neigbourhood of V[i] using global vector field bound M
            if abs(self.V[i] - Vstar) <= self.d[i] + M*(self.T[i+1]-self.T[i]):
                c.append(mustar)
            elif self.V[i] + self.d[i] + M*(self.T[i+1]-self.T[i]) < Vstar:
                J = Jac(0,[self.V[i] + self.d[i] + M*(self.T[i+1]-self.T[i]),0])
                c.append(max(real(eig(J+J.T)[0])))
            else:
                J = Jac(0,[self.V[i] - self.d[i] - M*(self.T[i+1]-self.T[i]),0])
                c.append(max(real(eig(J+J.T)[0])))
            # print 'c = ', c[i]
            
            
            # compute diameter of ball based on bound on expansion rate in neighbourhood of current state
            self.d.append(exp(c[i]*(self.T[i+1]-self.T[i]))*self.d[i]+self.tolerance)

   

##############################
##     TraceArray CLASS     ##
##############################


class TraceArray(list):
    # Contains a list of ReachTraces
    
    # Plot array of reach traces 
    def plotReachSet(self, NUM):
        fig = p.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlim(-0.1, 0.6)
        ax.set_ylim(-0.2, 0.6)
        for trace in self:
            for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                # print i
                e = Ellipse((trace.V[i],trace.I[i]), width=max(trace.d[i],0.002), height=max(trace.d[i],0.002), angle=0)
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_alpha(1)
                e.set_facecolor(p.rand(3))
        p.savefig('plotReachSet.pdf')

    # Plot all solutions of ODE
    def plotODE(self, NUM):
        fig = p.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlim(-0.1, 0.6)
        ax.set_ylim(-0.2, 0.6)
        for trace in self:
            Vplot = []
            Iplot = []
            for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                Vplot.append(trace.V[i])
                Iplot.append(trace.I[i])
            p.plot(Vplot, Iplot)
        p.savefig('plotODE.pdf')






    
param = [num_points, tolerance]

xmin = 0.45
xmax = 0.5
Earray = TraceArray()

i = 0
for x in linspace(xmax, xmin, 15):
    i += 1
    print 'Computing reach set from initial ball', i, '...'
    # Initial condition
    X0 = array([x, 0.1])
    trace = ReachTrace(f, Jac, X0, tmax, param)
    trace.algorithm3(e0, M)
    Earray.append(trace)

print 'Plotting reach sets...'
Earray.plotReachSet(400)

print 'Plotting ODE...' 
Earray.plotODE(400)

print 'Finished'




