###############################################################################
##
##  Figure2.py
##
##  Computes overapproximations of reachable set for tunnel diode oscillator
##  using algorithms described in "Reachability analysis of nonlinear
##  systems using matrix measures" by John Maidens and Murat Arcak
##
##  Used to produce Figure 2 from this paper
##
##  John Maidens
##  July 9, 2013
##
###############################################################################

import pylab as p
import json
from numpy import *
from scipy.integrate import ode
from scipy.linalg import eig, svd, inv
from matplotlib.patches import Ellipse
from scipy.optimize import minimize_scalar
from cvxpy import *
from numpy.ma import max



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
e0 = 1.5e-03

# Bound on vector field
M = 0







##############################
##     ReachTrace CLASS     ##
##############################

# Class to contain ellipsoids covering a reach trace
# Algorithms 3 and 4 are methods of this class

class ReachTrace:
    
    X = []  # Centers of ellipsoids
    V = []  # V coordinate of center
    I = []  # I coordinate of center
    T = []  # Time
    d1 = [] # Diameter of ellipsoid in first direction
    d2 = [] # Diameter of ellipsoid in second direction
    theta = [] # Rotation angle of ellipsoid
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
            muplot.append( 0.5*max(real(eig(J+J.T)[0])) )
        index = argmax(muplot)
        mustar = muplot[index]
        Vstar = Vtest[index]
        
        self.d1 = [e0]
        self.d2 = [e0]
        c = []
        for i in range(len(self.T)-1):
            # compute maximal expansion rate c_i in a neigbourhood of V[i] using global vector field bound M
            if abs(self.V[i] - Vstar) <= self.d1[i] + M*(self.T[i+1]-self.T[i]):
                c.append(mustar)
            elif self.V[i] + self.d1[i] + M*(self.T[i+1]-self.T[i]) < Vstar:
                J = Jac(0,[self.V[i] + self.d1[i] + M*(self.T[i+1]-self.T[i]),0])
                c.append(0.5*max(real(eig(J+J.T)[0])))
            else:
                J = Jac(0,[self.V[i] - self.d1[i] - M*(self.T[i+1]-self.T[i]),0])
                c.append(0.5*max(real(eig(J+J.T)[0])))            
            
            # compute diameter of ball based on bound on expansion rate in neighbourhood of current state
            self.d1.append(exp(c[i]*(self.T[i+1]-self.T[i]))*self.d1[i]+self.tolerance)
            self.d2.append(exp(c[i]*(self.T[i+1]-self.T[i]))*self.d2[i]+self.tolerance)
            self.theta.append(0)

                
    # Perform Algorithm 4
    def algorithm4(self, Gamma0, M):
        Gamma = Gamma0
        U, s, V = svd(inv(array(Gamma)))
        self.d1 = [sqrt(s[0])]
        self.d2 = [sqrt(s[1])]
        self.theta = [arccos(U[1,1])]
        for i in range(len(self.T)-1):
            
            def sdp(c):
                # For fixed c solve the semidefinite program for Algorithm 4
                
                # Variable Gamma_plus (for \Gamma_{i+1})
                Gamma_plus = variable(2,2,name='Gamma_plus')

                # Constraints
                c0 = belongs(Gamma_plus, semidefinite_cone)
                c1 = belongs(Gamma - Gamma_plus, semidefinite_cone)
                c2 = belongs(2*c*Gamma_plus - Gamma_plus*J - J.T*Gamma_plus, semidefinite_cone)

                # Objective function
                obj = -exp(-2*c*dT)*det_rootn(Gamma_plus)
                
                # Find solution
                p = program(minimize(obj), [c0, c1, c2])
                return p.solve(quiet = True)

            def f_Gamma(c):
                # Once the optimal c is found, find the ellipsoid shape matrix
                
                # Variable Gamma_plus (for \Gamma_{i+1})
                Gamma_plus = variable(2,2,name='Gamma_plus')
                
                # Constraints
                c0 = belongs(Gamma_plus, semidefinite_cone)
                c1 = belongs(Gamma - Gamma_plus, semidefinite_cone)
                c2 = belongs(2*c*Gamma_plus - Gamma_plus*J - J.T*Gamma_plus, semidefinite_cone)
                
                # Objective function
                obj = -exp(-2*c*dT)*det_rootn(Gamma_plus)
                
                # Find solution
                p = program(minimize(obj), [c0, c1, c2])
                p.solve(quiet = True)
                return Gamma_plus.value
                    
            # Search for c solving optimization problem using minimize_scalar to minimize function sdp
            J = matrix(Jac(self.T[i],self.X[i]))
            dT = self.T[i+1] - self.T[i]
            cmin = max(diag(J))
            cmax = cmin+1
            res = minimize_scalar(sdp, bounds=(cmin, cmax), method='bounded')
            cstar = res.x
            
            # Update Gamma
            Gamma = exp(-2*cstar*dT)*f_Gamma(cstar)
                    
            # Use Gamma to find width, height and angle of ellipsoid
            U, s, V = svd(inv(array(Gamma)))
            self.d1.append(sqrt(s[0]))
            self.d2.append(sqrt(s[1]))
            self.theta.append(arccos(U[1,1]))
                    
            # Update on our progress
            if i%50 == 0:
                print 'step', i, 'of', num_points

            





   

##############################
##     TraceArray CLASS     ##
##############################


class TraceArray(list):
    # Contains a list of ReachTraces
    
    # Plot array of reach traces 
    def plotReachSet(self, NUM, figname):
        fig = p.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlim(-0.1, 0.6)
        ax.set_ylim(-0.2, 0.6)
        for trace in self:
            for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                e = Ellipse((trace.V[i],trace.I[i]), width=trace.d1[i], height=trace.d2[i], angle=trace.theta[i])
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_alpha(1)
                e.set_facecolor(p.rand(3))
        for trace in self:
                e = Ellipse((trace.V[0],trace.I[0]), width=trace.d1[0], height=trace.d2[0], angle=trace.theta[0])
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_alpha(1)
                e.set_facecolor('r')
                e.set_edgecolor('r')
        p.savefig(figname)

    # Plot all solutions of ODE
    def plotODE(self, NUM, figname):
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
        p.savefig(figname)







#####################################
##     RUN ALGORITHMS AND PLOT     ##
#####################################

    
param = [num_points, tolerance]

# Set number minimum and maximum x values and number of traces to plot in between
xmin = 0.45
xmax = 0.5
num_traces = 15

num_plot = 400


# Perform Algorithm 3 for each initial ball
Earray_Alg3 = TraceArray()

i = 0
for x in linspace(xmax, xmin, num_traces):
    i += 1
    print 'Computing reach set with Alg 3 from initial ball', i, 'of', num_traces, '...'
    # Initial condition
    X0 = array([x, 0.1])
    trace = ReachTrace(f, Jac, X0, tmax, param)
    trace.algorithm3(e0, M)
    Earray_Alg3.append(trace)

# Save data in json format
data = []
for trace in Earray_Alg3:
    data.append([trace.V, trace.I, trace.d1, trace.d2, trace.theta])
output = open('data_plotReachSet_Algorithm_3.json', 'wb')
json.dump(data, output)
output.close()

# Plot the results
print 'Saving files...'
Earray_Alg3.plotODE(num_plot, 'plotReachSet_Brute_Force.pdf')
Earray_Alg3.plotReachSet(num_plot, 'plotReachSet_Algorithm_3.pdf')



# Perform Algorithm 4 for each initial ball
Earray_Alg4 = TraceArray()

i = 0
for x in linspace(xmax, xmin, num_traces):
    i += 1
    print 'Computing reach set with Alg 4 from initial ball', i, 'of', num_traces, '...'
    # Initial condition
    X0 = array([x, 0.1])
    trace = ReachTrace(f, Jac, X0, tmax, param)
    Gamma0 = 1/(e0*e0)*eye(2)

    trace.algorithm4(Gamma0, M)
    Earray_Alg4.append(trace)

# Save data in json format
data = []
for trace in Earray_Alg4:
    data.append([trace.V, trace.I, trace.d1, trace.d2, trace.theta])
output = open('data_plotReachSet_Algorithm_4.json', 'wb')
json.dump(data, output)
output.close()

# Plot the results
print 'Saving files...'
Earray_Alg4.plotReachSet(num_plot, 'plotReachSet_Algorithm_4.pdf')
print 'Finished'



