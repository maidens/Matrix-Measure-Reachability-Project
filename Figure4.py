###############################################################################
##
##  Figure4.py
##
##  Computes overapproximations of reachable set for transcription cascade
##  using algorithms described in "Reachability analysis of nonlinear
##  systems using matrix measures" by John Maidens and Murat Arcak
##
##  Used to produce Figure 4 from this paper
##
##  John Maidens
##  July 17, 2013
##
###############################################################################

import pylab as p
import json
import random
import time 
from numpy import *
from scipy.integrate import ode
from scipy.linalg import eig
from scipy import rand
from scipy.optimize import minimize_scalar
from matplotlib.patches import Ellipse, Polygon
from cvxpy import sum, less_equals, variable, program, geo_mean, maximize



#############################
##    SYSTEM PARAMETERS    ##
############################# 


# Other parameters
num_points = 50
tolerance = 1e-08
tmax = float(1)

# Function to generate random system equations for cascade of length n
def generate_system(n):
    
    # Generate random parameter values 
    delta = zeros(n)
    k1 = zeros(n)
    k2 = zeros(n)
    e1 = zeros(n)
    E = ones(n)
    e1[0] = 1
    A_c = zeros((n,n))
    for i in range(n):
        delta[i] = 10*random.random()
        k1[i] = 10*random.random()
        k2[i] = 10*random.random()
        if not i == 0:
            A_c[i,i-1] = 10
    print 'delta =', delta
    print 'k1 =', k1
    print 'k2 =', k2
    print 'A_c =', A_c

    # Define dynamics f
    def f(t,z):
        x = z[0:n]
        y = z[n:2*n]
        dx = dot(A_c,y) - delta*x + k1*y - k2*(E-y)*x + 3*(sin(10*t)+1)*e1
        dy = -k1*y + k2*(E-y)*x
        return hstack((dx,dy))

    # Define Jacobian of f
    def Jac(t,z):
        x = z[0:n]
        y = z[n:2*n]
        Jxx =  - diag(delta) - diag(k2*(E-y))
        Jxy = diag(k1) + diag(k2*x) + A_c
        Jyx = diag(k2*(E-y))
        Jyy = -diag(k1) - diag(k2*x)
        return vstack( ( hstack((Jxx,Jxy)),hstack((Jyx,Jyy)) ) )

    # Compute diagonal weightings of 1-norm in which f is contractive
        # c_tilde[n] is the contraction rate
    p_vec = []
    c_vec = []
    c_tilde = []
    q_vec = []
    for i in range(n):
        p_vec.append(1+delta[i]/(2*k2[i]))
        c_vec.append( min([(1-1/p_vec[i])*k1[i], delta[i]/2]) )
        if i == 0:
            c_tilde.append( c_vec[0] )
        else:
            c_tilde.append( min([0.9*c_tilde[i-1], c_vec[i]]) )
            q_vec.append(0.9*c_tilde[i-1])
    d = [1, p_vec[0]]
    for i in range(n-1):
        d.append(q_vec[i])
        d.append(q_vec[i]*p_vec[i+1])
    print 'd =', d
    x0 = ones(n)
    y0 = 0.8*ones(n)
    X0 = hstack((x0,y0))
    return f, Jac, d, c_tilde[n-1]


# function to compute 1-measure of square matrix 
def mu_1(J):
    n, m = J.shape
    assert n == m
    return max([J[i,i] + sum([abs(J[i,j]) for j in range(n)]) - abs(J[i,i]) for i in range(n)])







##############################
##     ReachTrace CLASS     ##
##############################

# Class to contain ellipsoids covering a reach trace
# Algorithm 3 is a method of this class

class ReachTrace:
    
    X = []  # Centers of ellipsoids
    x = []  # x coordinates of center
    y = []  # y coordinates of center
    T = []  # Time
    d1 = [] # Diameter of ellipsoid in first direction
    d2 = [] # Diameter of ellipsoid in second direction
    theta = [] # Rotation angle of ellipsoid
    
    d_norm1 = [] # Vector of 1-norm ball radii in each of d directions
    Jac = [] # Jacobian of system dynamics
    tolerance = 1e-08 # Tolerance for ODE solution
    
    # Solve ODE when object is initialized
    
    def __init__(self, f, Jac, X0, tmax, param):
        self.Jac = Jac
        num_points, self.tolerance = param
        r = ode(f, Jac).set_integrator('vode', method='bdf', with_jacobian=True, atol=tolerance)
        r.set_initial_value(X0, 0)
        dt = tmax/num_points
        t = 0
        self.X = [X0]
        self.x = [X0[0:n]]
        self.y = [X0[n:2*n]]
        self.T = [0]
        
        while r.successful() and r.t +dt < tmax:
            r.integrate(r.t+dt)
            t += 1
            self.X.append(r.y)
            self.x.append(r.y[0:n])
            self.y.append(r.y[n:2*n])
            self.T.append(r.t)


    # Perform Algorithm 3
    def algorithm3(self, d0, M):
        self.d_norm1 = [d0]
        c = []
        for i in range(len(self.T)-1):
            J = self.Jac(self.T[i], self.X[i])
            
            # define expansion rate c
            c.append(mu_1(diag(d)*J*diag([1/ent for ent in d])))
            
            
            # compute diameter of ball based on bound on expansion rate in neighbourhood of current state
            self.d_norm1.append( [exp(-c[i]*(self.T[i+1]-self.T[i]))*ent + self.tolerance for ent in self.d_norm1[i]])
        print 'c =', c


    # Perform Algorithm 4    
    def algorithm4(self, d0, M):
        d = d0
        # U, s, V = svd(inv(array(d)))
        self.d_norm1 = [1/d0]
        
        for i in range(len(self.T)-1):
            
            def sdp(c):
                # For fixed c solve the semidefinite program for Algorithm 4
                
                # Variable Gamma_plus (for \Gamma_{i+1})
                d_plus = variable(2*n, 1, name='d_plus')
                
                # Constraints
                constr = []
                for j in range(2*n):
                    ctt =  less_equals(d_plus[j,0], 0)
                    constr.append( less_equals(-d_plus[j,0], 0))
                    constr.append( less_equals( d_plus[j,0], d[j,0]))
                    constr.append( less_equals( J[j,j]*d_plus[j,0] + sum([abs(J[i,j])*d_plus[i,0] for i in range(2*n)]) - abs(J[j,j])*d_plus[j,0], c* d_plus[j,0]))
                
                # Objective function
                obj = geo_mean(d_plus)
                # Find solution
                p = program(maximize(obj), constr)
                
                return exp(c*dT)/p.solve(quiet = False)
                    
                    
            def f_Gamma(c):
                # For fixed c solve the semidefinite program for Algorithm 4
                
                # Variable Gamma_plus (for \Gamma_{i+1})
                d_plus = variable(2*n, 1, name='d_plus')
                
                # Constraints
                constr = []
                for j in range(2*n):
                    ctt =  less_equals(d_plus[j,0], 0)
                    constr.append( less_equals(-d_plus[j,0], 0))
                    constr.append( less_equals( d_plus[j,0], d[j,0]))
                    constr.append( less_equals( J[j,j]*d_plus[j,0] + sum([abs(J[i,j])*d_plus[i,0] for i in range(2*n)]) - abs(J[j,j])*d_plus[j,0], c* d_plus[j,0]))
                
                # Objective function
                obj = geo_mean(d_plus)
                # Find solution
                p = program(maximize(obj), constr)

                return d_plus.value
            
            # Search for c solving optimization problem using minimize_scalar to minimize function sdp
            J = matrix(Jac(self.T[i],self.X[i]))
            dT = self.T[i+1] - self.T[i]
            cmin = max(diag(J))
            cmax = cmin+1
            res = minimize_scalar(sdp, bounds=(cmin, cmax), method='bounded')
            cstar = res.x
            self.d_norm1.append(f_Gamma(cstar))










##############################
##     TraceArray CLASS     ##
##############################


class TraceArray(list):
    # Contains a list of ReachTraces
    
    # Plot array of reach traces for 2-norm data
    def plotReachSet_norm2(self, NUM, figname):
        fig = p.figure()
        for j in range(n):

            ax = fig.add_subplot(3,3,j)# , aspect='equal')
            ax.set_xlim(0, 4)
            ax.set_ylim(0, 1)
            for trace in self:
                for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                    
                    e = Ellipse((trace.x[i][j],trace.y[i][j]), width=trace.d1[i], height=trace.d2[i], angle=trace.theta[i])
                    ax.add_artist(e)
                    e.set_clip_box(ax.bbox)
                    e.set_alpha(1)
                    e.set_facecolor(p.rand(3))
            for trace in self:
                    e = Ellipse((trace.x[0][j],trace.y[0][j]), width=trace.d1[0], height=trace.d2[0], angle=trace.theta[0])
                    ax.add_artist(e)
                    e.set_clip_box(ax.bbox)
                    e.set_alpha(1)
                    e.set_facecolor('r')
                    e.set_edgecolor('r')
        p.savefig(figname)
    
    # Plot array of reach traces for 1-norm data 
    def plotReachSet_norm1(self, NUM, figname):
        fig = p.figure()
        for j in range(n):
            ax = fig.add_subplot(2,2,j+1 , aspect='equal')
            ax.set_xlim(0, 4)
            ax.set_ylim(0, 1)
            ax.set_xlabel('$x_'+str(j+1)+'$')
            ax.set_ylabel('$y_'+str(j+1)+'$')
            for trace in self:
                for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                    verts = [(trace.x[i][j] + 1/trace.d_norm1[i][2*j], trace.y[i][j]                          ),
                             (trace.x[i][j]                            , trace.y[i][j] - 1/trace.d_norm1[i][2*j+1]),
                             (trace.x[i][j] - 1/trace.d_norm1[i][2*j], trace.y[i][j]                          ),
                             (trace.x[i][j]                              , trace.y[i][j] + 1/trace.d_norm1[i][2*j+1])]
                    poly = Polygon(verts, facecolor='0.8', edgecolor='k')
                    ax.add_artist(poly)
                    poly.set_clip_box(ax.bbox)
                    poly.set_alpha(1)
                    if i==0:
                        poly.set_facecolor('r')
                    else:
                        poly.set_facecolor(p.rand(3))

        p.savefig(figname)

    # Plot all solutions of ODE
    def plotODE(self, NUM, figname):
        fig = p.figure()
        for j in range(n):
            ax = fig.add_subplot(3,3,j+1) #, aspect='equal')
            ax.set_xlim(0, 4)
            ax.set_ylim(0, 1)
            for trace in self:
                xplot = []
                yplot = []
                for i in [int(floor(k*len(trace.T)/NUM)) for k in range(NUM)]:
                    xplot.append(trace.x[i][j])
                    yplot.append(trace.y[i][j])
            p.plot(xplot, yplot)
        p.savefig(figname)







########################################################
##     RUN ALGORITHMS AND CREATE SCALABILITY PLOT     ##
########################################################


# list to contain computaton times
times_list = []

# list to contain state dimensions 
dimension_list = []

for i in range(50):
    n = i+1  # number of cascaded systems 
    dimension_list.append(2*n)
    
    # Dynamics
    f, Jac, d, c = generate_system(n)

    # Initial condition
    x0 = ones(n)
    y0 = 0.8*ones(n)
    X0 = hstack((x0,y0))
    
    # tic
    tic = time.clock()

    # Compute simulation trace
    param = [num_points, tolerance]
    trace = ReachTrace(f, Jac, X0, tmax, param)

    # Run Algorithm 3
    e0 = 5
    M = 0
    print 'n =', n
    d0 = [e0*ent for ent in d]
    trace.algorithm3(d0, M)

    # toc 
    toc = time.clock()
    times_list.append(toc - tic)

# plot the computation times 
p.plot(dimension_list, times_list)
p.xlabel('state dimension')
p.ylabel('computation time (s)')
p.savefig('scalability.pdf')
