from matplotlib.patches import Ellipse
import cPickle
import pylab as p
from numpy import floor



###################################
##    CREATE PHASE PLANE PLOT    ##
###################################

# Load ball diameter data d

print 'Starting...'

pkl_file = open('data_alg3_short.pkl', 'rb')
d = cPickle.load(pkl_file)
# print d
pkl_file.close()

print 'Ellipsoid data loaded...'


# Load ODE solution data V, I, and T
pkl_file = open('data_odesolve_short.pkl', 'rb')
V, I, T = cPickle.load(pkl_file)
# print V
pkl_file.close()

print 'ODE data loaded'



NUM = 1500
num_points = len(T)
tmax = T[len(T)-1]
step = tmax/NUM
M = 3.5


ells = [Ellipse((V[i],I[i]), width=d[i], height=d[i], angle=0) for i in [int(floor(k*num_points/NUM)) for k in range(NUM)]]

fig = p.figure()
ax = fig.add_subplot(111, aspect='equal')
for e in ells:
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_alpha(1)
    e.set_facecolor(p.rand(3))

ax.set_xlim(-0.1, 0.6)
ax.set_ylim(-0.1, 1.1)

p.show()

