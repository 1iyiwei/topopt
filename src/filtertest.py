'''
tester for topology optimization code
'''

import numpy as np
import math
import matplotlib.pyplot as plt

from loads import HalfBeam
from constraints import DensityConstraint
from fesolvers import LilFESolver, CooFESolver
from topopt import Topopt

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

if __name__ == "__main__":
    # material properties
    young = 1
    poisson = 0.3

    # constraints
    volfrac = 0.4
    xmin = 0.001
    xmax = 1.0

    # input parameters
    nelx = 180
    nely = 60

    penal = 3.0
    rmin = 5.4

    delta = 0.02
    loopy = 2  # math.inf

    # loading/problem
    load = HalfBeam(nelx, nely)

    # constraints
    density_constraint = DensityConstraint(volume_frac = volfrac, density_min = xmin, density_max = xmax)

    # optimizer
    verbose = True
    fesolver = CooFESolver(verbose = verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose = verbose)

    # compute new filter result
    history = False
    newfilt = True
    xn = optimizer.init(load, density_constraint)
    xn, loop, cnew = optimizer.layout(load, density_constraint, xn, penal, rmin, delta, loopy, history, newfilt)

    # compute old filter result
    history = False
    newfilt = False
    xo = optimizer.init(load, density_constraint)
    xo, loop, cold = optimizer.layout(load, density_constraint, xo, penal, rmin, delta, loopy, history, newfilt)    

    # print copmlicance difference
    print('Compliance old ', cold)
    print('Compliance new ', cnew)
    print('  Relative differnce ', 100*(cnew-cold)/(cold))

    x = xo - xn
    xabs = max(np.absolute(x).min(), x.max())

    # plot
    plt.figure()
    plt.imshow(x, cmap=plt.cm.bwr)
    plt.title(r'$x_{old}-x_{new}$, $r_{min}=%0.1f$, $nely,nelx=%d,%d$' % (rmin, nely, nelx))
    plt.colorbar()
    plt.clim(-xabs, xabs)
    plt.xticks([])
    plt.yticks([])
    plt.show()
