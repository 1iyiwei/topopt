'''
tester for topology optimization code
'''
import time
import math

from loads import HalfBeam, Canti, Michell
from constraints import DensityConstraint
from fesolvers import CooFESolver
from topopt import Topopt
from plotting import Plot

if __name__ == "__main__":
    t = time.time()
    # material properties
    young = 1
    poisson = 0.3

    # constraints
    volfrac = 0.3
    xmin = 1e-9
    xmax = 1.0

    # input parameters
    nelx = 60*7
    nely = 110*7

    penal = 3
    rmin = 2.5

    delta = 0.02
    loopy = 300  # math.inf

    # loading/problem
    load = Michell(nelx, nely)

    # constraints
    density_constraint = DensityConstraint(volume_frac=volfrac, density_min=xmin, density_max=xmax)

    # optimizer
    verbose = True
    fesolver = CooFESolver(verbose=verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose=verbose)

    # compute
    history = True
    x = optimizer.init(load, density_constraint)
    x, x_more = optimizer.layout(load, density_constraint, x, penal, rmin, delta, loopy, history)

    print('Elapsed time is: ', time.time() - t, 'seconds.')

    if history:
        x_history = x_more
        loop = len(x_history)
    else:
        loop = x_more
        x_history = None

    # save
    if x_history:
        import imageio
        imageio.mimsave('topopt.gif', x_history)

    # plot
    pl = Plot(x, load, nelx, nely)
    pl.figure(title='loop '+str(loop))
    pl.boundary()
    pl.loading()
    pl.show()
