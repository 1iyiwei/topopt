'''
tester for topology optimization code
'''
import time
import math

from loads import HalfBeam, Canti, Michell, BiAxial
from constraints import DensityConstraint
from fesolvers import CvxFEA, SciPyFEA
from topopt import Topopt
from plotting import Plot

if __name__ == "__main__":
    t = time.time()
    # material properties
    young = 1
    poisson = 0.3

    # constraints
    Emin = 1e-9

    # input parameters
    nelx = 360
    nely = 250

    volfrac = 0.2

    penal = 3.0
    rmin = 2.5

    delta = 0.02
    loopy = 70

    # loading/problem
    load = Canti(nelx, nely)

    # constraints
    density_constraint = DensityConstraint(volume_frac=volfrac, Emin=Emin)

    # optimizer
    verbose = True
    fesolver = CvxFEA(verbose=verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose=verbose)

    # compute
    filt = 'sensitivity'
    history = True
    x = optimizer.init(load, density_constraint)
    x, x_more = optimizer.layout(load, density_constraint, x, penal, rmin, delta, loopy, filt, history)

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
    pl.figure()
    pl.loading()
    pl.boundary()
    pl.show()
