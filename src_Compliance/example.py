'''
tester for topology optimization code
'''
import time
import math

from loads import HalfBeam, Beam, Canti, Michell, BiAxial
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
    volfrac = 0.5
    move = 1

    # input parameters
    nelx = 500
    nely = 100

    penal = 3.0
    rmin = 1.5

    loopy = 1000  # math.inf
    delta = 0.02

    # loading/problem
    load = Canti(nelx, nely)

    # constraints5
    den_con = DensityConstraint(load, move, volume_frac=volfrac, Emin=Emin)

    # optimizer
    verbose = True
    fesolver = CvxFEA(verbose=verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose=verbose)

    # compute
    filt = 'density'
    history = True
    x = optimizer.init(load, den_con)
    x, x_more = optimizer.layout(load, den_con, x, penal, rmin, delta, loopy, filt, history)

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
