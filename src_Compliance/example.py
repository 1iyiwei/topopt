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
    volfrac = 0.4
    move = 1

    # input parameters
    nelx = 50
    nely = 40

    # optimizer parameters
    penal = 3.0
    rmin = 1.5
    filt = 'density'
    loopy = 10  # math.inf
    delta = 0.001

    # plotting and printing options
    verbose = True
    plotting = True
    history = True

    # constraints object created
    den_con = DensityConstraint(nelx, nely, move, volume_frac=volfrac)

    # loading case object, other classes can be selected and created
    load = Canti(nelx, nely, young, Emin, poisson)

    # FEA object is generated, other solvers can be selected and created
    fesolver = CvxFEA(verbose=verbose)

    # create optimizer object and initialise the problem
    optimizer = Topopt(den_con, load, fesolver, verbose=verbose)

    # execute the optimization
    print('Elapsed time is: ', time.time() - t, 'seconds.')
    x, x_more = optimizer.layout(penal, rmin, delta, loopy, filt, history)
    print('Elapsed time is: ', time.time() - t, 'seconds.')

    if history:
        x_history = x_more
        loop = len(x_history)
    else:
        loop = x_more
        x_history = None

    # plot
    pl = Plot(x, load, nelx, nely, plotting)
    pl.figure()
    pl.loading()
    pl.boundary()

    # save
    if x_history:
        import imageio
        imageio.mimsave('topopt.gif', x_history)
