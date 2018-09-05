"""
This is the main progam code that sets up the topology optimisation problem.
This optimisation tries to minimize the stress intensity factor. Which should
maximise the fatigue crack propegation life and thus the increase damage
tolerance.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""

import time
import math

from loads import EdgeCrack, DoubleEdgeCrack, CompactTension
from constraints import DensityConstraint
from fesolvers import CvxFEA, SciPyFEA
from topopt import Topopt
from plotting import Plot

if __name__ == "__main__":
    # material properties
    young = 1
    poisson = 0.3
    ext_stiff = 0.0

    # constraints
    Emin = 1e-9
    volfrac = 1.05
    move = 0.5

    # mesh dimensions
    nelx = 200
    nely = 200
    crack_length = 60

    # optimization parameters
    penal = 1.0
    rmin = 1.1
    filt = 'sensitivity'
    loopy = 300  # math.inf
    delta = 0.01

    # plotting and printing options
    verbose = True
    plotting = True
    save_plot = False
    history = True

    # loading case object, other classes can be selected and created
    load = CompactTension(nelx, crack_length, young, Emin, poisson, ext_stiff)   

    # constraints object created
    den_con = DensityConstraint(load, move, volume_frac=volfrac, density_min=1, density_max=2)

    # FEA object is generated, other solvers can be selected and created
    fesolver = CvxFEA(verbose=verbose)

    # create optimizer object and initialise the problem
    optimizer = Topopt(den_con, load, fesolver, verbose=verbose)

    # compute
    t = time.time()
    x, x_history = optimizer.layout(penal, rmin, delta, loopy, filt, history)
    print('Elapsed time is: ', time.time() - t, 'seconds.')

    # plotting
    pl = Plot(load)
    pl.loading(load)
    pl.boundary(load)
    pl.add(x)

    if save_plot:
        pl.save('topopt')

    if history:
        for i in x_history:
            pl.add(i, animated=True)
        pl.save('topopt')

    if plotting:
        pl.show()
