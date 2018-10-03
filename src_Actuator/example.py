"""
This is the main progam code that sets up the topology optimisation problem.
This optimisation tries to maximize the displacement at a set location. A so-
called compiant design.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""

import time
import math

from loads import Inverter
from constraints import DensityConstraint
from fesolvers import CvxFEA, SciPyFEA
from topopt import Topopt
from plotting import Plot

if __name__ == "__main__":
    # material properties
    young = 1
    Emin = 1e-9
    poisson = 0.3
    ext_stiff = 0.01

    # constraints
    volfrac = 0.3
    move = 1

    # mesh dimensions
    nelx = 40*3
    nely = 20*3

    # optimizer parameters
    penal = 3.0
    rmin = 1.5
    filt = 'density'
    loopy = 200  # math.inf
    delta = 0.005

    # plotting and printing options
    verbose = True
    plotting = True
    save_plot = True
    history = True

    # constraints object created
    den_con = DensityConstraint(nelx, nely, move, volume_frac=volfrac)

    # loading case object, other classes can be selected and created
    load = Inverter(nelx, nely, young, Emin, poisson, ext_stiff)

    # FEA object is generated, other solvers can be selected and created
    fesolver = CvxFEA(verbose=verbose)

    # create optimizer object and initialise the problem
    optimizer = Topopt(den_con, load, fesolver, verbose=verbose)

    # compute
    t = time.time()
    x, x_history = optimizer.layout(penal, rmin, delta, loopy, filt, history)
    print('Elapsed time is: ', time.time() - t, 'seconds.')

    # plotting
    directory = 'CT0001/'
    pl = Plot(load, directory)
    pl.boundary(load)
    pl.loading(load)

    if history:
        for i in x_history:
            pl.add(i, animated=True)
        pl.save('video')

    pl.add(x, animated=False)

    if save_plot:
        pl.save('figure')

    if plotting:
        pl.show()

    if save_pointcloud:
        pl.saveXYZ(x, x_size=60, thickness=1)
