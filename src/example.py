'''
tester for topology optimization code
'''

import math
import matplotlib.pyplot as plt

from loads import HalfBeam
from fesolvers import LilFESolver, CooFESolver
from topopt import Topopt

if __name__ == "__main__":
    # material properties
    young = 1
    poisson = 0.3

    # default input parameters
    nelx = 180
    nely = 60
    volfrac = 0.4
    penal = 3.0
    rmin = 5.4
    delta = 0.01
    loopy = math.inf

    # loading/problem
    load = HalfBeam(nelx, nely)

    # optimizer
    verbose = True
    fesolver = CooFESolver(verbose = verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose = verbose)

    # compute
    x = optimizer.init(load, volfrac)
    x, loop = optimizer.layout(load, x, volfrac, penal, rmin, delta, loopy)

    # plot
    plt.figure()
    plt.imshow(x, cmap=plt.cm.gray)
    plt.title(str(loop) + ' loops')
    plt.xticks([])
    plt.yticks([])
    plt.show()
