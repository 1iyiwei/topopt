'''
tester for topology optimization code
'''

import numpy as np
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
    delta = 0.02
    loopy = math.inf

    # loading/problem
    load = HalfBeam(nelx, nely)

    # optimizer
    verbose = True
    fesolver = CooFESolver(verbose = verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose = verbose)

    # compute
    history = False
    x = optimizer.init(load, volfrac)
    if history:
        x, x_history = optimizer.layout(load, x, volfrac, penal, rmin, delta, loopy, history)
        loop = len(x_history)
    else:    
        x, loop = optimizer.layout(load, x, volfrac, penal, rmin, delta, loopy, history)
        x_history = None

    # save
    if x_history:
        # import scipy.misc
        # sequence = [scipy.misc.toimage(x, cmin=0, cmax=1) for x in x_history]
        import imageio
        imageio.mimsave('topopt.gif', x_history)

    # plot
    plt.figure()
    plt.imshow(x, cmap=plt.cm.gray)
    plt.title(str(loop) + ' loops')
    plt.xticks([])
    plt.yticks([])
    plt.show()
