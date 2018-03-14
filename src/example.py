'''
tester for topology optimization code
'''
import time
import math
import matplotlib.pyplot as plt

from loads import HalfBeam
from constraints import DensityConstraint
from fesolvers import CooFESolver
from topopt import Topopt

if __name__ == "__main__":
    # material properties
    t = time.time()
    young = 1
    poisson = 0.3

    # constraints
    volfrac = 0.4
    xmin = 0.001
    xmax = 1.0

    # input parameters
    nelx = 300
    nely = 100

    penal = 3.0
    rmin = 3.0

    delta = 0.02
    loopy = math.inf

    # loading/problem
    load = HalfBeam(nelx, nely)

    # constraints
    density_constraint = DensityConstraint(volume_frac=volfrac, density_min=xmin, density_max=xmax)

    # optimizer
    verbose = True
    fesolver = CooFESolver(verbose=verbose)
    optimizer = Topopt(fesolver, young, poisson, verbose=verbose)

    # compute
    history = False
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
