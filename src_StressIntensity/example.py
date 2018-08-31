'''
tester for topology compliant optimization code
'''
import time
import math

from loads import EdgeCrack, DoubleEdgeCrack, CompactTension
from constraints import DensityConstraint
from fesolvers import CvxFEA, SciPyFEA
from topopt import Topopt
from plotting import Plot

if __name__ == "__main__":
    t = time.time()
    # material properties
    young = 1
    poisson = 0.3
    ext_stiff = 0.0

    # constraints
    Emin = 1e-9
    volfrac = 1.1
    move = 1

    # mesh dimensions
    nelx = 300
    nely = 320
    crack_length = 200

    # optimization settings
    penal = 1.0
    rmin = 1.2
    loopy = 1000 # math.inf
    delta = 0.01

    # loading/problem
    load = CompactTension(nelx, crack_length, ext_stiff)   

    # constraints
    den_con = DensityConstraint(load, move, volume_frac=volfrac, density_min=1, density_max=2.0, Emin=Emin)

    # optimizer
    verbose = True
    fesolver = CvxFEA(verbose=verbose)
    optimizer = Topopt(fesolver, load, den_con, young, poisson, verbose=verbose)

    # compute
    filt = 'sensitivity'
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
    pl.figure(den_con)
    pl.loading()
    pl.boundary()
    pl.show()
