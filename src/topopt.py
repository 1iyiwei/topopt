'''
topology optimization

2D only for now
'''

import numpy as np
import cv2 as cv
import math


class Topopt(object):
    '''
    young: young's modulus
    poisson: poisson ratio
    '''
    def __init__(self, fesolver, young=1, poisson=0.3, verbose=False):
        self.fesolver = fesolver
        self.young = young
        self.poisson = poisson
        self.dim = 2
        self.verbose = verbose

    # topology optimization
    def layout(self, load, constraint, x, penal, rmin, delta, loopy, history=False):

        loop = 0  # number of loop iterations
        change = 1.0  # maximum density change from prior iteration

        if history:
            x_history = [x]

        while (change > delta) and (loop < loopy):
            loop = loop + 1
            x, change, c = self.iter(load, constraint, x, penal, rmin)
            if self.verbose: print('It.: {0:4d},  Obj.: {1:8.2f},  ch.: {2:0.3f}'.format(loop, c, change), flush=True)
            if history: x_history.append(x)

        # done
        if history:
            return x, x_history
        else:
            return x, loop

    # initialization
    def init(self, load, constraint):
        (nelx, nely) = load.shape()
        # mean density
        return np.ones((nely, nelx))*constraint.volume_frac()

    # iteration
    def iter(self, load, constraint, x, penal, rmin):

        xold = x.copy()

        # element stiffness matrix
        ke = self.lk(self.young, self.poisson)

        # displacement via finite element analysis
        u = self.fesolver.displace(load, x, ke, penal)

        # compliance and derivative
        c, dc = self.comp(load, x, u, ke, penal)

        # filter
        dc = self.filt(x, rmin, dc)

        # update
        x = self.update(constraint, x, dc)

        # how much has changed?
        change = np.amax(abs(x-xold))

        return x, change, c

    # updated compliance algorithm
    def comp(self, load, x, u, ke, penal):
        nely, nelx = x.shape
        xe = x.T.flatten()  # flat list wich desities

        edof = load.edof(nelx, nely)[0]
        ue = u[edof]  # list with the displacements of the nodes of that element

        # calculating the compliance in 3 steps
        dot = np.dot(ke, ue.reshape((nelx*nely, 8, 1)))
        ce = np.sum(ue.T*dot[:, :, 0], axis=0)  # element compliance
        c = np.sum(ce * xe.T**penal)  # total compliance

        dc = -penal * (xe ** (penal-1)) * ce  # compliance derivative
        dc = dc.reshape((nelx, nely)).T
        return c, dc

    # new filter based upon C++ accelerated code
    def filt(self, x, rmin, dc):
        rminf = math.floor(rmin)
        nely, nelx = x.shape

        # define normalized convolution kernel based upon rmin
        size = rminf*2+1
        kernel = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                dis = np.sqrt((rminf-i)**2 + (rminf-j)**2)
                kernel[i, j] = np.max((0, rmin - dis))
        kernel = kernel/np.sum(kernel)  # normalisation

        # elementwise multiplication of x and dc
        xdc = dc*x
        xdcn = cv.filter2D(xdc, -1, kernel, borderType=cv.BORDER_REFLECT)
        dcn = xdcn/x

        return dcn

    # optimality criteria update
    def update(self, constraint, x, dc):
        volfrac = constraint.volume_frac()
        xmin = constraint.density_min()
        xmax = constraint.density_max()

        # ugly hardwired constants to fix later
        move = 0.2 * xmax
        l1 = 0
        l2 = 100000
        lt = 1e-4

        nely, nelx = x.shape
        while (l2-l1 > lt):
            lmid = 0.5*(l2+l1)
            xnew = np.multiply(x, np.sqrt(-dc/lmid))

            x_below = np.maximum(xmin, x - move)
            x_above = np.minimum(xmax, x + move)
            xnew = np.maximum(x_below, np.minimum(x_above, xnew))

            if (np.sum(xnew) - volfrac*nelx*nely) > 0:
                l1 = lmid
            else:
                l2 = lmid

        return xnew

    # element (local) stiffness matrix
    def lk(self, young, poisson):
        e = young
        nu = poisson
        k = np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
        ke = e/(1-nu**2) * \
            np.array([[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
                      [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
                      [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
                      [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
                      [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
                      [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
                      [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
                      [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
        return ke
