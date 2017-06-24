'''
topology optimization

2D only for now
'''

import numpy as np
import math

class Topopt(object):
    
    '''
    young: young's modulus
    poisson: poisson ratio
    '''
    def __init__(self, fesolver, young = 1, poisson = 0.3, verbose = False):
        self.fesolver = fesolver
        self.young = young
        self.poisson = poisson
        self.dim = 2
        self.verbose = verbose

    # topology optimization
    def layout(self, load, x, volfrac, penal, rmin, delta, loopy, history = False):

        loop = 0 # number of loop iterations
        change = 1.0 # maximum density change from prior iteration
        
        if history:
            x_history = [x]

        while (change > delta) and (loop < loopy):
            loop = loop + 1
            x, change = self.iter(load, x, volfrac, penal, rmin)
            if self.verbose: print('iteration ', loop, ', change ', change, flush = True)
            if history: x_history.append(x)

        # done
        if history:
            return x, x_history
        else:
            return x, loop

    # initialization
    def init(self, load, volfrac):
        (nelx, nely) = load.shape()
        # mean density
        return np.ones((nely, nelx))*volfrac

    # iteration
    def iter(self, load, x, volfrac, penal, rmin):

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
        x = self.update(x, volfrac, dc)

        # how much has changed?
        change = np.amax(abs(x-xold))

        return x, change
    
    # compliance and its derivative
    def comp(self, load, x, u, ke, penal):
        c = 0
        dc = np.zeros(x.shape)

        nely, nelx = x.shape
        for ely in range(nely):
            for elx in range(nelx):
                ue = u[load.edof(elx, ely, nelx, nely)]
                ce = np.dot(ue.transpose(), np.dot(ke, ue))
                c = c + (x[ely,elx]**penal)*ce
                dc[ely,elx] = -penal*(x[ely,elx]**(penal-1))*ce

        return c, dc

    # filter
    def filt(self, x, rmin, dc):
        rminf = math.floor(rmin)

        dcn = np.zeros(x.shape)
        nely, nelx = x.shape

        for i in range(nelx):
            for j in range(nely):
                sum = 0.0
                for k in range(max(i-rminf, 0), min(i+rminf+1, nelx)):
                    for l in range(max(j-rminf, 0), min(j+rminf+1, nely)):
                        fac = max(0, rmin - math.sqrt((i-k)**2+(j-l)**2));
                        sum = sum + fac;
                        dcn[j,i] = dcn[j,i] + fac*x[l,k]*dc[l,k];
            
                dcn[j,i] = dcn[j,i]/(x[j,i]*sum);

        return dcn

    # optimality criteria update
    def update(self, x, volfrac, dc):
        # ugly hardwired constants to fix later
        l1 = 0
        l2 = 100000
        lt = 1e-4
        move = 0.2
        xmin = 0.001
        xmax = 1.0

        nely, nelx = x.shape
        while (l2-l1 > lt):
            lmid = 0.5*(l2+l1)
            xnew = np.multiply(x, np.sqrt(-dc/lmid))
            xnew = np.maximum(xmin, np.maximum(x-move, np.minimum(xmax, np.minimum(x+move, xnew))));
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
        ke = e/(1-nu**2)* \
            np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
                       [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
                       [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
                       [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
                       [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
                       [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
                       [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
                       [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);

        return ke
