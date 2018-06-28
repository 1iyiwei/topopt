'''
topology optimization

2D only for now
'''
import numpy as np
import math
from mma import mma
from scipy.ndimage import convolve


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
    def layout(self, load, constraint, x, penal, rmin, delta, loopy, filt, history=False):
        # check if an existing filter was selected
        if filt != 'sensitivity' and filt != 'density':
            raise ValueError('No valid filter was selected, density of sensitivity are the only options')

        itr = 0  # number of loop iterations
        change = 1.0  # maximum density change from prior iteration
        xold1 = constraint.volfrac*np.ones((load.nely*load.nelx))
        xold2 = constraint.volfrac*np.ones((load.nely*load.nelx))
        low = np.zeros((load.nely*load.nelx))
        upp = np.zeros((load.nely*load.nelx))

        if history:
            xf_history = [1-x]

        while (change > delta) and (itr < loopy):
            itr = itr + 1
            x, change, c, xold1, xold2, low, upp = self.iter(load, constraint, x, penal, rmin, filt, itr, xold1, xold2, low, upp)

            if self.verbose: print('It.: {0:4d},  Obj.: {1:8.2f},  ch.: {2:0.3f}'.format(itr, c, change), flush=True)
            if history:
                xf = self.densityfilt(x, rmin, filt)
                xf_history.append(1-xf)

        # the filtered density is the physical desity
        xf = self.densityfilt(x, rmin, filt)

        # done
        if history:
            return xf, xf_history
        else:
            return xf, itr

    # initialization
    def init(self, load, constraint):
        (nelx, nely) = load.shape()
        xlist, ylist, values = load.passive()

        x = np.ones((nely, nelx))*constraint.volume_frac()
        x[ylist, xlist] = values
        return x

    # iteration
    def iter(self, load, constraint, x, penal, rmin, filt, itr, xold1, xold2, low, upp):
        # element stiffness matrix
        Emin = constraint.Emin
        ke = self.lk(self.young, self.poisson)
        kmin = self.lk(Emin, self.poisson)

        # applying the density filter if required
        xf = self.densityfilt(x, rmin, filt)

        # displacement via FEA
        u = self.fesolver.displace(load, xf, ke, kmin, penal)

        # compliance, its derivative
        c, dc = self.comp(load, xf, u, ke, penal)

        # applying the sensitvity filter if required
        dcf = self.sensitivityfilt(xf, rmin, dc, filt)

        # Prepairing MMA update scheme
        m = 1  # amount of constraint functions
        n = load.nelx*load.nely  # amount of elements
        xmin = constraint.xmin(load, x)  # vector with min element density
        xmax = constraint.xmax(load, x)  # vector with max element density
        x = x.flatten()
        dcf = dcf.flatten()
        volcon = constraint.current_volconstrain(xf, load)  # value of constraint function
        dvolcondx = constraint.volume_derivative(load)  # constraint derivative
        a0 = 1
        a = np.zeros((m))
        c_ = 1000*np.ones((m))
        d = a

        # Execute MMA update scheme
        xnew, low, upp = mma(m, n, itr, x, xmin, xmax, xold1, xold2, c, dcf, volcon, dvolcondx, low, upp, a0, a, c_, d)

        # Update variables
        xold2 = xold1
        xold1 = x
        x = xnew.reshape((load.nely, load.nelx))

        # What is the maximum change
        change = np.amax(abs(xnew - xold1))

        return x, change, c, xold1, xold2, low, upp

    # updated compliance algorithm
    def comp(self, load, x, u, ke, penal):
        '''funcion calculates compliance and compliance density derivative'''

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

    # sensitivity filter
    def densityfilt(self, x, rmin, filt):
        ''' filters the density with a radius of rmin if:

            >>> filt=='density'  '''
        if filt == 'density':
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

            # apply convolution filter
            xf = convolve(x, kernel, mode='reflect')

        else:
            xf = x

        return xf

    # sensitivity filter
    def sensitivityfilt(self, x, rmin, dc, filt):
        ''' filters sensitivity with a radius of rmin if:

            >>> filt=='sensitivity'
        '''
        if filt == 'sensitivity':
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
            xdcn = convolve(xdc, kernel, mode='reflect')
            dcn = np.divide(xdcn, x, out=np.zeros_like(xdcn), where=x!=0)  # fix devision by 0

        else:
            dcn = dc

        return dcn

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
