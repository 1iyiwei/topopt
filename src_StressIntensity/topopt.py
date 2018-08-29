'''
topology optimization

2D only for now
'''
import csv
import numpy as np
from math import *
from mma import mma
from scipy.ndimage import convolve


class Topopt(object):
    '''
    young: young's modulus
    poisson: poisson ratio
    '''
    def __init__(self, fesolver, load, constraint, young=1, poisson=0.3, verbose=False):
        self.fesolver = fesolver
        self.young = young
        self.poisson = poisson
        self.dim = 2
        self.verbose = verbose
        Emin = constraint.Emin
        self.ke = self.lk(load, young, poisson)
        self.kmin = self.lk(load, Emin, poisson)

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

        while (change >= delta) and (itr < loopy):
            itr = itr + 1
            x, change, uout, xold1, xold2, low, upp = self.iter(load, constraint, x, penal, rmin, filt, itr, xold1, xold2, low, upp)

            if self.verbose: print('It.: {0:4d},  Obj.: {1:8.4f},  ch.: {2:0.3f}'.format(itr, uout, change), flush=True)
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
        ke = self.ke
        kmin = self.kmin

        # applying the density filter if required
        xf = self.densityfilt(x, rmin, filt)

        # Repair constant values after filtering of round off errors ect
        xlist, ylist, values = load.passive()
        x[ylist, xlist] = values

        # displacement via FEA
        u, lamba = self.fesolver.displace(load, xf, ke, kmin, penal)

        # compliance, its derivative
        uout, duout = self.disp(load, xf, u, lamba, ke, penal)

        # applying the sensitvity filter if required
        duoutf = self.sensitivityfilt(xf, rmin, duout, filt)

        # Prepairing MMA update scheme
        m = 1  # amount of constraint functions
        n = load.nelx*load.nely  # amount of elements
        xmin = constraint.xmin(load, x)  # vector with min element density
        xmax = constraint.xmax(load, x)  # vector with max element density
        x = x.flatten()
        duoutf = duoutf.flatten()
        volcon = constraint.current_volconstrain(xf, load)  # value of constraint function
        dvolcondx = constraint.volume_derivative(load)  # constraint derivative
        a0 = 1
        a = np.zeros((m))
        c_ = 1000*np.ones((m))
        d = a

        # Execute MMA update scheme
        xnew, low, upp = mma(m, n, itr, x, xmin, xmax, xold1, xold2, uout, duoutf, volcon, dvolcondx, low, upp, a0, a, c_, d)

        # Update variables
        xold2 = xold1
        xold1 = x
        x = xnew.reshape((load.nely, load.nelx))

        # What is the maximum change
        change = np.amax(abs(xnew - xold1))

        return x, change, uout, xold1, xold2, low, upp

    # updated compliance algorithm
    def disp(self, load, x, u, lamba, ke, penal):
        '''funcion calculates displacement at the actuator tip (dout) and
        the displacement density derivative'''

        # calculate displacment at actuator tip
        l = load.displaceloc()
        uout = np.dot(l.T, u)[0, 0]

        # calculating derivative
        nely, nelx = x.shape
        duout = np.zeros((nely, nelx))

        num = 0
        for elx in range(nelx):
            for ely in range(nely):
                ue = u[load.edof[num]]
                lambae = lamba[load.edof[num]]
                length = len(ue)
                unum = ue.reshape(length, 1)
                lambanum = lambae.reshape(length, 1)
                ce = np.dot(lambanum.T, np.dot(ke[num], unum))
                duout[ely, elx] = -penal * (x[ely, elx] ** (penal - 1)) * ce
                num += 1

        return uout, duout

    # sensitivity filter
    def densityfilt(self, x, rmin, filt):
        ''' filters the density with a radius of rmin if:

            >>> filt=='density'  '''
        if filt == 'density':
            rminf = floor(rmin)
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
            rminf = floor(rmin)
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
            dcn = np.divide(xdcn, x, out=np.zeros_like(xdcn), where=x!=0)

        else:
            dcn = dc

        return dcn

    # Importing the stiffness matrixes
    def import_stiffness(self, elementtype, E, nu):
        file = '../Stiffness Matrices/' + elementtype + '.csv'

        # import file
        with open(file, 'r') as f:
            reader = csv.reader(f)
            k = list(reader)
        f.close()

        # generate empy matrix with correct size
        K = np.zeros((len(k), len(k)))

        for i in range(len(k)):
            for j in range(len(k[i])):
                kij = k[i][j].strip()
                exec("Kij = {}".format(kij), locals(), globals())
                K[i, j] = Kij
        self.reset_Kij()
        return K

    # reset global Kij, fixes bug with exec??
    def reset_Kij(self):
        global Kij
        del Kij

    # list with all local stiffeness matrix of every element
    def lk(self, load, E, nu):
        k = []
        for elx in range(load.nelx):
            for ely in range(load.nely):
                if [elx, ely] in load.hoe:
                    index = [i for i, x in enumerate(load.hoe) if x == [elx, ely]][0]
                    if load.hoe_type[index] == '-1,-1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(-1:-1)', E, nu))
                    elif load.hoe_type[index] == '-1,1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(-1:1)', E, nu))
                    elif load.hoe_type[index] == '1,-1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(1:-1)', E, nu))
                    elif load.hoe_type[index] == '1,1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(1:1)', E, nu))
                    else:
                        raise Exception('The element type requested does not exist')
                else:
                    kk = np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
                    ke = E/(1-nu**2) * \
                        np.array([[kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6], kk[7]],
                                  [kk[1], kk[0], kk[7], kk[6], kk[5], kk[4], kk[3], kk[2]],
                                  [kk[2], kk[7], kk[0], kk[5], kk[6], kk[3], kk[4], kk[1]],
                                  [kk[3], kk[6], kk[5], kk[0], kk[7], kk[2], kk[1], kk[4]],
                                  [kk[4], kk[5], kk[6], kk[7], kk[0], kk[1], kk[2], kk[3]],
                                  [kk[5], kk[4], kk[3], kk[2], kk[1], kk[0], kk[7], kk[6]],
                                  [kk[6], kk[3], kk[4], kk[1], kk[2], kk[7], kk[0], kk[5]],
                                  [kk[7], kk[2], kk[1], kk[4], kk[3], kk[6], kk[5], kk[0]]])
                    k.append(ke)

        return k
