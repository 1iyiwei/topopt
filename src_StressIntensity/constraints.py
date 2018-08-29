'''
constraints
'''
import numpy as np


class DensityConstraint(object):
    def __init__(self, load, move, volume_frac, density_min=0.0, density_max=1.0, Emin=1e-9):
        self.volfrac = volume_frac
        self.move = move
        self.density_min = density_min
        self.density_max = density_max
        v = 1/(load.nelx*load.nely*self.volume_frac())
        self.dvdx = v*np.ones((1, load.nely*load.nelx))
        self.Emin = Emin

    def volume_frac(self):
        return self.volfrac

    def volume_derivative(self, load):
        return self.dvdx

    def xmin(self, load, x):
        elx, ely, value = load.passive()
        xmin = self.density_min*np.ones((load.nely, load.nelx))
        xmin = np.maximum(xmin, x - self.move)
        x[ely, elx] = value
        return xmin.flatten()

    def xmax(self, load, x):
        elx, ely, value = load.passive()
        xmax = self.density_max*np.ones((load.nely, load.nelx))
        xmax = np.minimum(xmax, x + self.move)
        x[ely, elx] = value
        return xmax.flatten()

    def zero_stiffness(self):
        return self.Emin

    def current_volconstrain(self, x, load):
        cur_vol = np.sum(x)/(load.nelx*load.nely*self.volfrac) - 1
        return cur_vol
