'''
finite element solvers for the displacement from stiffness matrix and force
'''

import numpy as np
# https://docs.scipy.org/doc/scipy-0.18.1/reference/sparse.html
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


class FESolver(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    # finite element computation for displacement
    def displace(self, load, x, ke, penal):
        f = load.force()

        fixdofs = np.array(load.fixdofs())
        freedofs = np.array(load.freedofs())

        nely, nelx = x.shape

        k_freedofs = self.gk_freedofs(load, x, ke, penal)

        u = np.zeros(load.dim*(nely+1)*(nelx+1))

        u[freedofs] = spsolve(k_freedofs, f[freedofs], permc_spec='COLAMD')
        u[fixdofs] = 0.0
        return u

    # global stiffness matrix
    def gk_freedofs(self, load, x, ke, penal):
        raise NotImplementedError


# coo_matrix should be faster
class CooFESolver(FESolver):
    def __init__(self, verbose=False):
        super().__init__(verbose)

    def gk_freedofs(self, load, x, ke, penal):
        nelx, nely = load.shape()

        edof, x_list, y_list = load.edof(nelx, nely)
        kd = x.T.reshape(nelx*nely, 1, 1) ** penal
        value_list = (np.tile(ke, (nelx*nely, 1, 1))*kd).flatten()

        # coo_matrix automatically sums duplicated entries, so it is handy
        dof = load.dim*(nelx+1)*(nely+1)
        k = coo_matrix((value_list, (y_list, x_list)), shape=(dof, dof)).tocsc()

        freedofs = load.freedofs()
        k_freedofs = k[freedofs,:][:,freedofs]
        return k_freedofs
