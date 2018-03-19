'''
finite element solvers for the displacement from stiffness matrix and force
'''

import numpy as np
# https://docs.scipy.org/doc/scipy-0.18.1/reference/sparse.html
from scipy.sparse import coo_matrix, lil_matrix, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

class FESolver(object):
    def __init__(self, verbose = False):
        self.verbose = verbose

    # finite element computation for displacement
    def displace(self, load, x, ke, penal):
        f = load.force()

        fixdofs = np.array(load.fixdofs())
        freedofs = np.array(load.freedofs())

        nely, nelx = x.shape

        k_freedofs = self.gk_freedofs(load, x, ke, penal)

        u = np.zeros(load.dim*(nely+1)*(nelx+1));

        u[freedofs] = spsolve(k_freedofs, f[freedofs])  

        u[fixdofs] = 0.0

        return u

    # global stiffness matrix
    def gk_freedofs(self, load, x, ke, penal):
        raise NotImplementedError

# Using lil_matrix is quite slow
class LilFESolver(FESolver):
    
    def __init__(self, verbose = False):
        super().__init__(verbose)

    def gk_freedofs(self, load, x, ke, penal):
        nelx, nely = load.shape()
        dof = load.dim*(nelx+1)*(nely+1)

        k = lil_matrix((dof, dof))

        for elx in range(nelx):
            for ely in range(nely):
                sel = load.edof(elx, ely, nelx, nely)
                k[np.ix_(sel, sel)] += ke*(x[ely, elx]**penal);

        freedofs = np.array(load.freedofs())

        k_freedofs = k[np.ix_(freedofs, freedofs)].tocsc()

        return k_freedofs

# coo_matrix should be faster
class CooFESolver(FESolver):

    def __init__(self, verbose = False):
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
