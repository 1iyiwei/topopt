'''
finite element solvers for the displacement from stiffness matrix and force
'''
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import cvxopt
import cvxopt.cholmod


class FESolver(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        freedofs = np.array(load.freedofs())
        nely, nelx = x.shape

        f = load.force()
        B_free = cvxopt.matrix(f[freedofs])

        k_free = self.gk_freedofs(load, x, ke, kmin, penal)
        k_free = cvxopt.spmatrix(k_free.data, k_free.row, k_free.col)

        u = np.zeros(load.dim*(nely+1)*(nelx+1))

        # setting up a fast cholesky decompositon solver
        cvxopt.cholmod.linsolve(k_free, B_free)
        u[freedofs] = np.array(B_free)[:, 0]
#        u[freedofs] = spsolve(k_free.tocsr(), f[freedofs])
        return u

    # global stiffness matrix
    def gk_freedofs(self, load, x, ke, penal):
        raise NotImplementedError


# coo_matrix should be faster
class CooFESolver(FESolver):
    def __init__(self, verbose=False):
        super().__init__(verbose)

    def gk_freedofs(self, load, x, ke, kmin, penal):
        nelx, nely = load.shape()
        freedofs = load.freedofs()

        edof, x_list, y_list = load.edof(nelx, nely)

        #  SIMP - Ee(xe) = Emin + x^p (E0-Emin)
        kd = x.T.reshape(nelx*nely, 1, 1) ** penal  # knockdown factor
        value_list = ((np.tile(kmin, (nelx*nely, 1, 1)) + np.tile(ke-kmin, (nelx*nely, 1, 1))*kd)).flatten()

        # coo_matrix sums duplicated entries and sipmlyies slicing
        dof = load.dim*(nelx+1)*(nely+1)
        k = coo_matrix((value_list, (y_list, x_list)), shape=(dof, dof)).tocsc()
        k_free = k[freedofs, :][:, freedofs].tocoo()
        return k_free

