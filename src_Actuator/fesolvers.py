'''
finite element solvers for the displacement from stiffness matrix and force
'''
import numpy as np
from scipy.sparse import coo_matrix

from scipy.sparse.linalg import spsolve

import cvxopt
import cvxopt.cholmod


# coo_matrix should be faster
class CSCStiffnessMatrix(object):
    def __init__(self, verbose=False):
        self.verbose = verbose

    def displace(self, load, x, ke, kmin, penal):
        raise NotImplementedError

    def gk_freedofs(self, load, x, ke, kmin, penal):
        freedofs = np.array(load.freedofs())
        nelx, nely = load.shape()

        edof, x_list, y_list = load.edof(nelx, nely)

        #  SIMP - Ee(xe) = Emin + x^p (E-Emin)
        kd = x.T.reshape(nelx*nely, 1, 1) ** penal  # knockdown factor
        value_list = ((np.tile(kmin, (nelx*nely, 1, 1)) + np.tile(ke-kmin, (nelx*nely, 1, 1))*kd)).flatten()

        # coo_matrix sums duplicated entries and sipmlyies slicing
        dof = load.dim*(nelx+1)*(nely+1)
        k = coo_matrix((value_list, (y_list, x_list)), shape=(dof, dof)).tocsc()

        # adding external spring stiffness to load and actuator locations
        loc_force = np.where(load.force() != 0)[0]
        loc_actuator = np.where(load.displaceloc() != 0)[0]
        loc = np.hstack((loc_force, loc_actuator))
        k[loc, loc] += load.ext_stiff*np.ones(len(loc))

        # selecting only the free directions of the siffness matirx
        k = k[freedofs, :][:, freedofs]

        return k


class CvxFEA(CSCStiffnessMatrix):
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        freedofs = np.array(load.freedofs())
        nely, nelx = x.shape

        f = load.force()
        l = load.displaceloc()
        B_free = cvxopt.matrix(np.hstack((f[freedofs], l[freedofs])))

        k_free = self.gk_freedofs(load, x, ke, kmin, penal).tocoo()
        k_free = cvxopt.spmatrix(k_free.data, k_free.row, k_free.col)

        u = np.zeros((load.dim*(nely+1)*(nelx+1), 1))
        lamba = np.zeros((load.dim*(nely+1)*(nelx+1), 1))

        # setting up a fast cholesky decompositon solver
        cvxopt.cholmod.linsolve(k_free, B_free)
        u[freedofs] = np.array(B_free[:, 0])
        lamba[freedofs] = np.array(B_free[:, 1])

        return u, lamba


class SciPyFEA(CSCStiffnessMatrix):
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        freedofs = np.array(load.freedofs())
        nely, nelx = x.shape

        f = load.force()
        l = load.displaceloc()

        f_free = np.hstack((f[freedofs], l[freedofs]))
        k_free = self.gk_freedofs(load, x, ke, kmin, penal)

        # solving the system f = Ku with scipy
        u = np.zeros((load.dim*(nely+1)*(nelx+1), 1))
        lamba = np.zeros((load.dim*(nely+1)*(nelx+1), 1))        
        res = spsolve(k_free, f_free)
        u[freedofs] = res[:, 0]
        lamba[freedofs] = res[:, 1]

        return u, lamba