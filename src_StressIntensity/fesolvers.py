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

        #  SIMP - Ee(xe) = Emin + x^p (E-Emin)
        kd = x.T.reshape(nelx*nely) ** penal  # knockdown factor
        value_list = [kmini + kdi*(kei - kmini) for kei, kmini, kdi in zip(ke, kmin, kd)]
        value_list = [item for sublist in value_list for subsublist in sublist for item in subsublist]
        value_list = np.array(value_list)

        # coo_matrix sums duplicated entries and sipmlyies slicing
        dof = load.maxnode + 1
        k = coo_matrix((value_list, (load.y_list, load.x_list)), shape=(dof, dof)).tocsc()

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

        max_node = load.maxnode + 1
        u = np.zeros((max_node, 1))
        lamba = np.zeros((max_node, 1))

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
        max_node = load.maxnode + 1
        u = np.zeros((max_node, 1))
        lamba = np.zeros((max_node, 1))

        res = spsolve(k_free, f_free)
        u[freedofs] = res[:, 0].reshape(len(freedofs), 1)
        lamba[freedofs] = res[:, 1].reshape(len(freedofs), 1)

        return u, lamba