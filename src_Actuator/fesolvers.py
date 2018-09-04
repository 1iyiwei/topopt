"""
Finite element solvers for the displacement from stiffness matrix and force
vector. This version of the code is meant for local compliant maximization.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""

# Importing modules for parent class
import numpy as np
from scipy.sparse import coo_matrix

# Importing linear algabra solver for the CvxFEA class
import cvxopt
import cvxopt.cholmod

# Importing linear algabra solver for the SciPyFEA class
from scipy.sparse.linalg import spsolve


# coo_matrix should be faster
class CSCStiffnessMatrix(object):
    """
    This parent FEA class can only assemble the global stiffness matrix and
    exclude all fixed degrees of freedom from it. This function, gk_freedofs
    is used in all FEA solvers classes. The displace function is not
    implemented in this parrent class as it does not contain a solver for the
    linear problem.

    Atributes
    --------
    verbose : bool
        False if the FEA should not print updates.

    Methods
    -------
    displace(load, x, ke, kmin, penal)
        This function is not implemented, see child classes for implemetations
        of this function.
    gk_freedogs(self, load, x, ke, kmin, penal)
        Generates the global stiffness matrix with deleted fixed degrees of
        freedom. This includes adding the external stiffness to the load
        introduction and displacement output.
    """
    def __init__(self, verbose=False):
        self.verbose = verbose

    def displace(self, load, x, ke, kmin, penal):
        raise NotImplementedError

    def gk_freedofs(self, load, x, ke, kmin, penal):
        """
        Generates the global stiffness matrix with deleted fixed degrees of
        freedom. It generates a list with stiffness values and their x and y
        indices in the global stiffness matrix. Some combination of x and y
        appear multiple times as the degree of freedom might apear in multiple
        elements of the FEA. The SciPy coo_matrix function adds them up at the
        background. At the location of the force introduction and displacement
        output an external stiffness is added due to stability reasons.

        Parameters
        --------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffnes matrix.
        kmin : 2-D array size(8, 8)
            Local stiffness matrix for an empty element.
        penal : float
            Material model penalisation (SIMP).

        Returns
        -------
        k : 2-D sparse csc matrix
            Global stiffness matrix without fixed degrees of freedom.
        """
        freedofs = np.array(load.freedofs())
        nelx = load.nelx
        nely = load.nely

        edof, x_list, y_list = load.edof()

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
    """
    This parent FEA class can assemble the global stiffness matrix and solve
    the FE problem with a Supernodal Sparse Cholesky Factorization. It solves
    for both the equalibrium and adjoint problem.

    Atributes
    --------
    verbose : bool
        False if the FEA should not print updates.

    Methods
    -------
    displace(load, x, ke, kmin, penal)
        FE solver based upon a Supernodal Sparse Cholesky Factorization.
    gk_freedogs(self, load, x, ke, kmin, penal)
        Generates the global stiffness matrix with deleted fixed degrees of
        freedom. Function inherented from parent.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        """
        FE solver based upon a Supernodal Sparse Cholesky Factorization. It
        requires the instalation of the cvx module. It solves both the FEA
        equalibrium and adjoint problems. [1]_

        Parameters
        -------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffnes matrix.
        kmin : 2-D array size(8, 8)
            Local stiffness matrix for an empty element.
        penal : float
            Material model penalisation (SIMP).

        Returns
        -------
        u : 1-D column array shape(max(edof), 1)
            The displacement vector.
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.

        References
        ---------
        .. [1] Y. Chen, T. A. Davis, W. W. Hager, S. Rajamanickam, "Algorithm
            887: CHOLMOD, Supernodal Sparse Cholesky Factorization and
            Update/Downdate", ACM Transactions on Mathematical Software, 35(3),
            22:1-22:14, 2008.
        """
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
    """
    This parent FEA class can assemble the global stiffness matrix and solve
    the FE problem with a sparse solver based upon umfpack. This solver is
    slowen than the CvxFEA solver. It is however more robust. It solves
    for both the equalibrium and adjoint problem.

    Atributes
    --------
    verbose : bool
        False if the FEA should not print updates.

    Methods
    -------
    displace(load, x, ke, kmin, penal)
        FE solver based upon a SciPy sparse sysems solver that uses umfpack.
    gk_freedogs(self, load, x, ke, kmin, penal)
        Generates the global stiffness matrix with deleted fixed degrees of
        freedom. Function inherented from parent.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        """
        FE solver based upon the sparse SciPy solver that uses umfpack.

        Parameters
        -------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffnes matrix.
        kmin : 2-D array size(8, 8)
            Local stiffness matrix for an empty element.
        penal : float
            Material model penalisation (SIMP).

        Returns
        -------
        u : 1-D column array shape(max(edof), 1)
            The displacement vector.
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.
        """
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
        u[freedofs] = res[:, 0].reshape((len(freedofs), 1))
        lamba[freedofs] = res[:, 1].reshape((len(freedofs), 1))

        return u, lamba