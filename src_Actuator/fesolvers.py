"""
Finite element solvers for the displacement from stiffness matrix, force and
adjoin vector. This version of the code is meant for local compliant
maximization.

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

# Imporning linear algabla conjugate gradient solver
from scipy.sparse.linalg import cg
from scipy.sparse import diags


# coo_matrix should be faster
class FESolver(object):
    """
    This parent FEA class can only assemble the global stiffness matrix and
    exclude all fixed degrees of freedom from it. This stiffness csc-sparse
    stiffness matrix is assembled in the gk_freedof method. This
    class solves the FE problem with a sparse LU-solver based upon umfpack.
    This solver is slow and inefficient. It is however more robust.

    For this local compliance (actuator) maximization this solver solves two
    problems, the equilibrium and the adjoint problem which will be
    required to compute the gradients.

    Parameters
    ----------
    verbose : bool, optional
        False if the FEA should not print updates

    Attributes
    ----------
    verbose : bool
        False if the FEA should not print updates.
    """
    def __init__(self, verbose=False):
        self.verbose = verbose

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        """
        FE solver based upon the sparse SciPy solver that uses umfpack.

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considered for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffness matrix.
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

    # sparce stiffness matrix assembly
    def gk_freedofs(self, load, x, ke, kmin, penal):
        """
        Generates the global stiffness matrix with deleted fixed degrees of
        freedom. It generates a list with stiffness values and their x and y
        indices in the global stiffness matrix. Some combination of x and y
        appear multiple times as the degree of freedom might appear in multiple
        elements of the FEA. The SciPy coo_matrix function adds them up at the
        background. At the location of the force introduction and displacement
        output an external stiffness is added due to stability reasons.

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considered for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffness matrix.
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
        k[loc_force, loc_force] = load.ext_stiff + np.float64(k[loc_force, loc_force])
        k[loc_actuator, loc_actuator] = load.ext_stiff + np.float64(k[loc_actuator, loc_actuator])

        # selecting only the free directions of the siffness matirx
        k = k[freedofs, :][:, freedofs]

        return k


class CvxFEA(FESolver):
    """
    This parent FEA class can assemble the global stiffness matrix and solve
    the FE problem with a Supernodal Sparse Cholesky Factorization. It solves
    for both the equilibrium and adjoin problems.

    Attributes
    ----------
    verbose : bool
        False if the FEA should not print updates.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        """
        FE solver based upon a Supernodal Sparse Cholesky Factorization. It
        requires the installation of the cvx module. It solves both the FEA
        equilibrium and adjoint problems. [1]_

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffness matrix.
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


class CGFEA(FESolver):
    """
    This parent FEA class can assemble the global stiffness matrix and solve
    the FE problem with a sparse solver based upon a preconditioned conjugate
    gradient solver. The preconditioning is based upon the inverse of the
    diagonal of the stiffness matrix.

    Recommendations

    - Make the tolerance change over the iterations, low accuracy is
      required for first iteration, more accuracy for the later ones.
    - Add more advanced preconditioned.
    - Add gpu acceleration.

    Attributes
    ----------
    verbose : bool
        False if the FEA should not print updates.
    ufree_old : array len(freedofs)
        Displacement field of previous iteration.
    lambafree_old : array len(freedofs)
        Ajoint equation result of previous iteration.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)
        self.ufree_old = None
        self.lambafree_old = None

    # finite element computation for displacement
    def displace(self, load, x, ke, kmin, penal):
        """
        FE solver based upon the sparse SciPy solver that uses a preconditioned
        conjugate gradient solver, preconditioning is based upon the inverse
        of the diagonal of the stiffness matrix. Currently the relative
        tolerance is hardcoded as 1e-5. It solves both the equilibrium and
        adjoint problems.

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considered for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        ke : 2-D array size(8, 8)
            Local fully dense stiffness matrix.
        kmin : 2-D array size(8, 8)
            Local stiffness matrix for an empty element.
        penal : float
            Material model penalisation (SIMP).

        Returns
        -------
        u : 1-D array len(max(edof)+1)
            Displacement of all degrees of freedom
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.
        """
        freedofs = np.array(load.freedofs())
        nely, nelx = x.shape

        f_free = load.force()[freedofs]
        l_free = load.displaceloc()[freedofs]
        k_free = self.gk_freedofs(load, x, ke, kmin, penal)

        # Preconditioning
        L = diags(1/k_free.diagonal())

        # solving the system f = Ku with a cg implementation
        u = np.zeros((load.dim*(nely+1)*(nelx+1), 1))
        u[freedofs, 0], info1 = cg(k_free, f_free, x0=self.ufree_old, tol=1e-5, M=L)

        # solving adjoint problem l = Klamba with cg
        lamba = np.zeros((load.dim*(nely+1)*(nelx+1), 1))
        lamba[freedofs, 0], info2 = cg(k_free, l_free, x0=self.lambafree_old, tol=1e-5, M=L)

        # update uold
        self.ufree_old = u[freedofs]
        self.lamdafree_old = lamba[freedofs]

        if self.verbose is True:
            if info1 > 0:
                print('Convergence tolerance FEA not achieved after ', info1, ' itrations')
            if info1 < 0:
                print('Illegal input or breakdown FEA', info1)
            if info2 > 0:
                print('Convergence tolerance adjoint problem not achieved after ', info2, ' itrations')
            if info2 < 0:
                print('Illegal input or breakdown adjoint problem', info2)

        return u, lamba
