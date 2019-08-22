"""
Finite element solvers for the displacement from stiffness matrix, force and
adjoint vector. This version of the code is meant for the fatigue crack growth
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
    exclude all fixed degrees of freedom from it. This function, gk_freedofs
    is used in all FEA solvers classes. The displace function is not
    implemented in this parrent class as it does not contain a solver for the
    linear problem.

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

    # finite element computation for displacement and adjoint
    def displace(self, load, x, penal, length):
        """
        FE solver based upon the sparse SciPy solver that uses umfpack.

        Parameters
        -------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        penal : float
            Material model penalisation (SIMP).
        length : int
            Length of the current crack conciderd.

        Returns
        -------
        u : 1-D column array shape(max(edof), 1)
            The displacement vector.
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.
        """
        freedofs = np.array(load.freedofs(length))
        nely, nelx = x.shape

        f = load.force(length)
        l = load.kiloc()

        f_free = np.hstack((f[freedofs], l[freedofs]))
        k_free = self.gk_freedofs(load, x, penal, length)

        # solving the system f = Ku with scipy
        num_dof = load.num_dofs
        u = np.zeros((num_dof, 1))
        lamba = np.zeros((num_dof, 1))

        res = spsolve(k_free, f_free)
        u[freedofs] = res[:, 0].reshape(len(freedofs), 1)
        lamba[freedofs] = res[:, 1].reshape(len(freedofs), 1)

        return u, lamba

    def gk_freedofs(self, load, x, penal, length):
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
        penal : float
            Material model penalisation (SIMP).
        length : int
            Length of the current crack conciderd.

        Returns
        -------
        k : 2-D sparse csc matrix
            Global stiffness matrix without fixed degrees of freedom.
        """
        # select propper stiffness lists from dictionary
        ke = load.k_list[str(length)]
        kmin = load.kmin_list[str(length)]
        x_list = load.x_list[str(length)]
        y_list = load.y_list[str(length)]

        freedofs = np.array(load.freedofs(length))
        nelx = load.nelx
        nely = load.nely

        #  SIMP - Ee(xe) = Emin + x^p (E-Emin)
        kd = x.T.reshape(nelx*nely) ** penal  # knockdown factor
        value_list = [kmini + kdi*(kei - kmini) for kei, kmini, kdi in zip(ke, kmin, kd)]
        value_list = [item for sublist in value_list for subsublist in sublist for item in subsublist]
        value_list = np.array(value_list)

        # coo_matrix sums duplicated entries and sipmlyies slicing
        dof = load.num_dofs
        k = coo_matrix((value_list, (y_list, x_list)), shape=(dof, dof)).tocsc()

        # adding external spring stiffness to load and actuator locations
        loc_force = np.where(load.force(length) != 0)[0]
        loc_ki = np.where(load.kiloc() != 0)[0]
        k[loc_force, loc_force] = load.ext_stiff + np.float64(k[loc_force, loc_force])
        k[loc_ki, loc_ki] = load.ext_stiff + np.float64(k[loc_ki, loc_ki])

        # selecting only the free directions of the siffness matirx
        k = k[freedofs, :][:, freedofs]

        return k


class CvxFEA(FESolver):
    """
    This parent FEA class can assemble the global stiffness matrix and solve
    the FE problem with a Supernodal Sparse Cholesky Factorization. It solves
    for both the equalibrium and adjoint problem.

    Attributes
    ----------
    verbose : bool
        False if the FEA should not print updates.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)

    # finite element computation for displacement
    def displace(self, load, x, penal, length):
        """
        FE solver based upon a Supernodal Sparse Cholesky Factorization. It
        requires the instalation of the cvx module. It solves both the FEA
        equalibrium and adjoint problems. [1]_

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        penal : float
            Material model penalisation (SIMP).
        length : int
            Length of the current crack conciderd.


        Returns
        -------
        u : 1-D column array shape(max(edof), 1)
            The displacement vector.
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.

        References
        ----------
        .. [1] Y. Chen, T. A. Davis, W. W. Hager, S. Rajamanickam, "Algorithm
            887: CHOLMOD, Supernodal Sparse Cholesky Factorization and
            Update/Downdate", ACM Transactions on Mathematical Software, 35(3),
            22:1-22:14, 2008.
        """
        freedofs = np.array(load.freedofs(length))
        nely, nelx = x.shape

        f = load.force(length)
        l = load.kiloc()
        B_free = cvxopt.matrix(np.hstack((f[freedofs], l[freedofs])))

        k_free = self.gk_freedofs(load, x, penal, length).tocoo()
        k_free = cvxopt.spmatrix(k_free.data, k_free.row, k_free.col)

        num_dof = load.num_dofs
        u = np.zeros((num_dof, 1))
        lamba = np.zeros((num_dof, 1))

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

    Recomendations

    - Make the tolerance change over the iterations, low accuracy is
      required for first itteration, more accuracy for the later ones.
    - Add more advanced preconditioner.
    - Add gpu accerelation.

    Attributes
    ----------
    verbose : bool
        False if the FEA should not print updates.
    ufree_old : dict
        Displacement field of previous iteration for every crack length,
        the keys are the related cracklengths.
    lambafree_old : array len(freedofs)
        Ajoint equation result of previos iteration for every crack length,
        the keys are the related cracklengths.
    """
    def __init__(self, verbose=False):
        super().__init__(verbose)
        self.keys = []
        self.ufree_old = {}
        self.lambafree_old = {}

    # finite element computation for displacement
    def displace(self, load, x, penal, length):
        """
        FE solver based upon the sparse SciPy solver that uses a preconditioned
        conjugate gradient solver, preconditioning is based upon the inverse
        of the diagonal of the stiffness matrix. Currently the relative
        tolerance is hardcoded as 1e-3.

        Parameters
        ----------
        load : object, child of the Loads class
            The loadcase(s) considerd for this optimisation problem.
        x : 2-D array size(nely, nelx)
            Current density distribution.
        penal : float
            Material model penalisation (SIMP).
        length : int
            Length of the current crack conciderd.

        Returns
        -------
        u : 1-D array len(max(edof)+1)
            Displacement of all degrees of freedom
        lamba : 1-D column array shape(max(edof), 1)
            Adjoint equation solution.
        """
        freedofs = np.array(load.freedofs(length))
        nely, nelx = x.shape
        num_dof = load.num_dofs

        f_free = load.force(length)[freedofs]
        l_free = load.kiloc()[freedofs]
        k_free = self.gk_freedofs(load, x, penal, length)

        # Preconditioning
        L = diags(1/k_free.diagonal())

        if str(length) in self.ufree_old:
            u0 = self.ufree_old[str(length)]
            lamba0 = self.lambafree_old[str(length)]
        else:
            u0 = None
            lamba0 = None

        # solving the system f = Ku with a cg implementation
        u = np.zeros((num_dof, 1))
        u[freedofs, 0], info1 = cg(k_free, f_free, x0=u0, tol=1e-3, M=L)

        # solving adjoint problem l = Klamba with cg
        lamba = np.zeros((num_dof, 1))
        lamba[freedofs, 0], info2 = cg(k_free, l_free, x0=lamba0, tol=1e-3, M=L)

        # update uold and lamdaold
        self.ufree_old[str(length)] = u[freedofs]
        self.lambafree_old[str(length)] = lamba[freedofs]

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
