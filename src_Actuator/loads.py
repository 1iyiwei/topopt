"""
This file containts the Load class that allows the generation of an object that
contains geometric, mesh, loads and boundary conditions that belong to the
load case. This version of the code is meant for local compliant maximization.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""

import numpy as np


class Load(object):
    """
    Load parent class that contains the basic functions used in all load cases.
    This class and its children do cantain information about the load case
    conciderd in the optimisation. The load case consists of the mesh, the
    loads, and the boundaries conditions. The class is constructed such that
    new load cases can be generated simply by adding a child and changing the
    function related to the geometry, loads and boundaries.

    Atributes
    -------
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    young : float
        Youngs modulus of the materias.
    Emin : float
        Artifical Youngs modulus of the material to ensure a stable FEA.
        It is used in the SIMP based material model.
    poisson : float
        Poisson ration of the material.
    dim : int
        Amount of dimensions conciderd in the problem, set at 2.

    Methods
    ------
    node(elx, ely)
        Returns the topleft node number of the element.
    nodes(elx, ely)
        Returns all nodes of the element.
    edof()
        Generats an array with the possitions af all degrees of freedom that
        belong to all elements.
    lk(E, nu)
        Calculates the element stiffness matrix.
    force()
        Returns an 2-D column array, the force vector of the loading condition.
    displaceloc()
        Returns a 2-D column array, a vector, l, locating the objective node,
        such that u·l = u_oute
    alldofs()
        Returns a list with all degrees of freedom.
    fixdofs()
        Returns a list with indices that are fixed by the boundary conditions.
    freedofs()
        Returns a list of arr indices that are not fixede
    passive()
        Retuns three lists containing the location and magnitude of fixed
        density valuese
    """
    def __init__(self, nelx, nely, young, Emin, poisson, ext_stiff):
        self.nelx = nelx
        self.nely = nely
        self.young = young
        self.Emin = Emin
        self.poisson = poisson
        self.dim = 2
        self.ext_stiff = ext_stiff

    # compute 1D index from 2D position for node (boundary of element)
    def node(self, elx, ely):
        """
        Calculates the topleft node number of the requested element

        Parameters
        ---------
        elx : int
            X position of the conciderd element.
        ely : int
            Y position of the conciderd element.

        Returns
        -------
        topleft : int
            The node number of the top left node
        """
        return (self.nely+1)*elx + ely

    # compute the 4 boundary nodes of an element
    def nodes(self, elx, ely):
        """
        Calculates all node numbers of the requested element

        Parameters
        ---------
        elx : int
            X position of the conciderd element.
        ely : int
            Y position of the conciderd element.

        Returns
        -------
        n1 : int
            The node number of the top left node.
        n2 : int
            The node number of the top right node.
        n3 : int
            The node number of the bottom right node.
        n4 : int
            The node number of the bottom left node.
        """
        n1 = self.node(elx,     ely    ) 
        n2 = self.node(elx + 1, ely,   ) 
        n3 = self.node(elx + 1, ely + 1)
        n4 = self.node(elx,     ely + 1)
        return n1, n2, n3, n4

    # edof that returns an array
    def edof(self):
        """
        Generates an array with the position of the nodes of each element in
        the global stiffness matrix.

        Results
        ------
        edof : 2-D array size(nelx*nely, 8)
            The list with all elements and their degree of freedom numbers.
        x_list : 1-D array len(nelx*nely*8*8)
            The list with the x indices of all ellements to be inserted into
            the global stiffniss matrix.
        y_list : 1-D array len(nelx*nely*8*8)
            The list with the y indices of all ellements to be inserted into
            the global stiffniss matrix.
        """
        # Creating list with element numbers
        elx = np.repeat(range(self.nelx), self.nely).reshape((self.nelx*self.nely, 1))  # x position of element
        ely = np.tile(range(self.nely), self.nelx).reshape((self.nelx*self.nely, 1))  # y position of element

        n1, n2, n3, n4 = self.nodes(elx, ely)
        edof = np.array([self.dim*n1,self.dim*n1+1, self.dim*n2,self.dim*n2+1,
                         self.dim*n3,self.dim*n3+1, self.dim*n4,self.dim*n4+1])
        edof = edof.T[0]

        x_list = np.repeat(edof, 8)  # flat list pointer of each node in an element
        y_list = np.tile(edof, 8).flatten()  # flat list pointer of each node in element
        return edof, x_list, y_list

    # element (local) stiffness matrix
    def lk(self, E, nu):
        """
        Calculates the local siffness matrix depending on E and nu.

        Parameters
        ---------
        E : float
            Youngs modulus of the material.
        nu : float
            Poisson ratio of the material.

        Returns
        -------
        ke : 2-D array size(8, 8)
            Local stiffness matrix.
        """
        k = np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
        ke = E/(1-nu**2) * \
            np.array([[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
                      [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
                      [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
                      [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
                      [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
                      [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
                      [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
                      [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
        return ke

    def force(self):
        """
        Returns an 1D array, the force vector of the loading condition.

        Returns
        -------
        f : 1-D column array length covering all degrees of freedom
            Empy force vector.
        """
        return np.zeros((self.dim*(self.nely+1)*(self.nelx+1), 1))

    def displaceloc(self):
        """
        Returns a zero vector, there is supposed to be an 1 implemented at the
        index where displacment output should be maximised, such that
        u·l = u_out

        Returns
        -------
        l : 1-D column array length covering all degrees of freedom
            Empty for the governing class.
        """
        return np.zeros((self.dim*(self.nely+1)*(self.nelx+1), 1))

    def alldofs(self):
        """
        Returns a list with all degrees of freedom.

        Returns
        -------
        all : 1-D list
            List with numbers from 0 to the maximum degree of freedom number.
        """
        return [x for x in range(self.dim*(self.nely+1)*(self.nelx+1))]

    def fixdofs(self):
        """
        Returns a list with indices that are fixed by the boundary conditions.

        Returns
        -------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom. This list is
            empty in this parrent class.
        """
        return []

    def freedofs(self):
        """
        Returns a list of arr indices that are not fixed

        Returns
        -------
        free : 1-D list
            List containing all elemens of alldogs except those that appear in
            the freedofs list.
        """
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        """
        Retuns three lists containing the location and magnitude of fixed
        density values


        Returns
        ------
        elx : 1-D list
            X coordinates of all passive elements, empty for the parrent class.
        ely : 1-D list
            Y ccordinates of all passive elements, empty for the parrent class.
        values : 1-D list
            Density values of all passive elements, empty for the parrent class.
        """
        return [], [], []


# example loading scenario, inverter with horizontal mirror axis
class Inverter(Load):
    """
    This child of the Loads class represents a top half of the symetric inverter
    design used for MEMS actuators. It contains an positive horizontal force at
    the bottom left corner which causes a negative displacement at the bottom
    right corner.

    Methods
    --------
    No methods are added compared to the parrent class. Only the force,
    displacementloc and fixdof equations are changed to contain the propper
    values for the boundary conditions and optimisation objective.
    """
    def __init__(self, nelx, nely, young, Emin, poisson, ext_stiff):
        super().__init__(nelx, nely, young, Emin, poisson, ext_stiff)

    def force(self):
        """
        The force vector containts a load in positive x direction at the bottom
        left corner node.
        
        Returns
        ------
        f : 1-D column array length covering all degrees of freedom
            Value of 1 at the index related to the bottom left node.           
        """
        f = super().force()
        f[self.dim*self.nely] = 1.0
        return f

    def displaceloc(self):
        """
        The maximisation should occur in negative x direction at the bottom
        right corner. Positive movement is thus in negative x direction.

        Returns
        -------
        l : 1-D column array length covering all degrees of freedom
            Value of -1 at the index related to the bottom right node.
        """
        dout = super().displaceloc()
        dout[-2] = -1
        return dout

    def fixdofs(self):
        """
        The boundary conditions of this problem fixes the bottom of the desing
        space in y direction (due to symetry). While the topleft element is
        fixed in both x and y direction on the left side.

        Returns
        -------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom.
        """
        # bottom fixed to the wall (y direction only), top left corner fixed to a point
        fix1 = np.arange(self.dim*self.nely + 1, self.dim*(self.nelx+1)*(self.nely+1), self.dim*(self.nely+1))
        return (np.hstack((fix1, [0, 1, 2, 3]))).tolist()
