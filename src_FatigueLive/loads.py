"""
This file containts the Load class that allows the generation of an object that
contains geometric, mesh, loads and boundary conditions that belong to the
load case. This version of the code is meant for the fatigue live maximization.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""

import numpy as np
from math import *
import csv

class Load(object):
    """
    Load parent class that contains the basic functions used in all load cases.
    This class and its children do cantain information about the load case
    conciderd in the optimisation. The load case consists of the mesh, the
    loads, and the boundaries conditions. The class is constructed such that
    new load cases can be generated simply by adding a child and changing the
    function related to the geometry, loads and boundaries.

    Parameters
    ---------
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
    ext_stiff : float
        Extra stiffness to be added to global stiffness matrix. Due to
        interactions with meganisms outside design domain.
    hoe : dict
        Dictionary with for every cracklength the x end y element locations
        that need to be enriched.

    Atributes
    -------
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    dim : int
        Amount of dimensions conciderd in the problem, set at 2.
    edof : dict
        Dictionary containing list with all elements and their degree of
        freedom numbers for all crack_lengtgs, str(length) is the key.
    x_list : dict
        Dictionary with a 1D list that contains the x indices of all degrees of
        freedom for all cracklengths, str(length) is the key.
    y_list : dict
        Dictionary with a 1D list that contains the y indices of all degrees of
        freedom for all cracklengths, str(length) is the key.
    num_dofs : int
        Amount of degrees of freedom.
    young : float
        Youngs modulus of the materias.
    Emin : float
        Artifical Youngs modulus of the material to ensure a stable FEA.
        It is used in the SIMP based material model.
    poisson : float
        Poisson ration of the material.
    k_list : dict
        Dictionary containing a list for every crack length, these lists
        contain the element stiffness matrices of full density for all
        elements, str(length) is the key.
    kmin_list : list len(nelx*nely)
        Dictionary containing a list for every crack length, these lists
        contain the epmy element stiffness matrices for all elements,
        str(length) is the key.
    ext_stiff : float
        Extra stiffness to be added to global stiffness matrix. Due to
        interactions with meganisms outside design domain.

    Methods
    ------
    node(elx, ely)
        Returns the topleft node number of the element.
    nodes(elx, ely)
        Returns all nodes of the element.
    edofcalc(hoe)
        Generats an array with the possitions af all degrees of freedom that
        belong to all elements.
    import_stiffness(elementtype, E, nu)
        Imports and reads stiffness matrices from .csv file
    reset_Kij()
        Bug fix that removes Kij as a global variable.
    lk(load, E, nu)
        Generates a list with all element stiffness matrices.
    force()
        Returns an 2-D column array, the force vector of the loading condition.
    kiloc()
        Returns a 2-D column array, a vector, l, locating the stress intensity
        node, such that u·l = KI
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
    def __init__(self, nelx, nely, young, Emin, poisson, ext_stiff, hoe):
        self.nelx = nelx
        self.nely = nely
        self.dim = 2
        self.young = young
        self.Emin = Emin
        self.poisson = poisson
        self.ext_stiff = ext_stiff
        self.edof = {}
        self.x_list = {}
        self.y_list = {}
        self.k_list = {}
        self.kmin_list = {}
        for length in hoe:
            edof, x_list, y_list, num_dofs = self.edofcalc(hoe[length])
            self.edof[length] = edof
            self.x_list[length] = x_list
            self.y_list[length] = y_list
            self.num_dofs = num_dofs
            k_list = self.lk(young, poisson, hoe[length])
            kmin_list = self.lk(Emin, poisson, hoe[length])
            self.k_list[length] = k_list
            self.kmin_list[length] = kmin_list

    # compute 1D index from 2D position for node (boundary of element)
    def node(self, elx, ely):
        """
        Calculates the topleft node number of the requested element. Does not
        toke Higher Order Elements in account.

        Parameters
        ---------
        elx : int
            X position of the conciderd element.
        ely : int
            Y position of the conciderd element.

        Returns
        -------
        topleft : int
            The node number of the top left node.
        """
        return (self.nely+1)*elx + ely

    # compute the 4 boundary nodes of an element
    def nodes(self, elx, ely):
        """
        Calculates all node numbers of the requested element. Does not take
        Higher Order Elements in account.

        Parameters
        ---------
        elx : int
            X position of the conciderd element.
        ely : int
            Y position of the conciderd element.

        Returns
        -------
        n0 : int
            The node number of the bottom left node.
        n1 : int
            The node number of the bottom right node.
        n2 : int
            The node number of the top left node.
        n3 : int
            The node number of the top right node.
        """
        n3 = self.node(elx,     ely    ) 
        n2 = self.node(elx + 1, ely    ) 
        n1 = self.node(elx + 1, ely + 1) 
        n0 = self.node(elx,     ely + 1)
        return n0, n1, n2, n3

    # edof that returns an array
    def edofcalc(self, hoe):
        """
        Generates an array with the position of the nodes of each element in
        the global stiffness matrix. This takes the Higher Order Elements in
        account.

        Parameters
        ----------
        hoe : list
            A list containing the x and y location of the higher order elemens
            for this crack length.

        Results
        ------
        edof : 2-D list size(nelx*nely, # degrees of freedom per element)
            The list with all elements and their degree of freedom numbers.
        x_list : 1-D array
            The list with the x indices of all ellements to be inserted into
            the global stiffniss matrix.
        y_list : 1-D array
            The list with the y indices of all ellements to be inserted into
            the global stiffniss matrix.
        num_dofs : int
            The amount of degrees of freedom.
        """
        # Creating list with element numbers
        elx = np.repeat(range(self.nelx), self.nely).reshape((self.nelx*self.nely, 1))  # x position of element
        ely = np.tile(range(self.nely), self.nelx).reshape((self.nelx*self.nely, 1))  # y position of element

        node_loc = np.array(self.nodes(elx, ely), dtype=np.int32)
        shape = list(np.shape(node_loc))
        shape[0] = shape[0]*2
        edof = np.zeros(shape, dtype=np.int32)
        edof[::2] = self.dim*node_loc
        edof[1::2] = self.dim*node_loc + 1
        edof = edof.T[0]

        # shift all elements such that proper space is create around higher order elements
        edofold = edof.copy()
        for i in range(len(hoe)):
            hoei = hoe[i]
            n0, n1, n2, n3 = self.nodes(hoei[0], hoei[1]) # dofs of hoe

            edof[np.where(edofold > self.dim*n3 + 1)] += self.dim * 4
            edof[np.where(edofold > self.dim*n0 + 1)] += self.dim * 2
            edof[np.where(edofold > self.dim*n2 + 1)] += self.dim * 2

            # check if higher order element exists to the left.
            if [hoei[0] - 1, hoei[1]] in hoe:
                edof[np.where(edofold > self.dim*n3 + 1)] -= self.dim * 2

            # check if higher order element exists to the top
            if [hoei[0], hoei[1] - 1] in hoe:
                edof[np.where(edofold > self.dim*n3 + 1)] -= self.dim * 2

        # add full version of edof of the higher order elements
        edofmax = np.max(edof)
        edof = edof.tolist()
        for i in range(len(hoe)):
            hoei = hoe[i]
            element = self.nely*hoei[0] + hoei[1]
            edof_ele = edof[element]  # old dofs of hoe

            n0 = edof_ele[0]
            n1 = edof_ele[0] + self.dim*1
            n2 = edof_ele[0] + self.dim*2
            n3 = edof_ele[2]
            n4 = edof_ele[2] - self.dim*1
            n5 = edof_ele[2] - self.dim*2
            n6 = edof_ele[4]
            n7 = edof_ele[0] - self.dim*3
            n8 = edof_ele[0] - self.dim*4
            n9 = edof_ele[6]
            n10 = edof_ele[0] - self.dim*2
            n11 = edof_ele[0] - self.dim*1
            kk = edofmax + 1

            edof[element] = [n0, n0+1, n1, n1+1, n2, n2+1, n3, n3+1, n4, n4+1,
                             n5, n5+1, n6, n6+1, n7, n7+1, n8, n8+1, n9, n9+1,
                             n10, n10+1, n11, n11+1, kk, kk+1]

        x_list = [[[i]*len(edofi) for i in edofi] for edofi in edof]
        x_list = [item for sublist in x_list for subsublist in sublist for item in subsublist]
        y_list = [edofi*len(edofi) for edofi in edof]
        y_list = [item for sublist in y_list for item in sublist]

        num_dofs = max(x_list) + 1

        return edof, np.array(x_list), np.array(y_list), num_dofs

    # Importing the stiffness matrixes
    def import_stiffness(self, elementtype, E, nu):
        """
        This function imports a matrix from a csv file that has variables to
        the material properties. The correct material properties are added.

        Parameters
        ----------
        elementtype : str
            Describes what .csv file should be used for the import.
        E : float
            Youngs modulus of the material.
        nu : float
            Poissons ratio of the material.

        Returns
        -------
        lk : array size(dofs, dofs)
            Element stiffness matrix
        """
        file = 'Stiffness Matrices/' + elementtype + '.csv'

        # import file
        with open(file, 'r') as f:
            reader = csv.reader(f)
            k = list(reader)
        f.close()

        # generate empy matrix with correct size
        lk = np.zeros((len(k), len(k)))

        for i in range(len(k)):
            for j in range(len(k[i])):
                kij = k[i][j].strip()
                exec("Kij = {}".format(kij), locals(), globals())
                lk[i, j] = Kij
        self.reset_Kij()
        return lk

    # reset global Kij, fixes bug with exec??
    def reset_Kij(self):
        """
        Resets the global variable Kij. This is neccesary as function
        import_stiffness will not clean up its local variables itself.
        """
        global Kij
        del Kij

    # list with all local stiffeness matrix of every element
    def lk(self, E, nu, hoe):
        """
        Generates a list with all element stiffness matrices. It differenciates
        between the element types used.

        Parameters
        ---------
        E : float
            Youngs modulus of the material.
        nu : float
            Poissons ratio of the material.

        Returns
        k : list len(nelx*nely)
            Returns a list with all local stiffness matrices.
        """
        k = []
        for elx in range(self.nelx):
            for ely in range(self.nely):
                if [elx, ely] in hoe:
                    index = [i for i, x in enumerate(hoe) if x == [elx, ely]][0]
                    if self.hoe_type[index] == '-1,-1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(-1:-1)', E, nu))
                    elif self.hoe_type[index] == '-1,1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(-1:1)', E, nu))
                    elif self.hoe_type[index] == '1,-1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(1:-1)', E, nu))
                    elif self.hoe_type[index] == '1,1':
                        k.append(self.import_stiffness('Stiffness_Cubic_PlaneStress_Enriched(1:1)', E, nu))
                    else:
                        raise Exception('The element type requested does not exist')
                else:
                    kk = np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
                    ke = E/(1-nu**2) * \
                        np.array([[kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6], kk[7]],
                                  [kk[1], kk[0], kk[7], kk[6], kk[5], kk[4], kk[3], kk[2]],
                                  [kk[2], kk[7], kk[0], kk[5], kk[6], kk[3], kk[4], kk[1]],
                                  [kk[3], kk[6], kk[5], kk[0], kk[7], kk[2], kk[1], kk[4]],
                                  [kk[4], kk[5], kk[6], kk[7], kk[0], kk[1], kk[2], kk[3]],
                                  [kk[5], kk[4], kk[3], kk[2], kk[1], kk[0], kk[7], kk[6]],
                                  [kk[6], kk[3], kk[4], kk[1], kk[2], kk[7], kk[0], kk[5]],
                                  [kk[7], kk[2], kk[1], kk[4], kk[3], kk[6], kk[5], kk[0]]])
                    k.append(ke)

        return k

    def force(self):
        """
        Returns an 1D array, the force vector of the loading condition.

        Returns
        -------
        f : 1-D column array length covering all degrees of freedom
            Empy force vector.
        """
        return np.zeros((self.num_dofs, 1))

    def kiloc(self):
        """
        The location of the stress intensity factor KI can be found at the
        second last index.

        Returns
        -------
        l : 1-D column array length covering all degrees of freedom
            Zeros except for the second last index.
        """
        l = np.zeros((self.num_dofs, 1))

        # higher order element
        l[-2] = 1
        return l

    def alldofs(self):
        """
        Returns a list with all degrees of freedom.

        Returns
        -------
        all : 1-D list
            List with numbers from 0 to the maximum degree of freedom number.
        """
        return [x for x in range(self.num_dofs)]

    def fixdofs(self, length_i):
        """
        Returns a list with indices that are fixed by the boundary conditions.

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom. This list is
            empty in this parrent class.
        """
        return []

    def freedofs(self, length_i):
        """
        Returns a list of arr indices that are not fixed

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        free : 1-D list
            List containing all elemens of alldogs except those that appear in
            the freedofs list.
        """
        return list(set(self.alldofs()) - set(self.fixdofs(length_i)))

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
        fix_ele : 1-D list
            List with all element numbers that are allowed to change.
        """
        return [], [], [], [x for x in range(self.nelx*self.nely)]


# example loading scenario, inverter with horizontal mirror axis
class EdgeCrack(Load):
    """
    This child class of Load class represents the symetric top half of an edge
    crack. The crack is positioned to the bottom left and propegates towards
    the right. Special elements are placed around the crack tip. The plate is
    subjected to a distributed tensile load (σ=1) on the top.

    For a perfectly flat plate analytical expressions for K_I are known. [1]_

    The stress intensity factors calculated can be be interperted in two ways:

    1. Without schaling. This means that all elements have a size of 2 length units.
    2. With schaling, comparison to reality should be based upon. ::

        Kreal = sigma_real * sqrt(real_crack_length) * Ksim/sqrt(2*crack_length)

    Parameters
    ---------
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    crack_length : array
        All crack lengs conciderd.
    young : float
        Youngs modulus of the materias.
    Emin : float
        Artifical Youngs modulus of the material to ensure a stable FEA.
        It is used in the SIMP based material model.
    poisson : float
        Poisson ration of the material.
    ext_stiff : float
        Extra stiffness to be added to global stiffness matrix. Due to
        interactions with meganisms outside design domain.

    Atributes
    ---------
    Two attributes are added with respect to the parrent class.

    crack_length : int
        Is the amount of elements that the crack is long.
    hoe : list len(2)
        List containing the x end y element locations that need to be enriched.

    Methods
    -------
    No methods are added compared to the parrent class. Only the force and
    fixdof are changed to contain the propper values for the boundary
    conditions.

    References
    ----------
    .. [1] Tada, H., Paris, P., & Irwin, G. (2000). "Part II 2.10-2.12 The
        Single Edge Notch Test Specimen", The stress analysis of cracks handbook
        (3rd ed.). New York: ASME Press, pp:52-54.
    """
    def __init__(self, nelx, nely, crack_length, young, Emin, poisson, ext_stiff):
        self.crack_length = crack_length
        hoe = {}
        for length in crack_length:
            hoe_i = [[length-1, nely-1], [length, nely-1]]
            hoe[str(length)] = hoe_i
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, nely, young, Emin, poisson, ext_stiff, hoe)

    def force(self, length_i):
        """
        The top of the design space is pulled upwards by 1MPa. This means that
        the nodal forces are 2 upwards, except for the top corners they have a
        load of 1 only.

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        f : 1-D column array length covering all degrees of freedom
            Force vector.
        """
        # select proper element from dictionary
        edof = self.edof[str(length_i)]

        f = super().force()
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        forceloc = [edof[i][5] for i in top_ele][:-1]
        f[forceloc] = 2
        f[1] = 1
        f[[edof[i][5] for i in top_ele][-1]] = 1
        return f

    def fixdofs(self, length_i):
        """
        The boundary conditions limmit y-translation at the bottom of the design
        space (due to symetry) and x-translations at the top (due to the clamps)

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        ------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom.
        """
        # select proper element from dictionary
        edof = self.edof[str(length_i)]

        # higher order element after crack
        ele = self.nely*(length_i+1)-1
        bottom = [1, 3, 5]
        fix = [edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(length_i+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [edof[i][1] for i in bot_ele] + [edof[bot_ele[-1]][3]]

        # top elements
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        # all elements in x direction fixed
        fix2 = [0] + [edof[i][4] for i in top_ele]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def passive(self):
        """
        Retuns three lists containing the location and magnitude of fixed
        density values. The elements around the crack tip are fixed at a
        density of one.

        Returns
        ------
        elx : 1-D list
            X coordinates of all passive elements, empty for the parrent class.
        ely : 1-D list
            Y ccordinates of all passive elements, empty for the parrent class.
        values : 1-D list
            Density values of all passive elements, empty for the parrent class.
        fix_ele : 1-D list
            List with all element numbers that are allowed to change.
        """
        elx = [x for x in range(np.min(self.crack_length)-1, np.max(self.crack_length)+1)]
        ely = [self.nely-1]*len(elx)
        values = [1]*len(elx)

        fixele = []
        for i in range(len(elx)):
            fixele.append(self.nelx*ely[i] + elx[i])
        free_ele = list(set(range(self.nelx*self.nely)) - set(fixele))

        return elx, ely, values, free_ele


class DoubleEdgeCrack(Load):
    """
    This child class of Load class represents the symetric top rigth quarter of
    an double edge crack plate. The crack is positioned to the bottom left and
    propegatestowards the right. Special elements are placed around the crack
    tip. The plate is subjected to a distributed tensile load (σ=1) on the top.

    For a perfectly flat plate analytical expressions for K_I are known. [1]_

    The stress intensity factors calculated can be be interperted in two ways:

    1. Without schaling. This means that all elements have a size of 2 length units.
    2. With schaling, comparison to reality should be based upon. ::

        Kreal = sigma_real * sqrt(real_crack_length) * Ksim/sqrt(2*crack_length)

    Parameters
    ---------
    nelx : int
        Number of elements in x direction.
    young : float
        Youngs modulus of the materias.
    Emin : float
        Artifical Youngs modulus of the material to ensure a stable FEA.
        It is used in the SIMP based material model.
    poisson : float
        Poisson ration of the material.
    ext_stiff : float
        Extra stiffness to be added to global stiffness matrix. Due to
        interactions with meganisms outside design domain.

    Atributes
    ---------
    Two attributes are added and one was changed with respect to the parrent
    class.

    nely : int
        Number of y elements, this is now a function of nelx.
    crack_length : int
        Is the amount of elements that the crack is long, this is a function of
        nelx.
    hoe : list len(2)
        List containing the x end y element locations that need to be enriched.

    Methods
    -------
    No methods are added compared to the parrent class. Only the force,
    fixdof and passive equations are changed to contain the
    propper values for the boundary conditions and
    passive elements.

    References
    ----------
    .. [1] Tada, H., Paris, P., & Irwin, G. (2000). "Part II 2.6-2.9a The
        Double Edge Notch Test Specimen", The stress analysis of cracks handbook
        (3rd ed.). New York: ASME Press, pp:46-51.
    """
    def __init__(self, nelx, young, Emin, poisson, ext_stiff):
        nelx = nelx
        nely = 4*nelx
        self.crack_length = [int(nelx/5*2)]
        hoe = {}
        for length in self.crack_length:
            hoe_i = [[length-1, nely-1], [length, nely-1]]
            hoe[str(length)] = hoe_i
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, nely, young, Emin, poisson, ext_stiff, hoe)

    def force(self, length_i):
        """
        The top of the design space is pulled upwards by 1MPa. This means that
        the nodal forces are 2 upwards, except for the top left corner has
        a load of 1 only.

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        f : 1-D column array length covering all degrees of freedom
            Force vector
        """
        # select proper element from dictionary
        edof = self.edof[str(length_i)]

        f = super().force()
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        forceloc = [edof[i][5] for i in top_ele]
        f[forceloc] = 2
        f[1] = 1
        return f

    def fixdofs(self, length_i):
        """
        The right side is fixed in x direction (symetry around the y axis) while
        the bottom side is fixed in y direction (symetry around the x axis).

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom.
        """
        # select proper element from dictionary
        edof = self.edof[str(length_i)]
        
        # higher order element after crack
        ele = self.nely*(length_i+1)-1
        bottom = [1, 3, 5]
        fix = [edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(length_i+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [edof[i][1] for i in bot_ele] + [edof[bot_ele[-1]][3]]

        # right side elements
        start_ele = (self.nelx-1)*self.nely
        right_ele = [x for x in range(start_ele, self.nelx*self.nely, 1)]
        # all elements in x direction fixed
        fix2 = [edof[i][4] for i in right_ele] + [edof[bot_ele[-1]][2]]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def passive(self):
        """
        Retuns three lists containing the location and magnitude of fixed
        density values. The elements around the crack tip are fixed at a
        density of one.

        Returns
        ------
        elx : 1-D list
            X coordinates of all passive elements, empty for the parrent class.
        ely : 1-D list
            Y ccordinates of all passive elements, empty for the parrent class.
        values : 1-D list
            Density values of all passive elements, empty for the parrent class.
        fix_ele : 1-D list
            List with all element numbers that are allowed to change.
        """
        elx = [x for x in range(np.min(self.crack_length)-1, np.max(self.crack_length)+1)]
        ely = [self.nely-1]*len(elx)
        values = [1]*len(elx)

        fixele = []
        for i in range(len(elx)):
            fixele.append(self.nelx*ely[i] + elx[i])
        free_ele = list(set(range(self.nelx*self.nely)) - set(fixele))

        return elx, ely, values, free_ele


class CompactTension(Load):
    """
    This child class of Load class represents the symetric top half of an
    compact tension specimen. The crack is positioned to the bottom left and
    propegatestowards the right. Special elements are placed around the crack
    tip. The plate is subjected to upwards load of one. The design follows the
    ASTM standard. [1]_

    For a perfectly flat plate analytical expressions for K_I do exist. [1]_
    [2]_

    The stress intensity factors calculated is for the case where the element
    size is the real dimensions times two and a load of 1.

    Parameters
    ---------
    nelx : int
        Number of elements in x direction.
    crack_length : array
        An array containing all crack lengths conciderd.
    young : float
        Youngs modulus of the materias.
    Emin : float
        Artifical Youngs modulus of the material to ensure a stable FEA.
        It is used in the SIMP based material model.
    poisson : float
        Poisson ration of the material.
    ext_stiff : float
        Extra stiffness to be added to global stiffness matrix. Due to
        interactions with meganisms outside design domain.

    Atributes
    ---------
    Two attributes are added and one was changed with respect to the parrent
    class.

    nely : int
        Number of y elements, this is now a function of nelx.
    crack_length : array
        Is for all cracks conciderd the crack_length.

    Methods
    -------
    No methods are added compared to the parrent class. Only the force and
    fixdof are changed to contain the propper values for the boundary
    conditions.

    References
    ----------
    .. [1] ASTM Standard E647-15e1, “Standard Test Method for Measurement of
        Fatigue Crack Growth Rates,” ASTM Book of Standards, vol. 0.30.1, 2015.
    .. [2] Tada, H., Paris, P., & Irwin, G. (2000). "Part II 2.19-2.21 The
        Compact Tension Test Specimen", The stress analysis of cracks handbook
        (3rd ed.). New York: ASME Press, pp:61-63.
    """
    def __init__(self, nelx, crack_length, young, Emin, poisson, ext_stiff):
        nelx = nelx
        nely = int(np.round(nelx/1.25*1.2/2))
        self.crack_length = crack_length
        hoe = {}
        for length in crack_length:
            hoe_i = [[length-1, nely-1], [length, nely-1]]
            hoe[str(length)] = hoe_i
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, nely, young, Emin, poisson, ext_stiff, hoe)

    def force(self, length_i):
        """
        The ASTM standard requires the force to be located approx. 1/5 of nelx
        and at 0.195 * nely from the top.

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        f : 1-D column array length covering all degrees of freedom
            Force vector
        """
        # select correct edof list from dictionary
        edof = self.edof[str(length_i)]

        f = super().force()
        load_ele = int(self.nely/0.6*0.325)+int(self.nelx/5)*self.nely
        forceloc = edof[load_ele][3]
        f[forceloc] = 1
        return f

    def fixdofs(self, length_i):
        """
        The bottom of the design space is fixed in y direction (due to symetry
        around the x axis). While at the location that the load is introduced
        x translations are constraint.

        Parameters
        ----------
        length_i : int
            Length of the crack for the current mesh

        Returns
        -------
        fix : 1-D list
            List with all the numbers of fixed degrees of freedom.
        """
        # select correct edof list from dictionary
        edof = self.edof[str(length_i)]

        # higher order element after crack
        ele = self.nely*(length_i+1)-1
        bottom = [1, 3, 5]
        fix = [edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(length_i+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [edof[i][1] for i in bot_ele] + [edof[bot_ele[-1]][3]]

        # load introduction in y direction rigid
        load_ele = int(self.nely/0.6*0.325)+int(self.nelx/5)*self.nely
        fix2 = edof[load_ele][2]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def passive(self):
        """
        Retuns three lists containing the location and magnitude of fixed
        density values. The elements around the crack tip are fixed at a
        density of one.

        Returns
        ------
        elx : 1-D list
            X coordinates of all passive elements, empty for the parrent class.
        ely : 1-D list
            Y ccordinates of all passive elements, empty for the parrent class.
        values : 1-D list
            Density values of all passive elements, empty for the parrent class.
        fix_ele : 1-D list
            List with all element numbers that are allowed to change.
        """
        # passive crack tip elements
        elx = [x for x in range(np.min(self.crack_length)-1, np.max(self.crack_length)+1)]
        ely = [self.nely-1]*len(elx)

        # passive elements at load introduction
        r = self.nely/30
        center_y = self.nelx/5
        center_x = self.nely/0.6*0.325
        x = np.vstack([np.arange(self.nelx)]*self.nely)
        y = np.hstack([np.arange(self.nely).reshape((self.nely, 1))]*self.nelx)
        distance2 = (x-center_x)**2 + (y-center_y)**2
        loc_x, loc_y = np.where(distance2 <= r**2)
        elx += loc_x.tolist()
        ely += loc_y.tolist()

        values = [1]*len(elx)

        fixele = []
        for i in range(len(elx)):
            fixele.append(self.nelx*ely[i] + elx[i])
        free_ele = list(set(range(self.nelx*self.nely)) - set(fixele))

        return elx, ely, values, free_ele
