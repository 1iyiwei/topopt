'''
loading scenario aka boundary condition
density at elements, force/displacement at nodes (element corners)
dofs are for displacements
'''

import numpy as np


class Load(object):
    def __init__(self, nelx, nely):
        self.nelx = nelx
        self.nely = nely
        self.dim = 2

    '''
    different convention from numpy array shape: x/y versus row/col
    '''
    def shape(self):
        return (self.nelx, self.nely)

    # compute 1D index from 2D position for node (boundary of element)
    def node(self, elx, ely, nelx, nely):
        return (nely+1)*elx + ely

    # compute the 4 boundary nodes of an element
    def nodes(self, elx, ely, nelx, nely):
        n1 = self.node(elx,     ely,     nelx, nely) 
        n2 = self.node(elx + 1, ely,     nelx, nely) 
        n3 = self.node(elx + 1, ely + 1, nelx, nely) 
        n4 = self.node(elx,     ely + 1, nelx, nely)
        return n1, n2, n3, n4

    # edof that returns an array
    def edof(self, nelx, nely):
        """
        Generates an array with the position of the nodes of each element in
        the global stiffness matrix
        """
        # Creating list with element numbers
        elx = np.repeat(range(nelx), nely).reshape((nelx*nely, 1))  # x position of element
        ely = np.tile(range(nely), nelx).reshape((nelx*nely, 1))  # y position of element

        n1, n2, n3, n4 = self.nodes(elx, ely, nelx, nely)
        edof = np.array([self.dim*n1,self.dim*n1+1, self.dim*n2,self.dim*n2+1,
                         self.dim*n3,self.dim*n3+1, self.dim*n4,self.dim*n4+1])
        edof = edof.T[0]

        x_list = np.repeat(edof, 8)  # flat list pointer of each node in an element
        y_list = np.tile(edof, 8).flatten()  # flat list pointer of each node in element
        return edof, x_list, y_list

    def force(self):
        return np.zeros(self.dim*(self.nely+1)*(self.nelx+1))

    def alldofs(self):
        return [x for x in range(self.dim*(self.nely+1)*(self.nelx+1))]

    def fixdofs(self):
        return []

    def freedofs(self):
        return self.alldofs()


# example loading scenario, half mbb-beam (symetry around y axis)
class HalfBeam(Load):
    def __init__(self, nelx, nely):
        super().__init__(nelx, nely)

    def force(self):
        f = super().force()
        # downward force at the upper left corner
        f[self.dim-1] = -1.0
        return f

    def fixdofs(self):
        # left side fixed to a wall (x direction only), lower right corner fixed to a point
        n1, n2, n3, n4 = self.nodes(self.nelx-1, self.nely-1, self.nelx, self.nely)
        return ([x for x in range(0, self.dim*(self.nely+1), self.dim)] + [self.dim*n3+1])

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = []  # np.arange(0, self.nelx/2, dtype=int)
        ely = []  # np.ones(np.shape(xlist), dtype=int)*int(self.nely/2)
        values = []  # np.ones(np.shape(xlist))*0.001
        return elx, ely, values


# cantilever beam load in the middle of the end
class Canti(Load):
    def __init__(self, nelx, nely):
        super().__init__(nelx, nely)
        if nely % 2 != 0:
            raise ValueError('nely needs to be even in a cantilever beam')

    def force(self):
        f = super().force()
        # downward force at the right side of the cantilever
        n1, n2, n3, n4 = self.nodes(self.nelx-1, int(self.nely/2), self.nelx, self.nely)
        f[self.dim*n2+1] = -1
        return f

    def fixdofs(self):
        # left side fixed to a wall (x, and y directions)
        return ([x for x in range(0, self.dim*(self.nely+1))])

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = []
        ely = []
        values = []
        return elx, ely, values


# the Michell structures wich analytical salutions (symetry arround y axsi)
class Michell(Load):
    def __init__(self, nelx, nely):
        super().__init__(nelx, nely)
        if nely % 2 != 0:
            raise ValueError('nely needs to be even in a michell strucure')

    def force(self):
        f = super().force()
        # downward force at the right side of the cantilever
        n1, n2, n3, n4 = self.nodes(0, int(self.nely/2), self.nelx, self.nely)
        f[self.dim*n2+1] = -1
        return f

    def fixdofs(self):
        # at the symety axis x direction is fixed, and point fix at middle end
        n1, n2, n3, n4 = self.nodes(self.nelx-1, int(self.nely/2), self.nelx, self.nely)
        return ([self.dim*n2+1]+[x for x in range(0, self.dim*(self.nely+1), self.dim)])

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = []
        ely = []
        values = []
        return elx, ely, values

# the biaxcal tension example
class BiAxial(Load):
    def __init__(self, nelx, nely):
        super().__init__(nelx, nely)

    def force(self):
        f = super().force()
        # uprard force at the top side
        loc_up = np.arange(self.dim*(self.nely+1)+1, self.dim*(self.nelx)*(self.nely+1), self.dim*(self.nely+1))
        f[loc_up] = 1
        # right force at the right hand side
        loc_right = np.arange(self.dim*(self.nely+1)*self.nelx+2, self.dim*(self.nely+1)*(self.nelx+1)-2, self.dim)
        f[loc_right] = 1
        # bottom force down
        loc_down = np.arange(2*self.dim*(self.nely+1)-1, self.dim*(self.nely+1)*(self.nelx), self.dim*(self.nely+1))
        f[loc_down] = -1
        # left force left
        loc_left = np.arange(2, self.dim*(self.nely+1)-2, self.dim)
        f[loc_left] = -1
        return f

    def fixdofs(self):
        return  [0, self.dim*(self.nely), self.dim*(self.nely)+1, self.dim*(self.nelx+1)*(self.nely+1)-1]

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        return [], [], []

