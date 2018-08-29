'''
loading scenario aka boundary condition
density at elements, force/displacement at nodes (element corners)
dofs are for displacements
'''

import numpy as np


class Load(object):
    def __init__(self, nelx, nely, young, Emin, poisson):
        self.nelx = nelx
        self.nely = nely
        self.young = young
        self.Emin = Emin
        self.poisson = poisson
        self.dim = 2

#    '''
#    different convention from numpy array shape: x/y versus row/col
#    '''
#    def shape(self):
#        return (self.nely, self.nelx)

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

    # element (local) stiffness matrix
    def lk(self, E, nu):
        """ Returns a siffness matrix depending on E and nu """
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
        return np.zeros(self.dim*(self.nely+1)*(self.nelx+1))

    def alldofs(self):
        return [x for x in range(self.dim*(self.nely+1)*(self.nelx+1))]

    def fixdofs(self):
        return []

    def freedofs(self):
        return self.alldofs()


# example loading scenario, half mbb-beam (symetry around y axis)
class HalfBeam(Load):
    def __init__(self, nelx, nely, young, Emin, poisson):
        super().__init__(nelx, nely, young, Emin, poisson)

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


# example loading scenario, half mbb-beam without use of symetry axis
class Beam(Load):
    def __init__(self, nelx, nely, young, Emin, poisson):
        super().__init__(nelx, nely, young, Emin, poisson)
        if nelx % 2 != 0:
            raise ValueError('nelx needs to be even in a mbb beam')

    def force(self):
        f = super().force()
        # downward force in the middle of the upper part
        n1, n2, n3, n4 = self.nodes(int(self.nelx/2), 0, self.nelx, self.nely)
        f[self.dim*n1+1] = -1.0
        return f

    def fixdofs(self):
        # left side fixed to a wall (x direction only), lower right corner fixed to a point
        n11, n12, n13, n14 = self.nodes(0, self.nely-1, self.nelx, self.nely)
        n21, n22, n23, n24 = self.nodes(self.nelx-1, self.nely-1, self.nelx, self.nely)
        return ([self.dim*n14, self.dim*n14+1, self.dim*n23+1])

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = []  # np.arange(0, self.nelx/2, dtype=int)
        ely = []  # np.ones(np.shape(xlist), dtype=int)*int(self.nely/2)
        values = []  # np.ones(np.shape(xlist))*0.001
        return elx, ely, values


# cantilever beam load in the middle of the end
class Canti(Load):
    def __init__(self, nelx, nely, young, Emin, poisson):
        super().__init__(nelx, nely, young, Emin, poisson)
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
    def __init__(self, nelx, nely, young, Emin, poisson):
        super().__init__(nelx, nely, young, Emin, poisson)
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
    def __init__(self, nelx, nely, young, Emin, poisson):
        super().__init__(nelx, nely, young, Emin, poisson)

    def force(self):
        f = super().force()
        # uprard force at the top side and half load in corner
        loc_up = np.arange(1, self.dim*(self.nelx+1)*(self.nely+1), self.dim*(self.nely+1))
        f[loc_up] = -1
        # right force at the right hand side
        loc_right = np.arange(self.dim*(self.nely+1)*self.nelx, self.dim*(self.nely+1)*(self.nelx+1), self.dim)
        f[loc_right] = -1
        # bottom force down
        loc_down = np.arange(self.dim*(self.nely+1)-1, self.dim*(self.nely+1)*(self.nelx+1), self.dim*(self.nely+1))
        f[loc_down] = 1
        # left force left
        loc_left = np.arange(0, self.dim*(self.nely+1), self.dim)
        f[loc_left] = 1
        return f

    def fixdofs(self):
        return  [0, self.dim*(self.nely), self.dim*(self.nely)+1, self.dim*(self.nelx+1)*(self.nely+1)-1]

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = 2*[x for x in range(self.nelx)] + [0 for y in range(self.nely)] \
                + [self.nelx-1 for y in range(self.nely)] + [0+1, self.nelx-2, 0+1, self.nelx-2]
        ely = [0 for x in range(self.nelx)] + [self.nely-1 for x in range(self.nelx)] \
                + 2*[y for y in range(self.nely)] + [0+1, 0+1, self.nely-2, self.nely-2]
        values = np.ones((len(elx)))
        return elx, ely, values
