'''
loading scenario aka boundary condition
density at elements, force/displacement at nodes (element corners)
dofs are for displacements
'''

import numpy as np


class Load(object):
    def __init__(self, nelx, nely, ext_stiff):
        self.nelx = nelx
        self.nely = nely
        self.dim = 2
        self.ext_stiff = ext_stiff

    '''
    different convention from numpy array shape: x/y versus row/col
    '''
    def shape(self):
        return (self.nelx, self.nely)

    # compute 1D index from 2D position for node (boundary of element)
    def node(self, elx, ely, nelx, nely):
        return (nely+1)*elx + ely

    # compute the 4 boundary nodes of an element
    # requires more testing, errors might exist here
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
        return np.zeros((self.dim*(self.nely+1)*(self.nelx+1), 1))

    def displaceloc(self):
        return np.zeros((self.dim*(self.nely+1)*(self.nelx+1), 1))

    def alldofs(self):
        return [x for x in range(self.dim*(self.nely+1)*(self.nelx+1))]

    def fixdofs(self):
        return []

    def freedofs(self):
        return self.alldofs()


# example loading scenario, inverter with horizontal mirror axis
class Inverter(Load):
    def __init__(self, nelx, nely, ext_stiff):
        super().__init__(nelx, nely, ext_stiff)

    def force(self):
        f = super().force()
        # positive horizontal force bottom left corner
        f[self.dim*self.nely] = 1.0
        return f

    def displaceloc(self):
        dout = super().displaceloc()
        # actuater location is horizontal in the bottom right
        # positive movement is in the negative x direction
        dout[-2] = -1
        return dout

    def fixdofs(self):
        # bottom fixed to the wall (y direction only), top left corner fixed to a point
        fix1 = np.arange(self.dim*self.nely + 1, self.dim*(self.nelx+1)*(self.nely+1), self.dim*(self.nely+1))
        return (np.hstack((fix1, [0, 1, 2, 3]))).tolist()

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = []  # np.arange(0, self.nelx/2, dtype=int)
        ely = []  # np.ones(np.shape(xlist), dtype=int)*int(self.nely/2)
        values = []  # np.ones(np.shape(xlist))*0.001
        return elx, ely, values
