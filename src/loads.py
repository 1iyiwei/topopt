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
        return (nely+1)*elx + ely; 

    # compute the 4 boundary nodes of an element
    def nodes(self, elx, ely, nelx, nely):
        n1 = self.node(elx,     ely,     nelx, nely) 
        n2 = self.node(elx + 1, ely,     nelx, nely) 
        n3 = self.node(elx + 1, ely + 1, nelx, nely) 
        n4 = self.node(elx,     ely + 1, nelx, nely)
        return n1, n2, n3, n4

    # edof
    def edofOld(self, elx, ely, nelx, nely):
        n1, n2, n3, n4 = self.nodes(elx, ely, nelx, nely)
        return np.array([self.dim*n1,self.dim*n1+1, self.dim*n2,self.dim*n2+1, self.dim*n3,self.dim*n3+1, self.dim*n4,self.dim*n4+1])

    # edof that returns an array
    def edof(self, nelx, nely):
        """
        Generates an array with the position of the nodes of each element in
        the global stiffness matrix

        The following nodal defenitions are used:
        _________________
        >>>
                       |                                       |
        --   2*el+2*elx , 2*el+2*elx+1 ---- 2*el+2*elx+2*nely+2 , 2*el+2*elx+2*nely+3 --
                       |                                       |
                       |                                       |
                       |                                       |
                       |                                       |
                       |                                       |
                       |                                       |
                       |                                       |
        -- 2*el+2*elx+1 , 2*el+2*elx+2 ---- 2*el+2*elx+2*nely+4 , 2*el+2*elx+2*nely+5 --
                       |                                       |
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

# example loading scenario, half mbb-beam
class HalfBeam(Load):
    def __init__(self, nelx, nely):
        super().__init__(nelx, nely)
        
    def force(self):
        f = super().force()
        # downward force at the upper left corner
        f[self.dim-1] = -1.0
        return f

    def fixdofs(self):
        # left side fixed to a wall, lower right corner fixed to a point
        return ([x for x in range(0, self.dim*(self.nely+1), self.dim)] + [self.dim*(self.nelx+1)*(self.nely+1)-1])

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))
