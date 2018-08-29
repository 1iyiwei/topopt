'''
loading scenario aka boundary condition
density at elements, force/displacement at nodes (element corners)
dofs are for displacements
'''

import numpy as np


class Load(object):
    def __init__(self, nelx, nely, ext_stiff, hoe):
        self.nelx = nelx
        self.nely = nely
        self.dim = 2
        self.ext_stiff = ext_stiff
        self.edof, self.x_list, self.y_list, self.maxnode = self.edof(nelx, nely, hoe)

    '''
    different convention from numpy array shape: x/y versus row/col
    '''
    def shape(self):
        return (self.nelx, self.nely)

    # compute 1D index from 2D position for node (boundary of element)
    def node(self, elx, ely):
        return (self.nely+1)*elx + ely

    # compute the 4 boundary nodes of an element
    # requires more testing, errors might exist here
    def nodes(self, elx, ely):
        n3 = self.node(elx,     ely    ) 
        n2 = self.node(elx + 1, ely    ) 
        n1 = self.node(elx + 1, ely + 1) 
        n0 = self.node(elx,     ely + 1)
        return n0, n1, n2, n3

    # edof that returns an array
    def edof(self, nelx, nely, hoe):
        """
        Generates an 2D list with the position of the nodes of each element in
        the global stiffness matrix
        """
        # Creating list with element numbers
        elx = np.repeat(range(nelx), nely).reshape((nelx*nely, 1))  # x position of element
        ely = np.tile(range(nely), nelx).reshape((nelx*nely, 1))  # y position of element

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
            element = nely*hoei[0] + hoei[1]
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
        maxnode = max(x_list)

        return edof, np.array(x_list), np.array(y_list), maxnode

    def force(self):
        max_node = self.maxnode + 1
        return np.zeros((max_node, 1))

    def displaceloc(self):
        max_node = self.maxnode + 1
        return np.zeros((max_node, 1))

    def alldofs(self):
        max_node = self.maxnode + 1
        return [x for x in range(max_node)]

    def fixdofs(self):
        return []

    def freedofs(self):
        return self.alldofs()


# example loading scenario, inverter with horizontal mirror axis
class EdgeCrack(Load):
    def __init__(self, nelx, nely, crack_length, ext_stiff):
        self.crack_length = crack_length
        self.hoe = [[crack_length-1, nely-1], [crack_length, nely-1]]
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, nely, ext_stiff, self.hoe)

    def force(self):
        f = super().force()
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        forceloc = [self.edof[i][5] for i in top_ele][:-1]
        f[forceloc] = 2
        f[1] = 1
        f[[self.edof[i][5] for i in top_ele][-1]] = 1
        return f

    def displaceloc(self):
        dout = super().displaceloc()
        # higher order element
        ele = self.hoe[0][0]*self.nely + self.hoe[0][1]
        # node with K_I
        node = self.edof[ele][-2]
        dout[node] = -1
        return dout

    def fixdofs(self):
        # higher order element after crack
        ele = self.nely*(self.crack_length+1)-1
        bottom = [1, 3, 5]
        fix = [self.edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(self.crack_length+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [self.edof[i][1] for i in bot_ele] + [self.edof[bot_ele[-1]][3]]

        # top elements
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        # all elements in x direction fixed
        fix2 = [0] + [self.edof[i][4] for i in top_ele]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = [self.crack_length-1, self.crack_length]
        ely = [self.nely-1, self.nely-1]
        values = [1 for x in elx]
        return elx, ely, values


class DoubleEdgeCrack(Load):
    def __init__(self, nelx, ext_stiff):
        self.nelx = nelx
        self.nely = 4*nelx
        self.crack_length = int(nelx/5*2)
        self.hoe = [[self.crack_length-1, self.nely-1], [self.crack_length, self.nely-1]]
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, self.nely, ext_stiff, self.hoe)

    def force(self):
        f = super().force()
        top_ele = [x for x in range(0, self.nelx*self.nely, self.nely)]
        forceloc = [self.edof[i][5] for i in top_ele]
        f[forceloc] = 2
        f[1] = 1
        return f

    def displaceloc(self):
        dout = super().displaceloc()
        # higher order element
        ele = self.hoe[0][0]*self.nely + self.hoe[0][1]
        # node with K_I
        node = self.edof[ele][-2]
        dout[node] = -1
        return dout

    def fixdofs(self):
        # higher order element after crack
        ele = self.nely*(self.crack_length+1)-1
        bottom = [1, 3, 5]
        fix = [self.edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(self.crack_length+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [self.edof[i][1] for i in bot_ele] + [self.edof[bot_ele[-1]][3]]

        # right side elements
        start_ele = (self.nelx-1)*self.nely
        right_ele = [x for x in range(start_ele, self.nelx*self.nely, 1)]
        # all elements in x direction fixed
        fix2 = [self.edof[i][4] for i in right_ele] + [self.edof[bot_ele[-1]][2]]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = [self.crack_length-1, self.crack_length]
        ely = [self.nely-1, self.nely-1]
        values = [1 for x in elx]
        return elx, ely, values


class CompactTension(Load):
    def __init__(self, nelx, crack_length, ext_stiff):
        self.nelx = nelx
        self.nely = int(np.round(nelx/1.25*1.2/2))
        self.crack_length = crack_length
        self.hoe = [[self.crack_length-1, self.nely-1], [self.crack_length, self.nely-1]]
        self.hoe_type = ['1,-1', '-1,-1']

        super().__init__(nelx, self.nely, ext_stiff, self.hoe)

    def force(self):
        f = super().force()
        load_ele = int(self.nely/0.6*0.325)+int(self.nelx/5)*self.nely
        forceloc = self.edof[load_ele][3]
        f[forceloc] = 1
        return f

    def displaceloc(self):
        dout = super().displaceloc()
        # higher order element
        ele = self.hoe[0][0]*self.nely + self.hoe[0][1]
        # node with K_I
        node = self.edof[ele][-2]
        dout[node] = -1
        return dout

    def fixdofs(self):
        # higher order element after crack
        ele = self.nely*(self.crack_length+1)-1
        bottom = [1, 3, 5]
        fix = [self.edof[ele][i] for i in bottom]

        # bottom elements fixed in y direction after the crack elements
        start_ele = self.nely*(self.crack_length+2)-1
        bot_ele = [x for x in range(start_ele, self.nelx*self.nely, self.nely)]
        # all bottem left y direction dof of all these elements
        fix1 = [self.edof[i][1] for i in bot_ele] + [self.edof[bot_ele[-1]][3]]

        # load introduction in y direction rigid
        load_ele = int(self.nely/0.6*0.325)+int(self.nelx/5)*self.nely
        fix2 = self.edof[load_ele][2]

        return (np.hstack((fix, fix1, fix2))).tolist()

    def freedofs(self):
        return list(set(self.alldofs()) - set(self.fixdofs()))

    def passive(self):
        elx = [self.crack_length-1, self.crack_length]
        ely = [self.nely-1, self.nely-1]
        values = [1 for x in elx]
        return elx, ely, values