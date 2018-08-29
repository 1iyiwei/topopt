"""
Plotting the simulated TopOpt geometry with boundery conditions and loads
"""
import matplotlib.pyplot as plt
import numpy as np

class Plot(object):
    def __init__(self, x, load, nelx, nely):
        self.x = x
        self.edof = load.edof
        self.fixdofs = load.fixdofs()
        self.force = load.force()
        self.nelx = load.nelx
        self.nely = load.nely
        self.dim = 2

    def figure(self, constraint, title=None):
        plt.figure()
        vmax = 1
        vmin = -1
        plt.imshow(1-self.x, vmin=vmin, vmax=vmax, cmap=plt.cm.gray)
        if title != None:
            plt.title(title)
        plt.xticks([])
        plt.yticks([])

    def find(self, dof):
        """ This function returns the location, x,y of any degree of freedom
        by corresponding it with the edof array """
        # find first location in edofs
        # check length element
        # place at correct location in elment
        elenum = 0
        for sublist in self.edof:
            subnum = 0
            for item in sublist:
                if item == dof:
                    break
                subnum += 1
            if item == dof:
                break
            elenum += 1

        # element center coordinates
        elx = int(elenum/self.nely)
        ely = elenum % self.nely

        # type of element, 8 or 26 dof
        ele_len = len(self.edof[elenum])

        # when element is 8 dofs
        if ele_len == 8:
            xlist = [-0.5, 0.5, 0.5, -0.5]
            ylist = [0.5, 0.5, -0.5, -0.5]
            x = elx + xlist[int(subnum/2)]
            y = ely + ylist[int(subnum/2)]

        # when element has 26 dofs
        if ele_len == 26:
            xlist = [-0.5, -0.5+1/3, -0.5+2/3, 0.5, 0.5, 0.5, 0.5, 0.5-1/3, 0.5-2/3, -0.5, -0.5, -0.5, 0]
            ylist = [0.5, 0.5, 0.5, 0.5, 0.5-1/3, 0.5-2/3, -0.5, -0.5, -0.5, -0.5, -0.5+1/3, -0.5+2/3, 0]
            x = elx + xlist[int(subnum/2)]
            y = ely + ylist[int(subnum/2)]

        return x, y

    def boundary(self):
        wedgepropsH = dict(arrowstyle="wedge, tail_width=1.", color='g')
        wedgepropsV = dict(arrowstyle="wedge, tail_width=1.", color='b')

        for i in self.fixdofs:
            nodex, nodey = self.find(i)
            if i % self.dim == 0:
                plt.annotate('', xy=(nodex, nodey), xytext=(-15, 0),
                             textcoords='offset points', arrowprops=wedgepropsH)
            if i % self.dim == 1:
                plt.annotate('', xy=(nodex, nodey), xytext=(0, -15),
                             textcoords='offset points', arrowprops=wedgepropsV)

    def loading(self):
        arrowprops = dict(arrowstyle="<-",fc="r", ec="r", mutation_scale=20)
        forceloc = np.nonzero(self.force)[0]

        for i in forceloc:
            force = self.force[i]
            nodex, nodey = self.find(i)
            if i % self.dim == 0:
                plt.annotate('', xy=(nodex, nodey), xytext=(30*force, 0),
                             textcoords='offset points', arrowprops=arrowprops)
            if i % self.dim == 1:
                plt.annotate('', xy=(nodex, nodey), xytext=(0, 30*force),
                             textcoords='offset points', arrowprops=arrowprops)

    def show(self, block=True):
        plt.show()
