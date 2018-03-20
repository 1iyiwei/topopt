"""
Plotting the simulated TopOpt geometry with boundery conditions and loads
"""
import matplotlib.pyplot as plt
import numpy as np

class Plot(object):
    def __init__(self, x, load, nelx, nely):
        self.x = x
        self.fixdofs = load.fixdofs()
        self.force = load.force()
        self.nelx = nelx
        self.nely = nely
        self.dim = 2

    def figure(self, title=None):
        plt.figure()
        plt.imshow(1-self.x, cmap=plt.cm.gray)
        plt.title(title)
        plt.xticks([])
        plt.yticks([])

    def boundary(self):
        wedgeprops = dict(arrowstyle="wedge, tail_width=1.", color='r')

        for i in self.fixdofs:
            if i % self.dim == 0:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(-15, 0),
                             textcoords='offset points', arrowprops=wedgeprops)

            if i % self.dim == 1:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(0, -15),
                             textcoords='offset points', arrowprops=wedgeprops)

    def loading(self):
        arrowprops = dict(arrowstyle="simple",fc="g", ec="g", mutation_scale=20)

        forceloc = np.nonzero(self.force)[0]

        for i in forceloc:
            force = self.force[i]

            if i % self.dim == 0:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(-60*force, 0),
                             textcoords='offset points', arrowprops=arrowprops)

            if i % self.dim == 1:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(0, -60*force),
                             textcoords='offset points', arrowprops=arrowprops)
