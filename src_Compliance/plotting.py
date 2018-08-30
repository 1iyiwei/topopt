"""
Plotting the simulated TopOpt geometry with boundery conditions and loads.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
"""
import matplotlib.pyplot as plt
import numpy as np


class Plot(object):
    """
    This class contains functions that allows the visualisation of the TopOpt
    algorithem. It can print the density distribution, the boundary conditions
    and the forces.

    Atributes
    --------
    x : 2-D array size(nely, nelx)
        Density distribution.
    load : object, child of the Loads class
        The loadcase(s) considerd for this optimisation problem.
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    plotting : bool
        Figure is only plotted if plotting == True

    Methods
    -------
    figure(title=None)
        Plotting the density distribution.
    boundary()
        Plotting the boundary conditions.
    loading()
        Plotting the forces acting on the problem.
    show()
        Displaying the generated figure.
    """
    def __init__(self, x, load, nelx, nely, plotting=False):
        self.x = x
        self.fixdofs = load.fixdofs()
        self.force = load.force()
        self.nelx = nelx
        self.nely = nely
        self.plotting = plotting
        self.dim = 2

    def figure(self, title=None):
        """
        Plotting the density distribution where a title can be added.

        Parameters
        ----------
        title : str
            A title that can be given to the plot.
        """
        if self.plotting is not True:
            return

        plt.figure()
        plt.imshow(1-self.x, vmin=0, vmax=1, cmap=plt.cm.gray)
        if title is not None:
            plt.title(title)
        plt.xticks([])
        plt.yticks([])

    def boundary(self):
        """
        Plotting the boundary conditions.
        """
        if self.plotting is not True:
            return
        wedgepropsH = dict(arrowstyle="wedge, tail_width=1.", color='g')
        wedgepropsV = dict(arrowstyle="wedge, tail_width=1.", color='b')

        for i in self.fixdofs:
            if i % self.dim == 0:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(-15, 0),
                             textcoords='offset points', arrowprops=wedgepropsH)

            if i % self.dim == 1:
                node = int((i)/self.dim)
                nodex = int(node/(self.nely + 1))-0.5
                nodey = node % (self.nely + 1)-0.5

                plt.annotate('', xy=(nodex, nodey), xytext=(0, -15),
                             textcoords='offset points', arrowprops=wedgepropsV)

    def loading(self):
        """
        Plotting the loading conditions.
        """
        if self.plotting is not True:
            return
        arrowprops = dict(arrowstyle="simple",fc="r", ec="r", mutation_scale=20)

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

    def show(self):
        """
        Displaying the generated figure
        """
        if self.plotting is not True:
            return
        plt.show()
