'''
Constraints class used to specify the density constraints of the topology
optimisation problem. It contains functions for minimum and maximum element
density in the upcomming iteration and the magnitude of the volume constraint
function itself of the current design.

Bram Lagerweij
Aerospace Structures and Materials Department TU Delft
2018
'''
import numpy as np


class DensityConstraint(object):
    """
    This object relates to the constraints used in this optimization.
    It can be used for the MMA updatescheme to derive what the limit is for all
    element densities at every itteration.
    The class itself is not changed by the itterations.

    Attributes
    -------
    nelx : int
        Number of elements in x direction.
    nely : int
        Number of elements in y direction.
    move : float
        Maximum change in density of an element over 1 itteration.
    volume_frac : float
        Maximum volume that can be filled with material.
    volume_derivative : 2D array size(1, nelx*nely)
        Sensityvity of the density constraint to the density in each element.
    density_min : float (optional)
        Minumum density, set at 0.0 if not specified.
    density_max : float (optional)
        Maximum density, set at 0.0 if not specified.

    Methods
    -------
    xmin(load, x)
        Returns the minimum density value of all ellements of this itteration.
    xmax(load, x)
        Returns the maximum density value of all ellements of this itteration.
    current_volconstrain(x)
        Returns the current magnitude of the volume constraint funcion.
    """
    def __init__(self, nelx, nely, move, volume_frac, density_min=0.0, density_max=1.0):
        self.nelx = nelx
        self.nely = nely
        self.move = move
        self.volume_frac = volume_frac
        self.volume_derivative = 1/(nelx*nely*volume_frac)*np.ones((1, nely*nelx))
        self.density_min = density_min
        self.density_max = density_max

    def xmin(self, load, x):
        """
        This function calculates the minimum density value of all ellements of
        this itteration.

        Parameters
        _______
        load : object, a child object of the Loads class
            Current load case object.
        x : 2D array size(nely, nelx)
            Density distribution of this itteration.

        Returns
        _______
        xmin : 2D array size(nely, nelx)
            Minimum density values of this itteration for the update scheme.
        """
        elx, ely, value = load.passive()
        xmin = self.density_min*np.ones((self.nely, self.nelx))
        xmin = np.maximum(xmin, x - self.move)
        xmin[ely, elx] = value
        return xmin

    def xmax(self, load, x):
        """
        This function calculates the maximum density value of all ellements of
        this itteration.

        Parameters
        _______
        load : object, a child object of the Loads class
            Current load case object.
        x : 2D array size(nely, nelx)
            Density distribution of this itteration.

        Returns
        _______
        xmax : 2D array size(nely, nelx)
            Maximum density values of this itteration after updating.
        """
        elx, ely, value = load.passive()
        xmax = self.density_max*np.ones((self.nely, self.nelx))
        xmax = np.minimum(xmax, x + self.move)
        xmax[ely, elx] = value
        return xmax

    def current_volconstrain(self, x):
        """
        Calculates the current magnitude of the volume constraint funcion: ::

                           ∑ x
          cur_vol = ────────────────── - 1
                    nelx*nelx*vol_frac
        Parameters
        _______
        x : 2D array size(nely, nelx)
            Density distribution of this itteration.

        Returns
        _______
        curvol : float
            Curent value of the density constraint function.
        """
        cur_vol = np.sum(x)/(self.nelx*self.nely*self.volume_frac) - 1
        return cur_vol
