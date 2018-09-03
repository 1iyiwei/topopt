# Topology optimization #

[![](https://rawgit.com/AJJLagerweij/topopt/master/img/Cantilever_Beam.svg)](https://rawgit.com/AJJLagerweij/topopt/master/img/Canti.mp4)

This is an implementation of the [classic topology optimization code](http://www.topopt.dtu.dk/) described in [A 99 line topology optimization code written in Matlab](http://www.topopt.dtu.dk/files/matlab.pdf) by Ole Sigmund. The main difference is that this code is written in python and that a method of moving asymptotes update scheme is used, [MMA](https://doi.org/10.1002/nme.1620240207) is developed by Krister Svanberg. Currently three versions of the code exist, a compilance minimalization, a actuator design version that maximizes the displacement at a certain point and a stress intensity minimisation code.
Start with [example.py](./src_Compliance/example.py).

## Prerequisites ##
Python 3 with NumPy, SciPy, matplotlib and cvxopt as the core. The ffmpeg packege is required only when .mp4 movies are created. To simplify the setup Anaconda enviroments (including Spyder) are avalible both for [Window](./anaconda/TopOpt_Windows.yml) and [Linux](./anaconda/TopOpt_Linux.yml).

## Generating a new load class ##
The folowing section will explain how to set up a simulation for a new geometry. One is expected to have some knowlege of FEM and topology optimization. No explanation on the simplation settings, such as the resolution, filter size or volume constrain will follow. For information on those topics I recomend reading "Topology Optimization" from M.P. Bends&#248;e and O. Sigmund.

Making a new load class is fairly simple, make the following steps:
 1. Open the [loads.py](./src_Compliance/loads.py) file
 2. Copy the [HalfBeam class](https://github.com/AJJLagerweij/topopt/blob/0343ef3eeac7ce95a76d9c00cfbf2ee66c383696/src_Compliance/loads.py#L66-L89) to the botom of the file and change the name of the class.
 3. Change the boundery conditions e.g. the load vector and fix certain degrees of freedom (fixdofs). The folowing section wil show how to do it for an example problem.
 4. Change the passive elements defenition, those that do not change in desity, for example with a density of 0.0001 as a hole is planned there.

To allow the load vector and fixdofs to change with your mesh size it is important to formulate it as equations of the amount of elements in x and y direction. The nodes are numbered from the top left corner and go down to the botom before starting on the next column of nodes. Funciton [nodes](https://github.com/AJJLagerweij/topopt/blob/0343ef3eeac7ce95a76d9c00cfbf2ee66c383696/src_Compliance/loads.py#L27) can helps with finding the nodal coodinates of point as a function of the number of elements in x and y direction. Take the HalfBeam example which among others is fixed in the y direction at the botom end of the beam. Then using `n1, n2, n3, n4 = self.nodes(self.nelx-1, self.nely-1, self.nelx, self.nely)` returns the nodal coordinates of the element at position nelx-1, nely-1. (the -1 is used as python starts counting from 0) The position of `n1` is in the top left the others are defined in clockwise order. Thus the botom right node, which is the one we need, is `n3`. To convert the nodal position to the position of the x orientation of the node in the displacement or load vector simply multiply by the dimensions used (`self.dim`). For the location of the y orientation simply add `+1` to the x location.

Sometimes a wole side is fixed (due to a wall or symetry) than using the `range(start, stop, step)` can be usefull. Take for exapmle the HalfBeam again. At the left side all movement in x direction is constrained due to symetry. As it is known that ranging numbers from 0 (the top left node its x orientation) to the botom of the last element (`self.dim*(self.nely+1)`) with steps of `self.dim` (2 in the 2D case) contains all the locations in the displacement vector that need to be fixed. Thus the code `[x for x in range(0, self.dim*(self.nely+1), self.dim)]` retuns a list with all these locations.

The [passive](https://github.com/AJJLagerweij/topopt/blob/0343ef3eeac7ce95a76d9c00cfbf2ee66c383696/src_Compliance/loads.py#L85-L89) elements defenition is fairly simple, three list (or arrays) need to be exported, elx contains the x coordinates of the fixed element, ely the y coordinates and lastly the values list contains the value at the given coordinate. The order of the elements should be the same in all three lists.

The same steps need no be taken to create a new actuator load case, but copy the [Inverter class](https://github.com/AJJLagerweij/topopt/blob/455e2723b66b0be2c2eea74a514b9ee37709f416/src_Actuator/loads.py#L72-L101) located in the src_Actuator folder. There is two extra steps however:
 1. Add the location and orientation of the actuator output. This is done with the [displaceloc function](https://github.com/AJJLagerweij/topopt/blob/455e2723b66b0be2c2eea74a514b9ee37709f416/src_Actuator/loads.py#L82-L87) which works similar as a force vector.
 2. Chose the [ext_stiffness](https://github.com/AJJLagerweij/topopt/blob/59266aa4b883c275de1d3175ea43ad3af0239c06/src_Actuator/example.py#L18) which is the external spring stiffness of the in and output.

![HalfBeamFBD](https://rawgit.com/AJJLagerweij/topopt/master/img/HalfBeamFBD.svg)

## To be implemented ##
Currently four changes are propose:
 1. Adding a faster algebraic multigrid preconditioner with Conjugate Gradient itterative solver.
 2. Fixing the division by zero error in the MMA combined with passive elements.
 3. Improof documentation of the Inverter and StressIntensity codes.
 4. Add multiloadcase optimisation capebilities.

## Special Thanks To ##
 1. Ole Sigmund and Martin Bends&#248;e from the Technical University of Denmark for there contributions to the field of topology optimization
 2. Li-Yi Wei from the University of Hong Kong for making a beautiful [python code](https://github.com/1iyiwei/topopt)
 3. Svanberg, K. (1987). The method of moving asymptotes — a new method for structural optimization. International Journal for Numerical Methods in Engineering, 24(2), 359–373. DOI: [10.1002/nme.162024027](https://doi.org/10.1002/nme.1620240207)
