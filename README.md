# Topology optimization #
[![](https://rawgit.com/AJJLagerweij/topopt/master/img/Cantilever_Beam.svg)](https://rawgit.com/AJJLagerweij/topopt/master/img/Canti.mp4)

## Introduction ##
This code is used in a master thesis that tries to increase damage tolerance with topology optimization. It tries to find a method that distributes density in such a way that a crack growth as slow as possible. This type of optimization is possible with the code in the [src_StressIntensity](./src_StressIntensity) folder. The other two folders, [src_Compliance](./src_Compliance) and [src_Actuator](./src_Actuator) contain versions used for the development of the code.

### src_Compliance ###
The goal if this part of the code is compliance minimization, e.g. stiffness maximization. This is an implementation of the [classic topology optimization code](http://www.topopt.dtu.dk/) described in [A 99 line topology optimization code written in Matlab](http://www.topopt.dtu.dk/files/matlab.pdf) by Ole Sigmund. The main difference is that this code is written in python and that a method of moving asymptotes update scheme is used, [MMA](https://doi.org/10.1002/nme.1620240207) is developed by Krister Svanberg. Start with [example.py](./src_Compliance/example.py).

### src_Actuator ###
This is an algorithm that maximizes the compliance on one or multiple degrees of freedom. It creates actuator like design with maximizes the displacement at a certain point, this can be used for MEMs (MicroElectroMechanical Systems) actuators. In basic, it is the same code as the src_Compliances however a different objective and linearisation is used. This is described in [Topology Optimization: Theory, Methods and Applications](www.doi.org/10.1007/978-3-662-05086-6) Sec: 2.6 'Synthesis of compliant mechanisms' by M.P. Bendsøe and O. Sigmund 2003. Start with [example.py](./src_Actuator/example.py).

### src_StressIntensity ###
Is the main goal of the project, this algorithm can minimize the stress intensity (KI) of a crack thas is assumed to be present. Minimizing (KI) for a given crack creates a design such that the given crack propegates as slow as possible within the set boundaries. Executing this optimization for multiple crack simultaniousy, with correctly set weighting factors, should maximize the time used for the crack to growth from the smallest to the longest length. The KI values are obtained with enriched cubic element as described in [Stress Intensity Factors by Enriched Finite Elements](www.doi.org/10.1016/0013-7944(78)90059-0) by N.L. Gifford and P.D. Hilton 1978. These crack tip elements are surounded by linear elements and connected with the hanging node method, resulting in a non-conform mesh. This enrichment method allows the FEA to calculate the stress intensity factors directly and without any post-prosesing steps simplifying the linearisation step which stays similar to that of the Actuator code. Start by using the [example.py](./src_StressIntensity/example.py).

## Prerequisites ##
Python 3 with NumPy, SciPy, matplotlib and cvxopt as the core. The ffmpeg packege is required only when .mp4 movies are created. To simplify the setup Anaconda enviroments (including Spyder) are avalible both for [Window](./anaconda/TopOpt_Windows.yml) and [Linux](./anaconda/TopOpt_Linux.yml).

## Generating a new load class ##
The folowing section will explain how to set up a simulation for a new geometry. One is expected to have some knowlege of FEM and topology optimization. No explanation on the simplation settings, such as the resolution, filter size or volume constrain will follow. For information on those topics I recomend reading "Topology Optimization" from M.P. Bends&#248;e and O. Sigmund.

Making a new load class is fairly simple and similar for all optimization types, make the following steps:
 1. Open the [loads.py](./src_Compliance/loads.py) file
 2. Copy the [HalfBeam class](https://github.com/AJJLagerweij/topopt/blob/49d6241989899567b59f3e474c6ad6799b4d0641/src_Compliance/loads.py#L241-L284) to the botom of the file and change the name of the class.
 3. Change the boundery conditions e.g. the load vector and fix certain degrees of freedom (fixdofs). The folowing section wil show how to do it for an example problem.
 4. Change the passive elements defenition, those that do not change in desity, an example cousd be setting the density of certain elements to 0 as a hole is planned at that location.

To allow the load vector and fixdofs to change with your mesh size it is important to formulate it as equations of the amount of elements in x and y direction. The nodes are numbered from the top left corner and go down to the botom before starting on the next column of nodes. Function [nodes](https://github.com/AJJLagerweij/topopt/blob/49d6241989899567b59f3e474c6ad6799b4d0641/src_Compliance/loads.py#L89-L116) can be called upon to finding the nodal coodinates of point as a function of the number of elements in x and y directly. Take the HalfBeam example which among others is fixed in the y direction at the botom end of the beam. Then using `n1, n2, n3, n4 = self.nodes(self.nelx-1, self.nely-1)` returns the nodal coordinates of the element at position nelx-1, nely-1. (the -1 is used as python starts counting from 0) The position of `n1` is in the top left the others are defined in clockwise order. Thus, the botom right node, which is the one we need, is `n3`. To convert the nodal position to the position of the x orientation of the node in the displacement or load vector simply multiply by the dimensions used (`self.dim`). For the location of the y orientation simply add `+1` to the x location.

Sometimes a wole side is fixed (due to a wall or symetry) than using the `range(start, stop, step)` funcion can be usefull. Take for exapmle the HalfBeam again. On the left side all movement in x direction is constrained due to symetry. It is known that ranging numbers from 0 (the top left node its x orientation) to the botom of the last element (`self.dim*(self.nely+1)`) with steps of `self.dim` (2 in the 2D case) contains all the locations in the displacement vector that need to be fixed. Therefore, the code `[x for x in range(0, self.dim*(self.nely+1), self.dim)]` retuns a list with all these locations.

An [example](https://github.com/AJJLagerweij/topopt/blob/2641a28770d6f0f6ed1ae9930a8e156444244b01/src_StressIntensity/loads.py#L910-L936) of passive elements is used in the StressIntensity minimization code. Defining these passive elements done by three list (or arrays) need to be exported, elx contains the x coordinates of the fixed element, ely the y coordinates and lastly the values list contains the value at the given coordinate. The order of the elements should be the same in all three lists. Then the element numbers of the free elements can be found when running:

```python
fixele = []
for i in range(len(elx)):
    fixele.append(self.nelx*ely[i] + elx[i])
free_ele = list(set(range(self.nelx*self.nely)) - set(fixele))
```
![HalfBeamFBD](https://rawgit.com/AJJLagerweij/topopt/master/img/HalfBeamFBD.svg)

The same steps need to be taken to create a new actuator load case, but copy the [Inverter class](https://github.com/AJJLagerweij/topopt/blob/49d6241989899567b59f3e474c6ad6799b4d0641/src_Actuator/loads.py#L263-L321) located in the src_Actuator folder. There is two extra steps however:
 1. Add the location and orientation of the actuator output. This is done with the [displaceloc function](https://github.com/AJJLagerweij/topopt/blob/49d6241989899567b59f3e474c6ad6799b4d0641/src_Actuator/loads.py#L297) which works similar as a force vector.
 2. Chose the [ext_stiffness](https://github.com/AJJLagerweij/topopt/blob/49d6241989899567b59f3e474c6ad6799b4d0641/src_Actuator/example.py#L25) which is the external spring stiffness of the in and output used to stabelize the problem.
 
 ![InverterFBD](https://rawgit.com/AJJLagerweij/topopt/master/img/Inverter.svg)
 
Using the StressIntensity load class works similar to the Compliance minimization one. Copy the [Compact Tension](https://github.com/AJJLagerweij/topopt/blob/2641a28770d6f0f6ed1ae9930a8e156444244b01/src_StressIntensity/loads.py#L784-L936) case and fix the loads and fixed degrees of freedom. Extra attention should be provided to the location of the Higher Order Elements (defined as `hoe`) and the passive elements. All elements that are used in the crack tip analysis must have a density value of 1 exactly as shown in the CompactTension example.

![CompactTension](https://rawgit.com/AJJLagerweij/topopt/master/img/CompactTension.svg)

## To be implemented ##
Four improvements are proposed:
 1. Adding a faster algebraic multigrid preconditioner with Conjugate Gradient itterative solver.
 2. Creating a conform cubic element mesh for the StressIntensity minimization with a multiresolution implementation as shown in [Higher‐order multi‐resolution topology optimization using the finite cell method](https://doi.org/10.1002/nme.5432) by J.P. Groen, M. Langelaar, O. Sigmund, and M. Ruess 2017. This can both improof the performance of the optimization and acuracy of the stress intensity factors calculated.
 3. Add multiloadcase optimisation capebilities to Compliance and Actuator code.

## Special Thanks To ##
 1. Ole Sigmund and Martin Bends&#248;e from the Technical University of Denmark for their contributions to the field of topology optimization
 2. Li-Yi Wei from the University of Hong Kong for making a beautiful [python code](https://github.com/1iyiwei/topopt)
 3. Svanberg, K. (1987). The method of moving asymptotes — a new method for structural optimization. International Journal for Numerical Methods in Engineering, 24(2), 359–373. DOI: [10.1002/nme.162024027](https://doi.org/10.1002/nme.1620240207)
