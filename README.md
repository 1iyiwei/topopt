[![Documentation Status](https://readthedocs.org/projects/topopt/badge/?version=latest)](https://topopt.readthedocs.io/en/latest/?badge=latest)
# Topology optimization #
[![](https://rawgit.com/AJJLagerweij/topopt/master/img/Cantilever_Beam.svg)](https://rawgit.com/AJJLagerweij/topopt/master/img/Canti.mp4)

## Introduction ##
This code is used in a master thesis that tries to increase damage tolerance with topology optimization. It tries to find a method that distributes density in such a way that a crack growth as slow as possible. This type of optimization is possible with the code in the [src_StressIntensity](./src_StressIntensity) folder. The other two folders, [src_Compliance](./src_Compliance) and [src_Actuator](./src_Actuator) contain versions used for the development of the code.

### src_Compliance ###
The goal if this part of the code is compliance minimization, e.g. stiffness maximization. This is an implementation of the [classic topology optimization code](http://www.topopt.dtu.dk/) described in [A 99 line topology optimization code written in Matlab](http://www.topopt.dtu.dk/files/matlab.pdf) by Ole Sigmund. The main difference is that this code is written in python and that a method of moving asymptotes update scheme is used, [MMA](https://doi.org/10.1002/nme.1620240207) is developed by Krister Svanberg. Start with [example.py](./src_Compliance/example.py).

### src_Actuator ###
This is an algorithm that maximizes the compliance on one or multiple degrees of freedom. It creates actuator like design with maximizes the displacement at a certain point, this can be used for MEMs (MicroElectroMechanical Systems) actuators. In basic, it is the same code as the src_Compliances however a different objective and linearisation is used. This is described in [Topology Optimization: Theory, Methods and Applications](www.doi.org/10.1007/978-3-662-05086-6) Sec: 2.6 'Synthesis of compliant mechanisms' by M.P. Bendsøe and O. Sigmund 2003. Start with [example.py](./src_Actuator/example.py).

### src_StressIntensity ###
This algorithm can minimize the stress intensity (KI) of a crack thas is assumed to be present. Minimizing (KI) for a given crack creates a design such that the given crack propegates as slow as possible within the set boundaries. Executing this optimization for multiple crack simultaniousy, with correctly set weighting factors, should maximize the time used for the crack to growth from the smallest to the longest length. The KI values are obtained with enriched cubic element as described in [Stress Intensity Factors by Enriched Finite Elements](www.doi.org/10.1016/0013-7944(78)90059-0) by N.L. Gifford and P.D. Hilton 1978. These crack tip elements are surounded by linear elements and connected with the hanging node method, resulting in a non-conform mesh. This enrichment method allows the FEA to calculate the stress intensity factors directly and without any post-prosesing steps simplifying the linearisation step which stays similar to that of the Actuator code. Start by using the [example.py](./src_StressIntensity/example.py).

### src_FatigueLive ###
This optimization maximizes the fatigue crack life, which is the amount of cycles that is required to grow from a staring length to a stop length. Both the starting and stop length are an input an are given in the numers of elements that the crack is long at the stard and the end. The algorihm requires the crackpath to be known before running the optimization. The algorihm calculates the stress intensity factor for all crack lengths between the starting and stop length. Than it uses the Paris rule to compute the crack growth rate. The fatigue crack growth life can be determined by integrating the inverse of the crack growth rate to the crack length.  Start by using the [example.py](./src_FatigueLive/example.py).

## Prerequisites ##
Python 3 with NumPy, SciPy, matplotlib and cvxopt as the core. The ffmpeg packege is required only when .mp4 movies are created. To simplify the setup Anaconda enviroments (including Spyder) are avalible both for [Window](./anaconda/TopOpt_Windows.yml) and [Linux](./anaconda/TopOpt_Linux.yml).

## Documentation ##
More detailed documentation can be found at [ReadTheDocs](https://topopt.readthedocs.io/en/latest/). This shows the docstrings of all modules, classes and expains how new optimization can be soleved. For that purpose is expains how to define a new losadcase and what the meaning of all the optimization settig are.

## To be implemented ##
Various improvements are proposed:
 1. Adding a faster algebraic multigrid preconditioner with Conjugate Gradient itterative solver.
 2. Creating a conform cubic element mesh for the StressIntensity minimization with a multiresolution implementation as shown in [Higher‐order multi‐resolution topology optimization using the finite cell method](https://doi.org/10.1002/nme.5432) by J.P. Groen, M. Langelaar, O. Sigmund, and M. Ruess 2017. This can both improof the performance of the optimization and acuracy of the stress intensity factors calculated.
 3. Coordinate system transformation to convert the local element coordinate system to the global one.
 4. Use more advanced elements that allow the crackpath to penetrate through elements and the cracktip to exist inside the element.
 5. A variable thickness cracktip.
 6. Add multiloadcase optimisation capebilities to Compliance and Actuator code.
 7. Adding a compliance constraint to the sress intensity factor minimization.
 8. Improved filtering stratigies to control minimum featuresize whithout the appearace sf artifacts.
 
## Special Thanks To ##
 1. Ole Sigmund and Martin Bends&#248;e from the Technical University of Denmark for their contributions to the field of topology optimization
 2. Li-Yi Wei from the University of Hong Kong for making a beautiful [python code](https://github.com/1iyiwei/topopt)
 3. Svanberg, K. (1987). The method of moving asymptotes — a new method for structural optimization. International Journal for Numerical Methods in Engineering, 24(2), 359–373. DOI: [10.1002/nme.162024027](https://doi.org/10.1002/nme.1620240207)
