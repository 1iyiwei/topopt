Background in Topology Optimization
===================================
+------------------------------------------------------------------+
| .. image:: nstatic/Canti_1.svg                                   |
|                                                                  |
+------------------------------------------------------------------+
| |beginfigref|                                                    |
|                                                                  |
| .. math::                                                        |
|    :label: Example_canti                                         |
|                                                                  |
|    \text{TO example, a cantilever beam with maximum stiffness.}  |
|                                                                  |
| |endfigref|                                                      |
+------------------------------------------------------------------+


This chapter will provide the reader with a basic insight into topology optimization (TO).
TO can alter the layout of the structure. Within a design space it tries to distribute a limited amount of material such that a certain objective is maximized or minimized.
This design space is limited by; the size of the design region, a material constrain, boundary conditions and others.

Here the formulation of a basic algorithm and the problems that can be encountered are disucessed.
It will provide the reader the basic grasp that is required before a change in optimization objective can be discussed.
For that purpose it will introduce a basic example of the TO algorithm that minimizes the global compliance, and thus maximizes stiffness.

This type of TO tries to minimize the global compliance.
It will be the main example algorithm as it has been researched and documented extensively among others by the TopOpt group at the Technical University of Denmark (DTU) [1]_, [2]_, [3]_, [4]_, [5]_.
The goal of the method is to minimize the compliance by distributing the assigned mass. It has to satisfy certain constraints, the volume constrain :math:`V` limits the amount of mass available and the structure should be in equilibrium.
If required, more constraints can be formulated.
One can limit the size of the finest features and take manufacturing limitations in account or introduce a local density constraint to create porous structures which ensures structural stability [6]_.

Different implementations of global compliance TO exist, the one discussed here is based on a gradient method.
Hence, it requires a continuous expression for the compliance as a function of the mass/density distribution.
Therefore, it must allow elements with density values that are between 0 and 1 and it uses a proportional stiffness with penalization method (SIMP) to approximate a discrete 0-1 problem.

.. contents::
   :local:
   :depth: 1

Continuum Formulation
---------------------
The linear elastic optimization for small deformation as presented by N. Olhoff and J.E. Taylor [7]_ is used.
It considers a design region :math:`\Omega` that is in :math:`\boldsymbol{\!R}^2` or :math:`\boldsymbol{\!R}^3` of which a subregion :math:`\Omega^m` is filled with material [1]_.
The optimal topology is reached when the optimal stiffness tensor :math:`\boldsymbol{E}_{ijkl}(\boldsymbol{x})` is found.

As all space within :math:`\Omega^m` is filled an equation of the mass distribution :math:`X` can be formulated as a discrete function,

.. math::

   X(\boldsymbol{x}) = \;\; \begin{cases} 1 \qquad \text{ if } \;\; \boldsymbol{x} \; \in \; \Omega^m \\ 0 \qquad \text{ if } \;\; \boldsymbol{x} \; \in \; \Omega\backslash\Omega^m \end{cases}

This can be used to define the stiffness tensor,

.. math::

   \boldsymbol{E}_{ijkl}(\boldsymbol{x}) = X(\boldsymbol{x})\boldsymbol{\overline{E}}_{ijkl}

in terms of this mass distribution function and the constant rigidity tensor :math:`\boldsymbol{\overline{E}}_{ijkl}`.
The constant rigidity tensor is function of the material properties only.
As :math:`X` is a discrete function all admissible tensors are discrete and thus the optimization problem has a discrete valued parameter function.

The amount of work due of the deformation :math:`\boldsymbol{u}` can be calculated by with a virtual work method.
With the standard linearized strain formulation this results in,

.. math::

   l(\boldsymbol{u}) = \int_{\Omega}\boldsymbol{fu}\text{ d}\Omega + \int_{\Gamma_T} \boldsymbol{tu} \text{ d}\Gamma_T

A bi-linear energy equation with virtual work :math:`a(\boldsymbol{u},\hat{\boldsymbol{u}})` is formulated,

.. math::

   a(\boldsymbol{u},\hat{\boldsymbol{u}}) =\int_{\Omega} \boldsymbol{E}_{ijkl}\boldsymbol{\varepsilon}_{kl}(\boldsymbol{u})\boldsymbol{\varepsilon}_{ij}(\hat{\boldsymbol{u}})\text{ d}\Omega

:math:`\hat{\boldsymbol{u}}` is an arbitrary kinematically admissible deformation.
Equilibrium is ensured when :math:`l(\hat{\boldsymbol{u}}) = a(\boldsymbol{u}, \hat{\boldsymbol{u}})` is satisfied for all admissible deformations :math:`\hat{\boldsymbol{u}}`.

As minimizing the work, due to the traction forces for a given load, minimizes the deformation of a structure the problem can be formulated as:

.. math::

    \min_{\Omega^m} \;\;& l(\boldsymbol{u}) \\
    &\begin{array}{llll}
    \text{s.t. :} & a(\boldsymbol{u},\hat{\boldsymbol{u}}) = l(\hat{\boldsymbol{u}}) \\
    & \int_{\Omega} X(\boldsymbol{x}) \text{d}\Omega \; = \; \text{ Vol}(\Omega^m) \; \leq \; V
    \end{array}


Discretization
---------------
To solve the continuum problem of the previous section it is discretized into a finite element analysis with :math:`N` elements:

.. math::

   \min_{X_1, X_2, \dots, X_N} \;\: & c = \boldsymbol{f}^T \boldsymbol{u}\\
   &\hspace{-0.6cm}\begin{array}{llll}
   \text{s.t. :} & \boldsymbol{Ku} = \boldsymbol{f} \\
   & \displaystyle\sum^N_{e=1} v_eX_e \; \leq \; V \\
   & X_e \in \{0, 1\} \;\;\; \forall \;\;\; e \in \{1, 2, \dots, N\}\\
   \text{where :} & \boldsymbol{K} = \displaystyle\sum_{e=1}^{N}\boldsymbol{K}_e(X_e, \overline{E})
   \end{array}

it shows that the element stiffness matrix :math:`\boldsymbol{K_e}` depends on the element material value :math:`X_e` and the material stiffness :math:`\overline{E}`.
The problem becomes unstable towards the element type and mesh when the discrete formulations of density are used.
Such a distribution problem generally has no solution [8]_, [9]_. Iterative search methods would not work because they require the calculation of gradients.
Therefore, the problem is changed so that the density becomes a continuous equation ranging from 0 to 1.

.. math::

   0 \leq X_e \leq 1

This method would result in a design with intermediate values.
Although this makes sense for variable thickness plate design, see the work of M.P. Rossow and J.E. Taylor [10]_, for discrete topology design loses its direct physical representation.
There is either material or there is not, intermediate values are meaningless.
Adding a penalization that reduces the effectiveness of intermediate values results in a formulation that suppresses these intermediate values.
The method used here, developed by E. Andreassen [5]_, is derived from the classical penalized proportional stiffness method (SIMP) [1]_, [3]_.
Here :math:`E_{\min}` is a small artificial stiffness used to avoid elements with zero stiffness as that could make the FEA unstable.

.. math::

   \boldsymbol{E}_{ijkl}(\boldsymbol{x}) = \boldsymbol{E}_{ijkl, \min} + X(\boldsymbol{x})^p\left(\boldsymbol{\overline{E}}_{ijkl} - \boldsymbol{E}_{ijkl, \min}\right)

When :math:`p > 1` the intermediate density values are less effective as there stiffness is low in comparison to the volume occupied. When :math:`p` is sufficiently large, generally :math:`p\geq3`, the design converges to a solution that is close to a discrete (0-1) design.

Sensitivity analysis and MMA
-----------------------------
The main focus on developing a robust and stable algorithm is the update scheme.
The MMA scheme was chosen as it proofed to be very effective for this type of optimization [3]_.
MMA is an efficient method meant for non-linear non-convex problems that approaches those problems by generating purely convex sub-problems, based on the gradient information.
It can be used to iterative solve the optimization problem.

The gradient of one element in the discretized form is :math:`\partial c/\partial X_e`.
This derivative does not have to be explicitly calculated as the problem is self adjoint.
This is used by  the following proof. It starts with a new formulation of the work, the difference is the zero term at the end.
Again :math:`\hat{\boldsymbol{u}}` is any arbitrary admissible deformation [3]_.

.. math::

   c = \boldsymbol{f}^T \boldsymbol{u} - \hat{\boldsymbol{u}}^T\left( \boldsymbol{Ku} - \boldsymbol{f} \right)

taking the derivative to the density leads to:

.. math::

   \frac{\partial c}{\partial X_e} = \left( \boldsymbol{f}^T - \hat{\boldsymbol{u}}^T\boldsymbol{K} \right) \frac{\partial \boldsymbol{u}}{\partial X_e} - \hat{\boldsymbol{u}}^T \frac{\partial\boldsymbol{K}}{\partial X_e}\boldsymbol{u}

when :math:`\hat{\boldsymbol{u}}` satisfies the adjoint equation it becomes:

.. math::

   \frac{\partial c}{\partial X_e} = & - \hat{\boldsymbol{u}}^T	\frac{\partial\boldsymbol{K}}{\partial X_e}\boldsymbol{u} \\
   & \text{when} \hspace{0.5cm} \boldsymbol{f}^T - \hat{\boldsymbol{u}}^T\boldsymbol{K} = 0

Satisfying this adjoint equation is simple, just choose :math:`\hat{\boldsymbol{u}} = \boldsymbol{u}`.
The derivative of the stiffness matrix to the density of an element can be derived leading to the final expression of the gradient:

.. math::

   \frac{\partial c}{\partial X_e} = - pX_e^{p-1}\boldsymbol{u}^T\boldsymbol{K}_e\boldsymbol{u}

MMA approaches the problem with multiple convex approximations around the expansion point (current iteration).
The goal here is to find the optimal density distribution of the current iteration where the influence of the densities is approximated with a convex function.
This approximation is based on the sensitivity and some information of previous iterations. Solving these convex equation can be done by various basic algorithms.
The obtained optimum is not the real optimum of the optimization problem as the convex function used is only an approximation of the real problem.
However, it is a step into the direction of the real optimum. The obtained density distribution is then used as an input of the next iteration [3]_ (pp. 19-21).
The optimization of this local problem must meet all the constraints. This means that the updated design has to meet the global volume constraint.

The MMA will approximate the compliance at iteration :math:`k`.
Here :math:`X^k` is a vector with the densities of all elements at the current iteration.
A description on the calculations of :math:`U_e` and :math:`L_e` follows later. The method was developed by K. Svansberg [11]_.

.. math::

   c &\approx c^k + \sum^{N}_{e =1}\left( \frac{r_e}{U_e- X_e} + \frac{s_e}{X_e - L_e} \right) \\
   &\begin{array}{ll}
   \text{where: } &  r_e = \begin{cases} 0 & \text{ if } \;  \frac{\partial c}{\partial X_e} \leq 0 \\ \left(U_e - X_e^{k}\right)^2\frac{\partial c}{\partial X_e}\phantom{-} & \text{ if } \; \frac{\partial c}{\partial X_e} > 0 \end{cases} \\
   & s_e = \begin{cases} 0 & \text{ if } \;  \frac{\partial c}{\partial X_e} \geq 0 \\ -\left(X_e^{k} - L_e\right)^2\frac{\partial c}{\partial X_e} & \text{ if } \; \frac{\partial c}{\partial X_e} < 0 \end{cases} \\
  \end{array}

That all the density sensitivities are negative can be derived from adjoint sensitivity equation. This simplifies the expression and resulted in:

.. math::

   c \approx c^k + \sum^{N}_{e =1}-\frac{\left( X_e^{k}-L_e\right)^2}{X_e-L_e}\frac{\partial c}{\partial X_e}

Then the optimization, on :math:`X_e`, used in this iteration is defined as:

.. math::

   \min_{X_1, X_2, \dots, X_N} \;\; & c^k - \sum^{N}_{e =1}\frac{\left(X_e^{k} - L_e\right)^2}{X_e- L_e}\frac{\partial c}{\partial X_e}\\
   &\begin{array}{llll}
   \text{s.t. :} & \displaystyle\sum^{N}_{e=1}v_eX_e \; \leq \; V \\
   & 0 \geq X_e \geq 1 \;\;\; \forall \;\;\; e \in \{1, 2, \dots, N\}
   \end{array}

here the moving asymptote, :math:`L_e`, can be varied and is chosen to improve convergence and stability, choosing this wisely is important.
In general the goal is to stabilize the process when it is oscillating, i.e. moving the asymptote closer.
Or to relax the problem when it is monotone, i.e. moving the asymptote further and thus causing larger steps to be taken at that iteration.
This can be done by including the behavior of previous iterations or calculating the second derivative of the optimization objective to the design variables.
Several implementations exist, they are tuned to work for specific problems [11]_, [12]_.

The update scheme minimizes the local approximation to decide on the new densities. Starting with the minimalization of the Lagrange function:

.. math::

   \mathcal{L} = c^k  - \sum^{N}_{e =1}\frac{\left(X_e^{k} - L_e\right)^2}{X_e- L_e}\frac{\partial c}{\partial X_e} + \Lambda\left( \sum_{e=1}^N v_eX_e -V \right) +  \sum_{e=1}^N \lambda^-_e\left(X_e - 0 \right) + \sum_{e=1}^N \lambda^+_e\left(1 - X_e \right)

This separable and purely convex problem can be solved by a range of algorithms. It can easily be changed into a formulation with other or more constraints.


Filtering Techniques
--------------------
Filtering the sensitivities was proposed by O. Sigmund [13]_ .
The method is derived from image processing and uses a normalized convolution filter to blur the figure.
The density distribution :math:`X_e` and the gradient can be interpreted as a figure with gray scale pixels.
The gradient itself is not filtered, but the gradient multiplied by the densities is filtered before the update scheme decides on the densities of the next iteration [14]_, [15]_.

The sensitivity filter can be described as,

.. math::

   \widehat{\frac{\partial C}{\partial X_k}} =& \dfrac{1}{X_k \sum_{i=1}^{N}H_i}\sum_{i=1}^{N} \; H_i \; X_i \; \frac{\partial l(\boldsymbol{u})}{\partial X_i} \\
   & H_i = \begin{cases} r_{min} - \text{dist}(k,i) & \text{if} \hspace{5mm} \text{dist}(k,i) < r_{min}\\
   0 &  \text{if} \hspace{5mm} \text{dist}(k,i) \geq r_{min}
   \end{cases}

where :math:`k` is the element to be filtered.
The value of the filtered compliance density gradient at element :math:`i` is depended on three main things, the density, density gradient and the distance to the surrounding nodes :math:`i`.
All nodes that fall within radius :math:`r_{min}` are contributing but the further the node is the lower its contribution. Note that the filter is normalized by dividing it by :math:`\sum\hat{H}_i`.
There is limited understanding why this filter works, there is no physical or theoretical basis for it. From experience, it was simply observed that it works well.

**Sensitivity filtering figure**

**Ref to figure** show the same simulations. The only difference is that the simulations is that the are filtered.
It was observed that scaling the filter size :math:`r_{min}` with the resolution results in similar designs.
The main difference between the designs is that higher resolution simulations result in a smoother structure.
But filtering this way leads to less discrete designs. Larger filters cause more pixels to have intermediate density values.
Three solutions do exist; lowering the filter size for the last couple of iterations, increasing the SIMP penalty factor or applying extra post processing steps.

Another filter that can be considered is the linear density filter which was proposed by T.E. Bruns, D.A. Tortorelli and B. Bourdin [16]_, [17]_. Here the blur filter,

.. math::

   \widehat{X_e} =& \dfrac{1}{\sum_{i=1}^{N}H_i}\sum_{i=1}^{N} \; H_i \; X_i \\
   & H_i = \begin{cases} r_{min} - \text{dist}(k,i) & \text{if} \hspace{5mm} \text{dist}(k,i) < r_{min}\\
   0 &  \text{if} \hspace{5mm} \text{dist}(k,i) \geq r_{min}
   \end{cases}

is applied directly on the densities.
These filtered densities, :math:`\widehat{X_e}`, are used in the FEA and SA.
This means that the design variables :math:`X_e` lose there physical meaning as the FEA gives it the relation to reality, therefore the final geometry should be based on the filtered densities [18]_.

A comparison between **Ref to figure** shows that filtering the densities suppresses the finer features well.
Comparing the performance difference of the sensitivity and density filters is difficult.
Many criteria can be used such as, computational effort, how discrete the final design is, the magnitude of the final compliance and whether the volume constrained is still maintained.
A small comparison was made by O. Sigmund [18]_.
The performance of the filters depends greatly on the design case used.
The paper clearly shows that better filters exist then those presented in this communication however as the density and sensitivity filters are computational efficient and simple to implement they were chosen as the basic filters used in the code.

**Density filtering figure***

Computational Implementation
----------------------------
The iterative implementation of topology optimization as proposed by M. Beckers, [19]_ or M.P. Bendsøe and O. Sigmund [3]_ are similar.
It exists out of three parts, initialization, optimization and post processing.
The flowchart for the methods used in this communication can be found in **Ref to Flowchart**.

**FLOWCHART!!!**

In the initialization phase the problem is set up.
It defines the design domain, the loading conditions, the initial design and generates the finite element mesh that will be used in the optimization phase.

The optimization phase is the iterative method that solves the topology problem.
It will analyze the current design with a FEA. After which it will calculate the sensitivity of the global compliance to the density of each element, this is the local gradient of which the calculation is discussed before
The Method of Moving Asymptotes (MMA), developed by K. Svanberg [11]_, is used to formulate a simplified convex approximation of the problem which is optimized to formulate the updated design.
These steps are performed in a loop until the design is converged, i.e. when the change in design between two iterations becomes negligible.

Post processing is required to remove the last elements with intermediate values and generate a shape out of the design, for example a CAD or STL file.
This algorithm will not contain any of the post processing steps.
The code used in this communication simply plots the final shape and load case.

Changing the Objective
----------------------
Topology optimization can be used for several objectives; classical examples are, truss structure design, antenna/microphone design, heat convection problems [3]_, [20]_ and MEMS actuator designs [2]_, [3]_, [21]_, [22]_. In general all these TO algorithms approach the optimization as a material distribution problem within a design space with a resource constraint witch is solved with an iterative gradient method.

When changing the objective and/or problem one should start with a formulation of the problem which consists of the objective, variables and constraints. Then the changes should be made in the calculation of the objective and sensitivity. Important therefor is the method used to link the optimization variables to the objective, in the case of compliance minimization it consists of the variables to density formulation (SIMP \cref{eq:SIMP_Lit}) and the FEA that links stiffness to compliance. Beneficial would be a (self) adjoint formulation because it allows for an efficient calculation of the sensitivities. The parts of the method that are unlikely to change are; the overall methodology, described in \cref{fig:Flowchart_Lit}, the method of moving asymptotes and its update scheme.

Sometimes optimization objectives are formulated in the form of several sub objectives resulting in multi objective optimization formulations.
Optimizing for multiple objectives or load cases at once is common. For most structures several considerations, such as costs, weight and strength are taken in account. In addition do most structures experience multiple load-cases during their life. Several TO algorithms have been developed for this purpose. The most basic methods will be discussed here.

**Flowchart 2**

The method sets up multiple FEA as shown in **Ref to Flowchart 2**.
Then the total objective will be linked to sub objectives.
For instance the goal might be to minimize the compliance due to :math:`n` load cases.
One could formulate the total objective (:math:`O`) as the weighted sum of the compliance of all load cases,

.. math::

   O = \sum_{i = 1}^{n} w_i c_i

resulting in a gradient function that can be formulated as,

.. math::

   \frac{\partial O}{\partial X_e} = \sum_{i = 1}^{n} w_i \frac{\partial c_i}{\partial X_e}

Another example can be made with a similar method. Assume that adding up the objective is not what is wanted but that the goal is to prohibit two different failure modes.
Hence, the design update is based on the most critical case resulting in objective,

.. math::

   O = \max \left( o_1, o_2, \dots, o_n \right)

An example of such a formulation can be found in the TO based damage tolerance optimization algorithm presented by Z. Kang, P. Liu and M. Li [23]_.
Where they optimize geometries for the most cricital crack in every iteration. The sensitivities can then be formulated as:

.. math::

   \begin{align}
   \frac{\partial O}{\partial X_e} =& \sum_{i = 1}^{n} s_i \frac{\partial o_i}{\partial X_e} \\
   & \text{where} \quad s_i = \begin{cases}
   1 \quad \text{if} \quad o_i  = O\\
   0 \quad \text{if} \quad o_i \neq O
   \end{cases}
   \end{align}

These basic multiple load case algorithms can be summarized in the flowchart shown in **Ref to Flowchart 2**.
In general the FEA requires most of the computational time therefore the method as shown here is computationally inefficient.
More advanced algorithms have been developed but these are outside the scope of this communication [24]_, [25]_.

References
----------

.. [1]  M. P. Bendsøe, “`Optimal shape design as a material distribution problem <https://www.doi.org/10.1007/BF01650949>`_,” Struct. Optim., vol. 1, no. 4, pp. 193–202, Dec. 1989.
.. [2]  O. Sigmund, “`A 99 line topology optimization code written in matlab <https://www.doi.org/10.1007/s001580050176>`_,” Struct. Multidiscip. Optim., vol. 21, no. 2, pp. 120–127, 2001.
.. [3]  M. P. Bendsøe and O. Sigmund, `Topology Optimization <https://www.doi.org/10.1007/978-3-662-05086-6>`_. Berlin, Heidelberg: Springer Berlin Heidelberg, 2004.
.. [4]  B. S. Lazarov and O. Sigmund, “`Filters in topology optimization based on Helmholtz-type differential equations <https://www.doi.org/10.1002/nme.3072>`_,” Int. J. Numer. Methods Eng., vol. 86, no. 6, pp. 765–781, May 2011.
.. [5]  E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov, and O. Sigmund, “`Efficient topology optimization in MATLAB using 88 lines of code <https://www.doi.org/10.1007/s00158-010-0594-7>`_,” Struct. Multidiscip. Optim., vol. 43, no. 1, pp. 1–16, Jan. 2011.
.. [6]  J. Wu, N. Aage, R. Westermann, and O. Sigmund, “`Infill Optimization for Additive Manufacturing—Approaching Bone-Like Porous Structures <https://www.doi.org/10.1109/TVCG.2017.2655523>`_,” IEEE Trans. Vis. Comput. Graph., vol. 24, no. 2, pp. 1127–1140, Feb. 2018.
.. [7]  N. Olhoff and J. E. Taylor, “`On Structural Optimization <https://www.doi.org/10.1115/1.3167196>`_,” J. Appl. Mech., vol. 50, no. 4b, p. 1139, 1983.
.. [8]  G. Strang and R. V. Kohn, “`Optimal design in elasticity and plasticity <https://www.doi.org/10.1002/nme.1620220113>`_,” Int. J. Numer. Methods Eng., vol. 22, no. 1, pp. 183–188, Jan. 1986.
.. [9]  R. V. Kohn and G. Strang, “`Optimal design and relaxation of variational problems, I <https://www.doi.org/10.1002/nme.1620220113>`_,” Commun. Pure Appl. Math., vol. 39, no. 1, pp. 113–137, 1986.
.. [10]  M. P. Rossow and J. E. Taylor, “`A Finite Element Method for the Optimal Design of Variable Thickness Sheets <https://www.doi.org/10.2514/3.50631>`_,” AIAA J., vol. 11, no. 11, pp. 1566–1569, Nov. 1973.
.. [11]  K. Svanberg, “`The method of moving asymptotes - a new method for structural optimization <https://www.doi.org/10.1002/nme.1620240207>`_,” Int. J. Numer. Methods Eng., vol. 24, no. 2, pp. 359–373, Feb. 1987.
.. [12]  K. Svanberg, “MMA and GCMMA – two methods for nonlinear optimization,” Stockholm, Sweden, 2007.
.. [13]  O. Sigmund, “Design of Material Structures Using Topology Optimization,” PHD thesis, 1994, pp. 72-75.
.. [14]  O. Sigmund, “`On the design of compliant mechanisms using topology optimization <https://www.doi.org/10.1080/08905459708945415>`_,” Mech. Struct. Mach., vol. 25, no. 4, pp. 493–524, 1997.
.. [15]  O. Sigmund and J. Petersson, “`Numerical instabilities in topology optimization: A survey on procedures dealing with checkerboards, mesh-dependencies and local minima <https://www.doi.org/10.1007/BF01214002>`_,” Struct. Optim., vol. 16, no. 1, pp. 68–75, Aug. 1998.
.. [16]  T. E. Bruns and D. A. Tortorelli, “`Topology optimization of non-linear elastic structures and compliant mechanisms <https://www.doi.org/10.1016/S0045-7825(00)00278-4>`_,” Comput. Methods Appl. Mech. Eng., vol. 190, no. 26–27, pp. 3443–3459, Mar. 2001.
.. [17]  B. Bourdin, “`Filters in topology optimization <https://www.doi.org/10.1002/nme.116>`_,” Int. J. Numer. Methods Eng., vol. 50, no. 9, pp. 2143–2158, Mar. 2001.
.. [18]  O. Sigmund, “`Morphology-based black and white filters for topology optimization <https://www.doi.org/10.1007/s00158-006-0087-x>`_,” Struct. Multidiscip. Optim., vol. 33, no. 4–5, pp. 401–424, Feb. 2007.
.. [19]  M. Beckers, “`Topology optimization using a dual method with discrete variables <https://www.doi.org/10.1007/BF01197709>`_,” Struct. Optim., vol. 17, no. 1, pp. 14–24, Feb. 1999.
.. [20]  S. Turteltaub, “`Functionally graded materials for prescribed field evolution <https://www.doi.org/10.1016/S0045-7825(01)00408-X>`_,” Comput. Methods Appl. Mech. Eng., vol. 191, no. 21–22, pp. 2283–2296, Mar. 2002.
.. [21]  O. Sigmund, “`Design of multiphysics actuators using topology optimization – Part I: One-material structures <https://www.doi.org/10.1016/S0045-7825(01)00251-1>`_,” Comput. Methods Appl. Mech. Eng., vol. 190, no. 49–50, pp. 6577–6604, Oct. 2001.
.. [22]  O. Sigmund, “`Design of multiphysics actuators using topology optimization – Part II: Two-material structures <https://www.doi.org/10.1016/S0045-7825(01)00252-3>`_,” Comput. Methods Appl. Mech. Eng., vol. 190, no. 49–50, pp. 6605–6627, Oct. 2001.
.. [23]  Z. Kang, P. Liu, and M. Li, “`Topology optimization considering fracture mechanics behaviors at specified locations <https://www.doi.org/10.1007/s00158-016-1623-y>`_,” Struct. Multidiscip. Optim., vol. 55, no. 5, pp. 1847–1864, May 2017.
.. [24]  K. A. James, J. S. Hansen, and J. R. R. A. Martins, “`Structural topology optimization for multiple load cases using a dynamic aggregation technique <https://www.doi.org/10.1080/03052150902926827>`_,” Eng. Optim., vol. 41, no. 12, pp. 1103–1118, 2009.
.. [25]  E. Nutu, “Multiple load case topology optimization based on bone mechanical adaptation theory,” UPB Sci. Bull. Ser. D Mech. Eng., vol. 77, no. 4, pp. 131–140, 2015.


.. |beginfigref| raw:: latex

                     \begin{minipage}{\textwidth}


.. |endfigref| raw:: latex

                   \end{minipage}

