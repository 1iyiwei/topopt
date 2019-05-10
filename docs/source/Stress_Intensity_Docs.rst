Stress Intensity Factor Minimization
====================================

In this stress intensity factor minimization a structure with crack is optimized to have minimal crack growth rate.
Thow this works is discussed in :ref:`Stress Intensity Explanation`
An example as how to use the optimization is shown in an example optimization example.py_

.. contents::
   :local:
   :depth: 1

Density Constraints
-------------------

.. automodule:: src_StressIntensity.constraints
    :members:

Load Cases
-------------------

.. automodule:: src_StressIntensity.loads

Parent Load Case
^^^^^^^^^^^^^^^^^

.. autoclass:: src_StressIntensity.loads.Load
    :members:

Child Load Cases
^^^^^^^^^^^^^^^^

.. autoclass:: src_StressIntensity.loads.EdgeCrack
    :members:
    :show-inheritance:

.. autoclass:: src_StressIntensity.loads.DoubleEdgeCrack
    :members:
    :show-inheritance:

.. autoclass:: src_StressIntensity.loads.CompactTension
    :members:
    :show-inheritance:

Finite Element Solvers
----------------------

.. automodule:: src_StressIntensity.fesolvers

Parent Solver
^^^^^^^^^^^^^^

.. autoclass:: src_StressIntensity.fesolvers.FESolver
    :members:

Child Solvers
^^^^^^^^^^^^^^

.. autoclass:: src_StressIntensity.fesolvers.CvxFEA
    :members:
    :show-inheritance:

.. autoclass:: src_StressIntensity.fesolvers.CGFEA
    :members:
    :show-inheritance:

Optimization Module
-------------------

.. automodule:: src_StressIntensity.topopt
    :members:

Plotting Module
-------------------
.. automodule:: src_StressIntensity.plotting

.. autoclass:: src_StressIntensity.plotting.Plot
    :members:

.. autoclass:: src_StressIntensity.plotting.FasterFFMpegWriter
    :members:
    :show-inheritance:
    :inherited-members:

.. _example.py: https://github.com/AJJLagerweij/topopt/blob/master/src_StressIntensity/example.py

