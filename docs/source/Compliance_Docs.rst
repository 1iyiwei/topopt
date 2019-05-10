Maximum Compliance Actuators
============================

Intro

.. contents::
   :local:
   :depth: 1

Density Constraints
-------------------

.. automodule:: src_Actuator.constraints
    :members:

Load Cases
-------------------

.. automodule:: src_Actuator.loads

Parent Load Case
^^^^^^^^^^^^^^^^^

.. autoclass:: src_Actuator.loads.Load
    :members:

Child Load Cases
^^^^^^^^^^^^^^^^

.. autoclass:: src_Actuator.loads.Inverter
    :members:
    :show-inheritance:

Finite Element Solvers
----------------------

.. automodule:: src_Actuator.fesolvers

Parent Solver
^^^^^^^^^^^^^^

.. autoclass:: src_Actuator.fesolvers.FESolver
    :members:

Child Solvers
^^^^^^^^^^^^^^

.. autoclass:: src_Actuator.fesolvers.CvxFEA
    :members:
    :show-inheritance:

.. autoclass:: src_Actuator.fesolvers.CGFEA
    :members:
    :show-inheritance:

Optimization Module
-------------------

.. automodule:: src_Actuator.topopt
    :members:

Plotting Module
-------------------
.. automodule:: src_Actuator.plotting

.. autoclass:: src_Actuator.plotting.Plot
    :members:

.. autoclass:: src_Actuator.plotting.FasterFFMpegWriter
    :members:
    :show-inheritance:
    :inherited-members:

.. _example.py: https://github.com/AJJLagerweij/topopt/blob/master/src_Actuator/example.py

