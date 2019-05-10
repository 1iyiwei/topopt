Fatigue Crack Growth Life Maximization
======================================

Intro

.. contents::
   :local:
   :depth: 1

Density Constraints
-------------------

.. automodule:: src_FatigueLive.constraints
    :members:

Load Cases
-------------------

.. automodule:: src_FatigueLive.loads

Parent Load Case
^^^^^^^^^^^^^^^^^

.. autoclass:: src_FatigueLive.loads.Load
    :members:

Child Load Cases
^^^^^^^^^^^^^^^^

.. autoclass:: src_FatigueLive.loads.EdgeCrack
    :members:
    :show-inheritance:

.. autoclass:: src_FatigueLive.loads.CompactTension
    :members:
    :show-inheritance:

Finite Element Solvers
----------------------

.. automodule:: src_FatigueLive.fesolvers

Parent Solver
^^^^^^^^^^^^^^

.. autoclass:: src_FatigueLive.fesolvers.FESolver
    :members:

Child Solvers
^^^^^^^^^^^^^^

.. autoclass:: src_FatigueLive.fesolvers.CvxFEA
    :members:
    :show-inheritance:

.. autoclass:: src_FatigueLive.fesolvers.CGFEA
    :members:
    :show-inheritance:

Optimization Module
-------------------

.. automodule:: src_FatigueLive.topopt
    :members:

Plotting Module
-------------------
.. automodule:: src_FatigueLive.plotting

.. autoclass:: src_FatigueLive.plotting.Plot
    :members:

.. autoclass:: src_FatigueLive.plotting.FasterFFMpegWriter
    :members:
    :show-inheritance:
    :inherited-members:

.. _example.py: https://github.com/AJJLagerweij/topopt/blob/master/src_FatigueLive/example.py

