Global Compliance Minimization Docs
===================================

The total compliance minimization does design structures with maximum stiffness as is discussed at :ref:`Global Compliance Minimization`.
It exists of five module.
An example as how to use the optimization is shown in an example optimization example.py_

.. contents::
   :local:
   :depth: 1


Density Constraints
-------------------

.. automodule:: src_Compliance.constraints
    :members:

Load Cases
-------------------

.. automodule:: src_Compliance.loads

Parent Load Case
^^^^^^^^^^^^^^^^^

.. autoclass:: src_Compliance.loads.Load
    :members:

Child Load Cases
^^^^^^^^^^^^^^^^

.. autoclass:: src_Compliance.loads.HalfBeam
    :members:
    :show-inheritance:
.. autoclass:: src_Compliance.loads.Beam
    :members:
    :show-inheritance:
.. autoclass:: src_Compliance.loads.Canti
    :members:
    :show-inheritance:
.. autoclass:: src_Compliance.loads.Michell
    :members:
    :show-inheritance:
.. autoclass:: src_Compliance.loads.BiAxial
    :members:
    :show-inheritance:

Finite Element Solvers
----------------------

.. automodule:: src_Compliance.fesolvers

Parent Solver
^^^^^^^^^^^^^^
.. autoclass:: src_Compliance.fesolvers.FESolver
    :members:

Child Solvers
^^^^^^^^^^^^^^

.. autoclass:: src_Compliance.fesolvers.CvxFEA
    :members:
    :show-inheritance:
.. autoclass:: src_Compliance.fesolvers.CGFEA
    :members:
    :show-inheritance:

Optimization Module
-------------------

.. automodule:: src_Compliance.topopt
    :members:

Plotting Module
-------------------

.. automodule:: src_Compliance.plotting
    :members:
    :show-inheritance:

.. _example.py: ../../../src_Compliance/example.py
