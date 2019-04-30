Compliance Minimization
=======================
intro 1


Density Constraints
-------------------
.. automodule:: src_Compliance.constraints
    :members:

Load Cases
-------------------
.. automodule:: src_Compliance.loads

Parrent Load Case
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

Parrent Solver
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
