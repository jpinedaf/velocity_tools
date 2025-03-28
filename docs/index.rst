Welcome to Velocity Tools's documentation!
==========================================
.. image:: https://img.shields.io/pypi/v/velocity-tools
   :target: https://pypi.org/project/velocity-tools/
   :alt: PyPI
.. image:: https://img.shields.io/pypi/pyversions/velocity-tools
   :target: https://pypi.org/project/velocity-tools/
   :alt: Python version
.. image:: https://img.shields.io/github/license/jpinedaf/velocity_tools

The ``velocity_tools`` is a Python package with the aim of aiding the 
analysis of velocity cubes and velocity fields. 
It is developed to support the different kinematical analysis of 
cubes and some sections of the package can improve the analysis 
of data cubes with `SpectralCube <https://spectral-cube.readthedocs.io>`_.


The package includes the following functionalities:

* Calculation of Keplerian rotation velocity map.
* Calculation of offsets and sky rotation for a given center.
* Calculation of average profiles for a given array of measurements.
* Masking of data cubes based on velocity and linewidth.
* Deprojection of relative coordinates for a given inclination and rotation angles and an arbitrary center. 
* Calculation of velocity gradient, assuming solid velocity rotation. 
* Streamline model using `Mendoza et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.393..579M>`_ formalism.


Getting Started
^^^^^^^^^^^^^^^
.. toctree::
     :maxdepth: 2
     :includehidden:

     installation.rst

Velocity gradient calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
     :maxdepth: 2
     :caption: Velocity gradient calculation

     vgrad.rst

Offsets and Sky rotation calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
     :maxdepth: 2
     :caption: Offsets and Sky rotation calculation

     offset_rotation.rst

Velocity fields
^^^^^^^^^^^^^^^
.. toctree::
     :maxdepth: 2
     :caption: Velocity fields

     Velocity_fields.rst

Streamline model
^^^^^^^^^^^^^^^^
.. toctree::
     :maxdepth: 2
     :caption: Streamline model

     streamline.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`