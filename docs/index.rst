Welcome to Velocity Tools's documentation!
===========================================

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

.. toctree::
     :maxdepth: 2
     :caption: Contents:

     installation
     usage


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`