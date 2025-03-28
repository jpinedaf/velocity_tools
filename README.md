# velocity_tools
![PyPI - License](https://img.shields.io/pypi/l/velocity_tools?color=red)
[![codecov](https://codecov.io/gh/jpinedaf/velocity_tools/graph/badge.svg?token=4JFPKTRSX0)](https://codecov.io/gh/jpinedaf/velocity_tools)

Repository with a compilation of tools for velocity analysis
It includes:
- Calculation of Keplerian rotation velocity map.
- Deprojection of relative coordinates for a given inclination and rotation angles and an arbitrary center. 
- Calculation of velocity gradient, assuming solid velocity rotation. 
- Streamline model using [Mendoza et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009MNRAS.393..579M) formalism.


## Credits
---

This is developed by:

- Jaime E Pineda ([@jpinedaf](http://github.com/jpinedaf))

with contibutions from 
- M. Teresa Valdivia-Mena ([@tere-valdivia](http://github.com/tere-valdivia))


## Dependencies
---

- astropy (>=5.0)
- scipy (>=1.7)
- numpy (>=1.21)
- matplotlib (>=3.4)
- scikit-image (>=0.9)
- spectral-cube (>=0.5)
- radio-beam (>=0.3)