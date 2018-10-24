#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import sys

# import ah_bootstrap
from setuptools import setup

PACKAGENAME = "velocity_tools"
DESCRIPTION = ('Tools to help analysis of velocity cubes/field')
AUTHOR = 'Jaime E. Pineda'
AUTHOR_EMAIL = 'jpineda@mpe.mpg.de'
LICENSE = "MIT"
URL = "https://github.com/jpinedaf/velocity_tools.git"

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.1.dev1'

here = os.path.abspath(os.path.dirname(__file__))
def read(fname):
    return open(os.path.join(here, fname)).read()

setup(
	name = PACKAGENAME,
    version = VERSION,
    description = DESCRIPTION,
    long_description=read('README.md'),
    packages = ['velocity_tools'],
    install_requires = ['astropy>=2.0'],
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    license = LICENSE,
    url = URL,
    )
      # scripts=scripts,
      # long_description=LONG_DESCRIPTION,
      # cmdclass=cmdclassd,
      # zip_safe=False,
      # use_2to3=True,
      # entry_points=entry_points,
      # **package_info
# )
