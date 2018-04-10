import os
from setuptools import setup, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
#Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
 
setup(name = "TASMANIAN-DEVIL",
      setup_requires=['numpy', 'pandas'],
      install_requires=['numpy', 'scipy', 'pandas', 'patsy', 'cobrapy==0.8.2', 'gurobipy'],
      packages = ["TASMANIAN-DEVIL"],
      entry_points = {"console_scripts": ['tas = TASMANIAN-DEVIL.TASMANIAN-DEVIL:main']},
      version = "0.1.0",
      description = "Something Something TASMANIAN-DEVIL",
      long_description=long_description,
      license='GPLv3',
      url = "www.nowhere.com",)

