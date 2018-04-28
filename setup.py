import os
from setuptools import setup, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
#Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
 
setup(name = "TASMANIAN_DEVIL",
      setup_requires=['numpy', 'scipy<1.0'],
      install_requires=['numpy', 'scipy<1.0', 'pandas', 'pycrypto', 'argparse', 'cobra==0.8.2', 'sympy', 'python-libsbml'],
      packages = ["TASMANIAN_DEVIL"],
      entry_points = {"console_scripts": ['tas = TASMANIAN_DEVIL.TASMANIAN_DEVIL:main']},
      version = "0.1.0",
      description = "TASMANIAN DEVIL: a software package for classifying gene activity from omics data sets, simplifying metabolic networks, and visualizing the estimated phenotypic fluxes of nutrients",
      long_description=long_description,
      license='GPLv3',
      url = "https://tlawrence3.github.io/Tasmanian-Devil/build/html/index.html",)

