import os
from setuptools import setup, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
#Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
 
setup(name = "TASMANIAN_DEVIL",
      setup_requires=['numpy', 'pandas'],
      install_requires=['numpy', 'scipy', 'pandas', 'patsy', 'cobra==0.8.2'],
      packages = ["TASMANIAN_DEVIL"],
      entry_points = {"console_scripts": ['tas = TASMANIAN_DEVIL.TASMANIAN_DEVIL:main']},
      version = "0.1.0",
      description = "Something Something TASMANIAN_DEVIL",
      long_description=long_description,
      license='GPLv3',
      url = "www.nowhere.com",)

