Installation
============


We use Gurobi as our mathematical programming solver. Because of this, TASMANIAN DEVIL requires Python 2.7 or Python 3.6. Currently, we only support and strongly recommend installing TASMANIAN DEVIL using Anaconda.

- Install the current version of Gurobi
	+ http://www.gurobi.com/downloads/gurobi-optimizer

* Obtain a Gurobi license (free academic licenses are available)
	+ http://www.gurobi.com/downloads/licenses/license-center

* Install Anaconda with Python 2.7 or Python 3.6:
	+ https://www.anaconda.com/download/

* After installing Anaconda open a terminal window (we will be using the terminal for the remainder of the installation process) and run the following commands to set up additional channels::

	conda config --add channels defaults
	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda config --add channels http://conda.anaconda.org/gurobi

* Create and activate a new Anaconda environment with the version of Python that you downloaded Anaconda with. Note that you must activate your conda environment every time you want to use TASMANIAN-DEVIL::

	conda create -n tas python=3.6 scipy=0.19.1 gurobi glpk
	source activate tas

* Activate your Gurobi installations using your license key (obtained through Gurobi's website) using the command below::

	grbketkey <your-license-key>

* Clone TASMANIAN-DEVIL from GitHub and install using setup script::

	git clone https://github.com/tlawrence3/Tasmanian-Devil
	cd Tasmanian-Devil
	python setup.py install
