INSTALLATION
============

This installation will only work with Python 2.7 or Python 3.6. We highly recommend using a conda environment, as this will manage the pip versions better than a virtual environment when running setup.py. 

To install COBRA, need to do the following first if you have a Linux machine (Ubuntu 14.04 tested):
	- First install pip
	- pip install numpy scipy
	- For Linux:
		* sudo apt-get install swig
		* For pythhon 2.7:
			+ sudo apt-get install libglpk-dev
			+ sudo apt-get install glpk-utils
	- For Mac OSX:
		* Need to install brew:
			+ https://docs.brew.sh/Installation
		* brew install glpk

Downlaod anaconda from the website: 
	- First install Anaconda with Python 2.7 or Python 3.6:
		* https://www.anaconda.com/download/
	- Create a conda environment:
		* https://conda.io/docs/user-guide/tasks/manage-environments.html	
	- Then install Gurobi using conda:	
		* conda install gurobipy
	- Then install a Gurobi academic key from their website usirng grbketkey

Then run test commands found in README.md 

