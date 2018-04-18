Different versions of gurobi give different optimum. 

INSTALLATION INSTRUCTIONS
#To install COBRA, need to do the following first if you have a Linux machine (Ubuntu 14.04 tested):
	For Linux:
		sudo apt-get install swig
	For Mac OSX:
		Need to include brew
		brew install glpk

#To install LibSBML for Ubuntu for the model module, enter from the terminal: 
#	sudo apt-get install python-dev
#	sudo apt-get install libxml2-dev
#	sudo apt-get install zlib1g-dev
#	sudo apt-get install libbz2-dev
#	sudo apt-get install python-pip
#	sudo pip install python-libsbml
	
Gurobi is the linear optimization solver of choice. It is free for academic users. Install the Gurobi MILP solver following the protocol on Gurobi's website: 
	http://user.gurobi.com/download/gurobi-optimizer
In brief, change into the Gurobi directory with the setup.py file and type:
	python setup.py install	
After installing Gurobi, set Gurobi's path accordingly as descibed in Gurobi's installation protocol. Here is an example of an appropriate file path:
	export GUROBI_HOME="/home/eddie/gurobi563/linux64"
	export PATH="${PATH}:${GUROBI_HOME}/bin"
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
Then go back to gurobi.com and download a free academic licencse. Enter your specific grbgetkey command for your license to activate Gurobi.


For people who do not know where the setyp.py is, it is \Library\gurobi***\mac64



Add the path of the install in MATLAB by Add with Suubfolders. To save the path, first you may have to change the read/write privileges for the folder that pathdef.m is in, such as:
		sudo chmod -R 777 /usr/local/MATLAB/R2015b/toolbox/local/

To install ILOG CPLEX:
	Sign up for the Academic Initiative for IBM: https://developer.ibm.com/academic/offers/
	Download IBM ILOG CPLEX Optimization Studio V12.6.1 for Linux x86-64 Multilingual (CNK2KWML)  
	To download Java if you don't already have it, which is needed for the install:
		sudo add-apt-repository ppa:webupd8team/java
		sudo apt-get update
		sudo apt-get install oracle-java7-installer
		sudo update-java-alternatives -s java-7-oracle
	After the download is complete, change into the downloaded directory and type:
		chmod +x CPLEX_OPTIM_STUDIO_12.6.1_LNXX86-.bin
		sudo ./CPLEX_OPTIM_STUDIO_12.6.1_LNXX86-.bin
	Add the path of the install in MATLAB by Add with Subfolders


To install cobrapy for Ubuntu, enter from the terminal:
	sudo apt-get install python-dev libglpk-dev
	sudo pip install cobra


To install mlabwrap:
	Downlaod mlabwrap:
		http://sourceforge.net/projects/mlabwrap/
		Unpack it and change into the direcotry
		If there is a discrepancy between the Matlab license user and root privileges for distribution packages for python, change priviliges for dist-packages:				
			sudo chmod -R ugo+rw /usr/local/lib/python2.7/dist-packages
		install mlabwrap by:
			python setup.py install
		Make the following changes to symbolic links if necessary in /usr/local/MATLAB/R2015b/bin/glnxa64/:
			sudo mv libssl.so.1.0.0 tmp-libssl.so.1.0.0
			sudo ln -s /lib/x86_64-linux-gnu/libssl.so.1.0.0 libssl.so.1.0.0
			sudo mv libcrypto.so.1.0.0 tmp-libcrypto.so.1.0.0
			sudo ln -s /lib/x86_64-linux-gnu/libcrypto.so.1.0.0 libcrypto.so.1.0.0


#To install the most recent version of cobrapy for Ubuntu (at the time cobra 0.4.0b6):
#	sudo pip install cython
#	clone the repository using git (set up an account on git if you do not have one):
#		sudo apt-get install git-all	
#		fork the repository on github.com
#		create a new folder on your local system and chage into that folder
#		from the command line run:
#			git clone https://github.com/YOUR-USERNAME/cobrapy
#		change into the cobrapy directory
#		from the command line run:
#			python setup.py develop --user


To install the COBRA Toolbox for MATLAB:
	First download the libSBML interface for MATLAB:
		http://sourceforge.net/projects/sbml/files/libsbml/5.12.0/stable/Linux/64-bit/MATLAB%20interface/
		I extracted it, and then saved it into Downloads.
		I then added the path for MATLAB.
	Download the COBRA Toolbox:
		https://opencobra.github.io/cobratoolbox/
		Extract the package and change into the following folder: 
			<YOUR_COBRA_ROOT_FOLDER_HERE>/external/toolboxes/SBMLToolbox-4.1.0/toolbox
		Then run:
			installSBMLToolbox			
		Add the path of the install in MATLAB by Add with subfolders.	

To install fastFVA:
	Download glpk from sourceforge:
		http://ftp.gnu.org/gnu/glpk/
	Unpack glpk, and run:
		./configure
		make
		make install
	Add the path for glpk.h in Matlab (It is most likely in usr/local/include)	
	Download fastFVA:
		http://wwwen.uni.lu/lcsb/research/mol_systems_physiology/fastfva
	Unpack fastFVA
	In Matlab, change into the fastFVA/code directory and enter:
		mex -largeArrayDims glpkFVAcc.cpp -lglpk -lm

	I need to figure out how to get CPLEX working in Linux. 
	
		 
To install python-matlab-bridge:
	sudo apt-get install libzmq3
    	sudo pip install pyzmq (installed pyzmq-15.1.0)
    	sudo pip install ipython
    	sudo pip install jupyter
	Download the .tar.gz file from the following website: http://arokem.github.io/python-matlab-bridge/
	Unpack the files, change into the new directory, and run:
		python setup.py install

To install matlab_wrapper:
	sudo pip install matlab_wrapper

Download and install the most recent R and RStudio packages from online for importing and analyzing microarray data sets.


FILE DESCRIPTION for EXAMO
The following are the description of the files:
git/EXAMO-reborn/iMM904_Testing/160202_Model_Conversion_Testing/
	This contains the files necessary to convert any .xml or .mat metabolic reconstruction in the necessary .pkl file format for EXAMO (and/or make desired changes and convert the model to a COBRA compatible .mat format). 
	-160415_import_cobra_model_46_iMM904.py imports the model and generates the output .pkl EXAMO and .mat COBRA files. Below is the program description upon entering python 160416_import_cobra_model_46_iMM904.py -h: 		
		Convert SBML model (.xml or .mat file) into EXAMO pickle file

		positional arguments:
		  s           Necessary variable: SBML file
		  b           Necessary variable: biomass rxn; Append "R_" to the front if the
    		              reaction name does not begin with this
  		  e           Necessary variable: extracellular compartment abbreviation.
              		      Instead of brackets, use underscores (ex: "_e")

		optional arguments:
		  -h, --help  show this help message and exit
		  -l L        lb file
		  -u U        ub file
		  -g G        gene2rxn file
		  -m M        metabolite mapping complexes .pkl file. Must first create .pkl
		              file if using this option
		  -n N        nucleotide conversions .pkl file. Must first create .pkl file if
		              using this option
		  -c C        metabolite to carbon mapping file

	To create a model with metabolite mappings, nucleotide conversion, and carbon mappings as an example, the command would be:
		python 151022_import_cobra_model_46_iMM904.py iMM904_NADcorrected_1127_FTHFLm.xml R_biomass_published _e -m 150723_iMM904_NADcorrected_1127_FTHFLm_metabolite_mappings.pkl -n 150723_iMM904_NADcorrected_1127_FTHFLm_nucleotide_conversions.pkl -c iMM904_NADcorrected_1127_FTHFLm_metabolite_dict.csv > test.txt
	- iMM904_NADcorrected_1127_FTHFLm_genes_genes2rxns.txt is an example if you wanted to change the gene2rx mappings for reactions. The lb and ub files wouuld be formatted the same way, with the boundaries on a new line for each respective reaction in the model.
	- 150723_iMM904_NADcorrected_1127_FTHFLm_metabolite_mappings.pkl is an example of a metabolite mappings dictionary. This was produced by 150723_iMM904_NADcorrected_1127_FTHFLm_metabolite_mappings.py
	- 150723_iMM904_NADcorrected_1127_FTHFLm_nucleotide_conversions.pkl is an example of a nucleotide conversions dictionary. This was produced by 150723_iMM904_NADcorrected_1127_FTHFLm_nucleotide_conversions.py.
	- iMM904_NADcorrected_1127_FTHFLm_metabolite_dict.csv contains the format needed to map the number of carbons to every metabolite and create a carbon balanced model, eliminating nonbalanced reactions if neccessary. 
		*To note, eliminating rxns that are not carbon balanced may result in the model not being solvable for a biomass flux if the metabolites involved are the only way to produce certain metabolites in the biomass reaction.


git/EXAMO-reborn/151113_iMM904_Testing
	This contains the files necessary to run EXAMO.  
	-160720_01_findZeroAndHighFreqRxns_1.py is the first script of the original EXAMO package and determines high frequency and zero frequency reactions. Below is the program description upon entering python 150720_01_findZeroAndHighFreqRxns_1.py:
		_01_findZeroAndHighFreqRxns.py

		positional arguments:
		  s           Necessary variable: SBML file
		  p           Necessary variable: pickle file
		  d           Necessary variable: condition description
		  b           Necessary variable: biomass reaction

		optional arguments:
		  -h, --help  show this help message and exit
		  -e E        Minimum value for fluxes forced to be non-zero (for the
		              implementation of the pathway algorithm created by Shlomi);
		              default value: 1E-1
		  -a A        Activity threshold (above which a flux is considered to be
		              larger than 0); default value: 1E-10
		  -o O        Defined biomass production

	To run the first script with a defined biomass for example for the Aerobic condition for iMM904, the command would be:	
		python 150720_01_findZeroAndHighFreqRxns_1.py iMM904_NADcorrected_1127_FTHFLm.xml iMM904_NADcorrected_1127_FTHFLm.pkl 151006_Aerobic_0.25 R_biomass_published -o 0.2879
	-150723_02_04.py is the second through fourth scirpts of the original EXAMO package. This runs the MBA profile, then adds back the most occurring reactions to the core high frequency reactions until biomass produced, and then the minimziation of fluxes is determiend. Below is the program description upon entering python 150723_02_04.py -h
		_02_04.py

		positional arguments:
		  s                     Necessary variable: SBML file
		  p                     Necessary variable: pickle file
		  d                     Necessary variable: condition description
		  b                     Necessary variable: biomass reaction
		  r                     Necessary variable: number of repetitions of the 2nd-
		                        4th scripts for final flux states

		optional arguments:
		  -h, --help            show this help message and exit
		  -e E                  Minimum value for active flux of reversible reactions
		                        in an irreversible model; default value: 1E-10
		  -a A                  Activity threshold (above which a flux is considered
		                        to be larger than 0); default value: 1E-10
		  -t T                  Activity threshold for finding active reactions;
		                        default value: 1E-10
		  -o O                  Defined biomass production
		  -EXrxns EXRXNS        Pickle file containing extracellular reactions list
		  -EXtrrxns EXTRRXNS    Pickle file containing extracellular transport
		                        reactions list
		  -Othertrrxns OTHERTRRXNS
		                        Pickle file containing other compartmental transport
		                        reactions list
	To run the second through fourth scripts with 5 repetitions defined biomass and adjusted miminum value for active flux, activity threshold, and activity threshold for finding active reactions, the command would be:
		python 150723_02_04.py iMM904_NADcorrected_1127_FTHFLm.xml iMM904_NADcorrected_1127_FTHFLm.pkl 151006_Aerobic_0.25 R_biomass_published 5 -o 0.2879 -e 1E-9 -a 1E-9 -t 1E-9
	-iMM904_EXtrrxns.pkl is aan example of extracellular transport reactions created by iMM904_EXtrrxns.py.
	-iMM904_EXrxns.pkl is an example of extracellular reactions created by iMM904_EXrxns.py.
	-iMM904_Othertrrxns.pkl is an example of other compartmental transport reactions created by iMM904_Othertrrxns.py. 
	-To note, the .pkl and .xml or .mat metabolic reconstructions need to be in the data subfolder, along with the gene rules.csv file. 



git/EXAMO-reborn/Gene_Rules
	This contains the files necessary to create the gene rules for each condition.
	-Gasch_Microarray_Analysis.R analyzes the microarray data from GEO and creates the necessary input file for the gene rules script. Run the script from RStudio.
	-Rintala_Micorarray_Analysis.R analyzes the microarray data from GEO and creates the ncessary input file for the gene rules script. Run the script from RStudio.
	-151007_Rintala_geneClassification.pl creates the gene rules for aerobic and anaerobic conditions. For a threshold of top and bottom 25% of the expressed genes being defined as on or off respectively, run from the command line:
		perl 151007_Rintala_geneClassification.pl 0.25
	-151012_Gasch_geneClassification.pl creates the gene rules for glucose and ethanol conditions. For a threshold of top and bottom 25% of the expressed genes being defined as on or off respectively, run from the command line:
		perl 151012_Gasch_geneClassification.pl 0.25
	-151015_negative_control_geneClassification.pl creates the negative control gene rules. For a threshold of top and bottom 25% of the expressed genes being defined as on or off respectively, run from the command line:
		perl 151015_negative_control_geneClassification.pl 0.25
