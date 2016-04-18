INSTALLATION INSTRUCTIONS
To install LibSBML for Ubuntu, enter from the terminal: (pip install python-libsbml)
	sudo apt-get install python-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libz-dev
	sudo apt-get install libbz2-dev
	sudo apt-get install python-pip
	sudo pip install python-libsbml
	
Install the Gurobi MILP solver following the protocol on Gurobi's website: 
	http://user.gurobi.com/download/gurobi-optimizer
In brief, change into the Gurobi directory with the setup.py file and type:
	python setup.py install	
After installing Gurobi, set Gurobi's path accordingly as descibed in Gurobi's installation protocol. Here is an example of an appropriate file path:
	export GUROBI_HOME="/home/eddie/gurobi563/linux64"
	export PATH="${PATH}:${GUROBI_HOME}/bin"
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
Then go back to gurobi.com and download a free academic licencse. Enter your specific grbgetkey command for your license to activate Gurobi.
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


FILE DESCRIPTION
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
	- iMM904_NADcorrected_1127_FTHFLm_genes_genes2rxns.txt is an example if you wanted to change the gene2rx mappings for reactions. The lb and ub files wouuld be formatted the same way, with the boundaries on a new line for each respective reaction in the model.
	- 150723_iMM904_NADcorrected_1127_FTHFLm_metabolite_mappings.pkl is an example of a metabolite mappings dictionary.
	- 150723_iMM904_NADcorrected_1127_FTHFLm_nucleotide_conversions.pkl is an example of a nucleotide conversions dictionary.
	- iMM904_NADcorrected_1127_FTHFLm_metabolite_dict.csv contains the format needed to map the number of carbons to every metabolite and create a carbon balanced model, eliminating nonbalanced reactions if neccessary. 
		*To note, eliminating rxns that are not carbon balanced may result in the model not being solvable for a biomass flux if the metabolites involved are the only way to produce certain metabolites in the biomass reaction .


Conversion_Script/
	-This contains the files necessary to convert .xml or .mat metabolic reonstructions into the necessary .pkl file format for EXAMO.
	-To create metabolite mappings or nucleotide conversions (which are options for the conversion script), the dictionaries must first be created and saved as .pkl files.
	-150722_Recon2.v04_metabolite_mappings.py: produces 150722_Recon2.v04_metabolite_mappings.pkl by:
		python 150722_Recon2.v04_metabolite_mappings.py
	-150722_Recon2.v04_nucleotide_conversions.py: produces 150722_Recon2.v04_nucleoitde_conversions.pkl by:
		python 150722_Recon2.v04_nucleotide_conversions.py
	-150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.py: produces 150723_iMM904_NADcorrected_1127_metabolite_mappings.pkl by:
		python 150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.py
	-150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.py: produces 150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.pkl by:
		python 150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.py
	To note, eliminating rxns that are not carbon balanced may result in the model not being solvable for a biomass flux if the metabolites involved are the only way to produce certain metabolites in the biomass reaction .
	-150722_import_cobra_model_46.py: the file conversion script. 
		-To create the .pkl file for Recon2 with metabolite mappings, nucleotide conversions, and a carbon balanced model:
			python 150723_import_cobra_model_46.py Recon2.v04.mat R_biomass_reaction _e_ -m 150722_Recon2.v04_metabolite_mappings.pkl -n 150722_Recon2.v04_nucleotide_conversions.pkl -c recon2_metabolite_dict4.csv
		-To create the .pkl file for iMM904 with metabolite mappings, nucleotide conversions, lower boundary adjustments, upper boundary adjustments, and changes to the gene2rxn mappings:
			python 150722_import_cobra_model_46.py iMM904_NADcorrected_1127_MTHFDi.xml R_biomass_published _e -l iMM904_NADcorrected_1127_MTHFDi_lb_AA_AC_B_Polyamine_Sterol.txt -u iMM904_NADcorrected_1127_MTHFDi_ub.txt -g iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.txt -m 150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.pkl -n 150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.pkl > 150723_iMM904_NADcorrected_1127_MTHFDi_rxns_l_u_g_m_n.txt
		-To create the .pkl file for iAF1260
			python 150722_import_cobra_model_46.py iAF1260.xml R_Ec_biomass_iAF1260_core_59p81M _e
			

150723_recon2_iMM904_iAF1260/
	-This contains updates to the EXAMO package allowing for the import of other metabolic reconstructions besides yeast. 
	-.pkl model files are placed in /data with a geneRules.csv file and the original .xml or .mat files. ##To note, iMM904_NADcorrected.xml was used instead of iMM904_NADcorrected_1127_MTHFDi.xml because the latter did not have default boundary constraints and produced too large of a biomass flux.
	-Flux states are written out as metabolicState*.csv files.
	-To run the first script for Recon2 for TM and TP SKCM:
		python 150720_01_findZeroAndHighFreqRxns_1.py Recon2.v04.mat 20150721210414_m_n_c_Recon2.v04.pkl 150708_SKCM_TP_20 R_biomass_reaction
		python 150720_01_findZeroAndHighFreqRxns_1.py Recon2.v04.mat 20150721210414_m_n_c_Recon2.v04.pkl 150708_SKCM_TM_20 R_biomass_reaction
	-To run the second-fourth scripts for Recon2 with 1 repetition:
		python 150723_02_04.py Recon2.v04.mat 20150721210414_m_n_c_Recon2.v04.pkl 150708_SKCM_TM_20 R_biomass_reaction 1
		python 150723_02_04.py Recon2.v04.mat 20150721210414_m_n_c_Recon2.v04.pkl 150708_SKCM_TP_20 R_biomass_reaction 1
	-To run the first script for iMM904:
		python 150720_01_findZeroAndHighFreqRxns_1.py iMM904_NADcorrected.xml 20150722132206_l_u_g_m_n_iMM904_NADcorrected_1127_MTHFDi.pkl 141013_GCR1A R_biomass_published
	-To run the second-fourth scripts for iMM904 wiht 1 repetition:
		python 150723_02_04.py iMM904_NADcorrected.xml 20150722132206_l_u_g_m_n_iMM904_NADcorrected_1127_MTHFDi.pkl 141013_GCR1A R_biomass_published 1
	-To run the first script for iAF1260:
		python 150720_01_findZeroAndHighFreqRxns_1.py iAF1260.xml 20150722125653_iAF1260.pkl 150616_iAF1260 R_Ec_biomass_iAF1260_core_59p81M
	-To run the second-fourth scripts for iAF1260 with 1 reptition:
		python 150723_02_04.py iAF1260.xml 20150722125653_iAF1260.pkl 150616_iAF1260 R_Ec_biomass_iAF1260_core_59p81M 1


150203_GCR1_Baker_1_lb_AA_AC_B_Polyamine_Sterol_Sugars_para_pass_noexport/
	-This contains the most recent work for the paraellilzation in the mba.py script in /examoModules. Other changes to mba.py have been made in 150723_recon2_iMM904_iAF1260/examoModules/mba.py to allow for other models to be used besides yeast. This parallelization version wasn't used to test other models since it ran approximately 7 times slower than the forking version. This is what I had planned to work on next.   
