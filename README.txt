INSTALLATION INSTRUCTIONS
To install LibSBML for Ubuntu, enter from the terminal: (pip install python-libsbml)
	sudo apt-get install python-dev
	sudo apt-get install libxml2-dev
	sudo apt-get install libz-dev
	sudo apt-get install libbz2-dev
	sudo apt-get install python-pip
	sudo pip install python-libsbml

To install cobrapy for Ubuntu, enter from the terminal:
	sudo pip install cobra

	
Install the Gurobi MILP solver following the protocol on Gurobi's website: 
	http://user.gurobi.com/download/gurobi-optimizer
In brief, change into the Gurobi directory with the setup.py file and type:
	python setup.py install	
After installing Gurobi, set Gurobi's path accordingly as descibed in Gurobi's installation protocol. Here is an example of an appropriate file path:
	export GUROBI_HOME="/home/eddie/gurobi563/linux64"
	export PATH="${PATH}:${GUROBI_HOME}/bin"
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
Then go back to gurobi.com and download a free academic licencse. Enter your specific grbgetkey command for your license to activate Gurobi.



FILE DESCRIPTION
The following are the description of the files:
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
