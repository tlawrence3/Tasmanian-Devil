# INTRODUCTION

  We have developed The Algorithm for Simplified Metabolic ANalysIs by Altering Networks and Deducing flux Estimates for VIsuaLization (TASMANIAN DEVIL). There are four separate modules that comprise this software package: gene activity determination, genome-scale metabolic model importation and simplification to reduce network complexity, robust heuristic model building and metabolic flux prediction, and flux visualization from a reference network topology. All modules can be utilized independently or in conjunction with one another. TASMANIAN DEVIL has the potential to be adopted by a wide spectrum of researchers using Linux or macOS platforms. More information for the individual modules can be found at https://tlawrence3.github.io/Tasmanian-Devil/build/html/index.html.

# INSTALLATION

This installation will only work with Python 2.7 or Python 3.6. We highly recommend using a conda environment, as this will manage the pip versions better than a virtual environment when running setup.py.

To install COBRA, need to do the following first:
* For Linux:
	* sudo apt-get install swig
	* For pythhon 2.7:
		* sudo apt-get install libglpk-dev
		* sudo apt-get install glpk-utils
* For macOS:
	* Need to install brew:
		* https://docs.brew.sh/Installation
	* brew install glpk

Downlaod anaconda and build a conda environment: 
* First install Anaconda with Python 2.7 or Python 3.6:
	* https://www.anaconda.com/download/
* Create a conda environment and activate it:
	* https://conda.io/docs/user-guide/tasks/manage-environments.html	
* Then install Gurobi using conda from within the conda environment:	
	* conda install gurobipy
* Then install a Gurobi academic key from their website usirng grbketkey
* conda uninstall scipy
	* We are currently relying on an older version of scipy installed in the setup.py script 
	
Install TASMANIAN-DEVIL:
* git clone https://github.com/tlawrence3/Tasmanian-Devil
* Change into the TASMANIAN-DEVIL folder
	* python setup.py install

# TUTORIAL

Each module's help page contains documentation about usage. The following are test examples and general advice for using each module.

**model module**:

Simplest command for yeast metabolism:
```bash
tas model test_data/iMM904_NADcorrected_1127_MTHFDi.xml _e \
-s > test_iMM904_1.txt
```
Adjusting for lower boundaries, gene2rxn mappings, and 1st attempt at adding adaptations to test fixing the biomass reaction for yeast metabolism:
```bash
tas model test_data/iMM904_NADcorrected_1127_MTHFDi.xml _e \
-l test_data/YPD_lb.csv \
-g test_data/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv \
-a test_data/iMM904_NADcorrected_1127_MTHFDi_adaptations_V1.csv \
-d test_data/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv \
-s > test_iMM904_2.txt
```
After inspecting output for which reactions are unbalanced, try a 2nd adaptation to correct for these reactions. <br />
Adjusting for lower boundaries, gene2rxn mapings, and 2nd attempt at adding adaptations to test fixing unbalanced reactions for yeast metabolism:
```bash
tas model test_data/iMM904_NADcorrected_1127_MTHFDi.xml _e \
-l test_data/YPD_lb.csv \
-g test_data/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv \
-a test_data/iMM904_NADcorrected_1127_MTHFDi_adaptations_V2.csv \
-d test_data/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv \
-s > test_iMM904_3.txt
```
After inspecting output for which reactions are unbalanced still after the adapations, force the model to be carbon balanced. The algorithm will identify any metabolites to be removed from the biomass reaction as needed.

Adjusting for lower boundaries, gene2rxn mapings, adding adaptations to fix unbalanced reactions, accounting for nucleoside conversions and metabolite mapping complexes, removing metabolites without carbons, and remove inactive reactions and carbon balance the model for yeast metabolism:
```bash
tas model test_data/iMM904_NADcorrected_1127_MTHFDi.xml _e \ 
-l test_data/YPD_lb.csv \ 
-g test_data/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv \ 
-a test_data/iMM904_NADcorrected_1127_MTHFDi_adaptations_V2.csv \ 
-m test_data/150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.csv \ 
-n test_data/150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.csv \ 
-d test_data/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv \
-s -z -r > test_iMM904_4.txt
```
Simplest command for human metabolism:
```bash
tas model test_data/Recon2_v04.xml _e -s > test_Recon2_1.txt
```
Model for human metabolism accounting for nucleoside conversions and metabolite mapping complexes, removing metabolites without carbons, and removing inactive reactions and carbon balancing the model:
```bash
tas model test_data/Recon2_v04.xml _e \ 
-m test_data/150722_Recon2.v04_metabolite_mappings.csv \ 
-n test_data/150721_Recon2.v04_nucleotide_conversions.csv \ 
-d test_data/Recon2_metabolite_carbon_dict4.csv -s -z -r > test_Recon2_2.txt
```
**gene module**:

```bash
tas gene test_data/151012_Gasch_glucose.txt \ 
-m test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
-o test_data/151012_glucose_0.25.csv -c
```
**flux module**:

simplest command:
```bash
tas flux test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
test_data/151012_ethanol_0.25.csv _e 1 1 1 -c
```
For multiple concurrent processes, repetitions of concurrent processes and repetitions with pruning reactions in a desired order by compartment:
```bash
tas flux test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
test_data/151012_ethanol_0.25.csv _e 2 2 2 -c -b 0.2879 \ 
-EXrxns test_data/EXrxns.csv \ 
-EXtrrxns test_data/EXtrrxns.csv \ 
-Othertrrxns test_data/Othertrrxns.csv
```
**visualization module**:

Visualizing 1 condition
```bash
tas visualization test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
test_data/151012_glucose_0.25.csv \ 
test_data/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened \ 
1 _e -c \ 
-c1 test_data/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl \
-b1 test_data/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl
```
Comparing 2 conditions
```bash
tas visualization test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
test_data/151012_glucose_0.25.csv \ 
test_data/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened \ 
1 _e -c \
-c1 test_data/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl \
-b1 test_data/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl \ 
-c2 test_data/RxnsClassifiedByExpression_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl \
-b2 test_data/freqBasedRxns_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl \ 
-m2 test_data/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat \ 
-g2 test_data/151012_ethanol_0.25.csv \ 
-f2 test_data/metabolicState_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi
```
