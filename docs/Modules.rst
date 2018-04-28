Modules
=======

A more detailed description of optional longer tags and variable inputs are available on the project's help pages for each module that comes with the installation.

model
~~~~~


Make modifications to SBML or COBRA models. Look at iMM904 example in test_data folder of installation for examples of formatting files

-model		Necessary variable: metabolic reconstruction file
-extracellular	Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e')
-o		Name of model to be written out	
-c		Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag
-s		Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag
-l		lower boundary constraints file
-u		upper boundary constraints file
-g		gene2rxn file
-d		Tab-delimited file to specify dicitonary mappings of number of carbons in every metabolite. This is to check whether the model is carbon balanced. See iMM904 example for documentation
-m		Tab-delimited metabolite mapping complexes file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this if metFormulas is not in model
-n		Tab-delimited nucleotide conversions file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this if metFormulas is not in model
-a		Tab delimited file to change specific reaction stoichiometries or to remove metabolites from reactions. See iMM904 example for documentation; must have -d argument as well to use this if metFormulas is not in model
-z		Flag to specify whether to remove metabolites wihtout any carbons. Must have -d argument as well to use this if metFormulas is not in model
-r		Flag to specify whether to remove inactive reactions from final model and remove carbon unbalanced reactions. Recommend using only after first inspecting reactions to be removed. Must have -d argument as well to use this if metFormulas is not in model

gene
~~~~

Classify gene activity in the context of a model

-m	Metabolic reconstruction file
-c	Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag
-s	Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag
-o	gene_classification.csv", help='Name of csv gene rule file to be written out
-u	upper threshold by which to define genes as being active; default value: upper 25% of genes mapped to reactions
-l	lower threshold by which to define genes as being inactive; default value: lower 25% of genes mapped to reactions


flux
~~~~

Predict condition-specific fluxes

-model					Necessary variable: metabolic reconstruction file
-description				Necessary variable: comma separated gene rule file
-extracellular				Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e')
-concurrentprocesses			Necessary variable: number of concurrent processes for creating pruned models
-repetitionsofconcurrentprocesses	Necessary variable: number of times to repeat the chosen number of processes to create the profile of models
-repetitionsoffluxstates		Necessary variable: number of repetitions to produce final flux states (csv files)
-c					Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag
-s					Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag
-e					Minimum value for fluxes forced to be non-zero (for the implementation of the pathway algorithm created by Shlomi); default value: 1E-1
-z					Activity threshold when classifying HFR and ZFR reactions from gene rules 
-i					Minimum value for active flux of reversible reactions in an irreversible model; default value: 1E-10
-t					Activity threshold for finding active reactions when pruning a model
-a					Activity threshold (above which a flux is considered to be larger than 0 when minimizing the sum of fluxes); default value: 1E-10
-b					Defined biomass production 
-EXrxns					Comma separated file containing extracellular reactions
-EXtrrxns				Comma separated file containing extracellular transport reactions
-Othertrrxns				Comma separated file containing other compartmental transport reactions

visualization
~~~~~~~~~~~~~

Visualize fluxes on predefined maps

-model1				Necessary variable: metabolic reconstruction file 1
-geneCalls1			Necessary variable: comma separated gene rule file 1
-fluxState1			Necessary variable: Name of flux state 1
-pathways			Necessary variable: name of pathway(s)
-repetitionsoffluxstates	Necessary variable: number of repetitions to produce final flux states (csv files)
-extracellular			Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e')
-c				Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag
-s				Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag
--rxnsClassifiedByExpression1	Reactions classified by expression pickle file 1 from flux module
--freqBasedRxns1		Reactions classified by frequency pickle file 1 from flux module
--model2			Metabolic reconstruction file 2
--geneCalls2			Comma separated gene rule file 2
--fluxState2			Name of flux state 2
--rxnsClassifiedByExpression2	Reactions classified by expression pickle file 2 from flux module
--freqBasedRxns2		Reactions classified by frequency pickle file 2 from flux module


When creating the predefined pathway map in CellDesigner first, make sure to select "Grid Snap" and "Grid Visible" under the "Edit" tab. Space all reactants at least 12 cells away from products for visualization purposes. Format metabolites as reactants or products based on their predefined state in the starting SBML model. This will ensure that reversibility is properly taken into account for a negative flux. If a reaction has 1 reactant and 2 products or 2 products and 1 reactant, make it a "Dissociation" or a "Heterodimer Association", respectively. Otherwise, make all reactions a "State Transition". If there are more than 1 product or 1 reactant that doesn't match the cases described above, add them as reactants or products to the reaction by opting for "Add Reactant" and "Add Product" with the reaction selected. Name the metabolites as you create them. Do not name the reactions as you are creating them. Rather, replace the reaction IDs by "Replace Reaction ID" under the "Edit" tab. Name the reactions as they were in the model, except without the "R\_" appended to the front for each reaction name. Save the model as an SBML file.
