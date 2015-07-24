# Copyright (C) 2012 Sergio Rossell
#
# This script is part of the EXAMO software
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
# 



#This package was modified by Eddie Gibb (EG). Further documentation of changes can be found in the paper.


"""
This script accepts a comma separated value file with gene expression calls and a genome-scale constraint based metabolic model as a dictionary with the following content:
 
'S' : [scipy.sparse coo_matrix] stoichiometry matrix
'idSp' : [list of strings] list of ids of the model's metabolic species (metabolites)
'idRs' : [list of strings] list of ids of the model's reactions
'lb' : [list of floats] lowwer flux boundaries for the model's reactions
'ub' : [list of floats] upper flux boundaries for the model's reactions
'gene2rxn' : [dict] keys are reaction ids, values are the reactions Boolean gene-to-reaction mapping
'genes' : [list of strings] list of all genes in the model 

The script's outputs are writen in the `data' directory. 

OUTPUTS

1. rxnsClassifiedByExpression_description.pkl [dict]. The keys of this dictionary are 'rH' (highly expressed reactions), 'rL' (lowly expressed reactions) and 'rU' (all other reactions --unclassified). The values of each key are lists of the reactions classified into each of these categories.

2. freqBasedRxns_description.pkl [dict]. Keys are 'hfr' (high-frequency-reactions) and 'zfr' (zero-frequency-reactions). THe values of each key are lists of the reactions classified into each of these categories.
"""

import csv
from sys import path; path.append('./examoModules/')

from examoModules import *

################################################################################
# _01_findZeroAndHighFreqRxns.py

################################################################################
# INPUTS
description = sys.argv[1]


fGeneCalls = 'data/geneCalls_%s.csv' % description#genes classified by expression

# model dictionary of the original model with blocked reactions deleted
fModelDict = 'data/150126_degapped_iMM904_with_metabolite_mapping_complexes_36_stoichiometry_lb_AA_AC_B_Polyamine_Sterol.pkl'

# id of the biomass reaction
biomassRxn = 'R_biomass_published'

# Algorithm variables
#threshold above which a flux is considered to be larger than zero
activityThreshold = 1E-10

# minimum value for fluxes forced to be non-zero
eps = 1E-1 #EG Using a higher threshold for the 1st script

# files for export
fOutRxnsByExpression = 'data/rxnsClassifiedByExprssion_%s.pkl'

fOutFreqBasedRxns = 'data/freqBasedRxns_%s.pkl'

################################################################################
# STATEMENTS

########################################
# 0. Importing model dictionary
md = importPickle(fModelDict)
m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])

#copy of rxns identifiers
idRs = m.idRs[:]
idRs.remove(biomassRxn)
#forcing non-zero biomass production
m.lb[m.idRs.index(biomassRxn)] = 0.2879

########################################
# 1. Classify reactions by expression
## Importing gene calls as a dictionary
geneCalls = {}
f = open(fGeneCalls, 'rU')
for line in csv.reader(f):
    geneCalls[line[0]] = int(line[1])
f.close()
## Classifying reactions by expression
rxnDict = classifyRxnsByExpression(geneCalls, m.gene2rxn, m.genes)
exportPickle(rxnDict, fOutRxnsByExpression % description)

# 2. Maximizing agreement score and exploring alternative optima
rH = rxnDict['rH']
rL = rxnDict['rL']
#EG added eps to arguments
model = MetabGeneExpModel_gurobi(m.idSp, m.idRs, m.S, m.lb, m.ub, rH, rL, eps)
scores, imgeSols = model.exploreAlternativeOptima(idRs)

# 3. identifying zero and high frequency reactions
zfr, hfr = getZeroAndHighFrequencyRxns(scores, imgeSols, idRs)
exportPickle({'zfr' : zfr, 'hfr' : hfr}, fOutFreqBasedRxns % description)

print 'zfr ', len(zfr)
print 'hfr ', len(hfr)
