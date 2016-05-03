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
from cobra.io import read_sbml_model
from cobra.core import ArrayBasedModel, Model
from cobra.io.mat import _cell
from decimal import Decimal, getcontext, ROUND_DOWN
import cobra
import argparse
import numpy as np
import scipy.io
from mlabwrap import mlab
import os

################################################################################
'''Run methods for condition-specific modeling in MATLAB'''

################################################################################
#Define argparse command line inputs
parser = argparse.ArgumentParser(description='_01_findZeroAndHighFreqRxns.py')
parser.add_argument('m', nargs="+", type=str, help='Necessary variable: MATLAB model file')
parser.add_argument('d', nargs="+", type=str, help='Necessary variable: condition description')
parser.add_argument('-b', nargs="+", type=str, help='Necessary variable: biomass reaction')
parser.add_argument('-o', type=float, help='Defined biomass production')

args = parser.parse_args()

model = str(args.m)
model = model[2:-2]
description = str(args.d)
description = description[2:-2]
if args.b is not None:
    biomass_reaction = str(args.b)
if args.o is not None:
    lb_biomass = args.o

########################################
#Classify reactions by expression
## Importing gene calls as a dictionary
fGeneCalls = 'data/geneRules/geneCalls_%s.csv' % description#genes classified by expression
geneCalls = {}
f = open(fGeneCalls, 'rU')
for line in csv.reader(f):
    geneCalls[line[0]] = int(line[1])
f.close()

data_list = []
genes_cobra_list = []

for i,j in geneCalls.items():
	if j == 1: 
		genes_cobra_list.append(i)
		data_list.append(j)
	if j == -1:
		genes_cobra_list.append(i)
		data_list.append(0)

#Create the expressionData object from the gene rules

data_matlab = np.zeros((len(data_list),), dtype=np.double)
data_matlab[:] = data_list
data_matlab = data_matlab[np.newaxis].T

genes_matlab = np.zeros((len(genes_cobra_list),), dtype=np.object)
genes_matlab[:] = genes_cobra_list
genes_matlab = genes_matlab[np.newaxis].T


expressionData = str('expressionData_' + description)
expressionData = expressionData.replace(".", "")

expressionData_object = {'Data': data_matlab, 'Locus': genes_matlab}

scipy.io.savemat('data/geneRules/%s.mat' % description, expressionData_object)

#Run pFBA, GIMME and iMAT
mlab.addpath(os.getcwd())
try:
	mlab.Cobra_Methods_test_2(model,description,lb_biomass,biomass_reaction)
except NameError:
	print "hi"
	mlab.Cobra_Methods_test_2(model_name,description)
