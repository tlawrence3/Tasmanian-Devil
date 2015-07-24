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
"""
Classifying reactions by expression.

Genes are classified by expression according to a high and a low threshold. Genes with intensities bellow the threshold are assigned a -1 call, genes with intesities above the high threshold are assigned a 1 call, all other genes ar assigned a 0 call.
"""

import re

################################################################################
# FUNCTIONS

########################################
# CLASSIFYING REACTIONS BY EXPRESSION

def orderGeneNamesByLength(geneList):
    """
    Group gene names by their lengths
    RETURNS byLengthDict [dict]: {stringLength : [list of gene ids that are
        that long]
    """
    byLengthDict = {}
    lengthSet = set()
    for gn in geneList:
        lengthSet.add(len(gn))
    for length in lengthSet:
        for gn in geneList:
            if length == len(gn):
                try:
                    byLengthDict[length].append(gn)
                except KeyError:
                    byLengthDict[length] = [gn]
    return byLengthDict

def substituteGeneNamesByCalls(g2r, geneList, geneCalls, callTranslationDict):
    """
    ACCEPTS
    g2r [str] Boolean gene-to-reaction mapping
    geneLists [list] list of genes who's calls will be subsitituted in g2r
    geneCalls [dict] {geneId : geneCall}
    callTranslationDict [dict] : {geneCall : rxnCall} rules for translating a
        geneCall to a rxnCall (see createRxnGeneCalls)
    RETURNS:
    gs [str] string with gene names substituted by their expression calls
    """    
    for gn in geneList:
        try:
            if geneCalls[gn] == -1:#lowly expressed
                g2r = g2r.replace(gn, callTranslationDict[-1])
            elif geneCalls[gn] == 0:#undecided expression
                g2r = g2r.replace(gn, callTranslationDict[0])
            elif geneCalls[gn] == 1:#highly expressed
                g2r = g2r.replace(gn, callTranslationDict[1])
        except KeyError:# gene expression unavailable
            g2r = g2r.replace(gn, callTranslationDict[0])#undecided expression
    return g2r

def createRxnGeneCalls(geneCalls, gene2rxn, modelGenes):
    """
    creates a dictionary of reaction calls {rxn : call}, where call may be
    -1 (lowly expressed), 0 (moderate expression) or 1 (highly expressed).
    ACCEPTS
    geneCalls [dict]: {geneName : call} where calls are as explained above
    gene2rxn [dict]: {rxn : Boolean gene to rxn mapping (string)} 
    modelGenes [set]: gene names included in the model
    RETURNS
    rxnGeneCalls [dict] : {rxn : calls}, where calls may be -1, 0 or 1
    """
    orphanRxns = set([rxn for rxn in gene2rxn if gene2rxn[rxn] == ''])
    #ordering gene names by their lenghts (some gene names are part of others
    # e.g. 'YCR024C' and 'YCR024C-A' in yeast). Hence long names should be 
    # replaced before short ones
    geneNamesByLength = orderGeneNamesByLength(modelGenes)
    orderedGeneLengths = geneNamesByLength.keys()
    orderedGeneLengths.sort(reverse = True)
    rxnGeneCalls = {}
    for rxn in set(gene2rxn.keys()) - orphanRxns:
        # lowly expressed genes are evaluated to -1 while higly expressed
        # genes are evaluated to 1. Lowly expressed reactions are those
        # for which their Boolean expression evaluates to -1
        g2r = gene2rxn[rxn]
        for length in orderedGeneLengths:
            g2r = substituteGeneNamesByCalls(g2r, geneNamesByLength[length], 
                    geneCalls, {-1 : '-1', 0 : '0', 1 : '1'})
	g2r = str(g2r)
	i = re.split('or',g2r)
	count = 0
	countdict = {}
	countdict_calls = {}
	for j in i:
		count += 1
		j = re.split('and', j)
		genedict = []
		for k in j:
			k = re.sub(' ', '', k)
			k = re.sub('\(', '', k)
			k = re.sub('\)', '', k)
			genedict.append(k)
		countdict[count] = genedict
	for i in countdict:
		one_count = 0
		negative_one_count = 0
		for j in countdict[i]:
			if int(j) == 1:
				one_count += 1
			if int(j) == -1:
				negative_one_count += 1
		if len(countdict[i]) == one_count:
			countdict_calls[i] = 1
		elif negative_one_count > 0:
			countdict_calls[i] = -1
		else:
			countdict_calls[i] = 0
	call_one = 0
	call_zero = 0
	call_negative_one = 0
	for i in countdict_calls:
		if countdict_calls[i] == 1:
			call_one += 1
		elif countdict_calls[i] == 0:
			call_zero += 1
		else:
			call_negative_one += 1
	if call_one >= 1:
		rxnGeneCalls[rxn] = 1
	elif call_zero >= 1:	
		rxnGeneCalls[rxn] = 0
	else:
		rxnGeneCalls[rxn] = -1 	
    return rxnGeneCalls

def classifyRxnsByExpression(geneCalls, gene2rxn, modelGenes):
    """
    Classifies reactions with gene associations as  highly, uncertain or 
    lowly expressed. 
    rxnDict [dict]: {rxn : rxnObject} cbModel.rxns or cbModel.reactions (old)
    geneCalls [dict]: {geneName : call} where calls are as explained above
    modelGenes [set]: gene names included in the model
    RETURNS
    d [dict]: {'rL' : lowlyExpressed, 'rU' : uncertainExpression, 
        'rH' : highlyExpressed}
    """
    rxnCalls = createRxnGeneCalls(geneCalls, gene2rxn, modelGenes)
    lowlyExpressed = set()
    highlyExpressed = set()
    uncertainExpression = set()
    for rxn in rxnCalls:
        if rxnCalls[rxn] == -1:
            lowlyExpressed.add(rxn)
        elif rxnCalls[rxn] == 0:
            uncertainExpression.add(rxn)
        elif rxnCalls[rxn] == 1:
            highlyExpressed.add(rxn)
    d = {'rL' : lowlyExpressed, 'rU' : uncertainExpression, 
            'rH' : highlyExpressed}
    return d



################################################################################
# TESTING 

if __name__ == '__main__':       
    import csv
    ########################################
    # INPUTS
    fGeneCalls = 'data/geneCalls/geneCalls_eth_15_85.csv'
    fModelDict = 'data/iMM904_blkRxnsDeleted_dict.pkl'
    
    ## Importing gene-to-reaction mapping
    gene2rxn = importPickle(fModelDict)['gene2rxn']
    modelGenes = importPickle(fModelDict)['genes']
    ## Importing gene calls as a dictionary
    geneCalls = {}
    f = open(fGeneCalls, 'rU')
    for line in csv.reader(f):
        geneCalls[line[0]] = int(line[1])
    f.close()
    #NOTE createRxnGeneCalls might be the best way to report rH, rU and RL
    #rxnCalls = createRxnGeneCalls(geneCalls, gene2rxn, modelGenes)
    
    rxnsByExpression = classifyRxnsByExpression(geneCalls, gene2rxn, modelGenes)
    exportPickle(rxnsByExpression, 'data/rxnDictTemp.pkl')

