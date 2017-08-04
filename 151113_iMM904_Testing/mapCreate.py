#This script was created by Eddie Gibb. 
#Script to visualize pathways. Please read the README file for details on creating pathway maps in CellDesigner.

import argparse
import csv
import numpy as np
from examoModules import *

# INPUTS
#Define argparse command line inputs
parser = argparse.ArgumentParser(description='mapCreate.py')
parser.add_argument('m1', nargs="+", type=str, help='Necessary variable: input model 1 file')
parser.add_argument('g1', nargs="+", type=str, help='Necessary variable: gene calls 1 file')
parser.add_argument('f1', nargs="+", type=str, help='Necessary variable: name of flux state 1')
parser.add_argument('r', nargs="+", type=int, help='Necessary variable: Number of repetitions run for the 2nd-4th scripts for final flux states')
parser.add_argument('-p', nargs = "+", type=str, help='Necessary variable: Pathway(s)')
parser.add_argument('-m2', nargs="+", type=str, help='Optional variable: input model 2 file')
parser.add_argument('-g2', nargs="+", type=str, help='Optional variable: gene calls 2 file')
parser.add_argument('-f2', nargs="+", type=str, help='Optional variable: name of flux state 2')


args = parser.parse_args()


model1 = str(args.m1)
model1 = model1[2:-2]
pickle_model1 = importPickle('data/%s' % model1)
md_model1 = CbModel(pickle_model1['S'], pickle_model1['idSp'], pickle_model1['idRs'], pickle_model1['lb'], pickle_model1['ub'], pickle_model1['rxns'], pickle_model1['genes'])


geneCallsFile1 = str(args.g1)
geneCallsFile1 = geneCallsFile1[2:-2]
geneCallsFile1 = open('data/%s' % geneCallsFile1, 'r')
csvreader1 = csv.reader(geneCallsFile1)
#Import the gene rules and create dictionaries
geneCalls1 = {}
for row in csvreader1:
	geneCalls1[row[0]] = int(row[1])
geneCallsFile1.close()



repetitions = str(args.r)
repetitions = int(repetitions[1:-1])



fluxstate1 = str(args.f1)
fluxstate1 = fluxstate1[2:-2]
fModelDict1_hfr = 'data/freqBasedRxns_%s.pkl' % fluxstate1
fbr1 = importPickle(fModelDict1_hfr)
fOutRxnsByExpression1 = 'data/rxnsClassifiedByExprssion_%s.pkl' % fluxstate1
gbr1 = importPickle(fOutRxnsByExpression1) 
#Classify the rxns as being in rH, rL or neither and hfr, zfr, or neither for both conditions
gbr1_rH = {}
fbr1_hfr = {}
for t in md_model1.rxns.keys():
	if t in gbr1['rH']:
		gbr1_rH[t] = 2
	elif t in gbr1['rL']:
		gbr1_rH[t] = 0
	else:
		gbr1_rH[t] = 1
	if t in fbr1['hfr']:
		fbr1_hfr[t] = 2
	elif t in fbr1['zfr']:
		fbr1_hfr[t] = 0
	else:
		fbr1_hfr[t] = 1
rxnFreq = {}
#Map raw data for every gene for every reaction 
rxngenesandor1 = {}
for t in md_model1.rxns.keys():
	genelist = str(md_model1.gene2rxn[t])
	i = re.split('or',genelist)
	count = 0
	or_count_1 = 0
	countdict_rules = {}
	countdict_calls = {}
	rxn_genes_and_or = [t]
	for j in i:
		count += 1
		j = re.split('and', j)
		genedict_rules = {}
		gene_list = []
		for k in j:
			k = re.sub(' ', '', k)
			k = re.sub('\(', '', k)
			k = re.sub('\)', '', k)
			if k != "":
				if k in geneCalls1:
					gene_list.append(k)
					countdict_rules[count] = gene_list 				
	for l in countdict_rules:
		one_count = 0
		negative_one_count = 0
		for m in countdict_rules[l]:
			if geneCalls1[m] == 1:
				one_count += 1
			if geneCalls1[m] == -1:
				negative_one_count += 1
		#Make it so that rules are 2, 1, and 0 instead of 1, 0, and -1 respectively, so that the descriptions can be directly given for CellDesigner
		if len(countdict_rules[l]) == one_count:
			countdict_calls[l] = 2
		elif negative_one_count > 0:
			countdict_calls[l] = 0
		else:
			countdict_calls[l] = 1
		if (len(countdict_rules[l]) == 1 and len(countdict_rules) > 1):
			or_count_1 += 1
		else:
			if len(countdict_rules[l]) == 1:
				rxn_genes_and_or.extend(["_",str(len(countdict_rules[l])),"ONLY",str(countdict_calls[l])])
			else:
				rxn_genes_and_or.extend(["_",str(len(countdict_rules[l])),"AND",str(countdict_calls[l])])
	if (or_count_1 > 0):
		rxn_genes_and_or.extend(["_",str(or_count_1),"OR",str(gbr1_rH[t])])
	rxn_genes_and_or_join = ''.join(rxn_genes_and_or)
	rxngenesandor1[t] = rxn_genes_and_or_join
# Retrieving MBA candidate reaction lists
for i in range(0,repetitions):
	mbaCandRxnsDirectory = 'data/mbaCandRxns/%s_%d/' % (fluxstate1, i)
	files = os.popen('ls %s | grep %s' % (mbaCandRxnsDirectory, fluxstate1)).read().splitlines()
	rxnSets = []
	for fn in files:
		rxnSets.append(importPickle(mbaCandRxnsDirectory + fn))
	# Quantifying the number of times a rxn is among the candidate models
	rxnFreq[i] = {}
	for rs in rxnSets:
		for rxn in rs:
			try:
				rxnFreq[i][rxn] += 1
			except KeyError:
				rxnFreq[i][rxn] = 1
	for t in md_model1.rxns.keys():
		if t not in rxnFreq[i]:
			rxnFreq[i][t] = 0
cond1freqavg = {}
for t in md_model1.rxns.keys():
	freqavg_list = []
	for i in range(0,repetitions):
		freqavg_list.append(float(rxnFreq[i][t]))
	freqavg = np.mean(freqavg_list)
	freqavg = str("%.2f" % freqavg)
	if freqavg == "0.00":
		freqtext = [t]
		freqtext.append("_")
		freqtext.append("0")
		freqtext.append("_")
		freqtext.append("0")
		freqtextjoin = ''.join(freqtext)
		cond1freqavg[t] = freqtextjoin[2:]
	else:
		frequencies = []
		for i in range(0,repetitions):
			frequencies.append(rxnFreq[i][t])
		if np.std(frequencies) != 0:
			freqavgsearch = re.search('\.00', freqavg)
			if (freqavgsearch is not None):
				freqavg = freqavg[:-3]
			else:
				freqavg = re.sub('\.','',freqavg)
				freqavg+='E_2'
			freqtext = [t]
			freqtext.append("_")
			freqtext.append(str(freqavg))
			freqtext.append("_")
			if np.std(frequencies) == 0:
				freqtext.append("0")
			else:
				std = round(np.std(frequencies),2)
				std = str("%.2f" % std)
				stdsearch = re.search('\.00', std)
				stdsearch0 = re.search('0\.', std)
				if (stdsearch is not None):
					std = std[:-3]
				elif (stdsearch0 is not None):
					std = std[2:]
					std+='E_2'
				else:
					std = re.sub('\.','',std)
					std+='E_2'
				freqtext.append(std)
		else:
			freqtext = [t]
			freqtext.append("_")
			freqtext.append(str(freqavg[:-3]))
			freqtext.append("_0")
		freqtextjoin = ''.join(freqtext)
		cond1freqavg[t] = freqtextjoin[2:]
fluxdict = {}
fluxavgdict1 = {}
fluxvardict1 = {}
fluxstddict1 = {}
idRsflux = {}
notincond = {}
for i in range(0,repetitions):
	metabolicState_file = open('data/metabolicState_%s_%d.csv' % (fluxstate1, i), 'r')
	csvreader = csv.reader(metabolicState_file)
	idRsflux[i] = {}
	idRscsv = []
	flux = []
	for row in csvreader:
		idRscsv.append(row[0])
		flux.append(float(row[1]))
	for j, item1 in enumerate(idRscsv):
		idRsflux[i][item1] = flux[j]
	for t in pickle_model1['rxns']:
		if t not in notincond:
			notincond[t] = 0
		if t not in fluxdict:
			fluxdict[t] = []
		if t not in idRsflux[i]: 
			notincond[t] += 1
			fluxdict[t].append(0)
		else:
			fluxdict[t].append(idRsflux[i][t])
cond1fluxavg = {}
for t in pickle_model1['rxns']:
	if notincond == repetitions: 
		fluxtext = [t]
		fluxtext.append("_NA")
		fluxtextjoin = ''.join(fluxtext)
		cond1fluxavg[t] = fluxtextjoin[2:]
		fluxavgdict1[t] = 0
		fluxvardict1[t] = 0
		fluxstddict1[t] = 0
	else:
		fluxavgdict1[t] = round(np.mean(fluxdict[t]),4)
		fluxvardict1[t] = round(np.var(fluxdict[t]),4)
		fluxstddict1[t] = round(np.std(fluxdict[t]),4)
		fluxavg = fluxavgdict1[t]
		fluxavg = str("%.4f" % fluxavg)
		negative_count = 0
		fluxavg_search = re.search('\-', fluxavg)
		if (fluxavg_search is not None):
			fluxavg = fluxavg[1:]
			negative_count = 1
		if fluxavg == "0.0000":
			fluxtext = [t]
			fluxtext.append("_")
			fluxtext.append("0")
			fluxtext.append("_")
			fluxtext.append("0")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
		else:
			if (len(fluxavg) == 6):
				fluxavgsearch = re.search('0\.', fluxavg)
				if (fluxavgsearch is not None):
					fluxavg = fluxavg[2:]
					fluxavg+='E_4'
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='E_4'
				if negative_count == 1:
					fluxavg = "_" + fluxavg
			else:
				fluxavg = re.sub('\.','',fluxavg)
				fluxavg+='E_4'
				if negative_count == 1:
					fluxavg = "_" + fluxavg
			fluxtext = [t]
			fluxtext.append("_")
			fluxtext.append(str(fluxavg))
			fluxtext.append("_")
			if fluxavgdict1[t] == 0:
				fluxtext.append("0")
			else:
				std = round(np.std(fluxdict[t]),4)
				std = str("%.4f" % std)
				if (len(std) == 6):
					stdsearch = re.search('0\.', std)
					if (stdsearch is not None):
						std = std[2:]
						std+='E_4'
					else:
						std = re.sub('\.','',std)
						std+='E_4'
				else:
					std = re.sub('\.','',std)
					std+='E_4'
				fluxtext.append(std)
			fluxtextjoin = ''.join(fluxtext)
			fluxtextjoin = re.sub('E_4_0000E_4','E_4_0',fluxtextjoin)				
			cond1fluxavg[t] = fluxtextjoin[2:]


pathway_list = []
for pathway in args.p:
	pathway_list.append(pathway)



required_together = ('m2','g2','f2')

if not all([getattr(args,x) for x in required_together]):
    raise RuntimeError("Must supply m2, g2, and f2 if going to use second condition")


if args.m2 is not None: 
	model2 = str(args.m2)
	model2 = model2[2:-2]
	pickle_model2 = importPickle('data/%s' % model2)
	md_model2 = CbModel(pickle_model2['S'], pickle_model2['idSp'], pickle_model2['idRs'], pickle_model2['lb'], pickle_model2['ub'], pickle_model2['rxns'], pickle_model2['genes'])


if args.g2 is not None:
	geneCallsFile2 = str(args.g2)
	geneCallsFile2 = geneCallsFile2[2:-2]
	geneCallsFile2 = open('data/%s' % geneCallsFile2, 'r')
	csvreader2 = csv.reader(geneCallsFile2)
	geneCalls2 = {}
	for row in csvreader2:
		geneCalls2[row[0]] = int(row[1])
	geneCallsFile2.close()
	

#Initialize cond2fluxavg
cond2fluxavg = None
if args.f2 is not None:
	fluxstate2 = str(args.f2)
	fluxstate2 = fluxstate2[2:-2]
	fModelDict2_hfr = 'data/freqBasedRxns_%s.pkl' % fluxstate2
	fbr2 = importPickle(fModelDict2_hfr)
	fOutRxnsByExpression2 = 'data/rxnsClassifiedByExprssion_%s.pkl' % fluxstate2
	gbr2 = importPickle(fOutRxnsByExpression2) 
	#Classify the rxns as being in rH, rL or neither and hfr, zfr, or neither for both conditions
	gbr2_rH = {}
	fbr2_hfr = {}
	for t in md_model2.rxns.keys():
		if t in gbr2['rH']:
			gbr2_rH[t] = 2
		elif t in gbr2['rL']:
			gbr2_rH[t] = 0
		else:
			gbr2_rH[t] = 1
		if t in fbr2['hfr']:
			fbr2_hfr[t] = 2
		elif t in fbr2['zfr']:
			fbr2_hfr[t] = 0
		else:
			fbr2_hfr[t] = 1
	rxnFreq = {}
	#Map raw data for every gene for every reaction 
	rxngenesandor2 = {}
	for t in md_model2.rxns.keys():
		genelist = str(md_model2.gene2rxn[t])
		i = re.split('or',genelist)
		count = 0
		or_count_1 = 0
		countdict_rules = {}
		countdict_calls = {}
		rxn_genes_and_or = [t]
		for j in i:
			count += 1
			j = re.split('and', j)
			genedict_rules = {}
			gene_list = []	
			for k in j:
				k = re.sub(' ', '', k)
				k = re.sub('\(', '', k)
				k = re.sub('\)', '', k)
				if k != "":
					if k in geneCalls2:
						gene_list.append(k)
						countdict_rules[count] = gene_list 				
		for l in countdict_rules:
			one_count = 0
			negative_one_count = 0
			for m in countdict_rules[l]:
				if geneCalls2[m] == 1:
					one_count += 1
				if geneCalls2[m] == -1:
					negative_one_count += 1
			#Make it so that rules are 2, 1, and 0 instead of 1, 0, and -1 respectively, so that the descriptions can be directly given for CellDesigner
			if len(countdict_rules[l]) == one_count:
				countdict_calls[l] = 2
			elif negative_one_count > 0:
				countdict_calls[l] = 0
			else:
				countdict_calls[l] = 1
			if (len(countdict_rules[l]) == 1 and len(countdict_rules) > 1):
				or_count_1 += 1
			else:
				if len(countdict_rules[l]) == 1:
					rxn_genes_and_or.extend(["_",str(len(countdict_rules[l])),"ONLY",str(countdict_calls[l])])
				else:
					rxn_genes_and_or.extend(["_",str(len(countdict_rules[l])),"AND",str(countdict_calls[l])])
		if (or_count_1 > 0):
			rxn_genes_and_or.extend(["_",str(or_count_1),"OR",str(gbr2_rH[t])])
		rxn_genes_and_or_join = ''.join(rxn_genes_and_or)
		rxngenesandor2[t] = rxn_genes_and_or_join
	# Retrieving MBA candidate reaction lists
	for i in range(0,repetitions):
		mbaCandRxnsDirectory = 'data/mbaCandRxns/%s_%d/' % (fluxstate2, i)
		files = os.popen('ls %s | grep %s' % (mbaCandRxnsDirectory, fluxstate2)).read().splitlines()
		rxnSets = []
		for fn in files:
			rxnSets.append(importPickle(mbaCandRxnsDirectory + fn))
		# Quantifying the number of times a rxn is among the candidate models
		rxnFreq[i] = {}
		for rs in rxnSets:
			for rxn in rs:
				try:
					rxnFreq[i][rxn] += 1
				except KeyError:
					rxnFreq[i][rxn] = 1
		for t in md_model2.rxns.keys():
			if t not in rxnFreq[i]:
				rxnFreq[i][t] = 0
	cond2freqavg = {}
	for t in md_model2.rxns.keys():
		freqavg_list = []
		for i in range(0,repetitions):
			freqavg_list.append(float(rxnFreq[i][t]))
		freqavg = np.mean(freqavg_list)
		freqavg = str("%.2f" % freqavg)
		if freqavg == "0.00":
			freqtext = [t]
			freqtext.append("_")
			freqtext.append("0")
			freqtext.append("_")
			freqtext.append("0")
			freqtextjoin = ''.join(freqtext)
			cond2freqavg[t] = freqtextjoin[2:]
		else:
			frequencies = []
			for i in range(0,repetitions):
				frequencies.append(rxnFreq[i][t])
			if np.std(frequencies) != 0:
				freqavgsearch = re.search('\.00', freqavg)
				if (freqavgsearch is not None):
					freqavg = freqavg[:-3]
				else:
					freqavg = re.sub('\.','',freqavg)
					freqavg+='E_2'
				freqtext = [t]
				freqtext.append("_")
				freqtext.append(str(freqavg))
				freqtext.append("_")
				if np.std(frequencies) == 0:
					freqtext.append("0")
				else:
					std = round(np.std(frequencies),2)
					std = str("%.2f" % std)
					stdsearch = re.search('\.00', std)
					stdsearch0 = re.search('0\.', std)
					if (stdsearch is not None):
						std = std[:-3]
					elif (stdsearch0 is not None):
						std = std[2:]
						std+='E_2'
					else:
						std = re.sub('\.','',std)
						std+='E_2'
					freqtext.append(std)
			else:
				freqtext = [t]
				freqtext.append("_")
				freqtext.append(str(freqavg[:-3]))
				freqtext.append("_0")
			freqtextjoin = ''.join(freqtext)
			cond2freqavg[t] = freqtextjoin[2:]
	fluxdict = {}
	fluxavgdict2 = {}
	fluxvardict2 = {}
	fluxstddict2 = {}
	idRsflux = {}
	notincond = {}
	for i in range(0,repetitions):
		metabolicState_file = open('data/metabolicState_%s_%d.csv' % (fluxstate2, i), 'r')
		csvreader = csv.reader(metabolicState_file)
		idRsflux[i] = {}
		idRscsv = []
		flux = []
		for row in csvreader:
			idRscsv.append(row[0])
			flux.append(float(row[1]))
		for j, item1 in enumerate(idRscsv):
			idRsflux[i][item1] = flux[j]
		for t in pickle_model1['rxns']:
			if t not in notincond:
				notincond[t] = 0
			if t not in fluxdict:
				fluxdict[t] = []
			if t not in idRsflux[i]: 
				notincond[t] += 1
				fluxdict[t].append(0)
			else:
				fluxdict[t].append(idRsflux[i][t])
	cond2fluxavg = {}
	for t in pickle_model2['rxns']:
		if notincond == repetitions: 
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond2fluxavg[t] = fluxtextjoin[2:]
			fluxavgdict2[t] = 0
			fluxvardict2[t] = 0
			fluxstddict2[t] = 0
		else:
			fluxavgdict2[t] = round(np.mean(fluxdict[t]),4)
			fluxvardict2[t] = round(np.var(fluxdict[t]),4)
			fluxstddict2[t] = round(np.std(fluxdict[t]),4)
			fluxavg = fluxavgdict2[t]
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond2fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='E_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='E_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='E_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxtext.append("_")
				if fluxavgdict1[t] == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxdict[t]),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='E_4'
						else:
							std = re.sub('\.','',std)
							std+='E_4'
					else:
						std = re.sub('\.','',std)
						std+='E_4'
					fluxtext.append(std)
				fluxtextjoin = ''.join(fluxtext)
				fluxtextjoin = re.sub('E_4_0000E_4','E_4_0',fluxtextjoin)				
				cond2fluxavg[t] = fluxtextjoin[2:]
	


#Scale the fluxes for the pathway maps.
if cond2fluxavg:
	cond1fluxavgscaled = {}
	cond2fluxavgscaled = {}
	max1 = max(fluxavgdict1, key=fluxavgdict1.get)
	max2 = max(fluxavgdict2, key=fluxavgdict2.get)
	if fluxavgdict1[max1] > fluxavgdict2[max2]:
		max_both = fluxavgdict1[max1]
	else:
		max_both = fluxavgdict2[max2] 
	for t in sorted(pickle_model1['rxns']):
		cond1fluxavgscaled[t] = float(round((abs(fluxavgdict1[t])/max_both)*9+1,3))
	for t in sorted(pickle_model2['rxns']):
		cond2fluxavgscaled[t] = float(round((abs(fluxavgdict2[t])/max_both)*9+1,3))
else:
	cond1fluxavgscaled = {}	
	max1 = max(fluxavgdict1, key=fluxavgdict1.get)
	for t in sorted(pickle_model1['rxns']):
		cond1fluxavgscaled[t] = float(round((abs(fluxavgdict1[t])/max1)*9+1,3))



for pathway in pathway_list:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_Rules_%s.xml' % fluxstate1,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
		
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(rxngenesandor1[rxn][2:]),reference[o])
			if (md_model1.lb[md_model1.idRs.index(rxn)] < 0 and md_model1.ub[md_model1.idRs.index(rxn)] > 0):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if gbr1_rH[rxn] == 2:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			if gbr1_rH[rxn] == 0:
				lines2_join = re.sub(color.group(0),'ff0000ff', lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

for pathway in pathway_list:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_Rules_%s.xml' % fluxstate2,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(rxngenesandor2[rxn][2:]),reference[o])
			if (md_model2.lb[md_model2.idRs.index(rxn)] < 0 and md_model2.ub[md_model2.idRs.index(rxn)] > 0):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if gbr2_rH[rxn] == 2:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			if gbr2_rH[rxn] == 0:
				lines2_join = re.sub(color.group(0),'ff0000ff', lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

for pathway in pathway_list:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_MBA_%s.xml' % fluxstate1,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(cond1freqavg[rxn]),reference[o])
			if (md_model1.lb[md_model1.idRs.index(rxn)] < 0 and md_model1.ub[md_model1.idRs.index(rxn)] > 0):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if fbr1_hfr[rxn] == 2:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			if fbr1_hfr[rxn] == 0:
				lines2_join = re.sub(color.group(0),'ff0000ff', lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

for pathway in pathway_list:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_MBA_%s.xml' % fluxstate2,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(cond2freqavg[rxn]),reference[o])
			if (md_model2.lb[md_model2.idRs.index(rxn)] < 0 and md_model2.ub[md_model2.idRs.index(rxn)] > 0):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if fbr2_hfr[rxn] == 2:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			if fbr2_hfr[rxn] == 0:
				lines2_join = re.sub(color.group(0),'ff0000ff', lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

for pathway in pathway_list:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_Flux_%s.xml' % fluxstate1,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(cond1fluxavg[rxn]),reference[o])
			rxn_NA = str(rxn[2:]) + '_NA'
			if (md_model1.lb[md_model1.idRs.index(rxn)] < 0 and md_model1.ub[md_model1.idRs.index(rxn)] > 0 and fluxavgdict1[rxn] == 0 and cond1fluxavg[rxn] != rxn_NA):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if cond1fluxavg[rxn] == rxn_NA:
				lines2_join = re.sub(color.group(0),'fff0f0f0',lines2_join)
			elif fluxavgdict1[rxn] != 0:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			else:
				continue
			width = re.search('(?<=celldesigner:line\swidth=")(.*?\s)(?=color)',lines2_join)
			lines2_join = re.sub(str(width.group(0)),str(cond1fluxavgscaled[rxn])+'" ',lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

pathway_list_flux_cond1 = []
for pathway in pathway_list:
	pathway_list_flux_cond1.append(pathway + '_Flux_%s' % fluxstate1)

for pathway in pathway_list_flux_cond1:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_edit.xml','w')

	for o in range(0, len(reference)):
		flip = re.search('__', reference[o])
		if flip is None:
			continue
		else:
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined, re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			rectangleindex = re.search('(?<=rectangleIndex=")(.*?)(?=")',lines2_join)
			annotationsubsection = re.findall('(?<=reversible="false">\n)(.*?)(?=<listOfReactants>)',lines2_join, re.DOTALL)
			annotationsubsection_join = ''.join(annotationsubsection)
			countnannotationsubsection_join = annotationsubsection_join.count("\n")
			reactantssubsection = re.findall('(?<=</celldesigner:reactionType>\n)(.*?)(?=<celldesigner:baseProducts>)',lines2_join, re.DOTALL)
			reactantssubsection_join = ''.join(reactantssubsection)
			reactantssubsection_join = re.sub("baseReactant","baseProduct",reactantssubsection_join)
			productssubsection = re.findall('(?<=</celldesigner:baseReactants>\n)(.*?</celldesigner:baseProducts>\n)(?=<)',lines2_join, re.DOTALL)
			productssubsection_join = ''.join(productssubsection)
			productssubsection_join = re.sub("baseProduct","baseReactant",productssubsection_join)
			productsandreactants_join = ''.join((productssubsection_join,reactantssubsection_join))
			countnproductsandreactants = productsandreactants_join.count("\n")
			productsandreactantssplit = productsandreactants_join.split('\n')
			for p in range(0, countnproductsandreactants-1):
				reference[o+4+p] = ''.join((productsandreactantssplit[p],"\n"))
			reactantlinks = re.findall('(?<=<celldesigner:listOfReactantLinks>\n)(.*?)(?=</celldesigner:listOfReactantLinks>)',lines2_join, re.DOTALL)
			reactantlinks_join = ''.join(reactantlinks)
			if len(reactantlinks_join) > 0:
				reactantlinks_join = ''.join(("<celldesigner:listOfProductLinks>\n",reactantlinks_join,"</celldesigner:listOfProductLinks>\n"))
			reactantlinks_join = re.sub("reactant","product",reactantlinks_join)
			productlinks = re.findall('(?<=<celldesigner:listOfProductLinks>\n)(.*?)(?=</celldesigner:listOfProductLinks>)',lines2_join, re.DOTALL)
			productlinks_join = ''.join(productlinks)
			if len(productlinks_join) > 0:
				productlinks_join = ''.join(("<celldesigner:listOfReactantLinks>\n",productlinks_join,"</celldesigner:listOfReactantLinks>\n"))
			productlinks_join = re.sub("product","reactant",productlinks_join)
			annotation_end = re.findall('(?<=rectangleIndex)(.*?)(?=</annotation>)',lines2_join,re.DOTALL)
			annotation_end_join = ''.join(annotation_end)			
			editPoints_reactants = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',reactantlinks_join, re.DOTALL)					
			editPoints_reactants_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',reactantlinks_join)
			if len(editPoints_reactants) > 0:
				editPoints_reactants_split = {}
				editPoints_reactants_split2 = {}
				editPoints_reactants_index = editPoints_reactants[-1]
				editPoints_reactants_split = re.split('\s|,',editPoints_reactants_index)
				for r in range(0, len(editPoints_reactants_split[1::2])):
					if r == len(editPoints_reactants_split[1::2]) - 1:
						editPoints_reactants_split2[2+r*4-1] = -float(editPoints_reactants_split[1+r*2])
						editPoints_reactants_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_reactants_split2[2] = -float(editPoints_reactants_split[1])
					else:
						editPoints_reactants_split2[2+r*4-1] = -float(editPoints_reactants_split[1+r*2])
						editPoints_reactants_split2[3+r*4-1] = ' ' 
				for s in range(0, len(editPoints_reactants_split[::2])):
					if s == 0:
						editPoints_reactants_split2[0] = 1-float(editPoints_reactants_split[0])
						editPoints_reactants_split2[1] = ','
						editPoints_reactants_split2[2] = -float(editPoints_reactants_split[1])
					else:
						editPoints_reactants_split2[s*4-1] = 1-float(editPoints_reactants_split[s*2])
						editPoints_reactants_split2[1+s*4-1] = ','
				editPoints_reactants_split2_list = []
				for i, v in enumerate(editPoints_reactants_split2.keys()):
					editPoints_reactants_split2_list.append(editPoints_reactants_split2[i])
				editPoints_reactants_split3_list = {}
				for r in range(0, len(editPoints_reactants_split2_list[1::2])):
					if len(editPoints_reactants_split2_list) == 3:
						rounded_r = 0
					if len(editPoints_reactants_split2_list) > 3:
						rounded_r = int((len(editPoints_reactants_split2_list)-3)/4)
					if rounded_r == 0: 
						editPoints_reactants_split3_list[0] = editPoints_reactants_split2_list[0]
						editPoints_reactants_split3_list[1] = editPoints_reactants_split2_list[1]
						editPoints_reactants_split3_list[2] = editPoints_reactants_split2_list[2]
					if rounded_r > 0:
						editPoints_reactants_split3_list[0] = editPoints_reactants_split2_list[3+4*(rounded_r-1)]
						editPoints_reactants_split3_list[1] = editPoints_reactants_split2_list[4+4*(rounded_r-1)]
						editPoints_reactants_split3_list[2] = editPoints_reactants_split2_list[5+4*(rounded_r-1)]
						editPoints_reactants_split3_list[3] = editPoints_reactants_split2_list[6+4*(rounded_r-1)]
						for i in range(0,rounded_r-1):
							editPoints_reactants_split3_list[4+i*4] = editPoints_reactants_split2_list[3+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[5+i*4] = editPoints_reactants_split2_list[4+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[6+i*4] = editPoints_reactants_split2_list[5+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[7+i*4] = editPoints_reactants_split2_list[6+4*(rounded_r-1-i)]	
						editPoints_reactants_split3_list[4+(rounded_r-1)*4] = editPoints_reactants_split2_list[0]
						editPoints_reactants_split3_list[5+(rounded_r-1)*4] = editPoints_reactants_split2_list[1]
						editPoints_reactants_split3_list[6+(rounded_r-1)*4] = editPoints_reactants_split2_list[2]
				editPoints_reactants_split4_list = []
				editPoints_reactants_split4_string = ""
				for i, v in enumerate(editPoints_reactants_split3_list.keys()):
					editPoints_reactants_split4_list.append(editPoints_reactants_split3_list[i])
					editPoints_reactants_split4_string += str(editPoints_reactants_split4_list[i])	
				reactantlinks_join = re.sub(editPoints_reactants_search.group(0),editPoints_reactants_split4_string,reactantlinks_join)	
			editPoints_products = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',productlinks_join, re.DOTALL)
			editPoints_products_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',productlinks_join)
			if len(editPoints_products) > 0:
				editPoints_products_split = {}
				editPoints_products_split2 = {}
				editPoints_products_index = editPoints_products[-1]
				editPoints_products_split = re.split('\s|,',editPoints_products_index)
				for r in range(0, len(editPoints_products_split[1::2])):
					if r == len(editPoints_products_split[1::2]) - 1:
						editPoints_products_split2[2+r*4-1] = -float(editPoints_products_split[1+r*2])
						editPoints_products_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_products_split2[2] = -float(editPoints_products_split[1])
					else:
						editPoints_products_split2[2+r*4-1] = -float(editPoints_products_split[1+r*2])
						editPoints_products_split2[3+r*4-1] = ' ' 
				for s in range(0, len(editPoints_products_split[::2])):
					if s == 0:
						editPoints_products_split2[0] = 1-float(editPoints_products_split[0])
						editPoints_products_split2[1] = ','
						editPoints_products_split2[2] = -float(editPoints_products_split[1])
					else:
						editPoints_products_split2[s*4-1] = 1-float(editPoints_products_split[s*2])
						editPoints_products_split2[1+s*4-1] = ','
				editPoints_products_split2_list = []
				for i, v in enumerate(editPoints_products_split2.keys()):
					editPoints_products_split2_list.append(editPoints_products_split2[i])
				editPoints_products_split3_list = {}
				for r in range(0, len(editPoints_products_split2_list[1::2])):
					if len(editPoints_products_split2_list) == 3:
						rounded_r = 0
					if len(editPoints_products_split2_list) > 3:
						rounded_r = int((len(editPoints_products_split2_list)-3)/4)
					if rounded_r == 0: 
						editPoints_products_split3_list[0] = editPoints_products_split2_list[0]
						editPoints_products_split3_list[1] = editPoints_products_split2_list[1]
						editPoints_products_split3_list[2] = editPoints_products_split2_list[2]
					if rounded_r > 0:
						editPoints_products_split3_list[0] = editPoints_products_split2_list[3+4*(rounded_r-1)]
						editPoints_products_split3_list[1] = editPoints_products_split2_list[4+4*(rounded_r-1)]
						editPoints_products_split3_list[2] = editPoints_products_split2_list[5+4*(rounded_r-1)]
						editPoints_products_split3_list[3] = editPoints_products_split2_list[6+4*(rounded_r-1)]
						for i in range(0,rounded_r-1):
							editPoints_products_split3_list[4+i*4] = editPoints_products_split2_list[3+4*(rounded_r-1-i)]
							editPoints_products_split3_list[5+i*4] = editPoints_products_split2_list[4+4*(rounded_r-1-i)]
							editPoints_products_split3_list[6+i*4] = editPoints_products_split2_list[5+4*(rounded_r-1-i)]
							editPoints_products_split3_list[7+i*4] = editPoints_products_split2_list[6+4*(rounded_r-1-i)]	
						editPoints_products_split3_list[4+(rounded_r-1)*4] = editPoints_products_split2_list[0]
						editPoints_products_split3_list[5+(rounded_r-1)*4] = editPoints_products_split2_list[1]
						editPoints_products_split3_list[6+(rounded_r-1)*4] = editPoints_products_split2_list[2]
				editPoints_products_split4_list = []
				editPoints_products_split4_string = ""
				for i, v in enumerate(editPoints_products_split3_list.keys()):
					editPoints_products_split4_list.append(editPoints_products_split3_list[i])
					editPoints_products_split4_string += str(editPoints_products_split4_list[i])
				productlinks_join = re.sub(editPoints_products_search.group(0),editPoints_products_split4_string,productlinks_join)
			editPoints_last = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',annotation_end_join, re.DOTALL)
			editPoints_last_total = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
			if len(editPoints_last) > 0:
				editPoints_last_split = {}
				editPoints_last_split2 = {}
				editPoints_last_index = editPoints_last[-1]
				editPoints_last_split = re.split('\s|,',editPoints_last_index)
				for r in range(0, len(editPoints_last_split[1::2])):
					if r == len(editPoints_last_split[1::2]) - 1:
						editPoints_last_split2[2+r*4-1] = -float(editPoints_last_split[1+r*2])
						editPoints_last_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_last_split2[2] = -float(editPoints_last_split[1])
					else:
						editPoints_last_split2[2+r*4-1] = -float(editPoints_last_split[1+r*2])
						editPoints_last_split2[3+r*4-1] = ' '
				for s in range(0, len(editPoints_last_split[::2])):
					if s == 0:
						editPoints_last_split2[0] = 1-float(editPoints_last_split[0])
						editPoints_last_split2[1] = ','
						editPoints_last_split2[2] = -float(editPoints_last_split[1])
					else:
						editPoints_last_split2[s*4-1] = 1-float(editPoints_last_split[s*2])
						editPoints_last_split2[1+s*4-1] = ','
				editPoints_last_split2_list = []
				for i, v in enumerate(editPoints_last_split2.keys()):
					editPoints_last_split2_list.append(editPoints_last_split2[i])
				editPoints_last_split3_list = {}
				for r in range(0, len(editPoints_last_split2_list[1::2])):
						if len(editPoints_last_split2_list) == 3:
							rounded_r = 0
						if len(editPoints_last_split2_list) > 3:
							rounded_r = int((len(editPoints_last_split2_list)-3)/4)
						if rounded_r == 0: 
							editPoints_last_split3_list[0] = editPoints_last_split2_list[0]
							editPoints_last_split3_list[1] = editPoints_last_split2_list[1]
							editPoints_last_split3_list[2] = editPoints_last_split2_list[2]
						if rounded_r > 0:
							editPoints_last_split3_list[0] = editPoints_last_split2_list[3+4*(rounded_r-1)]
							editPoints_last_split3_list[1] = editPoints_last_split2_list[4+4*(rounded_r-1)]
							editPoints_last_split3_list[2] = editPoints_last_split2_list[5+4*(rounded_r-1)]
							editPoints_last_split3_list[3] = editPoints_last_split2_list[6+4*(rounded_r-1)]
							for i in range(0,rounded_r-1):
								editPoints_last_split3_list[4+i*4] = editPoints_last_split2_list[3+4*(rounded_r-1-i)]
								editPoints_last_split3_list[5+i*4] = editPoints_last_split2_list[4+4*(rounded_r-1-i)]
								editPoints_last_split3_list[6+i*4] = editPoints_last_split2_list[5+4*(rounded_r-1-i)]
								editPoints_last_split3_list[7+i*4] = editPoints_last_split2_list[6+4*(rounded_r-1-i)]	
							editPoints_last_split3_list[4+(rounded_r-1)*4] = editPoints_last_split2_list[0]
							editPoints_last_split3_list[5+(rounded_r-1)*4] = editPoints_last_split2_list[1]
							editPoints_last_split3_list[6+(rounded_r-1)*4] = editPoints_last_split2_list[2]
				editPoints_last_split4_list = []
				editPoints_last_split4_string = ""
				for i, v in enumerate(editPoints_last_split3_list.keys()):
					editPoints_last_split4_list.append(editPoints_last_split3_list[i])
					editPoints_last_split4_string += str(editPoints_last_split4_list[i])
				editPoints_until = re.findall('(?<=reaction\smetaid=")(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
				editPoints_until_join = ''.join(editPoints_until)
				editPoints_until_more = re.findall('(?<=</celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
				editPoints_last_searchable = re.findall('(?<=rectangleIndex)(.*?</celldesigner:editPoints>)(?=\n)',lines2_join, re.DOTALL)
				editPoints_last_searchable_join_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',annotation_end_join)
				countnedit1 = editPoints_until_join.count("\n")
				if len(editPoints_last_total) > 1:
					countnedit2 = editPoints_until_more[0].count("\n")
					countnedit2_1 = countnedit1 + countnedit2 
				if len(editPoints_last_total) > 2:
					countnedit2 = editPoints_until_more[0].count("\n")
					countnedit3 = editPoints_until_more[1].count("\n")
					countnedit3_2_1 = countnedit1 + countnedit2 + countnedit3
				if len(editPoints_last_total) == 1:
					reference[o+countnedit1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit1])
				if len(editPoints_last_total) == 2:
					reference[o+countnedit2_1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit2_1])
				if len(editPoints_last_total) == 3:
					reference[o+countnedit3_2_1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit3_2_1])
			reactantsandproductlinks_join = ''.join((productlinks_join,reactantlinks_join))
			if len(reactantsandproductlinks_join) > 0:	
				countnreactantsandproductlinks = reactantsandproductlinks_join.count("\n")
				reactantsandproductlinkssplit = reactantsandproductlinks_join.split('\n')
				for q in range(0, countnreactantsandproductlinks):
					reference[o+4+countnproductsandreactants+q] = ''.join((reactantsandproductlinkssplit[q],"\n"))
			listofreactants = re.findall('(?<=<listOfReactants>\n)(.*?)(?=</listOfReactants)',lines2_join, re.DOTALL)
			listofreactants_join = ''.join(listofreactants)
			listofreactants_join_rename = ''.join(("<listOfReactants>\n",listofreactants_join)) 
			listofreactants_join_rename = re.sub("Reactant","Product",listofreactants_join_rename)
			listofproducts = re.findall('(?<=<listOfProducts>\n)(.*?)(?=</listOfProducts)',lines2_join,re.DOTALL)
			listofproducts_join = ''.join(listofproducts)
			listofproducts_join_rename = ''.join(("<listOfProducts>\n",listofproducts_join,"</listOfProducts>\n"))
			listofproducts_join_rename = re.sub("Product","Reactant",listofproducts_join_rename)
			listofreactantsandproducts_join = ''.join((listofproducts_join_rename,listofreactants_join_rename))
			countnlistofreactantsandproducts = listofreactantsandproducts_join.count("\n")
			listofreactantsandproductssplit = listofreactantsandproducts_join.split('\n')
			for r in range(0, countnlistofreactantsandproducts):
				reference[o+countnannotationsubsection_join+1+r] = ''.join((listofreactantsandproductssplit[r],"\n"))

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

for pathway in pathway_list: 
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_Flux_%s.xml' % fluxstate2,'w')

	for o in range(0, len(reference)):
		reaction_name = re.search('(?<=reaction\smetaid=")(.*?)(?="\sid)',reference[o])
		reaction_id = re.search('(?<="\sid=")(.*?)(?="\sreversible)', reference[o])
		reaction_id_name = re.search('(?<="\sname=")(.*?)(?="\sreversible)', reference[o])
		if (reaction_id is None) and (reaction_id_name is None):
			continue
		else:
			if reaction_id_name is None:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id.group(0))
				else:
					rxn = str(reaction_id.group(0))
			else:
				if str(reaction_id.group(0))[:2] != 'R_':
					rxn = 'R_' + str(reaction_id_name.group(0))
				else:
					rxn = str(reaction_id_name.group(0))
			reference[o] = re.sub(reaction_id.group(0),str(cond2fluxavg[rxn]),reference[o])
			rxn_NA = str(rxn[2:]) + '_NA'
			if (md_model1.lb[md_model2.idRs.index(rxn)] < 0 and md_model2.ub[md_model2.idRs.index(rxn)] > 0 and fluxavgdict2[rxn] == 0 and cond2fluxavg[rxn] != rxn_NA):
				reversible = re.search('(?<=reversible=")(.*?)(?=">)',reference[o])
				reference[o] = re.sub(reversible.group(0),'true',reference[o])
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=<reaction)',linesjoined, re.DOTALL)
			if len(reactionlength) == 0:
				reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined,re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			color = re.search('(?<=color=")(.*?)(?=")',lines2_join)
			if (cond2fluxavg[rxn] == rxn_NA):
				lines2_join = re.sub(color.group(0),'fff0f0f0',lines2_join)
			elif fluxavgdict2[rxn] != 0:
				lines2_join = re.sub(color.group(0),'ffff0000',lines2_join)
			else:
				continue
			width = re.search('(?<=celldesigner:line\swidth=")(.*?\s)(?=color)',lines2_join)
			lines2_join = re.sub(str(width.group(0)),str(cond2fluxavgscaled[rxn])+'" ',lines2_join)
			lines2_join_list = re.split('\n', lines2_join)
			lines2_join_list_remake = []
			for i in lines2_join_list:
				string = i + '\n'
				lines2_join_list_remake.append(string)
			count = -1 
			del lines2_join_list_remake[-1]
			for i in lines2_join_list_remake:
				count += 1
				reference[o+count] = re.sub(reference[o+count],i,reference[o+count])
	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()

pathway_list_flux_cond2 = []
for pathway in pathway_list:
	pathway_list_flux_cond2.append(pathway + '_Flux_%s' % fluxstate2)

for pathway in pathway_list_flux_cond2:
	file_name_reference = open(pathway + '.xml')
	reference = file_name_reference.readlines()
	file_name_reference_out = open(pathway + '_edit.xml','w')

	for o in range(0, len(reference)):
		flip = re.search('__', reference[o])
		if flip is None:
			continue
		else:
			linesjoin = reference[o:o+300]
			linesjoined = ''.join(linesjoin)
			reactionlength = re.findall('(?<=>)(.*?)(?=</reaction)',linesjoined, re.DOTALL)
			reactionlength_join = ''.join(reactionlength[0])
			countnreactionlength_join = reactionlength_join.count("\n")
			lines2 = reference[o:o+countnreactionlength_join]
			lines2_join = ''.join(lines2)
			rectangleindex = re.search('(?<=rectangleIndex=")(.*?)(?=")',lines2_join)
			annotationsubsection = re.findall('(?<=reversible="false">\n)(.*?)(?=<listOfReactants>)',lines2_join, re.DOTALL)
			annotationsubsection_join = ''.join(annotationsubsection)
			countnannotationsubsection_join = annotationsubsection_join.count("\n")
			reactantssubsection = re.findall('(?<=</celldesigner:reactionType>\n)(.*?)(?=<celldesigner:baseProducts>)',lines2_join, re.DOTALL)
			reactantssubsection_join = ''.join(reactantssubsection)
			reactantssubsection_join = re.sub("baseReactant","baseProduct",reactantssubsection_join)
			productssubsection = re.findall('(?<=</celldesigner:baseReactants>\n)(.*?</celldesigner:baseProducts>\n)(?=<)',lines2_join, re.DOTALL)
			productssubsection_join = ''.join(productssubsection)
			productssubsection_join = re.sub("baseProduct","baseReactant",productssubsection_join)
			productsandreactants_join = ''.join((productssubsection_join,reactantssubsection_join))
			countnproductsandreactants = productsandreactants_join.count("\n")
			productsandreactantssplit = productsandreactants_join.split('\n')
			for p in range(0, countnproductsandreactants-1):
				reference[o+4+p] = ''.join((productsandreactantssplit[p],"\n"))
			reactantlinks = re.findall('(?<=<celldesigner:listOfReactantLinks>\n)(.*?)(?=</celldesigner:listOfReactantLinks>)',lines2_join, re.DOTALL)
			reactantlinks_join = ''.join(reactantlinks)
			if len(reactantlinks_join) > 0:
				reactantlinks_join = ''.join(("<celldesigner:listOfProductLinks>\n",reactantlinks_join,"</celldesigner:listOfProductLinks>\n"))
			reactantlinks_join = re.sub("reactant","product",reactantlinks_join)
			productlinks = re.findall('(?<=<celldesigner:listOfProductLinks>\n)(.*?)(?=</celldesigner:listOfProductLinks>)',lines2_join, re.DOTALL)
			productlinks_join = ''.join(productlinks)
			if len(productlinks_join) > 0:
				productlinks_join = ''.join(("<celldesigner:listOfReactantLinks>\n",productlinks_join,"</celldesigner:listOfReactantLinks>\n"))
			productlinks_join = re.sub("product","reactant",productlinks_join)
			annotation_end = re.findall('(?<=rectangleIndex)(.*?)(?=</annotation>)',lines2_join,re.DOTALL)
			annotation_end_join = ''.join(annotation_end)			
			editPoints_reactants = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',reactantlinks_join, re.DOTALL)					
			editPoints_reactants_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',reactantlinks_join)
			if len(editPoints_reactants) > 0:
				editPoints_reactants_split = {}
				editPoints_reactants_split2 = {}
				editPoints_reactants_index = editPoints_reactants[-1]
				editPoints_reactants_split = re.split('\s|,',editPoints_reactants_index)
				for r in range(0, len(editPoints_reactants_split[1::2])):
					if r == len(editPoints_reactants_split[1::2]) - 1:
						editPoints_reactants_split2[2+r*4-1] = -float(editPoints_reactants_split[1+r*2])
						editPoints_reactants_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_reactants_split2[2] = -float(editPoints_reactants_split[1])
					else:
						editPoints_reactants_split2[2+r*4-1] = -float(editPoints_reactants_split[1+r*2])
						editPoints_reactants_split2[3+r*4-1] = ' ' 
				for s in range(0, len(editPoints_reactants_split[::2])):
					if s == 0:
						editPoints_reactants_split2[0] = 1-float(editPoints_reactants_split[0])
						editPoints_reactants_split2[1] = ','
						editPoints_reactants_split2[2] = -float(editPoints_reactants_split[1])
					else:
						editPoints_reactants_split2[s*4-1] = 1-float(editPoints_reactants_split[s*2])
						editPoints_reactants_split2[1+s*4-1] = ','
				editPoints_reactants_split2_list = []
				for i, v in enumerate(editPoints_reactants_split2.keys()):
					editPoints_reactants_split2_list.append(editPoints_reactants_split2[i])
				editPoints_reactants_split3_list = {}
				for r in range(0, len(editPoints_reactants_split2_list[1::2])):
					if len(editPoints_reactants_split2_list) == 3:
						rounded_r = 0
					if len(editPoints_reactants_split2_list) > 3:
						rounded_r = int((len(editPoints_reactants_split2_list)-3)/4)
					if rounded_r == 0: 
						editPoints_reactants_split3_list[0] = editPoints_reactants_split2_list[0]
						editPoints_reactants_split3_list[1] = editPoints_reactants_split2_list[1]
						editPoints_reactants_split3_list[2] = editPoints_reactants_split2_list[2]
					if rounded_r > 0:
						editPoints_reactants_split3_list[0] = editPoints_reactants_split2_list[3+4*(rounded_r-1)]
						editPoints_reactants_split3_list[1] = editPoints_reactants_split2_list[4+4*(rounded_r-1)]
						editPoints_reactants_split3_list[2] = editPoints_reactants_split2_list[5+4*(rounded_r-1)]
						editPoints_reactants_split3_list[3] = editPoints_reactants_split2_list[6+4*(rounded_r-1)]
						for i in range(0,rounded_r-1):
							editPoints_reactants_split3_list[4+i*4] = editPoints_reactants_split2_list[3+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[5+i*4] = editPoints_reactants_split2_list[4+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[6+i*4] = editPoints_reactants_split2_list[5+4*(rounded_r-1-i)]
							editPoints_reactants_split3_list[7+i*4] = editPoints_reactants_split2_list[6+4*(rounded_r-1-i)]	
						editPoints_reactants_split3_list[4+(rounded_r-1)*4] = editPoints_reactants_split2_list[0]
						editPoints_reactants_split3_list[5+(rounded_r-1)*4] = editPoints_reactants_split2_list[1]
						editPoints_reactants_split3_list[6+(rounded_r-1)*4] = editPoints_reactants_split2_list[2]
				editPoints_reactants_split4_list = []
				editPoints_reactants_split4_string = ""
				for i, v in enumerate(editPoints_reactants_split3_list.keys()):
					editPoints_reactants_split4_list.append(editPoints_reactants_split3_list[i])
					editPoints_reactants_split4_string += str(editPoints_reactants_split4_list[i])	
				reactantlinks_join = re.sub(editPoints_reactants_search.group(0),editPoints_reactants_split4_string,reactantlinks_join)	
			editPoints_products = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',productlinks_join, re.DOTALL)
			editPoints_products_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',productlinks_join)
			if len(editPoints_products) > 0:
				editPoints_products_split = {}
				editPoints_products_split2 = {}
				editPoints_products_index = editPoints_products[-1]
				editPoints_products_split = re.split('\s|,',editPoints_products_index)
				for r in range(0, len(editPoints_products_split[1::2])):
					if r == len(editPoints_products_split[1::2]) - 1:
						editPoints_products_split2[2+r*4-1] = -float(editPoints_products_split[1+r*2])
						editPoints_products_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_products_split2[2] = -float(editPoints_products_split[1])
					else:
						editPoints_products_split2[2+r*4-1] = -float(editPoints_products_split[1+r*2])
						editPoints_products_split2[3+r*4-1] = ' ' 
				for s in range(0, len(editPoints_products_split[::2])):
					if s == 0:
						editPoints_products_split2[0] = 1-float(editPoints_products_split[0])
						editPoints_products_split2[1] = ','
						editPoints_products_split2[2] = -float(editPoints_products_split[1])
					else:
						editPoints_products_split2[s*4-1] = 1-float(editPoints_products_split[s*2])
						editPoints_products_split2[1+s*4-1] = ','
				editPoints_products_split2_list = []
				for i, v in enumerate(editPoints_products_split2.keys()):
					editPoints_products_split2_list.append(editPoints_products_split2[i])
				editPoints_products_split3_list = {}
				for r in range(0, len(editPoints_products_split2_list[1::2])):
					if len(editPoints_products_split2_list) == 3:
						rounded_r = 0
					if len(editPoints_products_split2_list) > 3:
						rounded_r = int((len(editPoints_products_split2_list)-3)/4)
					if rounded_r == 0: 
						editPoints_products_split3_list[0] = editPoints_products_split2_list[0]
						editPoints_products_split3_list[1] = editPoints_products_split2_list[1]
						editPoints_products_split3_list[2] = editPoints_products_split2_list[2]
					if rounded_r > 0:
						editPoints_products_split3_list[0] = editPoints_products_split2_list[3+4*(rounded_r-1)]
						editPoints_products_split3_list[1] = editPoints_products_split2_list[4+4*(rounded_r-1)]
						editPoints_products_split3_list[2] = editPoints_products_split2_list[5+4*(rounded_r-1)]
						editPoints_products_split3_list[3] = editPoints_products_split2_list[6+4*(rounded_r-1)]
						for i in range(0,rounded_r-1):
							editPoints_products_split3_list[4+i*4] = editPoints_products_split2_list[3+4*(rounded_r-1-i)]
							editPoints_products_split3_list[5+i*4] = editPoints_products_split2_list[4+4*(rounded_r-1-i)]
							editPoints_products_split3_list[6+i*4] = editPoints_products_split2_list[5+4*(rounded_r-1-i)]
							editPoints_products_split3_list[7+i*4] = editPoints_products_split2_list[6+4*(rounded_r-1-i)]	
						editPoints_products_split3_list[4+(rounded_r-1)*4] = editPoints_products_split2_list[0]
						editPoints_products_split3_list[5+(rounded_r-1)*4] = editPoints_products_split2_list[1]
						editPoints_products_split3_list[6+(rounded_r-1)*4] = editPoints_products_split2_list[2]
				editPoints_products_split4_list = []
				editPoints_products_split4_string = ""
				for i, v in enumerate(editPoints_products_split3_list.keys()):
					editPoints_products_split4_list.append(editPoints_products_split3_list[i])
					editPoints_products_split4_string += str(editPoints_products_split4_list[i])
				productlinks_join = re.sub(editPoints_products_search.group(0),editPoints_products_split4_string,productlinks_join)
			editPoints_last = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',annotation_end_join, re.DOTALL)
			editPoints_last_total = re.findall('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
			if len(editPoints_last) > 0:
				editPoints_last_split = {}
				editPoints_last_split2 = {}
				editPoints_last_index = editPoints_last[-1]
				editPoints_last_split = re.split('\s|,',editPoints_last_index)
				for r in range(0, len(editPoints_last_split[1::2])):
					if r == len(editPoints_last_split[1::2]) - 1:
						editPoints_last_split2[2+r*4-1] = -float(editPoints_last_split[1+r*2])
						editPoints_last_split2[3+r*4-1] = ' '
					elif r == 0:
						editPoints_last_split2[2] = -float(editPoints_last_split[1])
					else:
						editPoints_last_split2[2+r*4-1] = -float(editPoints_last_split[1+r*2])
						editPoints_last_split2[3+r*4-1] = ' '
				for s in range(0, len(editPoints_last_split[::2])):
					if s == 0:
						editPoints_last_split2[0] = 1-float(editPoints_last_split[0])
						editPoints_last_split2[1] = ','
						editPoints_last_split2[2] = -float(editPoints_last_split[1])
					else:
						editPoints_last_split2[s*4-1] = 1-float(editPoints_last_split[s*2])
						editPoints_last_split2[1+s*4-1] = ','
				editPoints_last_split2_list = []
				for i, v in enumerate(editPoints_last_split2.keys()):
					editPoints_last_split2_list.append(editPoints_last_split2[i])
				editPoints_last_split3_list = {}
				for r in range(0, len(editPoints_last_split2_list[1::2])):
						if len(editPoints_last_split2_list) == 3:
							rounded_r = 0
						if len(editPoints_last_split2_list) > 3:
							rounded_r = int((len(editPoints_last_split2_list)-3)/4)
						if rounded_r == 0: 
							editPoints_last_split3_list[0] = editPoints_last_split2_list[0]
							editPoints_last_split3_list[1] = editPoints_last_split2_list[1]
							editPoints_last_split3_list[2] = editPoints_last_split2_list[2]
						if rounded_r > 0:
							editPoints_last_split3_list[0] = editPoints_last_split2_list[3+4*(rounded_r-1)]
							editPoints_last_split3_list[1] = editPoints_last_split2_list[4+4*(rounded_r-1)]
							editPoints_last_split3_list[2] = editPoints_last_split2_list[5+4*(rounded_r-1)]
							editPoints_last_split3_list[3] = editPoints_last_split2_list[6+4*(rounded_r-1)]
							for i in range(0,rounded_r-1):
								editPoints_last_split3_list[4+i*4] = editPoints_last_split2_list[3+4*(rounded_r-1-i)]
								editPoints_last_split3_list[5+i*4] = editPoints_last_split2_list[4+4*(rounded_r-1-i)]
								editPoints_last_split3_list[6+i*4] = editPoints_last_split2_list[5+4*(rounded_r-1-i)]
								editPoints_last_split3_list[7+i*4] = editPoints_last_split2_list[6+4*(rounded_r-1-i)]	
							editPoints_last_split3_list[4+(rounded_r-1)*4] = editPoints_last_split2_list[0]
							editPoints_last_split3_list[5+(rounded_r-1)*4] = editPoints_last_split2_list[1]
							editPoints_last_split3_list[6+(rounded_r-1)*4] = editPoints_last_split2_list[2]
				editPoints_last_split4_list = []
				editPoints_last_split4_string = ""
				for i, v in enumerate(editPoints_last_split3_list.keys()):
					editPoints_last_split4_list.append(editPoints_last_split3_list[i])
					editPoints_last_split4_string += str(editPoints_last_split4_list[i])
				editPoints_until = re.findall('(?<=reaction\smetaid=")(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
				editPoints_until_join = ''.join(editPoints_until)
				editPoints_until_more = re.findall('(?<=</celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',lines2_join, re.DOTALL)
				editPoints_last_searchable = re.findall('(?<=rectangleIndex)(.*?</celldesigner:editPoints>)(?=\n)',lines2_join, re.DOTALL)
				editPoints_last_searchable_join_search = re.search('(?<=<celldesigner:editPoints>)(.*?)(?=</celldesigner:editPoints>)',annotation_end_join)
				countnedit1 = editPoints_until_join.count("\n")
				if len(editPoints_last_total) > 1:
					countnedit2 = editPoints_until_more[0].count("\n")
					countnedit2_1 = countnedit1 + countnedit2 
				if len(editPoints_last_total) > 2:
					countnedit2 = editPoints_until_more[0].count("\n")
					countnedit3 = editPoints_until_more[1].count("\n")
					countnedit3_2_1 = countnedit1 + countnedit2 + countnedit3
				if len(editPoints_last_total) == 1:
					reference[o+countnedit1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit1])
				if len(editPoints_last_total) == 2:
					reference[o+countnedit2_1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit2_1])
				if len(editPoints_last_total) == 3:
					reference[o+countnedit3_2_1] = re.sub(editPoints_last_searchable_join_search.group(0),editPoints_last_split4_string,reference[o+countnedit3_2_1])
			reactantsandproductlinks_join = ''.join((productlinks_join,reactantlinks_join))
			if len(reactantsandproductlinks_join) > 0:	
				countnreactantsandproductlinks = reactantsandproductlinks_join.count("\n")
				reactantsandproductlinkssplit = reactantsandproductlinks_join.split('\n')
				for q in range(0, countnreactantsandproductlinks):
					reference[o+4+countnproductsandreactants+q] = ''.join((reactantsandproductlinkssplit[q],"\n"))
			listofreactants = re.findall('(?<=<listOfReactants>\n)(.*?)(?=</listOfReactants)',lines2_join, re.DOTALL)
			listofreactants_join = ''.join(listofreactants)
			listofreactants_join_rename = ''.join(("<listOfReactants>\n",listofreactants_join)) 
			listofreactants_join_rename = re.sub("Reactant","Product",listofreactants_join_rename)
			listofproducts = re.findall('(?<=<listOfProducts>\n)(.*?)(?=</listOfProducts)',lines2_join,re.DOTALL)
			listofproducts_join = ''.join(listofproducts)
			listofproducts_join_rename = ''.join(("<listOfProducts>\n",listofproducts_join,"</listOfProducts>\n"))
			listofproducts_join_rename = re.sub("Product","Reactant",listofproducts_join_rename)
			listofreactantsandproducts_join = ''.join((listofproducts_join_rename,listofreactants_join_rename))
			countnlistofreactantsandproducts = listofreactantsandproducts_join.count("\n")
			listofreactantsandproductssplit = listofreactantsandproducts_join.split('\n')
			for r in range(0, countnlistofreactantsandproducts):
				reference[o+countnannotationsubsection_join+1+r] = ''.join((listofreactantsandproductssplit[r],"\n"))

	reference2 = reference	
	for p, item3 in enumerate(reference2):
		file_name_reference_out.write(reference2[p])
	file_name_reference_out.close()
