from cobra.io import read_sbml_model
from cobra.core import ArrayBasedModel
from scipy.sparse import coo_matrix
import cobra

import os
import sys
import time
from numpy import array, delete, shape, object as np_object
import numpy as np
from scipy.io import loadmat
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix
import re
import operator
from libsbml import SBMLDocument, SpeciesReference, KineticLaw, Parameter
from libsbml import readSBML, writeSBML
import collections
import csv
import argparse
import cPickle as pickle

#Define argparse command line inputs
parser = argparse.ArgumentParser(description='Convert SBML model (.xml or .mat file) into EXAMO pickle file')
parser.add_argument('s', nargs="+", type=str, help='Necessary variable: SBML file')
parser.add_argument('b', nargs="+", type=str, help='Necessary variable: biomass rxn; Append "R_" to the front if the reaction name does not begin with this')
parser.add_argument('e', nargs="+", type=str, help='Necessary variable: extracellular compartment abbreviation. Instead of brackets, use underscores (ex: "_e")')
parser.add_argument('-l', type=str, help='lb file')
parser.add_argument('-u', type=str, help='ub file')
parser.add_argument('-g', type=str, help='gene2rxn file')
parser.add_argument('-m', type=str, help='metabolite mapping complexes .pkl file. Must first create .pkl file if using this option')
parser.add_argument('-n', type=str, help='nucleotide conversions .pkl file. Must first create .pkl file if using this option')
parser.add_argument('-c', type=str, help='metabolite to carbon mapping file')

args = parser.parse_args()

args_b = str(args.b)
args_b = args_b[2:-2]

#exportPickle function from EXAMO package
def exportPickle(obj, fileName, mode = 'wb', protocol = -1):
    import pickle
    f = open(fileName, mode)
    pickle.dump(obj, f, protocol = -1)
    f.close()

#Modify the file name
args_s = str(args.s)
args_s = args_s[2:-2]

#Name of exported pickle file
locTime = time.localtime()
test_model = '%i%02i%02i%02i%02i%02i' % locTime[:6]
if args.l is not None:
	test_model = test_model + str('_l')
if args.u is not None:
	test_model = test_model + str('_u')
if args.g is not None:
	test_model = test_model + str('_g')
if args.m is not None:
	test_model = test_model + str('_m')
if args.n is not None:
	test_model = test_model + str('_n')
if args.c is not None:
	test_model = test_model + str('_c')

test_model = test_model + '_' + args_s[:-4] + str('.pkl')

#Create necessary variables and import the model
if args_s[-4:] == '.xml':
	cobra_model = read_sbml_model(args_s)
if args_s[-4:] == '.mat':
	cobra_model = cobra.io.mat.load_matlab_model(args_s)
idRs = []
lb = []
ub = []
pathway = []
rxn2lb = {}
rxn2ub = {}
gene2rxn = {}
rxns = {}
rxns_original = {}
rxns_to_remove = []	
rxns_to_remove_2 = []
rxns_to_remove_3 = []
unbalanced_rxns = []

idSp = []
for i in cobra_model.metabolites:
	metabolite_name = re.sub("LPAREN","",str(i.id))
	metabolite_name = re.sub("RPAREN","",metabolite_name)
	metabolite_name = re.sub("\(","_",metabolite_name)
	metabolite_name = re.sub("\)","_",metabolite_name)
	metabolite_name = re.sub("\[","_",metabolite_name)
	metabolite_name = re.sub("\]","_",metabolite_name)
	metabolite_name = re.sub("\-","_",metabolite_name)
	metabolite_name = re.sub("__","_",metabolite_name)
	idSp.append(str("M_")+str(metabolite_name))

seen_set = set([x for x in idSp if idSp.count(x) > 1])
seen_list = list(seen_set)

cobra_model = cobra.core.ArrayBasedModel(cobra_model)
#cobra_model = cobra.core.ArrayBasedModel.SMatrix_lil(ArrayBasedModel(cobra_model))
S = coo_matrix(cobra_model.S)
S = sparse.lil_matrix(S)

#Import lower boundary adjustments if the argument is supplied from the command line. 
if args.l is not None:
	args_l = str(args.l)
	lb_file = open(args_l)
	lb = lb_file.readlines()
	for i, item1 in enumerate(lb):
		lb[i] = item1.rstrip()
		lb[i] = float(lb[i])
	lb_file.close()
	for i, item1 in enumerate(cobra_model.reactions):
		reaction_name = re.sub("LPAREN","",str(item1))
		reaction_name = re.sub("RPAREN","",reaction_name)
		reaction_name = re.sub("\(","_",reaction_name)
		reaction_name = re.sub("\)","_",reaction_name)
		reaction_name = re.sub("\[","_",reaction_name)
		reaction_name = re.sub("\]","_",reaction_name)
		reaction_name = re.sub("\-","_",reaction_name)
		reaction_name = re.sub("__","_",reaction_name)
		reaction_name = str("R_")+reaction_name
		rxn2lb[reaction_name] = lb[i]

#Import upper boundary adjustments if the argument is supplied from the command line. 
if args.u is not None:
	args_u = str(args.u)
	ub_file = open(args_u)
	ub = ub_file.readlines()
	for i, item1 in enumerate(ub):
		ub[i] = item1.strip()
		ub[i] = float(ub[i])
	ub_file.close()
	for i, item1 in enumerate(cobra_model.reactions):
		reaction_name = re.sub("LPAREN","",str(item1))
		reaction_name = re.sub("RPAREN","",reaction_name)
		reaction_name = re.sub("\(","_",reaction_name)
		reaction_name = re.sub("\)","_",reaction_name)
		reaction_name = re.sub("\[","_",reaction_name)
		reaction_name = re.sub("\]","_",reaction_name)
		reaction_name = re.sub("\-","_",reaction_name)
		reaction_name = re.sub("__","_",reaction_name)
		reaction_name = str("R_")+reaction_name
		rxn2ub[reaction_name] = ub[i]

#Import gene rule adjustments if the argument is supplied from the command line. 
if args.g is not None:
	args_g = str(args.g)
	genes_genes2rxn_file = open(args_g)
	genes_genes2rxn = genes_genes2rxn_file.readlines()
	for i, item1 in enumerate(genes_genes2rxn):
		genes_genes2rxn[i] = item1.strip()
	genes_genes2rxn_file.close()
	for i, item1 in enumerate(cobra_model.reactions):
		reaction_name = re.sub("LPAREN","",str(item1))
		reaction_name = re.sub("RPAREN","",reaction_name)
		reaction_name = re.sub("\(","_",reaction_name)
		reaction_name = re.sub("\)","_",reaction_name)
		reaction_name = re.sub("\[","_",reaction_name)
		reaction_name = re.sub("\]","_",reaction_name)
		reaction_name = re.sub("\-","_",reaction_name)
		reaction_name = re.sub("__","_",reaction_name)
		reaction_name = str("R_")+reaction_name
		gene2rxn[reaction_name] = genes_genes2rxn[i]

#Create the necesssary rxn dictionaries for EXAMO.
count = 0
b_met = []
for i in cobra_model.reactions:
	count += 1
	reaction_name = re.sub("LPAREN","",i.id)
	reaction_name = re.sub("RPAREN","",reaction_name)
	reaction_name = re.sub("\(","_",reaction_name)
	reaction_name = re.sub("\)","_",reaction_name)
	reaction_name = re.sub("\[","_",reaction_name)
	reaction_name = re.sub("\]","_",reaction_name)
	reaction_name = re.sub("\-","_",reaction_name)
	reaction_name = re.sub("__","_",reaction_name)
	reaction_name = str("R_")+reaction_name
	reactants = {}
	reactants_original = {}
	products = {}
	products_original = {}
	idRs.append(reaction_name)
	if args.l is None:
		rxn2lb[reaction_name] = i.lower_bound
		lb.append(float(i.lower_bound))
	if args.u is None:
		rxn2ub[reaction_name] = i.upper_bound
		ub.append(float(i.upper_bound))
	if args.g is None:
		gene2rxn[reaction_name] = i.gene_reaction_rule
	pathway.append(i.subsystem)
	#Now need rxn
	for j in i.metabolites:
		if i.metabolites[j] < 0:
			metabolite_name = re.sub("LPAREN","",str(j.id))
			metabolite_name = re.sub("RPAREN","",metabolite_name)
			metabolite_name = re.sub("\(","_",metabolite_name)
			metabolite_name = re.sub("\)","_",metabolite_name)
			metabolite_name = re.sub("\[","_",metabolite_name)
			metabolite_name = re.sub("\]","_",metabolite_name)
			metabolite_name = re.sub("\-","_",metabolite_name)
			metabolite_name = re.sub("__","_",metabolite_name)
			reactants[str("M_")+str(metabolite_name)] = -1*i.metabolites[j]
			reactants_original[str("M_")+str(metabolite_name)] = -1*i.metabolites[j]
		if i.metabolites[j] > 0:
			metabolite_name = re.sub("LPAREN","",str(j.id))
			metabolite_name = re.sub("RPAREN","",metabolite_name)
			metabolite_name = re.sub("\(","_",metabolite_name)
			metabolite_name = re.sub("\)","_",metabolite_name)
			metabolite_name = re.sub("\[","_",metabolite_name)
			metabolite_name = re.sub("\]","_",metabolite_name)
			metabolite_name = re.sub("\-","_",metabolite_name)
			metabolite_name = re.sub("__","_",metabolite_name)
			products[str("M_")+str(metabolite_name)] = i.metabolites[j]
			products_original[str("M_")+str(metabolite_name)] = i.metabolites[j]
	if len(products) == 0:
		for j in reactants:
			extracellular = str(args.e)
			extracellular = extracellular[2:-2]
			extracellular_string = extracellular + '\Z'
			idSp_e = re.search(extracellular_string, j)
			if idSp_e is not None:
				if extracellular[-1:] == '_':
					products[j[:-2]+str("b_")] = reactants[j]
					products_original[j[:-2]+str("b_")] = reactants[j]
					b_met.append(j[:-2]+str("b_"))
				else:
					products[j[:-1]+str("b")] = reactants[j]
					products_original[j[:-1]+str("b")] = reactants[j]
					b_met.append(j[:-1]+str("b"))
	if (rxn2lb[reaction_name] < 0 and rxn2ub[reaction_name] > 0):
		reversible = True
	else:
		reversible = False
	rxns[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants, 'products': products, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}
	rxns_original[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants_original, 'products': products_original, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}  


genes = set()
for i in cobra_model.genes:
	genes.add(i.id)

#Create dictionaries to look up reactions that have metabolite mappings if the argument is supplied from the command line. 
if args.m is not None:
	args_m = str(args.m)
	md = pickle.load(open(args_m, 'rb'))
	metabolite_mappings = md['metabolite_mappings']
	met_dict_rxns = {}
	met_dict_mets = {}
	rxns_mets_to_delete = {}
	reactant_count_dict = {}
	product_count_dict = {}
	for t in rxns:
		reactant_count = 0
		secondary_check = 0
		exception_reactant_count = 0
		product_count = 0
		exception_product_count = 0
		proceed = 0
		metabolite_mappings_rxn_reactant_list = []
		metabolite_mappings_rxn_product_list = []
		for i in metabolite_mappings:
			j_count = 0
			if i in rxns[t]['reactants']:
				proceed +=1
				if proceed == 1:
					for j in metabolite_mappings[i]:
						if j in rxns[t]['products']:
							metabolite_mappings_rxn_reactant_list.append(i)
							metabolite_mappings_rxn_product_list.append(j)						
							try:
								met_dict_mets[i].append(t)
							except KeyError:
								met_dict_mets[i] = [t]
							try:
								met_dict_mets[j].append(t)
							except KeyError:
								met_dict_mets[j] = [t]
							for reactant in rxns[t]['reactants']:
								if reactant != i:
									if reactant in metabolite_mappings:
										for k in metabolite_mappings[reactant]:
											if k in rxns[t]['products']:
												metabolite_mappings_rxn_reactant_list.append(reactant)
												metabolite_mappings_rxn_product_list.append(k)
												try:
													met_dict_mets[reactant].append(t)
												except KeyError:
													met_dict_mets[reactant] = [t]
												try:
													met_dict_mets[k].append(t)
												except KeyError:
													met_dict_mets[k] = [t]
												secondary_check += 1
									else:
								 		reactant_count += 1
									#if reactant in met_exception_list:
									#	exception_reactant_count += 1
							#for product in rxns[t]['products']:
							#	if product not in metabolite_mappings[i]:
							#		if product in met_exception_list:
							#			exception_product_count += 1
							#To account for sink reactions, have several conditions that may result in a metabolite mapping exception
							if (((len(rxns[t]['reactants']) > secondary_check + 1) and (len(rxns[t]['products']) > secondary_check + 1))  or ((len(rxns[t]['reactants']) > secondary_check + 1) and (len(rxns[t]['products']) == secondary_check + 1)) or ((len(rxns[t]['reactants']) == secondary_check + 1) and (len(rxns[t]['products']) > secondary_check + 1))): 						
								met_dict_mets_list = metabolite_mappings_rxn_reactant_list + metabolite_mappings_rxn_product_list
								met_dict_rxns[t] = met_dict_mets_list
						else: 
							j_count += 1
						if j_count == len(metabolite_mappings[i]):					
							proceed = 0

	#Remove mets from reactions with metabolite mappings
	rxn_index = 0	
	for i in idRs:
		rxn_index += 1
		met_index = 0
		rxn_index_0 = []
		met_index_0 = []
		if i in met_dict_rxns:
			rxn_index_0.append(rxn_index - 1)
			for j in idSp:
				met_index += 1
				if j in met_dict_rxns[i]:
					met_index_0.append(met_index - 1)
					if j in rxns[i]['reactants']:
						del rxns[i]['reactants'][j]
					if j in rxns[i]['products']:
						del rxns[i]['products'][j]
		S[met_index_0, rxn_index_0] = 0

#Map lb and ub to rxns based on the order in which they appear
rxn_lb = {}
rxn_ub = {}
rxn_gene2rxn = {}
count = 0 
for i in idRs:
	count += 1
	count_lb = 0
	count_ub = 0
	for j in lb:
		count_lb += 1
		if count == count_lb:
			rxn_lb[i] = j 
	for j in ub:
		count_ub += 1
		if count == count_ub:
			rxn_ub[i] = j

#Remove repetitive nucleotides that share a carbon source if the nucleotide conversion argument is supplied from the command line. 
if args.n is not None:
	args_n = str(args.n)
	md = pickle.load(open(args_n, 'rb'))
	for repetitivemet in md['nucleotide_conversions']:
		mets_to_delete = []
		for i in md['nucleotide_conversions'][repetitivemet]:
			met_index = 0
			mets_to_delete.append(i)
			met_index_to_change = []
			met_index_to_delete = []
			for j in idSp:
				met_index += 1
				metabolite = str(repetitivemet)
				if j == metabolite:
					met_index_to_change.append(met_index - 1) 
			met_index = 0	
			for j in idSp:
				num = 0
				met_index += 1		
				if j == i:
					met_index_to_delete.append(met_index - 1)
					rxn_index = 0
					for t in idRs:
						rxn_index_to_change = []
						rxn_index += 1
						if i in rxns[t]['reactants']:
							rxn_index_to_change.append(rxn_index - 1)
							if metabolite not in rxns[t]['reactants']:
								rxns[t]['reactants'].update({metabolite : rxns[t]['reactants'][i]})
								num = -1*rxns[t]['reactants'][metabolite]
								S[met_index_to_change, rxn_index_to_change] = num			
							else:											
								rxns[t]['reactants'].update({metabolite : rxns[t]['reactants'][metabolite] + rxns[t]['reactants'][i]})
								num = -1*rxns[t]['reactants'][metabolite]
								S[met_index_to_change, rxn_index_to_change] = num
							S[met_index_to_delete, rxn_index_to_change] = 0
							del rxns[t]['reactants'][i]
						if i in rxns[t]['products']:
							rxn_index_to_change.append(rxn_index - 1)
							if metabolite not in rxns[t]['products']:
								rxns[t]['products'].update({metabolite : rxns[t]['products'][i]})
								S[met_index_to_change, rxn_index_to_change] = rxns[t]['products'][metabolite]
							else:											
								rxns[t]['products'].update({metabolite : rxns[t]['products'][metabolite] + rxns[t]['products'][i]})
								S[met_index_to_change, rxn_index_to_change] = rxns[t]['products'][metabolite]
							S[met_index_to_delete, rxn_index_to_change] = 0
							del rxns[t]['products'][i]

#Remove metabolites that appear as both a product and a reactant. 
rxn_count = 0
for t in idRs:
	rxn_count += 1
	rxn_index = []
	met_count = 0
	for i in idSp:
		met_index = []
		met_count += 1
		if i in rxns[t]['reactants']:
			if i in rxns[t]['products']:
			 	reactants_amount = rxns[t]['reactants'][i]
				products_amount = rxns[t]['products'][i]
				num = 0
				met_index.append(met_count - 1)
				rxn_index.append(rxn_count - 1)
				if reactants_amount > products_amount:
					rxns[t]['reactants'][i] = rxns[t]['reactants'][i] - products_amount
					rxns[t]['products'][i] = 0
					del rxns[t]['products'][i]
					num = -1*rxns[t]['reactants'][i]
					S[met_index, rxn_index] = num				
				if products_amount > reactants_amount:
					rxns[t]['products'][i] = rxns[t]['products'][i] - reactants_amount
					rxns[t]['reactants'][i] = 0
					del rxns[t]['reactants'][i]
					S[met_index, rxn_index] = rxns[t]['products'][i]				
				if reactants_amount == products_amount:
					rxns[t]['reactants'][i] = rxns[t]['reactants'][i] - products_amount
					rxns[t]['products'][i] = rxns[t]['products'][i] - reactants_amount
					del rxns[t]['reactants'][i]
					del rxns[t]['products'][i]
					S[met_index, rxn_index] = 0				
	if len(rxns[t]['reactants']) == 0:
		if len(rxns[t]['products']) == 0:
			rxns_to_remove_2.append(t)

rxn_index = 0
rxn_index_list = []
for i in idRs:
	rxn_index += 1
	if i not in rxns_to_remove_2:
		rxn_index_list.append(rxn_index-1)

for i in rxns_to_remove_2:
	del rxns[i]
	del gene2rxn[i]
	idRs.remove(i)

S = sparse.lil_matrix(sparse.csr_matrix(S)[:,rxn_index_list])

last_string_list = []
for i in idSp:
	if i[-1] == '_':
		b_name = '_b_'
		if b_name not in last_string_list:
			last_string_list.append(b_name)
		last_string = i[-3:]
		if last_string not in last_string_list:
			last_string_list.append(last_string)
	else:
		b_name = '_b'
		if b_name not in last_string_list:
			last_string_list.append(b_name)
		last_string = i[-2:]
		if last_string not in last_string_list:
			last_string_list.append(last_string)

#Import metabolite dictionary mapped to carbons if the argument is supplied from the command line. 
if args.c is not None:
	args_c = str(args.c)
	metabolite_dict_file = open(args_c)
	csvreader1 = csv.reader(metabolite_dict_file)
	metabolite_dict = {}
	metabolite_dict_csv = []
	for row in csvreader1:
		metabolite_dict_csv.append(row[0])
	for i in metabolite_dict_csv:
		i = re.split('\s|,|\t', i)
		for j in last_string_list:
			met_name = i[0]
			met_name += j
			metabolite_dict[met_name] = i[1]

	#Identify metabolites that are not in the model and cannot be imported from extracellular sources
	metabolite_dict_to_delete = []
	for i in metabolite_dict:
		if i not in idSp:
			if i not in b_met:
				metabolite_dict_to_delete.append(i)
	for i in metabolite_dict_to_delete:
		del metabolite_dict[i]

	#Idenitfy metabolites with 0 carbons. They will be removed from the model later
	met_exception_list = []
	for i in metabolite_dict:
		if int(metabolite_dict[i]) == 0:
			met_exception_list.append(i)

	#Remove metabolites without carbons
	for t in rxns:
		reactant_count = 0
		product_count = 0	
		for reactant in rxns[t]['reactants']:
			if reactant in met_exception_list:
				reactant_count += 1
		if reactant_count == len(rxns[t]['reactants']):
			if t != args_b:	
				rxns_to_remove.append(t)
		for met in met_exception_list:
			if met in rxns[t]['reactants']:
				del rxns[t]['reactants'][met]
		for product in rxns[t]['products']:
			if product in met_exception_list:
				product_count += 1
		if product_count == len(rxns[t]['products']):
			if t not in rxns_to_remove:
				if t != args_b:
					rxns_to_remove.append(t)
		for met in met_exception_list:
			if met in rxns[t]['products']:
				del rxns[t]['products'][met]
								
	rxn_index = 0
	rxn_index_list = []
	for i in idRs:
		rxn_index += 1
		if i not in rxns_to_remove:
			rxn_index_list.append(rxn_index-1)
	met_index = 0
	met_index_list = []
	for i in idSp:
		met_index += 1
		if i not in met_exception_list:
			met_index_list.append(met_index-1)

	S = sparse.lil_matrix(sparse.csr_matrix(S)[met_index_list, :])
	S = sparse.lil_matrix(sparse.csr_matrix(S)[:,rxn_index_list])

	for met in met_exception_list:
		if ((met[-2:] != '_b') and (met[-3:] != '_b_')):
			idSp.remove(met)	
	#Delete reactions and gene2rxn dictionary entries for reactions that cannot be carried out due to cofactor presence
	for t in rxns_to_remove:
		del rxns[t]
		del gene2rxn[t]
		idRs.remove(t)

	#Identify unbalanced rxns
	metabolite_dict_rxns = {}
	metabolite_dict_rxns_original = {}
	for t in rxns:
		metabolite_dict_reactants = {}
		metabolite_dict_products = {}
		for i in rxns[t]['reactants']:
			metabolite_dict_reactants[i] = metabolite_dict[i]
		for i in rxns[t]['products']:
			metabolite_dict_products[i] = metabolite_dict[i]
		rxnmetlist = {'met_reactants': metabolite_dict_reactants, 'met_products': metabolite_dict_products}
		metabolite_dict_rxns[t] = rxnmetlist 
		metabolite_dict_reactants_original = {}
		metabolite_dict_products_original = {}
		for j in rxns_original[t]['reactants']:
			metabolite_dict_reactants_original[j] = metabolite_dict[j]
		for j in rxns_original[t]['products']:
			metabolite_dict_products_original[j] = metabolite_dict[j]
		rxnmetlist = {'met_reactants': metabolite_dict_reactants_original, 'met_products': metabolite_dict_products_original}
		metabolite_dict_rxns_original[t] = rxnmetlist

	count = 0
	for t in rxns:
		reactants_count = 0
		products_count = 0
		for i in rxns[t]['reactants']:
			reactants_count += int(metabolite_dict[i])*rxns[t]['reactants'][i]
		for i in rxns[t]['products']:
			products_count += int(metabolite_dict[i])*rxns[t]['products'][i]
		for i in rxns_original[t]['reactants']:
			if i in met_exception_list:
				products_count -= int(metabolite_dict[i])*rxns_original[t]['reactants'][i]
		for i in rxns_original[t]['products']:
			if i in met_exception_list:
				reactants_count -= int(metabolite_dict[i])*rxns_original[t]['products'][i]
		if round(reactants_count,3) != round(products_count,3):
			print "\n"		
			print t
			print "carbons in original model: %s" % metabolite_dict_rxns_original[t]
			print "stoichiometry in original model: %s" % rxns_original[t]
			print "carbons in adapted model: %s" % metabolite_dict_rxns[t]
			print "stoichiometry in adapted model: %s" % rxns[t]		
			print "carbons in reactants in adapted model: %s" % reactants_count
			print "cabons in products in adapted model: %s" % products_count
			count += 1
			if t != args_b:
				unbalanced_rxns.append(t)
	#Identify unbalanced rxns
	metabolite_dict_rxns = {}
	metabolite_dict_rxns_original = {}
	for t in rxns:
		metabolite_dict_reactants = {}
		metabolite_dict_products = {}
		for i in rxns[t]['reactants']:
			metabolite_dict_reactants[i] = metabolite_dict[i]
		for i in rxns[t]['products']:
			metabolite_dict_products[i] = metabolite_dict[i]
		rxnmetlist = {'met_reactants': metabolite_dict_reactants, 'met_products': metabolite_dict_products}
		metabolite_dict_rxns[t] = rxnmetlist 
		metabolite_dict_reactants_original = {}
		metabolite_dict_products_original = {}
		for j in rxns_original[t]['reactants']:
			metabolite_dict_reactants_original[j] = metabolite_dict[j]
		for j in rxns_original[t]['products']:
			metabolite_dict_products_original[j] = metabolite_dict[j]
		rxnmetlist = {'met_reactants': metabolite_dict_reactants_original, 'met_products': metabolite_dict_products_original}
		metabolite_dict_rxns_original[t] = rxnmetlist

	#If a metabolite is not in a balanced rxn and only in an unbalanced rxn, it will be identified and later removed from the model. The rxns identified as not being balanced need to be verified whether they are due to a discrepancy in the metabolite_dict or whether the reactions are really not balanced due to an error in original model. 
	potential_unbalanced_mets = []
	for i in unbalanced_rxns:
		for k in rxns[i]['reactants']:
			if k not in potential_unbalanced_mets:		
				potential_unbalanced_mets.append(k)
		for k in rxns[i]['products']:
			if k not in potential_unbalanced_mets:
				potential_unbalanced_mets.append(k)

	balanced_rxns = list(set(idRs) - set(unbalanced_rxns))

	mets_to_remove_from_potential_unbalanced_mets = []
	for i in potential_unbalanced_mets:
		for j in balanced_rxns:
			if j != args_b:
				if i in rxns[j]['reactants']:
					if i not in mets_to_remove_from_potential_unbalanced_mets:
						mets_to_remove_from_potential_unbalanced_mets.append(i)
				if i in rxns[j]['products']:
					if i not in mets_to_remove_from_potential_unbalanced_mets:
						mets_to_remove_from_potential_unbalanced_mets.append(i)

	mets_to_remove = list(set(potential_unbalanced_mets) - set(mets_to_remove_from_potential_unbalanced_mets))

	#Remove unbalanced reactions
	rxn_index = 0
	rxn_index_list = []
	for i in idRs:
		rxn_index += 1
		if i not in unbalanced_rxns:
			rxn_index_list.append(rxn_index-1)

	for i in unbalanced_rxns:
		del rxns[i]
		del gene2rxn[i]
		idRs.remove(i)

	S = sparse.lil_matrix(sparse.csr_matrix(S)[:,rxn_index_list])

#Remove reactions that do not have any reactants or products##Consider deleting
rxns_to_delete = []
#for i in rxns:
#	if i != args_b:
#		if len(rxns[i]['reactants']) == 0:
#			if i[:4] != 'R_EX':
#				rxns_to_delete.append(i)
#		if len(rxns[i]['products']) == 0:
#			if i[:4] != 'R_EX':
#				rxns_to_delete.append(i)

#rxn_index = 0
#rxn_index_list = []
#for i in idRs:
#	rxn_index += 1
#	if i not in rxns_to_delete:
#		rxn_index_list.append(rxn_index-1)

#for i in rxns_to_delete:
#	del rxns[i]
#	del gene2rxn[i]
#	idRs.remove(i)

#S = sparse.lil_matrix(sparse.csr_matrix(S)[:,rxn_index_list])



#Remove mets that do not appear in any rxns##Consider deleting
idSp_to_remove = []
met_index = 0
met_index_list = []
for i in idSp:
	met_index += 1
	met_occurrence = 0
	for j in rxns:
		if i in rxns[j]['reactants']:
			met_occurrence += 1
		if i in rxns[j]['products']:	
			met_occurrence += 1
	if met_occurrence == 0:
		idSp_to_remove.append(i)
	else:
		met_index_list.append(met_index - 1)

S = sparse.lil_matrix(sparse.csr_matrix(S)[met_index_list,:])

for i in idSp_to_remove:
	idSp.remove(i)

#Remove lb and ub entries for reactions that were removed from the model 
lb = []
ub = []
for i in idRs:
	print i
	for rxn_mapping in rxn_lb:
		if rxn_mapping == i:		
			if i not in rxns_to_remove:
				if i not in rxns_to_remove_2:
					if i not in rxns_to_remove_3:
						if i not in unbalanced_rxns:
							if i not in rxns_to_delete:
								lb.append(rxn_lb[i])
	for rxn_mapping in rxn_ub:
		if rxn_mapping == i:		
			if i not in rxns_to_remove:
				if i not in rxns_to_remove_2:
					if i not in rxns_to_remove_3:
						if i not in unbalanced_rxns:
							if i not in rxns_to_delete:
								ub.append(rxn_ub[i])


#Convert the matrix into a coo matrix
S = coo_matrix(S)

#Print all of the reactions in the final model
#for t in rxns:
	#print "%s -> %s" % (rxns[t]['reactants'], rxns[t]['products'])

#Export the model as a pickle file
exportPickle({'idSp' : idSp, 'idRs' : idRs, 'genes' : genes, 'lb' : lb, 'ub' : ub, 'gene2rxn': gene2rxn, 'S' : S, 'rxns' : rxns }, test_model)
