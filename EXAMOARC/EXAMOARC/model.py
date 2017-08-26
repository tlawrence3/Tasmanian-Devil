import re
import cobra
import scipy as sp
import numpy as np
import cPickle as pickle #Remove this if we are able to get rid of all pickling
import csv

def set_parameter(cobra_model, args_s, args_c, args_e, args_l, args_u, args_g, model_desc):
	#Import the model
	if args_s:
		cobra_model = cobra.io.read_sbml_model(cobra_model)
	if args_c:
		cobra_model = cobra.io.mat.load_matlab_model(cobra_model)

	idRs = []
	lb = []
	ub = []
	genes_genes2rxn = []
	pathway = []
	rxn2lb = {}
	rxn2ub = {}
	gene2rxn = {}
	rxns = {}
	rxns_original = {}

	idSp = []
	for i in cobra_model.metabolites:
		name = name_sub(i.id, "M_")
		idSp.append(name)

	cobra_model = cobra.core.ArrayBasedModel(cobra_model)
	S = sp.sparse.coo_matrix(cobra_model.S)
	S = sp.sparse.lil_matrix(S)

	#Import lower boundary adjustments if the argument is supplied from the command line. 
	if args_l:
		lb_dict = {}		
		csv_file = csv.reader(args_l)		
		for row in csv_file:
			name = row[0]
			name = name_sub(name, "R_")
			rxn2lb[name] = float(row[1])

	#Import upper boundary adjustments if the argument is supplied from the command line. 
	if args_u:
		ub_dict = {}		
		csv_file = csv.reader(args_u)		
		for row in csv_file:
			name = row[0]
			name = name_sub(name, "R_")
			rxn2ub[name] = float(row[1])

	#Import gene rule adjustments if the argument is supplied from the command line. 
	if args_g:
		genes_dict = {}
		csv_file = csv.reader(args_g)
		for row in csv_file:
			name = row[0]
			name = name_sub(name, "R_")
			gene2rxn[name] = row[1]

	#Create the necesssary rxn dictionaries for EXAMO.
	b_met = []	
	for i in cobra_model.reactions:
		reaction_name = name_sub(i.id, "R_")
		reactants = {}
		reactants_original = {}
		products = {}
		products_original = {}
		idRs.append(reaction_name)
		if not args_l:
			rxn2lb[reaction_name] = i.lower_bound
			lb.append(float(i.lower_bound))
		if not args_u:
			rxn2ub[reaction_name] = i.upper_bound
			ub.append(float(i.upper_bound))
		if not args_g:
			gene2rxn[reaction_name] = i.gene_reaction_rule
		pathway.append(i.subsystem)
		#Now need rxn
		for j in i.metabolites:
			if i.metabolites[j] < 0:
				metabolite_name = name_sub(j.id, "M_")
				reactants[metabolite_name] = -1*i.metabolites[j]
				reactants_original[metabolite_name] = -1*i.metabolites[j]
			if i.metabolites[j] > 0:
				metabolite_name = name_sub(j.id, "M_")
				products[metabolite_name] = i.metabolites[j]
				products_original[metabolite_name] = i.metabolites[j]
		if len(products) == 0:
			for j in reactants:
				extracellular_string = args_e + '\Z'
				idSp_e = re.search(extracellular_string, j)
				if idSp_e:
					if args_e[-1:] == '_':
						products[j[:-2]+'b_'] = reactants[j]
						products_original[j[:-2]+'b_'] = reactants[j]
						b_met.append(j[:-2]+str("b_"))
					else:
						products[j[:-1]+'b'] = reactants[j]
						products_original[j[:-1]+'b'] = reactants[j]
						b_met.append(j[:-1]+str("b"))
		if (rxn2lb[reaction_name] < 0 and rxn2ub[reaction_name] > 0):
			reversible = True
		else:
			reversible = False
		rxns[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants, 'products': products, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}
		rxns_original[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants_original, 'products': products_original, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}  

	#Create list/sets/matlab cells of genes using cobrapy and array based model for reacitons
	genes_cobra = cobra.io.mat._cell(cobra_model.genes.list_attr('id'))
	grRules = cobra.io.mat._cell(cobra_model.reactions.list_attr('gene_reaction_rule'))
	c = np.array(cobra_model.reactions.list_attr('objective_coefficient')) * 1
	subsystem = cobra.io.mat._cell(cobra_model.reactions.list_attr('subsystem'))
	metNames = cobra.io.mat._cell(cobra_model.metabolites.list_attr('name'))
	metFormulas = cobra.io.mat._cell([str(m.formula) for m in cobra_model.metabolites])
	b = np.array(cobra_model.metabolites.list_attr('_bound')) * 1.

	#Remove mets that do not appear in any rxns. This is a model cleanup step and may be needed for linear optimizations not to fail for EXAMO 
	#idSp_to_remove = []
	#met_index = 0
	#met_index_list = []
	#for i in idSp:
	#	met_index += 1
	#	met_occurrence = 0
	#	for j in rxns:
	#		if i in rxns[j]['reactants']:
	#			met_occurrence += 1
	#		if i in rxns[j]['products']:	
	#			met_occurrence += 1
	#	if met_occurrence == 0:
	#		idSp_to_remove.append(i)
	#	else:
	#		met_index_list.append(met_index - 1)


	#S = sp.sparse.lil_matrix(sp.sparse.csr_matrix(S)[met_index_list,:])

	#met_index = 0
	#met_index_to_delete = []
	#for i in idSp:
	#	met_index += 1
	#	if i in idSp_to_remove:
	#		met_index_to_delete.append(met_index-1)
	#metNames = np.delete(metNames, met_index_to_delete)
	#metFormulas = np.delete(metFormulas, met_index_to_delete)
	#b = np.delete(b, met_index_to_delete)

	#for i in idSp_to_remove:
	#	idSp.remove(i)
	
	#Create the model and cobra specific objects dictionaries to be able to export to other functions
	model = {'idSp' : idSp, 'idRs' : idRs, 'lb' : lb, 'ub' : ub, 'gene2rxn': gene2rxn, 'S' : S, 'rxns' : rxns }
	cobra_specific_objects = {'grRules': grRules, 'c': c, 'subsystem': subsystem, 'metNames': metNames, 'metFormulas': metFormulas, 'b': b}	
	return model, cobra_specific_objects, b_met, rxns_original


def name_sub(string, prepend):
	name = re.sub("LPAREN","",string)
	name = re.sub("RPAREN","",name)
	name = re.sub("\(","_",name)
	name = re.sub("\)","_",name)
	name = re.sub("\[","_",name)
	name = re.sub("\]","_",name)
	name = re.sub("\-","_",name)
	name = re.sub("__","_",name)
	if name[:2] != prepend:
		name = prepend+name
	return name


def modify(model, cobra_specific_objects, args_a):
	modifications = {}
	csv_file = csv.reader(args_a)
	for row in csv_file:
		rxn = name_sub(row[0], "R_")
		met = name_sub(row[1], "M_")
		stoich = float(row[2])
		if rxn not in modifications:
			modifications[rxn] = {}
		modifications[rxn][met] = stoich
	rxn_index = 0
	for t in model['idRs']:
		rxn_index += 1
		met_index = 0
		rxn_index_change = []
		met_index_0 = []
		if t in modifications:
			rxn_index_change = rxn_index - 1
			for j in model['idSp']:
				met_index += 1
				if j in modifications[t]:
					met_index_change = met_index - 1
					model['S'][met_index_change, rxn_index_change] = modifications[t][j]
					if modifications[t][j] < 0:
						model['rxns'][t]['reactants'][j] = -1*modifications[t][j]
						if ((j in model['rxns'][t]['products']) or (modifications[t][j] == 0 and j in model['rxns'][t]['reactants'])):
							del model['rxns'][t]['products'][j]
					if modifications[t][j] > 0:
						model['rxns'][t]['products'][j] = modifications[t][j]
						if ((j in model['rxns'][t]['reactants']) or (modifications[t][j] == 0 and j in model['rxns'][t]['products'])):
							del model['rxns'][t]['reactants'][j]
					if (modifications[t][j] == 0 and j in model['rxns'][t]['reactants']):
						del model['rxns'][t]['reactants'][j]
					if (modifications[t][j] == 0 and j  in model['rxns'][t]['products']):
						del model['rxns'][t]['products'][j]
					
	return model, cobra_specific_objects


def metabolite_mapping(model, cobra_specific_objects, args_m):
	#Adjust for reactions that have metbaolite mappings
	if args_m:
		metabolite_mappings = {}		
		csv_file = csv.reader(args_m)
		for row in csv_file:
			met1 = name_sub(row[0], "M_")
			met2 = name_sub(row[1], "M_")
			metabolite_mappings[met1] = met2				
		met_rxns = []
		for t in model['rxns']:
			reactant_count = 0
			secondary_check = 0
			product_count = 0
			for i in metabolite_mappings.keys():
				if ((i in model['rxns'][t]['reactants']) and (metabolite_mappings[i] in model['rxns'][t]['products'])):			
					reactant_count += 1
					product_count += 1
			if ((reactant_count == len(model['rxns'][t]['reactants'])) and (product_count == len(model['rxns'][t]['produts']))):
				met_rxns.append(t)						
		#Remove mets from reactions with metabolite mappings
		rxn_index = 0
		for i in model['idRs']:
			rxn_index += 1
			met_index = 0
			rxn_index_0 = []
			met_index_0 = []
			if i in met_rxns:
				rxn_index_0.append(rxn_index - 1)
				for j in model['idSp']:
					met_index += 1
					if j in met_dict_rxns[i]:
						met_index_0.append(met_index - 1)
						if j in model['rxns'][i]['reactants']:
							del model['rxns'][i]['reactants'][j]
						if j in model['rxns'][i]['products']:
							del model['rxns'][i]['products'][j]
			model['S'][met_index_0, rxn_index_0] = 0
	return model, cobra_specific_objects


def nucleotide_conversion(model, cobra_specific_objects, args_n):
	#Remove repetitive nucleotides that share a carbon source 
	if args_n:
		nucleotide_conversions = {}		
		csv_file = csv.reader(args_n)
		for row in csv_file:
			met1 = name_sub(row[0], "M_")
			nucleotide_conversions[met1] = []
			for i in range(1,len(row)):
				met = name_sub(row[i], "M_")
				if met != "":
					nucleotide_conversions[met1].append(met)
		for repetitivemet in nucleotide_conversions:
			mets_to_delete = []
			for i in nucleotide_conversions[repetitivemet]:
				met_index = 0
				mets_to_delete.append(i)
				met_index_to_change = []
				met_index_to_delete = []
				for j in model['idSp']:
					met_index += 1
					metabolite = repetitivemet
					if j == metabolite:
						met_index_to_change.append(met_index - 1) 
				met_index = 0	
				for j in model['idSp']:
					num = 0
					met_index += 1		
					if j == i:
						met_index_to_delete.append(met_index - 1)
						rxn_index = 0
						for t in model['idRs']:
							rxn_index_to_change = []
							rxn_index += 1
							if i in model['rxns'][t]['reactants']:
								rxn_index_to_change.append(rxn_index - 1)
								if metabolite not in model['rxns'][t]['reactants']:
									model['rxns'][t]['reactants'].update({metabolite : model['rxns'][t]['reactants'][i]})
									num = -1*model['rxns'][t]['reactants'][metabolite]
									model['S'][met_index_to_change, rxn_index_to_change] = num			
								else:											
									model['rxns'][t]['reactants'].update({metabolite : model['rxns'][t]['reactants'][metabolite] + model['rxns'][t]['reactants'][i]})
									num = -1*model['rxns'][t]['reactants'][metabolite]
									model['S'][met_index_to_change, rxn_index_to_change] = num
								model['S'][met_index_to_delete, rxn_index_to_change] = 0
								del model['rxns'][t]['reactants'][i]
							if i in model['rxns'][t]['products']:
								rxn_index_to_change.append(rxn_index - 1)
								if metabolite not in model['rxns'][t]['products']:
									model['rxns'][t]['products'].update({metabolite : model['rxns'][t]['products'][i]})
									model['S'][met_index_to_change, rxn_index_to_change] = model['rxns'][t]['products'][metabolite]
								else:											
									model['rxns'][t]['products'].update({metabolite : model['rxns'][t]['products'][metabolite] + model['rxns'][t]['products'][i]})
									model['S'][met_index_to_change, rxn_index_to_change] = model['rxns'][t]['products'][metabolite]
								model['S'][met_index_to_delete, rxn_index_to_change] = 0
								del model['rxns'][t]['products'][i]

	#Remove metabolites that appear as both a product and a reactant. 
	rxns_to_remove = []	
	rxn_count = 0
	for t in model['idRs']:
		rxn_count += 1
		rxn_index = []
		met_count = 0
		for i in model['idSp']:
			met_index = []
			met_count += 1
			if i in model['rxns'][t]['reactants']:
				if i in model['rxns'][t]['products']:
				 	reactants_amount = model['rxns'][t]['reactants'][i]
					products_amount = model['rxns'][t]['products'][i]
					num = 0
					met_index.append(met_count - 1)
					rxn_index.append(rxn_count - 1)
					if reactants_amount > products_amount:
						model['rxns'][t]['reactants'][i] = model['rxns'][t]['reactants'][i] - products_amount
						model['rxns'][t]['products'][i] = 0
						del model['rxns'][t]['products'][i]
						num = -1*model['rxns'][t]['reactants'][i]
						model['S'][met_index, rxn_index] = num				
					if products_amount > reactants_amount:
						model['rxns'][t]['products'][i] = model['rxns'][t]['products'][i] - reactants_amount
						model['rxns'][t]['reactants'][i] = 0
						del model['rxns'][t]['reactants'][i]
						model['S'][met_index, rxn_index] = model['rxns'][t]['products'][i]				
					if reactants_amount == products_amount:
						model['rxns'][t]['reactants'][i] = model['rxns'][t]['reactants'][i] - products_amount
						model['rxns'][t]['products'][i] = model['rxns'][t]['products'][i] - reactants_amount
						del model['rxns'][t]['reactants'][i]
						del model['rxns'][t]['products'][i]
						model['S'][met_index, rxn_index] = 0				
		if len(model['rxns'][t]['reactants']) == 0:
			if len(model['rxns'][t]['products']) == 0:
				rxns_to_remove.append(t)

	rxn_index = 0
	rxn_index_list = []
	rxn_index_list_to_delete = []
	for i in model['idRs']:
		rxn_index += 1
		if i not in rxns_to_remove:
			rxn_index_list.append(rxn_index-1)
		else:
			rxn_index_list_to_delete.append(rxn_index-1)
	cobra_specific_objects['grRules'] = np.delete(cobra_specific_objects['grRules'], rxn_index_list_to_delete)
	cobra_specific_objects['c'] = np.delete(cobra_specific_objects['c'], rxn_index_list_to_delete)
	cobra_specific_objects['subsystem'] = np.delete(cobra_specific_objects['subsystem'], rxn_index_list_to_delete)

	for t in rxns_to_remove:
		del model['rxns'][t]
		del model['gene2rxn'][t]
		model['idRs'].remove(t)

	model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[:,rxn_index_list])


	return model, cobra_specific_objects


def balance_reactions(model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original, args_biomass, args_c2m, args_z, args_b):
	#Import metabolite dictionary mapped to carbons if the argument is supplied from the command line. 
	if args_c2m:
		csvreader1 = csv.reader(args_c2m)
		metabolite_dict = {}
		metabolite_dict_csv = []
		for row in csvreader1:
			metabolite_dict_csv.append(row[0])

		#Account for possibility of every metabolite coming from any compartment, including transport into extracellular compartment
		last_string_list = []
		for i in model['idSp']:
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

		for i in metabolite_dict_csv:
			i = re.split('\s|,|\t', i)
			for j in last_string_list:
				met_name = i[0]
				met_name += j
				metabolite_dict[met_name] = i[1]
		
		#Remove metabolites in the metabolite dictionary that are not in the model
		metabolite_dict_to_delete = []
		for i in metabolite_dict:
			if i not in model['idSp']:
				if i not in mets_to_extracellular_comp:
					metabolite_dict_to_delete.append(i)
		for i in metabolite_dict_to_delete:
			del metabolite_dict[i]	
		
		#Idenitfy metabolites with 0 carbons. They will be removed from the model	
		if args_z:
			met_exception_list = []
			rxns_to_remove = []	
			for i in metabolite_dict:
				if int(metabolite_dict[i]) == 0:
					met_exception_list.append(i)
			#Remove reactions that only contain metabolites with 0 carbons for either reactants or products			
			for t in model['rxns']:
				reactant_count = 0
				product_count = 0	
				for reactant in model['rxns'][t]['reactants']:
					if reactant in met_exception_list:
						reactant_count += 1
				if reactant_count == len(model['rxns'][t]['reactants']):
					if t != args_biomass:	
						rxns_to_remove.append(t)
				for met in met_exception_list:
					if met in model['rxns'][t]['reactants']:
						del model['rxns'][t]['reactants'][met]
				for product in model['rxns'][t]['products']:
					if product in met_exception_list:
						product_count += 1
				if product_count == len(model['rxns'][t]['products']):
					if t not in rxns_to_remove:
						if t != args_biomass:
							rxns_to_remove.append(t)
				for met in met_exception_list:
					if met in model['rxns'][t]['products']:
						del model['rxns'][t]['products'][met]
								
			rxn_index = 0
			rxn_index_list = []
			rxn_index_list_to_delete = []
			for i in model['idRs']:
				rxn_index += 1
				if i not in rxns_to_remove:
					rxn_index_list.append(rxn_index-1)
				else:
					rxn_index_list_to_delete.append(rxn_index-1)
			cobra_specific_objects['grRules'] = np.delete(cobra_specific_objects['grRules'], rxn_index_list_to_delete)
			cobra_specific_objects['c'] = np.delete(cobra_specific_objects['c'], rxn_index_list_to_delete)
			cobra_specific_objects['subsystem'] = np.delete(cobra_specific_objects['subsystem'], rxn_index_list_to_delete)

			met_index = 0
			met_index_list = []
			met_index_list_to_delete = []
			for i in model['idSp']:
				met_index += 1
				if i not in met_exception_list:
					met_index_list.append(met_index-1)
				else:
					met_index_list_to_delete.append(met_index-1)
			cobra_specific_objects['metNames'] = np.delete(cobra_specific_objects['metNames'], met_index_list_to_delete)
			cobra_specific_objects['metFormulas'] = np.delete(cobra_specific_objects['metFormulas'], met_index_list_to_delete)
			cobra_specific_objects['b'] = np.delete(cobra_specific_objects['b'], met_index_list_to_delete)

			model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[met_index_list, :])
			model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[:,rxn_index_list])

			for met in met_exception_list:
				if ((met[-2:] != '_b') and (met[-3:] != '_b_')):
					model['idSp'].remove(met)
	
			#Delete reactions and gene2rxn dictionary entries for reactions that cannot be carried out due to cofactor presence
			for t in rxns_to_remove:
				del model['rxns'][t]
				del model['gene2rxn'][t]
				model['idRs'].remove(t)

		#Identify unbalanced rxns
		#Need to adapt to print to a file
		metabolite_dict_rxns = {}
		metabolite_dict_rxns_original = {}
		unbalanced_rxns = []
		for t in model['rxns']:
			metabolite_dict_reactants = {}
			metabolite_dict_products = {}
			for i in model['rxns'][t]['reactants']:
				metabolite_dict_reactants[i] = metabolite_dict[i]
			for i in model['rxns'][t]['products']:
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
		for t in model['rxns']:
			reactants_count = 0
			products_count = 0
			try:
				if met_exception_list:
					pass
			except:
				met_exception_list = []
			for i in model['rxns'][t]['reactants']:
				reactants_count += int(metabolite_dict[i])*model['rxns'][t]['reactants'][i]
			for i in model['rxns'][t]['products']:
				products_count += int(metabolite_dict[i])*model['rxns'][t]['products'][i]
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
				print "stoichiometry in adapted model: %s" % model['rxns'][t]		
				print "carbons in reactants in adapted model: %s" % reactants_count
				print "cabons in products in adapted model: %s" % products_count
				count += 1
				if t != args_biomass:
					unbalanced_rxns.append(t)
		#The rxns identified as not being balanced need to be verified whether they are due to a discrepancy in the metabolite_dict or whether the reactions are really not balanced due to an error in original model. 
		#If verified, then the rxns can be deleted, but other adjustments may need to be made to the model.
		#If a metabolite is not in a balanced rxn and only in an unbalanced rxn, it will be removed from the model.
		if args_b: 
			potential_unbalanced_mets = []
			for i in unbalanced_rxns:
				for k in model['rxns'][i]['reactants']:
					if k not in potential_unbalanced_mets:		
						potential_unbalanced_mets.append(k)
				for k in model['rxns'][i]['products']:
					if k not in potential_unbalanced_mets:
						potential_unbalanced_mets.append(k)

			balanced_rxns = list(set(model['idRs']) - set(unbalanced_rxns))

			mets_to_remove_from_potential_unbalanced_mets = []
			for i in potential_unbalanced_mets:
				for j in balanced_rxns:
					if j != args_biomass:
						if i in model['rxns'][j]['reactants']:
							if i not in mets_to_remove_from_potential_unbalanced_mets:
								mets_to_remove_from_potential_unbalanced_mets.append(i)
						if i in model['rxns'][j]['products']:
							if i not in mets_to_remove_from_potential_unbalanced_mets:
								mets_to_remove_from_potential_unbalanced_mets.append(i)

			mets_to_remove = list(set(potential_unbalanced_mets) - set(mets_to_remove_from_potential_unbalanced_mets))

			met_index = 0
			met_index_list = []
			met_index_list_to_delete = []
			for i in model['idSp']:
				met_index += 1
				if i not in mets_to_remove:
					met_index_list.append(met_index-1)
				else:
					met_index_list_to_delete.append(met_index-1)
			cobra_specific_objects['metNames'] = np.delete(cobra_specific_objects['metNames'], met_index_list_to_delete)
			cobra_specific_objects['metFormulas'] = np.delete(cobra_specific_objects['metFormulas'], met_index_list_to_delete)
			cobra_specific_objects['b'] = np.delete(cobra_specific_objects['b'], met_index_list_to_delete)

			model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[met_index_list, :])

			for met in mets_to_remove:
				if ((met[-2:] != '_b') and (met[-3:] != '_b_')):
					model['idSp'].remove(met)

			#Remove unbalanced reactions
			rxn_index = 0
			rxn_index_list = []
			rxn_index_list_to_delete = []
			for i in model['idRs']:
				rxn_index += 1
				if i not in unbalanced_rxns:
					rxn_index_list.append(rxn_index-1)
				else:
					rxn_index_list_to_delete.append(rxn_index-1)
			cobra_specific_objects['grRules'] = np.delete(cobra_specific_objects['grRules'], rxn_index_list_to_delete)
			cobra_specific_objects['c'] = np.delete(cobra_specific_objects['c'], rxn_index_list_to_delete)
			cobra_specific_objects['subsystem'] = np.delete(cobra_specific_objects['subsystem'], rxn_index_list_to_delete)

			for i in unbalanced_rxns:
				del model['rxns'][i]
				del model['gene2rxn'][i]
				model['idRs'].remove(i)

			model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[:,rxn_index_list])
	return model, cobra_specific_objects


def metabolite_cleanup(model, cobra_specific_objects):
	#Remove mets that do not appear in any rxns
	idSp_to_remove = []
	met_index = 0
	met_index_list = []
	for i in model['idSp']:
		met_index += 1
		met_occurrence = 0
		for t in model['rxns']:
			if i in model['rxns'][t]['reactants']:
				met_occurrence += 1
			if i in model['rxns'][t]['products']:	
				met_occurrence += 1
		if met_occurrence == 0:
			idSp_to_remove.append(i)
		else:
			met_index_list.append(met_index - 1)

	model['S'] = sp.sparse.lil_matrix(sp.sparse.csr_matrix(model['S'])[met_index_list,:])

	met_index = 0
	met_index_to_delete = []
	for i in model['idSp']:
		met_index += 1
		if i in idSp_to_remove:
			met_index_to_delete.append(met_index-1)
	cobra_specific_objects['metNames'] = np.delete(cobra_specific_objects['metNames'], met_index_to_delete)
	cobra_specific_objects['metFormulas'] = np.delete(cobra_specific_objects['metFormulas'], met_index_to_delete)
	cobra_specific_objects['b'] = np.delete(cobra_specific_objects['b'], met_index_to_delete)

	for i in idSp_to_remove:
		model['idSp'].remove(i)
	return model, cobra_specific_objects


def model_export(model, cobra_specific_objects, model_desc):
	#Convert EXAMO data structures into COBRA compliant data structures that can be ported to MATLAB
	rxns_matlab = np.zeros((len(model['idRs']),), dtype=np.object)
	rxns_matlab[:] = model['idRs']
	rxns_matlab = rxns_matlab[np.newaxis].T

	mets_matlab = np.zeros((len(model['idSp']),), dtype=np.object)
	mets_matlab[:] = model['idSp']
	mets_matlab = mets_matlab[np.newaxis].T

	#Create gene set (for EXAMO model) and list (for COBRA model)
	genes = set()
	genes_list = []
	for t in model['idRs']:
		values = model['rxns'][t]['genes']
	        values_split = values.split('or')
	        for count, j in enumerate(values_split):
			j = j.split('and')
			genelist = []
			for k in j:
				k = k.translate(None, ' ()')
				if not k:
					continue
				else:
					genelist.append(k)
					if k not in genes_list:
						genes.add(k)
						genes_list.append(k)

	#Create dictionary for gene occurrence to replace genes with their index in grRules for creation of rules object
	genes_list_dict = {}
	for i, j in enumerate(genes_list):
		genes_list_dict[j] = i

	#Create reaction gene matrix and rules object
	rxnGeneMat = sp.sparse.csr_matrix((len(model['idRs']),len(genes)))
	rxnGeneMat = sp.sparse.lil_matrix(rxnGeneMat)
	rxn_index = 0
	rules = cobra_specific_objects['grRules'].copy()
	for i in range(0,len(model['idRs'])):
		rxn_index += 1
		gene_index = 0
		rxn_index_0 = []
		gene_index_0 = []	
		for j in genes_list:
			gene_index += 1
			genes_search = re.search(str(j), str(cobra_specific_objects['grRules'][i]))
			if genes_search is not None: 
				rxn_index_0.append(rxn_index - 1)
				gene_index_0.append(gene_index - 1)
			rules[i] = re.sub(j,'x(%s)' % genes_list_dict[j],str(rules[i]))
			rules[i] = re.sub('or','|',str(rules[i]))
			rules[i] = re.sub('and','&',str(rules[i]))
		rxnGeneMat[rxn_index_0, gene_index_0] = 1

	genes_matlab = np.zeros((len(genes_list),), dtype=np.object)
	genes_matlab[:] = genes_list
	genes_matlab = genes_matlab[np.newaxis].T

	ub = []
	for t in model['idRs']:
		ub.append(model['rxns'][t]['ub'])
	ub_matlab = np.zeros((len(ub),), dtype=np.float64)
	ub_matlab[:] = ub
	ub_matlab = ub_matlab[np.newaxis].T

	lb = []
	for t in model['idRs']:
		lb.append(model['rxns'][t]['lb'])	
	lb_matlab = np.zeros((len(lb),), dtype=np.float64)
	lb_matlab[:] = lb
	lb_matlab = lb_matlab[np.newaxis].T

	rev = []
	for t in model['idRs']:
		rev.append(model['rxns'][t]['reversible'])
	rev_cobra = np.zeros((len(rev),), dtype=np.float64)
	rev_cobra[:] = rev
	rev_cobra = rev_cobra[np.newaxis].T

	grRules = cobra_specific_objects['grRules'][np.newaxis].T

	rules = rules[np.newaxis].T

	c = cobra_specific_objects['c'][np.newaxis].T

	subsystem = cobra_specific_objects['subsystem'][np.newaxis].T

	metNames = cobra_specific_objects['metNames'][np.newaxis].T

	metFormulas = cobra_specific_objects['metFormulas'][np.newaxis].T

	b = cobra_specific_objects['b'][np.newaxis].T

	#Convert the matrices into coo matrices
	S = sp.sparse.coo_matrix(model['S'])

	rxnGeneMat = sp.sparse.coo_matrix(rxnGeneMat)

	model_matlab = {'rxns': rxns_matlab, 'mets': mets_matlab, 'ub': ub_matlab, 'lb': lb_matlab, 'S': S, 'grRules': grRules, 'rules': rules, 'genes': genes_matlab, 'rxnGeneMat': rxnGeneMat, 'rev': rev_cobra, 'c': c, 'subsystem': subsystem, 'metNames': metNames, 'metFormulas': metFormulas, 'b': b, 'description': model_desc}
	#model_matlab = {'rxns': rxns_matlab, 'mets': mets_matlab, 'ub': ub_matlab, 'lb': lb_matlab, 'S': S, 'grRules': grRules, 'rules': rules, 'genes': genes_matlab, 'rev': rev_cobra, 'c': c, 'subsystem': subsystem, 'metNames': metNames, 'metFormulas': metFormulas, 'b': b, 'description': model_desc}
	#Need to add genes to model and export model for EXAMO
	model['genes'] = genes
	sp.io.savemat('%s' % model_desc[:-4], model_matlab)
