import re
import cobra
import scipy as sp
import numpy as np
import cPickle as pickle #Remove this if we are able to get rid of all pickling

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
	rxns_to_remove = []	
	rxns_to_remove_2 = []
	unbalanced_rxns = []

	idSp = []
	for i in cobra_model.metabolites:
		#metabolite_name = re.sub("LPAREN","",str(i.id))
		#metabolite_name = re.sub("RPAREN","",metabolite_name)
		#metabolite_name = re.sub("\(","_",metabolite_name)
		#metabolite_name = re.sub("\)","_",metabolite_name)
		#metabolite_name = re.sub("\[","_",metabolite_name)
		#metabolite_name = re.sub("\]","_",metabolite_name)
		#metabolite_name = re.sub("\-","_",metabolite_name)
		#metabolite_name = re.sub("__","_",metabolite_name)
		idSp.append(str("M_")+str(i.id))

	cobra_model = cobra.core.ArrayBasedModel(cobra_model)
	S = sp.sparse.coo_matrix(cobra_model.S)
	S = sp.sparse.lil_matrix(S)

	#Need to change the code to import as dictionaries
	#Import lower boundary adjustments if the argument is supplied from the command line. 
	if args_l:
		for x in args_l:
			lb.append(float(x))
		for i, item1 in enumerate(cobra_model.reactions):
			#reaction_name = re.sub("LPAREN","",str(item1))
			#reaction_name = re.sub("RPAREN","",reaction_name)
			#reaction_name = re.sub("\(","_",reaction_name)
			#reaction_name = re.sub("\)","_",reaction_name)
			#reaction_name = re.sub("\[","_",reaction_name)
			#reaction_name = re.sub("\]","_",reaction_name)
			#reaction_name = re.sub("\-","_",reaction_name)
			#reaction_name = re.sub("__","_",reaction_name)
			reaction_name = 'R_'+str(item1)
			rxn2lb[reaction_name] = lb[i]

	#Import upper boundary adjustments if the argument is supplied from the command line. 
	if args_u:
		for x in args_u:
			ub.append(float(x))
		for i, item1 in enumerate(cobra_model.reactions):
			#reaction_name = re.sub("LPAREN","",str(item1))
			#reaction_name = re.sub("RPAREN","",reaction_name)
			#reaction_name = re.sub("\(","_",reaction_name)
			#reaction_name = re.sub("\)","_",reaction_name)
			#reaction_name = re.sub("\[","_",reaction_name)
			#reaction_name = re.sub("\]","_",reaction_name)
			#reaction_name = re.sub("\-","_",reaction_name)
			#reaction_name = re.sub("__","_",reaction_name)
			reaction_name = 'R_'+str(item1)
			rxn2ub[reaction_name] = ub[i]

	#Import gene rule adjustments if the argument is supplied from the command line. 
	if args_g:
		for x in args_g:
			genes_genes2rxn.append(x)
		for i, item1 in enumerate(cobra_model.reactions):
			#reaction_name = re.sub("LPAREN","",str(item1))
			#reaction_name = re.sub("RPAREN","",reaction_name)
			#reaction_name = re.sub("\(","_",reaction_name)
			#reaction_name = re.sub("\)","_",reaction_name)
			#reaction_name = re.sub("\[","_",reaction_name)
			#reaction_name = re.sub("\]","_",reaction_name)
			#reaction_name = re.sub("\-","_",reaction_name)
			#reaction_name = re.sub("__","_",reaction_name)
			reaction_name = 'R_'+str(item1)
			gene2rxn[reaction_name] = genes_genes2rxn[i]

	#Create the necesssary rxn dictionaries for EXAMO.
	b_met = []
	for i in cobra_model.reactions:
		#reaction_name = re.sub("LPAREN","",i.id)
		#reaction_name = re.sub("RPAREN","",reaction_name)
		#reaction_name = re.sub("\(","_",reaction_name)
		#reaction_name = re.sub("\)","_",reaction_name)
		#reaction_name = re.sub("\[","_",reaction_name)
		#reaction_name = re.sub("\]","_",reaction_name)
		#reaction_name = re.sub("\-","_",reaction_name)
		#reaction_name = re.sub("__","_",reaction_name)
		reaction_name = 'R_'+i.id
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
				#metabolite_name = re.sub("LPAREN","",str(j.id))
				#metabolite_name = re.sub("RPAREN","",metabolite_name)
				#metabolite_name = re.sub("\(","_",metabolite_name)
				#metabolite_name = re.sub("\)","_",metabolite_name)
				#metabolite_name = re.sub("\[","_",metabolite_name)
				#metabolite_name = re.sub("\]","_",metabolite_name)
				#metabolite_name = re.sub("\-","_",metabolite_name)
				#metabolite_name = re.sub("__","_",metabolite_name)
				reactants['M_'+j.id] = -1*i.metabolites[j]
				reactants_original['M_'+j.id] = -1*i.metabolites[j]
			if i.metabolites[j] > 0:
				#metabolite_name = re.sub("LPAREN","",str(j.id))
				#metabolite_name = re.sub("RPAREN","",metabolite_name)
				#metabolite_name = re.sub("\(","_",metabolite_name)
				#metabolite_name = re.sub("\)","_",metabolite_name)
				#metabolite_name = re.sub("\[","_",metabolite_name)
				#metabolite_name = re.sub("\]","_",metabolite_name)
				#metabolite_name = re.sub("\-","_",metabolite_name)
				#metabolite_name = re.sub("__","_",metabolite_name)
				products['M_'+j.id] = i.metabolites[j]
				products_original['M_'+j.id] = i.metabolites[j]
		if len(products) == 0:
			for j in reactants:
				extracellular_string = args_e + '\Z'
				idSp_e = re.search(extracellular_string, j)
				if idSp_e is not None:
					if args_e[-1:] == '_':
						products[j[:-2]+'b_'] = reactants[j]
						products_original[j[:-2]+'b_'] = reactants[j]
						b_met.append(j[:-2]+'b_')
					else:
						products[j[:-1]+'b'] = reactants[j]
						products_original[j[:-1]+'b'] = reactants[j]
						b_met.append(j[:-1]+'b')
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

	#Create gene set (for EXAMO model) and list (for COBRA model)
	genes = set()
	genes_list = []
	for keys,values in list(gene2rxn.items()):
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
	rxnGeneMat = sp.sparse.csr_matrix((len(idRs),len(genes)))
	rxnGeneMat = sp.sparse.lil_matrix(rxnGeneMat)
	rxn_index = 0
	rules = grRules.copy()
	for i in range(0,len(cobra_model.reactions)):
		rxn_index += 1
		gene_index = 0
		rxn_index_0 = []
		gene_index_0 = []	
		for j in genes_list:
			gene_index += 1
			genes_search = re.search(str(j), str(grRules[i]))
			if genes_search is not None: 
				rxn_index_0.append(rxn_index - 1)
				gene_index_0.append(gene_index - 1)
			rules[i] = re.sub(j,'x(%s)' % genes_list_dict[j],str(rules[i]))
			rules[i] = re.sub('or','|',str(rules[i]))
			rules[i] = re.sub('and','&',str(rules[i]))
		rxnGeneMat[rxn_index_0, gene_index_0] = 1

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
		rxn_index_list_to_delete = []
		for i in idRs:
			rxn_index += 1
			if i not in rxns_to_remove_2:
				rxn_index_list.append(rxn_index-1)
			else:
				rxn_index_list_to_delete.append(rxn_index-1)
		grRules = np.delete(grRules, rxn_index_list_to_delete)
		rules = np.delete(rules, rxn_index_list_to_delete)
		c = np.delete(c, rxn_index_list_to_delete)
		subsystem = np.delete(subsystem, rxn_index_list_to_delete)

		for i in rxns_to_remove_2:
			del rxns[i]
			del gene2rxn[i]
			idRs.remove(i)

		S = sp.sparse.lil_matrix(sp.sparse.csr_matrix(S)[:,rxn_index_list])

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


	S = sp.sparse.lil_matrix(sp.sparse.csr_matrix(S)[met_index_list,:])

	met_index = 0
	met_index_to_delete = []
	for i in idSp:
		met_index += 1
		if i in idSp_to_remove:
			met_index_to_delete.append(met_index-1)
	metNames = np.delete(metNames, met_index_to_delete)
	metFormulas = np.delete(metFormulas, met_index_to_delete)
	b = np.delete(b, met_index_to_delete)

	for i in idSp_to_remove:
		idSp.remove(i)

	#Remove lb and ub entries for reactions that were removed from the model, and map reversibility 
	lb = []
	ub = []
	rev = []
	for i in idRs:
		for rxn_mapping in rxn_lb:
			if rxn_mapping == i:		
				if i not in rxns_to_remove:
					if i not in rxns_to_remove_2:
						if i not in unbalanced_rxns:
							lb.append(rxn_lb[i])
							ub.append(rxn_ub[i])
							if rxn_lb < 0 and rxn_ub > 0:
								rev.append(1)
							else:
								rev.append(0)

	#Create the mmodel dictionary to be able to export to other functions
	model = {'idSp' : idSp, 'idRs' : idRs, 'genes' : genes, 'lb' : lb, 'ub' : ub, 'gene2rxn': gene2rxn, 'S' : S, 'rxns' : rxns }	
	return model


def metabolite_mapping(model, args_m):
	md = pickle.load(args_m)
	metabolite_mappings = md['metabolite_mappings']
	met_dict_rxns = {}
	met_dict_mets = {}
	rxns_mets_to_delete = {}
	reactant_count_dict = {}
	product_count_dict = {}
	for t in model['rxns']:
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
			if i in model['rxns'][t]['reactants']:
				proceed +=1
				if proceed == 1:
					for j in metabolite_mappings[i]:
						if j in model['rxns'][t]['products']:
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
							for reactant in model['rxns'][t]['reactants']:
								if reactant != i:
									if reactant in metabolite_mappings:
										for k in metabolite_mappings[reactant]:
											if k in model['rxns'][t]['products']:
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
							if (((len(model['rxns'][t]['reactants']) > secondary_check + 1) and (len(model['rxns'][t]['products']) > secondary_check + 1))  or ((len(model['rxns'][t]['reactants']) > secondary_check + 1) and (len(model['rxns'][t]['products']) == secondary_check + 1)) or ((len(model['rxns'][t]['reactants']) == secondary_check + 1) and (len(model['rxns'][t]['products']) > secondary_check + 1))): 						
								met_dict_mets_list = metabolite_mappings_rxn_reactant_list + metabolite_mappings_rxn_product_list
								met_dict_rxns[t] = met_dict_mets_list
						else: 
							j_count += 1
						if j_count == len(metabolite_mappings[i]):					
							proceed = 0

	#Remove mets from reactions with metabolite mappings
	rxn_index = 0	
	for i in model['idRs']:
		rxn_index += 1
		met_index = 0
		rxn_index_0 = []
		met_index_0 = []
		if i in met_dict_rxns:
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
	return model


def nucleotide_conversion():
	print "hi"


def modify():
	print "hi"


def balance_reactions():
	print "hi"

def model_export(mdoel, model_desc):
	#Convert EXAMO data structures into COBRA compliant data structures that can be ported to MATLAB
	rxns_matlab = np.zeros((len(idRs),), dtype=np.object)
	rxns_matlab[:] = idRs
	rxns_matlab = rxns_matlab[np.newaxis].T

	mets_matlab = np.zeros((len(idSp),), dtype=np.object)
	mets_matlab[:] = idSp
	mets_matlab = mets_matlab[np.newaxis].T

	genes_matlab = np.zeros((len(genes_list),), dtype=np.object)
	genes_matlab[:] = genes_list
	genes_matlab = genes_matlab[np.newaxis].T

	ub_matlab = np.zeros((len(ub),), dtype=np.float64)
	ub_matlab[:] = ub
	ub_matlab = ub_matlab[np.newaxis].T

	lb_matlab = np.zeros((len(lb),), dtype=np.float64)
	lb_matlab[:] = lb
	lb_matlab = lb_matlab[np.newaxis].T

	rev_cobra = np.zeros((len(rev),), dtype=np.float64)
	rev_cobra[:] = rev
	rev_cobra = rev_cobra[np.newaxis].T

	grRules = grRules[np.newaxis].T

	rules = rules[np.newaxis].T

	c = c[np.newaxis].T

	subsystem = subsystem[np.newaxis].T

	metNames = metNames[np.newaxis].T

	metFormulas = metFormulas[np.newaxis].T

	b = b[np.newaxis].T

	#Convert the matrices into coo matrices
	S = sp.sparse.coo_matrix(S)

	rxnGeneMat = sp.sparse.coo_matrix(rxnGeneMat)

	model_matlab = {'rxns': rxns_matlab, 'mets': mets_matlab, 'ub': ub_matlab, 'lb': lb_matlab, 'S': S, 'grRules': grRules, 'rules': rules, 'genes': genes_matlab, 'rxnGeneMat': rxnGeneMat, 'rev': rev_cobra, 'c': c, 'subsystem': subsystem, 'metNames': metNames, 'metFormulas': metFormulas, 'b': b, 'description': model_desc}
