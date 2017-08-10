import re
import cobra
import scipy as sp
import numpy as np

def set_parameter(cobra_model, extracellular, args_l, args_u, args_g)
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
	S = sp.sparse.coo_matrix(cobra_model.S)
	S = sp.sparse.lil_matrix(S)

	#Import lower boundary adjustments if the argument is supplied from the command line. 
	if args_l is not None:
		lb_file = open('%s' % args_l)
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
	if args_l is not None:
		ub_file = open('%s' % args_u)
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
	if args_g is not None:
		genes_genes2rxn_file = open('%s' % args_g)
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
	b_met = []
	for i in cobra_model.reactions:
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
		if args_l is None:
			rxn2lb[reaction_name] = i.lower_bound
			lb.append(float(i.lower_bound))
		if args_u is None:
			rxn2ub[reaction_name] = i.upper_bound
			ub.append(float(i.upper_bound))
		if args_g is None:
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

	#Create list/sets/matlab cells of genes using cobrapy and array based model for reacitons
	genes_cobra = cobra.io.mat._cell(cobra_model.genes.list_attr("id"))
	grRules = cobra.io.mat._cell(cobra_model.reactions.list_attr("gene_reaction_rule"))
	c = np.array(cobra_model.reactions.list_attr("objective_coefficient")) * 1
	subsystem = cobra.io.mat._cell(cobra_model.reactions.list_attr("subsystem"))
	metNames = cobra.io.mat._cell(cobra_model.metabolites.list_attr("name"))
	metFormulas = cobra.io.mat._cell([str(m.formula) for m in cobra_model.metabolites])
	b = np.array(cobra_model.metabolites.list_attr("_bound")) * 1.

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
	rxnGeneMat = csr_matrix((len(idRs),len(genes)))
	rxnGeneMat = lil_matrix(rxnGeneMat)
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
			rules[i] = re.sub(j,"x(%s)" % genes_list_dict[j],str(rules[i]))
			rules[i] = re.sub("or","|",str(rules[i]))
			rules[i] = re.sub("and","&",str(rules[i]))
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



def metabolite_mapping()
	if args.m is not None:
		md = pickle.load(open('data/models/%s' % args_m, 'rb'))
		####Need to figure out how to import properly
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



def nucleotide_conversion()
	if args.n is not None:
		####Need to figure out how to import properly
		md = pickle.load(open('data/models/%s' % args_n, 'rb'))
		#Remove repetitive nucleotides
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


def modify()



def balance_reactions()
