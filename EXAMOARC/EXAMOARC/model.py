import re
import cobra
import scipy as sp
import numpy as np

def set_parameter(cobra_model, extracellular, args_l, args_u, args_g):
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
	if args_l:
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
	if args_l:
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
	if args_g:
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
		if args_l:
			rxn2lb[reaction_name] = i.lower_bound
			lb.append(float(i.lower_bound))
		if args_u:
			rxn2ub[reaction_name] = i.upper_bound
			ub.append(float(i.upper_bound))
		if args_g:
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



def metabolite_mapping():
	print "something"


def nucleotide_conversion():
	print "hi"


def modify():
	print "hi"


def balance_reactions():
	print "hi"
