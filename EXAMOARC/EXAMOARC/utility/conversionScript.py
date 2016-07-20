from __future__ import print_function
from cobra.io import read_sbml_model
from cobra.core import ArrayBasedModel, Model
from cobra.io.mat import _cell
#from cobra.mlab.mlab import python_list_to_matlab_cell, scipy_sparse_to_mlab_sparse
from scipy.sparse import coo_matrix
import scipy.io
import cobra

import os
import sys
import time
from numpy import array, delete, transpose, shape, object as np_object
import numpy as np
from scipy.io import loadmat
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
import re
import operator
import collections

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

def cobra_to_examo(model_file, c, s):
    if c:
        return cobra.io.mat.load_matlab_model(model_file)
    if s:
        return read_sbml_model(model_file)

def model_cleanup(cobra_model):
    for i in cobra_model.metabolites:
        metabolite_name = re.sub("LPAREN","",str(i.id))
        metabolite_name = re.sub("RPAREN","",metabolite_name)
        metabolite_name = re.sub("\(","_",metabolite_name)
        metabolite_name = re.sub("\)","_",metabolite_name)
        metabolite_name = re.sub("\[","_",metabolite_name)
        metabolite_name = re.sub("\]","_",metabolite_name)
        metabolite_name = re.sub("\-","_",metabolite_name)
        metabolite_name = re.sub("__","_",metabolite_name)
    for i in cobra_model.reactions:
        reaction_name = re.sub("LPAREN","",str(item1))
        reaction_name = re.sub("RPAREN","",reaction_name)
        reaction_name = re.sub("\(","_",reaction_name)
        reaction_name = re.sub("\)","_",reaction_name)
        reaction_name = re.sub("\[","_",reaction_name)
        reaction_name = re.sub("\]","_",reaction_name)
        reaction_name = re.sub("\-","_",reaction_name)
        reaction_name = re.sub("__","_",reaction_name)
        reaction_name = str("R_")+reaction_name


#Import lower boundary adjustments if the argument is supplied from the command line. 
def lowbound_set(lb_file, model):
    for i, item1 in enumerate(lb_file.readline()):
        lb[i] = item1.strip()
        lb[i] = float(lb[i])
    for i, item1 in enumerate(model.reactions):
        model.reactions.item1.lower_bound  = lb[i]
    return model

#Import upper boundary adjustments if the argument is supplied from the command line. 
def upperbound_set(ub_file, model):
    for i, item1 in enumerate(ub_file.readlines()):
        ub[i] = item1.strip()
        ub[i] = float(ub[i])
    for i, item1 in enumerate(model.reactions):
        model.reactions.item1.upper_bound = ub[i]
    return model
    
#Import gene rule adjustments if the argument is supplied from the command line. 
def gene2rxn_set(gene2rxn_file, model):
    for i, item1 in enumerate(genes2rxn_file.readlines()):
        genes_genes2rxn[i] = item1.strip()
    for i, item1 in enumerate(model.reactions):
        model.reactions.item1.gene_reaction_rule = genes_genes2rxn[i]
    return model

#cobra_model = cobra.core.ArrayBasedModel(cobra_model)
#S = coo_matrix(cobra_model.S)
#S = lil_matrix(S)
##Create the necesssary rxn dictionaries for EXAMO.
#count = 0
#b_met = []
#for i in cobra_model.reactions:
#    count += 1
#    reactants = {}
#    reactants_original = {}
#    products = {}
#    products_original = {}
#    idRs.append(reaction_name)
#    if args.l is None:
#        rxn2lb[reaction_name] = i.lower_bound
#        lb.append(float(i.lower_bound))
#    if args.u is None:
#        rxn2ub[reaction_name] = i.upper_bound
#        ub.append(float(i.upper_bound))
#    if args.g is None:
#        gene2rxn[reaction_name] = i.gene_reaction_rule
#    pathway.append(i.subsystem)
#    #Now need rxn
#    for j in i.metabolites:
#        if i.metabolites[j] < 0:
#            metabolite_name = re.sub("LPAREN","",str(j.id))
#            metabolite_name = re.sub("RPAREN","",metabolite_name)
#            metabolite_name = re.sub("\(","_",metabolite_name)
#            metabolite_name = re.sub("\)","_",metabolite_name)
#            metabolite_name = re.sub("\[","_",metabolite_name)
#            metabolite_name = re.sub("\]","_",metabolite_name)
#            metabolite_name = re.sub("\-","_",metabolite_name)
#            metabolite_name = re.sub("__","_",metabolite_name)
#            reactants[str("M_")+str(metabolite_name)] = -1*i.metabolites[j]
#            reactants_original[str("M_")+str(metabolite_name)] = -1*i.metabolites[j]
#        if i.metabolites[j] > 0:
#            metabolite_name = re.sub("LPAREN","",str(j.id))
#            metabolite_name = re.sub("RPAREN","",metabolite_name)
#            metabolite_name = re.sub("\(","_",metabolite_name)
#            metabolite_name = re.sub("\)","_",metabolite_name)
#            metabolite_name = re.sub("\[","_",metabolite_name)
#            metabolite_name = re.sub("\]","_",metabolite_name)
#            metabolite_name = re.sub("\-","_",metabolite_name)
#            metabolite_name = re.sub("__","_",metabolite_name)
#            products[str("M_")+str(metabolite_name)] = i.metabolites[j]
#            products_original[str("M_")+str(metabolite_name)] = i.metabolites[j]
#    if len(products) == 0:
#        for j in reactants:
#            extracellular = str(args.e)
#            extracellular = extracellular[2:-2]
#            extracellular_string = extracellular + '\Z'
#            idSp_e = re.search(extracellular_string, j)
#            if idSp_e is not None:
#                if extracellular[-1:] == '_':
#                    products[j[:-2]+str("b_")] = reactants[j]
#                    products_original[j[:-2]+str("b_")] = reactants[j]
#                    b_met.append(j[:-2]+str("b_"))
#                else:
#                    products[j[:-1]+str("b")] = reactants[j]
#                    products_original[j[:-1]+str("b")] = reactants[j]
#                    b_met.append(j[:-1]+str("b"))
#    if (rxn2lb[reaction_name] < 0 and rxn2ub[reaction_name] > 0):
#        reversible = True
#    else:
#        reversible = False
#    rxns[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants, 'products': products, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}
#    rxns_original[reaction_name] = {'name': i.name, 'id': reaction_name, 'reactants': reactants_original, 'products': products_original, 'reversible': reversible, 'genes': gene2rxn[reaction_name], 'lb': rxn2lb[reaction_name], 'ub': rxn2ub[reaction_name], 'pathway': i.subsystem}  
#
##Create list/sets/matlab cells of genes using cobrapy and array based model for reacitons
#genes_cobra = _cell(cobra_model.genes.list_attr("id"))
#grRules = _cell(cobra_model.reactions.list_attr("gene_reaction_rule"))
#c = array(cobra_model.reactions.list_attr("objective_coefficient")) * 1
#subsystem = _cell(cobra_model.reactions.list_attr("subsystem"))
#metNames = _cell(cobra_model.metabolites.list_attr("name"))
#metFormulas = _cell([str(m.formula) for m in cobra_model.metabolites])
#b = array(cobra_model.metabolites.list_attr("_bound")) * 1.
#
##Create gene set (for EXAMO model) and list (for COBRA model)
#genes = set()
#genes_list = []
#for keys,values in list(gene2rxn.items()):
#        values_split = values.split('or')
#        for count, j in enumerate(values_split):
#            j = j.split('and')
#            genelist = []
#            for k in j:
#                k = k.translate(None, ' ()')
#                if not k:
#                    continue
#                else:
#                    genelist.append(k)
#            if k not in genes_list:
#                genes.add(k)
#                genes_list.append(k)
#
##Create dictionary for gene occurrence to replace genes with their index in grRules for creation of rules object
#genes_list_dict = {}
#for i, j in enumerate(genes_list):
#    genes_list_dict[j] = i
#
##Create reaction gene matrix and rules object
#rxnGeneMat = csr_matrix((len(idRs),len(genes)))
#rxnGeneMat = lil_matrix(rxnGeneMat)
#rxn_index = 0
#rules = grRules.copy()
#for i in range(0,len(cobra_model.reactions)):
#    rxn_index += 1
#    gene_index = 0
#    rxn_index_0 = []
#    gene_index_0 = []    
#    for j in genes_list:
#        gene_index += 1
#        genes_search = re.search(str(j), str(grRules[i]))
#        if genes_search is not None: 
#            rxn_index_0.append(rxn_index - 1)
#            gene_index_0.append(gene_index - 1)
#        rules[i] = re.sub(j,"x(%s)" % genes_list_dict[j],str(rules[i]))
#        rules[i] = re.sub("or","|",str(rules[i]))
#        rules[i] = re.sub("and","&",str(rules[i]))
#    rxnGeneMat[rxn_index_0, gene_index_0] = 1
#
##Create dictionaries to look up reactions that have metabolite mappings if the argument is supplied from the command line. 
if args.m is not None:
    args_m = str(args.m)
    md = pickle.load(open('data/models/%s' % args_m, 'rb'))
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
                                    #    exception_reactant_count += 1
                            #for product in rxns[t]['products']:
                            #    if product not in metabolite_mappings[i]:
                            #        if product in met_exception_list:
                            #            exception_product_count += 1
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
    md = pickle.load(open('data/models/%s' % args_n, 'rb'))

    #The following block are changes are specific to iMM904 model
    #Change adp to a reactant in the biomass rxns so that it does not create loops that will feed back into the system, and modify the adp and atp amount in relation to amp as described in the ratios in Wilson et al's paper. 
    rxn_index = 0
    for i in idRs:
        rxn_index += 1
        met_index = 0
        rxn_index_0 = []
        met_index_0 = []
        if i == args_b:
            rxn_index_0.append(rxn_index - 1)
            for j in idSp:
                met_index += 1
                if j == 'M_adp_c':
                    met_index_0.append(met_index - 1)
        S[met_index_0, rxn_index_0] = -0.20

    rxn_index = 0
    for i in idRs:
        rxn_index += 1
        met_index = 0
        rxn_index_0 = []
        met_index_0 = []
        if i == args_b:
            rxn_index_0.append(rxn_index - 1)
            for j in idSp:
                met_index += 1
                if j == 'M_atp_c':
                    met_index_0.append(met_index - 1)
        S[met_index_0, rxn_index_0] = -1.058

    del rxns[args_b]['products']['M_adp_c']
    rxns[args_b]['reactants']['M_adp_c'] = 0.20
    rxns[args_b]['reactants']['M_atp_c'] = 1.058

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

S = lil_matrix(csr_matrix(S)[:,rxn_index_list])

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
    metabolite_dict_file = open('data/models/%s' % args_c)
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
    #The following block of code is specific to iMM904, made to balance the carbons in the model. 
    #Make changes to other reactions so that the carbons are balanced
    rxn_index = 0
    for i in idRs:
        if i != args_b:
            rxn_index += 1
            rxn_index_0 = []
            met_index = 0
            for j in idSp:
                met_index += 1
                met_index_0 = []
                if j in rxns[i]['reactants']:
                    if rxns[i]['reactants'][j] == 0.01:
                        rxn_index_0.append(rxn_index - 1)
                        met_index_0.append(met_index - 1)
                        S[met_index_0, rxn_index_0] = -1
                        rxns[i]['reactants'][j] = 1
                if j in rxns[i]['products']:
                    if rxns[i]['products'][j] == 0.01:
                        rxn_index_0.append(rxn_index - 1)
                        met_index_0.append(met_index - 1)
                        S[met_index_0, rxn_index_0] = 1
                        rxns[i]['products'][j] = 1
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
    rxn_index_list_to_delete = []
    for i in idRs:
        rxn_index += 1
        if i not in rxns_to_remove:
            rxn_index_list.append(rxn_index-1)
        else:
            rxn_index_list_to_delete.append(rxn_index-1)
    grRules = np.delete(grRules, rxn_index_list_to_delete)
    rules = np.delete(rules, rxn_index_list_to_delete)
    c = np.delete(c, rxn_index_list_to_delete)
    subsystem = np.delete(subsystem, rxn_index_list_to_delete)

    met_index = 0
    met_index_list = []
    met_index_list_to_delete = []
    for i in idSp:
        met_index += 1
        if i not in met_exception_list:
            met_index_list.append(met_index-1)
        else:
            if ((i[-2:] != '_b') and (i[-3:] != '_b_')):
                met_index_list_to_delete.append(met_index-1)
            else:
                met_index_list.append(met_index-1)
    metNames = np.delete(metNames, met_index_list_to_delete)
    metFormulas = np.delete(metFormulas, met_index_list_to_delete)
    b = np.delete(b, met_index_list_to_delete)

    S = lil_matrix(csr_matrix(S)[met_index_list, :])
    S = lil_matrix(csr_matrix(S)[:,rxn_index_list])

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
            print("\n")    
            print(t)
            print("carbons in original model: %s" % metabolite_dict_rxns_original[t])
            print("stoichiometry in original model: %s" % rxns_original[t])
            print("carbons in adapted model: %s" % metabolite_dict_rxns[t])
            print("stoichiometry in adapted model: %s" % rxns[t])
            print("carbons in reactants in adapted model: %s" % reactants_count)
            print("cabons in products in adapted model: %s" % products_count)
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
    rxn_index_list_to_delete = []
    for i in idRs:
        rxn_index += 1
        if i not in unbalanced_rxns:
            rxn_index_list.append(rxn_index-1)
        else:
            rxn_index_list_to_delete.append(rxn_index-1)
    grRules = np.delete(grRules, rxn_index_list_to_delete)
    rules = np.delete(rules, rxn_index_list_to_delete)
    c = np.delete(c, rxn_index_list_to_delete)
    subsystem = np.delete(subsystem, rxn_index_list_to_delete)

    for i in unbalanced_rxns:
        del rxns[i]
        del gene2rxn[i]
        idRs.remove(i)

    S = lil_matrix(csr_matrix(S)[:,rxn_index_list])

    #The following is specific to iMM904
    #Remove metabolites from the biomass reaction produced from unbalanced reactions. If there is no other way to produce these metabolites, then they must be removed. This will be determined investigatively.
    mets_to_remove_from_biomass = ['M_pa_SC_c', 'M_pc_SC_c', 'M_pe_SC_c', 'M_ps_SC_c', 'M_ptd1ino_SC_c', 'M_triglyc_SC_c']

    rxn_index = 0
    for i in idRs:
        rxn_index += 1
        met_index = 0
        rxn_index_0 = []
        met_index_0 = []
        if i == args_b:
            rxn_index_0.append(rxn_index - 1)
            for j in idSp:
                met_index += 1
                if j in mets_to_remove_from_biomass:
                    met_index_0.append(met_index - 1)
                    del rxns[args_b]['reactants'][j]
        S[met_index_0, rxn_index_0] = 0

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


S = lil_matrix(csr_matrix(S)[met_index_list,:])

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
S = coo_matrix(S)

rxnGeneMat = coo_matrix(rxnGeneMat)

model = {'rxns': rxns_matlab, 'mets': mets_matlab, 'ub': ub_matlab, 'lb': lb_matlab, 'S': S, 'grRules': grRules, 'rules': rules, 'genes': genes_matlab, 'rxnGeneMat': rxnGeneMat, 'rev': rev_cobra, 'c': c, 'subsystem': subsystem, 'metNames': metNames, 'metFormulas': metFormulas, 'b': b, 'description': test_model[:-4]}
#Export the model as a .mat file for COBRA MATLAB usage
scipy.io.savemat('data/models/%s.mat' % test_model[:-4], model)

#Print all of the reactions in the final model
#for t in rxns:
#    print "%s -> %s" % (rxns[t]['reactants'], rxns[t]['products'])

#Export the model as a pickle file for EXAMO usage
exportPickle({'idSp' : idSp, 'idRs' : idRs, 'genes' : genes, 'lb' : lb, 'ub' : ub, 'gene2rxn': gene2rxn, 'S' : S, 'rxns' : rxns }, ('data/models/%s' % test_model))
