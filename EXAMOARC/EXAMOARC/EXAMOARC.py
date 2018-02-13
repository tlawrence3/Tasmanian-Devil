from __future__ import print_function
import argparse
import cobra
import csv
import decimal
import os
import shutil
import numpy as np
import multiprocessing as mp
import Crypto.Random
import gene as gene_class
import model as model_class
import flux as flux_class


def gene(args):
    gene_class.gene_classify(args.expression_set, args.upper, args.lower, args.output)

def model(args):
    #Make sure there is a way to check if model is carbon balanced
    metFormulas_list = []
    if ((args.metabolitemappingcomplexes or args.nucleotideconversions or args.adaptation or args.zerocarbons or args.removeinactiverxnsandbalance) and not args.metabolite2carbon):
        if args.sbml:
            cobra_model = cobra.io.read_sbml_model(args.model)	
        if args.cobra:
            cobra_model = cobra.io.mat.load_matlab_model(args.model)
        metFormulas = cobra.io.mat._cell([str(m.formula) for m in cobra_model.metabolites])
        metFormulas_list = filter(None, metFormulas.tolist())
        if not metFormulas_list:
            raise RuntimeError("Must supply -d argument as well. Look at the help documetation.")

    #Import other arguments and change name of exported model file
    model_desc = ''
    if args.lowerbound:
        model_desc = model_desc + 'l'
    if args.upperbound:
        model_desc = model_desc + 'u'
    if args.gene2rxn:
        model_desc = model_desc + 'g'
    if args.metabolitemappingcomplexes:
        model_desc = model_desc + 'm'
    if args.nucleotideconversions:
        model_desc = model_desc + 'n'
    if args.removeinactiverxnsandbalance:
        model_desc = model_desc + 'c'
    if args.adaptation:
        model_desc = model_desc + 'mod'
    #Perhaps consider if OS is Windows or Linux for file structure
    name_split = args.model.split('/')
    if len(name_split) > 2:
        model_desc = '/'.join(name_split[0:-2]) + '/' + model_desc + name_split[-1]
    elif len(name_split) == 2:
        model_desc = name_split[0] + '/' + model_desc + name_split[1]
    elif len(name_split) == 1:
        model_desc = model_desc + name_split[0]
    else:
        model_desc = model_desc + args.model


    ##Make the changes to the model
    model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original, biomass_rxn = model_class.set_parameter(args.model, args.sbml, args.cobra, args.extracellular, args.lowerbound, args.upperbound, args.gene2rxn, model_desc)
    model, cobra_specific_objects = model_class.modify(model, cobra_specific_objects, args.adaptation)    
    model, cobra_specific_objects = model_class.metabolite_mapping(model, cobra_specific_objects, args.metabolitemappingcomplexes)
    model, cobra_specific_objects = model_class.nucleotide_conversion(model, cobra_specific_objects, args.nucleotideconversions)
    model, cobra_specific_objects, unbalanced_rxns_mets_unique_list, unbalanced_rxns_mets_potential_list = model_class.balance_reactions(model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original, biomass_rxn, args.metabolite2carbon, metFormulas_list, args.zerocarbons)
    model, cobra_specific_objects = model_class.metabolite_cleanup(model,cobra_specific_objects)
    model_class.model_export(model, cobra_specific_objects, model_desc)
    model_class.remove_inactive_rxns_and_account_for_biomass(model_desc,args.removeinactiverxnsandbalance, args.extracellular, args.metabolite2carbon, metFormulas_list) 
    

def flux(args):
    #Import the model, adjust the biomass production, and create the names for the exported files
    description = str(args.d)
    description = description[2:-6]
    if args.sbml:
        cobra_model = cobra.io.read_sbml_model(args.model)	
    if args.cobra:
        cobra_model = cobra.io.mat.load_matlab_model(args.model)
    
    if args.biomassprod:
        lb_biomass = args.biomassprod
    else:
        cobra_model.optimize(solver='gurobi')
        decimal.getcontext().rounding = decimal.ROUND_DOWN
        decimal.getcontext().prec = 4
        lb_biomass = decimal.Decimal(cobra_model.solution.f) + decimal.Decimal('0.0')
    #Perhaps consider if OS is Windows or Linux for file structure
    name_split = args.model.split('/')
    description_split = description.split('/')
    if len(description_split) > 2:
        fOutRxnsByExpression = '/'.join(description_split[0:-2]) + '/' + 'RxnsClassifiedByExpression_' + description_split[-1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = '/'.join(description_split[0:-2]) + '/' + 'freqBasedRxns_' + description_split[-1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
    elif (description_split) == 2:
        fOutRxnsByExpression = '/'.join(description_split[0]) + '/' + 'RxnsClassifiedByExpression_' + description_split[1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxsn = '/'.join(description_split[0]) + '/' + 'freqBasedRxns_' + description_split[1] + '_' + name_split[-1][:-4] + '.pkl'
    elif len(description_split) == 1:
        fOutRxnsByExpression = 'RxnsClassifiedByExpression_' + despriction_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = 'freqBasedRxns_' + description_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl' 

    ######Perform iMAT
    # 0. Import model dictionary
    md = importPickle(fModelDict)
    m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])

    #Identify biomass rxn
    c = np.array(cobra_model.reactions.list_attr('objective_coefficient')) * 1
    count = 0
    check_objective = 0	
    for i in c:
        count += 1
            if i == 1:
                check_objective += 1
                if check_objective == 1:
                    biomass_rxn = idRs[count-1]
                else:
                    print "More than one objective being optimized for. Change objective to be only one reaction before attempting to convert."
                    raise SystemExit

    #copy of rxns identifiers
    idRs = m.idRs[:]
    idRs.remove(biomass_rxn)
    #forcing non-zero biomass production
    m.lb[m.idRs.index(biomass_rxn)] = lb_biomass

    #Not including this would cause the model to fail the iMAT if reactions are included in the model that have the following boundary constraints 
    for i in m.idRs:
        if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
            m.ub[m.idRs.index(i)] = 1000

    # 1. Classify reactions by expression
    ## Importing gene calls as a dictionary
    geneCalls = {}
    f = open(args.d, 'rU')
    for line in csv.reader(f):
        geneCalls[line[0]] = int(line[1])
    f.close()

    ## Classifying reactions by expression
    rxnDict = classifyRxnsByExpression(geneCalls, m.gene2rxn, m.genes)
    exportPickle(rxnDict, fOutRxnsByExpression)

    # 2. Maximizing agreement score and exploring alternative optima
    rH = rxnDict['rH']
    rL = rxnDict['rL']

    #EG added eps to arguments
    model = MetabGeneExpModel_gurobi(m.idSp, m.idRs, m.S, m.lb, m.ub, rH, rL, eps)
    scores, imgeSols = model.exploreAlternativeOptima(idRs)

    # 3. identifying zero and high frequency reactions
    zfr, hfr = getZeroAndHighFrequencyRxns(scores, imgeSols, idRs, activityThreshold)
    exportPickle({'zfr' : zfr, 'hfr' : hfr}, fOutFreqBasedRxns)

    #Pint the number of HFR and ZFR genes
    print 'zfr ', len(zfr)
    print 'hfr ', len(hfr)


def main():
    parser = argparse.ArgumentParser(prog="EXAMO-ARC")
    subparsers = parser.add_subparsers(help="sub-command help")
    
    # gene subcommand parser
    ####Need to include the genes from the model as an additional argument to filter better. 
    parser_gene = subparsers.add_parser("gene", help="need to write")
    
    parser_gene.add_argument("expression_set", type=argparse.FileType("r"), help="need to write")
    parser_gene.add_argument("-o", "--output", type=argparse.FileType("w"), 
                             default="gene_classification.txt", help="need to write")
    parser_gene.add_argument("-u", "--upper", type=float, default=0.75,
                             help="high help")
    parser_gene.add_argument("-l", "--lower", type=float, default=0.25,
                             help="low help")
    parser_gene.set_defaults(func=gene)
    
    # model subcommand parser
    parser_model = subparsers.add_parser("model", help='Make modifications to models and adapt into flux module commpatible format.')
    model_group = parser_model.add_mutually_exclusive_group(required=True)
    parser_model.add_argument("model", type=str, 
                              help='Necessary variable: metabolic reconstruction file.')
    parser_model.add_argument("extracellular", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e').")
    parser_model.add_argument("-o", "--output", type=argparse.FileType("w"),
                              help='Name of model to be written out')
    model_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    model_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_model.add_argument("-l", "--lowerbound", type=argparse.FileType("r"), default=None,
                              help='lower boundary constraints file')
    parser_model.add_argument("-u", "--upperbound", type=argparse.FileType("r"), default=None,
                              help='upper boundary constraints file')
    parser_model.add_argument("-g", "--gene2rxn", type=argparse.FileType("r"), default=None,
                              help='gene2rxn file')
    parser_model.add_argument("-d", "--metabolite2carbon", type=str, default=None,
                              help='Tab-delimited file to specify dicitonary mappings of number of carbons in every metabolite. This is to check whether the model is carbon balanced. See iMM904 example for documentation.')
    parser_model.add_argument("-m", "--metabolitemappingcomplexes", type=argparse.FileType("r"), default=None,
                              help='Tab-delimited metabolite mapping complexes file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this if metFormulas is not in model.')
    parser_model.add_argument("-n", "--nucleotideconversions", type=argparse.FileType("r"), default=None,
                              help='Tab-delimited nucleotide conversions file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this if metFormulas is not in model.')
    parser_model.add_argument("-a", "--adaptation",  type=argparse.FileType("r"), default=None, 
                              help='Tab delimited file to change specific reaction stoichiometries or to remove metabolites from reactions. See iMM904 example for documentation; must have -d argument as well to use this if metFormulas is not in model.') 
    parser_model.add_argument("-z", "--zerocarbons", action="store_true",
                              help='Flag to specify whether to remove metabolites wihtout any carbons. Must have -d argument as well to use this if metFormulas is not in model.')
    parser_model.add_argument("-r", "--removeinactiverxnsandbalance", action="store_true",
                              help='Flag to specify whether to remove inactive reactions from final model and remove carbon unbalanced reactions. Recommend using only after first inspecting reactions to be removed. Must have -d argument as well to use this if metFormulas is not in model.')	
    parser_model.set_defaults(func=model)
    
    # flux subcommand parser
    parser_flux = subparsers.add_parser("flux", help='Predict condition-specific fluxes')
    flux_group = parser_flux.add_mutually_exclusive_group(required=True)
    parser_flux.add_argument("model", help='Necessary variable: metabolic reconstruction file.')
    parser_flux.add_argument("d", type=argparse.FileType=("r"), help='Necessary variable: Comma separated gene rule file')
    flux_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    flux_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_flux.add_argument("-e", "--epsearly", type=float, default=1E-1, help='Minimum value for fluxes forced to be non-zero (for the implementation of the pathway algorithm created by Shlomi); default value: 1E-1')
    parser_flux.add_argument("-a", "--activity", type=float, default=1E-10, help='Activity threshold (above which a flux is considered to be larger than 0); default value: 1E-10')
    parser_flux.add_argument("-b", "--biomassprod", type=float, help='Defined biomass production')
    parser_flux.set_defaults(func=flux)
    
    #parse args
    args = parser.parse_args()
    args.func(args)
