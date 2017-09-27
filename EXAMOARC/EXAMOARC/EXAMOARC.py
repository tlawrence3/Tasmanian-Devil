from __future__ import print_function
import argparse
import cobra
import gene as gene_class
import model as model_class


import cobra

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
    print("flux subcommand")

def main():
    parser = argparse.ArgumentParser(prog="EXAMO-ARC")
    subparsers = parser.add_subparsers(help="sub-command help")
    
    # gene subcommand parser
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
    parser_model = subparsers.add_parser("model", help='Make modifications to models and adapt into EXAMO commpatible format.')
    model_group = parser_model.add_mutually_exclusive_group(required=True)
    parser_model.add_argument("model", type=str, 
                              help='Necessary variable: metabolic reconstruction file.')
    parser_model.add_argument("extracellular", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets, use underscores (ex: '_e').")
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
    parser_flux = subparsers.add_parser("flux", help="need to write")
    flux_group = parser_flux.add_mutually_exclusive_group(required=True)
    parser_flux.add_argument("model", type=argparse.FileType("r"), help="need to write")
    parser_flux.add_argument("-o", "--output", type=str, help="need to write")
    flux_group.add_argument("-c", "--cobra", action="store_true",
                             help="cobra help")
    flux_group.add_argument("-s", "--sbml", action="store_true",
                             help="sbml help")
    parser_flux.add_argument("-b", "--biomass-rxn", type=str, help="biomass help")
    parser_flux.add_argument("-a", type=float, help="a")
    parser_flux.set_defaults(func=flux)
    
    #parse args
    args = parser.parse_args()
    args.func(args)
