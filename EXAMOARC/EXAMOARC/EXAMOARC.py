from __future__ import print_function
import argparse
import gene as gene_class
import model as model_class


import cobra

def gene(args):
    gene_class.gene_classify(args.expression_set, args.upper, args.lower, args.output)

def model(args):
    if ((args.metabolitemappingcomplexes or args.nucleotideconversions or args.adaptation or args.balance) and not args.metabolite2carbon):
        raise RuntimeError("Must supply -d argument as well. Look at the help documentation.")

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
    if args.balance:
        model_desc = model_desc + 'c'
    if args.adaptation:
        model_desc = model_desc + 'mod'
    #perhaps consider if OS is Windows or Linux for file structure
    name_split = args.model.name.split('/')
    if len(name_split) > 2:
        model_desc = '/'.join(name_split[0:-2]) + '/' + model_desc + name_split[-1]
    elif len(name_split) == 2:
        model_desc = name_split[0] + '/' + model_desc + name_split[1]
    elif len(name_split) == 1:
        model_desc = model_desc + name_split[0]
    else:
        model_desc = model_desc + args.model.name


    ##Make the changes to the model
    model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original = model_class.set_parameter(args.model, args.sbml, args.cobra, args.extracellular, args.lowerbound, args.upperbound, args.gene2rxn, model_desc)
    model, cobra_specific_objects = model_class.modify(model, cobra_specific_objects, args.adaptation)    
    model, cobra_specific_objects = model_class.metabolite_mapping(model, cobra_specific_objects, args.metabolitemappingcomplexes)
    model, cobra_specific_objects = model_class.nucleotide_conversion(model, cobra_specific_objects, args.nucleotideconversions)
    model, cobra_specific_objects = model_class.balance_reactions(model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original, args.biomassRxn, args.metabolite2carbon, args.zerocarbons, args.balance)
    model, cobra_specific_objects = model_class.metabolite_cleanup(model,cobra_specific_objects)
    model_matlab = model_class.model_export(model, cobra_specific_objects, model_desc)
    #Test functionality
    #print ("%s" % model_desc)
    #cobra_model = cobra.io.mat.load_matlab_model("lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat")
    #cobra_model.optimize(solver="gurobi")
    #print(cobra_model.solution.f)
    

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
    parser_model = subparsers.add_parser("model", help="Make modifications to models and adapt into EXAMO commpatible format")
    model_group = parser_model.add_mutually_exclusive_group(required=True)
    parser_model.add_argument("model", type=argparse.FileType("r"), 
                              help="Necessary variable: metabolic reconstruction file")
    parser_model.add_argument("biomassRxn", type=str,
                             help="Necessary variable: biomass rxn; Append 'R_' to the front if the reaction name does not begin with this.")
###Need to throw in a check here to make sure that the biomass rxn follows the naming convention.
    parser_model.add_argument("extracellular", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets, use underscores (ex: '_e').)")
    parser_model.add_argument("-o", "--output", type=argparse.FileType("w"),
                              help="Name of model to be written out")
    model_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    model_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_model.add_argument("-l", "--lowerbound", type=argparse.FileType("r"), default=None,
                              help="'lower boundary constraints file")
    parser_model.add_argument("-u", "--upperbound", type=argparse.FileType("r"), default=None,
                              help="upper boundary constraints file")
    parser_model.add_argument("-g", "--gene2rxn", type=argparse.FileType("r"), default=None,
                              help="'gene2rxn file'")
    parser_model.add_argument("-d", "--metabolite2carbon", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited file to specify dicitonary mappings of number of carbons in every metabolite. This is to check whether the model is carbon balanced. See iMM904 example for documentation.")
    parser_model.add_argument("-m", "--metabolitemappingcomplexes", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited metabolite mapping complexes file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this.")
    parser_model.add_argument("-n", "--nucleotideconversions", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited nucleotide conversions file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this.")
    parser_model.add_argument("-a", "--adaptation",  type=argparse.FileType("r"), default=None, 
                              help="Tab delimited file to change specific reaction stoichiometries or to remove metabolites from reactions. See iMM904 example for documentation; must have -d argument as well to use this.") 
    parser_model.add_argument("-z", "--zerocarbons", action="store_true",
                              help="Flag to specify whether to remove metabolites wihtout any carbons. Must have -d argument as well to use this.")
    parser_model.add_argument("-b", "--balance", action="store_true", 
                              help="Flag to specify whether to remove carbon unbalanced reactions. Recommend using only after first inspecting reactions to be removed. Must have -d argument as well to use this.")	
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
