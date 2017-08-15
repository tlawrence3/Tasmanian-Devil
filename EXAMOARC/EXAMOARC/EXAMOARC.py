from __future__ import print_function
import argparse
import gene as gene_class
import model as model_class

def gene(args):
    gene_class.gene_classify(args.expression_set, args.upper, args.lower, args.output)

def model(args):
    if ((args.metabolite-mapping-complexes or args.nucleotide-conversions or args.adaptation or args.balance) and not args.metabolite2carbon):
        raise RuntimeError("Must supply -d argument as well. Look at the help documentation.")

    #Import the metabolic reconstruction file name
    #args_s = args.model.name
    #model_file = args_s[2:-2]

    #Import biomass reaction
    #biomassRxn = str(args.biomassRxn)
    #biomassRxn = biomassRxn[2:-2]

    #Import extracellular compartment abbreviation
    #extracellular = str(args.extracellularcompartmentname)
    #extracellular = extracellular[2:-2]

    #Import other arguments and change name of exported model file
    #Prepend the code
    #test_model = model_file[:-4]
    test_model = ''
    if args.lower-bound:
        test_model = test_mdoel + 'l'
    if args.upper-bound:
        test_model = test_model + 'u'
    if args.gene2rxn:
        test_model = test_model + 'g'
    if args.metabolite-mapping-complexes:
        test_model = test_model + 'm'
    if args.nucleotide-conversions:
        test_model = test_model + 'n'
    if args.balance:
        test_model = test_model + 'c'
    if args.adaptation:
        test_model = test_model + 'mod'
    test_model = test_model + args.model.name

    #Import the model
    #Could just be read_sbml_model(args.model)
    if args.sbml:
        cobra_model = cobra.io.read_sbml_model(args.model)
    if args.cobra:
        cobra_model = cobra.io.mat.load_matlab_model(args.model)

    ##Make the changes to the model
    model = model_class.set_parameter(args.model, args.sbml, args.cobra, args.extracellular, args.lower-bound, args.upper-boound, args.gene2rxn)
    #model = model.metabolite_mapping(model, args_m)
    #model = model.nucleotide_conversion(model, args_n)
    #model = model.modfiy(model, args_mod)
    #model = mdoel.balance_reactions(model, args_c2m, args_c)

    

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
    parser_gene.add_argument("-l", "--lower", type=float, default=0.2,
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
    parser_model.add_argument("externalcompartmentname", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets, use underscores (ex: '_e').)")
    parser_model.add_argument("-o", "--output", type=argparse.FileType("w"),
                              help="Name of model to be written out")
    model_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    model_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_model.add_argument("-l", "--lower-bound", type=argparse.FileType("r"), default=None,
                              help="'lower boundary constraints file")
    parser_model.add_argument("-u", "--uppper-bound", type=argparse.FileType("r"), default=None,
                              help="upper boundary constraints file")
    parser_model.add_argument("-g", "--gene2rxn", type=argparse.FileType("r"), default=None,
                              help="'gene2rxn file'")
    parser_model.add_argument("-d", "--metabolite2carbon", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited file to specify dicitonary mappings of number of carbons in every metabolite. This is to check whether the model is carbon balanced. See iMM904 example for documentation.")
    parser_model.add_argument("-m", "--metabolite-mapping-complexes", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited metabolite mapping complexes file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this.")
    parser_model.add_argument("-n", "--nucleotide-conversions", type=argparse.FileType("r"), default=None,
                              help="Tab-delimited nucleotide conversions file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -d argument as well to use this.")
    parser_model.add_argument("-a", "--adaptation",  type=argparse.FileType("r"), default=None, 
                              help='Tab delimited file to change specific reaction stoichiometries or to remove metabolites from reactions. See iMM904 example for documentation; must have -d argument as well to use this.') 
    parser_model.add_argument("-b", "--balance", action="store_true", 
                              help='Flag to specify whether to remove carbon unbalanced reactions. Recommend using only after first inspecting reactions to be removed. Must have -d argument as well to use this.')	
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
