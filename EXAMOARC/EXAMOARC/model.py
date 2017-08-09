from conversion import *
from utility import *
import cobra

##Define argparse command line inputs
parser = argparse.ArgumentParser(description='Make modifications to models and adapt into EXAMO commpatible format')
parser.add_argument('s', nargs="+", type=str, help='Necessary variable: metabolic reconstruction file')
parser.add_argument('b', nargs="+", type=str, help='Necessary variable: biomass rxn; Append "R_" to the front if the reaction name does not begin with this')
parser.add_argument('e', nargs="+", type=str, help='Necessary variable: extracellular compartment abbreviation. Instead of brackets, use underscores (ex: "_e")')
parser.add_argument('-xml', type=str, help='Flag to specify whether model is .xml; must have either -xml or -mat flag')
parser.add_argument('-mat', type=str, help='Flag to specify whether model is .mat; must have either -xml or -mat flag')
parser.add_argument('-l', type=str, help='lb file')
parser.add_argument('-u', type=str, help='ub file')
parser.add_argument('-g', type=str, help='gene2rxn file')
parser.add_argument('-c2m', type=str, help='Tab delimited file to specify number of carbons in every metabolite. This is to check whether the model is carbon balanced. See iMM904 example for documentation.')
parser.add_argument('-m', type=str, help='Tab delimited metabolite mapping complexes file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -c2m argument as well to use this.')
parser.add_argument('-n', type=str, help='Tab delimited nucleotide conversions file. See iMM904 example for documentation. Make sure model is carbon balanced if you use this; must have -c2m argument as well to use this.')
parser.add_argument('-mod', type=str, help='Tab delimited file to change specific reaction stoichiometries or to remove metabolites from reactions. See iMM904 example for documentation; must have -c2m argument as well to use this.') 
parser.add_argument('-c', type=str, help='Flag to specify whether to remove carbon unbalanced reactions. Recommend using only after first inspecting reactions to be removed. Must have -c2m argument as well to use this.')
#Need to figure out if removing carbons with 0 metabolites on top of removing carbons that are not balanced. 


##Import argparse arguments
args = parser.parse_args()

if args.xml is None and args.mat is None:
    raise RuntimeError("Must supply -xml or -mat flag. Look at the help documentation.")

if args.m is not None and args.c2m is None:
    raise RuntimeError("Must supply -c argument as well. Look at the help documentation.")

if args.n is not None and args.c2m is None:
    raise RuntimeError("Must supply -c argument as well. Look at the help documentation.")

if args.mod is not None and args.c2m is None:
    raise RuntimeError("Must supply -c argument. Look at the help documentation.")

if args.c is not None and args.c2m is None:
    raise RuntimeError("Must supply -c argument. Look at the help documentation.")

#Import the metabolic reconstruction file name
args_s = str(args.s)
model_file = args_s[2:-2]

#Import biomass reaction
args_b = str(args.b)
biomass_rxn = args_b[2:-2]

#Import extracellular compartment abbreviation
extracellular = str(args.e)
extracellular = extracellular[2:-2]

#Name of exported modle file
test_model = model_file[:-4]
if args.l is not None:
	args_l = str(args.l)
	test_model = test_model + str('_l_') + args_l[:-4]
if args.u is not None:
	test_model = test_model + str('_u')
if args.g is not None:
	test_model = test_model + str('_g')
if args.m is not None:
	test_model = test_model + str('_m')
if args.n is not None:
	test_model = test_model + str('_n')
if args.mod is not None:
	args_mod = str(args.mod)
	test_model = test_model + str('_mod_') + args_mod[:-4]
if args.c is not None:
	test_model = test_model + str('_c')

#Import the model
if args.xml is not None:
	cobra_model = cobra.io.read_sbml_model('%s' % model_file)
if args_s.mat is not None:
	cobra_model = cobra.io.mat.load_matlab_model('%s' % model_file)

#Import the lower boundary file name if it exists
ars_l = None
if args.l is not None:
	args_l = str(args.l)

#Import the upper boundary file name if it exists
args_u = None
if args.u is not None:
	args_u = str(args.u)

#Import gene rule adjustments file name if it exists 
args_g = None
if args.g is not None:
	args_g = str(args.g)

#Import carbon to metabolite file name if it exists
args_c2m = None
if args.c2m is not None:
	args_c2m = str(args.c2m)

#Import metabolite mapping complexes file name if it exists
args_m = None
if args.m is not None:
	args_m = str(args.m)

#Import nucleotide conversions file name if it exists
args_n = None
if args.n is not None:
	args_n = str(args.n)

#Import model specific changes file name if it exists
args_mod = None
if args.mod is not None:
	args_mod = str(args.mod)

#Import whether to balance carbons or not
args_c = False
if args.c is not None:
	args_c = True


##Make the changes to the model
model = cobra_model.set_parameter(cobra_model, extracellular, args_l, args_u, args_g)
model = model.metabolite_mapping(model, args_m)
model = model.nucleotide_conversion(model, args_n)
model = model.modfiy(model, args_mod)
model = mdoel.balance_reactions(model, args_c2m, args_c)


