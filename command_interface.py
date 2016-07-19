import argparse

parser = argparse.ArgumentParser(prog="EXAMO-ARC")
subparsers = parser.add_subparsers(help="sub-command help")

# gene subcommand parser
parser_gene = subparsers.add_parser("gene", help="need to write")

parser_gene.add_argument("model", type=argparse.FileType("r"), help="need to write")
parser_gene.add_argument("expression_set", type=argparse.FileType("r"), help="need to write")
parser_gene.add_argument("-o", "--output", type=argparse.FileType("w"), help="need to write")
parser_gene.add_argument("-u", "--upper", type=float, default=0.2,
                         help="high help")
parser_gene.add_argument("-l", "--lower", type=float, default=0.2,
                         help="low help")
gene_group = parser_gene.add_mutually_exclusive_group(required=True)
gene_group.add_argument("-c", "--cobra", action="store_true",
                        help="cobra help")
gene_group.add_argument("-s", "--sbml", action="store_true",
                        help="sbml help")

# model subcommand parser
parser_model = subparsers.add_parser("model", help="need to write")
model_group = parser_model.add_mutually_exclusive_group(required=True)
parser_model.add_argument("model", type=argparse.FileType("r"), help="need to write")
parser_model.add_argument("-o", "--output", type=argparse.FileType("w"), help="need to write")
model_group.add_argument("-c", "--cobra", action="store_true",
                         help="cobra help")
model_group.add_argument("-s", "--sbml", action="store_true",
                         help="sbml help")
parser_model.add_argument("-e", "--external", type=str,
                          help="external help")
parser_model.add_argument("-l", "--lower-bound", type=argparse.FileType("r"),
                          help="lower bound help")
parser_model.add_argument("-u", "--uppper-bound", type=argparse.FileType("r"),
                          help="upper bound help")
parser_model.add_argument("-g", "--gene2rxn", type=argparse.FileType("r"),
                          help="gene2rxn help")
parser_model.add_argument("--metabolite-mapping-complexes", type=argparse.FileType("r"),
                          help="metab help")
parser_model.add_argument("--nucleotide-conversions", type=argparse.FileType("r"),
                          help="nuc help")
parser_model.add_argument("--metabolite2carbon", action="store_true", help="carbon help")

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

#parse args
args = parser.parse_args()
