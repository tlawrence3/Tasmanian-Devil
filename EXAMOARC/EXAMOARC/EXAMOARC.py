from __future__ import print_function
import argparse
import cobra
import csv
import decimal
import os
import shutil
import time
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
    description = args.d.name
    if args.sbml:
        cobra_model = cobra.io.read_sbml_model(args.model)	
    if args.cobra:
        cobra_model = cobra.io.mat.load_matlab_model(args.model)
    
    if args.biomassprod:
        lb_biomass = args.biomassprod
    else:
        solution = cobra_model.optimize()
        decimal.getcontext().rounding = decimal.ROUND_DOWN
        decimal.getcontext().prec = 4
        lb_biomass = decimal.Decimal(solution.objective_value) + decimal.Decimal('0.0')
    #Perhaps consider if OS is Windows or Linux for file structure
    name_split = args.model.split('/')
    description_split = description.split('/')
    model_desc = ''
    if args.EXrxns:
        model_desc = model_desc + 'EXrxns'
    if args.EXtrrxns:
        model_desc = model_desc + 'EXtrrxns'
    if args.Othertrrxns:
        model_desc = model_desc + 'Othertrrxns'
    if len(description_split) > 2:
        fOutRxnsByExpression = '/'.join(description_split[0:-2]) + '/' + 'RxnsClassifiedByExpression_' + model_desc + description_split[-1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = '/'.join(description_split[0:-2]) + '/' + 'freqBasedRxns_' + model_desc + description_split[-1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        model_desc = model_desc + description_split[-1][:-4] + '_' + name_split[-1][:-4]
        file_path = '/'.join(description_split[0:-2]) + '/'
    elif len(description_split) == 2:
        fOutRxnsByExpression = description_split[0] + '/' + 'RxnsClassifiedByExpression_' + model_desc + description_split[1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = description_split[0] + '/' + 'freqBasedRxns_' + model_desc + description_split[1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        model_desc = model_desc + description_split[-1][:-4] + '_' + name_split[-1][:-4]
        file_path = description_split[0] + '/'
    elif len(description_split) == 1:
        fOutRxnsByExpression = 'RxnsClassifiedByExpression_' + model_desc + description_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = 'freqBasedRxns_' + model_desc + description_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl' 
        model_desc = model_desc + description_split[0][:-4] + '_' + name_split[-1][:-4]
        file_path = ''

    ######Perform iMAT
    # 0. Import model dictionary
    args_lowerbound = None
    args_upperbound = None
    args_gene2rxn = None
    model, cobra_specific_objects, mets_to_extracellular_comp, rxns_original, biomass_rxn = model_class.set_parameter(args.model, args.sbml, args.cobra, args.extracellular, args_lowerbound, args_upperbound, args_gene2rxn)
    m = flux_class.CbModel(model['S'], model['idSp'], model['idRs'], model['lb'], model['ub'], model['rxns'], model['genes'])

    #copy of rxns identifiers
    idRs = m.idRs[:]
    idRs.remove(biomass_rxn)
    #forcing non-zero biomass production
    m.lb[m.idRs.index(biomass_rxn)] = lb_biomass

    #Not including this would cause the model to fail the iMAT if reactions are included in the model that have the following boundary constraints 
    for i in m.idRs:
        if ((model['rxns'][i]['lb'] <= 0) and (model['rxns'][i]['ub'] == 0)):
            m.ub[m.idRs.index(i)] = 1000

    # 1. Classify reactions by expression
    ## Importing gene calls as a dictionary
    geneCalls = {}
    csv_file = csv.reader(args.d)
    for line in csv_file:
        geneCalls[line[0]] = int(line[1])

    ## Classifying reactions by expression
    rxnDict = flux_class.classifyRxnsByExpression(geneCalls, m.gene2rxn, m.genes)
    flux_class.exportPickle(rxnDict, fOutRxnsByExpression)

    # 2. Maximizing agreement score and exploring alternative optima
    rH = rxnDict['rH']
    rL = rxnDict['rL']

    #EG added eps to arguments
    model_mgem = flux_class.MetabGeneExpModel_gurobi(m.idSp, m.idRs, m.S, m.lb, m.ub, rH, rL, args.epsearly)
    scores, imgeSols = model_mgem.exploreAlternativeOptima(idRs)

    # 3. identifying zero and high frequency reactions
    zfr, hfr = flux_class.getZeroAndHighFrequencyRxns(scores, imgeSols, idRs, args.epsidentifyingfrequencies)
    flux_class.exportPickle({'zfr' : zfr, 'hfr' : hfr}, fOutFreqBasedRxns)

    #Pint the number of HFR and ZFR genes
    print ('zfr ', len(zfr))
    print ('hfr ', len(hfr))

    ######Prune model
    repetitions = args.r
    eps = args.epsirreversiblerxns
    thresh = args.epspruning
    activityThreshold = args.espsolving

    #Create lists of extracellular reactions, extracellular transport reactions, and other compartmental transport reactions, so that the reactions can be pruned in that order first.
    EXrxns = []
    EXtrrxns = []
    Othertrrxns = []
    if args.EXrxns:
        csv_file = csv.reader(args.EXrxns)
	for line in csv_file:
            EXrxns.append(line[0])
    if args.EXtrrxns:
        csv_file = csv.reader(args.EXtrrxns)
        for line in csv_file:
            EXtrrxns.append(line[0])
    if args.Othertrrxns:
        csv_file = csv.reader(args.Othertrrxns)
        for line in csv_file:
            Othertrrxns.append(line[0])
    for repetition in range(args.r): 
    ################################################################################
    # originally _02_minimizeNetwork_part_A.py in EXAMO
    ################################################################################

        number_concurrent_processes = 1
        reps = 1

        #Making subdirectories for candidate reactions
        mbaCandRxnsDirectory = file_path + 'data/mbaCandRxns/%s_%s/' % (model_desc, str(repetition))
        if os.path.exists(mbaCandRxnsDirectory):
            shutil.rmtree(mbaCandRxnsDirectory)
            os.makedirs(mbaCandRxnsDirectory, 0777)
        else:
            os.makedirs(mbaCandRxnsDirectory, 0777)

        fOutMbaCandRxns = ''.join((mbaCandRxnsDirectory, "mbaCandRxns_%s.pkl"))

	################################################################################
        # STATEMENTS
        # Instantiating CbModel 
        m0 = flux_class.CbModel(model['S'], model['idSp'], model['idRs'], model['lb'], model['ub'], model['rxns'],model['genes'])
        #EG Changed the minimum biomass flux to be the maximum amount with default boundary constraints 
        m0.lb[m0.idRs.index(biomass_rxn)] = lb_biomass

        for i in m0.idRs:
            if ((model['rxns'][i]['lb'] <= 0) and (model['rxns'][i]['ub'] == 0)):
                m0.ub[m0.idRs.index(i)] = 1000

        biomass_set = {biomass_rxn}
        new_hfr = hfr.union(biomass_set)

        zfr_check = 0
        try:
            m = flux_class.deleteCbmRxns(m0, zfr)
            act = flux_class.findActiveRxns(m, thresh, new_hfr)
            cH = new_hfr & act
            act = flux_class.findActiveRxns(m, thresh, cH)
            cH2 = cH & act
            #######################################################################
            # STATEMENTS
            new_hfr2 = hfr & set(m.idRs)
            #minimizing the sum of fluxes
            mtry1result = flux_class.MipSeparateFwdRev_gurobi(m, new_hfr2, eps)
            mtry1result.initMipGurobi()
            mtry1result.minSumFluxes_gurobi()
            #EG Added activityThreshold and the m0.rxns dictionary to the
            #function, so that the reactants and products could be written out
            nz = flux_class.getNzRxnsGurobi(mtry1result, activityThreshold, m.rxns)[1]
        except:
            zfr_check = 1
	
        if zfr_check == 1:
            zfr_list = []
            for i in zfr:
                try:
                    m1 = flux_class.deleteCbmRxns(m0, i)
                    act = flux_class.findActiveRxns(m1, thresh, new_hfr)
                    #######################################################################
                    # STATEMENTS
                    new_hfr2 = hfr & set(m0.idRs)
                    #minimizing the sum of fluxes
                    mtry1result = flux_class.MipSeparateFwdRev_gurobi(m1, new_hfr2, eps)
                    mtry1result.initMipGurobi()
                    mtry1result.minSumFluxes_gurobi()
                    #Added activityThreshold and the m0.rxns dictionary to the
                    #function, so that the reactants and products could be written out
                    nz = flux_class.getNzRxnsGurobi(mtry1result, activityThreshold, m1.rxns)[1]
                    m0 = m1
                    zfr_list.append(i)
                except:
                    continue
            m = flux_class.deleteCbmRxns(m0, zfr_list)
            act = flux_class.findActiveRxns(m, thresh, new_hfr)
            cH = new_hfr & act
            act = flux_class.findActiveRxns(m, thresh, cH)
            cH2 = cH & act		

        #Make a directory for temporary files for every time a rxn is pruned
        mbaCandRxnsDirectorySubset = file_path + 'data/test/%s_%s/' % (model_desc, str(repetition))
        if not os.path.exists(mbaCandRxnsDirectorySubset):
            os.makedirs(mbaCandRxnsDirectorySubset, 0777)

        def pruneReps():
            locTime = time.localtime()
            pid = os.getpid()
            Crypto.Random.atfork()
            for x in xrange(reps):
                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
                tag = '%s_%s_%s_%s' % (model_desc, pid, x, timeStr)
                #Added despricription, repetition, and lists of compartmental reactions to the function
                cr = flux_class.iterativePrunning(i, m, cH2, fOutFreqBasedRxns, biomass_rxn, lb_biomass, repetition, thresh, eps, activityThreshold, EXrxns, EXtrrxns, Othertrrxns)
                flux_class.exportPickle(cr, fOutMbaCandRxns % tag)

        processes = []
        for _ in xrange(number_concurrent_processes):
            p = mp.Process(target = pruneReps)
            p.start()
            processes.append(p)
   
        for p in processes:
            p.join()

        #EG Delete the temporary files generated for every time a rxn is pruned
        shutil.rmtree(mbaCandRxnsDirectorySubset)

	################################################################################
        # _03_minimizeNetwork_part_B_new.py 
	################################################################################
        # STATEMENTS
        #Create file to export
        fOutModel = file_path + 'examo_%s_%s.pkl'

        # Instantiating CbModel 
        m = flux_class.CbModel(model['S'], model['idSp'], model['idRs'], model['lb'], model['ub'], model['rxns'], model['genes'])

        # Retrieving MBA candidate reaction lists
        files = os.popen('ls %s | grep %s' % (mbaCandRxnsDirectory, model_desc)).read().splitlines()

        rxnSets = []
        for fn in files:
            l = set(flux_class.importPickle(mbaCandRxnsDirectory + fn)) - hfr
            rxnSets.append(flux_class.importPickle(mbaCandRxnsDirectory + fn))

        # Quantifying the number of times a rxn is among the candidate models
        rxnFreq = {}
        for rs in rxnSets:
            for rxn in rs:
                try:
                    rxnFreq[rxn] += 1
                except KeyError:
                    rxnFreq[rxn] = 1

        freq = {}#{frequency : list of reactions}
        for rxn in rxnFreq:
            try:
                freq[rxnFreq[rxn]].add(rxn)
            except KeyError:
                freq[rxnFreq[rxn]] = set([rxn])
        orderedFreq = freq.keys()
        orderedFreq.sort(reverse = True)

        #Making sure that all hfr reactions are active.
        act = flux_class.findActiveRxns(m, thresh, hfr)
        mRxns = hfr & act
        #Rather than identifying which reactions need to be added to make all of the hfrs active, reactions will be added until the a flux can be achieved for the biomass reaction 
        for num in orderedFreq:
            mRxns.update(freq[num])
            excRxns = set(m.idRs) - mRxns
            try:
                m1 = flux_class.deleteCbmRxns(m, excRxns)
                flux_class.exportPickle(m1, fOutModel % (model_desc, str(repetition)))
                rev = [0 if val >= 0 else 1 for val in m1.lb]
                gene2rxn = {}
                rxns = {}
                for rxn in m1.rxns:
                    gene2rxn[rxn] = m1.rxns[rxn]['genes']
                    rxns[rxn] = m1.rxns[rxn]
                mDict = {
                    'S' : m1.S,
                    'idSp' : m1.idSp,
                    'idRs' : m1.idRs,
                    'lb' : m1.lb,
                    'ub' : m1.ub,
                    'rev' : rev,
                    'genes' : m1.genes,
                    'gene2rxn' : gene2rxn,
                    'rxns' : rxns,
	            'descrip' : 'examo %s' % model_desc}
                flux_class.exportPickle(mDict, fOutModel % (model_desc, str(repetition) + '_dict'))
                ################################################################################                
                # _04_predictMetabolicState.py
		################################################################################
                # OUTPUT
                fOutMetabState = file_path + 'metabolicState_%s_%s.csv'
		################################################################################
                # STATEMENTS

                md = flux_class.importPickle(fOutModel % (model_desc, str(repetition) + '_dict'))
                mtry = flux_class.CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
                hfr = hfr & set(mtry.idRs)
                #forcing biomass production
                mtry.lb[mtry.idRs.index(biomass_rxn)] = lb_biomass
                for i in mtry.idRs:
                    if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
                        mtry.ub[mtry.idRs.index(i)] = 1000

                #minimizing the sum of fluxes
                mprod = flux_class.MipSeparateFwdRev_gurobi(mtry, hfr, eps)
                mprod.initMipGurobi()
                mprod.minSumFluxes_gurobi()
                #Added activityThreshold and the md['rxns'] dictionary to the function, so that the reactants and products could be written out
                nz = flux_class.getNzRxnsGurobi(mprod, activityThreshold, md['rxns'])[1]

                # reporting the flux distribution obtained		
                f = open(fOutMetabState % (model_desc, str(repetition)), 'w')
                csv.writer(f).writerows(nz)
                f.close()
                break
            except:
                mprod = flux_class.MipSeparateFwdRev_gurobi(mtry, hfr, eps)
                mprod.initMipGurobi()
                mprod.minSumFluxes_gurobi()
                continue


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
    parser_flux.add_argument("d", type=argparse.FileType("r"), help='Necessary variable: comma separated gene rule file')
    parser_flux.add_argument("extracellular", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e').")
    parser_flux.add_argument("r", type=int, 
                             help='Necessary variable: number of repetitions of the 2nd-4th scripts for final flux states')
    flux_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    flux_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_flux.add_argument("-e", "--epsearly", type=float, default=1E-1, 
                             help='Minimum value for fluxes forced to be non-zero (for the implementation of the pathway algorithm created by Shlomi); default value: 1E-1')
    parser_flux.add_argument("-z", "--epsidentifyingfrequencies", type=float, default=1E-10,
                             help='Activity threshold when classifying HFR and ZFR reactions from gene rules') 
    parser_flux.add_argument("-i", "--epsirreversiblerxns", type=float, default=1E-10,
                             help='Minimum value for active flux of reversible reactions in an irreversible model; default value: 1E-10')
    parser_flux.add_argument("-t", "--epspruning", type=float, default=1E-10, 
                             help='Activity threshold for finding active reactions when pruning a model')
    parser_flux.add_argument("-a", "--espsolving", type=float, default=1E-10, 
                             help='Activity threshold (above which a flux is considered to be larger than 0 when minimizing the sum of fluxes); default value: 1E-10')
    parser_flux.add_argument("-b", "--biomassprod", type=float, help='Defined biomass production') 
    parser_flux.add_argument('-EXrxns', type=argparse.FileType("r"), help='Comma separated file containing extracellular reactions')
    parser_flux.add_argument('-EXtrrxns', type=argparse.FileType("r"), help='Comma separated file containing extracellular transport reactions')
    parser_flux.add_argument('-Othertrrxns', type=argparse.FileType("r"), help='Comma separated file containing other compartmental transport reactions')
    parser_flux.set_defaults(func=flux)
    
    #parse args
    args = parser.parse_args()
    args.func(args)
