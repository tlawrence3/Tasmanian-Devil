from __future__ import print_function
import argparse
import cobra
import csv
import decimal
import os
import shutil
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
        model_desc = model_desc + name_split[-1][:-4] + '_' + name_split[-1:-4]
    elif (description_split) == 2:
        fOutRxnsByExpression = '/'.join(description_split[0]) + '/' + 'RxnsClassifiedByExpression_' + model_desc + description_split[1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxsn = '/'.join(description_split[0]) + '/' + 'freqBasedRxns_' + model_desc + description_split[1][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        model_desc = model_desc + description_split[-1][:-4] + '_' + name_split[-1:-4]
    elif len(description_split) == 1:
        fOutRxnsByExpression = 'RxnsClassifiedByExpression_' + model_desc + description_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl'
        fOutFreqBasedRxns = 'freqBasedRxns_' + model_desc + description_split[0][:-4] + '_' + name_split[-1][:-4] + '.pkl' 
        model_desc = model_desc + description_split[0][:-4] + '_' + name_split[-1:-4]

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
        if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
            m.ub[m.idRs.index(i)] = 1000

    # 1. Classify reactions by expression
    ## Importing gene calls as a dictionary
    geneCalls = {}
    csv_file = csv.reader(args.d)
    for line in csv_file:
        geneCalls[line[0]] = int(line[1])
    f.close()

    ## Classifying reactions by expression
    rxnDict = flux_class.classifyRxnsByExpression(geneCalls, m.gene2rxn, m.genes)
    flux_class.exportPickle(rxnDict, fOutRxnsByExpression)

    # 2. Maximizing agreement score and exploring alternative optima
    rH = rxnDict['rH']
    rL = rxnDict['rL']

    #EG added eps to arguments
    model = flux_class.MetabGeneExpModel_gurobi(m.idSp, m.idRs, m.S, m.lb, m.ub, rH, rL, args.epsearly)
    scores, imgeSols = model.flux_class.exploreAlternativeOptima(idRs)

    # 3. identifying zero and high frequency reactions
    zfr, hfr = flux_class.getZeroAndHighFrequencyRxns(scores, imgeSols, idRs, args.epsidentifyingfrequencies)
    flux_class.exportPickle({'zfr' : zfr, 'hfr' : hfr}, fOutFreqBasedRxns)

    #Pint the number of HFR and ZFR genes
    print ('zfr ', len(zfr))
    print ('hfr ', len(hfr))

    ######Prune model
    repetitions = args.r
    eps = args.i
    thresh = args.t
    activityThreshold = args.a

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

        number_concurrent_processes = 10
        reps = 5

        #Making subdirectories for candidate reactions
        mbaCandRxnsDirectory = 'data/mbaCandRxns/%s_%s/' % (model_desc, str(repetition))
        if os.path.exists(mbaCandRxnsDirectory):
            shutil.rmtree(mbaCandRxnsDirectory)
            os.mkdir(mbaCandRxnsDirectory, 0777)
        else:
            os.mkdir(mbaCandRxnsDirectory, 0777)

        fOutMbaCandRxns = ''.join((mbaCandRxnsDirectory, "mbaCandRxns_%s.pkl"))

	################################################################################
        # STATEMENTS
        # Instantiating CbModel 
        m0 = flux_class.CbModel(model['S'], model['idSp'], model['idRs'], model['lb'], model['ub'], model['rxns'],
		model['genes'])
        #EG Changed the minimum biomass flux to be the maximum amount with default boundary constraints 
        m0.lb[m0.idRs.index(biomass_rxn)] = lb_biomass

        for i in m0.idRs:
            if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
                m0.ub[m0.idRs.index(i)] = 1000

        biomass_set = {biomassRxn}
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
            mtry1result.flux_class.initMipGurobi()
            mtry1result.flux_class.minSumFluxes_gurobi()
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
                    mtry1result.flux_class.initMipGurobi()
                    mtry1result.flux_class.minSumFluxes_gurobi()
                    #EG Added activityThreshold and the m0.rxns dictionary to the
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

        #EG Make a directory for temporary files for every time a rxn is pruned
        mbaCandRxnsDirectorySubset = '/data/test/%s_%s/' % (model_desc, str(repetition))
        if not os.path.exists(mbaCandRxnsDirectorySubset):
            os.mkdir(mbaCandRxnsDirectorySubset, 0777)

        def pruneReps():
            locTime = time.localtime()
            pid = os.getpid()
            Crypto.Random.atfork()
            for x in xrange(reps):
                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
                tag = '%s_%s_%s_%s' % (model_desc, pid, x, timeStr)
                try:
                    #EG Added despricription, repetition, and lists of compartmental reactions to the function
                    cr = iterativePrunning(i, m, cH2, description, pickle_model_name, biomassRxn, lb_biomass, repetition, thresh, eps, activityThreshold, EXrxns, EXtrrxns, Othertrrxns)
                    exportPickle(cr, fOutMbaCandRxns % tag)
                except:
                    print 'gurobi error, no solution found %s'  % description


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
        #and 
        # _04_predictMetabolicState.py
        #EG combined the 3rd and 4th scripts to make sure that biomass could actuall be produced for the model generated

	################################################################################
        # INPUTS
        md = importPickle(fModelDict)

        # Importing model information
        fbr = importPickle('data/freqBasedRxns_%s_%s.pkl' % (description, pickle_model_name))

        fOutModel = 'data/examo_%s_%s_%s.pkl'

	################################################################################
        # STATEMENTS
        # Instantiating CbModel 
	################################################################################
        # STATEMENTS
        # Instantiating CbModel 
        m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])

        # Retrieving MBA candidate reaction lists
        files = os.popen('ls %s | grep %s' % (mbaCandRxnsDirectory, description)).read().splitlines()


        candRxnsMinusHFR = []
        allRs = set()
        rxnSets = []
        for fn in files:
            l = set(importPickle(mbaCandRxnsDirectory + fn)) - fbr['hfr']
            rxnSets.append(importPickle(mbaCandRxnsDirectory + fn))
            allRs.update(l)
            candRxnsMinusHFR.append(l)

        # checking the degree of overlap in MBA candidate reaction lists
        overlap = []
        l = range(len(candRxnsMinusHFR))
        while l:
            ind = l.pop()
            for i in l:
                a = candRxnsMinusHFR[ind]
                b = candRxnsMinusHFR[i]
                overlap.append(float(len(a & b))/(max(len(a), len(b))))
        overlap = array(overlap)

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

        # EG Making sure that all hfr reactions are active.
        act = findActiveRxns(m, thresh, fbr['hfr'])
        mRxns = fbr['hfr'] & act

        #EG Rather than identifying which reactions need to be added to make all of the hfrs active, reactions will be added until the a flux can be achieved for the biomass reaction 
        for num in orderedFreq:
            mRxns.update(freq[num])
            excRxns = set(m.idRs) - mRxns
            try:
                m1 = deleteCbmRxns(m, excRxns)
                exportPickle(m1, fOutModel % (description, pickle_model_name_after, str(repetition)))
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
	            'descrip' : 'examo %s' % description}
                exportPickle(mDict, fOutModel % (description, pickle_model_name_after, str(repetition) + '_dict'))#EG End of the original 3rd script

                #EG Beginning of the 4th script
			################################################################################
                # INPUTS
                fModelExamo = 'data/examo_%s_%s_%s_dict.pkl'
                fFreqBasedRxns = 'data/freqBasedRxns_%s_%s.pkl'
                fOutMetabState = 'data/metabolicState_%s_%s_%s.csv'
			################################################################################
                # STATEMENTS

                hfr = importPickle(fFreqBasedRxns % (description, pickle_model_name))['hfr']
                md = importPickle(fModelExamo % (description, pickle_model_name_after, str(repetition)))
                mtry = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
                hfr = hfr & set(mtry.idRs)
                #forcing biomass production
                mtry.lb[mtry.idRs.index(biomassRxn)] = lb_biomass

                for i in mtry.idRs:
                    if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
                        mtry.ub[mtry.idRs.index(i)] = 1000

                #minimizing the sum of fluxes
                mprod = MipSeparateFwdRev_gurobi(mtry, hfr, eps)
                mprod.initMipGurobi()
                mprod.minSumFluxes_gurobi()
                #EG Added activityThreshold and the md['rxns'] dictionary to the function, so that the reactants and products could be written out
                nz = getNzRxnsGurobi(mprod, activityThreshold, md['rxns'])[1]

                # reporting the flux distribution obtained		
                f = open(fOutMetabState % (description, pickle_model_name_after, str(repetition)), 'w')
                csv.writer(f).writerows(nz)
                f.close()
                break
            except:
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
    parser_model.add_argument("extracellular", type=str,
                             help="Necessary variable: extracellular compartment abbreviation. Instead of brackets or parentheses, use underscores (ex: '_e').")
    parser_flux.add_argument("r", type=int, 
                             help='Necessary variable: number of repetitions of the 2nd-4th scripts for final flux states')
    flux_group.add_argument("-c", "--cobra", action="store_true",
                             help='Flag to specify whether model is a COBRA Toolbox (.mat) file type; must have either -xml or -mat flag.')
    flux_group.add_argument("-s", "--sbml", action="store_true",
                             help='Flag to specify whether model is a SBML (.xml) file type; must have either -xml or -mat flag.')
    parser_flux.add_argument("-e", "--epsearly", type=float, default=1E-1, 
                             help='Minimum value for fluxes forced to be non-zero (for the implementation of the pathway algorithm created by Shlomi); default value: 1E-1')
    parser_flux.add_argument("-z", "--epsidentifyingfreqencies", type=float, default=1E-10,
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
