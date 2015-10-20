# Copyright (C) 2012 Sergio Rossell
#
# This script is part of the EXAMO software
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
# 



#This package was modified by Eddie Gibb (EG). Further documentation of changes can be found in the paper.

import os
import shutil
import csv
import re
import sys
from copy import copy
from numpy import array
import multiprocessing as mp

from sys import path; path.append('./examoModules/')
from examoModules import *
from decimal import Decimal, getcontext, ROUND_DOWN
import cobra
import argparse
from Crypto.Random import atfork
import cPickle as pickle

################################################################################
# INPUTS
#Define argparse command line inputs
parser = argparse.ArgumentParser(description='_01_findZeroAndHighFreqRxns.py')
parser.add_argument('s', nargs="+", type=str, help='Necessary variable: SBML file')
parser.add_argument('p', nargs="+", type=str, help='Necessary variable: pickle file')
parser.add_argument('d', nargs="+", type=str, help='Necessary variable: condition description')
parser.add_argument('b', nargs="+", type=str, help='Necessary variable: biomass reaction')
parser.add_argument('r', nargs="+", type=int, help='Necessary variable: number of repetitions of the 2nd-4th scripts for final flux states')
parser.add_argument('-e', type=float, default=1E-10, help='Minimum value for active flux of reversible reactions in an irreversible model; default value: 1E-10')
parser.add_argument('-a', type=float, default=1E-10, help='Activity threshold (above which a flux is considered to be larger than 0); default value: 1E-10')
parser.add_argument('-t', type=float, default=1E-10, help='Activity threshold for finding active reactions; default value: 1E-10')
parser.add_argument('-EXrxns', type=str, help='Pickle file containing extracellular reactions list')
parser.add_argument('-EXtrrxns', type=str, help='Pickle file containing extracellular transport reactions list')
parser.add_argument('-Othertrrxns', type=str, help='Pickle file containing other compartmental transport reactions list')

args = parser.parse_args()

model = str(args.s)
model = model[2:-2]
pickle_model = str(args.p)
pickle_model = pickle_model[2:-2]
pickle_model_name = pickle_model[:-4]
description = str(args.d)
description = description[2:-2]
biomassRxn = str(args.b)
biomassRxn = biomassRxn[2:-2]
repetitions = str(args.r)
repetitions = int(repetitions[1:-1])
eps = args.e
activityThreshold = args.a
thresh = args.t
#EG Create lists of extracellular reactions, extracellular transport reactions, and other compartmental transport reactions, so that the reactions can be pruned in that order first.
if args.EXrxns is not None:
    EXrxns_file = str(args.EXrxns)
    md = pickle.load(open(EXrxns_file, 'rb'))
    EXrxns = md['EXrxns']
else:
    EXrxns = []
if args.EXtrrxns is not None:
    EXtrrxns_file = str(args.EXtrrxns)
    md = pickle.load(open(EXtrrxns_file, 'rb'))
    EXtrrxns = md['EXtrrxns']
else:
    EXtrrxns = []
if args.Othertrrxns is not None:
    Othertrrxns_file = str(args.Othertrrxns)
    md = pickle.load(open(Othertrrxns_file, 'rb'))
    Othertrrxns = md['Othertrrxns']
else:
    Othertrrxns = []

#Create necessary variables and import the model
if model[-4:] == '.xml':
    cobra_model = cobra.io.read_sbml_model('data/%s' % model)
if model[-4:] == '.mat':
    cobra_model = cobra.io.mat.load_matlab_model('data/%s' % model)

#model dictionary of the original model with blocked reactions deleted
fModelDict = 'data/%s' % pickle_model

cobra_model.optimize(solver='gurobi')
getcontext().rounding = ROUND_DOWN
getcontext().prec = 4
lb_biomass = Decimal(cobra_model.solution.f) + Decimal('0.0')

for repetition in range(repetitions):

	################################################################################
    # _02_minimizeNetwork_part_A.py

	################################################################################
    # INPUTS

    number_concurrent_processes = 1
    reps = 1

    md = importPickle(fModelDict)

    # Importing model information
    fbr = importPickle('data/freqBasedRxns_%s_%s.pkl' % (description, pickle_model_name))

    #EG Making subdirectories for candidate reactions
    mbaCandRxnsDirectory = 'data/mbaCandRxns/%s_%s_%s/' % (description, pickle_model_name, str(repetition))
    if os.path.exists(mbaCandRxnsDirectory):
        shutil.rmtree(mbaCandRxnsDirectory)
        os.mkdir(mbaCandRxnsDirectory, 0777)
    else:
        os.mkdir(mbaCandRxnsDirectory, 0777)

    fOutMbaCandRxns = ''.join((mbaCandRxnsDirectory, "mbaCandRxns_%s.pkl"))

	################################################################################
    # STATEMENTS
    # Instantiating CbModel 
    m0 = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'],
		md['genes'])
    #EG Changed the minimum biomass flux to be the maximum amount with default boundary constraints 
    m0.lb[m0.idRs.index(biomassRxn)] = lb_biomass

    for i in m0.idRs:
        if ((md['rxns'][i]['lb'] <= 0) and (md['rxns'][i]['ub'] == 0)):
            m0.ub[m0.idRs.index(i)] = 1000

    biomass_set = {biomassRxn}
    new_hfr = fbr['hfr'].union(biomass_set)

    zfr_check = 0
    try:
        m = deleteCbmRxns(m0, fbr['zfr'])
        act = findActiveRxns(m, thresh, new_hfr)
        cH = new_hfr & act
        act = findActiveRxns(m, thresh, cH)
        cH2 = cH & act
    except:
        zfr_check = 1
	
    if zfr_check == 1:
        zfr_list = []
        for i in fbr['zfr']:
            try:
                m1 = deleteCbmRxns(m0, i)
                act = findActiveRxns(m1, thresh, new_hfr)
                m0 = m1
                zfr_list.append(i)
            except:
                continue
        m = deleteCbmRxns(m0, zfr_list)
        act = findActiveRxns(m, thresh, new_hfr)
        cH = new_hfr & act
        act = findActiveRxns(m, thresh, cH)
        cH2 = cH & act		

    #EG Make a directory for temporary files for every time a rxn is pruned
    mbaCandRxnsDirectorySubset = 'examoModules/data/%s_%s_%s/' % (description, pickle_model_name, str(repetition))
    if not os.path.exists(mbaCandRxnsDirectorySubset):
        os.mkdir(mbaCandRxnsDirectorySubset, 0777)

    def pruneReps():
        locTime = time.localtime()
        pid = os.getpid()
        atfork()
        for x in xrange(reps):
            timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
            tag = '%s_%s_%s_%s_%s' % (description, pickle_model_name, pid, x, timeStr)
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
            exportPickle(m1, fOutModel % (description, pickle_model_name, str(repetition)))
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
            exportPickle(mDict, fOutModel % (description, pickle_model_name, str(repetition) + '_dict'))#EG End of the original 3rd script

            #EG Beginning of the 4th script
			################################################################################
            # INPUTS
            fModelExamo = 'data/examo_%s_%s_%s_dict.pkl'
            fFreqBasedRxns = 'data/freqBasedRxns_%s_%s.pkl'
            fOutMetabState = 'data/metabolicState_%s_%s_%s.csv'
			################################################################################
            # STATEMENTS

            hfr = importPickle(fFreqBasedRxns % (description, pickle_model_name))['hfr']
            md = importPickle(fModelExamo % (description, pickle_model_name, str(repetition)))
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
            f = open(fOutMetabState % (description, pickle_model_name, str(repetition)), 'w')
            csv.writer(f).writerows(nz)
            f.close()
            break
        except:
            continue
