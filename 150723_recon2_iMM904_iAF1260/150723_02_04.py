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

from sys import path; path.append('./examoModules/')
from examoModules import *
from decimal import Decimal, getcontext, ROUND_DOWN
import cobra

################################################################################
# INPUTS
model = sys.argv[1]
pickle_model = sys.argv[2]
description = sys.argv[3]
biomassRxn = sys.argv[4]
repetitions = int(sys.argv[5])

#Create necessary variables and import the model
if model[-4:] == '.xml':
    cobra_model = cobra.io.read_sbml_model('data/%s' % model)
if model[-4:] == '.mat':
    cobra_model = cobra.io.mat.load_matlab_model('data/%s' % model)

# model dictionary of the original model with blocked reactions deleted
fModelDict = 'data/%s' % pickle_model

cobra_model.optimize(solver='gurobi')
getcontext().rounding = ROUND_DOWN
getcontext().prec = 4
lb_biomass = Decimal(cobra_model.solution.f) + Decimal('0.0')

# Algorithm variables
#threshold above which a flux is considered to be larger than zero
activityThreshold = 1E-10

for repetition in range(repetitions):

	################################################################################
    # _02_minimizeNetwork_part_A.py

	################################################################################
    # INPUTS

    numProc = 1# 100	#number of parallel processes used
    numRep =1# 10		#number of times each process is repeated

    md = importPickle(fModelDict)

    # Importing model information
    fbr = importPickle('data/freqBasedRxns_%s.pkl' % description)

    #EG Making subdirectories for candidate reactions
    mbaCandRxnsDirectory = 'data/mbaCandRxns/%s_%s/' % (description, str(repetition))
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
        act = findActiveRxns(m, 1E-10, new_hfr)
        cH = new_hfr & act
        act = findActiveRxns(m, 1E-10, cH)
        cH2 = cH & act
    except:
        zfr_check = 1
	
    if zfr_check == 1:
        zfr_list = []
        for i in fbr['zfr']:
            try:
                m1 = deleteCbmRxns(m0, i)
                act = findActiveRxns(m1, 1E-10, new_hfr)
                m0 = m1
                zfr_list.append(i)
            except:
                continue
        m = deleteCbmRxns(m0, zfr_list)
        act = findActiveRxns(m, 1E-10, new_hfr)
        cH = new_hfr & act
        act = findActiveRxns(m, 1E-10, cH)
        cH2 = cH & act		

    #EG Create lists of extracellular reactions, extracellular transport reactions, and other compartmental transport reactions, so that the reactions can be pruned in that order first.
    EXrxns = []

    EXtrrxns = [] 

    Othertrrxns = []

    #EG Make a directory for temporary files for every time a rxn is pruned
    mbaCandRxnsDirectorySubset = 'examoModules/data/%s_%s/' % (description, str(repetition))
    if not os.path.exists(mbaCandRxnsDirectorySubset):
        os.mkdir(mbaCandRxnsDirectorySubset, 0777)

    #Run the MBA
    for i in range(numProc):
        pid = os.fork()
        if pid == 0:
            for j in range(numRep):
                locTime = time.localtime()
                pid = os.getpid()
                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
                tag = '%s_%s_%s' % (description, pid, timeStr)
                try:
                    #EG Added despricription, repetition, and lists of compartmental reactions to the function
                    cr = iterativePrunning(i, m, cH2, description, biomassRxn, lb_biomass, repetition, activityThreshold, EXrxns, EXtrrxns, Othertrrxns)
                    exportPickle(cr, fOutMbaCandRxns % tag)
                except:
                    print 'gurobi error, no solution found %s'  % description
            os._exit(0)

    #EG Count the number of exported pickle files, and do not proceed until all files have been written
    file_count = 0
    while file_count != numProc*numRep:
        for file_directory, dirs, files in os.walk(mbaCandRxnsDirectory):
            file_count = len(files)

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
    fbr = importPickle('data/freqBasedRxns_%s.pkl' % description)

    fOutModel = 'data/examo_%s_%s.pkl'

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
    act = findActiveRxns(m, 1E-10, fbr['hfr'])
    mRxns = fbr['hfr'] & act

    #EG Rather than identifying which reactions need to be added to make all of the hfrs active, reactions will be added until the a flux can be achieved for the biomass reaction 
    for num in orderedFreq:
        mRxns.update(freq[num])
        excRxns = set(m.idRs) - mRxns
        try:
            m1 = deleteCbmRxns(m, excRxns)
            exportPickle(m1, fOutModel % (description, str(repetition)))
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
            exportPickle(mDict, fOutModel % (description, str(repetition) + '_dict'))#EG End of the original 3rd script

            #EG Beginning of the 4th script
			################################################################################
            # INPUTS
            eps = 1E-10

            fModelExamo = 'data/examo_%s_%s_dict.pkl'
            fFreqBasedRxns = 'data/freqBasedRxns_%s.pkl'

            fOutMetabState = 'data/metabolicState_%s_%s.csv'

			################################################################################
            # STATEMENTS

            hfr = importPickle(fFreqBasedRxns % description)['hfr']
            md = importPickle(fModelExamo % (description, str(repetition)))
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
            f = open(fOutMetabState % (description, str(repetition)), 'w')
            csv.writer(f).writerows(nz)
            f.close()
            break
        except:
            continue
