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

"""
MBA implementation
"""
from numpy import array
from gurobipy import *
from multiprocessing import Pool
import itertools
import random 
import time 
import os
import csv
import re

from sys import path; path.append('./modules/')
from utilities import importPickle, exportPickle
from examo import *



################################################################################
# FUNCTIONS

#For importing the file, need to not feed in other string names in the file. 

def solvemax(rxn, md):
    try:
        thresh = 1E-10
        #md = importPickle(fOutModeldict)
        cbm = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
        cbm.initLp()
        #cbm.guro.setObjective(0)
        exec 'cbm.guro.setObjective(cbm.%s, GRB.MAXIMIZE)' % rxn
        cbm.guro.optimize()
        sol = abs(array([v.x for v in cbm.guro.getVars()]))
        indices = (sol > thresh).nonzero()[0]
        return indices
    except:
        print "exception 1"
        indices = []
        return indices

def solvemin(rxn, md):
    try:
        thresh = 1E-10
        #md = importPickle(fOutModeldict)
        cbm = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
        cbm.initLp()
        #cbm.guro.setObjective(0)
        exec 'cbm.guro.setObjective(cbm.%s, GRB.MINIMIZE)' % rxn
        cbm.guro.optimize()
        sol = abs(array([v.x for v in cbm.guro.getVars()]))
        indices = (sol > thresh).nonzero()[0]
        return indices
    except:
        print "exception 2"
        indices = []
        return indices

def solvemax_star(rxn_cbm):
    return solvemax(*rxn_cbm)

def solvemin_star(rxn_cbm):
    return solvemin(*rxn_cbm)

def findActiveRxns(cbm, thresh, mDict, rl = []):
    act = set()
    arrayIdRs = array(cbm.idRs[:])
    cbm.initLp()
    if rl:
        idRs = rl
    else:
        idRs = cbm.idRs[:]
    # maximizing all reactions at once
    # reseting the objective
    cbm.guro.setObjective(0)
    # setting the objective
    s = 'cbm.linobj = LinExpr([1.0] * len(cbm.idRs), ['
    for var in cbm.guro.getVars():
        s += 'cbm.%s, ' % var.varName
    s = s.rstrip(', ')
    s += '])'
    exec s
    #EG Initially set the objective to maximize
    cbm.guro.setObjective(cbm.linobj, 1)#1 for maximize
    cbm.guro.optimize()
    sol = abs(array([v.x for v in cbm.guro.getVars()]))
    indices = (sol > thresh).nonzero()[0]
    act.update(arrayIdRs[indices])
    idRs = list(set(idRs) - act)
    idRs_ub = []
    # maximizing
    #EG Reduce the number of reactions that need to be investigated based off of upper boundary constraints
    for i in idRs:
        if cbm.ub[cbm.idRs.index(i)] != 0:
	    idRs_ub.append(i)
    # reseting the objective
    workers = 7
    pool = Pool(processes = workers)
    results = pool.map(solvemax_star, itertools.izip(idRs, itertools.repeat(mDict)))
    pool.close()
    pool.join()
    for i in results:
	act.update(arrayIdRs[i])
    idRs = list(set(idRs) - act)
    idRs_lb = []
    # minimizing
    #EG Reduce the number of reactions that need to be investigated based off of lower boundary constraints
    for i in idRs:
	if cbm.lb[cbm.idRs.index(i)] != 0:
	    idRs_lb.append(i)
    # reseting the objective
    pool = Pool(processes = workers)
    results = pool.map(solvemin_star, itertools.izip(idRs, itertools.repeat(mDict)))
    pool.close()
    pool.join()
    for i in results:
        act.update(arrayIdRs[i])
    return act

def pruneRxn(cbm, cH, rxn, thresh, description, repetition):
    try:
        #EG Prune a reaction. If a flux soltuion cannot be obtained, or if the biomass flux becomes inactive, stop pruning.
        rxntodelete = rxn
        m0 = deleteCbmRxns(cbm, rxntodelete)
	fOutModel = 'examoModules/data/iMM904_examo_%s_%s.pkl'
	fOutModeldict = 'examoModules/data/iMM904_examo_test_dict.pkl'
	#exportPickle(m0, fOutModel % (description, os.getpid()))
	rev = [0 if val >= 0 else 1 for val in m0.lb]
	gene2rxn = {}
	rxns = {}
	for rxn in m0.rxns:
   	    gene2rxn[rxn] = m0.rxns[rxn]['genes']
    	    rxns[rxn] = m0.rxns[rxn]
	mDict = {
	'S' : m0.S,
	'idSp' : m0.idSp,
	'idRs' : m0.idRs,
	'lb' : m0.lb,
	'ub' : m0.ub,
	'rev' : rev,
	'genes' : m0.genes,
	'gene2rxn' : gene2rxn,
	'rxns' : rxns,
	'descrip' : 'iMM904'}
	#exportPickle(mDict, fOutModeldict)
        #NOTE the threshold for is set a bit higher for cH rxns
        act = findActiveRxns(m0, thresh, mDict, cH)
        actlist = set()
        for i in act: 
	    if i in cH:
	        actlist.add(i)

	#EG Beginning of the 4th script
	################################################################################
	# INPUTS
	eps = 1E-10
	activityThreshold = 1E-10
	fModelExamo = 'examoModules/data/iMM904_examo_test_dict.pkl'
	fFreqBasedRxns = 'data/freqBasedRxns_%s.pkl'
	################################################################################
	# STATEMENTS

	hfr = importPickle(fFreqBasedRxns % description)['hfr']
	md = importPickle(fModelExamo)
	mtry1 = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
	hfr = hfr & set(mtry1.idRs)
	#forcing biomass production
	mtry1.lb[mtry1.idRs.index('R_biomass_published')] = 0.2879
	#minimizingg the sum of fluxes
	mtry1result = MipSeparateFwdRev_gurobi(mtry1, hfr, eps)
	mtry1result.initMipGurobi()
	mtry1result.minSumFluxes_gurobi()
	#EG Added activityThreshold and the md['rxns'] dictionary to the function, so that the reactants and products could be written out
	nz = getNzRxnsGurobi(mtry1result, activityThreshold, md['rxns'])[1]
    except:
        print "exited early"
	return cbm
    if (len(cH - set(actlist)) != 0):#not all cH rxns are active
	print "exited 1"
	return cbm
    else:
	#EG Identify the reactions that became inactive after the reaction was delted. If extra deleted reactions cause the model to be unsolvable, or if extra deleted inactive reactions cause any of the hfrs to become inactive, or if a solution cannot be obtained with a biomass flux, only delete the one reaction. Otherwise, delete the inactive reactions.	
	try:
	    inact = set(m0.idRs) - act - cH	
	    m1 = deleteCbmRxns(m0, inact)	
	    fOutModel2 = 'examoModules/data/iMM904_examo_%s_%s2.pkl'
	    fOutModeldict2 = 'examoModules/data/iMM904_examo_test_dict2.pkl'
	    #exportPickle(m1, fOutModel2 % (description, os.getpid()))
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
	    'descrip' : 'iMM904'}
	    #exportPickle(mDict, fOutModeldict2)
	    act2 = findActiveRxns(m1, thresh, mDict, cH)
	    actlist2 = set()
	    for j in act2:
 	        if j in cH:
		    actlist2.add(j)

	    #EG Beginning of the 4th script
	    ################################################################################
	    # INPUTS
	    eps = 1E-10
	    activityThreshold = 1E-10

	    fModelExamo2 = 'examoModules/data/iMM904_examo_test_dict2.pkl'
	    fFreqBasedRxns = 'data/freqBasedRxns_%s.pkl'

	    ################################################################################
	    # STATEMENTS

	    hfr = importPickle(fFreqBasedRxns % description)['hfr']
	    md = importPickle(fModelExamo2)
	    mtry2 = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
	    hfr = hfr & set(mtry2.idRs)
	    #forcing biomass production
	    mtry2.lb[mtry2.idRs.index('R_biomass_published')] = 0.2879
	    #minimizing the sum of fluxes
	    mtry2result = MipSeparateFwdRev_gurobi(mtry2, hfr, eps)
	    mtry2result.initMipGurobi()
	    mtry2result.minSumFluxes_gurobi()
	    #EG Added activityThreshold and the md['rxns'] dictionary to the function, so that the reactants and products could be written out
	    nz = getNzRxnsGurobi(mtry2result, activityThreshold, md['rxns'])[1]
            print inact
	    return m1	
	except:
            print "exited 2"
	    m0 = deleteCbmRxns(cbm, rxntodelete)
	    return m0

	if (len(cH - set(actlist2)) != 0):#not all cH rxns are active
	   m0 = deleteCbmRxns(cbm, rxntodelete)
	   return m0

#EG 131112 Avoided creating sets for prunableRxns so that randomness would be preserverd, and first try pruning transport reactions before other reactions in the model
def iterativePrunning(i, m, cH, description, repetition, thresh = 1E-10, EXrxns = set(), EXtrrxns = set(), Othertrrxns = set()):
    """
    solver can be 'cplex', 'glpk' or 'gurobi'
    """
    semilla = int((time.time() * 1E6) * os.getpid())
    random.seed(semilla)
    EXrxnsprune = list(set(list(EXrxns)) - cH)
    for j in range(i+1):
	random.shuffle(EXrxnsprune)
    while EXrxnsprune:
	rxn1 = EXrxnsprune.pop()
        print rxn1
        print len(EXrxnsprune)
	try:
            mTemp1 = pruneRxn(mTemp1, cH, rxn1, thresh, description, repetition)
	    EXrxnsprune2 = []
	    for k in mTemp1.idRs:
		if k in EXrxnsprune:
			EXrxnsprune2.append(k)			
	    for l in range(len(EXrxnsprune2)):
		random.shuffle(EXrxnsprune2)
		EXrxnsprune = EXrxnsprune2
        except NameError:
            mTemp1 = pruneRxn(m, cH, rxn1, thresh, description, repetition)
            EXrxnsprune2 = []
	    for k in mTemp1.idRs:
		if k in EXrxnsprune:
			EXrxnsprune2.append(k)
	    for l in range(len(EXrxnsprune2)):
		random.shuffle(EXrxnsprune2)
		EXrxnsprune = EXrxnsprune2
    EXtrrxnsprune = list(set(list(EXtrrxns)) - cH)
    EXtrrxnsprunelist = []
    for j in EXtrrxnsprune:
	if j in mTemp1.idRs:
	    EXtrrxnsprunelist.append(j)
    for j in range(i+1):
	random.shuffle(EXtrrxnsprune)
    while EXtrrxnsprune:
	rxn2 = EXtrrxnsprune.pop()
        print rxn2
        print len(EXtrrxnsprune)
        mTemp1 = pruneRxn(mTemp1, cH, rxn2, thresh, description, repetition)
	EXtrrxnsprune2 = []
	for k in mTemp1.idRs:
	    if k in EXtrrxnsprune:
		EXtrrxnsprune2.append(k)			
	for l in range(len(EXtrrxnsprune2)):
	    random.shuffle(EXtrrxnsprune2)
	    EXtrrxnsprune = EXtrrxnsprune2
    prunableRxns = []
    for j in mTemp1.idRs:
	if j not in list(cH):
	    if j not in EXrxns:
	        if j not in EXtrrxns:
		    if j not in Othertrrxns:
		        prunableRxns.append(j)
    semilla = int((time.time() * 1E6) * os.getpid())
    random.seed(semilla)
    for j in range(i+1):
	random.shuffle(prunableRxns)
    while prunableRxns:
        rxn3 = prunableRxns.pop()
        print rxn3
	print len(prunableRxns)
        mTemp1 = pruneRxn(mTemp1, cH, rxn3, thresh, description, repetition)
	prunableRxns2 = []
	for k in mTemp1.idRs:
	    if k in prunableRxns:
		prunableRxns2.append(k)			
	for l in range(len(prunableRxns2)):
	    random.shuffle(prunableRxns2)
	    prunableRxns = prunableRxns2
    return mTemp1.idRs

if __name__ == '__main__':
    ########################################
    # INPUTS
    
    md = importPickle('../data/iMM904_blkRxnsDeleted_dict.pkl')
    rH = importPickle('../data/rxnsClassifiedByExprssion_eth_15_85.pkl')['rH']
    fOutMbaCandRxns = '../data/mbaCandRxns/mbaCandRxns_%s.pkl'

    activityThreshold = 1E-10
    descrip = 'eth_15_85'
    ########################################
    # STATEMENTS
    m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'],
        md['genes'])
    #act = findActiveRxns(m, activityThreshold, m.idRs)
    #mT = pruneRxn(m, rH, 'R_PFK', activityThreshold)
    #iterativePrunning(m, rH, activityThreshold)
    
    numProc = 2
    numRep = 2
    # making sure that all rH reactions are active to begin with
    act = findActiveRxns(m, 1E-10, rH)
    cH = rH & act
    cH.add('R_biomass_published')
    #print '%i processes and %i repetitions' % (numProc, numRep)
    for i in range(numProc):
        pid = os.fork()
        if pid == 0:
            for j in range(numRep):
                locTime = time.localtime()
                pid = os.getpid()
                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
                tag = '%s_%s_%s' % (descrip, pid, timeStr)
                cr = iterativePrunning(m, cH, activityThreshold, m.idRs[:10])
                exportPickle(cr, fOutMbaCandRxns % tag)
            os._exit(0)




