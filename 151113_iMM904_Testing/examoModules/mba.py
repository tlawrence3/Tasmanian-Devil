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
import time
import os
import shutil
import sys
import re

from sys import path; path.append('./modules/')
from utilities import importPickle, exportPickle
from examo import *
from decimal import Decimal, getcontext, ROUND_DOWN
from Crypto.Random import random, atfork



################################################################################
# FUNCTIONS

#@profile
def findActiveRxns(cbm, thresh, rl = []):
    act = set()
    arrayIdRs = array(cbm.idRs[:])
    init = cbm.initLp()
    init 
    if rl:
        idRs = rl
    else:
        idRs = cbm.idRs[:]
    # maximizing all reactions at once
    # reseting the objective
    #cbm.guro.setObjective(0)
    # setting the objective
    #s = 'cbm.linobj = LinExpr([1.0] * len(cbm.idRs), ['
    #for var in cbm.guro.getVars():
    #    s += 'cbm.%s, ' % var.varName
    #s = s.rstrip(', ')
    #s += '])'
    #exec s
    #EG Initially set the objective to maximize
    #cbm.guro.setObjective(cbm.linobj)
    #cbm.guro.optimize()
    #sol = abs(array([v.x for v in cbm.guro.getVars()]))
    #indices = (sol > thresh).nonzero()[0]
    #act.update(arrayIdRs[indices])
    #idRs = list(set(idRs) - act)
    #print len(idRs)
    # maximizing
    for rxn in idRs:
        if rxn not in act:
            # reseting the objective
            cbm.guro.setObjective(0)
            exec 'cbm.guro.setObjective(cbm.%s, GRB.MAXIMIZE)' % rxn
            cbm.guro.optimize()
            sol = abs(array([v.x for v in cbm.guro.getVars()]))
            indices = (sol > thresh).nonzero()[0]
            act.update(arrayIdRs[indices])
    idRs = list(set(idRs) - act)
    # minimizing
    for rxn in idRs:
        if rxn not in act:
            # reseting the objective
            cbm.guro.setObjective(0)
            exec 'cbm.guro.setObjective(cbm.%s, GRB.MINIMIZE)' % rxn
            cbm.guro.optimize()
            sol = abs(array([v.x for v in cbm.guro.getVars()]))
            indices = (sol > thresh).nonzero()[0]
            act.update(arrayIdRs[indices])
    return act

#@profile
def pruneRxn(cbm, cH, rxn, thresh, eps, activityThreshold, description, pickle_model_name, 
             repetition, biomassRxn, lb_biomass):
    try:
        #EG Prune a reaction. If a flux solution cannot be obtained
        #or if the biomass flux becomes inactive, stop pruning.
        rxntodelete = rxn
        m0 = deleteCbmRxns(cbm, rxntodelete)
        act = findActiveRxns(m0, thresh)
        cH_act = cH & act
        if (len(cH - cH_act) != 0):#not all cH rxns are active
            #print "not all active 1"
            return cbm
        #######################################################################
        # INPUTS
        fFreqBasedRxns = 'data/freqBasedRxns_%s_%s.pkl'
        #######################################################################
        # STATEMENTS
        hfr = importPickle(fFreqBasedRxns % (description, pickle_model_name))['hfr']
        hfr = hfr & set(m0.idRs)
        #forcing biomass production
        m0.lb[m0.idRs.index(biomassRxn)] = lb_biomass
        #minimizing the sum of fluxes
        mtry1result = MipSeparateFwdRev_gurobi(m0, hfr, eps)
        mtry1result.initMipGurobi()
        mtry1result.minSumFluxes_gurobi()
        #EG Added activityThreshold and the m0.rxns dictionary to the
        #function, so that the reactants and products could be written out
        #print "hello 8"
        nz = getNzRxnsGurobi(mtry1result, activityThreshold, m0.rxns)[1]
        #print "hello 9"
    except:
        #print "exception 1"
        return cbm
        #EG Identify the reactions that became inactive after the
        #reaction was deleted. If extra deleted reactions cause the
        #model to be unsolvable, or if extra deleted inactive reactions
        #cause any of the hfrs to become inactive, or if a solution
        #cannot be obtained with a biomass flux, only delete the one
        #reaction. Otherwise, delete the inactive reactions.
    try:
        inact = set(m0.idRs) - act - cH
        m1 = deleteCbmRxns(m0, inact)
        act2 = findActiveRxns(m1, thresh)
        cH_act2 = cH & act2
        if (len(cH - cH_act2) != 0):#not all cH rxns are active
            #print rxntodelete
            return m0
        # STATEMENTS
        hfr = importPickle(fFreqBasedRxns % (description, pickle_model_name))['hfr']
        hfr = hfr & set(m1.idRs)
        #forcing biomass production
        m1.lb[m1.idRs.index(biomassRxn)] = lb_biomass
        #minimizing the sum of fluxes
        mtry2result = MipSeparateFwdRev_gurobi(m1, hfr, eps)
        mtry2result.initMipGurobi()
        mtry2result.minSumFluxes_gurobi()
        #EG Added activityThreshold and the m1.rxns dictionary
        #to the function, so that the reactants and products could
        #be written out
        nz = getNzRxnsGurobi(mtry2result, activityThreshold, m1.rxns)[1]
        #print inact
        return m1
    except:
        #print "exception 2"
        return m0

#EG 131112 Avoided creating sets for prunableRxns so that randomness
#would be preserverd, and first try pruning transport reactions before
#other reactions in the model
#@profile
def iterativePrunning(i, m, cH, description, pickle_model_name, biomassRxn, lb_biomass,
                      repetition, thresh = 1E-10, eps = 1E-10, activityThreshold = 1E-10, EXrxns = [],
                      EXtrrxns = [], Othertrrxns = []):
    """
    solver can be 'cplex', 'glpk' or 'gurobi'
    """
    if len(EXrxns) > 0:
        EXrxnsprune = list(set(list(EXrxns)) - cH)
        random.shuffle(EXrxnsprune)
        while EXrxnsprune:
            rxn1 = EXrxnsprune.pop()
            #print rxn1
            try:
                mTemp1 = pruneRxn(mTemp1, cH, rxn1, thresh, eps, activityThreshold, description,
                                  pickle_model_name, repetition, biomassRxn, lb_biomass)
                EXrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXrxnsprune:
                        EXrxnsprune2.append(k)
                random.shuffle(EXrxnsprune2)
                EXrxnsprune = EXrxnsprune2
            except NameError:
                mTemp1 = pruneRxn(m, cH, rxn1, thresh, eps, activityThreshold, description,
                                  pickle_model_name, repetition, biomassRxn, lb_biomass)
                EXrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXrxnsprune:
                        EXrxnsprune2.append(k)
                random.shuffle(EXrxnsprune2)
                EXrxnsprune = EXrxnsprune2
    if len(EXtrrxns) > 0:
        EXtrrxnsprune = list(set(list(EXtrrxns)) - cH)
        EXtrrxnsprunelist = []
        for j in EXtrrxnsprune:
            if j in mTemp1.idRs:
                EXtrrxnsprunelist.append(j)
        random.shuffle(EXtrrxnsprune)
        while EXtrrxnsprune:
            rxn2 = EXtrrxnsprune.pop()
            #print rxn2
            try:
                mTemp1 = pruneRxn(mTemp1, cH, rxn2, thresh, eps, activityThreshold, description,
                                  pickle_model_name, repetition, biomassRxn, lb_biomass)
                EXtrrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXtrrxnsprune:
                        EXtrrxnsprune2.append(k)
                random.shuffle(EXtrrxnsprune2)
                EXtrrxnsprune = EXtrrxnsprune2
            except NameError:
                mTemp1 = pruneRxn(m, cH, rxn2, thresh, eps, activityThreshold, description,
                                  pickle_model_name, repetition, biomassRxn, lb_biomass)
                EXtrrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXtrrxnsprune:
                        EXtrrxnsprune2.append(k)
                random.shuffle(EXtrrxnsprune2)
                EXtrrxnsprune = EXtrrxnsprune2

    prunableRxns = []
    try:
        for j in mTemp1.idRs:
            if j not in list(cH):
                if j not in EXrxns:
                    if j not in EXtrrxns:
                        if j not in Othertrrxns:
                            prunableRxns.append(j)
    except NameError:
        for j in m.idRs:
            if j not in list(cH):
                if j not in EXrxns:
                    if j not in EXtrrxns:
                        if j not in Othertrrxns:
                            prunableRxns.append(j)
    random.shuffle(prunableRxns)
    while prunableRxns:
        rxn3 = prunableRxns.pop()
        try:
            mTemp1 = pruneRxn(mTemp1, cH, rxn3, thresh, eps, activityThreshold, description,
                              pickle_model_name, repetition, biomassRxn, lb_biomass)
            prunableRxns2 = []
            for k in mTemp1.idRs:
                if k in prunableRxns:
                    prunableRxns2.append(k)
            random.shuffle(prunableRxns2)
            prunableRxns = prunableRxns2
        except NameError:
            mTemp1 = pruneRxn(m, cH, rxn3, thresh, eps, activityThreshold, description,
                              pickle_model_name, repetition, biomassRxn, lb_biomass)
            prunableRxns2 = []
            for k in mTemp1.idRs:
                if k in prunableRxns:
                    prunableRxns2.append(k)
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
    print '%i processes and %i repetitions' % (numProc, numRep)
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




