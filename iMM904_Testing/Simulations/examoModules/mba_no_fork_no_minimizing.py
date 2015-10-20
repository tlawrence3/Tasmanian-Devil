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
import random
import time
import os
import shutil
import sys
import re

from sys import path; path.append('./modules/')
from utilities import importPickle, exportPickle
from examo import *
from decimal import Decimal, getcontext, ROUND_DOWN
import cobra



######################################################################
# FUNCTIONS

def findActivefromZfr(cbm, thresh, rl = []):
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
    return act

@profile
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
    cbm.guro.setObjective(0)
    # setting the objective
    s = 'cbm.linobj = LinExpr([1.0] * len(cbm.idRs), ['
    for var in cbm.guro.getVars():
        s += 'cbm.%s, ' % var.varName
    s = s.rstrip(', ')
    s += '])'
    exec s
    #EG Initially set the objective to maximize
    cbm.guro.setObjective(cbm.linobj)#1 for maximize
    cbm.guro.optimize()
    sol = abs(array([v.x for v in cbm.guro.getVars()]))
    indices = (sol > thresh).nonzero()[0]
    act.update(arrayIdRs[indices])
    idRs = list(set(idRs) - act)
    # maximizing
    for rxn in idRs:
        if rxn not in act:
            # reseting the objective
            cbm.guro.setObjective(0)
            exec 'cbm.guro.setObjective(cbm.%s, GRB.MAXIMIZE)' % rxn
            cbm.guro.optimize()
            sol = abs(array([v.x for v in cbm.guro.getVars()]))
            indices = (sol > thresh).nonzero()[0]
            #act2 = act.copy()
            act.update(arrayIdRs[indices])
            #if act2 != act:
                #print "added 1"
                #if rxn in act:
                    #print rxn
    return act

@profile
def pruneRxn(cbm, cH, rxn, thresh, description, repetition, biomassRxn,
             lb_biomass):
    try:
        #EG Prune a reaction. If a flux soltuion cannot be obtained
        #or if the biomass flux becomes inactive, stop pruning.
        rxntodelete = rxn
        m0 = deleteCbmRxns(cbm, rxntodelete)
        #NOTE the threshold for is set a bit higher for cH rxns
        act = findActiveRxns(m0, thresh, cH)
        cH_act = cH & act
        if (len(cH - cH_act) != 0):#not all cH rxns are active
            print "not all active 1"
            return cbm
        #######################################################################
        # INPUTS
        eps = 1E-10
        activityThreshold = 1E-10
        fFreqBasedRxns = '../data/freqBasedRxns_%s.pkl'
        #######################################################################
        # STATEMENTS
        hfr = importPickle(fFreqBasedRxns % description)['hfr']
        hfr = hfr & set(m0.idRs)
        #forcing biomass production
        m0.lb[m0.idRs.index(biomassRxn)] = lb_biomass
        #minimizingg the sum of fluxes
        mtry1result = MipSeparateFwdRev_gurobi(m0, hfr, eps)
        mtry1result.initMipGurobi()
        mtry1result.minSumFluxes_gurobi()
        #EG Added activityThreshold and the m0.rxns dictionary to the
        #function, so that the reactants and products could be written out
        nz = getNzRxnsGurobi(mtry1result, activityThreshold, m0.rxns)[1]
    except:
        print "exception 1"
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
        act2 = findActiveRxns(m1, thresh, cH)
        cH_act2 = cH & act2
        if (len(cH - cH_act2) != 0):#not all cH rxns are active
            print rxntodelete
            return m0
        ###################################################################
        # INPUTS
        eps = 1E-10
        activityThreshold = 1E-10
        fFreqBasedRxns = '../data/freqBasedRxns_%s.pkl'
        ###################################################################
        # STATEMENTS
        hfr = importPickle(fFreqBasedRxns % description)['hfr']
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
        print inact
        return m1
    except:
        print "exception 2"
        return m0

#EG 131112 Avoided creating sets for prunableRxns so that randomness
#would be preserverd, and first try pruning transport reactions before
#other reactions in the model
@profile
def iterativePrunning(i, m, cH, description, biomassRxn, lb_biomass,
                      repetition, thresh = 1E-10, EXrxns = [],
                      EXtrrxns = [], Othertrrxns = []):
    """
    solver can be 'cplex', 'glpk' or 'gurobi'
    """
    if len(EXrxns) > 0:
        EXrxnsprune = list(set(list(EXrxns)) - cH)
        random.shuffle(EXrxnsprune)
        while EXrxnsprune:
            rxn1 = EXrxnsprune.pop()
            try:
                mTemp1 = pruneRxn(mTemp1, cH, rxn1, thresh, description,
                                  repetition, biomassRxn, lb_biomass)
                EXrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXrxnsprune:
                        EXrxnsprune2.append(k)
                random.shuffle(EXrxnsprune2)
                EXrxnsprune = EXrxnsprune2
            except NameError:
                mTemp1 = pruneRxn(m, cH, rxn1, thresh, description,
                                  repetition, biomassRxn, lb_biomass)
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
            try:
                mTemp1 = pruneRxn(mTemp1, cH, rxn2, thresh, description,
                                  repetition, biomassRxn, lb_biomass)
                EXtrrxnsprune2 = []
                for k in mTemp1.idRs:
                    if k in EXtrrxnsprune:
                        EXtrrxnsprune2.append(k)
                random.shuffle(EXtrrxnsprune2)
                EXtrrxnsprune = EXtrrxnsprune2
            except NameError:
                mTemp1 = pruneRxn(m, cH, rxn2, thresh, description,
                                  repetition, biomassRxn, lb_biomass)
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
            mTemp1 = pruneRxn(mTemp1, cH, rxn3, thresh, description,
                              repetition, biomassRxn, lb_biomass)
            prunableRxns2 = []
            for k in mTemp1.idRs:
                if k in prunableRxns:
                    prunableRxns2.append(k)
            random.shuffle(prunableRxns2)
            prunableRxns = prunableRxns2
            #prunableRxns2_appended = prunableRxns.append()
            #prunableRxns2_appended k for k in mTemp1.idRs if k in prunableRnxs
            #for k in mTemp1.idRs:
            #    if k in prunableRxns:
            #        prunableRxns2.append(k)
            #random.shuffle(prunableRxns2_appended)
            #prunableRxns = prunableRxns2_appended
        except NameError:
            mTemp1 = pruneRxn(m, cH, rxn3, thresh, description,
                              repetition, biomassRxn, lb_biomass)
            prunableRxns2 = []
            for k in mTemp1.idRs:
                if k in prunableRxns:
                    prunableRxns2.append(k)
            random.shuffle(prunableRxns2)
            prunableRxns = prunableRxns2
            #prunableRxns2_appended = prunableRxns.append()
            #prunableRxns2_appended k for k in mTemp1.idRs if k in prunableRnxs
            #random.shuffle(prunableRxns2_appended)
            #prunableRxns = prunableRxns2_appended
            #for k in mTemp1.idRs:
            #    if k in prunableRxns:
            #        prunableRxns2.append(k)
            #random.shuffle(prunableRxns2)
            #prunableRxns = prunableRxns2
    return mTemp1.idRs

if __name__ == '__main__':
#    ########################################
#    # INPUTS
#    md = importPickle('../data/iMM904_blkRxnsDeleted_dict.pkl')
#    rH = importPickle('../data/rxnsClassifiedByExprssion_eth_15_85.pkl')['rH']
#    fOutMbaCandRxns = '../data/mbaCandRxns/mbaCandRxns_%s.pkl'
#
#    activityThreshold = 1E-10
#    descrip = 'eth_15_85'
#    ########################################
#    # STATEMENTS
#    m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'],
#        md['genes'])
#    #act = findActiveRxns(m, activityThreshold, m.idRs)
#    #mT = pruneRxn(m, rH, 'R_PFK', activityThreshold)
#    #iterativePrunning(m, rH, activityThreshold)
#    numProc = 2
#    numRep = 2
#    # making sure that all rH reactions are active to begin with
#    act = findActiveRxns(m, 1E-10, rH)
#    cH = rH & act
#    cH.add('R_biomass_published')
#    #print '%i processes and %i repetitions' % (numProc, numRep)
#    for i in range(numProc):
#        pid = os.fork()
#        if pid == 0:
#            for j in range(numRep):
#                locTime = time.localtime()
#                pid = os.getpid()
#                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
#                tag = '%s_%s_%s' % (descrip, pid, timeStr)
#                cr = iterativePrunning(m, cH, activityThreshold, m.idRs[:10])
#                exportPickle(cr, fOutMbaCandRxns % tag)
#            os._exit(0)


    ########################################
    # INPUTS
    model = sys.argv[1]
    pickle_model = sys.argv[2]
    description = sys.argv[3]
    biomassRxn = sys.argv[4]
    repetitions = int(sys.argv[5])
    #Create necessary variables and import the model
    if model[-4:] == '.xml':
        cobra_model = cobra.io.read_sbml_model('../data/%s' % model)
    if model[-4:] == '.mat':
        cobra_model = cobra.io.mat.load_matlab_model('../data/%s' % model)
    # model dictionary of the original model with blocked reactions deleted
    fModelDict = '../data/%s' % pickle_model

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
        fbr = importPickle('../data/freqBasedRxns_%s.pkl' % description)

        #EG Making subdirectories for candidate reactions
        mbaCandRxnsDirectory = '../data/mbaCandRxns/%s_%s/' % (description, str(repetition))
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
                    act = findActiveRnxs(m1, 1E-10, new_hfr)
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
        mbaCandRxnsDirectorySubset = '../data/%s_%s/' % (description, str(repetition))
        if not os.path.exists(mbaCandRxnsDirectorySubset):
            os.mkdir(mbaCandRxnsDirectorySubset, 0777)
        #Run the MBA
        for i in range(numProc):
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
