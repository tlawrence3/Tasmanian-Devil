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

from gurobipy import *
from scipy.sparse import coo_matrix


class CbModel():
    """
    Constraint-based model class
    """
    #{{{
    def __init__(self, S, idSp, idRs, lb, ub, rxns = {}, genes = set(),
                 descrip = ''):
        """
        S should be a sparse coo_matrix
        """
        #{{{2
        self.S = S
        self.idSp = idSp
        self.idRs = idRs
        self.lb = lb
        self.ub = ub
        self.rxns = rxns
        self.descrip = descrip
        self.genes = genes
        self.gene2rxn = {}
        for rxn in self.rxns:
            self.gene2rxn[rxn] = self.rxns[rxn]['genes']
        #}}}2

    def separateFwdRevRxns(self):
        """
        Returns
        modelSep [CbModel instance] reversible reactions separated into 
            fwd and rev
        """
        #{{{
        #from cbModel import CbModel
        #from scipy.sparse import coo_matrix
        from numpy import concatenate
        indRev = [i for i, rxn in enumerate(self.idRs) if (self.lb[i] < 0.) 
                and (self.ub[i] > 0.)]
        indIrrev = [i for i in range(self.S.shape[1]) if i not in indRev]
        # augmenting the stoichiometric matrix
        S = self.S.toarray()
        Srev = -S[:, indRev]
        Saug = concatenate((S, Srev), axis = 1)
        # updating idRs, lb and ub
        idRs = self.idRs + [rxn + '_rev' for i, rxn in enumerate(self.idRs) 
                if i in indRev]
        lb = [0. if val <= 0. else val for val in self.lb] + [0.]*len(indRev)
        ub = self.ub + [abs(val) for i, val in enumerate(self.lb) 
                if i in indRev]
        cbmNew = CbModel(coo_matrix(Saug), self.idSp, idRs, lb, ub,
                descrip = self.descrip + '_fwdAndRevSeparated')
        #species dictionary (unchanged)
        if 'species' in dir(self):
            cbmNew.species = self.species
        #genes set (unchanged)
        if 'genes' in dir(self):
            cbmNew.genes = self.genes
        #reactions dictionary
        if 'rxns' in dir(self):
            rxns = {}
            for rxn in cbmNew.idRs:
                if '_rev' in rxn[-4:]:
                    rxns[rxn] = {}
                    rxns[rxn].update(self.rxns[rxn[:-4]])
                    rxns[rxn]['id'] = rxn
                    rs = rxns[rxn]['reactants']
                    ps = rxns[rxn]['products']
                    rxns[rxn]['reactants'] = ps
                    rxns[rxn]['products'] = rs
                else:
                    rxns[rxn] = self.rxns[rxn]
            cbmNew.rxns = rxns
        return cbmNew 
        #}}}

    def initLp(self, name = 'unnamed'):
        from numpy import array
        self.guro = Model(name)
        #turning off the writing of the gurobi.log file
        self.guro.setParam('OutputFlag', 0) 
        for i, rxn in enumerate(self.idRs):
            exec 'self.{} = self.guro.addVar(lb = {}, ub = {}, vtype = GRB.CONTINUOUS, name = "{}")'.format(rxn, self.lb[i], self.ub[i], rxn)
        self.guro.update()
        # adding constraints
        for i, row in enumerate(self.S.toarray()):
            nz = row.nonzero()[0]
            pair = zip(row[nz], array(self.idRs)[nz])
            s = ''
            for p in pair:
                s += '({} * self.{}) + '.format(p[0], p[1])
            s = s.rstrip(' + ')
            s += ' == 0.'
            exec 'self.guro.addConstr( {}, "{}")'.format(s, self.idSp[i])

    def setObjective(self, obj, optSense = 'max'):
        s = ''
        s += 'self.{}'.format(obj)
        if optSense == 'max':
            exec 'self.guro.setObjective({}, GRB.MAXIMIZE)'.format(s)
        elif optSense == 'min':
            exec 'self.guro.setObjective({}, GRB.MINIMIZE)'.format(s)
    
    def solveLp(self):
        self.guro.optimize()
        if self.guro.getAttr('status') == 2:
            self.solObjective = self.guro.objVal
            self.solution = [v.x for v in self.guro.getVars()]
        else:
            self.solObjective = None
            self.solution = []

    def fva(self, rl = []):
        """
        Performs flux variability analysis (FVA)
        """
        if not rl:
            rl = self.idRs
        fv = []
        fvMax = []
        # maximizing each
        for rxn in rl:
            # reseting the objective
            self.guro.setObjective(0)
            exec 'self.guro.setObjective(self.{}, GRB.MAXIMIZE)'.format(rxn)
            self.guro.optimize()
            fvMax.append(self.guro.objVal)
        # minimizing each
        fvMin = []
        for rxn in rl:
            # reseting the objective
            self.guro.setObjective(0)
            exec 'self.guro.setObjective(self.{}, GRB.MINIMIZE)'.format(rxn)
            self.guro.optimize()
            fvMin.append(self.guro.objVal)
        # getting the lists together
        for i, rxn in enumerate(rl):
            fv.append([rxn, fvMin[i], fvMax[i]])
        return fv

    #}}}


class MetabGeneExpModel(object):
    #{{{
    def __init__(self, idSp, idRs, S, lb, ub, rH, rL, eps = 1E-10):
        """
        idSp [list] strings with the names of species
        idRs [list] strings with the names of reactions
        S [array or coo_matrix] stoichiometry matrix
        lb [list] floats with lower bound values
        ub [list] floats with upper bound values
        rH [list] names of reactions classified as highly expressed
        rL [list] names of reactions classified as lowly expressed
        """
        #{{{2
        self.S = coo_matrix(S)
        self.Srows = [int(i) for i in self.S.row]
        self.Scols = [int(i) for i in self.S.col]
        self.Svals = [float(i) for i in self.S.data]
        self.idSp = idSp
        self.idRs = idRs
        self.lb = lb
        self.lb = lb
        self.ub = ub
        self.eps = eps
        #copying the S for the construction of the lhs
        self.lhsRows = self.Srows[:]
        self.lhsCols = self.Scols[:]
        self.lhsVals = self.Svals[:]
        self.rowNames = self.idSp[:]
        self.colNames = self.idRs[:]
        self.rH = rH
        self.rL = rL
        self.rHrev = [rxn for i, rxn in enumerate(self.idRs) if 
                (self.lb[i] < 0.) and (self.ub[i] > 0.) and (rxn in self.rH)]
        #}}}2

    def rowAndColNames(self):
        """
        Adds names for the integer variables and constraints
        """
        #{{{2
        for rxn in self.rH:
            self.colNames.append('yp_{}'.format(rxn))
            self.rowNames.append('cp_{}'.format(rxn))
            if rxn in self.rHrev:
                self.colNames.append('ym_{}'.format(rxn))
                self.rowNames.append('cm_{}'.format(rxn))
                #yp + ym <= 1
                self.rowNames.append('cYsum_{}'.format(rxn))        
        for rxn in self.rL:
            self.colNames.append('y_{}'.format(rxn))
            self.rowNames.append('cl_{}'.format(rxn))
            self.rowNames.append('cu_{}'.format(rxn))
        #}}}2

    def buildLhsMatrix(self): 
        #{{{2
        for rxn in self.rH:
            self.lhsRows.append(self.rowNames.index('cp_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('yp_{}'.format(rxn)))
            self.lhsVals.append(self.lb[self.colNames.index(rxn)] - self.eps)   
            self.lhsRows.append(self.rowNames.index('cp_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn))
            self.lhsVals.append(1)
            if rxn in self.rHrev: 
                self.lhsRows.append(self.rowNames.index('cm_{}'.format(rxn)))
                self.lhsCols.append(self.colNames.index('ym_{}'.format(rxn)))
                self.lhsVals.append(self.ub[self.colNames.index(rxn)] + self.eps)
                self.lhsRows.append(self.rowNames.index('cm_{}'.format(rxn)))
                self.lhsCols.append(self.colNames.index(rxn))
                self.lhsVals.append(1)
                # yp + ym <= 1
                self.lhsRows.append(self.rowNames.index('cYsum_{}'.format(rxn)))
                self.lhsCols.append(self.colNames.index('yp_{}'.format(rxn)))
                self.lhsVals.append(1)
                self.lhsRows.append(self.rowNames.index('cYsum_{}'.format(rxn)))
                self.lhsCols.append(self.colNames.index('ym_{}'.format(rxn)))
                self.lhsVals.append(1)
        
        for rxn in self.rL:
            self.lhsRows.append(self.rowNames.index('cl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(self.lb[self.colNames.index(rxn)])
            self.lhsRows.append(self.rowNames.index('cl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn))
            self.lhsVals.append(1)
        
            self.lhsRows.append(self.rowNames.index('cu_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(self.ub[self.colNames.index(rxn)])
            self.lhsRows.append(self.rowNames.index('cu_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn))
            self.lhsVals.append(1)
        #}}}2

    def buildRhsAndSenses(self):
        #{{{2
        self.senses = 'E'*len(self.idSp)
        self.rhsRows = []
        self.rhsCols = []
        self.rhsVals = []
        for rxn in self.rH:
            self.rhsRows.append(self.rowNames.index('cp_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(self.lb[self.colNames.index(rxn)])
            self.senses += 'G'
            if rxn in self.rHrev: 
                self.rhsRows.append(self.rowNames.index('cm_{}'.format(rxn)))
                self.rhsCols.append(0)
                self.rhsVals.append(self.ub[self.colNames.index(rxn)])
                self.senses += 'L'        
                # yp + ym <= 1
                #NOTE cplex's L is less or equal, not strictly less!
                self.rhsRows.append(self.rowNames.index('cYsum_{}'.format(rxn)))
                self.rhsCols.append(0)
                self.rhsVals.append(1)
                self.senses += 'L'

        for rxn in self.rL:
            self.rhsRows.append(self.rowNames.index('cl_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(self.lb[self.colNames.index(rxn)])
            self.senses += 'G'
        
            self.rhsRows.append(self.rowNames.index('cu_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(self.ub[self.colNames.index(rxn)])
            self.senses += 'L'
        self.rhs = list(coo_matrix((self.rhsVals, (self.rhsRows, self.
            rhsCols)), shape = (len(self.rowNames), 1)).toarray().flatten())
        #}}}2
    #}}}

class MetabGeneExpModel_gurobi(MetabGeneExpModel):
    """
    Uses gurobi as its milp solver
    """
    #{{{
    def initGurobi(self):
        from numpy import array
        self.guro = Model('imge')
        #turning off the writing of the gurobi.log file
        self.guro.setParam('OutputFlag', 0) 
        # Initializing variables
        for i, rxn in enumerate(self.idRs):
            exec 'self.{} = self.guro.addVar(lb = {}, ub = {}, vtype = GRB.CONTINUOUS, name = "{}")'.format(rxn, self.lb[i], self.ub[i], rxn)
        for var in [v for v in self.colNames if v not in self.idRs]:
            exec 'self.{0} = self.guro.addVar(vtype = GRB.BINARY, name = "{0}")'.format(var)
        self.guro.update()
        # adding constraints
        lhs = coo_matrix((self.lhsVals, (self.lhsRows, self.lhsCols))).toarray()
        for i, row in enumerate(lhs):
            nz = row.nonzero()[0]
            pair = zip(row[nz], array(self.colNames)[nz])
            s = ''
            for p in pair:
                s += '({} * self.{}) + '.format(p[0], p[1]) 
            s = s.rstrip(' + ')
            if self.senses[i] == 'E':
                s += ' == {}'.format(self.rhs[i])
            elif self.senses[i] == 'G':
                s += ' >= {}'.format(self.rhs[i])
            elif self.senses[i] == 'L':
                s += ' <= {}'.format(self.rhs[i])
            exec 'self.guro.addConstr( {}, "{}")'.format(s, self.rowNames[i])
        #setting the objective
        s = ''
        for name in self.colNames:
            if name.startswith('y'):
                s += 'self.{} + '.format(name)
        s = s.rstrip(' + ')
        exec 'self.guro.setObjective({}, GRB.MAXIMIZE)'.format(s)

    def solve(self):
        self.rowAndColNames()
        self.buildLhsMatrix()
        self.buildRhsAndSenses()
        self.initGurobi()
        self.guro.optimize()
        self.initialized = 1

    def modulateEnzymeActivityAndDirection(self, rxn):
        """
        """
        #{{{2
        try:
            dummy =  self.initialized
        except AttributeError:
            self.solve()
        # 1. forcing rxn to be inactive reaction
        # The original bounds are saved in self.d before changing them
        self.d = {}
        exec 'self.d["{0}.lb"] = self.{0}.lb'.format(rxn)
        exec 'self.d["{0}.ub"] = self.{0}.ub'.format(rxn)
        exec "self.{}.setAttr('lb', 0.)".format(rxn)
        exec "self.{}.setAttr('ub', 0.)".format(rxn)
        if rxn in self.rH:
            exec 'self.d["yp_{0}.ub"] = self.yp_{0}.ub'.format(rxn)
            exec "self.yp_{}.setAttr('ub', 0)".format(rxn)
            if rxn in self.rHrev:
                exec 'self.d["ym_{0}.ub"] = self.ym_{0}.ub'.format(rxn)
                exec "self.ym_{}.setAttr('ub', 0)".format(rxn)
        elif rxn in self.rL:
            exec 'self.d["y_{0}.lb"] = self.y_{0}.lb'.format(rxn) #saving lb
            exec "self.y_{}.setAttr('lb', 1)".format(rxn) #changing lb
        else:
            pass# i.e. rxns without binary variables associated to them
        self.guro.update()
        self.guro.optimize()
        status = self.guro.getAttr('status')
        if status == 2:
            scoreZero = int(round(self.guro.objVal))
            solZero = [v.x for v in self.guro.getVars()]
        else:
            scoreZero = -1
            solZero = []
        # restoring the original bounds
        for key in self.d:
            try:
                l = key.split('.')
                exec "self.{}.setAttr('{}', {})".format(l[0], l[1], self.d[key])
            except AttributeError:
                pass
        self.guro.update()
        del self.d

        # 2. forcing forward direction
        exec 'self.fwd = self.{}.ub'.format(rxn) 
        if self.fwd == 0:#if irreversible in the rev dir
            scoreFwd = 0
            solFwd = []
            del self.fwd
        else:
            self.d = {}
            exec 'self.d["{0}.lb"] = self.{0}.lb'.format(rxn)
            exec "self.{}.setAttr('lb', self.eps)".format(rxn)
            if rxn in self.rH:
                exec 'self.d["yp_{0}.lb"] = self.yp_{0}.lb'.format(rxn)
                exec "self.yp_{}.setAttr('lb', 1)".format(rxn)
                if rxn in self.rHrev:
                    exec 'self.d["ym_{0}.ub"] = self.ym_{0}.ub'.format(rxn)
                    exec "self.ym_{}.setAttr('ub', 0)".format(rxn)
            elif rxn in self.rL:
                    exec 'self.d["y_{0}.ub"] = self.y_{0}.ub'.format(rxn)
                    exec "self.y_{}.setAttr('ub', 0)".format(rxn)
        self.guro.update()
        self.guro.optimize()
        status = self.guro.getAttr('status')
        if status == 2:
            scoreFwd = int(round(self.guro.objVal))
            solFwd = [v.x for v in self.guro.getVars()]
        else:
            scoreFwd = -1
            solFwd = []
        # restoring the original bounds
        for key in self.d:
            try:
                l = key.split('.')
                exec "self.{}.setAttr('{}', {})".format(l[0], l[1], self.d[key])
            except AttributeError:
                pass
        self.guro.update()
        del self.d

        # 3. forcing reverse direction
        # forcing reverse direction
        exec 'self.rev = self.{}.lb'.format(rxn)
        if self.rev == 0: #if irreversible
            scoreRev = 0
            solRev = []
            del self.rev
        else:
            self.d = {}
            exec 'self.d["{0}.ub"] = self.{0}.ub'.format(rxn)
            exec "self.{}.setAttr('ub', -self.eps)".format(rxn)
            if rxn in self.rH:
                exec 'self.d["yp_{0}.ub"] = self.yp_{0}.ub'.format(rxn)
                exec "self.yp_{}.setAttr('ub', 0)".format(rxn)
                if rxn in self.rHrev:
                    exec 'self.d["ym_{0}.lb"] = self.ym_{0}.lb'.format(rxn)
                    exec "self.ym_{}.setAttr('lb', 1)".format(rxn)
            elif rxn in self.rL:
                exec 'self.d["y_{0}.ub"] = self.y_{0}.ub'.format(rxn)
                exec "self.y_{}.setAttr('ub', 0)".format(rxn)
            else:
                pass
            self.guro.update()
            self.guro.optimize()
            status = self.guro.getAttr('status')
            if status == 2:
                scoreRev = int(round(self.guro.objVal))
                solRev = [v.x for v in self.guro.getVars()]
            else:
                scoreRev = -1
                solRev = []
            # restoring the original bounds
            for key in self.d:
                try:
                    l = key.split('.')
                    exec "self.{}.setAttr('{}', {})".format(l[0], l[1], self.d[key])
                except AttributeError:
                    pass
            self.guro.update()
            del self.d
        return scoreZero, scoreFwd, scoreRev, (solZero, solFwd, solRev)
    #}}}2

    def exploreAlternativeOptima(self, idRs):
        """
        """
        self.solve()
        wtSol = int(round(self.guro.objVal))
        scores = []
        imgeSols = {}
        scores.append(['wildType'] + [wtSol, '', ''])
        for rxn in idRs:
            try:
                s = self.modulateEnzymeActivityAndDirection(rxn)
                scores.append([rxn] + list(s)[:3])
                imgeSols[rxn] = s[3]
            except:
                scores.append([rxn] + [-1]*3)
                imgeSols[rxn] = ([], [], [])
        return scores, imgeSols

    #}}}

class MipSeparateFwdRev(object):
    """
    03 May 2012
    Creates a version of a cbModel in which reversible reactions into their
    fwd and rev components. Integer variables ensure that only one of the
    fwd, rev pair can carry a flux.
    """
    def __init__(self, cbm, mfr = [], eps = 1E-10):
        #{{{2
        self.m0 = cbm
        self.mfr = mfr
        self.eps = eps
        # identify mfr that are reversible
        self.mfrRev = [rxn for i, rxn in enumerate(self.m0.idRs) if 
                (self.m0.lb[i] < 0.) and (self.m0.ub[i] > 0.) and
                (rxn in self.mfr)]
        #NOTE assumes that irreversible rxns have lb = 0 (true in all
        #models I've come across)
        # Modify lb of irreversible mfr
        for rxn in set(self.mfr) - set(self.mfrRev):
            ind = self.m0.idRs.index(rxn)
            self.m0.lb[ind] = self.eps
        # 3. separate reversible rxns into fwd and rev reactions
        self.m = self.m0.separateFwdRevRxns()     
        #copying the S for the construction of the lhs
        self.lhsRows = [int(i) for i in self.m.S.row]   
        self.lhsCols = [int(i) for i in self.m.S.col]   
        self.lhsVals = [float(i) for i in self.m.S.data]
        self.rowNames = self.m.idSp[:]
        self.colNames = self.m.idRs[:]
        #}}}2
    
    def rowAndColNames(self):
        #{{{
        for rxn in self.mfrRev:
            self.rowNames.append('cfl_{}'.format(rxn)) #constraint fwd lb
            self.rowNames.append('cfu_{}'.format(rxn)) #constraint fwd ub
            self.rowNames.append('crl_{}'.format(rxn)) #constraint rev lb
            self.rowNames.append('cru_{}'.format(rxn)) #constraint rev ub
            self.colNames.append('y_{}'.format(rxn)) #integer variable    
        #}}}
    
    def buildLhs(self):
        #{{{
        for rxn in self.mfrRev:
            #fwd lowerbound constrain
            self.lhsRows.append(self.rowNames.index('cfl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn))
            self.lhsVals.append(1)
            self.lhsRows.append(self.rowNames.index('cfl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(-self.eps)
            #fwd upperbound constrain
            self.lhsRows.append(self.rowNames.index('cfu_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn))
            self.lhsVals.append(1)
            self.lhsRows.append(self.rowNames.index('cfu_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(-self.m.ub[self.colNames.index(rxn)])
            #rev lowerbound constrain
            self.lhsRows.append(self.rowNames.index('crl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn + '_rev'))
            self.lhsVals.append(1)
            self.lhsRows.append(self.rowNames.index('crl_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(self.eps)
            #rev upperbound constrain
            self.lhsRows.append(self.rowNames.index('cru_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index(rxn + '_rev'))
            self.lhsVals.append(1)
            self.lhsRows.append(self.rowNames.index('cru_{}'.format(rxn)))
            self.lhsCols.append(self.colNames.index('y_{}'.format(rxn)))
            self.lhsVals.append(self.m.ub[self.colNames.index(rxn + '_rev')])
        #}}}

    def buildRhsAndSenses(self):
        #{{{
        self.senses = 'E' * len(self.m.idSp)
        self.rhsRows = []
        self.rhsCols = []
        self.rhsVals = []
        for rxn in self.mfrRev:
            #fwd ub
            self.rhsRows.append(self.rowNames.index('cfl_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(0)
            self.senses += 'G'
            #fwd lb
            self.rhsRows.append(self.rowNames.index('cfu_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(0)
            self.senses += 'L'
            #
            #rev ub
            self.rhsRows.append(self.rowNames.index('crl_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(self.eps)
            self.senses += 'G'
            #rev lb
            self.rhsRows.append(self.rowNames.index('cru_{}'.format(rxn)))
            self.rhsCols.append(0)
            self.rhsVals.append(self.m.ub[self.colNames.index(rxn + '_rev')])
            self.senses += 'L'

        self.rhs = list(coo_matrix((self.rhsVals, (self.rhsRows, self.rhsCols)), 
                                   shape = (len(self.rowNames), 1)).toarray().flatten())    
        #}}}

class MipSeparateFwdRev_gurobi(MipSeparateFwdRev):
    def initMipGurobi(self):
        #{{{
        from numpy import array
        self.rowAndColNames()
        self.buildLhs()
        self.buildRhsAndSenses()
        self.guro = Model('minSumFluxes')
        #turning off the writing of the gurobi.log file
        self.guro.setParam('OutputFlag', 0) 
        # adding variables
        for i, rxn in enumerate(self.m.idRs):
            exec 'self.{0} = self.guro.addVar(lb = {}, ub = {}, vtype = GRB.CONTINUOUS, name = "{0}")'.format(rxn, self.m.lb[i], self.m.ub[i])  
        for var in [v for v in self.colNames if v not in self.m.idRs]:
            exec 'self.{0} = self.guro.addVar(lb = 0, ub = 1, vtype = GRB.BINARY, name = "{0}")'.format(var)
        self.guro.update()
        # adding constraints
        lhs = coo_matrix((self.lhsVals, (self.lhsRows, self.lhsCols))).toarray()
        for i, row in enumerate(lhs):
            nz = row.nonzero()[0]
            pair = zip(row[nz], array(self.colNames)[nz])
            s = ''
            for p in pair:
                s += '({} * self.{}) + '.format(p[0], p[1]) 
            s = s.rstrip(' + ')
            if self.senses[i] == 'E':
                s += ' == {}'.format(self.rhs[i])
            elif self.senses[i] == 'G':
                s += ' >= {}'.format(self.rhs[i])
            elif self.senses[i] == 'L':
                s += ' <= {}'.format(self.rhs[i])
            exec 'self.guro.addConstr( {}, "{}")'.format(s, self.rowNames[i])
            #}}}

    def minSumFluxes_gurobi(self):
        #{{{
        # setting the objective
        s = 'self.linobj = LinExpr([1.0] * len(self.m.idRs), ['
        for var in self.guro.getVars():
            if not var.varName.startswith('y_'): #excluding binary variables
                s += 'self.{}, '.format(var.varName)
        s = s.rstrip(', ')
        s += '])'
        exec s
        self.guro.setObjective(self.linobj, GRB.MINIMIZE)#1 for minimize
        self.guro.optimize()
        self.initialized = 1
        #}}}

################################################################################
# FUNCTIONS

######################################## 

def getZeroAndHighFrequencyRxns(scores, imgeSols, idRs, 
        activityThreshold = 1E-10):
    wtScore = scores[0][1]
    numSols = 0
    actDict = {}#counts for active reactions
    inactDict = {}#counts for inactive reactions
    for i, rxn in enumerate(idRs):
        row = scores[i+1][1:4]
        sol = imgeSols[rxn]
        for j, val in enumerate(row):
            if val >= wtScore:#select only max agreement scores
                numSols += 1
                for k, r in enumerate(idRs):
                    if abs(sol[j][k]) > activityThreshold:
                        try:
                            actDict[r] += 1
                        except KeyError:
                            actDict[r] = 1
                    else:
                        try:
                            inactDict[r] += 1
                        except KeyError:
                            inactDict[r] = 1
    
    zfr = set()#zero frequency reactions
    for rxn in inactDict:
        if inactDict[rxn] == numSols:# rxn inactive in all optimal solutions
            zfr.add(rxn)
    
    hfr = set()#high frequency reactions
    for rxn in actDict:
        if actDict[rxn] == numSols:# rxn active in all optimal solutiosn
            hfr.add(rxn)
    return zfr, hfr


def deleteCbmRxns(cbm, rl):
    """
    09 Dec 2011
    Deletes a list of reactions (or a single reaction) in rl from a
    CbModel instance.
    cbm [CbModel instance]
    rl [list or str] reaction(s) to be deleted
    RETURNS
    cbmNew [CbModel instance] model with the reaction(s) in rl deleted
    """
    if type(rl) == str:
        rl = [rl]
    indRs = [i for i, rxn in enumerate(cbm.idRs) if rxn not in rl]
    idRs =  [rxn for i, rxn in enumerate(cbm.idRs) if rxn not in rl]
    S = cbm.S.toarray()[:, indRs]
    actSpInd = [i for i, row in enumerate(S) if sum(abs(row)) > 1E-10]
    S = S[actSpInd]
    idSp = [met for i, met in enumerate(cbm.idSp) if i in actSpInd]
    lb = [val for i, val in enumerate(cbm.lb) if i in indRs]
    ub = [val for i, val in enumerate(cbm.ub) if i in indRs]
    cbmNew = CbModel(coo_matrix(S), idSp, idRs, lb, ub)
    #reactions dictionary
    if 'rxns' in dir(cbm):
        rxns = {}
        for key in cbm.rxns:
            if key in idRs:
                rxns[key] = cbm.rxns[key]
        cbmNew.rxns = rxns
    #species dictionary
    if 'species' in dir(cbm):
        species = {}
        for key in cbm.species:
            if key in idSp:
                species[key] = cbm.species[key]
        cbmNew.species = species
    #model genes
    if 'genes' in dir(cbm):
        import re
        reGenes = re.compile('[a-zA-Z0-9_\-\.]*')
        genes = set()
        for rxn in rxns:
            gl = reGenes.findall(rxns[rxn]['genes'])
            genes.update(gl)
        genes = genes - set(['', 'and', 'or'])
        cbmNew.genes = genes
    return cbmNew

######################################## 
# Minimizing the sum of fluxes

def getNzRxnsGurobi(mg, actThreshold, md):
    #computing the fluxes for reversible reactions
    revRxns = [i.varName for i in mg.guro.getVars() if i.varName.endswith('_rev')]
    revRxnsSols = []
    revRxnsIds = []
    for rxn in revRxns:
        fwd = rxn[:-4]
        exec 'valRev = mg.{}.x'.format(rxn)
        exec 'valFwd = mg.{}.x'.format(fwd)
        revRxnsSols.append(valFwd - valRev)
        revRxnsIds.append(fwd)
    sols = []
    for var in mg.guro.getVars():
        if not var.varName.startswith('y_'):# skip integer variables:
            if var.varName in revRxnsIds:
                sols.append([var.varName, revRxnsSols[revRxnsIds.index(var.varName)]])
            elif var.varName.endswith('_rev'):
                continue
            else:
                sols.append([var.varName, var.x])

    nzSols = [[i[0], i[1], abs(i[1]), md[i[0]]['reactants'],
               md[i[0]]['products'], md[i[0]]['lb'],
               md[i[0]]['ub']] for i in sols]

    return [i[0] for i in nzSols], nzSols

def altMinSumJs(m, mfr, eps, actThreshold):
    """
    Explores the space of flux distributions by minimizing the sum of fluxes
    """
    m1 = MipSeparateFwdRev_gurobi(m, mfr, eps)
    m1.initMipGurobi()
    m1.minSumFluxes_gurobi()
    sol0 = getNzRxnsGurobi(m1, actThreshold)
    maxUb = max(m.ub)
    minLb = min(m.lb)
    #check that all solutions are at least 10% below lb or ub
    for row in sol0[1]:
        if (row[1] > 0.9*maxUb) or (row[1] < 0.9*minLb):
            print 'warning: sol0 possible flux bound hit'
            print row

    cand0 = set(sol0[0]) - mfr
    essRxns = set()
    tried = set()
    sols = set()
    rep = []
    while cand0:
        cand1 = set()
        m1.d = {}
        for rxn in cand0:
            exec 'm1.d["{0}.lb"] = m1.{0}.lb'.format(rxn) #saving lb
            exec 'm1.d["{0}.ub"] = m1.{0}.ub'.format(rxn) #saving ub
            exec "m1.{}.setAttr('lb', 0.)".format(rxn) #changing lb
            exec "m1.{}.setAttr('ub', 0.)".format(rxn) #changing ub
            try:
                exec 'm1.d["{0}_rev.lb"] = m1.{0}_rev.lb'.format(rxn) #saving lb rev
                exec 'm1.d["{0}_rev.ub"] = m1.{0}_rev.ub'.format(rxn) #saving ub rev
                exec "m1.{}_rev.setAttr('lb', 0.)".format(rxn) #changing lb
                exec "m1.{}_rev.setAttr('ub', 0.)".format(rxn) #changing ub
            except AttributeError:
                pass
            m1.guro.update()
            m1.guro.optimize()
            status = m1.guro.getAttr('status')
            if status != 2:
                essRxns.add(rxn)
            else:
                tried.add(rxn)
                sol = getNzRxnsGurobi(m1)
                #check that all solutions are at least 10% below lb or ub
                for row in sol[1]:
                    if (row[1] > 0.9*maxUb) or (row[1] < 0.9*minLb):
                        print 'warning: sol possible flux bound hit'
                        print row
                sols.update(sol[0])
                cand1.update(sol[0])
            #restoring the original bounds
            for key in m1.d:
                try:
                    l = key.split('.')
                    exec "m1.{}.setAttr('{}', {})".format(l[0], l[1], m1.d[key])
                except AttributeError:
                    print key
            m1.guro.update()
        cand0 = cand1 - (tried | mfr)
        rep.append([len(cand0), len(tried | mfr), len(essRxns)])
        return sols




################################################################################
# TESTING
# Testing the MetabGeneExpModel_gurobi class and its methods
if __name__ == '__main__':
    # Example Figure 1 Rossell et al
    Srows = [0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8]
    Scols = [0, 1, 2, 1, 3, 8, 2, 3, 4,10, 3, 5, 4, 7, 8, 5, 6,11, 4, 6, 6, 7, 9,10]
    Svals = [1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1, 1,-1,-1,-1, 1, 1,-1, 1,-1]
    idRs = ['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'rBio']
    idSp = ['pep', 'ooa', 'acCoa', 'cit', 'mal', 'icit', 'gly', 'succ', 'acAld']
    rH = ['r3', 'r5', 'r6']
    rL = ['r2', 'r7']
    irrev = ['r0', 'r1', 'r2', 'r3', 'r4', 'r6', 'r9', 'r10', 'rBio']
    #creating lb and ub
    #ub = [10.]*len(idRs)
    ub = [cplex.infinity]*len(idRs)
    lb = []
    for i, rxn in enumerate(idRs):
        if rxn in irrev:
            lb.append(0.)
        else:
            #lb.append(-10.) 
            lb.append(-cplex.infinity)   
    S = coo_matrix((Svals, (Srows, Scols)))
    m = MetabGeneExpModel_gurobi(idSp, idRs, S, lb, ub, rH, rL)
    m.solve()
    for v in m.guro.getVars():
        print v.varName, v.x
    print 'Obj: ', m.guro.objVal
    
    solutionsTable = []
    for rxn in idRs:
        solutionsTable.append(m.exploreAltOptima(rxn))
    #printing out scores
    scores =[]
    for i, rxn in enumerate(idRs):
        scores.append([rxn, solutionsTable[i][:3]])
    print scores
