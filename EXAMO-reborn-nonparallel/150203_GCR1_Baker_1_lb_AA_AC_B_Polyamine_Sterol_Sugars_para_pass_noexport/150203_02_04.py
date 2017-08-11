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

################################################################################
# INPUTS
description = sys.argv[1]


# id of the biomass reaction
biomassRxn = 'R_biomass_published'

# Algorithm variables
#threshold above which a flux is considered to be larger than zero
activityThreshold = 1E-10

repetitions = int(sys.argv[2])

for repetition in range(repetitions):

	################################################################################
	# _02_minimizeNetwork_part_A.py

	################################################################################
	# INPUTS

	numProc = 1# 100	#number of parallel processes used

	md = importPickle('data/150126_degapped_iMM904_with_metabolite_mapping_complexes_36_stoichiometry_lb_AA_AC_B_Polyamine_Sterol.pkl')

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
	m0.lb[m0.idRs.index('R_biomass_published')] = 0.2879

	# Deleting zero frequency creactions
	fOutModeldict = 'examoModules/data/iMM904_examo_testactive_dict.pkl'	
	m = deleteCbmRxns(m0, fbr['zfr'])
	rev = [0 if val >= 0 else 1 for val in m0.lb]
	gene2rxn = {}
	rxns = {}
	for rxn in m.rxns:
   	    gene2rxn[rxn] = m.rxns[rxn]['genes']
    	    rxns[rxn] = m.rxns[rxn]
	mDict = {
	'S' : m.S,
	'idSp' : m.idSp,
	'idRs' : m.idRs,
	'lb' : m.lb,
	'ub' : m.ub,
	'rev' : rev,
	'genes' : m.genes,
	'gene2rxn' : gene2rxn,
	'rxns' : rxns,
	'descrip' : 'iMM904'}
	#exportPickle(mDict, fOutModeldict)

	# making sure that all hfr reactions are active to begin with
	act = findActiveRxns(m, 1E-10, mDict, fbr['hfr'])
	print "hello"
	cH = fbr['hfr'] & act

	# Making sure that biomass can be produced
	cH.add('R_biomass_published')

	#EG Create lists of extracellular reactions, extracellular transport reactions, and other compartmental transport reactions, so that the reactions can be pruned in that order first.
	EXrxns = ['R_EX_13BDglcn_e_', 'R_EX_2hb_e_', 'R_EX_2mbac_e_', 'R_EX_2mbald_e_', 'R_EX_2mbtoh_e_', 'R_EX_2mppal_e_', 'R_EX_2phetoh_e_', 'R_EX_3c3hmp_e_', 'R_EX_3mbald_e_', 'R_EX_3mop_e_', 'R_EX_4abut_e_', 'R_EX_4abz_e_', 'R_EX_5aop_e_', 'R_EX_Nbfortyr_e_', 'R_EX_abt_e_', 'R_EX_ac_e_', 'R_EX_acald_e_', 'R_EX_aces_e_', 'R_EX_ade_e_', 'R_EX_adn_e_', 'R_EX_akg_e_', 'R_EX_ala_L_e_', 'R_EX_alltn_e_', 'R_EX_alltt_e_', 'R_EX_amet_e_', 'R_EX_arab_L_e_', 'R_EX_arg_L_e_', 'R_EX_asn_L_e_', 'R_EX_asp_L_e_', 'R_EX_btd_RR_e_', 'R_EX_chol_e_', 'R_EX_cit_e_', 'R_EX_co2_e_', 'R_EX_csn_e_', 'R_EX_cys_L_e_', 'R_EX_cytd_e_', 'R_EX_dad_2_e_', 'R_EX_dca_e_', 'R_EX_dcyt_e_', 'R_EX_ddca_e_', 'R_EX_dgsn_e_', 'R_EX_din_e_', 'R_EX_dttp_e_', 'R_EX_duri_e_', 'R_EX_epist_e_', 'R_EX_epistest_SC_e_', 'R_EX_ergst_e_', 'R_EX_ergstest_SC_e_', 'R_EX_etha_e_', 'R_EX_etoh_e_', 'R_EX_fe2_e_', 'R_EX_fecost_e_', 'R_EX_fecostest_SC_e_', 'R_EX_fmn_e_', 'R_EX_for_e_', 'R_EX_fru_e_', 'R_EX_fum_e_', 'R_EX_g3pc_e_', 'R_EX_g3pi_e_', 'R_EX_gal_e_', 'R_EX_galur_e_', 'R_EX_gam6p_e_', 'R_EX_gcald_e_', 'R_EX_glc_e_', 'R_EX_gln_L_e_', 'R_EX_glu_L_e_', 'R_EX_glx_e_', 'R_EX_gly_e_', 'R_EX_glyc_e_', 'R_EX_gsn_e_', 'R_EX_gthox_e_', 'R_EX_gthrd_e_', 'R_EX_gua_e_', 'R_EX_h2o_e_', 'R_EX_h_e_', 'R_EX_hdca_e_', 'R_EX_hdcea_e_', 'R_EX_hexc_e_', 'R_EX_his_L_e_', 'R_EX_hxan_e_', 'R_EX_iamac_e_', 'R_EX_iamoh_e_', 'R_EX_ibutac_e_', 'R_EX_ibutoh_e_', 'R_EX_id3acald_e_', 'R_EX_ile_L_e_', 'R_EX_ind3eth_e_', 'R_EX_inost_e_', 'R_EX_ins_e_', 'R_EX_lac_D_e_', 'R_EX_lac_L_e_', 'R_EX_lanost_e_', 'R_EX_lanostest_SC_e_', 'R_EX_leu_L_e_', 'R_EX_lys_L_e_', 'R_EX_mal_L_e_', 'R_EX_malt_e_', 'R_EX_man_e_', 'R_EX_melib_e_', 'R_EX_met_L_e_', 'R_EX_nac_e_', 'R_EX_nadp_e_', 'R_EX_nh4_e_', 'R_EX_nmn_e_', 'R_EX_o2_e_', 'R_EX_oaa_e_', 'R_EX_ocdca_e_', 'R_EX_ocdcea_e_', 'R_EX_ocdcya_e_', 'R_EX_orn_e_', 'R_EX_pacald_e_', 'R_EX_pap_e_', 'R_EX_pc_SC_e_', 'R_EX_pectin_e_', 'R_EX_phe_L_e_', 'R_EX_pheac_e_', 'R_EX_pi_e_', 'R_EX_pnto_R_e_', 'R_EX_pro_L_e_', 'R_EX_ptd1ino_SC_e_', 'R_EX_ptrc_e_', 'R_EX_pyr_e_', 'R_EX_rib_D_e_', 'R_EX_ribflv_e_', 'R_EX_sbt_D_e_', 'R_EX_sbt_L_e_', 'R_EX_ser_L_e_', 'R_EX_so3_e_', 'R_EX_so4_e_', 'R_EX_spmd_e_', 'R_EX_sprm_e_', 'R_EX_srb_L_e_', 'R_EX_succ_e_', 'R_EX_sucr_e_', 'R_EX_thm_e_', 'R_EX_thmmp_e_', 'R_EX_thmpp_e_', 'R_EX_thr_L_e_', 'R_EX_thym_e_', 'R_EX_thymd_e_', 'R_EX_tre_e_', 'R_EX_trp_L_e_', 'R_EX_ttdca_e_', 'R_EX_tyr_L_e_', 'R_EX_ura_e_', 'R_EX_urea_e_', 'R_EX_uri_e_', 'R_EX_val_L_e_', 'R_EX_xan_e_', 'R_EX_xtsn_e_', 'R_EX_xyl_D_e_', 'R_EX_xylt_e_', 'R_EX_zymst_e_', 'R_EX_zymstest_SC_e_']

	EXtrrxns = ['R_2MBACt', 'R_2MBALDt', 'R_2MBTOHt', 'R_2MPPALt', 'R_2PHETOHt', 'R_3C3HMPt', 'R_3MBALDt', 'R_3MOPt', 'R_4ABZt', 'R_5AOPt2', 'R_ABTt', 'R_ABUTt2r', 'R_ACALDt', 'R_ACESt', 'R_ACt2r', 'R_ACtr', 'R_ADEt2', 'R_ADNt2', 'R_AKGMAL', 'R_AKGt2r', 'R_ALAt2r', 'R_ALLTNti', 'R_ALLTTti', 'R_AMETt2', 'R_ARAB_Lt', 'R_ARGt2r', 'R_ASNt2r', 'R_ASPt2r', 'R_ATPS', 'R_CHLt2', 'R_CITt2r', 'R_CO2t', 'R_CSNt2', 'R_CYSt2r', 'R_CYTDt2', 'R_DADNt2', 'R_DCYTt2', 'R_DGSNt2', 'R_DINSt2', 'R_DTTPt', 'R_DURIt2', 'R_D_LACt2', 'R_EPISTt', 'R_ERGSTt', 'R_ETHAt', 'R_ETOHt', 'R_FE2t', 'R_FECOSTt', 'R_FORt', 'R_FRUt2', 'R_FUMt2r', 'R_G3PCt', 'R_GALt2', 'R_GAM6Pt', 'R_GCALDt', 'R_GLCt1', 'R_GLNt2r', 'R_GLUt2r', 'R_GLXt', 'R_GLYCt', 'R_GLYCt2', 'R_GLYt2r', 'R_GSNt2', 'R_GTHOXti', 'R_GTHRDt2', 'R_GUAt2r', 'R_H2Ot', 'R_HDCAt', 'R_HDCEAt', 'R_HISt2r', 'R_HXANt2r', 'R_IAMACt', 'R_IAMOHt', 'R_IBUTACt', 'R_IBUTOHt', 'R_ID3ACALDt', 'R_ILEt2r', 'R_IND3ETHt', 'R_INSTt2', 'R_INSt2', 'R_LANOSTt', 'R_LEUt2r', 'R_LYSt2r', 'R_L_LACt2r', 'R_MALTt2', 'R_MALt2r', 'R_MANt2', 'R_MELIBt2', 'R_METt2r', 'R_NACt', 'R_NADPt', 'R_NFORTYRt', 'R_NH4t', 'R_NH4ti', 'R_NMNTP', 'R_O2t', 'R_OAAt', 'R_OCDCAt', 'R_OCDCEAt', 'R_OCDCYAt', 'R_ORNt2r', 'R_PACALDt', 'R_PAPt', 'R_PHEACt', 'R_PHEt2r', 'R_PIt2r', 'R_PNTOt2', 'R_PROt2r', 'R_PTRCt3i', 'R_PYRt', 'R_PYRt2', 'R_RIBFLVt2', 'R_RIBt2', 'R_SBT_Dt', 'R_SBT_Lt', 'R_SERt2r', 'R_SO3ti', 'R_SO4ti', 'R_SPMDt3i', 'R_SPRMt2i', 'R_SRB_Lt', 'R_SUCCt2r', 'R_THMDt2', 'R_THMt2', 'R_THRt2r', 'R_THYMt3r', 'R_TREt2', 'R_TRPt2r', 'R_TTDCAtr', 'R_TYRt2r', 'R_URAt2', 'R_UREA2t2', 'R_URIt2', 'R_VALt2r', 'R_XANt', 'R_XTSNt2', 'R_XYLTt', 'R_XYLt', 'R_ZYMSTt'] 

	Othertrrxns = ['R_ASNt6', 'R_ASNt7', 'R_CO2tv', 'R_GLCNtv', 'R_GLCtv', 'R_GLNt6', 'R_GLNt7', 'R_H2Otv', 'R_ILEt6', 'R_ILEt7', 'R_LEUt6', 'R_LEUt7', 'R_PEtv_SC', 'R_PStv_SC', 'R_TREt2v', 'R_TYRt6', 'R_TYRt7', 'R_ACRNtp', 'R_AKGtp', 'R_ASPGLUtp', 'R_ATP2tp_H', 'R_ATPtp_H', 'R_CITtap', 'R_CRNCARtp', 'R_CRNtp', 'R_CYSTtp', 'R_GLXtp', 'R_H2Otp', 'R_HCYSt2p', 'R_MALOAAtp', 'R_NH4tp', 'R_PIt2p', 'R_PYRt2p', 'R_TRDOXtp', 'R_TRDRDtp', 'R_ASPt2n', 'R_ASPt5n', 'R_CO2tn', 'R_H2Otn', 'R_HCO3tn', 'R_2DDA7Ptm', 'R_2DHPtm', 'R_2MBALDtm', 'R_2MBTOHtm', 'R_2MPPALtm', 'R_2OBUTtm', 'R_2OXOADPtim', 'R_2PHETOHtm', 'R_34HPPt2m', 'R_3C3HMPtm', 'R_3DH5HPBtm', 'R_3MBALDtm', 'R_3MOBtm', 'R_3MOPtm', 'R_3OPHB_5tm', 'R_4ABZtm', 'R_4HBZtm', 'R_5AOPtm', 'R_ACALDtm', 'R_ACRNtm', 'R_ACtm', 'R_AHCYStm', 'R_ALAtmi', 'R_AMETtm', 'R_ASPGLU2m', 'R_ASPt2m', 'R_ATPtm_H', 'R_CITtam', 'R_CITtbm', 'R_CITtcm', 'R_CO2tm', 'R_COAtim', 'R_CRNCARtm', 'R_CRNtim', 'R_CTPtm', 'R_DHAPtm', 'R_DHNPTtm', 'R_DHPTtm', 'R_D_LACt2m', 'R_D_LACtm', 'R_E4Ptm', 'R_ETOHtm', 'R_FE2utm', 'R_FORtm', 'R_FRDcm', 'R_GCALDtm', 'R_GLUt5m', 'R_GLUt7m', 'R_GLYC3Ptm', 'R_GLYt2m', 'R_GTPt2m', 'R_H2Otm', 'R_HMGCOAtm', 'R_IAMOHtm', 'R_IBUTOHtm', 'R_ID3ACALDtm', 'R_ILEtmi', 'R_IND3ETHtm', 'R_IPDPtm', 'R_MALtm', 'R_NH4tm', 'R_O2tm', 'R_OAAt2m', 'R_ORNt3m', 'R_PACALDtm', 'R_PAN4Ptm', 'R_PANTtm', 'R_PAPtm', 'R_PAtm_SC', 'R_PENDPtm', 'R_PEtm_SC', 'R_PIt2m', 'R_PIt5m', 'R_PPPG9tm', 'R_PROtm', 'R_PStm_SC', 'R_PYRt2m', 'R_SERt2m', 'R_SUCCtm', 'R_SUCFUMtm', 'R_THRt2m', 'R_TYRt2m', 'R_VALt2m',  'R_6PGLter', 'R_DOLPt2er', 'R_ERGSTter', 'R_ERGTETROLter', 'R_G6Pter', 'R_H2Oter', 'R_MANNANter', 'R_O2ter', 'R_SQ23EPXter', 'R_SQLter']

	#EG Make a directory for temporary files for every time a rxn is pruned
	mbaCandRxnsDirectorySubset = 'examoModules/data/%s_%s/' % (description, str(repetition))
	if not os.path.exists(mbaCandRxnsDirectorySubset):
		os.mkdir(mbaCandRxnsDirectorySubset, 0777)

	#Run the MBA
	for i in range(numProc):
		locTime = time.localtime()
		pid = os.getpid()
		timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
		tag = '%s_%s_%s' % (description, pid, timeStr)
		try:
		#EG Added despricription, repetition, and lists of compartmental reactions to the function
			cr = iterativePrunning(i, m, cH, description, repetition, activityThreshold, EXrxns, EXtrrxns, Othertrrxns)
			exportPickle(cr, fOutMbaCandRxns % tag)
		except:
			print 'gurobi error, no solution found %s'  % description

	print "after"

	#EG Count the number of exported pickle files, and do not proceed until all files have been written
	file_count = 0
	while file_count != numProc:
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
	md = importPickle('data/150126_degapped_iMM904_with_metabolite_mapping_complexes_36_stoichiometry_lb_AA_AC_B_Polyamine_Sterol.pkl')

	# Importing model information
	fbr = importPickle('data/freqBasedRxns_%s.pkl' % description)

	fOutModel = 'data/iMM904_examo_%s_%s.pkl'

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

	#EG Rather than identifying which reactions need to be added to make all of the hfrs active, reactions will be added until the a flux can be achieved for the biomass reaction 
	for num in orderedFreq:
		cH.update(freq[num])
		excRxns = set(m.idRs) - cH
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
			    'descrip' : 'iMM904 examo %s' % description}
			exportPickle(mDict, fOutModel % (description, str(repetition) + '_dict'))#EG End of the original 3rd script

			#EG Beginning of the 4th script
			################################################################################
			# INPUTS
			eps = 1E-10

			fModelExamo = 'data/iMM904_examo_%s_%s_dict.pkl'
			fFreqBasedRxns = 'data/freqBasedRxns_%s.pkl'

			fOutMetabState = 'data/metabolicState_%s_%s.csv'

			################################################################################
			# STATEMENTS

			hfr = importPickle(fFreqBasedRxns % description)['hfr']
			md = importPickle(fModelExamo % (description, str(repetition)))
			mtry = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'], md['genes'])
			hfr = hfr & set(mtry.idRs)
			#forcing biomass production
			mtry.lb[mtry.idRs.index('R_biomass_published')] = 0.2879

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
