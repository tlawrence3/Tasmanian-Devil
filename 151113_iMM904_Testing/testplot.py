#%matplotlib inline
import matplotlib.pyplot as plt
import scipy.io
import os
import re
import csv
import numpy as np
import math
from examoModules import *

#store file names and conditions for EXAMO and COBRA results
pickle_files = []
mat_files = []
condition_list = []
input_Models = {}
COBRA_Results = {}
EXAMO_Results = {}
EXAMO_Results_Ex = {}
for file in os.listdir(os.getcwd()+'/data/COBRAResults/'):
    file_condition_entry = re.findall('(?<=Results_\d{6}_)(.*?)(?=_0.25_iMM904_NAD)',file,re.DOTALL)
    file_input_entry = re.findall('(?<=FTHFLm_)(.*?)(?=.mat)',file,re.DOTALL)
    if len(file_input_entry) > 0:
        file_input_entry_copy = 'iMM904_NADcorrected_1127_FTHFLm_'+file_input_entry[0]
    else:
        file_input_entry_copy = 'iMM904_NADcorrected_1127_FTHFLm'
    file_input_entry_copy = file_input_entry_copy + '.pkl'
    if file_input_entry_copy not in pickle_files:
        pickle_files.append(file_input_entry_copy)
    mat_files.append(file)
    try:
        file_name_combine = file_condition_entry[0].title() + '_' + file_input_entry[0]
    except:
        file_name_combine = file_condition_entry[0].title()
    file_name_combine = re.sub('YPD','lb',file_name_combine)
    file_name_combine = re.sub('YPEtoh','lb',file_name_combine)
    file_name_combine = re.sub('Rintala_Anaerobic','lb',file_name_combine)
    file_name_combine = re.sub('Rintala_Aerobic','lb',file_name_combine)
    condition_list.append(file_name_combine)
    input_Models[file_name_combine] = file_input_entry_copy 
    COBRA_Results[file_name_combine] = {'file_name': file}
    file_name_combine_for_metabolic_state = re.sub('Results', 'metabolicState', file)
    EXAMO_Results[file_name_combine] = {'file_name': file_name_combine_for_metabolic_state[0:-4]}
    EXAMO_Results_Ex[file_name_combine] = {'file_name': file_name_combine_for_metabolic_state[0:-4] + '_EXrxns_EXtrrxns_Othertrrxns'}
    
#group conditions by medium
glucose_list = []
ethanol_list = []
aerobic_list = []
anaerobic_list = []
negative_control_list = []
for name in condition_list:
    if 'Glucose' in name:
        glucose_list.append(name)
    if 'Ethanol' in name:
        ethanol_list.append(name)
    if 'Aerobic' in name:
        aerobic_list.append(name)
    if 'Anaerobic' in name:
        anaerobic_list.append(name)
    if 'Negative_Control' in name:
        negative_control_list.append(name)
        
#import pickle models using cobrapy
pickle_models = {}
md_models = {}
for model in input_Models:
    pickle_models[model] = importPickle(os.getcwd()+'/data/models/'+input_Models[model])
    md_models[model] = CbModel(pickle_models[model]['S'], pickle_models[model]['idSp'],
                               pickle_models[model]['idRs'], pickle_models[model]['lb'],
                               pickle_models[model]['ub'], pickle_models[model]['rxns'],
                               pickle_models[model]['genes'])
        
#import list of tested genes for YPD
YPD_tested_genes_list = []
YPD_tested_genes_file = open(os.getcwd()+'/data/experimentalResults/'+'Giaver_Tested_Genes.txt')
YPD = YPD_tested_genes_file.readlines()
for i, item1 in enumerate(YPD):
    YPD_tested_genes_list.append(item1.strip())

#import list of essential genes for YPD
YPD_essential_genes_list = []
YPD_genes_file = open(os.getcwd()+'/data/experimentalResults/'+'Giaver_YPD_Essential_Genes.txt')
YPD = YPD_genes_file.readlines()
for i, item1 in enumerate(YPD):
    YPD_essential_genes_list.append(item1.strip())

#import list of tested genes for YPEtoh
YPEtoh_tested_genes_list = []
YPEtoh_tested_genes_file = open(os.getcwd()+'/data/experimentalResults/'+'Snitkin_Tested_Mutants.txt')
YPEtoh_tested_genes = YPEtoh_tested_genes_file.readlines()
for i, item1 in enumerate(YPEtoh_tested_genes):
    YPEtoh_tested_genes_list.append(item1.strip())

#import list of essential genes for YPEtoh
YPEtoh_essential_genes_list = []
YPEtoh_genes_file = open(os.getcwd()+'/data/experimentalResults/'+'Snitkin_YPEtoh_Essential_Genes.txt')
YPEtoh = YPEtoh_genes_file.readlines()
for i, item1 in enumerate(YPEtoh):
    YPEtoh_essential_genes_list.append(item1.strip())
    
#import list of fluxes
Rintala_Fluxes_file = open(os.getcwd()+'/data/experimentalResults/'+'Rintala_Fluxes.csv', 'r')
csvreader_experimental_fluxes = csv.reader(Rintala_Fluxes_file)
counter = 0
experimental_fluxes_rxns = {}
for row in csvreader_experimental_fluxes:
    counter += 1
    if counter > 1:
        rxnid = row[0]
        rxn1 = row[1]
        rxn2 = row[2]
        rxn3 = row[3]
        if rxn3 != '':
            rxn_list = [rxn1,rxn2,rxn3]
        elif rxn2 != '':
            rxn_list = [rxn1, rxn2]
        else:
            rxn_list = [rxn1]
        Rep1lb_Aerobic = float(row[4])
        Rep1ub_Aerobic = float(row[5])
        Rep2lb_Aerobic = float(row[6])
        Rep2ub_Aerobic = float(row[7])
        Rep1lb_Anaerobic = float(row[8])
        Rep1ub_Anaerobic = float(row[9])
        Rep2lb_Anaerobic = float(row[10])
        Rep2ub_Anaerobic = float(row[11])
                        
        Aerobic_flux_list = [Rep1lb_Aerobic,Rep1ub_Aerobic,Rep2lb_Aerobic,Rep2ub_Aerobic]
        Anaerobic_flux_list = [Rep1lb_Anaerobic,Rep1ub_Anaerobic,Rep2lb_Anaerobic,Rep2ub_Anaerobic]
        Aerobic_average = round(np.average(Aerobic_flux_list),4)
        Aerobic_std = round(np.std(Aerobic_flux_list),4)
        Anaerobic_average = round(np.average(Anaerobic_flux_list),4)
        Anaerobic_std = round(np.std(Anaerobic_flux_list),4)
        experimental_fluxes_rxns[rxnid] = {'rxns': rxn_list, 'Aerobic_average': Aerobic_average, 'Anaerobic_average': Anaerobic_average, 'Aerobic_std': Aerobic_std, 'Anaerobic_std': Anaerobic_std}

#determine the possible genes for each model
genes_possible_model_dict = {}
for state in md_models:
    gene_list = []
    for rxn in md_models[model].gene2rxn.keys():
        g2r = md_models[model].gene2rxn[rxn]
        g2r = str(g2r)
        i = g2r.split('or')
        for count, j in enumerate(i):
            j = j.split('and')
            for k in j:
                k = k.translate(None, ' ()')
                if not k:
                    continue
                else:
                    if k not in gene_list:
                        gene_list.append(k)
    genes_possible_model_dict[state] = gene_list

#Import fluxes for EXAMO
results = [EXAMO_Results, EXAMO_Results_Ex]
for result in results:
    fluxavgdict = {}
    fluxvardict = {}
    fluxstddict = {}
    rxninallreps = {}
    rxninnonereps = {}
    idRsfluxstate = {}
    allidRsfluxstate = {}
    for state in result:
        metabolicState_file_1 = open(os.getcwd()+'/data/%s_0.csv' % result[state]['file_name'], 'r')
        metabolicState_file_2 = open(os.getcwd()+'/data/%s_1.csv' % result[state]['file_name'], 'r')
        metabolicState_file_3 = open(os.getcwd()+'/data/%s_2.csv' % result[state]['file_name'], 'r')
        metabolicState_file_4 = open(os.getcwd()+'/data/%s_3.csv' % result[state]['file_name'], 'r')
        metabolicState_file_5 = open(os.getcwd()+'/data/%s_4.csv' % result[state]['file_name'], 'r')
        csvreader1 = csv.reader(metabolicState_file_1)
        csvreader2 = csv.reader(metabolicState_file_2)
        csvreader3 = csv.reader(metabolicState_file_3)
        csvreader4 = csv.reader(metabolicState_file_4)
        csvreader5 = csv.reader(metabolicState_file_5)

        idRsflux1 = {}
        idRscsv1 = []
        flux1 = []

        for row in csvreader1:
            idRscsv1.append(row[0])
            flux1.append(float(row[1]))
        for i, item1 in enumerate(idRscsv1):
            idRsflux1[item1] = flux1[i]

        idRsflux2 = {}
        idRscsv2 = []
        flux2 = []

        for row in csvreader2:
            idRscsv2.append(row[0])
            flux2.append(float(row[1]))
        for i, item1 in enumerate(idRscsv2):
            idRsflux2[item1] = flux2[i]

        idRsflux3 = {}
        idRscsv3 = []
        flux3 = []

        for row in csvreader3:
            idRscsv3.append(row[0])
            flux3.append(float(row[1]))
        for i, item1 in enumerate(idRscsv3):
            idRsflux3[item1] = flux3[i]

        idRsflux4 = {}
        idRscsv4 = []
        flux4 = []

        for row in csvreader4:
            idRscsv4.append(row[0])
            flux4.append(float(row[1]))
        for i, item1 in enumerate(idRscsv4):
            idRsflux4[item1] = flux4[i]

        idRsflux5 = {}
        idRscsv5 = []
        flux5 = []

        for row in csvreader5:
            idRscsv5.append(row[0])
            flux5.append(float(row[1]))
        for i, item1 in enumerate(idRscsv5):
            idRsflux5[item1] = flux5[i]

        allidRsflux1 = {}
        allidRsflux2 = {}
        allidRsflux3 = {}
        allidRsflux4 = {}
        allidRsflux5 = {}
        fluxavgdict1 = {}
        fluxvardict1 = {}
        fluxstddict1 = {}
        rxninallreps_list = []
        rxninnonereps_list = []
        for t in pickle_models[state]['rxns']:
            notincond1 = 0
            if t not in idRsflux1:
                allidRsflux1[t] = 0
                notincond1 += 1
            else:
                allidRsflux1[t] = idRsflux1[t]
            if t not in idRsflux2:
                allidRsflux2[t] = 0
                notincond1 += 1
            else:
                allidRsflux2[t] = idRsflux2[t]
            if t not in idRsflux3:
                allidRsflux3[t] = 0
                notincond1 += 1
            else:
                allidRsflux3[t] = idRsflux3[t]
            if t not in idRsflux4:
                allidRsflux4[t] = 0
                notincond1 += 1
            else:
                allidRsflux4[t] = idRsflux4[t]
            if t not in idRsflux5:
                allidRsflux5[t] = 0
                notincond1 += 1
            else:
                allidRsflux5[t] = idRsflux5[t]
            if notincond1 == 0:
                rxninallreps_list.append(t)
            if notincond1 == 5:
                fluxavgdict1[t] = 0
                fluxvardict1[t] = 0
                fluxstddict1[t] = 0
                rxninnonereps_list.append(t) 
            else:
                fluxavgdict[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
                fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]] 
                fluxvardict[t] = round(np.var(fluxes),4)
                fluxstddict[t] = round(np.std(fluxes),4)

        result[state].update({'fluxavgdict': fluxavgdict, 'fluxvardict': fluxvardict, 'fluxstddict': fluxstddict, 'rxninallreps': rxninallreps, 'rxninnonereps': rxninnonereps, 'idRsfluxstate': {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}, 'allidRsfluxstate': {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}})        

#determine the number of negative and positive genes present for the glucose, negative control, and ethanol data sets for EXAMO
gene_inclusion_conditions = [glucose_list, negative_control_list]
for result in results:
    for state in result:
        for condition in gene_inclusion_conditions:
            if state in condition:
                average_genes = []
                idRsflux_dict = {}
                idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
                for flux_state in idRsflux_list:
                    gene_list_true_positive = []
                    gene_list_false_positive = []
                    for rxn in result[state]['idRsfluxstate'][flux_state].keys():
                        g2r = md_models[state].gene2rxn[rxn]
                        g2r = str(g2r)
                        i = g2r.split('or')
                        for count, j in enumerate(i):
                            j = j.split('and')
                            for k in j:
                                k = k.translate(None, ' ()')
                                if not k:
                                    continue
                                else:
                                    if k in YPD_tested_genes_list:
                                        if k in YPD_essential_genes_list:
                                            if k not in gene_list_true_positive:
                                                gene_list_true_positive.append(k)
                                        else:
                                            if k not in gene_list_false_positive:
                                                gene_list_false_positive.append(k)
                    gene_list_true_negative = []
                    gene_list_false_negative = []
                    for rxn in md_models[state].gene2rxn.keys():
                        if rxn not in result[state]['idRsfluxstate'][flux_state].keys():
                            g2r = md_models[state].gene2rxn[rxn]
                            g2r = str(g2r)
                            i = g2r.split('or')
                            for count, j in enumerate(i):
                                j = j.split('and')
                                for k in j:
                                    k = k.translate(None, ' ()')
                                    if not k:
                                        continue
                                    else:
                                        if k in YPD_tested_genes_list:
                                            if k not in YPD_essential_genes_list:
                                                if k not in gene_list_true_negative:
                                                    gene_list_true_negative.append(k)
                                            else:
                                                if k not in gene_list_false_negative:
                                                    gene_list_false_negative.append(k)
                    idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
                sensitivity_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))]
                ppv_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))]
                sensitivity_average = round(np.average(sensitivity_list),4)
                sensitivity_std = round(np.std(sensitivity_list),4)
                ppv_average = round(np.average(ppv_list),4)
                ppv_std = round(np.std(ppv_list),4)
                result[state].update({'glucose_sensitivity_average': sensitivity_average, 'glucose_sensitivity_std': sensitivity_std, 'glucose_PPV_average': ppv_average, 'glucose_PPV_std': ppv_std})

gene_inclusion_conditions = [ethanol_list, negative_control_list]                
for result in results:
    for state in result:
        for condition in gene_inclusion_conditions:
            if state in condition:
                average_genes = []
                idRsflux_dict = {}
                idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
                for flux_state in idRsflux_list:
                    gene_list_true_positive = []
                    gene_list_false_positive = []
                    for rxn in result[state]['idRsfluxstate'][flux_state].keys():
                        g2r = md_models[state].gene2rxn[rxn]
                        g2r = str(g2r)
                        i = g2r.split('or')
                        for count, j in enumerate(i):
                            j = j.split('and')
                            for k in j:
                                k = k.translate(None, ' ()')
                                if not k:
                                    continue
                                else:
                                    if k in YPD_tested_genes_list:
                                        if k in YPD_essential_genes_list:
                                            if k not in gene_list_true_positive:
                                                gene_list_true_positive.append(k)
                                        else:
                                            if k not in gene_list_false_positive:
                                                gene_list_false_positive.append(k)
                    gene_list_true_negative = []
                    gene_list_false_negative = []
                    for rxn in md_models[state].gene2rxn.keys():
                        if rxn not in result[state]['idRsfluxstate'][flux_state].keys():
                            g2r = md_models[state].gene2rxn[rxn]
                            g2r = str(g2r)
                            i = g2r.split('or')
                            for count, j in enumerate(i):
                                j = j.split('and')
                                for k in j:
                                    k = k.translate(None, ' ()')
                                    if not k:
                                        continue
                                    else:
                                        if k in YPD_tested_genes_list:
                                            if k not in YPD_essential_genes_list:
                                                if k not in gene_list_true_negative:
                                                    gene_list_true_negative.append(k)
                                            else:
                                                if k not in gene_list_false_negative:
                                                    gene_list_false_negative.append(k)
                    idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
                sensitivity_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))]
                ppv_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))]
                sensitivity_average = round(np.average(sensitivity_list),4)
                sensitivity_std = round(np.std(sensitivity_list),4)
                ppv_average = round(np.average(ppv_list),4)
                ppv_std = round(np.std(ppv_list),4)
                result[state].update({'ethanol_sensitivity_average': sensitivity_average, 'ethanol_sensitivity_std': sensitivity_std, 'ethanol_PPV_average': ppv_average, 'ethanol_PPV_std': ppv_std})                

#Compare the measured flux data with the modeled flux data for EXAMO for aerobic, anaerobic, and negative control
flux_conditions = [aerobic_list, negative_control_list]
for result in results:
    for state in result:
        for condition in flux_conditions:
            if state in condition:
                rxn_dict = {}
                tot_abs_difference = 0
                experimental_total = 0
                for rxnid in experimental_fluxes_rxns.keys():
                    if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
                        rxn_sum = 0
                        rxn_var_sum = 0
                        for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
                            rxn_sum = rxn_sum + result[state]['fluxavgdict'][rxn]
                            rxn_var_sum = rxn_var_sum + (result[state]['fluxvardict'][rxn]**2)
                        rxn_var = math.sqrt(rxn_var_sum)
                    else:
                        for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
                            rxn_sum = result[state]['fluxavgdict'][rxn]
                            rxn_var = result[state]['fluxvardict'][rxn]
                    #experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_average']
                    abs_difference = abs(abs(rxn_sum) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_average']))
                    tot_abs_difference = tot_abs_difference + abs_difference
                    rxn_dict[rxnid] = {'Experimental_flux': experimental_fluxes_rxns[rxnid]['Aerobic_average'], 'Modeled_Flux': rxn_sum, 'abs_difference': abs_difference}
                #total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)             
                result[state].update({'Aerobic_tot_abs_difference': tot_abs_difference, 'Aerobic_rxn_dict': rxn_dict})

flux_conditions = [anaerobic_list, negative_control_list]
for result in results:
    for state in result:
        for condition in flux_conditions:
            if state in condition:
                rxn_dict = {}
                tot_abs_difference = 0
                experimental_total = 0
                for rxnid in experimental_fluxes_rxns.keys():
                    if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
                        rxn_sum = 0
                        rxn_var_sum = 0
                        for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
                            rxn_sum = rxn_sum + result[state]['fluxavgdict'][rxn]
                            rxn_var_sum = rxn_var_sum + (result[state]['fluxvardict'][rxn]**2)
                        rxn_var = math.sqrt(rxn_var_sum)
                    else:
                        for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
                            rxn_sum = result[state]['fluxavgdict'][rxn]
                            rxn_var = result[state]['fluxvardict'][rxn]
                    #experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_average']
                    abs_difference = abs(abs(rxn_sum) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_average']))
                    tot_abs_difference = tot_abs_difference + abs_difference
                    rxn_dict[rxnid] = {'Experimental_flux': experimental_fluxes_rxns[rxnid]['Anaerobic_average'], 'Modeled_Flux': rxn_sum, 'abs_difference': abs_difference}
                #total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)             
                result[state].update({'Anaerobic_tot_abs_difference': tot_abs_difference, 'Anaerobic_rxn_dict': rxn_dict})
                
print 'hi'        