#Need to include 
%matplotlib notebook

from __future__ import division
import os
import csv
import re
import math
import numpy as np
import scipy.spatial.distance
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import glob

from sys import path; path.append('./151113_iMM904_Testing/')
from examoModules import *


#First import all of the EXAMO-ARC V:1 data
pickle_models_list = ['iMM904_NADcorrected_1127_FTHFLm', 'iMM904_NADcorrected_1127_FTHFLm_g', 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', 'iMM904_NADcorrected_1127_FTHFLm_YPD', 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']

pickle_models_list_dataset_S2 = ['iMM904_NADcorrected_1127_FTHFLm']

pickle_models_list_dataset_S2_original_model = ['iMM904_blkRxnsDeleted_dict']

freqBased_list = ['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']

freqBased_list_dataset_S2 = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']

freqBased_list_dataset_S2_original_model = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']

freqBased_dict = {'151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm':'iMM904_NADcorrected_1127_FTHFLm' , '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'}

freqBased_dict_dataset_S2 = {'151006_Aerobic_0.25_0': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_1': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_2': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_3': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_4': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_0': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_1': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_2': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_3': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_4': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_0': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_1': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_2': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_3': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_4': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_0': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_1': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_2': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_3': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_4': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_0': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_1': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_2': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_3': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_4': 'iMM904_NADcorrected_1127_FTHFLm'}

freqBased_dict_dataset_S2_original_model = {'151006_Aerobic_0.25_0': 'iMM904_blkRxnsDeleted_dict', '151006_Aerobic_0.25_1': 'iMM904_blkRxnsDeleted_dict', '151006_Aerobic_0.25_2': 'iMM904_blkRxnsDeleted_dict', '151006_Aerobic_0.25_3': 'iMM904_blkRxnsDeleted_dict', '151006_Aerobic_0.25_4': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25_0': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25_1': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25_2': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25_3': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25_4': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25_0': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25_1': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25_2': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25_3': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25_4': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25_0': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25_1': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25_2': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25_3': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25_4': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25_0': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25_1': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25_2': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25_3': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25_4': 'iMM904_blkRxnsDeleted_dict'}

metabolicState_list = ['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', 
'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']

metabolicState_list_dataset_S2 = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_list_dataset_S2_eps = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_list_dataset_S2_original_model = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_list_dataset_S2_original_model_eps = ['151006_Aerobic_0.25', '151006_Anaerobic_0.25', '151012_ethanol_0.25', '151012_glucose_0.25', '151015_negative_control_0.25']


#Count the number of mbaCandRxns files
mbaCandRxns_dict = {}
for state in metabolicState_list:
	count0 = len(os.walk('./151113_iMM904_Testing/data/mbaCandRxns/%s_0' % state).next()[2])
	count1 = len(os.walk('./151113_iMM904_Testing/data/mbaCandRxns/%s_1' % state).next()[2])
	count2 = len(os.walk('./151113_iMM904_Testing/data/mbaCandRxns/%s_2' % state).next()[2])
	count3 = len(os.walk('./151113_iMM904_Testing/data/mbaCandRxns/%s_3' % state).next()[2])
	count4 = len(os.walk('./151113_iMM904_Testing/data/mbaCandRxns/%s_4' % state).next()[2])
	metabolicState_path0 = './151113_iMM904_Testing/data/metabolicState_%s_0.csv' % state
	metabolicState_path1 = './151113_iMM904_Testing/data/metabolicState_%s_1.csv' % state
	metabolicState_path2 = './151113_iMM904_Testing/data/metabolicState_%s_2.csv' % state
	metabolicState_path3 = './151113_iMM904_Testing/data/metabolicState_%s_3.csv' % state
	metabolicState_path4 = './151113_iMM904_Testing/data/metabolicState_%s_4.csv' % state
	if os.path.exists(metabolicState_path0):
		metabolicState0 = 1
	else:
		metabolicState0 = 0
	if os.path.exists(metabolicState_path1):
		metabolicState1 = 1
	else:
		metabolicState1 = 0	
	if os.path.exists(metabolicState_path2):
		metabolicState2 = 1
	else:
		metabolicState2 = 0
	if os.path.exists(metabolicState_path3):
		metabolicState3 = 1
	else:
		metabolicState3 = 0
	if os.path.exists(metabolicState_path4):
		metabolicState4 = 1
	else:
		metabolicState4 = 0
	mbaCandRxns_dict[state] = {'count0': count0, 'count1': count1, 'count2': count2, 'count3': count3, 'count4': count4, 'metabolicState0': metabolicState0, 'metabolicState1': metabolicState1, 'metabolicState2': metabolicState2, 'metabolicState3': metabolicState3, 'metabolicState4': metabolicState4}


mbaCandRxns_dict_dataset_S2 = {}
for state in metabolicState_list_dataset_S2:
	list_of_files = glob.glob('./dataset_S2/data/mbaCandRxns/*')
	count0 = 0
	count1 = 0
	count2 = 0
	count3 = 0
	count4 = 0
	for filename in list_of_files:
		if re.match('./dataset_S2/data/mbaCandRxns/mbaCandRxns_%s_0' % state, filename):
			count0 += 1
		if re.match('./dataset_S2/data/mbaCandRxns/mbaCandRxns_%s_1' % state, filename):
			count1 += 1
		if re.match('./dataset_S2/data/mbaCandRxns/mbaCandRxns_%s_2' % state, filename):
			count2 += 1
		if re.match('./dataset_S2/data/mbaCandRxns/mbaCandRxns_%s_3' % state, filename):
			count3 += 1
		if re.match('./dataset_S2/data/mbaCandRxns/mbaCandRxns_%s_4' % state, filename):
			count4 += 1			
	metabolicState_path0 = './dataset_S2/data/metabolicState_%s_0.csv' % state
	metabolicState_path1 = './dataset_S2/data/metabolicState_%s_1.csv' % state
	metabolicState_path2 = './dataset_S2/data/metabolicState_%s_2.csv' % state
	metabolicState_path3 = './dataset_S2/data/metabolicState_%s_3.csv' % state
	metabolicState_path4 = './dataset_S2/data/metabolicState_%s_4.csv' % state
	if os.path.exists(metabolicState_path0):
		metabolicState0 = 1
	else:
		metabolicState0 = 0
	if os.path.exists(metabolicState_path1):
		metabolicState1 = 1
	else:
		metabolicState1 = 0	
	if os.path.exists(metabolicState_path2):
		metabolicState2 = 1
	else:
		metabolicState2 = 0
	if os.path.exists(metabolicState_path3):
		metabolicState3 = 1
	else:
		metabolicState3 = 0
	if os.path.exists(metabolicState_path4):
		metabolicState4 = 1
	else:
		metabolicState4 = 0
	mbaCandRxns_dict_dataset_S2[state] = {'count0': count0, 'count1': count1, 'count2': count2, 'count3': count3, 'count4': count4, 'metabolicState0': metabolicState0, 'metabolicState1': metabolicState1, 'metabolicState2': metabolicState2, 'metabolicState3': metabolicState3, 'metabolicState4': metabolicState4}

mbaCandRxns_dict_dataset_S2_eps = {}
for state in metabolicState_list_dataset_S2_eps:
	list_of_files = glob.glob('./dataset_S2_eps/data/mbaCandRxns/*')
	count0 = 0
	count1 = 0
	count2 = 0
	count3 = 0
	count4 = 0
	for filename in list_of_files:
		if re.match('./dataset_S2_eps/data/mbaCandRxns/mbaCandRxns_%s_0' % state, filename):
			count0 += 1
		if re.match('./dataset_S2_eps/data/mbaCandRxns/mbaCandRxns_%s_1' % state, filename):
			count1 += 1
		if re.match('./dataset_S2_eps/data/mbaCandRxns/mbaCandRxns_%s_2' % state, filename):
			count2 += 1
		if re.match('./dataset_S2_eps/data/mbaCandRxns/mbaCandRxns_%s_3' % state, filename):
			count3 += 1
		if re.match('./dataset_S2_eps/data/mbaCandRxns/mbaCandRxns_%s_4' % state, filename):
			count4 += 1			
	metabolicState_path0 = './dataset_S2_eps/data/metabolicState_%s_0.csv' % state
	metabolicState_path1 = './dataset_S2_eps/data/metabolicState_%s_1.csv' % state
	metabolicState_path2 = './dataset_S2_eps/data/metabolicState_%s_2.csv' % state
	metabolicState_path3 = './dataset_S2_eps/data/metabolicState_%s_3.csv' % state
	metabolicState_path4 = './dataset_S2_eps/data/metabolicState_%s_4.csv' % state
	if os.path.exists(metabolicState_path0):
		metabolicState0 = 1
	else:
		metabolicState0 = 0
	if os.path.exists(metabolicState_path1):
		metabolicState1 = 1
	else:
		metabolicState1 = 0	
	if os.path.exists(metabolicState_path2):
		metabolicState2 = 1
	else:
		metabolicState2 = 0
	if os.path.exists(metabolicState_path3):
		metabolicState3 = 1
	else:
		metabolicState3 = 0
	if os.path.exists(metabolicState_path4):
		metabolicState4 = 1
	else:
		metabolicState4 = 0
	mbaCandRxns_dict_dataset_S2_eps[state] = {'count0': count0, 'count1': count1, 'count2': count2, 'count3': count3, 'count4': count4, 'metabolicState0': metabolicState0, 'metabolicState1': metabolicState1, 'metabolicState2': metabolicState2, 'metabolicState3': metabolicState3, 'metabolicState4': metabolicState4}

mbaCandRxns_dict_dataset_S2_original_model = {}
for state in metabolicState_list_dataset_S2_original_model:
	list_of_files = glob.glob('./dataset_S2_original_model/data/mbaCandRxns/*')
	count0 = 0
	count1 = 0
	count2 = 0
	count3 = 0
	count4 = 0
	for filename in list_of_files:
		if re.match('./dataset_S2_original_model/data/mbaCandRxns/mbaCandRxns_%s_0' % state, filename):
			count0 += 1
		if re.match('./dataset_S2_original_model/data/mbaCandRxns/mbaCandRxns_%s_1' % state, filename):
			count1 += 1
		if re.match('./dataset_S2_original_model/data/mbaCandRxns/mbaCandRxns_%s_2' % state, filename):
			count2 += 1
		if re.match('./dataset_S2_original_model/data/mbaCandRxns/mbaCandRxns_%s_3' % state, filename):
			count3 += 1
		if re.match('./dataset_S2_original_model/data/mbaCandRxns/mbaCandRxns_%s_4' % state, filename):
			count4 += 1			
	metabolicState_path0 = './dataset_S2_original_model/data/metabolicState_%s_0.csv' % state
	metabolicState_path1 = './dataset_S2_original_model/data/metabolicState_%s_1.csv' % state
	metabolicState_path2 = './dataset_S2_original_model/data/metabolicState_%s_2.csv' % state
	metabolicState_path3 = './dataset_S2_original_model/data/metabolicState_%s_3.csv' % state
	metabolicState_path4 = './dataset_S2_original_model/data/metabolicState_%s_4.csv' % state
	if os.path.exists(metabolicState_path0):
		metabolicState0 = 1
	else:
		metabolicState0 = 0
	if os.path.exists(metabolicState_path1):
		metabolicState1 = 1
	else:
		metabolicState1 = 0	
	if os.path.exists(metabolicState_path2):
		metabolicState2 = 1
	else:
		metabolicState2 = 0
	if os.path.exists(metabolicState_path3):
		metabolicState3 = 1
	else:
		metabolicState3 = 0
	if os.path.exists(metabolicState_path4):
		metabolicState4 = 1
	else:
		metabolicState4 = 0
	mbaCandRxns_dict_dataset_S2_original_model[state] = {'count0': count0, 'count1': count1, 'count2': count2, 'count3': count3, 'count4': count4, 'metabolicState0': metabolicState0, 'metabolicState1': metabolicState1, 'metabolicState2': metabolicState2, 'metabolicState3': metabolicState3, 'metabolicState4': metabolicState4}

mbaCandRxns_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_list_dataset_S2_original_model_eps:
	list_of_files = glob.glob('./dataset_S2_original_model_eps/data/mbaCandRxns/*')
	count0 = 0
	count1 = 0
	count2 = 0
	count3 = 0
	count4 = 0
	for filename in list_of_files:
		if re.match('./dataset_S2_original_model_eps/data/mbaCandRxns/mbaCandRxns_%s_0' % state, filename):
			count0 += 1
		if re.match('./dataset_S2_original_model_eps/data/mbaCandRxns/mbaCandRxns_%s_1' % state, filename):
			count1 += 1
		if re.match('./dataset_S2_original_model_eps/data/mbaCandRxns/mbaCandRxns_%s_2' % state, filename):
			count2 += 1
		if re.match('./dataset_S2_original_model_eps/data/mbaCandRxns/mbaCandRxns_%s_3' % state, filename):
			count3 += 1
		if re.match('./dataset_S2_original_model_eps/data/mbaCandRxns/mbaCandRxns_%s_4' % state, filename):
			count4 += 1			
	metabolicState_path0 = './dataset_S2_original_model_eps/data/metabolicState_%s_0.csv' % state
	metabolicState_path1 = './dataset_S2_original_model_eps/data/metabolicState_%s_1.csv' % state
	metabolicState_path2 = './dataset_S2_original_model_eps/data/metabolicState_%s_2.csv' % state
	metabolicState_path3 = './dataset_S2_original_model_eps/data/metabolicState_%s_3.csv' % state
	metabolicState_path4 = './dataset_S2_original_model_eps/data/metabolicState_%s_4.csv' % state
	if os.path.exists(metabolicState_path0):
		metabolicState0 = 1
	else:
		metabolicState0 = 0
	if os.path.exists(metabolicState_path1):
		metabolicState1 = 1
	else:
		metabolicState1 = 0	
	if os.path.exists(metabolicState_path2):
		metabolicState2 = 1
	else:
		metabolicState2 = 0
	if os.path.exists(metabolicState_path3):
		metabolicState3 = 1
	else:
		metabolicState3 = 0
	if os.path.exists(metabolicState_path4):
		metabolicState4 = 1
	else:
		metabolicState4 = 0
	mbaCandRxns_dict_dataset_S2_original_model_eps[state] = {'count0': count0, 'count1': count1, 'count2': count2, 'count3': count3, 'count4': count4, 'metabolicState0': metabolicState0, 'metabolicState1': metabolicState1, 'metabolicState2': metabolicState2, 'metabolicState3': metabolicState3, 'metabolicState4': metabolicState4}




metabolicState_dict = {'151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm':'iMM904_NADcorrected_1127_FTHFLm' , '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c',
'151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c',
'151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns':'iMM904_NADcorrected_1127_FTHFLm' , '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'}

metabolicState_dict_dataset_S2 = {'151006_Aerobic_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25': 'iMM904_NADcorrected_1127_FTHFLm'}

metabolicState_dict_dataset_S2_eps = {'151006_Aerobic_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25': 'iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25': 'iMM904_NADcorrected_1127_FTHFLm'}

metabolicState_dict_dataset_S2_original_model = {'151006_Aerobic_0.25': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25': 'iMM904_blkRxnsDeleted_dict'}

metabolicState_dict_dataset_S2_original_model_eps = {'151006_Aerobic_0.25': 'iMM904_blkRxnsDeleted_dict', '151006_Anaerobic_0.25': 'iMM904_blkRxnsDeleted_dict', '151012_ethanol_0.25': 'iMM904_blkRxnsDeleted_dict', '151012_glucose_0.25': 'iMM904_blkRxnsDeleted_dict', '151015_negative_control_0.25': 'iMM904_blkRxnsDeleted_dict'}

pickle_models = {}
md_models = {}
for model in pickle_models_list:
	pickle_models[model] = importPickle('./151113_iMM904_Testing/data/%s.pkl' % model)
	md_models[model] = CbModel(pickle_models[model]['S'], pickle_models[model]['idSp'], pickle_models[model]['idRs'], pickle_models[model]['lb'], pickle_models[model]['ub'], pickle_models[model]['rxns'], pickle_models[model]['genes'])

pickle_models_dataset_S2 = {}
md_models_dataset_S2 = {}
for model in pickle_models_list_dataset_S2:
	pickle_models_dataset_S2[model] = importPickle('./dataset_S2/data/%s.pkl' % model)
	md_models_dataset_S2[model] = CbModel(pickle_models_dataset_S2[model]['S'], pickle_models_dataset_S2[model]['idSp'], pickle_models_dataset_S2[model]['idRs'], pickle_models_dataset_S2[model]['lb'], pickle_models_dataset_S2[model]['ub'], pickle_models_dataset_S2[model]['rxns'], pickle_models_dataset_S2[model]['genes'])

pickle_models_dataset_S2_eps = {}
md_models_dataset_S2_eps = {}
for model in pickle_models_list_dataset_S2:
	pickle_models_dataset_S2_eps[model] = importPickle('./dataset_S2_eps/data/%s.pkl' % model)
	md_models_dataset_S2_eps[model] = CbModel(pickle_models_dataset_S2_eps[model]['S'], pickle_models_dataset_S2_eps[model]['idSp'], pickle_models_dataset_S2_eps[model]['idRs'], pickle_models_dataset_S2_eps[model]['lb'], pickle_models_dataset_S2_eps[model]['ub'], pickle_models_dataset_S2_eps[model]['rxns'], pickle_models_dataset_S2_eps[model]['genes'])

pickle_models_dataset_S2_original_model = {}
md_models_dataset_S2_original_model = {}
for model in pickle_models_list_dataset_S2_original_model:
	pickle_models_dataset_S2_original_model[model] = importPickle('./dataset_S2_original_model/data/%s.pkl' % model)
	md_models_dataset_S2_original_model[model] = CbModel(pickle_models_dataset_S2_original_model[model]['S'], pickle_models_dataset_S2_original_model[model]['idSp'], pickle_models_dataset_S2_original_model[model]['idRs'], pickle_models_dataset_S2_original_model[model]['lb'], pickle_models_dataset_S2_original_model[model]['ub'], pickle_models_dataset_S2_original_model[model]['rxns'], pickle_models_dataset_S2_original_model[model]['genes'])

pickle_models_dataset_S2_original_model_eps = {}
md_models_dataset_S2_original_model_eps = {}
for model in pickle_models_list_dataset_S2_original_model:
	pickle_models_dataset_S2_original_model_eps[model] = importPickle('./dataset_S2_original_model_eps/data/%s.pkl' % model)
	md_models_dataset_S2_original_model_eps[model] = CbModel(pickle_models_dataset_S2_original_model_eps[model]['S'], pickle_models_dataset_S2_original_model_eps[model]['idSp'], pickle_models_dataset_S2_original_model_eps[model]['idRs'], pickle_models_dataset_S2_original_model_eps[model]['lb'], pickle_models_dataset_S2_original_model_eps[model]['ub'], pickle_models_dataset_S2_original_model_eps[model]['rxns'], pickle_models_dataset_S2_original_model_eps[model]['genes'])

fbr_freqBased = {}
gbr_fOutRxnsByExpression = {}
for freq in freqBased_list:
	fbr_freqBased[freq] = importPickle('./151113_iMM904_Testing/data/freqBasedRxns_%s.pkl' % freq)
	gbr_fOutRxnsByExpression = importPickle('./151113_iMM904_Testing/data/rxnsClassifiedByExprssion_%s.pkl' %freq)

fbr_freqBased_dataset_S2 = {}
gbr_fOutRxnsByExpression_dataset_S2 = {}
for freq in freqBased_list_dataset_S2:
	try:
		fbr_freqBased_dataset_S2[freq] = importPickle('./dataset_S2/data/freqBasedRxns_%s.pkl' % freq)
	except IOError:
		fbr_freqBased_dataset_S2[freq] = 0
	try:
		gbr_fOutRxnsByExpression_dataset_S2 = importPickle('./dataset_S2/data/rxnsClassifiedByExprssion_%s.pkl' %freq)
	except IOError:
		fbr_freqBased_dataset_S2[freq] = 0

fbr_freqBased_dataset_S2_eps = {}
gbr_fOutRxnsByExpression_dataset_S2_eps = {}
for freq in freqBased_list_dataset_S2:
	try: 
		fbr_freqBased_dataset_S2_eps[freq] = importPickle('./dataset_S2_eps/data/freqBasedRxns_%s.pkl' % freq)
	except IOError:
		fbr_freqBased_dataset_S2_eps[freq] = 0
	try:
		gbr_fOutRxnsByExpression_dataset_S2_eps = importPickle('./dataset_S2_eps/data/rxnsClassifiedByExprssion_%s.pkl' %freq)
	except IOError:
		gbr_fOutRxnsByExpression_dataset_S2_eps = 0

fbr_freqBased_dataset_S2_original_model = {}
gbr_fOutRxnsByExpression_dataset_S2_original_model = {}
for freq in freqBased_list_dataset_S2_original_model:
	try:
		fbr_freqBased_dataset_S2_original_model[freq] = importPickle('./dataset_S2_original_model/data/freqBasedRxns_%s.pkl' % freq)
	except IOError:
		fbr_freqBased_dataset_S2_original_model[freq] = 0
	try:
		gbr_fOutRxnsByExpression_dataset_S2_original_model = importPickle('./dataset_S2_original_model/data/rxnsClassifiedByExprssion_%s.pkl' %freq)
	except IOError:
		gbr_fOutRxnsByExpression_dataset_S2_original_model = 0

fbr_freqBased_dataset_S2_original_model_eps = {}
gbr_fOutRxnsByExpression_dataset_S2_original_model_eps = {}
for freq in freqBased_list_dataset_S2_original_model:
	try:
		fbr_freqBased_dataset_S2_original_model_eps[freq] = importPickle('./dataset_S2_original_model_eps/data/freqBasedRxns_%s.pkl' % freq)
	except IOError:
		fbr_freqBased_dataset_S2_original_model_eps[freq] = 0
	try:
		gbr_fOutRxnsByExpression_dataset_S2_original_model_eps = importPickle('./dataset_S2_original_model_eps/data/rxnsClassifiedByExprssion_%s.pkl' %freq)
	except IOError: 
		gbr_fOutRxnsByExpression_dataset_S2_original_model_eps = 0

fluxdict = {}
fluxavgdict = {}
fluxvardict = {}
fluxstddict = {}
fluxavgdict_text = {}
rxninallreps = {}
rxninnonereps = {}
idRsfluxstate = {}
allidRsfluxstate = {}
for state in metabolicState_list:
	metabolicState_file_1 = open('./151113_iMM904_Testing/data/metabolicState_%s_0.csv' % state, 'r')
	metabolicState_file_2 = open('./151113_iMM904_Testing/data/metabolicState_%s_1.csv' % state, 'r')
	metabolicState_file_3 = open('./151113_iMM904_Testing/data/metabolicState_%s_2.csv' % state, 'r')
	metabolicState_file_4 = open('./151113_iMM904_Testing/data/metabolicState_%s_3.csv' % state, 'r')
	metabolicState_file_5 = open('./151113_iMM904_Testing/data/metabolicState_%s_4.csv' % state, 'r')
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
	cond1fluxavg = {}
	rxninallreps_list = []
	rxninnonereps_list = []
	for t in pickle_models[metabolicState_dict[state]]['rxns']:
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
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
			if not fluxdict:
				fluxdict[state] = {t: [0,0,0,0,0]}
			else:
				if state not in fluxdict:
					fluxdict[state] = {t: [0,0,0,0,0]}
				else:
					fluxdict[state][t] = [0,0,0,0,0]
			fluxavgdict1[t] = 0
			fluxvardict1[t] = 0
			fluxstddict1[t] = 0
			rxninnonereps_list.append(t) 
		else:
			fluxavgdict1[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			if not fluxdict:
				fluxdict[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
			else:
				if state not in fluxdict:
					fluxdict[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
				else:
					fluxdict[state][t] = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			fluxvardict1[t] = round(np.var(fluxes),4)
			fluxstddict1[t] = round(np.std(fluxes),4)
			fluxavg = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='e_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
				fluxtext.append("_")
				if round(np.std(fluxes),4) == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxes),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='e_4'
					else:
						std = re.sub('\.','',std)
						std+='e_4'
					fluxtext.append(std)		
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]

	fluxavgdict[state] = fluxavgdict1
	fluxvardict[state] = fluxvardict1
	fluxstddict[state] = fluxstddict1
	fluxavgdict_text[state] = cond1fluxavg
	rxninallreps[state] = rxninallreps_list
	rxninnonereps[state] = rxninnonereps_list
	idRsfluxstate[state] = {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}
	allidRsfluxstate[state] = {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}


fluxdict_dataset_S2 = {}
fluxavgdict_dataset_S2 = {}
fluxvardict_dataset_S2 = {}
fluxstddict_dataset_S2 = {}
fluxavgdict_text_dataset_S2 = {}
rxninallreps_dataset_S2 = {}
rxninnonereps_dataset_S2 = {}
idRsfluxstate_dataset_S2 = {}
allidRsfluxstate_dataset_S2 = {}
for state in metabolicState_list_dataset_S2:
	try:
		metabolicState_file_1 = open('./dataset_S2/data/metabolicState_%s_0.csv' % state, 'r')
		csvreader1 = csv.reader(metabolicState_file_1)

		idRsflux1 = {}
		idRscsv1 = []
		flux1 = []

		for row in csvreader1:
			idRscsv1.append(row[0])
			flux1.append(float(row[1]))
		for i, item1 in enumerate(idRscsv1):
			idRsflux1[item1] = flux1[i]
	except IOError:
		idRsflux1 = {}
	try: 
		metabolicState_file_2 = open('./dataset_S2/data/metabolicState_%s_1.csv' % state, 'r')
		csvreader2 = csv.reader(metabolicState_file_2)

		idRsflux2 = {}
		idRscsv2 = []
		flux2 = []

		for row in csvreader2:
			idRscsv2.append(row[0])
			flux2.append(float(row[1]))
		for i, item1 in enumerate(idRscsv2):
			idRsflux2[item1] = flux2[i]
	except IOError:
		idRsflux2 = {} 
	try:	
		metabolicState_file_3 = open('./dataset_S2/data/metabolicState_%s_2.csv' % state, 'r')
		csvreader3 = csv.reader(metabolicState_file_3)

		idRsflux3 = {}
		idRscsv3 = []
		flux3 = []

		for row in csvreader3:
			idRscsv3.append(row[0])
			flux3.append(float(row[1]))
		for i, item1 in enumerate(idRscsv3):
			idRsflux3[item1] = flux3[i]
	except IOError:
		idRsflux3 = {} 
	try:	
		metabolicState_file_4 = open('./dataset_S2/data/metabolicState_%s_3.csv' % state, 'r')
		csvreader4 = csv.reader(metabolicState_file_4)

		idRsflux4 = {}
		idRscsv4 = []
		flux4 = []

		for row in csvreader4:
			idRscsv4.append(row[0])
			flux4.append(float(row[1]))
		for i, item1 in enumerate(idRscsv4):
			idRsflux4[item1] = flux4[i]
	except IOError:
		idRsflux4 = {}
	try:	
		metabolicState_file_5 = open('./dataset_S2/data/metabolicState_%s_4.csv' % state, 'r')
		csvreader5 = csv.reader(metabolicState_file_5)

		idRsflux5 = {}
		idRscsv5 = []
		flux5 = []

		for row in csvreader5:
			idRscsv5.append(row[0])
			flux5.append(float(row[1]))
		for i, item1 in enumerate(idRscsv5):
			idRsflux5[item1] = flux5[i]
	except IOError:
		idRsflux5 = {}


	allidRsflux1 = {}
	allidRsflux2 = {}
	allidRsflux3 = {}
	allidRsflux4 = {}
	allidRsflux5 = {}
	fluxavgdict1 = {}
	fluxvardict1 = {}
	fluxstddict1 = {}
	cond1fluxavg = {}
	rxninallreps_list = []
	rxninnonereps_list = []
	for t in pickle_models_dataset_S2[metabolicState_dict_dataset_S2[state]]['rxns']:
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
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
			if not fluxdict_dataset_S2:
				fluxdict_dataset_S2[state] = {t: [0,0,0,0,0]}
			else:
				if state not in fluxdict_dataset_S2:
					fluxdict_dataset_S2[state] = {t: [0,0,0,0,0]}
				else:
					fluxdict_dataset_S2[state][t] = [0,0,0,0,0]
			fluxavgdict1[t] = 0
			fluxvardict1[t] = 0
			fluxstddict1[t] = 0
			rxninnonereps_list.append(t) 
		else:
			fluxavgdict1[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			if not fluxdict_dataset_S2:
				fluxdict_dataset_S2[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
			else:
				if state not in fluxdict_dataset_S2:
					fluxdict_dataset_S2[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
				else:
					fluxdict_dataset_S2[state][t] = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			fluxvardict1[t] = round(np.var(fluxes),4)
			fluxstddict1[t] = round(np.std(fluxes),4)
			fluxavg = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='e_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
				fluxtext.append("_")
				if round(np.std(fluxes),4) == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxes),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='e_4'
					else:
						std = re.sub('\.','',std)
						std+='e_4'
					fluxtext.append(std)		
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]

	fluxavgdict_dataset_S2[state] = fluxavgdict1
	fluxvardict_dataset_S2[state] = fluxvardict1
	fluxstddict_dataset_S2[state] = fluxstddict1
	fluxavgdict_text_dataset_S2[state] = cond1fluxavg
	rxninallreps_dataset_S2[state] = rxninallreps_list
	rxninnonereps_dataset_S2[state] = rxninnonereps_list
	idRsfluxstate_dataset_S2[state] = {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}
	allidRsfluxstate_dataset_S2[state] = {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}


fluxdict_dataset_S2_eps = {}
fluxavgdict_dataset_S2_eps = {}
fluxvardict_dataset_S2_eps = {}
fluxstddict_dataset_S2_eps = {}
fluxavgdict_text_dataset_S2_eps = {}
rxninallreps_dataset_S2_eps = {}
rxninnonereps_dataset_S2_eps = {}
idRsfluxstate_dataset_S2_eps = {}
allidRsfluxstate_dataset_S2_eps = {}
for state in metabolicState_list_dataset_S2_eps:
	try:
		metabolicState_file_1 = open('./dataset_S2_eps/data/metabolicState_%s_0.csv' % state, 'r')
		csvreader1 = csv.reader(metabolicState_file_1)

		idRsflux1 = {}
		idRscsv1 = []
		flux1 = []

		for row in csvreader1:
			idRscsv1.append(row[0])
			flux1.append(float(row[1]))
		for i, item1 in enumerate(idRscsv1):
			idRsflux1[item1] = flux1[i]
	except IOError:
		idRsflux1 = {}
	try: 
		metabolicState_file_2 = open('./dataset_S2_eps/data/metabolicState_%s_1.csv' % state, 'r')
		csvreader2 = csv.reader(metabolicState_file_2)

		idRsflux2 = {}
		idRscsv2 = []
		flux2 = []

		for row in csvreader2:
			idRscsv2.append(row[0])
			flux2.append(float(row[1]))
		for i, item1 in enumerate(idRscsv2):
			idRsflux2[item1] = flux2[i]
	except IOError:
		idRsflux2 = {} 
	try:	
		metabolicState_file_3 = open('./dataset_S2_eps/data/metabolicState_%s_2.csv' % state, 'r')
		csvreader3 = csv.reader(metabolicState_file_3)

		idRsflux3 = {}
		idRscsv3 = []
		flux3 = []

		for row in csvreader3:
			idRscsv3.append(row[0])
			flux3.append(float(row[1]))
		for i, item1 in enumerate(idRscsv3):
			idRsflux3[item1] = flux3[i]
	except IOError:
		idRsflux3 = {} 
	try:	
		metabolicState_file_4 = open('./dataset_S2_eps/data/metabolicState_%s_3.csv' % state, 'r')
		csvreader4 = csv.reader(metabolicState_file_4)

		idRsflux4 = {}
		idRscsv4 = []
		flux4 = []

		for row in csvreader4:
			idRscsv4.append(row[0])
			flux4.append(float(row[1]))
		for i, item1 in enumerate(idRscsv4):
			idRsflux4[item1] = flux4[i]
	except IOError:
		idRsflux4 = {}
	try:	
		metabolicState_file_5 = open('./dataset_S2_eps/data/metabolicState_%s_4.csv' % state, 'r')
		csvreader5 = csv.reader(metabolicState_file_5)

		idRsflux5 = {}
		idRscsv5 = []
		flux5 = []

		for row in csvreader5:
			idRscsv5.append(row[0])
			flux5.append(float(row[1]))
		for i, item1 in enumerate(idRscsv5):
			idRsflux5[item1] = flux5[i]
	except IOError:
		idRsflux5 = {}


	allidRsflux1 = {}
	allidRsflux2 = {}
	allidRsflux3 = {}
	allidRsflux4 = {}
	allidRsflux5 = {}
	fluxavgdict1 = {}
	fluxvardict1 = {}
	fluxstddict1 = {}
	cond1fluxavg = {}
	rxninallreps_list = []
	rxninnonereps_list = []
	for t in pickle_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]]['rxns']:
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
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
			if not fluxdict_dataset_S2_eps:
				fluxdict_dataset_S2_eps[state] = {t: [0,0,0,0,0]}
			else:
				if state not in fluxdict_dataset_S2_eps:
					fluxdict_dataset_S2_eps[state] = {t: [0,0,0,0,0]}
				else:
					fluxdict_dataset_S2_eps[state][t] = [0,0,0,0,0]
			fluxavgdict1[t] = 0
			fluxvardict1[t] = 0
			fluxstddict1[t] = 0
			rxninnonereps_list.append(t) 
		else:
			fluxavgdict1[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			if not fluxdict_dataset_S2_eps:
				fluxdict_dataset_S2_eps[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
			else:
				if state not in fluxdict_dataset_S2_eps:
					fluxdict_dataset_S2_eps[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
				else:
					fluxdict_dataset_S2_eps[state][t] = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			fluxvardict1[t] = round(np.var(fluxes),4)
			fluxstddict1[t] = round(np.std(fluxes),4)
			fluxavg = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='e_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
				fluxtext.append("_")
				if round(np.std(fluxes),4) == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxes),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='e_4'
					else:
						std = re.sub('\.','',std)
						std+='e_4'
					fluxtext.append(std)		
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]

	fluxavgdict_dataset_S2_eps[state] = fluxavgdict1
	fluxvardict_dataset_S2_eps[state] = fluxvardict1
	fluxstddict_dataset_S2_eps[state] = fluxstddict1
	fluxavgdict_text_dataset_S2_eps[state] = cond1fluxavg
	rxninallreps_dataset_S2_eps[state] = rxninallreps_list
	rxninnonereps_dataset_S2_eps[state] = rxninnonereps_list
	idRsfluxstate_dataset_S2_eps[state] = {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}
	allidRsfluxstate_dataset_S2_eps[state] = {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}


fluxdict_dataset_S2_original_model = {}
fluxavgdict_dataset_S2_original_model = {}
fluxvardict_dataset_S2_original_model = {}
fluxstddict_dataset_S2_original_model = {}
fluxavgdict_text_dataset_S2_original_model = {}
rxninallreps_dataset_S2_original_model = {}
rxninnonereps_dataset_S2_original_model = {}
idRsfluxstate_dataset_S2_original_model = {}
allidRsfluxstate_dataset_S2_original_model = {}
for state in metabolicState_list_dataset_S2_original_model:
	try:
		metabolicState_file_1 = open('./dataset_S2_original_model/data/metabolicState_%s_0.csv' % state, 'r')
		csvreader1 = csv.reader(metabolicState_file_1)

		idRsflux1 = {}
		idRscsv1 = []
		flux1 = []

		for row in csvreader1:
			idRscsv1.append(row[0])
			flux1.append(float(row[1]))
		for i, item1 in enumerate(idRscsv1):
			idRsflux1[item1] = flux1[i]
	except IOError:
		idRsflux1 = {}
	try: 
		metabolicState_file_2 = open('./dataset_S2_original_model/data/metabolicState_%s_1.csv' % state, 'r')
		csvreader2 = csv.reader(metabolicState_file_2)

		idRsflux2 = {}
		idRscsv2 = []
		flux2 = []

		for row in csvreader2:
			idRscsv2.append(row[0])
			flux2.append(float(row[1]))
		for i, item1 in enumerate(idRscsv2):
			idRsflux2[item1] = flux2[i]
	except IOError:
		idRsflux2 = {} 
	try:	
		metabolicState_file_3 = open('./dataset_S2_original_model/data/metabolicState_%s_2.csv' % state, 'r')
		csvreader3 = csv.reader(metabolicState_file_3)

		idRsflux3 = {}
		idRscsv3 = []
		flux3 = []

		for row in csvreader3:
			idRscsv3.append(row[0])
			flux3.append(float(row[1]))
		for i, item1 in enumerate(idRscsv3):
			idRsflux3[item1] = flux3[i]
	except IOError:
		idRsflux3 = {} 
	try:	
		metabolicState_file_4 = open('./dataset_S2_original_model/data/metabolicState_%s_3.csv' % state, 'r')
		csvreader4 = csv.reader(metabolicState_file_4)

		idRsflux4 = {}
		idRscsv4 = []
		flux4 = []

		for row in csvreader4:
			idRscsv4.append(row[0])
			flux4.append(float(row[1]))
		for i, item1 in enumerate(idRscsv4):
			idRsflux4[item1] = flux4[i]
	except IOError:
		idRsflux4 = {}
	try:	
		metabolicState_file_5 = open('./dataset_S2_original_model/data/metabolicState_%s_4.csv' % state, 'r')
		csvreader5 = csv.reader(metabolicState_file_5)

		idRsflux5 = {}
		idRscsv5 = []
		flux5 = []

		for row in csvreader5:
			idRscsv5.append(row[0])
			flux5.append(float(row[1]))
		for i, item1 in enumerate(idRscsv5):
			idRsflux5[item1] = flux5[i]
	except IOError:
		idRsflux5 = {}


	allidRsflux1 = {}
	allidRsflux2 = {}
	allidRsflux3 = {}
	allidRsflux4 = {}
	allidRsflux5 = {}
	fluxavgdict1 = {}
	fluxvardict1 = {}
	fluxstddict1 = {}
	cond1fluxavg = {}
	rxninallreps_list = []
	rxninnonereps_list = []
	for t in pickle_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]]['rxns']:
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
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
			if not fluxdict_dataset_S2_original_model:
				fluxdict_dataset_S2_original_model[state] = {t: [0,0,0,0,0]}
			else:
				if state not in fluxdict_dataset_S2_original_model:
					fluxdict_dataset_S2_original_model[state] = {t: [0,0,0,0,0]}
				else:
					fluxdict_dataset_S2_original_model[state][t] = [0,0,0,0,0]
			fluxavgdict1[t] = 0
			fluxvardict1[t] = 0
			fluxstddict1[t] = 0
			rxninnonereps_list.append(t) 
		else:
			fluxavgdict1[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			if not fluxdict_dataset_S2_original_model:
				fluxdict_dataset_S2_original_model[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
			else:
				if state not in fluxdict_dataset_S2_original_model:
					fluxdict_dataset_S2_original_model[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
				else:
					fluxdict_dataset_S2_original_model[state][t] = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			fluxvardict1[t] = round(np.var(fluxes),4)
			fluxstddict1[t] = round(np.std(fluxes),4)
			fluxavg = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='e_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
				fluxtext.append("_")
				if round(np.std(fluxes),4) == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxes),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='e_4'
					else:
						std = re.sub('\.','',std)
						std+='e_4'
					fluxtext.append(std)		
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]

	fluxavgdict_dataset_S2_original_model[state] = fluxavgdict1
	fluxvardict_dataset_S2_original_model[state] = fluxvardict1
	fluxstddict_dataset_S2_original_model[state] = fluxstddict1
	fluxavgdict_text_dataset_S2_original_model[state] = cond1fluxavg
	rxninallreps_dataset_S2_original_model[state] = rxninallreps_list
	rxninnonereps_dataset_S2_original_model[state] = rxninnonereps_list
	idRsfluxstate_dataset_S2_original_model[state] = {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}
	allidRsfluxstate_dataset_S2_original_model[state] = {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}


fluxdict_dataset_S2_original_model_eps = {}
fluxavgdict_dataset_S2_original_model_eps = {}
fluxvardict_dataset_S2_original_model_eps = {}
fluxstddict_dataset_S2_original_model_eps = {}
fluxavgdict_text_dataset_S2_original_model_eps = {}
rxninallreps_dataset_S2_original_model_eps = {}
rxninnonereps_dataset_S2_original_model_eps = {}
idRsfluxstate_dataset_S2_original_model_eps = {}
allidRsfluxstate_dataset_S2_original_model_eps = {}
for state in metabolicState_list_dataset_S2_original_model_eps:
	try:
		metabolicState_file_1 = open('./dataset_S2_original_model_eps/data/metabolicState_%s_0.csv' % state, 'r')
		csvreader1 = csv.reader(metabolicState_file_1)

		idRsflux1 = {}
		idRscsv1 = []
		flux1 = []

		for row in csvreader1:
			idRscsv1.append(row[0])
			flux1.append(float(row[1]))
		for i, item1 in enumerate(idRscsv1):
			idRsflux1[item1] = flux1[i]
	except IOError:
		idRsflux1 = {}
	try: 
		metabolicState_file_2 = open('./dataset_S2_original_model_eps/data/metabolicState_%s_1.csv' % state, 'r')
		csvreader2 = csv.reader(metabolicState_file_2)

		idRsflux2 = {}
		idRscsv2 = []
		flux2 = []

		for row in csvreader2:
			idRscsv2.append(row[0])
			flux2.append(float(row[1]))
		for i, item1 in enumerate(idRscsv2):
			idRsflux2[item1] = flux2[i]
	except IOError:
		idRsflux2 = {} 
	try:	
		metabolicState_file_3 = open('./dataset_S2_original_model_eps/data/metabolicState_%s_2.csv' % state, 'r')
		csvreader3 = csv.reader(metabolicState_file_3)

		idRsflux3 = {}
		idRscsv3 = []
		flux3 = []

		for row in csvreader3:
			idRscsv3.append(row[0])
			flux3.append(float(row[1]))
		for i, item1 in enumerate(idRscsv3):
			idRsflux3[item1] = flux3[i]
	except IOError:
		idRsflux3 = {} 
	try:	
		metabolicState_file_4 = open('./dataset_S2_original_model_eps/data/metabolicState_%s_3.csv' % state, 'r')
		csvreader4 = csv.reader(metabolicState_file_4)

		idRsflux4 = {}
		idRscsv4 = []
		flux4 = []

		for row in csvreader4:
			idRscsv4.append(row[0])
			flux4.append(float(row[1]))
		for i, item1 in enumerate(idRscsv4):
			idRsflux4[item1] = flux4[i]
	except IOError:
		idRsflux4 = {}
	try:	
		metabolicState_file_5 = open('./dataset_S2_original_model_eps/data/metabolicState_%s_4.csv' % state, 'r')
		csvreader5 = csv.reader(metabolicState_file_5)

		idRsflux5 = {}
		idRscsv5 = []
		flux5 = []

		for row in csvreader5:
			idRscsv5.append(row[0])
			flux5.append(float(row[1]))
		for i, item1 in enumerate(idRscsv5):
			idRsflux5[item1] = flux5[i]
	except IOError:
		idRsflux5 = {}


	allidRsflux1 = {}
	allidRsflux2 = {}
	allidRsflux3 = {}
	allidRsflux4 = {}
	allidRsflux5 = {}
	fluxavgdict1 = {}
	fluxvardict1 = {}
	fluxstddict1 = {}
	cond1fluxavg = {}
	rxninallreps_list = []
	rxninnonereps_list = []
	for t in pickle_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]]['rxns']:
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
			fluxtext = [t]
			fluxtext.append("_NA")
			fluxtextjoin = ''.join(fluxtext)
			cond1fluxavg[t] = fluxtextjoin[2:]
			if not fluxdict_dataset_S2_original_model_eps:
				fluxdict_dataset_S2_original_model_eps[state] = {t: [0,0,0,0,0]}
			else:
				if state not in fluxdict_dataset_S2_original_model_eps:
					fluxdict_dataset_S2_original_model_eps[state] = {t: [0,0,0,0,0]}
				else:
					fluxdict_dataset_S2_original_model_eps[state][t] = [0,0,0,0,0]
			fluxavgdict1[t] = 0
			fluxvardict1[t] = 0
			fluxstddict1[t] = 0
			rxninnonereps_list.append(t) 
		else:
			fluxavgdict1[t] = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			if not fluxdict_dataset_S2_original_model_eps:
				fluxdict_dataset_S2_original_model_eps[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
			else:
				if state not in fluxdict_dataset_S2_original_model_eps:
					fluxdict_dataset_S2_original_model_eps[state] = {t: [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]}
				else:
					fluxdict_dataset_S2_original_model_eps[state][t] = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
			fluxvardict1[t] = round(np.var(fluxes),4)
			fluxstddict1[t] = round(np.std(fluxes),4)
			fluxavg = (round((allidRsflux1[t] + allidRsflux2[t] + allidRsflux3[t] + allidRsflux4[t] + allidRsflux5[t])/5,4))
			fluxavg = str("%.4f" % fluxavg)
			negative_count = 0
			fluxavg_search = re.search('\-', fluxavg)
			if (fluxavg_search is not None):
				fluxavg = fluxavg[1:]
				negative_count = 1
			if fluxavg == "0.0000":
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtext.append("_")
				fluxtext.append("0")
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]
			else:
				if (len(fluxavg) == 6):
					fluxavgsearch = re.search('0\.', fluxavg)
					if (fluxavgsearch is not None):
						fluxavg = fluxavg[2:]
						fluxavg+='e_4'
					else:
						fluxavg = re.sub('\.','',fluxavg)
						fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				else:
					fluxavg = re.sub('\.','',fluxavg)
					fluxavg+='e_4'
					if negative_count == 1:
						fluxavg = "_" + fluxavg
				fluxtext = [t]
				fluxtext.append("_")
				fluxtext.append(str(fluxavg))
				fluxes = [allidRsflux1[t], allidRsflux2[t], allidRsflux3[t], allidRsflux4[t], allidRsflux5[t]]
				fluxtext.append("_")
				if round(np.std(fluxes),4) == 0:
					fluxtext.append("0")
				else:
					std = round(np.std(fluxes),4)
					std = str("%.4f" % std)
					if (len(std) == 6):
						stdsearch = re.search('0\.', std)
						if (stdsearch is not None):
							std = std[2:]
							std+='e_4'
					else:
						std = re.sub('\.','',std)
						std+='e_4'
					fluxtext.append(std)		
				fluxtextjoin = ''.join(fluxtext)
				cond1fluxavg[t] = fluxtextjoin[2:]

	fluxavgdict_dataset_S2_original_model_eps[state] = fluxavgdict1
	fluxvardict_dataset_S2_original_model_eps[state] = fluxvardict1
	fluxstddict_dataset_S2_original_model_eps[state] = fluxstddict1
	fluxavgdict_text_dataset_S2_original_model_eps[state] = cond1fluxavg
	rxninallreps_dataset_S2_original_model_eps[state] = rxninallreps_list
	rxninnonereps_dataset_S2_original_model_eps[state] = rxninnonereps_list
	idRsfluxstate_dataset_S2_original_model_eps[state] = {'idRsflux1': idRsflux1, 'idRsflux2': idRsflux2, 'idRsflux3': idRsflux3, 'idRsflux4': idRsflux4, 'idRsflux5': idRsflux5}
	allidRsfluxstate_dataset_S2_original_model_eps[state] = {'allidRsflux1': allidRsflux1, 'allidRsflux2': allidRsflux2, 'allidRsflux3': allidRsflux3, 'allidRsflxu4': allidRsflux4, 'allidRsflux5': allidRsflux5}



#import list of tested genes for YPD
YPD_tested_genes_list = []
YPD_tested_genes_file = open('Giaver_Tested_Genes.txt')
YPD = YPD_tested_genes_file.readlines()
for i, item1 in enumerate(YPD):
	YPD_tested_genes_list.append(item1.strip())

#import list of essential genes for YPD
YPD_essential_genes_list = []
YPD_genes_file = open('Giaver_YPD_Essential_Genes.txt')
YPD = YPD_genes_file.readlines()
for i, item1 in enumerate(YPD):
	YPD_essential_genes_list.append(item1.strip())

#import list of tested genes for YPEtoh
YPEtoh_tested_genes_list = []
YPEtoh_tested_genes_file = open('Snitkin_Tested_Mutants.txt')
YPEtoh_tested_genes = YPEtoh_tested_genes_file.readlines()
for i, item1 in enumerate(YPEtoh_tested_genes):
	YPEtoh_tested_genes_list.append(item1.strip())

#import list of essential genes for YPEtoh
YPEtoh_essential_genes_list = []
YPEtoh_genes_file = open('Snitkin_YPEtoh_Essential_Genes.txt')
YPEtoh = YPEtoh_genes_file.readlines()
for i, item1 in enumerate(YPEtoh):
	YPEtoh_essential_genes_list.append(item1.strip())

#determine the possible genes for each model
genes_possible_model_dict = {}
for state in metabolicState_list:
	gene_list = []
	for rxn in md_models[metabolicState_dict[state]].gene2rxn.keys():
		g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
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

genes_possible_model_dict_dataset_S2 = {}
for state in metabolicState_list_dataset_S2:
	gene_list = []
	for rxn in md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn.keys():
		g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
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
	genes_possible_model_dict_dataset_S2[state] = gene_list

genes_possible_model_dict_dataset_S2_eps = {}
for state in metabolicState_list_dataset_S2_eps:
	gene_list = []
	for rxn in md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn.keys():
		g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
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
	genes_possible_model_dict_dataset_S2_eps[state] = gene_list

genes_possible_model_dict_dataset_S2_original_model = {}
for state in metabolicState_list_dataset_S2_original_model:
	gene_list = []
	for rxn in md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn.keys():
		g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
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
	genes_possible_model_dict_dataset_S2_original_model[state] = gene_list

genes_possible_model_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_list_dataset_S2_original_model_eps:
	gene_list = []
	for rxn in md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn.keys():
		g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
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
	genes_possible_model_dict_dataset_S2_original_model_eps[state] = gene_list

#determine the genes never present in any of the pruned models for each state
genes_never_present_model_dict = {}
for state in metabolicState_list:
	gene_list = []
	for rxn in rxninnonereps[state]:
		if rxn in md_models[metabolicState_dict[state]].gene2rxn.keys():
			g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
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
	genes_never_present_model_dict[state] = gene_list

genes_never_present_model_dict_dataset_S2 = {}
for state in metabolicState_list_dataset_S2:
	gene_list = []
	for rxn in rxninnonereps_dataset_S2[state]:
		if rxn in md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
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
	genes_never_present_model_dict_dataset_S2[state] = gene_list

genes_never_present_model_dict_dataset_S2_eps = {}
for state in metabolicState_list_dataset_S2_eps:
	gene_list = []
	for rxn in rxninnonereps_dataset_S2_eps[state]:
		if rxn in md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
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
	genes_never_present_model_dict_dataset_S2_eps[state] = gene_list

genes_never_present_model_dict_dataset_S2_original_model = {}
for state in metabolicState_list_dataset_S2_original_model:
	gene_list = []
	for rxn in rxninnonereps_dataset_S2_original_model[state]:
		if rxn in md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
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
	genes_never_present_model_dict_dataset_S2_original_model[state] = gene_list

genes_never_present_model_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_list_dataset_S2_original_model_eps:
	gene_list = []
	for rxn in rxninnonereps_dataset_S2_original_model_eps[state]:
		if rxn in md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
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
	genes_never_present_model_dict_dataset_S2_original_model_eps[state] = gene_list					

#determine the genes always present in each of the pruned models for each state
genes_always_present_model_dict = {}
for state in metabolicState_list:
	gene_list = []
	for rxn in rxninallreps[state]:
		if rxn in md_models[metabolicState_dict[state]].gene2rxn.keys():
			g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
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
	genes_always_present_model_dict[state] = gene_list

genes_always_present_model_dict_dataset_S2 = {}
for state in metabolicState_list_dataset_S2:
	gene_list = []
	for rxn in rxninallreps_dataset_S2[state]:
		if rxn in md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
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
	genes_always_present_model_dict_dataset_S2[state] = gene_list

genes_always_present_model_dict_dataset_S2_eps = {}
for state in metabolicState_list_dataset_S2_eps:
	gene_list = []
	for rxn in rxninallreps_dataset_S2_eps[state]:
		if rxn in md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
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
	genes_always_present_model_dict_dataset_S2_eps[state] = gene_list

genes_always_present_model_dict_dataset_S2_original_model = {}
for state in metabolicState_list_dataset_S2_original_model:
	gene_list = []
	for rxn in rxninallreps_dataset_S2_original_model[state]:
		if rxn in md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
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
	genes_always_present_model_dict_dataset_S2_original_model[state] = gene_list

genes_always_present_model_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_list_dataset_S2_original_model_eps:
	gene_list = []
	for rxn in rxninallreps_dataset_S2_original_model_eps[state]:
		if rxn in md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn.keys():
			g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
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
	genes_always_present_model_dict_dataset_S2_original_model_eps[state] = gene_list					

#define the metabolic states for the ethanol growth condition (and negative control)
metabolicState_ethanol_list = ['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c',  '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', 
'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']

metabolicState_ethanol_list_dataset_S2 = ['151012_ethanol_0.25', '151015_negative_control_0.25']

metabolicState_ethanol_list_dataset_S2_eps = ['151012_ethanol_0.25', '151015_negative_control_0.25']

metabolicState_ethanol_list_dataset_S2_original_model = ['151012_ethanol_0.25', '151015_negative_control_0.25']

metabolicState_ethanol_list_dataset_S2_original_model_eps = ['151012_ethanol_0.25', '151015_negative_control_0.25']

#define the metabolic states for the glucose growth condition (and negative control)
metabolicState_glucose_list = ['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c',  '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']

metabolicState_glucose_list_dataset_S2 = ['151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_glucose_list_dataset_S2_eps = ['151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_glucose_list_dataset_S2_original_model = ['151012_glucose_0.25', '151015_negative_control_0.25']

metabolicState_glucose_list_dataset_S2_original_model_eps = ['151012_glucose_0.25', '151015_negative_control_0.25']

#determine the number of negative and positive genes present for the ethanol data
ethanol_positive_negative_dict = {}
for state in metabolicState_ethanol_list:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate[state][flux_state].keys():
			g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
			g2r = str(g2r)
			i = g2r.split('or')
			for count, j in enumerate(i):
				j = j.split('and')
				for k in j:
					k = k.translate(None, ' ()')
					if not k:
						continue
					else:
						if k in YPEtoh_tested_genes_list:
							#Tested with YPD background experimentally. Therefore exclude YPD essential genes
							if k not in YPD_essential_genes_list:
								if k in YPEtoh_essential_genes_list:
									if k not in gene_list_true_positive:
										gene_list_true_positive.append(k)
								if k not in YPEtoh_essential_genes_list:
									if k not in gene_list_false_positive:
										gene_list_false_positive.append(k)
						
		gene_list_true_negative = []
		gene_list_false_negative = []
		for rxn in md_models[metabolicState_dict[state]].gene2rxn.keys():
			if rxn not in idRsfluxstate[state][flux_state].keys():
				g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
				g2r = str(g2r)
				i = g2r.split('or')
				for count, j in enumerate(i):
					j = j.split('and')
					for k in j:
						k = k.translate(None, ' ()')
						if not k:
							continue
						else:
							if k in YPEtoh_tested_genes_list:
								#Tested with YPD background experimentally. Therefore exclude YPD essential genes
								if k not in YPD_essential_genes_list:
									if k not in YPEtoh_essential_genes_list:
										if k not in gene_list_true_negative:
											gene_list_true_negative.append(k)
									if k in YPEtoh_essential_genes_list:
										if k not in gene_list_false_negative:
											gene_list_false_negative.append(k)

		idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
	sensitivity_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))]
	ppv_list = [float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive'])), float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive'])),float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive'])),float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive'])),float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))]
	specificity_list = [float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive'])), float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive'])), float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive'])), float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive'])), float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))]
	npv_list = [float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative'])), float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative'])), float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative'])), float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative'])), float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))]
	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	ethanol_positive_negative_dict[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}



ethanol_positive_negative_dict_dataset_S2 = {}
for state in metabolicState_ethanol_list_dataset_S2:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2[state][flux_state].keys():
			g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
			g2r = str(g2r)
			i = g2r.split('or')
			for count, j in enumerate(i):
				j = j.split('and')
				for k in j:
					k = k.translate(None, ' ()')
					if not k:
						continue
					else:
						if k in YPEtoh_tested_genes_list:
							#Tested with YPD background experimentally. Therefore exclude YPD essential genes
							if k not in YPD_essential_genes_list:
								if k in YPEtoh_essential_genes_list:
									if k not in gene_list_true_positive:
										gene_list_true_positive.append(k)
								if k not in YPEtoh_essential_genes_list:
									if k not in gene_list_false_positive:
										gene_list_false_positive.append(k)
						
		gene_list_true_negative = []
		gene_list_false_negative = []
		for rxn in md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2[state][flux_state].keys():
					g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
					g2r = str(g2r)
					i = g2r.split('or')
					for count, j in enumerate(i):
						j = j.split('and')
						for k in j:
							k = k.translate(None, ' ()')
							if not k:
								continue
							else:
								if k in YPEtoh_tested_genes_list:
									#Tested with YPD background experimentally. Therefore exclude YPD essential genes
									if k not in YPD_essential_genes_list:
										if k not in YPEtoh_essential_genes_list:
											if k not in gene_list_true_negative:
												gene_list_true_negative.append(k)
										if k in YPEtoh_essential_genes_list:
											if k not in gene_list_false_negative:
												gene_list_false_negative.append(k)

		idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN	
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	ethanol_positive_negative_dict_dataset_S2[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}



ethanol_positive_negative_dict_dataset_S2_eps = {}
for state in metabolicState_ethanol_list_dataset_S2_eps:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_eps[state][flux_state].keys():
			g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
			g2r = str(g2r)
			i = g2r.split('or')
			for count, j in enumerate(i):
				j = j.split('and')
				for k in j:
					k = k.translate(None, ' ()')
					if not k:
						continue
					else:
						if k in YPEtoh_tested_genes_list:
							#Tested with YPD background experimentally. Therefore exclude YPD essential genes
							if k not in YPD_essential_genes_list:
								if k in YPEtoh_essential_genes_list:
									if k not in gene_list_true_positive:
										gene_list_true_positive.append(k)
								if k not in YPEtoh_essential_genes_list:
									if k not in gene_list_false_positive:
										gene_list_false_positive.append(k)
						
		gene_list_true_negative = []
		gene_list_false_negative = []
		for rxn in md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_eps[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_eps[state][flux_state].keys():
					g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
					g2r = str(g2r)
					i = g2r.split('or')
					for count, j in enumerate(i):
						j = j.split('and')
						for k in j:
							k = k.translate(None, ' ()')
							if not k:
								continue
							else:
								if k in YPEtoh_tested_genes_list:
									#Tested with YPD background experimentally. Therefore exclude YPD essential genes
									if k not in YPD_essential_genes_list:
										if k not in YPEtoh_essential_genes_list:
											if k not in gene_list_true_negative:
												gene_list_true_negative.append(k)
										if k in YPEtoh_essential_genes_list:
											if k not in gene_list_false_negative:
												gene_list_false_negative.append(k)

		idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	ethanol_positive_negative_dict_dataset_S2_eps[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


ethanol_positive_negative_dict_dataset_S2_original_model = {}
for state in metabolicState_ethanol_list_dataset_S2_original_model:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_original_model[state][flux_state].keys():
			g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
			g2r = str(g2r)
			i = g2r.split('or')
			for count, j in enumerate(i):
				j = j.split('and')
				for k in j:
					k = k.translate(None, ' ()')
					if not k:
						continue
					else:
						if k in YPEtoh_tested_genes_list:
							#Tested with YPD background experimentally. Therefore exclude YPD essential genes
							if k not in YPD_essential_genes_list:
								if k in YPEtoh_essential_genes_list:
									if k not in gene_list_true_positive:
										gene_list_true_positive.append(k)
								if k not in YPEtoh_essential_genes_list:
									if k not in gene_list_false_positive:
										gene_list_false_positive.append(k)
						
		gene_list_true_negative = []
		gene_list_false_negative = []
		for rxn in md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_original_model[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_original_model[state][flux_state].keys():
					g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
					g2r = str(g2r)
					i = g2r.split('or')
					for count, j in enumerate(i):
						j = j.split('and')
						for k in j:
							k = k.translate(None, ' ()')
							if not k:
								continue
							else:
								if k in YPEtoh_tested_genes_list:
									#Tested with YPD background experimentally. Therefore exclude YPD essential genes
									if k not in YPD_essential_genes_list:
										if k not in YPEtoh_essential_genes_list:
											if k not in gene_list_true_negative:
												gene_list_true_negative.append(k)
										if k in YPEtoh_essential_genes_list:
											if k not in gene_list_false_negative:
												gene_list_false_negative.append(k)

		idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN	
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	ethanol_positive_negative_dict_dataset_S2_original_model[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}

ethanol_positive_negative_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_ethanol_list_dataset_S2_original_model_eps:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_original_model_eps[state][flux_state].keys():
			g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
			g2r = str(g2r)
			i = g2r.split('or')
			for count, j in enumerate(i):
				j = j.split('and')
				for k in j:
					k = k.translate(None, ' ()')
					if not k:
						continue
					else:
						if k in YPEtoh_tested_genes_list:
							#Tested with YPD background experimentally. Therefore exclude YPD essential genes
							if k not in YPD_essential_genes_list:
								if k in YPEtoh_essential_genes_list:
									if k not in gene_list_true_positive:
										gene_list_true_positive.append(k)
								if k not in YPEtoh_essential_genes_list:
									if k not in gene_list_false_positive:
										gene_list_false_positive.append(k)
						
		gene_list_true_negative = []
		gene_list_false_negative = []
		for rxn in md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_original_model_eps[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_original_model_eps[state][flux_state].keys():
					g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
					g2r = str(g2r)
					i = g2r.split('or')
					for count, j in enumerate(i):
						j = j.split('and')
						for k in j:
							k = k.translate(None, ' ()')
							if not k:
								continue
							else:
								if k in YPEtoh_tested_genes_list:
									#Tested with YPD background experimentally. Therefore exclude YPD essential genes
									if k not in YPD_essential_genes_list:
										if k not in YPEtoh_essential_genes_list:
											if k not in gene_list_true_negative:
												gene_list_true_negative.append(k)
										if k in YPEtoh_essential_genes_list:
											if k not in gene_list_false_negative:
												gene_list_false_negative.append(k)


		idRsflux_dict[flux_state] = {'true_positive': len(gene_list_true_positive), 'false_positive': len(gene_list_false_positive), 'true_negative': len(gene_list_true_negative), 'false_negative': len(gene_list_false_negative)}
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN	
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	ethanol_positive_negative_dict_dataset_S2_original_model_eps[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}

#Write out the results for the new EXAMO implementation
ethanol_positive_negative_dict_names = {'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm': 'E', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'E_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'E_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'E_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'E_lb', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'E_lb_g', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'E_lb_g_m_n_c', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'E_lb_m_n_c',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm': 'N', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'N_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'N_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'N_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh': 'N_lb', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g': 'N_lb_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c': 'N_lb_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c': 'N_lb_m_n_c',  '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'E_Ex', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'E_g_Ex', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'E_g_m_n_c_Ex', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'E_m_n_c_Ex', 
'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns': 'E_lb_Ex',
'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns': 'E_lb_g_Ex',
'151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'E_lb_g_m_n_c_Ex', '151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'E_lb_m_n_c_Ex',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'N_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'N_g_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_g_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_g_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_g_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_m_n_c_Ex'}


f = open('170503_ethanol_positive_negative.txt', 'w')

for i in ethanol_positive_negative_dict:
	print >> f, "Condition\tSensitivity Average\tSensitivity Std\tPPV Average\tPPV Std\tSpecificity Average\tSpecificity STD\tNPV Average\tNPV STD"
	print >> f, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (ethanol_positive_negative_dict_names[i], ethanol_positive_negative_dict[i]['sensitivity_average'], ethanol_positive_negative_dict[i]['sensitivity_std'], ethanol_positive_negative_dict[i]['PPV_average'], ethanol_positive_negative_dict[i]['PPV_std'], ethanol_positive_negative_dict[i]['specificity_average'], ethanol_positive_negative_dict[i]['specificity_std'], ethanol_positive_negative_dict[i]['npv_average'], ethanol_positive_negative_dict[i]['npv_std'])
f.close()

#determine the number of negative and positive genes present for the glucose data
glucose_positive_negative_dict = {}
for state in metabolicState_glucose_list:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate[state][flux_state].keys():
			g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
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
		for rxn in md_models[metabolicState_dict[state]].gene2rxn.keys():
			if idRsfluxstate[state][flux_state]:
				if rxn not in idRsfluxstate[state][flux_state].keys():
					g2r = md_models[metabolicState_dict[state]].gene2rxn[rxn]
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
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	glucose_positive_negative_dict[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


glucose_positive_negative_dict_dataset_S2 = {}
for state in metabolicState_glucose_list_dataset_S2:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2[state][flux_state].keys():
			g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
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
		for rxn in md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2[state][flux_state].keys():
					g2r = md_models_dataset_S2[metabolicState_dict_dataset_S2[state]].gene2rxn[rxn]
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
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	glucose_positive_negative_dict_dataset_S2[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


glucose_positive_negative_dict_dataset_S2_eps = {}
for state in metabolicState_glucose_list_dataset_S2_eps:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_eps[state][flux_state].keys():
			g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
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
		for rxn in md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_eps[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_eps[state][flux_state].keys():
					g2r = md_models_dataset_S2_eps[metabolicState_dict_dataset_S2_eps[state]].gene2rxn[rxn]
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
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	glucose_positive_negative_dict_dataset_S2_eps[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


glucose_positive_negative_dict_dataset_S2_original_model = {}
for state in metabolicState_glucose_list_dataset_S2_original_model:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_original_model[state][flux_state].keys():
			g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
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
		for rxn in md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_original_model[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_original_model[state][flux_state].keys():
					g2r = md_models_dataset_S2_original_model[metabolicState_dict_dataset_S2_original_model[state]].gene2rxn[rxn]
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
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	glucose_positive_negative_dict_dataset_S2_original_model[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


glucose_positive_negative_dict_dataset_S2_original_model_eps = {}
for state in metabolicState_glucose_list_dataset_S2_original_model_eps:
	average_genes = []
	idRsflux_dict = {}
	idRsflux_list = ['idRsflux1', 'idRsflux2', 'idRsflux3', 'idRsflux4', 'idRsflux5']
	for flux_state in idRsflux_list:
		gene_list_true_positive = []
		gene_list_false_positive = []
		for rxn in idRsfluxstate_dataset_S2_original_model_eps[state][flux_state].keys():
			g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
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
		for rxn in md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn.keys():
			if idRsfluxstate_dataset_S2_original_model_eps[state][flux_state]:
				if rxn not in idRsfluxstate_dataset_S2_original_model_eps[state][flux_state].keys():
					g2r = md_models_dataset_S2_original_model_eps[metabolicState_dict_dataset_S2_original_model_eps[state]].gene2rxn[rxn]
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
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try: 
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	sensitivity_list = [num1,num2,num3,num4,num5]
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_positive'])/(float(idRsflux_dict['idRsflux1']['true_positive'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_positive'])/(float(idRsflux_dict['idRsflux2']['true_positive'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_positive'])/(float(idRsflux_dict['idRsflux3']['true_positive'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_positive'])/(float(idRsflux_dict['idRsflux4']['true_positive'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_positive'])/(float(idRsflux_dict['idRsflux5']['true_positive'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	ppv_list = [num1,num2,num3,num4,num5]
	
	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_positive']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_positive']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_positive']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_positive']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_positive']))
	except ZeroDivisionError:
		num5 = np.NaN
	specificity_list = [num1,num2,num3,num4,num5]

	try:
		num1 = float(idRsflux_dict['idRsflux1']['true_negative'])/(float(idRsflux_dict['idRsflux1']['true_negative'])+float(idRsflux_dict['idRsflux1']['false_negative']))
	except ZeroDivisionError:
		num1 = np.NaN
	try:
		num2 = float(idRsflux_dict['idRsflux2']['true_negative'])/(float(idRsflux_dict['idRsflux2']['true_negative'])+float(idRsflux_dict['idRsflux2']['false_negative']))
	except ZeroDivisionError:
		num2 = np.NaN
	try:
		num3 = float(idRsflux_dict['idRsflux3']['true_negative'])/(float(idRsflux_dict['idRsflux3']['true_negative'])+float(idRsflux_dict['idRsflux3']['false_negative']))
	except ZeroDivisionError:
		num3 = np.NaN
	try:
		num4 = float(idRsflux_dict['idRsflux4']['true_negative'])/(float(idRsflux_dict['idRsflux4']['true_negative'])+float(idRsflux_dict['idRsflux4']['false_negative']))
	except ZeroDivisionError:
		num4 = np.NaN
	try:
		num5 = float(idRsflux_dict['idRsflux5']['true_negative'])/(float(idRsflux_dict['idRsflux5']['true_negative'])+float(idRsflux_dict['idRsflux5']['false_negative']))
	except ZeroDivisionError:
		num5 = np.NaN
	npv_list = [num1,num2,num3,num4,num5]

	if np.isnan(np.nanmean(sensitivity_list)):
		sensitivity_average = 0
	else:
		sensitivity_average = round(np.nanmean(sensitivity_list),4)
	if np.isnan(np.nanstd(sensitivity_list)):
		sensitivity_std = 0
	else:
		sensitivity_std = round(np.nanstd(sensitivity_list),4)
	if np.isnan(np.nanmean(ppv_list)):
		ppv_average = 0
	else:
		ppv_average = round(np.nanmean(ppv_list),4)
	if np.isnan(np.nanstd(ppv_list)):
		ppv_std = 0
	else:
		ppv_std = round(np.nanstd(ppv_list),4)
	if np.isnan(np.nanmean(specificity_list)):
		specificity_average = 0
	else:
		specificity_average = round(np.nanmean(specificity_list),4)
	if np.isnan(np.nanstd(specificity_list)):
		specificity_std = 0
	else:	
		specificity_std = round(np.nanstd(specificity_list),4)
	if np.isnan(np.nanmean(npv_list)):
		npv_average = 0
	else:
		npv_average = round(np.nanmean(npv_list),4)
	if np.isnan(np.nanstd(npv_list)):
		npv_std = 0
	else:
		npv_std = round(np.nanstd(npv_list),4)
	glucose_positive_negative_dict_dataset_S2_original_model_eps[state] = {'sensitivity_average': sensitivity_average, 'sensitivity_std': sensitivity_std, 'PPV_average': ppv_average, 'PPV_std': ppv_std, 'specificity_average': specificity_average, 'specificity_std': specificity_std, 'npv_average': npv_average, 'npv_std': npv_std}


#Write out the results for the new EXAMO implementation
glucose_positive_negative_dict_names = {'151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm': 'G', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'G_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'G_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'G_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'G_lb', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'G_lb_g', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'G_lb_g_m_n_c', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'G_lb_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm': 'N', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g': 'N_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c': 'N_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c': 'N_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD': 'N_lb', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g': 'N_lb_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c': 'N_lb_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c': 'N_lb_m_n_c',  '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'G_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'G_g_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'G_g_m_n_c_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'G_m_n_c_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns': 'G_lb_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns': 'G_lb_g_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'G_lb_g_m_n_c_Ex', '151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'g_lb_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns': 'N_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns': 'N_g_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_g_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_g_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_g_m_n_c_Ex', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns': 'N_lb_m_n_c_Ex'}

f = open('170503_glucose_positive_negative.txt', 'w')

for i in glucose_positive_negative_dict:
	print >> f, "Condition\tSensitivity Average\tSensitivity Std\tPPV Average\tPPV Std\tSpecificity Average\tSpecificity STD\tNPV Average\tNPV STD"
	print >> f, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (glucose_positive_negative_dict_names[i], glucose_positive_negative_dict[i]['sensitivity_average'], glucose_positive_negative_dict[i]['sensitivity_std'], glucose_positive_negative_dict[i]['PPV_average'], glucose_positive_negative_dict[i]['PPV_std'], glucose_positive_negative_dict[i]['specificity_average'], glucose_positive_negative_dict[i]['specificity_std'], glucose_positive_negative_dict[i]['npv_average'], glucose_positive_negative_dict[i]['npv_std'])

f.close()






n_groups = 20 
index = np.arange(n_groups)

sensitivity_ethanol = [ethanol_positive_negative_dict_dataset_S2_original_model['151012_ethanol_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2['151012_ethanol_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2_eps['151012_ethanol_0.25']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],
ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average']]

sensitivity_ethanol_negative = [ethanol_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['sensitivity_average'],ethanol_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],
ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['sensitivity_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average']]


sensitivity_glucose = [glucose_positive_negative_dict_dataset_S2_original_model['151012_glucose_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2['151012_glucose_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2_eps['151012_glucose_0.25']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],
glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average']]

sensitivity_glucose_negative = [glucose_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['sensitivity_average'],glucose_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],
glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['sensitivity_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_average']]


sensitivity_ethanol_std = [ethanol_positive_negative_dict_dataset_S2_original_model['151012_ethanol_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2['151012_ethanol_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2_eps['151012_ethanol_0.25']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],
ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std']]

sensitivity_ethanol_negative_std = [ethanol_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['sensitivity_std'],ethanol_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],
ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['sensitivity_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std']]


sensitivity_glucose_std = [glucose_positive_negative_dict_dataset_S2_original_model['151012_glucose_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2['151012_glucose_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2_eps['151012_glucose_0.25']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],
glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std']]

sensitivity_glucose_negative_std = [glucose_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['sensitivity_std'],glucose_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],
glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['sensitivity_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['sensitivity_std']]





PPV_ethanol = [ethanol_positive_negative_dict_dataset_S2_original_model['151012_ethanol_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2['151012_ethanol_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2_eps['151012_ethanol_0.25']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],
ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average']]

PPV_ethanol_negative = [ethanol_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['PPV_average'],ethanol_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],
ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['PPV_average'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average']]


PPV_glucose = [glucose_positive_negative_dict_dataset_S2_original_model['151012_glucose_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2['151012_glucose_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2_eps['151012_glucose_0.25']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],
glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['PPV_average'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average']]

PPV_glucose_negative = [glucose_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['PPV_average'],glucose_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],
glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['PPV_average'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_average']]


PPV_ethanol_std = [ethanol_positive_negative_dict_dataset_S2_original_model['151012_ethanol_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2['151012_ethanol_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2_eps['151012_ethanol_0.25']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],
ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std']]

PPV_ethanol_negative_std = [ethanol_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['PPV_std'],ethanol_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],
ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['PPV_std'], ethanol_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std']]


PPV_glucose_std = [glucose_positive_negative_dict_dataset_S2_original_model['151012_glucose_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2['151012_glucose_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2_eps['151012_glucose_0.25']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],
glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['PPV_std'], glucose_positive_negative_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std']]

PPV_glucose_negative_std = [glucose_positive_negative_dict_dataset_S2_original_model['151015_negative_control_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2['151015_negative_control_0.25']['PPV_std'],glucose_positive_negative_dict_dataset_S2_eps['151015_negative_control_0.25']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],
glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'],glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['PPV_std'], glucose_positive_negative_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['PPV_std']]



metabolicState_count_ethanol_0 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['count0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['count0'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['count0'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count0'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['count0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0']]

metabolicState_count_ethanol_1 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['count1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['count1'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['count1'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count1'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['count1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1']]

metabolicState_count_ethanol_2 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['count2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['count2'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['count2'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count2'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['count2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2']]

metabolicState_count_ethanol_3 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['count3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['count3'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['count3'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count3'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['count3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3']]

metabolicState_count_ethanol_4 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['count4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['count4'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['count4'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count4'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['count4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4']]

metabolicState_count_glucose_0 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['count0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['count0'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['count0'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count0'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0']]

metabolicState_count_glucose_1 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['count1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['count1'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['count1'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count1'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1']]

metabolicState_count_glucose_2 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['count2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['count2'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['count2'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count2'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2']]

metabolicState_count_glucose_3 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['count3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['count3'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['count3'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count3'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3']]

metabolicState_count_glucose_4 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['count4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['count4'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['count4'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count4'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4']]

metabolicState_count_Aerobic_0 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count0'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['count0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0']]

metabolicState_count_Aerobic_1 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count1'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['count1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1']]

metabolicState_count_Aerobic_2 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count2'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['count2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2']]

metabolicState_count_Aerobic_3 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count3'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['count3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3']]

metabolicState_count_Aerobic_4 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count4'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['count4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4']]

metabolicState_count_Anaerobic_0 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['count0'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count0'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['count0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0']]

metabolicState_count_Anaerobic_1 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['count1'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count1'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['count1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1']]

metabolicState_count_Anaerobic_2 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['count2'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count2'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['count2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2']]

metabolicState_count_Anaerobic_3 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['count3'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count3'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['count3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3']]

metabolicState_count_Anaerobic_4 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['count4'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count4'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['count4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4']]

metabolicState_count_negative_control_0 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['count0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['count0'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['count0'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count0'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count0']]

metabolicState_count_negative_control_1 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['count1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['count1'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['count1'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count1'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count1']]

metabolicState_count_negative_control_2 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['count2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['count2'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['count2'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count2'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count2']]

metabolicState_count_negative_control_3 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['count3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['count3'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['count3'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count3'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count3']]

metabolicState_count_negative_control_4 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['count4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['count4'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['count4'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['count4'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['count4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['count4']]





metabolicState_state_ethanol_0 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0']]

metabolicState_state_ethanol_1 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1']]

metabolicState_state_ethanol_2 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2']]

metabolicState_state_ethanol_3 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3']]

metabolicState_state_ethanol_4 = [mbaCandRxns_dict_dataset_S2_original_model['151012_ethanol_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_ethanol_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2['151012_ethanol_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_eps['151012_ethanol_0.25']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],
mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_ethanol_0.25_iMM904_NADcorrected_1127_FTHFLm_YPEtoh_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4']]

metabolicState_state_glucose_0 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0']]

metabolicState_state_glucose_1 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1']]

metabolicState_state_glucose_2 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2']]

metabolicState_state_glucose_3 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3']]

metabolicState_state_glucose_4 = [mbaCandRxns_dict_dataset_S2_original_model['151012_glucose_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151012_glucose_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2['151012_glucose_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_eps['151012_glucose_0.25']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],
mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151012_glucose_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4']]

metabolicState_state_Aerobic_0 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0']]

metabolicState_state_Aerobic_1 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1']]

metabolicState_state_Aerobic_2 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2']]

metabolicState_state_Aerobic_3 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3']]

metabolicState_state_Aerobic_4 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Aerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2['151006_Aerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_eps['151006_Aerobic_0.25']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],
mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4']]

metabolicState_state_Anaerobic_0 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0']]

metabolicState_state_Anaerobic_1 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1']]

metabolicState_state_Anaerobic_2 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2']]

metabolicState_state_Anaerobic_3 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3']]

metabolicState_state_Anaerobic_4 = [mbaCandRxns_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2['151006_Anaerobic_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_eps['151006_Anaerobic_0.25']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],
mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4']]

metabolicState_state_negative_control_0 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['metabolicState0'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState0'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState0']]

metabolicState_state_negative_control_1 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['metabolicState1'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState1'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState1']]

metabolicState_state_negative_control_2 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['metabolicState2'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState2'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState2']]

metabolicState_state_negative_control_3 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['metabolicState3'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState3'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState3']]

metabolicState_state_negative_control_4 = [mbaCandRxns_dict_dataset_S2_original_model['151015_negative_control_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2['151015_negative_control_0.25']['metabolicState4'],mbaCandRxns_dict_dataset_S2_eps['151015_negative_control_0.25']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],
mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'],mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c']['metabolicState4'], mbaCandRxns_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_YPD_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['metabolicState4']]



error_bar_length_ethanol_0 = np.repeat(0,len(metabolicState_count_ethanol_0))
error_bar_std_ethanol_0 = np.repeat(0,len(metabolicState_count_ethanol_0))
error_bar_length_glucose_0 = np.repeat(0,len(metabolicState_count_glucose_0))
error_bar_std_glucose_0 = np.repeat(0,len(metabolicState_count_glucose_0))
error_bar_length_Aerobic_0 = np.repeat(0,len(metabolicState_count_Aerobic_0))
error_bar_std_Aerobic_0 = np.repeat(0,len(metabolicState_count_Aerobic_0))
error_bar_length_Anaerobic_0 = np.repeat(0,len(metabolicState_count_Anaerobic_0))
error_bar_std_Anaerobic_0 = np.repeat(0,len(metabolicState_count_Anaerobic_0))
error_bar_length_negative_control_0 = np.repeat(0,len(metabolicState_count_negative_control_0))
error_bar_std_negative_control_0 = np.repeat(0,len(metabolicState_count_negative_control_0))

error_bar_length_ethanol_1 = np.repeat(0,len(metabolicState_count_ethanol_1))
error_bar_std_ethanol_1 = np.repeat(0,len(metabolicState_count_ethanol_1))
error_bar_length_glucose_1 = np.repeat(0,len(metabolicState_count_glucose_1))
error_bar_std_glucose_1 = np.repeat(0,len(metabolicState_count_glucose_1))
error_bar_length_Aerobic_1 = np.repeat(0,len(metabolicState_count_Aerobic_1))
error_bar_std_Aerobic_1 = np.repeat(0,len(metabolicState_count_Aerobic_1))
error_bar_length_Anaerobic_1 = np.repeat(0,len(metabolicState_count_Anaerobic_1))
error_bar_std_Anaerobic_1 = np.repeat(0,len(metabolicState_count_Anaerobic_1))
error_bar_length_negative_control_1 = np.repeat(0,len(metabolicState_count_negative_control_1))
error_bar_std_negative_control_1 = np.repeat(0,len(metabolicState_count_negative_control_1))

error_bar_length_ethanol_2 = np.repeat(0,len(metabolicState_count_ethanol_2))
error_bar_std_ethanol_2 = np.repeat(0,len(metabolicState_count_ethanol_2))
error_bar_length_glucose_2 = np.repeat(0,len(metabolicState_count_glucose_2))
error_bar_std_glucose_2 = np.repeat(0,len(metabolicState_count_glucose_2))
error_bar_length_Aerobic_2 = np.repeat(0,len(metabolicState_count_Aerobic_2))
error_bar_std_Aerobic_2 = np.repeat(0,len(metabolicState_count_Aerobic_2))
error_bar_length_Anaerobic_2 = np.repeat(0,len(metabolicState_count_Anaerobic_2))
error_bar_std_Anaerobic_2 = np.repeat(0,len(metabolicState_count_Anaerobic_2))
error_bar_length_negative_control_2 = np.repeat(0,len(metabolicState_count_negative_control_2))
error_bar_std_negative_control_2 = np.repeat(0,len(metabolicState_count_negative_control_2))

error_bar_length_ethanol_3 = np.repeat(0,len(metabolicState_count_ethanol_3))
error_bar_std_ethanol_3 = np.repeat(0,len(metabolicState_count_ethanol_3))
error_bar_length_glucose_3 = np.repeat(0,len(metabolicState_count_glucose_3))
error_bar_std_glucose_3 = np.repeat(0,len(metabolicState_count_glucose_3))
error_bar_length_Aerobic_3 = np.repeat(0,len(metabolicState_count_Aerobic_3))
error_bar_std_Aerobic_3 = np.repeat(0,len(metabolicState_count_Aerobic_3))
error_bar_length_Anaerobic_3 = np.repeat(0,len(metabolicState_count_Anaerobic_3))
error_bar_std_Anaerobic_3 = np.repeat(0,len(metabolicState_count_Anaerobic_3))
error_bar_length_negative_control_3 = np.repeat(0,len(metabolicState_count_negative_control_3))
error_bar_std_negative_control_3 = np.repeat(0,len(metabolicState_count_negative_control_3))

error_bar_length_ethanol_4 = np.repeat(0,len(metabolicState_count_ethanol_4))
error_bar_std_ethanol_4 = np.repeat(0,len(metabolicState_count_ethanol_4))
error_bar_length_glucose_4 = np.repeat(0,len(metabolicState_count_glucose_4))
error_bar_std_glucose_4 = np.repeat(0,len(metabolicState_count_glucose_4))
error_bar_length_Aerobic_4 = np.repeat(0,len(metabolicState_count_Aerobic_4))
error_bar_std_Aerobic_4 = np.repeat(0,len(metabolicState_count_Aerobic_4))
error_bar_length_Anaerobic_4 = np.repeat(0,len(metabolicState_count_Anaerobic_4))
error_bar_std_Anaerobic_4 = np.repeat(0,len(metabolicState_count_Anaerobic_4))
error_bar_length_negative_control_4 = np.repeat(0,len(metabolicState_count_negative_control_4))
error_bar_std_negative_control_4 = np.repeat(0,len(metabolicState_count_negative_control_4))



fig, (ax0,ax1,ax2) = plt.subplots(3, sharex=True)

plt.close()
bar_width = 0.17
error_bar_length = np.repeat(0,len(sensitivity_ethanol))
rects1 = ax0.bar(index, sensitivity_ethanol, bar_width, alpha = 0.8, yerr = [error_bar_length, sensitivity_ethanol_std], capsize = 3, color = 'b', label = 'Ethanol')
rects2 = ax0.bar(index + bar_width, sensitivity_ethanol_negative, bar_width, alpha = 0.8, yerr = [error_bar_length, sensitivity_ethanol_negative_std], capsize = 3, color = 'r', label = 'Ethanol Neg') 
rects3 = ax0.bar(index + bar_width + bar_width, sensitivity_glucose, bar_width, alpha = 0.8, yerr = [error_bar_length,sensitivity_glucose_std], capsize = 3, color = 'g', label = 'Glucose') 
rects4  = ax0.bar(index + bar_width + bar_width + bar_width, sensitivity_glucose_negative, bar_width, alpha = 0.8, yerr = [error_bar_length, sensitivity_glucose_std], capsize = 3, color = 'm', label = 'Glucose Neg')
#plt.xlabel('Model Parameters',fontweight='bold',fontsize=10) 
ax0.set_ylabel('Sensitivity (TP/(TP+FN))',fontweight='bold',fontsize=9)
ax0.tick_params(labelsize=9)
#plt.title('Sensitivity of Conditions',fontweight='bold') 
#ax[0].xticks(index + 2*bar_width, ('C_orig', 'C_orig_eps', 'C_mod', 'C_mod_eps', 'C', 'C_Ex', 'C_g', 'C_g_Ex', 'C_lb', 'C_lb_EX', 'C_lb_g', 'C_lb_g_EX', 'C_m_n_c', 'C_m_n_c_EX', 'C_m_n_c_g', 'C_m_n_c_g_EX', 'C_m_n_c_lb', 'C_m_n_c_lb_EX', 'C_m_n_c_lb_g', 'C_m_n_c_lb_g_EX'), rotation=90)
a = ax0.legend(loc="lower right",fontsize=8)
a.get_frame().set_alpha(1)
a.set_title('Condition',prop={'weight':'bold','size':'9'})
#plt.tight_layout()
b = ax0.axvspan(-0.25,1.75,alpha=0.1,color='blue')
c = ax0.axvspan(1.75,3.75,alpha=0.1,color='orange')
d = ax0.axvspan(3.75,19.75,alpha=0.1,color='red')
e = ax0.legend([b,c,d],['EXAMO','EXAMO with Model\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],loc='lower center',fontsize=8)
e.get_frame().set_alpha(1)
e.set_title('Software',prop={'weight':'bold','size':'9'})
ax0.add_artist(a)
ax0.add_artist(e)
ax0.set_ylim(0,1)
#fig = plt.gcf()
#fig.set_size_inches(6.69,4.5)
#fig.savefig('Figure_2.png',dpi=300)


bar_width = 0.17
error_bar_length = np.repeat(0,len(PPV_ethanol))
rects1 = ax1.bar(index, PPV_ethanol, bar_width, alpha = 0.8, yerr = [error_bar_length, PPV_ethanol_std], capsize = 3, color = 'b', label = 'Ethanol')
rects2 = ax1.bar(index + bar_width, PPV_ethanol_negative, bar_width, alpha = 0.8, yerr = [error_bar_length, PPV_ethanol_negative_std], capsize = 3, color = 'r', label = 'Ethanol Neg') 
rects3 = ax1.bar(index + bar_width + bar_width, PPV_glucose, bar_width, alpha = 0.8, yerr = [error_bar_length,PPV_glucose_std], capsize = 3, color = 'g', label = 'Glucose') 
rects4  = ax1.bar(index + bar_width + bar_width + bar_width, PPV_glucose_negative, bar_width, alpha = 0.8, yerr = [error_bar_length, PPV_glucose_std], capsize = 3, color = 'm', label = 'Glucose Neg')
#plt.xlabel('Model Parameters',fontweight='bold',fontsize=10) 
ax1.set_ylabel('Precision (TP/(TP+FP))',fontweight='bold',fontsize=10)
ax1.tick_params(labelsize=9)
#plt.title('Precision of Conditions',fontweight='bold') 
#plt.xticks(index + 2*bar_width, ('C_orig', 'C_orig_eps', 'C_mod', 'C_mod_eps', 'C', 'C_Ex', 'C_g', 'C_g_Ex', 'C_lb', 'C_lb_EX', 'C_lb_g', 'C_lb_g_EX', 'C_m_n_c', 'C_m_n_c_EX', 'C_m_n_c_g', 'C_m_n_c_g_EX', 'C_m_n_c_lb', 'C_m_n_c_lb_EX', 'C_m_n_c_lb_g', 'C_m_n_c_lb_g_EX'), rotation=90) 
a = ax1.legend(loc="upper right",fontsize=8)
a.get_frame().set_alpha(1)
a.set_title('Condition',prop={'weight':'bold','size':'9'})
#plt.tight_layout()
b = ax1.axvspan(-0.25,1.75,alpha=0.1,color='blue')
c = ax1.axvspan(1.75,3.75,alpha=0.1,color='orange')
d = ax1.axvspan(3.75,19.75,alpha=0.1,color='red')
e = ax1.legend([b,c,d],['EXAMO','EXAMO with Model\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],loc='upper center',fontsize=8)
e.get_frame().set_alpha(1)
e.set_title('Software',prop={'weight':'bold','size':'9'})
ax1.add_artist(a)
ax1.add_artist(e)
ax1.set_ylim(0,0.3)
#fig = plt.gcf()
#fig.set_size_inches(6.69,4.5)
#fig.savefig('Figure_3.png',dpi=300)


bar_width = 0.136
#This top section is just for indexing so that the legend does not have hatching
rects1 = plt.bar(index, metabolicState_count_ethanol_0, bar_width, yerr = [error_bar_length_ethanol_0, error_bar_std_ethanol_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol')
rects2 = plt.bar(index+bar_width, metabolicState_count_glucose_0, bar_width, yerr = [error_bar_length_glucose_0, error_bar_std_glucose_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose')
rects3 = plt.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_0, bar_width, yerr = [error_bar_length_Aerobic_0, error_bar_std_Aerobic_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic')
rects4 = plt.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_0, bar_width, yerr = [error_bar_length_Anaerobic_0, error_bar_std_Anaerobic_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic')
rects5 = plt.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_0, bar_width, yerr = [error_bar_length_negative_control_0, error_bar_std_negative_control_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control')

#This next section is for creating a legend with hatching
hatch1 = plt.bar(index, metabolicState_count_ethanol_0, bar_width, yerr = [error_bar_length_ethanol_0, error_bar_std_ethanol_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'white', label = 'Solvable', edgecolor = 'black')
hatch2 = plt.bar(index+bar_width, metabolicState_count_glucose_0, bar_width, yerr = [error_bar_length_glucose_0, error_bar_std_glucose_0], error_kw=dict(capsize = 2, capthick = 0.25), hatch = "//////", alpha = 0.8, color = 'white', label = 'Unsolvable', edgecolor = 'black')

#Now create the plots with bars
#fig2 = ax2.figure()
#ax2_fig = ax2.add_subplot(111)
rects = ax2.bar(index, metabolicState_count_ethanol_0, bar_width, yerr = [error_bar_length_ethanol_0, error_bar_std_ethanol_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol') + \
    ax2.bar(index+bar_width, metabolicState_count_glucose_0, bar_width, yerr = [error_bar_length_glucose_0, error_bar_std_glucose_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose') + \
    ax2.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_0, bar_width, yerr = [error_bar_length_Aerobic_0, error_bar_std_Aerobic_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_0, bar_width, yerr = [error_bar_length_Anaerobic_0, error_bar_std_Anaerobic_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_0, bar_width, yerr = [error_bar_length_negative_control_0, error_bar_std_negative_control_0], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control') + \
    ax2.bar(index, metabolicState_count_ethanol_1, bar_width, bottom = metabolicState_count_ethanol_0, yerr = [error_bar_length_ethanol_1, error_bar_std_ethanol_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol') + \
    ax2.bar(index+bar_width, metabolicState_count_glucose_1, bar_width, bottom = metabolicState_count_glucose_0, yerr = [error_bar_length_glucose_1, error_bar_std_glucose_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose') + \
    ax2.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_1, bar_width, bottom = metabolicState_count_Aerobic_0, yerr = [error_bar_length_Aerobic_1, error_bar_std_Aerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_1, bar_width, bottom = metabolicState_count_Anaerobic_0, yerr = [error_bar_length_Anaerobic_1, error_bar_std_Anaerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_1, bar_width, bottom = metabolicState_count_negative_control_0, yerr = [error_bar_length_negative_control_1, error_bar_std_negative_control_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control') + \
    ax2.bar(index, metabolicState_count_ethanol_2, bar_width, bottom = np.sum((metabolicState_count_ethanol_0, metabolicState_count_ethanol_1),axis=0), yerr = [error_bar_length_ethanol_1, error_bar_std_ethanol_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol') + \
    ax2.bar(index+bar_width, metabolicState_count_glucose_2, bar_width, bottom = np.sum((metabolicState_count_glucose_0, metabolicState_count_glucose_1),axis=0), yerr = [error_bar_length_glucose_1, error_bar_std_glucose_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose') + \
    ax2.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_2, bar_width, bottom = np.sum((metabolicState_count_Aerobic_0, metabolicState_count_Aerobic_1),axis=0), yerr = [error_bar_length_Aerobic_1, error_bar_std_Aerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_2, bar_width, bottom = np.sum((metabolicState_count_Anaerobic_0, metabolicState_count_Anaerobic_1),axis=0), yerr = [error_bar_length_Anaerobic_1, error_bar_std_Anaerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_2, bar_width, bottom = np.sum((metabolicState_count_negative_control_0, metabolicState_count_negative_control_1),axis=0), yerr = [error_bar_length_negative_control_1, error_bar_std_negative_control_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control') + \
    ax2.bar(index, metabolicState_count_ethanol_3, bar_width, bottom = np.sum((metabolicState_count_ethanol_0, metabolicState_count_ethanol_1, metabolicState_count_ethanol_2),axis=0), yerr = [error_bar_length_ethanol_1, error_bar_std_ethanol_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol') + \
    ax2.bar(index+bar_width, metabolicState_count_glucose_3, bar_width, bottom = np.sum((metabolicState_count_glucose_0, metabolicState_count_glucose_1, metabolicState_count_glucose_2),axis=0), yerr = [error_bar_length_glucose_1, error_bar_std_glucose_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose') + \
    ax2.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_3, bar_width, bottom = np.sum((metabolicState_count_Aerobic_0, metabolicState_count_Aerobic_1, metabolicState_count_Aerobic_2),axis=0), yerr = [error_bar_length_Aerobic_1, error_bar_std_Aerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_3, bar_width, bottom = np.sum((metabolicState_count_Anaerobic_0, metabolicState_count_Anaerobic_1, metabolicState_count_Anaerobic_2),axis=0), yerr = [error_bar_length_Anaerobic_1, error_bar_std_Anaerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_3, bar_width, bottom = np.sum((metabolicState_count_negative_control_0, metabolicState_count_negative_control_1, metabolicState_count_negative_control_2),axis=0), yerr = [error_bar_length_negative_control_1, error_bar_std_negative_control_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control') + \
    ax2.bar(index, metabolicState_count_ethanol_4, bar_width, bottom = np.sum((metabolicState_count_ethanol_0, metabolicState_count_ethanol_1, metabolicState_count_ethanol_2, metabolicState_count_ethanol_3),axis=0), yerr = [error_bar_length_ethanol_1, error_bar_std_ethanol_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'b', label = 'Ethanol') + \
    ax2.bar(index+bar_width, metabolicState_count_glucose_4, bar_width, bottom = np.sum((metabolicState_count_glucose_0, metabolicState_count_glucose_1, metabolicState_count_glucose_2, metabolicState_count_glucose_3),axis=0), yerr = [error_bar_length_glucose_1, error_bar_std_glucose_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'g', label = 'Glucose') + \
    ax2.bar(index+bar_width+bar_width, metabolicState_count_Aerobic_4, bar_width, bottom = np.sum((metabolicState_count_Aerobic_0, metabolicState_count_Aerobic_1, metabolicState_count_Aerobic_2, metabolicState_count_Aerobic_3),axis=0), yerr = [error_bar_length_Aerobic_1, error_bar_std_Aerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'r', label = 'Aerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width, metabolicState_count_Anaerobic_4, bar_width, bottom = np.sum((metabolicState_count_Anaerobic_0, metabolicState_count_Anaerobic_1, metabolicState_count_Anaerobic_2, metabolicState_count_Anaerobic_3),axis=0), yerr = [error_bar_length_Anaerobic_1, error_bar_std_Anaerobic_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'purple', label = 'Anaerobic') + \
    ax2.bar(index+bar_width+bar_width+bar_width+bar_width, metabolicState_count_negative_control_4, bar_width, bottom = np.sum((metabolicState_count_negative_control_0, metabolicState_count_negative_control_1, metabolicState_count_negative_control_2, metabolicState_count_negative_control_3),axis=0), yerr = [error_bar_length_negative_control_1, error_bar_std_negative_control_1], error_kw=dict(capsize = 2, capthick = 0.25), alpha = 0.8, color = 'gray', label = 'Negative Control')
   
metabolicState_state_list = [metabolicState_state_ethanol_0,metabolicState_state_glucose_0,metabolicState_state_Aerobic_0,metabolicState_state_Anaerobic_0,metabolicState_state_negative_control_0,
                      metabolicState_state_ethanol_1,metabolicState_state_glucose_1,metabolicState_state_Aerobic_1,metabolicState_state_Anaerobic_1,metabolicState_state_negative_control_1,
                      metabolicState_state_ethanol_2,metabolicState_state_glucose_2,metabolicState_state_Aerobic_2,metabolicState_state_Anaerobic_2,metabolicState_state_negative_control_2,
                      metabolicState_state_ethanol_3,metabolicState_state_glucose_3,metabolicState_state_Aerobic_3,metabolicState_state_Anaerobic_3,metabolicState_state_negative_control_3,
                      metabolicState_state_ethanol_4,metabolicState_state_glucose_4,metabolicState_state_Aerobic_4,metabolicState_state_Anaerobic_4,metabolicState_state_negative_control_4]
metablicState_list_flattened = sum(metabolicState_state_list, [])
patterns = []
for state in metablicState_list_flattened:
    if state == 0:
        patterns.append("//////")
    if state == 1:
        patterns.append("")
    
for bar, pattern in zip(rects, patterns):
        bar.set_hatch(pattern)

ax2.set_xlabel('Model Parameters',fontweight='bold',fontsize=10) 
ax2.set_ylabel('Number of Successfully\nPruned Models',fontweight='bold',fontsize=10) 
ax2.tick_params(labelsize=9)
#plt.title('Comparison of Successfully Pruned Models and Optimizations',fontweight='bold') 
ax2.set_xticks(index + 2*bar_width)
ax2.set_xticklabels(('C_orig', 'C_orig_eps', 'C_mod', 'C_mod_eps', 'C', 'C_Ex', 'C_g', 'C_g_Ex', 'C_lb', 'C_lb_EX', 'C_lb_g', 'C_lb_g_EX', 'C_m_n_c', 'C_m_n_c_EX', 'C_m_n_c_g', 'C_m_n_c_g_EX', 'C_m_n_c_lb', 'C_m_n_c_lb_EX', 'C_m_n_c_lb_g', 'C_m_n_c_lb_g_EX'), rotation=90)
a = ax2.legend([rects1,rects2,rects3,rects4,rects5],['Ethanol','Glucose','Aerobic','Anaerobic','Negative Control'],loc="upper right",fontsize=8)
a.get_frame().set_alpha(1)
a.set_title('Condition',prop={'weight':'bold','size':'9'})
b = ax2.axvspan(-0.25,1.75,alpha=0.1,color='blue')
c = ax2.axvspan(1.75,3.75,alpha=0.1,color='orange')
d = ax2.axvspan(3.75,19.75,alpha=0.1,color='red')
e = ax2.legend([b,c,d],['EXAMO','EXAMO with Model\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],loc='upper center',fontsize=8)
e.get_frame().set_alpha(1)
e.set_title('Software',prop={'weight':'bold','size':'9'})
f = ax2.legend([hatch1,hatch2],['Solvable','Unsolvable'],loc='upper left',fontsize=8)
f.get_frame().set_alpha(1)
f.set_title('Solvability',prop={'weight':'bold','size':'9'})
ax2.add_artist(a)
ax2.add_artist(e)
ax2.set_ylim(0,540)

ax0.text(-3.7,1,'a',fontsize=22)
ax1.text(-3.7,0.3,'b',fontsize=22)
ax2.text(-3.7,540,'c',fontsize=22)


plt.tight_layout()
fig.subplots_adjust(left=0.11,bottom=0.155,top=0.98,right=0.99)
fig.set_size_inches(6.69,8.85)
plt.show()
fig.savefig('Figure_2.png',dpi=300)





#import list of fluxes
Rintala_Fluxes_file = open('Rintala_Fluxes.csv', 'r')
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
		Aerobic_Average = (Rep1lb_Aerobic+Rep1ub_Aerobic+Rep2lb_Aerobic+Rep2ub_Aerobic)/4
		Anaerobic_Average = (Rep1lb_Anaerobic+Rep1ub_Anaerobic+Rep2lb_Anaerobic+Rep2ub_Anaerobic)/4
		experimental_fluxes_rxns[rxnid] = {'rxns': rxn_list, 'Aerobic_Average': Aerobic_Average, 'Anaerobic_Average': Anaerobic_Average}

metabolicState_Aerobic_list = ['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns',
'151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']


metabolicState_Aerobic_list_dataset_S2 = ['151006_Aerobic_0.25', '151015_negative_control_0.25']

metabolicState_Aerobic_list_dataset_S2_eps = ['151006_Aerobic_0.25', '151015_negative_control_0.25']

metabolicState_Aerobic_list_dataset_S2_original_model = ['151006_Aerobic_0.25', '151015_negative_control_0.25']

metabolicState_Aerobic_list_dataset_S2_original_model_eps = ['151006_Aerobic_0.25', '151015_negative_control_0.25']


metabolicState_Anaerobic_list = [ '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c',  '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns',  '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns', '151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']


metabolicState_Anaerobic_list_dataset_S2 = ['151006_Anaerobic_0.25', '151015_negative_control_0.25']

metabolicState_Anaerobic_list_dataset_S2_eps = ['151006_Anaerobic_0.25', '151015_negative_control_0.25']

metabolicState_Anaerobic_list_dataset_S2_original_model = ['151006_Anaerobic_0.25', '151015_negative_control_0.25']

metabolicState_Anaerobic_list_dataset_S2_original_model_eps = ['151006_Anaerobic_0.25', '151015_negative_control_0.25']


rxn_Aerobic_metabolicState_dict = {}
rxn_Aerobic_tot_difference = {}
rxn_Aerobic_tot_difference_per = {}
for state in metabolicState_Aerobic_list:
	rxn_Aerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]			
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]
				rxn_average = fluxavgdict[state][rxn]
				rxn_var = fluxvardict[state][rxn]
				rxn_average_0 = fluxdict[state][rxn][0]
				rxn_average_1 = fluxdict[state][rxn][1]
				rxn_average_2 = fluxdict[state][rxn][2]
				rxn_average_3 = fluxdict[state][rxn][3]
				rxn_average_4 = fluxdict[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Aerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Aerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference, 'rxnid': rxn_from_model}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Aerobic_metabolicState_dict[state] = rxn_Aerobic_dict
	rxn_Aerobic_tot_difference[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Aerobic_tot_difference_per[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Aerobic_metabolicState_dict_dataset_S2 = {}
rxn_Aerobic_tot_difference_dataset_S2 = {}
rxn_Aerobic_tot_difference_per_dataset_S2 = {}
for state in metabolicState_Aerobic_list_dataset_S2:
	rxn_Aerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]			
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]
				rxn_average = fluxavgdict_dataset_S2[state][rxn]
				rxn_var = fluxvardict_dataset_S2[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Aerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Aerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference, 'rxnid': rxn_from_model}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Aerobic_metabolicState_dict_dataset_S2[state] = rxn_Aerobic_dict
	rxn_Aerobic_tot_difference_dataset_S2[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Aerobic_tot_difference_per_dataset_S2[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Aerobic_metabolicState_dict_dataset_S2_eps = {}
rxn_Aerobic_tot_difference_dataset_S2_eps = {}
rxn_Aerobic_tot_difference_per_dataset_S2_eps = {}
for state in metabolicState_Aerobic_list_dataset_S2_eps:
	rxn_Aerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]			
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_eps[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_eps[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_eps[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_eps[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_eps[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_eps[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_eps[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]
				rxn_average = fluxavgdict_dataset_S2_eps[state][rxn]
				rxn_var = fluxvardict_dataset_S2_eps[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_eps[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_eps[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_eps[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_eps[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_eps[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Aerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Aerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference, 'rxnid': rxn_from_model}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Aerobic_metabolicState_dict_dataset_S2_eps[state] = rxn_Aerobic_dict
	rxn_Aerobic_tot_difference_dataset_S2_eps[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Aerobic_tot_difference_per_dataset_S2_eps[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Aerobic_metabolicState_dict_dataset_S2_original_model = {}
rxn_Aerobic_tot_difference_dataset_S2_original_model = {}
rxn_Aerobic_tot_difference_per_dataset_S2_original_model = {}
for state in metabolicState_Aerobic_list_dataset_S2_original_model:
	rxn_Aerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]			
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_original_model[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_original_model[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_original_model[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_original_model[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_original_model[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_original_model[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_original_model[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]
				rxn_average = fluxavgdict_dataset_S2_original_model[state][rxn]
				rxn_var = fluxvardict_dataset_S2_original_model[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_original_model[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_original_model[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_original_model[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_original_model[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_original_model[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Aerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Aerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference, 'rxnid': rxn_from_model}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Aerobic_metabolicState_dict_dataset_S2_original_model[state] = rxn_Aerobic_dict
	rxn_Aerobic_tot_difference_dataset_S2_original_model[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Aerobic_tot_difference_per_dataset_S2_original_model[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps = {}
rxn_Aerobic_tot_difference_dataset_S2_original_model_eps = {}
rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps = {}
for state in metabolicState_Aerobic_list_dataset_S2_original_model_eps:
	rxn_Aerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]			
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_original_model_eps[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_original_model_eps[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_original_model_eps[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_original_model_eps[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_original_model_eps[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_original_model_eps[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_original_model_eps[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_from_model = experimental_fluxes_rxns[rxnid]['rxns'][0]
				rxn_average = fluxavgdict_dataset_S2_original_model_eps[state][rxn]
				rxn_var = fluxvardict_dataset_S2_original_model_eps[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_original_model_eps[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_original_model_eps[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_original_model_eps[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_original_model_eps[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_original_model_eps[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Aerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Aerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Aerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Aerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference, 'rxnid': rxn_from_model}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps[state] = rxn_Aerobic_dict
	rxn_Aerobic_tot_difference_dataset_S2_original_model_eps[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]



rxn_Anaerobic_metabolicState_dict = {}
rxn_Anaerobic_tot_difference = {}
rxn_Anaerobic_tot_difference_per = {}
for state in metabolicState_Anaerobic_list:
	rxn_Anaerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_average = fluxavgdict[state][rxn]
				rxn_var = fluxvardict[state][rxn]
				rxn_average_0 = fluxdict[state][rxn][0]
				rxn_average_1 = fluxdict[state][rxn][1]
				rxn_average_2 = fluxdict[state][rxn][2]
				rxn_average_3 = fluxdict[state][rxn][3]
				rxn_average_4 = fluxdict[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Anaerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Anaerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Anaerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Anaerobic_metabolicState_dict[state] = rxn_Anaerobic_dict
	rxn_Anaerobic_tot_difference[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Anaerobic_tot_difference_per[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Anaerobic_metabolicState_dict_dataset_S2 = {}
rxn_Anaerobic_tot_difference_dataset_S2 = {}
rxn_Anaerobic_tot_difference_per_dataset_S2 = {}
for state in metabolicState_Anaerobic_list_dataset_S2:
	rxn_Anaerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_average = fluxavgdict_dataset_S2[state][rxn]
				rxn_var = fluxvardict_dataset_S2[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Anaerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Anaerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Anaerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Anaerobic_metabolicState_dict_dataset_S2[state] = rxn_Anaerobic_dict
	rxn_Anaerobic_tot_difference_dataset_S2[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Anaerobic_tot_difference_per_dataset_S2[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Anaerobic_metabolicState_dict_dataset_S2_eps = {}
rxn_Anaerobic_tot_difference_dataset_S2_eps = {}
rxn_Anaerobic_tot_difference_per_dataset_S2_eps = {}
for state in metabolicState_Anaerobic_list_dataset_S2_eps:
	rxn_Anaerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_eps[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_eps[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_eps[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_eps[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_eps[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_eps[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_eps[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_average = fluxavgdict_dataset_S2_eps[state][rxn]
				rxn_var = fluxvardict_dataset_S2_eps[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_eps[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_eps[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_eps[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_eps[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_eps[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Anaerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Anaerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Anaerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Anaerobic_metabolicState_dict_dataset_S2_eps[state] = rxn_Anaerobic_dict
	rxn_Anaerobic_tot_difference_dataset_S2_eps[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Anaerobic_tot_difference_per_dataset_S2_eps[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model = {}
rxn_Anaerobic_tot_difference_dataset_S2_original_model = {}
rxn_Anaerobic_tot_difference_per_dataset_S2_original_model = {}
for state in metabolicState_Anaerobic_list_dataset_S2_original_model:
	rxn_Anaerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_original_model[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_original_model[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_original_model[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_original_model[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_original_model[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_original_model[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_original_model[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_average = fluxavgdict_dataset_S2_original_model[state][rxn]
				rxn_var = fluxvardict_dataset_S2_original_model[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_original_model[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_original_model[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_original_model[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_original_model[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_original_model[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Anaerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Anaerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Anaerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model[state] = rxn_Anaerobic_dict
	rxn_Anaerobic_tot_difference_dataset_S2_original_model[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Anaerobic_tot_difference_per_dataset_S2_original_model[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]


rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps = {}
rxn_Anaerobic_tot_difference_dataset_S2_original_model_eps = {}
rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps = {}
for state in metabolicState_Anaerobic_list_dataset_S2_original_model_eps:
	rxn_Anaerobic_dict = {}
	tot_abs_difference = 0
	tot_abs_difference_0 = 0
	tot_abs_difference_1 = 0
	tot_abs_difference_2 = 0
	tot_abs_difference_3 = 0
	tot_abs_difference_4 = 0
	experimental_total = 0
	for rxnid in experimental_fluxes_rxns.keys():
		if len(experimental_fluxes_rxns[rxnid]['rxns']) > 1:
			rxn_sum = 0
			rxn_sum_0 = 0
			rxn_sum_1 = 0
			rxn_sum_2 = 0
			rxn_sum_3 = 0
			rxn_sum_4 = 0
			rxn_var_sum = 0
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_sum = rxn_sum + fluxavgdict_dataset_S2_original_model_eps[state][rxn]
				rxn_sum_0 = rxn_sum_0 + fluxdict_dataset_S2_original_model_eps[state][rxn][0]
				rxn_sum_1 = rxn_sum_1 + fluxdict_dataset_S2_original_model_eps[state][rxn][1]
				rxn_sum_2 = rxn_sum_2 + fluxdict_dataset_S2_original_model_eps[state][rxn][2]
				rxn_sum_3 = rxn_sum_3 + fluxdict_dataset_S2_original_model_eps[state][rxn][3]
				rxn_sum_4 = rxn_sum_4 + fluxdict_dataset_S2_original_model_eps[state][rxn][4]
				rxn_var_sum = rxn_var_sum + (fluxvardict_dataset_S2_original_model_eps[state][rxn]**2)
			rxn_var = math.sqrt(rxn_var_sum)
			rxn_average = rxn_sum
			rxn_average_0 = rxn_sum_0
			rxn_average_1 = rxn_sum_1
			rxn_average_2 = rxn_sum_2
			rxn_average_3 = rxn_sum_3
			rxn_average_4 = rxn_sum_4
		else:
			for rxn in experimental_fluxes_rxns[rxnid]['rxns']:
				rxn_average = fluxavgdict_dataset_S2_original_model_eps[state][rxn]
				rxn_var = fluxvardict_dataset_S2_original_model_eps[state][rxn]
				rxn_average_0 = fluxdict_dataset_S2_original_model_eps[state][rxn][0]
				rxn_average_1 = fluxdict_dataset_S2_original_model_eps[state][rxn][1]
				rxn_average_2 = fluxdict_dataset_S2_original_model_eps[state][rxn][2]
				rxn_average_3 = fluxdict_dataset_S2_original_model_eps[state][rxn][3]
				rxn_average_4 = fluxdict_dataset_S2_original_model_eps[state][rxn][4]
		experimental_total = experimental_total + experimental_fluxes_rxns[rxnid]['Anaerobic_Average']
		abs_difference = abs(abs(rxn_average) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_0 = abs(abs(rxn_average_0) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_1 = abs(abs(rxn_average_1) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_2 = abs(abs(rxn_average_2) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_3 = abs(abs(rxn_average_3) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		abs_difference_4 = abs(abs(rxn_average_4) - abs(experimental_fluxes_rxns[rxnid]['Anaerobic_Average']))
		tot_abs_difference = tot_abs_difference + abs_difference
		tot_abs_difference_0 = tot_abs_difference_0 + abs_difference_0
		tot_abs_difference_1 = tot_abs_difference_1 + abs_difference_1
		tot_abs_difference_2 = tot_abs_difference_2 + abs_difference_2
		tot_abs_difference_3 = tot_abs_difference_3 + abs_difference_3
		tot_abs_difference_4 = tot_abs_difference_4 + abs_difference_4
		rxn_Anaerobic_dict[rxnid] = {'Experimental_Flux': experimental_fluxes_rxns[rxnid]['Anaerobic_Average'], 'Modeled_Flux': rxn_average, 'abs_difference': abs_difference}
	total_experimental_percentage_difference = round(tot_abs_difference/experimental_total,4)
	rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps[state] = rxn_Anaerobic_dict
	rxn_Anaerobic_tot_difference_dataset_S2_original_model_eps[state] = {'tot_abs_difference': tot_abs_difference}
	rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps[state] = [tot_abs_difference_0, tot_abs_difference_1, tot_abs_difference_2, tot_abs_difference_3, tot_abs_difference_4]





from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

C_all_Aerobic_df = pd.DataFrame()
C_all_Aerobic_df['tot_abs_difference'] =  [rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4]]
C_all_Aerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

C_all_Aerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_all_Aerobic = ols(design, C_all_Aerobic_df).fit()
aov_table_all_Aerobic = anova_lm(model_all_Aerobic, typ=2)
print aov_table_all_Aerobic


C_m_n_c_Aerobic_df = pd.DataFrame()
C_m_n_c_Aerobic_df['tot_abs_difference'] =  [rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4]]

C_m_n_c_Aerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

C_m_n_c_Aerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_m_n_c = ols(design, C_m_n_c_Aerobic_df).fit()
aov_table_m_n_c_Aerobic = anova_lm(model_m_n_c, typ=2)
print aov_table_m_n_c_Aerobic







C_all_Anaerobic_df = pd.DataFrame()
C_all_Anaerobic_df['tot_abs_difference'] =  [rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4]]
C_all_Anaerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

C_all_Anaerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_all_Anaerobic = ols(design, C_all_Anaerobic_df).fit()
aov_table_all_Anaerobic = anova_lm(model_all_Anaerobic, typ=2)
print aov_table_all_Anaerobic


C_m_n_c_Anaerobic_df = pd.DataFrame()
C_m_n_c_Anaerobic_df['tot_abs_difference'] =  [rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c'][4],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4]]

C_m_n_c_Anaerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

C_m_n_c_Anaerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_m_n_c = ols(design, C_m_n_c_Anaerobic_df).fit()
aov_table_m_n_c_Anaerobic = anova_lm(model_m_n_c, typ=2)
print aov_table_m_n_c_Anaerobic



C_eps_combined_Aerobic_df = pd.DataFrame()
C_eps_combined_Aerobic_df['tot_abs_difference'] = [rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][4],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][4],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][4],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][4]]

C_eps_combined_Aerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]

C_eps_combined_Aerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_eps = ols(design, C_eps_combined_Aerobic_df).fit()
aov_table_eps_combined_Aerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_combined_Aerobic

#Use eps_Aerobic_df for graph
C_eps_Aerobic_df = pd.DataFrame()
C_eps_Aerobic_df['tot_abs_difference'] = [rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][4],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][4]]

C_eps_Aerobic_df['Condition'] = [1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition)'
model_eps = ols(design, C_eps_Aerobic_df).fit()
aov_table_eps_Aerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_Aerobic


C_eps_original_model_Aerobic_df = pd.DataFrame()
C_eps_original_model_Aerobic_df['tot_abs_difference'] = [rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Aerobic_0.25'][4],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][4]]

C_eps_original_model_Aerobic_df['Condition'] = [1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition)'
model_eps = ols(design, C_eps_original_model_Aerobic_df).fit()
aov_table_eps_original_model_Aerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_original_model_Aerobic




C_eps_combined_Anaerobic_df = pd.DataFrame()
C_eps_combined_Anaerobic_df['tot_abs_difference'] = [rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][4],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][4],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][4],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][4]]

C_eps_combined_Anaerobic_df['Condition'] = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]

C_eps_combined_Anaerobic_df['Group'] = [1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_eps = ols(design, C_eps_combined_Anaerobic_df).fit()
aov_table_eps_combined_Anaerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_combined_Anaerobic

#Use eps_Anaerobic_df for graph
C_eps_Anaerobic_df = pd.DataFrame()
C_eps_Anaerobic_df['tot_abs_difference'] = [rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'][4],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'][4]]

C_eps_Anaerobic_df['Condition'] = [1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition)'
model_eps = ols(design, C_eps_Anaerobic_df).fit()
aov_table_eps_Anaerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_Anaerobic


C_eps_original_model_Anaerobic_df = pd.DataFrame()
C_eps_original_model_Anaerobic_df['tot_abs_difference'] = [rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151006_Anaerobic_0.25'][4],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][0],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][1],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][2],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][3],rxn_Anaerobic_tot_difference_per_dataset_S2_original_model_eps['151015_negative_control_0.25'][4]]

C_eps_original_model_Anaerobic_df['Condition'] = [1,1,1,1,1,2,2,2,2,2]

design = 'tot_abs_difference ~ C(Condition)'
model_eps = ols(design, C_eps_original_model_Anaerobic_df).fit()
aov_table_eps_original_model_Anaerobic = anova_lm(model_eps, typ=2)
print aov_table_eps_original_model_Anaerobic




C_all_Aerobic_vs_C_eps = pd.DataFrame()
C_all_Aerobic_vs_C_eps['tot_abs_difference'] = [rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][4]]


C_all_Aerobic_vs_C_eps['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2]

C_all_Aerobic_vs_C_eps['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,17,17,17,17,17]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_eps = ols(design, C_all_Aerobic_vs_C_eps).fit()
aov_table_all_Aerobic_vs_C_eps = anova_lm(model_eps, typ=2)
print aov_table_all_Aerobic_vs_C_eps



C_m_n_c_Aerobic_vs_C_eps = pd.DataFrame()
C_m_n_c_Aerobic_vs_C_eps['tot_abs_difference'] = [rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c'][4],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][0],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][1],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][2],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][3],rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'][4],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][0],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][1],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][2],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][3],rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'][4]]


C_m_n_c_Aerobic_vs_C_eps['Condition'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2]

C_m_n_c_Aerobic_vs_C_eps['Group'] = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9]

design = 'tot_abs_difference ~ C(Condition) + C(Group)'
model_eps = ols(design, C_m_n_c_Aerobic_vs_C_eps).fit()
aov_table_m_n_c_Aerobic_vs_C_eps = anova_lm(model_eps, typ=2)
print aov_table_m_n_c_Aerobic_vs_C_eps
#This is how you extract the p-value
aov_table_m_n_c_Aerobic_vs_C_eps.ix[0,3]





Aerobic_average_eps = np.mean(rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'])

Aerobic_average_m_n_c = np.mean(rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_average_all = np.mean(rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_negative_control_average_eps = np.mean(rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'])

Aerobic_negative_control_average_m_n_c = np.mean(rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_negative_control_average_all = np.mean(rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_average_eps = np.average(rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'])

Anaerobic_average_m_n_c = np.average(rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_average_all = np.average(rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_negative_control_average_eps = np.mean(rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'])

Anaerobic_negative_control_average_m_n_c = np.mean(rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_negative_control_average_all = np.mean(rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])





Aerobic_std_eps = np.std(rxn_Aerobic_tot_difference_per_dataset_S2_eps['151006_Aerobic_0.25'])

Aerobic_std_m_n_c = np.std(rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_std_all = np.std(rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_negative_control_std_eps = np.std(rxn_Aerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'])

Aerobic_negative_control_std_m_n_c = np.std(rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Aerobic_negative_control_std_all = np.std(rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']+rxn_Aerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_std_eps = np.std(rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151006_Anaerobic_0.25'])

Anaerobic_std_m_n_c = np.std(rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_std_all = np.std(rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_negative_control_std_eps = np.std(rxn_Anaerobic_tot_difference_per_dataset_S2_eps['151015_negative_control_0.25'])

Anaerobic_negative_control_std_m_n_c = np.std(rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])

Anaerobic_negative_control_std_all = np.std(rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']+rxn_Anaerobic_tot_difference_per['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns'])



n_groups = 3 
index = np.arange(n_groups)

Aerobic_average_individual = [Aerobic_average_eps, Aerobic_average_m_n_c, Aerobic_average_all]
Aerobic_negative_control_average_individual = [Aerobic_negative_control_average_eps, Aerobic_negative_control_average_m_n_c, Aerobic_negative_control_average_all]
Anaerobic_average_individual = [Anaerobic_average_eps, Anaerobic_average_m_n_c, Anaerobic_average_all]
Anaerobic_negative_control_average_individual = [Anaerobic_negative_control_average_eps, Anaerobic_negative_control_average_m_n_c, Anaerobic_negative_control_average_all]

Aerobic_std_individual = [Aerobic_std_eps, Aerobic_std_m_n_c, Aerobic_std_all]
Aerobic_negative_control_std_individual = [Aerobic_negative_control_std_eps, Aerobic_negative_control_std_m_n_c, Aerobic_negative_control_std_all]
Anaerobic_std_individual = [Anaerobic_std_eps, Anaerobic_std_m_n_c, Anaerobic_std_all]
Anaerobic_negative_control_std_individual = [Anaerobic_negative_control_std_eps, Anaerobic_negative_control_std_m_n_c, Anaerobic_negative_control_std_all]

fig, (ax0, ax1) = plt.subplots(2)

plt.close()
bar_width = 0.17
error_bar_length = np.repeat(0,len(Aerobic_average_individual))
rects1 = ax0.bar(index, Aerobic_average_individual, bar_width, alpha = 0.8, yerr = [error_bar_length, Aerobic_std_individual], capsize = 3, color = 'b', label = 'Aerobic')
rects2 = ax0.bar(index + bar_width, Aerobic_negative_control_average_individual, bar_width, alpha = 0.8, yerr = [error_bar_length, Aerobic_negative_control_std_individual], capsize = 3, color = 'r', label = 'Aerobic Neg') 
rects3 = ax0.bar(index + bar_width + bar_width, Anaerobic_average_individual, bar_width, alpha = 0.8, yerr = [error_bar_length, Anaerobic_std_individual], capsize = 3, color = 'g', label = 'Anaerobic') 
rects4  = ax0.bar(index + bar_width + bar_width + bar_width, Anaerobic_negative_control_average_individual, bar_width, alpha = 0.8, yerr = [error_bar_length, Anaerobic_negative_control_std_individual], capsize = 3, color = 'm', label = 'Anaerobic Neg')

props = {'connectionstyle':'bar','arrowstyle':'-','lw':2}
ax0.annotate('***', xy=(0+bar_width/4-0.01,23+max(Aerobic_average_individual[0]+Aerobic_std_individual[0],Aerobic_negative_control_average_individual[0]+Aerobic_negative_control_std_individual[0])),fontsize=9)
ax0.annotate('', xy=(0,10+max(Aerobic_average_individual[0]+Aerobic_std_individual[0],Aerobic_negative_control_average_individual[0]+Aerobic_negative_control_std_individual[0])), xytext=(0+bar_width,10+max(Aerobic_average_individual[0]+Aerobic_std_individual[0],Aerobic_negative_control_average_individual[0]+Aerobic_negative_control_std_individual[0])), arrowprops=props)

ax0.annotate('***', xy=(0+bar_width+bar_width+bar_width/4-0.01,23+max(Anaerobic_average_individual[0]+Anaerobic_std_individual[0],Anaerobic_negative_control_average_individual[0]+Anaerobic_negative_control_std_individual[0])),fontsize=9)
ax0.annotate('', xy=(0+bar_width+bar_width,10+max(Anaerobic_average_individual[0]+Anaerobic_std_individual[0],Anaerobic_negative_control_average_individual[0]+Anaerobic_negative_control_std_individual[0])), xytext=(0+bar_width+bar_width+bar_width,10+max(Anaerobic_average_individual[0]+Anaerobic_std_individual[0],Anaerobic_negative_control_average_individual[0]+Anaerobic_negative_control_std_individual[0])), arrowprops=props)

ax0.annotate('***', xy=(1+bar_width/4-0.01,23+max(Aerobic_average_individual[1]+Aerobic_std_individual[1],Aerobic_negative_control_average_individual[1]+Aerobic_negative_control_std_individual[1])),fontsize=9)
ax0.annotate('', xy=(1,10+max(Aerobic_average_individual[1]+Aerobic_std_individual[1],Aerobic_negative_control_average_individual[1]+Aerobic_negative_control_std_individual[1])), xytext=(1+bar_width,10+max(Aerobic_average_individual[1]+Aerobic_std_individual[1],Aerobic_negative_control_average_individual[1]+Aerobic_negative_control_std_individual[1])), arrowprops=props)

ax0.annotate('***', xy=(1+bar_width+bar_width+bar_width/4-0.01,23+max(Anaerobic_average_individual[1]+Anaerobic_std_individual[1],Anaerobic_negative_control_average_individual[1]+Anaerobic_negative_control_std_individual[1])),fontsize=9)
ax0.annotate('', xy=(1+bar_width+bar_width,10+max(Anaerobic_average_individual[1]+Anaerobic_std_individual[1],Anaerobic_negative_control_average_individual[1]+Anaerobic_negative_control_std_individual[1])), xytext=(1+bar_width+bar_width+bar_width,10+max(Anaerobic_average_individual[1]+Anaerobic_std_individual[1],Anaerobic_negative_control_average_individual[1]+Anaerobic_negative_control_std_individual[1])), arrowprops=props)

ax0.annotate('NS', xy=(2+bar_width/4-0.01,29+max(Aerobic_average_individual[2]+Aerobic_std_individual[2],Aerobic_negative_control_average_individual[2]+Aerobic_negative_control_std_individual[2])),fontsize=9)
ax0.annotate('', xy=(2,10+max(Aerobic_average_individual[2]+Aerobic_std_individual[2],Aerobic_negative_control_average_individual[2]+Aerobic_negative_control_std_individual[2])), xytext=(2+bar_width,10+max(Aerobic_average_individual[2]+Aerobic_std_individual[2],Aerobic_negative_control_average_individual[2]+Aerobic_negative_control_std_individual[2])), arrowprops=props)

ax0.annotate('NS', xy=(2+bar_width+bar_width+bar_width/4-0.01,29+max(Anaerobic_average_individual[2]+Anaerobic_std_individual[2],Anaerobic_negative_control_average_individual[2]+Anaerobic_negative_control_std_individual[2])),fontsize=9)
ax0.annotate('', xy=(2+bar_width+bar_width,10+max(Anaerobic_average_individual[2]+Anaerobic_std_individual[2],Anaerobic_negative_control_average_individual[2]+Anaerobic_negative_control_std_individual[2])), xytext=(2+bar_width+bar_width+bar_width,10+max(Anaerobic_average_individual[2]+Anaerobic_std_individual[2],Anaerobic_negative_control_average_individual[2]+Anaerobic_negative_control_std_individual[2])), arrowprops=props)


#label_diff(-0.17,0.17,'***',index+bar_width/2,Aerobic_average_individual)


#ax0.xlabel('Model Parameters',fontweight='bold',fontsize=10) 
ax0.set_ylabel('Total Absolute Flux\nDifference (mmol/h/gDW)',fontweight='bold',fontsize=10) 
ax0.tick_params(labelsize=9)
#plt.title('Comparison of Flux Predictions',fontweight='bold') 
ax0.set_xticks(index + 0.25)
ax0.set_xticklabels(('C_eps', 'C_m_n_c Cohort', 'C_all Cohort'), rotation=0)
a = ax0.legend(loc="upper right",fontsize=8)
a.get_frame().set_alpha(1)
a.set_title('Condition',prop={'weight':'bold','size':'9'})
b = ax0.axvspan(-0.25,0.75,alpha=0.1,color='orange')
c = ax0.axvspan(0.75,2.75,alpha=0.1,color='red')
d = ax0.legend([b,c],['EXAMO with Model\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],loc='upper left',fontsize=8)
d.get_frame().set_alpha(1)
d.set_title('Software',prop={'weight':'bold','size':'9'})

def create_proxy(label):
    line = matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black', mec='none', marker=r'$\mathregular{{{}}}$'.format(label))
    return line

labels = ['***','**','*','NS']
descriptions = ['p<0.001', 'p<0.01', 'p<0.05', 'p>0.05']
proxies = [create_proxy(item) for item in labels]
e = ax0.legend(proxies,descriptions,markerscale=2, loc='upper center',fontsize=8)
e.set_title('Significance',prop={'weight':'bold','size':'9'})


ax0.add_artist(a)
ax0.add_artist(d)
ax0.add_artist(e)
ax0.set_ylim(0,355)


n_groups = 2
index = np.arange(n_groups)
Aerobic_Anaerobic_average_comparing_software_eps = [Aerobic_average_eps, Anaerobic_average_eps] 
Aerobic_Anaerobic_average_comparing_software_m_n_c = [Aerobic_average_m_n_c, Anaerobic_average_m_n_c]
Aerobic_Anaerobic_average_comparing_software_all = [Aerobic_average_all, Anaerobic_average_all]

Aerobic_Anaerobic_std_comparing_software_eps = [Aerobic_std_eps, Anaerobic_std_eps] 
Aerobic_Anaerobic_std_comparing_software_m_n_c = [Aerobic_std_m_n_c, Anaerobic_std_m_n_c]
Aerobic_Anaerobic_std_comparing_software_all = [Aerobic_std_all, Anaerobic_std_all]

bar_width = 0.2
error_bar_length = np.repeat(0,len(Aerobic_Anaerobic_average_comparing_software_eps))
rects1 = ax1.bar(index, Aerobic_Anaerobic_average_comparing_software_eps, bar_width, alpha = 0.8, yerr = [error_bar_length, Aerobic_Anaerobic_std_comparing_software_eps], capsize = 3, color = 'orange', label = 'C_eps')
rects2 = ax1.bar(index + bar_width, Aerobic_Anaerobic_average_comparing_software_m_n_c, bar_width, alpha = 0.8, yerr = [error_bar_length, Aerobic_Anaerobic_std_comparing_software_m_n_c], capsize = 3, color = 'r', label = 'C_m_n_c Cohort') 
rects3 = ax1.bar(index + bar_width + bar_width, Aerobic_Anaerobic_average_comparing_software_all, bar_width, alpha = 0.8, yerr = [error_bar_length, Aerobic_Anaerobic_std_comparing_software_all], capsize = 3, color = 'darkred', label = 'C_all Cohort') 


props = {'connectionstyle':'bar','arrowstyle':'-', 'lw':2}
ax1.annotate('***', xy=(0+0.075,41.5+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_m_n_c[0]+Aerobic_Anaerobic_std_comparing_software_m_n_c[0])),fontsize=9)
ax1.annotate('', xy=(0,10+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_m_n_c[0]+Aerobic_Anaerobic_std_comparing_software_m_n_c[0])), xytext=(0+bar_width,10+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_m_n_c[0]+Aerobic_Anaerobic_std_comparing_software_m_n_c[0])), arrowprops=props)

ax1.annotate('***', xy=(0+bar_width-0.025,80+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_all[0]+Aerobic_Anaerobic_std_comparing_software_all[0])),fontsize=9)
ax1.annotate('', xy=(0,15+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_all[0]+Aerobic_Anaerobic_std_comparing_software_all[0])), xytext=(0+bar_width+bar_width,15+max(Aerobic_Anaerobic_average_comparing_software_eps[0]+Aerobic_Anaerobic_std_comparing_software_eps[0],Aerobic_Anaerobic_average_comparing_software_all[0]+Aerobic_Anaerobic_std_comparing_software_all[0])), arrowprops=props)

ax1.annotate('***', xy=(1+0.075,41.5+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_m_n_c[1]+Aerobic_Anaerobic_std_comparing_software_m_n_c[1])),fontsize=9)
ax1.annotate('', xy=(1,10+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_m_n_c[1]+Aerobic_Anaerobic_std_comparing_software_m_n_c[1])), xytext=(1+bar_width,10+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_m_n_c[1]+Aerobic_Anaerobic_std_comparing_software_m_n_c[1])), arrowprops=props)

ax1.annotate('***', xy=(1+bar_width-0.025,80+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_all[1]+Aerobic_Anaerobic_std_comparing_software_all[1])),fontsize=9)
ax1.annotate('', xy=(1,15+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_all[1]+Aerobic_Anaerobic_std_comparing_software_all[1])), xytext=(1+bar_width+bar_width,15+max(Aerobic_Anaerobic_average_comparing_software_eps[1]+Aerobic_Anaerobic_std_comparing_software_eps[1],Aerobic_Anaerobic_average_comparing_software_all[1]+Aerobic_Anaerobic_std_comparing_software_all[1])), arrowprops=props)

ax1.set_xlabel('Conditions',fontweight='bold',fontsize=10) 
ax1.set_ylabel('Total Absolute Flux\nDifference (mmol/h/gDW)',fontweight='bold',fontsize=10) 
ax1.tick_params(labelsize=9)
#plt.title('Comparison of Flux Predictions',fontweight='bold') 
ax1.set_xticks(index + 0.2)
ax1.set_xticklabels(('Aerobic', 'Anaerobic'), rotation=0)
a = ax1.legend(loc="upper right",fontsize=8)
a.get_frame().set_alpha(1)
a.set_title('Model Parameters',prop={'weight':'bold','size':'9'})

labels = ['***','**','*','NS']
descriptions = ['p<0.001', 'p<0.01', 'p<0.05', 'p>0.05']
proxies = [create_proxy(item) for item in labels]
e = ax1.legend(proxies,descriptions,markerscale=2, loc='upper center',fontsize=8)
e.set_title('Significance',prop={'weight':'bold','size':'9'})


ax1.add_artist(a)
ax1.add_artist(e)
ax1.set_ylim(0,355)


ax0.text(-0.74,355,'a',fontsize=22)
ax1.text(-0.365,355,'b',fontsize=22)


plt.tight_layout()
fig.subplots_adjust(left=0.11,bottom=0.1,top=0.95,right=0.99)
fig.set_size_inches(6.69,5.5)
fig.savefig('Figure_4.png',dpi=300)





#Create the data frames for the heatmap
rowvalues = ['GLUK;HEX1;SBTD_D2/x1','FBA/x4','PYK/x9','PDHm/x10','GND/x3','RPE/x5','TKT2/x6','TALA/x7','Csm/x11','ACONTm/x12','MDHm/x13','ME1m;ME2m/x14','PPCK/x15','ACS/x17','ALDD2y/x18','ALCD2ir/x19','G3PT/x20','PYRDC/x24','OAAt2m/x21','PYRt2m/x23','GLCt1/x25','ETOHt/x26','ACt2r/x27','GLYCt/x28','OAAt/x34','PYRt/x37','SUCCt2r/x38']
group_labels = ['C_orig', 'C_orig_eps', 'C_mod', 'C_mod_eps', 'C', 'C_Ex', 'C_g', 'C_g_Ex', 'C_lb', 'C_lb_EX', 'C_lb_g', 'C_lb_g_EX', 'C_m_n_c', 'C_m_n_c_EX', 'C_m_n_c_g', 'C_m_n_c_g_EX', 'C_m_n_c_lb', 'C_m_n_c_lb_EX', 'C_m_n_c_lb_g', 'C_m_n_c_lb_g_EX']
software_class = ['blue','blue','#eecb13','#eecb13','red','red','red','red','red','red','red','red','red','red','red','red','red','red','red','red']
Aerobic_df = pd.DataFrame()
Aerobic_df['C_orig'] = [rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151006_Aerobic_0.25']['x38']['abs_difference']]

Aerobic_df['C_orig_eps'] = [rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Aerobic_0.25']['x38']['abs_difference']] 

Aerobic_df['C_mod'] =  [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

Aerobic_df['C_mod_eps'] =  [rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151006_Aerobic_0.25']['x38']['abs_difference']]

Aerobic_df['C'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x38']['abs_difference']]


Aerobic_df['C_Ex'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_g'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x38']['abs_difference']]

Aerobic_df['C_g_Ex'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_lb'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x38']['abs_difference']]

Aerobic_df['C_lb_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_lb_g'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x38']['abs_difference']]

Aerobic_df['C_lb_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_m_n_c'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_g'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_lb'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_lb_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_lb_g'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x38']['abs_difference']]

Aerobic_df['C_m_n_c_lb_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]





Aerobic_negative_control_df = pd.DataFrame()
Aerobic_negative_control_df['C_orig'] = [rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x38']['abs_difference']]

Aerobic_negative_control_df['C_orig_eps'] = [rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x38']['abs_difference']] 

Aerobic_negative_control_df['C_mod'] =  [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

Aerobic_negative_control_df['C_mod_eps'] =  [rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x38']['abs_difference']]

Aerobic_negative_control_df['C'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x38']['abs_difference']]


Aerobic_negative_control_df['C_Ex'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_g'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x38']['abs_difference']]

Aerobic_negative_control_df['C_g_Ex'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_lb'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['x38']['abs_difference']]

Aerobic_negative_control_df['C_lb_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_lb_g'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['x38']['abs_difference']]

Aerobic_negative_control_df['C_lb_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_g'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_lb'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_lb_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_lb_g'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['x38']['abs_difference']]

Aerobic_negative_control_df['C_m_n_c_lb_g_EX'] =  [rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Aerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]




Anaerobic_df = pd.DataFrame()
Anaerobic_df['C_orig'] = [rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151006_Anaerobic_0.25']['x38']['abs_difference']]

Anaerobic_df['C_orig_eps'] = [rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['x38']['abs_difference']] 

Anaerobic_df['C_mod'] =  [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

Anaerobic_df['C_mod_eps'] =  [rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151006_Anaerobic_0.25']['x38']['abs_difference']]

Anaerobic_df['C'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['x38']['abs_difference']]


Anaerobic_df['C_Ex'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_g'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x38']['abs_difference']]

Anaerobic_df['C_g_Ex'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_lb'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x38']['abs_difference']]

Anaerobic_df['C_lb_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_lb_g'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x38']['abs_difference']]

Anaerobic_df['C_lb_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_g'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_lb'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_lb_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_lb_g'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x38']['abs_difference']]

Anaerobic_df['C_m_n_c_lb_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]





Anaerobic_negative_control_df = pd.DataFrame()
Anaerobic_negative_control_df['C_orig'] = [rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model['151015_negative_control_0.25']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_orig_eps'] = [rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_original_model_eps['151015_negative_control_0.25']['x38']['abs_difference']] 

Anaerobic_negative_control_df['C_mod'] =  [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]

Anaerobic_negative_control_df['C_mod_eps'] =  [rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict_dataset_S2_eps['151015_negative_control_0.25']['x38']['abs_difference']]

Anaerobic_negative_control_df['C'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['x38']['abs_difference']]


Anaerobic_negative_control_df['C_Ex'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_g'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_g_Ex'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_lb'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_lb_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_lb_g'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_lb_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_g'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_lb'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_lb_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_lb_g'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['x38']['abs_difference']]

Anaerobic_negative_control_df['C_m_n_c_lb_g_EX'] =  [rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x1']['abs_difference'], rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x4']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x9']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x10']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x3']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x5']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x6']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x7']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x11']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x12']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x13']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x14']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x15']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x17']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x18']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x19']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x20']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x24']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x21']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x23']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x25']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x26']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x27']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x28']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x34']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x37']['abs_difference'],
rxn_Anaerobic_metabolicState_dict['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['x38']['abs_difference']]


#Create the labels for the columns of the heatmap based on the total flux difference
#To achieve rounding with sig figs to a certian number of places, I used Randle Taylor's code from http://randlet.com/blog/python-significant-figures-format/
def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)


Aerobic_totals = [to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_original_model['151006_Aerobic_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_original_model_eps['151006_Aerobic_0.25']['tot_abs_difference']),3),
'NA',
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_eps['151006_Aerobic_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151006_Aerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3)]


Aerobic_negative_control_totals = [to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_original_model['151015_negative_control_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_original_model_eps['151015_negative_control_0.25']['tot_abs_difference']),3),
'NA',
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference_dataset_S2_eps['151015_negative_control_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Aerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Aerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3)]


Anaerobic_totals = [to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_original_model['151006_Anaerobic_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_original_model_eps['151006_Anaerobic_0.25']['tot_abs_difference']),3),
'NA',
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_eps['151006_Anaerobic_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151006_Anaerobic_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3)]



Anaerobic_negative_control_totals = [to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_original_model['151015_negative_control_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_original_model_eps['151015_negative_control_0.25']['tot_abs_difference']),3),
'NA',
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference_dataset_S2_eps['151015_negative_control_0.25']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c']['tot_abs_difference']),3),
to_precision('%s' % float('%.3g' % rxn_Anaerobic_tot_difference['151015_negative_control_0.25_iMM904_NADcorrected_1127_FTHFLm_Rintala_Anaerobic_g_m_n_c_EXrxns_EXtrrxns_Othertrrxns']['tot_abs_difference']),3)]



#This top section is just for indexing for creating the secondary legend later
#Now import seaborn for heatmap capabilities
import seaborn as sns
plt.close()
n_groups = 20
index = np.arange(n_groups)
color1 = plt.bar(index,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], color="blue")
color2 = plt.bar(index,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], color="#eecb13") 
color3 = plt.bar(index,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], color="red")

#Now create the heatmaps and the subplots
fig, ax = plt.subplots(2,2)
plot1 = sns.heatmap(Aerobic_df, cmap=plt.cm.YlOrRd, ax=ax[0,0], cbar=0, vmin=0, vmax=5, xticklabels = False)
plot2 = sns.heatmap(Aerobic_negative_control_df,cmap=plt.cm.YlOrRd, ax=ax[0,1], cbar=0, vmin=0, vmax=5, xticklabels = False, yticklabels=False)
plot3 = sns.heatmap(Anaerobic_df,cmap=plt.cm.YlOrRd, ax=ax[1,0], cbar=0, vmin=0, vmax=5)
plot4 = sns.heatmap(Anaerobic_negative_control_df,cmap=plt.cm.YlOrRd, ax=ax[1,1], cbar=1, vmin=0, vmax=5, yticklabels=False, cbar_ax = fig.add_axes([.92,.375,0.02,.4]))
#plot4 = sns.heatmap(Anaerobic_negative_control_df,cmap=plt.cm.YlOrRd, ax=ax[1,1], cbar=1, vmin=0, vmax=5, yticklabels=False, cbar_ax = fig.add_axes([.91,.375,0.03,.4]), cbar_kws={'label': 'Absolute Flux Difference (Experimental - Simulated) (mmol/h/gDW)'})

ax[1,0].axes.set_xticklabels(group_labels,rotation=90,fontsize=9)
for color,tick in zip(software_class,ax[1,0].xaxis.get_major_ticks()):
    tick.label1.set_color(color) #set the color property

ax[0,0].axes.set_yticklabels(reversed(rowvalues),rotation=0, fontsize=7)


xl, xh=ax[0,0].get_xlim()
left=xl-(xh-xl)*1.09
left_line = xl-(xh-xl)*1.1
right=xh

Lines=ax[0,0].hlines([7,9,10,11,14,16,17,19,23], left_line,right)
Lines.set_clip_on(False)
ax[0,0].set_xlim((xl,xh))


ax[0,0].text(left, 0.1, "Transport Extracullar", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 7.1, "Transport Mitochondrial", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 9.1, "Complex Alcohol Metab.", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 10.1, "Glycerolipid Metab.", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 11.1, "Pyruvate Metab.", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 14.1, "Anaplerotic Reactions", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 16.1, "Ox. Phosphorylation", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 17.1, "Citric Acid Cycle", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 19.1, "Pentose Phosphate", horizontalalignment="left", fontsize=9)
ax[0,0].text(left, 23.1, "Glycolysis", horizontalalignment="left", fontsize=9)

ax[0,0].text(left+3, 28, "Total Absolute Flux Difference:", horizontalalignment = "left", fontsize = 8, weight = "bold")


ax[1,0].axes.set_yticklabels(reversed(rowvalues),rotation=0, fontsize=9)

Lines=ax[1,0].hlines([7,9,10,11,14,16,17,19,23], left_line,right)
Lines.set_clip_on(False)
ax[1,0].set_xlim((xl,xh))

ax[1,0].text(left, 0.1, "Transport Extracullar", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 7.1, "Transport Mitochondrial", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 9.1, "Complex Alcohol Metab.", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 10.1, "Glycerolipid Metab.", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 11.1, "Pyruvate Metab.", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 14.1, "Anaplerotic Reactions", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 16.1, "Ox. Phosphorylation", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 17.1, "Citric Acid Cycle", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 19.1, "Pentose Phosphate", horizontalalignment="left", fontsize=9)
ax[1,0].text(left, 23.1, "Glycolysis", horizontalalignment="left", fontsize=9)

ax[1,0].text(left+3, 28, "Total Absolute Flux Difference:", horizontalalignment = "left", fontsize = 8, weight = "bold")

#Add hlines to other plots
ax[0,1].hlines([7,9,10,11,14,16,17,19,23], xl, xh)
ax[1,1].hlines([7,9,10,11,14,16,17,19,23], xl, xh)


ax[0,0].text(0, 27.5, Aerobic_totals[0], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(1, 27.5, Aerobic_totals[1], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(2, 27.5, Aerobic_totals[2], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(3, 27.5, Aerobic_totals[3], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(4, 27.5, Aerobic_totals[4], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(5, 27.5, Aerobic_totals[5], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(6, 27.5, Aerobic_totals[6], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(7, 27.5, Aerobic_totals[7], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(8, 27.5, Aerobic_totals[8], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(9, 27.5, Aerobic_totals[9], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(10, 27.5, Aerobic_totals[10], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(11, 27.5, Aerobic_totals[11], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(12, 27.5, Aerobic_totals[12], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(13, 27.5, Aerobic_totals[13], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(14, 27.5, Aerobic_totals[14], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(15, 27.5, Aerobic_totals[15], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(16, 27.5, Aerobic_totals[16], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(17, 27.5, Aerobic_totals[17], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(18, 27.5, Aerobic_totals[18], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,0].text(19, 27.5, Aerobic_totals[19], rotation=90, fontsize=7, verticalalignment='bottom')


ax[0,1].text(0, 27.5, Aerobic_negative_control_totals[0], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(1, 27.5, Aerobic_negative_control_totals[1], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(2, 27.5, Aerobic_negative_control_totals[2], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(3, 27.5, Aerobic_negative_control_totals[3], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(4, 27.5, Aerobic_negative_control_totals[4], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(5, 27.5, Aerobic_negative_control_totals[5], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(6, 27.5, Aerobic_negative_control_totals[6], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(7, 27.5, Aerobic_negative_control_totals[7], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(8, 27.5, Aerobic_negative_control_totals[8], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(9, 27.5, Aerobic_negative_control_totals[9], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(10, 27.5, Aerobic_negative_control_totals[10], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(11, 27.5, Aerobic_negative_control_totals[11], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(12, 27.5, Aerobic_negative_control_totals[12], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(13, 27.5, Aerobic_negative_control_totals[13], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(14, 27.5, Aerobic_negative_control_totals[14], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(15, 27.5, Aerobic_negative_control_totals[15], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(16, 27.5, Aerobic_negative_control_totals[16], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(17, 27.5, Aerobic_negative_control_totals[17], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(18, 27.5, Aerobic_negative_control_totals[18], rotation=90, fontsize=7, verticalalignment='bottom')
ax[0,1].text(19, 27.5, Aerobic_negative_control_totals[19], rotation=90, fontsize=7, verticalalignment='bottom')


ax[1,0].text(0, 27.5, Anaerobic_totals[0], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(1, 27.5, Anaerobic_totals[1], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(2, 27.5, Anaerobic_totals[2], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(3, 27.5, Anaerobic_totals[3], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(4, 27.5, Anaerobic_totals[4], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(5, 27.5, Anaerobic_totals[5], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(6, 27.5, Anaerobic_totals[6], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(7, 27.5, Anaerobic_totals[7], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(8, 27.5, Anaerobic_totals[8], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(9, 27.5, Anaerobic_totals[9], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(10, 27.5, Anaerobic_totals[10], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(11, 27.5, Anaerobic_totals[11], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(12, 27.5, Anaerobic_totals[12], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(13, 27.5, Anaerobic_totals[13], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(14, 27.5, Anaerobic_totals[14], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(15, 27.5, Anaerobic_totals[15], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(16, 27.5, Anaerobic_totals[16], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(17, 27.5, Anaerobic_totals[17], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(18, 27.5, Anaerobic_totals[18], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,0].text(19, 27.5, Anaerobic_totals[19], rotation=90, fontsize=7, verticalalignment='bottom')


ax[1,1].text(0, 27.5, Anaerobic_negative_control_totals[0], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(1, 27.5, Anaerobic_negative_control_totals[1], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(2, 27.5, Anaerobic_negative_control_totals[2], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(3, 27.5, Anaerobic_negative_control_totals[3], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(4, 27.5, Anaerobic_negative_control_totals[4], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(5, 27.5, Anaerobic_negative_control_totals[5], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(6, 27.5, Anaerobic_negative_control_totals[6], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(7, 27.5, Anaerobic_negative_control_totals[7], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(8, 27.5, Anaerobic_negative_control_totals[8], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(9, 27.5, Anaerobic_negative_control_totals[9], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(10, 27.5, Anaerobic_negative_control_totals[10], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(11, 27.5, Anaerobic_negative_control_totals[11], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(12, 27.5, Anaerobic_negative_control_totals[12], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(13, 27.5, Anaerobic_negative_control_totals[13], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(14, 27.5, Anaerobic_negative_control_totals[14], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(15, 27.5, Anaerobic_negative_control_totals[15], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(16, 27.5, Anaerobic_negative_control_totals[16], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(17, 27.5, Anaerobic_negative_control_totals[17], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(18, 27.5, Anaerobic_negative_control_totals[18], rotation=90, fontsize=7, verticalalignment='bottom')
ax[1,1].text(19, 27.5, Anaerobic_negative_control_totals[19], rotation=90, fontsize=7, verticalalignment='bottom')


ax[1,1].axes.set_xticklabels(group_labels,rotation=90,fontsize=8.5)
for color,tick in zip(software_class,ax[1,1].xaxis.get_major_ticks()):
    tick.label1.set_color(color) #set the color property
    
    
ax[1,0].set_xlabel('Gene Rules for Condition',fontsize=10,labelpad=15, fontweight='bold')
ax[1,1].set_xlabel('Negative Control Gene Rules\nfor Condition',fontsize=10,labelpad=15, fontweight='bold')
ax[0,1].set_ylabel('Aerobic', fontsize=10, rotation=270, labelpad=10, fontweight='bold')
ax[0,1].yaxis.set_label_position("right")
ax[1,1].set_ylabel('Anaerobic', fontsize=10, rotation=270, labelpad=10, fontweight='bold')
ax[1,1].yaxis.set_label_position("right")
ax[1,1].collections[0].colorbar.set_ticks([0, 1, 2, 3, 4, 5])
ax[1,1].collections[0].colorbar.set_ticklabels(['0', '1', '2', '3', '4', '>5'])
#ax[1,1].collections[0].colorbar.ax.set_title('Absolute Flux Difference (Experimental - Simulated) (mmol/h/gDW)',fontsize=10,ha='right',rotation=90)
ax[1,1].collections[0].colorbar.ax.tick_params(labelsize=8)

ax[1,1].text(26.5, 9.5, 'Absolute Flux Difference (Experimental - Simulated) (mmol/h/gDW)', rotation=90, fontsize=10, verticalalignment='bottom')

ax[0,0].tick_params(labelsize=7)
ax[0,1].tick_params(labelsize=7)
ax[1,0].tick_params(labelsize=7)
ax[1,1].tick_params(labelsize=7)

a = plt.legend([color1,color2,color3],['EXAMO','EXAMO with Model Converted\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],bbox_to_anchor=[-30.5,-0.88,0.03,.4],fontsize=9)
#a = plt.legend([color1,color2,color3],['EXAMO','EXAMO with model converted\nby EXAMO-ARC.V.1','EXAMO-ARC.V.1'],bbox_to_anchor=[-21.5,-0.54,0.03,.4])
a.set_title('Software',prop={'weight':'bold', 'size':10})
plt.gca().add_artist(a)

fig.subplots_adjust(hspace=7)

plt.tight_layout(rect=[0.13,-0.1,0.93,0.99])
#plt.tight_layout(rect=[0.09,-0.1,0.9,1])
#plt.subplots_adjust(top=1.04)
fig.set_size_inches(6.69,8)
fig.savefig('Figure_5.png',dpi=300)
