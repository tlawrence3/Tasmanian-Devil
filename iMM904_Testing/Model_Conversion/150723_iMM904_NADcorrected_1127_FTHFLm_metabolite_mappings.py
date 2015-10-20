#This is an example of how to create the dictionaries used for nucleotide conversions that are converted into a pickle file. 
import cPickle as pickle

metabolite_mappings = {'M_fad_m':['M_fadh2_m'], 'M_fadh2_m':['M_fad_m'], 'M_nad_c':['M_nadh_c'], 'M_nad_m':['M_nadh_m'], 'M_nad_x':['M_nadh_x'], 'M_nadh_c':['M_nad_c'], 'M_nadh_m':['M_nad_m'], 'M_nadh_x':['M_nad_x'], 'M_nadp_c':['M_nadph_c'], 'M_nadp_m':['M_nadph_m'], 'M_nadp_r':['M_nadph_r'], 'M_nadph_c':['M_nadp_c'], 'M_nadph_m':['M_nadp_m'], 'M_nadph_r':['M_nadp_r']}

pickle.dump({'metabolite_mappings': metabolite_mappings}, open('150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.pkl', 'wb'))
