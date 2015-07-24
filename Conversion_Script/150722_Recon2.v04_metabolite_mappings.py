#This is an example of how to create the dictionaries used for nucleotide conversions that are converted into a pickle file. 
import cPickle as pickle

metabolite_mappings = {'M_fad_c_':['M_fadh2_c_'], 'M_fad_m_':['M_fadh2_m_'], 'M_fad_x_':['M_fadh2_x_'], 'M_fad_r_': ['M_fadh2_r_'], 'M_fadh2_c_':['M_fad_c_'], 'M_fadh2_m_':['M_fad_m_'], 'M_fadh2_x_':['M_fad_x_'], 'M_fadh2_r_': ['M_fad_r_'], 'M_nad_c_':['M_nadh_c_'], 'M_nad_m_':['M_nadh_m_'], 'M_nad_r_':['M_nadh_r_'], 'M_nad_x_':['M_nadh_x_'], 'M_nadh_c_':['M_nad_c_'], 'M_nadh_m_':['M_nad_m_'],'M_nadh_r_':['M_nad_r_'], 'M_nadh_x_':['M_nad_x_'], 'M_nadp_c_':['M_nadph_c_'], 'M_nadp_m_':['M_nadph_m_'], 'M_nadp_r_':['M_nadph_r_'], 'M_nadp_x_':['M_nadph_x_'], 'M_nadp_l_': ['M_nadph_l_'], 'M_nadp_n_': ['M_nadph_n_'], 'M_nadph_c_':['M_nadp_c_'], 'M_nadph_m_':['M_nadp_m_'], 'M_nadph_r_':['M_nadp_r_'], 'M_nadph_x_':['M_nadp_x_'], 'M_nadph_l_': ['M_nadp_l_'], 'M_nadph_n_': ['M_nadp_n_']}

pickle.dump({'metabolite_mappings': metabolite_mappings}, open('150722_Recon2.v04_metabolite_mappings.pkl', 'wb'))
