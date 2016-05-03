#This is an example of how to create the dictionaries used for nucleotide conversions that are converted into a pickle file. 
import cPickle as pickle

M_amp_c = {'M_adp_c': ['M_amp_c'], 'M_atp_c': ['M_amp_c'], 'M_damp_c': ['M_amp_c'], 'M_dadp_c': ['M_amp_c'], 'M_datp_c': ['M_amp_c']}
M_amp_m = {'M_adp_m': ['M_amp_m'], 'M_atp_m': ['M_amp_m']}
M_amp_x = {'M_adp_x': ['M_amp_x'], 'M_atp_x': ['M_amp_x']}
M_gmp_c = {'M_gdp_c': ['M_gmp_c'], 'M_gtp_c': ['M_gmp_c'], 'M_dgmp_c': ['M_gmp_c'], 'M_dgdp_c': ['M_gmp_c'], 'M_dgtp_c': ['M_gmp_c']}
M_cmp_c = {'M_cdp_c': ['M_cmp_c'], 'M_ctp_c': ['M_cmp_c'], 'M_dcmp_c': ['M_cmp_c'], 'M_dcdp_c': ['M_cmp_c'], 'M_dctp_c': ['M_cmp_c']}
M_cmp_m = {'M_cdp_m': ['M_cmp_m'], 'M_ctp_m': ['M_cmp_m']}
M_dtmp_c = {'M_dtdp_c': ['M_dtmp_c'], 'M_dttp_c': ['M_dtmp_c']}
M_ump_c = {'M_udp_c': ['M_ump_c'], 'M_utp_c': ['M_ump_c'], 'M_dump_c': ['M_ump_c'], 'M_dudp_c': ['M_ump_c'], 'M_dutp_c': ['M_ump_c']}


nucleotide_conversions = {'M_amp_c': M_amp_c, 'M_amp_m': M_amp_m, 'M_amp_x': M_amp_x, 'M_gmp_c': M_gmp_c, 'M_cmp_c': M_cmp_c, 'M_cmp_m': M_cmp_m, 'M_dtmp_c': M_dtmp_c, 'M_ump_c': M_ump_c}

pickle.dump({'nucleotide_conversions': nucleotide_conversions}, open('150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.pkl', 'wb'))
