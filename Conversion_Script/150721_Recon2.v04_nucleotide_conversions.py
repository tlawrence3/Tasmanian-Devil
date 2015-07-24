#This is an example of how to create the dictionaries used for nucleotide conversions that are converted into a pickle file. 
import cPickle as pickle

M_amp_c_ = {'M_adp_c_': ['M_amp_c_'], 'M_atp_c_': ['M_amp_c_'], 'M_damp_c_': ['M_amp_c_'], 'M_dadp_c_': ['M_amp_c_'], 'M_datp_c_': ['M_amp_c_'], 'M_camp_c_': ['M_amp_c_']}
M_amp_e_ = {'M_adp_e_': ['M_amp_e_'], 'M_atp_e_': ['M_amp_e_'], 'M_camp_e_': ['M_amp_e_']}
M_amp_g_ = {'M_adp_g_': ['M_amp_g_'], 'M_atp_g_': ['M_amp_g_'], 'M_camp_g_': ['M_amp_g_']}
M_amp_l_ = {'M_adp_l_': ['M_amp_l_'], 'M_atp_l_': ['M_amp_l_'], 'M_damp_l_': ['M_amp_l_']}
M_amp_m_ = {'M_adp_m_': ['M_amp_m_'], 'M_atp_m_': ['M_amp_m_'], 'M_dadp_m_': ['M_amp_m_'], 'M_datp_m_': ['M_amp_m_']}
M_dadp_n_ = {'M_adp_n_': ['M_dadp_n_'], 'M_atp_n_': ['M_dadp_n_'], 'M_datp_n_': ['M_dadp_n_']}
M_amp_r_ = {'M_adp_r_': ['M_amp_r_'], 'M_atp_r_': ['M_amp_r_']}
M_amp_x_ = {'M_adp_x_': ['M_amp_x_'], 'M_atp_x_': ['M_amp_x_']}

M_cmp_c_ = {'M_cdp_c_': ['M_cmp_c_'], 'M_ctp_c_': ['M_cmp_c_'], 'M_dcmp_c_': ['M_cmp_c_'], 'M_dcdp_c_': ['M_cmp_c_'], 'M_dctp_c_': ['M_cmp_c_']}
M_cmp_e_ = {'M_cdp_e_': ['M_cmp_e_'], 'M_ctp_e_': ['M_cmp_e_']}
M_cmp_l_ = {'M_dcmp_l_': ['M_cmp_l_']}
M_cmp_m_ = {'M_cdp_m_': ['M_cmp_m_'], 'M_ctp_m_': ['M_cmp_m_'], 'M_dcmp_m_': ['M_cmp_m_'], 'M_dcdp_m_': ['M_cmp_m_'], 'M_dctp_m_': ['M_cmp_m_']}
M_cmp_n_ = {'M_cdp_n_': ['M_cmp_n_'], 'M_ctp_n_': ['M_cmp_n_'], 'M_dcmp_n_': ['M_cmp_n_'], 'M_dcdp_n_': ['M_cmp_n_'], 'M_dctp_n_': ['M_cmp_n_']}

M_dtmp_c_ = {'M_dtdp_c_': ['M_dtmp_c_'], 'M_dttp_c_': ['M_dtmp_c_']}
M_dtmp_e_ = {'M_dtdp_e_': ['M_dtmp_e_'], 'M_dttp_e_': ['M_dtmp_e_']}
M_dtmp_m_ = {'M_dtdp_m_': ['M_dtmp_m_'], 'M_dttp_m_': ['M_dtmp_m_']}
M_dtmp_n_ = {'M_dtdp_n_': ['M_dtmp_n_'], 'M_dttp_n_': ['M_dtmp_n_']}

M_gmp_c_ = {'M_gdp_c_': ['M_gmp_c_'], 'M_gtp_c_': ['M_gmp_c_'], 'M_dgmp_c_': ['M_gmp_c_'], 'M_dgdp_c_': ['M_gmp_c_'], 'M_dgtp_c_': ['M_gmp_c_'], 'M_35cgmp_c_': ['M_gmp_c_']}
M_gmp_e_ = {'M_gdp_e_': ['M_gmp_e_'], 'M_gtp_e_': ['M_gmp_e_'], 'M_dgmp_e_': ['M_gmp_e_'], 'M_dgtp_e_': ['M_gmp_e_'], 'M_35cgmp_e_': ['M_gmp_e_']}
M_gmp_g_ = {'M_gdp_g_': ['M_gmp_g_'], 'M_35cgmp_g_': ['M_gmp_g_']}
M_gmp_l_ = {'M_dgmp_l_': ['M_gmp_l_']}
M_gmp_m_ = {'M_gdp_m_': ['M_gmp_m_'], 'M_gtp_m_': ['M_gmp_m_'], 'M_dgmp_m_': ['M_gmp_m_'], 'M_dgdp_m_': ['M_gmp_m_'], 'M_dgtp_m_': ['M_gmp_m_']}
M_gmp_n_ = {'M_gdp_n_': ['M_gmp_n_'], 'M_gtp_n_': ['M_gmp_n_'], 'M_dgdp_n_': ['M_gmp_n_'], 'M_dgtp_n_': ['M_gmp_n_'], 'M_35cgmp_n_': ['M_gmp_n_']}

M_ump_c_ = {'M_udp_c_': ['M_ump_c_'], 'M_utp_c_': ['M_ump_c_'], 'M_dump_c_': ['M_ump_c_'], 'M_dudp_c_': ['M_ump_c_'], 'M_dutp_c_': ['M_ump_c_']}
M_ump_e_ = {'M_udp_e_': ['M_ump_e_'], 'M_utp_e_': ['M_ump_e_']}
M_ump_g_ = {'M_udp_g_': ['M_ump_g_']}
M_ump_l_ = {'M_udp_l_': ['M_ump_l_']}
M_ump_m_ = {'M_udp_m_': ['M_ump_m_'], 'M_utp_m_': ['M_ump_m_'], 'M_dump_m_': ['M_ump_m_'], 'M_dudp_m_': ['M_ump_m_'], 'M_dutp_m_': ['M_ump_m_']}
M_ump_n_ = {'M_udp_n_': ['M_ump_n_'], 'M_utp_n_': ['M_ump_n_'], 'M_dump_n_': ['M_ump_n_'], 'M_dudp_n_': ['M_ump_n_'], 'M_dutp_n_': ['M_ump_n_']}
M_ump_r_ = {'M_udp_r_': ['M_ump_r_']}


M_imp_c_ = {'M_idp_c_': ['M_imp_c_'], 'M_itp_c_': ['M_imp_c_'], 'M_dimp_c_': ['M_imp_c_'], 'M_didp_c_': ['M_imp_c_'], 'M_ditp_c_': ['M_imp_c_']}
M_imp_e_ = {'M_idp_e_': ['M_imp_e_'], 'M_itp_e_': ['M_imp_e_']}
M_imp_m_ = {'M_idp_m_': ['M_imp_m_'], 'M_itp_m_': ['M_imp_m_'], 'M_didp_m_': ['M_imp_m_'], 'M_ditp_m_': ['M_imp_m_']}
M_didp_n_ = {'M_idp_n_': ['M_didp_n_'], 'M_itp_n_': ['M_didp_n_'], 'M_ditp_n_': ['M_didp_n_']}


nucleotide_conversions = {'M_amp_c_': M_amp_c_, 'M_amp_e_': M_amp_e_, 'M_amp_g_': M_amp_g_, 'M_amp_l_': M_amp_l_, 'M_amp_m_': M_amp_m_, 'M_dadp_n_': M_dadp_n_, 'M_amp_r_': M_amp_r_, 'M_amp_x_': M_amp_x_, 'M_cmp_c_': M_cmp_c_, 'M_cmp_e_': M_cmp_e_, 'M_cmp_l_': M_cmp_l_, 'M_cmp_m_': M_cmp_m_, 'M_cmp_n_': M_cmp_n_, 'M_dtmp_c_': M_dtmp_c_, 'M_dtmp_e_': M_dtmp_e_, 'M_dtmp_m_': M_dtmp_m_, 'M_dtmp_n_': M_dtmp_n_, 'M_gmp_c_': M_gmp_c_, 'M_gmp_e_': M_gmp_e_, 'M_gmp_g_': M_gmp_g_, 'M_gmp_l_': M_gmp_l_, 'M_gmp_m_': M_gmp_m_, 'M_gmp_n_': M_gmp_n_, 'M_ump_c_': M_ump_c_, 'M_ump_e_': M_ump_e_, 'M_ump_g_': M_ump_g_, 'M_ump_l_': M_ump_l_, 'M_ump_m_': M_ump_m_, 'M_ump_n_': M_ump_n_, 'M_ump_r_': M_ump_r_, 'M_imp_c_': M_imp_c_, 'M_imp_e_': M_imp_e_, 'M_imp_m_': M_imp_m_, 'M_didp_n_': M_didp_n_}

pickle.dump({'nucleotide_conversions': nucleotide_conversions}, open('150721_Recon2.v04_nucleotide_conversions.pkl', 'wb'))
