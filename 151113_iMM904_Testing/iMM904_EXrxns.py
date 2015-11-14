#This is an example of how to create the list for the extracellular reactions that is converted into a pickle file. 
import cPickle as pickle

EXrxns = ['R_EX_13BDglcn_e_', 'R_EX_2hb_e_', 'R_EX_2mbac_e_', 'R_EX_2mbald_e_', 'R_EX_2mbtoh_e_', 'R_EX_2mppal_e_', 'R_EX_2phetoh_e_', 'R_EX_3c3hmp_e_', 'R_EX_3mbald_e_', 'R_EX_3mop_e_', 'R_EX_4abut_e_', 'R_EX_4abz_e_', 'R_EX_5aop_e_', 'R_EX_Nbfortyr_e_', 'R_EX_abt_e_', 'R_EX_ac_e_', 'R_EX_acald_e_', 'R_EX_aces_e_', 'R_EX_ade_e_', 'R_EX_adn_e_', 'R_EX_akg_e_', 'R_EX_ala_L_e_', 'R_EX_alltn_e_', 'R_EX_alltt_e_', 'R_EX_amet_e_', 'R_EX_arab_L_e_', 'R_EX_arg_L_e_', 'R_EX_asn_L_e_', 'R_EX_asp_L_e_', 'R_EX_btd_RR_e_', 'R_EX_chol_e_', 'R_EX_cit_e_', 'R_EX_co2_e_', 'R_EX_csn_e_', 'R_EX_cys_L_e_', 'R_EX_cytd_e_', 'R_EX_dad_2_e_', 'R_EX_dca_e_', 'R_EX_dcyt_e_', 'R_EX_ddca_e_', 'R_EX_dgsn_e_', 'R_EX_din_e_', 'R_EX_dttp_e_', 'R_EX_duri_e_', 'R_EX_epist_e_', 'R_EX_epistest_SC_e_', 'R_EX_ergst_e_', 'R_EX_ergstest_SC_e_', 'R_EX_etha_e_', 'R_EX_etoh_e_', 'R_EX_fe2_e_', 'R_EX_fecost_e_', 'R_EX_fecostest_SC_e_', 'R_EX_fmn_e_', 'R_EX_for_e_', 'R_EX_fru_e_', 'R_EX_fum_e_', 'R_EX_g3pc_e_', 'R_EX_g3pi_e_', 'R_EX_gal_e_', 'R_EX_galur_e_', 'R_EX_gam6p_e_', 'R_EX_gcald_e_', 'R_EX_glc_e_', 'R_EX_gln_L_e_', 'R_EX_glu_L_e_', 'R_EX_glx_e_', 'R_EX_gly_e_', 'R_EX_glyc_e_', 'R_EX_gsn_e_', 'R_EX_gthox_e_', 'R_EX_gthrd_e_', 'R_EX_gua_e_', 'R_EX_h2o_e_', 'R_EX_h_e_', 'R_EX_hdca_e_', 'R_EX_hdcea_e_', 'R_EX_hexc_e_', 'R_EX_his_L_e_', 'R_EX_hxan_e_', 'R_EX_iamac_e_', 'R_EX_iamoh_e_', 'R_EX_ibutac_e_', 'R_EX_ibutoh_e_', 'R_EX_id3acald_e_', 'R_EX_ile_L_e_', 'R_EX_ind3eth_e_', 'R_EX_inost_e_', 'R_EX_ins_e_', 'R_EX_lac_D_e_', 'R_EX_lac_L_e_', 'R_EX_lanost_e_', 'R_EX_lanostest_SC_e_', 'R_EX_leu_L_e_', 'R_EX_lys_L_e_', 'R_EX_mal_L_e_', 'R_EX_malt_e_', 'R_EX_man_e_', 'R_EX_melib_e_', 'R_EX_met_L_e_', 'R_EX_nac_e_', 'R_EX_nadp_e_', 'R_EX_nh4_e_', 'R_EX_nmn_e_', 'R_EX_o2_e_', 'R_EX_oaa_e_', 'R_EX_ocdca_e_', 'R_EX_ocdcea_e_', 'R_EX_ocdcya_e_', 'R_EX_orn_e_', 'R_EX_pacald_e_', 'R_EX_pap_e_', 'R_EX_pc_SC_e_', 'R_EX_pectin_e_', 'R_EX_phe_L_e_', 'R_EX_pheac_e_', 'R_EX_pi_e_', 'R_EX_pnto_R_e_', 'R_EX_pro_L_e_', 'R_EX_ptd1ino_SC_e_', 'R_EX_ptrc_e_', 'R_EX_pyr_e_', 'R_EX_rib_D_e_', 'R_EX_ribflv_e_', 'R_EX_sbt_D_e_', 'R_EX_sbt_L_e_', 'R_EX_ser_L_e_', 'R_EX_so3_e_', 'R_EX_so4_e_', 'R_EX_spmd_e_', 'R_EX_sprm_e_', 'R_EX_srb_L_e_', 'R_EX_succ_e_', 'R_EX_sucr_e_', 'R_EX_thm_e_', 'R_EX_thmmp_e_', 'R_EX_thmpp_e_', 'R_EX_thr_L_e_', 'R_EX_thym_e_', 'R_EX_thymd_e_', 'R_EX_tre_e_', 'R_EX_trp_L_e_', 'R_EX_ttdca_e_', 'R_EX_tyr_L_e_', 'R_EX_ura_e_', 'R_EX_urea_e_', 'R_EX_uri_e_', 'R_EX_val_L_e_', 'R_EX_xan_e_', 'R_EX_xtsn_e_', 'R_EX_xyl_D_e_', 'R_EX_xylt_e_', 'R_EX_zymst_e_', 'R_EX_zymstest_SC_e_']

pickle.dump({'EXrxns': EXrxns}, open('iMM904_EXrxns.pkl', 'wb'))

