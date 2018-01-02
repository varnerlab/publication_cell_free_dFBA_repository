# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-03-20T10:48:32
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# Setup default flux bounds array -
	default_bounds_array = [
		0	100.0	;	# 1 M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c
		0	100.0	;	# 2 M_g6p_c --> M_f6p_c
		0	100.0	;	# 3 M_f6p_c --> M_g6p_c
		0	100.0	;	# 4 M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0	100.0	;	# 5 M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0	100.0	;	# 6 M_fdp_c --> M_dhap_c+M_g3p_c
		0	100.0	;	# 7 M_dhap_c+M_g3p_c --> M_fdp_c
		0	100.0	;	# 8 M_dhap_c --> M_g3p_c
		0	100.0	;	# 9 M_g3p_c --> M_dhap_c
		0	100.0	;	# 10 M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0	100.0	;	# 11 M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0	100.0	;	# 12 M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0	100.0	;	# 13 M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0	100.0	;	# 14 M_3pg_c --> M_2pg_c
		0	100.0	;	# 15 M_2pg_c --> M_3pg_c
		0	100.0	;	# 16 M_2pg_c --> M_h2o_c+M_pep_c
		0	100.0	;	# 17 M_h2o_c+M_pep_c --> M_2pg_c
		0	100.0	;	# 18 M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0	100.0	;	# 19 M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0	100.0	;	# 20 M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0	100.0	;	# 21 M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0	100.0	;	# 22 M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0	100.0	;	# 23 M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0	100.0	;	# 24 M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0	100.0	;	# 25 M_6pgl_c+M_h2o_c --> M_6pgc_c
		0	100.0	;	# 26 M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0	100.0	;	# 27 M_ru5p_D_c --> M_xu5p_D_c
		0	100.0	;	# 28 M_xu5p_D_c --> M_ru5p_D_c
		0	100.0	;	# 29 M_r5p_c --> M_ru5p_D_c
		0	100.0	;	# 30 M_ru5p_D_c --> M_r5p_c
		0	100.0	;	# 31 M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0	100.0	;	# 32 M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 33 M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 34 M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0	100.0	;	# 35 M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0	100.0	;	# 36 M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0	100.0	;	# 37 M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0	100.0	;	# 38 M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0	100.0	;	# 39 M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0	100.0	;	# 40 M_cit_c --> M_icit_c
		0	100.0	;	# 41 M_icit_c --> M_cit_c
		0	100.0	;	# 42 M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0	100.0	;	# 43 M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0	100.0	;	# 44 M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0	100.0	;	# 45 M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0	100.0	;	# 46 M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0	100.0	;	# 47 M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0	100.0	;	# 48 M_fum_c+M_h2o_c --> M_mal_L_c
		0	100.0	;	# 49 M_mal_L_c --> M_fum_c+M_h2o_c
		0	100.0	;	# 50 M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0	100.0	;	# 51 M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c
		0	100.0	;	# 52 2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_h_e
		0	100.0	;	# 53 4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_h_e
		0	100.0	;	# 54 2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_h_e
		0	100.0	;	# 55 M_adp_c+M_pi_c+4.0*M_h_e --> M_atp_c+4.0*M_h_c+M_h2o_c
		0	100.0	;	# 56 3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_h_e
		0	100.0	;	# 57 M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0	100.0	;	# 58 M_nadh_c+M_nadp_c+2.0*M_h_e --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0	100.0	;	# 59 M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0	100.0	;	# 60 M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0	100.0	;	# 61 M_ppi_c+M_h2o_c --> 2.0*M_pi_c
		0	100.0	;	# 62 M_icit_c --> M_glx_c+M_succ_c
		0	100.0	;	# 63 M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0	100.0	;	# 64 M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0	100.0	;	# 65 M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0	100.0	;	# 66 M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0	100.0	;	# 67 M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0	100.0	;	# 68 M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0	100.0	;	# 69 M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0	100.0	;	# 70 M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0	100.0	;	# 71 M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0	100.0	;	# 72 M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0	100.0	;	# 73 M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0	100.0	;	# 74 M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0	100.0	;	# 75 M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0	100.0	;	# 76 M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0	100.0	;	# 77 M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c
		0	100.0	;	# 78 M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c
		0	100.0	;	# 79 M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0	100.0	;	# 80 M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c
		0	100.0	;	# 81 M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c
		0	100.0	;	# 82 M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c
		0	100.0	;	# 83 M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0	100.0	;	# 84 M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0	100.0	;	# 85 M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c
		0	100.0	;	# 86 M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c
		0	100.0	;	# 87 M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0	100.0	;	# 88 M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c
		0	100.0	;	# 89 M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c
		0	100.0	;	# 90 2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0	100.0	;	# 91 M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0	100.0	;	# 92 M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 93 M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0	100.0	;	# 94 M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0	100.0	;	# 95 M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0	100.0	;	# 96 M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0	100.0	;	# 97 M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c
		0	100.0	;	# 98 M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c
		0	100.0	;	# 99 2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0	100.0	;	# 100 M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c
		0	100.0	;	# 101 M_asp_L_c --> M_fum_c+M_nh3_c
		0	100.0	;	# 102 M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c
		0	100.0	;	# 103 M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c
		0	100.0	;	# 104 M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c
		0	100.0	;	# 105 M_ser_L_c --> M_nh3_c+M_pyr_c
		0	100.0	;	# 106 M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0	100.0	;	# 107 M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0	100.0	;	# 108 M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c
		0	100.0	;	# 109 M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c
		0	100.0	;	# 110 M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 111 M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c
		0	100.0	;	# 112 M_lys_L_c --> M_co2_c+M_cadav_c
		0	100.0	;	# 113 M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c
		0	100.0	;	# 114 M_glu_L_c --> M_co2_c+M_gaba_c
		0	100.0	;	# 115 M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c
		0	100.0	;	# 116 M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c
		0	100.0	;	# 117 M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0	100.0	;	# 118 M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c
		0	100.0	;	# 119 M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0	100.0	;	# 120 M_4adochor_c --> M_4abz_c+M_pyr_c
		0	100.0	;	# 121 M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c
		0	100.0	;	# 122 M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0	100.0	;	# 123 M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0	100.0	;	# 124 M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c
		0	100.0	;	# 125 M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c
		0	100.0	;	# 126 M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c
		0	100.0	;	# 127 M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0	100.0	;	# 128 M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c
		0	100.0	;	# 129 M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c
		0	100.0	;	# 130 M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0	100.0	;	# 131 M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c
		0	100.0	;	# 132 M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c
		0	100.0	;	# 133 2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c
		0	100.0	;	# 134 M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c
		0	100.0	;	# 135 M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0	100.0	;	# 136 M_omp_c --> M_ump_c+M_co2_c
		0	100.0	;	# 137 M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c
		0	100.0	;	# 138 M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c
		0	100.0	;	# 139 M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0	100.0	;	# 140 M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0	100.0	;	# 141 M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0	100.0	;	# 142 M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0	100.0	;	# 143 M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0	100.0	;	# 144 M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0	100.0	;	# 145 M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0	100.0	;	# 146 M_saicar_c --> M_fum_c+M_aicar_c
		0	100.0	;	# 147 M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0	100.0	;	# 148 M_faicar_c --> M_imp_c+M_h2o_c
		0	100.0	;	# 149 M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0	100.0	;	# 150 M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0	100.0	;	# 151 M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0	100.0	;	# 152 M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c
		0	100.0	;	# 153 M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c
		0	100.0	;	# 154 M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c
		0	100.0	;	# 155 M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c
		0	100.0	;	# 156 M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0	100.0	;	# 157 M_utp_c+M_h2o_c --> M_udp_c+M_pi_c
		0	100.0	;	# 158 M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c
		0	100.0	;	# 159 M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c
		0	100.0	;	# 160 M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0	100.0	;	# 161 M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0	100.0	;	# 162 M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0	100.0	;	# 163 M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0	100.0	;	# 164 M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0	100.0	;	# 165 M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0	100.0	;	# 166 M_amp_c+M_atp_c --> 2.0*M_adp_c
		0	100.0	;	# 167 GENE_CAT+RNAP --> OPEN_GENE_CAT
		0	100.0	;	# 168 OPEN_GENE_CAT+151.0*M_gtp_c+144.0*M_ctp_c+189.0*M_utp_c+176.0*M_atp_c+660.0*M_h2o_c --> mRNA_CAT+GENE_CAT+RNAP+660.0*M_ppi_c
		0	100.0	;	# 169 mRNA_CAT --> 151.0*M_gmp_c+144.0*M_cmp_c+189.0*M_ump_c+176.0*M_amp_c
		0	100.0	;	# 170 mRNA_CAT+RIBOSOME --> RIBOSOME_START_CAT
		0	100.0	;	# 171 RIBOSOME_START_CAT+438.0*M_gtp_c+438.0*M_h2o_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA --> RIBOSOME+mRNA_CAT+PROTEIN_CAT+438.0*M_gdp_c+438.0*M_pi_c+219.0*tRNA
		0	100.0	;	# 172 15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA+15.0*M_h2o_c --> 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+15.0*M_ppi_c
		0	100.0	;	# 173 5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0	100.0	;	# 174 10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0	100.0	;	# 175 12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 176 5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0	100.0	;	# 177 12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 178 13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0	100.0	;	# 179 10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0	100.0	;	# 180 12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_his_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 181 9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
		0	100.0	;	# 182 13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0	100.0	;	# 183 12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0	100.0	;	# 184 9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_met_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
		0	100.0	;	# 185 20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA+20.0*M_h2o_c --> 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+20.0*M_ppi_c
		0	100.0	;	# 186 7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA+7.0*M_h2o_c --> 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+7.0*M_ppi_c
		0	100.0	;	# 187 10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0	100.0	;	# 188 13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0	100.0	;	# 189 5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0	100.0	;	# 190 11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+11.0*M_ppi_c
		0	100.0	;	# 191 16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16.0*M_h2o_c --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c
		0	100.0	;	# 192 tRNA --> []
		0	100.0	;	# 193 [] --> tRNA
		0	100.0	;	# 194 PROTEIN_CAT --> []
		0	100.0	;	# 195 [] --> M_o2_c
		0	100.0	;	# 196 M_co2_c --> []
		0	100.0	;	# 197 M_h_c --> []
		0	100.0	;	# 198 [] --> M_h_c
		0	100.0	;	# 199 [] --> M_h2s_c
		0	100.0	;	# 200 M_h2s_c --> []
		0	100.0	;	# 201 [] --> M_h2o_c
		0	100.0	;	# 202 M_h2o_c --> []
		0	100.0	;	# 203 [] --> M_pi_c
		0	100.0	;	# 204 M_pi_c --> []
		0	100.0	;	# 205 [] --> M_nh3_c
		0	100.0	;	# 206 M_nh3_c --> []
		0	100.0	;	# 207 [] --> M_glc_D_c
		0	100.0	;	# 208 [] --> M_hco3_c
		0	100.0	;	# 209 M_hco3_c --> []
		0	100.0	;	# 210 M_pyr_c --> []
		0	100.0	;	# 211 [] --> M_pyr_c
		0	100.0	;	# 212 M_ac_c --> []
		0	100.0	;	# 213 M_lac_D_c --> []
		0	100.0	;	# 214 M_succ_c --> []
		0	100.0	;	# 215 M_mal_L_c --> []
		0	100.0	;	# 216 M_fum_c --> []
		0	100.0	;	# 217 M_etoh_c --> []
		0	100.0	;	# 218 M_mglx_c --> []
		0	100.0	;	# 219 M_prop_c --> []
		0	100.0	;	# 220 M_indole_c --> []
		0	100.0	;	# 221 M_cadav_c --> []
		0	100.0	;	# 222 M_gaba_c --> []
		0	100.0	;	# 223 M_glycoA_c --> []
		0	100.0	;	# 224 [] --> M_ala_L_c
		0	100.0	;	# 225 M_ala_L_c --> []
		0	100.0	;	# 226 [] --> M_arg_L_c
		0	100.0	;	# 227 M_arg_L_c --> []
		0	100.0	;	# 228 [] --> M_asn_L_c
		0	100.0	;	# 229 M_asn_L_c --> []
		0	100.0	;	# 230 [] --> M_asp_L_c
		0	100.0	;	# 231 M_asp_L_c --> []
		0	100.0	;	# 232 [] --> M_cys_L_c
		0	100.0	;	# 233 M_cys_L_c --> []
		0	100.0	;	# 234 [] --> M_glu_L_c
		0	100.0	;	# 235 M_glu_L_c --> []
		0	100.0	;	# 236 [] --> M_gln_L_c
		0	100.0	;	# 237 M_gln_L_c --> []
		0	100.0	;	# 238 [] --> M_gly_L_c
		0	100.0	;	# 239 M_gly_L_c --> []
		0	100.0	;	# 240 [] --> M_his_L_c
		0	100.0	;	# 241 M_his_L_c --> []
		0	100.0	;	# 242 [] --> M_ile_L_c
		0	100.0	;	# 243 M_ile_L_c --> []
		0	100.0	;	# 244 [] --> M_leu_L_c
		0	100.0	;	# 245 M_leu_L_c --> []
		0	100.0	;	# 246 [] --> M_lys_L_c
		0	100.0	;	# 247 M_lys_L_c --> []
		0	100.0	;	# 248 [] --> M_met_L_c
		0	100.0	;	# 249 M_met_L_c --> []
		0	100.0	;	# 250 [] --> M_phe_L_c
		0	100.0	;	# 251 M_phe_L_c --> []
		0	100.0	;	# 252 [] --> M_pro_L_c
		0	100.0	;	# 253 M_pro_L_c --> []
		0	100.0	;	# 254 [] --> M_ser_L_c
		0	100.0	;	# 255 M_ser_L_c --> []
		0	100.0	;	# 256 [] --> M_thr_L_c
		0	100.0	;	# 257 M_thr_L_c --> []
		0	100.0	;	# 258 [] --> M_trp_L_c
		0	100.0	;	# 259 M_trp_L_c --> []
		0	100.0	;	# 260 [] --> M_tyr_L_c
		0	100.0	;	# 261 M_tyr_L_c --> []
		0	100.0	;	# 262 [] --> M_val_L_c
		0	100.0	;	# 263 M_val_L_c --> []
		0	100.0	;	# 264 M_atp_c --> []
		0	100.0	;	# 265 [] --> M_atp_c
		0	100.0	;	# 266 M_adp_c --> []
		0	100.0	;	# 267 [] --> M_adp_c
		0	100.0	;	# 268 M_amp_c --> []
		0	100.0	;	# 269 [] --> M_amp_c
		0	100.0	;	# 270 M_gtp_c --> []
		0	100.0	;	# 271 [] --> M_gtp_c
		0	100.0	;	# 272 M_gdp_c --> []
		0	100.0	;	# 273 [] --> M_gdp_c
		0	100.0	;	# 274 M_gmp_c --> []
		0	100.0	;	# 275 [] --> M_gmp_c
		0	100.0	;	# 276 M_ctp_c --> []
		0	100.0	;	# 277 [] --> M_ctp_c
		0	100.0	;	# 278 M_cdp_c --> []
		0	100.0	;	# 279 [] --> M_cdp_c
		0	100.0	;	# 280 M_cmp_c --> []
		0	100.0	;	# 281 [] --> M_cmp_c
		0	100.0	;	# 282 M_utp_c --> []
		0	100.0	;	# 283 [] --> M_utp_c
		0	100.0	;	# 284 M_udp_c --> []
		0	100.0	;	# 285 [] --> M_udp_c
		0	100.0	;	# 286 M_ump_c --> []
		0	100.0	;	# 287 [] --> M_ump_c
	];

	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1 GENE_CAT
		0.0	0.0	;	# 2 M_10fthf_c
		0.0	0.0	;	# 3 M_13dpg_c
		0.0	0.0	;	# 4 M_2ddg6p_c
		0.0	0.0	;	# 5 M_2pg_c
		0.0	0.0	;	# 6 M_3pg_c
		0.0	0.0	;	# 7 M_4abz_c
		0.0	0.0	;	# 8 M_4adochor_c
		0.0	0.0	;	# 9 M_5mthf_c
		0.0	0.0	;	# 10 M_5pbdra
		0.0	0.0	;	# 11 M_6pgc_c
		0.0	0.0	;	# 12 M_6pgl_c
		0.0	0.0	;	# 13 M_78dhf_c
		0.0	0.0	;	# 14 M_78mdp_c
		0.0	0.0	;	# 15 M_ac_c
		0.0	0.0	;	# 16 M_accoa_c
		0.0	0.0	;	# 17 M_actp_c
		0.0	0.0	;	# 18 M_adp_c
		0.0	0.0	;	# 19 M_aicar_c
		0.0	0.0	;	# 20 M_air_c
		0.0	0.0	;	# 21 M_akg_c
		0.0	0.0	;	# 22 M_ala_L_c
		0.0	0.0	;	# 23 M_ala_L_c_tRNA
		0.0	0.0	;	# 24 M_amp_c
		0.0	0.0	;	# 25 M_arg_L_c
		0.0	0.0	;	# 26 M_arg_L_c_tRNA
		0.0	0.0	;	# 27 M_asn_L_c
		0.0	0.0	;	# 28 M_asn_L_c_tRNA
		0.0	0.0	;	# 29 M_asp_L_c
		0.0	0.0	;	# 30 M_asp_L_c_tRNA
		0.0	0.0	;	# 31 M_atp_c
		0.0	0.0	;	# 32 M_cadav_c
		0.0	0.0	;	# 33 M_cair_c
		0.0	0.0	;	# 34 M_cdp_c
		0.0	0.0	;	# 35 M_chor_c
		0.0	0.0	;	# 36 M_cit_c
		0.0	0.0	;	# 37 M_clasp_c
		0.0	0.0	;	# 38 M_cmp_c
		0.0	0.0	;	# 39 M_co2_c
		0.0	0.0	;	# 40 M_coa_c
		0.0	0.0	;	# 41 M_ctp_c
		0.0	0.0	;	# 42 M_cys_L_c
		0.0	0.0	;	# 43 M_cys_L_c_tRNA
		0.0	0.0	;	# 44 M_dhap_c
		0.0	0.0	;	# 45 M_dhf_c
		0.0	0.0	;	# 46 M_e4p_c
		0.0	0.0	;	# 47 M_etoh_c
		0.0	0.0	;	# 48 M_f6p_c
		0.0	0.0	;	# 49 M_faicar_c
		0.0	0.0	;	# 50 M_fdp_c
		0.0	0.0	;	# 51 M_fgam_c
		0.0	0.0	;	# 52 M_fgar_c
		0.0	0.0	;	# 53 M_for_c
		0.0	0.0	;	# 54 M_fum_c
		0.0	0.0	;	# 55 M_g3p_c
		0.0	0.0	;	# 56 M_g6p_c
		0.0	0.0	;	# 57 M_gaba_c
		0.0	0.0	;	# 58 M_gar_c
		0.0	0.0	;	# 59 M_gdp_c
		0.0	0.0	;	# 60 M_glc_D_c
		0.0	0.0	;	# 61 M_gln_L_c
		0.0	0.0	;	# 62 M_gln_L_c_tRNA
		0.0	0.0	;	# 63 M_glu_L_c
		0.0	0.0	;	# 64 M_glu_L_c_tRNA
		0.0	0.0	;	# 65 M_glx_c
		0.0	0.0	;	# 66 M_gly_L_c
		0.0	0.0	;	# 67 M_gly_L_c_tRNA
		0.0	0.0	;	# 68 M_glycoA_c
		0.0	0.0	;	# 69 M_gmp_c
		0.0	0.0	;	# 70 M_gtp_c
		0.0	0.0	;	# 71 M_h2o2_c
		0.0	0.0	;	# 72 M_h2o_c
		0.0	0.0	;	# 73 M_h2s_c
		0.0	0.0	;	# 74 M_h_c
		0.0	0.0	;	# 75 M_hco3_c
		0.0	0.0	;	# 76 M_his_L_c
		0.0	0.0	;	# 77 M_his_L_c_tRNA
		0.0	0.0	;	# 78 M_icit_c
		0.0	0.0	;	# 79 M_ile_L_c
		0.0	0.0	;	# 80 M_ile_L_c_tRNA
		0.0	0.0	;	# 81 M_imp_c
		0.0	0.0	;	# 82 M_indole_c
		0.0	0.0	;	# 83 M_lac_D_c
		0.0	0.0	;	# 84 M_leu_L_c
		0.0	0.0	;	# 85 M_leu_L_c_tRNA
		0.0	0.0	;	# 86 M_lys_L_c
		0.0	0.0	;	# 87 M_lys_L_c_tRNA
		0.0	0.0	;	# 88 M_mal_L_c
		0.0	0.0	;	# 89 M_met_L_c
		0.0	0.0	;	# 90 M_met_L_c_tRNA
		0.0	0.0	;	# 91 M_methf_c
		0.0	0.0	;	# 92 M_mglx_c
		0.0	0.0	;	# 93 M_mlthf_c
		0.0	0.0	;	# 94 M_mql8_c
		0.0	0.0	;	# 95 M_mqn8_c
		0.0	0.0	;	# 96 M_nad_c
		0.0	0.0	;	# 97 M_nadh_c
		0.0	0.0	;	# 98 M_nadp_c
		0.0	0.0	;	# 99 M_nadph_c
		0.0	0.0	;	# 100 M_nh3_c
		0.0	0.0	;	# 101 M_o2_c
		0.0	0.0	;	# 102 M_oaa_c
		0.0	0.0	;	# 103 M_omp_c
		0.0	0.0	;	# 104 M_or_c
		0.0	0.0	;	# 105 M_pep_c
		0.0	0.0	;	# 106 M_phe_L_c
		0.0	0.0	;	# 107 M_phe_L_c_tRNA
		0.0	0.0	;	# 108 M_pi_c
		0.0	0.0	;	# 109 M_ppi_c
		0.0	0.0	;	# 110 M_pro_L_c
		0.0	0.0	;	# 111 M_pro_L_c_tRNA
		0.0	0.0	;	# 112 M_prop_c
		0.0	0.0	;	# 113 M_prpp_c
		0.0	0.0	;	# 114 M_pyr_c
		0.0	0.0	;	# 115 M_q8_c
		0.0	0.0	;	# 116 M_q8h2_c
		0.0	0.0	;	# 117 M_r5p_c
		0.0	0.0	;	# 118 M_ru5p_D_c
		0.0	0.0	;	# 119 M_s7p_c
		0.0	0.0	;	# 120 M_saicar_c
		0.0	0.0	;	# 121 M_ser_L_c
		0.0	0.0	;	# 122 M_ser_L_c_tRNA
		0.0	0.0	;	# 123 M_succ_c
		0.0	0.0	;	# 124 M_succoa_c
		0.0	0.0	;	# 125 M_thf_c
		0.0	0.0	;	# 126 M_thr_L_c
		0.0	0.0	;	# 127 M_thr_L_c_tRNA
		0.0	0.0	;	# 128 M_trp_L_c
		0.0	0.0	;	# 129 M_trp_L_c_tRNA
		0.0	0.0	;	# 130 M_tyr_L_c
		0.0	0.0	;	# 131 M_tyr_L_c_tRNA
		0.0	0.0	;	# 132 M_udp_c
		0.0	0.0	;	# 133 M_ump_c
		0.0	0.0	;	# 134 M_utp_c
		0.0	0.0	;	# 135 M_val_L_c
		0.0	0.0	;	# 136 M_val_L_c_tRNA
		0.0	0.0	;	# 137 M_xmp_c
		0.0	0.0	;	# 138 M_xu5p_D_c
		0.0	0.0	;	# 139 OPEN_GENE_CAT
		0.0	0.0	;	# 140 PROTEIN_CAT
		0.0	0.0	;	# 141 RIBOSOME
		0.0	0.0	;	# 142 RIBOSOME_START_CAT
		0.0	0.0	;	# 143 RNAP
		0.0	0.0	;	# 144 mRNA_CAT
		0.0	0.0	;	# 145 tRNA
		-0.0001	0.0	;	# 146 M_h_e
	];

	# Setup the objective coefficient array -
	objective_coefficient_array = [
		0.0	;	# 1 R_glk_atp::M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c
		0.0	;	# 2 R_pgi::M_g6p_c --> M_f6p_c
		0.0	;	# 3 R_pgi_reverse::M_f6p_c --> M_g6p_c
		0.0	;	# 4 R_pfk::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c
		0.0	;	# 5 R_fdp::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0.0	;	# 6 R_fbaA::M_fdp_c --> M_dhap_c+M_g3p_c
		0.0	;	# 7 R_fbaA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c
		0.0	;	# 8 R_tpiA::M_dhap_c --> M_g3p_c
		0.0	;	# 9 R_tpiA_reverse::M_g3p_c --> M_dhap_c
		0.0	;	# 10 R_gapA::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0.0	;	# 11 R_gapA_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0.0	;	# 12 R_pgk::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0.0	;	# 13 R_pgk_reverse::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0.0	;	# 14 R_gpm::M_3pg_c --> M_2pg_c
		0.0	;	# 15 R_gpm_reverse::M_2pg_c --> M_3pg_c
		0.0	;	# 16 R_eno::M_2pg_c --> M_h2o_c+M_pep_c
		0.0	;	# 17 R_eno_reverse::M_h2o_c+M_pep_c --> M_2pg_c
		0.0	;	# 18 R_pyk::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c
		0.0	;	# 19 R_pck::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0.0	;	# 20 R_ppc::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c
		0.0	;	# 21 R_pdh::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c
		0.0	;	# 22 R_pps::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c
		0.0	;	# 23 R_zwf::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0.0	;	# 24 R_zwf_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0.0	;	# 25 R_pgl::M_6pgl_c+M_h2o_c --> M_6pgc_c
		0.0	;	# 26 R_gnd::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c
		0.0	;	# 27 R_rpe::M_ru5p_D_c --> M_xu5p_D_c
		0.0	;	# 28 R_rpe_reverse::M_xu5p_D_c --> M_ru5p_D_c
		0.0	;	# 29 R_rpi::M_r5p_c --> M_ru5p_D_c
		0.0	;	# 30 R_rpi_reverse::M_ru5p_D_c --> M_r5p_c
		0.0	;	# 31 R_talAB::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0.0	;	# 32 R_talAB_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0.0	;	# 33 R_tkt1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0.0	;	# 34 R_tkt1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0.0	;	# 35 R_tkt2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0.0	;	# 36 R_tkt2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0.0	;	# 37 R_edd::M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0.0	;	# 38 R_eda::M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0.0	;	# 39 R_gltA::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c
		0.0	;	# 40 R_acn::M_cit_c --> M_icit_c
		0.0	;	# 41 R_acn_reverse::M_icit_c --> M_cit_c
		0.0	;	# 42 R_icd::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c
		0.0	;	# 43 R_icd_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c
		0.0	;	# 44 R_sucAB::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c
		0.0	;	# 45 R_sucCD::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0.0	;	# 46 R_sdh::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0.0	;	# 47 R_frd::M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0.0	;	# 48 R_fum::M_fum_c+M_h2o_c --> M_mal_L_c
		0.0	;	# 49 R_fum_reverse::M_mal_L_c --> M_fum_c+M_h2o_c
		0.0	;	# 50 R_mdh::M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0.0	;	# 51 R_mdh_reverse::M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c
		0.0	;	# 52 R_cyd::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_h_e
		0.0	;	# 53 R_cyo::4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_h_e
		0.0	;	# 54 R_app::2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_h_e
		0.0	;	# 55 R_atp::M_adp_c+M_pi_c+4.0*M_h_e --> M_atp_c+4.0*M_h_c+M_h2o_c
		0.0	;	# 56 R_nuo::3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_h_e
		0.0	;	# 57 R_pnt1::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0.0	;	# 58 R_pnt2::M_nadh_c+M_nadp_c+2.0*M_h_e --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0.0	;	# 59 R_ndh1::M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0.0	;	# 60 R_ndh2::M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0.0	;	# 61 R_ppa::M_ppi_c+M_h2o_c --> 2.0*M_pi_c
		0.0	;	# 62 R_aceA::M_icit_c --> M_glx_c+M_succ_c
		0.0	;	# 63 R_aceB::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c
		0.0	;	# 64 R_maeA::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c
		0.0	;	# 65 R_maeB::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c
		0.0	;	# 66 R_pta::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0.0	;	# 67 R_pta_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0.0	;	# 68 R_ackA::M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0.0	;	# 69 R_ackA_reverse::M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0.0	;	# 70 R_acs::M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0.0	;	# 71 R_adhE::M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0.0	;	# 72 R_adhE_reverse::M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0.0	;	# 73 R_ldh::M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0.0	;	# 74 R_ldh_reverse::M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0.0	;	# 75 R_pflAB::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0.0	;	# 76 R_alaAC::M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0.0	;	# 77 R_alaAC_reverse::M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c
		0.0	;	# 78 R_arg::M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c
		0.0	;	# 79 R_aspC::M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0.0	;	# 80 R_asnB::M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c
		0.0	;	# 81 R_asnA::M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c
		0.0	;	# 82 R_cysEMK::M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c
		0.0	;	# 83 R_gltBD::M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0.0	;	# 84 R_gdhA::M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0.0	;	# 85 R_gdhA_reverse::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c
		0.0	;	# 86 R_glnA::M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c
		0.0	;	# 87 R_glyA::M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0.0	;	# 88 R_his::M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c
		0.0	;	# 89 R_ile::M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c
		0.0	;	# 90 R_leu::2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0.0	;	# 91 R_lys::M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0.0	;	# 92 R_met::M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c
		0.0	;	# 93 R_phe::M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0.0	;	# 94 R_pro::M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0.0	;	# 95 R_serABC::M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0.0	;	# 96 R_thr::M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0.0	;	# 97 R_trp::M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c
		0.0	;	# 98 R_tyr::M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c
		0.0	;	# 99 R_val::2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0.0	;	# 100 R_arg_deg::M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c
		0.0	;	# 101 R_asp_deg::M_asp_L_c --> M_fum_c+M_nh3_c
		0.0	;	# 102 R_asn_deg::M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c
		0.0	;	# 103 R_gly_deg::M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c
		0.0	;	# 104 R_mglx_deg::M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c
		0.0	;	# 105 R_ser_deg::M_ser_L_c --> M_nh3_c+M_pyr_c
		0.0	;	# 106 R_pro_deg::M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0.0	;	# 107 R_thr_deg1::M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0.0	;	# 108 R_thr_deg2::M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c
		0.0	;	# 109 R_thr_deg3::M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c
		0.0	;	# 110 R_trp_deg::M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c
		0.0	;	# 111 R_cys_deg::M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c
		0.0	;	# 112 R_lys_deg::M_lys_L_c --> M_co2_c+M_cadav_c
		0.0	;	# 113 R_gln_deg::M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c
		0.0	;	# 114 R_glu_deg::M_glu_L_c --> M_co2_c+M_gaba_c
		0.0	;	# 115 R_gaba_deg1::M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c
		0.0	;	# 116 R_gaba_deg2::M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c
		0.0	;	# 117 R_chor::M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0.0	;	# 118 R_fol_e::M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c
		0.0	;	# 119 R_fol_1::M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0.0	;	# 120 R_fol_2a::M_4adochor_c --> M_4abz_c+M_pyr_c
		0.0	;	# 121 R_fol_2b::M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c
		0.0	;	# 122 R_fol_3::M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0.0	;	# 123 R_fol_4::M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0.0	;	# 124 R_gly_fol::M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c
		0.0	;	# 125 R_gly_fol_reverse::M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c
		0.0	;	# 126 R_mthfd::M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c
		0.0	;	# 127 R_mthfd_reverse::M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0.0	;	# 128 R_mthfc::M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c
		0.0	;	# 129 R_mthfc_reverse::M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c
		0.0	;	# 130 R_mthfr2a::M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0.0	;	# 131 R_mthfr2b::M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c
		0.0	;	# 132 R_prpp_syn::M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c
		0.0	;	# 133 R_or_syn_1::2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c
		0.0	;	# 134 R_or_syn_2::M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c
		0.0	;	# 135 R_omp_syn::M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0.0	;	# 136 R_ump_syn::M_omp_c --> M_ump_c+M_co2_c
		0.0	;	# 137 R_ctp_1::M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c
		0.0	;	# 138 R_ctp_2::M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c
		0.0	;	# 139 R_A_syn_1::M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0.0	;	# 140 R_A_syn_2::M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0.0	;	# 141 R_A_syn_3::M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0.0	;	# 142 R_A_syn_4::M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0.0	;	# 143 R_A_syn_5::M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0.0	;	# 144 R_A_syn_6::M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0.0	;	# 145 R_A_syn_7::M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0.0	;	# 146 R_A_syn_8::M_saicar_c --> M_fum_c+M_aicar_c
		0.0	;	# 147 R_A_syn_9::M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0.0	;	# 148 R_A_syn_10::M_faicar_c --> M_imp_c+M_h2o_c
		0.0	;	# 149 R_A_syn_12::M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0.0	;	# 150 R_xmp_syn::M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0.0	;	# 151 R_gmp_syn::M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0.0	;	# 152 R_atp_amp::M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c
		0.0	;	# 153 R_utp_ump::M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c
		0.0	;	# 154 R_ctp_cmp::M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c
		0.0	;	# 155 R_gtp_gmp::M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c
		0.0	;	# 156 R_atp_adp::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c
		0.0	;	# 157 R_utp_adp::M_utp_c+M_h2o_c --> M_udp_c+M_pi_c
		0.0	;	# 158 R_ctp_adp::M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c
		0.0	;	# 159 R_gtp_adp::M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c
		0.0	;	# 160 R_udp_utp::M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0.0	;	# 161 R_cdp_ctp::M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0.0	;	# 162 R_gdp_gtp::M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0.0	;	# 163 R_atp_ump::M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0.0	;	# 164 R_atp_cmp::M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0.0	;	# 165 R_atp_gmp::M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0.0	;	# 166 R_adk_atp::M_amp_c+M_atp_c --> 2.0*M_adp_c
		0.0	;	# 167 transcriptional_initiation_CAT::GENE_CAT+RNAP --> OPEN_GENE_CAT
		0.0	;	# 168 transcription_CAT::OPEN_GENE_CAT+151.0*M_gtp_c+144.0*M_ctp_c+189.0*M_utp_c+176.0*M_atp_c+660.0*M_h2o_c --> mRNA_CAT+GENE_CAT+RNAP+660.0*M_ppi_c
		0.0	;	# 169 mRNA_degradation_CAT::mRNA_CAT --> 151.0*M_gmp_c+144.0*M_cmp_c+189.0*M_ump_c+176.0*M_amp_c
		0.0	;	# 170 translation_initiation_CAT::mRNA_CAT+RIBOSOME --> RIBOSOME_START_CAT
		0.0	;	# 171 translation_CAT::RIBOSOME_START_CAT+438.0*M_gtp_c+438.0*M_h2o_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA --> RIBOSOME+mRNA_CAT+PROTEIN_CAT+438.0*M_gdp_c+438.0*M_pi_c+219.0*tRNA
		0.0	;	# 172 tRNA_charging_M_ala_L_c_CAT::15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA+15.0*M_h2o_c --> 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+15.0*M_ppi_c
		0.0	;	# 173 tRNA_charging_M_arg_L_c_CAT::5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0.0	;	# 174 tRNA_charging_M_asn_L_c_CAT::10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0.0	;	# 175 tRNA_charging_M_asp_L_c_CAT::12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 176 tRNA_charging_M_cys_L_c_CAT::5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0.0	;	# 177 tRNA_charging_M_glu_L_c_CAT::12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 178 tRNA_charging_M_gln_L_c_CAT::13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0.0	;	# 179 tRNA_charging_M_gly_L_c_CAT::10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0.0	;	# 180 tRNA_charging_M_his_L_c_CAT::12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_his_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 181 tRNA_charging_M_ile_L_c_CAT::9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
		0.0	;	# 182 tRNA_charging_M_leu_L_c_CAT::13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0.0	;	# 183 tRNA_charging_M_lys_L_c_CAT::12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c
		0.0	;	# 184 tRNA_charging_M_met_L_c_CAT::9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_met_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c
		0.0	;	# 185 tRNA_charging_M_phe_L_c_CAT::20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA+20.0*M_h2o_c --> 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+20.0*M_ppi_c
		0.0	;	# 186 tRNA_charging_M_pro_L_c_CAT::7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA+7.0*M_h2o_c --> 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+7.0*M_ppi_c
		0.0	;	# 187 tRNA_charging_M_ser_L_c_CAT::10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c
		0.0	;	# 188 tRNA_charging_M_thr_L_c_CAT::13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c
		0.0	;	# 189 tRNA_charging_M_trp_L_c_CAT::5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c
		0.0	;	# 190 tRNA_charging_M_tyr_L_c_CAT::11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+11.0*M_ppi_c
		0.0	;	# 191 tRNA_charging_M_val_L_c_CAT::16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16.0*M_h2o_c --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c
		0.0	;	# 192 tRNA_exchange::tRNA --> []
		0.0	;	# 193 tRNA_exchange_reverse::[] --> tRNA
		0.0	;	# 194 PROTEIN_export_CAT::PROTEIN_CAT --> []
		0.0	;	# 195 M_o2_c_exchange::[] --> M_o2_c
		0.0	;	# 196 M_co2_c_exchange::M_co2_c --> []
		0.0	;	# 197 M_h_c_exchange::M_h_c --> []
		0.0	;	# 198 M_h_c_exchange_reverse::[] --> M_h_c
		0.0	;	# 199 M_h2s_c_exchange::[] --> M_h2s_c
		0.0	;	# 200 M_h2s_c_exchange_reverse::M_h2s_c --> []
		0.0	;	# 201 M_h2o_c_exchange::[] --> M_h2o_c
		0.0	;	# 202 M_h2o_c_exchange_reverse::M_h2o_c --> []
		0.0	;	# 203 M_pi_c_exchange::[] --> M_pi_c
		0.0	;	# 204 M_pi_c_exchange_reverse::M_pi_c --> []
		0.0	;	# 205 M_nh3_c_exchange::[] --> M_nh3_c
		0.0	;	# 206 M_nh3_c_exchange_reverse::M_nh3_c --> []
		0.0	;	# 207 M_glc_D_c_exchange::[] --> M_glc_D_c
		0.0	;	# 208 M_hco3_c_exchange::[] --> M_hco3_c
		0.0	;	# 209 M_hco3_c_exchange_reverse::M_hco3_c --> []
		0.0	;	# 210 M_pyr_c_exchange::M_pyr_c --> []
		0.0	;	# 211 M_pyr_c_exchange_reverse::[] --> M_pyr_c
		0.0	;	# 212 M_ac_c_exchange::M_ac_c --> []
		0.0	;	# 213 M_lac_D_c_exchange::M_lac_D_c --> []
		0.0	;	# 214 M_succ_c_exchange::M_succ_c --> []
		0.0	;	# 215 M_mal_L_c_exchange::M_mal_L_c --> []
		0.0	;	# 216 M_fum_c_exchange::M_fum_c --> []
		0.0	;	# 217 M_etoh_c_exchange::M_etoh_c --> []
		0.0	;	# 218 M_mglx_c_exchange::M_mglx_c --> []
		0.0	;	# 219 M_prop_c_exchange::M_prop_c --> []
		0.0	;	# 220 M_indole_c_exchange::M_indole_c --> []
		0.0	;	# 221 M_cadav_c_exchange::M_cadav_c --> []
		0.0	;	# 222 M_gaba_c_exchange::M_gaba_c --> []
		0.0	;	# 223 M_glycoA_c_exchange::M_glycoA_c --> []
		0.0	;	# 224 M_ala_L_c_exchange::[] --> M_ala_L_c
		0.0	;	# 225 M_ala_L_c_exchange_reverse::M_ala_L_c --> []
		0.0	;	# 226 M_arg_L_c_exchange::[] --> M_arg_L_c
		0.0	;	# 227 M_arg_L_c_exchange_reverse::M_arg_L_c --> []
		0.0	;	# 228 M_asn_L_c_exchange::[] --> M_asn_L_c
		0.0	;	# 229 M_asn_L_c_exchange_reverse::M_asn_L_c --> []
		0.0	;	# 230 M_asp_L_c_exchange::[] --> M_asp_L_c
		0.0	;	# 231 M_asp_L_c_exchange_reverse::M_asp_L_c --> []
		0.0	;	# 232 M_cys_L_c_exchange::[] --> M_cys_L_c
		0.0	;	# 233 M_cys_L_c_exchange_reverse::M_cys_L_c --> []
		0.0	;	# 234 M_glu_L_c_exchange::[] --> M_glu_L_c
		0.0	;	# 235 M_glu_L_c_exchange_reverse::M_glu_L_c --> []
		0.0	;	# 236 M_gln_L_c_exchange::[] --> M_gln_L_c
		0.0	;	# 237 M_gln_L_c_exchange_reverse::M_gln_L_c --> []
		0.0	;	# 238 M_gly_L_c_exchange::[] --> M_gly_L_c
		0.0	;	# 239 M_gly_L_c_exchange_reverse::M_gly_L_c --> []
		0.0	;	# 240 M_his_L_c_exchange::[] --> M_his_L_c
		0.0	;	# 241 M_his_L_c_exchange_reverse::M_his_L_c --> []
		0.0	;	# 242 M_ile_L_c_exchange::[] --> M_ile_L_c
		0.0	;	# 243 M_ile_L_c_exchange_reverse::M_ile_L_c --> []
		0.0	;	# 244 M_leu_L_c_exchange::[] --> M_leu_L_c
		0.0	;	# 245 M_leu_L_c_exchange_reverse::M_leu_L_c --> []
		0.0	;	# 246 M_lys_L_c_exchange::[] --> M_lys_L_c
		0.0	;	# 247 M_lys_L_c_exchange_reverse::M_lys_L_c --> []
		0.0	;	# 248 M_met_L_c_exchange::[] --> M_met_L_c
		0.0	;	# 249 M_met_L_c_exchange_reverse::M_met_L_c --> []
		0.0	;	# 250 M_phe_L_c_exchange::[] --> M_phe_L_c
		0.0	;	# 251 M_phe_L_c_exchange_reverse::M_phe_L_c --> []
		0.0	;	# 252 M_pro_L_c_exchange::[] --> M_pro_L_c
		0.0	;	# 253 M_pro_L_c_exchange_reverse::M_pro_L_c --> []
		0.0	;	# 254 M_ser_L_c_exchange::[] --> M_ser_L_c
		0.0	;	# 255 M_ser_L_c_exchange_reverse::M_ser_L_c --> []
		0.0	;	# 256 M_thr_L_c_exchange::[] --> M_thr_L_c
		0.0	;	# 257 M_thr_L_c_exchange_reverse::M_thr_L_c --> []
		0.0	;	# 258 M_trp_L_c_exchange::[] --> M_trp_L_c
		0.0	;	# 259 M_trp_L_c_exchange_reverse::M_trp_L_c --> []
		0.0	;	# 260 M_tyr_L_c_exchange::[] --> M_tyr_L_c
		0.0	;	# 261 M_tyr_L_c_exchange_reverse::M_tyr_L_c --> []
		0.0	;	# 262 M_val_L_c_exchange::[] --> M_val_L_c
		0.0	;	# 263 M_val_L_c_exchange_reverse::M_val_L_c --> []
		0.0	;	# 264 M_atp_exchange::M_atp_c --> []
		0.0	;	# 265 M_atp_exchange_reverse::[] --> M_atp_c
		0.0	;	# 266 M_adp_exchange::M_adp_c --> []
		0.0	;	# 267 M_adp_exchange_reverse::[] --> M_adp_c
		0.0	;	# 268 M_amp_exchange::M_amp_c --> []
		0.0	;	# 269 M_amp_exchange_reverse::[] --> M_amp_c
		0.0	;	# 270 M_gtp_exchange::M_gtp_c --> []
		0.0	;	# 271 M_gtp_exchange_reverse::[] --> M_gtp_c
		0.0	;	# 272 M_gdp_exchange::M_gdp_c --> []
		0.0	;	# 273 M_gdp_exchange_reverse::[] --> M_gdp_c
		0.0	;	# 274 M_gmp_exchange::M_gmp_c --> []
		0.0	;	# 275 M_gmp_exchange_reverse::[] --> M_gmp_c
		0.0	;	# 276 M_ctp_exchange::M_ctp_c --> []
		0.0	;	# 277 M_ctp_exchange_reverse::[] --> M_ctp_c
		0.0	;	# 278 M_cdp_exchange::M_cdp_c --> []
		0.0	;	# 279 M_cdp_exchange_reverse::[] --> M_cdp_c
		0.0	;	# 280 M_cmp_exchange::M_cmp_c --> []
		0.0	;	# 281 M_cmp_exchange_reverse::[] --> M_cmp_c
		0.0	;	# 282 M_utp_exchange::M_utp_c --> []
		0.0	;	# 283 M_utp_exchange_reverse::[] --> M_utp_c
		0.0	;	# 284 M_udp_exchange::M_udp_c --> []
		0.0	;	# 285 M_udp_exchange_reverse::[] --> M_udp_c
		0.0	;	# 286 M_ump_exchange::M_ump_c --> []
		0.0	;	# 287 M_ump_exchange_reverse::[] --> M_ump_c
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"R_glk_atp::M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c"
		"R_pgi::M_g6p_c --> M_f6p_c"
		"R_pgi_reverse::M_f6p_c --> M_g6p_c"
		"R_pfk::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c"
		"R_fdp::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c"
		"R_fbaA::M_fdp_c --> M_dhap_c+M_g3p_c"
		"R_fbaA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c"
		"R_tpiA::M_dhap_c --> M_g3p_c"
		"R_tpiA_reverse::M_g3p_c --> M_dhap_c"
		"R_gapA::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c"
		"R_gapA_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c"
		"R_pgk::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c"
		"R_pgk_reverse::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c"
		"R_gpm::M_3pg_c --> M_2pg_c"
		"R_gpm_reverse::M_2pg_c --> M_3pg_c"
		"R_eno::M_2pg_c --> M_h2o_c+M_pep_c"
		"R_eno_reverse::M_h2o_c+M_pep_c --> M_2pg_c"
		"R_pyk::M_adp_c+M_pep_c --> M_atp_c+M_pyr_c"
		"R_pck::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c"
		"R_ppc::M_co2_c+M_h2o_c+M_pep_c --> M_oaa_c+M_pi_c"
		"R_pdh::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c+M_h_c"
		"R_pps::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+M_pep_c+M_pi_c"
		"R_zwf::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c"
		"R_zwf_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c"
		"R_pgl::M_6pgl_c+M_h2o_c --> M_6pgc_c"
		"R_gnd::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c+M_h_c"
		"R_rpe::M_ru5p_D_c --> M_xu5p_D_c"
		"R_rpe_reverse::M_xu5p_D_c --> M_ru5p_D_c"
		"R_rpi::M_r5p_c --> M_ru5p_D_c"
		"R_rpi_reverse::M_ru5p_D_c --> M_r5p_c"
		"R_talAB::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c"
		"R_talAB_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c"
		"R_tkt1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c"
		"R_tkt1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c"
		"R_tkt2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c"
		"R_tkt2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c"
		"R_edd::M_6pgc_c --> M_2ddg6p_c+M_h2o_c"
		"R_eda::M_2ddg6p_c --> M_g3p_c+M_pyr_c"
		"R_gltA::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c"
		"R_acn::M_cit_c --> M_icit_c"
		"R_acn_reverse::M_icit_c --> M_cit_c"
		"R_icd::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c+M_h_c"
		"R_icd_reverse::M_akg_c+M_co2_c+M_nadph_c+M_h_c --> M_icit_c+M_nadp_c"
		"R_sucAB::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c+M_h_c"
		"R_sucCD::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c"
		"R_sdh::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c"
		"R_frd::M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c"
		"R_fum::M_fum_c+M_h2o_c --> M_mal_L_c"
		"R_fum_reverse::M_mal_L_c --> M_fum_c+M_h2o_c"
		"R_mdh::M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c"
		"R_mdh_reverse::M_oaa_c+M_h_c+M_nadh_c --> M_mal_L_c+M_nad_c"
		"R_cyd::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_h_e"
		"R_cyo::4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_h_e"
		"R_app::2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_h_e"
		"R_atp::M_adp_c+M_pi_c+4.0*M_h_e --> M_atp_c+4.0*M_h_c+M_h2o_c"
		"R_nuo::3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_h_e"
		"R_pnt1::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c"
		"R_pnt2::M_nadh_c+M_nadp_c+2.0*M_h_e --> 2.0*M_h_c+M_nad_c+M_nadph_c"
		"R_ndh1::M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c"
		"R_ndh2::M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c"
		"R_ppa::M_ppi_c+M_h2o_c --> 2.0*M_pi_c"
		"R_aceA::M_icit_c --> M_glx_c+M_succ_c"
		"R_aceB::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_mal_L_c"
		"R_maeA::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c+M_h_c"
		"R_maeB::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c+M_h_c"
		"R_pta::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c"
		"R_pta_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c"
		"R_ackA::M_actp_c+M_adp_c --> M_ac_c+M_atp_c"
		"R_ackA_reverse::M_ac_c+M_atp_c --> M_actp_c+M_adp_c"
		"R_acs::M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c"
		"R_adhE::M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c"
		"R_adhE_reverse::M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c"
		"R_ldh::M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c"
		"R_ldh_reverse::M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c"
		"R_pflAB::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c"
		"R_alaAC::M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c"
		"R_alaAC_reverse::M_ala_L_c+M_akg_c --> M_pyr_c+M_glu_L_c"
		"R_arg::M_accoa_c+2.0*M_glu_L_c+3.0*M_atp_c+M_nadph_c+M_h_c+M_h2o_c+M_nh3_c+M_co2_c+M_asp_L_c --> M_coa_c+2.0*M_adp_c+2.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c+M_arg_L_c"
		"R_aspC::M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c"
		"R_asnB::M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_ppi_c+M_amp_c"
		"R_asnA::M_asp_L_c+M_atp_c+M_nh3_c --> M_asn_L_c+M_ppi_c+M_amp_c"
		"R_cysEMK::M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_ac_c"
		"R_gltBD::M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c"
		"R_gdhA::M_akg_c+M_nadph_c+M_nh3_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c"
		"R_gdhA_reverse::M_glu_L_c+M_h2o_c+M_nadp_c --> M_akg_c+M_nadph_c+M_nh3_c+M_h_c"
		"R_glnA::M_glu_L_c+M_atp_c+M_nh3_c --> M_gln_L_c+M_adp_c+M_pi_c"
		"R_glyA::M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c"
		"R_his::M_gln_L_c+M_r5p_c+2.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_nadh_c+M_amp_c+M_pi_c+2.0*M_ppi_c+2.0*M_h_c"
		"R_ile::M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh3_c+M_co2_c+M_nadp_c+M_akg_c"
		"R_leu::2.0*M_pyr_c+M_glu_L_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c"
		"R_lys::M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c"
		"R_met::M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c+M_h2o_c+2.0*M_h_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh3_c+M_pyr_c"
		"R_phe::M_chor_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c"
		"R_pro::M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c"
		"R_serABC::M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c"
		"R_thr::M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c"
		"R_trp::M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+M_amp_c"
		"R_tyr::M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c+M_h_c"
		"R_val::2.0*M_pyr_c+M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c"
		"R_arg_deg::M_arg_L_c+4.0*M_h2o_c+M_nad_c+M_akg_c+M_succoa_c --> M_h_c+M_co2_c+2.0*M_glu_L_c+2.0*M_nh3_c+M_nadh_c+M_succ_c+M_coa_c"
		"R_asp_deg::M_asp_L_c --> M_fum_c+M_nh3_c"
		"R_asn_deg::M_asn_L_c+M_amp_c+M_ppi_c --> M_nh3_c+M_asp_L_c+M_atp_c"
		"R_gly_deg::M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c"
		"R_mglx_deg::M_mglx_c+M_nad_c+M_h2o_c --> M_pyr_c+M_nadh_c+M_h_c"
		"R_ser_deg::M_ser_L_c --> M_nh3_c+M_pyr_c"
		"R_pro_deg::M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c"
		"R_thr_deg1::M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c"
		"R_thr_deg2::M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh3_c+M_mglx_c+M_h_c"
		"R_thr_deg3::M_thr_L_c+M_pi_c+M_adp_c --> M_nh3_c+M_for_c+M_atp_c+M_prop_c"
		"R_trp_deg::M_trp_L_c+M_h2o_c --> M_indole_c+M_nh3_c+M_pyr_c"
		"R_cys_deg::M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh3_c+M_pyr_c"
		"R_lys_deg::M_lys_L_c --> M_co2_c+M_cadav_c"
		"R_gln_deg::M_gln_L_c+M_h2o_c --> M_nh3_c+M_glu_L_c"
		"R_glu_deg::M_glu_L_c --> M_co2_c+M_gaba_c"
		"R_gaba_deg1::M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadh_c"
		"R_gaba_deg2::M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+M_h_c+M_nadph_c"
		"R_chor::M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c+M_h_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c"
		"R_fol_e::M_gtp_c+4.0*M_h2o_c --> M_for_c+3.0*M_pi_c+M_glycoA_c+M_78mdp_c"
		"R_fol_1::M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c"
		"R_fol_2a::M_4adochor_c --> M_4abz_c+M_pyr_c"
		"R_fol_2b::M_4abz_c+M_78mdp_c --> M_78dhf_c+M_h2o_c"
		"R_fol_3::M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c"
		"R_fol_4::M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c"
		"R_gly_fol::M_gly_L_c+M_thf_c+M_nad_c --> M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c"
		"R_gly_fol_reverse::M_mlthf_c+M_nh3_c+M_co2_c+M_nadh_c+M_h_c --> M_gly_L_c+M_thf_c+M_nad_c"
		"R_mthfd::M_mlthf_c+M_nadp_c --> M_methf_c+M_nadph_c"
		"R_mthfd_reverse::M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c"
		"R_mthfc::M_h2o_c+M_methf_c --> M_10fthf_c+M_h_c"
		"R_mthfc_reverse::M_10fthf_c+M_h_c --> M_h2o_c+M_methf_c"
		"R_mthfr2a::M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c"
		"R_mthfr2b::M_mlthf_c+M_h_c+M_nadph_c --> M_5mthf_c+M_nadp_c"
		"R_prpp_syn::M_r5p_c+M_atp_c --> M_prpp_c+M_amp_c"
		"R_or_syn_1::2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_adp_c+M_glu_L_c+M_pi_c+M_clasp_c+2.0*M_h_c"
		"R_or_syn_2::M_clasp_c+M_asp_L_c+M_q8_c --> M_or_c+M_q8h2_c+M_h2o_c+M_pi_c"
		"R_omp_syn::M_prpp_c+M_or_c --> M_omp_c+M_ppi_c"
		"R_ump_syn::M_omp_c --> M_ump_c+M_co2_c"
		"R_ctp_1::M_utp_c+M_atp_c+M_nh3_c --> M_ctp_c+M_adp_c+M_pi_c"
		"R_ctp_2::M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c"
		"R_A_syn_1::M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c"
		"R_A_syn_2::M_atp_c+M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c"
		"R_A_syn_3::M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c"
		"R_A_syn_4::M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c"
		"R_A_syn_5::M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c"
		"R_A_syn_6::M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c"
		"R_A_syn_7::M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c"
		"R_A_syn_8::M_saicar_c --> M_fum_c+M_aicar_c"
		"R_A_syn_9::M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c"
		"R_A_syn_10::M_faicar_c --> M_imp_c+M_h2o_c"
		"R_A_syn_12::M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c"
		"R_xmp_syn::M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c"
		"R_gmp_syn::M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c"
		"R_atp_amp::M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c"
		"R_utp_ump::M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c"
		"R_ctp_cmp::M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c"
		"R_gtp_gmp::M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c"
		"R_atp_adp::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c"
		"R_utp_adp::M_utp_c+M_h2o_c --> M_udp_c+M_pi_c"
		"R_ctp_adp::M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c"
		"R_gtp_adp::M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c"
		"R_udp_utp::M_udp_c+M_atp_c --> M_utp_c+M_adp_c"
		"R_cdp_ctp::M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c"
		"R_gdp_gtp::M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c"
		"R_atp_ump::M_atp_c+M_ump_c --> M_adp_c+M_udp_c"
		"R_atp_cmp::M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c"
		"R_atp_gmp::M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c"
		"R_adk_atp::M_amp_c+M_atp_c --> 2.0*M_adp_c"
		"transcriptional_initiation_CAT::GENE_CAT+RNAP --> OPEN_GENE_CAT"
		"transcription_CAT::OPEN_GENE_CAT+151.0*M_gtp_c+144.0*M_ctp_c+189.0*M_utp_c+176.0*M_atp_c+660.0*M_h2o_c --> mRNA_CAT+GENE_CAT+RNAP+660.0*M_ppi_c"
		"mRNA_degradation_CAT::mRNA_CAT --> 151.0*M_gmp_c+144.0*M_cmp_c+189.0*M_ump_c+176.0*M_amp_c"
		"translation_initiation_CAT::mRNA_CAT+RIBOSOME --> RIBOSOME_START_CAT"
		"translation_CAT::RIBOSOME_START_CAT+438.0*M_gtp_c+438.0*M_h2o_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA --> RIBOSOME+mRNA_CAT+PROTEIN_CAT+438.0*M_gdp_c+438.0*M_pi_c+219.0*tRNA"
		"tRNA_charging_M_ala_L_c_CAT::15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA+15.0*M_h2o_c --> 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+15.0*M_ppi_c"
		"tRNA_charging_M_arg_L_c_CAT::5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c"
		"tRNA_charging_M_asn_L_c_CAT::10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c"
		"tRNA_charging_M_asp_L_c_CAT::12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_charging_M_cys_L_c_CAT::5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c"
		"tRNA_charging_M_glu_L_c_CAT::12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_charging_M_gln_L_c_CAT::13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c"
		"tRNA_charging_M_gly_L_c_CAT::10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c"
		"tRNA_charging_M_his_L_c_CAT::12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_his_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_charging_M_ile_L_c_CAT::9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c"
		"tRNA_charging_M_leu_L_c_CAT::13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c"
		"tRNA_charging_M_lys_L_c_CAT::12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA+12.0*M_h2o_c --> 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+12.0*M_ppi_c"
		"tRNA_charging_M_met_L_c_CAT::9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA+9.0*M_h2o_c --> 9.0*M_met_L_c_tRNA+9.0*M_amp_c+9.0*M_ppi_c"
		"tRNA_charging_M_phe_L_c_CAT::20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA+20.0*M_h2o_c --> 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+20.0*M_ppi_c"
		"tRNA_charging_M_pro_L_c_CAT::7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA+7.0*M_h2o_c --> 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+7.0*M_ppi_c"
		"tRNA_charging_M_ser_L_c_CAT::10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA+10.0*M_h2o_c --> 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+10.0*M_ppi_c"
		"tRNA_charging_M_thr_L_c_CAT::13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA+13.0*M_h2o_c --> 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+13.0*M_ppi_c"
		"tRNA_charging_M_trp_L_c_CAT::5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA+5.0*M_h2o_c --> 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+5.0*M_ppi_c"
		"tRNA_charging_M_tyr_L_c_CAT::11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA+11.0*M_h2o_c --> 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+11.0*M_ppi_c"
		"tRNA_charging_M_val_L_c_CAT::16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16.0*M_h2o_c --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c"
		"tRNA_exchange::tRNA --> []"
		"tRNA_exchange_reverse::[] --> tRNA"
		"PROTEIN_export_CAT::PROTEIN_CAT --> []"
		"M_o2_c_exchange::[] --> M_o2_c"
		"M_co2_c_exchange::M_co2_c --> []"
		"M_h_c_exchange::M_h_c --> []"
		"M_h_c_exchange_reverse::[] --> M_h_c"
		"M_h2s_c_exchange::[] --> M_h2s_c"
		"M_h2s_c_exchange_reverse::M_h2s_c --> []"
		"M_h2o_c_exchange::[] --> M_h2o_c"
		"M_h2o_c_exchange_reverse::M_h2o_c --> []"
		"M_pi_c_exchange::[] --> M_pi_c"
		"M_pi_c_exchange_reverse::M_pi_c --> []"
		"M_nh3_c_exchange::[] --> M_nh3_c"
		"M_nh3_c_exchange_reverse::M_nh3_c --> []"
		"M_glc_D_c_exchange::[] --> M_glc_D_c"
		"M_hco3_c_exchange::[] --> M_hco3_c"
		"M_hco3_c_exchange_reverse::M_hco3_c --> []"
		"M_pyr_c_exchange::M_pyr_c --> []"
		"M_pyr_c_exchange_reverse::[] --> M_pyr_c"
		"M_ac_c_exchange::M_ac_c --> []"
		"M_lac_D_c_exchange::M_lac_D_c --> []"
		"M_succ_c_exchange::M_succ_c --> []"
		"M_mal_L_c_exchange::M_mal_L_c --> []"
		"M_fum_c_exchange::M_fum_c --> []"
		"M_etoh_c_exchange::M_etoh_c --> []"
		"M_mglx_c_exchange::M_mglx_c --> []"
		"M_prop_c_exchange::M_prop_c --> []"
		"M_indole_c_exchange::M_indole_c --> []"
		"M_cadav_c_exchange::M_cadav_c --> []"
		"M_gaba_c_exchange::M_gaba_c --> []"
		"M_glycoA_c_exchange::M_glycoA_c --> []"
		"M_ala_L_c_exchange::[] --> M_ala_L_c"
		"M_ala_L_c_exchange_reverse::M_ala_L_c --> []"
		"M_arg_L_c_exchange::[] --> M_arg_L_c"
		"M_arg_L_c_exchange_reverse::M_arg_L_c --> []"
		"M_asn_L_c_exchange::[] --> M_asn_L_c"
		"M_asn_L_c_exchange_reverse::M_asn_L_c --> []"
		"M_asp_L_c_exchange::[] --> M_asp_L_c"
		"M_asp_L_c_exchange_reverse::M_asp_L_c --> []"
		"M_cys_L_c_exchange::[] --> M_cys_L_c"
		"M_cys_L_c_exchange_reverse::M_cys_L_c --> []"
		"M_glu_L_c_exchange::[] --> M_glu_L_c"
		"M_glu_L_c_exchange_reverse::M_glu_L_c --> []"
		"M_gln_L_c_exchange::[] --> M_gln_L_c"
		"M_gln_L_c_exchange_reverse::M_gln_L_c --> []"
		"M_gly_L_c_exchange::[] --> M_gly_L_c"
		"M_gly_L_c_exchange_reverse::M_gly_L_c --> []"
		"M_his_L_c_exchange::[] --> M_his_L_c"
		"M_his_L_c_exchange_reverse::M_his_L_c --> []"
		"M_ile_L_c_exchange::[] --> M_ile_L_c"
		"M_ile_L_c_exchange_reverse::M_ile_L_c --> []"
		"M_leu_L_c_exchange::[] --> M_leu_L_c"
		"M_leu_L_c_exchange_reverse::M_leu_L_c --> []"
		"M_lys_L_c_exchange::[] --> M_lys_L_c"
		"M_lys_L_c_exchange_reverse::M_lys_L_c --> []"
		"M_met_L_c_exchange::[] --> M_met_L_c"
		"M_met_L_c_exchange_reverse::M_met_L_c --> []"
		"M_phe_L_c_exchange::[] --> M_phe_L_c"
		"M_phe_L_c_exchange_reverse::M_phe_L_c --> []"
		"M_pro_L_c_exchange::[] --> M_pro_L_c"
		"M_pro_L_c_exchange_reverse::M_pro_L_c --> []"
		"M_ser_L_c_exchange::[] --> M_ser_L_c"
		"M_ser_L_c_exchange_reverse::M_ser_L_c --> []"
		"M_thr_L_c_exchange::[] --> M_thr_L_c"
		"M_thr_L_c_exchange_reverse::M_thr_L_c --> []"
		"M_trp_L_c_exchange::[] --> M_trp_L_c"
		"M_trp_L_c_exchange_reverse::M_trp_L_c --> []"
		"M_tyr_L_c_exchange::[] --> M_tyr_L_c"
		"M_tyr_L_c_exchange_reverse::M_tyr_L_c --> []"
		"M_val_L_c_exchange::[] --> M_val_L_c"
		"M_val_L_c_exchange_reverse::M_val_L_c --> []"
		"M_atp_exchange::M_atp_c --> []"
		"M_atp_exchange_reverse::[] --> M_atp_c"
		"M_adp_exchange::M_adp_c --> []"
		"M_adp_exchange_reverse::[] --> M_adp_c"
		"M_amp_exchange::M_amp_c --> []"
		"M_amp_exchange_reverse::[] --> M_amp_c"
		"M_gtp_exchange::M_gtp_c --> []"
		"M_gtp_exchange_reverse::[] --> M_gtp_c"
		"M_gdp_exchange::M_gdp_c --> []"
		"M_gdp_exchange_reverse::[] --> M_gdp_c"
		"M_gmp_exchange::M_gmp_c --> []"
		"M_gmp_exchange_reverse::[] --> M_gmp_c"
		"M_ctp_exchange::M_ctp_c --> []"
		"M_ctp_exchange_reverse::[] --> M_ctp_c"
		"M_cdp_exchange::M_cdp_c --> []"
		"M_cdp_exchange_reverse::[] --> M_cdp_c"
		"M_cmp_exchange::M_cmp_c --> []"
		"M_cmp_exchange_reverse::[] --> M_cmp_c"
		"M_utp_exchange::M_utp_c --> []"
		"M_utp_exchange_reverse::[] --> M_utp_c"
		"M_udp_exchange::M_udp_c --> []"
		"M_udp_exchange_reverse::[] --> M_udp_c"
		"M_ump_exchange::M_ump_c --> []"
		"M_ump_exchange_reverse::[] --> M_ump_c"
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"GENE_CAT"
		"M_10fthf_c"
		"M_13dpg_c"
		"M_2ddg6p_c"
		"M_2pg_c"
		"M_3pg_c"
		"M_4abz_c"
		"M_4adochor_c"
		"M_5mthf_c"
		"M_5pbdra"
		"M_6pgc_c"
		"M_6pgl_c"
		"M_78dhf_c"
		"M_78mdp_c"
		"M_ac_c"
		"M_accoa_c"
		"M_actp_c"
		"M_adp_c"
		"M_aicar_c"
		"M_air_c"
		"M_akg_c"
		"M_ala_L_c"
		"M_ala_L_c_tRNA"
		"M_amp_c"
		"M_arg_L_c"
		"M_arg_L_c_tRNA"
		"M_asn_L_c"
		"M_asn_L_c_tRNA"
		"M_asp_L_c"
		"M_asp_L_c_tRNA"
		"M_atp_c"
		"M_cadav_c"
		"M_cair_c"
		"M_cdp_c"
		"M_chor_c"
		"M_cit_c"
		"M_clasp_c"
		"M_cmp_c"
		"M_co2_c"
		"M_coa_c"
		"M_ctp_c"
		"M_cys_L_c"
		"M_cys_L_c_tRNA"
		"M_dhap_c"
		"M_dhf_c"
		"M_e4p_c"
		"M_etoh_c"
		"M_f6p_c"
		"M_faicar_c"
		"M_fdp_c"
		"M_fgam_c"
		"M_fgar_c"
		"M_for_c"
		"M_fum_c"
		"M_g3p_c"
		"M_g6p_c"
		"M_gaba_c"
		"M_gar_c"
		"M_gdp_c"
		"M_glc_D_c"
		"M_gln_L_c"
		"M_gln_L_c_tRNA"
		"M_glu_L_c"
		"M_glu_L_c_tRNA"
		"M_glx_c"
		"M_gly_L_c"
		"M_gly_L_c_tRNA"
		"M_glycoA_c"
		"M_gmp_c"
		"M_gtp_c"
		"M_h2o2_c"
		"M_h2o_c"
		"M_h2s_c"
		"M_h_c"
		"M_hco3_c"
		"M_his_L_c"
		"M_his_L_c_tRNA"
		"M_icit_c"
		"M_ile_L_c"
		"M_ile_L_c_tRNA"
		"M_imp_c"
		"M_indole_c"
		"M_lac_D_c"
		"M_leu_L_c"
		"M_leu_L_c_tRNA"
		"M_lys_L_c"
		"M_lys_L_c_tRNA"
		"M_mal_L_c"
		"M_met_L_c"
		"M_met_L_c_tRNA"
		"M_methf_c"
		"M_mglx_c"
		"M_mlthf_c"
		"M_mql8_c"
		"M_mqn8_c"
		"M_nad_c"
		"M_nadh_c"
		"M_nadp_c"
		"M_nadph_c"
		"M_nh3_c"
		"M_o2_c"
		"M_oaa_c"
		"M_omp_c"
		"M_or_c"
		"M_pep_c"
		"M_phe_L_c"
		"M_phe_L_c_tRNA"
		"M_pi_c"
		"M_ppi_c"
		"M_pro_L_c"
		"M_pro_L_c_tRNA"
		"M_prop_c"
		"M_prpp_c"
		"M_pyr_c"
		"M_q8_c"
		"M_q8h2_c"
		"M_r5p_c"
		"M_ru5p_D_c"
		"M_s7p_c"
		"M_saicar_c"
		"M_ser_L_c"
		"M_ser_L_c_tRNA"
		"M_succ_c"
		"M_succoa_c"
		"M_thf_c"
		"M_thr_L_c"
		"M_thr_L_c_tRNA"
		"M_trp_L_c"
		"M_trp_L_c_tRNA"
		"M_tyr_L_c"
		"M_tyr_L_c_tRNA"
		"M_udp_c"
		"M_ump_c"
		"M_utp_c"
		"M_val_L_c"
		"M_val_L_c_tRNA"
		"M_xmp_c"
		"M_xu5p_D_c"
		"OPEN_GENE_CAT"
		"PROTEIN_CAT"
		"RIBOSOME"
		"RIBOSOME_START_CAT"
		"RNAP"
		"mRNA_CAT"
		"tRNA"
		"M_h_e"
	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
