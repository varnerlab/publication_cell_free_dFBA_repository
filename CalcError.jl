using Interpolations
function CalcError(Upper,Lower,Mean,experimental_time,time_state_array,data_dictionary)

species = data_dictionary["list_of_metabolite_symbols"]
tstart = 0.0
tstop = 3.0
tstep_sim = 0.01
tstep_exp = (tstop-tstart)/(length(Upper["GENE_CAT"])-1)
t_array_sim = tstart:tstep_sim:tstop
t_array_exp = tstart:tstep_exp:tstop
species_to_consider = ["M_glc_D_c"; "M_g6p_c"; "M_f6p_c"; "M_fdp_c"; "M_dhap_c"; "M_g3p_c"; "M_13dpg_c"; "M_3pg_c"; "M_2pg_c"; "M_pep_c"; "M_pyr_c"; "M_lac_D_c"; "M_accoa_c"; "M_ac_c"; "M_cit_c"; "M_icit_c"; "M_glx_c"; "M_akg_c"; "M_succoa_c"; "M_succ_c"; "M_fum_c"; "M_mal_L_c"; "M_oaa_c"; "M_6pgc_c"; "M_6pgl_c"; "M_ru5p_D_c"; "M_r5p_c"; "M_xu5p_D_c"; "M_s7p_c"; "M_e4p_c"; "M_2ddg6p_c"]
error = zeros(length(species),1)
for i in 1:length(species)
    key = species[i]
    if in(key,species_to_consider)
        u = Upper[key]
        l = Lower[key]
        m = Mean[key]
        sim = time_state_array[i,:]
        itp = interpolate((t_array_sim,),sim,Gridded(Linear()))
        sim_itp = itp[t_array_exp]
        error[i] = sum((max(0,sim_itp-u)+max(0,l-sim_itp))./(m+0.001))
        if isnan(error[i])
          error[i] = 0
        end
    end
end
return error
end
