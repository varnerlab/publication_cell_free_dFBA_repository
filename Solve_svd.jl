tic()
include("Include.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)
TXTL_dictionary = TXTLDictionary(0,0,0)
number_of_fluxes = length(data_dictionary["default_flux_bounds_array"][:,1])
number_of_species = length(data_dictionary["species_bounds_array"][:,1])
#Set objective reaction
data_dictionary["objective_coefficient_array"][171] = -1;
#-------------------------------------------------------------------------------------------------#Define case number
# 1 = Amino Acid Uptake & Synthesis
# 2 = Amino Acid Uptake w/o Synthesis
# 3 = Amino Acid Synthesis w/o Uptake
case = 1
if case == 1
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if case == 2
  data_dictionary["AASyn"] = 0;
  data_dictionary["AAUptake"] = 30
  data_dictionary["AASecretion"] = 0;
end
if case == 3
  data_dictionary["AASyn"] = 100;
  data_dictionary["AAUptake"] = 0
  data_dictionary["AASecretion"] = 100;
end
#-------------------------------------------------------------------------------------------------
volume = TXTL_dictionary["volume"]
plasmid_concentration = 5;
gene_copy_number = (volume/1e9)*(6.02e23)*plasmid_concentration;
TXTL_dictionary["gene_copies"] = gene_copy_number
data_dictionary = DataDictionary(0,0,0)
#Set objective reaction
data_dictionary["objective_coefficient_array"][194] = -1;
data_dictionary["AASyn"] = 100;
data_dictionary["AAUptake"] = 30
data_dictionary["AASecretion"] = 0;
data_dictionary["GlcUptake"] = 30#*rand(1)[1];
data_dictionary["Oxygen"] = 30#+2*rand(1)[1];
TXTL_dictionary["RNAP_concentration_nM"] = 75 #+ (80-60)*rand(1)[1];
TXTL_dictionary["RNAP_elongation_rate"] = 25 #+ (30-20)*rand(1)[1];
TXTL_dictionary["RIBOSOME_concentration"] = 0.0016# + (0.0018-0.0012)*rand(1)[1];
TXTL_dictionary["RIBOSOME_elongation_rate"] = 2# + (3-1.5)*rand(1)[1];
data_dictionary["Oxygen"] = 100;
data_dictionary["AAUptake"] = 1 #+ (1.5-0.5)*rand(1)[1];
data_dictionary["GlcUptake"] = 30# + (38-28)*rand(1)[1];
data_dictionary["AC"] = 8# + (12-6)*rand(1)[1];
data_dictionary["SUCC"] = 2# + (3-1)*rand(1)[1];
data_dictionary["PYR"] = 7# + (10-5)*rand(1)[1];
data_dictionary["MAL"] = 3 #+ (5-2)*rand(1)[1];
data_dictionary["LAC"] = 8 #+ (12-8)*rand(1)[1];
data_dictionary["ENERGY"] = 1#1*rand(1)[1];
data_dictionary["ALA"] = 2 #+ (4-3)*rand(1)[1];
data_dictionary["ASP"] = 1#+rand(1)[1];
data_dictionary["GLN"] = 0.25#*rand(1)[1];
data_dictionary["LysUptake"] = 1.5#*rand(1)[1];
data_dictionary["LysSecretion"] = 0.5#1*rand(1)[1];
#TXTL_dictionary["RNAP_concentration_nM"] = 60 + (80-60)*rand(1)[1];
#TXTL_dictionary["RNAP_elongation_rate"] = 20 + (30-20)*rand(1)[1];
#TXTL_dictionary["RIBOSOME_concentration"] = 0.0012 + (0.0018-0.0012)*rand(1)[1];
#TXTL_dictionary["RIBOSOME_elongation_rate"] = 1.5 + (3-1.5)*rand(1)[1];
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]
#Mean = Unpack("./data/Mean")
#Std = Unpack("./data/Std")
#Upper = Unpack("./data/Upper")
#Lower = Unpack("./data/Lower")

species_list = readdlm("./data/species_list")
Mean_raw = readdlm("./data/Mean")
Std_raw = readdlm("./data/Std")
Upper_raw = readdlm("./data/Upper")
Lower_raw = readdlm("./data/Lower")
Mean = Dict()
Std = Dict()
Upper = Dict()
Lower = Dict()
for i in 1:length(species_list)
	rxn_idx = species_list[i]
	Mean[rxn_idx] = Mean_raw[:,i]
	Std[rxn_idx] = Std_raw[:,i]
	Upper[rxn_idx] = Upper_raw[:,i]
	Lower[rxn_idx] = Lower_raw[:,i]
end


#amplitude = 2, decay_rate = 1
data_dictionary["species_drift_amplitude"] = 2 # percent allowed to drift for unconstraint species for each iteration (0.01 hr)
data_dictionary["species_drift_decay_rate"] = 1
metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
t0 = 0
tf = 3
tstep = 0.01
experimental_time = t0:tstep:tf
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

#notes: two possible implementations ---
#  1. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data be open or a function of concentration
#  2. select subset of species to be constraint to synthetic data
#    have everything else that are not constrained to data constrained at steady state (sv=0)
tRNA_index = Int64[]
for i in 1:length(metabolite_list)
	rxn_idx = metabolite_list[i]
	if length(rxn_idx) > 5
		if rxn_idx[end-4:end] == "_tRNA"
			push!(tRNA_index,i)
		end
	end
end

initial_conditions = zeros(n_species,1)
for species_index in 1:n_species
	if in(species_index,tRNA_index)
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][2]
	else
		initial_conditions[species_index] = Mean[metabolite_list[species_index]][1]
	end
end
initial_conditions[142] = .001
initial_conditions[144] = .0001

# this determines which parameters are fitted to data
dataset_index_table = readdlm("FVA_svd/svd_measurement_constraints")
dataset_index = sort(dataset_index_table[:,1])
non_dataset_index = collect(1:number_of_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

#species_constraint_index = vcat(dataset_index,tRNA_index)
species_constraint_index = copy(dataset_index)
#species_constraint_index = dataset_index[1:10]

plot_color = "orangered"

data_dictionary["initial_conditions"] = initial_conditions
# Original values:
GENE_IC = 5.19e-6
RNAP_IC = 6.75e-5
RIBO_IC = 0.002
mRNA_IC = 0.0001
tRNA_IC = 1.6
TX_Vmax = 122.7
plasmid_saturation_coefficient = 3.5e-6
TL_Vmax = 163.6
mRNA_saturation_coefficient = 0.045
mRNA_degradation_rate = 20
params = vcat(GENE_IC, RNAP_IC, RIBO_IC, mRNA_IC, tRNA_IC, TX_Vmax, plasmid_saturation_coefficient, TL_Vmax, mRNA_saturation_coefficient, mRNA_degradation_rate)
params[8] = 40
RNAP_IC = 1e-3
RIBO_IC = 0.0021
RNAP_elongation_rate = 25
plasmid_saturation_coefficient = 116e-6
polysome_amplification = 5
RIBOSOME_elongation_rate = 2
mRNA_saturation_coefficient = 0.045
mRNA_degradation_rate = 10
params = vcat(RNAP_IC, RIBO_IC, RNAP_elongation_rate, plasmid_saturation_coefficient, polysome_amplification, RIBOSOME_elongation_rate, mRNA_saturation_coefficient, mRNA_degradation_rate)
data_dictionary["initial_conditions"][143] = params[1]
data_dictionary["initial_conditions"][141] = params[2]
TXTL_dictionary["RNAP_elongation_rate"] = params[3]
TXTL_dictionary["plasmid_saturation_coefficient"] = params[4]
TXTL_dictionary["polysome_amplification"] = params[5]
TXTL_dictionary["RIBOSOME_elongation_rate"] = params[6]
TXTL_dictionary["mRNA_saturation_coefficient"] = params[7]
TXTL_dictionary["mRNA_degradation_rate"] = params[8]

species_constraint_array  = dataset_index
Labels = "svd"
num_constraint_sets = 1

# set up initial conditions (use kinetic model data)
time_state_array = zeros(n_species,length_time)
time_state_array[:,1] = initial_conditions
time_flux_array = zeros(n_rxn,length_time)
FBA = data_dictionary["default_flux_bounds_array"]
SBA = data_dictionary["species_bounds_array"]

EXIT_FLAG = zeros(length_time-1,num_constraint_sets)

overall_state_array = zeros()
# loop for all constraint subsets
FVA_run = true
fva_flux_array_1_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_max = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_1_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_2_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_3_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)
fva_flux_array_4_min = zeros(number_of_fluxes,number_of_fluxes,num_constraint_sets)

constraint_index_array = collect(1:num_constraint_sets)

if !isdir("FVA_svd")
  mkdir("FVA_svd")
end
for constraint_index in 1 #constraint_index_array[1:end] ;println(Labels[constraint_index])
	constraint_subset = "svd"#Labels[constraint_index] #"base_case"

	if !isdir("FVA_svd/$constraint_subset")
		mkdir("FVA_svd/$constraint_subset")
	end
  writedlm("FVA_svd/$constraint_subset/constraint_set",dataset_index_table)
	# discretized dfba
	state_array = deepcopy(initial_conditions)
	index = 1

	species_constraint_index = Array{Int64}(species_constraint_array[constraint_index,:])
#	species_constraint_index = copy(dataset_index)

	syn_data_upper = zeros(n_species,length_time)
	syn_data_lower = zeros(n_species,length_time)
	syn_mean = zeros(n_species,length_time)
	for syn_index in species_constraint_index
		itp = interpolate((experimental_time_exp,),Mean[metabolite_list[syn_index]],Gridded(Linear()))
		Mean_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Std[metabolite_list[syn_index]],Gridded(Linear()))
		Std_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Lower[metabolite_list[syn_index]],Gridded(Linear()))
		Lower_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Upper[metabolite_list[syn_index]],Gridded(Linear()))
		Upper_itp = itp[experimental_time]
		syn_data_upper[syn_index,:] = Upper_itp
		syn_data_lower[syn_index,:] = Lower_itp
		syn_mean[syn_index,:] = Mean_itp
	end
	syn_data = Dict{AbstractString,Any}()
	syn_data["mean"] = syn_mean
	syn_data["upper"] = syn_data_upper
	syn_data["lower"] = syn_data_lower
    syn_data["lower"][3,:] = syn_data_lower[3,:]/5.3

    Exit_flag = Int64[]
    uptake_array = 0
    flux_array = 0

	FVA_min = Dict()
	FVA_max = Dict()
	Percentage_failed_tps_min = Dict()
	Percentage_failed_tps_max = Dict()
	for rxn_idx in 1:194
		FVA_min[rxn_idx] = zeros(length(t0:tstep:tf-tstep),194) # Only record 194 fluxes
		FVA_max[rxn_idx] = zeros(length(t0:tstep:tf-tstep),194) # Only record 194 fluxes
	end

	#dfba
	fva_time = experimental_time[1:end-1]
	for t in 1:length(fva_time)
		time_index = fva_time[t] ;println(time_index)
		# constraint
		data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
		data_dictionary["objective_coefficient_array"][171] = -1
		data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array)
		species_constraints = calculate_constraints(state_array,parameter,species_constraint_index,syn_data,index,tstep)
		data_dictionary["species_bounds_array"] = species_constraints
		data_dictionary["state_array"] = state_array

		# solve the lp problem -
		(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
		push!(Exit_flag,exit_flag)

		# mass balance for concentration
		state_array = state_array + (uptake_array * tstep) # cell free mass balanace
		state_array[101] = initial_conditions[101]
		Z0 = objective_value
		if FVA_run && in(time_index,collect(0.0:0.1:3.0))
			for flux_index in 1:194
				#constraint objective flux to value determined
				data_dictionary["default_flux_bounds_array"][171] = 1*Z0
				data_dictionary["default_flux_bounds_array"][171] = Z0

				# minimize flux
				data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
				data_dictionary["objective_coefficient_array"][flux_index] = -1
				(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

				rxn_idx = copy(flux_index)

				FVA_min[rxn_idx][t,:] = flux_array[1:194]
				Percentage_failed_tps_min[rxn_idx] = sum(5-exit_flag)/length(experimental_time[1:end-1])

#		            @show flux_index, constraint_index
#		            @show flux_array
#					eval(parse("fva_flux_array_"fva_flag*"_max"))[:,flux_index,constraint_index] = flux_array

				# maximize flux
				data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
				data_dictionary["objective_coefficient_array"][flux_index] = 1
				(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

				FVA_max[rxn_idx][t,:] = flux_array[1:194]
				Percentage_failed_tps_max[rxn_idx] = sum(5-exit_flag)/length(experimental_time[1:end-1])
			end
		end # if FVA_run

		index = index + 1
		time_state_array[:,index] = state_array
		time_flux_array[:,index] = flux_array

	end # for time_index in [t0:tstep:tf-tstep;]

	writedlm("FVA_svd/$constraint_subset/time_state_array.txt",time_state_array)
	for rxn_idx in 1:194
		if !isdir("FVA_svd/$constraint_subset/$rxn_idx")
			mkdir("FVA_svd/$constraint_subset/$rxn_idx")
		end
		writedlm("FVA_svd/$constraint_subset/$rxn_idx/FVA_min",FVA_min[rxn_idx])
		writedlm("FVA_svd/$constraint_subset/$rxn_idx/FVA_max",FVA_max[rxn_idx])
		writedlm("FVA_svd/$constraint_subset/$rxn_idx/Percentage_failed_tps_min",Percentage_failed_tps_min[rxn_idx])
		writedlm("FVA_svd/$constraint_subset/$rxn_idx/Percentage_failed_tps_max",Percentage_failed_tps_max[rxn_idx])

	end

end # for constraint_index in constraint_index_array

plot_flag = false
if plot_flag
	Plot(plot_color)
	PlotSeparate(plot_color)
end

time_elapsed = toc()
