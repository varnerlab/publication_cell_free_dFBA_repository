include("Include.jl")

data_dictionary, TXTL_dictionary, species_list, Mean, Upper, Lower = LoadDictionaries()

t0_exp = 0
tf_exp = 3
tstep_exp = 0.01
experimental_time_exp = [t0_exp:tstep_exp:tf_exp;]

metabolite_list = data_dictionary["list_of_metabolite_symbols"]
rxn_list = data_dictionary["list_of_reaction_strings"]

(n_species,n_rxn) = size(data_dictionary["stoichiometric_matrix"])
number_of_fluxes = find(rxn_list.=="tRNA_charging_M_val_L_c_CAT::16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA+16.0*M_h2o_c --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+16.0*M_ppi_c")[1]

t0 = 0
tf = 3
tstep = 0.01
experimental_time = t0:tstep:tf
length_time = length(experimental_time)
data_dictionary["tstep"] = tstep

# this determines which parameters are fitted to data
dataset_index = [15; 18; 22; 24; 27; 29; 31; 34; 38; 41; 42; 59; 60; 61; 66; 69; 70; 76; 79; 83; 84; 86; 88; 89; 106; 110; 114; 121; 123; 126; 128; 130; 132; 133; 134; 135; 140]
non_dataset_index = collect(1:n_species)
deleteat!(non_dataset_index,dataset_index)

species_list = data_dictionary["list_of_metabolite_symbols"]
species_list_dataset = species_list[dataset_index]
species_list_non_dataset = species_list[non_dataset_index]

species_constraint_index = copy(dataset_index)

plot_color = "black"

num_constraint_sets = 1#size(species_constraint_array,1)

# not currently used
parameter = ones(202,2)
parameter[:,1] = parameter[:,1]*110
parameter[:,2] = parameter[:,2]*0.03

initial_conditions = data_dictionary["initial_conditions"]

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

constraint_index_array = collect(1:num_constraint_sets)

Exit_flag = Int64[]
uptake_array = 0
flux_array = 0

if !isdir("Ens")
	mkdir("Ens")
end
if !isdir("Ens/Best")
	mkdir("Ens/Best")
end

syn_data = Dict{AbstractString,Any}()

fva_time_array = collect(0.0:0.25:3.0)

fluxes_to_optimize = [[2 3], [23 24], 26, 37, [73 74], 21, 75]

species_to_vary = readdlm("species_to_vary")
fixed_indices = Int64[]
for i in 1:length(species_list)
	if !in(species_list[i],species_to_vary)
		push!(fixed_indices,i)
	end
end

species_boolean = zeros(size(species_list))
for i in 1:length(species_list)
	species_boolean[i] = in(i,dataset_index)
end

switch_prob = .25

if ~isfile("Ens/num_dir")
    writedlm("Ens/num_dir",0)
end

alpha = 1
num_to_reset = 50

if isfile("Ens/Best/error")
	error_best = readdlm("Ens/Best/error")[1]
	species_boolean = readdlm("Ens/Best/species_boolean")
else
	error_best = readdlm("error_best")[1]
end
species_boolean_best = copy(species_boolean)
error = copy(error_best)

error_bc = 3608115.9638814344
flux_uncertainty_bc = 4000 # based on branch points only 20932.739199515625

num_dir = convert(Int64,readdlm("Ens/num_dir")[1]) # Number of pre-existing directories

num_iter = 10000
for i = num_dir+1:num_dir+num_iter
	# Reset periodically
	if rem(i,num_to_reset) == 0
		species_boolean_new = copy(species_boolean_best)
	else
		species_boolean_new = copy(species_boolean)
	end
	
	# Perturb constraint array
	for j in 1:length(species_boolean_new)
#		if !in(j,fixed_indices)
			if rand() < switch_prob
				species_boolean_new[j] = 1-species_boolean_new[j]
			end
#		end
	end
	
	species_boolean_new[fixed_indices] = 0
	
	num_diff = Int64(sum(abs(species_boolean_new-species_boolean)))
	
	# discretized dfba
	state_array = deepcopy(initial_conditions)
	index = 1
	
	species_constraint_index_new = find(species_boolean_new)
	
	syn_data_upper = zeros(n_species,length_time)
	syn_data_lower = zeros(n_species,length_time)
	syn_mean = zeros(n_species,length_time)
	for syn_index in species_constraint_index_new
		itp = interpolate((experimental_time_exp,),Mean[metabolite_list[syn_index]],Gridded(Linear()))
		Mean_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Lower[metabolite_list[syn_index]],Gridded(Linear()))
		Lower_itp = itp[experimental_time]
		itp = interpolate((experimental_time_exp,),Upper[metabolite_list[syn_index]],Gridded(Linear()))
		Upper_itp = itp[experimental_time]
		syn_data_upper[syn_index,:] = Upper_itp
		syn_data_lower[syn_index,:] = Lower_itp
		syn_mean[syn_index,:] = Mean_itp
	end
	syn_data["mean"] = syn_mean
	syn_data["upper"] = syn_data_upper
	syn_data["lower"] = syn_data_lower
	
	Exit_flag = Int64[]
	uptake_array = 0
	flux_array = 0
	
	FVA_min = Dict()
	FVA_max = Dict()
	Percentage_failed_tps_min = Dict()
	Percentage_failed_tps_max = Dict()
	for rxn_idx in fluxes_to_optimize
		key = split(rxn_list[rxn_idx[1]],':')[1]
		FVA_min[key] = zeros(length(fva_time_array),191) # Only record 191 fluxes
		FVA_max[key] = zeros(length(fva_time_array),191) # Only record 191 fluxes
	end
	
	#dfba
	fva_time = experimental_time[1:end-1]
	for j in 1:length(fva_time)
		time_index = fva_time[j]
		# constraint
		data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
		data_dictionary["objective_coefficient_array"][171] = -1
		data_dictionary = Bounds(data_dictionary,TXTL_dictionary,state_array)
		species_constraints = calculate_constraints(state_array, parameter, species_constraint_index_new, syn_data, index, tstep)
		data_dictionary["species_bounds_array"] = species_constraints
		#  data_dictionary["default_flux_bounds_array"] = flux_constraints
		data_dictionary["state_array"] = state_array

		# solve the lp problem -
		(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
		push!(Exit_flag,exit_flag)

		# mass balance for concentration
		state_array = state_array + (uptake_array * tstep) # cell free mass balanace
		state_array[101] = initial_conditions[101]
		Z0 = objective_value
		if FVA_run && in(time_index,fva_time_array)
			for rxn_idx in fluxes_to_optimize
				key = split(rxn_list[rxn_idx[1]],':')[1]
				if length(rxn_idx) == 1
					#constraint objective flux to value determined
					data_dictionary["default_flux_bounds_array"][171] = Z0

					# minimize flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][rxn_idx] = 1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
					FVA_min[key][Int64(time_index*4+1),:] = flux_array[1:191]

					# maximize flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][rxn_idx] = -1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
					FVA_max[key][Int64(time_index*4+1),:] = flux_array[1:191]
				else
					#constraint objective flux to value determined
					data_dictionary["default_flux_bounds_array"][171] = Z0

					# maximize reverse flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][rxn_idx[2]] = -1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
					FVA_min[key][Int64(time_index*4+1),:] = flux_array[1:191]

					# maximize forward flux
					data_dictionary["objective_coefficient_array"] = zeros(number_of_fluxes)
					data_dictionary["objective_coefficient_array"][rxn_idx[1]] = -1
					(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
					FVA_max[key][Int64(time_index*4+1),:] = flux_array[1:191]
				end
			end
		end # if FVA_run
		
		index = index + 1
		time_state_array[:,index] = state_array
		time_flux_array[:,index] = flux_array

	end # for time_index in [t0:tstep:tf-tstep;]
	
	if sum(5-Exit_flag) == 0
		error_vector = CalcError(Upper,Lower,Mean,experimental_time,time_state_array,data_dictionary)
		
		FVA_MIN = zeros(length(fva_time_array),length(fluxes_to_optimize))
		for j in 1:length(fluxes_to_optimize)
			rxn_idx = fluxes_to_optimize[j]
			key = split(rxn_list[rxn_idx[1]],':')[1]
			if length(rxn_idx) == 1
				FVA_MIN[:,j] = FVA_min[key][:,rxn_idx]
			else
				FVA_MIN[:,j] = FVA_min[key][:,rxn_idx[1]]-FVA_min[key][:,rxn_idx[2]]
			end
		end
		FVA_MAX = zeros(length(fva_time_array),length(fluxes_to_optimize))
		for j in 1:length(fluxes_to_optimize)
			rxn_idx = fluxes_to_optimize[j]
			key = split(rxn_list[rxn_idx[1]],':')[1]
			if length(rxn_idx) == 1
				FVA_MAX[:,j] = FVA_max[key][:,rxn_idx]
			else
				FVA_MAX[:,j] = FVA_max[key][:,rxn_idx[1]]-FVA_max[key][:,rxn_idx[2]]
			end
		end
		flux_uncertainty = norm(FVA_MAX-FVA_MIN)
	
		error_new = sum(error_vector)*flux_uncertainty
#		error_new = length(species_constraint_index_new)
		
		if sum(error_vector) > error_bc
			error_new *= 1e9
		end
		if flux_uncertainty > flux_uncertainty_bc
			error_new *= 1e9
		end
		if sum(species_boolean_new) == 0
			error_new *= 1e9
		end
		
		acc_prob = exp(alpha*(error-error_new)/error)
		
		error_round = round(error/10^floor(log(10,error)),3)*10^floor(log(10,error))
		error_new_round = round(error_new/10^floor(log(10,error_new)),3)*10^floor(log(10,error_new))
		error_best_round = round(error_best/10^floor(log(10,error_best)),3)*10^floor(log(10,error_best))
		acc_prob_round = round(acc_prob,2)
		
		# If a new best is achieved, overwrite parameters and cost and save the relevant information to Best directory
		if error_new < error_best
		    # Save to Best directory
		    writedlm("Ens/Best/species_boolean",species_boolean_new)
		    writedlm("Ens/Best/time_state_array",time_state_array)
		    writedlm("Ens/Best/exit_flag",Exit_flag)
		    writedlm("Ens/Best/state_error",sum(error_vector))
		    writedlm("Ens/Best/flux_uncertainty",flux_uncertainty)
		    writedlm("Ens/Best/error",error_new)
		    species_boolean_best = copy(species_boolean_new)
		    error_best = copy(error_new)
		    println("$i: error_new = $error_new_round, error = $error_round, best = $error_best_round, acc_prob = $acc_prob_round, num_diff = $num_diff, NEW BEST")
		else
		    println("$i: error_new = $error_new_round, error = $error_round, best = $error_best_round, acc_prob = $acc_prob_round, num_diff = $num_diff")
		end
		# If new cost is lower than previous cost, choose new params as reference point for perturbation
		if rand() < acc_prob # Accept all better sets and some worse sets
		    # Create directory for next sample and save the relevant information
		    num_dir += 1
			writedlm("Ens/num_dir",i) # Record new number of directories
		    mkdir("Ens/$num_dir")
		    writedlm("Ens/$num_dir/species_boolean",species_boolean_new)
		    writedlm("Ens/$num_dir/time_state_array",time_state_array)
		    writedlm("Ens/$num_dir/exit_flag",Exit_flag)
		    writedlm("Ens/$num_dir/state_error",sum(error_vector))
		    writedlm("Ens/$num_dir/flux_uncertainty",flux_uncertainty)
		    writedlm("Ens/$num_dir/error",error_new)
		    species_boolean = copy(species_boolean_new)
		    error = copy(error_new)
			for rxn_idx in fluxes_to_optimize
				key = split(rxn_list[rxn_idx[1]],':')[1]
				mkdir("Ens/$num_dir/$key")
				writedlm("Ens/$num_dir/$key/FVA_min",FVA_min[key])
				writedlm("Ens/$num_dir/$key/FVA_max",FVA_max[key])
			end
		end
	end
end
