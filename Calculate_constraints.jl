function calculate_constraints(state_array, parameters, constraint_index, syn_data,time_index,tstep)
  alpha = 1.1 # looseness/tightness of constraints
  beta = 10
  n_species = length(state_array[:,1])
  syn_upper = syn_data["upper"]
  syn_lower = syn_data["lower"]
  syn_mean = syn_data["mean"]
  state_array = max(0,state_array)

  species_drift_amplitude = data_dictionary["species_drift_amplitude"]
  species_drift_decay_rate = data_dictionary["species_drift_decay_rate"]
  species_drift = species_drift_amplitude*exp(-time_index*tstep*species_drift_decay_rate)

  lower_bound = 0
  upper_bound = 100000000
  species_lower_bound = max(min(state_array*(-species_drift)/tstep,-0.00000000001/tstep),(lower_bound-state_array)/tstep);
  species_upper_bound = min(max(state_array*(species_drift)/tstep,0.00000000001/tstep),(upper_bound-state_array)/tstep);

  # constrain the species to synthetic data
  for c_index in constraint_index
      species_lower_bound[c_index] = (syn_lower[c_index,time_index+1]-state_array[c_index]) / tstep
      species_upper_bound[c_index] = (syn_upper[c_index,time_index+1]-state_array[c_index]) / tstep
  end
  # oxygen
  species_lower_bound[101] = -10000
  species_upper_bound[101] = 10000
  #h2o
  species_lower_bound[72] = -10000
  species_upper_bound[72] = 10000




  species_constraints = [species_lower_bound species_upper_bound]
  #flux_constraints = [flux_lower_bound flux_upper_bound]

  return species_constraints#, flux_constraints
end
