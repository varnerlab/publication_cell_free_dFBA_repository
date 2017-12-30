function TXTLDictionary(time_start,time_stop,time_step)

#Blocked Reactions
degradation_rxn = [100:1:116;];

AA_synthesis_rxn = [81;87;88;89;90;91;92;93;94;96;97;98;99;];

# setup promoter models -
Promoter = Float64[]
push!(Promoter,10) #K1

# Setup the mRNA elongation rate, and global translation
RNAP_concentration_nM = 75;                          # 60-75nM (ACS SynBio Garamella 2016)
RNAP_elongation_rate = 25;                          # >5 NT/s (ACS SynBio Garamella 2016)
RIBOSOME_concentration = 0.0016;                      # 0.0016mM with 72% MaxActive (Underwood, Swartz, Puglisi 2005 Biotech Bioeng) & <0.0023mM (ACS SynBio Garamella 2016)
RIBOSOME_elongation_rate = 2;                      # >1 (ACS SynBio Garamella 2016) & 1.5 AA/sec (Underwood, Swartz, Puglisi 2005 Biotech Bioeng)
#number_of_cells = 2e8;                                # 1e9 cells/ml

# Protein of Interest Sequences
mRNA_length = 660;
protein_length = 219;
gene_copies = 3.125e10;
mRNA_degradation_rate = 10#5.2;
polysome_gain = 5#10;
volume = 15e-6
inducer = 35;
plasmid_saturation_coefficient = 3.5;
mRNA_saturation_coefficient = 0.045;

# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
TXTL_dictionary = Dict{AbstractString,Any}()
TXTL_dictionary["AA_degradation_rxn"] = degradation_rxn
TXTL_dictionary["AA_synthesis_rxn"] = AA_synthesis_rxn
TXTL_dictionary["Promoter"] = Promoter
TXTL_dictionary["RNAP_concentration_nM"] = RNAP_concentration_nM
TXTL_dictionary["RNAP_elongation_rate"] = RNAP_elongation_rate
TXTL_dictionary["RIBOSOME_concentration"] = RIBOSOME_concentration
TXTL_dictionary["RIBOSOME_elongation_rate"] = RIBOSOME_elongation_rate
TXTL_dictionary["mRNA_length"] = mRNA_length
TXTL_dictionary["protein_length"] = protein_length
TXTL_dictionary["gene_copies"] = gene_copies
TXTL_dictionary["mRNA_degradation_rate"] = mRNA_degradation_rate
TXTL_dictionary["polysome_gain"] = polysome_gain
TXTL_dictionary["volume"] = volume
TXTL_dictionary["inducer"] = inducer
TXTL_dictionary["plasmid_saturation_coefficient"] = plasmid_saturation_coefficient
TXTL_dictionary["mRNA_saturation_coefficient"] = mRNA_saturation_coefficient
# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
return TXTL_dictionary
end
