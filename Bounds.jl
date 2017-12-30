function Bounds(DF,TXTL,state_array)
FB = DF["default_flux_bounds_array"]

Vmax = 600 # 1000 s^-1 * (3600 s/h) * (0.000167 mM)

FB[1:191,1] = 0
FB[1:191,2] = Vmax

FB[192:end,1] = 0
FB[192:end,2] = 0
FB[197,2] = 100 # hydrogen transport

#==============================================TXTL=====================================================#
#RNAP_elongation_rate = TXTL["RNAP_elongation_rate"];
RIBOSOME_concentration = state_array[141]#TXTL["RIBOSOME_concentration"];
RIBOSOME_elongation_rate = TXTL["RIBOSOME_elongation_rate"];
#mRNA_degradation_rate = TXTL["mRNA_degradation_rate"]*2;
mRNA_length = TXTL["mRNA_length"];
protein_length = TXTL["protein_length"];
gene_copies = TXTL["gene_copies"];
volume = TXTL["volume"];
polysome_amplification = TXTL["polysome_gain"];
plasmid_saturation_coefficient = TXTL["plasmid_saturation_coefficient"]*1e-6;
mRNA_saturation_coefficient = TXTL["mRNA_saturation_coefficient"]
Promoter = TXTL["Promoter"]
inducer = TXTL["inducer"]
#====================================Transcription===================================================#
#Compute the promoter strength P -
P = 0.9 #0.8*(K1)/(1+K1);
RNAP_elongation_rate = TXTL["RNAP_elongation_rate"]
plasmid_saturation_coefficient = TXTL["plasmid_saturation_coefficient"]
# 1 GENE_CAT
TX_Vmax = RNAP_elongation_rate/mRNA_length*3600*P
TX = TX_Vmax*state_array[143]*state_array[1]/(plasmid_saturation_coefficient+state_array[1]);

#====================================Translation===================================================#
polysome_amplification = TXTL["polysome_amplification"]
RIBOSOME_elongation_rate = TXTL["RIBOSOME_elongation_rate"]
mRNA_saturation_coefficient = TXTL["mRNA_saturation_coefficient"]
mRNA_degradation_rate = TXTL["mRNA_degradation_rate"]
# 141 RIBOSOME
# 144 mRNA_CAT
TL_Vmax = polysome_amplification*(3*RIBOSOME_elongation_rate)*(1/mRNA_length)*3600
TL = TL_Vmax*state_array[141]*state_array[144]/(mRNA_saturation_coefficient+state_array[144])
DEG = state_array[144]*mRNA_degradation_rate
#===================================================================================================#
FB[167,1] = TX #transcriptional initiation
FB[168,1] = TX #transcriptional elongation
FB[167,2] = TX #transcriptional initiation
FB[168,2] = TX #transcriptional elongation
FB[169,1] = 0 #mRNA_degradation
FB[169,2] = DEG #mRNA_degradation
FB[170,1] = 0 #translation initiation
FB[171,1] = 0 #translation elongation
FB[170,2] = TL #translation initiation
FB[171,2] = TL #translation elongation
FB[172:191,1] = 0 #tRNA charging
FB[172:191,2] = 50*TL #tRNA charging
FB = max(0,FB)
DF["default_flux_bounds_array"] = FB
return DF
end
