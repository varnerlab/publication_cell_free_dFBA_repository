using PyPlot
using LaTeXStrings

#data = readdlm("data",'\t')
#data_13dpg = readdlm("combined_data_combo_13dpg")
#data = vcat(data,data_13dpg)
#data = readdlm("data_corrected",'\t')
#combined_data_combo_new = readdlm("combined_data_combo_new")
#combined_data_combo_13dpg = readdlm("13dpg")
#data = vcat(data,combined_data_combo_new,combined_data_combo_13dpg)

data = readdlm("FINAL",'\t')
data = readdlm("combined_data_SA_238ens",'\t')

labels = data[:,1]
error = data[:,2]
FVA_norm = data[:,3]

species_indices = Dict()
for i in 1:length(labels)
	species_indices[labels[i]] = i
end

error /= error[1]
FVA_norm /= FVA_norm[1]

m = 8
a = .3

no_glc = Array{Int64}(readdlm("no_glc"))
no_g6p = Array{Int64}(readdlm("no_g6p"))

#figure(figsize=(22.5,11.5))
figure(figsize=(12,8))

#semilogx(error[3:end],FVA_norm[3:end],"go",markersize=m,markeredgecolor="none",alpha=a)
#for i in 3:length(error)
#	if in(i,no_glc) && !in(i,no_g6p)
#		semilogx(error[i],FVA_norm[i],color="blue","o",markersize=m,markeredgecolor="none",alpha=a)
#	elseif in(i,no_g6p) && !in(i,no_glc)
#		semilogx(error[i],FVA_norm[i],"ro",markersize=m,markeredgecolor="none",alpha=a)
#	elseif in(i,no_glc) && in(i,no_g6p)
#		semilogx(error[i],FVA_norm[i],color="k","o",markersize=m,markeredgecolor="none",alpha=a)
#	else
#		semilogx(error[i],FVA_norm[i],"go",markersize=m,markeredgecolor="none",alpha=a)
#	end
#end
for i in 3:length(error)
	if in(i,no_glc)
		semilogx(error[i],FVA_norm[i],color="orange","o",markersize=m,markeredgecolor="none",alpha=a)
	else
		semilogx(error[i],FVA_norm[i],"go",markersize=m,markeredgecolor="none",alpha=a)
	end
end

#semilogx(error[2],FVA_norm[2],"ro",markersize=m,markeredgecolor="none")
semilogx(error[1],FVA_norm[1],"k*",markersize=m*2,markeredgecolor="none")
xlabel("error / base case error")
ylabel("flux uncertainty / base case flux uncertainty")
#xlabel("error")
#ylabel("flux uncertainty")

ax = collect(axis())
#ax[3] = 15000
#ax[4] = 25000
axis(ax)

#for i in 1:length(labels)
#	key = labels[i]
#	if in(key,["–GLC"; "–PYR"; "–ADP"; "–ATP"; "+CIT"; "–TYR"; "–UMP"; "+R5P"; "–MAL"])
#		annotate(key,xy=(error[species_indices[key]]+.005,FVA_norm[species_indices[key]]))
#	elseif in(key,["–AC"])
#		annotate(key,xy=(error[species_indices[key]],FVA_norm[species_indices[key]]-.003))
#	elseif in(key,["+AKG"])
#		annotate(key,xy=(error[species_indices[key]]-.082,FVA_norm[species_indices[key]]-.0015))
#	elseif in(key,["–PHE"])
#		annotate(key,xy=(error[species_indices[key]]-.035,FVA_norm[species_indices[key]]-.005))
#	elseif in(key,["+2PG"])
#		annotate(key,xy=(error[species_indices[key]]-.04,FVA_norm[species_indices[key]]-.005))
#	elseif in(key,["+F6P"])
#		annotate(key,xy=(error[species_indices[key]]+.007,FVA_norm[species_indices[key]]-.002))
#	elseif in(key,["+G6P"])
#		annotate(key,xy=(error[species_indices[key]]-.04,FVA_norm[species_indices[key]]-.005))
#	end
#end

savefig("New.pdf")

