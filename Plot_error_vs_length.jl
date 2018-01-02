using PyPlot
using LaTeXStrings

data = readdlm("error_vs_length",'\t')

error = data[:,1]
l = data[:,2]

species_indices = Dict()
for i in 1:length(labels)
	species_indices[labels[i]] = i
end

error /= error[1]
#l /= l[1]

m = 8
a = .3

no_glc = Array{Int64}(readdlm("no_glc"))
no_g6p = Array{Int64}(readdlm("no_g6p"))

figure(figsize=(12,8))
for i in 1:length(error)
	if in(i,no_glc)
		semilogx(error[i],l[i],color="orange","o",markersize=m,markeredgecolor="none",alpha=a)
	else
		semilogx(error[i],l[i],"go",markersize=m,markeredgecolor="none",alpha=a)
	end
end

semilogx(error[1],l[1],"k*",markersize=m*2,markeredgecolor="none")
xlabel("error / base case error")
ylabel("length")
savefig("New.pdf")

