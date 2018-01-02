include("DataDictionary.jl")
include("LoadDictionaries.jl")
include("LoadSynthData.jl")
include("FluxDriver.jl")
include("calculate_constraints.jl")
include("CalculateConstraints.jl")
include("Bounds.jl")
include("TXTLDictionary.jl")
include("Utility.jl")
include("CalcError.jl")
include("Plot.jl")
include("PlotSeparate.jl")
include("PlotErrorbar.jl")
include("dFBA.jl")
include("RunFBA.jl")
using PyPlot
using Interpolations
using LaTeXStrings
using GLPK
term_out(GLPK.OFF)
