using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
# Pkg.add("LaTeXStrings")
# Pkg.add("Random")
#Pkg.add("FFTW")
# Pkg.add("Plots")

using FFTW
using Plots
using LaTeXStrings
using Random

# delete old plots folder and create a new one
figdir = joinpath(@__DIR__, "Plots") # plotting path
rm(figdir,recursive = true)
mkpath(figdir)

# nicer plots 
 function plot_kwargs()
        fontsizes = (
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
            legendfontsize = 14,
        )
        (;
            #gridlinewidth = 2,            
            linewidth = 3,
            markersize = 6,
            markerstrokewidth = 0,
            fontsizes...,
            size = (600, 500),
        )
    end

include("pfc.jl")
include("gsh.jl")