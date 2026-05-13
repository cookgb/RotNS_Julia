import Pkg
Pkg.activate(".")

# Prevent CondaPkg from resolving
ENV["JULIA_CONDAPKG_OFFLINE"] = "true"

using Plots HDF5 SciMLBase OrdinaryDiffEq NonlinearSolve
using DelimitedFiles ArgParse

println("✓ Ready!")