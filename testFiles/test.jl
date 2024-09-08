# moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
cd(@__DIR__)


# import structs
include("../structsAndFunctions/structs.jl")
using .structs: Molecule, Lasers, GeneralSettings


# define molecule constants
include("../simulationSettings/moleculeVariables.jl")
using .moleculeVariables: SrF, CaF, BaF, MgF, CaOH, SrOH
mol = SrF


