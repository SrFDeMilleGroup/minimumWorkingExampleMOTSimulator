# 1) Go to directory and load external variables + functions
cd(@__DIR__) # moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well

include("../structsAndFunctions/structs.jl")
using .structs: Molecule, Lasers, GeneralSettings

include("../simulationSettings/moleculeVariables.jl")
using .moleculeVariables: SrF, CaF, BaF, MgF, CaOH, SrOH
mol = SrF

# 3) Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
include("../structsAndFunctions/generateGeneralSettings.jl")
using .generalSettings: generateGeneralSettings
general = generateGeneralSettings(mol)

# 5) User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).
include("../structsAndFunctions/generateLaserSettings.jl")
using .laserSettings: generateLaserSettings
lasers = generateLaserSettings(mol)






include("../structsAndFunctions/simulateIt.jl")
using .simulateIt: simulateOBE