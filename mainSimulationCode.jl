"""
unit of energy: hbar * Gamma
unit of velocity: Gamma / k, where k is the wavevector
unit of time: 1 / Gamma
unit of length: 1 / k (= wavelength / 2pi)
unit of force: 1e-3 * hbar * Gamma * k (??)
"""

using Distributed: @everywhere, pmap, addprocs
using BenchmarkTools: @time


## 1) Go to directory and load external variables + functions ##

# moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
cd(@__DIR__)

# import structs
include("./structsAndFunctions/structs.jl")
using .structs: Molecule, Lasers, GeneralSettings

# define molecule constants
include("./simulationSettings/moleculeVariables.jl")
using .moleculeVariables: SrF, CaF, BaF, MgF, CaOH, SrOH
mol = SrF

# User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).
include("./structsAndFunctions/generateLaserSettings.jl")
using .laserSettings: generateLaserSettings
lasers = generateLaserSettings(mol)

# Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
include("./structsAndFunctions/generateGeneralSettings.jl")
using .generalSettings: generateGeneralSettings
general = generateGeneralSettings(mol)

# actual OBE simulation functinos
include("./structsAndFunctions/simulateIt.jl")
using .simulateIt: simulateOBE

# save simulation results and settings
include("./structsAndFunctions/saveSimulation.jl")
using .saveSimulation: saveToCsv


## 2) Start OBE simulation, with user-defined parameters ##

displacements_list = repeat(general.displacementsInMM, inner=length(general.userSpeeds)*length(general.longSpeeds))
longSpeeds_list = repeat(general.longSpeeds, inner=length(general.userSpeeds), outer=length(general.displacementsInMM))
userSpeeds_list = repeat(general.userSpeeds, outer=length(general.longSpeeds)*length(general.displacementsInMM))
mol_list = fill(mol, length(displacements_list))
lasers_list = fill(lasers, length(displacements_list))
general_list = fill(general, length(displacements_list))

obeResults = pmap(simulateOBE, mol_list, lasers_list, general_list, displacements_list, userSpeeds_list, longSpeeds_list)


## 3) save simulation results and settings ##

saveToCsv(mol, lasers, general, displacements_list, longSpeeds_list, userSpeeds_list, obeResults)
        