using Distributed: @everywhere, pmap, addprocs

addprocs(4)


## 1) Go to directory and load external variables + functions ##

# moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
cd(@__DIR__)


# import structs
@everywhere include("../structsAndFunctions/structs.jl")
@everywhere using .structs: Molecule, Lasers, GeneralSettings

@everywhere begin
    # define molecule constants
    include("../simulationSettings/moleculeVariables.jl")
    using .moleculeVariables: SrF, CaF, BaF, MgF, CaOH, SrOH
    mol = SrF
end

# User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).
include("../structsAndFunctions/generateLaserSettings.jl")
using .laserSettings: generateLaserSettings
lasers = generateLaserSettings(mol)

# Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
include("../structsAndFunctions/generateGeneralSettings.jl")
using .generalSettings: generateGeneralSettings
general = generateGeneralSettings(mol)

# actual OBE simulation functinos
include("../structsAndFunctions/simulateIt.jl")
using .simulateIt: simulateOBE

# save simulation results and settings
include("../structsAndFunctions/saveSimulation.jl")
using .saveSimulation: saveToCsv
