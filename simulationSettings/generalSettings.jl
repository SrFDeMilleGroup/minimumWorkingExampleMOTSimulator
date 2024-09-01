"""
simulationType: string, e.g., redMOT, blueMOT, transCooling, etc.

numTrialsPerValueSet: number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds)

displacementsInMM: mm, initial displacements from either origin (if 3D) or else z-axis (if 2D) in mm

initDispDir: if "XY", will force initial r to go along (x+y)/sqrt(2). 
Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = "Same". If "Z", force go along Z

longSpeeds: longitudinal speeds in units of Gamma/k (=4.4m/s for SrF), doesn't matter for 3D sims, sets vel to 140 m/s

userSpeeds: in unit of Gamma/k (=4.4m/s for SrF), speeds in xy plane (for 2d force profile) or in 3D

velDirRelToR: relative direction of molecule velocity, w.r.t. initial displacement r. Options are ["Same", "Orthogonal", "Opposite", "Random"]

forceProfile: either "ThreeD", (forces calculated are (f \\dot r) / |r|, (f \\dot v) / |v|), 
or "TwoD" (f \\dot (rx,ry,0) / |(rx,ry,0)|, f \\dot (vx,vy,0) / |(vx,vy,0|, and fz are all calculated)

bGradReal: in units Gauss/cm, unless bFieldSetting == "Static", then this becomes the static field in Gauss

bFieldSetting: can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT"),
or static "Static" (2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).

"""


# A) parameters for quick test of restoring force
simulationType = "SrFRedMOTNormalValues"
numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [0.5, 1.5, 3.0, 4.5, 6.0, 7.5]
displacementsInMM::Vector{Float64} = [0.5, 3.0]
initDispDir::String = "XY"
longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [-4, -3, -2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2, 3, 4]
userSpeeds::Vector{Float64} = [0.1, 0.5]
velDirRelToR::String = "Same"
forceProfile::String = "ThreeD"
bGradReal::Float64 = 12.5
bFieldSetting::String = "ThreeD"

# B) typical choices for simulating red-MOT
# simulationType = "SrFRedMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [1, 2, 3, 5, 7, 9, 11, 14, 17]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [.05, .1, .2, .4, .6, 1, 1.5, 2, 2.5, 3, 3.5, 5, 6.5, 8]
# velDirRelToR::String = "Same"
# forceProfile::String = "ThreeD"
# bGradReal::Float64 = 12.5
# bFieldSetting::String = "ThreeD"

# C) typical choices for quick checks of red-det sub-doppler heating magnitude (ideally not too large) + magnitude of ~20 m/s de-celeration (should be high for red MOT)
# simulationType = "SrFRedMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [0.1]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [0.05, 1]
# velDirRelToR::String = "Same"
# forceProfile::String = "ThreeD"
# bGradReal::Float64 = 12.5
# bFieldSetting::String = "ThreeD"

# D) typical choices for quick checks of sub-doppler trasnverse cooling
# simulationType = "SrFTransCooling"
# numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [0.1]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [0.2, 0.5, 1.0, 1.5]
# velDirRelToR::String = "Same"
# forceProfile::String = "TwoD"
# bGradReal::Float64 = 12.5
# bFieldSetting::String = "Static"

# E) typical choices for simulating blue-MOT
# simulationType = "SrFblueMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [.04, .07, .1, .15, .2, .25, .3, .4, .5, .6, .8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3]
# velDirRelToR::String = "Same"
# forceProfile::String = "ThreeD"
# bGradReal::Float64 = 12.5
# bFieldSetting::String = "ThreeD"

# F) typical choices for simulating slowing
# simulationType = "SrFSlowing"
# numTrialsPerValueSet::Int64 = 100
# displacementsInMM::Vector{Float64} = [0.01]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [35] # longitudinal velocity from source, for 2D MOT/tranverse cooling 
# #longSpeeds::Vector{Float64} = vcat([-20, -15, -10, -5, -3, -1, -.5, .5, 1, 2], (3:1:51)) # longitudinal speed for 2d force profile (normalized units v/(gam/k))
# userSpeeds::Vector{Float64} = [0.4]
# velDirRelToR::String = "Same"
# forceProfile::String = "TwoD"
# bGradReal::Float64 = 12.5
# bFieldSetting::String = "Static"