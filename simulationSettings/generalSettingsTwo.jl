module generalSettingsTwo

export longSpeeds, displacementsInMM, userSpeeds, forceProfile, bFieldSetting

"""
longSpeeds: longitudinal speeds in units of Gamma/k (=4.4m/s for SrF), doesn't matter for 3D sims, sets vel to 140 m/s

displacementsInMM: mm, initial displacements from either origin (if 3D) or else z-axis (if 2D) in mm

userSpeeds: in unit of Gamma/k (=4.4m/s for SrF), speeds in xy plane (for 2d force profile) or in 3D

forceProfile: either "ThreeD", (forces calculated are (f \\dot r) / |r|, (f \\dot v) / |v|), 
    or "TwoD" (f \\dot (rx,ry,0) / |(rx,ry,0)|, f \\dot (vx,vy,0) / |(vx,vy,0|, and fz are all calculated)

bFieldSetting: can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT"),
    or static "Static" (2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).
    
"""

# A) parameters for quick test of restoring force
longSpeeds::Vector{Float64} = [32]
displacementsInMM::Vector{Float64} = [0.5, 1.5, 3.0, 4.5, 6.0, 7.5]
userSpeeds::Vector{Float64} = [-4, -3, -2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2, 3, 4]
forceProfile::String = "ThreeD"
bFieldSetting::String = "ThreeD"

# B) typical choices for simulating red-MOT
# longSpeeds::Vector{Float64} = [32]
# displacementsInMM::Vector{Float64} = [1, 2, 3, 5, 7, 9, 11, 14, 17]
# userSpeeds::Vector{Float64} = [.05, .1, .2, .4, .6, 1, 1.5, 2, 2.5, 3, 3.5, 5, 6.5, 8]
# forceProfile::String = "ThreeD"
# bFieldSetting::String = "ThreeD"

# C) typical choices for quick checks of red-det sub-doppler heating magnitude (ideally not too large) + magnitude of ~20 m/s de-celeration (should be high for red MOT)
# longSpeeds::Vector{Float64} = [32]
# displacementsInMM::Vector{Float64} = [0.1]
# userSpeeds::Vector{Float64} = [0.05, 1]
# forceProfile::String = "ThreeD"
# bFieldSetting::String = "ThreeD"

# D) typical choices for quick checks of sub-doppler trasnverse cooling
# longSpeeds::Vector{Float64} = [32]
# displacementsInMM::Vector{Float64} = [0.1]
# userSpeeds::Vector{Float64} = [0.2, 0.5, 1.0, 1.5]
# forceProfile::String = "TwoD"
# bFieldSetting::String = "Static"

# E) typical choices for simulating blue-MOT
# longSpeeds::Vector{Float64} = [32]
# displacementsInMM::Vector{Float64} = [.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]
# userSpeeds::Vector{Float64} = [.04, .07, .1, .15, .2, .25, .3, .4, .5, .6, .8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3]
# forceProfile::String = "ThreeD"
# bFieldSetting::String = "ThreeD"

# F) typical choices for simulating slowing
# longSpeeds::Vector{Float64} = [35] # longitudinal velocity from source, for 2D MOT/tranverse cooling 
# #longSpeeds::Vector{Float64} = vcat([-20, -15, -10, -5, -3, -1, -.5, .5, 1, 2], (3:1:51)) # longitudinal speed for 2d force profile (normalized units v/(gam/k))
# userSpeeds::Vector{Float64} = [0.4]
# displacementsInMM::Vector{Float64} = [0.01]
# forceProfile::String = "TwoD"
# bFieldSetting::String = "Static"

@assert forceProfile in ["ThreeD", "TwoD"] "Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'."
@assert bFieldSetting in ["ThreeD", "TwoD", "Static"] "Invalid bFieldSetting value: $bFieldSetting. It must be one of ['ThreeD', 'TwoD', 'Static']."

end