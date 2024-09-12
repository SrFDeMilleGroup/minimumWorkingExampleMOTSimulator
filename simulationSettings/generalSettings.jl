"""
simulationType: string, notes for yourself about simulation types, e.g., redMOT, blueMOT, transCooling, etc.

numTrialsPerValueSet: number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds).

forceProfile: Options are:
    'ThreeD': Used for 3D MOT and molasses simulations. 
                In this case, 'displacementsInMM' and 'userSpeeds' are treated as magnitudes of 3D displacements and velocities respectively. 'longSpeeds' is ignored. 
                The direction of initial displacements is indicated by 'initDispDir'. And the direction of velocity is set with respect to displacements by 'velDirRelToR'. 
                In order to account for laser field variation, the actual initial displacements in computation are randomly sampled from a cube of size of one wavelength round the values specified here. 
                The projections of 3D force on initial displacement f \\dot r / |r| and f \\dot v / |v| are calculated and returned. 
    'TwoD': Used for slowing and transverse cooling simulations. 
            In this case, 'displacementsInMM' and 'userSpeeds' are treated as magnitudes of 2D displacements and velocities in x-y plane. 
            Displacement along z is taken as zero (up to a random value within +/- 1/2 wavelength from 0). 'longSpeeds' is used as velocity along z direction. 
            'initDispDir' is ignored, and the direction of initial displacements in x-y plane is always random. The relative direction of velocity and displacement in x-y plane is set by 'velDirRelToR'. 
            The projection of 3D force on x-y plane displacement f \\dot r / |r| and f \\dot v / |v|, as well as its z component f_z are calculated and returned.

displacementsInMM: mm, magnitude of initial displacements. Simulation will iterate through the entire list.

initDispDir: Direction of initial displacememnt. Only used if forceProfile is 'ThreeD'. Valid options are 'XY' ((x+y)/sqrt(2) direction, where slowed molecules come into MOT region for most experiments), 'Z' and 'Random'.

longSpeeds: longitudinal speeds in units of Gamma/k (=4.4m/s for SrF), Only used if 'foreProfile' is 'TwoD'.

userSpeeds: in unit of Gamma/k (=4.4m/s for SrF), speeds in xy plane (for 2d force profile) or in 3D.

velDirRelToR: relative direction of molecule velocity, w.r.t. initial displacement r. Options are ["Same", "Orthogonal", "Opposite", "Random"].

bFieldSetting: can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT"),
or static "StaticXY" (in (x+y)/sqrt(2) direction) or "StaticZ" (in z direction). 
(Static B fields are used for 2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).

bGradReal: Gauss/cm for B field gradient, or Gauss for uniform B field, depending on 'bFieldSetting'.
"""


# A) parameters for quick test of restoring force
simulationType::String = "SrFRedMOTNormalValues"
numTrialsPerValueSet::Int64 = 2
forceProfile::String = "ThreeD"
# displacementsInMM::Vector{Float64} = [0.5, 1.5, 3.0, 4.5, 6.0, 7.5]
displacementsInMM::Vector{Float64} = [0.5, 3.0]
initDispDir::String = "XY"
longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [-4, -3, -2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2, 3, 4]
userSpeeds::Vector{Float64} = [0.1, 0.5]
velDirRelToR::String = "Same"
bFieldSetting::String = "ThreeD"
bGradReal::Float64 = 12.5

# B) typical choices for simulating red-MOT
# simulationType::String = "SrFRedMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# forceProfile::String = "ThreeD"
# displacementsInMM::Vector{Float64} = [1, 2, 3, 5, 7, 9, 11, 14, 17]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [.05, .1, .2, .4, .6, 1, 1.5, 2, 2.5, 3, 3.5, 5, 6.5, 8]
# velDirRelToR::String = "Same"
# bFieldSetting::String = "ThreeD"
# bGradReal::Float64 = 12.5

# C) typical choices for quick checks of red-det sub-doppler heating magnitude (ideally not too large) + magnitude of ~20 m/s de-celeration (should be high for red MOT)
# simulationType::String = "SrFRedMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# forceProfile::String = "ThreeD"
# displacementsInMM::Vector{Float64} = [0.1]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [0.05, 1]
# velDirRelToR::String = "Same"
# bFieldSetting::String = "ThreeD"
# bGradReal::Float64 = 12.5

# D) typical choices for quick checks of sub-doppler trasnverse cooling
# simulationType::String = "SrFTransCooling"
# numTrialsPerValueSet::Int64 = 100
# forceProfile::String = "TwoD"
# displacementsInMM::Vector{Float64} = [0.1]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [0.2, 0.5, 1.0, 1.5]
# velDirRelToR::String = "Same"
# bFieldSetting::String = "StaticXY"
# bGradReal::Float64 = 12.5

# E) typical choices for simulating blue-MOT
# simulationType::String = "SrFblueMOTNormalValues"
# numTrialsPerValueSet::Int64 = 100
# forceProfile::String = "ThreeD"
# displacementsInMM::Vector{Float64} = [.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [32]
# userSpeeds::Vector{Float64} = [.04, .07, .1, .15, .2, .25, .3, .4, .5, .6, .8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3]
# velDirRelToR::String = "Same"
# bFieldSetting::String = "ThreeD"
# bGradReal::Float64 = 12.5

# F) typical choices for simulating slowing
# simulationType::String = "SrFSlowing"
# numTrialsPerValueSet::Int64 = 100
# forceProfile::String = "TwoD"
# displacementsInMM::Vector{Float64} = [0.01]
# initDispDir::String = "XY"
# longSpeeds::Vector{Float64} = [35] # longitudinal velocity from source, for 2D MOT/tranverse cooling 
# #longSpeeds::Vector{Float64} = vcat([-20, -15, -10, -5, -3, -1, -.5, .5, 1, 2], (3:1:51)) # longitudinal speed for 2d force profile (normalized units v/(gam/k))
# userSpeeds::Vector{Float64} = [0.4]
# velDirRelToR::String = "Same"
# bFieldSetting::String = "StaticXY"
# bGradReal::Float64 = 12.5