module laserSettings

export s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps

"""
s0: single laser pass peak saturation parameter.  I_Sat ~ 3 mW/cm^2 for XA and ~ 4 mW/cm^2 for XB

laserEnergy: in unit of Gamma, relative to the energy difference E_{e}-E_{g}

polSign: -/+ determine sigma-/+. for other 'polType' these are unused

whichTransition: can be "XA" (couples X,v=0 to A,v=0), "XB" (couples X,v=0 to B,v=0), and "XARepump" (couples X,v=1 to A,v=0).
    if there is no "XARepump", vibrational branching IS TURNED OFF (obviously, or else all population would accumulate in v=1).
    "CouplingMatrices" and "laser masks" populate based on what values of "whichTransition" are chosen. 

polType: can be "3D" (sig +/-, with z-axis (quadrupole coil axis) reversed wrt other axes), "2DSS" (sig +/- but lasers only in x,y direction. if \\sig+ along +x then \\sig- along +y),  
    "2DPar"(lasers in x,y direction both polarized along z), "2DPerp" (x laser polarized along y, y polarized along z), "Slower" (z laser linearly polarized along x). 

sidebandFreqs: in unit of Gamma

sidebandAmps: in unit of radians

"""

# A) blue XB 2D/3D MOT params (note: forceProfile and bFieldSetting should both be "ThreeD" here)
# s0::Vector{Float64} = [20., 20., 20., 20.]
# detunings = [3, 3, 3, 3]
# laserEnergy::Vector{Float64} = -stateEnergiesGround .+ detunings
# polSign = [-1, 1, 1, -1]
# whichTransition = ["XB", "XB", "XB", "XB"]
# polType = ["2DSS", "2DSS", "2DSS", "2DSS"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# B) X->A transverse cooling params
# s0::Vector{Float64} = [20., 20., 20., 20., 800.]
# tcDetuning = 3
# laserEnergy::Vector{Float64} = vcat(-stateEnergiesGround .+ tcDetuning, -47)
# polSign = [1, 1, 1, 1, 1] # doesn't matter here
# whichTransition = ["XA", "XA", "XA", "XA", "XB"]
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0., 0.6]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0., 44.]

# C) X->A transverse cooling params with repump
# s0::Vector{Float64} = [20., 20., 20., 20., 20., 20., 20., 20., 800.] # last laser is slowing laser
# tcDetuning = 3
# laserEnergy::Vector{Float64} = vcat(-stateEnergiesGround .+ tcDetuning, -stateEnergiesGround, -47)
# polSign = [1, 1, 1, 1, 1, 1, 1, 1, 1] # doesn't matter here
# whichTransition = ["XA", "XA", "XA", "XA", "XARepump", "XARepump", "XARepump", "XARepump", "XB"]
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0., 0., 0., 0., 0., 0.6]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0., 0., 0., 0., 0., 44.]

# D) red XA 3D 5-laser MOT params
s0::Vector{Float64} = [10.4, 19.2, 10.4, 31.3, 8.7]
laserEnergy::Vector{Float64} = [-1.0, -9.8, -18.6, -26.8, -20.8]
polSign::Vector{Int64} = [1, 1, 1, -1, -1]
whichTransition::Vector{String} = ["XA", "XA", "XA", "XA", "XA"]
polType::Vector{String} = ["3D", "3D", "3D", "3D", "3D"]
sidebandFreqs::Vector{Float64} = [0., 0., 0., 0., 0.]
sidebandAmps::Vector{Float64} = [0., 0., 0., 0., 0.]

# D2) blue XA only 1 fiber eom
# s0::Vector{Float64} = [30., 30.]
# sideBandDriveFreq = 9.5
# carrierFreq = -6.6
# laserEnergy::Vector{Float64} = [carrierFreq, -22.0]
# polSign::Vector{Int64} = [1, -1]
# whichTransition::Vector{String} = ["XA", "XA"]
# polType::Vector{String} = ["3D", "3D"]
# sidebandFreqs::Vector{Float64} = [sideBandDriveFreq, 0.]
# sidebandAmps::Vector{Float64} = [1.9, s0.]

# D2) blue XA Only, both fiber eom (V1)
# s0::Vector{Float64} = [20., 20.]
# sideBandDriveFreq1 = 19.8
# sideBandDriveFreq2 = 13.5
# carrierFreq1 = 2.5
# carrierFreq2 = -22.
# laserEnergy::Vector{Float64} = [carrierFreq1, carrierFreq2]
# polSign::Vector{Int64} = [1, -1]
# whichTransition::Vector{String} = ["XA", "XA"]
# polType::Vector{String} = ["3D", "3D"]
# sidebandFreqs::Vector{Float64} = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps::Vector{Float64} = [1.6, 0.8]

# D2) blue XA Only, both fiber eom (V3)
# s0::Vector{Float64} = [20., 20.]
# sideBandDriveFreq1 = 25.9
# sideBandDriveFreq2 = 9.8
# carrierFreq1 = 7.5
# carrierFreq2 = -22.
# laserEnergy::Vector{Float64} = [carrierFreq1, carrierFreq2]
# polSign::Vector{Int64} = [1, -1]
# whichTransition::Vector{String} = ["XA", "XA"]
# polType::Vector{String} = ["3D", "3D"]
# sidebandFreqs::Vector{Float64} = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps::Vector{Float64} = [1.6, 1.8]

# D3) blue XA Only, try to find best single freq
# s0::Vector{Float64} = [30., 2., 8., 6., 4.]
# laserEnergy::Vector{Float64} = [4., -7.5, -16.6, -21.6, -24.4]
# polSign::Vector{Int64} = [1, 1, 1, -1, 1]
# whichTransition::Vector{String} = ["XA", "XA", "XA", "XA", "XA"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0., 0.]

# D4) blue XA real
# s0::Vector{Float64} = [36., 7., 7., 2.]
# laserEnergy::Vector{Float64} = [-21.2, -17.4, +1.1, -8.5]
# polSign::Vector{Int64} = [-1, 1, 1, -1]
# whichTransition::Vector{String} = ["XA", "XA", "XA", "XA"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# D5) blueXA in CaOH Style for CaF
# s0::Vector{Float64} = [1.3, 3.5, 1.2, .0]
# laserEnergy::Vector{Float64} = [+1.1, -13.7, -16.6, -16.7]
# polSign::Vector{Int64} = [1, -1, -1, -1]
# whichTransition::Vector{String} = ["XA", "XA", "XA", "XA"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# D5) blueXA in CaOH Style for SrF F=2
# s0::Vector{Float64} = [5., 10., 2., 2.]
# laserEnergy::Vector{Float64} = [+1.5, -19.6+2.5, -25.9+1.1, -25.9+1.0]
# polSign::Vector{Int64} = [1, -1, -1, 1]
# whichTransition::Vector{String} = ["XA", "XA", "XA", "XA"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# D5) blueXA in CaOH Style for SrF F=1Down
# s0::Vector{Float64} = [5., 5., 2., 2.]
# laserEnergy::Vector{Float64} = [-19.5+1.5, -25.9+2.5, +1.1, +1.0]
# polSign::Vector{Int64} = [1, -1, -1, 1]
# whichTransition::Vector{String} = ["XA", "XA", "XA", "XA"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# D6) red XA 2D/3D 4-laser MOT params
# s0::Vector{Float64} = [30., 30.] ./ 5
# laserEnergy::Vector{Float64} = [-2., -26.9]
# polSign::Vector{Int64} = [1, -1]
# whichTransition::Vector{String} = ["XA", "XA"]
# polType::Vector{String} = ["3D", "3D"]
# sidebandFreqs::Vector{Float64} = [14.6, 5.3]
# sidebandAmps::Vector{Float64} = [1.4, 1.4]

# E) Bichromatic MOT
# s0::Vector{Float64} = [30., 10., 30., 45.] ./ 50
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy::Vector{Float64} = -stateEnergiesGround .+ detunings
# laserEnergy::Vector{Float64} = [-0.5, -8.5, -25.1, -29.9]
# polSign::Vector{Int64} = [1, -1, 1, -1]
# whichTransition::Vector{String} = ["XA", "XB", "XA", "XB"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# F) Bichromatic (or not?) Blue MOT
# s0::Vector{Float64} = [3., 1., 1., 4.] .* 4
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy::Vector{Float64} = -stateEnergiesGround .+ detunings
# laserEnergy::Vector{Float64} = [2, -6.5, -18.6, -22.9]
# polSign::Vector{Int64} = [-1, -1, 1, -1]
# whichTransition::Vector{String} = ["XA", "XB", "XA", "XB"]
# polType::Vector{String} = ["3D", "3D", "3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0., 0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0., 0., 0.]

# F2) Lambda MOT
# s0::Vector{Float64} = [7., 23.] .* 1
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy::Vector{Float64} = -stateEnergiesGround .+ detunings
# laserEnergy::Vector{Float64} = [2.0, -24.0]
# polSign::Vector{Int64} = [1, 1]
# whichTransition::Vector{String} = ["XB", "XB"]
# polType::Vector{String} = ["3D", "3D"]
# sidebandFreqs::Vector{Float64} = [0., 0.]
# sidebandAmps::Vector{Float64} = [0., 0.]

# G) Slowing with push
# s0::Vector{Float64} = [280., 35.] .* 1.0
# laserEnergy::Vector{Float64} = [-40., -13.]
# polSign::Vector{Int64} = [1, 1] # doesn't matter here
# whichTransition::Vector{String} = ["XB", "XA"]
# polType::Vector{String} = ["Slower", "Push"]
# sidebandFreqs::Vector{Float64} = [0.6, 6.5]
# sidebandAmps::Vector{Float64} = [44., 2.5]

# H) Slowing without push
# s0::Vector{Float64} = [280.]
# laserEnergy::Vector{Float64} = [-40.]
# polSign::Vector{Int64} = [1] # doesn't matter here
# whichTransition::Vector{String} = ["XB"]
# polType::Vector{String} = ["Slower"]
# sidebandFreqs::Vector{Float64} = [0.6]
# sidebandAmps::Vector{Float64} = [44.]

@assert length(s0) == length(laserEnergy) == length(polSign) == length(whichTransition) == length(polType) == length(sidebandFreqs) == length(sidebandAmps) "All laser settings arrays must have the same length."
@assert all([i in [1, -1] for i in polSign]) "Invalid polSign value(s): $polSign. All values must be either 1 or -1."
@assert all([i in ["XA", "XB", "XARepump"] for i in whichTransition]) "Invalid whichTransition value(s): $whichTransition. All values must be either 'XA', 'XB', or 'XARepump'."
@assert all([i in ["3D", "2DSS", "2DPar", "2DPerp", "Slower", "Push"] for i in polType]) "Invalid polType value(s): $polType. All values must be one of ['3D', '2DSS', '2DPar', '2DPerp', 'Slower', 'Push']."

end