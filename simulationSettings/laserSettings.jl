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

beamWaistInMM: in unit of mm, only used if polType is 3D. Handles finite MOT beam waists
"""


# A) blue XB 2D/3D MOT params (note: forceProfile and bFieldSetting should both be "ThreeD" here)
# s0 = [20., 20., 20., 20.]
# detunings = [3, 3, 3, 3]
# laserEnergy = -stateEnergiesGround .+ detunings
# polSign = [-1, 1, 1, -1]
# whichTransition = ["XB", "XB", "XB", "XB"]
# polType = ["2DSS", "2DSS", "2DSS", "2DSS"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# B) X->A transverse cooling params
# s0 = [20., 20., 20., 20., 800.]
# tcDetuning = 3
# laserEnergy = vcat(-mol.stateEnergiesGround .+ tcDetuning, -47)
# polSign = [1, 1, 1, 1, 1] # doesn't matter here
# whichTransition = ["XA", "XA", "XA", "XA", "XB"]
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# sidebandFreqs = [0., 0., 0., 0., 0.6]
# sidebandAmps = [0., 0., 0., 0., 44.]
# beamWaistInMM = 7.0

# C) X->A transverse cooling params with repump
# s0 = [20., 20., 20., 20., 20., 20., 20., 20., 800.] # last laser is slowing laser
# tcDetuning = 3
# laserEnergy = vcat(-mol.stateEnergiesGround .+ tcDetuning, -mol.stateEnergiesGround, -47)
# polSign = [1, 1, 1, 1, 1, 1, 1, 1, 1] # doesn't matter here
# whichTransition = ["XA", "XA", "XA", "XA", "XARepump", "XARepump", "XARepump", "XARepump", "XB"]
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# sidebandFreqs = [0., 0., 0., 0., 0., 0., 0., 0., 0.6]
# sidebandAmps = [0., 0., 0., 0., 0., 0., 0., 0., 44.]
# beamWaistInMM = 7.0

# D) red XA 3D 5-laser MOT params
s0 = [10.4, 19.2, 10.4, 31.3, 8.7]
laserEnergy = [-1.0, -9.8, -18.6, -26.8, -20.8]
polSign = [1, 1, 1, -1, -1]
whichTransition = ["XA", "XA", "XA", "XA", "XA"]
polType = ["3D", "3D", "3D", "3D", "3D"]
sidebandFreqs = [0., 0., 0., 0., 0.]
sidebandAmps = [0., 0., 0., 0., 0.]
beamWaistInMM = 7.0

# D2) blue XA only 1 fiber eom
# s0 = [30., 30.]
# sideBandDriveFreq = 9.5
# carrierFreq = -6.6
# laserEnergy = [carrierFreq, -22.0]
# polSign = [1, -1]
# whichTransition = ["XA", "XA"]
# polType = ["3D", "3D"]
# sidebandFreqs = [sideBandDriveFreq, 0.]
# sidebandAmps = [1.9, 0.]
# beamWaistInMM = 7.0

# D2) blue XA Only, both fiber eom (V1)
# s0 = [20., 20.]
# sideBandDriveFreq1 = 19.8
# sideBandDriveFreq2 = 13.5
# carrierFreq1 = 2.5
# carrierFreq2 = -22.
# laserEnergy = [carrierFreq1, carrierFreq2]
# polSign = [1, -1]
# whichTransition = ["XA", "XA"]
# polType = ["3D", "3D"]
# sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps = [1.6, 0.8]
# beamWaistInMM = 7.0

# D2) blue XA Only, both fiber eom (V3)
# s0 = [20., 20.]
# sideBandDriveFreq1 = 25.9
# sideBandDriveFreq2 = 9.8
# carrierFreq1 = 7.5
# carrierFreq2 = -22.
# laserEnergy = [carrierFreq1, carrierFreq2]
# polSign = [1, -1]
# whichTransition = ["XA", "XA"]
# polType = ["3D", "3D"]
# sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps = [1.6, 1.8]
# beamWaistInMM = 7.0

# D3) blue XA Only, try to find best single freq
# s0 = [30., 2., 8., 6., 4.]
# laserEnergy = [4., -7.5, -16.6, -21.6, -24.4]
# polSign = [1, 1, 1, -1, 1]
# whichTransition = ["XA", "XA", "XA", "XA", "XA"]
# polType = ["3D", "3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0., 0.]
# beamWaistInMM = 7.0

# D4) blue XA real
# s0 = [36., 7., 7., 2.]
# laserEnergy = [-21.2, -17.4, +1.1, -8.5]
# polSign = [-1, 1, 1, -1]
# whichTransition = ["XA", "XA", "XA", "XA"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for CaF
# s0 = [1.3, 3.5, 1.2, .0]
# laserEnergy = [+1.1, -13.7, -16.6, -16.7]
# polSign = [1, -1, -1, -1]
# whichTransition = ["XA", "XA", "XA", "XA"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for SrF F=2
# s0 = [5., 10., 2., 2.]
# laserEnergy = [+1.5, -19.6+2.5, -25.9+1.1, -25.9+1.0]
# polSign = [1, -1, -1, 1]
# whichTransition = ["XA", "XA", "XA", "XA"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for SrF F=1Down
# s0 = [5., 5., 2., 2.]
# laserEnergy = [-19.5+1.5, -25.9+2.5, +1.1, +1.0]
# polSign = [1, -1, -1, 1]
# whichTransition = ["XA", "XA", "XA", "XA"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# D6) red XA 2D/3D 4-laser MOT params
# s0 = [30., 30.] ./ 5
# laserEnergy = [-2., -26.9]
# polSign = [1, -1]
# whichTransition = ["XA", "XA"]
# polType = ["3D", "3D"]
# sidebandFreqs = [14.6, 5.3]
# sidebandAmps = [1.4, 1.4]
# beamWaistInMM = 7.0

# E) Bichromatic MOT
# s0 = [30., 10., 30., 45.] ./ 50
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [-0.5, -8.5, -25.1, -29.9]
# polSign = [1, -1, 1, -1]
# whichTransition = ["XA", "XB", "XA", "XB"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# F) Bichromatic (or not?) Blue MOT
# s0 = [3., 1., 1., 4.] .* 4
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [2, -6.5, -18.6, -22.9]
# polSign = [-1, -1, 1, -1]
# whichTransition = ["XA", "XB", "XA", "XB"]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# beamWaistInMM = 7.0

# F2) Lambda MOT
# s0 = [7., 23.] .* 1
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [2.0, -24.0]
# polSign = [1, 1]
# whichTransition = ["XB", "XB"]
# polType = ["3D", "3D"]
# sidebandFreqs = [0., 0.]
# sidebandAmps = [0., 0.]
# beamWaistInMM = 7.0

# G) Slowing with push
# s0 = [280., 35.] .* 1.0
# laserEnergy = [-40., -13.]
# polSign = [1, 1] # doesn't matter here
# whichTransition = ["XB", "XA"]
# polType = ["Slower", "Push"]
# sidebandFreqs = [0.6, 6.5]
# sidebandAmps = [44., 2.5]
# beamWaistInMM = 7.0

# H) Slowing without push
# s0 = [280.]
# laserEnergy = [-40.]
# polSign = [1] # doesn't matter here
# whichTransition = ["XB"]
# polType = ["Slower"]
# sidebandFreqs = [0.6]
# sidebandAmps = [44.]
# beamWaistInMM = 7.0