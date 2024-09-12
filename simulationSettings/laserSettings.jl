"""
s0: Saturation parameter corresponding the peak laser intensity of a single laser pass. I_Sat ~ 3 mW/cm^2 for XA and ~ 4 mW/cm^2 for XB for SrF.

laserEnergy: in unit of Gamma, laser detuning measured from the lowest ground hyperfine level.

polType: lsaer polarization and propagation direction. Options are:
    '3D': 6 laser beams propagate in +/-x, +/-y, +/-z 6 directions respectively. 
            All lasers are circularly polarized, but x and y lasers have opposite polarization handedness from z lasers, as required by 3D MOT quadrupole B field. 
            Sigma +/- polarization is indicated by polSign. This is mainly used to simulate 3D MOT or molasses cooling.
    '2DSS': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. 
            All lasers are circularly polarized, but x lasers have opposite handedness from y lasers, as required by 2D MOT B field. 
            Sigma +/- polarization is indicated by polSign. This is mainly used to simulate 2D MOT.
    '2DPar': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. 
            All lasers are linear polarized along z axis. This is mainly used to simulate transverse cooling of a molecule beam, assuming it travels along z axis.
    '2DPerp': Same as '2DPar', except x lasers are polarized along y, and y lasers are polarized along z axis.
    'Slower': 1 laser beam propagate along -z direction, and is linearly polarized along x axis.
    'Push': Same as 'Slower' except the laser propagates along +z direction.
        
polSign: -/+ determine sigma-/+ for polType is 3D or 2DSS. For other 'polType' these are unused.

sidebandFreqs: in unit of Gamma, modulation frequencies of EO modulators.

sidebandAmps: in unit of radians, modulation depths of EO modulators.

whichTransition: can be "XA" (couples X,v=0 to A,v=0), "XB" (couples X,v=0 to B,v=0), and "XARepump" (couples X,v=1 to A,v=0).
    if there is no "XARepump", vibrational branching IS TURNED OFF (obviously, or else all population would accumulate in v=1).
    "CouplingMatrices" and "laser masks" populate based on what values of "whichTransition" are chosen. 

beamWaistInMM: mm, laser beam waist. Used to simulate the effects of the finite size of 3D MOT laser beams. 
    This is only used if polType is '3D', otherwise the laser beams are assumed infinitely large.
"""


# A) blue XB 2D/3D MOT params (note: forceProfile and bFieldSetting should both be "ThreeD" here)
# s0 = [20., 20., 20., 20.]
# detunings = [3, 3, 3, 3]
# laserEnergy = -stateEnergiesGround .+ detunings
# polType = ["2DSS", "2DSS", "2DSS", "2DSS"]
# polSign = [-1, 1, 1, -1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XB", "XB", "XB", "XB"]
# beamWaistInMM = 7.0

# B) X->A transverse cooling params
# s0 = [20., 20., 20., 20., 800.]
# tcDetuning = 3
# laserEnergy = vcat(-mol.stateEnergiesGround .+ tcDetuning, -47)
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# polSign = [1, 1, 1, 1, 1] # doesn't matter here
# sidebandFreqs = [0., 0., 0., 0., 0.6]
# sidebandAmps = [0., 0., 0., 0., 44.]
# whichTransition = ["XA", "XA", "XA", "XA", "XB"]
# beamWaistInMM = 7.0

# C) X->A transverse cooling params with repump
# s0 = [20., 20., 20., 20., 20., 20., 20., 20., 800.] # last laser is slowing laser
# tcDetuning = 3
# laserEnergy = vcat(-mol.stateEnergiesGround .+ tcDetuning, -mol.stateEnergiesGround, -47)
# polType = ["2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "Slower"]
# polSign = [1, 1, 1, 1, 1, 1, 1, 1, 1] # doesn't matter here
# sidebandFreqs = [0., 0., 0., 0., 0., 0., 0., 0., 0.6]
# sidebandAmps = [0., 0., 0., 0., 0., 0., 0., 0., 44.]
# whichTransition = ["XA", "XA", "XA", "XA", "XARepump", "XARepump", "XARepump", "XARepump", "XB"]
# beamWaistInMM = 7.0

# D) red XA 3D 5-laser MOT params
s0 = [10.4, 19.2, 10.4, 31.3, 8.7]
laserEnergy = [-1.0, -9.8, -18.6, -26.8, -20.8]
polType = ["3D", "3D", "3D", "3D", "3D"]
polSign = [1, 1, 1, -1, -1]
sidebandFreqs = [0., 0., 0., 0., 0.]
sidebandAmps = [0., 0., 0., 0., 0.]
whichTransition = ["XA", "XA", "XA", "XA", "XA"]
beamWaistInMM = 7.0

# D2) blue XA only 1 fiber eom
# s0 = [30., 30.]
# sideBandDriveFreq = 9.5
# carrierFreq = -6.6
# laserEnergy = [carrierFreq, -22.0]
# polType = ["3D", "3D"]
# polSign = [1, -1]
# sidebandFreqs = [sideBandDriveFreq, 0.]
# sidebandAmps = [1.9, 0.]
# whichTransition = ["XA", "XA"]
# beamWaistInMM = 7.0

# D2) blue XA Only, both fiber eom (V1)
# s0 = [20., 20.]
# sideBandDriveFreq1 = 19.8
# sideBandDriveFreq2 = 13.5
# carrierFreq1 = 2.5
# carrierFreq2 = -22.
# laserEnergy = [carrierFreq1, carrierFreq2]
# polType = ["3D", "3D"]
# polSign = [1, -1]
# sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps = [1.6, 0.8]
# whichTransition = ["XA", "XA"]
# beamWaistInMM = 7.0

# D2) blue XA Only, both fiber eom (V3)
# s0 = [20., 20.]
# sideBandDriveFreq1 = 25.9
# sideBandDriveFreq2 = 9.8
# carrierFreq1 = 7.5
# carrierFreq2 = -22.
# laserEnergy = [carrierFreq1, carrierFreq2]
# polType = ["3D", "3D"]
# polSign = [1, -1]
# sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2]
# sidebandAmps = [1.6, 1.8]
# whichTransition = ["XA", "XA"]
# beamWaistInMM = 7.0

# D3) blue XA Only, try to find best single freq
# s0 = [30., 2., 8., 6., 4.]
# laserEnergy = [4., -7.5, -16.6, -21.6, -24.4]
# polSign = [1, 1, 1, -1, 1]
# polType = ["3D", "3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0., 0.]
# whichTransition = ["XA", "XA", "XA", "XA", "XA"]
# beamWaistInMM = 7.0

# D4) blue XA real
# s0 = [36., 7., 7., 2.]
# laserEnergy = [-21.2, -17.4, +1.1, -8.5]
# polSign = [-1, 1, 1, -1]
# polType = ["3D", "3D", "3D", "3D"]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XA", "XA", "XA"]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for CaF
# s0 = [1.3, 3.5, 1.2, .0]
# laserEnergy = [+1.1, -13.7, -16.6, -16.7]
# polType = ["3D", "3D", "3D", "3D"]
# polSign = [1, -1, -1, -1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XA", "XA", "XA"]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for SrF F=2
# s0 = [5., 10., 2., 2.]
# laserEnergy = [+1.5, -19.6+2.5, -25.9+1.1, -25.9+1.0]
# polType = ["3D", "3D", "3D", "3D"]
# polSign = [1, -1, -1, 1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XA", "XA", "XA"]
# beamWaistInMM = 7.0

# D5) blueXA in CaOH Style for SrF F=1Down
# s0 = [5., 5., 2., 2.]
# laserEnergy = [-19.5+1.5, -25.9+2.5, +1.1, +1.0]
# polType = ["3D", "3D", "3D", "3D"]
# polSign = [1, -1, -1, 1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XA", "XA", "XA"]
# beamWaistInMM = 7.0

# D6) red XA 2D/3D 4-laser MOT params
# s0 = [30., 30.] ./ 5
# laserEnergy = [-2., -26.9]
# polType = ["3D", "3D"]
# polSign = [1, -1]
# sidebandFreqs = [14.6, 5.3]
# sidebandAmps = [1.4, 1.4]
# whichTransition = ["XA", "XA"]
# beamWaistInMM = 7.0

# E) Bichromatic MOT
# s0 = [30., 10., 30., 45.] ./ 50
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [-0.5, -8.5, -25.1, -29.9]
# polType = ["3D", "3D", "3D", "3D"]
# polSign = [1, -1, 1, -1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XB", "XA", "XB"]
# beamWaistInMM = 7.0

# F) Bichromatic (or not?) Blue MOT
# s0 = [3., 1., 1., 4.] .* 4
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [2, -6.5, -18.6, -22.9]
# polType = ["3D", "3D", "3D", "3D"]
# polSign = [-1, -1, 1, -1]
# sidebandFreqs = [0., 0., 0., 0.]
# sidebandAmps = [0., 0., 0., 0.]
# whichTransition = ["XA", "XB", "XA", "XB"]
# beamWaistInMM = 7.0

# F2) Lambda MOT
# s0 = [7., 23.] .* 1
# # detunings = [-2, -9.5, -23.6, -27.9]
# # laserEnergy = -mol.stateEnergiesGround .+ detunings
# laserEnergy = [2.0, -24.0]
# polType = ["3D", "3D"]
# polSign = [1, 1]
# sidebandFreqs = [0., 0.]
# sidebandAmps = [0., 0.]
# whichTransition = ["XB", "XB"]
# beamWaistInMM = 7.0

# G) Slowing with push
# s0 = [280., 35.] .* 1.0
# laserEnergy = [-40., -13.]
# polType = ["Slower", "Push"]
# polSign = [1, 1] # doesn't matter here
# sidebandFreqs = [0.6, 6.5]
# sidebandAmps = [44., 2.5]
# whichTransition = ["XB", "XA"]
# beamWaistInMM = 7.0

# H) Slowing without push
# s0 = [280.]
# laserEnergy = [-40.]
# polType = ["Slower"]
# polSign = [1] # doesn't matter here
# sidebandFreqs = [0.6]
# sidebandAmps = [44.]
# whichTransition = ["XB"]
# beamWaistInMM = 7.0