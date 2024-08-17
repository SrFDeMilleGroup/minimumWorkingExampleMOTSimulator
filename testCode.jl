# 1) Go to directory and load external variables + functions
cd(@__DIR__); # moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well
include("moleculeVariables/SrFVariables.jl") # change this to whatever molecule you care about
include("auxFunctions/auxFunctions.jl"); # supplementary functions


# 2) User choices with respect to saving output files
saveInRealUnits = 1; # if 1, save vel+accel in m/s, mm/ms^2.  If 0, save in normalized units (vel= v/(gam/k)), (force=1e-3*hbar*k*gam)
saveData = 1; # if you want to save the data
saveDataFolderTag = "SrFRedMOTNormalValues"; # If you want to put anything additional in "savefoldername" to tag it, see variable folderString after lasers are declared.
addHeaders = 1;


# 3) Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
bGradReal = 12.5; # in units Gauss/cm.  if bFieldSetting == "Static", this becomes the static field, in Gauss
waistInMM = 7; # only used if polType is 3D.  Handles finite MOT beam waists
numTrialsPerValueSet = 2; # number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds)
velDirRelToR = 0; # -1 will choose random values for direction of v, r. 0 will force them to be same direction. 1 forces orthogonal.  2 forces opposite.
forceXY = 1; # if '1', will force v, r to go along (x+y)/sqrt(2). Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = 0. If 2, force go along Z
if velDirRelToR == -1
    runType = "Random"; # goes into folder name.
elseif velDirRelToR == 0
    runType = "SameDir"; # goes into folder name.
elseif velDirRelToR == 1
    runType = "OrthoDir"; # goes into folder name.
elseif velDirRelToR == 2
    runType = "OppDir"; # goes into folder name.
else
    throw("Invalid velDirRelToR value: $velDirRelToR. Valid values are -1, 0, 1, 2.");
end
vRound = 0.02; # in unit of Gamma/k, nearest normalized unit to round all speeds to. Simulation will have periodicity of 2*pi/vRound, so if you round finely, will take a while.  Usually choose 0.02


# 4) User choices for what displacements from either origin (if 3D) or else z-axis (if 2D) and speeds
#displacements will always be in millimeters, user speeds always in normalized units Gamma/k (for SrF, 1 normalized unit is roughly 4.4 m/s)

# 4A) parameters for quick test of restoring force
displacementsInMM = [0.5, 1.5, 3.0, 4.5, 6.0, 7.5]; # mm
userSpeeds = [-4, -3, -2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2, 3, 4]; # speed in xy plane (for 2d force profile) or in 3D (in unit of Gamma/k)
forceProfile = "ThreeD"; # either "ThreeD", (forces calculated are (f \dot r) / |r|, (f \dot v) / |v|), or "TwoD" (f \dot (rx,ry,0) / |(rx,ry,0)|, f \dot (vx,vy,0) / |(vx,vy,0)|, and fz are all calculated)
bFieldSetting = "ThreeD"; # can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT") or static "Static" (2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).  


# 4B) typical choices for simulating red-MOT
#=
longSpeeds = [32]; # doesn't matter for 3D sims, sets vel to 140 m/s
displacementsInMM = [1, 2, 3, 5, 7, 9, 11, 14, 17];
userSpeeds = [.05, .1, .2, .4, .6, 1, 1.5, 2, 2.5, 3, 3.5, 5, 6.5, 8];
forceProfile = "ThreeD";
bFieldSetting = "ThreeD";
=#

# 4C) typical choices for quick checks of red-det sub-doppler heating magnitude (ideally not too large) + magnitude of ~20 m/s de-celeration (should be high for red MOT)
#=
longSpeeds = [32]; # doesn't matter for 3D sims, sets vel to 140 m/s
displacementsInMM = [0.1];
userSpeeds = [0.05, 1];
forceProfile = "ThreeD";
bFieldSetting = "ThreeD";
=#

# 4D) typical choices for quick checks of sub-doppler trasnverse cooling
#=
longSpeeds = [32]; # sets vel to 140 m/s (e.g., longitudinal velocity from source) for 2D MOT/tranverse cooling 
displacementsInMM = [0.1];
userSpeeds = [0.2, 0.5, 1.0, 1.5];
forceProfile = "TwoD";
bFieldSetting = "Static";
=#

# 4E) typical choices for simulating blue-MOT
#=
longSpeeds = [32]; # doesn't matter for 3D sims, sets vel to 140 m/s
displacementsInMM = [.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6];
userSpeeds = [.04, .07, .1, .15, .2, .25, .3, .4, .5, .6, .8, 1, 1.2, 1.4, 1.6, 2, 2.5, 3];
forceProfile = "ThreeD";
bFieldSetting = "ThreeD";
=#

# 4F) typical choices for simulating slowing
#=
longSpeeds = [35]; # longitudinal velocity from source, for 2D MOT/tranverse cooling 
#longSpeeds = vcat([-20, -15, -10, -5, -3, -1, -.5, .5, 1, 2], (3:1:51)); # longitudinal speed for 2d force profile (normalized units v/(gam/k))
userSpeeds = [0.4];
displacementsInMM = 0.01;
forceProfile = "TwoD";
bFieldSetting = "Static";
=#

if forceProfile == "ThreeD"
    headers = ["Speed" "av" "del_av" "ar"  "del_ar" "PF1Down" "PF0" "PF1Up" "PF2" "PExc"];
elseif forceProfile == "TwoD"
    headers = ["Speed" "av" "del_av" "ar"  "del_ar"  "LongSpeed" "az" "del_az" "PF1Down" "PF0" "PF1Up" "PF2" "PExc"];
else 
    throw("Invalid forceProfile value: $forceProfile. Valid values are 'ThreeD' or 'TwoD'.");
end


# 5) User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).

#= 
polType can be "3D" (sig +/-, with z-axis (quadrupole coil axis) reversed wrt other axes), "2DSS" (sig +/- but lasers only in x,y direction.  if \sig+ along +x then \sig- along +y).  
"2DPar"(lasers in x,y direction both polarized along z).  "2DPerp" (x laser polarized along y, y polarized along z).  "Slower" (z laser linearly polarized along x) 
=#

#= 
whichTransition can be "XA" (couples X,v=0 to A,v=0), "XB" (couples X,v=0 to B,v=0), and "XARepump" (couples X,v=1 to A,v=0).
if there is no "XARepump", vibrational branching IS TURNED OFF (obviously, or else all population would accumulate in v=1).
"CouplingMatrices" and "laser masks" populate based on what values of "whichTransition" are chosen. 
=#

# 5A) blue XB 2D/3D MOT params (note: forceProfile and bFieldSetting should both be "ThreeD" here)
#=
s0 = [20., 20., 20., 20.]; # single laser pass peak saturation parameter.  I_Sat ~ 3 mW/cm^2 for XA and ~ 4 mW/cm^2 for XB
detunings = [3, 3, 3, 3];
laserEnergy = -stateEnergiesGround .+ detunings; # relative to the energy difference E_{e}-E_{g}
polSign = [-1, 1, 1, -1]; # -/+ determine sigma-/+. for other 'polType' these are unused
whichTransition = ["XB", "XB", "XB", "XB"];
polType = ["2DSS", "2DSS", "2DSS", "2DSS"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5B) X->A transverse cooling params
#=
s0 = [20., 20., 20., 20., 800.]; # last laser is slowing laser
tcDetuning = 3;
laserEnergy = vcat(-stateEnergiesGround .+ tcDetuning, -47);
polSign = [1, 1, 1, 1, 1]; # doesn't matter here
whichTransition = ["XA", "XA", "XA", "XA", "XB"];
polType = ["2DPar", "2DPar", "2DPar", "2DPar", "Slower"];
sidebandFreqs = [0., 0., 0., 0., 0.6];
sidebandAmps = [0., 0., 0., 0., 44.];
=#

# 5C) X->A transverse cooling params with repump
#=
s0 = [20., 20., 20., 20., 20., 20., 20., 20., 800.]; # last laser is slowing laser
tcDetuning = 3;
laserEnergy = vcat(-stateEnergiesGround .+ tcDetuning, -stateEnergiesGround, -47);
polSign = [1, 1, 1, 1, 1, 1, 1, 1, 1]; # doesn't matter here
whichTransition = ["XA", "XA", "XA", "XA", "XARepump", "XARepump", "XARepump", "XARepump", "XB"];
polType = ["2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "2DPar", "Slower"];
sidebandFreqs = [0., 0., 0., 0., 0., 0., 0., 0., 0.6];
sidebandAmps = [0., 0., 0., 0., 0., 0., 0., 0., 44.];
=#

# 5D) red XA 3D 5-laser MOT params

s0 = [10.4, 19.2, 10.4, 31.3, 8.7];
detunings=[0, 0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [-1.0, -9.8, -18.6, -26.8, -20.8];
polSign = [1, 1, 1, -1, -1];
whichTransition = ["XA", "XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0., 0.];


# 5D2) blue XA only 1 fiber eom
#=
s0 = [30., 30.];
detunings=[0, 0]; # not used here, just write actual laser energies
sideBandDriveFreq = 9.5;
carrierFreq = -6.6;
laserEnergy = [carrierFreq, -22.0];
polSign = [1, -1];
whichTransition = ["XA", "XA"];
polType = ["3D", "3D"];
sidebandFreqs = [sideBandDriveFreq, 0.];
sidebandAmps = [1.9, s0.];
=#

# 5D2) blue XA Only, both fiber eom (V1)
#=
s0 = [20., 20.];
detunings=[0, 0]; # not used here, just write actual laser energies
sideBandDriveFreq1 = 19.8;
sideBandDriveFreq2 = 13.5;
carrierFreq1 = 2.5;
carrierFreq2 = -22.;
laserEnergy = [carrierFreq1, carrierFreq2];
polSign = [1, -1];
whichTransition = ["XA", "XA"];
polType = ["3D", "3D"];
sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2];
sidebandAmps = [1.6, 0.8];
=#

# 5D2) blue XA Only, both fiber eom (V3)
#=
s0 = [20., 20.];
detunings=[0, 0]; # not used here, just write actual laser energies
sideBandDriveFreq1 = 25.9;
sideBandDriveFreq2 = 9.8;
carrierFreq1 = 7.5;
carrierFreq2 = -22.;
laserEnergy = [carrierFreq1, carrierFreq2];
polSign = [1, -1];
whichTransition = ["XA", "XA"];
polType = ["3D", "3D"];
sidebandFreqs = [sideBandDriveFreq1, sideBandDriveFreq2];
sidebandAmps = [1.6, 1.8];
=#

# 5D3) blue XA Only, try to find best single freq
#=
s0 = [30., 2., 8., 6., 4.];
detunings=[0, 0]; # not used here, just write actual laser energies
laserEnergy = [4., -7.5, -16.6, -21.6, -24.4];
polSign = [1, 1, 1, -1, 1];
whichTransition = ["XA", "XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0., 0.];
=#

# 5D4) blue XA real
#=
s0 = [36., 7., 7., 2.];
detunings=[0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [-21.2, -17.4, +1.1, -8.5];
polSign = [-1, 1, 1, -1];
whichTransition = ["XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5D5) blueXA in CaOH Style for CaF
#=
s0 = [1.3, 3.5, 1.2, .0];
detunings=[0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [+1.1, -13.7, -16.6, -16.7];
polSign = [1, -1, -1, -1];
whichTransition = ["XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5D5) blueXA in CaOH Style for SrF F=2
#=
s0 = [5., 10., 2., 2.];
detunings=[0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [+1.5, -19.6+2.5, -25.9+1.1, -25.9+1.0];
polSign = [1, -1, -1, 1];
whichTransition = ["XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5D5) blueXA in CaOH Style for SrF F=1Down
#=
s0 = [5., 5., 2., 2.];
detunings=[0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [-19.5+1.5, -25.9+2.5, +1.1, +1.0];
polSign = [1, -1, -1, 1];
whichTransition = ["XA", "XA", "XA", "XA"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5D6) red XA 2D/3D 4-laser MOT params
#=
s0 = [30., 30.] ./ 5;
detunings=[0, 0, 0, 0, 0]; # not used here, just write actual laser energies
laserEnergy = [-2., -26.9];
polSign = [1, -1];
whichTransition = ["XA", "XA"];
polType = ["3D", "3D"];
sidebandFreqs = [14.6, 5.3];
sidebandAmps = [1.4, 1.4];
=#

# 5E) Bichromatic MOT
#=
s0 = [30., 10., 30., 45.] ./ 50;
# detunings = [-2, -9.5, -23.6, -27.9];
# laserEnergy = -stateEnergiesGround .+ detunings; # relative to the energy difference E_{e}-E_{g}
laserEnergy = [-0.5, -8.5, -25.1, -29.9];
polSign = [1, -1, 1, -1];
whichTransition = ["XA", "XB", "XA", "XB"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5F) Bichromatic (or not?) Blue MOT
#=
s0 = [3., 1., 1., 4.] .* 4;
# detunings = [-2, -9.5, -23.6, -27.9];
# laserEnergy = -stateEnergiesGround .+ detunings; # relative to the energy difference E_{e}-E_{g}
laserEnergy = [2, -6.5, -18.6, -22.9];
polSign = [-1, -1, 1, -1];
whichTransition = ["XA", "XB", "XA", "XB"];
polType = ["3D", "3D", "3D", "3D"];
sidebandFreqs = [0., 0., 0., 0.];
sidebandAmps = [0., 0., 0., 0.];
=#

# 5F2) Lambda MOT
#=
s0 = [7., 23.] .* 1;
# detunings = [-2, -9.5, -23.6, -27.9];
# laserEnergy = -stateEnergiesGround .+ detunings; # relative to the energy difference E_{e}-E_{g}
laserEnergy = [2.0, -24.0];
polSign = [1, 1];
whichTransition = ["XB", "XB"];
polType = ["3D", "3D"];
sidebandFreqs = [0., 0.];
sidebandAmps = [0., 0.];
=#

# 5G) Slowing with push
#=
s0 = [280., 35.] .* 1.0;
laserEnergy = [-40., -13.];
polSign = [1, 1]; # doesn't matter here
whichTransition = ["XB", "XA"];
polType = ["Slower", "Push"];
sidebandFreqs = [0.6, 6.5];
sidebandAmps = [44., 2.5];
=#

#5H) Slowing without push
#=
s0 = [280.];
laserEnergy = [-40.];
polSign = [1]; # doesn't matter here
whichTransition = ["XB"];
polType = ["Slower"];
sidebandFreqs = [0.6]; # units of \Gamma
sidebandAmps = [44.]; # radians
=#

checkErrors(bFieldSetting,forceProfile,whichTransition,polType);
#E nd of User Adjustable Values


#6) Stuff for setting up simulation based on user's choices

# stuff needed to determine minimum number of states, and which coupling terms to use, and which lasers actually 'use' a given coupling term (see 'laserMasks')
(couplingMatrices,bCouplingMatrices,stateEnergyMatrix,laserMasks,wavenumberRatios,numZeemanStatesGround,numZeemanStatesExcited) = createCouplingTermsandLaserMasks(whichTransition)

lasers = Lasers(s0,laserEnergy,polSign,whichTransition,polType,sidebandFreqs,sidebandAmps,wavenumberRatios,laserMasks); # define lasers structure, see auxFunctions

# set bGrad (units Gauss * wavevector) (or make "bGrad" static, in units Gauss)
if bFieldSetting != "Static"
    bGrad = (1 / kA * 1e2) * bGradReal; # prefactor converts Gauss/cm to Gauss * wavevector (or Gauss if "Static")
else
    bGrad = bGradReal;
end

numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited;

waist = waistInMM * 1e-3 * kA; # convert to unit of 1/k, waist only used in 3D MOT code

# everything below here is part of the initizliation block that you probably won't want to change

rInit = [0., 0., 0.]; # placehold not used
vInit = [0., 0., 0.]; # placehold not used
# note: p will include a lot of pre-allocated stuff.  This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below
# in retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. (actually you can't. Julia forces ODEProblem to have a variable 'p')
pPreInitialized = preInitializer(length(s0),numZeemanStatesGround,numZeemanStatesTotal)


p = [rInit, vInit, stateEnergyMatrix, lasers, waist, bGrad * normalizedBohrMag,
    couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3]];
append!(p,pPreInitialized);
push!(p,bFieldSetting);

# density matrix Initialization
pStart = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
pStart[1:12, 1:12] = Matrix(I, 12, 12) ./ 12; # molecules equally populate X N=1 states

# NOTE, when executing in VSCode using alt-enter, code will stop here.  To actually run simulation, hit alt-enter again with cursor below the double #
##

#saveFolder
if bFieldSetting == "Static"
    bString="BFieldGauss"
else
    bString="BGradGPerCM"
end
if saveData==1
    folderString = string(@__DIR__, "\\saveData\\", saveDataFolderTag, "bFieldSetting", bFieldSetting, bString, bGradReal, "Force", forceProfile, "NumLasers", length(s0), "Date", Dates.format(now(),"yyyymmdd_HHMM"))
    mkpath(folderString)
end


# OK, that's the setup, now for actually obtaining some acceleration curves via our OBE solver (note: the bulk of the work is 'under the hood' in auxFunctions)

# Simulation: force vs speed for various displacements
# 7) initialize a bunch of different storage variables
forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2);
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a \dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #a \dot r/|r|
if forceProfile == "TwoD"
    forceVsLong = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2); #az
end
pExcVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pF1DownVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pF0VsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pF1UpVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
pF2VsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2);
#initialize some 'masks' that zero out subset of population values...helpful for quick calculation of populations in various ground states
maskExc = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskExc[(numZeemanStatesGround+1):(numZeemanStatesTotal),(numZeemanStatesGround+1):(numZeemanStatesTotal)] .= ones(numZeemanStatesExcited, numZeemanStatesExcited);
maskF1Down = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskF1Down[1:3, 1:3] .= ones(3, 3);
maskF0 = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskF0[4, 4] = 1;
maskF1Up = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskF1Up[5:7, 5:7] .= ones(3, 3);
maskF2 = zeros(numZeemanStatesTotal, numZeemanStatesTotal);
maskF2[8:12, 8:12] .= ones(5, 5);

#8) Iterate over user choices for displacements and speeds

for l = 1:length(displacementsInMM)
    currDisp = displacementsInMM[l];
    for k=1:length(longSpeeds)
        currLongSpeed = longSpeeds[k];
        for j = 1:length(userSpeeds)
            currSpeed = userSpeeds[j];
            if abs(currSpeed)<0.04
                vRound = 0.002
            elseif abs(currSpeed)<0.1
                vRound=0.01;
            elseif abs(currSpeed)<0.5
                vRound = 0.02;
            else
                vRound = 0.05;
            end
            #8A) Set up and solve OBEs
            (randRxs,randRys,randRzs,randVxs,randVys,randVzs) = generateRandPosAndVel(forceProfile,numTrialsPerValueSet,velDirRelToR,currDisp,currSpeed,vRound,currLongSpeed,forceXY);
            tForSteadyState = maximum([10 / currSpeed, 270]); # obtained by trial and error. Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
            periodLength = 2 * pi / vRound;
            saveTimes = tForSteadyState:0.1:(tForSteadyState+periodLength) # times to record obe solution for force integration
            for i = 1:(numTrialsPerValueSet*2)
                forceVsTime[i] = zeros(length(saveTimes), 3)
            end
            prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + periodLength), p)#set up OBE problem to solve

            function prob_func(prob, i, repeat)#change position and velocity of funtion each 'ensemble' sample based on the randomly chosen values for pos/vel vector (magnitude fixed)
                prob.p[1][1] = randRxs[i]
                prob.p[1][2] = randRys[i]
                prob.p[1][3] = randRzs[i]
                prob.p[2][1] = randVxs[i]
                prob.p[2][2] = randVys[i]
                prob.p[2][3] = randVzs[i]
                remake(prob)
            end

            #these two lines here actually handle the parallized runs of the ode solver
            ens_prob = EnsembleProblem(prob, prob_func=prob_func)#solve obe problem for various initial conditions re-set by 'prob_func' each iteration
            @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=numTrialsPerValueSet * 2, saveat=saveTimes)#parallelized OBE solver, runs on amount of threads made available by CPU (Threads.nthreads())
            
            #8B) calculate forces (f\dot r/|r|, etc.) for each random R, V trial..
            @time for i = 1:(numTrialsPerValueSet*2)
                currSol = sol[i]
                makeForceVsTime!(forceVsTime[i], currSol.t, currSol.u, lasers,
                couplingMatrices, stateEnergyMatrix, waist, [randRxs[i], randRys[i], randRzs[i]], [randVxs[i], randVys[i], randVzs[i]])
                if forceProfile=="TwoD"
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    forceVsLong[j,i] = trapz(currSol.t, forceVsTime[i][:, 3]) / 1e-3 / (currSol.t[end] - currSol.t[1]);
                else
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                end
                pExcVsSpeed[j, i] = mean(real(tr.([maskExc .* v for v in currSol.u])))
                pF1DownVsSpeed[j, i] = mean(real(tr.([maskF1Down .* v for v in currSol.u])))
                pF0VsSpeed[j, i] = mean(real(tr.([maskF0 .* v for v in currSol.u])))
                pF1UpVsSpeed[j, i] = mean(real(tr.([maskF1Up .* v for v in currSol.u])))
                pF2VsSpeed[j, i] = mean(real(tr.([maskF2 .* v for v in currSol.u])))
            end#for all trials
        end#for speeds

        #8C) for given set of speeds, for current choices of longSpeed and displacement, average a\dot v, a\dot r, populations,etc. over runs
        forceVsSpeedAvg = mean(forceVsSpeed, dims=2)
        forceVsSpeedAvg = dropdims(forceVsSpeedAvg, dims=(2))#converts to vector
        forceVsSpeedUnc = std(forceVsSpeed, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
        forceVsSpeedUnc = dropdims(forceVsSpeedUnc, dims=(2))

        forceVsPosAvg = mean(forceVsPos, dims=2)
        forceVsPosAvg = dropdims(forceVsPosAvg, dims=(2))
        forceVsPosUnc = std(forceVsPos, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
        forceVsPosUnc = dropdims(forceVsPosUnc, dims=(2))
        if forceProfile=="TwoD"
            forceVsLongAvg = mean(forceVsLong, dims=2)
            forceVsLongAvg = dropdims(forceVsLongAvg, dims=(2))
            forceVsLongUnc = std(forceVsLong, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
            forceVsLongUnc = dropdims(forceVsLongUnc, dims=(2))
        end

        pExcVsSpeedAvg = mean(pExcVsSpeed, dims=2)
        pExcVsSpeedAvg = dropdims(pExcVsSpeedAvg, dims=(2))
        pF1DownVsSpeedAvg = mean(pF1DownVsSpeed, dims=2)
        pF1DownVsSpeedAvg = dropdims(pF1DownVsSpeedAvg, dims=(2))
        pF0VsSpeedAvg = mean(pF0VsSpeed, dims=2)
        pF0VsSpeedAvg = dropdims(pF0VsSpeedAvg, dims=(2))
        pF1UpVsSpeedAvg = mean(pF1UpVsSpeed, dims=2)
        pF1UpVsSpeedAvg = dropdims(pF1UpVsSpeedAvg, dims=(2))
        pF2VsSpeedAvg = mean(pF2VsSpeed, dims=2)
        pF2VsSpeedAvg = dropdims(pF2VsSpeedAvg, dims=(2))

        #8D) convert to real units if applicable and save data
        (forceVsSpeedAvgSaveVals,forceVsSpeedUncSaveVals,forceVsPosAvgSaveVals,forceVsPosUncSaveVals) = (forceVsSpeedAvg,forceVsSpeedUnc,forceVsPosAvg,forceVsPosUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
        userSpeedsSaveVals = userSpeeds.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
         if forceProfile=="TwoD"
            (forceVsLongAvgSaveVals,forceVsLongUncSaveVals) = (forceVsLongAvg,forceVsLongUnc).*(accelFactor*saveInRealUnits+1*(1-saveInRealUnits));
            currLongSpeedSaveVals = currLongSpeed.*(velFactor*saveInRealUnits+1*(1-saveInRealUnits));
        end
        if saveData==1
            open(string(folderString,"\\forceVsSpeedDisplacement",displacementsInMM[l],"MM",runType,".dat"),"a") do io
                if addHeaders==1 && k==1
                    if forceProfile=="TwoD"
                        writedlm(io,[headers ; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pF1DownVsSpeedAvg,pF0VsSpeedAvg,pF1UpVsSpeedAvg,pF2VsSpeedAvg,pExcVsSpeedAvg)]);
                    else
                        writedlm(io,[headers ; hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pF1DownVsSpeedAvg,pF0VsSpeedAvg,pF1UpVsSpeedAvg,pF2VsSpeedAvg, pExcVsSpeedAvg)]);
                    end
                else #if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                    if forceProfile=="TwoD"
                        writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,fill(currLongSpeedSaveVals,length(userSpeeds)),forceVsLongAvgSaveVals,forceVsLongUncSaveVals,pF1DownVsSpeedAvg,pF0VsSpeedAvg,pF1UpVsSpeedAvg,pF2VsSpeedAvg,pExcVsSpeedAvg));
                    else
                        writedlm(io,hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals,pF1DownVsSpeedAvg,pF0VsSpeedAvg,pF1UpVsSpeedAvg,pF2VsSpeedAvg, pExcVsSpeedAvg));
                    end
                end
            end
        end
    end#for longitudinal speeds

end#for displacements

laserVarHeaders = ["s0" "energy" "polSign" "whichTransition" "polType" "sidebandFreqs" "sidebandAmps"]
if saveData ==1
    open(string(folderString,"\\laserVariables.dat"),"w") do io
        writedlm(io,[laserVarHeaders ; hcat(s0,laserEnergy,polSign,whichTransition,polType,sidebandFreqs,sidebandAmps)]);
    end
end

