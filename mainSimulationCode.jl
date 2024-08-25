# unit of energy: hbar * Gamma
# unit of velocity: Gamma / k, where k is the wavevector
# unit of time: 1 / Gamma
# unit of length: 1 / k (= wavelength / 2pi)
# unit of force: 1e-3 * hbar * Gamma * k (??)


using DifferentialEquations: ODEProblem, EnsembleProblem, solve, remake, Tsit5, EnsembleThreads
using BenchmarkTools: @time
using LinearAlgebra: tr, I
using DelimitedFiles: writedlm
using Statistics: mean, std
using Dates: Dates, now
using SharedArrays: SharedArray
using Trapz: trapz


# 1) Go to directory and load external variables + functions
cd(@__DIR__) # moves julia terminal to directory where this file is.  This directory should have auxFunctions+SrF(or whatever)Variables files as well

include("simulationSettings/moleculeVariables.jl")
using .moleculeVariables: SrF, CaF, BaF, MgF, CaOH, SrOH, Molecule
mol = SrF

include("auxFunctions/auxFunctions.jl") # supplementary functions


# 2) User choices with respect to saving output files
include("simulationSettings/saveSettings.jl")
using .saveSettings: saveInRealUnits, saveData, saveDataFolderTag, addHeaders


# 3) Non Laser Detuning/Pol Simulation Variables (B-field, beam-waist etc.)
include("simulationSettings/generalSettingsOne.jl")
using .generalSettingsOne: bGradReal, waistInMM, numTrialsPerValueSet, velDirRelToR, initDispDir


# 4) User choices for what displacements and speeds
include("simulationSettings/generalSettingsTwo.jl")
using .generalSettingsTwo: longSpeeds, displacementsInMM, userSpeeds, forceProfile, bFieldSetting


# 5) User choices for laser parameters (detuning, polarization, etc) example laser values (these all work for SrF).
include("simulationSettings/laserSettings.jl")
using .laserSettings: s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps

 
#6) Stuff for setting up simulation based on user's choices
# stuff needed to determine minimum number of states, and which coupling terms to use, and which lasers actually 'use' a given coupling term (see 'laserMasks')
(couplingMatrices, bCouplingMatrices, stateEnergyMatrix, laserMasks, wavenumberRatios, numZeemanStatesGround, numZeemanStatesExcited) = createCouplingTermsandLaserMasks(whichTransition, mol)
numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited

# define lasers structure, see auxFunctions
lasers = Lasers(s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps, wavenumberRatios, laserMasks)

# set bGrad (units Gauss * wavevector) (or make "bGrad" static, in units Gauss)
bGrad = bFieldSetting == "Static" ? bGradReal : (1 / mol.kA * 1e2) * bGradReal

waist = waistInMM * 1e-3 * mol.kA # convert to unit of 1/k, waist only used in 3D MOT code

rInit = [0., 0., 0.] # placehold not used
vInit = [0., 0., 0.] # placehold not used
# note: p will include a lot of pre-allocated stuff.  This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below
# in retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. (actually you can't. Julia forces ODEProblem to have a variable 'p')
pPreInitialized = preInitializer(length(s0), numZeemanStatesGround, numZeemanStatesTotal)

p = [rInit, vInit, stateEnergyMatrix, lasers, waist, bGrad * mol.normalizedBohrMag,
    couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3]]
append!(p, pPreInitialized)
push!(p, bFieldSetting)

# initial value of density matrix
pStart = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
pStart[1:12, 1:12] = Matrix(I, 12, 12) ./ 12 # molecules equally populate X N=1 states

# initialize a bunch of different storage variables for simulation of force vs speed at various displacements
forceVsTime = Array{Array{ComplexF64,2},1}(undef, numTrialsPerValueSet * 2)
forceVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2) # a \dot v/|v|
forceVsPos = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2) # a \dot r/|r|

if forceProfile == "TwoD"
    forceVsLong = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2) # az
elseif forceProfile == "ThreeD"
else
    error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
end

pExcVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2)
pF1DownVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2)
pF0VsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2)
pF1UpVsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2)
pF2VsSpeed = SharedArray{Float64}(length(userSpeeds), numTrialsPerValueSet * 2)

# initialize some 'masks' that zero out subset of population values...helpful for quick calculation of populations in various ground states
maskExc = zeros(numZeemanStatesTotal, numZeemanStatesTotal)
maskExc[(numZeemanStatesGround+1) : (numZeemanStatesTotal), (numZeemanStatesGround+1) : (numZeemanStatesTotal)] .= ones(numZeemanStatesExcited, numZeemanStatesExcited)
maskF1Down = zeros(numZeemanStatesTotal, numZeemanStatesTotal)
maskF1Down[1:3, 1:3] .= ones(3, 3)
maskF0 = zeros(numZeemanStatesTotal, numZeemanStatesTotal)
maskF0[4, 4] = 1
maskF1Up = zeros(numZeemanStatesTotal, numZeemanStatesTotal)
maskF1Up[5:7, 5:7] .= ones(3, 3)
maskF2 = zeros(numZeemanStatesTotal, numZeemanStatesTotal)
maskF2[8:12, 8:12] .= ones(5, 5)

# NOTE, when executing in VSCode using alt-enter (Windows) or opt-enter (MacOS), code will stop here if the line below the next line ONLY contains two #'s. 
# In this case, hit alt-enter or opt-enter again with cursor below the double #.
# ##


# OK, that's the setup, now for actually obtaining some acceleration curves via our OBE solver (note: the bulk of the work is 'under the hood' in auxFunctions)
# 8) Iterate over user choices for displacements and speeds
if saveData
    bString = bFieldSetting == "Static" ? "BFieldGauss" : "BGradGPerCM"
    folderString = string(@__DIR__, "/savedData/", saveDataFolderTag, "bFieldSetting", bFieldSetting, bString, bGradReal, "Force", forceProfile, "NumLasers", length(s0), "Date", Dates.format(now(),"yyyymmdd_HHMMSS"))
    mkpath(folderString)
end

for currDisp in displacementsInMM
    for (k, currLongSpeed) in enumerate(longSpeeds)
        for (j, currSpeed) in enumerate(userSpeeds)
            if abs(currSpeed) < 0.04
                vRound = 0.002
            elseif abs(currSpeed) < 0.1
                vRound = 0.01
            elseif abs(currSpeed) < 0.5
                vRound = 0.02
            else
                vRound = 0.05
            end

            # 8A) Set up and solve OBEs
            (randRxs, randRys, randRzs, randVxs, randVys, randVzs) = generateRandPosAndVel(forceProfile, numTrialsPerValueSet, velDirRelToR, currDisp, currSpeed, vRound, currLongSpeed, initDispDir, mol)
            tForSteadyState = maximum([10 / currSpeed, 270]) # obtained by trial and error. Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
            hamiltonianPeriod = 2 * pi / vRound
            saveTimes = tForSteadyState : 0.1 : (tForSteadyState + hamiltonianPeriod) # times to record obe solution for force integration
            for i = 1 : (numTrialsPerValueSet * 2)
                forceVsTime[i] = zeros(length(saveTimes), 3) # 3 is for x, y, z three different directions
            end
            prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + hamiltonianPeriod), p)#set up OBE problem to solve

            # change initial position and velocity of each 'ensemble' sample based on the randomly chosen values for pos/vel vector (magnitude fixed)
            function prob_func(prob, i, repeat)
                prob.p[1][1] = randRxs[i]
                prob.p[1][2] = randRys[i]
                prob.p[1][3] = randRzs[i]
                prob.p[2][1] = randVxs[i]
                prob.p[2][2] = randVys[i]
                prob.p[2][3] = randVzs[i]
                return remake(prob) # makes a copy of the problem and return it, so different threads sees different problems
            end

            # these two lines here actually handle the parallized runs of the ode solver
            ens_prob = EnsembleProblem(prob, prob_func=prob_func) # solve obe problem for various initial conditions re-set by 'prob_func' each iteration
            @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=numTrialsPerValueSet * 2, saveat=saveTimes) # parallelized OBE solver, runs on amount of threads made available by CPU (Threads.nthreads())
            
            # 8B) calculate forces (f\dot r/|r|, etc.) for each random R, V trial..
            @time for i = 1 : (numTrialsPerValueSet*2)
                currSol = sol[i]
                makeForceVsTime!(forceVsTime[i], currSol.t, currSol.u, lasers,
                couplingMatrices, stateEnergyMatrix, waist, [randRxs[i], randRys[i], randRzs[i]], [randVxs[i], randVys[i], randVzs[i]])

                if forceProfile == "TwoD"
                    # forceVsSpeed = (f \dot v) / |v| with time average
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    # forceVsPos = (f \dot r) / |r| with time average
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) + randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                    # forceVsLong = fz with time average
                    forceVsLong[j,i] = trapz(currSol.t, forceVsTime[i][:, 3]) / 1e-3 / (currSol.t[end] - currSol.t[1])
                elseif forceProfile == "ThreeD"
                    # forceVsSpeed = (f \dot v) / |v| with time average 
                    forceVsSpeed[j, i] = (randVxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randVys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randVzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                    # forceVsPos = (f \dot r) / |r| with time average
                    forceVsPos[j, i] = (randRxs[i] * trapz(currSol.t, forceVsTime[i][:, 1]) +
                    randRys[i] * trapz(currSol.t, forceVsTime[i][:, 2]) +
                    randRzs[i] * trapz(currSol.t, forceVsTime[i][:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                else
                    error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
                end

                pExcVsSpeed[j, i] = mean(real(tr.([maskExc .* v for v in currSol.u])))
                pF1DownVsSpeed[j, i] = mean(real(tr.([maskF1Down .* v for v in currSol.u])))
                pF0VsSpeed[j, i] = mean(real(tr.([maskF0 .* v for v in currSol.u])))
                pF1UpVsSpeed[j, i] = mean(real(tr.([maskF1Up .* v for v in currSol.u])))
                pF2VsSpeed[j, i] = mean(real(tr.([maskF2 .* v for v in currSol.u])))

            end # for all trials
        end # for speeds

        # 8C) for given set of speeds, for current choices of longSpeed and displacement, average a \dot v, a \dot r, populations,etc. over runs
        forceVsSpeedAvg = mean(forceVsSpeed, dims=2)
        forceVsSpeedAvg = dropdims(forceVsSpeedAvg, dims=(2)) # converts to vector
        forceVsSpeedUnc = std(forceVsSpeed, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
        forceVsSpeedUnc = dropdims(forceVsSpeedUnc, dims=(2))

        forceVsPosAvg = mean(forceVsPos, dims=2)
        forceVsPosAvg = dropdims(forceVsPosAvg, dims=(2))
        forceVsPosUnc = std(forceVsPos, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
        forceVsPosUnc = dropdims(forceVsPosUnc, dims=(2))

        if forceProfile == "TwoD"
            forceVsLongAvg = mean(forceVsLong, dims=2)
            forceVsLongAvg = dropdims(forceVsLongAvg, dims=(2))
            forceVsLongUnc = std(forceVsLong, dims=2) ./ sqrt(numTrialsPerValueSet * 2)
            forceVsLongUnc = dropdims(forceVsLongUnc, dims=(2))
        elseif forceProfile == "ThreeD"
        else
            error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
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

        # 8D) convert to real units if applicable and save data
        (forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals) = (forceVsSpeedAvg, forceVsSpeedUnc, forceVsPosAvg, forceVsPosUnc) .* (saveInRealUnits ? mol.accelFactor : 1)
        userSpeedsSaveVals = userSpeeds .* (saveInRealUnits ? mol.velFactor : 1)
        
        if forceProfile == "TwoD"
            (forceVsLongAvgSaveVals, forceVsLongUncSaveVals) = (forceVsLongAvg,forceVsLongUnc) .* (saveInRealUnits ? mol.accelFactor : 1)
            currLongSpeedSaveVals = currLongSpeed .* (saveInRealUnits ? mol.velFactor : 1)
        elseif forceProfile == "ThreeD"
        else
            error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
        end
        
        if saveData
            open(string(folderString, "/forceVsSpeedDisplacement", currDisp, "MM", velDirRelToR, "Dir", ".dat"), "a") do io
                if addHeaders && k==1
                    if forceProfile == "TwoD"
                        headers = ["Speed" "av" "av_std" "ar"  "ar_std"  "LongSpeed" "az" "az_std" "PF1Down" "PF0" "PF1Up" "PF2" "PExc"]
                        writedlm(io, headers)
                    elseif forceProfile == "ThreeD"
                        headers = ["Speed" "av" "av_std" "ar"  "ar_std" "PF1Down" "PF0" "PF1Up" "PF2" "PExc"]
                        writedlm(io, headers)
                    else
                        throw("Invalid forceProfile value: $forceProfile. Valid values are 'ThreeD' or 'TwoD'.")
                    end
                end
                
                # if you've already added headers/don't want them, just append the current forceVsSpeed to the relevant file (so, if you have different longSpeeds, they'll all show up in same file since file is distinguished by displacement)
                if forceProfile == "TwoD"
                    writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals, fill(currLongSpeedSaveVals,length(userSpeeds)), forceVsLongAvgSaveVals, forceVsLongUncSaveVals, pF1DownVsSpeedAvg, pF0VsSpeedAvg, pF1UpVsSpeedAvg, pF2VsSpeedAvg, pExcVsSpeedAvg))
                elseif forceProfile == "ThreeD"
                    writedlm(io, hcat(userSpeedsSaveVals, forceVsSpeedAvgSaveVals, forceVsSpeedUncSaveVals, forceVsPosAvgSaveVals, forceVsPosUncSaveVals, pF1DownVsSpeedAvg, pF0VsSpeedAvg, pF1UpVsSpeedAvg, pF2VsSpeedAvg, pExcVsSpeedAvg))
                else
                    throw("Invalid forceProfile value: $forceProfile. Valid values are 'ThreeD' or 'TwoD'.")
                end
            end
        end

        # don't loop over longSpeeds since it doesn't matter for "ThreeD" forceProfile
        if forceProfile == "TwoD"
        elseif forceProfile == "ThreeD"
            break 
        else
            error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
        end

    end # for longitudinal speeds

end # for displacements

laserVarHeaders = ["s0" "energy" "polSign" "whichTransition" "polType" "sidebandFreqs" "sidebandAmps"]
if saveData
    open(string(folderString, "/laserVariables.dat"), "w") do io
        writedlm(io, [laserVarHeaders ; hcat(s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps)])
    end
end
