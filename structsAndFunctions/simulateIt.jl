module simulateIt

    """
    unit of energy: hbar * Gamma
    unit of velocity: Gamma / k, where k is the wavevector
    unit of time: 1 / Gamma
    unit of length: 1 / k (= wavelength / 2pi)
    unit of force: 1e-3 * hbar * Gamma * k (??)
    """

    using DifferentialEquations: ODEProblem, EnsembleProblem, solve, remake, Tsit5, EnsembleThreads
    using BenchmarkTools: @time
    using LinearAlgebra: tr, I
    using Statistics: mean, std
    using Trapz: trapz
    using Accessors: @reset

    using ..structs: Molecule, Lasers, GeneralSettings

    include("auxFunctions/obeInitialization.jl") # supplementary functions
    using .obeInitialization: createCouplingTermsandLaserMasks, preInitializer, generateRandPosAndVel

    include("auxFunctions/obeEvaluation.jl") # supplementary functions
    using .obeEvaluation: densityMatrixChangeTerms!

    include("auxFunctions/forceCalculation.jl") # supplementary functions
    using .forceCalculation: makeForceVsTime!


    export simulateOBE
    

    function simulateOBE(mol::Molecule, lasers::Lasers, general::GeneralSettings, currDisp::Float64, currSpeed::Float64, currLongSpeed::Float64)

        ## OBE initialization ##

        # Set vRound and freqRound for velocity and laser frequency/molecule energy to round to
        # The maximum common divisor of vRound and freqRound, say valRound, will set the hamiltonian period to be 2*pi/valRound
        if abs(currSpeed) < 0.04
            vRound = 0.002
        elseif abs(currSpeed) < 0.1
            vRound = 0.01
        elseif abs(currSpeed) < 0.5
            vRound = 0.02
        else
            vRound = 0.05
        end

        (randRxs, randRys, randRzs, randVxs, randVys, randVzs) = generateRandPosAndVel(general, currDisp, currSpeed, vRound, currLongSpeed, mol)

        freqRound = 0.1
        @reset mol.stateEnergiesGround = round.(mol.stateEnergiesGround ./ freqRound) .* freqRound
        @reset mol.stateEnergiesExcited = round.(mol.stateEnergiesExcited ./ freqRound) .* freqRound
        @reset lasers.laserEnergy = round.(lasers.laserEnergy ./ freqRound) .* freqRound
        @reset lasers.sidebandFreqs = round.(lasers.sidebandFreqs ./ freqRound) .* freqRound

        # stuff needed to determine minimum number of states, and which coupling terms to use
        (couplingMatrices, bCouplingMatrices, stateEnergyMatrix, numZeemanStatesGround, numZeemanStatesExcited) = createCouplingTermsandLaserMasks(lasers, mol)
        numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited

        rInit = [0., 0., 0.] # initial position, placehold not used
        vel = [0., 0., 0.] # velocity, placehold not used

        # note: p will include a lot of pre-allocated stuff. 
        # This is basically all of the stuff 'passed to' the obe solver, in addition to the initial condition of the density matrix defined below.
        # In retrospect p is not the best choice for the variable name but it's the julia house style...maybe replace later. 
        # (actually you can't. Julia forces ODEProblem to have a variable 'p')
        pPreInitialized = preInitializer(lasers, numZeemanStatesGround, numZeemanStatesTotal)

        p = [rInit, vel, stateEnergyMatrix, lasers, general, mol,
            couplingMatrices[1], couplingMatrices[2], couplingMatrices[3], bCouplingMatrices[1], bCouplingMatrices[2], bCouplingMatrices[3]]
        append!(p, pPreInitialized)
        p = Tuple(p)

        # initial value of density matrix
        pStart = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pStart[1:12, 1:12] .= Matrix(I, 12, 12) ./ 12 # molecules equally populate X N=1 states


        ## OBE evaluation ##

        # obtained by trial and error. Could potentially be handled more rigrorously (solve ode in steps of 'period length' until solution 'converges')
        tForSteadyState = maximum([10 / currSpeed, 270])
        hamiltonianPeriod = 2 * pi / vRound
        saveTimes = tForSteadyState : 0.1 : (tForSteadyState + hamiltonianPeriod) # times to record obe solution for force integration

        # set up OBE problem to solve
        prob = ODEProblem(densityMatrixChangeTerms!, pStart, (0.0, tForSteadyState + hamiltonianPeriod), p)

        # change initial position and velocity of each 'ensemble' sample based on the randomly chosen values for pos/vel vector (magnitude fixed)
        function prob_func(prob, i, repeat)
            prob.p[1][1] = randRxs[i]
            prob.p[1][2] = randRys[i]
            prob.p[1][3] = randRzs[i]
            prob.p[2][1] = randVxs[i]
            prob.p[2][2] = randVys[i]
            prob.p[2][3] = randVzs[i]

            # makes a copy of the problem and return it, so different threads sees different problems
            return remake(prob)
        end

        # define obe problem for various initial conditions re-set by 'prob_func' each iteration
        ens_prob = EnsembleProblem(prob, prob_func=prob_func) 

        # solve ensemble problem concurrently, runs on amount of threads made available by the current CPU core (Threads.nthreads())
        @time sol = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=general.numTrialsPerValueSet * 2, saveat=saveTimes)
        

        ## Force calculation ##

        forceVsTime = zeros(length(saveTimes), 3) # 3 is for x, y, z three different directions
        forceProjOnVel = Vector{Float64}(undef, general.numTrialsPerValueSet * 2) # a \dot v/|v|
        forceProjOnPos = Vector{Float64}(undef, general.numTrialsPerValueSet * 2) # a \dot r/|r|
        forceProjOnLong = zeros(Float64, general.numTrialsPerValueSet * 2) # a_z, won't be used if forceProfile is "ThreeD", so initialize to zero rather than uninitialized

        # population on each hyperfine level
        pExc = Vector{Float64}(undef, general.numTrialsPerValueSet * 2)
        pF1Down = Vector{Float64}(undef, general.numTrialsPerValueSet * 2)
        pF0 = Vector{Float64}(undef, general.numTrialsPerValueSet * 2)
        pF1Up = Vector{Float64}(undef, general.numTrialsPerValueSet * 2)
        pF2 = Vector{Float64}(undef, general.numTrialsPerValueSet * 2)

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

        @time for i = 1 : (general.numTrialsPerValueSet*2)
            currSol = sol[i]

            # Note that this function re-write the forceVsTime array everytime
            makeForceVsTime!(forceVsTime, currSol.t, currSol.u, lasers, couplingMatrices, stateEnergyMatrix, [randRxs[i], randRys[i], randRzs[i]], [randVxs[i], randVys[i], randVzs[i]])

            if general.forceProfile == "TwoD"
                # forceProjOnVel = (f \dot v) / |v| with time average
                forceProjOnVel[i] = (randVxs[i] * trapz(currSol.t, forceVsTime[:, 1]) + randVys[i] * trapz(currSol.t, forceVsTime[:, 2])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                # forceProjOnPos = (f \dot r) / |r| with time average
                forceProjOnPos[i] = (randRxs[i] * trapz(currSol.t, forceVsTime[:, 1]) + randRys[i] * trapz(currSol.t, forceVsTime[:, 2])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2) / (currSol.t[end] - currSol.t[1])
                # forceProjOnLong = fz with time average
                forceProjOnLong[i] = trapz(currSol.t, forceVsTime[:, 3]) / 1e-3 / (currSol.t[end] - currSol.t[1])
            elseif general.forceProfile == "ThreeD"
                # forceProjOnVel = (f \dot v) / |v| with time average 
                forceProjOnVel[i] = (randVxs[i] * trapz(currSol.t, forceVsTime[:, 1]) + randVys[i] * trapz(currSol.t, forceVsTime[:, 2]) +
                randVzs[i] * trapz(currSol.t, forceVsTime[:, 3])) / 1e-3 / sqrt(randVxs[i] .^ 2 + randVys[i] .^ 2 + randVzs[i] .^2) / (currSol.t[end] - currSol.t[1])
                # forceProjOnPos = (f \dot r) / |r| with time average
                forceProjOnPos[i] = (randRxs[i] * trapz(currSol.t, forceVsTime[:, 1]) + randRys[i] * trapz(currSol.t, forceVsTime[:, 2]) +
                randRzs[i] * trapz(currSol.t, forceVsTime[:, 3])) / 1e-3 / sqrt(randRxs[i] .^ 2 + randRys[i] .^ 2 + randRzs[i] .^2) / (currSol.t[end] - currSol.t[1])
            else
                error("Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'.")
            end

            pExc[i] = mean(real(tr.([maskExc .* v for v in currSol.u])))
            pF1Down[i] = mean(real(tr.([maskF1Down .* v for v in currSol.u])))
            pF0[i] = mean(real(tr.([maskF0 .* v for v in currSol.u])))
            pF1Up[i] = mean(real(tr.([maskF1Up .* v for v in currSol.u])))
            pF2[i] = mean(real(tr.([maskF2 .* v for v in currSol.u])))
        end

        forceProjOnVelAvg = mean(forceProjOnVel)
        forceProjOnVelUnc = std(forceProjOnVel) / sqrt(general.numTrialsPerValueSet * 2)

        forceProjOnPosAvg = mean(forceProjOnPos)
        forceProjOnPosUnc = std(forceProjOnPos) / sqrt(general.numTrialsPerValueSet * 2)

        # forceProjOnLong is just zero if forceProfile is "ThreeD"
        forceProjOnLongAvg = mean(forceProjOnLong)
        forceProjOnLongUnc = std(forceProjOnLong) / sqrt(general.numTrialsPerValueSet * 2)

        pExcAvg = mean(pExc)
        pF1DownAvg = mean(pF1Down)
        pF0Avg = mean(pF0)
        pF1UpAvg = mean(pF1Up)
        pF2Avg = mean(pF2)

        return (forceProjOnVelAvg, forceProjOnVelUnc, forceProjOnPosAvg, forceProjOnPosUnc, forceProjOnLongAvg, forceProjOnLongUnc, pExcAvg, pF1DownAvg, pF0Avg, pF1UpAvg, pF2Avg)
    end
end