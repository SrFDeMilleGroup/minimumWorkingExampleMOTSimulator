module obeInitialization

    using Random: Xoshiro

    using ..moleculeVariables: Molecule
    using ..laserSettings: Lasers
    using ..generalSettings: GeneralSettings

    export preInitializer, generateRandPosAndVel, createCouplingTermsandLaserMasks

    # Fix the random number generator seed for reproducibility
    # Don't use the default RNG, as it may be called in other background processes.
    myRNG = Xoshiro(123) 


    function preInitializer(lasers::Lasers, numZeemanStatesGround, numZeemanStatesTotal)
        # initializes a bunch of stuff used in the OBE solver.  Julia likes things pre-initialized if possible

        # holds the modified coupling matrices used in decay terms
        coupleMatEff1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        coupleMatEff2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        coupleMatEff3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)

        # convenient for fast evaluation of terms used in the 'decay' term of the density matrix evolution (second term in eq 1 of main writeup)
        decayMaskAllButTopLeft = zeros(Float64, numZeemanStatesTotal, numZeemanStatesTotal);
        decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1
        decayMaskAllButTopLeft[1:numZeemanStatesGround, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1 / 2
        decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, 1:numZeemanStatesGround] .= -1 / 2
        decayMaskForCalcTopLeft = zeros(Int64, numZeemanStatesTotal, numZeemanStatesTotal)
        decayMaskForCalcTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= 1

        # now we make a bunch of initializations.  This makes the julia code run much faster at the cost of some readability...
        r = Vector{Float64}(undef, 3)

        # fieldTerms[i] are the projections of the light field for laser[i] at a given position on the \sigma^-, \pi, \sigma^+ basis
        fieldTerms = [zeros(ComplexF64, 3) for i in 1:lasers.numLasers]

        # will eventually 'hold' the atom-light matrix term of the hamiltonian during the diff-eq solver (see densityMatrixChangeTerms! in auxFunctions)
        atomLightTerm = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)

        # bField terms are the projections of magnetic field at a given position on the \sigma^-,\pi,\sigma^+ basis. bFieldTermFull basically holds the 'mu' tensor
        bFieldTerms = Vector{ComplexF64}(undef, 3)
        bFieldTermFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)

        # will eventually hold the -\mu \cdot B (and hermitian conjugate) terms
        uProdBField = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        bFieldProdU = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)

        # initializations of some matrices used to speed up the decay term calculation
        decayFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pOnlyExcitedStates = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft1PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft2PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft3PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)
        pTopLeft3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal)

        # return pre-initialized stuff from here to join the rest of the pre-initialized stuff in the 'main' program
        pPreInitialized = [coupleMatEff1, coupleMatEff2, coupleMatEff3, decayMaskAllButTopLeft, 
        decayMaskForCalcTopLeft, r, fieldTerms, atomLightTerm, bFieldTerms, bFieldTermFull, uProdBField, bFieldProdU, decayFull,
        pOnlyExcitedStates, pTopLeft1PreMult, pTopLeft2PreMult, pTopLeft3PreMult, pTopLeft1, pTopLeft2, pTopLeft3]

        return pPreInitialized
    end


    function generateRandPosAndVel(general::GeneralSettings, currDisp, currSpeed, vRound, longSpeed, mol::Molecule)
        # Function generates a set of random positions and 'pseudo'-random velocities (direction determined by 'velDirRelToR' + whether 'force profile' is 2D or 3D.)
        # if forceProfile is TwoD: z position is assumed to not matter, z velocity is fixed to longSpeed, and direction of velocity relative to random choice of \phi where x=disp*(cos(\phi)), etc. determined by velDirRelToR
        # if forceProfile is ThreeD: longSpeed isn't used, and direction of velocity chosen relative to random x,y,z direction of position is determined by velDirRelToR

        forceProfile = general.forceProfile
        velDirRelToR = general.velDirRelToR
        initDispDir = general.initDispDir
        numTrialsPerSpeed = general.numTrialsPerValueSet

        if forceProfile == "TwoD"
            # randomize position direction
            randPhisPos = rand(myRNG, numTrialsPerSpeed) * 2 * pi
            randRxs = cos.(randPhisPos) .* currDisp .* 1e-3 .* mol.kA
            randRys = sin.(randPhisPos) .* currDisp .* 1e-3 .* mol.kA
            randRzs = rand(myRNG, numTrialsPerSpeed) * 2 * pi

            if velDirRelToR == "Random" #randomize phi
                randPhisVels = rand(myRNG, numTrialsPerSpeed) * 2 * pi
            elseif velDirRelToR == "Same"
                randPhisVels = randPhisPos
            elseif velDirRelToR == "Orthogonal"
                randPhisVels = randPhisPos .+ pi/2
            elseif velDirRelToR == "Opposite"
                randPhisVels = randPhisPos .+ pi
            else
                throw(ArgumentError(string("Invalid choice of velDirRelToR, ", velDirRelToR, ". Valid options are Same, Orthogonal, Opposite, Random.")))
            end
            randVxs = round.(currSpeed .* cos.(randPhisVels) ./ vRound) .* vRound
            randVys = round.(currSpeed .* sin.(randPhisVels) ./ vRound) .* vRound
            randVzs = fill(round(longSpeed/vRound)*vRound, numTrialsPerSpeed)

        elseif forceProfile == "ThreeD"
            # if initDispDir=="XY", force position to be along (x+y)/sqrt(2) (e.g., entering from slower); if initDispDir=="Z", then it forces along Z
            if initDispDir == "XY"
                randRxs = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* mol.kA .+ 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)
                randRys = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* mol.kA .+ 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)
                randRzs = 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)            
            elseif initDispDir == "Z"
                randRxs = 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)
                randRys = 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)
                randRzs = currDisp .* 1e-3 .* mol.kA .+ 2 .* pi .* (rand(myRNG, numTrialsPerSpeed) .- 0.5)
            else
                throw(ArgumentError(string("Invalid choice of initDispDir, ", initDispDir, ". Valid options are XY or Z.")))
            end
            normTerms = sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2)
            randX = randRxs ./ normTerms
            randY = randRys ./ normTerms
            randZ = randRzs ./ normTerms

            if velDirRelToR == "Random" # random velocity direction as wel
                randX = rand(myRNG, numTrialsPerSpeed) .- 0.5 # re-roll
                randY = rand(myRNG, numTrialsPerSpeed) .- 0.5
                randZ = rand(myRNG, numTrialsPerSpeed) .- 0.5
                normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2)
                randVxs = randX ./ normTerms .* currSpeed
                randVys = randY ./ normTerms .* currSpeed
                randVzs = randZ ./ normTerms .* currSpeed
            elseif velDirRelToR == "Same"
                randVxs = randX .* currSpeed
                randVys = randY .* currSpeed
                randVzs = randZ .* currSpeed
            elseif velDirRelToR == "Orthogonal"
                randX2 = rand(myRNG, numTrialsPerSpeed) .- 0.5
                randY2 = rand(myRNG, numTrialsPerSpeed) .- 0.5
                randZ2 = rand(myRNG, numTrialsPerSpeed) .- 0.5
                for i = 1:numTrialsPerSpeed
                    (randX2[i],randY2[i],randZ2[i]) =[randX2[i],randY2[i],randZ2[i]] - dot([randX[i],randY[i],randZ[i]],[randX2[i],randY2[i],randZ2[i]]) .* [randX[i],randY[i],randZ[i]]
                end
                normTerms = sqrt.(randX2.^2 .+ randY2.^2 .+ randZ2 .^2)
                randVxs = randX2 ./ normTerms .* currSpeed
                randVys = randY2 ./ normTerms .* currSpeed
                randVzs = randZ2 ./ normTerms .* currSpeed
            elseif velDirRelToR == "Opposite"
                randVxs = -randX .* currSpeed
                randVys = -randY .* currSpeed
                randVzs = -randZ .* currSpeed
            else
                throw(ArgumentError(string("Invalid choice of velDirRelToR, ", velDirRelToR, ". Valid options are Same, Orthogonal, Opposite, Random.")))
            end
            randVxs = round.(randVxs ./ vRound) .* vRound
            randVys = round.(randVys ./ vRound) .* vRound
            randVzs = round.(randVzs ./ vRound) .* vRound
        else
            throw(ArgumentError(string("Invalid choice of forceProfile, ", forceProfile, ". Valid options are TwoD or ThreeD.")))
        end

        # run for both +/- r and +/- v (better statistics)
        randRxs = [randRxs; -randRxs]
        randRys = [randRys; -randRys]
        randRzs = [randRzs; -randRzs]
        randVxs = [randVxs; -randVxs]
        randVys = [randVys; -randVys]
        if forceProfile == "TwoD"
            randVzs = [randVzs; randVzs] # for 2D, Vz is always "longSpeed"
        elseif forceProfile == "ThreeD"
            randVzs = [randVzs; -randVzs]
        else
            throw(ArgumentError(string("Invalid choice of forceProfile, ", forceProfile, ". Valid options are TwoD or ThreeD.")))
        end
        
        # velocity along any dimension cannot be zero (particle should have x,y,z all change throughout OBE evolution to ensure periodicity)
        for i = 1:numTrialsPerSpeed*2
            if randVxs[i] == 0
                randVxs[i] = vRound * sign(rand(myRNG) - 0.5)
            end
            if randVys[i] == 0
                randVys[i] = vRound * sign(rand(myRNG) - 0.5)
            end
            if randVzs[i] == 0
                randVzs[i] = vRound * sign(rand(myRNG) - 0.5)
            end
        end
        
        return randRxs, randRys, randRzs, randVxs, randVys, randVzs
    end


    function createCouplingTermsandLaserMasks(lasers::Lasers, mol::Molecule)
        # This function does a number of things

        # 1) determine how many ground and excited states are needed (12 ground if no lasers are "XARepumps", 24 if there are repumps. 4 excited if only one of "A" or "B" are used, 8 if both are)

        # 2) Based on this, write out "stateEnergyMatrix". Ultimately this is subtracted from the laser energy in the OBE solver exp(-i*t*(energyDiff)) like term. All columns are identical.  
        # 2 (cont)) each row (i) is the energy of |i> relative to |F=1,J=1/2> (if i is a ground state) or |E,F'=1> for |i> corresponding to either E=A\Pi or E=B\Sigma.

        # 3) Establish 'coupling' (C matrices, eq 12-14 of writeup, basically 'clebsch-gordan' like terms) and 'b-coupling' matrices (C_B matrices, eq 27-29 of writeup.  Basically a 'B-field' coupling matrix based on g_{F} terms)
        # 3 (cont)) Terms C_{i,j} are zero unless i=ground and j=excited.  size of matrix determined by number of excited states and ground states needed.  C_{i,j}[k] is the coupling from i->j for polarization k
        # 3 (cont)) Terms C_{B,i,j} are zero unless i and j are in same F manifold.  size of matrix determined by number of excited states and ground states needed.  C_{B,i,j}[k] is the coupling from i->j for <B\cdot p_{k}>/|B|, where p_{k} is the \sigma^-/+,\pi basis

        #1)
        bichrom = (("XA" in lasers.whichTransition) && ("XB" in lasers.whichTransition)) ? 1 : 0 # winds up 0 if only XA of XB are used, 1 if both are
        XToB = (!("XA" in lasers.whichTransition) && ("XB" in lasers.whichTransition)) ? 1 : 0 # winds up 1 if only lasers are XB
        repump = ("XARepump" in lasers.whichTransition) ? 1 : 0 # winds up 0 if no repump, 1 if there are repumps

        numZeemanStatesGround = 12 + 12 * repump
        numZeemanStatesExcited = 4 + 4 * bichrom
        numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited

        #2)
        # fill(groundStateEnergy, numZeemanStatesEachHyperfineLevel)
        stateEnergiesColumnFormat = [fill(mol.stateEnergiesGround[1], 3); fill(mol.stateEnergiesGround[2], 1); fill(mol.stateEnergiesGround[3], 3); fill(mol.stateEnergiesGround[4], 5)]
        if repump == 1
            # NOTE this assumes hyperfine splitting is the same in v=1 repump...not quite right but close enough
            stateEnergiesColumnFormat = vcat(stateEnergiesColumnFormat, stateEnergiesColumnFormat)
        end
        stateEnergiesColumnFormat = [stateEnergiesColumnFormat; fill(0, 4+4*bichrom)]

        stateEnergyMatrix = repeat(stateEnergiesColumnFormat, 1, numZeemanStatesTotal)
        if XToB == 1 || bichrom == 1
            stateEnergyMatrix[1:numZeemanStatesGround, end] = stateEnergyMatrix[1:numZeemanStatesGround, end] .- mol.stateEnergiesExcited[2] # handles excited state hyperfine splitting of |B\Sigma,F=0> level. 
            if bichrom == 1
                stateEnergyMatrix[1:numZeemanStatesGround, end-4] = stateEnergyMatrix[1:numZeemanStatesGround, end-4] .- mol.stateEnergiesExcited[1] # handles excited state hyperfine splitting of |A\Pi,F=0> level. 
            end
        else
            stateEnergyMatrix[1:numZeemanStatesGround, end] = stateEnergyMatrix[1:numZeemanStatesGround, end] .- mol.stateEnergiesExcited[1] # handles excited state hyperfine splitting of |A\Pi,F=0> level. 
        end

        #3)
        couplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)]

        makeCouplingMatrices!(couplingMatrices, XToB, repump, bichrom, mol)

        bCouplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)]

        makeBCouplingMatrices!(bCouplingMatrices, XToB, repump, bichrom, mol)

        return couplingMatrices, bCouplingMatrices, stateEnergyMatrix, numZeemanStatesGround, numZeemanStatesExcited
    end


    function makeCouplingMatrices!(couplingMatrices, XToB, repump, bichrom, mol::Molecule)
        # makes C_{i,j}[k] matrices.  What these look like depend on what ground/excited states are included
        # Choice 1) bichrom means that both A and B are 'spoken' to, and thus there are 8 excited states.
        # Choice 2) XToB=0 is true if no lasers 'talk' to B.  Thus, all 4 excited states are A states
        # Choice 3) Thus, if bichrom=0 and XToB=1, all 4 excited states are B states
        # In all cases, repump can be added, and thus there are 12 ground states. These can only be coupled to the A state

        # these are all hardcoded for the assumption of a SrF, CaF, etc. type molecule where the alkaline has no hyperfine structure and, in the X\Sigma,N=1 state there
        # is mixing between 'pure' |F=1,J=1/2> and |F=1,J=3/2> that can be parameterized by a,b where |F=1,J~3/2> = a|F=1,J=3/2>+b|F=1,J=1/2> and |F=1,J~1/2> = -b|F=1,J=3/2>+a|F=1,J=1/2>
        # See Appendix A in writeup

        a = mol.jMixingRatioA
        b = mol.jMixingRatioB

        if bichrom == 1 
            # note: 12*repump term in second index forces 'excited' index to start at appropriate place, e.g. 13 for no repump, 25 if there is repump
            couplingMatrices[1][1, 14+12*repump] = -sqrt(2) / 3 * a - b / 6
            couplingMatrices[1][1, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[1][2, 15+12*repump] = -sqrt(2) / 3 * a - b / 6
            couplingMatrices[1][4, 15+12*repump] = sqrt(2) / 3
            couplingMatrices[1][5, 14+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[1][5, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[1][6, 15+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[1][8, 13+12*repump] = -1 / sqrt(6)
            couplingMatrices[1][9, 14+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[1][10, 15+12*repump] = -1 / 6
            couplingMatrices[1][1, 14+4+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[1][1, 16+4+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[1][2, 15+4+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[1][4, 15+4+12*repump] = 1 / 3
            couplingMatrices[1][5, 14+4+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[1][5, 16+4+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[1][6, 15+4+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[1][8, 13+4+12*repump] = 1 / sqrt(3)
            couplingMatrices[1][9, 14+4+12*repump] = 1 / sqrt(6)
            couplingMatrices[1][10, 15+4+12*repump] = 1 / 3 / sqrt(2)

            couplingMatrices[2][1, 13+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[2][2, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[2][3, 15+12*repump] = -sqrt(2) / 3 * a - 1 / 6 * b
            couplingMatrices[2][4, 14+12*repump] = -sqrt(2) / 3
            couplingMatrices[2][5, 13+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[2][6, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[2][7, 15+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[2][9, 13+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[2][10, 14+12*repump] = -1 / 3
            couplingMatrices[2][11, 15+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[2][1, 13+4+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[2][2, 16+4+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[2][3, 15+4+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[2][4, 14+4+12*repump] = -1 / 3
            couplingMatrices[2][5, 13+4+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[2][6, 16+4+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[2][7, 15+4+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[2][9, 13+4+12*repump] = 1 / sqrt(6)
            couplingMatrices[2][10, 14+4+12*repump] = sqrt(2) / 3
            couplingMatrices[2][11, 15+4+12*repump] = 1 / sqrt(6)

            couplingMatrices[3][2, 13+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[3][3, 14+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[3][3, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[3][4, 13+12*repump] = sqrt(2) / 3
            couplingMatrices[3][6, 13+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[3][7, 14+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[3][7, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[3][10, 13+12*repump] = -1 / 6
            couplingMatrices[3][11, 14+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[3][12, 15+12*repump] = -1 / sqrt(6)
            couplingMatrices[3][2, 13+4+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[3][3, 14+4+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[3][3, 16+4+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[3][4, 13+4+12*repump] = 1 / 3
            couplingMatrices[3][6, 13+4+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[3][7, 14+4+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[3][7, 16+4+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[3][10, 13+4+12*repump] = 1 / 3 / sqrt(2)
            couplingMatrices[3][11, 14+4+12*repump] = 1 / sqrt(6)
            couplingMatrices[3][12, 15+4+12*repump] = 1 / sqrt(3)

            if repump == 1
                couplingMatrices[1][13:24, 25:28] = couplingMatrices[1][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)
                couplingMatrices[1][13:24, 29:32] = couplingMatrices[1][1:12, 29:32] .* sqrt(mol.v1BranchingRatioB)
                couplingMatrices[2][13:24, 25:28] = couplingMatrices[2][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)
                couplingMatrices[2][13:24, 29:32] = couplingMatrices[2][1:12, 29:32] .* sqrt(mol.v1BranchingRatioB)
                couplingMatrices[3][13:24, 25:28] = couplingMatrices[3][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)
                couplingMatrices[3][13:24, 29:32] = couplingMatrices[3][1:12, 29:32] .* sqrt(mol.v1BranchingRatioB)

                couplingMatrices[1][1:12, 25:28] = couplingMatrices[1][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
                couplingMatrices[1][1:12, 29:32] = couplingMatrices[1][1:12, 29:32] .* sqrt(1-mol.v1BranchingRatioB)
                couplingMatrices[2][1:12, 25:28] = couplingMatrices[2][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
                couplingMatrices[2][1:12, 29:32] = couplingMatrices[2][1:12, 29:32] .* sqrt(1-mol.v1BranchingRatioB)
                couplingMatrices[3][1:12, 25:28] = couplingMatrices[3][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
                couplingMatrices[3][1:12, 29:32] = couplingMatrices[3][1:12, 29:32] .* sqrt(1-mol.v1BranchingRatioB)
            end

        elseif XToB == 0 
            # excited states are all "A" states
            couplingMatrices[1][1, 14+12*repump] = -sqrt(2) / 3 * a - b / 6
            couplingMatrices[1][1, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[1][2, 15+12*repump] = -sqrt(2) / 3 * a - b / 6
            couplingMatrices[1][4, 15+12*repump] = sqrt(2) / 3
            couplingMatrices[1][5, 14+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[1][5, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[1][6, 15+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[1][8, 13+12*repump] = -1 / sqrt(6)
            couplingMatrices[1][9, 14+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[1][10, 15+12*repump] = -1 / 6
            
            couplingMatrices[2][1, 13+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[2][2, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[2][3, 15+12*repump] = -sqrt(2) / 3 * a - 1 / 6 * b
            couplingMatrices[2][4, 14+12*repump] = -sqrt(2) / 3
            couplingMatrices[2][5, 13+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[2][6, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[2][7, 15+12*repump] = a / 6 - sqrt(2) / 3 * b
            couplingMatrices[2][9, 13+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[2][10, 14+12*repump] = -1 / 3
            couplingMatrices[2][11, 15+12*repump] = -1 / (2 * sqrt(3))
            
            couplingMatrices[3][2, 13+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[3][3, 14+12*repump] = sqrt(2) / 3 * a + 1 / 6 * b
            couplingMatrices[3][3, 16+12*repump] = -sqrt(2) / 3 * a + b / 3
            couplingMatrices[3][4, 13+12*repump] = sqrt(2) / 3
            couplingMatrices[3][6, 13+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[3][7, 14+12*repump] = -a / 6 + sqrt(2) / 3 * b
            couplingMatrices[3][7, 16+12*repump] = -a / 3 - sqrt(2) / 3 * b
            couplingMatrices[3][10, 13+12*repump] = -1 / 6
            couplingMatrices[3][11, 14+12*repump] = -1 / (2 * sqrt(3))
            couplingMatrices[3][12, 15+12*repump] = -1 / sqrt(6)

            if repump == 1
                couplingMatrices[1][13:24, 25:28] = couplingMatrices[1][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)
                couplingMatrices[2][13:24, 25:28] = couplingMatrices[2][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)
                couplingMatrices[3][13:24, 25:28] = couplingMatrices[3][1:12, 25:28] .* sqrt(mol.v1BranchingRatioA)

                couplingMatrices[1][1:12, 25:28] = couplingMatrices[1][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
                couplingMatrices[2][1:12, 25:28] = couplingMatrices[2][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
                couplingMatrices[3][1:12, 25:28] = couplingMatrices[3][1:12, 25:28] .* sqrt(1-mol.v1BranchingRatioA)
            end

        else
            # excited states are all b states
            couplingMatrices[1][1, 14+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[1][1, 16+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[1][2, 15+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[1][4, 15+12*repump] = 1 / 3
            couplingMatrices[1][5, 14+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[1][5, 16+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[1][6, 15+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[1][8, 13+12*repump] = 1 / sqrt(3)
            couplingMatrices[1][9, 14+12*repump] = 1 / sqrt(6)
            couplingMatrices[1][10, 15+12*repump] = 1 / 3 / sqrt(2)
            
            couplingMatrices[2][1, 13+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[2][2, 16+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[2][3, 15+12*repump] = -a / 3 + b / 3 / sqrt(2)
            couplingMatrices[2][4, 14+12*repump] = -1 / 3
            couplingMatrices[2][5, 13+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[2][6, 16+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[2][7, 15+12*repump] = -a / 3 / sqrt(2) - b / 3
            couplingMatrices[2][9, 13+12*repump] = 1 / sqrt(6)
            couplingMatrices[2][10, 14+12*repump] = sqrt(2) / 3
            couplingMatrices[2][11, 15+12*repump] = 1 / sqrt(6)
            
            couplingMatrices[3][2, 13+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[3][3, 14+12*repump] = a / 3 - b / 3 / sqrt(2)
            couplingMatrices[3][3, 16+12*repump] = -a / 3 - sqrt(2) * b / 3
            couplingMatrices[3][4, 13+12*repump] = 1 / 3
            couplingMatrices[3][6, 13+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[3][7, 14+12*repump] = a / 3 / sqrt(2) + b / 3
            couplingMatrices[3][7, 16+12*repump] = sqrt(2) * a / 3 - b / 3
            couplingMatrices[3][10, 13+12*repump] = 1 / 3 / sqrt(2)
            couplingMatrices[3][11, 14+12*repump] = 1 / sqrt(6)
            couplingMatrices[3][12, 15+12*repump] = 1 / sqrt(3)

            if repump == 1
                # NOTE, there's really no reason this should ever execute...B and the vibrational repump are decoupled.  force this to not happen in main program.
                throw(ArgumentError("Repump and XToB are both 1. This is not allowed for now."))
                couplingMatrices[1][13: 24,25:28] = couplingMatrices[1][1:12, 25:28] .* sqrt(0)
                couplingMatrices[2][13: 24,25:28] = couplingMatrices[2][1:12, 25:28] .* sqrt(0)
                couplingMatrices[3][13: 24,25:28] = couplingMatrices[3][1:12, 25:28] .* sqrt(0)
            end
        end
    end


    function makeBCouplingMatrices!(bCouplingMatrices, XToB, repump, bichrom, mol::Molecule)
    # describes magnetic field induced larmor precession (for 'perpendicular' fields with-respect-to magnetic moment) and energy shifts (for parallel fields).  Depends on g factor for given hyperfine state
    
        gs = mol.gFactors

        bCouplingMatrices[1][2, 1] = gs[1]
        bCouplingMatrices[1][3, 2] = gs[1]
        bCouplingMatrices[1][6, 5] = gs[2]
        bCouplingMatrices[1][7, 6] = gs[2]
        bCouplingMatrices[1][9, 8] = sqrt(2) * gs[3]
        bCouplingMatrices[1][10, 9] = sqrt(3) * gs[3]
        bCouplingMatrices[1][11, 10] = sqrt(3) * gs[3]
        bCouplingMatrices[1][12, 11] = sqrt(2) * gs[3]

        if bichrom == 1
            bCouplingMatrices[1][14+12*repump, 13+12*repump] = gs[4]
            bCouplingMatrices[1][15+12*repump, 14+12*repump] = gs[4]
            bCouplingMatrices[1][14+12*repump+4, 13+12*repump+4] = gs[5]
            bCouplingMatrices[1][15+12*repump+4, 14+12*repump+4] = gs[5]
        elseif XToB == 1
            bCouplingMatrices[1][14+12*repump, 13+12*repump] = gs[5]
            bCouplingMatrices[1][15+12*repump, 14+12*repump] = gs[5]
        else
            bCouplingMatrices[1][14+12*repump, 13+12*repump] = gs[4]
            bCouplingMatrices[1][15+12*repump, 14+12*repump] = gs[4]
        end

        bCouplingMatrices[2][1, 1] = -gs[1]
        bCouplingMatrices[2][3, 3] = gs[1]
        bCouplingMatrices[2][5, 5] = -gs[2]
        bCouplingMatrices[2][7, 7] = gs[2]
        bCouplingMatrices[2][8, 8] = -2 * gs[3]
        bCouplingMatrices[2][9, 9] = -gs[3]
        bCouplingMatrices[2][11, 11] = gs[3]
        bCouplingMatrices[2][12, 12] = 2 * gs[3]

        if bichrom == 1
            bCouplingMatrices[2][13+12*repump, 13+12*repump] = -gs[4]
            bCouplingMatrices[2][15+12*repump, 15+12*repump] = gs[4]
            bCouplingMatrices[2][13+12*repump+4, 13+12*repump+4] = -gs[5]
            bCouplingMatrices[2][15+12*repump+4, 15+12*repump+4] = gs[5]
        elseif XToB == 1
            bCouplingMatrices[2][13+12*repump, 13+12*repump] = -gs[5]
            bCouplingMatrices[2][15+12*repump, 15+12*repump] = gs[5]
        else
            bCouplingMatrices[2][13+12*repump, 13+12*repump] = -gs[4]
            bCouplingMatrices[2][15+12*repump, 15+12*repump] = gs[4]
        end

        bCouplingMatrices[3][1, 2] = -gs[1]
        bCouplingMatrices[3][2, 3] = -gs[1]
        bCouplingMatrices[3][5, 6] = -gs[2]
        bCouplingMatrices[3][6, 7] = -gs[2]
        bCouplingMatrices[3][8, 9] = -sqrt(2) * gs[3]
        bCouplingMatrices[3][9, 10] = -sqrt(3) * gs[3]
        bCouplingMatrices[3][10, 11] = -sqrt(3) * gs[3]
        bCouplingMatrices[3][11, 12] = -sqrt(2) * gs[3]

        if bichrom == 1
            bCouplingMatrices[3][13+12*repump, 14+12*repump] = -gs[4]
            bCouplingMatrices[3][14+12*repump, 15+12*repump] = -gs[4]
            bCouplingMatrices[3][13+12*repump+4, 14+12*repump+4] = -gs[5]
            bCouplingMatrices[3][14+12*repump+4, 15+12*repump+4] = -gs[5]
        elseif XToB == 1
            bCouplingMatrices[3][13+12*repump, 14+12*repump] = -gs[5]
            bCouplingMatrices[3][14+12*repump, 15+12*repump] = -gs[5]
        else
            bCouplingMatrices[3][13+12*repump, 14+12*repump] = -gs[4]
            bCouplingMatrices[3][14+12*repump, 15+12*repump] = -gs[4]
        end
        
        if repump == 1
            bCouplingMatrices[1][13:24, 13:24] = bCouplingMatrices[1][1:12, 1:12]
            bCouplingMatrices[2][13:24, 13:24] = bCouplingMatrices[2][1:12, 1:12]
            bCouplingMatrices[3][13:24, 13:24] = bCouplingMatrices[3][1:12, 1:12]
        end
    end
end