module obeEvaluation

    using LinearAlgebra: mul!

    using ..moleculeVariables: Molecule
    using ..laserSettings: Lasers
    using ..generalSettings: GeneralSettings

    export densityMatrixChangeTerms!


    function makeFieldTerms!(fieldTerms::Vector{Vector{ComplexF64}}, r::Vector{Float64}, lasers::Lasers) 
        # These are different for 2D MOT
        # returns 'field terms' for all lasers. Field terms depend on the polarization type (and, if \sigma +/-, the sign)
        # This is basically the field for laser [i] due to the 6 (if 3D), or 1 (for Slower/push), or 4 (for all 2D lasers) passes of the beam expressed in the standard \sigma^- ([i][1]), \pi ([i][2]) and \sigma^+ ([i][3])
        # This is calculated in the way illustrated in JOSAB 6(11) 2023-2045 (1989) by Cohen-Tannoudji + Dalibard section 2.  See also Eq15-16 and subsequent expressions in my writeup for the 3D example.

        # IMPORTANT CAVEAT: all terms are 'pre-conjugated' since only the complex conjugate of this term is ever used (Eq 21 of my writeup).  Better to just express it pre-conjugated instead of 
        # repeatedly taking conjugates in the diff-eq solver

        polSign::Vector{Int64} = lasers.polSign
        polType::Vector{String} = lasers.polType
        wavenumberRatios::Vector{Float64} = lasers.wavenumberRatios
        waist::Float64 = lasers.beamWaist

        for i = 1:lasers.numLasers # iterate through all lasers
            if polType[i] == "Slower" 
                # polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
                fieldTerms[i][1] = 1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) + im * sin(r[3] * wavenumberRatios[i]))
                fieldTerms[i][2] = 0
                fieldTerms[i][3] = -1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) + im * sin(r[3] * wavenumberRatios[i]))

            elseif polType[i] == "Push"
                fieldTerms[i][1] = 1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) - im * sin(r[3] * wavenumberRatios[i]))
                fieldTerms[i][2] = 0
                fieldTerms[i][3] = -1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) - im * sin(r[3] * wavenumberRatios[i]))

            elseif polType[i] == "2DSS"
                fieldTerms[i][1] = polSign[i] * sin(r[1] * wavenumberRatios[i]) + im * cos(r[2] * wavenumberRatios[i])
                fieldTerms[i][2] = sqrt(2) * im * (cos(r[1] * wavenumberRatios[i]) - polSign[i]*sin(r[2] * wavenumberRatios[i]))
                fieldTerms[i][3] = polSign[i] * sin(r[1] * wavenumberRatios[i]) - im * cos(r[2] * wavenumberRatios[i])

            elseif polType[i] == "2DPerp"
                fieldTerms[i][1] = -sqrt(2) * im * cos(r[1] * wavenumberRatios[i])
                fieldTerms[i][2] = 2 * cos(r[2] * wavenumberRatios[i])
                fieldTerms[i][3] = -sqrt(2) * im * cos(r[1] * wavenumberRatios[i])

            elseif polType[i] == "2DPar"
                fieldTerms[i][1] = 0
                fieldTerms[i][2] = 2 * (cos(r[2] * wavenumberRatios[i])+cos(r[1] * wavenumberRatios[i]))
                fieldTerms[i][3] = 0

            elseif polType[i] == "3D"
                fieldTerms[i][1] = cos(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) + polSign[i] * sin(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) -
                im * (polSign[i] * sin(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) - cos(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2))
            
                fieldTerms[i][2] = sqrt(2) * im * (cos(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) +
                polSign[i] * sin(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2))
            
                fieldTerms[i][3] = cos(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) + polSign[i] * sin(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) +
                im * (polSign[i] * sin(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) - cos(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2))
            else
                throw(ArgumentError("Invalid polType value: $polType. It must be one of ['Slower', 'Push', '2DSS', '2DPerp', '2DPar', '3D']."))
            end
        end
    end

        
    function makeBFieldTerms!(bFieldTerms::Vector{ComplexF64}, r::Vector{Float64}, general::GeneralSettings)
        # expresses B field at position r in the \sigma^+/-, pi basis
        
        if general.bFieldSetting == "TwoD"
            bFieldTerms[1] = 1 / sqrt(2) * (r[1] - im * r[2])
            bFieldTerms[2] = -0 * (r[3])
            bFieldTerms[3] = 1 / sqrt(2) * (-r[1] - im * r[2])
        elseif general.bFieldSetting == "ThreeD"
            bFieldTerms[1] = 1 / sqrt(2) * (r[1] + im * r[2])
            bFieldTerms[2] = -1 * (r[3]) # note: really this should be -2r[3] for a quadropole field.  In practice, I prefer to run my f(r) for random direction at constant B.  So, assume \tilde{r}=(x,y,z/2).
            bFieldTerms[3] = 1 / sqrt(2) * (-r[1] + im * r[2])
        elseif general.bFieldSetting == "Static"
            bFieldTerms[1] = (im+1)/2
            bFieldTerms[2] = 0
            bFieldTerms[3] = (-1+im)/2
            # bFieldTerms[1] = 0
            # bFieldTerms[2] = 1
            # bFieldTerms[3] = 0
        else
            throw(ArgumentError("Invalid bFieldSetting value: $(general.bFieldSetting). It must be one of ['ThreeD', 'TwoD', 'Static']."))
        end
    end


    function densityMatrixChangeTerms!(du, u, p, t)
        # The meat of the program. Here's where the density matrix is actually evolved.

        # user inputs (these vary, things like initial position, velocity, laser params, etc.).  These are determined by the user-chosen parameters in the main program
        rInit = p[1]::Vector{Float64}
        vel = p[2]::Vector{Float64}
        stateEnergyMatrix = p[3]::Matrix{Float64}

        lasers = p[4]::Lasers
        general = p[5]::GeneralSettings
        mol = p[6]::Molecule

        # coupling matrices passed by user.  
        coupleMat1 = p[7]::Matrix{Float64}
        coupleMat2 = p[8]::Matrix{Float64}
        coupleMat3 = p[9]::Matrix{Float64}
        bCoupleMat1 = p[10]::Matrix{Float64}
        bCoupleMat2 = p[11]::Matrix{Float64}
        bCoupleMat3 = p[12]::Matrix{Float64}

        # coupling matrices used in decay calc
        coupleMatEff1 = p[13]::Matrix{ComplexF64}
        coupleMatEff2 = p[14]::Matrix{ComplexF64}
        coupleMatEff3 = p[15]::Matrix{ComplexF64}

        # decay 'masks' used in calculating the decay term.  
        decayMaskAllButTopLeft = p[16]::Matrix{Float64}
        decayMaskForCalcTopLeft = p[17]::Matrix{Int64}

        # pre-cached r Array
        r = p[18]::Vector{Float64}

        # pre-cached matrices for atom light term.  
        fieldTerms = p[19]::Vector{Vector{ComplexF64}}
        atomLightTerm = p[20]::Matrix{ComplexF64}
        atomLightTerm .= zeros(ComplexF64, size(coupleMat1,1), size(coupleMat1,2))
        
        # pre-cached matrices for b field term. 
        bFieldTerms = p[21]::Vector{ComplexF64}
        bFieldTermFull = p[22]::Matrix{ComplexF64}
        uProdBField = p[23]::Matrix{ComplexF64}
        bFieldProdU = p[24]::Matrix{ComplexF64}

        # pre-cached matrices for decay term
        decayFull = p[25]::Matrix{ComplexF64}
        pOnlyExcitedStates = p[26]::Matrix{ComplexF64}
        pTopLeft1PreMult = p[27]::Matrix{ComplexF64}
        pTopLeft2PreMult = p[28]::Matrix{ComplexF64}
        pTopLeft3PreMult = p[29]::Matrix{ComplexF64}
        pTopLeft1 = p[30]::Matrix{ComplexF64}
        pTopLeft2 = p[31]::Matrix{ComplexF64}
        pTopLeft3 = p[32]::Matrix{ComplexF64}

        # 1) evolve position
        @. r = rInit + vel * t # use .= to assign new values to each element of r (instead of creating a new array and allocating memory for it)

        # 2) Calculate field terms at new position
        makeFieldTerms!(fieldTerms, r, lasers)

        # 3)calculate -E dot D term (see Eq 21 of writeup)
        for i = 1:lasers.numLasers
            atomLightTerm .= atomLightTerm .+ sqrt(lasers.s0[i]/8) .* -exp(im * lasers.laserEnergy[i] * t + im * lasers.sidebandAmps[i] * sin(lasers.sidebandFreqs[i] * t)) .* 
            lasers.laserMasks[i] .* ((fieldTerms[i][1] .* coupleMat1) .+ (fieldTerms[i][2] .* coupleMat2) .+ (fieldTerms[i][3] .* coupleMat3))
        end
        @. atomLightTerm = atomLightTerm * exp(im * t * stateEnergyMatrix) # subtracts relevant hyperfine energies from 'laserEnergy'
        @. atomLightTerm = atomLightTerm + atomLightTerm' # needed here because, the way coupleMat is defined, 'atomLightTerm' up til now only has the top right half of the hermitian coupling matrix

        # 4) calculate -mu dot B term (see Eq 32 of writeup)
        makeBFieldTerms!(bFieldTerms, r, general)
        @. bFieldTermFull = general.bGrad * mol.normalizedBohrMag * (bFieldTerms[1] * bCoupleMat1 + bFieldTerms[2] * bCoupleMat2 + bFieldTerms[3] * bCoupleMat3) + atomLightTerm # 'bTermFull' also sums the -mu dot B term with the calculated -D dot E term
        
        # 5) take commutator of [H,u] where u is density matrix and H = -D dot E + -mu dot B
        mul!(uProdBField, u, bFieldTermFull) # uH
        mul!(bFieldProdU, bFieldTermFull, u) # Hu

        # 6) Take decay (aka 'coupling to reservoir') into account (Eq 46 of writeup)
        @. pOnlyExcitedStates = u * decayMaskForCalcTopLeft

        # 6A) these next 6 lines calculate the last term in eq 46 of writeup
        @. coupleMatEff1 = coupleMat1 * exp(-im * t * stateEnergyMatrix)
        mul!(pTopLeft1PreMult, coupleMatEff1, pOnlyExcitedStates)
        mul!(pTopLeft1, pTopLeft1PreMult, coupleMatEff1')

        @. coupleMatEff2 = coupleMat2 * exp(-im * t * stateEnergyMatrix)
        mul!(pTopLeft2PreMult, coupleMatEff2, pOnlyExcitedStates)
        mul!(pTopLeft2, pTopLeft2PreMult, coupleMatEff2')

        @. coupleMatEff3 = coupleMat3 * exp(-im * t * stateEnergyMatrix)
        mul!(pTopLeft3PreMult, coupleMatEff3, pOnlyExcitedStates)
        mul!(pTopLeft3, pTopLeft3PreMult, coupleMatEff3')

        @. decayFull = (u * decayMaskAllButTopLeft) + pTopLeft1 + pTopLeft2 + pTopLeft3 # u.*decayMask term represents 1st and 2nd term of eq 46 in writeup

        @. du = im * (uProdBField - bFieldProdU) + decayFull # finally, add the 'Liouville' term and the decay term (Eq 1 of writeup) to step the density matrix
    end
end