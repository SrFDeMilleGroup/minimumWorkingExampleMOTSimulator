module forceCalculation

    using LinearAlgebra: tr, mul!
    using Random: Xoshiro

    using ..moleculeVariables: Molecule
    using ..laserSettings: Lasers
    using ..generalSettings: GeneralSettings

    export makeForceVsTime!


    # everything here is used to calculate a force given a density matrix, see section 1.6 of writeup
    function makeDFieldTerms!(dFieldTerms, r::Vector{Float64}, lasers::Lasers)
        # dfieldterms (dE/dr) will be 3x3 matrix, first element is xyz second is sig+ pi sig-
        # fieldTerms = zeros(ComplexF64,3,1);
        # pre conjugated, just like 'makeFieldTerms'

        polSign::Vector{Int64} = lasers.polSign
        polType::Vector{String} = lasers.polType
        wavenumberRatios::Vector{Float64} = lasers.wavenumberRatios
        waist::Float64 = lasers.beamWaist

        for i = 1:lasers.numLasers 
            # iterates through all lasers
            # polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
            if polType[i] == "Slower" 
                dFieldTerms[i][1, 1] = 0
                dFieldTerms[i][1, 2] = 0
                dFieldTerms[i][1, 3] = 0

                dFieldTerms[i][2, 1] = 0
                dFieldTerms[i][2, 2] = 0
                dFieldTerms[i][2, 3] = 0

                dFieldTerms[i][3, 1] = 1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) + im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i]
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = -1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) + im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i]

            elseif polType[i] == "Push"
                dFieldTerms[i][1, 1] = 0
                dFieldTerms[i][1, 2] = 0
                dFieldTerms[i][1, 3] = 0

                dFieldTerms[i][2, 1] = 0
                dFieldTerms[i][2, 2] = 0
                dFieldTerms[i][2, 3] = 0

                dFieldTerms[i][3, 1] = 1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) - im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i]
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = -1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) - im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i]

            elseif polType[i] == "2DSS"
                dFieldTerms[i][1, 1] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][1, 2] = -sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][1, 3] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]

                dFieldTerms[i][2, 1] = -im * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][2, 2] = -sqrt(2) * im * (polSign[i] * cos(r[2] * wavenumberRatios[i])) * wavenumberRatios[i]
                dFieldTerms[i][2, 3] = im * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]

                dFieldTerms[i][3, 1] = 0
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = 0

            elseif polType[i] == "2DPerp"
                dFieldTerms[i][1, 1] = sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][1, 2] = 0
                dFieldTerms[i][1, 3] = sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]

                dFieldTerms[i][2, 1] = 0
                dFieldTerms[i][2, 2] = -2 * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][2, 3] = 0

                dFieldTerms[i][3, 1] = 0
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = 0

            elseif polType[i] == "2DPar"
                dFieldTerms[i][1, 1] = 0
                dFieldTerms[i][1, 2] = -2 * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][1, 3]=0

                dFieldTerms[i][2, 1] = 0
                dFieldTerms[i][2, 2] = -2 * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]
                dFieldTerms[i][2, 3] = 0

                dFieldTerms[i][3, 1] = 0
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = 0

            elseif polType[i] == "3D"
                dFieldTerms[i][1, 1] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i]
                dFieldTerms[i][1, 2] = -sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i]
                dFieldTerms[i][1, 3] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i]

                dFieldTerms[i][2, 1] = -im * sin(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2) * wavenumberRatios[i]
                dFieldTerms[i][2, 2] = sqrt(2) * im * (polSign[i] * cos(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2)) * wavenumberRatios[i]
                dFieldTerms[i][2, 3] = im * (sin(r[2] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[3]^2)/waist^2)) * wavenumberRatios[i]
        
                dFieldTerms[i][3, 1] = (-sin(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) - im * (polSign[i] * cos(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2))) * wavenumberRatios[i]
                dFieldTerms[i][3, 2] = 0
                dFieldTerms[i][3, 3] = (-sin(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2) + im * (polSign[i] * cos(r[3] * wavenumberRatios[i]) * exp(-2*(r[1]^2+r[2]^2)/waist^2))) * wavenumberRatios[i]
            else
                throw(ArgumentError("Invalid polType value: $polType. It must be one of ['Slower', 'Push', '2DSS', '2DPerp', '2DPar', '3D']."))
            end
        end
    end


    function forceCalc!(force, dFieldTerms::Vector{Matrix{ComplexF64}}, rho::Matrix{ComplexF64}, 
        lasers::Lasers, couplingMatrices::Vector{Matrix}, stateEnergyMatrix::Matrix{Float64}, t::Float64)
        # calculates force given position and lasers (used to calculate dFieldTerms) and density matrix \rho.  
        # Both r(t) and \rho(t) are recorded vs time by the OBE solver, so this runs afterwards to calculate what forces the particle experienced over the trajectory

        # force pre-factor is calculated for each laser [i]. Has the rotating-frame frequency exponent + phase modulation term + intensity term \sqrt(s0/8). Hyperfine energies are subtracted later
        forcePrefactor = zeros(ComplexF64, 1, lasers.numLasers)
        for i = 1:lasers.numLasers
            forcePrefactor[i] = sqrt(lasers.s0[i] / 8) * exp(im * lasers.laserEnergy[i] * t + im * lasers.sidebandAmps[i] * sin(lasers.sidebandFreqs[i] * t))
        end

        # calculate x force. Implements Eq 48 of main writeup
        dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho, 1), size(rho, 2))
        dRhoDPosTimesDensityMatContainer = zeros(ComplexF64, size(rho, 1), size(rho, 2))
        for i = 1:lasers.numLasers
            @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + forcePrefactor[i] * (dFieldTerms[i][1, 1] * couplingMatrices[1] + dFieldTerms[i][1, 2] * couplingMatrices[2] + dFieldTerms[i][1, 3] * couplingMatrices[3]) * lasers.laserMasks[i]
        end
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix * exp(im * t * stateEnergyMatrix)
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + dRhoDPosCalcMatrix'
        mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix) # multiplies dp_{x}/dt by density matrix \rho
        force[1] = real(tr(dRhoDPosTimesDensityMatContainer)) # takes trace of \rho*dp_{x}/dt to determine average force over enemble (Eq 49 of writeup)

        # similarly, calculate y and z force
        dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho, 1), size(rho, 2))
        for i = 1:lasers.numLasers
            @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + forcePrefactor[i] * (dFieldTerms[i][2, 1] * couplingMatrices[1] + dFieldTerms[i][2, 2] * couplingMatrices[2] + dFieldTerms[i][2, 3] * couplingMatrices[3]) * lasers.laserMasks[i]
        end
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix * exp(im * t * stateEnergyMatrix)
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + dRhoDPosCalcMatrix'
        mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix)
        force[2] = real(tr(dRhoDPosTimesDensityMatContainer))

        dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho, 1), size(rho, 2))
        for i = 1:lasers.numLasers
            @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + forcePrefactor[i] * (dFieldTerms[i][3, 1] * couplingMatrices[1] + dFieldTerms[i][3, 2] * couplingMatrices[2] + dFieldTerms[i][3, 3] * couplingMatrices[3]) * lasers.laserMasks[i]
        end
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix * exp(im * t * stateEnergyMatrix)
        @. dRhoDPosCalcMatrix = dRhoDPosCalcMatrix + dRhoDPosCalcMatrix'
        mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix)
        force[3] = real(tr(dRhoDPosTimesDensityMatContainer))
    end


    function makeForceVsTime!(forceVsTime, times::Vector{Float64}, rhos::Vector{Matrix{ComplexF64}}, lasers::Lasers, couplingMatrices::Vector{Matrix}, 
        stateEnergyMatrix::Matrix{Float64}, rInit::Vector{Float64}, v::Vector{Float64})
        # given a set of times, an initial position and velocity, the lasers used, and \rho(t), calculate force vs t

        # initialize some stuff
        dFieldContainer = Vector{Matrix{ComplexF64}}(undef, lasers.numLasers)
        for i = 1:lasers.numLasers
            dFieldContainer[i] = zeros(ComplexF64, 3, 3)
        end
        forceCalcContainer = zeros(ComplexF64, 3, 1)
        r = Vector{Float64}(undef, 3)

        # iterate through time, propegating r in the same way done in the OBEs.  Then determine force experienced given r(t), \rho(t), and the lasers used
        for i = 1:length(times)
            @. r = rInit + v * times[i]
            makeDFieldTerms!(dFieldContainer, r, lasers)
            forceCalc!(forceCalcContainer, dFieldContainer, rhos[i], lasers, couplingMatrices, stateEnergyMatrix, times[i])
            forceVsTime[i, :] = forceCalcContainer
        end
    end

end