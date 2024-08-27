module laserSettings

    using ..moleculeVariables: Molecule # each successive . leads to the parent of the current module

    # include("../simulationSettings/laserSettings.jl")
    # using .laserSettings: s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps, beamWaistInMM

    export Lasers, generateLaserSettings

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

    wavenumberRatios: ratio of k_{Laser} to k_{A} 

    laserMasks: make "Masks" for lasers based on what transition the laser corresponds to. This is multiplied element-wise with coupling matrix in the OBE solver (densityMatrixChangeTerms). 
        This is zero for terms that are not coupled together by the matrix (e.g., turns off X->A coupling for X->B laser, etc. and 1 for terms that are)

    numLasers: number of lasers
    
    beamWaistInMM: in unit of mm, only used if polType is 3D. Handles finite MOT beam waists

    laserBeamWaist: in unit of 1/k, converted from beamWaistInMM
    """

    @kwdef struct Lasers
        s0::Vector{Float64}
        laserEnergy::Vector{Float64}
        polSign::Vector{Int64}
        whichTransition::Vector{String}
        polType::Vector{String}
        sidebandFreqs::Vector{Float64}
        sidebandAmps::Vector{Float64}
        wavenumberRatios::Vector{Float64}
        laserMasks::Vector{Matrix{Float64}}
        numLasers::Int64
        beamWaistInMM::Float64
        beamWaist::Float64
    end

    function generateLaserSettings(mol::Molecule)::Lasers

        include("./simulationSettings/laserSettings.jl")

        @assert length(s0) == length(laserEnergy) == length(polSign) == length(whichTransition) == length(polType) == length(sidebandFreqs) == length(sidebandAmps) "All laser settings arrays must have the same length."
        @assert all([i in [1, -1] for i in polSign]) "Invalid polSign value(s): $polSign. All values must be either 1 or -1."
        @assert all([i in ["XA", "XB", "XARepump"] for i in whichTransition]) "Invalid whichTransition value(s): $whichTransition. All values must be either 'XA', 'XB', or 'XARepump'."
        @assert all([i in ["3D", "2DSS", "2DPar", "2DPerp", "Slower", "Push"] for i in polType]) "Invalid polType value(s): $polType. All values must be one of ['3D', '2DSS', '2DPar', '2DPerp', 'Slower', 'Push']."
        @assert beamWaistInMM > 0 "Invalid waistInMM value: $waistInMM. It must be greater than 0."

        bichrom = (("XA" in whichTransition) && ("XB" in whichTransition)) ? 1 : 0 # winds up 0 if only XA of XB are used, 1 if both are
        repump = ("XARepump" in whichTransition) ? 1 : 0 # winds up 0 if no repump, 1 if there are repumps

        numZeemanStatesGround = 12 + 12 * repump
        numZeemanStatesExcited = 4 + 4 * bichrom # 4 for XA or XB along, 8 for both
        numZeemanStatesTotal = numZeemanStatesGround + numZeemanStatesExcited

        wavenumberRatios = Vector{Float64}(undef, length(whichTransition))
        laserMasks = [zeros(numZeemanStatesTotal, numZeemanStatesTotal) for i=1:length(whichTransition)]
        for (i, currTransition) in enumerate(whichTransition)
            if currTransition == "XA"
                wavenumberRatios[i] = 1.0
                laserMasks[i][1:12, (13+12*repump):(16+12*repump)] .= 1
            elseif currTransition == "XB"
                wavenumberRatios[i] = mol.kB / mol.kA
                laserMasks[i][1:12, (13+12*repump+4*bichrom):(16+12*repump+4*bichrom)] .= 1
            elseif currTransition == "XARepump"
                wavenumberRatios[i] = mol.kRepump / mol.kA
                laserMasks[i][13:24, 25:28] .= 1
            else
                error("Invalid whichTransition value: $currTransition. Must be one of 'XA', 'XB', or 'XARepump'.")
            end
        end

        return Lasers(s0, laserEnergy, polSign, whichTransition, polType, sidebandFreqs, sidebandAmps, wavenumberRatios, laserMasks, length(s0), beamWaistInMM, beamWaistInMM*1e-3*mol.kA)

    end
end