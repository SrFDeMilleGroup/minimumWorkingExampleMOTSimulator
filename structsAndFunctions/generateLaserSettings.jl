module laserSettings

    using ..structs: Molecule, Lasers # each successive . leads to the parent of the current module

    export generateLaserSettings


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