module generalSettings

    using ..structs: Molecule, GeneralSettings

    export generateGeneralSettings


    function generateGeneralSettings(mol::Molecule)::GeneralSettings

        include("./simulationSettings/generalSettings.jl")

        @assert numTrialsPerValueSet > 0 "Invalid numTrialsPerValueSet value: $numTrialsPerValueSet. It must be greater than 0."
        @assert initDispDir in ["XY", "Z"] "Invalid initDispDir value: $initDispDir. It must be either 'XY' or 'Z'. "
        @assert velDirRelToR in ["Same", "Orthogonal", "Opposite", "Random"] "Invalid velDirRelToR value: $velDirRelToR. It must be one of ['Same', 'Orthogonal', 'Opposite', 'Random']."
        @assert forceProfile in ["ThreeD", "TwoD"] "Invalid forceProfile value: $forceProfile. It must be either 'ThreeD' or 'TwoD'."
        @assert bFieldSetting in ["ThreeD", "TwoD", "Static"] "Invalid bFieldSetting value: $bFieldSetting. It must be one of ['ThreeD', 'TwoD', 'Static']."

        bGrad = bFieldSetting == "Static" ? bGradReal : (1 / mol.kA * 1e2) * bGradReal
        return GeneralSettings(numTrialsPerValueSet, displacementsInMM, initDispDir, longSpeeds, userSpeeds, velDirRelToR, forceProfile, bGradReal, bGrad, bFieldSetting)
    
    end
end