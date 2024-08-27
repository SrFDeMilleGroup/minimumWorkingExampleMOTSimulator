module generalSettings

    using ..moleculeVariables: Molecule

    export generateGeneralSettings, GeneralSettings

    """
    numTrialsPerValueSet: number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds)

    displacementsInMM: mm, initial displacements from either origin (if 3D) or else z-axis (if 2D) in mm

    initDispDir: if "XY", will force initial r to go along (x+y)/sqrt(2). 
        Simulates slowing/trapping of molecules moving along slowing axis in tandem with velDirToR = "Same". If "Z", force go along Z

    longSpeeds: longitudinal speeds in units of Gamma/k (=4.4m/s for SrF), doesn't matter for 3D sims, sets vel to 140 m/s

    userSpeeds: in unit of Gamma/k (=4.4m/s for SrF), speeds in xy plane (for 2d force profile) or in 3D

    velDirRelToR: relative direction of molecule velocity, w.r.t. initial displacement r. Options are ["Same", "Orthogonal", "Opposite", "Random"]

    forceProfile: either "ThreeD", (forces calculated are (f \\dot r) / |r|, (f \\dot v) / |v|), 
        or "TwoD" (f \\dot (rx,ry,0) / |(rx,ry,0)|, f \\dot (vx,vy,0) / |(vx,vy,0|, and fz are all calculated)

    bGradReal: in units Gauss/cm, unless bFieldSetting == "Static", then this becomes the static field in Gauss

    bGrad: convert bGradReal to units Gauss * wavevector, unless bFieldSetting == "Static", then this remains unchanged in Gauss

    bFieldSetting: can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT"),
        or static "Static" (2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).
        
    """

    @kwdef struct GeneralSettings
        numTrialsPerValueSet::Int64
        displacementsInMM::Vector{Float64}
        initDispDir::String
        longSpeeds::Vector{Float64}
        userSpeeds::Vector{Float64}
        velDirRelToR::String
        forceProfile::String
        bGradReal::Float64
        bGrad::Float64
        bFieldSetting::String
    end

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