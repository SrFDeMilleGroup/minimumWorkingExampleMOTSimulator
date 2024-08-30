# unit of energy: hbar * Gamma
# unit of velocity: Gamma / k, where k is the wavevector
# unit of time: 1 / Gamma
# unit of length: 1 / k (= wavelength / 2pi)
# unit of force: 1e-3 * hbar * Gamma * k (??)

module structs

    export Molecule, Lasers, GeneralSettings

    @kwdef struct Molecule

        """
        lambdaA: m, XA trasition wavelength. Note: positions normalized to \\tilde{x}=k_{X->A}x
        lambdaB: m, XB transition wavelength
        lambdaRepump: m, v10 XA repump transition wavelength
        v1BranchingRatioA: ratio of population decay from A\\pi,v=0 into X\\Sigma,v=1
        v1BranchingRatioB: ratio of population decay from B\\Sigma,v=0 into X\\Sigma,v=1
        Gamma: Hz, linewidth (happens to be same for B and A). Haven't figure out a good way to implement differing gamma in bichromatic traps.

        mass: kg, mass of molecule
        jMixingRatioA: j mixing terms a and b, see john barry thesis chapt 2
        gFactors: g values. First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
        stateEnergiesGround: in unit of \\hbar\\Gamma, X\\Sigma hyperfine energies, with 0 corresponding to the F=1\\DOWN energy
        stateEnergiesExcited: in unit of \\hbar\\Gamma, Energy of F=0 relative to "0" (F=1). Entry 1 for A state, Entry 2 for B State. 
            (Splitting negligible in A state (probably not zero, update if we ever measure this)
        
        jMixingRatioB: j mixing terms a and b, see john barry thesis chapt 2
        normalizedBohrMag: in units \\hbar * \\Gamma / Gauss, bohr magneton = 1.4 MHz/G
        kA: wavevector
        kB: wavevector
        kRepump: wavevector
        velFactor: unit of velocity in this simulation, Gamma / kA
        hbar: SI units
        accelFactor: normalized force units in program are 1e-3\\hbar*k*\\gam. So the factor converts this to m/s^2
        """

        lamdaA::Float64
        lamdaB::Float64
        lamdaRepump::Float64
        v1BranchingRatioA::Float64
        v1BranchingRatioB::Float64
        Gamma::Float64
        
        mass::Float64
        jMixingRatioA::Float64
        gFactors::Vector{Float64}
        stateEnergiesGround::Vector{Float64}
        stateEnergiesExcited::Vector{Float64}

        jMixingRatioB::Float64 = sqrt(1 - jMixingRatioA^2)
        normalizedBohrMag::Float64 = 2 * pi * 1.39962449171e6 / Gamma
        kA::Float64 = 2 * pi / lamdaA
        kB::Float64 = 2 * pi / lamdaB
        kRepump::Float64 = 2 * pi / lamdaRepump
        velFactor::Float64 = Gamma / kA
        hbar::Float64 = 1.05e-34
        accelFactor::Float64 = 1e-3 * hbar * kA * Gamma / mass
    end


    @kwdef struct Lasers

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


    @kwdef struct GeneralSettings

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

end