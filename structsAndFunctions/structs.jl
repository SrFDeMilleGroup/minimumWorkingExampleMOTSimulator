module structs

    """
    - unit of frequency: Gamma, where Gamma is the linewidth of a molecular transition, usually X-A transition.
    - unit of energy: hbar * Gamma
    - unit of time: 1 / Gamma.
    - unit of length: 1 / k (= wavelength / 2pi), where k is the (angular) wavenumber of a molecular transition, usually X-A transition.
    - unit of velocity: Gamma / k.
    - unit of force: 1e-3 * hbar * Gamma * k.
    """

    export Molecule, Lasers, GeneralSettings

    @kwdef struct Molecule

        """
        name: name of molecule

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
        accelFactor: normalized force units in program are 1e-3\\hbar*k*\\gam. So the factor converts this to mm/ms^2
        """

        name::String

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
        hbar::Float64 = 1.05457182e-34
        accelFactor::Float64 = (1e-3 * hbar * kA * Gamma / mass) * 1e-3
    end


    @kwdef struct Lasers

        """
        s0: Saturation parameter corresponding the peak laser intensity of a single laser pass.  I_Sat ~ 3 mW/cm^2 for XA and ~ 4 mW/cm^2 for XB for SrF.

        laserEnergy: in unit of Gamma, laser detuning measured from the lowest ground hyperfine level.

        polType: lsaer polarization and propagation direction. Options are:
            '3D': 6 laser beams propagate in +/-x, +/-y, +/-z 6 directions respectively. 
                    All lasers are circularly polarized, but x and y lasers have opposite polarization handedness from z lasers, as required by 3D MOT quadrupole B field. 
                    Sigma +/- polarization is indicated by polSign. This is mainly used to simulate 3D MOT or molasses cooling.
            '2DSS': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. 
                    All lasers are circularly polarized, but x lasers have opposite handedness from y lasers, as required by 2D MOT B field. 
                    Sigma +/- polarization is indicated by polSign. This is mainly used to simulate 2D MOT.
            '2DPar': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. 
                    All lasers are linear polarized along z axis. This is mainly used to simulate transverse cooling of a molecule beam, assuming it travels along z axis.
            '2DPerp': Same as '2DPar', except x lasers are polarized along y, and y lasers are polarized along z axis.
            'Slower': 1 laser beam propagate along -z direction, and is linearly polarized along x axis.
            'Push': Same as 'Slower' except the laser propagates along +z direction.
                
        polSign: -/+ determine sigma-/+ for polType is 3D or 2DSS. For other 'polType' these are unused.

        sidebandFreqs: in unit of Gamma, modulation frequencies of EO modulators.

        sidebandAmps: in unit of radians, modulation depths of EO modulators.

        whichTransition: can be "XA" (couples X,v=0 to A,v=0), "XB" (couples X,v=0 to B,v=0), and "XARepump" (couples X,v=1 to A,v=0).
            if there is no "XARepump", vibrational branching IS TURNED OFF (obviously, or else all population would accumulate in v=1).
            "CouplingMatrices" and "laser masks" populate based on what values of "whichTransition" are chosen. 

        wavenumberRatios: ratio of k_{Laser} to k_{A}.

        laserMasks: make "Masks" for lasers based on what transition the laser corresponds to. This is multiplied element-wise with coupling matrix in the OBE solver (densityMatrixChangeTerms). 
            This is zero for terms that are not coupled together by the matrix (e.g., turns off X->A coupling for X->B laser, etc. and 1 for terms that are).

        numLasers: number of lasers.
        
        beamWaistInMM: mm, laser beam waist. Used to simulate the effects of the finite size of 3D MOT laser beams. 
            This is only used if polType is '3D', otherwise the laser beams are assumed infinitely large.

        laserBeamWaist: in unit of 1/k, converted from beamWaistInMM.
        """

        s0::Vector{Float64}
        laserEnergy::Vector{Float64}
        polType::Vector{String}
        polSign::Vector{Int64}
        sidebandFreqs::Vector{Float64}
        sidebandAmps::Vector{Float64}
        whichTransition::Vector{String}
        wavenumberRatios::Vector{Float64}
        laserMasks::Vector{Matrix{Float64}}
        numLasers::Int64
        beamWaistInMM::Float64
        beamWaist::Float64
    end


    @kwdef struct GeneralSettings

        """
        simulationType: string, notes for yourself about simulation types, e.g., redMOT, blueMOT, transCooling, etc.

        numTrialsPerValueSet: number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds).

        forceProfile: Options are:
            'ThreeD': Used for 3D MOT and molasses simulations. 
                        In this case, 'displacementsInMM' and 'userSpeeds' are treated as magnitudes of 3D displacements and velocities respectively. 'longSpeeds' is ignored. 
                        The direction of initial displacements is indicated by 'initDispDir'. And the direction of velocity is set with respect to displacements by 'velDirRelToR'. 
                        In order to account for laser field variation, the actual initial displacements in computation are randomly sampled from a cube of size of one wavelength round the values specified here. 
                        The projections of 3D force on initial displacement f \\dot r / |r| and f \\dot v / |v| are calculated and returned. 
            'TwoD': Used for slowing and transverse cooling simulations. 
                        In this case, 'displacementsInMM' and 'userSpeeds' are treated as magnitudes of 2D displacements and velocities in x-y plane. 
                        Displacement along z is taken as zero (up to a random value within +/- 1/2 wavelength from 0). 'longSpeeds' is used as velocity along z direction. 
                        'initDispDir' is ignored, and the direction of initial displacements in x-y plane is always random. The relative direction of velocity and displacement in x-y plane is set by 'velDirRelToR'. 
                        The projection of 3D force on x-y plane displacement f \\dot r / |r| and f \\dot v / |v|, as well as its z component f_z are calculated and returned.

        displacementsInMM: mm, magnitude of initial displacements. Simulation will iterate through the entire list.

        initDispDir: Direction of initial displacememnt. Only used if forceProfile is 'ThreeD'. Valid options are 'XY' ((x+y)/sqrt(2) direction, where slowed molecules come into MOT region for most experiments), 'Z' and 'Random'.

        longSpeeds: longitudinal speeds in units of Gamma/k (=4.4m/s for SrF), Only used if 'foreProfile' is 'TwoD'.

        userSpeeds: in unit of Gamma/k (=4.4m/s for SrF), speeds in xy plane (for 2d force profile) or in 3D.

        velDirRelToR: relative direction of molecule velocity, w.r.t. initial displacement r. Options are ["Same", "Orthogonal", "Opposite", "Random"].

        bFieldSetting: can set to 3D quadrupole "ThreeD" (e.g. 3D-MOT"), 2D quadrupole "TwoD" (e.g. 2D-MOT"),
        or static "StaticXY" (in (x+y)/sqrt(2) direction) or "StaticZ" (in z direction). 
        (Static B fields are used for 2D transverse slowing primarily, could also use to simulate e.g. lambda-cooling in 3D field).

        bGradReal: Gauss/cm for B field gradient, or Gauss for uniform B field, depending on 'bFieldSetting'.

        bGrad: convert bGradReal to units Gauss * wavevector, unless bFieldSetting is static, then this remains unchanged in Gauss.
        """

        simulationType::String
        numTrialsPerValueSet::Int64
        forceProfile::String
        displacementsInMM::Vector{Float64}
        initDispDir::String
        longSpeeds::Vector{Float64}
        userSpeeds::Vector{Float64}
        velDirRelToR::String
        bFieldSetting::String
        bGradReal::Float64
        bGrad::Float64
        
    end

end