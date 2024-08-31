module moleculeVariables

    """
    unit of energy: hbar * Gamma
    unit of velocity: Gamma / k, where k is the wavevector
    unit of time: 1 / Gamma
    unit of length: 1 / k (= wavelength / 2pi)
    unit of force: 1e-3 * hbar * Gamma * k (??)
    """

    using ..structs: Molecule

    export SrF, CaF, BaF, MgF, CaOH, SrOH


    SrF = Molecule(name = "SrF",
                   lamdaA = 663e-9, # m, note: positions normalized to \tilde{x}=k_{SrF,X->A}x 
                   lamdaB = 579e-9, # m
                   lamdaRepump = 685e-9, # m
                   v1BranchingRatioA = 1/50, # ratio of population decay from A\pi,v=0 into X\Sigma,v=1
                   v1BranchingRatioB = 3.866e-3, # ratio of population decay from B\Sigma,v=0 into X\Sigma,v=1
                   Gamma = 2 * pi * 6.63e6, # Hz, linewidth (happens to be same for B and A).  Haven't figure out a good way to implement differing gamma in bichromatic traps.
                   mass = (88 + 19) * 1.67e-27, # kg, mass of SrF
                   jMixingRatioA = 0.888, # j mixing terms a and b, see john barry thesis chapt 2
                   gFactors = [-0.47, 0.97, 0.5, -0.088, 1.088], # g values.  First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
                   stateEnergiesGround = [0.0, 7.5, 19.6, 25.9], # in unit of \hbar\Gamma, X\Sigma hyperfine energies, with 0 corresponding to the F=1\DOWN energy
                   stateEnergiesExcited = [0.0, -2.0] # in unit of \hbar\Gamma, Energy of F=0 relative to "0" (F=1).  Entry 1 for A state, Entry 2 for B State.  Splitting negligible in A state (probably not zero, update if we ever measure this
                   )

    CaF = Molecule(name = "CaF",
                   lamdaA = 606e-9,
                   lamdaB = 531e-9,
                   lamdaRepump = 628e-9,
                   v1BranchingRatioA = 1 - 0.978,
                   v1BranchingRatioB = 1 - 0.998,
                   Gamma = 2 * pi * 8.3e6,
                   mass = (40 + 19) * 1.67e-27,
                   jMixingRatioA = 0.772496,
                   gFactors = [-0.295, 0.795, 0.5, -0.02, 1.02],
                   stateEnergiesGround = [0.0, 9.3, 15.1, 18.1],
                   stateEnergiesExcited = [0.0, -3.1]
                   )

    BaF = Molecule(name = "BaF",
                   lamdaA = 860e-9,
                   lamdaB = 736.7e-9,
                   lamdaRepump = 896e-9,
                   v1BranchingRatioA = 1 - 0.9508,
                   v1BranchingRatioB = 1 - 0.81,
                   Gamma = 2 * pi * 3.0e6,
                   mass = (138 + 19) * 1.67e-27,
                   jMixingRatioA = 0.9593,
                   gFactors = [0.015, 0.485, 0.5, -0.202, 1.202],
                   stateEnergiesGround = [0.0, 9.3, 39.2, 50.6],
                   stateEnergiesExcited = [-5.6, -3.1]
                   )

    MgF = Molecule(name = "MgF",
                   lamdaA = 359.3e-9,
                   lamdaB = 268.9e-9,
                   lamdaRepump = 368.7e-9,
                   v1BranchingRatioA = 1 - 0.97,
                   v1BranchingRatioB = 1 - 0.998,
                   Gamma = 2 * pi * 22e6,
                   mass = (24 + 19) * 1.67e-27,
                   jMixingRatioA = 0.6990,
                   gFactors = [-0.21, 0.71, 0.5, -0.0002, 1.0002],
                   stateEnergiesGround = [0.0, 5.0, 10.5, 10.9],
                   stateEnergiesExcited = [0.0, -3.1]
                   )

    CaOH = Molecule(name = "CaOH",
                    lamdaA = 626.4e-9,
                    lamdaB = 555.2e-9,
                    lamdaRepump = 651e-9,
                    v1BranchingRatioA = 1 - 0.9521,
                    v1BranchingRatioB = 1 - 0.9711,
                    Gamma = 2 * pi * 6.4e6,
                    mass = (40 + 16 + 1) * 1.67e-27,
                    jMixingRatioA = 0.999633,
                    gFactors = [-0.33, 0.83, 0.5, -0.04, 1.04],
                    stateEnergiesGround = [0.0, 0.0, 8.0, 8.2],
                    stateEnergiesExcited = [0.0, -3.5]
                    )

    SrOH = Molecule(name = "SrOH",
                    lamdaA = 688e-9,
                    lamdaB = 611e-9,
                    lamdaRepump = 713.4e-9,
                    v1BranchingRatioA = 1 - 0.96,
                    v1BranchingRatioB = 1 - 0.98,
                    Gamma = 2 * pi * 7e6,
                    mass = (88 + 16 + 1) * 1.67e-27,
                    jMixingRatioA = 0.999633,
                    gFactors = [-0.33, 0.83, 0.5, -0.096, 1.096],
                    stateEnergiesGround = [0.0, 0.0, 15.6, 15.7],
                    stateEnergiesExcited = [0.0, 0.0]
                    )

end
