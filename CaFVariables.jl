const lamA=606e-9;#note: positions normalized to \tilde{x}=k_{SrF,X->A}x
const lamB=531e-9;
const lamRepump = 628e-9;
const v1BranchingRatioA = 1-.978;#ratio of population decay from A\pi,v=0 into X\Sigma,v=1
const v1BranchingRatioB = 1-.998;#ratio of population decay from B\Sigma,v=0 into X\Sigma,v=1
const gam = 2 * pi * 8.3e6;#linewidth (happens to be same for B and A for SrF).  Haven't figure out a good way to implement differing gamma in bichromatic traps.  For now, assume this is small effect?
const normalizedBohrMag = 1.758820e7 / 2 / gam;
const mass = (40 + 19) * 1.67e-27;#mass of CaF
const a = 0.772496;#j mixing terms a and b, from Eunmi Chae thesis
const b = sqrt(1 - a^2);
const gs = [-0.295, 0.795, 0.5, -0.02,1.02];#g values.  First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
const stateEnergiesGround = [0.0, 9.3, 15.1, 18.1]; #X\Sigma hyperfine energies, with 0 corresponding to the F=1\DOWN energy
const stateEnergiesExcited = [0.0,-3.1];#Energy of F=0 relative to "0" (F=1).  Entry 1 for A state, Entry 2 for B State.  Splitting negligible in A state (probably not zero, update if we ever measure this)

const kA = 2 * pi / lamA; #wavenumber
const kB = 2 * pi / lamB;
const kRepump = 2 * pi / lamRepump;
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2