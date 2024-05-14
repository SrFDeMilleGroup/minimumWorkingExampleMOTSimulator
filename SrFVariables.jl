const lamA=663e-9;#note: positions normalized to \tilde{x}=k_{SrF,X->A}x
const lamB=579e-9;
const lamRepump = 685e-9;
const v1BranchingRatioA = 1/50;#ratio of population decay from A\pi,v=0 into X\Sigma,v=1
const v1BranchingRatioB = 3.866e-3;#ratio of population decay from B\Sigma,v=0 into X\Sigma,v=1
const gam = 2 * pi * 6.63e6;#linewidth (happens to be same for B and A).  Haven't figure out a good way to implement differing gamma in bichromatic traps.  
const normalizedBohrMag = 0.2114; #\mu_{B}/\hbar\Gamma in units 1/Gauss
const mass = (88 + 19) * 1.67e-27;#mass of SrF
const a = 0.888;#j mixing terms a and b, see john barry thesis chapt 2
const b = sqrt(1 - a^2);
const gs = [-0.47, 0.97, 0.5, -0.088,1.088];#g values.  First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
const stateEnergiesGround = [0.0, 7.5, 19.6, 25.9]; #X\Sigma hyperfine energies, with 0 corresponding to the F=1\DOWN energy
const stateEnergiesExcited = [0.0,-2.0];#Energy of F=0 relative to "0" (F=1).  Entry 1 for A state, Entry 2 for B State.  Splitting negligible in A state (probably not zero, update if we ever measure this)

const kA = 2 * pi / lamA; #wavenumber
const kB = 2 * pi / lamB;
const kRepump = 2 * pi / lamRepump;
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2