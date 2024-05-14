const lamA=626.4e-9;#note: positions normalized to \tilde{x}=k_{SrF,X->A}x
const lamB=555.2-9;#note: this has 18 MHz linewidth.  Will try to add in different linewidths later, if possible.  
const lamRepump = 651e-9;
const v1BranchingRatioA = 1-.9521;#ratio of population decay from A\pi,v=0 into X\Sigma,v=1, from New J Phys. 21 (2019) 052002
const v1BranchingRatioB = 1-.9711;#ratio of population decay from B\Sigma,v=0 into X\Sigma,v=1
const gam = 2 * pi * 6.4e6;#linewidth (happens to be same for B and A for SrF).  Haven't figure out a good way to implement differing gamma in bichromatic traps.  For now, assume this is small effect?
const normalizedBohrMag = 1.758820e7 / 2 / gam;
const mass = (40 + 16 + 1) * 1.67e-27;#mass of CaOH
const a = 0.999633;#j mixing terms a and b, calculated by me and also from J Phys. B: At. Mol. Opt. Phys 50 (2017) 015001
const b = sqrt(1 - a^2);
const gs = [-.33, 0.83, 0.5, -0.04,1.04];#g values (note: for any B-field work, should use 'zeeman' code, since hyperfine splitting is small).  First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
const stateEnergiesGround = [0.0, 0.0, 8.0, 8.2]; #X\Sigma hyperfine energies, with 0 corresponding to the F=1\DOWN energy
const stateEnergiesExcited = [0,-3.5];#Energy of F=0 relative to "0" (F=1).  Entry 1 for A state, Entry 2 for B State (can't actually find this anywhere so just guessing based on CaF). 

const kA = 2 * pi / lamA; #wavenumber
const kB = 2 * pi / lamB;
const kRepump = 2 * pi / lamRepump;
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2