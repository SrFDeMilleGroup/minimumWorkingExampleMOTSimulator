const lamA=860e-9;#note: positions normalized to \tilde{x}=k_{SrF,X->A}x
const lamB=736.7e-9;#note: this has 18 MHz linewidth.  Will try to add in different linewidths later, if possible.  
const lamRepump = 896e-9;
const v1BranchingRatioA = 1-.9508;#ratio of population decay from A\pi,v=0 into X\Sigma,v=1, from Chin. J. Chem. Phys., Vol 35, No. 1, 58
const v1BranchingRatioB = 1-.81;#ratio of population decay from B\Sigma,v=0 into X\Sigma,v=1
const gam = 2 * pi * 3.0e6;#linewidth (happens to be same for B and A for SrF).  Haven't figure out a good way to implement differing gamma in bichromatic traps.  For now, assume this is small effect?
const normalizedBohrMag = 1.758820e7 / 2 / gam;
const mass = (138 + 19) * 1.67e-27;#mass of MgF
const a = 0.9593;#j mixing terms a and b, calculated by me and also from J Phys. B: At. Mol. Opt. Phys 50 (2017) 015001
const b = sqrt(1 - a^2);
const gs = [.015, 0.485, 0.5, -0.202,1.202];#g values.  First 3 are for X state F=1DOWN, F=1UP, F=2. 4th is for A state F=1, 5th is for B state F=1.
const stateEnergiesGround = [0.0, 9.3, 39.2, 50.6]; #X\Sigma hyperfine energies, with 0 corresponding to the F=1\DOWN energy
const stateEnergiesExcited = [-5.6,-3.1];#Energy of F=0 relative to "0" (F=1).  Entry 1 for A state (nl-edm group), Entry 2 for B State (can't actually find this anywhere). 

const kA = 2 * pi / lamA; #wavenumber
const kB = 2 * pi / lamB;
const kRepump = 2 * pi / lamRepump;
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2