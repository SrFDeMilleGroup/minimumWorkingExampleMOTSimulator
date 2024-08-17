using DifferentialEquations
#using Plots
using BenchmarkTools
using Octavian
using LinearAlgebra
using SharedArrays
#using Distributed
#using WignerSymbols
#using Coverage
using Trapz
using DelimitedFiles
using Statistics
using Dates

struct Lasers{T1<:Vector{Float64},T2<:Vector{Int64},T3<:Vector{String},T4<:Vector{Matrix{Float64}}} #structure 'defining' a laser
    s0::T1;#saturation intensity at laser center (single pass)
    laserEnergy::T1;#energy of laser (note: zero energy defined to be energy of transition from |X\Sigma,F=1,J=1/2> to |F'=1>)
    polSign::T2;#polarization sign (for configurations using \sigma+/- light.  Defines if x-axis, say, is +\sigma or -\sigma (and corresponding changes to other axes...))
    whichTransition::T3;#"XB", "XA", or "XARepump"
    polType::T3;#= polType can be "3D" (sig +/-, with z-axis (quadrupole coil axis) reversed wrt other axes), "2DSS" (sig +/- but lasers only in x,y direction.  if \sig+ along +x then \sig- along +y).  
    "2DPar"(lasers in x,y direction both polarized along z).  "2DPerp" (x laser polarized along y, y polarized along z).  "Slower" (z laser linearly polarized along x) =#
    sidebandFreqs::T1;#frequency at which sidebands are driven
    sidebandAmps::T1;#phase modulation depth in radians
    wavenumberRatios::T1;#ratio of k_{Laser} to k_{A} 
    laserMasks::T4;#used in calculation of density matrix evolution.  Turns off coupling terms corresponding to, for example, X->B and X(v=1)->A for a laser with 'whichTransition'="XA"
end

function preInitializer(numLasers,numZeemanStatesGround,numZeemanStatesTotal)#initializes a bunch of stuff used in the OBE solver.  Julia likes things pre-initialized if possible

    #holds the modified coupling matrices used in decay terms

    coupleMatEff1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    coupleMatEff2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    coupleMatEff3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #convenient for fast evaluation of terms used in the 'decay' term of the density matrix evolution (second term in eq 1 of main writeup)

    decayMaskAllButTopLeft = zeros(Float64, numZeemanStatesTotal, numZeemanStatesTotal);
    decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1;
    decayMaskAllButTopLeft[1:numZeemanStatesGround, (numZeemanStatesGround+1):numZeemanStatesTotal] .= -1 / 2;
    decayMaskAllButTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, 1:numZeemanStatesGround] .= -1 / 2;
    decayMaskForCalcTopLeft = zeros(Int64, numZeemanStatesTotal, numZeemanStatesTotal);
    decayMaskForCalcTopLeft[(numZeemanStatesGround+1):numZeemanStatesTotal, (numZeemanStatesGround+1):numZeemanStatesTotal] .= 1;

    #now we make a bunch of initializations.  This makes the julia code run much faster at the cost of some readability...
    r = Array{Float64,1}(undef, 3)

    #fieldTerms[i] are the projections of the light field for laser[i] at a given position on the \sigma^-,\pi,\sigma^+ basis
    fieldTerms = Array{Array{ComplexF64,1},1}(undef, numLasers);
    for i = 1:numLasers
        fieldTerms[i] = vec(zeros(ComplexF64, 1, 3))
    end

    #will eventually 'hold' the atom-light matrix term of the hamiltonian during the diff-eq solver (see densityMatrixChangeTerms! in auxFunctions)
    atomLightTerm = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #bField terms are the projections of magnetic field at a given position on the \sigma^-,\pi,\sigma^+ basis. bFieldTermFull basically holds the 'mu' tensor
    bFieldTerms = Array{ComplexF64,1}(undef, 3);
    bFieldTermFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #will eventually hold the -\mu\cdotB (and hermitian conjugate) terms
    uProdBField = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    bFieldProdU = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #initializations of some matrices used to speed up the decay term calculation
    decayFull = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pOnlyExcitedStates = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft1PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft2PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft3PreMult = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft1 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft2 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);
    pTopLeft3 = zeros(ComplexF64, numZeemanStatesTotal, numZeemanStatesTotal);

    #return pre-initialized stuff from here to join the rest of the pre-initialized stuff in the 'main' program
    pPreInitialized = [coupleMatEff1,coupleMatEff2,coupleMatEff3,decayMaskAllButTopLeft, 
    decayMaskForCalcTopLeft, r, fieldTerms, atomLightTerm, bFieldTerms, bFieldTermFull, uProdBField, bFieldProdU, decayFull,
    pOnlyExcitedStates, pTopLeft1PreMult, pTopLeft2PreMult, pTopLeft3PreMult, pTopLeft1, pTopLeft2, pTopLeft3];
    return pPreInitialized;
end

function generateRandPosAndVel(forceProfile,numTrialsPerSpeed,velDirRelToR,currDisp,currSpeed,vRound,longSpeed,forceXY);
    #Function generates a set of random positions and 'pseudo'-random velocities (direction determined by 'velDirRelToR' + whether 'force profile' is 2D or 3D.)
    #if forceProfile is TwoD: z position is assumed to not matter, z velocity is fixed to longSpeed, and direction of velocity relative to random choice of \phi where x=disp*(cos(\phi)), etc. determined by velDirRelToR
    #if forceProfile is ThreeD: longSpeed isn't used, and direction of velocity chosen relative to random x,y,z direction of position is determined by velDirRelToR
    #velDirRelToR=-1 gives random orientation
    #velDirRelToR=0 forces v parallel to r
    #velDirRelToR=1 forces v perpendicular to r (or, for 2D, perpendicular in the xy plane at least)
    #velDirRelToR=2 forces v anti-parallel to r
    if forceProfile=="TwoD"
        #randomize position direction
        randPhisPos = rand(numTrialsPerSpeed, 1) * 2 * pi
        randRxs = cos.(randPhisPos) .* currDisp .* 1e-3 .* kA
        randRys = sin.(randPhisPos) .* currDisp .* 1e-3 .* kA
        randRzs = rand(numTrialsPerSpeed, 1) * 2 * pi

        if velDirRelToR==-1#randomize phi
            randPhisVels = rand(numTrialsPerSpeed, 1) * 2 * pi
        else#adjust phis from positions to force v either parallel, orthogonal, or anti-parallel to the r choice
            randPhisVels = randPhisPos .+ pi/2*velDirRelToR;#pi/2 for ortho, pi for total reversal, 0 for same.
        end
        randVxs = round.(currSpeed .* cos.(randPhisVels) ./ vRound) .* vRound
        randVys = round.(currSpeed .* sin.(randPhisVels) ./ vRound) .* vRound
        randVzs = round.(currSpeed .* cos.(randPhisVels) ./ vRound) .* 0 .+ longSpeed#
    else #if 3D
        #random position direction
        randX = randn(numTrialsPerSpeed, 1);
        randY = randn(numTrialsPerSpeed, 1);
        randZ = randn(numTrialsPerSpeed, 1);
        normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        randRxs = randX ./ normTerms .* currDisp .* 1e-3 .* kA;
        randRys = randY ./ normTerms .* currDisp .* 1e-3 .* kA;
        randRzs = randZ ./ normTerms .* currDisp .* 1e-3 .* kA;
        #forceXY forces position to be along (x+y)/sqrt(2) (e.g., entering from slower)
        if forceXY == 1
            randRxs = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRys = 1 ./ sqrt(2) .* currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRzs = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randX = randRxs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randY = randRys ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randZ = randRzs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        elseif forceXY == 2
            randRxs = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRys = 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randRzs = currDisp .* 1e-3 .* kA .+ 2 .* pi .* randn(numTrialsPerSpeed, 1);
            randX = randRxs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randY = randRys ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            randZ = randRzs ./ sqrt.(randRxs.^2 .+ randRys.^2 .+ randRzs.^2);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
        end
        if velDirRelToR == -1#random velocity direction as wel
            randX = randn(numTrialsPerSpeed, 1);#re-roll
            randY = randn(numTrialsPerSpeed, 1);
            randZ = randn(numTrialsPerSpeed, 1);
            normTerms = sqrt.(randX.^2 .+ randY.^2 .+ randZ.^2);
            randVxs = randX ./ normTerms .* currSpeed;
            randVys = randY ./ normTerms .* currSpeed;
            randVzs = randZ ./ normTerms .* currSpeed;
        elseif velDirRelToR == 0#same dir
            randVxs = randX ./ normTerms .* currSpeed;
            randVys = randY ./ normTerms .* currSpeed;
            randVzs = randZ ./ normTerms .* currSpeed;
        elseif velDirRelToR == 1#ortho dir
	    randX2 = randn(numTrialsPerSpeed,1);
	    randY2 = randn(numTrialsPerSpeed,1);
	    randZ2 = randn(numTrialsPerSpeed,1);
	    for i=1:length(randX2)
	        (randX2[i],randY2[i],randZ2[i]) =[randX2[i],randY2[i],randZ2[i]]-dot([randX[i],randY[i],randZ[i]],[randX2[i],randY2[i],randZ2[i]])./dot([randX[i],randY[i],randZ[i]],[randX[i],randY[i],randZ[i]]) .* [randX[i],randY[i],randZ[i]];
	    end
	    normTerms = sqrt.(randX2.^2 .+ randY2.^2 .+ randZ2 .^2);
            randVxs = randX2 ./ normTerms .* currSpeed;
            randVys = randY2 ./ normTerms .* currSpeed;
            randVzs = randZ2 ./ normTerms .* currSpeed;
        elseif velDirRelToR == 2#negative dir
            randVxs = -randX ./ normTerms .* currSpeed;
            randVys = -randY ./ normTerms .* currSpeed;
            randVzs = -randZ ./ normTerms .* currSpeed;
        end
        randVxs = round.(randVxs ./ vRound) .* vRound;
        randVys = round.(randVys ./ vRound) .* vRound;
        randVzs = round.(randVzs ./ vRound) .* vRound;
    end
    #run for both +/- r and +/- v (better statistics)
    randRxs = [randRxs; -randRxs];
    randRys = [randRys; -randRys];
    randRzs = [randRzs; -randRzs];
    randVxs = [randVxs; -randVxs];
    randVys = [randVys; -randVys];
    if forceProfile=="ThreeD"
        randVzs = [randVzs; -randVzs];
    else
        randVzs = [randVzs; randVzs];#for 2D, Vz is always "longSpeed"
    end
    
    for i=1:length(randVzs)#velocity along any dimension cannot be zero (particle should have x,y,z all change throughout OBE evolution to ensure periodicity)
        if randVxs[i]==0
            randVxs[i] = vRound .* sign.(randn(Float64));
        end
        if randVys[i]==0
            randVys[i] = vRound .* sign.(randn(Float64));
        end
        if randVzs[i]==0
            randVzs[i] = vRound .* sign.(randn(Float64));
        end
    end
    
    return randRxs,randRys,randRzs,randVxs,randVys,randVzs;

end

function createCouplingTermsandLaserMasks(whichTransition)
    #This function does a number of things
    #1) determine how many ground and excited states are needed (12 ground if no lasers are "XARepumps", 24 if there are repumps.  4 excited if only one of "A" or "B" are used, 8 if both are)

    #2) Based on this, write out "stateEnergyMatrix".  Ultimately this is subtracted from the laser energy in the OBE solver exp(-i*t*(energyDiff)) like term.  All columns are identical.  
    #2 (cont)) each row (i) is the energy of |i> relative to |F=1,J=1/2> (if i is a ground state) or |E,F'=1> for |i> corresponding to either E=A\Pi or E=B\Sigma.

    #3) make "Masks" for lasers based on what transition the laser corresponds to.  This is multiplied element-wise with coupling matrix in the OBE solver (densityMatrixChangeTerms).  This is zero
    # for terms that are not coupled together by the matrix (e.g., turns off X->A coupling for X->B laser, etc. and 1 for terms that are)

    #4) similarly, record wavenumber ratio based on what transition laser corresponds to.  

    #5) Establish 'coupling' (C matrices, eq 12-14 of writeup, basically 'clebsch-gordan' like terms) and 'b-coupling' matrices (C_B matrices, eq 27-29 of writeup.  Basically a 'B-field' coupling matrix based on g_{F} terms)
    #5 (cont)) Terms C_{i,j} are zero unless i=ground and j=excited.  size of matrix determined by number of excited states and ground states needed.  C_{i,j}[k] is the coupling from i->j for polarization k
    #5 (cont)) Terms C_{B,i,j} are zero unless i and j are in same F manifold.  size of matrix determined by number of excited states and ground states needed.  C_{B,i,j}[k] is the coupling from i->j for <B\cdot p_{k}>/|B|, where p_{k} is the \sigma^-/+,\pi basis

    #1)
    bichrom=0;#winds up 0 if only XA of XB are used, 1 if both are
    repump=0;#winds up 0 if no repump, 1 if there are repumps
    XToB = 0;#winds up 1 if only lasers are XB
    if "XB" in whichTransition
        if "XA" in whichTransition
            bichrom=1;
        else
            XToB = 1;
        end
    end
    if "XARepump" in whichTransition
        repump=1;
    end
    numZeemanStatesGround = 12+12*repump;
    numZeemanStatesExcited = 4+4*bichrom;
    numZeemanStatesTotal = numZeemanStatesGround+numZeemanStatesExcited;
    #2)
    stateEnergiesColumnFormat = [fill(stateEnergiesGround[1], 3); fill(stateEnergiesGround[2], 1); fill(stateEnergiesGround[3], 3); fill(stateEnergiesGround[4], 5)];#; fill(0, 4)];
    if repump==1
        stateEnergiesColumnFormat = vcat(stateEnergiesColumnFormat,stateEnergiesColumnFormat);#NOTE this assumes hyperfine splitting is the same in v=1 repump...not quite right but close enough
    end
    stateEnergiesColumnFormat = [stateEnergiesColumnFormat;fill(0,4+4*bichrom)];


    stateEnergyMatrix = repeat(stateEnergiesColumnFormat, 1, numZeemanStatesTotal);
    if XToB == 1 || bichrom==1
        stateEnergyMatrix[1:numZeemanStatesGround, end] = stateEnergyMatrix[1:numZeemanStatesGround, end] .- stateEnergiesExcited[2]; #handles excited state hyperfine splitting of |B\Sigma,F=0> level. 
        if bichrom==1
            stateEnergyMatrix[1:numZeemanStatesGround, end-4] = stateEnergyMatrix[1:numZeemanStatesGround, end-4] .- stateEnergiesExcited[1]; #handles excited state hyperfine splitting of |A\Pi,F=0> level. 
        end
    else
        stateEnergyMatrix[1:numZeemanStatesGround, end] = stateEnergyMatrix[1:numZeemanStatesGround, end] .- stateEnergiesExcited[1];#handles excited state hyperfine splitting of |A\Pi,F=0> level. 
    end

    #3+4)
    laserMasks = [zeros(numZeemanStatesTotal,numZeemanStatesTotal) for i=1:(length(whichTransition))]
    wavenumberRatios = Vector{Float64}(undef,length(whichTransition));
    for i=1:length(whichTransition)
        currTransition = whichTransition[i];
        if currTransition=="XA"
            laserMasks[i][1:12,(13+12*repump):(16+12*repump)] .= 1;
            wavenumberRatios[i]=1.0;
        elseif currTransition=="XB"
            laserMasks[i][1:12,(13+12*repump+4*bichrom):(16+12*repump+4*bichrom)] .= 1;
            wavenumberRatios[i] = kB/kA;
        else
            laserMasks[i][13:24,25:28] .= 1;
            wavenumberRatios[i] = kRepump/kA;
        end
    end
    #5)
    couplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)];

    makeCouplingMatrices!(couplingMatrices, a, b, XToB,repump,bichrom,v1BranchingRatioA,v1BranchingRatioB);

    bCouplingMatrices = Matrix[zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal), zeros(numZeemanStatesTotal, numZeemanStatesTotal)];

    makeBCouplingMatrices!(bCouplingMatrices, gs, XToB,repump,bichrom)


    return couplingMatrices,bCouplingMatrices,stateEnergyMatrix,laserMasks,wavenumberRatios,numZeemanStatesGround,numZeemanStatesExcited;
end

function checkErrors(bFieldSetting,forceProfile,whichTransition,polType)
    #checks if user made an error (invalid choice for laser transition type, or polarization type, etc.)

    for i = 1:length(whichTransition)
        (whichTransition[i]=="XA" || whichTransition[i]=="XB" || whichTransition[i]=="XARepump") ||
        throw(ArgumentError(string("invalid choice of ",whichTransition[i]," in whichTransition element ",i,". Valid options are XA, XB, or XARepump")))

        (polType[i]=="3D" || polType[i]=="2DSS" || polType[i]=="2DPar" || polType[i]=="2DPerp" || polType[i]=="Slower" || polType[i]=="Push") ||
        throw(ArgumentError(string("invalid choice of ",polType[i]," in polType element ",i,".  Valid options are 3D, 2DSS, 2DPar, 2DPerp, Slower")))
    end
    (bFieldSetting=="ThreeD" || bFieldSetting=="TwoD" || bFieldSetting=="Static") ||
    throw(ArgumentError(string("invalid choice of bFieldSetting, ",bFieldSetting,".  Valid options are ThreeD, TwoD, or Static")))

    (forceProfile=="ThreeD" || forceProfile=="TwoD") ||
    throw(ArgumentError(string("invalid choice of forceProfile, ",forceProfile,".  Valid options are ThreeD or TwoD")))
end

function makeCouplingMatrices!(couplingMatrices,a,b,XToB,repump,bichrom,v1BranchingRatioA,v1BranchingRatioB)
    #makes C_{i,j}[k] matrices.  What these look like depend on what ground/excited states are included
    #Choice 1) bichrom means that both A and B are 'spoken' to, and thus there are 8 excited states.
    #Choice 2) XToB=0 is true if no lasers 'talk' to B.  Thus, all 4 excited states are A states
    #Choice 3) Thus, if bichrom=0 and XToB=1, all 4 excited states are B states
    #In all cases, repump can be added, and thus there are 12 ground states.  These can only be coupled to the A state

    #these are all hardcoded for the assumption of a SrF, CaF, etc. type molecule where the alkaline has no hyperfine structure and, in the X\Sigma,N=1 state there
    #is mixing between 'pure' |F=1,J=1/2> and |F=1,J=3/2> that can be parameterized by a,b where |F=1,J~3/2> = a|F=1,J=3/2>+b|F=1,J=1/2> and |F=1,J~1/2> = -b|F=1,J=3/2>+a|F=1,J=1/2>
    #See Appendix A in writeup

    if bichrom ==1 #note: 12*repump term in second index forces 'excited' index to start at appropriate place, e.g. 13 for no repump, 25 if there is repump
        couplingMatrices[1][1,14+12*repump] = -sqrt(2)/3*a-b/6;
        couplingMatrices[1][1,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[1][2,15+12*repump] = -sqrt(2)/3*a-b/6;
        couplingMatrices[1][4,15+12*repump] = sqrt(2)/3;
        couplingMatrices[1][5,14+12*repump] = (a/6-sqrt(2)/3*b);
        couplingMatrices[1][5,16+12*repump] =(-a/3-sqrt(2)/3*b);
        couplingMatrices[1][6,15+12*repump]=(a/6-sqrt(2)/3*b);
        couplingMatrices[1][8,13+12*repump] = -1/sqrt(6);
        couplingMatrices[1][9,14+12*repump]=-1/(2*sqrt(3));
        couplingMatrices[1][10,15+12*repump]=-1/6;
        couplingMatrices[1][1,14+4+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[1][1,16+4+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[1][2,15+4+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[1][4,15+4+12*repump] = 1/3;
        couplingMatrices[1][5,14+4+12*repump] = (-a/3/sqrt(2)-b/3);
        couplingMatrices[1][5,16+4+12*repump] =(sqrt(2)*a/3-b/3);
        couplingMatrices[1][6,15+4+12*repump]=(-a/3/sqrt(2)-b/3);
        couplingMatrices[1][8,13+4+12*repump] = 1/sqrt(3);
        couplingMatrices[1][9,14+4+12*repump]=1/sqrt(6);
        couplingMatrices[1][10,15+4+12*repump]=1/3/sqrt(2);

        couplingMatrices[2][1,13+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[2][2,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[2][3,15+12*repump] = -sqrt(2)/3*a-1/6*b;
        couplingMatrices[2][4,14+12*repump] = -sqrt(2)/3;
        couplingMatrices[2][5,13+12*repump] = (-a/6+sqrt(2)/3*b);
        couplingMatrices[2][6,16+12*repump] =  (-a/3-sqrt(2)/3*b);
        couplingMatrices[2][7,15+12*repump] = (a/6-sqrt(2)/3*b);
        couplingMatrices[2][9,13+12*repump] = -1/(2*sqrt(3));
        couplingMatrices[2][10,14+12*repump] = -1/3;
        couplingMatrices[2][11,15+12*repump] = -1/(2*sqrt(3));
        couplingMatrices[2][1,13+4+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[2][2,16+4+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[2][3,15+4+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[2][4,14+4+12*repump] = -1/3;
        couplingMatrices[2][5,13+4+12*repump] = (a/3/sqrt(2)+b/3);
        couplingMatrices[2][6,16+4+12*repump] =  (sqrt(2)*a/3-b/3);
        couplingMatrices[2][7,15+4+12*repump] = (-a/3/sqrt(2)-b/3);
        couplingMatrices[2][9,13+4+12*repump] = 1/sqrt(6);
        couplingMatrices[2][10,14+4+12*repump] = sqrt(2)/3;
        couplingMatrices[2][11,15+4+12*repump] = 1/sqrt(6);

        couplingMatrices[3][2,13+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[3][3,14+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[3][3,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[3][4,13+12*repump] = sqrt(2)/3;
        couplingMatrices[3][6,13+12*repump] =  (-a/6+sqrt(2)/3*b);
        couplingMatrices[3][7,14+12*repump] = (-a/6+sqrt(2)/3*b);
        couplingMatrices[3][7,16+12*repump] =  (-a/3-sqrt(2)/3*b);
        couplingMatrices[3][10,13+12*repump] = -1/6;
        couplingMatrices[3][11,14+12*repump] = -1/(2*sqrt(3));
        couplingMatrices[3][12,15+12*repump] = -1/sqrt(6);
        couplingMatrices[3][2,13+4+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[3][3,14+4+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[3][3,16+4+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[3][4,13+4+12*repump] = 1/3;
        couplingMatrices[3][6,13+4+12*repump] =  (a/3/sqrt(2)+b/3);
        couplingMatrices[3][7,14+4+12*repump] = (a/3/sqrt(2)+b/3);
        couplingMatrices[3][7,16+4+12*repump] =  (sqrt(2)*a/3-b/3);
        couplingMatrices[3][10,13+4+12*repump] = 1/3/sqrt(2);
        couplingMatrices[3][11,14+4+12*repump] = 1/sqrt(6);
        couplingMatrices[3][12,15+4+12*repump] = 1/sqrt(3);

        if repump==1
            couplingMatrices[1][13:24,25:28]=couplingMatrices[1][1:12,25:28] .* sqrt(v1BranchingRatioA);
            couplingMatrices[1][13:24,29:32]=couplingMatrices[1][1:12,29:32] .* sqrt(v1BranchingRatioB);
            couplingMatrices[2][13:24,25:28]=couplingMatrices[2][1:12,25:28] .* sqrt(v1BranchingRatioA);
            couplingMatrices[2][13:24,29:32]=couplingMatrices[2][1:12,29:32] .* sqrt(v1BranchingRatioB);
            couplingMatrices[3][13:24,25:28]=couplingMatrices[3][1:12,25:28] .* sqrt(v1BranchingRatioA);
            couplingMatrices[3][13:24,29:32]=couplingMatrices[3][1:12,29:32] .* sqrt(v1BranchingRatioB);

            couplingMatrices[1][1:12,25:28]=couplingMatrices[1][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
            couplingMatrices[1][1:12,29:32]=couplingMatrices[1][1:12,29:32] .* sqrt(1-v1BranchingRatioB);
            couplingMatrices[2][1:12,25:28]=couplingMatrices[2][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
            couplingMatrices[2][1:12,29:32]=couplingMatrices[2][1:12,29:32] .* sqrt(1-v1BranchingRatioB);
            couplingMatrices[3][1:12,25:28]=couplingMatrices[3][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
            couplingMatrices[3][1:12,29:32]=couplingMatrices[3][1:12,29:32] .* sqrt(1-v1BranchingRatioB);
        end

    elseif XToB==0 #excited states are all "A" states
        couplingMatrices[1][1,14+12*repump] = -sqrt(2)/3*a-b/6;
        couplingMatrices[1][1,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[1][2,15+12*repump] = -sqrt(2)/3*a-b/6;
        couplingMatrices[1][4,15+12*repump] = sqrt(2)/3;
        couplingMatrices[1][5,14+12*repump] = (a/6-sqrt(2)/3*b);
        couplingMatrices[1][5,16+12*repump] =(-a/3-sqrt(2)/3*b);
        couplingMatrices[1][6,15+12*repump]=(a/6-sqrt(2)/3*b);
        couplingMatrices[1][8,13+12*repump] = -1/sqrt(6);
        couplingMatrices[1][9,14+12*repump]=-1/(2*sqrt(3));
        couplingMatrices[1][10,15+12*repump]=-1/6;
        
        couplingMatrices[2][1,13+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[2][2,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[2][3,15+12*repump] = -sqrt(2)/3*a-1/6*b;
        couplingMatrices[2][4,14+12*repump] = -sqrt(2)/3;
        couplingMatrices[2][5,13+12*repump] = (-a/6+sqrt(2)/3*b);
        couplingMatrices[2][6,16+12*repump] =  (-a/3-sqrt(2)/3*b);
        couplingMatrices[2][7,15+12*repump] = (a/6-sqrt(2)/3*b);
        couplingMatrices[2][9,13+12*repump] = -1/(2*sqrt(3));
        couplingMatrices[2][10,14+12*repump] = -1/3;
        couplingMatrices[2][11,15+12*repump] = -1/(2*sqrt(3));
        
        couplingMatrices[3][2,13+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[3][3,14+12*repump] = sqrt(2)/3*a+1/6*b;
        couplingMatrices[3][3,16+12*repump] = -sqrt(2)/3*a+b/3;
        couplingMatrices[3][4,13+12*repump] = sqrt(2)/3;
        couplingMatrices[3][6,13+12*repump] =  (-a/6+sqrt(2)/3*b);
        couplingMatrices[3][7,14+12*repump] = (-a/6+sqrt(2)/3*b);
        couplingMatrices[3][7,16+12*repump] =  (-a/3-sqrt(2)/3*b);
        couplingMatrices[3][10,13+12*repump] = -1/6;
        couplingMatrices[3][11,14+12*repump] = -1/(2*sqrt(3));
        couplingMatrices[3][12,15+12*repump] = -1/sqrt(6);

        if repump==1
            couplingMatrices[1][13:24,25:28]=couplingMatrices[1][1:12,25:28] .* sqrt(v1BranchingRatioA);
            couplingMatrices[2][13:24,25:28]=couplingMatrices[2][1:12,25:28] .* sqrt(v1BranchingRatioA);
            couplingMatrices[3][13:24,25:28]=couplingMatrices[3][1:12,25:28] .* sqrt(v1BranchingRatioA);

            couplingMatrices[1][1:12,25:28]=couplingMatrices[1][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
            couplingMatrices[2][1:12,25:28]=couplingMatrices[2][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
            couplingMatrices[3][1:12,25:28]=couplingMatrices[3][1:12,25:28] .* sqrt(1-v1BranchingRatioA);
        end

   else#excited states are all b states
        couplingMatrices[1][1,14+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[1][1,16+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[1][2,15+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[1][4,15+12*repump] = 1/3;
        couplingMatrices[1][5,14+12*repump] = (-a/3/sqrt(2)-b/3);
        couplingMatrices[1][5,16+12*repump] =(sqrt(2)*a/3-b/3);
        couplingMatrices[1][6,15+12*repump]=(-a/3/sqrt(2)-b/3);
        couplingMatrices[1][8,13+12*repump] = 1/sqrt(3);
        couplingMatrices[1][9,14+12*repump]=1/sqrt(6);
        couplingMatrices[1][10,15+12*repump]=1/3/sqrt(2);
        
        couplingMatrices[2][1,13+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[2][2,16+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[2][3,15+12*repump] = -a/3+b/3/sqrt(2);
        couplingMatrices[2][4,14+12*repump] = -1/3;
        couplingMatrices[2][5,13+12*repump] = (a/3/sqrt(2)+b/3);
        couplingMatrices[2][6,16+12*repump] =  (sqrt(2)*a/3-b/3);
        couplingMatrices[2][7,15+12*repump] = (-a/3/sqrt(2)-b/3);
        couplingMatrices[2][9,13+12*repump] = 1/sqrt(6);
        couplingMatrices[2][10,14+12*repump] = sqrt(2)/3;
        couplingMatrices[2][11,15+12*repump] = 1/sqrt(6);
        
        couplingMatrices[3][2,13+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[3][3,14+12*repump] = a/3-b/3/sqrt(2);
        couplingMatrices[3][3,16+12*repump] = -a/3-sqrt(2)*b/3;
        couplingMatrices[3][4,13+12*repump] = 1/3;
        couplingMatrices[3][6,13+12*repump] =  (a/3/sqrt(2)+b/3);
        couplingMatrices[3][7,14+12*repump] = (a/3/sqrt(2)+b/3);
        couplingMatrices[3][7,16+12*repump] =  (sqrt(2)*a/3-b/3);
        couplingMatrices[3][10,13+12*repump] = 1/3/sqrt(2);
        couplingMatrices[3][11,14+12*repump] = 1/sqrt(6);
        couplingMatrices[3][12,15+12*repump] = 1/sqrt(3);

        if repump==1#NOTE, there's really no reason this should ever execute...B and the vibrational repump are decoupled.  force this to not happen in main program.
            couplingMatrices[1][13:24,25:28]=couplingMatrices[1][1:12,25:28] .* sqrt(0);
            couplingMatrices[2][13:24,25:28]=couplingMatrices[2][1:12,25:28] .* sqrt(0);
            couplingMatrices[3][13:24,25:28]=couplingMatrices[3][1:12,25:28] .* sqrt(0);

        end

   end
end

function makeBCouplingMatrices!(bCouplingMatrices, gs, XToB,repump,bichrom)
   #describes magnetic field induced larmor precession (for 'perpendicular' fields with-respect-to magnetic moment) and energy shifts (for parallel fields).  Depends on g factor for given hyperfine state
   
    bCouplingMatrices[1][2,1]=gs[1];
    bCouplingMatrices[1][3,2]=gs[1];
    bCouplingMatrices[1][6,5]=gs[2];
    bCouplingMatrices[1][7,6]=gs[2];
    bCouplingMatrices[1][9,8] = sqrt(2)*gs[3];
    bCouplingMatrices[1][10,9] = sqrt(3)*gs[3];
    bCouplingMatrices[1][11,10] = sqrt(3)*gs[3];
    bCouplingMatrices[1][12,11] = sqrt(2)*gs[3];
    if bichrom==1
        bCouplingMatrices[1][14+12*repump,13+12*repump] = gs[4];
        bCouplingMatrices[1][15+12*repump,14+12*repump] = gs[4];
        bCouplingMatrices[1][14+12*repump+4,13+12*repump+4] = gs[5];
        bCouplingMatrices[1][15+12*repump+4,14+12*repump+4] = gs[5];
    elseif XToB==1
        bCouplingMatrices[1][14+12*repump,13+12*repump] = gs[5];
        bCouplingMatrices[1][15+12*repump,14+12*repump] = gs[5];
    else
        bCouplingMatrices[1][14+12*repump,13+12*repump] = gs[4];
        bCouplingMatrices[1][15+12*repump,14+12*repump] = gs[4];
    end

    bCouplingMatrices[2][1,1]=-gs[1];
    bCouplingMatrices[2][3,3]=gs[1];
    bCouplingMatrices[2][5,5]=-gs[2];
    bCouplingMatrices[2][7,7]=gs[2];
    bCouplingMatrices[2][8,8]=-2*gs[3];
    bCouplingMatrices[2][9,9] = -gs[3];
    bCouplingMatrices[2][11,11] = gs[3];
    bCouplingMatrices[2][12,12] = 2*gs[3];
    if bichrom==1
        bCouplingMatrices[2][13+12*repump,13+12*repump] = -gs[4];
        bCouplingMatrices[2][15+12*repump,15+12*repump] = gs[4];
        bCouplingMatrices[2][13+12*repump+4,13+12*repump+4] = -gs[5];
        bCouplingMatrices[2][15+12*repump+4,15+12*repump+4] = gs[5];
    elseif XToB==1
        bCouplingMatrices[2][13+12*repump,13+12*repump] = -gs[5];
        bCouplingMatrices[2][15+12*repump,15+12*repump] = gs[5];
    else
        bCouplingMatrices[2][13+12*repump,13+12*repump] = -gs[4];
        bCouplingMatrices[2][15+12*repump,15+12*repump] = gs[4];
    end

    bCouplingMatrices[3][1,2]=-gs[1];
    bCouplingMatrices[3][2,3]=-gs[1];
    bCouplingMatrices[3][5,6]=-gs[2];
    bCouplingMatrices[3][6,7]=-gs[2];
    bCouplingMatrices[3][8,9] = -sqrt(2)*gs[3];
    bCouplingMatrices[3][9,10] = -sqrt(3)*gs[3];
    bCouplingMatrices[3][10,11] = -sqrt(3)*gs[3];
    bCouplingMatrices[3][11,12] = -sqrt(2)*gs[3];
    if bichrom==1
        bCouplingMatrices[3][13+12*repump,14+12*repump] = -gs[4];
        bCouplingMatrices[3][14+12*repump,15+12*repump] = -gs[4];
        bCouplingMatrices[3][13+12*repump+4,14+12*repump+4] = -gs[5];
        bCouplingMatrices[3][14+12*repump+4,15+12*repump+4] = -gs[5];
    elseif XToB==1
        bCouplingMatrices[3][13+12*repump,14+12*repump] = -gs[5];
        bCouplingMatrices[3][14+12*repump,15+12*repump] = -gs[5];
    else
        bCouplingMatrices[3][13+12*repump,14+12*repump] = -gs[4];
        bCouplingMatrices[3][14+12*repump,15+12*repump] = -gs[4];
    end
    
    if repump==1
        bCouplingMatrices[1][13:24,13:24]=bCouplingMatrices[1][1:12,1:12];
        bCouplingMatrices[2][13:24,13:24]=bCouplingMatrices[2][1:12,1:12];
        bCouplingMatrices[3][13:24,13:24]=bCouplingMatrices[3][1:12,1:12];
    end

end

function propR!(r, rInit::Array{Float64,1}, v::Array{Float64,1}, t::Float64)
    r[1] = rInit[1] + v[1] * t;
    r[2] = rInit[2] + v[2] * t;
    r[3] = rInit[3] + v[3] * t;
end

function makeFieldTerms!(fieldTerms, r::Array{Float64,1}, polSign::Array{Int64,1}, polType::Array{String,1}, wavenumberRatios::Array{Float64,1},waist::Float64) #These are different for 2D MOT
    #returns 'field terms' for all lasers.  Field terms depend on the polarization type (and, if \sigma +/-, the sign)
    #This is basically the field for laser [i] due to the 6 (if 3D), or 1 (for Slower/push), or 4 (for all 2D lasers) passes of the beam expressed in the standard \sigma^- ([i][1]), \pi ([i][2]) and \sigma^+ ([i][3])
    #This is calculated in the way illustrated in JOSAB 6(11) 2023-2045 (1989) by Cohen-Tannoudji + Dalibard section 2.  See also Eq15-16 and subsequent expressions in my writeup for the 3D example.

    #IMPORTANT CAVEAT: all terms are 'pre-conjugated' since only the complex conjugate of this term is ever used (Eq 21 of my writeup).  Better to just express it pre-conjugated instead of 
    #repeatedly taking conjugates in the diff-eq solver

    for i=1:length(polSign)#iterate through all lasers
        if polType[i]=="Slower"#polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
            fieldTerms[i][1] = 1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) + im * sin(r[3] * wavenumberRatios[i]));
            fieldTerms[i][2] = 0;
            fieldTerms[i][3] = -1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) + im * sin(r[3] * wavenumberRatios[i]));

        elseif polType[i]=="Push"
            fieldTerms[i][1] = 1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) - im * sin(r[3] * wavenumberRatios[i]));
            fieldTerms[i][2] = 0;
            fieldTerms[i][3] = -1/sqrt(2) * (cos(r[3] * wavenumberRatios[i]) - im * sin(r[3] * wavenumberRatios[i]));

        elseif polType[i]=="2DSS"
            fieldTerms[i][1] = polSign[i] * sin(r[1] * wavenumberRatios[i]) + im * cos(r[2] * wavenumberRatios[i]);
        
            fieldTerms[i][2] = sqrt(2) * im * (cos(r[1] * wavenumberRatios[i]) - polSign[i]*sin(r[2] * wavenumberRatios[i]));
        
            fieldTerms[i][3] = polSign[i] * sin(r[1] * wavenumberRatios[i]) - im * cos(r[2] * wavenumberRatios[i]);

        elseif polType[i]=="2DPerp"
            fieldTerms[i][1] = -sqrt(2) * im * cos(r[1] * wavenumberRatios[i]);
            fieldTerms[i][2] = 2 * cos(r[2] * wavenumberRatios[i]);
            fieldTerms[i][3] = -sqrt(2) * im * cos(r[1] * wavenumberRatios[i]);
        elseif polType[i]=="2DPar"
            fieldTerms[i][1] = 0;
            fieldTerms[i][2] = 2 * (cos(r[2] * wavenumberRatios[i])+cos(r[1] * wavenumberRatios[i]));
            fieldTerms[i][3]=0;

        else#add 2Dperp,2DPar later
            fieldTerms[i][1] = cos(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) .+ polSign[i] * sin(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) .-
             im * (polSign[i] * sin(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) .- cos(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2));
        
            fieldTerms[i][2] = sqrt(2) * im * (cos(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) .+
             polSign[i] * sin(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2));
        
            fieldTerms[i][3] = cos(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) .+ polSign[i] * sin(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) .+
             im * (polSign[i] * sin(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) - cos(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2));

        end
    end#end polSign (e.g. end iteration through lasers)
end

    
function makeBFieldTerms!(bFieldTerms, r::Array{Float64,1},bFieldSetting::String)
    #expresses B field at position r in the \sigma^+/-, pi basis
    # r=zeros(1,3);
    
    # BFieldPTerms = zeros(ComplexF64,1,3);
    if bFieldSetting=="TwoD"
        bFieldTerms[1] = 1 / sqrt(2) * (r[1] - im * r[2]);
        bFieldTerms[2] = -0 * (r[3]);
        bFieldTerms[3] = 1 / sqrt(2) * (-r[1] - im * r[2]);
    elseif bFieldSetting=="ThreeD"
        bFieldTerms[1] = 1 / sqrt(2) * (r[1] + im * r[2]);
        bFieldTerms[2] = -1 * (r[3]);#note: really this should be -2r[3] for a quadropole field.  In practice, I prefer to run my f(r) for random direction at constant B.  So, assume \tilde{r}=(x,y,z/2).
        bFieldTerms[3] = 1 / sqrt(2) * (-r[1] + im * r[2]);
    else#if Static
        #bFieldTerms[1] = im/2;
        #bFieldTerms[2] = 1/sqrt(2);
        #bFieldTerms[3] = im/2;
        bFieldTerms[1] = (im+1)/2;
        bFieldTerms[2] = 0;
        bFieldTerms[3] = (-1+im)/2;
    end
    # return BFieldPTerms
end

function densityMatrixChangeTerms!(du, u, p, t)
    #The meat of the program.  Here's where the density matrix is actually evolved.

    # user inputs (these vary, things like initial position, velocity, laser params, etc.).  These are determined by the user-chosen parameters in the main program
    rInit = p[1]::Vector{Float64};
    v = p[2]::Vector{Float64};
    stateEnergyMatrix = p[3]::Matrix{Float64};
    lasers = p[4]::Lasers{Vector{Float64}, Vector{Int64}, Vector{String}, Vector{Matrix{Float64}}};
    laserEnergy = lasers.laserEnergy;
    s0 = lasers.s0;
    polSign = lasers.polSign;
    polType = lasers.polType;
    sidebandFreqs = lasers.sidebandFreqs;
    sidebandAmps = lasers.sidebandAmps;
    laserMasks = lasers.laserMasks;
    wavenumberRatios = lasers.wavenumberRatios;
    waist = p[5]::Float64;
    bToHamConvert = p[6]::Float64;

    # coupling matrices passed by user.  
    coupleMat1 = p[7]::Matrix{Float64};
    coupleMat2 = p[8]::Matrix{Float64};
    coupleMat3 = p[9]::Matrix{Float64};
    bCoupleMat1 = p[10]::Matrix{Float64};
    bCoupleMat2 = p[11]::Matrix{Float64};
    bCoupleMat3 = p[12]::Matrix{Float64};
    # coupling matrices used in decay calc
    coupleMatEff1 = p[13]::Matrix{ComplexF64};
    coupleMatEff2 = p[14]::Matrix{ComplexF64};
    coupleMatEff3 = p[15]::Matrix{ComplexF64};

    # decay 'masks' used in calculating the decay term.  
    decayMaskAllButTopLeft = p[16]::Matrix{Float64};
    decayMaskForCalcTopLeft = p[17]::Matrix{Int64};

    # pre-cached r Array
    r = p[18]::Vector{Float64};
    # pre-cached matrices for atom light term.  
    fieldTerms = p[19]::Vector{Vector{ComplexF64}};
    atomLightTerm = p[20]::Matrix{ComplexF64};
    atomLightTerm = zeros(ComplexF64, size(coupleMat1,1), size(coupleMat1,2));
    
    # pre-cached matrices for b field term. 
    bFieldTerms = p[21]::Vector{ComplexF64};
    bFieldTermFull = p[22]::Matrix{ComplexF64} ;
    uProdBField = p[23]::Matrix{ComplexF64} ;
    bFieldProdU = p[24]::Matrix{ComplexF64} ;

    # pre-cached matrices for decay term
    decayFull = p[25]::Matrix{ComplexF64};
    pOnlyExcitedStates = p[26]::Matrix{ComplexF64};
    pTopLeft1PreMult = p[27]::Matrix{ComplexF64};
    pTopLeft2PreMult = p[28]::Matrix{ComplexF64};
    pTopLeft3PreMult = p[29]::Matrix{ComplexF64};
    pTopLeft1 = p[30]::Matrix{ComplexF64};
    pTopLeft2 = p[31]::Matrix{ComplexF64};
    pTopLeft3 = p[32]::Matrix{ComplexF64};

    #Whether B-field is 3D, 2D, or static

    bFieldSetting = p[33]::String;


    #1) evolve position
    propR!(r, rInit, v, t);
    #2) Calculate field terms at new position
    makeFieldTerms!(fieldTerms, r, polSign,polType,wavenumberRatios,waist);
    #3)calculate -E dot D term (see Eq 21 of writeup)
    for i=1:length(s0)
        atomLightTerm .= atomLightTerm.+ sqrt(s0[i] / 8) .* -exp(1im * laserEnergy[i] * t + 1im * sidebandAmps[i]*sin(sidebandFreqs[i] * t)) .* laserMasks[i] .* ((fieldTerms[i][1] .* coupleMat1) .+
         (fieldTerms[i][2] .* coupleMat2) .+ (fieldTerms[i][3] .* coupleMat3));
    end
    atomLightTerm .= atomLightTerm.*exp.(1im*t.*stateEnergyMatrix);#subtracts relevant hyperfine energies from 'laserEnergy'
    atomLightTerm .= atomLightTerm .+ atomLightTerm';#needed here because, the way coupleMat is defined, 'atomLightTerm' up til now only has the top right half of the hermitian coupling matrix

    #4) calculate -mu dot B term (see Eq 32 of writeup)
    makeBFieldTerms!(bFieldTerms, r, bFieldSetting);
    bFieldTermFull .= bToHamConvert .* (bFieldTerms[1] .* bCoupleMat1 .+ bFieldTerms[2] .* bCoupleMat2 .+ bFieldTerms[3] .* bCoupleMat3).+atomLightTerm; #'bTermFull' also sums the -mu dot B term with the calculated -D dot E term
    #5) take commutator of [H,u] where u is density matrix and H = -D dot E + -mu dot B
    mul!(uProdBField, u, bFieldTermFull);#uH
    mul!(bFieldProdU, bFieldTermFull, u);#Hu

    #6) Take decay (aka 'coupling to reservoir') into account (Eq 46 of writeup)
    pOnlyExcitedStates .= u .* decayMaskForCalcTopLeft;

    #6A) these next 6 lines calculate the last term in eq 46 of writeup
    coupleMatEff1 = coupleMat1 .* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft1PreMult, coupleMatEff1, pOnlyExcitedStates);
    mul!(pTopLeft1, pTopLeft1PreMult, coupleMatEff1')

    coupleMatEff2 = coupleMat2 .* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft2PreMult, coupleMatEff2, pOnlyExcitedStates);
    mul!(pTopLeft2, pTopLeft2PreMult, coupleMatEff2')

    coupleMatEff3 = coupleMat3 .* exp.(-1im*t.*stateEnergyMatrix);
    mul!(pTopLeft3PreMult, coupleMatEff3, pOnlyExcitedStates);
    mul!(pTopLeft3, pTopLeft3PreMult, coupleMatEff3')

    decayFull .= (u .* decayMaskAllButTopLeft) .+ pTopLeft1 .+ pTopLeft2 .+ pTopLeft3;#u.*decayMask term represents 1st and 2nd term of eq 46 in writeup

    du .= 1im .* (uProdBField .- bFieldProdU) .+ decayFull;#finally, add the 'Liouville' term and the decay term (Eq 1 of writeup) to step the density matrix
end

#everything else here is used to calculate a force given a density matrix, see section 1.6 of writeup

function makeDFieldTerms!(dFieldTerms, r::Array{Float64,1}, polSign::Array{Int64,1}, polType::Array{String,1}, wavenumberRatios::Array{Float64,1},waist::Float64)
    # dfieldterms (dE/dr) will be 3x3 matrix, first element is xyz second is sig+ pi sig-
    # fieldTerms = zeros(ComplexF64,3,1);
    # pre conjugated, just like 'makeFieldTerms'
    for i=1:length(polSign) #iterates through all lasers
        if polType[i]=="Slower" #polarization vs position for a given laser depends on whether it a "3D\sig\sig", "2D\sig\sig", slower, etc.
            dFieldTerms[i][1,1] = 0;
            dFieldTerms[i][1,2] = 0;
            dFieldTerms[i][1,3] = 0;

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = 0;
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) + im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i];
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = -1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) + im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i];

        elseif polType[i]=="Push"
            dFieldTerms[i][1,1] = 0;
            dFieldTerms[i][1,2] = 0;
            dFieldTerms[i][1,3] = 0;

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = 0;
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) - im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i];
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = -1/sqrt(2) * (-sin(r[3] * wavenumberRatios[i]) - im * cos(r[3] * wavenumberRatios[i])) * wavenumberRatios[i];

        elseif polType[i]=="2DSS"
            dFieldTerms[i][1,1] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * wavenumberRatios[i];
            dFieldTerms[i][1,2] = -sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i];
            dFieldTerms[i][1,3] = polSign[i] * cos(r[1] * wavenumberRatios[i]) * wavenumberRatios[i];

            dFieldTerms[i][2,1] = -im * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i];
            dFieldTerms[i][2,2] = -sqrt(2) * im * (polSign[i] * cos(r[2] * wavenumberRatios[i])) * wavenumberRatios[i];
            dFieldTerms[i][2,3] = im * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i];

            dFieldTerms[i][3,1] = 0;
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = 0;

        elseif polType[i]=="2DPerp"
            dFieldTerms[i][1,1] = sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i];
            dFieldTerms[i][1,2] = 0;
            dFieldTerms[i][1,3] = sqrt(2) * im * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i];

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = -2 * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 0;
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = 0;
        elseif polType[i]=="2DPar"
            dFieldTerms[i][1,1] = 0;
            dFieldTerms[i][1,2] = -2 * sin(r[1] * wavenumberRatios[i]) * wavenumberRatios[i]
            dFieldTerms[i][1,3]=0;

            dFieldTerms[i][2,1] = 0;
            dFieldTerms[i][2,2] = -2 * sin(r[2] * wavenumberRatios[i]) * wavenumberRatios[i]
            dFieldTerms[i][2,3] = 0;

            dFieldTerms[i][3,1] = 0;
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = 0;
        else
            dFieldTerms[i][1,1] = polSign[i] * cos(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i];
            dFieldTerms[i][1,2] = -sqrt(2) * im * sin(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i];
            dFieldTerms[i][1,3] = polSign[i] * cos(r[1] * wavenumberRatios[i])*exp(-2*(r[2]^2+r[3]^2)/waist^2) * wavenumberRatios[i];
    
            dFieldTerms[i][2,1] = -im * sin(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2) * wavenumberRatios[i];
            dFieldTerms[i][2,2] = sqrt(2) * im * (polSign[i] * cos(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2)) * wavenumberRatios[i];
            dFieldTerms[i][2,3] = im * (sin(r[2] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[3]^2)/waist^2)) * wavenumberRatios[i];
    
            dFieldTerms[i][3,1] = (-sin(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) - im * (polSign[i] * cos(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2))) * wavenumberRatios[i];
            dFieldTerms[i][3,2] = 0;
            dFieldTerms[i][3,3] = (-sin(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2) + im * (polSign[i] * cos(r[3] * wavenumberRatios[i])*exp(-2*(r[1]^2+r[2]^2)/waist^2))) * wavenumberRatios[i];
        end
    end#end polSign

    # return fieldTerms
end

function forceCalc!(force, dFieldTerms::Vector{Matrix{ComplexF64}}, rho::Matrix{ComplexF64}, 
    lasers::Lasers{Vector{Float64}, Vector{Int64}, Vector{String}, Vector{Matrix{Float64}}}, couplingMatrices::Vector{Matrix}, stateEnergyMatrix::Matrix{Float64}, t::Float64)
    #calculates force given position and lasers (used to calculate dFieldTerms) and density matrix \rho.  
    #Both r(t) and \rho(t) are recorded vs time by the OBE solver, so this runs afterwords to calculate what forces the particle experienced over the trajectory
    s0 = lasers.s0;
    sidebandFreqs = lasers.sidebandFreqs;
    sidebandAmps = lasers.sidebandAmps;
    laserMasks = lasers.laserMasks;
    laserEnergy = lasers.laserEnergy;

    #force pre-factor is calculated for each laser [i].  Has the rotating-frame frequency exponent + phase modulation term + intensity term \sqrt(s0/8).  Hyperfine energies are subtracted later
    forcePrefactor = zeros(ComplexF64,1,length(laserEnergy));
    for i=1:length(laserEnergy)
        forcePrefactor[i] = sqrt(s0[i] / 8) * exp(1im * laserEnergy[i] * t + 1im * sidebandAmps[i]*sin(sidebandFreqs[i] * t));
    end

    #calculate x force.  Implements Eq 48 of main writeup
    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    dRhoDPosTimesDensityMatContainer = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][1,1] * couplingMatrices[1] .+ dFieldTerms[i][1,2] * couplingMatrices[2] .+ dFieldTerms[i][1,3] * couplingMatrices[3]) .* laserMasks[i];
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix.*exp.(1im*t.*stateEnergyMatrix);
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);#multiplies dp_{x}/dt by density matrix \rho
    force[1] = real(tr(dRhoDPosTimesDensityMatContainer));#takes trace of \rho*dp_{x}/dt to determine average force over enemble (Eq 49 of writeup)

    #similarly, calculate y and z force
    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][2,1] * couplingMatrices[1] .+ dFieldTerms[i][2,2] * couplingMatrices[2] .+ dFieldTerms[i][2,3] * couplingMatrices[3]) .* laserMasks[i];
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix.*exp.(1im*t.*stateEnergyMatrix);
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);
    force[2] = real(tr(dRhoDPosTimesDensityMatContainer));

    dRhoDPosCalcMatrix = zeros(ComplexF64, size(rho,1), size(rho,2));
    for i=1:length(laserEnergy)
        dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ forcePrefactor[i] * (dFieldTerms[i][3,1] * couplingMatrices[1] .+ dFieldTerms[i][3,2] * couplingMatrices[2] .+ dFieldTerms[i][3,3] * couplingMatrices[3]) .* laserMasks[i];
    end
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix.*exp.(1im*t.*stateEnergyMatrix);
    dRhoDPosCalcMatrix .= dRhoDPosCalcMatrix .+ dRhoDPosCalcMatrix';
    mul!(dRhoDPosTimesDensityMatContainer, rho, dRhoDPosCalcMatrix);
    force[3] = real(tr(dRhoDPosTimesDensityMatContainer));
end#forceCalc!

function makeForceVsTime!(forceVsTime, times::Vector{Float64}, rhos::Vector{Matrix{ComplexF64}},
    lasers::Lasers{Vector{Float64}, Vector{Int64}, Vector{String}, Vector{Matrix{Float64}}}, couplingMatrices::Vector{Matrix}, 
    stateEnergyMatrix::Matrix{Float64},waist::Float64, rInit::Vector{Float64}, v::Vector{Float64})
    #given a set of times, an initial position and velocity, the lasers used, and \rho(t), calculate force vs t
    polSign = lasers.polSign;
    polType = lasers.polType;
    wavenumberRatios = lasers.wavenumberRatios;
    #initialize some stuff
    dFieldContainer = Array{Array{ComplexF64, 2},1}(undef,length(polSign));
    for i=1:length(polSign)
        dFieldContainer[i]=zeros(ComplexF64,3,3)
    end#end polSign
    forceCalcContainer = zeros(ComplexF64, 3, 1);
    r = Array{Float64,1}(undef, 3)
    #iterate through time, propegating r in the same way done in the OBEs.  Then determine force experienced given r(t), \rho(t), and the lasers used
    for i = 1:length(times)
        propR!(r, rInit, v, times[i]);
        makeDFieldTerms!(dFieldContainer, r, polSign,polType, wavenumberRatios,waist)
        forceCalc!(forceCalcContainer, dFieldContainer, rhos[i], lasers, couplingMatrices, stateEnergyMatrix,times[i]);
        forceVsTime[i,:] = forceCalcContainer;
    end#end times

end#makeForceVsTime!
