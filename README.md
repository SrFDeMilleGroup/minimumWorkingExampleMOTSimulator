# minimumWorkingExampleOBESimulator
Minimum working example of a master equation simulator for molecular slowing and magneto-optical traps (MOTs).

This README.md file provides a brief introduction to the programming side of the simulation. For the physics side and details of master equations, please see *mainOBEWriteup.pdf*.

If you find this repository helpful, please considering citing our paper: 
- Langin, T. K., & DeMille, D. (2023). [Toward improved loading, cooling, and trapping of molecules in magneto-optical traps.](https://iopscience.iop.org/article/10.1088/1367-2630/acc34d) New Journal of Physics, 25(4), 043005.


## Run the simulation
### Install julia on local computers
Follow the instructions [here](https://julialang.org/downloads/) to use terminals to install julia (top of the page). It willl also install [juliaup](https://github.com/JuliaLang/juliaup#mac-and-linux), which is a version manager of julia. Juliaup also adds julia to PATH by default.

Here, the code is developed in julia 1.10.4.

### Run on local computers
Run *mainSimulationCode.jl* file in julia. 

User settings can be modified in *simulationSettings/generalSettings.jl* and *simulationSettings/laserSettings.jl* (see next section for the meaning of each variable). Using 4 processors by default, it parallelizes the computation for different molecule initial displacements and velocities. Number of processors and molecule species can be modified in *mainSimulationCode.jl*. Simulation results and settings are saved in *savedData* folder.

### Run on clusters
To do


## Simulation setup
By "simulating slowing and MOTs for molecules", what we actually mean is to calculate the radiation force (and therefore acceleration) molecules experience in external laser and magnetic fields. Here we can calculate radiation forces under different conditions, including initial position, velocity, B-field, laser configurations, etc. But in this simulator, we don't evaluate the trajectory of molecules from the beam source to MOT (although it can be derived from radiation forces straightforwardly). To calculate radiation forces, we first solve master equations to obtain the [quasi-steady state](#quasi-steady-state) density matrix, and calculate forces as expectation values of time derivative of the momentum operator. 

Forces are generally three-dimensional, but for us, the two most interesting components are their projections on initial molecule displacement ($\vec{f}\cdot\vec{r}/\left|\vec{r}\right|$) and velocity ($\vec{f}\cdot\vec{v}/\left|\vec{v}\right|$). They are referred to as trapping and slowing/cooling forces respectively. And these are what we record at the end of simulation.

### Units in simulation
To make most variables in this simulation code dimensionless, we introduce the following unit system.
- $\hbar=1$
- unit of frequency: $\Gamma$, where $\Gamma$ is the linewidth of a molecular transition, usually X-A transition.
- unit of energy: $\hbar\Gamma$
- unit of time: $1/\Gamma$.
- unit of length: $1/k$ (= $\lambda/2\pi$), where $k$ is the (angular) wavenumber of molecular transition, usually X-A transition, and $\lambda$ is the corresponding wavelength.
- unit of velocity: $\Gamma/k$.
- unit of force: $10^{-3}\times\hbar\Gamma\times k$.

(Note that this is different from atomic units, where the unit of length is Bohr radius, and the unit of time is also different. The logic here is that we define three independent units: $\hbar$, unit of time, and unit of length first, and then derive the rest from them. Although not explicitely invloved, mass would also have a derived unit $\hbar k^2/\Gamma$.)

### Quasi-steady state
Because of molecule motion and laser fields, the Hamiltonian in the master equation is time-dependent and periodic, which therefore makes the solution time-dependent. The density matrix doesn't actually reach a steady state, but only a periodic quasi-steady state. The radiation forces are averaged over a period to generate final results. Period depends on two quantities, 
1. greatest common devisor (GCD) of three velocity components $v_x$, $v_y$, $v_z$, denoted as $v_0$;
2. GCD of all frequencies, including laser detuning and molecule energy level splittings, denoted as $\omega_0$.

The period will be the least common multiple (LCM) of $\lambda/v_0$ and $2\pi/\omega_0$. To express it under the unit system in ths simulation, it becomes the LCM of $2\pi/v_0$ and $2\pi/\omega_0$.

To make period deterministic, we round molecule velocities to a small number `vRound`, and all frequencies to `freqRound`, and use them as $v_0$ and $\omega_0$ respectively. (Although one can argue that they are not necessary the *greatest* common devisor, but except making us average forces over more than actual period and costing some computation resources, it does no harm to the final results.)

In the current setup of the simulation, `freqRound` and `vRound` are two constants hard-coded in *structsAndFunctions/simulateIt.jl*. `freqRound` is fixed at 0.1, and `vRound` has four possible values (0.002, 0.01, 0.02, 0.05) depending on molecule initial velocity. In this setup, the quasi-steady state period is $2\pi$/`vRound`. 

### Initial displacement and velocity
In ths evolution of master equations, molecules are considered to move in the external fields from initial displacement $\vec{r}_{\text{init}}$, with constant velocity $\vec{v}$. Molecule position is calculated at every step of the numerical evaluation of master equation as $\vec{r}\left(t\right)=\vec{r}_{\text{init}}+\vec{v}t$. 

Although users need to specify initial displacements as part of the simulation inputs, the actual $\vec{r}_{\text{init}}$ for computation could be randomly sampled from a cube of size of one wavelength around the user specified initial displacement to account for laser field variation. The velocities are also specified by users, but its direction could be randomly sampled up to some constraints. More details about the settings can be found in [Simulation settings II: general settings](#simulation-settings-ii-general-settings). Typically we repeat this process many times and average to get the final results. In julia, this forms a `EnsembleProblem` (from `DifferentialEquations` package).

### File hierarchy
Here we modularize the code as much as possible. Contents of most julia files (.jl) are wrapped up as modules, to avoid potential namespace conflicts. Julia has this inconvenient (or in my opinion, immature) module system that a local module can only be imported once, and has to be from the top-most level file that uses it. It poses limitations for us to organize files and modules.

Currently we organize the files in the following way.

- *mainSimulationCode.jl* (top-most level file) imports from
  - *structsAndFunctions/structs.jl*
  - *simulationSettings/moleculeVariables.jl* (defines molecular constants)
  - *structsAndFunctions/generateLaserSettings.jl*, which imports from
    - *simulationSettings/laserSettings.jl* (defines laser configurations)
  - *structsAndFunctions/generateGeneralSettings.jl*, which imports from
    - *simulationSettings/generalSettings.jl* (defines other configuratinos, like molecule initial displacements and velocities, etc.)
  - *structsAndFunctions/simulateIt.jl* (main master equation functions), which imports from
    - *structsAndFunctions/auxFunctions/obeInitialization.jl*
    - *structsAndFunctions/auxFunctions/obeEvaluation.jl*
    - *structsAndFunctions/auxFunctions/forceCalculation.jl*
  - *structsAndFunctions/saveSimulation.jl* (saves simulation results)

### Simulation settings I: laser settings
This simulation requries user inputs for laser settings, as listed here. They can be modified in *simulationSettings/laserSettings.jl*.
- `s0::Vector{Float64}`: Saturation parameter corresponding the peak laser intensity of a single laser pass. Every entry in the vector corresponds to one laser (same for all vector settings here). 
- `laserEnergy::Vector{Float64}`: In unit of $\Gamma$, laser detuning measured from the lowest ground hyperfine level.
- `polType::Vector{String}`: laser polarization and propagation direction, valid options are below.
  - '3D': 6 laser beams propagate in +/-x, +/-y, +/-z 6 directions respectively. All lasers are circularly polarized, but x and y lasers have opposite polarization handedness from z lasers, as required by 3D MOT quadrupole B field. Sigma +/- polarization is indicated by `polSign`. This is mainly used to simulate 3D MOT or molasses cooling.
  - '2DSS': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. All lasers are circularly polarized, but x lasers have opposite handedness from y lasers, as required by 2D MOT B field. Sigma +/- polarization is indicated by `polSign`. This is mainly used to simulate 2D MOT.
  - '2DPar': 4 laser beams propagate in +/-x, +/-y 4 directions respectively. All lasers are linear polarized along z axis. This is mainly used to simulate transverse cooling of a molecule beam, assuming it travels along z axis.
  - '2DPerp': Same as '2DPar', except x lasers are polarized along y, and y lasers are polarized along z axis.
  - 'Slower': 1 laser beam propagate along -z direction, and is linearly polarized along x axis.
  - 'Push': Same as 'Slower' except the laser propagates along +z direction.
- `polSign::Vector{Int64}`: If `polType` is 3D or 2DSS, -1/+1 determines sigma-/+ polarization for lasers, otherwise this settings is unused.
- `sidebandFreqs::Vector{Float64}`: In unit of $\Gamma$, modulation frequencies of EO modulators.
- `sidebandAmps::Vector{Float64}`: Radian, modulation depths of EO modulators.
- `whichTransition::Vector{String}`: The transition this laser addresses. Valid options are 'XA', 'XB' and 'XARepump'.
- `beamWaistInMM::Float64`: mm, laser beam waist. Used to simulate the effects of the finite size of 3D MOT laser beams. This is only used if `polType` is '3D', otherwise the laser beams are assumed infinitely large.

One thing we should point out is that here we are using different coordinate systems for 3D MOT/molasses and slowing/transverse cooling simulations. For 3D MOT, we assume z axis is the symmetry axis of the quadrupole B field, while for slowing/transverse cooling, we assume z axis is the molecule beam direction. They are usually perpendicular to each other, for most if not all molecule laser cooling expriments have been developed. This means that 3D MOT and slowing/transverse cooling can only be simulated separately ('3D' can't be present at the same time with other options in `polType`). Modification is needed if users want to simulate 3D MOT and slowing at the same time (e.g., to see slowing laser's effect on MOT trapped molecules). And note that, currently this simulation code doesn't prevent users from mixing `polType` together. It's users' responsibility to make sure input settings are sensible. 

Here we also present a typical set of values for 3D MOT simulation (SrF 5-frequency red-detuned DC MOT).
```julia
s0 = [10.4, 19.2, 10.4, 31.3, 8.7]
laserEnergy = [-1.0, -9.8, -18.6, -26.8, -20.8]
polType = ["3D", "3D", "3D", "3D", "3D"]
polSign = [1, 1, 1, -1, -1]
sidebandFreqs = [0., 0., 0., 0., 0.]
sidebandAmps = [0., 0., 0., 0., 0.]
whichTransition = ["XA", "XA", "XA", "XA", "XA"]
beamWaistInMM = 7.0
```

### Simulation settings II: general settings
All other simulation settings are summarized here. They can be modified in *simulationSettings/generalSettings.jl*.
- `simulationType::String`: A note to user themself about the simulation.
- `numTrialsPerValueSet::Int64`: Number of trials per set of values (displacementsInMM, userSpeeds, longSpeeds) to run and average. Averaging is needed because we randomly choose initial molecule displacements and velocities under user specified constraints.
- `forceProfile::String`: Valid options are below.
  - 'ThreeD': Used for 3D MOT and molasses simulations. In this case, `displacementsInMM` and `userSpeeds` are treated as magnitudes of 3D displacements and velocities respectively. `longSpeeds` is ignored. The direction of initial displacements is indicated by `initDispDir`. And the direction of velocity is set with respect to displacements by `velDirRelToR`. In order to account for laser field variation, the actual initial displacements in computation are randomly sampled from a cube of size of one wavelength round the values specified here. The projections of 3D force on initial displacement $\vec{f}\cdot\vec{r}/\left|\vec{r}\right|$ and velocity $\vec{f}\cdot\vec{v}/\left|\vec{v}\right|$ are calculated and returned. 
  - 'TwoD': Used for slowing and transverse cooling simulations. In this case, `displacementsInMM` and `userSpeeds` are treated as magnitudes of 2D displacements and velocities in x-y plane. Displacement along z is taken as zero (up to a random value within +/- 1/2 wavelength from 0). `longSpeeds` is used as velocity along z direction. `initDispDir` is ignored, and the direction of initial displacements in x-y plane is always random. The relative direction of velocity and displacement in x-y plane is set by `velDirRelToR`. The projection of 3D force on x-y plane displacement $\vec{f}\cdot\vec{r}/\left|\vec{r}\right|$ and x-y plane velocity $\vec{f}\cdot\vec{v}/\left|\vec{v}\right|$, as well as its z component $f_z$ are calculated and returned.
- `displacementsInMM::Vector{Float64}`: mm, magnitude of initial displacements. Simulation will iterate through the entire list.
- `initDispDir::String` : Direction of initial displacememnt. Only used if `forceProfile` is 'ThreeD'. Valid options are 'XY' ((x+y)/sqrt(2) direction, where slowed molecules come into MOT region for most experiments), 'Z' and 'Random'.
- `longSpeeds::Vector{Float64}`: In unit of $\Gamma/k$, a list of molecule beam longitudinal speeds to simulate through. Only used if `foreProfile` is 'TwoD'.
- `userSpeeds::Vector{Float64}`: In unit of $\Gamma/k$, a list of molecule speeds to simulate through.
- `velDirRelToR::String`: Relative direction of `userSpeeds` and `displacementsInMM`. Valid options are 'Same', 'Orthogonal', 'Opposite' and 'Random'.
- `bFieldSetting::String`: Set B field types. Valid options are
  - 'ThreeD': 3D quadrupole B field centered at the origin, as used for 3D MOT. `bGradReal` will have the unit Gauss/cm.
  - 'TwoD': 2D B field gradient centered at the origin, as used in 2D MOT. `bGradReal` will have the unit Gauss/cm.
  - 'StaticXY': Uniform B field (no gradient) along (x+y)/sqrt(2) direction, can be used for slowing, transverse cooling, 3D molasses simulation, etc. `bGradReal` will have the unit Gauss.
  - 'StaticZ': Same as 'StaticXY' but B field along z direction.
- `bGradReal::Float64`: Gauss/cm for B field gradient, or Gauss for uniform B field, depending on `bFieldSetting`.

Here we also present a typical set of values for 3D MOT simulation.
```julia
simulationType = "SrFRedMOTNormalValues"
numTrialsPerValueSet::Int64 = 100
forceProfile::String = "ThreeD"
displacementsInMM::Vector{Float64} = [0.5, 1.5, 3.0, 4.5, 6.0, 7.5]
initDispDir::String = "XY"
longSpeeds::Vector{Float64} = [32]
userSpeeds::Vector{Float64} = [-4, -3, -2, -1, -0.5, -0.1, -0.05, 0.05, 0.1, 0.5, 1, 2, 3, 4]
velDirRelToR::String = "Same"
bFieldSetting::String = "ThreeD"
bGradReal::Float64 = 12.5
```


## To-do
1. analysis code
2. extend to cluster
3. implement randomized laser phases
4. specify package version
5. implement machine learning to further optimize MOT and slowing