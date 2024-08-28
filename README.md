# minimumWorkingExampleMOTSimulator
minimum working example for molecular MOTs & amp; also example of how to run it on cluster.

See `gettingStartedWithJuliaMoleculeSimulator.pdf'

## Install Julia on local computers
Follow the instructions [here](https://julialang.org/downloads/) to use terminals to install Julia (top of the page). It willl also install [juliaup](https://github.com/JuliaLang/juliaup#mac-and-linux), which is a version manager of Julia. Juliaup also adds Julia to PATH by default.

## Todo
1. clarify in write about the molecule velocity direction 
2. clarify about dragging molecule through laser field with constant v, random initial position and average, Hamiltonian period
3. what we actually save: f dot v / |v|, etc
4. using .= in OBE solver iteration function to avoid re-allocating memory. So pre-cached array may have different values every iteration. Check its value before using it. @. faster
5. save all settings, save all numbers in one file, save units, save in .csv
6. analysis code
7. extend to cluster
8. parallelize
9. longitudinal static B field
10. frequency round to 0.1