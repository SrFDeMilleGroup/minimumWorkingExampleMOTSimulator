#!/bin/bash
#SBATCH --job-name=BLUEF1
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH -c 20
#SBATCH --mail-user=tlangin@uchicago.edu
#SBATCH --array=1-1

module load julia/1.7.2
molecules=(SrF);
aOrBOrABs=(A);
addRepumps=(no);
tests=(g);
for addRepump in ${addRepumps[@]};do
  for aOrBOrAB in ${aOrBOrABs[@]}; do
    for molecule in ${molecules[@]}; do
      for test in ${tests[@]}; do
        export JULIA_NUM_THREADS=20
        julia twoSigMoleculeSimulationHyperfineBasisXABLUE.jl $molecule $aOrBOrAB $addRepump $test
      done
    done
  done
done
