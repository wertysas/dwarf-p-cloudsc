#!/bin/bash
# The job name
#SBATCH --job-name=dwarf-cloudsc-gpu-double-block
# Set the error and output files to [job name]-[JobID].out
#SBATCH --output=%x-%j.out 
#SBATCH --error=%x-%j.out
# Wall clock time limit
#SBATCH --time=00:05:00
# Send an email on failure
#SBATCH --mail-type=FAIL
# REQUIRED node settings/stats
#SBATCH --mem=200G
#SBATCH --gpus-per-task=1
# Choose the queue
#SBATCH --qos=dg
# Set the initial working directory

#SBATCH --chdir=/scratch/ecm2953/code/dwarf-p-cloudsc/build-gpu

export APGI_ACC_CUDA_HEAPSIZE=8GB

gridpoints=( 40000 80000 120000 160000 )
for ngptot in "${gridpoints[@]}"
do
    echo "----------------------------------------------------------------- NGPTOT=$ngptot"
    srun dwarf-cloudsc-gpu-scc-hoist 1 $ngptot 128
done
