#!/bin/bash --garcia
#SBATCH --partition=ompi*
#SBATCH --job-name=testslurm
#SBATCH --output=test.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00


mpirun -np 1 ./Arepo param.txt