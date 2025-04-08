#!/bin/bash --garcia
#SBATCH --partition=ompi*
#SBATCH --job-name=testslurm
#SBATCH --output=test.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00

make clean
make CONFIG=./spherical1d/polytrope_1d_spherical/Config.sh
python3 ./spherical1d/polytrope_1d_spherical/create.py
cd spherical1d/polytrope_1d_spherical
mpirun -np 1 ./Arepo param.txt
python3 ./check.py