#!/bin/bash

#SBATCH --job-name=formew
#SBATCH --output=/z/users.../logfile-%j.log
#SBATCH --export=ALL
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --time 0-01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

module load gnu_comp/11.1.0
module load cmake/3.25.1

pip3 install --user meson
pip3 install --user ninja

module load compilers/gcc/11.2.0
module load mpi/openmpi/4.1.1/gcc-11.2.0

module load apps/anaconda3/2023.03-poetry

source /add path to form in cosma
export FORMTMP=/tmp

tform -w$SLURM_CPUS_PER_TASK -l my-form-script.frm

