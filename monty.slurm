#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --sockets-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00
#SBATCH --account=coursesuwashpc
#SBATCH --output=output.out
#SBATCH --job-name=Montecarlo

srun ./monty