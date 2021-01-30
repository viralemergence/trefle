#! /bin/bash
#SBATCH --time=02:00:00
#SBATCH --array=1-829
#SBATCH --mem-per-cpu=2300M
#SBATCH --cpus-per-task=40
#SBATCH --job-name=trefle-prediction
#SBATCH --output=%x-%a.out

module load StdEnv/2020 julia/1.5.2

julia --project -t 38 02_imputation.jl
