#!/bin/bash
#
#SBATCH --array=1-1000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=12:00:00
#SBATCH --job-name=z_linear
#SBATCH --output=./out/sim_%A_%a.out
#SBATCH --error=./err/sim_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=marco.morucci@nyu.edu

module purge
module load r/intel/4.0.4
module load cplex/20.1.0

cd /scratch/$USER/linear_sim
export LIBRARY_PATH=$LIBRARY_PATH:/home/mm7936/GLPK/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mm7936/GLPK/lib
export CPATH=$CPATH:/home/mm7936/GLPK/include

/share/apps/r/4.0.4/intel/bin/Rscript --vanilla ./R/run_sim.R z_comparison \
-o "data/z_linear" \ 
-d "gen_linear"\
-s "50 100 250 500"\
-c "20"\
-a "0.0 0.2 0.5 0.8 1.2 2"\
-p "0.2"\
-v "0.1"